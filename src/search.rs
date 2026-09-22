//! One bounded range worker shared by benchmarking and the checkpoint supervisor.

use crate::arith::mul_mod;
use crate::doubling_cycles::{Backend, DoublingCycleContext};
use crate::moments::{Canonical, DEFAULT_BATCH_SIZE, MAX_PRIME};
use crate::sieve::{integer_sqrt, segmented_primes, simple_primes};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, Read, Write};
use std::path::Path;
use std::time::Instant;

pub const FORMAT: &str = "lerch-range-v2";
pub const MAX_WIDTH: u64 = 100_001;
pub const MAX_THREADS: usize = 256;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Interval {
    pub start: u64,
    pub end: u64,
}

impl Interval {
    pub fn validate(self, threads: usize) -> Result<(), String> {
        if !(2..=MAX_PRIME).contains(&self.start) || self.end < self.start || self.end > MAX_PRIME {
            return Err(format!(
                "inclusive bounds must satisfy 2 <= start <= end <= {MAX_PRIME}"
            ));
        }
        if self.end - self.start >= MAX_WIDTH {
            return Err(format!(
                "range width must not exceed {MAX_WIDTH}; use search for multiple chunks"
            ));
        }
        if !(1..=MAX_THREADS).contains(&threads) {
            return Err(format!("threads must be in 1..={MAX_THREADS}"));
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Provenance {
    pub package_version: String,
    pub backend: Backend,
    pub kernel: String,
    pub max_prime: u64,
    pub binary_sha256: String,
}

impl Provenance {
    pub fn capture(backend: Backend) -> Result<Self, String> {
        backend.validate_platform()?;
        let path = std::env::current_exe().map_err(|error| format!("locate binary: {error}"))?;
        let mut file =
            File::open(path).map_err(|error| format!("open binary for provenance: {error}"))?;
        let mut hash = Sha256::new();
        let mut buffer = [0u8; 65_536];
        loop {
            let count = file
                .read(&mut buffer)
                .map_err(|error| format!("hash binary: {error}"))?;
            if count == 0 {
                break;
            }
            hash.update(&buffer[..count]);
        }
        Ok(Self {
            package_version: env!("CARGO_PKG_VERSION").into(),
            backend,
            kernel: backend.kernel().into(),
            max_prime: MAX_PRIME,
            binary_sha256: format!("{:x}", hash.finalize()),
        })
    }
}

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Summary {
    pub primes: usize,
    pub first_prime: Option<u64>,
    pub last_prime: Option<u64>,
    pub residue_terms: u64,
    pub odd_recurrence_terms: u64,
    pub primary_pair_steps: u64,
    pub hits: Vec<Canonical>,
    pub sha256: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Timing {
    pub seconds: f64,
    pub pool_sieve_seconds: f64,
    pub validation_setup_traversal_collection_seconds: f64,
    pub pool_teardown_seconds: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Sample {
    pub method: String,
    pub threads: usize,
    pub timing: Timing,
    pub summary: Summary,
    pub full_oracle_tuples_matched: bool,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Archive {
    pub format: String,
    pub interval: Interval,
    pub requested_batch: u64,
    pub sample: Sample,
    pub canonical: Vec<Canonical>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub provenance: Option<Provenance>,
}

pub fn canonical_sha256(results: &[Canonical]) -> String {
    let mut hash = Sha256::new();
    hash.update(b"lerch-canonical-v1\0");
    for row in results {
        for value in [row.p, row.primitive_root, row.q1, row.q2] {
            hash.update(value.to_le_bytes());
        }
        hash.update([u8::from(row.lerch_remainder.is_some())]);
        hash.update(row.lerch_remainder.unwrap_or(0).to_le_bytes());
        hash.update([u8::from(row.is_lerch)]);
    }
    format!("{:x}", hash.finalize())
}

fn summarize(results: &[Canonical], primes: &[u64]) -> Result<Summary, String> {
    if results.len() != primes.len()
        || results.iter().zip(primes).any(|(row, &p)| {
            row.p != p
                || row.primitive_root == 0
                || row.primitive_root >= p
                || row.q1 >= p
                || row.q2 >= p
                || if p == 2 {
                    *row != Canonical::two()
                } else {
                    let expression = (mul_mod(row.q1, row.q1, p) + row.q2 + 2 * (p - row.q1)) % p;
                    let remainder = mul_mod(expression, p.div_ceil(2), p);
                    row.lerch_remainder != Some(remainder) || row.is_lerch != (remainder == 0)
                }
        })
    {
        return Err("range coverage/order/canonical tuple or Lerch identity mismatch".into());
    }
    let odd_recurrence_terms = results
        .iter()
        .filter(|row| row.p != 2)
        .map(|row| row.p - 1)
        .sum::<u64>();
    Ok(Summary {
        primes: results.len(),
        first_prime: results.first().map(|row| row.p),
        last_prime: results.last().map(|row| row.p),
        residue_terms: results.iter().map(|row| row.p - 1).sum(),
        odd_recurrence_terms,
        primary_pair_steps: odd_recurrence_terms / 2,
        hits: results.iter().filter(|row| row.is_lerch).cloned().collect(),
        sha256: canonical_sha256(results),
    })
}

pub fn measure(interval: Interval, backend: Backend, threads: usize) -> Result<Archive, String> {
    interval.validate(threads)?;
    let provenance = Provenance::capture(backend)?;
    let timer = Instant::now();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|error| format!("thread pool: {error}"))?;
    let base = simple_primes(integer_sqrt(interval.end));
    let primes = segmented_primes(interval.start, interval.end, &base);
    let pool_sieve = timer.elapsed();
    let results = pool.install(|| {
        primes
            .par_iter()
            .map(|&p| {
                if p == 2 {
                    Ok(Canonical::two())
                } else {
                    DoublingCycleContext::from_sieved_prime(p, DEFAULT_BATCH_SIZE, &base)?
                        .check(backend)
                }
            })
            .collect::<Result<Vec<_>, String>>()
    })?;
    std::hint::black_box(&results);
    let collected = timer.elapsed();
    drop(pool);
    let elapsed = timer.elapsed();
    let summary = summarize(&results, &primes)?;
    Ok(Archive {
        format: FORMAT.into(),
        interval,
        requested_batch: DEFAULT_BATCH_SIZE,
        sample: Sample {
            method: backend.method().into(),
            threads,
            timing: Timing {
                seconds: elapsed.as_secs_f64(),
                pool_sieve_seconds: pool_sieve.as_secs_f64(),
                validation_setup_traversal_collection_seconds: (collected - pool_sieve)
                    .as_secs_f64(),
                pool_teardown_seconds: (elapsed - collected).as_secs_f64(),
            },
            summary,
            full_oracle_tuples_matched: false,
        },
        canonical: results,
        provenance: Some(provenance),
    })
}

/// Read-only comparison accepts original v1 archives; it never upgrades their files.
/// Equality to another SIMD archive is a regression check, not independent verification.
pub fn compare_archive(result: &mut Archive, path: &Path) -> Result<(), String> {
    let file =
        File::open(path).map_err(|error| format!("open oracle {}: {error}", path.display()))?;
    let oracle: Archive = serde_json::from_reader(BufReader::new(file))
        .map_err(|error| format!("read oracle {}: {error}", path.display()))?;
    oracle.interval.validate(oracle.sample.threads)?;
    let primes = segmented_primes(
        oracle.interval.start,
        oracle.interval.end,
        &simple_primes(integer_sqrt(oracle.interval.end)),
    );
    if !["lerch-range-v1", FORMAT].contains(&oracle.format.as_str())
        || oracle.interval != result.interval
        || oracle.requested_batch != DEFAULT_BATCH_SIZE
        || (oracle.format == FORMAT && oracle.provenance.is_none())
        || summarize(&oracle.canonical, &primes)? != oracle.sample.summary
        || !oracle.sample.timing.seconds.is_finite()
        || oracle.sample.timing.seconds <= 0.0
    {
        return Err("oracle must be an intact archive for this exact interval".into());
    }
    if oracle.canonical != result.canonical {
        let mismatch = result
            .canonical
            .iter()
            .zip(&oracle.canonical)
            .find(|(got, expected)| got != expected);
        return Err(format!(
            "canonical tuple mismatch: lengths {}/{}, first={mismatch:?}",
            result.canonical.len(),
            oracle.canonical.len()
        ));
    }
    result.sample.full_oracle_tuples_matched = true;
    Ok(())
}

pub fn create_output(path: &Path) -> Result<File, String> {
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent).map_err(|error| format!("create output directory: {error}"))?;
    }
    OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|error| format!("create new output {}: {error}", path.display()))
}

pub fn persist<T: Serialize>(output: &mut File, value: &T) -> Result<(), String> {
    serde_json::to_writer(&mut *output, value)
        .map_err(|error| format!("serialize output: {error}"))?;
    output
        .write_all(b"\n")
        .and_then(|()| output.flush())
        .and_then(|()| output.sync_all())
        .map_err(|error| format!("persist output: {error}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn range_bounds_are_checked_before_allocation() {
        for interval in [
            Interval { start: 0, end: 10 },
            Interval { start: 10, end: 9 },
            Interval {
                start: 2,
                end: MAX_PRIME + 1,
            },
            Interval {
                start: 2,
                end: 100_003,
            },
        ] {
            assert!(interval.validate(1).is_err());
        }
        assert!(
            Interval {
                start: 2,
                end: 100_002
            }
            .validate(1)
            .is_ok()
        );
        assert!(Interval { start: 2, end: 2 }.validate(0).is_err());
    }

    #[test]
    fn canonical_summary_rejects_changed_identity_and_incomplete_coverage() {
        let rows = vec![Canonical::two(), crate::reference::generic(3).unwrap()];
        let summary = summarize(&rows, &[2, 3]).unwrap();
        assert_eq!(summary.primes, 2);
        assert_eq!(summary.primary_pair_steps, 1);
        assert!(summarize(&rows, &[2, 3, 5]).is_err());
        let mut bad = rows.clone();
        bad[1].lerch_remainder = Some(1);
        assert!(summarize(&bad, &[2, 3]).is_err());
        bad = rows;
        bad[0].is_lerch = true;
        assert!(summarize(&bad, &[2, 3]).is_err());
    }
}
