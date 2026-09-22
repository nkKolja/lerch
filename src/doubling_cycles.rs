//! All-prime doubling-coset traversal, shared by the two production SIMD kernels.

use crate::arith::{
    distinct_prime_factors, fermat_quotient_mod_p, inverse_mod, pow_mod, prime_base,
    primitive_root_from_factors,
};
use crate::carry_words::CarryWord32;
use crate::moments::{
    Canonical, DEFAULT_BATCH_SIZE, MAX_PRIME, PairedSums, Segment, add, finish, square_batch_size,
};
use crate::reduction::Montgomery32;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Backend {
    Neon,
    Avx512,
}

impl Backend {
    pub fn detect() -> Result<Self, String> {
        for backend in [Self::Neon, Self::Avx512] {
            if backend.validate_platform().is_ok() {
                return Ok(backend);
            }
        }
        Err("no supported SIMD backend: need AArch64 NEON or x86_64 AVX-512F; use reference instead".into())
    }

    pub fn parse(name: &str) -> Result<Self, String> {
        let backend = match name {
            "auto" => return Self::detect(),
            "neon" => Self::Neon,
            "avx512" => Self::Avx512,
            _ => {
                return Err(format!(
                    "unknown backend {name}; choose auto, neon or avx512"
                ));
            }
        };
        backend.validate_platform()?;
        Ok(backend)
    }

    pub fn validate_platform(self) -> Result<(), String> {
        match self {
            Self::Neon => crate::neon::check_platform(),
            Self::Avx512 => crate::avx512::check_platform(),
        }
    }

    pub fn method(self) -> &'static str {
        match self {
            Self::Neon => "carry32",
            Self::Avx512 => "avx512-64",
        }
    }

    pub fn kernel(self) -> &'static str {
        match self {
            Self::Neon => "neon16-centered-carry32-division-r64",
            Self::Avx512 => "avx512-64-centered-paired-r64",
        }
    }
}

#[derive(Clone, Copy)]
pub struct DoublingCycleContext {
    p: u64,
    order: u64,
    cycles: u64,
    g: u64,
    inverse_g: u64,
    quotient_g: u64,
    quotient_two: u64,
    batch: u64,
    reducer: Montgomery32,
}

impl DoublingCycleContext {
    pub fn new(p: u64, requested_batch: u64) -> Result<Self, String> {
        if !(3..=MAX_PRIME).contains(&p) {
            return Err(format!("odd prime must be in 3..={MAX_PRIME}"));
        }
        let base = prime_base(p)?;
        Self::from_sieved_prime(p, requested_batch, &base)
    }

    /// Only the inclusive sieve or the checked public constructor may call this.
    pub(crate) fn from_sieved_prime(
        p: u64,
        requested_batch: u64,
        base: &[u64],
    ) -> Result<Self, String> {
        let batch = square_batch_size(p, requested_batch)?;
        let factors = distinct_prime_factors(p - 1, base);
        let mut order = p - 1;
        for &factor in &factors {
            while order.is_multiple_of(factor) && pow_mod(2, order / factor, p) == 1 {
                order /= factor;
            }
        }
        let g = primitive_root_from_factors(p, &factors);
        let reducer = Montgomery32::new(p);
        let quotient_two = fermat_quotient_mod_p(2, p);
        Ok(Self {
            p,
            order,
            cycles: (p - 1) / order,
            g,
            inverse_g: inverse_mod(g, p),
            quotient_g: if g == 2 {
                quotient_two
            } else {
                fermat_quotient_mod_p(g, p)
            },
            quotient_two: reducer.encode(quotient_two),
            batch,
            reducer,
        })
    }

    pub fn p(&self) -> u64 {
        self.p
    }
    pub fn order(&self) -> u64 {
        self.order
    }
    pub fn cycles(&self) -> u64 {
        self.cycles
    }
    pub fn primitive_root(&self) -> u64 {
        self.g
    }
    pub fn batch_size(&self) -> u64 {
        self.batch
    }

    fn representatives(&self) -> Representatives<'_> {
        Representatives {
            context: self,
            remaining: self.cycles,
            len: self.order,
            c: 1,
            v: 1,
            u: 0,
        }
    }

    fn paired_representatives(&self) -> Representatives<'_> {
        let mut representatives = self.representatives();
        if self.order.is_multiple_of(2) {
            representatives.len /= 2;
        } else {
            representatives.remaining /= 2;
        }
        representatives
    }

    fn segments<const LANES: usize>(&self, seed: Segment) -> [Segment; LANES] {
        let lanes = LANES as u64;
        std::array::from_fn(|lane| {
            let lane = lane as u64;
            let start = lane * (seed.len / lanes) + lane.min(seed.len % lanes);
            let len = seed.len / lanes + u64::from(lane < seed.len % lanes);
            if start == 0 || len == 0 {
                Segment { len, ..seed }
            } else {
                let c = seed.c * pow_mod(2, start, self.p) % self.p;
                Segment {
                    c,
                    v: self.reducer.encode(inverse_mod(c, self.p)),
                    u: self.reducer.encode(fermat_quotient_mod_p(c, self.p)),
                    len,
                }
            }
        })
    }

    fn traverse<const LANES: usize>(
        &self,
        mut kernel: impl FnMut(&[Segment; LANES], u64) -> Result<PairedSums, String>,
    ) -> Result<Canonical, String> {
        let mut total = PairedSums::default();
        let mut merge = |partial: PairedSums| {
            total.first = add(total.first, partial.first, self.p);
            total.square = add(total.square, partial.square, self.p);
            total.pairs += partial.pairs;
        };
        let mut representatives = self.paired_representatives();
        let len = representatives.len;
        for _ in 0..representatives.remaining / LANES as u64 {
            let seeds = std::array::from_fn(|_| representatives.next().expect("full coset group"));
            merge(kernel(&seeds, LANES as u64 * len)?);
        }
        // Short coset tails use the same SIMD kernel with balanced segments, not
        // another scalar production implementation or an unsupported-CPU fallback.
        for seed in representatives {
            merge(kernel(&self.segments(seed), len)?);
        }
        finish(self.p, self.g, self.reducer, total)
    }

    pub fn check(&self, backend: Backend) -> Result<Canonical, String> {
        backend.validate_platform()?;
        match backend {
            Backend::Neon => {
                let words = CarryWord32::new(self.p);
                self.traverse(|seeds, pairs| {
                    crate::neon::paired_sums(
                        words,
                        self.batch,
                        self.reducer,
                        self.quotient_two,
                        seeds,
                        pairs,
                    )
                })
            }
            Backend::Avx512 => self.traverse(|seeds, pairs| {
                crate::avx512::paired_sums(
                    self.p,
                    self.batch,
                    self.reducer,
                    self.quotient_two,
                    seeds,
                    pairs,
                )
            }),
        }
    }
}

pub fn check_prime(p: u64, backend: Backend) -> Result<Canonical, String> {
    backend.validate_platform()?;
    if p == 2 {
        return Ok(Canonical::two());
    }
    DoublingCycleContext::new(p, DEFAULT_BATCH_SIZE)?.check(backend)
}

/// Canonical generic-g recurrence across coset representatives only.
/// It avoids per-coset exponentiation when the doubling order is small.
struct Representatives<'a> {
    context: &'a DoublingCycleContext,
    remaining: u64,
    len: u64,
    c: u64,
    v: u64,
    u: u64,
}

impl Iterator for Representatives<'_> {
    type Item = Segment;

    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            return None;
        }
        let ctx = self.context;
        let seed = Segment {
            c: self.c,
            v: ctx.reducer.encode(self.v),
            u: ctx.reducer.encode(self.u),
            len: self.len,
        };
        self.remaining -= 1;
        if self.remaining != 0 {
            let product = ctx.g * self.c;
            let carry = product / ctx.p;
            self.c = product - carry * ctx.p;
            self.v = self.v * ctx.inverse_g % ctx.p;
            self.u = (self.u + ctx.quotient_g + carry * self.v) % ctx.p;
        }
        Some(seed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::{direct_invariants, generic};
    use crate::sieve::{integer_sqrt, segmented_primes, simple_primes};

    fn native_backends() -> Vec<Backend> {
        [Backend::Neon, Backend::Avx512]
            .into_iter()
            .filter(|backend| backend.validate_platform().is_ok())
            .collect()
    }

    #[test]
    fn all_primes_through_2000_match_definitions_and_cover_sign_pairs() {
        let base = simple_primes(integer_sqrt(2000));
        let (mut even_order, mut odd_order, mut root_two, mut other_root) = (0, 0, 0, 0);
        for p in segmented_primes(3, 2000, &base) {
            let direct = direct_invariants(p, false, false);
            let reference = generic(p).unwrap();
            for batch in [1, 4, 18, 31, 32, 33, 460, 461, DEFAULT_BATCH_SIZE] {
                let ctx = DoublingCycleContext::new(p, batch).unwrap();
                for backend in native_backends() {
                    let got = ctx.check(backend).unwrap();
                    assert_eq!(got, reference, "p={p}, batch={batch}, backend={backend:?}");
                    assert_eq!(
                        (got.q1, got.q2, got.lerch_remainder),
                        (direct.q1, direct.q2, Some(direct.lerch_remainder))
                    );
                }
            }
            let ctx = DoublingCycleContext::new(p, DEFAULT_BATCH_SIZE).unwrap();
            if ctx.order.is_multiple_of(2) {
                even_order += 1;
            } else {
                odd_order += 1;
            }
            if ctx.cycles == 1 {
                root_two += 1;
            } else {
                other_root += 1;
            }
            let mut seen = vec![false; p as usize];
            let mut pairs = 0;
            for seed in ctx.paired_representatives() {
                let (mut c, mut v, mut u) = (seed.c, seed.v, seed.u);
                for _ in 0..seed.len {
                    for residue in [c, p - c] {
                        assert!(!seen[residue as usize], "duplicate pair p={p}, c={c}");
                        seen[residue as usize] = true;
                    }
                    assert_eq!(ctx.reducer.decode(u), fermat_quotient_mod_p(c, p));
                    assert_eq!(
                        ctx.reducer.decode(add(u, v, p)),
                        fermat_quotient_mod_p(p - c, p)
                    );
                    assert_eq!(c * ctx.reducer.decode(v) % p, 1);
                    let carry = 2 * c >= p;
                    c = 2 * c % p;
                    v = crate::moments::half(v, p);
                    u = add(u, ctx.quotient_two, p);
                    if carry {
                        u = add(u, v, p);
                    }
                    pairs += 1;
                }
            }
            assert_eq!(2 * pairs, p - 1);
            assert!(seen[1..].iter().all(|&x| x));
            assert!(!seen[0]);
        }
        assert!(even_order > 0 && odd_order > 0 && root_two > 50 && other_root > 50);
    }

    #[test]
    fn high_coset_counts_and_known_hits() {
        for (p, order, cycles) in [
            (3, 2, 1),
            (7, 3, 2),
            (17, 8, 2),
            (31, 5, 6),
            (43, 14, 3),
            (127, 7, 18),
            (8191, 13, 630),
            (524_287, 19, 27_594),
        ] {
            let ctx = DoublingCycleContext::new(p, u64::MAX).unwrap();
            assert_eq!((ctx.order, ctx.cycles), (order, cycles));
            let reference = generic(p).unwrap();
            for backend in native_backends() {
                assert_eq!(ctx.check(backend).unwrap(), reference);
            }
        }
        for p in [3, 103, 839, 2237, 42_447_347] {
            for backend in native_backends() {
                let got = check_prime(p, backend).unwrap();
                assert!(got.is_lerch);
                assert_eq!(got.lerch_remainder, Some(0));
                if p == 42_447_347 {
                    assert_eq!(
                        (got.primitive_root, got.q1, got.q2),
                        (2, 34_227_565, 10_415_263)
                    );
                }
            }
        }
    }

    #[test]
    fn invalid_inputs_and_unsupported_backends_are_rejected() {
        for p in [0, 1, 2, MAX_PRIME + 1, u64::MAX] {
            assert!(DoublingCycleContext::new(p, DEFAULT_BATCH_SIZE).is_err());
        }
        for p in [
            4,
            9,
            15,
            341,
            561,
            1105,
            2047,
            1_000_000_000,
            1_000_000_001,
            MAX_PRIME,
        ] {
            assert!(matches!(DoublingCycleContext::new(p, DEFAULT_BATCH_SIZE),
                Err(error) if error == format!("{p} is not prime")));
        }
        assert!(DoublingCycleContext::new(7, 0).is_err());
        let ctx = DoublingCycleContext::new(1_999_999_973, DEFAULT_BATCH_SIZE).unwrap();
        assert_eq!(ctx.batch_size(), 4);
        for backend in [Backend::Neon, Backend::Avx512] {
            if backend.validate_platform().is_err() {
                assert!(ctx.check(backend).is_err());
            }
        }
        assert!(Backend::parse("avx2").is_err());
    }

    fn shifted<const LANES: usize>(
        kernel: impl Fn(&DoublingCycleContext, u64, &[Segment; LANES], u64) -> PairedSums,
    ) {
        for p in [
            103,
            2237,
            42_447_347,
            999_999_893,
            1_000_000_007,
            1_073_741_827,
            1_431_655_777,
            1_499_999_957,
            1_999_999_811,
            1_999_999_973,
        ] {
            let ctx = DoublingCycleContext::new(p, DEFAULT_BATCH_SIZE).unwrap();
            for len in [
                0, 1, 15, 16, 17, 31, 32, 33, 63, 64, 65, 287, 288, 289, 460, 461, 493, 511, 512,
                513, 529, 1041,
            ] {
                let len = len.min((p - 1) / 2);
                for start in [0, 1, 31, 32, 33, p / 2, p - 1 - len] {
                    if start + len > p - 1 {
                        continue;
                    }
                    let c = pow_mod(2, start, p);
                    let seed = Segment {
                        c,
                        v: ctx.reducer.encode(inverse_mod(c, p)),
                        u: ctx.reducer.encode(fermat_quotient_mod_p(c, p)),
                        len,
                    };
                    let segments = ctx.segments::<LANES>(seed);
                    let (mut first, mut square) = (0, 0);
                    for exponent in start..start + len {
                        let c = pow_mod(2, exponent, p);
                        let x = (2 * fermat_quotient_mod_p(c, p) + inverse_mod(c, p)) % p;
                        first = (first + x) % p;
                        square = (square + x * x) % p;
                    }
                    let expected = (ctx.reducer.encode(first), ctx.reducer.encode(square), len);
                    for requested in [1, 4, 18, 31, 32, 33, 460, 461, 8192] {
                        let batch = square_batch_size(p, requested).unwrap();
                        let got = kernel(&ctx, batch, &segments, len);
                        assert_eq!(
                            (got.first, got.square, got.pairs),
                            expected,
                            "p={p}, start={start}, len={len}, lanes={LANES}, batch={batch}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn shifted_native_segments_match_centered_definitions_through_two_billion() {
        for backend in native_backends() {
            match backend {
                Backend::Neon => shifted::<16>(|ctx, batch, segments, len| {
                    crate::neon::paired_sums(
                        CarryWord32::new(ctx.p),
                        batch,
                        ctx.reducer,
                        ctx.quotient_two,
                        segments,
                        len,
                    )
                    .unwrap()
                }),
                Backend::Avx512 => shifted::<64>(|ctx, batch, segments, len| {
                    crate::avx512::paired_sums(
                        ctx.p,
                        batch,
                        ctx.reducer,
                        ctx.quotient_two,
                        segments,
                        len,
                    )
                    .unwrap()
                }),
            }
        }
    }

    #[test]
    #[ignore = "explicit bounded native run: two complete large original-generic traversals"]
    fn two_billion_full_moments_match_original_generic() {
        let backend = Backend::detect().expect("requires a native supported SIMD CPU");
        for p in [1_499_999_957, 1_999_999_973] {
            assert!((2..=integer_sqrt(p)).all(|divisor| p % divisor != 0));
            let timer = std::time::Instant::now();
            let reference = generic(p).unwrap();
            println!(
                "original-generic p={p} tuple={reference:?} seconds={:.9}",
                timer.elapsed().as_secs_f64()
            );
            let timer = std::time::Instant::now();
            let ctx = DoublingCycleContext::new(p, DEFAULT_BATCH_SIZE).unwrap();
            let got = ctx.check(backend).unwrap();
            assert_eq!(got, reference);
            if p == 1_999_999_973 {
                assert_eq!(
                    (got.q1, got.q2, got.lerch_remainder),
                    (30_303_682, 1_772_269_194, Some(1_433_408_083))
                );
            }
            println!(
                "full-tuple-match p={p} backend={backend:?} order={} cycles={} batch={} seconds={:.9}",
                ctx.order(),
                ctx.cycles(),
                ctx.batch_size(),
                timer.elapsed().as_secs_f64()
            );
        }
    }
}
