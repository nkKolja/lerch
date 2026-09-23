use lerch_prime_search::arith::fermat_quotient_mod_p;
use lerch_prime_search::reference::{direct_invariants, generic};
use lerch_prime_search::reference_recurrence::recurrence_values;
use lerch_prime_search::search::{
    Interval, Provenance, compare_archive, create_output, measure, persist,
};
use lerch_prime_search::sieve::{integer_sqrt, segmented_primes, simple_primes};
use lerch_prime_search::verify::{direct_lerch_remainder_bigint, verify};
use lerch_prime_search::{Backend, DEFAULT_BATCH_SIZE, DoublingCycleContext, check_prime};
use serde::Serialize;
use serde_json::json;
use sha2::{Digest, Sha256};
use std::collections::HashSet;
use std::path::Path;
use std::process::Command;
use std::time::Instant;

const SEARCH_RUNNER: &str = include_str!("../scripts/search.py");

fn usage() -> &'static str {
    "lerch-prime-search <command> [options]\n\
     All production checks support primes through 2000000000.\n\
     Backend: auto|neon-inverse|neon|avx512. Auto selects neon-inverse on ARM, avx512 on x86.\n\
     Explicit neon retains the previous NEON16 carry32 kernel.\n\
     commands:\n\
       check --prime P [--backend B] [--output FILE]\n\
       search --start N --end N --output-dir DIR [--backend B] [--threads N]\n\
              [--chunk-size N] [--resume] [--deadline-utc RFC3339] [--max-chunks N]\n\
              [--source-sha SHA] [--progress-width N] [--min-free-bytes N]\n\
       benchmark --start N --end N [--backend B] [--threads N] [--output FILE]\n\
                 [--oracle FILE]  (inclusive width <=100001, one pass, no warmup)\n\
       verify --prime P [--backend B] [--generic] [--output FILE]\n\
       reference --prime P [--output FILE]  (original generic O(p), no SIMD required)\n\
       validate [--limit N] [--bigint-limit N] [--backend B]\n\
     search uses Python 3 stdlib on Linux/macOS and stores a deadline of at most 24h.\n\
     Above 1B an explicit absolute deadline is required. Resume never extends it.\n\
     Resume pins the resolved backend, source, binary and runner; upgrades need a new output directory.\n\
     Search records candidates, NOT independent proofs. Default verify is expensive:\n\
     definitions plus bigint p^3 for Lerch candidates. --generic is only a cross-check.\n\
     Outputs are create-only; archived v1 data is never rewritten."
}

fn value(args: &[String], name: &str) -> Result<Option<String>, String> {
    let Some(index) = args.iter().position(|arg| arg == name) else {
        return Ok(None);
    };
    args.get(index + 1)
        .cloned()
        .map(Some)
        .ok_or_else(|| format!("missing value for {name}"))
}

fn parsed<T: std::str::FromStr>(args: &[String], name: &str, default: T) -> Result<T, String> {
    value(args, name)?
        .map(|text| text.parse().map_err(|_| format!("invalid {name}: {text}")))
        .unwrap_or(Ok(default))
}

fn required<T: std::str::FromStr>(args: &[String], name: &str) -> Result<T, String> {
    let text = value(args, name)?.ok_or_else(|| format!("requires {name}"))?;
    text.parse().map_err(|_| format!("invalid {name}: {text}"))
}

fn reject_unknown(args: &[String], valued: &[&str], flags: &[&str]) -> Result<(), String> {
    let mut seen = HashSet::new();
    let mut i = 0;
    while i < args.len() {
        if !seen.insert(&args[i]) {
            return Err(format!("duplicate option: {}", args[i]));
        }
        if valued.contains(&args[i].as_str()) {
            if args.get(i + 1).is_none_or(|value| value.starts_with("--")) {
                return Err(format!("missing value for {}", args[i]));
            }
            i += 2;
        } else if flags.contains(&args[i].as_str()) {
            i += 1;
        } else {
            return Err(format!("unknown option: {}", args[i]));
        }
    }
    Ok(())
}

fn backend(args: &[String]) -> Result<Backend, String> {
    Backend::parse(value(args, "--backend")?.as_deref().unwrap_or("auto"))
}

fn default_threads() -> usize {
    std::thread::available_parallelism().map_or(1, |n| n.get().min(8))
}

fn emit<T: Serialize>(args: &[String], record: &T) -> Result<(), String> {
    if let Some(path) = value(args, "--output")? {
        persist(&mut create_output(Path::new(&path))?, record)?;
    }
    println!(
        "{}",
        serde_json::to_string_pretty(record).map_err(|error| error.to_string())?
    );
    Ok(())
}

fn check(args: &[String]) -> Result<(), String> {
    reject_unknown(args, &["--prime", "--backend", "--output"], &[])?;
    let p = required(args, "--prime")?;
    let backend = backend(args)?;
    let provenance = Provenance::capture(backend)?;
    let timer = Instant::now();
    let (canonical, setup) = if p == 2 {
        (check_prime(p, backend)?, None)
    } else {
        let ctx = DoublingCycleContext::new(p, DEFAULT_BATCH_SIZE)?;
        let setup = json!({
            "order": ctx.order(), "cycles": ctx.cycles(),
            "batch_size": ctx.kernel_batch_size(backend),
        });
        (ctx.check(backend)?, Some(setup))
    };
    emit(
        args,
        &json!({
            "format": "lerch-check-v2", "provenance": provenance, "setup": setup,
            "canonical": canonical, "seconds": timer.elapsed().as_secs_f64(),
        }),
    )
}

fn benchmark(args: &[String]) -> Result<(), String> {
    reject_unknown(
        args,
        &[
            "--start",
            "--end",
            "--backend",
            "--threads",
            "--output",
            "--oracle",
        ],
        &[],
    )?;
    let interval = Interval {
        start: required(args, "--start")?,
        end: required(args, "--end")?,
    };
    let backend = backend(args)?;
    let threads = parsed(args, "--threads", default_threads())?;
    interval.validate(threads)?;
    let mut output = value(args, "--output")?
        .map(|path| create_output(Path::new(&path)))
        .transpose()?;
    let mut archive = measure(interval, backend, threads)?;
    if let Some(path) = value(args, "--oracle")? {
        compare_archive(&mut archive, Path::new(&path))?;
    }
    if let Some(output) = &mut output {
        persist(output, &archive)?;
    }
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "format": archive.format, "interval": archive.interval,
            "provenance": archive.provenance, "sample": archive.sample,
        }))
        .map_err(|error| error.to_string())?
    );
    Ok(())
}

fn reference(args: &[String]) -> Result<(), String> {
    reject_unknown(args, &["--prime", "--output"], &[])?;
    let p = required(args, "--prime")?;
    let timer = Instant::now();
    let canonical = generic(p)?;
    emit(
        args,
        &json!({
            "format": "lerch-reference-v2", "method": "original-generic-recurrence",
            "canonical": canonical, "seconds": timer.elapsed().as_secs_f64(),
            "definition_checked": false,
        }),
    )
}

fn verification(args: &[String]) -> Result<(), String> {
    reject_unknown(args, &["--prime", "--backend", "--output"], &["--generic"])?;
    let p = required(args, "--prime")?;
    let backend = backend(args)?;
    let provenance = Provenance::capture(backend)?;
    let timer = Instant::now();
    let canonical = check_prime(p, backend)?;
    let transcript = verify(&canonical, args.iter().any(|arg| arg == "--generic"))?;
    let verified = transcript.verified;
    emit(
        args,
        &json!({
            "format": "lerch-verification-v2", "provenance": provenance,
            "canonical": canonical, "verification": transcript, "seconds": timer.elapsed().as_secs_f64(),
        }),
    )?;
    if verified {
        Ok(())
    } else {
        Err("independent verification failed".into())
    }
}

fn validate(args: &[String]) -> Result<(), String> {
    reject_unknown(args, &["--limit", "--bigint-limit", "--backend"], &[])?;
    let limit = parsed(args, "--limit", 2000u64)?;
    let bigint_limit = parsed(args, "--bigint-limit", 200u64)?;
    if !(3..=100_000).contains(&limit) || bigint_limit > 1000 {
        return Err("validation requires limit in 3..=100000 and bigint-limit <=1000".into());
    }
    let backend = backend(args)?;
    let base = simple_primes(integer_sqrt(limit));
    let primes = segmented_primes(3, limit, &base);
    let mut hits = Vec::new();
    for &p in &primes {
        let fast = check_prime(p, backend)?;
        let direct = direct_invariants(p, false, false);
        if fast != generic(p)?
            || fast.q1 != direct.q1
            || fast.q2 != direct.q2
            || fast.q1 != direct.wilson
            || fast.lerch_remainder != Some(direct.lerch_remainder)
        {
            return Err(format!("aggregate definition mismatch at p={p}"));
        }
        for (a, q) in recurrence_values(p, &base) {
            if q != fermat_quotient_mod_p(a, p) {
                return Err(format!("q_p(a) mismatch at p={p}, a={a}"));
            }
        }
        if p <= bigint_limit && fast.lerch_remainder != Some(direct_lerch_remainder_bigint(p)) {
            return Err(format!("bigint definition mismatch at p={p}"));
        }
        if fast.is_lerch {
            hits.push(p);
        }
    }
    println!(
        "{}",
        json!({
            "validation": "passed", "backend": backend, "inclusive_limit": limit,
            "odd_primes_checked": primes.len(), "bigint_limit": bigint_limit.min(limit), "lerch_hits": hits,
        })
    );
    Ok(())
}

fn search(args: &[String]) -> Result<(), String> {
    reject_unknown(
        args,
        &[
            "--start",
            "--end",
            "--output-dir",
            "--backend",
            "--threads",
            "--chunk-size",
            "--progress-width",
            "--deadline-utc",
            "--max-chunks",
            "--source-sha",
            "--min-free-bytes",
        ],
        &["--resume"],
    )?;
    let backend = backend(args)?;
    let backend_name = match backend {
        Backend::NeonInverse => "neon-inverse",
        Backend::Neon => "neon",
        Backend::Avx512 => "avx512",
    };
    let binary =
        std::env::current_exe().map_err(|error| format!("locate worker binary: {error}"))?;
    let mut command = Command::new("python3");
    command
        .args(["-B", "-c", SEARCH_RUNNER, "--binary"])
        .arg(binary)
        .args(["--backend", backend_name, "--runner-sha256"])
        .arg(format!("{:x}", Sha256::digest(SEARCH_RUNNER.as_bytes())));
    if value(args, "--threads")?.is_none() {
        command.args(["--threads", &default_threads().to_string()]);
    }
    let mut i = 0;
    while i < args.len() {
        if args[i] == "--backend" {
            i += 2;
        } else {
            command.arg(&args[i]);
            i += 1;
        }
    }
    #[cfg(unix)]
    {
        use std::os::unix::process::CommandExt;
        Err(format!(
            "start Python 3 search supervisor: {}",
            command.exec()
        ))
    }
    #[cfg(not(unix))]
    {
        Err("checkpoint supervision currently requires Linux or macOS".into())
    }
}

fn run(args: &[String]) -> Result<(), String> {
    if args.is_empty()
        || ["help", "--help", "-h"].contains(&args[0].as_str())
        || (args.len() == 2 && ["--help", "-h"].contains(&args[1].as_str()))
    {
        println!("{}", usage());
        return Ok(());
    }
    match args[0].as_str() {
        "check" => check(&args[1..]),
        "search" => search(&args[1..]),
        "benchmark" => benchmark(&args[1..]),
        "reference" => reference(&args[1..]),
        "verify" => verification(&args[1..]),
        "validate" => validate(&args[1..]),
        other => Err(format!("unknown command: {other}\n{}", usage())),
    }
}

fn main() {
    if let Err(error) = run(&std::env::args().skip(1).collect::<Vec<_>>()) {
        eprintln!("error: {error}");
        std::process::exit(2);
    }
}
