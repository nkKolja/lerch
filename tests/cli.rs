use lerch_prime_search::reference::generic;
use lerch_prime_search::search::Archive;
use lerch_prime_search::{Backend, Canonical};
use serde_json::{Value, json};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT_DIRECTORY: AtomicUsize = AtomicUsize::new(0);

struct Temporary(PathBuf);

impl Temporary {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "lerch-cli-{}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos(),
            NEXT_DIRECTORY.fetch_add(1, Ordering::Relaxed),
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for Temporary {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_lerch-prime-search"))
        .args(args)
        .output()
        .unwrap()
}

fn success(args: &[&str]) -> Value {
    let output = run(args);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap()
}

fn text(path: &Path) -> &str {
    path.to_str().unwrap()
}

fn native() -> Option<Backend> {
    Backend::detect().ok()
}

#[test]
fn reference_is_available_without_a_simd_backend_and_rejects_bad_inputs() {
    let result = success(&["reference", "--prime", "103"]);
    let canonical: Canonical = serde_json::from_value(result["canonical"].clone()).unwrap();
    assert_eq!(canonical, generic(103).unwrap());
    assert_eq!(result["definition_checked"], false);
    for args in [
        vec!["reference", "--prime", "9"],
        vec!["reference", "--prime", "2000000001"],
        vec!["reference", "--prime", "3", "--prime", "5"],
        vec!["reference", "--prime"],
        vec!["reference", "--prime", "3", "--unknown"],
        vec!["check", "--prime", "3", "--backend", "scalar"],
    ] {
        let output = run(&args);
        assert!(!output.status.success());
        assert!(!output.stderr.is_empty());
    }
}

#[test]
fn automatic_backend_really_runs_the_selected_kernel() {
    let Some(backend) = native() else {
        let output = run(&["check", "--prime", "103"]);
        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains("no supported SIMD backend"));
        return;
    };
    for p in [3, 103, 839, 2237, 42_447_347] {
        let result = success(&["check", "--prime", &p.to_string()]);
        let canonical: Canonical = serde_json::from_value(result["canonical"].clone()).unwrap();
        assert_eq!(canonical, generic(p).unwrap());
        assert!(canonical.is_lerch);
        assert_eq!(result["provenance"]["kernel"], backend.kernel());
        assert_eq!(
            result["provenance"]["binary_sha256"]
                .as_str()
                .unwrap()
                .len(),
            64
        );
    }
    let unsupported = match backend {
        Backend::NeonInverse | Backend::Neon => "avx512",
        Backend::Avx512 => "neon",
    };
    assert!(
        !run(&["check", "--prime", "103", "--backend", unsupported])
            .status
            .success()
    );
}

#[test]
fn neon_inverse_is_the_arm_default_and_explicit_neon_remains_the_baseline() {
    if Backend::NeonInverse.validate_platform().is_err() {
        let output = run(&["check", "--prime", "103", "--backend", "neon-inverse"]);
        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains("NEON requires"));
        return;
    }
    assert_eq!(Backend::detect().unwrap(), Backend::NeonInverse);
    for p in [2, 3, 5, 7, 17, 103, 257, 8191, 524_287] {
        let prime = p.to_string();
        let old = success(&["check", "--prime", &prime, "--backend", "neon"]);
        let new = success(&["check", "--prime", &prime, "--backend", "neon-inverse"]);
        let auto = success(&["check", "--prime", &prime]);
        assert_eq!(old["canonical"], new["canonical"], "p={p}");
        assert_eq!(auto["canonical"], new["canonical"], "p={p}");
        assert_eq!(old["provenance"]["kernel"], Backend::Neon.kernel());
        for result in [&auto, &new] {
            assert_eq!(result["provenance"]["backend"], "neon-inverse");
            assert_eq!(result["provenance"]["kernel"], "neon8-scaled-inverse-r64");
            assert_eq!(result["format"], "lerch-check-v2");
        }
    }
    let range = success(&[
        "benchmark",
        "--start",
        "97",
        "--end",
        "103",
        "--threads",
        "1",
        "--backend",
        "neon-inverse",
    ]);
    assert_eq!(range["sample"]["summary"]["primes"], 3);
    assert_eq!(range["sample"]["method"], "neon-inverse");
    assert_eq!(range["provenance"]["backend"], "neon-inverse");
    assert_eq!(range["provenance"]["kernel"], "neon8-scaled-inverse-r64");
    let validation = success(&[
        "validate",
        "--limit",
        "2000",
        "--bigint-limit",
        "200",
        "--backend",
        "neon-inverse",
    ]);
    assert_eq!(validation["validation"], "passed");
    assert_eq!(validation["backend"], "neon-inverse");
}

#[test]
fn benchmark_covers_every_prime_and_compares_archived_v1_read_only() {
    let Some(backend) = native() else {
        return;
    };
    let dir = Temporary::new();
    let output = dir.path("range.json");
    let summary = success(&[
        "benchmark",
        "--start",
        "2",
        "--end",
        "300",
        "--threads",
        "1",
        "--output",
        text(&output),
    ]);
    assert_eq!(summary["sample"]["method"], backend.method());
    let archive: Archive = serde_json::from_slice(&fs::read(&output).unwrap()).unwrap();
    assert_eq!(archive.format, "lerch-range-v2");
    assert_eq!(archive.sample.summary.primes, 62);
    assert_eq!(
        archive
            .sample
            .summary
            .hits
            .iter()
            .map(|row| row.p)
            .collect::<Vec<_>>(),
        [3, 103]
    );
    for row in &archive.canonical {
        assert_eq!(*row, generic(row.p).unwrap());
    }
    let mut old = serde_json::to_value(&archive).unwrap();
    old["format"] = json!("lerch-range-v1");
    old["sample"]["method"] = json!("generic");
    old.as_object_mut().unwrap().remove("provenance");
    let oracle = dir.path("original-v1.json");
    let original = serde_json::to_vec(&old).unwrap();
    fs::write(&oracle, &original).unwrap();
    let matched = success(&[
        "benchmark",
        "--start",
        "2",
        "--end",
        "300",
        "--threads",
        "1",
        "--oracle",
        text(&oracle),
    ]);
    assert_eq!(matched["sample"]["full_oracle_tuples_matched"], true);
    assert_eq!(fs::read(&oracle).unwrap(), original);
    old["canonical"][1]["q1"] = json!(0);
    fs::write(&oracle, serde_json::to_vec(&old).unwrap()).unwrap();
    assert!(
        !run(&[
            "benchmark",
            "--start",
            "2",
            "--end",
            "300",
            "--oracle",
            text(&oracle),
        ])
        .status
        .success()
    );
    let original = fs::read(&output).unwrap();
    assert!(
        !run(&[
            "benchmark",
            "--start",
            "2",
            "--end",
            "300",
            "--output",
            text(&output),
        ])
        .status
        .success()
    );
    assert_eq!(fs::read(&output).unwrap(), original);
}

#[test]
fn empty_and_prime_two_ranges_have_canonical_shapes() {
    if native().is_none() {
        return;
    }
    let two = success(&["benchmark", "--start", "2", "--end", "2", "--threads", "1"]);
    assert_eq!(two["sample"]["summary"]["primes"], 1);
    assert_eq!(two["sample"]["summary"]["primary_pair_steps"], 0);
    let empty = success(&[
        "benchmark",
        "--start",
        "14",
        "--end",
        "16",
        "--threads",
        "1",
    ]);
    assert_eq!(empty["sample"]["summary"]["primes"], 0);
    assert!(empty["sample"]["summary"]["first_prime"].is_null());
    assert!(
        !run(&["benchmark", "--start", "2", "--end", "100003"])
            .status
            .success()
    );
}

#[test]
fn generic_cross_check_cannot_claim_independent_candidate_verification() {
    if native().is_none() {
        return;
    }
    let generic = success(&["verify", "--prime", "103", "--generic"]);
    assert_eq!(generic["verification"]["verified"], true);
    assert_eq!(generic["verification"]["definition_checked"], false);
    assert_eq!(
        generic["verification"]["candidate_status"],
        "pending-independent-verification"
    );
    let definition = success(&["verify", "--prime", "103"]);
    assert_eq!(definition["verification"]["verified"], true);
    assert_eq!(definition["verification"]["definition_checked"], true);
    assert_eq!(
        definition["verification"]["candidate_status"],
        "verified-lerch"
    );
    assert_eq!(
        definition["verification"]["power_sum_minus_factorial_minus_p_mod_p3"],
        "0"
    );
    assert!(!run(&["verify", "--prime", "2"]).status.success());
}

#[test]
fn check_and_verify_outputs_never_overwrite_existing_files() {
    if native().is_none() {
        return;
    }
    let dir = Temporary::new();
    let output = dir.path("proof.json");
    fs::write(&output, "preserve this proof").unwrap();
    for command in ["check", "verify", "reference"] {
        assert!(
            !run(&[command, "--prime", "3", "--output", text(&output)])
                .status
                .success()
        );
        assert_eq!(fs::read_to_string(&output).unwrap(), "preserve this proof");
    }
}

#[test]
fn primary_search_uses_simd_and_resumes_without_resetting_budget_or_candidate_state() {
    let Some(backend) = native() else {
        return;
    };
    let dir = Temporary::new();
    let output = dir.path("search");
    let args = [
        "search",
        "--start",
        "2",
        "--end",
        "103",
        "--chunk-size",
        "30",
        "--threads",
        "1",
        "--min-free-bytes",
        "0",
        "--output-dir",
        text(&output),
    ];
    let mut paused = args.to_vec();
    paused.extend(["--max-chunks", "1"]);
    let result = run(&paused);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let path = output.join("manifest.json");
    let first: Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
    assert_eq!(first["format"], "lerch-campaign-v2");
    assert_eq!(first["status"], "paused");
    assert_eq!(first["next_start"], 32);
    assert_eq!(first["method"], backend.method());
    assert_eq!(
        first["config"]["backend"],
        serde_json::to_value(backend).unwrap()
    );
    assert_eq!(first["config"]["binary_sha256"].as_str().unwrap().len(), 64);
    assert_eq!(first["config"]["runner_sha256"].as_str().unwrap().len(), 64);
    assert_eq!(
        first["candidate_hits"][0]["verification_status"],
        "pending-independent-verification"
    );
    let before = fs::read(&path).unwrap();
    assert!(!run(&args).status.success());
    assert_eq!(fs::read(&path).unwrap(), before);
    let mut resume = args.to_vec();
    resume.push("--resume");
    let result = run(&resume);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let final_manifest: Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
    assert_eq!(final_manifest["status"], "complete");
    assert_eq!(final_manifest["next_start"], 104);
    assert_eq!(final_manifest["primes"], 27);
    assert_eq!(final_manifest["config"], first["config"]);
    let hits = final_manifest["candidate_hits"].as_array().unwrap();
    assert_eq!(
        hits.iter()
            .map(|row| row["p"].as_u64().unwrap())
            .collect::<Vec<_>>(),
        [3, 103]
    );
    assert!(
        hits.iter()
            .all(|row| row["verification_status"] == "pending-independent-verification")
    );
    let complete = fs::read(&path).unwrap();
    assert!(run(&resume).status.success());
    assert_eq!(fs::read(&path).unwrap(), complete);
}

#[test]
fn old_neon_checkpoint_cannot_silently_resume_with_the_new_auto_backend() {
    if Backend::NeonInverse.validate_platform().is_err() {
        return;
    }
    let dir = Temporary::new();
    let output = dir.path("carry32");
    let args = [
        "search",
        "--start",
        "2",
        "--end",
        "13",
        "--chunk-size",
        "5",
        "--threads",
        "1",
        "--min-free-bytes",
        "0",
        "--source-sha",
        "backend-resume-regression",
        "--output-dir",
        text(&output),
    ];
    let mut start = args.to_vec();
    start.extend(["--backend", "neon", "--max-chunks", "1"]);
    assert!(run(&start).status.success());
    let manifest = output.join("manifest.json");
    let before = fs::read(&manifest).unwrap();
    let first: Value = serde_json::from_slice(&before).unwrap();
    assert_eq!(first["config"]["backend"], "neon");
    assert_eq!(first["method"], "carry32");
    for backend in ["auto", "neon-inverse"] {
        let mut resume = args.to_vec();
        resume.extend(["--backend", backend, "--resume"]);
        let result = run(&resume);
        assert!(!result.status.success());
        assert!(String::from_utf8_lossy(&result.stderr).contains("refusing to resume different"));
        assert_eq!(fs::read(&manifest).unwrap(), before);
    }
    let mut resume = args.to_vec();
    resume.extend(["--backend", "neon", "--resume"]);
    assert!(run(&resume).status.success());
    let completed: Value = serde_json::from_slice(&fs::read(&manifest).unwrap()).unwrap();
    assert_eq!(completed["status"], "complete");
    assert_eq!(completed["primes"], 6);
    assert_eq!(completed["config"], first["config"]);
}

#[test]
fn archived_campaign_resume_is_refused_without_creating_any_files() {
    if native().is_none() {
        return;
    }
    let dir = Temporary::new();
    let output = dir.path("archive");
    fs::create_dir(&output).unwrap();
    let path = output.join("manifest.json");
    let old = b"{\"format\":\"lerch-campaign-v1\"}";
    fs::write(&path, old).unwrap();
    let result = run(&[
        "search",
        "--start",
        "3",
        "--end",
        "7",
        "--output-dir",
        text(&output),
        "--resume",
    ]);
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("read-only"));
    assert_eq!(fs::read(&path).unwrap(), old);
    assert_eq!(fs::read_dir(&output).unwrap().count(), 1);
}
