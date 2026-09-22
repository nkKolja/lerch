"""Bounded fake-worker tests: no remote host or large arithmetic campaign."""

import argparse
from datetime import datetime, timedelta, timezone
import fcntl
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import time
import unittest


RUNNER = Path(__file__).resolve().parents[1] / "scripts/search.py"
SPEC = importlib.util.spec_from_file_location("search", RUNNER)
search = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(search)


def archive(low, high, backend, threads, binary_sha):
    rows = []
    for p in search.prime_list(low, high):
        hit = p in (3, 103)
        rows.append({
            "p": p, "primitive_root": 1 if p == 2 else 2,
            "q1": int(hit), "q2": 0 if p == 2 else (1 if hit else 2),
            "lerch_remainder": None if p == 2 else int(not hit), "is_lerch": hit,
        })
    digest = hashlib.sha256(b"lerch-canonical-v1\0")
    for row in rows:
        digest.update(struct.pack(
            "<QQQQBQB", row["p"], row["primitive_root"], row["q1"], row["q2"],
            int(row["lerch_remainder"] is not None), row["lerch_remainder"] or 0,
            int(row["is_lerch"]),
        ))
    odd_terms = sum(row["p"] - 1 for row in rows if row["p"] != 2)
    return {
        "format": "lerch-range-v2", "interval": {"start": low, "end": high},
        "requested_batch": 8192, "canonical": rows,
        "provenance": {
            "backend": backend, "kernel": search.KERNELS[backend], "max_prime": 2_000_000_000,
            "binary_sha256": binary_sha, "package_version": "fake-worker-not-native",
        },
        "sample": {
            "method": search.METHODS[backend], "threads": threads,
            "timing": {
                "seconds": 0.001, "pool_sieve_seconds": 0.0001,
                "validation_setup_traversal_collection_seconds": 0.0008,
                "pool_teardown_seconds": 0.0001,
            },
            "full_oracle_tuples_matched": False,
            "summary": {
                "primes": len(rows), "first_prime": rows[0]["p"] if rows else None,
                "last_prime": rows[-1]["p"] if rows else None,
                "residue_terms": sum(row["p"] - 1 for row in rows),
                "odd_recurrence_terms": odd_terms, "primary_pair_steps": odd_terms // 2,
                "hits": [row for row in rows if row["is_lerch"]], "sha256": digest.hexdigest(),
            },
        },
    }


def fake_worker():
    control = Path(os.environ["LERCH_FAKE_CONTROL"])
    settings = json.loads(control.read_text())
    if sys.argv[1] != "benchmark":
        raise RuntimeError("supervisor did not call the primary bounded worker")
    flags = dict(zip(sys.argv[2::2], sys.argv[3::2]))
    low, high = int(flags["--start"]), int(flags["--end"])
    with control.with_suffix(".calls").open("a") as calls:
        calls.write(json.dumps({"pid": os.getpid(), "start": low}) + "\n")
    output = Path(flags["--output"])
    with output.open("x") as file:
        file.write('{"unfinished":')
    print("fake worker started", flush=True)
    time.sleep(settings.get("delay_by_start", {}).get(str(low), settings.get("delay", 0)))
    if settings.get("exit_code"):
        sys.exit(settings["exit_code"])
    result = archive(low, high, flags["--backend"], int(flags["--threads"]),
                     search.sha256_file(Path(sys.argv[0])))
    if settings.get("bad_identity"):
        result["canonical"][0]["lerch_remainder"] = 1
    output.write_text(json.dumps(result))


class SearchSupervisorTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="lerch-supervisor-test-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.output = self.root / "search"
        self.control = self.root / "control.json"
        self.control.write_text("{}")
        self.binary = self.root / "fake-worker"
        self.binary.write_text(
            f"#!{sys.executable}\nimport runpy\n"
            f"runpy.run_path({str(Path(__file__).resolve())!r})['fake_worker']()\n"
        )
        self.binary.chmod(0o755)
        self.env = dict(os.environ, LERCH_FAKE_CONTROL=str(self.control),
                        PYTHONDONTWRITEBYTECODE="1")

    def deadline(self, seconds=1.0):
        return (datetime.now(timezone.utc) + timedelta(seconds=seconds)).isoformat()

    def run_runner(self, *options, start=3, end=7, output=None):
        return subprocess.run([
            sys.executable, "-B", str(RUNNER),
            "--binary", str(self.binary), "--output-dir", str(output or self.output),
            "--source-sha", "fake-worker-not-a-native-result", "--backend", "neon",
            "--start", str(start), "--end", str(end), "--chunk-size", "2",
            "--progress-width", "10", "--threads", "1", "--min-free-bytes", "0",
            *options,
        ], env=self.env, capture_output=True, text=True, timeout=10)

    def manifest(self, output=None):
        return json.loads(((output or self.output) / "manifest.json").read_text())

    def calls(self):
        path = self.control.with_suffix(".calls")
        return [json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []

    def assert_ok(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_absolute_timezone_parsing(self):
        self.assertEqual(search.absolute_utc("2026-09-23T17:38:00+02:00"),
                         search.absolute_utc("2026-09-23T15:38:00Z"))
        for value in ["2026-09-23T15:38:00", "tomorrow", "NaN", ""]:
            with self.assertRaises(argparse.ArgumentTypeError):
                search.absolute_utc(value)

    def test_expired_startup_and_resume_never_reset_the_deadline(self):
        self.assert_ok(self.run_runner("--deadline-utc", self.deadline(-1)))
        first = self.manifest()
        self.assertEqual(first["status"], "deadline_reached")
        self.assertEqual((first["next_start"], first["chunks"]), (3, []))
        self.assertIsNone(first["search_started_utc"])
        self.assertEqual(self.calls(), [])
        self.assert_ok(self.run_runner("--resume"))
        self.assertEqual(self.manifest()["config"], first["config"])
        self.assertEqual(self.manifest()["deadline_reached_utc"], first["deadline_reached_utc"])
        self.assertEqual(self.calls(), [])
        before = (self.output / "manifest.json").read_bytes()
        rejected = self.run_runner("--resume", "--deadline-utc", self.deadline(60))
        self.assertNotEqual(rejected.returncode, 0)
        self.assertIn("refusing to resume different", rejected.stderr)
        self.assertEqual((self.output / "manifest.json").read_bytes(), before)

    def test_active_chunk_is_killed_reaped_and_never_checkpointed(self):
        self.control.write_text(json.dumps({"delay": 4}))
        started = time.monotonic()
        self.assert_ok(self.run_runner("--deadline-utc", self.deadline()))
        self.assertLess(time.monotonic() - started, 3.5)
        manifest = self.manifest()
        self.assertEqual(manifest["status"], "deadline_reached")
        self.assertEqual((manifest["next_start"], manifest["primes"], manifest["chunks"]), (3, 0, []))
        self.assertEqual(list((self.output / "chunks").glob("*.json.gz")), [])
        partials = list((self.output / "chunks").glob("*.partial-*.json"))
        self.assertEqual(len(partials), 1)
        self.assertEqual(partials[0].read_text(), '{"unfinished":')
        self.assertEqual(len(self.calls()), 1)
        with self.assertRaises(ProcessLookupError):
            os.kill(self.calls()[0]["pid"], 0)
        self.assertIn("binary timed out", manifest["deadline_reason"])

    def test_completed_prefix_survives_expiry_and_separate_authorized_continuation(self):
        self.control.write_text(json.dumps({"delay_by_start": {"5": 4}}))
        self.assert_ok(self.run_runner("--deadline-utc", self.deadline()))
        first = self.manifest()
        self.assertEqual((first["status"], first["next_start"]), ("deadline_reached", 5))
        self.assertEqual([(row["start"], row["end"]) for row in first["chunks"]], [(3, 4)])
        self.assert_ok(self.run_runner("--resume"))
        self.assertEqual(len(self.calls()), 2)
        self.control.write_text("{}")
        continuation = self.root / "explicit-new-budget"
        self.assert_ok(self.run_runner(
            "--deadline-utc", self.deadline(10), start=5, output=continuation,
        ))
        self.assertEqual((self.manifest(continuation)["status"], self.manifest(continuation)["next_start"]),
                         ("complete", 8))
        self.assertEqual(self.manifest()["next_start"], 5)

    def test_default_budget_is_persisted_at_24_hours_and_explicit_resume_is_required(self):
        self.assert_ok(self.run_runner("--max-chunks", "1"))
        first = self.manifest()
        self.assertEqual((first["status"], first["next_start"]), ("paused", 5))
        self.assertAlmostEqual(search.timestamp(first["config"]["deadline_utc"])
                               - search.timestamp(first["created_utc"]), 86_400)
        before = (self.output / "manifest.json").read_bytes()
        self.assertNotEqual(self.run_runner().returncode, 0)
        self.assertEqual((self.output / "manifest.json").read_bytes(), before)
        self.assert_ok(self.run_runner("--resume"))
        final = self.manifest()
        self.assertEqual(final["config"], first["config"])
        self.assertEqual((final["status"], final["next_start"]), ("complete", 8))
        self.assertEqual(final["candidate_hits"][0]["verification_status"],
                         "pending-independent-verification")

    def test_deadline_over_24_hours_is_rejected_before_creating_output(self):
        result = self.run_runner("--deadline-utc", self.deadline(25 * 3600))
        self.assertEqual(result.returncode, 2)
        self.assertIn("24-hour", result.stderr)
        self.assertFalse(self.output.exists())
        self.assertEqual(self.calls(), [])

    def test_expired_completed_search_stays_complete(self):
        deadline = self.deadline(0.6)
        self.assert_ok(self.run_runner("--deadline-utc", deadline, end=3))
        before = (self.output / "manifest.json").read_bytes()
        time.sleep(max(0, search.timestamp(deadline) - time.time()) + 0.02)
        self.assert_ok(self.run_runner("--resume", end=3))
        self.assertEqual((self.output / "manifest.json").read_bytes(), before)
        self.assertEqual(len(self.calls()), 1)

    def test_two_billion_neon_requires_explicit_deadline_and_is_not_legacy_capped(self):
        rejected = self.run_runner(start=1_999_999_973, end=1_999_999_973)
        self.assertEqual(rejected.returncode, 2)
        self.assertIn("require --deadline-utc", rejected.stderr)
        self.assertEqual(self.calls(), [])
        self.assert_ok(self.run_runner(
            "--deadline-utc", self.deadline(10), start=1_999_999_973, end=1_999_999_973,
        ))
        self.assertEqual(self.manifest()["status"], "complete")
        self.assertEqual(self.manifest()["primes"], 1)

    def test_invalid_bounds_deadlines_and_backends_fail_before_worker(self):
        for options in [
            ["--end", "2000000001"], ["--end", "18446744073709551615"],
            ["--deadline-utc", "2026-09-23T15:38:00"], ["--deadline-utc", "invalid"],
            ["--backend", "avx2"], ["--chunk-size", "100001"], ["--threads", "0"],
        ]:
            self.assertEqual(self.run_runner(*options).returncode, 2)
        self.assertEqual(self.calls(), [])

    def test_worker_failure_is_reported_not_disguised_as_deadline_success(self):
        self.control.write_text(json.dumps({"exit_code": 7}))
        result = self.run_runner("--deadline-utc", self.deadline(10))
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.manifest()["status"], "failed")
        self.assertIn("command exited 7", self.manifest()["error"])
        self.assertEqual(self.manifest()["next_start"], 3)

    def test_invalid_worker_tuple_cannot_become_a_checkpoint(self):
        self.control.write_text(json.dumps({"bad_identity": True}))
        result = self.run_runner("--deadline-utc", self.deadline(10))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("invalid canonical tuple", result.stderr)
        self.assertEqual(self.manifest()["chunks"], [])
        self.assertEqual(list((self.output / "chunks").glob("*.gz")), [])

    def test_lock_prevents_duplicate_workers(self):
        self.assert_ok(self.run_runner("--max-chunks", "1"))
        before = (self.output / "manifest.json").read_bytes()
        with (self.output / "search.lock").open("a+") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            result = self.run_runner("--resume")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual((self.output / "manifest.json").read_bytes(), before)
        self.assertEqual(len(self.calls()), 1)

    def test_corrupt_checkpoint_and_manifest_totals_fail_without_rewriting(self):
        self.assert_ok(self.run_runner("--max-chunks", "1"))
        manifest = self.manifest()
        path = self.output / "manifest.json"
        before = path.read_bytes()
        chunk = self.output / manifest["chunks"][0]["archive"]
        original = chunk.read_bytes()
        chunk.write_bytes(b"corrupt")
        result = self.run_runner("--resume")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("compressed checkpoint checksum mismatch", result.stderr)
        self.assertEqual(path.read_bytes(), before)
        chunk.write_bytes(original)
        manifest["primes"] += 1
        path.write_text(json.dumps(manifest))
        before = path.read_bytes()
        result = self.run_runner("--resume")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("cursor/totals/candidate", result.stderr)
        self.assertEqual(path.read_bytes(), before)
        self.assertEqual(len(self.calls()), 1)

    def test_committed_archive_is_recovered_after_manifest_interruption(self):
        self.assert_ok(self.run_runner("--max-chunks", "1"))
        manifest = self.manifest()
        archive_path = self.output / manifest["chunks"][0]["archive"]
        original = archive_path.read_bytes()
        manifest.update({"next_start": 3, "primes": 0, "residue_terms": 0,
                         "chunks": [], "candidate_hits": [], "progress": []})
        (self.output / "manifest.json").write_text(json.dumps(manifest))
        self.assert_ok(self.run_runner("--resume", "--max-chunks", "1"))
        recovered = self.manifest()
        self.assertTrue(recovered["chunks"][0]["recovered_after_interruption"])
        self.assertEqual(archive_path.read_bytes(), original)
        self.assertEqual(len(self.calls()), 1)
        self.assertEqual(recovered["next_start"], 5)

    def test_different_configuration_and_binary_are_rejected(self):
        self.assert_ok(self.run_runner("--max-chunks", "1"))
        path = self.output / "manifest.json"
        before = path.read_bytes()
        for options in [["--threads", "2"], ["--source-sha", "different"], ["--backend", "avx512"]]:
            result = self.run_runner("--resume", *options)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("refusing to resume different", result.stderr)
            self.assertEqual(path.read_bytes(), before)
        with self.binary.open("a") as file:
            file.write("# different worker bytes\n")
        self.assertNotEqual(self.run_runner("--resume").returncode, 0)
        self.assertEqual(path.read_bytes(), before)
        self.assertEqual(len(self.calls()), 1)

    def test_archived_v1_resume_is_strictly_read_only(self):
        self.output.mkdir()
        old = b'{"format":"lerch-campaign-v1"}'
        path = self.output / "manifest.json"
        path.write_bytes(old)
        result = self.run_runner("--resume")
        self.assertEqual(result.returncode, 2)
        self.assertIn("read-only", result.stderr)
        self.assertEqual(path.read_bytes(), old)
        self.assertEqual([entry.name for entry in self.output.iterdir()], ["manifest.json"])

    def test_stop_file_and_disk_reserve_prevent_new_work(self):
        self.assert_ok(self.run_runner("--max-chunks", "1"))
        (self.output / "STOP").touch()
        self.assert_ok(self.run_runner("--resume"))
        self.assertEqual(self.manifest()["status"], "paused")
        self.assertEqual(len(self.calls()), 1)
        (self.output / "STOP").unlink()
        result = self.run_runner("--resume", "--min-free-bytes", str(2 ** 63))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("free disk below", result.stderr)
        self.assertEqual(len(self.calls()), 1)

    def test_create_only_checkpoint_write_preserves_existing_bytes(self):
        path = self.root / "checkpoint.gz"
        original = gzip.compress(b"immutable", mtime=0)
        path.write_bytes(original)
        with self.assertRaises(FileExistsError):
            search.atomic_bytes(path, b"replacement", replace=False)
        self.assertEqual(path.read_bytes(), original)
        self.assertEqual(list(self.root.glob("checkpoint.gz.tmp-*")), [])


if __name__ == "__main__":
    unittest.main()
