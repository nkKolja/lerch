#!/usr/bin/env python3
"""Small definition-derived fixtures and deliberately damaged evidence."""

from __future__ import annotations

import copy
import gzip
import hashlib
import json
import math
from pathlib import Path
import struct
import tempfile
import unittest

import audit_results as auditor


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, sort_keys=True) + "\n")


def fixture_rows(low: int, high: int) -> list[dict]:
    rows = []
    for p in range(low, high + 1):
        if p < 2 or any(p % d == 0 for d in range(2, p)):
            continue
        root = next(g for g in range(1, p) if len({pow(g, j, p) for j in range(p - 1)}) == p - 1)
        quotients = [(pow(a, p - 1, p * p) - 1) // p for a in range(1, p)]
        q1 = sum(quotients) % p if p != 2 else 0
        q2 = sum(q * q for q in quotients) % p if p != 2 else 0
        # Compute L from its definition, independently of the moment identity.
        numerator = sum(pow(a, p - 1) for a in range(1, p)) - math.factorial(p - 1) - p
        remainder = (numerator // (p * p)) % p if p != 2 else None
        rows.append({
            "p": p, "primitive_root": root, "q1": q1, "q2": q2,
            "lerch_remainder": remainder, "is_lerch": remainder == 0,
        })
    return rows


def fixture_summary(rows: list[dict]) -> dict:
    digest = hashlib.sha256(b"lerch-canonical-v1\0")
    for row in rows:
        for key in ("p", "primitive_root", "q1", "q2"):
            digest.update(struct.pack("<Q", row[key]))
        digest.update(bytes([row["lerch_remainder"] is not None]))
        digest.update(struct.pack("<Q", row["lerch_remainder"] or 0))
        digest.update(bytes([row["is_lerch"]]))
    terms = sum(row["p"] - 1 for row in rows)
    odd_terms = sum(row["p"] - 1 for row in rows if row["p"] != 2)
    return {
        "primes": len(rows), "first_prime": rows[0]["p"] if rows else None,
        "last_prime": rows[-1]["p"] if rows else None, "residue_terms": terms,
        "odd_recurrence_terms": odd_terms, "primary_pair_steps": odd_terms // 2,
        "hits": [row for row in rows if row["is_lerch"]], "sha256": digest.hexdigest(),
    }


def fixture(root: Path) -> dict:
    (root / "chunks").mkdir()
    records = []
    # Ends with an empty composite singleton, as the historical 1B campaign does.
    for low, high in ((2, 5), (6, 9), (10, 13), (14, 14)):
        rows = fixture_rows(low, high)
        summary = fixture_summary(rows)
        archive = {
            "format": "lerch-range-v1", "interval": {"start": low, "end": high},
            "requested_batch": 8192,
            "sample": {
                "method": "avx512-64", "threads": 8, "summary": summary,
                "full_oracle_tuples_matched": False,
                "timing": {
                    "seconds": 0.5, "pool_sieve_seconds": 0.1,
                    "validation_setup_traversal_collection_seconds": 0.3,
                    "pool_teardown_seconds": 0.1,
                },
            },
            "canonical": rows,
        }
        packed = gzip.compress(json.dumps(archive).encode(), mtime=0)
        name = f"chunks/chunk-{low:010d}-{high:010d}.json.gz"
        (root / name).write_bytes(packed)
        records.append({
            "start": low, "end": high, "archive": name, "archive_bytes": len(packed),
            "archive_sha256": hashlib.sha256(packed).hexdigest(), "summary": summary,
            "compute_seconds": 0.5, "wall_seconds": 0.6, "recovered_after_interruption": False,
            "completed_utc": "2026-01-01T00:00:01+00:00",
        })
    progress = []
    for i, group in enumerate((records[:2], records[2:]), 1):
        progress.append({
            "start": group[0]["start"], "end": group[-1]["end"],
            "primes": sum(r["summary"]["primes"] for r in group),
            "residue_terms": sum(r["summary"]["residue_terms"] for r in group),
            "candidate_hits": [hit["p"] for r in group for hit in r["summary"]["hits"]],
            "utc": f"2026-01-01T00:00:0{i}+00:00", "completed_timestamp": 1767225600 + i,
            "wall_seconds": 1.0, "cumulative_wall_seconds": float(i),
        })
    manifest = {
        "format": "lerch-campaign-v1", "status": "complete", "method": "avx512-64",
        "config": {
            "start": 2, "end": 14, "chunk_width": 4, "progress_width": 8, "threads": 8,
            "source_sha": "a" * 40, "binary_sha256": "b" * 64, "runner_sha256": "c" * 64,
        },
        "chunks": records, "next_start": 15, "progress": progress,
        "search_started_utc": "2026-01-01T00:00:00+00:00",
        "primes": sum(r["summary"]["primes"] for r in records),
        "residue_terms": sum(r["summary"]["residue_terms"] for r in records),
        "candidate_hits": [{
            **hit, "chunk_start": r["start"], "chunk_end": r["end"],
            "verification_status": "pending-independent-verification",
        } for r in records for hit in r["summary"]["hits"]],
    }
    write_json(root / "manifest.json", manifest)
    write_progress(root, manifest)
    return manifest


def write_progress(root: Path, manifest: dict) -> None:
    (root / "progress.jsonl").write_text("".join(json.dumps(row) + "\n" for row in manifest["progress"]))
    (root / "progress.log").write_text("".join(
        f"{r['utc']} completed=[{r['start']},{r['end']}] "
        f"primes={r['primes']} terms={r['residue_terms']} "
        f"candidate_hits={r['candidate_hits']} block_wall_seconds={r['wall_seconds']:.3f} "
        f"cumulative_wall_seconds={r['cumulative_wall_seconds']:.3f}\n"
        for r in manifest["progress"]
    ))


class CampaignAuditTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.manifest = fixture(self.root)

    def save(self) -> None:
        write_json(self.root / "manifest.json", self.manifest)

    def rewrite_chunk(self, mutate, *, refresh_summary: bool = False) -> None:
        record = self.manifest["chunks"][0]
        path = self.root / record["archive"]
        archive = json.loads(gzip.decompress(path.read_bytes()))
        mutate(archive)
        if refresh_summary:
            summary = fixture_summary(archive["canonical"])
            archive["sample"]["summary"] = summary
            record["summary"] = summary
        packed = gzip.compress(json.dumps(archive).encode(), mtime=0)
        path.write_bytes(packed)
        record["archive_bytes"] = len(packed)
        record["archive_sha256"] = hashlib.sha256(packed).hexdigest()
        self.save()

    def test_complete_and_parallel_are_equivalent(self) -> None:
        digest = auditor.sha256_file(self.root / "manifest.json")
        before = {p.relative_to(self.root): p.read_bytes() for p in self.root.rglob("*") if p.is_file()}
        serial = auditor.audit(self.root, expected_sha256=digest)
        parallel = auditor.audit(self.root, workers=2, expected_sha256=digest)
        for field in ("audit", "primes", "residue_terms", "candidate_hits",
                      "archive_index_sha256", "canonical_index_sha256"):
            self.assertEqual(serial[field], parallel[field])
        self.assertEqual(serial["primes"], 6)
        self.assertEqual([hit["p"] for hit in serial["candidate_hits"]], [3])
        self.assertEqual(serial["chunks"], 4)
        self.assertFalse(serial["fermat_quotients_recomputed"])
        self.assertEqual(before, {
            p.relative_to(self.root): p.read_bytes() for p in self.root.rglob("*") if p.is_file()
        })

    def test_bad_manifest_anchor(self) -> None:
        with self.assertRaisesRegex(auditor.AuditError, "anchored manifest"):
            auditor.audit(self.root, expected_sha256="0" * 64)

    def test_missing_archive(self) -> None:
        (self.root / self.manifest["chunks"][0]["archive"]).unlink()
        with self.assertRaisesRegex(auditor.AuditError, "missing evidence"):
            auditor.audit(self.root)

    def test_corrupt_gzip_checksum(self) -> None:
        path = self.root / self.manifest["chunks"][0]["archive"]
        raw = bytearray(path.read_bytes())
        raw[-1] ^= 1
        path.write_bytes(raw)
        with self.assertRaisesRegex(auditor.AuditError, "archive SHA-256"):
            auditor.audit(self.root)

    def test_gaps_overlaps_duplicates_and_missing_endpoint(self) -> None:
        original = copy.deepcopy(self.manifest)
        for mutation in (
            lambda m: m["chunks"].pop(1),
            lambda m: m["chunks"].insert(1, copy.deepcopy(m["chunks"][0])),
            lambda m: m["chunks"][1].update(start=5),
            lambda m: m["chunks"].pop(),
            lambda m: m.update(next_start=14),
        ):
            with self.subTest(mutation=mutation):
                self.manifest = copy.deepcopy(original)
                mutation(self.manifest)
                self.save()
                with self.assertRaises(auditor.AuditError):
                    auditor.audit(self.root)

    def test_reordered_or_missing_primes_even_with_refreshed_hashes(self) -> None:
        self.rewrite_chunk(lambda a: a["canonical"].reverse(), refresh_summary=True)
        with self.assertRaisesRegex(auditor.AuditError, "prime order/coverage"):
            auditor.audit(self.root)

    def test_missing_prime_even_with_refreshed_hashes(self) -> None:
        self.rewrite_chunk(lambda a: a["canonical"].pop(), refresh_summary=True)
        with self.assertRaisesRegex(auditor.AuditError, "prime-count"):
            auditor.audit(self.root)

    def test_lerch_identity_even_with_refreshed_hashes(self) -> None:
        self.rewrite_chunk(lambda a: a["canonical"][2].update(q2=0), refresh_summary=True)
        with self.assertRaisesRegex(auditor.AuditError, "moment/Lerch identity"):
            auditor.audit(self.root)

    def test_integer_is_not_boolean(self) -> None:
        self.rewrite_chunk(lambda a: a["canonical"][1].update(is_lerch=1))
        with self.assertRaisesRegex(auditor.AuditError, "Lerch flag"):
            auditor.audit(self.root)

    def test_canonical_hash_mismatch(self) -> None:
        self.rewrite_chunk(lambda a: a["canonical"][2].update(primitive_root=3))
        with self.assertRaisesRegex(auditor.AuditError, "archive summary"):
            auditor.audit(self.root)

    def test_aggregate_and_progress_mismatch(self) -> None:
        self.manifest["primes"] += 1
        self.save()
        with self.assertRaisesRegex(auditor.AuditError, "campaign prime total"):
            auditor.audit(self.root)
        self.manifest["primes"] -= 1
        self.manifest["progress"][0]["residue_terms"] += 1
        self.save()
        write_progress(self.root, self.manifest)
        with self.assertRaisesRegex(auditor.AuditError, "progress residue terms"):
            auditor.audit(self.root)

    def test_paths_and_symlinks_are_rejected(self) -> None:
        original = self.manifest["chunks"][0]["archive"]
        self.manifest["chunks"][0]["archive"] = "../outside.json.gz"
        self.save()
        with self.assertRaisesRegex(auditor.AuditError, "archive path"):
            auditor.audit(self.root)
        self.manifest["chunks"][0]["archive"] = original
        self.save()
        path = self.root / original
        moved = path.with_suffix(".original")
        path.rename(moved)
        path.symlink_to(moved.name)
        with self.assertRaisesRegex(auditor.AuditError, "symlink"):
            auditor.audit(self.root)

    def test_duplicate_json_keys_are_rejected(self) -> None:
        path = self.root / "manifest.json"
        path.write_text(path.read_text().replace('"status": "complete"', '"status":"paused","status":"complete"'))
        with self.assertRaisesRegex(auditor.AuditError, "duplicate JSON key"):
            auditor.audit(self.root)


if __name__ == "__main__":
    unittest.main()
