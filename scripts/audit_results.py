#!/usr/bin/env python3
"""Read-only audit of frozen lerch-campaign-v1 / lerch-range-v1 evidence.

Regenerates every prime with a Python segmented sieve, not Fermat quotients.
The stored moments are checked for canonical ranges and algebraic consistency;
this is not an independent arithmetic repetition of the search.
"""

from __future__ import annotations

import argparse
from contextlib import nullcontext
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
import multiprocessing
from pathlib import Path
import re
import struct
import sys
import time
from typing import Any
import zlib


MAX_CHUNK_BYTES = 64 * 1024 * 1024
ROW_FIELDS = {"p", "primitive_root", "q1", "q2", "lerch_remainder", "is_lerch"}
_worker_settings: tuple[Path, str, int, list[int]] | None = None


class AuditError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def integer(value: Any, low: int, high: int, label: str) -> int:
    require(type(value) is int and low <= value <= high, f"invalid {label}")
    return value


def seconds(value: Any, label: str, *, positive: bool = False) -> float:
    require(
        type(value) in (int, float) and math.isfinite(value)
        and (value > 0 if positive else value >= 0),
        f"invalid {label}",
    )
    return float(value)


def json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def same(actual: Any, expected: Any, label: str) -> None:
    require(json_bytes(actual) == json_bytes(expected), f"{label} mismatch")


def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def read_json(raw: bytes) -> Any:
    return json.loads(raw, object_pairs_hook=unique_object)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def regular_file(root: Path, relative: str) -> Path:
    require(isinstance(relative, str), "invalid evidence path type")
    parts = Path(relative).parts
    require(
        bool(parts) and not Path(relative).is_absolute()
        and all(part not in (".", "..") for part in parts),
        "unsafe evidence path",
    )
    path = root
    for part in parts:
        path = path / part
        require(not path.is_symlink(), f"symlink is not evidence: {relative}")
    require(path.is_file(), f"missing evidence file: {relative}")
    return path


def base_primes(limit: int) -> list[int]:
    flags = bytearray(b"\1") * (limit + 1)
    flags[:2] = b"\0" * min(2, len(flags))
    for p in range(2, math.isqrt(limit) + 1):
        if flags[p]:
            flags[p * p::p] = b"\0" * ((limit - p * p) // p + 1)
    return [p for p, present in enumerate(flags) if present]


def interval_primes(low: int, high: int, bases: list[int]) -> list[int]:
    flags = bytearray(b"\1") * (high - low + 1)
    for p in bases:
        if p * p > high:
            break
        first = max(p * p, ((low + p - 1) // p) * p)
        if first <= high:
            flags[first - low::p] = b"\0" * ((high - first) // p + 1)
    return [low + i for i, present in enumerate(flags) if present and low + i >= 2]


def canonical_summary(rows: list[dict[str, Any]], primes: list[int]) -> dict[str, Any]:
    require(len(rows) == len(primes), "independent sieve prime-count mismatch")
    digest = hashlib.sha256(b"lerch-canonical-v1\0")
    hits = []
    for row, p in zip(rows, primes):
        require(type(row) is dict and set(row) == ROW_FIELDS, "invalid canonical row fields")
        integer(row["p"], p, p, "prime order/coverage")
        integer(row["primitive_root"], 1, p - 1, f"primitive-root range at {p}")
        q1 = integer(row["q1"], 0, p - 1, f"q1 at {p}")
        q2 = integer(row["q2"], 0, p - 1, f"q2 at {p}")
        remainder = row["lerch_remainder"]
        require(type(row["is_lerch"]) is bool, f"invalid Lerch flag at {p}")
        if p == 2:
            require(remainder is None and q1 == q2 == 0 and not row["is_lerch"],
                    "invalid special row for 2")
        else:
            integer(remainder, 0, p - 1, f"Lerch remainder at {p}")
            expected = ((q1 * q1 + q2 - 2 * q1) * ((p + 1) // 2)) % p
            require(remainder == expected, f"moment/Lerch identity mismatch at {p}")
            require(row["is_lerch"] == (remainder == 0), f"Lerch flag mismatch at {p}")
        digest.update(struct.pack(
            "<QQQQBQB", p, row["primitive_root"], q1, q2,
            remainder is not None, remainder if remainder is not None else 0, row["is_lerch"],
        ))
        if row["is_lerch"]:
            hits.append(row)
    odd_terms = sum(p - 1 for p in primes if p != 2)
    return {
        "primes": len(primes),
        "first_prime": primes[0] if primes else None,
        "last_prime": primes[-1] if primes else None,
        "residue_terms": sum(p - 1 for p in primes),
        "odd_recurrence_terms": odd_terms,
        "primary_pair_steps": odd_terms // 2,
        "hits": hits,
        "sha256": digest.hexdigest(),
    }


def initialize_worker(root: Path, method: str, threads: int, bases: list[int]) -> None:
    global _worker_settings
    _worker_settings = root, method, threads, bases


def audit_chunk(record: dict[str, Any]) -> dict[str, Any]:
    assert _worker_settings is not None
    root, method, threads, bases = _worker_settings
    name = record["archive"]
    path = regular_file(root, name)
    size = integer(record["archive_bytes"], 1, MAX_CHUNK_BYTES, "archive size")
    require(path.stat().st_size == size, f"archive size mismatch: {name}")
    packed = path.read_bytes()
    require(hashlib.sha256(packed).hexdigest() == record["archive_sha256"],
            f"archive SHA-256 mismatch: {name}")
    with gzip.GzipFile(fileobj=io.BytesIO(packed)) as source:
        raw = source.read(MAX_CHUNK_BYTES + 1)
    require(len(raw) <= MAX_CHUNK_BYTES, f"oversized uncompressed archive: {name}")
    archive = read_json(raw)
    require(archive["format"] == "lerch-range-v1", f"unsupported chunk format: {name}")
    same(archive["interval"], {"start": record["start"], "end": record["end"]}, "interval")
    integer(archive["requested_batch"], 8192, 8192, "requested batch")
    sample = archive["sample"]
    require(sample["method"] == method, f"method mismatch: {name}")
    integer(sample["threads"], threads, threads, "worker count")
    require(type(sample["full_oracle_tuples_matched"]) is bool, "invalid oracle flag")
    primes = interval_primes(record["start"], record["end"], bases)
    summary = canonical_summary(archive["canonical"], primes)
    same(sample["summary"], summary, f"archive summary: {name}")
    same(record["summary"], summary, f"manifest chunk summary: {name}")
    timing = sample["timing"]
    compute = seconds(timing["seconds"], "compute time", positive=True)
    require(compute == seconds(record["compute_seconds"], "manifest compute time"),
            f"compute time mismatch: {name}")
    parts = [seconds(timing[key], key) for key in (
        "pool_sieve_seconds", "validation_setup_traversal_collection_seconds",
        "pool_teardown_seconds",
    )]
    require(math.isclose(math.fsum(parts), compute, rel_tol=1e-12, abs_tol=1e-8),
            f"timing components mismatch: {name}")
    seconds(record["wall_seconds"], "chunk wall time")
    require(type(record["recovered_after_interruption"]) is bool, "invalid recovery flag")
    return {
        "summary": summary,
        "uncompressed_bytes": len(raw),
        "full_oracle_tuples_matched": sample["full_oracle_tuples_matched"],
    }


def utc_timestamp(value: str) -> float:
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    require(parsed.tzinfo is not None, "timestamp must include timezone")
    return parsed.timestamp()


def validate_progress(root: Path, manifest: dict[str, Any]) -> None:
    progress = manifest["progress"]
    saved = [read_json(line) for line in regular_file(root, "progress.jsonl").read_bytes().splitlines()]
    same(saved, progress, "progress.jsonl")
    expected_log = "".join(
        f"{row['utc']} completed=[{row['start']},{row['end']}] "
        f"primes={row['primes']} terms={row['residue_terms']} "
        f"candidate_hits={row['candidate_hits']} block_wall_seconds={row['wall_seconds']:.3f} "
        f"cumulative_wall_seconds={row['cumulative_wall_seconds']:.3f}\n"
        for row in progress
    )
    require(regular_file(root, "progress.log").read_text() == expected_log, "progress.log mismatch")
    config = manifest["config"]
    cursor = config["start"]
    chunk_index = 0
    started = utc_timestamp(manifest["search_started_utc"])
    previous_time = started
    for block in progress:
        high = min(cursor + config["progress_width"] - 1, config["end"])
        integer(block["start"], cursor, cursor, "progress start")
        integer(block["end"], high, high, "progress end")
        records = []
        while chunk_index < len(manifest["chunks"]):
            record = manifest["chunks"][chunk_index]
            if record["end"] > high:
                break
            records.append(record)
            chunk_index += 1
        require(records and records[0]["start"] == cursor and records[-1]["end"] == high,
                "progress block does not partition complete chunks")
        same(block["primes"], sum(row["summary"]["primes"] for row in records), "progress primes")
        same(block["residue_terms"], sum(row["summary"]["residue_terms"] for row in records),
             "progress residue terms")
        same(block["candidate_hits"], [hit["p"] for row in records for hit in row["summary"]["hits"]],
             "progress candidates")
        completed = seconds(block["completed_timestamp"], "progress timestamp")
        cumulative = seconds(block["cumulative_wall_seconds"], "cumulative wall time")
        wall = seconds(block["wall_seconds"], "progress wall time")
        require(completed >= previous_time and all((
            abs(utc_timestamp(block["utc"]) - completed) < 0.05,
            abs(completed - started - cumulative) < 0.05,
            abs(completed - previous_time - wall) < 0.05,
        )), "progress timestamp/elapsed-time mismatch")
        previous_time = completed
        cursor = high + 1
    require(cursor == config["end"] + 1 and chunk_index == len(manifest["chunks"]),
            "progress coverage incomplete")


def audit(root: Path, workers: int = 1, expected_sha256: str | None = None) -> dict[str, Any]:
    started = time.monotonic()
    root = root.resolve(strict=True)
    integer(workers, 1, 64, "audit workers")
    raw_manifest = regular_file(root, "manifest.json").read_bytes()
    manifest_hash = hashlib.sha256(raw_manifest).hexdigest()
    if expected_sha256 is not None:
        require(manifest_hash == expected_sha256, "anchored manifest SHA-256 mismatch")
    manifest = read_json(raw_manifest)
    require(manifest["format"] == "lerch-campaign-v1", "unsupported manifest format")
    require(manifest["status"] == "complete", "campaign is not complete")
    config = manifest["config"]
    low = integer(config["start"], 2, 2_000_000_000, "campaign start")
    high = integer(config["end"], low, 2_000_000_000, "campaign end")
    width = integer(config["chunk_width"], 1, 100_000, "chunk width")
    block_width = integer(config["progress_width"], 1, 2_000_000_000, "progress width")
    threads = integer(config["threads"], 1, 10, "search workers")
    require(isinstance(manifest["method"], str), "invalid method")
    for key, length in (("source_sha", 40), ("binary_sha256", 64), ("runner_sha256", 64)):
        require(re.fullmatch(f"[0-9a-f]{{{length}}}", config[key]) is not None,
                f"invalid provenance {key}")
    require(type(manifest["chunks"]) is list and manifest["chunks"], "empty chunk list")
    cursor = low
    for record in manifest["chunks"]:
        block_end = min(low + ((cursor - low) // block_width + 1) * block_width - 1, high)
        expected_end = min(cursor + width - 1, block_end, high)
        require(cursor <= high, "extra or duplicate completed interval")
        integer(record["start"], cursor, cursor, "contiguous chunk start")
        integer(record["end"], expected_end, expected_end, "complete chunk end")
        require(record["archive"] == f"chunks/chunk-{cursor:010d}-{expected_end:010d}.json.gz",
                "noncanonical archive path")
        require(re.fullmatch("[0-9a-f]{64}", record["archive_sha256"]) is not None,
                "invalid archive digest")
        cursor = expected_end + 1
    require(cursor == high + 1, "chunk coverage does not reach endpoint")
    integer(manifest["next_start"], cursor, cursor, "completed coverage cursor")

    bases = base_primes(math.isqrt(high))
    settings = (root, manifest["method"], threads, bases)
    initialize_worker(*settings)
    context = (
        nullcontext(None) if workers == 1 else multiprocessing.get_context("spawn").Pool(
            workers, initializer=initialize_worker, initargs=settings,
        )
    )
    primes = terms = odd_terms = pairs = raw_bytes = oracle_chunks = 0
    candidates = []
    archive_index = hashlib.sha256(b"lerch-archive-index-v1\0")
    canonical_index = hashlib.sha256(b"lerch-canonical-index-v1\0")
    with context as pool:
        results = map(audit_chunk, manifest["chunks"]) if pool is None else pool.imap(
            audit_chunk, manifest["chunks"], chunksize=1,
        )
        for index, (record, result) in enumerate(zip(manifest["chunks"], results), 1):
            summary = result["summary"]
            primes += summary["primes"]
            terms += summary["residue_terms"]
            odd_terms += summary["odd_recurrence_terms"]
            pairs += summary["primary_pair_steps"]
            raw_bytes += result["uncompressed_bytes"]
            oracle_chunks += result["full_oracle_tuples_matched"]
            candidates.extend({
                **hit, "chunk_start": record["start"], "chunk_end": record["end"],
                "verification_status": "pending-independent-verification",
            } for hit in summary["hits"])
            archive_index.update(f"{record['archive_sha256']}  {record['archive']}\n".encode())
            canonical_index.update(
                struct.pack("<QQ", record["start"], record["end"])
                + bytes.fromhex(summary["sha256"])
            )
            if index % 1000 == 0:
                print(f"audited {index}/{len(manifest['chunks'])} chunks; {primes} primes",
                      file=sys.stderr, flush=True)
    same(manifest["primes"], primes, "campaign prime total")
    same(manifest["residue_terms"], terms, "campaign residue total")
    same(manifest["candidate_hits"], candidates, "campaign candidate list")
    validate_progress(root, manifest)
    require(regular_file(root, "manifest.json").read_bytes() == raw_manifest,
            "manifest changed during audit")
    return {
        "format": "lerch-evidence-audit-v1",
        "audit": "passed",
        "audited_utc": datetime.now(timezone.utc).isoformat(),
        "auditor_sha256": sha256_file(Path(__file__)),
        "manifest_sha256": manifest_hash,
        "manifest_digest_anchored": expected_sha256 is not None,
        "campaign_format": manifest["format"],
        "interval": {"start": low, "end": high},
        "chunks": len(manifest["chunks"]),
        "progress_blocks": len(manifest["progress"]),
        "primes": primes,
        "residue_terms": terms,
        "odd_recurrence_terms": odd_terms,
        "primary_pair_steps": pairs,
        "candidate_hits": candidates,
        "archive_bytes": sum(row["archive_bytes"] for row in manifest["chunks"]),
        "uncompressed_json_bytes": raw_bytes,
        "archive_index_sha256": archive_index.hexdigest(),
        "canonical_index_sha256": canonical_index.hexdigest(),
        "sum_chunk_compute_seconds": math.fsum(row["compute_seconds"] for row in manifest["chunks"]),
        "sum_chunk_wall_seconds": math.fsum(row["wall_seconds"] for row in manifest["chunks"]),
        "search_started_utc": manifest["search_started_utc"],
        "last_progress_utc": manifest["progress"][-1]["utc"],
        "campaign_wall_seconds": manifest["progress"][-1]["cumulative_wall_seconds"],
        "method": manifest["method"],
        "search_workers": threads,
        "source_sha": config["source_sha"],
        "binary_sha256": config["binary_sha256"],
        "runner_sha256": config["runner_sha256"],
        "independent_prime_regeneration": True,
        "moment_remainder_identity_checked": True,
        "fermat_quotients_recomputed": False,
        "primitive_roots_certified": False,
        "chunks_marked_full_oracle_matched": oracle_chunks,
        "audit_workers": workers,
        "audit_elapsed_seconds": time.monotonic() - started,
        "python": sys.version.split()[0],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--manifest-sha256")
    parser.add_argument("--output", type=Path, help="new report file outside the input directory")
    args = parser.parse_args()
    try:
        if args.output is not None:
            output = args.output.resolve()
            require(not output.is_relative_to(args.directory.resolve()), "report must be outside input")
            require(not output.exists(), "report output already exists")
        report = audit(args.directory, args.workers, args.manifest_sha256)
        encoded = json.dumps(report, indent=2, allow_nan=False) + "\n"
        if args.output is not None:
            with args.output.open("x") as target:
                target.write(encoded)
        print(encoded, end="")
    except (OSError, ValueError, KeyError, TypeError, EOFError, zlib.error) as error:
        parser.exit(1, f"audit failed: {error}\n")


if __name__ == "__main__":
    main()
