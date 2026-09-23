#!/usr/bin/env python3
"""Portable checkpoint/deadline supervisor, also embedded in the Rust search command."""

import argparse
from datetime import datetime, timedelta, timezone
import fcntl
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import struct
import subprocess
import sys
import time
import uuid


FORMAT = "lerch-campaign-v2"
METHODS = {"neon-inverse": "neon-inverse", "neon": "carry32", "avx512": "avx512-64"}
KERNELS = {
    "neon-inverse": "neon8-scaled-inverse-r64",
    "neon": "neon16-centered-carry32-division-r64",
    "avx512": "avx512-64-centered-paired-r64",
}
STOP_REQUESTED = False


class DeadlineReached(RuntimeError):
    pass


def utc():
    return datetime.now(timezone.utc).isoformat()


def timestamp(value):
    return datetime.fromisoformat(value.replace("Z", "+00:00")).timestamp()


def absolute_utc(value):
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
        if parsed.tzinfo is None:
            raise ValueError("timezone required")
        return parsed.astimezone(timezone.utc).isoformat()
    except (ValueError, OverflowError) as error:
        raise argparse.ArgumentTypeError("deadline must be an absolute time with a timezone") from error


def remaining_seconds(deadline):
    remaining = timestamp(deadline) - time.time()
    if remaining <= 0:
        raise DeadlineReached(f"absolute UTC deadline reached: {deadline}")
    return remaining


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_bytes(path, content, replace=True):
    temporary = path.with_name(path.name + ".tmp-" + uuid.uuid4().hex)
    try:
        with temporary.open("xb") as output:
            output.write(content)
            output.flush()
            os.fsync(output.fileno())
        if replace:
            os.replace(temporary, path)
        else:
            # A science checkpoint is create-only, even if another writer races us.
            os.link(temporary, path)
        directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_json(path, value):
    atomic_bytes(path, (json.dumps(value, sort_keys=True) + "\n").encode())


def prime_list(low, high):
    """Independent stdlib sieve used to validate the native worker's coverage."""
    limit = math.isqrt(high)
    base = bytearray(b"\1") * (limit + 1)
    base[:2] = b"\0\0"
    for p in range(2, math.isqrt(limit) + 1):
        if base[p]:
            base[p * p::p] = b"\0" * ((limit - p * p) // p + 1)
    segment = bytearray(b"\1") * (high - low + 1)
    for p in range(2, limit + 1):
        if base[p]:
            first = max(p * p, ((low + p - 1) // p) * p)
            if first <= high:
                segment[first - low::p] = b"\0" * ((high - first) // p + 1)
    return [low + i for i, present in enumerate(segment) if present and low + i >= 2]


def validate_archive(archive, low, high, config):
    backend = config["backend"]
    if (
        archive["format"] != "lerch-range-v2"
        or archive["interval"] != {"start": low, "end": high}
        or archive["requested_batch"] != 8192
        or archive["sample"]["method"] != METHODS[backend]
        or archive["sample"]["threads"] != config["threads"]
    ):
        raise ValueError("chunk archive configuration mismatch")
    provenance = archive["provenance"]
    if (
        provenance["backend"] != backend
        or provenance["kernel"] != KERNELS[backend]
        or provenance["max_prime"] != 2_000_000_000
        or provenance["binary_sha256"] != config["binary_sha256"]
        or not isinstance(provenance["package_version"], str)
        or not provenance["package_version"]
    ):
        raise ValueError("chunk archive provenance mismatch")
    rows = archive["canonical"]
    primes = prime_list(low, high)
    if [row["p"] for row in rows] != primes:
        raise ValueError("independent sieve found missing, duplicate, or reordered primes")
    digest = hashlib.sha256(b"lerch-canonical-v1\0")
    for row in rows:
        p = row["p"]
        remainder = row["lerch_remainder"]
        if (
            any(type(row[key]) is not int for key in ("p", "primitive_root", "q1", "q2"))
            or not 0 < row["primitive_root"] < p
            or not 0 <= row["q1"] < p
            or not 0 <= row["q2"] < p
            or type(row["is_lerch"]) is not bool
            or (p == 2 and (
                row["primitive_root"] != 1 or remainder is not None
                or row["q1"] or row["q2"] or row["is_lerch"]
            ))
            or (p != 2 and (
                type(remainder) is not int or not 0 <= remainder < p
                or remainder != ((row["q1"] ** 2 + row["q2"] - 2 * row["q1"]) * ((p + 1) // 2)) % p
                or row["is_lerch"] != (remainder == 0)
            ))
        ):
            raise ValueError(f"invalid canonical tuple at {p}")
        digest.update(struct.pack("<QQQQBQB", p, row["primitive_root"], row["q1"], row["q2"],
                                  int(remainder is not None), remainder or 0, int(row["is_lerch"])))
    odd_terms = sum(p - 1 for p in primes if p != 2)
    summary = {
        "primes": len(primes),
        "first_prime": primes[0] if primes else None,
        "last_prime": primes[-1] if primes else None,
        "residue_terms": sum(p - 1 for p in primes),
        "odd_recurrence_terms": odd_terms,
        "primary_pair_steps": odd_terms // 2,
        "hits": [row for row in rows if row["is_lerch"]],
        "sha256": digest.hexdigest(),
    }
    if archive["sample"]["summary"] != summary:
        raise ValueError("chunk summary/hash does not match every canonical tuple")
    timing = archive["sample"]["timing"]
    values = [timing[key] for key in (
        "seconds", "pool_sieve_seconds", "validation_setup_traversal_collection_seconds",
        "pool_teardown_seconds",
    )]
    if any(type(value) not in (int, float) or not math.isfinite(value) or value < 0 for value in values):
        raise ValueError("invalid chunk compute time")
    if values[0] <= 0 or not math.isclose(values[0], sum(values[1:]), rel_tol=1e-9, abs_tol=1e-8):
        raise ValueError("inconsistent chunk timing components")
    return summary


def request_stop(_signum, _frame):
    global STOP_REQUESTED
    STOP_REQUESTED = True


def stopped(root):
    return STOP_REQUESTED or (root / "STOP").exists()


def save(root, manifest):
    manifest["updated_utc"] = utc()
    atomic_json(root / "manifest.json", manifest)
    atomic_bytes(root / "progress.jsonl", "".join(
        json.dumps(row, sort_keys=True) + "\n" for row in manifest["progress"]
    ).encode())


def pause_at_deadline(root, manifest, error):
    if manifest["status"] == "complete":
        print("deadline expired; campaign was already complete", flush=True)
        return
    manifest["status"] = "deadline_reached"
    manifest.pop("error", None)
    manifest.setdefault("deadline_reached_utc", utc())
    manifest["deadline_reason"] = str(error)
    save(root, manifest)
    print(f"deadline_reached: next_start={manifest['next_start']}; "
          f"unfinished output retained, not checkpointed; {error}", flush=True)


def checked_process(command, log, deadline):
    with log.open("xb") as output:
        try:
            result = subprocess.run(
                command, stdout=output, stderr=subprocess.STDOUT, check=False,
                timeout=remaining_seconds(deadline),
            )
        except subprocess.TimeoutExpired as error:
            # run() kills and reaps the worker process, including its native threads.
            raise DeadlineReached(f"binary timed out at {deadline}; inspect {log.name}") from error
        finally:
            output.flush()
            os.fsync(output.fileno())
    if result.returncode:
        raise RuntimeError(f"command exited {result.returncode}; inspect {log.name}")


def chunk_end(config, low):
    block_end = min(
        config["start"] + ((low - config["start"]) // config["progress_width"] + 1)
        * config["progress_width"] - 1, config["end"],
    )
    return min(low + config["chunk_width"] - 1, block_end), block_end


def archive_name(low, high):
    return f"chunks/chunk-{low:010d}-{high:010d}.json.gz"


def candidates(summary, low, high):
    return [{
        **hit, "chunk_start": low, "chunk_end": high,
        "verification_status": "pending-independent-verification",
    } for hit in summary["hits"]]


def progress_totals(records, start, end):
    return {
        "start": start, "end": end,
        "primes": sum(row["summary"]["primes"] for row in records),
        "residue_terms": sum(row["summary"]["residue_terms"] for row in records),
        "candidate_hits": [hit["p"] for row in records for hit in row["summary"]["hits"]],
    }


def validate_completed(root, manifest):
    config = manifest["config"]
    cursor, primes, terms, hits, expected_progress, block = config["start"], 0, 0, [], [], []
    for record in manifest["chunks"]:
        remaining_seconds(config["deadline_utc"])
        high, block_end = chunk_end(config, cursor)
        if record["start"] != cursor or record["end"] != high or high > config["end"]:
            raise ValueError("checkpoint intervals are not contiguous configured chunks")
        if record["archive"] != archive_name(cursor, high):
            raise ValueError("checkpoint archive name does not match its interval")
        path = root / record["archive"]
        if sha256_file(path) != record["archive_sha256"]:
            raise ValueError("compressed checkpoint checksum mismatch")
        archive = json.loads(gzip.decompress(path.read_bytes()))
        summary = validate_archive(archive, cursor, high, config)
        if (record["summary"] != summary or record["archive_bytes"] != path.stat().st_size
                or record["compute_seconds"] != archive["sample"]["timing"]["seconds"]):
            raise ValueError("checkpoint record does not match its archive")
        primes += summary["primes"]
        terms += summary["residue_terms"]
        hits.extend(candidates(summary, cursor, high))
        block.append(record)
        if high == block_end:
            expected_progress.append(progress_totals(block, block[0]["start"], high))
            block = []
        cursor = high + 1
    if (manifest["next_start"] != cursor or manifest["primes"] != primes
            or manifest["residue_terms"] != terms or manifest["candidate_hits"] != hits):
        raise ValueError("checkpoint cursor/totals/candidate verification state mismatch")
    if len(manifest["progress"]) != len(expected_progress) or any(
        any(record[key] != value for key, value in expected.items())
        for record, expected in zip(manifest["progress"], expected_progress)
    ):
        raise ValueError("checkpoint progress totals mismatch")
    if manifest["status"] == "complete" and cursor != config["end"] + 1:
        raise ValueError("incomplete coverage labeled complete")


def campaign(args, root, manifest):
    config = manifest["config"]
    deadline = config["deadline_utc"]
    remaining_seconds(deadline)
    if manifest["next_start"] > config["end"]:
        if manifest["status"] != "complete":
            manifest["status"] = "complete"
            save(root, manifest)
        print("campaign already complete", flush=True)
        return
    if not manifest["search_started_utc"]:
        manifest["search_started_utc"] = utc()
    manifest.pop("error", None)
    completed_this_run = 0
    while manifest["next_start"] <= config["end"]:
        remaining_seconds(deadline)
        if stopped(root) or (args.max_chunks is not None and completed_this_run >= args.max_chunks):
            manifest["status"] = "paused"
            save(root, manifest)
            print(f"paused at next_start={manifest['next_start']}", flush=True)
            return
        if shutil.disk_usage(root).free < args.min_free_bytes:
            raise RuntimeError("free disk below safety reserve; refusing a new chunk")
        manifest["status"] = "running"
        save(root, manifest)
        low = manifest["next_start"]
        high, block_end = chunk_end(config, low)
        compressed = root / archive_name(low, high)
        wall_start = time.monotonic()
        recovered = compressed.exists()
        if recovered:
            raw = gzip.decompress(compressed.read_bytes())
        else:
            partial = compressed.with_name(f"chunk-{low:010d}-{high:010d}.partial-{uuid.uuid4().hex}.json")
            checked_process([
                str(args.binary), "benchmark", "--start", str(low), "--end", str(high),
                "--backend", config["backend"], "--output", str(partial),
                "--threads", str(config["threads"]),
            ], partial.with_suffix(".log"), deadline)
            raw = partial.read_bytes()
        archive = json.loads(raw)
        summary = validate_archive(archive, low, high, config)
        remaining_seconds(deadline)
        if not recovered:
            packed = gzip.compress(raw, mtime=0)
            if gzip.decompress(packed) != raw:
                raise RuntimeError("gzip round-trip mismatch")
            atomic_bytes(compressed, packed, replace=False)
        record = {
            "start": low, "end": high, "completed_utc": utc(),
            "archive": archive_name(low, high), "archive_sha256": sha256_file(compressed),
            "archive_bytes": compressed.stat().st_size, "summary": summary,
            "compute_seconds": archive["sample"]["timing"]["seconds"],
            "wall_seconds": time.monotonic() - wall_start,
            "recovered_after_interruption": recovered,
        }
        manifest["chunks"].append(record)
        manifest["next_start"] = high + 1
        manifest["primes"] += summary["primes"]
        manifest["residue_terms"] += summary["residue_terms"]
        manifest["candidate_hits"].extend(candidates(summary, low, high))
        if high == block_end:
            block_start = config["start"] + ((low - config["start"]) // config["progress_width"]) * config["progress_width"]
            records = [row for row in manifest["chunks"] if row["start"] >= block_start]
            progress = progress_totals(records, block_start, high)
            progress.update({"utc": utc(), "cumulative_wall_seconds":
                             time.time() - timestamp(manifest["search_started_utc"])})
            manifest["progress"].append(progress)
            print(json.dumps({"progress": progress}), flush=True)
        if high == config["end"]:
            manifest["status"] = "complete"
        save(root, manifest)
        if not recovered:
            partial.unlink()
        completed_this_run += 1
        print(f"checkpoint [{low},{high}] backend={config['backend']} primes={summary['primes']} "
              f"seconds={record['compute_seconds']:.6f} next={high + 1} "
              f"candidates={[hit['p'] for hit in summary['hits']]}", flush=True)
    print(f"completed inclusive search [{config['start']},{config['end']}]; "
          "candidate hits require separate independent verification", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--start", type=int, default=2)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--chunk-size", type=int, default=100_000)
    parser.add_argument("--progress-width", type=int, default=10_000_000)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--backend", choices=METHODS, required=True)
    parser.add_argument("--source-sha", help="optional source revision recorded without inference")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--deadline-utc", type=absolute_utc,
                        help="absolute deadline, at most 24h after creation; inherited, never reset on resume")
    parser.add_argument("--max-chunks", type=int)
    parser.add_argument("--min-free-bytes", type=int, default=1_073_741_824)
    parser.add_argument("--runner-sha256", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if (
        not 2 <= args.start <= args.end <= 2_000_000_000
        or not 1 <= args.chunk_size <= 100_000 or not 1 <= args.progress_width
        or not 1 <= args.threads <= 256 or args.min_free_bytes < 0
        or (args.max_chunks is not None and args.max_chunks < 1)
    ):
        parser.error("invalid bounded search configuration")
    root = args.output_dir.resolve()
    path = root / "manifest.json"
    previous = None
    if args.resume:
        if not path.is_file():
            parser.error("--resume requires an existing v2 manifest")
        previous = json.loads(path.read_text())
        if previous["format"] != FORMAT:
            parser.error("archived v1 data is read-only; audit it separately, do not resume or rewrite it")
    elif root.exists() and any(root.iterdir()):
        parser.error("output directory is not empty; use --resume for an existing v2 search")
    created = previous["created_utc"] if previous else utc()
    if args.deadline_utc is None:
        if previous:
            args.deadline_utc = previous["config"]["deadline_utc"]
        elif args.end > 1_000_000_000:
            parser.error("searches above 1000000000 require --deadline-utc")
        else:
            args.deadline_utc = (datetime.fromisoformat(created) + timedelta(hours=24)).isoformat()
    if timestamp(args.deadline_utc) > timestamp(created) + 86_400.000001:
        parser.error("deadline may not exceed the original 24-hour budget")
    args.binary = args.binary.resolve(strict=True)
    runner_sha = args.runner_sha256 or sha256_file(Path(__file__).resolve())
    if len(runner_sha) != 64 or any(char not in "0123456789abcdef" for char in runner_sha):
        parser.error("invalid runner SHA-256")
    config = {
        "start": args.start, "end": args.end, "chunk_width": args.chunk_size,
        "progress_width": args.progress_width, "threads": args.threads,
        "backend": args.backend, "method": METHODS[args.backend],
        "source_sha": args.source_sha, "binary_sha256": sha256_file(args.binary),
        "runner_sha256": runner_sha, "deadline_utc": args.deadline_utc,
    }
    if previous is not None and previous["config"] != config:
        raise ValueError("refusing to resume different source, binary, backend, deadline, or configuration")
    root.mkdir(parents=True, exist_ok=True)
    with (root / "search.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if args.resume:
            manifest = json.loads(path.read_text())
            if manifest != previous:
                raise ValueError("manifest changed while acquiring the search lock; retry")
        else:
            if any(entry.name != "search.lock" for entry in root.iterdir()):
                raise ValueError("output directory changed while acquiring the search lock")
            manifest = {
                "format": FORMAT, "config": config, "created_utc": created,
                "status": "initializing", "method": config["method"],
                "next_start": args.start, "search_started_utc": None,
                "primes": 0, "residue_terms": 0, "candidate_hits": [],
                "chunks": [], "progress": [],
            }
        lock.seek(0)
        lock.truncate()
        lock.write(str(os.getpid()) + "\n")
        lock.flush()
        (root / "chunks").mkdir(exist_ok=True)
        signal.signal(signal.SIGTERM, request_stop)
        signal.signal(signal.SIGINT, request_stop)
        try:
            remaining_seconds(config["deadline_utc"])
            validate_completed(root, manifest)
        except DeadlineReached as error:
            pause_at_deadline(root, manifest, error)
            return
        try:
            campaign(args, root, manifest)
        except DeadlineReached as error:
            pause_at_deadline(root, manifest, error)
        except (OSError, ValueError, KeyError, TypeError, RuntimeError, subprocess.SubprocessError) as error:
            manifest["status"] = "paused" if STOP_REQUESTED else "failed"
            manifest["error"] = f"{type(error).__name__}: {error}"
            save(root, manifest)
            print(f"{manifest['status']}: {manifest['error']}", file=sys.stderr, flush=True)
            if not STOP_REQUESTED:
                raise


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, TypeError, RuntimeError, subprocess.SubprocessError) as error:
        raise SystemExit(f"search failed: {error}") from error
