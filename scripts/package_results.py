#!/usr/bin/env python3
"""Build a deterministic allowlisted tar without changing the original campaign."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile

from audit_results import read_json, regular_file, require, sha256_file, validate_progress


PREFIX = "lerch-campaign-v1"


def package(root: Path, audit_path: Path, license_path: Path, output: Path) -> dict:
    root = root.resolve(strict=True)
    require(not output.resolve().is_relative_to(root), "archive must be outside input")
    manifest_raw = regular_file(root, "manifest.json").read_bytes()
    manifest = read_json(manifest_raw)
    manifest_sha = hashlib.sha256(manifest_raw).hexdigest()
    report_raw = audit_path.read_bytes()
    report = read_json(report_raw)
    audit_script = Path(__file__).with_name("audit_results.py")
    require(
        report["format"] == "lerch-evidence-audit-v1" and report["audit"] == "passed"
        and report["manifest_digest_anchored"] is True
        and report["manifest_sha256"] == manifest_sha
        and report["auditor_sha256"] == sha256_file(audit_script)
        and report["independent_prime_regeneration"] is True
        and report["moment_remainder_identity_checked"] is True,
        "a matching, anchored full audit is required",
    )
    require(manifest["format"] == "lerch-campaign-v1" and manifest["status"] == "complete",
            "not a completed v1 campaign")
    require(report["chunks"] == len(manifest["chunks"]), "audit chunk count mismatch")
    config = manifest["config"]
    validate_progress(root, manifest)
    paths = ["manifest.json", "progress.jsonl", "progress.log"]
    for record in manifest["chunks"]:
        require(record["archive"] == (
            f"chunks/chunk-{record['start']:010d}-{record['end']:010d}.json.gz"
        ), "noncanonical archive member name")
        paths.append(record["archive"])
    require(len(set(paths)) == len(paths), "duplicate archive member")
    expected = {r["archive"]: (r["archive_sha256"], r["archive_bytes"]) for r in manifest["chunks"]}
    bundle = {
        "format": "lerch-campaign-bundle-v1",
        "original_manifest_sha256": manifest_sha,
        "canonical_archives_unmodified": True,
        "source_sha": config["source_sha"],
        "binary_sha256": config["binary_sha256"],
        "runner_sha256": config["runner_sha256"],
        "original_files": len(paths),
        "packaging": "uncompressed USTAR; fixed member order; uid/gid/mtime=0; mode=0644",
        "excluded": [
            "unlisted/orphan files", "partial checkpoints", "service/SSH logs",
            "locks", "morning partial snapshot", "deployment configuration",
        ],
        "audit_scope": "hashes, independent sieve, row bounds/identity, summaries/progress; no full Fermat-quotient recomputation",
    }
    extras = {
        "provenance/audit.json": report_raw,
        "provenance/audit_results.py": audit_script.read_bytes(),
        "provenance/package_results.py": Path(__file__).read_bytes(),
        "provenance/LICENSE": license_path.read_bytes(),
        "provenance/bundle.json": (json.dumps(bundle, sort_keys=True, indent=2) + "\n").encode(),
    }
    checksums = []

    def add(archive: tarfile.TarFile, relative: str, raw: bytes) -> None:
        digest = hashlib.sha256(raw).hexdigest()
        if relative in expected:
            require((digest, len(raw)) == expected[relative], f"changed archive: {relative}")
        info = tarfile.TarInfo(f"{PREFIX}/{relative}")
        info.size = len(raw)
        info.mode = 0o644
        info.uid = info.gid = info.mtime = 0
        archive.addfile(info, io.BytesIO(raw))
        checksums.append(f"{digest}  {relative}\n")

    with output.open("xb") as destination:
        try:
            with tarfile.open(fileobj=destination, mode="w", format=tarfile.USTAR_FORMAT) as archive:
                for relative in paths:
                    add(archive, relative, regular_file(root, relative).read_bytes())
                for relative, raw in extras.items():
                    add(archive, relative, raw)
                add(archive, "SHA256SUMS", "".join(checksums).encode())
            require((root / "manifest.json").read_bytes() == manifest_raw,
                    "manifest changed during packaging")
        except (OSError, ValueError, KeyError, TypeError, tarfile.TarError, KeyboardInterrupt):
            output.unlink()
            raise
    return {
        "format": "lerch-release-asset-v1",
        "name": output.name,
        "bytes": output.stat().st_size,
        "sha256": sha256_file(output),
        "members": len(paths) + len(extras) + 1,
        "original_manifest_sha256": manifest_sha,
        "original_chunk_archives": len(manifest["chunks"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--audit-report", type=Path, required=True)
    parser.add_argument("--license", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        print(json.dumps(package(args.directory, args.audit_report, args.license, args.output), indent=2))
    except (OSError, ValueError, KeyError, TypeError, tarfile.TarError) as error:
        parser.exit(1, f"packaging failed: {error}\n")


if __name__ == "__main__":
    main()
