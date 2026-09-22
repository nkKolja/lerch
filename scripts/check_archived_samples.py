#!/usr/bin/env python3
"""Nine fixed, stratified original-generic checks against the frozen 1B data."""

from concurrent.futures import ThreadPoolExecutor
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess

from audit_results import read_json, regular_file, require, same, sha256_file


def check_samples(root: Path, binary: Path, output: Path) -> dict:
    manifest = read_json(regular_file(root, "manifest.json").read_bytes())
    require(
        manifest["format"] == "lerch-campaign-v1" and manifest["status"] == "complete"
        and (manifest["config"]["start"], manifest["config"]["end"]) == (200_000_000, 1_000_000_000)
        and len(manifest["chunks"]) == 8001,
        "this fixed sample protocol requires the completed 200M-1B campaign",
    )
    require(sha256_file(binary) == manifest["config"]["binary_sha256"], "wrong historical binary")
    require(not output.resolve().is_relative_to(root.resolve()), "output must be outside input")
    output.mkdir()
    selections = [(index, 0) for index in range(0, 8000, 1000)] + [(7999, -1)]
    expected = []
    for index, row_index in selections:
        record = manifest["chunks"][index]
        packed = regular_file(root, record["archive"]).read_bytes()
        require(hashlib.sha256(packed).hexdigest() == record["archive_sha256"],
                "sample archive checksum mismatch")
        archive = read_json(gzip.decompress(packed))
        expected.append((record["archive"], archive["canonical"][row_index]))

    def check(item: tuple[str, dict]) -> dict:
        archive_name, row = item
        p = row["p"]
        result_path = output / f"prime-{p}.json"
        command = [str(binary), "range", str(p), str(p), "generic", str(result_path), "1"]
        with (output / f"prime-{p}.log").open("x") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=180)
        result = read_json(result_path.read_bytes())
        same(result["canonical"], [row], f"full generic tuple at {p}")
        return {
            "p": p, "campaign_archive": archive_name, "full_canonical_tuple_matched": True,
            "canonical": row, "generic_seconds": result["sample"]["timing"]["seconds"],
            "raw_result": result_path.name, "raw_sha256": sha256_file(result_path),
        }

    with ThreadPoolExecutor(max_workers=8) as executor:
        samples = list(executor.map(check, expected))
    report = {
        "format": "lerch-stratified-generic-checks-v1",
        "selection": "first prime in chunks 0,1000,...,7000; last prime in chunk 7999",
        "manifest_sha256": sha256_file(root / "manifest.json"),
        "binary_sha256": manifest["config"]["binary_sha256"],
        "source_sha": manifest["config"]["source_sha"],
        "script_sha256": sha256_file(Path(__file__)),
        "samples": samples,
        "scope": "Separate original generic recurrence, same historical binary and shared arithmetic; not definition-level or a full campaign rerun. Concurrent single-worker timings are not benchmark samples.",
    }
    with (output / "summary.json").open("x") as target:
        json.dump(report, target, indent=2)
        target.write("\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        print(json.dumps(check_samples(args.directory, args.binary.resolve(), args.output_dir), indent=2))
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
        parser.exit(1, f"sample check failed: {error}\n")


if __name__ == "__main__":
    main()
