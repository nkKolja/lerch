#!/usr/bin/env python3
"""Archive exact historical Rust packages and the 1B campaign's pinned runner.

Requires the selected historical Git objects, not a checkout/reset of a commit.
No evidence trees, operational configuration, Git history or credentials are
included. Every original file is byte-for-byte unchanged.
"""

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile

from audit_results import require, sha256_file


SOURCE = "616b5e3f406aba620e9d6780e93275619490c7da"
RUNNER = "ae192d18e877d1796503af3b114765ff5e3ed1a8"
RUNNER_SHA256 = "a93f7794d27ff23996f5a0ee52e71d5bd865c18dc8c92c57ee9dfad3f8d44eb1"
REVISIONS = [
    SOURCE,
    "dc4079671dc81f013c9a543f24b5e8cf63f6e6f0",
    "41236ecc9a3dc6dc7baf2168ce33d9b17c5ad10f",
    "8eaaab04e2402d70509fcc1938dd132c3e12fbef",
]
SOURCE_PATHS = [
    "Cargo.toml", "Cargo.lock", "src", "examples", "tests", "tools",
    "LICENSE", "CITATION.cff", "scripts/independent_verify.py",
    "scripts/audit_search_results.py", "scripts/audit_bernoulli_output.py",
    "scripts/build_bernoulli_verifier.sh", "scripts/test_bernoulli_verifier.sh",
]


def package_source(output: Path, revision: str = SOURCE) -> dict:
    require(revision in REVISIONS, "unknown historical source revision")
    prefix = f"lerch-producing-source-{revision[:7]}"
    packed = subprocess.check_output(["git", "archive", "--format=tar", revision, *SOURCE_PATHS])
    files = {}
    with tarfile.open(fileobj=io.BytesIO(packed)) as archive:
        for member in archive:
            if member.isdir():
                continue
            require(member.isfile() and not Path(member.name).is_absolute()
                    and ".." not in Path(member.name).parts, "unexpected source member")
            require(member.name not in files, "duplicate source member")
            files[member.name] = archive.extractfile(member).read()
    if revision == SOURCE:
        runner = subprocess.check_output(["git", "show", f"{RUNNER}:scripts/epyc_campaign.py"])
        require(hashlib.sha256(runner).hexdigest() == RUNNER_SHA256, "historical runner hash mismatch")
        files["controller/epyc_campaign.py"] = runner
    provenance = {
        "format": "lerch-producing-source-v1",
        "source_commit": revision,
        "controller_commit": RUNNER if revision == SOURCE else None,
        "source_path_selection": SOURCE_PATHS,
        "source_files_unmodified": True,
        "normalization": "sorted USTAR members; mtime/uid/gid=0; mode=0644; gzip level=9, mtime=0",
        "scope": (
            "Complete Rust package source, lockfile, examples and Rust tests; exact controller added from its later commit. Not the latest publication implementation."
            if revision == SOURCE else
            "Complete historical Rust package source, lockfile, examples and Rust tests. Not the latest publication implementation."
        ),
        "files": [{
            "path": name, "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(),
            "commit": RUNNER if name.startswith("controller/") else revision,
        } for name, raw in sorted(files.items())],
    }
    files["SOURCE_PROVENANCE.json"] = (json.dumps(provenance, sort_keys=True, indent=2) + "\n").encode()
    sums = "".join(f"{hashlib.sha256(raw).hexdigest()}  {name}\n" for name, raw in sorted(files.items()))
    files["SHA256SUMS"] = sums.encode()
    with output.open("xb") as destination:
        with gzip.GzipFile(filename="", fileobj=destination, mode="wb", compresslevel=9, mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w", format=tarfile.USTAR_FORMAT) as archive:
                for name, raw in sorted(files.items()):
                    member = tarfile.TarInfo(f"{prefix}/{name}")
                    member.size = len(raw)
                    member.mode = 0o644
                    member.uid = member.gid = member.mtime = 0
                    archive.addfile(member, io.BytesIO(raw))
    return {
        "format": "lerch-source-release-asset-v1", "name": output.name,
        "bytes": output.stat().st_size, "sha256": sha256_file(output),
        "members": len(files), "source_commit": revision,
        "controller_commit": RUNNER if revision == SOURCE else None,
        "cargo_lock_sha256": hashlib.sha256(files["Cargo.lock"]).hexdigest(),
        "license_sha256": hashlib.sha256(files["LICENSE"]).hexdigest(),
        "controller_sha256": RUNNER_SHA256 if revision == SOURCE else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--revision", choices=REVISIONS, default=SOURCE)
    args = parser.parse_args()
    try:
        print(json.dumps(package_source(args.output, args.revision), indent=2))
    except (OSError, ValueError, KeyError, tarfile.TarError, subprocess.SubprocessError) as error:
        parser.exit(1, f"source packaging failed: {error}\n")


if __name__ == "__main__":
    main()
