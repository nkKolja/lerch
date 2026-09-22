#!/usr/bin/env python3

import hashlib
import json
from pathlib import Path
import tarfile
import tempfile
import unittest

import audit_results
import package_results
from test_audit_results import fixture


class PackageTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.work = Path(self.temporary.name)
        self.root = self.work / "campaign"
        self.root.mkdir()
        self.manifest = fixture(self.root)
        self.report = audit_results.audit(
            self.root, expected_sha256=audit_results.sha256_file(self.root / "manifest.json"),
        )
        self.report_path = self.work / "audit.json"
        self.report_path.write_text(json.dumps(self.report) + "\n")
        self.license = Path(__file__).resolve().parent.parent / "LICENSE"
        self.output = self.work / "results.tar"

    def package(self, output: Path | None = None) -> dict:
        return package_results.package(
            self.root, self.report_path, self.license, output or self.output,
        )

    def test_deterministic_exact_allowlist_and_checksums(self) -> None:
        for name in ("service.log", "campaign.lock", "secret.env", "chunks/orphan.json.gz"):
            (self.root / name).write_text("MUST NOT BE PACKAGED")
        first = self.package()
        second_path = self.work / "again.tar"
        second = self.package(second_path)
        self.assertEqual(first["sha256"], second["sha256"])
        with tarfile.open(self.output) as archive:
            members = archive.getmembers()
            self.assertEqual(len(members), first["members"])
            for member in members:
                self.assertTrue(member.isfile())
                self.assertTrue(member.name.startswith("lerch-campaign-v1/"))
                self.assertNotIn("..", Path(member.name).parts)
                self.assertEqual((member.uid, member.gid, member.mtime, member.mode), (0, 0, 0, 0o644))
                self.assertNotIn(b"MUST NOT BE PACKAGED", archive.extractfile(member).read())
            sums = archive.extractfile("lerch-campaign-v1/SHA256SUMS").read().decode()
            for line in sums.splitlines():
                digest, relative = line.split("  ", 1)
                raw = archive.extractfile(f"lerch-campaign-v1/{relative}").read()
                self.assertEqual(hashlib.sha256(raw).hexdigest(), digest)
            for record in self.manifest["chunks"]:
                raw = archive.extractfile(f"lerch-campaign-v1/{record['archive']}").read()
                self.assertEqual(raw, (self.root / record["archive"]).read_bytes())

    def test_rejects_mutation_after_audit_and_removes_partial_package(self) -> None:
        (self.root / self.manifest["chunks"][0]["archive"]).write_bytes(b"changed")
        with self.assertRaisesRegex(audit_results.AuditError, "changed archive"):
            self.package()
        self.assertFalse(self.output.exists())

    def test_rejects_changed_manifest(self) -> None:
        path = self.root / "manifest.json"
        path.write_bytes(path.read_bytes() + b"\n")
        with self.assertRaisesRegex(audit_results.AuditError, "matching, anchored full audit"):
            self.package()
        self.assertFalse(self.output.exists())

    def test_rejects_unanchored_or_different_auditor(self) -> None:
        for key, value in (("manifest_digest_anchored", False), ("auditor_sha256", "0" * 64)):
            with self.subTest(key=key):
                altered = {**self.report, key: value}
                self.report_path.write_text(json.dumps(altered))
                with self.assertRaisesRegex(audit_results.AuditError, "matching, anchored full audit"):
                    self.package()

    def test_rejects_output_inside_input(self) -> None:
        with self.assertRaisesRegex(audit_results.AuditError, "outside input"):
            self.package(self.root / "output.tar")

    def test_does_not_overwrite(self) -> None:
        self.output.write_bytes(b"existing")
        with self.assertRaises(FileExistsError):
            self.package()
        self.assertEqual(self.output.read_bytes(), b"existing")


if __name__ == "__main__":
    unittest.main()
