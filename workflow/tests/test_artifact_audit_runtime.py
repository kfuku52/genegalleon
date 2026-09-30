"""Full-audit consistency, progress protocol and bounded reuse regressions."""
import concurrent.futures
import hashlib
import importlib
import json
import os
import subprocess
import sys
import threading
import time
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
provenance = importlib.import_module("artifact_provenance")
runtime = importlib.import_module("artifact_audit_runtime")
store_module = importlib.import_module("gene_family_output_store")
AuditDigests, AuditProgress = runtime.AuditDigests, runtime.AuditProgress
source_identity, validate_observation = runtime.source_identity, runtime.validate_observation
ArchiveStoreError, GeneFamilyOutputStore = store_module.ArchiveStoreError, store_module.GeneFamilyOutputStore
INDEX_EPOCH_FILE, INDEX_UPDATE_FILE = store_module.INDEX_EPOCH_FILE, store_module.INDEX_UPDATE_FILE


def fixture(root, families=8):
    logical = root / "output/orthogroup"
    manifests = logical / "artifact_provenance"
    manifests.mkdir(parents=True)
    shared = root / "shared.tsv"
    shared.write_bytes(b"shared input\n" * 1000)
    for index in range(families):
        family = f"OG{index:07d}"
        payload = {"schema_version": 1, "family_id": family, "step": "example", "inputs": [{
            "label": "shared", "scope": "workspace", "path": "shared.tsv", "artifact_type": "file",
            "sha256": hashlib.sha256(shared.read_bytes()).hexdigest(), "size_bytes": shared.stat().st_size}],
            "outputs": []}
        (manifests / (family + ".example.json")).write_text(json.dumps(payload))
    return logical


def arguments(root, logical, report, workers=4):
    return provenance.build_parser().parse_args(["audit", "--workspace-root", str(root),
        "--logical-root", str(logical), "--mode", "orthogroup", "--output-tsv", str(report),
        "--workers", str(workers), "--progress-interval", "0"])


def test_digest_reuse_single_flight_and_final_change_fence(tmp_path):
    path = tmp_path / "shared"
    path.write_bytes(b"abc")
    memo, calls, lock = AuditDigests(), [], threading.Lock()
    def compute():
        with lock:
            calls.append(1)
        time.sleep(0.01)
        data = path.read_bytes()
        return hashlib.sha256(data).hexdigest(), len(data)
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        rows = list(pool.map(lambda _: memo.read(path, compute), range(40)))
    assert len(set(rows)) == 1 and len(calls) == 1
    assert memo.metrics()["digest_bytes_read"] == 3
    memo.validate()
    previous = path.stat()
    path.write_bytes(b"xyz")
    os.utime(path, ns=(previous.st_atime_ns, previous.st_mtime_ns))
    with pytest.raises(ValueError, match="changed before publication"):
        memo.validate()


def test_directory_membership_and_distinct_member_count_constraints(tmp_path):
    path = tmp_path / "directory"
    path.mkdir()
    child = path / "file"
    child.write_text("abc")
    memo, calls = AuditDigests(), []
    def compute():
        calls.append(1)
        return hashlib.sha256(child.read_bytes()).hexdigest(), child.stat().st_size
    memo.read(path, compute, directory=True, constraint=1)
    memo.read(path, compute, directory=True, constraint=2)
    assert len(calls) == 2
    child.write_text("xyz")
    with pytest.raises(ValueError, match="changed before publication"):
        memo.validate()


def test_coarse_filesystem_times_require_content_revalidation(tmp_path, monkeypatch):
    path = tmp_path / "source"
    path.write_bytes(b"abc")
    original = runtime.signature
    def coarse(path):
        values = original(path)
        return values[:4] + tuple(value // 1_000_000_000 * 1_000_000_000 for value in values[4:])
    monkeypatch.setattr(runtime, "signature", coarse)
    memo = AuditDigests()
    def compute():
        return hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_size
    memo.read(path, compute)
    memo.validate()
    assert memo.metrics()["digest_revalidation_reads"] == 1
    path.write_bytes(b"xyz")
    with pytest.raises(ValueError, match="changed before publication"):
        memo.validate()


def test_guard_capacity_fails_closed_instead_of_dropping_source_checks():
    memo = AuditDigests(limit=1)
    memo.guard("one", 1, lambda: 1)
    with pytest.raises(ValueError, match="bounded capacity"):
        memo.guard("two", 2, lambda: 2)


def test_digest_capacity_fails_closed_even_for_concurrent_first_reads(tmp_path):
    memo = AuditDigests(limit=1)
    paths = []
    for index in range(256):
        path = tmp_path / str(index)
        key = str(path.absolute()), None, False, None
        if not paths or hash(key) % 256 != hash((str(paths[0].absolute()), None, False, None)) % 256:
            path.write_bytes(b"abc")
            paths.append(path)
        if len(paths) == 2:
            break
    barrier = threading.Barrier(2)
    def read(path):
        first = True
        def compute():
            nonlocal first
            if first:
                first = False
                barrier.wait(timeout=5)
            return hashlib.sha256(path.read_bytes()).hexdigest(), 3
        return memo.read(path, compute)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(read, path) for path in paths]
        outcomes = []
        for future in futures:
            try:
                outcomes.append(future.result())
            except ValueError as exc:
                outcomes.append(str(exc))
    assert sum(isinstance(value, str) and "bounded capacity" in value for value in outcomes) == 1
    assert len(memo.entries) == 1
    memo.validate()


def test_workspace_mutation_between_manifest_reads_cannot_rebind_digest(tmp_path, monkeypatch):
    logical = fixture(tmp_path, families=2)
    original = provenance.audit_manifest
    def changing(payload, *args):
        result = original(payload, *args)
        if payload["family_id"] == "OG0000000":
            shared = tmp_path / "shared.tsv"
            shared.write_bytes(shared.read_bytes() + b"changed\n")
            second = logical / "artifact_provenance/OG0000001.example.json"
            updated = json.loads(second.read_text())
            updated["inputs"][0]["sha256"] = hashlib.sha256(shared.read_bytes()).hexdigest()
            updated["inputs"][0]["size_bytes"] = shared.stat().st_size
            second.write_text(json.dumps(updated))
        return result
    monkeypatch.setattr(provenance, "audit_manifest", changing)
    report = tmp_path / "audit.tsv"
    assert provenance.audit(arguments(tmp_path, logical, report, workers=1)) == 1
    assert "snapshot_identity\taudit_error" in report.read_text()


@pytest.mark.parametrize("marker", [INDEX_EPOCH_FILE, INDEX_UPDATE_FILE])
def test_store_snapshot_rejects_metadata_publication_or_interrupted_update(tmp_path, marker):
    store = GeneFamilyOutputStore(tmp_path)
    store.archive_root.mkdir()
    with pytest.raises(ArchiveStoreError):
        with store.read_snapshot():
            (store.archive_root / marker).write_text("changed")
    assert store._read_snapshot_epoch is None
    assert not store._snapshot_zip_readers


def test_optional_absence_is_fenced(tmp_path):
    logical = fixture(tmp_path)
    memo = AuditDigests()
    entry = {"scope": "workspace", "path": "optional", "state": "absent", "artifact_type": "file"}
    with pytest.raises(FileNotFoundError):
        provenance.audit_entry_digest(entry, GeneFamilyOutputStore(logical), logical, tmp_path, memo)
    (tmp_path / "optional").write_text("new output")
    with pytest.raises(ValueError, match="Logical audit source changed"):
        memo.validate()


def test_threaded_audit_matches_serial_and_binds_machine_result(tmp_path):
    logical = fixture(tmp_path)
    reports = [tmp_path / f"audit-{workers}.tsv" for workers in (1, 4)]
    for workers, report in zip((1, 4), reports, strict=True):
        assert provenance.audit(arguments(tmp_path, logical, report, workers)) == 0
        result = json.loads(report.with_suffix(".tsv.result.json").read_text())
        assert result["checked"] == 8 and result["status_counts"] == {"current": 8}
        assert result["source_sha256"] == source_identity(SUPPORT)
        assert result["report_sha256"] == hashlib.sha256(report.read_bytes()).hexdigest()
        assert result["metrics"]["unique_digest_reads"] == 1
    assert reports[0].read_bytes() == reports[1].read_bytes()


def test_manifest_mutation_during_audit_fails_before_success_publication(tmp_path, monkeypatch):
    logical = fixture(tmp_path, 1)
    original = provenance.audit_manifest
    def changing(*args):
        result = original(*args)
        next((logical / "artifact_provenance").iterdir()).write_text("{}")
        return result
    monkeypatch.setattr(provenance, "audit_manifest", changing)
    report = tmp_path / "audit.tsv"
    assert provenance.audit(arguments(tmp_path, logical, report, 1)) == 1
    assert "snapshot_identity\taudit_error" in report.read_text()
    assert json.loads(report.with_suffix(".tsv.result.json").read_text())["exit_code"] == 1


def test_full_audit_does_not_reuse_persistent_cross_run_digests(tmp_path, monkeypatch):
    logical = fixture(tmp_path)
    provenance.configure_digest_cache(tmp_path / "old-cache.sqlite3")
    previous = provenance.digest_cache().database
    original = provenance.audit_manifest
    def inspect(*args):
        assert provenance.digest_cache() is None
        return original(*args)
    monkeypatch.setattr(provenance, "audit_manifest", inspect)
    try:
        assert provenance.audit(arguments(tmp_path, logical, tmp_path / "audit.tsv")) == 0
        assert provenance.digest_cache().database == previous
    finally:
        provenance.configure_digest_cache(None)


@pytest.mark.parametrize("interval", [-1, float("nan"), float("inf")])
def test_progress_rejects_invalid_interval(tmp_path, interval):
    with pytest.raises(ValueError):
        AuditProgress(tmp_path / "progress.json", interval=interval)


def test_server_stdout_remains_numeric_with_progress(tmp_path):
    logical = fixture(tmp_path)
    argv = ["audit", "--workspace-root", str(tmp_path), "--logical-root", str(logical),
            "--mode", "orthogroup", "--output-tsv", str(tmp_path / "audit.tsv"), "--progress-interval", "0.01"]
    result = subprocess.run([sys.executable, "-B", str(SUPPORT / "artifact_provenance.py"), "serve"],
                            input=("\0".join(argv) + "\0\0").encode(), capture_output=True)
    assert result.returncode == 0 and result.stdout == b"0\n", result.stderr
    progress = json.loads((tmp_path / "audit.tsv.progress.json").read_text())
    assert progress["state"] == "completed" and progress["completed"] == progress["total"] == 8


def test_progress_validation_rejects_false_counts_and_attempt(tmp_path):
    attempt = "a" * 32
    with AuditProgress(tmp_path / "progress.json", attempt_dir=tmp_path / attempt,
                       source_sha256="b" * 64) as progress:
        row = progress.payload()
        validate_observation(row, attempt)
        row["completed"] = True
        with pytest.raises(ValueError):
            validate_observation(row, attempt)
