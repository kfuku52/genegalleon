"""Lossless metadata retirement and public read/restart compatibility."""
import hashlib
import json
import os
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
if str(SUPPORT) not in sys.path:
    sys.path.insert(0, str(SUPPORT))
import metadata_store as store  # noqa: E402
import workflow_observation as observer  # noqa: E402


def closed_attempt(root, number=1, state="exited"):
    aid = f"{number:032x}"
    attempt = root / aid
    attempt.mkdir(parents=True)
    value = {"schema": observer.SCHEMA, "attempt_id": aid, "workflow": "gg_gene_evolution",
             "started_at_ns": 1, "execution_state": state, "completion_state": "unverified",
             "exit_code": 0 if state == "exited" else None, "accepted_exit_codes": [0]}
    if state == "exited":
        value.update(finished_at_ns=2, execution_accepted=True)
    observer.atomic_json(attempt / "run.json", value)
    key = observer.json_digest(["OG0000001", "mafft"])
    observer.atomic_json(attempt / ("contract-" + key + ".json"), {
        "schema": "genegalleon-contract-observation-v1", "attempt_id": aid,
        "family_id": "OG0000001", "step": "mafft", "operation": "record", "exit_code": 0,
        "observed_at_ns": 2, "argv": [], "manifest_sha256": None})
    return attempt


def pack(root, **kwargs):
    return store.pack(root, writers_stopped=True, delete_originals=True, **kwargs)


def test_pack_restore_and_idempotent_reimport(tmp_path):
    root = tmp_path / "resource_metrics"
    root.mkdir()
    path = root / "attempt-a.json"
    raw = b'{"workflow":"gg_gene_evolution","cpus":4,"wall_seconds":12.25}\n'
    path.write_bytes(raw)
    os.chmod(path, 0o640)
    os.utime(path, ns=(1234567890, 1234567890))
    original_mtime = path.stat().st_mtime_ns
    result = pack(root)
    assert result["removed_files"] == 1 and not path.exists()
    assert observer.read_regular_bytes(path) == raw
    assert store.archived_bytes(path) == raw
    with sqlite3.connect((root / store.DATABASE).as_uri() + "?mode=ro", uri=True) as db:
        assert db.execute("SELECT workflow,cpus,wall_seconds FROM records").fetchone() == ("gg_gene_evolution", 4, 12.25)
    assert pack(root)["removed_files"] == 0
    assert store.restore(root, writers_stopped=True)["restored_files"] == 1
    assert path.read_bytes() == raw and path.stat().st_mtime_ns == original_mtime
    assert path.stat().st_mode & 0o777 == 0o640
    assert store.restore(root, writers_stopped=True)["restored_files"] == 0
    path.write_bytes(b'{"cpus":8}')
    assert observer.read_regular_bytes(path) == b'{"cpus":8}'
    pack(root)
    assert store.archived_bytes(path) == b'{"cpus":8}'


def test_public_status_identical_with_live_and_archived_attempts(tmp_path):
    root = tmp_path / "observations"
    closed = closed_attempt(root)
    live = closed_attempt(root, 2, "started")
    argv = [sys.executable, str(SUPPORT / "workflow_api.py"), "status", "--directory", str(root),
            "--page-size", "512", "--max-records", "20"]
    before = json.loads(subprocess.run(argv, capture_output=True, text=True, check=True).stdout)
    result = pack(root)
    assert result["removed_files"] == 2
    assert not (closed / "run.json").exists() and (live / "run.json").exists()
    after = json.loads(subprocess.run(argv, capture_output=True, text=True, check=True).stdout)
    for key in ("records", "fingerprints", "inventory_count", "inventory_sha256", "root_identity"):
        assert after[key] == before[key]


def test_source_change_before_publication_keeps_originals(tmp_path, monkeypatch):
    root = tmp_path / "versions"
    root.mkdir()
    path = root / "run.test.log"
    path.write_bytes(b"original")
    original = store._read_source
    count = 0
    def changed(path, **kwargs):
        nonlocal count
        result = original(path, **kwargs)
        count += 1
        if count == 1:
            path.write_bytes(b"changed")
        return result
    monkeypatch.setattr(store, "_read_source", changed)
    with pytest.raises(ValueError, match="changed"):
        pack(root)
    assert path.read_bytes() == b"changed" and not (root / store.DATABASE).exists()


def test_source_change_after_publication_is_not_deleted(tmp_path, monkeypatch):
    root = tmp_path / "versions"
    root.mkdir()
    path = root / "run.test.log"
    path.write_bytes(b"original")
    original = store.os.replace
    def changed(source, destination):
        original(source, destination)
        if Path(destination).name == store.DATABASE:
            path.write_bytes(b"changed")
    monkeypatch.setattr(store.os, "replace", changed)
    with pytest.raises(ValueError, match="Source changed"):
        pack(root)
    assert path.read_bytes() == b"changed" and store.archived_bytes(path) == b"original"
    assert observer.read_regular_bytes(path) == b"changed"


def test_corrupt_body_stops_reads(tmp_path):
    root = tmp_path / "versions"
    root.mkdir()
    path = root / "run.test.log"
    path.write_bytes(b"exact")
    pack(root)
    with sqlite3.connect(root / store.DATABASE) as db:
        db.execute("UPDATE files SET sha256=?", ("0" * 64,))
    with pytest.raises(ValueError, match="corrupt"):
        store.archived_bytes(path)


def test_restore_does_not_overwrite_or_escape(tmp_path):
    root = tmp_path / "versions"
    root.mkdir()
    path = root / "run.test.log"
    path.write_bytes(b"exact")
    pack(root)
    path.write_bytes(b"user change")
    with pytest.raises(ValueError, match="overwrite"):
        store.restore(root, writers_stopped=True)
    assert path.read_bytes() == b"user change"
    path.unlink()
    with sqlite3.connect(root / store.DATABASE) as db:
        db.execute("UPDATE files SET path='../outside.log'")
    with pytest.raises(ValueError, match="Unsafe metadata key"):
        store.restore(root, writers_stopped=True)
    assert not (tmp_path / "outside.log").exists()


def test_symlink_hardlink_and_non_metadata_roots_rejected(tmp_path):
    root = tmp_path / "versions"
    root.mkdir()
    other = tmp_path / "owned.log"
    other.write_bytes(b"keep")
    path = root / "run.test.log"
    path.symlink_to(other)
    with pytest.raises(ValueError, match="regular"):
        pack(root)
    path.unlink()
    os.link(other, path)
    with pytest.raises(ValueError, match="regular"):
        pack(root)
    assert other.read_bytes() == b"keep"
    tsv = tmp_path / "stat_tree"
    tsv.mkdir()
    (tsv / "OG1.tsv").write_bytes(b"gene\tvalue\n")
    with pytest.raises(ValueError, match="registered metadata"):
        pack(tsv)


def test_offline_acknowledgement_is_required(tmp_path):
    root = tmp_path / "versions"
    root.mkdir()
    with pytest.raises(ValueError, match="Stop all workspace writers"):
        store.pack(root)


def test_observation_fence_detects_database_replacement(tmp_path):
    root = tmp_path / "observations"
    attempt = closed_attempt(root)
    pack(root)
    with pytest.raises(ValueError, match="changed"):
        with observer.observe_files():
            observer.read_regular_json(attempt / "run.json")
            info = (root / store.DATABASE).stat()
            os.utime(root / store.DATABASE, ns=(info.st_atime_ns, info.st_mtime_ns + 1_000_000_000))
            assert (root / store.DATABASE).stat().st_mtime_ns != info.st_mtime_ns


def test_inventory_journal_round_trip(tmp_path):
    root = tmp_path / "output-inventory"
    family = root / "ab" / ("a" * 64)
    family.mkdir(parents=True)
    path = family / "123.paths"
    raw = b"mafft/OG1.fa\0artifact_provenance/OG1.mafft.json\0"
    path.write_bytes(raw)
    pack(root)
    assert store.journal_paths(family) == [path]
    assert store.read_bytes(path) == raw
    (family / "inventory.paths").write_bytes(raw)
    assert store.journal_paths(family) == [family / "inventory.paths"]


def test_provenance_reader_and_managed_deletion_generations(tmp_path):
    import artifact_provenance as provenance
    from gene_family_output_store import GeneFamilyOutputStore
    root = tmp_path / "orthogroup"
    directory = root / "artifact_provenance"
    directory.mkdir(parents=True)
    path = directory / "OG0000001.mafft.json"
    raw = b'{"schema_version":1,"family_id":"OG0000001","step":"mafft"}\n'
    path.write_bytes(raw)
    pack(directory)
    assert provenance.declared_path_exists(path)
    assert provenance.load_manifest(path)["family_id"] == "OG0000001"
    logical = GeneFamilyOutputStore(root)
    assert logical.file_names("artifact_provenance") == [path.name]
    artifact = logical.source_artifact("artifact_provenance", path.name)
    assert artifact.metadata_path == directory / store.DATABASE
    assert provenance.read_manifest_from_store(logical, path.name, artifact)["step"] == "mafft"
    assert provenance.sha256_path(path) == (hashlib.sha256(raw).hexdigest(), len(raw), "file")
    logical.delete("artifact_provenance/" + path.name, family_id="OG0000001")
    assert logical.source_artifact("artifact_provenance", path.name) is None
    assert not provenance.declared_path_exists(path)
    with pytest.raises(provenance.ProvenanceError):
        provenance.load_manifest(path)
    assert store.restore(directory, writers_stopped=True)["restored_files"] == 0
    # A new publication after deletion receives a later native store generation.
    path.write_bytes(raw)
    pack(directory)
    logical = GeneFamilyOutputStore(root)
    assert logical.source_artifact("artifact_provenance", path.name) is not None
    logical.delete("artifact_provenance/" + path.name, family_id="OG0000001")
    assert logical.file_names("artifact_provenance") == []


def test_missing_legacy_metadata_namespace_stays_read_only(tmp_path):
    path = tmp_path / "artifact_provenance" / "OG1.mafft.json"
    assert store.archived_exists(path) is False
    assert store.archived_bytes(path) is None
    with pytest.raises(FileNotFoundError):
        observer.read_regular_json(path)
    assert not path.parent.exists()


def test_native_zip_index_sqlite_conversion_update_and_repair(tmp_path):
    from gene_family_output_store import (
        GeneFamilyOutputStore,
        archive_completed_outputs,
        orthogroup_id_from_name,
        repair_archive_index,
    )
    root = tmp_path / "orthogroup"
    path = root / "mafft" / "OG0000001.aln.fa"
    path.parent.mkdir(parents=True)
    path.write_bytes(b">gene\nACGT\n")
    assert archive_completed_outputs(root, "orthogroup", ["OG0000001"], orthogroup_id_from_name,
                                     include_incomplete=True)
    before = GeneFamilyOutputStore(root).source_artifact("mafft", path.name)
    for name in ("index", "index-by-subdir"):
        directory = root / ".gg_store" / name
        result = pack(directory)
        assert result["removed_files"] > 0 and not list(directory.glob("*.json"))
    logical = GeneFamilyOutputStore(root)
    assert logical.source_artifact("mafft", path.name) == before
    with logical.open_binary("mafft", path.name) as handle:
        assert handle.read() == b">gene\nACGT\n"
    with sqlite3.connect(root / ".gg_store/index" / store.DATABASE) as db:
        assert db.execute("SELECT logical_path,family_id FROM archive_entries").fetchall() == [("mafft/OG0000001.aln.fa", "OG0000001")]
    # Native archive writes retain SQLite storage and publish through the existing pending marker.
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b">gene\nTGCA\n")
    assert archive_completed_outputs(root, "orthogroup", ["OG0000001"], orthogroup_id_from_name,
                                     include_incomplete=True)
    for name in ("index", "index-by-subdir"):
        assert not list((root / ".gg_store" / name).glob("*.json"))
    with GeneFamilyOutputStore(root).open_binary("mafft", path.name) as handle:
        assert handle.read() == b">gene\nTGCA\n"
    repair_archive_index(root, remove_orphans=False)
    with GeneFamilyOutputStore(root).open_binary("mafft", path.name) as handle:
        assert handle.read() == b">gene\nTGCA\n"


def test_index_named_outside_native_store_is_not_a_metadata_root(tmp_path):
    root = tmp_path / "index"
    root.mkdir()
    (root / "scientific.json").write_bytes(b"{}")
    with pytest.raises(ValueError, match="registered metadata"):
        pack(root)


def test_same_operation_can_finish_only_its_original_retirement(tmp_path, monkeypatch):
    root = tmp_path / "versions"
    root.mkdir()
    for name in ("one.log", "two.log"):
        (root / name).write_bytes(name.encode())
    original = store._retire
    def interrupted(root, expected, operation_id):
        key = next(iter(expected))
        original(root, {key: expected[key]}, operation_id)
        raise RuntimeError("interrupted after one removal")
    monkeypatch.setattr(store, "_retire", interrupted)
    with pytest.raises(RuntimeError, match="interrupted"):
        pack(root, operation_id="test-operation")
    monkeypatch.setattr(store, "_retire", original)
    result = pack(root, operation_id="test-operation")
    assert result["idempotent_replay"] and result["removed_files"] == 2
    assert not list(root.glob("*.log"))
    (root / "one.log").write_bytes(b"a later writer")
    with pytest.raises(ValueError, match="Source changed"):
        pack(root, operation_id="test-operation")
    assert (root / "one.log").read_bytes() == b"a later writer"


def test_sqlite_version_cache_and_bounded_context_retention(tmp_path, monkeypatch):
    import version_report_identity as versions
    root = tmp_path / "version_inventory"
    versions.private_cache(root)
    keys = [f"{value:064x}" for value in range(4)]
    for i, key in enumerate(keys):
        raw = f"tool version {i}\n".encode()
        (root / (key + ".log")).write_bytes(raw)
        (root / (key + ".json")).write_text(json.dumps({"schema_version": 1, "key": key,
                                                       "sha256": hashlib.sha256(raw).hexdigest()}))
        os.utime(root / (key + ".log"), ns=((i + 1) * 1_000_000_000, (i + 1) * 1_000_000_000))
    pack(root)
    assert all(versions.cache_ready(root / (key + ".log"), key) for key in keys)
    store.pack(root, writers_stopped=True, cache_limit=2)
    assert store.cache_policy(root) == 2
    assert not versions.cache_ready(root / (keys[0] + ".log"), keys[0])
    assert versions.cache_ready(root / (keys[-1] + ".log"), keys[-1])
    before = versions.environment_identity()
    monkeypatch.setenv("GG_OBSERVATION_ATTEMPT_DIR", "/other/run-id")
    monkeypatch.setenv("GG_RESOURCE_EVENT_FILE", "/other/events-id")
    monkeypatch.setenv("GG_ARRAY_JOB_ID", "123456")
    assert versions.environment_identity() == before
    monkeypatch.setenv("PYTHONPATH", "/different/library")
    assert versions.environment_identity() != before


def test_indexed_execution_fields_preserve_original_json(tmp_path):
    root = tmp_path / "observations"
    attempt = closed_attempt(root)
    path = attempt / "run.json"
    value = json.loads(path.read_text())
    value["runtime"] = {"runtime_id": "image-identity", "scheduler": {"GG_ARRAY_JOB_ID": "123", "GG_ARRAY_TASK_ID": "4"}}
    raw = json.dumps(value).encode()
    path.write_bytes(raw)
    pack(root)
    assert observer.read_regular_bytes(path) == raw
    with sqlite3.connect(root / store.DATABASE) as db:
        assert db.execute("SELECT job_id,array_task_id,started_at,finished_at FROM records WHERE path=?",
                          (attempt.name + "/run.json",)).fetchone() == ("123", "4", 1e-9, 2e-9)


def test_cache_publish_and_copy_use_sqlite_policy(tmp_path):
    import version_report_identity as versions
    root = tmp_path / "version_inventory"
    versions.private_cache(root)
    pack(root)
    store.pack(root, writers_stopped=True, cache_limit=1)
    source = tmp_path / "versions.log"
    destination = tmp_path / "copied.log"
    for i in range(2):
        source.write_bytes(f"fixture tool {i}\n".encode())
        key = f"{i:064x}"
        path = root / (key + ".log")
        subprocess.run([sys.executable, str(SUPPORT / "version_report_identity.py"), "publish-cache",
                        "--cache-file", str(path), "--key", key, "--source", str(source)], check=True)
        assert not path.exists()
        assert versions.cache_ready(path, key)
        subprocess.run([sys.executable, str(SUPPORT / "version_report_identity.py"), "copy-cache",
                        "--cache-file", str(path), "--key", key, "--destination", str(destination)], check=True)
        assert destination.read_bytes() == source.read_bytes()
    assert not versions.cache_ready(root / (f"{0:064x}" + ".log"), f"{0:064x}")
    assert store.cache_policy(root) == 1


def test_interrupted_retirement_flushes_completed_directory_and_resumes(tmp_path, monkeypatch):
    root = tmp_path / "versions"
    root.mkdir()
    for name in ["one.log", "two.log"]:
        (root / name).write_bytes(name.encode())
    unlink, sync = Path.unlink, store._sync_dir
    events = []
    def interrupted(path, *args, **kwargs):
        if path == root / "two.log":
            raise RuntimeError("retirement interrupted")
        result = unlink(path, *args, **kwargs)
        if path == root / "one.log":
            events.append("removed-one")
        return result
    def synced(path):
        events.append("synced-" + path.name)
        return sync(path)
    monkeypatch.setattr(Path, "unlink", interrupted)
    monkeypatch.setattr(store, "_sync_dir", synced)
    with pytest.raises(RuntimeError, match="retirement interrupted"):
        pack(root, operation_id="interrupted-directory-flush")
    assert events.index("removed-one") < len(events) - 1 and events[-1] == "synced-versions"
    assert store.archived_bytes(root / "one.log") == b"one.log"
    assert (root / "two.log").read_bytes() == b"two.log"
    monkeypatch.setattr(Path, "unlink", unlink)
    assert pack(root, operation_id="interrupted-directory-flush")["removed_this_call"] == 1


def test_flat_record_collection_preserves_opaque_lock_namespaces(tmp_path):
    root = tmp_path / "version_inventory"
    root.mkdir()
    (root / "fixture.log").write_bytes(b"verified version body")
    control = root / "fixture.log.lock" / "gate"
    control.mkdir(parents=True)
    outside = tmp_path / "unrelated-control"
    outside.mkdir()
    (outside / "owner.json").write_bytes(b"coordination state")
    pointer = control / "owner"
    pointer.symlink_to(outside, target_is_directory=True)
    before = pointer.lstat().st_ino
    assert pack(root)["removed_files"] == 1
    assert store.archived_bytes(root / "fixture.log") == b"verified version body"
    assert pointer.is_symlink() and pointer.lstat().st_ino == before
    assert (outside / "owner.json").read_bytes() == b"coordination state"


def coarse_database_stamps(monkeypatch):
    original = store.signature
    def coarse(path):
        value = original(path)
        return (*value[:4], 0, 0, *value[6:]) if Path(path).name == store.DATABASE else value
    monkeypatch.setattr(store, "signature", coarse)


@pytest.mark.parametrize("mutation", ["corrupt_sha", "unsafe_path"])
def test_same_stamp_sqlite_commits_invalidate_cached_reads(tmp_path, monkeypatch, mutation):
    root = tmp_path / "versions"
    root.mkdir()
    path = root / "run.test.log"
    path.write_bytes(b"exact")
    coarse_database_stamps(monkeypatch)
    pack(root)
    assert store.archived_bytes(path) == b"exact"
    database = root / store.DATABASE
    before = store.signature(database)
    with sqlite3.connect(database) as db:
        if mutation == "corrupt_sha":
            db.execute("UPDATE files SET sha256=?", ("0" * 64,))
        else:
            db.execute("UPDATE files SET path='../outside.log'")
    assert store.signature(database) == before
    with pytest.raises(ValueError, match="corrupt" if mutation == "corrupt_sha" else "Unsafe"):
        if mutation == "corrupt_sha":
            store.archived_bytes(path)
        else:
            store.restore(root, writers_stopped=True)
    assert not path.exists() and not (tmp_path / "outside.log").exists()


def test_observation_fence_sees_sqlite_commit_with_coarse_stamps(tmp_path, monkeypatch):
    root = tmp_path / "observations"
    attempt = closed_attempt(root)
    coarse_database_stamps(monkeypatch)
    original = observer.path_signature
    def coarse(path):
        value = original(path)
        return (*value[:4], 0, 0) if value is not None and Path(path).name == store.DATABASE else value
    monkeypatch.setattr(observer, "path_signature", coarse)
    pack(root)
    with pytest.raises(ValueError, match="changed"):
        with observer.observe_files():
            observer.read_regular_json(attempt / "run.json")
            with sqlite3.connect(root / store.DATABASE) as db:
                db.execute("UPDATE files SET sha256=?", ("0" * 64,))


def test_tied_cache_timestamps_keep_latest_publication(tmp_path):
    import version_report_identity as versions
    root = tmp_path / "version_inventory"
    versions.private_cache(root)
    keys = [f"{i:064x}" for i in range(3)]
    for i, key in enumerate(keys):
        path = root / (key + ".log")
        raw = f"tool version {i}\n".encode()
        path.write_bytes(raw)
        (root / (key + ".json")).write_text(json.dumps({"schema_version": 1, "key": key,
                                                       "sha256": hashlib.sha256(raw).hexdigest()}))
        os.utime(path, ns=(0, 0))
        pack(root)
        store.pack(root, writers_stopped=True, cache_limit=1)
        assert versions.cache_ready(path, key)
        assert all(not versions.cache_ready(root / (old + ".log"), old) for old in keys[:i])
        with sqlite3.connect(root / store.DATABASE) as db:
            assert db.execute("SELECT mtime_ns FROM files WHERE path=?", (path.name,)).fetchone() == (0,)


def test_observation_inventory_overlaps_independent_run_reads(tmp_path, monkeypatch):
    import threading
    root = tmp_path / "observations"
    closed = [closed_attempt(root, i) for i in range(8)]
    live = closed_attempt(root, 8, "started")
    original = store._read_source
    lock = threading.Lock()
    overlap = threading.Event()
    active = 0
    def read(path, **kwargs):
        nonlocal active
        if Path(path).name != "run.json":
            return original(path, **kwargs)
        with lock:
            active += 1
            if active >= 2:
                overlap.set()
        try:
            assert overlap.wait(2), "run headers were serialized"
            return original(path, **kwargs)
        finally:
            with lock:
                active -= 1
    monkeypatch.setattr(store, "_read_source", read)
    expected = sorted(str(p.relative_to(root)) for directory in closed for p in directory.glob("*.json"))
    assert store._inventory(root, 100) == expected
    assert (live / "run.json").exists() and active == 0


def test_observation_prefetch_rejects_symlink_before_read(tmp_path):
    root = tmp_path / "observations"
    root.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    path = outside / "run.json"
    path.write_bytes(b"do not read or retire")
    (root / ("a" * 32)).symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError, match="Symlinked"):
        pack(root)
    assert path.read_bytes() == b"do not read or retire" and not (root / store.DATABASE).exists()


def test_snapshot_checks_only_bounded_database_headers(tmp_path, monkeypatch):
    root = tmp_path / "versions"
    root.mkdir()
    for i in range(100):
        (root / f"run.{i}.log").write_bytes(f"version {i}".encode())
    pack(root)
    original = os.fdopen
    sizes = []
    class HeaderReader:
        def __init__(self, stream):
            self.stream = stream
        def __enter__(self):
            self.stream.__enter__()
            return self
        def __exit__(self, *args):
            return self.stream.__exit__(*args)
        def fileno(self):
            return self.stream.fileno()
        def read(self, size=-1):
            sizes.append(size)
            return self.stream.read(size)
    monkeypatch.setattr(store.os, "fdopen", lambda *args, **kwargs: HeaderReader(original(*args, **kwargs)))
    monkeypatch.setattr(store, "_file_sha", lambda *_: pytest.fail("full database hashing during lookup"))
    with store.read_snapshot():
        for i in range(100):
            assert store.archived_bytes(root / f"run.{i}.log") == f"version {i}".encode()
    assert sizes and -1 not in sizes and sum(sizes) <= 400


def test_collection_prefetch_preserves_order_and_bounded_read_ahead(tmp_path, monkeypatch):
    import threading
    import time
    root = tmp_path / "versions"
    root.mkdir()
    bodies = {f"{i:03d}.log": json.dumps({"workflow": "gg_gene_evolution", "number": i}).encode() for i in range(64)}
    for key, raw in bodies.items():
        (root / key).write_bytes(raw)
    identities = {key: list(store.signature(root / key)) for key in bodies}
    original_read, original_fields = store._read_source, store._fields
    lock, overlap = threading.Lock(), threading.Event()
    active = peak = started = 0
    checked = False
    def read(path, **kwargs):
        nonlocal active, peak, started
        with lock:
            active += 1
            peak = max(peak, active)
            started += 1
            if active >= 2:
                overlap.set()
        try:
            assert overlap.wait(2), "metadata source bodies were serialized"
            return original_read(path, **kwargs)
        finally:
            with lock:
                active -= 1
    def fields(raw):
        nonlocal checked
        if not checked:
            checked = True
            # Completed future results must not accumulate for the full input.
            time.sleep(0.03)
            with lock:
                assert 2 <= started <= 16
        return original_fields(raw)
    monkeypatch.setattr(store, "_read_source", read)
    monkeypatch.setattr(store, "_fields", fields)
    result = pack(root, operation_id="bounded-source-prefetch")
    assert result["removed_files"] == len(bodies) and active == 0 and 2 <= peak <= 16
    with sqlite3.connect(root / store.DATABASE) as db:
        assert [r[0] for r in db.execute("SELECT path FROM files ORDER BY rowid")] == sorted(bodies)
        assert {key: json.loads(value) for key, value in db.execute("SELECT path,source_identity FROM operation_files")} == identities
    assert {key: store.archived_bytes(root / key) for key in bodies} == bodies


def test_collection_prefetch_failure_joins_workers_and_keeps_sources(tmp_path, monkeypatch):
    import threading
    import time
    root = tmp_path / "versions"
    root.mkdir()
    bodies = {f"{i:03d}.log": f"version {i}".encode() for i in range(8)}
    for key, raw in bodies.items():
        (root / key).write_bytes(raw)
    original = store._read_source
    other_started, lock = threading.Event(), threading.Lock()
    active = 0
    def read(path, **kwargs):
        nonlocal active
        with lock:
            active += 1
        try:
            if Path(path).name == "000.log":
                assert other_started.wait(2), "failed read did not overlap another worker"
                raise ValueError("prefetched source failed validation")
            other_started.set()
            time.sleep(0.02)
            return original(path, **kwargs)
        finally:
            with lock:
                active -= 1
    monkeypatch.setattr(store, "_read_source", read)
    with pytest.raises(ValueError, match="prefetched source failed validation"):
        pack(root)
    assert active == 0 and not (root / store.DATABASE).exists()
    assert not list(root.glob(store.DATABASE + ".partial.*"))
    assert {key: (root / key).read_bytes() for key in bodies} == bodies


@pytest.mark.parametrize("name", ["index", "index-by-subdir"])
def test_native_index_read_budget_stays_serial(tmp_path, monkeypatch, name):
    from gene_family_output_store import archive_completed_outputs, orthogroup_id_from_name
    root = tmp_path / "orthogroup"
    families = [f"OG{i:07d}" for i in range(16)]
    for subdir in ("mafft", "muscle"):
        directory = root / subdir
        directory.mkdir(parents=True)
        for family in families:
            (directory / (family + ".aln.fa")).write_bytes(b">gene\nACGT\n")
    assert archive_completed_outputs(root, "orthogroup", families, orthogroup_id_from_name, include_incomplete=True)
    directory = root / ".gg_store" / name
    assert len(list(directory.glob("*.json"))) >= 2
    def no_pool(*args, **kwargs):
        pytest.fail("128 MiB native index records must not be prefetched concurrently")
    monkeypatch.setattr(store, "ThreadPoolExecutor", no_pool)
    assert pack(directory)["removed_files"] >= 2


def test_prefetched_source_mutation_with_same_signature_prevents_publication(tmp_path, monkeypatch):
    root = tmp_path / "versions"
    root.mkdir()
    bodies = {f"{i:03d}.log": b"original" for i in range(8)}
    for key, raw in bodies.items():
        (root / key).write_bytes(raw)
    path = root / "000.log"
    identity = store.signature(path)
    signature, fields = store.signature, store._fields
    changed = False
    def same_stamp(candidate):
        return identity if Path(candidate) == path else signature(candidate)
    def mutate(raw):
        nonlocal changed
        if not changed:
            changed = True
            path.write_bytes(b"modified")
        return fields(raw)
    monkeypatch.setattr(store, "signature", same_stamp)
    monkeypatch.setattr(store, "_fields", mutate)
    with pytest.raises(ValueError, match="Metadata source changed before publication"):
        pack(root)
    assert not (root / store.DATABASE).exists()
    assert path.read_bytes() == b"modified"
    assert all((root / key).read_bytes() == value for key, value in bodies.items() if key != path.name)


def test_prefetch_caller_failure_joins_readers_before_cleanup(tmp_path, monkeypatch):
    import threading
    import time
    root = tmp_path / "versions"
    root.mkdir()
    bodies = {f"{i:03d}.log": b"original" for i in range(8)}
    for key, raw in bodies.items():
        (root / key).write_bytes(raw)
    original = store._read_source
    other_started, failed, lock = threading.Event(), threading.Event(), threading.Lock()
    active = 0
    def read(path, **kwargs):
        nonlocal active
        with lock:
            active += 1
        try:
            if Path(path).name == "000.log":
                assert other_started.wait(2)
            else:
                other_started.set()
                assert failed.wait(2)
                time.sleep(0.02)
            return original(path, **kwargs)
        finally:
            with lock:
                active -= 1
    def caller_failed(raw):
        failed.set()
        raise ValueError("collector caller failed")
    monkeypatch.setattr(store, "_read_source", read)
    monkeypatch.setattr(store, "_fields", caller_failed)
    with pytest.raises(ValueError, match="collector caller failed"):
        pack(root)
    assert active == 0 and not (root / store.DATABASE).exists()
    assert {key: (root / key).read_bytes() for key in bodies} == bodies
