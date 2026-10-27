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
    result = pack(root)
    assert result["removed_files"] == 1 and not path.exists()
    assert observer.read_regular_bytes(path) == raw
    assert store.archived_bytes(path) == raw
    with sqlite3.connect((root / store.DATABASE).as_uri() + "?mode=ro", uri=True) as db:
        assert db.execute("SELECT workflow,cpus,wall_seconds FROM records").fetchone() == ("gg_gene_evolution", 4, 12.25)
    assert pack(root)["removed_files"] == 0
    assert store.restore(root, writers_stopped=True)["restored_files"] == 1
    assert path.read_bytes() == raw and path.stat().st_mtime_ns == 1234567890
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
            os.utime(root / store.DATABASE, ns=(info.st_atime_ns, info.st_mtime_ns + 1))


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
        os.utime(root / (key + ".log"), ns=(i + 1, i + 1))
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
