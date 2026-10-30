"""Exact offline selection, recovery and refusal preserve scientific/control data."""

import json
import sqlite3
import subprocess
import sys
import time
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
if str(SUPPORT) not in sys.path:
    sys.path.insert(0, str(SUPPORT))
import metadata_maintenance as maintenance  # noqa: E402
import metadata_store as metadata  # noqa: E402
import workflow_observation as observer  # noqa: E402
from shared_namespace_lock import acquire, release  # noqa: E402


def snapshot(root):
    return {str(p.relative_to(root)): p.read_bytes() for p in root.rglob("*") if p.is_file()}


def observation(root, number, *, closed=True, age=172800):
    aid = f"{number:032x}"
    path = root / aid
    path.mkdir(parents=True)
    value = {
        "schema": observer.SCHEMA,
        "attempt_id": aid,
        "workflow": "gg_gene_evolution",
        "execution_state": "exited" if closed else "started",
        "completion_state": "unverified",
        "started_at_ns": 1,
        "exit_code": 0 if closed else None,
        "accepted_exit_codes": [0],
        "execution_accepted": closed,
        "finished_at_ns": time.time_ns() - int(age * 10**9),
    }
    (path / "run.json").write_text(json.dumps(value) + "\n")
    return path


def collected(root):
    metadata.pack(root, writers_stopped=True, delete_originals=True)


def maintain(root, work, kind):
    return maintenance.maintain(root, work, kind=kind, operation_id="fixture-operation", writers_stopped=True)


def test_completed_observation_directory_retirement_replay_and_exact_restore(tmp_path):
    root = tmp_path / "observations"
    old = observation(root, 1)
    live = observation(root, 2, closed=False)
    recent = observation(root, 3, age=60)
    collected(root)
    before = {name: metadata.archived_bytes(root / name / "run.json") for name in (old.name, recent.name)}
    original = old.stat()
    public_before = metadata.archived_children(root, directories=True)
    database = (root / metadata.DATABASE).read_bytes()
    work = tmp_path / "journal"
    result = maintain(root, work, "prune-observations")
    assert result["removed_entries"] == 1 and result["removed_this_call"] == 1
    assert not old.exists() and live.exists() and recent.exists()
    assert metadata.archived_children(root, directories=True) == public_before
    assert {name: metadata.archived_bytes(root / name / "run.json") for name in before} == before
    assert (root / metadata.DATABASE).read_bytes() == database
    assert maintain(root, work, "prune-observations")["removed_this_call"] == 0
    assert maintenance.verify(root, work)
    assert maintenance.restore(root, work, writers_stopped=True)["restored_entries"] == 1
    assert old.stat().st_mode == original.st_mode and old.stat().st_mtime_ns == original.st_mtime_ns
    assert metadata.archived_children(root, directories=True) == public_before


def test_missing_unarchived_nonempty_recent_and_live_observation_directories_remain(tmp_path):
    root = tmp_path / "observations"
    old = observation(root, 1)
    collected(root)
    (old / "new.txt").write_text("keep")
    unknown = root / ("a" * 32)
    unknown.mkdir()
    result = maintain(root, tmp_path / "journal", "prune-observations")
    assert result["removed_entries"] == 0 and result["protected_items"] == 2
    assert (old / "new.txt").read_text() == "keep" and unknown.exists()


def cache(root):
    root.mkdir()
    (root / ("a" * 64 + ".log")).write_text("version\n")
    collected(root)
    return (root / metadata.DATABASE).read_bytes()


def legacy(root, name):
    marker = root / (name + ".log.lock")
    marker.write_bytes(b"")
    marker.chmod(0o600)
    readers = root / (marker.name + ".namespace-v1") / "readers"
    readers.mkdir(parents=True)
    return marker, readers.parent, readers


def test_legacy_cache_coordination_removal_preserves_bounded_locks_database_and_restore(tmp_path):
    root = tmp_path / "version_inventory"
    database = cache(root)
    paths = legacy(root, "b" * 64)
    before = {str(p.relative_to(root)): p.stat() for p in paths}
    unrelated = root / "notes"
    unrelated.write_text("keep")
    work = tmp_path / "journal"
    result = maintain(root, work, "retire-version-locks")
    assert result["removed_entries"] == 3 and all(not p.exists() for p in paths)
    assert unrelated.read_text() == "keep" and (root / ".versions-cache.lock").exists()
    assert (root / metadata.DATABASE).read_bytes() == database
    assert maintain(root, work, "retire-version-locks")["removed_this_call"] == 0
    maintenance.restore(root, work, writers_stopped=True)
    for p in paths:
        original = before[str(p.relative_to(root))]
        assert p.stat().st_mode == original.st_mode and p.stat().st_mtime_ns == original.st_mtime_ns
    assert paths[0].read_bytes() == b""


def test_unknown_cache_owner_and_nonempty_marker_are_retained_independently(tmp_path):
    root = tmp_path / "version_inventory"
    cache(root)
    paths = legacy(root, "b" * 64)
    token = acquire(paths[0], exclusive=True)
    changed = legacy(root, "c" * 64)
    changed[0].write_bytes(b"keep")
    eligible = legacy(root, "d" * 64)
    try:
        result = maintain(root, tmp_path / "journal", "retire-version-locks")
        assert result["removed_entries"] == 3 and result["protected_items"] == 2
        assert all(p.exists() for p in paths + changed) and all(not p.exists() for p in eligible)
    finally:
        release(paths[0], token, exclusive=True)


@pytest.mark.parametrize("kind", ["prune-observations", "retire-version-locks"])
def test_offline_requirement_and_journal_escape_are_enforced(tmp_path, kind):
    root = tmp_path / ("observations" if kind == "prune-observations" else "version_inventory")
    if kind == "prune-observations":
        observation(root, 1)
        collected(root)
    else:
        cache(root)
    before = snapshot(root)
    with pytest.raises(ValueError, match="Stop all"):
        maintenance.maintain(root, tmp_path / "journal", kind=kind, operation_id="fixture")
    assert snapshot(root) == before
    with pytest.raises(ValueError, match="journal location"):
        maintain(root, root / "journal", kind)


def test_interruption_resumes_original_intent_without_adopting_replaced_directory(tmp_path, monkeypatch):
    root = tmp_path / "observations"
    old = observation(root, 1)
    second = observation(root, 2)
    collected(root)
    work = tmp_path / "journal"
    original = Path.rmdir
    calls = 0

    def interrupted(path):
        nonlocal calls
        if path.parent == root:
            calls += 1
            if calls == 2:
                raise OSError("fixture interruption")
        return original(path)

    with monkeypatch.context() as scope:
        scope.setattr(Path, "rmdir", interrupted)
        with pytest.raises(OSError, match="interruption"):
            maintain(root, work, "prune-observations")
    assert not old.exists() and second.exists()
    second.rmdir()
    second.mkdir()
    with pytest.raises(ValueError, match="replaced|metadata changed"):
        maintain(root, work, "prune-observations")
    assert second.exists()


def test_database_or_manifest_drift_prevents_further_retirement(tmp_path):
    root = tmp_path / "observations"
    observation(root, 1)
    collected(root)
    work = tmp_path / "journal"
    maintain(root, work, "prune-observations")
    with sqlite3.connect(root / metadata.DATABASE) as db:
        db.execute("UPDATE metadata SET value='changed' WHERE key='generation'")
    with pytest.raises(ValueError, match="database changed"):
        maintain(root, work, "prune-observations")


def test_symlink_and_unsealed_manifest_are_not_adopted(tmp_path):
    root = tmp_path / "observations"
    observation(root, 1)
    collected(root)
    outside = tmp_path / "outside"
    outside.mkdir()
    work = tmp_path / "journal"
    work.symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError, match="Symlinked"):
        maintain(root, work, "prune-observations")
    work.unlink()
    work.mkdir(mode=0o700)
    (work / "source.jsonl.gz").write_bytes(b"unknown")
    with pytest.raises(ValueError, match="Unsealed"):
        maintain(root, work, "prune-observations")


def test_public_observation_status_keeps_the_same_logical_records(tmp_path):
    root = tmp_path / "observations"
    observation(root, 1)
    collected(root)
    argv = [
        sys.executable,
        str(SUPPORT / "workflow_api.py"),
        "status",
        "--directory",
        str(root),
        "--page-size",
        "512",
        "--max-records",
        "20",
    ]
    before = json.loads(subprocess.run(argv, check=True, capture_output=True, text=True).stdout)
    maintain(root, tmp_path / "journal", "prune-observations")
    result = subprocess.run(argv, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    after = json.loads(result.stdout)
    for key in ("records", "fingerprints", "inventory_count", "inventory_sha256", "root_identity"):
        assert after[key] == before[key]


def test_collected_cache_uses_one_lock_across_many_contexts(tmp_path):
    import version_report_identity as versions

    root = tmp_path / "version_inventory"
    cache(root)
    selected = {versions.cache_lock_path(root / (f"{i:064x}" + ".log")) for i in range(10000)}
    assert selected == {root / ".versions-cache.lock"}
    cold = tmp_path / "cold"
    assert versions.cache_lock_path(cold / ("a" * 64 + ".log")) == cold / ("a" * 64 + ".log.lock")
    for _ in range(16):
        token = acquire(next(iter(selected)), exclusive=True)
        release(next(iter(selected)), token, exclusive=True)
    assert not list(root.glob("*.log.lock"))


def test_unarchived_absent_attempt_is_not_treated_as_virtual(tmp_path):
    root = tmp_path / "observations"
    observation(root, 1)
    collected(root)
    sys.path.insert(0, str(SUPPORT))
    import workflow_api

    with pytest.raises(FileNotFoundError):
        workflow_api.inventory(root / ("b" * 32), pattern="contract-")
