"""Public preview is read-only and agrees with independently authorized maintenance."""
import json
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
import metadata_maintenance as maintenance  # noqa: E402
import metadata_plan as planner  # noqa: E402
import metadata_store as metadata  # noqa: E402
from shared_namespace_lock import acquire, release  # noqa: E402
from test_metadata_maintenance import cache, collected, legacy, observation  # noqa: E402


def snapshot(root):
    return {str(p.relative_to(root)): (p.lstat().st_mode, p.lstat().st_mtime_ns,
                                     p.lstat().st_ctime_ns, p.read_bytes() if p.is_file() else None)
            for p in [root, *root.rglob("*")]}


def test_preview_never_creates_controls_and_collection_membership_matches(tmp_path):
    root = tmp_path / "observations"
    observation(root, 1)
    observation(root, 2, closed=False)
    scientific = root / "annotation.tsv"
    scientific.write_bytes(b"keep\n")
    before = snapshot(root)
    result = planner.plan(root)
    assert snapshot(root) == before
    assert result["selected_files"] == 1 and result["protected_items"] == 2
    assert result["eligible"] and result["database_sha256"] == ""
    actual = metadata.pack(root, writers_stopped=True, delete_originals=True)
    assert actual["removed_files"] == result["selected_files"]
    assert scientific.read_bytes() == b"keep\n"


def test_archived_directory_plan_agrees_with_maintenance_and_preserves_database(tmp_path):
    root = tmp_path / "observations"
    observation(root, 1)
    observation(root, 2, age=60)
    collected(root)
    before = snapshot(root)
    result = planner.plan(root, operation="prune-observations")
    assert snapshot(root) == before
    assert result["selected_directories"] == 1 and result["protected_items"] == 1
    assert len(result["database_sha256"]) == 64 and result["database_generation"] >= 0
    actual = maintenance.maintain(root, tmp_path / "journal", kind="prune-observations",
                                 operation_id="fixture", writers_stopped=True)
    assert actual["removed_entries"] == result["selected_directories"]


def test_version_plan_preserves_unknown_owner_and_is_bounded(tmp_path):
    root = tmp_path / "version_inventory"
    cache(root)
    legacy(root, "b" * 64)
    owned = legacy(root, "c" * 64)
    token = acquire(owned[0], exclusive=True)
    try:
        before = snapshot(root)
        result = planner.plan(root, operation="retire-version-locks")
        assert snapshot(root) == before
        assert result["selected_files"] == 1 and result["selected_directories"] == 2
        assert result["protected_items"] == 1
        with pytest.raises(ValueError, match="bound"):
            planner.plan(root, max_entries=1)
    finally:
        release(owned[0], token, exclusive=True)


def test_namespace_ownership_age_threshold_and_no_members_are_explicit(tmp_path):
    root = tmp_path / "resource_metrics"
    root.mkdir()
    (root / "metric.json").write_text("{}")
    token = acquire(root / ".gg-metadata.guard", exclusive=False)
    try:
        result = planner.plan(root, minimum_files=2, minimum_idle_seconds=60)
        assert not result["eligible"]
        assert set(result["deferred_reasons"]) == {"namespace-owned", "below-minimum-entries", "sources-too-recent"}
    finally:
        release(root / ".gg-metadata.guard", token, exclusive=False)
    (root / "metric.json").unlink()
    assert planner.plan(root)["deferred_reasons"] == ["no-eligible-entries"]


def test_database_and_source_drift_are_rejected(tmp_path, monkeypatch):
    root = tmp_path / "observations"
    path = observation(root, 1) / "run.json"
    original = metadata._inventory

    def change(root, limit):
        selected = original(root, limit)
        path.write_text(path.read_text() + " ")
        return selected

    monkeypatch.setattr(metadata, "_inventory", change)
    with pytest.raises(ValueError, match="changed"):
        planner.plan(root)


def test_public_cli_and_unsafe_root_do_not_write(tmp_path):
    root = tmp_path / "resource_metrics"
    root.mkdir()
    (root / "metric.json").write_text("{}")
    before = snapshot(root)
    result = subprocess.run([sys.executable, str(SUPPORT / "metadata_plan.py"), str(root)],
                            text=True, capture_output=True, check=True)
    assert json.loads(result.stdout)["selected_files"] == 1
    assert snapshot(root) == before
    link = tmp_path / "alias"
    link.symlink_to(root, target_is_directory=True)
    with pytest.raises(ValueError, match="physical"):
        planner.plan(link)
    unsafe = root / "outside.json"
    unsafe.symlink_to(root / "metric.json")
    with pytest.raises(ValueError, match="regular"):
        planner.plan(root)


def test_database_generation_changed_during_plan_is_rejected(tmp_path, monkeypatch):
    root = tmp_path / "observations"
    observation(root, 1)
    collected(root)
    original = planner._sha

    def change(path):
        value = original(path)
        with sqlite3.connect(path) as db:
            db.execute("UPDATE metadata SET value='changed' WHERE key='generation'")
        return value

    monkeypatch.setattr(planner, "_sha", change)
    with pytest.raises(ValueError, match="changed"):
        planner.plan(root)
