"""Completed runtime scripts retain bytes and resume fences after SQLite collection."""
import hashlib
import importlib
import os
import shlex
import sqlite3
import stat
import subprocess
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"


@pytest.fixture
def modules(monkeypatch):
    monkeypatch.syspath_prepend(str(SUPPORT))
    return importlib.import_module("metadata_store"), importlib.import_module("runtime_snapshot_storage")


def snapshots(tmp_path, runs=8):
    root = tmp_path / "workspace/output/runtime"
    pool = root / ".shared-code-v1"
    pool.mkdir(parents=True)
    files = dict(VERSION=b"0.8.test\n", gg_common_params_sh=b"original parameters\n",
                 gg_gene_evolution_core_sh=b"echo science\n", gg_gene_evolution_entrypoint_sh=b"entrypoint\n")
    result = {}
    for name, raw in files.items():
        name = name.replace("_sh", ".sh")
        obj = pool / hashlib.sha256(raw).hexdigest()
        obj.write_bytes(raw)
        obj.chmod(0o444)
        for i in range(runs):
            run = root / "gg_gene_evolution_entrypoint" / f"123_{i+1}"
            run.mkdir(parents=True, exist_ok=True)
            path = run / name
            os.link(obj, path)
            result[path.relative_to(root).as_posix()] = (raw, stat.S_IMODE(path.stat().st_mode), path.stat().st_mtime_ns)
    return root, result


def test_linked_runtime_roundtrip_prunes_runs_but_preserves_pool_and_unknown_files(tmp_path, modules):
    metadata, runtime = modules
    root, before = snapshots(tmp_path)
    unknown = root / "gg_gene_evolution_entrypoint/123_1/notes.txt"
    unknown.write_bytes(b"user notes\n")
    inode = next((root / ".shared-code-v1").iterdir()).stat().st_ino
    result = metadata.pack(root, delete_originals=True, writers_stopped=True, operation_id="roundtrip")
    assert result["selected_files"] == result["removed_files"] == 32
    with sqlite3.connect(root / metadata.DATABASE) as db:
        assert db.execute("SELECT count(*) FROM files").fetchone() == (32,)
        assert db.execute("SELECT count(*) FROM runtime_bodies").fetchone() == (4,)
        assert db.execute("SELECT sum(length(body)) FROM files").fetchone() == (0,)
    assert runtime.prune(root, writers_stopped=True) == 7
    assert unknown.read_bytes() == b"user notes\n"
    assert next((root / ".shared-code-v1").iterdir()).stat().st_ino == inode
    for name, (raw, mode, mtime) in before.items():
        assert not (root / name).exists()
        assert metadata.archived_bytes(root / name) == raw
        assert metadata.archived_info(root / name)["mtime_ns"] == mtime
    replay = metadata.pack(root, delete_originals=True, writers_stopped=True, operation_id="roundtrip")
    assert replay["idempotent_replay"] and replay["removed_this_call"] == 0
    assert runtime.prune(root, writers_stopped=True) == 0
    metadata.restore(root, writers_stopped=True)
    for name, (raw, mode, mtime) in before.items():
        path = root / name
        assert (path.read_bytes(), stat.S_IMODE(path.stat().st_mode), path.stat().st_mtime_ns) == (raw, mode, mtime)


def test_prepare_restores_only_same_run_and_rejects_changed_source_before_any_restore(tmp_path, modules):
    metadata, runtime = modules
    root, before = snapshots(tmp_path, runs=2)
    metadata.pack(root, delete_originals=True, writers_stopped=True, operation_id="prepare")
    runtime.prune(root, writers_stopped=True)
    run = root / "gg_gene_evolution_entrypoint/123_1"
    source = tmp_path / "source"
    source.mkdir()
    sources = []
    for name, (raw, _mode, _mtime) in before.items():
        if "/123_1/" in name:
            path = source / Path(name).name
            path.write_bytes(raw)
            sources.append(path)
    changed = next(p for p in sources if p.name.endswith("core.sh"))
    raw = changed.read_bytes()
    changed.write_bytes(b"changed source\n")
    with pytest.raises(ValueError, match="content changed"):
        runtime.prepare(run, sources)
    assert not run.exists()
    changed.write_bytes(raw)
    assert runtime.prepare(run, sources) == 4
    assert runtime.prepare(run, sources) == 0
    assert not (run.parent / "123_2").exists()
    for path in run.iterdir():
        expected = before[path.relative_to(root).as_posix()]
        assert (path.read_bytes(), stat.S_IMODE(path.stat().st_mode), path.stat().st_mtime_ns) == expected
    with pytest.raises(ValueError, match="membership changed"):
        runtime.prepare(run, sources[:-1])


def test_snapshot_collection_retains_strict_membership_permissions_and_writers_stop(tmp_path, modules):
    metadata, runtime = modules
    root, _before = snapshots(tmp_path, runs=1)
    with pytest.raises(ValueError, match="Stop all"):
        metadata.pack(root, delete_originals=True)
    with pytest.raises(ValueError, match="Stop all"):
        runtime.prune(root)
    assert not runtime.eligible("../outside/core.sh")
    assert not runtime.eligible("gg_gene_evolution_entrypoint/../VERSION")
    assert not runtime.eligible("gg_gene_evolution_entrypoint/123_1/unknown.sh")
    linked = root / "gg_gene_evolution_entrypoint/123_1/VERSION"
    linked.chmod(0o644)
    with pytest.raises(ValueError, match="unlinked regular"):
        metadata.pack(root, delete_originals=True, writers_stopped=True)
    assert linked.exists()
    linked.chmod(0o444)
    other = root / "bad_link"
    other.symlink_to(tmp_path, target_is_directory=True)
    with pytest.raises(ValueError, match="Symlinked metadata"):
        metadata.pack(root, delete_originals=True, writers_stopped=True)
    assert linked.exists()


def test_runtime_retirement_rejects_external_link_or_byte_change(tmp_path, modules):
    metadata, runtime = modules
    root, before = snapshots(tmp_path, runs=2)
    metadata.pack(root, delete_originals=False, writers_stopped=True, operation_id="fence")
    with sqlite3.connect(root / metadata.DATABASE) as db:
        import json
        expected = {p:(json.loads(i),sha) for p,i,sha in db.execute("SELECT path,source_identity,sha256 FROM operation_files")}
    selected = root / next(iter(before))
    external = tmp_path / "external"
    os.link(selected, external)
    with pytest.raises(ValueError, match="snapshot changed"):
        runtime.retire(root, expected)
    assert selected.exists()
    external.unlink()
    selected.chmod(0o644)
    selected.write_bytes(b"changed\n")
    selected.chmod(0o444)
    with pytest.raises(ValueError):
        runtime.retire(root, expected)
    assert selected.exists()


def test_real_shell_entrypoint_snapshot_reuses_archived_run_and_rejects_drift(tmp_path, modules):
    metadata, runtime = modules
    workflow = tmp_path / "workflow"
    workflow.mkdir()
    core = workflow / "gg_gene_evolution_core.sh"
    core.write_text("echo science\n")
    (workflow / "gg_gene_evolution_entrypoint.sh").write_text("entrypoint\n")
    (workflow / "gg_common_params.sh").write_text("parameters\n")
    (tmp_path / "VERSION").write_text("fixture\n")
    output = tmp_path / "workspace/output"
    prefix = (f"source {shlex.quote(str(SUPPORT / 'gg_util.sh'))}; "
              f"gg_workspace_output_dir={shlex.quote(str(output))}; gg_workflow_dir={shlex.quote(str(workflow))}; "
              f"GG_JOB_ID=123; GG_ARRAY_TASK_ID=1; GG_SCHEDULER_KIND=gridengine; GG_BOOTSTRAP_PYTHON={shlex.quote(sys.executable)}; ")
    cmd = prefix + "gg_prepare_entrypoint_runtime_snapshot gg_gene_evolution_entrypoint.sh " + shlex.quote(str(core))
    result = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    captured = Path(result.stdout.strip())
    original = captured.read_bytes()
    metadata.pack(output / "runtime", delete_originals=True, writers_stopped=True, operation_id="shell")
    runtime.prune(output / "runtime", writers_stopped=True)
    core.write_text("changed source\n")
    result = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    assert result.returncode != 0 and not captured.exists()
    core.write_bytes(original)
    result = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert Path(result.stdout.strip()) == captured and captured.read_bytes() == original


@pytest.mark.parametrize("failure", ["body", "reference"])
def test_runtime_body_corruption_never_becomes_source_or_resume_success(tmp_path, modules, failure):
    metadata, runtime = modules
    root, before = snapshots(tmp_path, runs=1)
    metadata.pack(root, delete_originals=True, writers_stopped=True, operation_id="corrupt")
    with sqlite3.connect(root / metadata.DATABASE) as db:
        if failure == "body":
            db.execute("UPDATE runtime_bodies SET body=x'00'")
        else:
            db.execute("DELETE FROM runtime_bodies")
    with pytest.raises((ValueError, __import__("zlib").error)):
        metadata.archived_bytes(root / next(iter(before)))
    with pytest.raises((ValueError, __import__("zlib").error)):
        metadata.pack(root, writers_stopped=True, operation_id="new")
