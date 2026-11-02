"""Bounded provenance ownership and byte-identical runtime snapshot sharing."""
import argparse
import fcntl
import importlib
import shlex
import subprocess
from pathlib import Path

SUPPORT = Path(__file__).resolve().parents[1] / "support"
UTIL = SUPPORT / "gg_util.sh"


def shell(source, destination, pool, *, replace=False):
    command = "source " + shlex.quote(str(UTIL)) + "; gg_share_runtime_file "
    command += shlex.join([str(source), str(destination), str(pool), str(int(replace))])
    return subprocess.run(["bash", "-c", command], capture_output=True, text=True)


def test_runtime_files_share_inodes_without_linking_the_mutable_source(tmp_path):
    source = tmp_path / "core.sh"
    source.write_bytes(b"echo original\n")
    source.chmod(0o755)
    pool = tmp_path / "pool"
    first, second = tmp_path / "first", tmp_path / "second"
    for target in (first, second):
        result = shell(source, target, pool)
        assert result.returncode == 0, result.stderr
    assert first.stat().st_ino == second.stat().st_ino
    assert first.stat().st_ino != source.stat().st_ino
    assert first.stat().st_mode & 0o222 == 0
    source.write_bytes(b"echo new source\n")
    assert first.read_bytes() == second.read_bytes() == b"echo original\n"
    result = shell(source, first, pool)
    assert result.returncode != 0
    assert first.read_bytes() == b"echo original\n"
    assert not list(pool.glob(".publish.*"))


def test_existing_runtime_copy_is_retained_until_explicit_migration(tmp_path):
    source, target, pool = tmp_path / "source", tmp_path / "target", tmp_path / "pool"
    source.write_bytes(b"immutable captured bytes\n")
    target.write_bytes(source.read_bytes())
    inode = target.stat().st_ino
    result = shell(source, target, pool)
    assert result.returncode == 0, result.stderr
    assert target.stat().st_ino == inode
    result = shell(target, target, pool, replace=True)
    assert result.returncode == 0, result.stderr
    shared_inode = target.stat().st_ino
    assert shared_inode != inode
    result = shell(target, target, pool, replace=True)
    assert result.returncode == 0, result.stderr
    assert target.stat().st_ino == shared_inode
    assert len(list(pool.glob("[0-9a-f]*"))) == 1


def test_snapshot_sharing_rejects_redirected_and_writable_pool_objects(tmp_path):
    source, target, pool = tmp_path / "source", tmp_path / "target", tmp_path / "pool"
    source.write_bytes(b"source\n")
    assert shell(source, target, pool).returncode == 0
    pooled = next(pool.glob("[0-9a-f]*"))
    pooled.chmod(0o644)
    assert shell(source, tmp_path / "other", pool).returncode != 0
    target.unlink()
    target.symlink_to(source)
    assert shell(source, target, pool).returncode != 0
    assert source.read_bytes() == b"source\n"


def test_provenance_ownership_is_bounded_and_independent_of_flock(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(SUPPORT))
    module = importlib.import_module("artifact_provenance")
    monkeypatch.setattr(fcntl, "flock", lambda *_: None)
    for index in range(1024):
        args = argparse.Namespace(workspace_root=tmp_path, manifest=tmp_path / f"{index}.json", dry_run=False)
        with module.artifact_manifest_lock(args):
            pass
    root = tmp_path / ".gg_cache" / "artifact_provenance_lock_stripes-v1"
    assert len(list(root.glob("*.lock"))) <= 256
    assert len(list(root.rglob("*"))) <= 256 * 3
    assert not list(root.rglob("owner.json"))
    assert not (tmp_path / ".gg_cache" / "artifact_provenance_locks").exists()


def test_dry_run_provenance_does_not_create_lock_namespace(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(SUPPORT))
    module = importlib.import_module("artifact_provenance")
    args = argparse.Namespace(workspace_root=tmp_path, manifest=tmp_path / "m.json", dry_run=True)
    with module.artifact_manifest_lock(args):
        pass
    assert not (tmp_path / ".gg_cache").exists()
