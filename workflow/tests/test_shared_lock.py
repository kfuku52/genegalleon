import os
import subprocess
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))

import shared_lock  # noqa: E402


def test_snapshot_cannot_mix_metadata_and_inode_during_owner_handoff(tmp_path, monkeypatch):
    lock = tmp_path / "shared.lock"
    old_token = shared_lock.write_lock_file(lock, pid=99999999)
    old_stat = lock.stat()
    original_fstat = os.fstat
    replaced = False

    def replace_before_fstat(fd):
        nonlocal replaced
        if not replaced:
            replaced = True
            lock.unlink()
            shared_lock.write_lock_file(lock)
        return original_fstat(fd)

    monkeypatch.setattr(os, "fstat", replace_before_fstat)
    metadata, stat_result = shared_lock.read_lock_snapshot(lock)
    assert metadata["token"] == old_token
    assert (stat_result.st_dev, stat_result.st_ino) == (old_stat.st_dev, old_stat.st_ino)
    assert shared_lock.read_lock_metadata(lock)["token"] != old_token
    assert shared_lock.reclaim_if_stale(lock, 86400) is None
    assert lock.exists()


def test_previous_owner_cannot_release_or_heartbeat_replacement(tmp_path):
    lock = tmp_path / "shared.lock"
    previous = shared_lock.acquire_lock(lock, 1, 1, 0.01, "previous")
    os.utime(lock, (0, 0))
    replacement = shared_lock.acquire_lock(lock, 1, 1, 0.01, "replacement")
    before = lock.stat().st_mtime_ns
    assert previous.token != replacement.token
    assert not shared_lock.update_owned_lock(lock, previous.token)
    assert lock.stat().st_mtime_ns == before
    assert not shared_lock.release_lock(lock, previous)
    assert shared_lock.read_lock_metadata(lock)["token"] == replacement.token
    assert shared_lock.release_lock(lock, replacement)
    assert not lock.exists()


def test_mutation_guard_blocks_other_process_until_handoff_finishes(tmp_path):
    lock = tmp_path / "shared.lock"
    script = """
import sys
from shared_lock import try_create_lock
print('attempting', flush=True)
assert try_create_lock(sys.argv[1])
print('acquired', flush=True)
"""
    with shared_lock.lock_mutation_guard(lock):
        child = subprocess.Popen(
            [sys.executable, "-c", script, str(lock)],
            env={**os.environ, "PYTHONPATH": str(SUPPORT)},
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        assert child.stdout.readline().strip() == "attempting"
        with pytest.raises(subprocess.TimeoutExpired):
            child.communicate(timeout=0.2)
        assert not lock.exists()
    stdout, stderr = child.communicate(timeout=5)
    assert child.returncode == 0, stderr
    assert "acquired" in stdout
    assert lock.exists()


def test_shell_old_owner_preserves_replacement_and_stops_heartbeat(tmp_path):
    lock = tmp_path / "shared.lock"
    script = r'''
set -euo pipefail
source "$1/gg_shared_lock.sh"
gg_shared_lock_acquire "$2" audit
gg_shared_lock_start_heartbeat "$2"
heartbeat_pid=$GG_SHARED_LOCK_HEARTBEAT_PID
python - "$1" "$2" "$$" <<'PY'
import os
import sys
sys.path.insert(0, sys.argv[1])
from shared_lock import reclaim_if_stale, write_lock_file
os.utime(sys.argv[2], (0, 0))
assert reclaim_if_stale(sys.argv[2], 1)
write_lock_file(sys.argv[2], pid=int(sys.argv[3]), token='replacement')
os.utime(sys.argv[2], (1, 1))
PY
sleep 1.2
gg_shared_lock_stop_heartbeat "$heartbeat_pid"
gg_shared_lock_release "$2"
'''
    result = subprocess.run(
        ["bash", "-c", script, "audit", str(SUPPORT), str(lock)],
        env={**os.environ, "GG_LOCK_HEARTBEAT_SECONDS": "1"},
        capture_output=True, text=True, timeout=15,
    )
    assert result.returncode == 0, result.stderr
    assert shared_lock.read_lock_metadata(lock)["token"] == "replacement"
    assert lock.stat().st_mtime == 1


def test_shell_reports_guard_errors_without_retrying_as_contention(tmp_path):
    lock = tmp_path / "shared.lock"
    unrelated = tmp_path / "unrelated"
    unrelated.write_text("preserve")
    (tmp_path / ".shared.lock.guard").symlink_to(unrelated)
    result = subprocess.run(
        ["bash", "-c", 'source "$1/gg_shared_lock.sh"; gg_shared_lock_acquire "$2" audit',
         "audit", str(SUPPORT), str(lock)],
        capture_output=True, text=True, timeout=5,
    )
    assert result.returncode == 2, result.stderr
    assert "Shared lock operation failed" in result.stderr
    assert not lock.exists()
    assert unrelated.read_text() == "preserve"
