import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from workflow.support.shared_namespace_lock import (
    NamespaceLockError,
    acquire,
    namespace_lock,
    release,
)

HELPER = Path(__file__).resolve().parents[1] / "support" / "shared_namespace_lock.py"
SHELL = HELPER.with_name("gg_shared_lock.sh")


def test_readers_share_but_exclude_writer_across_processes(tmp_path):
    path = tmp_path / "family.lock"
    tokens = []
    for _ in range(2):
        result = subprocess.run([sys.executable, str(HELPER), "acquire-shared", str(path)],
                                capture_output=True, text=True, check=True)
        tokens.append(result.stdout.strip())
    assert acquire(path, exclusive=True, nonblocking=True) is None
    release(path, tokens[0], exclusive=False)
    assert acquire(path, exclusive=True, nonblocking=True) is None
    release(path, tokens[1], exclusive=False)
    with namespace_lock(path, exclusive=True) as writer:
        assert writer
        assert acquire(path, exclusive=False, nonblocking=True) is None
        assert acquire(path, exclusive=True, nonblocking=True) is None


def test_abandoned_remote_owner_is_not_stolen_by_age(tmp_path):
    path = tmp_path / "family.lock"
    token = acquire(path, exclusive=False)
    owner = Path(str(path) + ".namespace-v1/readers") / token
    payload = json.loads(owner.read_text())
    payload.update(pid=99999999, host="another-node", created_ns=1)
    owner.write_text(json.dumps(payload))
    os.utime(owner, (1, 1))
    with pytest.raises(NamespaceLockError, match="Timed out"):
        acquire(path, exclusive=True, timeout=0.01)
    assert owner.exists()
    release(path, token, exclusive=False)


def test_abandoned_gate_is_fail_closed(tmp_path):
    path = tmp_path / "family.lock"
    token = acquire(path, exclusive=True)
    with pytest.raises(NamespaceLockError, match="Timed out"):
        acquire(path, exclusive=False, timeout=0.01)
    with pytest.raises(NamespaceLockError, match="ownership changed"):
        release(path, "0" * 32, exclusive=True)
    release(path, token, exclusive=True)


def test_nested_readers_and_exception_release(tmp_path):
    path = tmp_path / "family.lock"
    with pytest.raises(ValueError):
        with namespace_lock(path, exclusive=False):
            with namespace_lock(path, exclusive=False):
                raise ValueError("body failure")
    with namespace_lock(path, exclusive=True, nonblocking=True) as acquired:
        assert acquired


@pytest.mark.parametrize("part", ["marker", "root", "readers", "gate"])
def test_symlinked_coordination_state_is_rejected(tmp_path, part):
    path = tmp_path / "family.lock"
    target = tmp_path / "elsewhere"
    target.mkdir()
    root = Path(str(path) + ".namespace-v1")
    if part == "marker":
        path.symlink_to(target)
    elif part == "root":
        root.symlink_to(target)
    else:
        root.mkdir()
        (root / part).symlink_to(target)
    with pytest.raises((NamespaceLockError, OSError)):
        acquire(path, exclusive=True, nonblocking=True)
    assert not list(target.iterdir())


def test_shell_holds_same_protocol_until_explicit_release(tmp_path):
    path = tmp_path / "family.lock"
    script = f'''set -euo pipefail
source "{SHELL}"
gg_advisory_shared_lock_acquire "$1"
printf 'ready\\n'
read -r reply
gg_advisory_shared_lock_release
'''
    process = subprocess.Popen(["bash", "-c", script, "test", str(path)],
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, text=True)
    try:
        assert process.stdout.readline() == "ready\n"
        assert acquire(path, exclusive=True, nonblocking=True) is None
        _, stderr = process.communicate("release\n", timeout=10)
        assert process.returncode == 0, stderr
        with namespace_lock(path, exclusive=True, nonblocking=True) as writer:
            assert writer
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()
