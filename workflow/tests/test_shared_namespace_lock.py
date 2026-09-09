import json
import os
import shlex
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


@pytest.mark.parametrize("ending,retained", [("exit 0", False), ("exit 7", False), ("exit 137", True), ("kill -TERM $$", True)])
def test_input_core_lock_cleanup_retains_interrupted_ownership(tmp_path, ending, retained):
    core = HELPER.parents[1] / "core" / "gg_input_generation_core.sh"
    text = core.read_text()
    # Exercise the actual lifecycle definitions, without unrelated workflow stages.
    lifecycle = text[text.index("array_lock_paths=()"):text.index("prepare_input_generation_tmp_dirs()")]
    path = tmp_path / "worker.lock"
    script = f"""
set -euo pipefail
gg_support_dir={shlex.quote(str(HELPER.parent))}
write_gg_input_generation_summary_on_exit() {{ return 0; }}
{lifecycle}
input_generation_lock {shlex.quote(str(path))} exclusive
{ending}
"""
    result = subprocess.run(["bash", "-c", script], capture_output=True, text=True)
    if retained:
        assert result.returncode >= 128, result.stderr
        assert acquire(path, exclusive=True, nonblocking=True) is None
    else:
        with namespace_lock(path, exclusive=True, nonblocking=True) as held:
            assert held, result.stderr


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


@pytest.mark.parametrize("nonblocking", [False, True])
def test_owner_releases_gate_between_mkdir_and_inspection(tmp_path, monkeypatch, nonblocking):
    path = tmp_path / "family.lock"
    owner = acquire(path, exclusive=True)
    gate = Path(str(path) + ".namespace-v1/gate")
    original_mkdir = Path.mkdir
    released = False

    def release_after_collision(directory, *args, **kwargs):
        nonlocal released
        try:
            return original_mkdir(directory, *args, **kwargs)
        except FileExistsError:
            if directory == gate and not released:
                release(path, owner, exclusive=True)
                released = True
            raise

    monkeypatch.setattr(Path, "mkdir", release_after_collision)
    contender = acquire(path, exclusive=True, nonblocking=nonblocking, timeout=1)
    assert released
    if nonblocking:
        assert contender is None
        contender = acquire(path, exclusive=True, nonblocking=True)
    assert contender is not None
    release(path, contender, exclusive=True)


def test_regular_file_at_gate_remains_an_error(tmp_path):
    path = tmp_path / "family.lock"
    with namespace_lock(path, exclusive=True):
        pass
    gate = Path(str(path) + ".namespace-v1/gate")
    gate.write_text("not a lock directory")
    with pytest.raises(NamespaceLockError, match="Invalid lock gate"):
        acquire(path, exclusive=True, nonblocking=True)
    assert gate.read_text() == "not a lock directory"


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


def test_input_core_shared_lock_retries_registration_contention(tmp_path):
    core = HELPER.parents[1] / "core" / "gg_input_generation_core.sh"
    text = core.read_text()
    lifecycle = text[text.index("array_lock_paths=()"):text.index("prepare_input_generation_tmp_dirs()")]
    path = tmp_path / "phase.lock"
    holder = acquire(path, exclusive=True)
    script = f"""
set -euo pipefail
gg_support_dir={shlex.quote(str(HELPER.parent))}
write_gg_input_generation_summary_on_exit() {{ return 0; }}
{lifecycle}
printf 'starting\\n'
input_generation_lock "$1" shared
printf 'acquired\\n'
"""
    process = subprocess.Popen(["bash", "-c", script, "test", str(path)],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        assert process.stdout.readline() == "starting\n"
        # The old nonblocking reader exits while the registration gate is held.
        with pytest.raises(subprocess.TimeoutExpired):
            process.wait(timeout=0.5)
    finally:
        release(path, holder, exclusive=True)
        stdout, stderr = process.communicate(timeout=10)
    assert process.returncode == 0, stdout + stderr
    assert stdout == "acquired\n"
    with namespace_lock(path, exclusive=True, nonblocking=True) as available:
        assert available
