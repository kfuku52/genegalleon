import json
import os
import signal
import subprocess
import sys
import textwrap
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
GG_UTIL = REPO_ROOT / "workflow" / "support" / "gg_util.sh"
SUPPORT_DIR = REPO_ROOT / "workflow" / "support"
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))

from shared_lock import FIELD_SEPARATOR, read_lock_metadata, try_create_lock  # noqa: E402


def test_python_and_shell_shared_lock_helpers_use_same_metadata_format(tmp_path):
    lock_path = tmp_path / "shared.lock"
    assert try_create_lock(lock_path, pid=12345, hostname="lock-host", boot_id="boot-1")

    script = textwrap.dedent(
        f"""
        set -euo pipefail
        source "{GG_UTIL}"
        gg_shared_lock_read_metadata "{lock_path}"
        """
    )
    completed = subprocess.run(
        ["bash", "-lc", script],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    fields = completed.stdout.strip().split(FIELD_SEPARATOR)
    assert fields[:4] == ["shared-lock-v2", "12345", "lock-host", "boot-1"]
    assert fields[5] != ""
    assert fields[6] != ""
    assert fields[7] != ""

    metadata = read_lock_metadata(lock_path)
    assert metadata["format"] == "shared-lock-v2"
    assert metadata["pid"] == 12345
    assert metadata["hostname"] == "lock-host"
    assert metadata["boot_id"] == "boot-1"


def test_shell_shared_lock_try_create_writes_python_readable_metadata(tmp_path):
    lock_path = tmp_path / "created-by-shell.lock"
    script = textwrap.dedent(
        f"""
        set -euo pipefail
        source "{GG_UTIL}"
        gg_shared_lock_try_create "{lock_path}" 4242
        """
    )
    completed = subprocess.run(
        ["bash", "-lc", script],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    payload = json.loads(lock_path.read_text(encoding="utf-8"))
    assert payload["format"] == "shared-lock-v2"
    assert payload["pid"] == 4242
    assert payload["hostname"]
    assert "created_at" in payload


def test_malformed_shared_lock_still_uses_stat_fields_for_stale_recovery(tmp_path):
    lock_path = tmp_path / "uniref90.lock"
    script = textwrap.dedent(
        f"""
        set -euo pipefail
        source "{GG_UTIL}"
        lock_path="{lock_path}"
        python - "$lock_path" <<'PY'
import os
import sys
import time

path = sys.argv[1]
with open(path, "w", encoding="utf-8") as handle:
    handle.write("legacy-lock-format\\n")
old = time.time() - 10
os.utime(path, (old, old))
PY
        summary="$(gg_shared_lock_owner_summary "$lock_path")"
        printf '%s\\n' "$summary"
        gg_shared_lock_reclaim_if_stale "$lock_path" "MMseqs2 UniRef90 taxonomy DB"
        [[ ! -e "$lock_path" ]]
        """
    )
    completed = subprocess.run(
        ["bash", "-lc", script],
        cwd=REPO_ROOT,
        env={**os.environ, "GG_LOCK_STALE_SECONDS": "1"},
        text=True,
        capture_output=True,
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert "pid=unknown" in completed.stdout
    assert "host=unknown" in completed.stdout
    assert "heartbeat_age=" in completed.stdout
    assert "heartbeat_age=unknown" not in completed.stdout
    assert "Recovered stale shared lock" in completed.stderr


def test_shared_semaphore_limits_concurrency_and_releases_slots(tmp_path):
    semaphore_dir = tmp_path / "ncbi_semaphore"
    script = textwrap.dedent(
        f"""
        set -euo pipefail
        source "{GG_UTIL}"
        semaphore_dir="{semaphore_dir}"

        gg_shared_semaphore_acquire "$semaphore_dir" 2 "NCBI test"
        slot1="$GG_SHARED_SEMAPHORE_SLOT_LOCK_FILE"
        idx1="$GG_SHARED_SEMAPHORE_SLOT_INDEX"

        gg_shared_semaphore_acquire "$semaphore_dir" 2 "NCBI test"
        slot2="$GG_SHARED_SEMAPHORE_SLOT_LOCK_FILE"
        idx2="$GG_SHARED_SEMAPHORE_SLOT_INDEX"

        [[ -n "$slot1" ]]
        [[ -n "$slot2" ]]
        [[ "$slot1" != "$slot2" ]]
        [[ "$idx1" != "$idx2" ]]

        if gg_shared_semaphore_acquire "$semaphore_dir" 2 "NCBI test"; then
          echo "expected semaphore acquisition to time out" >&2
          exit 1
        fi

        gg_shared_semaphore_release "$slot1"
        gg_shared_semaphore_release "$slot2"

        gg_shared_semaphore_acquire "$semaphore_dir" 2 "NCBI test"
        slot3="$GG_SHARED_SEMAPHORE_SLOT_LOCK_FILE"
        [[ -n "$slot3" ]]
        gg_shared_semaphore_release "$slot3"
        """
    )
    completed = subprocess.run(
        ["bash", "-lc", script],
        cwd=REPO_ROOT,
        env={
            **os.environ,
            "GG_LOCK_ACQUIRE_TIMEOUT_SECONDS": "1",
            "GG_LOCK_POLL_SECONDS": "1",
        },
        text=True,
        capture_output=True,
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert "timed out waiting for shared semaphore slot" in completed.stderr


def test_shared_semaphore_releases_slot_immediately_on_term(tmp_path):
    semaphore_dir = tmp_path / "ncbi_semaphore_term"
    script = textwrap.dedent(
        f"""
        set -euo pipefail
        source "{GG_UTIL}"
        gg_run_with_shared_semaphore "{semaphore_dir}" 1 "NCBI test" bash -lc "sleep 30"
        """
    )
    completed = subprocess.Popen(
        ["bash", "-c", script],
        cwd=REPO_ROOT,
        env={**os.environ, "GG_LOCK_STALE_SECONDS": "900"},
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
    )
    lock_path = semaphore_dir / "slot.1.lock"
    try:
        for _ in range(50):
            if lock_path.exists():
                break
            time.sleep(0.1)
        assert lock_path.exists()

        os.killpg(completed.pid, signal.SIGTERM)
        completed.wait(timeout=5)

        for _ in range(50):
            if not lock_path.exists():
                break
            time.sleep(0.1)
        assert not lock_path.exists()
    finally:
        if completed.poll() is None:
            os.killpg(completed.pid, signal.SIGKILL)
            completed.wait(timeout=5)
