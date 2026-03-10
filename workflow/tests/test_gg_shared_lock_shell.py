import os
import subprocess
import textwrap
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
GG_UTIL = REPO_ROOT / "workflow" / "support" / "gg_util.sh"


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
