import shlex
import subprocess
from pathlib import Path

from shell_static_helpers import WORKFLOW_DIR

UTIL = WORKFLOW_DIR / "support" / "gg_util.sh"


def run_site_command(tmp_path: Path, bind: str) -> subprocess.CompletedProcess[str]:
    command = (
        f"source {shlex.quote(str(UTIL))}; "
        "hostname() { printf 'audrey1\\n'; }; "
        f"export GG_CONTAINER_PROJECT_ROOT_BIND={shlex.quote(bind)}; "
        "gg_site_container_shell_command singularity singularity_command; "
        "printf 'command=%s\\n' \"${singularity_command[*]}\"; "
        "printf 'bind=%s\\n' \"${SINGULARITY_BINDPATH:-}\""
    )
    return subprocess.run(
        ["bash", "-euo", "pipefail", "-c", command],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )


def test_audrey1_isolates_shm_and_keeps_project_paths_visible(tmp_path):
    root = tmp_path / "project"
    root.mkdir()
    result = run_site_command(tmp_path, f"{root}:{root}")
    assert result.returncode == 0, result.stderr
    assert "site profile = audrey1" in result.stdout
    assert "command=singularity exec --contain" in result.stdout
    assert f"bind={root}:{root}" in result.stdout


def test_audrey1_rejects_unreviewed_bind(tmp_path):
    root = tmp_path / "project"
    root.mkdir()
    result = run_site_command(tmp_path, f"{root}:/unrelated")
    assert result.returncode != 0
    assert "invalid GG_CONTAINER_PROJECT_ROOT_BIND" in result.stderr


def test_audrey1_existing_unbound_runs_keep_previous_runtime(tmp_path):
    result = run_site_command(tmp_path, "")
    assert result.returncode == 0, result.stderr
    assert "command=singularity exec\n" in result.stdout
