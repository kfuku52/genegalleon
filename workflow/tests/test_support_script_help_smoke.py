import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SUPPORT_DIR = REPO_ROOT / "workflow" / "support"

# Library modules have no CLI: running them with --help only repeats imports.
# Keep process-level checks for actual entrypoints, where parser construction
# can fail before any functional test reaches the command.
# These commands require kftools, which is not in the fast-lane environment.
# Their dedicated statistics/tree tests cover behavior; do not fake the import
# just to make --help succeed.
RUNTIME_HELP_SCRIPTS = {"iqtree2mapnh.py", "orthogroup_statistics.py"}
SMOKE_HELP_SCRIPTS = sorted(
    script.name
    for script in SUPPORT_DIR.glob("*.py")
    if script.name not in RUNTIME_HELP_SCRIPTS
    and re.search(r"if\s+__name__\s*==\s*['\"]__main__['\"]", script.read_text(encoding="utf-8"))
)


@pytest.mark.parametrize("script_name", SMOKE_HELP_SCRIPTS)
def test_support_script_help_smoke(script_name: str, tmp_path: Path):
    script_path = SUPPORT_DIR / script_name

    run_cwd = tmp_path / "run_cwd"
    run_cwd.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    mpl_config_dir = run_cwd / "mplconfig"
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    env["MPLCONFIGDIR"] = str(mpl_config_dir)

    before = {p.relative_to(run_cwd).as_posix() for p in run_cwd.rglob("*")}
    proc = subprocess.run(
        [sys.executable, str(script_path), "--help"],
        cwd=str(run_cwd),
        capture_output=True,
        text=True,
        env=env,
    )
    after = {p.relative_to(run_cwd).as_posix() for p in run_cwd.rglob("*")}
    created = sorted(after - before)
    created = [path for path in created if not path.startswith("mplconfig")]
    assert proc.returncode == 0, (
        f"--help failed for {script_name}\n"
        f"stdout:\n{proc.stdout}\n"
        f"stderr:\n{proc.stderr}"
    )
    assert "usage:" in proc.stdout.lower(), f"No CLI help for {script_name}: {proc.stdout}"
    assert created == [], f"--help created unexpected files for {script_name}: {created}"
