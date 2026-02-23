from pathlib import Path
import os
import subprocess
import sys

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SUPPORT_DIR = REPO_ROOT / "workflow" / "support"

SKIP_HELP_SCRIPTS = set()
SMOKE_HELP_SCRIPTS = sorted(
    script.name
    for script in SUPPORT_DIR.glob("*.py")
    if script.name not in SKIP_HELP_SCRIPTS
)

REQUIRES_KFTOOLS_STUB = {"iqtree2mapnh.py", "orthogroup_statistics.py"}


def _write_kftools_stub(stub_root: Path) -> None:
    pkg_dir = stub_root / "kftools"
    pkg_dir.mkdir(parents=True, exist_ok=True)
    (pkg_dir / "__init__.py").write_text("", encoding="utf-8")
    for module_name in ("kfphylo.py", "kfseq.py", "kfog.py"):
        (pkg_dir / module_name).write_text("__all__ = []\n", encoding="utf-8")


@pytest.mark.parametrize("script_name", SMOKE_HELP_SCRIPTS)
def test_support_script_help_smoke(script_name: str, tmp_path: Path):
    script_path = SUPPORT_DIR / script_name
    assert script_path.exists(), f"Support script not found: {script_path}"

    run_cwd = tmp_path / "run_cwd"
    run_cwd.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    mpl_config_dir = run_cwd / "mplconfig"
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    env["MPLCONFIGDIR"] = str(mpl_config_dir)
    if script_name in REQUIRES_KFTOOLS_STUB:
        stub_root = tmp_path / "stub_pkgs"
        _write_kftools_stub(stub_root)
        existing_pythonpath = env.get("PYTHONPATH", "")
        env["PYTHONPATH"] = (
            f"{stub_root}{os.pathsep}{existing_pythonpath}"
            if existing_pythonpath
            else str(stub_root)
        )

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
    assert created == [], f"--help created unexpected files for {script_name}: {created}"
