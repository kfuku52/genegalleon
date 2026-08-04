import os
import subprocess
from pathlib import Path

from shell_static_helpers import REPO_ROOT


RUN_IN_SIF = REPO_ROOT / "workflow" / "tests" / "run_in_sif.sh"


def _fake_runtime(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "#!/usr/bin/env bash\n"
        "printf '%s\\n' \"$@\" > \"${GG_CAPTURE_ARGS}\"\n",
        encoding="utf-8",
    )
    path.chmod(0o755)


def test_run_in_sif_discovers_versioned_runtime_and_adds_read_only_bind(tmp_path):
    package_root = tmp_path / "packages"
    runtime = package_root / "apptainer" / "1.4.5" / "bin" / "apptainer"
    _fake_runtime(runtime)
    sif = tmp_path / "genegalleon.sif"
    sif.touch()
    extra = tmp_path / "external-project"
    extra.mkdir()
    captured = tmp_path / "runtime-args.txt"
    env = os.environ.copy()
    env.update(
        {
            "PATH": "/usr/bin:/bin",
            "GENEGALLEON_SIF": str(sif),
            "GENEGALLEON_SIF_EXTRA_BINDS": str(extra),
            "GG_CONTAINER_RUNTIME_PACKAGE_ROOT": str(package_root),
            "GG_CONTAINER_RUNTIME_LEGACY_DIR": str(tmp_path / "legacy"),
            "GG_CAPTURE_ARGS": str(captured),
        }
    )

    completed = subprocess.run(
        ["/bin/bash", str(RUN_IN_SIF), "python", "--version"],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    args = captured.read_text(encoding="utf-8").splitlines()
    assert args[:2] == ["exec", "--cleanenv"]
    assert f"{REPO_ROOT}:{REPO_ROOT}" in args
    assert f"{extra}:{extra}:ro" in args
    assert args[-3:] == [str(sif), "python", "--version"]


def test_run_in_sif_rejects_relative_extra_bind(tmp_path):
    runtime_dir = tmp_path / "bin"
    _fake_runtime(runtime_dir / "apptainer")
    sif = tmp_path / "genegalleon.sif"
    sif.touch()
    env = os.environ.copy()
    env.update(
        {
            "PATH": f"{runtime_dir}:/usr/bin:/bin",
            "GENEGALLEON_SIF": str(sif),
            "GENEGALLEON_SIF_EXTRA_BINDS": "relative/path",
            "GG_CAPTURE_ARGS": str(tmp_path / "runtime-args.txt"),
        }
    )

    completed = subprocess.run(
        ["/bin/bash", str(RUN_IN_SIF), "true"],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 1
    assert "Invalid GENEGALLEON_SIF_EXTRA_BINDS path" in completed.stderr
