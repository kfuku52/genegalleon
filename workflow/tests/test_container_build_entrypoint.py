from pathlib import Path
import os
import stat
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[2]
ENTRYPOINT = REPO_ROOT / "gg_container_build_entrypoint.sh"
CHECK_ENV_COVERAGE = REPO_ROOT / "container" / "scripts" / "check_env_coverage.sh"


def _write_executable(path: Path, body: str) -> None:
    path.write_text(body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _prepare_stubbed_path(tmp_path: Path, runtime_name: str) -> tuple[str, Path]:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    runtime_log = tmp_path / "runtime.log"

    _write_executable(
        bin_dir / runtime_name,
        f"""#!/usr/bin/env bash
set -euo pipefail
printf '%s\\n' "$*" >> "{runtime_log}"
if [[ "$1" == "build" ]]; then
  touch "$2"
fi
""",
    )
    _write_executable(
        bin_dir / "docker",
        """#!/usr/bin/env bash
exit 127
""",
    )
    return f"{bin_dir}{os.pathsep}{os.environ['PATH']}", runtime_log


def _prepare_stubbed_path_with_buildx(tmp_path: Path, runtime_name: str) -> tuple[str, Path, Path]:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    runtime_log = tmp_path / "runtime.log"
    docker_log = tmp_path / "docker.log"

    _write_executable(
        bin_dir / runtime_name,
        f"""#!/usr/bin/env bash
set -euo pipefail
printf '%s\\n' "$*" >> "{runtime_log}"
if [[ "$1" == "build" ]]; then
  touch "$2"
fi
""",
    )
    _write_executable(
        bin_dir / "docker",
        f"""#!/usr/bin/env bash
set -euo pipefail
printf '%s\\n' "$*" >> "{docker_log}"
if [[ "$1" == "buildx" && "$2" == "ls" ]]; then
  exit 0
fi
if [[ "$1" == "buildx" && "$2" == "inspect" ]]; then
  printf '%s\\n' 'Driver: docker'
  exit 0
fi
if [[ "$1" == "buildx" && "$2" == "build" ]]; then
  exit 0
fi
exit 0
""",
    )
    return f"{bin_dir}{os.pathsep}{os.environ['PATH']}", runtime_log, docker_log


def _run_entrypoint(tmp_path: Path, runtime_name: str, extra_env: dict[str, str] | None = None):
    path_value, runtime_log = _prepare_stubbed_path(tmp_path, runtime_name)
    env = os.environ.copy()
    env["PATH"] = path_value
    env["OUT"] = str(tmp_path / "genegalleon.sif")
    env["BUILD_SIF"] = "1"
    if extra_env:
        env.update(extra_env)
    completed = subprocess.run(
        ["bash", str(ENTRYPOINT)],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )
    return completed, runtime_log, Path(env["OUT"])


def _run_entrypoint_with_buildx(tmp_path: Path, runtime_name: str, extra_env: dict[str, str] | None = None):
    path_value, runtime_log, docker_log = _prepare_stubbed_path_with_buildx(tmp_path, runtime_name)
    env = os.environ.copy()
    env["PATH"] = path_value
    env["OUT"] = str(tmp_path / "genegalleon.sif")
    env["BUILD_SIF"] = "1"
    if extra_env:
        env.update(extra_env)
    completed = subprocess.run(
        ["bash", str(ENTRYPOINT)],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )
    return completed, runtime_log, docker_log, Path(env["OUT"])


def test_container_build_entrypoint_uses_public_image_without_docker_with_apptainer(tmp_path: Path):
    completed, runtime_log, out_path = _run_entrypoint(
        tmp_path,
        "apptainer",
        {
            "IMAGE_SOURCE": "public",
            "IMAGE": "ghcr.io/example/genegalleon",
            "TAG": "20260306-test",
        },
    )

    assert completed.returncode == 0, completed.stderr
    assert out_path.is_file()
    assert "[gg_container_build] image_source=public" in completed.stdout
    assert "[gg_container_build] engine=apptainer" in completed.stdout
    assert "[gg_container_build] step 1/1: build SIF from registry image" in completed.stdout
    assert runtime_log.read_text(encoding="utf-8").strip().splitlines() == [
        f"build {out_path} docker://ghcr.io/example/genegalleon:20260306-test"
    ]


def test_container_build_entrypoint_auto_detects_singularity_without_docker(tmp_path: Path):
    completed, runtime_log, out_path = _run_entrypoint(
        tmp_path,
        "singularity",
        {
            "IMAGE_SOURCE": "public",
            "IMAGE": "ghcr.io/example/genegalleon",
            "TAG": "latest",
        },
    )

    assert completed.returncode == 0, completed.stderr
    assert out_path.is_file()
    assert "[gg_container_build] engine=singularity" in completed.stdout
    assert runtime_log.read_text(encoding="utf-8").strip().splitlines() == [
        f"build {out_path} docker://ghcr.io/example/genegalleon:latest"
    ]


def test_container_build_entrypoint_falls_back_to_official_registry_image_when_defaults_are_used(tmp_path: Path):
    completed, runtime_log, out_path = _run_entrypoint(tmp_path, "apptainer")

    assert completed.returncode == 0, completed.stderr
    assert out_path.is_file()
    assert "[gg_container_build] image_source=auto" in completed.stdout
    assert "falling back to published image ghcr.io/kfuku52/genegalleon:latest" in completed.stdout
    assert runtime_log.read_text(encoding="utf-8").strip().splitlines() == [
        f"build {out_path} docker://ghcr.io/kfuku52/genegalleon:latest"
    ]


def test_container_build_entrypoint_falls_back_to_latest_when_local_default_tag_is_explicit(tmp_path: Path):
    completed, runtime_log, out_path = _run_entrypoint(
        tmp_path,
        "apptainer",
        {
            "TAG": "dev",
        },
    )

    assert completed.returncode == 0, completed.stderr
    assert out_path.is_file()
    assert "falling back to published image ghcr.io/kfuku52/genegalleon:latest" in completed.stdout
    assert runtime_log.read_text(encoding="utf-8").strip().splitlines() == [
        f"build {out_path} docker://ghcr.io/kfuku52/genegalleon:latest"
    ]


def test_container_build_entrypoint_prefers_public_pull_for_remote_image_in_auto_mode_even_with_buildx(tmp_path: Path):
    completed, runtime_log, docker_log, out_path = _run_entrypoint_with_buildx(
        tmp_path,
        "apptainer",
        {
            "IMAGE": "ghcr.io/example/genegalleon",
            "TAG": "20260306-test",
        },
    )

    assert completed.returncode == 0, completed.stderr
    assert out_path.is_file()
    assert "[gg_container_build] resolved_image_source=public" in completed.stdout
    assert "skipping local Docker build" in completed.stdout
    assert runtime_log.read_text(encoding="utf-8").strip().splitlines() == [
        f"build {out_path} docker://ghcr.io/example/genegalleon:20260306-test"
    ]
    assert not docker_log.exists()


def test_container_build_entrypoint_uses_native_local_build_without_docker(tmp_path: Path):
    completed, runtime_log, out_path = _run_entrypoint(
        tmp_path,
        "apptainer",
        {
            "IMAGE_SOURCE": "local",
            "IMAGE": "local/genegalleon",
            "TAG": "dev",
        },
    )

    assert completed.returncode == 0, completed.stderr
    assert out_path.is_file()
    assert "step 1/1: native local build from repository" in completed.stdout
    assert "[apptainer_local_build] definition=" in completed.stdout
    runtime_log_lines = runtime_log.read_text(encoding="utf-8").strip().splitlines()
    assert len(runtime_log_lines) == 1
    assert runtime_log_lines[0].startswith(f"build {out_path} ")
    assert runtime_log_lines[0].endswith(".def")


def test_container_build_entrypoint_rejects_non_registry_image_for_public_source(tmp_path: Path):
    completed, runtime_log, out_path = _run_entrypoint(
        tmp_path,
        "apptainer",
        {
            "IMAGE_SOURCE": "public",
            "IMAGE": "local/genegalleon",
            "TAG": "dev",
        },
    )

    assert completed.returncode != 0
    assert not out_path.exists()
    assert not runtime_log.exists()
    assert "IMAGE_SOURCE=public requires a registry image reference." in completed.stdout


def test_container_env_coverage_preflight_matches_current_repo():
    completed = subprocess.run(
        ["bash", str(CHECK_ENV_COVERAGE), str(REPO_ROOT)],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "Conda environment coverage matches" in completed.stdout
