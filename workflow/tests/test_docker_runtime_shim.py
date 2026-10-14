import os
import shlex
import stat
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_DIR = REPO_ROOT / "workflow"
SUPPORT_DIR = WORKFLOW_DIR / "support"
GG_UTIL_PATH = SUPPORT_DIR / "gg_util.sh"
SHIM_PATH = SUPPORT_DIR / "gg_wrapper_bin" / "singularity"
PROGRESS_ENTRYPOINT = WORKFLOW_DIR / "gg_progress_summary_entrypoint.sh"
REPO_VERSION = REPO_ROOT.joinpath("VERSION").read_text(encoding="utf-8").splitlines()[0].strip()
FIXTURE_IMAGE_ID = "sha256:" + "a" * 64


def _write_executable(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _run_bash(command: str, cwd: Path, env: dict[str, str] | None = None):
    run_env = os.environ.copy()
    if env:
        run_env.update(env)
    return subprocess.run(
        ["bash", "-c", command],
        cwd=cwd,
        capture_output=True,
        text=True,
        env=run_env,
        check=False,
    )


def _prepare_stub_docker(tmp_path: Path) -> Path:
    bin_dir = tmp_path / "bin"
    _write_executable(
        bin_dir / "docker",
        f"""#!/usr/bin/env bash
set -euo pipefail
base_dir={shlex.quote(str(tmp_path))}
if [[ "${{1:-}} ${{2:-}} ${{3:-}}" == 'image inspect --format' ]]; then
  printf '%s\\n' {shlex.quote(FIXTURE_IMAGE_ID)}
  exit 0
fi
counter_file="${{base_dir}}/docker.counter"
count=0
if [[ -f "${{counter_file}}" ]]; then
  count=$(cat "${{counter_file}}")
fi
count=$((count + 1))
printf '%s' "${{count}}" > "${{counter_file}}"
printf '%s\\n' "$@" > "${{base_dir}}/call_${{count}}.args"
cat > "${{base_dir}}/call_${{count}}.stdin"
""",
    )
    return bin_dir


def _prepare_stub_docker_with_known_images(tmp_path: Path, known_images: list[str]) -> Path:
    bin_dir = tmp_path / "bin"
    known_cases = "\n".join(
        f"    {shlex.quote(image)}) exit 0 ;;" for image in known_images
    )
    _write_executable(
        bin_dir / "docker",
        f"""#!/usr/bin/env bash
set -euo pipefail
base_dir={shlex.quote(str(tmp_path))}
if [[ "${{1:-}} ${{2:-}} ${{3:-}}" == 'image inspect --format' ]]; then
  printf '%s\\n' {shlex.quote(FIXTURE_IMAGE_ID)}
  exit 0
fi
counter_file="${{base_dir}}/docker.counter"
count=0
if [[ -f "${{counter_file}}" ]]; then
  count=$(cat "${{counter_file}}")
fi
count=$((count + 1))
printf '%s' "${{count}}" > "${{counter_file}}"
printf '%s\\n' "$@" > "${{base_dir}}/call_${{count}}.args"
if [[ "${{1:-}}" == "image" && "${{2:-}}" == "inspect" ]]; then
  case "${{3:-}}" in
{known_cases}
    *) exit 1 ;;
  esac
fi
cat > "${{base_dir}}/call_${{count}}.stdin"
""",
    )
    return bin_dir


def test_default_docker_shim_uses_checked_in_executable():
    completed = _run_bash(
        f"source {shlex.quote(str(GG_UTIL_PATH))}; unset GG_WRAPPER_BIN; gg_docker_singularity_shim_path",
        cwd=REPO_ROOT,
    )
    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(SHIM_PATH)


def test_detect_container_runtime_binary_returns_docker_shim(tmp_path: Path):
    wrapper_bin = tmp_path / "gg_wrapper_bin"
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"gg_support_dir={shlex.quote(str(SUPPORT_DIR))}; "
        f"GG_WRAPPER_BIN={shlex.quote(str(wrapper_bin))}; "
        "GG_CONTAINER_RUNTIME=docker; "
        "gg_detect_container_runtime_binary"
    )

    completed = _run_bash(command, cwd=REPO_ROOT)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(wrapper_bin / "singularity")


def test_set_singularity_command_auto_falls_back_to_pulled_public_image_when_sif_is_missing(tmp_path: Path):
    wrapper_bin = tmp_path / "gg_wrapper_bin"
    bin_dir = _prepare_stub_docker_with_known_images(
        tmp_path, ["ghcr.io/kfuku52/genegalleon:latest"]
    )
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"gg_support_dir={shlex.quote(str(SUPPORT_DIR))}; "
        f"gg_workflow_dir={shlex.quote(str(WORKFLOW_DIR))}; "
        f"gg_container_image_path={shlex.quote(str(tmp_path / 'missing.sif'))}; "
        f"GG_WRAPPER_BIN={shlex.quote(str(wrapper_bin))}; "
        "set_singularity_command; "
        'runtime_bin="$(gg_container_shell_command_runtime_bin)"; '
        'printf "runtime=%s\\nimage=%s\\n" "${runtime_bin}" "${GG_CONTAINER_DOCKER_IMAGE:-}"'
    )

    completed = _run_bash(
        command,
        cwd=REPO_ROOT,
        env={
            "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}",
            "GG_CONTAINER_DOCKER_IMAGE": "",
            "GG_WRAPPER_IMAGE": "",
        },
    )

    assert completed.returncode == 0, completed.stderr
    assert f"runtime={wrapper_bin / 'singularity'}" in completed.stdout
    assert "image=ghcr.io/kfuku52/genegalleon:latest" in completed.stdout


def test_docker_singularity_shim_translates_bind_mounts_and_envs(tmp_path: Path):
    bin_dir = _prepare_stub_docker(tmp_path)
    env = {
        "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}",
        "GG_CONTAINER_DOCKER_IMAGE": "local/genegalleon:dev",
        "GG_CONTAINER_BIND_MOUNTS": "/host/workspace:/workspace,/host/workflow:/script",
        "SINGULARITYENV_FOO": "foo",
        "APPTAINERENV_BAR": "bar",
    }

    completed = subprocess.run(
        ["bash", str(SHIM_PATH), "exec", "--contain", "/tmp/genegalleon.sif", "bash", "-s", "--"],
        input="echo shim-test\n",
        capture_output=True,
        text=True,
        env={**os.environ, **env},
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    args = (tmp_path / "call_1.args").read_text(encoding="utf-8").splitlines()
    assert args[0] == "run"
    assert "--rm" in args
    assert "-i" in args
    assert args[args.index("--label") + 1] == "io.github.kfuku52.project=genegalleon"
    assert "--contain" not in args
    assert args[args.index("--user") + 1] == f"{os.getuid()}:{os.getgid()}"
    assert "/host/workspace:/workspace" in args
    assert "/host/workflow:/script" in args
    assert "BAR=bar" in args
    assert "FOO=foo" in args
    assert "HOME=/workspace" in args
    assert args[args.index("-w") + 1] == "/workspace"
    assert args[-4:] == ["local/genegalleon:dev", "bash", "-s", "--"]
    assert (tmp_path / "call_1.stdin").read_text(encoding="utf-8") == "echo shim-test\n"


@pytest.mark.parametrize("devices", ["0", "2", "3", "0,1", "2,4", "GPU-498ed493-27bc-b686-e464-e1c4087f0626"])
def test_docker_nv_exposes_only_explicit_scheduler_devices(tmp_path, devices):
    bin_dir = _prepare_stub_docker(tmp_path)
    env = {**os.environ, "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}",
           "GG_CONTAINER_DOCKER_IMAGE": "fixture", "CUDA_VISIBLE_DEVICES": devices,
           "GG_CONTAINER_DOCKER_GPU_DEVICES": "", "APPTAINERENV_CUDA_VISIBLE_DEVICES": devices,
           "SINGULARITYENV_CUDA_VISIBLE_DEVICES": devices}
    result = subprocess.run(["bash", str(SHIM_PATH), "exec", "--nv", "fixture.sif", "true"],
                            capture_output=True, text=True, env=env, input="", timeout=10)
    assert result.returncode == 0, result.stderr
    args = (tmp_path / "call_1.args").read_text().splitlines()
    expected = f'"device={devices}"' if ',' in devices else f'device={devices}'
    assert args[args.index('--gpus') + 1] == expected
    assert '--nv' not in args and 'all' not in args
    assert f'CUDA_VISIBLE_DEVICES={devices}' in args
    visibility = [value for value in args if value.startswith('CUDA_VISIBLE_DEVICES=')]
    assert visibility[-1] == 'CUDA_VISIBLE_DEVICES=' + ','.join(str(i) for i in range(len(devices.split(','))))


@pytest.mark.parametrize("devices", ["", "all", "none", "void", "-1", "0,,1", "0,", ",0", "0,0", "0;echo unsafe", "GPU-invalid"])
def test_docker_nv_rejects_missing_or_unbounded_gpu_selection(tmp_path, devices):
    bin_dir = _prepare_stub_docker(tmp_path)
    env = {**os.environ, "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}",
           "GG_CONTAINER_DOCKER_IMAGE": "fixture", "CUDA_VISIBLE_DEVICES": devices,
           "GG_CONTAINER_DOCKER_GPU_DEVICES": "", "APPTAINERENV_CUDA_VISIBLE_DEVICES": "",
           "SINGULARITYENV_CUDA_VISIBLE_DEVICES": ""}
    result = subprocess.run(["bash", str(SHIM_PATH), "exec", "--nv", "fixture.sif", "true"],
                            capture_output=True, text=True, env=env, input="", timeout=10)
    assert result.returncode != 0
    assert not (tmp_path / 'docker.counter').exists()


@pytest.mark.parametrize("conflict", ["GG_CONTAINER_DOCKER_GPU_DEVICES", "APPTAINERENV_CUDA_VISIBLE_DEVICES", "SINGULARITYENV_CUDA_VISIBLE_DEVICES"])
def test_docker_nv_refuses_allocator_and_forwarded_visibility_conflicts(tmp_path, conflict):
    bin_dir = _prepare_stub_docker(tmp_path)
    env = {**os.environ, "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}",
           "GG_CONTAINER_DOCKER_IMAGE": "fixture", "CUDA_VISIBLE_DEVICES": '0',
           "GG_CONTAINER_DOCKER_GPU_DEVICES": "", "APPTAINERENV_CUDA_VISIBLE_DEVICES": "",
           "SINGULARITYENV_CUDA_VISIBLE_DEVICES": "", conflict: '1'}
    result = subprocess.run(["bash", str(SHIM_PATH), "exec", "--nv", "fixture.sif", "true"],
                            capture_output=True, text=True, env=env, input="", timeout=10)
    assert result.returncode != 0
    assert 'must match' in result.stderr
    assert not (tmp_path / 'docker.counter').exists()


def test_docker_gpu_override_is_explicit_when_no_scheduler_visibility_exists(tmp_path):
    bin_dir = _prepare_stub_docker(tmp_path)
    env = {**os.environ, "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}",
           "GG_CONTAINER_DOCKER_IMAGE": "fixture", "CUDA_VISIBLE_DEVICES": '',
           "GG_CONTAINER_DOCKER_GPU_DEVICES": '2', "APPTAINERENV_CUDA_VISIBLE_DEVICES": "",
           "SINGULARITYENV_CUDA_VISIBLE_DEVICES": ""}
    result = subprocess.run(["bash", str(SHIM_PATH), "exec", "--nv", "fixture.sif", "true"],
                            capture_output=True, text=True, env=env, input="", timeout=10)
    assert result.returncode == 0, result.stderr
    args = (tmp_path / "call_1.args").read_text().splitlines()
    assert args[args.index('--gpus') + 1] == 'device=2'


def test_docker_without_nv_retains_cpu_behavior_even_when_gpu_env_is_set(tmp_path):
    bin_dir = _prepare_stub_docker(tmp_path)
    env = {**os.environ, "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}",
           "GG_CONTAINER_DOCKER_IMAGE": "fixture", "CUDA_VISIBLE_DEVICES": '0',
           "GG_CONTAINER_DOCKER_GPU_DEVICES": 'all'}
    result = subprocess.run(["bash", str(SHIM_PATH), "exec", "fixture.sif", "true"],
                            capture_output=True, text=True, env=env, input="", timeout=10)
    assert result.returncode == 0, result.stderr
    assert '--gpus' not in (tmp_path / 'call_1.args').read_text().splitlines()


@pytest.mark.parametrize("flag", ["--nv=true", "--nvccli", "--rocm"])
def test_docker_does_not_silently_drop_other_gpu_requests(tmp_path, flag):
    bin_dir = _prepare_stub_docker(tmp_path)
    env = {**os.environ, "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}",
           "GG_CONTAINER_DOCKER_IMAGE": "fixture", "CUDA_VISIBLE_DEVICES": '0'}
    result = subprocess.run(["bash", str(SHIM_PATH), "exec", flag, "fixture.sif", "true"],
                            capture_output=True, text=True, env=env, input="", timeout=10)
    assert result.returncode != 0
    assert 'Unsupported Docker GPU option' in result.stderr
    assert not (tmp_path / 'docker.counter').exists()


def test_print_version_summary_formats_container_version_with_tag_when_distinct():
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"gg_workflow_dir={shlex.quote(str(WORKFLOW_DIR))}; "
        "gg_container_image_path=/tmp/genegalleon.sif; "
        "inspect_text=$(printf '%s\\n' "
        "'org.opencontainers.image.version: 0.4.7' "
        "'io.genegalleon.local_image_tag: 20260304-abcd123'); "
        'gg_print_version_summary "" "" "$inspect_text"'
    )

    completed = _run_bash(command, cwd=REPO_ROOT)

    assert completed.returncode == 0, completed.stderr
    assert f"genegalleon version: {REPO_VERSION}" in completed.stdout
    assert "container version: 0.4.7 (20260304-abcd123)" in completed.stdout


def test_print_version_summary_warns_when_repo_and_container_versions_mismatch():
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"gg_workflow_dir={shlex.quote(str(WORKFLOW_DIR))}; "
        "gg_container_image_path=/tmp/genegalleon.sif; "
        "inspect_text=$(printf '%s\\n' "
        "'org.opencontainers.image.version: 0.4.6' "
        "'io.genegalleon.local_image_tag: 20260304-abcd123'); "
        'gg_print_version_summary "" "" "$inspect_text"'
    )

    completed = _run_bash(command, cwd=REPO_ROOT)

    assert completed.returncode == 0, completed.stderr
    assert f"genegalleon version: {REPO_VERSION}" in completed.stdout
    assert "container version: 0.4.6 (20260304-abcd123)" in completed.stdout
    assert (
        f"WARNING: genegalleon version ({REPO_VERSION}) does not match container version (0.4.6)."
        in completed.stdout
    )


def test_progress_summary_entrypoint_dispatches_to_docker_shim(tmp_path: Path):
    workspace_dir = tmp_path / "workspace"
    wrapper_bin = tmp_path / "gg_wrapper_bin"
    workspace_dir.mkdir()
    bin_dir = _prepare_stub_docker(tmp_path)
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}{os.pathsep}{env['PATH']}"
    env["GG_CONTAINER_RUNTIME"] = "docker"
    env["GG_CONTAINER_DOCKER_IMAGE"] = "local/genegalleon:dev"
    env["GG_WRAPPER_BIN"] = str(wrapper_bin)
    env["GG_SITE_PROFILE"] = "default"
    env["gg_workspace_dir"] = str(workspace_dir)

    completed = subprocess.run(
        ["bash", str(PROGRESS_ENTRYPOINT)],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert (tmp_path / "docker.counter").read_text(encoding="utf-8") == "2"

    first_args = (tmp_path / "call_1.args").read_text(encoding="utf-8").splitlines()
    second_args = (tmp_path / "call_2.args").read_text(encoding="utf-8").splitlines()
    assert "local/genegalleon:dev" in first_args
    assert FIXTURE_IMAGE_ID in second_args
    assert f"{workspace_dir}:/workspace" in first_args
    assert f"{WORKFLOW_DIR}:/script" in first_args
    assert (tmp_path / "call_1.stdin").read_text(encoding="utf-8")
    assert (tmp_path / "call_2.stdin").read_text(encoding="utf-8")
    assert f"genegalleon version: {REPO_VERSION}" in completed.stdout
    assert "container version: dev" in completed.stdout
    assert "WARNING: genegalleon version" not in completed.stdout

    versions_dir = workspace_dir / "output" / "versions"
    version_logs = sorted(versions_dir.glob("container.*.versions.log"))
    assert len(version_logs) == 1
    run_logs = list(versions_dir.glob("run.*.log"))
    assert len(run_logs) == 1
    assert "context=" in run_logs[0].read_text(encoding="utf-8")
    version_log_text = version_logs[0].read_text(encoding="utf-8")
    assert f"genegalleon version: {REPO_VERSION}" in version_log_text
    assert "container version: dev" in version_log_text
    assert "WARNING: genegalleon version" not in version_log_text


def test_progress_summary_entrypoint_auto_detects_pulled_public_image_without_runtime_envs(tmp_path: Path):
    workspace_dir = tmp_path / "workspace"
    wrapper_bin = tmp_path / "gg_wrapper_bin"
    workspace_dir.mkdir()
    bin_dir = _prepare_stub_docker_with_known_images(
        tmp_path, ["ghcr.io/kfuku52/genegalleon:latest"]
    )
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}{os.pathsep}{env['PATH']}"
    env["GG_WRAPPER_BIN"] = str(wrapper_bin)
    env["GG_SITE_PROFILE"] = "default"
    env["gg_workspace_dir"] = str(workspace_dir)
    env["gg_container_image_path"] = str(tmp_path / "missing.sif")
    env.pop("GG_CONTAINER_RUNTIME", None)
    env.pop("GG_CONTAINER_DOCKER_IMAGE", None)

    completed = subprocess.run(
        ["bash", str(PROGRESS_ENTRYPOINT)],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "container image ref: ghcr.io/kfuku52/genegalleon:latest" in completed.stdout
    run_call_files = sorted(tmp_path.glob("call_*.args"))
    assert run_call_files
    run_calls = [
        path.read_text(encoding="utf-8").splitlines()
        for path in run_call_files
        if path.read_text(encoding="utf-8").splitlines()[:1] == ["run"]
    ]
    assert len(run_calls) == 2
    assert any(FIXTURE_IMAGE_ID in call for call in run_calls)
    assert any("ghcr.io/kfuku52/genegalleon:latest" in call for call in run_calls)
