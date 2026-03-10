from pathlib import Path
import os
import shlex
import stat
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_DIR = REPO_ROOT / "workflow"
SUPPORT_DIR = WORKFLOW_DIR / "support"
GG_UTIL_PATH = SUPPORT_DIR / "gg_util.sh"
SHIM_PATH = SUPPORT_DIR / "gg_wrapper_bin" / "singularity"
PROGRESS_ENTRYPOINT = WORKFLOW_DIR / "gg_progress_summary_entrypoint.sh"
REPO_VERSION = REPO_ROOT.joinpath("VERSION").read_text(encoding="utf-8").splitlines()[0].strip()


def _write_executable(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _run_bash(command: str, cwd: Path, env: dict[str, str] | None = None):
    run_env = os.environ.copy()
    if env:
        run_env.update(env)
    return subprocess.run(
        ["bash", "-lc", command],
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
    assert "local/genegalleon:dev" in second_args
    assert f"{workspace_dir}:/workspace" in first_args
    assert f"{WORKFLOW_DIR}:/script" in first_args
    assert (tmp_path / "call_1.stdin").read_text(encoding="utf-8")
    assert (tmp_path / "call_2.stdin").read_text(encoding="utf-8")
    assert f"genegalleon version: {REPO_VERSION}" in completed.stdout
    assert "genegalleon.sif version: dev" in completed.stdout

    version_logs = sorted((workspace_dir / "output" / "versions").glob("*.log"))
    assert len(version_logs) == 1
    version_log_text = version_logs[0].read_text(encoding="utf-8")
    assert f"genegalleon version: {REPO_VERSION}" in version_log_text
    assert "genegalleon.sif version: dev" in version_log_text
