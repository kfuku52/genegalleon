import shlex
import subprocess
from pathlib import Path

from shell_static_helpers import WORKFLOW_DIR

SITE_RUNTIME_PATH = WORKFLOW_DIR / "support" / "gg_site_runtime.sh"


def run_bash(command: str, cwd: Path):
    return subprocess.run(
        ["bash", "-c", command],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )


def make_runtime(path: Path):
    path.parent.mkdir(parents=True)
    path.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
    path.chmod(0o755)


def test_nig_runtime_discovery_prepends_versioned_apptainer_once(tmp_path):
    package_root = tmp_path / "pkg"
    runtime_bin = package_root / "apptainer" / "1.4.5" / "bin" / "apptainer"
    make_runtime(runtime_bin)
    base_path = "/usr/bin:/bin"

    command = (
        f"export PATH={shlex.quote(base_path)}; "
        f"source {shlex.quote(str(SITE_RUNTIME_PATH))}; "
        f"gg_nig_prepend_container_runtime_path "
        f"{shlex.quote(str(package_root))} {shlex.quote(str(tmp_path / 'legacy'))}; "
        f"gg_nig_prepend_container_runtime_path "
        f"{shlex.quote(str(package_root))} {shlex.quote(str(tmp_path / 'legacy'))}; "
        'printf "runtime=%s\\npath=%s\\n" "$(command -v apptainer)" "$PATH"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        f"runtime={runtime_bin}",
        f"path={runtime_bin.parent}:{base_path}",
    ]


def test_nig_runtime_discovery_preserves_existing_runtime_path(tmp_path):
    existing_bin = tmp_path / "existing" / "bin"
    existing_runtime = existing_bin / "apptainer"
    make_runtime(existing_runtime)
    package_root = tmp_path / "pkg"
    make_runtime(package_root / "apptainer" / "1.4.5" / "bin" / "apptainer")
    base_path = f"{existing_bin}:/usr/bin:/bin"

    command = (
        f"export PATH={shlex.quote(base_path)}; "
        f"source {shlex.quote(str(SITE_RUNTIME_PATH))}; "
        f"gg_nig_prepend_container_runtime_path "
        f"{shlex.quote(str(package_root))} {shlex.quote(str(tmp_path / 'legacy'))}; "
        'printf "runtime=%s\\npath=%s\\n" "$(command -v apptainer)" "$PATH"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        f"runtime={existing_runtime}",
        f"path={base_path}",
    ]


def test_nig_runtime_discovery_falls_back_to_versioned_singularity(tmp_path):
    package_root = tmp_path / "pkg"
    runtime_bin = package_root / "singularity" / "4.3.0" / "bin" / "singularity"
    make_runtime(runtime_bin)
    base_path = "/usr/bin:/bin"

    command = (
        f"export PATH={shlex.quote(base_path)}; "
        f"source {shlex.quote(str(SITE_RUNTIME_PATH))}; "
        f"gg_nig_prepend_container_runtime_path "
        f"{shlex.quote(str(package_root))} {shlex.quote(str(tmp_path / 'legacy'))}; "
        'printf "runtime=%s\\npath=%s\\n" "$(command -v singularity)" "$PATH"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        f"runtime={runtime_bin}",
        f"path={runtime_bin.parent}:{base_path}",
    ]


def test_nig_scheduler_prelude_discovers_runtime_without_pbs_workdir(tmp_path):
    marker = tmp_path / "runtime-discovery-called"

    command = (
        "export GG_SITE_PROFILE=nig; "
        "unset PBS_O_WORKDIR; "
        f"source {shlex.quote(str(SITE_RUNTIME_PATH))}; "
        "gg_nig_prepend_container_runtime_path() { "
        f": > {shlex.quote(str(marker))}; "
        "}; "
        "gg_site_scheduler_prelude"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert marker.is_file()
