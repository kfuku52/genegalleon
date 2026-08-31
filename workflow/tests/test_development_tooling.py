import ast
import json
import os
import re
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest
from shell_static_helpers import REPO_ROOT

RUN_IN_RUNTIME = REPO_ROOT / "workflow" / "tests" / "run_in_runtime.sh"
SCHEMA_TOOL = REPO_ROOT / "workflow" / "support" / "entrypoint_config_schema.py"
VERSION_TOOL = REPO_ROOT / "workflow" / "support" / "bump_version.py"
BUILD_HASH_TOOL = REPO_ROOT / "container" / "scripts" / "compute_build_input_hash.sh"
RUNTIME_FRESHNESS_TOOL = REPO_ROOT / "container" / "scripts" / "check_runtime_freshness.sh"
FETCH_GIT_TOOL = REPO_ROOT / "container" / "scripts" / "fetch_git_repo.sh"
DOCKER_TO_SIF_TOOL = REPO_ROOT / "container" / "apptainer_from_docker.sh"
DOWNLOAD_URL_TOOL = REPO_ROOT / "container" / "scripts" / "download_url.sh"
SOURCE_SHA_VARS = (
    "KFU52_AMALGKIT_REPO_SHA",
    "KFU52_CDSKIT_REPO_SHA",
    "KFU52_CSUBST_REPO_SHA",
    "KFU52_NWKIT_REPO_SHA",
    "BUSCO_REPO_SHA",
    "PAML_REPO_SHA",
    "KFL1OU_REPO_SHA",
    "KFFRACTBIAS_REPO_SHA",
    "KFTOOLS_REPO_SHA",
    "RKFTOOLS_REPO_SHA",
    "RADTE_REPO_SHA",
)


def _run(*args: str, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(args),
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )


def test_generated_config_schema_matches_registry_and_effective_defaults():
    checked = _run("python", str(SCHEMA_TOOL), "--check")
    assert checked.returncode == 0, checked.stderr

    rendered = _run("python", str(SCHEMA_TOOL), "--format", "json")
    assert rendered.returncode == 0, rendered.stderr
    schema = json.loads(rendered.stdout)
    common = {item["name"]: item for item in schema["common_parameters"]}
    transcriptome = {
        item["name"]: item for item in schema["entrypoints"]["gg_transcriptome_generation_entrypoint.sh"]["parameters"]
    }

    assert common["GG_COMMON_BUSCO_LINEAGE"]["default"] == "eukaryota_odb12"
    assert transcriptome["assembly_method"]["default"] == "auto"
    assert transcriptome["assembly_method"]["environment"] == "GG_TRANSCRIPTOME_ASSEMBLY_METHOD"


def test_version_helper_previews_and_updates_semver_atomically(tmp_path: Path):
    version_file = tmp_path / "VERSION"
    version_file.write_text("1.2.3\n", encoding="utf-8")

    preview = _run(
        "python",
        str(VERSION_TOOL),
        "minor",
        "--dry-run",
        "--version-file",
        str(version_file),
    )
    assert preview.returncode == 0, preview.stderr
    assert preview.stdout.strip() == "1.3.0"
    assert version_file.read_text(encoding="utf-8") == "1.2.3\n"

    updated = _run(
        "python",
        str(VERSION_TOOL),
        "patch",
        "--version-file",
        str(version_file),
    )
    assert updated.returncode == 0, updated.stderr
    assert updated.stdout.strip() == "1.2.4"
    assert version_file.read_text(encoding="utf-8") == "1.2.4\n"


def test_runtime_runner_dispatches_to_requested_docker_image(tmp_path: Path):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "python3").symlink_to(sys.executable)
    docker_log = tmp_path / "docker.log"
    docker = bin_dir / "docker"
    docker.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        f"printf '%s\\n' \"$*\" >> {docker_log!s}\n"
        'case "${1:-}" in\n'
        "  info) exit 0 ;;\n"
        '  image) [[ "${2:-}" == inspect ]] && exit 0 ;;\n'
        "  run) exit 0 ;;\n"
        "esac\n"
        "exit 1\n",
        encoding="utf-8",
    )
    docker.chmod(docker.stat().st_mode | stat.S_IXUSR)
    env = os.environ.copy()
    env.update(
        {
            "PATH": f"{bin_dir}:/usr/bin:/bin",
            "GG_TEST_RUNTIME": "docker",
            "GG_CONTAINER_DOCKER_IMAGE": "local/genegalleon:test",
            "GG_RUNTIME_FRESHNESS": "off",
        }
    )

    completed = _run("bash", str(RUN_IN_RUNTIME), "python", "--version", env=env)

    assert completed.returncode == 0, completed.stderr
    calls = docker_log.read_text(encoding="utf-8")
    assert "image inspect local/genegalleon:test" in calls
    assert f"--user {os.getuid()}:{os.getgid()}" in calls
    assert f"--volume {REPO_ROOT}:{REPO_ROOT}" in calls
    assert f"--workdir {REPO_ROOT} local/genegalleon:test python --version" in calls


@pytest.mark.parametrize("target", [None, "runtime"])
def test_dev_build_selects_development_unless_overridden(tmp_path, target):
    shutil.copy2(REPO_ROOT / "dev", tmp_path / "dev")
    (tmp_path / "gg_container_build_entrypoint.sh").write_text(
        'printf "%s %s %s %s\\n" "$BUILD_TARGET" "$IMAGE_SOURCE" "$IMAGE:$TAG" "$BUILD_SIF"\n'
    )
    env = os.environ.copy()
    env.pop("BUILD_TARGET", None)
    if target:
        env["BUILD_TARGET"] = target
    completed = _run("bash", str(tmp_path / "dev"), "build", env=env)
    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == f"{target or 'development'} local local/genegalleon:dev 0"


@pytest.mark.parametrize("runtime", ["auto", "sif"])
@pytest.mark.parametrize("entrypoint", ["sif", "runtime", "dev"])
def test_validation_entrypoints_discover_versioned_hpc_runtime(tmp_path, runtime, entrypoint):
    package_root = tmp_path / "packages"
    engine = package_root / "apptainer" / "1.4.5" / "bin" / "apptainer"
    engine.parent.mkdir(parents=True)
    engine.write_text('#!/bin/bash\nprintf "%s\\n" "$@"\n')
    engine.chmod(0o755)
    sif = tmp_path / "genegalleon.sif"
    sif.touch()
    env = os.environ.copy()
    env.update({"PATH": "/usr/bin:/bin", "GENEGALLEON_SIF": str(sif), "GG_TEST_RUNTIME": runtime,
                "GG_CONTAINER_RUNTIME_PACKAGE_ROOT": str(package_root), "GG_CONTAINER_RUNTIME_LEGACY_DIR": str(tmp_path / "none"),
                "GG_RUNTIME_FRESHNESS": "off"})
    if entrypoint == "dev":
        command = ["bash", str(REPO_ROOT / "dev"), "check", "fast"]
    else:
        script = RUN_IN_RUNTIME if entrypoint == "runtime" else RUN_IN_RUNTIME.with_name("run_in_sif.sh")
        command = ["bash", str(script), "python", "--version"]
    completed = _run(*command, env=env)
    assert completed.returncode == 0, completed.stderr
    assert "--cleanenv\n" in completed.stdout
    assert str(sif) in completed.stdout
    assert "PYTHONDONTWRITEBYTECODE=1" in completed.stdout


def test_build_hash_ignores_import_and_editor_artifacts_but_tracks_copied_inputs(tmp_path):
    root = tmp_path / "repo"
    shutil.copytree(REPO_ROOT / "container", root / "container", ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copytree(REPO_ROOT / "workflow/support/treevis", root / "workflow/support/treevis")
    shutil.copy2(REPO_ROOT / ".dockerignore", root / ".dockerignore")
    tool = root / "container/scripts/compute_build_input_hash.sh"

    def compute():
        completed = _run("bash", str(tool), "--runtime", "2026-08-30")
        assert completed.returncode == 0, completed.stderr
        return completed.stdout

    before = compute()
    imported = _run("python3", "-c", "import importlib.util, sys; sys.dont_write_bytecode=False; "
                    "s=importlib.util.spec_from_file_location('extract',sys.argv[1]); "
                    "s.loader.exec_module(importlib.util.module_from_spec(s))",
                    str(root / "container/scripts/extract_notung_jar.py"))
    assert imported.returncode == 0, imported.stderr
    assert list((root / "container/scripts/__pycache__").glob("*.pyc"))
    (root / "container/env/.DS_Store").write_text("finder metadata")
    (root / "container/env/base.required.txt~").write_text("editor backup")
    editor_cache = root / "container/env/editor.swp/content.txt"
    editor_cache.parent.mkdir()
    editor_cache.write_text("ignored directory contents")
    (root / "container/scripts/not-copied-into-image.sh").write_text("not a Docker input")
    assert compute() == before
    (root / "container/env/new-untracked-input.txt").write_text("new runtime input")
    added = compute()
    assert added != before
    with (root / "container/env/base.required.txt").open("a") as handle:
        handle.write("\n# actual copied source edit\n")
    assert compute() != added
    staging = tmp_path / "native-context"
    staged = _run("python3", str(root / "container/scripts/list_build_inputs.py"),
                  "--stage-native-context", str(staging))
    assert staged.returncode == 0, staged.stderr
    assert (staging / "env/new-untracked-input.txt").is_file()
    assert (staging / "treevis/DESCRIPTION").is_file()
    assert (staging / "source_branches.env").read_bytes() == (root / "container/source_branches.env").read_bytes()
    assert (staging / "pip-compatibility.requirements.txt").read_bytes() == (
        root / "container/pip-compatibility.requirements.txt"
    ).read_bytes()
    assert not list(staging.rglob("*.pyc"))
    assert not list(staging.rglob(".DS_Store"))
    assert not (staging / "scripts/not-copied-into-image.sh").exists()


def test_docker_context_excludes_the_same_generated_files_as_build_hash():
    patterns = set((REPO_ROOT / ".dockerignore").read_text().splitlines())
    assert {"**/__pycache__", "**/.pytest_cache", "**/.ruff_cache", "**/*.pyc", "**/*.pyo",
            "**/.DS_Store", "**/*.swp", "**/*.swo", "**/*~"} <= patterns


def test_build_input_hash_is_deterministic_and_tracks_version():
    env = os.environ.copy()
    env["GG_BUILD_VERSION"] = "1.2.3"
    first = _run("bash", str(BUILD_HASH_TOOL), "2026-08-26", env=env)
    second = _run("bash", str(BUILD_HASH_TOOL), "2026-08-26", env=env)
    assert first.returncode == second.returncode == 0
    assert re.fullmatch(r"[0-9a-f]{64}", first.stdout.strip())
    assert first.stdout == second.stdout

    env["GG_BUILD_VERSION"] = "1.2.4"
    changed = _run("bash", str(BUILD_HASH_TOOL), "2026-08-26", env=env)
    assert changed.returncode == 0
    assert changed.stdout != first.stdout

    env["GG_BUILD_VERSION"] = "1.2.3"
    next_security_epoch = _run("bash", str(BUILD_HASH_TOOL), "2026-08-27", env=env)
    assert next_security_epoch.returncode == 0
    assert next_security_epoch.stdout != first.stdout


def test_runtime_input_hash_tracks_runtime_content_but_not_release_metadata():
    env = os.environ.copy()
    env.update(
        {
            "GG_BUILD_VERSION": "1.2.3",
            "GG_BUILD_VCS_REF": "a" * 40,
            "KFU52_NWKIT_REPO_SHA": "b" * 40,
        }
    )
    first = _run("bash", str(BUILD_HASH_TOOL), "--runtime", "2026-08-27", env=env)
    assert first.returncode == 0, first.stderr
    assert re.fullmatch(r"[0-9a-f]{64}", first.stdout.strip())

    env["GG_BUILD_VERSION"] = "9.9.9"
    env["GG_BUILD_VCS_REF"] = "c" * 40
    metadata_changed = _run("bash", str(BUILD_HASH_TOOL), "--runtime", "2026-08-27", env=env)
    assert metadata_changed.returncode == 0, metadata_changed.stderr
    assert metadata_changed.stdout == first.stdout

    env["KFU52_NWKIT_REPO_SHA"] = "d" * 40
    upstream_changed = _run("bash", str(BUILD_HASH_TOOL), "--runtime", "2026-08-27", env=env)
    assert upstream_changed.returncode == 0, upstream_changed.stderr
    assert upstream_changed.stdout != first.stdout

    env["KFU52_NWKIT_REPO_SHA"] = "b" * 40
    next_security_epoch = _run("bash", str(BUILD_HASH_TOOL), "--runtime", "2026-08-28", env=env)
    assert next_security_epoch.returncode == 0, next_security_epoch.stderr
    assert next_security_epoch.stdout != first.stdout


def test_runtime_input_hash_is_platform_agnostic():
    env = os.environ.copy()
    env["GG_BUILD_PLATFORMS"] = "linux/amd64"
    amd64 = _run("bash", str(BUILD_HASH_TOOL), "--runtime", "2026-08-27", env=env)
    assert amd64.returncode == 0, amd64.stderr

    env["GG_BUILD_PLATFORMS"] = "linux/amd64,linux/arm64"
    multi_arch = _run("bash", str(BUILD_HASH_TOOL), "--runtime", "2026-08-27", env=env)

    assert multi_arch.returncode == 0, multi_arch.stderr
    assert multi_arch.stdout == amd64.stdout


@pytest.mark.parametrize("hash_args", [[], ["--runtime"]])
def test_build_hash_distinguishes_targets_but_ignores_build_parallelism(hash_args):
    env = os.environ.copy()
    env.pop("GG_BUILD_TARGET", None)
    env["BUILD_TARGET"] = "runtime"
    runtime = _run("bash", str(BUILD_HASH_TOOL), *hash_args, "2026-08-27", env=env)
    env["BUILD_TARGET"] = "development"
    development = _run("bash", str(BUILD_HASH_TOOL), *hash_args, "2026-08-27", env=env)
    assert runtime.returncode == development.returncode == 0
    assert development.stdout != runtime.stdout
    env["GG_BUILD_JOBS"] = "7"
    more_jobs = _run("bash", str(BUILD_HASH_TOOL), *hash_args, "2026-08-27", env=env)
    assert more_jobs.returncode == 0, more_jobs.stderr
    assert more_jobs.stdout == development.stdout


@pytest.mark.parametrize("target", ["runtime", "development"])
def test_docker_runtime_freshness_uses_exact_runtime_hash_and_fails_closed(tmp_path: Path, target: str):
    env = os.environ.copy()
    for index, variable in enumerate(SOURCE_SHA_VARS, start=1):
        env[variable] = f"{index:040x}"
    env.update(
        {
            "GG_BUILD_PLATFORMS": "linux/amd64",
            "GG_BUILD_VCS_REF": "runtime-test",
            "GG_BUILD_VERSION": "runtime-test",
            "GG_BUILD_TARGET": target,
        }
    )
    expected = _run("bash", str(BUILD_HASH_TOOL), "--runtime", "2026-08-27", env=env)
    assert expected.returncode == 0, expected.stderr
    busco_sha = env.pop("BUSCO_REPO_SHA")
    paml_sha = env.pop("PAML_REPO_SHA")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "python3").symlink_to(sys.executable)
    docker = bin_dir / "docker"
    docker.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        'case "${1:-}" in\n'
        f"  image) printf '%s\\n' 'amd64|{expected.stdout.strip()}|2026-08-27|{target}' ;;\n"
        "  run)\n"
        "    printf 'source\\trevision\\n'\n"
        f"    printf 'BUSCO\\t%s\\n' '{busco_sha}'\n"
        f"    printf 'paml\\t%s\\n' '{paml_sha}'\n"
        "    ;;\n"
        "  *) exit 1 ;;\n"
        "esac\n",
        encoding="utf-8",
    )
    docker.chmod(docker.stat().st_mode | stat.S_IXUSR)
    env.update(
        {
            "PATH": f"{bin_dir}:/usr/bin:/bin",
            "GG_RUNTIME_FRESHNESS": "always",
            "GG_RUNTIME_FRESHNESS_CACHE_DIR": str(tmp_path / "cache"),
            # Freshness must use the image's target even if the caller requests
            # a different target for their next build.
            "GG_BUILD_TARGET": "development" if target == "runtime" else "runtime",
        }
    )

    current = _run("bash", str(RUNTIME_FRESHNESS_TOOL), "--docker", "local/test:dev", env=env)
    assert current.returncode == 0, current.stderr
    assert "matches the daily owned upstream snapshot" in current.stdout

    env["KFU52_NWKIT_REPO_SHA"] = "f" * 40
    stale = _run("bash", str(RUNTIME_FRESHNESS_TOOL), "--docker", "local/test:dev", env=env)
    assert stale.returncode == 1
    assert "runtime is stale" in stale.stderr
    assert expected.stdout.strip() in stale.stderr


def test_docker_runtime_freshness_daily_cache_honors_one_off_revision_overrides(
    tmp_path: Path,
):
    base_env = os.environ.copy()
    for index, variable in enumerate(SOURCE_SHA_VARS, start=1):
        base_env[variable] = f"{index:040x}"
    base_env.update(
        {
            "GG_BUILD_PLATFORMS": "linux/amd64",
            "GG_BUILD_VCS_REF": "runtime-test",
            "GG_BUILD_VERSION": "runtime-test",
        }
    )
    base_hash = _run(
        "bash", str(BUILD_HASH_TOOL), "--runtime", "2026-08-27", env=base_env
    )
    assert base_hash.returncode == 0, base_hash.stderr

    override_env = base_env.copy()
    override_env["KFU52_NWKIT_REPO_SHA"] = "f" * 40
    override_hash = _run(
        "bash", str(BUILD_HASH_TOOL), "--runtime", "2026-08-27", env=override_env
    )
    assert override_hash.returncode == 0, override_hash.stderr
    assert override_hash.stdout != base_hash.stdout

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "python3").symlink_to(sys.executable)
    docker = bin_dir / "docker"
    docker.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        '[[ "${1:-}" == image && "${2:-}" == inspect ]]\n'
        "printf 'amd64|%s|2026-08-27|runtime\\n' \"${MOCK_RUNTIME_HASH}\"\n",
        encoding="utf-8",
    )
    docker.chmod(docker.stat().st_mode | stat.S_IXUSR)
    shared = {
        "PATH": f"{bin_dir}:/usr/bin:/bin",
        "GG_RUNTIME_FRESHNESS": "daily",
        "GG_RUNTIME_FRESHNESS_SCOPE": "all",
        "GG_RUNTIME_FRESHNESS_CACHE_DIR": str(tmp_path / "cache"),
    }

    first_env = base_env.copy()
    first_env.update(shared)
    first_env["MOCK_RUNTIME_HASH"] = base_hash.stdout.strip()
    first = _run(
        "bash", str(RUNTIME_FRESHNESS_TOOL), "--docker", "local/test:dev", env=first_env
    )
    assert first.returncode == 0, first.stderr

    second_env = override_env.copy()
    second_env.update(shared)
    second_env["MOCK_RUNTIME_HASH"] = override_hash.stdout.strip()
    second = _run(
        "bash", str(RUNTIME_FRESHNESS_TOOL), "--docker", "local/test:dev", env=second_env
    )
    assert second.returncode == 0, second.stderr


def test_docker_runtime_freshness_daily_cache_honors_source_ref_overrides(
    tmp_path: Path,
):
    unresolved_env = os.environ.copy()
    for index, variable in enumerate(SOURCE_SHA_VARS, start=1):
        unresolved_env[variable] = f"{index:040x}"
    unresolved_env.pop("KFU52_NWKIT_REPO_SHA")
    unresolved_env.update(
        {
            "GG_BUILD_PLATFORMS": "linux/amd64",
            "GG_BUILD_VCS_REF": "runtime-test",
            "GG_BUILD_VERSION": "runtime-test",
            "KFU52_NWKIT_REPO_REF": "base",
        }
    )

    base_env = unresolved_env.copy()
    base_env["KFU52_NWKIT_REPO_SHA"] = "a" * 40
    base_hash = _run(
        "bash", str(BUILD_HASH_TOOL), "--runtime", "2026-08-27", env=base_env
    )
    assert base_hash.returncode == 0, base_hash.stderr

    override_unresolved_env = unresolved_env.copy()
    override_unresolved_env["KFU52_NWKIT_REPO_REF"] = "override"
    override_env = override_unresolved_env.copy()
    override_env["KFU52_NWKIT_REPO_SHA"] = "b" * 40
    override_hash = _run(
        "bash", str(BUILD_HASH_TOOL), "--runtime", "2026-08-27", env=override_env
    )
    assert override_hash.returncode == 0, override_hash.stderr
    assert override_hash.stdout != base_hash.stdout

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "python3").symlink_to(sys.executable)
    git = bin_dir / "git"
    git.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        'if [[ "${1:-}" == "-C" && "${3:-}" == "rev-parse" ]]; then\n'
        "  printf '.git\\n'\n"
        "  exit 0\n"
        "fi\n"
        '[[ "${1:-}" == "ls-remote" && "${2:-}" == "--exit-code" ]]\n'
        'case "${4:-}" in\n'
        "  refs/heads/base) revision='" + "a" * 40 + "' ;;\n"
        "  refs/heads/override) revision='" + "b" * 40 + "' ;;\n"
        "  *) exit 1 ;;\n"
        "esac\n"
        "printf '%s\\t%s\\n' \"${revision}\" \"${4}\"\n",
        encoding="utf-8",
    )
    git.chmod(git.stat().st_mode | stat.S_IXUSR)
    docker = bin_dir / "docker"
    docker.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        '[[ "${1:-}" == image && "${2:-}" == inspect ]]\n'
        "printf 'amd64|%s|2026-08-27|runtime\\n' \"${MOCK_RUNTIME_HASH}\"\n",
        encoding="utf-8",
    )
    docker.chmod(docker.stat().st_mode | stat.S_IXUSR)
    shared = {
        "PATH": f"{bin_dir}:/usr/bin:/bin",
        "GG_RUNTIME_FRESHNESS": "daily",
        "GG_RUNTIME_FRESHNESS_SCOPE": "all",
        "GG_RUNTIME_FRESHNESS_CACHE_DIR": str(tmp_path / "cache"),
    }

    first_env = unresolved_env.copy()
    first_env.update(shared)
    first_env["MOCK_RUNTIME_HASH"] = base_hash.stdout.strip()
    first = _run(
        "bash", str(RUNTIME_FRESHNESS_TOOL), "--docker", "local/test:dev", env=first_env
    )
    assert first.returncode == 0, first.stderr

    second_env = override_unresolved_env.copy()
    second_env.update(shared)
    second_env["MOCK_RUNTIME_HASH"] = override_hash.stdout.strip()
    second = _run(
        "bash", str(RUNTIME_FRESHNESS_TOOL), "--docker", "local/test:dev", env=second_env
    )
    assert second.returncode == 0, second.stderr


def test_runtime_freshness_fallback_cache_is_private_and_scoped_by_uid(tmp_path: Path):
    env = os.environ.copy()
    for index, variable in enumerate(SOURCE_SHA_VARS, start=1):
        env[variable] = f"{index:040x}"
    env.update(
        {
            "GG_BUILD_VCS_REF": "runtime-test",
            "GG_BUILD_VERSION": "runtime-test",
            "GG_RUNTIME_FRESHNESS": "always",
            "GG_RUNTIME_FRESHNESS_SCOPE": "all",
        }
    )
    env.pop("GG_RUNTIME_FRESHNESS_CACHE_DIR", None)
    expected = _run("bash", str(BUILD_HASH_TOOL), "--runtime", "2026-08-27", env=env)
    assert expected.returncode == 0, expected.stderr

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "python3").symlink_to(sys.executable)
    git = bin_dir / "git"
    git.write_text("#!/usr/bin/env bash\nexit 1\n", encoding="utf-8")
    git.chmod(git.stat().st_mode | stat.S_IXUSR)
    docker = bin_dir / "docker"
    docker.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        '[[ "${1:-}" == image && "${2:-}" == inspect ]]\n'
        f"printf 'amd64|{expected.stdout.strip()}|2026-08-27|runtime\\n'\n",
        encoding="utf-8",
    )
    docker.chmod(docker.stat().st_mode | stat.S_IXUSR)
    tmp_root = tmp_path / "shared-tmp"
    tmp_root.mkdir()
    env.update(
        {
            "PATH": f"{bin_dir}:/usr/bin:/bin",
            "TMPDIR": str(tmp_root),
        }
    )

    completed = _run(
        "bash", str(RUNTIME_FRESHNESS_TOOL), "--docker", "local/test:dev", env=env
    )

    assert completed.returncode == 0, completed.stderr
    private_parent = tmp_root / f"genegalleon_{os.getuid()}"
    assert private_parent.is_dir()
    assert stat.S_IMODE(private_parent.stat().st_mode) == 0o700
    assert any((private_parent / "gg-runtime-freshness").glob("*.env"))


def test_sif_runtime_identity_check_uses_embedded_exact_hash(tmp_path: Path):
    runtime_hash = "a" * 64
    engine = tmp_path / "singularity"
    engine.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "[[ \"${1:-}\" == inspect && \"${2:-}\" == --json ]]\n"
        "printf '%s\\n' "
        f"'{{\"data\":{{\"attributes\":{{\"labels\":{{\"io.genegalleon.runtime-input\":\"{runtime_hash}\",\"io.genegalleon.security-refresh-epoch\":\"2026-08-27\",\"io.genegalleon.build-target\":\"runtime\"}}}}}}}}'\n",
        encoding="utf-8",
    )
    engine.chmod(engine.stat().st_mode | stat.S_IXUSR)
    env = os.environ.copy()
    env["GG_RUNTIME_FRESHNESS"] = "off"

    current = _run(
        "bash",
        str(RUNTIME_FRESHNESS_TOOL),
        "--sif",
        str(tmp_path / "runtime.sif"),
        "--engine",
        str(engine),
        "--expected-hash",
        runtime_hash,
        env=env,
    )
    assert current.returncode == 0, current.stderr
    assert f"has exact runtime input {runtime_hash}" in current.stdout

    stale = _run(
        "bash",
        str(RUNTIME_FRESHNESS_TOOL),
        "--sif",
        str(tmp_path / "runtime.sif"),
        "--engine",
        str(engine),
        "--expected-hash",
        "b" * 64,
        env=env,
    )
    assert stale.returncode == 1
    assert "runtime identity mismatch" in stale.stderr



def test_python_cli_arguments_do_not_publish_empty_help_text():
    offenders: list[str] = []
    for path in sorted((REPO_ROOT / "workflow").rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
                continue
            if node.func.attr != "add_argument":
                continue
            for keyword in node.keywords:
                if keyword.arg == "help" and isinstance(keyword.value, ast.Constant) and keyword.value.value == "":
                    offenders.append(f"{path.relative_to(REPO_ROOT)}:{node.lineno}")
    assert offenders == []


def test_format_species_cli_has_no_wildcard_imports():
    path = REPO_ROOT / "workflow" / "support" / "format_species_inputs.py"
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    wildcard_imports = [
        node.lineno
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and any(alias.name == "*" for alias in node.names)
    ]
    assert wildcard_imports == []


def test_gene_family_store_cli_definition_is_a_bounded_module():
    store_path = REPO_ROOT / "workflow" / "support" / "gene_family_output_store.py"
    cli_path = REPO_ROOT / "workflow" / "support" / "gene_family_output_cli.py"

    assert cli_path.is_file()
    assert len(cli_path.read_text(encoding="utf-8").splitlines()) < 400
    assert "from gene_family_output_cli import build_parser as build_cli_parser" in store_path.read_text(
        encoding="utf-8"
    )
    completed = _run("python", str(store_path), "--help")
    assert completed.returncode == 0, completed.stderr
    assert "archive-completed" in completed.stdout
    assert "convert-storage" in completed.stdout


def test_genome_annotation_core_can_load_helpers_without_running_the_workflow():
    core = REPO_ROOT / "workflow" / "core" / "gg_genome_annotation_core.sh"
    completed = _run(
        "bash",
        "-c",
        f'GG_CORE_SOURCE_ONLY=1 source "{core}"; '
        "declare -F resolve_species_file >/dev/null; "
        "declare -F gg_genome_annotation_main >/dev/null",
    )
    assert completed.returncode == 0, completed.stderr


def test_notung_download_and_container_license_are_verified():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    definition = (REPO_ROOT / "container" / "apptainer_local_build.def.template").read_text(encoding="utf-8")
    expected_sha = "81cbff670ab4d2416c01eba503f81c454aa5a724b0982373dd17510113882ae6"

    assert f'ARG NOTUNG_ZIP_SHA256="{expected_sha}"' in dockerfile
    assert 'verify_sha256.sh "${tmp_zip}" "${NOTUNG_ZIP_SHA256}"' in dockerfile
    assert 'org.opencontainers.image.licenses="MIT"' in dockerfile
    assert "@@NOTUNG_ZIP_SHA256@@" in definition
    assert 'verify_sha256.sh "${tmp_zip}" "${notung_zip_sha256}"' in definition
    assert "org.opencontainers.image.licenses MIT" in definition


def test_container_image_rejects_invalid_input_hash_labels():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")

    assert '[[ "${BUILD_INPUT_HASH}" =~ ^[0-9a-f]{64}$ ]]' in dockerfile
    assert '[[ "${RUNTIME_INPUT_HASH}" =~ ^[0-9a-f]{64}$ ]]' in dockerfile


def test_fetch_git_repo_propagates_clone_failures_with_or_without_a_ref(tmp_path: Path):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "python3").symlink_to(sys.executable)
    git = bin_dir / "git"
    git.write_text("#!/usr/bin/env bash\nexit 1\n", encoding="utf-8")
    git.chmod(git.stat().st_mode | stat.S_IXUSR)
    env = os.environ.copy()
    env.update(
        {
            "PATH": f"{bin_dir}:/usr/bin:/bin",
            "FETCH_GIT_MAX_ATTEMPTS": "1",
            "FETCH_GIT_RETRY_SLEEP_SEC": "0",
        }
    )

    for suffix, repo_ref in (("default", ""), ("explicit", "main")):
        destination = tmp_path / suffix
        completed = _run(
            "bash",
            str(FETCH_GIT_TOOL),
            "invalid://repository",
            repo_ref,
            str(destination),
            env=env,
        )

        assert completed.returncode != 0
        assert "failed to fetch" in completed.stderr
        assert not destination.exists()


def test_docker_to_sif_conversion_preserves_existing_output_on_failure(tmp_path: Path):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "python3").symlink_to(sys.executable)
    engine = bin_dir / "apptainer"
    engine.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        '[[ "${1:-}" == build ]]\n'
        "printf 'partial\\n' > \"${2}\"\n"
        "exit 1\n",
        encoding="utf-8",
    )
    engine.chmod(engine.stat().st_mode | stat.S_IXUSR)
    output = tmp_path / "genegalleon.sif"
    output.write_text("known-good\n", encoding="utf-8")
    env = os.environ.copy()
    env.update(
        {
            "PATH": f"{bin_dir}:/usr/bin:/bin",
            "ENGINE": "apptainer",
            "OUT": str(output),
        }
    )

    completed = _run("bash", str(DOCKER_TO_SIF_TOOL), env=env)

    assert completed.returncode != 0
    assert output.read_text(encoding="utf-8") == "known-good\n"
    assert list(tmp_path.glob(".gg-sif-build.*")) == []


def test_sif_build_wrappers_publish_from_an_adjacent_staging_path():
    for relative_path in (
        "container/apptainer_from_docker.sh",
        "container/apptainer_local_build.sh",
    ):
        text = (REPO_ROOT / relative_path).read_text(encoding="utf-8")
        assert ".gg-sif-build.XXXXXX" in text
        assert '"${staged_out}"' in text
        assert 'mv -f -- "${staged_out}" "${OUT}"' in text
        assert 'build "${OUT}"' not in text


def test_download_url_preserves_existing_destination_when_transfers_fail(tmp_path: Path):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "python3").symlink_to(sys.executable)
    for command_name in ("curl", "wget"):
        command = bin_dir / command_name
        command.write_text(
            "#!/usr/bin/env bash\n"
            "set -euo pipefail\n"
            "output=''\n"
            "previous=''\n"
            "for argument in \"$@\"; do\n"
            "  if [[ \"${previous}\" == -o || \"${previous}\" == -O ]]; then output=${argument}; fi\n"
            "  previous=${argument}\n"
            "done\n"
            "[[ -n \"${output}\" ]] && printf 'partial\\n' > \"${output}\"\n"
            "exit 1\n",
            encoding="utf-8",
        )
        command.chmod(command.stat().st_mode | stat.S_IXUSR)
    destination = tmp_path / "archive.tar.gz"
    destination.write_text("known-good\n", encoding="utf-8")
    env = os.environ.copy()
    env.update(
        {
            "PATH": f"{bin_dir}:/usr/bin:/bin",
            "DOWNLOAD_URL_MAX_ATTEMPTS": "1",
            "DOWNLOAD_URL_RETRY_DELAY_SEC": "0",
        }
    )

    completed = _run(
        "bash", str(DOWNLOAD_URL_TOOL), "https://invalid.example/archive", str(destination), env=env
    )

    assert completed.returncode != 0
    assert destination.read_text(encoding="utf-8") == "known-good\n"
    assert list(tmp_path.glob("archive.tar.gz.tmp.*")) == []
