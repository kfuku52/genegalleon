import ast
import json
import os
import re
import stat
import subprocess
from pathlib import Path

from shell_static_helpers import REPO_ROOT

RUN_IN_RUNTIME = REPO_ROOT / "workflow" / "tests" / "run_in_runtime.sh"
SCHEMA_TOOL = REPO_ROOT / "workflow" / "support" / "entrypoint_config_schema.py"
VERSION_TOOL = REPO_ROOT / "workflow" / "support" / "bump_version.py"
BUILD_HASH_TOOL = REPO_ROOT / "container" / "scripts" / "compute_build_input_hash.sh"
RUNTIME_FRESHNESS_TOOL = REPO_ROOT / "container" / "scripts" / "check_runtime_freshness.sh"
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
    assert f"--volume {REPO_ROOT}:{REPO_ROOT}" in calls
    assert f"--workdir {REPO_ROOT} local/genegalleon:test python --version" in calls


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


def test_docker_runtime_freshness_uses_exact_runtime_hash_and_fails_closed(tmp_path: Path):
    env = os.environ.copy()
    for index, variable in enumerate(SOURCE_SHA_VARS, start=1):
        env[variable] = f"{index:040x}"
    env.update(
        {
            "GG_BUILD_PLATFORMS": "linux/amd64",
            "GG_BUILD_VCS_REF": "runtime-test",
            "GG_BUILD_VERSION": "runtime-test",
        }
    )
    expected = _run("bash", str(BUILD_HASH_TOOL), "--runtime", "2026-08-27", env=env)
    assert expected.returncode == 0, expected.stderr
    busco_sha = env.pop("BUSCO_REPO_SHA")
    paml_sha = env.pop("PAML_REPO_SHA")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    docker = bin_dir / "docker"
    docker.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        'case "${1:-}" in\n'
        f"  image) printf '%s\\n' 'amd64|{expected.stdout.strip()}|2026-08-27' ;;\n"
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


def test_sif_runtime_identity_check_uses_embedded_exact_hash(tmp_path: Path):
    runtime_hash = "a" * 64
    engine = tmp_path / "singularity"
    engine.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "[[ \"${1:-}\" == inspect && \"${2:-}\" == --json ]]\n"
        "printf '%s\\n' "
        f"'{{\"data\":{{\"attributes\":{{\"labels\":{{\"io.genegalleon.runtime-input\":\"{runtime_hash}\",\"io.genegalleon.security-refresh-epoch\":\"2026-08-27\"}}}}}}}}'\n",
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
