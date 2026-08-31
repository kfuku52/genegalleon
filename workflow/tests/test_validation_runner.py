import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import prepare_test_wheels
import pytest
import run_checks

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_focused_validation_keeps_explicit_file_and_pytest_arguments():
    target = "workflow/tests/test_shared_namespace_lock.py"
    result = subprocess.run(
        [sys.executable, str(Path(run_checks.__file__)), "fast", "--list",
         "--workers", "3", "--", target, "-k", "release and not timeout", "-x"],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    commands = json.loads(result.stdout)
    assert len(commands) == 1
    assert commands[0][-4:] == [target, "-k", "release and not timeout", "-x"]
    assert "workflow/tests" not in commands[0]
    assert commands[0][commands[0].index("-n") + 1] == "3"


@pytest.mark.parametrize("suite", ["static", "fast"])
def test_validation_loads_suite_options_without_an_explicit_test_path(suite):
    result = subprocess.run(
        [sys.executable, str(Path(run_checks.__file__)), suite, "--collect-only", "-p", "no:cacheprovider"],
        cwd=REPO_ROOT, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert ("workflow/tests/test_ci_static.py::" in result.stdout) == (suite == "static")
    assert ("workflow/tests/test_generate_orthogroup_database.py::" in result.stdout) == (suite == "fast")


def test_dev_forwards_focused_arguments_through_the_container(tmp_path):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    fake_docker = bin_dir / "docker"
    fake_docker.write_text('#!/bin/bash\nif [[ "$1" == run ]]; then printf "%s\\n" "$@"; fi\n')
    fake_docker.chmod(0o755)
    env = os.environ | {
        "PATH": f"{bin_dir}:/usr/bin:/bin", "GG_TEST_RUNTIME": "docker",
        "GG_CONTAINER_DOCKER_IMAGE": "local/genegalleon:test", "GG_RUNTIME_FRESHNESS": "off",
    }
    extra = ["workflow/tests/test_shared_namespace_lock.py", "-k", "release and not timeout", "-x"]
    result = subprocess.run(["bash", str(REPO_ROOT / "dev"), "check", "fast", *extra],
                            env=env, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    arguments = result.stdout.splitlines()
    assert arguments[-len(extra):] == extra
    assert str(REPO_ROOT / "workflow/tests/run_checks.py") in arguments
    assert "PYTHONDONTWRITEBYTECODE=1" in arguments


def test_full_runtime_and_r_checks_include_every_r_test():
    manifest = run_checks.load_manifest()
    declared = manifest["r_commands"]
    r_files = {str(path.relative_to(REPO_ROOT)) for path in (REPO_ROOT / "workflow/tests").glob("test_*.R")}
    assert {command[1] for command in declared if command[0] == "Rscript"} == r_files
    assert ["bash", "workflow/tests/check_treevis_package.sh"] in declared
    for suite in ("full", "runtime", "r"):
        commands = run_checks.commands_for(suite, "2", [])
        assert all(command in commands for command in declared)
    assert "test_fractionation_bias_integration.py" in manifest["runtime_python_files"]
    assert manifest["environment"]["KFFRACTBIAS_RUN_INTEGRATION"] == "1"


@pytest.mark.parametrize("workers", [None, "2"])
@pytest.mark.parametrize("behavior", ["skip", "skip_collection", "check_environment"])
def test_strict_runtime_detects_skips_and_enables_required_integrations(tmp_path, workers, behavior):
    test_dir = tmp_path / "workflow/tests"
    test_dir.mkdir(parents=True)
    for name in ("conftest.py", "validation_manifest.json"):
        shutil.copyfile(REPO_ROOT / "workflow/tests" / name, test_dir / name)
    body = {
        "skip": 'import pytest\ndef test_runtime(): pytest.skip("missing required executable")\n',
        "skip_collection": 'import pytest\npytest.skip("missing dependency", allow_module_level=True)\n',
        "check_environment": 'import os\ndef test_runtime(): assert os.environ["KFFRACTBIAS_RUN_INTEGRATION"] == "1"\n',
    }[behavior]
    (test_dir / "test_required.py").write_text(body)
    env = os.environ.copy()
    env.pop("KFFRACTBIAS_RUN_INTEGRATION", None)
    command = [sys.executable, "-m", "pytest", "-q", "--gg-strict-runtime",
               "-c", str(REPO_ROOT / "pyproject.toml"), "--rootdir", str(tmp_path),
               "--confcutdir", str(tmp_path), str(test_dir)]
    if workers:
        command += ["-n", workers]
    result = subprocess.run(command, cwd=tmp_path, env=env, capture_output=True, text=True, check=False)
    if behavior == "check_environment":
        assert result.returncode == 0, result.stdout + result.stderr
    else:
        assert result.returncode != 0
        assert "Required runtime checks were skipped" in result.stdout


def test_wheel_inputs_resolve_only_the_temporary_build_and_install_offline(tmp_path):
    requirements = REPO_ROOT / "workflow/tests/requirements.txt"
    before = requirements.read_bytes()
    source = "b" * 40
    prepare_test_wheels.prepare(tmp_path, source)
    assert requirements.read_bytes() == before
    build = (tmp_path / "build-requirements.txt").read_text()
    install = (tmp_path / "install-requirements.txt").read_text()
    assert f"csubst.git@{source}" in build
    assert "git+" not in install
    assert "csubst" in install.splitlines()
    assert json.loads((tmp_path / "source-identity.json").read_text())["csubst_sha"] == source
    with pytest.raises(ValueError, match="resolved"):
        prepare_test_wheels.prepare(tmp_path, "master")
