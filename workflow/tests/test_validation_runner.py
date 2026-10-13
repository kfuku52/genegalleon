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
@pytest.mark.parametrize("shard", [None, "1/1"])
def test_validation_loads_suite_options_without_an_explicit_test_path(tmp_path, suite, shard):
    # Exercise real testpaths/conftest discovery without importing the entire
    # scientific suite a second time inside this test.
    test_dir = tmp_path / "workflow/tests"
    test_dir.mkdir(parents=True)
    for name in ("conftest.py", "validation_manifest.json", "run_checks.py"):
        shutil.copyfile(REPO_ROOT / "workflow/tests" / name, test_dir / name)
    shutil.copyfile(REPO_ROOT / "pyproject.toml", tmp_path / "pyproject.toml")
    for lane in ("static", "fast"):
        (test_dir / f"test_sentinel_{lane}.py").write_text("def test_sentinel(): pass\n")
    command = [sys.executable, str(test_dir / "run_checks.py"), suite, "--collect-only", "-p", "no:cacheprovider"]
    if shard:
        command += ["--gg-shard", shard]
    result = subprocess.run(
        command,
        cwd=tmp_path, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    for lane in ("static", "fast"):
        assert (f"workflow/tests/test_sentinel_{lane}.py::" in result.stdout) == (suite == lane)
    assert "1 test collected" in result.stdout


@pytest.mark.parametrize("suite", ["fast", "smoke"])
@pytest.mark.parametrize("option", ["--ignore", "--ignore-glob"])
def test_suite_selection_honors_pytest_ignore_options(tmp_path, suite, option):
    test_dir = tmp_path / "workflow/tests"
    test_dir.mkdir(parents=True)
    for name in ("conftest.py", "validation_manifest.json"):
        shutil.copyfile(REPO_ROOT / "workflow/tests" / name, test_dir / name)
    if suite == "smoke":
        keep_name = "test_busco_hmmsearch_wrapper.py"
        keep_test = "test_hmmsearch_wrapper_creates_modified_fas_symlink_when_missing"
        ignore_name = "test_gg_input_generation_end_to_end.py"
    else:
        keep_name, keep_test, ignore_name = "test_keep.py", "test_keep", "test_ignore.py"
    (test_dir / keep_name).write_text(f"def {keep_test}(): pass\n")
    # An ignored module must not even be imported during collection.
    (test_dir / ignore_name).write_text('raise RuntimeError("ignored module was imported")\n')
    ignored = str(test_dir / ignore_name) if option == "--ignore" else f"*{ignore_name}"
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "--collect-only", f"--gg-suite={suite}",
         option, ignored, "-c", str(REPO_ROOT / "pyproject.toml"),
         "--rootdir", str(tmp_path), "--confcutdir", str(tmp_path), str(test_dir)],
        cwd=tmp_path, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert keep_test in result.stdout
    assert "1 test collected" in result.stdout


@pytest.mark.parametrize("suite, expected_count", [("fast", 4), ("runtime", 6)])
def test_promoter_unit_and_tool_checks_are_collected_in_separate_lanes(tmp_path, suite, expected_count):
    test_dir = tmp_path / "workflow/tests"
    test_dir.mkdir(parents=True)
    for name in ("conftest.py", "validation_manifest.json", "test_get_promoter_fasta.py",
                 "test_get_promoter_fasta_runtime.py"):
        shutil.copyfile(REPO_ROOT / "workflow/tests" / name, test_dir / name)
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "--collect-only", f"--gg-suite={suite}",
         "-c", str(REPO_ROOT / "pyproject.toml"), "--rootdir", str(tmp_path),
         "--confcutdir", str(tmp_path), str(test_dir)],
        cwd=tmp_path, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert ("test_normalize_ncpu_clamps_non_positive_values" in result.stdout) == (suite == "fast")
    assert ("test_promoter_cli_preserves_chromosome_identifiers" in result.stdout) == (suite == "runtime")
    assert f"{expected_count} tests collected" in result.stdout


def test_dev_forwards_focused_arguments_through_the_container(tmp_path):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    fake_docker = bin_dir / "docker"
    fake_docker.write_text('#!/bin/bash\nif [[ "$1" == run ]]; then printf "%s\\n" "$@"; '
                           'elif [[ "$1" == image ]]; then printf "sha256:%064d\\n" 0; fi\n')
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


def test_split_runtime_commands_preserve_the_complete_validation_contract():
    complete = run_checks.commands_for("runtime", "2", [])
    python = run_checks.commands_for("runtime-python", "2", [])
    extra = run_checks.commands_for("runtime-extra", "2", [])
    r = run_checks.commands_for("r", "2", [])
    assert python[0][python[0].index("-n"):python[0].index("-n") + 4] == ["-n", "2", "--dist", "load"]
    serial_python = [argument for argument in python[0] if argument not in ("-n", "2", "--dist", "load")]
    assert [serial_python] + extra + r == complete
    assert run_checks.commands_for("runtime-python", "2", ["--gg-shard", "3/8"])[0] == (
        python[0] + ["--gg-shard", "3/8"]
    )


def test_shards_partition_the_selected_lane_without_missing_or_repeating_cases(tmp_path):
    test_dir = tmp_path / "workflow/tests"
    test_dir.mkdir(parents=True)
    for name in ("conftest.py", "validation_manifest.json"):
        shutil.copyfile(REPO_ROOT / "workflow/tests" / name, test_dir / name)
    (test_dir / "test_gift_retrieval.py").write_text(
        'import pytest\n@pytest.mark.parametrize("case", range(64))\ndef test_case(case): pass\n'
    )
    (test_dir / "test_other.py").write_text('raise RuntimeError("wrong lane imported")\n')
    command = [sys.executable, "-m", "pytest", "-q", "--collect-only", "--gg-suite=runtime",
               "-c", str(REPO_ROOT / "pyproject.toml"), "--rootdir", str(tmp_path),
               "--confcutdir", str(tmp_path), str(test_dir)]

    def collect(extra):
        result = subprocess.run(command + extra, cwd=tmp_path, capture_output=True, text=True, check=False)
        assert result.returncode == 0, result.stdout + result.stderr
        return {line for line in result.stdout.splitlines() if "::test_case[" in line}

    complete = collect([])
    partitions = [collect(["--gg-shard", f"{index}/4"]) for index in range(1, 5)]
    assert len(complete) == 64
    assert set.union(*partitions) == complete
    assert sum(map(len, partitions)) == len(complete)
    assert collect(["--gg-shard", "2/4"]) == partitions[1]


@pytest.mark.parametrize("shard", ["0/8", "9/8", "1/0", "bad"])
def test_invalid_shard_cannot_silently_run_partial_validation(tmp_path, shard):
    test_dir = tmp_path / "workflow/tests"
    test_dir.mkdir(parents=True)
    for name in ("conftest.py", "validation_manifest.json"):
        shutil.copyfile(REPO_ROOT / "workflow/tests" / name, test_dir / name)
    (test_dir / "test_sentinel.py").write_text("def test_sentinel(): pass\n")
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "--gg-shard", shard,
         "-c", str(REPO_ROOT / "pyproject.toml"), "--rootdir", str(tmp_path),
         "--confcutdir", str(tmp_path), str(test_dir)],
        cwd=tmp_path, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 4
    assert "--gg-shard must be INDEX/COUNT" in result.stderr


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
        "check_environment": 'import os\ndef test_runtime(): assert os.environ["KFFRACTBIAS_RUN_INTEGRATION"] == "1"; assert os.environ["GG_TEST_CSUBST_3DI"] == "1"\n',
    }[behavior]
    (test_dir / "test_required.py").write_text(body)
    env = os.environ.copy()
    env.pop("KFFRACTBIAS_RUN_INTEGRATION", None)
    env.pop("GG_TEST_CSUBST_3DI", None)
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
    nwkit_source = "c" * 40
    cdskit_source = "d" * 40
    kffractbias_source = "e" * 40
    prepare_test_wheels.prepare(tmp_path, source, nwkit_source, cdskit_source, kffractbias_source)
    assert requirements.read_bytes() == before
    build = (tmp_path / "build-requirements.txt").read_text()
    install = (tmp_path / "install-requirements.txt").read_text()
    assert f"csubst.git@{source}" in build
    assert f"nwkit.git@{nwkit_source}" in build
    assert f"cdskit.git@{cdskit_source}" in build
    assert f"kfFractBias.git@{kffractbias_source}" in build
    assert "cdskit" in install.splitlines()
    assert "nwkit" in install.splitlines()
    assert "git+" not in install
    assert "csubst" in install.splitlines()
    assert "kffractbias" in install.splitlines()
    assert json.loads((tmp_path / "source-identity.json").read_text())["csubst_sha"] == source
    assert json.loads((tmp_path / "source-identity.json").read_text())["nwkit_sha"] == nwkit_source
    assert json.loads((tmp_path / "source-identity.json").read_text())["kffractbias_sha"] == kffractbias_source
    with pytest.raises(ValueError, match="resolved"):
        prepare_test_wheels.prepare(tmp_path, "master", nwkit_source, cdskit_source, kffractbias_source)
    with pytest.raises(ValueError, match="resolved"):
        prepare_test_wheels.prepare(tmp_path, source, "master", cdskit_source, kffractbias_source)
    with pytest.raises(ValueError, match="resolved"):
        prepare_test_wheels.prepare(tmp_path, source, nwkit_source, cdskit_source, "master")
