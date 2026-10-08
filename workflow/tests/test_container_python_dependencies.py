"""Guard mandatory source-wheel dependencies without a container build or network."""
import os
import re
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


def _wheel(directory, name, version, requires=()):
    path = directory / f"{name}-{version}-py3-none-any.whl"
    dist = f"{name}-{version}.dist-info"
    metadata = f"Metadata-Version: 2.1\nName: {name}\nVersion: {version}\n"
    metadata += "".join(f"Requires-Dist: {requirement}\n" for requirement in requires)
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(f"{name}/__init__.py", "")
        archive.writestr(f"{dist}/METADATA", metadata)
        archive.writestr(f"{dist}/WHEEL", "Wheel-Version: 1.0\nGenerator: fixture\nRoot-Is-Purelib: true\nTag: py3-none-any\n")
        archive.writestr(f"{dist}/RECORD", "")
    return path


def _pip(python, *arguments):
    env = {**os.environ, "PIP_CONFIG_FILE": os.devnull, "PIP_NO_INDEX": "1", "PYTHONNOUSERSITE": "1"}
    env.pop("PYTHONPATH", None)
    return subprocess.run(
        [str(python), "-m", "pip", "--disable-pip-version-check", *arguments],
        env=env, text=True, capture_output=True, timeout=60,
    )


@pytest.fixture
def isolated_python(tmp_path):
    target = tmp_path / "venv"
    completed = subprocess.run(
        [sys.executable, "-m", "venv", str(target)], text=True, capture_output=True, timeout=60,
    )
    assert completed.returncode == 0, completed.stderr
    return target / ("Scripts/python.exe" if os.name == "nt" else "bin/python")


def test_nwkit_required_excel_reader_is_declared_for_both_builds():
    requirements = (REPO_ROOT / "container/pip-compatibility.requirements.txt").read_text().splitlines()
    excel_requirements = [line.strip() for line in requirements if re.match(r"^xlrd(?:[<=>!~]|$)", line.strip())]
    assert excel_requirements == ["xlrd>=2"]
    for relative in ("container/Dockerfile", "container/apptainer_local_build.def.template"):
        build = (REPO_ROOT / relative).read_text()
        assert "-r /opt/pg/pip-compatibility.requirements.txt" in build
        assert "python -m pip check" in build


def test_source_wheels_keep_dependency_resolution_off_and_final_check_on():
    installer = (REPO_ROOT / "container/scripts/install_source_artifacts.sh").read_text()
    assert "--no-index --no-deps --force-reinstall" in installer
    assert installer.index("--no-index --no-deps") < installer.index("python -m pip check")
    builder = (REPO_ROOT / "container/scripts/build_source_artifact.sh").read_text()
    assert "--no-deps --no-build-isolation" in builder


def test_source_dependency_check_rejects_missing_and_obsolete_xlrd_then_accepts_install(
    isolated_python, tmp_path,
):
    # This packaging-only upstream fixture declares the exact mandatory Excel
    # dependency. It does not execute NWKIT or build a scientific runtime.
    upstream = _wheel(tmp_path, "nwkit", "0.43.43", ("xlrd>=2",))
    installed = _pip(isolated_python, "install", "--no-deps", "--no-index", str(upstream))
    assert installed.returncode == 0, installed.stderr
    missing = _pip(isolated_python, "check")
    assert missing.returncode != 0
    assert "requires xlrd, which is not installed" in missing.stdout.lower()

    obsolete = _wheel(tmp_path, "xlrd", "1.2.0")
    assert _pip(isolated_python, "install", "--no-deps", "--no-index", str(obsolete)).returncode == 0
    incompatible = _pip(isolated_python, "check")
    assert incompatible.returncode != 0
    assert "xlrd>=2" in incompatible.stdout.lower()
    assert "xlrd 1.2.0" in incompatible.stdout.lower()

    current = _wheel(tmp_path, "xlrd", "2.0.2")
    satisfied = _pip(isolated_python, "install", "--no-deps", "--no-index", str(current))
    assert satisfied.returncode == 0, satisfied.stderr
    valid = _pip(isolated_python, "check")
    assert valid.returncode == 0, valid.stdout + valid.stderr
    assert "No broken requirements found" in valid.stdout
