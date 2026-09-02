import json
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "container/scripts/ensure_aster_loader_compatibility.sh"
ACTIVATE = (
    "export OLD_LD_LIBRARY_PATH=${LD_LIBRARY_PATH}\n"
    "export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib/:${LD_LIBRARY_PATH}\n"
    "export SINGULARITYENV_LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"
)
DEACTIVATE = (
    "export LD_LIBRARY_PATH=${OLD_LD_LIBRARY_PATH}\n"
    "unset OLD_LD_LIBRARY_PATH\n"
    "unset SINGULARITYENV_LD_LIBRARY_PATH"
)


def make_prefix(tmp_path: Path, *, build_number: int = 0) -> Path:
    prefix = tmp_path / "conda"
    (prefix / "bin").mkdir(parents=True)
    os.symlink(sys.executable, prefix / "bin/python")
    activate = prefix / "etc/conda/activate.d/aster_activate.sh"
    deactivate = prefix / "etc/conda/deactivate.d/aster_deactivate.sh"
    activate.parent.mkdir(parents=True)
    deactivate.parent.mkdir(parents=True)
    activate.write_text(ACTIVATE, encoding="utf-8")
    deactivate.write_text(DEACTIVATE, encoding="utf-8")
    metadata = prefix / "conda-meta/aster-1.25-test.json"
    metadata.parent.mkdir(parents=True)
    metadata.write_text(
        json.dumps(
            {
                "name": "aster",
                "version": "1.25",
                "build": f"test_{build_number}",
                "build_number": build_number,
            }
        ),
        encoding="utf-8",
    )
    return prefix


def run_correction(prefix: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["bash", str(SCRIPT), str(prefix)],
        check=False,
        capture_output=True,
        text=True,
    )


def test_removes_only_verified_aster_1_25_build_0_hooks(tmp_path: Path):
    prefix = make_prefix(tmp_path)

    completed = run_correction(prefix)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert not (prefix / "etc/conda/activate.d/aster_activate.sh").exists()
    assert not (prefix / "etc/conda/deactivate.d/aster_deactivate.sh").exists()
    log = (tmp_path / "pg/logs/aster_dependency_correction.txt").read_text()
    assert "status=removed legacy loader hooks" in log
    assert "bioconda-recipes/pull/68643" in log


def test_hook_free_package_is_unchanged(tmp_path: Path):
    prefix = tmp_path / "conda"
    (prefix / "bin").mkdir(parents=True)
    os.symlink(sys.executable, prefix / "bin/python")

    completed = run_correction(prefix)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "no correction needed" in completed.stdout
    assert not (tmp_path / "pg/logs/aster_dependency_correction.txt").exists()


def test_rejects_changed_hook_content_without_removing_files(tmp_path: Path):
    prefix = make_prefix(tmp_path)
    changed = prefix / "etc/conda/activate.d/aster_activate.sh"
    changed.write_text(ACTIVATE + "\nexport UNEXPECTED=1\n", encoding="utf-8")

    completed = run_correction(prefix)

    assert completed.returncode != 0
    assert "Refusing to modify changed ASTER hook" in completed.stderr
    assert changed.exists()
    assert (prefix / "etc/conda/deactivate.d/aster_deactivate.sh").exists()


def test_rejects_unrecognized_package_build_without_removing_files(tmp_path: Path):
    prefix = make_prefix(tmp_path, build_number=1)

    completed = run_correction(prefix)

    assert completed.returncode != 0
    assert "unrecognized package identity" in completed.stderr
    assert (prefix / "etc/conda/activate.d/aster_activate.sh").exists()
    assert (prefix / "etc/conda/deactivate.d/aster_deactivate.sh").exists()
