import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

BIOCONDA_FIX_URL = "https://github.com/bioconda/bioconda-recipes/pull/68643"


def _conda_prefix() -> Path:
    return Path(os.environ.get("CONDA_PREFIX", sys.prefix)).resolve()


def test_aster_does_not_install_loader_path_activation_hooks():
    prefix = _conda_prefix()
    hooks = sorted(
        path
        for directory in ("activate.d", "deactivate.d")
        for path in (prefix / "etc" / "conda" / directory).glob("aster_*.sh")
    )

    assert not hooks, (
        "ASTER installed legacy activation hooks that override the runtime "
        "library search path and can break OR-Tools. Use the hook-free "
        f"Bioconda package tracked at {BIOCONDA_FIX_URL}; found: {hooks}"
    )


def test_aster_and_ortools_scip_load_in_the_same_runtime():
    astral_hybrid = shutil.which("astral-hybrid")
    assert astral_hybrid, "The GeneGalleon runtime must provide astral-hybrid."
    subprocess.run(
        [astral_hybrid, "--help"],
        check=True,
        capture_output=True,
        text=True,
        timeout=30,
    )

    try:
        from ortools.linear_solver import pywraplp
    except ImportError as exc:
        pytest.fail(
            "OR-Tools failed to import beside ASTER. Check for the ASTER 1.25 "
            "build-0 loader hook conflict before changing GeneGalleon or "
            f"kfFractBias dependencies: {BIOCONDA_FIX_URL}\n{exc}",
            pytrace=False,
        )

    solver = pywraplp.Solver.CreateSolver("SCIP")
    assert solver is not None, "OR-Tools could not create its bundled SCIP solver."
