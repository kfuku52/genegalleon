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


def test_native_wastral_concordant_trees_have_finite_branches(tmp_path: Path):
    import math

    from nwkit import util

    gene_trees = tmp_path / "genes.nwk"
    # Fully concordant cherries exposed upstream NaN coalescent estimates.
    gene_trees.write_text("((A:1,B:1)1:1,((C:1,D:1)1:1,(E:1,F:1)1:1)1:1);\n" * 20)
    output = tmp_path / "species.nwk"
    completed = subprocess.run(
        ["astral-hybrid", "--input", str(gene_trees), "--output", str(output),
         "--mode", "3", "--support", "2", "--thread", "2"],
        capture_output=True, text=True, timeout=60, check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert "Weighted ASTRAL" in completed.stderr, "A Java fallback changes the estimator."
    for label in ("q1", "pp1"):
        labeled = tmp_path / f"species.{label}.nwk"
        converter = Path(__file__).resolve().parents[1] / "support/extract_astral_support_label.py"
        converted = subprocess.run(
            [sys.executable, str(converter), "--infile", str(output), "--outfile", str(labeled),
             "--label_key", label], capture_output=True, text=True, timeout=30, check=False,
        )
        assert converted.returncode == 0, converted.stderr
        assert "did not contain key" not in converted.stderr
        tree = util.read_tree(str(labeled), "auto", True, quiet=True)
        assert set(tree.leaf_names()) == set("ABCDEF")
        assert all(node.dist is None or math.isfinite(node.dist) for node in tree.traverse())
        tips = frozenset(tree.leaf_names())
        clades = {frozenset(node.leaf_names()) for node in tree.traverse()}
        clades |= {tips - clade for clade in clades}
        for cherry in ("AB", "CD", "EF"):
            assert frozenset(cherry) in clades
