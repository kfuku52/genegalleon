"""Container integration for NWKIT exploratory orthogroup selection."""

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytestmark = [pytest.mark.runtime, pytest.mark.integration]
SUPPORT = Path(__file__).resolve().parents[1] / "support"


def selection_module():
    spec = importlib.util.spec_from_file_location(
        "copy_selection", SUPPORT / "orthogroup_copy_number_trait_selection.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_missing_species_preserves_original_brownian_covariance(tmp_path):
    from ete4 import Tree
    from nwkit.ordinary_regression import build_phylogenetic_covariance

    tree = Tree("(((a:1.123456789012345,b:1):1,(c:1,d:1):1):2,(e:1,f:1):3);", parser=1)
    leaves = ["a", "b", "c", "d"]
    expected = build_phylogenetic_covariance(tree, leaves)
    path = tmp_path / "subset.nwk"
    selection_module().write_species_subset(tree, set(leaves), path)
    actual = build_phylogenetic_covariance(Tree(str(path), parser=1), leaves)
    np.testing.assert_allclose(actual, expected, rtol=1e-15, atol=0)
    assert len(list(tree.leaves())) == 6


@pytest.mark.parametrize("text", ["species\tx\na\t1\t2\n", "species\t\na\t1\n", "species\tleaf_name\na\t1\n"])
def test_malformed_species_table_is_rejected(tmp_path, text):
    path = tmp_path / "bad.tsv"
    path.write_text(text)
    with pytest.raises(ValueError, match="fields|non-empty|conflicts"):
        selection_module().load_species_table(path)


def test_copy_selection_mixed_families_and_rollback(tmp_path):
    rng = np.random.default_rng(52)
    leaves = [f"sp{i}" for i in range(18)]
    tree = tmp_path / "tree.nwk"
    tree.write_text("[&R] (" + ",".join(f"{leaf}:1" for leaf in leaves) + ");")
    counts = rng.poisson(2, (3, 18))
    data = pd.DataFrame(counts, columns=leaves)
    data.insert(0, "Orthogroup", ["OG1", "OG2", "OG3"])
    data.insert(0, "besthit_0.95", ["hit1", "hit2", "hit3"])
    count_path = tmp_path / "copy.tsv"
    data.to_csv(count_path, sep="\t", index=False)
    traits = tmp_path / "traits.tsv"
    trait_data = pd.DataFrame(
        {
            "species": leaves,
            "size": counts[0] + rng.normal(size=18),
            "state": np.tile([0, 1], 9),
            "copies": rng.poisson(3, 18),
        }
    )
    trait_data.to_csv(traits, sep="\t", index=False)
    folds = tmp_path / "folds.tsv"
    pd.DataFrame({"leaf_name": leaves, "fold": np.repeat(["cladeA", "cladeB", "cladeC"], 6)}).to_csv(
        folds, sep="\t", index=False
    )
    output = tmp_path / "out"
    command = [
        sys.executable,
        str(SUPPORT / "orthogroup_copy_number_trait_selection.py"),
        "--copy-number",
        str(count_path),
        "--tree",
        str(tree),
        "--traits",
        str(traits),
        "--folds",
        str(folds),
        "--outdir",
        str(output),
        "--strengths",
        "0.5",
        "--l1-ratios",
        "0.5",
        "--response-families",
        "state=binomial,copies=negative-binomial",
    ]
    result = subprocess.run(command, capture_output=True, text=True, env=os.environ.copy())
    assert result.returncode == 0, result.stderr
    manifest = pd.read_csv(output / "manifest.tsv", sep="\t")
    assert manifest.response_family.tolist() == ["gaussian", "binomial", "negative-binomial"]
    assert manifest.predictor_transform.eq("log1p").all()
    assert len(list(output.glob("*.predictions.tsv"))) == 3
    for path in output.glob("*.coefficients.tsv"):
        coefficients = pd.read_csv(path, sep="\t")
        assert not any("pval" in name or "p_value" in name for name in coefficients)
        assert coefficients.predictor_transform.eq("log1p").all()
        selected = coefficients[coefficients.Orthogroup.notna()]
        np.testing.assert_allclose(selected.training_center, np.log1p(counts).mean(axis=1))
    before = {file.name: file.read_bytes() for file in output.iterdir()}
    # Last trait fails only after the first two fits: nothing may be published.
    trait_data["copies"] = -1
    trait_data.to_csv(traits, sep="\t", index=False)
    failed = subprocess.run(command, capture_output=True, text=True)
    assert failed.returncode != 0
    assert {file.name: file.read_bytes() for file in output.iterdir()} == before


def test_core_selection_only_uses_imported_artifacts_and_mixed_run_requires_sequences(tmp_path):
    import re
    import shlex

    workflow = SUPPORT.parent
    entry = (workflow / "gg_genome_evolution_entrypoint.sh").read_text()
    config = entry.split("### Start: Modify this block to tailor your analysis ###", 1)[1].split(
        "### End: Modify this block to tailor your analysis ###", 1
    )[0]
    rng = np.random.default_rng(28)
    leaves = [f"sp{i}" for i in range(18)]
    tree_dir = tmp_path / "output/species_tree/species_tree_summary"
    counts_dir = tmp_path / "output/orthofinder/Orthogroups_filtered"
    trait_dir = tmp_path / "input/species_trait"
    for directory in (tree_dir, counts_dir, trait_dir, tmp_path / "input/species_cds"):
        directory.mkdir(parents=True)
    clades = ["(" + ",".join(f"{s}:1" for s in leaves[j : j + 6]) + "):1" for j in (0, 6, 12)]
    (tree_dir / "dated_species_tree.nwk").write_text(f"({clades[0]},({clades[1]},{clades[2]}):1);")
    counts = rng.poisson(2, (3, 18))
    data = pd.DataFrame(counts, columns=leaves)
    data.insert(0, "Orthogroup", ["OG1", "OG2", "OG3"])
    data.insert(0, "besthit_0.95", "")
    data.to_csv(counts_dir / "Orthogroups.GeneCount.selected.tsv", sep="\t", index=False)
    pd.DataFrame({"species": leaves, "height": counts[0] + rng.normal(size=18)}).to_csv(
        trait_dir / "species_trait.tsv", sep="\t", index=False
    )
    pd.DataFrame({"leaf_name": leaves, "fold": np.repeat(["a", "b", "c"], 6)}).to_csv(
        tmp_path / "input/folds.tsv", sep="\t", index=False
    )
    disabled = "\n".join(f"{name}=0" for name in re.findall(r"^(run_\w+)=", entry, re.M))
    script = f"""
set -euo pipefail
set -a
{config}
{disabled}
run_orthogroup_copy_number_trait_selection=1
orthogroup_copy_number_trait=height
orthogroup_copy_number_trait_selection_folds=input/folds.tsv
orthogroup_copy_number_trait_selection_strengths=0.5
orthogroup_copy_number_trait_selection_l1_ratios=0.5
species_tree_output_storage=files
gg_workflow_dir={shlex.quote(str(workflow))}
gg_support_dir={shlex.quote(str(SUPPORT))}
gg_workspace_dir={shlex.quote(str(tmp_path))}
gg_container_image_path=/unused-container-already-running.sif
GG_TASK_CPUS=1
GG_COMMON_TMP_ROOT=workspace
${{GG_TEST_MIXED_COMMAND:-:}}
exec bash {shlex.quote(str(workflow / "core/gg_genome_evolution_core.sh"))}
"""
    result = subprocess.run(["bash", "-c", script], cwd=workflow.parent, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    manifest = tmp_path / "output/genome_evolution/orthogroup_copy_number/trait_selection/manifest.tsv"
    assert pd.read_csv(manifest, sep="\t").n_species.tolist() == [18]
    assert "Copy-number selection only" in result.stdout
    assert not list((tmp_path / "input/species_cds").iterdir())
    # Enabling a sequence-consuming stage must restore the normal prerequisite.
    mixed = script.replace("${GG_TEST_MIXED_COMMAND:-:}", "run_species_busco=1")
    failed = subprocess.run(["bash", "-c", mixed], cwd=workflow.parent, capture_output=True, text=True)
    assert failed.returncode != 0
    assert "No species_cds fasta files" in failed.stdout
    assert manifest.exists()
