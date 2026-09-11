"""Run the real in-place dating stage and downstream artifact contracts."""

import json
import os
import shlex
import subprocess
import time
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[2]


def _fixture(tmp_path, model, route, bad_codon=False, use_defaults=False, engine="native", iqtree_model="", iqtree_interface="auto"):

    from nwkit.radte_codon import CODONS, codon_matrix

    tree = "((A_1:0.12[&&NHX:S=A:D=N],A_2:0.15[&&NHX:S=A:D=N])nA:0.17[&&NHX:S=A:D=Y],(B_1:0.13[&&NHX:S=B:D=N],B_2:0.11[&&NHX:S=B:D=N])nB:0.18[&&NHX:S=B:D=Y])Root:0[&&NHX:S=0:D=N];"
    (tmp_path / "gene.nhx").write_text(tree)
    (tmp_path / "species.nwk").write_text("(A:10,B:10)Root;")
    (tmp_path / "generax.nwk").write_text("(B:0.1,A:0.2)0;")
    (tmp_path / "species_intervals.tsv").write_text(
        "node\tlower\tupper\tlevel\tkind\tsource\nRoot\t8\t12\t0.95\tconfidence\texternal-study\n"
    )
    (tmp_path / "species_map.tsv").write_text("leaf_name\tspecies_label\nA_1\tA\nA_2\tA\nB_1\tB\nB_2\tB\n")
    from nwkit.reconcile import build_reconciliation_table
    from nwkit.util import read_tree
    gene = read_tree(str(tmp_path / "gene.nhx"), "auto", True)
    species = read_tree(str(tmp_path / "species.nwk"), "auto", True)
    build_reconciliation_table(gene, species, {name: name.split("_")[0] for name in gene.leaf_names()},
                               tree_id="OG0000001").to_csv(tmp_path / "reconciliation.tsv", sep="\t", index=False)
    q, pi = codon_matrix(model, np.ones((4, 1), dtype=np.uint64), np.ones(1), 2, 0.5, "fq")
    rng = np.random.default_rng(91)
    root = rng.choice(61, 150, p=pi)

    def descend(parent, length):
        p = np.maximum(expm(q * length), 0)
        p /= p.sum(axis=1, keepdims=True)
        return np.array([rng.choice(61, p=p[i]) for i in parent])

    left, right = descend(root, 0.17), descend(root, 0.18)
    sequences = {
        name: "".join(CODONS[i] for i in descend(parent, length))
        for name, parent, length in [
            ("A_1", left, 0.12),
            ("A_2", left, 0.15),
            ("B_1", right, 0.13),
            ("B_2", right, 0.11),
        ]
    }
    if bad_codon:
        sequences["A_1"] = "TAA" + sequences["A_1"][3:]
    (tmp_path / "alignment.fa").write_text("".join(f">{name}\n{seq}\n" for name, seq in sequences.items()))
    defaults = dict(
        radte_sequence_engine=engine,
        radte_iqtree_model=iqtree_model,
        radte_iqtree_interface=iqtree_interface,
        radte_substitution_model=model,
        radte_codon_frequencies="fq",
        radte_kappa="2" if model == "gy94" else "",
        radte_omega="0.5" if model == "gy94" else "",
        radte_gamma_shape="1",
        radte_gamma_categories="1",
        radte_inference="auto",
        radte_likelihood="auto",
        radte_uncertainty="profile",
        radte_interval_level="0.95",
        radte_rate_sd="0.3",
        radte_maxiter="500",
        radte_seed="1",
        radte_species_intervals_tsv=str(tmp_path / "species_intervals.tsv"),
        radte_max_age="30",
        genetic_code="1",
        input_sequence_mode="cds",
        run_tree_dating="1",
        run_generax="1" if route == "generax" else "0",
        species_tree_basename="dated_species_tree",
        species_tree_pruned=str(tmp_path / "species.nwk"),
        species_tree_generax=str(tmp_path / "generax.nwk"),
        file_og_generax_nhx=str(tmp_path / "gene.nhx"),
        file_og_reconciliation=str(tmp_path / "reconciliation.tsv"),
        file_og_rooted_tree_analysis=str(tmp_path / "gene.nhx"),
        file_og_unrooted_tree_analysis=str(tmp_path / "gene.nhx"),
        file_og_trimmed_aln_analysis=str(tmp_path / "alignment.fa"),
        file_og_dated_tree=str(tmp_path / "out/dated_tree/OG0000001_dated.nwk"),
        file_og_dated_tree_log=str(tmp_path / "out/dated_tree_log/OG0000001_dated.log.txt"),
        file_og_radte_prefix=str(tmp_path / "out/dated_tree_native/OG0000001_radte"),
        species_label_parser="legacy",
        species_label_regex="",
        species_label_map_tsv=str(tmp_path / "species_map.tsv"),
        file_species_genetic_code=str(tmp_path / "absent.tsv"),
        dir_output_active=str(tmp_path / "out"),
        gg_workspace_dir=str(tmp_path),
        gg_support_dir=str(ROOT / "workflow/support"),
        og_id="OG0000001",
        GG_TASK_CPUS="1",
    )
    if use_defaults:
        defaults.update(
            radte_substitution_model="auto",
            radte_codon_frequencies="",
            radte_kappa="",
            radte_omega="",
            radte_gamma_shape="",
            radte_gamma_categories="4",
            radte_rate_sd="",
        )
    if iqtree_model:
        defaults.update(radte_codon_frequencies="", radte_kappa="", radte_omega="", radte_gamma_shape="")
    core = (ROOT / "workflow/core/gg_gene_evolution_core.sh").read_text()
    start = core.index('task="Species-tree-guided divergence time estimation"')
    end = core.index('task="Expression matrix preparation"', start)
    # Only scheduler/reporting are stubbed. NWKIT, seqkit and bundle publication run.
    script = "set -euo pipefail\nsource " + shlex.quote(str(ROOT / "workflow/support/gg_util.sh")) + "\n"
    script += "\n".join(k + "=" + shlex.quote(v) for k, v in defaults.items()) + "\n"
    script += """
disable_if_no_input_file() { :; }
gg_step_start() { :; }
gg_step_skip() { :; }
gg_extract_expected_zip_prefix() { :; }
gg_artifact_prepare_stage() { printf -v "$1" 1; }
gg_artifact_record() { printf '%s\\n' "$@" > recorded-provenance.txt; }
"""
    script += core[start:end]
    started = time.monotonic()
    result = subprocess.run(
        ["bash", "-c", script],
        cwd=tmp_path,
        text=True,
        capture_output=True,
        env={**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"},
        timeout=300,
    )
    (tmp_path / "stage_wall_seconds.txt").write_text(str(time.monotonic() - started))
    return result


@pytest.mark.parametrize("model,route", [("gy94", "generax"), ("ecmk07", "generax"), ("ecmrest", "lca")])
def test_native_dating_stage_publishes_complete_codon_bundle(tmp_path, model, route):
    result = _fixture(tmp_path, model, route)
    assert result.returncode == 0, result.stdout + result.stderr
    prefix = tmp_path / "out/dated_tree_native/OG0000001_radte"
    manifest = json.loads(Path(str(prefix) + ".manifest.json").read_text())
    assert manifest["sequence_model"]["model"] == model
    assert manifest["options"]["backend"] == "native"
    assert manifest["uncertainty"].startswith("conditional-profile")
    assert manifest["method"] in {"sequence-marginal-quadratic", "sequence-empirical-bayes-map"}
    assert (tmp_path / "out/dated_tree/OG0000001_dated.nwk").read_bytes() == Path(
        str(prefix) + ".dated.nwk"
    ).read_bytes()
    assert json.loads((tmp_path / "out/dated_tree_log/OG0000001_dated.log.txt").read_text()) == manifest
    assert Path(str(prefix) + ".pdf").read_bytes().startswith(b"%PDF")
    assert "engine=nwkit-native-sequence-v1" in (tmp_path / "recorded-provenance.txt").read_text()
    import pandas as pd

    species = pd.read_csv(str(prefix) + ".species.tsv", sep="\t")
    root = species[species.node == "Root"].iloc[0]
    assert root.age == root.age_min == root.age_max == 10


def test_failed_codon_dating_keeps_previous_results(tmp_path):
    old = tmp_path / "out/dated_tree/OG0000001_dated.nwk"
    old.parent.mkdir(parents=True)
    old.write_text("previous validated tree")
    result = _fixture(tmp_path, "gy94", "generax", bad_codon=True)
    assert result.returncode != 0
    assert "Stop codon" in result.stdout + result.stderr
    assert old.read_text() == "previous validated tree"
    assert not (tmp_path / "recorded-provenance.txt").exists()


def _assert_profile_output_contract(prefix, manifest):
    import pandas as pd

    nodes = pd.read_csv(str(prefix) + ".nodes.tsv", sep="\t")
    if manifest["log_rate_sd"] == 0:
        # A boundary optimum cannot be reported as a regular profile interval.
        assert manifest["uncertainty"] == "unavailable-strict-clock-limit"
        assert "estimated_rate_variance_at_zero_boundary" in manifest["diagnostics"]
        assert nodes[["interval_lower", "interval_upper"]].isna().all().all()
    else:
        assert manifest["uncertainty"].startswith("conditional-profile")


def test_default_cds_model_estimates_nuisance_parameters(tmp_path):
    result = _fixture(tmp_path, "gy94", "generax", use_defaults=True)
    assert result.returncode == 0, result.stdout + result.stderr
    manifest = json.loads((tmp_path / "out/dated_tree_native/OG0000001_radte.manifest.json").read_text())
    model = manifest["sequence_model"]
    assert model["model"] == "gy94"
    assert model["codon_frequencies"] == "f3x4"
    assert model["prefit"]["status"] == "estimated-unclocked-conditional-model"
    assert model["kappa"] > 0 and model["omega"] > 0 and model["gamma_shape"] > 0
    assert model["alignment_codon_sites"] == 150
    assert manifest["options"]["uncertainty"] == "profile"
    assert manifest["options"]["interval_level"] == 0.95
    _assert_profile_output_contract(tmp_path / "out/dated_tree_native/OG0000001_radte", manifest)
    if manifest["log_rate_sd"] > 0:
        assert "profile_quadratic_failed_validation_refitted_exact" in manifest["diagnostics"]


@pytest.mark.parametrize("iqtree_model", ["", "GY+F3X4+R4"])
def test_iqtree_stage_retains_profile_and_species_age_contract(tmp_path, iqtree_model):
    result = _fixture(tmp_path, "gy94", "generax", engine="iqtree", iqtree_model=iqtree_model, use_defaults=not iqtree_model, iqtree_interface="cli")
    assert result.returncode == 0, result.stdout + result.stderr
    prefix = tmp_path / "out/dated_tree_native/OG0000001_radte"
    manifest = json.loads(Path(str(prefix) + ".manifest.json").read_text())
    assert manifest["sequence_model"]["engine"] == "iqtree"
    assert manifest["sequence_model"]["iqtree_interface"] == "standard-cli-iq2mc"
    assert manifest["sequence_model"]["frozen_model"].startswith("GY{")
    _assert_profile_output_contract(tmp_path / "out/dated_tree_native/OG0000001_radte", manifest)
    assert manifest["options"]["backend"] == "native"
    assert manifest["options"]["sequence_engine"] == "iqtree"
    assert Path(str(prefix) + ".pdf").read_bytes().startswith(b"%PDF")
    assert "engine=nwkit-iqtree-sequence-v1" in (tmp_path / "recorded-provenance.txt").read_text()
    import pandas as pd

    species = pd.read_csv(str(prefix) + ".species.tsv", sep="\t")
    root = species[species.node == "Root"].iloc[0]
    assert root.age == root.age_min == root.age_max == 10


@pytest.mark.parametrize("iqtree_model", ["", "GY+F3X4+R4"])
def test_external_iqtree_library_stage_preserves_outputs_and_provenance(tmp_path, iqtree_model):
    from nwkit.iqtree_library import find_worker

    if find_worker() is None:
        pytest.skip("External IQ-TREE library worker is not installed")
    result = _fixture(tmp_path, "gy94", "generax", engine="iqtree", iqtree_model=iqtree_model, use_defaults=not iqtree_model, iqtree_interface="library")
    assert result.returncode == 0, result.stdout + result.stderr
    prefix = tmp_path / "out/dated_tree_native/OG0000001_radte"
    manifest = json.loads(Path(str(prefix) + ".manifest.json").read_text())
    assert manifest["sequence_model"]["iqtree_interface"] == "library-worker-v1"
    library = manifest["sequence_model"]["iqtree_library"]
    assert len(library["library_sha256"]) == len(library["worker_sha256"]) == 64
    _assert_profile_output_contract(tmp_path / "out/dated_tree_native/OG0000001_radte", manifest)
    assert Path(str(prefix) + ".pdf").read_bytes().startswith(b"%PDF")
    provenance = (tmp_path / "recorded-provenance.txt").read_text()
    assert "iqtree_interface=library" in provenance
    assert library["library_sha256"] in provenance
    import pandas as pd

    species = pd.read_csv(str(prefix) + ".species.tsv", sep="\t")
    root = species[species.node == "Root"].iloc[0]
    assert root.age == root.age_min == root.age_max == 10
