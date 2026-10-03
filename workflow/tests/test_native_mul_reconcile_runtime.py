import json
import os
import stat
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
from nwkit.util import read_tree

ROOT = Path(__file__).resolve().parents[2]
CORE = ROOT / "workflow/core/gg_genome_evolution_core.sh"


def run_wrapper(tmp_path, *, h1="Hybrid_three", fault=None, gene_text=None, empty=False, locus=False, bootstrap=0, null_calibration="plug-in", locus_model_path=None, locus_species_text=None, missing_locus_species=False):
    source = CORE.read_text()
    start = source.index("busco_grampa() {")
    end = source.index("\n}\n", start) + 3
    trees = tmp_path / "rooted"
    trees.mkdir(exist_ok=True)
    if not empty:
        (trees / "OG001.nwk").write_text(gene_text or "((Alpha_one_a,Hybrid_three_x1),(Beta_two_b,Hybrid_three_x2));\n")
    species = tmp_path / "species.nwk"
    species.write_text("[&R]((Alpha_one:1,Hybrid_three:1):1,Beta_two:2);\n" if locus else "[&R]((Alpha_one,Hybrid_three),Beta_two);\n")
    model_path = tmp_path / "locus_model.json"
    locus_species = tmp_path / "generations.nwk"
    if locus:
        locus_species.write_text(locus_species_text or species.read_text())
    if locus:
        model_path.write_text(json.dumps({
            "schema": "nwkit-mul-locus-mc-model-v1", "copy_role": "distinct-loci",
            "root_locus_count": 1, "species_time_unit": "generations", "ancestral_stem": 0.5,
            "detection": {"Alpha_one": 0.8, "Hybrid_three": 0.8, "Beta_two": 0.8},
            "max_observed_tips": 4,
            "parameter_grid": [{"duplication": 0, "loss": 0, "ne": 1, "hybridization_age": 0.5}],
            "samples": 64, "seed": 20261115, "confidence": 0.99,
            "max_attempts": 10000, "max_locus_nodes": 1000, "max_coalescent_states": 100000,
            "integration": "hybrid-rb",
        }))
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    forbidden = bin_dir / "grampa.py"
    forbidden.write_text("#!/bin/sh\nexit 91\n")
    forbidden.chmod(forbidden.stat().st_mode | stat.S_IXUSR)
    if fault == "native":
        forbidden = bin_dir / "nwkit"
        forbidden.write_text("#!/bin/sh\nexit 92\n")
        forbidden.chmod(forbidden.stat().st_mode | stat.S_IXUSR)
    elif fault == "late-publish":
        forbidden = bin_dir / "mv"
        forbidden.write_text(
            '#!/bin/bash\nstage=0\nfor arg in "$@"; do\n'
            'case "$arg" in */.locus.gg-stage.*) stage=1;; esac\ndone\n'
            'if [[ $stage == 1 && "${@: -1}" == */results/locus ]]; then exit 94; fi\n'
            'exec /bin/mv "$@"\n'
        )
        forbidden.chmod(forbidden.stat().st_mode | stat.S_IXUSR)
    elif fault in ("publish", "locus-publish"):
        target = "locus" if fault == "locus-publish" else "grampa_out.txt"
        forbidden = bin_dir / "mv"
        forbidden.write_text(
            '#!/bin/bash\nfor arg in "$@"; do\n'
            f'case "$arg" in */results/{target}) exit 93;; esac\ndone\nexec /bin/mv "$@"\n'
        )
        forbidden.chmod(forbidden.stat().st_mode | stat.S_IXUSR)
    script = "\n".join(
        [
            "set -euo pipefail",
            'source "${gg_support_dir}/gg_util.sh"',
            source[start:end],
            'busco_grampa "$INDIR" "$OUTDIR" "$OUTDIR/grampa_summary.tsv"',
        ]
    )
    return subprocess.run(
        ["bash", "-c", script],
        cwd=tmp_path,
        env={
            **os.environ,
            "PATH": str(bin_dir) + os.pathsep + os.environ["PATH"],
            "INDIR": "rooted",
            "OUTDIR": "results",
            "file_dated_species_tree": "species.nwk",
            "grampa_h1": h1,
            "gg_support_dir": str(ROOT / "workflow/support"),
            "GG_TASK_CPUS": "1",
            "species_label_parser": "legacy",
            "species_label_regex": "",
            "species_label_map_tsv": "",
            "grampa_locus_model": str(locus_model_path or model_path) if locus else "",
            "grampa_locus_species_tree": str(locus_species) if locus and not missing_locus_species else "",
            "grampa_locus_h2": "Alpha_one Beta_two" if locus else "",
            "grampa_locus_bootstrap": str(bootstrap),
            "grampa_locus_null_calibration": null_calibration,
        },
        capture_output=True,
        text=True,
        timeout=120,
    )


@pytest.mark.parametrize("bootstrap,calibration", [(0, "plug-in"), (2, "grid-supremum")])
def test_optional_locus_pipeline_has_distinct_schema_and_normalized_model(tmp_path, bootstrap, calibration):
    result = run_wrapper(tmp_path, locus=True, bootstrap=bootstrap, null_calibration=calibration,
                         gene_text="((Alpha_one_a,Hybrid_three_x),Beta_two_b);\n")
    assert result.returncode == 0, result.stdout + result.stderr
    output = tmp_path / "results"
    old = json.loads((output / "nwkit_mul_reconcile.json").read_text())
    saved = json.loads((output / "locus/results.json").read_text())
    assert old["method"] == "exact-MUL-LCA-DL-parsimony-v1"
    assert saved["schema"] == "nwkit-mul-locus-mc-v1"
    assert saved["model"]["integration"] == "hybrid-rb"
    assert set(saved["model"]["detection"]) == {"Alpha-one", "Hybrid-three", "Beta-two"}
    assert set(json.loads((tmp_path / "locus_model.json").read_text())["detection"]) == {"Alpha_one", "Hybrid_three", "Beta_two"}
    assert saved["num_gene_trees"] == 1
    assert saved["banks"][0]["interval_method"] == "chernoff-kl"
    assert (output / "locus/null_search.tsv").exists() == bool(bootstrap)
    assert not (output / "locus/best_mul_tree.nwk").exists()
    assert not list(tmp_path.glob("tmp.mul-reconcile.*"))


@pytest.mark.parametrize("fault", ["invalid-family", "locus-publish", "late-publish", "invalid-calibration", "input-alias", "missing-generations", "wrong-topology", "invalid-bootstrap"])
def test_locus_failure_preserves_both_bundles(tmp_path, fault):
    gene = "((Alpha_one_a,Hybrid_three_x),Beta_two_b);\n"
    initial = run_wrapper(tmp_path, locus=True, gene_text=gene)
    assert initial.returncode == 0, initial.stdout + initial.stderr
    output = tmp_path / "results"
    before = {str(p.relative_to(output)): p.read_bytes() for p in output.rglob("*") if p.is_file()}
    result = run_wrapper(
        tmp_path, locus=True,
        gene_text="((Alpha_one_a,Hybrid_three_x1),(Beta_two_b,Hybrid_three_x2));" if fault == "invalid-family" else gene,
        fault=fault if fault in ("locus-publish", "late-publish") else None,
        null_calibration="grid-supremum" if fault == "invalid-calibration" else "plug-in",
        locus_model_path=output / "locus/input_model.json" if fault == "input-alias" else None,
        missing_locus_species=fault == "missing-generations",
        locus_species_text="((Alpha_one:1,Beta_two:1):1,Hybrid_three:2);" if fault == "wrong-topology" else None,
        bootstrap=-1 if fault == "invalid-bootstrap" else 0,
    )
    assert result.returncode != 0
    assert {str(p.relative_to(output)): p.read_bytes() for p in output.rglob("*") if p.is_file()} == before
    assert not list(tmp_path.glob("tmp.mul-reconcile.*"))


def test_generations_tree_child_order_preserves_numeric_selectors(tmp_path):
    result = run_wrapper(tmp_path, locus=True, gene_text="((Alpha_one_a,Hybrid_three_x),Beta_two_b);",
                         locus_species_text="(Beta_two:2,(Hybrid_three:1,Alpha_one:1):1);")
    assert result.returncode == 0, result.stdout + result.stderr
    species = read_tree(str(tmp_path / "results/locus/species_tree.nwk"), "auto", True, quiet=True)
    assert list(species.leaf_names()) == ["Alpha-one", "Hybrid-three", "Beta-two"]
    assert species.children[0].dist == 1


def test_optional_locus_provenance_covers_every_stage_and_config_content():
    source = CORE.read_text()
    assert source.count("grampa_locus_provenance ") == 3
    assert '--input "locus_config=${grampa_locus_model}"' in source
    assert '--optional-output "locus_bundle=${output_root}/locus"' in source


def test_native_replacement_runs_existing_busco_wrapper_and_summary_without_grampa(tmp_path):
    result = run_wrapper(tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    output = tmp_path / "results"
    summary = pd.read_csv(output / "grampa_summary.tsv", sep="\t")
    assert set(summary["file_name"]) == {"OG001.nwk"}
    assert set(summary["total_score"]) == {0}
    assert set(summary["Hybrid-three"]) == {"x1,x2"}
    assert {"grampa_det.txt", "grampa_out.txt", "grampa_checknums.txt", "best_mul_tree.nwk"}.issubset(
        {path.name for path in output.iterdir()}
    )
    metadata = json.loads((output / "nwkit_mul_reconcile.json").read_text())
    assert metadata["method"] == "exact-MUL-LCA-DL-parsimony-v1"
    assert metadata["gene_filtering"].startswith("none")
    assert (output / "best_mul_tree.nwk").read_text().strip().endswith(";")
    assert not list(tmp_path.glob("tmp.mul-reconcile.*"))


def test_multiple_h1_selectors_and_existing_scratch_files(tmp_path):
    scratch = tmp_path / "grampa_out"
    scratch.mkdir()
    sentinel = scratch / "curated.txt"
    sentinel.write_text("unrelated existing content")
    trees = tmp_path / "rooted"
    trees.mkdir()
    (trees / ".ignored.nwk").write_text("an unrelated hidden file, not a gene tree")
    result = run_wrapper(tmp_path, h1="Hybrid_three\t1")
    assert result.returncode == 0, result.stdout + result.stderr
    scores = pd.read_csv(tmp_path / "results/grampa_out.txt", sep="\t", keep_default_na=False)
    assert set(scores["h1.node"]) == {"NA", "Hybrid-three", "<1>"}
    assert sentinel.read_text() == "unrelated existing content"
    assert (tmp_path / "results/busco_genetree_filenames.txt").read_text() == "OG001.nwk\n"


@pytest.mark.parametrize("fault", ["native", "publish", "invalid"])
def test_failed_analysis_or_bundle_publish_preserves_every_prior_result(tmp_path, fault):
    initial = run_wrapper(tmp_path)
    assert initial.returncode == 0, initial.stdout + initial.stderr
    output = tmp_path / "results"
    before = {path.name: path.read_bytes() for path in output.iterdir() if path.is_file()}
    result = run_wrapper(
        tmp_path,
        h1="1",
        fault=fault,
        gene_text="(Alpha_one_a,Unknown_species_b);" if fault == "invalid" else None,
    )
    assert result.returncode != 0
    assert {path.name: path.read_bytes() for path in output.iterdir() if path.is_file()} == before
    assert not list(tmp_path.glob("tmp.mul-reconcile.*"))


def test_empty_input_does_not_relabel_previous_results_as_current(tmp_path):
    result = run_wrapper(tmp_path, empty=True)
    assert result.returncode == 0, result.stderr
    initial = run_wrapper(tmp_path)
    assert initial.returncode == 0, initial.stderr
    summary = tmp_path / "results/grampa_summary.tsv"
    before = summary.read_bytes()
    (tmp_path / "rooted/OG001.nwk").unlink()
    result = run_wrapper(tmp_path, empty=True)
    assert result.returncode != 0
    assert summary.read_bytes() == before
    summary.unlink()
    result = run_wrapper(tmp_path, empty=True)
    assert result.returncode != 0
    assert (tmp_path / "results/nwkit_mul_reconcile.json").is_file()


def test_all_three_stage_contracts_invalidate_old_grampa_engine():
    source = CORE.read_text()
    assert source.count('--parameter "engine=nwkit-mul-reconcile-v1"') == 3
    assert 'grampa.py "${grampa_args[@]}"' not in source
    for prefix in ("busco_grampa_dna", "busco_grampa_pep", "orthogroup_grampa"):
        assert f'{prefix}_provenance_args+=(--input "species_map=${{species_label_map_tsv}}")' in source
    assert source.count('--optional-output "model=') == 3
    assert source.count('--optional-output "best_tree=') == 3
    assert 'rm -f -- "${file_busco_grampa' not in source
    assert 'rm -f -- "${file_orthogroup_grampa}"' not in source


@pytest.mark.parametrize("preset,use_map", [("legacy", False), ("taxonomic", False), ("legacy", True)])
def test_structured_preparation_preserves_topology_and_species_mapping(tmp_path, preset, use_map):
    from nwkit.mul_reconcile_model import validate_binary

    alpha = "Alpha_one_subsp_blue" if preset == "taxonomic" else "Alpha_one"
    species = tmp_path / "species.nwk"
    species.write_text(f"[&R](({alpha},Hybrid_three),Beta_two);\n")
    genes = tmp_path / "OG001.nwk"
    labels = [alpha + "_a", "Hybrid_three_x1", "Beta_two_b", "Hybrid_three_x2"]
    if use_map:
        labels = ["raw1", "raw2", "raw3", "raw4"]
    genes.write_text(f"[&R](({labels[0]}:0.1,{labels[1]}:0.2)99:0.3,({labels[2]}:0.4,{labels[3]}:0.5)88:0.6);\n")
    paths = {
        "species-out": tmp_path / "prepared_species.nwk",
        "genes-out": tmp_path / "prepared_genes.nwk",
        "names-out": tmp_path / "names.txt",
    }
    args = [
        sys.executable,
        str(ROOT / "workflow/support/prepare_mul_reconcile.py"),
        "--species-tree",
        str(species),
        "--gene-tree",
        str(genes),
        "--species-parser",
        preset,
    ]
    if use_map:
        mapping = tmp_path / "species.tsv"
        mapping.write_text(
            "leaf_name\tspecies_label\nraw1\tAlpha_one\nraw2\tHybrid_three\nraw3\tBeta_two\nraw4\tHybrid_three\n"
        )
        args += ["--species-map-tsv", str(mapping)]
    for flag, path in paths.items():
        args += ["--" + flag, str(path)]
    result = subprocess.run(args, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    prepared = read_tree(str(paths["genes-out"]), "auto", True, quiet=True)
    validate_binary(prepared, "Prepared gene tree")
    assert prepared.children[0].dist == pytest.approx(0.3)
    assert prepared.children[0].children[0].dist == pytest.approx(0.1)
    first = "raw1" if use_map else "a"
    assert list(prepared.children[0].leaf_names()) == [
        first + "_" + alpha.replace("_", "-"),
        ("raw2" if use_map else "x1") + "_Hybrid-three",
    ]
    assert paths["names-out"].read_text() == "OG001.nwk\n"


@pytest.mark.parametrize(
    "bad_input",
    [
        "[&U](Alpha_one_a,Beta_two_b);",
        "(Alpha_one_a,Unknown_two_b);",
        "(Alpha_one_a,Alpha_one_a);",
        "(Alpha_one_a,Beta_two_b,Hybrid_three_x);",
    ],
)
def test_invalid_preparation_preserves_existing_files(tmp_path, bad_input):
    species = tmp_path / "species.nwk"
    species.write_text("((Alpha_one,Hybrid_three),Beta_two);\n")
    genes = tmp_path / "gene.nwk"
    genes.write_text(bad_input)
    args = [
        sys.executable,
        str(ROOT / "workflow/support/prepare_mul_reconcile.py"),
        "--species-tree",
        str(species),
        "--gene-tree",
        str(genes),
    ]
    outputs = [tmp_path / name for name in ("prepared_species.nwk", "prepared_genes.nwk", "names.txt")]
    for flag, path in zip(("species-out", "genes-out", "names-out"), outputs, strict=True):
        path.write_text("prior output")
        args += ["--" + flag, str(path)]
    result = subprocess.run(args, capture_output=True, text=True, timeout=30)
    assert result.returncode != 0
    assert all(path.read_text() == "prior output" for path in outputs)


def test_duplicate_preparation_targets_fail_without_losing_a_file(tmp_path):
    species = tmp_path / "species.nwk"
    species.write_text("(Alpha_one,Beta_two);")
    genes = tmp_path / "gene.nwk"
    genes.write_text("(Alpha_one_a,Beta_two_b);")
    shared = tmp_path / "shared.nwk"
    shared.write_text("prior output")
    result = subprocess.run(
        [
            sys.executable, str(ROOT / "workflow/support/prepare_mul_reconcile.py"),
            "--species-tree", str(species), "--gene-tree", str(genes),
            "--species-out", str(shared), "--genes-out", str(shared),
            "--names-out", str(tmp_path / "names.txt"),
        ], capture_output=True, text=True, timeout=30,
    )
    assert result.returncode != 0
    assert shared.read_text() == "prior output"
    assert not (tmp_path / "names.txt").exists()


def test_filename_inventory_is_literal_not_csv_quoted_or_missing(tmp_path):
    initial = run_wrapper(tmp_path)
    assert initial.returncode == 0, initial.stderr
    original = tmp_path / "rooted/OG001.nwk"
    original.rename(tmp_path / 'rooted/"quoted".nwk')
    result = run_wrapper(tmp_path, empty=True)
    assert result.returncode == 0, result.stderr
    summary = pd.read_csv(tmp_path / "results/grampa_summary.tsv", sep="\t", keep_default_na=False)
    assert set(summary["file_name"]) == {'"quoted".nwk'}
