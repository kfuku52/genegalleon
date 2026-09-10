import importlib.util
import json
import shlex
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
from detect_ou_shift_native import aggregate_replicates, main  # noqa: E402

TREE = "(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);"


def test_replicate_means_error_of_mean_and_missing_coordinates():
    table = pd.DataFrame(
        {
            "gene": ["a", "b", "c"],
            "root_1": [2.0, 3.0, np.nan],
            "root_2": [4.0, np.nan, np.nan],
            "leaf_1": [5.0, 6.0, 7.0],
        }
    )
    output, traits, errors, audit = aggregate_replicates(table)
    assert traits == ["root", "leaf"]
    assert output.root[:2].tolist() == [3, 3]
    assert output.loc[0, errors[0]] ** 2 == pytest.approx(1)
    assert np.isnan(output.loc[2, "root"])
    assert audit[1]["sampling_variance_status"] == "unavailable_single_observation"
    assert audit[2]["sampling_variance_status"] == "missing"
    assert aggregate_replicates(table, "")[1] == ["root_1", "root_2", "leaf_1"]
    table["root_1"] = table["root_1"].astype(object)
    table.loc[0, "root_1"] = "invalid"
    with pytest.raises(ValueError):
        aggregate_replicates(table)


def inputs(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text(TREE)
    table = pd.DataFrame(
        np.random.default_rng(39).normal(size=(8, 4)), columns=["root_1", "root_2", "leaf_1", "leaf_2"]
    )
    table.insert(0, "gene", list("abcdefgh"))
    table.loc[1, ["leaf_1", "leaf_2"]] = np.nan
    table.to_csv(tmp_path / "traits.tsv", sep="\t", index=False, na_rep="NA")
    return [
        "--tree",
        str(tree),
        "--traits",
        str(tmp_path / "traits.tsv"),
        "--output-prefix",
        str(tmp_path / "result"),
        "--max-shifts",
        "1",
        "--calibration-replicates",
        "19",
        "--alpha",
        "0.5",
        "--process-tip-variance",
        "1",
        "--estimate-measurement-error",
        "no",
    ]


def test_adapter_real_native_search_artifacts_resume_and_topology_safe_summary(tmp_path):
    args = inputs(tmp_path)
    main(args)
    model = json.loads((tmp_path / "result.model.json").read_text())
    assert model["information_criterion"]["criterion"] == "AICc"
    assert model["selection_calibration"] is None
    assert model["search"]["strategy"] == "covariance_updated_optimum_path"
    assert model["trait_names"] == ["root", "leaf"]
    assert len(model["genegalleon_adapter"]["replicate_audit"]) == 16
    assert (tmp_path / "result.pdf").read_bytes().startswith(b"%PDF")
    assert len(pd.read_csv(tmp_path / "result.tips.tsv", sep="\t")) == 16
    summary = pd.read_csv(tmp_path / "result.branch-summary.tsv", sep="\t")
    assert len(summary) == 15
    assert summary.ou_native_selection_frequency.isna().all()
    args[args.index("--output-prefix") + 1] = str(tmp_path / "resumed")
    args.extend(["--resume-model", str(tmp_path / "result.model.json")])
    main(args)
    assert json.loads((tmp_path / "resumed.model.json").read_text()) == model

    spec = importlib.util.spec_from_file_location("native_ou_stats", SUPPORT / "orthogroup_statistics.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    mapped = module.load_native_ou_branch_table(tmp_path / "result.model.json", tmp_path / "tree.nwk")
    root = mapped.loc[mapped.ou_native_regime == model["branches"][0]["regime"]]
    assert 14 in root.branch_id.values  # GeneGalleon root is the largest clade rank.
    assert set(mapped.branch_id) == set(range(15))
    (tmp_path / "tree.nwk").write_text(TREE.replace("a:1", "z:1"))
    with pytest.raises(ValueError, match="does not match"):
        module.load_native_ou_branch_table(tmp_path / "result.model.json", tmp_path / "tree.nwk")


def test_failed_adapter_does_not_replace_any_previous_artifact(tmp_path):
    args = inputs(tmp_path)
    destination = tmp_path / "result.model.json"
    destination.write_text("previous model")
    (tmp_path / "tree.nwk").write_text(TREE.replace("a:1", "a:0"))
    with pytest.raises(ValueError, match="positive"):
        main(args)
    assert destination.read_text() == "previous model"
    assert not (tmp_path / "result.pdf").exists()


def test_native_core_stage_executes_real_adapter_and_publishes_complete_bundle(tmp_path):
    inputs(tmp_path)
    core = (SUPPORT.parent / "core" / "gg_gene_evolution_core.sh").read_text()
    defaults = core.split('run_native_ou="${run_native_ou:-0}"', 1)[1].split("treevis_query_marker=", 1)[0]
    stage = core.split('task="NWKIT OU shift detection"', 1)[1].split(
        'task="Expression-trait phylogenetic regression"', 1
    )[0]
    variables = {
        "run_native_ou": "1",
        "native_ou_max_shifts": "0",
        "native_ou_convergence": "0",
        "native_ou_estimate_measurement_error": "no",
        "gg_support_dir": str(SUPPORT),
        "file_og_expression": str(tmp_path / "traits.tsv"),
        "file_og_dated_tree_analysis": str(tmp_path / "tree.nwk"),
        "file_og_native_ou_prefix": str(tmp_path / "published" / "family"),
        "dir_output_active": str(tmp_path / "published"),
        "gg_workspace_dir": str(tmp_path),
        "og_id": "family",
    }
    script = "set -euo pipefail\n" + "\n".join(f"{name}={shlex.quote(value)}" for name, value in variables.items())
    script += "\n" + defaults + "\n"
    script += """
disable_if_no_input_file() { :; }
gg_artifact_prepare_stage() { printf -v "$1" '%s' 1; }
gg_step_start() { :; }
gg_step_skip() { :; }
gg_artifact_record() { printf '%s\\n' "$@" > manifest-arguments.txt; }
mv_out_bundle() {
  while (( $# )); do mkdir -p "$(dirname "$2")"; mv "$1" "$2"; shift 2; done
}
task="native test"
"""
    script += stage
    subprocess.run(["bash", "-c", script], cwd=tmp_path, check=True, capture_output=True, text=True)
    model = json.loads((tmp_path / "published" / "family.model.json").read_text())
    assert model["shift_branch_ids"] == []
    assert len(list((tmp_path / "published").glob("family.*"))) == 8
    manifest = (tmp_path / "manifest-arguments.txt").read_text()
    assert "nwkit_implementation=" + model["implementation_sha256"] in manifest
    assert "replicate_separator=_" in manifest
    assert "criterion=AICc" in manifest
    assert "search_strategy=native-path" in manifest


def test_adapter_preserves_numeric_internal_names_for_summary_alignment(tmp_path):
    args = inputs(tmp_path)
    (tmp_path / "tree.nwk").write_text(TREE.replace("):1", ")95:1"))
    args[args.index("--max-shifts") + 1] = "0"
    main(args)
    model = json.loads((tmp_path / "result.model.json").read_text())
    assert sum(row["name"] == "95" for row in model["branches"]) == 6
    spec = importlib.util.spec_from_file_location("native_ou_numeric_stats", SUPPORT / "orthogroup_statistics.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert len(module.load_native_ou_branch_table(tmp_path / "result.model.json", tmp_path / "tree.nwk")) == 15


@pytest.mark.parametrize("criterion,strategy,convergence", [("bootstrap", "auto", True), ("BIC", "lasso", False)])
def test_adapter_alternative_selection_modes(tmp_path, criterion, strategy, convergence):
    args = inputs(tmp_path)
    args.extend(["--criterion", criterion, "--search-strategy", strategy])
    if convergence:
        args.append("--convergence")
    main(args)
    model = json.loads((tmp_path / "result.model.json").read_text())
    assert model["completion_status"] == "complete"
    if criterion == "bootstrap":
        assert model["selection_calibration"] is not None
    else:
        assert model["information_criterion"]["criterion"] == criterion


def test_empty_traits_fail_without_publishing_model(tmp_path):
    args = inputs(tmp_path)
    (tmp_path / "traits.tsv").write_text("gene\na\nb\n")
    with pytest.raises(ValueError, match="at least one trait"):
        main(args)
    assert not (tmp_path / "result.model.json").exists()
