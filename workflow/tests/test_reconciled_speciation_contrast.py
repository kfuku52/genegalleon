from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "support" / "reconciled_speciation_contrast.py"


def load_module():
    spec = importlib.util.spec_from_file_location("reconciled_speciation_contrast", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    pandas.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_prepare_converts_wide_expression_replicates_without_false_pairing(tmp_path: Path):
    mod = load_module()
    expression = tmp_path / "expression.tsv"
    traits = tmp_path / "traits.tsv"
    write_tsv(
        expression,
        [
            {"gene_id": "Genus_a_g1", "root_1": 1.0, "root_2": 1.4, "leaf_1": 5.0, "leaf_2": 5.5},
            {"gene_id": "Genus_b_g1", "root_1": 2.0, "root_2": 2.4, "leaf_1": 6.0, "leaf_2": 6.5},
        ],
    )
    write_tsv(
        traits,
        [
            {"species": "Genus_a", "body_size": 1.0, "habitat": "wet"},
            {"species": "Genus_b", "body_size": 2.0, "habitat": "dry"},
        ],
    )
    expression_out = tmp_path / "expression.out.tsv"
    traits_out = tmp_path / "traits.out.tsv"
    plan = tmp_path / "plan.tsv"
    metadata = tmp_path / "metadata.tsv"

    assert (
        mod.main(
            [
                "prepare",
                "--expression",
                str(expression),
                "--species-traits",
                str(traits),
                "--responses",
                "all",
                "--predictors",
                "all",
                "--predictor-mode",
                "separate",
                "--expression-output",
                str(expression_out),
                "--species-traits-output",
                str(traits_out),
                "--analysis-plan-output",
                str(plan),
                "--metadata-output",
                str(metadata),
            ]
        )
        == 0
    )

    observed = pandas.read_csv(expression_out, sep="\t", keep_default_na=False)
    assert list(observed.columns) == ["leaf_name", "root", "leaf", "biological_id"]
    assert set(observed["biological_id"]) == {
        "auto:root:1",
        "auto:root:2",
        "auto:leaf:1",
        "auto:leaf:2",
    }
    assert observed.groupby(["leaf_name", "biological_id"]).size().max() == 1
    assert observed.loc[observed["biological_id"].str.startswith("auto:root"), "leaf"].eq("NA").all()

    trait_observed = pandas.read_csv(traits_out, sep="\t", keep_default_na=False)
    assert list(trait_observed.columns) == ["leaf_name", "body_size", "habitat"]
    plan_observed = pandas.read_csv(plan, sep="\t", keep_default_na=False)
    assert list(plan_observed["predictors"]) == ["body_size", "habitat"]
    assert plan_observed.loc[plan_observed["predictors"] == "habitat", "categorical_predictors"].iloc[0] == "habitat"

    metadata_observed = dict(
        pandas.read_csv(metadata, sep="\t", keep_default_na=False).itertuples(index=False, name=None)
    )
    assert metadata_observed["status"] == "ready"
    assert metadata_observed["responses"] == "root,leaf"
    assert metadata_observed["response_biological_id"] == "biological_id"


def test_prepare_supports_explicit_technical_replicates_and_batch(tmp_path: Path):
    mod = load_module()
    expression = tmp_path / "expression.tsv"
    traits = tmp_path / "traits.tsv"
    sample_metadata = tmp_path / "samples.tsv"
    write_tsv(
        expression,
        [
            {"gene": "Genus_a_g1", "s1": 1.0, "s2": 1.2},
            {"gene": "Genus_b_g1", "s1": 2.0, "s2": 2.2},
        ],
    )
    write_tsv(
        sample_metadata,
        [
            {
                "column": "s1",
                "response": "root",
                "biological_id": "bio1",
                "technical_id": "tech1",
                "batch": "A",
            },
            {
                "column": "s2",
                "response": "root",
                "biological_id": "bio1",
                "technical_id": "tech2",
                "batch": "A",
            },
        ],
    )
    write_tsv(
        traits,
        [
            {"species": "Genus_a", "size": 1.0},
            {"species": "Genus_b", "size": 2.0},
        ],
    )

    argv = [
        "prepare",
        "--expression",
        str(expression),
        "--species-traits",
        str(traits),
        "--sample-metadata",
        str(sample_metadata),
        "--expression-output",
        str(tmp_path / "expression.out.tsv"),
        "--species-traits-output",
        str(tmp_path / "traits.out.tsv"),
        "--analysis-plan-output",
        str(tmp_path / "plan.tsv"),
        "--metadata-output",
        str(tmp_path / "metadata.tsv"),
    ]
    assert mod.main(argv) == 0
    observed = pandas.read_csv(tmp_path / "expression.out.tsv", sep="\t")
    assert {"biological_id", "technical_id", "batch"}.issubset(observed.columns)
    assert set(observed["technical_id"]) == {"tech1", "tech2"}
    metadata = dict(
        pandas.read_csv(tmp_path / "metadata.tsv", sep="\t", keep_default_na=False).itertuples(index=False, name=None)
    )
    assert metadata["response_technical_id"] == "technical_id"
    assert metadata["response_batch"] == "batch"
    assert metadata["standard_error_columns"] == ""


def test_prepare_supports_summarized_expression_known_se(tmp_path: Path):
    mod = load_module()
    write_tsv(
        tmp_path / "expression.tsv",
        [
            {"gene": "Genus_a_g1", "root_mean": 1.0, "root_se": 0.1, "root_n": 4},
            {"gene": "Genus_b_g1", "root_mean": 2.0, "root_se": 0.2, "root_n": 5},
        ],
    )
    write_tsv(
        tmp_path / "samples.tsv",
        [
            {
                "column": "root_mean",
                "response": "root",
                "standard_error_column": "root_se",
                "sample_size_column": "root_n",
            }
        ],
    )
    write_tsv(
        tmp_path / "traits.tsv",
        [{"species": "Genus_a", "size": 1.0}, {"species": "Genus_b", "size": 2.0}],
    )
    assert (
        mod.main(
            [
                "prepare",
                "--expression",
                str(tmp_path / "expression.tsv"),
                "--species-traits",
                str(tmp_path / "traits.tsv"),
                "--sample-metadata",
                str(tmp_path / "samples.tsv"),
                "--within-variance",
                "known-se",
                "--expression-output",
                str(tmp_path / "expression.out.tsv"),
                "--species-traits-output",
                str(tmp_path / "traits.out.tsv"),
                "--analysis-plan-output",
                str(tmp_path / "plan.tsv"),
                "--metadata-output",
                str(tmp_path / "metadata.tsv"),
            ]
        )
        == 0
    )
    observed = pandas.read_csv(tmp_path / "expression.out.tsv", sep="\t")
    assert list(observed.columns) == ["leaf_name", "root", "root__standard_error", "root__sample_size"]
    assert observed.shape[0] == 2
    metadata = dict(
        pandas.read_csv(tmp_path / "metadata.tsv", sep="\t", keep_default_na=False).itertuples(index=False, name=None)
    )
    assert metadata["response_biological_id"] == ""
    assert metadata["standard_error_columns"] == "root__standard_error"
    assert metadata["sample_size_columns"] == "root__sample_size"


def test_prepare_routes_predictor_replicates_and_known_se_per_analysis(tmp_path: Path):
    mod = load_module()
    write_tsv(
        tmp_path / "expression.tsv",
        [
            {"gene": "Genus_a_g1", "expression": 1.0},
            {"gene": "Genus_b_g1", "expression": 2.0},
        ],
    )
    write_tsv(
        tmp_path / "traits.tsv",
        [
            {"species": "Genus_a", "size": 1.0, "height": 3.0, "size_se": 0.1, "height_se": 0.3},
            {"species": "Genus_b", "size": 2.0, "height": 4.0, "size_se": 0.2, "height_se": 0.4},
        ],
    )
    assert (
        mod.main(
            [
                "prepare",
                "--expression",
                str(tmp_path / "expression.tsv"),
                "--species-traits",
                str(tmp_path / "traits.tsv"),
                "--predictors",
                "size,height",
                "--predictor-within-variance",
                "known-se",
                "--predictor-standard-error-columns",
                "size_se,height_se",
                "--expression-output",
                str(tmp_path / "expression.out.tsv"),
                "--species-traits-output",
                str(tmp_path / "traits.out.tsv"),
                "--analysis-plan-output",
                str(tmp_path / "plan.tsv"),
                "--metadata-output",
                str(tmp_path / "metadata.tsv"),
            ]
        )
        == 0
    )
    plan = pandas.read_csv(tmp_path / "plan.tsv", sep="\t", keep_default_na=False)
    assert list(plan["predictor_standard_error_columns"]) == ["size_se", "height_se"]
    assert list(plan["predictor_sample_size_columns"]) == [".", "."]

    duplicated_traits = pandas.read_csv(tmp_path / "traits.tsv", sep="\t")
    pandas.concat([duplicated_traits, duplicated_traits.iloc[[0]]], ignore_index=True).to_csv(
        tmp_path / "traits.duplicated.tsv", sep="\t", index=False
    )

    with pytest.raises(SystemExit):
        mod.main(
            [
                "prepare",
                "--expression",
                str(tmp_path / "expression.tsv"),
                "--species-traits",
                str(tmp_path / "traits.duplicated.tsv"),
                "--predictors",
                "size",
                "--expression-output",
                str(tmp_path / "bad-expression.tsv"),
                "--species-traits-output",
                str(tmp_path / "bad-traits.tsv"),
                "--analysis-plan-output",
                str(tmp_path / "bad-plan.tsv"),
                "--metadata-output",
                str(tmp_path / "bad-metadata.tsv"),
            ]
        )


def test_prepare_allows_partial_missing_raw_predictor_replicates(tmp_path: Path):
    mod = load_module()
    write_tsv(
        tmp_path / "expression.tsv",
        [
            {"gene": "Genus_a_g1", "expression": 1.0},
            {"gene": "Genus_b_g1", "expression": 2.0},
        ],
    )
    write_tsv(
        tmp_path / "traits.tsv",
        [
            {"species": "Genus_a", "sample": "a1", "size": 1.0},
            {"species": "Genus_a", "sample": "a2", "size": "NA"},
            {"species": "Genus_b", "sample": "b1", "size": 2.0},
            {"species": "Genus_b", "sample": "b2", "size": 2.4},
        ],
    )
    assert (
        mod.main(
            [
                "prepare",
                "--expression",
                str(tmp_path / "expression.tsv"),
                "--species-traits",
                str(tmp_path / "traits.tsv"),
                "--predictors",
                "size",
                "--predictor-biological-id",
                "sample",
                "--expression-output",
                str(tmp_path / "expression.out.tsv"),
                "--species-traits-output",
                str(tmp_path / "traits.out.tsv"),
                "--analysis-plan-output",
                str(tmp_path / "plan.tsv"),
                "--metadata-output",
                str(tmp_path / "metadata.tsv"),
            ]
        )
        == 0
    )
    metadata = dict(
        pandas.read_csv(tmp_path / "metadata.tsv", sep="\t", keep_default_na=False).itertuples(
            index=False, name=None
        )
    )
    assert metadata["status"] == "ready"
    assert metadata["predictors"] == "size"
    prepared = pandas.read_csv(tmp_path / "traits.out.tsv", sep="\t", keep_default_na=False)
    assert prepared.loc[prepared["sample"] == "a2", "size"].tolist() == ["NA"]


def test_prepare_reports_not_estimable_for_incomplete_or_constant_responses(tmp_path: Path):
    mod = load_module()
    write_tsv(
        tmp_path / "expression.tsv",
        [
            {"gene": "Genus_a_g1", "constant": 1, "missing": 1},
            {"gene": "Genus_b_g1", "constant": 1, "missing": "NA"},
        ],
    )
    write_tsv(
        tmp_path / "traits.tsv",
        [{"species": "Genus_a", "size": 1}, {"species": "Genus_b", "size": 2}],
    )
    mod.main(
        [
            "prepare",
            "--expression",
            str(tmp_path / "expression.tsv"),
            "--species-traits",
            str(tmp_path / "traits.tsv"),
            "--expression-output",
            str(tmp_path / "expression.out.tsv"),
            "--species-traits-output",
            str(tmp_path / "traits.out.tsv"),
            "--analysis-plan-output",
            str(tmp_path / "plan.tsv"),
            "--metadata-output",
            str(tmp_path / "metadata.tsv"),
        ]
    )
    metadata = dict(
        pandas.read_csv(tmp_path / "metadata.tsv", sep="\t", keep_default_na=False).itertuples(index=False, name=None)
    )
    assert metadata["status"] == "not_estimable"
    assert "constant:constant_or_empty" in metadata["reason"]
    assert "missing:missing_leaf_values" in metadata["reason"]


def test_aggregate_adds_analysis_ids_to_every_bundle_table_and_summarizes(tmp_path: Path):
    mod = load_module()
    bundle_rows = []
    for analysis_id, predictor, coefficient, p_value in (
        ("p001_size", "size", 2.0, 0.04),
        ("p002_habitat", "habitat[dry]", -3.0, 0.01),
    ):
        prefix = tmp_path / analysis_id
        shared = pandas.DataFrame([{"tree_id": "OG1", "gene_clade_id": "a|b", "event_type": "speciation"}])
        for suffix in mod.REQUIRED_NWKIT_SUFFIXES:
            if suffix == ".regression.tsv":
                frame = pandas.DataFrame(
                    [
                        {
                            "model_id": f"OG1:{analysis_id}",
                            "tree_id": "OG1",
                            "response": "root",
                            "term": predictor,
                            "coefficient": coefficient,
                            "standard_error": 0.5,
                            "p_value": p_value,
                            "n_gene_contrasts": 5,
                            "n_species_events": 4,
                            "inference_status": "ok",
                        }
                    ]
                )
            elif suffix == ".random-effects.tsv":
                frame = pandas.DataFrame(columns=["model_id", "tree_id", "response", "term"])
            else:
                frame = shared
            frame.to_csv(f"{prefix}{suffix}", sep="\t", index=False)
        audit = tmp_path / f"{analysis_id}.audit.jsonl"
        audit.write_text(json.dumps({"command": "regress"}) + "\n", encoding="utf-8")
        bundle_rows.append({"analysis_id": analysis_id, "prefix": prefix, "audit": audit})
    write_tsv(tmp_path / "bundles.tsv", bundle_rows)

    assert (
        mod.main(
            [
                "aggregate",
                "--bundle-list",
                str(tmp_path / "bundles.tsv"),
                "--output-prefix",
                str(tmp_path / "combined"),
                "--status-output",
                str(tmp_path / "status.tsv"),
                "--audit-output",
                str(tmp_path / "audit.jsonl"),
                "--tree-id",
                "OG1",
            ]
        )
        == 0
    )
    reconciliation = pandas.read_csv(tmp_path / "combined.reconciliation.tsv", sep="\t")
    assert reconciliation.shape[0] == 2
    assert set(reconciliation["analysis_id"]) == {"p001_size", "p002_habitat"}
    results = pandas.read_csv(tmp_path / "combined.regression.tsv", sep="\t")
    assert list(results["analysis_id"]) == ["p001_size", "p002_habitat"]
    status = pandas.read_csv(tmp_path / "status.tsv", sep="\t")
    assert status.loc[0, "status"] == "ok"
    assert status.loc[0, "min_p_value_raw"] == 0.01
    assert status.loc[0, "min_p_value"] == 0.02
    assert status.loc[0, "best_analysis_id"] == "p002_habitat"
    assert len((tmp_path / "audit.jsonl").read_text().splitlines()) == 2

    tree_stats = mod.summarize_for_stat_tree(tmp_path / "combined.regression.tsv", tmp_path / "status.tsv")
    assert tree_stats["rsc_status"] == "ok"
    assert tree_stats["rsc_min_p_value_raw"] == 0.01
    assert tree_stats["rsc_min_p_value"] == 0.02
    assert tree_stats["rsc_best_analysis_id"] == "p002_habitat"
    assert tree_stats["rsc_best_coefficient"] == -3.0
    assert len(tree_stats) < 100


def test_prepare_combines_paired_responses_and_keeps_single_response_wide(tmp_path: Path):
    mod = load_module()
    write_tsv(
        tmp_path / "expression.tsv",
        [
            {"gene": "A_g1", "r1": 1.0, "l1": 3.0},
            {"gene": "B_g1", "r1": 2.0, "l1": 4.0},
        ],
    )
    write_tsv(
        tmp_path / "samples.tsv",
        [
            {"column": "r1", "response": "root", "biological_id": "sample1"},
            {"column": "l1", "response": "leaf", "biological_id": "sample1"},
        ],
    )
    write_tsv(
        tmp_path / "traits.tsv",
        [{"species": "A", "size": 1}, {"species": "B", "size": 2}],
    )
    assert (
        mod.main(
            [
                "prepare",
                "--expression",
                str(tmp_path / "expression.tsv"),
                "--species-traits",
                str(tmp_path / "traits.tsv"),
                "--sample-metadata",
                str(tmp_path / "samples.tsv"),
                "--expression-output",
                str(tmp_path / "expression.out.tsv"),
                "--species-traits-output",
                str(tmp_path / "traits.out.tsv"),
                "--analysis-plan-output",
                str(tmp_path / "plan.tsv"),
                "--metadata-output",
                str(tmp_path / "metadata.tsv"),
            ]
        )
        == 0
    )
    observed = pandas.read_csv(tmp_path / "expression.out.tsv", sep="\t")
    assert list(observed.columns) == ["leaf_name", "root", "leaf"]
    assert observed.shape == (2, 3)
    metadata = dict(
        pandas.read_csv(tmp_path / "metadata.tsv", sep="\t", keep_default_na=False).itertuples(index=False, name=None)
    )
    assert metadata["response_biological_id"] == ""
    assert metadata["response_sampling_uncertainty"] == "no"


def test_prepare_treats_na_optional_sample_roles_as_absent(tmp_path: Path):
    mod = load_module()
    write_tsv(
        tmp_path / "expression.tsv",
        [{"gene": "A_g1", "s1": 1.0}, {"gene": "B_g1", "s1": 2.0}],
    )
    write_tsv(
        tmp_path / "samples.tsv",
        [
            {
                "column": "s1",
                "response": "expression",
                "biological_id": "bio1",
                "technical_id": "NA",
                "batch": "NA",
            }
        ],
    )
    write_tsv(
        tmp_path / "traits.tsv",
        [{"species": "A", "size": 1}, {"species": "B", "size": 2}],
    )
    assert (
        mod.main(
            [
                "prepare",
                "--expression",
                str(tmp_path / "expression.tsv"),
                "--species-traits",
                str(tmp_path / "traits.tsv"),
                "--sample-metadata",
                str(tmp_path / "samples.tsv"),
                "--expression-output",
                str(tmp_path / "expression.out.tsv"),
                "--species-traits-output",
                str(tmp_path / "traits.out.tsv"),
                "--analysis-plan-output",
                str(tmp_path / "plan.tsv"),
                "--metadata-output",
                str(tmp_path / "metadata.tsv"),
            ]
        )
        == 0
    )
    assert list(pandas.read_csv(tmp_path / "expression.out.tsv", sep="\t")) == [
        "leaf_name",
        "expression",
    ]


def test_prepare_combines_paired_responses_with_mixed_replication(tmp_path: Path):
    mod = load_module()
    write_tsv(
        tmp_path / "expression.tsv",
        [
            {"gene": "A_g1", "r1": 1.0, "l1": 3.0, "r2": 1.5},
            {"gene": "B_g1", "r1": 2.0, "l1": 4.0, "r2": 2.5},
        ],
    )
    write_tsv(
        tmp_path / "samples.tsv",
        [
            {"column": "r1", "response": "root", "biological_id": "sample1"},
            {"column": "l1", "response": "leaf", "biological_id": "sample1"},
            {"column": "r2", "response": "root", "biological_id": "sample2"},
        ],
    )
    write_tsv(
        tmp_path / "traits.tsv",
        [{"species": "A", "size": 1}, {"species": "B", "size": 2}],
    )
    assert (
        mod.main(
            [
                "prepare",
                "--expression",
                str(tmp_path / "expression.tsv"),
                "--species-traits",
                str(tmp_path / "traits.tsv"),
                "--sample-metadata",
                str(tmp_path / "samples.tsv"),
                "--expression-output",
                str(tmp_path / "expression.out.tsv"),
                "--species-traits-output",
                str(tmp_path / "traits.out.tsv"),
                "--analysis-plan-output",
                str(tmp_path / "plan.tsv"),
                "--metadata-output",
                str(tmp_path / "metadata.tsv"),
            ]
        )
        == 0
    )
    observed = pandas.read_csv(tmp_path / "expression.out.tsv", sep="\t", keep_default_na=False)
    assert observed.shape[0] == 4
    paired = observed[observed["biological_id"] == "sample1"]
    assert paired["root"].ne("NA").all()
    assert paired["leaf"].ne("NA").all()
    unpaired = observed[observed["biological_id"] == "sample2"]
    assert unpaired["leaf"].eq("NA").all()


def test_prepare_header_only_expression_is_not_estimable(tmp_path: Path):
    mod = load_module()
    (tmp_path / "expression.tsv").write_text("gene\texpression\n", encoding="utf-8")
    write_tsv(
        tmp_path / "traits.tsv",
        [{"species": "A", "size": 1}, {"species": "B", "size": 2}],
    )
    assert (
        mod.main(
            [
                "prepare",
                "--expression",
                str(tmp_path / "expression.tsv"),
                "--species-traits",
                str(tmp_path / "traits.tsv"),
                "--expression-output",
                str(tmp_path / "expression.out.tsv"),
                "--species-traits-output",
                str(tmp_path / "traits.out.tsv"),
                "--analysis-plan-output",
                str(tmp_path / "plan.tsv"),
                "--metadata-output",
                str(tmp_path / "metadata.tsv"),
            ]
        )
        == 0
    )
    metadata = dict(
        pandas.read_csv(tmp_path / "metadata.tsv", sep="\t", keep_default_na=False).itertuples(index=False, name=None)
    )
    assert metadata["status"] == "not_estimable"
    assert "expression_has_no_rows" in metadata["reason"]


def test_status_excludes_nonconverged_rows_from_best_result(tmp_path: Path):
    mod = load_module()
    results = pandas.DataFrame(
        [
            {
                "analysis_id": "failed",
                "model_id": "m1",
                "response": "expression",
                "term": "size",
                "p_value": 1e-12,
                "inference_status": "ok",
                "optimizer_converged": "no",
            },
            {
                "analysis_id": "usable",
                "model_id": "m2",
                "response": "expression",
                "term": "size",
                "p_value": 0.02,
                "inference_status": "ok",
                "optimizer_converged": "yes",
            },
        ]
    )
    status = mod._status_from_results("OG1", results, 2)
    assert status.loc[0, "best_analysis_id"] == "usable"
    assert status.loc[0, "min_p_value"] == 0.02
    assert status.loc[0, "n_estimable_models"] == 1
    assert status.loc[0, "n_nonconverged_result_rows"] == 1


def test_status_and_stat_tree_summary_exclude_intercepts(tmp_path: Path):
    mod = load_module()
    results = pandas.DataFrame(
        [
            {
                "analysis_id": "p001_size",
                "model_id": "m1",
                "response": "expression",
                "predictor_type": "intercept",
                "term": "Intercept",
                "coefficient": 20.0,
                "p_value": 1e-30,
                "inference_status": "ok",
                "optimizer_converged": "yes",
            },
            {
                "analysis_id": "p001_size",
                "model_id": "m1",
                "response": "expression",
                "predictor_type": "continuous",
                "term": "size",
                "coefficient": 1.5,
                "p_value": 0.03,
                "inference_status": "ok",
                "optimizer_converged": "yes",
            },
        ]
    )
    status = mod._status_from_results("OG1", results, 4)
    assert status.loc[0, "status"] == "ok"
    assert status.loc[0, "best_term"] == "size"
    assert status.loc[0, "min_p_value"] == 0.03

    results.to_csv(tmp_path / "pgls.tsv", sep="\t", index=False)
    status.to_csv(tmp_path / "status.tsv", sep="\t", index=False)
    summary = mod.summarize_for_stat_tree(tmp_path / "pgls.tsv", tmp_path / "status.tsv")
    assert summary["rsc_best_term"] == "size"
    assert summary["rsc_min_p_value"] == 0.03

    intercept_only = mod._status_from_results("OG1", results.iloc[[0]], 4)
    assert intercept_only.loc[0, "status"] == "not_estimable"


def test_status_and_stat_tree_summary_adjust_all_usable_associations(tmp_path: Path):
    mod = load_module()
    results = pandas.DataFrame(
        [
            {
                "analysis_id": f"p{index}",
                "model_id": f"m{index}",
                "response": "expression",
                "predictor_type": "continuous",
                "term": "size",
                "p_value": p_value,
                "inference_status": "ok",
                "optimizer_converged": "yes",
            }
            for index, p_value in enumerate((0.01, 0.03, 0.2), start=1)
        ]
    )
    status = mod._status_from_results("OG1", results, 3)
    assert status.loc[0, "n_tested_associations"] == 3
    assert status.loc[0, "min_p_value_raw"] == pytest.approx(0.01)
    assert status.loc[0, "min_p_value_holm"] == pytest.approx(0.03)
    assert status.loc[0, "min_p_value_bh"] == pytest.approx(0.03)
    assert status.loc[0, "min_p_value"] == pytest.approx(0.03)

    results.to_csv(tmp_path / "pgls.tsv", sep="\t", index=False)
    status.to_csv(tmp_path / "status.tsv", sep="\t", index=False)
    summary = mod.summarize_for_stat_tree(
        tmp_path / "pgls.tsv", tmp_path / "status.tsv"
    )
    assert summary["rsc_num_tested_associations"] == 3
    assert summary["rsc_best_p_value_raw"] == pytest.approx(0.01)
    assert summary["rsc_best_p_value"] == pytest.approx(0.03)
    assert summary["rsc_best_p_value_adjustment"] == "holm"


def test_audit_value_error_is_classified_as_analysis_not_estimable(tmp_path: Path):
    mod = load_module()
    audit = tmp_path / "audit.jsonl"
    audit.write_text(
        json.dumps(
            {
                "status": "error",
                "error": {
                    "type": "ValueError",
                    "message": "predictor matrix is\nrank deficient",
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    assert (
        mod.main(
            [
                "inspect-audit-error",
                "--audit",
                str(audit),
                "--output",
                str(tmp_path / "error.tsv"),
            ]
        )
        == 0
    )
    observed = dict(
        pandas.read_csv(tmp_path / "error.tsv", sep="\t", keep_default_na=False).itertuples(index=False, name=None)
    )
    assert observed["status"] == "not_estimable"
    assert observed["reason"] == "predictor matrix is rank deficient"


def test_nhx_detection_and_reconciliation_preflight(tmp_path: Path):
    mod = load_module()
    nhx = tmp_path / "tree.nhx"
    nhx.write_text("(Genus_a_g1:1,Genus_b_g1:1)[&&NHX:D=N];\n", encoding="utf-8")
    plain = tmp_path / "tree.nwk"
    plain.write_text("(Genus_a_g1:1,Genus_b_g1:1);\n", encoding="utf-8")
    assert mod.main(["has-nhx", "--tree", str(nhx)]) == 0
    assert mod.main(["has-nhx", "--tree", str(plain)]) == 1

    write_tsv(
        tmp_path / "reconciliation.tsv",
        [
            {
                "event_type": "speciation",
                "eligible": "yes",
                "coverage_status": "complete",
                "species_event_id": "a|b",
            },
            {
                "event_type": "speciation",
                "eligible": "yes",
                "coverage_status": "complete",
                "species_event_id": "a|c",
            },
        ],
    )
    mod.main(
        [
            "inspect-reconciliation",
            "--reconciliation",
            str(tmp_path / "reconciliation.tsv"),
            "--min-species-events",
            "3",
            "--output",
            str(tmp_path / "preflight.tsv"),
        ]
    )
    preflight = dict(
        pandas.read_csv(tmp_path / "preflight.tsv", sep="\t", keep_default_na=False).itertuples(index=False, name=None)
    )
    assert preflight["status"] == "not_estimable"
    assert "eligible_species_events=2" in preflight["reason"]

    mod.main(
        [
            "empty-bundle",
            "--output-prefix",
            str(tmp_path / "empty"),
            "--status-output",
            str(tmp_path / "empty.status.tsv"),
            "--audit-output",
            str(tmp_path / "empty.audit.jsonl"),
            "--tree-id",
            "OG-empty",
            "--reason",
            "too_few_events",
            "--reconciliation",
            str(tmp_path / "reconciliation.tsv"),
        ]
    )
    for suffix in mod.REGRESSION_BUNDLE_SUFFIXES:
        assert pandas.read_csv(tmp_path / f"empty{suffix}", sep="\t").columns[0] == "analysis_id"
    empty_reconciliation = pandas.read_csv(tmp_path / "empty.reconciliation.tsv", sep="\t")
    assert set(empty_reconciliation["analysis_id"]) == {"preflight"}
    empty_status = pandas.read_csv(tmp_path / "empty.status.tsv", sep="\t")
    assert empty_status.loc[0, "status"] == "not_estimable"
    assert "too_few_events" in empty_status.loc[0, "reason"]
