from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas
import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support" / "reconciled_speciation_contrast.py"


def run(*args: object) -> None:
    subprocess.run([str(arg) for arg in args], check=True, text=True)


def test_reconciled_speciation_contrast_end_to_end_with_nwkit(tmp_path: Path):
    gene_names = [f"Genus_{letter}_g1" for letter in "abcde"]
    species_names = [f"Genus_{letter}" for letter in "abcde"]
    gene_tree = tmp_path / "gene.nwk"
    species_tree = tmp_path / "species.nwk"
    expression = tmp_path / "expression.tsv"
    traits = tmp_path / "traits.tsv"
    gene_tree.write_text(
        "((({}:1,{}:1):1,{}:2):1,({}:1,{}:1):2);".format(*gene_names),
        encoding="utf-8",
    )
    species_tree.write_text(
        "((({}:1,{}:1):1,{}:2):1,({}:1,{}:1):2);".format(*species_names),
        encoding="utf-8",
    )
    pandas.DataFrame(
        {
            "gene": gene_names,
            "expression_1": [1.8, 4.8, 6.8, 7.8, 11.8],
            "expression_2": [2.2, 5.2, 7.2, 8.2, 12.2],
        }
    ).to_csv(expression, sep="\t", index=False)
    pandas.DataFrame(
        {
            "species": species_names,
            "body_size": [1.0, 2.0, 4.0, 3.0, 7.0],
            "habitat": ["wet", "dry", "dry", "wet", "dry"],
        }
    ).to_csv(traits, sep="\t", index=False)

    prepared_expression = tmp_path / "prepared-expression.tsv"
    prepared_traits = tmp_path / "prepared-traits.tsv"
    plan = tmp_path / "plan.tsv"
    metadata = tmp_path / "metadata.tsv"
    run(
        sys.executable,
        SUPPORT,
        "prepare",
        "--expression",
        expression,
        "--species-traits",
        traits,
        "--responses",
        "all",
        "--predictors",
        "all",
        "--predictor-mode",
        "separate",
        "--expression-output",
        prepared_expression,
        "--species-traits-output",
        prepared_traits,
        "--analysis-plan-output",
        plan,
        "--metadata-output",
        metadata,
    )

    reconciliation = tmp_path / "preflight.reconciliation.tsv"
    preflight_audit = tmp_path / "preflight.audit.jsonl"
    preflight_metadata = tmp_path / "preflight.tsv"
    run(
        "nwkit",
        "reconcile",
        "--infile",
        gene_tree,
        "--species-tree",
        species_tree,
        "--tree-id",
        "OG1",
        "--event-source",
        "lca",
        "--species-parser",
        "legacy",
        "--outfile",
        reconciliation,
        "--audit",
        preflight_audit,
    )
    run(
        sys.executable,
        SUPPORT,
        "inspect-reconciliation",
        "--reconciliation",
        reconciliation,
        "--coverage",
        "complete",
        "--min-species-events",
        "2",
        "--output",
        preflight_metadata,
    )
    preflight = dict(
        pandas.read_csv(preflight_metadata, sep="\t", keep_default_na=False).itertuples(index=False, name=None)
    )
    assert preflight["status"] == "ready"

    pgls_common = [
        "--gene-tree",
        gene_tree,
        "--reconciliation-tree",
        gene_tree,
        "--species-tree",
        species_tree,
        "--expression",
        prepared_expression,
        "--species-traits",
        prepared_traits,
        "--responses",
        "expression",
        "--tree-id",
        "OG1",
        "--event-source",
        "lca",
        "--species-parser",
        "legacy",
        "--gene-evolution-model",
        "lambda",
        "--gene-evolution-parameter",
        "auto",
        "--species-evolution-model",
        "lambda",
        "--species-evolution-parameter",
        "auto",
        "--response-biological-id",
        "biological_id",
        "--response-within-variance",
        "pooled",
        "--event-weighting",
        "event",
        "--speciation-coverage",
        "complete",
        "--reconciled-model",
        "hierarchical",
        "--inference",
        "wald",
        "--predictor-categorical-replicate-policy",
        "latent",
        "--categorical-origin-diagnostics",
        "none",
    ]
    bundle_rows = []
    for analysis_id, predictor, categorical in (
        ("p001_body_size", "body_size", False),
        ("p002_habitat", "habitat", True),
    ):
        prefix = tmp_path / analysis_id
        audit = tmp_path / f"{analysis_id}.audit.jsonl"
        args = [
            "nwkit",
            "regress",
            *pgls_common,
            "--predictors",
            predictor,
            "--out-prefix",
            prefix,
            "--audit",
            audit,
        ]
        if categorical:
            args.extend(["--categorical-predictors", predictor])
        run(*args)
        bundle_rows.append({"analysis_id": analysis_id, "prefix": prefix, "audit": audit})

    bundle_list = tmp_path / "bundles.tsv"
    pandas.DataFrame(bundle_rows).to_csv(bundle_list, sep="\t", index=False)
    combined = tmp_path / "combined"
    status = tmp_path / "status.tsv"
    combined_audit = tmp_path / "combined.audit.jsonl"
    run(
        sys.executable,
        SUPPORT,
        "aggregate",
        "--bundle-list",
        bundle_list,
        "--output-prefix",
        combined,
        "--status-output",
        status,
        "--audit-output",
        combined_audit,
        "--source-audit",
        preflight_audit,
        "--tree-id",
        "OG1",
    )

    result = pandas.read_csv(tmp_path / "combined.regression.tsv", sep="\t")
    assert not result.empty
    assert set(result["analysis_id"]) == {"p001_body_size", "p002_habitat"}
    family_status = pandas.read_csv(status, sep="\t")
    assert family_status.loc[0, "n_analyses"] == 2
    assert (tmp_path / "combined.reconciliation.tsv").is_file()
    assert len(combined_audit.read_text(encoding="utf-8").splitlines()) == 3

    # NWKIT known-SE mode consumes one summarized mean/SE row per leaf and is
    # intentionally distinct from the raw biological-replicate mode above.
    summarized_expression = tmp_path / "summarized-expression.tsv"
    summarized_traits = tmp_path / "summarized-traits.tsv"
    summarized_metadata = tmp_path / "summarized-samples.tsv"
    pandas.DataFrame(
        {
            "gene": gene_names,
            "expression_mean": [2.0, 5.0, 7.0, 8.0, 12.0],
            "expression_se": [0.2] * 5,
            "expression_n": [2] * 5,
        }
    ).to_csv(summarized_expression, sep="\t", index=False)
    pandas.DataFrame(
        {
            "column": ["expression_mean"],
            "response": ["expression"],
            "standard_error_column": ["expression_se"],
            "sample_size_column": ["expression_n"],
        }
    ).to_csv(summarized_metadata, sep="\t", index=False)
    pandas.DataFrame(
        {
            "species": species_names,
            "body_size": [1.0, 2.0, 4.0, 3.0, 7.0],
            "body_size_se": [0.1] * 5,
        }
    ).to_csv(summarized_traits, sep="\t", index=False)

    known_expression = tmp_path / "known-expression.tsv"
    known_traits = tmp_path / "known-traits.tsv"
    run(
        sys.executable,
        SUPPORT,
        "prepare",
        "--expression",
        summarized_expression,
        "--species-traits",
        summarized_traits,
        "--sample-metadata",
        summarized_metadata,
        "--responses",
        "expression",
        "--predictors",
        "body_size",
        "--within-variance",
        "known-se",
        "--predictor-within-variance",
        "known-se",
        "--predictor-standard-error-columns",
        "body_size_se",
        "--expression-output",
        known_expression,
        "--species-traits-output",
        known_traits,
        "--analysis-plan-output",
        tmp_path / "known-plan.tsv",
        "--metadata-output",
        tmp_path / "known-metadata.tsv",
    )
    known_prefix = tmp_path / "known-se"
    run(
        "nwkit",
        "regress",
        "--gene-tree",
        gene_tree,
        "--reconciliation-tree",
        gene_tree,
        "--species-tree",
        species_tree,
        "--expression",
        known_expression,
        "--species-traits",
        known_traits,
        "--responses",
        "expression",
        "--predictors",
        "body_size",
        "--tree-id",
        "OG1-known-se",
        "--out-prefix",
        known_prefix,
        "--event-source",
        "lca",
        "--species-parser",
        "legacy",
        "--response-within-variance",
        "known-se",
        "--response-standard-error-columns",
        "expression__standard_error",
        "--response-sample-size-columns",
        "expression__sample_size",
        "--predictor-within-variance",
        "known-se",
        "--predictor-standard-error-columns",
        "body_size_se",
        "--gene-evolution-model",
        "lambda",
        "--gene-evolution-parameter",
        "auto",
        "--species-evolution-model",
        "lambda",
        "--species-evolution-parameter",
        "auto",
    )
    response_summary = pandas.read_csv(tmp_path / "known-se.response-tip-summary.tsv", sep="\t")
    predictor_summary = pandas.read_csv(tmp_path / "known-se.predictor-tip-summary.tsv", sep="\t")
    assert set(response_summary["variance_method"]) == {"known-se"}
    assert set(predictor_summary["variance_method"]) == {"known-se"}


def test_species_tree_comparators_use_matched_species_inputs(tmp_path: Path):
    species = [f"Genus_{letter}" for letter in "abcde"]
    outgroup = "Genus_f"
    species_tree = tmp_path / "species.nwk"
    species_tree.write_text(
        "(((({}:1,{}:1):1,{}:2):1,({}:1,{}:1):2):1,{}:4);".format(
            *species, outgroup
        ),
        encoding="utf-8",
    )
    genes = [f"{name}_g{copy}" for name in species for copy in (1, 2)]
    reconciliation = tmp_path / "reconciliation.tsv"
    pandas.DataFrame(
        {
            "node_class": ["tip"] * len(genes),
            "gene_name": genes,
            "species_name": [name for name in species for _copy in (1, 2)],
        }
    ).to_csv(reconciliation, sep="\t", index=False)

    expression_rows = []
    for species_index, name in enumerate(species):
        for copy, copy_shift in ((1, 0.1), (2, 0.6)):
            expression_rows.append(
                {
                    "leaf_name": f"{name}_g{copy}",
                    "expression": 1.0 + species_index + copy_shift,
                }
            )
    expression = tmp_path / "prepared-expression.tsv"
    pandas.DataFrame(expression_rows).to_csv(expression, sep="\t", index=False)
    traits = tmp_path / "prepared-traits.tsv"
    pandas.DataFrame(
        {
            "leaf_name": [*species, outgroup],
            "body_size": [1.0, 2.2, 3.7, 4.1, 6.4, 9.0],
            "habitat": ["land", "water", "land", "water", "land", "water"],
            "stage": [0, 1, 2, 1, 2, 0],
        }
    ).to_csv(traits, sep="\t", index=False)
    plan = tmp_path / "plan.tsv"
    pandas.DataFrame(
        [
            {
                "analysis_id": "p001_body_size",
                "predictors": "body_size",
                "categorical_predictors": ".",
                "ordered_predictors": ".",
                "factor_reference": ".",
                "predictor_standard_error_columns": ".",
                "predictor_sample_size_columns": ".",
                "origin_applicable": "no",
            },
            {
                "analysis_id": "p002_habitat",
                "predictors": "habitat",
                "categorical_predictors": "habitat",
                "ordered_predictors": ".",
                "factor_reference": ".",
                "predictor_standard_error_columns": ".",
                "predictor_sample_size_columns": ".",
                "origin_applicable": "yes",
            },
            {
                "analysis_id": "p003_stage",
                "predictors": "stage",
                "categorical_predictors": ".",
                "ordered_predictors": "stage=0|1|2",
                "factor_reference": ".",
                "predictor_standard_error_columns": ".",
                "predictor_sample_size_columns": ".",
                "origin_applicable": "yes",
            },
        ]
    ).to_csv(plan, sep="\t", index=False)
    metadata = tmp_path / "metadata.tsv"
    pandas.DataFrame(
        [
            ("responses", "expression"),
            ("response_biological_id", ""),
            ("response_technical_id", ""),
            ("response_batch", ""),
            ("standard_error_columns", ""),
        ],
        columns=["key", "value"],
    ).to_csv(metadata, sep="\t", index=False)

    prefix = tmp_path / "comparison"
    run(
        sys.executable,
        Path(__file__).resolve().parents[1] / "support" / "species_tree_pgls.py",
        "--methods",
        "species-nwkit,species-rphylopars",
        "--tree-id",
        "OG1",
        "--species-tree",
        species_tree,
        "--reconciliation",
        reconciliation,
        "--expression",
        expression,
        "--species-traits",
        traits,
        "--analysis-plan",
        plan,
        "--metadata",
        metadata,
        "--aggregation",
        "sum",
        "--expression-value-type",
        "identity",
        "--within-variance",
        "pooled",
        "--response-evolution-model",
        "brownian",
        "--predictor-evolution-model",
        "brownian",
        "--rphylopars-script",
        Path(__file__).resolve().parents[1] / "support" / "species_tree_rphylopars.R",
        "--native-out",
        f"{prefix}.native.tsv",
        "--rphylopars-out",
        f"{prefix}.rphylopars.tsv",
        "--comparison-out",
        f"{prefix}.comparison.tsv",
        "--status-out",
        f"{prefix}.status.tsv",
        "--audit-out",
        f"{prefix}.audit.jsonl",
        "--expression-summary-out",
        f"{prefix}.expression.tsv",
        "--expression-audit-out",
        f"{prefix}.expression-audit.tsv",
        "--response-tip-summary-out",
        f"{prefix}.response-summary.tsv",
        "--response-sampling-covariance-out",
        f"{prefix}.response-covariance.tsv",
        "--predictor-tip-summary-out",
        f"{prefix}.predictor-summary.tsv",
        "--predictor-sampling-covariance-out",
        f"{prefix}.predictor-covariance.tsv",
    )
    native = pandas.read_csv(f"{prefix}.native.tsv", sep="\t")
    rphylopars = pandas.read_csv(f"{prefix}.rphylopars.tsv", sep="\t")
    comparison = pandas.read_csv(f"{prefix}.comparison.tsv", sep="\t")
    status = pandas.read_csv(f"{prefix}.status.tsv", sep="\t")
    assert "body_size" in set(native["term"])
    assert native["term"].astype(str).str.contains("habitat").any()
    assert native["term"].astype(str).str.contains("stage").any()
    stage_rows = native.query("analysis_id == 'p003_stage' and predictor_type != 'intercept'")
    assert set(stage_rows["predictor_type"]) == {"ordered"}
    assert set(native["n_species"]) == {len(species)}
    assert "body_size" in set(rphylopars["term"])
    assert set(comparison["analysis_method"]) == {
        "species_nwkit",
        "species_rphylopars",
    }
    assert set(status.query("analysis_method == 'species_nwkit'")["status"]) == {"ok"}
    rphylopars_status = status.query("analysis_method == 'species_rphylopars'")
    assert rphylopars_status.query("analysis_id == 'p001_body_size'")["status"].tolist() == ["ok"]
    assert rphylopars_status.query("analysis_id == 'p002_habitat'")["status"].tolist() == ["not_estimable"]
    assert rphylopars_status.query("analysis_id == 'p003_stage'")["status"].tolist() == ["not_estimable"]
    assert (
        rphylopars_status.query("analysis_id == 'p002_habitat'")["reason"].str.contains("continuous predictors").all()
    )
    assert status.query("analysis_method == 'rsc'")["status"].tolist() == ["not_requested"]


def test_rsc_only_method_does_not_require_species_aggregation_inputs(tmp_path: Path):
    rsc_results = tmp_path / "rsc.tsv"
    pandas.DataFrame(
        [
            {
                "analysis_id": "p001_size",
                "tree_id": "OG1",
                "response": "expression",
                "term": "size",
                "coefficient": 1.2,
                "standard_error": 0.3,
                "p_value": 0.02,
                "inference_status": "ok",
            }
        ]
    ).to_csv(rsc_results, sep="\t", index=False)
    rsc_status = tmp_path / "rsc-status.tsv"
    pandas.DataFrame([{"tree_id": "OG1", "status": "ok", "reason": "", "n_result_rows": 1}]).to_csv(
        rsc_status, sep="\t", index=False
    )
    prefix = tmp_path / "rsc-only"
    run(
        sys.executable,
        Path(__file__).resolve().parents[1] / "support" / "species_tree_pgls.py",
        "--methods",
        "rsc",
        "--tree-id",
        "OG1",
        "--species-tree",
        tmp_path / "absent-tree.nwk",
        "--reconciliation",
        tmp_path / "absent-reconciliation.tsv",
        "--expression",
        tmp_path / "absent-expression.tsv",
        "--species-traits",
        tmp_path / "absent-traits.tsv",
        "--analysis-plan",
        tmp_path / "absent-plan.tsv",
        "--metadata",
        tmp_path / "absent-metadata.tsv",
        "--response-evolution-model",
        "brownian",
        "--predictor-evolution-model",
        "brownian",
        "--rphylopars-script",
        Path(__file__).resolve().parents[1] / "support" / "species_tree_rphylopars.R",
        "--rsc-results",
        rsc_results,
        "--rsc-status",
        rsc_status,
        "--native-out",
        f"{prefix}.native.tsv",
        "--rphylopars-out",
        f"{prefix}.rphylopars.tsv",
        "--comparison-out",
        f"{prefix}.comparison.tsv",
        "--status-out",
        f"{prefix}.status.tsv",
        "--audit-out",
        f"{prefix}.audit.jsonl",
        "--expression-summary-out",
        f"{prefix}.expression.tsv",
        "--expression-audit-out",
        f"{prefix}.expression-audit.tsv",
        "--response-tip-summary-out",
        f"{prefix}.response-summary.tsv",
        "--response-sampling-covariance-out",
        f"{prefix}.response-covariance.tsv",
        "--predictor-tip-summary-out",
        f"{prefix}.predictor-summary.tsv",
        "--predictor-sampling-covariance-out",
        f"{prefix}.predictor-covariance.tsv",
    )
    comparison = pandas.read_csv(f"{prefix}.comparison.tsv", sep="\t")
    status = pandas.read_csv(f"{prefix}.status.tsv", sep="\t")
    assert set(comparison["analysis_method"]) == {"rsc"}
    assert status.query("analysis_method == 'rsc'")["status"].tolist() == ["ok"]
    assert set(status.query("analysis_method != 'rsc'")["status"]) == {"not_requested"}


@pytest.mark.parametrize(
    (
        "model",
        "parameter",
        "predictor_model",
        "predictor_parameter",
        "parameter_name",
        "expected_parameter",
        "expected_parameter_status",
        "expected_status",
        "inference",
        "expected_reason",
        "predictor_sampling_variance",
    ),
    [
        ("lambda", "0.6", "lambda", "0.6", "lambda", 0.6, "fixed", "ok", "wald", "", 0.02),
        ("eb", "-0.1", "eb", "-0.1", "rate_change", -0.1, "fixed", "ok", "wald", "", 0.02),
        ("lambda", "auto", "lambda", "auto", "lambda", None, "estimated", "ok", "wald", "", 0.02),
        (
            "lambda",
            "0.6",
            "brownian",
            "auto",
            "",
            None,
            "",
            "not_estimable",
            "wald",
            "joint evolutionary model",
            0.02,
        ),
        (
            "brownian",
            "auto",
            "brownian",
            "auto",
            "",
            None,
            "",
            "not_estimable",
            "parametric-bootstrap",
            "does not implement requested inference",
            0.02,
        ),
        (
            "brownian",
            "auto",
            "brownian",
            "auto",
            "",
            None,
            "",
            "not_estimable",
            "wald",
            "cannot fit a singular sampling-error matrix",
            0.0,
        ),
    ],
)
def test_rphylopars_comparator_handles_shape_parameters_and_rejects_unmatched_models(
    tmp_path: Path,
    model: str,
    parameter: str,
    predictor_model: str,
    predictor_parameter: str,
    parameter_name: str,
    expected_parameter: float | None,
    expected_parameter_status: str,
    expected_status: str,
    inference: str,
    expected_reason: str,
    predictor_sampling_variance: float,
):
    species = list("ABCDEFGH")
    tree = tmp_path / "species.nwk"
    tree.write_text(
        "(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);",
        encoding="utf-8",
    )
    summary_rows = []
    for index, name in enumerate(species):
        summary_rows.extend(
            [
                {
                    "aggregation": "sum",
                    "source": "response",
                    "leaf_name": name,
                    "trait": "expression",
                    "value": 1.0 + 1.4 * index + (0.2 if index % 2 else -0.2),
                    "sampling_variance": 0.01,
                    "has_offdiagonal_sampling_covariance": "no",
                },
                {
                    "aggregation": "sum",
                    "source": "predictor:p001_size",
                    "leaf_name": name,
                    "trait": "size",
                    "value": float(index),
                    "sampling_variance": predictor_sampling_variance,
                    "has_offdiagonal_sampling_covariance": "no",
                },
            ]
        )
    summary = tmp_path / "summary.tsv"
    pandas.DataFrame(summary_rows).to_csv(summary, sep="\t", index=False)
    plan = tmp_path / "plan.tsv"
    pandas.DataFrame(
        [
            {
                "analysis_id": "p001_size",
                "predictors": "size",
                "categorical_predictors": ".",
                "ordered_predictors": ".",
            }
        ]
    ).to_csv(plan, sep="\t", index=False)
    result = tmp_path / "result.tsv"
    status = tmp_path / "status.tsv"
    run(
        "Rscript",
        Path(__file__).resolve().parents[1] / "support" / "species_tree_rphylopars.R",
        f"--tree={tree}",
        f"--summary={summary}",
        f"--plan={plan}",
        "--responses=expression",
        "--tree_id=OG1",
        f"--model={model}",
        f"--parameter={parameter}",
        f"--predictor_model={predictor_model}",
        f"--predictor_parameter={predictor_parameter}",
        "--branch_length=original",
        "--predictor_branch_length=original",
        "--reml=yes",
        "--confidence_level=0.95",
        f"--inference={inference}",
        "--sampling_covariance=require-diagonal",
        f"--outfile={result}",
        f"--status_out={status}",
    )
    fitted = pandas.read_csv(result, sep="\t")
    method_status = pandas.read_csv(status, sep="\t")
    assert set(method_status["status"]) == {expected_status}
    if expected_status == "not_estimable":
        assert fitted.empty
        assert method_status["reason"].str.contains(expected_reason).all()
        return
    assert set(fitted["evolution_parameter_name"]) == {parameter_name}
    assert set(fitted["evolution_parameter_status"]) == {expected_parameter_status}
    if expected_parameter is None:
        assert pandas.to_numeric(fitted["evolution_parameter"], errors="coerce").notna().all()
    else:
        assert fitted["evolution_parameter"].iloc[0] == pytest.approx(expected_parameter)
    assert set(fitted["optimizer_converged"]) == {"not_reported"}
