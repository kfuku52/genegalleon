from __future__ import annotations

import math

import pandas
import pytest
from scipy import sparse

from workflow.support.species_tree_pgls import (
    _aggregate_values,
    _covariance_diagonal_and_offdiagonal,
    _parse_aggregations,
    _parse_methods,
    aggregate_species_expression,
    summarize_for_stat_tree,
)


def test_method_and_aggregation_selection_are_explicit():
    assert _parse_methods("all") == ["rsc", "species-nwkit", "species-rphylopars"]
    assert _parse_methods("species-nwkit,rsc") == ["species-nwkit", "rsc"]
    assert _parse_aggregations("all") == ["sum", "mean", "max"]
    with pytest.raises(ValueError, match="Unknown PGLS"):
        _parse_methods("legacy")
    with pytest.raises(ValueError, match="duplicated"):
        _parse_methods("rsc,rsc")


@pytest.mark.parametrize(
    ("value_type", "values", "expected_sum", "expected_mean"),
    [
        ("identity", [2.0, 6.0], 8.0, 4.0),
        ("log", [math.log(2.0), math.log(6.0)], math.log(8.0), math.log(4.0)),
        ("log2", [1.0, math.log2(6.0)], 3.0, 2.0),
        ("log2p1", [1.0, math.log2(7.0)], 3.0, math.log2(4.5)),
    ],
)
def test_paralog_aggregation_occurs_on_linear_scale(value_type, values, expected_sum, expected_mean):
    total, _ = _aggregate_values(values, None, "sum", value_type)
    mean, _ = _aggregate_values(values, None, "mean", value_type)
    assert total == pytest.approx(expected_sum)
    assert mean == pytest.approx(expected_mean)


def test_known_standard_errors_are_propagated_for_sum_and_mean():
    total, total_se = _aggregate_values([2.0, 6.0], [0.3, 0.4], "sum", "identity")
    mean, mean_se = _aggregate_values([2.0, 6.0], [0.3, 0.4], "mean", "identity")
    assert total == 8.0
    assert total_se == pytest.approx(0.5)
    assert mean == 4.0
    assert mean_se == pytest.approx(0.25)


def test_sparse_sampling_covariance_is_inspected_without_dense_conversion():
    covariance = sparse.csr_matrix([[1.0, 0.0, 0.2], [0.0, 2.0, 0.0], [0.2, 0.0, 3.0]])
    diagonal, has_offdiagonal = _covariance_diagonal_and_offdiagonal(covariance, ["A", "B", "C"])
    assert diagonal.tolist() == [1.0, 2.0, 3.0]
    assert has_offdiagonal is True


def test_species_aggregation_never_treats_paralogs_as_replicates():
    expression = pandas.DataFrame(
        {
            "leaf_name": ["A_g1", "A_g2", "A_g1", "A_g2", "B_g1", "B_g1"],
            "expression": [1.0, 3.0, 2.0, 4.0, 10.0, 12.0],
            "biological_id": ["a1", "a1", "a2", "a2", "b1", "b2"],
        }
    )
    reconciliation = pandas.DataFrame(
        {
            "node_class": ["tip", "tip", "tip"],
            "gene_name": ["A_g1", "A_g2", "B_g1"],
            "species_name": ["A", "A", "B"],
        }
    )
    outputs, audit = aggregate_species_expression(
        expression,
        reconciliation,
        ["expression"],
        ["sum", "mean"],
        value_type="identity",
        missing_policy="error",
        tree_id="OG1",
    )
    summed = outputs["sum"].set_index(["leaf_name", "biological_id"])["expression"]
    averaged = outputs["mean"].set_index(["leaf_name", "biological_id"])["expression"]
    assert summed.loc[("A", "a1")] == 4.0
    assert summed.loc[("A", "a2")] == 6.0
    assert averaged.loc[("A", "a1")] == 2.0
    assert set(audit.loc[audit["species"] == "A", "expected_paralog_count"]) == {2}


def test_incomplete_paralog_measurement_is_not_silently_ignored():
    expression = pandas.DataFrame(
        {
            "leaf_name": ["A_g1", "A_g2"],
            "expression": [1.0, "NA"],
            "biological_id": ["a1", "a1"],
        }
    )
    reconciliation = pandas.DataFrame(
        {
            "node_class": ["tip", "tip"],
            "gene_name": ["A_g1", "A_g2"],
            "species_name": ["A", "A"],
        }
    )
    with pytest.raises(ValueError, match="Incomplete paralog expression coverage"):
        aggregate_species_expression(
            expression,
            reconciliation,
            ["expression"],
            ["sum"],
            value_type="identity",
            missing_policy="error",
            tree_id="OG1",
        )


def test_gene_tree_paralog_absent_from_expression_is_detected():
    expression = pandas.DataFrame({"leaf_name": ["A_g1"], "expression": [1.0], "biological_id": ["a1"]})
    reconciliation = pandas.DataFrame(
        {
            "node_class": ["tip", "tip"],
            "gene_name": ["A_g1", "A_g2"],
            "species_name": ["A", "A"],
        }
    )
    with pytest.raises(ValueError, match="1/2"):
        aggregate_species_expression(
            expression,
            reconciliation,
            ["expression"],
            ["sum"],
            value_type="identity",
            missing_policy="error",
            tree_id="OG1",
        )


def test_species_comparison_summary_is_bounded_and_method_specific(tmp_path):
    comparison = tmp_path / "comparison.tsv"
    status = tmp_path / "status.tsv"
    pandas.DataFrame(
        [
            {
                "analysis_method": "species_nwkit",
                "aggregation": "sum",
                "analysis_id": "p001_size",
                "response": "expression",
                "term": "size",
                "coefficient": 2.0,
                "standard_error": 0.2,
                "p_value": 0.01,
                "inference_status": "ok",
                "evolution_model": "brownian",
                "n_species": 20,
            },
            {
                "analysis_method": "species_rphylopars",
                "aggregation": "sum",
                "analysis_id": "p001_size",
                "response": "expression",
                "term": "size",
                "coefficient": 2.1,
                "standard_error": 0.25,
                "p_value": 0.02,
                "inference_status": "ok",
                "evolution_model": "brownian",
                "n_species": 20,
                "coefficient_difference_vs_species_nwkit": 0.1,
            },
        ]
    ).to_csv(comparison, sep="\t", index=False)
    pandas.DataFrame(
        [
            {"analysis_method": "species_nwkit", "status": "ok"},
            {"analysis_method": "species_rphylopars", "status": "ok"},
        ]
    ).to_csv(status, sep="\t", index=False)

    summary = summarize_for_stat_tree(comparison, status)
    assert summary["pgls_species_nwkit_num_ok"] == 1
    assert summary["pgls_species_nwkit_best_term"] == "size"
    assert summary["pgls_species_rphylopars_best_coefficient"] == pytest.approx(2.1)
    assert summary["pgls_nwkit_rphylopars_max_abs_coefficient_difference"] == pytest.approx(0.1)
