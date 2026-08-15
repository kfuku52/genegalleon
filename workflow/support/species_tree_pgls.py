#!/usr/bin/env python3
"""Run matched species-tree PGLS comparators for GeneGalleon RSC analyses.

This adapter deliberately lives in GeneGalleon.  NWKIT remains the numerical
implementation for the native species-tree fit; this module supplies the
gene-family-specific paralog aggregation, shared input preparation, method
orchestration, and comparison tables required by the workflow.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import math
import subprocess
import tempfile
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Iterable, Sequence

import numpy
import pandas
from scipy import sparse

METHODS = ("rsc", "species-nwkit", "species-rphylopars")
AGGREGATIONS = ("sum", "mean", "max")
MISSING_STRINGS = {"", "NA", "NaN", "nan", "?", "missing", "unknown", "."}

AUDIT_COLUMNS = [
    "tree_id",
    "aggregation",
    "response",
    "species",
    "biological_id",
    "technical_id",
    "batch",
    "expected_paralog_count",
    "observed_paralog_count",
    "coverage_fraction",
    "missing_policy",
    "paralog_covariance_mode",
    "paralog_covariance_pairs",
    "aggregated_value",
    "aggregated_standard_error",
]

STATUS_COLUMNS = [
    "tree_id",
    "analysis_method",
    "aggregation",
    "analysis_id",
    "response",
    "status",
    "reason",
    "n_species",
    "n_result_rows",
    "engine_version",
]

MINIMAL_RESULT_COLUMNS = [
    "analysis_method",
    "aggregation",
    "estimand",
    "analysis_id",
    "tree_id",
    "response",
    "term",
    "coefficient",
    "standard_error",
    "p_value",
    "inference_status",
]


def _is_missing(value: object) -> bool:
    if pandas.isna(value):
        return True
    return str(value).strip() in MISSING_STRINGS


def _read_tsv(path: Path, *, allow_empty: bool = False) -> pandas.DataFrame:
    frame = pandas.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        na_filter=False,
    )
    if frame.empty and not allow_empty:
        raise ValueError(f"TSV input has no rows: {path}")
    if frame.columns.duplicated().any():
        duplicated = frame.columns[frame.columns.duplicated()].tolist()
        raise ValueError(f"TSV input has duplicate columns: {', '.join(duplicated)}")
    return frame


def _read_metadata(path: Path) -> dict[str, str]:
    frame = _read_tsv(path)
    if list(frame.columns) != ["key", "value"]:
        raise ValueError("Preparation metadata must contain exactly key and value columns")
    if frame["key"].duplicated().any():
        raise ValueError("Preparation metadata contains duplicate keys")
    return {str(row.key): str(row.value) for row in frame.itertuples(index=False)}


def _csv(value: object) -> list[str]:
    return [item.strip() for item in str(value or "").split(",") if item.strip() and item.strip() != "."]


def _parse_methods(value: str) -> list[str]:
    methods = _csv(value)
    if methods == ["all"]:
        return list(METHODS)
    unknown = sorted(set(methods) - set(METHODS))
    if unknown:
        raise ValueError("Unknown PGLS method(s): " + ", ".join(unknown))
    if len(methods) != len(set(methods)):
        raise ValueError("PGLS methods must not be duplicated")
    if not methods:
        raise ValueError("At least one PGLS method is required")
    return methods


def _parse_aggregations(value: str) -> list[str]:
    if value == "all":
        return list(AGGREGATIONS)
    if value not in AGGREGATIONS:
        raise ValueError("Species-expression aggregation must be sum, mean, max, or all")
    return [value]


def _as_float(value: object, label: str) -> float:
    try:
        converted = float(str(value).strip())
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be numeric: {value!r}") from exc
    if not math.isfinite(converted):
        raise ValueError(f"{label} must be finite: {value!r}")
    return converted


def _to_linear(value: float, value_type: str) -> float:
    if value_type == "identity":
        return value
    if value_type == "log":
        return math.exp(value)
    if value_type == "log2":
        return 2.0**value
    if value_type == "log2p1":
        return 2.0**value - 1.0
    raise ValueError(f"Unsupported expression value type: {value_type}")


def _from_linear(value: float, value_type: str) -> float:
    if value_type == "identity":
        return value
    if value_type == "log":
        if value <= 0.0:
            raise ValueError("Log-scale paralog aggregation produced a non-positive value")
        return math.log(value)
    if value_type == "log2":
        if value <= 0.0:
            raise ValueError("Log2-scale paralog aggregation produced a non-positive value")
        return math.log2(value)
    if value_type == "log2p1":
        if value < 0.0:
            raise ValueError("Log2p1-scale paralog aggregation produced a negative value")
        return math.log2(value + 1.0)
    raise ValueError(f"Unsupported expression value type: {value_type}")


def _aggregate_values(
    values: Sequence[float],
    standard_errors: Sequence[float] | None,
    aggregation: str,
    value_type: str,
    sampling_covariance: numpy.ndarray | sparse.spmatrix | None = None,
) -> tuple[float, float | None]:
    linear = numpy.asarray([_to_linear(value, value_type) for value in values], dtype=float)
    if not numpy.isfinite(linear).all():
        raise ValueError("Paralog aggregation overflowed on the linear expression scale")
    if aggregation == "sum":
        combined_linear = float(linear.sum())
    elif aggregation == "mean":
        combined_linear = float(linear.mean())
    elif aggregation == "max":
        combined_linear = float(linear.max())
    else:
        raise ValueError(f"Unsupported paralog aggregation: {aggregation}")
    combined = _from_linear(combined_linear, value_type)
    if standard_errors is None:
        if sampling_covariance is not None:
            raise ValueError("Paralog sampling covariance requires aligned standard errors")
        return combined, None

    errors = numpy.asarray(standard_errors, dtype=float)
    if errors.shape != linear.shape or not numpy.isfinite(errors).all() or (errors < 0.0).any():
        raise ValueError("Paralog standard errors must be finite, non-negative, and aligned")
    validated_covariance = None
    if sampling_covariance is not None:
        validated_covariance = _validated_paralog_covariance(errors, sampling_covariance)

    if aggregation == "max":
        maxima = numpy.flatnonzero(linear == combined_linear)
        # The max function is not differentiable at a tie.  Using the largest
        # reported SE is deterministic under row reordering and is a
        # conservative approximation to selecting one of the tied paralogs.
        return combined, float(errors[maxima].max())

    if value_type == "identity":
        derivative = numpy.ones(len(linear), dtype=float)
        if aggregation == "mean":
            derivative /= len(linear)
    elif value_type in {"log", "log2"}:
        denominator = float(linear.sum())
        derivative = linear / denominator
    else:
        denominator = float(linear.sum()) + (len(linear) if aggregation == "mean" else 1.0)
        derivative = (linear + 1.0) / denominator
    diagonal_variance = float(numpy.square(derivative * errors).sum())
    if validated_covariance is None:
        return combined, math.sqrt(max(0.0, diagonal_variance))

    offdiagonal_variance = float(derivative @ validated_covariance @ derivative)
    variance = diagonal_variance + offdiagonal_variance
    return combined, math.sqrt(max(0.0, variance))


def _validated_paralog_covariance(
    errors: numpy.ndarray,
    sampling_covariance: numpy.ndarray | sparse.spmatrix,
) -> numpy.ndarray | sparse.csr_matrix:
    """Validate a zero-diagonal covariance update without forcing it dense."""

    shape = (len(errors), len(errors))
    if sparse.issparse(sampling_covariance):
        offdiagonal = sparse.csr_matrix(sampling_covariance, dtype=float)
        difference = offdiagonal - offdiagonal.T
        if (
            offdiagonal.shape != shape
            or not numpy.isfinite(offdiagonal.data).all()
            or (difference.nnz and float(numpy.max(numpy.abs(difference.data))) > 1e-12)
            or not numpy.allclose(offdiagonal.diagonal(), 0.0, rtol=0.0, atol=1e-12)
        ):
            raise ValueError(
                "Paralog sampling covariance must be a finite, symmetric, "
                "zero-diagonal matrix aligned to the reported standard errors"
            )
        full_covariance = offdiagonal + sparse.diags(numpy.square(errors), format="csr")
        scale = max(
            1.0,
            float(numpy.max(numpy.abs(full_covariance.data))) if full_covariance.nnz else 0.0,
        )
        tolerance = numpy.finfo(float).eps * scale * max(1, len(errors)) * 100.0
        if full_covariance.nnz == 0:
            minimum_eigenvalue = 0.0
        elif len(errors) <= 256:
            minimum_eigenvalue = float(numpy.linalg.eigvalsh(full_covariance.toarray()).min())
        else:
            minimum_eigenvalue = float(
                sparse.linalg.eigsh(
                    full_covariance,
                    k=1,
                    which="SA",
                    return_eigenvectors=False,
                    tol=1e-8,
                )[0]
            )
        validated: numpy.ndarray | sparse.csr_matrix = offdiagonal
    else:
        offdiagonal = numpy.asarray(sampling_covariance, dtype=float)
        if (
            offdiagonal.shape != shape
            or not numpy.isfinite(offdiagonal).all()
            or not numpy.allclose(offdiagonal, offdiagonal.T, rtol=1e-12, atol=1e-12)
            or not numpy.allclose(numpy.diag(offdiagonal), 0.0, rtol=0.0, atol=1e-12)
        ):
            raise ValueError(
                "Paralog sampling covariance must be a finite, symmetric, "
                "zero-diagonal matrix aligned to the reported standard errors"
            )
        full_covariance = offdiagonal + numpy.diag(numpy.square(errors))
        scale = max(1.0, float(numpy.max(numpy.abs(full_covariance))))
        tolerance = numpy.finfo(float).eps * scale * max(1, len(errors)) * 100.0
        minimum_eigenvalue = float(numpy.linalg.eigvalsh(full_covariance).min())
        validated = offdiagonal
    if minimum_eigenvalue < -tolerance:
        raise ValueError(
            "Paralog standard errors and sampling covariances do not form a positive-semidefinite covariance matrix"
        )
    return validated


def _prepare_paralog_sampling_covariance(
    frame: pandas.DataFrame | None,
    gene_species: dict[str, str],
    responses: Sequence[str],
    tree_id: str,
) -> dict[str, dict[tuple[str, str], float]]:
    """Validate optional within-species, cross-paralog known covariance."""
    if frame is None:
        return {}
    required = {"response", "gene_name_1", "gene_name_2", "sampling_covariance"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError("Paralog sampling covariance is missing columns: " + ", ".join(missing))
    selected = frame.copy()
    if "tree_id" in selected:
        selected = selected.loc[selected["tree_id"].astype(str) == str(tree_id)]
    response_set = set(responses)
    prepared: dict[str, dict[tuple[str, str], float]] = {}
    for record in selected.to_dict("records"):
        response = str(record["response"])
        first = str(record["gene_name_1"])
        second = str(record["gene_name_2"])
        known = (first in gene_species, second in gene_species)
        if not any(known):
            if "tree_id" in frame:
                raise ValueError(
                    "Paralog sampling covariance for tree {!r} contains unknown gene tips: {}, {}".format(
                        tree_id, first, second
                    )
                )
            continue
        if not all(known):
            raise ValueError(
                "Paralog sampling covariance pair mixes a mapped and unmapped gene tip: {}, {}".format(first, second)
            )
        if response not in response_set:
            raise ValueError("Paralog sampling covariance contains an unselected response: " + response)
        if first == second:
            raise ValueError(
                "Paralog sampling covariance diagonals come from standard errors; "
                "gene_name_1 and gene_name_2 must differ"
            )
        if gene_species[first] != gene_species[second]:
            raise ValueError(
                "Paralog sampling covariance cannot connect different species: {}, {}".format(first, second)
            )
        value = _as_float(record["sampling_covariance"], "paralog sampling covariance")
        pair = tuple(sorted((first, second)))
        by_pair = prepared.setdefault(response, {})
        if pair in by_pair:
            raise ValueError(
                "Paralog sampling covariance contains a duplicated unordered pair for response {!r}: {}, {}".format(
                    response, *pair
                )
            )
        by_pair[pair] = value
    return prepared


def _paralog_covariance_for_observations(
    leaf_names: Sequence[str],
    response: str,
    covariance_by_response: dict[str, dict[tuple[str, str], float]],
) -> tuple[sparse.csr_matrix | None, int]:
    by_pair = covariance_by_response.get(response, {})
    if not by_pair:
        return None, 0
    index_by_leaf = {str(leaf_name): index for index, leaf_name in enumerate(leaf_names)}
    row_indices: list[int] = []
    column_indices: list[int] = []
    values: list[float] = []
    for (first, second), value in by_pair.items():
        first_index = index_by_leaf.get(first)
        second_index = index_by_leaf.get(second)
        if first_index is None or second_index is None:
            continue
        row_indices.extend((first_index, second_index))
        column_indices.extend((second_index, first_index))
        values.extend((value, value))
    pair_count = len(values) // 2
    if not pair_count:
        return None, 0
    matrix = sparse.csr_matrix(
        (values, (row_indices, column_indices)),
        shape=(len(leaf_names), len(leaf_names)),
        dtype=float,
    )
    return matrix, pair_count


def _gene_species_map(reconciliation: pandas.DataFrame) -> dict[str, str]:
    required = {"node_class", "gene_name", "species_name"}
    missing = sorted(required - set(reconciliation.columns))
    if missing:
        raise ValueError("Reconciliation table is missing gene/species mapping columns: " + ", ".join(missing))
    tips = reconciliation.loc[reconciliation["node_class"].astype(str) == "tip", ["gene_name", "species_name"]].copy()
    tips = tips.loc[~tips["gene_name"].map(_is_missing) & ~tips["species_name"].map(_is_missing)]
    if tips.empty:
        raise ValueError("Reconciliation table contains no mapped gene tips")
    duplicated = tips.groupby("gene_name")["species_name"].nunique()
    conflicts = duplicated[duplicated > 1]
    if not conflicts.empty:
        raise ValueError("Reconciliation assigns a gene tip to multiple species: " + str(conflicts.index[0]))
    return dict(zip(tips["gene_name"].astype(str), tips["species_name"].astype(str), strict=True))


def _prune_tree_to_family_species(tree: Any, gene_species: dict[str, str]) -> list[str]:
    family_species = set(gene_species.values())
    tree_species = {str(name) for name in tree.leaf_names()}
    missing = sorted(family_species - tree_species)
    if missing:
        raise ValueError("Reconciliation contains species absent from the species tree: " + ", ".join(missing))
    if not family_species:
        raise ValueError("Reconciliation contains no family species")
    if family_species != tree_species:
        tree.prune(sorted(family_species), preserve_branch_length=True)
    return [str(name) for name in tree.leaf_names()]


def _subset_species_traits(frame: pandas.DataFrame, leaf_names: Sequence[str]) -> pandas.DataFrame:
    if "leaf_name" not in frame:
        raise ValueError("Prepared species traits require a leaf_name column")
    available = set(frame["leaf_name"].astype(str))
    missing = sorted(set(leaf_names) - available)
    if missing:
        raise ValueError("Prepared species traits lack family species: " + ", ".join(missing))
    return frame.loc[frame["leaf_name"].astype(str).isin(set(leaf_names))].copy()


def _write_tree(tree: Any, path: Path) -> None:
    from nwkit.util import write_tree

    write_tree(
        tree,
        SimpleNamespace(outfile=str(path), name_quote="auto"),
        format="auto",
        quiet=True,
    )


def _association_rows(frame: pandas.DataFrame) -> pandas.DataFrame:
    """Exclude intercept-only rows from association status and ranking."""
    if frame.empty:
        return frame
    mask = pandas.Series(True, index=frame.index, dtype=bool)
    intercept_labels = {"intercept", "(intercept)"}
    if "predictor_type" in frame:
        mask &= ~frame["predictor_type"].astype(str).str.strip().str.lower().isin(intercept_labels)
    for column in ("term", "source_term"):
        if column in frame:
            mask &= ~frame[column].astype(str).str.strip().str.lower().isin(intercept_labels)
    return frame.loc[mask]


def _usable_association_rows(frame: pandas.DataFrame) -> pandas.DataFrame:
    """Return converged, successful predictor-association rows."""
    selected = _association_rows(frame)
    if selected.empty:
        return selected
    mask = pandas.Series(True, index=selected.index, dtype=bool)
    if "inference_status" in selected:
        mask &= selected["inference_status"].astype(str).str.lower().eq("ok")
    false_values = {"no", "false", "0"}
    if "optimizer_converged" in selected:
        mask &= ~selected["optimizer_converged"].astype(str).str.lower().isin(false_values)
    for prefix in ("response_evolution", "predictor_evolution"):
        status_column = f"{prefix}_parameter_status"
        convergence_column = f"{prefix}_optimizer_converged"
        if status_column not in selected or convergence_column not in selected:
            continue
        estimated = selected[status_column].astype(str).str.lower().eq("estimated")
        failed = selected[convergence_column].astype(str).str.lower().isin(false_values)
        mask &= ~(estimated & failed)
    return selected.loc[mask]


def _adjust_association_p_values(
    values: pandas.Series,
) -> tuple[numpy.ndarray, numpy.ndarray]:
    """Return Holm family-wise and Benjamini-Hochberg adjusted p-values."""
    raw = numpy.asarray(values, dtype=float)
    count = len(raw)
    order = numpy.argsort(raw, kind="stable")
    ranked = raw[order]
    holm_ranked = numpy.minimum(
        1.0,
        numpy.maximum.accumulate((count - numpy.arange(count)) * ranked),
    )
    bh_ranked = numpy.minimum(
        1.0,
        numpy.minimum.accumulate((count / numpy.arange(count, 0, -1)) * ranked[::-1])[::-1],
    )
    holm = numpy.empty(count, dtype=float)
    bh = numpy.empty(count, dtype=float)
    holm[order] = holm_ranked
    bh[order] = bh_ranked
    return holm, bh


def aggregate_species_expression(
    expression: pandas.DataFrame,
    reconciliation: pandas.DataFrame,
    responses: Sequence[str],
    aggregations: Sequence[str],
    *,
    value_type: str,
    missing_policy: str,
    tree_id: str,
    paralog_sampling_covariance: pandas.DataFrame | None = None,
) -> tuple[dict[str, pandas.DataFrame], pandas.DataFrame]:
    if "leaf_name" not in expression:
        raise ValueError("Prepared expression requires a leaf_name column")
    mapping = _gene_species_map(reconciliation)
    covariance_by_response = _prepare_paralog_sampling_covariance(
        paralog_sampling_covariance,
        mapping,
        responses,
        tree_id,
    )
    unknown = sorted(set(expression["leaf_name"].astype(str)) - set(mapping))
    if unknown:
        raise ValueError("Prepared expression contains unmapped gene tips: " + ", ".join(unknown[:5]))
    working = expression.copy()
    working.insert(1, "species", working["leaf_name"].astype(str).map(mapping))
    expected = pandas.Series(list(mapping.values()), dtype=str).value_counts().to_dict()
    metadata_columns = [column for column in ("biological_id", "technical_id", "batch") if column in working]
    replicate_key = ["leaf_name", *metadata_columns]
    duplicate_mask = working.duplicated(replicate_key, keep=False)
    if duplicate_mask.any():
        duplicate = working.loc[duplicate_mask, replicate_key].iloc[0].to_dict()
        raise ValueError(
            "Prepared expression contains duplicate gene/replicate identity: "
            + ", ".join(f"{key}={value!r}" for key, value in duplicate.items())
        )
    group_columns = ["species", *metadata_columns]
    se_by_response = {
        response: f"{response}__standard_error" for response in responses if f"{response}__standard_error" in working
    }
    outputs: dict[str, pandas.DataFrame] = {}
    audit_rows: list[dict[str, object]] = []
    for aggregation in aggregations:
        rows: list[dict[str, object]] = []
        for key, group in working.groupby(group_columns, sort=False, dropna=False):
            key_values = key if isinstance(key, tuple) else (key,)
            identities = dict(zip(group_columns, key_values, strict=True))
            row: dict[str, object] = {"leaf_name": identities["species"]}
            row.update({column: identities[column] for column in metadata_columns})
            observed_any = False
            for response in responses:
                if response not in group:
                    raise ValueError(f"Prepared expression is missing response column {response!r}")
                observed = group.loc[~group[response].map(_is_missing)]
                if observed.empty:
                    row[response] = numpy.nan
                    continue
                observed_any = True
                expected_count = int(expected[str(identities["species"])])
                observed_count = int(observed["leaf_name"].nunique())
                if observed_count < expected_count and missing_policy == "error":
                    raise ValueError(
                        "Incomplete paralog expression coverage for response {!r}, species {!r}: {}/{}".format(
                            response, identities["species"], observed_count, expected_count
                        )
                    )
                values = [_as_float(value, f"expression response {response}") for value in observed[response]]
                errors = None
                se_column = se_by_response.get(response)
                if se_column:
                    if observed[se_column].map(_is_missing).any():
                        raise ValueError(f"Missing standard error while aggregating response {response!r}")
                    errors = [_as_float(value, f"standard error for {response}") for value in observed[se_column]]
                offdiagonal_covariance, covariance_pair_count = _paralog_covariance_for_observations(
                    observed["leaf_name"].astype(str).tolist(),
                    response,
                    covariance_by_response,
                )
                value, standard_error = _aggregate_values(
                    values,
                    errors,
                    aggregation,
                    value_type,
                    offdiagonal_covariance,
                )
                row[response] = value
                if se_column:
                    row[se_column] = standard_error
                audit_rows.append(
                    {
                        "tree_id": tree_id,
                        "aggregation": aggregation,
                        "response": response,
                        "species": identities["species"],
                        "biological_id": identities.get("biological_id", ""),
                        "technical_id": identities.get("technical_id", ""),
                        "batch": identities.get("batch", ""),
                        "expected_paralog_count": expected_count,
                        "observed_paralog_count": observed_count,
                        "coverage_fraction": observed_count / expected_count,
                        "missing_policy": missing_policy,
                        "paralog_covariance_mode": ("provided" if covariance_pair_count else "independent"),
                        "paralog_covariance_pairs": covariance_pair_count,
                        "aggregated_value": value,
                        "aggregated_standard_error": "" if standard_error is None else standard_error,
                    }
                )
            if observed_any:
                rows.append(row)
        columns = ["leaf_name", *responses, *metadata_columns, *se_by_response.values()]
        outputs[aggregation] = pandas.DataFrame(rows).reindex(columns=columns)
    return outputs, pandas.DataFrame(audit_rows, columns=AUDIT_COLUMNS)


def _parameter(value: str) -> float | None:
    return None if value in {"", "auto", "None", "none"} else _as_float(value, "evolution parameter")


def _ordered(value: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    if value in {"", "."}:
        return out
    for item in value.split(","):
        name, levels = item.split("=", 1)
        out[name] = levels.split("|")
    return out


def _references(value: str) -> dict[str, str]:
    out: dict[str, str] = {}
    if value in {"", "."}:
        return out
    for item in value.split(","):
        name, level = item.split("=", 1)
        out[name] = level
    return out


def _replicate_arguments(
    path: Path,
    *,
    biological_id: str,
    technical_id: str,
    batch: str,
    within_variance: str,
    technical_aggregation: str,
    standard_error_columns: str,
    sample_size_columns: str,
) -> SimpleNamespace:
    from nwkit.conventions import DEFAULT_TABLE_MISSING_VALUES_CSV

    return SimpleNamespace(
        trait=str(path),
        data=str(path),
        biological_id=biological_id or None,
        technical_id=technical_id or None,
        batch=batch or None,
        within_variance=within_variance,
        technical_aggregation=technical_aggregation,
        standard_error_columns=standard_error_columns or None,
        sample_size_columns=sample_size_columns or None,
        unmatched="error",
        missing_values=DEFAULT_TABLE_MISSING_VALUES_CSV,
    )


def _read_trait_inputs(
    path: Path,
    tree: Any,
    traits: Sequence[str],
    *,
    tree_id: str,
    biological_id: str,
    technical_id: str,
    batch: str,
    within_variance: str,
    technical_aggregation: str,
    standard_error_columns: str = "",
    sample_size_columns: str = "",
    categorical: Iterable[str] = (),
    ordered: Iterable[str] = (),
    categorical_policy: str = "latent",
):
    from nwkit.contrast import _read_mixed_replicate_traits
    from nwkit.ordinary_pgls import _ordinary_predictor_values
    from nwkit.replicates import TIP_SUMMARY_COLUMNS, ReplicateEstimates

    arguments = _replicate_arguments(
        path,
        biological_id=biological_id,
        technical_id=technical_id,
        batch=batch,
        within_variance=within_variance,
        technical_aggregation=technical_aggregation,
        standard_error_columns=standard_error_columns,
        sample_size_columns=sample_size_columns,
    )
    discrete = set(categorical) | set(ordered)
    replicate_requested = bool(biological_id) or within_variance == "known-se"
    if replicate_requested:
        return _read_mixed_replicate_traits(
            arguments,
            tree,
            list(traits),
            discrete,
            tree_id,
            categorical_policy=categorical_policy,
            option_name="GeneGalleon species-PGLS input",
        )
    return ReplicateEstimates(
        values_by_trait=_ordinary_predictor_values(
            arguments,
            tree,
            list(traits),
            duplicate_leaf_names="error",
            categorical=set(categorical),
            ordered=set(ordered),
        ),
        sampling_covariance_by_trait={},
        tip_summary=pandas.DataFrame(columns=TIP_SUMMARY_COLUMNS),
        model_by_trait={trait: "single-observation" for trait in traits},
    )


def _covariance_diagonal_and_offdiagonal(covariance: Any, leaf_names: Sequence[str]) -> tuple[numpy.ndarray, bool]:
    covariance_type = type(covariance)
    if covariance_type.__module__ == "nwkit.gaussian" and covariance_type.__name__ == "DiagonalLowRankCovariance":
        low_rank = covariance.low_rank
        if sparse.issparse(low_rank):
            low_rank = low_rank.tocsr()
            squared_norm = numpy.asarray(low_rank.multiply(low_rank).sum(axis=1)).ravel()
        else:
            low_rank = numpy.asarray(low_rank, dtype=float)
            squared_norm = numpy.sum(low_rank * low_rank, axis=1)
        diagonal = numpy.asarray(covariance.diagonal, dtype=float) + squared_norm
        # Inspect the low-rank Gram matrix in bounded blocks.  Materializing
        # the complete n-by-n product defeats the representation and used
        # hundreds of MiB for only 5,000 tips.
        n_rows = low_rank.shape[0]
        target_entries = 1_000_000
        block_size = max(1, min(n_rows, target_entries // max(1, n_rows)))
        for start in range(0, n_rows, block_size):
            stop = min(n_rows, start + block_size)
            gram = low_rank[start:stop] @ low_rank.T
            if sparse.issparse(gram):
                coordinates = gram.tocoo()
                offdiagonal = coordinates.col != coordinates.row + start
                has_offdiagonal = numpy.any(numpy.abs(coordinates.data[offdiagonal]) > 1e-12)
            else:
                gram[numpy.arange(stop - start), numpy.arange(start, stop)] = 0.0
                has_offdiagonal = numpy.any(numpy.abs(gram) > 1e-12)
            if has_offdiagonal:
                return diagonal, True
        return diagonal, False
    if isinstance(covariance, pandas.DataFrame):
        matrix = covariance.loc[list(leaf_names), list(leaf_names)].to_numpy(dtype=float)
    elif sparse.issparse(covariance):
        matrix = covariance.tocsr()
        diagonal = numpy.asarray(matrix.diagonal(), dtype=float)
        offdiagonal = matrix - sparse.diags(diagonal, format="csr")
        has_offdiagonal = bool(offdiagonal.nnz and numpy.any(numpy.abs(offdiagonal.data) > 1e-12))
        return diagonal, has_offdiagonal
    else:
        matrix = numpy.asarray(covariance, dtype=float)
    if matrix.ndim == 1:
        return matrix, False
    off = matrix.copy()
    numpy.fill_diagonal(off, 0.0)
    return numpy.diag(matrix), bool(numpy.any(numpy.abs(off) > 1e-12))


def _sampling_table(
    covariance_by_trait: dict[str, Any], leaf_names: Sequence[str], aggregation: str
) -> pandas.DataFrame:
    from nwkit.ordinary_pgls import _sampling_covariance_table

    table = _sampling_covariance_table(covariance_by_trait, list(leaf_names))
    table.insert(0, "aggregation", aggregation)
    return table


def _summary_long(
    aggregation: str,
    leaf_names: Sequence[str],
    response_estimates: Any,
    predictor_estimates_by_analysis: dict[str, Any],
) -> pandas.DataFrame:
    rows: list[dict[str, object]] = []
    sources = [("response", response_estimates)] + [
        (f"predictor:{analysis_id}", estimate) for analysis_id, estimate in predictor_estimates_by_analysis.items()
    ]
    seen: set[tuple[str, str]] = set()
    for source, estimates in sources:
        for trait, values in estimates.values_by_trait.items():
            key = (source if source == "response" else source, trait)
            if key in seen:
                continue
            seen.add(key)
            covariance = estimates.sampling_covariance_by_trait.get(trait)
            if covariance is None:
                diagonal = numpy.zeros(len(leaf_names), dtype=float)
                offdiagonal = False
            else:
                diagonal, offdiagonal = _covariance_diagonal_and_offdiagonal(covariance, leaf_names)
            for index, leaf_name in enumerate(leaf_names):
                rows.append(
                    {
                        "aggregation": aggregation,
                        "source": source,
                        "leaf_name": leaf_name,
                        "trait": trait,
                        "value": values[leaf_name],
                        "sampling_variance": float(diagonal[index]),
                        "has_offdiagonal_sampling_covariance": "yes" if offdiagonal else "no",
                    }
                )
    return pandas.DataFrame(rows)


def _native_status(
    tree_id: str,
    aggregation: str,
    analysis_id: str,
    responses: Sequence[str],
    frame: pandas.DataFrame | None,
    version: str,
    reason: str = "",
) -> list[dict[str, object]]:
    rows = []
    for response in responses:
        selected = pandas.DataFrame() if frame is None else frame.loc[frame["response"].astype(str) == response]
        usable = _usable_association_rows(selected)
        response_reason = reason
        if usable.empty and not response_reason:
            if selected.empty:
                response_reason = "nwkit_returned_no_response_rows"
            elif "inference_status" in selected:
                states = sorted(set(selected["inference_status"].dropna().astype(str)))
                if states and set(states) != {"ok"}:
                    response_reason = "nwkit_inference_status:" + ",".join(states)
                else:
                    response_reason = "nwkit_returned_no_usable_association_rows"
            else:
                response_reason = "nwkit_returned_no_usable_rows"
        rows.append(
            {
                "tree_id": tree_id,
                "analysis_method": "species_nwkit",
                "aggregation": aggregation,
                "analysis_id": analysis_id,
                "response": response,
                "status": "ok" if not usable.empty else "not_estimable",
                "reason": response_reason if usable.empty else "",
                "n_species": "" if selected.empty or "n_species" not in selected else selected["n_species"].max(),
                "n_result_rows": int(selected.shape[0]),
                "engine_version": version,
            }
        )
    return rows


def _run_rphylopars(
    args: argparse.Namespace,
    tree: Any,
    summary_long: pandas.DataFrame,
    plan: pandas.DataFrame,
    responses: Sequence[str],
) -> tuple[pandas.DataFrame, pandas.DataFrame, dict[str, object]]:
    with tempfile.TemporaryDirectory(prefix="genegalleon-rphylopars-") as temporary:
        directory = Path(temporary)
        summary_path = directory / "summary.tsv"
        plan_path = directory / "plan.tsv"
        tree_path = directory / "family-species-tree.nwk"
        results_path = directory / "results.tsv"
        status_path = directory / "status.tsv"
        summary_long.to_csv(summary_path, sep="\t", index=False, na_rep="NA")
        plan.to_csv(plan_path, sep="\t", index=False, na_rep="NA")
        _write_tree(tree, tree_path)
        command = [
            "Rscript",
            str(args.rphylopars_script),
            f"--tree={tree_path}",
            f"--summary={summary_path}",
            f"--plan={plan_path}",
            f"--responses={','.join(responses)}",
            f"--tree_id={args.tree_id}",
            f"--model={args.response_evolution_model}",
            f"--parameter={args.response_evolution_parameter}",
            f"--predictor_model={args.predictor_evolution_model}",
            f"--predictor_parameter={args.predictor_evolution_parameter}",
            f"--branch_length={args.branch_length}",
            f"--predictor_branch_length={args.predictor_branch_length}",
            f"--reml={args.reml}",
            f"--confidence_level={args.confidence_level}",
            f"--inference={args.inference}",
            f"--sampling_covariance={args.rphylopars_sampling_covariance}",
            f"--outfile={results_path}",
            f"--status_out={status_path}",
        ]
        completed = subprocess.run(command, check=False, text=True, capture_output=True)
        audit = {
            "analysis_method": "species_rphylopars",
            "command": command,
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
        }
        if completed.returncode != 0:
            raise RuntimeError(
                "Rphylopars adapter failed with exit status {}: {}".format(
                    completed.returncode,
                    " ".join(completed.stderr.split()) or "no diagnostic",
                )
            )
        results = pandas.read_csv(results_path, sep="\t", low_memory=False)
        status = pandas.read_csv(status_path, sep="\t", low_memory=False)
        if "engine_version" in status:
            audit["rphylopars_version"] = sorted(set(status["engine_version"].dropna().astype(str)))
        return results, status, audit


def _comparison_table(
    rsc_path: Path | None,
    native: pandas.DataFrame,
    rphylopars: pandas.DataFrame,
) -> pandas.DataFrame:
    frames: list[pandas.DataFrame] = []
    if rsc_path and rsc_path.is_file():
        rsc = pandas.read_csv(rsc_path, sep="\t", low_memory=False)
        if not rsc.empty:
            rsc.insert(0, "analysis_method", "rsc")
            rsc.insert(1, "aggregation", "gene_lineage_contrasts")
            rsc.insert(2, "estimand", "gene-lineage expression change per reconciled speciation event")
            rsc["directly_comparable_to_rsc"] = "reference"
            rsc["comparison_note"] = "RSC event-level estimand"
            frames.append(rsc)
    if not native.empty:
        current = native.copy()
        current["directly_comparable_to_rsc"] = "no"
        current["comparison_note"] = "species-tip estimand after within-sample paralog aggregation"
        frames.append(current)
    if not rphylopars.empty:
        current = rphylopars.copy()
        current["directly_comparable_to_rsc"] = "no"
        current["comparison_note"] = (
            "species-tip estimand; Rphylopars uses a joint trait model, and likelihood, "
            "parameter counting, and optimizer reporting are engine-specific"
        )
        frames.append(current)
    if not frames:
        return pandas.DataFrame(
            columns=[
                "analysis_method",
                "aggregation",
                "estimand",
                "analysis_id",
                "response",
                "term",
                "coefficient",
                "standard_error",
                "p_value",
                "directly_comparable_to_rsc",
                "comparison_note",
            ]
        )
    combined = pandas.concat(frames, ignore_index=True, sort=False)
    combined["coefficient_difference_vs_species_nwkit"] = numpy.nan
    keys = ["aggregation", "analysis_id", "response", "term"]
    if not native.empty and not rphylopars.empty and all(column in combined for column in keys):
        native_coefficients = native.drop_duplicates(keys, keep="first").set_index(keys)["coefficient"]
        for index, row in combined.loc[combined["analysis_method"] == "species_rphylopars"].iterrows():
            key = tuple(row[column] for column in keys)
            if key in native_coefficients.index:
                combined.at[index, "coefficient_difference_vs_species_nwkit"] = pandas.to_numeric(
                    row.get("coefficient"), errors="coerce"
                ) - pandas.to_numeric(native_coefficients.loc[key], errors="coerce")
    return combined


def _read_rsc_method_status(path: Path | None, nwkit_version: str = "") -> list[dict[str, object]]:
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return []
    status = pandas.read_csv(path, sep="\t", low_memory=False)
    if status.empty:
        return []
    row = status.iloc[0]
    raw_reason = row.get("reason", "")
    reason = "" if pandas.isna(raw_reason) else str(raw_reason)
    method_status = "not_requested" if "method_not_requested" in reason else row.get("status", "not_estimable")
    if method_status != "not_requested" and not nwkit_version:
        try:
            nwkit_version = importlib.metadata.version("nwkit")
        except importlib.metadata.PackageNotFoundError:
            nwkit_version = "unknown"
    return [
        {
            "tree_id": row.get("tree_id", ""),
            "analysis_method": "rsc",
            "aggregation": "gene_lineage_contrasts",
            "analysis_id": "all",
            "response": "all",
            "status": method_status,
            "reason": reason,
            "n_species": "",
            "n_result_rows": row.get("n_result_rows", 0),
            "engine_version": nwkit_version,
        }
    ]


def _append_unrequested_status(rows: list[dict[str, object]], methods: Sequence[str], tree_id: str) -> None:
    represented = {str(row.get("analysis_method", "")) for row in rows}
    for method in METHODS:
        normalized = method.replace("-", "_")
        if method in methods or normalized in represented:
            continue
        rows.append(
            {
                "tree_id": tree_id,
                "analysis_method": normalized,
                "aggregation": "",
                "analysis_id": "",
                "response": "",
                "status": "not_requested",
                "reason": "method_not_requested",
                "n_species": "",
                "n_result_rows": 0,
                "engine_version": "",
            }
        )


def _audit_status_record(row: dict[str, object]) -> dict[str, object]:
    normalized: dict[str, object] = {"event": "method_status"}
    for key, value in row.items():
        if pandas.isna(value):
            normalized[key] = None
        elif isinstance(value, numpy.generic):
            normalized[key] = value.item()
        else:
            normalized[key] = value
    return normalized


def _write_empty_outputs(args: argparse.Namespace, methods: Sequence[str], reason: str) -> int:
    native = pandas.DataFrame(columns=MINIMAL_RESULT_COLUMNS)
    rphylopars = pandas.DataFrame(columns=MINIMAL_RESULT_COLUMNS)
    status_rows = _read_rsc_method_status(args.rsc_status)
    for method in methods:
        if method == "rsc":
            continue
        status_rows.append(
            {
                "tree_id": args.tree_id,
                "analysis_method": method.replace("-", "_"),
                "aggregation": args.aggregation,
                "analysis_id": "",
                "response": "",
                "status": "not_estimable",
                "reason": reason,
                "n_species": "",
                "n_result_rows": 0,
                "engine_version": "",
            }
        )
    _append_unrequested_status(status_rows, methods, args.tree_id)
    for path in (
        args.native_out,
        args.rphylopars_out,
        args.comparison_out,
        args.status_out,
        args.audit_out,
        args.expression_summary_out,
        args.expression_audit_out,
        args.response_tip_summary_out,
        args.response_sampling_covariance_out,
        args.predictor_tip_summary_out,
        args.predictor_sampling_covariance_out,
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
    native.to_csv(args.native_out, sep="\t", index=False)
    rphylopars.to_csv(args.rphylopars_out, sep="\t", index=False)
    _comparison_table(args.rsc_results, native, rphylopars).to_csv(
        args.comparison_out, sep="\t", index=False, na_rep="NA"
    )
    pandas.DataFrame(status_rows, columns=STATUS_COLUMNS).to_csv(args.status_out, sep="\t", index=False, na_rep="NA")
    species_requested = any(method != "rsc" for method in methods)
    audit_records = [
        {
            "tree_id": args.tree_id,
            "analysis_method": "species_comparators",
            "status": "not_estimable" if species_requested else "not_requested",
            "reason": reason,
        },
        *(_audit_status_record(row) for row in status_rows),
    ]
    args.audit_out.write_text(
        "".join(json.dumps(record, sort_keys=True) + "\n" for record in audit_records),
        encoding="utf-8",
    )
    pandas.DataFrame(columns=["aggregation", "leaf_name", "response"]).to_csv(
        args.expression_summary_out, sep="\t", index=False
    )
    pandas.DataFrame(columns=AUDIT_COLUMNS).to_csv(args.expression_audit_out, sep="\t", index=False)
    pandas.DataFrame(columns=["aggregation", "tree_id", "leaf_name", "trait"]).to_csv(
        args.response_tip_summary_out, sep="\t", index=False
    )
    covariance_columns = [
        "aggregation",
        "tree_id",
        "trait",
        "leaf_name_1",
        "leaf_name_2",
        "sampling_covariance",
        "covariance_representation",
    ]
    pandas.DataFrame(columns=covariance_columns).to_csv(args.response_sampling_covariance_out, sep="\t", index=False)
    pandas.DataFrame(columns=["aggregation", "analysis_id", "tree_id", "leaf_name", "trait"]).to_csv(
        args.predictor_tip_summary_out, sep="\t", index=False
    )
    pandas.DataFrame(columns=["analysis_id", *covariance_columns]).to_csv(
        args.predictor_sampling_covariance_out, sep="\t", index=False
    )
    return 0


def summarize_for_stat_tree(comparison_path: str | Path, status_path: str | Path) -> dict[str, object]:
    """Return bounded species-comparator metrics for one ``stat_tree`` row."""
    out: dict[str, object] = {}
    status_file = Path(status_path)
    if status_file.is_file() and status_file.stat().st_size > 0:
        status = pandas.read_csv(status_file, sep="\t", low_memory=False)
        if {"analysis_method", "status"}.issubset(status.columns):
            for method in ("species_nwkit", "species_rphylopars"):
                selected = status.loc[status["analysis_method"].astype(str) == method]
                out[f"pgls_{method}_num_ok"] = int((selected["status"] == "ok").sum())
                out[f"pgls_{method}_num_not_estimable"] = int((selected["status"] == "not_estimable").sum())
                out[f"pgls_{method}_num_not_requested"] = int((selected["status"] == "not_requested").sum())

    comparison_file = Path(comparison_path)
    if not comparison_file.is_file() or comparison_file.stat().st_size == 0:
        return out
    comparison = pandas.read_csv(comparison_file, sep="\t", low_memory=False)
    out["pgls_comparison_num_result_rows"] = int(comparison.shape[0])
    if comparison.empty or "analysis_method" not in comparison:
        return out
    for method in ("species_nwkit", "species_rphylopars"):
        selected = comparison.loc[comparison["analysis_method"].astype(str) == method].copy()
        selected = _usable_association_rows(selected)
        if selected.empty or "p_value" not in selected:
            continue
        selected["_p_value"] = pandas.to_numeric(selected["p_value"], errors="coerce")
        selected = selected.loc[numpy.isfinite(selected["_p_value"]) & selected["_p_value"].between(0.0, 1.0)]
        if selected.empty:
            continue
        holm, bh = _adjust_association_p_values(selected["_p_value"])
        selected["_p_value_holm"] = holm
        selected["_p_value_bh"] = bh
        best = selected.loc[selected["_p_value"].idxmin()]
        for column in (
            "aggregation",
            "analysis_id",
            "response",
            "term",
            "coefficient",
            "standard_error",
            "p_value",
            "evolution_model",
            "evolution_parameter",
            "n_species",
        ):
            if column in best and not pandas.isna(best[column]):
                out[f"pgls_{method}_best_{column}"] = best[column]
        prefix = f"pgls_{method}"
        out[f"{prefix}_num_tested_associations"] = int(selected.shape[0])
        out[f"{prefix}_multiplicity_scope"] = "all_usable_family_associations_for_method"
        out[f"{prefix}_best_p_value_raw"] = float(best["_p_value"])
        out[f"{prefix}_best_p_value_holm"] = float(best["_p_value_holm"])
        out[f"{prefix}_best_p_value_bh"] = float(best["_p_value_bh"])
        # The unsuffixed field is safe for thresholding; the raw value remains
        # available under the explicit ``_raw`` name for descriptive ranking.
        out[f"{prefix}_best_p_value"] = float(best["_p_value_holm"])
        out[f"{prefix}_best_p_value_adjustment"] = "holm"
    if "coefficient_difference_vs_species_nwkit" in comparison:
        differences = pandas.to_numeric(comparison["coefficient_difference_vs_species_nwkit"], errors="coerce").abs()
        differences = differences.loc[numpy.isfinite(differences)]
        if not differences.empty:
            out["pgls_nwkit_rphylopars_max_abs_coefficient_difference"] = float(differences.max())
    return out


def run(args: argparse.Namespace) -> int:
    from nwkit.ordinary_pgls import fit_ordinary_pgls
    from nwkit.util import read_tree

    methods = _parse_methods(args.methods)
    aggregations = _parse_aggregations(args.aggregation)
    if args.empty_reason:
        return _write_empty_outputs(args, methods, args.empty_reason)
    if not any(method != "rsc" for method in methods):
        return _write_empty_outputs(args, methods, "species_methods_not_requested")
    for path in (
        args.native_out,
        args.rphylopars_out,
        args.comparison_out,
        args.status_out,
        args.audit_out,
        args.expression_summary_out,
        args.expression_audit_out,
        args.response_tip_summary_out,
        args.response_sampling_covariance_out,
        args.predictor_tip_summary_out,
        args.predictor_sampling_covariance_out,
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
    metadata = _read_metadata(args.metadata)
    responses = _csv(metadata.get("responses", ""))
    if not responses:
        raise ValueError("Prepared RSC metadata contains no usable responses")
    expression = _read_tsv(args.expression)
    reconciliation = _read_tsv(args.reconciliation)
    paralog_sampling_covariance = (
        None if args.paralog_sampling_covariance is None else _read_tsv(args.paralog_sampling_covariance)
    )
    plan = _read_tsv(args.analysis_plan)
    if plan.empty:
        raise ValueError("Prepared RSC analysis plan is empty")
    tree = read_tree(str(args.species_tree), "auto", True)
    leaf_names = _prune_tree_to_family_species(tree, _gene_species_map(reconciliation))
    species_traits = _subset_species_traits(_read_tsv(args.species_traits), leaf_names)
    aggregated, aggregation_audit = aggregate_species_expression(
        expression,
        reconciliation,
        responses,
        aggregations,
        value_type=args.expression_value_type,
        missing_policy=args.paralog_missing,
        tree_id=args.tree_id,
        paralog_sampling_covariance=paralog_sampling_covariance,
    )
    pandas.concat(
        [frame.assign(aggregation=name) for name, frame in aggregated.items()],
        ignore_index=True,
        sort=False,
    ).to_csv(args.expression_summary_out, sep="\t", index=False, na_rep="NA")
    aggregation_audit.to_csv(args.expression_audit_out, sep="\t", index=False, na_rep="NA")

    try:
        nwkit_version = importlib.metadata.version("nwkit")
    except importlib.metadata.PackageNotFoundError:
        nwkit_version = "unknown"
    native_frames: list[pandas.DataFrame] = []
    status_rows: list[dict[str, object]] = _read_rsc_method_status(args.rsc_status, nwkit_version)
    audit_records: list[dict[str, object]] = []
    response_summary_frames: list[pandas.DataFrame] = []
    response_covariance_frames: list[pandas.DataFrame] = []
    predictor_summary_frames: list[pandas.DataFrame] = []
    predictor_covariance_frames: list[pandas.DataFrame] = []
    rphylopars_summary_frames: list[pandas.DataFrame] = []

    for aggregation, frame in aggregated.items():
        with tempfile.TemporaryDirectory(prefix="genegalleon-species-pgls-") as temporary:
            expression_path = Path(temporary) / "species-expression.tsv"
            species_traits_path = Path(temporary) / "species-traits.tsv"
            frame.to_csv(expression_path, sep="\t", index=False, na_rep="NA")
            species_traits.to_csv(species_traits_path, sep="\t", index=False, na_rep="NA")
            response_estimates = _read_trait_inputs(
                expression_path,
                tree,
                responses,
                tree_id=args.tree_id,
                biological_id=metadata.get("response_biological_id", ""),
                technical_id=metadata.get("response_technical_id", ""),
                batch=metadata.get("response_batch", ""),
                within_variance=args.within_variance,
                technical_aggregation=args.technical_aggregation,
                standard_error_columns=metadata.get("standard_error_columns", ""),
                sample_size_columns="",
            )
            if not response_estimates.tip_summary.empty:
                response_summary_frames.append(response_estimates.tip_summary.assign(aggregation=aggregation))
            if response_estimates.sampling_covariance_by_trait:
                response_covariance_frames.append(
                    _sampling_table(response_estimates.sampling_covariance_by_trait, leaf_names, aggregation)
                )

            predictor_estimates_by_analysis: dict[str, Any] = {}
            for plan_row in plan.to_dict("records"):
                analysis_id = str(plan_row["analysis_id"])
                predictors = _csv(plan_row["predictors"])
                categorical = _csv(plan_row.get("categorical_predictors", ""))
                ordered = _ordered(str(plan_row.get("ordered_predictors", ".")))
                predictor_estimates = _read_trait_inputs(
                    species_traits_path,
                    tree,
                    predictors,
                    tree_id=args.tree_id,
                    biological_id=args.predictor_biological_id,
                    technical_id=args.predictor_technical_id,
                    batch=args.predictor_batch,
                    within_variance=args.predictor_within_variance,
                    technical_aggregation=args.predictor_technical_aggregation,
                    standard_error_columns=(
                        ""
                        if str(plan_row.get("predictor_standard_error_columns", ".")) == "."
                        else str(plan_row["predictor_standard_error_columns"])
                    ),
                    sample_size_columns=(
                        ""
                        if str(plan_row.get("predictor_sample_size_columns", ".")) == "."
                        else str(plan_row["predictor_sample_size_columns"])
                    ),
                    categorical=categorical,
                    ordered=ordered,
                    categorical_policy=args.categorical_replicate_policy,
                )
                predictor_estimates_by_analysis[analysis_id] = predictor_estimates
                if not predictor_estimates.tip_summary.empty:
                    predictor_summary_frames.append(
                        predictor_estimates.tip_summary.assign(aggregation=aggregation, analysis_id=analysis_id)
                    )
                if predictor_estimates.sampling_covariance_by_trait:
                    predictor_covariance_frames.append(
                        _sampling_table(
                            predictor_estimates.sampling_covariance_by_trait,
                            leaf_names,
                            aggregation,
                        ).assign(analysis_id=analysis_id)
                    )

                if "species-nwkit" not in methods:
                    continue
                try:
                    result = fit_ordinary_pgls(
                        tree,
                        response_estimates.values_by_trait,
                        predictor_estimates.values_by_trait,
                        responses,
                        predictors,
                        response_sampling_covariance=response_estimates.sampling_covariance_by_trait,
                        predictor_sampling_covariance=predictor_estimates.sampling_covariance_by_trait,
                        evolution_model=args.response_evolution_model,
                        evolution_parameter=_parameter(args.response_evolution_parameter),
                        branch_length=args.branch_length,
                        intercept=True,
                        confidence_level=args.confidence_level,
                        inference=args.inference,
                        bootstrap_replicates=args.bootstrap_replicates,
                        seed=args.seed,
                        reml=args.reml == "yes",
                        predictor_evolution_model=args.predictor_evolution_model,
                        predictor_evolution_parameter=_parameter(args.predictor_evolution_parameter),
                        predictor_branch_length=args.predictor_branch_length,
                        categorical_predictors=categorical,
                        ordered_predictors=ordered,
                        factor_references=_references(str(plan_row.get("factor_reference", "."))),
                        factor_coding=args.factor_coding,
                        allow_large_dense=args.allow_large_dense == "yes",
                    )
                    result.insert(0, "analysis_id", analysis_id)
                    result["tree_id"] = args.tree_id
                    result.insert(0, "estimand", f"species-level {aggregation} paralog expression")
                    result.insert(0, "aggregation", aggregation)
                    result.insert(0, "analysis_method", "species_nwkit")
                    native_frames.append(result)
                    native_status = _native_status(
                        args.tree_id,
                        aggregation,
                        analysis_id,
                        responses,
                        result,
                        nwkit_version,
                    )
                    status_rows.extend(native_status)
                    audit_status = "ok" if any(row["status"] == "ok" for row in native_status) else "not_estimable"
                    audit_records.append(
                        {
                            "analysis_method": "species_nwkit",
                            "aggregation": aggregation,
                            "analysis_id": analysis_id,
                            "status": audit_status,
                            "reason": ";".join(
                                str(row["reason"]) for row in native_status if row["status"] != "ok" and row["reason"]
                            ),
                            "nwkit_version": nwkit_version,
                        }
                    )
                except ValueError as exc:
                    reason = " ".join(str(exc).split())
                    status_rows.extend(
                        _native_status(
                            args.tree_id,
                            aggregation,
                            analysis_id,
                            responses,
                            None,
                            nwkit_version,
                            reason=reason,
                        )
                    )
                    audit_records.append(
                        {
                            "analysis_method": "species_nwkit",
                            "aggregation": aggregation,
                            "analysis_id": analysis_id,
                            "status": "not_estimable",
                            "error": {"type": "ValueError", "message": reason},
                            "nwkit_version": nwkit_version,
                        }
                    )
            if "species-rphylopars" in methods:
                rphylopars_summary_frames.append(
                    _summary_long(
                        aggregation,
                        leaf_names,
                        response_estimates,
                        predictor_estimates_by_analysis,
                    )
                )

    native = (
        pandas.concat(native_frames, ignore_index=True, sort=False)
        if native_frames
        else pandas.DataFrame(
            columns=[
                "analysis_method",
                "aggregation",
                "estimand",
                "analysis_id",
                "response",
                "term",
                "coefficient",
                "standard_error",
                "p_value",
            ]
        )
    )
    rphylopars = pandas.DataFrame(columns=native.columns)
    if "species-rphylopars" in methods:
        rphylopars, r_status, r_audit = _run_rphylopars(
            args,
            tree,
            pandas.concat(rphylopars_summary_frames, ignore_index=True, sort=False),
            plan,
            responses,
        )
        status_rows.extend(r_status.to_dict("records"))
        audit_records.append(r_audit)

    native.to_csv(args.native_out, sep="\t", index=False, na_rep="NA")
    rphylopars.to_csv(args.rphylopars_out, sep="\t", index=False, na_rep="NA")
    _comparison_table(args.rsc_results, native, rphylopars).to_csv(
        args.comparison_out, sep="\t", index=False, na_rep="NA"
    )
    _append_unrequested_status(status_rows, methods, args.tree_id)
    pandas.DataFrame(status_rows, columns=STATUS_COLUMNS).to_csv(args.status_out, sep="\t", index=False, na_rep="NA")
    audit_records.extend(_audit_status_record(row) for row in status_rows)
    args.audit_out.write_text(
        "".join(json.dumps(record, sort_keys=True) + "\n" for record in audit_records),
        encoding="utf-8",
    )

    def write_optional(frames: list[pandas.DataFrame], path: Path, columns: Sequence[str]) -> None:
        frame = pandas.concat(frames, ignore_index=True, sort=False) if frames else pandas.DataFrame(columns=columns)
        frame.to_csv(path, sep="\t", index=False, na_rep="NA")

    write_optional(
        response_summary_frames, args.response_tip_summary_out, ["aggregation", "tree_id", "leaf_name", "trait"]
    )
    write_optional(
        response_covariance_frames,
        args.response_sampling_covariance_out,
        [
            "aggregation",
            "tree_id",
            "trait",
            "leaf_name_1",
            "leaf_name_2",
            "sampling_covariance",
            "covariance_representation",
        ],
    )
    write_optional(
        predictor_summary_frames,
        args.predictor_tip_summary_out,
        ["aggregation", "analysis_id", "tree_id", "leaf_name", "trait"],
    )
    write_optional(
        predictor_covariance_frames,
        args.predictor_sampling_covariance_out,
        [
            "aggregation",
            "analysis_id",
            "tree_id",
            "trait",
            "leaf_name_1",
            "leaf_name_2",
            "sampling_covariance",
            "covariance_representation",
        ],
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--methods", required=True)
    parser.add_argument("--tree-id", required=True)
    parser.add_argument("--species-tree", required=True, type=Path)
    parser.add_argument("--reconciliation", required=True, type=Path)
    parser.add_argument("--expression", required=True, type=Path)
    parser.add_argument("--species-traits", required=True, type=Path)
    parser.add_argument("--analysis-plan", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--aggregation", choices=[*AGGREGATIONS, "all"], default="sum")
    parser.add_argument("--expression-value-type", choices=["identity", "log", "log2", "log2p1"], default="log2p1")
    parser.add_argument("--paralog-missing", choices=["error", "ignore"], default="error")
    parser.add_argument(
        "--paralog-sampling-covariance",
        type=Path,
        help=("Optional TSV of known within-species cross-paralog sampling covariances for standard-error propagation"),
    )
    parser.add_argument("--within-variance", choices=["pooled", "leaf", "known-se"], default="pooled")
    parser.add_argument("--technical-aggregation", choices=["error", "mean"], default="error")
    parser.add_argument("--predictor-biological-id", default="")
    parser.add_argument("--predictor-technical-id", default="")
    parser.add_argument("--predictor-batch", default="")
    parser.add_argument("--predictor-within-variance", choices=["pooled", "leaf", "known-se"], default="pooled")
    parser.add_argument("--predictor-technical-aggregation", choices=["error", "mean"], default="error")
    parser.add_argument("--categorical-replicate-policy", choices=["error", "latent"], default="latent")
    parser.add_argument("--factor-coding", choices=["treatment", "sum"], default="treatment")
    parser.add_argument("--branch-length", choices=["original", "unit"], default="original")
    parser.add_argument("--response-evolution-model", required=True)
    parser.add_argument("--response-evolution-parameter", default="auto")
    parser.add_argument("--predictor-evolution-model", required=True)
    parser.add_argument("--predictor-evolution-parameter", default="auto")
    parser.add_argument("--predictor-branch-length", choices=["original", "unit"], default="original")
    parser.add_argument("--inference", choices=["wald", "parametric-bootstrap"], default="wald")
    parser.add_argument("--bootstrap-replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--confidence-level", type=float, default=0.95)
    parser.add_argument("--reml", choices=["yes", "no"], default="yes")
    parser.add_argument("--allow-large-dense", choices=["yes", "no"], default="no")
    parser.add_argument(
        "--rphylopars-sampling-covariance", choices=["require-diagonal", "diagonalize"], default="require-diagonal"
    )
    parser.add_argument("--rphylopars-script", required=True, type=Path)
    parser.add_argument("--rsc-results", type=Path)
    parser.add_argument("--rsc-status", type=Path)
    parser.add_argument("--empty-reason", default="")
    parser.add_argument("--native-out", required=True, type=Path)
    parser.add_argument("--rphylopars-out", required=True, type=Path)
    parser.add_argument("--comparison-out", required=True, type=Path)
    parser.add_argument("--status-out", required=True, type=Path)
    parser.add_argument("--audit-out", required=True, type=Path)
    parser.add_argument("--expression-summary-out", required=True, type=Path)
    parser.add_argument("--expression-audit-out", required=True, type=Path)
    parser.add_argument("--response-tip-summary-out", required=True, type=Path)
    parser.add_argument("--response-sampling-covariance-out", required=True, type=Path)
    parser.add_argument("--predictor-tip-summary-out", required=True, type=Path)
    parser.add_argument("--predictor-sampling-covariance-out", required=True, type=Path)
    parser.set_defaults(func=run)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.bootstrap_replicates < 2:
        parser.error("--bootstrap-replicates must be at least 2")
    if args.seed < 0:
        parser.error("--seed must be non-negative")
    if not 0.0 < args.confidence_level < 1.0:
        parser.error("--confidence-level must lie strictly between zero and one")
    try:
        return int(args.func(args))
    except (OSError, RuntimeError, ValueError, pandas.errors.ParserError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
