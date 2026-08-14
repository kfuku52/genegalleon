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
        return combined, None

    errors = numpy.asarray(standard_errors, dtype=float)
    if errors.shape != linear.shape or not numpy.isfinite(errors).all() or (errors < 0.0).any():
        raise ValueError("Paralog standard errors must be finite, non-negative, and aligned")
    if aggregation == "max":
        maximum = int(numpy.argmax(linear))
        return combined, float(errors[maximum])

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
    variance = float(numpy.sum(numpy.square(derivative * errors)))
    return combined, math.sqrt(max(0.0, variance))


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


def aggregate_species_expression(
    expression: pandas.DataFrame,
    reconciliation: pandas.DataFrame,
    responses: Sequence[str],
    aggregations: Sequence[str],
    *,
    value_type: str,
    missing_policy: str,
    tree_id: str,
) -> tuple[dict[str, pandas.DataFrame], pandas.DataFrame]:
    if "leaf_name" not in expression:
        raise ValueError("Prepared expression requires a leaf_name column")
    mapping = _gene_species_map(reconciliation)
    unknown = sorted(set(expression["leaf_name"].astype(str)) - set(mapping))
    if unknown:
        raise ValueError("Prepared expression contains unmapped gene tips: " + ", ".join(unknown[:5]))
    working = expression.copy()
    working.insert(1, "species", working["leaf_name"].astype(str).map(mapping))
    expected = pandas.Series(list(mapping.values()), dtype=str).value_counts().to_dict()
    metadata_columns = [column for column in ("biological_id", "technical_id", "batch") if column in working]
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
                value, standard_error = _aggregate_values(
                    values,
                    errors,
                    aggregation,
                    value_type,
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
    replicate_requested = bool(biological_id) or within_variance == "known-se"
    if replicate_requested:
        return _read_mixed_replicate_traits(
            arguments,
            tree,
            list(traits),
            set(categorical),
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
        ),
        sampling_covariance_by_trait={},
        tip_summary=pandas.DataFrame(columns=TIP_SUMMARY_COLUMNS),
        model_by_trait={trait: "single-observation" for trait in traits},
    )


def _covariance_diagonal_and_offdiagonal(covariance: Any, leaf_names: Sequence[str]) -> tuple[numpy.ndarray, bool]:
    from nwkit.gaussian import DiagonalLowRankCovariance

    if isinstance(covariance, DiagonalLowRankCovariance):
        low_rank = numpy.asarray(covariance.low_rank, dtype=float)
        diagonal = numpy.asarray(covariance.diagonal, dtype=float) + numpy.sum(low_rank * low_rank, axis=1)
        gram = low_rank @ low_rank.T
        numpy.fill_diagonal(gram, 0.0)
        return diagonal, bool(numpy.any(numpy.abs(gram) > 1e-12))
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
        usable = selected
        if not usable.empty and "inference_status" in usable:
            usable = usable.loc[usable["inference_status"].astype(str) == "ok"]
        response_reason = reason
        if usable.empty and not response_reason:
            if selected.empty:
                response_reason = "nwkit_returned_no_response_rows"
            elif "inference_status" in selected:
                states = sorted(set(selected["inference_status"].dropna().astype(str)))
                response_reason = "nwkit_inference_status:" + ",".join(states or ["unknown"])
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
    summary_long: pandas.DataFrame,
    plan: pandas.DataFrame,
    responses: Sequence[str],
) -> tuple[pandas.DataFrame, pandas.DataFrame, dict[str, object]]:
    with tempfile.TemporaryDirectory(prefix="genegalleon-rphylopars-") as temporary:
        directory = Path(temporary)
        summary_path = directory / "summary.tsv"
        plan_path = directory / "plan.tsv"
        results_path = directory / "results.tsv"
        status_path = directory / "status.tsv"
        summary_long.to_csv(summary_path, sep="\t", index=False, na_rep="NA")
        plan.to_csv(plan_path, sep="\t", index=False, na_rep="NA")
        command = [
            "Rscript",
            str(args.rphylopars_script),
            f"--tree={args.species_tree}",
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
        if "inference_status" in selected:
            selected = selected.loc[selected["inference_status"].astype(str) == "ok"]
        if selected.empty or "p_value" not in selected:
            continue
        selected["_p_value"] = pandas.to_numeric(selected["p_value"], errors="coerce")
        selected = selected.loc[numpy.isfinite(selected["_p_value"])]
        if selected.empty:
            continue
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
    plan = _read_tsv(args.analysis_plan)
    if plan.empty:
        raise ValueError("Prepared RSC analysis plan is empty")
    tree = read_tree(str(args.species_tree), "auto", True)
    leaf_names = [str(name) for name in tree.leaf_names()]
    aggregated, aggregation_audit = aggregate_species_expression(
        expression,
        reconciliation,
        responses,
        aggregations,
        value_type=args.expression_value_type,
        missing_policy=args.paralog_missing,
        tree_id=args.tree_id,
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
            frame.to_csv(expression_path, sep="\t", index=False, na_rep="NA")
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
                predictor_estimates = _read_trait_inputs(
                    args.species_traits,
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
                        ordered_predictors=_ordered(str(plan_row.get("ordered_predictors", "."))),
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
