#!/usr/bin/env python3
"""GeneGalleon adapters and summaries for NWKIT reconciled PGLS.

The gene-evolution workflow stores expression in a gene-by-sample wide table,
whereas NWKIT consumes a leaf table and explicit replicate identifiers.  This
module performs that lossless format adaptation, combines one or more NWKIT
bundles, and provides compact metrics for ``stat_tree``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Iterable, Sequence

import pandas

MISSING_VALUES = {"", "NA", "NaN", "nan", "?", "missing", "unknown", "."}

PGLS_BUNDLE_SUFFIXES = (
    ".reconciliation.tsv",
    ".gene-contrasts.tsv",
    ".species-contrasts.tsv",
    ".response-sampling-covariance.tsv",
    ".response-tip-summary.tsv",
    ".predictor-sampling-covariance.tsv",
    ".predictor-tip-summary.tsv",
    ".random-effects.tsv",
    ".sensitivity.tsv",
    ".trait-origins.tsv",
    ".pgls.tsv",
)

REQUIRED_NWKIT_SUFFIXES = {
    ".reconciliation.tsv",
    ".gene-contrasts.tsv",
    ".species-contrasts.tsv",
    ".random-effects.tsv",
    ".pgls.tsv",
}

EMPTY_BUNDLE_COLUMNS = {
    ".reconciliation.tsv": [
        "tree_id",
        "gene_branch_id",
        "gene_clade_id",
        "event_type",
        "event_source",
        "species_event_id",
        "mapping_status",
        "eligible",
        "coverage_status",
        "reason",
    ],
    ".gene-contrasts.tsv": [
        "tree_id",
        "gene_clade_id",
        "lineage_clade_id",
        "event_type",
        "eligible",
        "coverage_status",
        "species_event_id",
        "trait",
        "evolution_model",
        "evolution_parameter_name",
        "evolution_parameter",
        "branch_length_mode",
        "raw_contrast",
        "contrast_variance",
        "standardized_contrast",
    ],
    ".species-contrasts.tsv": [
        "tree_id",
        "branch_clade_id",
        "descendant_taxa",
        "trait",
        "evolution_model",
        "evolution_parameter_name",
        "evolution_parameter",
        "branch_length_mode",
        "raw_contrast",
        "contrast_variance",
        "standardized_contrast",
    ],
    ".response-sampling-covariance.tsv": [
        "tree_id",
        "trait",
        "contrast_id_1",
        "contrast_id_2",
        "sampling_covariance",
        "covariance_representation",
    ],
    ".predictor-sampling-covariance.tsv": [
        "tree_id",
        "trait",
        "contrast_id_1",
        "contrast_id_2",
        "sampling_covariance",
        "covariance_representation",
    ],
    ".response-tip-summary.tsv": [
        "tree_id",
        "leaf_name",
        "trait",
        "n_biological",
        "n_technical",
        "mean",
        "within_sd",
        "standard_error",
        "variance_method",
        "batch_adjusted",
    ],
    ".predictor-tip-summary.tsv": [
        "tree_id",
        "leaf_name",
        "trait",
        "n_biological",
        "n_technical",
        "mean",
        "within_sd",
        "standard_error",
        "variance_method",
        "batch_adjusted",
    ],
    ".random-effects.tsv": [
        "analysis_id",
        "model_id",
        "tree_id",
        "response",
        "effect_type",
        "group_id",
        "term",
        "conditional_mode",
        "conditional_standard_error",
        "inference_status",
    ],
    ".sensitivity.tsv": [
        "analysis_id",
        "analysis_type",
        "model_id",
        "tree_id",
        "response",
        "group_id",
        "term",
        "full_coefficient",
        "omitted_coefficient",
        "coefficient_change",
        "inference_status",
        "message",
    ],
    ".trait-origins.tsv": [
        "analysis_id",
        "trait",
        "origin_id",
        "branch_id",
        "branch_clade_id",
        "from_state",
        "to_state",
        "posterior_frequency",
        "credible",
    ],
    ".pgls.tsv": [
        "analysis_id",
        "model_id",
        "tree_id",
        "response",
        "term",
        "source_term",
        "coefficient",
        "standard_error",
        "statistic",
        "degrees_of_freedom",
        "p_value",
        "confidence_interval_lower",
        "confidence_interval_upper",
        "n_gene_contrasts",
        "n_species_events",
        "inference_status",
    ],
}

STAT_TREE_RESULT_COLUMNS = (
    "coefficient",
    "standard_error",
    "statistic",
    "degrees_of_freedom",
    "p_value",
    "confidence_interval_lower",
    "confidence_interval_upper",
    "n_gene_contrasts",
    "n_species_events",
    "n_repeated_gene_contrasts",
    "n_lineages",
    "n_excluded_ineligible",
    "n_excluded_coverage",
    "response_evolution_model",
    "response_evolution_parameter_name",
    "response_evolution_parameter",
    "response_evolution_parameter_status",
    "response_evolution_optimizer_converged",
    "response_evolution_optimizer_message",
    "response_evolution_boundary_warning",
    "response_evolution_parameter_bootstrap_refit",
    "response_branch_length_mode",
    "predictor_evolution_model",
    "predictor_evolution_parameter_name",
    "predictor_evolution_parameter",
    "predictor_evolution_parameter_status",
    "predictor_evolution_optimizer_converged",
    "predictor_evolution_optimizer_message",
    "predictor_evolution_boundary_warning",
    "predictor_evolution_log_likelihood",
    "predictor_branch_length_mode",
    "coverage_policy",
    "small_sample_warning",
    "inference_status",
    "model",
    "inference_method",
    "reml",
    "optimizer_converged",
    "boundary_warning",
    "evolutionary_rate",
    "species_event_variance",
    "lineage_slope_variance",
    "mean_sampling_variance",
    "sampling_variance_fraction",
    "mean_predictor_sampling_variance",
    "mean_latent_predictor_variance",
    "predictor_uncertainty_fraction",
    "predictor_evolutionary_rate",
    "predictor_rate_optimizer_converged",
    "predictor_rate_optimizer_message",
    "predictor_rate_boundary_warning",
    "measurement_error_model",
    "log_likelihood",
    "event_random_effect",
    "lineage_random_slope",
    "ensemble_size",
    "between_tree_variance",
    "tree_support_fraction",
)


def _read_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"), None)
    if not header:
        raise ValueError(f"TSV input has no header: {path}")
    duplicates = sorted({name for name in header if header.count(name) > 1})
    if duplicates:
        raise ValueError(f"TSV input has duplicate column names ({', '.join(duplicates)}): {path}")
    if any(name == "" for name in header):
        raise ValueError(f"TSV input has an empty column name: {path}")
    return header


def read_tsv(path: Path, *, allow_empty: bool = False) -> pandas.DataFrame:
    _read_header(path)
    frame = pandas.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        na_filter=False,
    )
    if frame.empty and not allow_empty:
        raise ValueError(f"TSV input has no data rows: {path}")
    return frame


def _validate_generated_trait_names(names: Sequence[str], label: str) -> None:
    reserved = {"leaf_name", "biological_id", "technical_id", "batch"}
    conflicts = sorted(set(names) & reserved)
    generated = {
        generated_name for name in names for generated_name in (f"{name}__standard_error", f"{name}__sample_size")
    }
    generated_conflicts = sorted(set(names) & generated)
    if conflicts or generated_conflicts:
        invalid = sorted(set(conflicts + generated_conflicts))
        raise ValueError(f"{label} names conflict with reserved/generated NWKIT columns: " + ", ".join(invalid))


def _parse_csv_names(value: str, available: Sequence[str], option: str) -> list[str]:
    if str(value).strip().lower() == "all":
        selected = list(available)
    else:
        selected = [item.strip() for item in str(value).split(",")]
        if any(item == "" for item in selected):
            raise ValueError(f"{option} contains an empty column name")
    duplicates = sorted({name for name in selected if selected.count(name) > 1})
    if duplicates:
        raise ValueError(f"{option} contains duplicates: {', '.join(duplicates)}")
    missing = [name for name in selected if name not in available]
    if missing:
        raise ValueError(f"{option} selects columns not present in the input: {', '.join(missing)}")
    if not selected:
        raise ValueError(f"{option} did not select any columns")
    return selected


def _is_missing(value: object) -> bool:
    return str(value).strip() in MISSING_VALUES


def _metadata_value(row: dict[str, object], column: str) -> str:
    value = str(row.get(column, "")).strip()
    return "" if _is_missing(value) else value


def _numeric_nonmissing(values: Iterable[object], label: str) -> list[float]:
    numeric: list[float] = []
    invalid: list[str] = []
    for value in values:
        if _is_missing(value):
            continue
        try:
            converted = float(str(value).strip())
        except ValueError:
            invalid.append(str(value))
            continue
        if not math.isfinite(converted):
            invalid.append(str(value))
        else:
            numeric.append(converted)
    if invalid:
        examples = ", ".join(repr(value) for value in invalid[:3])
        raise ValueError(f"Continuous response {label!r} contains nonnumeric values: {examples}")
    return numeric


def _response_base(column: str, separator: str) -> tuple[str, str]:
    if separator and separator in column:
        base, suffix = column.rsplit(separator, 1)
        if base and suffix.isdigit():
            return base, suffix
    return column, "1"


def _slug(value: str, maximum: int = 64) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_.-")
    if not slug:
        slug = "item"
    if len(slug) <= maximum:
        return slug
    digest = hashlib.sha1(str(value).encode("utf-8")).hexdigest()[:10]
    return f"{slug[: maximum - 11]}_{digest}"


def _parse_keyed_setting(value: str, option: str) -> dict[str, str]:
    if str(value).strip() in {"", "."}:
        return {}
    out: dict[str, str] = {}
    for item in str(value).split(","):
        if "=" not in item:
            raise ValueError(f"{option} must use TRAIT=VALUE entries")
        key, setting = item.split("=", 1)
        key = key.strip()
        setting = setting.strip()
        if not key or not setting:
            raise ValueError(f"{option} contains an empty trait or value")
        if key in out:
            raise ValueError(f"{option} contains duplicate trait {key!r}")
        out[key] = setting
    return out


def _write_metadata(path: Path, values: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value"])
        for key, value in values.items():
            writer.writerow([key, value])


def _prepare_expression(args: argparse.Namespace) -> tuple[pandas.DataFrame, dict[str, object]]:
    expression = read_tsv(args.expression, allow_empty=True)
    leaf_column = expression.columns[0]
    expression[leaf_column] = expression[leaf_column].astype(str).str.strip()
    if (expression[leaf_column] == "").any():
        raise ValueError("Expression input contains an empty gene/leaf identifier")
    if expression[leaf_column].duplicated().any():
        duplicated = expression.loc[expression[leaf_column].duplicated(), leaf_column].iloc[0]
        raise ValueError(f"Expression input contains duplicate gene/leaf identifier: {duplicated}")

    source_columns = list(expression.columns[1:])
    if not source_columns:
        raise ValueError("Expression input contains no response columns")

    sample_metadata = None
    standard_error_by_column: dict[str, str] = {}
    sample_size_by_column: dict[str, str] = {}
    if args.sample_metadata:
        sample_metadata = read_tsv(args.sample_metadata)
        required = {"column", "response"}
        if args.within_variance != "known-se":
            required.add("biological_id")
        missing = sorted(required - set(sample_metadata.columns))
        if missing:
            raise ValueError("Expression sample metadata is missing required columns: " + ", ".join(missing))
        if sample_metadata["column"].duplicated().any():
            raise ValueError("Expression sample metadata contains duplicate 'column' entries")
        unknown = sorted(set(sample_metadata["column"]) - set(source_columns))
        if unknown:
            raise ValueError("Expression sample metadata refers to absent expression columns: " + ", ".join(unknown))
        mapping_rows = sample_metadata.to_dict("records")
        available_responses = list(dict.fromkeys(str(row["response"]).strip() for row in mapping_rows))
        if any(not response for response in available_responses):
            raise ValueError("Expression sample metadata contains an empty response")
    else:
        mapping_rows = []
        available_responses = []
        for source_column in source_columns:
            response, replicate = _response_base(source_column, args.replicate_separator)
            available_responses.append(response)
            mapping_rows.append(
                {
                    "column": source_column,
                    "response": response,
                    # Do not assume that replicate 1 of different responses is paired.
                    "biological_id": f"auto:{response}:{replicate}",
                    "technical_id": "",
                    "batch": "",
                    "standard_error_column": "",
                    "sample_size_column": "",
                }
            )
        available_responses = list(dict.fromkeys(available_responses))

    responses = _parse_csv_names(args.responses, available_responses, "--responses")
    _validate_generated_trait_names(responses, "Response")
    selected_rows = [row for row in mapping_rows if str(row["response"]).strip() in responses]

    if expression.empty:
        return pandas.DataFrame(columns=["leaf_name"]), {
            "status": "not_estimable",
            "reason": "expression_has_no_rows",
            "responses": "",
            "response_biological_id": "",
            "response_technical_id": "",
            "response_batch": "",
            "standard_error_columns": "",
            "sample_size_columns": "",
            "response_sampling_uncertainty": "no",
        }

    skipped: list[str] = []
    retained: list[str] = []
    for response in responses:
        response_columns = [str(row["column"]) for row in selected_rows if str(row["response"]).strip() == response]
        observed_by_leaf = expression[response_columns].apply(
            lambda row: any(not _is_missing(value) for value in row), axis=1
        )
        numeric = _numeric_nonmissing(expression[response_columns].to_numpy().ravel(), response)
        if not observed_by_leaf.all():
            skipped.append(f"{response}:missing_leaf_values")
        elif len(set(numeric)) < 2:
            skipped.append(f"{response}:constant_or_empty")
        else:
            retained.append(response)
    responses = retained
    selected_rows = [row for row in selected_rows if str(row["response"]).strip() in responses]

    if not responses:
        empty = pandas.DataFrame(columns=["leaf_name"])
        return empty, {
            "status": "not_estimable",
            "reason": ";".join(skipped) or "no_usable_responses",
            "responses": "",
            "response_biological_id": "",
            "response_technical_id": "",
            "response_batch": "",
            "standard_error_columns": "",
            "sample_size_columns": "",
            "response_sampling_uncertainty": "no",
        }

    if args.within_variance == "known-se":
        counts_by_response = {
            response: sum(str(row["response"]).strip() == response for row in selected_rows) for response in responses
        }
        repeated = [response for response, count in counts_by_response.items() if count != 1]
        if repeated:
            raise ValueError(
                "rsc_within_variance=known-se requires exactly one summarized "
                "measurement column per response: " + ", ".join(repeated)
            )
        if any(_metadata_value(row, role) for row in selected_rows for role in ("technical_id", "batch")):
            raise ValueError(
                "Known-SE expression input cannot also specify technical_id or batch; "
                "provide one summarized mean and standard error per gene/response"
            )

        se_columns: dict[str, str] = {}
        n_columns: dict[str, str] = {}
        output = expression[[leaf_column]].rename(columns={leaf_column: "leaf_name"}).copy()
        for mapping in selected_rows:
            response = str(mapping["response"]).strip()
            source = str(mapping["column"])
            se_source = _metadata_value(mapping, "standard_error_column")
            n_source = _metadata_value(mapping, "sample_size_column")
            if not se_source:
                raise ValueError(
                    f"rsc_within_variance=known-se requires standard_error_column for expression measurement {source!r}"
                )
            if se_source not in expression.columns:
                raise ValueError(f"Standard-error source column is absent: {se_source}")
            output[response] = expression[source]
            se_columns[response] = f"{response}__standard_error"
            output[se_columns[response]] = expression[se_source]
            if n_source:
                if n_source not in expression.columns:
                    raise ValueError(f"Sample-size source column is absent: {n_source}")
                n_columns[response] = f"{response}__sample_size"
                output[n_columns[response]] = expression[n_source]
        if n_columns and len(n_columns) != len(responses):
            raise ValueError("sample_size_column must be supplied for every selected response when it is used")
        return output, {
            "status": "ready",
            "reason": ";".join(skipped),
            "responses": ",".join(responses),
            "response_biological_id": "",
            "response_technical_id": "",
            "response_batch": "",
            "standard_error_columns": ",".join(se_columns[response] for response in responses),
            "sample_size_columns": ",".join(n_columns[response] for response in responses if response in n_columns),
            "response_sampling_uncertainty": "yes",
        }

    has_technical = any(_metadata_value(row, "technical_id") for row in selected_rows)
    has_batch = any(_metadata_value(row, "batch") for row in selected_rows)
    if has_technical and any(not _metadata_value(row, "technical_id") for row in selected_rows):
        raise ValueError("Expression sample metadata mixes empty and nonempty technical_id values")
    if has_batch and any(not _metadata_value(row, "batch") for row in selected_rows):
        raise ValueError("Expression sample metadata mixes empty and nonempty batch values")

    biological_counts: dict[str, int] = {}
    for response in responses:
        response_rows = [row for row in selected_rows if str(row["response"]).strip() == response]
        biological_ids = {_metadata_value(row, "biological_id") for row in response_rows}
        if "" in biological_ids:
            source = next(str(row["column"]) for row in response_rows if not _metadata_value(row, "biological_id"))
            raise ValueError(f"Empty biological_id for expression column {source!r}")
        observation_ids = {
            (
                _metadata_value(row, "biological_id"),
                _metadata_value(row, "technical_id") if has_technical else "",
            )
            for row in response_rows
        }
        if len(observation_ids) != len(response_rows):
            raise ValueError(
                "Expression sample metadata maps multiple columns to the same "
                f"observation/response for {response!r}; use distinct "
                "biological_id/technical_id values"
            )
        biological_counts[response] = len(biological_ids)
    replicate_mode = has_technical or has_batch or any(count > 1 for count in biological_counts.values())

    if not replicate_mode:
        if args.within_variance != "pooled":
            raise ValueError(f"rsc_within_variance={args.within_variance} requires biological expression replicates")
        rename = {leaf_column: "leaf_name"}
        for row in selected_rows:
            rename[str(row["column"])] = str(row["response"]).strip()
        output = expression[[leaf_column] + [str(row["column"]) for row in selected_rows]].rename(columns=rename)
        return output, {
            "status": "ready",
            "reason": ";".join(skipped),
            "responses": ",".join(responses),
            "response_biological_id": "",
            "response_technical_id": "",
            "response_batch": "",
            "standard_error_columns": "",
            "sample_size_columns": "",
            "response_sampling_uncertainty": "no",
        }

    se_columns: dict[str, str] = {}
    n_columns: dict[str, str] = {}
    for row in selected_rows:
        response = str(row["response"]).strip()
        source = str(row["column"])
        se_source = _metadata_value(row, "standard_error_column")
        n_source = _metadata_value(row, "sample_size_column")
        if se_source:
            if se_source not in expression.columns:
                raise ValueError(f"Standard-error source column is absent: {se_source}")
            standard_error_by_column[source] = se_source
            se_columns[response] = f"{response}__standard_error"
        if n_source:
            if n_source not in expression.columns:
                raise ValueError(f"Sample-size source column is absent: {n_source}")
            sample_size_by_column[source] = n_source
            n_columns[response] = f"{response}__sample_size"

    if standard_error_by_column or sample_size_by_column:
        raise ValueError(
            "Expression standard_error_column/sample_size_column metadata requires rsc_within_variance=known-se"
        )

    output_rows: list[dict[str, object]] = []
    output_columns = ["leaf_name"] + responses + ["biological_id"]
    if has_technical:
        output_columns.append("technical_id")
    if has_batch:
        output_columns.append("batch")
    output_columns.extend(se_columns[response] for response in responses if response in se_columns)
    output_columns.extend(n_columns[response] for response in responses if response in n_columns)

    observation_groups: dict[tuple[str, str], list[dict[str, object]]] = {}
    for mapping in selected_rows:
        biological = _metadata_value(mapping, "biological_id")
        technical = _metadata_value(mapping, "technical_id") if has_technical else ""
        key = (biological, technical)
        response = str(mapping["response"]).strip()
        existing = observation_groups.setdefault(key, [])
        if any(str(item["response"]).strip() == response for item in existing):
            raise ValueError(
                "Expression sample metadata maps multiple columns to the same "
                f"observation/response ({biological!r}, {technical!r}, {response!r}); "
                "use distinct biological_id/technical_id values"
            )
        if has_batch and existing:
            existing_batch = _metadata_value(existing[0], "batch")
            current_batch = _metadata_value(mapping, "batch")
            if existing_batch != current_batch:
                raise ValueError(
                    f"One biological/technical observation cannot span batches: {biological!r}, {technical!r}"
                )
        existing.append(mapping)

    for _, expression_row in expression.iterrows():
        leaf_name = str(expression_row[leaf_column]).strip()
        for (biological, technical), mappings in observation_groups.items():
            row: dict[str, object] = {column: "NA" for column in output_columns}
            row["leaf_name"] = leaf_name
            row["biological_id"] = biological
            if has_technical:
                row["technical_id"] = technical
            if has_batch:
                row["batch"] = _metadata_value(mappings[0], "batch")
            observed = False
            for mapping in mappings:
                source = str(mapping["column"])
                value = expression_row[source]
                if _is_missing(value):
                    continue
                observed = True
                response = str(mapping["response"]).strip()
                row[response] = value
            if observed:
                output_rows.append(row)
    output = pandas.DataFrame(output_rows, columns=output_columns)
    return output, {
        "status": "ready",
        "reason": ";".join(skipped),
        "responses": ",".join(responses),
        "response_biological_id": "biological_id",
        "response_technical_id": "technical_id" if has_technical else "",
        "response_batch": "batch" if has_batch else "",
        "standard_error_columns": ",".join(se_columns[response] for response in responses if response in se_columns),
        "sample_size_columns": ",".join(n_columns[response] for response in responses if response in n_columns),
        "response_sampling_uncertainty": "yes",
    }


def _prepare_predictors(args: argparse.Namespace) -> tuple[pandas.DataFrame, dict[str, object]]:
    traits = read_tsv(args.species_traits)
    leaf_column = traits.columns[0]
    traits[leaf_column] = traits[leaf_column].astype(str).str.strip()
    if (traits[leaf_column] == "").any():
        raise ValueError("Species-trait input contains an empty species identifier")

    identifier_columns = [
        value
        for value in (
            args.predictor_biological_id,
            args.predictor_technical_id,
            args.predictor_batch,
        )
        if value
    ]
    if len(identifier_columns) != len(set(identifier_columns)):
        raise ValueError("Predictor biological-ID, technical-ID, and batch columns must be distinct")
    if leaf_column in identifier_columns:
        raise ValueError("Predictor identifier/batch columns cannot also be the species column")
    auxiliary_columns = [
        item.strip()
        for value in (args.predictor_standard_error_columns, args.predictor_sample_size_columns)
        for item in str(value).split(",")
        if item.strip()
    ]
    if len(auxiliary_columns) != len(set(auxiliary_columns)):
        raise ValueError("Predictor standard-error/sample-size columns must be distinct")
    role_conflicts = sorted(set(identifier_columns) & set(auxiliary_columns))
    if role_conflicts:
        raise ValueError(
            "Predictor identifier/batch and standard-error/sample-size roles "
            "must use distinct columns: " + ", ".join(role_conflicts)
        )
    metadata_columns = set(identifier_columns + auxiliary_columns)
    for column in metadata_columns:
        if column not in traits.columns:
            raise ValueError(f"Predictor metadata column is absent from species traits: {column}")
    for column in identifier_columns:
        missing_ids = traits[column].map(_is_missing)
        if missing_ids.any():
            raise ValueError(f"Predictor identifier/batch column {column!r} contains missing values")
    if (args.predictor_technical_id or args.predictor_batch) and not args.predictor_biological_id:
        raise ValueError("Predictor technical-ID and batch columns require a biological-ID column")
    if args.predictor_biological_id:
        observation_key = [leaf_column, args.predictor_biological_id]
        repeated_observations = traits.duplicated(subset=observation_key, keep=False)
        if repeated_observations.any() and not args.predictor_technical_id:
            raise ValueError(
                "Repeated species/biological predictor observations require a technical-ID column or prior aggregation"
            )
        if (
            args.predictor_technical_id
            and traits.duplicated(subset=[*observation_key, args.predictor_technical_id], keep=False).any()
        ):
            raise ValueError("Predictor technical IDs must be unique within each species/biological observation")
    available = [column for column in traits.columns[1:] if column not in metadata_columns]
    predictors = _parse_csv_names(args.predictors, available, "--predictors")
    _validate_generated_trait_names(predictors, "Predictor")

    standard_error_columns = [
        item.strip() for item in str(args.predictor_standard_error_columns).split(",") if item.strip()
    ]
    sample_size_columns = [item.strip() for item in str(args.predictor_sample_size_columns).split(",") if item.strip()]
    if sample_size_columns and not standard_error_columns:
        raise ValueError("--predictor-sample-size-columns requires --predictor-standard-error-columns")
    if standard_error_columns and len(standard_error_columns) != len(predictors):
        raise ValueError("--predictor-standard-error-columns must contain one column for each selected predictor")
    if sample_size_columns and len(sample_size_columns) != len(predictors):
        raise ValueError("--predictor-sample-size-columns must contain one column for each selected predictor")
    predictor_standard_errors = (
        dict(zip(predictors, standard_error_columns, strict=True)) if standard_error_columns else {}
    )
    predictor_sample_sizes = dict(zip(predictors, sample_size_columns, strict=True)) if sample_size_columns else {}

    if args.predictor_within_variance == "known-se":
        if not standard_error_columns:
            raise ValueError("--predictor-within-variance known-se requires --predictor-standard-error-columns")
        if any((args.predictor_biological_id, args.predictor_technical_id, args.predictor_batch)):
            raise ValueError(
                "Known-SE predictor input requires one summarized row per species and "
                "cannot use predictor replicate IDs or batch"
            )
    elif standard_error_columns or sample_size_columns:
        raise ValueError("Predictor standard-error/sample-size columns require --predictor-within-variance known-se")
    elif args.predictor_within_variance != "pooled" and not args.predictor_biological_id:
        raise ValueError(
            f"--predictor-within-variance {args.predictor_within_variance} requires --predictor-biological-id"
        )

    if traits[leaf_column].duplicated().any() and not args.predictor_biological_id:
        raise ValueError(
            "Species-trait input has repeated species rows; set "
            "rsc_predictor_biological_id to identify biological replicates"
        )

    explicit_categorical = set(item.strip() for item in str(args.categorical_predictors).split(",") if item.strip())
    unknown_categorical = sorted(explicit_categorical - set(predictors))
    if unknown_categorical:
        raise ValueError("Categorical predictors were not selected: " + ", ".join(unknown_categorical))
    ordered = _parse_keyed_setting(args.ordered_predictors, "--ordered-predictors")
    references = _parse_keyed_setting(args.factor_reference, "--factor-reference")
    categorical_and_ordered = sorted(explicit_categorical & set(ordered))
    if categorical_and_ordered:
        raise ValueError("Predictors cannot be both unordered and ordered: " + ", ".join(categorical_and_ordered))
    referenced_ordered = sorted(set(references) & set(ordered))
    if referenced_ordered:
        raise ValueError("--factor-reference applies only to unordered predictors: " + ", ".join(referenced_ordered))
    for label, setting in (("ordered", ordered), ("factor-reference", references)):
        unknown = sorted(set(setting) - set(predictors))
        if unknown:
            raise ValueError(f"{label} predictors were not selected: {', '.join(unknown)}")

    retained: list[str] = []
    skipped: list[str] = []
    inferred_categorical: set[str] = set(explicit_categorical) | set(references)
    for predictor in predictors:
        values = [str(value).strip() for value in traits[predictor] if not _is_missing(value)]
        if len(set(values)) < 2:
            skipped.append(f"{predictor}:constant_or_empty")
            continue
        if len(values) != traits.shape[0]:
            skipped.append(f"{predictor}:missing_species_values")
            continue
        if predictor not in ordered:
            numeric = pandas.to_numeric(pandas.Series(values), errors="coerce")
            if numeric.isna().any():
                inferred_categorical.add(predictor)
        retained.append(predictor)

    retained_categorical = set(retained) & (inferred_categorical | set(ordered))
    if retained_categorical and args.predictor_within_variance == "known-se":
        raise ValueError(
            "Known standard errors do not apply to categorical predictors: " + ", ".join(sorted(retained_categorical))
        )
    if retained_categorical and args.predictor_batch:
        raise ValueError(
            "Batch adjustment is not supported for categorical predictors: " + ", ".join(sorted(retained_categorical))
        )

    columns = [leaf_column] + retained + sorted(metadata_columns, key=list(traits.columns).index)
    output = traits[columns].rename(columns={leaf_column: "leaf_name"})
    return output, {
        "status": "ready" if retained else "not_estimable",
        "reason": ";".join(skipped) if skipped else ("" if retained else "no_usable_predictors"),
        "predictors": ",".join(retained),
        "categorical_predictors": ",".join(predictor for predictor in retained if predictor in inferred_categorical),
        "ordered_predictors": ordered,
        "factor_reference": references,
        "standard_error_columns": predictor_standard_errors,
        "sample_size_columns": predictor_sample_sizes,
        "predictor_sampling_uncertainty": (
            "yes" if args.predictor_biological_id or args.predictor_within_variance == "known-se" else "no"
        ),
    }


def prepare(args: argparse.Namespace) -> int:
    expression, expression_meta = _prepare_expression(args)
    predictors, predictor_meta = _prepare_predictors(args)

    args.expression_output.parent.mkdir(parents=True, exist_ok=True)
    args.species_traits_output.parent.mkdir(parents=True, exist_ok=True)
    expression.to_csv(args.expression_output, sep="\t", index=False, na_rep="NA")
    predictors.to_csv(args.species_traits_output, sep="\t", index=False, na_rep="NA")

    status = "ready"
    reasons = [value for value in (expression_meta["reason"], predictor_meta["reason"]) if value]
    if expression_meta["status"] != "ready" or predictor_meta["status"] != "ready":
        status = "not_estimable"

    predictor_names = [name for name in str(predictor_meta["predictors"]).split(",") if name]
    categorical = set(name for name in str(predictor_meta["categorical_predictors"]).split(",") if name)
    ordered = predictor_meta["ordered_predictors"]
    references = predictor_meta["factor_reference"]
    predictor_standard_errors = predictor_meta["standard_error_columns"]
    predictor_sample_sizes = predictor_meta["sample_size_columns"]
    plan_rows: list[dict[str, str]] = []
    predictor_sets = [predictor_names] if args.predictor_mode == "joint" else [[name] for name in predictor_names]
    for index, selected in enumerate(predictor_sets, start=1):
        if not selected:
            continue
        analysis_id = "joint" if args.predictor_mode == "joint" else f"p{index:03d}_{_slug(selected[0])}"
        selected_set = set(selected)
        ordered_setting = ",".join(f"{name}={value}" for name, value in ordered.items() if name in selected_set)
        reference_setting = ",".join(f"{name}={value}" for name, value in references.items() if name in selected_set)
        selected_categorical = [name for name in selected if name in categorical]
        plan_rows.append(
            {
                "analysis_id": analysis_id,
                "predictors": ",".join(selected),
                "categorical_predictors": ",".join(selected_categorical) or ".",
                "ordered_predictors": ordered_setting or ".",
                "factor_reference": reference_setting or ".",
                "predictor_standard_error_columns": ",".join(
                    predictor_standard_errors[name] for name in selected if name in predictor_standard_errors
                )
                or ".",
                "predictor_sample_size_columns": ",".join(
                    predictor_sample_sizes[name] for name in selected if name in predictor_sample_sizes
                )
                or ".",
                "origin_applicable": "yes" if selected_categorical else "no",
            }
        )
    if status == "ready" and not plan_rows:
        status = "not_estimable"
        reasons.append("no_analysis_plan")

    args.analysis_plan_output.parent.mkdir(parents=True, exist_ok=True)
    pandas.DataFrame(
        plan_rows,
        columns=[
            "analysis_id",
            "predictors",
            "categorical_predictors",
            "ordered_predictors",
            "factor_reference",
            "predictor_standard_error_columns",
            "predictor_sample_size_columns",
            "origin_applicable",
        ],
    ).to_csv(args.analysis_plan_output, sep="\t", index=False)

    _write_metadata(
        args.metadata_output,
        {
            "status": status,
            "reason": ";".join(reasons),
            "responses": expression_meta["responses"],
            "predictors": predictor_meta["predictors"],
            "analysis_count": len(plan_rows),
            "response_biological_id": expression_meta["response_biological_id"],
            "response_technical_id": expression_meta["response_technical_id"],
            "response_batch": expression_meta["response_batch"],
            "standard_error_columns": expression_meta["standard_error_columns"],
            "sample_size_columns": expression_meta["sample_size_columns"],
            "response_sampling_uncertainty": expression_meta["response_sampling_uncertainty"],
            "predictor_sampling_uncertainty": predictor_meta["predictor_sampling_uncertainty"],
        },
    )
    return 0


def has_nhx(args: argparse.Namespace) -> int:
    text = args.tree.read_text(encoding="utf-8", errors="replace")
    return 0 if re.search(r"\[&&NHX:[^\]]*\bD=[YN01]", text) else 1


def inspect_reconciliation(args: argparse.Namespace) -> int:
    frame = read_tsv(args.reconciliation)
    required = {"event_type", "eligible", "coverage_status", "species_event_id"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError("Reconciliation table is missing columns: " + ", ".join(missing))
    selected = frame[(frame["event_type"] == "speciation") & (frame["eligible"] == "yes")]
    if args.coverage == "complete":
        selected = selected[selected["coverage_status"] == "complete"]
    event_ids = {
        str(value).strip() for value in selected["species_event_id"] if str(value).strip() not in MISSING_VALUES
    }
    status = "ready" if len(event_ids) >= args.min_species_events else "not_estimable"
    reason = (
        "" if status == "ready" else (f"eligible_species_events={len(event_ids)}<minimum={args.min_species_events}")
    )
    _write_metadata(
        args.output,
        {
            "status": status,
            "reason": reason,
            "eligible_species_events": len(event_ids),
        },
    )
    return 0


def inspect_audit_error(args: argparse.Namespace) -> int:
    records = []
    for line_number, line in enumerate(args.audit.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        try:
            records.append(json.loads(line))
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid NWKIT audit JSON at {args.audit}:{line_number}") from exc
    if not records:
        raise ValueError(f"NWKIT audit contains no records: {args.audit}")
    record = records[-1]
    error = record.get("error") or {}
    error_type = str(error.get("type", ""))
    message = " ".join(str(error.get("message", "")).split())
    if record.get("status") != "error" or not error_type:
        raise ValueError("NWKIT audit does not describe a failed command")
    _write_metadata(
        args.output,
        {
            "status": "not_estimable" if error_type == "ValueError" else "fatal",
            "error_type": error_type,
            "reason": message or error_type,
        },
    )
    return 0


def _read_bundle_list(path: Path) -> list[dict[str, str]]:
    frame = read_tsv(path)
    required = {"analysis_id", "prefix", "audit"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError("Bundle list is missing columns: " + ", ".join(missing))
    if frame["analysis_id"].duplicated().any():
        raise ValueError("Bundle list contains duplicate analysis_id values")
    return frame.to_dict("records")


def _runtime_empty_bundle_columns() -> dict[str, list[str]]:
    """Use the installed NWKIT schemas while retaining dependency-light adapters."""
    try:
        from nwkit.contrast import (
            BASE_CONTRAST_COLUMNS,
            RECONCILIATION_CONTEXT_COLUMNS,
            REPLICATE_CONTRAST_COLUMNS,
            SAMPLING_COVARIANCE_COLUMNS,
        )
        from nwkit.pgls import RANDOM_EFFECT_COLUMNS, RESULT_COLUMNS, SENSITIVITY_COLUMNS
        from nwkit.reconcile import RECONCILIATION_COLUMNS
        from nwkit.replicates import TIP_SUMMARY_COLUMNS
        from nwkit.rsc_diagnostics import ORIGIN_DIAGNOSTIC_COLUMNS
    except ImportError:
        return EMPTY_BUNDLE_COLUMNS

    return {
        ".reconciliation.tsv": list(RECONCILIATION_COLUMNS),
        ".gene-contrasts.tsv": (
            list(BASE_CONTRAST_COLUMNS[:4])
            + list(RECONCILIATION_CONTEXT_COLUMNS)
            + list(BASE_CONTRAST_COLUMNS[4:])
            + list(REPLICATE_CONTRAST_COLUMNS)
        ),
        ".species-contrasts.tsv": list(BASE_CONTRAST_COLUMNS),
        ".response-sampling-covariance.tsv": list(SAMPLING_COVARIANCE_COLUMNS),
        ".response-tip-summary.tsv": list(TIP_SUMMARY_COLUMNS),
        ".predictor-sampling-covariance.tsv": list(SAMPLING_COVARIANCE_COLUMNS),
        ".predictor-tip-summary.tsv": list(TIP_SUMMARY_COLUMNS),
        ".random-effects.tsv": ["analysis_id", *RANDOM_EFFECT_COLUMNS],
        ".sensitivity.tsv": ["analysis_id", *SENSITIVITY_COLUMNS],
        ".trait-origins.tsv": ["analysis_id", *ORIGIN_DIAGNOSTIC_COLUMNS],
        ".pgls.tsv": ["analysis_id", *RESULT_COLUMNS],
    }


def _empty_frame(suffix: str) -> pandas.DataFrame:
    columns = _runtime_empty_bundle_columns()[suffix]
    if "analysis_id" not in columns:
        columns = ["analysis_id", *columns]
    return pandas.DataFrame(columns=columns)


def _read_bundle_frame(path: Path, suffix: str, required: bool) -> pandas.DataFrame | None:
    if not path.is_file():
        if required:
            raise FileNotFoundError(f"NWKIT bundle member is missing: {path}")
        return None
    _read_header(path)
    return pandas.read_csv(path, sep="\t", low_memory=False)


def _usable_result_mask(results: pandas.DataFrame) -> pandas.Series:
    mask = pandas.Series(True, index=results.index, dtype=bool)
    if "inference_status" in results:
        mask &= results["inference_status"].astype(str).str.lower().eq("ok")
    false_values = {"no", "false", "0"}
    if "optimizer_converged" in results:
        mask &= ~results["optimizer_converged"].astype(str).str.lower().isin(false_values)
    for prefix in ("response_evolution", "predictor_evolution"):
        status_column = f"{prefix}_parameter_status"
        convergence_column = f"{prefix}_optimizer_converged"
        if status_column not in results or convergence_column not in results:
            continue
        estimated = results[status_column].astype(str).str.lower().eq("estimated")
        failed = results[convergence_column].astype(str).str.lower().isin(false_values)
        mask &= ~(estimated & failed)
    return mask


def _write_bundle_frames(
    output_prefix: Path,
    frames_by_suffix: dict[str, list[pandas.DataFrame]],
) -> dict[str, pandas.DataFrame]:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    combined: dict[str, pandas.DataFrame] = {}
    for suffix in PGLS_BUNDLE_SUFFIXES:
        frames = frames_by_suffix.get(suffix, [])
        if frames:
            frame = pandas.concat(frames, ignore_index=True, sort=False)
        else:
            frame = _empty_frame(suffix)
        frame.to_csv(f"{output_prefix}{suffix}", sep="\t", index=False, na_rep="NA")
        combined[suffix] = frame
    return combined


def _status_from_results(
    tree_id: str,
    results: pandas.DataFrame,
    analysis_count: int,
    reason: str = "",
) -> pandas.DataFrame:
    def add_reason(existing: str, value: str) -> str:
        return ";".join(item for item in (existing, value) if item)

    if results.empty:
        status = "not_estimable"
        reason = add_reason(reason, "nwkit_returned_no_model_rows")
        model_count = 0
        estimable_count = 0
    else:
        model_keys = [column for column in ("analysis_id", "model_id") if column in results]
        model_count = int(results[model_keys].drop_duplicates().shape[0]) if model_keys else len(results)
        estimable = results[_usable_result_mask(results)]
        estimable_count = int(estimable[model_keys].drop_duplicates().shape[0]) if model_keys else len(estimable)
        if estimable_count > 0:
            status = "ok"
        else:
            status = "not_estimable"
            reason = add_reason(reason, "nwkit_returned_no_estimable_models")

    usable_results = results[_usable_result_mask(results)]
    p_values = pandas.to_numeric(usable_results.get("p_value", pandas.Series(dtype=float)), errors="coerce")
    finite = p_values[p_values.notna() & pandas.Series([math.isfinite(v) for v in p_values], index=p_values.index)]
    best = usable_results.loc[finite.idxmin()] if not finite.empty else None
    row = {
        "tree_id": tree_id,
        "status": status,
        "reason": reason,
        "n_analyses": analysis_count,
        "n_models": model_count,
        "n_estimable_models": estimable_count,
        "n_result_rows": int(results.shape[0]),
        "n_usable_result_rows": int(usable_results.shape[0]),
        "n_nonconverged_result_rows": int(
            results.shape[0]
            - results[_usable_result_mask(results)].shape[0]
            - (
                0
                if "inference_status" not in results
                else (~results["inference_status"].astype(str).str.lower().eq("ok")).sum()
            )
        ),
        "min_p_value": "" if finite.empty else float(finite.min()),
        "best_analysis_id": "" if best is None else best.get("analysis_id", ""),
        "best_response": "" if best is None else best.get("response", ""),
        "best_term": "" if best is None else best.get("term", ""),
        "best_coefficient": "" if best is None else best.get("coefficient", ""),
        "max_n_species_events": (
            ""
            if results.empty or "n_species_events" not in results
            else pandas.to_numeric(results["n_species_events"], errors="coerce").max()
        ),
    }
    return pandas.DataFrame([row])


def aggregate(args: argparse.Namespace) -> int:
    bundle_rows = _read_bundle_list(args.bundle_list)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    output_columns: dict[str, list[str]] = {}
    result_frames: list[pandas.DataFrame] = []
    audit_records: list[str] = []
    analysis_reasons: list[str] = []
    if args.source_audit and args.source_audit.is_file():
        for line_number, line in enumerate(args.source_audit.read_text(encoding="utf-8").splitlines(), start=1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid NWKIT audit JSON at {args.source_audit}:{line_number}") from exc
            record["genegalleon_stage"] = "reconciliation_preflight"
            audit_records.append(json.dumps(record, sort_keys=True))
    for bundle in bundle_rows:
        analysis_id = str(bundle["analysis_id"])
        prefix = Path(str(bundle["prefix"]))
        for suffix in PGLS_BUNDLE_SUFFIXES:
            bundle_path = Path(f"{prefix}{suffix}")
            frame = _read_bundle_frame(bundle_path, suffix, suffix in REQUIRED_NWKIT_SUFFIXES)
            if frame is None:
                continue
            if "analysis_id" in frame.columns:
                frame = frame.drop(columns=["analysis_id"])
            frame.insert(0, "analysis_id", analysis_id)
            columns = list(frame.columns)
            if suffix in output_columns and columns != output_columns[suffix]:
                raise ValueError(f"NWKIT bundle member {suffix} has inconsistent columns across analyses")
            first_frame = suffix not in output_columns
            output_columns.setdefault(suffix, columns)
            frame.to_csv(
                f"{args.output_prefix}{suffix}",
                sep="\t",
                index=False,
                na_rep="NA",
                mode="w" if first_frame else "a",
                header=first_frame,
            )
            if suffix == ".pgls.tsv":
                result_frames.append(frame)

        audit_path = Path(str(bundle["audit"]))
        if not audit_path.is_file():
            raise FileNotFoundError(f"NWKIT audit file is missing: {audit_path}")
        for line_number, line in enumerate(audit_path.read_text(encoding="utf-8").splitlines(), start=1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid NWKIT audit JSON at {audit_path}:{line_number}") from exc
            record["genegalleon_analysis_id"] = analysis_id
            if record.get("status") == "error" and (record.get("error") or {}).get("type") == "ValueError":
                message = " ".join(str((record.get("error") or {}).get("message", "")).split())
                analysis_reasons.append(f"{analysis_id}:{message or 'not_estimable'}")
            audit_records.append(json.dumps(record, sort_keys=True))

    for suffix in PGLS_BUNDLE_SUFFIXES:
        if suffix not in output_columns:
            _empty_frame(suffix).to_csv(f"{args.output_prefix}{suffix}", sep="\t", index=False, na_rep="NA")
    results = (
        pandas.concat(result_frames, ignore_index=True, sort=False) if result_frames else _empty_frame(".pgls.tsv")
    )
    args.audit_output.parent.mkdir(parents=True, exist_ok=True)
    args.audit_output.write_text("\n".join(audit_records) + ("\n" if audit_records else ""), encoding="utf-8")
    status = _status_from_results(
        args.tree_id,
        results,
        len(bundle_rows),
        reason=";".join(value for value in (args.reason, *analysis_reasons) if value),
    )
    args.status_output.parent.mkdir(parents=True, exist_ok=True)
    status.to_csv(args.status_output, sep="\t", index=False, na_rep="NA")
    return 0


def empty_bundle(args: argparse.Namespace) -> int:
    frames = {suffix: [] for suffix in PGLS_BUNDLE_SUFFIXES}
    if args.reconciliation and args.reconciliation.is_file():
        reconciliation = pandas.read_csv(args.reconciliation, sep="\t", low_memory=False)
        if "analysis_id" in reconciliation.columns:
            reconciliation = reconciliation.drop(columns=["analysis_id"])
        reconciliation.insert(0, "analysis_id", "preflight")
        frames[".reconciliation.tsv"].append(reconciliation)
    _write_bundle_frames(args.output_prefix, frames)
    status = _status_from_results(
        args.tree_id,
        _empty_frame(".pgls.tsv"),
        0,
        reason=args.reason,
    )
    args.status_output.parent.mkdir(parents=True, exist_ok=True)
    status.to_csv(args.status_output, sep="\t", index=False, na_rep="NA")
    args.audit_output.parent.mkdir(parents=True, exist_ok=True)
    if args.source_audit and args.source_audit.is_file():
        args.audit_output.write_text(args.source_audit.read_text(encoding="utf-8"), encoding="utf-8")
    else:
        record = {
            "command": "genegalleon-rsc-preflight",
            "tree_id": args.tree_id,
            "status": "not_estimable",
            "reason": args.reason,
        }
        args.audit_output.write_text(json.dumps(record, sort_keys=True) + "\n", encoding="utf-8")
    return 0


def _sanitize_stat_component(value: object, maximum: int = 80) -> str:
    text = str(value or "").strip()
    normalized = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_").lower() or "none"
    if len(normalized) <= maximum:
        return normalized
    digest = hashlib.sha1(text.encode("utf-8")).hexdigest()[:10]
    return f"{normalized[: maximum - 11]}_{digest}"


def summarize_for_stat_tree(pgls_path: str | Path, status_path: str | Path) -> dict[str, object]:
    """Return stable, one-row RSC metrics for orthogroup ``stat_tree``."""
    out: dict[str, object] = {}
    status_file = Path(status_path)
    if status_file.is_file() and status_file.stat().st_size > 0:
        status = pandas.read_csv(status_file, sep="\t", nrows=1, low_memory=False)
        if not status.empty:
            for column, value in status.iloc[0].items():
                if column == "tree_id":
                    continue
                if pandas.isna(value):
                    continue
                out[f"rsc_{_sanitize_stat_component(column)}"] = value

    result_file = Path(pgls_path)
    if not result_file.is_file() or result_file.stat().st_size == 0:
        return out
    results = pandas.read_csv(result_file, sep="\t", low_memory=False)
    if results.empty:
        return out
    out["rsc_num_result_rows"] = int(results.shape[0])
    usable = results[_usable_result_mask(results)]
    out["rsc_num_usable_result_rows"] = int(usable.shape[0])
    out["rsc_num_nonusable_result_rows"] = int(results.shape[0] - usable.shape[0])
    if usable.empty or "p_value" not in usable:
        return out
    p_values = pandas.to_numeric(usable["p_value"], errors="coerce")
    finite = p_values.notna() & p_values.map(math.isfinite)
    if not finite.any():
        return out
    best = usable.loc[p_values[finite].idxmin()]
    identity_columns = (
        "analysis_id",
        "model_id",
        "response",
        "term",
        "source_term",
        "response_level",
        "predictor_level",
        "term_test",
    )
    for column in (*identity_columns, *STAT_TREE_RESULT_COLUMNS):
        if column not in results.columns:
            continue
        value = best[column]
        if pandas.isna(value):
            continue
        out[f"rsc_best_{_sanitize_stat_component(column)}"] = value
    return out


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--expression", required=True, type=Path)
    prepare_parser.add_argument("--species-traits", required=True, type=Path)
    prepare_parser.add_argument("--responses", default="all")
    prepare_parser.add_argument("--predictors", default="all")
    prepare_parser.add_argument("--predictor-mode", choices=["separate", "joint"], default="separate")
    prepare_parser.add_argument("--replicate-separator", default="_")
    prepare_parser.add_argument("--sample-metadata", type=Path)
    prepare_parser.add_argument("--within-variance", choices=["pooled", "leaf", "known-se"], default="pooled")
    prepare_parser.add_argument("--predictor-biological-id", default="")
    prepare_parser.add_argument("--predictor-technical-id", default="")
    prepare_parser.add_argument("--predictor-batch", default="")
    prepare_parser.add_argument(
        "--predictor-within-variance",
        choices=["pooled", "leaf", "known-se"],
        default="pooled",
    )
    prepare_parser.add_argument("--predictor-standard-error-columns", default="")
    prepare_parser.add_argument("--predictor-sample-size-columns", default="")
    prepare_parser.add_argument("--categorical-predictors", default="")
    prepare_parser.add_argument("--ordered-predictors", default="")
    prepare_parser.add_argument("--factor-reference", default="")
    prepare_parser.add_argument("--expression-output", required=True, type=Path)
    prepare_parser.add_argument("--species-traits-output", required=True, type=Path)
    prepare_parser.add_argument("--analysis-plan-output", required=True, type=Path)
    prepare_parser.add_argument("--metadata-output", required=True, type=Path)
    prepare_parser.set_defaults(func=prepare)

    nhx_parser = subparsers.add_parser("has-nhx")
    nhx_parser.add_argument("--tree", required=True, type=Path)
    nhx_parser.set_defaults(func=has_nhx)

    inspect_parser = subparsers.add_parser("inspect-reconciliation")
    inspect_parser.add_argument("--reconciliation", required=True, type=Path)
    inspect_parser.add_argument("--coverage", choices=["complete", "any"], default="complete")
    inspect_parser.add_argument("--min-species-events", type=int, default=2)
    inspect_parser.add_argument("--output", required=True, type=Path)
    inspect_parser.set_defaults(func=inspect_reconciliation)

    audit_error_parser = subparsers.add_parser("inspect-audit-error")
    audit_error_parser.add_argument("--audit", required=True, type=Path)
    audit_error_parser.add_argument("--output", required=True, type=Path)
    audit_error_parser.set_defaults(func=inspect_audit_error)

    aggregate_parser = subparsers.add_parser("aggregate")
    aggregate_parser.add_argument("--bundle-list", required=True, type=Path)
    aggregate_parser.add_argument("--output-prefix", required=True, type=Path)
    aggregate_parser.add_argument("--status-output", required=True, type=Path)
    aggregate_parser.add_argument("--audit-output", required=True, type=Path)
    aggregate_parser.add_argument("--source-audit", type=Path)
    aggregate_parser.add_argument("--tree-id", required=True)
    aggregate_parser.add_argument("--reason", default="")
    aggregate_parser.set_defaults(func=aggregate)

    empty_parser = subparsers.add_parser("empty-bundle")
    empty_parser.add_argument("--output-prefix", required=True, type=Path)
    empty_parser.add_argument("--status-output", required=True, type=Path)
    empty_parser.add_argument("--audit-output", required=True, type=Path)
    empty_parser.add_argument("--tree-id", required=True)
    empty_parser.add_argument("--reason", required=True)
    empty_parser.add_argument("--reconciliation", type=Path)
    empty_parser.add_argument("--source-audit", type=Path)
    empty_parser.set_defaults(func=empty_bundle)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "min_species_events", 2) < 1:
        parser.error("--min-species-events must be a positive integer")
    try:
        return int(args.func(args))
    except (OSError, ValueError, pandas.errors.ParserError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
