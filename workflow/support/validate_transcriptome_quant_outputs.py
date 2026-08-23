#!/usr/bin/env python3
"""Validate the complete amalgkit quant output set selected by metadata."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Iterable

REQUIRED_ABUNDANCE_COLUMNS = (
    "target_id",
    "length",
    "eff_length",
    "est_counts",
    "tpm",
)
NUMERIC_ABUNDANCE_COLUMNS = ("length", "eff_length", "est_counts", "tpm")
CONTROL_CHARACTER_PATTERN = re.compile(r"[\x00-\x1f\x7f]")


class QuantOutputValidationError(RuntimeError):
    """Raised when the selected quant output set is incomplete or invalid."""


def safe_run_component(value: str) -> str:
    """Apply the same single-component boundary used by amalgkit run IDs."""
    run_id = value.strip()
    if (
        not run_id
        or run_id in {".", ".."}
        or os.path.isabs(run_id)
        or os.path.basename(run_id) != run_id
        or "/" in run_id
        or "\\" in run_id
        or CONTROL_CHARACTER_PATTERN.search(run_id) is not None
    ):
        raise QuantOutputValidationError(f"unsafe metadata run ID: {run_id!r}")
    return run_id


def selected_runs(metadata_path: Path) -> list[str]:
    """Mirror amalgkit quant eligibility for exclusion and is_sampled fields."""
    try:
        with metadata_path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fieldnames = reader.fieldnames or []
            missing = [name for name in ("run", "scientific_name") if name not in fieldnames]
            if missing:
                raise QuantOutputValidationError(
                    "metadata is missing required column(s): " + ", ".join(missing)
                )
            rows = list(reader)
    except (OSError, csv.Error) as exc:
        raise QuantOutputValidationError(f"failed to read metadata: {exc}") from exc

    exclusion_active = "exclusion" in fieldnames and any(
        (row.get("exclusion") or "").strip() for row in rows
    )
    sampled_active = "is_sampled" in fieldnames and any(
        (row.get("is_sampled") or "").strip() for row in rows
    )
    if sampled_active:
        invalid = sorted(
            {
                (row.get("is_sampled") or "").strip().lower()
                for row in rows
                if (row.get("is_sampled") or "").strip().lower()
                not in {"", "yes", "no"}
            }
        )
        if invalid:
            raise QuantOutputValidationError(
                "metadata is_sampled contains invalid value(s): " + ", ".join(invalid)
            )

    selected: list[str] = []
    seen: set[str] = set()
    for row in rows:
        if exclusion_active and (row.get("exclusion") or "").strip().lower() != "no":
            continue
        if sampled_active and (row.get("is_sampled") or "").strip().lower() != "yes":
            continue
        run_id = safe_run_component(row.get("run") or "")
        scientific_name = (row.get("scientific_name") or "").strip()
        if not scientific_name:
            raise QuantOutputValidationError(
                f"metadata scientific_name is empty for run: {run_id}"
            )
        if run_id in seen:
            raise QuantOutputValidationError(f"metadata run ID is duplicated: {run_id}")
        seen.add(run_id)
        selected.append(run_id)
    if not selected:
        raise QuantOutputValidationError("metadata selected no eligible quant runs")
    return selected


def validate_abundance(path: Path, run_id: str) -> None:
    if not path.is_file() or path.is_symlink() or path.stat().st_size <= 0:
        raise QuantOutputValidationError(
            f"missing or empty regular abundance table for run: {run_id}"
        )
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fieldnames = reader.fieldnames or []
            missing = [name for name in REQUIRED_ABUNDANCE_COLUMNS if name not in fieldnames]
            if missing:
                raise QuantOutputValidationError(
                    f"abundance table for {run_id} is missing column(s): "
                    + ", ".join(missing)
                )
            seen_targets: set[str] = set()
            row_count = 0
            for row in reader:
                row_count += 1
                target_id = (row.get("target_id") or "").strip()
                if not target_id:
                    raise QuantOutputValidationError(
                        f"abundance table for {run_id} contains an empty target_id"
                    )
                if target_id in seen_targets:
                    raise QuantOutputValidationError(
                        f"abundance table for {run_id} contains duplicate target_id: {target_id}"
                    )
                seen_targets.add(target_id)
                for column in NUMERIC_ABUNDANCE_COLUMNS:
                    raw_value = (row.get(column) or "").strip()
                    try:
                        value = float(raw_value)
                    except ValueError as exc:
                        raise QuantOutputValidationError(
                            f"abundance table for {run_id} has invalid {column}"
                        ) from exc
                    if not math.isfinite(value) or value < 0:
                        raise QuantOutputValidationError(
                            f"abundance table for {run_id} has invalid {column}"
                        )
            if row_count == 0:
                raise QuantOutputValidationError(
                    f"abundance table for {run_id} contains no data rows"
                )
    except (OSError, csv.Error) as exc:
        raise QuantOutputValidationError(
            f"failed to read abundance table for {run_id}: {exc}"
        ) from exc


def validate_run_info(path: Path, run_id: str) -> None:
    if not path.is_file() or path.is_symlink() or path.stat().st_size <= 0:
        raise QuantOutputValidationError(
            f"missing or empty regular run-info JSON for run: {run_id}"
        )
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise QuantOutputValidationError(
            f"failed to read run-info JSON for {run_id}: {exc}"
        ) from exc
    if not isinstance(payload, dict) or "p_pseudoaligned" not in payload:
        raise QuantOutputValidationError(
            f"run-info JSON for {run_id} is missing p_pseudoaligned"
        )
    try:
        value = float(payload["p_pseudoaligned"])
    except (TypeError, ValueError) as exc:
        raise QuantOutputValidationError(
            f"run-info JSON for {run_id} has invalid p_pseudoaligned"
        ) from exc
    if not math.isfinite(value) or value < 0 or value > 100:
        raise QuantOutputValidationError(
            f"run-info JSON for {run_id} has invalid p_pseudoaligned"
        )


def validate_quant_outputs(metadata: Path, quant_root: Path) -> dict[str, object]:
    runs = selected_runs(metadata)
    if not quant_root.is_dir() or quant_root.is_symlink():
        raise QuantOutputValidationError("quant output root is not a regular directory")
    for run_id in runs:
        run_directory = quant_root / run_id
        if not run_directory.is_dir() or run_directory.is_symlink():
            raise QuantOutputValidationError(
                f"quant output directory is missing for run: {run_id}"
            )
        validate_abundance(run_directory / f"{run_id}_abundance.tsv", run_id)
        validate_run_info(run_directory / f"{run_id}_run_info.json", run_id)
    return {"status": "valid", "run_count": len(runs), "runs": runs}


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--quant-root", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        result = validate_quant_outputs(args.metadata, args.quant_root)
    except QuantOutputValidationError as exc:
        print(f"quant output validation error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(result, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
