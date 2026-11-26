#!/usr/bin/env python3
"""Validate the complete amalgkit quant output set selected by metadata."""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import os
import re
import sys
import zipfile
from pathlib import Path
from typing import Iterable, TextIO

from directory_storage_policy import is_transcriptome_quant_species
from safe_zip_extract import validated_members

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


def validate_abundance_stream(handle: TextIO, run_id: str) -> None:
    try:
        reader = csv.reader(handle, delimiter="\t")
        positions = {name: index for index, name in enumerate(next(reader, ()))}
        missing = [name for name in REQUIRED_ABUNDANCE_COLUMNS if name not in positions]
        if missing:
            raise QuantOutputValidationError(
                f"abundance table for {run_id} is missing column(s): "
                + ", ".join(missing)
            )
        target_position = positions["target_id"]
        numeric_positions = [(column, positions[column]) for column in NUMERIC_ABUNDANCE_COLUMNS]
        required_width = max(positions[column] for column in REQUIRED_ABUNDANCE_COLUMNS) + 1
        seen_targets: set[str] = set()
        for row in reader:
            if not row:
                continue
            # Match DictReader's blank defaults for omitted trailing fields.
            if len(row) < required_width:
                row.extend([""] * (required_width - len(row)))
            target_id = row[target_position].strip()
            if not target_id:
                raise QuantOutputValidationError(
                    f"abundance table for {run_id} contains an empty target_id"
                )
            if target_id in seen_targets:
                raise QuantOutputValidationError(
                    f"abundance table for {run_id} contains duplicate target_id: {target_id}"
                )
            seen_targets.add(target_id)
            for column, position in numeric_positions:
                raw_value = row[position].strip()
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
        if not seen_targets:
            raise QuantOutputValidationError(
                f"abundance table for {run_id} contains no data rows"
            )
    except (OSError, csv.Error) as exc:
        raise QuantOutputValidationError(
            f"failed to read abundance table for {run_id}: {exc}"
        ) from exc


def validate_run_info_stream(handle: TextIO, run_id: str) -> None:
    try:
        payload = json.load(handle)
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


def validate_quant_outputs(
    metadata: Path, quant_root: Path, *, require_transcript_map: bool = False,
) -> dict[str, object]:
    runs = selected_runs(metadata)
    archive_path = quant_root.with_name(quant_root.name + ".zip")
    required = [
        (f"{run}/{run}_{suffix}", validator, run)
        for run in runs
        for suffix, validator in (
            ("abundance.tsv", validate_abundance_stream),
            ("run_info.json", validate_run_info_stream),
        )
    ]
    if require_transcript_map:
        required.append(("transcript_gene.tsv", None, ""))
    if is_transcriptome_quant_species(quant_root) and os.path.lexists(archive_path):
        if os.path.lexists(quant_root) or archive_path.is_symlink() or not archive_path.is_file():
            raise QuantOutputValidationError("quantification ZIP has an unsafe or ambiguous path")
        with zipfile.ZipFile(archive_path) as archive:
            infos, _ = validated_members(archive, archive_path, quant_root.name)
            members = {info.filename: info for info in infos}
            for relative, validator, run in required:
                member = members.get(f"{quant_root.name}/{relative}")
                if member is None or member.file_size <= 0:
                    raise QuantOutputValidationError(f"missing or empty quantification file: {relative}")
                if validator is not None:
                    with io.TextIOWrapper(archive.open(member), encoding="utf-8", newline="") as handle:
                        validator(handle, run)
    else:
        if not quant_root.is_dir() or quant_root.is_symlink():
            raise QuantOutputValidationError("quant output root is not a regular directory")
        for relative, validator, run in required:
            path = quant_root / relative
            if path.parent.is_symlink() or not path.parent.is_dir():
                raise QuantOutputValidationError(f"quant output directory is missing for run: {run}")
            if not path.is_file() or path.is_symlink() or path.stat().st_size <= 0:
                raise QuantOutputValidationError(f"missing or empty regular quantification file: {relative}")
            if validator is not None:
                with path.open(encoding="utf-8", newline="") as handle:
                    validator(handle, run)
    return {"status": "valid", "run_count": len(runs), "runs": runs}


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--quant-root", required=True, type=Path)
    parser.add_argument("--require-transcript-map", action="store_true")
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        result = validate_quant_outputs(
            args.metadata, args.quant_root, require_transcript_map=args.require_transcript_map,
        )
    except (RuntimeError, OSError, UnicodeError, zipfile.BadZipFile) as exc:
        print(f"quant output validation error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(result, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
