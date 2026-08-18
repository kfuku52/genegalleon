#!/usr/bin/env python3
"""Aggregate per-species kfFractBias self-synteny retention outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import statistics
from collections import defaultdict
from pathlib import Path

FIELDNAMES = (
    "analysis_id",
    "species",
    "depth",
    "window_size",
    "denominator",
    "target_gene_count",
    "synteny_pair_count",
    "interchromosomal_pair_count",
    "intrachromosomal_pair_count",
    "analyzed_window_count",
    "mean_best_interchromosomal_retention_percent",
    "median_best_interchromosomal_retention_percent",
    "maximum_best_interchromosomal_retention_percent",
    "mean_intrachromosomal_retention_percent",
)


def _read_self_rows(table_path: Path) -> list[dict[str, str]]:
    with table_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Fractionation-bias table has no header: {table_path}")
        required = {"analysis_id", "target_species", "query_species", "quota"}
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise ValueError(f"Fractionation-bias table is missing required columns: {', '.join(missing)}")
        rows = []
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"Fractionation-bias table row {line_number} has too many fields")
            values = {key: (value or "").strip() for key, value in row.items()}
            if not any(values.values()):
                continue
            if values.get("mode", "compare") != "self":
                continue
            if values["target_species"] != values["query_species"]:
                raise ValueError(f"Self row {line_number} has different target and query species")
            quota_parts = values["quota"].split(":")
            if (
                len(quota_parts) != 2
                or quota_parts[0] != quota_parts[1]
                or not quota_parts[0].isdigit()
                or int(quota_parts[0]) < 1
            ):
                raise ValueError(f"Self row {line_number} does not have a symmetric positive quota")
            rows.append(values)
    if not rows:
        raise ValueError(f"No rows with mode=self were found in {table_path}")
    return rows


def _format_statistic(values: list[float], statistic) -> str:
    if not values:
        return ""
    return f"{statistic(values):.10g}"


def _window_statistics(path: Path) -> tuple[int, str, str, str, str]:
    best_interchromosomal: dict[tuple[str, str, str, str], list[float]] = defaultdict(list)
    intrachromosomal: list[float] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "target_seqid",
            "query_seqid",
            "window_index",
            "start_rank",
            "end_rank",
            "retention_percent",
        }
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(f"Invalid kfFractBias windows table: {path}")
        for row in reader:
            try:
                retention = float(row["retention_percent"])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Invalid retention_percent in {path}") from exc
            if not math.isfinite(retention) or not 0 <= retention <= 100:
                raise ValueError(f"Out-of-range retention_percent in {path}: {retention}")
            if row["target_seqid"] == row["query_seqid"]:
                intrachromosomal.append(retention)
                continue
            key = (row["target_seqid"], row["window_index"], row["start_rank"], row["end_rank"])
            best_interchromosomal[key].append(retention)

    best_values = [max(values) for values in best_interchromosomal.values()]
    return (
        len(best_values),
        _format_statistic(best_values, statistics.fmean),
        _format_statistic(best_values, statistics.median),
        _format_statistic(best_values, max),
        _format_statistic(intrachromosomal, statistics.fmean),
    )


def summarize(table_path: Path, results_dir: Path) -> list[dict[str, str | int]]:
    output_rows: list[dict[str, str | int]] = []
    seen_analysis_ids: set[str] = set()
    for config in _read_self_rows(table_path):
        analysis_id = config["analysis_id"]
        if not analysis_id:
            raise ValueError("A self row has an empty analysis_id")
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", analysis_id):
            raise ValueError(f"Unsafe self analysis_id: {analysis_id!r}")
        if analysis_id in seen_analysis_ids:
            raise ValueError(f"Duplicate self analysis_id: {analysis_id}")
        seen_analysis_ids.add(analysis_id)
        result_dir = results_dir / analysis_id
        summary_path = result_dir / f"{analysis_id}.summary.json"
        windows_path = result_dir / f"{analysis_id}.windows.tsv"
        if not summary_path.is_file() or not windows_path.is_file():
            raise ValueError(f"Self analysis outputs are incomplete for {analysis_id}: {result_dir}")
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if summary.get("analysis_mode") != "self_synteny_retention":
            raise ValueError(f"Unexpected analysis_mode in {summary_path}: {summary.get('analysis_mode')!r}")
        counts = summary.get("counts", {})
        parameters = summary.get("parameters", {})
        window_count, mean_best, median_best, maximum_best, mean_intrachromosomal = _window_statistics(windows_path)
        output_rows.append(
            {
                "analysis_id": analysis_id,
                "species": config["target_species"],
                "depth": config["quota"].split(":", 1)[0],
                "window_size": parameters.get("window_size", ""),
                "denominator": parameters.get("denominator", ""),
                "target_gene_count": counts.get("target_gene_count", ""),
                "synteny_pair_count": counts.get("synteny_pair_count", ""),
                "interchromosomal_pair_count": counts.get("interchromosomal_pair_count", ""),
                "intrachromosomal_pair_count": counts.get("intrachromosomal_pair_count", ""),
                "analyzed_window_count": window_count,
                "mean_best_interchromosomal_retention_percent": mean_best,
                "median_best_interchromosomal_retention_percent": median_best,
                "maximum_best_interchromosomal_retention_percent": maximum_best,
                "mean_intrachromosomal_retention_percent": mean_intrachromosomal,
            }
        )
    return output_rows


def _write_tsv(path: Path, rows: list[dict[str, str | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def _plot(path_pdf: Path, path_png: Path, rows: list[dict[str, str | int]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    labels = [str(row["analysis_id"]) for row in rows]
    values = [
        float(row["mean_best_interchromosomal_retention_percent"])
        if row["mean_best_interchromosomal_retention_percent"] != ""
        else 0.0
        for row in rows
    ]
    figure_height = max(3.5, min(24.0, 0.42 * len(rows) + 1.8))
    figure, axis = plt.subplots(figsize=(8.5, figure_height))
    positions = list(range(len(rows)))
    axis.barh(positions, values, color="#4472C4")
    axis.set_yticks(positions, labels=labels)
    axis.invert_yaxis()
    axis.set_xlim(0, 100)
    axis.set_xlabel("Mean best interchromosomal retention per target window (%)")
    axis.set_title("kfFractBias self-synteny retention summary")
    axis.grid(axis="x", alpha=0.25)
    figure.tight_layout()
    path_pdf.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path_pdf)
    figure.savefig(path_png, dpi=180)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", required=True, type=Path)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--output-tsv", required=True, type=Path)
    parser.add_argument("--output-pdf", required=True, type=Path)
    parser.add_argument("--output-png", required=True, type=Path)
    args = parser.parse_args()
    rows = summarize(args.table, args.results_dir)
    _write_tsv(args.output_tsv, rows)
    _plot(args.output_pdf, args.output_png, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
