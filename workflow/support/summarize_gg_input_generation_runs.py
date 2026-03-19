#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
import sys


DEFAULT_COLUMNS = [
    "started_utc",
    "exit_code",
    "provider",
    "input_generation_mode",
    "stage_format_status",
    "stage_validate_status",
    "stage_species_busco_status",
    "stage_busco_summary_status",
    "stage_trait_status",
    "num_species_cds",
    "num_species_gff",
    "num_species_genome",
    "num_species_busco_full",
    "num_species_busco_short",
    "num_busco_ids",
    "num_species_trait",
    "num_trait_columns",
    "cds_sequences_before",
    "cds_sequences_after",
    "cds_first_sequence_name",
]


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Summarize gg_input_generation run history from gg_input_generation_runs.tsv."
    )
    parser.add_argument(
        "--infile",
        default="workspace/output/input_generation/gg_input_generation_runs.tsv",
        help="Path to gg_input_generation run summary TSV.",
    )
    parser.add_argument(
        "--last-n",
        type=int,
        default=10,
        help="Number of most recent runs to show in the detailed section.",
    )
    parser.add_argument(
        "--outfile",
        default="-",
        help="Output path. Use '-' for stdout.",
    )
    return parser


def read_rows(path):
    with open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [row for row in reader]
    return rows


def parse_exit_code(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return 999


def sorted_rows(rows):
    return sorted(rows, key=lambda x: x.get("started_utc", ""))


def format_summary(rows, last_n):
    if len(rows) == 0:
        return "No runs found.\n"

    ordered = sorted_rows(rows)
    latest = ordered[-1]
    successful = sum(1 for row in ordered if parse_exit_code(row.get("exit_code", "")) == 0)
    failed = len(ordered) - successful

    provider_counts = {}
    for row in ordered:
        provider = row.get("provider", "")
        provider_counts[provider] = provider_counts.get(provider, 0) + 1

    lines = []
    lines.append("GG Input Generation Run Summary")
    lines.append("total_runs\t{}".format(len(ordered)))
    lines.append("successful_runs\t{}".format(successful))
    lines.append("failed_runs\t{}".format(failed))
    lines.append("latest_started_utc\t{}".format(latest.get("started_utc", "")))
    lines.append("latest_exit_code\t{}".format(latest.get("exit_code", "")))
    lines.append("latest_provider\t{}".format(latest.get("provider", "")))
    lines.append(
        "latest_stage_status\tformat={};validate={};species_busco={};busco_summary={};trait={}".format(
            latest.get("stage_format_status", ""),
            latest.get("stage_validate_status", ""),
            latest.get("stage_species_busco_status", ""),
            latest.get("stage_busco_summary_status", ""),
            latest.get("stage_trait_status", ""),
        )
    )
    lines.append("")
    lines.append("Runs By Provider")
    for provider in sorted(provider_counts.keys()):
        lines.append("{}\t{}".format(provider, provider_counts[provider]))
    lines.append("")

    detail_count = max(0, int(last_n))
    lines.append("Recent Runs (n={})".format(detail_count))
    lines.append("\t".join(DEFAULT_COLUMNS))
    for row in ordered[-detail_count:]:
        values = [row.get(col, "") for col in DEFAULT_COLUMNS]
        lines.append("\t".join(values))
    lines.append("")
    return "\n".join(lines)


def write_output(text, outfile):
    if outfile == "-":
        sys.stdout.write(text)
        return
    outpath = Path(outfile).expanduser().resolve()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "wt", encoding="utf-8") as handle:
        handle.write(text)


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    infile = Path(args.infile).expanduser().resolve()
    if not infile.exists():
        parser.error("Input file not found: {}".format(infile))

    rows = read_rows(infile)
    text = format_summary(rows, args.last_n)
    write_output(text, args.outfile)
    return 0


if __name__ == "__main__":
    sys.exit(main())
