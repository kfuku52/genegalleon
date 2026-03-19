#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
import sys

import format_species_inputs as fsi


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Merge gg_input_generation per-task shard outputs into canonical summaries."
    )
    parser.add_argument(
        "--species-summary-shard-dir",
        required=True,
        help="Directory containing per-task species summary TSV shards.",
    )
    parser.add_argument(
        "--species-summary-output",
        required=True,
        help="Canonical merged species summary TSV output path.",
    )
    parser.add_argument(
        "--task-stats-dir",
        required=True,
        help="Directory containing per-task formatting stats JSON files.",
    )
    parser.add_argument(
        "--aggregate-stats-output",
        required=True,
        help="Output JSON path for aggregated task stats.",
    )
    parser.add_argument(
        "--expected-task-count",
        type=int,
        default=0,
        help="Optional expected number of completed task shards/stats.",
    )
    return parser


def read_task_stats(path):
    with open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    species_summary_shard_dir = Path(args.species_summary_shard_dir).expanduser().resolve()
    species_summary_output = Path(args.species_summary_output).expanduser().resolve()
    task_stats_dir = Path(args.task_stats_dir).expanduser().resolve()
    aggregate_stats_output = Path(args.aggregate_stats_output).expanduser().resolve()

    shard_paths = sorted(species_summary_shard_dir.glob("*.tsv")) if species_summary_shard_dir.exists() else []
    stats_paths = sorted(task_stats_dir.glob("*.json")) if task_stats_dir.exists() else []

    if args.expected_task_count > 0 and len(stats_paths) != args.expected_task_count:
        parser.error(
            "Expected {} task stats shards but found {}".format(
                args.expected_task_count, len(stats_paths)
            )
        )

    merged_rows = fsi.retain_existing_species_summary_rows(
        fsi.read_species_summary_rows(species_summary_output)
    )
    shard_row_count = 0
    for shard_path in shard_paths:
        shard_rows = fsi.read_species_summary_rows(shard_path)
        shard_row_count += len(shard_rows)
        merged_rows.update(shard_rows)
    fsi.write_species_summary_rows(species_summary_output, merged_rows)

    species_processed = 0
    num_species_cds_files = 0
    num_species_gff_files = 0
    num_species_genome_files = 0
    cds_sequences_before = 0
    cds_sequences_after = 0
    cds_first_sequence_name = ""
    aggregated_cds_removed = 0

    for stats_path in stats_paths:
        payload = read_task_stats(stats_path)
        species_processed += int(payload.get("species_processed", 0) or 0)
        num_species_cds_files += int(payload.get("num_species_cds_files", 0) or 0)
        num_species_gff_files += int(payload.get("num_species_gff_files", 0) or 0)
        num_species_genome_files += int(payload.get("num_species_genome_files", 0) or 0)
        cds_sequences_before += int(payload.get("cds_sequences_before", 0) or 0)
        cds_sequences_after += int(payload.get("cds_sequences_after", 0) or 0)
        aggregated_cds_removed += int(payload.get("aggregated_cds_removed", 0) or 0)
        if cds_first_sequence_name == "":
            cds_first_sequence_name = str(payload.get("cds_first_sequence_name", "") or "")

    aggregate_stats_output.parent.mkdir(parents=True, exist_ok=True)
    with open(aggregate_stats_output, "wt", encoding="utf-8") as handle:
        json.dump(
            {
                "species_processed": species_processed,
                "num_species_cds_files": num_species_cds_files,
                "num_species_gff_files": num_species_gff_files,
                "num_species_genome_files": num_species_genome_files,
                "cds_sequences_before": cds_sequences_before,
                "cds_sequences_after": cds_sequences_after,
                "cds_first_sequence_name": cds_first_sequence_name,
                "aggregated_cds_removed": aggregated_cds_removed,
                "duplicate_cds_ids_skipped": aggregated_cds_removed,
                "merged_species_summary_rows": len(merged_rows),
                "task_stats_files": len(stats_paths),
                "task_species_summary_shards": shard_row_count,
            },
            handle,
            ensure_ascii=True,
            indent=2,
            sort_keys=True,
        )

    print(
        "Merged {} species summary shards and {} task stats shards.".format(
            len(shard_paths), len(stats_paths)
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
