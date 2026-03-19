#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
import sys

import format_species_inputs as fsi


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Run a single gg_input_generation species-formatting task from a task-plan JSON."
    )
    parser.add_argument("--task-plan", required=True, help="Task-plan JSON created by plan_input_generation_tasks.py.")
    parser.add_argument("--task-index", type=int, required=True, help="1-based task index to execute.")
    parser.add_argument(
        "--species-cds-dir",
        required=True,
        help="Output directory for formatted species CDS FASTA files.",
    )
    parser.add_argument(
        "--species-gff-dir",
        required=True,
        help="Output directory for formatted species GFF files.",
    )
    parser.add_argument(
        "--species-genome-dir",
        required=True,
        help="Output directory for formatted species genome FASTA files.",
    )
    parser.add_argument(
        "--species-summary-output",
        default="",
        help="Optional per-task species summary TSV path.",
    )
    parser.add_argument(
        "--stats-output",
        default="",
        help="Optional JSON path for task-level formatting stats.",
    )
    parser.add_argument(
        "--task-meta-output",
        default="",
        help="Optional JSON path for task metadata such as species_prefix and formatted output paths.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing formatted outputs.")
    parser.add_argument("--dry-run", action="store_true", help="Plan actions without writing formatted outputs.")
    return parser


def load_task_plan(path):
    with open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def deserialize_task(raw_task):
    task = dict(raw_task)
    for path_key in ("cds_path", "gff_path", "gbff_path", "genome_path"):
        value = str(task.get(path_key) or "").strip()
        task[path_key] = Path(value) if value != "" else None
    return task


def write_json(path_text, payload):
    if path_text == "":
        return
    path = Path(path_text).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wt", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=True, indent=2, sort_keys=True)


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    task_plan_path = Path(args.task_plan).expanduser().resolve()
    if not task_plan_path.exists():
        parser.error("Task plan not found: {}".format(task_plan_path))

    task_plan = load_task_plan(task_plan_path)
    tasks = task_plan.get("tasks") or []
    if args.task_index < 1 or args.task_index > len(tasks):
        parser.error(
            "--task-index {} is out of range for {} tasks".format(args.task_index, len(tasks))
        )

    task = deserialize_task(tasks[args.task_index - 1])
    output_cds_dir = Path(args.species_cds_dir).expanduser().resolve()
    output_gff_dir = Path(args.species_gff_dir).expanduser().resolve()
    output_genome_dir = Path(args.species_genome_dir).expanduser().resolve()
    output_cds_dir.mkdir(parents=True, exist_ok=True)
    output_gff_dir.mkdir(parents=True, exist_ok=True)
    output_genome_dir.mkdir(parents=True, exist_ok=True)

    cds_result = fsi.format_cds(task, output_cds_dir, args.overwrite, args.dry_run)
    gff_result = fsi.format_gff(task, output_gff_dir, args.overwrite, args.dry_run)
    genome_result = fsi.format_genome(task, output_genome_dir, args.overwrite, args.dry_run)

    run_started_utc = fsi.utc_now_iso()
    if fsi.format_task_succeeded(cds_result, gff_result, genome_result, args.dry_run):
        taxonomy_resolver = fsi.SpeciesTaxonomyMetadataResolver.from_environment()
        species_summary_row = fsi.build_species_summary_row(
            task,
            cds_result,
            gff_result,
            genome_result,
            run_started_utc=run_started_utc,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
            taxonomy_metadata=taxonomy_resolver.resolve(task["species_prefix"]),
        )
        if args.species_summary_output != "":
            row_key = fsi.species_row_key(
                task["provider"],
                task["species_key"],
                task["species_prefix"],
            )
            fsi.write_species_summary_rows(
                Path(args.species_summary_output).expanduser().resolve(),
                {row_key: species_summary_row},
            )

    stats = {
        "species_processed": 1,
        "num_species_cds_files": 1,
        "num_species_gff_files": 1,
        "num_species_genome_files": 0 if genome_result["status"] == "missing" else 1,
        "cds_sequences_before": cds_result.get("before_count", 0),
        "cds_sequences_after": cds_result.get("after_count", 0),
        "cds_first_sequence_name": cds_result.get("first_sequence_name", ""),
        "aggregated_cds_removed": cds_result.get("duplicates", 0),
        "duplicate_cds_ids_skipped": cds_result.get("duplicates", 0),
        "dry_run": int(args.dry_run),
        "species_prefix": task["species_prefix"],
        "species_key": task["species_key"],
        "provider": task["provider"],
    }
    write_json(args.stats_output, stats)
    write_json(
        args.task_meta_output,
        {
            "task_index": args.task_index,
            "task_count": len(tasks),
            "species_prefix": task["species_prefix"],
            "species_key": task["species_key"],
            "provider": task["provider"],
            "cds_output_path": str(cds_result.get("output_path") or ""),
            "gff_output_path": str(gff_result.get("output_path") or ""),
            "genome_output_path": str(genome_result.get("output_path") or ""),
            "cds_status": cds_result.get("status", ""),
            "gff_status": gff_result.get("status", ""),
            "genome_status": genome_result.get("status", ""),
        },
    )

    print(
        "[{}] {}: CDS={} ({}), GFF={} ({}), GENOME={} ({})".format(
            task["provider"],
            task["species_prefix"],
            task["cds_path"].name if task.get("cds_path") is not None else "DERIVED_FROM_GFF_GENOME",
            cds_result["status"],
            task["gff_path"].name,
            gff_result["status"],
            task["genome_path"].name if task.get("genome_path") is not None else "NA",
            genome_result["status"],
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
