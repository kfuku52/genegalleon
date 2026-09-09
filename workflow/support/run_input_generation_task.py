#!/usr/bin/env python3

import argparse
import csv
import json
import os
import sys
from pathlib import Path

import format_species_inputs as fsi
from format_species_provider_config import DEFAULT_INPUT_RELATIVE_DIRS
from input_generation_array_state import atomic_json, digest, load_plan


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Run a single gg_input_generation species-formatting task from a task-plan JSON."
    )
    parser.add_argument("--download-timeout", type=float, default=120)
    parser.add_argument("--http-header", action="append", default=[])
    parser.add_argument("--auth-bearer-token-env", default="")
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
    parser.add_argument(
        "--describe-only",
        action="store_true",
        help="Resolve manifest downloads if needed and write metadata without formatting outputs.",
    )
    parser.add_argument(
        "--reuse-existing",
        action="store_true",
        help="Reuse existing formatted outputs without applying input-audit regeneration checks.",
    )
    return parser


def load_task_plan(path):
    return load_plan(path)


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


def resolve_manifest_task(task, args):
    if "manifest_row" not in task:
        return task
    for path, expected in task.get("input_sha256", {}).items():
        if digest(path) != expected:
            raise ValueError("Local manifest input changed after planning: " + path)
    plan_sha256 = digest(args.task_plan)
    root = Path(str(args.task_plan) + ".tasks")
    root.mkdir(parents=True, exist_ok=True)
    resolved = root / (str(args.task_index) + ".json")
    if resolved.exists():
        cached = json.loads(resolved.read_text())
        if cached.get("plan_sha256") != plan_sha256 or cached.get("task_index") != args.task_index:
            raise ValueError("Resolved download cache belongs to another plan/task")
        actual = deserialize_task(cached["task"])
        if any(actual.get(key) != task.get(key) for key in ("species_prefix", "species_key", "provider")):
            raise ValueError("Resolved download cache species/provider mismatch")
        return actual
    if load_plan(args.task_plan).get("download_mode") == "staged":
        raise ValueError("Staged download receipt is missing; rerun array_prepare before workers")
    manifest = root / (str(args.task_index) + ".tsv")
    row = task["manifest_row"]
    with open(manifest, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    download_root = Path(task["download_dir"]) / "array" / plan_sha256 / str(args.task_index)
    report = fsi.download_from_manifest(
        manifest_path=manifest, download_root=download_root, provider_filter=task["provider"],
        overwrite=False, headers=fsi.parse_http_headers(args.http_header, args.auth_bearer_token_env),
        timeout=args.download_timeout, dry_run=False, jobs=int(os.environ.get("GG_TASK_CPUS", "1")),
        resolved_manifest_output_path=root / (str(args.task_index) + ".resolved.tsv"))
    if report["errors"]:
        raise ValueError("; ".join(report["errors"]))
    tasks, warnings, errors = fsi.discover_tasks(task["provider"], download_root / DEFAULT_INPUT_RELATIVE_DIRS[task["provider"]])
    tasks = [candidate for candidate in tasks if candidate["species_prefix"] == task["species_prefix"]]
    if errors or len(tasks) != 1:
        raise ValueError("Expected exactly one downloaded species task: " + repr(errors))
    actual = tasks[0]
    actual["input_sha256"] = {**task.get("input_sha256", {}), **{str(actual[key]): digest(actual[key]) for key in ("cds_path", "gff_path", "gbff_path", "genome_path") if actual.get(key)}}
    for key in ("gene_grouping_mode", "gff_repair_mode", "format_strict"):
        actual[key] = task[key]
    atomic_json(resolved, {"plan_sha256": plan_sha256, "task_index": args.task_index,
                          "task": {key: str(value) if isinstance(value, Path) else value for key, value in actual.items()}})
    return actual


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

    if args.dry_run and "manifest_row" in tasks[args.task_index - 1]:
        parser.error("Use the array submission helper for a download-free manifest preview")
    task = resolve_manifest_task(deserialize_task(tasks[args.task_index - 1]), args)
    for path, expected in task.get("input_sha256", {}).items():
        if digest(path) != expected:
            parser.error("Raw input changed after planning: " + path)
    output_cds_dir = Path(args.species_cds_dir).expanduser().resolve()
    output_gff_dir = Path(args.species_gff_dir).expanduser().resolve()
    output_genome_dir = Path(args.species_genome_dir).expanduser().resolve()
    if task.get("cds_path") is not None:
        cds_output_name = fsi.normalize_cds_output_basename(
            task["cds_path"].name, task["species_prefix"]
        )
    else:
        cds_output_name = fsi.build_derived_cds_output_basename(task)
    if task.get("gff_path") is not None:
        gff_output_name = fsi.normalize_gff_output_basename(
            task["gff_path"].name, task["species_prefix"]
        )
    elif task.get("gbff_path") is not None:
        gff_output_name = fsi.build_derived_gff_output_basename(task)
    else:
        gff_output_name = ""
    if task.get("genome_path") is not None:
        genome_output_name = fsi.normalize_genome_output_basename(
            task["genome_path"].name, task["species_prefix"]
        )
    elif task.get("gbff_path") is not None:
        genome_output_name = fsi.build_derived_genome_output_basename(task)
    else:
        genome_output_name = ""
    expected_paths = {
        "cds_output_path": str(output_cds_dir / cds_output_name),
        "gff_output_path": str(output_gff_dir / gff_output_name) if gff_output_name else "",
        "genome_output_path": (
            str(output_genome_dir / genome_output_name) if genome_output_name else ""
        ),
    }
    raw_paths = {
        key: str(task.get(key) or "")
        for key in ("cds_path", "gff_path", "gbff_path", "genome_path")
    }
    if args.describe_only:
        write_json(
            args.task_meta_output,
            {
                "task_index": args.task_index,
                "task_count": len(tasks),
                "species_prefix": task["species_prefix"],
                "species_key": task["species_key"],
                "provider": task["provider"],
                **raw_paths,
                **expected_paths,
            },
        )
        return 0
    output_cds_dir.mkdir(parents=True, exist_ok=True)
    output_gff_dir.mkdir(parents=True, exist_ok=True)
    output_genome_dir.mkdir(parents=True, exist_ok=True)

    cds_result = fsi.format_cds(
        task,
        output_cds_dir,
        args.overwrite,
        args.dry_run,
        reuse_existing=args.reuse_existing,
    )
    gff_result = fsi.format_gff(
        task,
        output_gff_dir,
        args.overwrite,
        args.dry_run,
        formatted_cds_path=cds_result.get("output_path"),
        reuse_existing=args.reuse_existing,
    )
    genome_result = fsi.format_genome(task, output_genome_dir, args.overwrite, args.dry_run)

    if int(gff_result.get("invalid_utf8_bytes", 0) or 0) > 0:
        sys.stderr.write(
            "Warning: [{}] {}: replaced {} invalid UTF-8 byte(s) in source GFF "
            "across {} line(s); source line sample={}.\n".format(
                task["provider"],
                task["species_prefix"],
                gff_result.get("invalid_utf8_bytes", 0),
                gff_result.get("invalid_utf8_line_count", 0),
                ",".join(str(value) for value in gff_result.get("invalid_utf8_lines", ())) or "NA",
            )
        )
    if int(cds_result.get("gff_mapping_fallback_tolerated", 0) or 0) > 0:
        sys.stderr.write(
            "Warning: [{}] {}: retained {} low-rate CDS/GFF fallback record(s); "
            "strict mode would reject this mismatch.\n".format(
                task["provider"],
                task["species_prefix"],
                cds_result.get("gff_unexpected_mapping_records", 0),
            )
        )

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
        "cds_grouping_source": cds_result.get("grouping_source", "header"),
        "cds_gff_records_mapped": cds_result.get("gff_records_mapped", 0),
        "cds_gff_records_unmapped": cds_result.get("gff_records_unmapped", 0),
        "cds_gff_records_ambiguous": cds_result.get("gff_records_ambiguous", 0),
        "cds_gff_fallback_records": cds_result.get("gff_unexpected_mapping_records", 0),
        "cds_gff_coordinate_rescued_transcripts": cds_result.get(
            "gff_coordinate_rescued_transcripts", 0
        ),
        "cds_gff_coordinate_rescued_groups": cds_result.get("gff_coordinate_rescued_groups", 0),
        "gff_repaired_gene_ids": gff_result.get("repair_gene_ids", 0),
        "gff_repaired_references": gff_result.get("repair_references", 0),
        "gff_repair_ambiguous": gff_result.get("repair_ambiguous", 0),
        "gff_repair_collisions": gff_result.get("repair_collisions", 0),
        "gff_repair_mode": gff_result.get("repair_mode", ""),
        "gff_normalized_bare_attribute_lines": gff_result.get("normalized_bare_attribute_lines", 0),
        "gff_invalid_utf8_bytes": gff_result.get("invalid_utf8_bytes", 0),
        "gff_invalid_utf8_sequences": gff_result.get("invalid_utf8_sequences", 0),
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
            **raw_paths,
            "cds_output_path": str(cds_result.get("output_path") or ""),
            "gff_output_path": str(gff_result.get("output_path") or ""),
            "genome_output_path": str(genome_result.get("output_path") or ""),
            "cds_status": cds_result.get("status", ""),
            "cds_grouping_source": cds_result.get("grouping_source", "header"),
            "cds_gff_grouping_audit_path": str(cds_result.get("gff_grouping_audit_path") or ""),
            "cds_gff_records_mapped": cds_result.get("gff_records_mapped", 0),
            "cds_gff_records_unmapped": cds_result.get("gff_records_unmapped", 0),
            "cds_gff_records_ambiguous": cds_result.get("gff_records_ambiguous", 0),
            "gff_status": gff_result.get("status", ""),
            "gff_repair_status": gff_result.get("repair_status", ""),
            "gff_repair_audit_path": str(gff_result.get("repair_audit_path") or ""),
            "gff_repaired_gene_ids": gff_result.get("repair_gene_ids", 0),
            "gff_repaired_references": gff_result.get("repair_references", 0),
            "genome_status": genome_result.get("status", ""),
        },
    )

    print(
        "[{}] {}: CDS={} ({}), GFF={} ({}), GENOME={} ({})".format(
            task["provider"],
            task["species_prefix"],
            task["cds_path"].name if task.get("cds_path") is not None else "DERIVED_FROM_GFF_GENOME",
            cds_result["status"],
            task["gff_path"].name if task.get("gff_path") else "DERIVED_FROM_GBFF",
            gff_result["status"],
            task["genome_path"].name if task.get("genome_path") is not None else "NA",
            genome_result["status"],
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
