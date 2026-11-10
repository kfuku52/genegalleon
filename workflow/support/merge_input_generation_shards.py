#!/usr/bin/env python3

import argparse
import csv
import hashlib
import io
import json
import os
import re
import stat
import sys
import tempfile
from pathlib import Path

import format_species_inputs as fsi
from format_species_summary import parse_species_summary_rows
from input_generation_array_state import parse_plan, read_verified_receipt


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
    parser.add_argument("--mapping-qc-output", default="", help="Optional per-species mapping and annotation QC TSV.")
    parser.add_argument(
        "--expected-task-count",
        type=int,
        default=0,
        help="Optional expected number of completed task shards/stats.",
    )
    parser.add_argument("--resolved-manifest-output", default="")
    parser.add_argument("--task-plan", default="", help="Require verified completion receipts for this exact plan.")
    return parser


def read_task_stats(path):
    with open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def read_payload(path):
    with os.fdopen(os.open(path, os.O_RDONLY | os.O_NONBLOCK), "rb") as handle:
        if not stat.S_ISREG(os.fstat(handle.fileno()).st_mode):
            raise ValueError("Expected a regular shard file: " + str(path))
        return handle.read()


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    species_summary_shard_dir = Path(args.species_summary_shard_dir).expanduser().resolve()
    species_summary_output = Path(args.species_summary_output).expanduser().resolve()
    task_stats_dir = Path(args.task_stats_dir).expanduser().resolve()
    aggregate_stats_output = Path(args.aggregate_stats_output).expanduser().resolve()

    shard_paths = sorted(species_summary_shard_dir.glob("*.tsv")) if species_summary_shard_dir.exists() else []
    all_json_paths = sorted(task_stats_dir.glob("*.json")) if task_stats_dir.exists() else []
    stats_paths = sorted(
        path for path in all_json_paths
        if re.fullmatch(r"[1-9][0-9]*\.json", path.name)
    )
    mapping_qc_paths = sorted(
        path for path in all_json_paths
        if re.fullmatch(r"[1-9][0-9]*\.mapping\.json", path.name)
    )
    ownership_qc_paths = sorted(
        path for path in all_json_paths
        if re.fullmatch(r"[1-9][0-9]*\.longest\.json", path.name)
    )
    if len(stats_paths) + len(mapping_qc_paths) + len(ownership_qc_paths) != len(all_json_paths):
        parser.error("Task stats directory contains an unexpected JSON file")

    if args.expected_task_count > 0 and len(stats_paths) != args.expected_task_count:
        parser.error(
            "Expected {} task stats shards but found {}".format(
                args.expected_task_count, len(stats_paths)
            )
        )

    receipt_files, stats_payloads, summary_rows = {}, {}, {}

    def bound_payload(path, label):
        index = int(path.name.split(".", 1)[0])
        if index not in receipt_files:
            parser.error(label + " index is outside the immutable task plan")
        expected = receipt_files[index].get(str(path.resolve()))
        if expected is None:
            parser.error(label + " is not bound to the worker completion receipt: " + str(path))
        payload = read_payload(path)
        if hashlib.sha256(payload).hexdigest() != expected:
            parser.error(label + " changed after worker completion verification: " + str(path))
        return payload

    def task_stats(path, label):
        if args.task_plan:
            return json.loads(bound_payload(path, label).decode("utf-8"))
        return read_task_stats(path)

    if args.task_plan:
        plan_payload = read_payload(args.task_plan)
        plan = parse_plan(plan_payload.decode("utf-8"))
        indices = range(1, plan["task_count"] + 1)
        plan_sha256 = hashlib.sha256(plan_payload).hexdigest()
        for index in indices:
            receipt = read_verified_receipt(args.task_plan, index, plan_sha256)
            if receipt is None:
                parser.error("Missing, stale, or changed worker completion receipts")
            receipt_files[index] = receipt["files"]
        expected_stats = {str(index) + ".json" for index in indices}
        expected_summaries = {str(index) + ".tsv" for index in indices}
        if {path.name for path in stats_paths} != expected_stats or {path.name for path in shard_paths} != expected_summaries:
            parser.error("Shard identities do not match the immutable task plan")
        for index, task in enumerate(plan["tasks"], start=1):
            stats_path = task_stats_dir / (str(index) + ".json")
            summary_path = species_summary_shard_dir / (str(index) + ".tsv")
            stats = task_stats(stats_path, "Task stats")
            payload = bound_payload(summary_path, "Species summary")
            rows = parse_species_summary_rows(io.StringIO(payload.decode("utf-8"), newline=""))
            if stats.get("species_prefix") != task["species_prefix"] or len(rows) != 1 or next(iter(rows.values()))["species_prefix"] != task["species_prefix"]:
                parser.error("Species identity mismatch in task " + str(index))
            stats_payloads[stats_path] = stats
            summary_rows[summary_path] = rows
    else:
        stats_payloads = {path: read_task_stats(path) for path in stats_paths}
    merged_rows = {} if args.task_plan else fsi.retain_existing_species_summary_rows(
        fsi.read_species_summary_rows(species_summary_output)
    )
    shard_row_count = 0
    for shard_path in shard_paths:
        shard_rows = summary_rows[shard_path] if args.task_plan else fsi.read_species_summary_rows(shard_path)
        shard_row_count += len(shard_rows)
        merged_rows.update(shard_rows)

    species_processed = 0
    num_species_cds_files = 0
    num_species_gff_files = 0
    num_species_genome_files = 0
    cds_sequences_before = 0
    cds_sequences_after = 0
    cds_first_sequence_name = ""
    aggregated_cds_removed = 0
    gff_repaired_gene_ids = 0
    gff_repaired_references = 0
    gff_repair_ambiguous = 0
    gff_repair_collisions = 0
    cds_gff_records_mapped = 0
    cds_gff_records_unmapped = 0
    cds_gff_records_ambiguous = 0
    cds_gff_coordinate_rescued_transcripts = 0
    cds_gff_coordinate_rescued_groups = 0
    phase_conflicts_total = 0
    utr_conflicts_total = 0
    source_gene_validation = {}

    for qc_path in ownership_qc_paths:
        index = qc_path.name.split(".", 1)[0]
        stats_path = task_stats_dir / (index + ".json")
        if stats_path not in stats_payloads:
            parser.error("Source ownership QC has no matching task stats")
        species = stats_payloads[stats_path].get("species_prefix")
        ownership = task_stats(qc_path, "Source ownership QC").get("source_gene_validation")
        if not species or not isinstance(ownership, dict) or set(ownership) != {species}:
            parser.error("Source ownership QC species does not match its task")
        source_gene_validation.update(ownership)

    qc_by_index = {}
    for qc_path in mapping_qc_paths:
        qc = task_stats(qc_path, "Mapping QC")
        if args.mapping_qc_output:
            qc_by_index[int(qc_path.name.split(".", 1)[0])] = qc
        phase_conflicts_total += int(qc.get("phase_conflicts_total", 0) or 0)
        utr_conflicts_total += int(qc.get("utr_conflicts_total", 0) or 0)

    for stats_path in stats_paths:
        payload = stats_payloads[stats_path]
        species_processed += int(payload.get("species_processed", 0) or 0)
        num_species_cds_files += int(payload.get("num_species_cds_files", 0) or 0)
        num_species_gff_files += int(payload.get("num_species_gff_files", 0) or 0)
        num_species_genome_files += int(payload.get("num_species_genome_files", 0) or 0)
        cds_sequences_before += int(payload.get("cds_sequences_before", 0) or 0)
        cds_sequences_after += int(payload.get("cds_sequences_after", 0) or 0)
        aggregated_cds_removed += int(payload.get("aggregated_cds_removed", 0) or 0)
        gff_repaired_gene_ids += int(payload.get("gff_repaired_gene_ids", 0) or 0)
        gff_repaired_references += int(payload.get("gff_repaired_references", 0) or 0)
        gff_repair_ambiguous += int(payload.get("gff_repair_ambiguous", 0) or 0)
        gff_repair_collisions += int(payload.get("gff_repair_collisions", 0) or 0)
        cds_gff_records_mapped += int(payload.get("cds_gff_records_mapped", 0) or 0)
        cds_gff_records_unmapped += int(payload.get("cds_gff_records_unmapped", 0) or 0)
        cds_gff_records_ambiguous += int(payload.get("cds_gff_records_ambiguous", 0) or 0)
        cds_gff_coordinate_rescued_transcripts += int(
            payload.get("cds_gff_coordinate_rescued_transcripts", 0) or 0
        )
        cds_gff_coordinate_rescued_groups += int(
            payload.get("cds_gff_coordinate_rescued_groups", 0) or 0
        )
        if cds_first_sequence_name == "":
            cds_first_sequence_name = str(payload.get("cds_first_sequence_name", "") or "")

    resolved_rows = []
    if args.task_plan and args.resolved_manifest_output and any("manifest_row" in task for task in plan["tasks"]):
        for index, task in enumerate(plan["tasks"], start=1):
            source = Path(args.task_plan + ".tasks") / (str(index) + ".resolved.tsv")
            payload = bound_payload(source, "Resolved manifest")
            with io.StringIO(payload.decode("utf-8"), newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            if len(rows) != 1 or rows[0]["species_key"] != task["species_key"]:
                parser.error("Resolved manifest identity mismatch in task " + str(index))
            resolved_rows.extend(rows)

    if args.mapping_qc_output:
        qc_rows = []
        for stats_path in sorted(stats_paths, key=lambda path: int(path.stem)):
            index = int(stats_path.stem)
            species = str(stats_payloads[stats_path].get("species_prefix") or "").strip()
            if not species:
                parser.error("Task stats is missing species_prefix for mapping QC")
            qc = qc_by_index.get(index)
            qc_rows.append({
                "task_index": index,
                "species_prefix": species,
                "qc_status": "available" if qc is not None else "not_recorded",
                "cds_ids_total": "" if qc is None else qc.get("cds_ids_total", ""),
                "mapped_ids_total": "" if qc is None else qc.get("mapped_ids_total", ""),
                "phase_conflicts_total": "" if qc is None else qc.get("phase_conflicts_total", ""),
                "utr_conflicts_total": "" if qc is None else qc.get("utr_conflicts_total", ""),
            })
        if not qc_rows:
            parser.error("Cannot publish mapping QC without task stats")
        target = Path(args.mapping_qc_output)
        target.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode="w", dir=target.parent, delete=False, newline="") as handle:
            staged_qc = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=list(qc_rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(qc_rows)
        os.replace(staged_qc, target)

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
                "gff_repaired_gene_ids": gff_repaired_gene_ids,
                "gff_repaired_references": gff_repaired_references,
                "gff_repair_ambiguous": gff_repair_ambiguous,
                "gff_repair_collisions": gff_repair_collisions,
                "cds_gff_records_mapped": cds_gff_records_mapped,
                "cds_gff_records_unmapped": cds_gff_records_unmapped,
                "cds_gff_records_ambiguous": cds_gff_records_ambiguous,
                "cds_gff_coordinate_rescued_transcripts": cds_gff_coordinate_rescued_transcripts,
                "cds_gff_coordinate_rescued_groups": cds_gff_coordinate_rescued_groups,
                "merged_species_summary_rows": len(merged_rows),
                "task_stats_files": len(stats_paths),
                "mapping_qc_files": len(mapping_qc_paths),
                "ownership_qc_files": len(ownership_qc_paths),
                "source_gene_validation": source_gene_validation,
                "phase_conflicts_total": phase_conflicts_total,
                "utr_conflicts_total": utr_conflicts_total,
                "task_species_summary_shards": shard_row_count,
            },
            handle,
            ensure_ascii=True,
            indent=2,
            sort_keys=True,
        )

    if resolved_rows:
        target = Path(args.resolved_manifest_output)
        target.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode="w", dir=target.parent, delete=False, newline="") as handle:
            staged_manifest = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=sorted({key for row in resolved_rows for key in row}), delimiter="\t")
            writer.writeheader()
            writer.writerows(resolved_rows)
        os.replace(staged_manifest, target)

    species_summary_output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=species_summary_output.parent, delete=False) as handle:
        staged = Path(handle.name)
    try:
        fsi.write_species_summary_rows(staged, merged_rows)
        os.replace(staged, species_summary_output)
    finally:
        staged.unlink(missing_ok=True)

    print(
        "Merged {} species summary shards and {} task stats shards.".format(
            len(shard_paths), len(stats_paths)
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
