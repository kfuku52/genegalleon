"""Persistent per-species run summary handling."""

import csv
import os
from pathlib import Path

from format_species_annotations import describe_task_genome_input, describe_task_gff_input, gene_grouping_mode_for_task
from format_species_constants import SPECIES_SUMMARY_COLUMNS
from format_species_discovery import utc_now_iso
from format_species_taxonomy import blank_species_taxonomy_metadata


def species_row_key(provider, species_key, species_prefix):
    stable_species = species_key if species_key != "" else species_prefix
    return "{}\t{}".format(provider, stable_species)


def read_species_summary_rows(path):
    rows = {}
    if not path.exists() or path.stat().st_size == 0:
        return rows
    with open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for raw in reader:
            provider = (raw.get("provider") or "").strip()
            species_key = (raw.get("species_key") or "").strip()
            species_prefix = (raw.get("species_prefix") or "").strip()
            if provider == "" or (species_key == "" and species_prefix == ""):
                continue
            key = species_row_key(provider, species_key, species_prefix)
            rows[key] = {col: str(raw.get(col) or "") for col in SPECIES_SUMMARY_COLUMNS}
    return rows


def species_summary_row_has_existing_outputs(row):
    cds_output = (row.get("cds_output_path") or "").strip()
    gff_output = (row.get("gff_output_path") or "").strip()
    genome_output = (row.get("genome_output_path") or "").strip()
    gff_status = (row.get("gff_status") or "").strip()
    genome_status = (row.get("genome_status") or "").strip()

    if cds_output == "":
        return False
    if not Path(cds_output).expanduser().exists():
        return False
    if gff_status != "missing":
        if gff_output == "":
            return False
        if not Path(gff_output).expanduser().exists():
            return False
    if genome_output != "" and genome_status != "missing" and not Path(genome_output).expanduser().exists():
        return False
    return True


def retain_existing_species_summary_rows(rows):
    retained = {}
    for key in sorted(rows.keys()):
        row = rows[key]
        if species_summary_row_has_existing_outputs(row):
            retained[key] = row
    return retained


def write_species_summary_rows(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = Path(str(path) + ".tmp.{}".format(os.getpid()))
    try:
        with open(tmp_path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=SPECIES_SUMMARY_COLUMNS, delimiter="\t")
            writer.writeheader()
            for key in sorted(rows.keys()):
                row = rows[key]
                payload = {col: str(row.get(col) or "") for col in SPECIES_SUMMARY_COLUMNS}
                writer.writerow(payload)
        tmp_path.replace(path)
    except Exception:
        try:
            tmp_path.unlink()
        except FileNotFoundError:
            pass
        except OSError:
            pass
        raise


def format_task_succeeded(cds_result, gff_result, genome_result, dry_run):
    if dry_run:
        return False
    if cds_result["status"] not in ("write", "skip"):
        return False
    if gff_result["status"] not in ("write", "skip", "missing"):
        return False
    if genome_result["status"] not in ("write", "skip", "missing"):
        return False
    return True


def result_output_name(result, fallback="NA"):
    output_path = result.get("output_path")
    if output_path is None:
        return fallback
    return output_path.name


def build_species_summary_row(
    task,
    cds_result,
    gff_result,
    genome_result,
    run_started_utc,
    overwrite,
    dry_run,
    taxonomy_metadata=None,
):
    metadata = blank_species_taxonomy_metadata()
    if taxonomy_metadata is not None:
        metadata.update({key: str(value or "") for key, value in taxonomy_metadata.items()})
    return {
        "updated_utc": utc_now_iso(),
        "run_started_utc": run_started_utc,
        "provider": task["provider"],
        "species_key": task["species_key"],
        "species_prefix": task["species_prefix"],
        "taxid": metadata["taxid"],
        "nuclear_genetic_code_id": metadata["nuclear_genetic_code_id"],
        "nuclear_genetic_code_name": metadata["nuclear_genetic_code_name"],
        "mitochondrial_genetic_code_id": metadata["mitochondrial_genetic_code_id"],
        "mitochondrial_genetic_code_name": metadata["mitochondrial_genetic_code_name"],
        "plastid_genetic_code_id": metadata["plastid_genetic_code_id"],
        "plastid_genetic_code_name": metadata["plastid_genetic_code_name"],
        "cds_input_path": str(cds_result.get("input_path") or ""),
        "gff_input_path": str(describe_task_gff_input(task) or ""),
        "genome_input_path": str(describe_task_genome_input(task) or ""),
        "cds_output_path": str(cds_result["output_path"]) if cds_result.get("output_path") is not None else "",
        "gff_output_path": str(gff_result["output_path"]) if gff_result.get("output_path") is not None else "",
        "genome_output_path": str(genome_result["output_path"]) if genome_result.get("output_path") is not None else "",
        "cds_status": cds_result["status"],
        "gff_status": gff_result["status"],
        "gff_repair_mode": gff_result.get("repair_mode", ""),
        "gff_repair_status": gff_result.get("repair_status", ""),
        "gff_repair_audit_path": str(gff_result.get("repair_audit_path") or ""),
        "gff_repaired_gene_ids": str(gff_result.get("repair_gene_ids", 0)),
        "gff_repaired_references": str(gff_result.get("repair_references", 0)),
        "gff_repair_ambiguous": str(gff_result.get("repair_ambiguous", 0)),
        "gff_repair_collisions": str(gff_result.get("repair_collisions", 0)),
        "genome_status": genome_result["status"],
        "cds_sequences_before": str(cds_result.get("before_count", "")),
        "cds_sequences_after": str(cds_result.get("after_count", "")),
        "cds_first_sequence_name": cds_result.get("first_sequence_name", "") or "NA",
        "aggregated_cds_removed": str(cds_result.get("duplicates", "")),
        "gene_grouping_mode": gene_grouping_mode_for_task(task),
        "overwrite": str(int(bool(overwrite))),
        "dry_run": str(int(bool(dry_run))),
    }
