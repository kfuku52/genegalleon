#!/usr/bin/env python3

import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import gzip
import json
from pathlib import Path
import re
import sys

import format_species_inputs as formatter


FASTA_EXTENSIONS = (
    ".fa",
    ".fas",
    ".fasta",
    ".fna",
    ".fa.gz",
    ".fas.gz",
    ".fasta.gz",
    ".fna.gz",
)

DERIVED_CDS_SUFFIX = " (derived CDS)"


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Validate that each formatted species CDS FASTA contains exactly the "
            "expected longest CDS representative per gene."
        )
    )
    parser.add_argument("--species-cds-dir", required=True, help="Directory containing formatted species CDS FASTA files.")
    parser.add_argument(
        "--species-summary",
        required=True,
        help="Path to gg_input_generation_species.tsv written by format_species_inputs.py.",
    )
    parser.add_argument(
        "--nthreads",
        type=int,
        default=None,
        help="Number of species validations to run concurrently.",
    )
    parser.add_argument(
        "--ncpu",
        type=int,
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--missing-limit",
        type=int,
        default=10,
        help="Maximum number of mismatching gene IDs to print per species.",
    )
    parser.add_argument(
        "--stats-output",
        default="",
        help="Optional JSON path for summary stats.",
    )
    return parser


def has_known_suffix(name, suffixes):
    lower = name.lower()
    return any(lower.endswith(suffix) for suffix in suffixes)


def list_nonhidden_files(directory, suffixes):
    out = []
    for path in sorted(Path(directory).iterdir()):
        if not path.is_file():
            continue
        if path.name.startswith("."):
            continue
        if has_known_suffix(path.name, suffixes):
            out.append(path)
    return out


def resolve_nthreads(args):
    if args.nthreads is not None:
        return max(1, int(args.nthreads))
    if args.ncpu is not None:
        return max(1, int(args.ncpu))
    return 1


def first_n(items, limit):
    return items[: max(0, int(limit))]


def read_species_summary_rows(path):
    rows = []
    with open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for index, row in enumerate(reader, start=2):
            clean = {key: str(value or "") for key, value in row.items()}
            clean["_summary_line"] = str(index)
            rows.append(clean)
    return rows


def index_summary_rows_by_output_name(rows):
    indexed = {}
    errors = []
    for row in rows:
        output_path = str(row.get("cds_output_path") or "").strip()
        if output_path == "":
            continue
        output_name = Path(output_path).name
        previous = indexed.get(output_name)
        if previous is not None:
            errors.append(
                "Duplicate species summary rows target the same CDS output '{}': lines {} and {}".format(
                    output_name,
                    previous.get("_summary_line", "?"),
                    row.get("_summary_line", "?"),
                )
            )
            continue
        indexed[output_name] = row
    return indexed, errors


def parse_derived_cds_description(description):
    parsed = {
        "gbff_path": None,
    }
    text = str(description or "").strip()
    if text == "" or not text.endswith(DERIVED_CDS_SUFFIX):
        return parsed
    core = text[: -len(DERIVED_CDS_SUFFIX)].strip()
    if core == "":
        return parsed
    for item in [part.strip() for part in core.split(" + ") if part.strip() != ""]:
        if formatter.is_gbff_filename(Path(item).name):
            parsed["gbff_path"] = Path(item).expanduser().resolve()
    return parsed


def build_task_from_summary_row(row):
    provider = str(row.get("provider") or "").strip()
    species_key = str(row.get("species_key") or row.get("species_prefix") or "").strip()
    species_prefix = str(row.get("species_prefix") or "").strip()
    if provider == "" or species_key == "" or species_prefix == "":
        raise ValueError(
            "species summary line {} is missing provider/species_key/species_prefix".format(
                row.get("_summary_line", "?")
            )
        )

    task = {
        "provider": provider,
        "species_key": species_key,
        "species_prefix": species_prefix,
        "gene_grouping_mode": formatter.normalize_gene_grouping_mode(row.get("gene_grouping_mode", "strict")),
        "cds_path": None,
        "gff_path": None,
        "gbff_path": None,
        "genome_path": None,
    }

    cds_input_path = str(row.get("cds_input_path") or "").strip()
    if cds_input_path != "" and not cds_input_path.endswith(DERIVED_CDS_SUFFIX):
        task["cds_path"] = Path(cds_input_path).expanduser().resolve()
    else:
        derived = parse_derived_cds_description(cds_input_path)
        if derived.get("gbff_path") is not None:
            task["gbff_path"] = derived["gbff_path"]

    gff_input_path = str(row.get("gff_input_path") or "").strip()
    if gff_input_path != "":
        task["gff_path"] = Path(gff_input_path).expanduser().resolve()

    genome_input_path = str(row.get("genome_input_path") or "").strip()
    if genome_input_path != "":
        task["genome_path"] = Path(genome_input_path).expanduser().resolve()

    return task


def validate_task_inputs_exist(task):
    if task.get("cds_path") is not None:
        if not task["cds_path"].exists():
            raise FileNotFoundError("raw CDS input not found: {}".format(task["cds_path"]))
        return
    if task.get("gff_path") is not None and task.get("genome_path") is not None:
        if not task["gff_path"].exists():
            raise FileNotFoundError("raw GFF input not found: {}".format(task["gff_path"]))
        if not task["genome_path"].exists():
            raise FileNotFoundError("raw genome input not found: {}".format(task["genome_path"]))
        return
    if task.get("gbff_path") is not None:
        if not task["gbff_path"].exists():
            raise FileNotFoundError("raw GBFF input not found: {}".format(task["gbff_path"]))
        if task.get("genome_path") is not None and not task["genome_path"].exists():
            raise FileNotFoundError("raw genome input not found: {}".format(task["genome_path"]))
        return
    raise ValueError(
        "species '{}' does not have enough source information in species summary to re-derive CDS".format(
            task.get("species_prefix", "")
        )
    )


def read_fasta_records(path):
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    records = {}
    duplicate_ids = []
    current_id = None
    chunks = []

    def flush():
        nonlocal current_id, chunks
        if current_id is None:
            return
        sequence = re.sub(r"\s+", "", "".join(chunks)).upper()
        if current_id in records:
            duplicate_ids.append(current_id)
        else:
            records[current_id] = sequence
        current_id = None
        chunks = []

    with opener(path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if line == "":
                continue
            if line.startswith(">"):
                flush()
                current_id = line[1:].strip().split()[0]
                chunks = []
                continue
            chunks.append(line)
    flush()
    return records, sorted(set(duplicate_ids))


def collect_expected_longest_records(task):
    expected = {}
    transcripts_by_gene = defaultdict(int)
    transcript_total = 0

    for header, sequence in formatter.iter_task_cds_records(task):
        transcript_total += 1
        transcript_id = formatter.build_formatted_cds_id(task, header)
        gene_id = formatter.build_gene_aggregate_id(task, header, transcript_id)
        seq = formatter.pad_to_codon_length(re.sub(r"\s+", "", sequence).upper())
        transcripts_by_gene[gene_id] += 1
        previous = expected.get(gene_id)
        if previous is None:
            expected[gene_id] = {
                "sequence": seq,
                "transcript_id": transcript_id,
            }
            continue
        previous_seq = previous["sequence"]
        previous_transcript_id = previous["transcript_id"]
        if len(seq) > len(previous_seq) or (
            len(seq) == len(previous_seq) and transcript_id < previous_transcript_id
        ):
            expected[gene_id] = {
                "sequence": seq,
                "transcript_id": transcript_id,
            }

    if transcript_total == 0:
        raise ValueError("no CDS records could be re-derived from source inputs")

    return expected, {
        "transcripts_total": transcript_total,
        "genes_total": len(expected),
        "multi_isoform_genes": sum(1 for count in transcripts_by_gene.values() if count > 1),
        "aggregated_cds_removed": max(0, transcript_total - len(expected)),
    }


def validate_single_species(task, missing_limit):
    cds_file = task["cds_file"]
    row = task["summary_row"]
    species_prefix = str(row.get("species_prefix") or "").strip()
    try:
        format_task = build_task_from_summary_row(row)
        validate_task_inputs_exist(format_task)
        expected_records, expected_stats = collect_expected_longest_records(format_task)
        actual_records, duplicate_output_ids = read_fasta_records(cds_file)

        expected_id_set = set(expected_records.keys())
        actual_id_set = set(actual_records.keys())
        missing_ids = sorted(expected_id_set - actual_id_set)
        extra_ids = sorted(actual_id_set - expected_id_set)
        shared_ids = sorted(expected_id_set.intersection(actual_id_set))
        length_mismatch_ids = sorted(
            gene_id
            for gene_id in shared_ids
            if len(actual_records[gene_id]) != len(expected_records[gene_id]["sequence"])
        )
        sequence_mismatch_ids = sorted(
            gene_id
            for gene_id in shared_ids
            if actual_records[gene_id] != expected_records[gene_id]["sequence"]
        )

        summary_before = str(row.get("cds_sequences_before") or "").strip()
        summary_after = str(row.get("cds_sequences_after") or "").strip()
        summary_aggregated = str(row.get("aggregated_cds_removed") or "").strip()
        count_mismatch_parts = []
        if summary_before != "" and int(summary_before) != expected_stats["transcripts_total"]:
            count_mismatch_parts.append(
                "summary_before={} expected_before={}".format(summary_before, expected_stats["transcripts_total"])
            )
        if summary_after != "" and int(summary_after) != expected_stats["genes_total"]:
            count_mismatch_parts.append(
                "summary_after={} expected_after={}".format(summary_after, expected_stats["genes_total"])
            )
        if summary_aggregated != "" and int(summary_aggregated) != expected_stats["aggregated_cds_removed"]:
            count_mismatch_parts.append(
                "summary_aggregated={} expected_aggregated={}".format(
                    summary_aggregated, expected_stats["aggregated_cds_removed"]
                )
            )

        stats = {
            "genes_total": expected_stats["genes_total"],
            "transcripts_total": expected_stats["transcripts_total"],
            "multi_isoform_genes": expected_stats["multi_isoform_genes"],
            "aggregated_cds_removed": expected_stats["aggregated_cds_removed"],
            "output_records": len(actual_records),
        }

        failure_parts = []
        if len(duplicate_output_ids) > 0:
            failure_parts.append(
                "duplicate_output_ids={} sample={}".format(
                    len(duplicate_output_ids), ",".join(first_n(duplicate_output_ids, missing_limit))
                )
            )
        if len(missing_ids) > 0:
            failure_parts.append(
                "missing={} sample={}".format(
                    len(missing_ids), ",".join(first_n(missing_ids, missing_limit))
                )
            )
        if len(extra_ids) > 0:
            failure_parts.append(
                "extra={} sample={}".format(
                    len(extra_ids), ",".join(first_n(extra_ids, missing_limit))
                )
            )
        if len(length_mismatch_ids) > 0:
            failure_parts.append(
                "length_mismatch={} sample={}".format(
                    len(length_mismatch_ids), ",".join(first_n(length_mismatch_ids, missing_limit))
                )
            )
        if len(sequence_mismatch_ids) > 0:
            failure_parts.append(
                "sequence_mismatch={} sample={}".format(
                    len(sequence_mismatch_ids), ",".join(first_n(sequence_mismatch_ids, missing_limit))
                )
            )
        failure_parts.extend(count_mismatch_parts)

        if len(failure_parts) > 0:
            return {
                "index": task["index"],
                "species_prefix": species_prefix,
                "ok": False,
                "stats_ready": True,
                "stats": stats,
                "error": "[{}] Longest CDS validation failed for {}: {}".format(
                    species_prefix,
                    cds_file.name,
                    "; ".join(failure_parts),
                ),
            }

        return {
            "index": task["index"],
            "species_prefix": species_prefix,
            "ok": True,
            "stats_ready": True,
            "stats": stats,
            "message": (
                "[{}] Longest CDS validation OK: genes={}, transcripts={}, multi_isoform_genes={}, aggregated={}"
            ).format(
                species_prefix,
                expected_stats["genes_total"],
                expected_stats["transcripts_total"],
                expected_stats["multi_isoform_genes"],
                expected_stats["aggregated_cds_removed"],
            ),
        }
    except Exception as exc:
        return {
            "index": task["index"],
            "species_prefix": species_prefix or cds_file.name,
            "ok": False,
            "stats_ready": False,
            "error": "[{}] Unexpected longest CDS validation error for {}: {}".format(
                species_prefix or cds_file.name,
                cds_file.name,
                exc,
            ),
        }


def run_validation_tasks(tasks, missing_limit, nthreads):
    if len(tasks) == 0:
        return []
    if nthreads <= 1 or len(tasks) <= 1:
        return [validate_single_species(task=task, missing_limit=missing_limit) for task in tasks]

    ordered_results = {}
    max_workers = min(nthreads, len(tasks))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_task = {
            executor.submit(validate_single_species, task, missing_limit): task
            for task in tasks
        }
        for future in as_completed(future_to_task):
            task = future_to_task[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {
                    "index": task["index"],
                    "species_prefix": task["summary_row"].get("species_prefix", ""),
                    "ok": False,
                    "stats_ready": False,
                    "error": "[{}] Unexpected executor error for {}: {}".format(
                        task["summary_row"].get("species_prefix", ""),
                        task["cds_file"].name,
                        exc,
                    ),
                }
            ordered_results[result["index"]] = result
    return [ordered_results[index] for index in sorted(ordered_results)]


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    nthreads = resolve_nthreads(args)

    cds_dir = Path(args.species_cds_dir).expanduser().resolve()
    species_summary = Path(args.species_summary).expanduser().resolve()
    if not cds_dir.is_dir():
        parser.error("--species-cds-dir not found: {}".format(cds_dir))
    if not species_summary.is_file():
        parser.error("--species-summary not found: {}".format(species_summary))

    cds_files = list_nonhidden_files(cds_dir, FASTA_EXTENSIONS)
    if len(cds_files) == 0:
        parser.error("No CDS FASTA files found in: {}".format(cds_dir))

    summary_rows = read_species_summary_rows(species_summary)
    if len(summary_rows) == 0:
        parser.error("No rows found in species summary: {}".format(species_summary))

    summary_by_output_name, index_errors = index_summary_rows_by_output_name(summary_rows)
    warnings = []
    errors = list(index_errors)
    stats = {
        "species_checked": 0,
        "species_passed": 0,
        "genes_total": 0,
        "transcripts_total": 0,
        "multi_isoform_genes_total": 0,
        "aggregated_cds_removed_total": 0,
        "nthreads": nthreads,
    }

    tasks = []
    for index, cds_file in enumerate(cds_files):
        row = summary_by_output_name.get(cds_file.name)
        if row is None:
            errors.append(
                "No matching species summary row found for formatted CDS output '{}' in {}".format(
                    cds_file.name,
                    species_summary,
                )
            )
            continue
        tasks.append(
            {
                "index": index,
                "cds_file": cds_file,
                "summary_row": row,
            }
        )

    current_output_names = {path.name for path in cds_files}
    for row in summary_rows:
        output_path = str(row.get("cds_output_path") or "").strip()
        if output_path == "":
            continue
        output_name = Path(output_path).name
        if output_name not in current_output_names:
            warnings.append(
                "Ignoring species summary row for '{}' because no current formatted CDS file exists in {}".format(
                    output_name,
                    cds_dir,
                )
            )

    for result in run_validation_tasks(tasks=tasks, missing_limit=args.missing_limit, nthreads=nthreads):
        if result.get("stats_ready", False):
            stats["species_checked"] += 1
            stats["genes_total"] += int(result["stats"]["genes_total"])
            stats["transcripts_total"] += int(result["stats"]["transcripts_total"])
            stats["multi_isoform_genes_total"] += int(result["stats"]["multi_isoform_genes"])
            stats["aggregated_cds_removed_total"] += int(result["stats"]["aggregated_cds_removed"])
        if result["ok"]:
            stats["species_passed"] += 1
            print(result["message"])
            continue
        errors.append(result["error"])

    if args.stats_output != "":
        stats_path = Path(args.stats_output).expanduser().resolve()
        stats_path.parent.mkdir(parents=True, exist_ok=True)
        with open(stats_path, "wt", encoding="utf-8") as handle:
            json.dump(stats, handle, ensure_ascii=True, indent=2, sort_keys=True)

    for warning in warnings:
        sys.stderr.write("Warning: {}\n".format(warning))
    for error in errors:
        sys.stderr.write("Error: {}\n".format(error))

    print(
        "Longest CDS validation summary: species_checked={}, species_passed={}, genes_total={}, transcripts_total={}, multi_isoform_genes_total={}, aggregated_cds_removed_total={}".format(
            stats["species_checked"],
            stats["species_passed"],
            stats["genes_total"],
            stats["transcripts_total"],
            stats["multi_isoform_genes_total"],
            stats["aggregated_cds_removed_total"],
        )
    )

    if len(errors) > 0:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
