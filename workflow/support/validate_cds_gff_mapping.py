#!/usr/bin/env python3

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


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

GFF_EXTENSIONS = (
    ".gff",
    ".gff3",
    ".gtf",
    ".gff.gz",
    ".gff3.gz",
    ".gtf.gz",
)

KNOWN_ALLOWED_MISSING_CDS_IDS = {
    # FernBase Azolla_filiculoides publishes two CDS-only entries whose
    # formatted IDs remain absent from the companion GFF.
    "Azolla_filiculoides": {
        "Azolla_filiculoides_Azfi_s0034.g025227",
        "Azolla_filiculoides_Azfi_s0093.g043301",
    },
}
TAXONOMIC_PROXIMITY_QUALIFIERS = frozenset(("cf", "aff", "nr"))
TAXONOMIC_INFRASPECIFIC_RANKS = frozenset(("subsp", "ssp", "var", "forma", "f"))


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Validate that formatted species CDS IDs can be resolved back to species GFF "
            "using the same gff2genestat logic used downstream."
        )
    )
    parser.add_argument("--species-cds-dir", required=True, help="Directory containing formatted species CDS FASTA files.")
    parser.add_argument("--species-gff-dir", required=True, help="Directory containing formatted species GFF files.")
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
        help="Maximum number of missing gene IDs to print per species.",
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


def species_prefix_from_name(name):
    parts = [part for part in Path(name).name.split("_") if part != ""]
    if len(parts) < 2:
        return ""
    second = parts[1].lower()
    third = parts[2].lower() if len(parts) >= 3 else ""
    count = 2
    if second == "sp":
        count = 3 if len(parts) >= 3 else 2
    elif second in TAXONOMIC_PROXIMITY_QUALIFIERS:
        count = 3 if len(parts) >= 3 else 2
    elif third in TAXONOMIC_PROXIMITY_QUALIFIERS:
        count = 3
    elif third in TAXONOMIC_INFRASPECIFIC_RANKS:
        count = 4 if len(parts) >= 4 else 3
    return "_".join(parts[:count]).split(".", 1)[0]


def read_fasta_ids(path):
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    ids = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            ids.append(line[1:].strip().split()[0])
    return ids


def read_gene_ids_from_tsv(path):
    if not path.exists() or path.stat().st_size == 0:
        return []
    ids = []
    with open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene_id = str(row.get("gene_id") or "").strip()
            if gene_id != "":
                ids.append(gene_id)
    return ids


def symlink_or_copy(src, dst):
    try:
        os.symlink(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def resolve_nthreads(args):
    if args.nthreads is not None:
        return max(1, int(args.nthreads))
    if args.ncpu is not None:
        return max(1, int(args.ncpu))
    return 1


def run_gff2genestat(cds_file, gff_file):
    script_dir = Path(__file__).resolve().parent
    gff2genestat = script_dir / "gff2genestat.py"
    with tempfile.TemporaryDirectory(prefix="validate_cds_gff_mapping_") as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        symlink_or_copy(str(gff_file), str(tmp_dir / gff_file.name))
        out_path = tmp_dir / "gff2genestat.tsv"
        cmd = [
            sys.executable,
            str(gff2genestat),
            "--dir_gff",
            str(tmp_dir),
            "--seqfile",
            str(cds_file),
            "--outfile",
            str(out_path),
            "--feature",
            "CDS",
            "--multiple_hits",
            "longest",
            "--ncpu",
            "1",
        ]
        completed = subprocess.run(cmd, capture_output=True, text=True, check=False)
        mapped_ids = read_gene_ids_from_tsv(out_path)
        return completed, mapped_ids


def first_n(items, limit):
    return items[: max(0, int(limit))]


def index_gff_files_by_species(gff_files):
    out = {}
    for path in gff_files:
        species_prefix = species_prefix_from_name(path.name)
        if species_prefix == "":
            continue
        out.setdefault(species_prefix, []).append(path)
    return out


def validate_single_species(task, missing_limit):
    species_prefix = task["species_prefix"]
    cds_file = task["cds_file"]
    gff_file = task["gff_file"]
    try:
        cds_ids = read_fasta_ids(cds_file)
        if len(cds_ids) == 0:
            return {
                "index": task["index"],
                "species_prefix": species_prefix,
                "ok": False,
                "stats_ready": False,
                "error": "[{}] No CDS IDs were found in {}".format(species_prefix, cds_file),
            }

        completed, mapped_ids = run_gff2genestat(cds_file=cds_file, gff_file=gff_file)
        if completed.returncode != 0:
            detail = "\n".join(
                text for text in (completed.stdout.strip(), completed.stderr.strip()) if text
            ).strip()
            if detail == "":
                detail = "no stderr/stdout"
            return {
                "index": task["index"],
                "species_prefix": species_prefix,
                "ok": False,
                "stats_ready": False,
                "error": "[{}] gff2genestat.py failed for CDS={} GFF={}: {}".format(
                    species_prefix, cds_file.name, gff_file.name, detail
                ),
            }

        cds_id_set = set(cds_ids)
        mapped_id_set = set(mapped_ids)
        missing_ids = sorted(cds_id_set - mapped_id_set)
        extra_ids = sorted(mapped_id_set - cds_id_set)
        allowed_missing_id_set = KNOWN_ALLOWED_MISSING_CDS_IDS.get(species_prefix, set())
        allowed_missing_ids = sorted(cds_id_set.intersection(allowed_missing_id_set) - mapped_id_set)
        unexpected_missing_ids = sorted(set(missing_ids) - set(allowed_missing_ids))
        stats = {
            "cds_ids": len(cds_ids),
            "mapped_ids": len(mapped_ids),
            "allowed_missing_ids": len(allowed_missing_ids),
        }

        if len(unexpected_missing_ids) > 0 or len(extra_ids) > 0:
            parts = []
            if len(unexpected_missing_ids) > 0:
                parts.append(
                    "missing={} sample={}".format(
                        len(unexpected_missing_ids), ",".join(first_n(unexpected_missing_ids, missing_limit))
                    )
                )
            if len(extra_ids) > 0:
                parts.append(
                    "extra={} sample={}".format(
                        len(extra_ids), ",".join(first_n(extra_ids, missing_limit))
                    )
                )
            return {
                "index": task["index"],
                "species_prefix": species_prefix,
                "ok": False,
                "stats_ready": True,
                "stats": stats,
                "error": "[{}] CDS-to-GFF mapping failed for CDS={} GFF={}: {}".format(
                    species_prefix, cds_file.name, gff_file.name, "; ".join(parts)
                ),
            }

        return {
            "index": task["index"],
            "species_prefix": species_prefix,
            "ok": True,
            "stats_ready": True,
            "stats": stats,
            "message": "[{}] CDS-to-GFF mapping OK: {}/{} IDs{}".format(
                species_prefix,
                len(mapped_ids),
                len(cds_ids),
                "" if len(allowed_missing_ids) == 0 else " (allowed_missing={})".format(len(allowed_missing_ids)),
            ),
        }
    except Exception as exc:
        return {
            "index": task["index"],
            "species_prefix": species_prefix,
            "ok": False,
            "stats_ready": False,
            "error": "[{}] Unexpected validation error for CDS={} GFF={}: {}".format(
                species_prefix, cds_file.name, gff_file.name, exc
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
                    "species_prefix": task["species_prefix"],
                    "ok": False,
                    "stats_ready": False,
                    "error": "[{}] Unexpected executor error for CDS={} GFF={}: {}".format(
                        task["species_prefix"],
                        task["cds_file"].name,
                        task["gff_file"].name,
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
    gff_dir = Path(args.species_gff_dir).expanduser().resolve()
    if not cds_dir.is_dir():
        parser.error("--species-cds-dir not found: {}".format(cds_dir))
    if not gff_dir.is_dir():
        parser.error("--species-gff-dir not found: {}".format(gff_dir))

    cds_files = list_nonhidden_files(cds_dir, FASTA_EXTENSIONS)
    gff_files = list_nonhidden_files(gff_dir, GFF_EXTENSIONS)
    if len(cds_files) == 0:
        parser.error("No CDS FASTA files found in: {}".format(cds_dir))
    if len(gff_files) == 0:
        parser.error("No GFF files found in: {}".format(gff_dir))

    warnings = []
    errors = []
    stats = {
        "species_checked": 0,
        "species_passed": 0,
        "cds_ids_total": 0,
        "mapped_ids_total": 0,
        "nthreads": nthreads,
    }

    gff_by_species = index_gff_files_by_species(gff_files)
    tasks = []
    for index, cds_file in enumerate(cds_files):
        species_prefix = species_prefix_from_name(cds_file.name)
        if species_prefix == "":
            errors.append("Could not parse species prefix from CDS file name: {}".format(cds_file.name))
            continue

        matching_gff = gff_by_species.get(species_prefix, [])
        if len(matching_gff) == 0:
            errors.append("[{}] No matching GFF file found for CDS file {}".format(species_prefix, cds_file.name))
            continue
        if len(matching_gff) > 1:
            warnings.append(
                "[{}] Multiple matching GFF files found. Using '{}'".format(
                    species_prefix, matching_gff[0].name
                )
            )
        tasks.append(
            {
                "index": index,
                "species_prefix": species_prefix,
                "cds_file": cds_file,
                "gff_file": matching_gff[0],
            }
        )

    for result in run_validation_tasks(tasks=tasks, missing_limit=args.missing_limit, nthreads=nthreads):
        if result.get("stats_ready", False):
            stats["species_checked"] += 1
            stats["cds_ids_total"] += int(result["stats"]["cds_ids"])
            stats["mapped_ids_total"] += int(result["stats"]["mapped_ids"])
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
        "CDS-to-GFF mapping validation summary: species_checked={}, species_passed={}, cds_ids_total={}, mapped_ids_total={}".format(
            stats["species_checked"],
            stats["species_passed"],
            stats["cds_ids_total"],
            stats["mapped_ids_total"],
        )
    )

    if len(errors) > 0:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
