#!/usr/bin/env python3
"""Verify the formatted outputs required by an input-generation dataset."""

import argparse
import csv
import gzip
from pathlib import Path

OUTPUT_COLUMNS = {
    "cds": "cds_output_path",
    "gff": "gff_output_path",
    "genome": "genome_output_path",
}
OUTPUT_LABELS = {"cds": "CDS", "gff": "GFF", "genome": "genome"}


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("rt", encoding="utf-8")


def has_fasta_sequence(path):
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        with open_text(path) as handle:
            first = handle.readline()
            sequence = handle.readline().strip()
            return first.startswith(">") and bool(sequence) and not sequence.startswith(">")
    except (OSError, UnicodeError):
        return False


def has_gff_feature(path):
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        with open_text(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                return len(fields) == 9 and bool(fields[0]) and bool(fields[2])
    except (OSError, UnicodeError):
        return False
    return False


def main(default_required=()):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-summary", required=True, type=Path)
    parser.add_argument("--expected-task-count", type=int, default=0)
    for output in OUTPUT_COLUMNS:
        parser.add_argument("--require-" + output, action="store_true")
    args = parser.parse_args()
    if args.expected_task_count < 0:
        parser.error("--expected-task-count must be nonnegative")
    required = [output for output in OUTPUT_COLUMNS if getattr(args, "require_" + output) or output in default_required]
    if not required:
        parser.error("At least one output must be required")
    if not args.species_summary.is_file():
        parser.error("Species summary is missing")
    with args.species_summary.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = {"species_prefix", *(OUTPUT_COLUMNS[output] for output in required)}
        if not columns.issubset(reader.fieldnames or ()):
            parser.error("Species summary lacks required output validation columns")
        rows = list(reader)
    if not rows or (args.expected_task_count and len(rows) != args.expected_task_count):
        parser.error("Species summary row count does not match the required species set")
    failures = []
    for output in required:
        missing = []
        validator = has_gff_feature if output == "gff" else has_fasta_sequence
        for row in rows:
            species = str(row.get("species_prefix") or "").strip()
            raw_path = str(row.get(OUTPUT_COLUMNS[output]) or "").strip()
            if not species or not raw_path or not validator(Path(raw_path)):
                missing.append(species or "<unknown>")
        if missing:
            failures.append("Required formatted {} is missing or invalid for: {}".format(
                OUTPUT_LABELS[output], ", ".join(missing[:20])))
    if failures:
        parser.error("; ".join(failures))
    if tuple(default_required) == ("genome",) and required == ["genome"]:
        print("Required genomes verified for {} species".format(len(rows)))
    else:
        print("Required {} verified for {} species".format(",".join(required), len(rows)))


if __name__ == "__main__":
    main()
