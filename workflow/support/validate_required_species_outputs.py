#!/usr/bin/env python3
"""Verify the formatted outputs required by an input-generation dataset."""

import argparse
import csv
import gzip
import math
import zlib
from pathlib import Path

OUTPUT_COLUMNS = {
    "cds": "cds_output_path",
    "gff": "gff_output_path",
    "genome": "genome_output_path",
}
OUTPUT_LABELS = {"cds": "CDS", "gff": "GFF", "genome": "genome"}
DNA_RESIDUES = frozenset("ACGTURYSWKMBDHVN")
DNA_SYMBOLS = DNA_RESIDUES | frozenset("-?.")


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("rt", encoding="utf-8")


def has_fasta_sequence(path):
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        with open_text(path) as handle:
            has_record = False
            has_residue = False
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if has_record and not has_residue:
                        return False
                    if not line[1:].strip():
                        return False
                    has_record = True
                    has_residue = False
                elif not has_record:
                    return False
                else:
                    sequence = line.upper()
                    if not set(sequence) <= DNA_SYMBOLS:
                        return False
                    has_residue = has_residue or any(base in DNA_RESIDUES for base in sequence)
            return has_record and has_residue
    except (OSError, UnicodeError, EOFError, zlib.error):
        return False


def has_gff_feature(path):
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        with open_text(path) as handle:
            has_feature = False
            in_fasta = False
            for line in handle:
                if in_fasta:
                    continue
                if line.startswith("##FASTA"):
                    in_fasta = True
                    continue
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 9 or any(not field.strip() for field in fields):
                    return False
                if fields[0] == "." or fields[2] == ".":
                    return False
                try:
                    start, end = int(fields[3]), int(fields[4])
                    score = None if fields[5] == "." else float(fields[5])
                except ValueError:
                    return False
                if start < 1 or end < start or (score is not None and not math.isfinite(score)):
                    return False
                if fields[6] not in {"+", "-", ".", "?"} or fields[7] not in {"0", "1", "2", "."}:
                    return False
                has_feature = True
            return has_feature
    except (OSError, UnicodeError, EOFError, zlib.error):
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
