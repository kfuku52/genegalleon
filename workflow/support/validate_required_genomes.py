#!/usr/bin/env python3
"""Verify that each formatted species has a nonempty genome FASTA."""

import argparse
import csv
import gzip
from pathlib import Path


def has_fasta_sequence(path):
    if not path.is_file() or path.stat().st_size == 0:
        return False
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8") as handle:
            first = handle.readline()
            sequence = handle.readline().strip()
            return first.startswith(">") and bool(sequence) and not sequence.startswith(">")
    except (OSError, UnicodeError):
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-summary", required=True, type=Path)
    parser.add_argument("--expected-task-count", type=int, default=0)
    args = parser.parse_args()
    if args.expected_task_count < 0:
        parser.error("--expected-task-count must be nonnegative")
    if not args.species_summary.is_file():
        parser.error("Species summary is missing")
    with args.species_summary.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"species_prefix", "genome_output_path"}.issubset(reader.fieldnames or ()):
            parser.error("Species summary lacks genome validation columns")
        rows = list(reader)
    if not rows or (args.expected_task_count and len(rows) != args.expected_task_count):
        parser.error("Species summary row count does not match the required species set")
    missing = []
    for row in rows:
        species = str(row.get("species_prefix") or "").strip()
        raw_path = str(row.get("genome_output_path") or "").strip()
        if not species or not raw_path:
            missing.append(species or "<unknown>")
            continue
        path = Path(raw_path)
        if not has_fasta_sequence(path):
            missing.append(species)
    if missing:
        parser.error("Required formatted genome is missing or invalid for: " + ", ".join(missing[:20]))
    print("Required genomes verified for {} species".format(len(rows)))


if __name__ == "__main__":
    main()
