#!/usr/bin/env python3
"""Small metadata-table helpers for explicit SRA accession workflows."""

import argparse
import csv
import os
import stat
import sys
import uuid
from contextlib import contextmanager
from pathlib import Path


def read_accessions(accession_path):
    requested = []
    seen = set()
    with Path(accession_path).open("rt", encoding="utf-8") as handle:
        for raw_line in handle:
            accession = raw_line.strip()
            if not accession or accession in seen:
                continue
            seen.add(accession)
            requested.append(accession)
    return requested


def metadata_runs(metadata_path):
    found = set()
    metadata_path = Path(metadata_path)
    if not metadata_path.exists() or metadata_path.stat().st_size == 0:
        return found
    with metadata_path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if "run" not in fieldnames:
            return found
        for row in reader:
            run = str(row.get("run", "") or "").strip()
            if run:
                found.add(run)
    return found


def print_missing_accessions(metadata_path, accession_path):
    found = metadata_runs(metadata_path)
    for accession in read_accessions(accession_path):
        if accession not in found:
            print(accession)


def row_is_transcriptomic(row):
    lib_source = str(row.get("lib_source", "") or "").strip().lower()
    lib_strategy = str(row.get("lib_strategy", "") or "").strip().lower()
    return lib_source == "transcriptomic" or lib_strategy in {"rna-seq", "est", "clone"}


@contextmanager
def metadata_writer(output_path, fieldnames, *, extrasaction="raise"):
    """Stream to a sibling file, publishing only after a complete successful read."""
    if output_path is None:
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter="\t", lineterminator="\n",
                                extrasaction=extrasaction)
        writer.writeheader()
        yield writer
        return
    output = Path(output_path).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    pending = output.with_name("." + output.name + "." + uuid.uuid4().hex + ".tmp")
    descriptor = os.open(pending, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            if output.exists():
                os.fchmod(handle.fileno(), stat.S_IMODE(output.stat().st_mode))
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n",
                                    extrasaction=extrasaction)
            writer.writeheader()
            yield writer
        os.replace(pending, output)
    finally:
        pending.unlink(missing_ok=True)


def filter_metadata_rows(metadata_path, output_path, required_column, predicate):
    with Path(metadata_path).open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames
        if not fieldnames:
            raise SystemExit("Metadata table is missing a header: {}".format(metadata_path))
        if fieldnames.count(required_column) != 1:
            raise SystemExit("Metadata table needs one required '{}' column: {}".format(required_column, metadata_path))
        with metadata_writer(output_path, fieldnames) as writer:
            for row in reader:
                if predicate(row):
                    writer.writerow(row)


def extract_transcriptomic_rows(metadata_path, accession_path, output_path):
    requested_set = set(read_accessions(accession_path))
    filter_metadata_rows(metadata_path, output_path, "run",
                         lambda row: str(row.get("run", "") or "").strip() in requested_set
                         and row_is_transcriptomic(row))


def extract_species_rows(metadata_path, species, output_path, *, qualify_unmarked=False):
    def normalized(value):
        return " ".join(str(value or "").replace("_", " ").split()).casefold()
    requested = normalized(species)
    if not requested:
        raise ValueError("A nonempty species name is required")
    def selected(row):
        if normalized(row.get("scientific_name")) != requested:
            return False
        if (qualify_unmarked and row.get("exclusion") == "no"
                and row.get("is_sampled") == "" and row.get("is_qualified") == ""):
            row["is_sampled"] = row["is_qualified"] = "yes"
        return True

    filter_metadata_rows(metadata_path, output_path, "scientific_name", selected)


def merge_metadata_tables_by_run(primary_table, extra_table, output_path):
    fieldnames = []
    seen_fieldnames = set()
    seen_runs = set()
    sources = []

    for path in (Path(primary_table), Path(extra_table)):
        if not path.exists() or path.stat().st_size == 0:
            continue
        with path.open("rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                continue
            sources.append((path, tuple(reader.fieldnames)))
            for fieldname in reader.fieldnames:
                if fieldname not in seen_fieldnames:
                    fieldnames.append(fieldname)
                    seen_fieldnames.add(fieldname)

    if len(fieldnames) == 0:
        raise SystemExit("Could not determine metadata header from {} or {}".format(primary_table, extra_table))

    with metadata_writer(output_path, fieldnames, extrasaction="ignore") as writer:
        for path, expected_fields in sources:
            with path.open("rt", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if tuple(reader.fieldnames or ()) != expected_fields:
                    raise ValueError("Metadata header changed during merge: " + str(path))
                for row in reader:
                    run = str(row.get("run", "") or "").strip()
                    key = run or tuple((name, row.get(name, "")) for name in reader.fieldnames)
                    if key not in seen_runs:
                        seen_runs.add(key)
                        writer.writerow({name: row.get(name, "") or "" for name in fieldnames})


def build_arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    missing_parser = subparsers.add_parser("missing")
    missing_parser.add_argument("metadata_table")
    missing_parser.add_argument("accession_file")

    extract_parser = subparsers.add_parser("extract-transcriptomic")
    extract_parser.add_argument("metadata_table")
    extract_parser.add_argument("accession_file")
    extract_parser.add_argument("output")

    species_parser = subparsers.add_parser("extract-species")
    species_parser.add_argument("metadata_table")
    species_parser.add_argument("species")
    species_parser.add_argument("output", help="Output TSV, or '-' for stdout")
    species_parser.add_argument("--qualify-unmarked", action="store_true",
                                help="Set blank is_sampled/is_qualified flags to yes only when exclusion=no")

    merge_parser = subparsers.add_parser("merge-by-run")
    merge_parser.add_argument("primary_table")
    merge_parser.add_argument("extra_table")
    merge_parser.add_argument("output")

    return parser


def main(argv=None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    if args.command == "missing":
        print_missing_accessions(args.metadata_table, args.accession_file)
        return 0
    if args.command == "extract-transcriptomic":
        extract_transcriptomic_rows(args.metadata_table, args.accession_file, args.output)
        return 0
    if args.command == "extract-species":
        extract_species_rows(args.metadata_table, args.species, None if args.output == "-" else args.output,
                             qualify_unmarked=args.qualify_unmarked)
        return 0
    if args.command == "merge-by-run":
        merge_metadata_tables_by_run(args.primary_table, args.extra_table, args.output)
        return 0
    parser.error("unknown command: {}".format(args.command))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
