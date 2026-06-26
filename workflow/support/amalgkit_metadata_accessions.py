#!/usr/bin/env python3
"""Small metadata-table helpers for explicit SRA accession workflows."""

import argparse
import csv
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


def extract_transcriptomic_rows(metadata_path, accession_path, output_path):
    requested_set = set(read_accessions(accession_path))
    metadata_path = Path(metadata_path)
    output_path = Path(output_path)
    with metadata_path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames
        if not fieldnames:
            raise SystemExit("Metadata table is missing a header: {}".format(metadata_path))
        if "run" not in fieldnames:
            raise SystemExit("Metadata table is missing required 'run' column: {}".format(metadata_path))

        rows = []
        for row in reader:
            run = str(row.get("run", "") or "").strip()
            if run in requested_set and row_is_transcriptomic(row):
                rows.append(row)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def merge_metadata_tables_by_run(primary_table, extra_table, output_path):
    fieldnames = []
    seen_fieldnames = set()
    rows = []
    seen_runs = set()

    for path in (Path(primary_table), Path(extra_table)):
        if not path.exists() or path.stat().st_size == 0:
            continue
        with path.open("rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                continue
            for fieldname in reader.fieldnames:
                if fieldname not in seen_fieldnames:
                    fieldnames.append(fieldname)
                    seen_fieldnames.add(fieldname)
            for row in reader:
                run = str(row.get("run", "") or "").strip()
                dedupe_key = run or tuple((name, row.get(name, "")) for name in reader.fieldnames)
                if dedupe_key in seen_runs:
                    continue
                seen_runs.add(dedupe_key)
                rows.append(row)

    if len(fieldnames) == 0:
        raise SystemExit("Could not determine metadata header from {} or {}".format(primary_table, extra_table))

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
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
    if args.command == "merge-by-run":
        merge_metadata_tables_by_run(args.primary_table, args.extra_table, args.output)
        return 0
    parser.error("unknown command: {}".format(args.command))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
