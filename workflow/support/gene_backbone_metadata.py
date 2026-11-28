#!/usr/bin/env python3
"""Adapt gene-family metadata to and from the CDSKIT backbone CLI."""

import argparse
import csv
from pathlib import Path

from Bio.SeqIO.FastaIO import SimpleFastaParser
from cdskit.sequence_io import open_text
from nwkit.species_parser import get_species_parser

FAMILY_PRIORITY = {"target_supported": 0, "competing_families": 2}


def prepare(args):
    labels = get_species_parser(args)
    statuses = {}
    if args.family_screen:
        with open_text(args.family_screen) as handle:
            rows = csv.DictReader(handle, delimiter="\t")
            required = {"gene_id", "status", "retained"}
            if not required <= set(rows.fieldnames or []):
                raise ValueError(f"Missing family-screen columns: {args.family_screen}")
            statuses = {row["gene_id"]: row["status"] for row in rows if row["retained"] == "1"}
    with open_text(args.fasta) as source, args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("gene_id", "species", "priority", "family_screen"))
        for title, _ in SimpleFastaParser(source):
            gene = title.split(None, 1)[0]
            status = statuses.get(gene, "")
            writer.writerow((gene, labels.parse(gene).species_label or "", FAMILY_PRIORITY.get(status, 1), status))


def report(args):
    with open_text(args.metadata) as handle:
        status = {row["gene_id"]: row["family_screen"] for row in csv.DictReader(handle, delimiter="\t")}
    with open_text(args.report) as handle, args.output.open("w", encoding="utf-8", newline="") as out:
        rows = csv.DictReader(handle, delimiter="\t")
        fields = ["family_screen" if name == "priority" else name for name in rows.fieldnames]
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            del row["priority"]
            row["family_screen"] = status[row["gene_id"]]
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(required=True)
    metadata = subparsers.add_parser("prepare")
    metadata.add_argument("--fasta", type=Path, required=True)
    metadata.add_argument("--output", type=Path, required=True)
    metadata.add_argument("--family-screen", type=Path)
    metadata.add_argument("--species-parser", default="extended-taxonomic")
    metadata.add_argument("--species-regex")
    metadata.add_argument("--species-map-tsv")
    metadata.set_defaults(handler=prepare)
    output = subparsers.add_parser("report")
    for arg in ("metadata", "report", "output"):
        output.add_argument("--" + arg, type=Path, required=True)
    output.set_defaults(handler=report)
    args = parser.parse_args()
    try:
        args.handler(args)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, f"Gene backbone metadata failed: {exc}\n")


if __name__ == "__main__":
    main()
