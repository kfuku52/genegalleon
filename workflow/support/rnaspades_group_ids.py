#!/usr/bin/env python3
"""Translate rnaSPAdes group assignments to GeneGalleon's public FASTA IDs."""
import argparse
import re
import sys

from cdskit.aggregate import read_group_table
from cdskit.regroup import transcript_key
from cdskit.tsvio import read_tsv_columns, write_tsv


def normalize(table, audit, species, output):
    if ("split_by_locus",) not in read_tsv_columns(audit, ("decision",)):
        return False

    def normalized_rows():
        for raw, group in read_group_table(table).items():
            old, isoform = transcript_key(raw)
            yield dict(transcript_id=f"{species}_{old}-{isoform}", group_id=f"{species}_{group}")

    write_tsv(output, normalized_rows(), ("transcript_id", "group_id"))
    return True


def rename_selected(table, fasta):
    groups = read_group_table(table)
    seen = set()
    with open(fasta, encoding="utf-8") as stream:
        for line in stream:
            if line.startswith(">"):
                group = groups[line[1:].split()[0]]
                if group in seen:
                    raise ValueError(f"More than one representative for {group}")
                seen.add(group)
                line = f">{group}\n"
            sys.stdout.write(line)


def restore_split_candidates(table, original):
    groups = read_group_table(table)
    keep = {name for name, group in groups.items() if re.search(r"\.split[0-9]+$", group)}
    for line in sys.stdin:
        if line.startswith(">"):
            keep.discard(line[1:].split()[0])
        sys.stdout.write(line)
    if not keep:
        return
    # The workflow sorts by name after padding. Restored records can be added
    # here without changing the order used for representative selection.
    with open(original, encoding="utf-8") as stream:
        selected = False
        for line in stream:
            if line.startswith(">"):
                selected = line[1:].split()[0] in keep
            if selected:
                sys.stdout.write(line)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    normalized = commands.add_parser("normalize")
    normalized.add_argument("groups")
    normalized.add_argument("audit")
    normalized.add_argument("species")
    normalized.add_argument("output")
    selected = commands.add_parser("rename-selected")
    selected.add_argument("groups")
    selected.add_argument("fasta")
    restore = commands.add_parser("restore-split-candidates")
    restore.add_argument("groups")
    restore.add_argument("original")
    args = parser.parse_args()
    if args.command == "normalize":
        print("yes" if normalize(args.groups, args.audit, args.species, args.output) else "no")
    elif args.command == "rename-selected":
        rename_selected(args.groups, args.fasta)
    else:
        restore_split_candidates(args.groups, args.original)


if __name__ == "__main__":
    main()
