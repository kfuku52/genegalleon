#!/usr/bin/env python3
"""Remove trimmed protein rows that contain only IQ-TREE missing data."""

import argparse
import re
from pathlib import Path

from Bio.SeqIO.FastaIO import SimpleFastaParser

# B, Z, and J are partial ambiguities. X, U, O, and gaps are missing data.
IQTREE_PROTEIN_STATE = re.compile(r"[ACDEFGHIKLMNPQRSTVWYBZJ]", re.ASCII | re.IGNORECASE)


def filter_alignment(source, output, excluded):
    kept = removed = 0
    with source.open() as fasta, output.open("w") as alignment, excluded.open("w") as audit:
        audit.write("gene_id\treason\n")
        for title, sequence in SimpleFastaParser(fasta):
            if IQTREE_PROTEIN_STATE.search(sequence):
                alignment.write(f">{title}\n{sequence}\n")
                kept += 1
            else:
                audit.write(f"{title.split()[0]}\tno_usable_residues_after_trimming\n")
                removed += 1
    return kept, removed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("excluded", type=Path)
    args = parser.parse_args()
    kept, removed = filter_alignment(args.source, args.output, args.excluded)
    print(f"Trimmed protein alignment: {kept} rows retained, {removed} rows with no usable residues removed.")


if __name__ == "__main__":
    main()
