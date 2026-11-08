#!/usr/bin/env python3
"""Validate the biological sequence contract of a FASTA artifact."""

from __future__ import annotations

import argparse
from pathlib import Path

from cdskit.fasta_contract import DNA_ALPHABET as DNA_ALPHABET
from cdskit.fasta_contract import PROTEIN_ALPHABET as PROTEIN_ALPHABET
from cdskit.fasta_contract import FastaSummary as FastaSummary
from cdskit.fasta_contract import SequenceContractError as SequenceContractError
from cdskit.fasta_contract import validate_fasta as _validate_fasta
from cdskit.fasta_contract import validate_fasta_stream as validate_fasta_stream


def validate_fasta(path: Path, expected: str) -> FastaSummary:
    if path.is_symlink():
        raise SequenceContractError(f"FASTA contract path must not be a symlink: {path}")
    return _validate_fasta(path, expected)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--expected", required=True, choices=("dna", "codon", "protein"))
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        summary = validate_fasta(args.input, args.expected)
    except (OSError, UnicodeError, SequenceContractError) as exc:
        print(f"FASTA sequence contract error: {exc}")
        return 1
    print(
        "FASTA sequence contract is valid: "
        f"type={args.expected} records={summary.records} residues={summary.residues} "
        f"aligned_length={summary.aligned_length} path={args.input}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
