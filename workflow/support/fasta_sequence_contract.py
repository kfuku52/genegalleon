#!/usr/bin/env python3
"""Validate the biological sequence contract of a FASTA artifact."""

from __future__ import annotations

import argparse
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

DNA_ALPHABET = frozenset("ACGTURYSWKMBDHVN-?.")
PROTEIN_ALPHABET = frozenset("ABCDEFGHIKLMNPQRSTVWXYZJUO*-?.")


class SequenceContractError(RuntimeError):
    """Raised when a FASTA artifact violates its declared sequence type."""


@dataclass(frozen=True)
class FastaSummary:
    records: int
    residues: int
    aligned_length: int


def _open_text(path: Path) -> TextIO:
    with path.open("rb") as handle:
        magic = handle.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="strict")
    return path.open("r", encoding="utf-8", errors="strict")


def validate_fasta(path: Path, expected: str) -> FastaSummary:
    """Validate one plain or gzip-compressed FASTA file.

    ``codon`` requires an aligned IUPAC nucleotide FASTA whose aligned and
    ungapped sequence lengths are multiples of three. ``dna`` validates only
    the alphabet, while ``protein`` accepts the IUPAC amino-acid alphabet.
    """

    if expected not in {"dna", "codon", "protein"}:
        raise SequenceContractError(f"Unsupported FASTA sequence type: {expected}")
    if path.is_symlink():
        raise SequenceContractError(f"FASTA contract path must not be a symlink: {path}")
    if not path.is_file():
        raise FileNotFoundError(path)

    allowed = PROTEIN_ALPHABET if expected == "protein" else DNA_ALPHABET
    record_lengths: list[int] = []
    ungapped_lengths: list[int] = []
    current_length: int | None = None
    current_ungapped = 0
    invalid: set[str] = set()
    residues = 0

    with _open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if len(line) == 1:
                    raise SequenceContractError(f"Empty FASTA header at line {line_number}: {path}")
                if current_length is not None:
                    record_lengths.append(current_length)
                    ungapped_lengths.append(current_ungapped)
                current_length = 0
                current_ungapped = 0
                continue
            if current_length is None:
                raise SequenceContractError(
                    f"FASTA sequence data appears before the first header at line {line_number}: {path}"
                )
            sequence = "".join(line.split()).upper()
            invalid.update(set(sequence) - allowed)
            current_length += len(sequence)
            current_ungapped += sum(character not in "-?." for character in sequence)
            residues += len(sequence)

    if current_length is not None:
        record_lengths.append(current_length)
        ungapped_lengths.append(current_ungapped)
    if not record_lengths:
        raise SequenceContractError(f"FASTA contains no records: {path}")
    if not residues or any(length == 0 for length in record_lengths):
        raise SequenceContractError(f"FASTA contains an empty sequence: {path}")
    if invalid:
        rendered = "".join(sorted(invalid))
        raise SequenceContractError(f"FASTA contains characters outside the {expected} alphabet: {rendered} ({path})")

    aligned_length = record_lengths[0]
    if expected == "codon":
        if len(set(record_lengths)) != 1:
            raise SequenceContractError(f"Codon FASTA is not an alignment: {path}")
        if aligned_length % 3:
            raise SequenceContractError(
                f"Codon FASTA aligned length is not a multiple of three: {aligned_length} ({path})"
            )
        invalid_ungapped = [length for length in ungapped_lengths if length % 3]
        if invalid_ungapped:
            raise SequenceContractError(
                f"Codon FASTA has ungapped sequence lengths that are not multiples of three: {path}"
            )

    return FastaSummary(
        records=len(record_lengths),
        residues=residues,
        aligned_length=aligned_length,
    )


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
