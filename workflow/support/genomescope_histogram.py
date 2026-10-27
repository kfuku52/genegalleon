#!/usr/bin/env python3
"""Validate count settings and preserve GenomeScope's terminal-bin convention."""

import argparse
import json
from pathlib import Path


def rows(path):
    previous = 0
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 2:
                raise ValueError(f"histogram line {number} must contain two integers")
            coverage, frequency = map(int, fields)
            if coverage <= previous or frequency < 0:
                raise ValueError("histogram counts must increase and frequencies must be nonnegative")
            previous = coverage
            yield coverage, frequency


def prepare(source, output, metadata, ceiling):
    if Path(source).resolve() == Path(output).resolve():
        raise ValueError("source and prepared histogram must differ")
    distinct = mass = maximum = boundary = 0
    for coverage, frequency in rows(source):
        distinct += frequency
        mass += coverage * frequency
        if frequency:
            maximum = coverage
        if coverage == ceiling:
            boundary = frequency
    if not distinct or maximum < 2:
        raise ValueError("histogram contains no usable spectrum")
    # GenomeScope discards its final row. End with an explicit empty bin so the
    # highest observed multiplicity contributes to both fitting and genome size.
    # Fill missing bins, retaining zero-frequency observations in the fit.
    with Path(output).open("w") as handle:
        previous = 0
        for coverage, frequency in rows(source):
            if coverage > maximum:
                break
            for missing in range(previous + 1, coverage):
                handle.write(f"{missing}\t0\n")
            handle.write(f"{coverage}\t{frequency}\n")
            previous = coverage
        handle.write(f"{maximum + 1}\t0\n")
    result = dict(distinct_kmers=distinct, total_kmer_observations=mass,
                  maximum_observed_multiplicity=maximum, counter_ceiling=ceiling,
                  ceiling_bin_distinct_kmers=boundary, terminal_zero_bin=maximum + 1)
    Path(metadata).write_text(json.dumps(result, indent=2) + "\n")
    return result


def validate(args):
    if not 1 <= args.kmer_length <= 256:
        raise ValueError("k-mer length must be between 1 and 256")
    if not 1 <= args.lower <= args.upper <= 2147483646:
        raise ValueError("require 1 <= lower <= upper <= 2147483646")
    if not 1 <= args.ploidy <= 6:
        raise ValueError("ploidy must be between 1 and 6")
    for value in (args.initial_coverage, args.max_kmercov):
        if value != -1 and value <= 0:
            raise ValueError("initial coverage and maximum coverage must be -1 or positive")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    check = commands.add_parser("validate")
    for name in ("kmer-length", "lower", "upper", "ploidy", "initial-coverage", "max-kmercov"):
        check.add_argument("--" + name, type=int, required=True)
    histogram = commands.add_parser("prepare")
    for name in ("input", "output", "metadata"):
        histogram.add_argument("--" + name, type=Path, required=True)
    histogram.add_argument("--ceiling", type=int, required=True)
    args = parser.parse_args()
    if args.command == "validate":
        validate(args)
    else:
        prepare(args.input, args.output, args.metadata, args.ceiling)


if __name__ == "__main__":
    main()
