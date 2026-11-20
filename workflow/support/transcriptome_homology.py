#!/usr/bin/env python3
"""Search transcript ORF candidates and supply protein hits to CDSKIT."""

import argparse
import csv
import re
import subprocess
import tempfile
from pathlib import Path

from Bio.Data.CodonTable import ambiguous_generic_by_id
from Bio.Seq import translate
from Bio.SeqIO.FastaIO import SimpleFastaParser
from cdskit.homology import HIT_COLUMNS, parse_hit, qualifies
from cdskit.longestcds import ALT_ORF_SUFFIX

DIAMOND_COLUMNS = ("qseqid", "sseqid", "evalue", "bitscore",
                   "qstart", "qend", "qlen", "sstart", "send", "slen")
ALT_PATTERN = re.compile(ALT_ORF_SUFFIX)


def fasta(path):
    with Path(path).open() as stream:
        yield from SimpleFastaParser(stream)


def prepare(baseline, alternates, peptides, code):
    """Keep baseline order and search each distinct translated sequence once."""
    baseline_lengths = {}
    seen = set()
    query_by_peptide = {}
    aliases = []
    codon_table = ambiguous_generic_by_id[code]
    with peptides.open("w") as aa_out:
        def emit(identifier, sequence):
            if identifier in seen:
                raise ValueError(f"Duplicate CDSKIT candidate ID: {identifier}")
            seen.add(identifier)
            if not sequence:
                return
            if len(sequence) % 3:
                raise ValueError(f"CDSKIT candidate is not in codon units: {identifier}")
            peptide = translate(sequence, table=codon_table, gap="-").removesuffix("*")
            if not peptide or "*" in peptide:
                return
            index = query_by_peptide.get(peptide)
            if index is None:
                index = len(aliases)
                query_by_peptide[peptide] = index
                aa_out.write(f">orf{index + 1}\n{peptide}\n")
                aliases.append([])
            aliases[index].append(identifier)

        for header, sequence in fasta(baseline):
            identifier = header.split()[0]
            if ALT_PATTERN.search(identifier):
                raise ValueError(f"Baseline ID uses a reserved alternate-ORF suffix: {identifier}")
            baseline_lengths[identifier] = len(sequence)
            emit(identifier, sequence)
        if alternates:
            for header, sequence in fasta(alternates):
                identifier = header.split()[0]
                source, count = ALT_PATTERN.subn("", identifier)
                if count != 1:
                    raise ValueError(f"Invalid CDSKIT alternate ORF ID: {identifier}")
                if source in baseline_lengths:
                    if len(sequence) > baseline_lengths[source]:
                        raise ValueError(f"Alternative is longer than its longest ORF: {identifier}")
                    emit(identifier, sequence)
    if not baseline_lengths:
        raise ValueError("No baseline CDSKIT candidates are available")
    return aliases


def search(database, peptides, hits, cpus):
    if not peptides.stat().st_size:
        hits.touch()
        return
    subprocess.run(
        ["diamond", "blastp", "--db", str(database), "--query", str(peptides), "--out", str(hits),
         "--threads", str(cpus), "--fast", "--evalue", "1e-5", "--max-target-seqs", "5",
         "--max-hsps", "1", "--outfmt", "6", *DIAMOND_COLUMNS],
        check=True,
    )


def expand_hits(raw, output, aliases):
    """Qualify each alignment once before expanding it to candidate IDs."""
    supported = set()
    with raw.open(newline="") as stream, output.open("w", newline="") as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(HIT_COLUMNS)
        for row in csv.reader(stream, delimiter="\t"):
            query, hit = parse_hit(row)
            index = int(query[3:]) - 1
            if not 0 <= index < len(aliases) or query != f"orf{index + 1}":
                raise ValueError("DIAMOND returned an invalid hit row or query ID")
            if qualifies(hit):
                for identifier in aliases[index]:
                    supported.add(identifier)
                    writer.writerow((identifier, *row[1:]))
    return supported


def write_eligible_candidates(baseline, alternates, supported, output):
    """Remove alternatives that cannot beat their length-first baseline."""
    with output.open("w") as destination:
        for header, sequence in fasta(baseline):
            destination.write(f">{header}\n{sequence}\n")
        for header, sequence in fasta(alternates):
            if header.split()[0] in supported:
                destination.write(f">{header}\n{sequence}\n")


def aggregate(candidates, hits, expressions, output, report, cpus, group_table=None):
    group_args = ["--group_table", str(group_table)] if group_table else []
    expression_args = ["--expression", *expressions] if expressions else []
    subprocess.run(
        ["cdskit", "aggregate", "--seq_file", str(candidates), "--out_file", str(output),
         "--threads", str(cpus), *expression_args, *group_args, "--hit_table", str(hits),
         "--selection_report", str(report)],
        check=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--alternates", type=Path)
    parser.add_argument("--database", required=True, type=Path)
    grouping = parser.add_mutually_exclusive_group(required=True)
    grouping.add_argument("--expression")
    grouping.add_argument("--group-table", type=Path)
    parser.add_argument("--genetic-code", type=int, required=True)
    parser.add_argument("--cpus", type=int, required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="orf-homology.", dir=args.output.parent) as temporary:
        root = Path(temporary)
        peptides = root / "queries.faa"
        raw, hits = root / "diamond.tsv", root / "hits.tsv"
        aliases = prepare(args.baseline, args.alternates, peptides, args.genetic_code)
        search(args.database, peptides, raw, args.cpus)
        supported = expand_hits(raw, hits, aliases)
        del aliases
        candidates = args.baseline
        expressions = [args.expression] if args.expression else []
        if args.alternates:
            candidates = root / "candidates.fna"
            write_eligible_candidates(args.baseline, args.alternates, supported, candidates)
            expressions.insert(0, ALT_ORF_SUFFIX)
        aggregate(candidates, hits, expressions, args.output, args.report, args.cpus, args.group_table)


if __name__ == "__main__":
    main()
