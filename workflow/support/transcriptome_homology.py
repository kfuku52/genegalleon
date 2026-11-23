#!/usr/bin/env python3
"""Search transcript ORF candidates and supply protein hits to CDSKIT."""

import argparse
import csv
import re
import subprocess
import tempfile
from contextlib import nullcontext
from heapq import merge
from pathlib import Path

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from cdskit.atomicio import atomic_output_paths
from cdskit.homology import HIT_COLUMNS, hit_values, parse_hit, qualifies
from cdskit.longestcds import ALT_ORF_SUFFIX
from fasta_sequence_store import open_text
from transcriptome_orf_translation import make_orf_peptide_translator

DIAMOND_COLUMNS = ("qseqid", "sseqid", "evalue", "bitscore",
                   "qstart", "qend", "qlen", "sstart", "send", "slen")
ALT_PATTERN = re.compile(ALT_ORF_SUFFIX)
TRINITY_TRANSCRIPT_ID = re.compile(r"(?:TRINITY_)?(DN\d+_c\d+_g\d+)_(i\d+)")


def fasta(path):
    with Path(path).open() as stream:
        yield from SimpleFastaParser(stream)


def prepare(baseline, alternates, peptides, code, extra_baseline=None):
    """Keep baseline order and search each distinct translated sequence once."""
    baseline_lengths = {}
    seen = set()
    query_by_peptide = {}
    aliases = []
    translate_orf = make_orf_peptide_translator(code)
    with peptides.open("w") as aa_out:
        def emit(identifier, sequence):
            if identifier in seen:
                raise ValueError(f"Duplicate CDSKIT candidate ID: {identifier}")
            seen.add(identifier)
            if not sequence:
                return
            if len(sequence) % 3:
                raise ValueError(f"CDSKIT candidate is not in codon units: {identifier}")
            peptide = translate_orf(sequence)
            if not peptide or "*" in peptide:
                return
            index = query_by_peptide.get(peptide)
            if index is None:
                index = len(aliases)
                query_by_peptide[peptide] = index
                aa_out.write(f">orf{index + 1}\n{peptide}\n")
                aliases.append([])
            aliases[index].append(identifier)

        for path in (baseline, extra_baseline):
            if not path:
                continue
            for header, sequence in fasta(path):
                identifier = header.split()[0]
                if path == extra_baseline and identifier in baseline_lengths:
                    continue
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


def retain_supported_longer_orfs(raw_hits, aliases, candidates, output, report):
    """Check provisional short replacements against local evidence on the same target."""
    with report.open(newline="") as stream:
        comparisons = [row for row in csv.DictReader(stream, delimiter="\t")
                       if row["reason"] == "supported_alternative"]
    if not comparisons:
        return

    requested = {}
    for row in comparisons:
        for candidate in (row["baseline_id"], row["selected_id"]):
            requested[candidate] = row["selected_target_id"]
    query_ids = {
        candidate: f"orf{index}"
        for index, identifiers in enumerate(aliases, 1)
        for candidate in identifiers
        if candidate in requested
    }
    targets_by_query = {}
    for candidate, query in query_ids.items():
        targets_by_query.setdefault(query, set()).add(requested[candidate])
    evidence = {}
    with raw_hits.open(newline="") as stream:
        # expand_hits already validated the DIAMOND TSV. Parse only alignments
        # that can inform a provisional replacement.
        for line in stream:
            query, _, rest = line.partition("\t")
            targets = targets_by_query.get(query)
            if targets:
                target, _, _ = rest.partition("\t")
                if target in targets:
                    row = line.rstrip("\r\n").split("\t")
                    _, hit = parse_hit(row)
                    evidence[query, target] = int(row[6]), hit

    retained = {}
    for row in comparisons:
        baseline, selected, target = row["baseline_id"], row["selected_id"], row["selected_target_id"]
        baseline_evidence = evidence.get((query_ids.get(baseline), target))
        selected_evidence = evidence.get((query_ids.get(selected), target))
        if not baseline_evidence or not selected_evidence:
            continue
        baseline_length, before = baseline_evidence
        selected_length, after = selected_evidence
        before_support = (before.bitscore, before.target_coverage)
        after_support = (after.bitscore, after.target_coverage)
        if (baseline_length > selected_length and before_support != after_support
                and all(old >= new for old, new in zip(before_support, after_support, strict=True))):
            retained[selected] = (baseline, before)
    if not retained:
        return

    baseline_ids = {baseline for baseline, _ in retained.values()}
    original = {record.id: record for record in SeqIO.parse(candidates, "fasta") if record.id in baseline_ids}
    fields = ("target_id", "evalue", "bitscore", "query_coverage", "target_coverage")
    with atomic_output_paths([output, report]) as (staged_fasta, staged_report):
        with output.open() as stream:
            SeqIO.write((original[retained[r.id][0]] if r.id in retained else r
                         for r in SeqIO.parse(stream, "fasta")), staged_fasta, "fasta")
        with report.open(newline="") as source, open(staged_report, "w", newline="") as target:
            rows = csv.DictReader(source, delimiter="\t")
            writer = csv.DictWriter(target, fieldnames=rows.fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                if row["selected_id"] in retained:
                    baseline, hit = retained[row["selected_id"]]
                    row.update(selected_id=baseline, selected_description=original[baseline].description,
                               reason="retain_dominant_baseline_hit")
                    for name, value in zip(fields, hit_values(hit), strict=True):
                        row[f"baseline_{name}"] = row[f"selected_{name}"] = value
                writer.writerow(row)


def write_eligible_candidates(
    baseline, alternates, supported, output, restored=(), sources=None, preview_alternates=None,
):
    """Remove alternatives that cannot beat their length-first baseline."""
    preview_context = preview_alternates.open("w") if preview_alternates else nullcontext()
    with output.open("w") as destination, preview_context as preview:
        records = merge(fasta(baseline), restored) if restored else fasta(baseline)
        for header, sequence in records:
            destination.write(f">{header}\n{sequence}\n")
        for header, sequence in fasta(alternates):
            identifier = header.split()[0]
            if identifier in supported:
                record = f">{header}\n{sequence}\n"
                if preview is not None:
                    preview.write(record)
                if sources is None or source_name(identifier) in sources:
                    destination.write(record)


def source_name(identifier):
    return ALT_PATTERN.sub("", identifier)


def target_links(candidates, target_fasta, root, cpus):
    accessions = {
        target
        for _, _, supports in candidates
        for per_orf in supports.values()
        for targets in per_orf
        for target in targets
    }
    if not accessions:
        return
    targets, database, alignments = (
        root / "split-targets.faa", root / "split-targets.dmnd", root / "split-target-hits.tsv",
    )
    with targets.open("w") as stream:
        for header, sequence in fasta(target_fasta):
            accession = header.split()[0]
            if accession in accessions:
                stream.write(f">{accession}\n{sequence}\n")
                accessions.remove(accession)
                if not accessions:
                    break
    if accessions:
        raise ValueError(f"Swiss-Prot FASTA is missing DIAMOND target: {min(accessions)}")
    subprocess.run(
        ["diamond", "makedb", "--in", str(targets), "--db", str(database), "--threads", str(cpus)],
        check=True,
    )
    subprocess.run(
        ["diamond", "blastp", "--db", str(database), "--query", str(targets),
         "--out", str(alignments), "--threads", str(cpus), "--sensitive",
         "--evalue", "1e-5", "--max-target-seqs", "0", "--max-hsps", "1",
         "--outfmt", "6", *DIAMOND_COLUMNS],
        check=True,
    )
    with alignments.open(newline="") as stream:
        for row in csv.reader(stream, delimiter="\t"):
            query, hit = parse_hit(row)
            if query != hit.target and qualifies(hit):
                yield query, hit.target


def split_input_transcripts(path, assembler):
    if assembler == "rnaspades":
        from cdskit.regroup import transcript_groups, transcript_key

        for group, members in transcript_groups(path).items():
            for raw in members:
                _, isoform = transcript_key(raw)
                yield raw, group, isoform
        return

    # Match the IDs emitted by the Trinity branch of longest-CDS extraction.
    seen = set()
    with open_text(path) as stream:
        for line in stream:
            if not line.startswith(">"):
                continue
            raw, *_ = line[1:].split(maxsplit=1) or [""]
            match = TRINITY_TRANSCRIPT_ID.fullmatch(raw)
            if not match:
                raise ValueError(f"Invalid Trinity transcript ID: {raw!r}")
            group, isoform = match[1].replace("_", "-"), match[2]
            if (group, isoform) in seen:
                raise ValueError(f"Duplicate Trinity transcript key: {group}/{isoform}")
            seen.add((group, isoform))
            yield raw, group, isoform
    if not seen:
        raise ValueError(f"No Trinity transcripts in {path}")


def split_groups(args, hits, root):
    from cdskit.homology_split import (
        SPLIT_COLUMNS,
        connected_targets,
        plan_splits,
        potential_splits,
        read_split_targets,
    )
    from cdskit.regroup import read_loci
    from cdskit.tsvio import write_tsv

    groups, raw_to_source = {}, {}
    for raw, group, isoform in split_input_transcripts(args.split_transcripts, args.split_assembler):
        normalized = f"{args.split_species}_{group}"
        source = f"{normalized}-{isoform}"
        groups.setdefault(normalized, []).append(source)
        raw_to_source[raw] = source
    loci = {}
    if args.split_loci:
        loci = {raw_to_source[raw]: locus for raw, locus in read_loci(args.split_loci, raw_to_source).items()}
    multi = {member for members in groups.values() if len(members) > 1 for member in members}
    support = read_split_targets(hits, source_name, multi)
    candidates = list(potential_splits(groups, support))
    links = list(target_links(candidates, args.split_targets, root, args.cpus))
    assignments, review = plan_splits(candidates, links, loci, apply=args.split_mode == "apply")
    if args.split_report:
        write_tsv(args.split_report, review, SPLIT_COLUMNS)
    if args.split_groups:
        write_tsv(args.split_groups, (
            dict(transcript_id=member, group_id=assignments.get(member, group))
            for group, members in groups.items() for member in members
        ), ("transcript_id", "group_id"))
    return assignments, (groups, review, connected_targets(links), loci)


def restrict_hits(hits, output, sources):
    with hits.open(newline="") as stream, output.open("w", newline="") as destination:
        reader = csv.reader(stream, delimiter="\t")
        writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
        writer.writerow(next(reader))
        writer.writerows(row for row in reader if source_name(row[0]) in sources)


def aggregate(candidates, hits, expressions, output, report, cpus, group_table=None,
              orf_transcripts=None, genetic_code=1):
    group_args = ["--group_table", str(group_table)] if group_table else []
    expression_args = ["--expression", *expressions] if expressions else []
    orf_args = ["--orf_transcripts", str(orf_transcripts), "--codon_table", str(genetic_code)] if orf_transcripts else []
    subprocess.run(
        ["cdskit", "aggregate", "--seq_file", str(candidates), "--out_file", str(output),
         "--threads", str(cpus), *expression_args, *group_args, *orf_args, "--hit_table", str(hits),
         "--selection_report", str(report)],
        check=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--alternates", type=Path)
    parser.add_argument("--orf-transcripts", type=Path)
    parser.add_argument("--database", required=True, type=Path)
    grouping = parser.add_mutually_exclusive_group(required=True)
    grouping.add_argument("--expression")
    grouping.add_argument("--group-table", type=Path)
    parser.add_argument("--genetic-code", type=int, required=True)
    parser.add_argument("--cpus", type=int, required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--split-mode", choices=("off", "report", "apply"), default="off")
    parser.add_argument("--split-transcripts", type=Path)
    parser.add_argument("--split-assembler", choices=("rnaspades", "trinity"), default="rnaspades")
    parser.add_argument("--split-species")
    parser.add_argument("--split-loci", type=Path)
    parser.add_argument("--split-targets", type=Path)
    parser.add_argument("--split-extra-baseline", type=Path)
    parser.add_argument("--split-report", type=Path)
    parser.add_argument("--split-groups", type=Path)
    parser.add_argument("--split-preview", type=Path)
    parser.add_argument("--rescue-report", type=Path)
    parser.add_argument("--rescue-groups", type=Path)
    args = parser.parse_args()
    split = args.split_mode != "off"
    rescue = bool(args.rescue_report)
    if rescue != bool(args.rescue_groups) or (rescue and split):
        parser.error("Catalogue rescue requires both rescue outputs and --split-mode off")
    if args.split_preview and args.split_mode != "report":
        parser.error("--split-preview requires --split-mode report")
    if split and not (args.split_report and args.split_groups):
        parser.error("Split review requires a group table and an audit output")
    if args.split_assembler == "trinity" and (split or args.split_loci):
        parser.error("Homology splitting and locus-table grouping require rnaSPAdes")
    if (split or rescue) and (args.group_table or not all((
        args.alternates, args.split_transcripts, args.split_species, args.split_targets,
        args.split_extra_baseline,
    ))):
        parser.error("Split review and catalogue rescue require assembly transcripts, alternate ORFs, and their inputs and outputs")
    with tempfile.TemporaryDirectory(prefix="orf-homology.", dir=args.output.parent) as temporary:
        root = Path(temporary)
        peptides = root / "queries.faa"
        raw, hits = root / "diamond.tsv", root / "hits.tsv"
        aliases = prepare(
            args.baseline, args.alternates, peptides, args.genetic_code,
            args.split_extra_baseline if split or rescue else None,
        )
        search(args.database, peptides, raw, args.cpus)
        supported = expand_hits(raw, hits, aliases)
        candidates = args.baseline
        expressions = [args.expression] if args.expression else []
        group_table = args.group_table
        restored, sources = [], None
        preview, all_hits = None, hits
        if split or rescue:
            assignments, preview = split_groups(args, hits, root)
            sources = {header.split()[0] for header, _ in fasta(args.baseline)}
            if assignments:
                restored = sorted(
                    (header, sequence) for header, sequence in fasta(args.split_extra_baseline)
                    if (identifier := header.split()[0]) in assignments and identifier not in sources
                )
                sources.update(header.split()[0] for header, _ in restored)
            filtered = root / "selection-hits.tsv"
            restrict_hits(hits, filtered, sources)
            hits = filtered
            if split:
                expressions = []
                group_table = args.split_groups
        preview_alternates = root / "supported-alternates.fa" if args.split_preview or rescue else None
        if args.alternates:
            candidates = root / "candidates.fna"
            write_eligible_candidates(
                args.baseline, args.alternates, supported, candidates, restored, sources, preview_alternates,
            )
            expressions.insert(0, ALT_ORF_SUFFIX)
        aggregate(candidates, hits, expressions, args.output, args.report, args.cpus, group_table,
                  args.orf_transcripts, args.genetic_code)
        retain_supported_longer_orfs(raw, aliases, candidates, args.output, args.report)
        del aliases
        if args.split_preview:
            from transcriptome_homology_split_review import write_split_preview

            write_split_preview(args, preview, all_hits, source_name, preview_alternates)
        if rescue:
            from transcriptome_catalog_rescue import rescue_representatives

            rescue_representatives(args, preview, all_hits, root, source_name, preview_alternates)


if __name__ == "__main__":
    main()
