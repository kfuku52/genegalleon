"""Add provisionally supported ORFs that are absent from the representative catalogue."""

import csv
import subprocess

from Bio.SeqIO.FastaIO import SimpleFastaParser
from cdskit.catalog_rescue import (
    CATALOG_DIAMOND_COLUMNS,
    RESCUE_COLUMNS,
    RescueCandidate,
    decide_catalog_rescues,
    read_catalog_hits,
)
from cdskit.tsvio import write_tsv
from transcriptome_homology_split_review import SUPPORTED_PREVIEW_DECISIONS, preview_selections
from transcriptome_orf_translation import make_orf_peptide_translator


def write_proteins(path, proteins):
    with path.open("w") as stream:
        for identifier, sequence in proteins.items():
            stream.write(f">{identifier}\n{sequence}\n")


def catalogue_hits(root, queries, catalogue, cpus):
    if not queries or not catalogue:
        return {}
    query, subject = root / "rescue-candidates.faa", root / "representatives.faa"
    database, output = root / "representatives.dmnd", root / "rescue-hits.tsv"
    write_proteins(query, queries)
    write_proteins(subject, catalogue)
    subprocess.run(
        ["diamond", "makedb", "--in", str(subject), "--db", str(database), "--threads", str(cpus)],
        check=True,
    )
    subprocess.run(
        ["diamond", "blastp", "--db", str(database), "--query", str(query), "--out", str(output),
         "--threads", str(cpus), "--sensitive", "--evalue", "1e-5",
         "--max-target-seqs", "0", "--max-hsps", "1", "--outfmt", "6", *CATALOG_DIAMOND_COLUMNS],
        check=True,
    )
    return read_catalog_hits(output, queries, catalogue)


def rescue_representatives(args, review, hits, root, source_name, alternates):
    with args.report.open(newline="") as stream:
        source_by_group = {
            row["group_id"]: source_name(row["selected_id"])
            for row in csv.DictReader(stream, delimiter="\t")
        }
    group_by_source = {source: group for group, source in source_by_group.items()}
    if len(group_by_source) != len(source_by_group):
        raise ValueError("Published CDS representatives must use different source transcripts")

    peptide = make_orf_peptide_translator(args.genetic_code)
    candidates, sequences = [], {}
    for row, selected in preview_selections(args, *review, hits, source_name, alternates):
        if row["decision"] in SUPPORTED_PREVIEW_DECISIONS:
            for child, record in selected.items():
                sequence = str(record.seq)
                candidates.append(RescueCandidate(
                    child, row["original_group_id"], record.id,
                    source_name(record.id), peptide(sequence),
                ))
                sequences[child] = sequence
    decisions = []
    if candidates:
        catalogue = {}
        with args.output.open() as stream:
            for title, nucleotide in SimpleFastaParser(stream):
                group = group_by_source[source_name(title.split()[0])]
                protein = peptide(nucleotide)
                if protein and "*" not in protein:
                    catalogue[group] = protein
        exact = set(catalogue.values())
        queries = {
            c.name: c.peptide for c in candidates
            if c.peptide and "*" not in c.peptide and c.peptide not in exact
            and c.parent in source_by_group and c.transcript not in group_by_source
        }
        alignments = catalogue_hits(root, queries, catalogue, args.cpus)
        decisions = decide_catalog_rescues(candidates, catalogue, source_by_group, alignments)
        with args.output.open("a") as stream:
            for row in decisions:
                if row["decision"] == "add":
                    child, source = row["candidate_cds_id"], row["transcript_id"]
                    stream.write(f">{row['selected_orf_id']}\n{sequences[child]}\n")
                    group_by_source[source] = child
    write_tsv(args.rescue_report, decisions, RESCUE_COLUMNS)
    write_tsv(
        args.rescue_groups,
        (dict(transcript_id=source, group_id=group) for source, group in group_by_source.items()),
        ("transcript_id", "group_id"),
    )
