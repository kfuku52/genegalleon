#!/usr/bin/env python3
"""Select a query-domain-informed MAFFT backbone and add every remaining sequence."""

import argparse
import csv
import math
import re
import subprocess
import sys
from collections import defaultdict
from contextlib import nullcontext
from dataclasses import dataclass
from itertools import compress
from pathlib import Path
from tempfile import TemporaryDirectory

from Bio import SeqIO
from nwkit.species_parser import get_species_parser

if __package__:
    from .fasta_sequence_store import open_text
else:
    from fasta_sequence_store import open_text


RPS_FIELDS = ("qacc", "sacc", "qstart", "qend", "sstart", "send", "bitscore")
GAPS = str.maketrans("", "", "-.")
RESIDUE_RUN = re.compile(r"[^.-]+")
FAMILY_PRIORITY = {"target_supported": 0, "competing_families": 2}


@dataclass(frozen=True)
class Domain:
    model: str
    start: int
    end: int
    first: int
    last: int
    score: float


def read_fasta(path):
    sequences = {}
    with open_text(path) as handle:
        for row in SeqIO.parse(handle, "fasta"):
            sequence = str(row.seq).upper().translate(GAPS)
            if not row.id or row.id in sequences or not sequence:
                raise ValueError(f"Expected unique IDs and nonempty sequences in {path}: {row.id!r}")
            sequences[row.id] = sequence
    if not sequences:
        raise ValueError(f"No FASTA records in {path}")
    return sequences


def read_domains(path, fieldnames=None):
    domains = defaultdict(list)
    with open_text(path) as handle:
        rows = csv.DictReader(handle, fieldnames=fieldnames, delimiter="\t")
        if not set(RPS_FIELDS) <= set(rows.fieldnames or ()):
            raise ValueError(f"Missing domain search columns in {path}")
        for row in rows:
            gene, model = row["qacc"], row["sacc"]
            hits = domains[gene]  # A no-hit row is evidence that this sequence was searched.
            if model:
                start, end = sorted((int(row["qstart"]), int(row["qend"])))
                first, last = sorted((int(row["sstart"]), int(row["send"])))
                score = float(row["bitscore"])
                if min(start, first) < 1 or not math.isfinite(score):
                    raise ValueError(f"Invalid RPS-BLAST domain for {gene} in {path}")
                hits.append(Domain(model, start, end, first, last, score))
    return domains


def write_fasta(path, seqs, ids):
    with path.open("w") as handle:
        for gene in ids:
            handle.write(f">{gene}\n{seqs[gene]}\n")


def annotate_queries(queries, candidates, domains, database, threads, scratch):
    shared = {gene for gene, aa in queries.items() if candidates.get(gene) == aa}
    annotated = {gene: domains[gene] for gene in shared if gene in domains}
    missing = sorted(queries.keys() - annotated.keys())
    if missing:
        fasta, hits = scratch / "queries.fa", scratch / "query_domains.tsv"
        write_fasta(fasta, queries, missing)
        subprocess.run(
            ["rpsblast", "-query", str(fasta), "-db", str(database), "-out", str(hits),
             "-evalue", "0.01", "-outfmt", "6 " + " ".join(RPS_FIELDS), "-num_threads", str(threads)],
            check=True,
        )
        forwarded = read_domains(hits, RPS_FIELDS)
        annotated.update((gene, forwarded.get(gene, ())) for gene in missing)
    return shared, annotated


def architecture(domains):
    """Collapse overlapping Pfam alternatives while retaining distinct domain copies."""
    chosen = []
    for hit in sorted(domains, key=lambda h: (-h.score, h.start, h.end, h.model, h.first, h.last)):
        if not any(
            2 * (min(hit.end, prev.end) - max(hit.start, prev.start) + 1)
            >= min(hit.end - hit.start + 1, prev.end - prev.start + 1)
            for prev in chosen
        ):
            chosen.append(hit)
    return sorted(chosen, key=lambda h: (h.start, h.end, h.model))


def architecture_match(query, candidate):
    """Count shared profile positions in the best ordered one-to-one domain match."""
    previous = [0] * (len(candidate) + 1)
    for q in query:
        current = [0]
        for i, hit in enumerate(candidate):
            overlap = 0
            if hit.model == q.model:
                overlap = max(0, min(hit.last, q.last) - max(hit.first, q.first) + 1)
            current.append(max(previous[i + 1], current[-1], previous[i] + overlap))
        previous = current
    return previous[-1]


def family_statuses(path):
    if path is None:
        return {}
    with open_text(path) as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        if not {"gene_id", "status", "retained"} <= set(rows.fieldnames or ()):
            raise ValueError(f"Missing family-screen columns in {path}")
        return {row["gene_id"]: row["status"] for row in rows if row["retained"] == "1"}


def select_backbone(candidates, query_domains, shared, domains, statuses, species_parser, size):
    if len(shared) > size:
        raise ValueError("mafft_backbone_size is smaller than the number of matching query sequences")
    profiles = {gene: architecture(hits) for gene, hits in sorted(query_domains.items()) if hits}
    profile_lengths = {gene: sum(hit.last - hit.first + 1 for hit in hits) for gene, hits in profiles.items()}
    models = {hit.model for hits in profiles.values() for hit in hits}

    evidence = {}
    for gene in candidates:
        hits = architecture(hit for hit in domains.get(gene, ()) if hit.model in models)
        score, matched = max(
            ((architecture_match(profile, hits) / profile_lengths[ref], ref) for ref, profile in profiles.items()),
            default=(0.0, ""),
        )
        evidence[gene] = (statuses.get(gene, ""), score, matched if score else "",
                          species_parser.parse(gene).species_label or "")

    def rank(gene):
        status, score, _, species = evidence[gene]
        return FAMILY_PRIORITY.get(status, 1), -score, not species, gene

    selected = dict.fromkeys(sorted(shared), "query")
    ranked = sorted(candidates.keys() - shared, key=rank)
    represented = {evidence[gene][3] for gene in selected if evidence[gene][3]}
    # Prefer one sequence per species before filling with additional copies.
    for gene in ranked:
        if len(selected) >= size:
            break
        species = evidence[gene][3]
        if species and species not in represented:
            selected[gene] = "species_representative"
            represented.add(species)
    for gene in ranked:
        if len(selected) >= size:
            break
        selected.setdefault(gene, "additional_representative")
    return selected, evidence, profiles


def align(source, output, threads, threadit, options):
    command = ["mafft", "--amino", "--anysymbol", "--quiet", "--thread", str(threads),
               "--threadit", str(threadit), *options, str(source)]
    with output.open("w") as handle:
        subprocess.run(command, stdout=handle, check=True)


def fasta_offsets(handle):
    """Index MAFFT's wrapped FASTA by record using byte offsets."""
    starts = []
    seen = set()
    offset = 0
    for line in handle:
        if line.startswith(b">"):
            fields = line[1:].split(None, 1)
            gene = fields[0].decode() if fields else ""
            if not gene or gene in seen:
                raise ValueError(f"Expected unique, nonempty MAFFT IDs: {gene!r}")
            starts.append((gene, offset))
            seen.add(gene)
        offset += len(line)
    ends = [start for _, start in starts[1:]] + [offset] if starts else []
    return {gene: (start, end - start) for (gene, start), end in zip(starts, ends, strict=True)}


def occupied_columns(rows):
    """Mark positions with a residue in at least one backbone row."""
    width = len(next(iter(rows.values())))
    occupied = bytearray(width)
    for row in rows.values():
        if len(row) != width:
            raise ValueError("MAFFT produced an unaligned backbone row")
        for match in RESIDUE_RUN.finditer(row):
            start, end = match.span()
            occupied[start:end] = b"\1" * (end - start)
    return occupied


def write_column_map(original, final, column_map):
    old_columns = occupied_columns(original)
    new_columns = occupied_columns(final)
    if old_columns.count(1) != new_columns.count(1):
        raise ValueError("MAFFT changed the number of occupied backbone columns")
    for gene, old_row in original.items():
        if "".join(compress(old_row, old_columns)) != "".join(compress(final[gene], new_columns)):
            raise ValueError(f"MAFFT changed the backbone alignment for {gene}")

    occupied = compress(range(1, len(old_columns) + 1), old_columns)
    with column_map.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("alignment_column", "backbone_column"))
        writer.writerows((i, next(occupied) if present else "") for i, present in enumerate(new_columns, 1))


def validate_and_write(core, alignment, candidates, selected, output, column_map):
    """Validate sequence identity and reference columns without loading the full MSA."""
    core_rows = {gene: str(record.seq).upper() for gene, record in SeqIO.to_dict(SeqIO.parse(core, "fasta")).items()}
    if core_rows.keys() != selected.keys():
        raise ValueError("MAFFT changed the backbone IDs")
    final_core = {}
    alignment_width = None
    with alignment.open("rb") as aligned:
        offsets = fasta_offsets(aligned)
        if offsets.keys() != candidates.keys():
            raise ValueError("MAFFT changed the family sequence IDs")
        with (nullcontext(sys.stdout) if output == Path("-") else output.open("w")) as handle:
            for gene, original in candidates.items():
                start, length = offsets[gene]
                aligned.seek(start)
                raw = aligned.read(length).partition(b"\n")[2]
                sequence = raw.translate(None, b" \t\r\n").upper().decode("ascii")
                if alignment_width is None:
                    alignment_width = len(sequence)
                if len(sequence) != alignment_width or sequence.translate(GAPS) != original:
                    raise ValueError(f"MAFFT changed residues or produced an unaligned row: {gene}")
                if gene in selected:
                    final_core[gene] = sequence
                handle.write(f">{gene}\n{sequence}\n")
    write_column_map(core_rows, final_core, column_map)


def run(args):
    if args.size < 2 or args.threads < 1 or args.threadit < 0:
        raise ValueError("Expected a backbone size of at least two and valid MAFFT thread counts")
    candidates, queries = read_fasta(args.input), read_fasta(args.query)
    domains = read_domains(args.domains)
    with TemporaryDirectory(prefix="mafft-backbone-", dir=".") as directory:
        scratch = Path(directory)
        shared, query_domains = annotate_queries(queries, candidates, domains, args.rps_db, args.threads, scratch)
        selected, evidence, profiles = select_backbone(candidates, query_domains, shared, domains,
                                                       family_statuses(args.family_screen), get_species_parser(args),
                                                       args.size)
        fasta, core, result = (scratch / name for name in ("backbone.fa", "backbone.aln.fa", "family.aln.fa"))
        write_fasta(fasta, candidates, selected)
        if len(selected) == 1:
            core = fasta
        else:
            options = ["--localpair", "--maxiterate", "1000"] if len(selected) < len(candidates) else ["--auto"]
            align(fasta, core, args.threads, args.threadit, options)
        if len(selected) < len(candidates):
            fragments = scratch / "remaining.fa"
            write_fasta(fragments, candidates, (gene for gene in candidates if gene not in selected))
            align(core, result, args.threads, args.threadit, ["--auto", "--addfragments", str(fragments)])
        else:
            result = core
        validate_and_write(core, result, candidates, selected, args.output, args.column_map)

    ranks = {gene: i for i, gene in enumerate(selected, 1)}
    with args.report.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("gene_id", "backbone_rank", "backbone_reason", "species", "family_screen", "best_domain_query",
                         "ordered_domain_coverage", "query_domain_profiles"))
        for gene, (status, score, matched, species) in evidence.items():
            writer.writerow((gene, ranks.get(gene, ""), selected.get(gene, "added"), species, status, matched,
                             f"{score:.6f}", len(profiles)))
    print(f"MAFFT backbone: {len(selected)} representatives, {len(candidates) - len(selected)} added, "
          f"{len(profiles)} queries with annotated domains", file=sys.stderr if args.output == Path("-") else sys.stdout,
          flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("input", "query", "domains", "rps-db", "report", "column-map"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="Aligned FASTA path, or '-' for standard output")
    parser.add_argument("--family-screen", type=Path)
    parser.add_argument("--size", type=int, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--threadit", type=int, default=0)
    parser.add_argument("--species-parser", default="extended-taxonomic")
    parser.add_argument("--species-regex")
    parser.add_argument("--species-map-tsv")
    args = parser.parse_args()
    try:
        run(args)
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f"MAFFT backbone failed: {exc}\n")


if __name__ == "__main__":
    main()
