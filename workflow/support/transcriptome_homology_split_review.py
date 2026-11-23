"""Preview which rnaSPAdes family splits survive representative selection."""

import csv
from collections import defaultdict
from itertools import combinations

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqRecord import SeqRecord
from cdskit.homology import HIT_COLUMNS, parse_hit, select_representative
from cdskit.homology_split import (
    FAMILY_SPLIT_DECISIONS,
    SPLIT_MIN_COVERAGE,
    connected_targets,
    propose_family_partition,
    review_selected_families,
)
from cdskit.tsvio import write_tsv
from transcriptome_orf_translation import make_orf_peptide_translator

PREVIEW_COLUMNS = (
    "original_group_id", "candidate_decision", "decision", "representatives", "transcript_groups",
)
SUPPORTED_PREVIEW_DECISIONS = frozenset({"distinct_representatives", "supported_distinct_loci"})
# Full-length similarity is a duplication flag without independent locus evidence.
MIN_PROTEIN_IDENTITY = 0.95


def similarity_veto(peptides):
    aligner = PairwiseAligner(match_score=2, mismatch_score=-1, open_gap_score=-4, extend_gap_score=-0.5)
    for left, right in combinations((peptide for peptide in peptides if peptide), 2):
        required = MIN_PROTEIN_IDENTITY * max(len(left), len(right))
        if min(len(left), len(right)) >= required:
            if aligner.align(left, right)[0].counts().identities >= required:
                return True
    return False


def fasta_records(path, source_name, assignment):
    with path.open() as stream:
        for title, sequence in SimpleFastaParser(stream):
            identifier = title.split()[0]
            source = source_name(identifier)
            if source in assignment:
                yield assignment[source], SeqRecord(Seq(sequence), id=identifier, description=title)


def selected_hits(path, query_ids):
    hits = defaultdict(dict)
    with path.open(newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        if tuple(next(reader, ())) != HIT_COLUMNS:
            raise ValueError("Invalid protein hit header for split preview")
        for row in reader:
            if row[0] in query_ids:
                query, hit = parse_hit(row)
                previous = hits[query].get(hit.target)
                if previous is None or hit.rank > previous.rank:
                    hits[query][hit.target] = hit
    return hits


def preview_selections(args, groups, audit, target_families, loci, hits, source_name, alternates):
    with args.report.open(newline="") as stream:
        originals = {
            row["group_id"]: source_name(row["selected_id"])
            for row in csv.DictReader(stream, delimiter="\t")
        }
    proposals, assignment = [], {}
    for row in audit:
        if row["decision"] not in FAMILY_SPLIT_DECISIONS:
            continue
        group = row["original_group_id"]
        partition, families = propose_family_partition(
            groups[group], row["transcript_target_families"], originals.get(group, "")
        )
        child_by_family = {family: f"{group}.split{i}" for i, family in enumerate(families, 1)}
        children = {member: child_by_family[family] for member, family in partition.items()}
        assignment.update(children)
        # Per-transcript ambiguous hits can connect previously separate target labels.
        components = connected_targets(row["transcript_target_families"].values())
        proposals.append((row, children, child_by_family, components))
    if not assignment:
        return
    baseline = sorted(
        fasta_records(args.split_extra_baseline, source_name, assignment), key=lambda item: item[1].id
    )
    alternatives = list(fasta_records(alternates, source_name, assignment))
    evidence = selected_hits(hits, {record.id for _, record in baseline + alternatives})
    records = defaultdict(list)
    for child, record in baseline:
        records[child].append(record)
    for child, record in alternatives:
        if record.id in evidence:
            records[child].append(record)
    translate_orf = make_orf_peptide_translator(args.genetic_code)
    for row, children, child_by_family, components in proposals:
        expected = {child: family for family, child in child_by_family.items()}
        observed, selected_loci, representatives, sequences, peptides = {}, {}, {}, {}, []
        for child, family in expected.items():
            _, record, _, _, hit = select_representative(records[child], evidence)
            sequences[child] = record
            source = source_name(record.id)
            if source in loci:
                selected_loci[child] = loci[source]
            peptide = translate_orf(str(record.seq))
            peptides.append(peptide)
            hit_family = ""
            if hit:
                canonical = target_families.get(hit.target, hit.target)
                hit_family = components.get(canonical, canonical)
                if min(hit.query_coverage, hit.target_coverage) >= SPLIT_MIN_COVERAGE:
                    observed[child] = hit_family
            representatives[child] = {
                "selected_id": record.id, "expected_family": family, "selected_family": hit_family,
                "target_id": hit.target if hit else "",
                "query_coverage": hit.query_coverage if hit else None,
                "target_coverage": hit.target_coverage if hit else None,
                "locus": selected_loci.get(child, ""),
            }
        yield dict(
            original_group_id=row["original_group_id"], candidate_decision=row["decision"],
            decision=review_selected_families(
                expected, observed, selected_loci, similar_proteins=similarity_veto(peptides)
            ),
            representatives=representatives, transcript_groups=children,
        ), sequences


def write_split_preview(args, review, hits, source_name, alternates):
    rows = (row for row, _ in preview_selections(args, *review, hits, source_name, alternates))
    write_tsv(args.split_preview, rows, PREVIEW_COLUMNS)
