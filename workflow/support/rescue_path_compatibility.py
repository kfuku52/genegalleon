"""Whole-protein support for alternative CDS paths with low genomic overlap.

Callers validate genomic ORFs, donor alignments and original ownership first.
Protein similarity supports predicted alternatives, not expressed isoforms or
orthology. Shared genomic codons must map to homologous protein positions.
"""
from functools import lru_cache

from Bio.Align import PairwiseAligner
from Bio.Data import CodonTable
from Bio.Seq import Seq

MINIMUM_PAIR_COVERAGE = .90
MINIMUM_PAIR_IDENTITY = .50
MAX_ALIGNMENT_CELLS = 25_000_000
MAX_OPTIMAL_ALIGNMENTS = 8


def _protein(model):
    code = model.get('genetic_code')
    sequence = model.get('sequence', '')
    if type(code) is not int or code not in CodonTable.unambiguous_dna_by_id or len(sequence) % 3:
        return None
    cumulative = 0
    for start, end, phase in sorted(model['cds'], reverse=model['strand'] == '-'):
        if end <= start or phase != (-cumulative) % 3:
            return None
        cumulative += end - start
    if cumulative != len(sequence):
        return None
    table = CodonTable.unambiguous_dna_by_id[code]
    if (not sequence or set(sequence) - set('ACGT') or sequence[:3] not in table.start_codons
            or sequence[-3:] not in table.stop_codons
            or any(sequence[i:i + 3] in table.forward_table for i in range(0, len(sequence), 3)
                   if sequence[i:i + 3] in table.stop_codons)):
        return None
    translated = str(Seq(sequence).translate(table=code)).removesuffix('*')
    return translated if translated and '*' not in translated else None


def _anchors(a, b, lengths):
    """Collect complete shared codons in transcript coordinates, on either strand."""
    shared = {}
    offset_a = 0
    for sa, ea, _ in sorted(a['cds'], reverse=a['strand'] == '-'):
        offset_b = 0
        for sb, eb, _ in sorted(b['cds'], reverse=b['strand'] == '-'):
            for position in range(max(sa, sb), min(ea, eb)):
                na = offset_a + (position - sa if a['strand'] == '+' else ea - 1 - position)
                nb = offset_b + (position - sb if b['strand'] == '+' else eb - 1 - position)
                if na % 3 != nb % 3:
                    return ()
                if na // 3 < lengths[0] and nb // 3 < lengths[1]:
                    key = na // 3, nb // 3
                    shared[key] = shared.get(key, 0) | (1 << (na % 3))
            offset_b += eb - sb
        offset_a += ea - sa
    return tuple(sorted(pair for pair, phases in shared.items() if phases == 7))


@lru_cache(maxsize=256)
def _alignment_metrics(protein_a, protein_b, anchors):
    aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1,
                              open_gap_score=-5, extend_gap_score=-.5)
    alignments = aligner.align(protein_a, protein_b)
    metrics = []
    for index in range(MAX_OPTIMAL_ALIGNMENTS):
        try:
            alignment = alignments[index]
        except IndexError:
            break
        pairs = {}
        for (a, b), (c, d) in zip(*alignment.aligned, strict=True):
            pairs.update((int(a + k), int(c + k)) for k in range(int(b - a)))
        if not all(pairs.get(a) == b for a, b in anchors):
            return 'contradictory_genomic_projection', ()
        aligned = len(pairs)
        identity = sum(protein_a[a] == protein_b[b] for a, b in pairs.items()) / max(1, aligned)
        metrics.append((aligned / len(protein_a), aligned / len(protein_b), identity))
    try:
        alignments[MAX_OPTIMAL_ALIGNMENTS]
    except IndexError:
        return 'evaluated', tuple(metrics)
    return 'optimal_alignment_budget_exceeded', ()


def protein_path_support(a, b, minimum_coverage=MINIMUM_PAIR_COVERAGE):
    """Return an auditable comparison; incomplete searches cannot justify a merge."""
    if a['seqid'] != b['seqid'] or a['strand'] != b['strand']:
        return {'supported': False, 'reason': 'different_sequence_or_strand'}
    if a.get('genetic_code') != b.get('genetic_code'):
        return {'supported': False, 'reason': 'different_genetic_codes'}
    protein_a, protein_b = _protein(a), _protein(b)
    if not protein_a or not protein_b:
        return {'supported': False, 'reason': 'unknown_or_unsupported_translation'}
    if len(protein_a) * len(protein_b) > MAX_ALIGNMENT_CELLS:
        return {'supported': False, 'reason': 'protein_alignment_cell_budget_exceeded'}
    anchors = _anchors(a, b, (len(protein_a), len(protein_b)))
    if not anchors:
        return {'supported': False, 'reason': 'no_shared_complete_codon_anchor'}
    reason, rows = _alignment_metrics(protein_a, protein_b, anchors)
    supported = reason == 'evaluated' and any(min(ca, cb) >= minimum_coverage
                                             and identity >= MINIMUM_PAIR_IDENTITY for ca, cb, identity in rows)
    return {'supported': supported, 'reason': reason if reason != 'evaluated' or supported else 'insufficient_pair_homology',
            'minimum_coverage': minimum_coverage, 'minimum_identity': MINIMUM_PAIR_IDENTITY,
            'shared_complete_codons': len(anchors),
            'alignments': [{'coverage_a': ca, 'coverage_b': cb, 'identity': identity} for ca, cb, identity in rows]}
