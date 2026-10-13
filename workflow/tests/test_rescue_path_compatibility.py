"""Alternative exons need homologous genomic projection, not common 5-prime ends."""
import copy
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
import rescue_coding_paths as paths
import rescue_path_compatibility as comparison

CODONS = ['ATG', 'GCT', 'CGT', 'AAT', 'GAT', 'TGT', 'CAA', 'GAA', 'GGT', 'CAT',
          'ATT', 'CTG', 'AAA', 'TTT', 'CCT', 'TCT', 'ACT', 'TGG', 'TAT', 'GTT']
DNA = ''.join(CODONS * 3) + 'TAA'


def model(blocks, sequence=DNA, strand='+'):
    return {'seqid': 'chr1', 'strand': strand, 'cds': blocks, 'sequence': sequence,
            'genetic_code': 1, 'coverage': 1., 'identity': .9, 'problems': [],
            'evidence': {'donor': 'Relative', 'query': 'gene'}}


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('shared_first', [True, False])
def test_low_overlap_paths_can_share_an_internal_or_terminal_exon(strand, shared_first):
    n = len(DNA)
    if shared_first:
        a = [[100, 118, 0], [300, 300 + n - 18, 0]]
        b = [[100, 118, 0], [600, 600 + n - 18, 0]]
    else:
        a = [[100, 118, 0], [500, 521, 0], [700, 700 + n - 39, 0]]
        b = [[300, 318, 0], [500, 521, 0], [900, 900 + n - 39, 0]]
    if strand == '-':
        a = [[2000 - e, 2000 - s, p] for s, e, p in a]
        b = [[2000 - e, 2000 - s, p] for s, e, p in b]
    first, second = model(a, strand=strand), model(b, strand=strand)
    assert paths.coding_overlap(first, second)[0] / n < .8
    evidence = comparison.protein_path_support(first, second)
    assert evidence['supported'] and evidence['shared_complete_codons'] > 0
    assert paths.compatible_paths(first, second)
    assert not paths.compatible_paths(first, second, owners_a=['gene_a'], owners_b=['gene_b'])


def test_tandem_like_shared_tail_head_is_not_a_homologous_codon_projection():
    n = len(DNA)
    a = model([[100, 100 + n, 0]])
    b = model([[100 + n - 18, 100 + 2 * n - 18, 0]])
    assert comparison.protein_path_support(a, b)['reason'] == 'contradictory_genomic_projection'
    assert not paths.compatible_paths(a, b)


def test_frame_strand_and_unknown_translation_cannot_be_merged():
    n = len(DNA)
    a = model([[100, 118, 0], [300, 300 + n - 18, 0]])
    b = model([[100, 118, 1], [600, 600 + n - 18, 1]])
    assert not paths.compatible_paths(a, b)
    b = copy.deepcopy(a)
    b['strand'] = '-'
    assert not paths.compatible_paths(a, b)
    b = model([[100, 118, 0], [600, 600 + n - 18, 0]])
    b.pop('genetic_code')
    assert not paths.compatible_paths(a, b)


def test_alignment_budget_abstains_without_changing_models(monkeypatch):
    n = len(DNA)
    a = model([[100, 118, 0], [300, 300 + n - 18, 0]])
    b = model([[100, 118, 0], [600, 600 + n - 18, 0]])
    before = copy.deepcopy([a, b])
    monkeypatch.setattr(comparison, 'MAX_ALIGNMENT_CELLS', 1)
    assert comparison.protein_path_support(a, b)['reason'] == 'protein_alignment_cell_budget_exceeded'
    assert [a, b] == before


def test_original_owner_guard_is_applied_before_unanchored_reassessment():
    from rescue_additional_candidates import reassess_unanchored_models
    n = len(DNA)
    a = model([[100, 118, 0], [300, 300 + n - 18, 0]])
    b = model([[100, 118, 0], [600, 600 + n - 18, 0]])
    for m in (a, b):
        m['problems'] = ['unanchored_genome_search']
    owners = paths.OwnershipIndex([
        {'seqid': 'chr1', 'strand': '+', 'start': 300, 'end': 300 + n - 18,
         'gene_id': 'original_a', 'cds': [[300, 300 + n - 18]]},
        {'seqid': 'chr1', 'strand': '+', 'start': 600, 'end': 600 + n - 18,
         'gene_id': 'original_b', 'cds': [[600, 600 + n - 18]]},
    ])
    reassess_unanchored_models([a, b], ownership=owners)
    assert all(m['placement_evidence']['reason'] == 'incompatible_unanchored_paths' for m in (a, b))
    assert all(m['problems'] == ['unanchored_genome_search'] for m in (a, b))
