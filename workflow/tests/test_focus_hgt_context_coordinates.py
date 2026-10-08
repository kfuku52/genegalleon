"""Indexed compression must retain the scalar genomic-coordinate geometry."""

import math
import random
import re
import struct
import sys
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / 'support'
sys.path.insert(0, str(SUPPORT))

from focus_hgt_context import GapCompressedCoordinates, GenomeCoordinates, model_span, render_context  # noqa: E402
from focus_hgt_gene_trees import write  # noqa: E402


def scalar_transform(display, position):
    reduction = 0
    for gap in display.gaps:
        covered = min(max(position - gap['genomic_start_bp'], 0), gap['original_gap_bp'])
        reduction += covered * gap['omitted_bp'] / gap['original_gap_bp']
    return position - reduction


def scalar_gap_type(gap, spans):
    left, right = gap['genomic_start_bp'], gap['genomic_end_exclusive_bp']
    within_gene = any(start <= left and right <= end for start, end in spans)
    mixed = any(start < right and end > left for start, end in spans)
    return 'intronic' if within_gene else 'mixed_noncoding' if mixed else 'intergenic'


def model(gene, start, end, **extra):
    return dict(gene_id=gene, chromosome='s1', start=str(start), end=str(end),
                feature_type='CDS', feature_blocks=f'{start}-{end}', utr_blocks='', strand='+', **extra)


@pytest.mark.parametrize('origin', [0, 10**17])
def test_indexed_transform_matches_scalar_bits_at_boundaries_midpoints_and_outer_positions(origin):
    rng = random.Random(82)
    models = [model(f'g{i}', origin + 1 + i * 30007, origin + 101 + i * 30007) for i in range(80)]
    focal = models[40]
    display = GapCompressedCoordinates(focal, models)
    positions = [origin - 1, origin, int(focal['start']), display.center, -0.0,
                 float('inf'), -float('inf'), float('nan')]
    for gap in display.gaps:
        for boundary in (gap['genomic_start_bp'], gap['genomic_end_exclusive_bp']):
            positions.extend([boundary - 1, boundary, boundary + 1, boundary + 0.5])
    positions.extend(origin + rng.randrange(2600000) + 0.5 for _ in range(200))
    for position in positions:
        expected, actual = scalar_transform(display, position), display.transform(position)
        if isinstance(expected, float) and math.isnan(expected):
            assert math.isnan(actual)
        else:
            assert type(actual) is type(expected)
            assert struct.pack('d', actual) == struct.pack('d', expected)
    assert display.anchor == scalar_transform(display, display.center)
    assert display.point(display.center) == 0
    assert display.point(int(models[0]['start'])) < 0


def test_scalar_numeric_behavior_and_empty_or_oversized_gap_cases_are_retained():
    left, right = model('left', 1, 100), model('right', 10000, 10100)
    display = GapCompressedCoordinates(left, [left, right])
    for position in [Decimal('5000.5'), Fraction(10001, 2), True]:
        expected, actual = scalar_transform(display, position), display.transform(position)
        assert type(actual) is type(expected) and actual == expected
    empty = GapCompressedCoordinates(None, [])
    assert type(empty.transform(1)) is int and empty.transform(1) == 1
    oversized = GapCompressedCoordinates(left, [left, model('distant', 10**400, 10**400 + 100)])
    assert oversized.transform(50.5) == scalar_transform(oversized, 50.5)
    with pytest.raises(OverflowError):
        oversized.transform(10**400)


@pytest.mark.parametrize('origin', [0, 10**17])
def test_indexed_gap_types_match_scalar_for_nested_partial_unknown_utr_and_trans_spliced_models(origin):
    models = [model('left', 1, 100), model('intronic', 7000, 20000),
              model('nested', 7000, 9000), model('partial', 24000, 35000),
              model('unknown', 42000, 49000), model('utr', 60000, 70000),
              model('exonic', 80000, 93000), model('trans', 100000, 111000)]
    models[1]['feature_blocks'] = '7000-7100;19900-20000'
    models[2]['feature_blocks'] = '7000-7100;8900-9000'
    models[3]['feature_blocks'] = '34900-35000'
    models[4]['feature_type'] = 'gene'
    models[5]['feature_blocks'] = '60000-60100;69900-70000'
    models[5]['utr_blocks'] = '59900-59999;70001-70100'
    models[6].update(feature_type='exon', feature_blocks='80000-80100;92900-93000')
    models[7].update(splice_mode='trans', feature_blocks='100000-100100;110900-111000')
    if origin:
        for row in models:
            row['start'], row['end'] = str(int(row['start']) + origin), str(int(row['end']) + origin)
            for field in ['feature_blocks', 'utr_blocks']:
                row[field] = ';'.join(f'{int(start)+origin}-{int(end)+origin}'
                                      for part in row[field].split(';') if part
                                      for start, end in [part.split('-')])
    focal = models[1]
    spans = [model_span(row) for row in models]
    display = GapCompressedCoordinates(focal, models)
    assert {gap['gap_type'] for gap in display.gaps} == {'intronic', 'mixed_noncoding', 'intergenic'}
    assert all(gap['gap_type'] == scalar_gap_type(gap, spans) for gap in display.gaps)
    random.Random(82).shuffle(models)
    reordered = GapCompressedCoordinates(focal, models)
    assert reordered.gaps == display.gaps
    assert reordered.audit() == display.audit()


def test_overlapping_loci_do_not_jointly_count_as_one_gap_spanning_gene():
    left, right = model('left', 1, 10000), model('right', 9000, 20000)
    left['feature_blocks'], right['feature_blocks'] = '1-100', '19900-20000'
    display = GapCompressedCoordinates(left, [left, right])
    assert len(display.gaps) == 1
    assert display.gaps[0]['gap_type'] == 'mixed_noncoding'
    assert display.gaps[0]['gap_type'] == scalar_gap_type(display.gaps[0], [model_span(left), model_span(right)])


def test_native_context_pdf_and_complete_audit_match_scalar_compression(tmp_path, monkeypatch):
    models = []
    for i, (start, end) in enumerate([(1, 100), (7000, 7100), (20000, 35000), (43000, 43100), (60000, 60100)]):
        row = model('placeholder' if i == 2 else f'neighbor{i}', start, end)
        if i == 0:
            row['feature_blocks'] = '1-50'
        if i == 2:
            row['feature_blocks'] = '20000-20100;34900-35000'
            row['utr_blocks'] = '19950-19999'
        models.append(row)
    root = tmp_path / 'gff'
    root.mkdir()
    for species in ['A', 'D']:
        rows = [dict(row, gene_id=species + ('_gene' if i == 2 else '_n' + str(i))) for i, row in enumerate(models)]
        write(root / (species + '.gff_info.tsv'), list(rows[0]), rows)
    stat = [dict(branch_id='2', node_name='root', child1='0', child2='1', support_generax_ufboot='45'),
            dict(branch_id='0', node_name='D_gene', child1='-999', child2='-999', support_generax_ufboot=''),
            dict(branch_id='1', node_name='A_gene', child1='-999', child2='-999', support_generax_ufboot='')]
    events = [dict(event_id='OG1:2:1', orthogroup='OG1', gene_tree_branch_id='2', gene_tree_node='root')]
    links = [dict(event_id='OG1:2:1', orthogroup='OG1', gene_id=species + '_gene', gene_species=species,
                  side=side, eligible_for_context='True', host_scaffold_status='measured', host_scaffold_id='s1',
                  host_scaffold_background_class_total_count='20', host_scaffold_background_class_compatible_count='9',
                  host_scaffold_background_class_incompatible_count='1', host_scaffold_background_class_unresolved_count='10',
                  host_scaffold_background_class_classified_fraction='0.5',
                  host_scaffold_background_class_compatible_fraction='0.9')
             for species, side in [('D', 'donor'), ('A', 'recipient')]]
    indexed = tmp_path / 'indexed.pdf'
    actual = render_context(indexed, stat, events, links, GenomeCoordinates(root), gene_tree_panel=False, max_genes_per_side=3)
    monkeypatch.setattr(GapCompressedCoordinates, 'transform', scalar_transform)
    indexed_init = GapCompressedCoordinates.__init__
    def scalar_classification_init(self, focal, neighbors):
        indexed_init(self, focal, neighbors)
        spans = [model_span(row) for row in neighbors]
        for gap in self.gaps:
            gap['gap_type'] = scalar_gap_type(gap, spans)
    monkeypatch.setattr(GapCompressedCoordinates, '__init__', scalar_classification_init)
    baseline = tmp_path / 'scalar.pdf'
    expected = render_context(baseline, stat, events, links, GenomeCoordinates(root), gene_tree_panel=False, max_genes_per_side=3)
    assert actual == expected
    def normalize(path):
        return re.sub(rb'/(CreationDate|ModDate) \([^)]*\)', b'', path.read_bytes())
    assert normalize(indexed) == normalize(baseline)
