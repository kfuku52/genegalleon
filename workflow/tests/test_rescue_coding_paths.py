"""Strand/coding ownership, copy conflicts and one-locus alternative export."""
import copy
import json
import sys
import tracemalloc
from importlib import import_module
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / 'support'
sys.path.insert(0, str(SUPPORT))
paths = import_module('rescue_coding_paths')
rescue = import_module('rescue_gene_models')
catalog_helper = import_module('gene_model_catalog')


def model(identifier, blocks, donors=('Relative',), *, strand='+', identity=1.0):
    evidence = {'donor': donors[0], 'query': 'gene'} if donors else {}
    return {'model_id': identifier, 'seqid': 'chr1', 'strand': strand, 'cds': blocks,
            'start': min(b[0] for b in blocks), 'end': max(b[1] for b in blocks),
            'sequence': 'A' * sum(b[1] - b[0] for b in blocks), 'coverage': 1.0, 'identity': identity,
            'status': 'accepted', 'problems': [], 'evidence': evidence,
            'support': [{'donor': donor, 'query': 'gene',
                         'alignment': {'coverage': 1.0, 'identity': identity, 'problems': []}} for donor in donors]}


@pytest.mark.parametrize('strand,collision', [('+', True), ('-', False), ('.', True), ('?', True)])
def test_unknown_strand_conservatively_owns_span_even_inside_its_coding_intron(strand, collision):
    original = {'seqid': 'chr1', 'start': 0, 'end': 100, 'gene_id': 'original', 'strand': strand,
                'cds': [[0, 12], [80, 100]]}
    prediction = model('new', [[40, 60, 0]])
    overlapping = paths.OwnershipIndex([original]).overlapping(prediction)
    assert bool(overlapping) is (collision and strand not in {'+', '-'})
    prediction = model('new', [[5, 23, 0]])
    assert bool(paths.OwnershipIndex([original]).overlapping(prediction)) is collision


def test_ownership_prefix_maximum_finds_nested_long_features_and_all_owners():
    features = [{'seqid': 'chr1', 'start': 0, 'end': 1000, 'gene_id': 'long', 'strand': '+'},
                {'seqid': 'chr1', 'start': 100, 'end': 120, 'gene_id': 'short', 'strand': '+'},
                {'seqid': 'chr1', 'start': 115, 'end': 200, 'gene_id': 'middle', 'strand': '+'}]
    assert {r['gene_id'] for r in paths.OwnershipIndex(features).overlapping(model('new', [[116, 122, 0]]))} == {'long', 'short', 'middle'}
    assert {r['gene_id'] for r in paths.OwnershipIndex(features).overlapping(model('new', [[900, 912, 0]]))} == {'long'}


def test_original_ownership_preserves_gene_ids_encoded_commas_and_noncoding_strand():
    text = ['chr1\tp\tgene\t1\t100\t.\t+\t.\tID=gene%2Ca\n',
            'chr1\tp\tmRNA\t1\t100\t.\t+\t.\tID=tx;Parent=gene%2Ca\n',
            'chr1\tp\tCDS\t1\t12\t.\t+\t0\tParent=tx\n',
            'chr1\tp\tCDS\t81\t100\t.\t+\t0\tParent=tx\n',
            'chr1\tp\tncRNA\t141\t170\t.\t?\t.\tID=nc\n']
    rows = paths.original_ownership(text, rescue.attributes)
    assert {r['gene_id'] for r in rows} == {'gene,a', 'nc'}
    coding = next(r for r in rows if r['gene_id'] == 'gene,a')
    assert coding['cds'] == [[0, 12], [80, 100]] and coding['strand'] == '+'
    assert next(r for r in rows if r['gene_id'] == 'nc')['strand'] == '.'


def test_exon_only_gtf_and_transcript_only_gff_bind_one_original_owner():
    gtf = ['chr1\tp\texon\t11\t20\t.\t-\t.\tgene_id "nc,g"; transcript_id "nc,t";\n',
           'chr1\tp\texon\t51\t70\t.\t-\t.\tgene_id "nc,g"; transcript_id "nc,t";\n']
    owner, = paths.original_ownership(gtf, rescue.attributes)
    assert (owner['gene_id'], owner['start'], owner['end'], owner['strand'], owner['cds']) == ('nc,g', 10, 70, '-', [])
    gff = ['chr1\tp\tmRNA\t1\t100\t.\t+\t.\tID=tx\n',
           'chr1\tp\tCDS\t1\t12\t.\t+\t0\tParent=tx\n',
           'chr1\tp\tCDS\t81\t100\t.\t+\t0\tParent=tx\n']
    owner, = paths.original_ownership(gff, rescue.attributes)
    assert owner['gene_id'] == 'tx' and owner['cds'] == [[0, 12], [80, 100]]


def test_connected_coding_components_keep_independent_intron_gene_out_of_conflict():
    a = model('a', [[0, 30, 0], [100, 130, 0]], ('Relative', 'Phylogenetic'))
    b = model('b', [[0, 24, 0], [100, 130, 0]], ('Relative',))
    independent = model('intron_gene', [[40, 70, 0]], ('Relative',))
    rows = [a, b, independent]
    assert paths.overlap_components(rows) == [[0, 1], [2]]
    paths.resolve_coding_paths(rows, nearest=['Relative'])
    assert a['status'] == 'accepted' and b['status'] == 'accepted_alternative_path'
    assert independent['status'] == 'accepted' and independent['problems'] == []
    assert len(a['alternative_coding_paths']) == 1


def test_bridge_prediction_cannot_fuse_two_independent_coding_loci():
    rows = [model('left', [[0, 60, 0]], ('Relative',)),
            model('bridge', [[30, 150, 0]], ('Relative', 'Phylogenetic')),
            model('right', [[120, 180, 0]], ('Relative',))]
    assert paths.overlap_components(rows) == [[0, 1, 2]]
    paths.resolve_coding_paths(rows)
    assert all(r['status'] == 'unresolved' and 'competing_new_models' in r['problems'] for r in rows)


@pytest.mark.parametrize('strand', ['+', '-'])
def test_shared_coding_frame_must_agree_on_both_strands(strand):
    a = model('a', [[0, 60, 0]], strand=strand)
    compatible = model('b', [[6, 60, 0]] if strand == '+' else [[0, 54, 0]], strand=strand)
    shifted = model('c', [[1, 61, 0]], strand=strand)
    assert paths.compatible_paths(a, compatible)
    assert not paths.compatible_paths(a, shifted)
    other_strand = {**copy.deepcopy(a), 'strand': '-' if strand == '+' else '+'}
    assert paths.coding_overlap(a, other_strand) == (0, 0)


def test_rank_counts_species_once_and_never_borrows_another_alignment_quality():
    a = model('a', [[0, 60, 0]], ('Relative',))
    a['support'].extend([{'donor': 'Relative', 'query': 'isoform2', 'alignment': {'coverage': 1, 'identity': .99, 'problems': []}},
                         {'donor': 'Unmeasured', 'query': 'gene'},
                         {'donor': 'Weak', 'query': 'gene', 'alignment': {'coverage': .3, 'identity': 1, 'problems': ['low_coverage']}}])
    assert paths.path_rank(a, ['Relative', 'Weak']) == (1, 1, 1.0, 1.0)
    # Reordering primary alignment cannot change path rank after deduplication.
    a['identity'] = .5
    assert paths.path_rank(a, ['Relative']) == (1, 1, 1.0, 1.0)


def test_ambiguous_same_locus_paths_remain_proposals_and_missing_provenance_cannot_break_tie():
    a, b = model('a', [[0, 60, 0]]), model('b', [[6, 60, 0]])
    paths.resolve_coding_paths([a, b], nearest=['Relative'])
    assert a['status'] == b['status'] == 'unresolved'
    assert a['problems'] == b['problems'] == ['ambiguous_coding_path_representative']
    a, b = model('a', [[0, 60, 0]], donors=(), identity=1), model('b', [[6, 60, 0]], donors=(), identity=.8)
    paths.resolve_coding_paths([a, b])
    assert a['status'] == b['status'] == 'unresolved'


def test_dense_ambiguous_component_does_not_allocate_a_pair_matrix():
    # Hundreds of competing unconfirmed paths must abstain without an
    # auxiliary all-pairs matrix, independent of the genomic span length.
    rows = [model(str(i), [[3 * i, 2400, 0]], donors=()) for i in range(600)]
    tracemalloc.start()
    try:
        paths.resolve_coding_paths(rows)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert all(row['status'] == 'unresolved' for row in rows)
    assert peak < 3_000_000


def test_finalize_exports_alternatives_under_one_gene_and_keeps_original_ids(tmp_path, monkeypatch):
    root = tmp_path / 'rescue'
    name = 'Species_target'
    worker = root / 'rescued' / name
    worker.mkdir(parents=True)
    prepared = root / 'prepared' / name
    prepared.mkdir(parents=True)
    (prepared / 'mapping.json').write_text(json.dumps({'feature': 'gene', 'attribute': 'ID'}))
    cds, gff, genome = [tmp_path / p for p in ('original.fa', 'original.gff3', 'genome.fa')]
    cds.write_text('>original\nATGAAACCCTAA\n')
    gff.write_text('##gff-version 3\nchr1\tp\tgene\t1\t12\t.\t+\t.\tID=original\n'
                   'chr1\tp\tmRNA\t1\t12\t.\t+\t.\tID=original_t;Parent=original\n'
                   'chr1\tp\tCDS\t1\t12\t.\t+\t0\tParent=original_t\n')
    genome.write_text('>chr1\nATGAAACCCTAA' + 'C' * 88 + 'ATGAAACCCGGGCCCTAA\n')
    primary = model(name + '_ggrescue_primary', [[100, 118, 0]], ('Relative', 'Phylogenetic'))
    primary['sequence'] = 'ATGAAACCCGGGCCCTAA'
    alternative = model(name + '_ggrescue_alternative', [[100, 106, 0], [112, 118, 0]])
    alternative['sequence'] = 'ATGAAACCCTAA'
    paths.resolve_coding_paths([primary, alternative])
    (worker / 'models.json').write_text(json.dumps([primary, alternative]))
    (worker / 'receipt.json').write_text('{}')
    plan = {'species': [name], 'common_references': [], 'request': {'sources': {name: {
        'species': name, 'fasta': str(cds), 'gff': str(gff), 'genome': str(genome), 'genetic_code': 1,
        'quality': {'complete_pct': 90}}}}}
    monkeypatch.setattr(rescue, 'plan_digest', lambda *_: 'frozen-plan')
    original_stage = rescue.stage
    def stage(output, destination, key, builder, guard, **kwargs):
        if Path(output) != root:
            return original_stage(output, destination, key, builder, guard, **kwargs)
        directory = output / destination
        directory.mkdir()
        builder(directory)
        return directory
    monkeypatch.setattr(rescue, 'stage', stage)
    def admission(source, directory, side, fraction, **kwargs):
        for suffix in ('json', 'tsv'):
            (directory / (side + '.anchor_admission.' + suffix)).write_text('{}')
        return [], {'anchor_admission': {'checked': True}}
    monkeypatch.setattr(rescue, 'prepare_rescue_genome', admission)
    result = rescue.finalize(root, plan)
    output_cds = result / 'species_cds' / (name + '.rescue.cds.fa')
    output_gff = result / 'species_gff' / (name + '.rescue.gff3')
    assert {i for i, _, _ in rescue.fasta_records(output_cds)} == {'original', primary['model_id']}
    text = output_gff.read_text()
    assert gff.read_text() in text
    rows = [line.split('\t') for line in text.splitlines() if '\tgenegalleon_rescue\t' in line]
    assert sum(row[2] == 'gene' for row in rows) == 1
    assert {rescue.attributes(row[8])['Parent'] for row in rows if row[2] == 'mRNA'} == {primary['model_id']}
    catalog = catalog_helper.build_catalog(name, output_cds, output_gff, genome)
    rescued = next(g for g in catalog['loci'] if g['gene_id'] == primary['model_id'])
    assert len(rescued['candidates']) == 2
    assert all(not c['quality']['sequence_mismatch'] for c in rescued['candidates'])
    alternate, = [c for c in rescued['candidates'] if c.get('rescue_alternative_coding_path')]
    assert alternate['origin'] == 'predicted' and not alternate['quality']['representative_eligible']
    assert catalog['summary']['rescue_alternative_coding_paths'] == 1
