"""Strand/coding ownership, copy conflicts and one-locus alternative export."""
import copy
import json
import sys
import tracemalloc
from importlib import import_module
from pathlib import Path

import pytest
from Bio.Seq import Seq

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


@pytest.mark.parametrize('reverse', [False, True])
def test_independent_locus_support_admits_compatible_paths_without_inflating_path_support(reverse):
    # The shorter path has the higher donor identity, but not decisively so.
    # Canonical fallback must therefore be independent of the donor rank.
    long = model('long', [[0, 90, 0]], ('Relative',), identity=.74)
    short = model('short', [[6, 90, 0]], ('Phylogenetic',), identity=.76)
    before = [copy.deepcopy(row) for row in (long, short)]
    paths.resolve_coding_paths([short, long] if reverse else [long, short], nearest=['Relative', 'Phylogenetic'])
    assert long['status'] == 'accepted' and short['status'] == 'accepted_alternative_path'
    assert short['parent_model_id'] == 'long'
    assert long['alternative_coding_paths'][0]['model_id'] == 'short'
    assert long['locus_support'] == {'counting_unit': 'independent_donor_species',
        'independent_donor_species': ['Phylogenetic', 'Relative'], 'minimum_species_for_fallback': 2,
        'orthology': 'unassigned', 'expected_copy': 'unassigned'}
    selection = long['path_selection']
    assert selection['reason'] == 'compatible_locus_canonical_fallback'
    assert selection['representative_status'] == 'ambiguous'
    assert selection['selection_policy'] == 'longest_cds_then_coding_shape'
    assert selection['ambiguity_reason'] == 'ambiguous_coding_path_representative'
    assert selection['rank'] == [1, 1, .74, 1.]
    assert selection['ranked_paths'][0]['model_id'] == 'short'
    for row, original in zip((long, short), before, strict=True):
        for field in ('sequence', 'cds', 'support', 'coverage', 'identity', 'evidence', 'problems'):
            assert row[field] == original[field]
        assert paths.path_rank(row, ['Relative', 'Phylogenetic'])[:2] == (1, 1)


@pytest.mark.parametrize('reverse', [False, True])
def test_canonical_equal_length_tie_is_coordinate_deterministic(reverse):
    a = model('a', [[0, 60, 0]], ('Relative',))
    b = model('b', [[3, 63, 0]], ('Phylogenetic',))
    paths.resolve_coding_paths([b, a] if reverse else [a, b])
    assert a['status'] == 'accepted' and b['status'] == 'accepted_alternative_path'
    assert a['path_selection']['representative_status'] == 'ambiguous'


def test_duplicate_queries_and_self_do_not_supply_independent_fallback_species():
    a = model('a', [[0, 60, 0]], ('Relative',))
    b = model('b', [[6, 60, 0]], ('Relative',))
    b['support'].extend([{'donor': 'Relative', 'query': 'isoform2',
                         'alignment': {'coverage': 1., 'identity': 1., 'problems': []}},
                        {'donor': 'Weak', 'query': 'isoform',
                         'alignment': {'coverage': .3, 'identity': 1., 'problems': ['low_coverage']}},
                        {'donor': 'Invalid', 'query': 'isoform',
                         'alignment': {'coverage': True, 'identity': 1., 'problems': []}}])
    paths.resolve_coding_paths([a, b])
    assert all(row['status'] == 'unresolved' for row in (a, b))
    a = model('a', [[0, 60, 0]], ('Self',))
    b = model('b', [[6, 60, 0]], ('Relative',))
    for row in (a, b):
        row['evidence']['target'] = 'Self'
    paths.resolve_coding_paths([a, b])
    assert all(row['status'] == 'unresolved' for row in (a, b))


@pytest.mark.parametrize('problem', ['frameshift', 'internal_stop', 'noncanonical_splice', 'invalid_phase',
                                     'low_coverage', 'overlap_existing_annotation', 'assembly_gap_or_ambiguity'])
def test_canonical_fallback_never_promotes_failed_path_or_masks_its_gate(problem):
    a = model('a', [[0, 60, 0]], ('Relative',))
    b = model('b', [[6, 60, 0]], ('Phylogenetic',))
    b['problems'].append(problem)
    paths.resolve_coding_paths([a, b])
    assert all(row['status'] == 'unresolved' for row in (a, b))
    assert problem in b['problems']
    assert all('alternative_coding_paths' not in row for row in (a, b))


def test_hard_splice_frame_conflict_is_not_rescued_by_locus_union_support():
    a = model('a', [[0, 60, 0]], ('Relative',))
    b = model('b', [[1, 61, 0]], ('Phylogenetic',))
    paths.resolve_coding_paths([a, b])
    assert all(row['status'] == 'unresolved' and 'competing_new_models' in row['problems'] for row in (a, b))


def validated_pair(strand='+'):
    coding = 'ATGAAAATGAAACCCGGGCCCTAA'
    dna = 'C' * 100 + coding + 'C' * 20
    blocks = ([[100, 124, 0]], [[106, 124, 0]])
    if strand == '-':
        blocks = tuple([[len(dna) - end, len(dna) - start, phase] for start, end, phase in path] for path in blocks)
        dna = str(Seq(dna).reverse_complement())
    class Genome:
        def fetch(self, seqid, start, end):
            assert seqid == 'chr1' and 0 <= start <= end <= len(dna)
            return dna[start:end]
        def get_reference_length(self, seqid):
            assert seqid == 'chr1'
            return len(dna)
    rows = [model('long', blocks[0], ('Relative',), strand=strand, identity=.74),
            model('short', blocks[1], ('Phylogenetic',), strand=strand, identity=.76)]
    checked = [rescue.validate_model({**row, 'frameshift': False}, Genome(), 1,
               {'minimum_coverage': .95, 'minimum_identity': .6}) for row in rows]
    assert all(row['problems'] == [] for row in checked)
    assert [row['sequence'] for row in checked] == [coding, coding[6:]]
    return checked


@pytest.mark.parametrize('strand', ['+', '-'])
def test_actual_genomic_orfs_admit_one_locus_and_preserve_original_ownership(strand):
    rows = validated_pair(strand)
    rescue.consolidate(rows, [], 'Species', nearest=['Relative', 'Phylogenetic'])
    primary, = [row for row in rows if row['status'] == 'accepted']
    assert len(primary['sequence']) == 24 and len(primary['alternative_coding_paths']) == 1
    assert primary['path_selection']['representative_status'] == 'ambiguous'
    assert primary['locus_support']['orthology'] == primary['locus_support']['expected_copy'] == 'unassigned'
    rows = validated_pair(strand)
    owner = {'seqid': 'chr1', 'strand': strand, 'start': min(b[0] for b in rows[0]['cds']),
             'end': max(b[1] for b in rows[0]['cds']), 'gene_id': 'healthy_original', 'cds': rows[0]['cds']}
    rescue.consolidate(rows, [owner], 'Species', nearest=['Relative', 'Phylogenetic'])
    assert all(row['status'] == 'unresolved' and row['revision_owner_ids'] == ['healthy_original'] for row in rows)
    assert all('overlap_existing_annotation' in row['problems'] and 'alternative_coding_paths' not in row for row in rows)


def test_unanchored_locus_corroboration_and_path_preference_remain_distinct():
    rows = validated_pair()
    for row in rows:
        row['evidence']['target'] = 'Species'
        row['evidence']['genome_only'] = True
        row['problems'] = ['unanchored_genome_search']
    placement = rescue.reassess_unanchored_models(rows, 2)
    assert placement['counts'] == {'supported_unanchored_annotation': 2}
    rescue.consolidate(rows, [], 'Species', nearest=['Relative', 'Phylogenetic'])
    primary, = [row for row in rows if row['status'] == 'accepted']
    assert primary['placement_evidence']['classification'] == 'unanchored'
    assert primary['path_selection']['representative_status'] == 'ambiguous'
    assert len(primary['support']) == 1 and len(primary['locus_support']['independent_donor_species']) == 2


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


@pytest.mark.parametrize('fallback', [False, True])
def test_finalize_exports_alternatives_under_one_gene_and_keeps_original_ids(tmp_path, monkeypatch, fallback):
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
    primary = model(name + '_ggrescue_primary', [[100, 118, 0]],
                    ('Relative',) if fallback else ('Relative', 'Phylogenetic'))
    primary['sequence'] = 'ATGAAACCCGGGCCCTAA'
    alternative = model(name + '_ggrescue_alternative', [[100, 106, 0], [112, 118, 0]],
                        ('Phylogenetic',) if fallback else ('Relative',))
    alternative['sequence'] = 'ATGAAACCCTAA'
    paths.resolve_coding_paths([primary, alternative])
    assert primary['path_selection']['representative_status'] == ('ambiguous' if fallback else 'supported_priority')
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
    gene_attributes, = [rescue.attributes(row[8]) for row in rows if row[2] == 'gene']
    assert gene_attributes['representative_status'] == ('ambiguous' if fallback else 'supported_priority')
    assert gene_attributes['locus_independent_donor_species_count'] == '2'
    assert gene_attributes['rescue_orthology'] == gene_attributes['rescue_expected_copy'] == 'unassigned'
    transcript_attributes = [rescue.attributes(row[8]) for row in rows if row[2] == 'mRNA']
    assert sorted(a['path_independent_donor_species_count'] for a in transcript_attributes) == (['1', '1'] if fallback else ['1', '2'])
    assert {rescue.attributes(row[8])['Parent'] for row in rows if row[2] == 'mRNA'} == {primary['model_id']}
    catalog = catalog_helper.build_catalog(name, output_cds, output_gff, genome)
    rescued = next(g for g in catalog['loci'] if g['gene_id'] == primary['model_id'])
    assert len(rescued['candidates']) == 2
    original_locus, = [g for g in catalog['loci'] if g['source_gene_id'] == 'original']
    assert len(original_locus['candidates']) == 1 and original_locus['candidates'][0]['cds'] == 'ATGAAACCCTAA'
    assert all(not c['quality']['sequence_mismatch'] for c in rescued['candidates'])
    assert all(c['rescue_path_selection']['representative_status'] == gene_attributes['representative_status']
               and c['rescue_path_selection']['locus_independent_donor_species_count'] == '2'
               for c in rescued['candidates'])
    alternate, = [c for c in rescued['candidates'] if c.get('rescue_alternative_coding_path')]
    assert alternate['origin'] == 'predicted' and not alternate['quality']['representative_eligible']
    assert catalog['summary']['rescue_alternative_coding_paths'] == 1
