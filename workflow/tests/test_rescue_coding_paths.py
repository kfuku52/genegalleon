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


def test_one_external_donor_admits_compatible_paths_but_missing_provenance_does_not():
    a, b = model('a', [[0, 60, 0]]), model('b', [[6, 60, 0]])
    paths.resolve_coding_paths([a, b], nearest=['Relative'])
    assert a['status'] == 'accepted' and b['status'] == 'accepted_alternative_path'
    assert a['path_selection']['representative_status'] == 'ambiguous'
    assert a['locus_support']['independent_donor_species'] == ['Relative']
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
        'independent_donor_species': ['Phylogenetic', 'Relative'], 'minimum_species_for_fallback': 1,
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


@pytest.mark.parametrize('engine', ['annevo', 'oriongeno'])
def test_native_support_preserves_valid_paths_with_ambiguous_representative(engine):
    a = model('a', [[0, 90, 0]], ('Relative', 'Phylogenetic'), identity=.74)
    b = model('b', [[6, 90, 0]], ('Relative', 'Phylogenetic'), identity=.76)
    b['support'][0][engine] = {'native_transcript_id': 'native.t'}
    paths.resolve_coding_paths([a, b])
    assert a['status'] == 'accepted' and b['status'] == 'accepted_alternative_path'
    assert a['path_selection']['representative_status'] == 'ambiguous'
    assert a['alternative_coding_paths'][0]['support'] == b['support']


@pytest.mark.parametrize('reverse', [False, True])
def test_native_duplicate_of_previous_long_path_cannot_remove_it_or_its_short_isoform(reverse):
    # Ancistrocladus Tripe_34502 regression: the new engine supports an
    # existing geometry, while donor identity does not choose a unique path.
    long = model('long', [[0, 30, 0], [90, 153, 0]], ('Relative',), identity=.82)
    short = model('short', [[93, 153, 0]], ('Phylogenetic',), identity=.826)
    native = copy.deepcopy(long['support'][0])
    native['annevo'] = {'native_transcript_id': 'native.long'}
    long['support'].append(native)
    original = copy.deepcopy([long, short])
    paths.resolve_coding_paths([short, long] if reverse else [long, short])
    assert long['status'] == 'accepted' and short['status'] == 'accepted_alternative_path'
    for row, before in zip([long, short], original, strict=True):
        for key in ['cds', 'sequence', 'support', 'evidence']:
            assert row[key] == before[key]


@pytest.mark.parametrize('bridge', [False, True])
@pytest.mark.parametrize('engine', ['annevo', 'oriongeno', 'boundary_candidate'])
def test_conflicting_new_addition_cannot_drop_or_fuse_transferred_cds(bridge, engine):
    left = model('left', [[0, 60, 0]])
    native = model('native', [[30, 150, 0]] if bridge else [[1, 61, 0]])
    native['support'][0][engine] = {'native_transcript_id': 'new.t'}
    native['evidence'][engine] = native['support'][0][engine]
    rows = [left, native]
    if bridge:
        rows.append(model('right', [[120, 180, 0]]))
    paths.resolve_coding_paths(rows)
    assert native['status'] == 'unresolved' and native['problems'] == ['competing_new_models']
    assert all(row['status'] == 'accepted' and row['problems'] == [] for row in rows if row is not native)
    assert all('alternative_coding_paths' not in row for row in rows)


def retained_path(row):
    return {paths.coding_shape(row): {'sequence_sha256': rescue.hashlib.sha256(row['sequence'].encode()).hexdigest(),
            'producer_receipt_sha256': '1' * 64, 'model_ids': ['previous'], 'source_directory': '/frozen/previous'}}


@pytest.mark.parametrize('prior_isoforms', [1, 2])
def test_new_self_only_alternatives_cannot_revoke_valid_prior_annotation(prior_isoforms):
    old = model('old', [[6, 96, 0]], ('Self',))
    previous = [old]
    if prior_isoforms == 2:
        previous.append(model('old_short', [[12, 96, 0]], ('Self',)))
    novel = model('new_longer', [[0, 96, 0]], ('Self',))
    independent = model('independent', [[120, 180, 0]], ('Relative',))
    for row in [*previous, novel]:
        row['evidence']['target'] = 'Self'
    retained = {key: value for row in previous for key, value in retained_path(row).items()}
    rescue.consolidate([*previous, novel, independent], [], 'Self', retained_paths=retained)
    assert old['status'] == independent['status'] == 'accepted'
    assert all(row['status'] == 'accepted_alternative_path' for row in previous[1:])
    assert novel['status'] == 'unresolved'
    assert novel['problems'] == ['ambiguous_coding_path_representative']
    assert novel['integration_conflict']['disposition'] == 'retain_previous_validated_annotation'
    assert old['locus_support']['independent_donor_species'] == []
    assert old['locus_support']['admission_basis'] == 'source_bound_currently_validated_prior_annotation'
    assert old['path_selection']['reason'] == 'retained_annotation_canonical_fallback'
    assert old['path_selection']['representative_status'] == 'ambiguous'
    assert len(old['alternative_coding_paths']) == prior_isoforms - 1


def test_prior_native_annotation_survives_conflicting_new_transfer_and_independent_locus():
    old = model('old_native', [[0, 60, 0]])
    old['evidence']['annevo'] = {'native_transcript_id': 'prior.t'}
    bridge = model('new_bridge', [[30, 150, 0]])
    independent = model('new_independent', [[120, 180, 0]])
    rescue.consolidate([old, bridge, independent], [], 'Target', retained_paths=retained_path(old))
    assert old['status'] == independent['status'] == 'accepted'
    assert bridge['status'] == 'unresolved' and bridge['problems'] == ['competing_new_models']
    assert bridge['integration_conflict']['disposition'] == 'retain_previous_validated_annotation'
    assert old['retained_coding_path']['model_ids'] == ['previous']
    assert len(old['support']) == 1  # History is not another donor or another alignment.


def test_existing_valid_annotation_does_not_require_a_new_unique_paralog_assignment():
    old, new = model('old', [[0, 60, 0]]), model('new', [[120, 180, 0]])
    prior = retained_path(old)
    for row in (old, new):
        row['problems'] = ['unanchored_genome_search']
    rescue.consolidate([old, new], [], 'Target', retained_paths=prior)
    assert old['status'] == 'accepted' and old['placement_evidence']['classification'] == 'retained_annotation'
    assert old['placement_evidence']['orthology'] == old['placement_evidence']['expected_copy'] == 'unassigned'
    assert new['status'] == 'unresolved' and new['problems'] == ['unanchored_genome_search']


def test_a_bad_new_path_does_not_disqualify_another_compatible_alternative():
    old = model('old', [[0, 90, 0]])
    compatible = model('compatible', [[6, 90, 0]])
    shifted = model('shifted', [[1, 91, 0]])
    rescue.consolidate([old, compatible, shifted], [], 'Target', retained_paths=retained_path(old))
    assert old['status'] == 'accepted' and compatible['status'] == 'accepted_alternative_path'
    assert shifted['status'] == 'unresolved'


def test_new_paths_that_compete_with_each_other_cannot_split_a_previous_locus():
    old = model('old', [[0, 180, 0]])
    a, b = model('a', [[0, 90, 0]]), model('b', [[90, 180, 0]])
    rescue.consolidate([old, a, b], [], 'Target', retained_paths=retained_path(old))
    assert old['status'] == 'accepted'
    assert a['status'] == b['status'] == 'unresolved'


@pytest.mark.parametrize('problem', ['missing_start', 'missing_stop', 'internal_stop', 'low_coverage', 'invalid_phase'])
def test_previous_adoption_cannot_override_current_failed_cds_or_homology_gate(problem):
    row = model('old', [[0, 60, 0]])
    prior = retained_path(row)
    row['problems'] = ['unanchored_genome_search', problem]
    rescue.consolidate([row], [], 'Target', retained_paths=prior)
    assert row['status'] == 'unresolved' and problem in row['problems']
    assert 'retained_coding_path' not in row


@pytest.mark.parametrize('change', ['sequence', 'coordinates', 'ownership'])
def test_retention_requires_same_sequence_and_coordinates_and_respects_original_owner(change):
    row = model('old', [[0, 60, 0]])
    prior = retained_path(row)
    owners = []
    if change == 'sequence':
        row['sequence'] = 'C' + row['sequence'][1:]
    elif change == 'coordinates':
        row['cds'] = [[3, 63, 0]]
    else:
        owners = [{'seqid': 'chr1', 'strand': '+', 'start': 0, 'end': 60,
                   'gene_id': 'original', 'cds': [[0, 60]]}]
    row['problems'] = ['unanchored_genome_search']
    rescue.consolidate([row], owners, 'Target', retained_paths=prior)
    assert row['status'] == 'unresolved'
    assert 'overlap_existing_annotation' in row['problems'] if change == 'ownership' else 'retained_coding_path' not in row


@pytest.mark.parametrize('legacy', [False, True])
def test_read_retained_paths_binds_primary_and_alternative_to_exact_source_receipt(tmp_path, legacy):
    from workflow.support import rescue_model_store as store
    from workflow.tests.test_rescue_model_store import receipt
    long, short = validated_pair()
    for i, row in enumerate((long, short)):
        row['model_id'] = 'old_' + str(i)
    long['alternative_coding_paths'] = [short]
    directory = tmp_path / 'previous'
    directory.mkdir()
    if legacy:
        (directory / 'models.json').write_text(json.dumps([long]))
    else:
        store.write_model_store(directory, [long])
    receipt(directory)
    expected = rescue.digest(directory / 'receipt.json')
    proof = rescue.FreshDigestBatch()
    retained = rescue.read_retained_coding_paths(directory, expected, verification=proof)
    assert set(retained) == {paths.coding_shape(long), paths.coding_shape(short)}
    assert all(row['producer_receipt_sha256'] == expected for row in retained.values())
    with pytest.raises(ValueError, match='receipt differs'):
        rescue.read_retained_coding_paths(directory, '0' * 64)
    member = directory / ('models.json' if legacy else 'model_store/accepted.sqlite')
    member.touch()
    with pytest.raises(OSError):
        proof.check()



@pytest.mark.parametrize('reverse', [False, True])
def test_old_bridge_cannot_fuse_recent_annotation_and_independent_historical_locus(reverse):
    left = model('left', [[0, 60, 0]])
    bridge = model('old_bridge', [[0, 60, 0], [120, 180, 0]])
    right = model('old_right', [[120, 180, 0]])
    for row, depth in [(left, 0), (bridge, 2), (right, 2)]:
        row['retained_coding_path'] = {'history_depth': depth}
    rows = [left, bridge, right]
    paths.resolve_coding_paths(list(reversed(rows)) if reverse else rows)
    assert left['status'] == right['status'] == 'accepted'
    assert bridge['status'] == 'unresolved'
    assert bridge['integration_conflict']['retained_model_ids'] == ['left', 'old_right']
    assert bridge['retained_coding_path']['history_depth'] == 2
    assert all('alternative_coding_paths' not in row for row in (left, right))


def test_currently_valid_older_frame_conflict_stays_a_review_proposal():
    recent, older = model('recent', [[0, 60, 0]]), model('older', [[1, 61, 0]])
    recent['retained_coding_path'] = {'history_depth': 0}
    older['retained_coding_path'] = {'history_depth': 3}
    paths.resolve_coding_paths([older, recent])
    assert recent['status'] == 'accepted' and older['status'] == 'unresolved'
    assert older['integration_conflict']['retained_model_ids'] == ['recent']


def test_conflicting_same_generation_prior_paths_are_preserved_for_review():
    a, b = model('a', [[0, 60, 0]]), model('b', [[1, 61, 0]])
    for row in (a, b):
        row['retained_coding_path'] = {'history_depth': 2, 'source_directory': '/verified/old'}
    paths.resolve_coding_paths([a, b])
    assert all(row['status'] == 'unresolved' and row['sequence'] for row in (a, b))
    assert all(row['retained_coding_path']['source_directory'] == '/verified/old' for row in (a, b))


@pytest.mark.parametrize('strand', ['+', '-'])
def test_missing_prior_path_is_reconstructed_with_current_donor_then_respects_ownership(strand):
    from workflow.tests.test_rescue_boundary_revalidation import PARAMS, checked
    dna = 'ATG' + 'GCT' * 30 + 'TAA'
    original, genome = checked(dna, [[0, len(dna), 0]], 'M' + 'A' * 30, strand)
    original['model_id'] = 'previous'
    retained = retained_path(original)
    metadata = retained[paths.coding_shape(original)]
    metadata['donor_queries'] = [['Relative', 'gene'], ['Gone', 'removed']]
    before = copy.deepcopy(retained)
    proteins = {'Relative': {'gene': 'M' + 'A' * 30}}
    candidates, audit = rescue.revalidate_retained_coding_paths(retained, [], genome, 1, PARAMS, proteins, 'Target')
    assert len(candidates) == 1 and audit['rows'][0]['emitted'] == 1
    candidate = candidates[0]
    assert candidate['cds'] == original['cds'] and candidate['sequence'] == original['sequence']
    assert candidate['coverage'] == candidate['identity'] == 1
    assert candidate['evidence']['donor'] == 'Relative'
    rescue.consolidate(candidates, [], 'Target', retained_paths=retained)
    assert candidate['status'] == 'accepted'
    assert retained == before
    repeated, _ = rescue.revalidate_retained_coding_paths(retained, candidates, genome, 1, PARAMS, proteins, 'Target')
    assert repeated == []
    owners = [{'seqid': 'chr', 'strand': strand, 'start': 0, 'end': len(dna),
               'gene_id': 'original_owner', 'cds': [[0, len(dna)]]}]
    fresh, _ = rescue.revalidate_retained_coding_paths(retained, [], genome, 1, PARAMS, proteins, 'Target')
    rescue.consolidate(fresh, owners, 'Target', retained_paths=retained)
    assert fresh[0]['status'] == 'unresolved' and fresh[0]['revision_owner_ids'] == ['original_owner']


@pytest.mark.parametrize('defect', ['DNA_change', 'internal_stop', 'donor_change', 'donor_missing', 'cell_budget'])
def test_historical_status_never_overrides_current_structure_or_homology(defect):
    from workflow.tests.test_rescue_boundary_revalidation import PARAMS, checked
    dna = 'ATG' + 'GCT' * 30 + 'TAA'
    old, genome = checked(dna, [[0, len(dna), 0]], 'M' + 'A' * 30)
    old['model_id'] = 'previous'
    retained = retained_path(old)
    retained[paths.coding_shape(old)]['donor_queries'] = [['Relative', 'gene']]
    donor = 'M' + 'A' * 30
    params = dict(PARAMS)
    if defect == 'DNA_change':
        genome.dna = genome.dna.replace('GCT', 'GCC', 1)
    elif defect == 'internal_stop':
        genome.dna = genome.dna.replace('GCT', 'TAA', 1)
    elif defect == 'donor_change':
        donor = 'M' + 'D' * 30
    elif defect == 'donor_missing':
        donor = ''
    else:
        params['terminal_max_alignment_cells'] = 10
    candidates, audit = rescue.revalidate_retained_coding_paths(
        retained, [], genome, 1, params, {'Relative': {'gene': donor}}, 'Target')
    assert candidates == [] and audit['rows'][0]['emitted'] == 0


def test_fixed_previous_path_homology_does_not_require_unique_boundary_traceback():
    from workflow.tests.test_rescue_boundary_revalidation import PARAMS, checked
    dna = 'ATG' + 'GCT' * 30 + 'TAA'
    donor = 'M' + 'A' * 28
    original, genome = checked(dna, [[0, len(dna), 0]], donor)
    original['model_id'] = 'previous'
    retained = retained_path(original)
    retained[paths.coding_shape(original)]['donor_queries'] = [['Relative', 'gene']]
    candidates, _ = rescue.revalidate_retained_coding_paths(
        retained, [], genome, 1, {**PARAMS, 'terminal_max_optimal_alignments': 1},
        {'Relative': {'gene': donor}}, 'Target')
    assert len(candidates) == 1
    assert candidates[0]['cds'] == original['cds'] and candidates[0]['coverage'] == 1


def retained_history_publication(root, records, parent=None, *, genome='same genome', legacy=False):
    from workflow.support import rescue_model_store as store
    from workflow.tests.test_rescue_model_store import receipt
    directory = root / 'rescued' / 'Target'
    directory.mkdir(parents=True)
    genome_path = root / 'genome.fa'
    genome_path.write_text(genome)
    request = {'sources': {'Target': {'genome': str(genome_path)}},
               'files': {str(genome_path): rescue.digest(genome_path)}}
    if parent is not None:
        request['prediction_cache'] = {'root': str(parent.parent.parent),
            'plan_sha256': rescue.digest(parent.parent.parent / 'plan.json'),
            'species': {'Target': {'receipt_sha256': rescue.digest(parent / 'receipt.json')}}}
    (root / 'plan.json').write_text(json.dumps({'request': request}))
    if legacy:
        (directory / 'models.json').write_text(json.dumps(records))
    else:
        store.write_model_store(directory, records)
    receipt(directory)
    value = json.loads((directory / 'receipt.json').read_text())
    value['key']['plan'] = rescue.digest(root / 'plan.json')
    (directory / 'receipt.json').write_text(json.dumps(value))
    return directory


@pytest.mark.parametrize('legacy', [False, True])
def test_retained_history_restores_a_path_lost_before_the_immediate_cache(tmp_path, legacy):
    long, short = validated_pair()
    long['model_id'], short['model_id'] = 'old_long', 'old_short'
    oldest = retained_history_publication(tmp_path / 'oldest', [long, short], legacy=legacy)
    middle = retained_history_publication(tmp_path / 'middle', [short], oldest, legacy=legacy)
    current = copy.deepcopy(short)
    current['model_id'] = 'current_short'
    latest = retained_history_publication(tmp_path / 'latest', [current], middle, legacy=legacy)
    proof, audit = rescue.FreshDigestBatch(), []
    retained = rescue.read_retained_coding_paths(latest, rescue.digest(latest / 'receipt.json'),
                                               verification=proof, history_audit=audit)
    assert set(retained) == {paths.coding_shape(long), paths.coding_shape(short)}
    assert retained[paths.coding_shape(long)]['source_directory'] == str(oldest)
    assert retained[paths.coding_shape(short)]['model_ids'] == ['current_short']
    assert len(audit) == 3 and all(row['status'] == 'verified' for row in audit)
    candidate = copy.deepcopy(long)
    candidate['problems'] = ['outside_expected_synteny_interval']
    rescue.consolidate([candidate], [], 'Target', retained_paths=retained)
    assert candidate['status'] == 'accepted'
    assert candidate['retained_coding_path']['producer_receipt_sha256'] == rescue.digest(oldest / 'receipt.json')
    failed = copy.deepcopy(long)
    failed['problems'] = ['low_coverage', 'outside_expected_synteny_interval']
    rescue.consolidate([failed], [], 'Target', retained_paths=retained)
    assert failed['status'] == 'unresolved' and 'retained_coding_path' not in failed
    member = oldest / ('models.json' if legacy else 'model_store/accepted.sqlite')
    member.touch()
    with pytest.raises(OSError):
        proof.check()


def test_retained_history_does_not_cross_genome_changes(tmp_path):
    long, short = validated_pair()
    long['model_id'], short['model_id'] = 'old_long', 'current_short'
    oldest = retained_history_publication(tmp_path / 'oldest', [long], genome='old genome')
    latest = retained_history_publication(tmp_path / 'latest', [short], oldest, genome='new genome')
    audit = []
    retained = rescue.read_retained_coding_paths(latest, rescue.digest(latest / 'receipt.json'), history_audit=audit)
    assert set(retained) == {paths.coding_shape(short)}
    assert audit[-1]['status'] == 'incompatible_genome_history'


@pytest.mark.parametrize('member', ['plan.json', 'accepted'])
def test_available_retained_history_rejects_corruption_instead_of_skipping_it(tmp_path, member):
    long, short = validated_pair()
    long['model_id'], short['model_id'] = 'old_long', 'current_short'
    oldest = retained_history_publication(tmp_path / 'oldest', [long])
    latest = retained_history_publication(tmp_path / 'latest', [short], oldest)
    path = oldest.parent.parent / 'plan.json' if member == 'plan.json' else oldest / 'model_store/accepted.sqlite'
    path.write_bytes(path.read_bytes() + b'changed')
    with pytest.raises(ValueError, match='changed'):
        rescue.read_retained_coding_paths(latest, rescue.digest(latest / 'receipt.json'))


def test_unavailable_archived_history_is_reported_without_consuming_it(tmp_path):
    long, short = validated_pair()
    long['model_id'], short['model_id'] = 'old_long', 'current_short'
    oldest = retained_history_publication(tmp_path / 'oldest', [long])
    latest = retained_history_publication(tmp_path / 'latest', [short], oldest)
    (oldest / 'receipt.json').unlink()
    audit = []
    retained = rescue.read_retained_coding_paths(latest, rescue.digest(latest / 'receipt.json'), history_audit=audit)
    assert set(retained) == {paths.coding_shape(short)}
    assert audit[-1] == {'directory': str(oldest), 'status': 'unavailable_history'}


def test_unbound_legacy_plan_cannot_extend_an_explicitly_verified_publication(tmp_path):
    long, short = validated_pair()
    long['model_id'], short['model_id'] = 'old_long', 'current_short'
    oldest = retained_history_publication(tmp_path / 'oldest', [long])
    latest = retained_history_publication(tmp_path / 'latest', [short], oldest)
    path = latest / 'receipt.json'
    value = json.loads(path.read_text())
    value['key'].pop('plan')
    path.write_text(json.dumps(value))
    audit = []
    retained = rescue.read_retained_coding_paths(latest, rescue.digest(path), history_audit=audit)
    assert set(retained) == {paths.coding_shape(short)}
    assert audit[-1]['status'] == 'unbound_plan_history'


def test_cyclic_retained_history_is_rejected(tmp_path):
    long, _ = validated_pair()
    long['model_id'] = 'old_long'
    directory = retained_history_publication(tmp_path / 'cycle', [long])
    plan = directory.parent.parent / 'plan.json'
    value = json.loads(plan.read_text())
    value['request']['prediction_cache'] = {'root': str(directory.parent.parent), 'plan_sha256': '0' * 64,
                                          'species': {'Target': {'receipt_sha256': '0' * 64}}}
    plan.write_text(json.dumps(value))
    path = directory / 'receipt.json'
    receipt = json.loads(path.read_text())
    receipt['key']['plan'] = rescue.digest(plan)
    path.write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match='Cyclic'):
        rescue.read_retained_coding_paths(directory, rescue.digest(path))


def test_duplicate_queries_do_not_inflate_species_and_self_cannot_supply_external_support():
    a = model('a', [[0, 60, 0]], ('Relative',))
    b = model('b', [[6, 60, 0]], ('Relative',))
    b['support'].extend([{'donor': 'Relative', 'query': 'isoform2',
                         'alignment': {'coverage': 1., 'identity': 1., 'problems': []}},
                        {'donor': 'Weak', 'query': 'isoform',
                         'alignment': {'coverage': .3, 'identity': 1., 'problems': ['low_coverage']}},
                        {'donor': 'Invalid', 'query': 'isoform',
                         'alignment': {'coverage': True, 'identity': 1., 'problems': []}}])
    paths.resolve_coding_paths([a, b])
    assert a['status'] == 'accepted' and b['status'] == 'accepted_alternative_path'
    assert a['locus_support']['independent_donor_species'] == ['Relative']
    a = model('a', [[0, 60, 0]], ('Self',))
    b = model('b', [[6, 60, 0]], ('Relative',))
    for row in (a, b):
        row['evidence']['target'] = 'Self'
    paths.resolve_coding_paths([a, b])
    assert a['status'] == 'accepted' and b['status'] == 'accepted_alternative_path'
    a = model('a', [[0, 60, 0]], ('Self',))
    b = model('b', [[6, 60, 0]], ('Self',))
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
    assert a['status'] == 'accepted' and b['status'] == 'unresolved'
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
    (worker / 'receipt.json').write_text(json.dumps({
        'key': {'plan': 'frozen-plan', 'species': name},
        'files': {'models.json': rescue.digest(worker / 'models.json')},
    }))
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
