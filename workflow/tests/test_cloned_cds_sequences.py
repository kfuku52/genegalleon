import csv
import gzip
import hashlib
import importlib
import json
import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
cloned = importlib.import_module('cloned_cds_models')
fsi = importlib.import_module('format_species_inputs')


def fixture(root, *, longer=False, transcriptome=False, ambiguity=False):
    root.mkdir(parents=True, exist_ok=True)
    rng = random.Random(31)
    sequence = 'ATG' + ''.join(rng.choice(('GCT', 'GGT', 'CTG', 'TTC', 'AAC', 'AGC')) for _ in range(300)) + 'TAA'
    source = root / 'source.fa'
    if transcriptome:
        source.write_text('>old\nCCC' + sequence + 'CCC\n')
    elif longer:
        source.write_text('>old\n' + sequence + '\n')
    else:
        source.write_text('>N_fragment\n' + sequence[:300] + '\n>C_fragment\n' + sequence[300:] + '\n')
    if ambiguity:
        changed = sequence[:400] + ('C' if sequence[400] != 'C' else 'A') + sequence[401:]
        if ambiguity == 'insertion':
            changed = sequence[:400] + 'GCG' + sequence[400:]
        elif ambiguity == 'deletion':
            changed = sequence[:400] + sequence[403:]
        assert changed != sequence
        source.write_text('>old1\n' + sequence + '\n>old2\n' + changed + '\n')
    incoming = root / 'clones.fa'
    incoming.write_text('>clone\n' + sequence + '\n')
    task = dict(provider='local', species_key='Test_species', species_prefix='Test_species', cds_path=source,
                gff_path=None, genome_path=None, genetic_code=1, gene_grouping_mode='strict', gff_repair_mode='safe',
                cds_input_type='transcriptome' if transcriptome else 'cds')
    row = dict(cloned_cds_fasta=str(incoming), cloned_cds_sha256=hashlib.sha256(incoming.read_bytes()).hexdigest(),
               cloned_cds_mode='sequence', cds_input_type=task['cds_input_type'])
    if longer:
        row['cloned_cds_longer_only_ids'] = '["clone"]'
    cloned.attach(task, row, root)
    return task, sequence, source, incoming, row


def test_sequence_search_replaces_both_fragments_with_one_clone_and_retains_audit(tmp_path):
    task, sequence, source, incoming, _ = fixture(tmp_path / 'inputs')
    before = source.read_bytes(), incoming.read_bytes()
    cloned.apply_to_task(task, tmp_path / 'curated')
    audit = json.loads(Path(task['cloned_model_audit']).read_text())
    assert audit['removed_CDS_ids'] == ['C_fragment', 'N_fragment']
    assert audit['decisions'][0]['status'] == 'replace'
    assert list(fsi.iter_fasta_records(task['cds_path']))[0][1] == sequence
    assert task['gff_path'] is None and task['genome_path'] is None
    assert before == (source.read_bytes(), incoming.read_bytes())
    metadata = dict(cloned_cds_mode='sequence', cloned_cds_path=str(incoming),
                    cloned_model_audit_path=task['cloned_model_audit'], cloned_source_cds_path=str(task['cds_path']))
    evidence = cloned.evidence_paths(metadata)
    assert any(path.endswith('removed_sequences.fa') for path in evidence)
    Path(next(path for path in evidence if path.endswith('removed_sequences.fa'))).unlink()
    with pytest.raises(ValueError, match='Missing cloned sequence artifact'):
        cloned.evidence_paths(metadata)


def test_longer_only_equal_length_is_explicitly_retained(tmp_path):
    task, sequence, _, _, _ = fixture(tmp_path / 'inputs', longer=True)
    cloned.apply_to_task(task, tmp_path / 'curated')
    audit = json.loads(Path(task['cloned_model_audit']).read_text())
    assert audit['decisions'][0]['status'] == 'not_longer'
    assert audit['removed_CDS_ids'] == []
    assert list(fsi.iter_fasta_records(task['cds_path'])) == [('old', sequence)]


@pytest.mark.parametrize('longer_only', [False, True])
def test_complete_exact_clone_in_extended_source_honors_length_policy(tmp_path, longer_only):
    task, sequence, source, _, _ = fixture(tmp_path / 'inputs', longer=longer_only)
    source.write_text('>extended\nATG' + 'GCT' * 9 + sequence + '\n')
    cloned.apply_to_task(task, tmp_path / 'curated')
    audit = json.loads(Path(task['cloned_model_audit']).read_text())
    decision = audit['decisions'][0]
    assert decision['status'] == ('not_longer' if longer_only else 'replace')
    assert decision['matches'][0]['match_basis'] == 'complete_clone_exact'
    assert decision['matches'][0]['query_coverage'] == 1.0
    assert decision['matches'][0]['target_coverage'] < .98
    if not longer_only:
        assert list(fsi.iter_fasta_records(task['cds_path']))[0][1] == sequence


def test_multiple_extended_models_are_retained_even_when_clone_is_exact(tmp_path):
    task, sequence, source, _, _ = fixture(tmp_path / 'inputs')
    source.write_text('>extended1\nATG' + 'GCT' * 9 + sequence + '\n>extended2\nATG' + 'AAC' * 9 + sequence + '\n')
    cloned.apply_to_task(task, tmp_path / 'curated')
    audit = json.loads(Path(task['cloned_model_audit']).read_text())
    assert audit['decisions'][0]['status'] == 'ambiguous_extended_targets'
    assert Path(task['cds_path']).read_text() == source.read_text()


@pytest.mark.parametrize('variant', [True, 'insertion', 'deletion'])
def test_divergent_overlapping_paralogs_are_preserved_for_review(tmp_path, variant):
    task, _, source, _, _ = fixture(tmp_path / 'inputs', ambiguity=variant)
    cloned.apply_to_task(task, tmp_path / 'curated')
    audit = json.loads(Path(task['cloned_model_audit']).read_text())
    assert audit['decisions'][0]['status'] == 'ambiguous_overlapping_targets'
    assert audit['removed_CDS_ids'] == []
    assert Path(task['cds_path']).read_text() == source.read_text()


def test_transcriptome_keeps_native_orf_processing_and_format_resume(tmp_path):
    task, sequence, _, _, _ = fixture(tmp_path / 'inputs', transcriptome=True)
    output = tmp_path / 'formatted'
    output.mkdir()
    first = fsi.format_cds(task, output, False, False)
    path = first['output_path']
    assert list(fsi.iter_fasta_records(path))[0][1] == sequence
    assert Path(str(path) + '.orf-extraction.json').is_file()
    assert task['cds_input_type'] == 'transcriptome'
    assert fsi.format_cds(task, output, False, False)['status'] == 'skip'


def test_declared_mrna_cds_coordinates_are_extracted_before_intact_cds_validation(tmp_path):
    task, sequence, _, incoming, row = fixture(tmp_path / 'inputs')
    incoming.write_text('>clone CDS=9-' + str(8 + len(sequence)) + '\nCCCCCCCC' + sequence + 'CCCC\n')
    row['cloned_cds_sha256'] = hashlib.sha256(incoming.read_bytes()).hexdigest()
    cloned.attach(task, row, tmp_path)
    assert cloned.cloned_records(incoming, ['clone'], 1) == {'clone': sequence}
    incoming.write_text('>clone CDS=9-999999\n' + sequence + '\n')
    with pytest.raises(ValueError, match='span outside'):
        cloned.cloned_records(incoming, ['clone'], 1)


@pytest.mark.parametrize('field,value', [('gff_url', 'file:///reference.gff'),
                                        ('genome_url', 'file:///reference.fa')])
def test_sequence_route_cannot_override_paired_genomic_inputs(tmp_path, field, value):
    task, _, _, _, row = fixture(tmp_path / 'inputs')
    with pytest.raises(ValueError, match='paired genomic'):
        cloned.attach(task, {**row, field: value}, tmp_path)


@pytest.mark.parametrize('value', [0, '0', 'NaN', '1.01'])
def test_sequence_thresholds_reject_invalid_fractions(tmp_path, value):
    task, _, _, _, row = fixture(tmp_path / 'inputs')
    with pytest.raises(ValueError, match='greater than zero'):
        cloned.attach(task, {**row, 'cloned_cds_min_identity': value}, tmp_path)


def test_genomic_route_cannot_silently_ignore_sequence_thresholds(tmp_path):
    task, _, _, _, row = fixture(tmp_path / 'inputs')
    with pytest.raises(ValueError, match='thresholds require'):
        cloned.attach(task, {**row, 'cloned_cds_mode': 'genomic', 'cloned_cds_min_identity': '.98'}, tmp_path)


def test_search_binding_rejects_changed_policies_and_original_sources(tmp_path):
    task, _, source, _, row = fixture(tmp_path / 'inputs')
    cloned.apply_to_task(task, tmp_path / 'curated')
    original = Path(task['cloned_model_audit'])
    other, *_ = fixture(tmp_path / 'inputs')
    cloned.apply_to_task(other, tmp_path / 'curated')
    assert Path(other['cloned_model_audit']) == original
    changed, *_ = fixture(tmp_path / 'inputs')
    cloned.attach(changed, {**row, 'cloned_cds_min_identity': '1.0'}, tmp_path)
    cloned.apply_to_task(changed, tmp_path / 'curated')
    assert Path(changed['cloned_model_audit']) != original
    source.write_text(source.read_text() + '>new_gene\nATGAAACCCTAA\n')
    next_task = {**task, 'cds_path': source}
    next_task.pop('_cloned_models_applied')
    cloned.apply_to_task(next_task, tmp_path / 'curated')
    assert Path(next_task['cloned_model_audit']) != original


def test_gzip_source_basename_survives_native_format_and_resume(tmp_path):
    task, sequence, source, _, _ = fixture(tmp_path / 'inputs')
    compressed = source.with_suffix('.fa.gz')
    with gzip.open(compressed, 'wt') as stream:
        stream.write(source.read_text())
    task['cds_path'] = compressed
    cloned.apply_to_task(task, tmp_path / 'curated')
    assert Path(task['cds_path']).name == compressed.name
    assert list(fsi.iter_fasta_records(task['cds_path']))[0][1] == sequence
    output = tmp_path / 'formatted'
    output.mkdir()
    assert fsi.format_cds(task, output, False, False)['status'] == 'write'
    assert fsi.format_cds(task, output, False, False)['status'] == 'skip'


def test_competing_clones_do_not_arbitrarily_choose_a_transcript(tmp_path):
    task, sequence, source, incoming, row = fixture(tmp_path / 'inputs', longer=True)
    incoming.write_text(incoming.read_text() + '>other_clone\n' + sequence + '\n')
    row['cloned_cds_sha256'] = hashlib.sha256(incoming.read_bytes()).hexdigest()
    row.pop('cloned_cds_longer_only_ids')
    cloned.attach(task, row, tmp_path)
    cloned.apply_to_task(task, tmp_path / 'curated')
    audit = json.loads(Path(task['cloned_model_audit']).read_text())
    assert {r['status'] for r in audit['decisions']} == {'competing_clones'}
    assert Path(task['cds_path']).read_text() == source.read_text()


def test_no_hit_preserves_sources_and_explicit_alignment_evidence(tmp_path):
    task, _, source, _, _ = fixture(tmp_path / 'inputs')
    source.write_text('>unrelated\nATG' + 'A' * 300 + 'TAA\n')
    cloned.apply_to_task(task, tmp_path / 'curated')
    audit = json.loads(Path(task['cloned_model_audit']).read_text())
    assert audit['decisions'][0]['status'] == 'no_supported_match'
    assert Path(task['cds_path']).read_text() == source.read_text()
    assert '# 0 hits found' in (Path(task['cloned_model_audit']).parent / audit['outputs']['alignments']).read_text()
    assert 'removed' not in audit['outputs']


def test_search_cap_counts_even_ineligible_subjects_and_checks_query_length(tmp_path):
    module = importlib.import_module('cloned_cds_sequences')
    path = tmp_path / 'hits.tsv'
    targets = {'target' + str(i): dict(id=str(i), length=100) for i in range(1000)}
    rows = ['\t'.join(['clone1', key, '50', '100', '1', '100', '1', '100', '300', '100',
                       '100', '1e-30', 'A' * 100, 'A' * 100]) for key in targets]
    path.write_text('\n'.join(rows) + '\n')
    grouped, capped = module.search_candidates(path, {'clone1': 'query'}, targets, .98, .98, {'query': 300})
    assert grouped == {'query': {}} and capped == {'query'}
    with pytest.raises(ValueError, match='coordinates changed'):
        module.search_candidates(path, {'clone1': 'query'}, targets, .98, .98, {'query': 301})


def group_row(row, clones, sources, *, unplaced=False):
    return {**row, 'cloned_cds_ids': json.dumps(clones),
            'cloned_cds_sequence_groups': json.dumps({'candidate': dict(clone_ids=clones,
                source_cds_ids=sources, unplaced=unplaced)})}


@pytest.mark.parametrize('reverse', [False, True])
def test_reviewed_competitors_keep_one_representative_and_every_experimental_alternative(tmp_path, reverse):
    task, sequence, source, incoming, row = fixture(tmp_path / 'inputs', longer=True)
    changed = sequence[:450] + ('GCT' if sequence[450:453] != 'GCT' else 'GGT') + sequence[453:]
    source.write_text('>old\n' + sequence[-306:] + '\n')
    incoming.write_text('>clone\n' + sequence + '\n>other_clone\n' + changed + '\n')
    row['cloned_cds_sha256'] = hashlib.sha256(incoming.read_bytes()).hexdigest()
    row.pop('cloned_cds_longer_only_ids')
    ids = ['other_clone', 'clone'] if reverse else ['clone', 'other_clone']
    cloned.attach(task, group_row(row, ids, ['old']), tmp_path)
    before = source.read_bytes(), incoming.read_bytes()
    output = tmp_path / 'formatted'
    output.mkdir()
    result = fsi.format_cds(task, output, False, False)
    records = list(fsi.iter_fasta_records(result['output_path']))
    assert len(records) == 1 and records[0][1] == (changed if reverse else sequence)
    audit_path = Path(task['cloned_model_audit'])
    audit = json.loads(audit_path.read_text())
    decisions = {x['clone_id']: x for x in audit['decisions']}
    assert decisions[ids[0]]['status'] == 'replace'
    assert decisions[ids[1]]['status'] == 'alternative_clone'
    assert decisions[ids[1]]['representative_clone'] == ids[0]
    assert decisions[ids[0]]['gene_identity_status'] == 'unresolved_candidate_group'
    assert audit['removed_CDS_ids'] == ['old']
    assert dict((h.split()[0], s) for h, s in fsi.iter_fasta_records(audit_path.parent / audit['outputs']['received_cloned_cds'])) == {'clone': sequence, 'other_clone': changed}
    assert list(fsi.iter_fasta_records(audit_path.parent / audit['outputs']['alternative_cloned_cds']))[0][1] == (sequence if reverse else changed)
    assert before == (source.read_bytes(), incoming.read_bytes())
    assert fsi.format_cds(task, output, False, False)['status'] == 'skip'


def test_explicit_unplaced_clone_is_added_without_removing_unknown_models_or_inventing_a_genome(tmp_path):
    task, sequence, source, _, row = fixture(tmp_path / 'inputs')
    unrelated = 'ATG' + 'AAA' * 100 + 'TAA'
    source.write_text('>unrelated\n' + unrelated + '\n')
    cloned.attach(task, group_row(row, ['clone'], [], unplaced=True), tmp_path)
    cloned.apply_to_task(task, tmp_path / 'curated')
    records = list(fsi.iter_fasta_records(task['cds_path']))
    assert len(records) == 2 and [s for _, s in records] == [unrelated, sequence]
    audit_path = Path(task['cloned_model_audit'])
    audit = json.loads(audit_path.read_text())
    decision = audit['decisions'][0]
    assert decision['status'] == 'add_unplaced' and decision['search_status'] == 'no_supported_match'
    assert decision['genomic_origin_status'] == 'not_assessed_no_genome'
    assert decision['genomic_coordinates'] is None
    assert decision['source_catalog_status'] == 'no_supported_match_in_sequence_catalog'
    assert decision['gene_identity_status'] == 'unresolved_unplaced'
    assert audit['removed_CDS_ids'] == [] and 'removed' not in audit['outputs']
    assert 'alternative_cloned_cds' not in audit['outputs']
    assert all((audit_path.parent / name).stat().st_size for name in audit['outputs'].values())
    assert task['genome_path'] is None and task['gff_path'] is None


@pytest.mark.parametrize('unplaced', [False, True])
def test_all_clones_remain_in_analysis_with_unknown_copy_number_and_unrelated_cds_unchanged(tmp_path, unplaced):
    task, sequence, source, incoming, row = fixture(tmp_path / 'inputs')
    changed = sequence[:450] + ('GCT' if sequence[450:453] != 'GCT' else 'GGT') + sequence[453:]
    neighbour = 'ATG' + 'AAA' * 100 + 'TAA'
    source.write_text('>neighbour\n' + neighbour + '\n' + ('' if unplaced else '>old\n' + sequence[-306:] + '\n'))
    incoming.write_text('>clone\n' + sequence + '\n>other_clone\n' + changed + '\n')
    row['cloned_cds_sha256'] = hashlib.sha256(incoming.read_bytes()).hexdigest()
    row = group_row(row, ['clone', 'other_clone'], [] if unplaced else ['old'], unplaced=unplaced)
    groups = json.loads(row['cloned_cds_sequence_groups'])
    groups['candidate']['analysis_policy'] = 'all_clones'
    row['cloned_cds_sequence_groups'] = json.dumps(groups)
    cloned.attach(task, row, tmp_path)
    before = source.read_bytes(), incoming.read_bytes()
    output = tmp_path / 'formatted'
    output.mkdir()
    result = fsi.format_cds(task, output, False, False)
    records = {h.split()[0]: s for h, s in fsi.iter_fasta_records(result['output_path'])}
    assert len(records) == 3 and set(records.values()) == {neighbour, sequence, changed}
    audit_path = Path(task['cloned_model_audit'])
    audit = json.loads(audit_path.read_text())
    assert audit['removed_CDS_ids'] == ([] if unplaced else ['old'])
    assert [d['status'] for d in audit['decisions']] == ['add_unplaced' if unplaced else 'replace', 'retain_unresolved_clone']
    for d in audit['decisions']:
        assert d['analysis_included'] and d['copy_number_status'] == 'unresolved'
        assert d['confirmed_gene_copy_count'] is None and d['gene_identity_status'] == 'unresolved_candidate_group'
        assert records[d['formatted_output_id']] == (sequence if d['clone_id'] == 'clone' else changed)
        assert d['genomic_coordinates'] is None
    group = audit['reviewed_sequence_groups'][0]
    assert group['analysis_clones'] == ['clone', 'other_clone'] and group['confirmed_gene_copy_count'] is None
    members = list(csv.DictReader((audit_path.parent / audit['outputs']['sequence_group_members']).open(), delimiter='\t'))
    assert {m['formatted_output_id'] for m in members} <= records.keys()
    assert all(m['analysis_included'] == 'True' and m['confirmed_gene_copy_count'] == '' for m in members)
    assert list(fsi.iter_fasta_records(audit_path.parent / audit['outputs']['alternative_cloned_cds']))[0][1] == changed
    assert before == (source.read_bytes(), incoming.read_bytes())


@pytest.mark.parametrize('groups,selected', [
    ({'group': dict(clone_ids=['clone'], source_cds_ids=[], unplaced=True)}, []),
    ({'group': dict(clone_ids=['clone'], source_cds_ids=[], unplaced=True)}, ['other']),
    ({'group': dict(clone_ids=['clone'], source_cds_ids=[], unplaced=False)}, ['clone']),
    ({'group': dict(clone_ids=['clone'], source_cds_ids=['old'], unplaced=True)}, ['clone']),
    ({'group': dict(clone_ids=['clone'], source_cds_ids=[], unplaced=1)}, ['clone']),
    ({'group': dict(clone_ids=['clone', 'clone'], source_cds_ids=['old'], unplaced=False)}, ['clone']),
    ({'group': dict(clone_ids=['clone'], source_cds_ids=['old'], unplaced=False),
      'other': dict(clone_ids=['other'], source_cds_ids=['old'], unplaced=False)}, ['clone', 'other']),
    ({'bad\nname': dict(clone_ids=['clone'], source_cds_ids=[], unplaced=True)}, ['clone']),
    ({'group': dict(clone_ids=['clone'], source_cds_ids=[], unplaced=True, analysis_policy='infer_alleles')}, ['clone']),
    ({'group': dict(clone_ids=['clone'], source_cds_ids=[], unplaced=True, analysis_policy=[])}, ['clone']),
])
def test_reviewed_sequence_group_schema_rejects_incomplete_or_overlapping_ownership(groups, selected):
    module = importlib.import_module('cloned_cds_sequences')
    with pytest.raises(ValueError, match='(Reviewed|Invalid)'):
        module.reviewed_groups({'cloned_cds_sequence_groups': json.dumps(groups)}, selected)


@pytest.mark.parametrize('status', ['search_limit_reached', 'not_longer', 'ambiguous_alignment',
                                    'ambiguous_extended_targets', 'ambiguous_overlapping_targets'])
def test_reviewed_group_does_not_override_unsafe_search_or_length_evidence(status):
    module = importlib.import_module('cloned_cds_sequences')
    decisions = [dict(clone_id='clone', status=status, matches=[dict(target_id='old')])]
    with pytest.raises(ValueError, match='current search evidence'):
        module.apply_reviewed_groups({'clone': 'ATGAAATAA'}, decisions,
            {'group': dict(clone_ids=['clone'], source_cds_ids=['old'], unplaced=False)}, {'old'})


def test_reviewed_unplaced_group_rejects_an_existing_supported_hit(tmp_path):
    task, _, _, _, row = fixture(tmp_path / 'inputs')
    cloned.attach(task, group_row(row, ['clone'], [], unplaced=True), tmp_path)
    with pytest.raises(ValueError, match='current search evidence'):
        cloned.apply_to_task(task, tmp_path / 'curated')


def test_reviewed_sequence_groups_cannot_be_used_for_genomic_inputs(tmp_path):
    task, _, _, _, row = fixture(tmp_path / 'inputs')
    with pytest.raises(ValueError, match='require cloned_cds_mode=sequence'):
        cloned.attach(task, {**group_row(row, ['clone'], [], unplaced=True), 'cloned_cds_mode': 'genomic'}, tmp_path)


def test_unplaced_sequence_cannot_bypass_a_longer_only_comparison(tmp_path):
    task, _, source, _, row = fixture(tmp_path / 'inputs', longer=True)
    source.write_text('>unrelated\nATG' + 'AAA' * 100 + 'TAA\n')
    cloned.attach(task, group_row(row, ['clone'], [], unplaced=True), tmp_path)
    with pytest.raises(ValueError, match='longer-only comparison'):
        cloned.apply_to_task(task, tmp_path / 'curated')


def test_reviewed_group_identity_and_representative_order_are_bound_to_the_native_plan(tmp_path):
    task, sequence, _, incoming, row = fixture(tmp_path / 'inputs')
    incoming.write_text(incoming.read_text() + '>other\n' + sequence + '\n')
    row['cloned_cds_sha256'] = hashlib.sha256(incoming.read_bytes()).hexdigest()
    cloned.attach(task, group_row(row, ['clone', 'other'], ['N_fragment', 'C_fragment']), tmp_path)
    first = cloned.binding_signature(task)
    forwarded = cloned.copy_binding({}, task)
    assert forwarded['cloned_cds_sequence_groups'] == task['cloned_cds_sequence_groups']
    cloned.attach(task, group_row(row, ['other', 'clone'], ['N_fragment', 'C_fragment']), tmp_path)
    assert cloned.binding_signature(task) != first
    reverse = cloned.binding_signature(task)
    changed = group_row(row, ['other', 'clone'], ['N_fragment', 'C_fragment'])
    groups = json.loads(changed['cloned_cds_sequence_groups'])
    groups['candidate']['analysis_policy'] = 'all_clones'
    changed['cloned_cds_sequence_groups'] = json.dumps(groups)
    cloned.attach(task, changed, tmp_path)
    assert cloned.binding_signature(task) != reverse
    cloned.attach(task, row, tmp_path)
    assert 'cloned_cds_sequence_groups' not in task
