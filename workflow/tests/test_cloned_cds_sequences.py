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
