import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
import cds_resolution as resolution

pytestmark = pytest.mark.fast


def resolved_view_fixture(tmp_path):
    sources = tmp_path / 'sources'
    sources.mkdir()
    source = sources / 'Plant_species.fa'
    source.write_text('>Plant_species_gene1\nATGAAATAA\n')
    target, _, _ = resolution.resolve(source, None, None, tmp_path / 'resolved')
    return sources, source, target, tmp_path / 'views'


def test_warm_view_hashes_resolved_target_once_per_invocation(tmp_path, monkeypatch):
    sources, _source, target, views = resolved_view_fixture(tmp_path)
    initial = resolution.resolved_view(sources, target.parent, views)
    state = sys.modules[resolution.FreshDigestBatch.__module__]
    original = state.digest
    reads = []
    def digest(path):
        reads.append(Path(path))
        return original(path)
    monkeypatch.setattr(state, 'digest', digest)
    assert resolution.resolved_view(sources, target.parent, views) == initial
    assert reads.count(target) == 1
    resolution.resolved_view(sources, target.parent, views)
    assert reads.count(target) == 2


def test_warm_view_rejects_target_replacement_while_alias_is_checked(tmp_path, monkeypatch):
    sources, _source, target, views = resolved_view_fixture(tmp_path)
    initial = resolution.resolved_view(sources, target.parent, views)
    state = sys.modules[resolution.FreshDigestBatch.__module__]
    original = state.digest
    def changed(path):
        value = original(path)
        if Path(path).parent == initial:
            replacement = tmp_path / 'replacement'
            replacement.write_bytes(target.read_bytes())
            replacement.replace(target)
        return value
    monkeypatch.setattr(state, 'digest', changed)
    with pytest.raises(OSError, match='changed'):
        resolution.resolved_view(sources, target.parent, views)


@pytest.mark.parametrize('change', ['rewrite', 'replace'])
def test_cold_view_revalidates_after_hardlink_writes_before_publication(tmp_path, monkeypatch, change):
    sources, _source, target, views = resolved_view_fixture(tmp_path)
    original = resolution.os.link
    def changed(source, destination):
        original(source, destination)
        if change == 'rewrite':
            target.write_text('>Plant_species_gene1\nATGCCCTAA\n')
        else:
            replacement = tmp_path / 'replacement'
            replacement.write_bytes(target.read_bytes())
            replacement.replace(target)
    monkeypatch.setattr(resolution.os, 'link', changed)
    with pytest.raises(ValueError, match='changed while creating view'):
        resolution.resolved_view(sources, target.parent, views)
    assert not list(views.iterdir())


def test_cold_view_keeps_valid_physical_aliases_and_unresolved_inputs(tmp_path):
    sources = tmp_path / 'sources'
    sources.mkdir()
    source = sources / 'Species_one.fa'
    source.write_text('>gene\nATGAAATAA\n')
    os.link(source, sources / 'Species_two.fa')
    view = resolution.resolved_view(sources, tmp_path / 'no-resolution', tmp_path / 'views')
    assert sorted(path.name for path in view.iterdir()) == ['Species_one.fa', 'Species_two.fa']
    assert all(os.path.samefile(path, source) for path in view.iterdir())


@pytest.mark.parametrize('change', ['extra', 'symlink-member', 'symlink-root', 'concurrent-extra'])
def test_view_rejects_unrecorded_members_and_redirects(tmp_path, monkeypatch, change):
    sources, _source, target, views = resolved_view_fixture(tmp_path)
    if change == 'concurrent-extra':
        def competing_view(staged, destination):
            shutil.copytree(staged, destination)
            (destination / 'Unexpected_species.fa').write_text('>extra\nATGAAATAA\n')
            raise FileExistsError('another publisher won')
        monkeypatch.setattr(resolution.os, 'rename', competing_view)
    else:
        view = resolution.resolved_view(sources, target.parent, views)
        if change == 'extra':
            (view / 'Unexpected_species.fa').write_text('>extra\nATGAAATAA\n')
        elif change == 'symlink-member':
            member = view / target.name
            member.unlink()
            member.symlink_to(target)
        else:
            moved = tmp_path / 'redirected'
            view.rename(moved)
            view.symlink_to(moved, target_is_directory=True)
    with pytest.raises(ValueError, match='CDS.*view'):
        resolution.resolved_view(sources, target.parent, views)


def test_phase_is_constraint_and_padding_does_not_mask_stops():
    assert resolution.evaluate('ATGAAATAA')['accepted']
    assert not resolution.evaluate('ATGTAATAA', phase=0)['accepted']
    rescued = resolution.evaluate('AATGAAATAA', phase=1)
    assert rescued['accepted'] and rescued['head_padding'] == 2
    assert rescued['sequence'] == 'NNAATGAAATAA'
    assert resolution.evaluate('ATGNNNTAA')['reason'] == 'internal_ambiguous_bases'
    # A stop before tail padding remains internal; masking is never a repair.
    assert resolution.evaluate('ATGTAAA', phase=0)['reason'] == 'premature_stop'


@pytest.mark.parametrize('phase', [0, 1, 2])
@pytest.mark.parametrize('position, reason', [
    (0, 'internal_ambiguous_bases'),
    (3, 'internal_ambiguous_bases'),
    (6, 'accepted'),
])
def test_ambiguity_check_preserves_internal_and_terminal_codon_boundaries(phase, position, reason):
    source = 'A' * phase + 'CCCCCCTAA'
    for ambiguous in 'RYSWKMBDHVN':
        offset = phase + position
        sequence = source[:offset] + ambiguous + source[offset + 1:]
        result = resolution.evaluate(sequence, phase=phase)
        assert result['reason'] == reason, (ambiguous, phase, position)
        assert result['internal_stop_count'] == 0
        assert result['sequence'] == 'N' * ((3 - phase) % 3) + sequence
        assert result['source_sha256'] == resolution.sequence_digest(sequence)


@pytest.mark.parametrize('phase', [1, 2])
def test_ambiguity_in_first_incomplete_codon_is_retained(phase):
    for ambiguous in 'RYSWKMBDHVN':
        sequence = ambiguous + 'A' * (phase - 1) + 'CCCCCCTAA'
        result = resolution.evaluate(sequence, phase=phase)
        assert result['accepted']
        assert result['sequence'] == 'N' * ((3 - phase) % 3) + sequence


def test_cds_rejection_precedence_is_retained():
    assert resolution.evaluate('CCCTAANNNTAA', phase=0)['reason'] == 'premature_stop'
    assert resolution.evaluate('CCCTGATAA', code=27, phase=0)['reason'] == 'translation_uncertain'
    assert resolution.evaluate('NCCTGATAA', code=27, phase=0)['reason'] == 'internal_ambiguous_bases'


def test_valid_source_priority_and_conditional_longest_selection():
    choice, reason, _ = resolution.select_candidate('ATGAAATAA', 'ATGCCCAAATAA', phase=0)
    assert choice == 'supplied' and 'disagreement' in reason
    choice, reason, _ = resolution.select_candidate('ATGAAATAA', 'CCCATGAAATAA', phase=0)
    assert choice == 'gff_genome' and reason == 'longer_in_frame_extension'
    choice, reason, _ = resolution.select_candidate('', 'ATGAAATAA', phase=0)
    assert choice == 'gff_genome' and reason == 'gff_genome_rescue'


def inputs(tmp_path, cds='ATGAAATAA', genome='ATGAAATAA', utr=True):
    source = tmp_path / 'input'
    source.mkdir()
    fasta = source / 'Plant_species.fa'
    fasta.write_text('>Plant_species_gene1\n' + cds + '\n')
    ref = tmp_path / 'genome.fa'
    ref.write_text('>chr1\n' + genome + '\n')
    gff = tmp_path / 'Plant_species.gff'
    gff.write_text('chr1\ts\tgene\t1\t9\t.\t+\t.\tID=gene1\n'
                   'chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=tx1;Parent=gene1\n'
                   'chr1\ts\tCDS\t1\t9\t.\t+\t0\tParent=tx1\n' +
                   ('chr1\ts\tfive_prime_UTR\t1\t3\t.\t+\t.\tParent=tx1\n' if utr else ''))
    return fasta, gff, ref


def test_resolution_preserves_sources_and_disables_only_conflicting_utr(tmp_path):
    fasta, gff, ref = inputs(tmp_path)
    before = [p.read_bytes() for p in (fasta, gff, ref)]
    output, traits, report = resolution.resolve(fasta, gff, ref, tmp_path / 'resolved')
    assert [p.read_bytes() for p in (fasta, gff, ref)] == before
    assert output.read_text() == fasta.read_text()
    row = pd.read_csv(traits, sep='\t').iloc[0]
    assert row.structure_status == 'sequence_verified'
    assert row.utr_status.startswith('conflicting:')
    assert pd.isna(row.utr_blocks)
    assert row.feature_size == 9
    payload = json.loads(report.read_text())
    assert payload['accepted_count'] == 1
    mtime = output.stat().st_mtime_ns
    resolution.resolve(fasta, gff, ref, tmp_path / 'resolved')
    assert output.stat().st_mtime_ns == mtime
    view = resolution.resolved_view(fasta.parent, output.parent, tmp_path / 'views')
    assert (view / fasta.name).is_file() and not (view / fasta.name).is_symlink()
    assert resolution.digest(view / fasta.name) == resolution.digest(output)
    ref.write_text('>chr1\nATGCCCTAA\n')
    with pytest.raises(ValueError, match='Stale CDS resolution'):
        resolution.resolved_view(fasta.parent, output.parent, tmp_path / 'views')


def test_resolution_uses_explicit_short_model_and_invalidates_changed_map(tmp_path):
    fasta, gff, genome = inputs(tmp_path, genome='CCCATGAAATAA', utr=False)
    gff.write_text('chr1\ts\tgene\t1\t12\t.\t+\t.\tID=gene1\n'
                   'chr1\ts\tmRNA\t4\t12\t.\t+\t.\tID=short;Parent=gene1\n'
                   'chr1\ts\tmRNA\t1\t12\t.\t+\t.\tID=long;Parent=gene1\n'
                   'chr1\ts\tCDS\t4\t12\t.\t+\t0\tParent=short\n'
                   'chr1\ts\tCDS\t1\t12\t.\t+\t0\tParent=long\n')
    selection = tmp_path / 'map.tsv'
    selection.write_text('species\tgene_id\tcandidate_id\tsource_transcript_id\tstatus\tscore\tmargin\treason\n'
                         'Plant_species\tPlant_species_gene1\tc1\tshort\tselected\t1\t0.2\tconserved\n')
    output, traits, report = resolution.resolve(fasta, gff, genome, tmp_path / 'resolved',
                                                representative_map=selection)
    row = pd.read_csv(traits, sep='\t').iloc[0]
    assert row.gff_transcript_id == 'short' and row.feature_size == 9 and row.start == 4
    assert output.read_text() == fasta.read_text()
    previous = json.loads(report.read_text())['contract']
    selection.write_text(selection.read_text().replace('\tshort\t', '\tlong\t'))
    resolution.resolve(fasta, gff, genome, tmp_path / 'resolved', representative_map=selection)
    changed = json.loads(report.read_text())['contract']
    assert changed['sources']['representative_map']['sha256'] != previous['sources']['representative_map']['sha256']
    assert pd.read_csv(traits, sep='\t').iloc[0].gff_transcript_id == 'long'
    assert output.read_text() == '>Plant_species_gene1\nCCCATGAAATAA\n'


def test_original_resolution_cannot_override_explicit_selected_short_coordinates(tmp_path):
    from gff2genestat import apply_cds_resolution

    fasta, gff, genome = inputs(tmp_path, genome='CCCATGAAATAA', utr=False)
    gff.write_text('chr1\ts\tgene\t1\t12\t.\t+\t.\tID=gene1\n'
                   'chr1\ts\tmRNA\t4\t12\t.\t+\t.\tID=short;Parent=gene1\n'
                   'chr1\ts\tmRNA\t1\t12\t.\t+\t.\tID=long;Parent=gene1\n'
                   'chr1\ts\tCDS\t4\t12\t.\t+\t0\tParent=short\n'
                   'chr1\ts\tCDS\t1\t12\t.\t+\t0\tParent=long\n')
    selection = tmp_path / 'map.tsv'
    selection.write_text('species\tgene_id\tcandidate_id\tsource_transcript_id\tstatus\tscore\tmargin\treason\n'
                         'Plant_species\tPlant_species_gene1\tc1\tshort\tselected\t1\t0.2\tconserved\n')
    output, original_traits, _report = resolution.resolve(fasta, gff, genome, tmp_path / 'resolved')
    assert output.read_text() == '>Plant_species_gene1\nCCCATGAAATAA\n'
    assert pd.read_csv(original_traits, sep='\t').iloc[0].gff_transcript_id == 'long'
    traits_path = tmp_path / 'selected.tsv'
    command = [sys.executable, str(Path(__file__).resolve().parents[1] / 'support/gff2genestat.py'),
               '--dir_gff', str(gff.parent), '--feature', 'CDS', '--multiple_hits', 'longest',
               '--seqfile', str(fasta), '--validate-cds-length', '--representative-map', str(selection),
               '--outfile', str(traits_path)]
    conflict = subprocess.run([*command, '--cds-resolution-dir', str(output.parent)],
                              capture_output=True, text=True, timeout=30)
    assert conflict.returncode != 0
    assert '--representative-map cannot be combined with --cds-resolution-dir' in conflict.stderr
    assert not traits_path.exists()
    selected = subprocess.run(command, capture_output=True, text=True, timeout=30)
    assert selected.returncode == 0, selected.stdout + selected.stderr
    traits = pd.read_csv(traits_path, sep='\t')
    row = traits.iloc[0]
    assert (row.gff_transcript_id, row.feature_size, row.start, row.end) == ('short', 9, 4, 12)
    assert row.phase_status == 'consistent'
    # Reproduce the original bug: a valid report for the raw longest CDS would
    # erase the selected transcript's coordinates because its sequence differs.
    apply_cds_resolution(traits, list(resolution.fasta_records(fasta)), output.parent)
    assert traits.iloc[0].structure_status == 'sequence_not_coordinate_matched'
    assert pd.isna(traits.iloc[0].feature_size)


def test_valid_disagreeing_cds_kept_without_intron_assignment(tmp_path):
    fasta, gff, ref = inputs(tmp_path, cds='ATGCCCAAATAA')
    output, traits, report = resolution.resolve(fasta, gff, ref, tmp_path / 'resolved')
    assert output.read_text() == fasta.read_text()
    row = pd.read_csv(traits, sep='\t').iloc[0]
    assert row.structure_status == 'sequence_not_coordinate_matched'
    assert pd.isna(row.feature_size) and pd.isna(row.intron_positions)
    assert json.loads(report.read_text())['decisions'][0]['selected_source'] == 'supplied'


def test_excluded_gene_does_not_abort_other_genes():
    records = [('good', 'good', 'ATGAAATAA'), ('bad', 'bad', 'NNNNNNNNNNNN')]
    output, decisions = resolution.resolve_records(records, pd.DataFrame(columns=resolution.TRAIT_COLUMNS), {}, 1)
    assert [r[0] for r in output] == ['good']
    assert decisions[1]['reason'] == 'no_valid_cds'


@pytest.mark.parametrize('index', [
    pd.Index([7, 19, 43]),
    pd.Index(['first', 'other', 'last']),
    pd.MultiIndex.from_tuples([('a', 7), ('b', 19), ('a', 43)]),
], ids=['nonconsecutive', 'string', 'multiindex'])
@pytest.mark.parametrize('dtype', ['object', 'string', 'category'])
def test_record_resolution_keeps_index_labels_and_last_duplicate_id(index, dtype):
    traits = pd.DataFrame({
        'gene_id': pd.array(['gene', 'other', 'gene'], dtype=dtype),
        'phase_status': ['consistent'] * 3,
        'cds_first_phase': [1, 0, 0],
        'structure_status': ['pending'] * 3,
    }, index=index)
    original = traits.copy(deep=True)
    sequence = 'ATGAAATAA'
    output, decisions = resolution.resolve_records(
        [('gene', 'gene descriptive header', sequence)], traits, {'gene': sequence}, 1)

    assert output == [('gene', 'gene descriptive header', sequence)]
    assert decisions[0]['reason'] == 'sources_agree'
    assert decisions[0]['structure_status'] == 'sequence_verified'
    assert 'phase_status' not in decisions[0]
    pd.testing.assert_frame_equal(traits.iloc[:2], original.iloc[:2])
    assert traits.iloc[2]['structure_status'] == 'sequence_verified'
    assert traits.iloc[2]['phase_status'] == 'consistent'
    assert traits.iloc[2]['cds_first_phase'] == 0
    pd.testing.assert_index_equal(traits.index, original.index)


@pytest.mark.parametrize('columns', [[], resolution.TRAIT_COLUMNS], ids=['no-columns', 'native-columns'])
def test_record_resolution_keeps_empty_traits(columns):
    traits = pd.DataFrame(columns=columns)
    original = traits.copy(deep=True)
    output, decisions = resolution.resolve_records([('gene', 'gene', 'ATGAAATAA')], traits, {}, 1)
    assert output == [('gene', 'gene', 'ATGAAATAA')]
    assert decisions[0]['structure_status'] == 'unverified_source'
    pd.testing.assert_frame_equal(traits, original)


def test_record_resolution_preserves_legacy_numeric_row_key_coercion():
    # iterrows rounds this integer in an all-numeric frame. Such helper inputs
    # must retain their existing unmatched-ID behavior instead of being updated.
    identifier = 2 ** 53 + 1
    traits = pd.DataFrame({'gene_id': [identifier], 'phase_status': [0.0], 'cds_first_phase': [0.0]})
    original = traits.copy(deep=True)
    output, decisions = resolution.resolve_records([(identifier, 'numeric', 'ATGAAATAA')], traits, {}, 1)
    assert output == [(identifier, 'numeric', 'ATGAAATAA')]
    assert decisions[0]['structure_status'] == 'unverified_source'
    pd.testing.assert_frame_equal(traits, original)


def test_record_resolution_preserves_mixed_identifier_values():
    identifier = 2 ** 53 + 1
    traits = pd.DataFrame({
        'gene_id': ['gene', identifier, None],
        'phase_status': ['consistent'] * 3,
        'cds_first_phase': [0] * 3,
        'structure_status': ['pending'] * 3,
    }, index=[7, 19, 43])
    sequence = 'ATGAAATAA'
    records = [('gene', 'gene', sequence), (identifier, 'numeric', sequence)]
    output, decisions = resolution.resolve_records(records, traits, {'gene': sequence, identifier: sequence}, 1)
    assert output == records
    assert [entry['structure_status'] for entry in decisions] == ['sequence_verified'] * 2
    assert traits['structure_status'].tolist() == ['sequence_verified', 'sequence_verified', 'pending']
    assert traits['gene_id'].tolist() == ['gene', identifier, None]


@pytest.mark.parametrize('traits, exception, message', [
    (pd.DataFrame([['gene', 'gene']], columns=['gene_id', 'gene_id']), TypeError, 'unhashable'),
    (pd.DataFrame({'other': ['gene']}), AttributeError, 'gene_id'),
    (pd.DataFrame([['gene', 'consistent', 0]], columns=pd.MultiIndex.from_tuples([
        ('gene_id', 'value'), ('phase_status', 'value'), ('cds_first_phase', 'value'),
    ])), TypeError, 'unhashable'),
], ids=['duplicate-columns', 'missing-column', 'hierarchical-columns'])
def test_record_resolution_preserves_malformed_trait_errors(traits, exception, message):
    original = traits.copy(deep=True)
    with pytest.raises(exception, match=message):
        resolution.resolve_records([], traits, {}, 1)
    pd.testing.assert_frame_equal(traits, original)


def test_negative_strand_reference_and_explicit_fragment_order(tmp_path):
    genome = tmp_path / 'genome.fa'
    genome.write_text('>a\nTTACCC\n>b\nCAT\n')
    traits = pd.DataFrame([{'gene_id': 'g', 'feature_blocks': '1-3;1-3',
                           'feature_block_sequences': 'b;a', 'feature_block_strands': '-;-'}])
    assert resolution.extract_genomic_candidates(traits, genome) == {'g': 'ATGTAA'}


@pytest.mark.parametrize('indexed', [False, True])
def test_resolution_traits_are_bound_to_exact_selected_sequence(tmp_path, indexed):
    from gff2genestat import apply_cds_resolution
    fasta, gff, ref = inputs(tmp_path)
    output, traits_path, _ = resolution.resolve(fasta, gff, ref, tmp_path / 'resolved')
    traits = pd.read_csv(traits_path, sep='\t')
    records = list(resolution.fasta_records(output))
    apply_cds_resolution(traits, records, output.parent, tmp_path / 'table-cache' if indexed else None)
    assert traits.iloc[0].structure_status == 'sequence_verified'
    records[0] = (records[0][0], records[0][1], 'ATGCCCTAA')
    apply_cds_resolution(traits, records, output.parent, tmp_path / 'table-cache' if indexed else None)
    assert traits.iloc[0].structure_status == 'sequence_not_coordinate_matched'
    assert pd.isna(traits.iloc[0].feature_size)


def test_genetic_code_and_missing_reference_do_not_fabricate_structure(tmp_path):
    assert resolution.evaluate('ATGTGATAA', code=1, phase=0)['reason'] == 'premature_stop'
    assert resolution.evaluate('ATGTGATAA', code=4, phase=0)['accepted']
    fasta, gff, _ = inputs(tmp_path)
    _, traits_path, _ = resolution.resolve(fasta, gff, None, tmp_path / 'resolved')
    row = pd.read_csv(traits_path, sep='\t').iloc[0]
    assert row.structure_status == 'sequence_not_coordinate_matched'
    assert pd.isna(row.num_intron)


def test_policy_rebuild_preserves_manifest_and_previous_output(tmp_path):
    workspace = tmp_path / 'workspace'
    manifests = workspace / 'output/artifact_provenance/genome_annotation'
    manifests.mkdir(parents=True)
    output = workspace / 'output/species_gff_info/Plant_species_gff_info.tsv'
    output.parent.mkdir()
    output.write_text('original\n')
    manifest = manifests / 'Plant_species.gff_info.json'
    manifest.write_text(json.dumps({'family_id': 'Plant_species', 'outputs': [
        {'scope': 'workspace', 'artifact_type': 'file', 'path': str(output.relative_to(workspace))}]}))
    backup = resolution.backup_family_outputs(workspace, 'Plant_species')
    assert (backup / output.relative_to(workspace)).read_text() == 'original\n'
    assert (backup / manifest.relative_to(workspace)).read_bytes() == manifest.read_bytes()
    assert resolution.backup_family_outputs(workspace, 'Plant_species') == backup
    output.write_text('regenerated\n')
    assert (backup / output.relative_to(workspace)).read_text() == 'original\n'
    output.unlink()
    output.symlink_to(tmp_path / 'outside')
    (tmp_path / 'outside').write_text('outside')
    with pytest.raises(ValueError, match='Unsafe output'):
        resolution.backup_family_outputs(workspace, 'Plant_species')


def test_resolution_refuses_to_overwrite_source(tmp_path):
    fasta, gff, ref = inputs(tmp_path)
    original = fasta.read_bytes()
    with pytest.raises(ValueError, match='must not overwrite'):
        resolution.resolve(fasta, gff, ref, fasta.parent)
    assert fasta.read_bytes() == original


def test_padding_does_not_manufacture_an_in_frame_extension():
    selected, reason, candidates = resolution.select_candidate('ATGAAA', 'ATGAAAC', phase=0)
    assert all(value['accepted'] for value in candidates.values())
    assert len(candidates['gff_genome']['sequence']) == 9  # includes two artificial Ns
    assert selected == 'supplied' and reason == 'supplied_valid_unresolved_source_disagreement'


def test_explicit_trans_spliced_cds_does_not_depend_on_utr_order(tmp_path):
    from format_species_annotation.genbank import derive_cds_records_from_gff_and_genome
    gff = tmp_path / 'input.gff'
    gff.write_text('a\ts\tgene\t1\t6\t.\t+\t.\tID=g\n'
                   'a\ts\tmRNA\t1\t6\t.\t+\t.\tID=tx;Parent=g\n'
                   'a\ts\tfive_prime_UTR\t1\t3\t.\t+\t.\tParent=tx\n'
                   'a\ts\tCDS\t1\t6\t.\t+\t0\tParent=tx;exception=trans-splicing;part=1/2\n'
                   'b\ts\tCDS\t1\t3\t.\t+\t0\tParent=tx;exception=trans-splicing;part=2/2\n')
    genome = tmp_path / 'genome.fa'
    genome.write_text('>a\nATGAAA\n>b\nTAA\n')
    records = list(derive_cds_records_from_gff_and_genome(
        {'gff_path': gff, 'genome_path': genome, 'provider': 'direct', 'species_key': 'Plant_species'}))
    assert len(records) == 1 and records[0][1] == 'ATGAAATAA'


def test_valid_supplied_frame_does_not_inherit_conflicting_gff_phase(tmp_path):
    from gff2genestat import apply_cds_resolution
    fasta, gff, ref = inputs(tmp_path, utr=False)
    gff.write_text(gff.read_text().replace('\t+\t0\tParent=tx1', '\t+\t1\tParent=tx1'))
    output, traits_path, report_path = resolution.resolve(fasta, gff, ref, tmp_path / 'resolved')
    traits = pd.read_csv(traits_path, sep='\t')
    row = traits.iloc[0]
    assert row.structure_status == 'sequence_verified' and row.feature_size == 9
    assert row.phase_status == 'conflicting' and pd.isna(row.cds_first_phase)
    decision = json.loads(report_path.read_text())['decisions'][0]
    assert decision['selected_source'] == 'supplied'
    assert decision['phase_status'] == 'selected_frame_disagrees_with_gff'
    # Previously generated coordinate-matched tables also lose the invalid
    # phase when consumed by the updated downstream workflow.
    traits.loc[0, 'cds_first_phase'] = 1
    traits.loc[0, 'phase_status'] = 'consistent'
    traits.to_csv(traits_path, sep='\t', index=False)
    meta_path = Path(str(output) + '.resolution.meta.json')
    meta = json.loads(meta_path.read_text())
    meta['traits_sha256'] = resolution.digest(traits_path)
    meta_path.write_text(json.dumps(meta))
    apply_cds_resolution(traits, list(resolution.fasta_records(output)), output.parent)
    assert pd.isna(traits.iloc[0].cds_first_phase) and traits.iloc[0].phase_status == 'conflicting'


def resolution_application_fixture(tmp_path):
    """Real receipts and tables, without invoking scientific resolution."""
    sequence = 'ATGAAATAA'
    identifiers = [f'Plant_species_gene{i}' for i in range(3)]
    resolved = pd.DataFrame({
        'gene_id': identifiers, 'feature_size': [9] * 3, 'num_intron': [2] * 3,
        'intron_positions': ['3;6'] * 3, 'feature_blocks': ['1-3;5-7;9-11'] * 3,
        'utr_blocks': [''] * 3, 'chromosome': ['chr1'] * 3, 'strand': ['+'] * 3,
        'start': [1, 21, 41], 'end': [11, 31, 51], 'cds_first_phase': [0, 1, 2],
        'structure_status': ['sequence_verified'] * 3, 'phase_status': ['consistent'] * 3,
        'cds_sequence_sha256': [resolution.sequence_digest(sequence)] * 3,
    })
    source = tmp_path / 'Plant_species.fa'
    source.write_text(''.join(f'>{identifier}\n{sequence}\n' for identifier in identifiers))
    directory = tmp_path / 'resolved'
    directory.mkdir()
    table = directory / 'Plant_species.fa.traits.tsv'
    resolved.to_csv(table, sep='\t', index=False)
    receipt = directory / 'Plant_species.fa.resolution.meta.json'
    receipt.write_text(json.dumps({'traits_sha256': resolution.digest(table), 'contract': {
        'sources': {'cds': {'path': str(source), 'sha256': resolution.digest(source)}}}}))
    target = resolved.copy(deep=True)
    target['feature_size'] = 6
    target['structure_status'] = 'pending'
    target['phase_status'] = 'pending'
    records = [(identifier, identifier, sequence) for identifier in identifiers]
    return target, records, directory, table, receipt, source


def compare_resolution_application_with_scalar(target, records, directory, monkeypatch):
    import gff2genestat as genestat
    expected, actual = target.copy(deep=True), target.copy(deep=True)
    def invoke(frame):
        try:
            genestat.apply_cds_resolution(frame, records, directory)
        except Exception as error:
            return type(error), str(error)
        return None
    with monkeypatch.context() as context:
        context.setattr(genestat, '_batch_cds_resolution_traits', lambda *args: False)
        expected_error = invoke(expected)
    actual_error = invoke(actual)
    assert actual_error == expected_error
    pd.testing.assert_frame_equal(actual, expected, check_exact=True, check_dtype=True)
    for column in expected:
        if expected[column].dtype == object:
            assert [type(value) for value in actual[column]] == [type(value) for value in expected[column]]
    return actual, actual_error


@pytest.mark.parametrize('index', [pd.Index([7, 19, 43, 59]),
    pd.Index(['first', 'other', 'third', 'last']),
    pd.MultiIndex.from_tuples([('a', 7), ('b', 19), ('c', 43), ('a', 59)])],
    ids=['integer', 'string', 'multiindex'])
@pytest.mark.parametrize('dtype', ['object', 'string', 'category'])
def test_apply_resolution_keeps_last_duplicate_id_phase_and_mismatch(tmp_path, monkeypatch, index, dtype):
    target, records, directory, *_ = resolution_application_fixture(tmp_path)
    target = pd.concat([target, target.iloc[[0]]], ignore_index=True)
    target.index = index
    target['gene_id'] = target.gene_id.astype(dtype)
    records[1] = (*records[1][:2], 'ATGCCCTAA')
    actual, error = compare_resolution_application_with_scalar(target, records, directory, monkeypatch)
    assert error is None
    assert actual.iloc[0].structure_status == 'pending'
    assert actual.iloc[3].structure_status == 'sequence_verified'
    assert actual.iloc[1].structure_status == 'sequence_not_coordinate_matched'
    assert pd.isna(actual.iloc[1].feature_size)
    assert actual.iloc[2].phase_status == 'conflicting'
    assert pd.isna(actual.iloc[2].cds_first_phase)


@pytest.mark.parametrize('layout', ['boolean-index', 'duplicate-index', 'nullable-phase',
    'missing-status', 'missing-phase', 'malformed-hash', 'array-hash', 'array-status', 'object-block',
    'object-numeric', 'numeric-structure-column'])
def test_apply_resolution_preserves_unusual_cell_errors_and_partial_writes(tmp_path, monkeypatch, layout):
    import gff2genestat as genestat
    target, records, directory, table, *_ = resolution_application_fixture(tmp_path)
    resolved = pd.read_csv(table, sep='\t')
    order = list(set(target.gene_id) & set(resolved.gene_id))
    second = resolved.index[resolved.gene_id == order[1]][0]
    if layout == 'boolean-index':
        target = target.iloc[:2]
        target.index = [False, True]
    elif layout == 'duplicate-index':
        target.index = [0, 0, 2]
    elif layout == 'nullable-phase':
        target['cds_first_phase'] = target.cds_first_phase.astype('Int64')
    elif layout == 'missing-phase':
        target = target.drop(columns='cds_first_phase')
    elif layout == 'object-numeric':
        target['feature_size'] = target.feature_size.astype(object)
    elif layout == 'numeric-structure-column':
        resolved['intron_positions'] = 3
        target['intron_positions'] = 0
        records[0] = (*records[0][:2], 'ATGCCCTAA')
        monkeypatch.setattr(genestat, 'load_resolution_traits', lambda path: resolved.copy(deep=True))
    else:
        column = ('cds_sequence_sha256' if layout in {'malformed-hash', 'array-hash'} else
                  'utr_blocks' if layout == 'object-block' else 'structure_status')
        resolved[column] = resolved[column].astype(object)
        target[column] = target[column].astype(object)
        value = {'missing-status': pd.NA, 'malformed-hash': pd.NA,
                 'array-hash': np.array(['bad', 'hash']),
                 'array-status': np.array(['sequence_verified', 'other']),
                 'object-block': ['1-3', '5-7']}[layout]
        resolved.at[second, column] = value
        monkeypatch.setattr(genestat, 'load_resolution_traits', lambda path: resolved.copy(deep=True))
    actual, error = compare_resolution_application_with_scalar(target, records, directory, monkeypatch)
    if layout in {'missing-status', 'malformed-hash', 'array-hash', 'array-status'}:
        assert error is not None
        first = actual.loc[actual.gene_id == order[0]].iloc[0]
        assert first.structure_status == 'sequence_verified'
        assert first.feature_size == 9


@pytest.mark.parametrize('change', ['traits-and-source', 'source', 'traits-replace-during-read', 'receipt-hash'])
def test_apply_resolution_retains_fresh_proofs_and_verification_precedence(tmp_path, monkeypatch, change):
    import gff2genestat as genestat
    target, records, directory, table, receipt, source = resolution_application_fixture(tmp_path)
    if change in {'traits-and-source', 'source'}:
        source.write_text(source.read_text().replace('AAA', 'CCC'))
    if change == 'traits-and-source':
        table.write_text(table.read_text().replace('chr1', 'chr2'))
    elif change == 'receipt-hash':
        report = json.loads(receipt.read_text())
        report['traits_sha256'] = 'changed'
        receipt.write_text(json.dumps(report))
    elif change == 'traits-replace-during-read':
        original = genestat.load_resolution_traits
        def replaced(path):
            result = original(path)
            replacement = table.with_suffix('.replacement')
            replacement.write_bytes(table.read_bytes())
            replacement.replace(table)
            return result
        monkeypatch.setattr(genestat, 'load_resolution_traits', replaced)
    expected = 'Stale CDS resolution source' if change == 'source' else (
        'changed' if change == 'traits-replace-during-read' else 'Changed CDS resolution traits')
    with pytest.raises((ValueError, OSError), match=expected):
        genestat.apply_cds_resolution(target, records, directory)
    assert target.structure_status.tolist() == ['pending'] * 3
    assert target.feature_size.tolist() == [6] * 3


def test_apply_resolution_retains_readonly_array_failure(tmp_path):
    import gff2genestat as genestat
    target, records, directory, *_ = resolution_application_fixture(tmp_path)
    columns = {column: target[column].to_numpy(copy=True) for column in target}
    columns['feature_size'].flags.writeable = False
    target = pd.DataFrame(columns, copy=False)
    with pytest.raises(ValueError, match='read-only'):
        genestat.apply_cds_resolution(target, records, directory)
    assert target.structure_status.tolist() == ['pending'] * 3
