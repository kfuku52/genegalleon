import json
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
import cds_resolution as resolution

pytestmark = pytest.mark.fast


def test_phase_is_constraint_and_padding_does_not_mask_stops():
    assert resolution.evaluate('ATGAAATAA')['accepted']
    assert not resolution.evaluate('ATGTAATAA', phase=0)['accepted']
    rescued = resolution.evaluate('AATGAAATAA', phase=1)
    assert rescued['accepted'] and rescued['head_padding'] == 2
    assert rescued['sequence'] == 'NNAATGAAATAA'
    assert resolution.evaluate('ATGNNNTAA')['reason'] == 'internal_ambiguous_bases'
    # A stop before tail padding remains internal; masking is never a repair.
    assert resolution.evaluate('ATGTAAA', phase=0)['reason'] == 'premature_stop'


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


def test_negative_strand_reference_and_explicit_fragment_order(tmp_path):
    genome = tmp_path / 'genome.fa'
    genome.write_text('>a\nTTACCC\n>b\nCAT\n')
    traits = pd.DataFrame([{'gene_id': 'g', 'feature_blocks': '1-3;1-3',
                           'feature_block_sequences': 'b;a', 'feature_block_strands': '-;-'}])
    assert resolution.extract_genomic_candidates(traits, genome) == {'g': 'ATGTAA'}


def test_resolution_traits_are_bound_to_exact_selected_sequence(tmp_path):
    from gff2genestat import apply_cds_resolution
    fasta, gff, ref = inputs(tmp_path)
    output, traits_path, _ = resolution.resolve(fasta, gff, ref, tmp_path / 'resolved')
    traits = pd.read_csv(traits_path, sep='\t')
    records = list(resolution.fasta_records(output))
    apply_cds_resolution(traits, records, output.parent)
    assert traits.iloc[0].structure_status == 'sequence_verified'
    records[0] = (records[0][0], records[0][1], 'ATGCCCTAA')
    apply_cds_resolution(traits, records, output.parent)
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
