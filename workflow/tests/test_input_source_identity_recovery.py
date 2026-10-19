"""Published RNA aliases and partial complementary phases preserve source models."""
import gzip
from pathlib import Path

import pytest

from workflow.tests.test_format_species_inputs import load_module, normalisation_bundle


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('provided', [True, False])
def test_explicit_rna_name_repairs_self_parented_rna_only(tmp_path, strand, provided):
    module = load_module()
    task = normalisation_bundle(tmp_path, [{'blocks': [('ATGA', 0), ('AATAA', 2)], 'strand': strand}], provided=provided)
    source = task['gff_path'].read_text().replace('ID=t1;Parent=g1', 'ID=g1;Name=t1;Parent=g1')
    task['gff_path'].write_text(source)
    output = tmp_path / 'output'
    output.mkdir()
    cds = module.format_cds(task, output, False, False)
    assert dict(module.iter_fasta_records(cds['output_path'])) == {'Test_species_g1': 'ATGAAATAA'}
    gff = module.format_gff(task, output, False, False, formatted_cds_path=cds['output_path'])
    with gzip.open(gff['output_path'], 'rt') as handle:
        rows = [line for line in handle if not line.startswith('#')]
    rna = next(line for line in rows if '\tmRNA\t' in line)
    assert 'ID=t1;' in rna and 'Parent=g1' in rna and 'gg_original_id=g1' in rna
    assert all('Parent=t1' in line for line in rows if '\tCDS\t' in line)
    assert task['gff_path'].read_text() == source
    from format_species_annotation.source_identity import locus_identity_audit
    audit = locus_identity_audit(task)
    assert any(m.get('reason') == 'explicit_rna_name_parent' for m in audit['mappings'])
    from validate_cds_gff_mapping import validate_single_species
    result = validate_single_species(dict(index=1, species_prefix='Test_species', cds_file=cds['output_path'],
        gff_file=gff['output_path'], genome_file=task['genome_path'], strict=True), 10)
    assert result['ok'], result.get('error')


@pytest.mark.parametrize('mutation', ['different_axis', 'outside', 'shared_child', 'existing_id', 'no_child', 'duplicate_name'])
def test_rna_name_without_unique_structural_evidence_is_not_repaired(tmp_path, mutation):
    load_module()
    task = normalisation_bundle(tmp_path, [{'blocks': [('ATGA', 0), ('AATAA', 2)]}])
    source = task['gff_path'].read_text().replace('ID=t1;Parent=g1', 'ID=g1;Name=t1;Parent=g1')
    if mutation == 'different_axis':
        lines = source.splitlines(True)
        index = next(i for i, line in enumerate(lines) if '\tCDS\t' in line)
        lines[index] = 'other' + lines[index][lines[index].index('\t'):]
        source = ''.join(lines)
    elif mutation == 'outside':
        lines = source.splitlines(True)
        index = next(i for i, line in enumerate(lines) if '\tCDS\t' in line)
        fields = lines[index].rstrip().split('\t')
        fields[4] = '999'
        lines[index] = '\t'.join(fields) + '\n'
        source = ''.join(lines)
    elif mutation == 'shared_child':
        source = source.replace('Parent=t1', 'Parent=t1,other')
    elif mutation == 'existing_id':
        source += 'chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=t1\n'
    elif mutation == 'no_child':
        source = source.replace('Parent=t1', 'Parent=missing')
    else:
        source += ('chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=g2\n'
                   'chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=g2;Name=t1;Parent=g2\n')
    task['gff_path'].write_text(source)
    from format_species_annotation.source_identity import locus_identity_audit
    audit = locus_identity_audit(task)
    assert not any(m.get('reason') == 'explicit_rna_name_parent' for m in audit['mappings'])
    assert task['gff_path'].read_text() == source


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('provided', [True, False])
def test_partial_cds_normalises_complementary_phases_without_deleting_bases(tmp_path, strand, provided):
    module = load_module()
    task = normalisation_bundle(tmp_path, [{'blocks': [('ATGA', 0), ('AACC', 1)], 'strand': strand}], provided=provided)
    original = {name: Path(task[name]).read_bytes() for name in ['gff_path', 'genome_path']}
    output = tmp_path / 'output'
    output.mkdir()
    cds = module.format_cds(task, output, False, False)
    assert dict(module.iter_fasta_records(cds['output_path'])) == {'Test_species_g1': 'ATGAAACCN'}
    gff = module.format_gff(task, output, False, False, formatted_cds_path=cds['output_path'])
    with gzip.open(gff['output_path'], 'rt') as handle:
        rows = [line.rstrip().split('\t') for line in handle if '\tCDS\t' in line]
    rows.sort(key=lambda row: int(row[3]), reverse=strand == '-')
    assert [row[7] for row in rows] == ['0', '2'], Path(str(cds['output_path']) + '.cds-normalisation.json').read_text()
    from validate_cds_gff_mapping import validate_single_species
    result = validate_single_species(dict(index=1, species_prefix='Test_species', cds_file=cds['output_path'],
        gff_file=gff['output_path'], genome_file=task['genome_path'], strict=True), 10)
    assert result['ok'], result.get('error')
    assert {name: Path(task[name]).read_bytes() for name in original} == original


@pytest.mark.parametrize('mutation', ['genomic_mismatch', 'shared_parent', 'missing_phase', 'exception', 'internal_stop'])
def test_partial_phase_normalisation_retains_unsupported_source_models(tmp_path, mutation):
    load_module()
    from format_species_annotation.cds_normalisation import iter_normalised_cds_records
    sequence = 'AATAACCCAA' if mutation == 'internal_stop' else 'AACC'
    task = normalisation_bundle(tmp_path, [{'blocks': [('ATGA', 0), (sequence, 1)]}])
    text = task['gff_path'].read_text()
    if mutation == 'genomic_mismatch':
        task['cds_path'].write_text('>t1\nATGCCACC\n')
    elif mutation == 'shared_parent':
        text = text.replace('Parent=t1', 'Parent=t1,t2')
    elif mutation == 'missing_phase':
        text = text.replace('\t1\tParent=t1', '\t.\tParent=t1')
    elif mutation == 'exception':
        text = text.replace('ID=t1;Parent=g1', 'ID=t1;Parent=g1;exception=ribosomal_slippage')
    task['gff_path'].write_text(text)
    state = {}
    list(iter_normalised_cds_records(task, state))
    assert state['gff_updates'] == {}


@pytest.mark.parametrize('kind,suffix,filename', [
    ('cds', '-translated-nucleotides.fa.gz', '.cds.fa.gz'),
    ('transcriptome', '-SOAPdenovo-Trans-assembly.fa.gz', '.transcripts.fa.gz'),
])
def test_onekp_explicit_transcriptome_selects_assembly_for_native_orfs(monkeypatch, kind, suffix, filename):
    load_module()
    from format_species_providers import onekp
    from format_species_providers.dispatch import resolve_provider_specific_download_urls_from_id
    sample = dict(code='YKSS', species='Osmunda regalis', is_mislabeled=False,
        assembly_url='https://ftp.cngb.org/pub/SciRAID/onekp/assemblies/YKSS-Osmunda_regalis/')
    monkeypatch.setattr(onekp, 'sample_index', lambda *_: {'YKSS': sample})
    result = resolve_provider_specific_download_urls_from_id('onekp', 'YKSS', 'Osmunda_regalis', 30, {},
        row={'cds_input_type': kind})
    assert result['cds_url'] == sample['assembly_url'] + 'YKSS' + suffix
    assert result['cds_filename'] == 'YKSS' + filename
    assert result['onekp_sample_code'] == 'YKSS'
    assert result['onekp_reported_species'] == 'Osmunda regalis'
    with pytest.raises(ValueError, match='cds_input_type'):
        onekp.resolve_onekp_download_urls_from_id('YKSS', 'Osmunda_regalis', 30, {}, row={'cds_input_type': 'protein'})
    sample['is_mislabeled'] = True
    with pytest.raises(ValueError, match='mislabeled'):
        onekp.resolve_onekp_download_urls_from_id('YKSS', 'Osmunda_regalis', 30, {}, row={'cds_input_type': kind})


def test_native_planner_accepts_onekp_sample_id_for_transcriptome(tmp_path):
    import json
    import subprocess
    import sys

    manifest = tmp_path / 'manifest.tsv'
    manifest.write_text('provider\tid\tspecies_key\tcds_input_type\n'
                        'onekp\tYKSS\tOsmunda_regalis\ttranscriptome\n')
    plan = tmp_path / 'plan.json'
    script = Path(__file__).resolve().parents[1] / 'support/plan_input_generation_tasks.py'
    result = subprocess.run([sys.executable, str(script), '--provider', 'all',
        '--download-manifest', str(manifest), '--download-dir', str(tmp_path / 'downloads'),
        '--stage-downloads', '--strict', '--outfile', str(plan)], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    task = json.loads(plan.read_text())['tasks'][0]
    assert task['cds_input_type'] == 'transcriptome'
    assert task['manifest_row']['id'] == 'YKSS'


@pytest.mark.parametrize('provider,code', [('local', 'YKSS'), ('unknown', 'YKSS'), ('onekp', ''), ('onekp', 'YKSS-extra')])
def test_transcript_source_still_requires_a_supported_explicit_identity(provider, code):
    load_module()
    from transcriptome_orfs import validate_manifest_row
    with pytest.raises(ValueError, match='explicit sequence URL/path'):
        validate_manifest_row(dict(provider=provider, id=code, cds_input_type='transcriptome'))


@pytest.mark.parametrize('recompute', ['changed_kind', 'unvalidated'])
def test_new_plan_recomputes_changed_input_kind_and_reuses_unchanged_species(tmp_path, recompute):
    import csv
    import hashlib
    import json
    import subprocess

    from workflow.tests.test_gg_input_generation_end_to_end import (
        CORE_PATH,
        REPO_ROOT,
        _core_env,
        _install_fake_toolchain,
        _write_direct_species_fixture,
        _write_minimal_ete_taxonomy_db,
        _write_runtime_busco_dataset,
    )
    raw = _write_direct_species_fixture(tmp_path)
    toolchain = _install_fake_toolchain(tmp_path)
    donor = tmp_path / 'donor'
    target = tmp_path / 'target'
    source_plan = donor / 'output/input_generation/tmp/task_plan.json'
    logs = {}
    for workspace in [donor, target]:
        _write_minimal_ete_taxonomy_db(workspace)
        _write_runtime_busco_dataset(workspace)
        manifest = workspace / 'manifest.tsv'
        with manifest.open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=['provider', 'id', 'species_key', 'cds_url', 'cds_input_type'], delimiter='\t')
            writer.writeheader()
            for species in ['Arabidopsis_thaliana', 'Oryza_sativa']:
                changed = recompute == 'changed_kind' and workspace == target and species == 'Arabidopsis_thaliana'
                writer.writerow(dict(provider='local', id=species, species_key=species,
                    cds_url=(raw / species / (species + '.cds.fa')).as_uri(),
                    cds_input_type='transcriptome' if changed else 'cds'))
        for mode, index in [('array_prepare', None), ('array_worker', 1), ('array_worker', 2), ('array_finalize', None)]:
            env = _core_env(workspace, None, toolchain, mode, index, require_cds=True)
            env.update(provider='all', download_manifest=str(manifest), download_dir=str(workspace / 'downloads'),
                       strict='1', overwrite='0', run_species_taxonomy='1')
            if workspace == target:
                env.update(resume_from_task_plan=str(source_plan),
                    resume_from_task_plan_sha256=hashlib.sha256(source_plan.read_bytes()).hexdigest(),
                    resume_from_input_generation_root=str(source_plan.parents[1]))
            result = subprocess.run(['bash', str(CORE_PATH)], cwd=REPO_ROOT, env=env,
                                    capture_output=True, text=True, timeout=180)
            assert result.returncode == 0, result.stdout + result.stderr
            logs[(workspace, mode, index)] = result.stdout
        if workspace == donor and recompute == 'unvalidated':
            # An interrupted/failed validator leaves formatting without proof
            # of validated inputs or worker completion. Reformat this species.
            (donor / 'output/input_generation/tmp/stage_checkpoints/validate.Arabidopsis_thaliana.json').unlink()
            Path(str(source_plan) + '.completed/1.json').unlink()
    assert 'Reused verified input formatting' not in logs[(target, 'array_worker', 1)]
    assert 'Reused verified input formatting' in logs[(target, 'array_worker', 2)]
    root = target / 'output/input_generation'
    receipt = json.loads((root / 'tmp/task_plan.json.completed/1.json').read_text())
    assert any(path.endswith('.orf-extraction.json') for path in receipt['files']) == (recompute == 'changed_kind')
    old = donor / 'output/input_generation/species_cds_busco_full/Oryza_sativa.busco.full.tsv'
    new = root / 'species_cds_busco_full/Oryza_sativa.busco.full.tsv'
    assert old.read_bytes() == new.read_bytes()
