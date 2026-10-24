"""Native cohort handoff keeps IDs, source evidence and independent QC."""
import errno
import gzip
import json
import sys
from importlib import import_module
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / 'support'
sys.path.insert(0, str(SUPPORT))
cohort = import_module('input_generation_cohort')
formatter = import_module('format_species_inputs')
state = import_module('input_generation_array_state')
atomic_json, digest, verify_receipt = state.atomic_json, state.digest, state.verify_receipt


def fixture(tmp_path, *, clones=False, transcriptome=False):
    selection = []
    for number, species in enumerate(['Plant_one', 'Plant_two'], 1):
        workspace = tmp_path / ('donor' + str(number))
        root = workspace / 'output/input_generation'
        directory = root / 'species_cds'
        directory.mkdir(parents=True)
        raw = workspace / 'source.fa'
        raw.write_text('>g1\nATGAAATAA\n>g2\nATGCCCTAA\n')
        task = dict(provider='local', species_key=species, species_prefix=species, cds_path=raw,
                    gff_path=None, genome_path=None, gbff_path=None, genetic_code=1,
                    gene_grouping_mode='strict', gff_repair_mode='safe')
        if transcriptome:
            task['cds_input_type'] = 'transcriptome'
        formatted = formatter.format_cds(task, directory, False, False)
        output = str(formatted['output_path'])
        row = dict(provider='local', species_key=species, species_prefix=species,
                   gene_grouping_mode='strict', gff_repair_mode='safe', cds_input_path=str(raw),
                   cds_output_path=output, gff_output_path='', genome_output_path='')
        files = [raw, Path(output)]
        if transcriptome:
            row['cds_input_type'] = 'transcriptome'
            audit = Path(output + '.orf-extraction.json')
            value = json.loads(audit.read_text())
            value['contract']['implementation'] = 'a' * 64
            value['contract'].pop('gene_grouping')
            audit.write_text(json.dumps(value))
            files.append(audit)
        if clones and number == 1:
            audit = workspace / 'cloned_models/audit.json'
            audit.parent.mkdir()
            received = audit.parent / 'received.fa'
            received.write_text('>DaIspS1\nATGAAATAA\n>DaIspS2\nATGCCCTAA\n')
            audit.write_text(json.dumps(dict(mode='sequence', outputs={'received_cloned_cds': received.name},
                decisions=[dict(clone_id='DaIspS1', formatted_output_id=species + '_g1', analysis_included=True,
                                confirmed_gene_copy_count=None, copy_number_status='unresolved'),
                           dict(clone_id='DaIspS2', formatted_output_id=species + '_g2', analysis_included=True,
                                confirmed_gene_copy_count=None, copy_number_status='unresolved')],
                removed_CDS_ids=['old_fragment'], original_sources_modified=False)))
            row.update(cloned_model_audit_path=str(audit), cloned_cds_mode='sequence')
            files.extend([audit, received])
        short = root / 'species_cds_busco_short' / (species + '.busco.short.txt')
        short.parent.mkdir()
        short.write_text('# BUSCO version is: 6.1.0\n# The lineage dataset is: eukaryota_odb12\n'
                         '# BUSCO was run in mode: transcriptome\nC:100.0%[S:100.0%,D:0.0%],F:0.0%,M:0.0%,n:1\n')
        files.append(short)
        full = root / 'species_cds_busco_full' / (species + '.busco.full.tsv')
        full.parent.mkdir()
        full.write_text('M1\tComplete\t' + species + '_g1\n')
        files.append(full)
        plan = root / 'tmp/task_plan.json'
        atomic_json(plan, dict(version=2, task_count=1, tasks=[dict(species_prefix=species)]))
        shard = root / 'tmp/species_summary_shards/1.tsv'
        cohort.write_table(shard, [row])
        files.append(shard)
        receipt = Path(str(plan) + '.completed/1.json')
        atomic_json(receipt, dict(task_index=1, species_prefix=species, plan_sha256=digest(plan),
                                 files={str(p): digest(p) for p in files}))
        assert verify_receipt(plan, 1)
        selection.append(dict(species=species, source_workspace=str(workspace), plan_path=str(plan),
                              plan_sha256=digest(plan), cds_path=output, gff_path='', genome_path='',
                              cloned_model_audit_path=row.get('cloned_model_audit_path', '')))
    table = tmp_path / 'selection.tsv'
    cohort.write_table(table, selection)
    plan = tmp_path / 'cohort_plan.json'
    atomic_json(plan, cohort.freeze(table, 2), immutable=True)
    return plan, tmp_path / 'view', table


def test_prepared_view_keeps_exact_cds_ids_and_cloned_history(tmp_path):
    plan, workspace, _ = fixture(tmp_path, clones=True)
    doc = cohort.load_cohort(plan)
    before = {Path(p): Path(p).read_bytes() for t in doc['tasks'] for p in t['files']}
    cohort.prepare(workspace, plan)
    cohort.prepare(workspace, plan)
    cohort.verify_view(workspace, plan)
    assert all(p.read_bytes() == data for p, data in before.items())
    for task in doc['tasks']:
        source = Path(task['summary']['cds_output_path'])
        view = workspace / 'input/species_cds' / source.name
        assert not view.is_symlink() and view.read_bytes() == source.read_bytes()
        with gzip.open(view, 'rt') as stream:
            assert [line.strip() for line in stream if line.startswith('>')] == ['>' + task['species'] + '_g1', '>' + task['species'] + '_g2']
    audit = Path(doc['tasks'][0]['summary']['cloned_model_audit_path'])
    assert (workspace / 'output/input_generation/cloned_model_history/Plant_one.json').read_bytes() == audit.read_bytes()
    assert json.loads((workspace / 'cohort_view.json').read_text())['sequence_identifiers_rewritten'] is False
    for index in (1, 2):
        assert cohort.qc_one(workspace, plan, index)['passed'] is True
    result = cohort.finalize(workspace, plan)
    assert result['passed'] and result['species'] == 2 and result['cloned_audits'] == 1


@pytest.mark.parametrize('problem', ['source_change', 'receipt_change', 'plan_change', 'missing_history'])
def test_mutated_native_evidence_cannot_prepare(tmp_path, problem):
    plan, workspace, _ = fixture(tmp_path, clones=True)
    doc = cohort.load_cohort(plan)
    task = doc['tasks'][0]
    path = {'source_change': task['summary']['cds_output_path'], 'receipt_change': task['receipt_path'],
            'plan_change': task['plan_path'], 'missing_history': task['summary']['cloned_model_audit_path']}[problem]
    Path(path).write_bytes(b'changed')
    with pytest.raises((ValueError, OSError)):
        cohort.prepare(workspace, plan)
    assert not (workspace / 'cohort_view.json').exists()


def test_finalizer_requires_every_exact_qc_result(tmp_path):
    plan, workspace, _ = fixture(tmp_path)
    cohort.prepare(workspace, plan)
    cohort.qc_one(workspace, plan, 1)
    with pytest.raises(FileNotFoundError):
        cohort.finalize(workspace, plan)
    cohort.qc_one(workspace, plan, 2)
    path = workspace / 'qc/2.json'
    result = json.loads(path.read_text())
    result['species'] = 'Wrong_species'
    path.write_text(json.dumps(result))
    with pytest.raises(ValueError, match='mismatched'):
        cohort.finalize(workspace, plan)


def test_source_selection_cannot_alias_or_choose_unverified_outputs(tmp_path):
    plan, workspace, table = fixture(tmp_path)
    rows = cohort.read_table(table)
    rows[1]['species'] = rows[0]['species']
    cohort.write_table(table, rows)
    with pytest.raises(ValueError, match='duplicate'):
        cohort.freeze(table, 2)
    with pytest.raises(ValueError, match='overlaps'):
        cohort.prepare(Path(cohort.load_cohort(plan)['tasks'][0]['source_workspace']), plan)


def test_cloned_snapshot_tampering_is_detected(tmp_path):
    plan, workspace, _ = fixture(tmp_path, clones=True)
    cohort.prepare(workspace, plan)
    (workspace / 'output/input_generation/cloned_model_history/Plant_one.json').write_text('{}')
    with pytest.raises(ValueError, match='hash mismatch'):
        cohort.verify_view(workspace, plan)


def test_unselected_species_file_is_not_silently_adopted(tmp_path):
    plan, workspace, _ = fixture(tmp_path)
    cohort.prepare(workspace, plan)
    (workspace / 'input/species_cds/Foreign_species.fa').write_text('>Foreign_species_g1\nATG\n')
    with pytest.raises(ValueError, match='Unexpected cohort species'):
        cohort.verify_view(workspace, plan)


def test_qc_cannot_publish_after_a_native_source_changes_during_validation(tmp_path, monkeypatch):
    plan, workspace, _ = fixture(tmp_path)
    cohort.prepare(workspace, plan)
    validator = import_module('validate_longest_cds_selection')
    original = validator.validate_single_species
    def change_input(task, limit):
        result = original(task, limit)
        raw = Path(task['summary_row']['cds_input_path'])
        raw.write_text(raw.read_text() + '>new\nATGTAA\n')
        return result
    monkeypatch.setattr(validator, 'validate_single_species', change_input)
    with pytest.raises(OSError, match='changed'):
        cohort.qc_one(workspace, plan, 1)
    assert not (workspace / 'qc/1.json').exists()


def test_verification_hashes_hardlinked_views_once_per_fresh_boundary(tmp_path, monkeypatch):
    plan, workspace, _ = fixture(tmp_path)
    cohort.prepare(workspace, plan)
    calls = []
    original = state.digest
    def counted(path):
        calls.append(Path(path).stat().st_ino)
        return original(path)
    monkeypatch.setattr(state, 'digest', counted)
    cohort.verify_view(workspace, plan)
    output = Path(cohort.load_cohort(plan)['tasks'][0]['summary']['cds_output_path'])
    assert calls.count(output.stat().st_ino) == 1
    output.write_bytes(b'changed')
    with pytest.raises(ValueError, match='hash mismatch'):
        cohort.verify_view(workspace, plan)


def test_raw_source_relocation_preserves_original_receipt_and_rejects_different_bytes(tmp_path):
    plan, workspace, selection = fixture(tmp_path)
    original = cohort.load_cohort(plan)['tasks'][0]
    raw = Path(original['summary']['cds_input_path'])
    replacement = tmp_path / 'recovered.fa'
    replacement.write_bytes(raw.read_bytes())
    raw.unlink()
    aliases = {str(raw): str(replacement)}
    doc = cohort.freeze(selection, 2, source_aliases=aliases)
    relocated = doc['tasks'][0]
    assert relocated['receipt_sha256'] == original['receipt_sha256']
    assert relocated['source_summary'] == original['source_summary']
    assert relocated['qc_summary']['cds_input_path'] == str(replacement)
    assert str(raw) not in relocated['files']
    relocated_plan = tmp_path / 'relocated_plan.json'
    atomic_json(relocated_plan, doc)
    cohort.prepare(workspace, relocated_plan)
    assert cohort.qc_one(workspace, relocated_plan, 1)['passed']
    replacement.write_text('>different\nATGTAA\n')
    with pytest.raises(ValueError, match='content differs'):
        cohort.freeze(selection, 2, source_aliases=aliases)
    with pytest.raises(ValueError, match='hash mismatch'):
        cohort.verify_view(workspace, relocated_plan)


def test_relocation_cannot_replace_unverified_evidence_or_symlinks(tmp_path):
    plan, _, selection = fixture(tmp_path)
    task = cohort.load_cohort(plan)['tasks'][0]
    raw = task['summary']['cds_input_path']
    with pytest.raises(ValueError, match='Unused/unverified'):
        cohort.freeze(selection, 2, source_aliases={task['receipt_path']: raw})
    link = tmp_path / 'link.fa'
    link.symlink_to(raw)
    with pytest.raises(ValueError, match='Invalid raw source relocation'):
        cohort.freeze(selection, 2, source_aliases={raw: str(link)})


def test_cohort_rederives_historical_orfs_and_preserves_the_original_audit(tmp_path):
    plan, workspace, _ = fixture(tmp_path, transcriptome=True)
    task = cohort.load_cohort(plan)['tasks'][0]
    audit = Path(task['summary']['cds_output_path'] + '.orf-extraction.json')
    before = audit.read_bytes()
    cohort.prepare(workspace, plan)
    result = cohort.qc_one(workspace, plan, 1)
    proof = result['ownership']['stats']['transcript_orf_revalidation']
    assert result['passed'] and proof['exact_sequences_rederived']
    assert proof['original_implementation'] == 'a' * 64
    assert proof['gene_grouping'] == 'historical-header-groups-v1'
    assert not proof['audit_rewritten'] and not proof['biological_gene_uniqueness_established']
    assert audit.read_bytes() == before


def test_prepare_hashes_each_physical_source_once_after_link_creation(tmp_path, monkeypatch):
    plan, workspace, _ = fixture(tmp_path)
    task = cohort.load_cohort(plan)['tasks'][0]
    sources = [Path(task['summary'][key]) for key in ['cds_input_path', 'cds_output_path']]
    calls = []
    original = state.digest
    def counted(path):
        calls.append(Path(path).stat().st_ino)
        return original(path)
    monkeypatch.setattr(state, 'digest', counted)
    cohort.prepare(workspace, plan)
    assert all(calls.count(source.stat().st_ino) == 1 for source in sources)
    assert Path(task['summary']['cds_output_path']).samefile(workspace / 'input/species_cds' / sources[1].name)


def test_failed_cross_filesystem_copy_cannot_publish_a_verified_view(tmp_path, monkeypatch):
    plan, workspace, _ = fixture(tmp_path)
    def cross_device(*args, **kwargs):
        raise OSError(errno.EXDEV, 'different filesystem')
    def corrupt_copy(source, destination, _buffer):
        destination.write(source.read() + b'corrupted')
    monkeypatch.setattr(cohort.os, 'link', cross_device)
    monkeypatch.setattr(cohort.shutil, 'copyfileobj', corrupt_copy)
    with pytest.raises(ValueError, match='hash mismatch'):
        cohort.prepare(workspace, plan)
    assert not (workspace / 'cohort_view.json').exists()


def test_partial_view_with_an_extra_species_cannot_be_published(tmp_path):
    plan, workspace, _ = fixture(tmp_path)
    foreign = workspace / 'input/species_cds/Foreign_species.fa'
    foreign.parent.mkdir(parents=True)
    foreign.write_text('>foreign\nATGTAA\n')
    with pytest.raises(ValueError, match='Unexpected cohort species'):
        cohort.prepare(workspace, plan)
    assert not (workspace / 'cohort_view.json').exists()


def test_cohort_plan_mutation_during_staging_prevents_publication(tmp_path, monkeypatch):
    plan, workspace, _ = fixture(tmp_path)
    original = cohort.os.link
    def change_plan(*args, **kwargs):
        plan.write_bytes(plan.read_bytes() + b' ')
        return original(*args, **kwargs)
    monkeypatch.setattr(cohort.os, 'link', change_plan)
    with pytest.raises(OSError, match='changed'):
        cohort.prepare(workspace, plan)
    assert not (workspace / 'cohort_view.json').exists()


def test_optional_qc_telemetry_keeps_native_results_and_failed_content_checks(tmp_path, monkeypatch):
    plan, workspace, _ = fixture(tmp_path, clones=True)
    cohort.prepare(workspace, plan)
    metrics = tmp_path / 'performance'
    monkeypatch.setenv('GG_PERFORMANCE_DIR', str(metrics))
    result = cohort.qc_one(workspace, plan, 1)
    assert result['passed'] and result['ownership']['ok']
    assert json.loads((workspace / 'qc/1.json').read_text()) == result
    task = cohort.load_cohort(plan)['tasks'][1]
    Path(task['summary']['cds_input_path']).write_text('>changed\nATGTAA\n')
    with pytest.raises(ValueError, match='hash mismatch'):
        cohort.qc_one(workspace, plan, 2)
    assert not (workspace / 'qc/2.json').exists()
    rows = [json.loads(line) for path in metrics.glob('*.jsonl') for line in path.read_text().splitlines()]
    content = [row for row in rows if row['phase'] == 'cohort.qc.content']
    selection = [row for row in rows if row['phase'] == 'cohort.qc.cds_selection']
    assert [row['status'] for row in content] == ['ok', 'failed']
    assert len(selection) == 1 and selection[0]['status'] == 'ok'
    assert all(row['seconds'] >= 0 and row['process_peak_rss_bytes'] > 0 for row in content + selection)
