"""Real genomic admission from frozen evidence, without another prediction run."""
import copy
import importlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / 'support'
sys.path.insert(0, str(SUPPORT))
rescue = importlib.import_module('rescue_gene_models')
replay = importlib.import_module('rescue_admission_replay')
store = importlib.import_module('rescue_model_store')


@pytest.fixture
def completed(tmp_path, request):
    """Actual prepared anchors, indexed genomic QC and ordinary stage receipts."""
    root, source = tmp_path / 'previous', tmp_path / 'producer'
    root.mkdir()
    tools = rescue.identities()
    bindings = {**replay.MODULES, **{name: '' for name in replay.FUNCTIONS}}
    for name in bindings:
        path = source / 'workflow/support' / name
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(SUPPORT / name, path)
    sources, files = {}, {}
    names = ['Target_species', 'Relative_species']
    for name in names:
        directory = tmp_path / name
        directory.mkdir()
        coding = 'ATGCCCTAA' if name == names[0] else 'ATGAAATAA'
        (directory / 'genome.fa').write_text('>owned\n' + coding + '\n>missing\nATGAAATAA\n>disrupted\nATGTAATAA\n')
        (directory / 'cds.fa').write_text(f'>{name}_g\n{coding}\n')
        (directory / 'busco.txt').write_text('# frozen synthetic pre-rescue summary\n')
        (directory / 'annotation.gff3').write_text(
            '##gff-version 3\nowned\ttest\tgene\t1\t9\t.\t+\t.\tID=g\n'
            'owned\ttest\tmRNA\t1\t9\t.\t+\t.\tID=t;Parent=g\n'
            'owned\ttest\tCDS\t1\t9\t.\t+\t0\tParent=t\n')
        sources[name] = {'species': name, 'fasta': str(directory / 'cds.fa'), 'gff': str(directory / 'annotation.gff3'),
                         'genome': str(directory / 'genome.fa'), 'mode': 'cds', 'genetic_code': 1,
                         'busco': str(directory / 'busco.txt'), 'quality': {'complete_pct': 90.0},
                         'feature': 'gene', 'attribute': 'ID'}
        files.update({str(p): rescue.digest(p) for p in directory.iterdir()})
    defaults = {a.dest: a.default for a in rescue.parser()._subparsers._group_actions[0].choices['plan']._actions}
    params = {key: defaults[key] for key in rescue.PARAMETERS}
    plan = {'output': str(root), 'species': names, 'common_references': [names[1]],
            'donors': {names[0]: [names[1]], names[1]: [names[0]]},
            'nearest_references': {names[0]: [names[1]], names[1]: [names[0]]},
            'synteny_jobs': [{'id': 'comparison_000001', 'index': 1, 'a': names[0], 'b': names[1], 'kind': 'pair'}],
            'request': {'schema': rescue.SCHEMA, 'sources': sources, 'files': files, 'parameters': params,
                        'tools': tools, 'species_profiles': {}, 'gemoma_jar': None,
                        'output_storage': {'format': 'compact', 'retain_search_inputs': False}}}
    history_case = getattr(request, 'param', None)
    if history_case in {'separate-history', 'separate-code', 'duplicate-history'}:
        ancestor = tmp_path / 'separate-annotation-history'
        ancestor.mkdir()
        ancestor_plan = {**copy.deepcopy(plan), 'output': str(ancestor)}
        if history_case == 'separate-code':
            ancestor_plan['request']['sources'][names[0]]['genetic_code'] = 4
        rescue.atomic_json(ancestor / 'plan.json', ancestor_plan)
        sha = rescue.digest(ancestor / 'plan.json')
        receipts = {}
        for name in names:
            directory = ancestor / 'rescued' / name
            rows = [{'model_id': 'historical_primary', 'seqid': 'missing', 'strand': '+', 'cds': [[0, 9, 0]],
                     'sequence': 'ATGAAATAA', 'problems': [], 'status': 'accepted',
                     'evidence': {'donor': names[1], 'query': names[1] + '_g'}}] if name == names[0] else []
            store.write_model_store(directory, rows, codec='gzip')
            rescue.atomic_json(directory / 'candidates.json', [])
            files = {str(p.relative_to(directory)): rescue.digest(p) for p in directory.rglob('*') if p.is_file()}
            rescue.atomic_json(directory / 'receipt.json', {'key': {'plan': sha, 'species': name}, 'files': files})
            receipts[name] = rescue.digest(directory / 'receipt.json')
        plan['request']['integration_reclassification'] = {'source_root': str(ancestor), 'source_plan_sha256': sha,
                                                           'source_worker_receipts': receipts}
        if history_case == 'duplicate-history':
            plan['request']['prediction_cache'] = rescue.frozen_prediction_cache_key(ancestor)
    for name, raw in {'selection.nwk': '(Target_species,Relative_species);', 'quality.tsv': 'quality\n',
                      'references.nwk': '(Target_species,Relative_species);', 'references.tsv': 'references\n'}.items():
        (root / name).write_text(raw)
    rescue.atomic_json(root / 'plan.json', plan)
    for name in names:
        rescue.prepared(root, plan, name)
    job = plan['synteny_jobs'][0]
    rescue.stage(root, Path('synteny') / job['id'], rescue.comparison_key(root, job),
                 lambda tmp: rescue.atomic_json(tmp / 'blocks.json', []))
    donor = next(rescue.fasta_records(root / 'prepared' / names[1] / 'genes.pep'))[0]
    evidence = {'id': 'transfer', 'target': names[0], 'donor': names[1], 'query': donor,
                'genome_only': True, 'placement': 'unanchored'}
    intact = {'query': 'transfer', 'id': 'raw', 'seqid': 'missing', 'strand': '+', 'cds': [[0, 9, 0]],
              'frameshift': False, 'coverage': 1.0, 'identity': 1.0, 'sequence': 'ATGAAATAA',
              'genetic_code': 1, 'problems': ['competing_new_models'], 'status': 'unresolved',
              'evidence': evidence, 'partial_evidence': {'partial': False},
              'query_start': 0, 'query_end': 2, 'query_length': 2,
              'unknown': {'number': -0.0, 'unicode': '植物'}}
    rejected = copy.deepcopy(intact)
    rejected.update(query='strict-rejection', seqid='disrupted', sequence='ATGTAATAA', problems=['internal_stop'])
    rows = [intact, rejected]
    for name in names:
        records = rows if name == names[0] else []
        def build(tmp, records=records, name=name):
            store.write_model_store(tmp, records, shard_bytes=400, codec='gzip')
            rescue.atomic_json(tmp / 'candidates.json', [])
            rescue.atomic_json(tmp / 'search_inputs.json', {'prepared_receipts': {}, 'schema': 1})
            rescue.write_tsv(tmp / 'audit.tsv', ['candidate', 'donor_gene', 'status', 'reasons', 'coverage', 'identity', 'model_id'],
                             [(m['query'], m['evidence']['query'], m['status'], ','.join(m['problems']), 1, 1, '') for m in records]
                             + [('no-result', donor, 'unresolved', 'no_alignment', '', '', '')])
            rescue.write_tsv(tmp / 'quality_flags.tsv', ['candidate', 'model_id', 'status', 'start_codon',
                             'donor_n_terminus_aligned', 'donor_c_terminus_aligned', 'donor_species', 'flags',
                             'donor_aligned_query_fraction', 'donor_internal_unaligned_query_fraction'],
                             [(m['query'], '', m['status'], '', '', '', '', '', '', '') for m in records])
            (tmp / 'unknown.txt').write_text('complete original publication\n')
        rescue.stage(root, Path('rescued') / name, rescue.rescue_key(root, plan, name), build)
    return root, source, plan, rows


def test_standard_replay_reuses_evidence_rechecks_admission_and_retains_rejections(completed, tmp_path, monkeypatch):
    previous, producer, old, rows = completed
    before = {str(p): p.read_bytes() for p in previous.rglob('*') if p.is_file()}
    output = tmp_path / 'new-generation'
    plan = replay.build_plan(previous, producer, output)
    assert plan['request']['parameters'] == old['request']['parameters']
    replay.prepare(output, plan, 2)
    def forbidden(*args, **kwargs):
        raise AssertionError('Prediction/alignment must not run in admission replay')
    monkeypatch.setattr(rescue, 'search_intervals', forbidden)
    monkeypatch.setattr(rescue, 'run', forbidden)
    native = importlib.import_module('annevo_gene_models')
    monkeypatch.setattr(native, 'nominate', forbidden)
    validator, calls = rescue.validate_model, []
    def fresh(*args):
        calls.append(args[0]['query'])
        return validator(*args)
    monkeypatch.setattr(rescue, 'validate_model', fresh)
    result = rescue.rescue(output, plan, 'Target_species', 2)
    current = list(store.iter_models(result))
    assert current[0]['status'] == 'accepted' and current[0]['problems'] == []
    assert current[1] == rows[1]
    assert current[0]['unknown'] == rows[0]['unknown']
    assert calls == ['transfer']
    assert list(store.iter_accepted_models(result)) == [current[0]]
    audit = (result / 'audit.tsv').read_text().splitlines()
    assert audit[-1].startswith('no-result\t') and len(audit) == 4
    assert (result / 'unknown.txt').read_text() == 'complete original publication\n'
    proof = json.loads((result / 'admission_replay.json').read_text())
    assert proof['full_stream_verified'] and proof['admitted_cds_freshly_validated']
    assert proof['counts'] == {'raw_records': 2, 'admission_rechecked': 1, 'unchanged_strict_rejections': 1}
    assert {str(p): p.read_bytes() for p in previous.rglob('*') if p.is_file()} == before


@pytest.mark.parametrize('fault', ['validator', 'raw-module', 'runtime', 'unbound-source', 'native-confirmation'])
def test_replay_abstains_when_reuse_contract_cannot_be_proven(completed, tmp_path, fault):
    previous, source, plan, _ = completed
    if fault == 'runtime':
        plan['request']['tools']['biopython'] = 'different'
    else:
        filename = {'validator': 'rescue_gene_models.py', 'raw-module': 'rescue_raw_validation.py',
                    'unbound-source': 'rescue_model_quality.py', 'native-confirmation': 'annevo_gene_models.py'}[fault]
        path = source / 'workflow/support' / filename
        old = path.read_text()
        if fault in {'validator', 'native-confirmation'}:
            name = 'validate_model' if fault == 'validator' else 'confirm_transferred_structure'
            old = old.replace('def ' + name + '(', 'def altered_' + name + '(')
        else:
            old += '\n# changed producer bytes\n'
        path.write_text(old)
        field = {**replay.MODULES, 'rescue_gene_models.py': 'implementation',
                 'annevo_gene_models.py': 'annevo_candidate_implementation'}[filename]
        if fault != 'unbound-source':
            plan['request']['tools'][field] = rescue.digest(path)
    with pytest.raises((ValueError, KeyError)):
        replay.prove_contract(plan, source)
    assert not (tmp_path / 'new').exists()


@pytest.mark.parametrize('fault', ['parameters', 'donors', 'profiles', 'raw-bytes', 'prepared-bytes', 'genome'])
def test_changed_inputs_or_frozen_scope_never_publish(completed, tmp_path, fault):
    previous, source, old, _ = completed
    output = tmp_path / 'new'
    plan = replay.build_plan(previous, source, output)
    if fault == 'parameters':
        plan['request']['parameters']['minimum_coverage'] = .1
    elif fault == 'donors':
        plan['donors']['Target_species'] = []
    elif fault == 'profiles':
        plan['request']['species_profiles'] = {'Target_species': {'minimum_coverage': .1}}
    elif fault == 'raw-bytes':
        path = previous / 'rescued/Target_species/model_store/bodies.sqlite'
        path.write_bytes(path.read_bytes() + b'bad')
    elif fault == 'prepared-bytes':
        path = previous / 'prepared/Target_species/genes.pep'
        path.write_text(path.read_text() + 'bad')
    else:
        path = Path(old['request']['sources']['Target_species']['genome'])
        path.write_text(path.read_text().replace('ATGAAATAA', 'ATGCCCTAA'))
    rescue.atomic_json(output / 'plan.json', plan)
    with pytest.raises((ValueError, OSError)):
        replay.prepare(output, plan, 1)
        rescue.rescue(output, plan, 'Target_species', 1)
    assert not (output / 'rescued/Target_species/receipt.json').exists()


def test_replay_requires_all_workers_and_new_output(completed, tmp_path):
    previous, source, _, _ = completed
    with pytest.raises(ValueError, match='new, disjoint'):
        replay.build_plan(previous, source, previous / 'child')
    (previous / 'rescued/Relative_species/receipt.json').unlink()
    with pytest.raises(ValueError, match='all completed'):
        replay.build_plan(previous, source, tmp_path / 'new')


def test_native_unconfirmed_and_structural_failures_are_not_cleared():
    for reason in ('annevo_unconfirmed_coding_structure', 'oriongeno_unconfirmed_coding_structure',
                   'frameshift', 'internal_stop', 'invalid_phase', 'assembly_gap_within_model_span'):
        row = {'problems': [reason, 'competing_new_models'], 'status': 'unresolved'}
        assert not replay.eligible(row)
        assert replay.reset_admission(row)['problems'] == [reason]


def test_standard_cli_owns_plan_import_worker_and_resume(completed, tmp_path):
    previous, source, _, _ = completed
    output = tmp_path / 'cli'
    command = [sys.executable, str(Path(rescue.__file__))]
    result = subprocess.check_output([*command, 'reclassify-plan', '--previous-output', str(previous),
                                      '--producer-source', str(source), '--output', str(output)], text=True)
    assert json.loads(result)['admission_replay'] is True
    subprocess.check_call([*command, 'reclassify-prepare', '--output', str(output), '--cpus', '2'])
    subprocess.check_call([*command, 'rescue', '--output', str(output), '--cpus', '2', '--task-index', '1'])
    worker = output / 'rescued/Target_species'
    before = {str(p): (p.stat().st_ino, p.read_bytes()) for p in worker.rglob('*') if p.is_file()}
    again = subprocess.check_output([*command, 'rescue', '--output', str(output), '--cpus', '2', '--task-index', '1'], text=True)
    assert 'Reused rescued/Target_species' in again
    assert {str(p): (p.stat().st_ino, p.read_bytes()) for p in worker.rglob('*') if p.is_file()} == before


def test_producer_source_change_during_contract_check_abstains(completed, monkeypatch):
    _, source, plan, _ = completed
    original = replay._definitions
    changed = False
    def replacing(path, names):
        nonlocal changed
        value = original(path, names)
        if not changed:
            changed = True
            target = source / 'workflow/support/rescue_raw_validation.py'
            target.write_text(target.read_text() + '\n# changed during proof\n')
        return value
    monkeypatch.setattr(replay, '_definitions', replacing)
    with pytest.raises(OSError, match='changed'):
        replay.prove_contract(plan, source)


def test_invalid_primary_cannot_survive_as_an_unselected_rejection(completed, tmp_path):
    previous, source, _, rows = completed
    worker = previous / 'rescued/Target_species'
    prior_receipt = json.loads((worker / 'receipt.json').read_text())
    shutil.rmtree(worker / 'model_store')
    records = copy.deepcopy(rows)
    records[1]['status'] = 'accepted'
    store.write_model_store(worker, records, codec='gzip')
    files = {str(p.relative_to(worker)): rescue.digest(p) for p in worker.rglob('*')
             if p.is_file() and p.name != 'receipt.json'}
    rescue.atomic_json(worker / 'receipt.json', {'key': prior_receipt['key'], 'files': files})
    output = tmp_path / 'new'
    plan = replay.build_plan(previous, source, output)
    replay.prepare(output, plan, 1)
    with pytest.raises(ValueError, match='Prior primary|Previous accepted'):
        rescue.rescue(output, plan, 'Target_species', 1)
    assert not (output / 'rescued/Target_species/receipt.json').exists()


@pytest.mark.parametrize('completed', ['separate-history'], indirect=True)
def test_separate_producer_annotation_source_is_bound_and_kept_across_replays(completed, tmp_path):
    previous, source, old, _ = completed
    output = tmp_path / 'new'
    plan = replay.build_plan(previous, source, output)
    parent = old['request']['integration_reclassification']['source_root']
    assert plan['request']['admission_cache']['annotation_history']['root'] == parent
    assert 'integration_reclassification' not in plan['request']
    replay.prepare(output, plan, 1)
    result = rescue.rescue(output, plan, 'Target_species', 1)
    primary, = store.iter_accepted_models(result)
    assert primary['retained_coding_path']['source_directory'] == parent + '/rescued/Target_species'
    assert primary['retained_coding_path']['history_depth'] == 1
    assert primary['retained_coding_path']['locus_model_ids'] == ['historical_primary']
    assert replay._legacy_annotation_history(plan) == plan['request']['admission_cache']['annotation_history']


@pytest.mark.parametrize('completed', ['separate-code'], indirect=True)
def test_separate_annotation_history_cannot_cross_genetic_code(completed, tmp_path):
    previous, source, _, _ = completed
    output = tmp_path / 'new'
    plan = replay.build_plan(previous, source, output)
    replay.prepare(output, plan, 1)
    with pytest.raises(ValueError, match='genome/code differs'):
        rescue.rescue(output, plan, 'Target_species', 1)
    assert not (output / 'rescued/Target_species/receipt.json').exists()


@pytest.mark.parametrize('completed', ['duplicate-history'], indirect=True)
def test_annotation_ancestor_already_in_prediction_history_is_not_read_twice(completed, tmp_path):
    previous, source, _, _ = completed
    output = tmp_path / 'new'
    plan = replay.build_plan(previous, source, output)
    replay.prepare(output, plan, 1)
    worker = rescue.rescue(output, plan, 'Target_species', 1)
    history = json.loads((worker / 'retained_annotation_sources.json').read_text())
    parent = plan['request']['admission_cache']['annotation_history']['root']
    rows = [row for row in history if row['directory'] == parent + '/rescued/Target_species' and row['status'] == 'verified']
    assert len(rows) == 1
