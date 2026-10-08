"""Both rescue publications retain scientific gates and frozen consumer inputs."""
import copy
import json
import shutil
import sqlite3
from pathlib import Path

import pytest

from workflow.support import rescue_model_evidence as evidence
from workflow.support import rescue_model_store as model_store
from workflow.support.input_generation_array_state import digest
from workflow.support.rescue_prediction_cache import frozen_prediction_cache_key, verify_prediction_cache
from workflow.tests.test_rescue_additional_candidates import PARAMS, cache_fixture
from workflow.tests.test_rescue_model_evidence import fixture as evidence_fixture
from workflow.tests.test_rescue_refinement_bridge import (
    Genome,
    classify,
    refinement,
    rescue_model,
    resolve,
)
from workflow.tests.test_rescue_refinement_bridge import (
    fixture as revision_fixture,
)


def receipt(directory, key):
    files = {str(path.relative_to(directory)): digest(path) for path in directory.rglob('*')
             if path.is_file() and path.name != 'receipt.json'}
    (directory / 'receipt.json').write_text(json.dumps({'key': key, 'files': files}))


def publish(directory, models, *, revisions=(), sharded=True, key=None):
    directory.mkdir(parents=True, exist_ok=True)
    if sharded:
        # Deliberately replace a synthetic generation for frozen-key tests.
        if (directory / 'model_store').exists():
            shutil.rmtree(directory / 'model_store')
        for name in ('models.json', 'partial_models.json', 'revision_candidates.json'):
            (directory / name).unlink(missing_ok=True)
        model_store.write_model_store(directory, models, revisions=list(revisions), shard_bytes=1024)
    else:
        for name, values in (('models.json', models), ('partial_models.json', []),
                             ('revision_candidates.json', list(revisions))):
            (directory / name).write_text(json.dumps(values))
    receipt(directory, key or {})


@pytest.mark.parametrize('sharded', [False, True])
def test_revision_publication_preserves_normal_trust_and_hard_qc(tmp_path, sharded):
    catalog, models, edges = revision_fixture()
    raw = rescue_model(models)
    directory = tmp_path / 'worker'
    publish(directory, [], revisions=[raw], sharded=sharded)
    if sharded:
        request = {'revision_model_stores': {'Species_target': {
            'directory': str(directory), 'key': model_store.frozen_model_store_key(directory, kind='revision')}}}
    else:
        request = {'revision_candidates': {'Species_target': str(directory / 'revision_candidates.json')}}
    assert refinement.revision_worker(request, 'Species_target') == directory
    imported, proposals = refinement.import_rescue_revisions(
        refinement.revision_models(request, 'Species_target'), catalog, edges, refinement.DEFAULTS,
        Genome(raw['sequence']), resolve)
    assert not proposals
    accepted, = classify(catalog, imported, edges)
    assert accepted['status'] == 'accepted'
    assert accepted['change_type'] == 'model_revision'
    assert accepted['donors'] == ['Species_donor1', 'Species_donor2']
    raw['problems'].append('internal_stop')
    publish(directory, [], revisions=[raw], sharded=sharded)
    if sharded:
        request['revision_model_stores']['Species_target']['key'] = model_store.frozen_model_store_key(
            directory, kind='revision')
    imported, _ = refinement.import_rescue_revisions(
        refinement.revision_models(request, 'Species_target'), catalog, edges, refinement.DEFAULTS,
        Genome(raw['sequence']), resolve)
    rejected, = classify(catalog, imported, edges)
    assert rejected['status'] == 'proposal' and 'internal_stop' in rejected['problems']


def test_revision_generation_change_is_not_silently_imported(tmp_path):
    _catalog, models, _edges = revision_fixture()
    directory = tmp_path / 'worker'
    publish(directory, [], revisions=[rescue_model(models)])
    request = {'revision_model_stores': {'Species_target': {
        'directory': str(directory), 'key': model_store.frozen_model_store_key(directory, kind='revision')}}}
    publish(directory, [], revisions=[])
    with pytest.raises((ValueError, OSError)):
        list(refinement.revision_models(request, 'Species_target'))


@pytest.mark.parametrize('sharded', [False, True])
def test_prediction_reuse_keeps_query_context_and_discards_old_acceptance(tmp_path, sharded):
    old, new, plan, region = cache_fixture(tmp_path)
    directory = old / 'rescued/T'
    rows = json.loads((directory / 'models.json').read_text())
    rows.append({**copy.deepcopy(rows[0]), 'status': 'rejected', 'model_id': 'old_other_decision',
                 'problems': ['internal_stop'], 'sequence': 'OLD_UNTRUSTED_DNA'})
    old_key = json.loads((directory / 'receipt.json').read_text())['key']
    publish(directory, rows, sharded=sharded, key=old_key)
    frozen = frozen_prediction_cache_key(old, ['T'])
    cache = verify_prediction_cache(old, new, plan, 'T', [region], PARAMS, frozen=frozen)
    reused = list(cache.iter_models())
    assert len(reused) == 2
    assert all(row['query'] == 'region1' and row['evidence'] == region for row in reused)
    assert all('status' not in row and 'model_id' not in row and 'sequence' not in row
               and 'problems' not in row and 'support' not in row for row in reused)


@pytest.mark.parametrize('fault', ['member', 'manifest', 'receipt', 'omitted_shard', 'accepted_scope', 'boolean_schema'])
def test_prediction_store_freeze_detects_content_and_membership_changes(tmp_path, fault):
    old, new, plan, region = cache_fixture(tmp_path)
    directory = old / 'rescued/T'
    rows = json.loads((directory / 'models.json').read_text())
    old_key = json.loads((directory / 'receipt.json').read_text())['key']
    publish(directory, rows, key=old_key)
    frozen = frozen_prediction_cache_key(old, ['T'])
    store_key = frozen['species']['T']['model_store']
    if fault == 'manifest':
        member = 'model_store/manifest.json'
        with (directory / member).open('ab') as handle:
            handle.write(b'\n')
    elif fault == 'receipt':
        with (directory / 'receipt.json').open('ab') as handle:
            handle.write(b'\n')
    elif fault == 'accepted_scope':
        store_key['kind'] = 'accepted'
    elif fault == 'boolean_schema':
        store_key['schema'] = True
    else:
        member = next(name for name in store_key['files'] if name.endswith(('.gz', '.zst')))
        if fault == 'member':
            with (directory / member).open('ab') as handle:
                handle.write(b'changed')
        else:
            # A forged narrow frozen key must not enable unchecked raw reads.
            store_key['files'].pop(member)
            frozen['species']['T']['files'].pop(member)
    with pytest.raises((ValueError, OSError)):
        verify_prediction_cache(old, new, plan, 'T', [region], PARAMS, frozen=frozen)


def test_accepted_evidence_reads_no_rejected_body_or_shard(tmp_path, monkeypatch):
    args, _genome = evidence_fixture(tmp_path)
    directory = args.rescue_output / 'rescued' / args.species
    models = json.loads((directory / 'models.json').read_text())
    models.extend({**copy.deepcopy(models[0]), 'status': 'rejected', 'model_id': f'rejected_{i}',
                   'query': f'query_{i}'} for i in range(30))
    key = json.loads((directory / 'receipt.json').read_text())['key']
    publish(directory, models, key=key)
    raw = model_store.frozen_model_store_key(directory, kind='models')
    accepted = model_store.frozen_model_store_key(directory, kind='accepted')
    forbidden = {str((directory / member).resolve()) for member in set(raw['files']) - set(accepted['files'])}
    open_path, connect = Path.open, sqlite3.connect
    def guarded_open(path, *a, **k):
        assert str(path.resolve()) not in forbidden, 'accepted-only consumer read rejected model storage'
        return open_path(path, *a, **k)
    def guarded_connect(database, *a, **k):
        assert not any(path in str(database) for path in forbidden), 'accepted-only consumer opened raw body database'
        return connect(database, *a, **k)
    monkeypatch.setattr(Path, 'open', guarded_open)
    monkeypatch.setattr(sqlite3, 'connect', guarded_connect)
    assert evidence.audit(args)['accepted_models'] == 1
    accepted_row, = json.loads((args.output / 'evidence.json').read_text())
    assert accepted_row['model_id'] == models[0]['model_id']
    assert accepted_row['decision'] == 'unchanged' and accepted_row['functional_status'] == 'not_established'


def test_repeat_join_binds_all_consumed_accepted_store_members(tmp_path):
    import gene_model_refinement_busco as busco
    args, _genome = evidence_fixture(tmp_path)
    directory = args.rescue_output / 'rescued' / args.species
    models = json.loads((directory / 'models.json').read_text())
    key = json.loads((directory / 'receipt.json').read_text())['key']
    publish(directory, models, key=key)
    args.output = tmp_path / 'evidence' / args.species
    evidence.audit(args)
    store_key = model_store.frozen_model_store_key(directory, kind='accepted')
    changes = {'rescue_reference_selection': {'plan_sha256': digest(args.rescue_output / 'plan.json')},
               'evidence': {args.species: {
                   'rescue_receipt_sha256': digest(directory / 'receipt.json'),
                   'rescue_model_members_sha256': store_key['files'],
                   'rescued_loci_support': {'gene': {'source_model_id': models[0]['model_id']}}}},
               'species': {args.species: {'refinement_status': 'analysed', 'prior_rescued_loci': 1}}}
    busco.collect_rescue_repeat_evidence(changes, args.output.parent)
    assert changes['species'][args.species]['rescue_repeat_groups']['not_assessed'] == 1
    member = next(iter(store_key['files']))
    changes['evidence'][args.species]['rescue_model_members_sha256'][member] = '0' * 64
    with pytest.raises(ValueError, match='different rescue models'):
        busco.collect_rescue_repeat_evidence(changes, args.output.parent)
