"""Lossless admission updates, immutable source proof and sparse decoding."""
import copy
import json
import os

import pytest

from workflow.support import rescue_model_store as store
from workflow.support.rescue_model_store_patch import patch_model_store
from workflow.tests.test_rescue_model_store import model, published, receipt


@pytest.mark.parametrize('codec', ['gzip', 'auto'])
def test_patch_keeps_raw_order_unknown_fields_partials_and_independent_bytes(tmp_path, codec):
    records = [model(0, 'accepted'), model(1), model(2, 'accepted'), model(3)]
    source, before = published(tmp_path, records, partial=[records[3], records[1], records[3]], codec=codec)
    original = {str(p.relative_to(source)): p.read_bytes() for p in source.rglob('*') if p.is_file()}
    changes = {0: {**records[0], 'status': 'unresolved', 'problems': ['competing_new_models']},
               1: {**records[1], 'status': 'accepted', 'new_admission_metadata': {'unicode': '植物'}}}
    appended = {**model(4, 'accepted'), 'sequence': 'ATGAAATAA', 'cds': [[6, 15, 0]]}
    revisions = [{**model(2), 'support': [{'new': 'proof'}]}]
    output = tmp_path / 'result'
    result = patch_model_store(source, output, changes, appended=[appended], revisions=revisions)
    receipt(output)
    expected = [changes.get(i, m) for i, m in enumerate(records)] + [appended]
    assert json.dumps(list(store.iter_models(output)), ensure_ascii=False) == json.dumps(expected, ensure_ascii=False)
    assert list(store.iter_accepted_models(output)) == [changes[1], records[2], appended]
    assert list(store.iter_partial_models(output)) == [records[3], changes[1], records[3]]
    assert list(store.iter_revision_models(output)) == revisions
    assert result['counts']['models'] == before['counts']['models'] + 1
    assert result['counts']['partial'] == before['counts']['partial']
    assert {str(p.relative_to(source)): p.read_bytes() for p in source.rglob('*') if p.is_file()} == original
    for name in before['files']:
        assert not os.path.samefile(source / 'model_store' / name, output / 'model_store' / name)


@pytest.mark.parametrize('field,value', [('sequence', 'ATGCCCTAA'), ('partial_evidence', {'partial': True}),
                                        ('cds', [[4, 37, 0]])])
def test_patch_rejects_changed_predictor_dna_and_partial_membership(tmp_path, field, value):
    row = model()
    source, _ = published(tmp_path, [row])
    with pytest.raises(ValueError, match='predictor/DNA'):
        patch_model_store(source, tmp_path / 'output', {0: {**row, field: value}})
    assert not (tmp_path / 'output/model_store').exists()


@pytest.mark.parametrize('fault', ['body', 'shard', 'partial', 'overlay', 'receipt'])
def test_patch_rejects_corrupted_source_members(tmp_path, fault):
    source, manifest = published(tmp_path, [model(0, 'accepted'), model(1)])
    frozen = store.frozen_model_store_key(source)
    name = {'body': 'model_store/bodies.sqlite', 'shard': 'model_store/' + manifest['shards'][0]['path'],
            'partial': 'model_store/partial.sqlite', 'overlay': 'model_store/accepted.sqlite',
            'receipt': 'receipt.json'}[fault]
    path = source / name
    path.write_bytes(path.read_bytes() + b'corrupt')
    with pytest.raises((ValueError, OSError, store.sqlite3.DatabaseError)):
        patch_model_store(source, tmp_path / 'output', {}, frozen_key=frozen)
    assert not (tmp_path / 'output/model_store').exists()


def test_context_selection_keeps_original_ordinals_without_decoding_rejections(tmp_path, monkeypatch):
    rows = [model(i) for i in range(10)]
    for i, m in enumerate(rows):
        m['problems'] = [] if i in {2, 7} else ['internal_stop']
    source, _ = published(tmp_path, rows)
    unpack, calls = store._unpack, []
    def record(envelope, *args, **kwargs):
        calls.append(envelope['ordinal'])
        return unpack(envelope, *args, **kwargs)
    monkeypatch.setattr(store, '_unpack', record)
    assert list(store.iter_model_items(source, context_filter=lambda c: not c['problems'])) == [(2, rows[2]), (7, rows[7])]
    assert calls == [2, 7]


def test_context_selection_finishes_generation_fence_after_last_match(tmp_path):
    source, _ = published(tmp_path, [model(0), model(1)])
    iterator = store.iter_model_items(source, context_filter=lambda c: c['query'] == 'query0')
    assert next(iterator)[0] == 0
    (source / 'receipt.json').write_text((source / 'receipt.json').read_text() + '\n')
    with pytest.raises(ValueError, match='generation changed'):
        list(iterator)


def test_patch_source_replacement_during_copy_is_not_published(tmp_path, monkeypatch):
    source, _ = published(tmp_path, [model(0), model(1)])
    from workflow.support import rescue_model_store_patch as patch
    original, changed = patch._copy, False
    def replace(*args):
        nonlocal changed
        original(*args)
        if not changed:
            changed = True
            path = source / 'model_store/bodies.sqlite'
            replacement = source / 'new-body'
            replacement.write_bytes(path.read_bytes())
            replacement.replace(path)
    monkeypatch.setattr(patch, '_copy', replace)
    with pytest.raises(ValueError, match='changed'):
        patch_model_store(source, tmp_path / 'output', {})
    assert not (tmp_path / 'output/model_store').exists()


def test_patch_rejects_overlap_and_bad_ordinals_before_write(tmp_path):
    source, _ = published(tmp_path, [model()])
    with pytest.raises(ValueError, match='overlaps'):
        patch_model_store(source, source / 'child', {})
    for ordinal in (True, -1, 1):
        with pytest.raises(ValueError, match='ordinal'):
            patch_model_store(source, tmp_path / 'output', {ordinal: copy.deepcopy(model())})


@pytest.mark.parametrize('alter_predictor', [False, True])
def test_legacy_metric_shared_bodies_are_normalized_only_with_exact_predictor_values(tmp_path, monkeypatch, alter_predictor):
    rows = [model(0), model(1, 'accepted')]
    # This is the retained public v1 layout used by the archived producer,
    # before donor metrics moved into per-candidate contexts. No old code runs.
    with monkeypatch.context() as prior:
        prior.setattr(store, 'BODY_FIELDS', store.BODY_FIELDS | {'coverage', 'identity', 'query_span_coverage',
                                                              'query_start', 'query_end', 'query_length'})
        source, _ = published(tmp_path, rows)
    assert list(store.iter_models(source)) == rows
    changed = {**rows[0], 'status': 'accepted'}
    if alter_predictor:
        changed['coverage'] = .2
        with pytest.raises(ValueError, match='predictor/DNA'):
            patch_model_store(source, tmp_path / 'output', {0: changed})
        return
    output = tmp_path / 'output'
    patch_model_store(source, output, {0: changed})
    receipt(output)
    assert list(store.iter_models(output)) == [changed, rows[1]]
    assert list(store.iter_accepted_models(output)) == [changed, rows[1]]
    assert list(store.iter_partial_models(output)) == [rows[1]]


@pytest.mark.parametrize('field,value', [('identity', .5), ('query', 'other'), ('coverage', .2)])
def test_current_layout_also_rejects_changed_predictor_context(tmp_path, field, value):
    row = model()
    source, _ = published(tmp_path, [row])
    with pytest.raises(ValueError, match='predictor/DNA'):
        patch_model_store(source, tmp_path / 'output', {0: {**row, field: value}})


@pytest.mark.parametrize('codec', ['gzip', 'auto'])
@pytest.mark.parametrize('explicit_order', [False, True])
def test_new_partial_candidates_append_without_changing_old_partial_references(tmp_path, codec, explicit_order):
    rows = [model(0), model(1), model(2), model(3)]
    partials = [rows[3], rows[1], rows[3]] if explicit_order else [rows[1], rows[3]]
    source, before = published(tmp_path, rows, partial=partials, codec=codec)
    original = {str(p.relative_to(source)): p.read_bytes() for p in source.rglob('*') if p.is_file()}
    appended = [model(4), model(5), model(6), model(7)]
    with pytest.raises(ValueError, match='partial reference stream'):
        patch_model_store(source, tmp_path/'unqualified', {}, appended=appended)
    output = tmp_path/'output'
    result = patch_model_store(source, output, {}, appended=appended, append_partial_evidence=True)
    receipt(output)
    assert list(store.iter_models(output)) == rows + appended
    assert list(store.iter_partial_models(output)) == partials + [appended[1], appended[3]]
    assert result['counts']['partial'] == before['counts']['partial'] + 2
    assert result['partial_ordered'] == before['partial_ordered']
    assert {str(p.relative_to(source)): p.read_bytes() for p in source.rglob('*') if p.is_file()} == original
