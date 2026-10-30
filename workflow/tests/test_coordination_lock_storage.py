"""Offline metadata retirement preserves science and every live/unknown owner."""
import hashlib
import json
from pathlib import Path

import pytest

from workflow.support import coordination_lock_storage as storage
from workflow.support import gene_family_output_store as store
from workflow.support.shared_namespace_lock import namespace_lock


def populated(tmp_path, count=3):
    root=tmp_path/'output'
    root.mkdir()
    (root/'.gg_store').mkdir()
    science=root/'scientific.tsv'
    science.write_bytes(b'unchanged original scientific data\n')
    markers=[]
    for number in range(count):
        digest=hashlib.sha256(f'OG{number:07d}'.encode()).hexdigest()
        paths=[root/'.gg_run_locks'/f'.family.{digest}.lock.guard',
               root/'.gg_store_locks/family-locks-v2'/digest[:2]/f'{digest}.lock']
        for path in paths:
            with namespace_lock(path,exclusive=True):
                pass
            markers.append(path)
    journal=tmp_path/'journal'
    journal.mkdir(mode=0o700)
    return root,journal,markers


def test_full_idle_retirement_and_original_journal_replay(tmp_path):
    root,journal,markers=populated(tmp_path)
    science=(root/'scientific.tsv').read_bytes()
    expected=storage.snapshot(root)
    assert expected['groups']==6 and expected['entries']==18
    result=storage.retire(root,journal,expected,writers_stopped=True)
    assert result['verified'] and result['retired_entries']==18
    assert storage.snapshot(root)['groups']==0
    assert all(not p.exists() and not Path(str(p)+'.namespace-v1').exists()for p in markers)
    assert (root/'scientific.tsv').read_bytes()==science
    assert storage.retire(root,journal,expected,writers_stopped=True)==result
    # Original runtimes can recreate the same per-family paths; no layout change.
    with namespace_lock(markers[0],exclusive=True):
        assert markers[0].exists()


def test_missing_offline_writer_fence_never_starts(tmp_path):
    root,journal,markers=populated(tmp_path)
    expected=storage.snapshot(root)
    with pytest.raises(storage.CoordinationLockError,match='fenced'):
        storage.retire(root,journal,expected)
    assert all(p.exists()for p in markers)
    assert not list(journal.iterdir())


def test_primary_run_lease_is_not_retired_even_when_old(tmp_path):
    root,journal,markers=populated(tmp_path)
    expected=storage.snapshot(root)
    primary=markers[0].with_name(markers[0].name[1:-len('.guard')])
    primary.write_text(json.dumps({'pid':123,'owner':'preserve'}))
    with pytest.raises(storage.CoordinationLockError,match='Primary producer lease'):
        storage.retire(root,journal,expected,writers_stopped=True)
    assert all(p.exists()for p in markers)
    assert primary.read_text()==json.dumps({'pid':123,'owner':'preserve'})


def test_active_namespace_is_preserved(tmp_path):
    root,journal,markers=populated(tmp_path)
    expected=storage.snapshot(root)
    with namespace_lock(markers[0],exclusive=True):
        with pytest.raises(storage.CoordinationLockError,match='occupied'):
            storage.retire(root,journal,expected,writers_stopped=True)
    assert all(p.exists()for p in markers)


def test_family_gate_excludes_even_a_waiting_family_operation(tmp_path):
    root,journal,markers=populated(tmp_path)
    expected=storage.snapshot(root)
    gate=store.family_gate_path(root/'.gg_store')
    with namespace_lock(gate,exclusive=False):
        with pytest.raises(storage.CoordinationLockError,match='Family operations'):
            storage.retire(root,journal,expected,writers_stopped=True)
    assert all(p.exists()for p in markers)


def test_replaced_marker_is_not_resealed(tmp_path):
    root,journal,markers=populated(tmp_path)
    expected=storage.snapshot(root)
    original=markers[0].read_bytes()
    replacement=markers[0].with_name('replacement')
    replacement.write_bytes(original)
    replacement.chmod(0o600)
    replacement.replace(markers[0])
    with pytest.raises(storage.CoordinationLockError,match='snapshot changed'):
        storage.retire(root,journal,expected,writers_stopped=True)
    assert all(p.exists()for p in markers)


def test_adjacent_unknown_files_and_current_global_gates_are_kept(tmp_path):
    root,journal,markers=populated(tmp_path)
    unknown=root/'.gg_run_locks'/'unclassified-result.json'
    unknown.write_bytes(b'user-owned adjacent bytes')
    storage.retire(root,journal,storage.snapshot(root),writers_stopped=True)
    assert unknown.read_bytes()==b'user-owned adjacent bytes'
    assert store.family_gate_path(root/'.gg_store').exists()
    assert (root/'.gg_store').is_dir()
    assert all(not p.exists()for p in markers)


def test_partial_same_journal_resumes_without_another_snapshot(tmp_path,monkeypatch):
    root,journal,markers=populated(tmp_path)
    expected=storage.snapshot(root)
    target=Path(str(markers[0])+'.namespace-v1')/'readers'
    original=Path.rmdir
    def interrupted(path):
        result=original(path)
        if path==target:
            raise RuntimeError('simulated interruption')
        return result
    with monkeypatch.context()as patch:
        patch.setattr(Path,'rmdir',interrupted)
        with pytest.raises(RuntimeError,match='interruption'):
            storage.retire(root,journal,expected,writers_stopped=True)
    manifest=journal/'coordination-guards-original.jsonl.gz'
    digest=hashlib.sha256(manifest.read_bytes()).hexdigest()
    result=storage.retire(root,journal,expected,writers_stopped=True)
    assert result['verified']
    assert hashlib.sha256(manifest.read_bytes()).hexdigest()==digest
    assert all(not p.exists()for p in markers)


def test_unrecorded_new_guard_stops_before_removing_any_original(tmp_path):
    root,journal,markers=populated(tmp_path)
    expected=storage.snapshot(root)
    fresh=root/'.gg_run_locks'/('.family.'+'f'*64+'.lock.guard')
    with namespace_lock(fresh,exclusive=True):
        pass
    with pytest.raises(storage.CoordinationLockError,match='snapshot changed'):
        storage.retire(root,journal,expected,writers_stopped=True)
    assert all(p.exists()for p in markers)
    assert fresh.exists()
