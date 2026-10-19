"""Command-local runtime proofs retain exact membership and mutation fences."""
import copy
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
import annevo_runtime_cache as cache
import input_generation_array_state as state


@pytest.fixture
def runtime_tree(tmp_path):
    root = tmp_path / 'runtime'
    nested = root / 'lib/package'
    nested.mkdir(parents=True)
    (nested / 'module.py').write_text('original runtime\n')
    (root / 'alias').symlink_to('lib/package/module.py')
    return root, cache.inspect_tree(root)


def test_repeated_boundaries_keep_all_member_checks_without_reenumerating(runtime_tree, monkeypatch):
    root, recorded = runtime_tree
    verification = cache.RuntimeVerification()
    verification.verify_tree(root, recorded)
    def no_scan(*a, **kw):
        pytest.fail('an already verified runtime was enumerated again')
    monkeypatch.setattr(Path, 'rglob', no_scan)
    verification.verify_tree(root, copy.deepcopy(recorded))
    verification.check()


@pytest.mark.parametrize('change', ['add-root', 'add-nested', 'delete', 'rewrite', 'replace',
                                   'permission', 'directory-replace', 'directory-link', 'link', 'inventory'])
def test_reuse_rejects_runtime_or_inventory_changes(runtime_tree, change):
    root, recorded = runtime_tree
    verification = cache.RuntimeVerification()
    verification.verify_tree(root, recorded)
    member = root / 'lib/package/module.py'
    if change.startswith('add-'):
        (root / 'extra' if change == 'add-root' else member.with_name('extra')).write_text('new')
    elif change == 'delete':
        member.unlink()
    elif change == 'rewrite':
        before = member.stat()
        member.write_text('modified runtime\n')
        os.utime(member, ns=(before.st_atime_ns, before.st_mtime_ns))
    elif change == 'replace':
        replacement = root.parent / 'replacement'
        replacement.write_bytes(member.read_bytes())
        replacement.replace(member)
    elif change == 'permission':
        member.chmod(0o400)
    elif change.startswith('directory-'):
        old = member.parent.with_name('previous')
        member.parent.rename(old)
        if change == 'directory-link':
            member.parent.symlink_to(old, target_is_directory=True)
        else:
            member.parent.mkdir()
            member.write_bytes((old / member.name).read_bytes())
    elif change == 'link':
        (root / 'alias').unlink()
        (root / 'alias').symlink_to('./lib/package/module.py')
    else:
        recorded['lib/package/module.py']['sha256'] = 'a' * 64
    with pytest.raises((RuntimeError, OSError)):
        verification.verify_tree(root, recorded)
    with pytest.raises(RuntimeError, match='previously failed'):
        verification.check()


def test_first_verification_rejects_incomplete_inventory(runtime_tree):
    root, recorded = runtime_tree
    del recorded['lib/package/module.py']
    with pytest.raises(RuntimeError, match='inventory changed'):
        cache.RuntimeVerification().verify_tree(root, recorded)


def test_change_during_initial_walk_is_not_absorbed(runtime_tree, monkeypatch):
    root, recorded = runtime_tree
    original = cache.signature
    def changed(path):
        result = original(path)
        if path == root / 'lib/package/module.py':
            (path.parent / 'extra').write_text('racing addition')
        return result
    monkeypatch.setattr(cache, 'signature', changed)
    with pytest.raises(RuntimeError, match='changed'):
        cache.RuntimeVerification().verify_tree(root, recorded)


def test_full_byte_hash_is_shared_only_inside_one_verification(tmp_path, monkeypatch):
    model = tmp_path / 'model.pt'
    model.write_bytes(b'model weights')
    reads = []
    original = state.digest
    def digest(path):
        reads.append(str(path))
        return original(path)
    monkeypatch.setattr(state, 'digest', digest)
    first = cache.RuntimeVerification()
    expected = first.file_sha(model)
    assert first.file_sha(model) == expected
    first.check()
    assert reads == [str(model)]
    assert cache.RuntimeVerification().file_sha(model) == expected
    assert reads == [str(model), str(model)]


@pytest.mark.parametrize('change', ['rewrite', 'replace', 'delete'])
def test_changed_model_cannot_refresh_its_proof(tmp_path, change):
    model = tmp_path / 'model.pt'
    model.write_bytes(b'original weights')
    verification = cache.RuntimeVerification()
    verification.file_sha(model)
    if change == 'rewrite':
        model.write_bytes(b'modified weights')
    elif change == 'replace':
        replacement = tmp_path / 'replacement'
        replacement.write_bytes(model.read_bytes())
        replacement.replace(model)
    else:
        model.unlink()
    with pytest.raises(OSError):
        verification.file_sha(model)
    with pytest.raises(RuntimeError, match='previously failed'):
        verification.file_sha(model)


def test_model_changed_during_initial_hash_is_rejected(tmp_path, monkeypatch):
    model = tmp_path / 'model.pt'
    model.write_bytes(b'original weights')
    original = state.digest
    def changed(path):
        result = original(path)
        model.write_bytes(b'modified weights')
        return result
    monkeypatch.setattr(state, 'digest', changed)
    with pytest.raises(OSError, match='changed'):
        cache.RuntimeVerification().file_sha(model)
