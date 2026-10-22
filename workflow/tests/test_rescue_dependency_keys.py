"""Pair-key reuse retains fresh bytes, exact keys and immutable generations."""
import json
import os
import sys
from collections import Counter
from importlib import import_module
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
rescue = import_module('rescue_gene_models')
ReceiptSnapshot = import_module('refinement_receipt_snapshot').ReceiptSnapshot


def fixture(tmp_path):
    previous = import_module('test_rescue_prepared_snapshot')
    root, plan, names = previous.fixture(tmp_path, count=4)
    plan['species'] = names
    plan['donors'] = {name: [other for other in names if other != name] for name in names}
    rescue.atomic_json(root / 'plan.json', plan)
    plan_hash = rescue.digest(root / 'plan.json')
    for name in names:
        directory = root / 'prepared' / name
        receipt = json.loads((directory / 'receipt.json').read_text())
        receipt['key']['plan'] = plan_hash
        rescue.atomic_json(directory / 'receipt.json', receipt)
    for job in plan['synteny_jobs']:
        directory = root / 'synteny' / job['id']
        receipt = json.loads((directory / 'receipt.json').read_text())
        receipt['key'] = rescue.comparison_key(root, job)
        rescue.atomic_json(directory / 'receipt.json', receipt)
    return root, plan, names


def test_same_complete_keys_and_each_prepared_publication_verified_once(tmp_path, monkeypatch):
    root, plan, names = fixture(tmp_path)
    expected = [rescue.rescue_key(root, plan, name) for name in names]
    original = rescue.verified
    calls = Counter()

    def counted(directory, key, **kwargs):
        calls[str(directory)] += 1
        return original(directory, key, **kwargs)

    monkeypatch.setattr(rescue, 'verified', counted)
    with ReceiptSnapshot(counted) as snapshot:
        for _ in range(3):
            assert [rescue.rescue_key(root, plan, name, dependency_snapshot=snapshot) for name in names] == expected
    assert {calls[str(root / 'prepared' / name)] for name in names} == {1}
    assert all(calls[str(root / 'synteny' / job['id'])] == 1 for job in plan['synteny_jobs'])
    first = calls.copy()
    calls.clear()
    with ReceiptSnapshot(counted) as snapshot:
        for name in names:
            rescue.rescue_key(root, plan, name, dependency_snapshot=snapshot)
    assert calls == first  # No content proof crosses a CLI boundary.


def test_comparison_key_without_snapshot_keeps_fresh_verification(tmp_path, monkeypatch):
    root, plan, names = fixture(tmp_path)
    original = rescue.verified
    calls = Counter()

    def counted(directory, key, **kwargs):
        calls[str(directory)] += 1
        return original(directory, key, **kwargs)

    monkeypatch.setattr(rescue, 'verified', counted)
    job = plan['synteny_jobs'][0]
    expected = rescue.comparison_key(root, job)
    assert rescue.comparison_key(root, job) == expected
    assert calls[str(root / 'prepared' / names[0])] == 2


def test_two_preparation_contexts_cannot_be_combined(tmp_path):
    root, plan, _names = fixture(tmp_path)
    with pytest.raises(ValueError, match='one preparation'):
        rescue.comparison_key(root, plan['synteny_jobs'][0], object(), dependency_snapshot=object())


@pytest.mark.parametrize('target', ['member', 'receipt'])
@pytest.mark.parametrize('mutation', ['rewrite', 'replace', 'remove', 'permission', 'symlink'])
def test_changed_preparation_cannot_refresh_inside_snapshot(tmp_path, target, mutation):
    root, plan, names = fixture(tmp_path)
    snapshot = ReceiptSnapshot(rescue.verified)
    rescue.rescue_key(root, plan, names[0], dependency_snapshot=snapshot)
    path = root / 'prepared' / names[0] / ('genes.pep' if target == 'member' else 'receipt.json')
    raw, info = path.read_bytes(), path.stat()
    if mutation == 'rewrite':
        path.write_bytes(bytes([raw[0] ^ 1]) + raw[1:])
        os.utime(path, ns=(info.st_atime_ns, info.st_mtime_ns))
    elif mutation == 'replace':
        alternate = path.with_suffix('.replacement')
        alternate.write_bytes(raw)
        os.utime(alternate, ns=(info.st_atime_ns, info.st_mtime_ns))
        alternate.replace(path)
    elif mutation == 'remove':
        path.unlink()
    elif mutation == 'permission':
        path.chmod(0o600 if info.st_mode & 0o077 else 0o644)
    else:
        alternate = path.with_suffix('.target')
        alternate.write_bytes(raw)
        path.unlink()
        path.symlink_to(alternate)
    with pytest.raises((OSError, ValueError)):
        rescue.rescue_key(root, plan, names[0], dependency_snapshot=snapshot)
    if path.is_symlink():
        path.unlink()
    path.write_bytes(raw)
    with pytest.raises((OSError, ValueError)):
        rescue.rescue_key(root, plan, names[0], dependency_snapshot=snapshot)
    snapshot.close()


def test_initial_corruption_is_never_admitted_as_a_cached_key(tmp_path):
    root, plan, names = fixture(tmp_path)
    (root / 'prepared' / names[0] / 'genes.pep').write_text('corrupted before first proof')
    with ReceiptSnapshot(rescue.verified) as snapshot:
        with pytest.raises(ValueError):
            rescue.rescue_key(root, plan, names[0], dependency_snapshot=snapshot)
