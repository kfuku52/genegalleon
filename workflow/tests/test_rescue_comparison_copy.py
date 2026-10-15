"""A copied cache keeps exact bytes and rejects changes until the copy ends."""
import hashlib
import json
import os
from pathlib import Path

import pytest

from workflow.support import rescue_gene_models as rescue


@pytest.fixture(params=[1, 4])
def comparison(tmp_path, request):
    source, destination = tmp_path / 'source', tmp_path / 'destination'
    source.mkdir()
    destination.mkdir()
    files = {'alignment': b'ACGT' * (600 * 1024), 'nested/empty': b'',
             'blocks.json': b'[[["gene_a","gene_b"]]]\n'}
    for name, raw in files.items():
        path = source / name
        path.parent.mkdir(exist_ok=True)
        path.write_bytes(raw)
    key = {'job': {'a': 'A', 'b': 'B', 'kind': 'pair'}}
    receipt = {'key': key, 'files': {name: hashlib.sha256(raw).hexdigest() for name, raw in files.items()}}
    (source / 'receipt.json').write_text(json.dumps(receipt))
    return source, destination, key, receipt, request.param


def test_comparison_copy_preserves_all_bytes_and_independent_files(comparison):
    source, destination, key, receipt, workers = comparison
    rescue.copy_verified_comparison(source, destination, key, hash_workers=workers)
    assert {str(p.relative_to(destination)): rescue.digest(p) for p in destination.rglob('*') if p.is_file()} == receipt['files']
    for name in receipt['files']:
        assert not os.path.samefile(source / name, destination / name)
    (destination / 'alignment').write_text('a later writer')
    assert rescue.digest(source / 'alignment') == receipt['files']['alignment']


@pytest.mark.parametrize('fault', ['key', 'content', 'missing', 'unsafe', 'receipt_member', 'empty'])
def test_comparison_copy_rejects_invalid_source(comparison, fault):
    source, destination, key, receipt, workers = comparison
    if fault == 'key':
        key = {'another': 'job'}
    elif fault == 'content':
        (source / 'alignment').write_text('corrupt')
    elif fault == 'missing':
        (source / 'alignment').unlink()
    else:
        receipt['files'] = ({'../outside': '0' * 64} if fault == 'unsafe' else
                            {'receipt.json': '0' * 64} if fault == 'receipt_member' else {})
        (source / 'receipt.json').write_text(json.dumps(receipt))
    with pytest.raises((ValueError, OSError)):
        rescue.copy_verified_comparison(source, destination, key, hash_workers=workers)


@pytest.mark.parametrize('fault', ['content', 'receipt', 'replacement', 'redirect'])
def test_comparison_copy_fences_previously_copied_sources(comparison, monkeypatch, fault):
    source, destination, key, _receipt, workers = comparison
    original, changed = rescue.digest, False

    def mutate_after_destination_hash(path):
        nonlocal changed
        value = original(path)
        if Path(path) == destination / 'alignment' and not changed:
            changed = True
            target = source / 'alignment'
            if fault == 'content':
                target.write_bytes(target.read_bytes()[:-1] + b'X')
            elif fault == 'receipt':
                (source / 'receipt.json').write_text((source / 'receipt.json').read_text() + '\n')
            else:
                replacement = source.parent / 'replacement'
                replacement.write_bytes(target.read_bytes())
                info = target.stat()
                os.utime(replacement, ns=(info.st_atime_ns, info.st_mtime_ns))
                if fault == 'replacement':
                    replacement.replace(target)
                else:
                    target.unlink()
                    target.symlink_to(replacement)
        return value

    monkeypatch.setattr(rescue, 'digest', mutate_after_destination_hash)
    with pytest.raises((ValueError, OSError), match='changed'):
        rescue.copy_verified_comparison(source, destination, key, hash_workers=workers)
    assert changed


def test_comparison_copy_rehashes_written_destination(comparison, monkeypatch):
    source, destination, key, _receipt, workers = comparison
    original = rescue.digest

    def corrupt_before_hash(path):
        if Path(path) == destination / 'alignment':
            Path(path).write_text('corrupt destination')
        return original(path)

    monkeypatch.setattr(rescue, 'digest', corrupt_before_hash)
    with pytest.raises(ValueError, match='changed during copying'):
        rescue.copy_verified_comparison(source, destination, key, hash_workers=workers)


def test_comparison_copy_rejects_nonregular_source_without_blocking(comparison):
    source, destination, key, _receipt, workers = comparison
    (source / 'alignment').unlink()
    os.mkfifo(source / 'alignment')
    with pytest.raises(ValueError, match='regular file'):
        rescue.copy_verified_comparison(source, destination, key, hash_workers=workers)


@pytest.mark.parametrize('alias', ['source', 'nested', 'symlink', 'hardlink'])
def test_comparison_copy_cannot_overwrite_source_aliases(comparison, alias):
    source, destination, key, receipt, workers = comparison
    if alias == 'source':
        destination = source
    elif alias == 'nested':
        destination = source / 'copy'
    elif alias == 'symlink':
        (destination / 'alignment').symlink_to(source / 'alignment')
    else:
        os.link(source / 'alignment', destination / 'alignment')
    with pytest.raises((ValueError, OSError)):
        rescue.copy_verified_comparison(source, destination, key, hash_workers=workers)
    assert rescue.digest(source / 'alignment') == receipt['files']['alignment']


def test_failed_copy_does_not_replace_previous_complete_stage(comparison, tmp_path, monkeypatch):
    source, _destination, key, _receipt, workers = comparison
    rescue.stage(tmp_path, 'published', key, lambda p: (p / 'old').write_text('previous'))
    original = rescue.digest

    def change_source(path):
        result = original(path)
        if Path(path).name == 'alignment' and Path(path).parent != source:
            (source / 'alignment').write_text('changed after copying')
        return result

    monkeypatch.setattr(rescue, 'digest', change_source)
    with pytest.raises((ValueError, OSError), match='changed'):
        rescue.stage(tmp_path, 'published', {'new': 'generation'},
                     lambda p: rescue.copy_verified_comparison(source, p, key, hash_workers=workers))
    assert (tmp_path / 'published/old').read_text() == 'previous'
    assert rescue.verified(tmp_path / 'published', key)
