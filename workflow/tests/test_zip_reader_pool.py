"""ZIP reader lifetime, freshness and compaction integrity contracts."""

import contextlib
import gc
import os
import threading
import weakref
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from workflow.support import gene_family_output_store as store_module


def write_zip(path, payload=b'old'):
    with zipfile.ZipFile(path, 'w') as archive:
        archive.writestr('entry', payload)


def test_alternating_zip_readers_reuse_and_reset_close_every_reader(tmp_path, monkeypatch):
    paths = [tmp_path / f'{number}.zip' for number in range(2)]
    for path in paths:
        write_zip(path)
    store = store_module.GeneFamilyOutputStore(tmp_path)
    original = store_module.zipfile.ZipFile
    opened = []

    def observed(*args, **kwargs):
        archive = original(*args, **kwargs)
        opened.append(archive)
        return archive

    monkeypatch.setattr(store_module.zipfile, 'ZipFile', observed)
    for path in paths * 4:
        with store._cached_zip_reader(path) as archive:
            assert archive.read('entry') == b'old'
    assert len(opened) == 2
    assert all(archive.fp is not None for archive in opened)
    store._reset_cache()
    assert all(archive.fp is None for archive in opened)


@pytest.mark.parametrize('atomic', [False, True])
@pytest.mark.parametrize('advance_epoch', [False, True])
def test_alternating_zip_readers_detect_same_epoch_replacement(tmp_path, atomic, advance_epoch):
    paths = [tmp_path / f'{number}.zip' for number in range(2)]
    for path in paths:
        write_zip(path)
    store = store_module.GeneFamilyOutputStore(tmp_path)
    opened = []
    for path in paths:
        with store._cached_zip_reader(path) as archive:
            assert archive.read('entry') == b'old'
            opened.append(archive)
    before = paths[0].stat()
    if atomic:
        replacement = tmp_path / 'replacement.zip'
        write_zip(replacement, b'new')
        os.replace(replacement, paths[0])
    else:
        write_zip(paths[0], b'new')
    os.utime(paths[0], ns=(before.st_atime_ns, before.st_mtime_ns))
    if advance_epoch:
        store._cache_epoch = store._read_index_epoch()
        store.archive_root.mkdir()
        (store.archive_root / store_module.INDEX_EPOCH_FILE).write_text('changed\n')
        store._refresh_if_index_changed()
        assert all(archive.fp is None for archive in opened)
    with store._cached_zip_reader(paths[0]) as archive:
        assert archive.read('entry') == b'new'
        assert archive is not opened[0]
    assert opened[0].fp is None
    store._reset_cache()


def test_corrupt_replacement_is_not_hidden_by_idle_cached_reader(tmp_path):
    first, other = tmp_path / 'first.zip', tmp_path / 'other.zip'
    for path in (first, other):
        write_zip(path)
    store = store_module.GeneFamilyOutputStore(tmp_path)
    for path in (first, other):
        with store._cached_zip_reader(path) as archive:
            assert archive.read('entry') == b'old'
    first.write_bytes(b'corrupt same-epoch archive')
    with pytest.raises(zipfile.BadZipFile):
        with store._cached_zip_reader(first):
            pass
    store._reset_cache()


def test_active_reader_survives_replacement_and_full_pool_overflow(tmp_path):
    paths = [tmp_path / f'{number}.zip' for number in range(9)]
    for path in paths:
        write_zip(path)
    store = store_module.GeneFamilyOutputStore(tmp_path)
    with contextlib.ExitStack() as active:
        readers = [active.enter_context(store._cached_zip_reader(path)) for path in paths[:8]]
        with store._cached_zip_reader(paths[8]) as overflow:
            assert overflow.read('entry') == b'old'
        assert overflow.fp is None
        assert all(reader.fp is not None for reader in readers)
        replacement = tmp_path / 'replacement.zip'
        write_zip(replacement, b'new')
        os.replace(replacement, paths[0])
        with store._cached_zip_reader(paths[0]) as current:
            assert current.read('entry') == b'new'
            assert readers[0].read('entry') == b'old'
        assert current.fp is None
    store._reset_cache()
    assert all(reader.fp is None for reader in readers)


def test_thread_pools_do_not_share_readers_and_exception_exit_releases_users(tmp_path):
    path = tmp_path / 'source.zip'
    write_zip(path)
    store = store_module.GeneFamilyOutputStore(tmp_path)
    barrier = threading.Barrier(2)

    def read():
        with store._cached_zip_reader(path) as archive:
            barrier.wait(timeout=10)
            assert archive.read('entry') == b'old'
            return archive

    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(read) for _ in range(2)]
        readers = [future.result(timeout=10) for future in futures]
        assert readers[0] is not readers[1]
        store._reset_cache()
        assert all(reader.fp is None for reader in readers)
    with pytest.raises(RuntimeError, match='consumer failed'):
        with store._cached_zip_reader(path):
            raise RuntimeError('consumer failed')
    readers = next(iter(store._zip_reader_cache.values()))
    assert readers[path][1] == 0
    store._reset_cache()


@pytest.mark.parametrize('stale', [False, True])
def test_failed_idle_close_remains_registered_for_cleanup(tmp_path, stale):
    paths = [tmp_path / f'{number}.zip' for number in range(9)]
    for path in paths:
        write_zip(path)
    store = store_module.GeneFamilyOutputStore(tmp_path)
    opened = []
    for path in paths[:8]:
        with store._cached_zip_reader(path) as archive:
            opened.append(archive)
    original_close = opened[0].close
    attempts = 0

    def fail_once():
        nonlocal attempts
        attempts += 1
        if attempts == 1:
            raise OSError('injected reader close failure')
        original_close()

    opened[0].close = fail_once
    if stale:
        write_zip(paths[0], b'new')
    with pytest.raises(OSError, match='injected reader close failure'):
        with store._cached_zip_reader(paths[0] if stale else paths[8]):
            pass
    readers = next(iter(store._zip_reader_cache.values()))
    assert readers[paths[0]][0] is opened[0]
    store._reset_cache()
    assert all(archive.fp is None for archive in opened)


def test_snapshot_readers_close_at_fence_independently_of_normal_pool(tmp_path):
    path = tmp_path / 'source.zip'
    write_zip(path)
    store = store_module.GeneFamilyOutputStore(tmp_path)
    with store._cached_zip_reader(path) as normal:
        assert normal.read('entry') == b'old'
    with store.read_snapshot():
        with store._snapshot_zip_reader(path) as snapshot:
            assert snapshot is not normal
            assert snapshot.read('entry') == b'old'
    assert snapshot.fp is None
    assert normal.fp is not None
    assert not store._snapshot_zip_readers
    store._reset_cache()
    assert normal.fp is None


def test_compaction_drops_evicted_reader_objects_and_preserves_payloads(tmp_path, monkeypatch):
    source_root = tmp_path / 'sources'
    source_root.mkdir()
    artifacts = []
    expected = {}
    shards = store_module.MAX_OPEN_COMPACTION_SOURCES + 1
    for shard in range(shards):
        path = source_root / f'{shard}.zip'
        with zipfile.ZipFile(path, 'w') as archive:
            for number in range(2):
                name = f'stat_branch/OG{number * shards + shard:07d}_stat.branch.tsv'
                payload = name.encode()
                archive.writestr(name, payload)
                info = archive.getinfo(name)
                expected[name] = payload
                artifacts.append(store_module.Artifact(logical_path=name, subdir='stat_branch',
                    name=name.split('/', 1)[1], generation=1, zip_path=path, member_name=name,
                    size=info.file_size, crc=info.CRC, family_id=store_module.orthogroup_id_from_name(name)))
    original = store_module.zipfile.ZipFile
    references = []
    max_alive = 0

    def observed(path, *args, **kwargs):
        nonlocal max_alive
        archive = original(path, *args, **kwargs)
        if Path(path).parent == source_root:
            references.append(weakref.ref(archive))
            max_alive = max(max_alive, sum(reference() is not None for reference in references))
        return archive

    monkeypatch.setattr(store_module.zipfile, 'ZipFile', observed)
    destination, compacted = store_module._compact_artifact_chunk(tmp_path / 'out', 'stat_branch',
        sorted(artifacts, key=lambda artifact: artifact.logical_path), 'orthogroup', 2)
    gc.collect()
    assert max_alive <= store_module.MAX_OPEN_COMPACTION_SOURCES + 1
    assert all(reference() is None for reference in references)
    assert len(compacted) == len(expected)
    with zipfile.ZipFile(destination) as archive:
        assert {name: archive.read(name) for name in expected} == expected


@pytest.mark.parametrize('failure', ['digest', 'crc', 'close'])
def test_failed_compaction_closes_sources_and_preserves_previous_archive(tmp_path, monkeypatch, failure):
    source_root = tmp_path / 'sources'
    source_root.mkdir()
    artifacts = []
    for number in range(2):
        name = f'stat_branch/OG{number:07d}_stat.branch.tsv'
        path = source_root / f'{number}.zip'
        with zipfile.ZipFile(path, 'w') as archive:
            archive.writestr(name, name.encode())
            info = archive.getinfo(name)
        artifacts.append(store_module.Artifact(logical_path=name, subdir='stat_branch', name=name.split('/')[1],
            generation=1, zip_path=path, member_name=name, size=info.file_size, crc=info.CRC,
            sha256='0' * 64 if failure == 'digest' else None))
    original = store_module.zipfile.ZipFile
    readers = []

    def observed(path, *args, **kwargs):
        archive = original(path, *args, **kwargs)
        if Path(path).parent == source_root:
            readers.append(archive)
            if failure == 'close' and len(readers) == 1:
                original_close = archive.close

                def failing_close():
                    archive.close = original_close
                    original_close()
                    raise OSError('injected source close failure')

                archive.close = failing_close
        return archive

    monkeypatch.setattr(store_module.zipfile, 'ZipFile', observed)
    if failure == 'crc':
        def fail_crc(*args, **kwargs):
            raise store_module.ArchiveStoreError('injected CRC failure')
        monkeypatch.setattr(store_module, '_verify_zip_crc', fail_crc)
    output = tmp_path / 'out/stat_branch'
    output.mkdir(parents=True)
    previous = output / 'stat_branch.pack-000002.zip'
    previous.write_bytes(b'previous complete archive')
    expected_error = OSError if failure == 'close' else store_module.ArchiveStoreError
    with pytest.raises(expected_error):
        store_module._compact_artifact_chunk(tmp_path / 'out', 'stat_branch', artifacts, 'orthogroup', 2)
    assert previous.read_bytes() == b'previous complete archive'
    assert not list(output.glob('.*.partial.*'))
    assert all(reader.fp is None for reader in readers)
