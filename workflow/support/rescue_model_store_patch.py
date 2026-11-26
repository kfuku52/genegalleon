"""Immutable admission updates preserving raw ordinals and partial references."""
import bisect
import copy
import hashlib
import json
import math
import os
import shutil
import tempfile
import zipfile
from collections import deque
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from contextlib import ExitStack, nullcontext
from multiprocessing import get_context
from pathlib import Path

if __package__:
    from . import rescue_model_store as store
    from .stage_output_hashes import hash_paths
else:
    import rescue_model_store as store
    from stage_output_hashes import hash_paths


def _copy(source, destination, expected, *, verify_destination=True):
    checksum = hashlib.sha256()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open('rb') as inp, destination.open('xb') as out:
        for block in iter(lambda: inp.read(1024 * 1024), b''):
            checksum.update(block)
            out.write(block)
    if checksum.hexdigest() != expected or (verify_destination and store._digest(destination) != expected):
        raise ValueError('Model-store changed during immutable copying')


def copy_legacy_model_archive(source, destination, frozen_key):
    """Copy a migrated legacy store into a new publication in native ZIP layout."""
    root, source_store = store._layout(source)
    if not os.path.lexists(source_store / store.LEGACY_SHARD_ARCHIVE):
        return False
    manifest, manifest_sha = store._manifest(source_store)
    archive = store._legacy_archive(source_store, manifest)
    if archive is None:
        return False
    store._check_key(source, frozen_key, 'models')
    receipt, _ = store._receipt(root)
    prefix = source_store.relative_to(root).as_posix()
    members = {f'{prefix}/{name}': sha for name, sha in manifest['files'].items()}
    if (receipt['files'].get(f'{prefix}/manifest.json') != manifest_sha
            or any(receipt['files'].get(name) != sha for name, sha in members.items())):
        raise ValueError('Copy input is not bound to its producer receipt')

    target = Path(destination) / 'model_store'
    target.mkdir()
    with store.verified_legacy_shards(root, members) as archived:
        for name, sha in members.items():
            if name not in archived:
                _copy(store._member(root, name), Path(destination) / name, sha)
        hashes = manifest['files']
        hashes[store.SHARD_ARCHIVE] = store._digest(archive)
        _copy(archive, target / store.SHARD_ARCHIVE, hashes[store.SHARD_ARCHIVE])
    manifest['shard_archive'] = store.SHARD_ARCHIVE
    (target / 'manifest.json').write_bytes(store._json(manifest) + b'\n')
    store._check_key(source, frozen_key, 'models')
    return True


def _same_json_value(a, b):
    """Compare immutable values before the final envelope's JSON validation.

    Equal built-in scalars need no temporary JSON strings. Containers and
    subclasses retain the serialized comparison, including field order and
    JSON's distinctions between numeric types and the signs of zero.
    """
    kind = type(a)
    if kind is type(b):
        if kind in (str, int, bool, type(None)):
            return a == b
        if kind is float and math.isfinite(a) and math.isfinite(b):
            return a == b and (a != 0.0 or math.copysign(1.0, a) == math.copysign(1.0, b))
    return store._json(a) == store._json(b)


def _patch_shard(source_store, temporary, shard, expected, codec, level, changes, shard_root=None):
    normalized = set()
    with ExitStack() as resources:
        source_bodies = None
        name = shard['path']
        if shard_root is None:
            shard_root = store._worker_shards
        seen = 0
        with store._writer(temporary / name, codec, level) as out:
            for envelope in store._lines(shard_root / name, codec, expected):
                ordinal = shard['first'] + seen
                if envelope.get('ordinal') != ordinal:
                    raise ValueError('Patch source candidate order changed')
                if ordinal in changes:
                    # Serialize only the current record. A cohort can have
                    # millions of changed candidates with large DNA bodies.
                    # The final write validates the envelope once, with
                    # its ordinal, before publishing any candidate bytes.
                    replacement_body, replacement = store._pack(changes[ordinal], _validate_envelope=False)
                    replacement['ordinal'] = ordinal
                    if replacement['context'].get('partial_evidence') != envelope['context'].get('partial_evidence'):
                        raise ValueError('Admission patch changed predictor/DNA or partial evidence')
                    if replacement['body'] != envelope['body']:
                        # Earlier v1 producers put donor metrics in shared
                        # bodies. Prove the same predictor values before
                        # normalizing just the changed candidate's layout.
                        if source_bodies is None:
                            db = store._readonly(source_store / 'bodies.sqlite')
                            resources.callback(db.close)
                            source_bodies = store._Bodies(db, codec)
                        raw_body = source_bodies.read(envelope['body'])
                        original = store._unpack(envelope, raw_body)
                        body_keys = json.loads(raw_body)
                        if (store._pack(original)[0] != replacement_body
                                or any(k not in changes[ordinal]
                                       or not _same_json_value(original[k], changes[ordinal][k])
                                       for k in body_keys)):
                            raise ValueError('Admission patch changed predictor/DNA values')
                        normalized.add(ordinal)
                    else:
                        immutable = ('coverage', 'identity', 'query_span_coverage', 'query_start', 'query_end',
                                     'query_length', 'query', 'id', 'paf', 'raw_prediction', 'genetic_code')
                        a, b = envelope['context'], replacement['context']
                        if any((k in a) != (k in b) or (k in a and not _same_json_value(a[k], b[k]))
                               for k in immutable):
                            raise ValueError('Admission patch changed predictor/DNA values')
                    envelope = replacement
                out.write(store._json(envelope) + b'\n')
                seen += 1
        if seen != shard['count']:
            raise ValueError('Patch source shard count changed')
    return normalized


def _initialize_patch_shards(source, archive):
    store._initialize_worker_shards(source, archive)


def _patch_shards(source, output, shards, manifest, ordered_changes, changes, workers, shard_root, archive):
    def args(shard, first, last):
        name = shard['path']
        return (source, output, shard, manifest['files'][name], manifest['codec'],
                manifest['compression_level'], {i: changes[i] for i in ordered_changes[first:last]})

    if workers == 1 or len(shards) < workers * 2:
        for shard, first, last in shards:
            yield from _patch_shard(*args(shard, first, last), shard_root=shard_root)
        return
    pending = deque()
    with ProcessPoolExecutor(max_workers=workers, mp_context=get_context('spawn'),
                             initializer=_initialize_patch_shards,
                             initargs=(source, archive)) as pool:
        for shard, first, last in shards:
            pending.append(pool.submit(_patch_shard, *args(shard, first, last)))
            if len(pending) >= workers * 2:
                yield from pending.popleft().result()
        while pending:
            yield from pending.popleft().result()


APPEND_SHARD_BYTES = store.SHARD_BYTES
BODY_DB_BATCH = store.BODY_DB_BATCH


def _append_models(directory, count, appended, manifest, body_db, body_cache, workers):
    """Pack once in order, overlapping bounded shard compression with SQLite writes."""
    codec, level = manifest['codec'], manifest['compression_level']
    extension = 'gz' if codec == 'gzip' else 'zst'
    shards, partials = [], []
    packed = bytearray()
    bodies = {}
    pending = deque()
    first, rows = count, 0
    writers = min(3, workers - 1)

    def flush(executor):
        nonlocal packed, first, rows
        path = f'models-append-{first:012d}.jsonl.{extension}'
        args = directory / path, packed, codec, level
        if executor is None:
            store._write_shard(*args)
        else:
            if len(pending) == writers:
                pending.popleft().result()
            pending.append(executor.submit(store._write_shard, *args))
        shards.append({'path': path, 'first': first, 'count': rows})
        packed = bytearray()
        first += rows
        rows = 0

    context = ThreadPoolExecutor(max_workers=writers) if writers else nullcontext()
    with context as executor:
        for ordinal, record in enumerate(appended, count):
            body, sha, row, _ = store._encode_model(record, ordinal)
            if bodies.setdefault(sha, body) != body:
                raise ValueError('Model-store body hash collision')
            if len(bodies) == BODY_DB_BATCH or len(packed) + len(row) > APPEND_SHARD_BYTES:
                store._insert_bodies(body_db, bodies, body_cache, level)
                bodies.clear()
            if packed and len(packed) + len(row) > APPEND_SHARD_BYTES:
                flush(executor)
            if record.get('partial_evidence', {}).get('partial'):
                partials.append((ordinal, len(manifest['shards']) + len(shards), rows))
            packed.extend(row)
            rows += 1
        store._insert_bodies(body_db, bodies, body_cache, level)
        if packed:
            flush(executor)
        for result in pending:
            result.result()
    return shards, partials


def patch_model_store(source, destination, changes, *, appended=(), revisions=(), frozen_key=None,
                      append_partial_evidence=False, publication_workers=1):
    """Amend decisions without rewriting unselected DNA bodies or candidates.

    Existing coordinates, sequences, predictor values and partial evidence must
    remain exact. New independently validated records may only be appended.
    Unknown fields and original shard/partial ordering survive unchanged. The
    enclosing ordinary rescue stage supplies the new producer receipt.
    Independent changed shards are encoded in bounded worker batches.
    """
    if type(publication_workers) is not int or publication_workers < 1:
        raise ValueError('Model-store publication workers must be a positive integer')
    root, source_store = store._layout(source)
    destination = Path(destination)
    if (destination.resolve() == root.resolve() or root.resolve() in destination.resolve().parents
            or destination.resolve() in root.resolve().parents):
        raise ValueError('Patch destination overlaps its source')
    key = frozen_key or store.frozen_model_store_key(source)
    store._check_key(source, key, 'models')
    manifest, manifest_sha = store._manifest(source_store)
    receipt, receipt_sha = store._receipt(root)
    prefix = source_store.relative_to(root).as_posix()
    members = {f'{prefix}/{name}': sha for name, sha in store._physical_files(manifest).items()}
    members[f'{prefix}/manifest.json'] = manifest_sha
    if any(receipt['files'].get(name) != sha for name, sha in members.items()):
        raise ValueError('Patch input is not bound to its producer receipt')
    legacy_archive = store._legacy_archive(source_store, manifest)
    archived = store._legacy_members(root, source_store, manifest, members) if legacy_archive else {}
    paths = {name: store._member(root, name) for name in members if name not in archived}
    if legacy_archive:
        paths[legacy_archive.name] = legacy_archive
    stats = {name: store._stat(path) for name, path in paths.items()}
    count = manifest['counts']['models']
    if any(type(i) is not int or not 0 <= i < count for i in changes):
        raise ValueError('Patch ordinal outside the original model publication')
    appended = list(appended)
    if not append_partial_evidence and any(m.get('partial_evidence', {}).get('partial') for m in appended):
        raise ValueError('Admission append cannot change the partial reference stream')
    target = Path(destination) / 'model_store'
    if target.exists() or target.is_symlink():
        raise ValueError('Patch destination already exists')
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix='.model-patch-', dir=target.parent))
    databases = []
    try:
        for name in ('accepted.sqlite', 'revision.sqlite'):
            if store._digest(source_store / name) != manifest['files'][name]:
                raise ValueError('Model-store overlay changed before patching')
        accepted = dict(store._iter(source, 'accepted', with_ordinals=True))
        for ordinal in changes:
            accepted.pop(ordinal, None)
        for ordinal, model in changes.items():
            if model.get('status') == 'accepted':
                accepted[ordinal] = model
        for offset, model in enumerate(appended):
            if model.get('status') == 'accepted':
                accepted[count + offset] = model
        output = copy.deepcopy(manifest)
        normalized = set()
        ordered_changes = sorted(changes)
        unchanged, affected_shards = [], []
        for shard in manifest['shards']:
            first = bisect.bisect_left(ordered_changes, shard['first'])
            last = bisect.bisect_left(ordered_changes, shard['first'] + shard['count'], lo=first)
            if first == last:
                unchanged.append(source_store / shard['path'])
            else:
                affected_shards.append((shard, first, last))

        archive = temporary / store.SHARD_ARCHIVE
        with store._open_shards(source_store, manifest, legacy_archive=legacy_archive) as shard_root:
            with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_STORED) as handle:
                for path in unchanged:
                    name = str(path.relative_to(source_store))
                    store._copy_zip_member(shard_root / name, handle, name, manifest['files'][name])
            normalized.update(_patch_shards(source_store, temporary, affected_shards, manifest,
                                            ordered_changes, changes, publication_workers, shard_root,
                                            legacy_archive.name if legacy_archive else manifest.get('shard_archive')))
        deferred = {}
        _copy(source_store / 'bodies.sqlite', temporary / 'bodies.sqlite', manifest['files']['bodies.sqlite'])
        body_db = store.sqlite3.connect(temporary / 'bodies.sqlite')
        databases.append(body_db)
        body_db.execute('PRAGMA journal_mode=OFF')
        body_db.execute('PRAGMA synchronous=OFF')
        body_cache = store._Bodies(body_db, manifest['body_codec'])
        for ordinal in sorted(normalized):
            body, _ = store._pack(changes[ordinal])
            store._insert_body(body_db, body, manifest['compression_level'], manifest['codec'], body_cache)
        new_shards, new_partials = _append_models(temporary, count, appended, manifest, body_db, body_cache,
                                                   publication_workers)
        output['counts']['bodies'] = body_db.execute('SELECT COUNT(*) FROM bodies').fetchone()[0]
        body_db.commit()
        body_db.close()
        databases.remove(body_db)
        del body_cache
        partial_expected = manifest['files']['partial.sqlite']
        _copy(source_store / 'partial.sqlite', temporary / 'partial.sqlite', partial_expected,
              verify_destination=bool(new_partials))
        if new_partials:
            partial_db = store.sqlite3.connect(temporary / 'partial.sqlite')
            databases.append(partial_db)
            original_count = partial_db.execute('SELECT COUNT(*) FROM refs').fetchone()[0]
            if original_count != manifest['counts']['partial']:
                raise ValueError('Original partial reference count changed')
            partial_db.executemany('INSERT INTO refs VALUES (?,?,?,?)',
                                   ((position, *ref) for position, ref in enumerate(new_partials, original_count)))
            partial_db.commit()
            partial_db.close()
            databases.remove(partial_db)
            output['counts']['partial'] += len(new_partials)
        else:
            deferred['partial.sqlite'] = partial_expected
        known = {'bodies.sqlite', 'partial.sqlite', 'accepted.sqlite', 'revision.sqlite',
                 *[s['path'] for s in manifest['shards']]}
        if 'shard_archive' in manifest:
            known.add(manifest['shard_archive'])
        for name in store._physical_files(manifest).keys() - known:
            expected = manifest['files'][name]
            _copy(source_store / name, temporary / name, expected, verify_destination=False)
            deferred[name] = expected
        output['shards'].extend(new_shards)
        if new_partials:
            previous = output['scopes']['partial'][1:]
            partial_shards = sorted({shard for _, shard, _ in new_partials})
            output['scopes']['partial'] = ['partial.sqlite', 'bodies.sqlite',
                *(name for name in previous if name != 'bodies.sqlite'),
                *(output['shards'][shard]['path'] for shard in partial_shards)]
        overlay_counts = {}
        for kind, rows in (('accepted', sorted(accepted.items())), ('revision', enumerate(revisions))):
            db = store._overlay(temporary / (kind + '.sqlite'))
            databases.append(db)
            overlay_counts[kind] = 0
            for ordinal, model in rows:
                body, envelope = store._pack(model)
                envelope['ordinal'] = ordinal
                store._add_overlay(db, ordinal, body, envelope, manifest['compression_level'], manifest['codec'])
                overlay_counts[kind] += 1
            db.commit()
            db.close()
            databases.remove(db)
        output['counts'].update(models=count + len(appended), **overlay_counts)
        output['scopes']['models'] = ['bodies.sqlite', *[s['path'] for s in output['shards']]]
        files = (p for p in temporary.rglob('*') if p.is_file() and p != archive)
        hashes = hash_paths(temporary, files, workers=publication_workers, hash_function=store._digest)
        if any(hashes.get(name) != expected for name, expected in deferred.items()):
            raise ValueError('Unchanged model-store member differs at publication')
        with zipfile.ZipFile(archive, 'a', compression=zipfile.ZIP_STORED) as handle:
            for shard in output['shards']:
                name = shard['path']
                path = temporary / name
                if name in hashes:
                    store._copy_zip_member(path, handle, name, hashes[name])
                    path.unlink()
                else:
                    hashes[name] = manifest['files'][name]
        if any(sha != hashes[name] for name, sha in store._archive_shard_hashes(temporary, output['shards']).items()):
            raise ValueError('Copied model-store ZIP member differs')
        hashes[archive.name] = store._digest(archive)
        output['shard_archive'] = archive.name
        output['files'] = dict(sorted(hashes.items()))
        (temporary / 'manifest.json').write_bytes(store._json(output) + b'\n')
        store._check_key(source, key, 'models')
        if (store._receipt(root)[1] != receipt_sha
                or any(path.is_symlink() or store._stat(path) != stats[name] for name, path in paths.items())):
            raise ValueError('Patch source generation changed during publication')
        store._sync_files((p for p in temporary.rglob('*') if p.is_file()), publication_workers)
        os.rename(temporary, target)
        store._sync_directory(target.parent)
        return output
    except BaseException:
        for db in databases:
            db.close()
        shutil.rmtree(temporary, ignore_errors=True)
        raise
