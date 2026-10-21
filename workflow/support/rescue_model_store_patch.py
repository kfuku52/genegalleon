"""Immutable admission updates preserving raw ordinals and partial references."""
import copy
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path

try:
    import rescue_model_store as store
except ImportError:
    from . import rescue_model_store as store


def _copy(source, destination, expected):
    checksum = hashlib.sha256()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open('rb') as inp, destination.open('xb') as out:
        for block in iter(lambda: inp.read(1024 * 1024), b''):
            checksum.update(block)
            out.write(block)
    if checksum.hexdigest() != expected or store._digest(destination) != expected:
        raise ValueError('Model-store changed during immutable copying')


def patch_model_store(source, destination, changes, *, appended=(), revisions=(), frozen_key=None,
                      append_partial_evidence=False):
    """Amend decisions without rewriting unselected DNA bodies or candidates.

    Existing coordinates, sequences, predictor values and partial evidence must
    remain exact. New independently validated records may only be appended.
    Unknown fields and original shard/partial ordering survive unchanged. The
    enclosing ordinary rescue stage supplies the new producer receipt.
    """
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
    members = {f'{prefix}/{name}': sha for name, sha in manifest['files'].items()}
    members[f'{prefix}/manifest.json'] = manifest_sha
    if any(receipt['files'].get(name) != sha for name, sha in members.items()):
        raise ValueError('Patch input is not bound to its producer receipt')
    stats = {name: store._stat(store._member(root, name)) for name in members}
    count = manifest['counts']['models']
    if any(type(i) is not int or not 0 <= i < count for i in changes):
        raise ValueError('Patch ordinal outside the original model publication')
    appended = list(appended)
    new_partials = [i for i, m in enumerate(appended) if m.get('partial_evidence', {}).get('partial')]
    if new_partials and not append_partial_evidence:
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
        source_db = store._readonly(source_store / 'bodies.sqlite')
        databases.append(source_db)
        source_bodies = store._Bodies(source_db, manifest['body_codec'])
        for ordinal in changes:
            accepted.pop(ordinal, None)
        updates = {**changes, **{count + i: m for i, m in enumerate(appended)}}
        packed = {}
        for ordinal, model in updates.items():
            body, envelope = store._pack(model)
            envelope['ordinal'] = ordinal
            packed[ordinal] = (body, envelope)
            if model.get('status') == 'accepted':
                accepted[ordinal] = model
        output = copy.deepcopy(manifest)
        normalized = set()
        for shard in manifest['shards']:
            name = shard['path']
            affected = {i for i in changes if shard['first'] <= i < shard['first'] + shard['count']}
            if not affected:
                _copy(source_store / name, temporary / name, manifest['files'][name])
                continue
            seen = 0
            with store._writer(temporary / name, manifest['codec'], manifest['compression_level']) as out:
                for envelope in store._lines(source_store / name, manifest['codec'], manifest['files'][name]):
                    ordinal = shard['first'] + seen
                    if envelope.get('ordinal') != ordinal:
                        raise ValueError('Patch source candidate order changed')
                    if ordinal in affected:
                        replacement = packed[ordinal][1]
                        if replacement['context'].get('partial_evidence') != envelope['context'].get('partial_evidence'):
                            raise ValueError('Admission patch changed predictor/DNA or partial evidence')
                        if replacement['body'] != envelope['body']:
                            # Earlier v1 producers put donor metrics in shared
                            # bodies. Prove the same predictor values before
                            # normalizing just the changed candidate's layout.
                            raw_body = source_bodies.read(envelope['body'])
                            original = store._unpack(envelope, raw_body)
                            body_keys = json.loads(raw_body)
                            if (store._pack(original)[0] != packed[ordinal][0]
                                    or any(k not in changes[ordinal]
                                           or store._json(original[k]) != store._json(changes[ordinal][k])
                                           for k in body_keys)):
                                raise ValueError('Admission patch changed predictor/DNA values')
                            normalized.add(ordinal)
                        else:
                            immutable = ('coverage', 'identity', 'query_span_coverage', 'query_start', 'query_end',
                                         'query_length', 'query', 'id', 'paf', 'raw_prediction', 'genetic_code')
                            a, b = envelope['context'], replacement['context']
                            if any((k in a) != (k in b) or (k in a and store._json(a[k]) != store._json(b[k]))
                                   for k in immutable):
                                raise ValueError('Admission patch changed predictor/DNA values')
                        envelope = replacement
                    out.write(store._json(envelope) + b'\n')
                    seen += 1
            if seen != shard['count']:
                raise ValueError('Patch source shard count changed')
        source_db.close()
        databases.remove(source_db)
        _copy(source_store / 'bodies.sqlite', temporary / 'bodies.sqlite', manifest['files']['bodies.sqlite'])
        body_db = store.sqlite3.connect(temporary / 'bodies.sqlite')
        databases.append(body_db)
        body_db.execute('PRAGMA journal_mode=OFF')
        body_db.execute('PRAGMA synchronous=OFF')
        for ordinal in [*sorted(normalized), *range(count, count + len(appended))]:
            store._insert_body(body_db, packed[ordinal][0], manifest['compression_level'], manifest['codec'])
        output['counts']['bodies'] = body_db.execute('SELECT COUNT(*) FROM bodies').fetchone()[0]
        body_db.commit()
        body_db.close()
        databases.remove(body_db)
        _copy(source_store / 'partial.sqlite', temporary / 'partial.sqlite', manifest['files']['partial.sqlite'])
        if new_partials:
            partial_db = store.sqlite3.connect(temporary / 'partial.sqlite')
            databases.append(partial_db)
            original_count = partial_db.execute('SELECT COUNT(*) FROM refs').fetchone()[0]
            if original_count != manifest['counts']['partial']:
                raise ValueError('Original partial reference count changed')
            for position, line in enumerate(new_partials, original_count):
                partial_db.execute('INSERT INTO refs VALUES (?,?,?,?)',
                                   (position, count + line, len(manifest['shards']), line))
            partial_db.commit()
            partial_db.close()
            databases.remove(partial_db)
            output['counts']['partial'] += len(new_partials)
        known = {'bodies.sqlite', 'partial.sqlite', 'accepted.sqlite', 'revision.sqlite',
                 *[s['path'] for s in manifest['shards']]}
        for name in manifest['files'].keys() - known:
            _copy(source_store / name, temporary / name, manifest['files'][name])
        if appended:
            suffix = 'zst' if manifest['codec'] == 'zstd' else 'gz'
            member = f'models-{len(manifest["shards"]):06d}.jsonl.{suffix}'
            with store._writer(temporary / member, manifest['codec'], manifest['compression_level']) as out:
                for ordinal in range(count, count + len(appended)):
                    out.write(store._json(packed[ordinal][1]) + b'\n')
            output['shards'].append({'path': member, 'first': count, 'count': len(appended)})
            if new_partials:
                previous = output['scopes']['partial'][1:]
                output['scopes']['partial'] = ['partial.sqlite', 'bodies.sqlite',
                    *(name for name in previous if name != 'bodies.sqlite'), member]
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
        output['files'] = {str(p.relative_to(temporary)): store._digest(p)
                           for p in sorted(temporary.rglob('*')) if p.is_file()}
        (temporary / 'manifest.json').write_bytes(store._json(output) + b'\n')
        store._check_key(source, key, 'models')
        if (store._receipt(root)[1] != receipt_sha
                or any(store._stat(store._member(root, name)) != before for name, before in stats.items())):
            raise ValueError('Patch source generation changed during publication')
        for path in temporary.rglob('*'):
            if path.is_file():
                with path.open('rb') as handle:
                    os.fsync(handle.fileno())
        os.rename(temporary, target)
        store._sync_directory(target.parent)
        return output
    except BaseException:
        for db in databases:
            db.close()
        shutil.rmtree(temporary, ignore_errors=True)
        raise
