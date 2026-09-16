"""Opt-in URL/member cache shared by independent input-generation plans."""
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def materialize_shared_download(download, arguments):
    # Import lazily: locking owns both the public transport and lock helpers.
    from .local import validate_gzip_with_cache
    from .locking import _fsync_directory, acquire_download_lock, release_download_lock
    root = Path(os.environ['GG_DOWNLOAD_SHARED_CACHE_DIR']).expanduser().resolve()
    if Path(arguments["destination"]).resolve().is_relative_to(root):
        raise ValueError("download destination must be outside the shared cache")
    root.mkdir(parents=True, exist_ok=True)
    # Include request headers so data requested under different credentials cannot
    # alias. Neither URLs nor header values are written to the receipt.
    identity = json.dumps([arguments['url'], arguments['archive_member'],
                           sorted((k.lower(), v) for k, v in (arguments['headers'] or {}).items())],
                          separators=(',', ':'))
    key = hashlib.sha256(identity.encode()).hexdigest()
    entry = root / key
    entry.mkdir(exist_ok=True)
    suffix = Path(arguments['destination']).suffix.lower()
    suffix = suffix if suffix in ('.gz', '.zip') else '.bin'
    payload = entry / ('payload' + suffix)
    receipt = entry / ('receipt' + suffix + '.json')
    lock = entry / 'publish.lock'
    ownership = acquire_download_lock(lock, arguments['lock_stale_seconds'], arguments['warnings'], arguments['lock_context'])
    try:
        valid = False
        fetched = False
        if payload.exists() and not arguments['overwrite']:
            try:
                evidence = json.loads(receipt.read_text())
                valid = (evidence['schema_version'] == 2
                         and type(evidence['size']) is int and evidence['size'] > 0
                         and evidence['size'] == payload.stat().st_size
                         and evidence['sha256'] == sha256_file(payload))
            except (OSError, ValueError, KeyError, TypeError):
                pass
        if not valid:
            # Keep incomplete payloads for the transport's identity-aware resume.
            # A bad completed cache must never be accepted by the skip path.
            options = dict(arguments, destination=payload, overwrite=arguments['overwrite'] or payload.exists(), validation_cache=None,
                           validation_key='', validation_relative_target='')
            download(**options)
            fetched = True
            evidence = {'schema_version': 2, 'size': payload.stat().st_size,
                        'sha256': sha256_file(payload)}
            temp = receipt.with_suffix('.tmp')
            with open(temp, 'w', encoding='utf-8') as out:
                json.dump(evidence, out)
                out.write('\n')
                out.flush()
                os.fsync(out.fileno())
            temp.replace(receipt)
            _fsync_directory(entry)
        destination = Path(arguments['destination'])
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() and not arguments['overwrite']:
            if (destination.stat().st_size == payload.stat().st_size
                    and sha256_file(destination) == evidence['sha256']):
                return False, fetched
        # Copy rather than hardlink: downstream edits must not poison the cache.
        descriptor, temporary = tempfile.mkstemp(prefix=destination.name + '.tmp.', dir=destination.parent)
        try:
            with os.fdopen(descriptor, 'wb') as out, open(payload, 'rb') as source:
                shutil.copyfileobj(source, out, 1024 * 1024)
                out.flush()
                os.fsync(out.fileno())
            if sha256_file(temporary) != evidence["sha256"]:
                raise OSError("shared download copy checksum mismatch")
            error = validate_gzip_with_cache(Path(temporary), expected_path=destination)
            if error is not None:
                raise OSError('shared download failed validation: {}'.format(error))
            os.replace(temporary, destination)
            _fsync_directory(destination.parent)
        finally:
            Path(temporary).unlink(missing_ok=True)
        return True, fetched
    finally:
        release_download_lock(lock, ownership)
