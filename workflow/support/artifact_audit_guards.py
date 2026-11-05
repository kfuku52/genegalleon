"""Private, bounded spill storage for one audit's source identity fences."""
from __future__ import annotations

import json
import sqlite3
import tempfile
from pathlib import Path


class AuditCapacityError(ValueError):
    """An audit-wide resource failure, never a corrupt family manifest."""


def encode(value):
    """Preserve identity types without loading executable serialized objects."""
    def tagged(item):
        if item is None or type(item) in (str, int, float, bool):
            return [type(item).__name__, item]
        if type(item) in (tuple, list):
            return [type(item).__name__, [tagged(child) for child in item]]
        if type(item) is dict:
            pairs = [[tagged(key), tagged(child)] for key, child in item.items()]
            return ['dict', sorted(pairs, key=lambda pair: json.dumps(pair[0], sort_keys=True))]
        raise TypeError('Unsupported audit identity type: ' + type(item).__name__)
    body = json.dumps(tagged(value), ensure_ascii=False, separators=(',', ':'), allow_nan=False).encode()
    if len(body) > 16 * 1024 * 1024:
        raise AuditCapacityError('Audit guard identity exceeds its bounded capacity')
    return body


def decode(body):
    """Read only the tagged JSON types written by this private run's index."""
    def restore(item):
        kind, value = item
        if kind in {'NoneType', 'str', 'int', 'float', 'bool'}:
            return value
        if kind == 'tuple':
            return tuple(restore(child) for child in value)
        if kind == 'list':
            return [restore(child) for child in value]
        if kind == 'dict':
            return {restore(key): restore(child) for key, child in value}
        raise ValueError('Unknown identity type')
    try:
        return restore(json.loads(body))
    except (ValueError, TypeError, RecursionError) as exc:
        raise AuditCapacityError('Audit spill record is invalid') from exc


class GuardSpill:
    """No cross-run reuse: keep every reader, stream identities from private disk.

    SQLite is only an ephemeral index, never workflow or execution history.
    Disk grows with the inventory in the selected temporary filesystem; the
    cache, per-record size and two-million-reader budget remain bounded.
    Callers serialize access through their existing audit lock.
    """
    label = 'source guard'
    def __init__(self):
        self.readers = []
        self.failure = None
        try:
            self.directory = tempfile.TemporaryDirectory(prefix='gg-audit-guards-')
            self.database = sqlite3.connect(Path(self.directory.name) / 'guards.sqlite3',
                                            isolation_level=None, check_same_thread=False)
            self.database.execute('PRAGMA page_size=4096')
            # Do not impose a dataset-size ceiling below SQLite/filesystem
            # limits. Both guards and digests can legitimately exceed 1 GiB.
            self.database.execute('PRAGMA journal_mode=OFF')
            self.database.execute('PRAGMA synchronous=OFF')
            self.database.execute('PRAGMA cache_size=-4096')
            # Keep the key index separate from potentially multi-megabyte
            # inventory identities. A WITHOUT ROWID index carries those large
            # payloads through unrelated key lookups and insertion paths.
            self.database.execute('CREATE TABLE guards (key BLOB PRIMARY KEY NOT NULL, identity BLOB NOT NULL, reader INTEGER NOT NULL UNIQUE)')
            # This connection is private to one audit and never published.
            # A single transaction avoids rewriting index pages per source;
            # the bounded page cache still spills dirty pages to disk.
            self.database.execute('BEGIN')
        except (OSError, sqlite3.Error) as exc:
            self.close()
            raise self._storage_failure(exc) from exc
        except BaseException:
            self.close()
            raise

    def _storage_failure(self, exc):
        directory = getattr(self, 'directory', None)
        location = directory.name if directory is not None else 'the configured temporary filesystem'
        self.failure = (f'Audit {self.label} spill failed in {location}: {exc}. '
                        'Check free space/quota or set TMPDIR to a larger writable filesystem before starting the audit')
        return AuditCapacityError(self.failure)

    def _check_available(self):
        if self.failure is not None:
            raise AuditCapacityError(self.failure)

    def guard(self, key, identity, reader):
        self._check_available()
        encoded_key, expected = encode(key), encode(identity)
        try:
            prior = self.database.execute('SELECT identity,reader FROM guards WHERE key=?', (encoded_key,)).fetchone()
            if prior is not None:
                if prior[0] != expected:
                    raise ValueError('Audit source changed during the audit')
                # The first reader owns this query's original source/snapshot.
                # A later, short-lived store must not replace that context.
                return
            if len(self.readers) >= 2_000_000:
                raise AuditCapacityError(f'Audit {self.label} spill exceeds its bounded capacity')
            self.database.execute('INSERT INTO guards VALUES (?,?,?)', (encoded_key, expected, len(self.readers)))
            self.readers.append(reader)
        except sqlite3.Error as exc:
            # SQLITE_FULL can roll back the whole private transaction. Never
            # accept a later retry with part of the original inventory lost.
            raise self._storage_failure(exc) from exc

    def _records(self):
        self._check_available()
        try:
            count, low, high = self.database.execute('SELECT count(*), min(reader), max(reader) FROM guards').fetchone()
            if (count != len(self.readers) or count and (low != 0 or high != count - 1)):
                raise AuditCapacityError(f'Audit {self.label} spill inventory changed')
            # A record can be 16 MiB. Count-bounded batches retained hundreds
            # of large identities (and two batches during the next fetch).
            # Let the validation submission window bound retained records.
            yield from self.database.execute('SELECT key,identity,reader FROM guards ORDER BY rowid')
        except sqlite3.Error as exc:
            raise self._storage_failure(exc) from exc

    def values(self):
        for _key, identity, index in self._records():
            yield identity, self.readers[index]

    def close(self):
        try:
            if hasattr(self, 'database'):
                self.database.close()
        finally:
            if hasattr(self, 'directory'):
                self.directory.cleanup()


class DigestSpill(GuardSpill):
    """One run's digest identities, with a two-million-reader bound.

    Callables stay in process; keys, signatures, results and directory inventories
    use tagged JSON on private disk. No executable serialization or cross-run reuse.
    """
    label = 'digest'
    def get(self, key):
        self._check_available()
        try:
            row = self.database.execute('SELECT identity,reader FROM guards WHERE key=?', (encode(key),)).fetchone()
        except sqlite3.Error as exc:
            raise self._storage_failure(exc) from exc
        if row is not None:
            before, value, members = decode(row[0])
            return before, value, members, self.readers[row[1]]
        return None

    def items(self):
        for key, identity, index in self._records():
            before, value, members = decode(identity)
            yield decode(key), (before, value, members, self.readers[index])
