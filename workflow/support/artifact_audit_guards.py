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

    SQLite is only an ephemeral index, never workflow or execution history. Both
    its 512 MiB disk budget and two-million-reader budget fail closed.
    Callers serialize access through their existing audit lock.
    """
    label = 'source guard'
    max_pages = 131072

    def __init__(self):
        self.directory = tempfile.TemporaryDirectory(prefix='gg-audit-guards-')
        self.readers = []
        try:
            self.database = sqlite3.connect(Path(self.directory.name) / 'guards.sqlite3',
                                            isolation_level=None, check_same_thread=False)
            self.database.execute('PRAGMA page_size=4096')
            self.database.execute(f'PRAGMA max_page_count={self.max_pages}')
            self.database.execute('PRAGMA journal_mode=OFF')
            self.database.execute('PRAGMA synchronous=OFF')
            self.database.execute('PRAGMA cache_size=-4096')
            # Keep the key index separate from potentially multi-megabyte
            # inventory identities. A WITHOUT ROWID index carries those large
            # payloads through unrelated key lookups and insertion paths.
            self.database.execute('CREATE TABLE guards (key BLOB PRIMARY KEY NOT NULL, identity BLOB NOT NULL, reader INTEGER NOT NULL UNIQUE)')
        except BaseException:
            self.close()
            raise

    def guard(self, key, identity, reader):
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
            raise AuditCapacityError(f'Audit {self.label} spill failed') from exc

    def _records(self):
        try:
            count, low, high = self.database.execute('SELECT count(*), min(reader), max(reader) FROM guards').fetchone()
            if (count != len(self.readers) or count and (low != 0 or high != count - 1)):
                raise AuditCapacityError(f'Audit {self.label} spill inventory changed')
            cursor = self.database.execute('SELECT key,identity,reader FROM guards ORDER BY rowid')
            while rows := cursor.fetchmany(512):
                yield from rows
        except sqlite3.Error as exc:
            raise AuditCapacityError(f'Audit {self.label} spill failed') from exc

    def values(self):
        for _key, identity, index in self._records():
            yield identity, self.readers[index]

    def close(self):
        if hasattr(self, 'database'):
            self.database.close()
        self.directory.cleanup()


class DigestSpill(GuardSpill):
    """One run's digest identities, with a 1 GiB disk and two-million-reader bound.

    Callables stay in process; keys, signatures, results and directory inventories
    use tagged JSON on private disk. No executable serialization or cross-run reuse.
    """
    label = 'digest'
    max_pages = 262144

    def get(self, key):
        try:
            row = self.database.execute('SELECT identity,reader FROM guards WHERE key=?', (encode(key),)).fetchone()
        except sqlite3.Error as exc:
            raise AuditCapacityError('Audit digest spill failed') from exc
        if row is not None:
            before, value, members = decode(row[0])
            return before, value, members, self.readers[row[1]]
        return None

    def items(self):
        for key, identity, index in self._records():
            before, value, members = decode(identity)
            yield decode(key), (before, value, members, self.readers[index])
