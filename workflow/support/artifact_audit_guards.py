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


class GuardSpill:
    """No cross-run reuse: keep every reader, stream identities from private disk.

    SQLite is only an ephemeral index, never workflow or execution history. Both
    its 512 MiB disk budget and two-million-reader budget fail closed.
    Callers serialize access through their existing audit lock.
    """
    def __init__(self):
        self.directory = tempfile.TemporaryDirectory(prefix='gg-audit-guards-')
        self.readers = []
        try:
            self.database = sqlite3.connect(Path(self.directory.name) / 'guards.sqlite3',
                                            isolation_level=None, check_same_thread=False)
            self.database.execute('PRAGMA page_size=4096')
            self.database.execute('PRAGMA max_page_count=131072')
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
                self.readers[prior[1]] = reader
                return
            if len(self.readers) >= 2_000_000:
                raise AuditCapacityError('Audit source guard spill exceeds its bounded capacity')
            self.database.execute('INSERT INTO guards VALUES (?,?,?)', (encoded_key, expected, len(self.readers)))
            self.readers.append(reader)
        except sqlite3.Error as exc:
            raise AuditCapacityError('Audit source guard spill failed') from exc

    def values(self):
        try:
            count, low, high = self.database.execute('SELECT count(*), min(reader), max(reader) FROM guards').fetchone()
            if (count != len(self.readers) or count and (low != 0 or high != count - 1)):
                raise AuditCapacityError('Audit source guard spill inventory changed')
            cursor = self.database.execute('SELECT identity,reader FROM guards')
            while rows := cursor.fetchmany(512):
                for identity, index in rows:
                    yield identity, self.readers[index]
        except sqlite3.Error as exc:
            raise AuditCapacityError('Audit source guard spill failed') from exc

    def close(self):
        if hasattr(self, 'database'):
            self.database.close()
        self.directory.cleanup()
