"""Typed, content-bound row retrieval for repeatedly queried species tables."""

import hashlib
import inspect
import json
import os
import sqlite3
import tempfile
import zlib
from contextlib import closing
from pathlib import Path

import pandas as pd

try:
    from input_generation_array_state import FreshDigestBatch, atomic_json
    from shared_namespace_lock import namespace_lock
except ImportError:  # package imports
    from .input_generation_array_state import FreshDigestBatch, atomic_json
    from .shared_namespace_lock import namespace_lock


def _payload(row):
    text = json.dumps(row, ensure_ascii=False)
    # Wide numeric tables otherwise duplicate large decimal strings in the
    # index. Independent row compression retains indexed access and cheap
    # fresh content verification; narrow records avoid compression overhead.
    return zlib.compress(text.encode('utf-8'), level=1) if len(text) >= 1024 else text


def _build(database, frame):
    with closing(sqlite3.connect(database)) as connection, connection:
        connection.execute('CREATE TABLE rows (ordinal INTEGER PRIMARY KEY, gene_id TEXT, index_value TEXT, payload TEXT)')
        connection.executemany('INSERT INTO rows VALUES (?, ?, ?, ?)',
            ((i, str(row[1]), json.dumps(row[0]), _payload(row[1:]))
             for i, row in enumerate(frame.itertuples(index=True, name=None))))
        connection.execute('CREATE INDEX gene_ids ON rows(gene_id)')
        connection.execute('CREATE TABLE metadata (payload TEXT)')
        connection.execute('INSERT INTO metadata VALUES (?)', (json.dumps({
            'columns': list(frame.columns), 'dtypes': [str(t) for t in frame.dtypes],
            'index_name': frame.index.name, 'index_dtype': str(frame.index.dtype),
        }),))


def read_subset(path, identifiers, loader, cache_dir=None, *, source_proof=None):
    """Preserve parser types, duplicate rows and original row order/index.

    Without a cache this retains the direct parser path. Each cached read freshly
    hashes the source and the published index, fencing changes through retrieval.
    Cache publication uses the workflow's lock and atomic receipt primitives.
    """
    path = Path(path)
    if not cache_dir:
        if source_proof is not None:
            source_proof.read([path])
        frame = loader(path)
        if source_proof is not None:
            source_proof.check()
        return frame.loc[frame.iloc[:, 0].astype(str).isin(identifiers)].copy()
    root = Path(cache_dir).expanduser().resolve()
    source = source_proof if source_proof is not None else FreshDigestBatch()
    digest = source.read([path])[str(path)]
    if root == path.resolve() or root in path.resolve().parents:
        raise ValueError('Trait index cache must not contain the source')
    key = {'schema': 1, 'source_sha256': digest, 'source_path': str(path.absolute()),
           'pandas': pd.__version__,
           'infer_string': pd.options.future.infer_string, 'string_storage': pd.options.mode.string_storage,
           'reader': hashlib.sha256(inspect.getsource(loader).encode()).hexdigest(),
           'indexer': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    token = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()
    root.mkdir(parents=True, exist_ok=True)
    database, receipt = root / (token + '.sqlite3'), root / (token + '.json')
    with namespace_lock(root / '.locks' / (token + '.lock'), exclusive=True):
        source.check()
        index_proof = FreshDigestBatch()
        try:
            saved = json.loads(receipt.read_text())
            valid = saved['key'] == key and index_proof.read([database])[str(database)] == saved['sha256']
        except (OSError, ValueError, KeyError, TypeError):
            valid = False
        if not valid:
            # The database is a single atomic publication. A crash between it
            # and the receipt leaves a mismatch that the next reader rebuilds.
            with tempfile.TemporaryDirectory(prefix='.trait-index-', dir=root) as temporary:
                staged = Path(temporary) / 'rows.sqlite3'
                _build(staged, loader(path))
                proof = FreshDigestBatch()
                output_hash = proof.read([staged])[str(staged)]
                source.check()
                proof.check()
                os.replace(staged, database)
                atomic_json(receipt, {'key': key, 'sha256': output_hash})
            index_proof = FreshDigestBatch()
            if index_proof.read([database])[str(database)] != output_hash:
                raise OSError('Trait index changed during publication')
        with closing(sqlite3.connect(database.as_uri() + '?mode=ro', uri=True)) as connection:
            metadata = json.loads(connection.execute('SELECT payload FROM metadata').fetchone()[0])
            # Avoid SQLite's parameter limit; retain each duplicate source row.
            connection.execute('CREATE TEMP TABLE wanted (gene_id TEXT PRIMARY KEY)')
            connection.executemany('INSERT OR IGNORE INTO wanted VALUES (?)', ((str(g),) for g in identifiers))
            # Drive the join from requested IDs. Ordering by ordinal must not
            # persuade SQLite to scan every payload in the species table.
            rows = connection.execute('SELECT rows.index_value, rows.payload FROM wanted CROSS JOIN rows INDEXED BY gene_ids ON rows.gene_id = wanted.gene_id ORDER BY rows.ordinal').fetchall()
        frame = pd.DataFrame([json.loads(zlib.decompress(row[1]) if isinstance(row[1], bytes) else row[1]) for row in rows], columns=metadata['columns'],
                             index=pd.Index([json.loads(row[0]) for row in rows], dtype=metadata['index_dtype'], name=metadata['index_name']))
        for column, dtype in zip(frame.columns, metadata['dtypes'], strict=True):
            frame[column] = frame[column].astype(dtype)
        frame = frame.loc[frame.iloc[:, 0].astype(str).isin(identifiers)].copy()
        source.check()
        index_proof.check()
    return frame
