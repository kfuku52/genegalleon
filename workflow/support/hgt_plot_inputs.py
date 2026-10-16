#!/usr/bin/env python3
"""Private, source-bound family index for HGT plotting (never scientific output).

Read each candidate table once with the normal whole-table pandas inference.
Keep exact typed frames in a private SQLite spool using data-only JSON.
Family readers verify each payload digest and restore the whole-table dtypes.
Original TSVs remain the declared provenance inputs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sqlite3
import stat
import tempfile
from collections import Counter, defaultdict
from contextlib import closing
from pathlib import Path

import pandas
from content_digest_cache import cached_sha256_file
from shared_namespace_lock import namespace_lock

SCHEMA = 1


def signature(path):
    path = Path(path).absolute()
    info = path.lstat()
    if not stat.S_ISREG(info.st_mode) or info.st_uid != os.getuid():
        raise ValueError(f'Expected an owned regular file: {path}')
    return [info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns]


def sources(branch, gene):
    result = {}
    for key, path in [('branch', branch), ('gene', gene)]:
        path = Path(path).absolute()
        before = signature(path)
        digest, size = cached_sha256_file(path)
        if signature(path) != before:
            raise ValueError(f'Candidate input changed while hashing: {path}')
        result[key] = dict(path=str(path), signature=before, sha256=digest, size=size)
    return result


def connect_readonly(path):
    signature(path)
    connection = sqlite3.connect(Path(path).absolute().as_uri() + '?mode=ro', uri=True)
    connection.execute('PRAGMA query_only=ON')
    return connection


def metadata(connection):
    return json.loads(connection.execute('SELECT value FROM metadata').fetchone()[0])


def validate_sources(proof, branch, gene):
    if proof['schema'] != SCHEMA:
        raise ValueError('Unsupported HGT candidate index schema')
    for key, path in [('branch', branch), ('gene', gene)]:
        expected = proof['sources'][key]
        if str(Path(path).absolute()) != expected['path'] or signature(path) != expected['signature']:
            raise ValueError(f'Candidate input changed since index preparation: {path}')


def majority_prefixes(path):
    # Preserve the existing awk nonempty-taxonomy count and lexical tie-break,
    # including strings such as NA that pandas normally recognizes as missing.
    csv.field_size_limit(100_000_000)
    counts = defaultdict(Counter)
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row.get('orthogroup') and row.get('gene_taxon'):
                counts[row['orthogroup']][row['gene_taxon']] += 1
    return {family: min(values, key=lambda key: (-values[key], key)).replace(' ', '_') + '_'
            for family, values in counts.items()}


def encode_frame(frame):
    def scalar(value):
        if isinstance(value, float) and not math.isfinite(value):
            return {'float_hex': value.hex()}
        return None if pandas.isna(value) else value
    payload = dict(columns=frame.columns.tolist(), dtypes=[str(t) for t in frame.dtypes],
                   index=frame.index.tolist(), index_dtype=str(frame.index.dtype),
                   data=[[scalar(v) for v in row] for row in frame.itertuples(index=False, name=None)])
    return json.dumps(payload, separators=(',', ':'), allow_nan=False).encode()


def decode_frame(data):
    payload = json.loads(data)
    def scalar(value):
        if isinstance(value, dict):
            if set(value) != {'float_hex'}:
                raise ValueError('Invalid typed HGT index scalar')
            return float.fromhex(value['float_hex'])
        return float('nan') if value is None else value
    frame = pandas.DataFrame({
        column: pandas.Series([scalar(row[i]) for row in payload['data']], dtype=dtype)
        for i, (column, dtype) in enumerate(zip(payload['columns'], payload['dtypes'], strict=True))})
    frame.index = pandas.Index(payload['index'], dtype=payload['index_dtype'])
    return frame


def prepare(path, branch, gene):
    path = Path(path).absolute()
    path.parent.mkdir(parents=True, exist_ok=True)
    with namespace_lock(path.with_suffix('.build.guard'), exclusive=True):
        proof = dict(schema=SCHEMA, pandas_version=pandas.__version__, sources=sources(branch, gene))
        if path.exists():
            with closing(connect_readonly(path)) as connection:
                if metadata(connection) == proof:
                    return
        fd, temporary = tempfile.mkstemp(prefix=path.name + '.', dir=path.parent)
        os.close(fd)
        try:
            with closing(sqlite3.connect(temporary)) as connection:
                connection.execute('CREATE TABLE metadata (value TEXT NOT NULL)')
                connection.execute('CREATE TABLE frames (kind TEXT, family TEXT, payload BLOB, sha256 TEXT, '
                                   'PRIMARY KEY(kind, family))')
                connection.execute('CREATE TABLE prefixes (family TEXT PRIMARY KEY, prefix TEXT)')
                for kind, source in proof['sources'].items():
                    try:
                        frame = pandas.read_csv(source['path'], sep='\t', low_memory=False)
                    except pandas.errors.EmptyDataError:
                        frame = pandas.DataFrame()
                    if not frame.empty and 'orthogroup' not in frame:
                        raise ValueError(f'Candidate table has no orthogroup column: {source["path"]}')
                    if 'orthogroup' in frame:
                        for family, subset in frame.groupby(frame['orthogroup'].astype(str), sort=False):
                            data = encode_frame(subset)
                            connection.execute('INSERT INTO frames VALUES (?, ?, ?, ?)',
                                               (kind, family, data, hashlib.sha256(data).hexdigest()))
                    # Store the schema for families absent from one candidate table.
                    data = encode_frame(frame.iloc[:0])
                    connection.execute('INSERT INTO frames VALUES (?, ?, ?, ?)',
                                       (kind, '', data, hashlib.sha256(data).hexdigest()))
                    del frame
                connection.executemany('INSERT INTO prefixes VALUES (?, ?)', majority_prefixes(gene).items())
                validate_sources(proof, branch, gene)
                connection.execute('INSERT INTO metadata VALUES (?)', (json.dumps(proof, sort_keys=True),))
                connection.commit()
            os.replace(temporary, path)
        finally:
            Path(temporary).unlink(missing_ok=True)


def read_family(path, branch, gene, family):
    with closing(connect_readonly(path)) as connection:
        validate_sources(metadata(connection), branch, gene)
        frames = []
        for kind in ('branch', 'gene'):
            row = connection.execute('SELECT payload, sha256 FROM frames WHERE kind=? AND family=?',
                                     (kind, family)).fetchone()
            if row is None:
                row = connection.execute('SELECT payload, sha256 FROM frames WHERE kind=? AND family=?',
                                         (kind, '')).fetchone()
            data, expected = row
            if hashlib.sha256(data).hexdigest() != expected:
                raise ValueError('Corrupt private HGT candidate index payload')
            frames.append(decode_frame(data))
        prefix = connection.execute('SELECT prefix FROM prefixes WHERE family=?', (family,)).fetchone()
        return *frames, prefix[0] if prefix else ''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--database', required=True)
    parser.add_argument('--branch-tsv', required=True)
    parser.add_argument('--gene-tsv', required=True)
    args = parser.parse_args()
    prepare(args.database, args.branch_tsv, args.gene_tsv)


if __name__ == '__main__':
    main()
