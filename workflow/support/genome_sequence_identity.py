"""Exact contig-name/DNA identity independent of FASTA wrapping and gzip."""
import gzip
import hashlib
import json
import os
from pathlib import Path

_DNA = str.maketrans('', '', 'ACGTRYSWKMBDHVN')

def genome_sequence_digest(path):
    path = Path(path)
    with path.open('rb') as stream:
        compressed = stream.read(2) == b'\x1f\x8b'
    entries, identifier, sequence, length = {}, None, None, 0
    opener = gzip.open if compressed else open
    fragments, fragment_length = [], 0
    def flush():
        nonlocal length, fragment_length
        if fragments:
            chunk = ''.join(fragments).upper()
            if chunk.translate(_DNA):
                raise ValueError('Invalid genomic DNA symbol')
            sequence.update(chunk.encode('ascii'))
            length += len(chunk)
            fragments.clear()
            fragment_length = 0
    def finish():
        flush()
        if identifier is not None:
            if length == 0:
                raise ValueError('Empty genome contig')
            entries[identifier] = (length, sequence.hexdigest())
    with opener(path, 'rt', encoding='utf-8') as stream:
        def stamp(info):
            return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns
        before = stamp(os.fstat(stream.fileno()))
        for line in stream:
            if line.startswith('>'):
                finish()
                tokens = line[1:].split()
                if not tokens or tokens[0] in entries:
                    raise ValueError('Missing or duplicate genome contig identifier')
                identifier, sequence, length = tokens[0], hashlib.sha256(), 0
            else:
                chunk = line.strip()
                if not chunk:
                    continue
                if identifier is None:
                    raise ValueError('Genome sequence precedes its contig identifier')
                fragments.append(chunk)
                fragment_length += len(chunk)
                if fragment_length >= 64 * 1024:
                    flush()
        flush()
        if before != stamp(os.fstat(stream.fileno())) or before != stamp(path.stat()):
            raise ValueError('Genome changed during sequence identity reading')
    finish()
    if not entries:
        raise ValueError('Genome has no contigs')
    return hashlib.sha256(json.dumps(sorted(entries.items()), separators=(',', ':')).encode()).hexdigest()
