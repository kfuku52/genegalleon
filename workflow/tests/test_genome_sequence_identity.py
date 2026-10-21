"""Exact sequence identity and rejection contracts across streaming boundaries."""
import gzip
import hashlib
import importlib
import json
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
identity = importlib.import_module('genome_sequence_identity')


def expected(sequences):
    entries = sorted((name, (len(dna), hashlib.sha256(dna.encode()).hexdigest()))
                     for name, dna in sequences.items())
    return hashlib.sha256(json.dumps(entries, separators=(',', ':')).encode()).hexdigest()


@pytest.mark.parametrize('compressed', [False, True])
def test_identity_preserves_names_wrapping_case_unicode_and_contig_order(tmp_path, compressed):
    path = tmp_path / 'genome.fa'  # Gzip remains signature-based.
    text = '>chrβ first description\n \tacgt\r\n\u00a0rYSwKmbdhvnßſ \n>chrA\n tt\n'
    path.write_bytes(gzip.compress(text.encode()) if compressed else text.encode())
    digest = expected({'chrβ': 'ACGTRYSWKMBDHVNSSS', 'chrA': 'TT'})
    assert identity.genome_sequence_digest(path) == digest
    path.write_text('>chrA changed description\nTT\n>chrβ\nACGTRYSWKMBDHVNSSS\n')
    assert identity.genome_sequence_digest(path) == digest


@pytest.mark.parametrize('text,message', [
    ('  \n\t\n', 'no contigs'),
    ('ACGT\n', 'precedes'),
    ('>a\nACGT\n>a\nACGT\n', 'duplicate'),
    ('>a\nACGT\n> \nTT\n', 'identifier'),
    ('>a\n\n>b\nAT\n', 'Empty'),
    ('>a\nAT\n>b\n', 'Empty'),
    ('>a\nAC GT\n', 'Invalid'),
    ('>a\nAC\x00GT\n', 'Invalid'),
])
def test_invalid_inputs_still_fail(tmp_path, text, message):
    path = tmp_path / 'invalid.fa'
    path.write_text(text)
    with pytest.raises(ValueError, match=message):
        identity.genome_sequence_digest(path)


@pytest.mark.parametrize('position', [0, 65535, 65536, 65537, 1024**2 - 1, 1024**2, 1024**2 + 1, None])
def test_large_contig_validates_every_base_and_exact_length(tmp_path, position):
    dna = 'ACGTN' * 300000
    if position is not None:
        dna = dna[:position] + 'X' + dna[position + 1:]
    path = tmp_path / 'large.fa'
    path.write_text('>long\n' + '\n'.join(dna[i:i + 60] for i in range(0, len(dna), 60)) + '\n')
    if position is None:
        assert identity.genome_sequence_digest(path) == expected({'long': dna})
    else:
        with pytest.raises(ValueError, match='Invalid'):
            identity.genome_sequence_digest(path)


@pytest.mark.parametrize('damage', ['truncated', 'crc'])
def test_compressed_stream_is_read_through_its_integrity_trailer(tmp_path, damage):
    path = tmp_path / 'broken.fa'
    data = gzip.compress(b'>a\n' + b'ACGT\n' * 300000)
    if damage == 'truncated':
        data = data[:-4]
    else:
        data = data[:-8] + bytes([data[-8] ^ 1]) + data[-7:]
    path.write_bytes(data)
    with pytest.raises((OSError, EOFError)):
        identity.genome_sequence_digest(path)


@pytest.mark.parametrize('change', ['append', 'replace', 'restore_mtime'])
def test_input_change_at_end_of_reading_cannot_publish_an_identity(tmp_path, monkeypatch, change):
    path = tmp_path / 'genome.fa'
    path.write_text('>a\nACGT\n')
    before = path.stat()
    real_open = open
    class ChangingReader:
        def __enter__(self):
            self.stream = real_open(path, 'rt', encoding='utf-8')
            return self
        def __exit__(self, *args):
            self.stream.close()
        def fileno(self):
            return self.stream.fileno()
        def __iter__(self):
            yield from self.stream
            if change == 'replace':
                other = path.with_name('replacement.fa')
                other.write_bytes(path.read_bytes())
                other.replace(path)
            elif change == 'restore_mtime':
                path.write_text('>a\nTGCA\n')
                os.utime(path, ns=(before.st_atime_ns, before.st_mtime_ns))
            else:
                with real_open(path, 'a') as output:
                    output.write('ACGT\n')
    monkeypatch.setattr(identity, 'open', lambda *a, **kw: ChangingReader(), raising=False)
    with pytest.raises(ValueError, match='changed'):
        identity.genome_sequence_digest(path)
