#!/usr/bin/env python3
"""Measure FASTA record formatting and indexed extraction with exact output bytes."""
import argparse
import contextlib
import hashlib
import io
import json
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from types import SimpleNamespace


def worker(args):
    sys.path.insert(0, str(args.support_root))
    import fasta_sequence_store as store

    records = []
    for i in range(args.records):
        motif = 'ACGTNacgtn' if i % 2 else 'MKWVTFISLLFLFSSAYS'
        if i % 11 == 0:
            motif = 'ACG-TN*'
        elif i % 23 == 0:
            motif = 'AC GT\tNN'
        sequence = (motif * (args.length // len(motif) + 1))[:args.length]
        records.append((f'Genus_species_GENE{i:07} description {i}', sequence))
    output = io.StringIO()
    started = time.perf_counter()
    for header, sequence in records:
        store.write_record(output, header, sequence)
    format_seconds = time.perf_counter() - started
    format_hash = hashlib.sha256(output.getvalue().encode()).hexdigest()
    output.close()
    with tempfile.TemporaryDirectory() as name:
        directory = Path(name)
        source = directory / 'Genus_species.fa'
        source.write_text(''.join(f'>{header}\n{sequence}\n' for header, sequence in records))
        database, manifest = directory / 'sequences.sqlite', directory / 'manifest.json'
        # Use the real schema/build and default storage budgets; never bypass
        # integrity checks or construct an alternative test-only store.
        with contextlib.redirect_stdout(io.StringIO()):
            store.build_store(database, manifest, [(source, 'Genus_species')],
                              store.DEFAULT_MAX_DATABASE_BYTES, store.DEFAULT_MINIMUM_FREE_BYTES)
        assert store.database_content_current(database)
        database_before = hashlib.sha256(database.read_bytes()).hexdigest()
        source_before = hashlib.sha256(source.read_bytes()).hexdigest()
        patterns = directory / 'patterns.txt'
        patterns.write_text(''.join(f'gene{i:07}\n' for i in range(args.records - 1, -1, -1)))
        target = directory / 'result.fa'
        extract_args = SimpleNamespace(database=database, pattern_file=patterns, output=target,
            query_variants=True, ignore_case=True, prefix_species=True, require_all=True)
        stdout, stderr = io.StringIO(), io.StringIO()
        started = time.perf_counter()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            result = store.extract(extract_args)
        extract_seconds = time.perf_counter() - started
        assert result == 0 and stderr.getvalue() == ''
        assert hashlib.sha256(database.read_bytes()).hexdigest() == database_before
        assert hashlib.sha256(source.read_bytes()).hexdigest() == source_before
        assert store.database_content_current(database)
        digest = hashlib.sha256(target.read_bytes())
        digest.update(stdout.getvalue().replace(str(target), '$OUTPUT').encode())
        digest.update(stderr.getvalue().encode())
        print(json.dumps({'format_seconds': format_seconds, 'extract_seconds': extract_seconds,
            'format_sha256': format_hash, 'extract_sha256': digest.hexdigest(),
            'output_bytes': target.stat().st_size,
            'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--records', type=int, default=4096)
    parser.add_argument('--length', type=int, default=1536)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 16 <= args.records <= 16384 or not 60 <= args.length <= 8192 or args.records * args.length > 67108864:
        parser.error('Require 16..16384 records, 60..8192 characters each, at most 64 MiB of sequence')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    result = {'support_root': str(args.support_root), 'records': args.records, 'length': args.length,
              'python': sys.version, 'platform': platform.platform(), 'samples': []}
    for trial in range(4):
        raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()), '--worker',
            '--support-root', str(args.support_root), '--records', str(args.records), '--length', str(args.length)],
            text=True)
        sample = json.loads(raw)
        if trial:
            result['samples'].append(sample)
        print('warmup' if trial == 0 else trial, sample['format_seconds'], sample['extract_seconds'], flush=True)
    for metric in ('format_sha256', 'extract_sha256'):
        assert len({sample[metric] for sample in result['samples']}) == 1
    for metric in ('format_seconds', 'extract_seconds', 'peak_rss_kib_linux'):
        result['median_' + metric] = statistics.median(sample[metric] for sample in result['samples'])
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('medians', {key: value for key, value in result.items() if key.startswith('median_')})


if __name__ == '__main__':
    main()
