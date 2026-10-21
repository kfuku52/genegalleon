#!/usr/bin/env python3
"""Compare exact genome identities in alternating fresh-process trials.

Synthetic warm-cache FASTA parsing only; no prediction, GPU, NAS or accuracy claim.
"""
import argparse
import contextlib
import gzip
import hashlib
import importlib.util
import json
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def fixture(path, bases, contigs, compressed):
    line = (b'ACGTNacgtRYSWKMBDHVN' * 4)[:60]
    block = (line + b'\n') * 16384
    with contextlib.ExitStack() as stack:
        stream = stack.enter_context(path.open('wb'))
        if compressed:
            stream = stack.enter_context(gzip.GzipFile(filename='', fileobj=stream, mode='wb', mtime=0))
        for index in range(contigs):
            stream.write(f'>chr{index:06} synthetic control\n'.encode())
            remaining = bases // contigs + (index < bases % contigs)
            while remaining >= len(line) * 16384:
                stream.write(block)
                remaining -= len(line) * 16384
            full, tail = divmod(remaining, len(line))
            stream.write((line + b'\n') * full)
            if tail:
                stream.write(line[:tail] + b'\n')


def worker(args):
    spec = importlib.util.spec_from_file_location('measured_identity', args.source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    started, cpu = time.perf_counter(), time.process_time()
    digest = module.genome_sequence_digest(args.input)
    return dict(seconds=time.perf_counter() - started, cpu_seconds=time.process_time() - cpu,
                peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss *
                    (1 if sys.platform == 'darwin' else 1024), output_sha256=digest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[1] / 'support/genome_sequence_identity.py')
    parser.add_argument('--baseline-source', type=Path)
    parser.add_argument('--mib', type=int, default=128)
    parser.add_argument('--trials', type=int, default=3)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--input', type=Path)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args)))
        return
    if not args.output or not 8 <= args.mib <= 2048 or not 3 <= args.trials <= 10:
        parser.error('--output, 8..2048 MiB and 3..10 trials are required')
    if args.output.exists():
        parser.error('output already exists')
    sources = {'current': args.source.resolve()}
    if args.baseline_source:
        sources = {'baseline': args.baseline_source.resolve(), **sources}
    result = dict(python=sys.version, platform=platform.platform(), mib=args.mib,
                  sources={name: {'path': str(path), 'sha256': sha(path)} for name, path in sources.items()}, cases={})
    with tempfile.TemporaryDirectory(prefix='gg-identity-benchmark-') as temporary:
        for name, bases, contigs, compressed in (
                ('wrapped', args.mib * 1024**2, 8, False),
                ('gzip', args.mib * 1024**2, 8, True),
                ('fragmented', 8 * 1024**2, 4096, False)):
            path = Path(temporary) / (name + '.fa' + ('.gz' if compressed else ''))
            fixture(path, bases, contigs, compressed)
            case = dict(bases=bases, contigs=contigs, columns=60, input_bytes=path.stat().st_size,
                        input_sha256=sha(path), samples={label: [] for label in sources})
            digests = set()
            for trial in range(args.trials + 1):
                labels = list(sources)
                if trial % 2:
                    labels.reverse()
                for label in labels:
                    command = [sys.executable, '-B', str(Path(__file__).resolve()), '--worker',
                               '--source', str(sources[label]), '--input', str(path)]
                    sample = json.loads(subprocess.check_output(command, text=True))
                    digests.add(sample['output_sha256'])
                    if trial:
                        case['samples'][label].append(sample)
                print(name, 'warmup' if not trial else trial, flush=True)
            if len(digests) != 1:
                raise RuntimeError('Genome identities differ between sources or trials')
            case['output_sha256'] = digests.pop()
            case['medians'] = {label: {metric: statistics.median(row[metric] for row in rows)
                                     for metric in ('seconds', 'cpu_seconds', 'peak_rss_bytes')}
                               for label, rows in case['samples'].items()}
            result['cases'][name] = case
    if result['cases']['wrapped']['output_sha256'] != result['cases']['gzip']['output_sha256']:
        raise RuntimeError('Compression changed the genome identity')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({name: case['medians'] for name, case in result['cases'].items()}, indent=2))


if __name__ == '__main__':
    main()
