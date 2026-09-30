#!/usr/bin/env python3
"""Measure exact seeded ortholog-decay calculations and summaries."""
import argparse
import hashlib
import json
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def worker(args):
    sys.path.insert(0, str(args.support_root))
    import numpy as np
    import single_copy_ortholog_decay_plot as decay

    counts = np.random.default_rng(811).integers(0, 4, size=(args.orthogroups, args.species), dtype=np.int64)
    counts[0, :] = 1
    counts[1, :] = 0
    counts[2::31, :] = 2
    selected = counts[::2, :].copy()
    sizes = sorted({1, 3, 8, args.species // 2, args.species})
    seconds, digest = {}, hashlib.sha256()
    for name, matrix in [('all', None), ('selected', selected)]:
        started = time.perf_counter()
        values = decay.calculate_decay(counts, sizes, args.replicates, seed=73, selected_counts=matrix)
        summary = decay.summarize_decay(sizes, values)
        seconds[name] = time.perf_counter() - started
        digest.update(values.tobytes())
        digest.update(json.dumps({'columns': summary.columns.tolist(), 'dtypes': [str(x) for x in summary.dtypes]}, sort_keys=True).encode())
        digest.update(summary.to_csv(sep='\t', index=False).encode())
    print(json.dumps({'seconds': seconds, 'output_sha256': digest.hexdigest(),
                      'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--orthogroups', type=int, default=8192)
    parser.add_argument('--species', type=int, default=64)
    parser.add_argument('--replicates', type=int, default=128)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', type=Path)
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 128 <= args.orthogroups <= 32768 or not 16 <= args.species <= 256 or not 1 <= args.replicates <= 512:
        parser.error('--orthogroups must be 128..32768; --species 16..256; --replicates 1..512')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'orthogroups': args.orthogroups, 'species': args.species,
               'replicates': args.replicates, 'python': sys.version, 'platform': platform.platform(), 'samples': []}
    with tempfile.TemporaryDirectory(prefix='gg-decay-benchmark-') as temporary:
        for trial in range(4):
            root = Path(temporary) / str(trial)
            root.mkdir()
            raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()), '--worker', str(root),
                '--support-root', str(args.support_root), '--orthogroups', str(args.orthogroups),
                '--species', str(args.species), '--replicates', str(args.replicates)], text=True)
            sample = json.loads(raw)
            if trial:
                results['samples'].append(sample)
            print('warmup' if trial == 0 else trial, sample['seconds'], flush=True)
    assert len({row['output_sha256'] for row in results['samples']}) == 1
    results['median_seconds'] = {name: statistics.median(row['seconds'][name] for row in results['samples'])
                                 for name in ('all', 'selected')}
    results['median_peak_rss_kib_linux'] = statistics.median(row['peak_rss_kib_linux'] for row in results['samples'])
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    print('medians', results['median_seconds'], results['median_peak_rss_kib_linux'])


if __name__ == '__main__':
    main()
