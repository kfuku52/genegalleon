#!/usr/bin/env python3
"""Measure CDS admission and two-source selection with exact decision records."""
import argparse
import hashlib
import json
import platform
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path


def worker(args):
    sys.path.insert(0, str(args.support_root))
    import cds_resolution as resolution

    cases = []
    codons = ('GCT', 'GAA', 'GGT', 'GCC', 'TTT', 'TAT', 'GTC')
    for i in range(args.genes):
        base = 'ATG' + ''.join(codons[(i + j) % len(codons)] for j in range(args.codons)) + 'TAA'
        supplied, derived, phase = base, base, 0
        if i % 8 == 1:
            supplied = ' \n' + base.lower() + '\t'
        elif i % 8 == 2:
            supplied = derived = base[:30] + 'NNN' + base[33:]
        elif i % 8 == 3:
            supplied = derived = base[:30] + 'TAA' + base[33:]
        elif i % 8 == 4:
            supplied = ''
        elif i % 8 == 5:
            derived = 'CCC' + base
        elif i % 8 == 6:
            supplied = derived = base[2:]
            phase = 1
        elif i % 8 == 7:
            derived = base[:30] + 'GCA' + base[33:]
        cases.append((supplied, derived, phase))
    digest = hashlib.sha256()
    sample = {}
    for name in ('unconstrained', 'phase_constrained', 'two_source_selection'):
        started = time.perf_counter()
        if name == 'unconstrained':
            results = [resolution.evaluate(supplied) for supplied, _derived, _phase in cases]
        elif name == 'phase_constrained':
            results = [resolution.evaluate(derived, phase=phase) for _supplied, derived, phase in cases]
        else:
            results = [resolution.select_candidate(supplied, derived, phase=phase)
                       for supplied, derived, phase in cases]
        sample[name + '_seconds'] = time.perf_counter() - started
        digest.update(name.encode())
        digest.update(json.dumps(results, sort_keys=True).encode())
        del results
    sample.update(output_sha256=digest.hexdigest(),
                  peak_rss_kib_linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    print(json.dumps(sample))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--genes', type=int, default=1024)
    parser.add_argument('--codons', type=int, default=500)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 32 <= args.genes <= 8192 or not 16 <= args.codons <= 2048:
        parser.error('Require 32..8192 genes and 16..2048 codons')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'genes': args.genes, 'codons': args.codons,
               'python': sys.version, 'platform': platform.platform(), 'samples': []}
    for trial in range(4):
        raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()), '--worker',
            '--support-root', str(args.support_root), '--genes', str(args.genes), '--codons', str(args.codons)], text=True)
        sample = json.loads(raw)
        if trial:
            results['samples'].append(sample)
        print('warmup' if trial == 0 else trial, {k: v for k, v in sample.items() if k.endswith('_seconds')}, flush=True)
    assert len({row['output_sha256'] for row in results['samples']}) == 1
    for metric in ('unconstrained_seconds', 'phase_constrained_seconds', 'two_source_selection_seconds', 'peak_rss_kib_linux'):
        results['median_' + metric] = statistics.median(row[metric] for row in results['samples'])
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    print('medians', {k: v for k, v in results.items() if k.startswith('median_')})


if __name__ == '__main__':
    main()
