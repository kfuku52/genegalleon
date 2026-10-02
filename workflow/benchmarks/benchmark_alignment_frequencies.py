#!/usr/bin/env python3
"""Measure alignment-derived F3X4 frequencies for a full tree and its two subroots."""
import argparse
import hashlib
import itertools
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
    import iqtree2mapnh as mapping

    codons = [''.join(codon) for codon in itertools.product('ACGT', repeat=3)]
    records, names = [], []
    for i in range(args.genes):
        body = []
        for j in range(args.codons):
            if j % 101 == 0:
                codon = 'NNN'
            elif j % 103 == 0:
                codon = 'A--'
            elif j % 107 == 0:
                codon = 'éAT'
            elif j % 5 == 0:
                codon = ('AAA', 'CCT', 'GGT')[i % 3]
            else:
                codon = codons[(i * 17 + j * 7) % len(codons)]
            body.append(codon)
        sequence = ''.join(body) + 'A' * (i % 3)
        if i % 3 == 0:
            sequence = sequence.replace('T', 'U')
        if i % 2:
            sequence = sequence.lower()
        records.append(f'>Host_species_g-{i:05d}\n{sequence}\n')
        names.append(f'Host_species_g_{i:05d}')
    with tempfile.TemporaryDirectory(prefix='gg-frequency-benchmark-') as temporary:
        path = Path(temporary) / 'alignment.fa'
        source = ''.join(records).encode()
        path.write_bytes(source)
        started = time.perf_counter()
        frequencies = [mapping.alignment_subset_nuc_freqs(path, 'F3X4+G4', subset)
                       for subset in (None, names[:args.genes // 2], names[args.genes // 2:])]
        seconds = time.perf_counter() - started
        thetas = [mapping.kfseq.nuc_freq2theta(nuc_freqs=values) for values in frequencies]
        assert path.read_bytes() == source
    digest = hashlib.sha256(json.dumps({'frequencies': frequencies, 'thetas': thetas},
                                       sort_keys=True).encode()).hexdigest()
    print(json.dumps({'seconds': seconds, 'output_sha256': digest,
                      'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--genes', type=int, default=1024)
    parser.add_argument('--codons', type=int, default=512)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 16 <= args.genes <= 2048 or not 16 <= args.codons <= 2048 or args.genes * args.codons > 1048576:
        parser.error('Require 16..2048 genes/codons, with at most 1048576 total gene/codon combinations')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'genes': args.genes, 'codons': args.codons,
               'python': sys.version, 'platform': platform.platform(), 'samples': []}
    for trial in range(4):
        raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()), '--worker',
            '--support-root', str(args.support_root), '--genes', str(args.genes),
            '--codons', str(args.codons)], text=True)
        sample = json.loads(raw)
        if trial:
            results['samples'].append(sample)
        print('warmup' if trial == 0 else trial, sample['seconds'], flush=True)
    assert len({row['output_sha256'] for row in results['samples']}) == 1
    for metric in ('seconds', 'peak_rss_kib_linux'):
        results['median_' + metric] = statistics.median(row[metric] for row in results['samples'])
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    print('medians', results['median_seconds'], results['median_peak_rss_kib_linux'])


if __name__ == '__main__':
    main()
