#!/usr/bin/env python3
"""Measure sequence/coordinate compatibility reporting and exact coordinate clearing."""
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
    import gff2genestat as gff
    import pandas as pd

    rows, records = [], []
    for i in range(args.genes):
        gene = f'Host_species_g{i:05d}'
        sequence = 'ATG' * 16
        if args.scenario == 'mixed' and i % 7 == 0:
            sequence += 'A'
        elif args.scenario == 'mixed' and i % 7 == 1:
            sequence += 'NN'
        elif args.scenario == 'mixed' and i % 7 == 2:
            sequence = sequence[:-2]
        rows.append(dict(gene_id=gene, feature_size=48, num_intron=1, cds_first_phase=0, start=i * 100 + 1,
                         end=i * 100 + 51, intron_positions='24', feature_blocks='1-24;28-51', utr_blocks='',
                         chromosome='chr1', strand='+', feature_block_sequences='chr1;chr1',
                         feature_block_strands='+;+', transcript_junction_positions='',
                         feature_type='CDS', structure_status='unchecked', unused_value=i))
        records.append((gene, gene, sequence))
    traits = pd.DataFrame(rows, index=[i * 3 + 7 for i in range(args.genes)])
    traits.index.name = 'source_row'
    source = tuple(records)
    started = time.perf_counter()
    gff.mark_incompatible_structures(traits, records)
    seconds = time.perf_counter() - started
    assert tuple(records) == source
    digest = hashlib.sha256(json.dumps({'columns': traits.columns.tolist(),
        'dtypes': [str(x) for x in traits.dtypes], 'index_name': traits.index.name}, sort_keys=True).encode())
    digest.update(traits.to_csv(sep='\t', index=True).encode())
    print(json.dumps({'seconds': seconds, 'output_sha256': digest.hexdigest(),
                      'compatible': int(traits.structure_status.eq('length_compatible').sum()),
                      'mismatches': int(traits.structure_status.eq('cds_length_mismatch').sum()),
                      'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--genes', type=int, default=4096)
    parser.add_argument('--scenario', choices=['mixed', 'compatible'], default='mixed')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 16 <= args.genes <= 8192:
        parser.error('--genes must be 16..8192')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'genes': args.genes, 'scenario': args.scenario, 'python': sys.version,
               'platform': platform.platform(), 'samples': []}
    for trial in range(4):
        raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()), '--worker',
            '--support-root', str(args.support_root), '--genes', str(args.genes), '--scenario', args.scenario], text=True)
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
