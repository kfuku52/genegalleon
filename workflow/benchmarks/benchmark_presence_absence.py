#!/usr/bin/env python3
"""Measure complete family/species long-table construction with incomplete families."""
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
    import gene_family_presence_absence as presence_absence
    import numpy as np
    import pandas as pd

    species = [f'Genus_species{i:03d}' for i in range(args.species)]
    families = [f'OG{i:07d}' for i in range(args.families)]
    columns, statuses = {}, {}
    for i, family in enumerate(families):
        incomplete = i % 7 == 0
        columns[family] = [np.nan] * args.species if incomplete else [(i + j) % 5 for j in range(args.species)]
        statuses[family] = 'missing_stat_branch' if incomplete else 'complete'
    counts = pd.DataFrame(columns, index=species)
    presence = (counts > 0).astype('Int64').where(~counts.isna(), pd.NA)
    originals = (counts.copy(deep=True), presence.copy(deep=True))
    # Exercise species/family ordering independently of the matrix layout.
    selected_species = species[::-1]
    family_order = families[::-1]
    started = time.perf_counter()
    result = presence_absence.build_long_table(counts, presence, selected_species, family_order,
                                               statuses, 'orthogroup')
    seconds = time.perf_counter() - started
    pd.testing.assert_frame_equal(counts, originals[0])
    pd.testing.assert_frame_equal(presence, originals[1])
    digest = hashlib.sha256(json.dumps({'columns': result.columns.tolist(),
        'dtypes': [str(x) for x in result.dtypes], 'index_name': result.index.name}, sort_keys=True).encode())
    digest.update(result.to_csv(sep='\t', index=True).encode())
    print(json.dumps({'seconds': seconds, 'output_sha256': digest.hexdigest(), 'rows': len(result),
                      'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--families', type=int, default=2048)
    parser.add_argument('--species', type=int, default=64)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 16 <= args.families <= 8192 or not 8 <= args.species <= 256:
        parser.error('--families must be 16..8192 and --species must be 8..256')
    if args.families * args.species > 524288:
        parser.error('at most 524288 family/species combinations are allowed')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'families': args.families, 'species': args.species,
               'python': sys.version, 'platform': platform.platform(), 'samples': []}
    for trial in range(4):
        raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()), '--worker',
            '--support-root', str(args.support_root), '--families', str(args.families),
            '--species', str(args.species)], text=True)
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
