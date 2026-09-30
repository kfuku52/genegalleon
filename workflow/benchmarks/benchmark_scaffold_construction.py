#!/usr/bin/env python3
"""Measure complete host-scaffold classification and locus aggregation."""
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
    import pandas as pd
    import scaffold_taxonomy as scaffold

    class LocalRanks:
        def ranks(self, taxid):
            if taxid not in (3, 5):
                return {}
            return {rank: (1 if i == 0 else taxid * 10 + i)
                    for i, rank in enumerate(scaffold.RANKS)}

    rows, taxa, loci = [], [], {}
    for i in range(args.genes):
        gene, transcript = f'g{i:05d}', f't{i:05d}'
        rows.append({'gene_id': gene, 'chromosome': '' if i % 131 == 0 else f's{i // 64}',
                     'gff_transcript_id': transcript,
                     'splice_mode': 'trans-splicing' if i % 211 == 0 else 'cis-splicing'})
        if i % 17:
            taxa.append({'gene_id': gene, 'lca_taxid': (0, 3, 5)[(i // 7) % 3]})
        if (i // 2) % 4 < 3:
            loci[transcript] = f'L{i // 2}'
    gff, taxonomy = pd.DataFrame(rows), pd.DataFrame(taxa)
    original = (gff.copy(deep=True), taxonomy.copy(deep=True))
    started = time.perf_counter()
    genes, summaries = scaffold.build_tables(gff, taxonomy, 'Host_alpha', 3, LocalRanks(), loci)
    seconds = time.perf_counter() - started
    pd.testing.assert_frame_equal(gff, original[0])
    pd.testing.assert_frame_equal(taxonomy, original[1])
    digest = hashlib.sha256()
    for frame in (genes, summaries):
        digest.update(json.dumps({'columns': frame.columns.tolist(), 'dtypes': [str(x) for x in frame.dtypes],
                                  'index_name': frame.index.name}, sort_keys=True).encode())
        digest.update(frame.to_csv(sep='\t', index=True).encode())
    print(json.dumps({'seconds': seconds, 'output_sha256': digest.hexdigest(),
                      'gene_rank_rows': len(genes), 'scaffold_rank_rows': len(summaries),
                      'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--genes', type=int, default=2048)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 64 <= args.genes <= 8192 or args.genes % 64:
        parser.error('--genes must be a multiple of 64 in 64..8192')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'genes': args.genes, 'python': sys.version,
               'platform': platform.platform(), 'samples': []}
    for trial in range(4):
        raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()), '--worker',
            '--support-root', str(args.support_root), '--genes', str(args.genes)], text=True)
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
