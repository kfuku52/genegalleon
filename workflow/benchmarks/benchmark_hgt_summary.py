#!/usr/bin/env python3
"""Measure HGT candidate summaries and aggregation with local synthetic evidence."""
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
    import pandas as pd
    import score_hgt_candidates as scorer

    identifiers = [f'Genus_sp{i % 64}_g{i:05d}' for i in range(args.genes)]
    leaves = pd.DataFrame({
        'node_name': identifiers, 'taxon': ['Arabidopsis thaliana' if i % 2 else 'Oryza sativa' for i in range(args.genes)],
        'num_intron': [i % 4 for i in range(args.genes)], 'intron_is_imputed': [i % 11 == 0 for i in range(args.genes)],
        'expression_root': ['NA' if i % 5 == 0 else i % 7 for i in range(args.genes)],
        'expression_leaf': ['NA' if i % 7 == 0 else i % 3 for i in range(args.genes)],
        'synteny_support_score': [float('nan') if i % 13 == 0 else i % 4 / 4 for i in range(args.genes)],
        'sprot_best': [f'P{i:06d}' for i in range(args.genes)],
        'organism': ['Escherichia coli' if i % 3 else 'Bacillus subtilis' for i in range(args.genes)],
    })
    contamination = pd.DataFrame({'gene_id': identifiers[::3], 'is_compatible_lineage': False,
                                   'lca_taxid': '562', 'lca_sciname': 'Escherichia coli'})
    branches = [pd.Series({'orthogroup': 'OGfixture', 'branch_id': i + 1, 'node_name': f'Node{i}',
                           'gene_labels': '; '.join(ids + ['absent_gene']), 'generax_event': 'H',
                           'generax_transfer': 'Y', 'clade_min_expression_pearsoncor': 0.7})
                 for i, ids in enumerate((identifiers[::-1], identifiers[::2]))]
    resolver = scorer.TaxonomyResolver('')
    branch_records, gene_records = [], []
    started = time.perf_counter()
    for branch in branches:
        record, genes = scorer.summarize_candidate_branch(branch, leaves, contamination,
                                                          ['expression_root', 'expression_leaf'], resolver)
        branch_records.append(record)
        gene_records.extend(genes)
    summarize_seconds = time.perf_counter() - started
    branch_frame = pd.DataFrame(branch_records, columns=scorer.BRANCH_OUTPUT_COLUMNS)
    raw_genes = pd.DataFrame(gene_records)
    started = time.perf_counter()
    gene_frame = scorer.aggregate_gene_records(raw_genes)
    orthogroup_frame = scorer.aggregate_orthogroup_records(branch_frame, gene_frame)
    aggregate_seconds = time.perf_counter() - started
    digest = hashlib.sha256()
    for frame in (branch_frame, raw_genes, gene_frame, orthogroup_frame):
        digest.update(json.dumps({'columns': frame.columns.tolist(), 'dtypes': [str(x) for x in frame.dtypes]}, sort_keys=True).encode())
        digest.update(frame.to_csv(sep='\t', index=False).encode())
    print(json.dumps({'summarize_seconds': summarize_seconds, 'aggregate_seconds': aggregate_seconds,
                      'output_sha256': digest.hexdigest(), 'genes': len(gene_frame),
                      'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--genes', type=int, default=2048)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', type=Path)
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 32 <= args.genes <= 4096:
        parser.error('--genes must be 32..4096')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'genes': args.genes, 'python': sys.version,
               'platform': platform.platform(), 'samples': []}
    with tempfile.TemporaryDirectory(prefix='gg-hgt-summary-benchmark-') as temporary:
        for trial in range(4):
            root = Path(temporary) / str(trial)
            root.mkdir()
            raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()), '--worker', str(root),
                '--support-root', str(args.support_root), '--genes', str(args.genes)], text=True)
            sample = json.loads(raw)
            if trial:
                results['samples'].append(sample)
            print('warmup' if trial == 0 else trial, sample['summarize_seconds'], sample['aggregate_seconds'], flush=True)
    assert len({row['output_sha256'] for row in results['samples']}) == 1
    for metric in ('summarize_seconds', 'aggregate_seconds', 'peak_rss_kib_linux'):
        results['median_' + metric] = statistics.median(row[metric] for row in results['samples'])
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    print('medians', results['median_summarize_seconds'], results['median_aggregate_seconds'], results['median_peak_rss_kib_linux'])


if __name__ == '__main__':
    main()
