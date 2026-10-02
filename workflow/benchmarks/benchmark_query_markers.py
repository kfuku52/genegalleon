#!/usr/bin/env python3
"""Measure direct query-ID matching and complete branch marker annotation."""
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


def worker(args):
    sys.path.insert(0, str(args.support_root))
    import annotate_stat_branch_query_markers as markers
    import pandas as pd

    query, fasta, blast = [args.worker / name for name in ('queries.txt', 'query.fa', 'blast.tsv')]
    branch, output = args.worker / 'branches.tsv', args.worker / 'annotated.tsv'
    queries = [f'q{i:05d}' for i in range(args.queries)]
    query.write_text('\n'.join(queries + [f'prefix-{value}' for value in queries[::8]] + ['NA', '-x', 'x']) + '\n')
    tips = []
    for i in range(args.tips):
        value = queries[i % len(queries)]
        tips.append(value if i % 19 == 0 else f'Sp{i}_prefix' + ('−' if i % 7 == 0 else '-') + value)
    tips += ['NA', 'sp.-x', 'sp_x', 'unmatched']
    fasta.write_text(''.join(f'>{tip}\nMAAA\n' for tip in tips[::17]))
    rows = []
    for i, tip in enumerate(tips):
        if i % 8 == 0:
            rows.append({'branch_id': f'internal{i}', 'node_name': f'Node{i}', 'so_event': 'S', 'query_marker': 'stale'})
        rows.append({'branch_id': str(i), 'node_name': tip, 'so_event': 'L', 'query_marker': 'stale'})
    pd.DataFrame(rows).to_csv(branch, sep='\t', index=False)
    hits = []
    for i in range(min(128, args.queries)):
        for target, coverage in ((tips[i], 0.9), (tips[i + 1], 0.8), ('outside_tree', 1.0)):
            hits.append({'qacc': f'external{i}', 'sacc': target, 'qjointcov': coverage,
                         'evalue': 'invalid;1e-9;1e-12', 'bitscore': '100;150'})
    pd.DataFrame(hits).to_csv(blast, sep='\t', index=False)
    started = time.perf_counter()
    direct = markers.direct_query_sources_by_node(tips, query, fasta)
    direct_seconds = time.perf_counter() - started
    captured = io.StringIO()
    started = time.perf_counter()
    with contextlib.redirect_stdout(captured):
        markers.annotate_stat_branch(branch, query, blast, output, fasta)
    annotation_seconds = time.perf_counter() - started
    payload = json.dumps(direct, sort_keys=True).encode() + output.read_bytes()
    payload += captured.getvalue().replace(str(output), 'OUTPUT').encode()
    print(json.dumps({'direct_seconds': direct_seconds, 'annotation_seconds': annotation_seconds,
                      'output_sha256': hashlib.sha256(payload).hexdigest(), 'rows': len(rows),
                      'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--tips', type=int, default=4096)
    parser.add_argument('--queries', type=int, default=512)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', type=Path)
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 256 <= args.tips <= 8192 or not 32 <= args.queries <= 2048:
        parser.error('--tips must be 256..8192; --queries must be 32..2048')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'tips': args.tips, 'queries': args.queries,
               'python': sys.version, 'platform': platform.platform(), 'samples': []}
    with tempfile.TemporaryDirectory(prefix='gg-query-marker-benchmark-') as temporary:
        for trial in range(4):
            root = Path(temporary) / str(trial)
            root.mkdir()
            raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()), '--worker', str(root),
                '--support-root', str(args.support_root), '--tips', str(args.tips), '--queries', str(args.queries)], text=True)
            sample = json.loads(raw)
            if trial:
                results['samples'].append(sample)
            print('warmup' if trial == 0 else trial, sample['direct_seconds'], sample['annotation_seconds'], flush=True)
    assert len({row['output_sha256'] for row in results['samples']}) == 1
    for metric in ('direct_seconds', 'annotation_seconds', 'peak_rss_kib_linux'):
        results['median_' + metric] = statistics.median(row[metric] for row in results['samples'])
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    print('medians', results['median_direct_seconds'], results['median_annotation_seconds'], results['median_peak_rss_kib_linux'])


if __name__ == '__main__':
    main()
