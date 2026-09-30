#!/usr/bin/env python3
"""Measure complete native intron-ASR table validation and branch-ID translation."""
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

from benchmark_support_mapping import make_tree


def worker(args):
    sys.path.insert(0, str(args.support_root))
    import orthogroup_statistics as stats
    import pandas as pd

    tree = make_tree(stats, args.tips, 'balanced', False)
    leaves = list(tree.leaves())
    for leaf, name in zip(leaves[:3], ['NA', 'None', '001'], strict=True):
        leaf.name = name
    nodes = list(tree.traverse(strategy='levelorder'))
    native_ids = {node: index for index, node in enumerate(nodes)}
    rows = []
    for index, node in enumerate(nodes):
        leaf, root = stats.node_is_leaf(node), stats.node_is_root(node)
        if not leaf:
            node.name = f'n{index}'
        missing = not leaf or index % 11 == 0
        count = 'NA' if missing else index % 5
        present = (0.4 + (index % 5) / 20) if missing else float(count > 0)
        rows.append({'branch_id': index, 'parent': -1 if root else native_ids[node.up],
                     'node_class': 'root' if root else 'leaf' if leaf else 'intnode', 'name': node.name,
                     'num_intron': count, 'is_imputed': leaf and missing,
                     'p_intron_present': present, 'p_intron_absent': 1 - present,
                     'additional_text': 'NA', 'additional_number': 0.125})
    tree_path, table_path = args.worker / 'dated.nwk', args.worker / 'asr.tsv'
    tree_path.write_text(tree.write(parser=1, format_root_node=True))
    pd.DataFrame(rows[::-1]).to_csv(table_path, sep='\t', index=False)
    inputs = (tree_path.read_bytes(), table_path.read_bytes())
    started = time.perf_counter()
    result = stats.load_asr_intron_branch_table(str(table_path), str(tree_path))
    seconds = time.perf_counter() - started
    assert inputs == (tree_path.read_bytes(), table_path.read_bytes())
    digest = hashlib.sha256(json.dumps({'columns': result.columns.tolist(),
        'dtypes': [str(x) for x in result.dtypes], 'index_name': result.index.name}, sort_keys=True).encode())
    digest.update(result.to_csv(sep='\t', index=True).encode())
    print(json.dumps({'seconds': seconds, 'output_sha256': digest.hexdigest(), 'rows': len(result),
                      'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--tips', type=int, default=2048)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', type=Path)
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 32 <= args.tips <= 8192:
        parser.error('--tips must be 32..8192')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'tips': args.tips, 'python': sys.version,
               'platform': platform.platform(), 'samples': []}
    with tempfile.TemporaryDirectory(prefix='gg-asr-loading-benchmark-') as temporary:
        for trial in range(4):
            root = Path(temporary) / str(trial)
            root.mkdir()
            raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()), '--worker', str(root),
                '--support-root', str(args.support_root), '--tips', str(args.tips)], text=True)
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
