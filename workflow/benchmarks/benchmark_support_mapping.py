#!/usr/bin/env python3
"""Measure full unrooted support mapping on balanced and deep gene trees."""
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


def make_tree(module, tips, shape, reverse):
    tree = module.ete4.PhyloTree()
    pending = [(tree, 0, tips)]
    while pending:
        node, start, count = pending.pop()
        if count == 1:
            node.name = f'gene_{start:05d}'
            continue
        size = count // 2 if shape == 'balanced' else 1
        ranges = [(start, size), (start + size, count - size)]
        if reverse:
            ranges.reverse()
        for offset, length in ranges:
            child = node.add_child(dist=1.0, support=90.0)
            pending.append((child, offset, length))
    for branch_id, node in enumerate(tree.traverse()):
        node.add_prop('branch_id', branch_id)
    return tree


def tree_signature(tree):
    return [(node.name, node.dist, node.support, node.props['branch_id'], len(node.children))
            for node in tree.traverse()]


def worker(args):
    sys.path.insert(0, str(args.support_root))
    import orthogroup_statistics as stats

    rooted = make_tree(stats, args.tips, args.shape, False)
    supported = make_tree(stats, args.tips, args.shape, True)
    signatures = (tree_signature(rooted), tree_signature(supported))
    started = time.perf_counter()
    result = stats.map_internal_support_by_split(rooted, supported, support_max=100, require_support=True)
    seconds = time.perf_counter() - started
    assert signatures == (tree_signature(rooted), tree_signature(supported))
    digest = hashlib.sha256(json.dumps(result, sort_keys=True).encode()).hexdigest()
    print(json.dumps({'seconds': seconds, 'output_sha256': digest,
                      'mapped_branches': len(result[0]), 'diagnostics': result[1],
                      'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--tips', type=int, default=2048)
    parser.add_argument('--shape', choices=['balanced', 'comb'], default='balanced')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', action='store_true')
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 32 <= args.tips <= (4096 if args.shape == 'balanced' else 1024):
        parser.error('--tips must be 32..4096 (balanced) or 32..1024 (comb)')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'tips': args.tips, 'shape': args.shape,
               'python': sys.version, 'platform': platform.platform(), 'samples': []}
    for trial in range(4):
        raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()), '--worker',
            '--support-root', str(args.support_root), '--tips', str(args.tips), '--shape', args.shape], text=True)
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
