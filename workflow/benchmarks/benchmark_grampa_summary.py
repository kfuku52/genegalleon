#!/usr/bin/env python3
"""Measure the GRAMPA parsing CLI with repeated maps and complete output fingerprints."""
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


def create_inputs(directory, args):
    species = [f'Genus_species{i:03}' for i in range(args.species)]
    trees = []
    for i in range(args.trees):
        genes = [f'{name}_gene{i}_{copy}' if (i + j) % 2 else f'gene{i}_{copy}_{name}'
                 for j, name in enumerate(species) for copy in range(1 + i % 2)]
        trees.append('(' + ','.join(genes[::-1]) + ');')
    det_rows = []
    # Interleave maps and include a detailed record without an input tree.
    for mul_tree in (2, 1):
        for i in range(args.trees + 1, 0, -1):
            if args.format == 'modern':
                det_rows.append(f'{mul_tree}\t{i}\t2\t3\t5\t1\n')
            else:
                det_rows.append(f'* GT-{i} to MT-{mul_tree}\t2\t3\t5\t1\n')
    det_header = ('mul.tree\tgene.tree\tdups\tlosses\ttotal.score\tmaps\n' if args.format == 'modern'
                  else '# GT/MT combo\tdups\tlosses\tTotal score\tMaps\n')
    species_tree = '(' + ','.join(species[::-1]) + ');'
    if args.format == 'modern':
        out_text = 'mul.tree\th1.node\th2.node\tscore\tlabeled.tree\n'
        out_text += ''.join(f'{i}\tH1\tH2\t7\t{species_tree}\n' for i in (2, 1))
    else:
        out_text = ''.join(f'MT-{i}\tH1\tH2\t{species_tree}\t7\n' for i in (2, 1))
    texts = {'det.tsv': det_header + ''.join(det_rows), 'out.tsv': out_text,
             'trees.nwk': '\n'.join(trees) + '\n', 'species.nwk': species_tree,
             'names.tsv': ''.join(f'family{i:07}.nwk\n' for i in range(args.trees))}
    for name, text in texts.items():
        (directory / name).write_text(text)
    return hashlib.sha256(json.dumps(texts, sort_keys=True).encode()).hexdigest()


def worker(args):
    import pandas as pd

    with tempfile.TemporaryDirectory() as name:
        directory = Path(name)
        input_hash = create_inputs(directory, args)
        command = [sys.executable, str(args.support_root / 'parse_grampa.py'),
                   '--grampa_det', str(directory / 'det.tsv'), '--grampa_out', str(directory / 'out.tsv'),
                   '--gene_trees', str(directory / 'trees.nwk'), '--species_tree', str(directory / 'species.nwk'),
                   '--sorted_gene_tree_file_names', str(directory / 'names.tsv'), '--ncpu', str(args.cpus)]
        started = time.perf_counter()
        completed = subprocess.run(command, cwd=directory, text=True, capture_output=True, check=True)
        seconds = time.perf_counter() - started
        path = directory / 'grampa_summary.tsv'
        frame = pd.read_csv(path, sep='\t')
        stdout = '\n'.join(line for line in completed.stdout.splitlines()
                           if not line.startswith(('Starting parse_grampa.py:', 'Ending parse_grampa.py:')))
        digest = hashlib.sha256(path.read_bytes())
        digest.update(json.dumps({'dtypes': [str(dtype) for dtype in frame.dtypes], 'stdout': stdout,
                                  'stderr': completed.stderr, 'input_sha256': input_hash}, sort_keys=True).encode())
        assert frame.shape[0] == (args.trees + 1) * 2
        print(json.dumps({'seconds': seconds, 'output_sha256': digest.hexdigest(), 'rows': len(frame),
                          'peak_child_rss_kib_linux': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--trees', type=int, default=2048)
    parser.add_argument('--species', type=int, default=32)
    parser.add_argument('--format', choices=['legacy', 'modern'], default='modern')
    parser.add_argument('--cpus', type=int, choices=[1, 2], default=1)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 16 <= args.trees <= 4096 or not 4 <= args.species <= 64 or args.trees * args.species > 131072:
        parser.error('Require 16..4096 trees, 4..64 species and at most 131072 combinations')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    result = {'support_root': str(args.support_root), 'trees': args.trees, 'species': args.species,
              'format': args.format, 'cpus': args.cpus, 'python': sys.version,
              'platform': platform.platform(), 'samples': []}
    for trial in range(4):
        raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()), '--worker',
            '--support-root', str(args.support_root), '--trees', str(args.trees), '--species', str(args.species),
            '--format', args.format, '--cpus', str(args.cpus)], text=True)
        sample = json.loads(raw)
        if trial:
            result['samples'].append(sample)
        print('warmup' if trial == 0 else trial, sample['seconds'], flush=True)
    assert len({sample['output_sha256'] for sample in result['samples']}) == 1
    for metric in ('seconds', 'peak_child_rss_kib_linux'):
        result['median_' + metric] = statistics.median(sample[metric] for sample in result['samples'])
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('medians', {key: value for key, value in result.items() if key.startswith('median_')})


if __name__ == '__main__':
    main()
