#!/usr/bin/env python3
"""Measure raw expression replicate preparation, without fitting a model."""
import argparse
import csv
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
    import reconciled_speciation_contrast as adapter

    expression = args.worker / 'expression.tsv'
    metadata = args.worker / 'samples.tsv'
    responses = [f'response:{i}' for i in range(8)] + ['constant', 'missing_leaf']
    columns = [f'{response}_{replicate}' for response in responses for replicate in range(1, 4)]
    with expression.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['gene id', *columns])
        for gene in range(args.genes):
            values = []
            for response in range(len(responses)):
                for replicate in range(1, 4):
                    value = str((gene % 17) + response * 0.5 + replicate * 0.1)
                    if (gene + response + replicate) % 11 == 0 or (gene % 17 == 0 and replicate == 2):
                        value = 'NA'
                    if response == 8:
                        value = '7'
                    if response == 9 and gene == 0:
                        value = 'NA'
                    values.append(value)
            writer.writerow([f' Genus_sp{gene % 64}_g{gene:05d} ', *values])
    if args.case == 'paired':
        with metadata.open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['column', 'response', 'biological_id', 'technical_id', 'batch'])
            for column in reversed(columns):
                response, replicate = column.rsplit('_', 1)
                writer.writerow([column, f' {response} ', f'bio{replicate}', f'tech{replicate}', f'batch{replicate}'])
    command = ['prepare', '--expression', str(expression), '--species-traits', str(args.worker / 'unused.tsv'),
               '--expression-output', str(args.worker / 'expression.out.tsv'),
               '--species-traits-output', str(args.worker / 'traits.out.tsv'),
               '--analysis-plan-output', str(args.worker / 'plan.tsv'), '--metadata-output', str(args.worker / 'meta.tsv')]
    if args.case == 'paired':
        command += ['--sample-metadata', str(metadata)]
    parameters = adapter.build_parser().parse_args(command)
    started = time.perf_counter()
    frame, status = adapter._prepare_expression(parameters)
    seconds = time.perf_counter() - started
    assert status['status'] == 'ready'
    assert 'constant:constant_or_empty' in status['reason']
    assert 'missing_leaf:missing_leaf_values' in status['reason']
    payload = json.dumps({'columns': frame.columns.tolist(), 'dtypes': [str(value) for value in frame.dtypes],
                          'metadata': status}, sort_keys=True).encode() + frame.to_csv(sep='\t', index=False).encode()
    print(json.dumps({'seconds': seconds, 'output_sha256': hashlib.sha256(payload).hexdigest(),
                      'rows': len(frame), 'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--genes', type=int, default=3000)
    parser.add_argument('--case', choices=['unpaired', 'paired'], default='unpaired')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', type=Path)
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 32 <= args.genes <= 8192:
        parser.error('--genes must be 32..8192')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'genes': args.genes, 'case': args.case,
               'python': sys.version, 'platform': platform.platform(), 'samples': []}
    with tempfile.TemporaryDirectory(prefix='gg-expression-benchmark-') as temporary:
        for trial in range(4):
            root = Path(temporary) / str(trial)
            root.mkdir()
            raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()), '--worker', str(root),
                '--support-root', str(args.support_root), '--genes', str(args.genes), '--case', args.case], text=True)
            sample = json.loads(raw)
            if trial:
                results['samples'].append(sample)
            print('warmup' if trial == 0 else trial, sample['seconds'], flush=True)
    assert len({row['output_sha256'] for row in results['samples']}) == 1
    results['median_seconds'] = statistics.median(row['seconds'] for row in results['samples'])
    results['median_peak_rss_kib_linux'] = statistics.median(row['peak_rss_kib_linux'] for row in results['samples'])
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    print('median', results['median_seconds'], results['median_peak_rss_kib_linux'], flush=True)


if __name__ == '__main__':
    main()
