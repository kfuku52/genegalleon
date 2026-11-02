#!/usr/bin/env python3
"""Fresh-process comparisons of RSC summaries and OMArk top-clade parsing."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import pickle
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from functools import partial
from pathlib import Path

import numpy
import pandas

SUPPORT = Path(__file__).resolve().parents[1] / 'support'


def load(source):
    sys.path.insert(0, str(SUPPORT))
    spec = importlib.util.spec_from_file_location('summary_benchmark_source', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare(path, mode, rows, clades, clade_order):
    path.mkdir(parents=True, exist_ok=True)
    if mode != 'omark':
        index = numpy.arange(rows)
        frame = pandas.DataFrame({
            'analysis_id': numpy.char.add('analysis_', (index % 32).astype(str)),
            'model_id': numpy.char.add('model_', (index // 4).astype(str)),
            'inference_status': numpy.where(index % 7 == 0, 'failed', 'ok'),
            'optimizer_converged': numpy.where(index % 11 == 0, 'no', 'yes'),
            'response_evolution_parameter_status': numpy.where(index % 3 == 0, 'estimated', 'fixed'),
            'response_evolution_optimizer_converged': numpy.where(index % 13 == 0, 'false', 'true'),
            'predictor_evolution_parameter_status': numpy.where(index % 5 == 0, 'estimated', 'fixed'),
            'predictor_evolution_optimizer_converged': numpy.where(index % 17 == 0, '0', '1'),
            'predictor_type': numpy.where(index % 4 == 0, 'intercept', 'continuous'),
            'term': numpy.where(index % 4 == 0, '(Intercept)', 'size'),
            'source_term': numpy.where(index % 4 == 0, 'Intercept', 'size'),
            'n_species_events': (index % 100) + 1,
            'p_value': numpy.resize(numpy.array([0.0, 1.0, .2, .03, -.1, 1.1, numpy.inf, -numpy.inf, numpy.nan]), rows),
        })
        with (path / 'results.pickle').open('wb') as handle:
            pickle.dump(frame, handle, protocol=5)
        frame.to_csv(path / 'regression.tsv', sep='\t', index=False)
        pandas.DataFrame([{'tree_id': 'OG1', 'status': 'ok', 'n_models': rows}]).to_csv(
            path / 'status.tsv', sep='\t', index=False)
    else:
        for species in range(rows):
            lines = ['#The selected clade was Test', '#Number of conserved HOGs is: 1000',
                     'S:900,D:50[U:20,E:30],M:50', 'S:90.0%,D:5.0%[U:2.0%,E:3.0%],M:5.0%',
                     '#On the whole proteome, there are 20000 proteins',
                     "#Clade\tNCBI taxid\tNumber of associated proteins\tPercentage of proteome's total"]
            for clade in range(clades):
                proteins = (clade * 7919) % 500 if clade_order == 'mixed' else clades - clade
                lines.append(f'Clade_{clade % 97:03}\t{clade}\t{proteins}\t{clade / 10}%')
            species_dir = path / f'species_{species:05}'
            species_dir.mkdir(exist_ok=True)
            (species_dir / 'species.sum').write_text('\n'.join(lines) + '\n', encoding='utf-8')


def worker(args):
    module = load(args.source)
    if args.mode == 'status':
        with (args.fixtures / 'results.pickle').open('rb') as handle:
            frame = pickle.load(handle)
        operation = partial(module._status_from_results, 'OG1', frame, 32)
    elif args.mode == 'stat-tree':
        operation = partial(module.summarize_for_stat_tree, args.fixtures / 'regression.tsv', args.fixtures / 'status.tsv')
    else:
        files = sorted(args.fixtures.rglob('*.sum'))
        def operation():
            return pandas.DataFrame([module.parse_sum_file(path) for path in files], columns=module.HEADER)
    cpu = time.process_time()
    start = time.perf_counter()
    result = operation()
    seconds, cpu_seconds = time.perf_counter() - start, time.process_time() - cpu
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    status = Path('/proc/self/status')
    if status.is_file():
        peak = next(int(line.split()[1]) for line in status.read_text().splitlines() if line.startswith('VmHWM:'))
    with args.result.open('wb') as handle:
        pickle.dump(result, handle, protocol=5)
    print(json.dumps({'wall_s': seconds, 'cpu_s': cpu_seconds, 'peak_rss_kib': peak}))


def compare(args):
    filename = 'summarize_omark.py' if args.mode == 'omark' else 'reconciled_speciation_contrast.py'
    sources = {'baseline': args.baseline / filename}
    if not args.baseline_only:
        sources['current'] = SUPPORT / filename
    args.fixtures = args.fixtures.resolve()
    prepare(args.fixtures, args.mode, args.rows, args.clades, args.clade_order)
    samples = {name: [] for name in sources}
    expected = None
    with tempfile.TemporaryDirectory(prefix='summary-benchmark-') as directory:
        result = Path(directory) / 'result.pickle'
        for iteration in range(args.repeats + 1):
            names = list(sources)
            if iteration % 2:
                names.reverse()
            for name in names:
                command = [sys.executable, str(Path(__file__).resolve()), '--worker', '--source', str(sources[name]),
                           '--mode', args.mode, '--fixtures', str(args.fixtures), '--result', str(result)]
                completed = subprocess.run(command, check=True, text=True, capture_output=True)
                sample = json.loads(completed.stdout)
                with result.open('rb') as handle:
                    value = pickle.load(handle)
                if expected is None:
                    expected = value
                elif isinstance(expected, pandas.DataFrame):
                    pandas.testing.assert_frame_equal(expected, value, check_exact=True)
                else:
                    assert expected == value, (expected, value)
                if iteration:
                    samples[name].append(sample)
    report = {'mode': args.mode, 'rows': args.rows, 'clades': args.clades, 'clade_order': args.clade_order, 'repeats': args.repeats,
              'platform': platform.platform(), 'python': sys.version, 'pandas': pandas.__version__,
              'numpy': numpy.__version__, 'samples': samples, 'equivalent': True,
              'sources': {name: {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                          for name, path in sources.items()},
              'medians': {name: {key: statistics.median(sample[key] for sample in values) for key in values[0]}
                          for name, values in samples.items()}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report['medians'], indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--mode', choices=['status', 'stat-tree', 'omark'], default='status')
    parser.add_argument('--rows', type=int, default=200000)
    parser.add_argument('--clades', type=int, default=200)
    parser.add_argument('--clade-order', choices=['mixed', 'descending'], default='mixed')
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--fixtures', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--baseline-only', action='store_true')
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--source', type=Path)
    parser.add_argument('--result', type=Path)
    args = parser.parse_args()
    for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[key] = '1'
    if args.worker:
        worker(args)
    else:
        compare(args)


if __name__ == '__main__':
    main()
