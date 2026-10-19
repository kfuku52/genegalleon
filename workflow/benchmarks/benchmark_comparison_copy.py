#!/usr/bin/env python3
"""Compare isolated cache-copy implementations on identical receipted files.

Use the runtime wrapper. Supply a historical rescue_gene_models.py as the
baseline; both functions execute with the same installed runtime and globals.
Only newly created scratch directories and the requested report are written.
"""
import argparse
import ast
import hashlib
import inspect
import json
import os
import platform
import resource
import stat
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

SUPPORT = Path(__file__).resolve().parents[1] / 'support'


def io_counts():
    path = Path('/proc/self/io')
    return {k: int(v) for k, v in (line.split(':') for line in path.read_text().splitlines())} if path.exists() else {}


def worker(args):
    sys.path.insert(0, str(args.support_dir))
    import rescue_gene_models as rescue
    if args.candidate_support_dir:
        sys.path.insert(0, str(args.candidate_support_dir))
    filename = args.baseline_source if args.worker == 'baseline' else args.candidate_source
    text = filename.read_text()
    function = next(node for node in ast.parse(text).body
                    if isinstance(node, ast.FunctionDef) and node.name == 'copy_verified_comparison')
    namespace = dict(vars(rescue), stat=stat)
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(filename), 'exec'), namespace)
    operation = namespace['copy_verified_comparison']
    receipt = json.loads((args.input / 'receipt.json').read_text())
    with tempfile.TemporaryDirectory(prefix='comparison-copy-', dir=args.scratch) as temporary:
        before = io_counts()
        cpu_start = time.process_time()
        start = time.perf_counter()
        options = {'hash_workers': args.workers} if 'hash_workers' in inspect.signature(operation).parameters else {}
        operation(args.input, Path(temporary), receipt['key'], **options)
        seconds = time.perf_counter() - start
        cpu_seconds = time.process_time() - cpu_start
        after = io_counts()
        # Independent whole-file equality check is outside the measured copy.
        hashes = {name: rescue.digest(Path(temporary) / name) for name in receipt['files']}
        if hashes != receipt['files']:
            raise RuntimeError('Copied cache bytes differ from the producer receipt')
    scale = 1 if sys.platform != 'darwin' else 1 / 1024
    return dict(seconds=seconds, cpu_seconds=cpu_seconds, effective_workers=args.workers if options else 1,
                peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale,
                io={k: after[k] - before[k] for k in before}, output_hashes=hashes,
                function_sha256=hashlib.sha256(ast.get_source_segment(text, function).encode()).hexdigest(),
                receipt_sha256=rescue.digest(args.input / 'receipt.json'),
                proof_sha256=(hashlib.sha256(Path(inspect.getfile(sys.modules['refinement_receipt_snapshot'])).read_bytes()).hexdigest()
                              if 'refinement_receipt_snapshot' in sys.modules else None))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--baseline-source', required=True, type=Path)
    parser.add_argument('--candidate-source', type=Path, default=SUPPORT / 'rescue_gene_models.py')
    parser.add_argument('--support-dir', type=Path, default=SUPPORT)
    parser.add_argument('--candidate-support-dir', type=Path)
    parser.add_argument('--scratch', required=True, type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--worker', choices=('baseline', 'candidate'))
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args)))
        return
    if args.repeats < 1 or args.workers < 1 or args.output is None:
        parser.error('--output and positive --repeats are required')
    args.scratch.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    receipt_path = args.input / 'receipt.json'
    raw_receipt = receipt_path.read_bytes()
    receipt = json.loads(raw_receipt)
    results = {'baseline': [], 'candidate': []}
    for repetition in range(args.repeats + 1):
        modes = ('baseline', 'candidate') if repetition % 2 == 0 else ('candidate', 'baseline')
        for mode in modes:
            command = [sys.executable, __file__, '--worker', mode, '--workers', str(args.workers)]
            for option in ('input', 'baseline_source', 'candidate_source', 'support_dir', 'candidate_support_dir', 'scratch'):
                value = getattr(args, option)
                if value is not None:
                    command += ['--' + option.replace('_', '-'), str(value)]
            row = json.loads(subprocess.check_output(command, text=True,
                                                    env={**os.environ, 'PYTHONHASHSEED': '0'}))
            if row['output_hashes'] != receipt['files'] or row['receipt_sha256'] != hashlib.sha256(raw_receipt).hexdigest():
                raise RuntimeError('Cache copy or producer receipt changed between trials')
            if repetition:
                results[mode].append(row)
    if receipt_path.read_bytes() != raw_receipt:
        raise RuntimeError('Producer receipt changed during benchmark')
    summary = {mode: dict(median_seconds=statistics.median(r['seconds'] for r in rows),
                          median_cpu_seconds=statistics.median(r['cpu_seconds'] for r in rows),
                          peak_rss_kib=max(r['peak_rss_kib'] for r in rows),
                          median_rchar=statistics.median(r['io'].get('rchar', 0) for r in rows))
               for mode, rows in results.items()}
    report = dict(input=str(args.input), files=len(receipt['files']),
                  source_bytes=sum((args.input / p).stat().st_size for p in receipt['files']),
                  receipt_sha256=hashlib.sha256(raw_receipt).hexdigest(),
                  python=sys.version, platform=platform.platform(), repeats=args.repeats, workers=args.workers,
                  warmups=1, equivalent=True, results=results, summary=summary)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
