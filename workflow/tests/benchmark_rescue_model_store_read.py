"""Synthetic full-record reader benchmark, including publication verification."""
import argparse
import hashlib
import importlib.util
import json
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def implementation(path):
    sys.path.insert(0, str(REPO / 'workflow/support'))
    spec = importlib.util.spec_from_file_location('benchmark_store', path)
    store = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(store)
    return store


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--implementation', required=True)
    p.add_argument('--baseline')
    p.add_argument('--records', type=int, default=20000)
    p.add_argument('--body-group-size', type=int, default=8)
    p.add_argument('--repeats', type=int, default=3)
    p.add_argument('--scope', choices=['models', 'partial', 'accepted', 'revision'], default='models')
    p.add_argument('--codec', choices=['gzip', 'zstd'], default='gzip')
    p.add_argument('--source')
    args = p.parse_args()
    if min(args.records, args.body_group_size, args.repeats) < 1:
        p.error('records, body-group-size and repeats must be positive')
    store = implementation(args.implementation)
    if args.source:
        started = time.perf_counter()
        rows = list(store._iter(args.source, args.scope))
        wall = time.perf_counter() - started
        maximum = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        checksum = hashlib.sha256()
        for row in rows:
            checksum.update(store._json(row) + b'\n')
        print(json.dumps({'wall_seconds': wall, 'peak_rss_bytes': maximum if sys.platform == 'darwin' else maximum * 1024,
                          'records': len(rows), 'full_records_sha256': checksum.hexdigest()}))
        return
    sys.path.insert(0, str(REPO))
    from workflow.tests.benchmark_rescue_model_store_patch import record
    from workflow.tests.test_rescue_model_store import model, receipt
    with tempfile.TemporaryDirectory(prefix='gg-model-read-') as name:
        source = Path(name) / 'source'
        def rows():
            for i in range(args.records):
                yield record(i, args.body_group_size, model)
        manifest = store.write_model_store(source, rows(), codec=args.codec, shard_bytes=65536)
        receipt(source)
        implementations = {'after': args.implementation}
        if args.baseline:
            implementations = {'before': args.baseline, **implementations}
        runs = {label: [] for label in implementations}
        for repeat in range(args.repeats):
            order = list(implementations)
            if repeat % 2:
                order.reverse()
            for label in order:
                command = [sys.executable, '-B', str(Path(__file__).resolve()), '--implementation', implementations[label],
                           '--source', str(source), '--scope', args.scope]
                result = subprocess.run(command, capture_output=True, text=True, check=True)
                runs[label].append(json.loads(result.stdout))
        if len({r['full_records_sha256'] for trials in runs.values() for r in trials}) != 1:
            raise ValueError('Full decoded records differ between benchmark runs')
        print(json.dumps({'records': args.records, 'body_group_size': args.body_group_size, 'scope': args.scope,
            'codec': manifest['codec'], 'sequence_nt': 12006, 'python': sys.version, 'measurements': runs,
            'implementation_sha256': hashlib.sha256(Path(args.implementation).read_bytes()).hexdigest(),
            'baseline_sha256': hashlib.sha256(Path(args.baseline).read_bytes()).hexdigest() if args.baseline else None,
            'timing_scope': 'Fully materialized verified reader; fixture creation and output hashing are excluded'}))


if __name__ == '__main__':
    main()
