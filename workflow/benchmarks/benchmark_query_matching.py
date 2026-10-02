#!/usr/bin/env python3
"""Measure query ownership matching and complete summaries on private fixtures."""
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
from types import SimpleNamespace


def worker(args):
    sys.path.insert(0, str(args.support_root))
    import gene_family_output_store as store
    import query2family_output_summary as summary
    ids = ['Q', 'Q_a', 'Q_a.b'] + [f'query{i:05d}_id.part' for i in range(args.families)]
    if args.case == 'matching':
        names = [name + suffix for name in ids for suffix in ('_stat.branch.tsv', '.tree_plot.pdf')]
        names += [f'unknown{i:05d}_stat.branch.tsv' for i in range(args.families)]
        started = time.perf_counter()
        matchers = store.query_id_matchers(ids)
        if hasattr(store, 'query_id_extractor'):
            extract = store.query_id_extractor(matchers)
        else:
            def extract(name):
                return store.query_id_from_name(name, matchers)
        result = [extract(name) for name in names]
        seconds = time.perf_counter() - started
        fingerprint = hashlib.sha256(json.dumps(result).encode()).hexdigest()
    else:
        query_dir = args.worker / 'input/query_gene'
        output_root = args.worker / 'output/query2family'
        query_dir.mkdir(parents=True)
        for name in ids:
            (query_dir / name).write_text('gene\n')
        for subdir, suffix in [('stat_branch', '_stat.branch.tsv'), ('tree_plot', '_tree_plot.pdf'),
                               ('alignment', '_alignment.fa')]:
            directory = output_root / subdir
            directory.mkdir(parents=True)
            for name in ids:
                (directory / (name + suffix)).write_text('fixture\n')
        output = args.worker / 'summary.tsv'
        started = time.perf_counter()
        with contextlib.redirect_stdout(io.StringIO()):
            summary.run(SimpleNamespace(dir_query2family=str(output_root), dir_query_gene=str(query_dir),
                                        out=str(output), ncpu=1))
        seconds = time.perf_counter() - started
        fingerprint = hashlib.sha256(output.read_bytes()).hexdigest()
    print(json.dumps({'seconds': seconds, 'output_sha256': fingerprint,
                      'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--families', type=int, default=3000)
    parser.add_argument('--case', choices=['matching', 'summary', 'all'], default='all')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', type=Path)
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if args.families < 1 or args.families > 10000:
        parser.error('--families must be 1..10000')
    if args.worker:
        worker(args)
        return
    if args.output is None:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'families': args.families + 3,
               'python': sys.version, 'platform': platform.platform(), 'cases': {}}
    with tempfile.TemporaryDirectory(prefix='gg-query-benchmark-') as temporary:
        for case in (['matching', 'summary'] if args.case == 'all' else [args.case]):
            samples = []
            for trial in range(4):
                root = Path(temporary) / f'{case}-{trial}'
                root.mkdir()
                output = subprocess.check_output([sys.executable, str(Path(__file__).resolve()),
                    '--worker', str(root), '--case', case, '--support-root', str(args.support_root),
                    '--families', str(args.families)], text=True)
                sample = json.loads(output)
                if trial:
                    samples.append(sample)
            assert len({row['output_sha256'] for row in samples}) == 1
            results['cases'][case] = {'samples': samples,
                'median_seconds': statistics.median(row['seconds'] for row in samples)}
            args.output.write_text(json.dumps(results, indent=2) + '\n')
            print(case, results['cases'][case]['median_seconds'], flush=True)


if __name__ == '__main__':
    main()
