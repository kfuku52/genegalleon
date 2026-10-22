#!/usr/bin/env python3
"""Compare command-local rescue dependency verification on synthetic files.

This measures key construction and its full-byte dependency checks, not rescue,
refinement, model admission, or whole-workflow throughput. No project data is read.
"""
import argparse
import hashlib
import importlib.util
import itertools
import json
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SUPPORT = REPO / 'workflow/support'


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def load_source(path):
    sys.path.insert(0, str(SUPPORT))
    spec = importlib.util.spec_from_file_location('measured_rescue_keys', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True) + '\n')


def fixture(root, count, payload_bytes, entries):
    names = [f'Fixture_{i:03d}' for i in range(count)]
    pairs = list(itertools.combinations(names, 2))
    pairs = pairs[:min(158, len(pairs))]
    jobs = [dict(a=a, b=b, id=f'pair_{i:03d}', index=i, kind='pair') for i, (a, b) in enumerate(pairs, 1)]
    plan = dict(species=names, synteny_jobs=jobs,
                donors={name: sorted({j['b'] if j['a'] == name else j['a'] for j in jobs if name in (j['a'], j['b'])}) for name in names},
                request={'files': {f'synthetic_input_{i}': hashlib.sha256(str(i).encode()).hexdigest() for i in range(entries)}})
    write_json(root / 'plan.json', plan)
    plan_hash = hashlib.sha256((root / 'plan.json').read_bytes()).hexdigest()
    prepared = {}
    chunk = b'synthetic dependency verification payload\n' * 25000
    for name in names:
        path = root / 'prepared' / name / 'payload.blob'
        path.parent.mkdir(parents=True)
        value = hashlib.sha256()
        with path.open('wb') as stream:
            remaining = payload_bytes
            while remaining:
                part = chunk[:min(len(chunk), remaining)]
                stream.write(part)
                value.update(part)
                remaining -= len(part)
        receipt = path.parent / 'receipt.json'
        write_json(receipt, {'key': {'plan': plan_hash, 'species': name}, 'files': {'payload.blob': value.hexdigest()}})
        prepared[name] = hashlib.sha256(receipt.read_bytes()).hexdigest()
    for job in jobs:
        directory = root / 'synteny' / job['id']
        directory.mkdir(parents=True)
        (directory / 'blocks.json').write_bytes(b'[]\n')
        key = {'plan': plan_hash, 'job': job, 'prepared': {n: prepared[n] for n in sorted({job['a'], job['b']})}}
        write_json(directory / 'receipt.json', {'key': key, 'files': {'blocks.json': hashlib.sha256(b'[]\n').hexdigest()}})
    return dict(species=count, comparisons=len(jobs), prepared_payload_bytes=payload_bytes,
                plan_entries=entries, plan_bytes=(root / 'plan.json').stat().st_size)


def measure(source, root, passes):
    module = load_source(source)
    state = sys.modules[module.digest.__module__]
    counts = Counter()
    state.count = lambda name, amount=1: counts.update({name: amount})
    plan = json.loads((root / 'plan.json').read_text())
    started = time.perf_counter()
    keys = []
    with module.ReceiptSnapshot(module.verified) as snapshot:
        for _ in range(passes):
            for name in plan['species']:
                keys.append(module.rescue_key(root, plan, name, dependency_snapshot=snapshot))
    elapsed = time.perf_counter() - started
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return dict(wall_seconds=elapsed, max_rss_bytes=rss if sys.platform == 'darwin' else rss * 1024,
                output_sha256=digest(keys), **counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compare', help='Exact previous commit; no branches or defaults are changed')
    parser.add_argument('--species', type=int, default=23)
    parser.add_argument('--payload-mib', type=int, default=8)
    parser.add_argument('--plan-entries', type=int, default=20000)
    parser.add_argument('--passes', type=int, default=2)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--warmups', type=int, default=1)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--fixture-root', type=Path)
    args = parser.parse_args()
    if args.source:
        print(json.dumps(measure(args.source, args.fixture_root, args.passes), sort_keys=True))
        return
    if min(args.species, args.payload_mib, args.passes, args.repeats) < 1 or args.warmups < 0:
        parser.error('Positive workload and repeat sizes are required')
    current = SUPPORT / 'rescue_gene_models.py'
    with tempfile.TemporaryDirectory(prefix='gg-dependency-keys-') as temporary:
        directory = Path(temporary)
        root = directory / 'fixture'
        workload = fixture(root, args.species, args.payload_mib * 1024**2, args.plan_entries)
        sources = {'after': current}
        if args.compare:
            before = directory / 'before.py'
            before.write_bytes(subprocess.check_output(
                ['git', '-c', 'safe.directory=' + str(REPO), 'show',
                 args.compare + ':workflow/support/rescue_gene_models.py'], cwd=REPO))
            sources = {'before': before, **sources}
        results = {label: [] for label in sources}
        expected = None
        for trial in range(args.warmups + args.repeats):
            order = list(sources)
            if trial % 2:
                order.reverse()
            for label in order:
                command = [sys.executable, '-B', str(Path(__file__).resolve()), '--source', str(sources[label]),
                           '--fixture-root', str(root), '--passes', str(args.passes)]
                value = json.loads(subprocess.check_output(command, text=True, cwd=REPO))
                expected = expected or value['output_sha256']
                if value['output_sha256'] != expected:
                    raise ValueError('Complete ordered dependency keys differ')
                if trial >= args.warmups:
                    results[label].append(value)
        summary = {label: {key: statistics.median(row[key] for row in values)
                           for key in ('wall_seconds', 'max_rss_bytes', 'sha256_bytes', 'sha256_reads')}
                   for label, values in results.items()}
        print(json.dumps(dict(workload=workload, passes=args.passes, warmups=args.warmups,
                              repeats=args.repeats, baseline_commit=args.compare,
                              output_sha256=expected, medians=summary, trials=results), sort_keys=True))


if __name__ == '__main__':
    main()
