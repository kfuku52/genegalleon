#!/usr/bin/env python3
"""Compare full dominant-partner ordering, preserving typed metadata/errors.

Inputs/imports/verification excluded; actual anchor reads and ordering included.
Fresh alternating workers discard one warmup pair. Optional target/query modes
are measured; the normal both_length default is an unchanged small control.

Run in the GeneGalleon runtime, with a saved support module as the baseline:
python workflow/benchmarks/benchmark_karyotype_order.py \\
    --baseline /mounted/baseline/pairwise_synteny.py
"""
import argparse
import hashlib
import importlib.util
import json
import os
import platform
import random
import statistics
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
PROFILES = {
    'small': (4, 16, 'target', 20),
    'chromosomes': (24, 2000, 'target', 1),
    'fragmented': (512, 128, 'target', 1),
    'fragmented_query': (512, 128, 'query', 1),
    'subset': (512, 128, 'subset', 1),
    'none': (24, 2000, 'none', 1),
    'both_length': (4, 16, 'both_length', 5),
}


def fingerprint(value):
    def typed(item):
        if isinstance(item, dict):
            return ['dict', [[typed(k), typed(v)] for k, v in item.items()]]
        if isinstance(item, (tuple, list)):
            return [type(item).__name__, [typed(v) for v in item]]
        return [type(item).__name__, item]
    return hashlib.sha256(json.dumps(typed(value), ensure_ascii=False).encode()).hexdigest()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(source):
    sys.path.insert(0, str(ROOT / 'workflow/support'))
    spec = importlib.util.spec_from_file_location('candidate_synteny_order', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fixture(root, profile):
    chromosomes, count, mode, batch = PROFILES[profile]
    selected = [[f'{prefix}{c}' for c in range(chromosomes)] for prefix in ('t', 'q')]
    genomes = [[], []]
    for side in range(2):
        for c, seqid in enumerate(selected[side]):
            for r in range(count):
                # Shared starts exercise end priority; exact coordinate ties
                # and natural gene-ID ordering are covered by contracts().
                genomes[side].append(SimpleNamespace(seqid=seqid, gene_id=f'{seqid}_g{r}',
                    start=(r // 2) * 100, end=(r // 2) * 100 + 60 + r % 3))
        random.Random(7231 + side).shuffle(genomes[side])
        random.Random(8291 + side).shuffle(selected[side])
    anchors = root / 'anchors'
    with anchors.open('w') as stream:
        stream.write('###\n\n')
        for c in range(chromosomes):
            for r in range(0, count, max(1, count // 16)):
                stream.write(f't{c}_g{r} q{(c * 7 + r) % chromosomes}_g{(r * 3) % count} 100\n')
                if r % 3 == 0:
                    stream.write(f't{c}_g{r} q{(c * 7 + r) % chromosomes}_g{(r * 3) % count} 99\n')
    simple = root / 'simple'
    simple.write_text('t0_g0 t0_g15 q1_g0 q1_g15 10 +\nt1_g0 t1_g15 q0_g0 q0_g15 10 +\n')
    if mode == 'subset':
        selected[1] = selected[1][:8]
        mode = 'target'
    return selected, genomes, anchors, simple, mode, batch


def contracts(module, root):
    selected, genomes, anchors, simple, _, _ = fixture(root, 'small')
    results = []
    for mode in ('target', 'query', 'none', 'both_length', 'target_length', 'query_length', 'bad'):
        try:
            result = module.order_karyotype(selected, genomes, anchors, mode, simple)
        except Exception as exc:
            result = [type(exc).__name__, str(exc)]
        results.append([mode, result])
    for line in ('missing q0_g0 100\n', 't0_g0 missing 100\n', 't0_g0\n', '', '###\n'):
        anchors.write_text(line)
        for mode in ('target', 'query'):
            try:
                result = module.order_karyotype(selected, genomes, anchors, mode)
            except Exception as exc:
                result = [type(exc).__name__, str(exc)]
            results.append([mode, result])
    # Direct helper calls can supply repeated IDs and repeated selected seqids.
    genomes[1].append(genomes[1][0])
    selected[1].append(selected[1][0])
    anchors.write_text('t0_g0 q0_g0 100\n' + f't1_g0 {genomes[1][0].gene_id} 100\n')
    results.append(module.order_karyotype(selected, genomes, anchors, 'target'))
    genomes[1].extend([SimpleNamespace(seqid='q0', gene_id=name, start=0, end=60)
                       for name in ('q0_tie10', 'q0_tie2')])
    anchors.write_text('t0_g0 q0_tie2 100\nt0_g1 q0_tie10 100\n')
    results.append(module.order_karyotype(selected, genomes, anchors, 'target'))
    return fingerprint(results)


def worker(args):
    module = load(args.source)
    selected, genomes, anchors, simple, mode, batch = fixture(args.work, args.profile)
    input_files = {'anchors': sha256(anchors), 'simple': sha256(simple)}
    before = fingerprint([selected, [[vars(g) for g in genes] for genes in genomes], input_files])
    start, cpu = time.perf_counter(), time.process_time()
    for _ in range(batch):
        result = module.order_karyotype(selected, genomes, anchors, mode, simple)
    wall = (time.perf_counter() - start) / batch
    cpu = (time.process_time() - cpu) / batch
    rss = next(int(line.split()[1]) * 1024 for line in Path('/proc/self/status').read_text().splitlines()
               if line.startswith('VmHWM:'))
    assert before == fingerprint([selected, [[vars(g) for g in genes] for genes in genomes],
                                  {'anchors': sha256(anchors), 'simple': sha256(simple)}])
    return dict(wall_seconds=wall, cpu_seconds=cpu, peak_rss_bytes=rss,
                input_sha256=before, output_sha256=fingerprint(result),
                input_files_sha256=input_files,
                contract_sha256=contracts(module, args.work),
                source=str(args.source), source_sha256=sha256(args.source))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, default=ROOT / 'workflow/support/pairwise_synteny.py')
    parser.add_argument('--profiles', nargs='+', choices=PROFILES, default=list(PROFILES))
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--warmups', type=int, default=1)
    parser.add_argument('--output', type=Path, default=Path('tmp/karyotype-order-benchmark.json'))
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--source', type=Path)
    parser.add_argument('--work', type=Path)
    parser.add_argument('--profile', choices=PROFILES)
    args = parser.parse_args()
    if args.repeats < 1 or args.warmups < 0:
        parser.error('--repeats must be positive and --warmups must be nonnegative')
    if args.worker:
        print(json.dumps(worker(args)))
        return
    work = args.output.with_suffix('')
    work.mkdir(parents=True, exist_ok=True)
    samples, summaries = [], {}
    for profile in args.profiles:
        records = []
        identity = None
        for pair in range(args.repeats + args.warmups):
            for variant in (('baseline', 'candidate') if pair % 2 == 0 else ('candidate', 'baseline')):
                source = getattr(args, variant).resolve()
                command = [sys.executable, str(Path(__file__).resolve()), '--worker', '--source', str(source),
                           '--baseline', str(args.baseline.resolve()),
                           '--profile', profile, '--work', str(work)]
                result = subprocess.run(command, check=True, capture_output=True, text=True,
                    env=dict(os.environ, PYTHONHASHSEED='0', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
                             MKL_NUM_THREADS='1'))
                row = json.loads(result.stdout)
                value = [row[key] for key in ('input_sha256', 'output_sha256', 'contract_sha256')]
                if identity is None:
                    identity = value
                assert value == identity, (profile, variant, value, identity)
                if pair >= args.warmups:
                    records.append(dict(row, variant=variant, pair=pair - args.warmups, profile=profile))
        summaries[profile] = {variant: {key: statistics.median(row[key] for row in records if row['variant'] == variant)
            for key in ('wall_seconds', 'cpu_seconds', 'peak_rss_bytes')} for variant in ('baseline', 'candidate')}
        summaries[profile]['paired_speedups'] = [
            next(row['wall_seconds'] for row in records if row['variant'] == 'baseline' and row['pair'] == pair)
            / next(row['wall_seconds'] for row in records if row['variant'] == 'candidate' and row['pair'] == pair)
            for pair in range(args.repeats)
        ]
        samples.extend(records)
        print(profile, json.dumps(summaries[profile]), flush=True)
    args.output.write_text(json.dumps(dict(environment=dict(python=sys.version, platform=platform.platform()),
        repeats=args.repeats, warmups=args.warmups,
        sources={name: {'path': str(getattr(args, name).resolve()), 'sha256': sha256(getattr(args, name))}
                 for name in ('baseline', 'candidate')},
        summaries=summaries, measurements=samples), indent=2) + '\n')


if __name__ == '__main__':
    main()
