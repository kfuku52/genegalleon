#!/usr/bin/env python3
"""Differential benchmark of complete BUSCO annotation loading and merge main.

Alternate fresh baseline/candidate workers after a discarded warmup pair. Compare
complete typed frames, serialized TSVs, stdout and exceptions outside timers.
Fixture construction and imports are excluded; Linux process peak RSS includes
imports. Use the GeneGalleon runtime, and do not run concurrent benchmarks.
"""
import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
import os
import pickle
import platform
import re
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def load(path):
    spec = importlib.util.spec_from_file_location('busco_candidate', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def worker(args):
    module = load(args.source)
    output, stderr = io.StringIO(), io.StringIO()
    result = error = None
    with contextlib.redirect_stdout(output), contextlib.redirect_stderr(stderr):
        start_cpu = time.process_time()
        start_wall = time.perf_counter()
        try:
            if args.cli:
                sys.argv = [args.source, '--cds_fasta', args.fasta, '--busco_tsv', args.fixture,
                            '--out_tsv', args.cli_output]
                module.main()
            else:
                result = module.load_busco(args.fixture)
        except Exception as exc:
            error = (type(exc).__name__, str(exc))
        wall = time.perf_counter() - start_wall
        cpu = time.process_time() - start_cpu
    rss = next(int(line.split()[1]) for line in Path('/proc/self/status').read_text().splitlines()
               if line.startswith('VmHWM:'))
    payload = {'result': result, 'error': error, 'stderr': stderr.getvalue(),
               'stdout': re.sub(r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}(?:\.\d+)?', '<time>', output.getvalue())}
    if args.cli and Path(args.cli_output).exists():
        payload['cli_bytes'] = Path(args.cli_output).read_bytes()
    with open(args.payload, 'wb') as handle:
        pickle.dump(payload, handle, protocol=5)
    print(json.dumps({'wall_seconds': wall, 'cpu_seconds': cpu, 'rss_kib': rss,
                      'source_sha256': hashlib.sha256(Path(args.source).read_bytes()).hexdigest()}))


def main(args, root):
    import pandas as pd
    profiles = []
    specs = [('small', 128, 2), ('singleton', 50000, 1), ('mixed', 50000, 0),
             ('duplicated', 10000, 2), ('dense', 10000, 8)]
    for name, genes, copies in specs:
        fixture = root / (name + '.tsv')
        if args.profiles and name not in args.profiles and not (name == 'duplicated' and 'cli' in args.profiles):
            continue
        with fixture.open('w') as handle:
            handle.write('# BUSCO synthetic annotation fixture\n')
            for copy in range(20 if copies == 0 else copies):
                # Interleave repeated genes to test source-order preservation.
                for gene in range(genes):
                    if copies == 0 and (copy > 0 and (gene % 20 != 0 or copy > 1)):
                        continue
                    handle.write(f'B{gene}_{copy}\tComplete\tg{gene}:cds{copy}\t{100 + copy}\t200\tu{gene}\td{gene}_{copy}\n')
        profiles.append((name, fixture, False))
    fasta = root / 'cds.fa'
    if not args.profiles or 'cli' in args.profiles:
        with fasta.open('w') as handle:
            for gene in range(10000):
                handle.write(f'>g{gene}\nATG\n')
        profiles.append(('cli', root / 'duplicated.tsv', True))
    edges = {
        'empty': b'', 'comments': b'# annotation\n',
        'missing': b'B1\tMissing\t\t\t\t\t\n',
        'all_filtered': b'B1\tComplete\t:one\t1\t2\tu\td\n',
        'lexical': b'B1\tComplete\t001:one\t100\t200\tu\td\nB2\tNA\tNA:two\t\t\tNA\tnan\nB3\tDuplicated\t001:three\t1.0\t2.0\tu\td\nB4\tComplete\tnan:four\t0\t0\t\t\n',
        'quoted': 'B1\tComplete\tgeneA:one\t100\t200\tu\t"説明 α\tfield\ncontinued"\nB2\tComplete\tgeneB:two\t\t\t\t\nB3\tDuplicated\tgeneA:three\t101\t201\tNA\tnan\n'.encode(),
        'duplicates_all': b'B1\tComplete\tg2:a\t1\t2\tu\td\nB2\tComplete\tg1:b\t2\t3\tu\td\nB3\tComplete\tg2:c\t3\t4\tu\td\nB4\tComplete\tg1:d\t4\t5\tu\td\n',
        'invalid_utf8': b'B1\tComplete\tg1\t1\t2\tu\td\nB2\tComplete\tg2\t1\t2\tu\t\xff\n',
        'extra_field': b'B1\tComplete\tg1\t1\t2\tu\td\nB2\tComplete\tg2\t1\t2\tu\td\textra\n',
        'bad_quote': b'B1\tComplete\tg1\t1\t2\tu\td\nB2\tComplete\tg2\t1\t2\tu\t"bad\n',
    }
    for name, data in edges.items():
        fixture = root / (name + '.tsv')
        fixture.write_bytes(data)
        profiles.append(('edge_' + name, fixture, False))
    if args.profiles:
        available = {name for name, *_ in profiles}
        assert set(args.profiles) <= available, set(args.profiles) - available
        profiles = [profile for profile in profiles if profile[0] in args.profiles]
    report = {'platform': platform.platform(), 'python': platform.python_version(), 'pandas': pd.__version__,
              'harness_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'sources': {label: {'path': getattr(args, label),
                                 'sha256': hashlib.sha256(Path(getattr(args, label)).read_bytes()).hexdigest()}
                          for label in ('baseline', 'candidate')}, 'repeats': args.repeats, 'profiles': {}}
    for name, fixture, cli in profiles:
        edge = name.startswith('edge_')
        runs = {'baseline': [], 'candidate': []}
        for trial in range(1 if edge else args.repeats + 1):
            payloads = {}
            for label in (['baseline', 'candidate'] if trial % 2 == 0 else ['candidate', 'baseline']):
                payload_path = root / (label + '.pkl')
                command = [sys.executable, str(Path(__file__).resolve()), '--worker', '--baseline', args.baseline,
                           '--candidate', args.candidate, '--source', getattr(args, label), '--fixture', str(fixture), '--payload', str(payload_path)]
                if cli:
                    cli_output = root / (label + '.tsv')
                    cli_output.unlink(missing_ok=True)
                    command += ['--cli', '--fasta', str(fasta), '--cli-output', str(cli_output)]
                env = os.environ.copy()
                env.update({'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'PYTHONHASHSEED': '0'})
                completed = subprocess.run(command, check=True, capture_output=True, text=True, env=env)
                measurement = json.loads(completed.stdout)
                with payload_path.open('rb') as handle:
                    payloads[label] = pickle.load(handle)
                if edge or trial:
                    runs[label].append(measurement)
            before, after = payloads['baseline'], payloads['candidate']
            assert before['error'] == after['error'], (name, before['error'], after['error'])
            assert before['stdout'] == after['stdout'], name
            assert before['stderr'] == after['stderr'], name
            if not edge:
                assert before['error'] is None, (name, before['error'])
                if cli:
                    assert before.get('cli_bytes'), name
                    assert len(before['cli_bytes'].splitlines()) == 10001, name
            if before['result'] is None:
                assert after['result'] is None, name
            else:
                pd.testing.assert_frame_equal(before['result'], after['result'], check_exact=True)
                assert before['result'].to_csv(sep='\t') == after['result'].to_csv(sep='\t'), name
            assert before.get('cli_bytes') == after.get('cli_bytes'), name
        report['profiles'][name] = runs
        if not edge:
            print(name, {label: statistics.median(row['wall_seconds'] for row in values)
                         for label, values in runs.items()}, flush=True)
    Path(args.output).write_text(json.dumps(report, indent=2) + '\n')
    print('All typed frames, TSV/CLI bytes, errors and logs match.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', required=True)
    parser.add_argument('--candidate', default=str(Path(__file__).resolve().parents[1] / 'support/merge_cds_annotation.py'))
    parser.add_argument('--output', default='tmp/busco-annotation-benchmark.json')
    parser.add_argument('--profiles', nargs='+')
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--source')
    parser.add_argument('--fixture')
    parser.add_argument('--payload')
    parser.add_argument('--cli', action='store_true')
    parser.add_argument('--fasta')
    parser.add_argument('--cli-output')
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    if args.worker:
        worker(args)
    else:
        args.baseline = str(Path(args.baseline).resolve(strict=True))
        args.candidate = str(Path(args.candidate).resolve(strict=True))
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='gg-busco-annotation-') as directory:
            main(args, Path(directory))
