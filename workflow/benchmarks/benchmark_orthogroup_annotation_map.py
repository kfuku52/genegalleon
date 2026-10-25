#!/usr/bin/env python3
"""Compare orthogroup annotation maps and the complete merge parser/main.

Both sources run in alternating fresh processes after discarded warmups. Imports,
fixture creation, output hashing and complete typed Series/TSV/CLI verification
are outside timers. Peak RSS includes imports. Thirteen parser boundary cases
also compare returned types, logs and errors. Run in the GeneGalleon runtime:

python workflow/benchmarks/benchmark_orthogroup_annotation_map.py \\
    --baseline /mounted/baseline/merge_cds_annotation.py
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

ROOT = Path(__file__).resolve().parents[2]


def load(path):
    spec = importlib.util.spec_from_file_location('mapping_candidate', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def normalize_log(text):
    return re.sub(r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}(?:\.\d+)?', '<timestamp>', text)


def worker(args):
    mod = load(args.source)
    out = io.StringIO()
    result = None
    error = None
    with contextlib.redirect_stdout(out):
        start_cpu = time.process_time()
        start = time.perf_counter()
        try:
            if args.cli:
                sys.argv = [args.source, '--scientific_name', 'Species', '--orthogroup_tsv', args.fixture,
                            '--cds_fasta', args.fasta, '--out_tsv', args.cli_output]
                mod.main()
            else:
                result = mod.load_orthogroup_map(args.fixture, args.species)
        except Exception as exc:
            error = (type(exc).__name__, str(exc))
        seconds = time.perf_counter() - start
        cpu_seconds = time.process_time() - start_cpu
    rss = next(int(line.split()[1]) for line in Path('/proc/self/status').read_text().splitlines()
               if line.startswith('VmHWM:'))
    payload = {'result': result, 'error': error, 'stdout': normalize_log(out.getvalue())}
    if args.cli and Path(args.cli_output).exists():
        payload['cli_bytes'] = Path(args.cli_output).read_bytes()
    with open(args.payload, 'wb') as handle:
        pickle.dump(payload, handle, protocol=5)
    print(json.dumps({'wall_seconds': seconds, 'cpu_seconds': cpu_seconds, 'rss_kib': rss,
                      'source_sha256': hashlib.sha256(Path(args.source).read_bytes()).hexdigest()}))


def main(args, root):
    import pandas as pd
    profiles = []
    for name, rows, copies, missing_every, species_columns in [
        ('small', 64, 1, 0, 1),
        ('singleton', 50000, 1, 0, 1),
        ('multicopy', 50000, 4, 0, 1),
        ('giant_members', 1000, 200, 0, 1),
        ('sparse', 50000, 1, 2, 1),
        ('wide', 20000, 2, 0, 20),
    ]:
        if args.profiles and name not in args.profiles and not (name == 'singleton' and 'cli' in args.profiles):
            continue
        fixture = root / (name+'.tsv')
        with fixture.open('w') as handle:
            handle.write('Orthogroup\tSpecies'+''.join('\tOther'+str(i) for i in range(species_columns-1))+'\n')
            for i in range(rows):
                # Duplicates cross group boundaries; first-group assignment must win.
                members = ', '.join('g'+str((i*copies+j) % max(1, rows*copies-17)) for j in range(copies))
                if missing_every and i % missing_every:
                    members = ''
                handle.write(('OG'+str(i) if i % 37 else '')+'\t'+members+
                             ''.join('\tOther_gene'+str(i) for _ in range(species_columns-1))+'\n')
        profiles.append((name, fixture, False, 'Species'))
    fasta = root/'cds.fa'
    if not args.profiles or 'cli' in args.profiles:
        with fasta.open('w') as handle:
            for i in range(50000):
                handle.write('>g'+str(i)+'\nATG\n')
        profiles.append(('cli', root/'singleton.tsv', True, 'Species'))
    # All differential cases use the real parser, including invalid tails.
    edge_cases = {
        'lexical': b'Orthogroup\tSpecies\n001\t001, NA, nan, None\nNA\t None , , other\n\tnan,000\n',
        'whitespace': 'Orthogroup\tSpecies\nA\t , \u2003, \t\nB\t\u2003a\u2003, a\n'.encode(),
        'empty': b'Orthogroup\tSpecies\n',
        'all_missing': b'Orthogroup\tSpecies\nA\t\n',
        'all_filtered': b'Orthogroup\tSpecies\nA\t , , \n',
        'missing_groups': b'Orthogroup\tSpecies\n\ta,b\nX\tb,c\n',
        'missing_header': b'Other\tSpecies\nA\ta\n',
        'duplicate_header': b'Orthogroup\tSpecies\tSpecies\nA\ta\tb\n',
        'quoted': b'Orthogroup\tSpecies\n"OG\tX"\t"a, b\ncontinued"\n',
        'invalid_utf8': b'Orthogroup\tSpecies\nA\ta\nB\t\xff\n',
        'bad_quote': b'Orthogroup\tSpecies\nA\ta\nB\t"bad\n',
        'blank': b'',
    }
    for name, data in edge_cases.items():
        path = root/(name+'.tsv')
        path.write_bytes(data)
        profiles.append(('edge_'+name, path, False, 'Species'))
    profiles.append(('edge_same_column', root/'lexical.tsv', False, 'Orthogroup'))
    if args.profiles:
        available = {name for name, *_ in profiles}
        unknown = set(args.profiles) - available
        if unknown:
            raise ValueError('Unknown profiles: ' + ', '.join(sorted(unknown)))
        profiles = [profile for profile in profiles if profile[0] in args.profiles]
    report = {'platform': platform.platform(), 'python': platform.python_version(), 'pandas': pd.__version__,
              'harness_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'sources': {label: {'path': getattr(args, label),
                                 'sha256': hashlib.sha256(Path(getattr(args, label)).read_bytes()).hexdigest()}
                          for label in ('baseline', 'candidate')},
              'repeats': args.repeats, 'profiles': {}}
    for name, fixture, cli, species in profiles:
        is_edge = name.startswith('edge_')
        runs = {'baseline': [], 'candidate': []}
        for trial in range(1 if is_edge else args.repeats+1):
            order = ['baseline', 'candidate'] if trial % 2 == 0 else ['candidate', 'baseline']
            payloads = {}
            for label in order:
                payload_path = root/('payload-'+label+'.pkl')
                cli_output = root/('cli-'+label+'.tsv')
                command = [sys.executable, str(Path(__file__).resolve()), '--worker', '--baseline', args.baseline, '--source',
                           args.baseline if label == 'baseline' else args.candidate,
                           '--fixture', str(fixture), '--payload', str(payload_path), '--species', species]
                if cli:
                    command += ['--cli', '--fasta', str(fasta), '--cli-output', str(cli_output)]
                env = os.environ.copy()
                env.update({'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
                            'PYTHONHASHSEED': '0'})
                result = subprocess.run(command, env=env, capture_output=True, text=True, check=True)
                measurement = json.loads(result.stdout)
                with payload_path.open('rb') as handle:
                    payloads[label] = pickle.load(handle)
                if is_edge or trial:
                    runs[label].append(measurement)
            before = payloads['baseline']
            after = payloads['candidate']
            if not is_edge:
                assert before['error'] is None, (name, before['error'])
                assert after['error'] is None, (name, after['error'])
                if cli:
                    assert before.get('cli_bytes'), name
                    assert len(before['cli_bytes'].splitlines()) == 50001, name
                else:
                    assert before['result'] is not None and len(before['result']) > 0, name
            assert before['error'] == after['error'], (name, before['error'], after['error'])
            assert before['stdout'] == after['stdout'], (name, before['stdout'], after['stdout'])
            if before['result'] is None:
                assert after['result'] is None, name
            else:
                pd.testing.assert_series_equal(before['result'], after['result'], check_exact=True)
                assert before['result'].to_csv(sep='\t') == after['result'].to_csv(sep='\t'), name
            assert before.get('cli_bytes') == after.get('cli_bytes'), name
        report['profiles'][name] = runs
        if not is_edge:
            print(name, {label: statistics.median(row['wall_seconds'] for row in values)
                         for label, values in runs.items()}, flush=True)
    Path(args.output).write_text(json.dumps(report, indent=2)+'\n')
    print('All complete typed outputs, errors, logs and CLI bytes match.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', required=True, help='Saved baseline merge_cds_annotation.py.')
    parser.add_argument('--candidate', default=str(ROOT / 'workflow/support/merge_cds_annotation.py'))
    parser.add_argument('--output', default='tmp/orthogroup-annotation-map-benchmark.json')
    parser.add_argument('--profiles', nargs='+', help='Optional subset of named normal or edge profiles.')
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--source')
    parser.add_argument('--fixture')
    parser.add_argument('--payload')
    parser.add_argument('--species', default='Species')
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
        with tempfile.TemporaryDirectory(prefix='gg-annotation-map-benchmark-') as temporary:
            main(args, Path(temporary))
