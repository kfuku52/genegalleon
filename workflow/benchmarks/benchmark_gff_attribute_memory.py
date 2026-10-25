#!/usr/bin/env python3
"""Compare GFF attribute identification with fresh, alternating workers.

Synthetic fixtures and SHA-256 calculations run outside workers. Timers include
implementation import and identify_attribute_key; peak RSS covers the complete
worker process. These warm-cache measurements do not measure JCVI or a complete
annotation workflow. Run without concurrent benchmarks, tests, or builds.
"""

import argparse
import gzip
import hashlib
import importlib.util
import json
import platform
import resource
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from types import SimpleNamespace

SUPPORT = Path(__file__).resolve().parents[1] / 'support'


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def worker(args):
    sys.path.insert(0, str(SUPPORT))
    started_wall, started_cpu = time.perf_counter(), time.process_time()
    spec = importlib.util.spec_from_file_location('gff_attribute_benchmark_target', args.implementation)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    try:
        value = module.identify_attribute_key(SimpleNamespace(
            fasta_file=str(args.fasta), gff_file=str(args.gff), gff_feature='mRNA'))
        outcome = {'returned': value}
    except Exception as error:
        outcome = {'exception': type(error).__module__ + '.' + type(error).__qualname__,
                   'message': str(error)}
    wall, cpu = time.perf_counter() - started_wall, time.process_time() - started_cpu
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(json.dumps({'wall_s': wall, 'cpu_s': cpu,
                      'peak_rss_bytes': peak if sys.platform == 'darwin' else peak * 1024,
                      'pandas': module.pandas.__version__, 'outcome': outcome}))


def fixtures(root, args):
    fasta = root / 'genes.fa'
    with fasta.open('w', encoding='utf-8', newline='') as handle:
        for index in range(args.fasta_records):
            handle.write(f'>Arabidopsis_thaliana_gene{index:06d}\nATGATGATG\n')
    cases = []
    sizes = [('large', args.rows), ('small', args.small_rows)]
    for size, rows in sizes:
        positions = args.positions + (['malformed'] if args.include_errors else [])
        for position in positions:
            plain = root / f'{size}-{position}.gff3'
            match_row = {'first': 0, 'last': rows - 1, 'missing': -1, 'malformed': 0}[position]
            with plain.open('w', encoding='utf-8', newline='') as handle:
                handle.write('##gff-version 3\n')
                for index in range(rows):
                    feature = 'mRNA' if index % 5 == 0 or index == match_row else 'CDS'
                    attributes = (f'ID=unmatched{index:09d};Parent=parent{index // 5};'
                                  f'Note=synthetic%20annotation%20{index}')
                    if index == match_row:
                        attributes = 'Name=gene000000;ID=gene000000;Note=first%20matching%20key'
                    handle.write(f'chr1\tsynthetic\t{feature}\t{index * 3 + 1}\t{index * 3 + 3}'
                                 f'\t.\t+\t.\t{attributes}\n')
                if position == 'malformed':
                    handle.write('chr1\tsynthetic\tmRNA\t1\t3\t.\t+\t.\tID=extra\textra\n')
            for codec in args.codecs:
                path = plain
                if codec == 'gzip':
                    path = plain.with_suffix(plain.suffix + '.gz')
                    with plain.open('rb') as source, path.open('wb') as destination:
                        with gzip.GzipFile(filename='', mode='wb', fileobj=destination, mtime=0) as compressed:
                            shutil.copyfileobj(source, compressed)
                cases.append({'case': f'{size}-{position}-{codec}', 'gff': path,
                              'rows': rows, 'codec': codec, 'position': position,
                              'expected': ({'exception': 'pandas.errors.ParserError'}
                                           if position == 'malformed' else
                                           {'returned': None if position == 'missing' else 'Name'})})
                if args.include_errors and position == 'first' and codec == 'gzip':
                    truncated = path.with_name(f'{size}-truncated.gff3.gz')
                    shutil.copyfile(path, truncated)
                    with truncated.open('r+b') as handle:
                        handle.truncate(truncated.stat().st_size - 5)
                    cases.append({'case': f'{size}-truncated-gzip', 'gff': truncated,
                                  'rows': rows, 'codec': codec, 'position': 'truncated',
                                  'expected': {'exception': 'builtins.EOFError'}})
    return fasta, cases


def main(args):
    sources = {'baseline': args.baseline.resolve()}
    if not args.baseline_only:
        sources['candidate'] = args.implementation.resolve()
    implementations = {key: {'path': str(path), 'sha256': sha256(path)}
                       for key, path in sources.items()}
    samples, summary = [], {}
    with tempfile.TemporaryDirectory(prefix='gg-gff-attribute-memory-') as temporary:
        fasta, cases = fixtures(Path(temporary), args)
        inputs = [{'name': fasta.name, 'bytes': fasta.stat().st_size, 'sha256': sha256(fasta)}]
        for case in cases:
            path, name = case['gff'], case['case']
            reference = None
            inputs.append({'case': name, 'name': path.name, 'rows': case['rows'],
                           'bytes': path.stat().st_size, 'sha256': sha256(path),
                           'codec': case['codec'], 'position': case['position']})
            for trial in range(-1, args.repeats):
                order = list(sources)
                if trial % 2:
                    order.reverse()
                for variant in order:
                    command = [sys.executable, str(Path(__file__).resolve()), '--worker',
                               '--implementation', str(sources[variant]), '--fasta', str(fasta),
                               '--gff', str(path)]
                    result = json.loads(subprocess.check_output(command, text=True))
                    expected = case['expected']
                    if any(key not in result['outcome'] or result['outcome'][key] != value
                           for key, value in expected.items()):
                        raise AssertionError((name, variant, trial, case['expected'], result['outcome']))
                    if reference is None:
                        reference = result['outcome']
                    elif result['outcome'] != reference:
                        raise AssertionError((name, variant, trial, reference, result['outcome']))
                    result.update(case=name, variant=variant, trial=trial, warmup=trial < 0)
                    samples.append(result)
                    print(json.dumps(result), flush=True)
            summary[name] = {}
            for variant in sources:
                measured = [row for row in samples if row['case'] == name
                            and row['variant'] == variant and not row['warmup']]
                summary[name][variant] = {metric: statistics.median(row[metric] for row in measured)
                                         for metric in ('wall_s', 'cpu_s', 'peak_rss_bytes')}
            if len(sources) == 2:
                before, after = summary[name]['baseline'], summary[name]['candidate']
                summary[name]['candidate_over_baseline'] = {
                    metric: after[metric] / before[metric]
                    for metric in ('wall_s', 'cpu_s', 'peak_rss_bytes')}
    report = {'environment': {'python': sys.version, 'platform': platform.platform(),
                              'pandas': samples[0]['pandas']},
              'measurement': {'timers': 'implementation import and identify_attribute_key',
                              'rss': 'whole fresh worker process', 'cache': 'warm filesystem cache',
                              'warmups_per_variant': 1},
              'configuration': {'rows': args.rows, 'small_rows': args.small_rows,
                                'fasta_records': args.fasta_records, 'positions': args.positions,
                                'codecs': args.codecs, 'repeats': args.repeats,
                                'include_errors': args.include_errors},
              'implementations': implementations, 'inputs': inputs,
              'summary': summary, 'trials': samples}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'output': str(args.output), 'summary': summary}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--implementation', type=Path, default=SUPPORT / 'gff_attribute_identifier.py')
    parser.add_argument('--baseline-only', action='store_true')
    parser.add_argument('--rows', type=int, default=500000)
    parser.add_argument('--small-rows', type=int, default=32)
    parser.add_argument('--fasta-records', type=int, default=10000)
    parser.add_argument('--positions', nargs='+', choices=['first', 'last', 'missing'],
                        default=['first', 'last', 'missing'])
    parser.add_argument('--codecs', nargs='+', choices=['plain', 'gzip'], default=['plain', 'gzip'])
    parser.add_argument('--include-errors', action='store_true',
                        help='Also verify malformed rows and truncated gzip after an early match.')
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--output', type=Path, default=Path('tmp/gff-attribute-memory.json'))
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--fasta', type=Path)
    parser.add_argument('--gff', type=Path)
    arguments = parser.parse_args()
    if arguments.worker:
        worker(arguments)
    else:
        if arguments.baseline is None or min(arguments.rows, arguments.small_rows,
                                             arguments.fasta_records, arguments.repeats) < 1:
            parser.error('--baseline and positive row/FASTA/repeat counts are required')
        if arguments.include_errors and 'first' not in arguments.positions:
            parser.error('--include-errors requires first in --positions')
        main(arguments)
