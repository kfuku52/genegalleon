"""Compare karyotype, scaffold plotting and FASTQ profiling in fresh workers.

Use --prepare --fixtures DIRECTORY once, then supply --baseline MODULE and
--source MODULE for one --case (count-large/small/sparse, render-large/small,
scatter-large/small, scaffold-render-large/small, resource-fastq-large/small,
resource-fastq-gzip, or resource-fasta-control). Each measured pair alternates
order after one discarded warmup; --output preserves all samples and hashes.

Count cases time only the actual BED-count statements extracted from the
renderer. Scatter cases time the real run prefix through scaffold aggregation
and the three jittered scatter helpers. The fourth scatter and complete figure
are covered by scaffold-render cases. Imports and fixture generation are outside
timers. Complete rendering includes reads and SVG/PDF writes; decoded profiling
includes its complete scan and mutation fence. RSS is Linux worker VmHWM.
"""
import argparse
import ast
import hashlib
import importlib.util
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / 'workflow/support'
sys.path.insert(0, str(SUPPORT))


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bed_count_function(source):
    # Extract the real count assignment and its optional local preparation.
    # BED parsing/drawing are intentionally outside this phase measurement.
    tree = ast.parse(source.read_text())
    target = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'render_karyotype')
    scope = next(n for n in ast.walk(target) if isinstance(n, ast.With) and
                 any(isinstance(k, ast.For) and ast.unparse(k.target) == '(entry, line)' for k in n.body))
    loop = next(n for n in scope.body if isinstance(n, ast.For) and ast.unparse(n.target) == '(entry, line)')
    nodes = [n for n in loop.body if isinstance(n, ast.Assign) and any(
        ast.unparse(t) in ('selected_seqids', 'chromosome_counts', 'entry.sizes') for t in n.targets)]
    code = 'def count(entry):\n' + '\n'.join('    ' + ast.unparse(n) for n in nodes) + '\n    return entry.sizes\n'
    namespace = {}
    from collections import Counter
    namespace['Counter'] = Counter
    exec(code, namespace)
    return namespace['count']


def sample(args):
    if args.case.startswith('scaffold-render'):
        from types import SimpleNamespace

        import matplotlib
        import matplotlib.pyplot as plt
        module = load(args.source, 'scaffold_candidate')
        options = SimpleNamespace(fx2tab=(args.fixtures / (args.case + '.tsv')).resolve(),
                                  min_scaffold_size=100)
        directory = args.fixtures / args.case
        directory.mkdir(exist_ok=True)
        os.chdir(directory)
        os.environ['SOURCE_DATE_EPOCH'] = '1780000000'
        with matplotlib.rc_context({'svg.hashsalt': 'gg-scaffold-performance-20261026'}):
            start, cpu = time.perf_counter(), time.process_time()
            module.run(options)
            wall, cpu = time.perf_counter() - start, time.process_time() - cpu
        output = {ext: hashlib.sha256(Path('scaffold_size_histogram.' + ext).read_bytes()).hexdigest()
                  for ext in ['svg', 'pdf']}
        plt.close('all')
    elif args.case.startswith('count'):
        from types import SimpleNamespace

        from jcvi.formats.bed import Bed
        function = bed_count_function(args.source)
        contigs = 32 if args.case == 'count-large' else 2
        bed_name = args.case + '.bed'
        entry = SimpleNamespace(bed=Bed(str(args.fixtures / bed_name)),
                                seqids=['c' + str(i) for i in range(contigs - 1, -1, -1)])
        start, cpu = time.perf_counter(), time.process_time()
        output = function(entry)
        wall, cpu = time.perf_counter() - start, time.process_time() - cpu
    elif args.case.startswith('resource'):
        module = load(args.source, 'resource_candidate')
        kind = 'fasta' if 'fasta' in args.case else 'fastq'
        path = args.fixtures / (args.case + ('.gz' if 'gzip' in args.case else '.txt'))
        start, cpu = time.perf_counter(), time.process_time()
        output = module.scan(path, kind)
        wall, cpu = time.perf_counter() - start, time.process_time() - cpu
    elif args.case.startswith('scatter'):
        import numpy as np
        source = args.source.read_text().replace('    matplotlib = get_matplotlib()', '    return df\n    matplotlib = get_matplotlib()')
        namespace = {'__name__': 'scatter_candidate'}
        exec(compile(source, str(args.source), 'exec'), namespace)
        size = 1000000 if args.case.endswith('large') else 100
        from types import SimpleNamespace
        options = SimpleNamespace(fx2tab=args.fixtures / (args.case + '.tsv'), min_scaffold_size=100)
        values = []
        class Axis:
            def scatter(self, x, y, **kwargs):
                values.append([x, y, kwargs])
        rng = np.random.default_rng(11)
        start, cpu = time.perf_counter(), time.process_time()
        df = namespace['run'](options)
        for y in ('length', 'GC', 'GC-Skew'):
            namespace['scatter_by_size'](Axis(), df, y, '≥1.0E+02 bp', '<1.0E+02 bp', rng)
        wall, cpu = time.perf_counter() - start, time.process_time() - cpu
        checksum = hashlib.sha256()
        for x, y, kwargs in values:
            checksum.update(x.tobytes())
            checksum.update(y.tobytes())
            checksum.update(json.dumps(kwargs, sort_keys=True).encode())
        assert len(df) == size
        output = checksum.hexdigest()
    elif args.case.startswith('render'):
        import matplotlib
        from kffractbias.io import read_bed
        module = load(args.source, 'karyotype_candidate')
        directory = args.fixtures / args.case
        seqids = ['c' + str(i) for i in range(32 if args.case.endswith('large') else 2)]
        genomes = [read_bed(directory / (side + '.bed')) for side in ('target', 'query')]
        colors = module.chromosome_colors([seqids, seqids], genomes, directory / 'unused')
        pair = {'target_species': 'Target_species', 'query_species': 'Query_species'}
        os.environ['SOURCE_DATE_EPOCH'] = '1780000000'
        with matplotlib.rc_context({'svg.hashsalt': 'gg-performance-20261026'}):
            start, cpu = time.perf_counter(), time.process_time()
            report = module.render_karyotype(directory, pair, 'svg', colors)
            wall, cpu = time.perf_counter() - start, time.process_time() - cpu
        output = {'report': report, 'svg': hashlib.sha256((directory / 'karyotype.svg').read_bytes()).hexdigest()}
    else:
        raise ValueError(args.case)
    rss = next(int(line.split()[1]) for line in Path('/proc/self/status').read_text().splitlines()
               if line.startswith('VmHWM:'))
    raw = json.dumps(output, sort_keys=True, allow_nan=True).encode()
    print(json.dumps({'wall': wall, 'cpu': cpu, 'rss_kib': rss,
                      'output_sha256': hashlib.sha256(raw).hexdigest()}))


def fixtures(path):
    path.mkdir(parents=True, exist_ok=True)
    (path / 'count-sparse.bed').write_text(''.join(
        f'c{i % 100000}\t{i * 10}\t{i * 10 + 5}\tg{i}\n' for i in range(200000)))
    for size, n, contigs in [('large', 200000, 32), ('small', 100, 2)]:
        bed = ''.join(f'c{i % contigs}\t{i * 10}\t{i * 10 + 5}\tg{i}\n' for i in range(n))
        (path / ('count-' + size + '.bed')).write_text(bed)
        directory = path / ('render-' + size)
        directory.mkdir(exist_ok=True)
        for side, prefix in [('target', 't'), ('query', 'q')]:
            (directory / (side + '.bed')).write_text(bed.replace('\tg', '\t' + prefix))
        seqids = ','.join('c' + str(i) for i in range(contigs))
        (directory / 'seqids').write_text(seqids + '\n' + seqids + '\n')
        (directory / 'layout').write_text(
            f'0.7,0.12,0.92,0,,Target species,top,{directory / "target.bed"},top\n'
            f'0.3,0.12,0.92,0,,Query species,bottom,{directory / "query.bed"},bottom\n'
            f'e,0,1,{directory / "colored.simple"}\n')
        (directory / 'colored.simple').write_text(f'#88afc4*t0 t{contigs} q0 q{contigs} 5 +\n')
    for size, n in [('large', 1000000), ('small', 100)]:
        with (path / ('scatter-' + size + '.tsv')).open('w') as stream:
            stream.write('#id\tlength\tGC\tGC-Skew\n')
            for i in range(n):
                stream.write(f'g{i}\t{50 if i % 8 else 200}\t{i % 100}\t{(i % 100 - 50)/100}\n')
    for case, n in [('scaffold-render-large', 20000), ('scaffold-render-small', 100)]:
        with (path / (case + '.tsv')).open('w') as stream:
            stream.write('#id\tlength\tGC\tGC-Skew\n')
            for i in range(n):
                stream.write(f'g{i}\t{50 if i % 8 else 200}\t{i % 100}\t{(i % 100 - 50)/100}\n')
    import gzip
    for case, n in [('resource-fastq-large', 500000), ('resource-fastq-small', 100),
                    ('resource-fastq-gzip', 500000), ('resource-fasta-control', 500000)]:
        is_gzip = 'gzip' in case
        target = path / (case + ('.gz' if is_gzip else '.txt'))
        opener = gzip.open if is_gzip else open
        with opener(target, 'wb') as stream:
            for i in range(n):
                if 'fasta' in case:
                    stream.write(b'>r%d\n%s\n' % (i, b'ACGT' * 37 + b'AC'))
                else:
                    stream.write(b'@r%d\r\n%s\r\n+\r\n%s\r\n' % (i, b'ACGT' * 37 + b'AC', b'I' * 150))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--case')
    parser.add_argument('--fixtures', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--prepare', action='store_true')
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    if args.prepare:
        fixtures(args.fixtures.resolve())
        return
    if args.worker:
        sample(args)
        return
    if not args.case or not args.output or not (args.baseline or args.source) or args.repeats < 1:
        parser.error('--case, --output, at least one source, and positive --repeats are required')
    sources = {'baseline': args.baseline.resolve()} if args.baseline else {}
    if args.source:
        sources['candidate'] = args.source.resolve()
    result = {'case': args.case, 'platform': platform.platform(), 'python': sys.version,
              'sources': {k: {'path': str(v), 'sha256': hashlib.sha256(v.read_bytes()).hexdigest()}
                          for k, v in sources.items()}, 'samples': []}
    for trial in range(args.repeats + 1):
        for name in list(sources)[::(-1 if trial % 2 else 1)]:
            command = [sys.executable, __file__, '--worker', '--source', str(sources[name]),
                       '--case', args.case, '--fixtures', str(args.fixtures)]
            output = subprocess.run(command, text=True, capture_output=True, check=True)
            value = json.loads(output.stdout.strip().splitlines()[-1])
            if trial:
                result['samples'].append({'implementation': name, **value})
            print(trial, name, value['wall'], value['rss_kib'], flush=True)
    result['output_equivalent'] = len({v['output_sha256'] for v in result['samples']}) == 1
    result['medians'] = {name: {metric: statistics.median(v[metric] for v in result['samples']
                                                      if v['implementation'] == name)
                               for metric in ['wall', 'cpu', 'rss_kib']} for name in sources}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    if not result['output_equivalent']:
        raise AssertionError('complete output hashes differ; inspect ' + str(args.output))
    print(result['medians'])


if __name__ == '__main__':
    main()
