#!/usr/bin/env python3
"""Compare complete GFF block interpretation and actual gff2genestat.main.

Sources are loaded as the real caller's structural-helper dependency. Imports,
fixtures, error corpus and full typed/output comparisons are outside timers.
Alternate fresh workers after one discarded warmup, retaining all measured pairs.

Run in the GeneGalleon runtime with a saved structural-helper baseline:
python workflow/benchmarks/benchmark_gff_block_attributes.py \\
    --baseline /mounted/baseline/gff_feature_structure.py
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

ROOT = Path(__file__).resolve()
while not (ROOT / 'workflow/support').is_dir() and ROOT != ROOT.parent:
    ROOT = ROOT.parent
SUPPORT = ROOT / 'workflow/support'
PROFILES = {'ordinary': (20000, 8), 'rich': (20000, 8), 'trans_splicing': (10000, 8),
            'overlap': (10000, 8), 'small': (64, 3), 'main': (2000, 8)}


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def model(kind, count):
    attrs = 'ID=c;Parent=t'
    if kind == 'rich':
        attrs += ';' + ';'.join(f'Note{i}=metadata%3Bvalue' for i in range(20))
    if kind == 'overlap':
        attrs += ';gg_source_overlap=confirmed'
    result = []
    for i in range(count):
        text = attrs + (f';exception=trans-splicing;part={i+1}/{count}' if kind == 'trans_splicing' else '')
        result.append(('chr1', '+', i * (8 if kind == 'overlap' else 15) + 1,
                       i * (8 if kind == 'overlap' else 15) + 9, text))
    return result


def corpus(module):
    values = []
    attrs = ['ID=c;Parent=t', 'ID=c;Parent=t;gg_source_overlap=confirmed',
             'ID=c;Parent=t; gg_source_overlap = confirmed ',
             'ID=c;Parent=t;gg_source_overlap=unconfirmed',
             'ID=c;Parent=t;gg_source_overlap=confirmed;gg_source_overlap=other',
             'ID=c;Parent=t;exception=ribosomal%20slippage',
             'ID=c;Parent=t;exception=trans-splicing;part=1/2',
             'ID=c;Parent=t;exception=other%2Ctrans-splicing;part=1/2',
             'ID=c;Parent=t;pseudo=true', 'ID=c;Parent=t;exception=trans-splicing,ribosomal slippage;part=1/2',
             'ID=c;Parent=t;exception=x;exception=y', 'ID=c;Parent=t;Note=a;Note=b',
             'ID=c;Parent=t;gg_source_overlap=confirmed;exception=trans-splicing;part=1/2',
             'ID=c;Parent=t;\u2003gg_source_overlap\u2003=CONFIRMED',
             'ID=c;Parent=t;exception=%ff', 'ID=c;Parent=t;exception=trans-splicing;part=2/2']
    for text in attrs:
        values.append(module.has_source_overlap_marker(text))
        for second in (text, text.replace('part=1/2', 'part=2/2')):
            try:
                value = module.ordered_annotated_blocks([('chr1', '+', 1, 9, text), ('chr1', '+', 8, 16, second)], 'g')
                values.append(('ok', value))
            except Exception as exc:
                values.append((type(exc).__name__, str(exc)))
    return values


def fixture(root, count, exons):
    directory = root / 'gff'
    directory.mkdir()
    fasta = root / 'queries.fa'
    with (directory / 'Genus_species.gff').open('w') as gff, fasta.open('w') as seq:
        gff.write('##gff-version 3\n')
        for i in range(count):
            start = i * 200 + 1
            gene = f'g{i:06}'
            seq.write(f'>Genus_species_{gene}\n' + 'ATG' * (exons * 3) + '\n')
            gff.write(f'chr1\tfixture\tgene\t{start}\t{start+exons*15-7}\t.\t+\t.\tID={gene}\n')
            gff.write(f'chr1\tfixture\tmRNA\t{start}\t{start+exons*15-7}\t.\t+\t.\tID={gene}.t;Parent={gene}\n')
            for exon in range(exons):
                pos = start + exon * 15
                gff.write(f'chr1\tfixture\tCDS\t{pos}\t{pos+8}\t.\t+\t0\tID={gene}.c;Parent={gene}.t\n')
    return fasta, directory


def normal(text):
    text = re.sub(r'\d{4}-\d\d-\d\d \d\d:\d\d:\d\d(?:\.\d+)?', '<time>', text)
    return re.sub(r'completed in [\d,]+ secs:', 'completed in <elapsed> secs:', text)


def worker(args):
    sys.path.insert(0, str(SUPPORT))
    module = load('gff_feature_structure', args.source)
    edge = corpus(module)
    count, exons = PROFILES[args.profile]
    caller = load('measured_gff_stats', SUPPORT / 'gff2genestat.py') if args.profile == 'main' else None
    if caller:
        assert caller.ordered_annotated_blocks is module.ordered_annotated_blocks
    rows = model(args.profile, exons)
    stdout, stderr = io.StringIO(), io.StringIO()
    output = Path(args.output_dir) / 'stats.tsv'
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
        wall_begin, cpu_begin = time.perf_counter(), time.process_time()
        if caller:
            sys.argv = ['gff2genestat.py', '--mode', 'gene_delim', '--dir_gff', args.gff_dir,
                        '--seqfile', args.fasta, '--outfile', str(output)]
            result = caller.main()
        else:
            result = [module.ordered_annotated_blocks(rows, 'g') for _ in range(count)]
        wall, cpu = time.perf_counter()-wall_begin, time.process_time()-cpu_begin
    rss = next(int(line.split()[1]) for line in Path('/proc/self/status').read_text().splitlines() if line.startswith('VmHWM:'))
    value = {'result': result, 'edge': edge, 'stdout': normal(stdout.getvalue()), 'stderr': normal(stderr.getvalue())}
    if caller:
        import pandas as pd
        value['bytes'] = output.read_bytes()
        value['frame'] = pd.read_csv(output, sep='\t')
        assert len(value['frame']) == count
    with Path(args.payload).open('wb') as handle:
        pickle.dump(value, handle, protocol=5)
    print(json.dumps({'wall_seconds': wall, 'cpu_seconds': cpu, 'rss_kib': rss}))


def main(args):
    import pandas as pd
    report = {'python': platform.python_version(), 'platform': platform.platform(), 'pandas': pd.__version__,
              'repeats': args.repeats, 'profiles': {},
              'sources': {label: {'path': str(getattr(args, label)),
                                  'sha256': hashlib.sha256(Path(getattr(args, label)).read_bytes()).hexdigest()}
                          for label in ('baseline', 'candidate')}}
    with tempfile.TemporaryDirectory(prefix='gg-block-attrs-') as temporary:
        root = Path(temporary)
        for profile in args.profiles:
            work = root / profile
            work.mkdir()
            fasta, directory = fixture(work, *PROFILES[profile]) if profile == 'main' else ('', '')
            inputs = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in work.rglob('*') if p.is_file()}
            runs = {'baseline': [], 'candidate': []}
            for trial in range(args.repeats+1):
                payloads = {}
                for label in ('baseline', 'candidate') if trial % 2 == 0 else ('candidate', 'baseline'):
                    out = root / f'{profile}-{trial}-{label}'
                    out.mkdir()
                    payload = out / 'result.pkl'
                    command = [sys.executable, str(Path(__file__).resolve()), '--worker', '--source', str(getattr(args, label)),
                               '--profile', profile, '--output-dir', str(out), '--payload', str(payload),
                               '--fasta', str(fasta), '--gff-dir', str(directory)]
                    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', PYTHONHASHSEED='0')
                    completed = subprocess.run(command, check=True, capture_output=True, text=True, env=env)
                    measured = json.loads(completed.stdout)
                    with payload.open('rb') as handle:
                        payloads[label] = pickle.load(handle)
                    if trial:
                        runs[label].append(measured)
                before, after = payloads['baseline'], payloads['candidate']
                for key in ('result', 'edge', 'stdout', 'stderr', 'bytes'):
                    assert before.get(key) == after.get(key), (profile, key)
                if profile == 'main':
                    pd.testing.assert_frame_equal(before['frame'], after['frame'], check_exact=True)
            assert inputs == {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in work.rglob('*') if p.is_file()}
            report['profiles'][profile] = runs
            print(profile, {label: statistics.median(row['wall_seconds'] for row in values) for label, values in runs.items()}, flush=True)
    Path(args.output).write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--candidate', type=Path, default=SUPPORT / 'gff_feature_structure.py')
    parser.add_argument('--output', type=Path, default=ROOT / 'tmp/gff-block-attributes-benchmark.json')
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--profiles', nargs='+', choices=PROFILES, default=list(PROFILES))
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--source')
    parser.add_argument('--profile')
    parser.add_argument('--fasta')
    parser.add_argument('--gff-dir')
    parser.add_argument('--output-dir')
    parser.add_argument('--payload')
    args = parser.parse_args()
    if args.worker:
        worker(args)
    else:
        if args.baseline is None:
            parser.error('--baseline is required')
        if args.repeats < 1:
            parser.error('--repeats must be positive')
        args.baseline = args.baseline.resolve(strict=True)
        args.candidate = args.candidate.resolve(strict=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        main(args)
