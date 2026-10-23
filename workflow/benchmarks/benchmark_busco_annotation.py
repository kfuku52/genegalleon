#!/usr/bin/env python3
"""Compare BUSCO annotation loading and the complete merge CLI.

Synthetic BUSCO full tables cover singleton, mixed and duplicate-rich genes.
Fresh alternating workers include a discarded warmup pair. Loader timers exclude
imports; CLI timers include imports and main. Fixture generation, output hashing
and schema verification are outside timers. All input bytes are still parsed.
"""

import argparse
import contextlib
import csv
import hashlib
import importlib.util
import io
import json
import platform
import re
import runpy
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def memory():
    return next(int(line.split()[1]) * 1024 for line in Path('/proc/self/status').read_text().splitlines()
                if line.startswith('VmHWM:'))


def load_module(path):
    spec = importlib.util.spec_from_file_location('busco_annotation_implementation', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def normalise_log(text):
    return re.sub(r'^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d(?:\.\d+)?: Processing:',
                  '<datetime>: Processing:', text, flags=re.MULTILINE)


def frame_signature(frame):
    return {'columns': frame.columns.tolist(), 'dtypes': list(map(str, frame.dtypes)),
            'rows': len(frame), 'index_dtype': str(frame.index.dtype), 'index_name': frame.index.name}


def worker(args):
    capture, errors = io.StringIO(), io.StringIO()
    if args.case == 'loader':
        module = load_module(args.implementation)
        start_wall, start_cpu = time.perf_counter(), time.process_time()
        with contextlib.redirect_stdout(capture), contextlib.redirect_stderr(errors):
            frame = module.load_busco(str(args.inputs / 'busco.tsv'))
        wall, cpu = time.perf_counter() - start_wall, time.process_time() - start_cpu
        peak = memory()
        signature = frame_signature(frame)
        frame.to_csv(args.worker_output, sep='\t', index=True)
        frame.to_pickle(args.worker_output.with_suffix('.pickle'))
    else:
        sys.argv = [str(args.implementation), '--cds_fasta', str(args.inputs / 'cds.fa'),
                    '--busco_tsv', str(args.inputs / 'busco.tsv'),
                    '--uniprot_tsv', str(args.inputs / 'uniprot.tsv'),
                    '--expression_tsv', str(args.inputs / 'expression.tsv'),
                    '--fx2tab', str(args.inputs / 'fx2tab.tsv'),
                    '--out_tsv', str(args.worker_output), '--ncpu', str(args.ncpu)]
        start_wall, start_cpu = time.perf_counter(), time.process_time()
        with contextlib.redirect_stdout(capture), contextlib.redirect_stderr(errors):
            runpy.run_path(str(args.implementation), run_name='__main__')
        wall, cpu = time.perf_counter() - start_wall, time.process_time() - start_cpu
        peak = memory()
        import pandas as pd
        frame = pd.read_csv(args.worker_output, sep='\t', dtype=str, keep_default_na=False)
        signature = frame_signature(frame)
    signature['sha256'] = sha256(args.worker_output)
    print(json.dumps({'wall_s': wall, 'cpu_s': cpu, 'peak_rss_bytes': peak,
                      'signature': signature, 'stdout': normalise_log(capture.getvalue()),
                      'stderr': errors.getvalue()}))


def fixtures(root, markers, profile, genes, cds_codons):
    root.mkdir()
    # The prime stride interleaves duplicate keys without changing their order.
    stride = 7919
    while markers % stride == 0:
        stride += 2
    order = sorted(range(markers), key=lambda value: (value * stride) % markers)
    with (root / 'busco.tsv').open('w', encoding='utf-8', newline='') as handle:
        handle.write('# Synthetic BUSCO full table\n')
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        for marker in order:
            if profile == 'singleton':
                gene = marker
            elif profile == 'mixed':
                gene = (marker // 10) * 9 + min(marker % 10, 8)
            else:
                gene = marker // 4
            sequence = f'gene_{gene:06d}:{10 + marker}-310'
            writer.writerow([f'{marker:06d}at9999', 'Complete', sequence,
                             '' if marker % 13 == 0 else 100 + marker % 5,
                             300 + marker % 7, '' if marker % 7 == 0 else f'https://example.invalid/{marker}',
                             '' if marker % 11 == 0 else f'説明 α {marker}\tfield\ncontinued'])
        writer.writerow(['missing', 'Missing', '', '', '', '', ''])
    with (root / 'cds.fa').open('w') as handle:
        for gene in range(genes):
            handle.write(f'>gene_{gene:06d}\n' + 'ATG' * cds_codons + '\n')
    for name, headers in [('uniprot.tsv', ['gene_id', 'protein_name']),
                          ('expression.tsv', ['Identifier', 'leaf', 'root']),
                          ('fx2tab.tsv', ['#id', 'length', 'GC'])]:
        with (root / name).open('w', encoding='utf-8', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(headers)
            for gene in range(genes):
                if name == 'uniprot.tsv':
                    values = [f'protein {gene}']
                elif name == 'expression.tsv':
                    values = [gene % 13, gene % 29]
                else:
                    values = [cds_codons * 3, 0.5]
                writer.writerow([f'gene_{gene:06d}', *values])


def main(args):
    sources = {'baseline': args.baseline.resolve()}
    if not args.baseline_only:
        sources['candidate'] = args.implementation.resolve()
    genes = args.genes or max(128, args.markers * 4)
    samples = []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='gg-busco-annotation-benchmark-') as temporary:
        root = Path(temporary)
        inputs = root / 'inputs'
        fixtures(inputs, args.markers, args.profile, genes, args.cds_codons)
        input_manifest = [{'name': path.name, 'bytes': path.stat().st_size, 'sha256': sha256(path)}
                          for path in sorted(inputs.iterdir())]
        references = {}
        reference_frames = {}
        for case in args.cases:
            for trial in range(-1, args.repeats):
                order = list(sources) if trial % 2 == 0 else list(reversed(sources))
                for variant in order:
                    output = root / f'{case}-{variant}-{trial}.tsv'
                    command = [sys.executable, str(Path(__file__).resolve()), '--worker', '--case', case,
                               '--implementation', str(sources[variant]), '--inputs', str(inputs),
                               '--worker-output', str(output), '--ncpu', str(args.ncpu)]
                    result = json.loads(subprocess.check_output(command, text=True))
                    if case == 'loader':
                        import pandas as pd
                        frame = pd.read_pickle(output.with_suffix('.pickle'))
                        reference_frame = reference_frames.setdefault(case, frame)
                        pd.testing.assert_frame_equal(frame, reference_frame, check_exact=True,
                                                      check_index_type=True, check_column_type=True)
                    log = result['stdout'] if args.ncpu == 1 else sorted(result['stdout'].splitlines())
                    evidence = (result['signature'], output.read_bytes(), log, result['stderr'])
                    reference = references.setdefault(case, evidence)
                    if evidence != reference:
                        raise AssertionError((case, variant, trial, result['signature']))
                    result.update(case=case, variant=variant, trial=trial, warmup=trial < 0)
                    samples.append(result)
                    print(json.dumps({key: result[key] for key in
                                      ('case', 'variant', 'trial', 'wall_s', 'cpu_s', 'peak_rss_bytes')}), flush=True)
    summary = {}
    for case in args.cases:
        summary[case] = {}
        for variant in sources:
            rows = [row for row in samples if row['case'] == case and row['variant'] == variant
                    and not row['warmup']]
            summary[case][variant] = {metric: statistics.median(row[metric] for row in rows)
                                      for metric in ('wall_s', 'cpu_s', 'peak_rss_bytes')}
        if len(sources) == 2:
            summary[case]['paired_speedups'] = [
                next(row['wall_s'] for row in samples if row['case'] == case
                     and row['variant'] == 'baseline' and row['trial'] == trial)
                / next(row['wall_s'] for row in samples if row['case'] == case
                       and row['variant'] == 'candidate' and row['trial'] == trial)
                for trial in range(args.repeats)]
    import pandas as pd
    report = {'environment': {'python': sys.version, 'pandas': pd.__version__, 'platform': platform.platform()},
              'configuration': {'markers': args.markers, 'profile': args.profile, 'genes': genes,
                                'cds_codons': args.cds_codons, 'ncpu': args.ncpu,
                                'repeats': args.repeats, 'cases': args.cases},
              'implementations': {key: {'path': str(path), 'sha256': sha256(path)} for key, path in sources.items()},
              'inputs': input_manifest, 'summary': summary, 'trials': samples}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'output': str(args.output), 'summary': summary}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--implementation', type=Path, default=Path('workflow/support/merge_cds_annotation.py'))
    parser.add_argument('--baseline-only', action='store_true')
    parser.add_argument('--markers', type=int, default=5000)
    parser.add_argument('--genes', type=int, default=0)
    parser.add_argument('--cds-codons', type=int, default=100)
    parser.add_argument('--profile', choices=['singleton', 'mixed', 'duplicated'], default='mixed')
    parser.add_argument('--ncpu', type=int, default=1)
    parser.add_argument('--cases', nargs='+', choices=['loader', 'cli'], default=['loader', 'cli'])
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--output', type=Path, default=Path('tmp/busco-annotation-benchmark.json'))
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--case', choices=['loader', 'cli'])
    parser.add_argument('--inputs', type=Path)
    parser.add_argument('--worker-output', type=Path)
    arguments = parser.parse_args()
    if arguments.worker:
        worker(arguments)
    else:
        if arguments.baseline is None or min(arguments.markers, arguments.ncpu,
                                             arguments.repeats, arguments.cds_codons) < 1:
            parser.error('--baseline and positive marker/cpu/repeat/codon counts are required')
        if arguments.genes < 0:
            parser.error('--genes cannot be negative')
        main(arguments)
