#!/usr/bin/env python3
"""Compare BUSCO metadata reduction and the complete collector CLI.

Use a saved baseline source with --baseline. Inputs are synthetic full BUSCO
tables, including duplicated markers, missing genes, complementary metadata,
literal NA-like text and quoted Unicode/tab/newline descriptions. Both the
metadata kernel (compiled unchanged from main's AST) and the complete CLI run
in alternating fresh workers, with a discarded warmup pair. The CLI timer
includes imports and main execution; fixture generation and output hashing are
outside timers. VmHWM includes imports and, for the kernel, fixture loading.
"""

import argparse
import ast
import contextlib
import csv
import hashlib
import importlib.util
import io
import json
import platform
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
    spec = importlib.util.spec_from_file_location('busco_benchmark_implementation', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def metadata_kernel(path, module):
    tree = ast.parse(path.read_text(), filename=str(path))
    main = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == 'main')
    block = next(node for node in main.body if isinstance(node, ast.If)
                 and ast.unparse(node.test) == 'len(metadata_frames) > 0')
    unit = ast.parse('def metadata_kernel(metadata_frames):\n    pass\n')
    unit.body[0].body = [block, ast.Return(value=ast.Name(id='merged_metadata', ctx=ast.Load()))]
    ast.fix_missing_locations(unit)
    namespace = {'pd': module.pd, 'first_non_missing': module.first_non_missing}
    exec(compile(unit, str(path), 'exec'), namespace)
    return namespace['metadata_kernel']


def signature(frame, output):
    payload = frame.to_csv(sep='\t', index=False).encode()
    output.write_bytes(payload)
    return {'columns': frame.columns.tolist(), 'dtypes': list(map(str, frame.dtypes)),
            'rows': len(frame), 'index_dtype': str(frame.index.dtype), 'index_name': frame.index.name,
            'sha256': hashlib.sha256(payload).hexdigest()}


def worker(args):
    sys.path.insert(0, str(args.support.resolve()))
    capture = io.StringIO()
    if args.case == 'kernel':
        module = load_module(args.implementation)
        frames = []
        for path in sorted(args.inputs.glob('*.tsv')):
            table = module.read_busco_table(path)
            frames.append(table.loc[:, ['busco_id', 'orthodb_url', 'description']].drop_duplicates(
                subset=['busco_id'], keep='first', ignore_index=True))
        function = metadata_kernel(args.implementation, module)
        start_wall, start_cpu = time.perf_counter(), time.process_time()
        result = function(frames)
        wall, cpu = time.perf_counter() - start_wall, time.process_time() - start_cpu
        peak = memory()
        result_signature = signature(result, args.worker_output)
    else:
        sys.argv = [str(args.implementation), '--busco_outdir', str(args.inputs),
                    '--outfile', str(args.worker_output), '--ncpu', str(args.ncpu)]
        start_wall, start_cpu = time.perf_counter(), time.process_time()
        with contextlib.redirect_stdout(capture):
            runpy.run_path(str(args.implementation), run_name='__main__')
        wall, cpu = time.perf_counter() - start_wall, time.process_time() - start_cpu
        peak = memory()
        import pandas as pd
        result = pd.read_csv(args.worker_output, sep='\t', dtype=str, keep_default_na=False)
        result_signature = {'columns': result.columns.tolist(), 'dtypes': list(map(str, result.dtypes)),
                            'rows': len(result), 'sha256': sha256(args.worker_output)}
    print(json.dumps({'wall_s': wall, 'cpu_s': cpu, 'peak_rss_bytes': peak,
                      'signature': result_signature, 'stdout': capture.getvalue()}))


def fixtures(root, markers, species):
    root.mkdir()
    for sp in range(species):
        path = root / f'Species_{sp:03d}.tsv'
        with path.open('w', encoding='utf-8', newline='') as handle:
            handle.write('# Synthetic BUSCO full table\n')
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            for offset in range(markers):
                marker = (offset + sp * 13) % markers
                identifier = f'{marker:06d}at9999'
                missing = marker % 13 == 0
                sequence = '-' if missing else f'gene_{sp:03d}_{marker:06d}:1-300'
                url = '-' if (marker + sp) % 3 else f'https://example.invalid/{marker}'
                description = '-' if (marker + sp) % 5 else f'説明 α {marker}\tfield\ncontinued'
                if marker % 71 == 0:
                    url, description = 'NA', 'nan'
                writer.writerow([identifier, 'Missing' if missing else 'Complete', sequence,
                                 '100', '300', url, description])
                if marker % 17 == 0:
                    writer.writerow([identifier, 'Duplicated', f'copy_{sp}_{marker}:1-300',
                                     '99', '300', 'later-url', 'later-description'])


def main(args):
    sources = {'baseline': args.baseline.resolve()}
    if not args.baseline_only:
        sources['candidate'] = args.implementation.resolve()
    samples = []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='gg-busco-metadata-benchmark-') as temporary:
        root = Path(temporary)
        inputs = root / 'inputs'
        fixtures(inputs, args.markers, args.species)
        input_manifest = [{'name': path.name, 'bytes': path.stat().st_size, 'sha256': sha256(path)}
                          for path in sorted(inputs.glob('*.tsv'))]
        references = {}
        for case in args.cases:
            for trial in range(-1, args.repeats):
                order = list(sources) if trial % 2 == 0 else list(reversed(sources))
                for variant in order:
                    output = root / f'{case}-{variant}-{trial}.tsv'
                    command = [sys.executable, str(Path(__file__).resolve()), '--worker', '--case', case,
                               '--implementation', str(sources[variant]), '--support', str(args.support.resolve()),
                               '--inputs', str(inputs), '--worker-output', str(output), '--ncpu', str(args.ncpu)]
                    result = json.loads(subprocess.check_output(command, text=True))
                    content = output.read_bytes()
                    log = result['stdout'] if args.ncpu == 1 else sorted(result['stdout'].splitlines())
                    evidence = (result['signature'], content, log)
                    reference = references.setdefault(case, evidence)
                    assert evidence == reference, (case, variant, trial)
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
              'configuration': {'markers': args.markers, 'species': args.species, 'ncpu': args.ncpu,
                                'repeats': args.repeats, 'cases': args.cases},
              'implementations': {key: {'path': str(path), 'sha256': sha256(path)} for key, path in sources.items()},
              'inputs': input_manifest, 'summary': summary, 'trials': samples}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'output': str(args.output), 'summary': summary}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--implementation', type=Path, default=Path('workflow/support/collect_common_BUSCO_genes.py'))
    parser.add_argument('--support', type=Path, default=Path('workflow/support'))
    parser.add_argument('--baseline-only', action='store_true')
    parser.add_argument('--markers', type=int, default=2000)
    parser.add_argument('--species', type=int, default=10)
    parser.add_argument('--ncpu', type=int, default=1)
    parser.add_argument('--cases', nargs='+', choices=['kernel', 'cli'], default=['kernel', 'cli'])
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--output', type=Path, default=Path('tmp/busco-metadata-benchmark.json'))
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--case', choices=['kernel', 'cli'])
    parser.add_argument('--inputs', type=Path)
    parser.add_argument('--worker-output', type=Path)
    arguments = parser.parse_args()
    if arguments.worker:
        worker(arguments)
    else:
        if arguments.baseline is None or min(arguments.markers, arguments.species, arguments.repeats) < 1:
            parser.error('--baseline and positive marker/species/repeat counts are required')
        main(arguments)
