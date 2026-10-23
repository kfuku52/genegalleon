#!/usr/bin/env python3
"""Compare UniProt reference-name loading with an unchanged baseline module.

Run through workflow/tests/run_in_runtime.sh with --baseline /path/to/old.py and
--output tmp/uniprot-performance.json. Fixtures use mixed official/literal IDs,
180-aa sequences, sparse requested accessions, duplicates, raw/gzip input, and
dense/unfiltered controls. --profiles raw-cli gz-cli measures the complete
formatter on the same sparse family, excluding module import and input creation.

Each case uses one warmup and alternating fresh processes. Record wall/CPU time,
process peak RSS before output hashing, and exact canonical dictionary/TSV hashes.
Fixture bytes and implementation hashes make repeated measurements inspectable.
This isolates reference parsing; it does not establish whole-workflow speed.
"""
import argparse
import gzip
import hashlib
import importlib.metadata
import importlib.util
import json
import os
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def load(path):
    spec = importlib.util.spec_from_file_location('uniprot_under_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def accession(index):
    return f'P{index:07d}'


def token(index):
    return f'sp|{accession(index)}|ENTRY_{index}' if index % 5 else accession(index)


def fixture(root, rows):
    root.mkdir(parents=True, exist_ok=True)
    raw = root / f'uniprot-{rows}.fa'
    compressed = raw.with_suffix('.fa.gz')
    if not raw.exists():
        with raw.open('w') as output:
            for index in range(rows):
                output.write(f'>{token(index)} Example protein {index} with a realistic product name '
                             'OS=Example organism OX=123 GN=gene PE=1 SV=1\n')
                output.write('ACDEFGHIKLMNPQRSTVWY' * 9 + '\n')
            for index in range(0, rows, 1000):
                output.write(f'>{token(index)} Last product for {index} OS=Last organism OX=321\nM\n')
    if not compressed.exists():
        with raw.open('rb') as source, gzip.open(compressed, 'wb', compresslevel=3) as output:
            while block := source.read(1024 * 1024):
                output.write(block)
    hits = root / f'hits-{rows}.tsv'
    queries = root / f'queries-{rows}.fa'
    if not hits.exists() or not queries.exists():
        subjects = [accession(0)] + [accession(index) for index in range(1, rows, 1000)]
        subjects.extend(token(index) for index in range(2, rows, 1000))
        hits.write_text(''.join(f'gene{index}\t{subject}\t90\t180\t1e-30\t200\t180\n'
                                for index, subject in enumerate(subjects)))
        queries.write_text(''.join(f'>gene{index}\nM\n' for index in range(len(subjects))))
    return raw, compressed


def worker(args):
    module = load(args.script)
    mode, kind = args.profile.split('-')
    path = args.fixtures / f'uniprot-{args.rows}.fa'
    if mode == 'gz':
        path = path.with_suffix('.fa.gz')
    if kind == 'cli':
        with tempfile.TemporaryDirectory(prefix='gg-uniprot-benchmark-') as temporary:
            output = Path(temporary) / 'annotations.tsv'
            sys.argv = [str(args.script), '--diamond_tsv', str(args.fixtures / f'hits-{args.rows}.tsv'),
                        '--query_fasta', str(args.fixtures / f'queries-{args.rows}.fa'),
                        '--uniprot_fasta', str(path), '--outfile', str(output)]
            started = time.perf_counter()
            started_cpu = time.process_time()
            module.main()
            seconds = time.perf_counter() - started
            cpu_seconds = time.process_time() - started_cpu
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024)
            payload = output.read_bytes()
            proof = hashlib.sha256(payload).hexdigest()
            keys = payload.count(b'\n') - 1
    else:
        selected = selected_accessions(kind, args.rows)
        started = time.perf_counter()
        started_cpu = time.process_time()
        result = module.read_uniprot_descriptions(str(path), selected)
        seconds = time.perf_counter() - started
        cpu_seconds = time.process_time() - started_cpu
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024)
        proof = hashlib.sha256(json.dumps(result, sort_keys=True, ensure_ascii=False).encode()).hexdigest()
        keys = len(result)
    print(json.dumps(dict(seconds=seconds, cpu_seconds=cpu_seconds, peak_rss_bytes=rss,
                          sha256=proof, keys=keys)))


def selected_accessions(kind, rows):
    if kind == 'unfiltered':
        selected = None
    else:
        stride = 1000 if kind == 'sparse' else 1
        selected = {accession(index) for index in range(1, rows, stride)}
        # Include literal official tokens in addition to their saved accession.
        selected.update(token(index) for index in range(2, rows, stride))
        selected.add(accession(0))
    return selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--current', type=Path, default=Path(__file__).resolve().parents[1]
                        / 'support' / 'reformat_uniprot_diamond.py')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--fixtures', type=Path, default=Path('tmp/uniprot-description-benchmark'))
    parser.add_argument('--rows', type=int, default=300000)
    parser.add_argument('--trials', type=int, default=3)
    parser.add_argument('--profiles', nargs='+', choices=[f'{mode}-{kind}' for mode in ('raw', 'gz')
                        for kind in ('sparse', 'dense', 'unfiltered', 'cli')],
                        default=['raw-sparse', 'gz-sparse', 'raw-dense',
                                 'gz-dense', 'raw-unfiltered', 'gz-unfiltered'])
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--script', type=Path)
    parser.add_argument('--profile')
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    if args.baseline is None or args.output is None:
        parser.error('--baseline and --output are required')
    if args.rows < 1 or not 1 <= args.trials <= 20:
        parser.error('--rows must be positive and --trials must be 1..20')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    inputs = fixture(args.fixtures, args.rows)
    implementations = [('baseline', args.baseline)]
    if args.current:
        implementations.append(('current', args.current))
    results = dict(environment=dict(python=sys.version, platform=platform.platform(),
                                   pandas=importlib.metadata.version('pandas')),
                   rows=args.rows, trials=args.trials, warmup=1,
                   cache='warm filesystem caches; fresh process each sample', python_hash_seed='0',
                   scripts={label: dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                            for label, path in implementations},
                   inputs={str(path): dict(bytes=path.stat().st_size,
                                           sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for path in
                           (*inputs, args.fixtures / f'hits-{args.rows}.tsv',
                            args.fixtures / f'queries-{args.rows}.fa')},
                   samples={})
    for profile in args.profiles:
        groups = {label: [] for label, _ in implementations}
        for trial in range(args.trials + 1):
            for label, path in implementations if trial % 2 == 0 else implementations[::-1]:
                command = [sys.executable, str(Path(__file__).resolve()), '--worker', '--script', str(path.resolve()),
                           '--fixtures', str(args.fixtures.resolve()), '--rows', str(args.rows), '--profile', profile]
                sample = json.loads(subprocess.check_output(command, text=True,
                                    env=dict(os.environ, PYTHONHASHSEED='0')))
                print(profile, label, 'warmup' if trial == 0 else trial, sample, flush=True)
                if trial:
                    groups[label].append(sample)
        assert len({sample['sha256'] for group in groups.values() for sample in group}) == 1, profile
        results['samples'][profile] = groups
        results['medians'] = {profile: {label: {metric: statistics.median(sample[metric] for sample in samples)
                                               for metric in ('seconds', 'cpu_seconds', 'peak_rss_bytes')}
                                       for label, samples in groups.items()}
                              for profile, groups in results['samples'].items()}
        args.output.write_text(json.dumps(results, indent=2) + '\n')
    print(json.dumps(results['medians'], indent=2))


if __name__ == '__main__':
    main()
