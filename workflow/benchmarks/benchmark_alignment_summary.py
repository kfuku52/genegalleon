#!/usr/bin/env python3
"""Measure named alignment-statistics summaries on private raw/ZIP fixtures."""
import argparse
import contextlib
import hashlib
import importlib
import io
import json
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def worker(args):
    sys.path.insert(0, str(args.support_root))
    import pandas
    from gene_family_output_store import GeneFamilyOutputStore, convert_storage_to_zip

    module = importlib.import_module(args.reader + '_output_summary')
    families = [(f'HOG{i:07d}' if args.reader == 'orthogroup' else f'query{i:07d}_id.part')
                for i in range(args.families)]
    root = args.worker / 'output'
    subdir = 'alignment_stats_original'
    directory = root / subdir
    directory.mkdir(parents=True)
    columns = module._alignment_stats_columns()
    header = '\t'.join(columns) + '\n'
    for offset, family in enumerate(families):
        row = [7 + offset, 101, 707, 3, 0.2, 19, 11, 0.43]
        (directory / f'{family}_alignment_stats.original.tsv').write_text(
            header + '\t'.join(map(str, row)) + '\n')
    matchers = module._query_id_matchers(families) if args.reader == 'query2family' else None
    if args.storage == 'zip':
        identify = (lambda name: module._extract_query_id(name, matchers)) if matchers is not None else module._extract_orthogroup_id
        with contextlib.redirect_stdout(io.StringIO()):
            convert_storage_to_zip(root, args.reader, families, identify)
    kwargs = {'query_id_matchers': matchers} if matchers is not None else {}
    if args.storage == 'zip':
        kwargs.update(store=GeneFamilyOutputStore(root), logical_subdir=subdir)
    frame = pandas.DataFrame({'Total': [2] * len(families)}, index=families)
    started = time.perf_counter()
    with contextlib.redirect_stdout(io.StringIO()):
        result = module.get_alignment_stats(frame, str(directory), 'original', ncpu=args.ncpu, **kwargs)
    seconds = time.perf_counter() - started
    assert result.shape == (args.families, 9)
    for offset, family in enumerate(families):
        assert result.loc[family, 'No_of_taxa_original'] == 7 + offset
    fingerprint = hashlib.sha256(result.to_csv(sep='\t').encode()).hexdigest()
    print(json.dumps({'seconds': seconds, 'output_sha256': fingerprint,
                      'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--families', type=int, default=10000)
    parser.add_argument('--ncpu', type=int, default=4)
    parser.add_argument('--reader', choices=['orthogroup', 'query2family'], default='orthogroup')
    parser.add_argument('--storage', choices=['files', 'zip'], default='files')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', type=Path)
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 1 <= args.families <= 10000 or not 1 <= args.ncpu <= 32:
        parser.error('--families must be 1..10000 and --ncpu must be 1..32')
    if args.worker:
        worker(args)
        return
    if args.output is None:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'families': args.families,
               'reader': args.reader, 'storage': args.storage, 'ncpu': args.ncpu,
               'python': sys.version, 'platform': platform.platform(), 'samples': []}
    with tempfile.TemporaryDirectory(prefix='gg-stats-benchmark-') as temporary:
        for trial in range(4):
            root = Path(temporary) / str(trial)
            root.mkdir()
            output = subprocess.check_output([sys.executable, str(Path(__file__).resolve()),
                '--worker', str(root), '--support-root', str(args.support_root),
                '--families', str(args.families), '--reader', args.reader,
                '--storage', args.storage, '--ncpu', str(args.ncpu)], text=True)
            sample = json.loads(output)
            if trial:
                results['samples'].append(sample)
            print('warmup' if trial == 0 else trial, sample['seconds'], flush=True)
    assert len({row['output_sha256'] for row in results['samples']}) == 1
    results['median_seconds'] = statistics.median(row['seconds'] for row in results['samples'])
    results['median_peak_rss_kib_linux'] = statistics.median(row['peak_rss_kib_linux'] for row in results['samples'])
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    print('median', results['median_seconds'], results['median_peak_rss_kib_linux'], flush=True)


if __name__ == '__main__':
    main()
