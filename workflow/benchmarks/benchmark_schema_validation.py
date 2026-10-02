#!/usr/bin/env python3
"""Measure strict scan-header validation on private wide TSVs."""
import argparse
import hashlib
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
    import generate_orthogroup_database as database
    from gene_family_output_store import GeneFamilyOutputStore, convert_storage_to_zip

    root = args.worker / 'families'
    directory = root / 'csubst_scan'
    directory.mkdir(parents=True)
    columns = sorted(database.CSUBST_SCAN_BASELINE_COLUMNS[database.AA_CHANGE_TABLE])
    columns += [f'optional_metric_{i:05d}' for i in range(args.columns - len(columns))]
    ids = [f'OG{i:07d}' for i in range(args.families)]
    header = '\t'.join(columns) + '\n'
    for family in ids:
        (directory / f'{family}_csubst_scan.tsv').write_text(header)
    invalid = directory / f'{ids[-1]}_csubst_scan.tsv'
    invalid_header = '\t'.join([*columns, 'duplicate_z', 'duplicate_a', 'duplicate_z', 'duplicate_a']) + '\n'
    invalid.write_text(invalid_header)
    if args.storage == 'zip':
        convert_storage_to_zip(root, 'orthogroup', ids, lambda name: name.split('_')[0])
    store = GeneFamilyOutputStore(root) if args.storage == 'zip' else None
    started = time.perf_counter()
    try:
        database.validate_csubst_scan_schemas([(database.AA_CHANGE_TABLE, str(directory))], store=store)
    except ValueError as exc:
        diagnostic = str(exc).replace(str(args.worker), '$FIXTURE')
    else:
        raise AssertionError('Duplicate columns must remain fatal')
    seconds = time.perf_counter() - started
    assert 'duplicate columns: duplicate_a, duplicate_z' in diagnostic
    print(json.dumps({'seconds': seconds, 'diagnostic_sha256': hashlib.sha256(diagnostic.encode()).hexdigest(),
                      'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--families', type=int, default=128)
    parser.add_argument('--columns', type=int, default=2048)
    parser.add_argument('--storage', choices=['files', 'zip'], default='files')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', type=Path)
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 1 <= args.families <= 1024 or not 16 <= args.columns <= 8192:
        parser.error('--families must be 1..1024 and --columns must be 16..8192')
    if args.worker:
        worker(args)
        return
    if args.output is None:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'families': args.families,
               'columns': args.columns, 'storage': args.storage,
               'python': sys.version, 'platform': platform.platform(), 'samples': []}
    with tempfile.TemporaryDirectory(prefix='gg-schema-benchmark-') as temporary:
        for trial in range(4):
            root = Path(temporary) / str(trial)
            root.mkdir()
            output = subprocess.check_output([sys.executable, str(Path(__file__).resolve()),
                '--worker', str(root), '--support-root', str(args.support_root),
                '--families', str(args.families), '--columns', str(args.columns),
                '--storage', args.storage], text=True)
            sample = json.loads(output)
            if trial:
                results['samples'].append(sample)
            print('warmup' if trial == 0 else trial, sample['seconds'], flush=True)
    assert len({row['diagnostic_sha256'] for row in results['samples']}) == 1
    results['median_seconds'] = statistics.median(row['seconds'] for row in results['samples'])
    results['median_peak_rss_kib_linux'] = statistics.median(row['peak_rss_kib_linux'] for row in results['samples'])
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    print('median', results['median_seconds'], results['median_peak_rss_kib_linux'], flush=True)


if __name__ == '__main__':
    main()
