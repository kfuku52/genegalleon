#!/usr/bin/env python3
"""Measure verified restoration and idle BUSCO storage in equivalent fresh workers.

The baseline and candidate consume the same sealed ZIP. Fixture generation,
copying and byte/metadata comparisons are outside the timed operation. One
warmup pair is discarded; measured pairs alternate their execution order.
"""

import argparse
import hashlib
import json
import os
import platform
import resource
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def fingerprint(root):
    digest = hashlib.sha256()
    count = 0
    for path in sorted(root.rglob('*')):
        if path.is_file():
            digest.update(path.relative_to(root).as_posix().encode())
            digest.update(path.read_bytes())
            value = path.stat()
            digest.update(str((value.st_mode & 0o7777, value.st_mtime_ns)).encode())
            count += 1
    return {'sha256': digest.hexdigest(), 'files': count}


def zip_fingerprint(path):
    digest = hashlib.sha256()
    with zipfile.ZipFile(path) as archive:
        data = json.loads(archive.read('__gg_directory_manifest__.json'))
        for row in sorted(data['files'], key=lambda row: Path(row['path'])):
            digest.update(row['path'].encode())
            digest.update(archive.read(Path(data['directory']).name + '/' + row['path']))
            digest.update(str((row['mode'], row['mtime_ns'])).encode())
    return {'sha256': digest.hexdigest(), 'files': len(data['files'])}


def worker(args):
    sys.path.insert(0, str(args.support))
    import directory_output_store as store
    name = 'busco_mafft' if args.case == 'idle' else 'tree_plot_input'
    with tempfile.TemporaryDirectory(prefix='gg-directory-worker-') as temporary:
        root = Path(temporary) / ('genome_evolution' if args.case == 'idle' else 'hgt')
        root.mkdir()
        shutil.copyfile(args.fixture, root / (name + '.zip'))
        if args.case == 'idle':
            source = (args.support.parent / 'core/gg_genome_evolution_core.sh').read_text()
            script = 'set -euo pipefail\n'
            for function in ('genome_evolution_storage_require_native_paths', 'genome_evolution_storage_prepare', 'genome_evolution_storage_finish'):
                if function + '() {' not in source:
                    continue
                script += function + '() {' + source.split(function + '() {', 1)[1].split('\n}\n', 1)[0] + '\n}\n'
            script += 'genome_evolution_storage_mode=zip\ngenome_evolution_storage_lock_token=""\n'
            script += 'run_busco_dupaware_mafft=0\ngenome_evolution_storage_prepare\ngenome_evolution_storage_finish\n'
            env = {key: value for key, value in os.environ.items() if not key.startswith('run_busco_dupaware_')}
            env.update(gg_support_dir=str(args.support), gg_workspace_output_dir=str(root.parent),
                       dir_genome_evolution=str(root))
        before_usage = [resource.getrusage(kind) for kind in (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN)]
        before_cpu = sum(row.ru_utime + row.ru_stime for row in before_usage)
        started = time.perf_counter()
        if args.case == 'idle':
            subprocess.run(['bash', '-c', script], env=env, check=True, capture_output=True, text=True)
        else:
            store.convert(root, kind='hgt', to='raw', directories=[name])
        elapsed = time.perf_counter() - started
        after_usage = [resource.getrusage(kind) for kind in (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN)]
        cpu = sum(row.ru_utime + row.ru_stime for row in after_usage) - before_cpu
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        child_rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        multiplier = 1024 if sys.platform != 'darwin' else 1
        return {'seconds': elapsed, 'cpu_seconds': cpu, 'peak_rss_bytes': int(rss * multiplier),
                'peak_child_process_rss_bytes': int(child_rss * multiplier),
                'fingerprint': zip_fingerprint(root / (name + '.zip')) if args.case == 'idle' else fingerprint(root / name)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=('restore', 'idle'), default='restore')
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--candidate', type=Path, default=ROOT / 'workflow/support')
    parser.add_argument('--files', type=int, default=2000)
    parser.add_argument('--bytes-per-file', type=int, default=4096)
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--support', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--fixture', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.support is not None:
        print(json.dumps(worker(args)))
        return
    if args.baseline is None or min(args.files, args.bytes_per_file, args.repeats) < 1:
        parser.error('a baseline support directory and positive workload sizes are required')
    sys.path.insert(0, str(args.baseline))
    import directory_output_store as baseline
    with tempfile.TemporaryDirectory(prefix='gg-directory-fixture-') as temporary:
        root = Path(temporary)
        name = 'busco_mafft' if args.case == 'idle' else 'tree_plot_input'
        raw = root / name
        raw.mkdir()
        for index in range(args.files):
            seed = hashlib.sha256(str(index).encode()).hexdigest().encode()
            path = raw / f'OG{index:07d}.tsv'
            path.write_bytes((seed * ((args.bytes_per_file + 63) // 64))[:args.bytes_per_file])
            path.chmod(0o640)
            os.utime(path, ns=(1700000000123456789, 1700000000123456789))
        expected = fingerprint(raw)
        baseline.convert(root, kind='genome_evolution' if args.case == 'idle' else 'hgt', directories=[name])
        fixture = root / (name + '.zip')
        samples = {'baseline': [], 'candidate': []}
        for trial in range(args.repeats + 1):
            for label in (('baseline', 'candidate') if trial % 2 == 0 else ('candidate', 'baseline')):
                support = getattr(args, label).resolve()
                command = [sys.executable, str(Path(__file__).resolve()), '--support', str(support), '--fixture', str(fixture), '--case', args.case]
                result = json.loads(subprocess.check_output(command, text=True))
                if result['fingerprint'] != expected:
                    raise RuntimeError(f'{label} changed restored contents or metadata')
                if trial:
                    samples[label].append(result)
        report = {'environment': {'python': sys.version, 'platform': platform.platform()},
                  'workload': {'case': args.case, 'files': args.files, 'bytes_per_file': args.bytes_per_file, 'repeats': args.repeats,
                               'discarded_warmup_pairs': 1, 'fixture_sha256': hashlib.sha256(fixture.read_bytes()).hexdigest()},
                  'equivalent': True, 'fingerprint': expected, 'samples': samples,
                  'source_sha256': {label: {str(path.name): hashlib.sha256(path.read_bytes()).hexdigest()
                                                   for path in [getattr(args, label) / 'orthofinder_output_store.py',
                                                                getattr(args, label).parent / 'core/gg_genome_evolution_core.sh']}
                                    for label in samples},
                  'median': {label: {'seconds': statistics.median(row['seconds'] for row in rows),
                                     'cpu_seconds': statistics.median(row['cpu_seconds'] for row in rows),
                                     'peak_rss_bytes': statistics.median(row['peak_rss_bytes'] for row in rows),
                                     'peak_child_process_rss_bytes': statistics.median(row['peak_child_process_rss_bytes'] for row in rows)}
                             for label, rows in samples.items()}}
        report['speedup'] = report['median']['baseline']['seconds'] / report['median']['candidate']['seconds']
        text = json.dumps(report, indent=2) + '\n'
        if args.output is not None:
            args.output.write_text(text)
        print(text)


if __name__ == '__main__':
    main()
