#!/usr/bin/env python3
"""Compare input-plan/staging reads and exact receipts with a baseline checkout."""
import argparse
import contextlib
import gzip
import hashlib
import io
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


def worker(args):
    sys.path.insert(0, str(args.support_root))
    import input_generation_array_state as state
    import plan_input_generation_tasks as planner
    import run_input_generation_task as runner
    import stage_input_generation_downloads as staging

    root = args.fixture_root
    plan = root / 'plan.json'
    plan.unlink(missing_ok=True)
    shutil.rmtree(str(plan) + '.tasks', ignore_errors=True)
    roles = {'cds': root / 'cds.fa', 'gff': root / 'annotation.gff.gz', 'genome': root / 'genome.fa'}
    if args.alias_roles:
        roles['cds'] = roles['genome']
    manifest = root / 'manifest.tsv'
    manifest.write_text('provider\tid\tspecies_key\tbind_local_sources\tcds_url\tgff_url\tgenome_url\n'
                        + 'direct\tfixture\tExample_species\t1\t'
                        + '\t'.join(roles[role].as_uri() for role in ('cds', 'gff', 'genome')) + '\n')
    if args.staged_reuse:
        donor = root / 'donor/output/input_generation/tmp/task_plan.json'
        with contextlib.redirect_stdout(io.StringIO()):
            sys.argv = ['planner', '--provider', 'all', '--download-manifest', str(manifest),
                        '--download-dir', str(root / 'downloads'), '--stage-downloads', '--outfile', str(donor)]
            assert planner.main() == 0
            staging.stage_downloads(donor)
        settings = Path(str(donor) + '.settings.json')
        settings.write_text('{}\n')
        state.atomic_json(Path(str(donor) + '.prepared.json'), {
            'plan_sha256': state.digest(donor), 'settings_sha256': state.digest(settings),
            'files': state.digest_paths(Path(str(donor) + '.tasks').glob('1.*')),
        })
    original = state.digest
    reads = []

    def counted_digest(path):
        reads.append(Path(path).stat().st_size)
        return original(path)

    state.digest = staging.digest = runner.digest = counted_digest
    if hasattr(planner, 'digest'):
        planner.digest = counted_digest
    phases = {}
    for phase in ('plan', 'stage', 'resume', 'describe'):
        reads.clear()
        started = time.perf_counter()
        with contextlib.redirect_stdout(io.StringIO()):
            if phase == 'plan':
                if args.staged_reuse:
                    state.export_staged_manifest(donor, manifest)
                sys.argv = ['planner', '--provider', 'all', '--download-manifest', str(manifest),
                            '--download-dir', str(root / 'downloads'), '--stage-downloads', '--outfile', str(plan)]
                assert planner.main() == 0
            elif phase == 'describe':
                sys.argv = ['runner', '--task-plan', str(plan), '--task-index', '1', '--describe-only',
                            '--species-cds-dir', str(root / 'cds'), '--species-gff-dir', str(root / 'gff'),
                            '--species-genome-dir', str(root / 'genomes'), '--task-meta-output', str(root / 'metadata.json')]
                assert runner.main() == 0
            else:
                staging.stage_downloads(plan, require_gff=True, require_genome=True)
        phases[phase] = {'seconds': time.perf_counter() - started,
                         'sha256_calls': len(reads), 'sha256_bytes': sum(reads)}
    proof = hashlib.sha256()
    if args.compare_staged_reuse:
        # The new plan adds sealed staging evidence. Compare all staged input
        # roles, identities, parameters and source hashes, excluding that added
        # evidence field; no scientific input/result is excluded.
        actual = json.loads((Path(str(plan) + '.tasks') / '1.json').read_text())['task']
        actual.pop('staged_input_reuse', None)
        proof.update(json.dumps(actual, sort_keys=True).encode())
    else:
        for path in [plan, root / 'metadata.json', *sorted(Path(str(plan) + '.tasks').iterdir())]:
            proof.update(path.name.encode())
            proof.update(path.read_bytes())
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(json.dumps({'phases': phases, 'receipt_sha256': proof.hexdigest(),
                      'peak_rss_bytes': peak_rss * (1 if sys.platform == 'darwin' else 1024)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--baseline-support', type=Path)
    parser.add_argument('--genome-mib', type=int, default=256)
    parser.add_argument('--trials', type=int, default=3)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--alias-roles', action='store_true')
    parser.add_argument('--compare-staged-reuse', action='store_true', help='Compare ordinary bound inputs with native sealed staging export.')
    parser.add_argument('--staged-reuse', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--fixture-root', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not 1 <= args.genome_mib <= 4096 or not 1 <= args.trials <= 10:
        parser.error('Require 1..4096 genome MiB and 1..10 trials')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    roots = {'current': args.support_root.resolve()}
    if args.compare_staged_reuse:
        roots = {'baseline': args.support_root.resolve(), 'current': args.support_root.resolve()}
    if args.baseline_support:
        roots['baseline'] = args.baseline_support.resolve()
    result = {'python': sys.version, 'platform': platform.platform(), 'genome_mib': args.genome_mib,
              'alias_roles': args.alias_roles, 'cache_mode': 'warm; no filesystem cache eviction', 'samples': {}}
    with tempfile.TemporaryDirectory(prefix='gg-staging-benchmark-') as temporary:
        root = Path(temporary)
        with (root / 'genome.fa').open('wb') as handle:
            handle.write(b'>chr1\n')
            chunk = (b'ATG' * 1024 + b'\n') * 341
            for _ in range(args.genome_mib):
                handle.write(chunk)
        (root / 'cds.fa').write_bytes(b'>gene1\n' + b'ATG' * 1024 * 341 + b'\n')
        (root / 'annotation.gff.gz').write_bytes(gzip.compress(
            b'##gff-version 3\nchr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1\n', mtime=0))
        result['raw_bytes'] = sum(path.stat().st_size for path in root.iterdir())
        for trial in range(args.trials + 1):
            items = list(roots.items())
            if trial % 2:
                items.reverse()
            for label, support_root in items:
                command = [sys.executable, str(Path(__file__).resolve()), '--worker', '--support-root', str(support_root),
                           '--fixture-root', str(root)]
                if args.alias_roles:
                    command.append('--alias-roles')
                if args.compare_staged_reuse:
                    command.append('--compare-staged-reuse')
                    if label == 'current':
                        command.append('--staged-reuse')
                sample = json.loads(subprocess.check_output(command, text=True))
                print(label, 'warmup' if trial == 0 else trial, sample['phases'], flush=True)
                if trial:
                    result['samples'].setdefault(label, []).append(sample)
    assert len({sample['receipt_sha256'] for samples in result['samples'].values() for sample in samples}) == 1
    result['medians'] = {}
    for label, samples in result['samples'].items():
        result['medians'][label] = {
            phase: {metric: statistics.median(sample['phases'][phase][metric] for sample in samples)
                    for metric in ('seconds', 'sha256_calls', 'sha256_bytes')}
            for phase in ('plan', 'stage', 'resume', 'describe')}
        result['medians'][label]['peak_rss_bytes'] = statistics.median(
            sample['peak_rss_bytes'] for sample in samples)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('receipt_sha256', next(iter(result['samples'].values()))[0]['receipt_sha256'])
    print('medians', result['medians'])


if __name__ == '__main__':
    main()
