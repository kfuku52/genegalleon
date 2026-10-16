#!/usr/bin/env python3
"""Prepare and verify persistent public reads with the native transcriptome workflow.

Metadata stays with its owning project; this route stops before assembly,
BUSCO and expression analysis. A separately scheduled assembly can consume the
verified cache after an explicit data handoff.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

from annevo_runtime_download import atomic_json, file_sha
from artifact_provenance import ProvenanceError
from genome_annotation_cohort import audit, entrypoint, environment, physical
from shared_namespace_lock import namespace_lock

INPUT_MODE = 'metadata-read-preparation'
DISABLED = ('ASSEMBLY', 'LONGESTCDS', 'LONGESTCDS_FX2TAB', 'LONGESTCDS_MMSEQS2TAXONOMY',
            'LONGESTCDS_CONTAMINATION_REMOVAL', 'BUSCO_ISOFORMS', 'BUSCO_LONGEST_CDS',
            'BUSCO_CONTAMINATION_REMOVED_LONGEST_CDS', 'ASSEMBLY_STAT', 'AMALGKIT_QUANT',
            'AMALGKIT_MERGE', 'MULTISPECIES_SUMMARY', 'AMALGKIT_METADATA_OR_INTEGRATE')


def inspect_metadata(path, expected_sha, expected_runs):
    path = physical(path)
    if file_sha(path) != expected_sha:
        raise ValueError('Canonical selected metadata has changed')
    with path.open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    runs = [row['run'] for row in rows]
    if not runs or len(runs) != len(set(runs)) or sorted(runs) != sorted(expected_runs):
        raise ValueError('Selected metadata RUN membership has changed')
    return dict(metadata_sha256=expected_sha, run_ids=sorted(runs), files=[path.name])


def native_environment(manifest, *, generate):
    workspace = physical(manifest['workspace_root'], directory=True)
    native = {**manifest, 'work_root': str(workspace),
              'runtime_root': manifest['runtime_release_root'],
              'runtime_manifest_sha256': manifest['runtime_release_manifest_sha256']}
    env = environment(native, workspace)
    for key in list(env):
        if key.startswith('GG_TRANSCRIPTOME_'):
            del env[key]
    for flag in DISABLED:
        env['GG_TRANSCRIPTOME_RUN_' + flag] = '0'
    env.update(GG_TRANSCRIPTOME_MODE_TRANSCRIPTOME_ASSEMBLY='metadata',
               GG_TRANSCRIPTOME_RUN_AMALGKIT_GETFASTQ='1' if generate else '0',
               GG_TRANSCRIPTOME_REMOVE_AMALGKIT_FASTQ_AFTER_COMPLETION='0',
               GG_TRANSCRIPTOME_REPRESENTATIVE_INPUTS='', GG_COMMON_REPRESENTATIVE_INPUTS='',
               GG_TRANSCRIPTOME_TRANSCRIPTOME_GETFASTQ_CACHE_DIR=manifest['fastq_source_root'],
               GG_TRANSCRIPTOME_GETFASTQ_CACHE_DIR=manifest['fastq_source_root'])
    return native, env


def container_preflight(manifest, row, manifest_sha):
    metadata = Path(manifest['workspace_root']) / 'input/amalgkit_metadata' / (row['species'] + '_metadata.tsv')
    source = row['files'][0]
    inspect_metadata(metadata, source['sha256'], row['runs'])
    if file_sha(physical(manifest['image_path'])) != manifest['image_sha256']:
        raise ValueError('Read preparation SIF has changed')
    command = container_command(manifest) + ['python', '-B',
               '/script/support/transcriptome_read_preparation.py', '--inspect-json',
               json.dumps(dict(path=str(metadata), sha256=source['sha256'], runs=row['runs']))]
    result = subprocess.run(command, check=True, capture_output=True, text=True, timeout=120)
    actual = json.loads(result.stdout)
    expected = inspect_metadata(metadata, source['sha256'], row['runs'])
    if actual != expected:
        raise ValueError('Container metadata view differs from the frozen input')
    core = 'runtime/workflow/core/gg_transcriptome_generation_core.sh'
    return dict(schema_version=2, input_mode=INPUT_MODE, manifest_sha256=manifest_sha,
                image_sha256=manifest['image_sha256'], core_sha256=manifest['managed_files'][core]['sha256'],
                task_id=row['task_id'], species=row['species'], **actual)


def container_command(manifest):
    root = Path(manifest['root'])
    return [manifest['container_binary'], 'exec', '--bind', manifest['project_root'] + ':' + manifest['project_root'],
            '--bind', str(root / 'runtime/workflow') + ':/script:ro', '--bind', manifest['workspace_root'] + ':/workspace',
            '--bind', manifest['fastq_source_root'] + ':/gg_transcriptome_getfastq_cache', manifest['image_path']]


def verify_native(manifest, row):
    workspace = Path(manifest['workspace_root'])
    proof = physical(workspace / 'output/transcriptome_assembly/artifact_provenance' / (row['species'] + '.getfastq.json'))
    completion = physical(Path(manifest['fastq_source_root']) / row['species'] / 'getfastq_completion.json')
    expected = dict(native_provenance_sha256=file_sha(proof), completion_manifest_sha256=file_sha(completion))
    command = container_command(manifest) + ['python', '-B', '/script/support/transcriptome_read_preparation.py',
                                           '--verify-json', json.dumps(dict(species=row['species']))]
    result = subprocess.run(command, check=True, capture_output=True, text=True, timeout=120)
    if json.loads(result.stdout) != expected or file_sha(proof) != expected['native_provenance_sha256'] or file_sha(completion) != expected['completion_manifest_sha256']:
        raise ValueError('Native read preparation changed during verification')
    return expected


def prepare(manifest, row):
    root = Path(manifest['root'])
    native, env = native_environment(manifest, generate=True)
    env['GG_ARRAY_TASK_ID'] = str(row['task_id'])
    entrypoint(native, 'gg_transcriptome_generation', env)
    native, env = native_environment(manifest, generate=False)
    env['GG_ARRAY_TASK_ID'] = str(row['task_id'])
    entrypoint(native, 'gg_transcriptome_generation', env)
    proof = verify_native(manifest, row)
    completion = physical(Path(manifest['fastq_source_root']) / row['species'] / 'getfastq_completion.json')
    return dict(status='read-preparation-completed', species=row['species'], runs=row['runs'],
                metadata_sha256=row['files'][0]['sha256'], completion_manifest=str(completion),
                **proof,
                root=str(root), assembly_complete=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package-root', type=Path)
    parser.add_argument('--manifest-sha256')
    parser.add_argument('--task-id', type=int)
    parser.add_argument('--preflight-only', action='store_true')
    parser.add_argument('--inspect-json')
    parser.add_argument('--verify-json')
    args = parser.parse_args()
    if args.verify_json:
        species = json.loads(args.verify_json)['species']
        workspace = Path('/workspace')
        proof = workspace / 'output/transcriptome_assembly/artifact_provenance' / (species + '.getfastq.json')
        audit(proof, workspace, 'transcriptome_getfastq', species)
        completion = physical(Path('/gg_transcriptome_getfastq_cache') / species / 'getfastq_completion.json')
        print(json.dumps(dict(native_provenance_sha256=file_sha(proof), completion_manifest_sha256=file_sha(completion)), sort_keys=True))
        return
    if args.inspect_json:
        value = json.loads(args.inspect_json)
        print(json.dumps(inspect_metadata(value['path'], value['sha256'], value['runs']), sort_keys=True))
        return
    if args.package_root is None or args.task_id is None or args.manifest_sha256 is None:
        parser.error('A frozen package, manifest hash and task index are required')
    root = physical(args.package_root, directory=True)
    path = physical(root / 'manifest.json')
    if file_sha(path) != args.manifest_sha256:
        parser.error('Read preparation manifest changed')
    manifest = json.loads(path.read_text())
    if manifest['root'] != str(root) or manifest['input_mode'] != INPUT_MODE or manifest['busco_lineage'] != 'embryophyta_odb12':
        parser.error('Read preparation package identity changed')
    rows = manifest['species']
    if [r['task_id'] for r in rows] != list(range(1, len(rows) + 1)) or len({r['species'] for r in rows}) != len(rows):
        parser.error('Frozen species mapping changed')
    row = next(r for r in rows if r['task_id'] == args.task_id)
    for name, value in manifest['managed_files'].items():
        relative = Path(name)
        if relative.is_absolute() or '..' in relative.parts or file_sha(physical(root / name)) != value['sha256']:
            parser.error('Native runtime changed: ' + name)
    preflight = container_preflight(manifest, row, args.manifest_sha256)
    if args.preflight_only:
        print(json.dumps(dict(status='preflight-valid', container_input_preflight=preflight), sort_keys=True))
        return
    if not os.environ.get('SLURM_JOB_ID', '').isdigit() or os.environ.get('SLURM_ARRAY_TASK_ID') != str(args.task_id):
        parser.error('An exact scheduler task allocation is required')
    receipt = root / 'receipts' / (row['species'] + '.json')
    with namespace_lock(root / 'task-locks' / (row['species'] + '.lock'), exclusive=True, timeout=1):
        if receipt.exists() or receipt.is_symlink():
            parser.error('A prior preparation attempt exists; reconcile it first')
        result = prepare(manifest, row)
        atomic_json(receipt, {**result, 'manifest_sha256': args.manifest_sha256,
                             'job_id': os.environ.get('SLURM_ARRAY_JOB_ID', os.environ['SLURM_JOB_ID']),
                             'task_id': args.task_id, 'container_input_preflight': preflight})


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, ProvenanceError, subprocess.SubprocessError) as exc:
        print('Read preparation: ' + str(exc), file=sys.stderr)
        raise SystemExit(1) from exc
