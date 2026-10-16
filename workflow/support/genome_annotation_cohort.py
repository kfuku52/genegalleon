#!/usr/bin/env python3
"""Run independently scheduled acquisition, prediction and CDS BUSCO for a frozen cohort.

The existing entrypoints own downloading, prediction, formatting and resumption.
This host-side driver freezes task membership and hands verified physical inputs
between those stages; it never submits scheduler jobs or replaces curated inputs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

from annevo_runtime_download import atomic_json, file_sha
from artifact_provenance import audit_manifest, load_manifest, resolve_reference, runtime_support_root
from format_species_providers.common import provider_raw_dir
from gene_family_output_store import GeneFamilyOutputStore
from shared_namespace_lock import namespace_lock


def physical(path, *, directory=False):
    path = Path(path)
    if (not path.is_absolute() or path.resolve() != path or path.is_symlink()
            or not (path.is_dir() if directory else path.is_file())
            or path.stat().st_uid != os.getuid()):
        raise ValueError(f"Expected an owned physical {'directory' if directory else 'file'}: {path}")
    return path


def select(manifest, task):
    if manifest.get('schema_version') != 1 or manifest.get('busco_lineage') != 'embryophyta_odb12':
        raise ValueError('Unsupported cohort schema or BUSCO lineage')
    members = manifest['members']
    if (not members or [m['task'] for m in members] != list(range(1, len(members) + 1))
            or len({m['species'] for m in members}) != len(members)):
        raise ValueError('Cohort membership must have unique species and consecutive frozen tasks')
    for member in members:
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_.-]+', member['species']):
            raise ValueError('Invalid species label')
        if member['filename'] not in (member['species'] + '.genome.fa', member['species'] + '.genome.fa.gz'):
            raise ValueError('Genome filename must bind the frozen species and role')
    if type(task) is not int or not 1 <= task <= len(members):
        raise ValueError('Task is outside the frozen cohort')
    return members[task - 1]


def bind(source, target):
    source = physical(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    physical(target.parent, directory=True)
    if target.exists() or target.is_symlink():
        if target.is_symlink() or not target.samefile(source):
            raise ValueError(f'Conflicting native input: {target}')
        return
    os.link(source, target)


def environment(manifest, workspace):
    project = physical(manifest['project_root'], directory=True)
    cache = physical(manifest['downloads_root'], directory=True)
    runtime = physical(manifest['runtime_root'], directory=True)
    image = physical(manifest['image_path'])
    if any(not p.is_relative_to(project) for p in (workspace, cache, runtime, image)):
        raise ValueError('Cohort paths must stay in the bound project')
    env = dict(os.environ)
    # Do not let another analysis's scoped overrides leak into this phase.
    for key in list(env):
        if key.startswith(('GG_INPUT_', 'GG_GENOME_ANNOTATION_')):
            del env[key]
    env.update(gg_workspace_dir=str(workspace), gg_container_image_path=str(image),
               KFAUTO_RUNTIME_RELEASE_ROOT=str(runtime), GG_ARRAY_TASK_ID='1',
               GG_COMMON_DOWNLOADS_ROOT=str(cache), GG_COMMON_TMP_ROOT='workspace',
               GG_COMMON_BUSCO_LINEAGE=manifest['busco_lineage'], GG_COMMON_GENETIC_CODE='1',
               GG_CONTAINER_PROJECT_ROOT_BIND=f'{project}:{project}')
    return env


def entrypoint(manifest, name, env):
    runtime = Path(manifest['runtime_root'])
    subprocess.run(['bash', str(runtime / 'workflow' / (name + '_entrypoint.sh'))],
                   cwd=runtime, env=env, check=True)


def acquire(manifest, member, workspace):
    source = member['source']
    if source.get('local_path'):
        genome = physical(source['local_path'])
        if file_sha(genome) != source['sha256']:
            raise ValueError('Previously acquired genome has changed')
    else:
        row_path = workspace / 'genome.download.tsv'
        row = dict(provider='direct', id=member['species'], species_key=member['species'],
                   genome_url=source['url'], genome_filename=member['filename'],
                   genome_archive_member=source.get('archive_member', ''))
        with row_path.open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row), delimiter='\t')
            writer.writeheader()
            writer.writerow(row)
        env = environment(manifest, workspace)
        env.update(GG_INPUT_PROVIDER='direct', GG_INPUT_DOWNLOAD_MANIFEST=str(row_path),
                   GG_INPUT_DOWNLOAD_DIR=str(workspace / 'raw'), GG_INPUT_DOWNLOAD_ONLY='1',
                   GG_INPUT_DOWNLOAD_LIMIT_DIR=str(Path(manifest['downloads_root']) / 'input_download_limits'),
                   GG_INPUT_INPUT_GENERATION_MODE='single')
        for flag in ('VALIDATE_INPUTS', 'CDS_FX2TAB', 'SPECIES_BUSCO', 'MULTISPECIES_SUMMARY',
                     'SPECIES_TAXONOMY', 'GENERATE_SPECIES_TRAIT', 'GENE_MODEL_RESCUE', 'GENE_MODEL_REFINEMENT'):
            env['GG_INPUT_RUN_' + flag] = '0'
        entrypoint(manifest, 'gg_input_generation', env)
        genome = physical(provider_raw_dir('direct', workspace / 'raw', member['species']) / member['filename'])
        if source.get('md5'):
            digest = hashlib.md5(usedforsecurity=False)
            with genome.open('rb') as stream:
                for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
                    digest.update(chunk)
            if source.get('archive_member') or digest.hexdigest() != source['md5']:
                raise ValueError('Downloaded genome does not match the published MD5')
    target = workspace / 'input/species_genome' / member['filename']
    bind(genome, target)
    atomic_json(workspace / 'genome.acquired.json',
                dict(species=member['species'], source=source, path=str(target), sha256=file_sha(target)))


def audit(path, workspace, step, species):
    payload = load_manifest(path)
    if payload.get('step') != step or payload.get('family_id') != species:
        raise ValueError('Native artifact identity differs from the cohort member')
    logical = workspace / 'output/.gg_global_artifacts'
    with runtime_support_root(Path(__file__).parent):
        status, reason = audit_manifest(payload, GeneFamilyOutputStore(logical), logical, workspace)
    if status != 'current':
        raise ValueError(f'Native artifact is not current: {status}: {reason}')
    return payload


def annotate(manifest, member, workspace, phase):
    receipt = json.loads(physical(workspace / 'genome.acquired.json').read_text())
    genome = physical(workspace / 'input/species_genome' / member['filename'])
    if (receipt['source'] != member['source'] or receipt['species'] != member['species']
            or receipt['path'] != str(genome) or receipt['sha256'] != file_sha(genome)):
        raise ValueError('Acquisition receipt or genome has changed')
    species = member['species']
    prediction = workspace / 'output/genome_prediction/oriongeno' / species
    phase_workspace = workspace
    if phase == 'busco':
        audit(prediction / 'prediction.provenance.json', workspace, 'genome_annotation_oriongeno_prediction', species)
        formatted = audit(prediction / 'format.provenance.json', workspace, 'genome_annotation_oriongeno_format', species)
        entry = next(e for e in formatted['outputs'] if e['label'] == 'formatted')
        directory = resolve_reference(entry, workspace / 'output/.gg_global_artifacts', workspace)
        phase_workspace = workspace / 'busco'
        for role in ('cds', 'gff', 'genome'):
            files = [p for p in (directory / ('species_' + role)).iterdir() if p.is_file() and not p.name.startswith('.')]
            if len(files) != 1:
                raise ValueError('Expected one verified formatted file per role')
            bind(files[0], phase_workspace / 'input' / ('species_' + role) / files[0].name)
    env = environment(manifest, phase_workspace)
    for flag in ('CDSKIT_LOCALIZE', 'SCAFFOLD_TAXONOMY', 'MULTISPECIES_SUMMARY'):
        env['GG_GENOME_ANNOTATION_RUN_' + flag] = '0'
    env.update(GG_GENOME_ANNOTATION_GENOME_ANNOTATION_BACKEND='oriongeno' if phase == 'predict' else 'none',
               GG_GENOME_ANNOTATION_ORIONGENO_MODEL_PROFILE='plants',
               GG_GENOME_ANNOTATION_ORIONGENO_BATCH_SIZE='1',
               GG_GENOME_ANNOTATION_RUN_BUSCO_CDS='1' if phase == 'busco' else '0',
               GG_GENOME_ANNOTATION_BUSCO_LINEAGE=manifest['busco_lineage'])
    entrypoint(manifest, 'gg_genome_annotation', env)
    if phase == 'predict':
        audit(prediction / 'prediction.provenance.json', workspace, 'genome_annotation_oriongeno_prediction', species)
        audit(prediction / 'format.provenance.json', workspace, 'genome_annotation_oriongeno_format', species)
    else:
        proof = phase_workspace / 'output/artifact_provenance/genome_annotation' / (species + '.busco_cds.json')
        result = audit(proof, phase_workspace, 'genome_annotation_busco_cds', species)
        if result['parameters'].get('lineage') != manifest['busco_lineage']:
            raise ValueError('BUSCO artifact has the wrong lineage')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--task', type=int, required=True)
    parser.add_argument('--phase', choices=('acquire', 'predict', 'busco'), required=True)
    args = parser.parse_args()
    path = physical(args.manifest)
    if file_sha(path) != args.manifest_sha256:
        parser.error('Cohort manifest changed')
    manifest = json.loads(path.read_text())
    member = select(manifest, args.task)
    work_root = physical(manifest['work_root'], directory=True)
    workspace = work_root / member['species']
    environment(manifest, workspace)
    workspace.mkdir(parents=True, exist_ok=True, mode=0o700)
    physical(workspace, directory=True)
    with namespace_lock(workspace / '.cohort.lock', exclusive=True, timeout=60):
        if args.phase == 'acquire':
            acquire(manifest, member, workspace)
        else:
            annotate(manifest, member, workspace, args.phase)


if __name__ == '__main__':
    main()
