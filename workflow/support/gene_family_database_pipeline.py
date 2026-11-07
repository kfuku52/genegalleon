#!/usr/bin/env python3
"""Audit, build and record a database behind one source-content publication fence."""
from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import sys
import tempfile
from pathlib import Path

import artifact_provenance as provenance
import generate_orthogroup_database as database
from artifact_audit_runtime import AuditProgress, source_identity
from artifact_audit_runtime import signature as source_signature
from gene_family_output_store import GeneFamilyOutputStore, orthogroup_id_from_name, orthogroup_ids_from_genecount
from shared_namespace_lock import namespace_lock
from workflow_observation import observe_files, record_contract_result, report_error


class CatalogStore(GeneFamilyOutputStore):
    """A read-only enumeration scope; old raw and ZIP material stays intact."""
    def __init__(self, root, families):
        super().__init__(root)
        self.families = frozenset(families)

    def _file_names_unlocked(self, subdir, *, family_ids=None, family_from_name=None):
        selected = self.families if family_ids is None else self.families.intersection(family_ids)
        return super()._file_names_unlocked(subdir, family_ids=selected,
                                            family_from_name=orthogroup_id_from_name)


def catalog_store(root, record, mode):
    inputs = dict(provenance.parse_path_pairs(record.input, '--input'))
    parameters = dict(provenance.parse_path_pairs(record.parameter, '--parameter'))
    scope = parameters.get('family_scope')
    if scope is None:
        return GeneFamilyOutputStore(root), None
    if scope != 'selected_orthogroups_v1' or mode != 'orthogroup' or record.family_id != 'orthogroup':
        raise ValueError('Invalid selected orthogroup database scope')
    catalog = inputs.get('family_catalog')
    if catalog is None:
        raise ValueError('Selected orthogroup database requires its bound family catalog')
    path = Path(catalog)
    def identity():
        before = source_signature(path)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if source_signature(path) != before:
            raise ValueError('Family catalog changed while reading')
        return before, digest
    before = identity()
    families = orthogroup_ids_from_genecount(path)
    if not families or any(orthogroup_id_from_name(family) != family for family in families):
        raise ValueError('Family catalog must contain nonempty, exact orthogroup identities')
    if identity() != before:
        raise ValueError('Family catalog changed before selection')
    return CatalogStore(root, families), ((str(path.absolute()), 'family_catalog'), before, identity)


def validate_destinations(audit, record, destination):
    report = audit.output_tsv
    paths = [destination, record.manifest, report, report.with_suffix(report.suffix + '.progress.json'),
             report.with_suffix(report.suffix + '.result.json')]
    for index, path in enumerate(paths):
        if path.is_symlink() or path.exists() and not path.is_file():
            raise ValueError('Database publication requires regular output paths: ' + str(path))
        for other in paths[:index]:
            if path.resolve() == other.resolve() or path.exists() and other.exists() and os.path.samefile(path, other):
                raise ValueError('Database, provenance and audit output paths must be distinct')


def publish_database(build, record, payload):
    """Publish the validated pair, restoring both files on a caught failure.

    The caller holds the database and manifest writer locks. Back up the large
    closed SQLite file with a hard link; copy the small manifest so metadata
    readers never observe a multiply-linked source. Packed old manifests need
    no copy: removing a failed live override exposes the unchanged packed body.
    """
    private, final, mode, publish = build
    final, manifest = Path(final), record.manifest
    directories, backups, installed = [], {}, []
    preserve = False
    with provenance.metadata_writer(manifest):
        try:
            for target in ([final] if publish else []) + [manifest]:
                if target.is_symlink() or target.exists() and not target.is_file():
                    raise ValueError('Database publication requires regular output paths: ' + str(target))
                target.parent.mkdir(parents=True, exist_ok=True)
                directory = Path(tempfile.mkdtemp(prefix='.' + target.name + '.publish-', dir=target.parent))
                directories.append(directory)
                backups[target] = None
                if target.exists():
                    backup = directory / 'original'
                    if target == final:
                        os.link(target, backup, follow_symlinks=False)
                    else:
                        shutil.copy2(target, backup)
                    backups[target] = backup
            try:
                if publish:
                    if mode == 'create':
                        os.link(private, final)
                    else:
                        os.replace(private, final)
                    installed.append(final)
                # The writer can fail after replacing its target (e.g. during
                # lock cleanup); include the manifest in rollback before entry.
                installed.append(manifest)
                provenance.write_manifest_atomic(manifest, payload)
            except BaseException as exc:
                failures = []
                for target in reversed(installed):
                    try:
                        backup = backups[target]
                        if backup is None:
                            target.unlink(missing_ok=True)
                        else:
                            os.replace(backup, target)
                    except OSError as failure:
                        failures.append(str(target) + ': ' + str(failure))
                if failures:
                    preserve = True
                    raise RuntimeError('Database publication rollback needs recovery; backups retained in '
                                       + ', '.join(map(str, directories)) + '; ' + '; '.join(failures)) from exc
                raise
        finally:
            if not preserve:
                for directory in directories:
                    shutil.rmtree(directory)


def run(audit_argv, database_argv, record_argv):
    audit = provenance.build_parser().parse_args(['audit', *audit_argv])
    parser = database.build_parser()
    args = parser.parse_args(database_argv)
    record = provenance.build_parser().parse_args(['record', *record_argv])
    destination = Path(args.dbpath).absolute()
    root = audit.logical_root.resolve(strict=True)
    stores = provenance.parse_path_pairs(record.input_gene_family_store, '--input-gene-family-store')
    outputs = provenance.parse_path_pairs(record.output, '--output')
    if (not args.dbpath or not args.dir_gene_family or Path(args.dir_gene_family).resolve() != root
            or len(stores) != 1 or Path(stores[0][1]).resolve() != root
            or len(outputs) != 1 or Path(outputs[0][1]).absolute() != destination
            or record.optional_output or record.output_logical_directory
            or audit.workspace_root.resolve() != record.workspace_root.resolve()):
        raise ValueError('Audit, database and record must describe the same workspace, store and sole database output')
    validate_destinations(audit, record, destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    inventory = hashlib.sha256()
    build = None
    memo = provenance.AuditDigests()
    with namespace_lock(destination.with_name('.' + destination.name + '.build.guard'), exclusive=True), \
            provenance.artifact_manifest_lock(record), \
            AuditProgress(audit.output_tsv.with_suffix(audit.output_tsv.suffix + '.progress.json'),
                          interval=audit.progress_interval, attempt_dir=os.environ.get('GG_OBSERVATION_ATTEMPT_DIR'),
                          workers=audit.workers, source_sha256=source_identity(provenance.SCRIPT_DIR)) as progress:
        try:
            store, catalog_guard = catalog_store(root, record, audit.mode)
            with store.read_snapshot(), observe_files(), provenance.digest_observation(
                    workers=audit.workers, progress=progress, spill_guards=True) as memo:
                bound_inputs = {str(Path(path).resolve()) for _, path in provenance.parse_path_pairs(record.input, '--input')}
                for name, identity in database.reconciliation_input_identities(args).items():
                    if identity[0] not in bound_inputs:
                        raise ValueError('Reconciliation input missing from provenance record: ' + name)
                    memo.guard((identity[0], name), identity,
                               lambda name=name: database.reconciliation_input_identities(args)[name])
                if catalog_guard is not None:
                    memo.guard(*catalog_guard)
                provenance._collect_audit_rows(audit, memo, progress, audit.workers, store, rows, inventory, revalidate=False)
                if isinstance(store, CatalogStore):
                    if any(row['family_id'] not in store.families and row['family_id'] != '-' for row in rows):
                        raise ValueError('Manifest family differs from the selected catalog')
                    for subdir in ('stat_tree', 'stat_branch'):
                        actual = {orthogroup_id_from_name(name) for name in store.file_names(subdir)}
                        if actual != store.families:
                            raise ValueError('Selected family statistics are incomplete: ' + subdir)
                # Structural or stale-policy failures are evaluated before any DB build.
                if provenance.audit_failures(audit, rows):
                    raise ValueError('Gene-family artifact provenance audit failed')
                if all(store.file_names(subdir) for subdir in ('stat_tree', 'stat_branch')):
                    progress.phase('database_build')
                    build = database._populate_database(args, parser, store)
                    private, final, mode, publish = build
                    payload = provenance.build_contract(record, include_diagnostics=True,
                        output_sources={outputs[0][0]: str(private)})
                    if os.environ.get('GG_OBSERVATION_ATTEMPT_DIR'):
                        payload['diagnostics']['observation_attempt_id'] = Path(os.environ['GG_OBSERVATION_ATTEMPT_DIR']).name
                # Content and archive epoch fences both precede publication.
                store._close_snapshot_zip_readers()
            provenance._publish_audit(audit, memo, progress, rows, inventory)
            if build is None:
                print('Skipping database prep because required logical statistics are absent.')
                return 0
            publish_database(build, record, payload)
            if os.environ.get('GG_OBSERVATION_ATTEMPT_DIR'):
                record_contract_result(record, 0)
            print(f'Recorded audited database: {destination}')
            return 0
        except Exception as exc:
            rows.append({'family_id':'-', 'step':'snapshot_identity', 'status':'audit_error',
                         'reason':str(exc), 'manifest':''})
            provenance._publish_audit(audit, memo, progress, rows, inventory)
            if os.environ.get('GG_OBSERVATION_ATTEMPT_DIR'):
                record_contract_result(record, provenance.ERROR)
                report_error('provenance_error', step=record.step, detail=str(exc))
            raise
        finally:
            if build and build[0] != build[1]:
                Path(build[0]).unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('audit', 'database', 'record'):
        parser.add_argument('--' + name, action='append', required=True, help='one exact argv token; use --' + name + '=VALUE')
    args = parser.parse_args()
    return run(args.audit, args.database, args.record)


if __name__ == '__main__':
    sys.exit(main())
