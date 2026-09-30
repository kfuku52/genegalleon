#!/usr/bin/env python3
"""Audit, build and record a database behind one source-content publication fence."""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
from pathlib import Path

import artifact_provenance as provenance
import generate_orthogroup_database as database
from artifact_audit_runtime import AuditProgress, source_identity
from gene_family_output_store import GeneFamilyOutputStore
from shared_namespace_lock import namespace_lock
from workflow_observation import observe_files, record_contract_result, report_error


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
    destination.parent.mkdir(parents=True, exist_ok=True)
    store = GeneFamilyOutputStore(root)
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
            with store.read_snapshot(), observe_files(), provenance.digest_observation(workers=audit.workers, progress=progress) as memo:
                provenance._collect_audit_rows(audit, memo, progress, audit.workers, store, rows, inventory, revalidate=False)
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
            provenance._publish_audit(audit, memo, progress, rows, inventory)
            if build is None:
                print('Skipping database prep because required logical statistics are absent.')
                return 0
            if publish:
                if mode == 'create':
                    os.link(private, final)
                    os.unlink(private)
                else:
                    os.replace(private, final)
            provenance.write_manifest_atomic(record.manifest, payload)
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
