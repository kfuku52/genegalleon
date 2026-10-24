#!/usr/bin/env python3
"""Finalize a selected cohort of native formatted inputs without rewriting IDs."""

import argparse
import collections
import contextlib
import csv
import errno
import hashlib
import io
import json
import os
import shutil
import stat
from pathlib import Path

from input_generation_array_state import (
    FreshDigestBatch,
    atomic_json,
    digest,
    namespace_path,
    receipt_path,
    safe_component,
)
from shared_namespace_lock import namespace_lock

KIND = 'native-formatted-input-cohort-v1'
OUTPUTS = {'cds': 'species_cds', 'gff': 'species_gff', 'genome': 'species_genome'}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read_table(path):
    with Path(path).open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        rows = list(reader)
    require(rows and all(None not in r and all(v is not None for v in r.values()) for r in rows), 'Malformed/empty table')
    return rows


def write_table(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name('.' + path.name + '.tmp-' + str(os.getpid()))
    temporary.write_text(stream.getvalue())
    os.replace(temporary, path)


def freeze(selection, expected_species, *, source_aliases=None):
    """Seal native donor metadata; large-file verification belongs to prepare/QC."""
    rows = read_table(selection)
    require(len(rows) == expected_species, 'Selected species count differs')
    species, tasks, documents = set(), [], {}
    source_aliases, used_aliases = dict(source_aliases or {}), set()
    for old, new in source_aliases.items():
        require(all(isinstance(p, str) and Path(p).is_absolute() and '..' not in Path(p).parts for p in (old, new))
                and old != new and Path(new).is_file() and not Path(new).is_symlink(), 'Invalid raw source relocation')
    def document(path):
        path = Path(path)
        if str(path) not in documents:
            before = path.stat()
            data = path.read_bytes()
            after = path.stat()
            require((before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns) ==
                    (after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns), 'Donor metadata changed')
            documents[str(path)] = (json.loads(data), hashlib.sha256(data).hexdigest())
        return documents[str(path)]
    for row in rows:
        sp = row['species']
        require(safe_component(sp) and sp not in species, 'Unsafe/duplicate selected species')
        species.add(sp)
        workspace, plan_path = Path(row['source_workspace']), Path(row['plan_path'])
        require(workspace.is_absolute() and plan_path == workspace / 'output/input_generation/tmp/task_plan.json', 'Donor plan/root differ')
        plan, plan_sha = document(plan_path)
        require(plan_sha == row['plan_sha256'], 'Selected donor plan hash changed: ' + sp)
        matching = [i for i, task in enumerate(plan['tasks'], 1) if task['species_prefix'] == sp]
        require(len(matching) == 1 and len(plan['tasks']) == plan['task_count'], 'Donor species is ambiguous')
        index = matching[0]
        receipt = receipt_path(plan_path, index)
        proof, receipt_sha = document(receipt)
        require(proof['plan_sha256'] == plan_sha and proof['task_index'] == index and proof['species_prefix'] == sp,
                'Native worker receipt differs: ' + sp)
        expected = {}
        for name, sha in proof['files'].items():
            path = str(namespace_path(name, workspace))
            require(Path(path).is_absolute() and '..' not in Path(path).parts and
                    (path not in expected or expected[path] == sha), 'Unsafe/conflicting native receipt member')
            expected[path] = sha
        shard = workspace / 'output/input_generation/tmp/species_summary_shards' / (str(index) + '.tsv')
        require(str(shard) in expected and digest(shard) == expected[str(shard)], 'Native summary shard is unverified')
        summaries = read_table(shard)
        require(len(summaries) == 1 and summaries[0]['species_prefix'] == sp, 'Native summary species differs')
        summary = {key: str(namespace_path(value, workspace)) if value.startswith('/workspace/') else value
                   for key, value in summaries[0].items()}
        for role in OUTPUTS:
            value = summary.get(role + '_output_path', '')
            require((not value and not row.get(role + '_path')) or
                    (value == row.get(role + '_path') and value in expected), 'Selected formatted ' + role + ' differs: ' + sp)
        require(summary['cds_output_path'], 'Selected CDS missing: ' + sp)
        audit = summary.get('cloned_model_audit_path', '')
        require(not row.get('cloned_model_audit_path') or row['cloned_model_audit_path'] == audit, 'Selected clone audit differs')
        if audit:
            from cloned_cds_models import evidence_paths
            for path in evidence_paths(summary, summary=True):
                require(path in expected, 'Clone evidence absent from native receipt: ' + path)
        expected.update({str(plan_path): plan_sha, str(receipt): receipt_sha})
        require(all(isinstance(sha, str) and len(sha) == 64 and all(c in '0123456789abcdef' for c in sha)
                    for sha in expected.values()), 'Invalid native content hash')
        task = dict(species=sp, source_summary=summaries[0], source_workspace=str(workspace), plan_path=str(plan_path),
                          plan_sha256=plan_sha, task_index=index, receipt_path=str(receipt), receipt_sha256=receipt_sha,
                          summary=summary, files=expected)
        relocated, qc_summary = {}, dict(summaries[0])
        for field in ('cds_input_path', 'gff_input_path', 'genome_input_path'):
            old = summary.get(field, '')
            if old not in source_aliases:
                continue
            new = source_aliases[old]
            require(old in expected and digest(new) == expected[old], 'Relocated raw source content differs')
            require(new not in expected or expected[new] == expected[old], 'Conflicting relocated source hash')
            relocated[old] = dict(path=new, sha256=expected[old])
            expected[new] = expected.pop(old)
            summary[field] = qc_summary[field] = new
            used_aliases.add(old)
        if relocated:
            task.update(source_relocations=relocated, qc_summary=qc_summary)
        tasks.append(task)
    require(used_aliases == set(source_aliases), 'Unused/unverified raw source relocation')
    return dict(kind=KIND, schema_version=1, task_count=len(tasks), selection_sha256=digest(selection), tasks=tasks)


def load_cohort(path):
    path = Path(path)
    doc = json.loads(path.read_text())
    require(doc['kind'] == KIND and doc['schema_version'] == 1 and doc['task_count'] == len(doc['tasks']) > 0,
            'Invalid cohort plan')
    require(doc['task_count'] <= 10000 and all(safe_component(t['species']) for t in doc['tasks']) and
            len({t['species'] for t in doc['tasks']}) == doc['task_count'], 'Duplicate/unsafe cohort species')
    return doc


def verify_tasks(tasks, *, batch=None):
    batch = batch or FreshDigestBatch()
    expected = {}
    for task in tasks:
        for name, sha in task['files'].items():
            require(name not in expected or expected[name] == sha, 'Conflicting donor content hashes')
            expected[name] = sha
    batch.read(expected, expected=expected)
    batch.check()
    return batch


def links_for(task):
    row, sp = task['summary'], task['species']
    links = {}
    for role, directory in OUTPUTS.items():
        path = row.get(role + '_output_path')
        if path:
            require(Path(path).is_file() and not Path(path).is_symlink(), 'Generated donor output is not regular')
            links['input/' + directory + '/' + Path(path).name] = path
            links['output/input_generation/' + directory + '/' + Path(path).name] = path
    base = Path(task['source_workspace']) / 'output/input_generation'
    for directory, filename in [('species_cds_fx2tab', sp + '_fx2tab_cds.tsv'),
                                ('species_cds_busco_full', sp + '.busco.full.tsv'),
                                ('species_cds_busco_short', sp + '.busco.short.txt'),
                                ('species_cds_busco_full/single_copy', sp + '.json.gz')]:
        path = str(base / directory / filename)
        if path in task['files']:
            require(Path(path).is_file() and not Path(path).is_symlink(), 'Generated donor output is not regular')
            links['output/input_generation/' + directory + '/' + filename] = path
    for step in ('busco', 'fx2tab'):
        path = str(base / 'artifact_provenance' / (step + '.' + sp + '.json'))
        if path in task['files']:
            links['output/input_generation/artifact_provenance/' + Path(path).name] = path
    return links


def verify_view(workspace, plan, *, indices=None, batch=None):
    batch = batch or FreshDigestBatch()
    doc = load_cohort(plan)
    chosen = doc['tasks'] if indices is None else [doc['tasks'][i - 1] for i in indices]
    proof = workspace / 'cohort_view.json'
    saved = json.loads(proof.read_text())
    require(saved['kind'] == KIND and saved['plan_sha256'] == digest(plan) and saved['task_count'] == doc['task_count'],
            'Cohort view identity changed')
    expected_links = {k: v for t in doc['tasks'] for k, v in links_for(t).items()}
    require(saved['links'] == expected_links, 'Cohort view membership changed')
    for task in chosen:
        require(saved['summaries'][task['species']] == task['summary'], 'Cohort source summary changed')
    verify_prepared_files(workspace, doc, chosen, expected_links, batch=batch)
    return doc


def verify_prepared_files(workspace, doc, chosen, expected_links, *, batch=None):
    """Verify sources, copied/linked outputs and membership before publication."""
    batch = batch or FreshDigestBatch()
    require([{k: v for k, v in r.items() if v} for r in read_table(workspace / 'output/input_generation/gg_input_generation_species.tsv')] ==
            [{k: v for k, v in t['summary'].items() if v} for t in doc['tasks']], 'Cohort summary changed')
    for role in OUTPUTS.values():
        for base in ['input/', 'output/input_generation/']:
            directory = workspace / (base + role)
            if directory.exists():
                require(not directory.is_symlink() and {str(p.relative_to(workspace)) for p in directory.iterdir()} ==
                        {k for k in expected_links if Path(k).parent == Path(base + role)}, 'Unexpected cohort species files')
    view_hashes = {}
    for task in chosen:
        for relative, source in links_for(task).items():
            target = workspace / relative
            require(target.is_file() and not target.is_symlink(), 'Cohort view file changed: ' + str(target))
            view_hashes[str(target)] = task['files'][source]
        audit = task['summary'].get('cloned_model_audit_path')
        if audit:
            snapshot = workspace / 'output/input_generation/cloned_model_history' / (task['species'] + '.json')
            require(snapshot.is_file() and not snapshot.is_symlink(), 'Cohort clone history changed')
            view_hashes[str(snapshot)] = task['files'][audit]
    verify_tasks(chosen, batch=batch)
    batch.read(view_hashes, expected=view_hashes)
    batch.check()


def prepare(workspace, plan):
    """Make a byte-preserving input view; donors remain required and unchanged."""
    workspace, plan = Path(workspace), Path(plan)
    batch = FreshDigestBatch()
    plan_sha = batch.read([plan])[str(plan)]
    doc = load_cohort(plan)
    batch.check()
    require(workspace.is_absolute() and workspace.resolve() == workspace and not workspace.is_symlink(), 'Unsafe cohort workspace')
    require(all(not workspace.is_relative_to(t['source_workspace']) and
                not Path(t['source_workspace']).is_relative_to(workspace) for t in doc['tasks']), 'Cohort overlaps a donor')
    workspace.mkdir(parents=True, exist_ok=True)
    with namespace_lock(workspace / '.cohort.lock', exclusive=True, nonblocking=True) as acquired:
        require(acquired, 'Cohort has an active owner')
        if (workspace / 'cohort_view.json').exists():
            verify_view(workspace, plan, batch=batch)
            return
        for path in workspace.rglob('*'):
            require(not path.is_symlink(), 'Unsafe existing cohort symlink')
        # Hold every donor namespace while staging the unpublished view. All
        # sources and outputs must pass the fresh final content boundary before
        # publishing cohort_view.json, including cross-filesystem copies.
        with contextlib.ExitStack() as stack:
            for root in sorted({t['source_workspace'] for t in doc['tasks']}):
                acquired = stack.enter_context(namespace_lock(Path(root) / 'output/input_generation/.array-phase.lock',
                                                             exclusive=True, nonblocking=True))
                require(acquired, 'Native donor has an active owner: ' + root)
            all_links = {}
            for task in doc['tasks']:
                for relative, source in links_for(task).items():
                    require(relative not in all_links, 'Duplicate cohort destination')
                    all_links[relative] = source
            for relative, source in all_links.items():
                target = workspace / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                if target.exists():
                    require(stat.S_ISREG(target.lstat().st_mode), 'Partial cohort file is not regular')
                else:
                    try:
                        os.link(source, target, follow_symlinks=False)
                    except OSError as exc:
                        if exc.errno not in {errno.EXDEV, errno.EPERM, errno.EACCES, errno.EMLINK}:
                            raise
                        with Path(source).open('rb') as src, target.open('xb') as dst:
                            shutil.copyfileobj(src, dst, 8 * 1024 * 1024)
            summaries = [t['summary'] for t in doc['tasks']]
            write_table(workspace / 'output/input_generation/gg_input_generation_species.tsv', summaries)
            history = workspace / 'output/input_generation/cloned_model_history'
            for task in doc['tasks']:
                source = task['summary'].get('cloned_model_audit_path')
                if source:
                    history.mkdir(parents=True, exist_ok=True)
                    target = history / (task['species'] + '.json')
                    data = Path(source).read_bytes()
                    require(hashlib.sha256(data).hexdigest() == task['files'][source], 'Clone history changed')
                    target.write_bytes(data)
            # Link creation changes donor ctime. Hash every physical source and
            # output once after staging; never reuse a digest across writes.
            verify_prepared_files(workspace, doc, doc['tasks'], all_links, batch=batch)
            atomic_json(workspace / 'cohort_view.json', dict(kind=KIND, plan_sha256=plan_sha, task_count=doc['task_count'],
                        summaries={t['species']: t['summary'] for t in doc['tasks']}, links=all_links,
                        donor_retention_required=True, sequence_identifiers_rewritten=False), immutable=True)


def qc_one(workspace, plan, index):
    from busco_quality_metadata import read_short_summary
    from validate_cds_gff_mapping import validate_single_species as mapping
    from validate_longest_cds_selection import validate_single_species as ownership
    workspace, plan = Path(workspace), Path(plan)
    doc = load_cohort(plan)
    require(1 <= index <= doc['task_count'], 'Cohort task index is out of bounds')
    task = doc['tasks'][index - 1]
    lock = workspace / 'qc' / (str(index) + '.lock')
    with namespace_lock(workspace / '.cohort.lock', exclusive=False, timeout=30) as reader, \
            namespace_lock(lock, exclusive=True, nonblocking=True) as acquired:
        require(reader and acquired, 'Cohort QC has an active owner')
        batch = FreshDigestBatch()
        verify_view(workspace, plan, indices=[index], batch=batch)
        row, sp = task.get('qc_summary', task['source_summary']), task['species']
        if any(value.startswith('/workspace/') for value in row.values()):
            require(Path('/workspace').is_dir() and Path('/workspace').samefile(task['source_workspace']),
                    'Native donor /workspace namespace is required for QC')
        shared = dict(index=index, species_prefix=sp, cds_file=Path(row['cds_output_path']), summary_row=row,
                      revalidate_historical_orfs=True)
        result = ownership(shared, 10)
        mapping_result = None
        if row.get('gff_output_path'):
            mapping_result = mapping({**shared, 'gff_file': Path(row['gff_output_path']),
                                      'genome_file': Path(row['genome_output_path']) if row.get('genome_output_path') else None,
                                      'strict': True}, 10)
        quality_path = Path(task['source_workspace']) / 'output/input_generation/species_cds_busco_short' / (sp + '.busco.short.txt')
        quality = read_short_summary(quality_path) if str(quality_path) in task['files'] else None
        full = Path(task['source_workspace']) / 'output/input_generation/species_cds_busco_full' / (sp + '.busco.full.tsv')
        statuses = {}
        if quality is not None:
            require(str(full) in task['files'], 'BUSCO full table is not receipt bound')
            groups = collections.defaultdict(list)
            for line in full.read_text().splitlines():
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.split('\t')
                require(len(fields) >= 2 and fields[1] in {'Complete', 'Duplicated', 'Fragmented', 'Missing'}, 'Invalid BUSCO status')
                groups[fields[0]].append(fields[1])
            require(groups and all(len(set(v)) == 1 for v in groups.values()), 'Contradictory BUSCO statuses')
            labels = dict(Complete='S', Duplicated='D', Fragmented='F', Missing='M')
            statuses = {k: ('?' if v[0] == 'Complete' and len(v) != 1 else labels[v[0]]) for k, v in groups.items()}
        record = dict(kind=KIND, plan_sha256=digest(plan), index=index, species=sp, ownership=result, mapping=mapping_result,
                      busco=quality, marker_statuses=statuses, passed=result.get('ok') is True and (mapping_result is None or mapping_result.get('ok') is True))
        # Native QC only reads input files. Keep this one fresh content
        # boundary fenced until publication instead of hashing it a second time.
        batch.check()
        atomic_json(workspace / 'qc' / (str(index) + '.json'), record)
        require(record['passed'], 'Native cohort QC failed: ' + sp)
        return record


def finalize(workspace, plan):
    workspace, plan = Path(workspace), Path(plan)
    with namespace_lock(workspace / '.cohort.lock', exclusive=True, nonblocking=True) as acquired:
        require(acquired, 'Cohort QC workers are still active')
        doc = verify_view(workspace, plan)
        results = []
        for index, task in enumerate(doc['tasks'], 1):
            value = json.loads((workspace / 'qc' / (str(index) + '.json')).read_text())
            require(value['kind'] == KIND and value['plan_sha256'] == digest(plan) and value['index'] == index
                    and value['species'] == task['species'] and value['passed'] is True, 'Missing/failed/mismatched native QC')
            results.append(value)
        report = dict(kind=KIND, plan_sha256=digest(plan), species=doc['task_count'], passed=True,
                      species_with_gff=sum(bool(t['summary'].get('gff_output_path')) for t in doc['tasks']),
                      species_with_genome=sum(bool(t['summary'].get('genome_output_path')) for t in doc['tasks']),
                      cloned_audits=sum(bool(t['summary'].get('cloned_model_audit_path')) for t in doc['tasks']),
                      results=results, sequence_identifiers_rewritten=False, donor_retention_required=True)
        marker_sets = {tuple(sorted(r['marker_statuses'])) for r in results if r['marker_statuses']}
        lineages = {r['busco']['lineage'] for r in results if r['busco']}
        require(len(marker_sets) <= 1 and len(lineages) <= 1, 'Cohort BUSCO marker/lineage inventories differ')
        markers = next(iter(marker_sets), ())
        write_table(workspace / 'marker_occupancy.tsv', [dict(species=r['species'], **{m: r['marker_statuses'].get(m, '.') for m in markers}) for r in results])
        report.update(busco_markers=len(markers), busco_lineages=sorted(lineages),
                      zero_single_copy=[r['species'] for r in results if 'S' not in r['marker_statuses'].values()],
                      automatic_species_exclusion=False)
        atomic_json(workspace / 'cohort_qc.json', report)
        return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase', choices=['plan', 'prepare', 'verify', 'qc', 'finalize'])
    p.add_argument('--selection', type=Path)
    p.add_argument('--expected-species', type=int)
    p.add_argument('--source-aliases', type=Path, help='JSON map of original raw paths to byte-identical physical replacements')
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--workspace', type=Path)
    p.add_argument('--task-index', type=int)
    a = p.parse_args()
    if a.phase == 'plan':
        aliases = json.loads(a.source_aliases.read_text()) if a.source_aliases else None
        atomic_json(a.plan, freeze(a.selection, a.expected_species, source_aliases=aliases), immutable=True)
        result = dict(planned=True, species=load_cohort(a.plan)['task_count'], plan_sha256=digest(a.plan))
    elif a.phase == 'prepare':
        prepare(a.workspace, a.plan)
        result = dict(prepared=True, plan_sha256=digest(a.plan))
    elif a.phase == 'verify':
        doc = verify_view(a.workspace, a.plan)
        result = dict(verified=True, species=doc['task_count'], plan_sha256=digest(a.plan))
    elif a.phase == 'qc':
        result = qc_one(a.workspace, a.plan, a.task_index)
    else:
        result = finalize(a.workspace, a.plan)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
