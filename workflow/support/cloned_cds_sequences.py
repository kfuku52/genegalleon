"""Audited clone-priority sequence replacement for inputs without genomic annotation."""

import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
from contextlib import nullcontext
from pathlib import Path

from input_generation_array_state import FreshDigestBatch, atomic_json, digest

GROUP_FIELD = 'cloned_cds_sequence_groups'


def reviewed_groups(row, selected):
    """Bind reviewed transcript ownership without asserting genomic loci."""
    raw = row.get(GROUP_FIELD)
    if not raw:
        return {}
    groups = json.loads(str(raw))
    if not selected or not isinstance(groups, dict) or not groups:
        raise ValueError('Reviewed sequence groups require explicit selected clone IDs')
    clones, sources = set(), set()
    for name, group in groups.items():
        if (not isinstance(name, str) or re.fullmatch(r'[A-Za-z][A-Za-z0-9_.-]*', name) is None
                or not isinstance(group, dict) or not {'clone_ids', 'source_cds_ids', 'unplaced'} <= set(group)
                or set(group) - {'clone_ids', 'source_cds_ids', 'unplaced', 'analysis_policy'}
                or type(group['unplaced']) is not bool
                or not isinstance(group.get('analysis_policy', 'representative_only'), str)
                or group.get('analysis_policy', 'representative_only') not in {'representative_only', 'all_clones'}):
            raise ValueError('Invalid reviewed sequence group')
        for field, seen in (('clone_ids', clones), ('source_cds_ids', sources)):
            values = group[field]
            if (not isinstance(values, list) or field == 'clone_ids' and not values
                    or any(not isinstance(v, str) or not v or any(c.isspace() for c in v) for v in values)
                    or len(set(values)) != len(values) or seen.intersection(values)):
                raise ValueError('Reviewed sequence groups require unique clone and source IDs')
            seen.update(values)
        if bool(group['source_cds_ids']) == group['unplaced']:
            raise ValueError('Reviewed sequence group must name source CDS or explicitly declare unplaced')
    if clones != set(selected):
        raise ValueError('Reviewed sequence groups must cover exactly the selected cloned CDS')
    return groups


def apply_reviewed_groups(clones, decisions, groups, seen):
    """Apply reviewed analysis inclusion without deciding allele versus paralog."""
    rows = {r['clone_id']: r for r in decisions}
    if groups and {c for g in groups.values() for c in g['clone_ids']} != set(clones):
        raise ValueError('Reviewed sequence group selection changed')
    evidence = []
    for name, group in groups.items():
        ids, source_ids = group['clone_ids'], set(group['source_cds_ids'])
        if source_ids - seen:
            raise ValueError('Reviewed source CDS is absent: ' + ','.join(sorted(source_ids - seen)))
        for clone in ids:
            row = rows[clone]
            targets = {m['target_id'] for m in row['matches']}
            if (group['unplaced'] and row['status'] != 'no_supported_match'
                    or not group['unplaced'] and (row['status'] not in {'replace', 'competing_clones'}
                                                  or targets != source_ids)):
                raise ValueError('Reviewed sequence group conflicts with current search evidence: ' + clone)
        representative = ids[0]
        policy = group.get('analysis_policy', 'representative_only')
        analysis_clones = ids if policy == 'all_clones' else ids[:1]
        for clone in ids:
            row = rows[clone]
            row.update(search_status=row['status'], candidate_group=name, representative_clone=representative,
                       clone_sha256=hashlib.sha256(clones[clone].encode()).hexdigest(),
                       genomic_origin_status='not_assessed_no_genome', genomic_coordinates=None,
                       analysis_policy=policy, analysis_included=clone in analysis_clones,
                       copy_number_status='unresolved', confirmed_gene_copy_count=None,
                       source_catalog_status='no_supported_match_in_sequence_catalog' if group['unplaced'] else 'matched_sequence',
                       gene_identity_status='unresolved_candidate_group' if len(ids) > 1 else 'unresolved_unplaced' if group['unplaced'] else 'sequence_only')
            row['status'] = ('add_unplaced' if group['unplaced'] else 'replace') if clone == representative else 'alternative_clone'
            if clone != representative and policy == 'all_clones':
                row['status'] = 'retain_unresolved_clone'
        evidence.append(dict(candidate_group=name, representative_clone=representative,
                             alternative_clones=ids[1:], source_cds_ids=sorted(source_ids), unplaced=group['unplaced'],
                             analysis_policy=policy, analysis_clones=analysis_clones,
                             copy_number_status='unresolved', confirmed_gene_copy_count=None,
                             representative_selection='explicit_clone_order; no functional or genomic priority inferred'))
    return evidence


def _run(argv, directory, name):
    with (directory / (name + '.stdout')).open('w') as stdout, (directory / (name + '.stderr')).open('w') as stderr:
        subprocess.run([str(v) for v in argv], stdout=stdout, stderr=stderr, check=True)


def search_candidates(path, query_names, targets, identity, coverage, query_lengths):
    """Find source fragments or exact complete clones within extended source CDS."""
    grouped = {name: {} for name in query_names.values()}
    subjects = {name: set() for name in query_names.values()}
    for line in path.read_text().splitlines():
        if line.startswith('#'):
            continue
        fields = line.split('\t')
        if len(fields) != 14 or fields[0] not in query_names or fields[1] not in targets:
            raise ValueError('Unexpected cloned sequence search row')
        query = query_names[fields[0]]
        subject = fields[1]
        subjects[query].add(subject)
        qstart, qend, sstart, send, qlen, slen = map(int, fields[4:10])
        if not (1 <= qstart <= qend <= qlen and 1 <= min(sstart, send) <= max(sstart, send) <= slen
                and slen == targets[subject]['length'] and qlen == query_lengths[query]):
            raise ValueError('Cloned sequence alignment coordinates changed')
        target_coverage = (abs(send - sstart) + 1) / slen
        complete_clone_exact = (qstart == 1 and qend == qlen and len(fields[12]) == qlen
                                and fields[12].upper() == fields[13].upper())
        if float(fields[2]) / 100 < identity or (target_coverage < coverage and not complete_clone_exact):
            continue
        if (len(fields[12]) != int(fields[3]) or len(fields[13]) != int(fields[3])
                or len(fields[13].replace('-', '')) != abs(send - sstart) + 1):
            raise ValueError('Cloned sequence alignment length changed')
        bases = {}
        insertions = {}
        position = qstart
        pending = ''
        for qbase, sbase in zip(fields[12], fields[13], strict=True):
            if qbase == '-':
                pending += sbase.upper()
            else:
                bases[position] = sbase.upper()
                if position > qstart or pending:
                    insertions[position] = pending
                pending = ''
                position += 1
        if pending:
            insertions[position] = pending
        if position != qend + 1:
            raise ValueError('Cloned sequence alignment length changed')
        record = dict(subject=subject, target_id=targets[subject]['id'], target_length=slen,
                      identity=float(fields[2]) / 100, target_coverage=target_coverage,
                      query_coverage=(qend - qstart + 1) / qlen,
                      match_basis='source_CDS_coverage' if target_coverage >= coverage else 'complete_clone_exact',
                      query_start=qstart, query_end=qend, strand='+' if send >= sstart else '-',
                      bases=bases, insertions=insertions)
        grouped[query].setdefault(subject, []).append(record)
    return grouped, {name for name, hits in subjects.items() if len(hits) >= 1000}


def choose_replacements(clones, grouped, longer_only):
    decisions = {}
    for name, sequence in clones.items():
        hits = grouped[name]
        row = dict(clone_id=name, clone_cds_length=len(sequence), matches=[])
        decisions[name] = row
        if not hits:
            row['status'] = 'no_supported_match'
            continue
        if any(len(records) != 1 for records in hits.values()):
            row['status'] = 'ambiguous_alignment'
            continue
        records = [records[0] for records in hits.values()]
        row['matches'] = [{k: v for k, v in r.items() if k not in {'bases', 'insertions'}} for r in records]
        if name in longer_only and len(sequence) <= max(r['target_length'] for r in records):
            row['status'] = 'not_longer'
            continue
        if len(records) > 1 and any(r['match_basis'] == 'complete_clone_exact' for r in records):
            row['status'] = 'ambiguous_extended_targets'
            continue
        conflict = False
        for field in ('bases', 'insertions'):
            supported = {}
            for record in records:
                for position, base in record[field].items():
                    if position in supported and supported[position] != base:
                        conflict = True
                    supported[position] = base
        row['status'] = 'ambiguous_overlapping_targets' if conflict else 'replace'
    ownership = {}
    for name, row in decisions.items():
        if row['status'] == 'replace':
            for match in row['matches']:
                ownership.setdefault(match['subject'], []).append(name)
    for names in ownership.values():
        if len(names) > 1:
            # An incoming reference/clone mixture cannot silently choose between
            # alternative sequences competing for the same existing transcript.
            for name in names:
                decisions[name]['status'] = 'competing_clones'
    return list(decisions.values())


def build_sequence(task, directory):
    from cloned_cds_models import cloned_records
    from format_species_annotation.common import iter_fasta_records
    from format_species_annotation.tasks import build_gene_aggregate_id, iter_task_cds_records
    from transcriptome_orfs import write_audit

    clones = cloned_records(Path(task['cloned_cds_path']), task.get('cloned_cds_ids', []), int(task.get('genetic_code', 1)))
    longer = set(task.get('cloned_cds_longer_only_ids', []))
    if longer - clones.keys():
        raise ValueError('Longer-only policy names an unselected cloned CDS')
    source_task = dict(task)
    source_task.pop('_transcriptome_orf_audit', None)
    targets = {}
    seen = set()
    source = directory / 'search_subjects.fa'
    with source.open('w') as stream:
        for i, (header, sequence) in enumerate(iter_task_cds_records(source_task), 1):
            identifier = header.split()[0]
            if identifier in seen:
                raise ValueError('Duplicate source CDS identifier: ' + identifier)
            seen.add(identifier)
            key = 'target' + str(i)
            targets[key] = dict(id=identifier, length=len(sequence))
            stream.write('>' + key + '\n' + sequence + '\n')
    if not targets:
        raise ValueError('No source coding sequences available for clone search')
    write_audit(source_task, source, atomic_json)
    queries = directory / 'search_queries.fa'
    names = {'clone' + str(i): name for i, name in enumerate(clones, 1)}
    queries.write_text(''.join('>' + key + '\n' + clones[name] + '\n' for key, name in names.items()))
    _run(['makeblastdb', '-in', source, '-dbtype', 'nucl', '-parse_seqids', '-out', directory / 'source_db'], directory, 'makeblastdb')
    alignments = directory / 'alignments.tsv'
    _run(['blastn', '-task', 'blastn', '-query', queries, '-db', directory / 'source_db', '-strand', 'plus',
          '-dust', 'no', '-evalue', '1e-20', '-max_target_seqs', '1000',
          '-num_threads', str(max(1, int(os.environ.get('GG_TASK_CPUS', '1')))),
          '-outfmt', '7 qseqid sseqid pident length qstart qend sstart send qlen slen bitscore evalue qseq sseq',
          '-out', alignments], directory, 'blastn')
    grouped, capped = search_candidates(alignments, names, targets, task['cloned_cds_min_identity'],
                                        task['cloned_cds_min_target_coverage'],
                                        {name: len(sequence) for name, sequence in clones.items()})
    decisions = choose_replacements(clones, grouped, longer)
    # Reaching the search-result cap cannot establish unique target membership.
    for row in decisions:
        if row['clone_id'] in capped:
            row['status'] = 'search_limit_reached'
    groups = task.get(GROUP_FIELD, {})
    if any(g['unplaced'] and longer.intersection(g['clone_ids']) for g in groups.values()):
        raise ValueError('Unplaced sequence cannot establish a longer-only comparison')
    group_evidence = apply_reviewed_groups(clones, decisions, groups, seen)
    accepted = [r for r in decisions if r['status'] in {'replace', 'add_unplaced', 'retain_unresolved_clone'}]
    removed = {m['target_id'] for r in accepted for m in r['matches']}
    additions = {'gg_clone_' + hashlib.sha256(r['clone_id'].encode()).hexdigest()[:16] + '.t1': r for r in accepted}
    if additions.keys() & seen:
        raise ValueError('Cloned sequence ID collides with source input')
    output = directory / Path(task['cds_path']).name
    removed_path = directory / 'removed_sequences.fa'
    opener = gzip.open if output.suffix == '.gz' else open
    with opener(output, 'wt') as retained, (removed_path.open('w') if removed else nullcontext()) as discarded:
        for header, sequence in iter_fasta_records(task['cds_path']):
            stream = discarded if header.split()[0] in removed else retained
            stream.write('>' + header + '\n' + sequence + '\n')
        for identifier, row in additions.items():
            retained.write('>' + identifier + '\n' + clones[row['clone_id']] + '\n')
            row['output_id'] = identifier
            row['formatted_output_id'] = build_gene_aggregate_id(task, identifier, identifier)
    outputs = dict(cds=output.name, alignments=alignments.name)
    if groups:
        received = directory / 'received_cloned_cds.fa'
        received.write_text(''.join('>' + name + '\n' + sequence + '\n' for name, sequence in clones.items()))
        outputs['received_cloned_cds'] = received.name
        atomic_json(directory / 'reviewed_sequence_groups.json', groups)
        outputs['reviewed_sequence_groups'] = 'reviewed_sequence_groups.json'
        members = directory / 'sequence_group_members.tsv'
        fields = ['clone_id', 'candidate_group', 'output_id', 'formatted_output_id', 'analysis_policy',
                  'analysis_included', 'gene_identity_status', 'copy_number_status', 'confirmed_gene_copy_count']
        with members.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            writer.writerows(r for r in decisions if 'candidate_group' in r)
        outputs['sequence_group_members'] = members.name
        alternatives = [r for r in decisions if r['status'] in {'alternative_clone', 'retain_unresolved_clone'}]
        if alternatives:
            path = directory / 'alternative_cloned_cds.fa'
            path.write_text(''.join('>' + r['clone_id'] + ' candidate_group=' + r['candidate_group']
                                   + ' representative_clone=' + r['representative_clone'] + '\n'
                                   + clones[r['clone_id']] + '\n' for r in alternatives))
            outputs['alternative_cloned_cds'] = path.name
    if removed:
        outputs['removed'] = removed_path.name
    if Path(str(source) + '.orf-extraction.json').exists():
        outputs['source_orf_audit'] = source.name + '.orf-extraction.json'
    atomic_json(directory / 'audit.json', dict(schema_version=1, mode='sequence', species=task['species_prefix'],
        coordinate_system='sequence only; no genomic coordinates inferred',
        minimum_identity=task['cloned_cds_min_identity'], minimum_target_coverage=task['cloned_cds_min_target_coverage'],
        decisions=decisions, removed_CDS_ids=sorted(removed), original_sources_modified=False, outputs=outputs,
        **(dict(reviewed_sequence_groups=group_evidence) if groups else {})))


def apply_sequence_task(task, root):
    from rescue_gene_models import stage
    from transcriptome_orfs import input_type

    if not task.get('cds_path') or any(task.get(role + '_path') for role in ('gff', 'gbff', 'genome')):
        raise ValueError('Sequence replacement requires a sequence-only CDS/transcriptome source')
    if task.get('cloned_sequence_implementation') != digest(Path(__file__)):
        raise ValueError('Cloned sequence implementation changed; prepare a new native task plan')
    tools = {name: shutil.which(name) for name in ('blastn', 'makeblastdb')}
    if not all(tools.values()):
        raise ValueError('blastn and makeblastdb are required for cloned sequence replacement')
    sources = [Path(task[k]).resolve() for k in ('cds_path', 'cloned_cds_path')]
    batch = FreshDigestBatch()
    hashes = batch.read(sources)
    if hashes[str(sources[-1])] != task['cloned_cds_sha256']:
        raise ValueError('Cloned CDS changed after input qualification')
    key = dict(schema=1, inputs=hashes, species=task['species_prefix'],
               selection=task.get('cloned_cds_ids', []), longer_only=task.get('cloned_cds_longer_only_ids', []),
               reviewed_sequence_groups=task.get(GROUP_FIELD, {}),
               identity=task['cloned_cds_min_identity'], coverage=task['cloned_cds_min_target_coverage'],
               input_type=input_type(task), genetic_code=int(task.get('genetic_code', 1)),
               implementation=digest(Path(__file__)), tools={name: digest(path) for name, path in tools.items()},
               dependencies={name: digest(Path(__file__).with_name(name)) for name in
                             ('cloned_cds_models.py', 'transcriptome_orfs.py')})
    if input_type(task) == 'transcriptome':
        from cdskit import longestcds
        key['orf_implementation'] = digest(Path(longestcds.__file__))
    token = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()
    destination = stage(Path(root), task['species_prefix'] + '/' + token, key,
                        lambda folder: build_sequence(task, folder), guard=batch.check, publication_guard=batch.check)
    audit = json.loads((destination / 'audit.json').read_text())
    task['cds_path'] = destination / audit['outputs']['cds']
    task['cloned_model_audit'] = str(destination / 'audit.json')
    task['_cloned_models_applied'] = True
