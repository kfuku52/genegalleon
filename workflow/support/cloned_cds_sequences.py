"""Audited clone-priority sequence replacement for inputs without genomic annotation."""

import gzip
import hashlib
import json
import os
import shutil
import subprocess
from contextlib import nullcontext
from pathlib import Path

from input_generation_array_state import FreshDigestBatch, atomic_json, digest


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
    from format_species_annotation.tasks import iter_task_cds_records
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
    accepted = [r for r in decisions if r['status'] == 'replace']
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
    outputs = dict(cds=output.name, alignments=alignments.name)
    if removed:
        outputs['removed'] = removed_path.name
    if Path(str(source) + '.orf-extraction.json').exists():
        outputs['source_orf_audit'] = source.name + '.orf-extraction.json'
    atomic_json(directory / 'audit.json', dict(schema_version=1, mode='sequence', species=task['species_prefix'],
        coordinate_system='sequence only; no genomic coordinates inferred',
        minimum_identity=task['cloned_cds_min_identity'], minimum_target_coverage=task['cloned_cds_min_target_coverage'],
        decisions=decisions, removed_CDS_ids=sorted(removed), original_sources_modified=False, outputs=outputs))


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
