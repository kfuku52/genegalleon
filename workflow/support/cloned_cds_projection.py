"""Publish preferred experimental CDS with explicitly reviewed reference projections.

Reviewed targets come from the native genomic candidate audit or independent
coordinate verification. No approximate best hit or reference edit is implicit.
"""
import hashlib
import json
import shutil
from pathlib import Path
from urllib.parse import quote

from experimental_cds_evidence import PREFIX, proof
from input_generation_array_state import FreshDigestBatch, atomic_json, digest

FIELDS = ('cloned_cds_projection_targets', 'cloned_cds_priority_ids')


def implementation():
    files = {name: digest(Path(__file__).with_name(name)) for name in
             ('cloned_cds_projection.py', 'experimental_cds_evidence.py', 'cds_model_normalisation.py')}
    return hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()


def options(row, selected):
    targets = json.loads(str(row.get(FIELDS[0]) or '{}'))
    priority = json.loads(str(row.get(FIELDS[1]) or '[]'))
    if not isinstance(targets, dict) or not targets:
        raise ValueError('Experimental cloned CDS requires reviewed cloned_cds_projection_targets')
    if (not isinstance(priority, list) or any(not isinstance(x, str) or not x for x in priority)
            or len(set(priority)) != len(priority) or set(priority) - targets.keys()):
        raise ValueError('cloned_cds_priority_ids must select unique reviewed clone IDs in priority order')
    if selected and set(targets) != set(selected):
        raise ValueError('Reviewed projection targets must cover exactly the selected cloned CDS')
    allowed = {'seqid', 'strand', 'cds', 'source_gene_ids', 'new_locus', 'partial_5prime', 'partial_3prime'}
    for name, target in targets.items():
        if (not isinstance(name, str) or not name or any(c.isspace() for c in name)
                or not isinstance(target, dict) or set(target) - allowed
                or not {'seqid', 'strand', 'cds', 'source_gene_ids'}.issubset(target)):
            raise ValueError('Invalid reviewed cloned CDS projection')
        if (not isinstance(target['seqid'], str) or not target['seqid'] or any(c.isspace() for c in target['seqid'])
                or target['strand'] not in {'+', '-'}):
            raise ValueError('Invalid reviewed projection contig/strand')
        blocks = target['cds']
        if (not isinstance(blocks, list) or not blocks
                or any(not isinstance(b, list) or len(b) != 3 or any(type(x) is not int for x in b)
                       or not 0 <= b[0] < b[1] or b[2] not in (0, 1, 2) for b in blocks)):
            raise ValueError('Projection CDS blocks require zero-based half-open coordinates and GFF phases')
        ordered = sorted(blocks)
        if any(a[1] > b[0] for a, b in zip(ordered, ordered[1:], strict=False)):
            raise ValueError('Overlapping reviewed projection CDS blocks')
        cumulative = 0
        for a, b, phase in sorted(blocks, reverse=target['strand'] == '-'):
            if phase != (3 - cumulative % 3) % 3:
                raise ValueError('Reviewed projection has inconsistent reference CDS phases')
            cumulative += b - a
        ids = target['source_gene_ids']
        if (not isinstance(ids, list) or any(not isinstance(x, str) or not x or any(c in x for c in '\r\n\t') for x in ids)
                or len(set(ids)) != len(ids)
                or any(type(target.get(k, False)) is not bool for k in ('new_locus', 'partial_5prime', 'partial_3prime'))
                or bool(ids) == target.get('new_locus', False)):
            raise ValueError('Projection must name source genes or explicitly declare one new locus')
    return dict(cloned_cds_projection_targets=targets, cloned_cds_priority_ids=priority)


def overlapping(blocks, other):
    return any(a < d and c < b for a, b, _ in blocks for c, d in other)


def source_owners(target, genes, coding):
    named = set(target['source_gene_ids'])
    if named - genes.keys():
        raise ValueError('Reviewed projection source gene is absent: ' + ','.join(sorted(named - genes.keys())))
    actual = set()
    for gene, value in genes.items():
        if (value['seqid'], value['strand']) != (target['seqid'], target['strand']):
            if gene in named:
                raise ValueError('Reviewed projection source gene has different contig/strand')
            continue
        if overlapping(target['cds'], coding[gene]):
            actual.add(gene)
        if gene in named:
            lo, hi = min(b[0] for b in target['cds']), max(b[1] for b in target['cds'])
            for c, d in coding[gene]:
                missing = [(c, d)]
                for a, b, _ in target['cds']:
                    missing = [(x, y) for u, v in missing for x, y in ((u, min(v, a)), (max(u, b), v)) if x < y]
                for u, v in missing:
                    end = ('5prime' if target['strand'] == '+' else '3prime') if v <= lo else (
                        ('3prime' if target['strand'] == '+' else '5prime') if u >= hi else None)
                    if end is None or not target.get('partial_' + end):
                        raise ValueError('Projection truncates source coding evidence without its matching partial declaration')
    if actual != named:
        raise ValueError('Reviewed projection does not cover exactly its named source gene owners')
    return sorted(actual)


def choose(candidates, clones, priority):
    """Resolve connected same-gene/coding-overlap groups, never duplicate a locus."""
    groups = []
    for name, target in candidates.items():
        connected = [g for g in groups if any(
            set(target['source_gene_ids']) & set(candidates[n]['source_gene_ids'])
            or ((target['seqid'], target['strand']) == (candidates[n]['seqid'], candidates[n]['strand'])
                and overlapping(target['cds'], [(a, b) for a, b, _ in candidates[n]['cds']])) for n in g)]
        group = {name}.union(*connected)
        groups = [g for g in groups if g not in connected] + [group]
    selected, skipped, aliases = [], [], {}
    for group in groups:
        preferred = [name for name in priority if name in group]
        if len(group) > 1 and not preferred and len({clones[name] for name in group}) > 1:
            raise ValueError('Different cloned CDS compete for one gene; provide cloned_cds_priority_ids')
        winner = preferred[0] if preferred else min(group)
        selected.append(winner)
        aliases[winner] = sorted(group - {winner})
        skipped.extend(dict(clone_id=name, status='same_gene_lower_priority', preferred_clone=winner)
                       for name in sorted(group - {winner}))
    return sorted(selected), skipped, aliases


def build_pair(task, directory):
    import format_species_inputs as fsi
    from Bio.Seq import Seq
    from cloned_cds_models import cloned_records, removed_cds_gene_tokens, replacement_gff, source_graph
    from format_species_annotation.genome_intervals import genome_intervals
    clones = cloned_records(Path(task['cloned_cds_path']), task.get('cloned_cds_ids', []), int(task.get('genetic_code', 1)))
    targets = task['cloned_cds_projection_targets']
    if set(targets) != set(clones):
        raise ValueError('Reviewed projections and supplied clone IDs differ')
    records, genes, coding, ownership = source_graph(task['gff_path'])
    for target in targets.values():
        source_owners(target, genes, coding)
    prepared = fsi.prepare_cds_identifier_task(task)
    longer = set(task.get('cloned_cds_longer_only_ids', []))
    lengths = {}
    if longer:
        for header, sequence in fsi.iter_task_cds_records(prepared):
            identifier = fsi.build_formatted_cds_id(prepared, header)
            _, match = fsi.resolve_gene_aggregate_id(prepared, header, identifier)
            if match.get('status') == 'mapped':
                token = match['gene_token']
                lengths[token] = max(lengths.get(token, 0), len(sequence))
    candidates, skipped = {}, []
    for name, target in targets.items():
        if name in longer:
            tokens = removed_cds_gene_tokens(prepared, set(target['source_gene_ids']), genes)
            found = [lengths[t] for t in tokens if t in lengths]
            if not found or len(clones[name]) <= max(found):
                skipped.append(dict(clone_id=name, status='not_longer' if found else 'no_target_CDS',
                                    clone_cds_length=len(clones[name]), target_cds_lengths=found))
                continue
        candidates[name] = target
    selected, duplicates, aliases = choose(candidates, clones, task['cloned_cds_priority_ids'])
    skipped += duplicates
    removed = {gene for name in selected for gene in targets[name]['source_gene_ids']}
    if any(owned & removed and owned - removed for owned in ownership):
        raise ValueError('Replaced GFF feature is shared with a retained source gene')
    identifiers = {attr['ID'] for _, fields, attr, _ in records if fields and attr.get('ID')}
    additions, replacements, reference_records = [], [], []
    regions = [(t['seqid'], a, b) for t in targets.values() for a, b, _ in t['cds']]
    with genome_intervals(task['genome_path'], regions, scratch_dir=directory) as genome:
        for name in selected:
            target = targets[name]
            pieces = [genome.fetch(target['seqid'], a, b).upper() for a, b, _ in sorted(target['cds'], reverse=target['strand'] == '-')]
            reference = ''.join(str(Seq(p).reverse_complement()) if target['strand'] == '-' else p for p in pieces)
            if set(reference) - set('ACGT'):
                raise ValueError('Reviewed experimental projection crosses ambiguous reference bases')
            evidence, differences = proof(reference, clones[name], target['cds'], target['strand'], target['seqid'],
                                           partial_5prime=target.get('partial_5prime'), partial_3prime=target.get('partial_3prime'))
            gene, transcript, text = replacement_gff(task['species_prefix'], name, target, target['source_gene_ids'])
            new_ids = {gene, transcript} | {f'{transcript}.{kind}{i}' for kind in ('exon', 'CDS') for i in range(1, len(target['cds']) + 1)}
            if new_ids & identifiers:
                raise ValueError('Experimental model identifiers collide with source annotation')
            identifiers.update(new_ids)
            attributes = ';gg_sequence_origin=experimental_clone' + ''.join(
                ';' + PREFIX + key + '=' + quote(value, safe='._-:') for key, value in evidence.items())
            attributes += ';gg_clone_aliases=' + ','.join(quote(n, safe='._-:') for n in aliases[name])
            lines = []
            for line in text.splitlines():
                fields = line.split('\t')
                if fields[2] in {'gene', 'mRNA', 'CDS'}:
                    fields[8] += attributes
                lines.append('\t'.join(fields).replace('gg_transcript_extent=CDS_only', 'gg_transcript_extent=reference_projection'))
            additions.append((gene, transcript, clones[name], '\n'.join(lines) + '\n'))
            reference_records.append('>' + transcript + '.reference_projection\n' + reference + '\n')
            replacements.append(dict(clone_id=name, gene_id=gene, transcript_id=transcript,
                                     replaced_gene_ids=target['source_gene_ids'], seqid=target['seqid'], strand=target['strand'],
                                     CDS_intervals_0based_halfopen=target['cds'], new_locus=target.get('new_locus', False),
                                     projection=evidence, differences=differences, aliases=aliases[name],
                                     reference_cds_length=len(reference), clone_cds_length=len(clones[name]),
                                     genome_identical=reference == clones[name]))
    removed_tokens = removed_cds_gene_tokens(prepared, removed, genes)
    import gzip
    def output(path):
        return gzip.open(path, 'wt') if path.suffix == '.gz' else path.open('w')
    cds, gff = directory / Path(task['cds_path']).name, directory / Path(task['gff_path']).name
    removed_cds = []
    removed_path = directory / 'removed_source_cds.fa'
    with output(cds) as stream, removed_path.open('w') as archived:
        for header, sequence in fsi.iter_task_cds_records(prepared):
            identifier = fsi.build_formatted_cds_id(prepared, header)
            _, match = fsi.resolve_gene_aggregate_id(prepared, header, identifier)
            if match.get('status') != 'mapped':
                raise ValueError('Unresolved source CDS owner prevents experimental replacement')
            if match['gene_token'] in removed_tokens:
                removed_cds.append(header.split()[0])
                archived.write('>' + header + '\n' + sequence + '\n')
            else:
                stream.write('>' + header + '\n' + sequence + '\n')
        for gene, transcript, sequence, _ in additions:
            stream.write('>' + transcript + ' gene:' + gene + '\n' + sequence + '\n')
    archived_gff = directory / 'removed_source_models.gff3'
    with output(gff) as stream, archived_gff.open('w') as archived:
        archived.write('##gff-version 3\n')
        for (line, _, _, _), owned in zip(records, ownership, strict=True):
            if not owned & removed:
                stream.write(line)
            else:
                archived.write(line)
        for _, _, _, text in additions:
            stream.write(text)
    genome = directory / 'genome_source' / Path(task['genome_path']).name
    genome.parent.mkdir()
    shutil.copyfile(task['genome_path'], genome)
    (directory / 'reference_projected_cds.fa').write_text(''.join(reference_records))
    (directory / 'received_cloned_cds.fa').write_text(''.join('>' + name + '\n' + sequence + '\n' for name, sequence in clones.items()))
    archive_outputs = {}
    # Empty FASTA has no applicable records and is not a completion artifact.
    # Keep the native nonempty-output guard, rather than inventing a sequence.
    for role, path, size in (('removed_cds', removed_path, len(removed_cds)),
                             ('reference_projected_cds', directory / 'reference_projected_cds.fa', len(reference_records))):
        if size:
            archive_outputs[role] = path.name
        else:
            path.unlink()
    atomic_json(directory / 'audit.json', dict(schema_version=2, mode='experimental', species=task['species_prefix'],
                genomic_CDS_identity_required=False, reference_genome_modified=False, original_sources_modified=False,
                projections_reviewed=True, replacements=replacements, skipped=skipped, removed_CDS_ids=removed_cds,
                archive_record_counts=dict(removed_cds=len(removed_cds), reference_projected_cds=len(reference_records),
                                           received_cloned_cds=len(clones)),
                removed_CDS_gene_tokens=sorted(removed_tokens), outputs=dict(cds=cds.name, gff=gff.name,
                genome=str(genome.relative_to(directory)), removed_gff=archived_gff.name,
                received_cloned_cds='received_cloned_cds.fa', **archive_outputs)))


def apply(task, root):
    from Bio.Align import _pairwisealigner
    from rescue_gene_models import stage
    paths = [Path(task[k]).resolve() for k in ('cds_path', 'gff_path', 'genome_path', 'cloned_cds_path')]
    batch = FreshDigestBatch()
    hashes = batch.read(paths)
    if hashes[str(paths[-1])] != task['cloned_cds_sha256']:
        raise ValueError('Experimental cloned CDS changed after qualification')
    key = dict(schema=1, mode='experimental', inputs=hashes, species=task['species_prefix'],
               targets=task['cloned_cds_projection_targets'], priority=task['cloned_cds_priority_ids'],
               longer_only=task.get('cloned_cds_longer_only_ids', []), genetic_code=task.get('genetic_code', 1),
               implementation=task['cloned_model_implementation'], aligner=digest(_pairwisealigner.__file__),
               dependencies={name: digest(Path(__file__).with_name(name)) for name in
                             ('cloned_cds_projection.py', 'experimental_cds_evidence.py', 'gene_model_refinement.py',
                              'cds_model_normalisation.py', 'format_species_inputs.py')})
    token = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()
    destination = stage(Path(root), task['species_prefix'] + '/' + token, key,
                        lambda folder: build_pair(task, folder), guard=batch.check, publication_guard=batch.check)
    audit = json.loads((destination / 'audit.json').read_text())
    for role in ('cds', 'gff', 'genome'):
        task[role + '_path'] = destination / audit['outputs'][role]
    task['cloned_model_audit'] = str(destination / 'audit.json')
    task['_cloned_models_applied'] = True
