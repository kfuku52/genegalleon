#!/usr/bin/env python3
"""Curate source models or sequence-only inputs with explicitly supplied cloned CDS.

This is a coupled source curation step before ordinary native formatting. The
explicit experimental mode records a reviewed reference projection and its
nucleotide differences rather than asserting genome identity. Genomes are unchanged.
"""
import argparse
import bisect
import csv
import functools
import hashlib
import json
import os
import re
import shutil
import sys
from pathlib import Path
from urllib.parse import quote, unquote, urlparse

from input_generation_array_state import FreshDigestBatch, atomic_json, digest

FIELDS = ('cloned_cds_fasta', 'cloned_cds_sha256', 'cloned_cds_ids', 'cloned_cds_mode',
          'cloned_cds_longer_only_ids', 'cloned_cds_min_identity', 'cloned_cds_min_target_coverage',
          'cloned_cds_unresolved_policy', 'cloned_cds_index_batch_bases',
          'cloned_cds_projection_targets', 'cloned_cds_priority_ids', 'cloned_cds_sequence_groups')
TASK_FIELDS = ('cloned_cds_path', 'cloned_cds_sha256', 'cloned_cds_ids', 'cloned_model_implementation',
               'cloned_cds_mode', 'cloned_cds_longer_only_ids', 'cloned_cds_min_identity',
               'cloned_cds_min_target_coverage', 'cloned_sequence_implementation', 'cloned_cds_unresolved_policy',
               'cloned_cds_index_batch_bases', 'cloned_genomic_search_implementation',
               'cloned_cds_projection_targets', 'cloned_cds_priority_ids', 'cloned_projection_implementation',
               'cloned_cds_sequence_groups')
EVIDENCE_FIELDS = ('cloned_cds_path', 'cloned_model_audit_path', 'cloned_source_cds_path',
                   'cloned_source_gff_path', 'cloned_source_genome_path')


def evidence_paths(record, *, summary=False):
    mode = record.get('cloned_cds_mode') or 'genomic'
    if mode not in {'genomic', 'sequence', 'experimental'}:
        raise ValueError('Invalid cloned model evidence mode')
    fields = ('cloned_model_audit_path', 'cds_input_path') if summary else EVIDENCE_FIELDS[:3]
    if mode in {'genomic', 'experimental'}:
        fields += ('gff_input_path', 'genome_input_path') if summary else EVIDENCE_FIELDS[3:]
    paths = []
    for field in fields:
        value = record.get(field)
        if not value or any(c in value for c in '\r\n') or not Path(value).is_file() or not Path(value).stat().st_size:
            raise ValueError('Missing or unsafe cloned model evidence: ' + field)
        paths.append(value)
    audit_path = Path(record['cloned_model_audit_path'])
    audit = json.loads(audit_path.read_text())
    if (audit.get('mode') or 'genomic') != mode:
        raise ValueError('Cloned model evidence mode changed')
    if mode in {'sequence', 'experimental'}:
        for value in audit['outputs'].values():
            relative = Path(value)
            if relative.is_absolute() or '..' in relative.parts or any(c in value for c in '\r\n'):
                raise ValueError('Unsafe cloned sequence artifact path')
            path = audit_path.parent / relative
            if not path.is_file() or path.is_symlink():
                raise ValueError('Missing cloned sequence artifact: ' + str(relative))
            paths.append(str(path))
    return list(dict.fromkeys(paths))


def attach(task, row, parent):
    """Bind only a local, hash-qualified FASTA; a mixed FASTA needs explicit IDs."""
    if not any(str(row.get(k) or '').strip() for k in FIELDS):
        return task
    reference = str(row.get('cloned_cds_fasta') or '').strip()
    expected = str(row.get('cloned_cds_sha256') or '').strip()
    if not reference or not re.fullmatch(r'[0-9a-f]{64}', expected):
        raise ValueError('Cloned CDS requires cloned_cds_fasta and its exact cloned_cds_sha256')
    mode = str(row.get('cloned_cds_mode') or 'genomic').strip()
    if mode not in {'genomic', 'sequence', 'experimental'}:
        raise ValueError('cloned_cds_mode must be genomic, sequence or experimental')
    if mode != 'sequence' and any(row.get(k) for k in ('cloned_cds_min_identity', 'cloned_cds_min_target_coverage')):
        raise ValueError('Sequence search thresholds require cloned_cds_mode=sequence')
    if row.get('paired_curation') or (mode != 'sequence' and str(row.get('cds_input_type') or 'cds').strip().lower() != 'cds'):
        raise ValueError('Cloned model replacement requires genomic CDS/GFF inputs, without paired_curation')
    if mode == 'sequence' and any(row.get(prefix + role + suffix) for role in ('gff', 'gbff', 'genome')
                                  for prefix, suffix in (('', '_url'), ('', '_path'), ('local_', '_path'))):
        raise ValueError('Sequence replacement cannot override a paired genomic annotation')
    parsed = urlparse(reference)
    if parsed.scheme and (parsed.scheme != 'file' or parsed.netloc or parsed.query or parsed.fragment):
        raise ValueError('cloned_cds_fasta must be a local path or file URL')
    path = Path(unquote(parsed.path) if parsed.scheme else reference).expanduser()
    path = (Path(parent) / path).resolve() if not path.is_absolute() else path.resolve()
    if digest(path) != expected:
        raise ValueError('Cloned CDS SHA256 mismatch: ' + str(path))
    selected = str(row.get('cloned_cds_ids') or '').strip()
    ids = json.loads(selected) if selected else []
    if (not isinstance(ids, list) or any(not isinstance(v, str) or not v or re.search(r'\s', v) for v in ids)
            or len(set(ids)) != len(ids)):
        raise ValueError('cloned_cds_ids must be a JSON array of unique FASTA IDs')
    longer = json.loads(str(row.get('cloned_cds_longer_only_ids') or '[]'))
    if (not isinstance(longer, list) or any(not isinstance(v, str) or not v or re.search(r'\s', v) for v in longer)
            or len(set(longer)) != len(longer) or ids and set(longer) - set(ids)):
        raise ValueError('cloned_cds_longer_only_ids must select unique cloned CDS IDs')
    task.update(cloned_cds_path=path, cloned_cds_sha256=expected, cloned_cds_ids=ids,
                cloned_model_implementation=digest(Path(__file__)), cloned_cds_mode=mode,
                cloned_cds_longer_only_ids=longer)
    from cloned_cds_projection import FIELDS as PROJECTION_FIELDS
    if mode == 'experimental':
        from cloned_cds_projection import implementation, options
        task.update(options(row, ids), cloned_projection_implementation=implementation())
    elif any(row.get(k) for k in PROJECTION_FIELDS):
        raise ValueError('Reviewed projection options require cloned_cds_mode=experimental')
    if mode == 'genomic':
        policy = str(row.get('cloned_cds_unresolved_policy') or 'error')
        if policy not in {'error', 'retain_source'}:
            raise ValueError('cloned_cds_unresolved_policy must be error or retain_source')
        from cloned_cds_search import implementation
        raw = str(row.get('cloned_cds_index_batch_bases') or 'auto')
        if raw != 'auto' and (not raw.isdecimal() or int(raw) < 1):
            raise ValueError('cloned_cds_index_batch_bases must be auto or a positive integer')
        task.update(cloned_cds_unresolved_policy=policy, cloned_cds_index_batch_bases='auto' if raw == 'auto' else int(raw),
                    cloned_genomic_search_implementation=implementation())
    elif any(row.get(k) for k in ('cloned_cds_unresolved_policy', 'cloned_cds_index_batch_bases')):
        raise ValueError('Genomic cloned CDS options require genomic mode')
    if mode == 'sequence':
        from cloned_cds_sequences import GROUP_FIELD, reviewed_groups
        groups = reviewed_groups(row, ids)
        if groups:
            task[GROUP_FIELD] = groups
        else:
            task.pop(GROUP_FIELD, None)
        task['cloned_sequence_implementation'] = digest(Path(__file__).with_name('cloned_cds_sequences.py'))
        for field in ('cloned_cds_min_identity', 'cloned_cds_min_target_coverage'):
            raw = row.get(field)
            value = float(.98 if raw is None or raw == '' else raw)
            if not 0 < value <= 1:
                raise ValueError(field + ' must be greater than zero and at most one')
            task[field] = value
    elif row.get('cloned_cds_sequence_groups'):
        raise ValueError('Reviewed sequence groups require cloned_cds_mode=sequence')
    return task


def copy_binding(task, source):
    for key in TASK_FIELDS:
        if key in source:
            task[key] = source[key]
    return task


def binding_signature(task):
    if not task.get('cloned_cds_path'):
        return ''
    binding = {k: str(task[k]) if isinstance(task[k], Path) else task[k] for k in TASK_FIELDS if k in task}
    return hashlib.sha256(json.dumps(binding, sort_keys=True).encode()).hexdigest()


def verify_manifest(path):
    from format_species_manifest import read_download_manifest
    entries = []
    for row in read_download_manifest(path):
        task = attach({}, row, path.parent)
        if task:
            entries.append({k: str(v) if isinstance(v, Path) else v for k, v in task.items()})
    return hashlib.sha256(json.dumps(entries, sort_keys=True).encode()).hexdigest() if entries else ''


def cloned_records(path, selected, code):
    from Bio.Data import CodonTable
    from Bio.Seq import Seq
    from format_species_annotation.common import iter_fasta_records
    from format_species_annotation.rna import extract_input_cds
    table = CodonTable.unambiguous_dna_by_id[code]
    records = {}
    seen = set()
    for header, sequence in iter_fasta_records(path):
        name = header.split()[0]
        if name in seen:
            raise ValueError('Duplicate cloned FASTA ID: ' + name)
        seen.add(name)
        if selected and name not in selected:
            continue
        sequence = extract_input_cds({'genetic_code': code}, header, sequence)
        if (len(sequence) < 9 or len(sequence) % 3 or set(sequence) - set('ACGT')
                or sequence[:3] not in table.start_codons or sequence[-3:] not in table.stop_codons
                or '*' in str(Seq(sequence[:-3]).translate(table=code))):
            raise ValueError('Cloned CDS must be a complete, unambiguous, intact CDS including stop: ' + name)
        records[name] = sequence
    if not records or set(selected) - records.keys():
        raise ValueError('Missing selected cloned CDS IDs: ' + ', '.join(sorted(set(selected) - records.keys())))
    return records


def source_graph(gff):
    from fasta_sequence_store import open_text
    from gene_model_refinement import _gff_source_graph, _selected_gff_source_rows
    with open_text(gff) as handle:
        graph = _gff_source_graph(handle)
    records, parents, implicit = graph
    known = {attrs['ID'] for _, fields, attrs, _ in records if fields and attrs.get('ID')}
    missing_parents = any(tid not in known or info[0] not in known for tid, info in implicit.items())
    missing_parents |= any(parent not in known for _, fields, _, pids in records
                           if fields and fields[2].lower().endswith(('rna', 'transcript')) for parent in pids)
    if missing_parents:
        # Parse only the synthetic parents. Original records were normalized
        # during their first read; never serialize/reparse the full annotation.
        prefix, prefix_parents, _ = _gff_source_graph(_selected_gff_source_rows(
            graph, set(), retain_all=True, synthetic_parents_only=True))
        for identifier, parent_ids in prefix_parents.items():
            parents.setdefault(identifier, set()).update(parent_ids)
        del prefix_parents
    else:
        prefix = [('##gff-version 3\n', None, {}, [])]
    for _, fields, attrs, pids in records:
        if fields and pids:
            attrs.setdefault('Parent', pids[0])
    records = prefix + [row for row in records if not row[0].startswith('##gff-version')]
    del graph, implicit, known, prefix
    genes = {}
    for _line, fields, attrs, _parents in records:
        if fields and fields[2].lower().endswith('gene'):
            name = attrs.get('ID')
            if not name or name in genes or fields[6] not in ('+', '-'):
                raise ValueError('Missing, duplicate or unstranded source gene ID')
            genes[name] = dict(seqid=fields[0], strand=fields[6], start=int(fields[3]) - 1, end=int(fields[4]))

    active = set()

    @functools.lru_cache(None)
    def owners(identifier):
        if identifier in active:
            raise ValueError('Cyclic GFF parent graph')
        if identifier in genes:
            return frozenset([identifier])
        if identifier not in parents:
            raise ValueError('Unbound GFF parent: ' + identifier)
        active.add(identifier)
        try:
            inherited = [owners(parent) for parent in parents[identifier]]
            return inherited[0] if len(inherited) == 1 else frozenset().union(*inherited)
        finally:
            active.remove(identifier)

    coding = {g: [] for g in genes}
    ownership = []
    for _line, fields, attrs, parent_ids in records:
        owned = frozenset()
        if fields:
            if attrs.get('ID'):
                owned = owners(attrs['ID'])
            elif len(parent_ids) == 1:
                owned = owners(parent_ids[0])
            else:
                owned = frozenset(g for p in parent_ids for g in owners(p))
            if fields[2] == 'CDS':
                if not owned:
                    raise ValueError('CDS has no source gene owner')
                for gene in owned:
                    if (fields[0], fields[6]) != (genes[gene]['seqid'], genes[gene]['strand']):
                        raise ValueError('Source gene/CDS coordinate ownership conflict')
                    coding[gene].append((int(fields[3]) - 1, int(fields[4])))
        ownership.append(owned)
    owners.cache_clear()
    parents.clear()
    return records, genes, coding, ownership


def overlapping_genes(model, genes, coding):
    intervals = [(s, e) for s, e, _ in model['cds']]
    matches = []
    for name, gene in genes.items():
        if (gene['seqid'], gene['strand']) != (model['seqid'], model['strand']):
            continue
        blocks = coding[name]
        if not blocks:
            if any(gene['start'] < d and c < gene['end'] for c, d in intervals):
                raise ValueError('Overlapping gene lacks coding evidence for replacement: ' + name)
            continue
        if not any(a < d and c < b for a, b in blocks for c, d in intervals):
            continue
        # A split model is explained by the same coding path. An intersecting
        # independent model with coding bases outside that path is not deleted.
        if any(not any(c <= a < b <= d for c, d in intervals) for a, b in blocks):
            raise ValueError('Overlapping model has unexplained coding coordinates: ' + name)
        matches.append(name)
    if not matches:
        raise ValueError('No existing gene model overlaps the cloned coding path')
    return sorted(matches)


def overlap_gene_candidates(genes, coding, models):
    """Query coding extents once per clone, retaining original admission order."""
    if len(models) < 2:
        return None
    groups = {(model['seqid'], model['strand']): [] for model in models.values()}
    for rank, (name, gene) in enumerate(genes.items()):
        group = groups.get((gene['seqid'], gene['strand']))
        if group is not None:
            blocks = coding[name]
            start = min(a for a, _ in blocks) if blocks else gene['start']
            end = max(b for _, b in blocks) if blocks else gene['end']
            group.append((start, end, rank, name))
    indexes = {}
    for key, group in groups.items():
        group.sort()
        maximum = -1
        ends = []
        for _, end, _, _ in group:
            maximum = max(maximum, end)
            ends.append(maximum)
        indexes[key] = (group, [row[0] for row in group], ends)
    result = {}
    for name, model in models.items():
        matches = []
        if model['cds']:
            start = min(a for a, _, _ in model['cds'])
            end = max(b for _, b, _ in model['cds'])
            rows, starts, ends = indexes[(model['seqid'], model['strand'])]
            i = bisect.bisect_left(starts, end) - 1
            while i >= 0 and ends[i] > start:
                if rows[i][1] > start:
                    matches.append(rows[i])
                i -= 1
        # The original first conflicting model must still be the first error.
        matches.sort(key=lambda row: row[2])
        result[name] = {row[3]: genes[row[3]] for row in matches}
    return result


def replacement_gff(name, clone_id, model, replaced):
    gene = 'gg_clone_' + hashlib.sha256(clone_id.encode()).hexdigest()[:16]
    transcript = gene + '.t1'
    seqid, strand = model['seqid'], model['strand']
    start = min(s for s, _, _ in model['cds']) + 1
    end = max(e for _, e, _ in model['cds'])
    evidence = ';gg_cloned_cds=' + quote(clone_id, safe='._-:') + ';gg_replaces=' + ','.join(quote(g, safe='._-:') for g in replaced)
    evidence += ';gg_transcript_extent=CDS_only;gg_species=' + quote(name, safe='._-:')
    def row(kind, s, e, phase, attributes):
        return '\t'.join(map(str, (seqid, 'GeneGalleon_cloned_CDS', kind, s, e, '.', strand, phase, attributes))) + '\n'
    lines = [row('gene', start, end, '.', 'ID=' + gene + evidence),
             row('mRNA', start, end, '.', 'ID=' + transcript + ';Parent=' + gene + evidence)]
    for i, (s, e, phase) in enumerate(model['cds'], 1):
        lines += [row('exon', s + 1, e, '.', f'ID={transcript}.exon{i};Parent={transcript}'),
                  row('CDS', s + 1, e, phase, f'ID={transcript}.CDS{i};Parent={transcript}')]
    return gene, transcript, ''.join(lines)


def removed_cds_gene_tokens(prepared, removed, genes):
    """Resolve source feature IDs through the same native CDS ownership index."""
    aliases = prepared['_gff_cds_grouping_index']['alias_to_gene_tokens']
    tokens = set()
    for gene in removed:
        candidates = aliases.get(gene, ())
        if len(candidates) != 1:
            raise ValueError('Replaced source gene has no unique native CDS owner: ' + gene)
        tokens.update(candidates)
    for gene in genes.keys() - removed:
        if tokens.intersection(aliases.get(gene, ())):
            raise ValueError('Replaced and retained genes share a native CDS owner: ' + gene)
    return tokens


def build_pair(task, directory, *, search_root=None, search_key=None, verification=None):
    import format_species_inputs as fsi
    from Bio.Seq import Seq
    from cloned_cds_search import alignments
    from gene_model_catalog import indexed_genome
    from rescue_gene_models import read_miniprot, validate_model
    code = int(task.get('genetic_code', 1))
    clones = cloned_records(Path(task['cloned_cds_path']), task.get('cloned_cds_ids', []), code)
    queries = directory / 'clone_proteins.fa'
    # Internal simple IDs avoid treating spaces, punctuation or GFF escapes in
    # collaborator-provided IDs as aligner/GFF syntax.
    names = {f'clone{i}': name for i, name in enumerate(sorted(clones), 1)}
    queries.write_text(''.join('>' + key + '\n' + str(Seq(clones[name][:-3]).translate(table=code)) + '\n'
                               for key, name in names.items()))
    cpus = max(1, int(os.environ.get('GG_TASK_CPUS', '1')))
    with indexed_genome(task['genome_path']) as genome:
        candidates = {name: {} for name in clones}
        diagnostics = {name: [] for name in clones}
        capped = set()
        for alignment in alignments(genome, queries, directory, code, cpus,
                                    task.get('cloned_cds_index_batch_bases', 'auto'), checkpoint_root=search_root,
                                    identity=search_key, verification=verification):
            before = {name: len(rows) for name, rows in diagnostics.items()}
            for raw in read_miniprot(alignment):
                if raw['query'] not in names:
                    raise ValueError('Unexpected miniprot query')
                name = names[raw['query']]
                model = validate_model(raw, genome, code, dict(minimum_coverage=1.0, minimum_identity=1.0,
                                                             minimum_broad_coverage=1.0, minimum_lower_coverage_identity=1.0))
                if model['sequence'] != clones[name]:
                    model['problems'].append('cloned_CDS_genome_mismatch')
                diagnostics[name].append({k: model[k] for k in ('seqid', 'strand', 'cds', 'coverage', 'identity', 'problems')})
                if not model['problems']:
                    key = (model['seqid'], model['strand'], tuple(map(tuple, model['cds'])))
                    candidates[name][key] = model
            capped.update(name for name, rows in diagnostics.items() if len(rows) - before[name] >= 100)
    audit = dict(schema_version=1, species=task['species_prefix'], genomic_CDS_identity_required=True,
                 candidates=diagnostics, replacements=[], skipped=[],
                 unresolved_policy=task.get('cloned_cds_unresolved_policy', 'error'),
                 search_limit_reached=sorted(capped))
    atomic_json(directory / 'audit.json', audit)
    selected_models = {}
    for name in sorted(clones):
        models = list(candidates[name].values())
        if len(models) != 1 or name in capped:
            if task.get('cloned_cds_unresolved_policy', 'error') == 'retain_source':
                audit['skipped'].append(dict(clone_id=name, status='unresolved_genomic_path',
                                             complete_genome_identical_paths=len(models),
                                             candidate_limit_reached=name in capped))
                continue
            raise ValueError(f'Cloned CDS {name} has {len(models)} complete genome-identical paths; inspect audit.json')
        selected_models[name] = models[0]
    records, genes, coding, ownership = source_graph(task['gff_path'])
    scoped_genes = overlap_gene_candidates(genes, coding, selected_models)
    identifiers = {attrs['ID'] for _, fields, attrs, _ in records if fields and attrs.get('ID')}
    removed = set()
    additions = []
    prepared = fsi.prepare_cds_identifier_task(task)
    longer = set(task.get('cloned_cds_longer_only_ids', []))
    if longer - clones.keys():
        raise ValueError('Longer-only policy names an unselected cloned CDS')
    source_lengths = {}
    if longer:
        for header, sequence in fsi.iter_task_cds_records(prepared):
            transcript = fsi.build_formatted_cds_id(prepared, header)
            _gene, match = fsi.resolve_gene_aggregate_id(prepared, header, transcript)
            if match.get('status') == 'mapped':
                token = match['gene_token']
                source_lengths[token] = max(source_lengths.get(token, 0), len(sequence))
    for name, model in sorted(selected_models.items()):
        scope = scoped_genes[name] if scoped_genes is not None else genes
        replaced = overlapping_genes(model, scope, coding)
        if name in longer:
            tokens = removed_cds_gene_tokens(prepared, set(replaced), genes)
            lengths = [source_lengths[t] for t in tokens if t in source_lengths]
            if not lengths or len(clones[name]) <= max(lengths):
                audit['skipped'].append(dict(clone_id=name, status='not_longer' if lengths else 'no_target_CDS',
                                             clone_cds_length=len(clones[name]), target_cds_lengths=lengths))
                continue
        if removed.intersection(replaced):
            raise ValueError('Competing cloned CDS would replace the same source gene')
        gene, transcript, text = replacement_gff(task['species_prefix'], name, model, replaced)
        new_ids = {gene, transcript} | {f'{transcript}.{kind}{i}' for kind in ('exon', 'CDS')
                                      for i in range(1, len(model['cds']) + 1)}
        if new_ids.intersection(identifiers):
            raise ValueError('Cloned model ID collides with existing annotation')
        identifiers.update(new_ids)
        additions.append((gene, transcript, clones[name], text))
        removed.update(replaced)
        audit['replacements'].append(dict(clone_id=name, gene_id=gene, transcript_id=transcript,
                                          replaced_gene_ids=replaced, seqid=model['seqid'], strand=model['strand'],
                                          CDS_intervals_0based_halfopen=model['cds']))
    for owned in ownership:
        if owned & removed and owned - removed:
            raise ValueError('A replaced GFF feature is shared with an unrelated retained gene')
    removed_tokens = removed_cds_gene_tokens(prepared, removed, genes)
    removed_cds = []
    cds = directory / Path(task['cds_path']).name
    gff = directory / Path(task['gff_path']).name
    import gzip
    def output(path):
        return gzip.open(path, 'wt') if path.name.endswith('.gz') else path.open('w')
    with output(cds) as handle:
        for header, sequence in fsi.iter_task_cds_records(prepared):
            transcript = fsi.build_formatted_cds_id(prepared, header)
            _gene_id, match = fsi.resolve_gene_aggregate_id(prepared, header, transcript)
            if match.get('status') != 'mapped':
                raise ValueError('Unresolved source CDS identity prevents coupled model replacement: ' + header.split()[0])
            if match.get('gene_token') in removed_tokens:
                removed_cds.append(header.split()[0])
            else:
                handle.write('>' + header + '\n' + sequence + '\n')
        for gene, transcript, sequence, _text in additions:
            handle.write('>' + transcript + ' gene:' + gene + '\n' + sequence + '\n')
    with output(gff) as handle:
        for (line, _fields, _attrs, _parents), owned in zip(records, ownership, strict=True):
            if not owned.intersection(removed):
                handle.write(line)
        for _gene, _transcript, _sequence, text in additions:
            handle.write(text)
    # A single-run download cache can be cleaned after publication. Retain the
    # exact reference with the coupled source pair so later native QC never
    # depends on that temporary cache. Copy in bounded buffers; a hard link
    # would change the original input's inode metadata and violate its fence.
    genome_source = directory / 'genome_source' / Path(task['genome_path']).name
    genome_source.parent.mkdir()
    shutil.copyfile(task['genome_path'], genome_source)
    # Every accepted path is genome-identical. Keep old annotations and raw CDS
    # outside this derived directory, bound in the stage receipt.
    audit.update(removed_CDS_ids=removed_cds, removed_CDS_gene_tokens=sorted(removed_tokens), original_sources_modified=False,
                 outputs={'cds': cds.name, 'gff': gff.name,
                          'genome': str(genome_source.relative_to(directory))})
    atomic_json(directory / 'audit.json', audit)


def apply_to_task(task, root, *, dry_run=False):
    if not task.get('cloned_cds_path'):
        return
    if task.get('_cloned_models_applied'):
        return
    if dry_run:
        raise ValueError('Cloned model discovery requires execution; use array_prepare for a download-free preview')
    if task.get('cloned_model_implementation') != digest(Path(__file__)):
        raise ValueError('Cloned model implementation changed; prepare a new native task plan')
    if task.get('cloned_cds_mode') == 'sequence':
        from cloned_cds_sequences import apply_sequence_task
        return apply_sequence_task(task, root)
    if task.get('cloned_cds_mode') == 'experimental':
        from cloned_cds_projection import apply, implementation
        if task.get('cloned_projection_implementation') != implementation():
            raise ValueError('Experimental projection implementation changed; prepare a new native task plan')
        if not all(task.get(k) for k in ('cds_path', 'gff_path', 'genome_path')):
            raise ValueError('Experimental cloned CDS requires source CDS, GFF and genome together')
        return apply(task, root)
    from cloned_cds_search import implementation
    if task.get('cloned_genomic_search_implementation') != implementation():
        raise ValueError('Cloned genomic search implementation changed; prepare a new native task plan')
    if not all(task.get(k) for k in ('cds_path', 'gff_path', 'genome_path')):
        raise ValueError('Cloned model replacement requires source CDS, GFF and genome together')
    from rescue_gene_models import stage
    tool = shutil.which('miniprot')
    if not tool:
        raise ValueError('miniprot is required for cloned CDS model replacement')
    sources = [Path(task[k]).resolve() for k in ('cds_path', 'gff_path', 'genome_path', 'cloned_cds_path')]
    batch = FreshDigestBatch()
    hashes = batch.read(sources)
    if hashes[str(sources[-1])] != task['cloned_cds_sha256']:
        raise ValueError('Cloned CDS changed after input qualification')
    key = dict(schema=1, inputs=hashes, species=task['species_prefix'],
               selected_ids=task.get('cloned_cds_ids', []), genetic_code=int(task.get('genetic_code', 1)),
               longer_only_ids=task.get('cloned_cds_longer_only_ids', []),
               unresolved_policy=task.get('cloned_cds_unresolved_policy', 'error'),
               index_batch_bases=task.get('cloned_cds_index_batch_bases', 'auto'),
               implementation=digest(Path(__file__)), miniprot=digest(tool),
               native_dependencies={name: digest(Path(__file__).parent / name) for name in (
                   'rescue_gene_models.py', 'gene_model_catalog.py', 'gene_model_refinement.py',
                   'rescue_homology_policy.py', 'cloned_cds_search.py', 'cloned_cds_resources.py',
                   'cloned_cds_checkpoints.py')})
    token = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()
    destination = stage(Path(root), task['species_prefix'] + '/' + token, key,
                        lambda folder: build_pair(task, folder, search_root=Path(root) / '.cloned_search' / token,
                                                  search_key=key, verification=batch),
                        guard=batch.check, publication_guard=batch.check)
    audit = json.loads((destination / 'audit.json').read_text())
    for role in ('cds', 'gff', 'genome'):
        task[role + '_path'] = destination / audit['outputs'][role]
    task['cloned_model_audit'] = str(destination / 'audit.json')
    task['_cloned_models_applied'] = True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['verify-manifest', 'evidence-files'])
    parser.add_argument('--manifest', type=Path)
    parser.add_argument('--metadata', type=Path)
    parser.add_argument('--species-summary', type=Path)
    args = parser.parse_args()
    try:
        if args.command == 'verify-manifest':
            if args.manifest is None:
                parser.error('--manifest is required')
            print(verify_manifest(args.manifest.resolve()))
        else:
            if (args.metadata is None) == (args.species_summary is None):
                parser.error('supply exactly one of --metadata and --species-summary')
            if args.metadata is not None:
                meta = json.loads(args.metadata.read_text())
                rows = [meta] if meta.get('cloned_cds_path') else []
            else:
                with args.species_summary.open() as stream:
                    rows = [row for row in csv.DictReader(stream, delimiter='\t') if row.get('cloned_model_audit_path')]
            seen = set()
            for row in rows:
                for value in evidence_paths(row, summary=args.species_summary is not None):
                    if value not in seen:
                        print(value)
                        seen.add(value)
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))


if __name__ == '__main__':
    sys.exit(main())
