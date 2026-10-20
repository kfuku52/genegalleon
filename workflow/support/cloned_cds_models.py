#!/usr/bin/env python3
"""Replace coordinate-supported source models with explicitly supplied cloned CDS.

This is a coupled source curation step before ordinary native formatting. It
never changes a genome or accepts a CDS-only override of genomic annotation.
"""
import argparse
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

FIELDS = ('cloned_cds_fasta', 'cloned_cds_sha256', 'cloned_cds_ids')
TASK_FIELDS = ('cloned_cds_path', 'cloned_cds_sha256', 'cloned_cds_ids', 'cloned_model_implementation')
EVIDENCE_FIELDS = ('cloned_cds_path', 'cloned_model_audit_path', 'cloned_source_cds_path',
                   'cloned_source_gff_path', 'cloned_source_genome_path')


def attach(task, row, parent):
    """Bind only a local, hash-qualified FASTA; a mixed FASTA needs explicit IDs."""
    if not any(str(row.get(k) or '').strip() for k in FIELDS):
        return task
    reference = str(row.get('cloned_cds_fasta') or '').strip()
    expected = str(row.get('cloned_cds_sha256') or '').strip()
    if not reference or not re.fullmatch(r'[0-9a-f]{64}', expected):
        raise ValueError('Cloned CDS requires cloned_cds_fasta and its exact cloned_cds_sha256')
    if row.get('paired_curation') or str(row.get('cds_input_type') or 'cds').strip().lower() != 'cds':
        raise ValueError('Cloned model replacement requires genomic CDS/GFF inputs, without paired_curation')
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
    task.update(cloned_cds_path=path, cloned_cds_sha256=expected, cloned_cds_ids=ids,
                cloned_model_implementation=digest(Path(__file__)))
    return task


def copy_binding(task, source):
    for key in TASK_FIELDS:
        if key in source:
            task[key] = source[key]
    return task


def binding_signature(task):
    if not task.get('cloned_cds_path'):
        return ''
    binding = {k: str(task[k]) if isinstance(task[k], Path) else task[k] for k in TASK_FIELDS}
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
        sequence = sequence.upper()
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
        # GTF can require synthesized gene/transcript rows. Ordinary GFF3 is
        # already normalized by the native parser; avoid rendering and parsing
        # its entire annotation a second time.
        records, parents, _ = _gff_source_graph(_selected_gff_source_rows(graph, set(), retain_all=True))
    else:
        for _, fields, attrs, pids in records:
            if fields and pids:
                attrs.setdefault('Parent', pids[0])
        records = [('##gff-version 3\n', None, {}, [])] + [
            row for row in records if not row[0].startswith('##gff-version')]
    del graph, implicit, known
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
            owned = owners(attrs['ID']) if attrs.get('ID') else frozenset(g for p in parent_ids for g in owners(p))
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


def build_pair(task, directory):
    import format_species_inputs as fsi
    from Bio.Seq import Seq
    from gene_model_catalog import indexed_genome
    from rescue_gene_models import read_miniprot, run, validate_model
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
        run(['miniprot', '-u', '-T', code, '-t', cpus, '-N', 100, '--outs=0', '--gff', genome.filename.decode(), queries],
            directory, 'cloned_cds_miniprot', directory / 'models.gff')
        candidates = {name: {} for name in clones}
        diagnostics = {name: [] for name in clones}
        for raw in read_miniprot(directory / 'models.gff'):
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
    audit = dict(schema_version=1, species=task['species_prefix'], genomic_CDS_identity_required=True,
                 candidates=diagnostics, replacements=[])
    atomic_json(directory / 'audit.json', audit)
    records, genes, coding, ownership = source_graph(task['gff_path'])
    identifiers = {attrs['ID'] for _, fields, attrs, _ in records if fields and attrs.get('ID')}
    removed = set()
    additions = []
    for name in sorted(clones):
        models = list(candidates[name].values())
        if len(models) != 1 or len(diagnostics[name]) >= 100:
            raise ValueError(f'Cloned CDS {name} has {len(models)} complete genome-identical paths; inspect audit.json')
        model = models[0]
        replaced = overlapping_genes(model, genes, coding)
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
    prepared = fsi.prepare_cds_identifier_task(task)
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
    if not all(task.get(k) for k in ('cds_path', 'gff_path', 'genome_path')):
        raise ValueError('Cloned model replacement requires source CDS, GFF and genome together')
    if task.get('cloned_model_implementation') != digest(Path(__file__)):
        raise ValueError('Cloned model implementation changed; prepare a new native task plan')
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
               implementation=digest(Path(__file__)), miniprot=digest(tool),
               native_dependencies={name: digest(Path(__file__).parent / name) for name in (
                   'rescue_gene_models.py', 'gene_model_catalog.py', 'gene_model_refinement.py',
                   'rescue_homology_policy.py')})
    token = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()
    destination = stage(Path(root), task['species_prefix'] + '/' + token, key,
                        lambda folder: build_pair(task, folder), guard=batch.check, publication_guard=batch.check)
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
                fields = EVIDENCE_FIELDS
            else:
                with args.species_summary.open() as stream:
                    rows = [row for row in csv.DictReader(stream, delimiter='\t') if row.get('cloned_model_audit_path')]
                fields = ('cloned_model_audit_path', 'cds_input_path', 'gff_input_path', 'genome_input_path')
            seen = set()
            for row in rows:
                for field in fields:
                    value = row.get(field)
                    if (not value or not Path(value).is_file() or not Path(value).stat().st_size
                            or any(c in value for c in '\r\n')):
                        raise ValueError('Missing or unsafe cloned model evidence: ' + field)
                    if value not in seen:
                        print(value)
                        seen.add(value)
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))


if __name__ == '__main__':
    sys.exit(main())
