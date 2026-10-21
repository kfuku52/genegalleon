"""Re-evaluate UTRs in a new publication without repeating coding selection.

Only consumed members of the completed source view are read and checked against
its frozen receipt. No source file, CDS, protein or BUSCO result is rewritten.
"""
import csv
import io
import json
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote, unquote

import annevo_gene_models
import rescue_gene_models as rescue
import transcript_utrs as utrs
from gene_model_catalog import indexed_genome
from input_generation_array_state import atomic_json, digest
from rescue_prepared_snapshot import _Proof


def _coding_rows(text):
    return [line for line in text.splitlines() if not line.startswith('#') and len(line.split('\t')) == 9
            and line.split('\t')[2].lower() == 'cds']


def _coding_seed(text, selected_ids=None):
    """Remove a preceding completion's derived ends before evaluating donors."""
    models = utrs.read_models(text, strict=False)
    replaced = {tid: model for tid, model in models.items()
                if model['attrs'].get('utr_policy') == ('projected_longest_per_side',)
                and (selected_ids is None or tid in selected_ids)}
    lines, embedded = [], False
    for line in text.splitlines(keepends=True):
        if embedded:
            lines.append(line)
            continue
        if line.startswith('##FASTA'):
            embedded = True
            lines.append(line)
            continue
        if line.startswith('#') or not line.strip():
            lines.append(line)
            continue
        fields = line.rstrip('\r\n').split('\t')
        attrs = utrs.parse_gff_attributes(fields[8])
        tid = next(iter(attrs.get('ID', ())), '')
        if fields[2].lower() in {'exon', 'utr', 'five_prime_utr', 'three_prime_utr', '5utr', '3utr'}:
            parents = attrs.get('Parent', ())
            remaining = [p for p in parents if p not in replaced]
            if len(remaining) != len(parents):
                if not remaining:
                    continue
                fields[8] = ';'.join('Parent='+','.join(quote(p, safe='._-:') for p in remaining)
                    if item.startswith('Parent=') else item for item in fields[8].split(';'))
        if tid in replaced and fields[2].lower().endswith(('rna', 'transcript')):
            cds = replaced[tid]['cds']
            fields[3], fields[4] = str(min(b[0] for b in cds)+1), str(max(b[1] for b in cds))
            fields[8] = ';'.join(item for item in fields[8].split(';') if not item.startswith('utr_policy='))
        lines.append('\t'.join(fields)+'\n')
    result = ''.join(lines)
    if _coding_rows(result) != _coding_rows(text):
        raise ValueError('UTR coding seed changed CDS rows')
    return result


def reproject(source_output, output, expected_receipt, species=None):
    """Use a receipted completed refinement as immutable coding-model input."""
    source_output, output = Path(source_output).resolve(), Path(output).resolve()
    source = source_output / 'effective'
    if output == source_output or output.is_relative_to(source_output) or source_output.is_relative_to(output):
        raise ValueError('UTR reassessment requires a separate output directory')
    proof = _Proof()
    receipt, receipt_sha = proof.json(source / 'receipt.json')
    if receipt_sha != expected_receipt or not isinstance(receipt.get('files'), dict) or not receipt['files']:
        raise ValueError('Completed source receipt changed or is invalid')
    plan, plan_sha = proof.json(source_output / 'plan.json')
    if receipt.get('key', {}).get('plan') != plan_sha:
        raise ValueError('Source publication does not belong to its refinement plan')

    consumed = {}
    def member(path):
        path = Path(path).absolute()
        if not path.is_relative_to(source) or path.resolve(strict=True) != path:
            raise ValueError('Source view member escaped its publication')
        relative = str(path.relative_to(source))
        expected = receipt['files'].get(relative)
        if expected is None:
            raise ValueError('Source view member is not receipted: ' + relative)
        actual = proof.read([path])[str(path)]
        if actual != expected:
            raise ValueError('Source view member changed: ' + relative)
        consumed[str(path)] = actual
        return path

    manifest = member(source / 'inputs.tsv')
    rows = list(csv.DictReader(io.StringIO(manifest.read_text()), delimiter='\t'))
    by_name = {row['species']: row for row in rows}
    names = list(plan['species']) if species is None else list(species)
    if (not names or len(set(names)) != len(names) or len(by_name) != len(rows)
            or not set(names) <= set(by_name) or not set(names) <= set(plan['species'])):
        raise ValueError('UTR species must be unique members of the source view')
    mapping = list(csv.DictReader(io.StringIO(member(source / 'representative_map.tsv').read_text()), delimiter='\t'))
    gene_ids = {}
    for row in mapping:
        key = row['species'], unquote(row['source_transcript_id'])
        if key in gene_ids and gene_ids[key] != row['gene_id']:
            raise ValueError('Source representative transcript has conflicting gene owners')
        gene_ids[key] = row['gene_id']
    inputs = {}
    for name in names:
        row = by_name[name]
        paths = {role: member(row[role]) for role in ('gff', 'genome', 'cds', 'analysis_gff')}
        paths.update({role: member(source / role / (name+'.gff3')) for role in ('full_annotation', 'source_annotation')})
        paths['utr_audit'] = member(source / 'utr_audit' / (name+'.json'))
        prediction = plan['request'].get('oriongeno', {}).get('predictions', {}).get(name)
        if prediction:
            for role in ('gff', 'receipt'):
                path = Path(prediction[role]).absolute()
                expected = plan['request']['files'].get(str(path))
                actual = proof.read([path])[str(path)]
                if actual != expected:
                    raise ValueError('Frozen native UTR input changed: ' + role)
                consumed[str(path)] = actual
            identity = annevo_gene_models.prediction_identity(prediction['gff'], prediction['receipt'],
                consumed[str(paths['genome'])], int(row['genetic_code']), str(paths['genome']), engine='oriongeno')
            for field in ('receipt_sha256', 'GFF_sha256', 'source_commit', 'model_sha256'):
                if identity[field] != prediction['identity'][field]:
                    raise ValueError('Native UTR prediction provenance changed: ' + field)
        inputs[name] = paths, prediction
    proof.check()
    code = {name: digest(Path(__file__).with_name(name)) for name in
            ('transcript_utr_reassessment.py', 'transcript_utrs.py')}
    key = dict(schema='genegalleon-utr-reassessment-v1', source_receipt=receipt_sha,
               source_plan=plan_sha, species=names, files=consumed, implementation=code)
    existing = output / 'effective' / 'receipt.json'
    if existing.exists() and json.loads(existing.read_text()).get('key') != key:
        raise ValueError('UTR output belongs to a different generation; use a new directory')

    def build(tmp):
        for role in ('species_gff', 'analysis_gff', 'full_annotation', 'species_transcript', 'utr_audit'):
            (tmp / role).mkdir()
        summaries = []
        for name in names:
            paths, prediction = inputs[name]
            selected_text = paths['gff'].read_text()
            native = utrs.load_models(prediction['gff']) if prediction else None
            with indexed_genome(str(paths['genome'])) as genome:
                text, audit, sequences = utrs.complete(_coding_seed(selected_text), paths['source_annotation'].read_text(), genome, name, native)
            selected = utrs.read_models(text, strict=False)
            updates = {r['transcript_id']: dict(model=selected[r['transcript_id']],
                        exons=selected[r['transcript_id']]['exons'], choices={'5': r['utr5'], '3': r['utr3']})
                       for r in audit if 'utr5' in r}
            for role, original, rendered in (
                    ('species_gff', selected_text, text),
                    ('analysis_gff', paths['analysis_gff'].read_text(), None),
                    ('full_annotation', paths['full_annotation'].read_text(), None)):
                rendered = utrs.render(_coding_seed(original, selected), updates) if rendered is None else rendered
                if _coding_rows(rendered) != _coding_rows(original):
                    raise ValueError('UTR reassessment changed CDS rows: ' + role)
                (tmp / role / (name+'.gff3')).write_text(rendered)
            with (tmp / 'species_transcript' / (name+'.fa')).open('w') as stream:
                for tid, sequence in sorted(sequences.items()):
                    stream.write('>'+quote(tid, safe='._-:')+' utr_policy=projected_longest_per_side\n'+sequence+'\n')
            audited = {r['transcript_id'] for r in audit}
            audit += [r for r in json.loads(paths['utr_audit'].read_text()) if r['transcript_id'] not in audited]
            final_genes = [line[1:].split()[0] for line in paths['cds'].read_text().splitlines() if line.startswith('>')]
            if len(final_genes) != len(set(final_genes)):
                raise ValueError('Duplicate final coding-model identity')
            final_genes = set(final_genes)
            with5, with3, normalized = set(), set(), set()
            reasons = defaultdict(set)
            for record in audit:
                gene = gene_ids.get((name, unquote(record['transcript_id'])), record['gene_id'])
                record['final_gene_id'], record['final_exported'] = gene, gene in final_genes
                if gene not in final_genes:
                    continue
                if record.get('utr5', {}).get('length'):
                    with5.add(gene)
                if record.get('utr3', {}).get('length'):
                    with3.add(gene)
                if (record.get('exon_path_basis') == 'coge_coding_plus_terminal_noncoding_segments'
                        or record.get('source_exon_interpretations')):
                    normalized.add(gene)
                reasons[record.get('native_correspondence', {}).get('status', record.get('status', 'unknown'))].add(gene)
            atomic_json(tmp / 'utr_audit' / (name+'.json'), audit)
            summaries.append(dict(species=name, final_gene_models=len(final_genes), with_5prime_utr=len(with5),
                with_3prime_utr=len(with3), with_both_utrs=len(with5 & with3), with_any_utr=len(with5 | with3),
                interpreted_coge_models=len(normalized), correspondence={k: len(v) for k, v in reasons.items()},
                cds_unchanged=True, cds_sha256=consumed[str(paths['cds'])]))
            print(json.dumps(summaries[-1], sort_keys=True), flush=True)
        atomic_json(tmp / 'summary.json', summaries)
        atomic_json(tmp / 'provenance.json', key)
    output.mkdir(parents=True, exist_ok=True)
    result = rescue.stage(output, 'effective', key, build, guard=proof.check, publication_guard=proof.check)
    proof.check()
    return result
