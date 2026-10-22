"""Explicit UTRs, compressed audits and validation of native prediction outputs.

UTRs are exonic intervals outside the unchanged CDS endpoints. This export
does not predict sequence, extend a transcript, or choose among source models.
"""
import gzip
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import quote

try:
    from annevo_runtime_download import file_sha
    from format_species_annotation.common import parse_gff_attributes
    from transcript_utrs import read_models
except ImportError:
    from .annevo_runtime_download import file_sha
    from .format_species_annotation.common import parse_gff_attributes
    from .transcript_utrs import read_models

UTR_POLICY = 'exon_cds_boundaries_v1'
UTR_KINDS = {'utr', 'five_prime_utr', 'three_prime_utr', '5utr', '3utr'}
AUDIT_SUFFIXES = {'hmm_boundaries': '.hmm_boundaries.jsonl', 'utr_repairs': '.utr_repairs.jsonl'}


def companion_paths(output):
    output = Path(output)
    return {'native_gff': output.with_name(output.stem + '.native.gff3'),
            'repeat_gff': output.with_name(output.stem + '.repeat.gff3'),
            **{name: output.with_name(output.stem + suffix + '.gz') for name, suffix in AUDIT_SUFFIXES.items()}}


def compress_audit_jsonl(source, output):
    """Preserve exact JSONL bytes in deterministic gzip, retaining one line at a time."""
    source = Path(source)
    if source.is_symlink() or not source.is_file():
        raise ValueError('OrionGeno audit log is missing or nonregular: ' + source.name)
    sha, size, records = hashlib.sha256(), 0, 0
    with source.open('rb') as stream, Path(output).open('wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=6) as compressed:
            for line in stream:
                try:
                    value = json.loads(line)
                except (ValueError, UnicodeError) as error:
                    raise ValueError('Invalid OrionGeno audit JSONL: ' + source.name) from error
                if not isinstance(value, dict):
                    raise ValueError('OrionGeno audit JSONL records must be objects: ' + source.name)
                compressed.write(line)
                sha.update(line)
                size += len(line)
                records += 1
    return dict(compression='gzip', uncompressed_sha256=sha.hexdigest(), uncompressed_bytes=size, jsonl_records=records)


def _terminal_blocks(model):
    start, end = model['cds'][0][0], model['cds'][-1][1]
    low = [(a, min(b, start)) for a, b in model['exons'] if a < start]
    high = [(max(a, end), b) for a, b in model['exons'] if b > end]
    return {'five_prime_UTR': low, 'three_prime_UTR': high} if model['strand'] == '+' else {
        'five_prime_UTR': high, 'three_prime_UTR': low}


def _uncovered(block, covered):
    start, end = block
    for left, right in sorted(covered):
        if right <= start:
            continue
        if left >= end:
            break
        if left > start:
            yield start, left
        start = max(start, right)
        if start >= end:
            return
    if start < end:
        yield start, end


def write_explicit_utrs(source, output, *, lengths=None):
    """Preserve every source row and add only missing terminal exonic UTR rows."""
    text = Path(source).read_text()
    directive = re.search(r'^##FASTA(?:\r?\n|$)', text, flags=re.MULTILINE)
    annotation, separator, embedded = (text[:directive.start()], '##FASTA', text[directive.start()+7:]) if directive else (text, '', '')
    models = read_models(annotation)
    expected, covered, used, encoded_ids = {}, {}, set(), {}
    for tid, model in models.items():
        if lengths is not None and (model['seqid'] not in lengths or any(
                b > lengths[model['seqid']] for _, b in model['exons'])):
            raise ValueError('OrionGeno transcript lies outside the genome')
        expected[tid] = _terminal_blocks(model)
    for line in annotation.splitlines():
        if not line or line.startswith('#'):
            continue
        row = line.split('\t')
        attrs = parse_gff_attributes(row[8])
        used.update(attrs.get('ID', ()))
        original_attrs = dict(item.split('=', 1) for item in row[8].split(';') if '=' in item)
        if len(attrs.get('ID', ())) == 1 and original_attrs.get('ID'):
            encoded_ids[attrs['ID'][0]] = original_attrs['ID']
        if row[2].lower() not in UTR_KINDS:
            continue
        block = int(row[3]) - 1, int(row[4])
        for tid in attrs.get('Parent', ()):
            if tid not in models:
                continue
            kind = next((kind for kind, blocks in expected[tid].items()
                         if any(a <= block[0] < block[1] <= b for a, b in blocks)), None)
            declared = {'5utr': 'five_prime_UTR', 'five_prime_utr': 'five_prime_UTR',
                        '3utr': 'three_prime_UTR', 'three_prime_utr': 'three_prime_UTR'}.get(row[2].lower())
            if kind is None or declared is not None and kind != declared:
                raise ValueError('Explicit UTR conflicts with transcript exon/CDS boundaries')
            covered.setdefault((tid, kind), []).append(block)
    additions, counts = [], {'five_prime_UTR': 0, 'three_prime_UTR': 0}
    for tid, model in models.items():
        for kind, blocks in expected[tid].items():
            index = 0
            for block in blocks:
                for start, end in _uncovered(block, covered.get((tid, kind), ())):
                    index += 1
                    identifier = tid + '.ggutr.' + ('5.' if kind == 'five_prime_UTR' else '3.') + str(index)
                    while identifier in used:
                        identifier += '_'
                    used.add(identifier)
                    attrs = ('ID=' + quote(identifier, safe='._-:') + ';Parent=' + encoded_ids.get(tid, quote(tid, safe='._-:'))
                             + ';annotation_basis=' + UTR_POLICY)
                    additions.append('\t'.join([model['seqid'], 'GeneGalleon', kind, str(start + 1),
                                                str(end), '.', model['strand'], '.', attrs]) + '\n')
                    counts[kind] += 1
    rendered = annotation + ('\n' if additions and annotation and not annotation.endswith('\n') else '') + ''.join(additions)
    if separator:
        rendered += separator + embedded
    Path(output).write_text(rendered)
    return {'policy': UTR_POLICY, 'transcripts': len(models), 'added_features': counts}


def validate_repeat_gff(path, lengths):
    """Allow an empty prediction, but require GFF3 and genome-bounded intervals."""
    identifiers, count, header = set(), 0, False
    with Path(path).open() as stream:
        for line in stream:
            if line.rstrip() == '##gff-version 3':
                header = True
            if line.startswith('#') or not line.strip():
                continue
            row = line.rstrip('\n').split('\t')
            if (len(row) != 9 or row[2] != 'repeat_region' or row[0] not in lengths
                    or row[6] not in {'+', '-', '.'} or row[7] != '.'
                    or not 1 <= int(row[3]) <= int(row[4]) <= lengths[row[0]]):
                raise ValueError('Invalid OrionGeno repeat GFF coordinates or feature')
            attrs = parse_gff_attributes(row[8])
            ids = attrs.get('ID', ())
            if len(ids) != 1 or ids[0] in identifiers or attrs.get('Parent'):
                raise ValueError('Invalid OrionGeno repeat GFF identity')
            identifiers.add(ids[0])
            count += 1
    if not header:
        raise ValueError('OrionGeno repeat output lacks its GFF3 header')
    return {'repeat_region': count}


def artifact_record(staged, destination):
    return {'path': Path(destination).name, 'sha256': file_sha(staged),
            'bytes': Path(staged).stat().st_size}


def verify_companions(value, gff):
    """New receipts bind companions; historical receipts retain their old scope."""
    outputs = value.get('annotation_outputs')
    if outputs is None:
        if any(key in value.get('inference_settings', {}) for key in ('output_utrs', 'output_repeat')):
            raise ValueError('OrionGeno annotation output contract is missing')
        return
    if (not isinstance(outputs, dict) or type(outputs.get('schema_version')) is not int or outputs['schema_version'] not in (1, 2)
            or type(outputs.get('output_repeat')) is not bool or type(outputs.get('output_utrs')) is not bool
            or outputs.get('utr_policy') != (UTR_POLICY if outputs['output_utrs'] else 'native')):
        raise ValueError('Invalid OrionGeno annotation output contract')
    settings = value.get('inference_settings', {})
    if any(settings.get(key) is not outputs[key] for key in ('output_utrs', 'output_repeat')):
        raise ValueError('OrionGeno annotation output flags differ from inference settings')
    paths = companion_paths(gff)
    names = {'native_gff'} | ({'repeat_gff'} if outputs['output_repeat'] else set())
    if outputs['schema_version'] == 2:
        if outputs.get('audit_log_compression') != 'gzip' or settings.get('audit_log_compression') != 'gzip':
            raise ValueError('OrionGeno audit output contract is missing or inconsistent')
        names.update(AUDIT_SUFFIXES)
    elif 'audit_log_compression' in outputs or 'audit_log_compression' in settings:
        raise ValueError('OrionGeno audit outputs require their versioned contract')
    if not isinstance(outputs.get('files'), dict) or set(outputs['files']) != names:
        raise ValueError('OrionGeno annotation companions are incomplete')
    for name, record in outputs['files'].items():
        path = paths[name]
        if (not isinstance(record, dict) or record.get('path') != path.name
                or not path.is_file() or path.is_symlink() or record.get('bytes') != path.stat().st_size
                or record.get('sha256') != file_sha(path)):
            raise ValueError('OrionGeno annotation companion changed: ' + name)
        if name in AUDIT_SUFFIXES and (
                record.get('compression') != 'gzip'
                or not re.fullmatch('[0-9a-f]{64}', str(record.get('uncompressed_sha256')))
                or any(type(record.get(key)) is not int or record[key] < 0 for key in ('uncompressed_bytes', 'jsonl_records'))):
            raise ValueError('Invalid OrionGeno audit output metadata: ' + name)
