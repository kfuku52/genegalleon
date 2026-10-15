"""Project transcript ends onto selected CDSs and choose the longest UTR per side.

Coordinates are zero-based half-open. Source exon paths and outer transcript
ends stay fixed; a boundary change never translates a UTR along the genome.
"""
import bisect
import re
from collections import defaultdict
from urllib.parse import quote

from Bio.Seq import Seq
from format_species_annotation.common import parse_gff_attributes
from format_species_writers import open_text


def merge_blocks(blocks):
    result = []
    for a, b in sorted(set(tuple(block[:2]) for block in blocks)):
        if a < 0 or b <= a:
            raise ValueError('Invalid transcript interval')
        if result and a <= result[-1][1]:
            result[-1][1] = max(result[-1][1], b)
        else:
            result.append([a, b])
    return result


def read_models(text, *, strict=True):
    """Read explicit GFF3 or GTF transcript relationships, including shared rows."""
    models, children, genes = {}, defaultdict(list), set()
    for line in text.splitlines():
        if line == '##FASTA':
            break
        if not line or line.startswith('#'):
            continue
        f = line.split('\t')
        if len(f) != 9:
            raise ValueError('Invalid transcript GFF/GTF row')
        attrs = parse_gff_attributes(f[8])
        def scalar(key, attrs=attrs):
            return next(iter(attrs.get(key, ())), '')
        kind = f[2].lower()
        if kind.endswith('gene'):
            genes.add(scalar('ID'))
        if kind.endswith(('rna', 'transcript')):
            tid = scalar('ID') or scalar('transcript_id')
            parents = attrs.get('Parent', ()) or attrs.get('gene_id', ())
            if not tid or len(parents) != 1 or tid in models:
                raise ValueError('Transcript requires unique identity and one gene owner')
            models[tid] = dict(id=tid, gene=parents[0], seqid=f[0], strand=f[6],
                               cds=[], exons=[], utr=[], attrs=attrs)
        elif kind in {'exon', 'cds', 'utr', 'five_prime_utr', 'three_prime_utr', '5utr', '3utr'}:
            owners = attrs.get('Parent', ()) or attrs.get('transcript_id', ())
            for tid in owners:
                children[tid].append((kind, f, attrs))
    for tid, rows in children.items():
        if not any(kind == 'cds' for kind, _, _ in rows):
            continue
        if tid not in models:
            gene = next(iter(rows[0][2].get('gene_id', ())), '') or (tid if tid in genes else '')
            if not gene:
                raise ValueError('Unresolved transcript parent: ' + tid)
            models[tid] = dict(id=tid, gene=gene, seqid=rows[0][1][0], strand=rows[0][1][6],
                               cds=[], exons=[], utr=[], attrs={})
        model = models[tid]
        for kind, f, _ in rows:
            if (f[0], f[6]) != (model['seqid'], model['strand']) or f[6] not in {'+', '-'}:
                if strict:
                    raise ValueError('Transcript features disagree on contig or strand')
                model['problem'] = 'multiple_contigs_or_strands'
                continue
            block = [int(f[3]) - 1, int(f[4])]
            if block[0] < 0 or block[1] <= block[0]:
                raise ValueError('Invalid transcript coordinates')
            if kind == 'cds':
                model['cds'].append(block + [int(f[7]) if f[7] in {'0', '1', '2'} else None])
            else:
                model['exons' if kind == 'exon' else 'utr'].append(block)
        model['cds'] = sorted(set(tuple(b) for b in model['cds']), key=lambda b: (b[0], b[1], -1 if b[2] is None else b[2]))
        if any(a[1] > b[0] for a, b in zip(model['cds'], model['cds'][1:], strict=False)):
            if strict:
                raise ValueError('Overlapping transcript CDS blocks')
            model['problem'] = 'overlapping_cds'
        model['exons'] = merge_blocks(model['exons'] or [*model['cds'], *model['utr']])
        for a, b, _ in model['cds']:
            if not any(x <= a < b <= y for x, y in model['exons']):
                if strict:
                    raise ValueError('CDS is not contained in its transcript exons')
                model['problem'] = 'cds_outside_exons'
        if any(not any(x <= a < b <= y for x, y in model['exons']) for a, b in model['utr']):
            if strict:
                raise ValueError('UTR is not contained in its transcript exons')
            model['problem'] = 'utr_outside_exons'
    return {tid: m for tid, m in models.items() if m['cds']}


def load_models(path):
    with open_text(path, 'rt') as stream:
        return read_models(stream.read())


def project(donor, selected, side, genome):
    """Retain donor exons outside selected CDS; inspect any new terminal splice."""
    if donor['seqid'] not in genome.references or any(b > genome.get_reference_length(donor['seqid']) for _, b in donor['exons']):
        return dict(blocks=[], length=0, reason='source_outside_genome', boundary_shift=0, new_junction=False)
    low = (side == '5') == (selected['strand'] == '+')
    old_boundary = min(b[0] for b in donor['cds']) if low else max(b[1] for b in donor['cds'])
    boundary = min(b[0] for b in selected['cds']) if low else max(b[1] for b in selected['cds'])
    outer = donor['exons'][0][0] if low else donor['exons'][-1][1]
    if (old_boundary <= outer if low else old_boundary >= outer):
        return dict(blocks=[], length=0, reason='no_source_utr', boundary_shift=0, new_junction=False)
    blocks = [[a, min(b, boundary)] for a, b in donor['exons'] if a < boundary] if low else [
        [max(a, boundary), b] for a, b in donor['exons'] if b > boundary]
    blocks = [b for b in blocks if b[0] < b[1]]
    reason, new_junction = 'projected', False
    if not blocks:
        reason = 'cds_reaches_transcript_end'
    else:
        edge = blocks[-1][1] if low else blocks[0][0]
        if edge != boundary:
            # The CDS boundary lies in a donor intron or outside its exon path.
            # A synthetic splice is allowed only with a genomic splice motif.
            a, b = sorted((edge, boundary))
            intron = genome.fetch(selected['seqid'], a, b).upper()
            if selected['strand'] == '-':
                intron = str(Seq(intron).reverse_complement())
            if len(intron) < 4 or (intron[:2], intron[-2:]) not in {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')}:
                blocks, reason = [], 'terminal_connection_unresolved'
            else:
                new_junction = True
    return dict(blocks=blocks, length=sum(b-a for a, b in blocks), reason=reason,
                boundary_shift=boundary-old_boundary, new_junction=new_junction)


def _overlap(a, b):
    return sum(max(0, min(y, v)-max(x, u)) for x, y, *_ in a for u, v, *_ in b)


def complete(selected_text, original_text, genome, species, native=None):
    """Return selected GFF, per-transcript audit and spliced transcript sequences."""
    selected, original = read_models(selected_text, strict=False), read_models(original_text, strict=False)
    native_available = native is not None
    native = native or {}
    donors, original_by_gene, axes = defaultdict(list), defaultdict(list), defaultdict(list)
    for m in original.values():
        if not m.get('problem'):
            original_by_gene[m['gene']].append(m)
    for m in selected.values():
        axes[m['seqid'], m['strand']].append(m)
    tracks = {}
    for axis, rows in axes.items():
        rows.sort(key=lambda m: min(b[0] for b in m['cds']))
        starts, ends = [], []
        for m in rows:
            starts.append(min(b[0] for b in m['cds']))
            ends.append(max(ends[-1] if ends else 0, max(b[1] for b in m['cds'])))
        tracks[axis] = starts, ends, rows
    for model in native.values():
        if model['seqid'] not in genome.references or any(b > genome.get_reference_length(model['seqid']) for _, b in model['exons']):
            raise ValueError('OrionGeno transcript lies outside target genome')
        starts, ends, rows = tracks.get((model['seqid'], model['strand']), ([], [], []))
        left = bisect.bisect_right(ends, min(b[0] for b in model['cds']))
        right = bisect.bisect_left(starts, max(b[1] for b in model['cds']))
        overlapping = [m for m in rows[left:right] if _overlap(m['cds'], model['cds'])]
        shape = [b[:2] for b in model['cds']]
        exact = {m['gene'] for m in overlapping if [b[:2] for b in m['cds']] == shape}
        owners = exact or {m['gene'] for m in overlapping}
        if len(owners) == 1:
            donors[next(iter(owners))].append(model)
    updates, audit, sequences = {}, [], {}
    for tid, model in selected.items():
        if model.get('problem'):
            audit.append(dict(species=species, transcript_id=tid, gene_id=model['gene'],
                              status='unsupported_structure', reason=model['problem']))
            continue
        candidates = [('existing', m) for m in original_by_gene[model['gene']]
                      if (m['seqid'], m['strand']) == (model['seqid'], model['strand'])]
        candidates += [('oriongeno', m) for m in donors[model['gene']]]
        # A newly selected native/predicted transcript may not exist in source GFF.
        candidates += [('selected', model)]
        choices, evaluations = {}, {}
        for side in ('5', '3'):
            rows = []
            for source, donor in candidates:
                p = project(donor, model, side, genome)
                rows.append(dict(source=source, source_transcript=donor['id'], source_gene=donor['gene'], **p))
            ranked = sorted(rows, key=lambda r: (-r['length'], {'existing': 0, 'selected': 1, 'oriongeno': 2}[r['source']], r['source_transcript']))
            choices[side] = ranked[0] if ranked else dict(blocks=[], length=0, reason='no_candidate')
            evaluations[side] = rows
        exons = merge_blocks([*model['cds'], *choices['5']['blocks'], *choices['3']['blocks']])
        pieces = [genome.fetch(model['seqid'], a, b) for a, b in exons]
        sequence = ''.join(pieces).upper()
        if model['strand'] == '-':
            sequence = str(Seq(sequence).reverse_complement())
        sequences[tid] = sequence
        updates[tid] = dict(model=model, choices=choices, exons=exons)
        audit.append(dict(species=species, transcript_id=tid, gene_id=model['gene'],
                          utr5=choices['5'], utr3=choices['3'], candidates=evaluations,
                          native_status='mapped' if donors[model['gene']] else 'no_unique_coding_overlap' if native_available else 'missing_prediction',
                          annotation_basis='projected_longest_per_side'))
    return render(selected_text, updates), audit, sequences


def render(text, updates):
    """Update structural rows only; retain CDS IDs, coordinates, phases and attributes."""
    directive = re.search(r'^##FASTA(?:\r?\n|$)', text, flags=re.MULTILINE)
    text, marker, embedded = (text[:directive.start()], '##FASTA', text[directive.start()+7:]) if directive else (text, '', '')
    lines, gene_bounds, retained_features = [], defaultdict(list), set()
    for value in updates.values():
        gene_bounds[value['model']['gene']].extend(value['exons'])
    for line in text.splitlines(keepends=True):
        if not line.strip() or line.startswith('#'):
            lines.append(line)
            continue
        f = line.rstrip('\r\n').split('\t')
        attrs = parse_gff_attributes(f[8])
        tid = next(iter(attrs.get('ID', ())), '')
        parents = attrs.get('Parent', ())
        kind = f[2].lower()
        if kind in {'exon', 'utr', 'five_prime_utr', 'three_prime_utr', '5utr', '3utr'} and set(parents) & updates.keys():
            block = [int(f[3])-1, int(f[4])]
            remaining = []
            for parent in parents:
                if parent not in updates:
                    remaining.append(parent)
                    continue
                value = updates[parent]
                side = '5' if kind in {'five_prime_utr', '5utr'} else '3' if kind in {'three_prime_utr', '3utr'} else None
                expected = (value['exons'] if kind == 'exon' else value['choices'][side]['blocks'] if side
                            else [b for c in value['choices'].values() for b in c['blocks']])
                if block in expected:
                    remaining.append(parent)
                    retained_features.add((parent, 'exon' if kind == 'exon' else 'utr', *block))
            if not remaining:
                continue
            f[8] = ';'.join('Parent='+','.join(quote(p, safe='._-:') for p in remaining) if p.startswith('Parent=') else p for p in f[8].split(';'))
        if kind.endswith('gene') and tid in gene_bounds:
            # Full annotation can retain other isoforms outside selected bounds.
            f[3], f[4] = str(min(int(f[3]), min(a for a, _ in gene_bounds[tid])+1)), str(max(int(f[4]), max(b for _, b in gene_bounds[tid])))
        if tid in updates and kind.endswith(('rna', 'transcript')):
            value = updates[tid]
            f[3], f[4] = str(value['exons'][0][0]+1), str(value['exons'][-1][1])
            f[8] = ';'.join(p for p in f[8].split(';') if not p.startswith('utr_policy=')) + ';utr_policy=projected_longest_per_side'
        lines.append('\t'.join(f)+'\n')
    used = {next(iter(parse_gff_attributes(line.split('\t')[8]).get('ID', ())), '')
            for line in lines if line.strip() and not line.startswith('#')}
    for tid, value in updates.items():
        model = value['model']
        def emit(kind, block, suffix, extra='', tid=tid, model=model):
            identifier = tid+'.ggutr.'+suffix
            while identifier in used:
                identifier += '_'
            used.add(identifier)
            attrs = 'ID='+quote(identifier, safe='._-:')+';Parent='+quote(tid, safe='._-:')+extra
            lines.append(f"{model['seqid']}\tGeneGalleon\t{kind}\t{block[0]+1}\t{block[1]}\t.\t{model['strand']}\t.\t{attrs}\n")
        for i, block in enumerate(value['exons']):
            if (tid, 'exon', *block) not in retained_features:
                emit('exon', block, 'exon'+str(i+1))
        for side, choice in value['choices'].items():
            extra = ';annotation_basis=projected_longest;utr_source='+choice.get('source', 'none')+';utr_source_transcript='+quote(choice.get('source_transcript', ''), safe='._-:')
            for i, block in enumerate(choice['blocks']):
                if (tid, 'utr', *block) not in retained_features:
                    emit('five_prime_UTR' if side == '5' else 'three_prime_UTR', block, side+'.'+str(i+1), extra)
    return ''.join(lines) + (marker + embedded if marker else '')
