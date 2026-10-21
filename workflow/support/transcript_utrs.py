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


def _coge_terminal_segments(model, rows):
    """Recognize CoGe's separate coding and terminal noncoding segments.

    This is a source-specific interpretation, not repair of arbitrary incomplete
    exon paths. CoGe feature identity, unique Parent links and declared RNA bounds
    must independently support the union. Internal or overlapping segments fail.
    """
    fid = model['attrs'].get('coge_fid', ())
    if (model.get('problem') or model.get('source', '').lower() != 'coge'
            or len(fid) != 1 or not fid[0].isdecimal() or not model['exons']
            or model['attrs'].get('ID') != (model['id'],)
            or model['attrs'].get('Parent') != (model['gene'],)):
        return None
    low, high = model['bounds']
    coding_low, coding_high = model['cds'][0][0], model['cds'][-1][1]
    for kind, fields, attrs in rows:
        if (fields[1].lower() != 'coge' or attrs.get('Parent') != (model['id'],)
                or len(attrs.get('coge_fid', ())) != 1 or not attrs['coge_fid'][0].isdecimal()):
            return None
        if kind == 'exon' and attrs['coge_fid'] != fid:
            return None
    if any(not (b <= coding_low or a >= coding_high) for a, b in model['exons']):
        return None
    blocks = [*model['exons'], *[b[:2] for b in model['cds']], *model['utr']]
    if any(not low <= a < b <= high for a, b in blocks):
        return None
    result = merge_blocks(blocks)
    return result if (result[0][0], result[-1][1]) == (low, high) else None


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
            # Standalone noncoding RNAs need no coding-locus owner. Validate
            # ownership only when this transcript actually has CDS children.
            if not tid:
                continue
            if tid in models:
                models[tid]['problem'] = 'duplicate_transcript_identity'
                continue
            models[tid] = dict(id=tid, gene=parents[0] if len(parents) == 1 else '', seqid=f[0], strand=f[6],
                               cds=[], exons=[], utr=[], attrs=attrs, source=f[1],
                               bounds=(int(f[3])-1, int(f[4])))
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
        if not model['gene']:
            model['problem'] = 'unresolved_or_ambiguous_gene_owner'
        if strict and model.get('problem'):
            raise ValueError('Transcript requires unique identity and one gene owner: ' + tid)
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
        if any(not any(x <= a < b <= y for x, y in model['exons']) for a, b, _ in model['cds']):
            interpreted = _coge_terminal_segments(model, rows)
            if interpreted is not None:
                model['declared_exons'] = model['exons']
                model['exons'] = interpreted
                model['exon_path_basis'] = 'coge_coding_plus_terminal_noncoding_segments'
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


def project(donor, selected, side, genome, *, reference_lengths=None):
    """Retain donor exons outside selected CDS; audit any new terminal splice."""
    if reference_lengths is None:
        length = genome.get_reference_length(donor['seqid']) if donor['seqid'] in genome.references else None
    else:
        length = reference_lengths.get(donor['seqid'])
    if length is None or any(b > length for _, b in donor['exons']):
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
    reason, new_junction, motif, canonical = 'projected', False, None, None
    if not blocks:
        reason = 'cds_reaches_transcript_end'
    else:
        edge = blocks[-1][1] if low else blocks[0][0]
        if edge != boundary:
            # The CDS boundary lies in a donor intron or outside its exon path.
            # This is the requested coordinate projection, including novel
            # terminal connections. Record their motifs without imposing an
            # additional UTR-loss gate or claiming biological splice support.
            a, b = sorted((edge, boundary))
            intron = (genome.fetch(selected['seqid'], a, b) if b-a < 4 else
                      genome.fetch(selected['seqid'], a, a+2)+genome.fetch(selected['seqid'], b-2, b)).upper()
            if selected['strand'] == '-':
                intron = str(Seq(intron).reverse_complement())
            new_junction, motif = True, intron[:2]+'-'+intron[-2:]
            canonical = b-a >= 4 and (intron[:2], intron[-2:]) in {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')}
    return dict(blocks=blocks, length=sum(b-a for a, b in blocks), reason=reason,
                boundary_shift=boundary-old_boundary, new_junction=new_junction,
                junction_motif=motif, junction_canonical=canonical)


def _overlap(a, b):
    return sum(max(0, min(y, v)-max(x, u)) for x, y, *_ in a for u, v, *_ in b)


def _coding_index(models):
    axes, tracks = defaultdict(list), {}
    for model in models:
        axes[model['seqid'], model['strand']].append(model)
    for axis, rows in axes.items():
        rows.sort(key=lambda m: min(b[0] for b in m['cds']))
        starts, ends = [], []
        for model in rows:
            starts.append(min(b[0] for b in model['cds']))
            ends.append(max(ends[-1] if ends else 0, max(b[1] for b in model['cds'])))
        tracks[axis] = starts, ends, rows
    return tracks


def _coding_overlaps(tracks, model, *, opposite=False):
    strand = ('-' if model['strand'] == '+' else '+') if opposite else model['strand']
    starts, ends, rows = tracks.get((model['seqid'], strand), ([], [], []))
    left = bisect.bisect_right(ends, min(b[0] for b in model['cds']))
    right = bisect.bisect_left(starts, max(b[1] for b in model['cds']))
    return [m for m in rows[left:right] if _overlap(m['cds'], model['cds'])]


def _coding_introns(model):
    blocks = merge_blocks(model['cds'])
    return tuple((a[1], b[0]) for a, b in zip(blocks, blocks[1:], strict=False))


def native_correspondence(selected, native):
    """Resolve exact CDS, then unique full intron chain, then unique overlap.

    A nonempty identical chain is structural evidence for one locus even if a
    terminal coding extension also touches a neighbour. Shared chains, single
    exons and genuinely different multi-locus paths remain ambiguous.
    """
    tracks = _coding_index(selected.values())
    donors, evidence = defaultdict(list), defaultdict(list)
    for model in native.values():
        overlapping = _coding_overlaps(tracks, model)
        exact = {m['gene'] for m in overlapping if [b[:2] for b in m['cds']] == [b[:2] for b in model['cds']]}
        introns = _coding_introns(model)
        chain = {m['gene'] for m in overlapping if introns and _coding_introns(m) == introns}
        # Multiple exact owners must not be resolved by weaker evidence.
        owners = exact or (chain if len(chain) == 1 else {m['gene'] for m in overlapping})
        basis = ('exact_cds' if exact else 'identical_coding_intron_chain' if len(chain) == 1
                 else 'unique_coding_overlap')
        if len(owners) == 1:
            donors[next(iter(owners))].append(model)
        record = dict(transcript_id=model['id'], basis=basis if len(owners) == 1 else 'ambiguous_coding_overlap',
                      assigned_gene=next(iter(owners)) if len(owners) == 1 else None,
                      overlapping_genes=sorted({m['gene'] for m in overlapping}))
        for gene in record['overlapping_genes']:
            evidence[gene].append(record)
    native_tracks = _coding_index(native.values())
    native_contigs = {m['seqid'] for m in native.values()}
    diagnostics = {}
    for tid, model in selected.items():
        candidates = evidence[model['gene']]
        reason = ('mapped' if donors[model['gene']] else
                  ('ambiguous_coding_overlap' if any(c['assigned_gene'] is None for c in candidates)
                   else 'overlap_assigned_to_other_gene') if candidates
                  else 'opposite_strand_coding_overlap' if _coding_overlaps(native_tracks, model, opposite=True)
                  else 'no_coding_overlap' if model['seqid'] in native_contigs else 'no_native_prediction_on_contig')
        diagnostics[tid] = dict(status=reason, candidates=candidates)
    return donors, diagnostics


def complete(selected_text, original_text, genome, species, native=None):
    """Return selected GFF, per-transcript audit and spliced transcript sequences."""
    selected, original = read_models(selected_text, strict=False), read_models(original_text, strict=False)
    native_available = native is not None
    native = native or {}
    # pysam exposes all reference names as a tuple. Repeated membership scans
    # for every donor/end are costly on fragmented assemblies. Keep this index
    # local to one completion and retain the same reference/bounds checks.
    references = set(genome.references)
    used_references = {m['seqid'] for models in (selected, original, native) for m in models.values()}
    reference_lengths = {name: genome.get_reference_length(name) for name in used_references & references}
    original_by_gene = defaultdict(list)
    for m in original.values():
        if not m.get('problem'):
            original_by_gene[m['gene']].append(m)
    for model in native.values():
        if (model['seqid'] not in reference_lengths
                or any(b > reference_lengths[model['seqid']] for _, b in model['exons'])):
            raise ValueError('OrionGeno transcript lies outside target genome')
    donors, correspondence = native_correspondence(selected, native)
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
                p = project(donor, model, side, genome, reference_lengths=reference_lengths)
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
                          native_correspondence=correspondence[tid] if native_available else dict(status='missing_prediction', candidates=[]),
                          annotation_basis='projected_longest_per_side',
                          exon_path_basis=model.get('exon_path_basis', 'declared_or_cds_derived'),
                          source_exon_interpretations=[dict(transcript_id=m['id'], basis=m['exon_path_basis'],
                                                           declared_exons=m['declared_exons'])
                                                       for m in original_by_gene[model['gene']]
                                                       if m.get('exon_path_basis')]))
    return render(selected_text, updates), audit, sequences


def render(text, updates):
    """Update structural rows only; retain CDS IDs, coordinates, phases and attributes."""
    # Completion still emits a transcript sequence and audit for a CDS-only
    # model, but must not synthesize exons or label an unchanged source GFF.
    updates = {tid: value for tid, value in updates.items()
               if value['model'].get('utr') or any(c['blocks'] for c in value['choices'].values())
               or value['exons'] != value['model'].get('exons')}
    if not updates:
        return text
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
        if tid in updates and kind.endswith(('gene', 'rna', 'transcript')):
            value = updates[tid]
            if not kind.endswith('gene'):
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
