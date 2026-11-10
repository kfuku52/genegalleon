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

UTR_FEATURE_TYPES = {'5': 'five_prime_UTR', '3': 'three_prime_UTR'}


from cdskit.transcript_coordinates import TranscriptGroups as TranscriptGroups
from cdskit.transcript_coordinates import _coge_terminal_segments as _coge_terminal_segments
from cdskit.transcript_coordinates import _utr_low_side as _utr_low_side
from cdskit.transcript_coordinates import coding_boundary as coding_boundary
from cdskit.transcript_coordinates import exons_within_reference as exons_within_reference
from cdskit.transcript_coordinates import load_models as load_models
from cdskit.transcript_coordinates import merge_blocks as merge_blocks
from cdskit.transcript_coordinates import project as project
from cdskit.transcript_coordinates import read_model_rows as read_model_rows
from cdskit.transcript_coordinates import read_models as read_models
from cdskit.transcript_coordinates import terminal_utr_blocks as terminal_utr_blocks


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
    selected = read_models(selected_text, strict=False)
    original = selected if original_text == selected_text else read_models(original_text, strict=False)
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
        if not exons_within_reference(model, reference_lengths.get(model['seqid'])):
            raise ValueError('OrionGeno transcript lies outside target genome')
    donors, correspondence = native_correspondence(selected, native) if native_available else ({}, {})
    updates, audit, sequences = {}, [], {}
    for tid, model in selected.items():
        if model.get('problem'):
            audit.append(dict(species=species, transcript_id=tid, gene_id=model['gene'],
                              status='unsupported_structure', reason=model['problem']))
            continue
        candidates = [('existing', m) for m in original_by_gene[model['gene']]
                      if (m['seqid'], m['strand']) == (model['seqid'], model['strand'])]
        native_donors = donors.get(model['gene'], ())
        candidates += [('oriongeno', m) for m in native_donors]
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
                          native_status='mapped' if native_donors else 'no_unique_coding_overlap' if native_available else 'missing_prediction',
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
    groups = TranscriptGroups(updates)
    for value in updates.values():
        gene_bounds[value['model']['gene']].extend(value['exons'])
    for line in text.splitlines(keepends=True):
        if line.rstrip('\r\n') == '###':
            groups.close(len(lines))
        if not line.strip() or line.startswith('#'):
            lines.append(line)
            continue
        f = line.rstrip('\r\n').split('\t')
        attrs = parse_gff_attributes(f[8])
        groups.add(attrs)
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
    groups.close(None)
    used = {next(iter(parse_gff_attributes(line.split('\t')[8]).get('ID', ())), '')
            for line in lines if line.strip() and not line.startswith('#')}
    additions = defaultdict(list)
    for tid, value in updates.items():
        model = value['model']
        def emit(kind, block, suffix, extra='', tid=tid, model=model):
            identifier = tid+'.ggutr.'+suffix
            while identifier in used:
                identifier += '_'
            used.add(identifier)
            attrs = 'ID='+quote(identifier, safe='._-:')+';Parent='+quote(tid, safe='._-:')+extra
            additions[groups.placements.get(tid)].append(
                f"{model['seqid']}\tGeneGalleon\t{kind}\t{block[0]+1}\t{block[1]}\t.\t{model['strand']}\t.\t{attrs}\n")
        for i, block in enumerate(value['exons']):
            if (tid, 'exon', *block) not in retained_features:
                emit('exon', block, 'exon'+str(i+1))
        for side, choice in value['choices'].items():
            extra = ';annotation_basis=projected_longest;utr_source='+choice.get('source', 'none')+';utr_source_transcript='+quote(choice.get('source_transcript', ''), safe='._-:')
            for i, block in enumerate(choice['blocks']):
                if (tid, 'utr', *block) not in retained_features:
                    emit(UTR_FEATURE_TYPES[side], block, side+'.'+str(i+1), extra)
    def rendered_rows():
        last = '\n'
        for number, line in enumerate(lines):
            for addition in additions.get(number, ()):
                if not last.endswith('\n'):
                    yield '\n'
                yield addition
                last = addition
            yield line
            last = line
        for addition in additions.get(None, ()):
            if not last.endswith('\n'):
                yield '\n'
            yield addition
            last = addition
    return ''.join(rendered_rows()) + (marker + embedded if marker else '')
