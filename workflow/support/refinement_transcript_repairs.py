"""Attribute a coding-path repair to its source transcript(s) at one genomic locus."""
import hashlib
from collections import defaultdict
from itertools import pairwise

from Bio.Seq import Seq

try:
    import refinement_native_splits as native
    from rescue_coding_paths import coding_overlap, coding_shape
    from rescue_gene_models import validate_model
    from rescue_homology_policy import homology_supported
except ImportError:
    from . import refinement_native_splits as native
    from .rescue_coding_paths import coding_overlap, coding_shape
    from .rescue_gene_models import validate_model
    from .rescue_homology_policy import homology_supported


def junctions(blocks):
    return {(a[1], b[0]) for a, b in pairwise(sorted(blocks))}


def protected(model):
    return any(model['quality'].get(flag) for flag in (
        'annotated_exception', 'structure_problem', 'internal_stop', 'ambiguous',
        'invalid_base', 'translation_uncertain', 'phase_unresolved', 'sequence_mismatch'))


class TranscriptPaths:
    """Identical CDSs with different transcript/UTR identities share one comparison."""

    def __init__(self, models):
        self.groups = defaultdict(list)
        for m in models:
            self.groups[coding_shape({**m, 'cds': m['blocks']})].append(m)
        self.features = [(group, {**group[0], 'cds': group[0]['blocks']}, junctions(group[0]['blocks']))
                         for group in self.groups.values() if group[0]['blocks']]

    def match(self, path):
        if coding_shape(path) in self.groups:
            return None
        new = junctions(path['cds'])
        scores = []
        for models, source, old in self.features:
            _, frame = coding_overlap(source, path)
            if not frame:
                continue
            # Shared sequence alone is insufficient when mutually exclusive
            # internal exons match distinct transcripts equally well.
            a, b = source['cds'], path['cds']
            ends = frozenset(i for i, j in (
                ((0, 0), (-1, 1)) if path['strand'] == '+' else ((0, 1), (-1, 0))) if a[i][j] == b[i][j])
            scores.append(((old & new, frame, ends), models))
        winners = [models for score, models in scores if not any(
            other != score and all(a >= b for a,b in zip(other, score, strict=True))
            for other, _ in scores)]
        if len(winners) != 1:
            return None
        models = winners[0]
        if any(protected(m) or not m['quality'].get('usable') for m in models):
            return None
        return models


def nominate_projected(rows, paths, ownership, genome, code, params, name):
    results = {}
    for row in rows:
        if row['kind'] != 'coding_path':
            continue
        owner = row['owner_overlap'][0]['gene_id']
        support = [h for h in row['support'] if h['donor'] != name and
                   owner in (h.get('syntenic_locus'), h.get('reciprocal_locus'),
                             h.get('internal_gap_locus'), h.get('terminal_locus'))
                   and homology_supported(h['coverage'], h['identity'], params)]
        if len({h['donor'] for h in support}) < max(2, params['min_support']):
            continue
        blocks = row['cds']
        path = {**row, 'start': min(b[0] for b in blocks), 'end': max(b[1] for b in blocks)}
        if path['end'] - path['start'] > params['max_interval']:
            continue
        models = paths.match(path)
        if not models:
            continue
        owners = ownership.overlapping(path)
        if len(owners) != 1 or owners[0]['gene_id'] != owner or owners[0]['strand'] != row['strand']:
            continue
        original = models[0]
        _, frame = coding_overlap({**original, 'cds': original['blocks']}, path)
        if frame < (.90 if native.complete_source_orf(original) else 2/3) * len(original['cds']):
            continue
        checked = validate_model({**path, 'coverage': 1., 'identity': 1., 'frameshift': False},
                                 genome, code, params)
        if (checked['problems'] or checked['cds'] != blocks or
                hashlib.sha256(checked['sequence'].encode()).hexdigest() != row['cds_sha256']):
            continue
        protein = str(Seq(checked['sequence']).translate(table=code)).removesuffix('*')
        if protein == original['protein'] and junctions(blocks) == junctions(original['blocks']):
            continue
        results[coding_shape(path)] = dict(models=models, checked=checked, protein=protein,
                                           review={**row, 'support': support, 'adoption': 'automatic'})
    return results


def nominate_native(paths, native_paths, owner, genome, code, params, name, rna_available, projected):
    results = {}
    for path in native_paths:
        shape = coding_shape(path)
        if shape in projected:
            projected[shape]['review'].setdefault('oriongeno_transcripts', []).append(path['native_transcript_id'])
            continue
        if shape in results:
            results[shape]['review']['oriongeno_transcripts'].append(path['native_transcript_id'])
            continue
        models = paths.match(path)
        if not models:
            continue
        candidate = native.nominate_extensions({owner: [path]}, {owner: models[0]}, genome, code, params,
                                               name, rna_available=rna_available).get(shape)
        if candidate:
            results[shape] = {**candidate, 'models': models}
    return results


def changed_blocks(old, new):
    """Keep the source-to-revision change compact in the existing audit."""
    a, b = set(map(tuple, old)), set(map(tuple, new))
    return {'removed': sorted(a-b), 'added': sorted(b-a)}


def register(by_transcript, model):
    """Key the publication by transcript ID within its species and frozen annotation."""
    tid = model['source_transcript_id']
    if tid in by_transcript:
        raise ValueError('Two repairs would publish the same transcript ID: ' + tid)
    by_transcript[tid] = model
