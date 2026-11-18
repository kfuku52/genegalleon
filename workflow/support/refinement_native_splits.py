"""Use native coding paths and cross-species evidence to resolve locus boundaries."""
import csv
import hashlib
import re
from collections import defaultdict
from itertools import combinations, pairwise
from pathlib import Path

from Bio import Phylo
from Bio.Seq import Seq

try:
    from fasta_sequence_store import fasta_records
    from gene_model_store import iter_loci
    from rescue_coding_paths import OwnershipIndex, coding_overlap, coding_shape
    from rescue_gene_models import run, validate_model
    from rescue_homology_policy import homology_supported
except ImportError:
    from .fasta_sequence_store import fasta_records
    from .gene_model_store import iter_loci
    from .rescue_coding_paths import OwnershipIndex, coding_overlap, coding_shape
    from .rescue_gene_models import run, validate_model
    from .rescue_homology_policy import homology_supported


def complete_source_orf(model):
    quality = model['quality']
    return quality['valid_orf'] or (quality.get('usable') and quality.get('has_start')
            and quality.get('has_stop') and not quality.get('translation_offset')
            and not quality.get('incomplete_codon'))


def _in_frame_codons(fragment, path, protein_length):
    """Project complete source codons onto the native translated path."""
    reverse = path['strand'] == '-'
    boundaries = defaultdict(int)
    codons = set()
    offset = -path.get('quality', {}).get('translation_offset', 0)
    for start, end, phase in path['cds']:
        for a, b, source_phase in fragment['blocks']:
            left, right = max(start, a), min(end, b)
            if left >= right:
                continue
            if (end - b - phase + source_phase if reverse else a - start - phase + source_phase) % 3:
                continue
            begin = offset + (end - right if reverse else left - start)
            limit = min(3 * protein_length, begin + right - left)
            begin = max(0, begin)
            first, stop = (begin + 2) // 3, limit // 3
            codons.update(range(first, stop))
            for position in range(begin, min(first * 3, limit)):
                boundaries[position // 3] |= 1 << (position % 3)
            for position in range(max(first * 3, stop * 3, begin), limit):
                boundaries[position // 3] |= 1 << (position % 3)
        offset += end - start
    codons.update(i for i, mask in boundaries.items() if mask == 7)
    return codons


def nominate(paths, originals, ownership, genome, code, params, species, excluded_shapes, *, rna_available=False, intact_only=False, comparative=False, projected_shapes=()):
    """Validate native CDS and unambiguous source geometry before protein search."""
    candidates = {}
    for path in paths:
        shape = coding_shape(path)
        if shape in excluded_shapes or path['end'] - path['start'] > params['max_interval']:
            continue
        existing = ownership.overlapping(path)
        if len(existing) < 2 or any(owner['strand'] != path['strand']
                                    or owner['gene_id'] not in originals for owner in existing):
            continue
        fragments = [originals[owner['gene_id']] for owner in existing]
        intact = [fragment for fragment in fragments if complete_source_orf(fragment)]
        needs_rna = len(intact) > 1
        confirmed = comparative and shape in projected_shapes
        if needs_rna and not (rna_available or confirmed):
            continue
        bounds = sorted((min(a for a, _, _ in fragment['blocks']), max(b for _, b, _ in fragment['blocks']))
                        for fragment in fragments)
        members = {fragment['source_gene_id'] for fragment in fragments}
        for (_, left), (right, _) in pairwise(bounds):
            if left > right:
                break
            if left < right and any(owner['gene_id'] not in members for owner in
                                    ownership.overlapping({**path, 'start': left, 'end': right, 'cds': [[left, right]]})):
                break
        else:
            geometry = []
            low_retention = False
            for fragment in fragments:
                overlap, frame = coding_overlap({**fragment, 'cds': fragment['blocks']}, path)
                size = len(fragment['cds'])
                low_retention |= frame < .80 * size
                if frame < size * (2 / 3 if confirmed or not needs_rna else .80):
                    break
                geometry.append({'gene_id': fragment['source_gene_id'], 'cds_bases': size, 'overlap_bases': overlap})
            else:
                assess = comparative and (confirmed and needs_rna or low_retention and not needs_rna)
                if (intact_only and not (needs_rna or assess)) or (low_retention and not comparative):
                    continue
                if shape in candidates:
                    if 'native_transcript_id' in path:
                        candidates[shape]['review']['oriongeno_transcripts'].append(path['native_transcript_id'])
                    continue
                checked = validate_model({**path, 'coverage': 1., 'identity': 1., 'frameshift': False},
                                         genome, code, params)
                if checked['problems'] or coding_shape(checked) != shape:
                    continue
                protein = str(Seq(checked['sequence']).translate(table=code)).removesuffix('*')
                windows = {fragment['source_gene_id']: _in_frame_codons(fragment, path, len(protein))
                           for fragment in fragments}
                if not all(windows.values()):
                    continue
                key = (shape, tuple(sorted(members)))
                identifier = species + '_ggmerge_' + hashlib.sha256(repr(key).encode()).hexdigest()[:16]
                candidates[shape] = {'review': {
                    'candidate_id': identifier, 'seqid': path['seqid'], 'strand': path['strand'], 'cds': checked['cds'],
                    'cds_sha256': hashlib.sha256(checked['sequence'].encode()).hexdigest(),
                    'owner_overlap': sorted(geometry, key=lambda row: row['gene_id']),
                    'oriongeno_transcripts': ([path['native_transcript_id']] if 'native_transcript_id' in path else []),
                    'adoption': 'automatic',
                    'basis': ('oriongeno_cross_fragment_homology_and_target_rna' if needs_rna else
                              'oriongeno_and_cross_fragment_protein_homology')},
                    'checked': checked, 'protein': protein, 'windows': windows,
                    'protect_rna_paths': intact if needs_rna else []}
                if assess:
                    candidate = candidates[shape]
                    candidate['source_fragments'] = [
                        (f['source_gene_id'], f['protein'], _in_frame_codons(
                            {'blocks': path['cds']}, {**f, 'cds': f['blocks']}, len(f['protein'])))
                        for f in sorted(fragments, key=lambda f: min(windows[f['source_gene_id']]))]
                    candidate['protect_rna_paths'] = fragments
                    candidate['review']['basis'] = 'conserved_coding_bridge_and_complementary_fragments'
    return candidates


def nominate_extensions(groups, originals, genome, code, params, species, *, rna_available=False):
    """Recover missing terminal CDS with in-frame overlap and donor coverage of the added residues."""
    candidates = {}
    for owner, paths in groups.items():
        original = originals.get(owner)
        if not original:
            continue
        quality, source_cds = original['quality'], original['blocks']
        intact = complete_source_orf(original)
        if (intact and not rna_available) or (not intact and not quality['partial']):
            continue
        old_size = len(original['cds'])
        source = {**original, 'cds': source_cds}
        for path in paths:
            size = sum(b - a for a, b, _ in path['cds'])
            if size - old_size < max(90, .20 * old_size) or path['end'] - path['start'] > params['max_interval']:
                continue
            overlap, frame = coding_overlap(source, path)
            if frame < .90 * old_size:
                continue
            if not intact:
                start, stop = (0, 1) if source['strand'] == '+' else (1, 0)
                if ((quality['has_start'] and source_cds[0][start] != path['cds'][0][start])
                        or (quality['has_stop'] and source_cds[-1][stop] != path['cds'][-1][stop])):
                    continue
                if not (quality['has_start'] or quality['has_stop']) and frame < .95 * old_size:
                    continue
            shape = coding_shape(path)
            if shape in candidates:
                candidates[shape]['review']['oriongeno_transcripts'].append(path['native_transcript_id'])
                continue
            checked = validate_model({**path, 'coverage': 1., 'identity': 1., 'frameshift': False}, genome, code, params)
            if checked['problems'] or coding_shape(checked) != shape:
                continue
            protein = str(Seq(checked['sequence']).translate(table=code)).removesuffix('*')
            retained = _in_frame_codons(original, path, len(protein))
            added = set(range(len(protein))) - retained
            if not retained or not added:
                continue
            identifier = species + '_ggrepair_' + hashlib.sha256(repr((owner, shape)).encode()).hexdigest()[:16]
            candidates[shape] = {
                'review': {'candidate_id': identifier, 'kind': 'coding_path', 'seqid': path['seqid'],
                           'strand': path['strand'], 'cds': checked['cds'], 'adoption': 'automatic',
                           'cds_sha256': hashlib.sha256(checked['sequence'].encode()).hexdigest(),
                           'owner_overlap': [{'gene_id': owner, 'cds_bases': old_size, 'overlap_bases': overlap}],
                           'oriongeno_transcripts': [path['native_transcript_id']],
                           'basis': ('oriongeno_cds_extension_with_target_rna' if intact else
                                     'oriongeno_partial_cds_extension')},
                'checked': checked, 'protein': protein, 'windows': {'retained': retained, 'added': added},
                'protect_rna_paths': [original] if intact else []}
    return candidates


def collinear_gap(first, second):
    if not first or not second:
        return None
    a_id, a_seq, a_strand, a_start, a_end = first
    b_id, b_seq, b_strand, b_start, b_end = second
    if a_id == b_id or a_seq != b_seq or a_strand != b_strand or a_strand not in {'+', '-'}:
        return None
    gap = b_start - a_end if a_strand == '+' else a_start - b_end
    return gap if gap >= 0 else None


def native_paths_by_ownership(paths, ownership, *, include_merges, extensions=None, single_paths=None):
    merge_paths = []
    singles = defaultdict(list)
    source_sizes = {}
    for path in paths:
        owners = ownership.overlapping(path)
        if include_merges and len(owners) >= 2:
            merge_paths.append(path)
        elif len(owners) == 1 and owners[0]['strand'] == path['strand']:
            member = owners[0]
            key = member['gene_id']
            singles[key].append(path)
            if extensions is not None:
                if key not in source_sizes:
                    source_sizes[key] = sum(b - a for a, b in set(map(tuple, member['cds'])))
                old = source_sizes[key]
                if old and sum(b - a for a, b, _ in path['cds']) - old >= max(90, .20 * old):
                    extensions[key].append(path)
    if single_paths is not None:
        single_paths.update(singles)
    return merge_paths, {owner: paths for owner, paths in singles.items()
                         if len({path['native_gene_id'] for path in paths}) >= 2}


def fission_overlaps(candidate, original, blocked):
    blocks = original['blocks']
    source = {**original, 'cds': blocks, 'start': min(b[0] for b in blocks), 'end': max(b[1] for b in blocks)}
    return blocked.overlapping(source) or any(blocked.overlapping(path) for path in candidate['paths'])


def nominate_fissions(groups, originals, genome, code, params, blocked):
    """Find two complete ORFs that reconstruct a source CDS across a splice."""
    candidates = {}
    for owner, paths in groups.items():
        original = originals.get(owner)
        if not original or not original['quality']['valid_orf'] or len(original['blocks']) < 2:
            continue
        blocks = original['blocks']
        source = {**original, 'cds': blocks, 'start': min(b[0] for b in blocks), 'end': max(b[1] for b in blocks)}
        if blocked.overlapping(source):
            continue
        source_size = len(original['cds'])
        pieces = []
        for path in paths:
            size = sum(b - a for a, b, _ in path['cds'])
            if (size < 303 or path['end'] - path['start'] > params['max_interval']
                    or blocked.overlapping(path)):
                continue
            _, in_frame = coding_overlap(source, path)
            if in_frame < max(.80 * size, .10 * source_size):
                continue
            checked = validate_model({**path, 'coverage': 1., 'identity': 1., 'frameshift': False},
                                     genome, code, params)
            if checked['problems'] or coding_shape(checked) != coding_shape(path):
                continue
            protein = str(Seq(checked['sequence']).translate(table=code)).removesuffix('*')
            pieces.append((checked, in_frame, protein))
        reverse = source['strand'] == '-'
        pieces.sort(key=lambda piece: piece[0]['start'], reverse=reverse)
        for (first, a_size, a_protein), (second, b_size, b_protein) in combinations(pieces, 2):
            if (first['native_gene_id'] == second['native_gene_id']
                    or (second['end'] > first['start'] if reverse else first['end'] > second['start'])
                    or a_size + b_size < .85 * source_size):
                continue
            bridges = []
            for left, right in pairwise(source['cds']):
                if not (coding_overlap({**source, 'cds': [left]}, first)[1]
                        and coding_overlap({**source, 'cds': [right]}, second)[1]):
                    continue
                bridges.append([right[1], left[0]] if reverse else [left[1], right[0]])
            if len(bridges) != 1:
                continue
            windows = tuple(_in_frame_codons({'blocks': piece['cds']}, source, len(original['protein']))
                            for piece in (first, second))
            if not all(windows):
                continue
            key = owner, coding_shape(first), coding_shape(second)
            candidates.setdefault(key, {
                'owner': owner, 'source_protein': original['protein'], 'source_windows': windows,
                'proteins': (a_protein, b_protein), 'paths': (first, second),
                'bridge': bridges[0], 'source_in_frame_bases': (a_size, b_size),
                'source_span': {k: source[k] for k in ('seqid', 'start', 'end', 'strand')},
            })
    return candidates


def sister_clades(path, cohorts, *, midpoint_root=False):
    """Group donors by where their lineage meets the target, nearest first."""
    tree = Phylo.read(path, 'newick')
    leaves = tree.get_terminals()
    names = {leaf.name for leaf in leaves}
    expected = set(cohorts) | {donor for cohort in cohorts.values() for donor in cohort}
    if len(names) != len(leaves) or not expected <= names:
        raise ValueError('Refinement species tree must have distinct leaves and include all native donor cohorts')
    if midpoint_root:
        tree.root_at_midpoint()
    groups = {}
    for target, donors in cohorts.items():
        lineage = [tree.root, *tree.get_path(target)]
        allowed = set(donors)
        groups[target] = [
            [members for sibling in parent.clades if sibling is not child
             if (members := sorted(allowed.intersection(leaf.name for leaf in sibling.get_terminals())))]
            for parent, child in reversed(list(pairwise(lineage)))
        ]
    return groups


def fission_decision(sister_groups, split, fused, minimum_species):
    """Require support on distinct sister branches and within the nearest informative group."""
    informative = set(split) | set(fused)
    split_groups = [sorted(set(group) & split.keys()) for level in sister_groups for group in level
                    if set(group) & split.keys()]
    closest = next((sorted({name for group in level for name in group} & informative)
                    for level in sister_groups if any(informative.intersection(group) for group in level)), [])
    if len(split) < minimum_species:
        reason = 'insufficient_collinear_donor_species'
    elif not sister_groups:
        reason = 'species_tree_unavailable'
    elif len(split_groups) < 2:
        reason = 'split_support_on_one_sister_branch'
    elif any(fused[name].get('rna_witnesses', 0) >= 2 for name in closest if name in fused):
        reason = 'closest_donor_fusion_supported_by_rna'
    elif any(hit.get('syntenic_rna_witnesses', 0) >= 2 for hit in fused.values()):
        reason = 'syntenic_donor_fusion_supported_by_rna'
    elif not set(closest) & split.keys():
        reason = 'closest_informative_donors_support_fusion'
    elif any(name in fused for name in closest):
        reason = 'mixed_closest_donor_architectures'
    elif len(split) == 2 and any(hit.get('within_syntenic_flanks') for hit in fused.values()):
        reason = 'syntenic_fusion_with_minimal_split_support'
    else:
        reason = 'supported_split'
    return {'reason': reason, 'closest_informative_species': closest, 'split_support_groups': split_groups}


def _paired_runs(traceback, start):
    """Yield query spans paired to donor residues in a DIAMOND BTOP."""
    for token in re.findall(r'\d+|[^\d]{2}', traceback):
        if token.isdecimal():
            count = int(token)
            yield start, count
            start += count
        elif token[0] != '-':
            if token[1] != '-':
                yield start, 1
            start += 1


def _paired_positions(traceback, start):
    paired = set()
    for offset, count in _paired_runs(traceback, start):
        paired.update(range(offset, offset + count))
    return paired


def flanking_anchor_index(blocks, target_positions, donor_positions, max_interval):
    """Index bounded inter-anchor intervals in both genomes, including inversions."""
    def interval(left, right):
        if left is None or right is None or left[1] != right[1]:
            return None
        a, b = sorted((left, right), key=lambda x: x[3])
        return (a[1], a[4], b[3]) if 0 < b[3] - a[4] <= max_interval else None

    bounds = set()
    for block in blocks:
        for (ta, da), (tb, db) in pairwise(block):
            target = interval(target_positions.get(ta), target_positions.get(tb))
            donor = interval(donor_positions.get(da), donor_positions.get(db))
            if target and donor:
                bounds.add((target, donor))
    return OwnershipIndex({'seqid': seqid, 'start': start, 'end': end,
                           'gene_id': ((seqid, start, end), donor), 'donor': donor}
                          for (seqid, start, end), donor in bounds)


def uniquely_flanked_locus(index, span):
    windows = {row['donor'] for row in index.overlapping(span)
               if row['start'] <= span['start'] and span['end'] <= row['end']}
    return next(iter(windows)) if len(windows) == 1 else None


def outer_anchor_blocks(blocks, target, donor):
    """Index collinear blocks even when a candidate itself contains an anchor."""
    features = []
    for block in blocks:
        groups = defaultdict(list)
        for t, d in block:
            if t in target and d in donor:
                a, b = target[t], donor[d]
                groups[a[1], b[1]].append((a[3], a[4], b[3], b[4]))
        for (chrom, other), anchors in groups.items():
            if len(anchors) > 1:
                features.append({'seqid': chrom, 'start': min(x[0] for x in anchors),
                                 'end': max(x[1] for x in anchors), 'gene_id': len(features)+1,
                                 'other': other, 'anchors': anchors})
    return OwnershipIndex(features)


def outer_regions(index, span, limit):
    regions = set()
    for block in index.overlapping(span):
        anchors = block['anchors']
        left = [a for a in anchors if 0 <= span['start']-a[1] <= limit]
        right = [a for a in anchors if 0 <= a[0]-span['end'] <= limit]
        if not left or not right:
            continue
        lo, hi = max(a[1] for a in left), min(a[0] for a in right)
        for a in left:
            if a[1] != lo:
                continue
            for b in right:
                if b[0] != hi:
                    continue
                if a[3] <= b[2]:
                    start, end, strand = a[3], b[2], span['strand']
                elif b[3] <= a[2]:
                    start, end, strand = b[3], a[2], '-' if span['strand'] == '+' else '+'
                else:
                    continue
                if 0 < end-start <= 2*limit:
                    regions.add((block['other'], start, end, strand))
    return regions


def merge_query_proteins(candidates):
    for candidate in candidates.values():
        yield candidate['protein']
        for _, protein, _ in candidate.get('source_fragments', ()):
            yield protein


def bridges_paired(windows, paired):
    spans = sorted((min(window), max(window)) for window in windows.values())
    for (_, left), (right, _) in pairwise(spans):
        if left >= right:
            return False
        bridge = range(max(0, left - 2), right + 3)
        if left not in paired or right not in paired or sum(i in paired for i in bridge) < .8 * len(bridge):
            return False
    return True


class DonorSearch:
    """Share one lazily built protein database across the finalization targets."""

    def __init__(self, panel, db, directory, cpus, rna_counts=None, *, flanking_blocks=None):
        self.panel = panel
        self.db = db
        self.directory = Path(directory)
        self.cpus = cpus
        self.subjects = []
        self.positions = {}
        self.rna_counts = defaultdict(dict)
        for (species, seqid, strand, blocks), count in (rna_counts or {}).items():
            if count >= 2:
                self.rna_counts[species][seqid, strand, blocks] = count
        self.observed_fusions = {}
        self.database = self.directory / 'donors.dmnd'
        self._queries = ()
        self.flanking_blocks = flanking_blocks
        self.outer_target = None
        self.outer_blocks = {}
        self.ranks = {}

    def _proteins(self, name):
        path = self.panel.get('proteins', {}).get(name)
        if path:
            for identifier, _, protein in fasta_records(Path(path)):
                if protein:
                    yield identifier, protein
        else:
            for locus in iter_loci(self.db, name):
                for model in locus['candidates']:
                    if model['origin'] == 'original' and model['quality']['valid_orf']:
                        yield model['candidate_id'], model['protein']

    def _build(self):
        self.directory.mkdir()
        fasta = self.directory / 'donors.fa'
        with fasta.open('w') as handle:
            for donor in sorted({donor for cohort in self.panel['donors'].values() for donor in cohort}):
                for identifier, protein in self._proteins(donor):
                    handle.write(f'>d{len(self.subjects)}\n{protein}\n')
                    self.subjects.append((donor, identifier))
        if self.subjects:
            run([self.panel['executable']['path'], 'makedb', '--in', fasta, '--db', self.database,
                 '--threads', self.cpus], self.directory, 'makedb')
        fasta.unlink()

    def search(self, proteins):
        """Write one temporary DIAMOND result and return its query index."""
        self._queries = ()
        sequences = tuple(dict.fromkeys(proteins))
        if not sequences:
            return ()
        if not self.directory.exists():
            self._build()
        if not self.subjects:
            return ()
        query = self.directory / 'queries.fa'
        output = self.directory / 'hits.tsv'
        with query.open('w') as handle:
            for index, protein in enumerate(sequences):
                handle.write(f'>q{index}\n{protein}\n')
        run([self.panel['executable']['path'], 'blastp', '--query', query, '--db', self.database,
             '--threads', self.cpus, '--very-sensitive', '--evalue', '1e-5', '--max-target-seqs', '0',
             '--max-hsps', '1', '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'qstart', 'qend', 'qlen',
             'sstart', 'send', 'slen', 'bitscore', 'btop'], self.directory, 'blastp', output)
        self._queries = sequences
        return sequences

    def alignments(self, species, proteins, *, batch=None, subject_intervals=False):
        if not proteins:
            return
        sequences = self.search(proteins) if batch is None else batch
        if not sequences:
            return
        if sequences is not self._queries:
            raise ValueError('DIAMOND result was replaced by another search')
        query_ids = {index for index, protein in enumerate(sequences) if protein in proteins}
        output = self.directory / 'hits.tsv'
        allowed = set(self.panel['donors'][species])
        with output.open() as handle:
            for line in handle:
                q, s, identity, qstart, qend, qlen, sstart, send, slen, bitscore, traceback = line.split()
                index = int(q[1:])
                if index not in query_ids:
                    continue
                donor, identifier = self.subjects[int(s[1:])]
                if donor not in allowed:
                    continue
                qstart, qend, qlen, sstart, send, slen = map(int, (qstart, qend, qlen, sstart, send, slen))
                query_coverage = (qend - qstart + 1) / qlen
                donor_coverage = (send - sstart + 1) / slen
                yield sequences[index], {
                    'donor': donor, 'query': identifier, 'identity': float(identity) / 100,
                    'coverage': min(query_coverage, donor_coverage), 'query_coverage': query_coverage,
                    'donor_coverage': donor_coverage, 'bitscore': float(bitscore),
                    **({'subject_start': sstart - 1, 'subject_end': send} if subject_intervals else {}),
                }, traceback, qstart - 1, max(qlen, slen)

    def support(self, species, candidates, params, *, batch=None):
        guarded = {shape: c for shape, c in candidates.items() if c.get('source_fragments')}
        if guarded and batch is None:
            batch = self.search(merge_query_proteins(candidates))

        queries = defaultdict(list)
        for shape, candidate in candidates.items():
            queries[candidate['protein']].append((shape, candidate['windows']))
        best = defaultdict(dict)
        for protein, hit, traceback, start, length in self.alignments(species, queries, batch=batch):
            coverage, identity = hit['coverage'], hit['identity']
            if identity < .50 or coverage < .90 or not homology_supported(coverage, identity, params):
                continue
            paired = _paired_positions(traceback, start)
            paired_coverage = len(paired) / length
            if paired_coverage < params['minimum_coverage']:
                continue
            donor = hit['donor']
            for shape, windows in queries[protein]:
                fragment_coverage = {gene: len(window & paired) / len(window) for gene, window in windows.items()}
                if min(fragment_coverage.values()) < .80:
                    continue
                if shape in guarded and not bridges_paired(windows, paired):
                    continue
                previous = best[shape].get(donor)
                if previous is None or (hit['bitscore'], hit['query']) > (previous['bitscore'], previous['query']):
                    best[shape][donor] = {**hit, 'paired_coverage': paired_coverage,
                                         'fragment_coverage': fragment_coverage, 'flags': []}
        if guarded:
            self._check_source_architecture(species, guarded, best, batch)
        return {shape: [hits[donor] for donor in sorted(hits)] for shape, hits in best.items()
                if len(hits) >= max(2, params['min_support'])}

    def _check_source_architecture(self, species, candidates, best, batch):
        fragments = defaultdict(list)
        for shape, candidate in candidates.items():
            if best.get(shape):
                for gene, protein, retained in candidate['source_fragments']:
                    fragments[protein].append((shape, gene, retained))
        aligned, separate = {}, defaultdict(dict)
        for protein, hit, traceback, start, _ in self.alignments(species, fragments, batch=batch, subject_intervals=True):
            if hit['identity'] < .50 or hit['query_coverage'] < 2 / 3:
                continue
            donor, query = hit['donor'], hit['query']
            paired = _paired_positions(traceback, start)
            for shape, gene, retained in fragments[protein]:
                if hit['coverage'] >= .85 and len(paired) >= .80 * len(protein):
                    rank = separate[shape, donor]
                    if gene not in rank or hit['bitscore'] > rank[gene][0]:
                        rank[gene] = (hit['bitscore'], {query})
                    elif hit['bitscore'] == rank[gene][0]:
                        rank[gene][1].add(query)
                if query != best[shape].get(donor, {}).get('query'):
                    continue
                kept = len(paired & retained)
                if kept >= .90 * len(paired) and kept >= 15:
                    aligned[shape, donor, gene] = hit['subject_start'], hit['subject_end']
        for shape, candidate in candidates.items():
            genes = [gene for gene, _, _ in candidate['source_fragments']]
            for donor in list(best[shape]):
                spans = [aligned.get((shape, donor, gene)) for gene in genes]
                if not all(spans) or any(left[1] > right[0] for left, right in pairwise(spans)):
                    del best[shape][donor]
            for donor in self.panel['donors'][species]:
                ranks = separate.get((shape, donor), {})
                if not all(gene in ranks for gene in genes):
                    continue
                positions = self._positions(donor)
                # Tied paralogs cannot establish a separate, ordered architecture.
                if any(len(ranks[g][1]) != 1 for g in genes):
                    continue
                loci = [positions.get(next(iter(ranks[g][1]))) for g in genes]
                if all((gap := collinear_gap(a, b)) is not None and gap <= 10_000
                       for a, b in pairwise(loci)):
                    # An adjacent paralogue pair outside every flanked interval is
                    # informative only when the fused orthologue itself is localized.
                    projections = candidate.get('projected_queries', {})
                    span = candidate.get('checked')
                    if projections and span and self.relative_locus(
                            species, donor, span, [next(iter(ranks[g][1])) for g in genes]) == 'outside':
                        corroborated = [d for d, h in best[shape].items() if
                                        h['query'] in projections.get(d, ()) and self.relative_locus(
                                            species, d, span, [h['query']]) == 'inside']
                        if len(corroborated) >= 2 and sum(
                                h['query'] in projections.get(d, ()) for d, h in best[shape].items()) >= 3:
                            for d in corroborated:
                                best[shape][d]['flags'].append('separate_paralogs_outside_flanks:' + donor)
                            continue
                    best[shape].clear()
                    break

    def _positions(self, donor):
        if donor not in self.positions:
            bed = self.panel.get('positions', {}).get(donor)
            positions = {}
            if bed:
                with Path(bed).open() as handle:
                    positions = read_positions(handle)
            else:
                for locus in iter_loci(self.db, donor):
                    for model in locus['candidates']:
                        if model['origin'] == 'original' and model['quality']['valid_orf'] and model['blocks']:
                            blocks = model['blocks']
                            positions[model['candidate_id']] = (
                                model['source_gene_id'], model['seqid'], model['strand'],
                                min(b[0] for b in blocks), max(b[1] for b in blocks))
            self.positions[donor] = positions
        return self.positions[donor]

    def relative_locus(self, species, donor, span, queries):
        """Locate homologues between two outer anchors within 500 kb of each target edge."""
        if not self.flanking_blocks:
            return None
        if self.outer_target != species:
            self.outer_blocks.clear()
            self.outer_target = species
        positions = self._positions(donor)
        if donor not in self.outer_blocks:
            pair = self.flanking_blocks(species, donor)
            self.outer_blocks[donor] = outer_anchor_blocks(*pair, positions) if pair else None
        blocks = self.outer_blocks[donor]
        regions = outer_regions(blocks, span, 500_000) if blocks else ()
        if not regions or not queries or any(q not in positions for q in queries):
            return None
        if all(all(within_flanks(r[:3], positions[q]) and positions[q][2] == r[3] for q in queries)
               for r in regions):
            return 'inside'
        if all(outside_flanks(r[:3], positions, queries) for r in regions):
            return 'outside'
        return None

    def _consecutive(self, donor, a, b):
        if donor not in self.ranks:
            chromosomes = defaultdict(list)
            for key, (_, chrom, _, start, end) in self._positions(donor).items():
                chromosomes[chrom].append((start, end, key))
            self.ranks[donor] = {key: i for loci in chromosomes.values()
                                 for i, (_, _, key) in enumerate(sorted(loci))}
        return abs(self.ranks[donor][a] - self.ranks[donor][b]) == 1

    def _rna_supported_donor_models(self, donor):
        if donor not in self.observed_fusions:
            counts = self.rna_counts.get(donor, {})
            direct, supported_loci = {}, defaultdict(dict)
            if counts:
                for locus in iter_loci(self.db, donor):
                    for candidate in locus['candidates']:
                        if candidate['origin'] != 'original' or not candidate['quality']['valid_orf']:
                            continue
                        blocks = tuple(tuple(b[:2]) for b in candidate['blocks'])
                        count = counts.get((candidate['seqid'], candidate['strand'], blocks))
                        if count and len(blocks) > 1:
                            direct[candidate['candidate_id']] = count
                            proteins = supported_loci[candidate['source_gene_id']]
                            protein = candidate['protein']
                            proteins[protein] = max(proteins.get(protein, 0), count)
            mapping = self.panel.get('id_maps', {}).get(donor)
            if mapping and supported_loci:
                with Path(mapping).open(newline='') as handle:
                    selected = {row['jcvi_id']: row['locus_id'] for row in csv.DictReader(handle, delimiter='\t')
                                if row['status'] == 'selected' and row['locus_id'] in supported_loci}
                direct = {identifier: supported_loci[selected[identifier]][protein]
                          for identifier, protein in self._proteins(donor)
                          if identifier in selected and protein in supported_loci[selected[identifier]]}
            elif mapping:
                direct = {}
            self.observed_fusions[donor] = direct
        return self.observed_fusions[donor]

    def fission_support(self, species, candidates, params, *, batch=None, accession_keys=()):
        proteins = {protein for candidate in candidates.values()
                    for protein in (candidate['source_protein'], *candidate['proteins'])}
        source_windows = defaultdict(dict)
        for key, candidate in candidates.items():
            source_windows[candidate['source_protein']][key] = [
                (sum(1 << i for i in window), .80 * len(window)) for window in candidate['source_windows']]
        best = defaultdict(dict)
        spanning = defaultdict(dict)
        spanning_loci = defaultdict(lambda: defaultdict(set))
        spanning_rna_donors = defaultdict(lambda: defaultdict(set))
        for protein, hit, traceback, start, length in self.alignments(species, proteins, batch=batch):
            if hit['identity'] < .45 or hit['coverage'] < .65:
                continue
            windows_by_candidate = source_windows.get(protein) if hit['identity'] >= .50 else None
            count = mask = 0
            for offset, size in _paired_runs(traceback, start):
                count += size
                if windows_by_candidate:
                    mask |= ((1 << size) - 1) << offset
            proportion = count / length
            if proportion < .60:
                continue
            donor = hit['donor']
            if windows_by_candidate and proportion >= .80:
                for key, windows in windows_by_candidate.items():
                    if not all((mask & window).bit_count() >= minimum for window, minimum in windows):
                        continue
                    if self.flanking_blocks:
                        spanning_loci[key][donor].add(hit['query'])
                    if donor in self.rna_counts:
                        spanning_rna_donors[key][donor].add(hit['query'])
                    previous = spanning[key].get(donor)
                    if previous is None or (hit['bitscore'], hit['query']) > (previous['bitscore'], previous['query']):
                        spanning[key][donor] = {**hit, 'paired_coverage': proportion}
            previous = best[protein].get(donor)
            if previous is None or hit['bitscore'] > previous['bitscore']:
                best[protein][donor] = {**hit, 'paired_coverage': proportion, 'tied': False}
            elif hit['bitscore'] == previous['bitscore'] and hit['query'] != previous['query']:
                previous['tied'] = True

        result = {}
        flank_indices = {}
        groups = self.panel.get('sister_clades', {}).get(species, [])
        minimum = max(2, params['min_support'])
        for key, candidate in candidates.items():
            opposed = spanning[key]
            first, second = (best[protein] for protein in candidate['proteins'])
            support, extended = {}, {}
            for donor in sorted(first.keys() & second.keys()):
                a, b = first[donor], second[donor]
                if a['tied'] or b['tied']:
                    continue
                positions = self._positions(donor)
                gap = collinear_gap(positions.get(a['query']), positions.get(b['query']))
                if gap is None or gap > 100_000:
                    continue
                bucket = support
                if gap > 10_000:
                    if not self.flanking_blocks or not self._consecutive(donor, a['query'], b['query']):
                        continue
                    bucket = extended
                bucket[donor] = {'species': donor, 'gap_bases': gap,
                                      'first': {k: v for k, v in a.items() if k != 'tied'},
                                      'second': {k: v for k, v in b.items() if k != 'tied'}}
            # Treat donor species with both architectures as ambiguous at this locus.
            unambiguous = {name: hit for name, hit in support.items() if name not in opposed}
            decision = fission_decision(groups, unambiguous, opposed, minimum)
            if (extended and not opposed and len(support) + len(extended) >= max(3, minimum)
                    and decision['reason'] in ('insufficient_collinear_donor_species', 'split_support_on_one_sister_branch')):
                combined = {**support, **extended}
                revised = fission_decision(groups, combined, opposed, minimum)
                if revised['reason'] == 'supported_split':
                    localized_pairs = []
                    for d, h in combined.items():
                        a, b = h['first']['query'], h['second']['query']
                        if self._consecutive(d, a, b) and self.relative_locus(
                                species, d, candidate['source_span'], [a, b]) == 'inside':
                            localized_pairs.append(d)
                    if len(localized_pairs) >= 2:
                        support = unambiguous = combined
                        decision = {**revised, 'syntenic_adjacent_pairs': sorted(localized_pairs)}
            excluded = []
            localized = {}
            if self.flanking_blocks and (len(support) >= minimum or key in accession_keys) and opposed:
                for donor in sorted(opposed):
                    if donor not in flank_indices:
                        inputs = self.flanking_blocks(species, donor)
                        flank_indices[donor] = (flanking_anchor_index(*inputs, self._positions(donor), params['max_interval'])
                                                if inputs else None)
                    flank = flank_indices[donor]
                    if not flank:
                        continue
                    region = uniquely_flanked_locus(flank, candidate['source_span'])
                    positions = self._positions(donor)
                    if region:
                        localized[donor] = {q for q in spanning_loci[key][donor] if q in positions
                                            and within_flanks(region, positions[q])}
                        if localized[donor]:
                            opposed[donor]['within_syntenic_flanks'] = True
                    if region and donor in spanning_rna_donors[key]:
                        observed = self._rna_supported_donor_models(donor)
                        count, query = max(((observed[q], q) for q in spanning_rna_donors[key][donor] & observed.keys()
                                             if q in positions and within_flanks(region, positions[q])), default=(0, ''))
                        if count >= 2:
                            opposed[donor].update(syntenic_rna_witnesses=count, syntenic_rna_query=query)
                    if region and outside_flanks(region, positions, spanning_loci[key][donor]):
                        excluded.append(donor)
                        del opposed[donor]
                        if donor in support and all(within_flanks(region, positions[h['query']])
                                                    for h in (support[donor]['first'], support[donor]['second'])):
                            unambiguous[donor] = support[donor]
                decision = fission_decision(groups, unambiguous, opposed, minimum)
                if excluded:
                    decision['off_locus_fusion_donors'] = excluded
            rna_donors = set(decision['closest_informative_species']) & opposed.keys() & spanning_rna_donors[key].keys()
            for donor in rna_donors:
                representative = opposed[donor]['query']
                verified = self._rna_supported_donor_models(donor)
                witness = representative
                count = verified.get(representative, 0)
                if donor in localized:
                    count, witness = max(((verified[q], q) for q in localized[donor] & verified.keys()),
                                         default=(0, representative))
                elif count < 2 and verified:
                    # A lower-scoring isoform of the same locus can supply the
                    # RNA witness. Evidence from other paralogous loci cannot.
                    positions = self._positions(donor)
                    position = positions.get(representative)
                    if position is not None:
                        alternatives = [(verified[candidate], candidate)
                                        for candidate in spanning_rna_donors[key][donor] & verified.keys()
                                        if (other := positions.get(candidate)) is not None and other[:3] == position[:3]]
                        count, witness = max(alternatives, default=(count, witness))
                opposed[donor]['rna_witnesses'] = count
                if witness != representative:
                    opposed[donor]['rna_supported_query'] = witness
            if rna_donors:
                decision.update(fission_decision(groups, unambiguous, opposed, minimum))
            result[key] = {**decision, 'support': list(unambiguous.values()),
                           'fusion_annotation_support': [opposed[donor] for donor in sorted(opposed)],
                           'ambiguous_donor_species': sorted(support.keys() & opposed.keys())}
        return result


def read_positions(lines):
    positions = {}
    for line in lines:
        seqid, a, b, identifier, _, strand = line.rstrip().split('\t')
        positions[identifier] = (identifier, seqid, strand, int(a), int(b))
    return positions


def outside_flanks(region, positions, hits):
    seqid, start, end = region
    return all((pos := positions.get(hit)) is not None
               and (pos[1] != seqid or pos[4] <= start or pos[3] >= end) for hit in hits)


def within_flanks(region, position):
    seqid, start, end = region
    return position[1] == seqid and start <= position[3] and position[4] <= end
