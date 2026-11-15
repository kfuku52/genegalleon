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
    from rescue_coding_paths import coding_overlap, coding_shape
    from rescue_gene_models import run, validate_model
    from rescue_homology_policy import homology_supported
except ImportError:
    from .fasta_sequence_store import fasta_records
    from .gene_model_store import iter_loci
    from .rescue_coding_paths import coding_overlap, coding_shape
    from .rescue_gene_models import run, validate_model
    from .rescue_homology_policy import homology_supported


def _in_frame_codons(fragment, path, protein_length):
    """Project complete source codons onto the native translated path."""
    reverse = path['strand'] == '-'
    bases = defaultdict(int)
    offset = 0
    for start, end, phase in path['cds']:
        for a, b, source_phase in fragment['blocks']:
            left, right = max(start, a), min(end, b)
            if left >= right:
                continue
            if (end - b - phase + source_phase if reverse else a - start - phase + source_phase) % 3:
                continue
            begin = offset + (end - right if reverse else left - start)
            for position in range(begin, begin + right - left):
                bases[position // 3] |= 1 << (position % 3)
        offset += end - start
    return {i for i, mask in bases.items() if mask == 7 and i < protein_length}


def nominate(paths, originals, ownership, genome, code, params, species, excluded_shapes, *, rna_available=False, intact_only=False):
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
        intact = [fragment for fragment in fragments if fragment['quality']['valid_orf']]
        needs_rna = len(intact) > 1
        if (needs_rna and not rna_available) or (intact_only and not needs_rna):
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
            for fragment in fragments:
                overlap, frame = coding_overlap({**fragment, 'cds': fragment['blocks']}, path)
                size = len(fragment['cds'])
                if frame < size * (.80 if needs_rna else 2 / 3):
                    break
                geometry.append({'gene_id': fragment['source_gene_id'], 'cds_bases': size, 'overlap_bases': overlap})
            else:
                if shape in candidates:
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
                    'oriongeno_transcripts': [path['native_transcript_id']], 'adoption': 'automatic',
                    'basis': ('oriongeno_cross_fragment_homology_and_target_rna' if needs_rna else
                              'oriongeno_and_cross_fragment_protein_homology')},
                    'checked': checked, 'protein': protein, 'windows': windows,
                    'protect_rna_paths': intact if needs_rna else []}
    return candidates


def nominate_extensions(groups, originals, genome, code, params, species, *, rna_available=False):
    """Recover missing terminal CDS with in-frame overlap and donor coverage of the added residues."""
    candidates = {}
    for owner, paths in groups.items():
        original = originals.get(owner)
        if not original:
            continue
        quality, source_cds = original['quality'], original['blocks']
        intact = quality['valid_orf']
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


def native_paths_by_ownership(paths, ownership, *, include_merges, extensions=None):
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
    return merge_paths, {owner: paths for owner, paths in singles.items()
                         if len({path['native_gene_id'] for path in paths}) >= 2}


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
    elif not set(closest) & split.keys():
        reason = 'closest_informative_donors_support_fusion'
    else:
        reason = 'supported_split'
    return {'reason': reason, 'closest_informative_species': closest, 'split_support_groups': split_groups}


def _paired_positions(traceback, start):
    """Read DIAMOND BTOP and omit query residues aligned to donor gaps."""
    paired = set()
    for token in re.findall(r'\d+|[^\d]{2}', traceback):
        if token.isdecimal():
            count = int(token)
            paired.update(range(start, start + count))
            start += count
        elif token[0] != '-':
            if token[1] != '-':
                paired.add(start)
            start += 1
    return paired


class DonorSearch:
    """Share one lazily built protein database across the finalization targets."""

    def __init__(self, panel, db, directory, cpus, rna_counts=None):
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

    def alignments(self, species, proteins):
        if not proteins:
            return
        if not self.directory.exists():
            self._build()
        if not self.subjects:
            return
        query = self.directory / 'queries.fa'
        output = self.directory / 'hits.tsv'
        sequences = list(dict.fromkeys(proteins))
        with query.open('w') as handle:
            for index, protein in enumerate(sequences):
                handle.write(f'>q{index}\n{protein}\n')
        run([self.panel['executable']['path'], 'blastp', '--query', query, '--db', self.database,
             '--threads', self.cpus, '--very-sensitive', '--evalue', '1e-5', '--max-target-seqs', '0',
             '--max-hsps', '1', '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'qstart', 'qend', 'qlen',
             'sstart', 'send', 'slen', 'bitscore', 'btop'], self.directory, 'blastp', output)
        allowed = set(self.panel['donors'][species])
        with output.open() as handle:
            for line in handle:
                q, s, identity, qstart, qend, qlen, sstart, send, slen, bitscore, traceback = line.split()
                donor, identifier = self.subjects[int(s[1:])]
                if donor not in allowed:
                    continue
                qstart, qend, qlen, sstart, send, slen = map(int, (qstart, qend, qlen, sstart, send, slen))
                query_coverage = (qend - qstart + 1) / qlen
                donor_coverage = (send - sstart + 1) / slen
                yield sequences[int(q[1:])], {
                    'donor': donor, 'query': identifier, 'identity': float(identity) / 100,
                    'coverage': min(query_coverage, donor_coverage), 'query_coverage': query_coverage,
                    'donor_coverage': donor_coverage, 'bitscore': float(bitscore),
                }, traceback, qstart - 1, max(qlen, slen)

    def support(self, species, candidates, params):
        queries = defaultdict(list)
        for shape, candidate in candidates.items():
            queries[candidate['protein']].append((shape, candidate['windows']))
        best = defaultdict(dict)
        for protein, hit, traceback, start, length in self.alignments(species, queries):
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
                previous = best[shape].get(donor)
                if previous is None or (hit['bitscore'], hit['query']) > (previous['bitscore'], previous['query']):
                    best[shape][donor] = {**hit, 'paired_coverage': paired_coverage,
                                         'fragment_coverage': fragment_coverage, 'flags': []}
        return {shape: [hits[donor] for donor in sorted(hits)] for shape, hits in best.items()
                if len(hits) >= max(2, params['min_support'])}

    def _positions(self, donor):
        if donor not in self.positions:
            bed = self.panel.get('positions', {}).get(donor)
            positions = {}
            if bed:
                with Path(bed).open() as handle:
                    for line in handle:
                        seqid, a, b, identifier, _, strand = line.rstrip().split('\t')
                        positions[identifier] = (identifier, seqid, strand, int(a), int(b))
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

    def fission_support(self, species, candidates, params):
        proteins = {protein for candidate in candidates.values()
                    for protein in (candidate['source_protein'], *candidate['proteins'])}
        source_windows = defaultdict(dict)
        for key, candidate in candidates.items():
            source_windows[candidate['source_protein']][key] = candidate['source_windows']
        best = defaultdict(dict)
        spanning = defaultdict(dict)
        for protein, hit, traceback, start, length in self.alignments(species, proteins):
            if hit['identity'] < .45 or hit['coverage'] < .65:
                continue
            paired = _paired_positions(traceback, start)
            proportion = len(paired) / length
            if proportion < .60:
                continue
            for key, windows in source_windows.get(protein, {}).items():
                if hit['identity'] >= .50 and proportion >= .80 and all(
                        len(paired & window) >= .80 * len(window) for window in windows):
                    donor = hit['donor']
                    previous = spanning[key].get(donor)
                    if previous is None or (hit['bitscore'], hit['query']) > (previous['bitscore'], previous['query']):
                        spanning[key][donor] = {**hit, 'paired_coverage': proportion}
            donor = hit['donor']
            previous = best[protein].get(donor)
            if previous is None or hit['bitscore'] > previous['bitscore']:
                best[protein][donor] = {**hit, 'paired_coverage': proportion, 'tied': False}
            elif hit['bitscore'] == previous['bitscore'] and hit['query'] != previous['query']:
                previous['tied'] = True

        result = {}
        groups = self.panel.get('sister_clades', {}).get(species, [])
        minimum = max(2, params['min_support'])
        for key, candidate in candidates.items():
            opposed = spanning[key]
            first, second = (best[protein] for protein in candidate['proteins'])
            support = {}
            for donor in sorted(first.keys() & second.keys()):
                a, b = first[donor], second[donor]
                if a['tied'] or b['tied']:
                    continue
                positions = self._positions(donor)
                gap = collinear_gap(positions.get(a['query']), positions.get(b['query']))
                if gap is not None and gap <= 10_000:
                    support[donor] = {'species': donor, 'gap_bases': gap,
                                      'first': {k: v for k, v in a.items() if k != 'tied'},
                                      'second': {k: v for k, v in b.items() if k != 'tied'}}
            # Treat donor species with both architectures as ambiguous at this locus.
            unambiguous = {name: hit for name, hit in support.items() if name not in opposed}
            decision = fission_decision(groups, unambiguous, opposed, minimum)
            for donor in set(decision['closest_informative_species']) & opposed.keys() & self.rna_counts.keys():
                count = self._rna_supported_donor_models(donor).get(opposed[donor]['query'], 0)
                opposed[donor]['rna_witnesses'] = count
            decision = fission_decision(groups, unambiguous, opposed, minimum)
            result[key] = {**decision, 'support': list(unambiguous.values()),
                           'fusion_annotation_support': [opposed[donor] for donor in sorted(opposed)],
                           'ambiguous_donor_species': sorted(support.keys() & opposed.keys())}
        return result
