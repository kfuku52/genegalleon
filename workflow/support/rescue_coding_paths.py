"""Reconcile genomic coding paths without conflating loci or inventing sequence."""
import bisect
import math
from collections import defaultdict

try:
    from format_species_annotation.common import parse_gff_attributes
    from rescue_path_compatibility import protein_path_support
except ImportError:
    from .format_species_annotation.common import parse_gff_attributes
    from .rescue_path_compatibility import protein_path_support


def coding_shape(model):
    return model['seqid'], model['strand'], tuple(tuple(b) for b in model['cds'])


def coding_overlap(a, b):
    if a['seqid'] != b['seqid'] or a['strand'] != b['strand']:
        return 0, 0
    shared = in_frame = 0
    for left in a['cds']:
        for right in b['cds']:
            start, end = max(left[0], right[0]), min(left[1], right[1])
            if start >= end:
                continue
            shared += end - start
            if len(left) < 3 or len(right) < 3 or left[2] not in (0, 1, 2) or right[2] not in (0, 1, 2):
                continue
            first = start if a['strand'] == '+' else end - 1
            phase_a = (first - left[0] - left[2]) % 3 if a['strand'] == '+' else (left[1] - 1 - first - left[2]) % 3
            phase_b = (first - right[0] - right[2]) % 3 if b['strand'] == '+' else (right[1] - 1 - first - right[2]) % 3
            if phase_a == phase_b:
                in_frame += end - start
    return shared, in_frame


def compatible_paths(a, b, *, owners_a=(), owners_b=()):
    # Similar proteins and shared exons cannot override distinct original loci.
    if owners_a and owners_b and set(owners_a).isdisjoint(owners_b):
        return False
    shared, in_frame = coding_overlap(a, b)
    shorter = min(sum(e - s for s, e, *_ in a['cds']), sum(e - s for s, e, *_ in b['cds']))
    # A complete sequence and donor provenance are required before a new path
    # can resolve a conflict. Bare coordinate records remain review proposals.
    if not (a.get('sequence') and b.get('sequence') and alignment_support(a) and alignment_support(b)
            and shorter > 0 and shared > 0 and in_frame == shared):
        return False
    return shared / shorter >= .80 or protein_path_support(a, b)['supported']


def overlap_components(models):
    """Sweep coding intervals so intron-nested independent loci stay separate."""
    grouped = defaultdict(list)
    parent = list(range(len(models)))
    def find(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index
    def join(a, b):
        a, b = find(a), find(b)
        if a != b:
            parent[max(a, b)] = min(a, b)
    for index, model in enumerate(models):
        for start, end, *_ in model['cds']:
            grouped[model['seqid'], model['strand']].append((start, end, index))
    for records in grouped.values():
        active = []
        for start, end, index in sorted(records):
            active = [(stop, other) for stop, other in active if stop > start]
            for _, other in active:
                join(index, other)
            active.append((end, index))
    components = defaultdict(list)
    for index in range(len(models)):
        components[find(index)].append(index)
    return [components[key] for key in sorted(components)]


class OwnershipIndex:
    """Strand-aware original features; unknown strand remains conservative."""
    def __init__(self, features):
        grouped = defaultdict(list)
        for feature in features:
            grouped[feature['seqid']].append(feature)
        self.tracks = {}
        for seqid, rows in grouped.items():
            rows.sort(key=lambda r: (r['start'], r['end']))
            maxima, high = [], 0
            for row in rows:
                high = max(high, row['end'])
                maxima.append(high)
            self.tracks[seqid] = ([r['start'] for r in rows], rows, maxima)

    def overlapping(self, model, owner=None):
        starts, rows, maxima = self.tracks.get(model['seqid'], ([], [], []))
        index = bisect.bisect_left(starts, model['end']) - 1
        found = {}
        while index >= 0 and maxima[index] > model['start']:
            row = rows[index]
            strand = row.get('strand', '.')
            if (row['end'] > model['start'] and (strand not in {'+', '-'} or strand == model['strand'])
                    and (owner is None or row.get('gene_id') != owner)):
                # Intron overlap alone does not establish coding ownership.
                if strand in {'+', '-'} and row.get('cds') and not any(max(a[0], b[0]) < min(a[1], b[1])
                                             for a in row['cds'] for b in model['cds']):
                    index -= 1
                    continue
                key = row.get('gene_id') or (row['start'], row['end'])
                found[key] = row
            index -= 1
        return list(found.values())


def original_ownership(lines, attribute_parser):
    """Keep original locus IDs, strand and CDS rather than merging all spans."""
    features, parents, genes = [], {}, {}
    for line in lines:
        if line.strip() == '##FASTA':
            break
        if line.startswith('#') or not line.strip():
            continue
        fields = line.rstrip().split('\t')
        if len(fields) != 9:
            continue
        kind = fields[2]
        parsed = parse_gff_attributes(fields[8])
        attr = attribute_parser(fields[8])
        declared_genes = parsed.get('gene_id', ())
        if not (kind.lower().endswith(('gene', 'rna', 'transcript')) or kind in {'CDS', 'exon'} or declared_genes):
            continue
        identifier = next(iter(parsed.get('ID', ())), attr.get('ID', ''))
        transcripts = parsed.get('transcript_id', ())
        declared_parents = parsed.get('Parent', ())
        is_gene = kind.lower().endswith('gene')
        if not identifier and (is_gene or kind.lower().endswith(('rna', 'transcript'))):
            identifier = next(iter(declared_genes if is_gene else transcripts), '')
        lineage = declared_parents or (() if is_gene else transcripts or declared_genes)
        row = {'seqid': fields[0], 'start': int(fields[3]) - 1, 'end': int(fields[4]),
               'strand': fields[6] if fields[6] in {'+', '-'} else '.', 'id': identifier,
               'parents': lineage, 'direct_owner_ids': declared_genes, 'kind': kind}
        features.append(row)
        if is_gene and identifier:
            genes[identifier] = {**row, 'gene_id': identifier, 'cds': []}
        elif identifier:
            parents[identifier] = lineage
        for transcript in transcripts:
            if declared_genes:
                parents[transcript] = declared_genes
        # GTF can declare a locus entirely through CDS/exon attributes.
        for gene in declared_genes:
            if gene not in genes:
                genes[gene] = {**row, 'gene_id': gene, 'cds': [], 'implicit': True}
            elif genes[gene].get('implicit'):
                existing = genes[gene]
                if existing['seqid'] != row['seqid']:
                    raise ValueError('Implicit gene spans multiple contigs: ' + gene)
                existing['start'], existing['end'] = min(existing['start'], row['start']), max(existing['end'], row['end'])
                if existing['strand'] != row['strand']:
                    existing['strand'] = '.'
    # Transcript-only providers still declare a real locus. Bind its children
    # to that original identifier rather than inventing one owner per CDS.
    for row in features:
        if (row['id'] and row['kind'].lower().endswith(('rna', 'transcript'))
                and not row['parents'] and not row['direct_owner_ids']):
            genes.setdefault(row['id'], {**row, 'gene_id': row['id'], 'cds': []})
    def gene_ids(row):
        pending = [row['id']] if row['id'] in genes else [*row['parents'], *row['direct_owner_ids']]
        seen, found = set(), set()
        while pending:
            value = pending.pop()
            if not value or value in seen:
                continue
            seen.add(value)
            if value in genes:
                found.add(value)
            elif value in parents:
                pending.extend(parents[value])
        return found
    orphaned = []
    for row in features:
        owners = gene_ids(row)
        if not owners:
            orphaned.append({**row, 'gene_id': row['id'] or f"unbound_annotation:{row['seqid']}:{row['start']}:{row['end']}:{row['strand']}:{row['kind']}"})
        elif row['kind'] == 'CDS':
            for owner in owners:
                genes[owner]['cds'].append([row['start'], row['end']])
                if genes[owner]['strand'] != row['strand']:
                    genes[owner]['strand'] = '.'
    return [*genes.values(), *orphaned]


def alignment_support(model):
    """Retain each donor's checked metrics; isoforms are not extra species."""
    primary = model.get('evidence', {})
    support = model.get('support') or [primary]
    result = []
    for row in support:
        donor = row.get('donor')
        if not donor:
            continue
        metrics = row.get('alignment')
        if metrics is None:
            if (donor, row.get('query')) != (primary.get('donor'), primary.get('query')):
                continue
            metrics = model
        if 'coverage' not in metrics or 'identity' not in metrics:
            continue
        coverage, identity = metrics.get('coverage', 0), metrics.get('identity', 0)
        if (metrics.get('problems') or not all(not isinstance(v, bool) and isinstance(v, (int, float)) and math.isfinite(v) and 0 <= v <= 1
                                             for v in (coverage, identity))):
            continue
        result.append((donor, identity, coverage))
    return result


def path_rank(model, nearest=()):
    support = alignment_support(model)
    donors = {donor for donor, _, _ in support}
    identity, coverage = max(((identity, coverage) for _, identity, coverage in support), default=(0, 0))
    return len(donors & set(nearest)), len(donors), identity, coverage


def transferred_support(model):
    """A native predictor flag cannot erase a checked donor genomic mapping."""
    evidence = model.get('support') or [model.get('evidence', {})]
    return any(not (row.get('annevo') or row.get('oriongeno') or row.get('boundary_candidate'))
               and alignment_support({**model, 'support': [row]}) for row in evidence)


def compatible_components(models):
    """Hold conflicting native additions without invalidating transferred CDS.

    Reassess the transferred paths as their own coding components. This also
    prevents a native bridge from joining two otherwise independent loci. A
    conflict among transferred paths themselves still requires review.
    """
    for component in overlap_components(models):
        rows = [models[i] for i in component]
        valid = []
        for row in rows:
            if row.get('status') == 'accepted' and not row.get('problems'):
                valid.append(row)
            else:
                row['status'] = 'unresolved'
        def coherent(group):
            return all(compatible_paths(a, b) for i, a in enumerate(group) for b in group[i + 1:])
        if coherent(valid):
            yield valid
            continue
        retained = [row for row in valid if row.get('retained_coding_path')]
        if retained:
            retained_groups = [[retained[i] for i in group] for group in overlap_components(retained)]
            if not all(coherent(group) for group in retained_groups):
                raise ValueError('Previously adopted coding paths conflict after current validation')
            additions = [[] for _ in retained_groups]
            independent = []
            def hold(row, previous):
                row['status'] = 'unresolved'
                row['problems'] = sorted(set(row['problems']) | {'competing_new_models'})
                row['integration_conflict'] = {'disposition': 'retain_previous_validated_annotation',
                    'retained_model_ids': sorted(r['model_id'] for r in previous)}
            for row in valid:
                if row.get('retained_coding_path'):
                    continue
                owners = [i for i, group in enumerate(retained_groups)
                          if any(coding_overlap(row, other)[0] for other in group)]
                if not owners:
                    independent.append(row)
                    continue
                if len(owners) == 1 and all(compatible_paths(row, other) for other in retained_groups[owners[0]]):
                    additions[owners[0]].append(row)
                else:
                    hold(row, [other for i in owners for other in retained_groups[i]])
            # Existing independent genes remain separate even when a new
            # bridging prediction connected them in the combined graph.
            for previous, proposed in zip(retained_groups, additions, strict=True):
                conflicting = {id(row) for row in proposed for other in proposed
                               if row is not other and not compatible_paths(row, other)}
                for row in proposed:
                    if id(row) in conflicting:
                        hold(row, previous)
                yield [*previous, *(row for row in proposed if id(row) not in conflicting)]
            yield from compatible_components(independent)
            continue
        transferred = [row for row in valid if transferred_support(row)]
        if transferred and len(transferred) < len(valid):
            for row in valid:
                if not transferred_support(row):
                    row['status'] = 'unresolved'
                    row['problems'] = sorted(set(row['problems']) | {'competing_new_models'})
                    row['integration_conflict'] = {'disposition': 'retain_transferred_coding_paths',
                        'transferred_model_ids': sorted(r['model_id'] for r in transferred)}
            groups = [[transferred[i] for i in group] for group in overlap_components(transferred)]
        else:
            groups = [valid]
        for group in groups:
            if coherent(group):
                yield group
            else:
                for row in group:
                    row['status'] = 'unresolved'
                    row['problems'] = sorted(set(row['problems']) | {'competing_new_models'})


def resolve_coding_paths(models, nearest=()):
    """Admit a supported locus separately from ranking its coding paths.

    Callers must first validate each path's genomic CDS, donor alignment and
    original annotation ownership. A canonical fallback establishes neither
    translation initiation nor orthology; its support remains path-specific.
    """
    for rows in compatible_components(models):
        if len(rows) < 2:
            continue
        ranked = sorted(rows, key=lambda r: (path_rank(r, nearest), repr(coding_shape(r))), reverse=True)
        best, second = ranked[:2]
        first_rank, second_rank = path_rank(best, nearest), path_rank(second, nearest)
        decisive = first_rank[:2] > second_rank[:2] or first_rank[2] - second_rank[2] >= .10
        # One external species can support a complete, coherent genomic locus.
        # Locus admission does not establish a biologically preferred isoform.
        targets = {evidence['target'] for row in rows
                   for evidence in [row.get('evidence', {}), *row.get('support', [])]
                   if evidence.get('target')}
        donors = sorted({donor for row in rows for donor, _, _ in alignment_support(row)} - targets)
        retained_annotation_fallback = False
        if not decisive:
            if not donors:
                previous = [row for row in rows if row.get('retained_coding_path')]
                for row in rows:
                    if row.get('retained_coding_path'):
                        continue
                    row['status'] = 'unresolved'
                    row['problems'].append('ambiguous_coding_path_representative')
                    if previous:
                        row['integration_conflict'] = {'disposition': 'retain_previous_validated_annotation',
                            'retained_model_ids': sorted(r['model_id'] for r in previous)}
                if not previous:
                    continue
                # New self-only alternatives cannot revoke a source-verified
                # annotation. Keep only the currently valid prior paths; this
                # supplies no external evidence for the novel alternatives.
                rows = previous
                ranked = [row for row in ranked if row.get('retained_coding_path')]
                retained_annotation_fallback = True
            # A reproducible CDS choice is needed by downstream one-sequence
            # inputs. Prefer the longest validated path without claiming it
            # is the biologically preferred isoform. Tie-break on coordinates.
            best = min(rows, key=lambda row: (-sum(e - s for s, e, *_ in row['cds']),
                                              repr(coding_shape(row)), row['model_id']))
        alternatives = []
        for row in ranked:
            if row is best:
                continue
            row['status'] = 'accepted_alternative_path'
            row['parent_model_id'] = best['model_id']
            alternatives.append({k: row[k] for k in ('model_id', 'seqid', 'strand', 'cds', 'sequence', 'support', 'coverage', 'identity',
                                                    'evidence', 'placement_evidence', 'problems', 'boundary_candidate',
                                                    'boundary_revalidation', 'partial_evidence', 'retained_coding_path') if k in row})
        best['alternative_coding_paths'] = alternatives
        best['locus_support'] = {'counting_unit': 'independent_donor_species',
                                 'independent_donor_species': donors, 'minimum_species_for_fallback': 1,
                                 'orthology': 'unassigned', 'expected_copy': 'unassigned'}
        if retained_annotation_fallback:
            best['locus_support']['admission_basis'] = 'source_bound_currently_validated_prior_annotation'
        best['path_selection'] = {
            'reason': ('retained_annotation_canonical_fallback' if retained_annotation_fallback else
                       'supported_same_locus_coding_paths' if decisive else 'compatible_locus_canonical_fallback'),
            'rank': list(path_rank(best, nearest)), 'alternative_paths': len(alternatives),
            'representative_status': 'supported_priority' if decisive else 'ambiguous',
            'selection_policy': 'validated_donor_rank' if decisive else 'longest_cds_then_coding_shape',
            'ambiguity_reason': None if decisive else 'ambiguous_coding_path_representative',
            'decisive_identity_margin': .10,
            'ranked_paths': [{'model_id': row['model_id'], 'rank': list(path_rank(row, nearest)),
                              'donor_species': sorted({donor for donor, _, _ in alignment_support(row)})}
                             for row in ranked]}
