"""Review modeled origins using saved focal classifications and species branches.

These are attention flags, never a hard taxonomy filter or a contamination
diagnosis. Best-hit annotations and neighbor classifications are not substituted
for the exact event-linked focal gene. A reusable evidence scope reads each
species source once, retains only requested genes, and verifies its input hashes.
"""

import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path

from focus_hgt_direction import key, known, read_taxonomy
from focus_hgt_gene_trees import number, validate_link_identity

SIDES = ('donor', 'recipient')
STATUSES = ('compatible', 'incompatible', 'unresolved', 'missing', 'conflicting')
RANKS = ('domain', 'kingdom', 'class')
SCAFFOLD_RANKS = {'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'}
RANK_ORDER = ('domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
RANK_POSITIONS = {rank: index for index, rank in enumerate(RANK_ORDER)}
PAIR_IDENTITY_FIELDS = ('event_id', 'orthogroup', 'event_index', 'generax_transfer', 'generax_donor_node',
                        'generax_recipient_node', 'branch_id', 'gene_tree_branch_id', 'node_name', 'gene_tree_node')
EVENT_FIELDS = ['origin_review_status', 'origin_review_flags', 'origin_endpoint_relation',
                'origin_donor_taxonomy_status', 'origin_recipient_taxonomy_status', 'origin_passing_pair_count',
                *(f'origin_{side}_class_{status}_pair_count' for side in SIDES for status in STATUSES),
                *(f'origin_{side}_higher_taxon_mismatch_pair_count' for side in SIDES)]
PAIR_FIELDS = [*(f'{side}_origin_host_class_status' for side in SIDES),
               *(f'{side}_origin_highest_mismatch_rank' for side in SIDES),
               *(f'{side}_origin_attention_flags' for side in SIDES), 'origin_pair_review_flags']
GENE_FIELDS = ['origin_gene_species', 'origin_gene_species_status', 'origin_mmseqs2_status',
               'origin_mmseqs2_lca_taxid', 'origin_mmseqs2_lca_rank', 'origin_mmseqs2_lca_name',
               'origin_mmseqs2_lineage_taxids', 'origin_mmseqs2_source', 'origin_mmseqs2_source_sha256',
               *(f'origin_host_{rank}_taxid' for rank in RANKS),
               *(f'origin_host_{rank}_taxid_candidates' for rank in RANKS),
               *(f'origin_host_{rank}_status' for rank in RANKS),
               'origin_scaffold_domain_label', 'origin_scaffold_class_label',
               'origin_scaffold_taxonomy_source', 'origin_scaffold_taxonomy_source_sha256',
               'origin_highest_mismatch_rank', 'origin_gene_attention_flags',
               'origin_expression_evidence_status', 'origin_intron_evidence_status', 'origin_synteny_evidence_status']


def numeric_taxid(value, zero=False):
    value = str(value).strip()
    if not known(value):
        return None
    try:
        number = Decimal(value)
    except InvalidOperation:
        raise ValueError('Invalid saved origin taxid: ' + value) from None
    if not number.is_finite() or number < (0 if zero else 1) or number != number.to_integral_value():
        raise ValueError('Invalid saved origin taxid: ' + value)
    return int(number)


def taxid_values(value):
    """Read singleton or JSON-array cells emitted by species_taxonomy.pack_values."""
    value = str(value).strip()
    if not known(value):
        return ()
    if value.startswith('['):
        try:
            values = json.loads(value)
        except ValueError:
            raise ValueError('Malformed packed origin taxid values: ' + value) from None
        if not isinstance(values, list):
            raise ValueError('Invalid packed origin taxid values: ' + value)
        result = tuple(numeric_taxid(item) for item in values)
        if any(item is None for item in result):
            raise ValueError('Missing packed origin taxid value: ' + value)
        return result
    return (numeric_taxid(value),)


def rank_names(value):
    value = str(value).strip()
    if not known(value):
        return ()
    if re.match(r'^\[\s*"', value) or value == '[]':
        try:
            values = json.loads(value)
        except ValueError:
            raise ValueError('Malformed packed origin rank names: ' + value) from None
        if not isinstance(values, list) or any(not isinstance(item, str) for item in values):
            raise ValueError('Invalid packed origin rank names: ' + value)
        return tuple(item.strip() for item in values if known(item))
    # Scientific names can begin with square brackets; the producer preserves
    # singleton names verbatim rather than JSON-encoding them.
    return (value,)


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()


def fingerprint(path):
    stat = Path(path).stat()
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


def rank_below_or_equal(rank, reference):
    """Only explicit comparable ranks establish a negative rank comparison."""
    rank = {'superkingdom': 'domain', 'species group': 'species', 'species subgroup': 'species',
            'subspecies': 'species', 'strain': 'species'}.get(rank, rank)
    return rank in RANK_POSITIONS and RANK_POSITIONS[rank] >= RANK_POSITIONS[reference]


class OriginEvidence:
    """One analysis scope; no process-global cache or network/database lookup.

    Construct with the complete requested event-gene links before processing
    family batches. Source files are loaded lazily and retained records are
    restricted to these exact genes. Call ``finalize_report`` after the final
    batch to check every consumed source again.
    """

    def __init__(self, links=(), mmseqs2_taxonomy_dir='', scaffold_taxonomy_dir='', species_taxonomy=''):
        self.sequence_root = Path(mmseqs2_taxonomy_dir).resolve() if mmseqs2_taxonomy_dir else None
        self.scaffold_root = Path(scaffold_taxonomy_dir).resolve() if scaffold_taxonomy_dir else None
        self.sources, self.fingerprints, self.missing_sources = {}, {}, set()
        self.wanted, self.cache, self.gene_cache = defaultdict(set), {}, {}
        self.host_cache, self.branch_cache, self.relation_cache = {}, {}, {}
        self.host_rank_candidates, self.ambiguous_host_ranks = {}, {}
        self.taxonomy, self.rank_taxids = {}, defaultdict(set)
        gene_species = {}
        for link in links:
            species, gene = key(link.get('gene_species', '')), str(link.get('gene_id', '')).strip()
            if species and gene:
                if not re.fullmatch(r'[A-Za-z0-9_.-]+', species) or species in {'.', '..'}:
                    raise ValueError('Unsafe origin taxonomy species identifier')
                if gene in gene_species and gene_species[gene] != species:
                    raise ValueError('Ambiguous exact origin gene/species identity')
                gene_species[gene] = species
                self.wanted[species].add(gene)
        if species_taxonomy:
            path = Path(species_taxonomy).resolve()
            before = fingerprint(path)
            self.taxonomy, sha = read_taxonomy(path)
            if fingerprint(path) != before:
                raise ValueError('Origin species taxonomy changed during reading')
            self.sources[str(path)], self.fingerprints[str(path)] = sha, before
            for row in self.taxonomy.values():
                if row['resolution_status'] != 'resolved':
                    continue
                for rank in RANKS:
                    self.rank_taxids[rank].update(taxid_values(row.get(rank + '_taxid', '')))

    def _rows(self, path):
        """Hash the bytes used by a streaming CSV reader, with a stat fence."""
        before, result = fingerprint(path), hashlib.sha256()
        with path.open('rb') as handle:
            def lines():
                for index, line in enumerate(handle):
                    result.update(line)
                    yield line.decode('utf-8-sig' if index == 0 else 'utf-8')
            yield from csv.reader(lines(), delimiter='\t')
        if fingerprint(path) != before:
            raise ValueError('Origin taxonomy source changed during reading: ' + str(path))
        self.sources[str(path)], self.fingerprints[str(path)] = result.hexdigest(), before

    def _source(self, root, species, suffix):
        if root is None:
            return None, 'source_not_configured'
        path = root / (species + suffix)
        if not path.is_file():
            self.missing_sources.add(str(path))
            return None, 'source_file_unavailable'
        return path, 'gene_record_unavailable'

    def _load(self, species):
        if species in self.cache:
            return self.cache[species]
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', species) or species in {'.', '..'}:
            raise ValueError('Unsafe origin taxonomy species identifier')
        assignments, labels = {}, {}
        wanted = self.wanted[species]
        sequence, status = self._source(self.sequence_root, species, '_mmseqs2taxonomy.tsv')
        if sequence:
            for row in self._rows(sequence):
                if len(row) < 4 or not row[0]:
                    raise ValueError('Malformed saved origin MMseqs2 classification row')
                if row[0] not in wanted:
                    continue
                gene = row[0]
                if gene in assignments:
                    raise ValueError('Duplicate requested origin MMseqs2 query gene ID')
                assigned = numeric_taxid(row[1], zero=True)
                if assigned is None or not row[2].strip() or not row[3].strip():
                    raise ValueError('Missing saved origin MMseqs2 classification')
                lineage = tuple(numeric_taxid(value) for value in row[8].split(';')) if len(row) >= 9 and known(row[8]) else ()
                if lineage and (None in lineage or len(lineage) != len(set(lineage)) or lineage[-1] != assigned):
                    raise ValueError('Saved origin MMseqs2 lineage disagrees with query LCA')
                assignments[gene] = (assigned, row[2], row[3], lineage)
        scaffold, scaffold_status = self._source(self.scaffold_root, species, '_gene_taxonomy.tsv')
        if scaffold:
            rows = iter(self._rows(scaffold))
            fields = next(rows, [])
            required = {'species', 'gene_id', 'scaffold', 'locus_id', 'count_unit', 'rank', 'host_taxid', 'label'}
            if len(fields) != len(set(fields)) or not required <= set(fields):
                raise ValueError('Malformed origin scaffold taxonomy columns')
            indices = {field: fields.index(field) for field in required}
            host_ids, parsed_host_ids = {}, {}
            for values in rows:
                if len(values) != len(fields):
                    raise ValueError('Malformed origin scaffold taxonomy row')
                if key(values[indices['species']]) != species:
                    raise ValueError('Origin scaffold taxonomy has a different species')
                rank, label = values[indices['rank']], values[indices['label']]
                if rank not in SCAFFOLD_RANKS or label not in {'compatible', 'incompatible', 'unresolved'}:
                    raise ValueError('Invalid saved origin scaffold taxonomy rank or label')
                raw_host = values[indices['host_taxid']]
                if raw_host not in parsed_host_ids:
                    parsed_host_ids[raw_host] = numeric_taxid(raw_host)
                host = parsed_host_ids[raw_host]
                if host is None and label != 'unresolved':
                    raise ValueError('Resolved origin scaffold label without host taxid')
                if rank in host_ids and host_ids[rank] != host:
                    raise ValueError('Conflicting origin scaffold host taxids')
                host_ids[rank] = host
                gene = values[indices['gene_id']]
                if gene not in wanted:
                    continue
                if any(not values[indices[field]] for field in required - {'host_taxid'}):
                    raise ValueError('Empty requested origin scaffold gene identity')
                if values[indices['count_unit']] not in {'gff_locus', 'cds_id'}:
                    raise ValueError('Invalid origin scaffold counting unit')
                identity = tuple(values[indices[field]] for field in ('scaffold', 'locus_id', 'count_unit'))
                records = labels.setdefault(gene, {})
                if rank in records:
                    raise ValueError('Duplicate requested origin scaffold gene rank')
                if records and next(iter(records.values()))[2] != identity:
                    raise ValueError('Conflicting requested origin scaffold gene identity')
                records[rank] = (host, label, identity)
            for records in labels.values():
                if set(records) != SCAFFOLD_RANKS:
                    raise ValueError('Incomplete requested origin scaffold gene ranks')
        result = assignments, labels, sequence, status, scaffold, scaffold_status
        self.cache[species] = result
        return result

    def _rank_status(self, rank, host, host_lineage, assignment, saved):
        if assignment is None:
            return 'missing'
        assigned, query_rank, _, lineage = assignment
        if assigned == 0:
            return 'unresolved'
        if host is None:
            direct = 'unresolved'
        elif assigned == host or host in lineage:
            direct = 'compatible'
        elif assigned in host_lineage:
            direct = 'unresolved'
        elif not self.rank_taxids[rank].isdisjoint(lineage) or (
                query_rank in {rank, 'superkingdom' if rank == 'domain' else rank}):
            direct = 'incompatible'
        elif lineage and rank_below_or_equal(query_rank, rank):
            direct = 'incompatible'
        else:
            direct = 'unresolved'
        if saved and saved[1] != 'unresolved':
            if host is not None and saved[0] != host:
                return 'conflicting'
            if direct in {'compatible', 'incompatible'} and direct != saved[1]:
                return 'conflicting'
            if direct == 'unresolved' and assigned in host_lineage and assigned != host:
                # A focal LCA above the tested host rank cannot establish that
                # rank. A resolved saved label here describes different evidence.
                return 'conflicting'
            return saved[1]
        return direct

    def _host(self, species):
        if species not in self.host_cache:
            row = self.taxonomy.get(species, {})
            resolved = row.get('resolution_status') == 'resolved'
            lineage = {taxid for field, value in row.items() if field.endswith('_taxid')
                       for taxid in taxid_values(value)} if resolved else set()
            candidates = {rank: tuple(dict.fromkeys(taxid_values(row.get(rank + '_taxid', '')))) if resolved else ()
                          for rank in RANKS}
            ambiguous = {rank for rank, values in candidates.items() if len(values) > 1
                         or resolved and len(set(rank_names(row.get(rank, '')))) > 1}
            ranks = {rank: values[0] if len(values) == 1 and rank not in ambiguous else None
                     for rank, values in candidates.items()}
            self.host_rank_candidates[species] = candidates
            self.ambiguous_host_ranks[species] = ambiguous
            self.host_cache[species] = lineage, ranks
        return self.host_cache[species]

    def gene(self, gene, species):
        """Return the same immutable-by-contract evidence mapping on reuse."""
        identity = gene, species
        if identity in self.gene_cache:
            return self.gene_cache[identity]
        result = dict.fromkeys(GENE_FIELDS, '')
        result.update(origin_gene_species=species, origin_gene_species_status='matched_exact_link' if species else 'gene_species_unavailable')
        if not species:
            result.update(origin_mmseqs2_status='gene_species_unavailable', origin_gene_attention_flags='focal_gene_species_unavailable')
            for rank in RANKS:
                result[f'origin_host_{rank}_status'] = 'missing'
            self.gene_cache[identity] = result
            return result
        assignments, labels, sequence, source_status, scaffold, _ = self._load(species)
        assignment, saved = assignments.get(gene), labels.get(gene, {})
        result['origin_mmseqs2_status'] = source_status if assignment is None else 'assigned' if assignment[0] else 'unclassified'
        if sequence:
            result.update(origin_mmseqs2_source=str(sequence), origin_mmseqs2_source_sha256=self.sources[str(sequence)])
        if scaffold:
            result.update(origin_scaffold_taxonomy_source=str(scaffold),
                          origin_scaffold_taxonomy_source_sha256=self.sources[str(scaffold)])
        if assignment is not None:
            result.update(origin_mmseqs2_lca_taxid=str(assignment[0]), origin_mmseqs2_lca_rank=assignment[1],
                          origin_mmseqs2_lca_name=assignment[2], origin_mmseqs2_lineage_taxids=';'.join(map(str, assignment[3])))
        host_lineage, host_ranks = self._host(species)
        flags = []
        for rank in RANKS:
            host = host_ranks[rank]
            rank_saved = saved.get(rank)
            ambiguous = rank in self.ambiguous_host_ranks[species]
            if host is None and rank_saved and not ambiguous:
                host = rank_saved[0]
            result[f'origin_host_{rank}_taxid'] = str(host) if host else ''
            result[f'origin_host_{rank}_taxid_candidates'] = ';'.join(map(str, self.host_rank_candidates[species][rank]))
            result[f'origin_host_{rank}_status'] = ('missing' if assignment is None else 'unresolved') if ambiguous else \
                self._rank_status(rank, host, host_lineage, assignment, rank_saved)
            if ambiguous:
                flags.append('host_' + rank + '_rank_ambiguous')
            if rank in {'domain', 'class'}:
                result[f'origin_scaffold_{rank}_label'] = rank_saved[1] if rank_saved else ''
            if result[f'origin_host_{rank}_status'] == 'conflicting':
                flags.append('focal_' + rank + '_classification_sources_conflict')
        mismatches = [rank for rank in RANKS if result[f'origin_host_{rank}_status'] == 'incompatible']
        result['origin_highest_mismatch_rank'] = mismatches[0] if mismatches else ''
        if result['origin_highest_mismatch_rank']:
            flags.append('focal_' + result['origin_highest_mismatch_rank'] + '_mismatch')
        if result['origin_host_class_status'] in {'missing', 'unresolved'}:
            flags.append('focal_class_classification_' + result['origin_host_class_status'])
        if assignment and assignment[0] and not assignment[3]:
            flags.append('saved_focal_mmseqs2_lineage_unavailable')
        result['origin_gene_attention_flags'] = '; '.join(flags)
        self.gene_cache[identity] = result
        return result

    def branch_taxonomy(self, tips):
        if tips is None or not tips:
            return 'species_branch_unmapped'
        identity = tuple(tips)
        if identity not in self.branch_cache:
            self.branch_cache[identity] = self._branch_taxonomy(tips)
        return self.branch_cache[identity]

    def _branch_taxonomy(self, tips):
        if not self.taxonomy:
            return 'taxonomy_source_unavailable'
        rows = [self.taxonomy.get(key(tip)) for tip in tips]
        if any(not row or row.get('resolution_status') != 'resolved' for row in rows):
            return 'unresolved_taxonomy'
        for rank in RANKS:
            names = [set(rank_names(row.get(rank, ''))) for row in rows]
            ids = [set(taxid_values(row.get(rank + '_taxid', ''))) for row in rows]
            if any(len(values) > 1 for values in names + ids):
                return 'ambiguous_' + rank + '_clade'
            values = set().union(*names)
            id_values = set().union(*ids)
            if len(values) > 1 or len(id_values) > 1:
                return 'mixed_' + rank + '_clade'
            if values and any(not value for value in names) or id_values and any(not value for value in ids):
                return 'partially_unresolved_' + rank + '_clade'
        return 'resolved_uniform_clade' if all(rank_names(row.get('domain', '')) for row in rows) else 'unresolved_taxonomy'

    def endpoint_relation(self, donor, recipient):
        identity = tuple(donor or ()), tuple(recipient or ())
        if identity not in self.relation_cache:
            self.relation_cache[identity] = endpoint_relation(donor, recipient)
        return self.relation_cache[identity]

    def finalize_report(self):
        for path, expected in self.sources.items():
            if fingerprint(path) != self.fingerprints[path] or digest(path) != expected:
                raise ValueError('Origin taxonomy input changed during generation: ' + path)
        if any(Path(path).is_file() for path in self.missing_sources):
            raise ValueError('Unavailable origin taxonomy input appeared during generation')
        return dict(source_sha256=dict(self.sources), missing_sources=sorted(self.missing_sources),
                    classification_scope='exact_event_linked_focal_gene_not_best_hit_or_neighbor',
                    attention_flags_are_exclusion_criteria=False, source_snapshot_verified=True,
                    loaded_species_count=len(self.cache), cached_gene_count=len(self.gene_cache))


def endpoint_relation(donor, recipient):
    if not donor or not recipient:
        return 'species_branch_unmapped'
    donor, recipient = frozenset(donor), frozenset(recipient)
    if donor == recipient:
        return 'same_species_clade'
    if recipient < donor:
        return 'donor_ancestor_of_recipient'
    if donor < recipient:
        return 'recipient_ancestor_of_donor'
    return 'incomparable_species_clades' if donor.isdisjoint(recipient) else 'overlapping_species_clades'


def auxiliary_status(link, field):
    value = link.get(field, '')
    if not known(value):
        return 'unavailable'
    if field == 'expression_measured':
        return ('measured' if str(value).strip().lower() in {'true', '1'} else 'not_measured'
                if str(value).strip().lower() in {'false', '0'} else 'availability_not_established')
    if field == 'intron_supported':
        counts = [number(link.get(name, '')) for name in ('num_intron', 'intron_count')]
        if any(value is not None and value >= 0 and value.is_integer() for value in counts):
            return 'observed_count_available'
        return 'recorded_support_flag' if str(value).strip().lower() in {'true', 'false', '1', '0'} else 'availability_not_established'
    return 'recorded_score'


def annotate_origin(events, pairs, genes, nodes, links=(), mmseqs2_taxonomy_dir='', scaffold_taxonomy_dir='',
                    species_taxonomy='', evidence=None):
    """Add review fields in place to owned rows, retaining every input event.

    Pair counts and flags use only exact pairs passing the combined pair filter
    (legacy ``passes_pfam_filter`` when the new field is absent). A passing copy
    never changes another copy's focal classification. With a reusable evidence
    scope the caller must invoke ``evidence.finalize_report()`` after all batches.
    """
    ids = {row['event_id']: row for row in events}
    if len(ids) != len(events) or any(not str(value).strip() for value in ids):
        raise ValueError('Duplicate or empty origin event identity')
    for row in events:
        if row['generax_transfer'] != f"Y@{row['generax_donor_node']}@{row['generax_recipient_node']}":
            raise ValueError('Origin event transfer token disagrees with species branches')
    for rows, fields in ((events, EVENT_FIELDS), (pairs, PAIR_FIELDS), (genes, GENE_FIELDS)):
        if any(any(field in row for field in fields) for row in rows):
            raise ValueError('Reserved origin audit fields already exist')
    linked, identities = {}, set()
    for link in links:
        if link['event_id'] not in ids:
            continue
        identity = link['event_id'], link['side'], link['gene_id']
        if identity in identities or link['side'] not in SIDES:
            raise ValueError('Duplicate or invalid origin event-gene link')
        identities.add(identity)
        validate_link_identity(ids[link['event_id']], link)
        linked[identity] = link
    gene_rows = {}
    for row in genes:
        identity = row['event_id'], row['side'], row['gene_id']
        if row['event_id'] not in ids or row['side'] not in SIDES or identity in gene_rows:
            raise ValueError('Duplicate, invalid or unmapped origin gene audit row')
        validate_link_identity(ids[row['event_id']], row)
        gene_rows[identity] = row
    pair_identities = set()
    for row in pairs:
        if row['event_id'] not in ids:
            raise ValueError('Unmapped origin pair event identity')
        pair_identity = row['event_id'], row['donor_gene_id'], row['recipient_gene_id']
        if pair_identity in pair_identities:
            raise ValueError('Duplicate origin event-gene pair')
        pair_identities.add(pair_identity)
        pair_identity_row = {field: row[field] for field in PAIR_IDENTITY_FIELDS if field in row}
        pair_identity_row['gene_id'] = row['donor_gene_id']
        validate_link_identity(ids[row['event_id']], pair_identity_row)
        for side in SIDES:
            identity = row['event_id'], side, row[side + '_gene_id']
            if identity not in gene_rows:
                raise ValueError('Origin pair has no exact event-side gene audit row')
            link = linked.get(identity, {})
            if link.get('lineage_status', 'retained') != 'retained':
                raise ValueError('Origin pair uses a non-retained event gene')
            species = key(link.get('gene_species', ''))
            branch = key(ids[row['event_id']][f'generax_{side}_node'])
            if species and branch in nodes and species not in nodes[branch]:
                raise ValueError('Origin pair gene is outside its modeled species branch')
        passing = row.get('passes_pair_filter', row.get('passes_pfam_filter', ''))
        if str(passing).lower() not in {'true', 'false', '1', '0'}:
            raise ValueError('Missing or invalid origin pair filtering decision')
    scope = evidence or OriginEvidence(links, mmseqs2_taxonomy_dir, scaffold_taxonomy_dir, species_taxonomy)
    summaries, event_flags = defaultdict(Counter), defaultdict(set)
    assessments = {}
    for identity, row in gene_rows.items():
        link = linked.get(identity, {})
        species = key(link.get('gene_species', row.get('gene_species', '')))
        # A reusable scope must have received these requested exact identities.
        if species and row['gene_id'] not in scope.wanted[species]:
            raise ValueError('Origin evidence scope lacks the requested species/gene')
        assessment = scope.gene(row['gene_id'], species)
        assessments[identity] = assessment
        row.update(assessment)
        for field, target in [('expression_measured', 'expression'), ('intron_supported', 'intron'),
                              ('synteny_support_score', 'synteny')]:
            row[f'origin_{target}_evidence_status'] = auxiliary_status(link, field)
    for row in pairs:
        flags = set()
        passing = str(row.get('passes_pair_filter', row.get('passes_pfam_filter', ''))).lower() in {'true', '1'}
        event_id = row['event_id']
        summaries[event_id]['passing'] += int(passing)
        for side in SIDES:
            assessment = assessments[event_id, side, row[side + '_gene_id']]
            status, mismatch = assessment['origin_host_class_status'], assessment['origin_highest_mismatch_rank']
            row[side + '_origin_host_class_status'] = status
            row[side + '_origin_highest_mismatch_rank'] = mismatch
            row[side + '_origin_attention_flags'] = assessment['origin_gene_attention_flags']
            flags.update(side + '_' + flag.strip() for flag in assessment['origin_gene_attention_flags'].split(';') if flag.strip())
            if passing:
                summaries[event_id][side, status] += 1
                summaries[event_id][side, 'higher_mismatch'] += int(mismatch in {'domain', 'kingdom'})
                event_flags[event_id].update(side + '_' + flag.strip()
                                             for flag in assessment['origin_gene_attention_flags'].split(';') if flag.strip())
        row['origin_pair_review_flags'] = '; '.join(sorted(flags))
    for row in events:
        donor, recipient = (nodes.get(key(row[f'generax_{side}_node'])) for side in SIDES)
        relation = scope.endpoint_relation(donor, recipient)
        flags = event_flags[row['event_id']]
        if relation != 'incomparable_species_clades':
            flags.add(relation + '_origin_review')
        count = summaries[row['event_id']]
        row.update(origin_endpoint_relation=relation, origin_passing_pair_count=count['passing'])
        for side, tips in (('donor', donor), ('recipient', recipient)):
            branch_status = scope.branch_taxonomy(tips)
            row[f'origin_{side}_taxonomy_status'] = branch_status
            if branch_status != 'resolved_uniform_clade':
                flags.add(side + '_' + branch_status + '_origin_review')
            for status in STATUSES:
                row[f'origin_{side}_class_{status}_pair_count'] = count[side, status]
            row[f'origin_{side}_higher_taxon_mismatch_pair_count'] = count[side, 'higher_mismatch']
            if count['passing'] and count[side, 'incompatible'] == count['passing']:
                flags.add('all_passing_' + side + '_pairs_class_incompatible')
            elif count[side, 'incompatible']:
                flags.add('some_passing_' + side + '_pairs_class_incompatible')
        if not count['passing']:
            flags.add('focal_classification_not_assessed_on_a_passing_pair')
        row.update(origin_review_status='review_required' if flags else 'no_origin_attention_flag_detected',
                   origin_review_flags='; '.join(sorted(flags)))
    report = scope.finalize_report() if evidence is None else dict(source_sha256=dict(scope.sources),
                                                                  source_snapshot_verified=False,
                                                                  attention_flags_are_exclusion_criteria=False)
    report.update(event_count=len(events), passing_pair_count=sum(row['origin_passing_pair_count'] for row in events),
                  review_required_event_count=sum(row['origin_review_status'] == 'review_required' for row in events))
    return events, pairs, genes, report
