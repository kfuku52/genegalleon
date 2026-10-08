"""Per-gene existing annotations for focused genomic context pages."""

import csv
import hashlib
import io
import re
import sqlite3
from contextlib import closing
from decimal import Decimal, InvalidOperation
from functools import lru_cache
from pathlib import Path

RANKS = ('kingdom', 'phylum', 'class', 'order', 'family', 'genus')
TABLE_EDGES = (0, .14, .36, .70, 1)
TABLE_WIDTH_PT = .435 * 22 * 72
SEQUENCE_FIELDS = ('mmseqs2_status', 'mmseqs2_lca_taxid', 'mmseqs2_lca_rank', 'mmseqs2_lca_name',
                   'mmseqs2_source', 'mmseqs2_source_sha256', 'mmseqs2_host_class_label',
                   'mmseqs2_host_species_label', 'mmseqs2_host_species_taxid',
                   'mmseqs2_host_label_status', 'mmseqs2_host_label_source', 'mmseqs2_host_label_source_sha256',
                   'mmseqs2_lineage_taxids', 'mmseqs2_lineage_status', 'mmseqs2_taxonomy_source',
                   'mmseqs2_taxonomy_source_sha256', *(f'mmseqs2_{rank}' for rank in RANKS))
FIELDS = ('gene_id', 'orthogroup', 'protein_product_name', 'protein_product_status',
          'protein_product_source', 'protein_product_source_sha256', 'protein_product_feature_ids',
          'protein_product_mapping_scope', 'gene_description', 'swissprot_best_hit_protein_name',
          'swissprot_name_source', 'swissprot_name_source_sha256',
          'besthit_accession', 'besthit_organism', 'besthit_taxid', 'besthit_source',
          'besthit_source_sha256', 'besthit_status', 'annotation_validation_status', 'besthit_coverage_percent', 'besthit_identity_percent',
          'besthit_evalue', 'taxonomy_source', 'taxonomy_source_sha256',
          *(f'besthit_{rank}' for rank in RANKS), *SEQUENCE_FIELDS)


def available(value):
    return str(value).strip() if value is not None and str(value).strip().lower() not in {
        '', '.', 'na', 'nan', 'none', 'null', 'unavailable', 'annotation unavailable'} else ''


def taxid(value):
    """Saved numeric TSVs may have '.0'; fractional/invalid IDs are never rounded."""
    value = available(value)
    if not value:
        return None
    try:
        number = Decimal(value)
    except InvalidOperation as exc:
        raise ValueError('Invalid best-hit taxid: ' + value) from exc
    if not number.is_finite() or number <= 0 or number != number.to_integral_value():
        raise ValueError('Invalid best-hit taxid: ' + value)
    return int(number)


def source_digest(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


class ContextAnnotations:
    """Exact gene IDs only; neighbor genes never inherit focal annotations."""

    def __init__(self, path='', store=None, mmseqs2_taxonomy_dir='', scaffold_taxonomy_dir='', taxonomy_dbfile=''):
        self.rows, self.sources, self.display_audit = {}, {}, []
        self.store, self.leaf_cache, self.family_sources = store, {}, {}
        self.sequence_root = Path(mmseqs2_taxonomy_dir) if mmseqs2_taxonomy_dir else None
        self.host_root = Path(scaffold_taxonomy_dir) if scaffold_taxonomy_dir else None
        self.sequence_cache = {}
        self._annotation_source = None
        self.taxonomy_path = Path(taxonomy_dbfile).resolve() if taxonomy_dbfile else None
        if self.taxonomy_path is not None and not self.taxonomy_path.is_file():
            raise FileNotFoundError(self.taxonomy_path)
        if path:
            path = Path(path).resolve()
            self._annotation_source = str(path)
            raw = path.read_bytes()
            self.sources[str(path)] = hashlib.sha256(raw).hexdigest()
            reader = csv.DictReader(io.StringIO(raw.decode('utf-8-sig')), delimiter='\t')
            if reader.fieldnames and len(set(reader.fieldnames)) != len(reader.fieldnames):
                raise ValueError('Duplicate context annotation columns')
            required = {'gene_id', 'orthogroup', 'besthit_accession', 'besthit_organism',
                        *[f'besthit_{r}' for r in RANKS]}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError('Context annotations lack required per-gene columns')
            for row in reader:
                if None in row or any(value is None for value in row.values()):
                    raise ValueError('Malformed context annotation TSV row')
                gene = available(row['gene_id'])
                if not gene or gene in self.rows:
                    raise ValueError('Duplicate or empty context annotation gene ID')
                family = available(row['orthogroup'])
                if not re.fullmatch(r'[A-Za-z0-9_.-]+', family) or family in {'.', '..'}:
                    raise ValueError('Empty or unsafe context annotation orthogroup')
                taxid(row.get('besthit_taxid'))
                if (available(row['besthit_organism']) or available(row.get('swissprot_best_hit_protein_name'))
                        or available(row.get('besthit_taxid'))
                        or any(available(row[f'besthit_{r}']) for r in RANKS)) \
                        and not available(row['besthit_accession']):
                    raise ValueError('Best-hit organism/taxonomy requires the same hit accession')
                self.rows[gene] = row

    @lru_cache(maxsize=None)
    def query_lineage(self, assigned, saved_lineage):
        result = dict.fromkeys((f'mmseqs2_{rank}' for rank in RANKS), '')
        result['mmseqs2_lineage_status'] = 'unclassified' if not assigned else 'taxonomy_source_unavailable'
        if not assigned or self.taxonomy_path is None:
            return result
        path = self.taxonomy_path
        if str(path) not in self.sources:
            self.sources[str(path)] = source_digest(path)
        result.update(mmseqs2_taxonomy_source=str(path), mmseqs2_taxonomy_source_sha256=self.sources[str(path)])
        with closing(sqlite3.connect(path.as_uri() + '?mode=ro', uri=True)) as conn:
            node = conn.execute('SELECT spname,rank,track FROM species WHERE taxid=?', (assigned,)).fetchone()
            if node is None:
                merged = conn.execute('SELECT taxid_new FROM merged WHERE taxid_old=?', (assigned,)).fetchone()
                node = conn.execute('SELECT spname,rank,track FROM species WHERE taxid=?', (merged[0],)).fetchone() if merged else None
            if node is None:
                result['mmseqs2_lineage_status'] = 'lca_taxid_unresolved_in_existing_database'
                return result
            expected = tuple(reversed([int(value) for value in node[2].split(',')]))
            lineage = expected
            missing = False
            if saved_lineage:
                canonical = []
                for ancestor in saved_lineage:
                    present = conn.execute('SELECT taxid FROM species WHERE taxid=?', (ancestor,)).fetchone()
                    if present is None:
                        merged = conn.execute('SELECT taxid_new FROM merged WHERE taxid_old=?', (ancestor,)).fetchone()
                        present = conn.execute('SELECT taxid FROM species WHERE taxid=?', (merged[0],)).fetchone() if merged else None
                    if present is None:
                        missing = True
                    else:
                        canonical.append(present[0])
                # Never borrow ranks from a foreign lineage ending in the same LCA ID.
                positions = [expected.index(ancestor) for ancestor in canonical if ancestor in expected]
                if (len(positions) != len(canonical) or positions != sorted(set(positions))
                        or not canonical or canonical[-1] != expected[-1]):
                    result['mmseqs2_lineage_status'] = 'saved_lineage_conflicts_existing_database'
                    return result
                lineage = tuple(canonical)
            for ancestor in lineage:
                row = conn.execute('SELECT spname,rank FROM species WHERE taxid=?', (ancestor,)).fetchone()
                if row is None:
                    missing = True
                    continue
                if row[1] in RANKS:
                    field = 'mmseqs2_' + row[1]
                    if result[field] and result[field] != row[0]:
                        raise ValueError('Conflicting ranks in saved MMseqs2 lineage')
                    result[field] = row[0]
            result['mmseqs2_lineage_status'] = ('saved_lineage_partially_unresolved' if missing else
                                               'saved_lineage_resolved' if saved_lineage else 'existing_database_lineage')
        return result

    def sequence_taxonomy(self, gene, species):
        """Saved query classification and saved per-gene host labels, never a best hit."""
        result = dict.fromkeys(SEQUENCE_FIELDS, '')
        result.update(mmseqs2_status='source_unavailable', mmseqs2_host_label_status='source_unavailable')
        if not species:
            return result
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', species) or species in {'.', '..'}:
            raise ValueError('Unsafe sequence taxonomy species identifier')
        if species not in self.sequence_cache:
            assignments, labels, metadata = {}, {}, {}
            path = self.sequence_root / (species + '_mmseqs2taxonomy.tsv') if self.sequence_root else None
            if path and path.is_file():
                raw = path.read_bytes()
                sha = hashlib.sha256(raw).hexdigest()
                self.sources[str(path.resolve())] = sha
                metadata.update(mmseqs2_source=str(path.resolve()), mmseqs2_source_sha256=sha)
                for row in csv.reader(io.StringIO(raw.decode('utf-8-sig')), delimiter='\t'):
                    if len(row) < 4 or not all(available(v) for v in row[:4]):
                        raise ValueError('Malformed saved MMseqs2 classification row')
                    if row[0] in assignments:
                        raise ValueError('Duplicate MMseqs2 query gene ID')
                    assigned = 0 if row[1] in {'0', '0.0'} else taxid(row[1])
                    saved = tuple(taxid(v) for v in row[8].split(';')) if len(row) >= 9 and available(row[8]) else ()
                    if saved and (None in saved or len(saved) != len(set(saved)) or saved[-1] != assigned):
                        raise ValueError('Saved MMseqs2 lineage disagrees with query LCA')
                    assignments[row[0]] = dict(mmseqs2_lca_taxid=str(assigned), mmseqs2_lca_rank=row[2],
                                               mmseqs2_lca_name=row[3],
                                               mmseqs2_lineage_taxids=';'.join(map(str, saved)),
                                               mmseqs2_status='assigned' if assigned else 'unclassified')
            path = self.host_root / (species + '_gene_taxonomy.tsv') if self.host_root else None
            if path and path.is_file():
                import pandas
                from scaffold_taxonomy import validate_gene_table

                raw = path.read_bytes()
                sha = hashlib.sha256(raw).hexdigest()
                self.sources[str(path.resolve())] = sha
                metadata.update(mmseqs2_host_label_source=str(path.resolve()), mmseqs2_host_label_source_sha256=sha)
                data = pandas.read_csv(io.StringIO(raw.decode('utf-8-sig')), sep='\t', dtype=str, keep_default_na=False)
                validate_gene_table(data)
                if set(data.species) - {species}:
                    raise ValueError('Scaffold taxonomy file contains a different species')
                for row in data.to_dict('records'):
                    labels.setdefault(row['gene_id'], {})[row['rank']] = row
            self.sequence_cache[species] = assignments, labels, metadata
        assignments, labels, metadata = self.sequence_cache[species]
        result.update(metadata)
        if metadata.get('mmseqs2_source'):
            result['mmseqs2_status'] = 'gene_record_unavailable'
        result.update(assignments.get(gene, {}))
        if gene in assignments:
            assignment = assignments[gene]
            saved = tuple(int(v) for v in assignment['mmseqs2_lineage_taxids'].split(';')) if assignment['mmseqs2_lineage_taxids'] else ()
            result.update(self.query_lineage(int(assignment['mmseqs2_lca_taxid']), saved))
        if metadata.get('mmseqs2_host_label_source'):
            result['mmseqs2_host_label_status'] = 'gene_record_unavailable'
        if gene in labels:
            result.update(mmseqs2_host_label_status='measured',
                          mmseqs2_host_class_label=labels[gene]['class']['label'],
                          mmseqs2_host_species_label=labels[gene]['species']['label'],
                          mmseqs2_host_species_taxid=labels[gene]['species']['host_taxid'])
        return result

    def family_leaves(self, family):
        """Read a neighbor's own existing family, including archived store members."""
        if family not in self.leaf_cache:
            name = family + '_stat.branch.tsv'
            try:
                with self.store.open_binary('stat_branch', name) as handle:
                    raw = handle.read()
            except FileNotFoundError:
                self.leaf_cache[family] = None
            else:
                self.family_sources['stat_branch/' + name] = hashlib.sha256(raw).hexdigest()
                reader = csv.DictReader(io.StringIO(raw.decode('utf-8-sig')), delimiter='\t')
                fields = reader.fieldnames or []
                if len(set(fields)) != len(fields) or not {'node_name', 'child1', 'child2'} <= set(fields):
                    raise ValueError('Missing or duplicate neighbor family columns')
                rows = list(reader)
                if any(None in r or any(v is None for v in r.values()) for r in rows):
                    raise ValueError('Malformed neighbor family row')
                leaves = {r['node_name']: r for r in rows if r['child1'] == r['child2'] == '-999'}
                if len(leaves) != sum(r['child1'] == r['child2'] == '-999' for r in rows):
                    raise ValueError('Duplicate family leaf in context annotations')
                self.leaf_cache[family] = leaves
        return self.leaf_cache[family]

    def get(self, gene, family='', leaf=None, species=''):
        row = self.rows.get(gene)
        validation = 'supplemental_only' if row is not None else 'annotation_unavailable'
        if row is not None and self.store is not None and leaf is None:
            family = family or available(row['orthogroup'])
            leaves = self.family_leaves(family)
            if leaves is None:
                validation = 'family_source_unavailable'
            elif gene not in leaves:
                raise ValueError('Context annotation gene is absent from its own family: ' + gene)
            else:
                leaf = leaves[gene]
        if leaf is not None:
            if leaf.get('node_name') != gene or leaf.get('child1') != leaf.get('child2') or leaf.get('child1') != '-999':
                raise ValueError('Context annotation requires the exact family leaf: ' + gene)
            validation = 'exact_family_leaf_verified'
        if row is not None:
            if family and available(row['orthogroup']) != family:
                raise ValueError('Context annotation family/gene mapping disagrees: ' + gene)
            if leaf and 'sprot_best' in leaf and available(row['besthit_accession']) != available(leaf['sprot_best']):
                raise ValueError('Context annotation best hit disagrees with the exact family leaf: ' + gene)
            if leaf and available(leaf.get('organism')) and available(row['besthit_organism']) \
                    and available(row['besthit_organism']) != available(leaf['organism']):
                raise ValueError('Context annotation hit organism disagrees with the exact family leaf: ' + gene)
            if leaf and available(leaf.get('taxid_y')) and available(row.get('besthit_taxid')) \
                    and taxid(row['besthit_taxid']) != taxid(leaf['taxid_y']):
                raise ValueError('Context annotation hit taxid disagrees with the exact family leaf: ' + gene)
            if leaf and available(leaf.get('sprot_recname')) and available(row.get('swissprot_best_hit_protein_name')) \
                    and available(row['swissprot_best_hit_protein_name']) != available(leaf['sprot_recname']):
                raise ValueError('Context annotation hit protein name disagrees with the exact family leaf: ' + gene)
            result = {k: available(row.get(k)) for k in FIELDS}
            if leaf and result['besthit_accession']:
                for field, native in [('besthit_organism', 'organism'), ('besthit_taxid', 'taxid_y'),
                                      ('swissprot_best_hit_protein_name', 'sprot_recname')]:
                    result[field] = result[field] or available(leaf.get(native))
            result['annotation_validation_status'] = validation
            result.update(self.sequence_taxonomy(gene, species))
            return result
        result = dict.fromkeys(FIELDS, '')
        result.update(gene_id=gene, orthogroup=family, protein_product_status='annotation_unavailable',
                      besthit_status='annotation_unavailable', annotation_validation_status=validation)
        # The family's exact saved leaf row is useful even without a supplemental table.
        if leaf and leaf.get('node_name') == gene and leaf.get('child1') == leaf.get('child2') == '-999':
            result.update(besthit_accession=available(leaf.get('sprot_best')),
                          besthit_organism=available(leaf.get('organism')),
                          besthit_taxid=available(leaf.get('taxid_y')),
                          swissprot_best_hit_protein_name=available(leaf.get('sprot_recname'))
                          if available(leaf.get('sprot_best')) else '',
                          besthit_source='existing exact family stat.branch.tsv leaf',
                          besthit_status='existing_family_leaf_hit' if available(leaf.get('sprot_best')) else 'no_existing_hit')
        result.update(self.sequence_taxonomy(gene, species))
        return result

    def verify(self, records=None):
        """Verify one page's exact dependencies, or every accumulated input.

        Returned annotation records retain source paths even for cached query
        classifications and lineages. The canonical supplemental table is always
        checked; its per-row upstream provenance is not a substitute for it.
        Exporters must also call the unrestricted check before completing.
        """
        paths, families = None, None
        if records is not None:
            paths = {self._annotation_source} if self._annotation_source else set()
            families = set()
            for row in records:
                paths.update(value for field, value in row.items()
                             if field.endswith('_source') and value in self.sources)
                family = available(row.get('orthogroup'))
                if family:
                    families.add('stat_branch/' + family + '_stat.branch.tsv')
        for path, expected in self.sources.items():
            if paths is not None and path not in paths:
                continue
            if source_digest(path) != expected:
                raise ValueError('Context annotation input changed during rendering')
        for logical, expected in self.family_sources.items():
            if families is not None and logical not in families:
                continue
            with self.store.open_binary(*logical.split('/', 1)) as handle:
                if hashlib.sha256(handle.read()).hexdigest() != expected:
                    raise ValueError('Neighbor family input changed during rendering')


@lru_cache(maxsize=4096)
def text_width(text):
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextToPath

    return TextToPath().get_text_width_height_descent(text, FontProperties(size=8), ismath=False)[0]


def wrap_cell(text, width):
    """Wrap against actual font metrics so long products/ranks cannot cross columns."""
    lines = []
    for paragraph in text.split('\n'):
        line = ''
        for word in paragraph.split():
            trial = (line + ' ' + word).strip()
            if text_width(trial) <= width:
                line = trial
                continue
            if line:
                lines.append(line)
                line = ''
            while text_width(word) > width:
                n = 1
                while n < len(word) and text_width(word[:n+1]) <= width:
                    n += 1
                lines.append(word[:n])
                word = word[n:]
            line = word
        lines.append(line)
    return '\n'.join(lines)


def annotation_cells(row, species, label):
    gene = row['gene_id'].removeprefix(species + '_')
    accession, organism = available(row.get('besthit_accession')), available(row.get('besthit_organism'))
    product = available(row.get('swissprot_best_hit_protein_name'))
    product = product + '\n[best-hit prediction]' if product and accession else 'Unavailable'
    hit = (organism or 'Organism unavailable') + '\n[' + accession + ']' if accession else 'Unavailable'
    taxonomy = []
    for i in range(0, len(RANKS), 2):
        taxonomy.append(' | '.join(f'{rank.capitalize()}: {available(row.get("besthit_" + rank)) or "unavailable"}'
                                   for rank in RANKS[i:i+2]))
    assigned = available(row.get('mmseqs2_lca_name'))
    classification = (assigned + '\nLCA rank: ' + (available(row.get('mmseqs2_lca_rank')) or 'unavailable')
                      + ' | taxid: ' + (available(row.get('mmseqs2_lca_taxid')) or 'unavailable')) if assigned else 'Unavailable'
    classification += '\nHost class: ' + (available(row.get('mmseqs2_host_class_label')) or 'unavailable')
    classification += '\nHost species: ' + (available(row.get('mmseqs2_host_species_label')) or 'unavailable')
    lineage_status = available(row.get('mmseqs2_lineage_status'))
    missing_rank = 'unresolved' if lineage_status and 'unavailable' not in lineage_status else 'unavailable'
    for i in range(0, len(RANKS), 2):
        classification += '\n' + ' | '.join(f'{rank.capitalize()}: {available(row.get("mmseqs2_" + rank)) or missing_rank}'
                                               for rank in RANKS[i:i+2])
    cells = [label + '\n' + gene, product, hit + '\n' + '\n'.join(taxonomy), classification]
    widths = [(b-a)*TABLE_WIDTH_PT-10 for a, b in zip(TABLE_EDGES[:-1], TABLE_EDGES[1:], strict=True)]
    return [wrap_cell(cell, width) for cell, width in zip(cells, widths, strict=True)]


def context_annotation_rows(entry, annotations, family, leaves):
    """Keep the genomic track's exact left-to-right neighbor numbering."""
    focal_id, species = entry['link']['gene_id'], entry['link'].get('gene_species', '')
    neighbors = sorted(entry['neighbors'], key=lambda r: int(r['start']))
    if not neighbors:
        neighbors = [dict(gene_id=focal_id)]
    result, number = [], 0
    from focus_hgt_context import neighbor_relation
    for gene in neighbors:
        focal = gene['gene_id'] == focal_id
        if not focal:
            number += 1
        label = 'Focal' if focal else str(number)
        leaf = leaves.get(gene['gene_id'])
        annotation = annotations.get(gene['gene_id'], family if focal or leaf is not None else '', leaf, species=species)
        relation = 'focal' if focal else neighbor_relation(gene, entry.get('focal'))
        cells = annotation_cells(annotation, species, label if focal else label + ' (' + relation + ')')
        result.append(dict(annotation, context_focal_gene_id=focal_id, side=entry['side'],
                           context_role='focal' if focal else 'neighbor', neighbor_label=label,
                           context_neighbor_relation=relation,
                           context_genomic_start_bp=gene.get('start', ''), context_genomic_end_bp=gene.get('end', ''),
                           context_genomic_strand=gene.get('strand', ''), context_scaffold=gene.get('chromosome', ''),
                           event_ids='; '.join(sorted(entry['event_ids'])), cells=cells,
                           height_pt=max(cell.count('\n')+1 for cell in cells)*10 + 9))
    return result


def draw_annotation_table(ax, records, color):
    from matplotlib.patches import Rectangle

    height = 25 + sum(r['height_pt'] for r in records)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, height)
    ax.axis('off')
    edges = TABLE_EDGES
    headings = ['Track label / gene ID', 'Protein product\n(Swiss-Prot best-hit)',
                'Swiss-Prot best hit / accession\nBest-hit taxonomic ranks',
                'MMseqs2 classification\nQuery taxonomic ranks / host match']
    ax.add_patch(Rectangle((0, height-25), 1, 25, color='#eef1f4', zorder=0))
    for x, text in zip(edges[:-1], headings, strict=True):
        ax.text(x+.006, height-8, text, fontsize=8, weight='bold', va='top')
    top = height-25
    for record in records:
        focal = record['context_role'] == 'focal'
        if focal:
            ax.add_patch(Rectangle((0, top-record['height_pt']), 1, record['height_pt'],
                                   color='#f3e5db' if record['side'] == 'recipient' else '#e8f0f6', zorder=0))
        for column, cell in enumerate(record['cells']):
            ax.text(edges[column]+.006, top-5, cell, fontsize=8, linespacing=1.22, va='top',
                    color=color if focal and column == 0 else '#333333',
                    weight='bold' if focal and column == 0 else 'normal')
        top -= record['height_pt']
        ax.plot([0, 1], [top, top], color='#d6dce1', lw=.5)
