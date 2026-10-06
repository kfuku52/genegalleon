"""Per-gene existing annotations for focused genomic context pages."""

import csv
import hashlib
import io
from functools import lru_cache
from pathlib import Path

RANKS = ('kingdom', 'phylum', 'class', 'order', 'family', 'genus')
TABLE_EDGES = (0, .16, .41, .66, 1)
TABLE_WIDTH_PT = .435 * 22 * 72
FIELDS = ('gene_id', 'orthogroup', 'protein_product_name', 'protein_product_status',
          'protein_product_source', 'protein_product_source_sha256', 'protein_product_feature_ids',
          'protein_product_mapping_scope', 'gene_description', 'swissprot_best_hit_protein_name',
          'swissprot_name_source', 'swissprot_name_source_sha256',
          'besthit_accession', 'besthit_organism', 'besthit_taxid', 'besthit_source',
          'besthit_source_sha256', 'besthit_status', 'besthit_coverage_percent', 'besthit_identity_percent',
          'besthit_evalue', 'taxonomy_source', 'taxonomy_source_sha256',
          *(f'besthit_{rank}' for rank in RANKS))


def available(value):
    return str(value).strip() if value is not None and str(value).strip().lower() not in {'', 'na', 'nan', 'none'} else ''


class ContextAnnotations:
    """Exact gene IDs only; neighbor genes never inherit focal annotations."""

    def __init__(self, path=''):
        self.rows, self.sources, self.display_audit = {}, {}, []
        if path:
            path = Path(path).resolve()
            raw = path.read_bytes()
            self.sources[str(path)] = hashlib.sha256(raw).hexdigest()
            reader = csv.DictReader(io.StringIO(raw.decode()), delimiter='\t')
            required = {'gene_id', 'orthogroup', 'besthit_accession', 'besthit_organism',
                        *[f'besthit_{r}' for r in RANKS]}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError('Context annotations lack required per-gene columns')
            for row in reader:
                gene = available(row['gene_id'])
                if not gene or gene in self.rows:
                    raise ValueError('Duplicate or empty context annotation gene ID')
                if (available(row['besthit_organism']) or available(row.get('swissprot_best_hit_protein_name'))
                        or any(available(row[f'besthit_{r}']) for r in RANKS)) \
                        and not available(row['besthit_accession']):
                    raise ValueError('Best-hit organism/taxonomy requires the same hit accession')
                self.rows[gene] = row

    def get(self, gene, family='', leaf=None):
        row = self.rows.get(gene)
        if row is not None:
            if family and available(row['orthogroup']) != family:
                raise ValueError('Context annotation family/gene mapping disagrees: ' + gene)
            if leaf and available(leaf.get('sprot_best')) and available(row['besthit_accession']) != available(leaf['sprot_best']):
                raise ValueError('Context annotation best hit disagrees with the exact family leaf: ' + gene)
            if leaf and available(leaf.get('organism')) and available(row['besthit_organism']) \
                    and available(row['besthit_organism']) != available(leaf['organism']):
                raise ValueError('Context annotation hit organism disagrees with the exact family leaf: ' + gene)
            if leaf and available(leaf.get('taxid_y')) and available(row.get('besthit_taxid')) \
                    and int(float(row['besthit_taxid'])) != int(float(leaf['taxid_y'])):
                raise ValueError('Context annotation hit taxid disagrees with the exact family leaf: ' + gene)
            return {k: available(row.get(k)) for k in FIELDS}
        result = dict.fromkeys(FIELDS, '')
        result.update(gene_id=gene, orthogroup=family, protein_product_status='annotation_unavailable',
                      besthit_status='annotation_unavailable')
        # The family's exact saved leaf row is useful even without a supplemental table.
        if leaf and leaf.get('node_name') == gene and leaf.get('child1') == leaf.get('child2') == '-999':
            result.update(besthit_accession=available(leaf.get('sprot_best')),
                          besthit_organism=available(leaf.get('organism')),
                          besthit_taxid=available(leaf.get('taxid_y')),
                          swissprot_best_hit_protein_name=available(leaf.get('sprot_recname')),
                          besthit_source='existing exact family stat.branch.tsv leaf',
                          besthit_status='existing_family_leaf_hit' if available(leaf.get('sprot_best')) else 'no_existing_hit')
        return result

    def verify(self):
        for path, expected in self.sources.items():
            if hashlib.sha256(Path(path).read_bytes()).hexdigest() != expected:
                raise ValueError('Context annotation input changed during rendering')


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
    cells = [label + '\n' + gene, product, hit, '\n'.join(taxonomy)]
    widths = [(b-a)*TABLE_WIDTH_PT-10 for a, b in zip(TABLE_EDGES[:-1], TABLE_EDGES[1:], strict=True)]
    return [wrap_cell(cell, width) for cell, width in zip(cells, widths, strict=True)]


def context_annotation_rows(entry, annotations, family, leaves):
    """Keep the genomic track's exact left-to-right neighbor numbering."""
    focal_id, species = entry['link']['gene_id'], entry['link'].get('gene_species', '')
    neighbors = sorted(entry['neighbors'], key=lambda r: int(r['start']))
    if not neighbors:
        neighbors = [dict(gene_id=focal_id)]
    result, number = [], 0
    for gene in neighbors:
        focal = gene['gene_id'] == focal_id
        if not focal:
            number += 1
        label = 'Focal' if focal else str(number)
        annotation = annotations.get(gene['gene_id'], family if focal else '', leaves.get(gene['gene_id']))
        cells = annotation_cells(annotation, species, label)
        result.append(dict(annotation, context_focal_gene_id=focal_id, side=entry['side'],
                           context_role='focal' if focal else 'neighbor', neighbor_label=label,
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
    headings = ['Track label / gene ID', 'Protein product (best-hit)', 'Best-hit organism / accession', 'Best-hit taxonomic ranks']
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
