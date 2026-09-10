"""Validate GFF feature order, including explicitly ordered trans-spliced CDS."""
import re
from urllib.parse import unquote


def has_trans_splicing_exception(text):
    return any('trans-splicing' in [unquote(x).strip() for x in field.split('=', 1)[1].split(',')]
               for field in str(text).split(';') if field.startswith('exception='))


def ordered_annotated_blocks(rows, gene_id):
    """Return (sequence, strand, start, end) blocks and their splice mode.

    Rows contain those four coordinate fields followed by GFF attributes.
    Trans-splicing is accepted only when every row declares the exception and
    every distinct block has an unambiguous, complete 1..N part order.
    """
    rows = list(rows)
    coordinates = [(str(seq), str(strand), int(start), int(end))
                   for seq, strand, start, end, _attrs in rows]
    attributes = []
    for *_coords, text in rows:
        fields = {}
        for field in str(text).split(';'):
            if '=' in field:
                key, value = field.split('=', 1)
                if key in fields and fields[key] != value:
                    raise ValueError(f'Conflicting GFF attribute {key} for {gene_id}')
                fields[key] = value
        attributes.append(fields)
    declared = ['trans-splicing' in [unquote(x).strip() for x in a.get('exception', '').split(',')]
                for a in attributes]
    if not any(declared):
        return ordered_feature_blocks(coordinates, gene_id), 'cis'
    if not all(declared):
        raise ValueError(f'Incomplete trans-splicing annotation for {gene_id}')
    parts = {}
    coordinate_parts = {}
    totals = set()
    for coordinate, attr in zip(coordinates, attributes, strict=True):
        match = re.fullmatch(r'([1-9][0-9]*)(?:/([1-9][0-9]*))?', attr.get('part', ''))
        if match is None:
            raise ValueError(f'Trans-splicing requires explicit part order for {gene_id}')
        part = int(match[1])
        if match[2]:
            totals.add(int(match[2]))
        if part in parts and parts[part] != coordinate:
            raise ValueError(f'Conflicting trans-splicing part for {gene_id}')
        if coordinate in coordinate_parts and coordinate_parts[coordinate] != part:
            raise ValueError(f'Repeated trans-splicing coordinate for {gene_id}')
        parts[part] = coordinate
        coordinate_parts[coordinate] = part
    if set(parts) != set(range(1, len(parts) + 1)) or (totals and totals != {len(parts)}):
        raise ValueError(f'Incomplete trans-splicing part order for {gene_id}')
    # Validate coordinate geometry within each genomic strand, without inventing
    # an order between different strands/contigs or overriding annotated parts.
    systems = {}
    for block in parts.values():
        systems.setdefault(block[:2], []).append(block)
    for blocks in systems.values():
        ordered_feature_blocks(blocks, gene_id)
    return [parts[i] for i in sorted(parts)], 'trans-splicing'


def ordered_feature_blocks(rows, gene_id):
    blocks = {(str(sequence), str(strand), int(start), int(end))
              for sequence, strand, start, end in rows}
    if len({(sequence, strand) for sequence, strand, _start, _end in blocks}) != 1:
        raise ValueError(f'Conflicting GFF coordinate systems for {gene_id}')
    if any(strand not in ('+', '-') or start < 1 or end < start
           for _sequence, strand, start, end in blocks):
        raise ValueError(f'Invalid GFF coordinates or strand for {gene_id}')
    blocks = sorted(blocks, key=lambda block: (block[2], block[3]))
    if any(right[2] <= left[3] for left, right in zip(blocks, blocks[1:], strict=False)):
        raise ValueError(f'Overlapping GFF feature blocks for {gene_id}')
    return blocks[::-1] if blocks[0][1] == '-' else blocks
