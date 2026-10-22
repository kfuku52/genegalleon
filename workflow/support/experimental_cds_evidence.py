"""Portable proofs for experimental CDS projected onto an unchanged reference.

M is an exact nucleotide match, X a substitution, I experimental-only bases,
and D reference-only bases. Coordinates remain reference coordinates; the
experimental CDS is never asserted to be a genome-identical reconstruction.
"""
import hashlib
import json
import re

PREFIX = 'gg_projection_'
FIELDS = ('version', 'clone_sha256', 'reference_sha256', 'blocks_sha256', 'cigar', 'partial_5prime', 'partial_3prime')


def sha(sequence):
    return hashlib.sha256(sequence.encode()).hexdigest()


def blocks_digest(blocks, strand, seqid):
    return sha(json.dumps({'seqid': seqid, 'strand': strand, 'blocks': sorted([list(b) for b in blocks], reverse=strand == '-')},
                          sort_keys=True, separators=(',', ':')))


def metadata(rows):
    found = []
    for row in rows:
        attr = row['attributes']
        if any(key.startswith(PREFIX) for key in attr):
            if attr.get('gg_sequence_origin') != 'experimental_clone':
                raise ValueError('Experimental projection lacks its explicit sequence origin')
            value = {key: attr.get(PREFIX + key, '') for key in FIELDS}
            if value not in found:
                found.append(value)
    if len(found) > 1:
        raise ValueError('Conflicting experimental CDS projection proofs')
    return found[0] if found else None


def alignment(reference, clone):
    from Bio.Align import PairwiseAligner
    aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-3,
                             open_gap_score=-8, extend_gap_score=-1)
    coordinates = aligner.align(reference, clone)[0].coordinates.T.tolist()
    runs, differences = [], []
    def append(op, size):
        if runs and runs[-1][1] == op:
            runs[-1][0] += size
        else:
            runs.append([size, op])
    for (r0, c0), (r1, c1) in zip(coordinates, coordinates[1:], strict=False):
        if r1 > r0 and c1 > c0:
            if r1 - r0 != c1 - c0:
                raise ValueError('Invalid nucleotide alignment coordinates')
            for offset, (a, b) in enumerate(zip(reference[r0:r1], clone[c0:c1], strict=True)):
                append('M' if a == b else 'X', 1)
                if a != b:
                    differences.append(dict(kind='substitution', reference_offset=r0 + offset,
                                            clone_offset=c0 + offset, reference=a, clone=b))
        else:
            append('D' if r1 > r0 else 'I', max(r1 - r0, c1 - c0))
            differences.append(dict(kind='deletion' if r1 > r0 else 'insertion', reference_offset=r0,
                                    clone_offset=c0, reference=reference[r0:r1], clone=clone[c0:c1]))
    return ''.join(str(size) + op for size, op in runs), differences


def proof(reference, clone, blocks, strand, seqid, *, partial_5prime=False, partial_3prime=False):
    cigar, differences = alignment(reference, clone)
    value = dict(version='1', clone_sha256=sha(clone), reference_sha256=sha(reference),
                 blocks_sha256=blocks_digest(blocks, strand, seqid), cigar=cigar,
                 partial_5prime=str(bool(partial_5prime)).lower(), partial_3prime=str(bool(partial_3prime)).lower())
    verify(reference, clone, blocks, strand, value, seqid)
    return value, differences


def verify(reference, clone, blocks, strand, value, seqid):
    if (set(value) != set(FIELDS) or value['version'] != '1'
            or value['clone_sha256'] != sha(clone) or value['reference_sha256'] != sha(reference)
            or value['blocks_sha256'] != blocks_digest(blocks, strand, seqid)
            or any(value[k] not in {'true', 'false'} for k in ('partial_5prime', 'partial_3prime'))):
        raise ValueError('Experimental CDS/reference/coordinate proof changed')
    cigar = value['cigar']
    if not re.fullmatch(r'(?:[1-9]\d*[MXID])+', cigar):
        raise ValueError('Invalid experimental CDS difference CIGAR')
    r = c = 0
    for number, op in re.findall(r'(\d+)([MXID])', cigar):
        size = int(number)
        a, b = reference[r:r + size], clone[c:c + size]
        if op in 'MX':
            if len(a) != size or len(b) != size or any((x == y) != (op == 'M') for x, y in zip(a, b, strict=True)):
                raise ValueError('Experimental CDS difference CIGAR disagrees with sequence')
            r += size
            c += size
        elif op == 'D':
            r += size
        else:
            c += size
        if r > len(reference) or c > len(clone):
            raise ValueError('Experimental CDS difference CIGAR exceeds sequence')
    if (r, c) != (len(reference), len(clone)):
        raise ValueError('Experimental CDS difference CIGAR does not cover both sequences')
    return True


def validate_pair(cds, gff, genome, species):
    """Reverify exported proofs against actual CDS and reference coordinates."""
    import contextlib
    from pathlib import Path

    from Bio.Seq import Seq
    from fasta_sequence_store import open_text
    from format_species_annotation.common import parse_gff_attributes
    from format_species_annotation.genome_intervals import genome_intervals
    projected, declared, nodes = {}, set(), {}
    with open_text(Path(gff)) as stream:
        for line in stream:
            if line.strip() == '##FASTA':
                break
            if line.startswith('#') or not line.strip():
                continue
            if PREFIX not in line and 'gg_sequence_origin' not in line:
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < 9:
                raise ValueError('Invalid GFF row in experimental projection validation')
            attr = {key: ','.join(values) for key, values in parse_gff_attributes(fields[8]).items()}
            if not any(key.startswith(PREFIX) for key in attr) and attr.get('gg_sequence_origin') != 'experimental_clone':
                continue
            row = dict(attributes=attr, start=int(fields[3]) - 1, end=int(fields[4]),
                       seqid=fields[0], strand=fields[6], phase=fields[7])
            value = metadata([row])
            if value is None:
                raise ValueError('Experimental model is missing its projection proof')
            if fields[2] == 'gene' or fields[2].lower().endswith(('rna', 'transcript')):
                if attr.get('ID') in nodes or not attr.get('ID'):
                    raise ValueError('Experimental gene/transcript has missing or duplicate ID')
                nodes[attr['ID']] = (value, attr.get('Parent', ''))
            if fields[2].lower().endswith(('rna', 'transcript')):
                declared.add(attr['ID'])
            if fields[2] == 'CDS':
                parents = attr.get('Parent', '').split(',')
                if len(parents) != 1 or not parents[0]:
                    raise ValueError('Experimental projection CDS must have one transcript owner')
                projected.setdefault(parents[0], []).append(row)
    if declared != projected.keys():
        raise ValueError('Experimental transcript/projection CDS ownership changed')
    for transcript, rows in projected.items():
        value, parent = nodes[transcript]
        if parent not in nodes or nodes[parent][0] != value or metadata(rows) != value:
            raise ValueError('Experimental gene/transcript/CDS projection proofs disagree')
    if set(nodes) != declared | {nodes[t][1] for t in declared}:
        raise ValueError('Experimental projection ownership changed')
    if not projected:
        if nodes:
            raise ValueError('Experimental gene is missing its projected transcript')
        return 0
    if genome is None:
        raise ValueError('Experimental CDS projection validation requires the unchanged reference genome')
    import format_species_inputs as fsi
    task = fsi.prepare_cds_identifier_task(dict(species_prefix=species, species_key=species, provider='local',
        cds_path=Path(cds), gff_path=Path(gff), genome_path=Path(genome), gene_grouping_mode='strict', gff_repair_mode='safe'))
    grouping = task['_gff_cds_grouping_index']
    by_gene = {}
    for transcript in projected:
        token = grouping['transcript_gene_tokens'].get(transcript)
        if not token or token in by_gene:
            raise ValueError('Experimental gene has missing or multiple projected transcripts')
        by_gene[token] = transcript
    sequences = {}
    for header, sequence in fsi.iter_task_cds_records(task):
        identifier = fsi.build_formatted_cds_id(task, header)
        _, match = fsi.resolve_gene_aggregate_id(task, header, identifier)
        transcript = by_gene.get(match.get('gene_token'))
        if transcript:
            if transcript in sequences:
                raise ValueError('Experimental gene has duplicate CDS records')
            sequences[transcript] = ''.join(sequence.split()).upper()
    if sequences.keys() != projected.keys():
        raise ValueError('Experimental transcript is missing from the CDS FASTA')
    regions = [(r['seqid'], r['start'], r['end']) for rows in projected.values() for r in rows]
    with contextlib.ExitStack() as stack:
        reader = stack.enter_context(genome_intervals(Path(genome), regions))
        for transcript, rows in projected.items():
            if len({(r['seqid'], r['strand']) for r in rows}) != 1:
                raise ValueError('Experimental CDS projection spans inconsistent contig/strand')
            strand = rows[0]['strand']
            blocks = [[r['start'], r['end'], int(r['phase'])] for r in rows]
            pieces = [reader.fetch(rows[0]['seqid'], a, b) for a, b, _ in sorted(blocks, reverse=strand == '-')]
            reference = ''.join(str(Seq(p).reverse_complement()) if strand == '-' else p for p in pieces)
            verify(reference, sequences[transcript], blocks, strand, metadata(rows), rows[0]['seqid'])
    return len(projected)
