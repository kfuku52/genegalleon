"""Extract complete, stranded coding-path witnesses from genome-aligned RNA BAMs."""
import hashlib
from contextlib import ExitStack
from pathlib import Path

import pysam

try:
    from gene_model_catalog import indexed_genome
except ImportError:
    from .gene_model_catalog import indexed_genome


def discover(directory, species):
    """Accept species-named BAMs or per-species subdirectories."""
    directory = Path(directory).resolve()
    if not directory.is_dir():
        raise ValueError('RNA BAM directory does not exist: ' + str(directory))
    result = {}
    for name in species:
        paths = {directory / (name + '.bam'), *directory.glob(name + '.*.bam'),
                 *(directory / name).glob('*.bam')}
        entries = []
        for path in sorted(p.resolve() for p in paths if p.is_file()):
            indices = [Path(str(path) + ext) for ext in ('.csi', '.bai')] + [path.with_suffix('.bai')]
            index = next((str(p.resolve()) for p in indices if p.is_file()), None)
            entries.append({'bam': str(path), 'index': index})
        if entries:
            result[name] = entries
    return result


def verified_reference(bam, genome, checked):
    """Verify the coordinate dictionary and each available reference checksum."""
    expected = dict(zip(genome.references, genome.lengths, strict=True))
    rows = bam.header.to_dict().get('SQ', [])
    observed = {r['SN']: r['LN'] for r in rows}
    if len(rows) != len(observed) or expected != observed:
        raise ValueError('RNA BAM reference dictionary differs from refinement genome: ' + str(bam.filename))
    for row in rows:
        if 'M5' not in row:
            continue
        name = row['SN']
        if name not in checked:
            proof = hashlib.md5(usedforsecurity=False)
            for start in range(0, row['LN'], 1024 * 1024):
                proof.update(genome.fetch(name, start, min(row['LN'], start + 1024 * 1024)).upper().encode())
            checked[name] = proof.hexdigest()
        if checked[name] != row['M5'].lower():
            raise ValueError('RNA BAM reference checksum differs for ' + name)


def splice_strand(genome, seqid, blocks):
    plus, minus = True, True
    for (_, left), (right, _) in zip(blocks, blocks[1:], strict=False):
        ends = (genome.fetch(seqid, left, left + 2).upper(), genome.fetch(seqid, right - 2, right).upper())
        plus &= ends in {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')}
        minus &= ends in {('CT', 'AC'), ('CT', 'GC'), ('GT', 'AT')}
    return '+' if plus and not minus else '-' if minus and not plus else None


def coding_blocks(read, start, end):
    blocks = []
    for a, b in read.get_blocks():
        a, b = max(a, start), min(b, end)
        if a >= b:
            continue
        if blocks and a == blocks[-1][1]:
            blocks[-1][1] = b
        else:
            blocks.append([a, b])
    return blocks


def evidence_rows(species, genome_path, entries, candidates, scratch):
    """Each reported path is spanned by one primary read on the same strand."""
    if not entries:
        return []
    paths = {(c['seqid'], c['strand'], tuple(tuple(b[:2]) for b in c.get('blocks', c.get('cds', ()))))
             for c in candidates if c.get('blocks', c.get('cds')) and c['strand'] in {'+', '-'}}
    if not paths:
        return []
    result = []
    with ExitStack() as stack:
        genome = stack.enter_context(indexed_genome(genome_path))
        bams, checked = [], {}
        for i, item in enumerate(entries):
            index = item['index']
            if not index:
                index = str(scratch / f'rna_{i}.csi')
                pysam.index('-c', '-o', index, item['bam'])
            bam = stack.enter_context(pysam.AlignmentFile(item['bam'], 'rb', index_filename=index))
            verified_reference(bam, genome, checked)
            bams.append((bam, Path(item['bam']).name))
        for seqid, strand, blocks in sorted(paths):
            genomic = sorted(map(list, blocks))
            start, end = genomic[0][0], genomic[-1][1]
            inferred = splice_strand(genome, seqid, genomic)
            for bam, name in bams:
                for read in bam.fetch(seqid, start, end):
                    if (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate
                            or read.is_qcfail or read.mapping_quality < 20):
                        continue
                    observed = read.get_tag('XS') if read.has_tag('XS') else inferred
                    if observed != strand or coding_blocks(read, start, end) != genomic:
                        continue
                    result.append({'species': species, 'seqid': seqid, 'strand': strand,
                                   'cds_blocks': [list(b) for b in blocks],
                                   'transcript_id': name + ':' + read.query_name, 'count': 1})
                    break
    return result
