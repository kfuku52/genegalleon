"""Extract complete, stranded coding-path witnesses from genome-aligned RNA BAMs."""
import hashlib
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
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
    result, owners = {}, {}
    for name in species:
        paths = {directory / (name + '.bam'), *directory.glob(name + '.*.bam'),
                 *(directory / name).glob('*.bam')}
        entries = {}
        for path in sorted({p.resolve() for p in paths if p.is_file()}, key=str):
            stat = path.stat()
            identity = stat.st_dev, stat.st_ino
            if owners.setdefault(identity, name) != name:
                raise ValueError('RNA BAM is assigned to multiple species: ' + str(path))
            indices = [Path(str(path) + ext) for ext in ('.csi', '.bai')] + [path.with_suffix('.bai')]
            index = next((str(p.resolve()) for p in indices
                          if p.is_file() and p.stat().st_mtime_ns >= stat.st_mtime_ns), None)
            if identity not in entries or (index and not entries[identity]['index']):
                entries[identity] = {'bam': str(path), 'index': index}
        if entries:
            result[name] = sorted(entries.values(), key=lambda entry: Path(entry['bam']))
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
    """Clip aligned blocks to the CDS, accepting only CIGAR N as a splice."""
    blocks, position = [], read.reference_start
    for operation, length in read.cigartuples or ():
        if operation in (1, 4, 5, 6):
            continue
        if operation not in (0, 2, 3, 7, 8):
            return []
        stop = position + length
        a, b = max(position, start), min(stop, end)
        position = stop
        if a >= b or operation == 3:
            continue
        if operation == 2:
            return []
        if blocks and a == blocks[-1][1]:
            blocks[-1][1] = b
        else:
            blocks.append([a, b])
    return blocks


def overlapping_windows(paths):
    """Group overlapping spans and leave distinct loci as separate indexed queries."""
    contigs = defaultdict(list)
    for path in paths:
        seqid, _, blocks = path
        genomic = sorted(map(list, blocks))
        contigs[seqid].append((genomic[0][0], genomic[-1][1], path, genomic))
    for seqid, models in sorted(contigs.items()):
        window, end = [], -1
        for model in sorted(models):
            if model[0] >= end and window:
                yield seqid, window, end
                window = []
            window.append(model)
            end = max(end, model[1])
        if window:
            yield seqid, window, end


def bam_witnesses(bam, genome, windows):
    """Find the first complete spanning read per candidate with one fetch per window."""
    witnesses, inferred = {}, {}
    for seqid, models, window_end in windows:
        starts = [model[0] for model in models]
        min_span = min(end - start for start, end, _, _ in models)
        remaining = len(models)
        for read in bam.fetch(seqid, starts[0], window_end):
            if (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate
                    or read.is_qcfail or read.mapping_quality < 20):
                continue
            read_start, read_end = read.reference_start, read.reference_end
            if read_end is None or read_end - read_start < min_span:
                continue
            first = bisect_left(starts, read_start)
            last = bisect_right(starts, read_end - min_span)
            if first >= last:
                continue
            declared = read.get_tag('XS') if read.has_tag('XS') else None
            for i in range(first, last):
                start, end, path, genomic = models[i]
                if end > read_end or path in witnesses or coding_blocks(read, start, end) != genomic:
                    continue
                if declared in ('+', '-'):
                    observed = declared
                else:
                    if path not in inferred:
                        inferred[path] = splice_strand(genome, seqid, genomic)
                    observed = inferred[path]
                if observed == path[1]:
                    witnesses[path] = read.query_name
                    remaining -= 1
            if not remaining:
                break
    return witnesses


def evidence_rows(species, genome_path, entries, candidates, scratch):
    """Each reported path is spanned by one primary read on the same strand."""
    if not entries:
        return []
    paths = {(c['seqid'], c['strand'], tuple(tuple(b[:2]) for b in c.get('blocks', c.get('cds', ()))))
             for c in candidates if c.get('blocks', c.get('cds')) and c['strand'] in {'+', '-'}}
    if not paths:
        return []
    windows = list(overlapping_windows(paths))
    names = Counter(Path(item['bam']).name for item in entries)
    results = defaultdict(list)
    with indexed_genome(genome_path) as genome:
        checked = {}
        for i, item in enumerate(entries):
            index = item['index']
            if not index:
                index = str(scratch / f'rna_{i}.csi')
                pysam.index('-c', '-o', index, item['bam'])
            name = Path(item['bam']).name
            if names[name] > 1:
                name = item['bam']
            with pysam.AlignmentFile(item['bam'], 'rb', index_filename=index) as bam:
                verified_reference(bam, genome, checked)
                for path, read in bam_witnesses(bam, genome, windows).items():
                    seqid, strand, blocks = path
                    results[path].append({'species': species, 'seqid': seqid, 'strand': strand,
                                          'cds_blocks': [list(b) for b in blocks],
                                          'transcript_id': name + ':' + read, 'count': 1})
    return [row for path in sorted(results) for row in results[path]]
