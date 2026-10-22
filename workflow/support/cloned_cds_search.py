"""Bound miniprot indexing memory without cutting genomic contigs or changing seeds."""
import json
import tempfile
from pathlib import Path


def alignments(genome, queries, directory, code, cpus, batch_bases=250_000_000):
    from rescue_gene_models import run
    groups, current, bases = [], [], 0
    for seqid, length in zip(genome.references, genome.lengths, strict=True):
        if current and bases + length > batch_bases:
            groups.append(current)
            current, bases = [], 0
        current.append((seqid, length))
        bases += length
    if current:
        groups.append(current)
    directory = Path(directory)
    (directory / 'index_batches.json').write_text(json.dumps({
        'batch_bases': batch_bases, 'contigs_split': False,
        'batches': [dict(contigs=[r[0] for r in group], bases=sum(r[1] for r in group)) for group in groups]
    }, sort_keys=True) + '\n')
    with tempfile.TemporaryDirectory(prefix='.index-', dir=directory) as temporary:
        part = Path(temporary) / 'genome.fa'
        for index, group in enumerate(groups, 1):
            if len(groups) == 1:
                source = genome.filename.decode()
            else:
                with part.open('w') as handle:
                    for seqid, length in group:
                        handle.write('>' + seqid + '\n')
                        for start in range(0, length, 1_000_000):
                            handle.write(genome.fetch(seqid, start, min(length, start + 1_000_000)) + '\n')
                source = str(part)
            label = 'cloned_cds_miniprot' if len(groups) == 1 else f'cloned_cds_miniprot.{index}'
            target = directory / ('models.gff' if len(groups) == 1 else f'models.{index}.gff')
            run(['miniprot', '-u', '-T', code, '-t', cpus, '-N', 100, '--outs=0', '--gff', source, queries],
                directory, label, target)
            yield target
