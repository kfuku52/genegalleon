#!/usr/bin/env python3
"""Measure live GFF transcript selection, phase validation and feature summaries."""
import argparse
import contextlib
import hashlib
import io
import json
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def frame_fingerprint(frame):
    schema = [(str(column), str(dtype)) for column, dtype in frame.dtypes.items()]
    return hashlib.sha256(json.dumps(schema).encode() + frame.to_csv(sep='\t', index=False).encode()).hexdigest()


def worker(args):
    sys.path.insert(0, str(args.support_root))
    import gff2genestat as module
    import pandas

    rows = []
    genes = []
    for index in range(args.genes):
        gene = f'gene{index:07d}'
        genes.append('Species_a_' + gene)
        span = (args.exons - 1) * 200 + 110
        start = index * max(1000, span + 100)
        strand = '+' if index % 2 else '-'
        rows.append(('chr1', 'fixture', 'gene', start+1, start+span, '.', strand, '.', f'ID={gene}'))
        for isoform in range(2 if index % 10 == 0 else 1):
            transcript = f'{gene}.t{isoform}'
            rows.append(('chr1', 'fixture', 'mRNA', start+1, start+span, '.', strand, '.',
                         f'ID={transcript};Parent={gene}'))
            for block in reversed(range(args.exons - 1 if isoform else args.exons)):
                position = start + block * 200 + 11
                row = ('chr1', 'fixture', 'CDS', position, position+89, '.', strand, '0', f'Parent={transcript}')
                rows.append(row)
                if block == 0 and index % 25 == 0:
                    rows.append(row)  # Duplicate coordinate/phase records are still checked.
            position = start + (1 if strand == '+' else span - 9)
            rows.append(('chr1', 'fixture', 'UTR', position, position+9, '.', strand, '.', f'Parent={transcript}'))
    frame = pandas.DataFrame(rows, columns=['sequence', 'source', 'feature', 'start', 'end',
                                          'score', 'strand', 'phase', 'attributes'])
    names = pandas.Series(genes)
    columns = ['gene_id', 'feature_size', 'num_intron', 'intron_positions', 'chromosome', 'start', 'end',
               'strand', 'feature_blocks', 'feature_type', 'gff_transcript_id', 'utr_status', 'utr_blocks',
               'cds_first_phase', 'phase_status', 'cds_partial', 'splice_mode']
    warnings = io.StringIO()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(warnings):
        begin = time.perf_counter()
        selected = module.extract_by_ids(frame, names, 'CDS', 'longest')
        selection_done = time.perf_counter()
        annotated = module.attach_transcript_structure(selected, frame)
        structure_done = time.perf_counter()
        result = module.summarize_gene_features(annotated, columns)
        finished = time.perf_counter()
    assert len(result) == args.genes
    assert set(result['feature_size']) == {90 * args.exons}
    assert set(result['phase_status']) == {'consistent'}
    assert set(result['num_intron']) == {args.exons - 1}
    print(json.dumps({'seconds': finished-begin, 'selection_seconds': selection_done-begin,
        'structure_seconds': structure_done-selection_done, 'summary_seconds': finished-structure_done,
        'selected_sha256': frame_fingerprint(selected), 'annotated_sha256': frame_fingerprint(annotated),
        'result_sha256': frame_fingerprint(result), 'warnings_sha256': hashlib.sha256(warnings.getvalue().encode()).hexdigest(),
        'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'pandas': pandas.__version__}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--genes', type=int, default=2000)
    parser.add_argument('--exons', type=int, default=3)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', type=Path)
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 1 <= args.genes <= 10000:
        parser.error('--genes must be 1..10000')
    if not 2 <= args.exons <= 1024 or args.genes * args.exons > 200000:
        parser.error('--exons must be 2..1024, with at most 200000 total gene/exon blocks')
    if args.worker:
        worker(args)
        return
    if args.output is None:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'genes': args.genes, 'exons': args.exons,
               'python': sys.version, 'platform': platform.platform(), 'samples': []}
    with tempfile.TemporaryDirectory(prefix='gg-gff-benchmark-') as temporary:
        for trial in range(4):
            root = Path(temporary) / str(trial)
            root.mkdir()
            output = subprocess.check_output([sys.executable, str(Path(__file__).resolve()),
                '--worker', str(root), '--support-root', str(args.support_root),
                '--genes', str(args.genes), '--exons', str(args.exons)], text=True)
            sample = json.loads(output)
            if trial:
                results['samples'].append(sample)
            print('warmup' if trial == 0 else trial, sample['seconds'], flush=True)
    for key in ('selected_sha256', 'annotated_sha256', 'result_sha256', 'warnings_sha256'):
        assert len({row[key] for row in results['samples']}) == 1
    for key in ('seconds', 'selection_seconds', 'structure_seconds', 'summary_seconds', 'peak_rss_kib_linux'):
        results['median_'+key] = statistics.median(row[key] for row in results['samples'])
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    print('median', results['median_seconds'], results['median_peak_rss_kib_linux'], flush=True)


if __name__ == '__main__':
    main()
