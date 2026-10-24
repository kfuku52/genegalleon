"""Measure CDS phase/UTR checks with identical selected GFF rows in one runtime."""
import argparse
import ast
import cProfile
import hashlib
import json
import multiprocessing
import platform
import pstats
import resource
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
import gff2genestat as gff
import pandas


def fixture(genes, blocks):
    rows, utr = [], []
    for index in range(genes):
        gene, transcript, base = f'g{index}', f't{index}', index * (blocks * 30 + 100)
        strand = '-' if index % 2 else '+'
        for block in range(blocks):
            start = base + 10 + block * 30
            rows.append((gene, transcript, 'chr1', strand, start, start + 8,
                         '.' if index % 7 == 0 else '0', f'Parent={transcript}', 'CDS'))
        start = base + 1 if strand == '+' else base + 19 + (blocks - 1) * 30
        utr.append((gene, transcript, 'chr1', strand, start, start + 8, '.', f'Parent={transcript}', 'UTR'))
    columns = ['gene_id', 'selected_transcript', 'sequence', 'strand', 'start', 'end', 'phase', 'attributes', 'feature']
    selected = pandas.DataFrame(rows, columns=columns)
    ancestors = []
    for index in range(genes):
        gene, transcript, base = f'g{index}', f't{index}', index * (blocks * 30 + 100)
        strand = '-' if index % 2 else '+'
        ancestors.extend([(gene, transcript, 'chr1', strand, base + 1, base + blocks * 30 + 30, '.',
                           f'ID={gene}', 'gene'),
                          (gene, transcript, 'chr1', strand, base + 1, base + blocks * 30 + 30, '.',
                           f'ID={transcript};Parent={gene}', 'mRNA')])
    return selected, pandas.concat([selected, pandas.DataFrame(utr + ancestors, columns=columns)], ignore_index=True)


def measure(source, frames, profile, pipeline, connection):
    namespace = dict(vars(gff))
    function = next(node for node in ast.parse(source).body
                    if isinstance(node, ast.FunctionDef) and node.name == 'attach_transcript_structure')
    exec(compile(ast.Module(body=[function], type_ignores=[]), '<measured-structure>', 'exec'), namespace)
    timer = cProfile.Profile() if profile else None
    if timer:
        timer.enable()
    started, cpu = time.perf_counter(), time.process_time()
    if pipeline:
        identifiers = pandas.Series(['Plant_species_' + value for value in frames[0].gene_id.unique()])
        selected = gff.extract_by_ids(frames[1], identifiers, 'CDS', 'longest')
        annotated = namespace['attach_transcript_structure'](selected, frames[1])
        result = gff.summarize_gene_features(annotated, [
            'gene_id', 'feature_size', 'num_intron', 'intron_positions', 'chromosome', 'start', 'end',
            'strand', 'feature_blocks', 'feature_type', 'gff_transcript_id', 'utr_blocks', 'cds_first_phase',
            'splice_mode', 'feature_block_sequences', 'feature_block_strands', 'transcript_junction_positions',
            'phase_status', 'cds_partial', 'utr_status', 'structure_status'])
    else:
        result = namespace['attach_transcript_structure'](*frames)
    seconds, cpu_seconds = time.perf_counter() - started, time.process_time() - cpu
    if timer:
        timer.disable()
    signature = hashlib.sha256(pandas.util.hash_pandas_object(result, index=True).values.tobytes()).hexdigest()
    value = dict(seconds=seconds, cpu_seconds=cpu_seconds,
                 peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                 checksum=signature, rows=len(result), dtypes=list(map(str, result.dtypes)))
    if timer:
        stats = pstats.Stats(timer)
        value['profile'] = [dict(file=k[0], line=k[1], function=k[2], calls=v[1], cumulative=v[3])
                            for k, v in sorted(stats.stats.items(), key=lambda item: item[1][3], reverse=True)[:22]]
    connection.send(value)
    connection.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--genes', type=int, default=10000)
    parser.add_argument('--blocks', type=int, default=4)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--profile', action='store_true')
    parser.add_argument('--baseline-only', action='store_true')
    parser.add_argument('--pipeline', action='store_true', help='Include locus matching, selection and feature summaries.')
    args = parser.parse_args()
    assert args.genes > 0 and args.blocks > 0 and args.repeats > 0
    sources = {'baseline': args.baseline.read_text()}
    if not args.baseline_only:
        sources['candidate'] = Path(gff.__file__).read_text()
    frames = fixture(args.genes, args.blocks)
    result = {key: [] for key in sources}
    profiles = {}
    for repeat in range(args.repeats + 1):
        for key in (list(sources) if repeat % 2 else list(reversed(sources))):
            parent, child = multiprocessing.Pipe(duplex=False)
            process = multiprocessing.get_context('fork').Process(
                target=measure, args=(sources[key], frames, args.profile and repeat == 0, args.pipeline, child))
            process.start()
            child.close()
            value = parent.recv()
            process.join()
            parent.close()
            assert process.exitcode == 0
            if repeat:
                result[key].append(value)
            elif args.profile:
                profiles[key] = value
    signatures = {(v['checksum'], v['rows'], tuple(v['dtypes'])) for values in result.values() for v in values}
    assert len(signatures) == 1, 'Complete output frames differ'
    print(json.dumps(dict(genes=args.genes, blocks=args.blocks, repeats=args.repeats, warmups=1, pipeline=args.pipeline,
        environment=dict(python=platform.python_version(), platform=platform.platform(), pandas=pandas.__version__),
        source_sha256={key: hashlib.sha256(source.encode()).hexdigest() for key, source in sources.items()},
        equivalent=True, profiles=profiles, runs=result, medians={key: {metric: statistics.median(v[metric] for v in values)
        for metric in ('seconds', 'cpu_seconds', 'peak_rss_kib')} for key, values in result.items()})))


if __name__ == '__main__':
    main()
