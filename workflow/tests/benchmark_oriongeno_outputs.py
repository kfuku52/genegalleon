#!/usr/bin/env python3
"""Compare native GFF export/reading in fresh workers with an optional old support snapshot."""
import argparse
import hashlib
import json
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def fixture(path, genes):
    with path.open('w') as stream:
        stream.write('##gff-version 3\n')
        for index in range(genes):
            start, strand = index * 1000, '+' if index % 2 == 0 else '-'
            gene, tid = 'g' + str(index), 't' + str(index)
            def row(kind, left, right, attrs, phase='.', start=start, strand=strand):
                stream.write(f'chr\tsynthetic\t{kind}\t{start+left+1}\t{start+right}\t.\t{strand}\t{phase}\t{attrs}\n')
            row('gene', 0, 940, 'ID=' + gene)
            row('mRNA', 0, 940, 'ID=' + tid + ';Parent=' + gene)
            for i, (left, right) in enumerate([(0, 30), (80, 200), (300, 420), (500, 620), (700, 820), (900, 940)]):
                row('exon', left, right, f'ID={tid}.e{i};Parent={tid}')
            for i, (left, right) in enumerate([(100, 200), (300, 420), (500, 620), (700, 800)]):
                row('CDS', left, right, f'ID={tid}.c{i};Parent={tid}', '0')


def worker(args):
    root = Path(__file__).resolve().parents[1] / 'support'
    sys.path[:0] = [str(args.support), str(root)]
    import annevo_runtime as native
    import oriongeno_outputs as outputs
    import transcript_utrs as utrs
    started = time.perf_counter()
    if args.case == 'validated-export':
        native.gff_counts(args.input)
    value = (outputs.write_explicit_utrs(args.input, args.output, lengths={'chr': args.genes * 1000})
             if args.case != 'read' else utrs.load_models(args.input))
    if args.case == 'validated-export':
        value = {'export': value, 'counts': native.gff_counts(args.output)}
    seconds = time.perf_counter() - started
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    digest = hashlib.sha256()
    if args.case != 'read':
        with args.output.open('rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(block)
        digest.update(json.dumps(value, sort_keys=True).encode())
    else:
        for key in sorted(value):
            digest.update(json.dumps([key, value[key]], sort_keys=True).encode())
    sources = {Path(m.__file__).name: hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest() for m in (outputs, utrs, native)}
    print(json.dumps(dict(seconds=seconds, peak_rss_kib=rss, digest=digest.hexdigest(), sources=sources)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-support', type=Path)
    parser.add_argument('--genes', type=int, default=20000)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--support', type=Path)
    parser.add_argument('--input', type=Path)
    parser.add_argument('--output', type=Path)
    cases = ['export', 'read', 'validated-export']
    parser.add_argument('--case', choices=cases)
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    if platform.system() != 'Linux' or args.genes < 1 or args.repeats < 1:
        parser.error('Use Linux in the GeneGalleon runtime, positive genes and positive repeats')
    current = Path(__file__).resolve().parents[1] / 'support'
    supports = {'current': current}
    if args.baseline_support:
        supports = {'baseline': args.baseline_support, **supports}
    records = []
    with tempfile.TemporaryDirectory(prefix='gg-utr-benchmark-') as temporary:
        temporary = Path(temporary)
        source, output = temporary / 'input.gff3', temporary / 'output.gff3'
        fixture(source, args.genes)
        for case in cases:
            for repetition in range(args.repeats + 1):
                order = list(supports) if repetition % 2 == 0 else list(reversed(supports))
                for version in order:
                    argv = [sys.executable, '-B', str(Path(__file__).resolve()), '--worker', '--case', case,
                            '--support', str(supports[version]), '--input', str(source), '--output', str(output),
                            '--genes', str(args.genes)]
                    result = subprocess.run(argv, capture_output=True, text=True, check=True, timeout=300)
                    records.append(dict(case=case, version=version, warmup=repetition == 0,
                                        repetition=repetition, **json.loads(result.stdout)))
        summary = {}
        for case in cases:
            assert len({r['digest'] for r in records if r['case'] == case}) == 1, 'Benchmark outputs differ'
            summary[case] = {}
            for version in supports:
                rows = [r for r in records if r['case'] == case and r['version'] == version and not r['warmup']]
                summary[case][version] = dict(seconds=statistics.median(r['seconds'] for r in rows),
                    peak_rss_kib=statistics.median(r['peak_rss_kib'] for r in rows), digest=rows[0]['digest'])
        print(json.dumps(dict(schema_version=1, platform=platform.platform(), python=sys.version,
            genes=args.genes, input_bytes=source.stat().st_size, repeats=args.repeats, warmups=1,
            records=records, summary=summary), sort_keys=True))


if __name__ == '__main__':
    main()
