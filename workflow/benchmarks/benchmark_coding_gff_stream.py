"""Compare a supplied historical gff2genestat source with the bounded stream.

Run inside the GeneGalleon runtime; --baseline points to the old helper source.
The fixture contains linked alignment ancestors and compressed repetitive GFF.
"""
import argparse
import ast
import gzip
import hashlib
import json
import multiprocessing
import resource
import statistics
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
import gff2genestat as gff
import pandas


def measure(source, path, connection):
    namespace = dict(vars(gff))
    spool = []
    original_reader = gff.read_gff_table
    def reader(value):
        if isinstance(value, (str, Path)):
            spool.append(Path(value).stat().st_size)
        return original_reader(value)
    namespace['read_gff_table'] = reader
    function = next(n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef) and n.name == 'read_coding_gff_table')
    exec(compile(ast.Module(body=[function], type_ignores=[]), '<benchmark-helper>', 'exec'), namespace)
    started = time.perf_counter()
    table = namespace['read_coding_gff_table'](str(path))
    elapsed = time.perf_counter() - started
    checksum = hashlib.sha256(pandas.util.hash_pandas_object(table, index=True).values.tobytes()).hexdigest()
    connection.send(dict(seconds=elapsed, peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                         spool_bytes=sum(spool), rows=len(table), checksum=checksum, dtypes=list(map(str, table.dtypes))))
    connection.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--rows', type=int, default=250000)
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    sources = {'baseline': args.baseline.read_text(), 'stream': Path(gff.__file__).read_text()}
    results = {key: [] for key in sources}
    with tempfile.TemporaryDirectory(prefix='gg-gff-benchmark-') as temporary:
        path = Path(temporary) / 'fixture.gff.gz'
        with gzip.open(path, 'wt') as out:
            out.write('chr1\ts\tmatch\t1\t99999999\t.\t+\t.\tID=ancestor\n')
            for index in range(args.rows):
                out.write(f'chr1\ts\tCDS\t{index * 3 + 1}\t{index * 3 + 3}\t.\t+\t0\tID=c{index};Parent=ancestor\n')
        for repeat in range(args.repeats + 1):
            for key in (list(sources) if repeat % 2 else list(reversed(sources))):
                parent, child = multiprocessing.Pipe(duplex=False)
                process = multiprocessing.get_context('fork').Process(target=measure, args=(sources[key], path, child))
                process.start()
                child.close()
                value = parent.recv()
                process.join()
                parent.close()
                assert process.exitcode == 0
                if repeat:
                    results[key].append(value)
        signatures = {(v['checksum'], v['rows'], tuple(v['dtypes'])) for values in results.values() for v in values}
        assert len(signatures) == 1, 'Outputs differ'
        print(json.dumps(dict(rows=args.rows, repeats=args.repeats, warmups=1, equivalent=True, runs=results,
            medians={key: {metric: statistics.median(v[metric] for v in values)
                          for metric in ('seconds', 'peak_rss_kib', 'spool_bytes')} for key, values in results.items()})))


if __name__ == '__main__':
    main()
