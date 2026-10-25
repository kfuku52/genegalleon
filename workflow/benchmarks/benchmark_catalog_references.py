#!/usr/bin/env python3
"""Compare native indexed-genome CDS reconstruction and complete catalog calls.

Every timed helper uses pysam.FastaFile opened through indexed_genome. Catalog
profiles include source hashing, GFF parsing, scratch indexing and all publication
writes. Imports and fixture creation are outside timers. Fresh workers alternate
after a warmup pair; full sequence lists, catalog data, publication bytes, native
errors and input hashes are compared outside the timers. Peak RSS is Linux VmHWM.
"""
import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
import os
import pickle
import platform
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve()
while not (ROOT / 'workflow/support').is_dir() and ROOT != ROOT.parent:
    ROOT = ROOT.parent
sys.path.insert(0, str(ROOT / 'workflow/support'))
PROFILES = {
    'small': ('helper', 1, 128, 4, '+'),
    'chromosomal': ('helper', 24, 20000, 4, '+'),
    'fragmented': ('helper', 10000, 20000, 4, '+'),
    'reverse': ('helper', 10000, 10000, 4, '-'),
    'exon_rich': ('helper', 24, 1000, 64, '+'),
    'catalog_small': ('catalog', 1, 16, 1, '+'),
    'catalog_chromosomal': ('catalog', 24, 3000, 1, '+'),
    'catalog_fragmented': ('catalog', 10000, 3000, 1, '+'),
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    spec = importlib.util.spec_from_file_location('catalog_candidate', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_fixture(root, profile):
    kind, contigs, models, _blocks, _strand = PROFILES[profile]
    root.mkdir(parents=True, exist_ok=True)
    genome, gff, cds = (root/name for name in ('genome.fa', 'source.gff', 'source.cds.fa'))
    with genome.open('w') as handle:
        for i in range(contigs):
            handle.write(f'>chr{i:06}\n' + 'ATGAAACCCGGGCCCTAA' * 96 + '\n')
    if kind == 'catalog':
        with gff.open('w') as annotations, cds.open('w') as sequences:
            annotations.write('##gff-version 3\n')
            for i in range(models):
                # Uniformly include the full reference-name search range.
                seqid = f'chr{(i * 7919) % contigs:06}'
                gene, transcript = f'g{i:06}', f't{i:06}'
                for feature, identifier, parent, phase in [('gene', gene, '', '.'),
                                                          ('mRNA', transcript, gene, '.'),
                                                          ('CDS', f'c{i:06}', transcript, '0')]:
                    attributes = 'ID=' + identifier + (';Parent=' + parent if parent else '')
                    annotations.write(f'{seqid}\tsource\t{feature}\t1\t18\t.\t+\t{phase}\t{attributes}\n')
                sequences.write(f'>{gene}\nATGAAACCCGGGCCCTAA\n')
    return genome, gff, cds


def error_cases(module, genome):
    class StringName(str):
        pass

    result = []
    cases = [([], 'chr000000', '+'), ([], 'absent', '+'), ([[0, 18, 0]], 'absent', '+'),
             ([[0, 18, 0]], 'chr000000', '?'), ([[-1, 18, 0]], 'chr000000', '+'),
             ([[18, 18, 0]], 'chr000000', '+'), ([[0, 99999, 0]], 'chr000000', '+'),
             ([[0, 18]], 'chr000000', '+'), ([[0, 18, 0]], None, '+'),
             ([[0, 18, 0]], ['chr000000'], '+'),
             ([[0, 18, 0]], StringName('chr000000'), '+'),
             ([[0, 18, 0]], StringName('absent'), '+')]
    for blocks, seqid, strand in cases:
        try:
            result.append(('value', module.reconstruct_sequence(blocks, seqid, strand, genome)))
        except Exception as error:
            result.append(('error', type(error).__name__, str(error)))
    return result


def protocol_cases(module, native_genome):
    class DuckGenome:
        references = ['chr000000']

        def __contains__(self, _seqid):
            raise AssertionError('Duck genome __contains__ must not be called')

        def get_reference_length(self, _seqid):
            return 18

        def fetch(self, _seqid, _start, _end):
            return 'ATGAAACCCGGGCCCTAA'

    # No native-method assumption is imposed on duck-typed callers, even when
    # their __contains__ disagrees with the original references protocol.
    result = error_cases(module, DuckGenome())
    # A caller using only the references protocol need not import pysam.
    native_module = sys.modules.pop('pysam', None)
    try:
        result.extend(error_cases(module, DuckGenome()))
        result.extend(error_cases(module, native_genome))
        assert 'pysam' not in sys.modules, 'Reconstruction imported a new dependency'
    finally:
        if native_module is not None:
            sys.modules['pysam'] = native_module
    return result


def worker(args):
    module = load(args.source)
    kind, contigs, models, block_count, strand = PROFILES[args.profile]
    root = Path(args.fixture)
    paths = [root/name for name in ('genome.fa', 'source.gff', 'source.cds.fa')]
    inputs = {p.name: digest(p) for p in paths if p.exists()}
    output = io.StringIO()
    native_references_type = native_references_doc = None
    with contextlib.redirect_stdout(output):
        if kind == 'helper':
            blocks = [[i*24, i*24+18, i % 3] for i in range(block_count)]
            if strand == '-':
                blocks.reverse()
            identifiers = [f'chr{(i*7919) % contigs:06}' for i in range(models)]
            with module.indexed_genome(paths[0]) as genome:
                native_references_type = type(genome.references).__module__ + '.' + type(genome.references).__qualname__
                native_references_doc = type(genome).references.__doc__
                start_cpu, start = time.process_time(), time.perf_counter()
                value = [module.reconstruct_sequence(blocks, seqid, strand, genome) for seqid in identifiers]
                wall, cpu = time.perf_counter()-start, time.process_time()-start_cpu
                rss = next(int(line.split()[1]) for line in Path('/proc/self/status').read_text().splitlines()
                           if line.startswith('VmHWM:'))
                errors = error_cases(module, genome) + protocol_cases(module, genome)
            with module.indexed_genome(paths[0]) as closed_genome:
                pass
            errors.extend(error_cases(module, closed_genome))
            publications = None
        else:
            start_cpu, start = time.process_time(), time.perf_counter()
            value = module.build_catalog('Species_one', paths[2], paths[1], paths[0])
            module.write_catalog(value, args.publication)
            wall, cpu = time.perf_counter()-start, time.process_time()-start_cpu
            rss = next(int(line.split()[1]) for line in Path('/proc/self/status').read_text().splitlines()
                       if line.startswith('VmHWM:'))
            errors = None
            publications = {p.name: p.read_bytes() for p in Path(args.publication).iterdir() if p.is_file()}
            assert value['summary']['loci'] == models and value['summary']['usable_candidates'] == models
    assert inputs == {p.name: digest(p) for p in paths if p.exists()}
    assert not any((root/(name+'.fai')).exists() for name in ('genome.fa', 'source.cds.fa'))
    with open(args.payload, 'wb') as handle:
        pickle.dump((value, publications, errors, output.getvalue()), handle, protocol=5)
    print(json.dumps({'wall_seconds': wall, 'cpu_seconds': cpu, 'rss_kib': rss,
                      'source_sha256': digest(args.source),
                      'native_references_type': native_references_type,
                      'native_references_doc': native_references_doc}))


def main(args):
    import pysam
    report = {'platform': platform.platform(), 'python': platform.python_version(), 'pysam': pysam.__version__,
              'repeats': args.repeats, 'harness_sha256': digest(__file__),
              'sources': {label: {'path': getattr(args, label), 'sha256': digest(getattr(args, label))}
                          for label in ('baseline', 'candidate')}, 'profiles': {}}
    with tempfile.TemporaryDirectory(prefix='gg-catalog-references-') as directory:
        root = Path(directory)
        for profile in args.profiles:
            fixture = root/profile
            make_fixture(fixture, profile)
            runs = {'baseline': [], 'candidate': []}
            for trial in range(args.repeats+1):
                payloads = {}
                for label in (['baseline', 'candidate'] if trial % 2 == 0 else ['candidate', 'baseline']):
                    payload = root/(label+'.pkl')
                    publication = root/(label+'-output-'+profile+'-'+str(trial))
                    command = [sys.executable, str(Path(__file__).resolve()), '--worker', '--source',
                               getattr(args, label), '--baseline', args.baseline, '--profile', profile, '--fixture', str(fixture),
                               '--payload', str(payload), '--publication', str(publication)]
                    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', PYTHONHASHSEED='0')
                    env.pop('GG_GENOME_INDEX_CACHE', None)
                    process = subprocess.run(command, capture_output=True, text=True, check=True, env=env)
                    measurement = json.loads(process.stdout)
                    with payload.open('rb') as handle:
                        payloads[label] = pickle.load(handle)
                    if trial:
                        runs[label].append(measurement)
                assert payloads['baseline'] == payloads['candidate'], profile
            report['profiles'][profile] = runs
            print(profile, {label: statistics.median(row['wall_seconds'] for row in values)
                            for label, values in runs.items()}, flush=True)
    Path(args.output).write_text(json.dumps(report, indent=2)+'\n')
    print('All sequences, complete catalog/publication bytes, errors and source hashes match.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', required=True)
    parser.add_argument('--candidate', default=str(ROOT / 'workflow/support/gene_model_catalog.py'))
    parser.add_argument('--output', default='tmp/catalog-references-benchmark.json')
    parser.add_argument('--profiles', nargs='+', choices=PROFILES, default=list(PROFILES))
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--source')
    parser.add_argument('--profile')
    parser.add_argument('--fixture')
    parser.add_argument('--payload')
    parser.add_argument('--publication')
    args = parser.parse_args()
    if args.worker:
        worker(args)
    else:
        args.baseline = str(Path(args.baseline).resolve(strict=True))
        args.candidate = str(Path(args.candidate).resolve(strict=True))
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        main(args)
