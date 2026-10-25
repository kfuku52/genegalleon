#!/usr/bin/env python3
"""Compare complete GFF line normalization, including repair/error audit output.

Imports and fixture generation are outside operation timers. The file controls
include reading, normalise_line and output writing with raw or gzip transport.
Every fresh worker executes a differential contract corpus outside its timer.
Synthetic warm-cache fixtures isolate normalization and file I/O, rather than
whole-workflow throughput. Linux VmHWM includes imports and timed I/O buffers;
output hashing and the contract corpus run after the memory observation.
"""
import argparse
import gzip
import hashlib
import importlib.util
import io
import json
import os
import platform
import random
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROFILES = ['small', 'structural', 'gzip', 'metadata', 'annotated', 'repairs', 'late_error']


def source_identity(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def load(path):
    spec = importlib.util.spec_from_file_location('gff_syntax_measured', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def contract(module):
    cases = [
        ('', '.', 'gene'), ('.', '.', 'gene'), ('ID=one', 'x', 'gene'),
        ('ID=one;Parent=two;Name=three;', 'maker', 'mRNA'),
        ('ID=one;;Parent=two;', 'x', 'gene'), (' ID=one; Parent=two ', 'x', 'gene'),
        ('gene_id "one"; transcript_id "two";', 'gtf', 'CDS'),
        ('gene_id "one" transcript_id "two"', 'gtf', 'CDS'),
        ('ID=one;Name:ABC=DEF', 'AUGUSTUS', 'gene'),
        ('ID=one;Blast2Go:ABC=DEF', 'AUGUSTUS', 'mRNA'),
        ('ID=one;Name=ABC;1', 'funannotate', 'gene'),
        ('ID=one;product=ABC;2, variant 1', 'funannotate', 'mRNA'),
        ('ID=one;Note=A RecName: Full=protein; Short=alpha', 'x', 'mRNA'),
        ('ID=one;Note=This (a; b) protein', 'maker', 'mRNA'),
        ('ID=one;Note=This (a; b protein', 'maker', 'mRNA'),
        ('ID=one;Note=Cation transporter HKT1;1 (OsHKT1;1)', 'maker', 'mRNA'),
        ('ID=one;Note=nucleosome assembly protein 1;1-like [Species example]', 'x', 'mRNA'),
        ('ID=one;description=endo-1,3;1,4-beta protein', 'x', 'gene'),
        ('ID=one;Accession=GWHGABCD00001.1;HC', 'EVM', 'gene'),
        ('ID=one;note=Geneious type: Repeat Region; IRa', '.', 'misc_feature'),
        ('ID=one;orphan', 'x', 'gene'), ('ID=one; =invalid', 'x', 'gene'),
        ('ID=one;Parent=', 'x', 'CDS'), ('ID=', 'x', 'gene'),
        ('Name=a=b\nsecond=line;Parent=c', 'x', 'CDS'),
    ]
    rng = random.Random(61471)
    values = ['a', '', 'αβ', '%3B%2C', 'a=b', 'x\ty', 'one two', '\n', 'Name:ABC=DEF', ' Note=foo']
    for _ in range(300):
        keys = rng.choices(['ID', 'Parent', 'Name', 'gene_id', 'transcript_id'], k=rng.randrange(1, 6))
        text = ';'.join(key+'='+rng.choice(values) for key in keys) + rng.choice(['', ';'])
        cases.append((text, rng.choice(['maker', 'AUGUSTUS', 'funannotate', '.', 'x']), rng.choice(['gene', 'mRNA', 'CDS'])))
    result = []
    for text, source, feature in cases:
        for bare in (False, True):
            try:
                value = module.normalise_attributes(text, source, feature, allow_bare=bare)
                result.append(['ok', value])
            except Exception as exc:
                result.append([type(exc).__name__, str(exc)])
    lines = ['# comment\n', '\n', 'chr1\tx\tCDS\t1\t9\t.\t+\t0\tID=x;Parent=y\t\t\n',
             'chr1\tx\tCDS\t1\t9\t.\t+\t0\tID=x;orphan\n', 'bad\trow\n']
    for number, line in enumerate(lines, 1):
        changes = []
        try:
            value = module.normalise_line(line, 'contract.gff3', number, changes)
            result.append(['line', value, changes])
        except Exception as exc:
            result.append(['line-error', type(exc).__name__, str(exc), changes])
    return hashlib.sha256(json.dumps(result, ensure_ascii=False, sort_keys=True).encode()).hexdigest(), len(result)


def worker(args):
    implementation = source_identity(args.implementation)
    module = load(args.implementation)
    assert source_identity(args.implementation) == implementation
    is_gzip = str(args.fixture).endswith('.gz')
    opener = gzip.open if is_gzip else open
    error = None
    changes = []
    wall, cpu = time.perf_counter(), time.process_time()
    try:
        with opener(args.fixture, 'rt', encoding='utf-8') as source:
            with open(args.worker_output, 'wb') as raw:
                if is_gzip:
                    binary = gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0)
                else:
                    binary = raw
                with io.TextIOWrapper(binary, encoding='utf-8') as target:
                    for number, line in enumerate(source, 1):
                        target.write(module.normalise_line(line, 'input.gff3', number, changes))
    except Exception as exc:
        error = [type(exc).__name__, str(exc)]
    cpu, wall = time.process_time()-cpu, time.perf_counter()-wall
    peak = next(int(line.split()[1]) for line in Path('/proc/self/status').read_text().splitlines() if line.startswith('VmHWM:'))
    result = Path(args.worker_output).read_bytes()
    signature = hashlib.sha256(result + json.dumps(changes, ensure_ascii=False, sort_keys=True).encode()
                               + json.dumps(error).encode()).hexdigest()
    corpus_hash, corpus_cases = contract(module)
    print(json.dumps(dict(wall_seconds=wall, cpu_seconds=cpu, peak_rss_kib=peak,
                          output_sha256=signature, contract_sha256=corpus_hash, contract_cases=corpus_cases,
                          implementation=implementation,
                          error=error, changes=len(changes), output_bytes=len(result))))


def fixture(path, profile):
    sizes = {'small': 64, 'structural': 20000, 'gzip': 20000, 'metadata': 10000,
             'annotated': 10000, 'repairs': 3000, 'late_error': 1000}
    genes = sizes[profile]
    opener = gzip.open if path.name.endswith('.gz') else open
    with opener(path, 'wt', encoding='utf-8') as out:
        out.write('##gff-version 3\n')
        for gene in range(genes):
            position = gene*10000+1
            annotation = (f';Dbxref=GeneID:{gene};Note=An ordinary protein description {gene}'
                          if profile == 'annotated' else '')
            if profile == 'repairs':
                source, attr = 'funannotate', f'ID=G{gene};Name=Protein;1'
            elif profile == 'metadata':
                source, attr = 'fixture', f'ID=G{gene};Name=G{gene};Note=An ordinary protein description {gene};Dbxref=GeneID:{gene}'
            else:
                source, attr = 'fixture', f'ID=G{gene};Name=G{gene}'
            out.write(f'chr1\t{source}\tgene\t{position}\t{position+999}\t.\t+\t.\t{attr}{annotation}\n')
            out.write(f'chr1\tfixture\tmRNA\t{position}\t{position+999}\t.\t+\t.\tID=T{gene};Parent=G{gene}{annotation}\n')
            for exon in range(6):
                start = position+exon*150
                out.write(f'chr1\tfixture\tCDS\t{start}\t{start+89}\t.\t+\t0\tParent=T{gene}{annotation}\n')
        if profile == 'late_error':
            out.write('chr1\tfixture\tCDS\t1\t9\t.\t+\t0\tID=x;unrecoverable\n')


def main(args):
    if args.worker:
        worker(args)
        return
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    sources = {label: source_identity(getattr(args, label)) for label in ('baseline', 'candidate')}
    report = dict(platform=platform.platform(), python=platform.python_version(), repeats=args.repeats,
                  warmup_pairs=1, sources=sources, benchmark=source_identity(__file__), profiles={})
    with tempfile.TemporaryDirectory(prefix='gg-gff-attributes-') as temporary:
        root = Path(temporary)
        for profile in args.profiles:
            source = root / (profile+'.gff3'+('.gz' if profile == 'gzip' else ''))
            fixture(source, profile)
            samples = []
            for trial in range(-1, args.repeats):
                pair = {}
                for label in (['baseline', 'candidate'] if trial % 2 == 0 else ['candidate', 'baseline']):
                    cmd = [sys.executable, __file__, '--baseline', sources['baseline']['path'],
                           '--worker', '--implementation', sources[label]['path'],
                           '--fixture', str(source), '--worker-output', str(root / (label+'.out'))]
                    run = subprocess.run(cmd, check=True, text=True, capture_output=True,
                                         env=dict(os.environ, PYTHONHASHSEED='0'))
                    pair[label] = json.loads(run.stdout)
                    assert pair[label]['implementation'] == sources[label], pair[label]
                    if profile == 'late_error':
                        assert pair[label]['error'] and pair[label]['error'][0] == 'ValueError', pair[label]
                    else:
                        assert pair[label]['error'] is None, pair[label]
                    assert pair[label]['output_bytes'] > 0, pair[label]
                for key in ['output_sha256', 'contract_sha256', 'contract_cases', 'error', 'changes', 'output_bytes']:
                    assert pair['baseline'][key] == pair['candidate'][key], (profile, key, pair)
                if trial >= 0:
                    samples.append(pair)
            report['profiles'][profile] = dict(samples=samples,
                input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(), input_bytes=source.stat().st_size,
                before=statistics.median(pair['baseline']['wall_seconds'] for pair in samples),
                after=statistics.median(pair['candidate']['wall_seconds'] for pair in samples))
            print(profile, report['profiles'][profile]['before'], report['profiles'][profile]['after'], flush=True)
            output.write_text(json.dumps(report, indent=2)+'\n')
    print('Complete outputs, audits, errors and contract corpus match.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', required=True, help='Saved baseline gff_attribute_syntax.py')
    parser.add_argument('--candidate', default=str(ROOT / 'workflow/support/gff_attribute_syntax.py'))
    parser.add_argument('--profiles', nargs='+', choices=PROFILES, default=PROFILES)
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--output', default='tmp/gff-structural-attributes-benchmark.json')
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--implementation')
    parser.add_argument('--fixture')
    parser.add_argument('--worker-output')
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    main(args)
