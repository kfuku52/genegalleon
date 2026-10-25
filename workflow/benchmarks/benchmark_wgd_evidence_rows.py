#!/usr/bin/env python3
"""Compare complete WGD family-evidence scans with exact typed outputs/errors.

Imports, fixtures and verification are outside operation timers. Complete
positions/events/anchors are parsed by the actual loader, including rows outside
the selected family. Every worker also compares the full contract corpus.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import pickle
import platform
import statistics
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / 'workflow/support'
FILES = ['gene_positions.tsv', 'wgd_events.tsv', 'anchor_evidence.tsv']
PROFILES = {'small': (256, 512, 32, 0), 'sparse': (200000, 400000, 32, 0),
            'dense': (20000, 40000, 20000, 0), 'wide': (100000, 200000, 32, 16)}


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load(path):
    sys.path.insert(0, str(SUPPORT))
    spec = importlib.util.spec_from_file_location('wgd_loader_measured', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def typed(value):
    if isinstance(value, dict):
        return [type(value).__name__, [[typed(key), typed(item)] for key, item in value.items()]]
    if isinstance(value, (tuple, list)):
        return [type(value).__name__, [typed(item) for item in value]]
    return [type(value).__name__, value]


def fingerprint(value):
    return hashlib.sha256(json.dumps(typed(value), ensure_ascii=False, separators=(',', ':')).encode()).hexdigest()


def write_table(path, fields, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(fields)
        writer.writerows(rows)


def fixture(directory, shape):
    directory.mkdir(parents=True, exist_ok=True)
    positions, anchors, selected, extras = shape
    extra_names = [f'annotation_{index}' for index in range(extras)]
    notes = ['quoted\tvalue', 'continued\nvalue', 'α β description', '']
    write_table(directory / FILES[0], ['gene_id', 'species', 'locus_id', 'seqid', 'rank', 'start', 'end', 'note', *extra_names],
                ([f'gene_{index:07d}', f'Species_{index % 32:03d}', f'locus_{index}',
                  f'chr{index % 8}', str(index), str(index * 100), str(index * 100 + 90),
                  notes[index % 4], *[f'{column}_{index}' for column in extra_names]] for index in range(positions)))
    write_table(directory / FILES[1], ['species_event_id', 'event_support', 'note'],
                [[f'event{index}', 'WGD-supported' if index % 3 else 'unresolved', notes[index % 4]] for index in range(32)])
    write_table(directory / FILES[2], ['species', 'block_id', 'gene_a', 'gene_b', 'ks',
                                     'ks_status', 'species_event_id', 'placement_status', 'note', *extra_names],
                ([f'Species_{index % 32:03d}', f'block{index}', f'gene_{index % positions:07d}',
                  f'gene_{(index + 32) % positions:07d}', '0.3', 'ok', f'event{index % 32}',
                  'interval_supported', notes[index % 4],
                  *[f'{column}_{index}' for column in extra_names]] for index in range(anchors)))
    # Sparse families retain two copies in each selected species. Dense control
    # retains every gene and every anchor, exposing dictionary-construction cost.
    selected_indices = (range(positions) if selected == positions
                        else [index // 2 + (index % 2) * 32 for index in range(selected)])
    (directory / 'genes.json').write_text(json.dumps([f'gene_{index:07d}' for index in selected_indices]))


def contracts(directory):
    directory.mkdir(parents=True, exist_ok=True)
    base = {
        FILES[0]: b'gene_id\tspecies\nkeep_a\tA\nkeep_b\tA\noutside\tB\n',
        FILES[1]: b'species_event_id\tevent_support\nevent_A\tunresolved\n',
        FILES[2]: b'gene_a\tgene_b\tnote\nkeep_a\tkeep_b\tfirst\nkeep_b\tkeep_a\tsecond\noutside\tkeep_a\tignored\n',
    }
    cases = [('ordinary', {}), ('blank_rows', {
        FILES[0]: b'gene_id\tspecies\n\nkeep_a\tA\n\nkeep_b\tA\n',
        FILES[2]: b'gene_a\tgene_b\n\nkeep_a\tkeep_b\n\n'}),
        ('quoted_empty', {FILES[0]: b'gene_id\n""\nkeep_a\n'}),
        ('missing_positions_key_empty', {FILES[0]: b'other\tspecies\n'}),
        ('missing_positions_key', {FILES[0]: b'other\tspecies\nx\tA\n'}),
        ('missing_anchor_a_empty', {FILES[2]: b'other\tgene_b\n'}),
        ('missing_anchor_a', {FILES[2]: b'other\tgene_b\nx\tkeep_b\n'}),
        ('missing_anchor_b', {FILES[2]: b'gene_a\tother\nkeep_a\tx\n'}),
        ('missing_both_anchor_keys', {FILES[2]: b'other\tx\na\tb\n'}),
        ('width_before_missing_key', {FILES[2]: b'other\tgene_b\na\tb\tc\n'}),
        ('global_duplicate', {FILES[0]: base[FILES[0]] + b'outside\tC\n'}),
        ('global_duplicate_event', {FILES[1]: base[FILES[1]] + b'outside\tunresolved\noutside\tunresolved\n'}),
        ('duplicate_before_bad_anchor', {FILES[0]: base[FILES[0]] + b'outside\tC\n', FILES[2]: b''}),
        ('missing_event_key', {FILES[1]: b'other\tevent_support\nx\tunresolved\n'}),
        ('missing_event_key_empty', {FILES[1]: b'other\tevent_support\n'}),
        ('quoted_unicode', {FILES[0]: 'gene_id\tspecies\tnote\nkeep_a\tA\t"α\tβ\nγ"\nkeep_b\tA\t""\n'.encode()}),
    ]
    for table in (FILES[0], FILES[2]):
        cases.extend((f'{table}_{name}', {table: value}) for name, value in [
            ('empty', b''), ('blank_header', b'\n'), ('duplicate_header', b'gene_id\tgene_id\n'),
            ('empty_header_name', b'gene_a\t\n'), ('short_tail', base[table] + b'short\n'),
            ('long_tail', base[table] + b'a\tb\tc\td\n'),
            ('bad_quote_tail', base[table] + b'"unterminated'),
            ('invalid_utf8_tail', base[table] + b'bad\t\xff\n'),
        ])
    for name, changes in cases:
        target = directory / name
        target.mkdir(exist_ok=True)
        for file, data in base.items():
            (target / file).write_bytes(changes.get(file, data))
    return [name for name, _ in cases]


def call(module, directory, genes):
    try:
        return ['ok', module._load_family_evidence(directory, genes)]
    except Exception as error:
        return ['error', type(error).__name__, str(error)]


def worker(args):
    module = load(args.implementation)
    genes = set(json.loads((args.fixture / 'genes.json').read_text()))
    wall, cpu = time.perf_counter(), time.process_time()
    result = module._load_family_evidence(args.fixture, genes)
    wall, cpu = time.perf_counter() - wall, time.process_time() - cpu
    peak = next(int(line.split()[1]) for line in Path('/proc/self/status').read_text().splitlines()
                if line.startswith('VmHWM:'))
    assert isinstance(result[0], dict) and isinstance(result[1], dict)
    assert isinstance(result[2], defaultdict) and result[2].default_factory is list
    contract_results = [(directory.name, call(module, directory, {'keep_a', 'keep_b'}))
                        for directory in sorted(args.contracts.iterdir())]
    payload = [result, contract_results]
    with args.payload.open('wb') as stream:
        pickle.dump(payload, stream, protocol=5)
    print(json.dumps(dict(wall_seconds=wall, cpu_seconds=cpu, peak_rss_kib=peak,
                          output_sha256=fingerprint(result), contracts_sha256=fingerprint(contract_results),
                          source_sha256=sha(args.implementation), selected_positions=len(result[0]),
                          selected_anchors=sum(map(len, result[2].values())))))


def controller(args):
    work = args.output.with_suffix('')
    work.mkdir(parents=True, exist_ok=True)
    names = contracts(work / 'contracts')
    sources = {'baseline': args.baseline.resolve(), 'candidate': args.candidate.resolve()}
    env = dict(os.environ, PYTHONHASHSEED='0', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    report = dict(platform=platform.platform(), python=sys.version, repeats=args.repeats, warmups=1,
                  contract_cases=names, sources={name: {'path': str(path), 'sha256': sha(path)} for name, path in sources.items()},
                  profiles={})
    for profile in args.profiles:
        directory = work / profile
        fixture(directory, PROFILES[profile])
        before = {name: sha(directory / name) for name in [*FILES, 'genes.json']}
        samples = []
        for trial in range(args.repeats + 1):
            pair, outputs = {}, {}
            for name in (list(sources) if trial % 2 == 0 else list(reversed(sources))):
                payload = work / (name + '.pickle')
                command = [sys.executable, str(Path(__file__).resolve()), '--baseline', str(args.baseline),
                           '--worker', '--implementation', str(sources[name]),
                           '--fixture', str(directory), '--contracts', str(work / 'contracts'), '--payload', str(payload)]
                raw = subprocess.run(command, check=True, capture_output=True, text=True, env=env)
                pair[name] = json.loads(raw.stdout)
                with payload.open('rb') as stream:
                    outputs[name] = pickle.load(stream)
            assert pair['baseline']['output_sha256'] == pair['candidate']['output_sha256'], (profile, pair)
            assert pair['baseline']['contracts_sha256'] == pair['candidate']['contracts_sha256'], (profile, pair)
            assert typed(outputs['baseline']) == typed(outputs['candidate']), profile
            assert pair['baseline']['selected_positions'] == PROFILES[profile][2], pair
            if trial:
                samples.append(pair)
        assert before == {name: sha(directory / name) for name in before}
        summary = {name: {key: statistics.median(row[name][key] for row in samples)
                          for key in ['wall_seconds', 'cpu_seconds', 'peak_rss_kib']} for name in sources}
        summary['paired_speedups'] = [row['baseline']['wall_seconds'] / row['candidate']['wall_seconds'] for row in samples]
        report['profiles'][profile] = dict(shape=PROFILES[profile], input_sha256=before, samples=samples, summary=summary)
        print(profile, summary, flush=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print('Complete typed outputs, field/row/key order, errors and unchanged inputs match.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, default=SUPPORT / 'wgd_ssd.py')
    parser.add_argument('--profiles', nargs='+', choices=PROFILES, default=list(PROFILES))
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--output', type=Path, default=Path('tmp/wgd-evidence-rows-benchmark.json'))
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--implementation', type=Path)
    parser.add_argument('--fixture', type=Path)
    parser.add_argument('--contracts', type=Path)
    parser.add_argument('--payload', type=Path)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    worker(args) if args.worker else controller(args)
