#!/usr/bin/env python3
"""Compare full ortholog summaries and archive compaction in fresh workers.

Fixtures are generated with a saved baseline before timing. Each worker copies
its fixture, then times the real summary or compact command including its I/O,
verification and publication. Imports, copying and output comparison are outside
timers. Paired workers alternate; one warmup pair is discarded.
"""

import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
import weakref
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / 'workflow/support'
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SUPPORT))


def load_store(source):
    spec = importlib.util.spec_from_file_location('gene_family_output_store', source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def source_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def peak_rss():
    return next(int(line.split()[1]) * 1024 for line in Path('/proc/self/status').read_text().splitlines()
                if line.startswith('VmHWM:'))


def prepare_fixture(root, case, store):
    families_root = root / 'families'
    if case.startswith('reader'):
        from workflow.tests.test_query_gene_orthologs import weak_duplication_fixture

        seed = root / 'seed'
        seed.mkdir(parents=True)
        fixture = weak_duplication_fixture.__wrapped__(seed)
        count = {'reader-small': 8, 'reader-raw': 256, 'reader-normal': 256, 'reader-large': 1024}[case]
        query = root / 'queries'
        query.mkdir()
        families = [f'FAM{i:05d}' for i in range(count)]
        for subdir, suffix in [('stat_branch', '_stat.branch.tsv'), ('cds_fasta', '_cds.fasta')]:
            (families_root / subdir).mkdir(parents=True, exist_ok=True)
            for family in families:
                shutil.copyfile(fixture['output_root'] / subdir / ('FAM' + suffix),
                                families_root / subdir / (family + suffix))
        for family in families:
            shutil.copyfile(fixture['query_dir'] / 'FAM', query / family)
        if case != 'reader-raw':
            store.archive_completed_outputs(families_root, 'query2family', families,
                                            store.query_id_extractor(families), include_incomplete=True)
        return {'families': count, 'source_subdirs': 2, 'storage': 'raw' if case == 'reader-raw' else 'zip'}
    shards, members = {'compact-small': (8, 16), 'compact-normal': (64, 32),
                       'compact-large': (64, 64)}[case]
    raw = families_root / 'stat_branch'
    raw.mkdir(parents=True)
    artifacts = {}
    for shard in range(shards):
        paths = []
        for member in range(members):
            family = f'OG{member * shards + shard:07d}'
            path = raw / (family + '_stat.branch.tsv')
            path.write_text('gene\tannotation\n' + family + '\t' + 'annotation α ' * 16 + '\n')
            os.utime(path, ns=(1780000000000000000, 1780000000000000000))
            paths.append(path)
        _, records, _ = store._archive_chunk(families_root, families_root / 'archives', 'stat_branch',
                                             paths, 'orthogroup', shard + 1, store.orthogroup_id_from_name)
        artifacts.update((artifact.logical_path, artifact) for artifact in records)
    shutil.rmtree(raw)
    store._write_store_metadata(families_root, 'orthogroup', [a.family_id for a in artifacts.values()])
    logical = store.GeneFamilyOutputStore(families_root)
    logical._write_index(artifacts)
    logical._write_generation_counter(shards)
    return {'source_shards': shards, 'members_per_shard': members, 'total_members': shards * members,
            'source_order': 'interleaved family IDs across source shards'}


def fingerprint_archive(root, store):
    logical = store.GeneFamilyOutputStore(root)
    logical.verify(deep=True)
    result = []
    for subdir in logical.logical_subdirs():
        for artifact in logical.artifacts(subdir):
            with logical.open_binary(subdir, artifact.name) as handle:
                digest = hashlib.sha256(handle.read()).hexdigest()
            result.append({'logical_path': artifact.logical_path, 'generation': artifact.generation,
                           'size': artifact.size, 'crc': artifact.crc, 'sha256': artifact.sha256,
                           'mtime_ns': artifact.mtime_ns, 'family_id': artifact.family_id,
                           'payload_sha256': digest})
    logical._reset_cache()
    manifests = []
    zip_metadata = []
    for path in store._physical_archive_paths(root):
        manifest = store.GeneFamilyOutputStore(root)._read_manifest(path)
        manifests.append(manifest)
        with store.zipfile.ZipFile(path) as archive:
            infos = [{key: getattr(info, key) for key in
                      ('filename', 'compress_type', 'file_size', 'compress_size', 'CRC', 'date_time',
                       'external_attr', 'internal_attr', 'flag_bits', 'create_system', 'create_version',
                       'extract_version', 'volume')}
                     for info in archive.infolist()]
        zip_metadata.append({'relative_path': str(path.relative_to(root)),
                             'sha256': source_hash(path), 'members': infos})
    indexes = {}
    for directory in (store.INDEX_DIR_NAME, store.SUBDIR_INDEX_DIR_NAME):
        for path in sorted((logical.archive_root / directory).glob('*.json')):
            indexes[str(path.relative_to(logical.archive_root))] = json.loads(path.read_bytes())
    catalog = logical.archive_root / store.INDEX_CATALOG_FILE
    indexes[store.INDEX_CATALOG_FILE] = json.loads(catalog.read_bytes())
    epoch = (logical.archive_root / store.INDEX_EPOCH_FILE).read_text().strip()
    assert re.fullmatch(r'\d+-\d+-[0-9a-f]{32}', epoch)
    return {'artifacts': result, 'manifests': manifests, 'zip_metadata': zip_metadata, 'indexes': indexes}


@contextlib.contextmanager
def fixed_wall_clock():
    original_time, original_gmtime = time.time, time.gmtime
    time.time = lambda: 1780000000.0
    time.gmtime = lambda seconds=None: original_gmtime(1780000000 if seconds is None else seconds)
    try:
        yield
    finally:
        time.time, time.gmtime = original_time, original_gmtime


def worker(args):
    store = load_store(args.source)
    workspace = args.worker / 'fixture'
    shutil.copytree(args.fixture, workspace)
    counts = {'source_zip_opens': 0, 'max_live_source_readers': 0, 'max_fds': 0}
    references = []
    original = store.zipfile.ZipFile

    def observed(path, *positional, **keywords):
        archive = original(path, *positional, **keywords)
        mode = keywords.get('mode', positional[0] if positional else 'r')
        if mode == 'r' and '/archives/' in str(path) and '.part-' in str(path):
            counts['source_zip_opens'] += 1
            references.append(weakref.ref(archive))
            if counts['source_zip_opens'] <= 16 or counts['source_zip_opens'] % 16 == 0:
                counts['max_live_source_readers'] = max(counts['max_live_source_readers'],
                                                       sum(reference() is not None for reference in references))
                counts['max_fds'] = max(counts['max_fds'], len(list(Path('/proc/self/fd').iterdir())))
        return archive

    if args.instrument:
        store.zipfile.ZipFile = observed
    log = io.StringIO()
    if args.case.startswith('reader'):
        import query_gene_orthologs as summary

        output = workspace / 'output'
        options = SimpleNamespace(basis='reference_species', dir_gene_family=str(workspace / 'families'),
            dir_query_gene=str(workspace / 'queries'), family_file='', reference_species='Reference_species',
            dup_conf_score_threshold=.05,
            **{'out_' + name: str(output / (name + '.tsv')) for name in
               ('columns', 'glyphs', 'tree', 'synteny', 'ufboot', 'dup_conf', 'query_map', 'selection')})
        started, cpu_started = time.perf_counter(), time.process_time()
        with contextlib.redirect_stdout(log):
            summary.run(options)
        wall, cpu = time.perf_counter() - started, time.process_time() - cpu_started
        rss = peak_rss()
        result = {path.name: path.read_bytes().hex() for path in sorted(output.iterdir())}
    else:
        epoch_path = workspace / 'families' / store.STORE_DIR_NAME / store.INDEX_EPOCH_FILE
        before_epoch = epoch_path.read_bytes()
        with fixed_wall_clock():
            started, cpu_started = time.perf_counter(), time.process_time()
            with contextlib.redirect_stdout(log):
                store.compact_archives(workspace / 'families', 'orthogroup')
            wall, cpu = time.perf_counter() - started, time.process_time() - cpu_started
        rss = peak_rss()
        assert epoch_path.read_bytes() != before_epoch
        result = fingerprint_archive(workspace / 'families', store)
    store.zipfile.ZipFile = original
    (args.worker / 'output.json').write_text(json.dumps(result, ensure_ascii=False, sort_keys=True) + '\n')
    (args.worker / 'measurement.json').write_text(json.dumps({'wall_s': wall, 'cpu_s': cpu,
        'peak_rss_bytes': rss, 'counts': counts,
        'output_sha256': source_hash(args.worker / 'output.json')}) + '\n')


def main(args):
    sources = {'baseline': args.before.resolve()}
    if not args.baseline_only:
        sources['candidate'] = args.after.resolve()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result = {'environment': {'python': sys.version, 'platform': platform.platform()},
              'sources': {name: source_hash(path) for name, path in sources.items()},
              'benchmark_sha256': source_hash(Path(__file__)), 'pairs': args.pairs, 'cases': {}}
    result['instrumented'] = args.instrument
    result['compact_wall_clock_epoch'] = 1780000000
    with tempfile.TemporaryDirectory(prefix='gg-zip-readers-') as temporary:
        temporary = Path(temporary)
        store = load_store(sources['baseline'])
        for case in args.cases:
            fixture = temporary / ('fixture-' + case)
            configuration = prepare_fixture(fixture, case, store)
            samples = []
            reference = None
            for trial in range(-1, args.pairs):
                for name in (list(sources) if trial % 2 == 0 else list(reversed(sources))):
                    root = temporary / f'{case}-{trial}-{name}'
                    root.mkdir()
                    command = [sys.executable, str(Path(__file__).resolve()), '--worker', str(root),
                               '--source', str(sources[name]), '--fixture', str(fixture), '--case', case]
                    if args.instrument:
                        command.append('--instrument')
                    logfile = args.output.with_name(f'{args.output.stem}-{case}-{trial}-{name}.log')
                    with logfile.open('w') as log:
                        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
                    evidence = (root / 'output.json').read_bytes()
                    if reference is None:
                        reference = evidence
                    assert evidence == reference, (case, trial, name)
                    sample = json.loads((root / 'measurement.json').read_text())
                    sample.update(variant=name, trial=trial, warmup=trial < 0)
                    samples.append(sample)
                    print(json.dumps({'case': case, **sample}), flush=True)
            summary = {}
            for name in sources:
                rows = [row for row in samples if row['variant'] == name and not row['warmup']]
                summary[name] = {key: statistics.median(row[key] for row in rows)
                                 for key in ('wall_s', 'cpu_s', 'peak_rss_bytes')}
            result['cases'][case] = {'configuration': configuration, 'samples': samples, 'summary': summary}
            args.output.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before', type=Path)
    parser.add_argument('--after', type=Path, default=SUPPORT / 'gene_family_output_store.py')
    parser.add_argument('--baseline-only', action='store_true')
    parser.add_argument('--instrument', action='store_true', help='Separate structural reader/FD lifetime proof')
    parser.add_argument('--pairs', type=int, default=5)
    parser.add_argument('--cases', nargs='+', default=['reader-small', 'reader-raw', 'reader-normal', 'reader-large',
                                                     'compact-small', 'compact-normal', 'compact-large'])
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', type=Path)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--fixture', type=Path)
    parser.add_argument('--case')
    arguments = parser.parse_args()
    if arguments.worker:
        worker(arguments)
    else:
        if arguments.before is None or arguments.output is None or arguments.pairs < 1:
            parser.error('--before, --output and positive --pairs are required')
        main(arguments)
