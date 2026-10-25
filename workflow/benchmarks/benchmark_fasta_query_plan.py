#!/usr/bin/env python3
"""Compare the real indexed FASTA CLI with sparse/dense requested keys.

Fixtures, source/index building, imports, query-plan inspection and verification
are outside timers. Workers call the real parser/main/extract in fresh processes.
All comparisons use the same source/index fixture and verify full FASTA bytes,
exit codes and logs. Linux VmHWM includes imports and output buffers. Synthetic
warm-cache cases measure indexed extraction, not whole-workflow performance.
"""
import argparse
import ast
import contextlib
import hashlib
import importlib.util
import io
import json
import os
import platform
import random
import sqlite3
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / 'workflow/support'
PROFILES = ['small', 'small_all', 'sparse', 'requests512', 'requests1024', 'requests1025', 'normalized1024', 'requests4096',
            'dense', 'all', 'variants', 'missing', 'empty']


def load(path):
    sys.path.insert(0, str(SUPPORT))
    spec = importlib.util.spec_from_file_location('fasta_join_measured', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def query_plan(module, args):
    """Inspect the actual extract SQL, with actual requested keys, after timing."""
    source = ast.parse(Path(module.__file__).read_text())
    node = next(item for item in ast.walk(source) if isinstance(item, ast.JoinedStr)
                and 'SELECT s.identifier, s.header, s.sequence, s.species' in ast.unparse(item))
    column = 'identifier_lower' if '--ignore-case' in args.flags else 'identifier'
    extract_node = next(item for item in source.body if isinstance(item, ast.FunctionDef) and item.name == 'extract')
    clauses = [item for item in ast.walk(extract_node) if isinstance(item, ast.Assign)
               and any(isinstance(target, ast.Name) and target.id == 'from_clause' for target in item.targets)]
    assert len(clauses) <= 1
    requested = module.read_patterns(Path(args.patterns))
    with sqlite3.connect(f'file:{args.database}?mode=ro', uri=True) as connection:
        connection.execute('CREATE TEMP TABLE requested(pattern TEXT PRIMARY KEY)')
        if '--query-variants' in args.flags:
            species = [str(row[0]) for row in connection.execute('SELECT DISTINCT species FROM sequences')]
            requested += [f'{sp}_{pattern}' for sp in species for pattern in list(requested) if sp]
        if '--ignore-case' in args.flags:
            requested = [pattern.lower() for pattern in requested]
        connection.executemany('INSERT OR IGNORE INTO requested(pattern) VALUES (?)', ((pattern,) for pattern in requested))
        requested_count = connection.execute('SELECT COUNT(*) FROM requested').fetchone()[0]
        namespace = dict(vars(module), identifier_column=column, requested_count=requested_count,
                         requested=requested, connection=connection)
        if clauses:
            exec(compile(ast.Module(body=clauses, type_ignores=[]), '<actual-from-clause>', 'exec'), namespace)
        query = eval(compile(ast.Expression(node), '<actual-extract-query>', 'eval'), namespace)
        return list(connection.execute('EXPLAIN QUERY PLAN ' + query)), requested_count


def worker(args):
    source_sha256 = digest(args.implementation)
    module = load(args.implementation)
    assert digest(args.implementation) == source_sha256
    output = Path(args.worker_output)
    sys.argv = [str(args.implementation), 'extract', '--database', args.database,
                '--pattern-file', args.patterns, '--output', str(output), *args.flags]
    stdout, stderr = io.StringIO(), io.StringIO()
    wall, cpu = time.perf_counter(), time.process_time()
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
        code = module.main()
    wall, cpu = time.perf_counter() - wall, time.process_time() - cpu
    peak = next(int(line.split()[1]) for line in Path('/proc/self/status').read_text().splitlines()
                if line.startswith('VmHWM:'))
    value = dict(code=code, stdout=stdout.getvalue().replace(str(output), '$OUTPUT'),
                 stderr=stderr.getvalue().replace(str(output), '$OUTPUT'),
                 output_sha256=digest(output) if output.exists() else None,
                 output_bytes=output.stat().st_size if output.exists() else None)
    plan, requested_count = (query_plan(module, args)
                             if code == 0 and Path(args.patterns).stat().st_size else ([], None))
    print(json.dumps(dict(wall_seconds=wall, cpu_seconds=cpu, peak_rss_kib=peak,
                          source_sha256=source_sha256, result=value,
                          requested_entries=requested_count, plan=plan)))


def fixture(module, directory, count):
    directory.mkdir(parents=True, exist_ok=True)
    sources = []
    rng = random.Random(12017)
    # Reuse a small pool of random DNA bodies, keeping fixture creation bounded.
    # SQLite still retains one compressed record per identifier as in production.
    bodies = [''.join(rng.choices('ACGT', k=900)) for _ in range(64)]
    per_species = max(1, count // 8)
    for species_no in range(8):
        species = f'Genus_species{species_no}'
        source = directory / (species + '.fa')
        with source.open('w') as handle:
            for index in range(per_species):
                identifier = f'{species}_gene{index:07d}'
                handle.write(f'>{identifier} synthetic record {index}\n')
                sequence = bodies[(index + species_no) % len(bodies)]
                for offset in range(0, len(sequence), 80):
                    handle.write(sequence[offset:offset + 80] + '\n')
        sources.append((source, species))
    database, manifest = directory / 'sequences.sqlite3', directory / 'manifest.json'
    with contextlib.redirect_stdout(io.StringIO()):
        module.build_store(database, manifest, sources, module.DEFAULT_MAX_DATABASE_BYTES,
                           module.DEFAULT_MINIMUM_FREE_BYTES)
    assert module.database_content_current(database, refresh_state=False)
    return database, manifest, sources, per_species


def controller(args):
    module = load(args.baseline)
    root = Path(args.output).with_suffix('')
    root.mkdir(parents=True, exist_ok=True)
    fixtures = {}
    records = []
    sources = {'baseline': args.baseline.resolve(), 'candidate': args.candidate.resolve()}
    source_identities = {name: {'path': str(path), 'sha256': digest(path)} for name, path in sources.items()}
    env = dict(os.environ, PYTHONHASHSEED='0', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    for profile in args.profiles:
        count = 64 if profile in {'small', 'small_all'} else args.records
        if count not in fixtures:
            fixtures[count] = fixture(module, root / str(count), count)
        database, manifest, source_files, per_species = fixtures[count]
        before = {str(path): digest(path) for path in [database, manifest, *[path for path, _ in source_files]]}
        patterns = root / (profile + '.txt')
        flags = []
        if profile in {'all', 'small_all'}:
            patterns.write_text(''.join(f'Genus_species{species}_gene{index:07d}\n'
                                       for species in range(8) for index in range(per_species)))
            flags = ['--require-all']
        elif profile == 'normalized1024':
            identifiers = [f'Genus_species{index % 8}_gene{index:07d}' for index in range(1024)]
            patterns.write_text(''.join(identifier + '\n' + identifier.upper() + '\n' for identifier in identifiers))
            flags = ['--ignore-case']
        elif profile in {'small', 'sparse', 'dense', 'requests512', 'requests1024', 'requests1025', 'requests4096'}:
            requested_number = (int(profile.removeprefix('requests')) if profile.startswith('requests')
                                else (4 if profile == 'small' else (32 if profile == 'sparse' else per_species)))
            number = min(per_species, requested_number)
            patterns.write_text(''.join(f'Genus_species{index % 8}_gene{index:07d}\n' for index in range(number)))
            flags = ['--require-all']
        elif profile == 'variants':
            patterns.write_text('gene0000000\ngene0000003\ngene0000007\n')
            flags = ['--ignore-case', '--query-variants', '--prefix-species']
        elif profile == 'missing':
            patterns.write_text('absent\nGenus_species0_gene0000001\n')
            flags = ['--require-all']
        elif profile == 'empty':
            patterns.write_text('')
            flags = ['--query-variants']
        else:
            raise ValueError(profile)
        pairs = []
        for trial in range(args.repeats + 1):
            pair = {}
            for name in (list(sources) if trial % 2 == 0 else list(reversed(sources))):
                output = root / f'{profile}-{name}.fa'
                output.unlink(missing_ok=True)
                cmd = [sys.executable, str(Path(__file__).resolve()), '--baseline', str(sources['baseline']),
                       '--worker', '--implementation', str(sources[name]),
                       '--database', str(database), '--patterns', str(patterns), '--worker-output', str(output),
                       '--flags', *flags]
                raw = subprocess.run(cmd, check=True, capture_output=True, text=True, env=env)
                pair[name] = json.loads(raw.stdout)
                assert pair[name]['source_sha256'] == source_identities[name]['sha256'], pair[name]
            assert pair['baseline']['result'] == pair['candidate']['result'], (profile, pair)
            if profile != 'missing':
                assert pair['baseline']['result']['code'] == 0, (profile, pair)
            else:
                assert pair['baseline']['result']['code'] == 3, pair
            assert (root / f'{profile}-baseline.fa').read_bytes() == (root / f'{profile}-candidate.fa').read_bytes()
            if trial:
                pairs.append(pair)
        after = {path: digest(path) for path in before}
        assert before == after
        assert module.database_content_current(database, refresh_state=False)
        records.append(dict(profile=profile, pairs=pairs, input_hashes=before,
                            pattern_sha256=digest(patterns), flags=flags, indexed_records=per_species * 8))
        print(profile, {name: statistics.median(row[name]['wall_seconds'] for row in pairs) for name in sources}, flush=True)
    summary = {}
    for record in records:
        pairs = record['pairs']
        summary[record['profile']] = {
            name: {key: statistics.median(row[name][key] for row in pairs)
                   for key in ['wall_seconds', 'cpu_seconds', 'peak_rss_kib']} for name in sources}
        summary[record['profile']]['paired_speedups'] = [row['baseline']['wall_seconds'] / row['candidate']['wall_seconds'] for row in pairs]
    report = dict(python=sys.version, platform=platform.platform(), repeats=args.repeats,
                  warmup_pairs=1, sources=source_identities,
                  benchmark={'path': str(Path(__file__).resolve()), 'sha256': digest(__file__)},
                  records=args.records, summaries=summary, measurements=records)
    Path(args.output).write_text(json.dumps(report, indent=2) + '\n')
    print('Complete FASTA bytes, return codes, stdout/stderr and unmodified inputs match.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline', type=Path, required=True, help='Saved baseline fasta_sequence_store.py')
    parser.add_argument('--candidate', type=Path, default=SUPPORT / 'fasta_sequence_store.py')
    parser.add_argument('--records', type=int, default=100000)
    parser.add_argument('--profiles', nargs='+', choices=PROFILES, default=PROFILES)
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--output', type=Path, default=Path('tmp/fasta-query-plan-benchmark.json'))
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--implementation', type=Path)
    parser.add_argument('--database')
    parser.add_argument('--patterns')
    parser.add_argument('--worker-output')
    parser.add_argument('--flags', nargs=argparse.REMAINDER, default=[])
    args = parser.parse_args()
    if args.records < 8 or args.repeats < 1:
        parser.error('--records must be at least 8 and --repeats must be positive')
    if args.worker:
        worker(args)
    else:
        controller(args)
