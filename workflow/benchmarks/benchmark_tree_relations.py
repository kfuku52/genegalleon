#!/usr/bin/env python3
"""Compare relation helpers and complete real orthogroup-statistics CLI processes.

Helper profiles time only relation annotation including branch-ID assignment.
cli_* profiles run the actual source with its real runtime dependencies and time
process startup, imports, main, and both output writes. Fixtures, full byte/error/
log comparisons, and input/source invariance checks remain outside the timers.
Use --profiles cli_balanced_512 cli_ladder_512 cli_ladder_1536 for CLI confirmation.
"""
import argparse
import ast
import hashlib
import json
import os
import pickle
import platform
import re
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = next(parent for parent in Path(__file__).resolve().parents
            if (parent/'workflow/support/orthogroup_statistics.py').is_file())


def functions(path):
    import pandas
    parsed = ast.parse(Path(path).read_text(), filename=str(path))
    main = next(n for n in parsed.body if isinstance(n, ast.FunctionDef) and n.name == 'main')
    names = {'iter_ancestors', 'iter_descendants', 'iter_leaves', 'node_is_leaf'}
    nested = {'get_node_prop', 'set_node_prop', 'get_node_label', 'ensure_branch_ids', 'get_ancestor_and_descendant'}
    selected = [n for n in parsed.body if isinstance(n, ast.FunctionDef) and n.name in names]
    selected += [n for n in main.body if isinstance(n, ast.FunctionDef) and n.name in nested]
    namespace = {'pandas': pandas}
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(path), 'exec'), namespace)
    return namespace


def fixture(profile):
    import ete4
    import pandas as pd
    kind, size = profile.rsplit('_', 1)
    size = int(size)
    tree = ete4.PhyloTree()
    if kind == 'single':
        tree.name = 'NA'
    elif kind in {'star', 'duplicate'}:
        for i in range(size):
            tree.add_child(name='NA' if kind == 'duplicate' else f'g{i:06}')
    elif kind in {'ladder', 'unary'}:
        current = tree
        for i in range(size):
            if kind == 'ladder':
                current.add_child(name=f'g{i:06}')
            current = current.add_child()
        current.name = f'g{size:06}'
    else:
        pending = [(tree, 0, size)]
        while pending:
            node, start, stop = pending.pop()
            if stop-start == 1:
                node.name = f'g{start:06}'
            else:
                mid = (start+stop)//2
                left, right = node.add_child(), node.add_child()
                pending.extend([(right, mid, stop), (left, start, mid)])
    frame = pd.DataFrame({'branch_id': [*range(sum(1 for _ in tree.traverse())), -1, 0],
                          'existing': 'value'})
    frame.index = pd.Index(range(len(frame)-1, -1, -1), name='retained-index')
    return tree, frame


def source_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def newick_text(tree):
    # Iterative serialization also covers ladders beyond Python's recursion limit.
    tokens = []
    pending = [';', tree]
    while pending:
        item = pending.pop()
        if isinstance(item, str):
            tokens.append(item)
        elif item.is_leaf:
            tokens.append(item.name + ':1')
        else:
            tokens.append('(')
            pending.append(')' + str(item.name or '') + ':1')
            children = item.get_children()
            for index in range(len(children)-1, -1, -1):
                pending.append(children[index])
                if index:
                    pending.append(',')
    return ''.join(tokens) + '\n'


def normalize_cli_log(text, source):
    text = text.replace(str(Path(source).resolve()), '<source>')
    text = re.sub(r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}(?:\.\d+)?', '<timestamp>', text)
    return text


def cli_worker(args):
    before = {'source': source_hash(args.source), 'fixture': source_hash(args.fixture)}
    output = Path(args.cli_directory)
    output.mkdir()
    env = os.environ.copy()
    env['PYTHONPATH'] = os.pathsep.join(filter(None, [args.support_directory, env.get('PYTHONPATH', '')]))
    command = [sys.executable, str(Path(args.source).resolve()), '--rooted_tree', args.fixture]
    start = time.perf_counter()
    completed = subprocess.run(command, cwd=output, capture_output=True, text=True, env=env)
    wall = time.perf_counter()-start
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    payload = {
        'returncode': completed.returncode,
        'stdout': normalize_cli_log(completed.stdout, args.source),
        'stderr': normalize_cli_log(completed.stderr, args.source),
        'outputs': {str(p.relative_to(output)): p.read_bytes() for p in sorted(output.rglob('*.tsv'))},
    }
    assert before == {'source': source_hash(args.source), 'fixture': source_hash(args.fixture)}, 'source/input mutation'
    with open(args.payload, 'wb') as handle:
        pickle.dump(payload, handle, protocol=5)
    print(json.dumps({'wall_seconds': wall, 'cpu_seconds': usage.ru_utime+usage.ru_stime,
                      'rss_kib': usage.ru_maxrss, 'source_sha256': before['source'],
                      'fixture_sha256': before['fixture'], 'timer': 'complete_cli_process_including_imports'}))


def worker(args):
    if args.profile.startswith('cli_'):
        cli_worker(args)
        return
    before_source = source_hash(args.source)
    ns = functions(args.source)
    tree, frame = fixture(args.profile)
    before = pickle.dumps(frame, protocol=5)
    start_cpu = time.process_time()
    start = time.perf_counter()
    result = ns['get_ancestor_and_descendant'](tree, frame)
    wall = time.perf_counter()-start
    cpu = time.process_time()-start_cpu
    rss = next(int(line.split()[1]) for line in Path('/proc/self/status').read_text().splitlines()
               if line.startswith('VmHWM:'))
    assert pickle.dumps(frame, protocol=5) == before, 'input frame mutated'
    assert source_hash(args.source) == before_source, 'source mutated'
    state = [(n.name, ns['get_node_label'](n)) for n in tree.traverse()]
    with open(args.payload, 'wb') as handle:
        pickle.dump((result, state), handle, protocol=5)
    print(json.dumps({'wall_seconds': wall, 'cpu_seconds': cpu, 'rss_kib': rss,
                      'source_sha256': before_source, 'timer': 'relation_function_including_branch_ids'}))


def main(args):
    import pandas as pd
    report = {'platform': platform.platform(), 'python': platform.python_version(), 'pandas': pd.__version__,
              'repeats': args.repeats,
              'harness_sha256': source_hash(__file__),
              'sources': {name: {'path': getattr(args, name),
                                 'sha256': hashlib.sha256(Path(getattr(args, name)).read_bytes()).hexdigest()}
                          for name in ('baseline', 'candidate')}, 'profiles': {}}
    with tempfile.TemporaryDirectory(prefix='gg-tree-relations-') as tmp:
        for profile in args.profiles:
            cli = profile.startswith('cli_')
            if cli:
                tree, _frame = fixture(profile.removeprefix('cli_'))
                # Normal CLI inputs name internal nodes; missing names fail the
                # existing completeness guard before publication.
                for index, node in enumerate(tree.traverse()):
                    if not node.is_leaf and not node.name:
                        node.name = f'internal{index:06}'
                tree_path = Path(tmp)/(profile+'.nwk')
                tree_path.write_text(newick_text(tree))
                branch_rows = sum(1 for _ in tree.traverse())
            runs = {'baseline': [], 'candidate': []}
            for trial in range(args.repeats+1):
                payloads = {}
                for label in (['baseline', 'candidate'] if trial % 2 == 0 else ['candidate', 'baseline']):
                    payload = Path(tmp)/(label+'.pkl')
                    command = [sys.executable, str(Path(__file__).resolve()), '--worker', '--source',
                               getattr(args, label), '--profile', profile, '--payload', str(payload)]
                    if cli:
                        command += ['--fixture', str(tree_path), '--cli-directory', str(Path(tmp)/f'{profile}-{trial}-{label}'),
                                    '--support-directory', args.support_directory]
                    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', PYTHONHASHSEED='0')
                    process = subprocess.run(command, capture_output=True, text=True, check=True, env=env)
                    measured = json.loads(process.stdout)
                    assert measured['source_sha256'] == report['sources'][label]['sha256'], 'measured source changed'
                    with payload.open('rb') as handle:
                        payloads[label] = pickle.load(handle)
                    if trial:
                        runs[label].append(measured)
                if cli:
                    before, after = payloads['baseline'], payloads['candidate']
                    assert before['returncode'] == 0, (profile, before)
                    assert after['returncode'] == before['returncode'], (profile, after)
                    assert before['stdout'] == after['stdout'], (profile, before['stdout'], after['stdout'])
                    assert before['stderr'] == after['stderr'], (profile, before['stderr'], after['stderr'])
                    assert set(before['outputs']) == {'orthogroup.branch.tsv', 'orthogroup.tree.tsv'}, profile
                    assert before['outputs'] == after['outputs'], profile
                    assert len(before['outputs']['orthogroup.branch.tsv'].splitlines()) == branch_rows+1, profile
                    assert len(before['outputs']['orthogroup.tree.tsv'].splitlines()) == 2, profile
                    continue
                before, before_state = payloads['baseline']
                after, after_state = payloads['candidate']
                pd.testing.assert_frame_equal(before, after, check_exact=True)
                assert before.to_csv(sep='\t') == after.to_csv(sep='\t'), profile
                assert before_state == after_state, profile
            report['profiles'][profile] = runs
            print(profile, {label: statistics.median(row['wall_seconds'] for row in values)
                            for label, values in runs.items()}, flush=True)
    Path(args.output).write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline')
    parser.add_argument('--candidate', default=str(ROOT/'workflow/support/orthogroup_statistics.py'))
    parser.add_argument('--output', default='tree-relations.json')
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--profiles', nargs='+', default=['single_1', 'balanced_32', 'balanced_512', 'balanced_2048',
                                                       'ladder_512', 'ladder_1536', 'star_2048', 'unary_512', 'duplicate_32',
                                                       'cli_balanced_512', 'cli_ladder_512', 'cli_ladder_1536'])
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--source')
    parser.add_argument('--profile')
    parser.add_argument('--payload')
    parser.add_argument('--fixture')
    parser.add_argument('--cli-directory')
    parser.add_argument('--support-directory', default=str(ROOT/'workflow/support'))
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    if args.worker:
        worker(args)
    else:
        args.baseline = str(Path(args.baseline).resolve(strict=True))
        args.candidate = str(Path(args.candidate).resolve(strict=True))
        args.support_directory = str(Path(args.support_directory).resolve(strict=True))
        main(args)
