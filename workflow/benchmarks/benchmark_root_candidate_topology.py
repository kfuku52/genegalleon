#!/usr/bin/env python3
"""Compare candidate topology validation and the complete root adapter.

Run through workflow/tests/run_in_runtime.sh without concurrent tests or builds::

    python workflow/benchmarks/benchmark_root_candidate_topology.py \
        --baseline /mounted/baseline/species_tree_guided_gene_tree_rooting.py \
        --output tmp/root-candidate-comparison.json

One warmup and alternating fresh processes measure wall/CPU time and Linux peak
RSS. Topology workloads use 128 candidate validations at 128/512 balanced tips
and 512 pectinate tips; the full adapter uses 128 actual rerooted candidates.
The full call includes lazy NWKIT imports, native rooting, rendering, and bundle
publication. Input generation and verification are outside all timed sections.
Topology output is restored to canonical taxon partitions before hashing. Full
output equality covers stdout, selected Newick, TSV, and PDF; only temporary
paths in stdout and PDF CreationDate/ModDate metadata are normalized.

The saved baseline's nested function is compiled from its original AST with
its original closure variables. Current implementations expose the same check
as a module helper. These workloads isolate candidate-validation costs and do
not establish whole-workflow performance or SIF compatibility.
"""
import argparse
import ast
import contextlib
import hashlib
import importlib.metadata
import importlib.util
import io
import json
import platform
import re
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'workflow/support'))


def newick(tips, shape):
    nodes = [f'T{i:04d}:{1 + (i % 11) / 10}' for i in range(tips)]
    if shape == 'pectinate':
        result = nodes[0]
        for node in nodes[1:]:
            result = f'({result},{node}):0.5'
    else:
        while len(nodes) > 1:
            nodes = [f'({nodes[i]},{nodes[i + 1]}):0.5' if i + 1 < len(nodes) else nodes[i]
                     for i in range(0, len(nodes), 2)]
        result = nodes[0]
    return result.removesuffix(':0.5') + ';'


def topology_worker(args):
    from nwkit.clade_mapping import canonical_split
    from nwkit.util import read_tree

    tree = read_tree(newick(args.tips, args.shape), 'auto', True, quiet=True)
    leaves = frozenset(tree.leaf_names())
    labels = sorted(leaves)
    leaf_bits = {name: 1 << index for index, name in enumerate(labels)}
    source = ast.parse(args.source.read_text())
    function = next(node for node in ast.walk(source) if isinstance(node, ast.FunctionDef) and node.name in ('topology_splits', '_topology_splits'))
    namespace = dict(canonical_split=canonical_split, leaves=leaves, leaf_bits=leaf_bits,
                     leaf_mask=(1 << len(labels)) - 1)
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(args.source), 'exec'), namespace)
    cpu_started = time.process_time()
    started = time.perf_counter()
    for _ in range(args.candidates):
        result = (namespace['_topology_splits'](tree, leaf_bits, (1 << len(labels)) - 1)
                  if function.name == '_topology_splits' else namespace['topology_splits'](tree))
    elapsed = time.perf_counter() - started
    cpu_elapsed = time.process_time() - cpu_started
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if result and isinstance(next(iter(result)), int):
        result = {canonical_split(frozenset(label for index, label in enumerate(labels) if mask & (1 << index)),
                                  frozenset(label for index, label in enumerate(labels) if not mask & (1 << index)))
                  for mask in result}
    splits = sorted((sorted(left), sorted(right)) for left, right in result)
    payload = json.dumps(splits, sort_keys=True).encode()
    return dict(seconds=elapsed, cpu_seconds=cpu_elapsed, sha256=hashlib.sha256(payload).hexdigest(), splits=len(result),
                peak_rss_kib_linux=peak_rss)


def full_worker(args):
    from ete4 import Tree

    with tempfile.TemporaryDirectory(prefix='gg-root-benchmark-') as temporary:
        work = Path(temporary)
        source = work / 'source.nwk'
        source.write_text(newick(args.tips, args.shape))
        tree = Tree(source.read_text(), parser=1)
        candidates = work / 'candidates.nwk'
        texts = []
        for index in range(1, min(args.candidates + 1, len(list(tree.traverse())))):
            copy = tree.copy(method='deepcopy')
            copy.set_outgroup(list(copy.traverse())[index])
            texts.append(copy.write(parser=1, format_root_node=True))
        candidates.write_text('\n'.join(texts) + '\n')
        options = SimpleNamespace(in_tree=source, candidate_trees=candidates, out_tree=work/'selected.nwk',
            comparison_table=work/'comparison.tsv', comparison_plot=work/'comparison.pdf', species_parser='taxonomic')
        spec = importlib.util.spec_from_file_location('root_adapter_under_test', args.source)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        capture = io.StringIO()
        with contextlib.redirect_stdout(capture):
            cpu_started = time.process_time()
            children_started = resource.getrusage(resource.RUSAGE_CHILDREN)
            started = time.perf_counter()
            module.select_root(options)
            elapsed = time.perf_counter() - started
        stdout = re.sub(re.escape(str(work)) + r'/nwkit-root-[^/\s]+', '{work}/nwkit-root-{temporary}', capture.getvalue())
        stdout = stdout.replace(str(work), '{work}')
        payload = stdout.encode() + options.out_tree.read_bytes() + options.comparison_table.read_bytes() + re.sub(
            rb'/(CreationDate|ModDate) \([^)]*\)', b'', options.comparison_plot.read_bytes())
        children = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu_elapsed = time.process_time() - cpu_started + children.ru_utime + children.ru_stime - children_started.ru_utime - children_started.ru_stime
        return dict(seconds=elapsed, cpu_seconds=cpu_elapsed, sha256=hashlib.sha256(payload).hexdigest(), candidates=len(texts),
                    peak_rss_kib_linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                    child_peak_rss_kib_linux=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, default=ROOT / 'workflow/support/species_tree_guided_gene_tree_rooting.py')
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--worker', choices=['topology', 'full'])
    parser.add_argument('--shape', choices=['balanced', 'pectinate'], default='balanced')
    parser.add_argument('--tips', type=int, default=128)
    parser.add_argument('--candidates', type=int, default=128)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--full-repeats', type=int, default=1)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--cases', nargs='+', choices=['topology', 'full'], default=['topology', 'full'])
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(topology_worker(args) if args.worker == 'topology' else full_worker(args)))
        return
    results = dict(environment=dict(python=sys.version, platform=platform.platform(),
                   nwkit_version=importlib.metadata.version('nwkit')), runs=[], sources={})
    workloads = []
    if 'topology' in args.cases:
        workloads += [('topology', shape, tips, candidates) for shape, tips, candidates in
                      [('balanced', 128, 128), ('balanced', 512, 128), ('pectinate', 512, 128)]]
    if 'full' in args.cases:
        workloads += [('full', 'balanced', 128, 128)]
    implementations = [('before', args.baseline), ('after', args.source)] if args.baseline else [('before', args.source)]
    for label, source in implementations:
        results['sources'][label] = dict(path=str(source), sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    for case, shape, tips, candidates in workloads:
        repeats = args.full_repeats if case == 'full' else args.repeats
        for repeat in range(-1, repeats):
            for label, source in implementations if repeat % 2 else implementations[::-1]:
                command = [sys.executable, str(Path(__file__)), '--source', str(source), '--worker', case,
                           '--shape', shape, '--tips', str(tips), '--candidates', str(candidates)]
                completed = subprocess.run(command, capture_output=True, text=True)
                if completed.returncode:
                    print(completed.stdout, completed.stderr, file=sys.stderr)
                    completed.check_returncode()
                row = json.loads(completed.stdout.splitlines()[-1])
                row.update(implementation=label, case=case, shape=shape, tips=tips,
                           candidate_count=candidates, repeat=repeat)
                results['runs'].append(row)
                print(json.dumps(row), flush=True)
                if args.output:
                    args.output.write_text(json.dumps(results, indent=2) + '\n')
        measured = [row for row in results['runs'] if row['case'] == case and row['shape'] == shape and row['tips'] == tips]
        assert len({row['sha256'] for row in measured}) == 1, 'Output equivalence failed'
    results['summary'] = []
    for case, shape, tips, candidates in workloads:
        summary = dict(case=case, shape=shape, tips=tips, candidates=candidates)
        for label, _ in implementations:
            rows = [row for row in results['runs'] if row['repeat'] >= 0 and row['implementation'] == label
                    and row['case'] == case and row['shape'] == shape and row['tips'] == tips]
            summary[label] = {metric: statistics.median(row[metric] for row in rows)
                              for metric in ('seconds', 'cpu_seconds', 'peak_rss_kib_linux')}
        results['summary'].append(summary)
    if args.output:
        args.output.write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()
