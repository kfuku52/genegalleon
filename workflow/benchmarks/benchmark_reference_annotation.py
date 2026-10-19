#!/usr/bin/env python3
"""Compare WGD indexing/pair selection and complete Swiss-Prot annotation output."""

import argparse
import contextlib
import hashlib
import io
import json
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CASES = ('ks-placement', 'cross-child-pairs', 'uniprot')
SIZES = dict(species=128, anchors=4000, copies=4, families=8, references=200000, queries=1000)


def balanced_tree(names):
    if len(names) == 1:
        return names[0] + ':1'
    middle = len(names) // 2
    return '(' + balanced_tree(names[:middle]) + ',' + balanced_tree(names[middle:]) + '):1'


def worker(args):
    sys.path.insert(0, str(args.support))
    counts = {}
    with tempfile.TemporaryDirectory(prefix='gg-reference-annotation-') as temporary:
        root = Path(temporary)
        if args.case == 'ks-placement':
            import nwkit.clade_index
            import wgd_evidence as module
            from nwkit.util import read_tree

            path = root / 'species.nwk'
            path.write_text(balanced_tree([f'Species_{i:04d}' for i in range(args.species)]) + ';\n')
            tree = read_tree(str(path), 'auto', True, quiet=True)
            original_index = nwkit.clade_index.CladeIndex
            index = original_index(tree)
            bounds = {}
            leaves = list(tree.leaves())
            for leaf in leaves:
                node, depth = leaf, 1
                while node.up is not None:
                    node = node.up
                    bounds[leaf.name, index.clade_id_for_node(node)] = dict(
                        status='ok', interval_status='ok', monotone_from_younger_node='yes',
                        ci_lower=depth - .1, ci_upper=depth + .1)
                    depth += 1
            # Include early, ancestral, overlapping and unresolved placements.
            values = (0, .5, .95, 1.5, 2.5, 100)
            requests = [(leaves[i % len(leaves)].name, values[i % len(values)]) for i in range(args.anchors)]
            counts['clade_index_builds'] = 0

            def counted_index(*a, **kw):
                counts['clade_index_builds'] += 1
                return original_index(*a, **kw)

            nwkit.clade_index.CladeIndex = counted_index
            started = time.perf_counter()
            if hasattr(module, 'KsBranchPlacer'):
                place = module.KsBranchPlacer(tree, bounds).place
                result = [place(focal, ks) for focal, ks in requests]
            else:
                result = [module.branch_for_ks(focal, ks, tree, bounds) for focal, ks in requests]
            seconds = time.perf_counter() - started
            assert len(result) == args.anchors
            payload = json.dumps(result).encode()
        elif args.case == 'cross-child-pairs':
            import wgd_ssd as module
            from nwkit.util import read_tree

            names = [[f'Species_{species:04d}_{side}_{copy}' for copy in range(args.copies)
                      for species in range(args.species)] for side in ('a', 'b')]
            path = root / 'gene.nwk'
            path.write_text('(' + ','.join(balanced_tree(side) for side in names) + ');\n')
            tree = read_tree(str(path), 'auto', True, quiet=True)
            positions = {name: dict(species='_'.join(name.split('_')[:2])) for side in names for name in side}
            # Unmapped leaves remain excluded without changing pair order.
            positions.pop(names[0][-1])
            positions.pop(names[1][0])
            started = time.perf_counter()
            if hasattr(module, 'same_species_cross_child_pairs'):
                result = [module.same_species_cross_child_pairs(tree, positions) for _ in range(args.families)]
            else:
                result = [[(a.name, b.name) for a in tree.children[0].leaves() for b in tree.children[1].leaves()
                           if a.name in positions and b.name in positions
                           and positions[a.name]['species'] == positions[b.name]['species']]
                          for _ in range(args.families)]
            seconds = time.perf_counter() - started
            counts['pairs'] = sum(map(len, result))
            payload = json.dumps(result).encode()
        else:
            import reformat_uniprot_diamond as module

            reference, metadata = root / 'reference.fa', root / 'metadata.tsv'
            with reference.open('w') as fasta, metadata.open('w') as table:
                table.write('accession\t' + '\t'.join(module.METADATA_COLUMNS) + '\n')
                for i in range(args.references):
                    accession = f'P{i:07d}'
                    fasta.write(f'>{accession} Reference protein {i}\n' + ('M' + 'ACDEFGHIKL' * 29) + '\n')
                    table.write(accession + '\t' + '\t'.join(f'{column}_{i}' for column in module.METADATA_COLUMNS) + '\n')
            queries, hits, output = root / 'query.fa', root / 'hits.tsv', root / 'annotation.tsv'
            with queries.open('w') as fasta, hits.open('w') as table:
                for i in range(args.queries):
                    accession = f'P{i * args.references // args.queries:07d}'
                    fasta.write(f'>gene{i}\nMKT\n')
                    table.write(f'gene{i}\t{accession}\t90\t290\t1e-20\t200\t300\n')
            counts['reference_bytes'] = reference.stat().st_size + metadata.stat().st_size
            sys.argv = [str(args.support / 'reformat_uniprot_diamond.py'), '--query_fasta', str(queries),
                        '--diamond_tsv', str(hits), '--uniprot_fasta', str(reference),
                        '--uniprot_meta_tsv', str(metadata), '--outfile', str(output)]
            with contextlib.redirect_stdout(io.StringIO()):
                started = time.perf_counter()
                module.main()
                seconds = time.perf_counter() - started
            payload = output.read_bytes()
            counts['output_rows'] = len(payload.splitlines()) - 1
            assert counts['output_rows'] == args.queries
    return dict(seconds=seconds, peak_rss_kib_linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                output_sha256=hashlib.sha256(payload).hexdigest(), counts=counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-support', type=Path)
    parser.add_argument('--baseline-only', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cases', choices=CASES, nargs='+', default=list(CASES))
    parser.add_argument('--repeats', type=int, default=3)
    for name, default in SIZES.items():
        parser.add_argument('--' + name, type=int, default=default)
    parser.add_argument('--support', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--case', choices=CASES, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if min(args.repeats, *(getattr(args, name) for name in SIZES)) < 1 or args.species < 2:
        parser.error('sizes and repeats must be positive; species must be at least two')
    if args.queries > args.references:
        parser.error('queries must not exceed references')
    if args.support:
        print(json.dumps(worker(args)))
        return
    if args.baseline_only and not args.baseline_support:
        parser.error('--baseline-only requires --baseline-support')
    supports = {} if args.baseline_only else dict(current=ROOT / 'workflow/support')
    if args.baseline_support:
        supports = dict(baseline=args.baseline_support.resolve(), **supports)
    result = dict(environment=dict(platform=platform.platform(), python=sys.version),
                  configuration={name: getattr(args, name) for name in (*SIZES, 'repeats')}, warmups=1,
                  workloads={}, source_sha256={label: {name: hashlib.sha256((support / name).read_bytes()).hexdigest()
                      for name in ('wgd_evidence.py', 'wgd_ssd.py', 'reformat_uniprot_diamond.py')}
                      for label, support in supports.items()})
    for case in args.cases:
        rows = {label: [] for label in supports}
        fingerprints = set()
        for trial in range(args.repeats + 1):
            items = list(supports.items())
            if trial % 2:
                items.reverse()
            for label, support in items:
                command = [sys.executable, str(Path(__file__).resolve()), '--support', str(support),
                           '--case', case, '--output', str(args.output)]
                for name in SIZES:
                    command += ['--' + name, str(getattr(args, name))]
                completed = subprocess.run(command, text=True, capture_output=True)
                if completed.returncode:
                    raise RuntimeError(completed.stdout + completed.stderr)
                row = json.loads(completed.stdout)
                fingerprints.add(row['output_sha256'])
                if trial:
                    rows[label].append(row)
        if len(fingerprints) != 1:
            raise ValueError('Different outputs: ' + case)
        result['workloads'][case] = {label: dict(trials=values,
            median_seconds=statistics.median(row['seconds'] for row in values),
            median_peak_rss_kib_linux=statistics.median(row['peak_rss_kib_linux'] for row in values))
            for label, values in rows.items()}
        print(case, {label: result['workloads'][case][label]['median_seconds'] for label in supports}, flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
