#!/usr/bin/env python3
"""Compare strict declaration, wide trait-table and species-tree validation."""
import argparse
import csv
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


def worker(args):
    sys.path.insert(0, str(args.support_root))
    import artifact_provenance as provenance
    import prepare_orthogroup_copy_number as copy_number
    import reconciled_speciation_contrast as traits

    root = args.worker
    (root / 'inputs').mkdir()
    source = root / 'inputs/source.tsv'
    source.write_text('gene_id\tmetric\ngeneA\t7\n')
    declarations = [f'input_{i:05d}={source}' for i in range(args.size)]
    contract_args = provenance.build_parser().parse_args([
        'record', '--manifest', str(root / 'manifest.json'), '--step', 'benchmark', '--family-id', 'test',
        '--logical-root', str(root), '--workspace-root', str(root), '--output', f'result={source}',
    ])
    contract_args.input = declarations
    columns = [f'trait_{i:05d}' for i in range(args.size)]
    table = root / 'traits.tsv'
    with table.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(columns)
        for row in range(4):
            writer.writerow([str((row + col) % 7) for col in range(args.size)])
    tree_path = root / 'species.nwk'
    leaves = [f'species_{i:05d}' for i in range(args.size)]
    tree_path.write_text('(' + ','.join(f'{name}:1' for name in leaves) + ');\n')
    duplicate_tree = root / 'duplicate-species.nwk'
    duplicate_tree.write_text('(' + ','.join(f'{name}:1' for name in [*leaves, leaves[-1], leaves[0]]) + ');\n')
    timings = {}
    started = time.perf_counter()
    contract = provenance.build_contract(contract_args, include_diagnostics=False)
    # Cross-kind duplicate declarations must remain fatal before any I/O.
    contract_args.input_gene_family_store = [f'input_00000={root / "missing-store"}']
    try:
        provenance.build_contract(contract_args, include_diagnostics=False)
    except provenance.ProvenanceError as exc:
        contract_error = str(exc)
    else:
        raise AssertionError('Duplicate input label accepted')
    timings['provenance_contract'] = time.perf_counter() - started
    started = time.perf_counter()
    frame = traits.read_tsv(table)
    selected = traits._parse_csv_names('all', frame.columns.tolist(), '--traits')
    try:
        traits._parse_csv_names(','.join([*columns, columns[-1], columns[0]]), columns, '--traits')
    except ValueError as exc:
        trait_error = str(exc)
    else:
        raise AssertionError('Duplicate trait selection accepted')
    timings['trait_read_and_select'] = time.perf_counter() - started
    started = time.perf_counter()
    tree, names = copy_number.load_species_tree(str(tree_path))
    try:
        copy_number.load_species_tree(str(duplicate_tree))
    except SystemExit as exc:
        tree_error = str(exc)
    else:
        raise AssertionError('Duplicate species accepted')
    timings['species_tree_validation'] = time.perf_counter() - started
    assert names == leaves and selected == columns
    assert contract_error == 'Duplicate input key(s): input_00000'
    assert trait_error == f'--traits contains duplicates: {columns[0]}, {columns[-1]}'
    assert tree_error == f'ERROR: Dated species tree has duplicate leaf label(s): {leaves[0]}, {leaves[-1]}'
    payload = {'contract': contract, 'trait_tsv': frame.to_csv(sep='\t', index=False),
               'selected': selected, 'tree': tree.write(), 'names': names,
               'diagnostics': [contract_error, trait_error, tree_error]}
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    print(json.dumps({'seconds': timings, 'output_sha256': digest,
                      'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--size', type=int, default=4096)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', type=Path)
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 32 <= args.size <= 8192:
        parser.error('--size must be 32..8192')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'size': args.size, 'python': sys.version,
               'platform': platform.platform(), 'samples': []}
    with tempfile.TemporaryDirectory(prefix='gg-declaration-benchmark-') as temporary:
        for trial in range(4):
            root = Path(temporary) / str(trial)
            root.mkdir()
            raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()),
                '--worker', str(root), '--support-root', str(args.support_root), '--size', str(args.size)], text=True)
            sample = json.loads(raw)
            if trial:
                results['samples'].append(sample)
            print('warmup' if trial == 0 else trial, sample['seconds'], flush=True)
    assert len({row['output_sha256'] for row in results['samples']}) == 1
    results['median_seconds'] = {key: statistics.median(row['seconds'][key] for row in results['samples'])
                               for key in results['samples'][0]['seconds']}
    results['median_peak_rss_kib_linux'] = statistics.median(row['peak_rss_kib_linux'] for row in results['samples'])
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    print('median', results['median_seconds'], results['median_peak_rss_kib_linux'], flush=True)


if __name__ == '__main__':
    main()
