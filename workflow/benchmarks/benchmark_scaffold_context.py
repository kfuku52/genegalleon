#!/usr/bin/env python3
"""Measure complete scaffold-context attachment using per-species local tables.

Optionally compare a baseline support directory in alternating fresh processes.
Fixture generation and output serialization are outside the operation timer;
Linux VmHWM is the worker's peak through attach_context, including fixtures.
"""
import argparse
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
    import pandas as pd
    import scaffold_taxonomy as scaffold

    root = args.worker
    candidates = []
    for species in ('Host_alpha', 'Host_beta'):
        rows = []
        for i in range(args.genes):
            gene = f'g{i:05d}'
            unit = 'gff_locus' if i % 4 < 2 else 'cds_id'
            locus = f'L{i // 2}' if unit == 'gff_locus' else gene
            for j, rank in enumerate(scaffold.RANKS):
                rows.append({'species': species, 'gene_id': gene, 'scaffold': f's{i // 64}',
                             'locus_id': locus, 'count_unit': unit, 'rank': rank, 'host_taxid': '3',
                             'label': ('compatible', 'incompatible', 'unresolved')[(i // 2 + j) % 3]})
            if i % 2 == 0:
                candidates.append({'orthogroup': f'OG{i % 32:02d}', 'gene_id': gene,
                                   'gene_taxon': species.replace('_', ' ')})
        pd.DataFrame(rows).to_csv(root / f'{species}_gene_taxonomy.tsv', sep='\t', index=False)
    candidates.extend([{'orthogroup': 'OG00', 'gene_id': 'absent', 'gene_taxon': 'Host alpha'},
                       {'orthogroup': 'OG00', 'gene_id': 'unknown', 'gene_taxon': 'Unknown taxon'}])
    genes = pd.DataFrame(candidates)
    genes.index = pd.Index([3 * i + 7 for i in range(len(genes))], name='input_gene')
    branches = []
    for orthogroup, group in genes.groupby('orthogroup', sort=True):
        for recipient in ('Host_alpha', 'recipient'):
            branches.append({'orthogroup': orthogroup, 'generax_transfer': f'Y@Donor@{recipient}',
                             'candidate_genes': '; '.join(dict.fromkeys(group.gene_id))})
    branches.append({'orthogroup': 'OG00', 'generax_transfer': 'Y', 'candidate_genes': 'absent'})
    branches = pd.DataFrame(branches)
    branches = pd.concat([branches] * args.branch_repeats, ignore_index=True)
    branches.index = pd.Index([5 * i + 11 for i in range(len(branches))], name='input_branch')
    tree = root / 'species.nwk'
    tree.write_text('((Host_alpha:1,Host_beta:1)recipient:1,Donor:1)root;\n')
    started = time.perf_counter()
    started_cpu = time.process_time()
    branch_out, gene_out = scaffold.attach_context(branches, genes, root, tree)
    seconds = time.perf_counter() - started
    cpu_seconds = time.process_time() - started_cpu
    peak_rss = next(int(line.split()[1]) for line in Path('/proc/self/status').read_text().splitlines()
                    if line.startswith('VmHWM:')) if sys.platform == 'linux' else resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    payload = bytearray()
    for frame in (branch_out, gene_out):
        payload.extend(json.dumps({'columns': frame.columns.tolist(), 'dtypes': [str(x) for x in frame.dtypes],
                                  'index_name': frame.index.name}, sort_keys=True).encode())
        payload.extend(frame.to_csv(sep='\t', index=True).encode())
    if args.worker_output:
        args.worker_output.write_bytes(payload)
    print(json.dumps({'seconds': seconds, 'cpu_seconds': cpu_seconds,
                      'output_sha256': hashlib.sha256(payload).hexdigest(), 'genes': len(gene_out),
                      'branches': len(branch_out), 'peak_rss_kib': peak_rss, 'peak_rss_kib_linux': peak_rss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--genes', type=int, default=2048, help='Taxonomy genes per species (two species)')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', type=Path)
    parser.add_argument('--worker-output', type=Path)
    parser.add_argument('--baseline-support', type=Path)
    parser.add_argument('--branch-repeats', type=int, default=1)
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 64 <= args.genes <= 8192 or args.genes % 64:
        parser.error('--genes must be a multiple of 64 in 64..8192')
    if args.branch_repeats < 1 or args.repeats < 1:
        parser.error('--branch-repeats and --repeats must be positive')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    roots = {'candidate': args.support_root}
    if args.baseline_support:
        roots = {'baseline': args.baseline_support.resolve(), **roots}
    results = {'support_root': str(args.support_root), 'genes_per_species': args.genes,
               'branch_repeats': args.branch_repeats, 'python': sys.version,
               'platform': platform.platform(), 'arguments': sys.argv, 'implementations': {
                   name: {'root': str(root), 'sha256': hashlib.sha256((root / 'scaffold_taxonomy.py').read_bytes()).hexdigest()}
                   for name, root in roots.items()}, 'paired_samples': {name: [] for name in roots}}
    reference = None
    with tempfile.TemporaryDirectory(prefix='gg-scaffold-context-benchmark-') as temporary:
        for trial in range(args.repeats + 1):
            order = list(roots) if trial % 2 == 0 else list(reversed(roots))
            for name in order:
                root = Path(temporary) / f'{name}-{trial}'
                root.mkdir()
                output = root / 'output'
                raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()), '--worker', str(root),
                    '--support-root', str(roots[name]), '--genes', str(args.genes), '--branch-repeats',
                    str(args.branch_repeats), '--worker-output', str(output)], text=True)
                sample = json.loads(raw)
                payload = output.read_bytes()
                if reference is None:
                    reference = payload
                else:
                    assert payload == reference, (name, trial)
                if trial:
                    results['paired_samples'][name].append(sample)
                print(name, 'warmup' if trial == 0 else trial, sample, flush=True)
    results['medians'] = {name: {metric: statistics.median(row[metric] for row in samples)
                               for metric in ('seconds', 'cpu_seconds', 'peak_rss_kib')}
                          for name, samples in results['paired_samples'].items()}
    # Retain the original single-implementation result fields for saved readers.
    results['samples'] = results['paired_samples']['candidate']
    results['median_seconds'] = results['medians']['candidate']['seconds']
    results['median_peak_rss_kib_linux'] = results['medians']['candidate']['peak_rss_kib']
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    print('medians', results['medians'])


if __name__ == '__main__':
    main()
