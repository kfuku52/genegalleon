#!/usr/bin/env python3
"""Measure complete scaffold-context attachment using per-species local tables."""
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
    branches.index = pd.Index([5 * i + 11 for i in range(len(branches))], name='input_branch')
    tree = root / 'species.nwk'
    tree.write_text('((Host_alpha:1,Host_beta:1)recipient:1,Donor:1)root;\n')
    started = time.perf_counter()
    branch_out, gene_out = scaffold.attach_context(branches, genes, root, tree)
    seconds = time.perf_counter() - started
    digest = hashlib.sha256()
    for frame in (branch_out, gene_out):
        digest.update(json.dumps({'columns': frame.columns.tolist(), 'dtypes': [str(x) for x in frame.dtypes],
                                  'index_name': frame.index.name}, sort_keys=True).encode())
        digest.update(frame.to_csv(sep='\t', index=True).encode())
    print(json.dumps({'seconds': seconds, 'output_sha256': digest.hexdigest(), 'genes': len(gene_out),
                      'branches': len(branch_out), 'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--genes', type=int, default=2048, help='Taxonomy genes per species (two species)')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', type=Path)
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 64 <= args.genes <= 8192 or args.genes % 64:
        parser.error('--genes must be a multiple of 64 in 64..8192')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    results = {'support_root': str(args.support_root), 'genes_per_species': args.genes, 'python': sys.version,
               'platform': platform.platform(), 'samples': []}
    with tempfile.TemporaryDirectory(prefix='gg-scaffold-context-benchmark-') as temporary:
        for trial in range(4):
            root = Path(temporary) / str(trial)
            root.mkdir()
            raw = subprocess.check_output([sys.executable, str(Path(__file__).resolve()), '--worker', str(root),
                '--support-root', str(args.support_root), '--genes', str(args.genes)], text=True)
            sample = json.loads(raw)
            if trial:
                results['samples'].append(sample)
            print('warmup' if trial == 0 else trial, sample['seconds'], flush=True)
    assert len({row['output_sha256'] for row in results['samples']}) == 1
    for metric in ('seconds', 'peak_rss_kib_linux'):
        results['median_' + metric] = statistics.median(row[metric] for row in results['samples'])
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    print('medians', results['median_seconds'], results['median_peak_rss_kib_linux'])


if __name__ == '__main__':
    main()
