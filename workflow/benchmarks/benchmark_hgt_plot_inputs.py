#!/usr/bin/env python3
"""Paired whole-TSV versus native HGT index benchmark with exact TSV equality.

Run through workflow/tests/run_in_runtime.sh. Each repetition builds a new
index and includes its full preparation cost; warm reuse is not the comparison.
The synthetic wide candidate table follows the existing public HGT columns.
"""
import argparse
import hashlib
import importlib.util
import json
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pandas

SUPPORT = Path(__file__).resolve().parents[1] / 'support'
sys.path.insert(0, str(SUPPORT))
import hgt_plot_inputs  # noqa: E402


def load(path):
    spec = importlib.util.spec_from_file_location('hgt_annotation_benchmark', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fixture(root, families, genes):
    branch = root / 'branch.tsv'
    gene = root / 'gene.tsv'
    pandas.DataFrame([dict(orthogroup=f'OG{f:05d}', branch_id=0, candidate_gene_count=genes)
                      for f in range(families)]).to_csv(branch, sep='\t', index=False)
    pandas.DataFrame([dict(orthogroup=f'OG{f:05d}', gene_id=f'g{f}_{g}', gene_taxon='Plant species',
                          candidate_branch_count=1, intron_supported=bool(g % 2),
                          annotation='ACGT' * 256) for f in range(families) for g in range(genes)]).to_csv(
                              gene, sep='\t', index=False)
    return branch, gene


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-ref', required=True)
    parser.add_argument('--families', type=int, default=32)
    parser.add_argument('--genes-per-family', type=int, default=1000)
    parser.add_argument('--mode', choices=['baseline', 'indexed'], required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        branch, gene = fixture(root, args.families, args.genes_per_family)
        old = root / 'baseline.py'
        old.write_bytes(subprocess.check_output(['git', '-c', 'safe.directory=' + str(SUPPORT.parents[1]), 'show', args.baseline_ref +
                                                 ':workflow/support/annotate_hgt_tree_plot.py']))
        module = load(old)
        database = root / 'index.sqlite3'
        started = time.monotonic()
        if args.mode == 'indexed':
            hgt_plot_inputs.prepare(database, branch, gene)
        digest = hashlib.sha256()
        for f in range(args.families):
            family = f'OG{f:05d}'
            stat = pandas.DataFrame([dict(branch_id=0, node_name=f'g{f}_{g}', so_event='L')
                                     for g in range(args.genes_per_family)])
            if args.mode == 'baseline':
                branches, genes = module.safe_read_tsv(branch), module.safe_read_tsv(gene)
            else:
                branches, genes, _ = hgt_plot_inputs.read_family(database, branch, gene, family)
            out = module.annotate_leaf_rows(module.annotate_branch_rows(stat, branches, family), genes, family)
            digest.update(out.to_csv(sep='\t', index=False).encode())
        result = dict(mode=args.mode, families=args.families, genes_per_family=args.genes_per_family,
                      source_bytes=branch.stat().st_size + gene.stat().st_size,
                      seconds=time.monotonic() - started, peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                      output_sha256=digest.hexdigest(), baseline_ref=args.baseline_ref)
        print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
