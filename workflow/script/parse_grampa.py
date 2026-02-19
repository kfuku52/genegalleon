#!/usr/bin/env python
# coding: utf-8

import argparse
import datetime
from concurrent.futures import ProcessPoolExecutor
import pandas
import re
import sys
import os
import ete4


def load_tree(newick_or_path, parser=0):
    if isinstance(newick_or_path, str) and os.path.exists(newick_or_path):
        with open(newick_or_path, 'r', encoding='utf-8') as handle:
            newick_or_path = handle.read().strip()
    return ete4.PhyloTree(newick_or_path, parser=parser)

_WORKER_SPECIES_NAMES = None


def _init_worker(species_names):
    global _WORKER_SPECIES_NAMES
    _WORKER_SPECIES_NAMES = species_names


def summarize_gene_tree(task, species_names=None):
    idx, gt_txt = task
    if species_names is None:
        species_names = _WORKER_SPECIES_NAMES
    gt_id = 'GT-' + str(idx + 1)
    gt = load_tree(newick_or_path=gt_txt, parser=0)
    gene_names = list(gt.leaf_names())
    species_set = set(species_names)
    species_gene_lists = {species_name: [] for species_name in species_names}
    for gene_name in gene_names:
        matched_species = None
        gene_id = None
        parts = gene_name.rsplit('_', 2)
        if len(parts) >= 3:
            candidate_species = parts[-2] + '_' + parts[-1]
            if candidate_species in species_set:
                matched_species = candidate_species
                gene_id = gene_name[:-(len(candidate_species) + 1)]
        if matched_species is None:
            for species_name in species_names:
                suffix = '_' + species_name
                if gene_name.endswith(suffix):
                    matched_species = species_name
                    gene_id = gene_name[:-len(suffix)]
                    break
        if matched_species is not None:
            species_gene_lists[matched_species].append(gene_id)
    species_gene_map = {species_name: ','.join(species_gene_lists[species_name]) for species_name in species_names}
    return gt_id, species_gene_map


def update_det_with_species_genes(det, gt_id, species_names, species_gene_map, row_indices_by_gt):
    row_indices = row_indices_by_gt.get(gt_id)
    if row_indices is None:
        return
    det.loc[row_indices, species_names] = [species_gene_map[species_name] for species_name in species_names]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--grampa_det', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--grampa_out', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--gene_trees', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--species_tree', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--sorted_gene_tree_file_names', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--ncpu', metavar='INT', default=1, type=int, help='Number of CPU processes.')
    # Backward-compatible alias.
    parser.add_argument('--threads', dest='ncpu', metavar='INT', type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.ncpu = max(1, int(args.ncpu))
    print('Starting parse_grampa.py: {}'.format(datetime.datetime.now()))
    outfile = 'grampa_summary.tsv'

    print('Processing grampa det file')
    det = pandas.read_csv(args.grampa_det, sep='\t', header=0, low_memory=False)
    det.columns = det.columns.str.replace('^# ', '', regex=True)
    det.columns = det.columns.str.replace('Total score', 'total_score', regex=False)
    det.columns = det.columns.str.replace('Maps', 'maps', regex=False)
    det = det.loc[(det['GT/MT combo'].str.startswith('*')), :]
    det = det.loc[(~det['dups'].astype(str).str.endswith('maps found!')), :].reset_index(drop=True)
    det['gene_tree'] = det['GT/MT combo'].str.replace('* ', '', regex=False).str.replace(' to.*', '', regex=True)
    det['mul_tree'] = det['GT/MT combo'].str.replace('.* to ', '', regex=True)
    det = det.loc[:, ['gene_tree', 'mul_tree', 'dups', 'losses', 'total_score', 'maps']]
    print('{} Grampa maps were found.'.format(det.shape[0]))
    if det.shape[0] == 0:
        print('No Grampa maps were dound. A placeholder output file will be generated.')
        with open(outfile, 'w') as f:
            f.write('This is a placeholder file. No Grampa maps were found.')
        print('Ending parse_grampa.py: {}'.format(datetime.datetime.now()))
        sys.exit(0)

    print('Processing grampa input species trees')
    st = load_tree(newick_or_path=args.species_tree, parser=0)
    species_names = sorted(list(st.leaf_names()))
    det.loc[:, species_names] = ''
    print('{} species were found.'.format(len(species_names)))

    print('Processing grampa input gene trees')
    with open(args.gene_trees, 'r') as f:
        gt_txts = f.readlines()

    row_indices_by_gt = det.groupby('gene_tree').indices
    tasks = [(i, gt_txt) for i, gt_txt in enumerate(gt_txts)]
    if args.ncpu > 1 and len(tasks) > 1:
        max_workers = min(args.ncpu, len(tasks))
        with ProcessPoolExecutor(
            max_workers=max_workers,
            initializer=_init_worker,
            initargs=(species_names,),
        ) as executor:
            for gt_id, species_gene_map in executor.map(summarize_gene_tree, tasks):
                update_det_with_species_genes(det, gt_id, species_names, species_gene_map, row_indices_by_gt)
    else:
        for task in tasks:
            gt_id, species_gene_map = summarize_gene_tree(task, species_names=species_names)
            update_det_with_species_genes(det, gt_id, species_names, species_gene_map, row_indices_by_gt)

    print('Processing grampa out file')
    with open(args.grampa_out, 'r') as f:
        out_txts = f.readlines()
    out_txts = [l for l in out_txts if l.startswith('MT-')]
    out_txts = [re.split('\t', l) for l in out_txts]
    out = pandas.DataFrame(out_txts, columns=['mul_tree', 'H1_node', 'H2_node', 'mul_tree_string', 'multree_score'])
    out['multree_score'] = out['multree_score'].replace('\n', '', regex=False).astype(int)
    print('{} MUL trees were found.'.format(out.shape[0]))

    print('Adding the original file names of gene trees')
    gtname = pandas.read_csv(args.sorted_gene_tree_file_names, sep='\t', header=None)
    gtname.columns = ['file_name']
    gtname['gene_tree'] = 'GT-' + pandas.Series([str(i + 1) for i in range(gtname.shape[0])])

    print('Writing output table')
    df = pandas.merge(det, out, on='mul_tree', how='left')
    df = pandas.merge(gtname, df, on='gene_tree', how='right')
    df.to_csv(outfile, sep='\t', index=False, quoting=None)

    print('Ending parse_grampa.py: {}'.format(datetime.datetime.now()))


if __name__ == '__main__':
    main()
