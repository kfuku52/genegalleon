#!/usr/bin/env python3
# coding: utf-8

import argparse
import datetime
import os
import sys
import time

import ete4
import numpy
import pandas

REQUIRED_GENECOUNT_COLUMNS = ['besthit_0.95', 'Orthogroup']


def fail(message):
    raise SystemExit(f'ERROR: {message}')


def require_input_file(path, label):
    if not path:
        fail(f'{label} path is empty.')
    if not os.path.isfile(path):
        fail(f'{label} was not found: {path}')
    if os.path.getsize(path) == 0:
        fail(f'{label} is empty: {path}')


def load_tree(newick_or_path, parser=1):
    if isinstance(newick_or_path, str) and os.path.exists(newick_or_path):
        with open(newick_or_path, 'r', encoding='utf-8') as handle:
            newick_or_path = handle.read().strip()
    return ete4.PhyloTree(newick_or_path, parser=parser)


def load_gene_count_table(path):
    require_input_file(path, 'Orthogroup gene-count table')
    try:
        genecount_df = pandas.read_csv(path, sep='\t')
    except Exception as exc:
        fail(f'Could not read orthogroup gene-count table {path}: {exc}')
    if genecount_df.empty:
        fail(f'Orthogroup gene-count table has no rows: {path}')
    missing_columns = [col for col in REQUIRED_GENECOUNT_COLUMNS if col not in genecount_df.columns]
    if missing_columns:
        fail(
            'Orthogroup gene-count table is missing required column(s): '
            + ', '.join(missing_columns)
        )
    orthogroup_ids = genecount_df['Orthogroup'].astype(str).str.strip()
    empty_ids = genecount_df['Orthogroup'].isna() | (orthogroup_ids.str.len() == 0)
    if empty_ids.any():
        fail('Orthogroup gene-count table contains empty Orthogroup IDs.')
    duplicated_ids = genecount_df.loc[genecount_df['Orthogroup'].duplicated(), 'Orthogroup'].astype(str).unique()
    if len(duplicated_ids):
        fail('Orthogroup gene-count table contains duplicate Orthogroup IDs: ' + ', '.join(duplicated_ids[:20]))
    return genecount_df


def load_species_tree(path):
    require_input_file(path, 'Dated species tree')
    try:
        tree = load_tree(path, parser=1)
    except Exception as exc:
        fail(f'Could not read dated species tree {path}: {exc}')
    leaf_names = list(tree.leaf_names())
    if not leaf_names:
        fail(f'Dated species tree has no leaf labels: {path}')
    duplicated_leaf_names = sorted({leaf for leaf in leaf_names if leaf_names.count(leaf) > 1})
    if duplicated_leaf_names:
        fail('Dated species tree has duplicate leaf label(s): ' + ', '.join(duplicated_leaf_names[:20]))
    return tree, leaf_names


def get_species_count_table(genecount_df, leaf_names):
    missing_species = [leaf for leaf in leaf_names if leaf not in genecount_df.columns]
    if missing_species:
        fail(
            'Orthogroup gene-count table is missing species column(s) from the dated species tree: '
            + ', '.join(missing_species[:20])
        )
    species_counts = genecount_df.loc[:, leaf_names].apply(pandas.to_numeric, errors='coerce')
    invalid_cols = [
        col for col in leaf_names
        if species_counts[col].isna().any()
    ]
    if invalid_cols:
        fail(
            'Orthogroup gene-count table contains non-numeric copy numbers in species column(s): '
            + ', '.join(invalid_cols[:20])
        )
    if (species_counts < 0).any().any():
        fail('Orthogroup gene-count table contains negative copy numbers.')
    return species_counts

def get_pyplot():
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.rcParams['font.size'] = 8
    matplotlib.rcParams['font.family'] = 'Helvetica'
    matplotlib.rcParams['svg.fonttype'] = 'none' # none, path, or svgfont
    return plt

def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genecount', metavar='PATH', default='', type=str, help='', required=True)
    parser.add_argument('--dated_species_tree', metavar='PATH', default='', type=str, help='', required=True)
    parser.add_argument('--output_dir', metavar='PATH', default='', type=str, help='', required=True)
    parser.add_argument('--max_size_differential', metavar='INT', default=int(1e9), type=int, help='')
    return parser

def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    start = time.time()
    print('Starting {} at {}'.format(sys.argv[0], datetime.datetime.now()))

    genecount_df = load_gene_count_table(args.genecount)
    _tree, leaf_names = load_species_tree(args.dated_species_tree)
    species_counts = get_species_count_table(genecount_df, leaf_names)
    genecount_df['size_differentials'] = species_counts.max(axis=1) - species_counts.min(axis=1)

    # Plot size differentials for all orthogroups
    print(f'Plotting size differentials for all orthogroups: {args.output_dir}/size_differential_histogram.svg, {args.output_dir}/size_differential_histogram.pdf')
    os.makedirs(args.output_dir, exist_ok=True)
    plt = get_pyplot()
    fig, ax = plt.subplots(figsize=(4, 4))
    max_size_differential = genecount_df['size_differentials'].max()
    bin_stop = max(10, int(numpy.ceil(max_size_differential / 10.0)) * 10 + 10)
    bins = numpy.arange(0, bin_stop + 1, 10)
    ax.hist(genecount_df['size_differentials'], bins=bins, color='black')
    if args.max_size_differential < genecount_df['size_differentials'].max():
        ax.axvline(x=args.max_size_differential, color='red', linestyle='--', linewidth=1, label=f'Size differential <= {args.max_size_differential}')
        ax.legend()
    ax.set_xlabel('Size differential')
    ax.set_ylabel('Number of orthogroups')
    ax.set_yscale('log')
    fig.savefig(os.path.join(args.output_dir, 'size_differential_histogram.svg'), format='svg', dpi=300)
    fig.savefig(os.path.join(args.output_dir, 'size_differential_histogram.pdf'), format='pdf', dpi=300)

    # Filter orthogroups by size differential and write orthogroup copy-number table
    print(f'Removing orthogroups with size differentials > {args.max_size_differential}: {args.output_dir}/removed_orthogroups.tsv')
    genecount_df[genecount_df['size_differentials'] > args.max_size_differential].to_csv(os.path.join(args.output_dir, 'removed_orthogroups.tsv'), sep='\t', index=False)
    genecount_df = genecount_df[genecount_df['size_differentials'] <= args.max_size_differential]

    # Write orthogroup copy-number table
    print(f'Writing orthogroup copy-number table: {args.output_dir}/orthogroup_copy_number.tsv')
    copy_number_df = pandas.concat(
        [
            genecount_df[['besthit_0.95', 'Orthogroup']].reset_index(drop=True),
            species_counts.loc[genecount_df.index].reset_index(drop=True),
        ],
        axis=1,
    )
    copy_number_df.to_csv(os.path.join(args.output_dir, 'orthogroup_copy_number.tsv'), sep='\t', index=False)

    print('Ending {} at {}. Elapsed time: {:,} sec'.format(sys.argv[0], datetime.datetime.now(), int(time.time() - start)))


if __name__ == '__main__':
    main()
