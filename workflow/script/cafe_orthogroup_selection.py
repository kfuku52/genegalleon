#!/usr/bin/env python
# coding: utf-8

import argparse
import datetime
import pandas
import numpy
import os
import sys
import time
import ete4


def load_tree(newick_or_path, parser=1):
    if isinstance(newick_or_path, str) and os.path.exists(newick_or_path):
        with open(newick_or_path, 'r', encoding='utf-8') as handle:
            newick_or_path = handle.read().strip()
    return ete4.PhyloTree(newick_or_path, parser=parser)

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

    genecount_df = pandas.read_csv(args.genecount, sep='\t')
    tree = load_tree(args.dated_species_tree, parser=1)
    leaf_names = list(tree.leaf_names())
    genecount_df['size_differentials'] = genecount_df[leaf_names].max(axis=1) - genecount_df[leaf_names].min(axis=1)

    # Plot size differentials for all orthogroups
    print(f'Plotting size differentials for all orthogroups: {args.output_dir}/size_differential_histogram.svg, {args.output_dir}/size_differential_histogram.pdf')
    os.makedirs(args.output_dir, exist_ok=True)
    plt = get_pyplot()
    fig, ax = plt.subplots(figsize=(4, 4))
    bins = numpy.arange(0, genecount_df['size_differentials'].max() + 10, 10)
    ax.hist(genecount_df['size_differentials'], bins=bins, color='black')
    if args.max_size_differential < genecount_df['size_differentials'].max():
        ax.axvline(x=args.max_size_differential, color='red', linestyle='--', linewidth=1, label=f'Size differential <= {args.max_size_differential}')
        ax.legend()
    ax.set_xlabel('Size differential')
    ax.set_ylabel('Number of orthogroups')
    ax.set_yscale('log')
    fig.savefig(os.path.join(args.output_dir, 'size_differential_histogram.svg'), format='svg', dpi=300)
    fig.savefig(os.path.join(args.output_dir, 'size_differential_histogram.pdf'), format='pdf', dpi=300)

    # Filter orthogroups by size differential and write CAFE input file
    print(f'Removing orthogroups with size differentials > {args.max_size_differential}: {args.output_dir}/removed_orthogroups.tsv')
    genecount_df[genecount_df['size_differentials'] > args.max_size_differential].to_csv(os.path.join(args.output_dir, 'removed_orthogroups.tsv'), sep='\t', index=False)
    genecount_df = genecount_df[genecount_df['size_differentials'] <= args.max_size_differential]

    # Write CAFE input file
    print(f'Writing CAFE input file: {args.output_dir}/cafe_input.tsv')
    cafe_input_df = genecount_df[['besthit_0.95', 'Orthogroup'] + leaf_names]
    cafe_input_df.to_csv(os.path.join(args.output_dir, 'cafe_input.tsv'), sep='\t', index=False)

    print('Ending {} at {}. Elapsed time: {:,} sec'.format(sys.argv[0], datetime.datetime.now(), int(time.time() - start)))


if __name__ == '__main__':
    main()
