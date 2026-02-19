#!/usr/bin/env python
# coding: utf-8

import argparse
import datetime
import os
import pandas
import sys
import time

pandas.set_option("display.max_columns", None)

def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--orthofinder_og_genecount', metavar='PATH', default='', type=str, help='', required=True)
    parser.add_argument('--orthofinder_hog_genecount', metavar='PATH', default='', type=str, help='', required=True)
    parser.add_argument('--sonicparanoid_genecount', metavar='PATH', default='', type=str, help='', required=True)
    return parser


def get_pyplot():
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.rcParams['font.size'] = 8
    matplotlib.rcParams['font.family'] = 'Helvetica'
    matplotlib.rcParams['svg.fonttype'] = 'none' # none, path, or svgfont
    return plt


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    start = time.time()
    cwd = os.getcwd()
    print('Working at: {}'.format(cwd))
    print('Starting {} at {}'.format(sys.argv[0], datetime.datetime.now()))

    dfs = {}
    dfs['OrthoFinder Orthogroup'] = pandas.read_csv(args.orthofinder_og_genecount, sep='\t', header=0, low_memory=False)
    dfs['OrthoFinder Hierarchical orthogroup'] = pandas.read_csv(args.orthofinder_hog_genecount, sep='\t', header=0, low_memory=False)
    dfs['SonicParanoid'] = pandas.read_csv(args.sonicparanoid_genecount, sep='\t', header=0, low_memory=False)
    for key in dfs.keys():
        dfs[key] = dfs[key].sort_values(by='Total')
        dfs[key]['cumulative_num_gene'] = dfs[key]['Total'].cumsum()

    colors = {
        'OrthoFinder Orthogroup':'#991574',
        'OrthoFinder Hierarchical orthogroup':'#749915',
        'SonicParanoid':'#157499',
    }

    plt = get_pyplot()
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(7.2, 6.4), sharey=False, sharex=False)
    axes = axes.flat

    max_ngenes = max([dfs[key]['Total'].max() for key in dfs.keys()])
    xmax = min(max_ngenes, 1000)
    bins = range(0, xmax, 1)

    ax = axes[0]
    for key in dfs.keys():
        ax.hist(x=dfs[key]['Total'], bins=bins, color=colors[key], histtype='step', cumulative=True, label=key, linewidth=1.0, alpha=0.5)
    ax.set_xlim(0, xmax)
    ax.set_xlabel('Number of genes per orthogroup')
    ax.set_ylabel('Cumulative number of orthogroups')
    ax.legend(loc='lower right')

    ax = axes[1]
    for key in dfs.keys():
        ax.hist(x=dfs[key]['Total'], bins=bins, color=colors[key], histtype='step', label=key, linewidth=1.0, alpha=0.5)
    ax.set_xlim(0, xmax)
    ax.set_xlabel('Number of genes per orthogroup')
    ax.set_ylabel('Number of orthogroups')
    ax.legend(loc='upper right')

    ax = axes[2]
    for key in dfs.keys():
        ax.plot(dfs[key]['Total'], dfs[key]['cumulative_num_gene'], color=colors[key], label=key, linewidth=1.0, alpha=0.5)
    ax.set_xlim(0, xmax)
    ax.set_xlabel('Number of genes per orthogroup')
    ax.set_ylabel('Cumulative number of genes')
    ax.legend(loc='lower right')

    ax = axes[3]
    for key in dfs.keys():
        tmp = dfs[key].copy()
        tmp['Total2'] = tmp['Total'].copy()
        tmp = tmp.groupby('Total2')['Total'].sum().reset_index()
        tmp2 = pandas.DataFrame({'Total2': range(0, xmax+1)})
        tmp = pandas.merge(tmp2, tmp, on='Total2', how='left')
        tmp['Total'] = tmp['Total'].fillna(0)
        ax.plot(tmp['Total2'], tmp['Total'], color=colors[key], label=key, linewidth=1.0, alpha=0.5)
    ax.set_xlim(0, xmax)
    ax.set_xlabel('Number of genes per orthogroup')
    ax.set_ylabel('Number of genes')
    ax.legend(loc='upper right')

    outbase = 'orthogroup_histogram'
    fig.tight_layout(pad=0.25, w_pad=1.5, h_pad=0.5)
    for ext in ['svg', 'pdf']:
        outpath = os.path.join(outbase + '.' + ext)
        fig.savefig(outpath, format=ext)

    print('Ending {} at {}. Elapsed time: {:,} sec'.format(sys.argv[0], datetime.datetime.now(), int(time.time() - start)))


if __name__ == '__main__':
    main()
