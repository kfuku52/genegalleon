#!/usr/bin/env python3

import argparse
import datetime
import os
import sys
import time

import numpy
import pandas

pandas.options.mode.chained_assignment = None
pandas.set_option('display.max_columns', None)


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--min_scaffold_size', metavar='INT', default=1000000, type=int, help='')
    parser.add_argument('--fx2tab', metavar='PATH', default='', type=str, help='')
    return parser


def get_matplotlib():
    import matplotlib

    matplotlib.rcParams['font.size'] = 8
    matplotlib.rcParams['font.family'] = 'Helvetica'
    matplotlib.rcParams['svg.fonttype'] = 'none'  # none, path, or svgfont
    return matplotlib


def add_size_legend(ax, matplotlib, label_large, label_small):
    handles = [
        matplotlib.lines.Line2D([0], [0], marker='o', linestyle='', color='black', markersize=4, label=label_large),
        matplotlib.lines.Line2D([0], [0], marker='o', linestyle='', color='grey', markersize=4, label=label_small),
    ]
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.05), shadow=True, ncol=1, frameon=False)


def scatter_by_size(ax, df, ycol, label_large, label_small, rng):
    for label, color, xoffset in ((label_large, 'black', -0.035), (label_small, 'grey', 0.035)):
        sub = df.loc[(df['size'] == label), ycol]
        if sub.shape[0] == 0:
            continue
        jitter = rng.normal(loc=0.0, scale=0.01, size=sub.shape[0])
        xpos = xoffset + jitter
        ax.scatter(xpos, sub.values, alpha=0.75, s=8, color=color, edgecolors='none')


def run(args):
    sys.stderr.write('{} started: {}\n'.format(sys.argv[0], datetime.datetime.now()))

    label_large = '≥{:.1E} bp'.format(args.min_scaffold_size)
    label_small = '<{:.1E} bp'.format(args.min_scaffold_size)

    df = pandas.read_csv(args.fx2tab, sep='\t', header=0, low_memory=False)
    df['size'] = label_large
    df.loc[(df['length'] < args.min_scaffold_size), 'size'] = label_small
    df['color'] = 'black'
    df['log10_length'] = numpy.log10(df['length'])

    df2 = df.loc[(df['length'] >= args.min_scaffold_size), :].reset_index(drop=True)
    sm = df.loc[(df['length'] < args.min_scaffold_size), :].reset_index(drop=True)
    if sm.shape[0] > 0:
        merge_id = '{}\n({:,} sequences)'.format(label_small, sm.shape[0])
        merge_length = sm['length'].sum()
        merge_gc = ((sm['GC'] * sm['length']) / merge_length).sum()
        merge_gcskew = ((sm['GC-Skew'] * sm['length']) / merge_length).sum()
        new_row = {
            '#id': [merge_id],
            'length': [merge_length],
            'GC': [merge_gc],
            'GC-Skew': [merge_gcskew],
            'color': 'grey',
            'size': label_small,
        }
        sm2 = pandas.DataFrame(new_row)
        df2 = pandas.concat([df2, sm2], axis=0, ignore_index=True).reset_index(drop=True)

    skew_max = df['GC-Skew'].abs().max() * 1.1

    matplotlib = get_matplotlib()
    rng = numpy.random.default_rng(11)
    fig = matplotlib.pyplot.figure(figsize=(7.2, 4.8))

    ax = fig.add_subplot(2, 5, (1, 5))
    xpos = numpy.arange(df2.shape[0])
    ax.bar(xpos, df2['length'].values, color=df2['color'].values, width=0.8)
    ax.set_ylabel('Scaffold size (bp)')
    ax.set_xlabel('')
    ax.set_xticks(xpos)
    ax.set_xticklabels(df2['#id'].values, rotation=90)

    ax = fig.add_subplot(2, 5, 6)
    scatter_by_size(ax=ax, df=df, ycol='length', label_large=label_large, label_small=label_small, rng=rng)
    ax.set_ylabel('Scaffold size (bp)')
    ax.set_xlabel('')
    ax.set_xticks(ticks=[0], labels=[''])
    ax.set_xlim(-0.08, 0.08)
    ax.set_ylim(0, df['length'].max() * 1.1)
    ax.axhline(y=args.min_scaffold_size, linestyle='dashed', color='black', linewidth=0.5)
    add_size_legend(ax=ax, matplotlib=matplotlib, label_large=label_large, label_small=label_small)

    ax = fig.add_subplot(2, 5, 7)
    scatter_by_size(ax=ax, df=df, ycol='GC', label_large=label_large, label_small=label_small, rng=rng)
    ax.set_ylabel('GC content (%)')
    ax.set_xlabel('')
    ax.set_xticks(ticks=[0], labels=[''])
    ax.set_xlim(-0.08, 0.08)
    ax.set_ylim(0, 100)
    ax.axhline(y=50, linestyle='dashed', color='black', linewidth=0.5)
    add_size_legend(ax=ax, matplotlib=matplotlib, label_large=label_large, label_small=label_small)

    ax = fig.add_subplot(2, 5, 8)
    scatter_by_size(ax=ax, df=df, ycol='GC-Skew', label_large=label_large, label_small=label_small, rng=rng)
    ax.set_ylabel('GC skew')
    ax.set_xlabel('')
    ax.set_xticks(ticks=[0], labels=[''])
    ax.set_xlim(-0.08, 0.08)
    ax.set_ylim(-skew_max, skew_max)
    ax.axhline(y=0, linestyle='dashed', color='black', linewidth=0.5)
    add_size_legend(ax=ax, matplotlib=matplotlib, label_large=label_large, label_small=label_small)

    ax = fig.add_subplot(2, 5, (9, 10))
    for label, color in ((label_large, 'black'), (label_small, 'grey')):
        sub = df.loc[(df['size'] == label), :]
        if sub.shape[0] == 0:
            continue
        ax.scatter(sub['log10_length'].values, sub['GC'].values, alpha=0.75, s=8, color=color, label=label)
    ax.axvline(x=numpy.log10(args.min_scaffold_size), linestyle='dashed', color='black', linewidth=0.5)
    ax.set_ylabel('GC content (%)')
    ax.set_xlabel('Scaffold length (log10 bp)')
    ax.set_ylim(0, 100)
    ax.legend(shadow=True, ncol=1, frameon=False)

    fig.tight_layout(pad=0.25, w_pad=1.0, h_pad=1.0)
    outbase = 'scaffold_size_histogram'
    for ext in ['svg', 'pdf']:
        outpath = os.path.join(outbase + '.' + ext)
        fig.savefig(outpath, format=ext)

    sys.stderr.write('{} ended: {}\n'.format(sys.argv[0], datetime.datetime.now()))


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    args.min_scaffold_size = max(1, int(args.min_scaffold_size))
    run(args)


if __name__ == '__main__':
    main()
