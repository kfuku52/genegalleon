#!/usr/bin/env python
# coding: utf-8

import argparse
import datetime
import os
import sys
import time

import numpy
import pandas


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file_orthogroup_table', metavar='PATH', default='', type=str, help='', required=True)
    parser.add_argument('--mode', metavar='STR', default='hog2og', choices=['hog2og', 'sp2og'], type=str, help='', required=True)
    parser.add_argument('--dir_out', metavar='PATH', default='', type=str, help='', required=True)
    return parser


def run(args):
    start = time.time()
    print('Starting {} at {}'.format(sys.argv[0], datetime.datetime.now()))

    if os.path.exists(args.dir_out):
        print('--out_dir exists already: {}'.format(args.dir_out))
    else:
        print('Creating --out_dir: {}'.format(args.dir_out))
        os.makedirs(args.dir_out)

    df = pandas.read_csv(args.file_orthogroup_table, sep='\t', header=0, low_memory=False)
    if args.mode == 'hog2og':
        df = df.drop(['OG', 'Gene Tree Parent Clade'], axis=1)
        df.columns = df.columns.str.replace('^HOG$', 'Orthogroup', regex=True)
        pref = os.path.basename(args.file_orthogroup_table).replace('tsv', '')
        df['Orthogroup'] = df['Orthogroup'].str.replace(pref, '', regex=False)
    elif args.mode == 'sp2og':
        df.columns = df.columns.str.replace('^group_id$', 'Orthogroup', regex=True)
        df.columns = df.columns.str.replace(r'\..*', '', regex=True)
        df['Orthogroup'] = 'SP' + df['Orthogroup'].astype(str).str.zfill(7)
        species_cols = [c for c in df.columns if c != 'Orthogroup']
        df.loc[:, species_cols] = df.loc[:, species_cols].apply(
            lambda col: col.str.replace('*', '', regex=False).str.replace(',', ', ', regex=False)
        )

    df = df.fillna('')
    outpath_og = os.path.join(args.dir_out, 'Orthogroups.tsv')
    print('Writing: {}'.format(outpath_og))
    df.to_csv(outpath_og, sep='\t', index=False)

    sp_cols = df.columns[df.columns != 'Orthogroup']
    gc = pandas.DataFrame({'Orthogroup': df['Orthogroup']})
    if len(sp_cols) > 0:
        values = df.loc[:, sp_cols].to_numpy(dtype=str, copy=False)
        comma_counts = numpy.char.count(values, ',')
        gene_counts = numpy.where(values == '', 0, comma_counts + 1).astype(int, copy=False)
        gc.loc[:, sp_cols] = gene_counts
        gc['Total'] = gene_counts.sum(axis=1)
    else:
        gc['Total'] = 0
    print('Number of Orthogroups: {:,}'.format(gc.shape[0]), flush=True)
    print('Largest Orthogroup size: {:,}'.format(gc['Total'].max()), flush=True)
    print('Smallest Orthogroup size: {:,}'.format(gc['Total'].min()), flush=True)

    outpath_gc = os.path.join(args.dir_out, 'Orthogroups.GeneCount.tsv')
    print('Writing: {}'.format(outpath_gc))
    gc.to_csv(outpath_gc, sep='\t', index=False)

    txt_readme = 'The files in this directory were created by {}. '.format(sys.argv[0])
    txt_readme += 'The OrthoFinder outputs in the obsoleted "Orthogroups" directory were '
    txt_readme += 'mimicked by parsing {}\n'.format(args.file_orthogroup_table)

    outpath_readme = os.path.join(args.dir_out, 'README.txt')
    with open(outpath_readme, 'w') as f:
        f.write(txt_readme)

    for tmp_file in [f for f in os.listdir(os.getcwd()) if f.startswith('tmp.')]:
        print('Removing: {}'.format(tmp_file))
        os.remove(tmp_file)

    print('Ending {} at {}. Elapsed time: {:,} sec'.format(sys.argv[0], datetime.datetime.now(), int(time.time() - start)))


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
