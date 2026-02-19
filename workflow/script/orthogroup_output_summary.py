#!/usr/bin/env python
# coding: utf-8

import argparse
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy
import pandas


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir_og', metavar='PATH', type=str, required=True, help='')
    parser.add_argument('--genecount', metavar='PATH', type=str, required=True, help='')
    parser.add_argument('--out', metavar='PATH', type=str, required=True, help='')
    parser.add_argument('--ncpu', metavar='INT', default=1, type=int, help='Number of worker threads.')
    # Backward-compatible alias.
    parser.add_argument('--threads', dest='ncpu', metavar='INT', type=int, help=argparse.SUPPRESS)
    return parser


def _amas_columns():
    return [
        'No_of_taxa',
        'Alignment_length',
        'Total_matrix_cells',
        'Undetermined_characters',
        'Missing_percent',
        'No_variable_sites',
        'Parsimony_informative_sites',
        'GC_content',
    ]


def _read_amas_file(file_path, amas_cols):
    og_id = re.sub(r'\..*', '', os.path.basename(file_path))
    tmp = pandas.read_csv(file_path, sep='\t', header=0)
    return og_id, tmp.loc[0, amas_cols].values


def get_amas_stats(df, dir_amas, extension, ncpu):
    amas_cols = _amas_columns()
    amas_new_cols = [f'{col}_{extension}' for col in amas_cols]
    for ncol in amas_new_cols:
        if ncol not in df.columns:
            df.loc[:, ncol] = numpy.nan

    if not os.path.isdir(dir_amas):
        print(f'{extension}: {dir_amas} was not found. Skipping.', flush=True)
        return df

    is_prefilled = ~df[f'No_of_taxa_{extension}'].isna()
    prefilled_ogs = set(df.index[is_prefilled])
    files = sorted(os.listdir(dir_amas))
    files = [
        file
        for file in files
        if (file.startswith('OG') or file.startswith('HOG') or file.startswith('SP'))
        and re.sub(r'\..*', '', file) not in prefilled_ogs
    ]

    counter = 0
    if ncpu > 1 and len(files) > 1:
        max_workers = min(ncpu, len(files))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(_read_amas_file, os.path.join(dir_amas, file), amas_cols)
                for file in files
            ]
            for future in as_completed(futures):
                og_id, values = future.result()
                df.loc[og_id, amas_new_cols] = values
                counter += 1
    else:
        for file in files:
            file_path = os.path.join(dir_amas, file)
            og_id, values = _read_amas_file(file_path, amas_cols)
            df.loc[og_id, amas_new_cols] = values
            counter += 1

    idx_total = numpy.where(df.columns == 'Total')[0][0]
    idx_added = numpy.arange(idx_total + 1, df.columns.shape[0])
    original_cols = df.columns[numpy.arange(idx_total + 1)].tolist()
    sorted_amas_cols = df.columns[idx_added].sort_values().tolist()
    df = df.loc[:, original_cols + sorted_amas_cols]
    print(f'{extension}: {counter} AMAS results were appended.', flush=True)
    return df


def run(args):
    print('args:', vars(args))

    start = time.time()
    updated_genecount = args.genecount.replace('.tsv', '') + '.amas.tsv'
    if os.path.exists(updated_genecount):
        print('Updated --genecount file was detected. Reading.', flush=True)
        df = pandas.read_csv(updated_genecount, sep='\t', index_col=0, header=0)
    else:
        print('Updated --genecount file was not detected. Reading the original.', flush=True)
        df = pandas.read_csv(args.genecount, sep='\t', index_col=0, header=0)

    dir_amas = os.path.join(args.dir_og, 'amas.original')
    df = get_amas_stats(df, dir_amas, 'original', args.ncpu)
    dir_amas = os.path.join(args.dir_og, 'amas.cleaned')
    df = get_amas_stats(df, dir_amas, 'clean', args.ncpu)
    df.to_csv(updated_genecount, index=True, sep='\t')

    df.loc[:, 'SGE_TASK_ID'] = numpy.arange(df.shape[0]) + 1

    subdirs = os.listdir(args.dir_og)
    subdirs = [sd for sd in subdirs if os.path.isdir(os.path.join(args.dir_og, sd))]
    subdirs = sorted(subdirs)
    df = pandas.concat([df, pandas.DataFrame(data=0, index=df.index, columns=subdirs, dtype=int)], axis=1)
    for subdir in subdirs:
        subdir_path = os.path.join(args.dir_og, subdir)
        subdir_path = os.path.realpath(subdir_path)
        files = os.listdir(subdir_path)
        og_ids = [re.sub(r'\..*', '', f) for f in files if f.startswith('OG') or f.startswith('HOG') or f.startswith('SP')]
        df.loc[og_ids, subdir] = 1
        num_missing = (df.loc[:, subdir] == 0).sum()
        txt = 'Subdirectory {}: {:,} / {:,} files are missing.'
        print(txt.format(subdir, num_missing, df.shape[0]))

    print('Writing output file:', args.out, flush=True)
    df.to_csv(args.out, index=True, sep='\t')
    print('Done. Elapsed time: {:,} sec'.format(int(time.time() - start)))


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    args.ncpu = max(1, int(args.ncpu))
    run(args)


if __name__ == '__main__':
    main()
