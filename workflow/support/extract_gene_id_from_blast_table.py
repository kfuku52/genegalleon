#!/usr/bin/env python3

import argparse
import datetime
import sys
import time

import numpy
import pandas

pandas.options.mode.chained_assignment = None

def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--outfile', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--min_query_blast_coverage', metavar='FLOAT', default=0.25, type=float, help='')
    parser.add_argument('--max_num_gene_blast_hit_retrieval', metavar='INT', default=5000, type=int, help='')
    return parser


def parse_min_evalue(value):
    if pandas.isna(value):
        return pandas.NA
    text = str(value).strip()
    if text == '':
        return pandas.NA
    candidates = []
    for token in text.split(';'):
        token = token.strip()
        if token == '':
            continue
        try:
            candidates.append(float(token))
        except ValueError:
            continue
    if len(candidates) == 0:
        return pandas.NA
    return min(candidates)


def read_required_columns(path):
    header = pandas.read_csv(path, sep="\t", nrows=0)
    columns = header.columns.tolist()
    base_columns = ["sacc", "qjointcov"]
    missing_base = [c for c in base_columns if c not in columns]
    if missing_base:
        raise ValueError(f"Required column(s) missing: {', '.join(missing_base)}")

    if "min_evalue" in columns:
        usecols = base_columns + ["min_evalue"]
        evalue_source = "min_evalue"
    elif "evalue" in columns:
        usecols = base_columns + ["evalue"]
        evalue_source = "evalue"
    else:
        raise ValueError("Required column 'evalue' (or 'min_evalue') is missing.")

    df = pandas.read_csv(path, sep="\t", usecols=usecols)
    return df, evalue_source


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    start_time = time.time()
    print(f'{sys.argv[0]} started: {datetime.datetime.now()}', flush=True)

    print(f'Input file: {args.infile}', flush=True)
    df, evalue_source = read_required_columns(args.infile)

    df['qjointcov'] = pandas.to_numeric(df['qjointcov'], errors='coerce')
    if evalue_source == 'min_evalue':
        df['evalue_sort'] = pandas.to_numeric(df['min_evalue'], errors='coerce')
    else:
        df['evalue_sort'] = df['evalue'].map(parse_min_evalue)

    # Downstream filtering is by qjointcov threshold, so rows below this can be removed early.
    df = df.loc[df['qjointcov'] >= args.min_query_blast_coverage, ['sacc', 'qjointcov', 'evalue_sort']]
    if df.empty:
        gene_ids = pandas.Series(dtype=object)
        gene_ids.to_csv(args.outfile, sep='\t', header=False, index=False)
        print('Number of BLAST-hit gene IDs before filtering: 0', flush=True)
        print('Number of BLAST-hit gene IDs after filtering: 0', flush=True)
        print(f'{sys.argv[0]} completed: {datetime.datetime.now()}', flush=True)
        print(f'Elapsed time: {int(time.time() - start_time)} sec', flush=True)
        return

    # Choose the best hit per gene_id without globally sorting all rows:
    # 1) max qjointcov per sacc, 2) min evalue among ties.
    max_qjointcov = df.groupby('sacc', sort=False)['qjointcov'].transform('max')
    candidates = df.loc[df['qjointcov'].eq(max_qjointcov), ['sacc', 'qjointcov', 'evalue_sort']]
    evalue_rank = candidates['evalue_sort'].fillna(numpy.inf)
    min_evalue_rank = evalue_rank.groupby(candidates['sacc'], sort=False).transform('min')
    df2 = candidates.loc[evalue_rank.eq(min_evalue_rank), ['sacc', 'evalue_sort', 'qjointcov']]
    df2 = df2.drop_duplicates(subset=['sacc'], keep='first').reset_index(drop=True)

    df2 = df2.sort_values(
        by=['qjointcov', 'evalue_sort'],
        ascending=[False, True],
        na_position='last',
        kind='mergesort',
    ).reset_index(drop=True)
    df2.columns = ['gene_id', 'min_evalue', 'qjointcov']
    print(f'Number of BLAST-hit gene IDs before filtering: {df2.shape[0]}', flush=True)
    df3 = df2.head(args.max_num_gene_blast_hit_retrieval).reset_index(drop=True)
    gene_ids = df3.loc[:,'gene_id']
    gene_ids.to_csv(args.outfile, sep='\t', header=False, index=False)
    print(f'Number of BLAST-hit gene IDs after filtering: {gene_ids.shape[0]}', flush=True)
    print(f'{sys.argv[0]} completed: {datetime.datetime.now()}', flush=True)
    print(f'Elapsed time: {int(time.time() - start_time)} sec', flush=True)


if __name__ == '__main__':
    main()
