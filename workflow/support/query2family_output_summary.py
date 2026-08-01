#!/usr/bin/env python3
# coding: utf-8

import argparse
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy
import pandas

SUPPORT_DIR = Path(__file__).resolve().parent
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))

from gene_family_output_store import GeneFamilyOutputStore


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir_query2family', metavar='PATH', type=str, required=True, help='')
    parser.add_argument('--dir_query_gene', metavar='PATH', type=str, required=True, help='')
    parser.add_argument('--out', metavar='PATH', type=str, required=True, help='')
    parser.add_argument('--ncpu', metavar='INT', default=1, type=int, help='Number of worker threads.')
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


def _visible_entries(path):
    return [entry for entry in os.listdir(path) if not entry.startswith('.')]


def _query_ids_from_input_dir(path):
    query_dir = Path(path)
    if not query_dir.is_dir():
        raise FileNotFoundError(f'Input query_gene directory was not found: {query_dir}')
    return sorted(entry.name for entry in query_dir.iterdir() if entry.is_file() and not entry.name.startswith('.'))


def _query_id_matchers(query_ids):
    return sorted([str(query_id) for query_id in query_ids], key=lambda x: (-len(x), x))


def _extract_query_id(file_name, query_id_matchers):
    basename = os.path.basename(file_name)
    for query_id in query_id_matchers:
        if basename == query_id:
            return query_id
        if basename.startswith(query_id + '_') or basename.startswith(query_id + '.'):
            return query_id
    return None


def _read_amas_file(file_path, query_id, amas_cols):
    tmp = pandas.read_csv(file_path, sep='\t', header=0, usecols=amas_cols, nrows=1, low_memory=False)
    return query_id, tmp.iloc[0].to_list()


def _read_amas_store_file(store, subdir, file_name, query_id, amas_cols):
    with store.open_binary(subdir, file_name) as handle:
        tmp = pandas.read_csv(handle, sep='\t', header=0, usecols=amas_cols, nrows=1, low_memory=False)
    return query_id, tmp.iloc[0].to_list()


def get_amas_stats(
    df,
    dir_amas,
    extension,
    query_id_matchers,
    ncpu,
    store=None,
    logical_subdir=None,
):
    if store is None and not os.path.isdir(dir_amas):
        print(f'{extension}: {dir_amas} was not found. Skipping.', flush=True)
        return df
    if store is not None and logical_subdir not in store.logical_subdirs():
        print(f'{extension}: logical subdirectory {logical_subdir} was not found. Skipping.', flush=True)
        return df

    amas_cols = _amas_columns()
    amas_new_cols = [f'{col}_{extension}' for col in amas_cols]
    for ncol in amas_new_cols:
        if ncol not in df.columns:
            df.loc[:, ncol] = numpy.nan

    is_prefilled = ~df[f'No_of_taxa_{extension}'].isna()
    prefilled_query_ids = set(df.index[is_prefilled])
    files = (
        sorted(_visible_entries(dir_amas))
        if store is None
        else store.file_names(logical_subdir)
    )
    queued = []
    seen_query_ids = set()
    valid_query_ids = set(df.index.astype(str))
    for file in files:
        file_path = os.path.join(dir_amas, file)
        if store is None and not os.path.isfile(file_path):
            continue
        query_id = _extract_query_id(file, query_id_matchers)
        if query_id is None:
            continue
        if query_id in seen_query_ids:
            continue
        if query_id in prefilled_query_ids:
            continue
        if query_id not in valid_query_ids:
            continue
        seen_query_ids.add(query_id)
        queued.append((file, query_id))

    counter = 0
    result_rows = []
    if ncpu > 1 and len(queued) > 1:
        max_workers = min(ncpu, len(queued))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(
                    _read_amas_file,
                    os.path.join(dir_amas, file),
                    query_id,
                    amas_cols,
                )
                if store is None
                else executor.submit(
                    _read_amas_store_file,
                    store,
                    logical_subdir,
                    file,
                    query_id,
                    amas_cols,
                )
                for file, query_id in queued
            ]
            for future in as_completed(futures):
                query_id, values = future.result()
                result_rows.append((query_id, values))
                counter += 1
    else:
        for file, query_id in queued:
            if store is None:
                query_id, values = _read_amas_file(os.path.join(dir_amas, file), query_id, amas_cols)
            else:
                query_id, values = _read_amas_store_file(
                    store,
                    logical_subdir,
                    file,
                    query_id,
                    amas_cols,
                )
            result_rows.append((query_id, values))
            counter += 1

    if result_rows:
        result_df = pandas.DataFrame.from_records(
            [values for _, values in result_rows],
            index=[query_id for query_id, _ in result_rows],
            columns=amas_new_cols,
        )
        df.loc[result_df.index, amas_new_cols] = result_df

    print(f'{extension}: {counter} AMAS results were appended.', flush=True)
    return df


def run(args):
    print('args:', vars(args))

    start = time.time()
    query_ids = _query_ids_from_input_dir(args.dir_query_gene)
    if len(query_ids) == 0:
        raise ValueError(f'Input query_gene directory is empty: {args.dir_query_gene}')

    query_id_matchers = _query_id_matchers(query_ids)
    df = pandas.DataFrame(index=query_ids)
    df.index.name = 'query'
    store = GeneFamilyOutputStore(args.dir_query2family)

    df = get_amas_stats(
        df,
        os.path.join(args.dir_query2family, 'amas_original'),
        'original',
        query_id_matchers,
        args.ncpu,
        store=store,
        logical_subdir='amas_original',
    )
    df = get_amas_stats(
        df,
        os.path.join(args.dir_query2family, 'amas_cleaned'),
        'clean',
        query_id_matchers,
        args.ncpu,
        store=store,
        logical_subdir='amas_cleaned',
    )

    df.insert(0, 'GG_ARRAY_TASK_ID', numpy.arange(df.shape[0]) + 1)

    subdirs = store.logical_subdirs()
    df = pandas.concat([df, pandas.DataFrame(data=0, index=df.index, columns=subdirs, dtype=int)], axis=1)
    valid_query_ids = set(df.index.astype(str))
    for subdir in subdirs:
        files = store.file_names(subdir)
        query_ids_in_files = sorted(
            {
                query_id
                for query_id in (_extract_query_id(f, query_id_matchers) for f in files)
                if query_id is not None and query_id in valid_query_ids
            }
        )
        if query_ids_in_files:
            df.loc[query_ids_in_files, subdir] = 1
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
