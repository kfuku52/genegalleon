#!/usr/bin/env python3
# coding: utf-8

import argparse
import os
import sys
import time
from pathlib import Path

import numpy
import pandas

SUPPORT_DIR = Path(__file__).resolve().parent
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))

from alignment_statistics import _alignment_stats_columns as _alignment_stats_columns
from alignment_statistics import _read_alignment_stats_file as _read_alignment_stats_file
from alignment_statistics import _read_alignment_stats_store_file as _read_alignment_stats_store_file
from alignment_statistics import _visible_entries as _visible_entries
from alignment_statistics import append_alignment_stats
from gene_family_output_store import SHARED_OUTPUT_SUBDIRS, GeneFamilyOutputStore, query_id_extractor


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dir_query2family", metavar="PATH", type=str, required=True, help="Path used by --dir_query2family."
    )
    parser.add_argument(
        "--dir_query_gene", metavar="PATH", type=str, required=True, help="Path used by --dir_query_gene."
    )
    parser.add_argument("--out", metavar="PATH", type=str, required=True, help="Path used by --out.")
    parser.add_argument("--ncpu", metavar="INT", default=1, type=int, help="Number of worker threads.")
    return parser




def _query_ids_from_input_dir(path):
    query_dir = Path(path)
    if not query_dir.is_dir():
        raise FileNotFoundError(f"Input query_gene directory was not found: {query_dir}")
    return sorted(entry.name for entry in query_dir.iterdir() if entry.is_file() and not entry.name.startswith("."))


def _query_id_matchers(query_ids):
    return sorted([str(query_id) for query_id in query_ids], key=lambda x: (-len(x), x))


def _extract_query_id(file_name, query_id_matchers):
    basename = os.path.basename(file_name)
    for query_id in query_id_matchers:
        if basename == query_id:
            return query_id
        if basename.startswith(query_id + "_") or basename.startswith(query_id + "."):
            return query_id
    return None


def get_alignment_stats(df, dir_alignment_stats, extension, query_id_matchers, ncpu, store=None, logical_subdir=None):
    df, _found = append_alignment_stats(df, dir_alignment_stats, extension, ncpu, query_id_extractor(query_id_matchers),
                                  store, logical_subdir)
    return df


def run(args):
    print("args:", vars(args))

    start = time.time()
    query_ids = _query_ids_from_input_dir(args.dir_query_gene)
    if len(query_ids) == 0:
        raise ValueError(f"Input query_gene directory is empty: {args.dir_query_gene}")

    query_id_matchers = _query_id_matchers(query_ids)
    extract_query_id = query_id_extractor(query_id_matchers)
    df = pandas.DataFrame(index=query_ids)
    df.index.name = "query"
    store = GeneFamilyOutputStore(args.dir_query2family)

    df = get_alignment_stats(
        df,
        os.path.join(args.dir_query2family, "alignment_stats_original"),
        "original",
        query_id_matchers,
        args.ncpu,
        store=store,
        logical_subdir="alignment_stats_original",
    )
    df = get_alignment_stats(
        df,
        os.path.join(args.dir_query2family, "alignment_stats_cleaned"),
        "clean",
        query_id_matchers,
        args.ncpu,
        store=store,
        logical_subdir="alignment_stats_cleaned",
    )

    df.insert(0, "GG_ARRAY_TASK_ID", numpy.arange(df.shape[0]) + 1)

    subdirs = [subdir for subdir in store.logical_subdirs() if subdir not in SHARED_OUTPUT_SUBDIRS]
    df = pandas.concat([df, pandas.DataFrame(data=0, index=df.index, columns=subdirs, dtype=int)], axis=1)
    valid_query_ids = set(df.index.astype(str))
    for subdir in subdirs:
        files = store.file_names(subdir)
        query_ids_in_files = sorted(
            {
                query_id
                for query_id in (extract_query_id(f) for f in files)
                if query_id is not None and query_id in valid_query_ids
            }
        )
        if query_ids_in_files:
            df.loc[query_ids_in_files, subdir] = 1
        num_missing = (df.loc[:, subdir] == 0).sum()
        txt = "Subdirectory {}: {:,} / {:,} files are missing."
        print(txt.format(subdir, num_missing, df.shape[0]))

    print("Writing output file:", args.out, flush=True)
    df.to_csv(args.out, index=True, sep="\t")
    print("Done. Elapsed time: {:,} sec".format(int(time.time() - start)))


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    args.ncpu = max(1, int(args.ncpu))
    run(args)


if __name__ == "__main__":
    main()
