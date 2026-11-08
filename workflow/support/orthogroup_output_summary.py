#!/usr/bin/env python3
# coding: utf-8

import argparse
import os
import re
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
from gene_family_output_store import SHARED_OUTPUT_SUBDIRS, GeneFamilyOutputStore

ORTHOGROUP_ID_RE = re.compile(r"^(OG\d+|HOG\d+|SP\d+)(?=$|[._-])")




def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir_og", metavar="PATH", type=str, required=True, help="Path used by --dir_og.")
    parser.add_argument("--genecount", metavar="PATH", type=str, required=True, help="Path used by --genecount.")
    parser.add_argument("--out", metavar="PATH", type=str, required=True, help="Path used by --out.")
    parser.add_argument(
        "--updated-genecount-out",
        metavar="PATH",
        type=str,
        help="Write the alignment-statistics-augmented gene-count table here; the default is beside --out.",
    )
    parser.add_argument("--ncpu", metavar="INT", default=1, type=int, help="Number of worker threads.")
    return parser


def _extract_orthogroup_id(file_name):
    match = ORTHOGROUP_ID_RE.match(os.path.basename(file_name))
    if match is None:
        return None
    return match.group(1)




def get_alignment_stats(df, dir_alignment_stats, extension, ncpu, store=None, logical_subdir=None):
    df, found = append_alignment_stats(df, dir_alignment_stats, extension, ncpu, _extract_orthogroup_id,
                                store, logical_subdir, columns_before_check=True)
    if not found:
        return df
    idx_total = numpy.where(df.columns == "Total")[0][0]
    original_cols = df.columns[:idx_total + 1].tolist()
    sorted_alignment_stats_cols = df.columns[idx_total + 1:].sort_values().tolist()
    return df.loc[:, original_cols + sorted_alignment_stats_cols]


def run(args):
    print("args:", vars(args))

    start = time.time()
    updated_genecount = getattr(args, "updated_genecount_out", None)
    if updated_genecount is None:
        out_path = Path(args.out)
        updated_genecount = str(out_path.with_name(f"{out_path.stem}.genecount.alignment_stats.tsv"))
    df_original = pandas.read_csv(args.genecount, sep="\t", index_col=0, header=0)
    df_original.index = df_original.index.astype(str)
    if os.path.exists(updated_genecount):
        print("Updated --genecount file was detected. Reading and aligning to the original index.", flush=True)
        df_updated = pandas.read_csv(updated_genecount, sep="\t", index_col=0, header=0)
        df_updated.index = df_updated.index.astype(str)
        df = df_original.copy()
        extra_cols = [col for col in df_updated.columns if col not in df.columns]
        if extra_cols:
            df = pandas.concat([df, df_updated.reindex(df.index)[extra_cols]], axis=1)
    else:
        print("Updated --genecount file was not detected. Reading the original.", flush=True)
        df = df_original

    store = GeneFamilyOutputStore(args.dir_og)
    dir_alignment_stats = os.path.join(args.dir_og, "alignment_stats_original")
    df = get_alignment_stats(
        df,
        dir_alignment_stats,
        "original",
        args.ncpu,
        store=store,
        logical_subdir="alignment_stats_original",
    )
    dir_alignment_stats = os.path.join(args.dir_og, "alignment_stats_cleaned")
    df = get_alignment_stats(
        df,
        dir_alignment_stats,
        "clean",
        args.ncpu,
        store=store,
        logical_subdir="alignment_stats_cleaned",
    )
    df.to_csv(updated_genecount, index=True, sep="\t")

    df.loc[:, "GG_ARRAY_TASK_ID"] = numpy.arange(df.shape[0]) + 1

    subdirs = [subdir for subdir in store.logical_subdirs() if subdir not in SHARED_OUTPUT_SUBDIRS]
    df = pandas.concat([df, pandas.DataFrame(data=0, index=df.index, columns=subdirs, dtype=int)], axis=1)
    valid_og_ids = set(df.index.astype(str))
    for subdir in subdirs:
        files = store.file_names(subdir)
        og_ids = sorted(
            {
                og_id
                for og_id in (_extract_orthogroup_id(f) for f in files)
                if og_id is not None and og_id in valid_og_ids
            }
        )
        if og_ids:
            df.loc[og_ids, subdir] = 1
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
