#!/usr/bin/env python3
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
    parser.add_argument(
        "--file_orthogroup_table",
        metavar="PATH",
        default="",
        type=str,
        help="Path used by --file_orthogroup_table.",
        required=True,
    )
    parser.add_argument(
        "--mode",
        metavar="STR",
        default="hog2og",
        choices=["hog2og", "sp2og"],
        type=str,
        help="Value used by --mode.",
        required=True,
    )
    parser.add_argument(
        "--dir_out", metavar="PATH", default="", type=str, help="Path used by --dir_out.", required=True
    )
    return parser


def build_gene_count_table(df, sp_cols):
    columns = {"Orthogroup": df["Orthogroup"]}
    total = numpy.zeros(len(df), dtype=int)

    for sp_col in sp_cols:
        # Materialize one object column to avoid per-scalar Arrow lookups,
        # without creating fixed-width Unicode arrays for giant cells.
        series = df.loc[:, sp_col].fillna("").astype(str)
        counts = numpy.fromiter(
            (value.count(",") + bool(value) for value in series.to_numpy(dtype=object)),
            dtype=int,
            count=len(series),
        )
        columns[sp_col] = counts
        total += counts

    columns["Total"] = total
    return pandas.DataFrame(columns, index=df.index)


def run(args):
    start = time.time()
    print("Starting {} at {}".format(sys.argv[0], datetime.datetime.now()))

    df = pandas.read_csv(
        args.file_orthogroup_table,
        sep="\t",
        header=0,
        low_memory=False,
        dtype=str,
        keep_default_na=False,
    )
    if args.mode == "hog2og":
        df = df.drop(["OG", "Gene Tree Parent Clade"], axis=1)
        df.columns = df.columns.str.replace("^HOG$", "Orthogroup", regex=True)
        pref = os.path.splitext(os.path.basename(args.file_orthogroup_table))[0] + "."
        has_pref = df["Orthogroup"].str.startswith(pref, na=False)
        df.loc[has_pref, "Orthogroup"] = df.loc[has_pref, "Orthogroup"].str[len(pref) :]
    elif args.mode == "sp2og":
        df.columns = df.columns.str.replace("^group_id$", "Orthogroup", regex=True)
        df.columns = df.columns.str.replace(r"\..*", "", regex=True)
        df["Orthogroup"] = "SP" + df["Orthogroup"].astype(str).str.zfill(7)
        species_cols = [c for c in df.columns if c != "Orthogroup"]
        df.loc[:, species_cols] = df.loc[:, species_cols].apply(
            lambda col: col.astype(str).str.replace("*", "", regex=False).str.replace(",", ", ", regex=False)
        )

    if df["Orthogroup"].eq("").any():
        raise ValueError("Orthogroup table contains an empty family identity")
    duplicates = df.loc[df["Orthogroup"].duplicated(), "Orthogroup"].unique()
    if len(duplicates):
        raise ValueError("Orthogroup table contains duplicate family identities: " + ", ".join(duplicates[:5]))

    if os.path.exists(args.dir_out):
        print("--out_dir exists already: {}".format(args.dir_out))
    else:
        print("Creating --out_dir: {}".format(args.dir_out))
        os.makedirs(args.dir_out)

    outpath_og = os.path.join(args.dir_out, "Orthogroups.tsv")
    print("Writing: {}".format(outpath_og))
    df.to_csv(outpath_og, sep="\t", index=False)

    sp_cols = df.columns[df.columns != "Orthogroup"]
    gc = build_gene_count_table(df, sp_cols)
    print("Number of Orthogroups: {:,}".format(gc.shape[0]), flush=True)
    print("Largest Orthogroup size: {:,}".format(gc["Total"].max()), flush=True)
    print("Smallest Orthogroup size: {:,}".format(gc["Total"].min()), flush=True)

    outpath_gc = os.path.join(args.dir_out, "Orthogroups.GeneCount.tsv")
    print("Writing: {}".format(outpath_gc))
    gc.to_csv(outpath_gc, sep="\t", index=False)

    txt_readme = "The files in this directory were created by {}. ".format(sys.argv[0])
    txt_readme += 'The OrthoFinder outputs in the obsoleted "Orthogroups" directory were '
    txt_readme += "mimicked by parsing {}\n".format(args.file_orthogroup_table)

    outpath_readme = os.path.join(args.dir_out, "README.txt")
    with open(outpath_readme, "w") as f:
        f.write(txt_readme)

    print(
        "Ending {} at {}. Elapsed time: {:,} sec".format(sys.argv[0], datetime.datetime.now(), int(time.time() - start))
    )


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
