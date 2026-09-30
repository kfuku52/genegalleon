#!/usr/bin/env python3
# coding: utf-8

import argparse
import datetime
import os
import re
import sys
import time
from pathlib import Path

import pandas

pandas.set_option("display.max_columns", None)


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--orthofinder_og_genecount",
        metavar="PATH",
        default="",
        type=str,
        help="Flat OG gene-count TSV or native Orthogroups.txt membership file.",
        required=True,
    )
    parser.add_argument(
        "--orthofinder_hog_genecount",
        metavar="PATH",
        default="",
        type=str,
        help="Path used by --orthofinder_hog_genecount.",
        required=True,
    )
    return parser


def read_gene_counts(path):
    """Read the totals used by the comparison without changing OG membership."""
    if Path(path).name != "Orthogroups.txt":
        return pandas.read_csv(path, sep="\t", header=0, low_memory=False)

    rows = []
    seen_groups = set()
    seen_genes = set()
    with open(path) as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            match = re.fullmatch(r"(OG[0-9]+):\s*(\S+(?:\s+\S+)*)\s*", line.strip())
            if match is None:
                raise ValueError(f"Invalid Orthogroups.txt membership at line {number}")
            group, members = match.groups()
            genes = members.split()
            if group in seen_groups:
                raise ValueError(f"Duplicate orthogroup {group} at line {number}")
            if len(set(genes)) != len(genes) or seen_genes.intersection(genes):
                raise ValueError(f"Duplicate gene membership at line {number}")
            seen_groups.add(group)
            seen_genes.update(genes)
            rows.append((group, len(genes)))
    if not rows:
        raise ValueError("Orthogroups.txt has no orthogroups")
    return pandas.DataFrame(rows, columns=["Orthogroup", "Total"])


def get_pyplot():
    import matplotlib
    import matplotlib.pyplot as plt

    matplotlib.rcParams["font.size"] = 8
    matplotlib.rcParams["font.family"] = "Helvetica"
    matplotlib.rcParams["svg.fonttype"] = "none"  # none, path, or svgfont
    return plt


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    start = time.time()
    cwd = os.getcwd()
    print("Working at: {}".format(cwd))
    print("Starting {} at {}".format(sys.argv[0], datetime.datetime.now()))

    dfs = {}
    dfs["OrthoFinder Orthogroup"] = read_gene_counts(args.orthofinder_og_genecount)
    dfs["OrthoFinder Hierarchical orthogroup"] = read_gene_counts(args.orthofinder_hog_genecount)
    for key in dfs.keys():
        dfs[key] = dfs[key].sort_values(by="Total")
        dfs[key]["cumulative_num_gene"] = dfs[key]["Total"].cumsum()

    colors = {
        "OrthoFinder Orthogroup": "#991574",
        "OrthoFinder Hierarchical orthogroup": "#749915",
    }

    plt = get_pyplot()
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(7.2, 6.4), sharey=False, sharex=False)
    axes = axes.flat

    max_ngenes = max([dfs[key]["Total"].max() for key in dfs.keys()])
    xmax = min(max_ngenes, 1000)
    bins = range(0, xmax, 1)

    ax = axes[0]
    for key in dfs.keys():
        ax.hist(
            x=dfs[key]["Total"],
            bins=bins,
            color=colors[key],
            histtype="step",
            cumulative=True,
            label=key,
            linewidth=1.0,
            alpha=0.5,
        )
    ax.set_xlim(0, xmax)
    ax.set_xlabel("Number of genes per orthogroup")
    ax.set_ylabel("Cumulative number of orthogroups")
    ax.legend(loc="lower right")

    ax = axes[1]
    for key in dfs.keys():
        ax.hist(x=dfs[key]["Total"], bins=bins, color=colors[key], histtype="step", label=key, linewidth=1.0, alpha=0.5)
    ax.set_xlim(0, xmax)
    ax.set_xlabel("Number of genes per orthogroup")
    ax.set_ylabel("Number of orthogroups")
    ax.legend(loc="upper right")

    ax = axes[2]
    for key in dfs.keys():
        ax.plot(
            dfs[key]["Total"], dfs[key]["cumulative_num_gene"], color=colors[key], label=key, linewidth=1.0, alpha=0.5
        )
    ax.set_xlim(0, xmax)
    ax.set_xlabel("Number of genes per orthogroup")
    ax.set_ylabel("Cumulative number of genes")
    ax.legend(loc="lower right")

    ax = axes[3]
    for key in dfs.keys():
        tmp = dfs[key].copy()
        tmp["Total2"] = tmp["Total"].copy()
        tmp = tmp.groupby("Total2")["Total"].sum().reset_index()
        tmp2 = pandas.DataFrame({"Total2": range(0, xmax + 1)})
        tmp = pandas.merge(tmp2, tmp, on="Total2", how="left")
        tmp["Total"] = tmp["Total"].fillna(0)
        ax.plot(tmp["Total2"], tmp["Total"], color=colors[key], label=key, linewidth=1.0, alpha=0.5)
    ax.set_xlim(0, xmax)
    ax.set_xlabel("Number of genes per orthogroup")
    ax.set_ylabel("Number of genes")
    ax.legend(loc="upper right")

    outbase = "orthogroup_histogram"
    fig.tight_layout(pad=0.25, w_pad=1.5, h_pad=0.5)
    for ext in ["svg", "pdf"]:
        outpath = os.path.join(outbase + "." + ext)
        fig.savefig(outpath, format=ext)

    print(
        "Ending {} at {}. Elapsed time: {:,} sec".format(sys.argv[0], datetime.datetime.now(), int(time.time() - start))
    )


if __name__ == "__main__":
    main()
