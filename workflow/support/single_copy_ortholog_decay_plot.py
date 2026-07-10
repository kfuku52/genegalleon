#!/usr/bin/env python3
# coding: utf-8

import argparse
import datetime
import os
import sys
import time

import numpy
import pandas

METRIC_LABELS = {
    "strict_single_copy": "Strictly single-copy orthogroups",
    "non_missing": "Non-missing orthogroups",
    "all_observed": "All observed orthogroups",
}


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthogroup-genecount", metavar="PATH", required=True, help="Orthogroups.GeneCount.tsv input.")
    parser.add_argument("--outdir", metavar="PATH", required=True, help="Output directory.")
    parser.add_argument("--replicates", metavar="INT", default=1000, type=int, help="Random species subsets per X value.")
    parser.add_argument("--species-counts", metavar="STR", default="auto", help="Comma-separated species counts, ranges like 1-10 or 1-10:2, or auto.")
    parser.add_argument("--seed", metavar="INT", default=1, type=int, help="Random seed.")
    parser.add_argument("--plot-basename", metavar="STR", default="single_copy_ortholog_decay_plot", help="Output plot basename.")
    parser.add_argument("--summary-name", metavar="STR", default="single_copy_ortholog_decay_summary.tsv", help="Output summary TSV name.")
    parser.add_argument("--formats", metavar="STR", default="pdf,svg", help="Comma-separated plot formats.")
    return parser


def get_pyplot():
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    matplotlib.rcParams["font.size"] = 8
    matplotlib.rcParams["font.family"] = "Helvetica"
    matplotlib.rcParams["svg.fonttype"] = "none"
    return plt


def parse_species_counts(species_counts, num_species):
    value = str(species_counts).strip()
    if value == "" or value.lower() == "auto":
        return list(range(1, num_species + 1))

    counts = []
    for raw_token in value.split(","):
        token = raw_token.strip()
        if token == "":
            continue
        if "-" in token:
            range_part, sep, step_part = token.partition(":")
            start_txt, dash, end_txt = range_part.partition("-")
            if dash == "":
                raise ValueError("Invalid species-count range: {}".format(token))
            start = int(start_txt)
            end = int(end_txt)
            step = int(step_part) if sep else 1
            if step <= 0:
                raise ValueError("Species-count range step must be positive: {}".format(token))
            if end < start:
                raise ValueError("Species-count range end is smaller than start: {}".format(token))
            counts.extend(range(start, end + 1, step))
        else:
            counts.append(int(token))

    counts = sorted(set(counts))
    if len(counts) == 0:
        raise ValueError("--species-counts did not contain any values.")
    invalid = [count for count in counts if count < 1 or count > num_species]
    if invalid:
        raise ValueError(
            "Species counts must be between 1 and {:,}; invalid values: {}".format(
                num_species, ", ".join(str(x) for x in invalid)
            )
        )
    return counts


def _is_species_count_column(col):
    if col in {"Orthogroup", "Total"}:
        return False
    if col.startswith("geneid_") or col.startswith("besthit_"):
        return False
    return True


def load_gene_count_table(path):
    df = pandas.read_csv(path, sep="\t", header=0, low_memory=False)
    if "Orthogroup" not in df.columns:
        raise ValueError("Input must contain an Orthogroup column: {}".format(path))
    species_cols = [col for col in df.columns if _is_species_count_column(col)]
    if len(species_cols) == 0:
        raise ValueError("No species count columns were detected in: {}".format(path))

    counts_df = df.loc[:, species_cols].apply(pandas.to_numeric, errors="raise")
    if counts_df.isna().any().any():
        raise ValueError("Species count columns contain missing values: {}".format(path))
    counts = counts_df.to_numpy(dtype=numpy.int64, copy=True)
    if (counts < 0).any():
        raise ValueError("Species count columns contain negative values: {}".format(path))
    return species_cols, counts


def calculate_decay(counts, species_counts, replicates, seed):
    if replicates <= 0:
        raise ValueError("--replicates must be a positive integer.")
    num_orthogroups, num_species = counts.shape
    max_species_count = max(species_counts)
    species_count_to_index = {count: idx for idx, count in enumerate(species_counts)}
    values = numpy.zeros((replicates, len(species_counts), 3), dtype=numpy.float64)
    rng = numpy.random.default_rng(seed)

    for replicate in range(replicates):
        species_order = rng.permutation(num_species)
        present_seen = numpy.zeros(num_orthogroups, dtype=bool)
        present_count = numpy.zeros(num_orthogroups, dtype=numpy.int32)
        single_count = numpy.zeros(num_orthogroups, dtype=numpy.int32)

        for position in range(1, max_species_count + 1):
            selected_counts = counts[:, species_order[position - 1]]
            present = selected_counts >= 1
            present_seen |= present
            present_count += present
            single_count += selected_counts == 1

            summary_index = species_count_to_index.get(position)
            if summary_index is None:
                continue
            strict_single_copy = int((single_count == position).sum())
            non_missing = int((present_count == position).sum())
            all_observed = int(present_seen.sum())
            values[replicate, summary_index, :] = [strict_single_copy, non_missing, all_observed]

    return values


def summarize_decay(species_counts, values):
    rows = []
    metrics = [
        ("strict_single_copy", values[:, :, 0]),
        ("non_missing", values[:, :, 1]),
        ("all_observed", values[:, :, 2]),
    ]
    for metric, metric_values in metrics:
        means = metric_values.mean(axis=0)
        sds = metric_values.std(axis=0, ddof=1) if metric_values.shape[0] > 1 else numpy.zeros(len(species_counts))
        for species_count, mean, sd in zip(species_counts, means, sds):
            rows.append(
                {
                    "species_count": species_count,
                    "metric": metric,
                    "metric_label": METRIC_LABELS[metric],
                    "mean": mean,
                    "sd": sd,
                    "replicates": metric_values.shape[0],
                }
            )
    return pandas.DataFrame(rows)


def _metric_array(summary, metric, value_col):
    metric_rows = summary.loc[summary["metric"] == metric, :].sort_values("species_count")
    return metric_rows[value_col].to_numpy(dtype=float)


def plot_decay(summary, outdir, basename, formats):
    plt = get_pyplot()
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    summary = summary.copy()
    x = (
        summary.loc[summary["metric"] == "strict_single_copy", :]
        .sort_values("species_count")["species_count"]
        .to_numpy(dtype=float)
    )
    strict_mean = _metric_array(summary, "strict_single_copy", "mean")
    strict_sd = _metric_array(summary, "strict_single_copy", "sd")
    non_missing_mean = _metric_array(summary, "non_missing", "mean")
    non_missing_sd = _metric_array(summary, "non_missing", "sd")
    all_observed_mean = _metric_array(summary, "all_observed", "mean")
    all_observed_sd = _metric_array(summary, "all_observed", "sd")

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5.6, 3.6), sharey=False, sharex=False)

    line_specs = [
        ("strict_single_copy", strict_mean, strict_sd, "#0072B2"),
        ("non_missing", non_missing_mean, non_missing_sd, "#009E73"),
        ("all_observed", all_observed_mean, all_observed_sd, "#D55E00"),
    ]
    for metric, mean, sd, color in line_specs:
        lower = numpy.maximum(mean - sd, 0)
        upper = mean + sd
        finite_band = numpy.isfinite(lower) & numpy.isfinite(upper)
        if finite_band.any():
            ax.fill_between(x, lower, upper, where=finite_band, color=color, alpha=0.08, linewidth=0)
        ax.plot(x, mean, color=color, linewidth=1.2, label=METRIC_LABELS[metric])

    if x.size == 1:
        ax.set_xlim(float(x[0]) - 0.5, float(x[0]) + 0.5)
    else:
        ax.set_xlim(float(x.min()), float(x.max()))
    y_upper = max(all_observed_mean + all_observed_sd) * 1.15
    if y_upper <= 0:
        y_upper = 1.0
    ax.set_ylim(0, y_upper)
    ax.set_xlabel("Number of selected species")
    ax.set_ylabel("Number of orthogroups")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: "{:,.0f}".format(y)))
    ax.grid(True, which="major", color="#cccccc", linewidth=0.4, alpha=0.7)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc="best", frameon=False)
    fig.tight_layout()

    os.makedirs(outdir, exist_ok=True)
    output_paths = []
    for ext in formats:
        outpath = os.path.join(outdir, "{}.{}".format(basename, ext))
        fig.savefig(outpath, format=ext)
        output_paths.append(outpath)
    plt.close(fig)
    return output_paths


def parse_formats(formats):
    out = []
    for raw_token in str(formats).split(","):
        token = raw_token.strip().lower()
        if token == "":
            continue
        out.append(token)
    if len(out) == 0:
        raise ValueError("--formats did not contain any plot formats.")
    return out


def run(args):
    if args.replicates <= 0:
        raise ValueError("--replicates must be a positive integer.")
    start = time.time()
    print("Starting {} at {}".format(sys.argv[0], datetime.datetime.now()))
    print("Loading: {}".format(args.orthogroup_genecount))

    species_cols, counts = load_gene_count_table(args.orthogroup_genecount)
    species_counts = parse_species_counts(args.species_counts, len(species_cols))
    formats = parse_formats(args.formats)
    print("Number of orthogroups: {:,}".format(counts.shape[0]))
    print("Number of species: {:,}".format(len(species_cols)))
    print("Species-count values: {}".format(",".join(str(x) for x in species_counts)))
    print("Replicates per species-count value: {:,}".format(args.replicates))
    print("Random seed: {}".format(args.seed))

    values = calculate_decay(counts, species_counts, args.replicates, args.seed)
    summary = summarize_decay(species_counts, values)
    os.makedirs(args.outdir, exist_ok=True)
    summary_path = os.path.join(args.outdir, args.summary_name)
    print("Writing: {}".format(summary_path))
    summary.to_csv(summary_path, sep="\t", index=False)
    for outpath in plot_decay(summary, args.outdir, args.plot_basename, formats):
        print("Writing: {}".format(outpath))

    print("Ending {} at {}. Elapsed time: {:,} sec".format(sys.argv[0], datetime.datetime.now(), int(time.time() - start)))


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
