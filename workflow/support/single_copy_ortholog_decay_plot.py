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
    "selected_observed": "Selected orthogroups",
    "all_observed": "All orthogroups",
}


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--orthogroup-genecount", metavar="PATH", required=True, help="Orthogroups.GeneCount.tsv input."
    )
    parser.add_argument(
        "--selected-orthogroup-genecount",
        metavar="PATH",
        required=True,
        help="Selected Orthogroups.GeneCount.selected.tsv input.",
    )
    parser.add_argument("--outdir", metavar="PATH", required=True, help="Output directory.")
    parser.add_argument(
        "--replicates", metavar="INT", default=1000, type=int, help="Random species subsets per X value."
    )
    parser.add_argument(
        "--species-counts",
        metavar="STR",
        default="auto",
        help="Comma-separated species counts, ranges like 1-10 or 1-10:2, or auto.",
    )
    parser.add_argument("--seed", metavar="INT", default=1, type=int, help="Random seed.")
    parser.add_argument(
        "--plot-basename", metavar="STR", default="single_copy_ortholog_decay_plot", help="Output plot basename."
    )
    parser.add_argument(
        "--summary-name",
        metavar="STR",
        default="single_copy_ortholog_decay_summary.tsv",
        help="Output summary TSV name.",
    )
    parser.add_argument("--formats", metavar="STR", default="pdf,svg", help="Comma-separated plot formats.")
    return parser


def get_pyplot():
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    matplotlib.rcParams.update(
        {
            "font.family": "Helvetica",
            "font.size": 8,
            "axes.titlesize": 8,
            "axes.labelsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "figure.titlesize": 8,
            "svg.fonttype": "none",
        }
    )
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


def calculate_decay(counts, species_counts, replicates, seed, selected_counts=None):
    if replicates <= 0:
        raise ValueError("--replicates must be a positive integer.")
    num_orthogroups, num_species = counts.shape
    if selected_counts is not None and selected_counts.shape[1] != num_species:
        raise ValueError("All and selected gene-count matrices must have the same number of species columns.")
    max_species_count = max(species_counts)
    species_count_to_index = {count: idx for idx, count in enumerate(species_counts)}
    num_metrics = 4 if selected_counts is not None else 3
    values = numpy.zeros((replicates, len(species_counts), num_metrics), dtype=numpy.float64)
    rng = numpy.random.default_rng(seed)

    for replicate in range(replicates):
        species_order = rng.permutation(num_species)
        present_seen = numpy.zeros(num_orthogroups, dtype=bool)
        present_count = numpy.zeros(num_orthogroups, dtype=numpy.int32)
        single_count = numpy.zeros(num_orthogroups, dtype=numpy.int32)
        if selected_counts is not None:
            selected_present_seen = numpy.zeros(selected_counts.shape[0], dtype=bool)

        for position in range(1, max_species_count + 1):
            current_counts = counts[:, species_order[position - 1]]
            present = current_counts >= 1
            present_seen |= present
            present_count += present
            single_count += current_counts == 1
            if selected_counts is not None:
                selected_present_seen |= selected_counts[:, species_order[position - 1]] >= 1

            summary_index = species_count_to_index.get(position)
            if summary_index is None:
                continue
            strict_single_copy = int((single_count == position).sum())
            non_missing = int((present_count == position).sum())
            all_observed = int(present_seen.sum())
            if selected_counts is None:
                values[replicate, summary_index, :] = [strict_single_copy, non_missing, all_observed]
            else:
                selected_observed = int(selected_present_seen.sum())
                values[replicate, summary_index, :] = [
                    strict_single_copy,
                    non_missing,
                    selected_observed,
                    all_observed,
                ]

    return values


def summarize_decay(species_counts, values):
    rows = []
    if values.shape[2] == 3:
        metrics = [
            ("strict_single_copy", values[:, :, 0]),
            ("non_missing", values[:, :, 1]),
            ("all_observed", values[:, :, 2]),
        ]
    elif values.shape[2] == 4:
        metrics = [
            ("strict_single_copy", values[:, :, 0]),
            ("non_missing", values[:, :, 1]),
            ("selected_observed", values[:, :, 2]),
            ("all_observed", values[:, :, 3]),
        ]
    else:
        raise ValueError("Decay values must contain three or four metrics.")
    for metric, metric_values in metrics:
        means = metric_values.mean(axis=0)
        sds = metric_values.std(axis=0, ddof=1) if metric_values.shape[0] > 1 else numpy.zeros(len(species_counts))
        for species_count, mean, sd in zip(species_counts, means, sds, strict=True):
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


def truncate_at_mean_floor(x, mean, sd, mean_floor):
    """Stop a series at its first mean-floor crossing, interpolating in log space."""
    at_or_below = numpy.flatnonzero(mean <= mean_floor)
    if at_or_below.size == 0:
        return x, mean, sd

    stop = int(at_or_below[0])
    if mean[stop] == mean_floor:
        return x[: stop + 1], mean[: stop + 1], sd[: stop + 1]
    if stop == 0:
        return x[:0], mean[:0], sd[:0]

    previous_mean = float(mean[stop - 1])
    current_mean = float(mean[stop])
    if current_mean > 0:
        fraction = (numpy.log(mean_floor) - numpy.log(previous_mean)) / (
            numpy.log(current_mean) - numpy.log(previous_mean)
        )
    else:
        fraction = (mean_floor - previous_mean) / (current_mean - previous_mean)
    crossing_x = float(x[stop - 1]) + fraction * float(x[stop] - x[stop - 1])
    crossing_sd = float(sd[stop - 1]) + fraction * float(sd[stop] - sd[stop - 1])
    return (
        numpy.append(x[:stop], crossing_x),
        numpy.append(mean[:stop], mean_floor),
        numpy.append(sd[:stop], crossing_sd),
    )


def plot_decay(summary, outdir, basename, formats, selected_total):
    plt = get_pyplot()
    from matplotlib.ticker import FuncFormatter, LogLocator, MaxNLocator, NullFormatter

    summary = summary.copy()
    x = (
        summary.loc[summary["metric"] == "strict_single_copy", :]
        .sort_values("species_count")["species_count"]
        .to_numpy(dtype=float)
    )
    line_specs = [
        ("all_observed", "#D55E00", "-"),
        ("selected_observed", "#CC79A7", "-."),
        ("non_missing", "#009E73", "--"),
        ("strict_single_copy", "#0072B2", ":"),
    ]
    log_floor = 1.0

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3.6, 4.4), sharey=False, sharex=False)
    for metric, color, linestyle in line_specs:
        mean = _metric_array(summary, metric, "mean")
        sd = _metric_array(summary, metric, "sd")
        x_plot, mean_plot, sd_plot = truncate_at_mean_floor(x, mean, sd, log_floor)
        lower = numpy.maximum(mean_plot - sd_plot, log_floor)
        upper = mean_plot + sd_plot
        finite_band = numpy.isfinite(lower) & numpy.isfinite(upper) & (upper >= log_floor)
        if finite_band.any():
            ax.fill_between(
                x_plot,
                lower,
                upper,
                where=finite_band,
                color=color,
                alpha=0.11,
                linewidth=0,
            )
        ax.plot(
            x_plot,
            mean_plot,
            color=color,
            linestyle=linestyle,
            linewidth=1.6,
            label=METRIC_LABELS[metric],
        )

    if x.size == 1:
        ax.set_xlim(float(x[0]) - 0.5, float(x[0]) + 0.5)
    else:
        ax.set_xlim(float(x.min()), float(x.max()))
    all_upper = _metric_array(summary, "all_observed", "mean") + _metric_array(
        summary, "all_observed", "sd"
    )
    y_upper = max(log_floor * 1.35, float(all_upper.max()) * 1.35)
    ax.set_yscale("log", base=10)
    ax.set_ylim(log_floor, y_upper)
    ax.set_xlabel("Number of selected species")
    ax.set_ylabel("Number of orthogroups (log10 scale)")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=numpy.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda value, _: "{:,.0f}".format(value) if value >= log_floor else "")
    )
    ax.grid(True, which="major", color="#D0D0D0", linewidth=0.45, alpha=0.75)
    ax.grid(True, which="minor", axis="y", color="#E4E4E4", linewidth=0.3, alpha=0.45)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.08, 0.815),
        frameon=False,
        ncol=1,
        labelspacing=0.35,
        handlelength=2.0,
    )
    replicates = int(summary["replicates"].iloc[0])
    fig.suptitle(
        "Orthogroup counts across random species subsets",
        x=0.08,
        y=0.985,
        ha="left",
        fontsize=8,
        fontweight="semibold",
    )
    fig.text(
        0.08,
        0.95,
        "Mean +/- 1 SD across {:,} random subsets".format(replicates),
        ha="left",
        va="top",
        fontsize=8,
        color="#4A4A4A",
    )
    fig.text(
        0.08,
        0.92,
        "Selected set: {:,} orthogroups".format(selected_total),
        ha="left",
        va="top",
        fontsize=8,
        color="#4A4A4A",
    )
    fig.text(
        0.08,
        0.89,
        "Log10 y-axis; lines stop at mean = 1",
        ha="left",
        va="top",
        fontsize=8,
        color="#4A4A4A",
    )
    fig.text(
        0.08,
        0.86,
        "SD lower bounds are clipped at 1",
        ha="left",
        va="top",
        fontsize=8,
        color="#4A4A4A",
    )
    fig.subplots_adjust(left=0.20, right=0.97, bottom=0.11, top=0.62)

    os.makedirs(outdir, exist_ok=True)
    output_paths = []
    for ext in formats:
        outpath = os.path.join(outdir, "{}.{}".format(basename, ext))
        fig.savefig(outpath, format=ext, metadata={"Creator": "GeneGalleon"})
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
    print("Loading all orthogroups: {}".format(args.orthogroup_genecount))
    print("Loading selected orthogroups: {}".format(args.selected_orthogroup_genecount))

    species_cols, counts = load_gene_count_table(args.orthogroup_genecount)
    selected_species_cols, selected_counts = load_gene_count_table(args.selected_orthogroup_genecount)
    if species_cols != selected_species_cols:
        raise ValueError("All and selected tables must have identical species columns in the same order.")
    species_counts = parse_species_counts(args.species_counts, len(species_cols))
    formats = parse_formats(args.formats)
    print("Number of orthogroups: {:,}".format(counts.shape[0]))
    print("Number of selected orthogroups: {:,}".format(selected_counts.shape[0]))
    print("Number of species: {:,}".format(len(species_cols)))
    print("Species-count values: {}".format(",".join(str(x) for x in species_counts)))
    print("Replicates per species-count value: {:,}".format(args.replicates))
    print("Random seed: {}".format(args.seed))

    values = calculate_decay(counts, species_counts, args.replicates, args.seed, selected_counts=selected_counts)
    summary = summarize_decay(species_counts, values)
    os.makedirs(args.outdir, exist_ok=True)
    summary_path = os.path.join(args.outdir, args.summary_name)
    print("Writing: {}".format(summary_path))
    summary.to_csv(summary_path, sep="\t", index=False)
    for outpath in plot_decay(
        summary,
        args.outdir,
        args.plot_basename,
        formats,
        selected_total=selected_counts.shape[0],
    ):
        print("Writing: {}".format(outpath))

    print(
        "Ending {} at {}. Elapsed time: {:,} sec".format(sys.argv[0], datetime.datetime.now(), int(time.time() - start))
    )


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
