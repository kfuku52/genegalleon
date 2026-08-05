#!/usr/bin/env python3
# coding: utf-8

import argparse
import os
import re
import sqlite3

import numpy as np
import pandas as pd

FDR_PRIORITY = [
    "q_rate_enrichment_empirical_maxT_global",
    "q_rate_enrichment_empirical_global",
    "q_rate_enrichment_global",
    "q_rate_enrichment_empirical_maxT",
    "q_rate_enrichment_empirical",
    "q_rate_enrichment",
]
P_PRIORITY = [
    "p_rate_enrichment_empirical_maxT",
    "p_rate_enrichment_empirical",
    "p_rate_enrichment",
]
PVALUE_QVALUE_METHODS = [
    {
        "label": "Analytical",
        "short_label": "Analytical",
        "p_column": "p_rate_enrichment",
        "q_column": "q_rate_enrichment_global",
        "color": "#2F6B9A",
        "linestyle": "-",
    },
    {
        "label": "Empirical (candidate-level)",
        "short_label": "Empirical",
        "p_column": "p_rate_enrichment_empirical",
        "q_column": "q_rate_enrichment_empirical_global",
        "color": "#D17A22",
        "linestyle": "--",
    },
    {
        "label": "Empirical maxT (full-scan)",
        "short_label": "Empirical maxT",
        "p_column": "p_rate_enrichment_empirical_maxT",
        "q_column": "q_rate_enrichment_empirical_maxT_global",
        "color": "#6B7F2A",
        "linestyle": ":",
    },
]
PROBABILITY_COUNT_THRESHOLDS = (0.05, 0.01, 0.001)
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")
DISPLAY_COLUMNS = [
    "orthogroup",
    "trait",
    "state_change",
    "target_class",
    "scan_match",
    "scan_rate_exposure",
    "codon_site_alignment",
    "site",
    "site_rate",
    "site_rate_categorized",
    "site_rate_quantile",
    "from_state",
    "to_state",
    "support_unit_count",
    "unit_total",
    "support_fraction",
    "support_pp_sum",
    "support_pp_mean",
    "support_unit_ids",
    "support_branch_ids",
    "rate_ratio",
    "p_rate_enrichment",
    "p_rate_enrichment_empirical",
    "p_rate_enrichment_empirical_maxT",
    "q_rate_enrichment_global",
    "q_rate_enrichment_empirical_global",
    "q_rate_enrichment_empirical_maxT_global",
    "q_rate_enrichment_empirical",
    "q_rate_enrichment_empirical_by_trait",
    "q_rate_enrichment_empirical_by_trait_match",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Write interpretable CSUBST scan aa_change summaries from a "
            "GeneGalleon SQLite database."
        )
    )
    parser.add_argument("--dbpath", metavar="PATH", required=True, help="Gene-family SQLite database path.")
    parser.add_argument("--out_tsv", metavar="PATH", required=True, help="Output ranked TSV path.")
    parser.add_argument(
        "--out_prefix",
        metavar="PATH_PREFIX",
        required=True,
        help=(
            "Prefix for plot PDFs. Writes *_evidence_density.pdf, "
            "*_substitution_spectrum.pdf, *_foreground_unit_support_matrix.pdf, "
            "and *_pvalue_qvalue_distributions.pdf."
        ),
    )
    parser.add_argument("--table", metavar="NAME", default="aa_change", help="Database table name. Default: aa_change.")
    parser.add_argument(
        "--top_n",
        metavar="INT",
        default=30,
        type=int,
        help="Number of top candidates to show in the foreground-unit support matrix.",
    )
    return parser.parse_args()


def validate_table_name(table):
    if re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", str(table)) is None:
        raise ValueError(f"Unsafe table name: {table}")


def table_exists(conn, table):
    validate_table_name(table)
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
        (table,),
    ).fetchone()
    return row is not None


def read_table(conn, table):
    validate_table_name(table)
    columns = [row[1] for row in conn.execute(f'PRAGMA table_info("{table}")').fetchall()]
    selected = [col for col in DISPLAY_COLUMNS if col in columns]
    for col in FDR_PRIORITY + P_PRIORITY:
        if col in columns and col not in selected:
            selected.append(col)
    selected.extend(col for col in columns if col not in selected)
    if not selected:
        return pd.DataFrame()
    query = "SELECT {} FROM {}".format(
        ", ".join(f'"{col}"' for col in selected),
        f'"{table}"',
    )
    return pd.read_sql_query(query, conn)


def choose_score_column(df):
    for col in FDR_PRIORITY:
        if col in df.columns and pd.to_numeric(df[col], errors="coerce").notna().any():
            return col, "FDR"
    for col in P_PRIORITY:
        if col in df.columns and pd.to_numeric(df[col], errors="coerce").notna().any():
            return col, "P"
    return None, ""


def negative_log10(series):
    values = pd.to_numeric(series, errors="coerce")
    positive = values[values > 0]
    floor_value = positive.min() / 10 if not positive.empty else 1e-300
    return -np.log10(values.fillna(1.0).clip(lower=floor_value, upper=1.0))


def numeric_column(df, column, default=np.nan):
    if column in df.columns:
        return pd.to_numeric(df[column], errors="coerce")
    return pd.Series(default, index=df.index, dtype="float64")


def ranked_candidates(df):
    if df.empty:
        return df, None, ""
    score_col, score_kind = choose_score_column(df)
    if score_col is None:
        return df, None, ""
    out = df.copy()
    out["_score"] = pd.to_numeric(out[score_col], errors="coerce")
    out["_p_rate_enrichment"] = numeric_column(out, "p_rate_enrichment")
    out["_support_fraction"] = numeric_column(out, "support_fraction")
    out = out.sort_values(
        ["_score", "_p_rate_enrichment", "_support_fraction"],
        ascending=[True, True, False],
        na_position="last",
        kind="mergesort",
    )
    internal_cols = ["_score", "_p_rate_enrichment", "_support_fraction"]
    return out.drop(columns=[col for col in internal_cols if col in out.columns]), score_col, score_kind


def prepare_plot_data(df, score_col):
    out = df.copy()
    out["_score"] = pd.to_numeric(out[score_col], errors="coerce")
    out["_neg_log10"] = negative_log10(out["_score"])
    out["_support_fraction"] = numeric_column(out, "support_fraction", 0).fillna(0)
    out["_support_unit_count"] = numeric_column(out, "support_unit_count", 0).fillna(0)
    return out


def plot_paths(out_prefix):
    return {
        "evidence_density": f"{out_prefix}_evidence_density.pdf",
        "substitution_spectrum": f"{out_prefix}_substitution_spectrum.pdf",
        "foreground_unit_support_matrix": f"{out_prefix}_foreground_unit_support_matrix.pdf",
        "pvalue_qvalue_distributions": f"{out_prefix}_pvalue_qvalue_distributions.pdf",
    }


def ensure_parent(path):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)


def write_empty_plot(out_pdf, message):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ensure_parent(out_pdf)
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


def write_empty_plot_set(paths, message):
    for out_pdf in paths.values():
        write_empty_plot(out_pdf, message)


def finite_probability_values(df, column):
    if column not in df.columns:
        return np.array([], dtype=np.float64)
    values = pd.to_numeric(df[column], errors="coerce").to_numpy(dtype=np.float64)
    return values[np.isfinite(values) & (values >= 0.0) & (values <= 1.0)]


def probability_count_label(method, values):
    counts = [int((values <= threshold).sum()) for threshold in PROBABILITY_COUNT_THRESHOLDS]
    return f"{method['short_label']}: " + " / ".join(f"{count:,}" for count in counts)


def probability_series(df, column_key):
    series = []
    for method in PVALUE_QVALUE_METHODS:
        column = method[column_key]
        values = finite_probability_values(df, column)
        if values.size == 0:
            continue
        series.append((method, column, values))
    return series


def probability_floor(series):
    positive_minima = [values[values > 0].min() for _, _, values in series if (values > 0).any()]
    if not positive_minima:
        return 1e-6
    return max(float(min(positive_minima)) / 10.0, 1e-300)


def style_probability_axis(ax):
    ax.grid(True, which="major", color="#E5E7EB", linewidth=0.75)
    ax.grid(True, which="minor", color="#F3F4F6", linewidth=0.4)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)


def plot_probability_cdf(ax, series, value_label, panel_label):
    if not series:
        ax.axis("off")
        ax.text(0.5, 0.5, f"No finite {value_label} values were available.", ha="center", va="center")
        return
    floor_value = probability_floor(series)
    x_min = 10 ** np.floor(np.log10(floor_value))
    thresholds = np.geomspace(x_min, 1.0, 700)
    min_fraction = 1.0
    point_masses = []
    for method, _, values in series:
        ordered = np.sort(np.clip(values, x_min, 1.0))
        cumulative = np.searchsorted(ordered, thresholds, side="right") / ordered.size
        positive = cumulative > 0
        min_fraction = min(min_fraction, float(cumulative[positive].min()))
        ax.step(
            thresholds[positive],
            cumulative[positive] * 100.0,
            where="post",
            color=method["color"],
            linestyle=method["linestyle"],
            linewidth=2.0,
            label=probability_count_label(method, values),
        )
        if positive.sum() == 1:
            point_masses.append(
                (
                    method,
                    float(thresholds[positive][0]),
                    float(cumulative[positive][0] * 100.0),
                )
            )
    ax.axvline(0.05, color="#6B7280", linestyle=(0, (4, 3)), linewidth=1.0)
    ax.set_xscale("log")
    ax.set_xlim(x_min, 1.18)
    if min_fraction >= 0.01:
        ax.set_ylim(0.0, 105.0)
    else:
        ax.set_yscale("log")
        ax.set_ylim(max(min_fraction * 80.0, 1e-5), 125.0)
    for method, point_x, point_y in point_masses:
        ax.vlines(
            point_x,
            ax.get_ylim()[0],
            point_y,
            color=method["color"],
            linestyle=method["linestyle"],
            linewidth=2.7,
            zorder=4,
        )
    ax.set_xlabel(f"{value_label} threshold")
    ax.set_ylabel("Candidates at or below threshold (%)")
    ax.set_title(f"{panel_label}  Cumulative {value_label} distribution", loc="left", fontsize=12, weight="bold")
    ax.legend(
        loc="upper left",
        frameon=False,
        fontsize=7.8,
        handlelength=2.7,
        title="Candidate n <= 0.05 / 0.01 / 0.001",
        title_fontsize=7.8,
    )
    style_probability_axis(ax)


def plot_probability_histogram(ax, series, value_label, panel_label):
    if not series:
        ax.axis("off")
        ax.text(0.5, 0.5, f"No finite {value_label} values were available.", ha="center", va="center")
        return
    floor_value = probability_floor(series)
    max_evidence = max(
        float((-np.log10(np.clip(values, floor_value, 1.0))).max())
        for _, _, values in series
    )
    upper = max(3.2, min(16.0, np.ceil(max_evidence * 5.0) / 5.0))
    bins = np.linspace(0.0, upper, 86)
    minimum_nonzero_percent = 100.0
    for method, _, values in series:
        evidence = -np.log10(np.clip(values, floor_value, 1.0))
        evidence = np.clip(evidence, 0.0, np.nextafter(upper, 0.0))
        counts, _ = np.histogram(evidence, bins=bins)
        nonzero = counts[counts > 0]
        if nonzero.size > 0:
            minimum_nonzero_percent = min(
                minimum_nonzero_percent,
                float(nonzero.min()) / values.size * 100.0,
            )
        ax.hist(
            evidence,
            bins=bins,
            weights=np.full(evidence.shape, 100.0 / evidence.size),
            histtype="step",
            linewidth=1.9,
            color=method["color"],
            linestyle=method["linestyle"],
            label=method["label"],
        )
    ax.axvline(-np.log10(0.05), color="#6B7280", linestyle=(0, (4, 3)), linewidth=1.0)
    ax.set_yscale("log")
    ax.set_xlim(0.0, upper)
    ax.set_ylim(max(minimum_nonzero_percent * 0.7, 1e-5), 125.0)
    ax.set_xlabel(f"Evidence strength, -log10({value_label})")
    ax.set_ylabel("Candidates per bin (%)")
    ax.set_title(f"{panel_label}  {value_label} evidence strength", loc="left", fontsize=12, weight="bold")
    ax.legend(loc="upper right", frameon=False, fontsize=8.2, handlelength=3.0)
    style_probability_axis(ax)


def write_pvalue_qvalue_distributions(df, out_pdf):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ensure_parent(out_pdf)
    p_series = probability_series(df, "p_column")
    q_series = probability_series(df, "q_column")
    if not p_series and not q_series:
        write_empty_plot(out_pdf, "No finite P-value or global q-value columns were available.")
        return

    with plt.rc_context(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.edgecolor": "#4B5563",
            "axes.linewidth": 0.8,
            "xtick.color": "#374151",
            "ytick.color": "#374151",
            "text.color": "#1F2937",
        }
    ):
        fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.4))
        plot_probability_cdf(axes[0, 0], p_series, "P-value", "A")
        plot_probability_cdf(axes[0, 1], q_series, "global q-value", "B")
        plot_probability_histogram(axes[1, 0], p_series, "P", "C")
        plot_probability_histogram(axes[1, 1], q_series, "global q", "D")
        fig.suptitle(
            "CSUBST scan P-value and global q-value distributions",
            x=0.07,
            y=0.985,
            ha="left",
            fontsize=16,
            weight="bold",
            color="#111827",
        )
        fig.text(
            0.07,
            0.952,
            (
                f"All {df.shape[0]:,} candidate state-change rows; "
                "global q-values use BH-FDR across all gene families"
            ),
            ha="left",
            va="top",
            fontsize=9.5,
            color="#4B5563",
        )
        fig.text(
            0.07,
            0.018,
            (
                "Colors identify the P/q method pair. Empirical values are discrete at the permutation resolution; "
                "the maxT method compares against the strongest null candidate per permutation."
            ),
            ha="left",
            va="bottom",
            fontsize=8.2,
            color="#4B5563",
        )
        fig.subplots_adjust(left=0.08, right=0.985, top=0.90, bottom=0.095, hspace=0.34, wspace=0.25)
        fig.savefig(out_pdf)
        plt.close(fig)


def write_evidence_density(df, score_col, out_pdf):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ensure_parent(out_pdf)
    plot_df = prepare_plot_data(df, score_col)
    plot_df = plot_df[plot_df["_support_fraction"].notna() & plot_df["_neg_log10"].notna()].copy()
    if plot_df.empty:
        write_empty_plot(out_pdf, "No finite support/FDR values for evidence-density plot.")
        return

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    if plot_df.shape[0] >= 20:
        hb = ax.hexbin(
            plot_df["_support_fraction"],
            plot_df["_neg_log10"],
            gridsize=30,
            mincnt=1,
            cmap="Blues",
        )
        fig.colorbar(hb, ax=ax, pad=0.01).set_label("candidate count")
    else:
        ax.scatter(
            plot_df["_support_fraction"],
            plot_df["_neg_log10"],
            s=40,
            color="#5AA6CF",
            edgecolor="white",
            linewidth=0.45,
            alpha=0.85,
        )

    top_df = plot_df.sort_values("_neg_log10", ascending=False).head(16)
    ax.scatter(
        top_df["_support_fraction"],
        top_df["_neg_log10"],
        s=34,
        color="#D95F02",
        edgecolor="white",
        linewidth=0.45,
        zorder=3,
    )
    ax.set_xlabel("Foreground support fraction")
    ax.set_ylabel(f"-log10({score_col})")
    ax.set_title("CSUBST scan evidence density", loc="left", fontsize=13, weight="bold")
    ax.grid(color="#dddddd", linestyle=":", linewidth=0.7)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


def state_order(values):
    observed = [str(v) for v in values if pd.notna(v) and str(v) != ""]
    extra = sorted(v for v in set(observed) if v not in AA_ORDER)
    return extra + [aa for aa in AA_ORDER if aa in set(observed)]


def write_substitution_spectrum(df, out_pdf):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ensure_parent(out_pdf)
    if "from_state" not in df.columns or "to_state" not in df.columns:
        write_empty_plot(out_pdf, "from_state/to_state columns are unavailable.")
        return
    plot_df = df[["from_state", "to_state"]].dropna().copy()
    plot_df["from_state"] = plot_df["from_state"].astype(str)
    plot_df["to_state"] = plot_df["to_state"].astype(str)
    plot_df = plot_df[(plot_df["from_state"] != "") & (plot_df["to_state"] != "")]
    if plot_df.empty:
        write_empty_plot(out_pdf, "No state-change values were available.")
        return

    rows = state_order(plot_df["from_state"])
    cols = state_order(plot_df["to_state"])
    count = (
        plot_df.groupby(["from_state", "to_state"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=rows, columns=cols, fill_value=0)
    )

    width = max(5.8, min(10.0, 0.36 * len(cols) + 2.6))
    height = max(4.6, min(10.0, 0.34 * len(rows) + 2.5))
    fig, ax = plt.subplots(figsize=(width, height))
    im = ax.imshow(count.values, cmap="YlGnBu", aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows)
    ax.set_xlabel("To state")
    ax.set_ylabel("From state")
    ax.set_title("CSUBST scan substitution spectrum", loc="left", fontsize=13, weight="bold")
    threshold = max(2, int(np.nanpercentile(count.values[count.values > 0], 80))) if (count.values > 0).any() else 2
    for i in range(len(rows)):
        for j in range(len(cols)):
            value = int(count.values[i, j])
            if value >= threshold:
                ax.text(j, i, str(value), ha="center", va="center", fontsize=7, color="black")
    fig.colorbar(im, ax=ax, pad=0.01).set_label("candidate count")
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


def parse_unit_ids(value):
    if pd.isna(value):
        return []
    parts = re.split(r"[,;|\s]+", str(value).strip())
    return [part for part in parts if part]


def candidate_label(row):
    orthogroup = str(row["orthogroup"]) if "orthogroup" in row and pd.notna(row["orthogroup"]) else "family"
    state_change = str(row["state_change"]) if "state_change" in row and pd.notna(row["state_change"]) else "candidate"
    trait = str(row["trait"]) if "trait" in row and pd.notna(row["trait"]) else ""
    label = f"{orthogroup} {state_change}"
    if trait:
        label = f"{label} ({trait})"
    return label


def write_foreground_unit_support_matrix(df, score_col, out_pdf, top_n):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    ensure_parent(out_pdf)
    if "support_unit_ids" not in df.columns:
        write_empty_plot(out_pdf, "support_unit_ids column is unavailable.")
        return
    plot_df = prepare_plot_data(df, score_col)
    plot_df["_unit_list"] = plot_df["support_unit_ids"].map(parse_unit_ids)
    plot_df = plot_df[plot_df["_unit_list"].map(len) > 0].copy()
    if plot_df.empty:
        write_empty_plot(out_pdf, "No foreground-unit support values were available.")
        return

    plot_df = plot_df.sort_values(["_neg_log10", "_support_fraction"], ascending=[False, False]).head(max(1, top_n))
    unit_ids = sorted(
        {unit for units in plot_df["_unit_list"] for unit in units},
        key=lambda item: (0, int(item)) if str(item).isdigit() else (1, str(item)),
    )
    mat = np.zeros((plot_df.shape[0], len(unit_ids)), dtype=float)
    unit_index = {unit: idx for idx, unit in enumerate(unit_ids)}
    for row_idx, units in enumerate(plot_df["_unit_list"]):
        for unit in units:
            if unit in unit_index:
                mat[row_idx, unit_index[unit]] = 1.0

    if mat.shape[1] > 1:
        weights = np.linspace(1.0, 2.0, mat.shape[0])
        col_order = np.argsort(mat.T @ weights)[::-1]
        mat = mat[:, col_order]
        unit_ids = [unit_ids[idx] for idx in col_order]

    labels = [candidate_label(row) for _, row in plot_df.iterrows()]
    height = max(4.4, min(16.0, 1.4 + 0.24 * len(labels)))
    width = max(6.0, min(14.0, 2.8 + 0.4 * len(unit_ids)))
    fig, ax = plt.subplots(figsize=(width, height))
    ax.imshow(mat, aspect="auto", cmap=ListedColormap(["#f4f4f4", "#2166ac"]))
    ax.set_xticks(range(len(unit_ids)))
    ax.set_xticklabels(unit_ids)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Foreground unit")
    ax.set_title("CSUBST scan foreground-unit support", loc="left", fontsize=13, weight="bold")
    ax.tick_params(length=0)
    for x in np.arange(-0.5, len(unit_ids), 1):
        ax.axvline(x, color="white", linewidth=0.6)
    for y in np.arange(-0.5, len(labels), 1):
        ax.axhline(y, color="white", linewidth=0.6)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)


def main():
    args = parse_args()
    paths = plot_paths(args.out_prefix)
    ensure_parent(args.out_tsv)

    if not os.path.exists(args.dbpath):
        pd.DataFrame().to_csv(args.out_tsv, sep="\t", index=False)
        write_empty_plot_set(paths, f"Database not found: {args.dbpath}")
        return 0

    with sqlite3.connect(args.dbpath) as conn:
        if not table_exists(conn, args.table):
            pd.DataFrame().to_csv(args.out_tsv, sep="\t", index=False)
            write_empty_plot_set(paths, f"Table not found: {args.table}")
            return 0
        df = read_table(conn, args.table)

    ranked, score_col, score_kind = ranked_candidates(df)
    ranked.to_csv(args.out_tsv, sep="\t", index=False)
    if ranked.empty:
        write_empty_plot_set(paths, "No CSUBST scan candidates in aa_change.")
    elif score_col is None:
        write_empty_plot_set(paths, "No CSUBST scan P-value or FDR columns were available.")
    else:
        write_evidence_density(ranked, score_col, paths["evidence_density"])
        write_substitution_spectrum(ranked, paths["substitution_spectrum"])
        write_foreground_unit_support_matrix(ranked, score_col, paths["foreground_unit_support_matrix"], args.top_n)
        write_pvalue_qvalue_distributions(ranked, paths["pvalue_qvalue_distributions"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
