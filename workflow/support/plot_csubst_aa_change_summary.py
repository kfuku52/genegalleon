#!/usr/bin/env python3
# coding: utf-8

import argparse
import math
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
            "*_substitution_spectrum.pdf, and *_foreground_unit_support_matrix.pdf."
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
