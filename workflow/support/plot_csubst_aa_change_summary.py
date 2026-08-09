#!/usr/bin/env python3
# coding: utf-8

import argparse
import os
import re
import sqlite3
import sys
from pathlib import Path

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
SUPPORT_SIGNIFICANCE_THRESHOLD = 0.05
SUPPORT_BIN_COUNT = 10
PRIMARY_MIN_SUPPORT = 2
MIN_SUPPORT_SENSITIVITY_START = 3
ORTHOGROUP_BESTHIT_COLUMNS = [
    "besthit_0.05",
    "besthit_0.25",
    "besthit_0.5",
    "besthit_0.75",
    "besthit_0.95",
]
GLOBAL_QVALUE_COLUMNS = {
    "p_rate_enrichment": "q_rate_enrichment_global",
    "p_rate_enrichment_empirical": "q_rate_enrichment_empirical_global",
    "p_rate_enrichment_empirical_maxT": "q_rate_enrichment_empirical_maxT_global",
}
GROUPED_QVALUE_COLUMNS = {
    "p_rate_enrichment": {
        "q_rate_enrichment": ("orthogroup",),
        "q_rate_enrichment_by_trait": ("orthogroup", "trait"),
        "q_rate_enrichment_by_trait_match": ("orthogroup", "trait", "scan_match"),
    },
    "p_rate_enrichment_empirical": {
        "q_rate_enrichment_empirical": ("orthogroup",),
        "q_rate_enrichment_empirical_by_trait": ("orthogroup", "trait"),
        "q_rate_enrichment_empirical_by_trait_match": ("orthogroup", "trait", "scan_match"),
    },
}
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
            "Prefix for plot PDFs. Writes min_support_2 primary outputs plus "
            "min_support_3 and higher sensitivity outputs in the same directory."
        ),
    )
    parser.add_argument(
        "--orthogroup_annotation_tsv",
        metavar="PATH",
        help=(
            "Optional Orthogroups.GeneCount.annotated.tsv path. When supplied, "
            "the five representative besthit columns are joined by orthogroup."
        ),
    )
    parser.add_argument("--table", metavar="NAME", default="aa_change", help="Database table name. Default: aa_change.")
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


def attach_orthogroup_besthits(df, annotation_path):
    if not annotation_path:
        return df

    annotation_path = Path(annotation_path)
    if not annotation_path.is_file():
        print(
            f"Warning: orthogroup annotation table was not found: {annotation_path}. "
            "Writing CSUBST summaries without besthit columns.",
            file=sys.stderr,
            flush=True,
        )
        return df
    if "orthogroup" not in df.columns:
        raise ValueError(
            "Cannot join orthogroup annotations because the CSUBST summary has no "
            "'orthogroup' column."
        )

    annotation_key = "Orthogroup"
    required_columns = [annotation_key, *ORTHOGROUP_BESTHIT_COLUMNS]
    header = pd.read_csv(annotation_path, sep="\t", nrows=0).columns.tolist()
    missing_columns = [column for column in required_columns if column not in header]
    if missing_columns:
        raise ValueError(
            f"Orthogroup annotation table is missing required columns: "
            f"{', '.join(missing_columns)}"
        )

    annotation = pd.read_csv(
        annotation_path,
        sep="\t",
        usecols=required_columns,
        dtype={annotation_key: "string"},
        low_memory=False,
    )
    annotation[annotation_key] = annotation[annotation_key].str.strip()
    missing_keys = annotation[annotation_key].isna() | annotation[annotation_key].eq("")
    if missing_keys.any():
        raise ValueError(
            "Orthogroup annotation table contains "
            f"{int(missing_keys.sum()):,} blank Orthogroup key(s)."
        )
    duplicate_keys = annotation[annotation_key].duplicated(keep=False)
    if duplicate_keys.any():
        examples = annotation.loc[duplicate_keys, annotation_key].drop_duplicates().head(5)
        raise ValueError(
            "Orthogroup annotation table contains duplicate Orthogroup keys: "
            f"{', '.join(examples.astype(str))}"
        )

    join_key = "_orthogroup_annotation_key"
    summary_keys = df["orthogroup"].astype("string").str.strip()
    unique_summary_keys = summary_keys.dropna()
    unique_summary_keys = unique_summary_keys[unique_summary_keys.ne("")].drop_duplicates()
    if not unique_summary_keys.empty:
        matched = unique_summary_keys.isin(annotation[annotation_key])
        matched_count = int(matched.sum())
        total_count = int(unique_summary_keys.shape[0])
        if matched_count == 0:
            raise ValueError(
                "Orthogroup annotation join matched 0 of "
                f"{total_count:,} CSUBST orthogroups. Check for an OG/HOG identifier mismatch."
            )
        coverage = matched_count / total_count
        if matched_count < total_count:
            print(
                "Warning: orthogroup annotation join matched "
                f"{matched_count:,} of {total_count:,} unique CSUBST orthogroups "
                f"({coverage:.1%}); unmatched rows retain blank besthit values.",
                file=sys.stderr,
                flush=True,
            )
        else:
            print(
                "Orthogroup annotation join matched "
                f"{matched_count:,} of {total_count:,} unique CSUBST orthogroups (100.0%).",
                flush=True,
            )

    output = df.drop(columns=[col for col in ORTHOGROUP_BESTHIT_COLUMNS if col in df.columns]).copy()
    output[join_key] = summary_keys
    annotation = annotation.rename(columns={annotation_key: join_key})
    row_count = output.shape[0]
    output = output.merge(
        annotation,
        how="left",
        on=join_key,
        sort=False,
        validate="many_to_one",
    ).drop(columns=join_key)
    if output.shape[0] != row_count:
        raise RuntimeError(
            "Orthogroup annotation join changed the number of CSUBST summary rows: "
            f"{row_count:,} -> {output.shape[0]:,}."
        )

    ordered_columns = ["orthogroup", *ORTHOGROUP_BESTHIT_COLUMNS]
    ordered_columns.extend(column for column in output.columns if column not in ordered_columns)
    return output.loc[:, ordered_columns]


def choose_score_column(df):
    for col in FDR_PRIORITY:
        if col in df.columns and pd.to_numeric(df[col], errors="coerce").notna().any():
            return col, "FDR"
    for col in P_PRIORITY:
        if col in df.columns and pd.to_numeric(df[col], errors="coerce").notna().any():
            return col, "P"
    return None, ""


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


def calculate_bh_fdr(pvalues):
    pvalues = pd.to_numeric(pd.Series(pvalues), errors="coerce").to_numpy(dtype=float)
    qvalues = np.full(shape=pvalues.shape, fill_value=np.nan, dtype=float)
    finite = np.isfinite(pvalues)
    if not finite.any():
        return qvalues
    finite_index = np.flatnonzero(finite)
    finite_p = np.clip(pvalues[finite], 0.0, 1.0)
    order = np.argsort(finite_p, kind="mergesort")
    ranked = finite_p[order]
    ranks = np.arange(1, ranked.shape[0] + 1, dtype=float)
    ranked_q = ranked * ranked.shape[0] / ranks
    ranked_q = np.minimum.accumulate(ranked_q[::-1])[::-1]
    ranked_q = np.clip(ranked_q, 0.0, 1.0)
    qvalues[finite_index[order]] = ranked_q
    return qvalues


def calculate_grouped_bh_fdr(df, p_column, group_columns):
    qvalues = np.full(shape=df.shape[0], fill_value=np.nan, dtype=float)
    grouped_positions = df.groupby(list(group_columns), sort=False, dropna=False).indices.values()
    for positions in grouped_positions:
        positions = np.asarray(positions, dtype=np.int64)
        qvalues[positions] = calculate_bh_fdr(df.iloc[positions][p_column])
    return qvalues


def recalculate_sensitivity_qvalues(df):
    recalculated = []
    for p_column, q_column in GLOBAL_QVALUE_COLUMNS.items():
        if p_column not in df.columns:
            continue
        df[q_column] = calculate_bh_fdr(df[p_column])
        recalculated.append(q_column)
    for p_column, q_columns in GROUPED_QVALUE_COLUMNS.items():
        if p_column not in df.columns:
            continue
        for q_column, group_columns in q_columns.items():
            if not set(group_columns).issubset(df.columns):
                continue
            df[q_column] = calculate_grouped_bh_fdr(df, p_column, group_columns)
            recalculated.append(q_column)
    return recalculated


def min_support_sensitivity_thresholds(df, start=MIN_SUPPORT_SENSITIVITY_START):
    support = numeric_column(df, "support_unit_count")
    finite = support[np.isfinite(support)]
    if finite.empty:
        return []
    maximum = int(np.floor(finite.max()))
    if maximum < start:
        return []
    return list(range(start, maximum + 1))


def min_support_sensitivity_paths(out_prefix, threshold=None):
    prefix = Path(str(out_prefix))
    output_dir = prefix.parent
    stem = prefix.name
    paths = {
        "output_dir": output_dir,
        "manifest": output_dir / f"{stem}_min_support_manifest.tsv",
    }
    if threshold is not None:
        threshold_stem = f"{stem}_min_support_{threshold}"
        paths.update(
            {
                "summary_tsv": output_dir / f"{threshold_stem}_summary.tsv",
                "plot_pdf": output_dir / f"{threshold_stem}_pvalue_qvalue_distributions.pdf",
            }
        )
    return paths


def qvalue_manifest_metrics(values):
    values = np.asarray(values, dtype=np.float64)
    finite = values[np.isfinite(values)]
    metrics = {"min": float(finite.min()) if finite.size else np.nan}
    for threshold in PROBABILITY_COUNT_THRESHOLDS:
        metrics[f"le_{threshold}"] = int((finite <= threshold).sum())
    return metrics


def remove_stale_min_support_sensitivity_outputs(out_prefix, thresholds):
    paths = min_support_sensitivity_paths(out_prefix)
    output_dir = paths["output_dir"]
    if not output_dir.is_dir():
        return
    primary_paths = min_support_sensitivity_paths(out_prefix, PRIMARY_MIN_SUPPORT)
    expected = {primary_paths["summary_tsv"], primary_paths["plot_pdf"]}
    for threshold in thresholds:
        threshold_paths = min_support_sensitivity_paths(out_prefix, threshold)
        expected.update((threshold_paths["summary_tsv"], threshold_paths["plot_pdf"]))
    stem = Path(str(out_prefix)).name
    patterns = (
        f"{stem}_min_support_*_summary.tsv",
        f"{stem}_min_support_*_pvalue_qvalue_distributions.pdf",
    )
    for pattern in patterns:
        for candidate in output_dir.glob(pattern):
            if candidate not in expected:
                candidate.unlink()
    paths["manifest"].unlink(missing_ok=True)


def remove_legacy_min_support_output_layout(out_prefix):
    prefix = Path(str(out_prefix))
    legacy_files = (
        Path(f"{prefix}_summary.tsv"),
        Path(f"{prefix}_support_significance_rate.pdf"),
        Path(f"{prefix}_substitution_spectrum.pdf"),
        Path(f"{prefix}_pvalue_qvalue_distributions.pdf"),
    )
    for legacy_file in legacy_files:
        legacy_file.unlink(missing_ok=True)

    legacy_dir = prefix.parent / "min_support_sensitivity"
    if not legacy_dir.is_dir():
        return
    patterns = (
        f"{prefix.name}_min_support_*_summary.tsv",
        f"{prefix.name}_min_support_*_pvalue_qvalue_distributions.pdf",
        f"{prefix.name}_min_support_manifest.tsv",
    )
    for pattern in patterns:
        for candidate in legacy_dir.glob(pattern):
            candidate.unlink()
    try:
        legacy_dir.rmdir()
    except OSError:
        pass


def write_min_support_sensitivity(df, out_prefix):
    thresholds = min_support_sensitivity_thresholds(df)
    remove_stale_min_support_sensitivity_outputs(out_prefix, thresholds)
    if not thresholds:
        print(
            "Skipping min_support sensitivity outputs because support_unit_count "
            f"has no finite value >= {MIN_SUPPORT_SENSITIVITY_START}.",
            flush=True,
        )
        return None

    paths = min_support_sensitivity_paths(out_prefix)
    paths["output_dir"].mkdir(parents=True, exist_ok=True)
    support = numeric_column(df, "support_unit_count")
    manifest_rows = []
    for threshold in thresholds:
        threshold_paths = min_support_sensitivity_paths(out_prefix, threshold)
        subset = df.loc[support >= threshold].copy()
        q_columns = recalculate_sensitivity_qvalues(subset)
        subset.to_csv(threshold_paths["summary_tsv"], sep="\t", index=False)
        write_pvalue_qvalue_distributions(subset, threshold_paths["plot_pdf"])
        row = {
            "min_support": threshold,
            "candidate_rows": int(subset.shape[0]),
            "summary_tsv": threshold_paths["summary_tsv"].name,
            "plot_pdf": threshold_paths["plot_pdf"].name,
        }
        for q_column in q_columns:
            metrics = qvalue_manifest_metrics(subset[q_column])
            for metric, value in metrics.items():
                row[f"{q_column}_{metric}"] = value
        manifest_rows.append(row)
        print(f"min_support={threshold}: {subset.shape[0]:,} candidate rows", flush=True)

    pd.DataFrame(manifest_rows).to_csv(paths["manifest"], sep="\t", index=False)
    return paths["manifest"]


def plot_paths(out_prefix):
    primary_prefix = f"{out_prefix}_min_support_{PRIMARY_MIN_SUPPORT}"
    return {
        "support_significance_rate": f"{primary_prefix}_support_significance_rate.pdf",
        "substitution_spectrum": f"{primary_prefix}_substitution_spectrum.pdf",
        "pvalue_qvalue_distributions": f"{primary_prefix}_pvalue_qvalue_distributions.pdf",
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


def support_significance_data(
    df,
    bin_count=SUPPORT_BIN_COUNT,
    threshold=SUPPORT_SIGNIFICANCE_THRESHOLD,
):
    support = numeric_column(df, "support_fraction")
    support_valid = support.notna() & np.isfinite(support) & support.between(0.0, 1.0)
    edges = np.linspace(0.0, 1.0, bin_count + 1)
    centers = (edges[:-1] + edges[1:]) / 2.0
    support_bins = pd.cut(
        support,
        bins=edges,
        labels=centers,
        include_lowest=True,
        right=True,
    )
    candidate_counts = (
        support_bins[support_valid]
        .value_counts(sort=False)
        .reindex(centers, fill_value=0)
        .to_numpy(dtype=np.int64)
    )

    method_series = []
    for method in PVALUE_QVALUE_METHODS:
        column = method["q_column"]
        if column not in df.columns:
            continue
        q_values = pd.to_numeric(df[column], errors="coerce")
        q_valid = q_values.notna() & np.isfinite(q_values) & q_values.between(0.0, 1.0)
        valid = support_valid & q_valid
        if not valid.any():
            continue
        method_frame = pd.DataFrame(
            {
                "support_bin": support_bins[valid],
                "significant": q_values[valid] <= threshold,
            }
        )
        grouped = method_frame.groupby("support_bin", observed=False)["significant"].agg(["sum", "count"])
        grouped = grouped.reindex(centers, fill_value=0)
        denominators = grouped["count"].to_numpy(dtype=np.int64)
        significant_counts = grouped["sum"].to_numpy(dtype=np.int64)
        percentages = np.full(centers.shape, np.nan, dtype=np.float64)
        populated = denominators > 0
        percentages[populated] = significant_counts[populated] / denominators[populated] * 100.0
        method_series.append(
            {
                "method": method,
                "centers": centers.copy(),
                "percentages": percentages,
                "significant_counts": significant_counts,
                "candidate_counts": denominators,
            }
        )
    return centers, candidate_counts, method_series


def write_support_significance_rate(df, out_pdf):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ensure_parent(out_pdf)
    centers, candidate_counts, method_series = support_significance_data(df)
    if not method_series:
        write_empty_plot(out_pdf, "No finite foreground-support/global-q pairs were available.")
        return

    marker_shapes = ["o", "s", "^"]
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
        fig, (rate_ax, count_ax) = plt.subplots(
            2,
            1,
            figsize=(8.2, 6.2),
            sharex=True,
            gridspec_kw={"height_ratios": [3.3, 1.0], "hspace": 0.08},
        )
        maximum_rate = 0.0
        any_significant = False
        for index, item in enumerate(method_series):
            method = item["method"]
            populated = item["candidate_counts"] > 0
            if not populated.any():
                continue
            x_values = item["centers"][populated]
            y_values = item["percentages"][populated]
            significant_total = int(item["significant_counts"].sum())
            candidate_total = int(item["candidate_counts"].sum())
            any_significant = any_significant or significant_total > 0
            maximum_rate = max(maximum_rate, float(np.nanmax(y_values)))
            rate_ax.plot(
                x_values,
                y_values,
                color=method["color"],
                linestyle=method["linestyle"],
                marker=marker_shapes[index % len(marker_shapes)],
                markersize=4.5,
                linewidth=1.9,
                label=f"{method['label']} ({significant_total:,}/{candidate_total:,})",
            )

        upper_limit = 1.0 if maximum_rate <= 0 else min(100.0, maximum_rate * 1.18 + 0.1)
        rate_ax.set_ylim(0.0, upper_limit)
        rate_ax.set_ylabel(f"Candidates with global q <= {SUPPORT_SIGNIFICANCE_THRESHOLD:g} (%)")
        rate_ax.grid(True, color="#E5E7EB", linewidth=0.75)
        rate_ax.legend(
            loc="upper left",
            frameon=False,
            fontsize=8.3,
            handlelength=3.0,
            title="Significant / finite-q candidates",
            title_fontsize=8.3,
        )
        if not any_significant:
            rate_ax.text(
                0.98,
                0.08,
                f"No candidates with global q <= {SUPPORT_SIGNIFICANCE_THRESHOLD:g}",
                ha="right",
                va="bottom",
                transform=rate_ax.transAxes,
                color="#6B7280",
                fontsize=8.5,
            )

        populated_counts = candidate_counts > 0
        count_ax.bar(
            centers[populated_counts],
            candidate_counts[populated_counts],
            width=0.085,
            color="#D1D5DB",
            edgecolor="#6B7280",
            linewidth=0.6,
        )
        count_ax.set_ylabel("Candidate n")
        count_ax.set_xlabel("Foreground support fraction (0.1-wide bins)")
        count_ax.set_xlim(0.0, 1.0)
        count_ax.set_xticks(np.linspace(0.0, 1.0, 6))
        count_ax.grid(True, axis="y", color="#E5E7EB", linewidth=0.7)

        for axis in (rate_ax, count_ax):
            for spine in ["top", "right"]:
                axis.spines[spine].set_visible(False)
        fig.suptitle(
            "CSUBST scan foreground support and significance",
            x=0.09,
            y=0.985,
            ha="left",
            fontsize=14,
            weight="bold",
            color="#111827",
        )
        fig.text(
            0.09,
            0.945,
            "Global BH-FDR significance rate within each foreground-support bin",
            ha="left",
            va="top",
            fontsize=9,
            color="#4B5563",
        )
        fig.subplots_adjust(left=0.12, right=0.985, top=0.89, bottom=0.11)
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


def main():
    args = parse_args()
    remove_legacy_min_support_output_layout(args.out_prefix)
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

    df = attach_orthogroup_besthits(df, args.orthogroup_annotation_tsv)
    ranked, score_col, score_kind = ranked_candidates(df)
    ranked.to_csv(args.out_tsv, sep="\t", index=False)
    if ranked.empty:
        remove_stale_min_support_sensitivity_outputs(args.out_prefix, [])
        write_empty_plot_set(paths, "No CSUBST scan candidates in aa_change.")
    elif score_col is None:
        write_empty_plot_set(paths, "No CSUBST scan P-value or FDR columns were available.")
    else:
        write_support_significance_rate(ranked, paths["support_significance_rate"])
        write_substitution_spectrum(ranked, paths["substitution_spectrum"])
        write_pvalue_qvalue_distributions(ranked, paths["pvalue_qvalue_distributions"])
    if not ranked.empty:
        write_min_support_sensitivity(ranked, args.out_prefix)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
