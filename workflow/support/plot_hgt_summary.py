#!/usr/bin/env python3

import argparse
import math
import os
import textwrap
from typing import Dict, Iterable, List, Tuple

import numpy
import pandas

from score_hgt_candidates import TaxonomyResolver, normalize_sci_name


OVERVIEW_NUMERIC_COLUMNS: List[Tuple[str, str, str]] = [
    ("candidate_gene_count", "Cand", "Number of genes descending from the candidate GeneRax HGT branch."),
    ("matched_leaf_count", "Tips", "Number of leaf rows in the database that were matched back to the candidate branch."),
    ("besthit_gene_count", "HitGenes", "Number of candidate genes with a UniProt/Swiss-Prot best hit."),
    ("besthit_taxid_count", "HitTaxID", "Number of candidate genes whose best hit could be resolved to a taxonomy-aware label."),
    ("besthit_same_superkingdom_fraction", "SameSK", "Fraction of candidate genes whose best hit falls in the same superkingdom as the focal lineage."),
    ("intron_support_fraction", "Intron", "Fraction of candidate genes with intron support recorded in stat_branch."),
    ("expression_measured_fraction", "Expr", "Fraction of candidate genes with any expression measurement."),
    ("clade_min_expression_pearsoncor", "ExprCor", "Minimum clade-level Pearson correlation among measured expression profiles."),
    ("synteny_support_fraction", "SynFrac", "Fraction of candidate genes with positive synteny support."),
    ("synteny_mean_support_score", "SynMean", "Mean synteny support score across candidate genes."),
    ("contamination_incompatible_fraction", "Contam", "Fraction of candidate genes flagged as lineage-incompatible by contamination QC."),
]

OVERVIEW_TEXT_COLUMNS: List[Tuple[str, str, int, str]] = [
    ("besthit_lca_rank_mode", "HitLCA", 18, "Most frequent best-hit lineage relationship label among candidate genes."),
    ("contamination_top_lca_sciname", "TopContam", 28, "Most frequent contamination LCA scientific name among candidate genes."),
]

FLOW_FALLBACK_LABEL = "Unresolved"
FLOW_OTHER_LABEL = "Other"


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Plot overview and taxonomy-flow summaries for gg_hgt outputs."
    )
    parser.add_argument("--branch_tsv", metavar="PATH", required=True, type=str)
    parser.add_argument("--gene_tsv", metavar="PATH", required=True, type=str)
    parser.add_argument("--overview_pdf", metavar="PATH", required=True, type=str)
    parser.add_argument("--taxonomy_flow_pdf", metavar="PATH", required=True, type=str)
    parser.add_argument("--taxonomy_dbfile", metavar="PATH", default="", type=str)
    parser.add_argument("--flow_rank", metavar="TEXT", default="phylum", type=str)
    parser.add_argument("--flow_max_categories", metavar="INT", default=12, type=int)
    return parser


def get_pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch, Rectangle

    matplotlib.rcParams["font.size"] = 8
    matplotlib.rcParams["font.family"] = "Helvetica"
    matplotlib.rcParams["svg.fonttype"] = "none"
    return plt, PdfPages, Path, PathPatch, Rectangle


def ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def write_overview_readme(out_pdf: str) -> None:
    out_dir = os.path.dirname(out_pdf)
    if not out_dir:
        return
    ensure_parent_dir(out_pdf)
    readme_path = os.path.join(out_dir, "README.md")
    lines = [
        "# HGT Overview Metrics",
        "",
        "`hgt_branch_overview.pdf` shows one row per `orthogroup:branch_id` candidate.",
        "Numeric columns are normalized independently within each PDF page before coloring, so color intensity is only comparable within the same column.",
        "",
        "## Numeric Columns",
        "",
    ]
    for _col, short_label, description in OVERVIEW_NUMERIC_COLUMNS:
        lines.append(f"- `{short_label}`: {description}")
    lines.extend(["", "## Text Columns", ""])
    for _col, short_label, _width, description in OVERVIEW_TEXT_COLUMNS:
        lines.append(f"- `{short_label}`: {description}")
    lines.extend(
        [
            "",
            "## Row Labels",
            "",
            "- `OGXXXX:branch_id`: orthogroup ID and the branch ID carried through `stat_branch` / `gg_orthogroup.db`.",
        ]
    )
    with open(readme_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def safe_read_tsv(path: str) -> pandas.DataFrame:
    if not path or not os.path.exists(path):
        return pandas.DataFrame()
    try:
        return pandas.read_csv(path, sep="\t", low_memory=False)
    except pandas.errors.EmptyDataError:
        return pandas.DataFrame()


def shorten_text(value, max_nchar: int) -> str:
    if pandas.isna(value):
        return ""
    text = str(value).strip().replace("_", " ")
    if text.lower() == "nan":
        return ""
    text = " ".join(text.split())
    if len(text) <= max_nchar:
        return text
    if max_nchar <= 3:
        return text[:max_nchar]
    return text[: max_nchar - 3] + "..."


def wrap_label(value: str, width: int = 14) -> str:
    text = shorten_text(value, 40)
    if text == "":
        return FLOW_FALLBACK_LABEL
    return "\n".join(textwrap.wrap(text, width=width))


def value_to_boolish(value):
    if pandas.isna(value):
        return pandas.NA
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y"}:
        return True
    if text in {"0", "false", "f", "no", "n"}:
        return False
    return pandas.NA


def normalize_numeric_frame(df: pandas.DataFrame) -> pandas.DataFrame:
    out = df.copy()
    for col in out.columns:
        vals = pandas.to_numeric(out[col], errors="coerce")
        finite = vals[numpy.isfinite(vals)]
        if finite.empty:
            out[col] = numpy.nan
            continue
        vmin = float(finite.min())
        vmax = float(finite.max())
        if math.isclose(vmin, vmax):
            norm = vals.copy() * numpy.nan
            norm[vals.notna()] = 0.5
        else:
            norm = (vals - vmin) / (vmax - vmin)
        out[col] = norm
    return out


def blank_pdf(path: str, title: str, message: str) -> None:
    plt, PdfPages, _, _, _ = get_pyplot()
    ensure_parent_dir(path)
    with PdfPages(path) as pdf:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.axis("off")
        ax.text(0.5, 0.65, title, ha="center", va="center", fontsize=12, fontweight="bold")
        ax.text(0.5, 0.40, message, ha="center", va="center", fontsize=9)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)


def plot_overview(branch_df: pandas.DataFrame, out_pdf: str) -> None:
    plt, PdfPages, _, _, _ = get_pyplot()
    ensure_parent_dir(out_pdf)
    write_overview_readme(out_pdf)
    if branch_df.empty:
        blank_pdf(out_pdf, "HGT Branch Overview", "No HGT candidate branches were found.")
        return

    plot_df = branch_df.copy()
    plot_df["row_label"] = (
        plot_df["orthogroup"].astype(str).str.strip()
        + ":"
        + plot_df["branch_id"].astype(str).str.strip()
    )
    plot_df = plot_df.sort_values(["orthogroup", "branch_id"], ascending=[True, True], kind="mergesort").reset_index(drop=True)
    numeric_cols = [col for col, _label, _desc in OVERVIEW_NUMERIC_COLUMNS if col in plot_df.columns]
    text_specs = [(col, label, width) for col, label, width, _desc in OVERVIEW_TEXT_COLUMNS if col in plot_df.columns]
    chunk_size = 60

    with PdfPages(out_pdf) as pdf:
        for start in range(0, plot_df.shape[0], chunk_size):
            chunk = plot_df.iloc[start:start + chunk_size, :].copy()
            numeric_df = chunk.loc[:, numeric_cols].apply(pandas.to_numeric, errors="coerce")
            normalized = normalize_numeric_frame(numeric_df) if len(numeric_cols) > 0 else pandas.DataFrame(index=chunk.index)
            row_labels = chunk["row_label"].astype(str).tolist()

            fig_height = max(4.0, 0.24 * len(row_labels) + 1.8)
            text_width_units = max(1, len(text_specs))
            fig, axes = plt.subplots(
                1,
                1 + text_width_units,
                figsize=(max(12.0, 0.7 * max(1, len(numeric_cols)) + 2.8 * text_width_units), fig_height),
                gridspec_kw={"width_ratios": [max(2.3, 0.78 * max(1, len(numeric_cols)))] + [2.6] * text_width_units},
                squeeze=False,
            )
            ax_heat = axes[0, 0]
            if len(numeric_cols) > 0:
                matrix = normalized.to_numpy(dtype=float)
                masked = numpy.ma.masked_invalid(matrix)
                cmap = plt.get_cmap("viridis").copy()
                cmap.set_bad(color="#f2f2f2")
                im = ax_heat.imshow(masked, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
                ax_heat.figure.colorbar(im, ax=ax_heat, fraction=0.025, pad=0.02, label="Column-normalized")
                ax_heat.set_xticks(range(len(numeric_cols)))
                ax_heat.set_xticklabels(
                    [label for col, label, _desc in OVERVIEW_NUMERIC_COLUMNS if col in numeric_cols],
                    rotation=45,
                    ha="right",
                )
            else:
                ax_heat.imshow(numpy.zeros((len(row_labels), 1)), aspect="auto", interpolation="nearest", cmap="Greys")
                ax_heat.set_xticks([0])
                ax_heat.set_xticklabels(["No numeric columns"])
            ax_heat.set_yticks(range(len(row_labels)))
            ax_heat.set_yticklabels(row_labels)
            ax_heat.set_title("HGT Branch Overview")
            ax_heat.set_xlabel("Evidence")
            ax_heat.tick_params(axis="y", labelsize=7)

            for idx, (col, label, width) in enumerate(text_specs, start=1):
                ax_text = axes[0, idx]
                ax_text.set_xlim(0, 1)
                ax_text.set_ylim(len(row_labels) - 0.5, -0.5)
                ax_text.axis("off")
                ax_text.set_title(label)
                values = chunk[col].map(lambda x: shorten_text(x, width)).fillna("")
                for y, value in enumerate(values.tolist()):
                    ax_text.text(0.0, y, value, va="center", ha="left", fontsize=7, clip_on=False)

            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


def resolve_rank_label(
    tax_name: str,
    taxid_value,
    resolver: TaxonomyResolver,
    preferred_rank: str,
) -> str:
    preferred_rank = str(preferred_rank).strip().lower()
    taxid = 0
    if taxid_value not in ("", None) and not pandas.isna(taxid_value):
        try:
            taxid = int(float(taxid_value))
        except (TypeError, ValueError):
            taxid = 0
    if taxid <= 0 and resolver.enabled:
        taxid = resolver.resolve_name_taxid(tax_name)
    if taxid > 0 and resolver.enabled:
        lineage = resolver.lineage(taxid)
        if len(lineage) > 0:
            ordered_ranks = [preferred_rank, "subphylum", "phylum", "class", "order", "kingdom", "superkingdom"]
            chosen_taxid = 0
            for rank_name in ordered_ranks:
                chosen_taxid = resolver.rank_taxid_from_lineage(lineage, [rank_name])
                if chosen_taxid > 0:
                    break
            if chosen_taxid == 0:
                chosen_taxid = lineage[-1]
            if chosen_taxid > 0 and resolver.ncbi is not None:
                try:
                    names = resolver.ncbi.get_taxid_translator([chosen_taxid])
                    if chosen_taxid in names:
                        return str(names[chosen_taxid])
                except Exception:
                    pass
    normalized = normalize_sci_name(tax_name)
    if normalized == "":
        return FLOW_FALLBACK_LABEL
    parts = normalized.split()
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    return parts[0]


def collapse_to_top_categories(count_df: pandas.DataFrame, left_col: str, right_col: str, max_categories: int) -> pandas.DataFrame:
    out = count_df.copy()
    left_totals = out.groupby(left_col, sort=False)["count"].sum().sort_values(ascending=False)
    right_totals = out.groupby(right_col, sort=False)["count"].sum().sort_values(ascending=False)
    keep_left = set(left_totals.head(max_categories).index.tolist())
    keep_right = set(right_totals.head(max_categories).index.tolist())
    out[left_col] = out[left_col].where(out[left_col].isin(keep_left), FLOW_OTHER_LABEL)
    out[right_col] = out[right_col].where(out[right_col].isin(keep_right), FLOW_OTHER_LABEL)
    out = out.groupby([left_col, right_col], as_index=False, sort=False)["count"].sum()
    return out


def compute_stack_positions(totals: pandas.Series, gap: float = 0.02) -> Dict[str, Tuple[float, float]]:
    positions: Dict[str, Tuple[float, float]] = {}
    n = len(totals)
    if n == 0:
        return positions
    total_sum = float(totals.sum())
    if total_sum <= 0:
        return positions
    top_margin = 0.04
    bottom_margin = 0.04
    usable = 1.0 - top_margin - bottom_margin - gap * max(0, n - 1)
    y_top = 1.0 - top_margin
    for category, value in totals.items():
        height = usable * (float(value) / total_sum)
        positions[str(category)] = (y_top - height, y_top)
        y_top -= height + gap
    return positions


def add_ribbon(ax, path_cls, patch_cls, x0, x1, y0_low, y0_high, y1_low, y1_high, color):
    ctrl = (x1 - x0) * 0.35
    vertices = [
        (x0, y0_high),
        (x0 + ctrl, y0_high),
        (x1 - ctrl, y1_high),
        (x1, y1_high),
        (x1, y1_low),
        (x1 - ctrl, y1_low),
        (x0 + ctrl, y0_low),
        (x0, y0_low),
        (x0, y0_high),
    ]
    codes = [
        path_cls.MOVETO,
        path_cls.CURVE4,
        path_cls.CURVE4,
        path_cls.CURVE4,
        path_cls.LINETO,
        path_cls.CURVE4,
        path_cls.CURVE4,
        path_cls.CURVE4,
        path_cls.CLOSEPOLY,
    ]
    patch = patch_cls(path_cls(vertices, codes), facecolor=color, edgecolor="none", alpha=0.55)
    ax.add_patch(patch)


def plot_taxonomy_flow(
    gene_df: pandas.DataFrame,
    out_pdf: str,
    resolver: TaxonomyResolver,
    preferred_rank: str,
    max_categories: int,
) -> None:
    plt, PdfPages, Path, PathPatch, Rectangle = get_pyplot()
    ensure_parent_dir(out_pdf)
    if gene_df.empty:
        blank_pdf(out_pdf, "HGT Taxonomy Flow", "No HGT candidate genes were found.")
        return

    plot_df = gene_df.copy()
    plot_df["recipient_label"] = plot_df.apply(
        lambda row: resolve_rank_label(row.get("gene_taxon", ""), numpy.nan, resolver, preferred_rank),
        axis=1,
    )
    plot_df["besthit_label"] = plot_df.apply(
        lambda row: resolve_rank_label(row.get("besthit_organism", ""), row.get("besthit_taxid", numpy.nan), resolver, preferred_rank),
        axis=1,
    )
    plot_df["recipient_label"] = plot_df["recipient_label"].replace("", FLOW_FALLBACK_LABEL).fillna(FLOW_FALLBACK_LABEL)
    plot_df["besthit_label"] = plot_df["besthit_label"].replace("", FLOW_FALLBACK_LABEL).fillna(FLOW_FALLBACK_LABEL)

    flow_df = (
        plot_df.groupby(["recipient_label", "besthit_label"], as_index=False, sort=False)
        .size()
        .rename(columns={"size": "count"})
    )
    if flow_df.empty:
        blank_pdf(out_pdf, "HGT Taxonomy Flow", "No taxonomy flow records were available.")
        return
    flow_df = collapse_to_top_categories(flow_df, "recipient_label", "besthit_label", max_categories=max_categories)

    left_totals = flow_df.groupby("recipient_label", sort=False)["count"].sum().sort_values(ascending=False)
    right_totals = flow_df.groupby("besthit_label", sort=False)["count"].sum().sort_values(ascending=False)
    left_order = left_totals.index.tolist()
    right_order = right_totals.index.tolist()
    left_pos = compute_stack_positions(left_totals)
    right_pos = compute_stack_positions(right_totals)

    total_flow = float(flow_df["count"].sum())
    if total_flow <= 0:
        blank_pdf(out_pdf, "HGT Taxonomy Flow", "Taxonomy flow counts were all zero.")
        return

    left_cursor = {key: value[1] for key, value in left_pos.items()}
    right_cursor = {key: value[1] for key, value in right_pos.items()}
    flow_df["recipient_label"] = pandas.Categorical(flow_df["recipient_label"], categories=left_order, ordered=True)
    flow_df["besthit_label"] = pandas.Categorical(flow_df["besthit_label"], categories=right_order, ordered=True)
    flow_df = flow_df.sort_values(["recipient_label", "besthit_label"], ascending=[True, True], kind="mergesort").reset_index(drop=True)

    plt, PdfPages, Path, PathPatch, Rectangle = get_pyplot()
    with PdfPages(out_pdf) as pdf:
        fig, ax = plt.subplots(figsize=(10, max(5.0, 0.32 * max(len(left_order), len(right_order)) + 1.5)))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        left_x0, left_x1 = 0.12, 0.18
        right_x0, right_x1 = 0.82, 0.88
        cmap = plt.get_cmap("tab20")
        color_map = {cat: cmap(idx % 20) for idx, cat in enumerate(left_order)}

        usable = 1.0 - 0.04 - 0.04 - 0.02 * max(0, len(left_order) - 1)
        unit_height = usable / total_flow
        for row in flow_df.itertuples(index=False):
            left_cat = str(row.recipient_label)
            right_cat = str(row.besthit_label)
            height = float(row.count) * unit_height
            left_high = left_cursor[left_cat]
            left_low = left_high - height
            right_high = right_cursor[right_cat]
            right_low = right_high - height
            add_ribbon(ax, Path, PathPatch, left_x1, right_x0, left_low, left_high, right_low, right_high, color_map[left_cat])
            left_cursor[left_cat] = left_low
            right_cursor[right_cat] = right_low

        for category in left_order:
            y0, y1 = left_pos[category]
            ax.add_patch(Rectangle((left_x0, y0), left_x1 - left_x0, y1 - y0, facecolor=color_map[category], edgecolor="black", linewidth=0.4))
            ax.text(left_x0 - 0.02, (y0 + y1) / 2, wrap_label(category), ha="right", va="center", fontsize=7)
            ax.text(left_x1 + 0.01, (y0 + y1) / 2, str(int(left_totals[category])), ha="left", va="center", fontsize=7)

        for category in right_order:
            y0, y1 = right_pos[category]
            ax.add_patch(Rectangle((right_x0, y0), right_x1 - right_x0, y1 - y0, facecolor="#d9d9d9", edgecolor="black", linewidth=0.4))
            ax.text(right_x1 + 0.02, (y0 + y1) / 2, wrap_label(category), ha="left", va="center", fontsize=7)
            ax.text(right_x0 - 0.01, (y0 + y1) / 2, str(int(right_totals[category])), ha="right", va="center", fontsize=7)

        ax.text((left_x0 + left_x1) / 2, 1.01, "Recipient lineage", ha="center", va="bottom", fontsize=9, fontweight="bold")
        ax.text((right_x0 + right_x1) / 2, 1.01, "Best-hit lineage", ha="center", va="bottom", fontsize=9, fontweight="bold")
        ax.text(
            0.5,
            1.04,
            f"HGT taxonomy flow ({preferred_rank})",
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
        )
        ax.text(0.5, -0.03, f"Total candidate genes: {int(total_flow)}", ha="center", va="top", fontsize=8)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    branch_df = safe_read_tsv(args.branch_tsv)
    gene_df = safe_read_tsv(args.gene_tsv)
    resolver = TaxonomyResolver(args.taxonomy_dbfile)

    plot_overview(branch_df, args.overview_pdf)
    plot_taxonomy_flow(
        gene_df=gene_df,
        out_pdf=args.taxonomy_flow_pdf,
        resolver=resolver,
        preferred_rank=args.flow_rank,
        max_categories=max(1, int(args.flow_max_categories)),
    )


if __name__ == "__main__":
    main()
