#!/usr/bin/env python3

import argparse
import math
import os
import textwrap
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Set, Tuple

import numpy
import pandas
from score_hgt_candidates import TaxonomyResolver, normalize_sci_name

OVERVIEW_NUMERIC_COLUMNS: List[Tuple[str, str, str]] = [
    ("candidate_gene_count", "Cand", "Number of genes descending from the candidate GeneRax HGT branch."),
    (
        "matched_leaf_count",
        "Tips",
        "Number of leaf rows in the database that were matched back to the candidate branch.",
    ),
    ("besthit_gene_count", "HitGenes", "Number of candidate genes with a UniProt/Swiss-Prot best hit."),
    (
        "besthit_taxid_count",
        "HitTaxID",
        "Number of candidate genes whose best hit could be resolved to a taxonomy-aware label.",
    ),
    (
        "besthit_same_superkingdom_fraction",
        "SameSK",
        "Fraction of candidate genes whose best hit falls in the same superkingdom as the focal lineage.",
    ),
    ("intron_support_fraction", "Intron", "Fraction of candidate genes with intron support recorded in stat_branch."),
    ("expression_measured_fraction", "Expr", "Fraction of candidate genes with any expression measurement."),
    (
        "clade_min_expression_pearsoncor",
        "ExprCor",
        "Minimum clade-level Pearson correlation among measured expression profiles.",
    ),
    ("synteny_support_fraction", "SynFrac", "Fraction of candidate genes with positive synteny support."),
    ("synteny_mean_support_score", "SynMean", "Mean synteny support score across candidate genes."),
    (
        "contamination_incompatible_fraction",
        "Contam",
        "Fraction of candidate genes flagged as lineage-incompatible by contamination QC.",
    ),
]

OVERVIEW_TEXT_COLUMNS: List[Tuple[str, str, int, str]] = [
    ("besthit_lca_rank_mode", "HitLCA", 18, "Most frequent best-hit lineage relationship label among candidate genes."),
    (
        "contamination_top_lca_sciname",
        "TopContam",
        28,
        "Most frequent contamination LCA scientific name among candidate genes.",
    ),
]

FLOW_FALLBACK_LABEL = "Unresolved"
FLOW_OTHER_LABEL = "Other"
TRANSFER_EDGE_COLUMNS = [
    "donor_node",
    "recipient_node",
    "hgt_event_count",
    "orthogroup_count",
    "event_fraction",
    "mapped_to_species_tree",
    "display_rank",
    "displayed",
    "phylogenetic_distance",
    "distance_metric",
    "selection_reason",
]


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Plot overview, taxonomy-flow, and transfer-tree summaries for gg_hgt outputs."
    )
    parser.add_argument("--branch_tsv", metavar="PATH", required=True, type=str)
    parser.add_argument("--gene_tsv", metavar="PATH", required=True, type=str)
    parser.add_argument("--overview_pdf", metavar="PATH", required=True, type=str)
    parser.add_argument("--taxonomy_flow_pdf", metavar="PATH", required=True, type=str)
    parser.add_argument("--taxonomy_dbfile", metavar="PATH", default="", type=str)
    parser.add_argument("--flow_rank", metavar="TEXT", default="phylum", type=str)
    parser.add_argument("--flow_max_categories", metavar="INT", default=12, type=int)
    parser.add_argument("--species_trait", default="", help="Optional species_trait TSV; numeric/binary traits appear beside species-tree tips.")
    parser.add_argument(
        "--transfer_tree_pdf",
        metavar="PATH",
        default="",
        type=str,
        help="Optional PDF path for directed HGT links over the species tree.",
    )
    parser.add_argument(
        "--transfer_edges_tsv",
        metavar="PATH",
        default="",
        type=str,
        help="Optional TSV path for all parsed donor-to-recipient edge counts.",
    )
    parser.add_argument(
        "--species_tree",
        metavar="PATH",
        default="",
        type=str,
        help="Species-tree Newick file used by the transfer-tree plot.",
    )
    parser.add_argument(
        "--transfer_tree_max_edges",
        metavar="INT",
        default=200,
        type=int,
        help="Initial mapped-direction selection limit; reverse directions are then included on shared curves. 0 selects all.",
    )
    return parser


def get_pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.patches import PathPatch, Rectangle
    from matplotlib.path import Path

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
            "",
            "## Directed Transfer Tree",
            "",
            "`hgt_transfer_tree.pdf` overlays directed donor/source-to-recipient/target links on the species tree.",
            "Optional species_trait columns show observed numeric/binary tip traits, with all text at 8 pt. The workflow automatically reads input/species_trait/species_trait.tsv; hgt_summary_species_trait accepts a path or none. Binary 1 is orange and 0 gray; numeric colors are scaled independently per column and values are printed. Missing or unmatched values show NA, never zero. Shared schema/metadata contracts apply. No ancestral states are reconstructed.",
            "Links attach to the midpoint of the horizontal branch entering each labelled node. These positions are display conventions, not estimated transfer times. Root endpoints use a dashed display-only stem; zero-length branches coincide with their nodes. Color and ranking retain endpoint-node path distance as a lineage-separation proxy, not distance between inferred transfer locations.",
            "Link width is max(0.35, 5 * count / maximum_count) points, using the maximum across all parsed pairs. The visibility floor preserves rare distant events; counts below the floor share a width. It is not a probability score.",
            "Arrowheads point to the recipient/target. `hgt_transfer_edges.tsv` contains every parseable pair, including edges not drawn in the PDF.",
            "Both directions share one curve: each arrow-end half encodes the count toward that endpoint. A one-way link has a thin source half with no source arrow. Existing reverse directions are added after initial selection (`selection_reason=reciprocal`), so displayed direction counts can exceed the limit without adding more connections. TSV rows remain directional. All internal branch names are drawn above their incoming branch midpoints, regardless of HGT participation. Species names are to the right of terminal branches. All transfer-tree text is 8 pt, including title, legend and colorbar.",
            "The PDF selects up to 200 mapped pairs by alternating event-count and distance rankings (0 selects all). Darker blue means greater tree distance; width scales with event count with the visibility floor above. Distant links are drawn last.",
            "`phylogenetic_distance` is the path length between labelled endpoint nodes, not transfer time. `distance_metric` is branch_length when every non-root branch has a finite nonnegative length and at least one is positive; otherwise the entire tree uses topology_edges. `selection_reason` records count, distance, all, or reciprocal; unselected pairs are not_displayed. Unmapped distances are missing. Rankings break ties by count, distance, and endpoint labels deterministically. The color scale uses all mapped pairs, including hidden pairs.",
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


def parse_generax_transfer(value) -> Optional[Tuple[str, str]]:
    """Return (donor/source, recipient/target) from GeneRax's Y@donor@recipient field."""
    if value is None or pandas.isna(value):
        return None
    text = str(value).strip()
    if text == "" or text.lower() in {"nan", "none", "na"}:
        return None
    parts = text.split("@", 2)
    if len(parts) != 3 or parts[0].strip().upper() != "Y":
        return None
    donor = parts[1].strip()
    recipient = parts[2].strip()
    if donor == "" or recipient == "":
        return None
    return donor, recipient


def normalize_tree_label(value) -> str:
    if value is None or pandas.isna(value):
        return ""
    text = str(value).strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in {"'", '"'}:
        text = text[1:-1].strip()
    return text


def build_species_tree_label_map(tree) -> Dict[str, str]:
    """Map exact and underscore/space aliases to unique labels in a Bio.Phylo tree."""
    label_to_canonical: Dict[str, str] = {}
    ambiguous: Set[str] = set()
    for clade in tree.find_clades(order="preorder"):
        label = normalize_tree_label(getattr(clade, "name", ""))
        if label == "":
            continue
        canonical_aliases = {label, label.replace("_", " "), label.replace(" ", "_")}
        for alias in canonical_aliases:
            if alias in ambiguous:
                continue
            previous = label_to_canonical.get(alias)
            if previous is not None and previous != label:
                label_to_canonical.pop(alias, None)
                ambiguous.add(alias)
            else:
                label_to_canonical[alias] = label
    return label_to_canonical


def load_species_tree_layout(tree_path: str):
    """Read a Newick tree and return tree, normalized x/y coordinates, and label aliases."""
    if not tree_path or not os.path.isfile(tree_path):
        return None, {}, {}, {}
    from hgt_species_tree import read_species_tree
    tree = read_species_tree(tree_path)

    terminals = tree.get_terminals()
    if len(terminals) == 0:
        return None, {}, {}, {}

    y_by_id: Dict[int, float] = {}
    for index, clade in enumerate(terminals):
        y_by_id[id(clade)] = float(len(terminals) - index - 1)

    def assign_internal_y(clade) -> float:
        if clade.is_terminal():
            return y_by_id[id(clade)]
        child_y = [assign_internal_y(child) for child in clade.clades]
        y_by_id[id(clade)] = float(sum(child_y) / len(child_y))
        return y_by_id[id(clade)]

    assign_internal_y(tree.root)

    length_x_by_id: Dict[int, float] = {}
    depth_x_by_id: Dict[int, float] = {}

    def assign_x(clade, length_x: float, depth_x: float) -> None:
        length_x_by_id[id(clade)] = length_x
        depth_x_by_id[id(clade)] = depth_x
        for child in clade.clades:
            branch_length = child.branch_length
            try:
                branch_length = float(branch_length)
            except (TypeError, ValueError):
                branch_length = 0.0
            if not math.isfinite(branch_length) or branch_length < 0:
                branch_length = 0.0
            assign_x(child, length_x + branch_length, depth_x + 1.0)

    assign_x(tree.root, 0.0, 0.0)
    max_length_x = max(length_x_by_id.values())
    source_x = length_x_by_id if max_length_x > 0 else depth_x_by_id
    max_x = max(source_x.values())
    if max_x <= 0:
        max_x = 1.0
    x_by_id = {clade_id: 0.04 + 0.46 * (value / max_x) for clade_id, value in source_x.items()}
    return tree, x_by_id, y_by_id, build_species_tree_label_map(tree)


def resolve_tree_endpoint(endpoint: str, tree_labels) -> str:
    if isinstance(tree_labels, dict):
        return str(tree_labels.get(endpoint, ""))
    if tree_labels is not None and endpoint in tree_labels:
        return endpoint
    return ""


def empty_transfer_edge_table() -> pandas.DataFrame:
    return pandas.DataFrame(columns=TRANSFER_EDGE_COLUMNS)


def build_transfer_edge_table(
    branch_df: pandas.DataFrame, tree_labels=None, max_edges: int = 200
) -> pandas.DataFrame:
    """Aggregate branch-level GeneRax HGT transfers into directed donor/recipient edges."""
    if branch_df.empty or "generax_transfer" not in branch_df.columns:
        return empty_transfer_edge_table()

    pair_counts: Counter = Counter()
    pair_orthogroups = defaultdict(set)
    total_events = 0
    orthogroup_values = branch_df["orthogroup"] if "orthogroup" in branch_df.columns else pandas.Series(
        [""] * len(branch_df), index=branch_df.index
    )
    for transfer, orthogroup in zip(
        branch_df["generax_transfer"].tolist(), orthogroup_values.tolist(), strict=True
    ):
        parsed = parse_generax_transfer(transfer)
        if parsed is None:
            continue
        donor, recipient = parsed
        pair_counts[(donor, recipient)] += 1
        total_events += 1
        if orthogroup is not None and not pandas.isna(orthogroup):
            orthogroup_text = str(orthogroup).strip()
            if orthogroup_text and orthogroup_text.lower() != "nan":
                pair_orthogroups[(donor, recipient)].add(orthogroup_text)

    if not pair_counts:
        return empty_transfer_edge_table()

    max_edges = max(0, int(max_edges))
    records = []
    mapped_rank = 0
    for donor_recipient, event_count in sorted(
        pair_counts.items(), key=lambda item: (-item[1], item[0][0], item[0][1])
    ):
        donor, recipient = donor_recipient
        donor_tree_label = resolve_tree_endpoint(donor, tree_labels)
        recipient_tree_label = resolve_tree_endpoint(recipient, tree_labels)
        mapped = int(donor_tree_label != "" and recipient_tree_label != "")
        if mapped:
            mapped_rank += 1
            displayed = int(max_edges == 0 or mapped_rank <= max_edges)
            display_rank = mapped_rank
        else:
            displayed = 0
            display_rank = 0
        records.append(
            {
                "donor_node": donor,
                "recipient_node": recipient,
                "hgt_event_count": int(event_count),
                "orthogroup_count": len(pair_orthogroups[donor_recipient]),
                "event_fraction": float(event_count) / float(total_events),
                "mapped_to_species_tree": mapped,
                "display_rank": display_rank,
                "displayed": displayed,
            }
        )
    return pandas.DataFrame.from_records(records, columns=TRANSFER_EDGE_COLUMNS)


def write_transfer_edges_tsv(edge_df: pandas.DataFrame, out_tsv: str) -> None:
    ensure_parent_dir(out_tsv)
    edge_df.to_csv(out_tsv, sep="\t", index=False, float_format="%.6f")


def blank_pdf(path: str, title: str, message: str, fontsize: int = 8) -> None:
    plt, PdfPages, _, _, _ = get_pyplot()
    ensure_parent_dir(path)
    with PdfPages(path) as pdf:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.axis("off")
        ax.text(0.5, 0.65, title, ha="center", va="center", fontsize=fontsize, fontweight="bold")
        ax.text(0.5, 0.40, message, ha="center", va="center", fontsize=fontsize)
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
        plot_df["orthogroup"].astype(str).str.strip() + ":" + plot_df["branch_id"].astype(str).str.strip()
    )
    plot_df = plot_df.sort_values(["orthogroup", "branch_id"], ascending=[True, True], kind="mergesort").reset_index(
        drop=True
    )
    numeric_cols = [col for col, _label, _desc in OVERVIEW_NUMERIC_COLUMNS if col in plot_df.columns]
    text_specs = [(col, label, width) for col, label, width, _desc in OVERVIEW_TEXT_COLUMNS if col in plot_df.columns]
    chunk_size = 60

    with PdfPages(out_pdf) as pdf:
        for start in range(0, plot_df.shape[0], chunk_size):
            chunk = plot_df.iloc[start : start + chunk_size, :].copy()
            numeric_df = chunk.loc[:, numeric_cols].apply(pandas.to_numeric, errors="coerce")
            normalized = (
                normalize_numeric_frame(numeric_df) if len(numeric_cols) > 0 else pandas.DataFrame(index=chunk.index)
            )
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
                values = chunk[col].map(lambda x, _width=width: shorten_text(x, _width)).fillna("")
                for y, value in enumerate(values.tolist()):
                    ax_text.text(0.0, y, value, va="center", ha="left", fontsize=7, clip_on=False)

            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


def draw_species_tree(ax, tree, x_by_id: Dict[int, float], y_by_id: Dict[int, float], referenced_labels: Set[str]) -> None:
    """Draw a compact rectangular species tree in the left part of an axes."""
    for parent in tree.find_clades(order="preorder"):
        children = list(parent.clades)
        if not children:
            continue
        parent_x = x_by_id[id(parent)]
        child_y = [y_by_id[id(child)] for child in children]
        ax.plot(
            [parent_x, parent_x],
            [min(child_y), max(child_y)],
            color="#777777",
            linewidth=0.65,
            solid_capstyle="round",
            zorder=3,
        )
        for child in children:
            child_x = x_by_id[id(child)]
            y_value = y_by_id[id(child)]
            ax.plot(
                [parent_x, child_x],
                [y_value, y_value],
                color="#777777",
                linewidth=0.65,
                solid_capstyle="round",
                zorder=3,
            )

    anchors = species_branch_anchors(tree, x_by_id, y_by_id)
    for clade in tree.find_clades(order="preorder"):
        label = normalize_tree_label(getattr(clade, "name", ""))
        if label == "":
            continue
        terminal = clade.is_terminal()
        ax.annotate(
            label,
            xy=(x_by_id[id(clade)], y_by_id[id(clade)]) if terminal else anchors[id(clade)],
            xytext=(4, 0) if terminal else (0, 3),
            textcoords="offset points",
            ha="left" if terminal else "center",
            va="center" if terminal else "bottom",
            fontsize=8,
            color="#555555",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 0.2},
            zorder=7,
        )


def add_transfer_distances(edge_df, tree, label_map, max_edges):
    """Measure node-to-node paths and interleave count/distance priority lists."""
    edges = edge_df.copy()
    nodes = {normalize_tree_label(c.name): c for c in tree.find_clades() if c.name}
    branches = [c for c in tree.find_clades() if c is not tree.root]
    lengths = [c.branch_length for c in branches]
    use_lengths = (all(v is not None and math.isfinite(v) and v >= 0 for v in lengths)
                   and any(v > 0 for v in lengths if v is not None))
    metric = "branch_length" if use_lengths else "topology_edges"
    paths = {id(c): tree.get_path(c) for c in tree.find_clades()}
    distances = []
    for row in edges.itertuples():
        if not row.mapped_to_species_tree:
            distances.append(float("nan"))
            continue
        a = nodes[resolve_tree_endpoint(row.donor_node, label_map)]
        b = nodes[resolve_tree_endpoint(row.recipient_node, label_map)]
        unique = {id(c): c for c in paths[id(a)]}
        for c in paths[id(b)]:
            if id(c) in unique:
                del unique[id(c)]
            else:
                unique[id(c)] = c
        distances.append(sum(c.branch_length if use_lengths else 1 for c in unique.values()))
    edges["phylogenetic_distance"] = distances
    edges["distance_metric"] = metric
    edges["selection_reason"] = "not_displayed"
    edges["displayed"] = 0
    edges["display_rank"] = 0
    mapped = edges.loc[edges.mapped_to_species_tree.eq(1)]
    count_order = mapped.sort_values(
        ["hgt_event_count", "phylogenetic_distance", "donor_node", "recipient_node"],
        ascending=[False, False, True, True]).index.tolist()
    distance_order = mapped.sort_values(
        ["phylogenetic_distance", "hgt_event_count", "donor_node", "recipient_node"],
        ascending=[False, False, True, True]).index.tolist()
    limit = len(mapped) if max_edges <= 0 else min(max_edges, len(mapped))
    chosen = set()
    for count_idx, distance_idx in zip(count_order, distance_order, strict=True):
        for idx, reason in ((count_idx, "count"), (distance_idx, "distance")):
            if idx in chosen or len(chosen) >= limit:
                continue
            chosen.add(idx)
            edges.loc[idx, ["displayed", "display_rank", "selection_reason"]] = [
                1, len(chosen), "all" if max_edges <= 0 else reason]
    # Complete the reverse direction of every selected connection so that a
    # bidirectional link never looks unidirectional because of the display cap.
    pair_indices = {(r.donor_node, r.recipient_node): r.Index for r in mapped.itertuples()}
    for idx in sorted(chosen):
        row = edges.loc[idx]
        reverse = pair_indices.get((row.recipient_node, row.donor_node))
        if reverse is not None and not edges.at[reverse, "displayed"]:
            edges.loc[reverse, ["displayed", "display_rank", "selection_reason"]] = [
                1, int(edges.display_rank.max()) + 1, "reciprocal"]
    return edges


def transfer_connections(display_df):
    """Group displayed directions by their unordered endpoint pair."""
    groups = {}
    for row in display_df.itertuples(index=False):
        key = tuple(sorted((str(row.donor_node), str(row.recipient_node))))
        groups.setdefault(key, []).append(row)
    return sorted(groups.items(), key=lambda item: (
        max(r.phylogenetic_distance for r in item[1]), item[0]))


def transfer_half_paths(ax, start, end):
    """Split one shallow quadratic in display space, returning midpoint-to-tip paths."""
    from matplotlib.path import Path

    a, b = ax.transData.transform([start, end])
    inverse = ax.transData.inverted()
    if numpy.allclose(a, b, rtol=0, atol=1e-8):
        # Same-branch transfers (or coincident zero-length branch anchors) must
        # remain visible. The loop is a display convention, not a time estimate.
        radius = 12 * ax.figure.dpi / 72
        middle = a + [2 * radius, 0]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        return (Path(inverse.transform([middle, middle + [0, radius], a + [0, radius], a]), codes),
                Path(inverse.transform([middle, middle - [0, radius], b - [0, radius], b]), codes))
    delta = b - a
    control = (a + b) / 2 + 0.10 * numpy.array([delta[1], -delta[0]])
    left, right = (a + control) / 2, (control + b) / 2
    middle = (left + right) / 2
    codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
    return (Path(inverse.transform([middle, left, a]), codes),
            Path(inverse.transform([middle, right, b]), codes))


def species_branch_anchors(tree, x_by_id, y_by_id):
    """Midpoints of incoming horizontal branches; root uses a display-only stem."""
    anchors = {id(tree.root): (x_by_id[id(tree.root)] / 2, y_by_id[id(tree.root)])}
    for parent in tree.find_clades():
        for child in parent.clades:
            anchors[id(child)] = (
                (x_by_id[id(parent)] + x_by_id[id(child)]) / 2,
                y_by_id[id(child)],
            )
    return anchors


def read_transfer_traits(path):
    """Use the shared trait contract; preserve unknown values and reject ambiguous IDs."""
    if not path:
        return pandas.DataFrame()
    from species_trait_contract import select_analysis_traits

    frame, _ = select_analysis_traits(path)
    identifiers = frame.iloc[:, 0].map(lambda s: normalize_tree_label(s).replace(" ", "_"))
    if identifiers.eq("").any() or identifiers.duplicated().any():
        raise ValueError("Species traits require unique nonempty species IDs (including space/underscore aliases)")
    values = frame.iloc[:, 1:].replace({"": numpy.nan, "NA": numpy.nan, "NaN": numpy.nan, "nan": numpy.nan})
    values = values.apply(pandas.to_numeric, errors="raise")
    if numpy.isinf(values.to_numpy(dtype=float)).any():
        raise ValueError("Species traits must be finite or missing")
    values.index = identifiers
    return values


def draw_transfer_traits(ax, tree, y_by_id, traits):
    """Draw observed tip traits only; no ancestral reconstruction or missing-to-zero conversion."""
    from matplotlib import colormaps
    from matplotlib.colors import Normalize
    from matplotlib.patches import Rectangle

    for column_index, column in enumerate(traits.columns):
        series = traits[column]
        observed = series.dropna()
        binary = set(observed.unique()) <= {0, 1}
        low, high = (0, 1) if binary or observed.empty else (float(observed.min()), float(observed.max()))
        scale = Normalize(low, high if high > low else low + 1)
        x = 0.81 + column_index * 0.11
        ax.text(x + 0.025, len(tree.get_terminals()) + 0.2, column,
                rotation=45, ha="left", va="bottom", fontsize=8)
        for tip in tree.get_terminals():
            key = normalize_tree_label(tip.name).replace(" ", "_")
            value = series.get(key, numpy.nan)
            missing = pandas.isna(value)
            color = "#ffffff" if missing else ("#e69f00" if value == 1 else "#eeeeee") if binary else colormaps["YlOrBr"](scale(value))
            ax.add_patch(Rectangle((x, y_by_id[id(tip)] - 0.4), 0.05, 0.8,
                                   facecolor=color, edgecolor="#aaaaaa", linewidth=0.3, zorder=5))
            ax.text(x + 0.025, y_by_id[id(tip)], "NA" if missing else f"{value:g}",
                    ha="center", va="center", fontsize=8, zorder=6)
    ax.text(0.81, -0.9, "Tip traits; NA = missing", fontsize=8, ha="left", va="top")


def plot_transfer_tree(
    branch_df: pandas.DataFrame,
    out_pdf: str,
    species_tree_path: str = "",
    edges_tsv: str = "",
    max_edges: int = 200,
    species_trait_path: str = "",
) -> None:
    """Plot directed GeneRax HGT event counts over a species tree."""
    tree, x_by_id, y_by_id, tree_label_map = load_species_tree_layout(species_tree_path)
    traits = read_transfer_traits(species_trait_path)
    edge_df = build_transfer_edge_table(branch_df, tree_labels=tree_label_map, max_edges=max_edges)
    if tree is not None and not edge_df.empty:
        edge_df = add_transfer_distances(edge_df, tree, tree_label_map, max_edges)
    if edges_tsv:
        write_transfer_edges_tsv(edge_df, edges_tsv)
    if not out_pdf:
        return

    if tree is None:
        message = "Species tree was not found or could not be parsed."
        if species_tree_path:
            message += f"\nRequested path: {species_tree_path}"
        blank_pdf(out_pdf, "HGT Transfer Tree", message)
        return
    if edge_df.empty:
        blank_pdf(out_pdf, "HGT Transfer Tree", "No parseable Y@donor@recipient transfers were found.")
        return

    display_df = edge_df.loc[edge_df["displayed"].astype(int).eq(1)].copy()
    if display_df.empty:
        blank_pdf(
            out_pdf,
            "HGT Transfer Tree",
            "Transfer records were found, but no donor/recipient pair matched the species-tree labels.",
        )
        return

    plt, PdfPages, _, _, _ = get_pyplot()
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyArrowPatch

    terminal_count = len(tree.get_terminals())
    branch_anchors = species_branch_anchors(tree, x_by_id, y_by_id)
    fig_height = max(8.0, min(28.0, 0.16 * terminal_count + 2.8))
    edge_color = "#2b6ca3"
    cmap = LinearSegmentedColormap.from_list("hgt_distance", ["#c6dbef", "#08306b"])
    max_distance = float(edge_df["phylogenetic_distance"].max())
    norm = Normalize(0, max_distance if max_distance > 0 else 1)
    referenced_labels = set()
    for row in display_df.itertuples(index=False):
        donor_label = resolve_tree_endpoint(str(row.donor_node), tree_label_map)
        recipient_label = resolve_tree_endpoint(str(row.recipient_node), tree_label_map)
        if donor_label:
            referenced_labels.add(donor_label)
        if recipient_label:
            referenced_labels.add(recipient_label)
    clade_by_label = {}
    for clade in tree.find_clades(order="preorder"):
        label = normalize_tree_label(getattr(clade, "name", ""))
        if label:
            clade_by_label[label] = clade

    ensure_parent_dir(out_pdf)
    with PdfPages(out_pdf) as pdf:
        fig, ax = plt.subplots(figsize=(14.0, fig_height))
        ax.set_xlim(0.0, max(1.04, 0.92 + 0.11 * len(traits.columns)))
        ax.set_ylim(-1.2, float(terminal_count) + 2.0)
        ax.axis("off")

        display_df = display_df.sort_values(
            ["phylogenetic_distance", "hgt_event_count", "donor_node", "recipient_node"],
            ascending=[True, True, True, True],
            kind="mergesort",
        )
        max_count = max(1, int(edge_df["hgt_event_count"].max()))
        connections = transfer_connections(display_df)
        for (a, b), rows in connections:
            a_point = branch_anchors[id(clade_by_label[resolve_tree_endpoint(a, tree_label_map)])]
            b_point = branch_anchors[id(clade_by_label[resolve_tree_endpoint(b, tree_label_map)])]
            halves = transfer_half_paths(ax, a_point, b_point)
            by_target = {str(r.recipient_node): r for r in rows}
            for half_index, (target, path) in enumerate(zip((a, b), halves, strict=True)):
                row = None if a == b and half_index == 0 else by_target.get(target)
                # One-way connections keep a thin source half without an arrow.
                width = max(0.35, 5.0 * int(row.hgt_event_count) / max_count) if row else 0.35
                ax.add_patch(FancyArrowPatch(
                    path=path, arrowstyle="-|>" if row else "-",
                    mutation_scale=7.0, linewidth=width,
                    color=cmap(norm(rows[0].phylogenetic_distance)),
                    capstyle="butt", zorder=2,
                ))

        draw_species_tree(ax, tree, x_by_id, y_by_id, referenced_labels)
        if len(traits.columns):
            draw_transfer_traits(ax, tree, y_by_id, traits)
        if normalize_tree_label(tree.root.name):
            ax.plot([0, x_by_id[id(tree.root)]], [y_by_id[id(tree.root)]] * 2,
                    color="#777777", linewidth=0.65, linestyle="--")

        mapped_events = int(edge_df.loc[edge_df["mapped_to_species_tree"].astype(int).eq(1), "hgt_event_count"].sum())
        displayed_events = int(display_df["hgt_event_count"].sum())
        total_events = int(edge_df["hgt_event_count"].sum())
        handles = [
            Line2D([0], [0], color="#777777", linewidth=0.65, label="Species-tree branch"),
        ]
        for count in sorted({1, max(1, max_count // 10), max_count}):
            handles.append(Line2D([0], [0], color=edge_color,
                                  linewidth=max(0.35, 5.0 * count / max_count), label=f"{count} HGT events"))
        color_ax = ax.inset_axes([0.63, 0.12, 0.025, 0.40])
        colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=color_ax)
        colorbar.set_label(f"Endpoint-node distance ({edge_df['distance_metric'].iloc[0]})")
        ax.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(0.01, 1.01),
            frameon=False,
            fontsize=8,
            ncol=2,
            handlelength=2.5,
            columnspacing=1.0,
        )
        ax.text(
            0.5,
            float(terminal_count) + 1.8,
            "HGT events mapped on the species tree",
            ha="center",
            va="bottom",
            fontsize=8,
            fontweight="bold",
        )
        ax.text(
            0.5,
            float(terminal_count) + 0.88,
            "Arrow-end half width: directional count (0.35 pt floor) | darker: greater distance | shared curve for both directions",
            ha="center",
            va="bottom",
            fontsize=8,
            color="#444444",
        )
        ax.text(
            0.5,
            -0.04,
            f"Parsed events: {total_events:,} | mapped events: {mapped_events:,} | displayed events: {displayed_events:,} | directions: {len(display_df):,} | connections: {len(connections):,}",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=8,
            color="#444444",
        )
        from matplotlib.text import Text
        for text in fig.findobj(match=Text):
            text.set_fontsize(8)
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


def collapse_to_top_categories(
    count_df: pandas.DataFrame, left_col: str, right_col: str, max_categories: int
) -> pandas.DataFrame:
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
    preferred_rank_normalized = str(preferred_rank).strip().lower()
    precomputed_phylum_columns = {"recipient_phylum", "donor_phylum"}.issubset(plot_df.columns)
    precomputed_phylum_values = pandas.DataFrame(index=plot_df.index)
    if precomputed_phylum_columns:
        for column in ["recipient_phylum", "donor_phylum"]:
            precomputed_phylum_values[column] = plot_df[column].map(
                lambda value: "" if pandas.isna(value) else str(value).strip()
            )
    use_precomputed_phylum = (
        preferred_rank_normalized == "phylum"
        and precomputed_phylum_columns
        and (
            precomputed_phylum_values["recipient_phylum"].ne("").any()
            or precomputed_phylum_values["donor_phylum"].ne("").any()
        )
    )
    if use_precomputed_phylum:
        plot_df["recipient_label"] = precomputed_phylum_values["recipient_phylum"]
        plot_df["besthit_label"] = precomputed_phylum_values["donor_phylum"]
    else:
        plot_df["recipient_label"] = plot_df.apply(
            lambda row: resolve_rank_label(row.get("gene_taxon", ""), numpy.nan, resolver, preferred_rank),
            axis=1,
        )
        plot_df["besthit_label"] = plot_df.apply(
            lambda row: resolve_rank_label(
                row.get("besthit_organism", ""), row.get("besthit_taxid", numpy.nan), resolver, preferred_rank
            ),
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
    flow_df = flow_df.sort_values(
        ["recipient_label", "besthit_label"], ascending=[True, True], kind="mergesort"
    ).reset_index(drop=True)

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
            add_ribbon(
                ax, Path, PathPatch, left_x1, right_x0, left_low, left_high, right_low, right_high, color_map[left_cat]
            )
            left_cursor[left_cat] = left_low
            right_cursor[right_cat] = right_low

        for category in left_order:
            y0, y1 = left_pos[category]
            ax.add_patch(
                Rectangle(
                    (left_x0, y0),
                    left_x1 - left_x0,
                    y1 - y0,
                    facecolor=color_map[category],
                    edgecolor="black",
                    linewidth=0.4,
                )
            )
            ax.text(left_x0 - 0.02, (y0 + y1) / 2, wrap_label(category), ha="right", va="center", fontsize=7)
            ax.text(left_x1 + 0.01, (y0 + y1) / 2, str(int(left_totals[category])), ha="left", va="center", fontsize=7)

        for category in right_order:
            y0, y1 = right_pos[category]
            ax.add_patch(
                Rectangle(
                    (right_x0, y0), right_x1 - right_x0, y1 - y0, facecolor="#d9d9d9", edgecolor="black", linewidth=0.4
                )
            )
            ax.text(right_x1 + 0.02, (y0 + y1) / 2, wrap_label(category), ha="left", va="center", fontsize=7)
            ax.text(
                right_x0 - 0.01, (y0 + y1) / 2, str(int(right_totals[category])), ha="right", va="center", fontsize=7
            )

        ax.text(
            (left_x0 + left_x1) / 2, 1.01, "Recipient lineage", ha="center", va="bottom", fontsize=9, fontweight="bold"
        )
        ax.text(
            (right_x0 + right_x1) / 2, 1.01, "Best-hit lineage", ha="center", va="bottom", fontsize=9, fontweight="bold"
        )
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
    if args.transfer_tree_pdf or args.transfer_edges_tsv:
        if not args.transfer_tree_pdf:
            transfer_pdf = os.path.join(os.path.dirname(args.transfer_edges_tsv), "hgt_transfer_tree.pdf")
        else:
            transfer_pdf = args.transfer_tree_pdf
        plot_transfer_tree(
            branch_df=branch_df,
            out_pdf=transfer_pdf,
            species_tree_path=args.species_tree,
            edges_tsv=args.transfer_edges_tsv,
            max_edges=max(0, int(args.transfer_tree_max_edges)),
            species_trait_path=args.species_trait,
        )


if __name__ == "__main__":
    main()
