"""Gene-model summary from explicit, reconciled reporting counts.

Rescue counts distinct loci, refinement counts coding paths, and UTR percentages
use final exported gene models. Plotting never runs predictors or BUSCO.
"""

import csv
import math
from dataclasses import dataclass
from functools import cache
from io import StringIO
from pathlib import Path

METHODS = ("synteny_based", "miniprot", "oriongeno")
METHOD_COLOURS = ("#9ecae1", "#4292c6", "#238b45")
MINIPROT_STACK_COLOURS = ("#9ecae1", "#4292c6", "#08519c")
SWISS_KEYS = ("te_only", "other_only", "both", "no_informative_hit", "not_assessed")
SWISS_COLOURS = ("#b2182b", "#e3a64d", "#6b3529", "#d4d4d4", "#f4f4f4")
SWISS_LABELS = ("TE-related only", "Other proteins only", "TE-related + other", "No informative hit", "Not assessed")
UTR_KEYS = ("five_only", "three_only", "both", "none")
UTR_COLOURS = ("#f1b67a", "#e17c3a", "#9c3d10", "#e5e5e5")
ORIGIN_KEYS = ("original_gff_only", "oriongeno_only", "both")
ORIGIN_COLOURS = ("#a6cbe5", "#4b91c4", "#154f79")
PATH_CHANGE_COLOURS = ("#a35400", "#7d320e")
PATH_METHOD_LABELS = ("Miniprot only", "OrionGeno only", "Both methods")
PHASES = ("pre_rescue", "before", "after")
OMARK_STAGES = ("input", "rescued", "refined")
OMARK_HOG_COLOURS = ("#395967", "#dca256")
OMARK_PROTEIN_COLOURS = ("#385d8e", "#dfe4ea")
ARCHITECTURE_COLOURS = ("#287262", "#794c91")
ROW_PITCH_IN = 0.35
BAR_HEIGHT_PT = 4.4
BAR_STEP_PT = 8.0
BAR_HEIGHT = BAR_HEIGHT_PT / (72 * ROW_PITCH_IN)


def bar_offsets(count):
    if count not in (1, 2, 3):
        raise ValueError("The figure uses one, two or three bars per species")
    step = BAR_STEP_PT / (72 * ROW_PITCH_IN)
    return tuple((index - (count - 1) / 2) * step for index in range(count))


def layout(rows, omark=False, labelled=False, show_gene_totals=False):
    if type(rows) is not int or rows < 1:
        raise ValueError("The figure requires at least one species")
    width, left, right, gap = 7.2, 1.13, 0.14, 0.12
    panel = (width - left - right - 4 * gap) / 5
    columns = 6 if omark else 4
    positions = [(left + i * (panel + gap), panel) for i in range(columns)]
    extra = 0.72 if omark and show_gene_totals else 0
    last_x, last_width = positions[-1]
    positions[-1] = (last_x, last_width + extra)
    width = last_x + last_width + extra + right
    return dict(
        width=width,
        height=ROW_PITCH_IN * rows + 4.80 + (0.25 if labelled else 0),
        bottom=4.05,
        plot_height=ROW_PITCH_IN * rows,
        positions=positions,
    )


def validate(data):
    rows = data["species"]
    groups = data["groups"]
    if (
        data.get("schema_version", 1) != 1
        or data["mode"] not in {"proteins", "transcriptome"}
        or type(groups) is not int
        or groups < 1
        or not isinstance(rows, list)
        or not rows
        or any(
            not isinstance(r.get("species"), str) or not r["species"] or type(r.get("analysed")) is not bool
            for r in rows
        )
        or len({r["species"] for r in rows}) != len(rows)
    ):
        raise ValueError("The figure requires unique species and a common positive BUSCO total")
    if data["coding_path_legend"] != {"labels": list(PATH_METHOD_LABELS), "colours": list(ORIGIN_COLOURS)}:
        raise ValueError("Refined path methods must be the three exclusive nomination groups")
    dataset = data.get("busco_dataset")
    if dataset is not None and (not isinstance(dataset, str) or not dataset.strip() or "\n" in dataset):
        raise ValueError("BUSCO dataset must be a nonempty single-line name")
    # Reporting metadata must never silently replace missing thresholds with defaults.
    method_notes(data)
    for row in rows:
        required_phases = PHASES if row["analysed"] else PHASES[:1]
        # Historical CDS-only reports may include unplotted placeholders.
        phases = required_phases if not row.get("busco_refined_assessed") else (*required_phases, "after")
        for phase in phases:
            result = row["busco"][phase]
            if (
                any(
                    type(result[k]) is not int or result[k] < 0
                    for k in ["single", "duplicated", "fragmented", "missing"]
                )
                or sum(result[k] for k in ["single", "duplicated", "fragmented", "missing"]) != groups
                or result.get("total", groups) != groups
                or result.get("complete", result["single"] + result["duplicated"])
                != result["single"] + result["duplicated"]
            ):
                raise ValueError("Invalid BUSCO group counts")
        if not row["analysed"]:
            if any(row.get(k) is not None for k in ["rescue", "coding_paths", "utr"]):
                raise ValueError("Unavailable model panels cannot contain invented counts")
            continue
        utr = row["utr"]
        n = utr["final_genes"]
        if (
            type(n) is not int
            or n <= 0
            or sum(utr["sides"][k] for k in UTR_KEYS) != n
            or sum(utr["origins"][k] for k in (*ORIGIN_KEYS, "none")) != n
            or utr["sides"]["none"] != utr["origins"]["none"]
        ):
            raise ValueError("UTR bars must share the final gene denominator")
        for counts in [utr["sides"], utr["origins"]]:
            if any(type(n) is not int or n < 0 for n in counts.values()):
                raise ValueError("Invalid UTR gene counts")
        paths = row["coding_paths"]
        if (
            any(type(paths[k]) is not int or paths[k] < 0 for k in ("repairs", "isoforms"))
            or ("counts" in paths) != ("colours" in paths)
            or ("counts" in paths and (
                paths["colours"] != data["coding_path_legend"]["colours"]
                or len(paths["counts"]) != len(paths["colours"])
                or any(type(n) is not int or n < 0 for n in paths["counts"])
                or sum(paths["counts"]) != paths["repairs"] + paths["isoforms"]
            ))
            or ("architecture" in paths and (
                set(paths["architecture"]) != {"merges", "splits"}
                or any(type(n) is not int or n < 0 for n in paths["architecture"].values())
            ))
        ):
            raise ValueError("Coding-path categories differ")
        rescue = row["rescue"]
        support = rescue.get("swissprot_support")
        if support is not None and (
            not isinstance(support, dict)
            or set(support) != set(SWISS_KEYS)
            or any(type(n) is not int or n < 0 for n in support.values())
            or sum(support.values()) != rescue["unique_loci"]
        ):
            raise ValueError("Swiss-Prot support must partition the accepted rescue loci")
        if any(
            type(rescue[key]) is not int or rescue[key] < 0
            for key in ("unique_loci", "loci_supported_by_multiple_methods")
        ):
            raise ValueError("Rescue locus totals must be nonnegative integers")
        for membership, count in rescue["exclusive_intersections"].items():
            names = membership.split("+")
            if (
                not names
                or len(set(names)) != len(names)
                or not set(names) <= set(METHODS)
                or type(count) is not int
                or count < 0
            ):
                raise ValueError("Invalid exclusive rescue locus intersection")
        for method, amount in rescue["method_counts"].items():
            if (
                type(amount) is not int
                or amount < 0
                or amount
                != sum(count for key, count in rescue["exclusive_intersections"].items() if method in key.split("+"))
            ):
                raise ValueError("Rescue methods disagree with actual locus intersections")
        if rescue["loci_supported_by_multiple_methods"] != sum(
            count for key, count in rescue["exclusive_intersections"].items() if "+" in key
        ):
            raise ValueError("Rescue overlap count differs")
        if (
            set(rescue["method_counts"]) != set(METHODS)
            or sum(rescue["exclusive_intersections"].values()) != rescue["unique_loci"]
        ):
            raise ValueError("Missing complete rescue origin/intersection evidence")
    return rows


def miniprot_partition(rescue):
    # Only is relative to the other Miniprot route; OrionGeno may also support it.
    both = sum(
        count
        for membership, count in rescue["exclusive_intersections"].items()
        if {"synteny_based", "miniprot"} <= set(membership.split("+"))
    )
    guided = rescue["method_counts"]["synteny_based"] - both
    independent = rescue["method_counts"]["miniprot"] - both
    if min(guided, independent, both) < 0:
        raise ValueError("Miniprot route intersection exceeds its route total")
    return guided, independent, both


def method_notes(data):
    """Describe only recorded reporting thresholds; absence is not a default."""
    rescue = data.get("rescue_parameters", {})
    evidence = data.get("rescue_swissprot_evidence", {})
    if not isinstance(evidence, dict):
        raise ValueError("Swiss-Prot reporting evidence must be a dictionary")
    swiss = evidence.get("parameters", {})
    if not isinstance(rescue, dict) or not isinstance(swiss, dict):
        raise ValueError("Reporting thresholds must be dictionaries")
    fractions = ("minimum_coverage", "minimum_identity", "minimum_broad_coverage",
                 "minimum_lower_coverage_identity", "cscore")
    for key in fractions:
        if key in rescue and (type(rescue[key]) not in (int, float) or not math.isfinite(rescue[key])
                              or not 0 <= rescue[key] <= 1):
            raise ValueError("Invalid Miniprot reporting threshold: " + key)
    for key in ("min_anchors", "max_interval"):
        if key in rescue and (type(rescue[key]) is not int or rescue[key] <= 0):
            raise ValueError("Invalid synteny reporting threshold: " + key)
    for key in ("evalue", "minimum_alignment", "query_coverage", "target_coverage", "score_fraction"):
        if key in swiss and (type(swiss[key]) not in (int, float) or not math.isfinite(swiss[key]) or swiss[key] <= 0
                            or key in {"query_coverage", "target_coverage", "score_fraction"} and swiss[key] > 1):
            raise ValueError("Invalid Swiss-Prot reporting threshold: " + key)
    notes = []
    homology = fractions[:4]
    if all(key in rescue for key in homology):
        coverage, identity, broad, lower = (rescue[k] for k in homology)
        lower = max(identity, lower)
        if broad > coverage and lower != identity:
            rule = (f"identity >= {100 * lower:g}% at {100 * coverage:g}-<{100 * broad:g}% coverage, "
                    f"or >= {100 * identity:g}% at >= {100 * broad:g}% coverage.")
        else:
            rule = f"identity >= {100 * identity:g}%."
        notes.append(f"Miniprot (both routes): protein coverage >= {100 * coverage:g}%; " + rule)
    else:
        notes.append("Miniprot homology thresholds were not recorded in this reporting input.")
    if all(key in rescue for key in ("min_anchors", "cscore", "max_interval")):
        notes.append(f'Synteny guidance: blocks >= {rescue["min_anchors"]} anchors; C-score >= {rescue["cscore"]:g}; '
                     f'candidate interval <= {rescue["max_interval"] / 1000:g} kb.')
    notes.append('Miniprot "only" excludes the other Miniprot route.')
    if all(key in swiss for key in ("evalue", "minimum_alignment", "query_coverage", "target_coverage", "score_fraction")):
        notes.append(f'Swiss-Prot: E <= {swiss["evalue"]:g}; paired residues >= {swiss["minimum_alignment"]:g}; '
                     f'query/target coverage >= {100 * swiss["query_coverage"]:g}%/{100 * swiss["target_coverage"]:g}%; '
                     f'bit score >= {100 * swiss["score_fraction"]:g}% of best qualifying hit.')
    else:
        notes.append("Swiss-Prot thresholds were not recorded; loci without supplied support counts are not assessed.")
    notes.extend([
        'TE-related/other only: >=1 informative hit in that class and 0 in the other. TE flags: explicit TE protein name,',
        'Transposable element keyword or transposase-activity GO term.',
    ])
    return notes


@dataclass(frozen=True)
class OMArkStage:
    proteins: int
    consistent_full: int
    conserved_hogs: int
    missing_hogs: int

    def counts(self):
        return (
            (self.conserved_hogs - self.missing_hogs, self.conserved_hogs),
            (self.consistent_full, self.proteins),
        )


def omark_comparison_counts(tsv):
    """Read within-stage numerators and denominators from the OMArk comparison."""
    reader = csv.DictReader(StringIO(tsv), delimiter="\t")
    required = {
        "species", "stage", "representative_proteins", "consistent_full", "conserved_hogs", "missing_hogs",
    }
    if not required <= set(reader.fieldnames or ()):
        raise ValueError("Expected a native OMArk comparison TSV")
    stages = {}
    for row in reader:
        species, stage = row["species"], row["stage"]
        if not species or stage not in OMARK_STAGES:
            raise ValueError("Invalid species or stage in OMArk comparison TSV")
        counts = stages.setdefault(species, {})
        if stage in counts:
            raise ValueError(f"Duplicate OMArk comparison row: {species}, {stage}")
        if bool(row["conserved_hogs"]) != bool(row["missing_hogs"]):
            raise ValueError(f"Incomplete OMArk HOG counts: {species}, {stage}")
        value = OMArkStage(
            proteins=int(row["representative_proteins"]),
            consistent_full=int(row["consistent_full"]),
            conserved_hogs=int(row["conserved_hogs"] or 0),
            missing_hogs=int(row["missing_hogs"] or 0),
        )
        if not (0 <= value.consistent_full <= value.proteins and 0 <= value.missing_hogs <= value.conserved_hogs):
            raise ValueError(f"Invalid OMArk counts: {species}, {stage}")
        counts[stage] = value.counts()
    return stages


def _count_axis(ax, maximum, bar_fraction=0.66):
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    ax.set_xlim(0, max(1, maximum) / bar_fraction)
    ticks = MaxNLocator(nbins=3, integer=True, steps=[1, 2, 2.5, 5, 10]).tick_values(0, max(1, maximum))
    ax.set_xticks([tick for tick in ticks if 0 <= tick <= maximum])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x / 1000:g}k" if x >= 1000 else f"{x:g}"))
    ax.spines["bottom"].set_bounds(0, max(1, maximum))


def _annotate_percent(ax, y, count, total, gene_total=None):
    label = f"{100 * count / total:.1f}%"
    if gene_total is not None:
        label += f" (n={gene_total:,})"
    ax.text(0.99, y, label, transform=ax.get_yaxis_transform(), ha="right", va="center")


def _bar(ax, y, width, colour, left=0):
    from matplotlib.patches import Rectangle

    # Both limits are already fixed, so these segments need no autoscaling.
    ax.add_artist(Rectangle((left, y - BAR_HEIGHT / 2), width, BAR_HEIGHT,
                            transform=ax.transData, facecolor=colour, edgecolor="none"))


def _plot_omark_counts(ax, names, counts, metric, colours, gene_totals):
    maximum = max(
        (stages[metric][1] for species in names for stages in counts.get(species, {}).values()),
        default=0,
    )
    _count_axis(ax, maximum, 0.40 if gene_totals is not None else 0.66)
    for i, species in enumerate(names):
        stages = counts.get(species, {})
        for phase, offset in zip(OMARK_STAGES, bar_offsets(3), strict=True):
            if phase not in stages:
                continue
            count, total = stages[phase][metric]
            if not total:
                continue
            _bar(ax, i + offset, count, colours[0])
            _bar(ax, i + offset, total - count, colours[1], left=count)
            gene_total = None if gene_totals is None else gene_totals.get((species, phase))
            _annotate_percent(ax, i + offset, count, total, gene_total)


def _place_legends(fig, geometry, data, *, omark=False, show_method_notes=True, _expanded=False):
    """Fit each bar-group legend using both PNG and PDF font metrics."""
    import io

    from matplotlib.backends.backend_pdf import PdfFile, RendererPdf
    from matplotlib.font_manager import FontProperties
    from matplotlib.patches import Patch

    font = FontProperties(family="sans-serif", size=7.7)
    bold = FontProperties(family="sans-serif", size=7.7, weight="bold")
    legend_start, text_start = len(fig.legends), len(fig.texts)
    # At 72 dpi, the figure transforms and PDF renderer both use points. Agg
    # alone underestimates long legend stacks when the PDF uses core fonts.
    pdf = PdfFile(io.BytesIO())
    renderer = RendererPdf(pdf, 72, geometry["height"], geometry["width"])
    renderers = (renderer, fig.canvas.get_renderer())

    @cache
    def width(text, prop):
        return max(r.get_text_width_height_descent(text, prop, False)[0] for r in renderers)

    def wrap(text, points, prop):
        lines = []
        for paragraph in text.split("\n"):
            line = ""
            for word in paragraph.replace("-", "- ").split():
                candidate = (line + ("" if line.endswith("-") else " ") + word).strip()
                if line and width(candidate, prop) > points:
                    lines.append(line)
                    line = word
                else:
                    line = candidate
                # Long lineage identifiers and arbitrary labels must also fit.
                while width(line, prop) > points:
                    end = len(line) - 1
                    while end > 1 and width(line[:end], prop) > points:
                        end -= 1
                    lines.append(line[:end])
                    line = line[end:]
            lines.append(line)
        return "\n".join(lines)

    busco = list(zip(("#000000", "#B22222", "#666666", "#CCCCCC"),
                     ("Single-copy", "Duplicated", "Fragmented", "Missing"), strict=True))
    groups = [
        [("Miniprot", list(zip(MINIPROT_STACK_COLOURS,
                              ("Synteny-guided only", "Synteny-independent only", "Both routes"), strict=True))),
         ("OrionGeno", [(METHOD_COLOURS[2], "Supported loci")]),
         ("Swiss-Prot annotation of the Miniprot/OrionGeno union", list(zip(SWISS_COLOURS, SWISS_LABELS, strict=True)))],
        [("Model changes", list(zip(PATH_CHANGE_COLOURS, ("Model revision", "Isoform addition"), strict=True)))],
        [("UTR ends", list(zip(UTR_COLOURS, ("5' only", "3' only", "Both ends", "No UTR"), strict=True))),
         ("UTR sources", list(zip((*ORIGIN_COLOURS, UTR_COLOURS[-1]),
                                  ("Existing/selected only" if data.get("utr_include_selected") else "Original GFF only",
                                   "OrionGeno only", "Both sources", "No UTR"), strict=True)))],
        [("BUSCO categories", busco)],
    ]
    path_rows = [r["coding_paths"] for r in data["species"] if r["analysed"]]
    if any("architecture" in paths for paths in path_rows):
        groups[1][0] = (
            "CDS repairs", list(zip(PATH_CHANGE_COLOURS, ("Representative", "Other isoform"), strict=True))
        )
        groups[1].append(("Model structure", list(zip(ARCHITECTURE_COLOURS, ("Merged", "Split"), strict=True))))
    if omark:
        groups.extend([
            [("Conserved HOGs", list(zip(OMARK_HOG_COLOURS, ("Detected", "Missing"), strict=True)))],
            [("Representative proteins", list(zip(OMARK_PROTEIN_COLOURS, ("Consistent and full-length", "Remaining"), strict=True)))],
        ])
    bottom = geometry["height"]
    try:
        for column, ((x, size), blocks) in enumerate(zip(geometry["positions"], groups, strict=True)):
            top = geometry["bottom"] - 0.50
            for title, items in blocks:
                legend = fig.legend(
                    handles=[Patch(facecolor=colour, label=wrap(label, size * 72 - 12, font))
                             for colour, label in items],
                    title=wrap(title, size * 72 - 1.5, bold),
                    loc="upper left", bbox_to_anchor=(x / geometry["width"], top / geometry["height"]),
                    frameon=False, prop=font, title_fontproperties=bold, handlelength=.7,
                    handletextpad=.4, borderpad=0, borderaxespad=0, labelspacing=.30,
                    alignment="left",
                )
                bounds = [legend.get_window_extent(r).transformed(fig.dpi_scale_trans.inverted()) for r in renderers]
                if any(b.x0 < x - .001 or b.x1 > x + size + .001 for b in bounds):
                    raise ValueError("Legend exceeds its column: " + title)
                bottom = min(bottom, *(b.y0 for b in bounds))
                top = min(b.y0 for b in bounds) - .09
            if column == 3 and show_method_notes:
                caption = fig.text(x / geometry["width"], top / geometry["height"],
                                   wrap("Top: baseline\nMiddle: post-rescue\nBottom: post-refinement", size * 72 - 1.5, font),
                                   va="top", fontproperties=font, linespacing=1.35)
                bottom = min(bottom, *(caption.get_window_extent(r).y0 / 72 for r in renderers))
        notes_font = FontProperties(family="sans-serif", size=7.1)
        footer_top = .15
        if show_method_notes:
            notes = "\n".join(wrap(line, (geometry["width"] - .3) * 72, notes_font) for line in method_notes(data))
            footer = fig.text(.15 / geometry["width"], .15 / geometry["height"], notes,
                              va="bottom", fontproperties=notes_font, linespacing=1.35)
            footer_top = max(footer.get_window_extent(r).y1 / 72 for r in renderers)
        if footer_top >= bottom - .07:
            extra = footer_top + .09 - bottom
            if _expanded or not math.isfinite(extra):
                raise ValueError("Figure legends overlap the recorded-threshold footer")
            # Text metrics vary across font installations and renderers. Grow
            # only the legend reserve; keep panel widths, row pitch and bars in
            # physical units, then verify the complete layout again.
            for artist in [*fig.legends[legend_start:], *fig.texts[text_start:]]:
                artist.remove()
            geometry = {**geometry, "height": geometry["height"] + extra,
                        "bottom": geometry["bottom"] + extra}
            fig.set_size_inches(geometry["width"], geometry["height"])
            for ax, (x, width) in zip(fig.axes, geometry["positions"], strict=True):
                ax.set_position([x / geometry["width"], geometry["bottom"] / geometry["height"],
                                 width / geometry["width"], geometry["plot_height"] / geometry["height"]])
            return _place_legends(
                fig, geometry, data, omark=omark, show_method_notes=show_method_notes, _expanded=True
            )
    finally:
        pdf.close()


def _refinement_counts(paths):
    coding = paths["repairs"], paths["isoforms"]
    structure = paths.get("architecture")
    return (*coding, structure["merges"], structure["splits"]) if structure is not None else coding


def plot_method_summary(data, output, *, omark=None, omark_label="", show_method_notes=False, gene_totals=None):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import MaxNLocator

    rows = validate(data)
    names = [row["species"] for row in rows]
    if omark is not None and not set(names).intersection(omark):
        raise ValueError("No species in the figure match the OMArk comparison")
    g = layout(len(rows), omark=omark is not None, labelled=bool(omark_label),
               show_gene_totals=gene_totals is not None)
    if not show_method_notes:
        g["height"] -= .4
        g["bottom"] -= .4
    groups = data["groups"]
    with plt.rc_context(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Nimbus Sans", "DejaVu Sans"],
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "pdf.use14corefonts": True,
            "axes.linewidth": 0.5,
            "svg.fonttype": "none",
        }
    ):
        fig = plt.figure(figsize=(g["width"], g["height"]), dpi=72)
        axes = []
        for x, w in g["positions"]:
            ax = fig.add_axes(
                [x / g["width"], g["bottom"] / g["height"], w / g["width"], g["plot_height"] / g["height"]],
                sharey=axes[0] if axes else None,
            )
            axes.append(ax)
            ax.set_ylim(len(rows) - 0.55, -0.55)
            ax.spines[["top", "right", "left"]].set_visible(False)
            ax.tick_params(axis="y", length=0, pad=4)
            ax.tick_params(axis="x", length=2, width=0.5, pad=3)
        rescue_ax, path_ax, utr_ax, busco_ax, *omark_axes = axes
        axes[0].set_yticks(range(len(rows)), [r["species"].replace("_", "\n", 1) for r in rows])
        for label in axes[0].get_yticklabels():
            label.set_fontstyle("italic")
        for ax in axes[1:]:
            ax.tick_params(axis="y", labelleft=False)
        maximum_rescue = max(
            (max([*r["rescue"]["method_counts"].values(), r["rescue"]["unique_loci"]]) for r in rows if r["analysed"]),
            default=0,
        )
        maximum_paths = max(
            (sum(_refinement_counts(row["coding_paths"])) for row in rows if row["analysed"]),
            default=0,
        )
        rescue_ax.set_xlim(0, max(1, maximum_rescue) * 2.3)
        path_ax.set_xlim(0, max(1, maximum_paths) * 1.05)
        for ax in (rescue_ax, path_ax):
            ax.xaxis.set_major_locator(MaxNLocator(nbins=2, integer=True))
        _count_axis(busco_ax, groups)
        _count_axis(utr_ax, max((r["utr"]["final_genes"] for r in rows if r["analysed"]), default=0))
        for ax, label in ((utr_ax, "With UTR"), (busco_ax, "Complete")):
            ax.text(0.99, 1.015, label, transform=ax.transAxes, ha="right", va="bottom", fontsize=7)
        for i, row in enumerate(rows):
            if not row["analysed"]:
                # Zero baselines are drawn only for analysed rows below. These
                # placeholders have neither a zero line nor a numeric bar.
                for ax in (rescue_ax, path_ax, utr_ax):
                    ax.add_patch(
                        Rectangle(
                            (0.035, i - 0.38),
                            0.93,
                            0.76,
                            transform=ax.get_yaxis_transform(),
                            facecolor="#eef0f3",
                            edgecolor="none",
                        )
                    )
                    ax.text(0.07, i, "Not analyzed", transform=ax.get_yaxis_transform(),
                            va="center", color="#657585")
            phases = PHASES if row["analysed"] else PHASES[:1]
            if not row["analysed"] and row.get("busco_refined_assessed"):
                phases = ("pre_rescue", "after")
            for phase, offset in zip(PHASES, bar_offsets(3), strict=True):
                if phase not in phases:
                    continue
                result = row["busco"][phase]
                left = 0
                for key, colour in zip(
                    ["single", "duplicated", "fragmented", "missing"],
                    ["#000000", "#B22222", "#666666", "#CCCCCC"],
                    strict=True,
                ):
                    count = result[key]
                    _bar(busco_ax, i + offset, count, colour, left=left)
                    left += count
                _annotate_percent(busco_ax, i + offset, result["single"] + result["duplicated"], groups)
            if not row["analysed"]:
                continue
            rescue = row["rescue"]
            top, middle, bottom = bar_offsets(3)
            upper, lower = bar_offsets(2)
            left = 0
            for amount, colour in zip(miniprot_partition(rescue), MINIPROT_STACK_COLOURS, strict=True):
                _bar(rescue_ax, i + top, amount, colour, left=left)
                left += amount
            rescue_ax.annotate(str(left), (left, i + top), xytext=(2, 0), textcoords="offset points", va="center")
            amount = rescue["method_counts"]["oriongeno"]
            _bar(rescue_ax, i + middle, amount, METHOD_COLOURS[2])
            rescue_ax.annotate(
                str(amount), (amount, i + middle), xytext=(2, 0), textcoords="offset points", va="center"
            )
            support = rescue.get("swissprot_support")
            if support is None:
                support = {**dict.fromkeys(SWISS_KEYS, 0), "not_assessed": rescue["unique_loci"]}
            left = 0
            for key, colour in zip(SWISS_KEYS, SWISS_COLOURS, strict=True):
                _bar(rescue_ax, i + bottom, support[key], colour, left=left)
                left += support[key]
            rescue_ax.annotate(
                str(left), (left, i + bottom), xytext=(2, 0), textcoords="offset points", va="center"
            )
            counts = _refinement_counts(row["coding_paths"])
            colours = (*PATH_CHANGE_COLOURS, *ARCHITECTURE_COLOURS)
            left = 0
            for count, colour, label_x in zip(counts, colours, (0, 0.31, 0.54, 0.77), strict=False):
                _bar(path_ax, i + upper, count, colour, left=left)
                path_ax.text(label_x, i + lower, f"{count:,}", transform=path_ax.get_yaxis_transform(),
                             color=colour, va="center", fontsize=7.5)
                left += count
            value = row["utr"]
            n = value["final_genes"]
            left = 0
            for key, colour in zip(UTR_KEYS, UTR_COLOURS, strict=True):
                count = value["sides"][key]
                _bar(utr_ax, i + upper, count, colour, left=left)
                left += count
            _annotate_percent(utr_ax, i, n - value["sides"]["none"], n)
            left = 0
            for key, colour in zip((*ORIGIN_KEYS, "none"), (*ORIGIN_COLOURS, UTR_COLOURS[-1]), strict=True):
                count = value["origins"][key]
                _bar(utr_ax, i + lower, count, colour, left=left)
                left += count
        titles = [
            "Missing gene\nrescue",
            "Refined coding\npaths",
            "UTRs in final\ngene models",
            "BUSCO groups\n" + data["busco_dataset"] if data.get("busco_dataset") else "BUSCO\ngroups",
        ]
        labels = ["Gene loci", "Coding paths", "Final transcripts", "Groups"]
        if any(r["analysed"] and "architecture" in r["coding_paths"] for r in rows):
            titles[1] = "Gene model\nrefinement"
            labels[1] = "Paths or events"
        if omark_axes:
            titles.extend(("OMArk conserved\nHOGs", "OMArk consistent\nand full-length"))
            labels.extend(("HOGs", "Proteins"))
            for metric, (ax, colours) in enumerate(zip(omark_axes, (OMARK_HOG_COLOURS, OMARK_PROTEIN_COLOURS), strict=True)):
                _plot_omark_counts(ax, names, omark, metric, colours, gene_totals if metric == 1 else None)
        for ax, title, label in zip(axes, titles, labels, strict=True):
            ax.set_title(title, loc="left", pad=32)
            ax.set_xlabel(label, labelpad=5)
        _place_legends(fig, g, data, omark=omark is not None, show_method_notes=show_method_notes)
        if omark_label:
            fig.text(.15 / g["width"], .992, omark_label, va="top", color="#465a6a")
        output = Path(output)
        output.mkdir(parents=True, exist_ok=True)
        try:
            for suffix in ("pdf", "png", "svg"):
                fig.savefig(output / ("busco_comparison." + suffix), dpi=200)
        finally:
            plt.close(fig)
