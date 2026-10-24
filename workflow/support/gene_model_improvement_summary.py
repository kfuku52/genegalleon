"""Five-column gene-model summary from explicit, reconciled reporting counts.

Rescue counts distinct loci, refinement counts coding paths, and UTR percentages
use final exported gene models. Plotting never runs predictors or BUSCO.
"""

import math
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
PATH_CHANGE_COLOURS = (UTR_COLOURS[0], UTR_COLOURS[2])
PATH_METHOD_LABELS = ("Miniprot only", "OrionGeno only", "Both methods")
PHASES = ("pre_rescue", "before", "after")
ROW_PITCH_IN = 0.35
BAR_HEIGHT_PT = 4.4
BAR_STEP_PT = 8.0
BAR_HEIGHT = BAR_HEIGHT_PT / (72 * ROW_PITCH_IN)


def bar_offsets(count):
    if count not in (1, 2, 3):
        raise ValueError("The figure uses one, two or three bars per species")
    step = BAR_STEP_PT / (72 * ROW_PITCH_IN)
    return tuple((index - (count - 1) / 2) * step for index in range(count))


def layout(rows):
    if type(rows) is not int or rows < 1:
        raise ValueError("The figure requires at least one species")
    width, left, right, gap = 7.2, 1.13, 0.14, 0.12
    panel = (width - left - right - 4 * gap) / 5
    return dict(
        width=width,
        height=ROW_PITCH_IN * rows + 4.80,
        bottom=4.05,
        plot_height=ROW_PITCH_IN * rows,
        positions=[(left + i * (panel + gap), panel) for i in range(5)],
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
        for phase in PHASES if row["analysed"] else PHASES[:1]:
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
            len(paths["counts"]) != len(paths["colours"])
            or paths["colours"] != data["coding_path_legend"]["colours"]
            or len(paths["counts"]) != len(data["coding_path_legend"]["labels"])
            or any(type(n) is not int or n < 0 for n in paths["counts"])
            or any(type(paths[k]) is not int or paths[k] < 0 for k in ["repairs", "isoforms"])
            or sum(paths["counts"]) != paths["repairs"] + paths["isoforms"]
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


def _place_legends(fig, geometry, data):
    """Fit each bar-group legend using both PNG and PDF font metrics."""
    import io

    from matplotlib.backends.backend_pdf import PdfFile, RendererPdf
    from matplotlib.font_manager import FontProperties
    from matplotlib.patches import Patch

    font = FontProperties(family="sans-serif", size=7.7)
    bold = FontProperties(family="sans-serif", size=7.7, weight="bold")
    # At 72 dpi, the figure transforms and PDF renderer both use points. Agg
    # alone underestimates long legend stacks when the PDF uses core fonts.
    pdf = PdfFile(io.BytesIO())
    renderer = RendererPdf(pdf, 72, geometry["height"], geometry["width"])
    fig.canvas.draw()
    renderers = (renderer, fig.canvas.get_renderer())

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
        [("Model changes", list(zip(PATH_CHANGE_COLOURS, ("Model revision", "Isoform addition"), strict=True))),
         ("Nomination method", list(zip(ORIGIN_COLOURS, PATH_METHOD_LABELS, strict=True)))],
        [("BUSCO categories", busco)],
        [("Complete-group gain", [("#37956f", "Rescue"), ("#5275b5", "Refinement")])],
        [("UTR ends", list(zip(UTR_COLOURS, ("5' only", "3' only", "Both ends", "No UTR"), strict=True))),
         ("UTR sources", list(zip((*ORIGIN_COLOURS, UTR_COLOURS[-1]),
                                  ("Original GFF only", "OrionGeno only", "Both sources", "No UTR"), strict=True)))],
    ]
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
            if column == 2:
                fig.text(x / geometry["width"], top / geometry["height"],
                         wrap("Top: baseline\nMiddle: post-rescue\nBottom: post-refinement", size * 72 - 1.5, font),
                         va="top", fontproperties=font, linespacing=1.35)
        notes_font = FontProperties(family="sans-serif", size=7.1)
        notes = "\n".join(wrap(line, (geometry["width"] - .3) * 72, notes_font) for line in method_notes(data))
        footer = fig.text(.15 / geometry["width"], .15 / geometry["height"], notes,
                          va="bottom", fontproperties=notes_font, linespacing=1.35)
        if any(footer.get_window_extent(r).y1 / 72 >= bottom - .07 for r in renderers):
            raise ValueError("Figure legends overlap the recorded-threshold footer")
    finally:
        pdf.close()


def plot_method_summary(data, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import MaxNLocator, StrMethodFormatter

    rows = validate(data)
    g = layout(len(rows))
    deltas = [0]
    groups = data["groups"]
    for row in rows:
        if row["analysed"]:
            complete = [sum(row["busco"][p][k] for k in ["single", "duplicated"]) for p in PHASES]
            gain, refine = (complete[1] - complete[0]) * 100 / groups, (complete[2] - complete[1]) * 100 / groups
            start = gain if gain * refine >= 0 else 0
            deltas.extend([0, gain, start, start + refine])
    low, high = min(deltas), max(deltas)
    span = max(1, high - low)
    low = low - 0.04 * span if low < 0 else 0
    high = high + 0.12 * span if high else 1
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
        rescue_ax, path_ax, busco_ax, delta_ax, utr_ax = axes
        axes[0].set_yticks(range(len(rows)), [r["species"].replace("_", "\n", 1) for r in rows])
        for label in axes[0].get_yticklabels():
            label.set_fontstyle("italic")
        for ax in axes[1:]:
            ax.tick_params(axis="y", labelleft=False)
        maximum_rescue = max(
            (max([*r["rescue"]["method_counts"].values(), r["rescue"]["unique_loci"]]) for r in rows if r["analysed"]),
            default=0,
        )
        maximum_paths = max((sum(r["coding_paths"]["counts"]) for r in rows if r["analysed"]), default=0)
        rescue_ax.set_xlim(0, max(1, maximum_rescue) * 2.3)
        path_ax.set_xlim(0, max(1, maximum_paths) * 2.15)
        busco_ax.set_xlim(0, 150)
        busco_ax.set_xticks([0, 50, 100])
        utr_ax.set_xlim(0, 160)
        utr_ax.set_xticks([0, 50, 100])
        delta_ax.set_xlim(low, high)
        for ax in (rescue_ax, path_ax):
            ax.xaxis.set_major_locator(MaxNLocator(nbins=2, integer=True))
        delta_ax.xaxis.set_major_locator(MaxNLocator(nbins=3, steps=[1, 2, 2.5, 5, 10]))
        delta_ax.xaxis.set_major_formatter(StrMethodFormatter("{x:g}"))
        for ax in (busco_ax, delta_ax):
            top = ax.secondary_xaxis(
                "top", functions=(lambda percent: percent * groups / 100, lambda count: count * 100 / groups)
            )
            top.tick_params(axis="x", length=2, width=0.5, pad=2, labelsize=8)
            top.spines["top"].set_linewidth(0.5)
            top.set_xlabel("Groups", labelpad=3, fontsize=8)
            if ax is busco_ax:
                ticks = MaxNLocator(nbins=3, integer=True, steps=[1, 2, 2.5, 5, 10]).tick_values(0, groups)
                top.set_xticks([v for v in ticks if 0 <= v <= groups])
                top.spines["top"].set_bounds(0, groups)
                ax.spines["bottom"].set_bounds(0, 100)
            else:
                top.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
        for i, row in enumerate(rows):
            if not row["analysed"]:
                # Zero baselines are drawn only for analysed rows below. These
                # placeholders have neither a zero line nor a numeric bar.
                for ax in (rescue_ax, path_ax, delta_ax, utr_ax):
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
                    ax.text(0.07, i, "Not analysed", transform=ax.get_yaxis_transform(), va="center", color="#657585")
            phases = PHASES if row["analysed"] else PHASES[:1]
            for phase, offset in zip(phases, bar_offsets(3)[:len(phases)], strict=True):
                result = row["busco"][phase]
                left = 0
                for key, colour in zip(
                    ["single", "duplicated", "fragmented", "missing"],
                    ["#000000", "#B22222", "#666666", "#CCCCCC"],
                    strict=True,
                ):
                    amount = 100 * result[key] / groups
                    busco_ax.barh(i + offset, amount, left=left, height=BAR_HEIGHT, color=colour)
                    left += amount
                busco_ax.text(102, i + offset, str(result["single"] + result["duplicated"]), va="center")
            if not row["analysed"]:
                continue
            rescue = row["rescue"]
            top, middle, bottom = bar_offsets(3)
            upper, lower = bar_offsets(2)
            left = 0
            for amount, colour in zip(miniprot_partition(rescue), MINIPROT_STACK_COLOURS, strict=True):
                rescue_ax.barh(i + top, amount, left=left, height=BAR_HEIGHT, color=colour)
                left += amount
            rescue_ax.annotate(str(left), (left, i + top), xytext=(2, 0), textcoords="offset points", va="center")
            amount = rescue["method_counts"]["oriongeno"]
            rescue_ax.barh(i + middle, amount, height=BAR_HEIGHT, color=METHOD_COLOURS[2])
            rescue_ax.annotate(
                str(amount), (amount, i + middle), xytext=(2, 0), textcoords="offset points", va="center"
            )
            support = rescue.get("swissprot_support")
            if support is None:
                support = {**dict.fromkeys(SWISS_KEYS, 0), "not_assessed": rescue["unique_loci"]}
            left = 0
            for key, colour in zip(SWISS_KEYS, SWISS_COLOURS, strict=True):
                rescue_ax.barh(i + bottom, support[key], left=left, height=BAR_HEIGHT, color=colour)
                left += support[key]
            rescue_ax.annotate(
                str(left), (left, i + bottom), xytext=(2, 0), textcoords="offset points", va="center"
            )
            paths = row["coding_paths"]
            left = 0
            for amount, colour in zip([paths["repairs"], paths["isoforms"]], PATH_CHANGE_COLOURS, strict=True):
                path_ax.barh(i + upper, amount, left=left, height=BAR_HEIGHT, color=colour)
                left += amount
            path_ax.annotate(
                str(paths["repairs"]) + " / " + str(paths["isoforms"]),
                (left, i + upper),
                xytext=(2, 0),
                textcoords="offset points",
                va="center",
            )
            left = 0
            for amount, colour in zip(paths["counts"], paths["colours"], strict=True):
                path_ax.barh(i + lower, amount, left=left, height=BAR_HEIGHT, color=colour)
                left += amount
            path_ax.annotate(str(left), (left, i + lower), xytext=(2, 0), textcoords="offset points", va="center")
            complete = [sum(row["busco"][p][k] for k in ["single", "duplicated"]) for p in PHASES]
            gain, refine = (complete[1] - complete[0]) * 100 / groups, (complete[2] - complete[1]) * 100 / groups
            start = gain if gain * refine >= 0 else 0
            delta_ax.barh(i, gain, height=BAR_HEIGHT, color="#37956f")
            delta_ax.barh(i, refine, left=start, height=BAR_HEIGHT, color="#5275b5")
            delta_ax.plot([0, 0], [i - 0.42, i + 0.42], color="#777777", linewidth=0.5)
            value = row["utr"]
            n = value["final_genes"]
            left = 0
            for key, colour in zip(UTR_KEYS, UTR_COLOURS, strict=True):
                count = value["sides"][key]
                amount = 100 * count / n
                utr_ax.barh(i + upper, amount, left=left, height=BAR_HEIGHT, color=colour)
                left += amount
            utr_ax.text(102, i, "n=" + str(n), va="center", color="#333333")
            left = 0
            for key, colour in zip((*ORIGIN_KEYS, "none"), (*ORIGIN_COLOURS, UTR_COLOURS[-1]), strict=True):
                amount = 100 * value["origins"][key] / n
                utr_ax.barh(i + lower, amount, left=left, height=BAR_HEIGHT, color=colour)
                left += amount
        titles = [
            "Missing gene\nrescue",
            "Refined coding\npaths",
            "BUSCO groups\n" + data["busco_dataset"] if data.get("busco_dataset") else "BUSCO\ngroups",
            "Change in complete\ngroups (S+D)",
            "UTRs in final\ngene models",
        ]
        labels = ["Gene loci", "Coding paths", "Groups (%)", "Groups (%)", "Final genes (%)"]
        for ax, title, label in zip(axes, titles, labels, strict=True):
            ax.set_title(title, loc="left", pad=32)
            ax.set_xlabel(label, labelpad=5)
        _place_legends(fig, g, data)
        output = Path(output)
        output.mkdir(parents=True, exist_ok=True)
        try:
            for suffix in ("pdf", "png", "svg"):
                fig.savefig(output / ("busco_comparison." + suffix), dpi=200)
        finally:
            plt.close(fig)
