"""Five-column gene-model summary from explicit, reconciled reporting counts.

Rescue counts distinct loci, refinement counts coding paths, and UTR percentages
use final exported gene models. Plotting never runs predictors or BUSCO.
"""

import textwrap
from pathlib import Path

METHODS = ("synteny_based", "miniprot", "oriongeno")
METHOD_COLOURS = ("#5275b5", "#37956f", "#9a6ab2")
MINIPROT_STACK_COLOURS = ("#5275b5", "#37956f", "#6b7987")
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
        height=ROW_PITCH_IN * rows + 3.35,
        bottom=2.6,
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


def plot_method_summary(data, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch, Rectangle
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
        fig = plt.figure(figsize=(g["width"], g["height"]))
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
            for phase, offset in zip(phases, bar_offsets(len(phases)), strict=True):
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
            amount = rescue["unique_loci"]
            rescue_ax.barh(i + bottom, amount, height=BAR_HEIGHT, color="#111111")
            rescue_ax.annotate(
                str(amount), (amount, i + bottom), xytext=(2, 0), textcoords="offset points", va="center"
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
            "BUSCO\ngroups",
            "Change in complete\ngroups (S+D)",
            "UTRs in final\ngene models",
        ]
        labels = ["Gene loci", "Coding paths", "Groups (%)", "Groups (%)", "Final genes (%)"]
        for ax, title, label in zip(axes, titles, labels, strict=True):
            ax.set_title(title, loc="left", pad=32)
            ax.set_xlabel(label, labelpad=5)
        legends = [
            list(
                zip(
                    MINIPROT_STACK_COLOURS,
                    ["Miniprot (synteny-guided)", "Miniprot (synteny-independent)", "Both Miniprot routes"],
                    strict=True,
                )
            )
            + [(METHOD_COLOURS[2], "OrionGeno"), ("#111111", "Total accepted loci")],
            list(zip(PATH_CHANGE_COLOURS, ["Model revision", "Isoform addition"], strict=True))
            + list(zip(data["coding_path_legend"]["colours"], data["coding_path_legend"]["labels"], strict=True)),
            list(
                zip(
                    ["#000000", "#B22222", "#666666", "#CCCCCC"],
                    ["Single-copy", "Duplicated", "Fragmented", "Missing"],
                    strict=True,
                )
            ),
            [("#37956f", "Rescue"), ("#5275b5", "Refinement")],
            list(zip(UTR_COLOURS, ["5' only", "3' only", "Both ends", "No UTR"], strict=True))
            + list(zip(ORIGIN_COLOURS, ["Original GFF only", "OrionGeno only", "Both sources"], strict=True)),
        ]
        for (x, w), items in zip(g["positions"], legends, strict=True):
            handles = [Patch(facecolor=colour, label="\n".join(textwrap.wrap(label, 17))) for colour, label in items]
            fig.legend(
                handles=handles,
                loc="upper left",
                bbox_to_anchor=(x / g["width"], 2.1 / g["height"]),
                frameon=False,
                fontsize=8,
                handlelength=0.7,
                handletextpad=0.4,
                borderaxespad=0,
                labelspacing=0.4,
            )
        fig.text(
            0.03,
            0.08 / g["height"],
            "BUSCO: baseline, rescue, refinement (top to bottom); unanalysed: baseline only.\nMiniprot route segments are exclusive; OrionGeno can overlap them. Black = total accepted loci.\nUTR: ends (top), sources (bottom); n = final gene models; grey = no UTR.",
            fontsize=8,
            va="bottom",
        )
        output = Path(output)
        output.mkdir(parents=True, exist_ok=True)
        try:
            for suffix in ("pdf", "png", "svg"):
                fig.savefig(output / ("busco_comparison." + suffix), dpi=200)
        finally:
            plt.close(fig)
