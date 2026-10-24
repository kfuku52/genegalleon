"""Publication-width, three-stage gene-model improvement figure.

Inputs are the verified BUSCO rows and model-change summary. Rescue counts
gene loci; refinement counts accepted coding paths. Swiss-Prot protein
support stays separate from genomic repeat-overlap evidence.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from gene_model_improvement_summary import plot_method_summary
from matplotlib.text import Text

STATUS_COLOURS = ("#000000", "#B22222", "#666666", "#CCCCCC")
DONOR_KEYS = ("self_only", "relative_only", "phylogenetic_only", "multiple")
DONOR_COLOURS = ("#b59b5b", "#5275b5", "#9a6ab2", "#37956f")
DONOR_LABELS = ("Self only", "Nearest only", "Balanced only", "Multiple groups")
SWISS_KEYS = ("te_only", "other_only", "both", "no_informative_hit", "not_assessed")
SWISS_COLOURS = ("#b45158", "#477b80", "#9a6ab2", "#c7cdd4", "#edf0f3")
SWISS_LABELS = ("TE-related only", "Other protein only", "TE + other", "No informative hit", "Not assessed")


def _wrap(text, width_points):
    from matplotlib._afm import AFM
    with (Path(matplotlib.get_data_path()) / "fonts/pdfcorefonts/Helvetica.afm").open("rb") as handle:
        font = AFM(handle)
    lines, line = [], ""
    for word in text.split():
        proposed = (line + " " + word).strip()
        if line and font.string_width_height(proposed)[0] * 8 / 1000 > width_points:
            lines.append(line)
            line = word
        else:
            line = proposed
    if line:
        lines.append(line)
    return lines


def _donors(value, field, path=False):
    keys = (*DONOR_KEYS, "other_interspecies") if path else DONOR_KEYS
    grouped = value.get(field + "_groups")
    if grouped is not None:
        return keys, grouped
    original = value.get(field + "_counts")
    if original is None:
        return (), {}
    if path:
        return keys, dict(zip(keys, (original[k] for k in
                                    ("self_only", "nearest_only", "balanced_only", "both", "other_interspecies")), strict=True))
    return keys, dict(zip(keys, (value.get("rescue_self_only_loci", 0), original["nearest_only"],
                                original["balanced_only"], original["both"]), strict=True))


def plot_staged_comparison(rows, output, changes):
    """Render four equal panels; preserve signed changes and unavailable rows."""
    totals = {r[p + "_result"]["total"] for r in rows for p in ("pre_rescue", "before", "after")}
    if len(totals) != 1 or next(iter(totals)) <= 0:
        raise ValueError("Count axes require the same positive BUSCO group total")
    total = next(iter(totals))
    stats = changes["species"]
    with plt.rc_context({"font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Nimbus Sans", "DejaVu Sans"],
                         "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "svg.fonttype": "none",
                         "pdf.use14corefonts": True, "axes.linewidth": .5, "hatch.linewidth": .4}):
        width, height = 7.2, max(8.6, .3291666667 * len(rows) + 6.25)
        fig = plt.figure(figsize=(width, height))
        panel_width = (width - 1.15 - .15 - 3 * .28) / 4
        positions = [(1.15 + i * (panel_width + .28), panel_width) for i in range(4)]
        bottom, plot_height = 5.0, height - 6.25
        axes = []
        for x, size in positions:
            ax = fig.add_axes([x / width, bottom / height, size / width, plot_height / height],
                              sharey=axes[0] if axes else None)
            axes.append(ax)
            ax.spines[["top", "right", "left"]].set_visible(False)
            ax.tick_params(axis="y", length=0, pad=3)
            ax.tick_params(axis="x", length=2, width=.5, pad=3)
            ax.grid(axis="x", alpha=.15, linewidth=.4)
            ax.set_axisbelow(True)
            ax.set_ylim(len(rows) - .6, -.6)
        rescue_axis, path_axis, busco_axis, delta_axis = axes
        for ax in axes[1:]:
            ax.tick_params(axis="y", labelleft=False)
        axes[0].set_yticks(range(len(rows)), [r["species"].replace("_", "\n", 1) for r in rows])
        for label in axes[0].get_yticklabels():
            label.set_fontstyle("italic")
        change_min, change_max = 0, 0
        for i, row in enumerate(rows):
            value = stats[row["species"]]
            unavailable = row["refinement_status"] == "not_analysed"
            if unavailable:
                for ax in axes:
                    ax.axhspan(i - .48, i + .48, color="#eef0f3", zorder=0)
            for phase, offset in (("pre_rescue", -.26), ("before", 0), ("after", .26)):
                if unavailable and phase != "pre_rescue":
                    continue
                result, left = row[phase + "_result"], 0
                for category, colour in zip(("single", "duplicated", "fragmented", "missing"), STATUS_COLOURS, strict=True):
                    amount = 100 * result[category] / total
                    busco_axis.barh(i + offset, amount, left=left, height=.18, color=colour)
                    left += amount
                busco_axis.text(101, i + offset, f'{100 * result["complete"] / total:.2f}', va="center")
            if unavailable:
                for ax in (rescue_axis, path_axis, delta_axis):
                    ax.text(.03, i, "Not analysed", transform=ax.get_yaxis_transform(), va="center", color="#657585")
                continue
            gain, refinement = row["delta_rescue_complete_pp"], row["delta_complete_pp"]
            # Same-sign gains/losses stack. Opposite signs extend independently
            # from zero, so a loss cannot be hidden behind a preceding gain.
            start = gain if gain * refinement >= 0 else 0
            delta_axis.barh(i, gain, height=.65, color="#37956f")
            delta_axis.barh(i, refinement, left=start, height=.65, color="#5275b5")
            change_min = min(change_min, gain, start, start + refinement)
            change_max = max(change_max, gain, start, start + refinement)
            rescued = value["prior_rescued_loci"]
            keys, donors = _donors(value, "rescue_support")
            left = 0
            if keys:
                for key, colour in zip(keys, DONOR_COLOURS, strict=True):
                    rescue_axis.barh(i - .20, donors[key], left=left, height=.32, color=colour)
                    left += donors[key]
            else:
                rescue_axis.barh(i - .20, rescued, height=.32, color="#5275b5")
            rescue_axis.annotate(str(rescued), (rescued, i - .20), xytext=(2, 0), textcoords="offset points", va="center")
            proteins = value.get("rescue_swissprot_groups")
            if proteins is None:
                proteins = {**dict.fromkeys(SWISS_KEYS, 0), "not_assessed": rescued}
            left = 0
            for key, colour in zip(SWISS_KEYS, SWISS_COLOURS, strict=True):
                rescue_axis.barh(i + .20, proteins[key], left=left, height=.32, color=colour,
                                 hatch="///" if key == "not_assessed" else None)
                left += proteins[key]
            rescue_axis.annotate(str(rescued), (rescued, i + .20), xytext=(2, 0), textcoords="offset points", va="center")
            repair, isoform = value["accepted_repair_paths"], value["accepted_isoform_paths"]
            path_axis.barh(i - .20, repair, height=.32, color="#187d97")
            path_axis.barh(i - .20, isoform, left=repair, height=.32, color="#d38b21")
            path_axis.annotate(f"{repair} / {isoform}", (repair + isoform, i - .20), xytext=(2, 0), textcoords="offset points", va="center")
            keys, donors = _donors(value, "accepted_path_support", path=True)
            left = 0
            for key, colour in zip(keys, (*DONOR_COLOURS, "#8a929b")[:len(keys)], strict=True):
                path_axis.barh(i + .20, donors[key], left=left, height=.32, color=colour)
                left += donors[key]
            if keys:
                path_axis.annotate(str(left), (left, i + .20), xytext=(2, 0), textcoords="offset points", va="center")
        titles = ("Missing-gene\nrescue", "Gene-model\nrefinement",
                  "BUSCO: before / after rescue /\nafter refinement", "Change in\ncomplete groups")
        for ax, title in zip(axes, titles, strict=True):
            ax.set_title(title, pad=30, fontweight="normal")
        busco_axis.set_xlim(0, 116)
        busco_axis.set_xticks([0, 50, 100])
        busco_axis.set_xlabel("BUSCO completeness (%)", labelpad=5)
        delta_axis.set_xlim(min(0, change_min * 1.15), max(1, change_max * 1.15))
        delta_axis.set_xlabel("BUSCO completeness\nchange (% points)", labelpad=5)
        for ax in (busco_axis, delta_axis):
            count_axis = ax.secondary_xaxis("top", functions=(lambda pp: pp * total / 100, lambda n: n * 100 / total))
            count_axis.set_xlabel("BUSCO group count", labelpad=5)
            count_axis.tick_params(length=2, width=.5, pad=3)
            from matplotlib.ticker import MaxNLocator
            count_axis.xaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
        rescue_max = max(v["prior_rescued_loci"] or 0 for v in stats.values())
        path_max = max((v["accepted_repair_paths"] or 0) + (v["accepted_isoform_paths"] or 0) for v in stats.values())
        rescue_axis.set_xlim(0, max(1, rescue_max) * 1.42)
        path_axis.set_xlim(0, max(1, path_max) * 1.95)
        for ax in (rescue_axis, path_axis):
            ax.xaxis.set_major_locator(MaxNLocator(nbins=2, integer=True))
        rescue_axis.set_xlabel("Rescued gene loci\n(count)", labelpad=5)
        path_axis.set_xlabel("Accepted coding paths\n(count)", labelpad=5)
        fig.text(.15 / width, (height - .20) / height, "BUSCO completeness and gene-model improvement", va="top")
        identity = rows[0]["before_result"]
        mode = identity.get("mode", "transcriptome")
        fig.text(.15 / width, (height - .39) / height,
                 f'BUSCO {identity["busco_version"]}; {identity["lineage"]}; n = {total}; {mode} mode', va="top")

        def legend(index, top, title, keys, colours, labels):
            x = positions[index][0]
            y = top
            for line in title.split("\n"):
                fig.text(x / width, y / height, line, va="top")
                y -= 10 / 72
            y -= 3 / 72
            for key, colour, label in zip(keys, colours, labels, strict=True):
                hatch = "///" if key == "not_assessed" else None
                fig.patches.append(plt.Rectangle((x / width, (y - 6 / 72) / height), .085 / width, 5 / 72 / height,
                                                 transform=fig.transFigure, facecolor=colour, edgecolor="#555555" if hatch else "none",
                                                 linewidth=.2, hatch=hatch))
                for line in _wrap(label, (panel_width - .13) * 72):
                    fig.text((x + .13) / width, y / height, line, va="top")
                    y -= 10 / 72
                y -= 2 / 72
            return y

        top = 4.33
        end = legend(0, top, "Donor groups\n(upper bars)", DONOR_KEYS, DONOR_COLOURS, DONOR_LABELS)
        swiss_end = legend(0, end - .13, "Swiss-Prot support\n(lower bars)", SWISS_KEYS, SWISS_COLOURS, SWISS_LABELS)
        end = legend(1, top, "Coding-path type\n(upper bars)", ("repair", "isoform"), ("#187d97", "#d38b21"), ("Repair path", "Added isoform"))
        path_end = legend(1, end - .13, "Donor groups\n(lower bars)", (*DONOR_KEYS, "other_interspecies"),
                          (*DONOR_COLOURS, "#8a929b"), (*DONOR_LABELS, "Other species only"))
        legend(2, top, "BUSCO status", ("single", "duplicated", "fragmented", "missing"), STATUS_COLOURS,
               ("Single-copy", "Duplicated", "Fragmented", "Missing"))
        legend(3, top, "Complete groups\n(gain)", ("rescue", "refinement"), ("#37956f", "#5275b5"), ("Rescue", "Refinement"))
        notes = ["BUSCO bars: before rescue, after rescue, after refinement (top to bottom); labels are Complete (%).",
                 "Grey row: rescue/refinement not analysed; only the pre-rescue BUSCO bar is shown.",
                 f"Change bars stack rescue and refinement; both count axes use the same {total:,} BUSCO groups.",
                 "Rescue counts added loci. Refinement counts accepted repair / isoform paths, which may share a locus.",
                 "Donor groups: self, nearest relatives, balanced references or multiple groups; target RNA is separate.",
                 "Multiple = at least two donor groups; a donor in both reference lists supplies both groups."]
        parameters = changes.get("swissprot_evidence", {}).get("parameters", {})
        if all(k in parameters for k in ("evalue", "minimum_alignment", "query_coverage", "target_coverage", "score_fraction")):
            notes.append(f'Swiss-Prot: E <= {parameters["evalue"]:g}; paired length >= {parameters["minimum_alignment"]:g} aa; '
                         f'query/target coverage >= {100 * parameters["query_coverage"]:g}%/{100 * parameters["target_coverage"]:g}%; '
                         f'score >= {100 * parameters["score_fraction"]:g}% of best qualifying hit.')
        elif any(v.get("rescue_swissprot_groups") is not None for v in stats.values()):
            notes.append("Swiss-Prot thresholds are unavailable in the supplied historical audit.")
        else:
            notes.append("Swiss-Prot audit not supplied: protein support is not assessed; repeat overlap is a separate measurement.")
        notes.append("Swiss-Prot support is advisory; other support does not establish host function; no support does not exclude TE origin.")
        lines = [line for text in notes for line in _wrap(text, (width - .30) * 72)]
        note_top = .16 + len(lines) * 10 / 72
        if min(swiss_end, path_end) <= note_top + .10:
            raise ValueError("Figure legends overlap the methods footer")
        for i, line in enumerate(lines):
            fig.text(.15 / width, (note_top - i * 10 / 72) / height, line, va="top")
        for text in fig.findobj(match=Text):
            text.set_fontsize(8)
        output = Path(output)
        output.mkdir(parents=True, exist_ok=True)
        for suffix in ("png", "svg", "pdf"):
            fig.savefig(output / ("busco_comparison." + suffix), dpi=180, facecolor="white")
        plt.close(fig)


def main():
    """Redraw the five-column method/UTR view without running scientific stages."""
    import argparse
    import hashlib
    import json

    from input_generation_array_state import FreshDigestBatch, atomic_json

    parser = argparse.ArgumentParser(description="Draw a compact gene-model method and UTR summary")
    parser.add_argument("--figure-data", required=True, type=Path,
                        help="Explicit rescue/path/BUSCO/UTR reporting counts, not a model-admission input")
    parser.add_argument("--output", required=True, type=Path, help="Separate figure output directory")
    parser.add_argument("--refinement-output", type=Path,
                        help="Reclassify path methods from receipt-bound refinement predictions")
    args = parser.parse_args()
    output = args.output.resolve()
    if args.refinement_output is not None:
        root = args.refinement_output.resolve()
        if root == output or root in output.parents or output in root.parents:
            parser.error("Figure output must be separate from the immutable refinement tree")
    source = args.figure_data.resolve()
    targets = ("improvement_figure_data.json", "improvement_method_provenance.json",
               "busco_comparison.pdf", "busco_comparison.png", "busco_comparison.svg")
    if source in {output / name for name in targets}:
        parser.error("Figure output must not overwrite its reporting input")
    raw = source.read_bytes()
    source_sha = hashlib.sha256(raw).hexdigest()
    data = json.loads(raw)
    verification = FreshDigestBatch()
    verification.read([source], expected={str(source): source_sha})
    implementations = [Path(__file__), Path(__file__).with_name("gene_model_improvement_summary.py")]
    if args.refinement_output is not None:
        implementations.append(Path(__file__).with_name("refinement_path_methods.py"))
    source_hashes = verification.read(implementations)
    provenance = dict(figure_data_sha256=source_sha, predictors_executed=False, busco_executed=False,
                      renderer_sha256={Path(p).name: sha for p, sha in source_hashes.items()},
                      matplotlib_version=matplotlib.__version__)
    if args.refinement_output is not None:
        from refinement_path_methods import collect_path_methods
        data, methods = collect_path_methods(data, args.refinement_output, verification=verification)
        provenance["path_methods"] = methods
    plot_method_summary(data, output)
    verification.check()
    atomic_json(output / "improvement_figure_data.json", data)
    atomic_json(output / "improvement_method_provenance.json", provenance)


if __name__ == "__main__":
    main()
