"""Generate compact RNA-seq QC from existing AMALGKIT merge outputs."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import tempfile
from collections import Counter, defaultdict
from hashlib import blake2s
from pathlib import Path

METRICS = (
    ("mapping_rate", "Mapping\nrate", "%", 1),
    ("fastp_duplication_rate", "Duplication\nrate (fastp)", "%", 1),
    ("fastp_insert_size_peak", "Insert size\npeak (fastp)", "bp", 1),
    ("total_spots", "Total\nspots", "million", 10**6),
    ("total_bases", "Total\nbases", "Gbp", 10**9),
)
METRIC_KEYS = tuple(metric[0] for metric in METRICS)
IDENTIFIERS = ("run", "scientific_name", "sample_group", "lib_layout")
SUMMARY_COLUMNS = (*IDENTIFIERS, *METRIC_KEYS, "mean_tpm")
TISSUE_COLORS = {"flower": "#cc79a7", "leaf": "#009e73", "root": "#0072b2", "other": "#53687b"}
LAYOUT_COLORS = {"paired": "#405d75", "single": "#aebfd0", "unknown": "#d8dbe0"}
INK = "#152a3a"
MUTED = "#465a6a"


def read_table(path, required):
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        header = reader.fieldnames or ()
        if len(header) != len(set(header)) or not set(required).issubset(header):
            raise ValueError(f"{path}: expected columns {', '.join(required)}")
        rows = list(reader)
        if any(None in row or None in row.values() for row in rows):
            raise ValueError(f"{path}: incomplete or extra TSV fields")
        return rows


def by_run(rows, path):
    records = {row["run"]: row for row in rows}
    if "" in records or len(records) != len(rows):
        raise ValueError(f"{path}: empty or duplicate RUN")
    return records


def mean_detected_tpm(path, runs):
    """Compute the per-RUN sum and positive count in one bounded, vectorized pass."""
    import numpy as np
    import pandas as pd

    sums = np.zeros(len(runs))
    detected = np.zeros(len(runs), dtype=np.int64)
    with pd.read_csv(path, sep="\t", usecols=runs, dtype=float, chunksize=8192) as chunks:
        for chunk in chunks:
            values = chunk.to_numpy(copy=False) if list(chunk) == runs else chunk[runs].to_numpy(copy=False)
            if not np.isfinite(values).all() or (values < 0).any():
                raise ValueError(f"{path}: TPM must be finite and nonnegative")
            sums += values.sum(axis=0)
            detected += np.count_nonzero(values, axis=0)
    return {run: total / count if count else None for run, total, count in zip(runs, sums, detected, strict=True)}


def collect_merge_qc(root, retained=None):
    """Read existing species-level metadata and TPM, without opening counts or FASTQs."""
    rows = []
    for metadata in sorted(root.glob("*/*_metadata.tsv")):
        source = by_run(read_table(metadata, ("run", "scientific_name")), metadata)
        if retained is not None and not retained.intersection(source):
            continue
        tpm = metadata.with_name(metadata.name.removesuffix("_metadata.tsv") + "_tpm.tsv")
        with tpm.open(encoding="utf-8", newline="") as stream:
            header = next(csv.reader(stream, delimiter="\t"), [])
        if len(header) < 2 or header[0] != "target_id" or len(header) != len(set(header)):
            raise ValueError(f"{tpm}: expected target_id and unique RUN columns")
        if missing := set(header[1:]) - source.keys():
            raise ValueError(f"{tpm}: RUNs absent from metadata: {', '.join(sorted(missing))}")
        runs = [run for run in header[1:] if source[run].get("exclusion", "no").strip().lower() in ("", "no")]
        if not runs:
            continue
        means = mean_detected_tpm(tpm, runs)
        for run in runs:
            row = {key: source[run].get(key, "") for key in SUMMARY_COLUMNS}
            row["mean_tpm"] = means[run]
            rows.append(row)
    if not rows:
        raise ValueError(f"{root}: no quantified RUNs found in species merge directories")
    return rows


def metric_value(row, field):
    value = row[field]
    if value is None or str(value).strip().lower() in ("", "na", "nan"):
        return None
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise ValueError(f"{row['run']}: invalid {field} value {value!r}")
    if field in ("mapping_rate", "fastp_duplication_rate") and number > 100:
        raise ValueError(f"{row['run']}: {field} must be in [0, 100]")
    return number


def write_summary(path, rows):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=path.parent, prefix=f".{path.name}.", mode="w", encoding="utf-8", newline="", delete=False,
        ) as stream:
            temporary = Path(stream.name)
            writer = csv.DictWriter(stream, SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def draw_qc(original, retained):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.font_manager import FontProperties
    from matplotlib.patches import Circle, Patch, Wedge
    from matplotlib.ticker import MaxNLocator

    names = sorted({row["scientific_name"] for row in retained})
    if not all(name.strip() for name in names):
        raise ValueError("Quantified RUNs must have a scientific_name")
    tissues = {
        row["run"]: row["sample_group"].lower() if row["sample_group"].lower() in TISSUE_COLORS else "other"
        for row in retained
    }
    tissue_names = [name for name in TISSUE_COLORS if name != "other" or name in tissues.values()]
    counts = Counter((row["scientific_name"], tissues[row["run"]]) for row in retained)
    layouts = Counter((
        row["scientific_name"], row["lib_layout"].lower() if row["lib_layout"].lower() in LAYOUT_COLORS else "unknown",
    ) for row in retained)
    means = defaultdict(list)
    for row in original:
        value = metric_value(row, "mean_tpm")
        if value:
            means[row["scientific_name"]].append(value)

    layout_names = [layout for layout in LAYOUT_COLORS if any(key[1] == layout for key in layouts)]
    width = 518.4  # 7.2 inches
    row_height, top = 10.5, 59
    bottom = 48 + 7 * max(0, len(layout_names) - 2)
    plot_height = row_height * len(names)
    height = bottom + plot_height + top
    fig = Figure(figsize=(width / 72, height / 72), dpi=160, facecolor="white")
    FigureCanvasAgg(fig)
    renderer = fig.canvas.get_renderer()
    centers = {name: (index + .5) * row_height for index, name in enumerate(names)}
    label_right = 131.5
    count_width = 20
    count_centers = [width - 14.4 - count_width * (len(tissue_names) - i - 1) for i in range(len(tissue_names))]
    mean_left = count_centers[0] - 61
    mean_width = 44
    pie_left, pie_width = 137.5, 28
    first, gap = 170, 10.5
    panel_width = (mean_left - first - len(METRICS) * gap) / len(METRICS)

    def text(x, y, label, size=6.6, color=MUTED, **style):
        return fig.text(x / width, y / height, label, fontsize=size, color=color, **style)

    def heading(x, label, unit=""):
        text(x, bottom + plot_height + 18, label, 6.45, INK, ha="center", va="bottom", weight="bold")
        text(x, bottom + plot_height + 5.8, unit, 5.75, ha="center", va="bottom")

    def axis(x, span, domain, ticks, ticklabels=None, log=False):
        ax = fig.add_axes([x / width, bottom / height, span / width, plot_height / height])
        if log:
            ax.set_xscale("log")
        ax.set(xlim=domain, ylim=(plot_height, 0), yticks=[])
        ax.set_xticks(ticks, labels=ticklabels)
        ax.minorticks_off()
        ax.tick_params(axis="x", labelsize=5.5, length=2, pad=2, labelcolor=MUTED, color="#8996a2")
        ax.spines[["top", "left", "right"]].set_visible(False)
        ax.spines["bottom"].set_color("#a0adb6")
        ax.spines["bottom"].set_linewidth(.6)
        ax.set_axisbelow(True)
        ax.grid(axis="x", color="#dce4ea", linewidth=.45)
        ax.hlines([index * row_height for index in range(len(names))], *domain, lw=.4, color="#dce4ea", zorder=0)
        return ax

    text(12, height - 12, "RNA-seq quality control", 13, INK, weight="bold", va="top")
    text(width - 12, height - 14.5, f"{len(names)} species", 7.6, ha="right", va="top")
    text(label_right, bottom + plot_height + 5.8, "Species", 7.2, INK, va="bottom", ha="right", weight="bold")
    for name, center in centers.items():
        parts = re.split(r"(?<=\s)(var\.|subsp\.|f\.)(?=\s)", name)
        styles = ["normal" if index % 2 else "italic" for index in range(len(parts))]
        widths = [
            renderer.get_text_width_height_descent(part, FontProperties(size=7.4, style=style), False)[0]
            * 72 / fig.dpi for part, style in zip(parts, styles, strict=True)
        ]
        scale = min(1, (label_right - 10) / sum(widths))
        x = label_right - sum(widths) * scale
        for part, style, span in zip(parts, styles, widths, strict=True):
            text(x, bottom + plot_height - center, part, 7.4 * scale, INK, fontstyle=style, va="center")
            x += span * scale

    text(sum(count_centers) / len(count_centers), bottom + plot_height + 20.3,
         "Runs by tissue", 7, INK, ha="center", weight="bold")
    for x, tissue in zip(count_centers, tissue_names, strict=True):
        text(x, bottom + plot_height + 5.8, tissue.capitalize(), 5.85, TISSUE_COLORS[tissue],
             ha="center", weight="bold")
        for name, center in centers.items():
            text(x, bottom + plot_height - center, str(counts[name, tissue]), 7.35, INK, va="center", ha="center")

    pie = axis(pie_left, pie_width, (0, pie_width), [])
    pie.set_aspect("equal")
    pie.set_axis_off()
    heading(pie_left + pie_width / 2, "Library\nlayout")
    for name, center in centers.items():
        total = sum(layouts[name, layout] for layout in layout_names)
        pieces = [(layout, layouts[name, layout]) for layout in layout_names if layouts[name, layout]]
        xy = (pie_width / 2, center)
        if len(pieces) == 1:
            pie.add_patch(Circle(xy, 3.8, color=LAYOUT_COLORS[pieces[0][0]], lw=0))
        else:
            start = 90.0
            for layout, count in pieces:
                end = start + count / total * 360
                pie.add_patch(Wedge(xy, 3.8, start, end, facecolor=LAYOUT_COLORS[layout],
                                    edgecolor="white", linewidth=.4))
                start = end

    # Stable RUN offsets across panels. These are display coordinates, never metric offsets.
    offsets = {}
    for row in retained:
        run = row["run"]
        value = int.from_bytes(blake2s(run.encode(), digest_size=4).digest(), "big")
        index = tissue_names.index(tissues[run])
        y = (index - (len(tissue_names) - 1) / 2) * (6.2 / max(2, len(tissue_names) - 1))
        offsets[run] = centers[row["scientific_name"]] + y + value / (2**32 - 1) - .5
    for i, (field, label, unit, divisor) in enumerate(METRICS):
        points = [(value / divisor, row["run"]) for row in retained if (value := metric_value(row, field)) is not None]
        maximum = max((value for value, _ in points), default=0)
        if unit == "%":
            limit, ticks = 100, [0, 50, 100]
        elif field == "fastp_insert_size_peak":
            limit = max(300, math.ceil(maximum / 50) * 50)
            ticks = [0, limit / 2, limit]
        else:
            limit = (maximum or 1) * 1.06
            ticks = [
                value for value in MaxNLocator(nbins=3, steps=[1, 2.5, 5, 10]).tick_values(0, limit)
                if value <= limit
            ]
        x = first + i * (panel_width + gap)
        ax = axis(x, panel_width, (0, limit), ticks)
        heading(x + panel_width / 2, label, unit)
        ax.scatter(
            [value for value, _ in points], [offsets[run] for _, run in points],
            c=[TISSUE_COLORS[tissues[run]] for _, run in points], s=1.8**2,
            edgecolors="white", linewidths=.16, zorder=4, clip_on=False,
        )

    group_names = [name for name in names if means[name]]
    groups = [means[name] for name in group_names]
    values = [value for group in groups for value in group]
    low, high = (min(values) * .9, max(values) * 1.15) if values else (1, 10)
    log_scale = high / low >= 5
    if log_scale:
        raw_ticks = [
            multiplier * 10.0**decade
            for decade in range(math.floor(math.log10(low)), math.ceil(math.log10(high)) + 1)
            for multiplier in (1, 3) if low <= multiplier * 10.0**decade <= high
        ]
        ticks = raw_ticks[::max(1, math.ceil(len(raw_ticks) / 3))]
    else:
        ticks = MaxNLocator(nbins=2).tick_values(low, high).tolist()
        low, high = ticks[0], ticks[-1]
    labels = [f"{tick / 1000:g}k" if tick >= 1000 else f"{tick:g}" for tick in ticks]
    expression_axis = axis(mean_left, mean_width, (low, high), ticks, labels, log=log_scale)
    heading(mean_left + mean_width / 2, "Mean\nexpression", "TPM (log scale)" if log_scale else "TPM")
    if groups:
        expression_axis.boxplot(
            groups, positions=[centers[name] for name in group_names], orientation="horizontal",
            widths=5.2, patch_artist=True, manage_ticks=False, showfliers=False,
            boxprops={"facecolor": "white", "edgecolor": "#527884", "linewidth": .65},
            whiskerprops={"color": "#527884"}, capprops={"color": "#527884"}, medianprops={"color": "#527884"},
        )
    layout_labels = {"paired": "Paired-end", "single": "Single-end", "unknown": "Unknown"}
    fig.legend(
        handles=[Patch(label=layout_labels[name], color=LAYOUT_COLORS[name]) for name in layout_names],
        loc="upper center", bbox_to_anchor=((pie_left + pie_width / 2) / width, (bottom - 13) / height),
        frameon=False, fontsize=6, labelcolor=MUTED, handlelength=1, labelspacing=.25,
        borderaxespad=0, borderpad=0,
    )
    original_count = sum(row["scientific_name"] in centers for row in original)
    text(12, 13.3, f"Mean TPM uses transcripts with TPM > 0 in {original_count:,} original runs.")
    text(12, 4.6, f"Other columns use {len(retained):,} retained runs. Unavailable metrics are omitted.")
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--merge-dir", type=Path, help="amalgkit_merge directory with paired metadata and TPM files")
    source.add_argument("--summary", type=Path, help="Reuse a TSV produced with --write-summary")
    parser.add_argument("--retained-runs", type=Path, help="Optional TSV with a run column for the other panels")
    parser.add_argument("--write-summary", type=Path, help="Save the per-RUN summary for future plots")
    parser.add_argument("--output", type=Path, required=True, help="Output PDF, SVG or PNG")
    args = parser.parse_args()
    if args.output.suffix.lower() not in (".pdf", ".svg", ".png"):
        parser.error("--output must end in .pdf, .svg or .png")
    try:
        retained = (
            set(by_run(read_table(args.retained_runs, ("run",)), args.retained_runs)) if args.retained_runs else None
        )
        if retained == set():
            raise ValueError("--retained-runs is empty")
        original = (
            read_table(args.summary, SUMMARY_COLUMNS) if args.summary
            else collect_merge_qc(args.merge_dir, retained)
        )
        records = by_run(original, args.summary or args.merge_dir)
        if retained is not None and (missing := retained - records.keys()):
            raise ValueError(f"Retained RUNs not found: {', '.join(sorted(missing))}")
        selected = original if retained is None else [row for row in original if row["run"] in retained]
        if not selected:
            raise ValueError("No quantified RUNs to plot")
        fig = draw_qc(original, selected)
        if args.write_summary:
            write_summary(args.write_summary, original)
        from matplotlib import rc_context

        with rc_context({"pdf.fonttype": 42, "svg.fonttype": "none"}):
            fig.savefig(args.output, dpi=240)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
