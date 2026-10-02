"""Configurable, source-backed comparison figures from signed statistics.

Display filters never change the original Benjamini–Hochberg test family.
Local subgenome labels are not given persistent shapes or parental identities.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import json
import math
from pathlib import Path

METRICS = ("retention_difference", "expression_log2_ratio", "expression_detection_difference")
DEFAULTS = {"font_family": "DejaVu Sans", "font_size": 8, "font_files": [],
            "width_pt": None, "height_pt": None, "point_colour": "black", "error_bar_colour": "black",
            "individual_width_pt": None, "individual_height_pt": None,
            "significance_colour": "#D55E00", "point_size": 3.4, "error_bar_width": 0.8, "capsize": 1.5,
            "dpi": 180, "formats": ["png", "svg", "pdf"], "significance_threshold": 0.05,
            "show_significance": True, "species_order": [], "analysis_ids": [],
            "panel_labels": True,
            "reference": "", "pair_set": "", "tissue": "", "metrics": list(METRICS),
            "contrast_anchors": {}, "group_labels": {}, "x_limits": {}, "absolute_x_limits": {}}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_config(config):
    from matplotlib.colors import is_color_like

    unknown = set(config) - set(DEFAULTS)
    if unknown:
        raise ValueError(f"Unknown subgenome plot settings: {sorted(unknown)}")
    result = {**DEFAULTS, **config}
    for key in ("font_size", "point_size", "error_bar_width", "dpi"):
        if not isinstance(result[key], (int, float)) or not math.isfinite(result[key]) or result[key] <= 0:
            raise ValueError(f"{key} must be positive and finite")
    for key in ("width_pt", "height_pt", "individual_width_pt", "individual_height_pt"):
        if result[key] is not None and (not isinstance(result[key], (int, float))
                                       or not math.isfinite(result[key]) or result[key] <= 0):
            raise ValueError(f"{key} must be positive and finite")
    if not isinstance(result["capsize"], (int, float)) or not math.isfinite(result["capsize"]) or result["capsize"] < 0:
        raise ValueError("capsize must be finite and nonnegative")
    if not isinstance(result["significance_threshold"], (int, float)) or not 0 < result["significance_threshold"] < 1:
        raise ValueError("significance_threshold must be between zero and one")
    for key in ("show_significance", "panel_labels"):
        if not isinstance(result[key], bool):
            raise ValueError(f"{key} must be boolean")
    for key in ("font_family", "reference", "pair_set", "tissue"):
        if not isinstance(result[key], str) or (key == "font_family" and not result[key].strip()):
            raise ValueError(f"{key} must be a string")
    for key in ("species_order", "analysis_ids", "font_files", "metrics", "formats"):
        items = result[key]
        if not isinstance(items, list) or any(not isinstance(x, str) or not x for x in items) or len(set(items)) != len(items):
            raise ValueError(f"{key} must be a list of unique nonempty strings")
    if not result["metrics"] or not set(result["metrics"]) <= set(METRICS):
        raise ValueError("Unsupported or empty metric selection")
    if not result["formats"] or not set(result["formats"]) <= {"png", "svg", "pdf"}:
        raise ValueError("Formats must be png, svg or pdf")
    for key in ("point_colour", "error_bar_colour", "significance_colour"):
        if not is_color_like(result[key]):
            raise ValueError(f"Invalid {key}")
    for key in ("contrast_anchors", "group_labels", "x_limits", "absolute_x_limits"):
        if not isinstance(result[key], dict):
            raise ValueError(f"{key} must be an object")
    if any(not isinstance(v, str) or not v for v in result["contrast_anchors"].values()):
        raise ValueError("contrast_anchors maps species to a subgenome label")
    if any(not isinstance(labels, dict) or any(not isinstance(v, str) for v in labels.values())
           for labels in result["group_labels"].values()):
        raise ValueError("group_labels maps species to group/label objects")
    for key in ("x_limits", "absolute_x_limits"):
        for metric, limits in result[key].items():
            if metric not in METRICS or not isinstance(limits, list) or len(limits) != 2 or any(
                    not isinstance(x, (int, float)) or not math.isfinite(x) for x in limits) or limits[0] >= limits[1]:
                raise ValueError(f"{key} requires metric: [finite lower, finite upper]")
    return result


def load_config(path=None):
    inputs = {}
    config = {}
    if path:
        path = Path(path).resolve()
        config = json.loads(path.read_text())
        if not isinstance(config, dict):
            raise ValueError("Plot config must be a JSON object")
        inputs[str(path)] = sha256(path)
        config = validate_config(config)
        config["font_files"] = [str((path.parent / font).resolve()) for font in config["font_files"]]
        inputs.update({font: sha256(font) for font in config["font_files"]})
    return validate_config(config), inputs


def display_interval(row, absolute):
    scale = 100 if row["metric"] != "expression_log2_ratio" else 1
    effect = row["effect"] * scale
    low, high = row.get("ci_low"), row.get("ci_high")
    if (low is None) != (high is None):
        raise ValueError("Both confidence bounds are required together")
    if low is not None:
        low, high = low * scale, high * scale
        if not all(math.isfinite(v) for v in (effect, low, high)) or not low <= effect <= high:
            raise ValueError("Invalid signed confidence interval")
    elif not math.isfinite(effect):
        raise ValueError("Nonfinite effect")
    if absolute:
        effect = abs(effect)
        if low is not None:
            low, high = (0 if low <= 0 <= high else min(abs(low), abs(high))), max(abs(low), abs(high))
    return effect, low, high


def select_analyses(analyses, config, apply_filters=True):
    if len({a["analysis_id"] for a in analyses}) != len(analyses):
        raise ValueError("Comparison analysis IDs must be unique")
    selected = []
    for analysis in analyses:
        if apply_filters and any(config[key] and analysis.get(key, "") != config[key] for key in ("reference", "pair_set")):
            continue
        if apply_filters and config["analysis_ids"] and analysis["analysis_id"] not in config["analysis_ids"]:
            continue
        anchor = config["contrast_anchors"].get(analysis["species"], "") if apply_filters else ""
        rows = [r for r in analysis["statistics"] if r["metric"] in config["metrics"] and r["effect"] is not None
                and (not apply_filters or not config["tissue"] or not r["tissue"] or r["tissue"] == config["tissue"])
                and (not anchor or anchor in (r["subgenome_a"], r["subgenome_b"]))]
        identities = [(r["metric"], r["group_id"], r["subgenome_a"], r["subgenome_b"], r["tissue"]) for r in rows]
        if len(set(identities)) != len(identities):
            raise ValueError("Repeated comparison identity in a statistics table")
        if rows:
            selected.append({**analysis, "statistics": rows})
    if not selected:
        raise ValueError("No estimable contrasts match the display filters")
    if apply_filters and config["analysis_ids"] and set(config["analysis_ids"]) - {a["analysis_id"] for a in selected}:
        raise ValueError("Requested analysis_ids are missing or have no selected contrasts")
    if apply_filters and config["species_order"]:
        if set(config["species_order"]) != {a["species"] for a in selected}:
            raise ValueError("species_order must contain every selected species exactly once")
        selected.sort(key=lambda a: config["species_order"].index(a["species"]))
    return selected


def comparison_plot(analyses, output, *, config=None, absolute=True, stem=None, apply_filters=True, input_hashes=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib.legend_handler import HandlerBase
    from matplotlib.text import Text

    config = validate_config(config or {})
    analyses = select_analyses(analyses, config, apply_filters)
    metrics = [m for m in config["metrics"] if any(r["metric"] == m for a in analyses for r in a["statistics"])]
    for font in config["font_files"]:
        font_manager.fontManager.addfont(font)
    selector = font_manager.fontManager
    if config["font_files"]:
        # A preinstalled face of the same family must not supersede supplied fonts.
        selector = copy.copy(font_manager.fontManager)
        supplied = {str(Path(font).resolve()) for font in config["font_files"]}
        selector.ttflist = [entry for entry in selector.ttflist if str(Path(entry.fname).resolve()) in supplied]
    selector.findfont(font_manager.FontProperties(family=config["font_family"]), fallback_to_default=False)
    size = config["font_size"]
    ncols = len(analyses)
    max_groups = max(len({(r["group_id"], r["tissue"]) for r in a["statistics"] if r["metric"] == m})
                     for a in analyses for m in metrics)
    width = config["width_pt"] or 300 * ncols
    height = config["height_pt"] or len(metrics) * max(140, max_groups * (size + 7) + 65) + 25
    style = {"font.family": config["font_family"], "font.size": size, "axes.labelsize": size,
             "axes.titlesize": size, "xtick.labelsize": size, "ytick.labelsize": size,
             "legend.fontsize": size, "pdf.fonttype": 42, "svg.fonttype": "none"}
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    stem = stem or ("comparison_absolute" if absolute else "comparison")
    plotted, limits_by_metric = [], {}
    with plt.rc_context(style):
        fig, axes = plt.subplots(len(metrics), ncols, figsize=(width / 72, height / 72),
                                 squeeze=False, layout="constrained")
        try:
            for i, metric in enumerate(metrics):
                bounds = [v for a in analyses for r in a["statistics"] if r["metric"] == metric
                          for v in display_interval(r, absolute) if v is not None]
                lower, upper = min(0, min(bounds)), max(0, max(bounds))
                span = upper - lower or 1
                limits = config["absolute_x_limits" if absolute else "x_limits"].get(
                    metric, [lower - span * .04, upper + span * .13])
                if limits[0] > lower or limits[1] < upper:
                    raise ValueError(f"x_limits clips {metric} estimates, intervals or zero")
                limits_by_metric[metric] = limits
                for j, analysis in enumerate(analyses):
                    axis = axes[i, j]
                    selected = sorted((r for r in analysis["statistics"] if r["metric"] == metric),
                                      key=lambda r: (r["group_id"], r["tissue"], r["subgenome_a"], r["subgenome_b"]))
                    buckets = {}
                    for row in selected:
                        buckets.setdefault((row["group_id"], row["tissue"]), []).append(row)
                    for y, ((group, tissue), rows) in enumerate(buckets.items()):
                        for k, row in enumerate(rows):
                            position = y + ((k / (len(rows) - 1) - .5) * .6 if len(rows) > 1 else 0)
                            effect, low, high = display_interval(row, absolute)
                            axis.plot(effect, position, "o", color=config["point_colour"], ms=config["point_size"])
                            if low is not None:
                                axis.errorbar(effect, position, xerr=[[effect - low], [high - effect]], fmt="none",
                                              ecolor=config["error_bar_colour"], elinewidth=config["error_bar_width"],
                                              capsize=config["capsize"])
                            significant = row.get("q_value") is not None and row["q_value"] < config["significance_threshold"]
                            if config["show_significance"] and significant:
                                axis.annotate("*", (effect, position), xytext=(3, 1), textcoords="offset points",
                                              color=config["significance_colour"], fontsize=size)
                            plotted.append({**{key: value for key, value in analysis.items() if key != "statistics"},
                                            **row, "panel": f"{i + 1}:{j + 1}", "plotted_effect": effect,
                                            "plotted_ci_low": low, "plotted_ci_high": high, "plotted_y": position,
                                            "plotted_unit": "log2_ratio" if metric == "expression_log2_ratio" else "percentage_points"})
                    labels = config["group_labels"].get(analysis["species"], {})
                    axis.set_yticks(range(len(buckets)), [(labels.get(g, g) + (" " + t if t else "")) for g, t in buckets])
                    axis.set_ylim(max(0, len(buckets) - 1) + .65, -.65)
                    axis.set_xlim(limits)
                    axis.axvline(0, color="0.6", lw=.6, zorder=0)
                    axis.spines[["top", "right"]].set_visible(False)
                    label = {"retention_difference": "Retention / syntelog\ndetection difference\n(percentage points)",
                             "expression_detection_difference": "Expression detection\ndifference (percentage points)",
                             "expression_log2_ratio": "Mean log2\nexpression ratio (A/B)"}[metric]
                    if metric == "retention_difference" and analysis.get("retention_status", "").startswith("exploratory"):
                        label = "Syntelog\ndetection difference\n(percentage points)"
                    if absolute:
                        label = "Absolute " + label[0].lower() + label[1:]
                        if metric == "expression_log2_ratio":
                            label = label.replace(" (A/B)", "")
                    axis.set_xlabel(label)
                    if config["panel_labels"]:
                        index = i * ncols + j
                        axis.text(-.04, 1.03, chr(97 + index) if index < 26 else str(index + 1),
                                  transform=axis.transAxes, va="bottom", ha="right", fontweight="bold", fontsize=size)
                    if i == 0:
                        title = analysis["species"].replace("_", " ")
                        if sum(a["species"] == analysis["species"] for a in analyses) > 1:
                            title += "\n" + analysis["analysis_id"]
                        axis.set_title(title, fontstyle="italic", fontsize=size)
                    if not selected:
                        axis.text(.5, .5, "Not estimable", transform=axis.transAxes, ha="center", fontsize=size)
            if config["show_significance"]:
                class StarHandler(HandlerBase):
                    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
                        return [Text(x=width / 2 - xdescent, y=height / 2 - ydescent, text="*", ha="center", va="center",
                                     fontfamily=config["font_family"], fontsize=size,
                                     color=config["significance_colour"], transform=trans)]
                handle = object()
                fig.legend([handle], [f"Benjamini–Hochberg q < {config['significance_threshold']:g}"],
                           handler_map={handle: StarHandler()}, loc="outside lower center", frameon=False)
            if config["font_files"]:
                for text in fig.findobj(Text):
                    properties = text.get_fontproperties().copy()
                    properties.set_file(selector.findfont(properties, fallback_to_default=False))
                    text.set_fontproperties(properties)
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            fonts, text_count = {}, 0
            # Matplotlib retains unused tick artists outside the view limits.
            # They are not drawn, so their off-canvas bounds are not figure defects.
            unused_tick_texts = set()
            for axis in axes.flat:
                for coordinate, limits in ((axis.xaxis, axis.get_xlim()), (axis.yaxis, axis.get_ylim())):
                    low, high = sorted(limits)
                    for tick in [*coordinate.get_major_ticks(), *coordinate.get_minor_ticks()]:
                        if not low - 1e-12 <= tick.get_loc() <= high + 1e-12:
                            unused_tick_texts.update((id(tick.label1), id(tick.label2)))
            for text in fig.findobj(Text):
                if not text.get_visible() or not text.get_text() or id(text) in unused_tick_texts:
                    continue
                if not math.isclose(text.get_fontsize(), size):
                    raise ValueError("Figure contains an unexpected text size")
                path = font_manager.findfont(text.get_fontproperties(), fallback_to_default=False)
                if font_manager.FontProperties(fname=path).get_name() != config["font_family"]:
                    raise ValueError("Figure resolved to a different font family")
                fonts[path] = sha256(path)
                box = text.get_window_extent(renderer)
                if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width + 1 or box.y1 > fig.bbox.height + 1:
                    raise ValueError(f"Text escapes the figure canvas: {text.get_text()!r}; increase width_pt or height_pt")
                text_count += 1
            for extension in config["formats"]:
                fig.savefig(output / f"{stem}.{extension}", dpi=config["dpi"])
        finally:
            plt.close(fig)
    fields = list(dict.fromkeys(key for row in plotted for key in row))
    with (output / f"{stem}_points.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(plotted)
    provenance = {"schema_version": 1, "config": config, "absolute": absolute,
                  "figure_size_pt": [width, height], "point_count": len(plotted), "text_count": text_count,
                  "point_unit": "group-specific subgenome contrast", "resampling_unit": "nonoverlapping block",
                  "ci": "signed 95% block bootstrap confidence set mapped under abs" if absolute else "signed 95% block bootstrap interval",
                  "tests": "signed tests; Benjamini–Hochberg families unchanged by display filters",
                  "local_labels": "identities across local groups are unresolved; all points use circles",
                  "shared_x_limits": limits_by_metric, "font_hashes": fonts, "input_hashes": input_hashes or {},
                  "renderer_sha256": sha256(__file__),
                  "matplotlib_version": matplotlib.__version__,
                  "output_hashes": {f"{stem}.{ext}": sha256(output / f"{stem}.{ext}") for ext in config["formats"]}}
    provenance["output_hashes"][f"{stem}_points.tsv"] = sha256(output / f"{stem}_points.tsv")
    (output / f"{stem}_provenance.json").write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
    return provenance


def report(manifest, output, config_path=None):
    config, inputs = load_config(config_path)
    manifest = Path(manifest).resolve()
    inputs[str(manifest)] = sha256(manifest)
    with manifest.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or not {"analysis_id", "species", "statistics_file"} <= set(reader.fieldnames):
            raise ValueError("Report manifest requires analysis_id, species, statistics_file")
        entries = list(reader)
    analyses = []
    numeric = {"effect", "ci_low", "ci_high", "p_value", "q_value", "mc_p_ci_low", "mc_p_ci_high"}
    integers = {"n_loci", "n_blocks", "n_opportunities", "retained_a", "retained_b", "null_draws",
                "bootstrap_replicates", "n_nonzero_blocks", "inference_version", "multiple_testing_n"}
    for entry in entries:
        if any(not entry.get(key, "") for key in ("analysis_id", "species", "statistics_file")):
            raise ValueError("Incomplete comparison manifest")
        path = (manifest.parent / entry["statistics_file"]).resolve()
        inputs[str(path)] = sha256(path)
        with path.open(newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if not reader.fieldnames or not {"metric", "group_id", "subgenome_a", "subgenome_b", "tissue",
                                              "effect", "ci_low", "ci_high", "q_value", "n_loci", "n_blocks"} <= set(reader.fieldnames):
                raise ValueError("Incomplete statistics table")
            rows = [{k: (float(v) if v else None) if k in numeric else (int(v) if v else None) if k in integers else v
                     for k, v in row.items()} for row in reader]
        analyses.append({**entry, "statistics_file": str(path), "statistics": rows})
    results = [comparison_plot(analyses, output, config=config, absolute=absolute, input_hashes=inputs)
               for absolute in (False, True)]
    for path, expected in inputs.items():
        if sha256(path) != expected:
            raise ValueError(f"Report input changed during rendering: {path}")
    return results
