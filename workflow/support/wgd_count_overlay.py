"""Validate and display a completed count scan alongside gene-tree WGD evidence."""

import json
from math import isclose, isfinite
from types import SimpleNamespace

try:
    from wgd_evidence import number, read_table
except ImportError:
    from .wgd_evidence import number, read_table

COUNT_OVERLAY_FIELDS = ("species_event_id", "branch_id", "num_families", "num_bootstrap", "likelihood_ratio_statistic", "p_value",
                        "burst_minus_event_aic", "count_support", "overlay_status")
STYLE = {
    "count_supported": ("D", "#286447", "Conditional count support"),
    "count_bound": ("x", "#9E5131", "Nuisance bound"),
    "ssd_preferred": ("v", "#9E5131", "Branch burst preferred"),
    "background_compatible": (".", "#69717B", "Count background compatible"),
    "not_calibrated": ("+", "#69717B", "Count scan uncalibrated"),
}


def load_count_overlay(tree, table, model_file):
    from nwkit.wgd_count import _count_tree
    from nwkit.wgd_count_fit import bootstrap_exceedances
    from nwkit.wgd_tree import model_mapping

    design, _, _ = _count_tree(tree)
    model = json.loads(model_file.read_text(encoding="utf-8"))
    if model.get("schema_version") != 1 or model.get("method") != "native-linear-DL-retained-WGM-count-v1":
        raise ValueError("Unrecognized native WGD count overlay model")
    mapping, _ = model_mapping(model, design)
    expected = {}
    for item in model["candidates"]:
        event, burst = item["event"], item["branch_burst"]
        node = mapping[event["event"]["node"]]
        if node == 0 or design.lengths[node] <= 0 or design.clade_ids[node] in expected:
            raise ValueError("Count overlay requires unique non-root model candidates")
        expected[design.clade_ids[node]] = (node, event, burst)
    if not expected:
        raise ValueError("Count overlay has no candidates")
    p_values = model.get("bootstrap_max_statistics", [])
    if any(type(v) not in (float, int) or number(v) is None or v < 0 for v in p_values):
        raise ValueError("Invalid count overlay bootstrap statistics")
    alpha = float(model["alpha"])
    if not 0 < alpha < 1:
        raise ValueError("Invalid count overlay alpha")
    calibrated = bool(p_values)
    if model["calibration"] != ("plugin-parametric-bootstrap-search-maximum" if calibrated else "not-run"):
        raise ValueError("Count overlay calibration does not match its bootstrap draws")
    n = int(model["num_families"])
    if n < 2:
        raise ValueError("Count overlay needs at least two families")
    rows = []
    fit = model["background"]
    for row in read_table(table):
        event_id = row.get("species_event_id")
        if event_id not in expected:
            raise ValueError("Count overlay table contains a duplicate or unexpected event")
        node, event, burst = expected.pop(event_id)
        lr = number(row.get("likelihood_ratio_statistic"))
        aic = float(row["burst_minus_event_aic"])
        diff = max(0., 2 * (event["log_likelihood"] - fit["log_likelihood"]))
        aic_diff = burst["aic"] - event["aic"]
        probability = number(row.get("p_value"))
        pv = ((1 + int(bootstrap_exceedances(p_values, [SimpleNamespace(improvement=diff)])[0])) / (1 + len(p_values))
              if calibrated else None)
        bounds = (event["nuisance_bound_reached"], fit["nuisance_bound_reached"], burst["nuisance_bound_reached"])
        status = ("not_calibrated" if pv is None else "background_compatible" if pv > alpha else
                  "branch_burst_preferred" if aic_diff <= 0 else "count_supported_conditional")
        if (lr is None or not isfinite(aic) or not isclose(lr, diff, rel_tol=1e-10, abs_tol=1e-11)
                or not isclose(aic, aic_diff, rel_tol=1e-10, abs_tol=1e-11)
                or int(row["num_families"]) != n or int(row["num_bootstrap"]) != len(p_values)
                or row["p_value_method"] != model["calibration"] or row["count_support"] != status
                or (pv is None and row["p_value"] != "NA") or (pv is not None and
                    (probability is None or not isclose(probability, pv, rel_tol=1e-10)))
                or any(row[key] != str(bound) for key, bound in zip(
                    ("nuisance_bound_reached", "background_nuisance_bound_reached", "branch_burst_nuisance_bound_reached"),
                    bounds, strict=True))):
            raise ValueError("Count overlay TSV and JSON do not describe the same fit")
        display = ("count_bound" if any(bounds) else "count_supported" if status == "count_supported_conditional" else
                   "ssd_preferred" if status == "branch_burst_preferred" else status)
        values = (event_id, design.branch_ids[node], n, len(p_values), lr, probability, aic, status, display)
        rows.append(dict(zip(COUNT_OVERLAY_FIELDS, values, strict=True)))
    if expected:
        raise ValueError("Count overlay TSV does not cover the complete model candidate search")
    return sorted(rows, key=lambda row: row["branch_id"]), alpha


def draw_count_overlay(ax, positions, rows, alpha):
    from matplotlib.lines import Line2D

    calibrated = any(row["p_value"] is not None for row in rows)
    axis_name = "Count model: search-wide P" if calibrated else "Count model: likelihood ratio"
    n = rows[0]["num_families"]
    draws = rows[0]["num_bootstrap"]
    detail = f"n={n:,} families, {draws:,} draws" if calibrated else f"n={n:,} families"
    ax.set_xlabel(f"{axis_name}\n({detail})", fontsize=9)
    ax.set_axisbelow(True)
    ax.grid(axis="x", color="#E5E8EA", linewidth=.7)
    ax.tick_params(axis="both", labelsize=8, left=False, labelleft=False)
    for side in ("top", "left", "right"):
        ax.spines[side].set_visible(False)
    if calibrated:
        ax.axvline(alpha, color="#A3ABB2", ls=":", lw=.85)
        ax.set_xlim(-.04, 1.04)
    for row in rows:
        if row["species_event_id"] not in positions:
            continue
        x = row["p_value"] if calibrated else row["likelihood_ratio_statistic"]
        marker, color, _ = STYLE[row["overlay_status"]]
        ax.plot(x, positions[row["species_event_id"]], marker=marker, color=color, ls="none", ms=4.2)
    states = {row["overlay_status"] for row in rows}
    handles = [Line2D([], [], marker=marker, color=color, ls="none", label=label)
               for status, (marker, color, label) in STYLE.items() if status in states]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=7, frameon=False)
