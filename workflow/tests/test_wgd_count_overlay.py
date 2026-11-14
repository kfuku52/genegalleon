import csv
import json

import pytest

from workflow.support.wgd_count_overlay import draw_count_overlay, load_count_overlay


@pytest.fixture
def count_bundle(tmp_path):
    from nwkit.util import read_tree
    from nwkit.wgd_count import _count_tree

    tree = read_tree("((A:1,B:1):1,C:2);", "auto", True)
    design, _, _ = _count_tree(tree)
    node = 1
    event_id = design.clade_ids[node]
    background = {"log_likelihood": -20., "nuisance_bound_reached": False}
    event = {"log_likelihood": -18., "aic": 44., "nuisance_bound_reached": False, "event": {"node": node}}
    burst = {"log_likelihood": -19., "aic": 48., "nuisance_bound_reached": False}
    method = "plugin-parametric-bootstrap-search-maximum"
    model = dict(schema_version=1, method="native-linear-DL-retained-WGM-count-v1", num_families=50, alpha=.1,
                 tip_names=design.tip_names, species_event_ids=design.clade_ids, branch_rate_groups=[0] * len(design.parents),
                 species_tree={"parents": design.parents, "lengths": design.lengths}, background=background,
                 candidates=[{"event": event, "branch_burst": burst}], calibration=method,
                 bootstrap_max_statistics=[1., 4. - 5e-10] + [2.] * 37)
    row = dict(species_event_id=event_id, num_families=50, likelihood_ratio_statistic=4, burst_minus_event_aic=4,
               num_bootstrap=39, p_value_method=method, count_support="count_supported_conditional", p_value=.05,
               nuisance_bound_reached="False", background_nuisance_bound_reached="False", branch_burst_nuisance_bound_reached="False")

    def bundle():
        table_path, model_path = tmp_path / "counts.tsv", tmp_path / "model.json"
        with table_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, row, delimiter="\t")
            writer.writeheader()
            writer.writerow(row)
        model_path.write_text(json.dumps(model))
        return tree, table_path, model_path

    return model, row, bundle


def test_overlay_reuses_bootstrap_threshold_and_prior_alpha(count_bundle):
    _, _, bundle = count_bundle
    (row,), alpha = load_count_overlay(*bundle())
    assert row["p_value"] == .05 and alpha == .1
    assert row["num_bootstrap_exceedances"] == 1
    assert row["mc_p_lower"] == pytest.approx(1 - .975 ** (1 / 39))
    assert .10 < row["mc_p_upper"] < .15
    assert row["overlay_status"] == "count_supported"


def test_overlay_rejects_mismatched_fit_and_tree(count_bundle):
    model, row, bundle = count_bundle
    row["p_value"] = .025
    with pytest.raises(ValueError, match="same fit"):
        load_count_overlay(*bundle())
    row["p_value"] = .05
    model["species_tree"]["lengths"] = tuple(v * 2 for v in model["species_tree"]["lengths"])
    with pytest.raises(ValueError, match="branch lengths"):
        load_count_overlay(*bundle())


def test_overlay_displays_negative_aic_and_distinguishes_a_nuisance_bound(count_bundle):
    model, row, bundle = count_bundle
    model["candidates"][0]["branch_burst"]["aic"] = 41.
    row.update(burst_minus_event_aic=-3., count_support="branch_burst_preferred")
    assert load_count_overlay(*bundle())[0][0]["overlay_status"] == "ssd_preferred"
    model["candidates"][0]["branch_burst"]["nuisance_bound_reached"] = True
    row["branch_burst_nuisance_bound_reached"] = "True"
    assert load_count_overlay(*bundle())[0][0]["overlay_status"] == "count_bound"


def test_count_plot_labels_significant_branches_at_their_unchanged_p_values(count_bundle):
    from matplotlib.figure import Figure

    _, _, bundle = count_bundle
    (row,), alpha = load_count_overlay(*bundle())
    other = {**row, "species_event_id": "other", "branch_id": 99, "p_value": .5,
             "overlay_status": "background_compatible", "mc_p_lower": .4, "mc_p_upper": .6}
    axis = Figure().subplots()
    draw_count_overlay(axis, {row["species_event_id"]: 2.5, "other": 2.84}, [row, other], alpha)
    annotation, = axis.texts
    assert annotation.get_text() == f"B{row['branch_id']}"
    assert annotation.xy == (.05, 2.5)
    assert axis.lines[-1].get_xdata()[0] == .5
    assert axis.collections[-1].get_segments()[0].tolist() == [[.4, 2.84], [.6, 2.84]]


@pytest.mark.parametrize(("draw", "p_value", "status", "boundary"), [
    (0., .025, "count_supported_conditional", "mc_p_lower"),
    (5., 1., "background_compatible", "mc_p_upper"),
])
def test_mc_interval_includes_zero_and_one_at_binomial_boundaries(count_bundle, draw, p_value, status, boundary):
    model, row, bundle = count_bundle
    model["bootstrap_max_statistics"] = [draw] * 39
    row.update(p_value=p_value, count_support=status)
    (result,), _ = load_count_overlay(*bundle())
    assert result[boundary] == float(draw > 0)
