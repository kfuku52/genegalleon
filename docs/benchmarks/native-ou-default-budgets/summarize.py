"""Summarize completed searches and retain failed/timed-out measurements."""

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
rows = []
for run in sorted((ROOT / "runs").iterdir()):
    measurement = run / "measurement.json"
    if not measurement.exists():
        continue
    row = {"run": run.name, **json.loads(measurement.read_text())}
    model_file = run / "model.json"
    if model_file.exists():
        model = json.loads(model_file.read_text())
        traits = len(model["traits"])
        truth = json.loads((ROOT / "data" / f"1000tips-{traits}traits" / "truth.json").read_text())
        actual = {tuple(clade) for clade in truth["shift_clades"]}
        selected = {tuple(clade) for clade in model["shift_clades"]}
        errors = [
            (value - truth["tip_mean"][tip][trait]) ** 2
            for tip, values in model["predicted"].items()
            for trait, value in enumerate(values)
        ]
        score = model["information_criterion"]
        p, n = score["parameter_count"], score["sample_size"]
        expected = -2 * model["log_likelihood"] + 2 * p + 2 * p * (p + 1) / (n - p - 1)
        assert math.isclose(expected, score["score"], rel_tol=1e-12)
        assert score["score"] == min(
            candidate["information_criterion"]["score"]
            for candidate in model["candidates"]
            if candidate["information_criterion"]["status"] == "ok"
        )
        row.update(
            traits=traits,
            configuration=model["configuration"],
            selected_shifts=len(selected),
            selected_groups=len(model["groups"]),
            true_positive=len(actual & selected),
            false_positive=len(selected - actual),
            missed_shifts=len(actual - selected),
            tip_mean_rmse=math.sqrt(sum(errors) / len(errors)),
            aicc=score["score"],
            resolved_cap=model["metadata"]["shift_limit"]["resolved"],
            largest_fitted_shifts=max(len(candidate["shift_branch_ids"]) for candidate in model["candidates"]),
            refits=len(model["candidates"]),
            quick_evaluations=model["metadata"]["quick_evaluations"],
        )
    rows.append(row)
comparisons = []
for traits in (1, 2):
    suffix = "1trait" if traits == 1 else "2traits"
    paths = [ROOT / "runs" / f"{setting}-{suffix}" / "model.json" for setting in ("baseline", "expanded")]
    if not all(path.exists() for path in paths):
        continue

    def layouts(path):
        return {
            (tuple(candidate["shift_branch_ids"]), tuple(tuple(group) for group in candidate["groups"])): candidate
            for candidate in json.loads(path.read_text())["candidates"]
        }

    baseline, expanded = [layouts(path) for path in paths]
    shared = baseline.keys() & expanded.keys()
    delta = max(abs(baseline[key]["log_likelihood"] - expanded[key]["log_likelihood"]) for key in shared)
    assert delta < 1e-9
    comparisons.append({"traits": traits, "shared_layouts": len(shared), "max_log_likelihood_difference": delta})
(ROOT / "overlap-check.json").write_text(json.dumps(comparisons, indent=2) + "\n")
(ROOT / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")
for row in rows:
    print(json.dumps(row))
