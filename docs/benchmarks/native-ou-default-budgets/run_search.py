"""Time a complete AICc/convergence search with the GeneGalleon fit settings."""

import argparse
import csv
import json
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from nwkit.shift_native_fit import NativeFitOptions
from nwkit.shift_native_model import ShiftData
from nwkit.shift_native_search import NativeLayoutEvaluator
from nwkit.shift_native_selection import NativeSearchRunner
from nwkit.util import read_tree

parser = argparse.ArgumentParser()
parser.add_argument("directory", type=Path)
parser.add_argument("output", type=Path)
parser.add_argument("--candidate-pool", type=int, required=True)
parser.add_argument("--refit-budget", type=int, required=True)
parser.add_argument("--screening-budget", type=int, required=True)
args = parser.parse_args()
rows = list(csv.DictReader((args.directory / "traits.tsv").open(), delimiter="\t"))
columns = [name for name in rows[0] if name != "taxon"]
lookup = {r["taxon"]: [float(r[c]) for c in columns] for r in rows}
tree = read_tree(str(args.directory / "tree.nwk"), 1, True, quiet=True)
data = ShiftData.build(tree, np.array([lookup[name] for name in tree.leaf_names()]), columns)
configuration = dict(
    max_shifts="auto",
    convergence=True,
    criterion="AICc",
    search_strategy="auto",
    exhaustive_max_configurations=5000,
    candidate_pool=args.candidate_pool,
    refit_budget=args.refit_budget,
    screening_budget=args.screening_budget,
    beam_width=2,
    lasso_iterations=150,
    search_memory_mb=512,
)
start = time.perf_counter()
original = NativeLayoutEvaluator.evaluate


def evaluate(self, layout):
    value = original(self, layout)
    if len(self.records) % 10 == 0:
        print(
            json.dumps(
                {
                    "elapsed": time.perf_counter() - start,
                    "refits": len(self.records),
                    "last_shifts": len(layout.shifts),
                    "last_groups": len(layout.groups),
                }
            ),
            flush=True,
        )
    return value


NativeLayoutEvaluator.evaluate = evaluate
result = NativeSearchRunner(
    data, SimpleNamespace(**configuration), {"options": NativeFitOptions(estimate_measurement_error=True)}
)(data)
selected = result.best_information
if selected is None:
    raise ValueError("No finite AICc candidate")
nodes = list(tree.traverse("levelorder"))
predicted = data.centers + data.scales * np.column_stack([f.predicted for f in selected["fits"]])
record = {
    "status": "complete",
    "configuration": configuration,
    "fit_configuration": {
        "estimate_measurement_error": True,
        "estimate_alpha": True,
        "estimate_process_variance": True,
    },
    "elapsed_search_seconds": time.perf_counter() - start,
    "metadata": result.metadata,
    "information_criterion": selected["information_criterion"],
    "log_likelihood": selected["log_likelihood"],
    "shift_clades": [sorted(nodes[i].leaf_names()) for i in selected["layout"].shifts],
    "groups": [list(g) for g in selected["layout"].groups],
    "predicted": dict(zip(data.tree.leaf_names, predicted.tolist(), strict=True)),
    "candidates": result.records,
    "traits": selected["traits"],
}
args.output.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
print(
    json.dumps(
        {
            "complete": True,
            "seconds": record["elapsed_search_seconds"],
            "selected_shifts": len(record["shift_clades"]),
            "largest_fitted_shifts": max(len(r["shift_branch_ids"]) for r in result.records),
            "refits": len(result.records),
            "quick_evaluations": result.metadata.get("quick_evaluations"),
        }
    ),
    flush=True,
)
