"""Independent OU recursion for budget-tuning workloads, not an adoption study."""

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def generate(traits):
    directory = ROOT / "data" / f"1000tips-{traits}traits"
    directory.mkdir(parents=True, exist_ok=False)
    nodes = []

    def build(start, count):
        node = {"name": f"t{start}" if count == 1 else "", "children": []}
        if count > 1:
            node["children"] = [build(start, count // 2), build(start + count // 2, count - count // 2)]
        node["height"] = max((c["height"] for c in node["children"]), default=-1) + 1
        node["tips"] = sorted(sum((c["tips"] for c in node["children"]), [])) if count > 1 else [node["name"]]
        nodes.append(node)
        return node

    tree = build(0, 1000)
    height = tree["height"]

    def newick(node):
        return (
            node["name"]
            if not node["children"]
            else "("
            + ",".join(
                newick(c) + ":" + format((node["height"] - c["height"]) / height, ".17g") for c in node["children"]
            )
            + ")"
        )

    rng = np.random.default_rng(32101)
    candidates = [n for n in nodes if 5 <= len(n["tips"]) <= 10]
    chosen = [candidates[i] for i in rng.choice(len(candidates), 100, replace=False)]
    assert len(set(sum((n["tips"] for n in chosen), []))) == sum(len(n["tips"]) for n in chosen)
    # Five repeated nonbaseline optimum regimes; locations shared across traits.
    group_ids = rng.integers(1, 6, size=100)
    optima = np.array([[0, 6, -6, 9, -9, 12], [0, -6, 9, 6, -12, -9]], dtype=float)[:traits]
    assigned = {id(n): int(g) for n, g in zip(chosen, group_ids, strict=True)}
    noise = rng.normal(size=(len(nodes), 2))[:, :traits]
    tip_noise = rng.normal(size=(1000, 2))[:, :traits]
    node_index = {id(n): i for i, n in enumerate(nodes)}
    alphas = np.array([3.0, 1.5])[:traits]
    rows = []

    def simulate(node, state, mean, regime=0, parent_height=None):
        regime = assigned.get(id(node), regime)
        if parent_height is not None:
            length = (parent_height - node["height"]) / height
            decay = np.exp(-alphas * length)
            variance = -np.expm1(-2 * alphas * length) / -np.expm1(-2 * alphas)
            state = decay * state + np.sqrt(variance) * noise[node_index[id(node)]]
            mean = decay * mean + (1 - decay) * optima[:, regime]
        for child in node["children"]:
            simulate(child, state, mean, regime, node["height"])
        if not node["children"]:
            observed = state + mean + math.sqrt(0.2) * tip_noise[int(node["name"][1:])]
            rows.append((node["name"], observed, mean))

    simulate(tree, np.zeros(traits), np.zeros(traits))
    (directory / "tree.nwk").write_text(newick(tree) + ";\n")
    with (directory / "traits.tsv").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxon", *[f"x{i + 1}" for i in range(traits)]])
        writer.writerows([[name, *values] for name, values, _ in rows])
    truth = {
        "tips": 1000,
        "traits": traits,
        "seed": 32101,
        "true_shifts": 100,
        "true_nonbaseline_regimes": 5,
        "alpha_height": alphas.tolist(),
        "process_tip_variance": 1.0,
        "measurement_variance": 0.2,
        "shift_clades": [n["tips"] for n in chosen],
        "shift_groups": group_ids.tolist(),
        "tip_mean": {name: mean.tolist() for name, _, mean in rows},
        "sha256": {
            name: hashlib.sha256((directory / name).read_bytes()).hexdigest() for name in ["tree.nwk", "traits.tsv"]
        },
    }
    (directory / "truth.json").write_text(json.dumps(truth, indent=2) + "\n")


if __name__ == "__main__":
    for count in (1, 2):
        generate(count)
