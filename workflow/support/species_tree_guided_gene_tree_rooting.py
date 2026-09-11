#!/usr/bin/env python3
"""Select MAD/midpoint roots among NWKIT optimal reconciliation candidates.

Priority: compatible MAD, compatible midpoint, first canonical optimal root.
An empty or invalid candidate collection is an error, not a fallback.
"""

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path


def select_root(args):
    from nwkit.clade_mapping import canonical_split, projected_root_split
    from nwkit.file_paths import validate_outputs_do_not_replace_inputs
    from nwkit.output_transaction import output_transaction, validate_output_targets
    from nwkit.util import read_tree, read_tree_strings

    candidate_strings = read_tree_strings(str(args.candidate_trees))
    if not candidate_strings:
        raise ValueError("Reconciliation candidate tree collection is empty.")
    outputs = [args.out_tree, args.comparison_table, args.comparison_plot]
    validate_outputs_do_not_replace_inputs(
        [("input tree", args.in_tree), ("reconciliation candidates", args.candidate_trees)],
        [("output", path) for path in outputs],
    )
    validate_output_targets(outputs)
    args.out_tree.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="nwkit-root-", dir=args.out_tree.parent) as temporary:
        work = Path(temporary)
        candidates = []
        for index, text in enumerate(candidate_strings):
            path = work / f"candidate.{index}.nwk"
            path.write_text(text + "\n", encoding="utf-8")
            candidates.append(path)
        roots = {}
        for method in ("mad", "midpoint"):
            roots[method] = work / f"{method}.nwk"
            subprocess.run([
                "nwkit", "root", "--method", method, "--infile", str(args.in_tree),
                "--outfile", str(roots[method]),
            ], check=True)
        native_trees = {method: read_tree(str(path), "auto", True) for method, path in roots.items()}
        leaves = frozenset(native_trees["mad"].leaf_names())
        if any(len(tree) != len(leaves) or frozenset(tree.leaf_names()) != leaves for tree in native_trees.values()):
            raise ValueError("NWKIT root trees must have matching, unique tip labels.")
        native_splits = {method: projected_root_split(tree, leaves) for method, tree in native_trees.items()}
        if any(split is None for split in native_splits.values()):
            raise ValueError("NWKIT did not produce a resolved root edge.")
        def topology_splits(tree):
            descendants = tree.get_cached_content()
            return {canonical_split(frozenset(leaf.name for leaf in tips),
                                    leaves - frozenset(leaf.name for leaf in tips))
                    for node, tips in descendants.items() if not node.is_root}

        expected_topology = topology_splits(native_trees["mad"])
        reconciliation_splits = []
        for path in candidates:
            candidate = read_tree(str(path), "auto", True)
            if len(candidate) != len(leaves) or frozenset(candidate.leaf_names()) != leaves:
                raise ValueError(f"Reconciliation candidate tips are duplicated or differ from the input tree: {path}")
            if topology_splits(candidate) != expected_topology:
                raise ValueError(f"Reconciliation candidate has a different unrooted topology: {path}")
            split = projected_root_split(candidate, leaves)
            if split is None:
                raise ValueError(f"Reconciliation candidate does not have a resolved root edge: {path}")
            reconciliation_splits.append(split)
        compatible = {method: split in reconciliation_splits for method, split in native_splits.items()}
        if compatible["mad"]:
            selected, reason = roots["mad"], "mad_compatible_with_reconciliation"
        elif compatible["midpoint"]:
            selected, reason = roots["midpoint"], "midpoint_compatible_with_reconciliation"
        else:
            selected, reason = candidates[0], "first_reconciliation_candidate"
        print(f"Number of optimal reconciliation roots: {len(candidates)}", flush=True)
        for method in ("mad", "midpoint"):
            print(f"is_{method}_compatible_with_reconciliation: {compatible[method]}", flush=True)
        print(f"Selected root: {reason}", flush=True)
        subprocess.run([
            "nwkit", "drop", "--infile", str(selected), "--outfile", str(work / "selected.nwk"),
            "--target", "intnode", "--name", "yes", "--support", "yes",
        ], check=True)
        subprocess.run([
            "nwkit", "rootcompare", "--infile", str(args.in_tree),
            "--methods", "midpoint,mad", "--species-parser", args.species_parser,
            "--outfile", str(work / "comparison.tsv"),
            "--figure-out", str(work / "comparison.pdf"),
        ], check=True)
        with output_transaction(outputs, create_parents=True) as staged:
            for source, target in zip(("selected.nwk", "comparison.tsv", "comparison.pdf"), outputs, strict=True):
                shutil.copyfile(work / source, staged[target])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-tree", required=True, type=Path)
    parser.add_argument("--candidate-trees", required=True, type=Path)
    parser.add_argument("--out-tree", required=True, type=Path)
    parser.add_argument("--comparison-table", required=True, type=Path)
    parser.add_argument("--comparison-plot", required=True, type=Path)
    parser.add_argument("--species-parser", default="taxonomic")
    select_root(parser.parse_args())


if __name__ == "__main__":
    main()
