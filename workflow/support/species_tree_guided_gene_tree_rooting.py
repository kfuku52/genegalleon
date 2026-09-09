#!/usr/bin/env python3
"""Select NWKIT MAD/midpoint roots against existing NOTUNG candidate roots.

NOTUNG inference and the selection priority are unchanged: compatible MAD,
compatible midpoint, first NOTUNG candidate, then MAD when there are no candidates.
Root compatibility means the same unrooted edge, independent of its split ratio.
"""

import argparse
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


def natural_key(path):
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", path.name)]


def select_root(args):
    from nwkit.clade_mapping import projected_root_split
    from nwkit.file_paths import validate_outputs_do_not_replace_inputs
    from nwkit.output_transaction import output_transaction, validate_output_targets
    from nwkit.util import read_tree

    if args.notung_root_dir is not None and not args.notung_root_dir.is_dir():
        raise ValueError(f"NOTUNG candidate directory does not exist or is not a directory: {args.notung_root_dir}")
    candidates = sorted(
        (path for path in args.notung_root_dir.glob("*.rooting.*")
         if path.is_file() and re.search(r"\.rooting\.\d+$", path.name)),
        key=natural_key,
    ) if args.notung_root_dir is not None else []
    outputs = [args.out_tree, args.comparison_table, args.comparison_plot]
    validate_outputs_do_not_replace_inputs(
        [("input tree", args.in_tree), *(("NOTUNG candidate", path) for path in candidates)],
        [("output", path) for path in outputs],
    )
    validate_output_targets(outputs)
    args.out_tree.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="nwkit-root-", dir=args.out_tree.parent) as temporary:
        work = Path(temporary)
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
        notung_splits = []
        for path in candidates:
            candidate = read_tree(str(path), "auto", True)
            if len(candidate) != len(leaves) or frozenset(candidate.leaf_names()) != leaves:
                raise ValueError(f"NOTUNG candidate tips are duplicated or differ from the input tree: {path}")
            split = projected_root_split(candidate, leaves)
            if split is None:
                raise ValueError(f"NOTUNG candidate does not have a resolved root edge: {path}")
            notung_splits.append(split)
        compatible = {method: split in notung_splits for method, split in native_splits.items()}
        if compatible["mad"]:
            selected, reason = roots["mad"], "mad_compatible_with_notung"
        elif compatible["midpoint"]:
            selected, reason = roots["midpoint"], "midpoint_compatible_with_notung"
        elif candidates:
            selected, reason = candidates[0], "first_notung_candidate"
        else:
            selected, reason = roots["mad"], "mad_without_notung_candidates"
        print(f"Number of rooted NOTUNG trees: {len(candidates)}", flush=True)
        for method in ("mad", "midpoint"):
            print(f"is_{method}_compatible_with_notung: {compatible[method]}", flush=True)
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
    parser.add_argument("--notung-root-dir", type=Path)
    parser.add_argument("--out-tree", required=True, type=Path)
    parser.add_argument("--comparison-table", required=True, type=Path)
    parser.add_argument("--comparison-plot", required=True, type=Path)
    parser.add_argument("--species-parser", default="taxonomic")
    select_root(parser.parse_args())


if __name__ == "__main__":
    main()
