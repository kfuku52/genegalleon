#!/usr/bin/env python3
"""Draw a public-unit dated tree and its supplied age intervals with NWKIT."""

import argparse
import math
import subprocess
from pathlib import Path


def draw_dated_tree(infile, outfile):
    from nwkit.file_paths import validate_outputs_do_not_replace_inputs
    from nwkit.util import read_tree

    validate_outputs_do_not_replace_inputs([("tree", infile)], [("plot", outfile)])
    tree = read_tree(str(infile), "auto", True)
    if any(node.dist is None or not math.isfinite(node.dist) or node.dist < 0
           for node in tree.traverse() if node is not tree):
        raise ValueError("Dated-tree branch lengths must be finite and nonnegative.")
    height = max(3.0, 1.2 + 0.24 * len(tree))
    Path(outfile).parent.mkdir(parents=True, exist_ok=True)
    subprocess.run([
        "nwkit", "draw", "--infile", str(infile), "--outfile", str(outfile),
        "--species-overlap-node-plot", "no", "--support-labels", "no",
        "--time-constraints", "no", "--time-credible-intervals", "auto",
        "--scale-bar", "auto", "--branch-length-unit", "Ma",
        "--figure-width", "7.0", "--figure-height", str(height),
        "--tip-label-wrap", "auto", "--tip-spacing", "label-aware",
    ], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--infile", required=True, type=Path)
    parser.add_argument("--outfile", required=True, type=Path)
    args = parser.parse_args()
    draw_dated_tree(args.infile, args.outfile)


if __name__ == "__main__":
    main()
