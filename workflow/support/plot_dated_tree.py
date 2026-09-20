#!/usr/bin/env python3
"""Draw a public-unit dated tree and its supplied age intervals with NWKIT."""

import argparse
import math
import subprocess
import sys
import tempfile
from pathlib import Path


def _draw_command(infile, outfile, height, time_credible_intervals):
    return [
        "nwkit", "draw", "--infile", str(infile), "--outfile", str(outfile),
        "--species-overlap-node-plot", "no", "--support-labels", "no",
        "--time-constraints", "no", "--time-credible-intervals", time_credible_intervals,
        "--scale-bar", "auto", "--branch-length-unit", "Ma",
        "--figure-width", "7.0", "--figure-height", str(height),
        "--tip-label-wrap", "auto", "--tip-spacing", "label-aware",
    ]


def _run_checked(command):
    return subprocess.run(command, check=True, text=True, capture_output=True)


def _replay_failure(result):
    if result.stdout:
        sys.stdout.write(result.stdout)
    if result.stderr:
        sys.stderr.write(result.stderr)


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
    draw_command = _draw_command(infile, outfile, height, "auto")
    try:
        _run_checked(draw_command)
        return
    except subprocess.CalledProcessError as error:
        combined = f"{error.stdout or ''}\n{error.stderr or ''}"
        if "Dated-tree branch lengths are not ultrametric at an internal node." not in combined:
            _replay_failure(error)
            raise

    # MCMCtree writes rounded posterior means. A tiny rounding discrepancy can
    # make CI-aware parsing reject an otherwise valid dated tree. Keep the
    # point estimates, drop only the intervals, and make the limitation clear.
    print(
        "Warning: rounded dated-tree branch lengths are slightly non-ultrametric; "
        "rendering the point estimate without credible intervals.",
        file=sys.stderr,
    )
    with tempfile.TemporaryDirectory(prefix="gg-dated-tree-") as temporary:
        fallback = Path(temporary) / "dated.nwk"
        _run_checked([
            "nwkit", "convert", "--infile", str(infile), "--outfile", str(fallback),
            "--to", "newick", "--properties", "drop",
        ])
        _run_checked(_draw_command(fallback, outfile, height, "no"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--infile", required=True, type=Path)
    parser.add_argument("--outfile", required=True, type=Path)
    args = parser.parse_args()
    draw_dated_tree(args.infile, args.outfile)


if __name__ == "__main__":
    main()
