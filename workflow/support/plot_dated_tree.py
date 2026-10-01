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
    parser.add_argument(
        "--busco-summary", type=Path, help="Annotation summary TSV with BUSCO counts and lineage metadata."
    )
    parser.add_argument(
        "--busco-results",
        type=Path,
        help="BUSCO full/short-result directory for legacy summaries without lineage metadata.",
    )
    parser.add_argument("--busco-prefix", default="busco_cds")
    parser.add_argument(
        "--geological-background",
        choices=["none", "period"],
        help="Use the publication renderer with an optional ICS geological background.",
    )
    parser.add_argument("--figure-width", type=float)
    parser.add_argument("--figure-height", type=float)
    parser.add_argument("--font-family")
    parser.add_argument("--font-size", type=float)
    parser.add_argument("--tip-order", type=Path, help="TSV with species_id in top-to-bottom order.")
    parser.add_argument("--tip-annotations", type=Path, help="TSV: species_id, colour, font_weight.")
    parser.add_argument("--node-ages", choices=["none", "root", "all"], default="none")
    parser.add_argument(
        "--age-clades", type=Path, help="TSV: descendant_species, comma-separated exact clades to label."
    )
    parser.add_argument("--layout-report", type=Path)
    args = parser.parse_args()
    if args.busco_results is not None and args.busco_summary is None:
        parser.error("--busco-results requires --busco-summary")
    if (
        args.geological_background is not None
        or args.busco_summary is not None
        or args.node_ages != "none"
        or args.age_clades is not None
        or any(
            value is not None
            for value in (
                args.figure_width,
                args.figure_height,
                args.font_family,
                args.font_size,
                args.tip_order,
                args.tip_annotations,
                args.layout_report,
            )
        )
    ):
        from dated_tree_presentation import render_dated_tree

        render_dated_tree(
            args.infile,
            args.outfile,
            busco_summary=args.busco_summary,
            busco_results=args.busco_results,
            busco_prefix=args.busco_prefix,
            geological_background=args.geological_background or "period",
            figure_width=args.figure_width if args.figure_width is not None else 7.2,
            figure_height=args.figure_height,
            font_family=args.font_family if args.font_family is not None else "Helvetica",
            font_size=args.font_size if args.font_size is not None else 8,
            tip_order=args.tip_order,
            tip_annotations=args.tip_annotations,
            node_ages=args.node_ages,
            age_clades=args.age_clades,
            layout_report=args.layout_report,
        )
    else:
        draw_dated_tree(args.infile, args.outfile)


if __name__ == "__main__":
    main()
