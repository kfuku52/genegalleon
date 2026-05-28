#!/usr/bin/env python3
"""Write a Newick tree as topology with leaf names only."""

import argparse
from pathlib import Path

from ete4 import Tree


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Strip internal labels, support values, and branch lengths from a Newick tree."
    )
    parser.add_argument("--infile", required=True, help="Input Newick file.")
    parser.add_argument("--outfile", required=True, help="Output topology-only Newick file.")
    parser.add_argument(
        "--parser",
        type=int,
        default=0,
        help="ETE Newick parser format for the input tree. Default: 0.",
    )
    args = parser.parse_args()

    tree_text = Path(args.infile).read_text(encoding="utf-8").strip()
    if not tree_text:
        raise SystemExit(f"Input Newick file is empty: {args.infile}")

    tree = Tree(tree_text, parser=args.parser)
    Path(args.outfile).write_text(tree.write(parser=9) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
