#!/usr/bin/env python3
"""Select MAD/midpoint roots among NWKIT optimal reconciliation candidates.

Priority: compatible MAD, compatible midpoint, first canonical optimal root.
An empty or invalid candidate collection is an error, not a fallback.
"""

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from input_generation_array_state import FreshDigestBatch
from nwkit.root_candidates import _topology_splits as _topology_splits


def select_root(args):
    from nwkit.root_candidates import select_eligible_root
    inputs = FreshDigestBatch()
    inputs.read([args.in_tree, args.candidate_trees])
    return select_eligible_root(args, before_publish=inputs.check, priority=("mad", "midpoint", "candidate"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-tree", required=True, type=Path)
    parser.add_argument("--candidate-trees", required=True, type=Path)
    parser.add_argument("--out-tree", required=True, type=Path)
    parser.add_argument("--comparison-table", required=True, type=Path)
    parser.add_argument("--comparison-plot", required=True, type=Path)
    parser.add_argument("--species-parser", default="extended-taxonomic")
    select_root(parser.parse_args())


if __name__ == "__main__":
    main()
