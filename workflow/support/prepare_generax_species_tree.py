#!/usr/bin/env python3
"""Serialize a rooted species tree for GeneRax without changing its source.

GeneRax reads plain Newick but rejects NHX and leading rooting declarations.
The root stays encoded by the topology of this temporary tool input; the
annotated source remains the input tracked by artifact provenance.
"""

import argparse
from pathlib import Path

from nwkit.tool_trees import prepare_species_tree as prepare_species_tree


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    prepare_species_tree(args.input, args.output)


if __name__ == "__main__":
    main()
