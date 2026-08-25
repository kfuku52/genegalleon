#!/usr/bin/env python3
"""Write an unrooted, annotation-free GeneRax support-mapping target."""

import argparse
from pathlib import Path

import ete4


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, metavar="PATH")
    parser.add_argument("--output", required=True, metavar="PATH")
    return parser


def prepare_target(input_path, output_path):
    text = Path(input_path).read_text(encoding="utf-8").strip()
    if not text:
        raise ValueError(f"GeneRax tree is empty: {input_path}")
    tree = ete4.PhyloTree(text, parser=1)
    leaf_names = [str(leaf.name) for leaf in tree.leaves()]
    if len(leaf_names) < 4:
        raise ValueError("An unrooted UFBoot target requires at least four tips.")
    if len(leaf_names) != len(set(leaf_names)):
        raise ValueError("GeneRax UFBoot target contains duplicate leaf labels.")

    for node in tree.traverse():
        if node.is_leaf:
            continue
        node.name = ""
        node.props.pop("support", None)

    if len(tree.get_children()) == 2:
        tree.unroot()
    if len(tree.get_children()) < 3:
        raise ValueError(
            "Failed to suppress the GeneRax root for unrooted support mapping."
        )

    tree.write(outfile=str(output_path), parser=1)
    roundtrip = ete4.PhyloTree(
        Path(output_path).read_text(encoding="utf-8").strip(), parser=1
    )
    roundtrip_names = [str(leaf.name) for leaf in roundtrip.leaves()]
    if set(roundtrip_names) != set(leaf_names) or len(roundtrip_names) != len(leaf_names):
        raise ValueError("Unrooting changed the GeneRax target tip set.")
    if len(roundtrip.get_children()) < 3:
        raise ValueError("Prepared GeneRax support target is still rooted.")


def main():
    args = build_arg_parser().parse_args()
    prepare_target(args.input, args.output)


if __name__ == "__main__":
    main()
