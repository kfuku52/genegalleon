#!/usr/bin/env python3
"""Serialize a rooted species tree for GeneRax without changing its source.

GeneRax reads plain Newick but rejects NHX and leading rooting declarations.
The root stays encoded by the topology of this temporary tool input; the
annotated source remains the input tracked by artifact provenance.
"""

import argparse
from pathlib import Path


def prepare_species_tree(input_path, output_path):
    from ete4.parser import newick
    from nwkit.file_paths import validate_outputs_do_not_replace_inputs
    from nwkit.output_transaction import output_transaction
    from nwkit.rooting_state import get_rooting_info
    from nwkit.util import read_tree

    validate_outputs_do_not_replace_inputs(
        [("species tree", input_path)], [("GeneRax input", output_path)],
    )
    tree = read_tree(str(input_path), "auto", True, quiet=True)
    if get_rooting_info(tree).rooted is not True:
        raise ValueError("GeneRax species tree must be rooted.")
    tips = list(tree.leaf_names())
    if len(tips) < 2 or len(tips) != len(set(tips)):
        raise ValueError("GeneRax species tree requires at least two unique tips.")
    if any(len(node.children) != 2 for node in tree.traverse() if not node.is_leaf):
        raise ValueError("GeneRax requires a rooted, bifurcating species tree.")

    # ETE's default distance formatter rounds to six significant digits.
    parser = {kind: [dict(field) for field in fields]
              for kind, fields in newick.PARSERS[1].items()}
    for fields in parser.values():
        fields[1]["write"] = lambda value: format(value, ".17g")
    text = tree.write(parser=parser, props=[], format_root_node=True) + "\n"

    def signature(candidate):
        descendants = candidate.get_cached_content()
        return {frozenset(tip.name for tip in leaves):
                (node.name or "", node.props.get("dist"))
                for node, leaves in descendants.items()}

    # Verify the root, every descendant clade, labels and branch lengths before
    # publishing the plain copy. Supports and other annotations are tool-irrelevant.
    from ete4 import Tree

    if signature(Tree(text, parser=1)) != signature(tree):
        raise ValueError("GeneRax serialization changed the rooted species tree.")
    with output_transaction([output_path], create_parents=True) as staged:
        Path(staged[output_path]).write_text(text, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    prepare_species_tree(args.input, args.output)


if __name__ == "__main__":
    main()
