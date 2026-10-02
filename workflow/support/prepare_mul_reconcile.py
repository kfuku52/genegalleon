#!/usr/bin/env python3
"""Prepare legacy gene_species labels without editing Newick syntax as text."""

import argparse
from io import StringIO
from pathlib import Path
from types import SimpleNamespace

from nwkit.mul_reconcile_model import validate_binary
from nwkit.output_transaction import output_transaction, validate_output_targets
from nwkit.species_parser import get_species_parser
from nwkit.util import read_tree, validate_outputs_do_not_replace_inputs, validate_unique_named_leaves, write_tree


def input_text(tree):
    for node in tree.traverse():
        if not node.is_leaf:
            node.name = None
    output = StringIO()
    write_tree(tree, SimpleNamespace(outfile=output), 1, quiet=True, props=[])
    return output.getvalue()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-tree", required=True)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--gene-tree", action="append")
    inputs.add_argument("--gene-tree-dir")
    parser.add_argument("--species-parser", default="legacy")
    parser.add_argument("--species-regex", default=None)
    parser.add_argument("--species-map-tsv", default=None)
    parser.add_argument("--species-out", required=True)
    parser.add_argument("--genes-out", required=True)
    parser.add_argument("--names-out", required=True)
    args = parser.parse_args()
    if args.gene_tree_dir is not None:
        args.gene_tree = [
            str(path) for path in sorted(Path(args.gene_tree_dir).glob("*.nwk"))
            if not path.name.startswith(".")
        ]
    if not args.gene_tree:
        raise ValueError("At least one gene tree is required.")
    validate_output_targets([args.species_out, args.genes_out, args.names_out])
    validate_outputs_do_not_replace_inputs(
        [
            ("--species-tree", args.species_tree),
            ("--species-map-tsv", args.species_map_tsv),
            *(("--gene-tree", path) for path in args.gene_tree),
        ],
        [("--species-out", args.species_out), ("--genes-out", args.genes_out), ("--names-out", args.names_out)],
    )
    species = read_tree(args.species_tree, "auto", True, quiet=True)
    validate_binary(species, "Species tree")
    names = set(species.leaf_names())
    for leaf in species.leaves():
        leaf.name = leaf.name.replace("_", "-")
    validate_unique_named_leaves(species, "Normalized species tree")
    species_parser = get_species_parser(args)
    genes, filenames = [], []
    for path in args.gene_tree:
        gene = read_tree(path, "auto", True, quiet=True)
        validate_binary(gene, f"Gene tree {path}")
        for leaf in gene.leaves():
            name = species_parser.parse(leaf.name).species_label
            if name not in names:
                raise ValueError(f"Gene tip does not match a species: {leaf.name} ({path})")
            prefix = name + "_"
            gene_id = leaf.name[len(prefix) :] if leaf.name.startswith(prefix) else leaf.name
            if not gene_id:
                raise ValueError(f"Gene tip has no gene identifier: {leaf.name}")
            leaf.name = gene_id.replace("_", "-") + "_" + name.replace("_", "-")
        validate_unique_named_leaves(gene, f"Normalized gene tree {path}")
        genes.append(input_text(gene))
        filename = Path(path).name
        if any(character in filename for character in "\t\r\n"):
            raise ValueError(f"Gene-tree filename contains a table delimiter: {path!r}")
        filenames.append(filename)
    contents = {
        args.species_out: input_text(species) + "\n",
        args.genes_out: "\n".join(genes) + "\n",
        args.names_out: "\n".join(filenames) + "\n",
    }
    with output_transaction(contents) as staged:
        for path, content in contents.items():
            staged.write_text(path, lambda handle, content=content: handle.write(content))


if __name__ == "__main__":
    main()
