#!/usr/bin/env python3
"""Prepare legacy gene_species labels without editing Newick syntax as text."""

import argparse
import json
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


def rooted_clades(tree):
    return {frozenset(node.leaf_names()) for node in tree.traverse() if not node.is_leaf}


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
    parser.add_argument("--locus-model")
    parser.add_argument("--locus-model-out")
    parser.add_argument("--locus-species-tree")
    parser.add_argument("--locus-species-out")
    args = parser.parse_args()
    locus_options = (args.locus_model, args.locus_model_out, args.locus_species_tree, args.locus_species_out)
    if any(locus_options) and not all(locus_options):
        parser.error("All four --locus-model/--locus-species input/output options are required together.")
    if args.gene_tree_dir is not None:
        args.gene_tree = [
            str(path) for path in sorted(Path(args.gene_tree_dir).glob("*.nwk"))
            if not path.name.startswith(".")
        ]
    if not args.gene_tree:
        raise ValueError("At least one gene tree is required.")
    outputs = [("--species-out", args.species_out), ("--genes-out", args.genes_out), ("--names-out", args.names_out)]
    if args.locus_model_out:
        outputs.append(("--locus-model-out", args.locus_model_out))
        outputs.append(("--locus-species-out", args.locus_species_out))
    validate_output_targets([path for _, path in outputs])
    validate_outputs_do_not_replace_inputs(
        [
            ("--species-tree", args.species_tree),
            ("--species-map-tsv", args.species_map_tsv),
            ("--locus-model", args.locus_model),
            ("--locus-species-tree", args.locus_species_tree),
            *(("--gene-tree", path) for path in args.gene_tree),
        ],
        outputs,
    )
    species = read_tree(args.species_tree, "auto", True, quiet=True)
    validate_binary(species, "Species tree")
    names = set(species.leaf_names())
    clades = rooted_clades(species)
    child_order = {
        frozenset(node.leaf_names()): {frozenset(child.leaf_names()): i for i, child in enumerate(node.children)}
        for node in species.traverse() if not node.is_leaf
    }
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
    if args.locus_model:
        from nwkit.mul_locus_cli import json_pairs
        from nwkit.mul_locus_mc import validate_model

        with open(args.locus_model) as handle:
            model = json.load(handle, object_pairs_hook=json_pairs)
        locus_species = read_tree(args.locus_species_tree, "auto", True, quiet=True)
        validate_binary(locus_species, "Locus species tree")
        if set(locus_species.leaf_names()) != names or rooted_clades(locus_species) != clades:
            raise ValueError("Locus species tree must have the same species and rooted topology.")
        # Numeric H1/H2 selectors use postorder IDs, so match the original order.
        for node in locus_species.traverse():
            if not node.is_leaf:
                order = child_order[frozenset(node.leaf_names())]
                node.children.sort(key=lambda child: order[frozenset(child.leaf_names())])
        for leaf in locus_species.leaves():
            leaf.name = leaf.name.replace("_", "-")
        if not isinstance(model, dict) or not isinstance(model.get("detection"), dict) or set(model["detection"]) != names:
            raise ValueError("Locus detection must match the original species labels exactly.")
        model["detection"] = {name.replace("_", "-"): p for name, p in model["detection"].items()}
        validate_model(model, locus_species)
        contents[args.locus_model_out] = json.dumps(model, indent=2, allow_nan=False) + "\n"
        contents[args.locus_species_out] = input_text(locus_species) + "\n"
    with output_transaction(contents) as staged:
        for path, content in contents.items():
            staged.write_text(path, lambda handle, content=content: handle.write(content))


if __name__ == "__main__":
    main()
