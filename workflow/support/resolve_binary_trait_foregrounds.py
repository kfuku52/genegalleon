#!/usr/bin/env python3

import argparse
import csv
import math
from collections import Counter
from pathlib import Path

from Bio import Phylo
from species_labeling import extract_species_label

MISSING_VALUES = frozenset(("", ".", "na", "nan", "none", "null"))
NOT_BINARY = object()


def normalize_species_name(value):
    value = str(value or "").strip()
    return extract_species_label(value) or value


def match_species_name(value, tree_species):
    raw_name = str(value or "").strip()
    if raw_name in tree_species:
        return raw_name
    return normalize_species_name(raw_name)


def binary_value(value):
    text = "" if value is None else str(value).strip()
    if text.lower() in MISSING_VALUES:
        return None
    try:
        number = float(text)
    except ValueError:
        return NOT_BINARY
    if not math.isfinite(number) or number not in (0.0, 1.0):
        return NOT_BINARY
    return int(number)


def binary_trait_values(rows, column_index):
    parsed = [binary_value(row[column_index]) for row in rows]
    observed = [value for value in parsed if value is not None]
    if not observed or any(value is NOT_BINARY for value in observed):
        return None
    return parsed


def validate_tree(tree):
    terminal_names = [str(terminal.name or "").strip() for terminal in tree.get_terminals()]
    empty_names = sum(name == "" for name in terminal_names)
    if empty_names:
        raise ValueError(f"Species tree contains {empty_names} unnamed terminal(s).")
    duplicate_names = sorted(name for name, count in Counter(terminal_names).items() if count > 1)
    if duplicate_names:
        raise ValueError("Species tree contains duplicate terminal names: " + ", ".join(duplicate_names))
    return set(terminal_names)


def foreground_components(tree, foreground_species):
    parent_by_clade = {
        child: parent
        for parent in tree.find_clades(order="level")
        for child in parent.clades
    }
    all_foreground = {}
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            all_foreground[clade] = clade.name in foreground_species
        else:
            all_foreground[clade] = bool(clade.clades) and all(all_foreground[child] for child in clade.clades)

    component_leaves = []
    for clade in tree.find_clades(order="preorder"):
        parent = parent_by_clade.get(clade)
        if all_foreground[clade] and (parent is None or not all_foreground[parent]):
            component_leaves.append(tuple(sorted(terminal.name for terminal in clade.get_terminals())))
    return sorted(component_leaves)


def resolve_binary_foregrounds(header, rows, tree):
    if len(header) < 2:
        raise ValueError("Trait table must contain a species column and at least one trait column.")
    if any(len(row) != len(header) for row in rows):
        raise ValueError("Trait table contains rows whose field count differs from the header.")

    tree_species = validate_tree(tree)
    normalized_species = [match_species_name(row[0], tree_species) for row in rows]
    duplicate_species = sorted(name for name, count in Counter(normalized_species).items() if count > 1)
    if duplicate_species:
        raise ValueError("Trait table contains duplicate species after normalization: " + ", ".join(duplicate_species))

    resolved_traits = []
    for column_index, trait_name in enumerate(header[1:], start=1):
        parsed = binary_trait_values(rows, column_index)
        if parsed is None:
            continue

        foreground_species = {
            species
            for species, value in zip(normalized_species, parsed, strict=True)
            if value == 1
        }
        missing_foreground = sorted(foreground_species - tree_species)
        if missing_foreground:
            raise ValueError(
                f"Trait {trait_name!r} has foreground species absent from the species tree: "
                + ", ".join(missing_foreground)
            )

        components = foreground_components(tree, foreground_species)
        lineage_by_species = {
            species: lineage_id
            for lineage_id, component in enumerate(components, start=1)
            for species in component
        }
        for row, species, value in zip(rows, normalized_species, parsed, strict=True):
            if value == 1:
                row[column_index] = str(lineage_by_species[species])
        resolved_traits.append((trait_name, len(components)))
    return resolved_traits


def read_trait_table(path):
    with Path(path).open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError("Trait table is empty.") from exc
        return header, [row for row in reader if row]


def write_trait_table(path, header, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Assign stable lineage IDs to maximal foreground clades in binary species-trait columns."
    )
    parser.add_argument("--input", required=True, help="Input species_trait TSV.")
    parser.add_argument("--species-tree", required=True, help="Species tree in Newick format.")
    parser.add_argument("--output", required=True, help="Output TSV with resolved binary foreground IDs.")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    header, rows = read_trait_table(args.input)
    tree = Phylo.read(args.species_tree, "newick")
    resolved_traits = resolve_binary_foregrounds(header, rows, tree)
    write_trait_table(args.output, header, rows)
    for trait_name, lineage_count in resolved_traits:
        print(f"Resolved binary foreground trait {trait_name}: {lineage_count} lineage(s).", flush=True)


if __name__ == "__main__":
    main()
