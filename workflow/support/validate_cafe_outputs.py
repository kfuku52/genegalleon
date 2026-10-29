#!/usr/bin/env python3
"""Check saved CAFE reconstructions against the current taxa and extant counts.

CAFE can omit input families (for example root-absent families). Their absence
is reported rather than converted to zero or rejected as an identity mismatch.
This does not recover the settings of an unrecorded historical fit.
"""
import argparse
import csv
import io
import json
import re
from decimal import Decimal, InvalidOperation
from pathlib import Path

from Bio import Phylo
from cafe_branch_map import branch_inventory, first_newick


def read_table(path, family_column="FamilyID"):
    with Path(path).open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if not header or len(header) != len(set(header)) or any(not h.strip() for h in header):
            raise ValueError(f"Missing or duplicate CAFE table columns: {path}")
        family_headers = [i for i, name in enumerate(header) if name in {"FamilyID", "#FamilyID"}]
        if len(family_headers) != 1:
            raise ValueError(f"Missing or ambiguous {family_column} column: {path}")
        index = family_headers[0]
        header[index] = "FamilyID"
        rows = {}
        for row in reader:
            if not row:
                continue
            if len(row) != len(header) or not row[index].strip() or row[index] in rows:
                raise ValueError(f"Malformed or duplicate CAFE family row: {path}")
            rows[row[index]] = dict(zip(header, row, strict=True))
    return header, rows


def count(value, signed=False):
    try:
        number = Decimal(value)
    except InvalidOperation as exc:
        raise ValueError(f"Invalid CAFE copy number: {value!r}") from exc
    if not number.is_finite() or (not signed and number < 0) or number != number.to_integral_value():
        raise ValueError(f"Invalid CAFE copy number: {value!r}")
    return int(number)


def clades(tree):
    return {frozenset(t.name for t in node.get_terminals()) for node in tree.find_clades()}


def validate(copy_number, tree_path, output_dir):
    columns, inputs = read_table(copy_number)
    input_species = set(columns) - {"Desc", "FamilyID"}
    tree = Phylo.read(tree_path, "newick")
    tips = [t.name for t in tree.get_terminals()]
    if any(not tip for tip in tips) or len(tips) != len(set(tips)) or set(tips) != input_species:
        raise ValueError("CAFE copy-number species do not match the dated species tree")
    for row in inputs.values():
        for species in input_species:
            count(row[species])
    root = Path(output_dir)
    species, inventory = branch_inventory(root / "Gamma_asr.tre")
    if species != input_species:
        raise ValueError("Saved CAFE taxa do not match the current input species")
    reconstructed_tree = Phylo.read(io.StringIO(first_newick(root / "Gamma_asr.tre")), "newick")
    for tip in reconstructed_tree.get_terminals():
        tip.name = tip.name.replace("*", "").split("<", 1)[0]
    if clades(reconstructed_tree) != clades(tree):
        raise ValueError("Saved CAFE topology does not match the dated species tree")
    branch_keys = {r["native_key"] for r in inventory}
    count_columns, counts = read_table(root / "Gamma_count.tab")
    change_columns, changes = read_table(root / "Gamma_change.tab")
    if set(count_columns) - {"FamilyID"} != branch_keys or set(change_columns) - {"FamilyID"} != branch_keys:
        raise ValueError("Saved CAFE count/change branches do not match the reconstructed tree")
    if not counts or set(counts) != set(changes):
        raise ValueError("Saved CAFE count/change family identities disagree or are empty")
    leaf_keys = {r["descendants"][0]: r["native_key"] for r in inventory if len(r["descendants"]) == 1}
    if not set(counts) <= set(inputs):
        raise ValueError("Saved CAFE reconstruction contains families outside the current copy-number input")
    tree_families = set()
    with (root / "Gamma_asr.tre").open() as handle:
        for line in handle:
            match = re.match(r"^\s*TREE\s+(\S+)\s*=", line, re.IGNORECASE)
            if match:
                family = match[1]
                if family in tree_families or family not in counts:
                    raise ValueError("Saved CAFE ancestral tree contains duplicate or foreign family identities")
                tree_families.add(family)
    if not tree_families:
        raise ValueError("Saved CAFE ancestral tree has no named family reconstructions")
    for family, row in counts.items():
        for key in branch_keys:
            count(row[key])
            count(changes[family][key], signed=True)
        for taxon, key in leaf_keys.items():
            if count(row[key]) != count(inputs[family][taxon]):
                raise ValueError(f"Saved CAFE extant count does not match input: {family}/{taxon}")
    for name, field in [("Gamma_branch_probabilities.tab", "FamilyID"),
                        ("Gamma_family_results.txt", "#FamilyID")]:
        _, rows = read_table(root / name, field)
        if not set(rows) <= set(inputs):
            raise ValueError(f"Saved CAFE output contains foreign family identities: {name}")
    return {"validated": True, "species": len(species), "reconstructed_families": len(counts),
            "input_families_without_reconstruction": len(set(inputs) - set(counts)),
            "scope": "current taxa, topology, unique families and observed extant counts; historical fit settings require provenance"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--copy-number", required=True, type=Path)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        result = validate(args.copy_number, args.tree, args.output_dir)
    except (ValueError, OSError) as exc:
        parser.exit(1, f"CAFE semantic validation failed: {exc}\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
