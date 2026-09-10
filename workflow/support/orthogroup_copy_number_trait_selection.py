#!/usr/bin/env python3
"""Prepare orthogroup predictors for NWKIT exploratory nested-CV selection.

The manifest identifies the current result bundle. All traits must finish before
any result is published; existing unrelated or historical files are preserved.
"""

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from ete4.parser.newick import PARSERS
from nwkit.file_paths import validate_outputs_do_not_replace_inputs
from nwkit.output_transaction import output_transaction
from nwkit.rooting_state import require_rooted
from nwkit.util import read_tree

SUPPORT = Path(__file__).resolve().parent
if str(SUPPORT) not in sys.path:
    sys.path.insert(0, str(SUPPORT))

from species_trait_contract import metadata_path, select_analysis_traits, write_derived_contract
from species_trait_schema import schema_path

FAMILIES = {"gaussian", "binomial", "poisson", "negative-binomial"}


def tokens(text):
    return [token for token in re.split(r"[,\s]+", text.strip()) if token]


def family_mapping(text, traits):
    mapping = dict.fromkeys(traits, "gaussian")
    seen = set()
    for token in tokens(text):
        pair = token.split("=")
        if len(pair) != 2 or pair[0] not in mapping or pair[1] not in FAMILIES or pair[0] in seen:
            raise ValueError(f"Invalid or duplicate response family mapping: {token}")
        mapping[pair[0]] = pair[1]
        seen.add(pair[0])
    return mapping


def read_table(path):
    with open(path, encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if not header or any(not name.strip() for name in header) or len(header) != len(set(header)):
            raise ValueError(f"Table headers must be unique and non-empty: {path}")
        for line, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise ValueError(f"Table row {line} has the wrong number of fields: {path}")
    return pd.read_csv(path, sep="\t", keep_default_na=False, dtype=str)


def load_species_table(path, table=None):
    if table is None:
        table = read_table(path)
    if "leaf_name" in table.columns[1:]:
        raise ValueError(f"Species key conflicts with a leaf_name data column: {path}")
    table = table.rename(columns={table.columns[0]: "leaf_name"})
    table["leaf_name"] = table["leaf_name"].str.replace(" ", "_", regex=False)
    if table.leaf_name.duplicated().any() or (table.leaf_name == "").any():
        raise ValueError(f"Species keys must be unique and non-empty: {path}")
    return table.set_index("leaf_name")


def write_species_subset(tree, matched, path):
    """Retain the original root and branch precision when traits are missing."""
    use_tree = tree.copy()
    retained = set()
    for leaf in use_tree.leaves():
        if leaf.name in matched:
            node = leaf
            while node is not None:
                retained.add(node)
                node = node.up
    for node in list(use_tree.traverse("postorder")):
        if node not in retained and node.up is not None:
            node.detach()
    # Keep unary nodes, including the root: collapsing the root loses shared
    # Brownian covariance when all remaining species belong to one clade.
    parser = {key: [dict(item) for item in fields] for key, fields in PARSERS[1].items()}
    for fields in parser.values():
        fields[1]["write"] = lambda value: format(float(value), ".17g")
    use_tree.write(outfile=str(path), parser=parser)


def select_families(available, ids, file, maximum):
    requested = tokens(ids)
    if file:
        lines = Path(file).read_text().splitlines()
        file_ids = [tokens(line)[0] for line in lines if line.strip() and not line.lstrip().startswith("#")]
        if file_ids and file_ids[0].lower() in {"family_id", "query", "orthogroup", "gene_family_id", "id"}:
            file_ids = file_ids[1:]
        requested += file_ids
    if requested:
        requested = list(dict.fromkeys(requested))
        if not set(requested) <= set(available):
            raise ValueError("Requested orthogroup IDs are absent from the copy-number table.")
        return requested
    if maximum.lower() in {"all", "auto", "0"}:
        return list(available)
    value = int(maximum)
    if value < 0:
        raise ValueError("max-families must be non-negative or all.")
    return list(available)[:value]


def run(args):
    tree = read_tree(args.tree, "auto", True)
    require_rooted(tree, "Copy-number selection requires a rooted dated species tree.")
    for leaf in tree.leaves():
        leaf.name = leaf.name.replace(" ", "_")
    leaves = [leaf.name for leaf in tree.leaves()]
    if len(set(leaves)) != len(leaves):
        raise ValueError("Duplicate normalized species-tree labels.")
    counts = read_table(args.copy_number)
    if "Orthogroup" not in counts or counts.Orthogroup.duplicated().any() or (counts.Orthogroup == "").any():
        raise ValueError("Copy-number table requires unique non-empty Orthogroup IDs.")
    counts = counts.set_index("Orthogroup")
    counts.columns = counts.columns.str.replace(" ", "_", regex=False)
    if counts.columns.duplicated().any():
        raise ValueError("Copy-number columns collide after species-label normalization.")
    chosen = select_families(counts.index, args.family_ids, args.family_file, args.max_families)
    if not chosen:
        raise ValueError("No orthogroups selected.")
    copy_matrix = counts.loc[chosen, leaves].T.apply(pd.to_numeric, errors="raise")
    values = copy_matrix.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any():
        raise ValueError("Predictor copy numbers must be finite non-negative integers.")
    # A fixed transform uses no held-out statistics; NWKIT standardizes within folds.
    copy_matrix = np.log1p(copy_matrix)
    selected_traits, trait_audit = select_analysis_traits(args.traits, args.trait)
    trait_table = load_species_table(args.traits, selected_traits)
    folds = load_species_table(args.folds)
    traits = list(trait_table.columns)
    trait_selection = trait_audit["type_selection"]
    if not traits or len(set(traits)) != len(traits) or not set(traits) <= set(trait_table.columns):
        raise ValueError("Trait selection must contain unique existing trait columns.")
    if "fold" not in folds or not set(leaves) <= set(folds.index):
        raise ValueError("Fold file must supply a fold for every tree species.")
    mapping = family_mapping(args.response_families, traits)
    # Stable internal column names prevent gene/trait/key namespace collisions.
    predictor_names = [f"predictor_{i:06d}" for i in range(len(chosen))]
    copy_matrix.columns = predictor_names
    executable = shutil.which("nwkit")
    if not executable:
        raise ValueError("nwkit regress-select is required.")
    with tempfile.TemporaryDirectory(prefix="gg-copy-selection-") as temp:
        work = Path(temp)
        result = work / "results"
        result.mkdir()
        selected_path = result / "selected_species_traits.tsv"
        selected_traits.to_csv(selected_path, sep="\t", index=False)
        write_derived_contract(selected_path, trait_audit)
        (result / "species_trait_input.json").write_text(json.dumps(trait_audit, indent=2) + "\n")
        pd.DataFrame({"predictor": predictor_names, "Orthogroup": chosen, "predictor_transform": "log1p"}).to_csv(
            result / "predictors.tsv", sep="\t", index=False
        )
        (work / "predictors.txt").write_text("\n".join(predictor_names) + "\n")
        pd.DataFrame(trait_selection).to_csv(result / "trait_selection.tsv", sep="\t", index=False)
        manifest = []
        for i, trait in enumerate(traits):
            raw_response = trait_table[trait].reindex(leaves)
            missing = raw_response.isna() | raw_response.isin(["", "NA", "NaN", "nan"])
            numeric = pd.to_numeric(raw_response[~missing], errors="raise")
            if not np.isfinite(numeric.to_numpy()).all():
                raise ValueError(f"Trait {trait} contains non-finite values.")
            matched = [leaf for leaf in leaves if leaf in numeric.index]
            frame = copy_matrix.loc[matched].copy()
            frame["response"] = numeric.loc[matched]
            frame.index.name = "leaf_name"
            frame.to_csv(work / "data.tsv", sep="\t")
            folds.loc[matched, ["fold"]].to_csv(work / "folds.tsv", sep="\t")
            write_species_subset(tree, set(matched), work / "tree.nwk")
            name = f"trait_{i + 1:04d}"
            subprocess.run(
                [
                    executable,
                    "regress-select",
                    "--input-rooted",
                    "yes",
                    "--tree",
                    str(work / "tree.nwk"),
                    "--data",
                    str(work / "data.tsv"),
                    "--folds",
                    str(work / "folds.tsv"),
                    "--response",
                    "response",
                    "--predictor-file",
                    str(work / "predictors.txt"),
                    "--family",
                    mapping[trait],
                    "--strengths",
                    args.strengths,
                    "--l1-ratios",
                    args.l1_ratios,
                    "--prediction",
                    args.prediction,
                    "--out-prefix",
                    str(result / name),
                ],
                check=True,
            )
            for suffix in ("coefficients.tsv", "path.tsv", "stability.tsv"):
                result_path = result / f"{name}.{suffix}"
                table = pd.read_csv(result_path, sep="\t", keep_default_na=False)
                table.insert(
                    1, "Orthogroup", table["term"].map(dict(zip(predictor_names, chosen, strict=True))).fillna("")
                )
                table["predictor_transform"] = "log1p"
                table.to_csv(result_path, sep="\t", index=False)
            result_metadata_path = result / f"{name}.metadata.json"
            metadata = json.loads(result_metadata_path.read_text())
            metadata["trait_input_audit"] = trait_audit
            metadata["response_trait"] = trait
            metadata["predictor_transform"] = "log1p"
            metadata["predictor_units"] = "natural_log_of_one_plus_copy_number"
            result_metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
            manifest.append(
                dict(
                    trait=trait,
                    response_family=mapping[trait],
                    predictor_transform="log1p",
                    n_species=len(matched),
                    n_missing_response=len(leaves) - len(matched),
                    prefix=name,
                    inference_status="exploratory_no_post_selection_inference",
                )
            )
        pd.DataFrame(manifest).to_csv(result / "manifest.tsv", sep="\t", index=False)
        sources = sorted(result.iterdir())
        targets = [Path(args.outdir) / source.name for source in sources]
        inputs = [
            ("copy_number", args.copy_number),
            ("tree", args.tree),
            ("traits", args.traits),
            ("folds", args.folds),
        ]
        if schema_path(args.traits).exists():
            inputs.append(("trait_schema", str(schema_path(args.traits))))
        if metadata_path(Path(args.traits)).exists():
            inputs.append(("trait_metadata", str(metadata_path(Path(args.traits)))))
        if args.family_file:
            inputs.append(("family_file", args.family_file))
        validate_outputs_do_not_replace_inputs(inputs, [(str(i), str(path)) for i, path in enumerate(targets)])
        with output_transaction(targets, create_parents=True) as staged:
            for source, target in zip(sources, targets, strict=True):
                shutil.copyfile(source, staged[target])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("copy-number", "tree", "traits", "folds", "outdir"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--trait", default="all")
    parser.add_argument("--response-families", default="")
    parser.add_argument("--family-ids", default="")
    parser.add_argument("--family-file", default="")
    parser.add_argument("--max-families", default="all")
    parser.add_argument("--strengths", default="1,0.1,0.01")
    parser.add_argument("--l1-ratios", default="1,0.5")
    parser.add_argument("--prediction", choices=("conditional", "fixed"), default="conditional")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
