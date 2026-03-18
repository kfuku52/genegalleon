#!/usr/bin/env python3

import argparse
import os
import re

import numpy
import pandas


BRANCH_PREFIX_COLUMNS = [
    "candidate_gene_count",
    "matched_leaf_count",
    "besthit_gene_count",
    "besthit_taxid_count",
    "besthit_taxonomy_method",
    "besthit_same_superkingdom_fraction",
    "besthit_lca_rank_mode",
    "intron_support_fraction",
    "expression_measured_fraction",
    "clade_min_expression_pearsoncor",
    "synteny_support_fraction",
    "synteny_mean_support_score",
    "contamination_incompatible_fraction",
    "contamination_top_lca_sciname",
]


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Annotate a stat_branch table with HGT-specific branch and tip evidence columns."
    )
    parser.add_argument("--stat_branch", metavar="PATH", required=True, type=str)
    parser.add_argument("--branch_tsv", metavar="PATH", required=True, type=str)
    parser.add_argument("--gene_tsv", metavar="PATH", required=True, type=str)
    parser.add_argument("--orthogroup", metavar="TEXT", required=True, type=str)
    parser.add_argument("--outfile", metavar="PATH", required=True, type=str)
    return parser


def ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def safe_read_tsv(path: str) -> pandas.DataFrame:
    if not path or not os.path.exists(path):
        return pandas.DataFrame()
    try:
        return pandas.read_csv(path, sep="\t", low_memory=False)
    except pandas.errors.EmptyDataError:
        return pandas.DataFrame()


def parse_boolish(value):
    if pandas.isna(value):
        return pandas.NA
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y"}:
        return True
    if text in {"0", "false", "f", "no", "n"}:
        return False
    return pandas.NA


def bool_to_numeric(value):
    parsed = parse_boolish(value)
    if parsed is pandas.NA:
        return numpy.nan
    return 1.0 if parsed else 0.0


def normalize_text(value) -> str:
    if pandas.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() == "nan":
        return ""
    text = text.replace("_", " ")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def abbreviate_species_name(value) -> str:
    text = normalize_text(value)
    if text == "":
        return ""
    parts = text.split()
    if len(parts) >= 2 and re.match(r"^[A-Za-z-]+$", parts[0]):
        text = f"{parts[0][0]}. {parts[1]}"
        if len(parts) >= 3:
            text = f"{text} {parts[2]}"
    return text


def shorten_text(value, max_nchar: int = 18) -> str:
    text = normalize_text(value)
    if text == "":
        return ""
    if len(text) <= max_nchar:
        return text
    if max_nchar <= 3:
        return text[:max_nchar]
    return text[: max_nchar - 3] + "..."


def display_text(value, fallback: str = "-", abbreviate: bool = False) -> str:
    text = abbreviate_species_name(value) if abbreviate else normalize_text(value)
    if text == "":
        return fallback
    return text


def annotate_branch_rows(stat_branch: pandas.DataFrame, branch_df: pandas.DataFrame, orthogroup: str) -> pandas.DataFrame:
    branch_subset = branch_df.loc[branch_df.get("orthogroup", pandas.Series(dtype=object)).astype(str) == orthogroup, :].copy()
    if branch_subset.empty:
        stat_branch["hgtbranch_candidate_flag"] = 0
        return stat_branch
    merge_cols = ["branch_id"] + [col for col in BRANCH_PREFIX_COLUMNS if col in branch_subset.columns]
    branch_subset = branch_subset.loc[:, merge_cols].copy()
    rename_map = {col: f"hgtbranch_{col}" for col in merge_cols if col != "branch_id"}
    branch_subset = branch_subset.rename(columns=rename_map)
    out = stat_branch.merge(branch_subset, on="branch_id", how="left", sort=False)
    out["hgtbranch_candidate_flag"] = out["hgtbranch_candidate_gene_count"].notna().astype(int)
    return out


def annotate_leaf_rows(stat_branch: pandas.DataFrame, gene_df: pandas.DataFrame, orthogroup: str) -> pandas.DataFrame:
    gene_subset = gene_df.loc[gene_df.get("orthogroup", pandas.Series(dtype=object)).astype(str) == orthogroup, :].copy()
    if gene_subset.empty:
        if "so_event" in stat_branch.columns:
            leaf_mask = stat_branch["so_event"].astype(str).eq("L")
            stat_branch.loc[leaf_mask, "hgt_Cand"] = 0
        return stat_branch

    gene_subset["node_name"] = gene_subset["gene_id"].astype(str)
    gene_subset["hgt_Cand"] = pandas.to_numeric(gene_subset.get("candidate_branch_count", numpy.nan), errors="coerce")
    gene_subset["hgt_SameSK"] = pandas.to_numeric(gene_subset.get("besthit_same_superkingdom", numpy.nan), errors="coerce")
    gene_subset["hgt_ContamOK"] = gene_subset.get("contamination_is_compatible_lineage", pandas.Series(dtype=object)).map(bool_to_numeric)
    gene_subset["hgt_Expr"] = gene_subset.get("expression_measured", pandas.Series(dtype=object)).map(bool_to_numeric)
    gene_subset["hgt_Intron"] = gene_subset.get("intron_supported", pandas.Series(dtype=object)).map(bool_to_numeric)
    gene_subset["hgt_Syn"] = pandas.to_numeric(gene_subset.get("synteny_support_score", numpy.nan), errors="coerce")
    gene_subset["besthit_lca_rank_short"] = gene_subset.get("besthit_lca_rank", pandas.Series(dtype=object)).map(
        lambda x: shorten_text(x, 18)
    )
    gene_subset["besthit_lca_rank_display"] = gene_subset.get("besthit_lca_rank", pandas.Series(dtype=object)).map(
        lambda x: display_text(x, fallback="-", abbreviate=False)
    )
    gene_subset["besthit_organism_short"] = gene_subset.get("besthit_organism", pandas.Series(dtype=object)).map(
        lambda x: shorten_text(abbreviate_species_name(x), 28)
    )
    gene_subset["besthit_organism_display"] = gene_subset.get("besthit_organism", pandas.Series(dtype=object)).map(
        lambda x: display_text(x, fallback="-", abbreviate=True)
    )
    gene_subset["contamination_lca_sciname_short"] = gene_subset.get("contamination_lca_sciname", pandas.Series(dtype=object)).map(
        lambda x: shorten_text(abbreviate_species_name(x), 28)
    )
    gene_subset["contamination_lca_sciname_display"] = gene_subset.get("contamination_lca_sciname", pandas.Series(dtype=object)).map(
        lambda x: display_text(x, fallback="-", abbreviate=True)
    )

    keep_cols = [
        "node_name",
        "candidate_branch_ids",
        "hgt_Cand",
        "hgt_SameSK",
        "hgt_ContamOK",
        "hgt_Expr",
        "hgt_Intron",
        "hgt_Syn",
        "besthit_lca_rank_short",
        "besthit_lca_rank_display",
        "besthit_organism_short",
        "besthit_organism_display",
        "contamination_lca_sciname_short",
        "contamination_lca_sciname_display",
    ]
    keep_cols = [col for col in keep_cols if col in gene_subset.columns]
    gene_subset = gene_subset.loc[:, keep_cols].drop_duplicates(subset=["node_name"], keep="first")
    out = stat_branch.merge(gene_subset, on="node_name", how="left", sort=False)

    if "so_event" in out.columns:
        leaf_mask = out["so_event"].astype(str).eq("L")
        out.loc[leaf_mask & out["hgt_Cand"].isna(), "hgt_Cand"] = 0
    return out


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    stat_branch = safe_read_tsv(args.stat_branch)
    if stat_branch.empty:
        ensure_parent_dir(args.outfile)
        pandas.DataFrame().to_csv(args.outfile, sep="\t", index=False)
        return

    branch_df = safe_read_tsv(args.branch_tsv)
    gene_df = safe_read_tsv(args.gene_tsv)

    out = stat_branch.copy()
    out = annotate_branch_rows(out, branch_df, args.orthogroup)
    out = annotate_leaf_rows(out, gene_df, args.orthogroup)

    ensure_parent_dir(args.outfile)
    out.to_csv(args.outfile, sep="\t", index=False)


if __name__ == "__main__":
    main()
