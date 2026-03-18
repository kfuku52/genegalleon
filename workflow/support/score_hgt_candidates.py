#!/usr/bin/env python3

import argparse
from collections import Counter
import gzip
import os
import re
import sqlite3
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy
import pandas

try:
    from ete4 import NCBITaxa
except Exception:  # pragma: no cover - optional runtime dependency
    NCBITaxa = None


BRANCH_OUTPUT_COLUMNS = [
    "orthogroup",
    "branch_id",
    "node_name",
    "generax_event",
    "generax_transfer",
    "generax_event_parent",
    "taxon",
    "spnode_coverage",
    "candidate_gene_count",
    "matched_leaf_count",
    "candidate_genes",
    "besthit_gene_count",
    "besthit_taxid_count",
    "besthit_taxonomy_method",
    "besthit_same_superkingdom_fraction",
    "besthit_lca_rank_mode",
    "intron_measured_gene_count",
    "intron_supported_gene_count",
    "intron_support_fraction",
    "expression_measured_gene_count",
    "expression_measured_fraction",
    "clade_min_expression_pearsoncor",
    "synteny_measured_gene_count",
    "synteny_supported_gene_count",
    "synteny_support_fraction",
    "synteny_mean_support_score",
    "contamination_measured_gene_count",
    "contamination_incompatible_gene_count",
    "contamination_incompatible_fraction",
    "contamination_top_lca_taxid",
    "contamination_top_lca_sciname",
    "contamination_top_lca_fraction",
]

GENE_OUTPUT_COLUMNS = [
    "orthogroup",
    "gene_id",
    "gene_taxon",
    "candidate_branch_count",
    "candidate_branch_ids",
    "besthit_accession",
    "besthit_organism",
    "besthit_taxid",
    "besthit_taxonomy_method",
    "besthit_lca_rank",
    "besthit_same_superkingdom",
    "intron_supported",
    "expression_measured",
    "synteny_support_score",
    "contamination_lca_taxid",
    "contamination_lca_sciname",
    "contamination_is_compatible_lineage",
]

ORTHOGROUP_OUTPUT_COLUMNS = [
    "orthogroup",
    "hgt_branch_count",
    "hgt_gene_count",
    "candidate_branch_ids",
]


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Summarize GeneRax-first HGT candidate evidence from the gg_orthogroup SQLite database."
    )
    parser.add_argument("--dbpath", metavar="PATH", required=True, type=str)
    parser.add_argument("--branch_out", metavar="PATH", required=True, type=str)
    parser.add_argument("--gene_out", metavar="PATH", required=True, type=str)
    parser.add_argument("--orthogroup_out", metavar="PATH", required=True, type=str)
    parser.add_argument("--dir_contamination_tsv", metavar="PATH", default="", type=str)
    parser.add_argument("--taxonomy_dbfile", metavar="PATH", default=None, type=str)
    return parser


def ensure_parent_dir(path: str) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def empty_frame(columns: Sequence[str]) -> pandas.DataFrame:
    return pandas.DataFrame(columns=list(columns))


def write_tsv(df: pandas.DataFrame, path: str, columns: Sequence[str]) -> None:
    ensure_parent_dir(path)
    out = df.copy()
    for col in columns:
        if col not in out.columns:
            out[col] = pandas.NA
    out = out.loc[:, list(columns)]
    out.to_csv(path, sep="\t", index=False)


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def sqlite_table_exists(conn: sqlite3.Connection, table_name: str) -> bool:
    query = "SELECT name FROM sqlite_master WHERE type='table' AND name=?"
    row = conn.execute(query, (table_name,)).fetchone()
    return row is not None


def sqlite_table_columns(conn: sqlite3.Connection, table_name: str) -> List[str]:
    rows = conn.execute(f'PRAGMA table_info("{table_name}")').fetchall()
    return [str(row[1]) for row in rows]


def quote_ident(name: str) -> str:
    return '"' + str(name).replace('"', '""') + '"'


def first_present(available: Iterable[str], candidates: Sequence[str]) -> str:
    available_set = set(available)
    for candidate in candidates:
        if candidate in available_set:
            return candidate
    return ""


def read_candidate_orthogroups(conn: sqlite3.Connection, branch_columns: Sequence[str]) -> List[str]:
    event_col = first_present(branch_columns, ["generax_event"])
    transfer_col = first_present(branch_columns, ["generax_transfer"])
    if event_col == "" and transfer_col == "":
        return []

    clauses = []
    if event_col:
        clauses.append(f'COALESCE({quote_ident(event_col)}, "") = "H"')
    if transfer_col:
        clauses.append(f'COALESCE({quote_ident(transfer_col)}, "") LIKE "Y%"')
    query = "SELECT DISTINCT orthogroup FROM branch WHERE " + " OR ".join(clauses)
    rows = conn.execute(query).fetchall()
    return [str(row[0]) for row in rows if row[0] not in (None, "")]


def read_branch_subset(
    conn: sqlite3.Connection,
    orthogroups: Sequence[str],
    columns: Sequence[str],
    chunk_size: int = 500,
) -> pandas.DataFrame:
    if len(orthogroups) == 0:
        return empty_frame(columns)
    frames = []
    select_cols = ", ".join(quote_ident(col) for col in columns)
    for start in range(0, len(orthogroups), chunk_size):
        chunk = orthogroups[start:start + chunk_size]
        placeholders = ", ".join("?" for _ in chunk)
        query = f"SELECT {select_cols} FROM branch WHERE orthogroup IN ({placeholders})"
        frames.append(pandas.read_sql_query(query, conn, params=list(chunk)))
    if len(frames) == 0:
        return empty_frame(columns)
    return pandas.concat(frames, ignore_index=True)


def split_gene_labels(value) -> List[str]:
    if pandas.isna(value):
        return []
    text = str(value).strip()
    if text == "":
        return []
    return [token.strip() for token in text.split(";") if token.strip() != ""]


def parse_boolish(value):
    if pandas.isna(value):
        return pandas.NA
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y"}:
        return True
    if text in {"0", "false", "f", "no", "n"}:
        return False
    return pandas.NA


def normalize_sci_name(name: str) -> str:
    text = str(name).strip().replace("_", " ")
    text = re.sub(r"\([^)]*\)", "", text)
    text = re.sub(r"\s+", " ", text).strip()
    if text == "":
        return ""
    parts = text.split()
    if len(parts) >= 2 and re.match(r"^[A-Z][a-zA-Z-]+$", parts[0]):
        return parts[0] + " " + parts[1]
    return text


def safe_float(value, default=numpy.nan):
    try:
        if pandas.isna(value):
            return default
    except TypeError:
        pass
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


class TaxonomyResolver:
    def __init__(self, dbfile: str = ""):
        self.dbfile = str(dbfile).strip()
        self.enabled = False
        self.ncbi = None
        self.name_cache: Dict[str, int] = {}
        self.lineage_cache: Dict[int, List[int]] = {}
        self.rank_cache: Dict[int, str] = {}
        if NCBITaxa is None or self.dbfile == "":
            return
        try:
            self.ncbi = NCBITaxa(dbfile=self.dbfile)
            self.enabled = True
        except Exception:
            self.ncbi = None
            self.enabled = False

    def resolve_name_taxid(self, name: str) -> int:
        normalized = normalize_sci_name(name)
        if normalized == "":
            return 0
        if normalized in self.name_cache:
            return self.name_cache[normalized]
        candidates = [normalized]
        if " " in normalized:
            genus = normalized.split()[0]
            candidates.append(genus)
        taxid = 0
        if self.enabled and self.ncbi is not None:
            for candidate in dict.fromkeys(candidates):
                try:
                    translated = self.ncbi.get_name_translator([candidate])
                except Exception:
                    translated = {}
                if candidate in translated and len(translated[candidate]) > 0:
                    taxid = int(translated[candidate][0])
                    break
        self.name_cache[normalized] = taxid
        return taxid

    def lineage(self, taxid: int) -> List[int]:
        taxid = int(taxid)
        if taxid <= 0:
            return []
        if taxid in self.lineage_cache:
            return self.lineage_cache[taxid]
        lineage = []
        if self.enabled and self.ncbi is not None:
            try:
                lineage = [int(item) for item in self.ncbi.get_lineage(taxid)]
            except Exception:
                lineage = []
        self.lineage_cache[taxid] = lineage
        return lineage

    def rank(self, taxid: int) -> str:
        taxid = int(taxid)
        if taxid <= 0:
            return ""
        if taxid in self.rank_cache:
            return self.rank_cache[taxid]
        rank_name = ""
        if self.enabled and self.ncbi is not None:
            try:
                rank_name = str(self.ncbi.get_rank([taxid]).get(taxid, ""))
            except Exception:
                rank_name = ""
        self.rank_cache[taxid] = rank_name
        return rank_name

    def rank_taxid_from_lineage(self, lineage: Sequence[int], target_ranks: Sequence[str]) -> int:
        target_set = {rank.lower() for rank in target_ranks}
        for taxid in lineage:
            rank_name = self.rank(taxid).lower()
            if rank_name in target_set:
                return int(taxid)
        return 0

    def compare(self, query_name: str, hit_taxid, hit_name: str = "") -> Optional[Dict[str, object]]:
        if not self.enabled:
            return None
        query_taxid = self.resolve_name_taxid(query_name)
        hit_taxid_val = 0
        if hit_taxid not in ("", None) and not pandas.isna(hit_taxid):
            hit_taxid_val = int(safe_float(hit_taxid, default=0))
        if hit_taxid_val <= 0 and hit_name:
            hit_taxid_val = self.resolve_name_taxid(hit_name)
        if query_taxid <= 0 or hit_taxid_val <= 0:
            return None
        query_lineage = self.lineage(query_taxid)
        hit_lineage = self.lineage(hit_taxid_val)
        if len(query_lineage) == 0 or len(hit_lineage) == 0:
            return None
        hit_lineage_set = set(hit_lineage)
        common = [taxid for taxid in query_lineage if taxid in hit_lineage_set]
        if len(common) == 0:
            return {
                "method": "taxonomy_db",
                "query_taxid": query_taxid,
                "hit_taxid": hit_taxid_val,
                "lca_rank": "none",
                "same_superkingdom": 0,
            }
        lca_taxid = int(common[-1])
        lca_rank = self.rank(lca_taxid).strip().lower() or "no_rank"
        query_super = self.rank_taxid_from_lineage(query_lineage, ["superkingdom", "domain"])
        hit_super = self.rank_taxid_from_lineage(hit_lineage, ["superkingdom", "domain"])
        same_superkingdom = int(query_super > 0 and query_super == hit_super)
        return {
            "method": "taxonomy_db",
            "query_taxid": query_taxid,
            "hit_taxid": hit_taxid_val,
            "lca_rank": lca_rank,
            "same_superkingdom": same_superkingdom,
        }


def compare_species_name_heuristic(query_name: str, hit_name: str) -> Optional[Dict[str, object]]:
    query_norm = normalize_sci_name(query_name)
    hit_norm = normalize_sci_name(hit_name)
    if query_norm == "" or hit_norm == "":
        return None
    qparts = query_norm.lower().split()
    hparts = hit_norm.lower().split()
    if len(qparts) >= 2 and len(hparts) >= 2 and qparts[0] == hparts[0] and qparts[1] == hparts[1]:
        return {
            "method": "name_heuristic",
            "lca_rank": "species",
            "same_superkingdom": pandas.NA,
        }
    if len(qparts) >= 1 and len(hparts) >= 1 and qparts[0] == hparts[0]:
        return {
            "method": "name_heuristic",
            "lca_rank": "genus",
            "same_superkingdom": pandas.NA,
        }
    return {
        "method": "name_heuristic",
        "lca_rank": "genus_mismatch",
        "same_superkingdom": pandas.NA,
    }


def read_contamination_table(path: str) -> pandas.DataFrame:
    try:
        header = pandas.read_csv(path, sep="\t", header=0, nrows=0)
    except Exception:
        return empty_frame(["gene_id", "is_compatible_lineage"])
    gene_col = first_present(header.columns, ["query", "gene_id"])
    compat_col = first_present(header.columns, ["is_compatible_lineage"])
    if gene_col == "" or compat_col == "":
        return empty_frame(["gene_id", "is_compatible_lineage"])
    usecols = [gene_col, compat_col]
    optional_cols = [col for col in ["aligned_taxid", "lca_taxid", "lca_sciname"] if col in header.columns]
    usecols.extend(optional_cols)
    try:
        df = pandas.read_csv(path, sep="\t", header=0, usecols=usecols, low_memory=False)
    except Exception:
        return empty_frame(["gene_id", "is_compatible_lineage"])
    df = df.rename(columns={gene_col: "gene_id", compat_col: "is_compatible_lineage"})
    df["gene_id"] = df["gene_id"].astype(str)
    df["is_compatible_lineage"] = df["is_compatible_lineage"].map(parse_boolish)
    return df.drop_duplicates(subset=["gene_id"], keep="first")


def read_contamination_dir(dir_path: str) -> pandas.DataFrame:
    dir_path = str(dir_path).strip()
    if dir_path == "" or not os.path.isdir(dir_path):
        return empty_frame(["gene_id", "is_compatible_lineage"])
    frames = []
    for name in sorted(os.listdir(dir_path)):
        if name.startswith("."):
            continue
        path = os.path.join(dir_path, name)
        if not os.path.isfile(path):
            continue
        frames.append(read_contamination_table(path))
    if len(frames) == 0:
        return empty_frame(["gene_id", "is_compatible_lineage"])
    df = pandas.concat(frames, ignore_index=True)
    if df.empty:
        return df
    return df.drop_duplicates(subset=["gene_id"], keep="first")


def is_leaf_row(df: pandas.DataFrame) -> pandas.Series:
    if "num_leaf" in df.columns:
        num_leaf = pandas.to_numeric(df["num_leaf"], errors="coerce")
        return num_leaf.eq(1)
    if "so_event" in df.columns:
        return df["so_event"].astype(str).eq("L")
    if "generax_event" in df.columns:
        return df["generax_event"].astype(str).eq("L")
    return pandas.Series(False, index=df.index)


def intron_support_from_leaf_rows(leaf_rows: pandas.DataFrame) -> Dict[str, object]:
    if leaf_rows.empty:
        return {
            "measured_count": 0,
            "supported_count": 0,
            "support_fraction": numpy.nan,
        }
    measured_mask = pandas.Series(False, index=leaf_rows.index)
    supported_mask = pandas.Series(False, index=leaf_rows.index)
    if "num_intron" in leaf_rows.columns:
        num_intron = pandas.to_numeric(leaf_rows["num_intron"], errors="coerce")
        measured_mask = measured_mask | num_intron.notna()
        supported_mask = supported_mask | num_intron.fillna(0).gt(0)
    if "intron_present" in leaf_rows.columns:
        intron_present = pandas.to_numeric(leaf_rows["intron_present"], errors="coerce")
        measured_mask = measured_mask | intron_present.notna()
        supported_mask = supported_mask | intron_present.fillna(0).gt(0)
    measured_count = int(measured_mask.sum())
    supported_count = int((supported_mask & measured_mask).sum())
    support_fraction = supported_count / measured_count if measured_count > 0 else numpy.nan
    return {
        "measured_count": measured_count,
        "supported_count": supported_count,
        "support_fraction": support_fraction,
    }


def expression_support_from_leaf_rows(
    leaf_rows: pandas.DataFrame,
    expression_cols: Sequence[str],
    branch_row: pandas.Series,
) -> Dict[str, object]:
    result = {
        "measured_count": 0,
        "measured_fraction": numpy.nan,
        "clade_corr": safe_float(branch_row.get("clade_min_expression_pearsoncor", numpy.nan)),
    }
    if leaf_rows.empty or len(expression_cols) == 0:
        return result
    expr = leaf_rows.loc[:, list(expression_cols)].apply(pandas.to_numeric, errors="coerce")
    measured_mask = expr.notna().any(axis=1)
    measured_count = int(measured_mask.sum())
    result["measured_count"] = measured_count
    result["measured_fraction"] = measured_count / leaf_rows.shape[0] if leaf_rows.shape[0] > 0 else numpy.nan
    return result


def synteny_support_from_leaf_rows(leaf_rows: pandas.DataFrame) -> Dict[str, object]:
    result = {
        "measured_count": 0,
        "supported_count": 0,
        "support_fraction": numpy.nan,
        "mean_support_score": numpy.nan,
    }
    if leaf_rows.empty:
        return result
    if "synteny_support_score" not in leaf_rows.columns:
        return result
    synteny_score = pandas.to_numeric(leaf_rows["synteny_support_score"], errors="coerce")
    measured_mask = synteny_score.notna()
    supported_mask = synteny_score.fillna(0).gt(0)
    measured_count = int(measured_mask.sum())
    supported_count = int((supported_mask & measured_mask).sum())
    result["measured_count"] = measured_count
    result["supported_count"] = supported_count
    result["support_fraction"] = supported_count / measured_count if measured_count > 0 else numpy.nan
    if measured_count > 0:
        result["mean_support_score"] = float(synteny_score.loc[measured_mask].mean())
    return result


def besthit_support_from_leaf_rows(
    leaf_rows: pandas.DataFrame,
    taxonomy_resolver: TaxonomyResolver,
) -> Dict[str, object]:
    result = {
        "gene_count": 0,
        "taxid_count": 0,
        "method": "",
        "same_superkingdom_fraction": numpy.nan,
        "lca_rank_mode": "",
        "per_gene": {},
    }
    if leaf_rows.empty:
        return result
    hit_taxid_col = first_present(leaf_rows.columns, ["taxid_y", "sprot_taxid", "besthit_taxid", "uniprot_taxid"])
    hit_org_col = first_present(leaf_rows.columns, ["organism", "sprot_organism"])
    hit_acc_col = first_present(leaf_rows.columns, ["sprot_best"])
    taxon_col = first_present(leaf_rows.columns, ["taxon"])
    comparisons = []
    methods = []
    per_gene = {}
    for row in leaf_rows.itertuples(index=False):
        row_dict = row._asdict()
        gene_id = str(row_dict.get("node_name", "")).strip()
        hit_org = str(row_dict.get(hit_org_col, "")).strip() if hit_org_col else ""
        hit_acc = str(row_dict.get(hit_acc_col, "")).strip() if hit_acc_col else ""
        hit_taxid = row_dict.get(hit_taxid_col, "") if hit_taxid_col else ""
        query_name = str(row_dict.get(taxon_col, "")).strip() if taxon_col else ""
        has_hit = (hit_acc != "") or (hit_org != "") or (not pandas.isna(hit_taxid) and str(hit_taxid).strip() != "")
        if not has_hit:
            continue
        result["gene_count"] += 1
        comparison = taxonomy_resolver.compare(query_name, hit_taxid, hit_name=hit_org)
        if comparison is None:
            comparison = compare_species_name_heuristic(query_name, hit_org)
        if comparison is None:
            per_gene[gene_id] = {"besthit_accession": hit_acc, "besthit_organism": hit_org, "besthit_taxid": hit_taxid}
            continue
        if hit_taxid_col and str(hit_taxid).strip() not in {"", "nan"}:
            result["taxid_count"] += 1
        comparisons.append(comparison)
        methods.append(str(comparison.get("method", "")).strip())
        per_gene[gene_id] = {
            "besthit_accession": hit_acc,
            "besthit_organism": hit_org,
            "besthit_taxid": "" if pandas.isna(hit_taxid) else str(hit_taxid),
            "besthit_taxonomy_method": comparison.get("method", ""),
            "besthit_lca_rank": comparison.get("lca_rank", ""),
            "besthit_same_superkingdom": comparison.get("same_superkingdom", pandas.NA),
        }
    result["per_gene"] = per_gene
    if len(comparisons) == 0:
        return result
    lca_ranks = [str(c.get("lca_rank", "")).strip() for c in comparisons if str(c.get("lca_rank", "")).strip() != ""]
    same_super = [safe_float(c.get("same_superkingdom", numpy.nan)) for c in comparisons if not pandas.isna(c.get("same_superkingdom", numpy.nan))]
    if len(methods) > 0:
        result["method"] = Counter(methods).most_common(1)[0][0]
    if len(same_super) > 0:
        result["same_superkingdom_fraction"] = float(numpy.nanmean(same_super))
    if len(lca_ranks) > 0:
        result["lca_rank_mode"] = Counter(lca_ranks).most_common(1)[0][0]
    return result


def contamination_support_for_genes(
    gene_ids: Sequence[str],
    contamination_by_gene: pandas.DataFrame,
) -> Dict[str, object]:
    result = {
        "measured_count": 0,
        "incompatible_count": 0,
        "incompatible_fraction": numpy.nan,
        "top_lca_taxid": pandas.NA,
        "top_lca_sciname": "",
        "top_lca_fraction": numpy.nan,
        "per_gene": {},
    }
    if len(gene_ids) == 0 or contamination_by_gene.empty:
        return result
    subset = contamination_by_gene.loc[contamination_by_gene["gene_id"].isin(gene_ids), :].copy()
    if subset.empty:
        return result
    subset = subset.drop_duplicates(subset=["gene_id"], keep="first")
    measured_mask = subset["is_compatible_lineage"].isin([True, False])
    measured_count = int(measured_mask.sum())
    incompatible_mask = measured_mask & subset["is_compatible_lineage"].eq(False)
    incompatible_count = int(incompatible_mask.sum())
    result["measured_count"] = measured_count
    result["incompatible_count"] = incompatible_count
    result["incompatible_fraction"] = incompatible_count / measured_count if measured_count > 0 else numpy.nan
    if incompatible_count > 0 and "lca_sciname" in subset.columns:
        incompatible_subset = subset.loc[incompatible_mask, :].copy()
        incompatible_subset["lca_taxid_str"] = incompatible_subset.get("lca_taxid", pandas.Series(index=incompatible_subset.index, dtype=object)).fillna("").astype(str)
        incompatible_subset["lca_sciname_str"] = incompatible_subset.get("lca_sciname", pandas.Series(index=incompatible_subset.index, dtype=object)).fillna("").astype(str)
        incompatible_subset = incompatible_subset.loc[
            (incompatible_subset["lca_taxid_str"] != "") | (incompatible_subset["lca_sciname_str"] != ""),
            :
        ]
        if not incompatible_subset.empty:
            top_taxon = (
                incompatible_subset.groupby(["lca_taxid_str", "lca_sciname_str"], sort=False)
                .size()
                .reset_index(name="count")
                .sort_values(["count", "lca_sciname_str", "lca_taxid_str"], ascending=[False, True, True], kind="mergesort")
                .iloc[0]
            )
            result["top_lca_taxid"] = top_taxon["lca_taxid_str"]
            result["top_lca_sciname"] = top_taxon["lca_sciname_str"]
            result["top_lca_fraction"] = float(top_taxon["count"] / incompatible_count)
    per_gene = {}
    for row in subset.itertuples(index=False):
        row_dict = row._asdict()
        gene_id = str(row_dict.get("gene_id", "")).strip()
        per_gene[gene_id] = {
            "contamination_is_compatible_lineage": row_dict.get("is_compatible_lineage", pandas.NA),
            "contamination_lca_taxid": row_dict.get("lca_taxid", pandas.NA),
            "contamination_lca_sciname": row_dict.get("lca_sciname", ""),
        }
    result["per_gene"] = per_gene
    return result


def summarize_candidate_branch(
    branch_row: pandas.Series,
    leaf_rows: pandas.DataFrame,
    contamination_by_gene: pandas.DataFrame,
    expression_cols: Sequence[str],
    taxonomy_resolver: TaxonomyResolver,
) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    candidate_genes = split_gene_labels(branch_row.get("gene_labels", ""))
    if len(candidate_genes) == 0:
        node_name = str(branch_row.get("node_name", "")).strip()
        if node_name != "":
            candidate_genes = [node_name]
    matched_leaf_rows = leaf_rows.loc[leaf_rows["node_name"].isin(candidate_genes), :].copy()
    matched_leaf_rows = matched_leaf_rows.drop_duplicates(subset=["node_name"], keep="first")

    evidence = {
        "besthit": besthit_support_from_leaf_rows(matched_leaf_rows, taxonomy_resolver),
        "intron": intron_support_from_leaf_rows(matched_leaf_rows),
        "expression": expression_support_from_leaf_rows(matched_leaf_rows, expression_cols, branch_row),
        "synteny": synteny_support_from_leaf_rows(matched_leaf_rows),
        "contamination": contamination_support_for_genes(candidate_genes, contamination_by_gene),
    }

    branch_record = {
        "orthogroup": branch_row.get("orthogroup", ""),
        "branch_id": branch_row.get("branch_id", ""),
        "node_name": branch_row.get("node_name", ""),
        "generax_event": branch_row.get("generax_event", ""),
        "generax_transfer": branch_row.get("generax_transfer", ""),
        "generax_event_parent": branch_row.get("generax_event_parent", ""),
        "taxon": branch_row.get("taxon", ""),
        "spnode_coverage": branch_row.get("spnode_coverage", ""),
        "candidate_gene_count": len(candidate_genes),
        "matched_leaf_count": matched_leaf_rows.shape[0],
        "candidate_genes": "; ".join(candidate_genes),
        "besthit_gene_count": evidence["besthit"]["gene_count"],
        "besthit_taxid_count": evidence["besthit"]["taxid_count"],
        "besthit_taxonomy_method": evidence["besthit"]["method"],
        "besthit_same_superkingdom_fraction": evidence["besthit"]["same_superkingdom_fraction"],
        "besthit_lca_rank_mode": evidence["besthit"]["lca_rank_mode"],
        "intron_measured_gene_count": evidence["intron"]["measured_count"],
        "intron_supported_gene_count": evidence["intron"]["supported_count"],
        "intron_support_fraction": evidence["intron"]["support_fraction"],
        "expression_measured_gene_count": evidence["expression"]["measured_count"],
        "expression_measured_fraction": evidence["expression"]["measured_fraction"],
        "clade_min_expression_pearsoncor": evidence["expression"]["clade_corr"],
        "synteny_measured_gene_count": evidence["synteny"]["measured_count"],
        "synteny_supported_gene_count": evidence["synteny"]["supported_count"],
        "synteny_support_fraction": evidence["synteny"]["support_fraction"],
        "synteny_mean_support_score": evidence["synteny"]["mean_support_score"],
        "contamination_measured_gene_count": evidence["contamination"]["measured_count"],
        "contamination_incompatible_gene_count": evidence["contamination"]["incompatible_count"],
        "contamination_incompatible_fraction": evidence["contamination"]["incompatible_fraction"],
        "contamination_top_lca_taxid": evidence["contamination"]["top_lca_taxid"],
        "contamination_top_lca_sciname": evidence["contamination"]["top_lca_sciname"],
        "contamination_top_lca_fraction": evidence["contamination"]["top_lca_fraction"],
    }

    gene_records = []
    besthit_per_gene = evidence["besthit"]["per_gene"]
    contamination_per_gene = evidence["contamination"]["per_gene"]
    intron_supported = {}
    if not matched_leaf_rows.empty:
        intron_df = matched_leaf_rows.loc[:, ["node_name"]].copy()
        intron_flags = intron_support_from_leaf_rows(matched_leaf_rows)
        del intron_flags
        supported_mask = pandas.Series(False, index=matched_leaf_rows.index)
        if "num_intron" in matched_leaf_rows.columns:
            supported_mask = supported_mask | pandas.to_numeric(matched_leaf_rows["num_intron"], errors="coerce").fillna(0).gt(0)
        if "intron_present" in matched_leaf_rows.columns:
            supported_mask = supported_mask | pandas.to_numeric(matched_leaf_rows["intron_present"], errors="coerce").fillna(0).gt(0)
        intron_supported = dict(zip(matched_leaf_rows["node_name"], supported_mask))
    expression_measured = {}
    if not matched_leaf_rows.empty and len(expression_cols) > 0:
        expr_df = matched_leaf_rows.loc[:, ["node_name"] + list(expression_cols)].copy()
        expr_vals = expr_df.loc[:, expression_cols].apply(pandas.to_numeric, errors="coerce")
        expr_mask = expr_vals.notna().any(axis=1)
        expression_measured = dict(zip(expr_df["node_name"], expr_mask))
    synteny_by_gene = {}
    if "synteny_support_score" in matched_leaf_rows.columns:
        synteny_vals = pandas.to_numeric(matched_leaf_rows["synteny_support_score"], errors="coerce")
        synteny_by_gene = dict(zip(matched_leaf_rows["node_name"], synteny_vals))

    for gene_id in candidate_genes:
        leaf_match = matched_leaf_rows.loc[matched_leaf_rows["node_name"] == gene_id, :]
        gene_taxon = ""
        if not leaf_match.empty and "taxon" in leaf_match.columns:
            gene_taxon = str(leaf_match.iloc[0]["taxon"])
        besthit_info = besthit_per_gene.get(gene_id, {})
        contamination_info = contamination_per_gene.get(gene_id, {})
        gene_records.append(
            {
                "orthogroup": branch_row.get("orthogroup", ""),
                "gene_id": gene_id,
                "gene_taxon": gene_taxon,
                "candidate_branch_id": branch_row.get("branch_id", ""),
                "besthit_accession": besthit_info.get("besthit_accession", ""),
                "besthit_organism": besthit_info.get("besthit_organism", ""),
                "besthit_taxid": besthit_info.get("besthit_taxid", ""),
                "besthit_taxonomy_method": besthit_info.get("besthit_taxonomy_method", ""),
                "besthit_lca_rank": besthit_info.get("besthit_lca_rank", ""),
                "besthit_same_superkingdom": besthit_info.get("besthit_same_superkingdom", pandas.NA),
                "intron_supported": bool(intron_supported.get(gene_id, False)),
                "expression_measured": bool(expression_measured.get(gene_id, False)),
                "synteny_support_score": synteny_by_gene.get(gene_id, numpy.nan),
                "contamination_lca_taxid": contamination_info.get("contamination_lca_taxid", pandas.NA),
                "contamination_lca_sciname": contamination_info.get("contamination_lca_sciname", ""),
                "contamination_is_compatible_lineage": contamination_info.get("contamination_is_compatible_lineage", pandas.NA),
            }
        )
    return branch_record, gene_records


def aggregate_gene_records(gene_records: pandas.DataFrame) -> pandas.DataFrame:
    if gene_records.empty:
        return empty_frame(GENE_OUTPUT_COLUMNS)
    gene_records = gene_records.sort_values(
        ["orthogroup", "gene_id", "candidate_branch_id"],
        ascending=[True, True, True],
        kind="mergesort",
    ).reset_index(drop=True)

    rows = []
    for (_orthogroup, _gene_id), group in gene_records.groupby(["orthogroup", "gene_id"], sort=False):
        top = group.iloc[0]
        rows.append(
            {
                "orthogroup": top["orthogroup"],
                "gene_id": top["gene_id"],
                "gene_taxon": top.get("gene_taxon", ""),
                "candidate_branch_count": int(group["candidate_branch_id"].nunique()),
                "candidate_branch_ids": "; ".join(dict.fromkeys(group["candidate_branch_id"].astype(str).tolist())),
                "besthit_accession": top.get("besthit_accession", ""),
                "besthit_organism": top.get("besthit_organism", ""),
                "besthit_taxid": top.get("besthit_taxid", ""),
                "besthit_taxonomy_method": top.get("besthit_taxonomy_method", ""),
                "besthit_lca_rank": top.get("besthit_lca_rank", ""),
                "besthit_same_superkingdom": top.get("besthit_same_superkingdom", pandas.NA),
                "intron_supported": top.get("intron_supported", False),
                "expression_measured": top.get("expression_measured", False),
                "synteny_support_score": top.get("synteny_support_score", numpy.nan),
                "contamination_lca_taxid": top.get("contamination_lca_taxid", pandas.NA),
                "contamination_lca_sciname": top.get("contamination_lca_sciname", ""),
                "contamination_is_compatible_lineage": top.get("contamination_is_compatible_lineage", pandas.NA),
            }
        )
    return pandas.DataFrame(rows, columns=GENE_OUTPUT_COLUMNS)


def aggregate_orthogroup_records(branch_records: pandas.DataFrame, gene_records: pandas.DataFrame) -> pandas.DataFrame:
    if branch_records.empty:
        return empty_frame(ORTHOGROUP_OUTPUT_COLUMNS)
    sorted_branch = branch_records.sort_values(["orthogroup", "branch_id"], ascending=[True, True], kind="mergesort")
    gene_count_by_orthogroup = {}
    if not gene_records.empty:
        gene_count_by_orthogroup = gene_records.groupby("orthogroup")["gene_id"].nunique().to_dict()
    rows = []
    for orthogroup, group in sorted_branch.groupby("orthogroup", sort=False):
        rows.append(
            {
                "orthogroup": orthogroup,
                "hgt_branch_count": int(group["branch_id"].nunique()),
                "hgt_gene_count": int(gene_count_by_orthogroup.get(orthogroup, 0)),
                "candidate_branch_ids": "; ".join(dict.fromkeys(group["branch_id"].astype(str).tolist())),
            }
        )
    return pandas.DataFrame(rows, columns=ORTHOGROUP_OUTPUT_COLUMNS)


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    taxonomy_dbfile = os.environ.get("GG_TAXONOMY_DBFILE", "")
    if args.taxonomy_dbfile is not None:
        taxonomy_dbfile = args.taxonomy_dbfile
    taxonomy_resolver = TaxonomyResolver(taxonomy_dbfile)

    if not os.path.exists(args.dbpath):
        write_tsv(empty_frame(BRANCH_OUTPUT_COLUMNS), args.branch_out, BRANCH_OUTPUT_COLUMNS)
        write_tsv(empty_frame(GENE_OUTPUT_COLUMNS), args.gene_out, GENE_OUTPUT_COLUMNS)
        write_tsv(empty_frame(ORTHOGROUP_OUTPUT_COLUMNS), args.orthogroup_out, ORTHOGROUP_OUTPUT_COLUMNS)
        return

    with sqlite3.connect(args.dbpath) as conn:
        if not sqlite_table_exists(conn, "branch"):
            write_tsv(empty_frame(BRANCH_OUTPUT_COLUMNS), args.branch_out, BRANCH_OUTPUT_COLUMNS)
            write_tsv(empty_frame(GENE_OUTPUT_COLUMNS), args.gene_out, GENE_OUTPUT_COLUMNS)
            write_tsv(empty_frame(ORTHOGROUP_OUTPUT_COLUMNS), args.orthogroup_out, ORTHOGROUP_OUTPUT_COLUMNS)
            return
        branch_columns = sqlite_table_columns(conn, "branch")
        orthogroups = read_candidate_orthogroups(conn, branch_columns)
        if len(orthogroups) == 0:
            write_tsv(empty_frame(BRANCH_OUTPUT_COLUMNS), args.branch_out, BRANCH_OUTPUT_COLUMNS)
            write_tsv(empty_frame(GENE_OUTPUT_COLUMNS), args.gene_out, GENE_OUTPUT_COLUMNS)
            write_tsv(empty_frame(ORTHOGROUP_OUTPUT_COLUMNS), args.orthogroup_out, ORTHOGROUP_OUTPUT_COLUMNS)
            return
        selected_columns = list(dict.fromkeys([
            "orthogroup",
            "branch_id",
            "node_name",
            "gene_labels",
            "num_leaf",
            "so_event",
            "taxon",
            "spnode_coverage",
            "generax_event",
            "generax_transfer",
            "generax_event_parent",
            "clade_min_expression_pearsoncor",
            "num_intron",
            "intron_present",
            "synteny_support_score",
            "sprot_best",
            "organism",
            "taxid_y",
        ] + [col for col in branch_columns if col.startswith("expression_")]))
        selected_columns = [col for col in selected_columns if col in set(branch_columns)]
        branch_df = read_branch_subset(conn, orthogroups, selected_columns)

    if branch_df.empty:
        write_tsv(empty_frame(BRANCH_OUTPUT_COLUMNS), args.branch_out, BRANCH_OUTPUT_COLUMNS)
        write_tsv(empty_frame(GENE_OUTPUT_COLUMNS), args.gene_out, GENE_OUTPUT_COLUMNS)
        write_tsv(empty_frame(ORTHOGROUP_OUTPUT_COLUMNS), args.orthogroup_out, ORTHOGROUP_OUTPUT_COLUMNS)
        return

    contamination_by_gene = read_contamination_dir(args.dir_contamination_tsv)
    expression_cols = [col for col in branch_df.columns if col.startswith("expression_")]

    candidate_mask = pandas.Series(False, index=branch_df.index)
    if "generax_event" in branch_df.columns:
        candidate_mask = candidate_mask | branch_df["generax_event"].astype(str).eq("H")
    if "generax_transfer" in branch_df.columns:
        candidate_mask = candidate_mask | branch_df["generax_transfer"].astype(str).str.startswith("Y")

    leaf_mask = is_leaf_row(branch_df)
    branch_records = []
    gene_records = []
    for orthogroup, og_df in branch_df.groupby("orthogroup", sort=False):
        candidate_rows = og_df.loc[candidate_mask.loc[og_df.index], :].copy()
        if candidate_rows.empty:
            continue
        leaf_rows = og_df.loc[leaf_mask.loc[og_df.index], :].copy()
        leaf_rows["node_name"] = leaf_rows["node_name"].astype(str)
        leaf_rows = leaf_rows.drop_duplicates(subset=["node_name"], keep="first")
        for _, branch_row in candidate_rows.iterrows():
            branch_record, branch_gene_records = summarize_candidate_branch(
                branch_row=branch_row,
                leaf_rows=leaf_rows,
                contamination_by_gene=contamination_by_gene,
                expression_cols=expression_cols,
                taxonomy_resolver=taxonomy_resolver,
            )
            branch_records.append(branch_record)
            gene_records.extend(branch_gene_records)

    branch_out = pandas.DataFrame(branch_records, columns=BRANCH_OUTPUT_COLUMNS)
    if not branch_out.empty:
        branch_out = branch_out.sort_values(
            ["orthogroup", "branch_id"],
            ascending=[True, True],
            kind="mergesort",
        ).reset_index(drop=True)
    else:
        branch_out = empty_frame(BRANCH_OUTPUT_COLUMNS)

    gene_out = aggregate_gene_records(pandas.DataFrame(gene_records))
    orthogroup_out = aggregate_orthogroup_records(branch_out, gene_out)

    write_tsv(branch_out, args.branch_out, BRANCH_OUTPUT_COLUMNS)
    write_tsv(gene_out, args.gene_out, GENE_OUTPUT_COLUMNS)
    write_tsv(orthogroup_out, args.orthogroup_out, ORTHOGROUP_OUTPUT_COLUMNS)


if __name__ == "__main__":
    main()
