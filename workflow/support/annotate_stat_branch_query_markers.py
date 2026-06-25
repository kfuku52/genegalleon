#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path

import pandas


MARKER_COLUMNS = [
    "query_marker",
    "query_marker_source",
    "query_marker_best_evalue",
    "query_marker_best_qjointcov",
    "query_marker_best_bitscore",
]
BEST_HIT_MARKER = "Best hit"


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stat_branch", metavar="PATH", required=True)
    parser.add_argument("--query_gene", metavar="PATH", required=True)
    parser.add_argument("--query_aa_fasta", metavar="PATH", default="")
    parser.add_argument("--query_blast", metavar="PATH", required=True)
    parser.add_argument("--outfile", metavar="PATH", required=True)
    parser.add_argument("--min_query_blast_coverage", metavar="FLOAT", default=0.25, type=float)
    return parser


def _open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def normalize_id(value):
    return str(value).strip().replace("−", "-")


def read_fasta_headers(path):
    headers = []
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return headers
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]
                if header:
                    headers.append(header)
    return headers


def read_query_gene_ids(path):
    if not Path(path).exists() or Path(path).stat().st_size == 0:
        return [], False
    with _open_text(path) as handle:
        first = handle.read(1)
        handle.seek(0)
        if first == ">":
            return read_fasta_headers(path), True
        ids = []
        for line in handle:
            value = line.strip()
            if value:
                ids.append(value)
        return ids, False


def unique_preserve_order(values):
    out = []
    seen = set()
    for value in values:
        value = str(value)
        if value and value not in seen:
            out.append(value)
            seen.add(value)
    return out


def format_float(value):
    if pandas.isna(value):
        return ""
    return f"{float(value):.6g}"


def parse_semicolon_min_float(value):
    if pandas.isna(value):
        return pandas.NA
    values = []
    for token in str(value).split(";"):
        token = token.strip()
        if not token:
            continue
        try:
            values.append(float(token))
        except ValueError:
            continue
    if not values:
        return pandas.NA
    return min(values)


def parse_semicolon_max_float(value):
    if pandas.isna(value):
        return pandas.NA
    values = []
    for token in str(value).split(";"):
        token = token.strip()
        if not token:
            continue
        try:
            values.append(float(token))
        except ValueError:
            continue
    if not values:
        return pandas.NA
    return max(values)


def node_matches_query(node_name, query_id, allow_suffix_match=False):
    node = normalize_id(node_name)
    query = normalize_id(query_id)
    if not node or not query:
        return False
    if node == query:
        return True
    if allow_suffix_match:
        return any(node.endswith(sep + query) for sep in ("_", "-", "."))
    return False


def direct_query_sources_by_node(tip_names, query_gene_path, query_aa_fasta_path=""):
    query_ids, query_gene_is_fasta = read_query_gene_ids(query_gene_path)
    query_aa_headers = read_fasta_headers(query_aa_fasta_path)
    exact_ids = unique_preserve_order([normalize_id(x) for x in query_ids + query_aa_headers])
    gene_list_ids = unique_preserve_order([normalize_id(x) for x in query_ids]) if not query_gene_is_fasta else []

    direct = {}
    for tip in tip_names:
        sources = []
        for query_id in exact_ids:
            if node_matches_query(tip, query_id, allow_suffix_match=False):
                sources.append(query_id)
        for query_id in gene_list_ids:
            if node_matches_query(tip, query_id, allow_suffix_match=True):
                sources.append(query_id)
        sources = unique_preserve_order(sources)
        if sources:
            direct[tip] = sources
    return direct


def read_query_blast(path):
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return pandas.DataFrame()
    df = pandas.read_csv(path, sep="\t", header=0, dtype=str)
    if df.empty:
        return df
    required = {"qacc", "sacc", "qjointcov"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"query_blast table is missing required columns: {', '.join(sorted(missing))}")

    df = df.copy()
    df["qacc"] = df["qacc"].astype(str).str.strip()
    df["sacc"] = df["sacc"].astype(str).str.strip()
    df["qjointcov_num"] = pandas.to_numeric(df["qjointcov"], errors="coerce")
    if "min_evalue" in df.columns:
        df["evalue_num"] = pandas.to_numeric(df["min_evalue"], errors="coerce")
    elif "evalue" in df.columns:
        df["evalue_num"] = df["evalue"].map(parse_semicolon_min_float)
    else:
        df["evalue_num"] = pandas.NA
    if "max_bitscore" in df.columns:
        df["bitscore_num"] = pandas.to_numeric(df["max_bitscore"], errors="coerce")
    elif "bitscore" in df.columns:
        df["bitscore_num"] = df["bitscore"].map(parse_semicolon_max_float)
    else:
        df["bitscore_num"] = pandas.NA
    return df


def best_blast_sources_by_node(tip_names, query_blast_path, min_query_blast_coverage=0.25):
    df = read_query_blast(query_blast_path)
    if df.empty:
        return {}

    tip_by_norm = {normalize_id(tip): tip for tip in tip_names}
    df = df.loc[
        (df["qacc"] != "")
        & (df["sacc"] != "")
        & (df["qjointcov_num"] >= float(min_query_blast_coverage)),
        :,
    ].copy()
    if df.empty:
        return {}

    df["sacc_norm"] = df["sacc"].map(normalize_id)
    df = df.loc[df["sacc_norm"].isin(tip_by_norm), :].copy()
    if df.empty:
        return {}
    df["node_name"] = df["sacc_norm"].map(tip_by_norm)

    df["evalue_sort"] = pandas.to_numeric(df["evalue_num"], errors="coerce").fillna(float("inf"))
    df["bitscore_sort"] = pandas.to_numeric(df["bitscore_num"], errors="coerce").fillna(float("-inf"))
    df = df.sort_values(
        by=["qacc", "qjointcov_num", "evalue_sort", "bitscore_sort"],
        ascending=[True, False, True, False],
        kind="mergesort",
    )
    best = df.drop_duplicates(subset=["qacc"], keep="first").copy()

    out = {}
    for node_name, group in best.groupby("node_name", sort=False):
        group = group.sort_values(
            by=["qacc", "qjointcov_num", "evalue_sort", "bitscore_sort"],
            ascending=[True, False, True, False],
            kind="mergesort",
        )
        out[node_name] = {
            "query_ids": unique_preserve_order(group["qacc"].astype(str).tolist()),
            "evalues": [format_float(x) for x in group["evalue_num"].tolist()],
            "qjointcovs": [format_float(x) for x in group["qjointcov_num"].tolist()],
            "bitscores": [format_float(x) for x in group["bitscore_num"].tolist()],
        }
    return out


def annotate_stat_branch(
    stat_branch_path,
    query_gene_path,
    query_blast_path,
    outfile,
    query_aa_fasta_path="",
    min_query_blast_coverage=0.25,
):
    df = pandas.read_csv(stat_branch_path, sep="\t", header=0, dtype=str, keep_default_na=False)
    for column in MARKER_COLUMNS:
        if column in df.columns:
            df = df.drop(columns=[column])
    if "node_name" not in df.columns:
        raise ValueError("stat_branch table is missing required column: node_name")

    is_tip = df["so_event"].astype(str).eq("L") if "so_event" in df.columns else pandas.Series([True] * len(df), index=df.index)
    tip_names = df.loc[is_tip, "node_name"].astype(str).tolist()
    direct_by_node = direct_query_sources_by_node(tip_names, query_gene_path, query_aa_fasta_path)
    best_by_node = best_blast_sources_by_node(tip_names, query_blast_path, min_query_blast_coverage)

    values = {column: [""] * len(df) for column in MARKER_COLUMNS}
    values["query_marker"] = ["-"] * len(df)
    for idx, row in df.iterrows():
        if not bool(is_tip.loc[idx]):
            continue
        node_name = str(row["node_name"])
        direct_sources = direct_by_node.get(node_name, [])
        best = best_by_node.get(node_name, {})
        best_sources = best.get("query_ids", [])
        marker = BEST_HIT_MARKER if best_sources else "-"
        values["query_marker"][idx] = marker
        source_tokens = []
        if direct_sources:
            source_tokens.append("direct:" + ";".join(direct_sources))
        if best_sources:
            source_tokens.append("best:" + ";".join(best_sources))
        values["query_marker_source"][idx] = "|".join(source_tokens)
        values["query_marker_best_evalue"][idx] = ";".join(best.get("evalues", []))
        values["query_marker_best_qjointcov"][idx] = ";".join(best.get("qjointcovs", []))
        values["query_marker_best_bitscore"][idx] = ";".join(best.get("bitscores", []))

    for column in MARKER_COLUMNS:
        df[column] = values[column]

    df.to_csv(outfile, sep="\t", index=False)
    marked = sum(value != "-" for value in values["query_marker"])
    print(f"Query marker annotation completed: marked_tips={marked}, output={outfile}", flush=True)


def main():
    args = build_arg_parser().parse_args()
    annotate_stat_branch(
        stat_branch_path=args.stat_branch,
        query_gene_path=args.query_gene,
        query_aa_fasta_path=args.query_aa_fasta,
        query_blast_path=args.query_blast,
        outfile=args.outfile,
        min_query_blast_coverage=args.min_query_blast_coverage,
    )


if __name__ == "__main__":
    main()
