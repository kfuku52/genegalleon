#!/usr/bin/env python3
# coding: utf-8

import argparse
import datetime
import gzip
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas

try:
    from busco_table_utils import SOURCE_COORDINATE_PATTERN
    from cds_annotation_scope import QC_COLUMNS, canonical, load_scope
    from line_comment_stream import LineCommentStream
except ImportError:  # package import used by tests and library callers
    from workflow.support.busco_table_utils import SOURCE_COORDINATE_PATTERN
    from workflow.support.cds_annotation_scope import QC_COLUMNS, canonical, load_scope
    from workflow.support.line_comment_stream import LineCommentStream


def identifier_text(value):
    """Keep numeric-looking and NA-like identifiers literal, while empty stays missing."""
    return value if value else float("nan")


def read_fasta_ids(path):
    ids = []
    opener = gzip.open if path.endswith(".gz") else open
    mode = "rt" if path.endswith(".gz") else "r"
    with opener(path, mode) as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def dedupe_gene_id_index(tmp):
    if tmp is None:
        return None
    if "gene_id" not in tmp.columns:
        return tmp
    return tmp.drop_duplicates(subset=["gene_id"], keep="first").set_index("gene_id")


def load_orthogroup_map(path, species_col):
    if not os.path.exists(path):
        print("File not found: {}".format(path), flush=True)
        return None
    print("{}: Processing: {}".format(datetime.datetime.now(), path), flush=True)

    header = pandas.read_csv(path, sep="\t", nrows=0)
    if ("Orthogroup" not in header.columns) or (species_col not in header.columns):
        print("Required columns were not found in orthogroup file: Orthogroup, {}".format(species_col), flush=True)
        return None

    tmp = pandas.read_csv(path, sep="\t", usecols=["Orthogroup", species_col], low_memory=False,
                          converters={"Orthogroup": identifier_text, species_col: identifier_text})
    tmp = tmp.dropna(subset=[species_col]).copy()
    if tmp.empty:
        return None

    # Preserve the first orthogroup for each member in source order without
    # materializing exploded copies of the input columns.
    member_groups = {}
    for group, members in zip(tmp["Orthogroup"].to_numpy(), tmp[species_col].to_numpy(), strict=True):
        for member in str(members).split(","):
            member = member.strip()
            if member:
                member_groups.setdefault(member, group)
    return pandas.Series(
        list(member_groups.values()),
        index=pandas.Index(list(member_groups), name="gene_id", dtype=tmp[species_col].astype(str).dtype),
        name="Orthogroup",
        dtype=tmp["Orthogroup"].dtype,
    )


def load_uniprot(path):
    if not os.path.exists(path):
        print("File not found: {}".format(path), flush=True)
        return None
    print("{}: Processing: {}".format(datetime.datetime.now(), path), flush=True)
    return dedupe_gene_id_index(pandas.read_csv(path, sep="\t", header=0, index_col=None, low_memory=False,
                                               converters={"gene_id": identifier_text}))


def load_cdskit_localize(path):
    if not os.path.exists(path):
        print("File not found: {}".format(path), flush=True)
        return None
    print("{}: Processing: {}".format(datetime.datetime.now(), path), flush=True)
    tmp = pandas.read_csv(path, sep="\t", header=0, index_col=None, low_memory=False,
                          converters={"seq_id": identifier_text})
    if "seq_id" not in tmp.columns:
        print("seq_id column was not found in cdskit localize table: {}".format(path), flush=True)
        return None
    tmp = tmp.rename(columns={"seq_id": "gene_id"})
    rename_map = {col: "cdskit_localize_{}".format(col) for col in tmp.columns if col != "gene_id"}
    tmp = tmp.rename(columns=rename_map)
    return dedupe_gene_id_index(tmp)


def load_busco(path):
    if not os.path.exists(path):
        print("File not found: {}".format(path), flush=True)
        return None
    print("{}: Processing: {}".format(datetime.datetime.now(), path), flush=True)
    colnames = [
        "busco_id",
        "busco_status",
        "busco_sequence",
        "busco_score",
        "busco_length",
        "busco_orthodb_url",
        "busco_description",
    ]
    with LineCommentStream(path, quoted=True) as source:
        tmp = pandas.read_csv(source, sep="\t", header=None, index_col=None, low_memory=False, names=colnames,
                              converters={"busco_sequence": identifier_text})
    tmp = tmp.dropna(subset=["busco_sequence"]).copy()
    tmp["gene_id"] = tmp["busco_sequence"].astype(str).str.replace(SOURCE_COORDINATE_PATTERN, "", regex=True)
    tmp = tmp.loc[tmp["gene_id"].notna() & (tmp["gene_id"] != ""), :].copy()
    if tmp.empty:
        return None
    agg_cols = [c for c in tmp.columns if c != "gene_id"]
    # Make the serialized missing-value spelling explicit before string joins;
    # pandas string dtype can retain a missing scalar after astype(str).
    tmp_agg = tmp.loc[:, ["gene_id"] + agg_cols].astype(str).fillna("nan")
    # Most BUSCO genes have one row. Only repeated genes need string joins.
    repeated = tmp_agg["gene_id"].duplicated(keep=False)
    if not repeated.any():
        return dedupe_gene_id_index(tmp_agg)
    # Reuse row positions for every annotation column instead of allocating a
    # Series and invoking an aggregation callback for every gene and column.
    duplicate_rows = tmp_agg if repeated.all() else tmp_agg.loc[repeated]
    grouped = duplicate_rows.groupby("gene_id", sort=False).indices
    duplicates = pandas.DataFrame({
        "gene_id": pandas.Series(list(grouped), dtype=tmp_agg["gene_id"].dtype),
    })
    for column in agg_cols:
        values = duplicate_rows[column].to_numpy(dtype=object)
        duplicates[column] = pandas.Series(
            ["; ".join(values.take(indices)) for indices in grouped.values()],
            dtype=tmp_agg[column].dtype,
        )
    if repeated.all():
        return dedupe_gene_id_index(duplicates)
    gene_order = tmp_agg["gene_id"].drop_duplicates()
    tmp_agg = pandas.concat([tmp_agg.loc[~repeated], duplicates], ignore_index=True)
    return dedupe_gene_id_index(tmp_agg).reindex(gene_order)


def load_fx2tab(path):
    if not os.path.exists(path):
        print("File not found: {}".format(path), flush=True)
        return None
    print("{}: Processing: {}".format(datetime.datetime.now(), path), flush=True)
    tmp = pandas.read_csv(path, sep="\t", header=0, index_col=None, low_memory=False,
                          converters={"#id": identifier_text, "gene_id": identifier_text})
    tmp = tmp.rename(columns={"length": "cds_length", "#id": "gene_id"})
    return dedupe_gene_id_index(tmp)


def load_gff_info(path):
    if not os.path.exists(path):
        print("File not found: {}".format(path), flush=True)
        return None
    print("{}: Processing: {}".format(datetime.datetime.now(), path), flush=True)
    return dedupe_gene_id_index(pandas.read_csv(path, sep="\t", header=0, index_col=None, low_memory=False,
                                               converters={"gene_id": identifier_text}))


def load_expression(path):
    if not os.path.exists(path):
        print("File not found: {}".format(path), flush=True)
        return None
    print("{}: Processing: {}".format(datetime.datetime.now(), path), flush=True)
    tmp = pandas.read_csv(path, sep="\t", header=0, index_col=None, low_memory=False,
                          converters={0: identifier_text})
    if "Unnamed: 0" in tmp.columns:
        tmp = tmp.rename(columns={"Unnamed: 0": "gene_id"})
    else:
        tmp = tmp.rename(columns={tmp.columns[0]: "gene_id"})
    return dedupe_gene_id_index(tmp)


def load_mmseqs(path):
    if not os.path.exists(path):
        print("File not found: {}".format(path), flush=True)
        return None
    print("{}: Processing: {}".format(datetime.datetime.now(), path), flush=True)
    tmp = pandas.read_csv(path, sep="\t", header=None, index_col=None, low_memory=False,
                          converters={0: identifier_text})
    tmp.columns = [
        "gene_id",
        "lca_taxid",
        "lca_rank",
        "lca_sciname",
        "num_assigned_protein_fragment",
        "num_labeled_protein_fragment",
        "num_taxid_supporting_protein_fragment",
        "fraction_evalue_support",
        "lineage_taxids",
    ]  # https://github.com/soedinglab/MMseqs2/wiki#taxonomy-format
    return dedupe_gene_id_index(tmp)


def join_if_available(df, tmp):
    if tmp is None:
        return df
    if isinstance(tmp, pandas.Series):
        return df.join(tmp, how="left")
    if getattr(tmp.index, "name", None) != "gene_id":
        if "gene_id" not in tmp.columns:
            return df
        tmp = dedupe_gene_id_index(tmp)
    return df.join(tmp, how="left")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--scientific_name", metavar="PATH", default="", type=str, help="Path used by --scientific_name."
    )
    parser.add_argument("--cds_fasta", metavar="PATH", default="", type=str, help="Path used by --cds_fasta.")
    parser.add_argument("--resolved_cds_fasta", default="", help="QC-accepted CDS FASTA, separate from the original catalogue.")
    parser.add_argument("--cds_resolution_json", default="", help="Existing CDS-resolution decisions; does not rerun QC.")
    parser.add_argument("--scope_json", default="", help="Write the separate catalogue/eligible-sequence contract.")
    parser.add_argument("--existing_annotation_tsv", default="", help="Reuse exact accepted annotation fields without rerunning predictors.")
    parser.add_argument("--orthogroup_tsv", metavar="PATH", default="", type=str, help="Path used by --orthogroup_tsv.")
    parser.add_argument("--uniprot_tsv", metavar="PATH", default="", type=str, help="Path used by --uniprot_tsv.")
    parser.add_argument(
        "--cdskit_localize_tsv", metavar="PATH", default="", type=str, help="Path used by --cdskit_localize_tsv."
    )
    parser.add_argument("--busco_tsv", metavar="PATH", default="", type=str, help="Path used by --busco_tsv.")
    parser.add_argument("--expression_tsv", metavar="PATH", default="", type=str, help="Path used by --expression_tsv.")
    parser.add_argument("--gff_info", metavar="PATH", default="", type=str, help="Path used by --gff_info.")
    parser.add_argument(
        "--mmseqs2taxonomy_tsv", metavar="PATH", default="", type=str, help="Path used by --mmseqs2taxonomy_tsv."
    )
    parser.add_argument("--fx2tab", metavar="PATH", default="", type=str, help="Path used by --fx2tab.")
    parser.add_argument(
        "--out_tsv", metavar="PATH", default="cds_annotation.tsv", type=str, help="Path used by --out_tsv."
    )
    parser.add_argument("--ncpu", metavar="INT", default=1, type=int, help="Number of worker threads.")
    args = parser.parse_args()
    args.ncpu = max(1, int(args.ncpu))

    scoped = bool(args.resolved_cds_fasta or args.cds_resolution_json or args.scope_json)
    if scoped and not all((args.resolved_cds_fasta, args.cds_resolution_json, args.scope_json)):
        parser.error("--resolved_cds_fasta, --cds_resolution_json and --scope_json are required together")
    if args.existing_annotation_tsv and not scoped:
        parser.error("Saved annotation reuse requires a validated CDS-resolution scope")
    seqnames = read_fasta_ids(args.cds_fasta)
    if len(seqnames) != len(set(seqnames)):
        raise ValueError("Original CDS catalogue contains duplicate IDs")
    accepted, qc_rows, scope = set(seqnames), {}, None
    if scoped:
        seqnames, accepted, qc_rows, scope = load_scope(args.cds_fasta, args.resolved_cds_fasta, args.cds_resolution_json)
    df = pandas.DataFrame(index=pandas.Index(seqnames, name="gene_id"))

    orthogroup_map = load_orthogroup_map(args.orthogroup_tsv, args.scientific_name)
    if orthogroup_map is not None:
        df.loc[:, "orthogroup"] = df.index.to_series().map(orthogroup_map).fillna("")

    tasks = {
        "uniprot": (load_uniprot, args.uniprot_tsv),
        "cdskit_localize": (load_cdskit_localize, args.cdskit_localize_tsv),
        "busco": (load_busco, args.busco_tsv),
        "fx2tab": (load_fx2tab, args.fx2tab),
        "gff_info": (load_gff_info, args.gff_info),
        "expression": (load_expression, args.expression_tsv),
        "mmseqs": (load_mmseqs, args.mmseqs2taxonomy_tsv),
    }

    if args.existing_annotation_tsv:
        if any(path for _, path in tasks.values()):
            parser.error("Saved annotation reuse cannot be mixed with prediction provider tables")
        saved = pandas.read_csv(args.existing_annotation_tsv, sep="\t", dtype=str, na_filter=False)
        if ("gene_id" not in saved or saved.gene_id.duplicated().any()
                or set(saved.gene_id) != accepted or set(QC_COLUMNS) & set(saved.columns)):
            raise ValueError("Saved annotation must represent exactly the accepted CDS IDs, without prior QC columns")
        saved = saved.set_index("gene_id")
        if "orthogroup" in df and "orthogroup" in saved:
            if not df.loc[saved.index, "orthogroup"].equals(saved["orthogroup"]):
                raise ValueError("Saved annotation orthogroups disagree with the admitted catalogue")
            saved = saved.drop(columns="orthogroup")
        df = df.join(saved, how="left").fillna("")

    loaded = {k: None for k in tasks}
    if args.existing_annotation_tsv:
        pass
    elif args.ncpu == 1:
        for key, (func, path) in tasks.items():
            loaded[key] = func(path)
    else:
        max_workers = min(args.ncpu, len(tasks))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(func, path): key for key, (func, path) in tasks.items()}
            for future in as_completed(futures):
                key = futures[future]
                loaded[key] = future.result()

    for key in ["uniprot", "cdskit_localize", "busco", "fx2tab", "gff_info", "expression", "mmseqs"]:
        df = join_if_available(df, loaded[key])

    if scoped:
        # Predictions remain restricted to resolved CDS; excluded catalogue rows
        # retain their identity and evidence, without borrowing functional hits.
        excluded = ~df.index.isin(accepted)
        df.loc[excluded, df.columns.difference(["orthogroup"])] = ""
        for column in QC_COLUMNS:
            df[column] = [qc_rows[gene][column] for gene in df.index]

    df = df.reset_index()
    df.to_csv(args.out_tsv, sep="\t", index=False)
    if scope is not None:
        Path(args.scope_json).write_text(canonical(scope) + "\n")
    print("merge_cds_annotation.py completed!", flush=True)


if __name__ == "__main__":
    main()
