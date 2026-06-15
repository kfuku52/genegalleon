#!/usr/bin/env python3
# coding: utf-8

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy
import pandas

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from species_labeling import scientific_name_from_label


STAT_BRANCH_SUFFIX = "_stat.branch.tsv"


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("query2family", "orthogroup"), required=True)
    parser.add_argument("--dir_gene_family", metavar="PATH", type=str, required=True)
    parser.add_argument("--dir_query_gene", metavar="PATH", type=str, default="")
    parser.add_argument("--orthogroup_genecount", metavar="PATH", type=str, default="")
    parser.add_argument("--species_tree", metavar="PATH", type=str, required=True)
    parser.add_argument("--out_presence", metavar="PATH", type=str, required=True)
    parser.add_argument("--out_copy_number", metavar="PATH", type=str, required=True)
    parser.add_argument("--out_long", metavar="PATH", type=str, required=True)
    parser.add_argument("--out_plot_presence", metavar="PATH", type=str, required=True)
    parser.add_argument("--out_plot_copy_number", metavar="PATH", type=str, required=True)
    parser.add_argument("--out_plot_long", metavar="PATH", type=str, required=True)
    parser.add_argument("--out_selection", metavar="PATH", type=str, required=True)
    parser.add_argument("--include_incomplete", metavar="0|1", default=0, type=int)
    parser.add_argument("--max_families", metavar="N|auto|0", default="auto")
    parser.add_argument("--family_ids", metavar="ID[,ID...]", default="")
    parser.add_argument("--family_file", metavar="PATH", default="")
    return parser


def visible_query_ids(path):
    query_dir = Path(path)
    if not query_dir.is_dir():
        raise FileNotFoundError(f"Input query_gene directory was not found: {query_dir}")
    return sorted(entry.name for entry in query_dir.iterdir() if entry.is_file() and not entry.name.startswith("."))


def read_orthogroup_ids_from_genecount(path):
    genecount_path = Path(path)
    if not genecount_path.is_file() or genecount_path.stat().st_size == 0:
        return []
    with genecount_path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        family_col = None
        for candidate in ("Orthogroup", "HOG", "OG", "family_id", "query"):
            if candidate in fieldnames:
                family_col = candidate
                break
        if family_col is None and fieldnames:
            family_col = fieldnames[0]
        if family_col is None:
            return []
        return [str(row.get(family_col) or "").strip() for row in reader if str(row.get(family_col) or "").strip()]


def stat_branch_family_id(path):
    name = Path(path).name
    if name.endswith(STAT_BRANCH_SUFFIX):
        return name[: -len(STAT_BRANCH_SUFFIX)]
    return Path(path).stem


def visible_stat_branch_ids(dir_gene_family):
    stat_branch_dir = Path(dir_gene_family) / "stat_branch"
    if not stat_branch_dir.is_dir():
        raise FileNotFoundError(f"stat_branch directory was not found: {stat_branch_dir}")
    return sorted(
        stat_branch_family_id(path)
        for path in stat_branch_dir.glob(f"*{STAT_BRANCH_SUFFIX}")
        if path.is_file() and not path.name.startswith(".")
    )


def resolve_family_ids(mode, dir_gene_family, dir_query_gene="", orthogroup_genecount=""):
    if mode == "query2family":
        family_ids = visible_query_ids(dir_query_gene)
        if len(family_ids) == 0:
            raise ValueError(f"Input query_gene directory is empty: {dir_query_gene}")
        return family_ids, "query_gene"

    family_ids = read_orthogroup_ids_from_genecount(orthogroup_genecount)
    if family_ids:
        return family_ids, "orthogroup_genecount"

    family_ids = visible_stat_branch_ids(dir_gene_family)
    if family_ids:
        return family_ids, "stat_branch"
    raise ValueError(f"No orthogroup IDs were found from genecount or stat_branch under: {dir_gene_family}")


def read_newick_tip_labels(path):
    tree_path = Path(path)
    if not tree_path.is_file() or tree_path.stat().st_size == 0:
        raise FileNotFoundError(f"Species tree was not found or empty: {tree_path}")
    text = tree_path.read_text(encoding="utf-8", errors="replace").strip()
    tips = []
    i = 0
    n = len(text)

    def read_label(start):
        if start >= n:
            return "", start
        if text[start] == "'":
            out = []
            j = start + 1
            while j < n:
                ch = text[j]
                if ch == "'":
                    if (j + 1) < n and text[j + 1] == "'":
                        out.append("'")
                        j += 2
                        continue
                    return "".join(out), j + 1
                out.append(ch)
                j += 1
            return "".join(out), j
        j = start
        while j < n and text[j] not in ":,();":
            j += 1
        return text[start:j].strip(), j

    while i < n:
        ch = text[i]
        if ch not in "(,":
            i += 1
            continue
        i += 1
        while i < n and text[i].isspace():
            i += 1
        if i >= n or text[i] in "(),;":
            continue
        label, i = read_label(i)
        label = label.strip()
        if label != "":
            tips.append(label)
    return tips


def read_stat_branch_copy_numbers(path):
    counts = {}
    stat_path = Path(path)
    if not stat_path.is_file() or stat_path.stat().st_size == 0:
        return counts

    with stat_path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if "spnode_coverage" not in fieldnames:
            raise ValueError(f"stat_branch is missing spnode_coverage: {stat_path}")
        for row in reader:
            is_tip = False
            if row.get("so_event") == "L":
                is_tip = True
            elif row.get("num_leaf") == "1":
                is_tip = True
            if not is_tip:
                continue
            species = str(row.get("spnode_coverage") or "").strip()
            if species == "":
                continue
            counts[species] = counts.get(species, 0) + 1
    return counts


def build_species_aliases(species_order):
    aliases = {}
    for species in species_order:
        aliases[species] = species
        aliases[species.replace(" ", "_")] = species
        aliases[species.replace("_", " ")] = species
    return aliases


def resolve_counts_to_species_order(raw_counts, species_order):
    aliases = build_species_aliases(species_order)
    resolved = {species: 0 for species in species_order}
    unresolved = {}
    for raw_species, count in raw_counts.items():
        species = aliases.get(raw_species)
        if species is None:
            unresolved[raw_species] = count
            continue
        resolved[species] += int(count)
    return resolved, unresolved


def write_matrix(path, matrix):
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(out_path, sep="\t", index=True, na_rep="NA")


def split_family_ids(text):
    if text is None or str(text).strip() == "":
        return []
    return [item for item in re.split(r"[\s,]+", str(text).strip()) if item]


def read_family_file(path):
    if path is None or str(path).strip() == "":
        return []
    family_file = Path(path)
    if not family_file.is_file() or family_file.stat().st_size == 0:
        raise FileNotFoundError(f"gene_summary_family_file was not found or empty: {family_file}")
    out = []
    header_tokens = {"family_id", "query", "orthogroup", "gene_family_id", "id"}
    with family_file.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped == "" or stripped.startswith("#"):
                continue
            token = re.split(r"[\t, ]+", stripped, maxsplit=1)[0]
            if not out and token.lower() in header_tokens:
                continue
            out.append(token)
    return out


def deduplicate_ordered(values):
    seen = set()
    out = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def parse_max_families(value, mode):
    value = str(value or "auto").strip().lower()
    if value in ("", "auto"):
        return 100 if mode == "orthogroup" else 0
    if value in ("all", "none", "no", "unlimited"):
        return 0
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"--max_families must be auto, all, or a non-negative integer: {value}") from exc
    if parsed < 0:
        raise ValueError(f"--max_families must be non-negative: {value}")
    return parsed


def select_plot_families(available_ids, query_status, mode, max_families, family_ids="", family_file=""):
    explicit_ids = deduplicate_ordered(split_family_ids(family_ids) + read_family_file(family_file))
    available_set = set(available_ids)
    if explicit_ids:
        missing = [family_id for family_id in explicit_ids if family_id not in available_set]
        if missing:
            raise ValueError(
                "Requested gene family ID(s) were not available in the matrix: "
                + ", ".join(missing[:20])
                + (" ..." if len(missing) > 20 else "")
            )
        selected_ids = explicit_ids
        source = "explicit"
    else:
        max_count = parse_max_families(max_families, mode)
        selected_ids = list(available_ids if max_count == 0 else available_ids[:max_count])
        source = "all" if max_count == 0 else f"first_{max_count}"

    base_order = {family_id: idx + 1 for idx, family_id in enumerate(available_ids)}
    rows = []
    for idx, family_id in enumerate(selected_ids, start=1):
        rows.append(
            {
                "family_id": family_id,
                "plot_order": idx,
                "base_order": base_order.get(family_id, pandas.NA),
                "status": query_status.get(family_id, "unknown"),
                "selection_source": source,
            }
        )
    return selected_ids, pandas.DataFrame(rows)


def build_long_table(copy_number, presence, species_order, family_order, query_status, mode):
    long_rows = []
    query_order = {family_id: idx + 1 for idx, family_id in enumerate(family_order)}
    species_order_map = {species: idx + 1 for idx, species in enumerate(species_order)}
    for family_id in copy_number.columns:
        status = query_status.get(family_id, "unknown")
        for species in species_order:
            copy_value = copy_number.loc[species, family_id]
            presence_value = presence.loc[species, family_id]
            long_rows.append(
                {
                    "species": species,
                    "species_display": scientific_name_from_label(species),
                    "species_order": species_order_map[species],
                    "query": family_id,
                    "query_order": query_order[family_id],
                    "family_id": family_id,
                    "family_order": query_order[family_id],
                    "mode": mode,
                    "copy_number": copy_value,
                    "presence": presence_value,
                    "status": status,
                }
            )
    return pandas.DataFrame(long_rows)


def build_matrices(mode, dir_gene_family, species_tree, dir_query_gene="", orthogroup_genecount="", include_incomplete=False):
    family_ids, family_source = resolve_family_ids(
        mode=mode,
        dir_gene_family=dir_gene_family,
        dir_query_gene=dir_query_gene,
        orthogroup_genecount=orthogroup_genecount,
    )

    species_order = read_newick_tip_labels(species_tree)
    if len(species_order) == 0:
        raise ValueError(f"No tip labels were detected in species tree: {species_tree}")

    stat_branch_dir = Path(dir_gene_family) / "stat_branch"
    if not stat_branch_dir.is_dir():
        raise FileNotFoundError(f"{mode} stat_branch directory was not found: {stat_branch_dir}")

    query_status = {}
    count_columns = {}
    unresolved_rows = []
    emitted_family_order = []
    for family_id in family_ids:
        stat_branch = stat_branch_dir / f"{family_id}{STAT_BRANCH_SUFFIX}"
        raw_counts = read_stat_branch_copy_numbers(stat_branch)
        if stat_branch.is_file() and stat_branch.stat().st_size > 0:
            query_status[family_id] = "complete"
        else:
            query_status[family_id] = "missing_stat_branch"
        if query_status[family_id] != "complete" and not include_incomplete:
            continue
        emitted_family_order.append(family_id)
        if query_status[family_id] == "complete":
            resolved_counts, unresolved = resolve_counts_to_species_order(raw_counts, species_order)
            count_columns[family_id] = [resolved_counts[species] for species in species_order]
            for species, count in unresolved.items():
                unresolved_rows.append({"query": family_id, "family_id": family_id, "species": species, "copy_number": count})
        else:
            count_columns[family_id] = [numpy.nan for _ in species_order]

    copy_number = pandas.DataFrame(count_columns, index=species_order)
    copy_number.index.name = "species"
    presence = copy_number.copy()
    for col in presence.columns:
        presence[col] = pandas.to_numeric(presence[col], errors="coerce")
    presence = (presence > 0).astype("Int64").where(~copy_number.isna(), pandas.NA)

    long_df = build_long_table(
        copy_number=copy_number,
        presence=presence,
        species_order=species_order,
        family_order=emitted_family_order,
        query_status=query_status,
        mode=mode,
    )
    unresolved_df = pandas.DataFrame(unresolved_rows)
    return presence, copy_number, long_df, unresolved_df, query_status, family_source


def write_table(path, df):
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False, na_rep="NA")


def run(args):
    include_incomplete = bool(int(args.include_incomplete))
    presence, copy_number, long_df, unresolved_df, query_status, family_source = build_matrices(
        mode=args.mode,
        dir_gene_family=args.dir_gene_family,
        dir_query_gene=args.dir_query_gene,
        orthogroup_genecount=args.orthogroup_genecount,
        species_tree=args.species_tree,
        include_incomplete=include_incomplete,
    )
    if presence.shape[1] == 0:
        raise ValueError(
            f"No completed {args.mode} stat_branch files were detected. "
            "Use --include_incomplete 1 to emit placeholder columns."
        )

    plot_ids, selection_df = select_plot_families(
        available_ids=list(copy_number.columns),
        query_status=query_status,
        mode=args.mode,
        max_families=args.max_families,
        family_ids=args.family_ids,
        family_file=args.family_file,
    )
    if len(plot_ids) == 0:
        raise ValueError("No gene family IDs were selected for plotting.")

    plot_presence = presence.loc[:, plot_ids]
    plot_copy_number = copy_number.loc[:, plot_ids]
    plot_long = long_df[long_df["query"].isin(plot_ids)].copy()
    plot_order = {family_id: idx + 1 for idx, family_id in enumerate(plot_ids)}
    plot_long["query_order"] = plot_long["query"].map(plot_order)
    plot_long["family_order"] = plot_long["query_order"]
    plot_long = plot_long.sort_values(["query_order", "species_order"]).reset_index(drop=True)

    write_matrix(args.out_presence, presence)
    write_matrix(args.out_copy_number, copy_number)
    write_table(args.out_long, long_df)
    write_matrix(args.out_plot_presence, plot_presence)
    write_matrix(args.out_plot_copy_number, plot_copy_number)
    write_table(args.out_plot_long, plot_long)
    write_table(args.out_selection, selection_df)

    num_complete = sum(1 for status in query_status.values() if status == "complete")
    print(f"{args.mode} presence/absence matrix: {num_complete} completed gene family file(s).")
    print(f"Gene family ID source: {family_source}")
    print(f"Presence matrix: {args.out_presence}")
    print(f"Copy-number matrix: {args.out_copy_number}")
    print(f"Long table: {args.out_long}")
    print(f"Plot selection: {args.out_selection}")
    print(f"Plot families: {len(plot_ids)} / {presence.shape[1]}")
    print(f"Plot long table: {args.out_plot_long}")
    if unresolved_df.shape[0] > 0:
        unresolved_path = Path(args.out_long).with_suffix(".unresolved_species.tsv")
        unresolved_df.to_csv(unresolved_path, sep="\t", index=False)
        print(f"Warning: {unresolved_df.shape[0]} unresolved species/count row(s) were written to: {unresolved_path}")


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
