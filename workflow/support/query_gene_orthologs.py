#!/usr/bin/env python3
"""Summarize anchor genes and their ortholog copies.

A gene-tree tip belongs to a reference-gene column when its MRCA with that
reference-species tip is a speciation node (or when it is the reference tip
itself). A tip assigned to more than one adjacent column is represented as a
shared ancestral copy predating the duplication that separated those
reference genes. Tips are grouped by species and reference-gene set so their
group size is the plotted copy number. Plot labels are verified against the
family CDS FASTA. The tree table retains every original duplication and records
which displayed genes descend from each child, so bars count only D nodes with
displayed genes in both children without changing the original S/D calls.
A fourth table compares
each plotted copy with each covered reference gene using the family's local
synteny neighborhoods. Two or more distinct shared neighbor-similarity groups
provide local-synteny support; a single shared group is retained as
single-anchor evidence but is not called supported.
A fifth table retains candidate/reference pair provenance for Gene tree UFBoot,
while requiring every pair represented by one glyph to resolve to the same
orthology-defining speciation branch and support value.

The default reference-species basis preserves the historical output fields,
with additive displayed-descendant provenance in the tree table.
The query-gene basis coalesces query records that select the same gene-tree
tip, retains every original record in a query-to-anchor mapping table, and
writes semantically explicit ``anchor_*`` fields.
"""

import argparse
import csv
import gzip
import io
import math
import re
import warnings
from collections import defaultdict
from pathlib import Path

if __package__:
    from .gene_family_output_store import GeneFamilyOutputStore
else:
    from gene_family_output_store import GeneFamilyOutputStore

STAT_BRANCH_SUFFIX = "_stat.branch.tsv"
SYNTENY_SUFFIX = "_synteny.tsv"
SYNTENY_SUPPORT_MIN_ANCHORS = 2
DUP_CONF_FIELDS = [
    "family_id", "family_order", "species", "reference_cds_fasta_id",
    "candidate_cds_fasta_id", "mrca_branch_id", "mrca_event",
    "shared_species_count", "union_species_count", "dup_conf_score",
    "dup_conf_score_threshold", "branch_ufboot", "branch_ufboot_source",
]
QUERY_DUP_CONF_FIELDS = [field.replace("reference_", "anchor_", 1) for field in DUP_CONF_FIELDS]
COLUMN_FIELDS = [
    "column_order",
    "family_id",
    "family_order",
    "reference_species",
    "cds_fasta_id",
    "gene_id",
    "plot_label",
    "reference_tip_branch_id",
]
GLYPH_FIELDS = [
    "species",
    "family_id",
    "family_order",
    "reference_species",
    "relation",
    "reference_cds_fasta_ids",
    "reference_gene_ids",
    "reference_gene_count",
    "copy_number",
    "gene_ids",
    "start_order",
    "end_order",
    "is_contiguous",
    "lane_index",
    "lane_count",
]
TREE_FIELDS = [
    "family_id",
    "family_order",
    "reference_species",
    "node_id",
    "parent_node_id",
    "is_tip",
    "event",
    "cds_fasta_id",
    "gene_id",
    "column_order",
    "node_height",
    "plot_order",
    "mapped_species_node",
    "duplication_index",
    "in_reference_tree",
    "displayed_gene_ids",
    "displayed_child1_gene_ids",
    "displayed_child2_gene_ids",
]
SYNTENY_FIELDS = [
    "family_id",
    "family_order",
    "species",
    "reference_species",
    "relation",
    "reference_cds_fasta_id",
    "reference_gene_id",
    "column_order",
    "candidate_cds_fasta_id",
    "glyph_copy_number",
    "glyph_start_order",
    "glyph_end_order",
    "glyph_lane_index",
    "glyph_lane_count",
    "synteny_status",
    "support_min_anchor_count",
    "synteny_window_radius",
    "reference_neighbor_count",
    "candidate_neighbor_count",
    "flank_coverage",
    "shared_anchor_count",
    "local_synteny_score",
    "collinear_anchor_count",
    "collinearity_ratio",
    "collinear_orientation",
    "shared_group_ids",
]
UFBOOT_FIELDS = [
    "family_id",
    "family_order",
    "species",
    "reference_species",
    "relation",
    "reference_cds_fasta_id",
    "reference_gene_id",
    "column_order",
    "candidate_cds_fasta_id",
    "glyph_copy_number",
    "glyph_start_order",
    "glyph_end_order",
    "glyph_lane_index",
    "glyph_lane_count",
    "orthology_mrca_branch_id",
    "orthology_mrca_event",
    "ufboot_support_source",
    "decisive_branch_ufboot",
    "orthology_ufboot_status",
    "orthology_ufboot_unavailable_reason",
]
QUERY_COLUMN_FIELDS = [
    "column_order",
    "family_id",
    "family_order",
    "basis",
    "query_ids",
    "query_labels",
    "query_count",
    "anchor_source",
    "anchor_species",
    "anchor_cds_fasta_id",
    "anchor_gene_id",
    "plot_label",
    "anchor_tip_branch_id",
]
QUERY_GLYPH_FIELDS = [
    "species",
    "family_id",
    "family_order",
    "basis",
    "relation",
    "anchor_cds_fasta_ids",
    "anchor_gene_ids",
    "anchor_query_ids",
    "anchor_count",
    "copy_number",
    "gene_ids",
    "start_order",
    "end_order",
    "is_contiguous",
    "lane_index",
    "lane_count",
]
QUERY_TREE_FIELDS = [
    "family_id",
    "family_order",
    "basis",
    "node_id",
    "parent_node_id",
    "is_tip",
    "event",
    "anchor_cds_fasta_id",
    "anchor_gene_id",
    "column_order",
    "node_height",
    "plot_order",
    "mapped_species_node",
    "duplication_index",
    "in_anchor_tree",
    "displayed_gene_ids",
    "displayed_child1_gene_ids",
    "displayed_child2_gene_ids",
]
QUERY_SYNTENY_FIELDS = [
    field.replace("reference_", "anchor_") if field.startswith("reference_") else field
    for field in SYNTENY_FIELDS
]
QUERY_SYNTENY_FIELDS[QUERY_SYNTENY_FIELDS.index("anchor_species")] = "basis"
QUERY_UFBOOT_FIELDS = [
    field.replace("reference_", "anchor_") if field.startswith("reference_") else field
    for field in UFBOOT_FIELDS
]
QUERY_UFBOOT_FIELDS[QUERY_UFBOOT_FIELDS.index("anchor_species")] = "basis"
QUERY_MAP_FIELDS = [
    "family_id",
    "family_order",
    "query_order",
    "query_id",
    "query_label",
    "marker_source",
    "anchor_species",
    "anchor_cds_fasta_id",
    "anchor_gene_id",
    "anchor_tip_branch_id",
    "column_order",
    "merged_query_count",
    "source_species",
    "hog_ids",
]
SELECTION_FIELDS = [
    "family_id", "query_id", "source_species", "tree_species", "distance",
    "anchor_cds_fasta_id", "decision", "reason", "replaced_by",
]


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir_gene_family", metavar="PATH", required=True)
    parser.add_argument("--dir_query_gene", metavar="PATH", default="")
    parser.add_argument("--family_file", metavar="PATH", default="")
    parser.add_argument("--family_manifest", metavar="TSV", default="",
                        help="Ordered saved-family sources; see docs/presence-absence.md")
    parser.add_argument("--query_metadata", metavar="TSV", default="")
    parser.add_argument("--query_selection", choices=("all", "closest"), default="all")
    parser.add_argument("--species_tree", metavar="NEWICK", default="")
    parser.add_argument("--target_species", default="")
    parser.add_argument("--selection_species", default="",
                        help="Comma-separated species whose exact ortholog union must be preserved; empty=all")
    parser.add_argument("--query_label", choices=("id", "label"), default="id")
    parser.add_argument("--out_selection", metavar="TSV", default="")
    parser.add_argument("--out_overlap", metavar="TSV", default="")
    parser.add_argument("--out_long", metavar="TSV", default="")
    parser.add_argument(
        "--basis",
        choices=("reference_species", "query_gene"),
        default="reference_species",
    )
    parser.add_argument("--reference_species", metavar="SPECIES", default="")
    parser.add_argument("--out_columns", metavar="PATH", required=True)
    parser.add_argument("--out_glyphs", metavar="PATH", required=True)
    parser.add_argument("--out_tree", metavar="PATH", required=True)
    parser.add_argument("--out_synteny", metavar="PATH", required=True)
    parser.add_argument("--out_ufboot", metavar="PATH", required=True)
    parser.add_argument("--out_query_map", metavar="PATH", default="")
    parser.add_argument("--dup_conf_score_threshold", "--dup-conf-score-threshold",
                        type=validate_dup_conf_threshold, default=0,
                        help="0 disables extra candidates; positive values flag cross-species D-MRCA pairs with Jaccard score <= threshold")
    parser.add_argument("--out_dup_conf", metavar="TSV", default="",
                        help="Pairwise provenance for weak-duplication candidates")
    return parser


def validate_dup_conf_threshold(value):
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("dup_conf_score_threshold must be a finite number between 0 and 1") from exc
    if not math.isfinite(number) or not 0 <= number <= 1:
        raise ValueError("dup_conf_score_threshold must be a finite number between 0 and 1")
    return number


def species_overlap_counts(by_id, children, root):
    """Recompute Jaccard counts on the full saved tree without changing events."""
    species_sets, counts = {}, {}
    stack = [(root, False)]
    while stack:
        node, visited = stack.pop()
        node_children = children.get(node, [])
        if not node_children:
            species = normalize_species_label(by_id[node].get("spnode_coverage", ""))
            if not species:
                raise ValueError(f"Weak-duplication candidates require a species for every tree tip: branch={node}")
            species_sets[node] = {species}
            continue
        if not visited:
            stack.append((node, True))
            stack.extend((child, False) for child in reversed(node_children))
            continue
        if len(node_children) != 2:
            raise ValueError("Weak-duplication candidates require a strictly binary saved gene tree")
        left, right = (species_sets[child] for child in node_children)
        species_sets[node] = left | right
        counts[node] = (len(left & right), len(species_sets[node]))
        expected_event = "D" if counts[node][0] else "S"
        if by_id[node].get("so_event") != expected_event:
            raise ValueError(
                f"Saved species-overlap event disagrees with descendant species: "
                f"branch={node}, saved={by_id[node].get('so_event')!r}, expected={expected_event}"
            )
    return counts


def add_weak_duplication_candidates(store, columns, glyphs, threshold):
    """Add separate cross-species candidate glyphs; retain all strict calls and D nodes."""
    threshold = validate_dup_conf_threshold(threshold)
    if threshold == 0:
        return glyphs, []
    columns_by_family = defaultdict(list)
    for column in columns:
        columns_by_family[str(column["family_id"])].append(column)
    additions, evidence = [], []
    for family_id, family_columns in columns_by_family.items():
        rows = read_stat_branch(store, family_id)
        if not rows:
            raise ValueError(f"Weak-duplication candidates require stat_branch rows: family={family_id}")
        by_id, children, root = build_tree_index(rows)
        counts = species_overlap_counts(by_id, children, root)
        cds_ids = set(read_family_cds_fasta_ids(store, family_id))
        ancestors = {node: ancestor_chain(by_id, node) for node in by_id}
        support, support_source = normalized_ufboot_by_branch(by_id, family_id)
        family_columns.sort(key=lambda column: int(column["column_order"]))
        grouped = defaultdict(list)
        for tip, row in by_id.items():
            if row.get("so_event") != "L":
                continue
            species = str(row.get("spnode_coverage") or "").strip()
            candidate_id = str(row.get("node_name") or "").strip()
            matches_by_mrca = defaultdict(list)
            for column in family_columns:
                anchor = int(column["reference_tip_branch_id"])
                if normalize_species_label(species) == normalize_species_label(by_id[anchor].get("spnode_coverage", "")):
                    continue
                mrca = mrca_node(ancestors, tip, anchor)
                if by_id[mrca].get("so_event") != "D":
                    continue
                shared, union = counts[mrca]
                score = shared / union
                if score > threshold:
                    continue
                if not candidate_id or candidate_id not in cds_ids:
                    raise ValueError(
                        f"Additional ortholog candidate is absent from CDS FASTA: "
                        f"family={family_id}, candidate={candidate_id!r}"
                    )
                matches_by_mrca[mrca].append(column)
                evidence.append(dict(
                    family_id=family_id, family_order=column["family_order"], species=species,
                    reference_cds_fasta_id=column["cds_fasta_id"], candidate_cds_fasta_id=candidate_id,
                    mrca_branch_id=mrca, mrca_event="D", shared_species_count=shared,
                    union_species_count=union, dup_conf_score=score, dup_conf_score_threshold=threshold,
                    branch_ufboot=support.get(mrca) if mrca != root and support.get(mrca) is not None else "",
                    branch_ufboot_source=support_source,
                ))
            for mrca, matches in matches_by_mrca.items():
                # Split gapped anchor sets so no unassigned column is painted.
                runs = []
                for column in matches:
                    if not runs or int(column["column_order"]) != int(runs[-1][-1]["column_order"]) + 1:
                        runs.append([])
                    runs[-1].append(column)
                for run in runs:
                    grouped[(species, mrca, tuple(int(column["column_order"]) for column in run))].append(candidate_id)
        column_by_order = {int(column["column_order"]): column for column in family_columns}
        for (species, _mrca, orders), genes in sorted(grouped.items()):
            covered = [column_by_order[order] for order in orders]
            glyph = dict(
                species=species, family_id=family_id, family_order=covered[0]["family_order"],
                reference_species=covered[0]["reference_species"], relation="weak_duplication",
                reference_cds_fasta_ids=";".join(column["cds_fasta_id"] for column in covered),
                reference_gene_ids=";".join(column["gene_id"] for column in covered),
                reference_gene_count=len(covered), copy_number=len(genes), gene_ids=";".join(sorted(genes)),
                start_order=min(orders), end_order=max(orders), is_contiguous=1, lane_index=1, lane_count=1,
            )
            if "query_ids" in covered[0]:
                glyph["anchor_query_ids"] = ";".join(column["query_ids"] for column in covered)
            additions.append(glyph)
    result = [dict(glyph) for glyph in glyphs] + additions
    assign_lanes(result)
    result.sort(key=lambda row: (str(row["species"]), int(row["start_order"]), int(row["end_order"]), str(row["relation"])))
    evidence.sort(key=lambda row: (int(row["family_order"]), str(row["species"]), str(row["candidate_cds_fasta_id"]), str(row["reference_cds_fasta_id"])))
    return result, evidence


def _open_query_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _query_label_from_header(header, query_id):
    pipe_fields = [field.strip() for field in header.split("|")]
    if len(pipe_fields) > 1 and pipe_fields[1]:
        return pipe_fields[1]
    return query_id


def read_query_definitions(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with _open_query_text(path) as handle:
        first = handle.read(1)
        handle.seek(0)
        definitions = []
        if first == ">":
            for line in handle:
                if not line.startswith(">"):
                    continue
                header = line[1:].strip()
                query_id = header.split()[0] if header else ""
                if query_id:
                    definitions.append(
                        {
                            "query_id": query_id,
                            "query_label": _query_label_from_header(header, query_id),
                            "source_species": _source_species_from_header(header),
                        }
                    )
        else:
            for line in handle:
                query_id = line.strip()
                if query_id:
                    definitions.append({"query_id": query_id, "query_label": query_id})
    seen = set()
    unique = []
    for definition in definitions:
        query_id = definition["query_id"]
        if query_id in seen:
            continue
        seen.add(query_id)
        unique.append(definition)
    return unique


def _source_species_from_header(header):
    # Explicit metadata only: never infer source from the best-hit tip or ID prefixes.
    match = re.search(r"(?:^|[|\s])species=(.*?)(?=\s+[A-Za-z_]+=[^=]|\||$)", header)
    return normalize_species_label(match.group(1)) if match else ""


def read_query_metadata(path):
    if not path:
        return {}
    metadata = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"family_id", "query_id", "source_species"}.issubset(reader.fieldnames or []):
            raise ValueError("Query metadata requires family_id, query_id, source_species")
        for row in reader:
            key = (row["family_id"].strip(), row["query_id"].strip())
            if not all(key) or key in metadata:
                raise ValueError(f"Empty or duplicate query metadata key: {key}")
            metadata[key] = {field: normalize_species_label(row.get(field))
                             for field in ("source_species", "tree_species")}
    return metadata


def select_closest_queries(rows, definitions, family_id, species_tree, target_species,
                           selection_species="", policy="all"):
    """Conservatively replace only farther, orthologous, coverage-redundant queries."""
    assignments = select_query_tip_assignments(rows, definitions)
    by_id, _children, _root = build_tree_index(rows)
    ancestors = {node: ancestor_chain(by_id, node) for node in by_id}
    distances = {}
    if policy == "closest":
        if not species_tree or not target_species:
            raise ValueError("closest requires --species_tree and --target_species")
        from Bio import Phylo

        tree = Phylo.read(species_tree, "newick")
        tips = {normalize_species_label(tip.name): tip for tip in tree.get_terminals()}
        if len(tips) != len(tree.get_terminals()) or "" in tips:
            raise ValueError("Species tree must have unique non-empty tip labels")
        target = normalize_species_label(target_species)
        if target not in tips:
            raise ValueError(f"Target species is absent from species tree: {target}")
        # Rank shared ancestry, not total path length: a densely sampled source
        # clade must not make a genuinely nearer relative appear farther away.
        # Count target-lineage edges back to the source/target MRCA, independent
        # of missing/non-comparable branch lengths.
        paths = {name: [tree.root] + tree.get_path(tip) for name, tip in tips.items()}
        for name, path in paths.items():
            common = sum(a is b for a, b in zip(path, paths[target], strict=False))
            distances[name] = len(paths[target]) - common
    selected_species = {normalize_species_label(value) for value in selection_species.split(",")
                        if value.strip()}
    if selected_species:
        if not species_tree:
            raise ValueError("--selection_species requires --species_tree")
        from Bio import Phylo

        available_species = {normalize_species_label(tip.name)
                             for tip in Phylo.read(species_tree, "newick").get_terminals()}
        if selected_species - available_species:
            raise ValueError(f"Selection species absent from species tree: {sorted(selected_species - available_species)}")
        # A tree-tip species with zero family members is valid biological non-detection.
    coverage = {}
    for definition in definitions if policy == "closest" else []:
        query_id = definition["query_id"]
        tip = assignments[query_id]["tip"]
        coverage[query_id] = {
            str(row["node_name"]) for node, row in by_id.items()
            if row.get("so_event") == "L"
            and (not selected_species or normalize_species_label(row.get("spnode_coverage"))
                 in selected_species)
            and (node == tip or by_id[mrca_node(ancestors, node, tip)].get("so_event") == "S")
        }
    distance_by_query = {
        definition["query_id"]: distances.get(
            definition.get("tree_species") or definition.get("source_species", ""))
        for definition in definitions
    }
    kept = []
    decisions = {}
    # Nearer candidates are finalized first. Equal distances never eliminate each other.
    ordered = sorted(definitions, key=lambda d: (
        distance_by_query[d["query_id"]] is None,
        distance_by_query[d["query_id"]] or 0,
    ))
    for definition in ordered:
        query_id = definition["query_id"]
        distance = distance_by_query[query_id]
        replacement = ""
        if policy == "closest" and distance is not None:
            tip = assignments[query_id]["tip"]
            for candidate in kept:
                candidate_id = candidate["query_id"]
                candidate_distance = distance_by_query[candidate_id]
                if candidate_distance is None or candidate_distance >= distance:
                    continue
                other_tip = assignments[candidate_id]["tip"]
                same_lineage = tip == other_tip or by_id[
                    mrca_node(ancestors, tip, other_tip)].get("so_event") == "S"
                if same_lineage and coverage[query_id].issubset(coverage[candidate_id]):
                    replacement = candidate_id
                    break
        if not replacement:
            kept.append(definition)
        decisions[query_id] = (replacement, distance)
    retained_ids = {definition["query_id"] for definition in kept}
    before = set().union(*coverage.values()) if coverage else set()
    after = set().union(*(coverage[q] for q in retained_ids)) if retained_ids and coverage else set()
    if before != after:
        raise ValueError(f"Query selection changed the target ortholog union: family={family_id}")
    audit = []
    for definition in definitions:
        query_id = definition["query_id"]
        replacement, distance = decisions[query_id]
        audit.append(dict(
            family_id=family_id, query_id=query_id,
            source_species=definition.get("source_species", ""),
            tree_species=definition.get("tree_species") or definition.get("source_species", ""),
            distance="NA" if distance is None else distance,
            anchor_cds_fasta_id=by_id[assignments[query_id]["tip"]]["node_name"],
            decision="excluded" if replacement else "retained",
            reason=("nearer_ortholog_covers_target_union" if replacement else
                    "all_queries" if policy == "all" else
                    "unknown_source_or_tree_species" if distance is None else
                    "no_strictly_nearer_redundant_ortholog"),
            replaced_by=replacement,
        ))
    return [definition for definition in definitions if definition["query_id"] in retained_ids], audit


def read_family_manifest(path):
    """Read explicit saved tree sources without copying or pruning their artifacts."""
    base = Path(path).resolve().parent
    records = []
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"family_id", "source_dir", "source_family_id"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Family manifest requires {sorted(required)}")
        for row in reader:
            row = {key: (value or "").strip() for key, value in row.items()}
            if not all(row.get(key) for key in required):
                raise ValueError("Family manifest has empty required values")
            for key in ("source_dir", "query_file", "hog_table"):
                if row.get(key):
                    row[key] = str((base / row[key]).resolve())
            if bool(row.get("query_file")) == bool(row.get("anchor_species")):
                raise ValueError("Each manifest row requires query_file OR anchor_species")
            if bool(row.get("hog_ids")) != bool(row.get("hog_table")):
                raise ValueError("hog_ids and hog_table must be supplied together")
            records.append(row)
    ids = [row["family_id"] for row in records]
    if not ids or len(ids) != len(set(ids)):
        raise ValueError("Family manifest requires unique, non-empty family_id values")
    return records


class ManifestOutputStore:
    """Route logical plot IDs to their original live-or-ZIP family artifacts."""

    def __init__(self, records):
        self.records = {row["family_id"]: row for row in records}
        self.stores = {row["source_dir"]: GeneFamilyOutputStore(row["source_dir"])
                       for row in records}

    def _resolve(self, name):
        for family_id in sorted(self.records, key=len, reverse=True):
            if name.startswith(family_id + "_"):
                row = self.records[family_id]
                return self.stores[row["source_dir"]], row["source_family_id"] + name[len(family_id):]
        raise ValueError(f"No manifest source for artifact: {name}")

    def artifact(self, subdir, name):
        store, original = self._resolve(name)
        return store.artifact(subdir, original)

    def open_binary(self, subdir, name):
        store, original = self._resolve(name)
        return store.open_binary(subdir, original)


def manifest_query_definitions(record, rows, cds_fasta_ids):
    definitions = read_query_definitions(record["query_file"]) if record.get("query_file") else []
    members = defaultdict(set)
    if record.get("hog_table"):
        selected = set(record["hog_ids"].split(";"))
        observed = set()
        with _open_query_text(record["hog_table"]) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"HOG", "OG"}.issubset(reader.fieldnames or []):
                raise ValueError("HOG table requires HOG and OG columns")
            tip_species = {str(row["node_name"]): normalize_species_label(row.get("spnode_coverage"))
                           for row in rows if row.get("so_event") == "L"}
            tip_ids = set(tip_species)
            for hog in reader:
                if hog["HOG"] not in selected:
                    continue
                if hog["HOG"] in observed or hog["OG"] != record["source_family_id"]:
                    raise ValueError("Selected HOG is duplicated or belongs to another source OG")
                observed.add(hog["HOG"])
                for species, text in hog.items():
                    if species in ("HOG", "OG", "Gene Tree Parent Clade"):
                        continue
                    for gene in re.split(r"[,;]\s*", text or ""):
                        if not gene.strip():
                            continue
                        gene = gene.strip()
                        if gene not in tip_ids:
                            matches = [str(row["node_name"]) for row in rows
                                       if row.get("so_event") == "L"
                                       and normalize_species_label(row.get("spnode_coverage"))
                                       == normalize_species_label(species)
                                       and gene_id_from_cds_fasta_id(row["node_name"], species) == gene]
                            if len(matches) != 1:
                                raise ValueError(f"HOG member not uniquely present in saved tree: {gene}")
                            gene = matches[0]
                        if gene not in cds_fasta_ids:
                            raise ValueError(f"HOG member missing from saved CDS FASTA: {gene}")
                        if tip_species[gene] != normalize_species_label(species):
                            raise ValueError(f"HOG member is in the wrong species column: {gene}")
                        members[gene].add(hog["HOG"])
        if selected != observed:
            raise ValueError(f"Selected HOGs missing from table: {sorted(selected - observed)}")
    if record.get("anchor_species"):
        species = normalize_species_label(record["anchor_species"])
        for row in rows:
            if row.get("so_event") != "L" or normalize_species_label(row.get("spnode_coverage")) != species:
                continue
            gene = str(row["node_name"])
            if record.get("hog_table") and gene not in members:
                continue
            definitions.append(dict(query_id=gene,
                                    query_label=";".join(sorted(members[gene])) + " " + gene
                                    if members[gene] else gene,
                                    source_species=species))
            # Mark an exact saved tip; do not search or rebuild the full OG tree.
            row["query_marker_source"] = str(row.get("query_marker_source") or "") + "|direct:" + gene
        if not definitions:
            raise ValueError(f"No selected anchor-species genes in saved tree: {record['family_id']}")
    for definition in definitions:
        if members:
            # HOG membership remains annotation, never a replacement tree or orthology definition.
            definition["hog_ids"] = ";".join(sorted(members.get(definition["query_id"], set())))
    return definitions


def read_family_ids(query_dir, family_file=""):
    query_dir = Path(query_dir)
    available = sorted(
        path.name for path in query_dir.iterdir() if path.is_file() and not path.name.startswith(".")
    )
    if not family_file:
        return available
    selected = []
    with Path(family_file).open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames and "family_id" in reader.fieldnames:
            selected = [str(row.get("family_id") or "").strip() for row in reader]
        else:
            handle.seek(0)
            selected = [line.strip().split("\t", 1)[0] for line in handle if line.strip()]
    available_set = set(available)
    seen = set()
    selected_available = []
    for family_id in selected:
        if family_id not in available_set or family_id in seen:
            continue
        seen.add(family_id)
        selected_available.append(family_id)
    return selected_available


def parse_query_marker_sources(value):
    sources = []
    for group in str(value or "").split("|"):
        group = group.strip()
        if not group or ":" not in group:
            continue
        source_type, query_text = group.split(":", 1)
        source_type = source_type.strip().lower()
        for query_id in query_text.split(";"):
            query_id = query_id.strip()
            if query_id:
                sources.append((source_type, query_id))
    return sources


def _int_value(value, default=-999):
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return default


def read_stat_branch(store, family_id):
    name = f"{family_id}{STAT_BRANCH_SUFFIX}"
    artifact = store.artifact("stat_branch", name)
    if artifact is None or not artifact.size:
        return []
    with store.open_binary("stat_branch", name) as binary_handle:
        with io.TextIOWrapper(binary_handle, encoding="utf-8", errors="replace", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))


def read_family_synteny(store, family_id):
    """Read one family's normalized local-neighborhood table.

    An absent or header-only table means that local synteny is unavailable for
    this family. A present malformed table is a hard error because silently
    treating corrupt evidence as biological absence would be misleading.
    """

    name = f"{family_id}{SYNTENY_SUFFIX}"
    artifact = store.artifact("synteny", name)
    if artifact is None or not artifact.size:
        return []
    with store.open_binary("synteny", name) as binary_handle:
        with io.TextIOWrapper(binary_handle, encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"node_name", "offset", "group_id"}
            observed = set(reader.fieldnames or [])
            if not required.issubset(observed):
                raise ValueError(
                    "Synteny table is missing required columns: "
                    f"family={family_id}, missing={sorted(required - observed)}"
                )
            rows = []
            seen = set()
            for row_number, row in enumerate(reader, start=2):
                node_name = str(row.get("node_name") or "").strip()
                group_id = str(row.get("group_id") or "").strip()
                offset_text = str(row.get("offset") or "").strip()
                if not node_name or not group_id or not offset_text:
                    continue
                try:
                    offset_number = float(offset_text)
                except ValueError as exc:
                    raise ValueError(
                        "Synteny table contains a non-numeric offset: "
                        f"family={family_id}, row={row_number}, value={offset_text!r}"
                    ) from exc
                if not offset_number.is_integer() or int(offset_number) == 0:
                    raise ValueError(
                        "Synteny offsets must be non-zero integers: "
                        f"family={family_id}, row={row_number}, value={offset_text!r}"
                    )
                offset = int(offset_number)
                neighbor_gene = str(row.get("neighbor_gene") or "").strip()
                key = (node_name, offset, neighbor_gene, group_id)
                if key in seen:
                    continue
                seen.add(key)
                rows.append(
                    {
                        "node_name": node_name,
                        "offset": offset,
                        "neighbor_gene": neighbor_gene,
                        "group_id": group_id,
                    }
                )
    return rows


def build_synteny_neighborhoods(rows):
    neighborhoods = defaultdict(list)
    window_radius = 0
    offsets_by_node = defaultdict(dict)
    for row in rows:
        node_name = row["node_name"]
        offset = int(row["offset"])
        group_id = row["group_id"]
        previous_group = offsets_by_node[node_name].get(offset)
        if previous_group is not None and previous_group != group_id:
            raise ValueError(
                "Synteny table maps one focal-gene offset to multiple similarity groups: "
                f"node={node_name}, offset={offset}, groups={[previous_group, group_id]}"
            )
        offsets_by_node[node_name][offset] = group_id
        neighborhoods[node_name].append(row)
        window_radius = max(window_radius, abs(offset))
    for node_name in neighborhoods:
        neighborhoods[node_name].sort(
            key=lambda row: (int(row["offset"]), str(row["group_id"]), str(row["neighbor_gene"]))
        )
    return neighborhoods, window_radius


def _representative_group_offsets(rows):
    offsets_by_group = defaultdict(list)
    for row in rows:
        offsets_by_group[str(row["group_id"])].append(int(row["offset"]))
    return {
        group_id: min(offsets, key=lambda offset: (abs(offset), offset))
        for group_id, offsets in offsets_by_group.items()
    }


def _longest_strict_monotonic_subsequence(values, increasing=True):
    if not values:
        return 0
    lengths = [1] * len(values)
    for right in range(len(values)):
        for left in range(right):
            ordered = values[left] < values[right] if increasing else values[left] > values[right]
            if ordered:
                lengths[right] = max(lengths[right], lengths[left] + 1)
    return max(lengths)


def local_synteny_metrics(reference_rows, candidate_rows, window_radius):
    """Return transparent anchor-count and order summaries for one gene pair.

    Each neighbor-similarity group contributes at most one independent anchor,
    so tandem/proximal copies from one broad group cannot inflate support.
    """

    reference_offsets = _representative_group_offsets(reference_rows)
    candidate_offsets = _representative_group_offsets(candidate_rows)
    shared_groups = sorted(set(reference_offsets).intersection(candidate_offsets))
    shared_anchor_count = len(shared_groups)
    ordered_groups = sorted(
        shared_groups,
        key=lambda group_id: (reference_offsets[group_id], group_id),
    )
    candidate_order = [candidate_offsets[group_id] for group_id in ordered_groups]
    forward_count = _longest_strict_monotonic_subsequence(candidate_order, increasing=True)
    reverse_count = _longest_strict_monotonic_subsequence(candidate_order, increasing=False)
    collinear_anchor_count = max(forward_count, reverse_count)
    if collinear_anchor_count == 0:
        orientation = "none"
    elif forward_count > reverse_count:
        orientation = "forward"
    elif reverse_count > forward_count:
        orientation = "reverse"
    else:
        orientation = "ambiguous"

    expected_neighbor_count = max(0, 2 * int(window_radius))
    reference_neighbor_count = len({int(row["offset"]) for row in reference_rows})
    candidate_neighbor_count = len({int(row["offset"]) for row in candidate_rows})
    if expected_neighbor_count > 0:
        local_synteny_score = shared_anchor_count / expected_neighbor_count
        flank_coverage = min(reference_neighbor_count, candidate_neighbor_count) / expected_neighbor_count
    else:
        local_synteny_score = ""
        flank_coverage = ""
    collinearity_ratio = (
        collinear_anchor_count / shared_anchor_count if shared_anchor_count > 0 else 0.0
    )
    return {
        "reference_neighbor_count": reference_neighbor_count,
        "candidate_neighbor_count": candidate_neighbor_count,
        "flank_coverage": flank_coverage,
        "shared_anchor_count": shared_anchor_count,
        "local_synteny_score": local_synteny_score,
        "collinear_anchor_count": collinear_anchor_count,
        "collinearity_ratio": collinearity_ratio,
        "collinear_orientation": orientation,
        "shared_group_ids": ";".join(ordered_groups),
    }


def read_fasta_ids_from_store(store, subdir, name):
    with store.open_binary(subdir, name) as binary_handle:
        if name.endswith(".gz"):
            sequence_handle = gzip.GzipFile(fileobj=binary_handle, mode="rb")
        else:
            sequence_handle = binary_handle
        with io.TextIOWrapper(
            sequence_handle, encoding="utf-8", errors="replace", newline=""
        ) as text_handle:
            return [
                line[1:].strip().split()[0]
                for line in text_handle
                if line.startswith(">") and line[1:].strip()
            ]


def read_family_cds_fasta_ids(store, family_id):
    candidates = [
        f"{family_id}_cds.fa.gz",
        f"{family_id}_cds.fasta.gz",
        f"{family_id}_cds.fasta",
        f"{family_id}_cds.fa",
    ]
    for name in candidates:
        artifact = store.artifact("cds_fasta", name)
        if artifact is not None and artifact.size:
            ids = read_fasta_ids_from_store(store, "cds_fasta", name)
            if ids:
                seen_ids = set()
                duplicate_ids = set()
                for cds_fasta_id in ids:
                    if cds_fasta_id in seen_ids:
                        duplicate_ids.add(cds_fasta_id)
                    seen_ids.add(cds_fasta_id)
                if duplicate_ids:
                    raise ValueError(
                        "Family CDS FASTA contains duplicate sequence IDs: "
                        f"family={family_id}, ids={sorted(duplicate_ids)}"
                    )
                return ids
    raise FileNotFoundError(
        f"No non-empty CDS FASTA was found for reference-gene labels: family={family_id}, "
        f"checked={candidates}"
    )


def gene_id_from_cds_fasta_id(cds_fasta_id, species):
    cds_fasta_id = str(cds_fasta_id)
    species = str(species or "").strip()
    for separator in ("_", ".", "-"):
        prefix = f"{species}{separator}"
        if species and cds_fasta_id.startswith(prefix) and len(cds_fasta_id) > len(prefix):
            return cds_fasta_id[len(prefix) :]
    return cds_fasta_id


def resolve_query_cds_definitions(
    by_id,
    query_tip_by_id,
    definitions,
    cds_fasta_ids,
    family_id,
    preserve_query_label=False,
):
    definition_by_id = {definition["query_id"]: dict(definition) for definition in definitions}
    cds_fasta_id_set = set(cds_fasta_ids)
    query_ids_by_tip = defaultdict(list)
    for query_id, tip in query_tip_by_id.items():
        query_ids_by_tip[tip].append(query_id)
    duplicated_tips = {
        tip: query_ids for tip, query_ids in query_ids_by_tip.items() if len(query_ids) > 1
    }
    if duplicated_tips:
        raise ValueError(
            f"Multiple query records map to the same gene-tree tip: family={family_id}, "
            f"mappings={duplicated_tips}"
        )

    missing = []
    for query_id, tip in query_tip_by_id.items():
        cds_fasta_id = str(by_id[tip].get("node_name") or "").strip()
        if cds_fasta_id not in cds_fasta_id_set:
            missing.append((query_id, tip, cds_fasta_id))
            continue
        species = str(by_id[tip].get("spnode_coverage") or "").strip()
        definition = definition_by_id[query_id]
        definition["cds_fasta_id"] = cds_fasta_id
        definition["gene_id"] = gene_id_from_cds_fasta_id(cds_fasta_id, species)
        if not preserve_query_label:
            definition["query_label"] = definition["gene_id"]
    if missing:
        raise ValueError(
            f"Query gene-tree tips were not found in the family CDS FASTA: family={family_id}, "
            f"missing={missing}"
        )
    return definition_by_id


def select_query_tip_assignments(rows, definitions):
    definition_ids = {definition["query_id"] for definition in definitions}
    candidates = defaultdict(list)
    for row in rows:
        if str(row.get("so_event") or "") != "L":
            continue
        branch_id = _int_value(row.get("branch_id"))
        for source_type, query_id in parse_query_marker_sources(row.get("query_marker_source")):
            if query_id in definition_ids:
                priority = 0 if source_type == "direct" else 1
                candidates[query_id].append((priority, branch_id, source_type))

    selected = {}
    for definition in definitions:
        query_id = definition["query_id"]
        options = candidates.get(query_id, [])
        if not options:
            continue
        best_priority = min(priority for priority, _, _ in options)
        nodes = sorted(
            {node for priority, node, _ in options if priority == best_priority}
        )
        if len(nodes) != 1:
            raise ValueError(
                f"Query marker maps to multiple equally preferred gene-tree tips: query={query_id}, tips={nodes}"
            )
        selected_sources = sorted(
            {
                source_type
                for priority, node, source_type in options
                if priority == best_priority and node == nodes[0]
            }
        )
        selected[query_id] = {
            "tip": nodes[0],
            "marker_source": selected_sources[0],
        }
    missing_query_ids = [
        definition["query_id"]
        for definition in definitions
        if definition["query_id"] not in selected
    ]
    if missing_query_ids:
        raise ValueError(
            "Query records have no gene-tree tip marker: "
            f"queries={missing_query_ids}"
        )
    return selected


def select_query_tip_nodes(rows, definitions):
    """Return the historical query-to-tip mapping without marker provenance."""

    return {
        query_id: assignment["tip"]
        for query_id, assignment in select_query_tip_assignments(rows, definitions).items()
    }


def normalize_species_label(value):
    return "_".join(str(value or "").strip().replace("_", " ").split())


def select_reference_tip_nodes(rows, reference_species):
    reference_species = normalize_species_label(reference_species)
    selected = {}
    definitions = []
    for row in rows:
        if str(row.get("so_event") or "") != "L":
            continue
        if normalize_species_label(row.get("spnode_coverage")) != reference_species:
            continue
        branch_id = _int_value(row.get("branch_id"))
        cds_fasta_id = str(row.get("node_name") or "").strip()
        if branch_id < 0 or not cds_fasta_id:
            continue
        if cds_fasta_id in selected:
            raise ValueError(
                "Reference-species CDS FASTA ID maps to multiple gene-tree tips: "
                f"reference_species={reference_species}, cds_fasta_id={cds_fasta_id}"
            )
        selected[cds_fasta_id] = branch_id
        definitions.append({"query_id": cds_fasta_id, "query_label": cds_fasta_id})
    return selected, definitions


def build_tree_index(rows):
    if not rows:
        raise ValueError("stat_branch must contain at least one row")

    def required_integer(row, column, branch_context):
        value = row.get(column)
        try:
            text = str(value).strip()
            if not text or text.lower() == "nan":
                raise ValueError
            return int(text)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"stat_branch contains an invalid integer in {column}: "
                f"branch={branch_context}, value={value!r}"
            ) from exc

    by_id = {}
    parent_by_id = {}
    for row_number, row in enumerate(rows, start=2):
        branch_id = required_integer(row, "branch_id", f"row {row_number}")
        if branch_id < 0:
            raise ValueError(
                f"stat_branch branch_id must be non-negative: row={row_number}, branch_id={branch_id}"
            )
        if branch_id in by_id:
            raise ValueError(f"stat_branch contains duplicate branch_id: {branch_id}")
        by_id[branch_id] = row
        parent_by_id[branch_id] = required_integer(row, "parent", branch_id)

    roots = [branch_id for branch_id, parent in parent_by_id.items() if parent < 0]
    if len(roots) != 1:
        raise ValueError(f"stat_branch must contain exactly one root; observed roots={roots}")

    children_from_parent = defaultdict(list)
    for branch_id, parent in parent_by_id.items():
        if parent < 0:
            continue
        if parent == branch_id:
            raise ValueError(f"stat_branch branch cannot be its own parent: branch_id={branch_id}")
        if parent not in by_id:
            raise ValueError(
                f"stat_branch parent does not exist: branch_id={branch_id}, parent={parent}"
            )
        children_from_parent[parent].append(branch_id)

    child_columns = ("child1", "child2")
    has_any_child_column = any(
        any(column in row for column in child_columns) for row in rows
    )
    has_complete_child_schema = all(
        all(column in row for column in child_columns) for row in rows
    )
    if has_any_child_column and not has_complete_child_schema:
        raise ValueError(
            "stat_branch child columns are incomplete; provide both child1 and child2 "
            "for every row, or omit both columns"
        )

    children = defaultdict(list)
    if has_complete_child_schema:
        for branch_id, row in by_id.items():
            explicit = [
                required_integer(row, column, branch_id) for column in child_columns
            ]
            explicit = [child for child in explicit if child >= 0]
            if len(explicit) != len(set(explicit)):
                raise ValueError(
                    f"stat_branch lists the same child more than once: "
                    f"branch_id={branch_id}, children={explicit}"
                )
            missing_children = [child for child in explicit if child not in by_id]
            if missing_children:
                raise ValueError(
                    f"stat_branch child does not exist: branch_id={branch_id}, "
                    f"children={missing_children}"
                )
            implied = children_from_parent.get(branch_id, [])
            if set(explicit) != set(implied):
                raise ValueError(
                    "stat_branch child columns disagree with parent links: "
                    f"branch_id={branch_id}, explicit_children={explicit}, "
                    f"parent_link_children={implied}"
                )
            if explicit:
                children[branch_id] = explicit
    else:
        for branch_id, implied in children_from_parent.items():
            children[branch_id] = sorted(implied)

    root = roots[0]
    visit_state = {}

    def visit(node):
        state = visit_state.get(node, 0)
        if state == 1:
            raise ValueError(f"stat_branch contains a cycle involving branch_id={node}")
        if state == 2:
            return
        visit_state[node] = 1
        for child in children.get(node, []):
            visit(child)
        visit_state[node] = 2

    visit(root)
    unreachable = sorted(set(by_id) - set(visit_state))
    if unreachable:
        raise ValueError(
            "stat_branch contains branches disconnected from its root: "
            f"root={root}, branches={unreachable}"
        )
    return by_id, children, root


def depth_first_tip_order(by_id, children, root):
    ordered = []

    def visit(node):
        row = by_id[node]
        if str(row.get("so_event") or "") == "L" or not children.get(node):
            ordered.append(node)
            return
        for child in children[node]:
            visit(child)

    visit(root)
    return ordered


def ancestor_chain(by_id, node):
    chain = []
    seen = set()
    while node in by_id and node not in seen:
        chain.append(node)
        seen.add(node)
        parent = _int_value(by_id[node].get("parent"))
        if parent < 0:
            break
        node = parent
    return chain


def mrca_node(ancestor_by_node, left, right):
    right_ancestors = set(ancestor_by_node[right])
    for node in ancestor_by_node[left]:
        if node in right_ancestors:
            return node
    raise ValueError(f"No MRCA was found for branch IDs {left} and {right}")


def mapped_species_node_for_gene_node(row, is_query_tip=False):
    primary_field = "spnode_coverage" if is_query_tip else "spnode_generax"
    mapped_species_node = str(row.get(primary_field) or "").strip()
    if not mapped_species_node or mapped_species_node.lower() == "nan":
        mapped_species_node = str(row.get("spnode_coverage") or "").strip()
    if mapped_species_node.lower() == "nan":
        return ""
    return mapped_species_node


def annotate_displayed_duplications(by_id, children, root, glyphs, tree_nodes):
    """Record displayed descendants of original D nodes; deduplicate gene identities."""
    displayed = {
        gene for glyph in glyphs for gene in str(glyph.get("gene_ids") or "").split(";")
        if gene
    }
    tip_by_gene = {}
    for node, row in by_id.items():
        if str(row.get("so_event") or "") != "L":
            continue
        gene = str(row.get("node_name") or node)
        if gene in tip_by_gene:
            raise ValueError(f"Displayed duplication counts require unique saved tip gene IDs: {gene}")
        tip_by_gene[gene] = node
    missing = displayed - tip_by_gene.keys()
    if missing:
        raise ValueError(f"Displayed glyph genes are absent from the saved tree: {sorted(missing)}")
    displayed_by_node = {}
    stack = [(root, False)]
    while stack:
        node, visited = stack.pop()
        node_children = children.get(node, [])
        if not visited and node_children:
            stack.append((node, True))
            stack.extend((child, False) for child in reversed(node_children))
            continue
        if node_children:
            displayed_by_node[node] = set().union(*(displayed_by_node[child] for child in node_children))
        else:
            gene = str(by_id[node].get("node_name") or node)
            displayed_by_node[node] = {gene} if gene in displayed else set()
    compact_roots = [row for row in tree_nodes
                     if row["in_reference_tree"] == 1 and row["parent_node_id"] == ""]
    if len(compact_roots) != 1:
        raise ValueError("Displayed duplication counts require exactly one compact tree root")
    for row in tree_nodes:
        row["displayed_gene_ids"] = ""
        row["displayed_child1_gene_ids"] = ""
        row["displayed_child2_gene_ids"] = ""
        if row["event"] != "D":
            continue
        node_children = children.get(int(row["node_id"]), [])
        if len(node_children) != 2:
            raise ValueError(f"Displayed duplication counts require two children at D node {row['node_id']}")
        for number, child in enumerate(node_children, start=1):
            row[f"displayed_child{number}_gene_ids"] = ";".join(sorted(displayed_by_node[child]))
    compact_roots[0]["displayed_gene_ids"] = ";".join(sorted(displayed))


def refresh_displayed_duplications(store, glyphs, tree_nodes):
    """Refresh provenance after adding candidates, including aliased manifest sources."""
    glyphs_by_family = defaultdict(list)
    nodes_by_family = defaultdict(list)
    for glyph in glyphs:
        glyphs_by_family[str(glyph["family_id"])].append(glyph)
    for node in tree_nodes:
        nodes_by_family[str(node["family_id"])].append(node)
    for family_id, nodes in nodes_by_family.items():
        by_id, children, root = build_tree_index(read_stat_branch(store, family_id))
        annotate_displayed_duplications(by_id, children, root, glyphs_by_family[family_id], nodes)


def build_query_tree_nodes(
    by_id,
    children,
    query_tip_by_id,
    ordered_query_ids,
    definition_by_id,
    global_order,
    family_id,
    family_order,
):
    if not ordered_query_ids:
        return []
    query_tips = [query_tip_by_id[query_id] for query_id in ordered_query_ids]
    ancestor_by_tip = {tip: ancestor_chain(by_id, tip) for tip in query_tips}
    common_ancestors = set(ancestor_by_tip[query_tips[0]])
    for tip in query_tips[1:]:
        common_ancestors.intersection_update(ancestor_by_tip[tip])
    subtree_root = next(
        node for node in ancestor_by_tip[query_tips[0]] if node in common_ancestors
    )

    included = {subtree_root}
    for tip in query_tips:
        for node in ancestor_by_tip[tip]:
            included.add(node)
            if node == subtree_root:
                break
    included_children = {
        node: [child for child in children.get(node, []) if child in included]
        for node in included
    }
    retained = set(query_tips)
    retained.update(
        node for node, node_children in included_children.items() if len(node_children) >= 2
    )
    retained.add(subtree_root)

    query_id_by_tip = {query_tip_by_id[query_id]: query_id for query_id in ordered_query_ids}
    all_duplication_nodes = sorted(
        node for node, row in by_id.items() if str(row.get("so_event") or "") == "D"
    )
    duplication_index_by_node = {
        node: duplication_index
        for duplication_index, node in enumerate(all_duplication_nodes, start=1)
    }
    tree_nodes = []
    for node in sorted(retained):
        parent = _int_value(by_id[node].get("parent"))
        while parent >= 0 and parent not in retained:
            parent = _int_value(by_id[parent].get("parent"))
        query_id = query_id_by_tip.get(node, "")
        mapped_species_node = mapped_species_node_for_gene_node(
            by_id[node], is_query_tip=bool(query_id)
        )
        tree_nodes.append(
            {
                "family_id": family_id,
                "family_order": family_order,
                "node_id": node,
                "parent_node_id": parent if parent in retained else "",
                "is_tip": int(bool(query_id)),
                "event": str(by_id[node].get("so_event") or ""),
                "cds_fasta_id": definition_by_id[query_id]["cds_fasta_id"] if query_id else "",
                "gene_id": definition_by_id[query_id]["gene_id"] if query_id else "",
                "column_order": global_order[query_id] if query_id else "",
                "node_height": 0,
                "plot_order": global_order[query_id] if query_id else "",
                "mapped_species_node": mapped_species_node,
                "duplication_index": duplication_index_by_node.get(node, ""),
                "in_reference_tree": 1,
            }
        )
    tree_node_by_id = {tree_node["node_id"]: tree_node for tree_node in tree_nodes}
    tree_children = defaultdict(list)
    for tree_node in tree_nodes:
        parent = tree_node["parent_node_id"]
        if parent != "":
            tree_children[parent].append(tree_node["node_id"])

    def set_layout(node):
        tree_node = tree_node_by_id[node]
        node_children = tree_children.get(node, [])
        if not node_children:
            return int(tree_node["node_height"]), float(tree_node["plot_order"])
        child_layout = [set_layout(child) for child in node_children]
        tree_node["node_height"] = max(height for height, _ in child_layout) + 1
        tree_node["plot_order"] = sum(order for _, order in child_layout) / len(child_layout)
        return int(tree_node["node_height"]), float(tree_node["plot_order"])

    set_layout(subtree_root)
    for node in all_duplication_nodes:
        if node in retained:
            continue
        parent = _int_value(by_id[node].get("parent"))
        tree_nodes.append(
            {
                "family_id": family_id,
                "family_order": family_order,
                "node_id": node,
                "parent_node_id": parent if parent >= 0 else "",
                "is_tip": 0,
                "event": "D",
                "cds_fasta_id": "",
                "gene_id": "",
                "column_order": "",
                "node_height": "",
                "plot_order": "",
                "mapped_species_node": mapped_species_node_for_gene_node(by_id[node]),
                "duplication_index": duplication_index_by_node[node],
                "in_reference_tree": 0,
            }
        )
    return tree_nodes


def assign_lanes(glyphs):
    by_species_family = defaultdict(list)
    for glyph in glyphs:
        key = (glyph["species"], glyph.get("family_id", ""))
        by_species_family[key].append(glyph)
    for species_glyphs in by_species_family.values():
        lane_ends = []
        ordered = sorted(
            species_glyphs,
            key=lambda glyph: (glyph["start_order"], glyph["end_order"], glyph["family_order"]),
        )
        for glyph in ordered:
            lane_index = None
            for idx, end_order in enumerate(lane_ends):
                if end_order < glyph["start_order"]:
                    lane_index = idx
                    lane_ends[idx] = glyph["end_order"]
                    break
            if lane_index is None:
                lane_index = len(lane_ends)
                lane_ends.append(glyph["end_order"])
            glyph["lane_index"] = lane_index + 1
        lane_count = max(1, len(lane_ends))
        for glyph in species_glyphs:
            glyph["lane_count"] = lane_count


def collect_reference_synteny_evidence(store, columns, glyphs):
    column_by_reference = {
        (str(column["family_id"]), str(column["cds_fasta_id"])): column
        for column in columns
    }
    family_cache = {}
    evidence_rows = []
    for glyph in glyphs:
        family_id = str(glyph["family_id"])
        if family_id not in family_cache:
            synteny_rows = read_family_synteny(store, family_id)
            neighborhoods, window_radius = build_synteny_neighborhoods(synteny_rows)
            family_cache[family_id] = (neighborhoods, window_radius)
        neighborhoods, window_radius = family_cache[family_id]
        candidate_ids = [
            value for value in str(glyph.get("gene_ids") or "").split(";") if value
        ]
        reference_ids = [
            value
            for value in str(glyph.get("reference_cds_fasta_ids") or "").split(";")
            if value
        ]
        for reference_cds_fasta_id in reference_ids:
            column_key = (family_id, reference_cds_fasta_id)
            if column_key not in column_by_reference:
                raise ValueError(
                    "Ortholog glyph references a CDS FASTA ID absent from the column table: "
                    f"family={family_id}, reference={reference_cds_fasta_id}"
                )
            column = column_by_reference[column_key]
            reference_rows = neighborhoods.get(reference_cds_fasta_id, [])
            for candidate_cds_fasta_id in candidate_ids:
                row = {
                    "family_id": family_id,
                    "family_order": glyph["family_order"],
                    "species": glyph["species"],
                    "reference_species": glyph["reference_species"],
                    "relation": glyph["relation"],
                    "reference_cds_fasta_id": reference_cds_fasta_id,
                    "reference_gene_id": column["gene_id"],
                    "column_order": column["column_order"],
                    "candidate_cds_fasta_id": candidate_cds_fasta_id,
                    "glyph_copy_number": glyph["copy_number"],
                    "glyph_start_order": glyph["start_order"],
                    "glyph_end_order": glyph["end_order"],
                    "glyph_lane_index": glyph["lane_index"],
                    "glyph_lane_count": glyph["lane_count"],
                    "synteny_status": "not_evaluable",
                    "support_min_anchor_count": SYNTENY_SUPPORT_MIN_ANCHORS,
                    "synteny_window_radius": window_radius if window_radius > 0 else "",
                    "reference_neighbor_count": len(
                        {int(value["offset"]) for value in reference_rows}
                    ),
                    "candidate_neighbor_count": 0,
                    "flank_coverage": "",
                    "shared_anchor_count": "",
                    "local_synteny_score": "",
                    "collinear_anchor_count": "",
                    "collinearity_ratio": "",
                    "collinear_orientation": "not_evaluable",
                    "shared_group_ids": "",
                }
                if candidate_cds_fasta_id == reference_cds_fasta_id:
                    row["synteny_status"] = "reference_self"
                    row["collinear_orientation"] = "reference_self"
                else:
                    candidate_rows = neighborhoods.get(candidate_cds_fasta_id, [])
                    row["candidate_neighbor_count"] = len(
                        {int(value["offset"]) for value in candidate_rows}
                    )
                    if window_radius > 0 and reference_rows and candidate_rows:
                        metrics = local_synteny_metrics(
                            reference_rows=reference_rows,
                            candidate_rows=candidate_rows,
                            window_radius=window_radius,
                        )
                        row.update(metrics)
                        shared_anchor_count = int(metrics["shared_anchor_count"])
                        if shared_anchor_count >= SYNTENY_SUPPORT_MIN_ANCHORS:
                            row["synteny_status"] = "supported"
                        elif shared_anchor_count == 1:
                            row["synteny_status"] = "single_anchor"
                        else:
                            row["synteny_status"] = "no_support"
                evidence_rows.append(row)
    evidence_rows.sort(
        key=lambda row: (
            int(row["family_order"]),
            str(row["species"]),
            int(row["glyph_start_order"]),
            int(row["column_order"]),
            str(row["candidate_cds_fasta_id"]),
        )
    )
    return evidence_rows


def normalized_ufboot_by_branch(by_id, family_id):
    """Return branch UFBoot values normalized to the 0-100 percentage scale."""

    def has_value(row, field):
        value_text = str(row.get(field) or "").strip()
        return bool(value_text) and value_text.lower() not in {"na", "nan"}

    support_field = (
        "support_generax_ufboot"
        if any(has_value(row, "support_generax_ufboot") for row in by_id.values())
        else "support_unrooted"
    )
    raw_values = {}
    observed = []
    for branch_id, row in by_id.items():
        value_text = str(row.get(support_field) or "").strip()
        if not value_text or value_text.lower() in {"na", "nan"}:
            raw_values[branch_id] = None
            continue
        try:
            value = float(value_text)
        except ValueError as exc:
            raise ValueError(
                f"stat_branch contains a non-numeric {support_field} value: "
                f"family={family_id}, branch_id={branch_id}, value={value_text!r}"
            ) from exc
        if not math.isfinite(value) or value < 0:
            raise ValueError(
                f"stat_branch {support_field} values must be finite and non-negative: "
                f"family={family_id}, branch_id={branch_id}, value={value_text!r}"
            )
        raw_values[branch_id] = value
        observed.append(value)

    # The explicit GeneRax field is produced by IQ-TREE as a percentage.  Do
    # not reinterpret a legitimate 1% value as a 0-1 proportion.  Retain the
    # legacy heuristic only for the generic support_unrooted field.
    multiplier = (
        100.0
        if support_field == "support_unrooted" and observed and max(observed) <= 1.0001
        else 1.0
    )
    normalized = {}
    for branch_id, value in raw_values.items():
        if value is None:
            normalized[branch_id] = None
            continue
        percent = value * multiplier
        if percent > 100.0001:
            raise ValueError(
                f"stat_branch {support_field} values exceed the UFBoot percentage range: "
                f"family={family_id}, branch_id={branch_id}, normalized_value={percent}"
            )
        normalized[branch_id] = round(min(100.0, percent), 10)
    return normalized, support_field


def validate_glyph_ufboot_evidence(rows, glyph):
    """Require one orthology-defining speciation branch per plotted glyph."""

    family_id = str(glyph["family_id"])
    species = str(glyph["species"])
    glyph_span = f"{glyph['start_order']}-{glyph['end_order']}"
    if not rows:
        raise ValueError(
            "Ortholog glyph produced no UFBoot evidence rows: "
            f"family={family_id}, species={species}, glyph_span={glyph_span}"
        )

    statuses = {str(row["orthology_ufboot_status"]) for row in rows}
    if glyph.get("relation") == "weak_duplication":
        if statuses != {"not_evaluable"} or any(
            row["orthology_mrca_event"] != "D"
            or row["orthology_ufboot_unavailable_reason"] != "weak_duplication"
            or row["decisive_branch_ufboot"] != ""
            for row in rows
        ):
            raise ValueError("Weak-duplication glyphs cannot report speciation-based orthology support")
        if len({row["orthology_mrca_branch_id"] for row in rows}) != 1:
            raise ValueError("Weak-duplication glyph pairs must share one duplication MRCA")
        return
    if "reference_self" in statuses:
        if statuses != {"reference_self"}:
            raise ValueError(
                "Ortholog glyph cannot mix reference-self and non-self UFBoot evidence: "
                f"family={family_id}, species={species}, glyph_span={glyph_span}"
            )
        return

    mrca_branch_ids = {
        int(row["orthology_mrca_branch_id"])
        for row in rows
        if row["orthology_mrca_branch_id"] != ""
    }
    if len(mrca_branch_ids) != 1:
        raise ValueError(
            "Ortholog glyph candidate/reference pairs do not share one "
            "orthology-defining speciation branch: "
            f"family={family_id}, species={species}, glyph_span={glyph_span}, "
            f"mrca_branch_ids={sorted(mrca_branch_ids)}"
        )

    unavailable_reasons = {
        str(row["orthology_ufboot_unavailable_reason"]) for row in rows
    }
    if len(statuses) != 1 or len(unavailable_reasons) != 1:
        raise ValueError(
            "Ortholog glyph candidate/reference pairs do not share one UFBoot "
            "availability state: "
            f"family={family_id}, species={species}, glyph_span={glyph_span}"
        )

    if statuses == {"evaluated"}:
        support_values = {
            float(row["decisive_branch_ufboot"])
            for row in rows
            if row["decisive_branch_ufboot"] != ""
        }
        if len(support_values) != 1:
            raise ValueError(
                "Ortholog glyph candidate/reference pairs do not share one "
                "orthology-defining branch UFBoot value: "
                f"family={family_id}, species={species}, glyph_span={glyph_span}, "
                f"ufboot_values={sorted(support_values)}"
            )


def collect_reference_ufboot_evidence(
    store,
    columns,
    glyphs,
    require_one_branch_per_glyph=True,
):
    """Collect orthology-defining branch UFBoot for each plotted glyph.

    The current orthology assignment includes a non-reference tip in a reference
    column when their MRCA is reconciled as a speciation. The branch entering
    that MRCA is therefore the gene-tree branch most directly associated with
    the assignment. Every pair represented by one glyph must resolve to that
    same branch; pairwise rows are retained for provenance. Root MRCAs have no
    corresponding unrooted bipartition and are retained as not evaluable.
    """

    column_by_reference = {
        (str(column["family_id"]), str(column["cds_fasta_id"])): column
        for column in columns
    }
    family_cache = {}
    evidence_rows = []
    for glyph in glyphs:
        family_id = str(glyph["family_id"])
        if family_id not in family_cache:
            rows = read_stat_branch(store, family_id)
            if not rows:
                raise ValueError(
                    "Orthology UFBoot evidence requires stat_branch rows: "
                    f"family={family_id}"
                )
            by_id, _children, root = build_tree_index(rows)
            ancestor_by_node = {
                node: ancestor_chain(by_id, node) for node in by_id
            }
            leaf_by_name = {}
            for branch_id, row in by_id.items():
                if str(row.get("so_event") or "") != "L":
                    continue
                node_name = str(row.get("node_name") or "").strip()
                if not node_name:
                    continue
                if node_name in leaf_by_name:
                    raise ValueError(
                        "stat_branch contains duplicate leaf node_name values: "
                        f"family={family_id}, node_name={node_name!r}"
                    )
                leaf_by_name[node_name] = branch_id
            ufboot_by_branch, ufboot_support_source = normalized_ufboot_by_branch(
                by_id, family_id
            )
            family_cache[family_id] = (
                by_id,
                root,
                ancestor_by_node,
                leaf_by_name,
                ufboot_by_branch,
                ufboot_support_source,
            )
        (
            by_id,
            root,
            ancestor_by_node,
            leaf_by_name,
            ufboot_by_branch,
            ufboot_support_source,
        ) = family_cache[family_id]
        candidate_ids = [
            value for value in str(glyph.get("gene_ids") or "").split(";") if value
        ]
        reference_ids = [
            value
            for value in str(glyph.get("reference_cds_fasta_ids") or "").split(";")
            if value
        ]
        glyph_evidence_rows = []
        for reference_cds_fasta_id in reference_ids:
            column_key = (family_id, reference_cds_fasta_id)
            if column_key not in column_by_reference:
                raise ValueError(
                    "Ortholog glyph references a CDS FASTA ID absent from the column table: "
                    f"family={family_id}, reference={reference_cds_fasta_id}"
                )
            column = column_by_reference[column_key]
            reference_branch_id = int(column["reference_tip_branch_id"])
            if reference_branch_id not in by_id:
                raise ValueError(
                    "Reference-gene branch is absent from stat_branch: "
                    f"family={family_id}, branch_id={reference_branch_id}"
                )
            for candidate_cds_fasta_id in candidate_ids:
                row = {
                    "family_id": family_id,
                    "family_order": glyph["family_order"],
                    "species": glyph["species"],
                    "reference_species": glyph["reference_species"],
                    "relation": glyph["relation"],
                    "reference_cds_fasta_id": reference_cds_fasta_id,
                    "reference_gene_id": column["gene_id"],
                    "column_order": column["column_order"],
                    "candidate_cds_fasta_id": candidate_cds_fasta_id,
                    "glyph_copy_number": glyph["copy_number"],
                    "glyph_start_order": glyph["start_order"],
                    "glyph_end_order": glyph["end_order"],
                    "glyph_lane_index": glyph["lane_index"],
                    "glyph_lane_count": glyph["lane_count"],
                    "orthology_mrca_branch_id": "",
                    "orthology_mrca_event": "",
                    "ufboot_support_source": ufboot_support_source,
                    "decisive_branch_ufboot": "",
                    "orthology_ufboot_status": "not_evaluable",
                    "orthology_ufboot_unavailable_reason": "missing_support",
                }
                if candidate_cds_fasta_id == reference_cds_fasta_id:
                    row["orthology_ufboot_status"] = "reference_self"
                    row["orthology_ufboot_unavailable_reason"] = "reference_self"
                else:
                    candidate_branch_id = leaf_by_name.get(candidate_cds_fasta_id)
                    if candidate_branch_id is None:
                        raise ValueError(
                            "Ortholog glyph candidate is absent from stat_branch leaves: "
                            f"family={family_id}, candidate={candidate_cds_fasta_id!r}"
                        )
                    mrca = mrca_node(
                        ancestor_by_node,
                        candidate_branch_id,
                        reference_branch_id,
                    )
                    mrca_event = str(by_id[mrca].get("so_event") or "")
                    row["orthology_mrca_branch_id"] = mrca
                    row["orthology_mrca_event"] = mrca_event
                    is_weak = glyph.get("relation") == "weak_duplication"
                    if is_weak and mrca_event != "D":
                        raise ValueError("Weak-duplication glyph pair must retain its original D MRCA")
                    if not is_weak and mrca_event != "S":
                        raise ValueError(
                            "Ortholog glyph pair does not have a speciation MRCA: "
                            f"family={family_id}, candidate={candidate_cds_fasta_id!r}, "
                            f"reference={reference_cds_fasta_id!r}, mrca={mrca}, "
                            f"event={mrca_event!r}"
                        )
                    if is_weak:
                        row["orthology_ufboot_unavailable_reason"] = "weak_duplication"
                    elif mrca == root:
                        row["orthology_ufboot_unavailable_reason"] = "mrca_is_root"
                    else:
                        support = ufboot_by_branch.get(mrca)
                        if support is not None:
                            row["decisive_branch_ufboot"] = support
                            row["orthology_ufboot_status"] = "evaluated"
                            row["orthology_ufboot_unavailable_reason"] = ""
                glyph_evidence_rows.append(row)
        if require_one_branch_per_glyph or glyph.get("relation") == "weak_duplication":
            validate_glyph_ufboot_evidence(glyph_evidence_rows, glyph)
        evidence_rows.extend(glyph_evidence_rows)
    evidence_rows.sort(
        key=lambda row: (
            int(row["family_order"]),
            str(row["species"]),
            int(row["glyph_start_order"]),
            int(row["column_order"]),
            str(row["candidate_cds_fasta_id"]),
        )
    )
    return evidence_rows


def _collect_family_orthologs_for_anchors(
    rows,
    cds_fasta_ids,
    query_tip_by_id,
    definitions,
    anchor_basis_label,
    family_id,
    family_order,
    first_column_order,
    preserve_query_label=False,
    enforce_identity_species="",
):
    by_id, children, root = build_tree_index(rows)
    if not query_tip_by_id:
        return [], [], []

    definition_order = {definition["query_id"]: idx for idx, definition in enumerate(definitions)}
    tip_position = {tip: idx for idx, tip in enumerate(depth_first_tip_order(by_id, children, root))}
    ordered_query_ids = sorted(
        query_tip_by_id,
        key=lambda query_id: (tip_position.get(query_tip_by_id[query_id], 10**12), definition_order[query_id]),
    )
    definition_by_id = resolve_query_cds_definitions(
        by_id=by_id,
        query_tip_by_id=query_tip_by_id,
        definitions=definitions,
        cds_fasta_ids=cds_fasta_ids,
        family_id=family_id,
        preserve_query_label=preserve_query_label,
    )
    local_order = {query_id: idx + 1 for idx, query_id in enumerate(ordered_query_ids)}
    global_order = {
        query_id: first_column_order + idx for idx, query_id in enumerate(ordered_query_ids)
    }

    columns = []
    for query_id in ordered_query_ids:
        definition = definition_by_id[query_id]
        columns.append(
            {
                "column_order": global_order[query_id],
                "family_id": family_id,
                "family_order": family_order,
                "cds_fasta_id": definition["cds_fasta_id"],
                "gene_id": definition["gene_id"],
                "plot_label": definition["query_label"],
                "reference_tip_branch_id": query_tip_by_id[query_id],
            }
        )

    tree_nodes = build_query_tree_nodes(
        by_id=by_id,
        children=children,
        query_tip_by_id=query_tip_by_id,
        ordered_query_ids=ordered_query_ids,
        definition_by_id=definition_by_id,
        global_order=global_order,
        family_id=family_id,
        family_order=family_order,
    )

    ancestor_by_node = {node: ancestor_chain(by_id, node) for node in by_id}
    grouped = defaultdict(list)
    for branch_id, row in by_id.items():
        if str(row.get("so_event") or "") != "L":
            continue
        species = str(row.get("spnode_coverage") or "").strip()
        if not species:
            continue
        query_set = []
        for query_id in ordered_query_ids:
            query_tip = query_tip_by_id[query_id]
            ancestor = mrca_node(ancestor_by_node, branch_id, query_tip)
            if branch_id == query_tip or str(by_id[ancestor].get("so_event") or "") == "S":
                query_set.append(query_id)
        if query_set:
            grouped[(species, tuple(query_set))].append(str(row.get("node_name") or branch_id))

    glyphs = []
    for (species, query_set), gene_ids in grouped.items():
        local_positions = [local_order[query_id] for query_id in query_set]
        start_local = min(local_positions)
        end_local = max(local_positions)
        is_contiguous = (end_local - start_local + 1) == len(query_set)
        relation = "specific" if len(query_set) == 1 else "shared_ancestral"
        if not is_contiguous:
            relation = "ambiguous"
        glyphs.append(
            {
                "species": species,
                "family_id": family_id,
                "family_order": family_order,
                "relation": relation,
                "reference_cds_fasta_ids": ";".join(
                    definition_by_id[query_id]["cds_fasta_id"] for query_id in query_set
                ),
                "reference_gene_ids": ";".join(
                    definition_by_id[query_id]["gene_id"] for query_id in query_set
                ),
                "reference_gene_count": len(query_set),
                "copy_number": len(gene_ids),
                "gene_ids": ";".join(sorted(gene_ids)),
                "start_order": min(global_order[query_id] for query_id in query_set),
                "end_order": max(global_order[query_id] for query_id in query_set),
                "is_contiguous": int(is_contiguous),
                "lane_index": 1,
                "lane_count": 1,
            }
        )
    normalized_anchor_basis_label = normalize_species_label(anchor_basis_label)
    if enforce_identity_species:
        normalized_reference_species = normalize_species_label(enforce_identity_species)
        reference_glyphs = [
            glyph
            for glyph in glyphs
            if normalize_species_label(glyph["species"]) == normalized_reference_species
        ]
        observed_reference_cds_ids = []
        reference_glyph_is_one_to_one = len(reference_glyphs) == len(ordered_query_ids)
        for glyph in reference_glyphs:
            glyph_reference_ids = [
                value
                for value in str(glyph["reference_cds_fasta_ids"]).split(";")
                if value
            ]
            glyph_gene_ids = [
                value for value in str(glyph["gene_ids"]).split(";") if value
            ]
            observed_reference_cds_ids.extend(glyph_reference_ids)
            reference_glyph_is_one_to_one = reference_glyph_is_one_to_one and (
                glyph["relation"] == "specific"
                and glyph["reference_gene_count"] == 1
                and glyph["copy_number"] == 1
                and len(glyph_reference_ids) == 1
                and glyph_gene_ids == glyph_reference_ids
            )
        expected_reference_cds_ids = sorted(ordered_query_ids)
        reference_glyph_is_one_to_one = reference_glyph_is_one_to_one and (
            sorted(observed_reference_cds_ids) == expected_reference_cds_ids
        )
        if not reference_glyph_is_one_to_one:
            raise ValueError(
                "Reference-species genes did not produce a one-to-one identity row; "
                "the reconciled gene tree is inconsistent with reference-gene columns: "
                f"family={family_id}, reference_species={normalized_reference_species}"
            )
    annotate_displayed_duplications(by_id, children, root, glyphs, tree_nodes)
    for output_rows in (columns, glyphs, tree_nodes):
        for output_row in output_rows:
            output_row["reference_species"] = normalized_anchor_basis_label
    return columns, glyphs, tree_nodes


def collect_family_orthologs(
    rows,
    cds_fasta_ids,
    reference_species,
    family_id,
    family_order,
    first_column_order,
):
    query_tip_by_id, definitions = select_reference_tip_nodes(rows, reference_species)
    return _collect_family_orthologs_for_anchors(
        rows=rows,
        cds_fasta_ids=cds_fasta_ids,
        query_tip_by_id=query_tip_by_id,
        definitions=definitions,
        anchor_basis_label=reference_species,
        family_id=family_id,
        family_order=family_order,
        first_column_order=first_column_order,
        enforce_identity_species=reference_species,
    )


def write_tsv(path, fieldnames, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def collect_query_gene_orthologs(
    dir_gene_family,
    dir_query_gene,
    reference_species,
    family_file="",
):
    query_dir = Path(dir_query_gene)
    family_ids = read_family_ids(query_dir, family_file=family_file)
    store = GeneFamilyOutputStore(dir_gene_family)
    columns = []
    glyphs = []
    tree_nodes = []
    next_column_order = 1
    next_family_order = 1
    for family_id in family_ids:
        rows = read_stat_branch(store, family_id)
        if not rows:
            continue
        cds_fasta_ids = read_family_cds_fasta_ids(store, family_id)
        family_columns, family_glyphs, family_tree_nodes = collect_family_orthologs(
            rows=rows,
            cds_fasta_ids=cds_fasta_ids,
            reference_species=reference_species,
            family_id=family_id,
            family_order=next_family_order,
            first_column_order=next_column_order,
        )
        if not family_columns:
            continue
        columns.extend(family_columns)
        glyphs.extend(family_glyphs)
        tree_nodes.extend(family_tree_nodes)
        next_column_order += len(family_columns)
        next_family_order += 1

    assign_lanes(glyphs)
    columns.sort(key=lambda row: int(row["column_order"]))
    glyphs.sort(
        key=lambda row: (
            str(row["species"]),
            int(row["start_order"]),
            int(row["end_order"]),
            str(row["relation"]),
        )
    )
    tree_nodes.sort(key=lambda row: (int(row["family_order"]), int(row["node_id"])))
    if family_ids and not columns:
        raise ValueError(
            "No genes from the selected reference species were found in the selected "
            f"gene families: reference_species={normalize_species_label(reference_species)}"
        )
    return columns, glyphs, tree_nodes


def collect_family_query_anchor_orthologs(
    rows,
    cds_fasta_ids,
    query_definitions,
    family_id,
    family_order,
    first_column_order,
    first_query_order,
    query_label="id",
):
    """Collect query-anchored orthologs, coalescing records on one tree tip."""

    assignments = select_query_tip_assignments(rows, query_definitions)
    by_id, _children, _root = build_tree_index(rows)
    definitions_by_id = {
        definition["query_id"]: definition for definition in query_definitions
    }
    query_ids_by_tip = defaultdict(list)
    for definition in query_definitions:
        query_id = definition["query_id"]
        query_ids_by_tip[assignments[query_id]["tip"]].append(query_id)

    anchor_tip_by_id = {}
    anchor_definitions = []
    anchor_metadata_by_id = {}
    for tip, query_ids in query_ids_by_tip.items():
        anchor_cds_fasta_id = str(by_id[tip].get("node_name") or "").strip()
        if not anchor_cds_fasta_id:
            raise ValueError(
                "A query marker selected a gene-tree tip without node_name: "
                f"family={family_id}, branch_id={tip}, queries={query_ids}"
            )
        if anchor_cds_fasta_id in anchor_tip_by_id:
            raise ValueError(
                "Query anchors have duplicate CDS FASTA IDs: "
                f"family={family_id}, cds_fasta_id={anchor_cds_fasta_id}"
            )
        query_labels = [definitions_by_id[query_id]["query_label"] for query_id in query_ids]
        marker_sources = [assignments[query_id]["marker_source"] for query_id in query_ids]
        distinct_sources = sorted(set(marker_sources))
        anchor_source = distinct_sources[0] if len(distinct_sources) == 1 else "mixed"
        anchor_species = normalize_species_label(by_id[tip].get("spnode_coverage"))
        if not anchor_species:
            raise ValueError(
                "A query marker selected a gene-tree tip without species coverage: "
                f"family={family_id}, branch_id={tip}, queries={query_ids}"
            )
        plot_label = query_labels[0] if query_label == "label" else query_ids[0]
        if len(query_ids) > 1:
            plot_label = f"{plot_label} (+{len(query_ids) - 1})"
        anchor_tip_by_id[anchor_cds_fasta_id] = tip
        anchor_definitions.append(
            {
                "query_id": anchor_cds_fasta_id,
                "query_label": plot_label,
            }
        )
        anchor_metadata_by_id[anchor_cds_fasta_id] = {
            "query_ids": list(query_ids),
            "query_labels": query_labels,
            "marker_sources": marker_sources,
            "anchor_source": anchor_source,
            "anchor_species": anchor_species,
            "anchor_tip_branch_id": tip,
        }

    columns, glyphs, tree_nodes = _collect_family_orthologs_for_anchors(
        rows=rows,
        cds_fasta_ids=cds_fasta_ids,
        query_tip_by_id=anchor_tip_by_id,
        definitions=anchor_definitions,
        anchor_basis_label="query_gene",
        family_id=family_id,
        family_order=family_order,
        first_column_order=first_column_order,
        preserve_query_label=True,
    )
    column_by_anchor = {}
    for column in columns:
        anchor_cds_fasta_id = str(column["cds_fasta_id"])
        metadata = anchor_metadata_by_id[anchor_cds_fasta_id]
        column.update(
            {
                "query_ids": ";".join(metadata["query_ids"]),
                "query_labels": ";".join(metadata["query_labels"]),
                "query_count": len(metadata["query_ids"]),
                "anchor_source": metadata["anchor_source"],
                "anchor_species": metadata["anchor_species"],
            }
        )
        column_by_anchor[anchor_cds_fasta_id] = column

    for glyph in glyphs:
        anchor_ids = [
            value
            for value in str(glyph["reference_cds_fasta_ids"]).split(";")
            if value
        ]
        glyph["anchor_query_ids"] = ";".join(
            query_id
            for anchor_id in anchor_ids
            for query_id in anchor_metadata_by_id[anchor_id]["query_ids"]
        )

    query_map = []
    for query_index, definition in enumerate(query_definitions, start=first_query_order):
        query_id = definition["query_id"]
        assignment = assignments[query_id]
        tip = assignment["tip"]
        anchor_cds_fasta_id = str(by_id[tip].get("node_name") or "").strip()
        column = column_by_anchor[anchor_cds_fasta_id]
        metadata = anchor_metadata_by_id[anchor_cds_fasta_id]
        query_map.append(
            {
                "family_id": family_id,
                "family_order": family_order,
                "query_order": query_index,
                "query_id": query_id,
                "query_label": definition["query_label"],
                "marker_source": assignment["marker_source"],
                "anchor_species": metadata["anchor_species"],
                "anchor_cds_fasta_id": anchor_cds_fasta_id,
                "anchor_gene_id": column["gene_id"],
                "anchor_tip_branch_id": tip,
                "column_order": column["column_order"],
                "merged_query_count": len(metadata["query_ids"]),
                "source_species": definition.get("source_species", ""),
                "hog_ids": definition.get("hog_ids", ""),
            }
        )
    return columns, glyphs, tree_nodes, query_map


def collect_query_anchor_orthologs(
    dir_gene_family,
    dir_query_gene,
    family_file="",
    manifest_records=None,
    query_metadata=None,
    query_selection="all",
    species_tree="",
    target_species="",
    selection_species="",
    query_label="id",
    selection_audit=None,
    long_rows=None,
):
    query_dir = Path(dir_query_gene)
    if manifest_records is None and not dir_query_gene:
        raise ValueError("--dir_query_gene or --family_manifest is required")
    family_ids = ([record["family_id"] for record in manifest_records] if manifest_records
                  else read_family_ids(query_dir, family_file=family_file))
    store = (ManifestOutputStore(manifest_records) if manifest_records
             else GeneFamilyOutputStore(dir_gene_family))
    manifest_by_id = {record["family_id"]: record for record in manifest_records or []}
    columns = []
    glyphs = []
    tree_nodes = []
    query_map = []
    next_column_order = 1
    next_family_order = 1
    next_query_order = 1
    for family_id in family_ids:
        rows = read_stat_branch(store, family_id)
        if not rows:
            if manifest_records:
                raise ValueError(f"Saved family has no stat_branch artifact: {family_id}")
            continue
        cds_fasta_ids = read_family_cds_fasta_ids(store, family_id)
        if manifest_records:
            tip_ids = [str(row["node_name"]) for row in rows if row.get("so_event") == "L"]
            if len(tip_ids) != len(set(tip_ids)):
                raise ValueError(f"Saved tree repeats CDS FASTA IDs: {family_id}")
            missing = set(tip_ids) - set(cds_fasta_ids)
            if missing:
                raise ValueError(f"Saved tree tips missing from CDS FASTA: {family_id}: {sorted(missing)}")
            query_definitions = manifest_query_definitions(manifest_by_id[family_id], rows, cds_fasta_ids)
        else:
            query_definitions = read_query_definitions(query_dir / family_id)
        if not query_definitions:
            if manifest_records:
                raise ValueError(f"Manifest query file has no records: {family_id}")
            continue
        known_query_ids = {definition["query_id"] for definition in query_definitions}
        unknown_metadata_ids = {query_id for metadata_family, query_id in query_metadata or {}
                                if metadata_family == family_id and query_id not in known_query_ids}
        if unknown_metadata_ids:
            raise ValueError(f"Query metadata IDs absent from selected family {family_id}: {sorted(unknown_metadata_ids)}")
        for definition in query_definitions:
            metadata = (query_metadata or {}).get((family_id, definition["query_id"]), {})
            header_species = definition.get("source_species", "")
            if header_species and metadata.get("source_species") and header_species != metadata["source_species"]:
                raise ValueError(f"Conflicting source species metadata: {family_id}: {definition['query_id']}")
            definition.update({key: value for key, value in metadata.items() if value})
        query_definitions, audit = select_closest_queries(
            rows, query_definitions, family_id, species_tree, target_species,
            selection_species, query_selection,
        )
        if selection_audit is not None:
            selection_audit.extend(audit)
        if long_rows is not None:
            counts = defaultdict(int)
            for row in rows:
                if row.get("so_event") == "L":
                    counts[normalize_species_label(row.get("spnode_coverage"))] += 1
            for species, count in counts.items():
                long_rows.append(dict(species=species, species_display=species.replace("_", " "),
                                      query=family_id, query_order=next_family_order,
                                      presence=int(count > 0), copy_number=count, status="complete"))
        (
            family_columns,
            family_glyphs,
            family_tree_nodes,
            family_query_map,
        ) = collect_family_query_anchor_orthologs(
            rows=rows,
            cds_fasta_ids=cds_fasta_ids,
            query_definitions=query_definitions,
            family_id=family_id,
            family_order=next_family_order,
            first_column_order=next_column_order,
            first_query_order=next_query_order,
            query_label=query_label,
        )
        columns.extend(family_columns)
        glyphs.extend(family_glyphs)
        tree_nodes.extend(family_tree_nodes)
        query_map.extend(family_query_map)
        next_column_order += len(family_columns)
        next_query_order += len(family_query_map)
        next_family_order += 1

    assign_lanes(glyphs)
    columns.sort(key=lambda row: int(row["column_order"]))
    glyphs.sort(
        key=lambda row: (
            str(row["species"]),
            int(row["start_order"]),
            int(row["end_order"]),
            str(row["relation"]),
        )
    )
    tree_nodes.sort(key=lambda row: (int(row["family_order"]), int(row["node_id"])))
    query_map.sort(key=lambda row: int(row["query_order"]))
    if family_ids and not columns:
        raise ValueError(
            "No query records with gene-tree tip markers were found in the selected gene families"
        )
    return columns, glyphs, tree_nodes, query_map


def query_columns_for_output(columns):
    return [
        {
            "column_order": row["column_order"],
            "family_id": row["family_id"],
            "family_order": row["family_order"],
            "basis": "query_gene",
            "query_ids": row["query_ids"],
            "query_labels": row["query_labels"],
            "query_count": row["query_count"],
            "anchor_source": row["anchor_source"],
            "anchor_species": row["anchor_species"],
            "anchor_cds_fasta_id": row["cds_fasta_id"],
            "anchor_gene_id": row["gene_id"],
            "plot_label": row["plot_label"],
            "anchor_tip_branch_id": row["reference_tip_branch_id"],
        }
        for row in columns
    ]


def query_glyphs_for_output(glyphs):
    return [
        {
            "species": row["species"],
            "family_id": row["family_id"],
            "family_order": row["family_order"],
            "basis": "query_gene",
            "relation": row["relation"],
            "anchor_cds_fasta_ids": row["reference_cds_fasta_ids"],
            "anchor_gene_ids": row["reference_gene_ids"],
            "anchor_query_ids": row["anchor_query_ids"],
            "anchor_count": row["reference_gene_count"],
            "copy_number": row["copy_number"],
            "gene_ids": row["gene_ids"],
            "start_order": row["start_order"],
            "end_order": row["end_order"],
            "is_contiguous": row["is_contiguous"],
            "lane_index": row["lane_index"],
            "lane_count": row["lane_count"],
        }
        for row in glyphs
    ]


def query_tree_for_output(tree_nodes):
    output = []
    for row in tree_nodes:
        converted = dict(row)
        converted.pop("reference_species", None)
        converted["basis"] = "query_gene"
        converted["anchor_cds_fasta_id"] = converted.pop("cds_fasta_id")
        converted["anchor_gene_id"] = converted.pop("gene_id")
        converted["in_anchor_tree"] = converted.pop("in_reference_tree")
        output.append(converted)
    return output


def query_evidence_for_output(rows):
    output = []
    for row in rows:
        converted = {}
        for field, value in row.items():
            if field == "reference_species":
                converted["basis"] = "query_gene"
            elif field.startswith("reference_"):
                converted[field.replace("reference_", "anchor_", 1)] = value
            else:
                converted[field] = value
        if converted.get("synteny_status") == "reference_self":
            converted["synteny_status"] = "anchor_self"
        if converted.get("collinear_orientation") == "reference_self":
            converted["collinear_orientation"] = "anchor_self"
        if converted.get("orthology_ufboot_status") == "reference_self":
            converted["orthology_ufboot_status"] = "anchor_self"
        if converted.get("orthology_ufboot_unavailable_reason") == "reference_self":
            converted["orthology_ufboot_unavailable_reason"] = "anchor_self"
        output.append(converted)
    return output


def run(args):
    dup_threshold = validate_dup_conf_threshold(getattr(args, "dup_conf_score_threshold", 0))
    if dup_threshold > 0 and not getattr(args, "out_dup_conf", ""):
        raise ValueError("--out_dup_conf is required when --dup_conf_score_threshold is positive")
    manifest_path = getattr(args, "family_manifest", "")
    manifest_records = read_family_manifest(manifest_path) if manifest_path else None
    store = ManifestOutputStore(manifest_records) if manifest_records else GeneFamilyOutputStore(args.dir_gene_family)
    basis = getattr(args, "basis", "reference_species")
    if basis == "reference_species":
        if manifest_records or getattr(args, "query_selection", "all") != "all":
            raise ValueError("Family manifests and query selection require --basis=query_gene")
        if not str(args.reference_species).strip():
            raise ValueError("--reference_species is required when --basis=reference_species")
        columns, glyphs, tree_nodes = collect_query_gene_orthologs(
            dir_gene_family=args.dir_gene_family,
            dir_query_gene=args.dir_query_gene,
            reference_species=args.reference_species,
            family_file=args.family_file,
        )
        glyphs, dup_evidence = add_weak_duplication_candidates(store, columns, glyphs, dup_threshold)
        if dup_evidence:
            refresh_displayed_duplications(store, glyphs, tree_nodes)
        synteny_evidence = collect_reference_synteny_evidence(
            store=store,
            columns=columns,
            glyphs=glyphs,
        )
        ufboot_evidence = collect_reference_ufboot_evidence(
            store=store,
            columns=columns,
            glyphs=glyphs,
        )
        write_tsv(args.out_columns, COLUMN_FIELDS, columns)
        write_tsv(args.out_glyphs, GLYPH_FIELDS, glyphs)
        write_tsv(args.out_tree, TREE_FIELDS, tree_nodes)
        write_tsv(args.out_synteny, SYNTENY_FIELDS, synteny_evidence)
        write_tsv(args.out_ufboot, UFBOOT_FIELDS, ufboot_evidence)
        if getattr(args, "out_dup_conf", ""):
            write_tsv(args.out_dup_conf, DUP_CONF_FIELDS, dup_evidence)
        print(
            "Reference-species ortholog summary: "
            f"reference_species={normalize_species_label(args.reference_species)}, "
            f"columns={len(columns)}, glyphs={len(glyphs)}, tree_nodes={len(tree_nodes)}, "
            f"shared_ancestral={sum(row['relation'] == 'shared_ancestral' for row in glyphs)}, "
            f"synteny_supported={sum(row['synteny_status'] == 'supported' for row in synteny_evidence)}, "
            f"synteny_single_anchor={sum(row['synteny_status'] == 'single_anchor' for row in synteny_evidence)}, "
            f"ufboot_evaluated={sum(row['orthology_ufboot_status'] == 'evaluated' for row in ufboot_evidence)}, "
            f"ufboot_not_evaluable={sum(row['orthology_ufboot_status'] == 'not_evaluable' for row in ufboot_evidence)}",
            flush=True,
        )
        return

    out_query_map = str(getattr(args, "out_query_map", "")).strip()
    if not out_query_map:
        raise ValueError("--out_query_map is required when --basis=query_gene")
    selection_audit = []
    long_rows = []
    selection_output = getattr(args, "out_selection", "")
    if getattr(args, "query_selection", "all") == "closest" and not selection_output:
        raise ValueError("closest requires --out_selection to record excluded queries")
    columns, glyphs, tree_nodes, query_map = collect_query_anchor_orthologs(
        dir_gene_family=args.dir_gene_family,
        dir_query_gene=args.dir_query_gene,
        family_file=args.family_file,
        manifest_records=manifest_records,
        query_metadata=read_query_metadata(getattr(args, "query_metadata", "")),
        query_selection=getattr(args, "query_selection", "all"),
        species_tree=getattr(args, "species_tree", ""),
        target_species=getattr(args, "target_species", ""),
        selection_species=getattr(args, "selection_species", ""),
        query_label=getattr(args, "query_label", "id"),
        selection_audit=selection_audit,
        long_rows=long_rows,
    )
    glyphs, dup_evidence = add_weak_duplication_candidates(store, columns, glyphs, dup_threshold)
    if dup_evidence:
        refresh_displayed_duplications(store, glyphs, tree_nodes)
    synteny_evidence = collect_reference_synteny_evidence(
        store=store,
        columns=columns,
        glyphs=glyphs,
    )
    ufboot_evidence = collect_reference_ufboot_evidence(
        store=store,
        columns=columns,
        glyphs=glyphs,
        require_one_branch_per_glyph=False,
    )
    write_tsv(args.out_columns, QUERY_COLUMN_FIELDS, query_columns_for_output(columns))
    write_tsv(args.out_glyphs, QUERY_GLYPH_FIELDS, query_glyphs_for_output(glyphs))
    write_tsv(args.out_tree, QUERY_TREE_FIELDS, query_tree_for_output(tree_nodes))
    write_tsv(
        args.out_synteny,
        QUERY_SYNTENY_FIELDS,
        query_evidence_for_output(synteny_evidence),
    )
    write_tsv(
        args.out_ufboot,
        QUERY_UFBOOT_FIELDS,
        query_evidence_for_output(ufboot_evidence),
    )
    write_tsv(out_query_map, QUERY_MAP_FIELDS, query_map)
    if getattr(args, "out_dup_conf", ""):
        write_tsv(args.out_dup_conf, QUERY_DUP_CONF_FIELDS, query_evidence_for_output(dup_evidence))
    if selection_output:
        write_tsv(selection_output, SELECTION_FIELDS, selection_audit)
    if getattr(args, "out_long", ""):
        # Complete absent-family rows for the species-tree tips as well as observed taxa.
        species = {row["species"] for row in long_rows}
        if getattr(args, "species_tree", ""):
            from Bio import Phylo

            species.update(normalize_species_label(tip.name)
                           for tip in Phylo.read(args.species_tree, "newick").get_terminals())
        by_key = {(row["query"], row["species"]): row for row in long_rows}
        families = sorted({(row["query_order"], row["query"]) for row in long_rows})
        complete_long = []
        for order, family in families:
            for sp in sorted(species):
                complete_long.append(by_key.get((family, sp), dict(
                    species=sp, species_display=sp.replace("_", " "), query=family,
                    query_order=order, presence=0, copy_number=0, status="complete")))
        write_tsv(args.out_long, ["species", "species_display", "query", "query_order",
                                  "presence", "copy_number", "status"], complete_long)
    gene_families = defaultdict(set)
    for glyph in glyphs:
        for gene in glyph["gene_ids"].split(";"):
            gene_families[(glyph["species"], gene)].add(glyph["family_id"])
    overlaps = [dict(species=sp, cds_fasta_id=gene, family_ids=";".join(sorted(families)),
                     family_count=len(families))
                for (sp, gene), families in sorted(gene_families.items()) if len(families) > 1]
    if getattr(args, "out_overlap", ""):
        write_tsv(args.out_overlap, ["species", "cds_fasta_id", "family_ids", "family_count"], overlaps)
    if overlaps:
        warnings.warn(f"{len(overlaps)} genes occur in multiple plot blocks; do not sum blocks as independent counts",
                      stacklevel=2)
    print(
        "Query-gene ortholog summary: "
        f"query_records={len(query_map)}, anchors={len(columns)}, "
        f"coalesced_records={len(query_map) - len(columns)}, glyphs={len(glyphs)}, "
        f"tree_nodes={len(tree_nodes)}, "
        f"shared_across_anchors={sum(row['relation'] == 'shared_ancestral' for row in glyphs)}, "
        f"synteny_supported={sum(row['synteny_status'] == 'supported' for row in synteny_evidence)}, "
        f"ufboot_evaluated={sum(row['orthology_ufboot_status'] == 'evaluated' for row in ufboot_evidence)}",
        flush=True,
    )


def main():
    run(build_arg_parser().parse_args())


if __name__ == "__main__":
    main()
