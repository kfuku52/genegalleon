#!/usr/bin/env python3
"""Summarize reference-species genes and their ortholog copies.

A gene-tree tip belongs to a reference-gene column when its MRCA with that
reference-species tip is a speciation node (or when it is the reference tip
itself). A tip assigned to more than one adjacent column is represented as a
shared ancestral copy predating the duplication that separated those
reference genes. Tips are grouped by species and reference-gene set so their
group size is the plotted copy number. Plot labels are verified against the
family CDS FASTA. The tree table also retains every reconciled duplication
outside the compact reference-gene subtree so the plot can map the complete
family's duplication history onto the species tree. A fourth table compares
each plotted copy with each covered reference gene using the family's local
synteny neighborhoods. Two or more distinct shared neighbor-similarity groups
provide local-synteny support; a single shared group is retained as
single-anchor evidence but is not called supported.
A fifth table retains candidate/reference pair provenance for Gene tree UFBoot,
while requiring every pair represented by one glyph to resolve to the same
orthology-defining speciation branch and support value.
"""

import argparse
import csv
import gzip
import io
import math
from collections import defaultdict
from pathlib import Path

from gene_family_output_store import GeneFamilyOutputStore

STAT_BRANCH_SUFFIX = "_stat.branch.tsv"
SYNTENY_SUFFIX = "_synteny.tsv"
SYNTENY_SUPPORT_MIN_ANCHORS = 2
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


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir_gene_family", metavar="PATH", required=True)
    parser.add_argument("--dir_query_gene", metavar="PATH", required=True)
    parser.add_argument("--family_file", metavar="PATH", default="")
    parser.add_argument("--reference_species", metavar="SPECIES", required=True)
    parser.add_argument("--out_columns", metavar="PATH", required=True)
    parser.add_argument("--out_glyphs", metavar="PATH", required=True)
    parser.add_argument("--out_tree", metavar="PATH", required=True)
    parser.add_argument("--out_synteny", metavar="PATH", required=True)
    parser.add_argument("--out_ufboot", metavar="PATH", required=True)
    return parser


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


def resolve_query_cds_definitions(by_id, query_tip_by_id, definitions, cds_fasta_ids, family_id):
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
        definition["query_label"] = definition["gene_id"]
    if missing:
        raise ValueError(
            f"Query gene-tree tips were not found in the family CDS FASTA: family={family_id}, "
            f"missing={missing}"
        )
    return definition_by_id


def select_query_tip_nodes(rows, definitions):
    definition_ids = {definition["query_id"] for definition in definitions}
    candidates = defaultdict(list)
    for row in rows:
        if str(row.get("so_event") or "") != "L":
            continue
        branch_id = _int_value(row.get("branch_id"))
        for source_type, query_id in parse_query_marker_sources(row.get("query_marker_source")):
            if query_id in definition_ids:
                priority = 0 if source_type == "direct" else 1
                candidates[query_id].append((priority, branch_id))

    selected = {}
    for definition in definitions:
        query_id = definition["query_id"]
        options = candidates.get(query_id, [])
        if not options:
            continue
        best_priority = min(priority for priority, _ in options)
        nodes = sorted({node for priority, node in options if priority == best_priority})
        if len(nodes) != 1:
            raise ValueError(
                f"Query marker maps to multiple equally preferred gene-tree tips: query={query_id}, tips={nodes}"
            )
        selected[query_id] = nodes[0]
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


def collect_reference_ufboot_evidence(store, columns, glyphs):
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
                    if mrca_event != "S":
                        raise ValueError(
                            "Ortholog glyph pair does not have a speciation MRCA: "
                            f"family={family_id}, candidate={candidate_cds_fasta_id!r}, "
                            f"reference={reference_cds_fasta_id!r}, mrca={mrca}, "
                            f"event={mrca_event!r}"
                        )
                    if mrca == root:
                        row["orthology_ufboot_unavailable_reason"] = "mrca_is_root"
                    else:
                        support = ufboot_by_branch.get(mrca)
                        if support is not None:
                            row["decisive_branch_ufboot"] = support
                            row["orthology_ufboot_status"] = "evaluated"
                            row["orthology_ufboot_unavailable_reason"] = ""
                glyph_evidence_rows.append(row)
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


def collect_family_orthologs(
    rows,
    cds_fasta_ids,
    reference_species,
    family_id,
    family_order,
    first_column_order,
):
    by_id, children, root = build_tree_index(rows)
    query_tip_by_id, definitions = select_reference_tip_nodes(rows, reference_species)
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
    normalized_reference_species = normalize_species_label(reference_species)
    reference_glyphs = [
        glyph
        for glyph in glyphs
        if normalize_species_label(glyph["species"]) == normalized_reference_species
    ]
    observed_reference_cds_ids = []
    reference_glyph_is_one_to_one = len(reference_glyphs) == len(ordered_query_ids)
    for glyph in reference_glyphs:
        glyph_reference_ids = [
            value for value in str(glyph["reference_cds_fasta_ids"]).split(";") if value
        ]
        glyph_gene_ids = [value for value in str(glyph["gene_ids"]).split(";") if value]
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
    for output_rows in (columns, glyphs, tree_nodes):
        for output_row in output_rows:
            output_row["reference_species"] = normalized_reference_species
    return columns, glyphs, tree_nodes


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


def run(args):
    columns, glyphs, tree_nodes = collect_query_gene_orthologs(
        dir_gene_family=args.dir_gene_family,
        dir_query_gene=args.dir_query_gene,
        reference_species=args.reference_species,
        family_file=args.family_file,
    )
    store = GeneFamilyOutputStore(args.dir_gene_family)
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


def main():
    run(build_arg_parser().parse_args())


if __name__ == "__main__":
    main()
