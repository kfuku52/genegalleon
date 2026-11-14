#!/usr/bin/env python3
"""Plan sparse cross-species synteny comparisons on a rooted species tree."""

import argparse
import hashlib
import heapq
import itertools
from collections import defaultdict
from pathlib import Path

try:
    from busco_reference_quality import safe_token
    from wgd_cross_synteny import read_pairs
    from wgd_evidence import read_table
    from wgd_ssd import protect_output, species_tree, write_table
except ImportError:
    from .busco_reference_quality import safe_token
    from .wgd_cross_synteny import read_pairs
    from .wgd_evidence import read_table
    from .wgd_ssd import protect_output, species_tree, write_table

PAIR_FIELDS = ("analysis_id", "target_species", "query_species")
AUDIT_FIELDS = ("species_event_id", "branch_id", "descendant_taxa", "num_daughters", "num_daughters_with_self_synteny",
                "num_daughters_with_existing_pairs", "num_daughters_with_assessed_pairs",
                "num_daughters_after_planned_pairs", "num_planned_pairs", "status")


def closest_untried(first, second, used):
    # Enumerate shortest sums lazily, stopping before unneeded species pairs.
    if len(first) > len(second):
        first, second = second, first

    def candidate(i, j):
        da, a = first[i]
        db, b = second[j]
        return da + db, tuple(sorted((a, b))), i, j

    queue = [candidate(i, 0) for i in range(len(first))]
    heapq.heapify(queue)
    while queue:
        distance, pair, i, j = heapq.heappop(queue)
        if pair not in used:
            return distance, pair
        if j + 1 < len(second):
            heapq.heappush(queue, candidate(i, j + 1))
    return None


def plan_pairs(tree, pairs, assessed=None):
    from nwkit.clade_index import CladeIndex, LcaIndex
    from nwkit.util import assign_branch_ids

    clades, lcas, ids = CladeIndex(tree), LcaIndex(tree), assign_branch_ids(tree)
    self_species = {a for a, b in pairs if a == b}
    leaves = {leaf.name: leaf for leaf in tree.leaves()}
    depth = {tree: 0.0}
    for node in tree.descendants("preorder"):
        depth[node] = depth[node.up] + node.dist
    previous = defaultdict(set)
    for a, b in pairs:
        if a != b:
            previous[lcas.common_ancestor(leaves[a], leaves[b])].add(tuple(sorted((a, b))))
    rows, audit = [], []
    for node in tree.traverse("postorder"):
        if node.is_root or node.is_leaf:
            continue
        groups = [[(depth[leaves[name]] - depth[node], name) for name in clades.names_for_node(child)
                   if name in self_species] for child in node.children]
        groups = sorted((sorted(group) for group in groups if group), key=lambda group: group[0][1])
        group_for = {name: i for i, group in enumerate(groups) for _, name in group}
        used = previous[node]
        eligible = {pair for pair in used if all(name in group_for for name in pair)}
        known = {group_for[name] for pair in eligible for name in pair}
        observed = known if assessed is None else {group_for[name] for pair in eligible & assessed for name in pair}
        covered = set(observed)
        count = 0
        options = {}
        for i, j in itertools.combinations(range(len(groups)), 2):
            if i not in covered or j not in covered:
                choice = closest_untried(groups[i], groups[j], used)
                if choice:
                    options[i, j] = choice
        while options and len(covered) < len(groups):
            viable = ((-(len({i, j} - covered)), distance, pair, i, j)
                      for (i, j), (distance, pair) in options.items() if not {i, j} <= covered)
            best = min(viable, default=None)
            if best is None:
                break
            _, _, (a, b), i, j = best
            del options[i, j]
            covered.update((i, j))
            count += 1
            key = hashlib.sha256(f"{a}\0{b}".encode()).hexdigest()[:20]
            rows.append(dict(zip(PAIR_FIELDS, (f"wgd-{key}", a, b), strict=True)))
        all_groups = len(node.children)
        if len(groups) < all_groups:
            status = "missing_self_synteny"
        elif len(covered) < all_groups:
            status = "unassessed_pairs_exhausted"
        elif count:
            status = "comparisons_planned"
        else:
            status = "covered_by_assessed_pairs" if assessed is not None else "covered_by_existing_pairs"
        audit.append(dict(zip(AUDIT_FIELDS, (clades.clade_id_for_node(node), ids[node], clades.csv_for_node(node),
                                             all_groups, len(groups), len(known), len(observed) if assessed is not None else "NA",
                                             len(covered), count, status), strict=True)))
    return rows, audit


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-tree", required=True, type=Path)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--cross-synteny", required=True, type=Path, help="Existing self/cross anchor manifest")
    parser.add_argument("--pair-evidence", type=Path, help="Optional cross_synteny_pairs.tsv for retrying unassessed branches")
    parser.add_argument("--pairs-out", required=True, type=Path, help="New cross-species pairs in pairwise synteny input format")
    parser.add_argument("--audit-out", required=True, type=Path)
    args = parser.parse_args(argv)
    workspace = args.workspace.resolve()
    tree = species_tree(args.species_tree)
    names = {safe_token(name, "species") for name in tree.leaf_names()}
    inputs = {args.species_tree.resolve(), args.cross_synteny.resolve()}
    manifest = read_pairs(args.cross_synteny, lambda path: (workspace / path).resolve(), names, inputs.add)
    for path in inputs:
        if not path.is_file():
            raise ValueError(f"Missing synteny plan input: {path}")
    pairs = {tuple(sorted((row["species_a"], row["species_b"]))) for row in manifest}
    assessed = None
    if args.pair_evidence:
        inputs.add(args.pair_evidence.resolve())
        assessed, seen = set(), set()
        for row in read_table(args.pair_evidence):
            if not {"species_a", "species_b", "status"} <= row.keys():
                raise ValueError("Pair evidence requires species_a, species_b and status")
            pair = tuple(sorted((row["species_a"], row["species_b"])))
            if pair[0] == pair[1] or pair not in pairs or pair in seen:
                raise ValueError(f"Duplicate or unknown cross pair in pair evidence: {pair}")
            seen.add(pair)
            if row["status"] == "assessed":
                assessed.add(pair)
            elif row["status"] not in {"insufficient_cross_anchors", "insufficient_archive_anchors", "missing_self_comparison"}:
                raise ValueError(f"Unrecognized pair evidence status: {row['status']}")
    for output, other in ((args.pairs_out, args.audit_out), (args.audit_out, args.pairs_out)):
        protect_output(output, inputs | {other.resolve()}, workspace)
        if output.exists():
            raise ValueError(f"Refusing to replace existing planning output: {output}")
    rows, audit = plan_pairs(tree, pairs, assessed)
    write_table(args.pairs_out, rows, PAIR_FIELDS)
    write_table(args.audit_out, audit, AUDIT_FIELDS)


if __name__ == "__main__":
    main()
