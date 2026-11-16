"""Match genome-stage syntenic paths to exact daughter copies in a family tree."""

import csv

try:
    from wgd_cross_phase import trace_has_unique_phase
except ImportError:
    from .wgd_cross_phase import trace_has_unique_phase

CROSS_PATH_FIELDS = ("species_event_id", "family_id", "species_a", "species_b", "contig", "window", "gene_a", "gene_b",
                     "gene_a2", "gene_b2", "unique_phase", "both_orientations", "used_in_candidate_mode")
NODE_PATH_FIELDS = ("family_id", "gene_clade_id", "species_event_id", "species_a", "species_b", "contig", "window",
                    "gene_a", "gene_b", "gene_a2", "gene_b2", "unique_phase", "both_orientations",
                    "used_in_candidate_mode", "branch_event_support", "mapping_status")


def cross_paths(hits, phases, shared, *, require_unique=False):
    seen = set()
    for hit in sorted(hits):
        event, a, b, (contig, window), family, genes, *_ = hit
        if len(genes) != 4:
            continue
        key = hit[:6]
        if key in seen:
            continue
        seen.add(key)
        unique_flag = int(trace_has_unique_phase(phases, hit))
        shared_flag = int(hit in shared) if shared is not None else None
        selected = shared_flag if shared is not None else unique_flag if require_unique else 1
        yield dict(zip(CROSS_PATH_FIELDS, (event, family, a, b, contig, window, *genes, unique_flag, shared_flag,
                                          selected), strict=True))


def load_family_paths(path, family):
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t", strict=True)
        if tuple(reader.fieldnames or ()) != CROSS_PATH_FIELDS:
            raise ValueError("Invalid cross-synteny node path header")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Invalid cross-synteny node path row width")
            genes = tuple(row[key] for key in CROSS_PATH_FIELDS[6:10])
            if len(set(genes)) != 4 or any(not gene for gene in genes):
                raise ValueError("Cross-synteny node path must identify four distinct genes")
            if (row["unique_phase"] not in {"0", "1"} or row["both_orientations"] not in {"", "0", "1"}
                    or row["used_in_candidate_mode"] not in {"0", "1"}):
                raise ValueError("Invalid cross-synteny node path flags")
            if row["family_id"] == family:
                yield row


def daughter_below(node, ancestor):
    while node.up is not ancestor:
        node = node.up
        if node is None or node is ancestor:
            return None
    return node


def cross_node_rows(tree, clades, leaf_species, events, reconciliation, paths):
    leaves = {leaf.name: leaf for leaf in tree.leaves()}
    duplication = {row["gene_clade_id"]: row for row in reconciliation if row["event_type"] == "duplication"}
    for path in paths:
        event = path["species_event_id"]
        a, b = path["species_a"], path["species_b"]
        genes = tuple(path[key] for key in CROSS_PATH_FIELDS[6:10])
        clade = None
        status = "genes_not_in_family_tree"
        if all(gene in leaves for gene in genes):
            if any(leaf_species.get(gene) != name for gene, name in zip(genes, (a, b, a, b), strict=True)):
                status = "species_mapping_mismatch"
            else:
                x, y, xx, yy = (leaves[gene] for gene in genes)
                node = tree.common_ancestor([x, xx])
                clade = clades.clade_id_for_node(node)
                if (len(node.children) != 2 or daughter_below(x, node) is not daughter_below(y, node)
                        or daughter_below(xx, node) is not daughter_below(yy, node)):
                    status = "copy_split_mismatch"
                elif (record := duplication.get(clade)) is None or record["mapping_status"] != "mapped" or record["event_status"] != "resolved":
                    status = "duplication_not_resolved"
                elif record["species_event_id"] != event:
                    status = "species_branch_mismatch"
                else:
                    status = "exact_node_match"
                    if path["used_in_candidate_mode"] == "1":
                        node.add_prop("cross_synteny_matched_paths", node.props.get("cross_synteny_matched_paths", 0) + 1)
        yield {**{k: path[k] for k in CROSS_PATH_FIELDS}, "gene_clade_id": clade,
               "branch_event_support": events.get(event, {}).get("event_support", "unresolved"), "mapping_status": status}
