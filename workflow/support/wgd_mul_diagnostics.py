"""Join native MUL assignments to origin evidence without changing its rules."""

import json
import sys
from collections import defaultdict

DIAGNOSTIC_FIELDS = (
    "family_id", "gene_node", "gene_clade_id", "gene_topology_id", "species_topology_id",
    "reconciliation_event_type", "species_event_id", "classification", "reason",
    "mul_mapping_status", "mul_dl_duplication", "mul_best_hypotheses", "mul_optimal_mappings",
    "mul_candidate_ids", "mul_num_states", "mul_states", "meaning",
)


def join_nodes(tree, species, rows, model, reconciliation, origins, family_id):
    from nwkit.clade_index import CladeIndex
    from nwkit.mul_reconcile_nodes import NODE_FIELDS, rooted_topology_id
    from nwkit.util import read_tree

    if (model.get("method") != "exact-MUL-LCA-DL-parsimony-v1"
            or model.get("num_gene_trees") != 1
            or model.get("node_diagnostics", {}).get("schema") != "nwkit-mul-node-assignments-v1"):
        raise ValueError("Unsupported MUL node diagnostic model")
    best = model.get("best_hypotheses", [])
    if not best or len(best) != len(set(best)) or any(type(c) is not int or c < 0 for c in best):
        raise ValueError("Invalid tied-best MUL candidates")
    score_rows = {r["mul.tree"]: r for r in model["scores"]}
    if len(score_rows) != len(model["scores"]) or not set(best).issubset(score_rows):
        raise ValueError("Incomplete MUL candidate scores")
    if (any(type(r["score"]) is not int or r["score"] < 0 for r in score_rows.values())
            or set(best) != {c for c, r in score_rows.items() if r["score"] == min(x["score"] for x in score_rows.values())}):
        raise ValueError("MUL tied-best candidates do not match minimum scores")
    gene_id, species_id = rooted_topology_id(tree), rooted_topology_id(species)
    clades = CladeIndex(tree)
    nodes = list(tree.traverse("postorder"))
    by_clade = {clades.clade_id_for_node(n): n for n in nodes}
    events = {r["gene_clade_id"]: r for r in reconciliation}
    origin_by_clade = {r["gene_clade_id"]: r for r in origins}
    if len(events) != len(reconciliation) or set(events) != set(by_clade):
        raise ValueError("Reconciliation nodes do not match the diagnostic gene tree")
    duplication_ids = {k for k, r in events.items() if r["event_type"] == "duplication"}
    if len(origin_by_clade) != len(origins) or set(origin_by_clade) != duplication_ids:
        raise ValueError("Origin classifications do not match reconciliation duplications")
    allowed_species = set(species.leaf_names())
    candidate_nodes = {
        c: list(read_tree(score_rows[c]["labeled.tree"], "auto", True, quiet=True).traverse("preorder"))
        for c in best
    }
    groups, limits, states = defaultdict(dict), {}, defaultdict(list)
    for row in rows:
        if not set(NODE_FIELDS).issubset(row):
            raise ValueError("Incomplete MUL node diagnostic columns")
        if (row["gene_topology_id"] != gene_id or row["species_topology_id"] != species_id
                or row["gene.tree"] != "1"):
            raise ValueError("MUL diagnostic tree identity mismatch")
        candidate, mapping = int(row["mul.tree"]), int(row["mapping.id"])
        clade = row["gene_clade_id"]
        if candidate not in best or mapping < 1 or clade not in by_clade:
            raise ValueError("Unknown MUL candidate, mapping or gene node")
        score = score_rows[candidate]
        if any(row[k] != str(score[k]) for k in ("h1.node", "h2.node", "hypothesis.kind")):
            raise ValueError("MUL candidate metadata mismatch")
        node = by_clade[clade]
        index, leaf, count = int(row["gene_node"]), int(row["is_leaf"]), int(row["optimal.mappings"])
        if (not 0 <= index < len(nodes) or leaf != int(node.is_leaf)
                or row["gene_label"] != (node.name or "") or count < 1):
            raise ValueError("MUL gene node metadata mismatch")
        if limits.setdefault(candidate, count) != count:
            raise ValueError("Inconsistent optimal mapping count")
        key = candidate, mapping
        if clade in groups[key]:
            raise ValueError("Duplicate MUL gene node assignment")
        mul_node = int(row["mul_node"])
        if not 0 <= mul_node < len(candidate_nodes[candidate]):
            raise ValueError("Unknown mapped MUL node")
        target = candidate_nodes[candidate][mul_node]
        tips = json.loads(row["mul_descendant_tips"])
        mapped_species = json.loads(row["mapped_species"])
        expected_tips = sorted(target.leaf_names())
        expected_species = sorted({t if t in allowed_species else t[:-1] for t in expected_tips})
        if (tips != expected_tips or mapped_species != expected_species
                or not set(expected_species).issubset(allowed_species) or row["mul_label"] != target.name):
            raise ValueError("MUL mapped clade metadata mismatch")
        if node.is_leaf and (events[clade]["mapping_status"] != "mapped"
                             or mapped_species != [events[clade]["species_name"]]):
            raise ValueError("MUL and reconciliation leaf species differ")
        duplication, edge_losses, root_losses = (int(row[k]) for k in (
            "duplication", "child_edge_losses", "root_losses"))
        if (duplication not in (0, 1) or min(edge_losses, root_losses) < 0
                or (node.is_leaf and (duplication or edge_losses))
                or (node is not tree and root_losses) or int(row["total.score"]) != score["score"]):
            raise ValueError("Invalid MUL assignment score")
        state = (tuple(tips), tuple(mapped_species), duplication, edge_losses, root_losses)
        groups[key][clade] = (index, duplication + edge_losses + root_losses, mul_node)
        states[clade].append(state)
    for candidate in best:
        count = limits.get(candidate, 0)
        if not count or {k[1] for k in groups if k[0] == candidate} != set(range(1, count + 1)):
            raise ValueError("Incomplete optimal MUL mappings")
    for (candidate, _), assignments in groups.items():
        if (set(assignments) != set(by_clade)
                or {p[0] for p in assignments.values()} != set(range(len(nodes)))
                or sum(p[1] for p in assignments.values()) != score_rows[candidate]["score"]):
            raise ValueError("Incomplete or inconsistent MUL node assignment")
    for candidate in best:
        signatures = [tuple(v[2] for _, v in sorted(assignments.items()))
                      for (c, _), assignments in groups.items() if c == candidate]
        if len(signatures) != len(set(signatures)):
            raise ValueError("Repeated optimal MUL assignment")
    output = []
    for index, node in enumerate(nodes):
        if node.is_leaf:
            continue
        clade = clades.clade_id_for_node(node)
        options = sorted(set(states[clade]))
        duplication_options = {s[2] for s in options}
        status = "consistent" if len(options) == 1 else "ambiguous"
        dl_duplication = "all" if duplication_options == {1} else "none" if duplication_options == {0} else "ambiguous"
        origin = origin_by_clade.get(clade, {})
        node.add_prop("mul_mapping_status", status)
        node.add_prop("mul_dl_duplication", dl_duplication)
        node.add_prop("mul_gene_node", index)
        node.add_prop("mul_best_hypotheses", len(best))
        node.add_prop("mul_optimal_mappings", len(groups))
        output.append({
            "family_id": family_id, "gene_node": index, "gene_clade_id": clade,
            "gene_topology_id": gene_id, "species_topology_id": species_id,
            "reconciliation_event_type": events[clade]["event_type"],
            "species_event_id": events[clade]["species_event_id"],
            "classification": origin.get("classification", "NA"), "reason": origin.get("reason", "NA"),
            "mul_mapping_status": status, "mul_dl_duplication": dl_duplication,
            "mul_best_hypotheses": len(best), "mul_optimal_mappings": len(groups),
            "mul_candidate_ids": json.dumps(best), "mul_num_states": len(options),
            "mul_states": json.dumps([{"mul_descendant_tips": s[0], "mapped_species": s[1],
                                        "duplication": s[2], "child_edge_losses": s[3], "root_losses": s[4]}
                                       for s in options], sort_keys=True),
            "meaning": "cooptimal_DL_assignments_not_origin_probabilities",
        })
    return output


def run_diagnostics(args, tree, species, reconciliation, origins, run_tool, read_table, write_table):
    output = args.output / "mul"
    output.mkdir()
    command = [sys.executable, "-m", "nwkit", "mul-reconcile", "--infile", args.gene_tree,
               "--species-tree", args.species_tree, "--species-parser", args.species_parser,
               "--h1", args.mul_h1, "--outfile", output / "scores.tsv", "--model-out", output / "results.json",
               "--node-out", output / "nodes.tsv", "--max-candidates", args.mul_max_candidates,
               "--max-state-pairs", args.mul_max_state_pairs, "--max-maps", args.mul_max_maps]
    if args.mul_h2:
        command += ["--h2", args.mul_h2]
    if args.species_regex:
        command += ["--species-regex", args.species_regex]
    if args.species_map:
        command += ["--species-map-tsv", args.species_map]
    run_tool(command, args.output, "mul-reconcile")
    model = json.loads((output / "results.json").read_text())
    diagnostics = join_nodes(tree, species, read_table(output / "nodes.tsv"), model, reconciliation, origins, args.family_id)
    write_table(args.output / "node_diagnostics.tsv", diagnostics, DIAGNOSTIC_FIELDS)
    plot_diagnostics(tree, args.output / "duplication_origins_mul.pdf")
    return {"enabled": True, "score_model": "dl", "h1": args.mul_h1, "h2": args.mul_h2,
            "best_hypotheses": model["best_hypotheses"],
            "num_consistent_nodes": sum(r["mul_mapping_status"] == "consistent" for r in diagnostics),
            "num_ambiguous_nodes": sum(r["mul_mapping_status"] == "ambiguous" for r in diagnostics),
            "origin_rule_changed": False, "meaning": "cooptimal assignments, not probabilities; not locus-MC node inference"}


def plot_diagnostics(tree, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    colors = {"WGD-supported": "#197c53", "SSD-supported": "#be4d31", "unresolved": "#777777"}
    leaves = list(tree.leaves())
    y = {n: i for i, n in enumerate(reversed(leaves))}
    x = {tree: 0}
    for n in tree.traverse("preorder"):
        for child in n.children:
            x[child] = x[n] + 1
    depth = max(1, max(x.values()))
    fig, ax = plt.subplots(figsize=(10, max(4, len(leaves) * 0.3 + 2)))
    tip_labels, node_labels = [], []
    for n in tree.traverse("postorder"):
        if n.children:
            y[n] = sum(y[c] for c in n.children) / len(n.children)
            ax.plot([x[n], x[n]], [min(y[c] for c in n.children), max(y[c] for c in n.children)],
                    color="#aaaaaa", linewidth=0.8)
        color = colors.get(n.props.get("duplication_origin"), "#333333")
        if n is not tree:
            ax.plot([x[n.up], x[n]], [y[n], y[n]], color=color, linewidth=1.4)
        if n.is_leaf:
            tip_labels.append(ax.annotate(n.name, (x[n], y[n]), xytext=(7, 0),
                                          textcoords="offset points", va="center", fontsize=10))
        else:
            ax.plot(x[n], y[n], "^" if n.props["mul_mapping_status"] == "ambiguous" else "o",
                    color=color, markersize=7)
            node_labels.append(ax.annotate(f"N{n.props['mul_gene_node']}", (x[n], y[n]), xytext=(5, 8),
                                           textcoords="offset points", fontsize=9))
    ax.set_ylim(-1, len(leaves))
    ax.axis("off")
    fig.legend(handles=[Line2D([], [], color=c, marker="o", label=s) for s, c in colors.items()]
               + [Line2D([], [], color="#333333", marker="o", label="Other events")],
               loc="upper center", bbox_to_anchor=(0.5, 0.99), ncol=4, frameon=False, fontsize=10)
    fig.legend(handles=[Line2D([], [], color="#333333", marker=m, linestyle="none", label=s)
                        for m, s in (("o", "MUL mapping consistent"), ("^", "MUL mapping ambiguous"))],
               loc="upper center", bbox_to_anchor=(0.5, 0.91), ncol=2, frameon=False, fontsize=10)
    fig.text(0.5, 0.035, "Colors: origin evidence. Shapes: tied-best D+L mappings. Neither is an origin posterior.",
             ha="center", fontsize=10)
    fig.subplots_adjust(top=0.78, bottom=0.14, left=0.05, right=0.95)
    # Measure actual glyphs, retaining space per depth level on comb-like trees.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    tip_width = max(t.get_window_extent(renderer).width for t in tip_labels) / fig.dpi
    node_width = max((t.get_window_extent(renderer).width for t in node_labels), default=0) / fig.dpi
    unit = max(0.45, node_width + 0.15)
    right_margin = (tip_width + 0.3) / unit
    fig.set_size_inches(max(10, (1.05 * depth * unit + tip_width + 0.3) / 0.9), fig.get_figheight())
    ax.set_xlim(-0.05 * depth, depth + right_margin)
    fig.savefig(path)
    plt.close(fig)
