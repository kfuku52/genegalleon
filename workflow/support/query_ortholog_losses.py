"""Link explicit GeneRax terminal losses to displayed ortholog anchors."""

from bisect import bisect_left
from collections import defaultdict


def terminal_loss_links(species_taxa, records, anchors, plotted_tips, species_key):
    """Return (species, anchor) pairs supported by a GeneRax speciation.

    An explicit terminal loss must branch off at the target's species split,
    with a surviving lineage on the other side. Link the loss to anchors in
    that sister lineage or across a deeper, compatible speciation. Duplications
    do not establish orthology between an anchor and a lost lineage.
    """
    # The XML clades are nested. The last containing clade seen for any tip
    # is the immediate parent when visiting from largest to smallest.
    parent_by_species, containing_clade = {}, {}
    for clade in sorted(species_taxa, key=lambda clade: -len(species_taxa[clade])):
        tips = species_taxa[clade]
        parent_by_species[clade] = containing_clade.get(tips[0])
        containing_clade.update(dict.fromkeys(tips, clade))

    def side(branch, split):
        while branch is not None and branch != split:
            parent = parent_by_species[branch]
            if parent == split:
                return branch
            branch = parent
        return None

    # NWKit emits preorder records. Subtrees are contiguous, so a loss only
    # examines anchors in the surviving sister and in opposite speciation
    # subtrees, even when the family has many unrelated paralogs.
    children = defaultdict(list)
    leaves = {}
    subtree_end = list(range(1, len(records) + 1))
    has_leaf = [False] * len(records)
    previous_split = [None] * len(records)
    for node, row in enumerate(records):
        parent = row["parent"]
        if parent is not None:
            children[parent].append(node)
            previous_split[node] = (
                (parent, node) if records[parent]["events"][-1]["type"] == "speciation"
                else previous_split[parent]
            )
        for event in row["events"]:
            if event["type"] == "leaf":
                leaves[row["label"]] = (node, event["species"])
                has_leaf[node] = True
    for node in reversed(range(len(records))):
        parent = records[node]["parent"]
        if parent is not None:
            subtree_end[parent] = max(subtree_end[parent], subtree_end[node])
            has_leaf[parent] |= has_leaf[node]

    mismatches = [
        anchor for anchor, species in anchors.items()
        if anchor not in leaves or species_key(leaves[anchor][1]) != species
    ]
    if mismatches:
        raise ValueError(f"Displayed anchors missing from GeneRax XML or assigned to other species: {mismatches}")
    positions = sorted(leaves[anchor][0] for anchor in anchors)
    links = set()

    def link_subtree(subtree, split, lost_side, target):
        start = bisect_left(positions, subtree)
        stop = bisect_left(positions, subtree_end[subtree], lo=start)
        for index in range(start, stop):
            anchor = records[positions[index]]["label"]
            anchor_side = side(leaves[anchor][1], split)
            if anchor_side is not None and anchor_side != lost_side:
                links.add((target, anchor))

    for node, row in enumerate(records):
        events = row["events"]
        if len(events) != 1 or events[0]["type"] != "loss" or has_leaf[node]:
            continue
        target = events[0]["species"]
        if target not in plotted_tips:
            continue
        owner, split = row["parent"], parent_by_species[target]
        if owner is None or split is None or len(children[owner]) != 2:
            continue
        event = records[owner]["events"][-1]
        if event["type"] != "speciation" or event["species"] != split:
            continue
        left, right = children[owner]
        sister = right if left == node else left
        sister_side = side(records[sister]["events"][0]["species"], split)
        if not has_leaf[sister] or sister_side is None or sister_side == target:
            continue
        link_subtree(sister, split, target, target)
        ancestor = previous_split[owner]
        while ancestor is not None:
            parent, child = ancestor
            split = records[parent]["events"][-1]["species"]
            lost_side = side(target, split)
            if lost_side is not None:
                for other in children[parent]:
                    if other != child:
                        link_subtree(other, split, lost_side, target)
            ancestor = previous_split[parent]
    return links
