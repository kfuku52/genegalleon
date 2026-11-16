"""Self-anchor checks for whether a copy split falls on a terminal branch."""

import itertools
import random
from collections import defaultdict, deque
from dataclasses import dataclass
from math import isfinite
from statistics import median

TERMINAL_PAIR_FIELDS = ("species", "species_event_id", "branch_id", "family_id", "block", "gene_a", "gene_b",
                        "contig_a", "contig_b", "tree_status", "sister_species", "sister_gene", "outgroup_species",
                        "outgroup_gene", "quartet_status", "copy_ds", "copy_pair_status",
                        "sister_a_status", "outgroup_a_status", "sister_b_status", "outgroup_b_status", "exterior_status",
                        "signed_terminal_margin",
                        "other_sums_difference")
TERMINAL_FIELDS = ("species", "species_event_id", "branch_id", "status", "same_family_anchors", "sampled_families",
                   "topology_compatible_families", "six_pair_families", "terminal_minimum_families",
                   "non_tied_families", "original_blocks", "median_copy_ds", "copy_lower", "copy_upper",
                   "median_signed_margin", "margin_lower", "margin_upper")


@dataclass(frozen=True)
class SelfAnchor:
    species: str
    family: str
    block: int
    genes: tuple
    contigs: tuple


def sample_anchors(species, positions, links, membership, maximum, seed):
    blocks = defaultdict(list)
    for (a, b), ids in links.items():
        family = membership.get(a)
        if a < b and family is not None and family == membership.get(b):
            block = min(ids)
            blocks[block].append(SelfAnchor(species, family, block, (a, b), (positions[a][0], positions[b][0])))
    rng = random.Random(seed)
    queues = []
    for _, anchors in sorted(blocks.items()):
        anchors.sort(key=lambda anchor: anchor.genes)
        rng.shuffle(anchors)
        queues.append(deque(anchors))
    rng.shuffle(queues)
    active = deque(queues)
    used, selected = set(), []
    while active and len(selected) < maximum:
        queue = active.popleft()
        while queue:
            anchor = queue.popleft()
            if anchor.family not in used:
                used.add(anchor.family)
                selected.append(anchor)
                break
        if queue:
            active.append(queue)
    return sum(map(len, blocks.values())), selected


class TerminalKs:
    def __init__(self, tree, members, key, sources, inputs, maximum, seed):
        from nwkit.clade_index import CladeIndex
        from nwkit.util import assign_branch_ids

        beds, self_blocks = inputs
        self.codes = {name: row["genetic_code"] for name, row in sources.items()}
        self.tips = {leaf.name: leaf for leaf in tree.leaves()}
        index, ids = CladeIndex(tree), assign_branch_ids(tree)
        self.identifiers = {sp: (index.clade_id_for_node(tip), ids[tip]) for sp, tip in self.tips.items()}
        self.available = self.codes.keys() & self_blocks.keys()
        self.inventory, self.sampled = {}, defaultdict(list)
        self.by_family, self.records, self.requested = defaultdict(list), {}, {}
        for sp in sorted(self.available):
            membership = {gene.strip(): row[key] for row in members if row[sp].strip() for gene in row[sp].split(",")}
            total, selected = sample_anchors(sp, beds[sp], self_blocks[sp][0], membership, maximum, seed)
            self.inventory[sp] = total
            if selected:
                self.sampled[sp] = selected
                for anchor in selected:
                    self.by_family[anchor.family].append(anchor)
        self.nearest = {}
        for sp in self.sampled:
            tip, code = self.tips[sp], self.codes[sp]
            sibling = set(tip.up.leaf_names()) - {sp}
            distances = {name: (tree.get_distance(tip, leaf), name)
                         for name, leaf in self.tips.items() if name != sp and self.codes.get(name) == code}
            self.nearest[sp] = ({n: d for n, d in distances.items() if n in sibling},
                                {n: d for n, d in distances.items() if n not in sibling})

    def observe(self, family, tree, species_by_gene):
        anchors = self.by_family.get(family)
        if not anchors:
            return
        leaves = {leaf.name: leaf for leaf in tree.leaves()}
        for anchor in anchors:
            self.records[anchor] = self._place(tree, leaves, species_by_gene, anchor)

    def _place(self, tree, leaves, species_by_gene, anchor):
        def select(names, allowed):
            options = [(*allowed[sp], gene) for gene in names if (sp := species_by_gene[gene]) in allowed]
            return min(options)[1:] if options else (None, None)

        a, b = anchor.genes
        if a not in leaves or b not in leaves:
            return "anchor_missing_from_tree", (None,) * 4, ()
        dup = tree.common_ancestor([leaves[a], leaves[b]])
        if len(dup.children) != 2:
            return "unresolved_copy_split", (None,) * 4, ()
        if any(species_by_gene[gene] != anchor.species for gene in dup.leaf_names()):
            return "copy_ancestor_contains_other_species", (None,) * 4, ()
        if dup.up is None:
            return "no_gene_sister", (None,) * 4, ()
        if len(dup.up.children) != 2:
            return "unresolved_gene_sister", (None,) * 4, ()
        nearby, exterior = self.nearest[anchor.species]
        sister, s = select((gene for child in dup.up.children if child is not dup for gene in child.leaf_names()), nearby)
        if s is None:
            return "no_same_code_nearest_sister", (None,) * 4, ()
        outside, o = select(leaves.keys() - set(dup.up.leaf_names()), exterior)
        if o is None:
            return "no_same_code_rooted_outgroup", (sister, s, None, None), ()
        pairs = tuple(tuple(sorted(pair)) for pair in itertools.combinations((a, b, s, o), 2))
        for pair in pairs:
            self.requested[pair] = (*map(species_by_gene.__getitem__, pair), anchor.family, "terminal_self_anchor")
        return "topology_compatible", (sister, s, outside, o), pairs

    def summarize(self, results):
        from nwkit.ksrate_model import median_order_statistic_interval

        def interval(values):
            if not values:
                return None, None, None
            lo, hi = median_order_statistic_interval(values, .05)
            return median(values), lo if isfinite(lo) else None, hi if isfinite(hi) else None

        pairs, summaries = [], []
        for sp, branch in sorted(self.identifiers.items()):
            margins, by_block = [], defaultdict(list)
            selected = self.sampled.get(sp, ())
            compatible = 0
            for anchor in selected:
                topology, neighbors, selected_pairs = self.records.get(anchor, ("family_tree_unavailable", (None,) * 4, ()))
                compatible += topology == "topology_compatible"
                status, copy, margin, other = "not_requested", None, None, None
                pair_states = ("not_requested",) * 6
                if selected_pairs:
                    values = [results[pair] for pair in selected_pairs]
                    pair_states = tuple(s for s, _ in values)
                    status = "dS_incomplete"
                    if values[0][0] == "ok":
                        copy = float(values[0][1])
                    if all(s == "ok" for s, _ in values):
                        copy, a_s, a_o, b_s, b_o, s_o = (float(value) for _, value in values)
                        u, v = a_s + b_o, a_o + b_s
                        margin = min(u, v) - (copy + s_o)
                        other = abs(u - v)
                        status = "terminal_minimum" if margin > 0 else "tie" if margin == 0 else "not_terminal_minimum"
                        by_block[anchor.block].append((copy, margin))
                        margins.append(margin)
                values = (sp, *branch, anchor.family, anchor.block, *anchor.genes, *anchor.contigs, topology,
                          *neighbors, status, copy, *pair_states, margin, other)
                pairs.append(dict(zip(TERMINAL_PAIR_FIELDS, values, strict=True)))
            status = ("unavailable_self_or_cds" if sp not in self.available else "no_same_family_anchor" if not selected else
                      "no_complete_quartets" if not margins else "measured")
            averages = tuple(tuple(map(median, zip(*records, strict=True))) for records in by_block.values())
            ages, gaps = zip(*averages, strict=True) if averages else ((), ())
            values = (sp, *branch, status, self.inventory.get(sp, 0), len(selected), compatible,
                      len(margins), sum(v > 0 for v in margins), sum(v != 0 for v in margins), len(by_block),
                      *interval(ages), *interval(gaps))
            summaries.append(dict(zip(TERMINAL_FIELDS, values, strict=True)))
        return pairs, summaries
