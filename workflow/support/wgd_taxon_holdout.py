"""Leave-one-taxon sensitivity for exact branch-matched duplication scores."""

from collections import Counter, defaultdict
from functools import lru_cache
from types import MappingProxyType

HOLDOUT_FIELDS = ("species_event_id", "branch_id", "excluded_taxon", "branch_identifiable", "eligible_families",
                  "supporting_families", "gained_votes", "lost_votes_among_eligible", "duplication_score",
                  "baseline_eligible_families", "baseline_supporting_families", "baseline_duplication_score")


def selected_bits(mask):
    while mask:
        bit = mask & -mask
        yield bit
        mask ^= bit


def lone_taxa(groups):
    result = 0
    for group in groups:
        if not group & (group - 1):
            result |= group
    return result


class DuplicationTaxonHoldout:
    def __init__(self, taxa, branches):
        # Match the scan's sorted taxon order and use its original sister and
        # daughter groups, not an adjacent branch after pruning a singleton.
        self.taxa = {1 << i: name for i, name in enumerate(taxa)}
        self.universe = sum(self.taxa)
        self.groups = {}
        self.by_event = {}
        self.projections = defaultdict(list)
        for mask, (_, sisters, daughters, event) in branches.items():
            groups = sisters, *daughters
            identifiable = self.universe ^ lone_taxa(groups)
            self.groups[mask] = groups
            self.by_event[event] = mask, identifiable
            for drop in selected_bits(identifiable):
                keep = self.universe ^ drop
                projected = mask & keep
                rule = (mask, drop, keep, tuple(group & keep for group in daughters))
                # Dropping a descendant can also expose an already smaller
                # gene-tree clade; dropping an external taxon that is absent
                # from the node leaves an original vote unchanged.
                if drop & mask:
                    self.projections[projected].append(rule)
                self.projections[projected | drop].append(rule)
        self.ineligible_after_exclusion = Counter()
        self.votes = Counter()
        self.gained = Counter()
        self.lost = Counter()
        self._eligibility = lru_cache(maxsize=512)(self._eligible_after_exclusion)

    def _eligible_after_exclusion(self, full):
        eligible, ineligible = {}, []
        for mask, groups in self.groups.items():
            present = tuple(full & group for group in groups)
            if all(present):
                lost_taxa = lone_taxa(present)
                eligible[mask] = self.universe ^ lost_taxa
                ineligible.extend((mask, taxon) for taxon in selected_bits(lost_taxa))
        return MappingProxyType(eligible), tuple(ineligible)

    def observe(self, tree, observed, baseline_votes):
        eligible, ineligible = self._eligibility(observed[tree])
        self.ineligible_after_exclusion.update(ineligible)
        # A removed child makes a unary node disappear and may make a
        # multifurcation binary. Other changes need no new gene-tree parse.
        votes = defaultdict(int, {mask: eligible.get(mask, 0) & ~mask for mask in baseline_votes})
        for node, union in observed.items():
            if node.is_leaf or union not in self.projections:
                continue
            children = tuple(observed[c] for c in node.children)
            for mask, drop, keep, daughters in self.projections[union]:
                if not (eligible.get(mask, 0) & drop) or votes[mask] & drop:
                    continue
                copies = [m & keep for m in children if m & keep]
                if len(copies) == 2:
                    left, right = copies
                    if left & right and all(left & group and right & group for group in daughters):
                        votes[mask] |= drop
        for mask, drops in votes.items():
            self.votes.update((mask, drop) for drop in selected_bits(drops))
            if mask not in baseline_votes:
                self.gained.update((mask, drop) for drop in selected_bits(drops))
        for mask in baseline_votes:
            lost = eligible.get(mask, 0) & ~votes.get(mask, 0)
            self.lost.update((mask, drop) for drop in selected_bits(lost))

    def rows(self, candidates):
        for row in candidates:
            event = row["species_event_id"]
            mask, identifiable = self.by_event[event]
            for bit, name in self.taxa.items():
                key = mask, bit
                total = row["num_eligible_families"] - self.ineligible_after_exclusion[key]
                support = self.votes[key]
                yield dict(zip(HOLDOUT_FIELDS, (event, row["branch_id"], name, int(bool(identifiable & bit)), total,
                                                 support, self.gained[key], self.lost[key], support / total if total else "NA",
                                                 row["num_eligible_families"], row["num_supporting_families"],
                                                 row["duplication_score"]), strict=True))
