"""Matched-family sensitivity for dS boundaries with common external species."""

import random
from collections import defaultdict
from dataclasses import dataclass
from statistics import median

import numpy as np

OUTGROUP_FIELDS = ("species_event_id", "branch_id", "focal", "younger_sister", "older_sister", "common_outgroup",
                   "num_eligible_families", "num_sampled_families", "num_measured_all_outgroups", "younger_median_ds",
                   "older_median_ds", "median_signed_gap", "gap_lower", "gap_upper", "status")


@dataclass(frozen=True)
class CommonOutgroups:
    event: str
    focal: str
    younger: str
    older: str
    outgroups: tuple
    eligible: int
    families: tuple

    def pair_keys(self, outside):
        return tuple(tuple(sorted(pair)) for pair in (
            (self.focal, self.younger), (self.focal, self.older), (self.younger, outside),
            (self.older, outside), (self.focal, outside),
        ))


def select_common_outgroups(focals, trios, nodes, index, archived, codes, maximum, seed):
    by_taxon = defaultdict(set)
    for family, present in archived.items():
        for sp in present:
            by_taxon[sp].add(family)
    by_trio = {(trio.focal, trio.species_event_id): trio for trio in trios}
    for event, focal in sorted(set(focals)):
        parent = nodes[event].up
        young = by_trio.get((focal, event))
        old = by_trio.get((focal, index.clade_id_for_node(parent)))
        result = CommonOutgroups(event, focal, young.sister if young else None, old.sister if old else None, (), 0, ())
        if young is None or old is None:
            yield result
            continue
        candidates = set.intersection(*(by_taxon[sp] for sp in (focal, young.sister, old.sister)))
        seen = set(parent.leaf_names())
        options = []
        for depth, ancestor in enumerate(parent.ancestors()):
            descendants = set(ancestor.leaf_names())
            for name in sorted(descendants - seen):
                if codes.get(name) == codes[focal] and (families := candidates & by_taxon[name]):
                    options.append((depth, name, families))
            seen = descendants
        if not options:
            yield result
            continue
        _, first, common = min(options, key=lambda item: (item[0], -len(item[2]), item[1]))
        second = min((item for item in options if item[1] != first and item[2] & common),
                     key=lambda item: (item[0], -len(item[2] & common), item[1]), default=None)
        names = (first,)
        if second:
            names += (second[1],)
            common &= second[2]
        families = sorted(common)
        random.Random(seed).shuffle(families)
        yield CommonOutgroups(event, focal, young.sister, old.sister, names, len(common), tuple(families[:maximum]))


def common_outgroup_requests(selections, archived):
    for selection in selections:
        for outside in selection.outgroups:
            for a, b in selection.pair_keys(outside):
                for family in selection.families:
                    genes = archived[family]
                    yield (genes[a], genes[b]), (a, b, family, "common_outgroup")


def common_outgroup_rows(selections, sampled, branch_ids):
    from nwkit.ksrate_model import median_order_statistic_interval

    denominator = sum(len(row.outgroups) for row in selections)
    for row in selections:
        if not row.outgroups:
            values = (row.event, branch_ids[row.event], row.focal, row.younger, row.older, None, 0, 0, 0,
                      None, None, None, None, None, "no_common_outgroup" if row.younger and row.older else "missing_trio")
            yield dict(zip(OUTGROUP_FIELDS, values, strict=True))
            continue
        by_outgroup = {}
        for outside in row.outgroups:
            distances = [sampled.get(key, {}) for key in row.pair_keys(outside)]
            by_outgroup[outside] = {f: tuple(d[f] for d in distances) for f in row.families if all(f in d for d in distances)}
        complete = sorted(set.intersection(*(set(families) for families in by_outgroup.values())))
        for outside, distances in by_outgroup.items():
            young = old = gap = lo = hi = None
            status = "no_complete_families"
            if complete:
                fy, fo, yg, og, fg = np.array([distances[f] for f in complete], dtype=float).T
                signed = (fo - fy) + (yg - og)  # The focal-outgroup term cancels exactly.
                young = median((np.maximum(fy, fg) - yg) + np.minimum(fy, fg))
                old = median((np.maximum(fo, fg) - og) + np.minimum(fo, fg))
                gap = median(signed)
                low, high = median_order_statistic_interval(signed, .05 / denominator)
                status = "interval_unavailable"
                if np.isfinite(low):
                    lo, hi = low, high
                    status = "ordered" if lo > 0 else "reversed" if hi < 0 else "order_unresolved"
            values = (row.event, branch_ids[row.event], row.focal, row.younger, row.older, outside,
                      row.eligible, len(row.families), len(complete), young, old, gap, lo, hi, status)
            yield dict(zip(OUTGROUP_FIELDS, values, strict=True))
