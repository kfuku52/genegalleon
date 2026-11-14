"""Inspect additive dS contrasts for exact four-copy synteny paths."""

import random
from collections import defaultdict

SLOTS = ((0, 1), (2, 3), (0, 2), (1, 3), (0, 3), (1, 2))
QUARTET_FIELDS = ("species_event_id", "branch_id", "species_a", "species_b", "contig", "window", "family_id",
                  "gene_a", "gene_b", "gene_a2", "gene_b2", "pair_statuses", "six_ds", "paired_sum",
                  "within_species_sum", "crossed_sum", "paired_margin", "other_sums_difference", "status")


def quartet_requests(traced_hits, candidates, codes, limit, seed):
    eligible = {r["species_event_id"]: r["branch_id"] for r in candidates
                if r["event_support"] == "duplication-and-synteny-candidate"}
    units = defaultdict(dict)
    for event, a, b, window, family, genes, *_ in traced_hits:
        if (event not in eligible or len(genes) != 4 or a not in codes or b not in codes
                or codes[a] != codes[b]):
            continue
        unit = window, family
        current = units[event, a, b].get(unit)
        if current is None or genes < current:
            units[event, a, b][unit] = genes
    requests = {}
    selected = []
    rng = random.Random(seed)
    for (event, a, b), pool in sorted(units.items()):
        names = a, b, a, b
        keys = sorted(pool)
        if len(keys) > limit:
            keys = sorted(rng.sample(keys, limit))
        for (contig, window), family in keys:
            genes = pool[(contig, window), family]
            pairs = []
            for i, j in SLOTS:
                x, y = sorted((i, j), key=lambda k: (names[k], genes[k]))
                pair = genes[x], genes[y]
                pairs.append(pair)
                requests.setdefault(pair, (names[x], names[y], family, "quartet"))
            selected.append(((event, eligible[event], a, b, contig, window, family, *genes), pairs))
    return selected, requests


def quartet_results(selected, results):
    for fields, pairs in selected:
        metrics = [results[pair] for pair in pairs]
        statuses, distances = zip(*metrics, strict=True)
        expected = within = crossed = margin = difference = None
        status = "incomplete"
        if all(state == "ok" for state in statuses):
            expected, within, crossed = (sum(map(float, distances[i:i + 2])) for i in (0, 2, 4))
            margin = min(within, crossed) - expected
            difference = abs(within - crossed)
            status = "paired_minimum" if margin > 0 else "tied_minimum" if margin == 0 else "paired_not_minimum"
        values = (*fields, ",".join(statuses), ",".join(map(str, distances)), expected, within, crossed,
                  margin, difference, status)
        yield dict(zip(QUARTET_FIELDS, values, strict=True))
