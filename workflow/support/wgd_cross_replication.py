"""Compare exact focal copy pairs across separately supplied anchor comparisons."""

from collections import defaultdict
from itertools import combinations

REPLICATION_FIELDS = ("species_event_id", "branch_id", "focal", "comparator_a", "comparator_b", "phase_filter",
                      "copy_pairs_a", "copy_pairs_b", "shared_copy_pairs", "shared_families")


def replicated_copy_rows(candidates, comparisons, traced_hits, unique_hits=None):
    branch_ids = {row["species_event_id"]: row["branch_id"] for row in candidates}
    assessable = defaultdict(set)
    for row in comparisons:
        if row["status"] != "assessed":
            continue
        event, a, b = (row[key] for key in ("species_event_id", "species_a", "species_b"))
        assessable[event, a].add(b)
        assessable[event, b].add(a)
    replicated = {key: names for key, names in assessable.items() if len(names) >= 2}
    if not replicated:
        return
    for label, paths in (("all_paired", traced_hits), ("unique_phase", unique_hits)):
        if paths is None:
            continue
        by_comparison = defaultdict(set)
        for event, a, b, _, family, genes, *_ in paths:
            if len(genes) != 4:
                continue
            for i, (focal, other) in enumerate(((a, b), (b, a))):
                if other in replicated.get((event, focal), ()):
                    pair = tuple(sorted((genes[i], genes[i + 2])))
                    by_comparison[event, focal, other].add((family, pair))
        for (event, focal), comparators in sorted(replicated.items()):
            for a, b in combinations(sorted(comparators), 2):
                first, second = (by_comparison.get((event, focal, other), set()) for other in (a, b))
                shared = first & second
                yield dict(zip(REPLICATION_FIELDS, (event, branch_ids[event], focal, a, b, label,
                               len(first), len(second), len(shared), len({family for family, _ in shared})), strict=True))
