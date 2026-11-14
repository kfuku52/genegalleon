"""Audit whether direct cross blocks identify a unique four-copy phase."""

from collections import defaultdict

PHASE_STATES = ("concordant_only", "discordant_only", "ambiguous")
PHASE_FIELDS = ("species_event_id", "branch_id", "species_a", "species_b",
                *(f"{state}_{unit}" for state in PHASE_STATES for unit in ("quartets", "families", "windows")))


def observe_phase(source, other, genes):
    x, y, xx, yy = genes
    if xx not in other:
        return 0
    if y in source and yy in other:
        return 1
    if yy in source and y in other:
        return 2
    return 0


def record_phase(phases, species_pair, window, family, genes, phase):
    """Merge alternative paths using the same four genes, not just the family."""
    first, second, first_copy, second_copy = genes
    key = window, family, first, first_copy, *sorted((second, second_copy))
    by_quartet = phases[species_pair]
    by_quartet[key] = by_quartet.get(key, 0) | phase


def phase_rows(candidates, comparisons, phases):
    branches = {row["species_event_id"]: row["branch_id"] for row in candidates}
    for row in comparisons:
        event, a, b = (row[field] for field in ("species_event_id", "species_a", "species_b"))
        families, windows = defaultdict(set), defaultdict(set)
        counts = defaultdict(int)
        for key, mask in phases.get((event, a, b), {}).items():
            window, family = key[:2]
            state = PHASE_STATES[mask - 1]
            counts[state] += 1
            families[state].add(family)
            windows[state].add(window)
        yield dict(zip(PHASE_FIELDS, (event, branches[event], a, b,
                     *(value for state in PHASE_STATES for value in (counts[state], len(families[state]), len(windows[state])))),
                       strict=True))


def unique_phase_windows(rows):
    by_branch = defaultdict(int)
    for row in rows:
        by_branch[row["species_event_id"]] += row["concordant_only_windows"]
    return by_branch
