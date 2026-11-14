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


def unique_phase_support(quartets, floor):
    windows, families = set(), set()
    for (window, family, *_), phase in quartets.items():
        if phase == 1:
            windows.add(window)
            families.add(family)
            if len(windows) >= floor and len(families) >= floor:
                return True
    return False


def trace_has_unique_phase(phases, hit):
    event, a, b, window, family, genes, *_ = hit
    if len(genes) != 4:
        return False
    first, second, first_copy, second_copy = genes
    key = window, family, first, first_copy, *sorted((second, second_copy))
    return phases.get((event, a, b), {}).get(key) == 1


ORIENTATION_PHASE_FIELDS = ("species_event_id", "branch_id", "species_a", "species_b", "mode", "status",
                            "a_quartets", "b_quartets", "shared_quartets", "shared_families",
                            "shared_windows_a", "shared_windows_b")


def orientation_phase_overlap(candidates, first_rows, second_rows, forward, reverse, *, unique=False):
    branches = {row["species_event_id"]: row["branch_id"] for row in candidates}
    reverse_status = {(r["species_event_id"], r["species_b"], r["species_a"]): r["status"] for r in second_rows}
    forward_status = {(r["species_event_id"], r["species_a"], r["species_b"]): r["status"] for r in first_rows}

    def observed(scanner, invert, unique_phase):
        result = defaultdict(dict)
        for hit in scanner.orientation_hits if invert else scanner.traced_hits:
            event, a, b, window, family, genes, *_ = hit
            if len(genes) != 4 or (unique_phase and not trace_has_unique_phase(scanner.phases, hit)):
                continue
            first, other_a, second, other_b = genes
            copies = ((other_a, first), (other_b, second)) if invert else ((first, other_a), (second, other_b))
            pair = (event, b, a) if invert else (event, a, b)
            key = family, *sorted(copies)
            result[pair].setdefault(key, set()).add(window)
        return result

    for want_unique in ((False, True) if unique else (False,)):
        first, second = observed(forward, False, want_unique), observed(reverse, True, want_unique)
        for (event, a, b), status in forward_status.items():
            pair = event, a, b
            x, y = first.get(pair, {}), second.get(pair, {})
            assessed_a, assessed_b = status == "assessed", reverse_status[pair] == "assessed"
            available = assessed_a and assessed_b
            intersection = x.keys() & y.keys() if available else ()
            values = (event, branches[event], a, b, "unique_phase" if want_unique else "all_concordant",
                      "assessed_both" if available else "not_assessed_both",
                      len(x) if assessed_a else None, len(y) if assessed_b else None,
                      len(intersection) if available else None,
                      len({key[0] for key in intersection}) if available else None,
                      len({w for key in intersection for w in x[key]}) if available else None,
                      len({w for key in intersection for w in y[key]}) if available else None)
            yield dict(zip(ORIENTATION_PHASE_FIELDS, values, strict=True))
