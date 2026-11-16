"""Count how much the duplication score reuses synteny-matched tree families."""

from collections import defaultdict

OVERLAP_FIELDS = ("species_event_id", "branch_id", "eligible_families", "supporting_families",
                  "matched_paired_families", "matched_partial_only_families", "matched_either",
                  "eligible_after_exclusion", "support_after_exclusion", "score_after_exclusion")


def duplication_synteny_overlap(candidates, windows):
    paired, partial = defaultdict(set), defaultdict(set)
    for row in windows:
        event = row["species_event_id"]
        for groups, key in ((paired, "matched_paired_families"), (partial, "matched_partial_only_families")):
            if families := row[key]:
                groups[event].update(families.split(","))
    for row in candidates:
        event = row["species_event_id"]
        a, b = paired[event], partial[event]
        matched = len(a | b)
        eligible, support = row["num_eligible_families"], row["num_supporting_families"]
        if not 0 <= matched <= support <= eligible:
            raise ValueError(f"Cross-synteny matches must be a subset of voting tree families: {event}")
        rest, residual = eligible - matched, support - matched
        values = (event, row["branch_id"], eligible, support, len(a), len(b - a), matched,
                  rest, residual, residual / rest if rest else None)
        yield dict(zip(OVERLAP_FIELDS, values, strict=True))
