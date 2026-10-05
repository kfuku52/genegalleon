"""Support filters and BH correction for the resulting candidate family."""

import numpy as np
import pandas as pd

ANALYTICAL_P_COLUMN = "p_rate_enrichment_asymptotic"
SUPPORT_Q_COLUMN = "q_rate_enrichment_asymptotic_support_filtered"
SUPPORT_BH_POLICY = "bh_after_unit_and_lineage_support_filter_v1"
SUPPORT_BH_SCOPE = "all_orthogroups_traits_matches_after_support_filters"
SUPPORT_BH_COLUMNS = [
    "support_bh_candidate_count", "support_bh_test_count", "support_bh_undefined_count",
    "support_bh_min_unit_support", "support_bh_min_lineage_support",
]


def calculate_bh_fdr(pvalues):
    """BH over finite analytical tests; undefined tests remain undefined."""
    values = pd.to_numeric(pd.Series(pvalues), errors="raise").to_numpy(dtype=float, na_value=np.nan)
    invalid = np.isinf(values) | (np.isfinite(values) & ((values < 0) | (values > 1)))
    if invalid.any():
        raise ValueError("Analytical P values must be finite probabilities in [0, 1] or missing.")
    qvalues = np.full(values.shape, np.nan)
    positions = np.flatnonzero(np.isfinite(values))
    if positions.size:
        order = np.argsort(values[positions], kind="stable")
        ranked_positions = positions[order]
        adjusted = values[ranked_positions] * positions.size / np.arange(1, positions.size + 1)
        qvalues[ranked_positions] = np.minimum(1.0, np.minimum.accumulate(adjusted[::-1])[::-1])
    return qvalues


def rate_testable_mask(values):
    # SQLite integer flags become floats when a combined table contains NULL.
    testable = values.astype(str).str.strip().str.lower().eq("true") | pd.to_numeric(values, errors="coerce").eq(1)
    return testable.fillna(False)


def validated_analytical_pvalues(frame, source):
    if ANALYTICAL_P_COLUMN not in frame:
        raise ValueError(f"{source}: missing {ANALYTICAL_P_COLUMN}; BH requires analytical P values.")
    raw = frame[ANALYTICAL_P_COLUMN]
    values = pd.to_numeric(raw, errors="coerce")
    invalid = (raw.notna() & values.isna()) | (values.notna() & (~np.isfinite(values) | ~values.between(0, 1)))
    if invalid.any():
        raise ValueError(f"{source}: invalid probabilities in {ANALYTICAL_P_COLUMN}.")
    if "scan_rate_testable" in frame:
        testable = rate_testable_mask(frame["scan_rate_testable"])
        if (values.notna() & ~testable).any():
            raise ValueError(f"{source}: an untestable candidate has a finite probability.")
    return values


def support_filtered_bh(frame, min_unit_support, min_lineage_support, source):
    """Filter both bounds, then correct all finite P before any report selection."""
    keep = support_mask(frame, min_unit_support, min_lineage_support, source)
    pvalues = validated_analytical_pvalues(frame, source)
    selected = frame.loc[keep].copy()
    selected[SUPPORT_Q_COLUMN] = calculate_bh_fdr(pvalues.loc[keep])
    metadata = {
        "probability_policy": SUPPORT_BH_POLICY,
        "inference_scope": SUPPORT_BH_SCOPE,
        "support_bh_candidate_count": len(selected),
        "support_bh_test_count": int(selected[SUPPORT_Q_COLUMN].notna().sum()),
        "support_bh_undefined_count": int(selected[SUPPORT_Q_COLUMN].isna().sum()),
        "support_bh_min_unit_support": int(min_unit_support),
        "support_bh_min_lineage_support": int(min_lineage_support),
    }
    for column in SUPPORT_BH_COLUMNS:
        selected[column] = metadata[column]
    selected.attrs["support_bh_metadata"] = metadata
    return selected


def validate_minimum(value, name):
    if not isinstance(value, (int, np.integer)) or value < 0:
        raise ValueError(f"{name} must be an integer >= 0 (0 disables the condition).")


def validated_counts(frame, column, source):
    if column not in frame:
        raise ValueError(f"{source}: missing {column}; rebuild summaries from a scan containing this column.")
    counts = pd.to_numeric(frame[column], errors="coerce")
    values = counts.to_numpy(dtype=float, na_value=np.nan)
    invalid = ~np.isfinite(values) | (values < 0) | (values != np.trunc(values))
    if invalid.any():
        examples = []
        for position in np.flatnonzero(invalid)[:5]:
            row = frame.iloc[position]
            identity = ", ".join(f"{key}={row[key]}" for key in ("orthogroup", "trait") if key in row)
            examples.append(f"row {position + 2}" + (f" ({identity})" if identity else ""))
        raise ValueError(
            f"{source}: {column} must contain finite nonnegative integers for every candidate; "
            f"{int(invalid.sum())} invalid/missing row(s): {'; '.join(examples)}."
        )
    return counts


def support_mask(frame, min_unit_support, min_lineage_support, source):
    validate_minimum(min_unit_support, "min_unit_support")
    validate_minimum(min_lineage_support, "min_lineage_support")
    keep = pd.Series(True, index=frame.index)
    if min_unit_support:
        keep &= validated_counts(frame, "support_unit_count", source) >= min_unit_support
    if min_lineage_support:
        keep &= validated_counts(frame, "support_lineage_count", source) >= min_lineage_support
    return keep


def support_view_prefix(out_prefix, min_unit_support, min_lineage_support):
    return f"{out_prefix}_min_unit_support_{min_unit_support}_min_lineage_support_{min_lineage_support}"
