"""Shared association filtering and Holm/BH numeric kernels; callers own test scopes."""

import numpy
import pandas


def association_rows(frame: pandas.DataFrame) -> pandas.DataFrame:
    """Exclude intercept-only rows from association status and ranking."""
    if frame.empty:
        return frame
    mask = pandas.Series(True, index=frame.index, dtype=bool)
    intercept_labels = {"intercept", "(intercept)"}
    if "predictor_type" in frame:
        mask &= ~frame["predictor_type"].astype(str).str.strip().str.lower().isin(intercept_labels)
    for column in ("term", "source_term"):
        if column in frame:
            mask &= ~frame[column].astype(str).str.strip().str.lower().isin(intercept_labels)
    return frame.loc[mask]


def adjust_association_p_values(
    values: pandas.Series | numpy.ndarray,
) -> tuple[numpy.ndarray, numpy.ndarray]:
    """Return Holm family-wise and Benjamini-Hochberg adjusted p-values."""
    raw = numpy.asarray(values, dtype=float)
    count = len(raw)
    order = numpy.argsort(raw, kind="stable")
    ranked = raw[order]
    holm_ranked = numpy.minimum(
        1.0,
        numpy.maximum.accumulate((count - numpy.arange(count)) * ranked),
    )
    bh_ranked = numpy.minimum(
        1.0,
        numpy.minimum.accumulate((count / numpy.arange(count, 0, -1)) * ranked[::-1])[::-1],
    )
    holm = numpy.empty(count, dtype=float)
    bh = numpy.empty(count, dtype=float)
    holm[order] = holm_ranked
    bh[order] = bh_ranked
    return holm, bh
