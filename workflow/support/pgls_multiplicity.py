"""Expression association multiplicity, grouped by response and source predictor.

The comparison bundle is the single input: it already includes RSC rows.
Never combine it with the native RSC bundle (which would count tests twice).
"""

import numpy as np
import pandas as pd
from species_tree_pgls import _adjust_association_p_values, _association_rows, _usable_association_rows

PAIR_SCOPE = "all_families_methods_aggregations_by_response_predictor"
REQUIRED = {"tree_id", "analysis_method", "aggregation", "analysis_id", "response", "source_term", "term", "term_test", "p_value"}


def adjust_associations(frame):
    """Keep failed association rows; assign them p=1 internally for global BH.

    Counts refer to represented association rows, not unrun or missing models.
    Family-local adjustments retain the previous usable-row/method scope.
    """
    missing = REQUIRED.difference(frame.columns)
    if missing:
        raise ValueError(f"PGLS comparison is missing columns: {', '.join(sorted(missing))}")
    out = _association_rows(frame).copy().reset_index(drop=True)
    identity = ["tree_id", "analysis_method", "aggregation", "analysis_id", "response", "source_term", "term", "term_test"]
    identity += [c for c in ("model_id", "response_level", "predictor_level") if c in out]
    if out.duplicated(identity).any():
        raise ValueError("Duplicate PGLS association identity in comparison inputs")
    for column in ("tree_id", "analysis_method", "aggregation", "analysis_id", "response", "source_term", "term", "term_test"):
        if out[column].isna().any() or out[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"Missing PGLS association identity: {column}")
    raw = pd.to_numeric(out["p_value"], errors="coerce")
    usable = out.index.isin(_usable_association_rows(out).index)
    valid = usable & np.isfinite(raw) & raw.between(0, 1)
    out["p_value"] = raw
    out["p_value_usable"] = valid
    for column in ("p_value_family_holm", "p_value_family_bh", "p_value_global_bh"):
        out[column] = np.nan
    out["family_multiplicity_scope"] = "all_usable_family_associations_for_method"
    out["global_multiplicity_scope"] = PAIR_SCOPE
    out["global_n_associations"] = 0
    out["global_n_usable"] = 0
    for _, group in out.loc[valid].groupby(["tree_id", "analysis_method"], sort=False):
        holm, bh = _adjust_association_p_values(group["p_value"])
        out.loc[group.index, "p_value_family_holm"] = holm
        out.loc[group.index, "p_value_family_bh"] = bh
    for _, group in out.groupby(["response", "source_term"], sort=False):
        idx = group.index
        values = raw.loc[idx].where(valid.loc[idx], 1.0)
        _, bh = _adjust_association_p_values(values)
        out.loc[idx, "p_value_global_bh"] = np.where(valid.loc[idx], bh, np.nan)
        out.loc[idx, "global_n_associations"] = len(idx)
        out.loc[idx, "global_n_usable"] = int(valid.loc[idx].sum())
    return out


def write_association_table(engine, store):
    """Read live/ZIP comparison bundles and replace the DB's long-form table."""
    frames = []
    if store is not None:
        for name in sorted(store.file_names("pgls_comparison")):
            if not name.endswith(".tsv"):
                continue
            with store.open_binary("pgls_comparison", name) as handle:
                frame = pd.read_csv(handle, sep="\t", low_memory=False)
            if not frame.empty:
                frames.append(frame)
    frame = pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame(columns=sorted(REQUIRED))
    adjusted = adjust_associations(frame)
    with engine.begin() as conn:
        adjusted.to_sql("pgls_association", conn, if_exists="replace", index=False)
        conn.exec_driver_sql('CREATE INDEX pgls_association_pair ON pgls_association(response, source_term)')
        conn.exec_driver_sql('CREATE INDEX pgls_association_family ON pgls_association(tree_id)')
