from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sqlite3

import pandas


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "plot_csubst_aa_change_summary.py"


def load_module():
    spec = spec_from_file_location("plot_csubst_aa_change_summary", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_read_table_retains_current_csubst_scan_rate_and_empirical_q_columns(tmp_path):
    mod = load_module()
    db_path = tmp_path / "scan.db"
    source = pandas.DataFrame(
        [
            {
                "orthogroup": "OG0001",
                "trait": "traitA",
                "state_change": "A10V",
                "scan_rate_exposure": "state_aware",
                "site_rate": 0.125,
                "site_rate_categorized": 3.0,
                "site_rate_quantile": 0.75,
                "p_rate_enrichment": 0.01,
                "q_rate_enrichment_global": 0.02,
                "q_rate_enrichment_empirical": 0.03,
                "q_rate_enrichment_empirical_by_trait": 0.04,
                "q_rate_enrichment_empirical_by_trait_match": 0.05,
            }
        ]
    )
    with sqlite3.connect(db_path) as conn:
        source.to_sql("aa_change", conn, index=False)
        observed = mod.read_table(conn, "aa_change")

    expected_columns = {
        "scan_rate_exposure",
        "site_rate",
        "site_rate_categorized",
        "site_rate_quantile",
        "q_rate_enrichment_empirical",
        "q_rate_enrichment_empirical_by_trait",
        "q_rate_enrichment_empirical_by_trait_match",
    }
    assert expected_columns.issubset(observed.columns)
    assert observed.loc[0, "scan_rate_exposure"] == "state_aware"
    assert observed.loc[0, "site_rate"] == 0.125

    ranked, score_column, score_kind = mod.ranked_candidates(observed)
    assert ranked.shape[0] == 1
    assert score_column == "q_rate_enrichment_global"
    assert score_kind == "FDR"
