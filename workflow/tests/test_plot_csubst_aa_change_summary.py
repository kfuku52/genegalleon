import sqlite3
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

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
                "future_csubst_metric": 42.0,
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
        "future_csubst_metric",
    }
    assert expected_columns.issubset(observed.columns)
    assert observed.loc[0, "scan_rate_exposure"] == "state_aware"
    assert observed.loc[0, "site_rate"] == 0.125
    assert observed.loc[0, "future_csubst_metric"] == 42.0
    assert observed.columns[-1] == "future_csubst_metric"

    ranked, score_column, score_kind = mod.ranked_candidates(observed)
    assert ranked.shape[0] == 1
    assert score_column == "q_rate_enrichment_global"
    assert score_kind == "FDR"


def test_plot_paths_include_pvalue_qvalue_distribution_pdf(tmp_path):
    mod = load_module()
    paths = mod.plot_paths(tmp_path / "orthogroup_csubst_aa_change")

    assert paths["pvalue_qvalue_distributions"].endswith(
        "orthogroup_csubst_aa_change_pvalue_qvalue_distributions.pdf"
    )


def test_probability_count_label_reports_three_threshold_counts():
    mod = load_module()
    values = pandas.Series([0.0005, 0.001, 0.01, 0.02, 0.05, 0.2]).to_numpy()

    observed = mod.probability_count_label(mod.PVALUE_QVALUE_METHODS[0], values)

    assert observed == "Analytical: 5 / 3 / 2"


def test_write_pvalue_qvalue_distributions(tmp_path):
    mod = load_module()
    frame = pandas.DataFrame(
        {
            "p_rate_enrichment": [0.001, 0.01, 0.2, 1.0],
            "p_rate_enrichment_empirical": [0.002, 0.02, 0.3, 1.0],
            "p_rate_enrichment_empirical_maxT": [0.05, 0.4, 0.9, 1.0],
            "q_rate_enrichment_global": [0.004, 0.03, 0.4, 1.0],
            "q_rate_enrichment_empirical_global": [0.006, 0.04, 0.5, 1.0],
            "q_rate_enrichment_empirical_maxT_global": [1.0, 1.0, 1.0, 1.0],
        }
    )
    out_pdf = tmp_path / "pvalue_qvalue_distributions.pdf"

    mod.write_pvalue_qvalue_distributions(frame, out_pdf)

    assert out_pdf.is_file()
    assert out_pdf.stat().st_size > 1000


def test_main_writes_pvalue_qvalue_distribution_by_default(tmp_path, monkeypatch):
    mod = load_module()
    db_path = tmp_path / "scan.db"
    source = pandas.DataFrame(
        {
            "orthogroup": ["OG0001", "OG0001", "OG0002", "OG0002"],
            "trait": ["traitA"] * 4,
            "state_change": ["10V", "20L", "30A", "40G"],
            "from_state": ["A", "I", "V", "S"],
            "to_state": ["V", "L", "A", "G"],
            "support_fraction": [0.75, 0.5, 0.4, 0.25],
            "support_unit_count": [3, 2, 2, 1],
            "support_unit_ids": ["1,2,3", "1,2", "2,3", "3"],
            "p_rate_enrichment": [0.001, 0.01, 0.2, 1.0],
            "p_rate_enrichment_empirical": [0.002, 0.02, 0.3, 1.0],
            "p_rate_enrichment_empirical_maxT": [0.05, 0.4, 0.9, 1.0],
            "q_rate_enrichment_global": [0.004, 0.03, 0.4, 1.0],
            "q_rate_enrichment_empirical_global": [0.006, 0.04, 0.5, 1.0],
            "q_rate_enrichment_empirical_maxT_global": [1.0, 1.0, 1.0, 1.0],
        }
    )
    with sqlite3.connect(db_path) as conn:
        source.to_sql("aa_change", conn, index=False)

    out_prefix = tmp_path / "orthogroup_csubst_aa_change"
    out_tsv = tmp_path / "orthogroup_csubst_aa_change_summary.tsv"
    monkeypatch.setattr(
        "sys.argv",
        [
            str(SCRIPT_PATH),
            "--dbpath",
            str(db_path),
            "--out_prefix",
            str(out_prefix),
            "--out_tsv",
            str(out_tsv),
        ],
    )

    assert mod.main() == 0
    assert out_tsv.is_file()
    for path in mod.plot_paths(out_prefix).values():
        assert Path(path).is_file()
