import sqlite3
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas
import pytest

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


def test_attach_orthogroup_besthits_is_many_to_one_and_orders_columns(tmp_path):
    mod = load_module()
    frame = pandas.DataFrame(
        {
            "orthogroup": ["OG0001", "OG0001", "OG0002"],
            "site": [1, 2, 3],
        }
    )
    annotation_path = tmp_path / "Orthogroups.GeneCount.annotated.tsv"
    pandas.DataFrame(
        {
            "Orthogroup": ["OG0001", "OG0002"],
            "besthit_0.05": ["hit-a", "hit-b"],
            "besthit_0.25": ["hit-c", "hit-d"],
            "besthit_0.5": ["hit-e", "hit-f"],
            "besthit_0.75": ["hit-g", "hit-h"],
            "besthit_0.95": ["hit-i", "hit-j"],
            "unused_species_count": [4, 5],
        }
    ).to_csv(annotation_path, sep="\t", index=False)

    observed = mod.attach_orthogroup_besthits(frame, annotation_path)

    assert observed.shape[0] == frame.shape[0]
    assert observed.columns[:6].tolist() == ["orthogroup", *mod.ORTHOGROUP_BESTHIT_COLUMNS]
    assert observed.loc[observed["orthogroup"].eq("OG0001"), "besthit_0.05"].tolist() == [
        "hit-a",
        "hit-a",
    ]
    assert "unused_species_count" not in observed.columns


def test_attach_orthogroup_besthits_warns_on_partial_coverage_and_keeps_unmatched_na(
    tmp_path,
    capsys,
):
    mod = load_module()
    frame = pandas.DataFrame({"orthogroup": ["OG0001", "OG0002"], "site": [1, 2]})
    annotation_path = tmp_path / "Orthogroups.GeneCount.annotated.tsv"
    pandas.DataFrame(
        {
            "Orthogroup": ["OG0001"],
            **{column: [f"value-{column}"] for column in mod.ORTHOGROUP_BESTHIT_COLUMNS},
        }
    ).to_csv(annotation_path, sep="\t", index=False)

    observed = mod.attach_orthogroup_besthits(frame, annotation_path)

    assert pandas.isna(observed.loc[1, "besthit_0.5"])
    assert "matched 1 of 2 unique CSUBST orthogroups (50.0%)" in capsys.readouterr().err


def test_attach_orthogroup_besthits_rejects_duplicate_annotation_keys(tmp_path):
    mod = load_module()
    frame = pandas.DataFrame({"orthogroup": ["OG0001"]})
    annotation_path = tmp_path / "Orthogroups.GeneCount.annotated.tsv"
    pandas.DataFrame(
        {
            "Orthogroup": ["OG0001", " OG0001 "],
            **{column: ["hit-a", "hit-b"] for column in mod.ORTHOGROUP_BESTHIT_COLUMNS},
        }
    ).to_csv(annotation_path, sep="\t", index=False)

    with pytest.raises(ValueError, match="duplicate Orthogroup keys: OG0001"):
        mod.attach_orthogroup_besthits(frame, annotation_path)


def test_attach_orthogroup_besthits_rejects_zero_coverage(tmp_path):
    mod = load_module()
    frame = pandas.DataFrame({"orthogroup": ["OG0001"]})
    annotation_path = tmp_path / "Orthogroups.GeneCount.annotated.tsv"
    pandas.DataFrame(
        {
            "Orthogroup": ["HOG0001"],
            **{column: ["hit"] for column in mod.ORTHOGROUP_BESTHIT_COLUMNS},
        }
    ).to_csv(annotation_path, sep="\t", index=False)

    with pytest.raises(ValueError, match="matched 0 of 1 CSUBST orthogroups"):
        mod.attach_orthogroup_besthits(frame, annotation_path)


def test_attach_orthogroup_besthits_missing_optional_file_preserves_summary(tmp_path, capsys):
    mod = load_module()
    frame = pandas.DataFrame({"orthogroup": ["OG0001"], "site": [1]})

    observed = mod.attach_orthogroup_besthits(frame, tmp_path / "missing.tsv")

    pandas.testing.assert_frame_equal(observed, frame)
    assert "writing csubst summaries without besthit columns" in capsys.readouterr().err.lower()


def test_attach_orthogroup_besthits_without_annotation_argument_is_noop():
    mod = load_module()
    frame = pandas.DataFrame({"orthogroup": ["query-family-1"], "site": [1]})

    assert mod.attach_orthogroup_besthits(frame, None) is frame


def test_attach_orthogroup_besthits_requires_all_five_columns(tmp_path):
    mod = load_module()
    frame = pandas.DataFrame({"orthogroup": ["OG0001"]})
    annotation_path = tmp_path / "Orthogroups.GeneCount.annotated.tsv"
    pandas.DataFrame(
        {
            "Orthogroup": ["OG0001"],
            "besthit_0.05": ["hit"],
        }
    ).to_csv(annotation_path, sep="\t", index=False)

    with pytest.raises(ValueError, match="missing required columns: besthit_0.25"):
        mod.attach_orthogroup_besthits(frame, annotation_path)


def test_plot_paths_include_support_rate_and_pvalue_qvalue_distribution_pdfs(tmp_path):
    mod = load_module()
    paths = mod.plot_paths(tmp_path / "orthogroup_csubst_aa_change")

    assert "evidence_density" not in paths
    assert "foreground_unit_support_matrix" not in paths
    assert paths["support_significance_rate"].endswith(
        "orthogroup_csubst_aa_change_min_support_2_support_significance_rate.pdf"
    )
    assert paths["pvalue_qvalue_distributions"].endswith(
        "orthogroup_csubst_aa_change_min_support_2_pvalue_qvalue_distributions.pdf"
    )


def test_support_significance_data_retains_flat_max_t_as_zero_rate():
    mod = load_module()
    frame = pandas.DataFrame(
        {
            "support_fraction": [0.05, 0.15, 0.55, 0.95],
            "q_rate_enrichment_global": [0.01, 0.04, 0.2, 1.0],
            "q_rate_enrichment_empirical_global": [0.02, 0.08, 0.4, 1.0],
            "q_rate_enrichment_empirical_maxT_global": [1.0, 1.0, 1.0, 1.0],
        }
    )

    centers, candidate_counts, series = mod.support_significance_data(frame)

    assert centers.shape == (10,)
    assert candidate_counts.sum() == 4
    observed = {item["method"]["short_label"]: item for item in series}
    assert observed["Analytical"]["significant_counts"].sum() == 2
    assert observed["Empirical"]["significant_counts"].sum() == 1
    assert observed["Empirical maxT"]["significant_counts"].sum() == 0
    assert observed["Empirical maxT"]["percentages"][0] == 0.0


def test_probability_count_label_reports_three_threshold_counts():
    mod = load_module()
    values = pandas.Series([0.0005, 0.001, 0.01, 0.02, 0.05, 0.2]).to_numpy()

    observed = mod.probability_count_label(mod.PVALUE_QVALUE_METHODS[0], values)

    assert observed == "Analytical: 5 / 3 / 2"


def test_recalculate_sensitivity_qvalues_updates_global_and_grouped_fdr():
    mod = load_module()
    frame = pandas.DataFrame(
        {
            "orthogroup": ["OG1", "OG1", "OG2", "OG2"],
            "trait": ["A", "A", "A", "B"],
            "scan_match": ["m1", "m1", "m1", "m2"],
            "p_rate_enrichment": [0.01, 0.04, 0.03, 0.002],
            "p_rate_enrichment_empirical": [0.02, 0.08, 0.06, 0.004],
            "p_rate_enrichment_empirical_maxT": [0.1, 0.4, 0.3, 0.02],
            "q_rate_enrichment_global": [1.0] * 4,
            "q_rate_enrichment": [1.0] * 4,
        }
    )

    q_columns = mod.recalculate_sensitivity_qvalues(frame)

    assert "q_rate_enrichment_global" in q_columns
    assert "q_rate_enrichment_by_trait_match" in q_columns
    assert frame["q_rate_enrichment_global"].tolist() == [0.02, 0.04, 0.04, 0.008]
    assert frame["q_rate_enrichment"].tolist() == [0.02, 0.04, 0.03, 0.004]
    assert frame["q_rate_enrichment_by_trait_match"].tolist() == [0.02, 0.04, 0.03, 0.002]


def test_write_min_support_sensitivity_writes_threshold_series(tmp_path):
    mod = load_module()
    frame = pandas.DataFrame(
        {
            "orthogroup": ["OG1", "OG1", "OG2", "OG2"],
            "trait": ["A", "A", "A", "B"],
            "scan_match": ["m1", "m1", "m1", "m2"],
            "support_unit_count": [2, 3, 4, 5],
            "p_rate_enrichment": [0.001, 0.01, 0.03, 0.5],
            "p_rate_enrichment_empirical": [0.002, 0.02, 0.04, 0.8],
            "p_rate_enrichment_empirical_maxT": [0.05, 0.1, 0.2, 1.0],
            "q_rate_enrichment_global": [0.004, 0.03, 0.04, 0.5],
            "q_rate_enrichment_empirical_global": [0.008, 0.04, 0.06, 0.8],
            "q_rate_enrichment_empirical_maxT_global": [0.2, 0.3, 0.4, 1.0],
            "besthit_0.05": ["hit-a", "hit-b", "hit-c", "hit-d"],
            "besthit_0.25": ["hit-e", "hit-f", "hit-g", "hit-h"],
            "besthit_0.5": ["hit-i", "hit-j", "hit-k", "hit-l"],
            "besthit_0.75": ["hit-m", "hit-n", "hit-o", "hit-p"],
            "besthit_0.95": ["hit-q", "hit-r", "hit-s", "hit-t"],
        }
    )
    out_prefix = tmp_path / "orthogroup_csubst_aa_change"
    stale_paths = mod.min_support_sensitivity_paths(out_prefix, 6)
    stale_paths["output_dir"].mkdir(parents=True, exist_ok=True)
    stale_paths["summary_tsv"].write_text("stale\n", encoding="utf-8")
    stale_paths["plot_pdf"].write_text("stale\n", encoding="utf-8")

    manifest_path = mod.write_min_support_sensitivity(frame, out_prefix)

    assert manifest_path == tmp_path / "orthogroup_csubst_aa_change_min_support_manifest.tsv"
    manifest = pandas.read_csv(manifest_path, sep="\t")
    assert manifest["min_support"].tolist() == [3, 4, 5]
    assert manifest["candidate_rows"].tolist() == [3, 2, 1]
    assert manifest["q_rate_enrichment_global_le_0.05"].tolist() == [2, 0, 0]
    assert not stale_paths["summary_tsv"].exists()
    assert not stale_paths["plot_pdf"].exists()
    assert not (tmp_path / "min_support_sensitivity").exists()
    for threshold, expected_rows in ((3, 3), (4, 2), (5, 1)):
        paths = mod.min_support_sensitivity_paths(out_prefix, threshold)
        subset = pandas.read_csv(paths["summary_tsv"], sep="\t")
        assert subset.shape[0] == expected_rows
        assert (subset["support_unit_count"] >= threshold).all()
        assert set(mod.ORTHOGROUP_BESTHIT_COLUMNS).issubset(subset.columns)
        assert paths["plot_pdf"].is_file()
        assert paths["plot_pdf"].stat().st_size > 1000


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

    annotation_path = tmp_path / "Orthogroups.GeneCount.annotated.tsv"
    pandas.DataFrame(
        {
            "Orthogroup": ["OG0001", "OG0002"],
            **{
                column: [f"{column}-og1", f"{column}-og2"]
                for column in mod.ORTHOGROUP_BESTHIT_COLUMNS
            },
        }
    ).to_csv(annotation_path, sep="\t", index=False)

    out_prefix = tmp_path / "orthogroup_csubst_aa_change"
    out_tsv = tmp_path / "orthogroup_csubst_aa_change_min_support_2_summary.tsv"
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
            "--orthogroup_annotation_tsv",
            str(annotation_path),
        ],
    )

    assert mod.main() == 0
    assert out_tsv.is_file()
    primary = pandas.read_csv(out_tsv, sep="\t")
    assert primary.columns[:6].tolist() == ["orthogroup", *mod.ORTHOGROUP_BESTHIT_COLUMNS]
    assert primary.loc[primary["orthogroup"].eq("OG0002"), "besthit_0.5"].eq(
        "besthit_0.5-og2"
    ).all()
    for path in mod.plot_paths(out_prefix).values():
        assert Path(path).is_file()
    sensitivity_paths = mod.min_support_sensitivity_paths(out_prefix, 3)
    assert sensitivity_paths["manifest"].is_file()
    assert sensitivity_paths["summary_tsv"].is_file()
    assert sensitivity_paths["plot_pdf"].is_file()
    sensitivity = pandas.read_csv(sensitivity_paths["summary_tsv"], sep="\t")
    assert sensitivity.columns[:6].tolist() == ["orthogroup", *mod.ORTHOGROUP_BESTHIT_COLUMNS]
    assert not (tmp_path / "min_support_sensitivity").exists()


def test_remove_legacy_min_support_output_layout_removes_only_generated_files(tmp_path):
    mod = load_module()
    out_prefix = tmp_path / "orthogroup_csubst_aa_change"
    legacy_primary = tmp_path / "orthogroup_csubst_aa_change_summary.tsv"
    legacy_primary.write_text("legacy\n", encoding="utf-8")
    legacy_dir = tmp_path / "min_support_sensitivity"
    legacy_dir.mkdir()
    legacy_generated = legacy_dir / "orthogroup_csubst_aa_change_min_support_3_summary.tsv"
    legacy_generated.write_text("legacy\n", encoding="utf-8")
    unrelated = legacy_dir / "keep.txt"
    unrelated.write_text("keep\n", encoding="utf-8")

    mod.remove_legacy_min_support_output_layout(out_prefix)

    assert not legacy_primary.exists()
    assert not legacy_generated.exists()
    assert unrelated.read_text(encoding="utf-8") == "keep\n"
