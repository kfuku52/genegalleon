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


def support_view_source():
    return pandas.DataFrame({
        "orthogroup": ["OG1", "OG2", "OG3", "OG4"],
        "trait": ["aquatic"] * 4,
        "support_unit_count": [1, 5, 6, 7],
        "support_lineage_count": [1, 4, 3, 4],
        "support_fraction": [0.1, 0.5, 0.6, 0.7],
        "score_rate_enrichment": [4.0, 3.0, 2.0, 1.0],
        "p_rate_enrichment_asymptotic": [0.001, 0.002, 0.003, 0.004],
        "q_rate_enrichment_asymptotic_global": [0.01, 0.02, 0.03, 0.04],
    })


@pytest.mark.parametrize("units,lineages,expected", [
    (6, 4, ["OG4"]), (5, 4, ["OG2", "OG4"]), (6, 0, ["OG3", "OG4"]),
    (0, 4, ["OG2", "OG4"]), (0, 0, ["OG1", "OG2", "OG3", "OG4"]),
    (8, 4, []), (0, 5, []),
])
def test_summary_support_views_apply_bh_after_both_bounds_and_preserve_source_columns(tmp_path, units, lineages, expected):
    mod = load_module()
    source = support_view_source()
    prefix = tmp_path / "scan"
    manifest_path = mod.write_support_views(source, prefix, units, lineages)
    manifest = pandas.read_csv(manifest_path, sep="\t")
    primary_path = tmp_path / manifest.loc[0, "summary_tsv"]
    primary = pandas.read_csv(primary_path, sep="\t")
    assert primary["orthogroup"].tolist() == expected
    expected_frame = source.loc[source["orthogroup"].isin(expected)].reset_index(drop=True)
    pandas.testing.assert_frame_equal(primary.loc[:, source.columns], expected_frame, check_dtype=False)
    assert primary[mod.SUPPORT_Q_COLUMN].tolist() == pytest.approx([0.004] * len(expected))
    assert manifest.loc[0, "support_bh_test_count"] == len(expected)
    assert manifest["min_lineage_support"].eq(lineages).all()
    assert manifest["min_unit_support"].iloc[0] == units
    if units == 0:
        assert manifest["min_unit_support"].tolist() == [0]
    assert manifest["probability_policy"].eq("bh_after_unit_and_lineage_support_filter_v1").all()
    for row in manifest.itertuples():
        assert (tmp_path / row.summary_tsv).is_file()
        assert (tmp_path / row.plot_pdf).is_file()


@pytest.mark.parametrize("bad_count", [None, float("inf"), -1, 1.5, "unknown", "missing_column"])
def test_main_rejects_missing_lineage_counts_before_replacing_outputs(tmp_path, monkeypatch, bad_count):
    mod = load_module()
    source = support_view_source()
    if bad_count == "missing_column":
        source = source.drop(columns="support_lineage_count")
    else:
        source["support_lineage_count"] = source["support_lineage_count"].astype(object)
        source.loc[0, "support_lineage_count"] = bad_count
    db = tmp_path / "scan.db"
    with sqlite3.connect(db) as conn:
        source.to_sql("aa_change", conn, index=False)
    prefix = tmp_path / "scan"
    outputs = [tmp_path / "scan_min_support_2_summary.tsv", tmp_path / "scan_summary.tsv",
               tmp_path / "scan_all_candidates_summary.tsv", tmp_path / "scan_min_support_manifest.tsv"]
    for path in outputs:
        path.write_text("preserve\n")
    monkeypatch.setattr("sys.argv", [str(SCRIPT_PATH), "--dbpath", str(db), "--out_prefix", str(prefix),
                                    "--out_tsv", str(outputs[0]), "--min_unit_support", "6",
                                    "--min_lineage_support", "4"])
    with pytest.raises(ValueError, match="support_lineage_count"):
        mod.main()
    assert all(path.read_text() == "preserve\n" for path in outputs)


def test_custom_summary_preserves_all_candidates_and_legacy_views(tmp_path, monkeypatch):
    mod = load_module()
    source = support_view_source()
    db = tmp_path / "scan.db"
    with sqlite3.connect(db) as conn:
        source.to_sql("aa_change", conn, index=False)
    prefix = tmp_path / "scan"
    legacy = tmp_path / "scan_min_support_2_summary.tsv"
    monkeypatch.setattr("sys.argv", [str(SCRIPT_PATH), "--dbpath", str(db), "--out_prefix", str(prefix),
                                    "--out_tsv", str(legacy), "--min_unit_support", "6",
                                    "--min_lineage_support", "4"])
    assert mod.main() == 0
    pandas.testing.assert_frame_equal(pandas.read_csv(legacy, sep="\t"),
                                     pandas.read_csv(f"{prefix}_all_candidates_summary.tsv", sep="\t"))
    assert len(pandas.read_csv(legacy, sep="\t")) == len(source)
    assert len(pandas.read_csv(f"{prefix}_min_support_6_summary.tsv", sep="\t")) == 2
    filtered = Path(f"{prefix}_min_unit_support_6_min_lineage_support_4_summary.tsv")
    assert pandas.read_csv(filtered, sep="\t")["orthogroup"].tolist() == ["OG4"]
    mod.remove_stale_min_support_sensitivity_outputs(prefix, [])
    assert filtered.is_file()


def test_support_view_cleanup_is_limited_to_the_selected_lineage_bound(tmp_path):
    mod = load_module()
    source = support_view_source()
    prefix = tmp_path / "scan"
    stale = Path(f"{prefix}_min_unit_support_9_min_lineage_support_4_summary.tsv")
    other_lineage = Path(f"{prefix}_min_unit_support_9_min_lineage_support_3_summary.tsv")
    unit_only = Path(f"{prefix}_min_support_9_summary.tsv")
    user_file = Path(f"{prefix}_min_unit_support_9_min_lineage_support_4_notes.txt")
    for path in (stale, other_lineage, unit_only, user_file):
        path.write_text("preserve")
    mod.write_support_views(source, prefix, 6, 4)
    assert not stale.exists()
    assert all(path.read_text() == "preserve" for path in (other_lineage, unit_only, user_file))


def test_nullable_missing_count_is_rejected_and_disabled_lineage_condition_accepts_legacy_input():
    mod = load_module()
    source = support_view_source()
    source["support_lineage_count"] = pandas.Series([pandas.NA, 4, 3, 4], dtype="Int64")
    with pytest.raises(ValueError, match="OG1"):
        mod.support_mask(source, 6, 4, "mixed scan")
    legacy = source.drop(columns="support_lineage_count")
    assert mod.support_mask(legacy, 6, 0, "old scan").tolist() == [False, False, True, True]


@pytest.mark.parametrize("argument", ["--min_unit_support", "--min_lineage_support"])
def test_main_rejects_negative_support_bounds_before_outputs(tmp_path, monkeypatch, argument):
    mod = load_module()
    output = tmp_path / "output.tsv"
    monkeypatch.setattr("sys.argv", [str(SCRIPT_PATH), "--dbpath", str(tmp_path / "missing.db"),
                                    "--out_prefix", str(tmp_path / "scan"), "--out_tsv", str(output), argument, "-1"])
    with pytest.raises(ValueError, match="integer >= 0"):
        mod.main()
    assert not output.exists()


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
                "p_rate_enrichment_asymptotic": 0.01,
                "score_rate_enrichment": 0.03,
                "q_rate_enrichment_asymptotic_global": 0.04,
                "q_rate_enrichment_empirical_by_trait_match": 0.05,
                "lineage_total": 4,
                "support_lineage_count": 3,
                "support_lineage_fraction": 0.75,
                "support_lineage_ids": "7,11,13",
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
        "score_rate_enrichment",
        "q_rate_enrichment_asymptotic_global",
        "q_rate_enrichment_empirical_by_trait_match",
        "lineage_total",
        "support_lineage_count",
        "support_lineage_fraction",
        "support_lineage_ids",
        "future_csubst_metric",
    }
    assert expected_columns.issubset(observed.columns)
    assert observed.loc[0, "scan_rate_exposure"] == "state_aware"
    assert observed.loc[0, "site_rate"] == 0.125
    assert observed.loc[0, "future_csubst_metric"] == 42.0
    assert observed.columns[-1] == "future_csubst_metric"

    ranked, score_column, score_kind = mod.ranked_candidates(observed)
    assert ranked.shape[0] == 1
    assert score_column == "q_rate_enrichment_asymptotic_global"
    assert score_kind == "BH-FDR"
    assert ranked.loc[0, "lineage_total"] == 4
    assert ranked.loc[0, "support_lineage_count"] == 3
    assert ranked.loc[0, "support_lineage_fraction"] == 0.75
    assert ranked.loc[0, "support_lineage_ids"] == "7,11,13"


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
            "q_rate_enrichment_asymptotic_global": [0.01, 0.04, 0.2, 1.0],
            "q_rate_enrichment_empirical_by_trait_match": [0.02, 0.08, 0.4, 1.0],
            "p_rate_enrichment_empirical_maxT": [1.0, 1.0, 1.0, 1.0],
        }
    )

    centers, candidate_counts, series = mod.support_significance_data(frame)

    assert centers.shape == (10,)
    assert candidate_counts.sum() == 4
    observed = {item["method"]["short_label"]: item for item in series}
    assert observed["Analytical"]["significant_counts"].sum() == 2


def test_probability_count_label_reports_three_threshold_counts():
    mod = load_module()
    values = pandas.Series([0.0005, 0.001, 0.01, 0.02, 0.05, 0.2]).to_numpy()

    observed = mod.probability_count_label(mod.PVALUE_QVALUE_METHODS[0], values)

    assert observed == "Analytical: 5 / 3 / 2"


def test_write_min_support_sensitivity_writes_threshold_series(tmp_path):
    mod = load_module()
    frame = pandas.DataFrame(
        {
            "orthogroup": ["OG1", "OG1", "OG2", "OG2"],
            "trait": ["A", "A", "A", "B"],
            "scan_match": ["m1", "m1", "m1", "m2"],
            "support_unit_count": [2, 3, 4, 5],
            "p_rate_enrichment_asymptotic": [0.001, 0.01, 0.03, 0.5],
            "p_rate_enrichment_empirical": [0.002, 0.02, 0.04, 0.8],
            "p_rate_enrichment_empirical_maxT": [0.05, 0.1, 0.2, 1.0],
            "q_rate_enrichment_asymptotic_global": [0.004, 0.03, 0.04, 0.5],
            "q_rate_enrichment_empirical_by_trait_match": [0.008, 0.04, 0.06, 0.8],
            "p_rate_enrichment_bootstrap_maxT": [0.2, 0.3, 0.4, 1.0],
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
    assert manifest["q_rate_enrichment_asymptotic_global_le_0.05"].tolist() == [2, 1, 0]
    assert manifest[f"{mod.SUPPORT_Q_COLUMN}_le_0.05"].tolist() == [2, 0, 0]
    assert not stale_paths["summary_tsv"].exists()
    assert not stale_paths["plot_pdf"].exists()
    assert not (tmp_path / "min_support_sensitivity").exists()
    for threshold, expected_rows in ((3, 3), (4, 2), (5, 1)):
        paths = mod.min_support_sensitivity_paths(out_prefix, threshold)
        subset = pandas.read_csv(paths["summary_tsv"], sep="\t")
        assert subset.shape[0] == expected_rows
        expected = frame.loc[frame["support_unit_count"] >= threshold].reset_index(drop=True)
        pandas.testing.assert_frame_equal(subset.loc[:, frame.columns], expected, check_dtype=False)
        expected_q = {3: [0.03, 0.045, 0.5], 4: [0.06, 0.5], 5: [0.5]}
        assert subset[mod.SUPPORT_Q_COLUMN].tolist() == pytest.approx(expected_q[threshold])
        assert {column for column in subset if column.endswith("_global")} == {"q_rate_enrichment_asymptotic_global"}
        assert (subset["support_unit_count"] >= threshold).all()
        assert set(mod.ORTHOGROUP_BESTHIT_COLUMNS).issubset(subset.columns)
        assert paths["plot_pdf"].is_file()
        assert paths["plot_pdf"].stat().st_size > 1000


def test_write_pvalue_qvalue_distributions(tmp_path):
    mod = load_module()
    frame = pandas.DataFrame(
        {
            "p_rate_enrichment_asymptotic": [0.001, 0.01, 0.2, 1.0],
            "p_rate_enrichment_empirical": [0.002, 0.02, 0.3, 1.0],
            "p_rate_enrichment_empirical_maxT": [0.05, 0.4, 0.9, 1.0],
            "q_rate_enrichment_asymptotic_global": [0.004, 0.03, 0.4, 1.0],
            "q_rate_enrichment_empirical_by_trait_match": [0.006, 0.04, 0.5, 1.0],
            "p_rate_enrichment_bootstrap_maxT": [1.0, 1.0, 1.0, 1.0],
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
            "score_rate_enrichment": [3.0, 2.0, 0.7, 0.0],
            "support_fraction": [0.75, 0.5, 0.4, 0.25],
            "support_unit_count": [3, 2, 2, 1],
            "support_unit_ids": ["1,2,3", "1,2", "2,3", "3"],
            "p_rate_enrichment_asymptotic": [0.001, 0.01, 0.2, 1.0],
            "p_rate_enrichment_empirical": [0.002, 0.02, 0.3, 1.0],
            "p_rate_enrichment_empirical_maxT": [0.05, 0.4, 0.9, 1.0],
            "q_rate_enrichment_asymptotic_global": [0.004, 0.03, 0.4, 1.0],
            "q_rate_enrichment_empirical_by_trait_match": [0.006, 0.04, 0.5, 1.0],
            "p_rate_enrichment_bootstrap_maxT": [1.0, 1.0, 1.0, 1.0],
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
    assert len(primary) == 4  # Legacy all-candidate output remains complete.
    filtered = pandas.read_csv(f"{out_prefix}_min_unit_support_2_min_lineage_support_0_summary.tsv", sep="\t")
    assert len(filtered) == 3
    assert filtered["support_unit_count"].ge(2).all()
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


def test_summary_rejects_legacy_database_instead_of_using_old_global_q(tmp_path):
    mod = load_module()
    with sqlite3.connect(tmp_path / "legacy.sqlite3") as conn:
        pandas.DataFrame({"p_rate_enrichment": [0.01], "q_rate_enrichment_global": [0.02]}).to_sql("aa_change", conn, index=False)
        with pytest.raises(ValueError, match="rerun scan and rebuild"):
            mod.read_table(conn, "aa_change")


def test_ranking_uses_stable_score_even_if_calibration_is_unavailable():
    mod = load_module()
    frame = pandas.DataFrame({"score_rate_enrichment": [3.0, 8.0, float("nan")],
                              "p_rate_enrichment_asymptotic": [0.001, 1e-8, float("nan")],
                              "p_rate_enrichment_empirical_maxT": [float("nan")] * 3})
    ranked, column, kind = mod.ranked_candidates(frame)
    assert (column, kind) == ("score_rate_enrichment", "Score")
    assert ranked.index.tolist() == [1, 0, 2]
    assert ranked["p_rate_enrichment_empirical_maxT"].isna().all()


@pytest.mark.parametrize('value', [-0.01, 1.01, float('inf'), 'invalid'])
def test_summary_invalid_p_preflight_preserves_outputs(tmp_path, monkeypatch, value):
    mod = load_module()
    source = support_view_source().astype({'p_rate_enrichment_asymptotic': object})
    source.loc[0, 'p_rate_enrichment_asymptotic'] = value
    db = tmp_path / 'scan.db'
    with sqlite3.connect(db) as conn:
        source.to_sql('aa_change', conn, index=False)
    prefix = tmp_path / 'scan'
    output = tmp_path / 'scan_min_support_2_summary.tsv'
    output.write_text('preserve\n')
    monkeypatch.setattr('sys.argv', [str(SCRIPT_PATH), '--dbpath', str(db), '--out_prefix', str(prefix),
                                   '--out_tsv', str(output), '--min_unit_support', '6', '--min_lineage_support', '4'])
    with pytest.raises(ValueError, match='invalid probabilities'):
        mod.main()
    assert output.read_text() == 'preserve\n'


def test_summary_filtered_plots_use_filtered_q(tmp_path):
    mod = load_module()
    source = support_view_source()
    source['q_rate_enrichment_asymptotic_global'] = 1.0
    selected = mod.support_filtered_bh(source, 5, 4, 'source')
    assert mod.choose_score_column(selected)[0] == mod.SUPPORT_Q_COLUMN
    series = mod.probability_series(selected, 'q_column')
    assert series[0][1] == mod.SUPPORT_Q_COLUMN
    assert series[0][2].tolist() == [0.004, 0.004]
    assert '2 finite tests' in mod.bh_caption(selected)
    assert 'unit >= 5 and lineage >= 4' in mod.bh_caption(selected)
    selected.attrs.clear()  # A reread TSV still has enough metadata for captions.
    assert '2 finite tests' in mod.bh_caption(selected)


@pytest.mark.parametrize('bad_count', [-1, 1.5, 'invalid', None])
def test_disabled_unit_bound_validates_counts_before_output_replacement(tmp_path, monkeypatch, bad_count):
    mod = load_module()
    source = support_view_source().astype({'support_unit_count': object})
    source.loc[0, 'support_unit_count'] = bad_count
    db = tmp_path / 'scan.db'
    with sqlite3.connect(db) as conn:
        source.to_sql('aa_change', conn, index=False)
    prefix = tmp_path / 'scan'
    output = tmp_path / 'scan_min_support_2_summary.tsv'
    output.write_text('preserve\n')
    monkeypatch.setattr('sys.argv', [str(SCRIPT_PATH), '--dbpath', str(db), '--out_prefix', str(prefix),
                                   '--out_tsv', str(output), '--min_unit_support', '0'])
    with pytest.raises(ValueError, match='support_unit_count'):
        mod.main()
    assert output.read_text() == 'preserve\n'


@pytest.mark.parametrize('units,lineages', [(0, 0), (5, 4), (6, 0), (0, 4), (9, 6)])
def test_filtered_bh_matches_independent_scipy_reference(units, lineages):
    import numpy as np
    from scipy.stats import false_discovery_control

    mod = load_module()
    rng = np.random.default_rng(545)
    counts = rng.integers(1, 9, 400)
    source = pandas.DataFrame({
        'support_unit_count': counts,
        'support_lineage_count': [rng.integers(1, min(count, 5) + 1) for count in counts],
        'orthogroup': [f'OG{i}' for i in range(400)],
        'trait': [f'trait{i % 3}' for i in range(400)],
        'scan_match': [f'm{i % 2}' for i in range(400)],
        'p_rate_enrichment_asymptotic': rng.random(400),
        'q_rate_enrichment_asymptotic_global': 0.2,
    })
    source.loc[:19, 'p_rate_enrichment_asymptotic'] = np.nan
    source.loc[20:29, 'p_rate_enrichment_asymptotic'] = 0.0
    source.loc[30:39, 'p_rate_enrichment_asymptotic'] = 1.0
    source.loc[40:49, 'p_rate_enrichment_asymptotic'] = 0.001  # Ties.
    source.loc[50:99, 'p_rate_enrichment_asymptotic'] = np.geomspace(1e-280, 1e-5, 50)
    original = source.copy()
    output = mod.support_filtered_bh(source, units, lineages, 'randomized reference')
    finite = output['p_rate_enrichment_asymptotic'].notna()
    expected = false_discovery_control(output.loc[finite, 'p_rate_enrichment_asymptotic'].to_numpy(), method='bh') if finite.any() else []
    np.testing.assert_allclose(output.loc[finite, mod.SUPPORT_Q_COLUMN], expected, rtol=1e-13, atol=0)
    assert output.loc[~finite, mod.SUPPORT_Q_COLUMN].isna().all()
    assert output.attrs['support_bh_metadata']['support_bh_test_count'] == int(finite.sum())
    pandas.testing.assert_frame_equal(source, original)
