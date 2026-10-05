import argparse
import io
import json
import os
import subprocess
import sys
import zipfile
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas as pd
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "csubst_scan_candidate_sites.py"


def load_module():
    spec = spec_from_file_location("csubst_scan_candidate_sites", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("units,lineages,expected", [
    (6, 4, ["OG0003"]), (5, 4, ["OG0002", "OG0003"]), (6, 0, ["OG0001", "OG0003"]),
    (0, 4, ["OG0002", "OG0003"]), (0, 0, ["OG0002", "OG0001", "OG0003"]), (8, 4, []),
])
def test_candidate_selection_applies_filtered_bh_and_preserves_source_probabilities(tmp_path, units, lineages, expected):
    mod = load_module()
    source = candidate_rows()
    source["support_lineage_count"] = [4, 3, 4]
    source[mod.DEFAULT_PROBABILITY_COLUMN] = [0.02, 0.01, 0.03]
    summary = tmp_path / "all_candidates.tsv"
    write_summary(summary, source)
    selected = mod.load_threshold_candidates(summary, units, mod.DEFAULT_PROBABILITY_COLUMN, 0.05, 0,
                                             "no", "none", lineages, filter_unit_support=True)
    assert selected["orthogroup"].tolist() == expected
    for row in selected.itertuples():
        original = source.loc[source["orthogroup"].eq(row.orthogroup)].iloc[0]
        assert row.q_rate_enrichment_asymptotic_global == original["q_rate_enrichment_asymptotic_global"]
        assert row.p_rate_enrichment_asymptotic == original["p_rate_enrichment_asymptotic"]
    assert selected["_selection_min_lineage_support"].eq(lineages).all()
    expected_q = {(6, 4): [0.02], (5, 4): [0.002, 0.02], (6, 0): [0.004, 0.02],
                  (0, 4): [0.002, 0.02], (0, 0): [0.003, 0.003, 0.02], (8, 4): []}
    assert selected[mod.DEFAULT_PROBABILITY_COLUMN].tolist() == pytest.approx(expected_q[(units, lineages)])
    assert selected.attrs["support_bh_metadata"]["support_bh_test_count"] == len(expected)


def test_candidate_cap_applies_after_lineage_filter_and_preserves_analysis_identity(tmp_path):
    mod = load_module()
    source = candidate_rows().iloc[:2].copy()
    source["support_lineage_count"] = [4, 3]
    summary = tmp_path / "summary.tsv"
    write_summary(summary, source)
    unfiltered = mod.load_threshold_candidates(summary, 5, mod.DEFAULT_PROBABILITY_COLUMN, 0.05, 0, "no", "none")
    selected = mod.load_threshold_candidates(summary, 5, mod.DEFAULT_PROBABILITY_COLUMN, 0.05, 1, "no", "none", 4)
    assert selected["orthogroup"].tolist() == ["OG0002"]
    same_candidate = unfiltered.loc[unfiltered["orthogroup"].eq("OG0002")].iloc[0]
    assert selected.iloc[0]["_analysis_key"] == same_candidate["_analysis_key"]


@pytest.mark.parametrize("count", [None, float("inf"), -1, 1.5, "unknown", "missing_column"])
def test_lineage_filter_rejects_invalid_counts_even_outside_selected_pool(tmp_path, count):
    mod = load_module()
    source = candidate_rows()
    source["support_lineage_count"] = pd.Series([4, 4, 4], dtype=object)
    if count == "missing_column":
        source = source.drop(columns="support_lineage_count")
    else:
        source.loc[0, "support_lineage_count"] = count
    summary = tmp_path / "summary.tsv"
    write_summary(summary, source)
    with pytest.raises(ValueError, match="support_lineage_count"):
        mod.load_threshold_candidates(summary, 6, mod.DEFAULT_PROBABILITY_COLUMN, 0.001, 1, "no", "none", 4,
                                      filter_unit_support=True)


def test_full_source_discovery_supports_disabled_unit_filter_and_ignores_filtered_views(tmp_path):
    mod = load_module()
    prefix = tmp_path / "scan"
    full = Path(f"{prefix}_all_candidates_summary.tsv")
    write_summary(full, candidate_rows())
    Path(f"{prefix}_min_unit_support_6_min_lineage_support_4_summary.tsv").write_text("irrelevant\n")
    assert mod.discover_summary_tables(prefix, 0) == {0: full}
    assert list(mod.discover_summary_tables(prefix, 5)) == [7, 6, 5]
    assert list(mod.discover_summary_tables(prefix, 1)) == list(range(7, 0, -1))
    assert mod.discover_summary_tables(prefix, 8) == {}
    # Even an empty selected pool cannot silently conceal a missing lineage column.
    with pytest.raises(ValueError, match="support_lineage_count"):
        mod.discover_summary_tables(prefix, 8, 4)


def test_complete_unit_series_retains_legacy_sources_for_archive_reuse(tmp_path):
    mod = load_module()
    prefix = tmp_path / "orthogroup_csubst_aa_change"
    source = candidate_rows()
    source["support_lineage_count"] = 4
    full = Path(f"{prefix}_all_candidates_summary.tsv")
    write_summary(full, source)
    write_run_summaries(tmp_path)
    assert mod.discover_summary_tables(prefix, 5) == {
        value: Path(f"{prefix}_min_support_{value}_summary.tsv") for value in (7, 6, 5)
    }
    assert mod.discover_summary_tables(prefix, 0, 4) == {0: full}
    below_legacy_minimum = mod.discover_summary_tables(prefix, 1)
    assert below_legacy_minimum[1] == full
    assert below_legacy_minimum[5] == Path(f"{prefix}_min_support_5_summary.tsv")


def test_legacy_lineage_preflight_validates_counts_when_unit_bound_selects_no_threshold(tmp_path):
    mod = load_module()
    write_run_summaries(tmp_path)
    prefix = tmp_path / "orthogroup_csubst_aa_change"
    assert mod.discover_summary_tables(prefix, 8) == {}
    with pytest.raises(ValueError, match="support_lineage_count"):
        mod.discover_summary_tables(prefix, 8, 4)


def test_lineage_preflight_preserves_existing_archives_and_manifests(tmp_path):
    mod = load_module()
    args = make_run_args(tmp_path)
    args.min_lineage_support = 4
    write_run_summaries(tmp_path)
    for threshold in (6, 7):
        path = tmp_path / f"orthogroup_csubst_aa_change_min_support_{threshold}_summary.tsv"
        frame = pd.read_csv(path, sep="\t")
        frame["support_lineage_count"] = 4
        write_summary(path, frame)
    output = Path(args.out_dir)
    output.mkdir()
    suffix = mod.output_suffix(args.probability_column, args.probability_threshold, 0, "no", "none", 4)
    existing = [output / f"orthogroup_csubst_aa_change_candidate_sites_{suffix}_manifest.tsv",
                output / f"orthogroup_csubst_aa_change_candidate_sites_min_support_5_{suffix}.zip"]
    for path in existing:
        path.write_bytes(b"preserve")
    with pytest.raises(ValueError, match="support_lineage_count"):
        mod.run(args)
    assert all(path.read_bytes() == b"preserve" for path in existing)
    assert not list(output.glob("*.work"))


@pytest.mark.parametrize("minimum,expected_thresholds", [(0, [0]), (5, [7, 6, 5])])
def test_run_full_candidate_source_uses_independent_lineage_bound(tmp_path, monkeypatch, minimum, expected_thresholds):
    mod = load_module()
    args = make_run_args(tmp_path)
    args.min_support = minimum
    args.min_lineage_support = 4
    source = candidate_rows()
    source["support_lineage_count"] = [4, 3, 4]
    source[args.probability_column] = 0.01
    write_summary(Path(f"{args.summary_prefix}_all_candidates_summary.tsv"), source)
    selected_by_threshold = {}
    def fake_package(candidates, threshold, minimum_lineage_support, **kwargs):
        assert minimum_lineage_support == 4
        selected_by_threshold[threshold] = candidates["orthogroup"].tolist()
    monkeypatch.setattr(mod, "package_threshold", fake_package)
    monkeypatch.setattr(mod, "archive_matches_source", lambda *args, **kwargs: False)
    monkeypatch.setattr(mod, "inspect_required_report_inputs", available_input_states)
    monkeypatch.setattr(mod, "ensure_candidate_analyses", lambda **kwargs: [])
    manifest = pd.read_csv(mod.run(args), sep="\t")
    assert manifest["min_support"].tolist() == expected_thresholds
    assert manifest["min_lineage_support"].eq(4).all()
    assert all("_min_lineage_support_4" in name for name in manifest["archive_zip"])
    if minimum == 0:
        assert selected_by_threshold == {0: ["OG0002", "OG0003"]}
    else:
        assert selected_by_threshold == {7: ["OG0003"], 6: ["OG0003"], 5: ["OG0002", "OG0003"]}


def test_lineage_filtered_zip_records_and_validates_both_bounds(tmp_path):
    mod = load_module()
    source = candidate_rows().iloc[[1]].copy()
    source["support_lineage_count"] = 4
    source["support_lineage_ids"] = "001,7,11,13"
    summary = tmp_path / "summary.tsv"
    write_summary(summary, source)
    candidates = mod.load_threshold_candidates(summary, 5, mod.DEFAULT_PROBABILITY_COLUMN, 0.05, 0, "no", "none", 4)
    candidates["_required_input_signature"] = "input-signature"
    row = candidates.iloc[0]
    cache = tmp_path / "cache"
    make_candidate_cache(mod, cache, row)
    archive = mod.archive_path_for_threshold(tmp_path / "scan", tmp_path, 5, mod.DEFAULT_PROBABILITY_COLUMN,
                                             0.05, 0, "no", "none", 4)
    assert archive.name.endswith("_min_lineage_support_4.zip")
    mod.package_threshold(candidates, 5, summary, archive, tmp_path / "packages", cache,
                          mod.DEFAULT_PROBABILITY_COLUMN, 0.05, minimum_lineage_support=4)
    assert mod.archive_matches_source(archive, summary, minimum_support=5, minimum_lineage_support=4)
    assert not mod.archive_matches_source(archive, summary, minimum_support=5, minimum_lineage_support=3)
    assert not mod.archive_matches_source(archive, summary, minimum_support=6, minimum_lineage_support=4)
    with zipfile.ZipFile(archive) as zipped:
        root = archive.stem
        metadata = pd.read_csv(zipped.open(f"{root}/package_metadata.tsv"), sep="\t").iloc[0]
        manifest = pd.read_csv(zipped.open(f"{root}/candidate_manifest.tsv"), sep="\t").iloc[0]
        candidate = pd.read_csv(zipped.open(f"{root}/{manifest['candidate_tsv']}"), sep="\t").iloc[0]
        assert metadata["min_support"] == manifest["selection_min_support"] == candidate["selection_min_support"] == 5
        assert metadata["min_lineage_support"] == manifest["selection_min_lineage_support"] == candidate["selection_min_lineage_support"] == 4
        assert candidate[mod.DEFAULT_PROBABILITY_COLUMN] == 0.002
        assert candidate["q_rate_enrichment_asymptotic_global"] == 0.01
        assert metadata["support_bh_test_count"] == 1
        assert metadata["support_bh_candidate_count"] == 1
        assert metadata["probability_policy"] == "bh_after_unit_and_lineage_support_filter_v1"
        assert "support_lineage_count >= 4" in zipped.read(f"{root}/README.txt").decode()
    assert "min_lineage_support: 4" in mod.candidate_annotation_text(row, mod.DEFAULT_PROBABILITY_COLUMN, 0.05)


@pytest.mark.parametrize("parameter", ["min_support", "min_lineage_support"])
def test_validate_args_rejects_negative_support_bounds(tmp_path, parameter):
    mod = load_module()
    args = make_run_args(tmp_path)
    setattr(args, parameter, -1)
    with pytest.raises(ValueError, match="integer >= 0"):
        mod.validate_args(args)


def candidate_rows():
    return pd.DataFrame(
        [
            {
                "orthogroup": "OG0002",
                "trait": "aquatic",
                "state_change": "A>V",
                "codon_site_alignment": 9,
                "support_unit_count": 5,
                "support_unit_ids": "1,2,3,4,5",
                "support_branch_ids": "7, 3,7",
                "p_rate_enrichment_asymptotic": 0.001,
                "scan_calibration_status": "conditional_assignment",
                "q_rate_enrichment_asymptotic_global": 0.02,
                "besthit_0.05": "protein B",
            },
            {
                "orthogroup": "OG0001",
                "trait": "aquatic",
                "state_change": "G>D",
                "codon_site_alignment": 4,
                "support_unit_count": 6,
                "support_unit_ids": "1,2,3,4,5,6",
                "support_branch_ids": "8,9",
                "p_rate_enrichment_asymptotic": 0.002,
                "scan_calibration_status": "conditional_assignment",
                "q_rate_enrichment_asymptotic_global": 0.01,
                "besthit_0.05": "protein A",
            },
            {
                "orthogroup": "OG0003",
                "trait": "aquatic",
                "state_change": "L>F",
                "codon_site_alignment": 12,
                "support_unit_count": 7,
                "support_unit_ids": "1,2,3,4,5,6,7",
                "support_branch_ids": "10,11",
                "p_rate_enrichment_asymptotic": 0.02,
                "scan_calibration_status": "conditional_assignment",
                "q_rate_enrichment_asymptotic_global": 0.2,
                "besthit_0.05": "protein C",
            },
        ]
    )


@pytest.mark.parametrize("identifier_dtype", ["object", "string"])
@pytest.mark.parametrize("duplicate_metadata", [False, True])
def test_input_state_annotation_preserves_exact_cache_keys_frame_types_and_order(identifier_dtype, duplicate_metadata):
    mod = load_module()
    frame = pd.DataFrame({
        "orthogroup": pd.Series(["OG1", "OG2", "OG1"], dtype=identifier_dtype),
        "_analysis_key": pd.Series(["keyA", "keyB", "keyC"], dtype=identifier_dtype),
        "_candidate_id": pd.Series(["candA", "candB", "candC"], dtype=identifier_dtype),
        "metadata": pd.Series([1, pd.NA, 3], dtype="Int64"),
    })
    if duplicate_metadata:
        frame = pd.concat([frame, frame[["metadata"]]], axis=1)
    frame.index = pd.Index([8, 8, 2], name="source_row")
    original = frame.copy(deep=True)
    states = {"OG1": {"missing_required_inputs": ["a.tsv", "b.nwk"], "required_input_signature": "signatureA"},
              "OG2": {"missing_required_inputs": [], "required_input_signature": "signatureB"}}
    output = mod.annotate_candidate_input_state(frame, states)
    expected_keys = ["27da8d6f06fff3753dfb091053cfd4449be6a1b174b28b68c6f75eb8aad4cc39",
                     "73453a6f3812c1be649f825128d05806b4614accec8450fe867bc2df7ed3f5f5",
                     "4e3f2f6309aed8ba2d72c2af748a827b506d9ef0af09c587ab3d913e5ad1dafc"]
    assert output["_analysis_key"].tolist() == expected_keys
    assert output["_cache_name"].tolist() == [f"{name}_{key[:16]}" for name, key in zip(["candA", "candB", "candC"], expected_keys, strict=True)]
    assert output["_missing_required_inputs"].tolist() == ["a.tsv;b.nwk", "", "a.tsv;b.nwk"]
    assert output["_required_input_signature"].tolist() == ["signatureA", "signatureB", "signatureA"]
    assert output.columns.tolist() == frame.columns.tolist() + ["_missing_required_inputs", "_required_input_signature", "_cache_name"]
    pd.testing.assert_frame_equal(output.drop(columns=["_analysis_key", "_missing_required_inputs", "_required_input_signature", "_cache_name"]),
                                  frame.drop(columns="_analysis_key"))
    pd.testing.assert_frame_equal(frame, original)


def test_input_state_annotation_keeps_numeric_legacy_row_coercion():
    mod = load_module()
    frame = pd.DataFrame({"orthogroup": [1], "_analysis_key": [2.0], "_candidate_id": [3]})
    original = frame.copy(deep=True)
    output = mod.annotate_candidate_input_state(frame, {"1.0": {"missing_required_inputs": [], "required_input_signature": "signatureN"}})
    assert output["_analysis_key"].tolist() == ["6890eddd8869214a17e0d6f1d0346476386df22acfce6a8a3e8b888c3b55cbf0"]
    assert output["_cache_name"].tolist() == ["3.0_6890eddd8869214a"]
    pd.testing.assert_frame_equal(frame, original)


def test_input_state_annotation_keeps_empty_schema_and_lookup_error_order():
    mod = load_module()
    empty = pd.DataFrame()
    output = mod.annotate_candidate_input_state(empty, {})
    assert output.columns.tolist() == ["_missing_required_inputs", "_required_input_signature"]
    assert output.empty and empty.empty and len(empty.columns) == 0
    frame = pd.DataFrame({"orthogroup": ["unknown"]})
    original = frame.copy(deep=True)
    with pytest.raises(KeyError, match="unknown"):
        mod.annotate_candidate_input_state(frame, {})
    pd.testing.assert_frame_equal(frame, original)
    frame = pd.DataFrame({"orthogroup": ["known"], "_analysis_key": ["key"]})
    with pytest.raises(KeyError, match="_candidate_id"):
        mod.annotate_candidate_input_state(frame, {"known": {"missing_required_inputs": [], "required_input_signature": "signature"}})


@pytest.mark.parametrize("numeric_dtype", ["float32", "float64", "Float32"])
@pytest.mark.parametrize("duplicate_metadata", [False, True])
def test_candidate_identity_keeps_scalar_precision_nullable_types_and_exact_hashes(numeric_dtype, duplicate_metadata):
    mod = load_module()
    frame = pd.DataFrame({
        "orthogroup": ["OG_001", "NA"], "trait": pd.Series(["aquatic", pd.NA], dtype="string"),
        "state_change": ["A>V", "G>D"], "codon_site_alignment": pd.Series([9, 4], dtype="Int64"),
        "from_state": pd.Series([0.1, float("nan")], dtype=numeric_dtype),
        "to_state": pd.Series([1, 2], dtype="UInt64"), "_canonical_support_branch_ids": ["3,7", None],
        "unused": pd.Categorical(["u", "v"]),
    })
    if duplicate_metadata:
        frame = pd.concat([frame, frame[["unused"]]], axis=1)
    frame.index = pd.Index([8, 8], name="source_row")
    original = frame.copy(deep=True)
    output = mod.assign_candidate_ids(frame, "no", "none")
    if numeric_dtype == "float64":
        first_id = "OG_001_site9_A_V_bffa70544fc05999"
        first_key = "c5176013feda1aab58d957f5fdec9cc2507e576c433f599401f73de51f719f5c"
    else:
        first_id = "OG_001_site9_A_V_3bc868909033769b"
        first_key = "c3494a1c7b77f5d57fb3cac19bf0663c0cd4c06ecead8449ba6b69b1f3b7e100"
    expected_ids = [first_id, "NA_site4_G_D_8faacf8f38690d6d"]
    expected_keys = [first_key, "af336d6fb3bf5090ac8403c06eb2d1d977e47d3f6cfed19b38caaac3360aec39"]
    assert output["_candidate_id"].tolist() == expected_ids
    assert output["_analysis_key"].tolist() == expected_keys
    assert output["_cache_name"].tolist() == [f"{identifier}_{key[:16]}" for identifier, key in zip(expected_ids, expected_keys, strict=True)]
    pd.testing.assert_frame_equal(output.drop(columns=["_candidate_id", "_analysis_key", "_cache_name"]), original)
    pd.testing.assert_frame_equal(frame, original)


def test_candidate_identity_keeps_numeric_legacy_row_coercion():
    mod = load_module()
    frame = pd.DataFrame({"orthogroup": [1], "state_change": [2], "codon_site_alignment": [3], "unused": [0.5]})
    output = mod.assign_candidate_ids(frame, "no", "none")
    assert output["_candidate_id"].tolist() == ["1.0_site3_2.0_115c165d73c91d2d"]
    assert output["_analysis_key"].tolist() == ["848090a1a470c3c8b6a0d18122aa9fdbe36a90f369df6f87888ad51bbffdb7a0"]


def test_candidate_identity_keeps_duplicate_rejection_empty_schema_and_error_order():
    mod = load_module()
    empty = mod.assign_candidate_ids(pd.DataFrame(), "no", "none")
    assert empty.empty and empty.columns.tolist() == ["_candidate_id", "_analysis_key", "_cache_name"]
    with pytest.raises(KeyError, match="codon_site_alignment"):
        mod.assign_candidate_ids(pd.DataFrame({"orthogroup": ["OG1"]}), "no", "none")
    with pytest.raises(ValueError, match="cannot convert float NaN to integer"):
        mod.assign_candidate_ids(pd.DataFrame({"orthogroup": ["OG1"], "codon_site_alignment": [float("nan")]}), "no", "none")
    frame = pd.DataFrame({"orthogroup": ["OG1", "OG1"], "codon_site_alignment": [1, 1], "state_change": ["A>V", "A>V"]})
    original = frame.copy(deep=True)
    with pytest.raises(ValueError, match="Candidate IDs are not unique"):
        mod.assign_candidate_ids(frame, "no", "none")
    pd.testing.assert_frame_equal(frame, original)


def write_summary(path, frame=None):
    if frame is None:
        frame = candidate_rows()
    frame.to_csv(path, sep="\t", index=False)


def test_discover_summary_tables_returns_contiguous_thresholds_in_descending_order(tmp_path):
    mod = load_module()
    prefix = tmp_path / "orthogroup_csubst_aa_change"
    for threshold in (5, 6, 7):
        write_summary(tmp_path / f"{prefix.name}_min_support_{threshold}_summary.tsv")

    discovered = mod.discover_summary_tables(prefix, 5)

    assert list(discovered) == [7, 6, 5]


def test_discover_summary_tables_rejects_a_gap(tmp_path):
    mod = load_module()
    prefix = tmp_path / "orthogroup_csubst_aa_change"
    for threshold in (5, 7):
        write_summary(tmp_path / f"{prefix.name}_min_support_{threshold}_summary.tsv")

    with pytest.raises(FileNotFoundError, match="Missing threshold.*6"):
        mod.discover_summary_tables(prefix, 5)


def test_discover_summary_tables_returns_empty_when_observed_max_is_below_minimum(tmp_path):
    mod = load_module()
    prefix = tmp_path / "orthogroup_csubst_aa_change"
    for threshold in (2, 3, 4):
        write_summary(tmp_path / f"{prefix.name}_min_support_{threshold}_summary.tsv")

    assert mod.discover_summary_tables(prefix, 5) == {}


def test_load_threshold_candidates_filters_q_and_canonicalizes_branches(tmp_path):
    mod = load_module()
    summary = tmp_path / "summary.tsv"
    write_summary(summary)

    selected = mod.load_threshold_candidates(
        summary_path=summary,
        minimum_support=5,
        probability_column="q_rate_enrichment_asymptotic_global",
        probability_threshold=0.05,
        max_candidates=0,
        csubst_nonsyn_recode="no",
        pdb="none",
    )

    assert selected["orthogroup"].tolist() == ["OG0001", "OG0002"]
    assert selected["_canonical_support_branch_ids"].tolist() == ["8,9", "3,7"]
    assert selected["_candidate_rank"].tolist() == [1, 2]
    assert selected["_candidate_id"].str.contains(r"_[0-9a-f]{16}$").all()


@pytest.mark.parametrize("identifier", ["001", "9007199254740993", "18446744073709551615"])
def test_candidate_reader_preserves_single_foreground_id_as_text(tmp_path, identifier):
    mod = load_module()
    source = candidate_rows().iloc[[1]].copy()
    source["support_lineage_ids"] = identifier
    source["support_lineage_count"] = 1
    source["lineage_total"] = 2
    source["support_lineage_fraction"] = 0.5
    summary = tmp_path / "summary.tsv"
    write_summary(summary, source)
    selected = mod.load_threshold_candidates(
        summary, 5, "q_rate_enrichment_asymptotic_global", 0.05, 0, "no", "none",
    )
    assert selected.loc[0, "support_lineage_ids"] == identifier
    row = selected.iloc[0]
    assert f"Support foreground lineage IDs: {identifier}" in mod.candidate_annotation_text(
        row, "q_rate_enrichment_asymptotic_global", 0.05,
    )
    output = mod.candidate_output_frame(row, "q_rate_enrichment_asymptotic_global", 0.05)
    assert output.loc[0, "support_lineage_ids"] == identifier


def test_candidate_analysis_identity_is_stable_across_recalculated_probability_values(tmp_path):
    mod = load_module()
    first_path = tmp_path / "first.tsv"
    second_path = tmp_path / "second.tsv"
    first = candidate_rows().iloc[[0]].copy()
    second = first.copy()
    second["q_rate_enrichment_asymptotic_global"] = 0.03
    write_summary(first_path, first)
    write_summary(second_path, second)

    loaded = [
        mod.load_threshold_candidates(
            summary_path=path,
            minimum_support=5,
            probability_column="q_rate_enrichment_asymptotic_global",
            probability_threshold=0.05,
            max_candidates=0,
            csubst_nonsyn_recode="no",
            pdb="none",
        )
        for path in (first_path, second_path)
    ]

    assert loaded[0].loc[0, "_analysis_key"] == loaded[1].loc[0, "_analysis_key"]
    assert loaded[0].loc[0, "_candidate_id"] == loaded[1].loc[0, "_candidate_id"]


def test_candidate_analysis_identity_ignores_tool_versions_but_tracks_parameters(
    monkeypatch, tmp_path
):
    mod = load_module()
    summary = tmp_path / "summary.tsv"
    write_summary(summary, candidate_rows().iloc[[0]].copy())

    first = mod.load_threshold_candidates(
        summary_path=summary,
        minimum_support=5,
        probability_column="q_rate_enrichment_asymptotic_global",
        probability_threshold=0.05,
        max_candidates=0,
        csubst_nonsyn_recode="no",
        pdb="none",
    )
    monkeypatch.setattr(mod, "analysis_engine_signature", lambda: "new-engine")
    monkeypatch.setattr(mod, "csubst_version", lambda: "new-version")
    second = mod.load_threshold_candidates(
        summary_path=summary,
        minimum_support=5,
        probability_column="q_rate_enrichment_asymptotic_global",
        probability_threshold=0.05,
        max_candidates=0,
        csubst_nonsyn_recode="no",
        pdb="none",
    )
    changed_parameter = mod.load_threshold_candidates(
        summary_path=summary,
        minimum_support=5,
        probability_column="q_rate_enrichment_asymptotic_global",
        probability_threshold=0.05,
        max_candidates=0,
        csubst_nonsyn_recode="dayhoff6",
        pdb="none",
    )

    assert first.loc[0, "_analysis_key"] == second.loc[0, "_analysis_key"]
    assert first.loc[0, "_analysis_key"] != changed_parameter.loc[0, "_analysis_key"]


def test_trait_color_paths_do_not_collide_after_filename_sanitizing(tmp_path):
    mod = load_module()
    trait_file = tmp_path / "traits.tsv"
    pd.DataFrame(
        {
            "species": ["sp1", "sp2"],
            "trait/a": [1, 0],
            "trait?a": [0, 1],
        }
    ).to_csv(trait_file, sep="\t", index=False)

    paths = mod.write_trait_color_tables(
        trait_file,
        ["trait/a", "trait?a"],
        tmp_path / "colors",
    )

    assert paths["trait/a"] != paths["trait?a"]
    first = pd.read_csv(paths["trait/a"], sep="\t")
    second = pd.read_csv(paths["trait?a"], sep="\t")
    assert first["color"].tolist() == ["firebrick", "black"]
    assert second["color"].tolist() == ["black", "firebrick"]


def test_exclusive_run_lock_recovers_a_dead_same_host_owner(tmp_path):
    mod = load_module()
    lock_path = tmp_path / ".candidate.lock"
    lock_directory = Path(f"{lock_path}.d")
    lock_directory.mkdir()
    (lock_directory / "owner.json").write_text(
        json.dumps(
            {
                "hostname": mod.socket.gethostname(),
                "pid": 2_000_000_000,
                "token": "dead-owner",
                "created_unix": 0,
            }
        ),
        encoding="utf-8",
    )

    with mod.exclusive_run_lock(lock_path):
        current_owner = json.loads(
            (lock_directory / "owner.json").read_text(encoding="utf-8")
        )
        assert current_owner["pid"] == os.getpid()
        assert current_owner["token"] != "dead-owner"

    assert not lock_directory.exists()


def test_run_lock_uses_heartbeat_to_detect_a_stale_foreign_owner(tmp_path):
    mod = load_module()
    lock_directory = tmp_path / ".candidate.lock.d"
    lock_directory.mkdir()
    owner_path = lock_directory / "owner.json"
    owner_path.write_text(
        json.dumps(
            {
                "hostname": "another-container",
                "pid": 1,
                "token": "foreign-owner",
                "created_unix": mod.time.time(),
            }
        ),
        encoding="utf-8",
    )

    assert not mod.run_lock_is_stale(lock_directory)
    stale_time = mod.time.time() - mod.RUN_LOCK_STALE_SECONDS - 1
    os.utime(owner_path, (stale_time, stale_time))
    assert mod.run_lock_is_stale(lock_directory)


def test_load_threshold_candidates_rejects_invalid_branch_ids(tmp_path):
    mod = load_module()
    summary = tmp_path / "summary.tsv"
    frame = candidate_rows().iloc[[0]].copy()
    frame["support_branch_ids"] = "3,bad"
    write_summary(summary, frame)

    with pytest.raises(ValueError, match="Invalid branch ID"):
        mod.load_threshold_candidates(
            summary_path=summary,
            minimum_support=5,
            probability_column="q_rate_enrichment_asymptotic_global",
            probability_threshold=0.05,
            max_candidates=0,
            csubst_nonsyn_recode="no",
            pdb="none",
        )


def test_required_report_input_preflight_reports_and_recovers_missing_stat_branch(
    tmp_path,
):
    mod = load_module()
    family_dir = tmp_path / "orthogroup"
    for subdir in ("iqtree_anc", "stat_branch", "clipkit"):
        (family_dir / subdir).mkdir(parents=True)
    for orthogroup in ("OG0001", "OG0002"):
        (family_dir / "iqtree_anc" / f"{orthogroup}_iqtree.anc.zip").write_bytes(
            b"placeholder"
        )
        (family_dir / "clipkit" / f"{orthogroup}_cds.clipkit.fa").write_text(
            ">sp1\nAAA\n", encoding="utf-8"
        )
    (family_dir / "stat_branch" / "OG0001_stat.branch.tsv").write_text(
        "branch_id\n0\n", encoding="utf-8"
    )
    missing_stat = family_dir / "stat_branch" / "OG0002_stat.branch.tsv"
    missing_stat.write_bytes(b"")

    first = mod.inspect_required_report_inputs(
        family_dir, ["OG0001", "OG0002"]
    )

    assert first["OG0001"]["missing_required_inputs"] == []
    assert first["OG0002"]["missing_required_inputs"] == [
        "stat_branch/OG0002_stat.branch.tsv"
    ]
    first_signature = first["OG0002"]["required_input_signature"]

    missing_stat.write_text("branch_id\n0\n", encoding="utf-8")
    second = mod.inspect_required_report_inputs(family_dir, ["OG0002"])

    assert second["OG0002"]["missing_required_inputs"] == []
    assert second["OG0002"]["required_input_signature"] != first_signature


def make_candidate_cache(mod, cache_root, row):
    cache_dir = cache_root / row["_cache_name"]
    cache_dir.mkdir(parents=True)
    focused = cache_dir / f"{row['_candidate_id']}.focused_tree_site.pdf"
    mod.site_wrapper.create_pdf("Focused site tree", str(focused))
    site_dir = cache_dir / "csubst_sites" / "csubst.branch_id8,9"
    site_dir.mkdir(parents=True)
    (site_dir / "csubst.tsv").write_text(
        "codon_site_alignment\tOCNany2spe\n4\t2.0\n",
        encoding="utf-8",
    )
    mod.site_wrapper.create_pdf("Raw CSUBST sites summary", str(site_dir / "csubst.pdf"))
    pd.DataFrame(
        [
            {
                "output_kind": "site_table_tsv",
                "output_file": "csubst.tsv",
                "output_path": str((site_dir / "csubst.tsv").resolve()),
                "file_exists": "Y",
                "file_size_bytes": (site_dir / "csubst.tsv").stat().st_size,
            },
            {
                "output_kind": "site_summary_pdf",
                "output_file": "csubst.pdf",
                "output_path": str((site_dir / "csubst.pdf").resolve()),
                "file_exists": "Y",
                "file_size_bytes": (site_dir / "csubst.pdf").stat().st_size,
            },
            {
                "output_kind": "output_manifest",
                "output_file": "csubst.outputs.tsv",
                "output_path": str((site_dir / "csubst.outputs.tsv").resolve()),
                "file_exists": "Y",
                "file_size_bytes": 0,
            },
        ]
    ).to_csv(site_dir / "csubst.outputs.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "analysis_key": row["_analysis_key"],
                "required_input_signature": row["_required_input_signature"],
                "candidate_id": row["_candidate_id"],
            }
        ]
    ).to_csv(cache_dir / "analysis.complete.tsv", sep="\t", index=False)


def test_package_threshold_writes_self_contained_zip(monkeypatch, tmp_path):
    mod = load_module()
    summary = tmp_path / "summary.tsv"
    source = candidate_rows().iloc[[1]].copy()
    lineage_support = {
        "lineage_total": 4,
        "support_lineage_count": 3,
        "support_lineage_fraction": 0.75,
        "support_lineage_ids": "7,11,13",
    }
    for column, value in lineage_support.items():
        source[column] = value
    write_summary(summary, source)
    candidates = mod.load_threshold_candidates(
        summary_path=summary,
        minimum_support=5,
        probability_column="q_rate_enrichment_asymptotic_global",
        probability_threshold=0.05,
        max_candidates=0,
        csubst_nonsyn_recode="no",
        pdb="none",
    )
    candidates["_required_input_signature"] = "test-input-signature"
    row = candidates.iloc[0]
    annotation = mod.candidate_annotation_text(row, "q_rate_enrichment_asymptotic_global", 0.05)
    assert "Support unit count: 6" in annotation
    assert "Foreground lineage total: 4" in annotation
    assert "Support lineage count (grouped by foreground ID): 3" in annotation
    assert "Support lineage fraction: 0.75" in annotation
    assert "Support foreground lineage IDs: 7,11,13" in annotation
    cache_root = tmp_path / "cache"
    make_candidate_cache(mod, cache_root, row)
    archive = tmp_path / "candidate_sites_min_support_5.zip"

    mod.package_threshold(
        candidates=candidates,
        threshold=5,
        source_summary=summary,
        archive_path=archive,
        packages_root=tmp_path / "packages",
        cache_root=cache_root,
        probability_column="q_rate_enrichment_asymptotic_global",
        probability_threshold=0.05,
    )

    assert archive.is_file()
    with zipfile.ZipFile(archive) as zipped:
        names = zipped.namelist()
        roots = {name.split("/", 1)[0] for name in names if name}
        assert roots == {archive.stem}
        assert f"{archive.stem}/candidate_manifest.tsv" in names
        assert f"{archive.stem}/skipped_candidates.tsv" in names
        assert f"{archive.stem}/package_metadata.tsv" in names
        candidate_prefix = f"{archive.stem}/candidate_0001_{row['_candidate_id']}"
        assert f"{candidate_prefix}/candidate.tsv" in names
        assert f"{candidate_prefix}/{row['_candidate_id']}.focused_tree_site.pdf" in names
        assert f"{candidate_prefix}/{row['_candidate_id']}.report.pdf" in names
        candidate_table = pd.read_csv(zipped.open(f"{candidate_prefix}/candidate.tsv"), sep="\t")
        assert candidate_table.loc[0, "selection_min_support"] == 5
        assert candidate_table.loc[0, "selection_probability_column"] == "q_rate_enrichment_asymptotic_global"
        assert candidate_table.loc[0, "besthit_0.05"] == "protein A"
        candidate_manifest = pd.read_csv(zipped.open(f"{archive.stem}/candidate_manifest.tsv"), sep="\t")
        for column, value in lineage_support.items():
            assert candidate_table.loc[0, column] == value
            assert candidate_manifest.loc[0, column] == value
        output_manifest = pd.read_csv(
            zipped.open(f"{candidate_prefix}/csubst_sites/csubst.branch_id8,9/csubst.outputs.tsv"),
            sep="\t",
        )
        assert output_manifest["output_path"].tolist() == [
            "csubst.tsv",
            "csubst.pdf",
            "csubst.outputs.tsv",
        ]
        assert output_manifest["file_exists"].tolist() == ["Y", "Y", "Y"]
        self_row = output_manifest.loc[output_manifest["output_kind"] == "output_manifest"].iloc[0]
        zipped_manifest_info = zipped.getinfo(
            f"{candidate_prefix}/csubst_sites/csubst.branch_id8,9/csubst.outputs.tsv"
        )
        assert int(self_row["file_size_bytes"]) == zipped_manifest_info.file_size
    assert mod.archive_matches_source(archive, summary)

    with monkeypatch.context() as context:
        context.setattr(mod, "analysis_engine_signature", lambda: "changed-engine")
        context.setattr(mod, "csubst_version", lambda: "changed-csubst")
        context.setattr(mod, "runtime_dependency_versions", lambda: "changed-runtime")
        assert mod.archive_matches_source(archive, summary)

    original_summary = summary.read_text(encoding="utf-8")
    summary.write_text(original_summary + "\n", encoding="utf-8")
    assert not mod.archive_matches_source(archive, summary)
    summary.write_text(original_summary, encoding="utf-8")
    assert mod.archive_matches_source(archive, summary)

    original_archive = archive.read_bytes()
    raw_damaged = tmp_path / "raw-damaged.zip"
    with zipfile.ZipFile(archive, "r") as source_zip, zipfile.ZipFile(
        raw_damaged, "w"
    ) as target_zip:
        for member in source_zip.infolist():
            if member.filename.endswith("/csubst.tsv"):
                continue
            target_zip.writestr(member, source_zip.read(member.filename))
    raw_damaged.replace(archive)
    assert not mod.archive_matches_source(archive, summary)
    archive.write_bytes(original_archive)
    assert mod.archive_matches_source(archive, summary)

    damaged = tmp_path / "damaged.zip"
    with zipfile.ZipFile(archive, "r") as source_zip, zipfile.ZipFile(damaged, "w") as target_zip:
        for member in source_zip.infolist():
            if member.filename.endswith(f"{row['_candidate_id']}.report.pdf"):
                continue
            target_zip.writestr(member, source_zip.read(member.filename))
    damaged.replace(archive)
    assert not mod.archive_matches_source(archive, summary)


def test_package_threshold_records_candidates_skipped_for_missing_inputs(tmp_path):
    mod = load_module()
    summary = tmp_path / "summary.tsv"
    write_summary(summary, candidate_rows().iloc[[0]].copy())
    selected = mod.load_threshold_candidates(
        summary_path=summary,
        minimum_support=5,
        probability_column="q_rate_enrichment_asymptotic_global",
        probability_threshold=0.05,
        max_candidates=0,
        csubst_nonsyn_recode="no",
        pdb="none",
    )
    selected["_missing_required_inputs"] = (
        "stat_branch/OG0002_stat.branch.tsv"
    )
    selected["_required_input_signature"] = "missing-stat-branch"
    skipped = mod.skipped_candidate_frame(selected)
    eligible = selected.iloc[0:0].copy()
    archive = tmp_path / "skipped.zip"

    mod.package_threshold(
        candidates=eligible,
        threshold=5,
        source_summary=summary,
        archive_path=archive,
        packages_root=tmp_path / "packages",
        cache_root=tmp_path / "cache",
        probability_column="q_rate_enrichment_asymptotic_global",
        probability_threshold=0.05,
        skipped_candidates=skipped,
    )

    assert mod.archive_matches_source(
        archive,
        summary,
        expected_candidates=eligible,
        expected_skipped_candidates=skipped,
    )
    with zipfile.ZipFile(archive) as zipped:
        root = archive.stem
        candidate_manifest = pd.read_csv(
            zipped.open(f"{root}/candidate_manifest.tsv"), sep="\t"
        )
        skipped_manifest = pd.read_csv(
            zipped.open(f"{root}/skipped_candidates.tsv"), sep="\t"
        )
        metadata = pd.read_csv(
            zipped.open(f"{root}/package_metadata.tsv"), sep="\t"
        )
        assert candidate_manifest.empty
        assert skipped_manifest.loc[0, "orthogroup"] == "OG0002"
        assert skipped_manifest.loc[0, "reason_code"] == "missing_required_input"
        assert (
            skipped_manifest.loc[0, "missing_required_inputs"]
            == "stat_branch/OG0002_stat.branch.tsv"
        )
        assert metadata.loc[0, "selected_candidate_count"] == 1
        assert metadata.loc[0, "packaged_candidate_count"] == 0
        assert metadata.loc[0, "skipped_candidate_count"] == 1
        assert metadata.loc[0, "skipped_gene_family_count"] == 1


def test_archive_names_record_selection_and_optional_analysis_modes(tmp_path):
    mod = load_module()

    path = mod.archive_path_for_threshold(
        summary_prefix=tmp_path / "orthogroup_csubst_aa_change",
        out_dir=tmp_path,
        threshold=7,
        probability_column="q_rate_enrichment_asymptotic_global",
        probability_threshold=0.05,
        max_candidates=20,
        nonsyn_recode="dayhoff6",
        pdb="besthit",
    )

    assert path.name == (
        "orthogroup_csubst_aa_change_candidate_sites_min_support_7_"
        "q_rate_enrichment_asymptotic_global_le_0.05_top20_nonsynRecode-dayhoff6_pdb-besthit.zip"
    )


def test_validate_args_rejects_invalid_probability_threshold(tmp_path):
    mod = load_module()
    trait = tmp_path / "trait.tsv"
    trait.write_text("species\taquatic\nsp1\t1\n", encoding="utf-8")
    families = tmp_path / "orthogroup"
    families.mkdir()
    args = argparse.Namespace(
        min_support=5,
        probability_threshold=1.1,
        max_candidates=0,
        probability_column="q_rate_enrichment_asymptotic_global",
        ncpu=1,
        dir_orthogroup=str(families),
        file_trait=str(trait),
        summary_prefix=str(tmp_path / "summary"),
        out_dir=str(tmp_path / "out"),
    )

    with pytest.raises(ValueError, match="between 0 and 1"):
        mod.validate_args(args)


def test_safe_extract_zip_rejects_parent_path_traversal(tmp_path):
    mod = load_module()
    archive = tmp_path / "unsafe.zip"
    with zipfile.ZipFile(archive, "w") as zipped:
        zipped.writestr("../escaped.txt", "unsafe")

    with pytest.raises(ValueError, match="Unsafe ZIP archive"):
        mod.safe_extract_zip(archive, tmp_path / "extract", "OG0001.iqtree.anc")

    assert not (tmp_path / "escaped.txt").exists()


def test_make_csubst_manifests_portable_rejects_external_files(tmp_path):
    mod = load_module()
    candidate_dir = tmp_path / "candidate"
    site_dir = candidate_dir / "csubst_sites" / "csubst.branch_id1"
    site_dir.mkdir(parents=True)
    outside = tmp_path / "outside.tsv"
    outside.write_text("x\n", encoding="utf-8")
    pd.DataFrame(
        [
            {
                "output_kind": "site_table_tsv",
                "output_file": "../../../outside.tsv",
                "output_path": str(outside),
                "file_exists": "Y",
                "file_size_bytes": outside.stat().st_size,
            }
        ]
    ).to_csv(site_dir / "csubst.outputs.tsv", sep="\t", index=False)

    mod.make_csubst_manifests_portable(candidate_dir)

    manifest = pd.read_csv(site_dir / "csubst.outputs.tsv", sep="\t", dtype=str)
    assert pd.isna(manifest.loc[0, "output_path"])
    assert manifest.loc[0, "file_exists"] == "N"
    assert manifest.loc[0, "file_size_bytes"] == "-1"


def test_package_threshold_writes_valid_empty_zip(tmp_path):
    mod = load_module()
    summary = tmp_path / "summary.tsv"
    write_summary(summary)
    archive = tmp_path / "empty.zip"
    empty = mod.load_threshold_candidates(
        summary_path=summary,
        minimum_support=5,
        probability_column="q_rate_enrichment_asymptotic_global",
        probability_threshold=0.0,
        max_candidates=0,
        csubst_nonsyn_recode="no",
        pdb="none",
    )

    mod.package_threshold(
        candidates=empty,
        threshold=5,
        source_summary=summary,
        archive_path=archive,
        packages_root=tmp_path / "packages",
        cache_root=tmp_path / "cache",
        probability_column="q_rate_enrichment_asymptotic_global",
        probability_threshold=0.0,
    )

    assert mod.archive_matches_source(archive, summary)
    # Existing archives imply disabled lineage filtering when the additive
    # selection columns and metadata field are absent.
    with zipfile.ZipFile(archive) as zipped:
        manifest = pd.read_csv(zipped.open(f"{archive.stem}/candidate_manifest.tsv"), sep="\t")
        assert manifest.empty
        assert manifest.columns.tolist() == mod.CANDIDATE_MANIFEST_COLUMNS
    temporary = archive.with_suffix(".legacy.zip")
    with zipfile.ZipFile(archive) as original, zipfile.ZipFile(temporary, "w") as legacy:
        for member in original.namelist():
            data = original.read(member)
            drops = {
                "package_metadata.tsv": ["min_lineage_support"],
                "candidate_manifest.tsv": ["selection_min_support", "selection_min_lineage_support"],
                "skipped_candidates.tsv": ["min_lineage_support"],
            }
            if Path(member).name in drops:
                frame = pd.read_csv(io.BytesIO(data), sep="\t").drop(columns=drops[Path(member).name])
                data = frame.to_csv(sep="\t", index=False).encode()
            legacy.writestr(member, data)
    os.replace(temporary, archive)
    assert mod.archive_matches_source(archive, summary, minimum_support=5)
    assert not mod.archive_matches_source(archive, summary, minimum_support=5, minimum_lineage_support=4)


def make_run_args(tmp_path):
    family_dir = tmp_path / "orthogroup"
    family_dir.mkdir(exist_ok=True)
    trait_file = tmp_path / "species_trait.tsv"
    trait_file.write_text("species\taquatic\nsp1\t1\n", encoding="utf-8")
    return argparse.Namespace(
        summary_prefix=str(tmp_path / "orthogroup_csubst_aa_change"),
        dir_orthogroup=str(family_dir),
        file_trait=str(trait_file),
        out_dir=str(tmp_path / "out"),
        min_support=5,
        probability_column="q_rate_enrichment_asymptotic_global",
        probability_threshold=0.05,
        max_candidates=0,
        ncpu=2,
        csubst_nonsyn_recode="no",
        pdb="none",
    )


def write_run_summaries(tmp_path):
    frame = candidate_rows().copy()
    frame["q_rate_enrichment_asymptotic_global"] = 0.01
    for threshold in (5, 6, 7):
        write_summary(
            tmp_path
            / f"orthogroup_csubst_aa_change_min_support_{threshold}_summary.tsv",
            frame.loc[frame["support_unit_count"] >= threshold].copy(),
        )


def available_input_states(_dir_orthogroup, orthogroups):
    return {
        orthogroup: {
            "missing_required_inputs": [],
            "required_input_signature": f"available-{orthogroup}",
        }
        for orthogroup in orthogroups
    }


def test_run_processes_thresholds_in_descending_order(monkeypatch, tmp_path):
    mod = load_module()
    write_run_summaries(tmp_path)
    args = make_run_args(tmp_path)
    analyzed = []
    packaged = []

    def fake_ensure(candidates, **kwargs):
        analyzed.append(int(candidates["_selection_min_support"].iloc[0]))
        return []

    def fake_package(candidates, threshold, archive_path, **kwargs):
        packaged.append((threshold, candidates.shape[0]))
        Path(archive_path).parent.mkdir(parents=True, exist_ok=True)
        Path(archive_path).write_bytes(b"placeholder")

    monkeypatch.setattr(mod, "ensure_candidate_analyses", fake_ensure)
    monkeypatch.setattr(mod, "package_threshold", fake_package)
    monkeypatch.setattr(mod, "archive_matches_source", lambda *args, **kwargs: False)
    monkeypatch.setattr(mod, "inspect_required_report_inputs", available_input_states)

    manifest_path = mod.run(args)

    manifest = pd.read_csv(manifest_path, sep="\t")
    assert analyzed == [7, 6, 5]
    assert packaged == [(7, 1), (6, 2), (5, 3)]
    assert manifest["min_support"].tolist() == [7, 6, 5]
    assert manifest["status"].tolist() == ["completed", "completed", "completed"]


def test_run_skips_missing_inputs_and_packages_remaining_candidates(
    monkeypatch, tmp_path
):
    mod = load_module()
    write_run_summaries(tmp_path)
    args = make_run_args(tmp_path)
    analyzed = []
    packaged = []

    def input_states(_dir_orthogroup, orthogroups):
        return {
            orthogroup: {
                "missing_required_inputs": (
                    ["stat_branch/OG0002_stat.branch.tsv"]
                    if orthogroup == "OG0002"
                    else []
                ),
                "required_input_signature": f"state-{orthogroup}",
            }
            for orthogroup in orthogroups
        }

    def fake_ensure(candidates, **kwargs):
        analyzed.extend(candidates["orthogroup"].tolist())
        return []

    def fake_package(candidates, threshold, archive_path, skipped_candidates, **kwargs):
        packaged.append(
            (
                threshold,
                candidates["orthogroup"].tolist(),
                skipped_candidates["orthogroup"].tolist(),
            )
        )
        Path(archive_path).parent.mkdir(parents=True, exist_ok=True)
        Path(archive_path).write_bytes(b"placeholder")

    monkeypatch.setattr(mod, "inspect_required_report_inputs", input_states)
    monkeypatch.setattr(mod, "ensure_candidate_analyses", fake_ensure)
    monkeypatch.setattr(mod, "package_threshold", fake_package)
    monkeypatch.setattr(mod, "archive_matches_source", lambda *args, **kwargs: False)

    manifest_path = mod.run(args)

    manifest = pd.read_csv(manifest_path, sep="\t")
    assert packaged == [
        (7, ["OG0003"], []),
        (6, ["OG0001", "OG0003"], []),
        (5, ["OG0001", "OG0003"], ["OG0002"]),
    ]
    assert "OG0002" not in analyzed
    assert manifest["candidate_count"].tolist() == [1, 2, 3]
    assert manifest["packaged_candidate_count"].tolist() == [1, 2, 2]
    assert manifest["skipped_candidate_count"].tolist() == [0, 0, 1]
    assert manifest["skipped_gene_family_count"].tolist() == [0, 0, 1]
    assert manifest["status"].tolist() == [
        "completed",
        "completed",
        "completed_with_skips",
    ]
    skipped_path = tmp_path / "out" / manifest.loc[0, "skipped_candidates_tsv"]
    skipped = pd.read_csv(skipped_path, sep="\t")
    assert skipped["orthogroup"].tolist() == ["OG0002"]
    assert skipped["reason_code"].tolist() == ["missing_required_input"]
    assert skipped["missing_required_inputs"].tolist() == [
        "stat_branch/OG0002_stat.branch.tsv"
    ]


def test_run_records_batch_setup_failure_in_manifest(monkeypatch, tmp_path):
    mod = load_module()
    write_run_summaries(tmp_path)
    args = make_run_args(tmp_path)

    def fail_analysis(**kwargs):
        raise RuntimeError("materialization failed")

    monkeypatch.setattr(mod, "ensure_candidate_analyses", fail_analysis)
    monkeypatch.setattr(mod, "archive_matches_source", lambda *args, **kwargs: False)
    monkeypatch.setattr(mod, "inspect_required_report_inputs", available_input_states)

    with pytest.raises(RuntimeError, match="packaging failed"):
        mod.run(args)

    manifests = list((tmp_path / "out").glob("*_manifest.tsv"))
    assert len(manifests) == 1
    manifest = pd.read_csv(manifests[0], sep="\t")
    assert manifest.loc[0, "min_support"] == 7
    assert manifest.loc[0, "status"] == "failed"
    assert manifest.loc[0, "error"] == "materialization failed"
    assert manifest.loc[1:, "status"].tolist() == ["pending", "pending"]


def test_run_keeps_candidate_analysis_errors_as_hard_failures(monkeypatch, tmp_path):
    mod = load_module()
    write_run_summaries(tmp_path)
    args = make_run_args(tmp_path)
    packaged = []

    def failed_analysis(candidates, **kwargs):
        return [
            {
                "candidate_id": candidates.iloc[0]["_candidate_id"],
                "status": "failed",
                "error": "CSUBST/stat_branch branch identity mismatch",
            }
        ]

    monkeypatch.setattr(mod, "inspect_required_report_inputs", available_input_states)
    monkeypatch.setattr(mod, "ensure_candidate_analyses", failed_analysis)
    monkeypatch.setattr(
        mod, "package_threshold", lambda **kwargs: packaged.append(kwargs)
    )
    monkeypatch.setattr(mod, "archive_matches_source", lambda *args, **kwargs: False)

    with pytest.raises(RuntimeError, match="packaging failed"):
        mod.run(args)

    assert packaged == []
    manifest_path = next((tmp_path / "out").glob("*_manifest.tsv"))
    manifest = pd.read_csv(manifest_path, sep="\t")
    assert manifest.loc[0, "status"] == "failed"
    assert "branch identity mismatch" in manifest.loc[0, "error"]


def test_run_writes_empty_manifest_when_no_threshold_reaches_minimum(tmp_path):
    mod = load_module()
    for threshold in (2, 3, 4):
        write_summary(
            tmp_path
            / f"orthogroup_csubst_aa_change_min_support_{threshold}_summary.tsv"
        )
    args = make_run_args(tmp_path)

    manifest_path = mod.run(args)

    manifest = pd.read_csv(manifest_path, sep="\t")
    assert manifest.empty
    assert manifest.columns.tolist() == mod.ARCHIVE_MANIFEST_COLUMNS
    assert not list((tmp_path / "out").glob("*.zip"))
    assert not any((tmp_path / "out").glob(".*.work"))


def test_archived_family_inputs_are_materialized_once_per_orthogroup(
    monkeypatch, tmp_path
):
    mod = load_module()
    family_dir = tmp_path / "orthogroup"
    (family_dir / ".gg_store").mkdir(parents=True)
    materialized_dir = tmp_path / "materialized" / "OG0001"
    materialized_dir.mkdir(parents=True)
    events = []

    class FakeMaterializationDirectory:
        def __init__(self, parent, orthogroup):
            events.append(("lock", str(parent), orthogroup))
            self.name = str(materialized_dir)

        def cleanup(self):
            events.append(("cleanup",))

    def fake_materialize(dir_og, og, destination_root):
        events.append(("materialize", dir_og, og, destination_root))
        return []

    def fake_analyze(record, effective_dir_orthogroup, **kwargs):
        events.append(("analyze", record["_candidate_id"], effective_dir_orthogroup))
        return {"candidate_id": record["_candidate_id"], "status": "completed", "error": ""}

    monkeypatch.setattr(
        mod.site_wrapper,
        "LockedMaterializationDirectory",
        FakeMaterializationDirectory,
    )
    monkeypatch.setattr(mod.site_wrapper, "materialize_csubst_site_inputs", fake_materialize)
    monkeypatch.setattr(mod, "analyze_candidate", fake_analyze)
    records = [
        {"_candidate_id": "candidate1"},
        {"_candidate_id": "candidate2"},
    ]

    results = mod.analyze_orthogroup_batch(
        orthogroup="OG0001",
        records=records,
        cache_root=tmp_path / "cache",
        dir_orthogroup=str(family_dir),
        materialization_parent=tmp_path / "materialization_parent",
        trait_color_paths={},
        nonsyn_recode="no",
        pdb="none",
    )

    assert [result["candidate_id"] for result in results] == ["candidate1", "candidate2"]
    assert sum(event[0] == "materialize" for event in events) == 1
    assert sum(event[0] == "analyze" for event in events) == 2
    assert events[-1] == ("cleanup",)


def test_cli_end_to_end_reuses_analysis_across_threshold_zips(tmp_path):
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    csubst_log = tmp_path / "csubst.log"
    rscript_log = tmp_path / "rscript.log"
    fake_csubst = fake_bin / "csubst"
    fake_csubst.write_text(
        """#!/usr/bin/env python3
import os
import sys
from pathlib import Path
from reportlab.pdfgen import canvas

args = sys.argv[1:]
def value(flag):
    return args[args.index(flag) + 1]

branch_ids = value("--branch_id")
outdir = Path(value("--outdir")) / (value("--output_prefix") + ".branch_id" + branch_ids)
outdir.mkdir(parents=True, exist_ok=True)
site_tsv = outdir / "csubst.tsv"
site_pdf = outdir / "csubst.pdf"
site_tsv.write_text("codon_site_alignment\\tOCNany2spe\\n1\\t0.0\\n2\\t2.0\\n3\\t0.0\\n", encoding="utf-8")
pdf = canvas.Canvas(str(site_pdf))
pdf.drawString(72, 720, "Fake csubst sites summary")
pdf.save()
manifest = outdir / "csubst.outputs.tsv"
manifest.write_text(
    "output_kind\\toutput_file\\toutput_path\\tfile_exists\\tfile_size_bytes\\n"
    + f"site_table_tsv\\tcsubst.tsv\\t{site_tsv.resolve()}\\tY\\t{site_tsv.stat().st_size}\\n"
    + f"site_summary_pdf\\tcsubst.pdf\\t{site_pdf.resolve()}\\tY\\t{site_pdf.stat().st_size}\\n"
    + f"output_manifest\\tcsubst.outputs.tsv\\t{manifest.resolve()}\\tY\\t0\\n",
    encoding="utf-8",
)
with open(os.environ["FAKE_CSUBST_LOG"], "a", encoding="utf-8") as handle:
    handle.write(" ".join(args) + "\\n")
""",
        encoding="utf-8",
    )
    fake_rscript = fake_bin / "Rscript"
    fake_rscript.write_text(
        """#!/usr/bin/env python3
import os
import sys
from reportlab.pdfgen import canvas

pdf = canvas.Canvas("stat_branch2tree_plot.pdf")
pdf.drawString(72, 720, "Fake focused tree")
pdf.save()
with open(os.environ["FAKE_RSCRIPT_LOG"], "a", encoding="utf-8") as handle:
    handle.write(" ".join(sys.argv[1:]) + "\\n")
""",
        encoding="utf-8",
    )
    fake_csubst.chmod(0o755)
    fake_rscript.chmod(0o755)

    family_dir = tmp_path / "orthogroup"
    for subdir in ("iqtree_anc", "clipkit", "stat_branch", "rpsblast"):
        (family_dir / subdir).mkdir(parents=True)
    for og in ("OG0001", "OG0002"):
        iqtree_zip = family_dir / "iqtree_anc" / f"{og}_iqtree.anc.zip"
        with zipfile.ZipFile(iqtree_zip, "w") as archive:
            for filename, content in {
                "csubst.fasta": ">sp1\nGCTGTTGAT\n>sp2\nGCTGTTGAT\n>sp3\nGCTGTTGAT\n",
                "csubst.nwk": "(sp1:0.1,(sp2:0.1,sp3:0.1)RootedBC:0.1)RootedRoot;\n",
                "csubst.treefile": "(sp1:0.1,(sp2:0.1,sp3:0.1)IqtreeBC:0.1)IqtreeRoot;\n",
                "csubst.state": "placeholder\n",
                "csubst.rate": "placeholder\n",
                "csubst.iqtree": "placeholder\n",
                "csubst.log": f"Converting to codon sequences with genetic code {1 if og == 'OG0001' else 2} ...\n",
            }.items():
                archive.writestr(f"{og}.iqtree.anc/{filename}", content)
        (family_dir / "clipkit" / f"{og}_cds.clipkit.fa").write_text(
            ">sp1\nGCTGTTGAT\n>sp2\nGCTGTTGAT\n>sp3\nGCTGTTGAT\n",
            encoding="utf-8",
        )
        (family_dir / "stat_branch" / f"{og}_stat.branch.tsv").write_text(
            "branch_id\tnode_name\tgene_labels\n"
            "0\tsp1\tsp1\n"
            "1\tsp2\tsp2\n"
            "2\tsp3\tsp3\n"
            "3\tN1\tsp2; sp3\n"
            "4\tRoot\tsp1; sp2; sp3\n",
            encoding="utf-8",
        )
        (family_dir / "rpsblast" / f"{og}_rpsblast.tsv").write_text(
            "query\nsp1\n",
            encoding="utf-8",
        )

    summary_row = pd.DataFrame(
        [
            {
                "orthogroup": "OG0001",
                "trait": "aquatic",
                "state_change": "2V",
                "codon_site_alignment": 2,
                "support_unit_count": 6,
                "support_unit_ids": "1,2,3,4,5,6",
                "support_branch_ids": "0,1",
                "p_rate_enrichment_asymptotic": 0.001,
                "scan_calibration_status": "conditional_assignment",
                "q_rate_enrichment_asymptotic_global": 0.01,
                "besthit_0.05": "annotated protein",
            },
            {
                "orthogroup": "OG0002",
                "trait": "aquatic",
                "state_change": "3D",
                "codon_site_alignment": 3,
                "support_unit_count": 6,
                "support_unit_ids": "1,2,3,4,5,6",
                "support_branch_ids": "2,3",
                "p_rate_enrichment_asymptotic": 0.002,
                "scan_calibration_status": "conditional_assignment",
                "q_rate_enrichment_asymptotic_global": 0.02,
                "besthit_0.05": "second annotated protein",
            },
        ]
    )
    prefix = tmp_path / "orthogroup_csubst_aa_change"
    for threshold in (5, 6):
        summary_row.to_csv(
            tmp_path / f"{prefix.name}_min_support_{threshold}_summary.tsv",
            sep="\t",
            index=False,
        )
    trait_file = tmp_path / "species_trait.tsv"
    trait_file.write_text("species\taquatic\nsp1\t1\nsp2\t0\nsp3\t0\n", encoding="utf-8")
    output_dir = tmp_path / "out"
    env = os.environ.copy()
    env["PATH"] = str(fake_bin) + os.pathsep + env["PATH"]
    env["FAKE_CSUBST_LOG"] = str(csubst_log)
    env["FAKE_RSCRIPT_LOG"] = str(rscript_log)

    command = [
        sys.executable,
        str(SCRIPT_PATH),
        "--summary_prefix",
        str(prefix),
        "--dir_orthogroup",
        str(family_dir),
        "--file_trait",
        str(trait_file),
        "--out_dir",
        str(output_dir),
        "--min_support",
        "5",
        "--ncpu",
        "2",
        "--pdb",
        "none",
    ]
    processes = [
        subprocess.Popen(
            command,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        for _ in range(2)
    ]
    results = []
    for process in processes:
        stdout, stderr = process.communicate(timeout=60)
        results.append((process.returncode, stdout, stderr))

    assert all(returncode == 0 for returncode, _, _ in results), results
    combined_stdout = "\n".join(stdout for _, stdout, _ in results)
    assert "candidates packaged" in combined_stdout
    assert "existing ZIP retained" in combined_stdout
    assert len(csubst_log.read_text(encoding="utf-8").splitlines()) == 2
    assert "--genetic_code 1" in csubst_log.read_text(encoding="utf-8")
    assert "--genetic_code 2" in csubst_log.read_text(encoding="utf-8")
    assert "--pdb" not in csubst_log.read_text(encoding="utf-8")
    rscript_lines = rscript_log.read_text(encoding="utf-8").splitlines()
    assert len(rscript_lines) == 2
    assert any("amino_acid_site,1,2," in line for line in rscript_lines)
    assert any("amino_acid_site,2,3," in line for line in rscript_lines)
    run_manifest = pd.read_csv(next(output_dir.glob("*_manifest.tsv")), sep="\t")
    assert run_manifest["min_support"].tolist() == [6, 5]
    assert run_manifest["status"].tolist() == ["existing", "existing"]
    archives = sorted(output_dir.glob("*.zip"))
    assert len(archives) == 2
    from pypdf import PdfReader

    for archive_path in archives:
        assert archive_path.is_file()
        with zipfile.ZipFile(archive_path) as archive:
            assert archive.testzip() is None
            candidate_manifest_path = f"{archive_path.stem}/candidate_manifest.tsv"
            candidate_manifest = pd.read_csv(archive.open(candidate_manifest_path), sep="\t")
            assert candidate_manifest.shape[0] == 2
            for _, candidate in candidate_manifest.iterrows():
                report_path = f"{archive_path.stem}/{candidate['report_pdf']}"
                assert len(PdfReader(io.BytesIO(archive.read(report_path))).pages) == 3
                raw_manifest_path = (
                    f"{archive_path.stem}/{candidate['csubst_sites_dir']}"
                    f"/csubst.branch_id{candidate['support_branch_ids']}/csubst.outputs.tsv"
                )
                raw_manifest = pd.read_csv(archive.open(raw_manifest_path), sep="\t")
                assert raw_manifest["output_path"].tolist() == [
                    "csubst.tsv",
                    "csubst.pdf",
                    "csubst.outputs.tsv",
                ]
                self_row = raw_manifest.loc[
                    raw_manifest["output_kind"] == "output_manifest"
                ].iloc[0]
                assert int(self_row["file_size_bytes"]) == archive.getinfo(
                    raw_manifest_path
                ).file_size
    assert not any(output_dir.glob(".*.work"))
    assert not any(output_dir.glob(".*.lock*"))


@pytest.mark.parametrize("value", [-0.01, 1.01, float("inf"), "invalid"])
def test_candidate_selection_rejects_invalid_source_probabilities(tmp_path, value):
    mod = load_module()
    frame = candidate_rows().astype({"q_rate_enrichment_asymptotic_global": object})
    frame.loc[0, "q_rate_enrichment_asymptotic_global"] = value
    source = tmp_path / "summary.tsv"
    write_summary(source, frame)
    with pytest.raises(ValueError, match="invalid probabilities"):
        mod.load_threshold_candidates(source, 5, "q_rate_enrichment_asymptotic_global", 0.05, 0, "no", "none")


def test_missing_fdr_never_uses_asymptotic_values(tmp_path):
    mod = load_module()
    frame = candidate_rows()
    frame["q_rate_enrichment_asymptotic_global"] = float("nan")
    frame["scan_calibration_status"] = "unavailable_failed_trials"
    source = tmp_path / "summary.tsv"
    write_summary(source, frame)
    selected = mod.load_threshold_candidates(source, 5, "q_rate_enrichment_asymptotic_global", 0.05, 0, "no", "none")
    assert selected.empty
    frame = frame.drop(columns="q_rate_enrichment_asymptotic_global")
    write_summary(source, frame)
    with pytest.raises(ValueError, match="missing required candidate column"):
        mod.load_threshold_candidates(source, 5, "q_rate_enrichment_asymptotic_global", 0.05, 0, "no", "none")




def test_empirical_probability_columns_are_not_accepted(tmp_path):
    mod = load_module()
    source = tmp_path / "summary.tsv"
    write_summary(source)
    with pytest.raises(ValueError, match="Unsupported scan probability column"):
        mod.load_threshold_candidates(source, 5, "p_rate_enrichment_empirical_maxT", 0.05, 0, "no", "none")


def test_candidate_alphabet_must_match_requested_3di_mode():
    mod = load_module()
    frame = candidate_rows()
    with pytest.raises(ValueError, match="source scan alphabet"):
        mod.assign_candidate_ids(frame, "3di20", "none")
    frame["nonsyn_recode"] = "no"
    with pytest.raises(ValueError, match="differs from requested"):
        mod.assign_candidate_ids(frame, "3di20", "none")


def test_missing_3di_source_alphabet_is_not_treated_as_a_match():
    mod = load_module()
    frame = candidate_rows()
    frame["nonsyn_recode"] = None
    with pytest.raises(ValueError, match="differs from requested"):
        mod.assign_candidate_ids(frame, "3di20", "none")


def test_filtered_bh_pools_traits_and_matches_before_cutoff_and_cap(tmp_path):
    mod = load_module()
    source = candidate_rows()
    source['trait'] = ['aquatic', 'terrestrial', 'aquatic']
    source['scan_match'] = ['m1', 'm2', 'm2']
    source['support_lineage_count'] = [4, 4, 4]
    source['p_rate_enrichment_asymptotic'] = [0.01, 0.03, float('nan')]
    source[mod.DEFAULT_PROBABILITY_COLUMN] = 0.00001  # Stale view q must be recomputed.
    summary = tmp_path / 'summary.tsv'
    write_summary(summary, source)
    selected = mod.load_threshold_candidates(summary, 5, mod.DEFAULT_PROBABILITY_COLUMN,
                                             0.025, 1, 'no', 'none', 4)
    assert selected['orthogroup'].tolist() == ['OG0002']
    assert selected.iloc[0][mod.DEFAULT_PROBABILITY_COLUMN] == 0.02
    metadata = selected.attrs['support_bh_metadata']
    assert metadata['support_bh_candidate_count'] == 3
    assert metadata['support_bh_test_count'] == 2
    assert metadata['support_bh_undefined_count'] == 1
    assert selected.iloc[0]['q_rate_enrichment_asymptotic_global'] == 0.02


def test_filtered_bh_empty_selection_keeps_full_family_metadata(tmp_path):
    mod = load_module()
    source = candidate_rows()
    source['p_rate_enrichment_asymptotic'] = [0.04, 0.1, 0.6]
    source['support_lineage_count'] = [4, 3, 3]
    summary = tmp_path / 'summary.tsv'
    write_summary(summary, source)
    selected = mod.load_threshold_candidates(summary, 5, mod.DEFAULT_PROBABILITY_COLUMN,
                                             0.05, 1, 'no', 'none', 0)
    assert selected.empty  # q for the smallest P is .12, even with cap 1.
    assert selected.attrs['support_bh_metadata']['support_bh_test_count'] == 3
    selected = mod.load_threshold_candidates(summary, 5, mod.DEFAULT_PROBABILITY_COLUMN,
                                             0.05, 1, 'no', 'none', 4)
    assert selected['orthogroup'].tolist() == ['OG0002']
    assert selected.iloc[0][mod.DEFAULT_PROBABILITY_COLUMN] == 0.04


@pytest.mark.parametrize('value', [-0.01, 1.01, float('inf'), 'invalid'])
def test_filtered_bh_rejects_invalid_p_even_outside_support_pool(tmp_path, value):
    mod = load_module()
    source = candidate_rows().astype({'p_rate_enrichment_asymptotic': object})
    source.loc[0, 'p_rate_enrichment_asymptotic'] = value
    source['support_lineage_count'] = [3, 4, 4]
    summary = tmp_path / 'summary.tsv'
    write_summary(summary, source)
    with pytest.raises(ValueError, match='invalid probabilities'):
        mod.load_threshold_candidates(summary, 6, mod.DEFAULT_PROBABILITY_COLUMN,
                                      0.05, 0, 'no', 'none', 4, filter_unit_support=True)


def test_filtered_bh_requires_p_and_does_not_use_missing_or_global_q(tmp_path):
    mod = load_module()
    source = candidate_rows()
    source['p_rate_enrichment_asymptotic'] = float('nan')
    summary = tmp_path / 'summary.tsv'
    write_summary(summary, source)
    selected = mod.load_threshold_candidates(summary, 5, mod.DEFAULT_PROBABILITY_COLUMN, 0.05, 0, 'no', 'none')
    assert selected.empty
    assert selected.attrs['support_bh_metadata']['support_bh_undefined_count'] == 3
    source = source.drop(columns='p_rate_enrichment_asymptotic')
    write_summary(summary, source)
    with pytest.raises(ValueError, match='missing required candidate column'):
        mod.load_threshold_candidates(summary, 5, mod.DEFAULT_PROBABILITY_COLUMN, 0.05, 0, 'no', 'none')


def test_filtered_bh_family_survives_cap_missing_inputs_empty_zip_and_archive_reuse(tmp_path, monkeypatch):
    mod = load_module()
    args = make_run_args(tmp_path)
    args.min_support = 0
    args.min_lineage_support = 4
    args.probability_column = mod.DEFAULT_PROBABILITY_COLUMN
    args.max_candidates = 1
    source = candidate_rows()
    source['support_lineage_count'] = 4
    write_summary(Path(f'{args.summary_prefix}_all_candidates_summary.tsv'), source)
    def missing_inputs(_directory, orthogroups):
        return {orthogroup: {'missing_required_inputs': ['iqtree_anc/missing'],
                            'required_input_signature': f'missing-{orthogroup}'} for orthogroup in orthogroups}
    monkeypatch.setattr(mod, 'inspect_required_report_inputs', missing_inputs)
    monkeypatch.setattr(mod, 'ensure_candidate_analyses', lambda **kwargs: [])
    manifest_path = mod.run(args)
    manifest = pd.read_csv(manifest_path, sep='\t').iloc[0]
    assert manifest['support_bh_test_count'] == manifest['support_bh_candidate_count'] == 3
    assert manifest['candidate_count'] == manifest['skipped_candidate_count'] == 1
    assert manifest['packaged_candidate_count'] == 0
    archive_path = Path(args.out_dir) / manifest['archive_zip']
    with zipfile.ZipFile(archive_path) as archive:
        members = {name: archive.read(name) for name in archive.namelist()}
        metadata_name = f'{archive_path.stem}/package_metadata.tsv'
        metadata = pd.read_csv(io.BytesIO(members[metadata_name]), sep='\t')
        assert metadata.loc[0, 'support_bh_test_count'] == 3
        assert metadata.loc[0, 'packaged_candidate_count'] == 0
    assert pd.read_csv(mod.run(args), sep='\t').loc[0, 'status'] == 'existing_with_skips'
    for column, value in [('probability_policy', 'old_policy'), ('support_bh_test_count', 2)]:
        corrupted = metadata.copy()
        corrupted.loc[0, column] = value
        with zipfile.ZipFile(archive_path, 'w') as archive:
            for name, content in members.items():
                archive.writestr(name, corrupted.to_csv(sep='\t', index=False) if name == metadata_name else content)
        assert not mod.archive_matches_source(archive_path, Path(f'{args.summary_prefix}_all_candidates_summary.tsv'),
                                              minimum_support=0, minimum_lineage_support=4)


def test_filtered_bh_accepts_nullable_sqlite_testable_flags_and_rejects_untestable_p(tmp_path):
    mod = load_module()
    assert mod.rate_testable_mask(pd.Series([1, pd.NA, 0], dtype="Int64")).tolist() == [True, False, False]
    source = candidate_rows()
    # SQLite integer flags with NULL become float64 on a pandas read.
    source['scan_rate_testable'] = [1.0, 1.0, float('nan')]
    source['p_rate_enrichment_asymptotic'] = [0.01, 0.03, float('nan')]
    summary = tmp_path / 'summary.tsv'
    write_summary(summary, source)
    selected = mod.load_threshold_candidates(summary, 5, mod.DEFAULT_PROBABILITY_COLUMN, 1.0, 0, 'no', 'none')
    assert selected[mod.DEFAULT_PROBABILITY_COLUMN].tolist() == pytest.approx([0.02, 0.03])
    source.loc[0, 'scan_rate_testable'] = 0.0
    write_summary(summary, source)
    with pytest.raises(ValueError, match='untestable'):
        mod.load_threshold_candidates(summary, 5, mod.DEFAULT_PROBABILITY_COLUMN, 1.0, 0, 'no', 'none')


def test_filtered_bh_run_uses_full_family_instead_of_stale_threshold_views(tmp_path, monkeypatch):
    mod = load_module()
    args = make_run_args(tmp_path)
    args.probability_column = mod.DEFAULT_PROBABILITY_COLUMN
    args.min_lineage_support = 4
    source = candidate_rows().iloc[:2].copy()
    source['support_unit_count'] = 5
    source['support_lineage_count'] = 4
    source['p_rate_enrichment_asymptotic'] = [0.03, 0.9]
    write_summary(Path(f'{args.summary_prefix}_all_candidates_summary.tsv'), source)
    # The stale table omits a nonsignificant hypothesis from the BH family.
    write_summary(Path(f'{args.summary_prefix}_min_support_5_summary.tsv'), source.iloc[:1])
    observed = {}
    def package(candidates, threshold, **kwargs):
        observed[threshold] = (len(candidates), candidates.attrs['support_bh_metadata']['support_bh_test_count'])
    monkeypatch.setattr(mod, 'package_threshold', package)
    monkeypatch.setattr(mod, 'archive_matches_source', lambda *args, **kwargs: False)
    monkeypatch.setattr(mod, 'inspect_required_report_inputs', available_input_states)
    monkeypatch.setattr(mod, 'ensure_candidate_analyses', lambda **kwargs: [])
    mod.run(args)
    assert observed == {5: (0, 2)}  # q=.06: both P values must enter BH.


@pytest.mark.parametrize('member_name,column,value', [
    ('package_metadata.tsv', 'probability_column', 'q_rate_enrichment_asymptotic_global'),
    ('package_metadata.tsv', 'probability_threshold', 0.01),
    ('candidate_manifest.tsv', 'probability_value', 0.003),
])
def test_filtered_bh_archive_does_not_reuse_a_different_probability_policy(tmp_path, member_name, column, value):
    mod = load_module()
    source = candidate_rows().iloc[[1]].copy()
    source['support_lineage_count'] = 4
    summary = tmp_path / 'source.tsv'
    write_summary(summary, source)
    candidates = mod.load_threshold_candidates(summary, 5, mod.DEFAULT_PROBABILITY_COLUMN, 0.05, 0, 'no', 'none', 4)
    candidates['_required_input_signature'] = 'input-signature'
    row = candidates.iloc[0]
    cache = tmp_path / 'cache'
    make_candidate_cache(mod, cache, row)
    archive = mod.archive_path_for_threshold(tmp_path / 'scan', tmp_path, 5, mod.DEFAULT_PROBABILITY_COLUMN,
                                             0.05, 0, 'no', 'none', 4)
    mod.package_threshold(candidates, 5, summary, archive, tmp_path / 'packages', cache,
                          mod.DEFAULT_PROBABILITY_COLUMN, 0.05, minimum_lineage_support=4)
    with zipfile.ZipFile(archive) as zipped:
        members = {name: zipped.read(name) for name in zipped.namelist()}
    metadata_name = f'{archive.stem}/{member_name}'
    metadata = pd.read_csv(io.BytesIO(members[metadata_name]), sep='\t')
    metadata.loc[0, column] = value
    with zipfile.ZipFile(archive, 'w') as zipped:
        for name, content in members.items():
            zipped.writestr(name, metadata.to_csv(sep='\t', index=False) if name == metadata_name else content)
    assert not mod.archive_matches_source(archive, summary, expected_candidates=candidates,
                                          minimum_support=5, minimum_lineage_support=4,
                                          probability_column=mod.DEFAULT_PROBABILITY_COLUMN, probability_threshold=0.05)


def test_filtered_bh_legacy_fallback_uses_broadest_pool_for_every_threshold(tmp_path):
    mod = load_module()
    prefix = tmp_path / 'scan'
    source = candidate_rows()
    source['support_lineage_count'] = [4, 4, 4]
    source['p_rate_enrichment_asymptotic'] = [0.01, 0.03, 0.9]
    broadest = Path(f'{prefix}_min_support_2_summary.tsv')
    write_summary(broadest, source)
    write_summary(Path(f'{prefix}_min_support_6_summary.tsv'), source.iloc[[1]])
    discovered = mod.discover_summary_tables(prefix, 5, 4, support_bh_family=True)
    assert discovered == {7: broadest, 6: broadest, 5: broadest}
    candidates = mod.load_threshold_candidates(discovered[6], 6, mod.DEFAULT_PROBABILITY_COLUMN,
                                               0.05, 0, 'no', 'none', 4, filter_unit_support=True)
    assert candidates.empty  # Both .03 and .9 enter BH; q=.06.
    assert candidates.attrs['support_bh_metadata']['support_bh_test_count'] == 2


def test_filtered_bh_source_preflight_checks_p_when_unit_bound_selects_no_threshold(tmp_path):
    mod = load_module()
    prefix = tmp_path / 'scan'
    source = candidate_rows().astype({'p_rate_enrichment_asymptotic': object})
    source.loc[0, 'p_rate_enrichment_asymptotic'] = 'invalid'
    write_summary(Path(f'{prefix}_all_candidates_summary.tsv'), source)
    with pytest.raises(ValueError, match='invalid probabilities'):
        mod.discover_summary_tables(prefix, 99, 0, support_bh_family=True)
