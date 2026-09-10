"""Acquisition/geometry contracts, including negative controls for bias."""

import csv
import json
import math
import sys
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pandas
import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
import gbif_observations as gbif  # noqa: E402
from species_trait_contract import (  # noqa: E402
    metadata_path,
    select_foreground_traits,
    select_species_traits,
    write_trait_bundle,
)


def row(key=1, lat=10, lon=20, **kwargs):
    return gbif.normalize_record({"key": key, "speciesKey": 100, "decimalLatitude": lat,
                                  "decimalLongitude": lon, "occurrenceStatus": "PRESENT", "issues": [],
                                  "countryCode": "JP", "basisOfRecord": "HUMAN_OBSERVATION", **kwargs})


def search(rows, config=None, change=None):
    config = gbif.effective_config({"gbif_page_size": 2, **(config or {})})
    calls = []

    def fetch_json(url, timeout):
        params = parse_qs(urlparse(url).query)
        offset, limit = int(params["offset"][0]), int(params["limit"][0])
        calls.append((offset, limit))
        payload = {"count": len(rows), "limit": limit, "offset": offset,
                   "endOfRecords": offset + limit >= len(rows), "results": rows[offset:offset + limit]}
        return change(payload, len(calls)) if change else payload

    result = gbif.fetch_gbif_occurrence_records("https://example.test/v1/", "100", config, 1, fetch_json)
    return result, calls


@pytest.mark.parametrize("count,cap,status", [(0, 4, "complete_search"), (3, 4, "complete_search"),
                                            (4, 4, "complete_search"), (5, 4, "capped_partial")])
def test_search_boundaries(count, cap, status):
    result, calls = search([row(key=index + 1) for index in range(count)], {"gbif_max_occurrences_per_species": cap})
    assert result["status"] == status
    assert result["reported_count"] == count
    assert result["raw_fetched"] == min(count, cap)
    assert calls[0] == (0, 0)
    assert sum(limit for _offset, limit in calls) == min(count, cap)
    assert result["started_utc"] <= result["ended_utc"]


@pytest.mark.parametrize("mutation,status,reason", [
    (lambda p: {**p, "results": []}, "interrupted", "short_or_empty_page"),
    (lambda p: {**p, "results": p["results"][:1]}, "interrupted", "short_or_empty_page"),
    (lambda p: {**p, "count": 9}, "inconsistent_search", "count_changed_or_missing"),
    (lambda p: {**p, "endOfRecords": True}, "inconsistent_search", "end_of_records_mismatch"),
    (lambda p: {**p, "offset": 99}, "invalid_response", "page_position_mismatch"),
    (lambda p: {**p, "results": {}}, "invalid_response", "invalid_page"),
    (lambda p: {**p, "results": [None]}, "invalid_response", "invalid_page_rows"),
    (lambda p: {**p, "results": [p["results"][0], p["results"][0]]}, "inconsistent_search", "duplicate_or_missing_record_id"),
    (lambda p: {**p, "results": [{**r, "gbifID": None} for r in p["results"]]}, "inconsistent_search", "duplicate_or_missing_record_id"),
])
def test_search_never_marks_incomplete_pages_complete(mutation, status, reason):
    result, _calls = search([row(key=index + 1) for index in range(6)], change=lambda p, n: mutation(p) if n == 3 else p)
    assert result["status"] == status
    assert result["termination_reason"] == reason
    assert result["raw_fetched"] >= 2
    _, quality, _ = gbif.summarize_records(result, gbif.effective_config({}))
    assert quality["analysis_eligible"] is False


def test_search_request_failure_retains_partial_audit():
    def fail(payload, call):
        if call == 3:
            raise TimeoutError("request failed")
        return payload

    result, _ = search([row(key=index + 1) for index in range(6)], change=fail)
    assert result["status"] == "interrupted"
    assert result["raw_fetched"] == 2
    assert len(result["rows"]) == 2
    assert result["error_type"] == "TimeoutError"


@pytest.mark.parametrize("payload", [[], {}, {"count": -1}, {"count": True}, {"count": "9"}])
def test_invalid_count_is_not_zero(payload):
    result, _ = search([], change=lambda p, n: payload)
    assert result["status"] == "invalid_response"
    assert result["reported_count"] is None


def test_identity_dedup_is_distinct_from_location_and_filter_counts():
    records = [row(1), row(1), row(2), row(3, lat=999)]
    acquisition = {"rows": records, "status": "complete_download", "raw_fetched": 4}
    metrics, quality, retained = gbif.summarize_records(acquisition, gbif.effective_config({}))
    assert quality["retained_count"] == 2
    assert quality["unique_location_count"] == 1
    assert quality["unique_gbif_ids"] == 3
    assert quality["exclusions"] == {"duplicate_gbif_id": 1, "invalid_coordinates": 1}
    assert quality["raw_fetched"] == quality["retained_count"] + quality["excluded_count"]
    assert quality["analysis_eligible"] is False
    assert len(retained) == 2
    assert metrics["gbif_observed_occupied_grid_area_km2"] == gbif.occupied_grid_area_km2([(10, 20, "JP")], 1)


def test_zero_coordinates_parse_as_numbers_but_null_is_not_a_coordinate():
    config = gbif.effective_config({})
    assert gbif.gbif_parse_float(0) == 0
    assert gbif.gbif_row_to_point(row(lat=0, lon=20), config) == (0, 20, "JP")
    assert gbif.record_exclusion_reason(row(lat=0, lon=0), config) == "zero_coordinate"
    for invalid in [None, "", "NaN", "inf", float("inf")]:
        assert gbif.record_exclusion_reason(row(lat=invalid), config) == "invalid_coordinates"


def test_known_centroid_filter_excludes_nearby_not_distant_records():
    config = gbif.effective_config({"gbif_min_distance_from_known_centroid_m": 2000})
    for value, expected in [(1999.999, "distanceFromCentroidInMeters_threshold"), (2000, ""), (2000.001, ""), (None, "")]:
        assert gbif.record_exclusion_reason(row(distanceFromCentroidInMeters=value), config) == expected
    config["gbif_missing_centroid_distance"] = "exclude"
    assert gbif.record_exclusion_reason(row(), config) == "distanceFromCentroidInMeters_unknown"
    with pytest.raises(ValueError, match="removed"):
        gbif.effective_config({"gbif_max_distance_from_centroid_m": 2000})


def test_date_windows_missingness_and_establishment_are_explicit():
    config = gbif.effective_config({"gbif_year_min": 2000, "gbif_year_max": 2020,
                                  "gbif_include_establishment_means": "native"})
    assert gbif.record_exclusion_reason(row(year=2010), config) == "establishment_not_included_or_unknown"
    assert gbif.record_exclusion_reason(row(establishmentMeans="native"), config) == "date_unknown"
    assert gbif.record_exclusion_reason(row(eventDate="2000/2020", establishmentMeans="native"), config) == ""
    assert gbif.record_exclusion_reason(row(eventDate="1999/2001", year=2000, establishmentMeans="native"), config) == "date_outside_window"
    config["gbif_missing_date"] = "keep"
    assert gbif.record_exclusion_reason(row(establishmentMeans="native"), config) == ""


def test_coordinate_uncertainty_threshold_does_not_claim_missing_is_precise():
    config = gbif.effective_config({"gbif_max_coordinate_uncertainty_m": 100})
    assert gbif.record_exclusion_reason(row(coordinateUncertaintyInMeters=100), config) == ""
    assert gbif.record_exclusion_reason(row(coordinateUncertaintyInMeters=101), config) == "coordinateUncertaintyInMeters_threshold"
    config["gbif_missing_uncertainty"] = "exclude"
    assert gbif.record_exclusion_reason(row(), config) == "coordinateUncertaintyInMeters_unknown"


def test_download_issue_column_is_interpreted_like_search_issues():
    record = gbif.normalize_record({**row(), "issues": None, "issue": "COUNTRY_COORDINATE_MISMATCH;OTHER_ISSUE"})
    assert gbif.record_exclusion_reason(record, gbif.effective_config({})) == "geospatial_issue"


def test_antimeridian_poles_and_degenerate_circular_mean():
    assert gbif.minimal_longitude_interval([179, -179]) == (179, -179, 2)
    assert gbif.minimal_longitude_interval([180, -180]) == (-180, -180, 0)
    assert gbif.circular_mean_longitude([0, 180]) is None
    assert gbif.circular_mean_longitude([90, -90]) is None
    assert gbif.circular_mean_longitude([179, -179]) == -180
    for lat in [0, 90, -90]:
        assert gbif.occupied_grid_area_km2([(lat, -180, "")], 1) == gbif.occupied_grid_area_km2([(lat, -180, ""), (lat, 180, "")], 1)
    assert gbif.occupied_grid_area_km2([(90, 0, ""), (90, 55, "")], 1) == gbif.occupied_grid_area_km2([(90, -180, "")], 1)


@pytest.mark.parametrize("grid", [1, 7, 100, 180])
def test_complete_grid_tiles_sphere_with_clipped_end_cells(grid):
    # Independent oracle: sum over every cell must equal the sphere surface.
    points = [(min(-90 + (i + .5) * grid, 89.99), min(-180 + (j + .5) * grid, 179.99), "")
              for i in range(math.ceil(180 / grid)) for j in range(math.ceil(360 / grid))]
    area = gbif.occupied_grid_area_km2(points, grid)
    assert area == pytest.approx(4 * math.pi * gbif.EARTH_RADIUS_KM ** 2, rel=1e-10)


@pytest.mark.parametrize("grid", [0, -1, float("nan"), float("inf"), 181])
def test_invalid_grid_rejected_even_for_empty_data(grid):
    with pytest.raises(ValueError):
        gbif.occupied_grid_area_km2([], grid)


def test_sampling_effort_changes_record_mean_not_duplicate_cell_area():
    balanced = [(10, 10, ""), (30, 30, "")]
    clustered = [(10, 10, "")] * 99 + [(30, 30, "")]
    a, b = [gbif.build_gbif_distribution_metrics(points, 1) for points in (balanced, clustered)]
    assert a["gbif_observed_record_mean_lat"] == 20
    assert b["gbif_observed_record_mean_lat"] == pytest.approx(10.2)
    assert a["gbif_observed_occupied_grid_area_km2"] == b["gbif_observed_occupied_grid_area_km2"]
    truncated = gbif.build_gbif_distribution_metrics(clustered[:50], 1)
    assert truncated["gbif_observed_northern_limit_lat"] == 10
    assert truncated["gbif_observed_occupied_grid_area_km2"] < a["gbif_observed_occupied_grid_area_km2"]


def make_local_files(tmp_path, records, verified=True):
    path = tmp_path / "occurrences.tsv"
    fields = list(gbif.RECORD_FIELDS)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows([{**record, "issues": ";".join(record["issues"] or [])} for record in records])
    mapping = tmp_path / "map.tsv"
    mapping.write_text("species\ttaxon_key\tscientific_name\nTest_species\t100\tTest species\n")
    config = {"gbif_occurrence_file": str(path), "gbif_taxon_map": str(mapping)}
    if verified:
        meta = tmp_path / "download.json"
        meta.write_text(json.dumps({"status": "SUCCEEDED", "key": "test-download", "doi": "10.example/test",
                                    "totalRecords": len(records), "request": {"predicate": {"type": "equals", "key": "TAXON_KEY", "value": "100"}}}))
        config["gbif_download_metadata"] = str(meta)
    return config


def acquire(tmp_path, config, fetch_json=None):
    def no_network(**kwargs):
        pytest.fail("Local download must not make network requests")
    return gbif.fetch_gbif_distribution_table("gbif", config, ["Test_species"], tmp_path / "cache", 1, False, fetch_json or no_network)


def test_local_download_snapshot_reuse_and_refilter(tmp_path):
    config = make_local_files(tmp_path, [row(1, year=1990), row(2, lat=20, year=2020)])
    frame = acquire(tmp_path, config)
    bundle = frame.attrs["gbif_bundle"]
    assert bundle["quality"][0]["status"] == "complete_download"
    assert bundle["quality"][0]["retained_count"] == 2
    snapshot = Path(bundle["records_path"])
    assert gbif.file_sha256(snapshot) == bundle["records_sha256"]
    saved = list(gbif.read_snapshot(snapshot))
    assert len(saved) == 3
    again = acquire(tmp_path, config)
    assert again.attrs["gbif_bundle"]["records_path"] == str(snapshot)
    newer = acquire(tmp_path, {**config, "gbif_year_min": 2000})
    assert newer.iloc[0]["gbif_observed_northern_limit_lat"] == 20
    assert newer.attrs["gbif_bundle"]["quality"][0]["retained_count"] == 1
    assert snapshot.is_file()  # immutable prior acquisition


def test_unverified_local_file_is_only_descriptive_and_strict_preserves_output(tmp_path):
    config = make_local_files(tmp_path, [row()], verified=False)
    frame = acquire(tmp_path, config)
    assert frame.iloc[0]["gbif_observed_northern_limit_lat"] is None
    bundle = frame.attrs["gbif_bundle"]
    assert bundle["observations"][0]["gbif_observed_northern_limit_lat"] == 10
    assert bundle["quality"][0]["status"] == "complete_local_file"
    with pytest.raises(ValueError, match="incomplete"):
        acquire(tmp_path, {**config, "gbif_require_complete": "yes"})


def test_download_count_verification_rejects_an_extract(tmp_path):
    config = make_local_files(tmp_path, [row()])
    meta = Path(config["gbif_download_metadata"])
    payload = json.loads(meta.read_text())
    payload["totalRecords"] = 2
    meta.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="matching totalRecords"):
        acquire(tmp_path, config)


@pytest.mark.parametrize("match", [{}, {"usageKey": 100, "rank": "GENUS", "confidence": 100, "matchType": "HIGHERRANK"},
                                  {"speciesKey": 100, "rank": "SPECIES", "matchType": "EXACT"},
                                  {"speciesKey": 100, "rank": "SPECIES", "confidence": 89, "matchType": "EXACT"}])
def test_unapproved_species_match_never_queries_occurrences(tmp_path, match):
    def fetch_json(url, timeout):
        assert "species/match" in url
        return match
    frame = acquire(tmp_path, {}, fetch_json)
    assert frame.attrs["gbif_bundle"]["quality"][0]["status"] == "taxon_unresolved"
    assert frame.iloc[0, 1:].isna().all()


def test_input_identity_includes_indirect_file_contents(tmp_path):
    config = make_local_files(tmp_path, [row()])
    first = gbif.input_identity(config)
    mapping = Path(config["gbif_taxon_map"])
    mapping.write_text(mapping.read_text().replace("Test species", "Reviewed species"))
    assert gbif.input_identity(config) != first
    assert gbif.input_identity({**config, "gbif_use_cache": "no"}) != gbif.input_identity({**config, "gbif_use_cache": "no"})


def test_gbif_contract_requires_explicit_selection_masks_partial_and_tracks_alias(tmp_path):
    metric = "gbif_observed_northern_limit_lat"
    path = tmp_path / "traits.tsv"
    frame = pandas.DataFrame({"species": ["A", "B"], "height": [1, 2], "observed_north": [10, 30]})
    definitions = {"height": {"role": "trait", "source": "user"},
                   "observed_north": {"role": "observation", "source": "gbif", "source_column": metric}}
    bundle = {"quality": [{"species": "A", "status": "complete_search", "analysis_eligible": True},
                          {"species": "B", "status": "capped_partial", "analysis_eligible": False, "termination_reason": "record_limit"}],
              "observations": [{"species": "A", metric: 10}, {"species": "B", metric: 30}]}
    write_trait_bundle(frame, path, definitions, bundle, [])
    automatic, audit = select_species_traits(path)
    assert automatic.columns.tolist() == ["species", "height"]
    assert audit["excluded_columns"] == ["observed_north"]
    explicit, audit = select_species_traits(path, "observed_north")
    assert explicit["observed_north"].tolist() == ["10", ""]
    assert audit["masked_species"]["observed_north"][0]["reason"] == "record_limit"
    with pytest.raises(ValueError, match="not foreground"):
        select_foreground_traits(path, "observed_north")
    path.write_text(path.read_text().replace("30", "40"))
    with pytest.raises(ValueError, match="differs"):
        select_species_traits(path, "observed_north")


def test_legacy_quality_columns_are_not_automatically_tested(tmp_path):
    path = tmp_path / "traits.tsv"
    path.write_text("species\theight\tgbif_occurrence_truncated\tgbif_northern_limit_lat\nA\t2\t1\t30\n")
    frame, audit = select_species_traits(path)
    assert frame.columns.tolist() == ["species", "height"]
    with pytest.raises(ValueError, match="metadata"):
        select_species_traits(path, "gbif_occurrence_truncated")


def test_contract_bundle_does_not_overwrite_inputs(tmp_path):
    path = tmp_path / "traits.tsv"
    path.write_text("prior input\n")
    with pytest.raises(ValueError):
        write_trait_bundle(pandas.DataFrame({"species": ["A"], "height": [1]}), path, {}, None, [path])
    assert path.read_text() == "prior input\n"
    assert not metadata_path(path).exists()


def test_trait_bundle_rolls_back_a_late_publication_failure(tmp_path, monkeypatch):
    import species_trait_contract as contract
    outputs = [tmp_path / "data.tsv", tmp_path / "quality.tsv", tmp_path / "metadata.json"]
    for index, path in enumerate(outputs):
        path.write_text(f"original {index}")
    original_replace = contract.os.replace
    failed = False

    def fail_once(source, target):
        nonlocal failed
        if Path(target) == outputs[2] and not failed:
            failed = True
            raise OSError("injected commit failure")
        return original_replace(source, target)

    monkeypatch.setattr(contract.os, "replace", fail_once)
    with pytest.raises(OSError, match="injected"):
        with contract.trait_output_transaction(outputs) as staged:
            for path in staged.values():
                path.write_text("replacement")
    assert [path.read_text() for path in outputs] == [f"original {index}" for index in range(3)]
    assert not list(tmp_path.glob(".trait-bundle-*"))


def test_sensitivity_replay_is_reproducible_and_insufficient_counts_are_missing(tmp_path):
    from analyze_gbif_sensitivity import main
    config = make_local_files(tmp_path, [row(i + 1, lat=10 + i, year=1990 + i * 10,
                                          establishmentMeans="native", datasetKey="dataset-a") for i in range(4)])
    frame = acquire(tmp_path, config)
    metadata = tmp_path / "acquisition.json"
    metadata.write_text(json.dumps(frame.attrs["gbif_bundle"]))
    output = tmp_path / "sensitivity.tsv"
    arguments = ["--metadata", str(metadata), "--output", str(output), "--grid-degrees", "1,2",
                 "--year-windows", "source,2000:2020", "--selection", "records,one_per_cell",
                 "--countries", "source;US", "--missing-date", "keep,exclude",
                 "--sample-sizes", "all,2,9", "--replicates", "3", "--seed", "37"]
    assert main(arguments) == 0
    first = output.read_bytes()
    assert main(arguments) == 0
    assert output.read_bytes() == first
    result = pandas.read_csv(output, sep="\t")
    insufficient = result[result["sample_size_requested"] == "9"]
    assert len(insufficient) > 0
    assert insufficient["status"].eq("insufficient_records").all()
    assert insufficient[list(gbif.METRIC_DEFINITIONS)].isna().all().all()
    absent_region = result[result["countries"].eq("US")]
    assert absent_region["selected_records"].eq(0).all()
    assert absent_region[list(gbif.METRIC_DEFINITIONS)].isna().all().all()
    assert result["selected_ids_sha256"].str.len().eq(64).all()
    provenance = json.loads(Path(str(output) + ".metadata.json").read_text())
    assert provenance["records_sha256"] == frame.attrs["gbif_bundle"]["records_sha256"]
    snapshot = Path(provenance["records_path"])
    snapshot.write_bytes(snapshot.read_bytes() + b"corrupted")
    with pytest.raises(ValueError, match="hash"):
        main(arguments)
    assert output.read_bytes() == first


def test_explicit_cli_defaults_override_database_source_settings():
    from generate_species_trait import apply_gbif_cli_overrides, build_arg_parser
    args = build_arg_parser().parse_args(["--gbif-page-size", "300", "--gbif-min-match-confidence", "90"])
    resolved = apply_gbif_cli_overrides({"gbif_page_size": "1", "gbif_min_match_confidence": "99"}, args)
    assert resolved["gbif_page_size"] == "300"
    assert resolved["gbif_min_match_confidence"] == "90.0"


@pytest.mark.parametrize("species_key", [999, None, ""])
def test_search_rejects_missing_or_different_species_key(species_key):
    result, _ = search([row(speciesKey=species_key)])
    assert result["status"] == "invalid_response"
    assert result["termination_reason"] == "record_taxon_mismatch"


@pytest.mark.parametrize("update", [{"totalRecords": True}, {"request": {}}, {"request": {"predicate": {}}}])
def test_local_download_rejects_incomplete_verification(tmp_path, update):
    config = make_local_files(tmp_path, [row()])
    path = Path(config["gbif_download_metadata"])
    path.write_text(json.dumps({**json.loads(path.read_text()), **update}))
    with pytest.raises(ValueError, match="metadata"):
        acquire(tmp_path, config)


@pytest.mark.parametrize("format_name", ["gz", "zip"])
def test_compressed_local_download_matches_plain_table(tmp_path, format_name):
    import gzip
    import zipfile
    config = make_local_files(tmp_path, [row(), row(2, lat=30)])
    source = Path(config["gbif_occurrence_file"])
    compressed = tmp_path / ("records." + format_name)
    if format_name == "gz":
        with gzip.open(compressed, "wb") as handle:
            handle.write(source.read_bytes())
    else:
        with zipfile.ZipFile(compressed, "w") as archive:
            archive.write(source, "occurrence.csv")
    frame = acquire(tmp_path, {**config, "gbif_occurrence_file": str(compressed)})
    assert frame.iloc[0]["gbif_observed_northern_limit_lat"] == 30
    assert frame.attrs["gbif_bundle"]["quality"][0]["status"] == "complete_download"


def test_refresh_and_corrupt_snapshots_reacquire_while_partial_is_not_reused(tmp_path, monkeypatch):
    config = make_local_files(tmp_path, [row()])
    calls = []
    original = gbif.read_local_download

    def counted(*args):
        calls.append(1)
        return original(*args)

    monkeypatch.setattr(gbif, "read_local_download", counted)
    frame = acquire(tmp_path, config)
    acquire(tmp_path, config)
    assert len(calls) == 1
    snapshot = Path(frame.attrs["gbif_bundle"]["records_path"])
    snapshot.write_bytes(b"corrupt")
    acquire(tmp_path, config)
    assert len(calls) == 2
    fresh = {**config, "gbif_use_cache": "no"}
    acquire(tmp_path, fresh)
    acquire(tmp_path, fresh)
    assert len(calls) == 4
    unverified = {**config, "gbif_download_metadata": ""}
    acquire(tmp_path, unverified)
    acquire(tmp_path, unverified)
    assert len(calls) == 6


@pytest.mark.parametrize("issue", ["COORDINATE_REPROJECTION_FAILED", "COORDINATE_REPROJECTION_SUSPICIOUS"])
def test_unsafe_reprojection_is_excluded_from_local_records(issue):
    config = gbif.effective_config({})
    assert gbif.record_exclusion_reason(row(issues=[issue]), config) == "geospatial_issue"
    assert gbif.record_exclusion_reason(row(issues=["COORDINATE_REPROJECTED"]), config) == ""


@pytest.mark.parametrize("date", ["2020-13", "2020-02-30", "2020-invalid", "2020-01-01T25:00:00", "0000"])
def test_malformed_dates_are_unknown_even_when_the_year_is_parseable(date):
    assert gbif.year_interval(row(eventDate=date, year=2020)) is None


@pytest.mark.parametrize("date", ["2020", "2020-02", "2020-02-29", "2020-02-29T12:00:00Z"])
def test_valid_date_precision_preserves_year(date):
    assert gbif.year_interval(row(eventDate=date)) == (2020, 2020)


def test_booleans_are_not_numeric_coordinates():
    assert gbif.gbif_parse_float(True) is None


def test_poles_do_not_supply_artificial_longitude_boundaries_or_means():
    polar = gbif.build_gbif_distribution_metrics([(90, 55, ""), (-90, 10, "")], 1)
    for name in gbif.METRIC_DEFINITIONS:
        if name.endswith("_lon") or "longitudinal_breadth" in name:
            assert polar[name] is None
    mixed = gbif.build_gbif_distribution_metrics([(90, 55, ""), (20, 40, "")], 1)
    assert mixed["gbif_observed_record_circular_mean_lon"] == pytest.approx(40)
    assert mixed["gbif_observed_longitudinal_breadth_deg"] == 0
    assert gbif.gbif_row_to_point(row(lat=90, lon=55), gbif.effective_config({}))[:2] == (90, -180)


def test_foreground_contract_preserves_boolean_colors(tmp_path):
    from csubst_scan_candidate_sites import write_trait_color_tables
    path = tmp_path / "traits.tsv"
    path.write_text("species\ttrait\nA\tTrue\nB\tFalse\nC\t\n")
    outputs = write_trait_color_tables(path, ["trait"], tmp_path / "colors")
    colors = pandas.read_csv(outputs["trait"], sep="\t")
    assert colors["color"].tolist() == ["firebrick", "black", "black"]


def test_taxonomic_aliases_cannot_become_independent_species(tmp_path):
    def fetch_json(url, timeout):
        if "species/match" in url:
            return {"speciesKey": 100, "rank": "SPECIES", "matchType": "EXACT", "confidence": 100}
        return {"count": 0}
    frame = gbif.fetch_gbif_distribution_table("gbif", {}, ["Name_a", "Name_b"], tmp_path, 1, False, fetch_json)
    quality = frame.attrs["gbif_bundle"]["quality"]
    assert all(record["status"] == "complete_search" and not record["analysis_eligible"] for record in quality)
    assert all(record["ineligibility_reason"] == "shared_species_key" for record in quality)


def test_outer_workflow_cache_cannot_freeze_incomplete_acquisitions(tmp_path):
    output = tmp_path / "traits.tsv"
    sidecar = metadata_path(output)
    sidecar.write_text(json.dumps({"gbif": {"quality": [{"status": "interrupted", "analysis_eligible": False}]}}))
    assert gbif.input_identity({}, output) != gbif.input_identity({}, output)
    sidecar.write_text(json.dumps({"gbif": {"quality": [{"status": "complete_search", "analysis_eligible": True}]}}))
    assert gbif.input_identity({}, output) == gbif.input_identity({}, output)


@pytest.mark.parametrize("text", ["species\tx\na\t1\t2\n", "species\tx\na\n"])
def test_trait_contract_rejects_wrong_field_counts(tmp_path, text):
    path = tmp_path / "malformed.tsv"
    path.write_text(text)
    with pytest.raises(ValueError, match="number of fields"):
        select_species_traits(path)


def test_trait_transaction_rejects_hardlinked_outputs(tmp_path):
    import os

    from species_trait_contract import trait_output_transaction
    first, second = tmp_path / "a", tmp_path / "b"
    first.write_text("original")
    os.link(first, second)
    with pytest.raises(ValueError, match="hard links"):
        with trait_output_transaction([first, second]):
            pytest.fail("Aliased outputs must be rejected before writing")
    assert first.read_text() == second.read_text() == "original"


@pytest.mark.parametrize("malformation", ["short_row", "duplicate_header"])
def test_local_download_rejects_ambiguous_or_short_rows(tmp_path, malformation):
    config = make_local_files(tmp_path, [row()])
    path = Path(config["gbif_occurrence_file"])
    lines = path.read_text().splitlines()
    if malformation == "short_row":
        lines[1] = lines[1].rsplit("\t", 1)[0]
    else:
        lines[0] += "\tdecimalLatitude"
        lines[1] += "\t89"
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="fields|headers"):
        acquire(tmp_path, config)


@pytest.mark.parametrize("bad", ["{", "null", "[]", '{"records_path": []}', '{"value": NaN}'])
def test_corrupt_cache_manifest_is_reacquired(tmp_path, bad):
    config = make_local_files(tmp_path, [row()])
    original = acquire(tmp_path, config)
    cache = next((tmp_path / "cache/gbif").glob("*.json"))
    cache.write_text(bad)
    restored = acquire(tmp_path, config)
    assert restored.to_dict() == original.to_dict()
    assert json.loads(cache.read_text())["quality"][0]["analysis_eligible"] is True


def test_generator_stats_and_sidecars_publish_as_one_bundle(tmp_path, monkeypatch):
    import generate_species_trait as generator
    from species_trait_contract import sidecar_paths
    config = make_local_files(tmp_path, [row()])
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("provider\tid\tspecies_key\nlocal\t1\tTest_species\n")
    output, stats = tmp_path / "traits.tsv", tmp_path / "stats.json"
    args = ["--download-manifest", str(manifest), "--databases", "gbif",
            "--trait-plan", str(tmp_path / "unused-plan"), "--database-sources", str(tmp_path / "unused-sources"),
            "--downloads-dir", str(tmp_path / "cache"), "--output", str(output), "--stats-output", str(stats)]
    for key, value in config.items():
        args += ["--" + key.replace("_", "-"), value]
    assert generator.main(args) == 0
    files = [output, stats, *sidecar_paths(output).values()]
    original = {path: path.read_bytes() for path in files}
    replace = generator.os.replace
    failed = False
    def fail_once(source, target):
        nonlocal failed
        if Path(target) == metadata_path(output) and not failed:
            failed = True
            raise OSError("publication failure")
        return replace(source, target)
    monkeypatch.setattr(generator.os, "replace", fail_once)
    with pytest.raises(OSError, match="publication failure"):
        generator.main(args)
    assert {path: path.read_bytes() for path in files} == original
    monkeypatch.setattr(generator.os, "replace", replace)
    for destination in [metadata_path(output), Path(config["gbif_occurrence_file"])]:
        with pytest.raises(SystemExit):
            generator.main(args[:-2] + ["--stats-output", str(destination)])
        assert {path: path.read_bytes() for path in files} == original



def test_cache_summary_corruption_cannot_change_verified_observations(tmp_path):
    config = make_local_files(tmp_path, [row(lat=12)])
    original = acquire(tmp_path, config)
    cache = next((tmp_path / "cache/gbif").glob("*.json"))
    value = json.loads(cache.read_text())
    value["analysis_rows"][0]["gbif_observed_northern_limit_lat"] = 89
    cache.write_text(json.dumps(value))
    restored = acquire(tmp_path, config)
    assert restored.to_dict() == original.to_dict()



def test_gbif_preflight_skips_unused_files_and_lists_selected_indirect_inputs(tmp_path, capsys):
    import generate_species_trait as generator
    config = make_local_files(tmp_path, [row()])
    sources = tmp_path / "sources.tsv"
    pandas.DataFrame([{"database": "gbif", **config}]).to_csv(sources, sep="\t", index=False)
    base = ["--database-sources", str(sources), "--trait-plan", str(tmp_path / "unused-plan")]
    assert generator.main([*base, "--databases", "gbif", "--print-gbif-input-files"]) == 0
    assert set(capsys.readouterr().out.splitlines()) == set(config.values())
    Path(config["gbif_occurrence_file"]).unlink()
    assert generator.main([*base, "--databases", "austraits", "--print-gbif-input-identity"]) == 0
    assert capsys.readouterr().out.strip() == "not_requested"
    with pytest.raises(FileNotFoundError):
        generator.main([*base, "--databases", "gbif", "--print-gbif-input-identity"])


def test_rsc_preparation_rolls_back_sidecar_and_protects_inputs(tmp_path, monkeypatch):
    import reconciled_speciation_contrast as rsc
    import species_trait_contract as contract
    expression, traits = tmp_path / "expression.tsv", tmp_path / "traits.tsv"
    expression.write_text("gene\tvalue\nGenus_a_g1\t1\nGenus_b_g1\t2\n")
    traits.write_text("species\tx\nGenus_a\t1\nGenus_b\t2\n")
    outputs = {name: tmp_path / (name + ".tsv") for name in ("expression-output", "species-traits-output", "analysis-plan-output", "metadata-output")}
    base = ["prepare", "--expression", str(expression), "--species-traits", str(traits)]
    args = base + [token for name, path in outputs.items() for token in ("--" + name, str(path))]
    assert rsc.main(args) == 0
    sidecar = contract.metadata_path(outputs["species-traits-output"])
    originals = {path: path.read_bytes() for path in [*outputs.values(), sidecar]}
    replace = contract.os.replace
    failed = False
    def fail_once(source, target):
        nonlocal failed
        if Path(target) == sidecar and not failed:
            failed = True
            raise OSError("sidecar publication failed")
        return replace(source, target)
    monkeypatch.setattr(contract.os, "replace", fail_once)
    with pytest.raises(SystemExit) as error:
        rsc.main(args)
    assert error.value.code == 2
    assert {path: path.read_bytes() for path in originals} == originals
    monkeypatch.setattr(contract.os, "replace", replace)
    with pytest.raises(SystemExit) as error:
        rsc.main(args + ["--metadata-output", str(traits)])
    assert error.value.code == 2
    assert traits.read_text() == "species\tx\nGenus_a\t1\nGenus_b\t2\n"



def test_core_binarization_preserves_unknowns_and_excludes_observations(tmp_path):
    import os
    import subprocess
    core = (SUPPORT.parent / "core/gg_gene_evolution_core.sh").read_text()
    function = core[core.index("binarize_species_trait() {"):core.index("normalize_mapnh_params_for_mapnh_v1() {")]
    source, output = tmp_path / "source.tsv", tmp_path / "binary.tsv"
    source.write_text("species\tbinary\tsize\tbool\tgbif_observed_record_mean_lat\na\t0\t1\ttrue\t10\nb\t1\t3\tfalse\t20\nc\tNA\tNA\tunknown\t30\n")
    subprocess.run(["bash", "-c", function + '\nbinarize_species_trait "$1" "$2"', "test", str(source), str(output)],
                   env={**os.environ, "gg_support_dir": str(SUPPORT)}, check=True)
    frame = pandas.read_csv(output, sep="\t", keep_default_na=False)
    assert frame.columns.tolist() == ["species", "binary", "size", "bool"]
    assert frame.iloc[0].tolist() == ["a", "0", "0", "1"]
    assert frame.iloc[1].tolist() == ["b", "1", "1", "0"]
    assert frame.iloc[2].tolist() == ["c", "", "", ""]



def test_custom_database_name_cannot_strip_gbif_observation_roles(tmp_path):
    import generate_species_trait as generator
    with pytest.raises(ValueError, match="must use database=gbif"):
        generator.load_database_table("custom", {"acquisition_mode": "gbif_distribution"}, [], ["Test_species"], tmp_path, 1, False)


def test_numeric_and_observation_contracts_compose_and_preserve_derived_types(tmp_path):
    from species_trait_contract import select_analysis_traits, write_derived_contract
    from species_trait_schema import schema_path, schema_payload, select_traits

    path = tmp_path / "traits.tsv"
    metric = "gbif_observed_northern_limit_lat"
    frame = pandas.DataFrame({"species": ["A", "B"], "height": [1, 2], "habit": ["tree", "herb"],
                              "code": [1, 2], "observed_north": [10, 30]})
    definitions = {"observed_north": {"role": "observation", "source": "gbif", "source_column": metric}}
    bundle = {"quality": [{"species": "A", "status": "complete_search", "analysis_eligible": True},
                          {"species": "B", "status": "capped_partial", "analysis_eligible": False}],
              "observations": []}
    write_trait_bundle(frame, path, definitions, bundle, [])
    schema_path(path).write_bytes(schema_payload(path.read_bytes(),
                                 {"height": "numeric", "habit": "text", "code": "categorical", "observed_north": "numeric"}))
    selected, audit = select_analysis_traits(path)
    assert selected.columns.tolist() == ["species", "height"]
    reasons = {row["trait"]: row["reason"] for row in audit["type_selection"]}
    assert reasons == {"height": "", "habit": "declared_text", "code": "declared_categorical",
                       "observed_north": "observation_contract_excluded"}
    selected, audit = select_analysis_traits(path, "observed_north")
    assert selected.observed_north.tolist() == ["10", ""]
    derived = tmp_path / "selected.tsv"
    selected.to_csv(derived, sep="\t", index=False)
    write_derived_contract(derived, audit)
    assert select_traits(derived)[0]["trait"] == "observed_north"
    assert select_analysis_traits(derived, "observed_north")[0].observed_north.tolist() == ["10", ""]
    with pytest.raises(ValueError, match="explicitly encode"):
        select_analysis_traits(path, "code")
    schema_path(path).write_text(schema_path(path).read_text().replace('"numeric"', '"text"', 1))
    assert select_analysis_traits(path, "observed_north")[0].columns.tolist() == ["species", "observed_north"]
    stale_schema = json.loads(schema_path(path).read_text())
    stale_schema["table_sha256"] = "0" * 64
    schema_path(path).write_text(json.dumps(stale_schema))
    with pytest.raises(ValueError, match="does not match"):
        select_analysis_traits(path)
