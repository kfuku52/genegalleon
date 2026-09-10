"""Auditable descriptions of GBIF occurrence records, not estimates of species ranges.

The same normalized records and metric definitions are used for search, local
downloads and sensitivity analyses. No occurrence download is requested here.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
import tempfile
import zipfile
from collections import Counter
from contextlib import contextmanager
from datetime import datetime, timezone
from io import TextIOWrapper
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence
from urllib.parse import urlencode, urlsplit

import pandas

SCHEMA_VERSION = 2
DEFAULT_GBIF_API = "https://api.gbif.org/v1/"
DEFAULT_GBIF_PAGE_SIZE = 300
DEFAULT_GBIF_MAX_OCCURRENCES_PER_SPECIES = 100000
GBIF_SEARCH_HARD_LIMIT = 100000
DEFAULT_GBIF_GRID_DEGREES = 1.0
DEFAULT_GBIF_MIN_MATCH_CONFIDENCE = 90.0
EARTH_RADIUS_KM = 6371.0088
COMPLETE_STATES = {"complete_search", "complete_download"}

METRIC_DEFINITIONS = {
    "gbif_observed_northern_limit_lat": ("degree", "Maximum retained record latitude"),
    "gbif_observed_southern_limit_lat": ("degree", "Minimum retained record latitude"),
    "gbif_observed_latitudinal_breadth_deg": ("degree", "Maximum minus minimum retained latitude"),
    "gbif_observed_western_limit_lon": ("degree", "Start of the shortest retained nonpolar longitude arc"),
    "gbif_observed_eastern_limit_lon": ("degree", "End of the shortest retained nonpolar longitude arc"),
    "gbif_observed_longitudinal_breadth_deg": ("degree", "Length of the shortest retained nonpolar longitude arc"),
    "gbif_observed_occupied_grid_area_km2": ("km2", "Sum of spherical areas of cells containing retained records; not AOO or habitat area"),
    "gbif_observed_record_mean_lat": ("degree", "Record-weighted arithmetic mean latitude; not an area centroid"),
    "gbif_observed_record_circular_mean_lon": ("degree", "Record-weighted circular mean nonpolar longitude; undefined with no nonpolar records or a near-zero resultant"),
    "gbif_observed_country_count": ("count", "Distinct nonmissing country codes in retained records"),
}
DEFAULT_GBIF_DISTRIBUTION_TRAITS = tuple((name, "numeric") for name in METRIC_DEFINITIONS)

# Retain fields needed to audit/refilter the observation process and cite sources.
# Identifying collectors is unnecessary for these metrics.
RECORD_FIELDS = (
    "gbifID", "speciesKey", "taxonKey", "acceptedTaxonKey", "scientificName",
    "decimalLatitude", "decimalLongitude", "countryCode", "geodeticDatum",
    "coordinateUncertaintyInMeters", "coordinatePrecision", "distanceFromCentroidInMeters",
    "occurrenceStatus", "hasGeospatialIssues", "issues", "basisOfRecord",
    "year", "month", "day", "eventDate", "datasetKey", "occurrenceID",
    "eventID", "parentEventID", "samplingProtocol", "samplingEffort",
    "sampleSizeValue", "sampleSizeUnit", "establishmentMeans", "degreeOfEstablishment",
    "license", "references", "dataGeneralizations", "informationWithheld",
)
GEOSPATIAL_ISSUES = {
    "ZERO_COORDINATE", "COUNTRY_COORDINATE_MISMATCH", "COORDINATE_INVALID", "COORDINATE_OUT_OF_RANGE",
    "COORDINATE_REPROJECTION_FAILED", "COORDINATE_REPROJECTION_SUSPICIOUS",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def file_sha256(path: Path | str) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_text(value: object) -> str:
    return json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False, separators=(",", ":"))


def tokens(value: object) -> set[str]:
    return {item.strip().upper() for item in str(value or "").split(",") if item.strip()}


def gbif_parse_float(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        parsed = float(value)
    except (ValueError, TypeError):
        return None
    return parsed if math.isfinite(parsed) else None


def normalize_longitude(value: float) -> float:
    """Canonical half-open interval; +180 and -180 are the same location."""
    if not math.isfinite(value):
        raise ValueError("Longitude must be finite")
    return ((float(value) + 180.0) % 360.0) - 180.0


def circular_mean_longitude(longitudes: Sequence[float]) -> float | None:
    if not longitudes:
        return None
    mean_sin = math.fsum(math.sin(math.radians(lon)) for lon in longitudes) / len(longitudes)
    mean_cos = math.fsum(math.cos(math.radians(lon)) for lon in longitudes) / len(longitudes)
    if math.hypot(mean_sin, mean_cos) <= 1e-12:
        return None
    return normalize_longitude(math.degrees(math.atan2(mean_sin, mean_cos)))


def minimal_longitude_interval(longitudes: Sequence[float]) -> tuple[float | None, float | None, float | None]:
    if not longitudes:
        return None, None, None
    values = sorted(set(float(lon) % 360.0 for lon in longitudes))
    gaps = [values[(i + 1) % len(values)] - lon for i, lon in enumerate(values)]
    gaps[-1] += 360.0
    index = max(range(len(gaps)), key=gaps.__getitem__)
    return normalize_longitude(values[(index + 1) % len(values)]), normalize_longitude(values[index]), 360.0 - gaps[index]


def validate_grid(grid_degrees: float) -> None:
    if not math.isfinite(grid_degrees) or not 0 < grid_degrees <= 180:
        raise ValueError("gbif_grid_degrees must be finite and in (0, 180]")


def occupied_cell(point: tuple[float, float, str], grid_degrees: float) -> tuple[int, int]:
    lat, lon, _country = point
    if not math.isfinite(lat) or not -90 <= lat <= 90 or not math.isfinite(lon) or not -180 <= lon <= 180:
        raise ValueError("Invalid occurrence coordinates")
    # Shift before nextafter: adding 90 after clamping can round back to 180.
    latitude_offset = min(lat + 90.0, math.nextafter(180.0, 0.0))
    longitude_offset = normalize_longitude(lon) + 180.0
    # Longitude at a pole is undefined; put each pole in one canonical cell.
    if abs(lat) == 90:
        longitude_offset = 0.0
    return math.floor(latitude_offset / grid_degrees), math.floor(longitude_offset / grid_degrees)


def occupied_grid_area_km2(points: Sequence[tuple[float, float, str]], grid_degrees: float) -> float | None:
    validate_grid(grid_degrees)
    if not points:
        return None
    cells = {occupied_cell(point, grid_degrees) for point in points}
    areas = []
    for lat_index, lon_index in cells:
        lat_min = -90.0 + lat_index * grid_degrees
        lat_max = min(lat_min + grid_degrees, 90.0)
        lon_min = -180.0 + lon_index * grid_degrees
        lon_max = min(lon_min + grid_degrees, 180.0)
        areas.append(EARTH_RADIUS_KM ** 2 * math.radians(lon_max - lon_min)
                     * (math.sin(math.radians(lat_max)) - math.sin(math.radians(lat_min))))
    return math.fsum(areas)


def build_gbif_distribution_metrics(points: Sequence[tuple[float, float, str]], grid_degrees: float) -> dict:
    validate_grid(grid_degrees)
    if not points:
        return {name: None for name in METRIC_DEFINITIONS}
    latitudes = [point[0] for point in points]
    # A pole has no defined longitude. It contributes latitude and one cell,
    # but cannot define a longitude boundary or a directional mean.
    longitudes = [point[1] for point in points if abs(point[0]) != 90]
    west, east, breadth = minimal_longitude_interval(longitudes)
    return {
        "gbif_observed_northern_limit_lat": max(latitudes),
        "gbif_observed_southern_limit_lat": min(latitudes),
        "gbif_observed_latitudinal_breadth_deg": max(latitudes) - min(latitudes),
        "gbif_observed_western_limit_lon": west,
        "gbif_observed_eastern_limit_lon": east,
        "gbif_observed_longitudinal_breadth_deg": breadth,
        "gbif_observed_occupied_grid_area_km2": occupied_grid_area_km2(points, grid_degrees),
        "gbif_observed_record_mean_lat": math.fsum(latitudes) / len(points),
        "gbif_observed_record_circular_mean_lon": circular_mean_longitude(longitudes),
        "gbif_observed_country_count": len({point[2] for point in points if point[2]}),
    }


def effective_config(config: Mapping) -> dict:
    if str(config.get("gbif_max_distance_from_centroid_m", "")).strip():
        raise ValueError("gbif_max_distance_from_centroid_m was removed: use gbif_min_distance_from_known_centroid_m to exclude proximity to georeferencing centroids")
    result = {
        "uri": str(config.get("uri") or DEFAULT_GBIF_API).rstrip("/") + "/",
        "gbif_page_size": DEFAULT_GBIF_PAGE_SIZE,
        "gbif_max_occurrences_per_species": DEFAULT_GBIF_MAX_OCCURRENCES_PER_SPECIES,
        "gbif_grid_degrees": DEFAULT_GBIF_GRID_DEGREES,
        "gbif_min_match_confidence": DEFAULT_GBIF_MIN_MATCH_CONFIDENCE,
        "gbif_max_coordinate_uncertainty_m": None,
        "gbif_min_distance_from_known_centroid_m": None,
        "gbif_year_min": None,
        "gbif_year_max": None,
        "gbif_countries": "",
        "gbif_include_basis_of_record": "",
        "gbif_exclude_basis_of_record": "",
        "gbif_include_establishment_means": "",
        "gbif_missing_date": "exclude",
        "gbif_missing_uncertainty": "keep",
        "gbif_missing_centroid_distance": "keep",
        "gbif_use_cache": "yes",
        "gbif_require_complete": "no",
        "gbif_occurrence_file": "",
        "gbif_taxon_map": "",
        "gbif_download_metadata": "",
    }
    for key in result:
        value = config.get(key)
        if value is not None and str(value).strip() != "":
            result[key] = value
    api_parts = urlsplit(str(result["uri"]))
    if api_parts.scheme not in {"http", "https"} or not api_parts.hostname or api_parts.username or api_parts.password or api_parts.query or api_parts.fragment:
        raise ValueError("gbif_api must be an HTTP(S) base URI without credentials, query parameters or fragments")
    result["uri"] = str(result["uri"]).rstrip("/") + "/"
    for key in ("gbif_page_size", "gbif_max_occurrences_per_species", "gbif_year_min", "gbif_year_max"):
        if result[key] is not None:
            value = gbif_parse_float(result[key])
            if value is None or value < 1 or value != int(value):
                raise ValueError(f"{key} must be a positive integer")
            result[key] = int(value)
    result["gbif_page_size"] = min(result["gbif_page_size"], DEFAULT_GBIF_PAGE_SIZE)
    for key in ("gbif_grid_degrees", "gbif_min_match_confidence", "gbif_max_coordinate_uncertainty_m", "gbif_min_distance_from_known_centroid_m"):
        if result[key] is not None:
            value = gbif_parse_float(result[key])
            if value is None or value < 0:
                raise ValueError(f"{key} must be a finite nonnegative number")
            result[key] = value
    validate_grid(result["gbif_grid_degrees"])
    if not 0 <= result["gbif_min_match_confidence"] <= 100:
        raise ValueError("gbif_min_match_confidence must be in [0, 100]")
    if result["gbif_year_min"] and result["gbif_year_max"] and result["gbif_year_min"] > result["gbif_year_max"]:
        raise ValueError("gbif_year_min must not exceed gbif_year_max")
    for key in ("gbif_missing_date", "gbif_missing_uncertainty", "gbif_missing_centroid_distance"):
        result[key] = str(result[key]).lower()
        if result[key] not in {"keep", "exclude"}:
            raise ValueError(f"{key} must be keep or exclude")
    for key in ("gbif_use_cache", "gbif_require_complete"):
        value = str(result[key]).lower()
        if value not in {"yes", "no", "true", "false", "1", "0"}:
            raise ValueError(f"{key} must be yes or no")
        result[key] = "yes" if value in {"yes", "true", "1"} else "no"
    for key in ("gbif_countries", "gbif_include_basis_of_record", "gbif_exclude_basis_of_record", "gbif_include_establishment_means"):
        result[key] = ",".join(sorted(tokens(result[key])))
    for key in ("gbif_occurrence_file", "gbif_taxon_map", "gbif_download_metadata"):
        if result[key]:
            result[key] = str(Path(result[key]).expanduser().resolve())
    if bool(result["gbif_occurrence_file"]) != bool(result["gbif_taxon_map"]):
        raise ValueError("Local GBIF records require both gbif_occurrence_file and gbif_taxon_map")
    if result["gbif_download_metadata"] and not result["gbif_occurrence_file"]:
        raise ValueError("gbif_download_metadata requires gbif_occurrence_file")
    return result


def input_identity(config: Mapping, previous_output: Path | None = None) -> str:
    """Fingerprint local files referred to indirectly by a database-source map."""
    resolved = effective_config(config)
    identity = {"configuration": resolved, "producer_sha256": file_sha256(Path(__file__)),
                "files": {key: file_sha256(Path(resolved[key]))
                          for key in ("gbif_occurrence_file", "gbif_taxon_map", "gbif_download_metadata") if resolved[key]}}
    if resolved["gbif_use_cache"] == "no":
        identity["fresh_acquisition_requested_utc"] = utc_now()
    if previous_output is not None:
        sidecar = Path(str(previous_output) + ".metadata.json")
        if sidecar.is_file():
            try:
                quality = json.loads(sidecar.read_text()).get("gbif", {}).get("quality", [])
                retry = any(row.get("status") not in COMPLETE_STATES or not row.get("analysis_eligible") for row in quality)
            except (OSError, ValueError, AttributeError, TypeError):
                retry = True
            if retry:
                # The outer workflow cache must not suppress recovery of a
                # failed acquisition merely because an all-NA table exists.
                identity["incomplete_previous_acquisition_retry_utc"] = utc_now()
    return hashlib.sha256(json_text(identity).encode()).hexdigest()


def normalize_record(row: Mapping) -> dict:
    # SIMPLE_CSV and DWCA headers can use lower case or namespace URIs.
    folded = {str(key).rsplit("/", 1)[-1].lower(): value for key, value in row.items()}
    record = {key: folded.get(key.lower()) for key in RECORD_FIELDS}
    record["issues"] = record["issues"] or folded.get("issue") or []
    record["gbifID"] = record["gbifID"] or folded.get("key") or folded.get("gbifid")
    for key in ("gbifID", "speciesKey", "taxonKey", "acceptedTaxonKey"):
        if record[key] is not None:
            record[key] = str(record[key]).strip()
    for key, value in record.items():
        if isinstance(value, float) and not math.isfinite(value):
            record[key] = str(value)
    return record


def year_interval(row: Mapping) -> tuple[int, int] | None:
    date = str(row.get("eventDate") or "").strip()
    if date:
        parts = date.split("/")
        if len(parts) <= 2 and all(re.match(r"^\d{4}(?:-|$)", part) for part in parts):
            try:
                for part in parts:
                    if re.fullmatch(r"\d{4}", part):
                        datetime(int(part), 1, 1)
                    elif re.fullmatch(r"\d{4}-\d{2}", part):
                        datetime(int(part[:4]), int(part[5:]), 1)
                    else:
                        datetime.fromisoformat(part.replace("Z", "+00:00"))
            except ValueError:
                return None
            years = [int(part[:4]) for part in parts]
            if 0 < years[0] <= years[-1]:
                return years[0], years[-1]
        return None
    year = gbif_parse_float(row.get("year"))
    return (int(year), int(year)) if year is not None and year > 0 and year == int(year) else None


def record_exclusion_reason(row: Mapping, config: Mapping) -> str:
    if str(row.get("occurrenceStatus") or "").upper() != "PRESENT":
        return "not_present_or_unknown"
    issues = row.get("issues") or []
    if isinstance(issues, str):
        issues = re.split(r"[;,]", issues)
    has_issue = str(row.get("hasGeospatialIssues") or "").lower() == "true"
    if has_issue or GEOSPATIAL_ISSUES.intersection(issues):
        return "geospatial_issue"
    lat, lon = gbif_parse_float(row.get("decimalLatitude")), gbif_parse_float(row.get("decimalLongitude"))
    if lat is None or lon is None or not -90 <= lat <= 90 or not -180 <= lon <= 180:
        return "invalid_coordinates"
    if lat == 0 and lon == 0:
        return "zero_coordinate"
    country = str(row.get("countryCode") or "").upper()
    if config["gbif_countries"] and country not in tokens(config["gbif_countries"]):
        return "country"
    basis = str(row.get("basisOfRecord") or "").upper()
    if config["gbif_include_basis_of_record"] and basis not in tokens(config["gbif_include_basis_of_record"]):
        return "basis_not_included"
    if basis in tokens(config["gbif_exclude_basis_of_record"]):
        return "basis_excluded"
    establishment = str(row.get("establishmentMeans") or "").upper()
    if config["gbif_include_establishment_means"] and establishment not in tokens(config["gbif_include_establishment_means"]):
        return "establishment_not_included_or_unknown"
    interval = year_interval(row)
    if config["gbif_year_min"] or config["gbif_year_max"]:
        if interval is None:
            if config["gbif_missing_date"] == "exclude":
                return "date_unknown"
        elif ((config["gbif_year_min"] and interval[0] < config["gbif_year_min"])
              or (config["gbif_year_max"] and interval[1] > config["gbif_year_max"])):
            return "date_outside_window"
    for field, limit_key, missing_key, direction in (
        ("coordinateUncertaintyInMeters", "gbif_max_coordinate_uncertainty_m", "gbif_missing_uncertainty", "maximum"),
        ("distanceFromCentroidInMeters", "gbif_min_distance_from_known_centroid_m", "gbif_missing_centroid_distance", "minimum"),
    ):
        limit = config[limit_key]
        value = gbif_parse_float(row.get(field))
        if limit is not None:
            if value is None or value < 0:
                if config[missing_key] == "exclude":
                    return field + "_unknown"
            elif (direction == "maximum" and value > limit) or (direction == "minimum" and value < limit):
                return field + "_threshold"
    return ""


def gbif_row_to_point(row: Mapping, config: Mapping) -> tuple[float, float, str] | None:
    if record_exclusion_reason(row, config):
        return None
    lat = float(row["decimalLatitude"])
    lon = -180.0 if abs(lat) == 90 else normalize_longitude(float(row["decimalLongitude"]))
    return lat, lon, str(row.get("countryCode") or "").upper()


def gbif_occurrence_query_params(taxon_key: str, limit: int, offset: int, config: Mapping) -> dict:
    # Keep date/country/basis restrictions local too so saved snapshots can be
    # reused to compare them. The acquisition universe is recorded explicitly.
    return {"taxonKey": taxon_key, "hasCoordinate": "true", "hasGeospatialIssue": "false",
            "occurrenceStatus": "PRESENT", "limit": limit, "offset": offset}


def query_url(api: str, resource: str, params: dict) -> str:
    return api.rstrip("/") + "/" + resource + "?" + urlencode(params, doseq=True)


def fetch_gbif_occurrence_records(api: str, taxon_key: str, config: Mapping, timeout: float, fetch_json: Callable) -> dict:
    limit = min(config["gbif_max_occurrences_per_species"], GBIF_SEARCH_HARD_LIMIT)
    result = {"started_utc": utc_now(), "ended_utc": None, "status": "invalid_response",
              "termination_reason": "invalid_count", "reported_count": None, "raw_fetched": 0,
              "fetch_limit": limit, "page_count": 0, "rows": [], "counts_observed": [],
              "query": gbif_occurrence_query_params(taxon_key, 0, 0, config)}
    try:
        payload = fetch_json(url=query_url(api, "occurrence/search", result["query"]), timeout=timeout)
        count = payload.get("count") if isinstance(payload, dict) else None
        if isinstance(count, bool) or not isinstance(count, int) or count < 0:
            return result
        result["reported_count"] = count
        result["counts_observed"] = [count]
        target = min(count, limit)
        result["status"] = "complete_search" if count <= limit else "capped_partial"
        result["termination_reason"] = "zero_records" if count == 0 else ("count_reached" if count <= limit else "record_limit")
        offset = 0
        seen_ids = set()
        while offset < target:
            page_size = min(config["gbif_page_size"], target - offset)
            params = gbif_occurrence_query_params(taxon_key, page_size, offset, config)
            payload = fetch_json(url=query_url(api, "occurrence/search", params), timeout=timeout)
            result["page_count"] += 1
            if not isinstance(payload, dict) or not isinstance(payload.get("results"), list):
                result.update(status="invalid_response", termination_reason="invalid_page")
                break
            rows = payload["results"]
            if len(rows) > page_size or any(not isinstance(row, dict) for row in rows):
                result.update(status="invalid_response", termination_reason="invalid_page_rows")
                break
            result["raw_fetched"] += len(rows)
            result["rows"].extend(normalize_record(row) for row in rows)
            if any(str(row.get("speciesKey") or "") != taxon_key for row in rows):
                result.update(status="invalid_response", termination_reason="record_taxon_mismatch")
                break
            page_count = payload.get("count")
            if type(page_count) is not int or page_count != count:
                result["counts_observed"].append(page_count)
                result.update(status="inconsistent_search", termination_reason="count_changed_or_missing")
                break
            if (type(payload.get("offset")) is not int or type(payload.get("limit")) is not int
                    or payload["offset"] != offset or payload["limit"] != page_size):
                result.update(status="invalid_response", termination_reason="page_position_mismatch")
                break
            ids = [str(row.get("key") or row.get("gbifID") or "") for row in rows]
            if any(not key for key in ids) or len(set(ids)) != len(ids) or seen_ids.intersection(ids):
                result.update(status="inconsistent_search", termination_reason="duplicate_or_missing_record_id")
                break
            seen_ids.update(ids)
            offset += len(rows)
            if len(rows) < page_size:
                result.update(status="interrupted", termination_reason="short_or_empty_page")
                break
            if not isinstance(payload.get("endOfRecords"), bool) or payload["endOfRecords"] != (offset >= count):
                result.update(status="inconsistent_search", termination_reason="end_of_records_mismatch")
                break
    except Exception as exc:
        # Preserve the partial observation audit; never promote it to complete.
        result.update(status="interrupted", termination_reason="request_failed", error_type=type(exc).__name__)
    finally:
        result["ended_utc"] = utc_now()
    return result


def summarize_records(acquisition: Mapping, config: Mapping) -> tuple[dict, dict, list[dict]]:
    retained, points, seen = [], [], set()
    exclusions = Counter()
    missing_id = 0
    for row in acquisition["rows"]:
        key = str(row.get("gbifID") or "")
        if not key:
            missing_id += 1
        elif key in seen:
            exclusions["duplicate_gbif_id"] += 1
            continue
        if key:
            seen.add(key)
        reason = record_exclusion_reason(row, config)
        if reason:
            exclusions[reason] += 1
            continue
        retained.append(row)
        points.append(gbif_row_to_point(row, config))
    metrics = build_gbif_distribution_metrics(points, config["gbif_grid_degrees"])
    intervals = [interval for row in retained if (interval := year_interval(row)) is not None]
    quality = {key: value for key, value in acquisition.items() if key != "rows"}
    quality.update(
        unique_gbif_ids=len(seen), missing_gbif_ids=missing_id, retained_count=len(retained),
        excluded_count=sum(exclusions.values()), exclusions=dict(sorted(exclusions.items())),
        analysis_eligible=acquisition["status"] in COMPLETE_STATES and missing_id == 0 and exclusions["duplicate_gbif_id"] == 0,
        observed_year_min=min((item[0] for item in intervals), default=None),
        observed_year_max=max((item[1] for item in intervals), default=None),
        unknown_date_count=len(retained) - len(intervals),
        unknown_uncertainty_count=sum((value := gbif_parse_float(row.get("coordinateUncertaintyInMeters"))) is None or value < 0 for row in retained),
        unknown_centroid_distance_count=sum((value := gbif_parse_float(row.get("distanceFromCentroidInMeters"))) is None or value < 0 for row in retained),
        unknown_establishment_count=sum(not row.get("establishmentMeans") for row in retained),
        unknown_effort_count=sum(not row.get("samplingEffort") for row in retained),
        unique_location_count=len(set((point[0], point[1]) for point in points)),
        occupied_cell_count=len({occupied_cell(point, config["gbif_grid_degrees"]) for point in points}),
        circular_mean_status=("no_records" if not points else "undefined_at_poles" if all(abs(point[0]) == 90 for point in points)
                              else "undefined_resultant" if metrics["gbif_observed_record_circular_mean_lon"] is None else "defined"),
    )
    for field in ("datasetKey", "basisOfRecord", "countryCode", "establishmentMeans", "degreeOfEstablishment", "samplingProtocol", "geodeticDatum"):
        quality[field + "_counts"] = dict(sorted(Counter(str(row.get(field) or "unknown") for row in retained).items()))
    quality["unique_event_count"] = len({(row.get("datasetKey"), row["eventID"]) for row in retained if row.get("eventID")})
    return metrics, quality, retained


@contextmanager
def open_occurrence_table(path: Path):
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as archive:
            members = [name for name in archive.namelist() if not name.endswith("/") and name.lower().endswith((".csv", ".tsv"))]
            if len(members) != 1:
                raise ValueError("GBIF zip must contain exactly one SIMPLE_CSV .csv/.tsv table; extract other formats explicitly")
            with archive.open(members[0]) as raw, TextIOWrapper(raw, encoding="utf-8-sig", newline="") as handle:
                yield handle
    else:
        opener = gzip.open if path.suffix.lower() == ".gz" else open
        with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
            yield handle


def read_local_download(config: Mapping, species: Sequence[str]) -> tuple[dict, dict]:
    mapping = {}
    with Path(config["gbif_taxon_map"]).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"species", "taxon_key", "scientific_name"}.issubset(reader.fieldnames or []):
            raise ValueError("GBIF taxon map requires species, taxon_key and scientific_name (reviewed species-level mapping)")
        for row in reader:
            name = row["species"].strip()
            if name in mapping or not name or not re.fullmatch(r"[1-9]\d*", row["taxon_key"].strip()) or not row["scientific_name"].strip():
                raise ValueError("Invalid or duplicate GBIF taxon-map entry")
            mapping[name] = {"speciesKey": row["taxon_key"].strip(), "scientificName": row["scientific_name"].strip(),
                             "rank": "SPECIES", "matchType": "MANUAL", "confidence": None}
    missing = sorted(set(species) - set(mapping))
    if missing:
        raise ValueError("GBIF taxon map is missing target species: " + ", ".join(missing))
    keys = {mapping[name]["speciesKey"] for name in species}
    if len(keys) != len(species):
        raise ValueError("Each target species must have a distinct reviewed GBIF speciesKey; aliases cannot act as independent species")
    acquisitions = {key: {"rows": [], "started_utc": utc_now(), "status": "complete_local_file",
                          "termination_reason": "local_file_without_download_verification", "query": "unknown",
                          "fetch_limit": None, "page_count": 0, "reported_count": None} for key in keys}
    total_rows = 0
    with open_occurrence_table(Path(config["gbif_occurrence_file"])) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        normalized_headers = [name.rsplit("/", 1)[-1].lower() for name in (reader.fieldnames or [])]
        if len(normalized_headers) != len(set(normalized_headers)) or any(not name.strip() for name in normalized_headers):
            raise ValueError("GBIF occurrence headers must be distinct and nonempty after namespace normalization")
        headers = set(normalized_headers)
        if not {"gbifid", "specieskey", "decimallatitude", "decimallongitude", "occurrencestatus", "issue"}.issubset(headers) and not {"gbifid", "specieskey", "decimallatitude", "decimallongitude", "occurrencestatus", "issues"}.issubset(headers):
            raise ValueError("GBIF local input must be a tab-delimited interpreted occurrence table with gbifID, speciesKey, coordinates, occurrenceStatus and issue(s)")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Malformed GBIF row: wrong number of fields")
            total_rows += 1
            record = normalize_record(row)
            if record["speciesKey"] in acquisitions:
                acquisitions[record["speciesKey"]]["rows"].append(record)
    verification = {}
    if config["gbif_download_metadata"]:
        metadata_path = Path(config["gbif_download_metadata"])
        verification = json.loads(metadata_path.read_text())
        # A saved official /occurrence/download/{key} response supplies the
        # download universe/count. A truncated extract must not pass as complete.
        if (verification.get("status") != "SUCCEEDED" or not verification.get("key")
                or not verification.get("doi") or type(verification.get("totalRecords")) is not int
                or verification.get("totalRecords") != total_rows
                or not isinstance(verification.get("request"), dict)
                or not isinstance(verification["request"].get("predicate"), dict)
                or not verification["request"]["predicate"].get("type")):
            raise ValueError("GBIF download metadata must describe a SUCCEEDED download with key, DOI, request and matching totalRecords")
        for acquisition in acquisitions.values():
            acquisition.update(status="complete_download", termination_reason="download_count_verified",
                               query=verification["request"].get("predicate", "unknown"),
                               download_key=verification["key"], download_doi=verification["doi"],
                               download_created=verification.get("created"))
    for acquisition in acquisitions.values():
        acquisition.update(raw_fetched=len(acquisition["rows"]), ended_utc=utc_now(),
                           reported_count=len(acquisition["rows"]), source_file_total_records=total_rows)
    return mapping, acquisitions


def read_snapshot(path: Path) -> Iterable[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            yield json.loads(line)


def fetch_gbif_distribution_table(database: str, config: Mapping, species: Sequence[str], downloads_dir: Path,
                                  timeout: float, dry_run: bool, fetch_json: Callable) -> pandas.DataFrame | None:
    config = effective_config(config)
    producer_sha256 = file_sha256(Path(__file__))
    inputs = {key: {"path": config[key], "sha256": file_sha256(Path(config[key]))}
              for key in ("gbif_occurrence_file", "gbif_taxon_map", "gbif_download_metadata") if config[key]}
    identity = {"schema_version": SCHEMA_VERSION, "producer_sha256": producer_sha256,
                "species": sorted(species), "configuration": config, "inputs": inputs}
    cache_id = hashlib.sha256(json_text(identity).encode()).hexdigest()
    cache_path = downloads_dir / "gbif" / ("observations_" + cache_id + ".json")
    if dry_run:
        print("[dry-run] GBIF observation acquisition: " + (config["gbif_occurrence_file"] or config["uri"]), flush=True)
        return None
    # gg-cache-guard: audited - identity includes all effective parameters, source hashes, schema and producer hash; cached records are hash-checked.
    if config["gbif_use_cache"] == "yes" and cache_path.exists():
        try:
            bundle = json.loads(cache_path.read_text())
        except (ValueError, OSError):
            bundle = {}
        if not isinstance(bundle, dict):
            bundle = {}
        try:
            records_path = Path(bundle.get("records_path") or "")
            valid_bundle_hash = bundle.get("bundle_sha256") == hashlib.sha256(json_text({key: value for key, value in bundle.items() if key != "bundle_sha256"}).encode()).hexdigest()
        except (ValueError, TypeError):
            records_path, valid_bundle_hash = Path(), False
        if (valid_bundle_hash and bundle.get("identity") == identity and records_path.is_file() and bundle.get("records_sha256") == file_sha256(records_path)
                and sorted(row.get("species", "") for row in bundle.get("quality", [])) == sorted(species)
                and all(row["status"] in COMPLETE_STATES and row.get("analysis_eligible") is True for row in bundle.get("quality", []))):
            print(f"[gbif] Reusing observation snapshot acquired at {bundle['acquired_utc']}: {records_path}", flush=True)
            frame = pandas.DataFrame.from_records(bundle["analysis_rows"])
            frame.attrs["gbif_bundle"] = bundle
            return frame
    local_mapping, local_acquisitions = ({}, {})
    if config["gbif_occurrence_file"]:
        local_mapping, local_acquisitions = read_local_download(config, species)
    bundle = {"identity": identity, "acquired_utc": utc_now(), "quality": [], "observations": [],
              "analysis_rows": [], "metric_definitions": METRIC_DEFINITIONS,
              "meaning": "Descriptions of retained GBIF presence records, not unbiased estimates of biological ranges",
              "taxon_identity_rule": "Different target species must not share one GBIF speciesKey for analysis",
              "rejected_geospatial_flags": sorted(GEOSPATIAL_ISSUES),
              "geometry": {"crs": "GBIF interpreted WGS84 coordinates", "earth_radius_km": EARTH_RADIUS_KM,
                           "grid_origin_lat_lon": [-90, -180], "grid_degrees": config["gbif_grid_degrees"],
                           "longitude_domain": "[-180,180)", "end_cells": "clipped at +90/+180",
                           "pole_cell_longitude": "-180", "longitude_summary_poles": "excluded: longitude is undefined",
                           "weighting": "one per distinct gbifID; separate events may share coordinates",
                           "area_definition": "spherical cell areas; neither IUCN AOO nor habitat area"}}
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    # Snapshot is written before its hash-bearing cache manifest. Partial runs
    # remain auditable, but only complete acquisitions are eligible for reuse.
    with tempfile.NamedTemporaryFile(dir=cache_path.parent, suffix=".records.tmp", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        with gzip.open(temporary, "wt", encoding="utf-8") as handle:
            for name in species:
                match = local_mapping.get(name)
                if match is None:
                    try:
                        match = fetch_json(url=query_url(config["uri"], "species/match", {"name": name.replace("_", " "), "verbose": "false"}), timeout=timeout)
                    except Exception:
                        match = {}
                if not isinstance(match, dict):
                    match = {}
                key = str(match.get("speciesKey") or match.get("acceptedUsageKey") or match.get("usageKey") or "")
                confidence = gbif_parse_float(match.get("confidence"))
                matched = bool(re.fullmatch(r"[1-9]\d*", key) and (name in local_mapping or (str(match.get("rank") or "").upper() == "SPECIES"
                               and str(match.get("matchType") or "").upper() in {"EXACT", "FUZZY"}
                               and confidence is not None and confidence >= config["gbif_min_match_confidence"])))
                if not matched:
                    acquisition = {"rows": [], "status": "taxon_unresolved", "termination_reason": "species_match_not_accepted",
                                   "reported_count": None, "raw_fetched": 0, "started_utc": utc_now(), "ended_utc": utc_now()}
                elif local_acquisitions:
                    acquisition = local_acquisitions[key]
                else:
                    acquisition = fetch_gbif_occurrence_records(config["uri"], key, config, timeout, fetch_json)
                handle.write(json_text({"species": name, "match": match, "acquisition": {key: value for key, value in acquisition.items() if key != "rows"}}) + "\n")
                for row in acquisition["rows"]:
                    handle.write(json_text({"species": name, "record": row}) + "\n")
                metrics, quality, _retained = summarize_records(acquisition, config)
                quality.update(species=name, taxon_key=key or None, match=match)
                bundle["quality"].append(quality)
                bundle["observations"].append({"species": name, **metrics})
                bundle["analysis_rows"].append({"species": name, **(metrics if quality["analysis_eligible"] else {key: None for key in METRIC_DEFINITIONS})})
                print(f"[gbif] {name}: status={quality['status']} fetched={quality['raw_fetched']} retained={quality['retained_count']}", flush=True)
        matched_keys = Counter(row["taxon_key"] for row in bundle["quality"] if row["taxon_key"] and row["status"] != "taxon_unresolved")
        for quality, analysis_row in zip(bundle["quality"], bundle["analysis_rows"], strict=True):
            if matched_keys[quality["taxon_key"]] > 1:
                quality.update(analysis_eligible=False, ineligibility_reason="shared_species_key")
                analysis_row.update({key: None for key in METRIC_DEFINITIONS})
        bundle["records_sha256"] = file_sha256(temporary)
        records_path = cache_path.with_name(cache_path.stem + "." + bundle["records_sha256"] + ".records.jsonl.gz")
        temporary.replace(records_path)
        bundle["records_path"] = str(records_path.resolve())
    finally:
        temporary.unlink(missing_ok=True)
    bundle["bundle_sha256"] = hashlib.sha256(json_text(bundle).encode()).hexdigest()
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=cache_path.parent, suffix=".cache.tmp", delete=False) as handle:
        cache_tmp = Path(handle.name)
        handle.write(json_text(bundle) + "\n")
    try:
        cache_tmp.replace(cache_path)
    finally:
        cache_tmp.unlink(missing_ok=True)
    if config["gbif_require_complete"] == "yes" and any(not row["analysis_eligible"] for row in bundle["quality"]):
        raise ValueError(f"GBIF acquisition is incomplete or unresolved; audit preserved in {cache_path}")
    frame = pandas.DataFrame.from_records(bundle["analysis_rows"])
    frame.attrs["gbif_bundle"] = bundle
    return frame
