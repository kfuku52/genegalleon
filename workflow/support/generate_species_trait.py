#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import sys
import tempfile
import zipfile
from dataclasses import dataclass
from io import StringIO, TextIOWrapper
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple
from urllib.parse import quote, unquote, urlparse
from urllib.request import Request

import pandas

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from format_species_network import guarded_urlopen as urlopen
from gbif_observations import DEFAULT_GBIF_DISTRIBUTION_TRAITS, METRIC_DEFINITIONS
from gbif_observations import effective_config as gbif_effective_config
from gbif_observations import fetch_gbif_distribution_table as _fetch_gbif_distribution_table
from gbif_observations import input_identity as gbif_input_identity
from gift_retrieval import GiftRetrieval, load_reviewed_mappings
from species_labeling import base_species_label, species_label_from_taxonomic_text
from species_trait_contract import sidecar_paths, trait_bundle_payloads
from species_trait_schema import schema_path, schema_payload

try:
    from openpyxl import load_workbook
except Exception:  # pragma: no cover - exercised when openpyxl is unavailable
    load_workbook = None


SUPPORTED_DATABASES = {
    "try": {
        "acquisition_mode": "bulk",
        "scope": "all_species_snapshot",
        "notes": "Use locally provided TRY export due to access restrictions.",
    },
    "gift": {
        "acquisition_mode": "gift_api",
        "scope": "trait_subset_api",
        "notes": "Resolve target work_IDs via species lookup, then fetch requested trait IDs.",
    },
    "gbif": {
        "acquisition_mode": "gbif_distribution",
        "scope": "target_species_occurrence_search",
        "notes": "Resolve species in the GBIF backbone and summarize no-login occurrence-search coordinates.",
    },
    "bien": {
        "acquisition_mode": "species_api",
        "scope": "target_species_only",
        "notes": "Species-level retrieval via API/package layer is expected.",
    },
    "eol_traitbank": {
        "acquisition_mode": "species_api",
        "scope": "target_species_only",
        "notes": "Species-level structured API retrieval is expected.",
    },
    "austraits": {
        "acquisition_mode": "bulk",
        "scope": "all_species_snapshot",
        "notes": "Download full release snapshot then subset target species.",
    },
    "eltontraits": {
        "acquisition_mode": "bulk",
        "scope": "all_species_snapshot",
        "notes": "Download full release snapshot then subset target species.",
    },
    "combine": {
        "acquisition_mode": "bulk",
        "scope": "all_species_snapshot",
        "notes": "Download full release snapshot then subset target species.",
    },
    "birdbase": {
        "acquisition_mode": "bulk",
        "scope": "all_species_snapshot",
        "notes": "Download full release snapshot then subset target species.",
    },
    "amniote": {
        "acquisition_mode": "bulk",
        "scope": "all_species_snapshot",
        "notes": "Download full release snapshot then subset target species.",
    },
    "amphibio": {
        "acquisition_mode": "bulk",
        "scope": "all_species_snapshot",
        "notes": "Download full release snapshot then subset target species.",
    },
    "animaltraits": {
        "acquisition_mode": "bulk",
        "scope": "all_species_snapshot",
        "notes": "Download full release snapshot then subset target species.",
    },
    "fishbase": {
        "acquisition_mode": "bulk",
        "scope": "all_species_snapshot",
        "notes": "Cache bulk tables/snapshot then subset target species.",
    },
}

FASTA_EXTENSIONS = (
    ".fa",
    ".fas",
    ".fasta",
    ".fna",
    ".fa.gz",
    ".fas.gz",
    ".fasta.gz",
    ".fna.gz",
)
SPECIES_COLUMN_CANDIDATES = (
    "species",
    "species_name",
    "scientific_name",
    "scientificname",
    "binomial",
    "taxon_name",
    "taxon",
)
DEFAULT_MANIFEST_PATH = Path("workspace/input/input_generation/download_plan.xlsx")
DEFAULT_TRAIT_PLAN_PATH = Path("workspace/input/input_generation/trait_plan.tsv")
DEFAULT_DB_SOURCES_PATH = Path("workspace/input/input_generation/trait_database_sources.tsv")
DEFAULT_DOWNLOADS_DIR = Path("workspace/downloads/trait_datasets")
DEFAULT_OUTPUT_PATH = Path("workspace/input/species_trait/species_trait.tsv")
DEFAULT_GIFT_API = "https://gift.uni-goettingen.de/api/extended/"
DEFAULT_GIFT_PAGE_SIZE = 10000
GIFT_TRAIT_ID_PATTERN = re.compile(r"^\d+(?:\.\d+)+$")



@dataclass
class TraitPlanRow:
    database: str
    source_column: str
    output_trait: str
    value_type: str
    aggregation: str
    positive_values: Set[str]
    trait_key: str
    trait_key_column: str


def _log(message: str) -> None:
    print(message, flush=True)


def builtin_gbif_trait_plan_rows() -> List[TraitPlanRow]:
    return [
        TraitPlanRow(
            database="gbif",
            source_column=source_column,
            output_trait=source_column,
            value_type=value_type,
            aggregation="any" if value_type == "binary" else "median",
            positive_values={"1"} if value_type == "binary" else set(),
            trait_key="",
            trait_key_column="",
        )
        for source_column, value_type in DEFAULT_GBIF_DISTRIBUTION_TRAITS
    ]


def add_builtin_trait_plan_rows(
    plan_rows: Sequence[TraitPlanRow],
    requested_databases: Sequence[str],
) -> List[TraitPlanRow]:
    out = list(plan_rows)
    requested = {database.strip().lower() for database in requested_databases}
    if "gbif" in requested and not any(row.database == "gbif" for row in out):
        out.extend(builtin_gbif_trait_plan_rows())
    return out


def normalize_species_name(value: object) -> str:
    raw = str(value or "").strip()
    if raw == "":
        return ""
    return species_label_from_taxonomic_text(raw)


def parse_species_from_id_label(value: str) -> str:
    text = str(value or "").strip()
    if text == "":
        return ""
    match = re.search(r"\(([^()]+)\)", text)
    if match:
        normalized = normalize_species_name(match.group(1))
        if normalized != "":
            return normalized
    return normalize_species_name(text)


def parse_delimiter(path: Path, explicit: str) -> str:
    explicit_norm = str(explicit or "").strip().lower()
    if explicit_norm in ("tab", "\\t", "tsv"):
        return "\t"
    if explicit_norm in ("comma", ",", "csv"):
        return ","
    if explicit_norm == "pipe":
        return "|"
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith(".tsv") or suffixes.endswith(".tsv.gz"):
        return "\t"
    if suffixes.endswith(".csv") or suffixes.endswith(".csv.gz"):
        return ","
    return "\t"


def strip_known_extension(path: Path) -> str:
    name = path.name
    lowered = name.lower()
    for ext in FASTA_EXTENSIONS:
        if lowered.endswith(ext):
            return name[: len(name) - len(ext)]
    return path.stem


def read_manifest_rows(path: Path) -> List[Dict[str, str]]:
    suffix = "".join(path.suffixes).lower()
    if suffix.endswith(".xlsx"):
        if load_workbook is None:
            raise RuntimeError("openpyxl is required to parse XLSX manifest files.")
        workbook = load_workbook(path, read_only=True, data_only=True)
        sheet = workbook.active
        row_iter = sheet.iter_rows(values_only=True)
        header_values = next(row_iter, None)
        if header_values is None:
            return []
        header = [str(cell or "").strip() for cell in header_values]
        rows: List[Dict[str, str]] = []
        for cells in row_iter:
            if cells is None:
                continue
            row = {header[i]: str(cells[i] or "").strip() for i in range(min(len(header), len(cells)))}
            if any(value.strip() != "" for value in row.values()):
                rows.append(row)
        return rows
    delimiter = parse_delimiter(path, "")
    with open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return [{key: str(value or "").strip() for key, value in row.items()} for row in reader]


def species_from_manifest(path: Path) -> Set[str]:
    rows = read_manifest_rows(path)
    species: Set[str] = set()
    for row in rows:
        provider = str(row.get("provider", "") or "").strip()
        if provider.startswith("#"):
            continue
        candidate = normalize_species_name(row.get("species_key", ""))
        if candidate == "":
            candidate = parse_species_from_id_label(row.get("id", ""))
        if candidate != "":
            species.add(candidate)
    return species


def species_from_cds_dir(path: Path) -> Set[str]:
    species: Set[str] = set()
    for file_path in sorted(path.iterdir()):
        if not file_path.is_file():
            continue
        if file_path.name.startswith("."):
            continue
        stem = strip_known_extension(file_path)
        candidate = normalize_species_name(stem)
        if candidate != "":
            species.add(candidate)
    return species


def read_config_rows(path: Path, required: Set[str]):
    with path.open("rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        columns = reader.fieldnames or []
        if len(columns) != len(set(columns)) or any(not col or col != col.strip() for col in columns):
            raise ValueError("Duplicate, empty or padded columns in " + str(path))
        if not required.issubset(columns):
            raise ValueError("Missing required columns in {}: {}".format(path, sorted(required - set(columns))))
        rows = []
        for row in reader:
            if None in row:
                raise ValueError("Too many fields in " + str(path))
            if any(not str(row.get(key, "") or "").strip() for key in required):
                raise ValueError("Empty required field in " + str(path))
            rows.append({key: str(value or "").strip() for key, value in row.items()})
        return rows


def validate_trait_plan_row(row):
    allowed = {"numeric": {"median", "mean", "min", "max"},
               "binary": {"any", "all", "min", "max", "sum", "mean"},
               "categorical": {"first", "mode"}, "text": {"unique", "first"}}
    if row.value_type not in allowed or row.aggregation not in allowed[row.value_type]:
        raise ValueError("Unsupported trait type/aggregation: {}/{}".format(row.value_type, row.aggregation))
    if any(char in row.output_trait for char in "\t\r\n"):
        raise ValueError("Output trait names must not contain tabs or newlines")
    if row.output_trait in {"species", "__species_norm"}:
        raise ValueError("Reserved output trait name: " + row.output_trait)
    if row.positive_values and row.value_type != "binary":
        raise ValueError("positive_values requires a binary trait")


def read_trait_plan(path: Path) -> List[TraitPlanRow]:
    rows = read_config_rows(path, {"database", "source_column", "output_trait"})
    out: List[TraitPlanRow] = []
    for raw in rows:
        database = str(raw.get("database", "") or "").strip().lower()
        source_column = str(raw.get("source_column", "") or "").strip()
        output_trait = str(raw.get("output_trait", "") or "").strip()
        if database == "" or source_column == "" or output_trait == "":
            continue
        value_type = str(raw.get("value_type", "numeric") or "numeric").strip().lower()
        aggregation = str(raw.get("aggregation", "") or "").strip().lower()
        if aggregation == "":
            aggregation = {"binary": "any", "categorical": "mode", "text": "unique"}.get(value_type, "median")
        positive_raw = str(raw.get("positive_values", "") or "").strip()
        positive_values = {token.strip().lower() for token in positive_raw.split(",") if token.strip() != ""}
        trait_key = str(raw.get("trait_key", "") or "").strip()
        trait_key_column = str(raw.get("trait_key_column", "") or "").strip()
        out.append(
            TraitPlanRow(
                database=database,
                source_column=source_column,
                output_trait=output_trait,
                value_type=value_type,
                aggregation=aggregation,
                positive_values=positive_values,
                trait_key=trait_key,
                trait_key_column=trait_key_column,
            )
        )
    for row in out:
        validate_trait_plan_row(row)
    return out


def read_database_sources(path: Path) -> Dict[str, Dict[str, str]]:
    rows = read_config_rows(path, {"database"})
    out: Dict[str, Dict[str, str]] = {}
    for raw in rows:
        database = str(raw.get("database", "") or "").strip().lower()
        if database == "":
            continue
        if database in out:
            raise ValueError("Duplicate database source: " + database)
        normalized = {str(k or "").strip(): str(v or "").strip() for k, v in raw.items()}
        uri = normalized.get("uri", "")
        if uri != "":
            parsed = urlparse(uri)
            if parsed.scheme == "" and not Path(uri).is_absolute():
                normalized["uri"] = str((path.parent / uri).resolve())
        mapping_path = normalized.get("gift_species_mapping_file", "")
        if mapping_path:
            candidate = Path(mapping_path).expanduser()
            normalized["gift_species_mapping_file"] = str((path.parent / candidate).resolve() if not candidate.is_absolute() else candidate.resolve())
        out[database] = normalized
    return out


def resolve_requested_databases(
    databases_arg: str,
    plan_rows: Sequence[TraitPlanRow],
    source_rows: Dict[str, Dict[str, str]],
) -> List[str]:
    requested = str(databases_arg or "").strip().lower()
    if requested == "" or requested == "auto":
        return sorted({row.database for row in plan_rows})
    if requested == "all":
        return sorted(SUPPORTED_DATABASES.keys())
    dbs = [token.strip().lower() for token in requested.split(",") if token.strip() != ""]
    if len(dbs) == 0:
        return sorted({row.database for row in plan_rows})
    unknown = [db for db in dbs if db not in SUPPORTED_DATABASES and db not in source_rows]
    if unknown:
        raise ValueError("Unknown database(s): {}".format(", ".join(sorted(unknown))))
    return dbs


def detect_species_column(df: pandas.DataFrame, requested: str) -> str:
    requested_norm = str(requested or "").strip()
    if requested_norm != "" and requested_norm in df.columns:
        return requested_norm
    lowered = {col.lower(): col for col in df.columns}
    for candidate in SPECIES_COLUMN_CANDIDATES:
        if candidate in lowered:
            return lowered[candidate]
    return str(df.columns[0])


def read_table(path: Path, delimiter: str) -> pandas.DataFrame:
    suffix = "".join(path.suffixes).lower()
    if suffix.endswith(".xlsx"):
        return pandas.read_excel(path, dtype=str, keep_default_na=False)
    sep = parse_delimiter(path, delimiter)
    return pandas.read_csv(path, sep=sep, dtype=str, keep_default_na=False)


def read_table_from_text(text: str, delimiter: str) -> pandas.DataFrame:
    sep = parse_delimiter(Path("response.tsv"), delimiter)
    return pandas.read_csv(StringIO(text), sep=sep, dtype=str, keep_default_na=False)


def split_uri_list(uri_raw: str) -> List[str]:
    return [token.strip() for token in str(uri_raw or "").split(",") if token.strip() != ""]


def choose_archive_member(members: Sequence[str], requested: str) -> Optional[str]:
    if len(members) == 0:
        return None
    requested_norm = str(requested or "").strip()
    if requested_norm != "":
        if requested_norm in members:
            return requested_norm
        return None
    preferred_suffixes = (".tsv", ".csv", ".txt", ".xlsx")
    for member in members:
        lowered = member.lower()
        for suffix in preferred_suffixes:
            if lowered.endswith(suffix):
                return member
    return members[0]


def read_table_from_zip(path: Path, delimiter: str, archive_member: str) -> pandas.DataFrame:
    with zipfile.ZipFile(path) as archive:
        members = [name for name in archive.namelist() if not name.endswith("/")]
        selected_member = choose_archive_member(members=members, requested=archive_member)
        if selected_member is None:
            raise ValueError("archive_member not found in {}: {}".format(path, archive_member))
        if selected_member.lower().endswith(".xlsx"):
            with (
                archive.open(selected_member) as source,
                tempfile.SpooledTemporaryFile(
                    max_size=8 * 1024 * 1024,
                ) as temporary,
            ):
                shutil.copyfileobj(source, temporary, length=1024 * 1024)
                temporary.seek(0)
                return pandas.read_excel(temporary, dtype=str)
        sep = parse_delimiter(Path(selected_member), delimiter)
        for encoding in ("utf-8", "latin1"):
            try:
                with archive.open(selected_member) as source:
                    with TextIOWrapper(source, encoding=encoding) as text:
                        return pandas.read_csv(text, sep=sep, dtype=str)
            except UnicodeDecodeError:
                continue
    raise UnicodeDecodeError("utf-8", b"", 0, 1, "failed to decode ZIP member")


def read_bulk_table_from_path(path: Path, config: Dict[str, str]) -> pandas.DataFrame:
    delimiter = config.get("delimiter", "")
    archive_member = str(config.get("archive_member", "") or "").strip()
    suffix = "".join(path.suffixes).lower()
    if suffix.endswith(".zip"):
        return read_table_from_zip(path=path, delimiter=delimiter, archive_member=archive_member)
    return read_table(path=path, delimiter=delimiter)


def copy_or_download_file(
    uri: str,
    destination: Path,
    timeout: float,
    dry_run: bool,
) -> Optional[Path]:
    if uri == "":
        return None
    parsed = urlparse(uri)
    if parsed.scheme in ("",):
        src = Path(uri).expanduser().resolve()
        if not src.exists():
            return None
        if src == destination:
            return src
        if dry_run:
            _log("[dry-run] copy {} -> {}".format(src, destination))
            return destination
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, destination)
        return destination
    if parsed.scheme == "file":
        src = Path(unquote(parsed.path)).resolve()
        if not src.exists():
            return None
        if src == destination:
            return src
        if dry_run:
            _log("[dry-run] copy {} -> {}".format(src, destination))
            return destination
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, destination)
        return destination
    if parsed.scheme in ("http", "https", "ftp"):
        if dry_run:
            _log("[dry-run] download {} -> {}".format(uri, destination))
            return destination
        destination.parent.mkdir(parents=True, exist_ok=True)
        request = Request(uri, headers={"User-Agent": "genegalleon-trait-generator"})
        with urlopen(request, timeout=timeout) as response, open(destination, "wb") as handle:
            shutil.copyfileobj(response, handle)
        return destination
    raise ValueError("Unsupported URI scheme in: {}".format(uri))


def detect_existing_bulk_file(database: str, downloads_dir: Path) -> Optional[Path]:
    db_dir = downloads_dir / database
    candidates = [
        db_dir / "{}.tsv".format(database),
        db_dir / "{}.csv".format(database),
        db_dir / "{}.txt".format(database),
        db_dir / "{}.zip".format(database),
        db_dir / "{}.tsv.gz".format(database),
        db_dir / "{}.csv.gz".format(database),
        db_dir / "source.tsv",
        db_dir / "source.csv",
        db_dir / "source.txt",
        db_dir / "source.zip",
        db_dir / "source.tsv.gz",
        db_dir / "source.csv.gz",
        db_dir / "source.xlsx",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def load_bulk_database(
    database: str,
    config: Dict[str, str],
    downloads_dir: Path,
    timeout: float,
    dry_run: bool,
) -> Optional[pandas.DataFrame]:
    uri_values = split_uri_list(config.get("uri", ""))
    db_dir = downloads_dir / database
    table_paths: List[Path] = []

    if len(uri_values) > 0:
        for index, uri in enumerate(uri_values):
            parsed = urlparse(uri)
            basename = Path(parsed.path).name or "{}.{}.tsv".format(database, index + 1)
            target = db_dir / basename
            source_path = copy_or_download_file(uri=uri, destination=target, timeout=timeout, dry_run=dry_run)
            if source_path is None:
                continue
            if not dry_run:
                table_paths.append(source_path)
        if dry_run:
            return None
    else:
        existing = detect_existing_bulk_file(database=database, downloads_dir=downloads_dir)
        if existing is None:
            return None
        table_paths.append(existing)

    if len(table_paths) == 0:
        return None
    frames: List[pandas.DataFrame] = []
    for table_path in table_paths:
        frames.append(read_bulk_table_from_path(path=table_path, config=config))
    if len(frames) == 1:
        return frames[0]
    return pandas.concat(frames, ignore_index=True)


def fetch_species_api_table(
    database: str,
    config: Dict[str, str],
    species: Sequence[str],
    timeout: float,
    dry_run: bool,
) -> Optional[pandas.DataFrame]:
    uri_template = str(config.get("uri", "") or "").strip()
    if uri_template == "":
        return None
    response_format = str(config.get("response_format", "tsv") or "tsv").strip().lower()
    delimiter = config.get("delimiter", "\t")
    frames: List[pandas.DataFrame] = []
    for species_name in species:
        safe_species = quote(species_name.replace("_", " "))
        url = (
            uri_template.replace("{species}", species_name)
            .replace("{species_space}", species_name.replace("_", " "))
            .replace("{species_urlencoded}", safe_species)
        )
        if dry_run:
            _log("[dry-run] {} request: {}".format(database, url))
            continue
        request = Request(url, headers={"User-Agent": "genegalleon-trait-generator"})
        with urlopen(request, timeout=timeout) as response:
            payload = response.read().decode("utf-8")
        if response_format in ("tsv", "csv"):
            sep = "\t" if response_format == "tsv" else ","
            if str(delimiter or "").strip() != "":
                sep = parse_delimiter(Path("dummy.{}".format(response_format)), delimiter)
            frames.append(read_table_from_text(payload, sep))
        else:
            raise ValueError("Unsupported species_api response_format '{}' for '{}'.".format(response_format, database))
    if len(frames) == 0:
        return None
    return pandas.concat(frames, ignore_index=True)


def parse_bool_option(value: object, key_name: str, default: bool) -> bool:
    text = str(value or "").strip().lower()
    if text == "":
        return default
    if text in ("1", "true", "yes", "y", "on"):
        return True
    if text in ("0", "false", "no", "n", "off"):
        return False
    raise ValueError("Invalid boolean value for {}: {}".format(key_name, value))


def parse_positive_int_option(value: object, key_name: str, default: int) -> int:
    text = str(value or "").strip()
    if text == "":
        return default
    parsed = int(text)
    if parsed <= 0:
        raise ValueError("{} must be a positive integer: {}".format(key_name, value))
    return parsed


def parse_optional_float_option(value: object, key_name: str) -> Optional[float]:
    text = str(value or "").strip()
    if text == "":
        return None
    parsed = float(text)
    if math.isnan(parsed):
        raise ValueError("{} must not be NaN.".format(key_name))
    return parsed


def split_csv_tokens(value: object) -> List[str]:
    return [token.strip() for token in str(value or "").split(",") if token.strip() != ""]


def normalize_base_uri(uri: str, fallback: str) -> str:
    base = str(uri or "").strip()
    if base == "":
        base = fallback
    if not base.endswith("/"):
        base = base + "/"
    return base


def fetch_json_payload(url: str, timeout: float) -> object:
    request = Request(
        url,
        headers={
            "User-Agent": "genegalleon-trait-generator",
            "Accept": "application/json",
        },
    )
    with urlopen(request, timeout=timeout) as response:
        payload = response.read().decode("utf-8")
    try:
        return json.loads(payload)
    except json.JSONDecodeError as exc:
        raise ValueError("Failed to parse JSON response from {}: {}".format(url, exc)) from exc


def json_payload_to_rows(payload: object) -> List[Dict[str, object]]:
    if payload is None:
        return []
    if isinstance(payload, list):
        return [row for row in payload if isinstance(row, dict)]
    if isinstance(payload, dict):
        for candidate_key in ("data", "results", "items"):
            value = payload.get(candidate_key)
            if isinstance(value, list):
                return [row for row in value if isinstance(row, dict)]
        return [payload]
    return []


def resolve_gift_version(
    data_api_base: str,
    requested_version: str,
    timeout: float,
    versions_api: str,
    fetcher=None,
) -> str:
    version = str(requested_version or "").strip()
    if version == "":
        version = "latest"
    if version != "latest":
        return version
    if str(versions_api or "").strip() != "":
        versions_base = normalize_base_uri(versions_api, DEFAULT_GIFT_API)
    elif "/extended/" in data_api_base:
        versions_base = normalize_base_uri(data_api_base.replace("/extended/", "/", 1), DEFAULT_GIFT_API)
    else:
        versions_base = data_api_base
    versions_url = "{}index.php?query=versions".format(versions_base)
    rows = json_payload_to_rows((fetcher or fetch_json_payload)(url=versions_url, timeout=timeout))
    if len(rows) == 0:
        raise ValueError("GIFT versions endpoint returned no rows: {}".format(versions_url))
    resolved = str(rows[-1].get("version", "") or "").strip()
    if resolved == "":
        raise ValueError("GIFT versions endpoint did not include 'version': {}".format(versions_url))
    return resolved


def build_gift_index_url(data_api_base: str, version: str) -> str:
    if version == "beta":
        index_name = "index.php"
    else:
        index_name = "index{}.php".format(version)
    return "{}{}".format(data_api_base, index_name)


def split_genus_epithet(species_name: str) -> Optional[Sequence[str]]:
    base_label = base_species_label(species_name)
    tokens = [token for token in str(base_label or "").strip().split("_") if token]
    if len(tokens) < 2:
        return None
    if tokens[1].lower() == "sp":
        return None
    return (tokens[0], tokens[1])


def remap_species_name_to_target(
    species_name: str,
    target_species: Set[str],
    target_species_by_base: Dict[str, Set[str]],
) -> str:
    normalized = str(species_name or "").strip()
    if normalized == "":
        return ""
    if normalized in target_species:
        return normalized
    base_label = base_species_label(normalized)
    if "_" not in base_label:
        return ""
    candidates = target_species_by_base.get(base_label, set())
    if len(candidates) == 1:
        return next(iter(candidates))
    return ""


def resolve_gift_species_map(
    index_url: str,
    species: Sequence[str],
    timeout: float,
    fetcher=None,
    reviewed_mappings=None,
    report=None,
) -> pandas.DataFrame:
    records = []
    mappings = reviewed_mappings or {}
    fetch = fetcher or fetch_json_payload
    for species_name in species:
        target = base_species_label(species_name)
        reviewed = mappings.get(target)
        audit = {"species": species_name, "status": "unmatched", "candidates": []}
        if reviewed:
            audit["reviewed_mapping"] = reviewed
        if normalize_species_name(species_name) != target:
            audit["status"] = "qualified_taxon_requires_review"
        elif reviewed and reviewed["decision"] == "exclude":
            audit["status"] = "excluded_taxonomic_scope"
        else:
            lookup_name = normalize_species_name(reviewed["work_species"]) if reviewed else target
            genus_epithet = split_genus_epithet(lookup_name)
            if genus_epithet:
                lookup_url = "{}?query=names_matched_unique&genus={}&epithet={}".format(
                    index_url, quote(genus_epithet[0]), quote(genus_epithet[1]))
                rows = json_payload_to_rows(fetch(url=lookup_url, timeout=timeout))
                candidates = {(str(row.get("work_ID", "") or "").strip(),
                               str(row.get("work_species", "") or "").strip()) for row in rows}
                audit["candidates"] = [{"work_ID": key, "work_species": name} for key, name in sorted(candidates)]
                # Never select an arbitrary high-score synonym or an infraspecific hit.
                eligible = {(str(row.get("work_ID", "") or "").strip(), str(row.get("work_species", "") or "").strip())
                            for row in rows if all(str(row.get(flag, "1")) == "1" for flag in ("matched", "resolved", "accepted"))}
                exact = {(key, name) for key, name in eligible
                         if key.isdigit() and int(key) > 0 and normalize_species_name(name) == lookup_name}
                if reviewed:
                    exact = {(key, name) for key, name in exact if key == reviewed["work_ID"]}
                if len(exact) == 1:
                    key, name = next(iter(exact))
                    record = {"species": species_name, "work_ID": key, "work_species": name}
                    records.append(record)
                    audit.update(record, status="reviewed_synonym" if reviewed else "exact")
                elif reviewed:
                    audit["status"] = "reviewed_mapping_mismatch"
                elif candidates:
                    audit["status"] = "ambiguous_or_synonym_requires_review"
        if audit["status"] not in {"exact", "reviewed_synonym"}:
            _log("WARNING: [gift] {}: {}".format(species_name, audit["status"]))
        if report is not None:
            report.append(audit)
    return pandas.DataFrame(records, columns=["species", "work_ID", "work_species"])


def dedupe_keep_order(values: Sequence[str]) -> List[str]:
    deduped: List[str] = []
    seen: Set[str] = set()
    for value in values:
        key = str(value or "").strip()
        if key == "" or key in seen:
            continue
        seen.add(key)
        deduped.append(key)
    return deduped


def collect_gift_trait_tokens(
    database: str,
    config: Dict[str, str],
    plan_rows: Sequence[TraitPlanRow],
) -> List[str]:
    out: List[str] = []
    config_trait_ids = str(config.get("gift_trait_ids", "") or config.get("trait_ids", "") or "").strip()
    if config_trait_ids != "":
        out.extend([token.strip() for token in config_trait_ids.split(",") if token.strip() != ""])
    for row in plan_rows:
        if row.database != database:
            continue
        trait_id = str(row.trait_key or "").strip()
        if trait_id != "":
            out.append(trait_id)
    return dedupe_keep_order(out)


def is_gift_trait_id(value: str) -> bool:
    return bool(GIFT_TRAIT_ID_PATTERN.fullmatch(str(value or "").strip()))


def parse_float_or_default(value: object, default: float) -> float:
    try:
        return float(str(value or "").strip())
    except Exception:
        return default


def strip_string_series(series: pandas.Series, lower: bool = False) -> pandas.Series:
    stripped = series.fillna("").astype(str).str.strip()
    if lower:
        stripped = stripped.str.lower()
    return stripped


def aggregate_categorical_mode(
    values: pandas.Series,
    group_keys: pandas.Series,
) -> pandas.Series:
    valid_mask = values != ""
    if not valid_mask.any():
        return pandas.Series(dtype=object)
    counts = (
        pandas.DataFrame(
            {
                "__species_norm": group_keys.loc[valid_mask].to_numpy(),
                "__value": values.loc[valid_mask].to_numpy(),
            }
        )
        .groupby(["__species_norm", "__value"], observed=True, sort=False)
        .size()
        .rename("__count")
        .reset_index()
    )
    counts = counts.sort_values(
        ["__species_norm", "__count", "__value"],
        ascending=[True, False, True],
        kind="mergesort",
    )
    return counts.drop_duplicates(subset="__species_norm", keep="first").set_index("__species_norm")["__value"]


def fetch_gift_traits_meta_rows(index_url: str, timeout: float, fetcher=None) -> List[Dict[str, object]]:
    traits_meta_url = "{}?query=traits_meta".format(index_url)
    rows = json_payload_to_rows((fetcher or fetch_json_payload)(url=traits_meta_url, timeout=timeout))
    return [row for row in rows if isinstance(row, dict)]


def resolve_gift_trait_token_map(
    index_url: str,
    trait_tokens: Sequence[str],
    timeout: float,
    fetcher=None,
) -> Dict[str, str]:
    tokens = dedupe_keep_order([str(token or "").strip() for token in trait_tokens])
    if len(tokens) == 0:
        return {}
    resolved_map: Dict[str, str] = {}
    direct_ids = [token for token in tokens if is_gift_trait_id(token)]
    name_tokens = [token for token in tokens if not is_gift_trait_id(token)]
    for token in direct_ids:
        resolved_map[token] = token
    if len(name_tokens) == 0:
        return resolved_map

    traits_meta_rows = fetch_gift_traits_meta_rows(index_url=index_url, timeout=timeout, fetcher=fetcher)
    for token in name_tokens:
        token_norm = token.casefold()
        candidates = [row for row in traits_meta_rows if str(row.get("Trait2", "")).strip().casefold() == token_norm]
        if not candidates:
            candidates = [row for row in traits_meta_rows if str(row.get("Trait1", "")).strip().casefold() == token_norm]
        identities = {str(row.get("Lvl3", "")).strip() for row in candidates}
        if not identities:
            raise ValueError("Unknown GIFT trait name '{}'; use an exact Trait2 name or trait ID".format(token))
        if len(identities) != 1 or not all(is_gift_trait_id(value) for value in identities):
            raise ValueError("GIFT trait name '{}' is ambiguous; choose one trait ID from {}".format(token, sorted(identities)))
        resolved_map[token] = next(iter(identities))
    return resolved_map


def print_gift_traits(
    data_api_base: str,
    requested_version: str,
    versions_api: str,
    timeout: float,
    search: str,
    limit: int,
) -> None:
    resolved_version = resolve_gift_version(
        data_api_base=data_api_base,
        requested_version=requested_version,
        timeout=timeout,
        versions_api=versions_api,
    )
    index_url = build_gift_index_url(data_api_base=data_api_base, version=resolved_version)
    rows = fetch_gift_traits_meta_rows(index_url=index_url, timeout=timeout)
    if len(rows) == 0:
        print("trait_id\ttrait_name\ttrait_group\tvalue_type\tunits\tcount")
        return

    search_norm = str(search or "").strip().lower()
    filtered: List[Dict[str, object]] = []
    for row in rows:
        row_text = "\t".join(
            [
                str(row.get("Lvl3", "") or ""),
                str(row.get("Trait2", "") or ""),
                str(row.get("Trait1", "") or ""),
                str(row.get("Category", "") or ""),
            ]
        ).lower()
        if search_norm != "" and search_norm not in row_text:
            continue
        filtered.append(row)
    sorted_rows = sorted(
        filtered,
        key=lambda row: (
            parse_float_or_default(row.get("count", 0), 0.0),
            str(row.get("Lvl3", "") or ""),
        ),
        reverse=True,
    )
    if limit > 0:
        sorted_rows = sorted_rows[:limit]

    print("trait_id\ttrait_name\ttrait_group\tvalue_type\tunits\tcount")
    for row in sorted_rows:
        print(
            "{}\t{}\t{}\t{}\t{}\t{}".format(
                str(row.get("Lvl3", "") or "").strip(),
                str(row.get("Trait2", "") or "").strip(),
                str(row.get("Trait1", "") or "").strip(),
                str(row.get("type", "") or "").strip(),
                str(row.get("Units", "") or "").strip(),
                str(row.get("count", "") or "").strip(),
            )
        )


def fetch_gift_api_table(
    database: str,
    config: Dict[str, str],
    plan_rows: Sequence[TraitPlanRow],
    species: Sequence[str],
    timeout: float,
    dry_run: bool,
    downloads_dir: Optional[Path] = None,
) -> Optional[pandas.DataFrame]:
    data_api_base = normalize_base_uri(config.get("uri", ""), DEFAULT_GIFT_API)
    requested_version = str(config.get("gift_version", "latest") or "latest").strip().lower()
    trait_tokens = collect_gift_trait_tokens(database=database, config=config, plan_rows=plan_rows)
    if len(trait_tokens) == 0:
        _log(
            "WARNING: [gift] no trait IDs resolved. Set trait_key in trait_plan or gift_trait_ids in database_sources."
        )
        return None
    page_size = parse_positive_int_option(
        value=config.get("gift_page_size", ""),
        key_name="gift_page_size",
        default=DEFAULT_GIFT_PAGE_SIZE,
    )
    max_pages = 0
    max_pages_raw = str(config.get("gift_max_pages_per_trait", "") or "").strip()
    if max_pages_raw != "":
        max_pages = parse_positive_int_option(
            value=max_pages_raw,
            key_name="gift_max_pages_per_trait",
            default=0,
        )
    if page_size > 10000:
        raise ValueError("gift_page_size must be between 1 and 10000")
    bias_ref = parse_bool_option(config.get("gift_bias_ref", ""), key_name="gift_bias_ref", default=True)
    bias_deriv = parse_bool_option(config.get("gift_bias_deriv", ""), key_name="gift_bias_deriv", default=True)
    agreement_min: Optional[float] = None
    agreement_min_raw = str(config.get("gift_agreement_min", "") or "").strip()
    if agreement_min_raw != "":
        agreement_min = float(agreement_min_raw)
        if not math.isfinite(agreement_min) or not 0 <= agreement_min <= 1:
            raise ValueError("gift_agreement_min must be finite and between 0 and 1")

    mode = str(config.get("gift_cache_mode", "reuse") or "reuse").strip().lower()
    retries = int(config.get("gift_retries", "2") or "2")
    if mode not in {"reuse", "refresh", "offline"} or not 0 <= retries <= 5:
        raise ValueError("Invalid gift_cache_mode or gift_retries")

    if dry_run:
        display_version = requested_version if requested_version != "latest" else "<latest>"
        display_index_url = build_gift_index_url(data_api_base=data_api_base, version=display_version)
        _log("[dry-run] GIFT names below are candidates; reviewed mappings are applied after release resolution.")
        for species_name in species:
            genus_epithet = split_genus_epithet(species_name)
            if genus_epithet is None:
                continue
            genus, epithet = genus_epithet[0], genus_epithet[1]
            _log(
                "[dry-run] {} request: {}?query=names_matched_unique&genus={}&epithet={}".format(
                    database,
                    display_index_url,
                    quote(genus),
                    quote(epithet),
                )
            )
        for trait_token in trait_tokens:
            _log(
                "[dry-run] {} request: {}?query=traits&traitid={}&biasref={}&biasderiv={}&startat=0&limit={}".format(
                    database,
                    display_index_url,
                    quote(trait_token),
                    int(bias_ref),
                    int(bias_deriv),
                    page_size,
                )
            )
        return None

    directory = downloads_dir / "gift" if downloads_dir is not None else None
    with GiftRetrieval(directory, fetch_json_payload, mode=mode, retries=retries, logger=_log) as client:
        return _fetch_gift_api_table(database, config, plan_rows, species, timeout, data_api_base,
                                     requested_version, trait_tokens, page_size, max_pages,
                                     bias_ref, bias_deriv, agreement_min, client)


def _fetch_gift_api_table(database, config, plan_rows, species, timeout, data_api_base,
                          requested_version, trait_tokens, page_size, max_pages,
                          bias_ref, bias_deriv, agreement_min, client):
    resolved_version = resolve_gift_version(
        data_api_base=data_api_base,
        requested_version=requested_version,
        timeout=timeout,
        versions_api=str(config.get("gift_versions_api", "") or "").strip(),
        fetcher=client.fetch,
    )
    index_url = build_gift_index_url(data_api_base=data_api_base, version=resolved_version)
    trait_token_map = resolve_gift_trait_token_map(
        index_url=index_url,
        trait_tokens=trait_tokens,
        timeout=timeout,
        fetcher=client.fetch,
    )
    trait_ids = dedupe_keep_order(list(trait_token_map.values()))
    if len(trait_ids) == 0:
        _log("WARNING: [gift] no trait IDs resolved after applying traits_meta lookup.")
        return None
    mappings, mapping_receipts = load_reviewed_mappings(
        resolved_version, str(config.get("gift_species_mapping_file", "") or ""),
        include_bundled=urlparse(index_url).hostname == "gift.uni-goettingen.de")
    client.report.update(resolved_version=resolved_version, mapping_sources=mapping_receipts,
                         requested_species=list(species), trait_token_map=trait_token_map,
                         agreement_min=agreement_min)
    species_map = resolve_gift_species_map(index_url=index_url, species=species, timeout=timeout,
        fetcher=client.fetch, reviewed_mappings=mappings, report=client.report["species"])
    if species_map.shape[0] == 0:
        return pandas.DataFrame(columns=["species", "work_ID", "work_species", "trait_ID", "trait_value"])
    work_ids = set(species_map["work_ID"].astype(str).tolist())
    required_columns: Set[str] = {"agreement"}
    default_trait_key_column = str(config.get("trait_key_column", "") or "").strip()
    for plan_row in plan_rows:
        if plan_row.database != database:
            continue
        required_columns.add(plan_row.source_column)
        if plan_row.trait_key == "":
            continue
        trait_key_column = plan_row.trait_key_column or default_trait_key_column or "trait_name"
        if trait_key_column not in ("trait_ID", "trait_token"):
            required_columns.add(trait_key_column)

    frames: List[pandas.DataFrame] = []
    for trait_id in trait_ids:
        trait_token = next(token for token in trait_tokens if trait_token_map[token] == trait_id)
        page_index = 0
        start_at = 0
        trait_row_count = 0
        seen_pages = set()
        while True:
            traits_url = "{}?query=traits&traitid={}&biasref={}&biasderiv={}&startat={}&limit={}".format(
                index_url,
                quote(trait_id),
                int(bias_ref),
                int(bias_deriv),
                start_at,
                page_size,
            )
            rows = json_payload_to_rows(client.fetch(url=traits_url, timeout=timeout))
            if len(rows) == 0:
                break
            page_signature = hashlib.sha256(json.dumps(rows, sort_keys=True).encode()).hexdigest()
            if page_signature in seen_pages:
                raise ValueError("GIFT returned a repeated page for trait " + trait_id)
            seen_pages.add(page_signature)
            matched_rows: List[Dict[str, object]] = []
            for row in rows:
                if not isinstance(row, dict):
                    continue
                work_id = str(row.get("work_ID", "") or "").strip()
                if work_id == "" or work_id not in work_ids:
                    continue
                matched_row: Dict[str, object] = {"work_ID": work_id, "trait_ID": trait_id, "trait_token": trait_token}
                for column in sorted(required_columns):
                    matched_row[column] = row.get(column, "")
                matched_rows.append(matched_row)
            if len(matched_rows) > 0:
                frames.append(pandas.DataFrame.from_records(matched_rows))
                trait_row_count += len(matched_rows)
            page_index += 1
            start_at += len(rows)
            if len(rows) < page_size:
                break
            if max_pages > 0 and page_index >= max_pages:
                raise ValueError("GIFT acquisition incomplete: gift_max_pages_per_trait={} for trait {}".format(max_pages, trait_id))
        client.report["traits"].append({"trait_ID": trait_id, "trait_token": trait_token,
                                        "pages": page_index, "matched_rows": trait_row_count})

    if len(frames) == 0:
        return pandas.DataFrame(columns=["species", "work_ID", "work_species", "trait_ID", "trait_value"])

    merged = pandas.concat(frames, ignore_index=True)
    # Multiple input labels can name the same taxon; retain all instead of keeping the first.
    merged = merged.merge(species_map, on="work_ID", how="inner", validate="many_to_many")

    if agreement_min is not None and "agreement" in merged.columns:
        agreement_numeric = pandas.to_numeric(merged["agreement"], errors="coerce")
        merged = merged.loc[(agreement_numeric >= agreement_min) | agreement_numeric.isna(), :]
    merged.attrs["gift_trait_token_map"] = trait_token_map
    return merged


def fetch_gbif_distribution_table(database, config, species, downloads_dir, timeout, dry_run):
    return _fetch_gbif_distribution_table(
        database, config, species, downloads_dir, timeout, dry_run, fetch_json=fetch_json_payload,
    )


def load_database_table(
    database: str,
    config: Dict[str, str],
    plan_rows: Sequence[TraitPlanRow],
    species: Sequence[str],
    downloads_dir: Path,
    timeout: float,
    dry_run: bool,
) -> Optional[pandas.DataFrame]:
    default_mode = SUPPORTED_DATABASES.get(database, {}).get("acquisition_mode", "bulk")
    acquisition_mode = str(config.get("acquisition_mode", default_mode) or default_mode).strip().lower()
    if acquisition_mode == "gbif_distribution" and database != "gbif":
        raise ValueError("GBIF occurrence acquisition must use database=gbif so observation roles and quality cannot be lost")
    if database == "gbif" and acquisition_mode != "gbif_distribution":
        raise ValueError("GBIF observations require acquisition_mode=gbif_distribution; use gbif_occurrence_file for local exports")
    if acquisition_mode == "bulk":
        return load_bulk_database(
            database=database,
            config=config,
            downloads_dir=downloads_dir,
            timeout=timeout,
            dry_run=dry_run,
        )
    if acquisition_mode == "species_api":
        return fetch_species_api_table(
            database=database,
            config=config,
            species=species,
            timeout=timeout,
            dry_run=dry_run,
        )
    if acquisition_mode == "gift_api":
        return fetch_gift_api_table(
            database=database,
            config=config,
            plan_rows=plan_rows,
            species=species,
            timeout=timeout,
            dry_run=dry_run,
            downloads_dir=downloads_dir,
        )
    if acquisition_mode == "gbif_distribution":
        return fetch_gbif_distribution_table(
            database=database,
            config=config,
            species=species,
            downloads_dir=downloads_dir,
            timeout=timeout,
            dry_run=dry_run,
        )
    raise ValueError("Unsupported acquisition_mode '{}' for '{}'".format(acquisition_mode, database))


def format_output_value(value: object) -> str:
    if value is pandas.NA:
        return ""
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        if value.is_integer():
            return str(int(value))
        if not math.isfinite(value):
            raise ValueError("Cannot publish a non-finite trait")
        return repr(value)
    if isinstance(value, (int,)):
        return str(value)
    text = str(value).strip()
    return text


def aggregate_trait_column(
    db_df: pandas.DataFrame,
    plan_row: TraitPlanRow,
) -> pandas.Series:
    validate_trait_plan_row(plan_row)
    group_keys = db_df["__species_norm"]
    values = db_df[plan_row.source_column]
    if plan_row.value_type == "text":
        valid = strip_string_series(values)
        mask = values.notna() & valid.ne("")
        grouped_text = valid[mask].groupby(group_keys[mask], observed=True, sort=False)
        if plan_row.aggregation == "first":
            return grouped_text.first()
        return grouped_text.agg(lambda items: json.dumps(sorted(set(items)), ensure_ascii=False))
    if plan_row.value_type == "binary":
        if plan_row.positive_values:
            normalized = strip_string_series(values, lower=True)
            missing = values.isna() | normalized.isin(["", "na", "nan", "n/a", "null", "none", "unknown"])
            mapped = normalized.isin(plan_row.positive_values).astype(float).mask(missing)
        else:
            numeric = finite_trait_values(values)
            mapped = (numeric > 0).astype(float).mask(numeric.isna())
        grouped = mapped.groupby(group_keys, observed=True, sort=False)
        if plan_row.aggregation in ("all", "min"):
            return grouped.min()
        if plan_row.aggregation == "sum":
            return grouped.sum(min_count=1)
        if plan_row.aggregation == "mean":
            return grouped.mean()
        return grouped.max()

    if plan_row.value_type == "categorical":
        normalized = strip_string_series(values)
        normalized = normalized.mask(normalized.str.lower().isin(["na", "nan", "n/a", "null", "none", "unknown"]), "")
        if plan_row.aggregation == "first":
            valid_mask = normalized != ""
            if not valid_mask.any():
                return pandas.Series(dtype=object)
            return normalized.loc[valid_mask].groupby(group_keys.loc[valid_mask], observed=True, sort=False).first()
        return aggregate_categorical_mode(values=normalized, group_keys=group_keys)

    numeric = finite_trait_values(values)
    grouped = numeric.groupby(group_keys, observed=True, sort=False)
    result = getattr(grouped, plan_row.aggregation)()
    if result.dropna().map(lambda value: not math.isfinite(float(value))).any():
        raise ValueError("Numeric trait aggregation overflowed")
    return result


def finite_trait_values(values):
    normalized = strip_string_series(values, lower=True)
    missing = values.isna() | normalized.isin(["", "na", "nan", "n/a", "null", "none", "unknown"])
    numeric = pandas.to_numeric(values.mask(missing), errors="coerce")
    invalid = ~missing & (numeric.isna() | numeric.map(lambda value: not math.isfinite(float(value))))
    if invalid.any():
        raise ValueError("Non-numeric or non-finite trait values: " + repr(values[invalid].head(3).tolist()))
    return numeric


def validate_species_source(
    species_source: str,
    manifest_path: Path,
    species_cds_dir: Path,
) -> Set[str]:
    if species_source == "download_manifest":
        if not manifest_path.exists():
            raise FileNotFoundError("Download manifest not found: {}".format(manifest_path))
        species = species_from_manifest(manifest_path)
        if len(species) == 0:
            raise ValueError("No species could be parsed from manifest: {}".format(manifest_path))
        return species
    if species_source == "species_cds":
        if not species_cds_dir.exists():
            raise FileNotFoundError("species_cds directory not found: {}".format(species_cds_dir))
        species = species_from_cds_dir(species_cds_dir)
        if len(species) == 0:
            raise ValueError("No species could be parsed from species_cds directory: {}".format(species_cds_dir))
        return species
    raise ValueError("Unsupported species-source: {}".format(species_source))


def print_supported_databases() -> None:
    print("database\tacquisition_mode\tretrieval_scope\tnotes")
    for database in sorted(SUPPORTED_DATABASES.keys()):
        info = SUPPORTED_DATABASES[database]
        print(
            "{}\t{}\t{}\t{}".format(
                database,
                info.get("acquisition_mode", ""),
                info.get("scope", ""),
                info.get("notes", ""),
            )
        )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate workspace/input/species_trait/species_trait.tsv from selected species and trait sources. "
            "Species can be resolved from download_plan manifest (default) or species_cds filenames."
        )
    )
    parser.add_argument(
        "--species-source",
        choices=("download_manifest", "species_cds"),
        default="download_manifest",
        help="Source of target species set.",
    )
    parser.add_argument(
        "--download-manifest",
        default=str(DEFAULT_MANIFEST_PATH),
        help="Path to download_plan manifest (XLSX/TSV/CSV). Used when --species-source=download_manifest.",
    )
    parser.add_argument(
        "--species-cds-dir",
        default="workspace/output/input_generation/species_cds",
        help="species_cds directory used when --species-source=species_cds.",
    )
    parser.add_argument(
        "--trait-plan",
        default=str(DEFAULT_TRAIT_PLAN_PATH),
        help="Trait extraction plan TSV path.",
    )
    parser.add_argument(
        "--database-sources",
        default=str(DEFAULT_DB_SOURCES_PATH),
        help="Trait database source map TSV path.",
    )
    parser.add_argument(
        "--databases",
        default="auto",
        help="Comma-separated database IDs, or 'auto' (from trait_plan), or 'all'.",
    )
    parser.add_argument(
        "--downloads-dir",
        default=str(DEFAULT_DOWNLOADS_DIR),
        help="Directory for cached/downloaded trait database files.",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT_PATH),
        help="Output species trait TSV path.",
    )
    parser.add_argument(
        "--download-timeout",
        type=float,
        default=120.0,
        help="Network timeout in seconds for database retrieval.",
    )
    parser.add_argument("--print-gift-mapping-inputs", action="store_true",
                        help="Print configured reviewed mapping files for local provenance; no network requests.")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail when required inputs are missing or any requested output trait has no observed values.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Resolve and report actions without writing files.",
    )
    parser.add_argument(
        "--stats-output",
        default="",
        help="Optional JSON stats output path.",
    )
    parser.add_argument(
        "--print-supported-databases",
        action="store_true",
        help="Print supported database IDs with retrieval mode and exit.",
    )
    parser.add_argument(
        "--print-gift-traits",
        action="store_true",
        help="Print GIFT traits_meta table and exit.",
    )
    parser.add_argument(
        "--gift-api",
        default=DEFAULT_GIFT_API,
        help="Base URI for GIFT API (used with --print-gift-traits).",
    )
    parser.add_argument(
        "--gift-version",
        default="latest",
        help="GIFT version (latest|beta|stable version; used with --print-gift-traits).",
    )
    parser.add_argument(
        "--gift-versions-api",
        default="",
        help="Optional base URI for GIFT versions endpoint (used with --print-gift-traits).",
    )
    parser.add_argument(
        "--gift-trait-search",
        default="",
        help="Case-insensitive keyword filter for --print-gift-traits.",
    )
    parser.add_argument(
        "--gift-trait-limit",
        type=int,
        default=0,
        help="Max rows for --print-gift-traits (0 means no limit).",
    )
    parser.add_argument(
        "--gbif-api",
        default=None,
        help="Base URI for the GBIF API used by gbif_distribution.",
    )
    parser.add_argument(
        "--gbif-page-size",
        type=int,
        default=None,
        help="Occurrence search page size for gbif_distribution.",
    )
    parser.add_argument(
        "--gbif-max-occurrences-per-species",
        type=int,
        default=None,
        help="Maximum occurrence-search records to fetch per species without GBIF login.",
    )
    parser.add_argument(
        "--gbif-grid-degrees",
        type=float,
        default=None,
        help="Grid size in degrees for the observed occupied-cell area (not IUCN AOO).",
    )
    parser.add_argument(
        "--gbif-min-match-confidence",
        type=float,
        default=None,
        help="Minimum GBIF species-match confidence used for occurrence retrieval.",
    )
    parser.add_argument(
        "--gbif-max-coordinate-uncertainty-m",
        default="",
        help="Optional maximum coordinateUncertaintyInMeters for GBIF occurrence points.",
    )
    parser.add_argument(
        "--gbif-min-distance-from-known-centroid-m",
        default="",
        help="Exclude points nearer than this distance to known georeferencing centroids (not the species mean).",
    )
    parser.add_argument("--print-gbif-input-files", action="store_true", help="Print local GBIF source files for workflow provenance and exit.")
    parser.add_argument("--print-gbif-input-identity", action="store_true", help="Print the effective GBIF config/source identity for workflow provenance and exit.")
    for name, help_text in (
        ("gbif-year-min", "Minimum event year; date intervals must fit entirely inside the window."),
        ("gbif-year-max", "Maximum event year."),
        ("gbif-countries", "Comma-separated country codes to retain."),
        ("gbif-include-basis-of-record", "Comma-separated basisOfRecord values to retain."),
        ("gbif-exclude-basis-of-record", "Comma-separated basisOfRecord values to exclude."),
        ("gbif-include-establishment-means", "Retain these establishmentMeans values; unknown is not native."),
        ("gbif-missing-date", "keep|exclude for unknown dates with an active year filter (default exclude)."),
        ("gbif-missing-uncertainty", "keep|exclude for unknown uncertainty with an active threshold (default keep)."),
        ("gbif-missing-centroid-distance", "keep|exclude for unknown known-centroid distance (default keep)."),
        ("gbif-use-cache", "yes|no; reuse a verified saved acquisition or acquire a new snapshot (default yes)."),
        ("gbif-require-complete", "yes|no; fail instead of publishing NA for incomplete/unresolved acquisitions (default no)."),
        ("gbif-occurrence-file", "Local GBIF SIMPLE_CSV TSV, TSV.gz, or single-table ZIP; no API requests."),
        ("gbif-taxon-map", "Reviewed TSV with species, taxon_key and scientific_name for local records."),
        ("gbif-download-metadata", "Saved official download metadata JSON; required to establish complete_download."),
    ):
        parser.add_argument("--" + name, default="", help=help_text)
    return parser


def apply_gbif_cli_overrides(config: Dict[str, str], args: argparse.Namespace) -> Dict[str, str]:
    merged = dict(config)
    for name, value in vars(args).items():
        if not name.startswith("gbif_") or value is None or str(value).strip() == "":
            continue
        merged["uri" if name == "gbif_api" else name] = str(value).strip()
    return merged


def validate_trait_output_paths(outputs, input_paths, protected_directories):
    outputs = [Path(path) for path in outputs if path is not None]
    for index, path in enumerate(outputs):
        if path.is_symlink() or (path.exists() and not path.is_file()):
            raise ValueError("Trait output must be a regular file: " + str(path))
        for other in list(input_paths) + outputs[:index]:
            other = Path(other)
            if path.resolve() == other.resolve() or (path.exists() and other.exists() and os.path.samefile(path, other)):
                raise ValueError("Trait output aliases an input or another output: " + str(path))
        if any(path.resolve().is_relative_to(Path(directory).resolve()) for directory in protected_directories):
            raise ValueError("Trait output is inside an input/cache directory: " + str(path))


def publish_trait_outputs(payloads):
    """Prepare every file before installation and roll back on publication errors."""
    staged = {}
    backups = {}
    installed = []
    try:
        for path, payload in payloads.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            if path.is_symlink() or (path.exists() and not path.is_file()):
                raise ValueError("Trait output must be a regular file: " + str(path))
            with tempfile.NamedTemporaryFile(dir=path.parent, prefix="." + path.name, mode="wb", delete=False) as handle:
                staged[path] = Path(handle.name)
                handle.write(payload)
                handle.flush()
                os.fsync(handle.fileno())
            if path.exists():
                with tempfile.NamedTemporaryFile(dir=path.parent, prefix="." + path.name + ".backup", delete=False) as handle:
                    backups[path] = Path(handle.name)
                shutil.copy2(path, backups[path])
                shutil.copymode(path, staged[path])
        for path, temporary in staged.items():
            os.replace(temporary, path)
            installed.append(path)
    except BaseException as original_error:
        recovery_errors = []
        for path in reversed(installed):
            try:
                if path in backups:
                    backup = backups[path]
                    os.replace(backup, path)
                    del backups[path]
                else:
                    path.unlink()
            except OSError as recovery_error:
                backup = backups.pop(path, None)
                recovery_errors.append("{}: {}; retained backup: {}".format(path, recovery_error, backup))
        if recovery_errors:
            raise RuntimeError("Trait publication rollback needs recovery: " + "; ".join(recovery_errors)) from original_error
        raise
    finally:
        for temporary in list(staged.values()) + list(backups.values()):
            temporary.unlink(missing_ok=True)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.print_gift_mapping_inputs:
        path = Path(args.database_sources).expanduser().resolve()
        if path.exists():
            for config in read_database_sources(path).values():
                if config.get("gift_species_mapping_file"):
                    print(config["gift_species_mapping_file"])
        return 0
    if args.print_gbif_input_identity or args.print_gbif_input_files:
        sources = read_database_sources(Path(args.database_sources)) if Path(args.database_sources).exists() else {}
        plan = read_trait_plan(Path(args.trait_plan)) if Path(args.trait_plan).exists() else []
        requested = resolve_requested_databases(args.databases, plan, sources)
        if "gbif" not in requested:
            if args.print_gbif_input_identity:
                print("not_requested")
            return 0
        config = apply_gbif_cli_overrides(sources.get("gbif", {}), args)
        if args.print_gbif_input_files:
            config = gbif_effective_config(config)
            for key in ("gbif_occurrence_file", "gbif_taxon_map", "gbif_download_metadata"):
                if config[key]:
                    print(config[key])
        else:
            print(gbif_input_identity(config, Path(args.output).expanduser().resolve()))
        return 0

    if args.print_supported_databases:
        print_supported_databases()
        return 0
    if args.print_gift_traits:
        if int(args.gift_trait_limit) < 0:
            parser.error("--gift-trait-limit must be >= 0")
        gift_api_base = normalize_base_uri(args.gift_api, DEFAULT_GIFT_API)
        gift_versions_api = str(args.gift_versions_api or "").strip()
        if gift_versions_api != "":
            gift_versions_api = normalize_base_uri(gift_versions_api, DEFAULT_GIFT_API)
        print_gift_traits(
            data_api_base=gift_api_base,
            requested_version=str(args.gift_version or "latest").strip().lower(),
            versions_api=gift_versions_api,
            timeout=float(args.download_timeout),
            search=str(args.gift_trait_search or "").strip(),
            limit=int(args.gift_trait_limit),
        )
        return 0

    manifest_path = Path(args.download_manifest).expanduser().resolve()
    species_cds_dir = Path(args.species_cds_dir).expanduser().resolve()
    trait_plan_path = Path(args.trait_plan).expanduser().resolve()
    db_sources_path = Path(args.database_sources).expanduser().resolve()
    downloads_dir = Path(args.downloads_dir).expanduser().resolve()
    output_path = Path(args.output).expanduser().absolute()
    schema_output_path = schema_path(output_path)
    stats_output_path = Path(args.stats_output).expanduser().absolute() if args.stats_output else None

    warnings: List[str] = []
    errors: List[str] = []

    try:
        target_species = validate_species_source(
            species_source=args.species_source,
            manifest_path=manifest_path,
            species_cds_dir=species_cds_dir,
        )
    except Exception as exc:
        parser.error(str(exc))

    _log("Target species source: {}".format(args.species_source))
    _log("Resolved target species: {}".format(len(target_species)))

    trait_plan_missing = not trait_plan_path.exists()
    if trait_plan_missing:
        message = "Trait plan not found: {}".format(trait_plan_path)
        warnings.append(message)
        plan_rows: List[TraitPlanRow] = []
    else:
        plan_rows = read_trait_plan(trait_plan_path)

    if db_sources_path.exists():
        source_rows = read_database_sources(db_sources_path)
    else:
        source_rows = {}
        warnings.append("Database source map not found: {}".format(db_sources_path))

    input_paths = [manifest_path, trait_plan_path, db_sources_path, Path(__file__),
                   SCRIPT_DIR / "gift_retrieval.py", SCRIPT_DIR / "gift_species_mappings.tsv",
                   SCRIPT_DIR / "species_trait_schema.py"]
    for config in source_rows.values():
        if config.get("gift_species_mapping_file"):
            input_paths.append(Path(config["gift_species_mapping_file"]))
        for uri in split_uri_list(config.get("uri", "")):
            parsed = urlparse(uri)
            if parsed.scheme in {"", "file"}:
                input_paths.append(Path(unquote(parsed.path)))
    gbif_config = apply_gbif_cli_overrides(source_rows.get("gbif", {}), args)
    input_paths.extend(Path(gbif_config[key]).expanduser() for key in ("gbif_occurrence_file", "gbif_taxon_map", "gbif_download_metadata") if gbif_config.get(key))
    input_paths.extend([SCRIPT_DIR / "gbif_observations.py", SCRIPT_DIR / "species_trait_contract.py"])
    protected_directories = [downloads_dir / "gift", downloads_dir / "gbif"]
    if args.species_source == "species_cds":
        protected_directories.append(species_cds_dir)
    try:
        validate_trait_output_paths([output_path, schema_output_path, *sidecar_paths(output_path).values(), stats_output_path], input_paths, protected_directories)
    except ValueError as exc:
        parser.error(str(exc))

    requested_databases = resolve_requested_databases(
        databases_arg=args.databases,
        plan_rows=plan_rows,
        source_rows=source_rows,
    )
    if len(requested_databases) == 0 and len(plan_rows) > 0:
        requested_databases = sorted({row.database for row in plan_rows})
    plan_rows = add_builtin_trait_plan_rows(plan_rows=plan_rows, requested_databases=requested_databases)
    if trait_plan_missing and args.strict and len(plan_rows) == 0:
        parser.error("Trait plan not found: {}".format(trait_plan_path))

    species_sorted = sorted(target_species)
    target_species_by_base: Dict[str, Set[str]] = {}
    for species_name in species_sorted:
        base_label = base_species_label(species_name)
        if "_" not in base_label:
            continue
        target_species_by_base.setdefault(base_label, set()).add(species_name)
    result = pandas.DataFrame({"species": species_sorted}).set_index("species")
    trait_types = {}
    for row in plan_rows:
        if row.database in requested_databases:
            validate_trait_plan_row(row)
            if row.output_trait in trait_types and trait_types[row.output_trait] != row.value_type:
                raise ValueError("Conflicting value types for output trait: " + row.output_trait)
            trait_types[row.output_trait] = row.value_type
            result[row.output_trait] = ""
    db_frames: Dict[str, pandas.DataFrame] = {}
    db_configs: Dict[str, Dict[str, str]] = {}
    trait_key_series_cache: Dict[Tuple[str, str], pandas.Series] = {}
    gbif_bundle = None
    trait_definitions = {}

    for database in requested_databases:
        config = source_rows.get(database, {})
        default_mode = SUPPORTED_DATABASES.get(database, {}).get("acquisition_mode", "bulk")
        if "acquisition_mode" not in config:
            config = dict(config)
            config["acquisition_mode"] = default_mode
        if database == "gbif":
            config = apply_gbif_cli_overrides(config=config, args=args)
        try:
            db_table = load_database_table(
                database=database,
                config=config,
                plan_rows=plan_rows,
                species=species_sorted,
                downloads_dir=downloads_dir,
                timeout=float(args.download_timeout),
                dry_run=bool(args.dry_run),
            )
        except Exception as exc:
            message = "[{}] failed to load source: {}".format(database, exc)
            if args.strict or database == "gbif":
                errors.append(message)
            else:
                warnings.append(message)
            continue

        if db_table is None:
            warnings.append("[{}] source table is unavailable.".format(database))
            continue
        if db_table.shape[0] == 0:
            warnings.append("[{}] source table is empty.".format(database))
            continue

        if database == "gbif":
            gbif_bundle = db_table.attrs.get("gbif_bundle")
        db_configs[database] = config
        species_column = detect_species_column(db_table, config.get("species_column", ""))
        db_table = db_table.copy()
        db_table["__species_norm"] = db_table[species_column].map(
            lambda value: str(value) if str(value) in target_species else normalize_species_name(value))
        db_table["__species_norm"] = db_table["__species_norm"].map(
            lambda value: remap_species_name_to_target(
                species_name=value,
                target_species=target_species,
                target_species_by_base=target_species_by_base,
            )
        )
        db_table = db_table.loc[db_table["__species_norm"] != "", :]
        if db_table.shape[0] == 0:
            warnings.append("[{}] no rows matched target species.".format(database))
            continue
        db_frames[database] = db_table
        _log("[{}] matched rows: {}".format(database, db_table.shape[0]))

    for plan_row in plan_rows:
        if plan_row.database not in requested_databases:
            continue
        if plan_row.database == "gbif":
            if plan_row.source_column not in METRIC_DEFINITIONS:
                errors.append("GBIF quality/legacy fields cannot be exported as traits: " + plan_row.source_column)
                continue
            if plan_row.value_type != "numeric" or plan_row.trait_key:
                errors.append("GBIF observation metrics must remain numeric without trait-key filtering")
                continue
            unit, meaning = METRIC_DEFINITIONS[plan_row.source_column]
            definition = {"source": "gbif", "source_column": plan_row.source_column, "role": "observation", "unit": unit, "meaning": meaning}
        else:
            definition = {"source": plan_row.database, "source_column": plan_row.source_column, "role": "trait"}
        prior_definition = trait_definitions.get(plan_row.output_trait)
        if prior_definition and prior_definition != definition and (definition["source"] == "gbif" or prior_definition["source"] == "gbif"):
            errors.append("A GBIF observation cannot share an output column with another definition: " + plan_row.output_trait)
            continue
        trait_definitions[plan_row.output_trait] = definition
        if args.dry_run:
            continue
        db_df = db_frames.get(plan_row.database)
        if db_df is None:
            message = "[{}] no data frame available for trait '{}'.".format(
                plan_row.database,
                plan_row.output_trait,
            )
            if args.strict:
                errors.append(message)
            else:
                warnings.append(message)
            continue
        db_config = db_configs.get(plan_row.database, {})
        db_filtered = db_df
        if plan_row.trait_key != "":
            trait_key_column = plan_row.trait_key_column or str(db_config.get("trait_key_column", "") or "").strip()
            if trait_key_column == "":
                trait_key_column = "trait_name"
            if trait_key_column not in db_filtered.columns:
                message = "[{}] missing trait_key_column '{}'.".format(plan_row.database, trait_key_column)
                if args.strict:
                    errors.append(message)
                else:
                    warnings.append(message)
                continue
            cache_key = (plan_row.database, trait_key_column)
            trait_key_values = trait_key_series_cache.get(cache_key)
            if trait_key_values is None:
                trait_key_values = strip_string_series(db_df[trait_key_column])
                trait_key_series_cache[cache_key] = trait_key_values
            resolved_key = db_df.attrs.get("gift_trait_token_map", {}).get(plan_row.trait_key, plan_row.trait_key)
            if plan_row.database == "gift" and trait_key_column in {"trait_ID", "trait_token"}:
                db_filtered = db_filtered.loc[strip_string_series(db_df["trait_ID"]) == resolved_key, :]
            else:
                db_filtered = db_filtered.loc[trait_key_values == plan_row.trait_key, :]
            if (
                db_filtered.shape[0] == 0
                and plan_row.database == "gift"
                and trait_key_column == "trait_ID"
                and "trait_token" in db_df.columns
            ):
                # Allow GIFT plan rows to use trait names in trait_key while source uses resolved trait_ID.
                cache_key = (plan_row.database, "trait_token")
                trait_token_values = trait_key_series_cache.get(cache_key)
                if trait_token_values is None:
                    trait_token_values = strip_string_series(db_df["trait_token"])
                    trait_key_series_cache[cache_key] = trait_token_values
                db_filtered = db_df.loc[trait_token_values == plan_row.trait_key, :]
            if db_filtered.shape[0] == 0:
                message = "[{}] no rows matched trait_key '{}' in '{}'.".format(
                    plan_row.database,
                    plan_row.trait_key,
                    trait_key_column,
                )
                if args.strict:
                    errors.append(message)
                else:
                    warnings.append(message)
                continue

        if plan_row.source_column not in db_filtered.columns:
            message = "[{}] missing source column '{}'.".format(plan_row.database, plan_row.source_column)
            if args.strict:
                errors.append(message)
            else:
                warnings.append(message)
            continue

        try:
            aggregated = aggregate_trait_column(db_df=db_filtered, plan_row=plan_row)
        except ValueError as exc:
            message = "[{}] invalid trait '{}': {}".format(plan_row.database, plan_row.output_trait, exc)
            (errors if args.strict else warnings).append(message)
            continue
        colname = plan_row.output_trait
        if colname not in result.columns:
            result[colname] = ""
        formatted = aggregated.map(format_output_value)
        formatted = formatted.loc[formatted != ""]
        if formatted.shape[0] > 0:
            empty_mask = result[colname] == ""
            if empty_mask.any():
                result.loc[empty_mask, colname] = formatted.reindex(result.index[empty_mask], fill_value="")

    if len(errors) > 0:
        for message in errors:
            _log("ERROR: {}".format(message))
        return 1

    output_df = result.reset_index()
    trait_columns = [col for col in output_df.columns if col != "species"]
    num_species_with_any_trait = 0
    if len(trait_columns) > 0:
        has_any = output_df.loc[:, trait_columns].ne("").any(axis=1)
        num_species_with_any_trait = int(has_any.sum())

    if args.strict and len(trait_columns) == 0:
        _log("ERROR: No trait columns were generated.")
        return 1

    empty_traits = [col for col in trait_columns if not output_df[col].ne("").any()]
    if args.strict and empty_traits and not args.dry_run:
        _log("ERROR: No observed values for requested traits: {}".format(", ".join(empty_traits)))
        return 1

    payloads = {}
    if args.dry_run:
        _log("[dry-run] trait schema would be written to: {}".format(schema_output_path))
        _log("[dry-run] species_trait output would be written to: {}".format(output_path))

    if stats_output_path is not None:
        stats = {
            "species_source": args.species_source,
            "num_target_species": len(species_sorted),
            "num_trait_columns": len(trait_columns),
            "num_species_with_any_trait": num_species_with_any_trait,
            "num_observed_by_trait": {col: int(output_df[col].ne("").sum()) for col in trait_columns},
            "num_requested_databases": len(requested_databases),
            "num_loaded_databases": len(db_frames),
            "output_path": str(output_path),
            "dry_run": int(bool(args.dry_run)),
        }
        payloads[stats_output_path] = json.dumps(stats, indent=2, ensure_ascii=False, allow_nan=False).encode("utf-8")
        if args.dry_run:
            _log("[dry-run] stats output would be written to: {}".format(stats_output_path))

    if not args.dry_run:
        validate_trait_output_paths([output_path, *sidecar_paths(output_path).values(), stats_output_path], input_paths, protected_directories)
        if gbif_bundle is not None:
            input_paths.append(Path(gbif_bundle["records_path"]))
        bundle_payloads = trait_bundle_payloads(output_df, output_path, {key: value for key, value in trait_definitions.items() if key in trait_columns}, gbif_bundle)
        # Publish stats and all sidecars together, with metadata last.
        payloads = {**payloads, schema_output_path: schema_payload(bundle_payloads[output_path], trait_types), **bundle_payloads}
        validate_trait_output_paths(list(payloads), input_paths, protected_directories)
        publish_trait_outputs(payloads)
        _log("species_trait.tsv written: {}".format(output_path))

    for message in warnings:
        _log("WARNING: {}".format(message))
    _log("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
