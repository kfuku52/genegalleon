"""Validated, resumable public GIFT JSON requests and reviewed name mappings."""

import csv
import hashlib
import json
import math
import os
import re
import tempfile
import time
import uuid
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from http.client import IncompleteRead, RemoteDisconnected
from io import StringIO
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import parse_qs, urlencode, urlparse, urlunparse

CACHE_SCHEMA = "genegalleon-gift-page-v2"
MAPPING_FILE = Path(__file__).with_name("gift_species_mappings.tsv")


def canonical_json(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix="." + path.name, mode="wb", delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(canonical_json(value))
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def validate_gift_payload(url, payload):
    """An API error object is never an empty result or a reusable page."""
    rows = payload
    if isinstance(payload, dict):
        rows = next((payload[key] for key in ("data", "results", "items") if isinstance(payload.get(key), list)), None)
    if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
        raise ValueError("GIFT returned a non-table JSON response: " + url)
    query = parse_qs(urlparse(url).query).get("query", [""])[0]
    required = {"traits": "work_ID", "traits_meta": "Lvl3", "versions": "version"}.get(query)
    if required and any(str(row.get(required, "") or "").strip() == "" for row in rows):
        raise ValueError("GIFT response lacks {}: {}".format(required, url))
    if query.startswith("names_matched") and any("work_ID" not in row for row in rows):
        raise ValueError("GIFT species response lacks work_ID: " + url)
    if query == "traits":
        for row in rows:
            if not re.fullmatch(r"[1-9][0-9]*", str(row["work_ID"])) or "trait_value" not in row:
                raise ValueError("GIFT trait row requires a positive work_ID and trait_value: " + url)
            if isinstance(row["trait_value"], (dict, list)):
                raise ValueError("GIFT trait_value must be scalar: " + url)
            agreement = row.get("agreement")
            if str(agreement or "").strip().lower() not in {"", "na", "nan", "n/a", "null", "none", "unknown"}:
                try:
                    valid_agreement = math.isfinite(float(agreement)) and 0 <= float(agreement) <= 1
                except (TypeError, ValueError):
                    valid_agreement = False
                if not valid_agreement:
                    raise ValueError("Invalid GIFT agreement value: " + url)
    canonical_json(payload)  # reject non-finite JSON numbers before writing
    return payload


class GiftRetrieval:
    """Cache keys include the exact endpoint, release, bias flags and page bounds.

    Stable releases are reused until refresh is requested. Latest-version
    discovery and beta pages are always fetched online; offline requires a
    concrete release. A failed run retains only individually validated pages.
    """

    def __init__(self, directory, fetcher, mode="reuse", retries=2, logger=print):
        if mode not in {"reuse", "refresh", "offline"}:
            raise ValueError("gift_cache_mode must be reuse, refresh or offline")
        if not 0 <= retries <= 5:
            raise ValueError("gift_retries must be between 0 and 5")
        self.directory = Path(directory) if directory is not None else None
        if mode == "offline" and self.directory is None:
            raise ValueError("offline GIFT retrieval requires a downloads directory")
        self.fetcher, self.mode, self.retries, self.logger = fetcher, mode, retries, logger
        self.report = {
            "schema": "genegalleon-gift-retrieval-v1",
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "cache_mode": mode,
            "requests": [],
            "species": [],
            "traits": [],
        }
        self.report_path = self.directory / "runs" / (uuid.uuid4().hex + ".json") if self.directory else None
        self.seen = set()
        self.generations = {}

    def __enter__(self):
        self.report["status"] = "running"
        if self.report_path:
            atomic_json(self.report_path, self.report)
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.report.update(status="failed" if exc else "complete", finished_utc=datetime.now(timezone.utc).isoformat())
        if exc:
            self.report["error"] = str(exc)
        if self.report_path:
            atomic_json(self.report_path, self.report)
            self.logger("[gift] retrieval report: " + str(self.report_path))
        return False

    def cache_path(self, url, mutable):
        if self.directory is None:
            return None
        parsed = urlparse(url)
        params = parse_qs(parsed.query, keep_blank_values=True)
        if params.get("query") != ["traits"]:
            return self.directory / "pages" / (hashlib.sha256(url.encode()).hexdigest() + ".json")
        # Pin a generation for the whole paginated query. A refresh starts a new
        # namespace before fetching anything; interrupted refreshes resume there,
        # and concurrent readers continue using their original generation.
        params.pop("startat", None)
        scope = urlunparse(parsed._replace(query=urlencode(sorted(params.items()), doseq=True), fragment=""))
        scope_key = hashlib.sha256(scope.encode()).hexdigest()
        if scope_key not in self.generations:
            pointer = self.directory / "scopes" / (scope_key + ".json")
            generation = None
            # gg-cache-guard: audited - query scope/schema and the pinned generation UUID are validated.
            if self.mode != "refresh" and not mutable and pointer.exists():
                try:
                    state = json.loads(pointer.read_text(encoding="utf-8"))
                    if (
                        not isinstance(state, dict)
                        or state.get("schema") != CACHE_SCHEMA
                        or state.get("scope") != scope
                        or not re.fullmatch(r"[0-9a-f]{32}", str(state.get("generation", "")))
                    ):
                        raise ValueError("invalid cache generation")
                    generation = state["generation"]
                except (OSError, ValueError, TypeError) as exc:
                    if self.mode == "offline":
                        raise ValueError("Invalid offline GIFT cache generation: " + scope) from exc
                    self.logger("WARNING: [gift] invalid cache generation; starting again: " + scope)
            if generation is None:
                if self.mode == "offline":
                    raise FileNotFoundError("Missing offline GIFT cache generation: " + scope)
                generation = uuid.uuid4().hex
                atomic_json(pointer, {"schema": CACHE_SCHEMA, "scope": scope, "generation": generation})
            self.generations[scope_key] = generation
            self.report.setdefault("generations", {})[scope] = generation
        return (
            self.directory
            / "pages"
            / scope_key
            / self.generations[scope_key]
            / (hashlib.sha256(url.encode()).hexdigest() + ".json")
        )

    def fetch(self, url, timeout):
        try:
            return self._fetch(url, timeout)
        except Exception as exc:
            self.report["requests"].append({"url": url, "source": "failed", "error": str(exc)})
            raise

    def _fetch(self, url, timeout):
        query = parse_qs(urlparse(url).query).get("query", [""])[0]
        mutable = query == "versions" or urlparse(url).path.endswith("/index.php")
        if mutable and self.mode == "offline":
            raise ValueError("Offline GIFT retrieval requires an explicit stable gift_version, not latest or beta")
        path = self.cache_path(url, mutable)
        reuse = not mutable and (self.mode != "refresh" or url in self.seen)
        # gg-cache-guard: audited - exact URL, schema and payload SHA-256 are checked before reuse.
        if path and reuse and path.exists():
            try:
                record = json.loads(path.read_text(encoding="utf-8"))
                if not isinstance(record, dict) or record.get("schema") != CACHE_SCHEMA or record.get("url") != url:
                    raise ValueError("cache identity mismatch")
                timestamp = datetime.fromisoformat(record["retrieved_utc"])
                if timestamp.tzinfo is None:
                    raise ValueError("cache timestamp lacks timezone")
                payload = record["payload"]
                digest = hashlib.sha256(canonical_json(payload)).hexdigest()
                if record["sha256"] != digest:
                    raise ValueError("cache checksum mismatch")
                validate_gift_payload(url, payload)
                self.report["requests"].append(
                    {"url": url, "source": "cache", "sha256": digest, "retrieved_utc": record["retrieved_utc"]}
                )
                self.seen.add(url)
                return payload
            except (OSError, ValueError, KeyError, TypeError) as exc:
                if self.mode == "offline":
                    raise ValueError("Invalid offline GIFT cache {}: {}".format(path, exc)) from exc
                self.logger("WARNING: [gift] invalid cache; fetching again: {} ({})".format(path, exc))
        if self.mode == "offline":
            raise FileNotFoundError("Missing offline GIFT cache page: " + url)
        for attempt in range(self.retries + 1):
            try:
                payload = validate_gift_payload(url, self.fetcher(url=url, timeout=timeout))
                break
            except (HTTPError, URLError, TimeoutError, ConnectionError, IncompleteRead, RemoteDisconnected) as exc:
                transient = not isinstance(exc, HTTPError) or exc.code in {408, 429, 500, 502, 503, 504}
                if not transient or attempt == self.retries:
                    raise
                delay = 2**attempt
                if isinstance(exc, HTTPError):
                    retry_after = exc.headers.get("Retry-After", "") if exc.headers else ""
                    if retry_after:
                        try:
                            if retry_after.isdigit():
                                server_delay = int(retry_after)
                            else:
                                when = parsedate_to_datetime(retry_after)
                                if when.tzinfo is None:
                                    when = when.replace(tzinfo=timezone.utc)
                                server_delay = max(0, math.ceil((when - datetime.now(timezone.utc)).total_seconds()))
                        except (TypeError, ValueError, OverflowError):
                            raise exc from None
                        if server_delay > 30:
                            raise
                        delay = max(delay, server_delay)
                self.logger(
                    "WARNING: [gift] transient request failure; retry {}/{} in {}s: {}".format(
                        attempt + 1, self.retries, delay, url
                    )
                )
                time.sleep(delay)
        digest = hashlib.sha256(canonical_json(payload)).hexdigest()
        record = {
            "schema": CACHE_SCHEMA,
            "url": url,
            "sha256": digest,
            "payload": payload,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        }
        if path:
            atomic_json(path, record)
        self.report["requests"].append(
            {key: record[key] for key in ("url", "sha256", "retrieved_utc")} | {"source": "network"}
        )
        self.seen.add(url)
        return payload


def load_reviewed_mappings(version, custom_path="", include_bundled=True):
    mappings = {}
    receipts = []
    paths = ([MAPPING_FILE] if include_bundled else []) + ([Path(custom_path)] if custom_path else [])
    for path in paths:
        raw = path.read_bytes()
        receipts.append({"path": str(path), "sha256": hashlib.sha256(raw).hexdigest()})
        with StringIO(raw.decode("utf-8-sig"), newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t", strict=True)
            required = {"species", "gift_version", "work_ID", "work_species", "decision", "evidence_url"}
            columns = reader.fieldnames or []
            if len(columns) != len(set(columns)) or any(not column for column in columns):
                raise ValueError("Duplicate/empty GIFT species mapping columns: " + str(path))
            if not required.issubset(columns):
                raise ValueError("GIFT species mapping is missing columns: " + str(path))
            seen = set()
            for raw_row in reader:
                if None in raw_row or any(value is None for value in raw_row.values()):
                    raise ValueError("Ragged GIFT species mapping: " + str(path))
                row = {key: value.strip() for key, value in raw_row.items()}
                key = (row["gift_version"], row["species"])
                if (
                    key in seen
                    or not all(key)
                    or row["decision"] not in {"map", "exclude"}
                    or urlparse(row["evidence_url"]).scheme != "https"
                    or not urlparse(row["evidence_url"]).hostname
                    or not re.fullmatch(r"[A-Z][A-Za-z.-]*_[a-z][A-Za-z.-]*", row["species"])
                ):
                    raise ValueError("Invalid or duplicate reviewed GIFT species mapping: " + str(key))
                seen.add(key)
                if row["decision"] == "map" and (not row["work_ID"].isdigit() or len(row["work_species"].split()) != 2):
                    raise ValueError("Reviewed GIFT mappings require a numeric work_ID and binomial work_species")
                if row["gift_version"] == version:
                    mappings[row["species"]] = row
    return mappings, receipts
