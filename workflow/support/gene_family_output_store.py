#!/usr/bin/env python3
"""ZIP-backed logical storage for GeneGalleon gene-family outputs.

Live files always take precedence over archived members.  Archives are
immutable ZIP shards grouped by logical output subdirectory.  Deletions are
recorded as tombstones so removing a live override does not unexpectedly
reveal an older archived version.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import contextlib
import csv
import hashlib
import json
import os
import re
import shutil
import stat
import sys
import threading
import time
import uuid
import weakref
import zipfile
from collections import OrderedDict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import BinaryIO, Callable, Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from content_digest_cache import cached_sha256_file
from shared_namespace_lock import NamespaceLockError, namespace_lock

STORE_DIR_NAME = ".gg_store"
ACTIVE_ARCHIVE_DIR_NAME = "archives"
LEGACY_ARCHIVE_DIR_NAME = ".gg_archives"
# ARCHIVE_DIR_NAME remains the metadata/locking root name for callers that
# import it.  ZIP payload lives separately in ACTIVE_ARCHIVE_DIR_NAME.
ARCHIVE_DIR_NAME = STORE_DIR_NAME
MANIFEST_MEMBER = ".genegalleon-manifest.json"
TOMBSTONE_FILE = "tombstones.jsonl"
FAMILY_STATE_FILE = "family-state.jsonl"
FAMILY_STATE_DIR_NAME = "family-state"
FAMILY_STATE_LOCK_DIR_NAME = "family-state-locks"
FAMILY_LOCK_DIR_NAME = "family-locks"
INDEX_FILE = "index.json"
INDEX_DIR_NAME = "index"
SUBDIR_INDEX_DIR_NAME = "index-by-subdir"
INDEX_CATALOG_FILE = "index-catalog.json"
INDEX_EPOCH_FILE = "index.epoch"
INDEX_UPDATE_FILE = "index-update.pending"
STORAGE_CONVERSION_FILE = "storage-conversion.pending"
STORAGE_CONVERSION_LOCK_FILE = "storage-conversion.lock"
STORE_METADATA_FILE = "store.json"
PURE_RAW_RETIRED_DIR_NAME = ".gg_store.pure-raw-retired"
GENERATION_FILE = "generation.counter"
LOCK_FILE = "archive.lock"
PRODUCER_LOCK_FILE = "producer.lock"
ARCHIVE_SCHEMA_VERSION = 1
INDEX_SCHEMA_VERSION = 1
STORE_METADATA_SCHEMA_VERSION = 1
FAMILY_LOCK_STRIPES = 16
MAX_REFERENCED_SHARDS_PER_SUBDIR = 8
MAX_OPEN_COMPACTION_SOURCES = 32
TOMBSTONE_LOG_COMPACT_BYTES = 4 * 1024 * 1024
DEFAULT_LARGE_ZIP_WARNING_BYTES = 20 * 1024 * 1024 * 1024
MATERIALIZATION_RECEIPT_NAME = ".gg_materialized.jsonl"
ORTHOGROUP_ID_RE = re.compile(r"^(OG\d+|HOG\d+|SP\d+)(?=$|[._-])")
EXCLUDED_SUBDIRS = {"tmp", ACTIVE_ARCHIVE_DIR_NAME}
PRECOMPRESSED_SUFFIXES = (
    ".gz",
    ".zip",
    ".bz2",
    ".xz",
    ".zst",
    ".pdf",
    ".png",
    ".jpg",
    ".jpeg",
    ".gif",
)
ArchivedSourceSignature = Tuple[int, int, int, int, int, str]

# GeneGalleon releases predating the underscore-based output layout used dots
# in both subdirectory and filename components.  Keep these mappings at the
# logical reader boundary so an archived historical workspace can be resumed
# without rewriting every ZIP up front.  Newly materialized/produced files use
# the current names and therefore migrate incrementally during normal reruns.
LEGACY_OUTPUT_PATH_RULES: Tuple[Tuple[str, str, str, str], ...] = (
    ("amas.cleaned", ".amas.cleaned.tsv", "amas_cleaned", "_amas.cleaned.tsv"),
    ("amas.original", ".amas.original.tsv", "amas_original", "_amas.original.tsv"),
    ("cds.fasta", ".cds.fasta", "cds_fasta", "_cds.fasta"),
    ("character.gff", ".gff.tsv", "character_gff_info", "_gff.tsv"),
    ("clipkit.log", ".cds.clipkit.log", "clipkit_log", "_cds.clipkit.log"),
    ("clipkit", ".cds.clipkit.fasta", "clipkit", "_cds.clipkit.fasta"),
    ("generax.nwk", ".generax.nwk", "generax_nwk", "_generax.nwk"),
    ("generax.tree", ".generax.nhx", "generax_tree", "_generax.nhx"),
    ("generax.xml", ".generax.xml", "generax_xml", "_generax.xml"),
    ("iqtree.model", ".model.gz", "iqtree_model", "_model.gz"),
    ("iqtree.tree", ".iqtree.nwk", "iqtree_tree", "_iqtree.nwk"),
    ("mafft", ".cds.aln.fasta", "mafft", "_cds.aln.fasta"),
    ("mapdNdS.dN.tree", ".mapdNdS.dN.nwk", "mapdnds_dn_tree", "_mapdNdS.dN.nwk"),
    ("mapdNdS.dS.tree", ".mapdNdS.dS.nwk", "mapdnds_ds_tree", "_mapdNdS.dS.nwk"),
    ("rpsblast", ".rpsblast.tsv", "rpsblast", "_rpsblast.tsv"),
    ("stat.branch", ".stat.branch.tsv", "stat_branch", "_stat.branch.tsv"),
    ("stat.tree", ".stat.tree.tsv", "stat_tree", "_stat.tree.tsv"),
    ("tree_plot", ".tree_plot.pdf", "tree_plot", "_tree_plot.pdf"),
)
LEGACY_SUBDIR_ALIASES = {
    legacy_subdir: current_subdir for legacy_subdir, _, current_subdir, _ in LEGACY_OUTPUT_PATH_RULES
}
SHARED_OUTPUT_SUBDIRS = frozenset({"parameters"})


def _canonical_output_path(subdir: str, name: str) -> Tuple[str, str]:
    for legacy_subdir, legacy_suffix, current_subdir, current_suffix in LEGACY_OUTPUT_PATH_RULES:
        if subdir != legacy_subdir or not name.endswith(legacy_suffix):
            continue
        family_id = name[: -len(legacy_suffix)]
        if family_id:
            return current_subdir, family_id + current_suffix
    return subdir, name


def _legacy_output_candidates(subdir: str, name: str) -> List[Tuple[str, str]]:
    candidates: List[Tuple[str, str]] = []
    for legacy_subdir, legacy_suffix, current_subdir, current_suffix in LEGACY_OUTPUT_PATH_RULES:
        if subdir != current_subdir or not name.endswith(current_suffix):
            continue
        family_id = name[: -len(current_suffix)]
        if family_id:
            candidates.append((legacy_subdir, family_id + legacy_suffix))
    return candidates


def _equivalent_output_logical_paths(subdir: str, name: str) -> Set[str]:
    logical_paths = {_safe_logical_path(subdir, name)}
    canonical_subdir, canonical_name = _canonical_output_path(subdir, name)
    logical_paths.add(_safe_logical_path(canonical_subdir, canonical_name))
    return logical_paths


class ArchiveStoreError(RuntimeError):
    """Raised when the logical archive store is inconsistent or corrupt."""


@dataclass(frozen=True)
class Artifact:
    logical_path: str
    subdir: str
    name: str
    generation: int
    live_path: Optional[Path] = None
    zip_path: Optional[Path] = None
    member_name: Optional[str] = None
    size: Optional[int] = None
    crc: Optional[int] = None
    sha256: Optional[str] = None
    mtime_ns: Optional[int] = None
    family_id: Optional[str] = None

    @property
    def is_live(self) -> bool:
        return self.live_path is not None


def _safe_logical_path(subdir: str, name: str) -> str:
    if not subdir or subdir.startswith(".") or "/" in subdir or "\\" in subdir:
        raise ArchiveStoreError(f"Unsafe logical subdirectory: {subdir!r}")
    if not name or name.startswith(".") or "/" in name or "\\" in name:
        raise ArchiveStoreError(f"Unsafe logical filename: {name!r}")
    return f"{subdir}/{name}"


def _compression_for(path: Path, compression: str = "adaptive") -> int:
    if compression == "store":
        return zipfile.ZIP_STORED
    if compression == "deflate":
        return zipfile.ZIP_DEFLATED
    if compression != "adaptive":
        raise ValueError(f"Unsupported ZIP compression mode: {compression}")
    lower_name = path.name.lower()
    if lower_name.endswith(PRECOMPRESSED_SUFFIXES):
        return zipfile.ZIP_STORED
    return zipfile.ZIP_DEFLATED


def _stat_signature(path: Path) -> Tuple[int, int, int, int, int]:
    stat_result = path.stat()
    return (
        int(stat_result.st_dev),
        int(stat_result.st_ino),
        int(stat_result.st_size),
        int(stat_result.st_mtime_ns),
        int(stat_result.st_ctime_ns),
    )


def _signature_matches(path: Path, signature: Sequence[int]) -> bool:
    try:
        return _stat_signature(path) == tuple(int(value) for value in signature)
    except FileNotFoundError:
        return False


def _fsync_directory(path: Path) -> None:
    try:
        descriptor = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    except OSError:
        pass
    finally:
        os.close(descriptor)


def _validate_archive_root(archive_root: Path) -> None:
    if archive_root.is_symlink():
        raise ArchiveStoreError(f"Symlinked GeneGalleon archive roots are not supported: {archive_root}")


def _archive_state_root(root: Path | str) -> Path:
    """Return the metadata root, retaining read support for legacy stores."""

    resolved = Path(root).resolve()
    current = resolved / STORE_DIR_NAME
    legacy = resolved / LEGACY_ARCHIVE_DIR_NAME
    if current.exists() or current.is_symlink() or not (legacy.exists() or legacy.is_symlink()):
        return current
    return legacy


def _archive_payload_root(root: Path | str) -> Path:
    """Return the immutable-shard root for the selected store layout."""

    resolved = Path(root).resolve()
    state_root = _archive_state_root(resolved)
    if state_root.name == LEGACY_ARCHIVE_DIR_NAME:
        return state_root
    return resolved / ACTIVE_ARCHIVE_DIR_NAME


def _zip_has_genegalleon_manifest(path: Path) -> bool:
    if not path.is_file() or path.is_symlink():
        return False
    try:
        with zipfile.ZipFile(path, "r") as archive:
            return MANIFEST_MEMBER in archive.namelist()
    except (OSError, zipfile.BadZipFile):
        return False


def _physical_archive_paths(root: Path | str) -> List[Path]:
    """List active shards and finalized, user-facing ZIP files."""

    resolved = Path(root).resolve()
    payload_root = _archive_payload_root(resolved)
    paths = list(payload_root.glob("*/*.zip")) if payload_root.is_dir() else []
    if _archive_state_root(resolved).name != LEGACY_ARCHIVE_DIR_NAME:
        paths.extend(path for path in resolved.glob("*.zip") if _zip_has_genegalleon_manifest(path))
    return sorted(set(paths))


def _unrelated_final_zip_candidates(root: Path | str) -> List[Path]:
    """List user ZIPs whose names could collide with managed final ZIPs."""

    resolved = Path(root).resolve()
    return sorted(
        path
        for path in resolved.glob("*.zip")
        if path.is_file()
        and not path.is_symlink()
        and not _zip_has_genegalleon_manifest(path)
        and (resolved / path.stem).is_dir()
        and not (resolved / path.stem).is_symlink()
    )


def _conflicting_final_zip_paths(
    root: Path | str,
    managed_subdirs: Set[str],
) -> List[Path]:
    return [path for path in _unrelated_final_zip_candidates(root) if path.stem in managed_subdirs]


def _catalog_conflicting_final_zip_paths(
    root: Path | str,
    valid_family_ids: Set[str],
    family_from_name: Callable[[str], Optional[str]],
) -> List[Path]:
    """Find collisions without rescanning unrelated output subdirectories."""

    conflicts: List[Path] = []
    for zip_path in _unrelated_final_zip_candidates(root):
        live_dir = zip_path.with_suffix("")
        with os.scandir(live_dir) as entries:
            for entry in entries:
                if entry.name.startswith(".") or entry.is_symlink() or not entry.is_file(follow_symlinks=False):
                    continue
                family_id = family_from_name(entry.name)
                if family_id is not None and family_id in valid_family_ids:
                    conflicts.append(zip_path)
                    break
    return conflicts


def _final_archive_path(root: Path | str, subdir: str) -> Path:
    _safe_logical_path(subdir, "__archive__")
    return Path(root).resolve() / f"{subdir}.zip"


def _family_index_bucket(family_id: Optional[str]) -> str:
    key = family_id if family_id not in {None, ""} else "__legacy__"
    return hashlib.sha256(str(key).encode("utf-8")).hexdigest()[:2]


def _family_lock_bucket(family_id: Optional[str]) -> str:
    key = family_id if family_id not in {None, ""} else "__legacy__"
    bucket = int(hashlib.sha256(str(key).encode("utf-8")).hexdigest()[:8], 16)
    return f"{bucket % FAMILY_LOCK_STRIPES:02x}"


def _subdir_index_name(subdir: str) -> str:
    _safe_logical_path(subdir, "__index__")
    digest = hashlib.sha256(subdir.encode("utf-8")).hexdigest()[:16]
    return f"{digest}.json"


def _store_lock_path(path: Path) -> Path:
    # Store migration and pure-raw conversion move/remove the metadata tree
    # while locks are held. Coordination must outlive both representations.
    for parent in path.parents:
        if parent.name in {STORE_DIR_NAME, LEGACY_ARCHIVE_DIR_NAME}:
            return parent.parent / ".gg_store_locks" / path.relative_to(parent)
    return path


def family_lock_path(archive_root: Path, family_id: str) -> Path:
    return _store_lock_path(archive_root / FAMILY_LOCK_DIR_NAME / f"{_family_lock_bucket(family_id)}.lock")


@contextlib.contextmanager
def _bucket_lock(
    lock_path: Path,
    *,
    exclusive: bool,
    nonblocking: bool = False,
) -> Iterator[bool]:
    try:
        with namespace_lock(_store_lock_path(lock_path), exclusive=exclusive, nonblocking=nonblocking) as acquired:
            yield acquired
    except NamespaceLockError as exc:
        raise ArchiveStoreError(str(exc)) from exc


def family_bucket_lock(
    archive_root: Path,
    family_id: str,
    *,
    exclusive: bool,
    nonblocking: bool = False,
) -> contextlib.AbstractContextManager[bool]:
    _validate_archive_root(archive_root)
    return _bucket_lock(
        family_lock_path(archive_root, family_id),
        exclusive=exclusive,
        nonblocking=nonblocking,
    )


def state_bucket_lock(
    archive_root: Path,
    family_id: str,
) -> contextlib.AbstractContextManager[bool]:
    _validate_archive_root(archive_root)
    path = archive_root / FAMILY_STATE_LOCK_DIR_NAME / f"{_family_lock_bucket(family_id)}.lock"
    return _bucket_lock(path, exclusive=True)


@contextlib.contextmanager
def lock_available_family_ids(
    archive_root: Path,
    family_ids: Iterable[str],
    *,
    nonblocking: bool,
) -> Iterator[Set[str]]:
    family_ids_set = {str(family_id) for family_id in family_ids}
    representative_by_bucket: Dict[str, str] = {}
    for family_id in sorted(family_ids_set):
        representative_by_bucket.setdefault(_family_lock_bucket(family_id), family_id)
    acquired_buckets: Set[str] = set()
    with contextlib.ExitStack() as stack:
        for bucket, representative in sorted(representative_by_bucket.items()):
            acquired = stack.enter_context(
                family_bucket_lock(
                    archive_root,
                    representative,
                    exclusive=True,
                    nonblocking=nonblocking,
                )
            )
            if acquired:
                acquired_buckets.add(bucket)
        yield {family_id for family_id in family_ids_set if _family_lock_bucket(family_id) in acquired_buckets}


@contextlib.contextmanager
def all_family_bucket_locks(archive_root: Path) -> Iterator[None]:
    _validate_archive_root(archive_root)
    with contextlib.ExitStack() as stack:
        for bucket_number in range(FAMILY_LOCK_STRIPES):
            bucket = f"{bucket_number:02x}"
            stack.enter_context(
                _bucket_lock(
                    archive_root / FAMILY_LOCK_DIR_NAME / f"{bucket}.lock",
                    exclusive=True,
                )
            )
        yield


def orthogroup_id_from_name(file_name: str) -> Optional[str]:
    match = ORTHOGROUP_ID_RE.match(os.path.basename(file_name))
    return None if match is None else match.group(1)


def query_id_matchers(query_ids: Iterable[str]) -> List[str]:
    return sorted((str(query_id) for query_id in query_ids), key=lambda value: (-len(value), value))


def query_id_from_name(file_name: str, matchers: Sequence[str]) -> Optional[str]:
    basename = os.path.basename(file_name)
    for query_id in matchers:
        if basename == query_id:
            return query_id
        if basename.startswith(query_id + "_") or basename.startswith(query_id + "."):
            return query_id
    return None


def query_ids_from_input_dir(path: Path) -> List[str]:
    if not path.is_dir():
        raise FileNotFoundError(f"Input query_gene directory was not found: {path}")
    return sorted(entry.name for entry in path.iterdir() if entry.is_file() and not entry.name.startswith("."))


def orthogroup_ids_from_genecount(path: Path) -> List[str]:
    if not path.is_file():
        raise FileNotFoundError(f"Orthogroup gene-count table was not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        return [row[0] for row in reader if row and row[0]]


def family_context(
    mode: str,
    query_dir: Optional[Path] = None,
    genecount: Optional[Path] = None,
) -> Tuple[List[str], Callable[[str], Optional[str]]]:
    if mode == "query2family":
        if query_dir is None:
            raise ValueError("--query-dir is required for query2family mode")
        family_ids = query_ids_from_input_dir(query_dir)
        matchers = query_id_matchers(family_ids)
        return family_ids, lambda name: query_id_from_name(name, matchers)
    if mode == "orthogroup":
        if genecount is None:
            raise ValueError("--genecount is required for orthogroup mode")
        family_ids = orthogroup_ids_from_genecount(genecount)
        return family_ids, orthogroup_id_from_name
    raise ValueError(f"Unsupported gene-family mode: {mode}")


def family_context_with_supplement(
    mode: str,
    query_dir: Optional[Path] = None,
    genecount: Optional[Path] = None,
    family_id_file: Optional[Path] = None,
) -> Tuple[List[str], Callable[[str], Optional[str]]]:
    if mode == "query2family":
        family_ids = query_ids_from_input_dir(query_dir) if query_dir is not None else []
    elif mode == "orthogroup":
        family_ids = orthogroup_ids_from_genecount(genecount) if genecount is not None else []
    else:
        raise ValueError(f"Unsupported gene-family mode: {mode}")
    if family_id_file is not None:
        if not family_id_file.is_file():
            raise FileNotFoundError(f"Family ID file was not found: {family_id_file}")
        with family_id_file.open("r", encoding="utf-8") as handle:
            family_ids.extend(family_id for line in handle if (family_id := line.rstrip("\r\n")))
    if not family_ids:
        required = "--query-dir or --family-id-file" if mode == "query2family" else "--genecount or --family-id-file"
        raise ValueError(f"{required} is required to identify family-owned files")
    family_ids = sorted(set(family_ids))
    if mode == "query2family":
        matchers = query_id_matchers(family_ids)
        return family_ids, lambda name: query_id_from_name(name, matchers)
    return family_ids, orthogroup_id_from_name


def _family_catalog_digest(family_ids: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for family_id in sorted({str(value) for value in family_ids}):
        digest.update(family_id.encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def _read_store_metadata(root: Path) -> Optional[dict]:
    metadata_path = _archive_state_root(root) / STORE_METADATA_FILE
    if metadata_path.is_symlink():
        raise ArchiveStoreError(f"Symlinked archive store metadata is not supported: {metadata_path}")
    if not metadata_path.is_file():
        return None
    try:
        with metadata_path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        if not isinstance(payload, dict):
            raise ArchiveStoreError(f"Archive store metadata is not a JSON object: {metadata_path}")
        if int(payload.get("schema_version", -1)) != STORE_METADATA_SCHEMA_VERSION:
            raise ArchiveStoreError(
                f"Unsupported archive store metadata schema in {metadata_path}: {payload.get('schema_version')!r}"
            )
        if payload.get("mode") not in {"query2family", "orthogroup"}:
            raise ArchiveStoreError(f"Archive store metadata has an invalid mode: {metadata_path}")
        return payload
    except ArchiveStoreError:
        raise
    except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
        raise ArchiveStoreError(f"Failed to read archive store metadata {metadata_path}: {exc}") from exc


def _write_store_metadata(
    root: Path,
    mode: str,
    family_ids: Iterable[str],
    *,
    compression: str = "adaptive",
    compression_level: int = 6,
    catalog_sources: Optional[Sequence[Path]] = None,
    preserve_existing_catalog: bool = False,
) -> dict:
    if mode not in {"query2family", "orthogroup"}:
        raise ValueError(f"Unsupported gene-family mode: {mode}")
    archive_root = _archive_state_root(root)
    _validate_archive_root(archive_root)
    archive_root.mkdir(parents=True, exist_ok=True)
    existing = _read_store_metadata(root)
    if existing is not None and existing.get("mode") != mode:
        raise ArchiveStoreError(f"Archive store metadata uses a different gene-family mode: {existing.get('mode')}")
    family_ids_set = sorted({str(value) for value in family_ids})
    payload = dict(existing or {})
    payload.update(
        {
            "schema_version": STORE_METADATA_SCHEMA_VERSION,
            "mode": mode,
            "compression": compression,
            "compression_level": int(compression_level),
            "family_lock_stripes": FAMILY_LOCK_STRIPES,
        }
    )
    if not preserve_existing_catalog or existing is None:
        payload.update(
            {
                "catalog_family_count": len(family_ids_set),
                "catalog_family_ids_sha256": _family_catalog_digest(family_ids_set),
                "catalog_updated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
        )
    if catalog_sources and not preserve_existing_catalog:
        payload["catalog_sources"] = [str(Path(path).resolve()) for path in catalog_sources]
    metadata_path = archive_root / STORE_METADATA_FILE
    temporary = metadata_path.with_name(f".{metadata_path.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, sort_keys=True, separators=(",", ":"))
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, metadata_path)
        _fsync_directory(archive_root)
    finally:
        if temporary.exists():
            temporary.unlink()
    _write_archive_readme(root)
    return payload


def _write_archive_readme(root: Path | str) -> None:
    resolved = Path(root).resolve()
    readme_path = resolved / "README_GENE_FAMILY_OUTPUTS.txt"
    content = (
        "GeneGalleon ZIP-backed gene-family outputs\n"
        "\n"
        "Finalized output sets are ordinary <subdirectory>.zip files in this directory.\n"
        "New or manually changed files remain in <subdirectory>/ and override ZIP members.\n"
        "While a run is active, immutable ZIP parts are visible below archives/<subdirectory>/.\n"
        "ARCHIVE_STATUS.tsv is a snapshot of the physical location of every logical output set.\n"
        "Internal indexes and deletion records are under .gg_store/; locks are under .gg_store_locks/.\n"
        "Do not edit or remove either internal tree.\n"
        "After manual file changes, run gg_gene_family_archive.sh refresh-status --root THIS_DIRECTORY.\n"
        "Use gg_gene_family_archive.sh list, verify, restore, delete, or finalize to manage files.\n"
    )
    temporary = readme_path.with_name(f".{readme_path.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, readme_path)
        _fsync_directory(resolved)
    finally:
        if temporary.exists():
            temporary.unlink()


def _manifest_modes(root: Path) -> Set[str]:
    modes: Set[str] = set()
    for zip_path in _physical_archive_paths(root):
        if zip_path.is_symlink() or zip_path.parent.is_symlink():
            raise ArchiveStoreError(f"Symlinked ZIP shards are not supported: {zip_path}")
        try:
            with zipfile.ZipFile(zip_path, "r") as archive:
                with archive.open(MANIFEST_MEMBER, "r") as handle:
                    manifest = json.load(handle)
            mode = str(manifest.get("mode", ""))
        except (OSError, KeyError, ValueError, zipfile.BadZipFile, json.JSONDecodeError) as exc:
            raise ArchiveStoreError(f"Failed to infer the gene-family mode from {zip_path}: {exc}") from exc
        if mode not in {"query2family", "orthogroup"}:
            raise ArchiveStoreError(f"ZIP shard has an invalid gene-family mode: {zip_path}")
        modes.add(mode)
    return modes


def resolve_archive_mode(
    root: Path,
    explicit_mode: Optional[str] = None,
    *,
    required: bool = True,
) -> Optional[str]:
    if explicit_mode is not None and explicit_mode not in {"query2family", "orthogroup"}:
        raise ValueError(f"Unsupported gene-family mode: {explicit_mode}")
    metadata = _read_store_metadata(root)
    observed_mode = None if metadata is None else str(metadata["mode"])
    marker = _read_storage_conversion_marker(Path(root).resolve())
    if observed_mode is None and marker is not None:
        observed_mode = str(marker["mode"])
    if observed_mode is None:
        modes = _manifest_modes(root)
        if len(modes) > 1:
            raise ArchiveStoreError("ZIP shards use mixed gene-family modes: " + ", ".join(sorted(modes)))
        if modes:
            observed_mode = next(iter(modes))
    if explicit_mode is not None and observed_mode is not None and explicit_mode != observed_mode:
        raise ArchiveStoreError(
            f"Archive metadata uses a different gene-family mode: requested={explicit_mode}, archived={observed_mode}"
        )
    resolved = explicit_mode or observed_mode
    if required and resolved is None:
        raise ValueError("--mode is required because this raw store has no archive metadata")
    return resolved


class ProgressReporter:
    """Emit line-buffered conversion progress and periodic heartbeats."""

    def __init__(self, interval_seconds: float = 30.0):
        self.interval_seconds = max(0.0, float(interval_seconds))
        self.started = time.monotonic()
        self._state: Dict[str, object] = {"phase": "starting"}
        self._lock = threading.Lock()
        self._stop = threading.Event()
        self._thread: Optional[threading.Thread] = None

    def _emit(self) -> None:
        with self._lock:
            state = dict(self._state)
        fields = [
            f"elapsed_seconds={int(time.monotonic() - self.started)}",
            *(f"{key}={state[key]}" for key in sorted(state)),
        ]
        print("progress\t" + "\t".join(fields), file=sys.stderr, flush=True)

    def start(self) -> None:
        self._emit()
        if self.interval_seconds <= 0:
            return

        def heartbeat() -> None:
            while not self._stop.wait(self.interval_seconds):
                self._emit()

        self._thread = threading.Thread(target=heartbeat, daemon=True)
        self._thread.start()

    def update(self, *, force: bool = False, **fields: object) -> None:
        with self._lock:
            self._state.update(fields)
        if force:
            self._emit()

    def close(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=max(1.0, self.interval_seconds + 1.0))
        self._emit()


@contextlib.contextmanager
def producer_read_lock(
    archive_root: Path,
    *,
    nonblocking: bool = False,
) -> Iterator[bool]:
    """Keep archive maintenance from replacing a source while it is open."""

    _validate_archive_root(archive_root)
    if not archive_root.is_dir():
        yield True
        return
    with _bucket_lock(
        archive_root / PRODUCER_LOCK_FILE,
        exclusive=False,
        nonblocking=nonblocking,
    ) as acquired:
        yield acquired


class GeneFamilyOutputStore:
    def __init__(self, root: Path | str, family_filter: Optional[str] = None):
        if family_filter is not None and not isinstance(family_filter, str):
            raise TypeError(
                f"family_filter must be a single family ID string or None, not {type(family_filter).__name__}"
            )
        self.root = Path(root).resolve()
        self.archive_root = _archive_state_root(self.root)
        self.payload_root = _archive_payload_root(self.root)
        _validate_archive_root(self.archive_root)
        _validate_archive_root(self.payload_root)
        self.family_filter = family_filter
        self._archives_loaded = False
        self._archived: Dict[str, Artifact] = {}
        self._tombstones: Dict[str, Tuple[int, str]] = {}
        self._tombstones_loaded = False
        self._family_states: Dict[str, Tuple[int, str, str]] = {}
        self._states_loaded = False
        self._zip_manifests: Dict[Path, dict] = {}
        self._index_generation = 0
        self._index_loaded = False
        self._subdir_archived: Dict[str, Dict[str, Artifact]] = {}
        self._index_catalog: Optional[dict] = None
        self._cache_epoch: Optional[str] = None
        self._index_update_active = False
        self._cache_lock = threading.RLock()
        self._zip_reader_cache: weakref.WeakKeyDictionary[
            threading.Thread,
            Tuple[Path, zipfile.ZipFile, int],
        ] = weakref.WeakKeyDictionary()

    def _assert_no_pending_index_update(self) -> None:
        marker = self.archive_root / INDEX_UPDATE_FILE
        if marker.is_symlink():
            raise ArchiveStoreError(f"Symlinked archive index update markers are not supported: {marker}")
        if marker.is_file() and not self._index_update_active:
            raise ArchiveStoreError(
                "An interrupted archive index update was detected; run the repair "
                "command before reading or modifying this store"
            )

    def _begin_index_update(self, *, recover_pending: bool = False) -> None:
        self.archive_root.mkdir(parents=True, exist_ok=True)
        marker = self.archive_root / INDEX_UPDATE_FILE
        if marker.is_symlink():
            raise ArchiveStoreError(f"Symlinked archive index update markers are not supported: {marker}")
        if marker.is_file() and not recover_pending:
            raise ArchiveStoreError(
                "An interrupted archive index update was detected; run the repair "
                "command before reading or modifying this store"
            )
        temporary = marker.with_name(f".{marker.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
        payload = {
            "schema_version": INDEX_SCHEMA_VERSION,
            "pid": os.getpid(),
            "hostname": os.uname().nodename,
            "created_ns": time.time_ns(),
        }
        try:
            with temporary.open("w", encoding="utf-8") as handle:
                json.dump(payload, handle, sort_keys=True, separators=(",", ":"))
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, marker)
            _fsync_directory(self.archive_root)
            self._index_update_active = True
        finally:
            if temporary.exists():
                temporary.unlink()

    def _finish_index_update(self) -> None:
        marker = self.archive_root / INDEX_UPDATE_FILE
        if marker.is_file() and not marker.is_symlink():
            marker.unlink()
            _fsync_directory(self.archive_root)
        self._index_update_active = False
        self._cache_epoch = self._read_index_epoch()

    @contextlib.contextmanager
    def _index_update(self, *, recover_pending: bool = False) -> Iterator[None]:
        self._begin_index_update(recover_pending=recover_pending)
        try:
            yield
        except BaseException:
            # Keep the durable marker so readers fail closed until repair has
            # rebuilt all denormalized index views from the ZIP manifests.
            self._index_update_active = False
            raise
        else:
            self._finish_index_update()

    def _load_tombstones(self) -> None:
        if self._tombstones_loaded:
            return
        tombstone_path = self.archive_root / TOMBSTONE_FILE
        if not tombstone_path.is_file():
            self._tombstones_loaded = True
            return
        if tombstone_path.is_symlink():
            raise ArchiveStoreError(f"Symlinked tombstone logs are not supported: {tombstone_path}")
        with tombstone_path.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                stripped = line.strip()
                if not stripped:
                    continue
                try:
                    record = json.loads(stripped)
                    logical_path = str(record["logical_path"])
                    generation = int(record["generation"])
                    operation = str(record["operation"])
                except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
                    raise ArchiveStoreError(
                        f"Invalid tombstone record at {tombstone_path}:{line_number}: {exc}"
                    ) from exc
                try:
                    tombstone_subdir, tombstone_name = logical_path.split("/", 1)
                    _safe_logical_path(tombstone_subdir, tombstone_name)
                except (ValueError, ArchiveStoreError) as exc:
                    raise ArchiveStoreError(
                        f"Unsafe tombstone path at {tombstone_path}:{line_number}: {logical_path!r}"
                    ) from exc
                if generation < 0:
                    raise ArchiveStoreError(f"Negative tombstone generation at {tombstone_path}:{line_number}")
                if operation not in {"delete", "undelete"}:
                    raise ArchiveStoreError(
                        f"Unsupported tombstone operation at {tombstone_path}:{line_number}: {operation}"
                    )
                previous = self._tombstones.get(logical_path)
                if previous is None or generation >= previous[0]:
                    self._tombstones[logical_path] = (generation, operation)
        self._tombstones_loaded = True

    def _load_legacy_family_states(self) -> Dict[str, Tuple[int, str, str]]:
        states: Dict[str, Tuple[int, str, str]] = {}
        state_path = self.archive_root / FAMILY_STATE_FILE
        if not state_path.is_file():
            return states
        if state_path.is_symlink():
            raise ArchiveStoreError(f"Symlinked family-state logs are not supported: {state_path}")
        with state_path.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                stripped = line.strip()
                if not stripped:
                    continue
                try:
                    record = json.loads(stripped)
                    family_id = str(record["family_id"])
                    generation = int(record["generation"])
                    status = str(record["status"])
                    run_token = str(record.get("run_token", ""))
                except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
                    raise ArchiveStoreError(
                        f"Invalid family-state record at {state_path}:{line_number}: {exc}"
                    ) from exc
                if not family_id:
                    raise ArchiveStoreError(f"Empty family ID in state record at {state_path}:{line_number}")
                if status not in {"running", "complete", "failed"}:
                    raise ArchiveStoreError(f"Unsupported family state at {state_path}:{line_number}: {status}")
                previous = states.get(family_id)
                if previous is None or generation >= previous[0]:
                    states[family_id] = (generation, status, run_token)
        return states

    def _state_bucket_path(self, family_id: str) -> Path:
        return self.archive_root / FAMILY_STATE_DIR_NAME / f"{_family_index_bucket(family_id)}.json"

    def _read_state_bucket(self, path: Path) -> Dict[str, Tuple[int, str, str]]:
        if not path.is_file():
            return {}
        if path.is_symlink():
            raise ArchiveStoreError(f"Symlinked family-state indexes are not supported: {path}")
        try:
            with path.open("r", encoding="utf-8") as handle:
                payload = json.load(handle)
            if not isinstance(payload, dict):
                raise ArchiveStoreError(f"Family-state index is not a JSON object: {path}")
            if int(payload.get("schema_version", -1)) != INDEX_SCHEMA_VERSION:
                raise ArchiveStoreError(
                    f"Unsupported family-state index schema in {path}: {payload.get('schema_version')!r}"
                )
            records = payload.get("families")
            if not isinstance(records, dict):
                raise ArchiveStoreError(f"Family-state index has no family mapping: {path}")
            states: Dict[str, Tuple[int, str, str]] = {}
            for family_id, record in records.items():
                if not isinstance(record, dict):
                    raise ArchiveStoreError(f"Family-state index has an invalid record: {family_id}")
                status = str(record["status"])
                if not str(family_id):
                    raise ArchiveStoreError(f"Family-state index has an empty family ID: {path}")
                if status not in {"running", "complete", "failed"}:
                    raise ArchiveStoreError(f"Family-state index has an invalid status for {family_id}: {status}")
                generation = int(record["generation"])
                if generation < 0:
                    raise ArchiveStoreError(f"Family-state index has a negative generation for {family_id}: {path}")
                states[str(family_id)] = (
                    generation,
                    status,
                    str(record.get("run_token", "")),
                )
            return states
        except ArchiveStoreError:
            raise
        except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
            raise ArchiveStoreError(f"Failed to read family-state index {path}: {exc}") from exc

    def _write_state_bucket(
        self,
        path: Path,
        states: Dict[str, Tuple[int, str, str]],
    ) -> None:
        if path.parent.is_symlink():
            raise ArchiveStoreError(f"Symlinked family-state index directories are not supported: {path.parent}")
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "schema_version": INDEX_SCHEMA_VERSION,
            "families": {
                family_id: {
                    "generation": generation,
                    "status": status,
                    "run_token": run_token,
                }
                for family_id, (generation, status, run_token) in sorted(states.items())
            },
        }
        temporary = path.with_name(f".{path.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
        try:
            with temporary.open("w", encoding="utf-8") as handle:
                json.dump(payload, handle, sort_keys=True, separators=(",", ":"))
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, path)
            _fsync_directory(path.parent)
        finally:
            if temporary.exists():
                temporary.unlink()

    def _ensure_state_index(self) -> None:
        state_dir = self.archive_root / FAMILY_STATE_DIR_NAME
        if state_dir.is_symlink():
            raise ArchiveStoreError(f"Symlinked family-state index directories are not supported: {state_dir}")
        marker = state_dir / ".complete"
        if marker.is_symlink():
            raise ArchiveStoreError(f"Symlinked family-state completion markers are not supported: {marker}")
        if marker.is_file():
            return
        with archive_lock(self.archive_root):
            if marker.is_file():
                return
            legacy_states = self._load_legacy_family_states()
            by_bucket: Dict[str, Dict[str, Tuple[int, str, str]]] = {}
            for family_id, state in legacy_states.items():
                by_bucket.setdefault(_family_index_bucket(family_id), {})[family_id] = state
            state_dir.mkdir(parents=True, exist_ok=True)
            for bucket, states in by_bucket.items():
                self._write_state_bucket(state_dir / f"{bucket}.json", states)
            temporary = marker.with_name(f".{marker.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
            try:
                with temporary.open("w", encoding="utf-8") as handle:
                    handle.write("1\n")
                    handle.flush()
                    os.fsync(handle.fileno())
                os.replace(temporary, marker)
                _fsync_directory(state_dir)
            finally:
                if temporary.exists():
                    temporary.unlink()

    def _load_family_states(self) -> None:
        if self._states_loaded:
            return
        state_dir = self.archive_root / FAMILY_STATE_DIR_NAME
        if state_dir.is_symlink():
            raise ArchiveStoreError(f"Symlinked family-state index directories are not supported: {state_dir}")
        marker = state_dir / ".complete"
        if marker.is_symlink():
            raise ArchiveStoreError(f"Symlinked family-state completion markers are not supported: {marker}")
        if marker.is_file():
            if self.family_filter is None:
                state_paths = sorted(state_dir.glob("*.json"))
            else:
                state_paths = [self._state_bucket_path(self.family_filter)]
            for state_path in state_paths:
                for family_id, state in self._read_state_bucket(state_path).items():
                    if self.family_filter is None or family_id == self.family_filter:
                        self._family_states[family_id] = state
        else:
            legacy_states = self._load_legacy_family_states()
            if self.family_filter is None:
                self._family_states.update(legacy_states)
            elif self.family_filter in legacy_states:
                self._family_states[self.family_filter] = legacy_states[self.family_filter]
        self._states_loaded = True
        self._cache_epoch = self._read_index_epoch()

    def _artifact_from_index_record(self, logical_path: str, record: dict) -> Artifact:
        try:
            subdir, name = logical_path.split("/", 1)
        except ValueError as exc:
            raise ArchiveStoreError(f"Archive index contains an invalid logical path: {logical_path!r}") from exc
        _safe_logical_path(subdir, name)
        relative_zip = Path(str(record["zip_path"]))
        zip_location = str(record.get("zip_location", "active"))
        unsafe = relative_zip.is_absolute() or ".." in relative_zip.parts
        if zip_location == "active":
            unsafe = unsafe or len(relative_zip.parts) != 2 or relative_zip.parent.name != subdir
            zip_path = self.payload_root / relative_zip
        elif zip_location == "final":
            unsafe = unsafe or len(relative_zip.parts) != 1 or relative_zip.name != f"{subdir}.zip"
            zip_path = self.root / relative_zip
        else:
            unsafe = True
            zip_path = self.root / relative_zip
        if unsafe or relative_zip.suffix.lower() != ".zip":
            raise ArchiveStoreError(f"Archive index contains an unsafe ZIP path: {relative_zip}")
        artifact = Artifact(
            logical_path=logical_path,
            subdir=subdir,
            name=name,
            generation=int(record["generation"]),
            zip_path=zip_path,
            member_name=str(record["member_name"]),
            size=int(record["size"]),
            crc=int(record["crc"]),
            sha256=str(record.get("sha256", "")) or None,
            mtime_ns=(int(record["mtime_ns"]) if record.get("mtime_ns") is not None else None),
            family_id=(str(record["family_id"]) if record.get("family_id") not in {None, ""} else None),
        )
        if artifact.member_name != logical_path:
            raise ArchiveStoreError(f"Archive index member path mismatch: {artifact.member_name!r} != {logical_path!r}")
        if artifact.generation < 0 or (artifact.size or 0) < 0 or (artifact.crc or 0) < 0:
            raise ArchiveStoreError(f"Archive index contains invalid numeric metadata: {logical_path}")
        if artifact.sha256 is not None and re.fullmatch(r"[0-9a-f]{64}", artifact.sha256) is None:
            raise ArchiveStoreError(f"Archive index contains an invalid SHA256 digest: {logical_path}")
        return artifact

    def _read_index_epoch(self) -> str:
        epoch_path = self.archive_root / INDEX_EPOCH_FILE
        if epoch_path.is_symlink():
            raise ArchiveStoreError(f"Symlinked archive index epochs are not supported: {epoch_path}")
        epoch = ""
        if epoch_path.is_file():
            try:
                epoch = epoch_path.read_text(encoding="utf-8").strip()
            except OSError as exc:
                raise ArchiveStoreError(f"Failed to read archive index epoch {epoch_path}: {exc}") from exc
        if not epoch:
            index_dir = self.archive_root / INDEX_DIR_NAME
            try:
                epoch = f"legacy:{index_dir.stat().st_mtime_ns}" if index_dir.is_dir() else "none"
            except OSError:
                epoch = "none"
        # The extra signatures close the small crash window between replacing
        # metadata and replacing index.epoch. Inode is included because the
        # JSON files are atomically replaced and can retain size/mtime values.
        signatures: List[str] = []
        for name in (
            INDEX_CATALOG_FILE,
            TOMBSTONE_FILE,
            INDEX_UPDATE_FILE,
            FAMILY_STATE_DIR_NAME,
        ):
            path = self.archive_root / name
            try:
                stat_result = path.lstat()
            except FileNotFoundError:
                signatures.append(f"{name}:missing")
            except OSError as exc:
                raise ArchiveStoreError(f"Failed to stat archive metadata {path}: {exc}") from exc
            else:
                signatures.append(
                    f"{name}:{stat_result.st_dev}:{stat_result.st_ino}:{stat_result.st_size}:{stat_result.st_mtime_ns}"
                )
        return "|".join([epoch, *signatures])

    def _load_index_catalog(self) -> Optional[dict]:
        self._assert_no_pending_index_update()
        if self._index_catalog is not None:
            return self._index_catalog
        catalog_path = self.archive_root / INDEX_CATALOG_FILE
        if not catalog_path.is_file():
            return None
        if catalog_path.is_symlink():
            raise ArchiveStoreError(f"Symlinked archive index catalogs are not supported: {catalog_path}")
        index_dir = self.archive_root / INDEX_DIR_NAME
        subdir_index_dir = self.archive_root / SUBDIR_INDEX_DIR_NAME
        if index_dir.is_symlink() or subdir_index_dir.is_symlink():
            raise ArchiveStoreError(
                "Symlinked archive index directories are not supported: "
                f"{index_dir if index_dir.is_symlink() else subdir_index_dir}"
            )
        try:
            with catalog_path.open("r", encoding="utf-8") as handle:
                catalog = json.load(handle)
            if not isinstance(catalog, dict):
                raise ArchiveStoreError(f"Archive index catalog is not a JSON object: {catalog_path}")
            if int(catalog.get("schema_version", -1)) != INDEX_SCHEMA_VERSION:
                raise ArchiveStoreError(
                    f"Unsupported archive index catalog schema in {catalog_path}: {catalog.get('schema_version')!r}"
                )
            buckets = catalog.get("buckets")
            if (
                not isinstance(buckets, list)
                or not all(isinstance(bucket, str) and re.fullmatch(r"[0-9a-f]{2}", bucket) for bucket in buckets)
                or len(buckets) != len(set(buckets))
            ):
                raise ArchiveStoreError(f"Archive index catalog has invalid buckets: {catalog_path}")
            for bucket in buckets:
                bucket_path = self.archive_root / INDEX_DIR_NAME / f"{bucket}.json"
                if not bucket_path.is_file():
                    raise ArchiveStoreError(
                        f"ZIP members are absent from the archive index because the "
                        f"catalog references a missing bucket: {bucket_path}"
                    )
            physical_buckets = {
                path.stem for path in (self.archive_root / INDEX_DIR_NAME).glob("*.json") if path.is_file()
            }
            if physical_buckets != set(buckets):
                raise ArchiveStoreError(f"Archive index catalog and bucket files disagree: {catalog_path}")
            subdir_indexes = catalog.get("subdir_indexes", {})
            if not isinstance(subdir_indexes, dict):
                raise ArchiveStoreError(f"Archive index catalog has invalid subdirectory indexes: {catalog_path}")
            subdirs = catalog.get("subdirs")
            if (
                not isinstance(subdirs, list)
                or not all(isinstance(subdir, str) for subdir in subdirs)
                or len(subdirs) != len(set(subdirs))
            ):
                raise ArchiveStoreError(f"Archive index catalog has invalid subdirectories: {catalog_path}")
            for subdir in subdirs:
                _safe_logical_path(subdir, "__catalog__")
            if set(subdir_indexes) != set(subdirs):
                raise ArchiveStoreError(f"Archive index catalog subdirectory mappings disagree: {catalog_path}")
            subdir_counts = catalog.get("subdir_counts")
            if (
                not isinstance(subdir_counts, dict)
                or set(subdir_counts) != set(subdirs)
                or not all(
                    isinstance(count, int) and not isinstance(count, bool) and count >= 0
                    for count in subdir_counts.values()
                )
            ):
                raise ArchiveStoreError(f"Archive index catalog has invalid subdirectory counts: {catalog_path}")
            artifact_count = catalog.get("artifact_count")
            if (
                not isinstance(artifact_count, int)
                or isinstance(artifact_count, bool)
                or artifact_count < 0
                or artifact_count != sum(subdir_counts.values())
            ):
                raise ArchiveStoreError(f"Archive index catalog has an invalid artifact count: {catalog_path}")
            for subdir, file_name in subdir_indexes.items():
                if file_name != _subdir_index_name(str(subdir)):
                    raise ArchiveStoreError(f"Archive index catalog has an invalid subdirectory index for {subdir}")
                subdir_path = self.archive_root / SUBDIR_INDEX_DIR_NAME / str(file_name)
                if not subdir_path.is_file():
                    raise ArchiveStoreError(
                        f"Archive index catalog references a missing subdirectory index: {subdir_path}"
                    )
            physical_subdir_indexes = {
                path.name for path in (self.archive_root / SUBDIR_INDEX_DIR_NAME).glob("*.json") if path.is_file()
            }
            if physical_subdir_indexes != set(subdir_indexes.values()):
                raise ArchiveStoreError(f"Archive catalog and subdirectory index files disagree: {catalog_path}")
            self._index_catalog = catalog
            return catalog
        except ArchiveStoreError:
            raise
        except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
            raise ArchiveStoreError(f"Failed to read archive index catalog {catalog_path}: {exc}") from exc

    def _write_index_catalog(self, artifacts: Dict[str, Artifact]) -> None:
        catalog_path = self.archive_root / INDEX_CATALOG_FILE
        buckets = sorted({_family_index_bucket(artifact.family_id) for artifact in artifacts.values()})
        subdirs = sorted({artifact.subdir for artifact in artifacts.values()})
        subdir_counts: Dict[str, int] = {}
        for artifact in artifacts.values():
            subdir_counts[artifact.subdir] = subdir_counts.get(artifact.subdir, 0) + 1
        payload = {
            "schema_version": INDEX_SCHEMA_VERSION,
            "artifact_count": len(artifacts),
            "buckets": buckets,
            "subdirs": subdirs,
            "subdir_indexes": {subdir: _subdir_index_name(subdir) for subdir in subdirs},
            "subdir_counts": subdir_counts,
        }
        temporary = catalog_path.with_name(f".{catalog_path.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
        try:
            with temporary.open("w", encoding="utf-8") as handle:
                json.dump(payload, handle, sort_keys=True, separators=(",", ":"))
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, catalog_path)
        finally:
            if temporary.exists():
                temporary.unlink()
        self._index_catalog = payload

    def _write_index_epoch(self) -> None:
        epoch_path = self.archive_root / INDEX_EPOCH_FILE
        temporary = epoch_path.with_name(f".{epoch_path.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
        epoch = f"{time.time_ns()}-{os.getpid()}-{uuid.uuid4().hex}"
        try:
            with temporary.open("w", encoding="utf-8") as handle:
                handle.write(epoch + "\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, epoch_path)
            _fsync_directory(self.archive_root)
        finally:
            if temporary.exists():
                temporary.unlink()
        self._cache_epoch = self._read_index_epoch()

    def _write_json_index_payload(self, path: Path, payload: dict) -> None:
        if path.parent.is_symlink():
            raise ArchiveStoreError(f"Symlinked archive index directories are not supported: {path.parent}")
        encoded = (json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")
        try:
            if path.is_file() and path.read_bytes() == encoded:
                return
        except OSError:
            pass
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(f".{path.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
        try:
            with temporary.open("wb") as handle:
                handle.write(encoded)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, path)
        finally:
            if temporary.exists():
                temporary.unlink()

    def _read_family_bucket_artifacts(self, bucket: str) -> Dict[str, Artifact]:
        if re.fullmatch(r"[0-9a-f]{2}", bucket) is None:
            raise ArchiveStoreError(f"Invalid archive index bucket: {bucket!r}")
        index_path = self.archive_root / INDEX_DIR_NAME / f"{bucket}.json"
        if not index_path.is_file():
            return {}
        if index_path.is_symlink():
            raise ArchiveStoreError(f"Symlinked archive indexes are not supported: {index_path}")
        try:
            with index_path.open("r", encoding="utf-8") as handle:
                payload = json.load(handle)
            if not isinstance(payload, dict):
                raise ArchiveStoreError(f"Archive index bucket is not a JSON object: {index_path}")
            if int(payload.get("schema_version", -1)) != INDEX_SCHEMA_VERSION:
                raise ArchiveStoreError(f"Unsupported GeneGalleon archive index schema in {index_path}")
            records = payload.get("artifacts")
            if not isinstance(records, dict):
                raise ArchiveStoreError(f"Archive index has no artifact mapping: {index_path}")
            artifacts: Dict[str, Artifact] = {}
            for logical_path, record in records.items():
                if not isinstance(record, dict):
                    raise ArchiveStoreError(f"Archive index has an invalid artifact record: {logical_path}")
                artifact = self._artifact_from_index_record(
                    str(logical_path),
                    record,
                )
                if _family_index_bucket(artifact.family_id) != bucket:
                    raise ArchiveStoreError(
                        f"Archive index record is in the wrong family bucket: {logical_path} in {index_path}"
                    )
                artifacts[str(logical_path)] = artifact
            return artifacts
        except ArchiveStoreError:
            raise
        except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
            raise ArchiveStoreError(f"Failed to read archive index bucket {index_path}: {exc}") from exc

    def _merge_index_subdirs(
        self,
        updated_subdirs: Dict[str, Dict[str, Artifact]],
    ) -> None:
        if not updated_subdirs:
            return
        with self._index_update():
            self._merge_index_subdirs_uncommitted(updated_subdirs)

    def _merge_index_subdirs_uncommitted(
        self,
        updated_subdirs: Dict[str, Dict[str, Artifact]],
    ) -> None:
        catalog = self._load_index_catalog()
        if catalog is None:
            catalog = {
                "schema_version": INDEX_SCHEMA_VERSION,
                "artifact_count": 0,
                "buckets": [],
                "subdirs": [],
                "subdir_indexes": {},
                "subdir_counts": {},
            }
        catalog = dict(catalog)
        subdir_counts = {str(subdir): int(count) for subdir, count in dict(catalog.get("subdir_counts", {})).items()}
        affected_by_bucket: Dict[str, Dict[str, Artifact]] = {}
        affected_buckets: Set[str] = set()
        subdir_index_dir = self.archive_root / SUBDIR_INDEX_DIR_NAME
        for subdir, artifacts in updated_subdirs.items():
            previous_artifacts = self._load_subdir_artifacts(subdir)
            affected_buckets.update(
                _family_index_bucket(artifact.family_id) for artifact in previous_artifacts.values()
            )
            if artifacts:
                subdir_counts[subdir] = len(artifacts)
                payload = {
                    "schema_version": INDEX_SCHEMA_VERSION,
                    "subdir": subdir,
                    "artifacts": {
                        logical_path: self._index_record(artifact)
                        for logical_path, artifact in sorted(artifacts.items())
                    },
                }
                self._write_json_index_payload(
                    subdir_index_dir / _subdir_index_name(subdir),
                    payload,
                )
            else:
                subdir_counts.pop(subdir, None)
                subdir_path = subdir_index_dir / _subdir_index_name(subdir)
                if subdir_path.is_file() and not subdir_path.is_symlink():
                    subdir_path.unlink()
            for logical_path, artifact in artifacts.items():
                bucket = _family_index_bucket(artifact.family_id)
                affected_buckets.add(bucket)
                affected_by_bucket.setdefault(bucket, {})[logical_path] = artifact
        index_dir = self.archive_root / INDEX_DIR_NAME
        nonempty_affected_buckets: Set[str] = set()
        for bucket in sorted(affected_buckets):
            replacements = affected_by_bucket.get(bucket, {})
            bucket_artifacts = self._read_family_bucket_artifacts(bucket)
            changed_subdirs = set(updated_subdirs)
            bucket_artifacts = {
                logical_path: artifact
                for logical_path, artifact in bucket_artifacts.items()
                if artifact.subdir not in changed_subdirs
            }
            bucket_artifacts.update(replacements)
            bucket_path = index_dir / f"{bucket}.json"
            if not bucket_artifacts:
                if bucket_path.is_file() and not bucket_path.is_symlink():
                    bucket_path.unlink()
                continue
            nonempty_affected_buckets.add(bucket)
            payload = {
                "schema_version": INDEX_SCHEMA_VERSION,
                "generation": max(
                    (artifact.generation for artifact in bucket_artifacts.values()),
                    default=0,
                ),
                "artifacts": {
                    logical_path: self._index_record(artifact)
                    for logical_path, artifact in sorted(bucket_artifacts.items())
                },
            }
            self._write_json_index_payload(bucket_path, payload)
        buckets = set(str(bucket) for bucket in catalog.get("buckets", []))
        buckets.difference_update(affected_buckets)
        buckets.update(nonempty_affected_buckets)
        subdirs = set(str(subdir) for subdir in catalog.get("subdirs", []))
        for subdir, artifacts in updated_subdirs.items():
            if artifacts:
                subdirs.add(subdir)
            else:
                subdirs.discard(subdir)
        catalog.update(
            {
                "schema_version": INDEX_SCHEMA_VERSION,
                "artifact_count": sum(subdir_counts.values()),
                "buckets": sorted(buckets),
                "subdirs": sorted(subdirs),
                "subdir_indexes": {subdir: _subdir_index_name(subdir) for subdir in sorted(subdirs)},
                "subdir_counts": subdir_counts,
            }
        )
        catalog_path = self.archive_root / INDEX_CATALOG_FILE
        self._write_json_index_payload(catalog_path, catalog)
        _fsync_directory(index_dir)
        _fsync_directory(subdir_index_dir)
        self._index_catalog = catalog
        self._write_index_epoch()
        self._reset_cache()

    def _load_index(self) -> bool:
        self._assert_no_pending_index_update()
        index_dir = self.archive_root / INDEX_DIR_NAME
        legacy_index_path = self.archive_root / INDEX_FILE
        if index_dir.is_dir():
            if index_dir.is_symlink():
                raise ArchiveStoreError(f"Symlinked archive index directories are not supported: {index_dir}")
            catalog = self._load_index_catalog()
            if self.family_filter is None:
                index_paths = sorted(index_dir.glob("*.json"))
            else:
                buckets = {
                    _family_index_bucket(self.family_filter),
                    _family_index_bucket(None),
                }
                index_paths = [
                    index_dir / f"{bucket}.json"
                    for bucket in sorted(buckets)
                    if (index_dir / f"{bucket}.json").is_file()
                ]
            if catalog is not None and self.family_filter is not None:
                catalog_buckets = set(catalog["buckets"])
                expected = {
                    _family_index_bucket(self.family_filter),
                    _family_index_bucket(None),
                }
                index_paths = [index_dir / f"{bucket}.json" for bucket in sorted(expected & catalog_buckets)]
            if not index_paths and catalog is None and _physical_archive_paths(self.root):
                return False
        elif legacy_index_path.is_file():
            index_paths = [legacy_index_path]
        else:
            return False
        try:
            for index_path in index_paths:
                if index_path.is_symlink():
                    raise ArchiveStoreError(f"Symlinked archive indexes are not supported: {index_path}")
                with index_path.open("r", encoding="utf-8") as handle:
                    payload = json.load(handle)
                if not isinstance(payload, dict):
                    raise ArchiveStoreError(f"Archive index is not a JSON object: {index_path}")
                if int(payload.get("schema_version", -1)) != INDEX_SCHEMA_VERSION:
                    raise ArchiveStoreError(
                        f"Unsupported GeneGalleon archive index schema in {index_path}: "
                        f"{payload.get('schema_version')!r}"
                    )
                self._index_generation = max(
                    self._index_generation,
                    int(payload.get("generation", 0)),
                )
                artifacts = payload.get("artifacts")
                if not isinstance(artifacts, dict):
                    raise ArchiveStoreError(f"Archive index has no artifact mapping: {index_path}")
                for logical_path, record in artifacts.items():
                    if not isinstance(record, dict):
                        raise ArchiveStoreError(f"Archive index has an invalid artifact record: {logical_path}")
                    artifact = self._artifact_from_index_record(str(logical_path), record)
                    if index_path.parent == index_dir and index_path.stem != _family_index_bucket(artifact.family_id):
                        raise ArchiveStoreError(
                            f"Archive index record is in the wrong family bucket: {logical_path} in {index_path}"
                        )
                    if (
                        self.family_filter is not None
                        and artifact.family_id is not None
                        and artifact.family_id != self.family_filter
                    ):
                        continue
                    if artifact.logical_path in self._archived:
                        raise ArchiveStoreError(
                            f"Archive index contains a duplicate logical path: {artifact.logical_path}"
                        )
                    self._archived[artifact.logical_path] = artifact
            self._index_loaded = True
            return True
        except ArchiveStoreError:
            raise
        except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
            raise ArchiveStoreError(f"Failed to read GeneGalleon archive index: {exc}") from exc

    def family_state(self, family_id: str) -> Optional[str]:
        self._refresh_if_index_changed()
        if self.family_filter is not None:
            state_dir = self.archive_root / FAMILY_STATE_DIR_NAME
            marker = state_dir / ".complete"
            if state_dir.is_symlink() or marker.is_symlink():
                raise ArchiveStoreError(
                    f"Symlinked family-state metadata is not supported: "
                    f"{state_dir if state_dir.is_symlink() else marker}"
                )
            if marker.is_file():
                state = self._read_state_bucket(self._state_bucket_path(family_id)).get(family_id)
                return None if state is None else state[1]
        self._load_family_states()
        state = self._family_states.get(family_id)
        return None if state is None else state[1]

    def _read_manifest(self, zip_path: Path) -> dict:
        try:
            with zipfile.ZipFile(zip_path, "r") as archive:
                try:
                    manifest_handle = archive.open(MANIFEST_MEMBER, "r")
                except KeyError as exc:
                    raise ArchiveStoreError(
                        f"Registered GeneGalleon manifest was not found in ZIP: {zip_path}"
                    ) from exc
                with manifest_handle:
                    manifest = json.load(manifest_handle)
                if not isinstance(manifest, dict):
                    raise ArchiveStoreError(f"Archive manifest is not a JSON object: {zip_path}")
                if int(manifest.get("schema_version", -1)) != ARCHIVE_SCHEMA_VERSION:
                    raise ArchiveStoreError(
                        f"Unsupported GeneGalleon archive schema in {zip_path}: {manifest.get('schema_version')!r}"
                    )
                try:
                    manifest_generation = int(manifest["generation"])
                except (KeyError, TypeError, ValueError) as exc:
                    raise ArchiveStoreError(f"Archive manifest has an invalid generation: {zip_path}") from exc
                if manifest_generation < 0:
                    raise ArchiveStoreError(f"Archive manifest has a negative generation: {zip_path}")
                if manifest.get("mode") not in {"query2family", "orthogroup"}:
                    raise ArchiveStoreError(
                        f"Archive manifest has an invalid mode in {zip_path}: {manifest.get('mode')!r}"
                    )
                members = manifest.get("members")
                if not isinstance(members, list):
                    raise ArchiveStoreError(f"Archive manifest has no member list: {zip_path}")
                if not all(isinstance(member, dict) for member in members):
                    raise ArchiveStoreError(f"Archive manifest contains an invalid member record: {zip_path}")
                declared_subdir = str(manifest.get("subdir", ""))
                active_location = declared_subdir == zip_path.parent.name
                final_location = zip_path.parent.resolve() == self.root and zip_path.name == f"{declared_subdir}.zip"
                if not active_location and not final_location:
                    raise ArchiveStoreError(
                        f"Archive location differs from its manifest in {zip_path}: {declared_subdir!r}"
                    )
                zip_infos = archive.infolist()
                zip_member_names = [info.filename for info in zip_infos]
                if len(zip_member_names) != len(set(zip_member_names)):
                    raise ArchiveStoreError(f"ZIP contains duplicate member names: {zip_path}")
                manifest_member_names = [str(member.get("member_name", "")) for member in members]
                if len(manifest_member_names) != len(set(manifest_member_names)):
                    raise ArchiveStoreError(f"Manifest contains duplicate member names: {zip_path}")
                zip_info_by_name = {info.filename: info for info in zip_infos}
                expected_member_names = set(manifest_member_names)
                expected_member_names.add(MANIFEST_MEMBER)
                unexpected_member_names = sorted(set(zip_info_by_name) - expected_member_names)
                missing_member_names = sorted(expected_member_names - set(zip_info_by_name))
                if unexpected_member_names or missing_member_names:
                    raise ArchiveStoreError(
                        f"ZIP member inventory differs from its manifest in {zip_path}: "
                        f"unexpected={unexpected_member_names}, missing={missing_member_names}"
                    )
                for member in members:
                    logical_path = str(member.get("logical_path", ""))
                    member_name = str(member.get("member_name", ""))
                    if member_name != logical_path:
                        raise ArchiveStoreError(
                            f"Manifest member path mismatch in {zip_path}: {member_name!r} != {logical_path!r}"
                        )
                    try:
                        logical_subdir, logical_name = logical_path.split("/", 1)
                    except ValueError as exc:
                        raise ArchiveStoreError(
                            f"Manifest contains an invalid logical path in {zip_path}: {logical_path!r}"
                        ) from exc
                    _safe_logical_path(logical_subdir, logical_name)
                    if logical_subdir != declared_subdir:
                        raise ArchiveStoreError(
                            f"Manifest member is in the wrong logical subdirectory in {zip_path}: {logical_path}"
                        )
                    info = zip_info_by_name.get(member_name)
                    if info is None:
                        raise ArchiveStoreError(f"Manifest member is missing from ZIP {zip_path}: {member_name}")
                    if int(member.get("size", -1)) != int(info.file_size):
                        raise ArchiveStoreError(f"Manifest size mismatch in {zip_path}: {member_name}")
                    if int(member.get("crc", -1)) != int(info.CRC):
                        raise ArchiveStoreError(f"Manifest CRC mismatch in {zip_path}: {member_name}")
                    member_generation = int(member.get("generation", manifest_generation))
                    if not 0 <= member_generation <= manifest_generation:
                        raise ArchiveStoreError(
                            f"Manifest member has an invalid generation in {zip_path}: {member_name}"
                        )
                    sha256 = str(member.get("sha256", ""))
                    if re.fullmatch(r"[0-9a-f]{64}", sha256) is None:
                        raise ArchiveStoreError(
                            f"Manifest member has an invalid SHA256 digest in {zip_path}: {member_name}"
                        )
                return manifest
        except (
            OSError,
            zipfile.BadZipFile,
            json.JSONDecodeError,
            KeyError,
            TypeError,
            UnicodeError,
            ValueError,
        ) as exc:
            raise ArchiveStoreError(f"Failed to read GeneGalleon ZIP archive {zip_path}: {exc}") from exc

    def _load_archives(self) -> None:
        with self._cache_lock:
            self._load_archives_unlocked()

    def _load_archives_unlocked(self) -> None:
        if self._archives_loaded:
            return
        self._load_tombstones()
        index_loaded = self._load_index()
        if not index_loaded:
            zip_paths = _physical_archive_paths(self.root)
            for zip_path in zip_paths:
                if zip_path.is_symlink() or zip_path.parent.is_symlink():
                    raise ArchiveStoreError(f"Symlinked ZIP shards are not supported: {zip_path}")
                manifest = self._read_manifest(zip_path)
                self._zip_manifests[zip_path] = manifest
                generation = int(manifest["generation"])
                for member in manifest["members"]:
                    logical_path = str(member["logical_path"])
                    subdir, name = logical_path.split("/", 1)
                    _safe_logical_path(subdir, name)
                    source_signature = member.get("source_signature")
                    mtime_ns = member.get("mtime_ns")
                    if mtime_ns is None and isinstance(source_signature, list) and len(source_signature) >= 4:
                        mtime_ns = source_signature[3]
                    candidate = Artifact(
                        logical_path=logical_path,
                        subdir=subdir,
                        name=name,
                        generation=int(member.get("generation", generation)),
                        zip_path=zip_path,
                        member_name=str(member["member_name"]),
                        size=int(member["size"]),
                        crc=int(member["crc"]),
                        sha256=str(member.get("sha256", "")) or None,
                        mtime_ns=int(mtime_ns) if mtime_ns is not None else None,
                        family_id=(str(member["family_id"]) if member.get("family_id") not in {None, ""} else None),
                    )
                    previous = self._archived.get(logical_path)
                    if previous is None or candidate.generation >= previous.generation:
                        self._archived[logical_path] = candidate
            if self._archived:
                self._index_generation = max(
                    self._index_generation,
                    max(artifact.generation for artifact in self._archived.values()),
                )
        if self._tombstones:
            self._index_generation = max(
                self._index_generation,
                max(generation for generation, _ in self._tombstones.values()),
            )
        self._cache_epoch = self._read_index_epoch()
        self._archives_loaded = True

    def _load_subdir_artifacts(self, subdir: str) -> Dict[str, Artifact]:
        with self._cache_lock:
            cached = self._subdir_archived.get(subdir)
            if cached is not None:
                return cached
            self._load_tombstones()
            catalog = self._load_index_catalog()
            if catalog is None:
                self._load_archives_unlocked()
                artifacts = {
                    logical_path: artifact
                    for logical_path, artifact in self._archived.items()
                    if artifact.subdir == subdir
                }
                self._subdir_archived[subdir] = artifacts
                return artifacts
            file_name = catalog.get("subdir_indexes", {}).get(subdir)
            if file_name is None:
                self._subdir_archived[subdir] = {}
                return self._subdir_archived[subdir]
            index_path = self.archive_root / SUBDIR_INDEX_DIR_NAME / str(file_name)
            if index_path.is_symlink():
                raise ArchiveStoreError(f"Symlinked subdirectory indexes are not supported: {index_path}")
            try:
                with index_path.open("r", encoding="utf-8") as handle:
                    payload = json.load(handle)
                if not isinstance(payload, dict):
                    raise ArchiveStoreError(f"Subdirectory index is not a JSON object: {index_path}")
                if (
                    int(payload.get("schema_version", -1)) != INDEX_SCHEMA_VERSION
                    or str(payload.get("subdir", "")) != subdir
                ):
                    raise ArchiveStoreError(f"Invalid GeneGalleon subdirectory index: {index_path}")
                records = payload.get("artifacts")
                if not isinstance(records, dict):
                    raise ArchiveStoreError(f"Subdirectory index has no artifact mapping: {index_path}")
                artifacts: Dict[str, Artifact] = {}
                for logical_path, record in records.items():
                    if not isinstance(record, dict):
                        raise ArchiveStoreError(f"Subdirectory index has an invalid artifact record: {logical_path}")
                    artifact = self._artifact_from_index_record(
                        str(logical_path),
                        record,
                    )
                    if artifact.subdir != subdir:
                        raise ArchiveStoreError(
                            f"Subdirectory index record is in the wrong subdirectory: {logical_path} in {index_path}"
                        )
                    artifacts[str(logical_path)] = artifact
            except ArchiveStoreError:
                raise
            except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
                raise ArchiveStoreError(f"Failed to read GeneGalleon subdirectory index {index_path}: {exc}") from exc
            self._subdir_archived[subdir] = artifacts
            self._cache_epoch = self._read_index_epoch()
            return artifacts

    def _logical_subdirs_unlocked(self) -> List[str]:
        self._refresh_if_index_changed()
        catalog = self._load_index_catalog()
        if catalog is None:
            self._load_archives()
            subdirs: Set[str] = {artifact.subdir for artifact in self._archived.values()}
        else:
            subdirs = set(str(subdir) for subdir in catalog.get("subdirs", []))
        if self.root.is_dir():
            subdirs.update(
                entry.name
                for entry in self.root.iterdir()
                if entry.is_dir()
                and not entry.is_symlink()
                and not entry.name.startswith(".")
                and entry.name not in EXCLUDED_SUBDIRS
            )
        return sorted(subdirs)

    def logical_subdirs(self) -> List[str]:
        with producer_read_lock(self.archive_root):
            return sorted({LEGACY_SUBDIR_ALIASES.get(subdir, subdir) for subdir in self._logical_subdirs_unlocked()})

    def _live_artifact(self, subdir: str, name: str) -> Optional[Artifact]:
        logical_path = _safe_logical_path(subdir, name)
        live_path = self.root / subdir / name
        if not live_path.is_file() or live_path.is_symlink() or live_path.name.startswith("."):
            return None
        return Artifact(
            logical_path=logical_path,
            subdir=subdir,
            name=name,
            generation=sys.maxsize,
            live_path=live_path,
            size=live_path.stat().st_size,
            mtime_ns=live_path.stat().st_mtime_ns,
            sha256=_sha256_path(live_path),
        )

    def _archived_artifact_is_deleted(self, artifact: Artifact) -> bool:
        for equivalent_path in _equivalent_output_logical_paths(
            artifact.subdir,
            artifact.name,
        ):
            tombstone = self._tombstones.get(equivalent_path)
            if tombstone is not None and tombstone[1] == "delete" and tombstone[0] >= artifact.generation:
                return True
        return False

    def _archived_artifact_unchecked(
        self,
        subdir: str,
        name: str,
    ) -> Optional[Artifact]:
        logical_path = _safe_logical_path(subdir, name)
        if self.family_filter is None:
            archived = self._load_subdir_artifacts(subdir).get(logical_path)
        else:
            self._load_archives()
            archived = self._archived.get(logical_path)
        if archived is None:
            return None
        if self._archived_artifact_is_deleted(archived):
            return None
        return archived

    def _physical_artifact_unchecked(
        self,
        subdir: str,
        name: str,
    ) -> Optional[Artifact]:
        live_artifact = self._live_artifact(subdir, name)
        if live_artifact is not None:
            return live_artifact
        return self._archived_artifact_unchecked(subdir, name)

    def _artifact_unchecked(self, subdir: str, name: str) -> Optional[Artifact]:
        artifact = self._physical_artifact_unchecked(subdir, name)
        if artifact is not None:
            return artifact
        requested_logical_path = _safe_logical_path(subdir, name)
        for legacy_subdir, legacy_name in _legacy_output_candidates(subdir, name):
            artifact = self._physical_artifact_unchecked(
                legacy_subdir,
                legacy_name,
            )
            if artifact is None:
                continue
            tombstone = self._tombstones.get(requested_logical_path)
            if tombstone is not None and tombstone[1] == "delete" and tombstone[0] >= artifact.generation:
                return None
            return replace(
                artifact,
                logical_path=requested_logical_path,
                subdir=subdir,
                name=name,
            )
        return None

    def artifact(self, subdir: str, name: str) -> Optional[Artifact]:
        with producer_read_lock(self.archive_root):
            self._refresh_if_index_changed()
            return self._artifact_unchecked(subdir, name)

    def _physical_file_names_unlocked(
        self,
        subdir: str,
        *,
        family_ids: Optional[Set[str]] = None,
        family_from_name: Optional[Callable[[str], Optional[str]]] = None,
    ) -> List[str]:
        self._refresh_if_index_changed()

        def selected_name(name: str) -> bool:
            if family_ids is None:
                return True
            assert family_from_name is not None
            _, canonical_name = _canonical_output_path(subdir, name)
            return family_from_name(canonical_name) in family_ids

        live_names: Set[str] = set()
        live_dir = self.root / subdir
        if live_dir.is_dir():
            live_names.update(
                entry.name
                for entry in live_dir.iterdir()
                if not entry.name.startswith(".")
                and selected_name(entry.name)
                and entry.is_file()
                and not entry.is_symlink()
            )
        if self.family_filter is None:
            archived_values = self._load_subdir_artifacts(subdir).values()
        else:
            self._load_archives()
            archived_values = (artifact for artifact in self._archived.values() if artifact.subdir == subdir)
        archived_names = {
            artifact.name
            for artifact in archived_values
            if (
                family_ids is None
                or (artifact.family_id in family_ids if artifact.family_id is not None else selected_name(artifact.name))
            )
            and not self._archived_artifact_is_deleted(artifact)
        }
        return sorted(live_names | archived_names)

    def _file_names_unlocked(
        self,
        subdir: str,
        *,
        family_ids: Optional[Set[str]] = None,
        family_from_name: Optional[Callable[[str], Optional[str]]] = None,
    ) -> List[str]:
        self._refresh_if_index_changed()
        physical_subdirs = {subdir}
        physical_subdirs.update(
            legacy_subdir for legacy_subdir, current_subdir in LEGACY_SUBDIR_ALIASES.items() if current_subdir == subdir
        )
        names: Set[str] = set()
        for physical_subdir in physical_subdirs:
            for physical_name in self._physical_file_names_unlocked(
                physical_subdir, family_ids=family_ids, family_from_name=family_from_name
            ):
                canonical_subdir, canonical_name = _canonical_output_path(
                    physical_subdir,
                    physical_name,
                )
                if canonical_subdir == subdir:
                    names.add(canonical_name)
        return sorted(names)

    def file_names(self, subdir: str) -> List[str]:
        with producer_read_lock(self.archive_root):
            return self._file_names_unlocked(subdir)

    def artifacts(self, subdir: str) -> List[Artifact]:
        with producer_read_lock(self.archive_root):
            self._refresh_if_index_changed()
            return [
                artifact
                for name in self._file_names_unlocked(subdir)
                if (artifact := self._artifact_unchecked(subdir, name)) is not None
            ]

    def logical_exists(self, logical_path: str) -> bool:
        subdir, name = logical_path.split("/", 1)
        return self.artifact(subdir, name) is not None

    @contextlib.contextmanager
    def open_binary(
        self,
        subdir: str,
        name: str,
        *,
        _producer_locked: bool = False,
    ) -> Iterator[BinaryIO]:
        lock_context = contextlib.nullcontext() if _producer_locked else producer_read_lock(self.archive_root)
        with lock_context:
            with self._cache_lock:
                self._refresh_if_index_changed()
                artifact = self._artifact_unchecked(subdir, name)
            if artifact is None:
                raise FileNotFoundError(self.root / subdir / name)
            if artifact.live_path is not None:
                with artifact.live_path.open("rb") as handle:
                    yield handle
                return
            assert artifact.zip_path is not None
            assert artifact.member_name is not None
            if artifact.zip_path.is_symlink() or artifact.zip_path.parent.is_symlink():
                raise ArchiveStoreError(f"Symlinked ZIP shards are not supported: {artifact.zip_path}")
            thread = threading.current_thread()
            cached_reader = False
            with self._cache_lock:
                cached = self._zip_reader_cache.get(thread)
                if cached is not None and cached[0] == artifact.zip_path:
                    archive = cached[1]
                    self._zip_reader_cache[thread] = (
                        cached[0],
                        cached[1],
                        cached[2] + 1,
                    )
                    cached_reader = True
                elif cached is None or cached[2] == 0:
                    if cached is not None:
                        cached[1].close()
                    archive = zipfile.ZipFile(artifact.zip_path, "r")
                    self._zip_reader_cache[thread] = (
                        artifact.zip_path,
                        archive,
                        1,
                    )
                    cached_reader = True
                else:
                    archive = zipfile.ZipFile(artifact.zip_path, "r")
            try:
                with archive.open(artifact.member_name, "r") as handle:
                    yield handle
            finally:
                if cached_reader:
                    with self._cache_lock:
                        cached = self._zip_reader_cache.get(thread)
                        if cached is not None and cached[1] is archive:
                            self._zip_reader_cache[thread] = (
                                cached[0],
                                cached[1],
                                max(0, cached[2] - 1),
                            )
                else:
                    archive.close()

    def materialize(
        self,
        subdir: str,
        name: str,
        destination_root: Optional[Path] = None,
        overwrite: bool = False,
        *,
        _producer_locked: bool = False,
    ) -> Path:
        destination_root = self.root if destination_root is None else Path(destination_root).resolve()
        destination = destination_root / subdir / name
        if destination.is_symlink() or destination.parent.is_symlink():
            raise ArchiveStoreError(f"Symlinked materialization destinations are not supported: {destination}")
        if destination.is_file() and not overwrite:
            return destination
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_name(f".{destination.name}.materialize.{os.getpid()}.{uuid.uuid4().hex}")
        try:
            lock_context = contextlib.nullcontext() if _producer_locked else producer_read_lock(self.archive_root)
            with lock_context:
                with self._cache_lock:
                    self._refresh_if_index_changed()
                    artifact = self._artifact_unchecked(subdir, name)
                if artifact is None:
                    raise FileNotFoundError(self.root / subdir / name)
                with (
                    self.open_binary(
                        subdir,
                        name,
                        _producer_locked=True,
                    ) as source,
                    temporary.open("wb") as target,
                ):
                    shutil.copyfileobj(source, target, length=1024 * 1024)
                    target.flush()
                    os.fsync(target.fileno())
            if artifact.mtime_ns is not None:
                os.utime(temporary, ns=(artifact.mtime_ns, artifact.mtime_ns))
            os.replace(temporary, destination)
            _fsync_directory(destination.parent)
        finally:
            if temporary.exists():
                temporary.unlink()
        return destination

    def materialize_family(
        self,
        family_id: str,
        family_from_name: Callable[[str], Optional[str]],
        destination_root: Optional[Path] = None,
        subdirs: Optional[Set[str]] = None,
        receipt_path: Optional[Path] = None,
        run_token: str = "",
    ) -> List[Path]:
        return self.materialize_families(
            {family_id},
            family_from_name,
            destination_root=destination_root,
            subdirs=subdirs,
            receipt_path=receipt_path,
            run_token=run_token,
        )

    def materialize_families(
        self,
        family_ids: Set[str],
        family_from_name: Callable[[str], Optional[str]],
        destination_root: Optional[Path] = None,
        subdirs: Optional[Set[str]] = None,
        receipt_path: Optional[Path] = None,
        run_token: str = "",
    ) -> List[Path]:
        destination_root_resolved = self.root if destination_root is None else Path(destination_root).resolve()
        selected: List[Tuple[Artifact, str]] = []
        with producer_read_lock(self.archive_root):
            self._refresh_if_index_changed()
            logical_subdirs = sorted(
                {LEGACY_SUBDIR_ALIASES.get(subdir, subdir) for subdir in self._logical_subdirs_unlocked()}
            )
            for subdir in logical_subdirs:
                if subdirs is not None and subdir not in subdirs:
                    continue
                # Filter names before stat/hash: another family (or a shared
                # lockfile) can disappear while its producer is running.
                for name in self._file_names_unlocked(
                    subdir, family_ids=family_ids, family_from_name=family_from_name
                ):
                    artifact = self._artifact_unchecked(subdir, name)
                    if artifact is None:
                        continue
                    artifact_family_id = (
                        artifact.family_id if artifact.family_id is not None else family_from_name(name)
                    )
                    if artifact_family_id not in family_ids:
                        continue
                    if artifact.is_live and destination_root_resolved == self.root:
                        continue
                    selected.append((artifact, str(artifact_family_id)))

            if destination_root_resolved == self.root:
                restored_in_place: List[Path] = []
                for artifact, artifact_family_id in selected:
                    if receipt_path is not None and artifact.zip_path is not None:
                        self._append_materialization_receipt(
                            Path(receipt_path),
                            artifact,
                            family_id=artifact_family_id,
                            run_token=run_token,
                        )
                    restored_in_place.append(
                        self.materialize(
                            artifact.subdir,
                            artifact.name,
                            destination_root=destination_root_resolved,
                        )
                    )
                return restored_in_place

            restored: List[Path] = []
            touched_directories: Set[Path] = set()
            by_zip: Dict[Optional[Path], List[Tuple[Artifact, str]]] = {}
            for artifact, artifact_family_id in selected:
                by_zip.setdefault(artifact.zip_path, []).append((artifact, artifact_family_id))

            def write_artifact(
                artifact: Artifact,
                source: BinaryIO,
            ) -> Path:
                destination = destination_root_resolved / artifact.subdir / artifact.name
                if destination.is_symlink() or destination.parent.is_symlink():
                    raise ArchiveStoreError(f"Symlinked materialization destinations are not supported: {destination}")
                if destination.is_file():
                    return destination
                destination.parent.mkdir(parents=True, exist_ok=True)
                temporary = destination.with_name(f".{destination.name}.materialize.{os.getpid()}.{uuid.uuid4().hex}")
                try:
                    with temporary.open("wb") as target:
                        shutil.copyfileobj(source, target, length=1024 * 1024)
                        target.flush()
                        os.fsync(target.fileno())
                    if artifact.mtime_ns is not None:
                        os.utime(
                            temporary,
                            ns=(artifact.mtime_ns, artifact.mtime_ns),
                        )
                    os.replace(temporary, destination)
                    touched_directories.add(destination.parent)
                finally:
                    if temporary.exists():
                        temporary.unlink()
                return destination

            live_entries = by_zip.pop(None, [])
            for artifact, _ in live_entries:
                if artifact.live_path is None:
                    raise ArchiveStoreError(f"Logical artifact has no readable source: {artifact.logical_path}")
                with artifact.live_path.open("rb") as source:
                    restored.append(write_artifact(artifact, source))
            for zip_path, entries in sorted(by_zip.items(), key=lambda item: str(item[0])):
                assert zip_path is not None
                if zip_path.is_symlink() or zip_path.parent.is_symlink():
                    raise ArchiveStoreError(f"Symlinked ZIP shards are not supported: {zip_path}")
                with zipfile.ZipFile(zip_path, "r") as archive:
                    for artifact, _ in sorted(
                        entries,
                        key=lambda item: item[0].logical_path,
                    ):
                        if artifact.member_name is None:
                            raise ArchiveStoreError(f"Archive member is missing a member name: {artifact.logical_path}")
                        with archive.open(artifact.member_name, "r") as source:
                            restored.append(write_artifact(artifact, source))
            for directory in sorted(touched_directories):
                _fsync_directory(directory)
            return sorted(restored)

    def _append_materialization_receipt(
        self,
        receipt_path: Path,
        artifact: Artifact,
        *,
        family_id: str,
        run_token: str,
    ) -> None:
        if artifact.zip_path is None:
            return
        if receipt_path.is_symlink() or receipt_path.parent.is_symlink():
            raise ArchiveStoreError(f"Symlinked materialization receipts are not supported: {receipt_path}")
        receipt_path.parent.mkdir(parents=True, exist_ok=True)
        receipt_existed = receipt_path.exists()
        record = {
            "schema_version": ARCHIVE_SCHEMA_VERSION,
            "root": str(self.root),
            "family_id": family_id,
            "logical_path": artifact.logical_path,
            "generation": artifact.generation,
            "size": int(artifact.size or 0),
            "sha256": artifact.sha256 or "",
            "run_token": run_token,
        }
        with receipt_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        if not receipt_existed:
            _fsync_directory(receipt_path.parent)

    def materialize_all(
        self,
        destination_root: Path,
        overwrite: bool = False,
    ) -> List[Path]:
        destination_root = Path(destination_root).resolve()
        restored: List[Path] = []
        for subdir in self.logical_subdirs():
            for name in self.file_names(subdir):
                restored.append(
                    self.materialize(
                        subdir,
                        name,
                        destination_root=destination_root,
                        overwrite=overwrite,
                    )
                )
        return restored

    def _next_generation(self) -> int:
        counter_path = self.archive_root / GENERATION_FILE
        if counter_path.is_symlink():
            raise ArchiveStoreError(f"Symlinked generation counters are not supported: {counter_path}")
        counter_generation = 0
        if counter_path.is_file():
            try:
                counter_generation = int(counter_path.read_text(encoding="utf-8").strip() or "0")
            except (OSError, ValueError) as exc:
                raise ArchiveStoreError(f"Invalid GeneGalleon generation counter {counter_path}: {exc}") from exc
        else:
            self._load_archives()
        self._index_generation = max(self._index_generation, counter_generation) + 1
        temporary = counter_path.with_name(f".{GENERATION_FILE}.partial.{os.getpid()}.{uuid.uuid4().hex}")
        try:
            with temporary.open("w", encoding="utf-8") as handle:
                handle.write(f"{self._index_generation}\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, counter_path)
            _fsync_directory(self.archive_root)
        finally:
            if temporary.exists():
                temporary.unlink()
        return self._index_generation

    def _index_record(self, artifact: Artifact) -> dict:
        if artifact.zip_path is None or artifact.member_name is None:
            raise ArchiveStoreError(f"Cannot index a live artifact as an archive member: {artifact.logical_path}")
        final_path = _final_archive_path(self.root, artifact.subdir)
        if artifact.zip_path.resolve() == final_path.resolve():
            relative_zip = Path(final_path.name)
            zip_location = "final"
        else:
            try:
                relative_zip = artifact.zip_path.relative_to(self.payload_root)
            except ValueError as exc:
                raise ArchiveStoreError(f"Archive member is outside the payload root: {artifact.zip_path}") from exc
            zip_location = "active"
        return {
            "zip_path": relative_zip.as_posix(),
            "zip_location": zip_location,
            "member_name": artifact.member_name,
            "generation": artifact.generation,
            "size": int(artifact.size or 0),
            "crc": int(artifact.crc or 0),
            "sha256": artifact.sha256 or "",
            "mtime_ns": artifact.mtime_ns,
            "family_id": artifact.family_id,
        }

    def _write_index(
        self,
        artifacts: Optional[Dict[str, Artifact]] = None,
        *,
        recover_pending: bool = False,
    ) -> None:
        with self._index_update(recover_pending=recover_pending):
            self._write_index_uncommitted(artifacts)

    def _write_index_uncommitted(
        self,
        artifacts: Optional[Dict[str, Artifact]] = None,
    ) -> None:
        artifacts = self._archived if artifacts is None else artifacts
        self.archive_root.mkdir(parents=True, exist_ok=True)
        index_dir = self.archive_root / INDEX_DIR_NAME
        if index_dir.is_symlink():
            raise ArchiveStoreError(f"Symlinked archive index directories are not supported: {index_dir}")
        index_dir.mkdir(parents=True, exist_ok=True)
        by_bucket: Dict[str, Dict[str, Artifact]] = {}
        by_subdir: Dict[str, Dict[str, Artifact]] = {}
        for logical_path, artifact in artifacts.items():
            bucket = _family_index_bucket(artifact.family_id)
            by_bucket.setdefault(bucket, {})[logical_path] = artifact
            by_subdir.setdefault(artifact.subdir, {})[logical_path] = artifact
        expected_paths: Set[Path] = set()
        for bucket, bucket_artifacts in sorted(by_bucket.items()):
            index_path = index_dir / f"{bucket}.json"
            expected_paths.add(index_path)
            bucket_generation = max(
                (artifact.generation for artifact in bucket_artifacts.values()),
                default=0,
            )
            payload = {
                "schema_version": INDEX_SCHEMA_VERSION,
                "generation": bucket_generation,
                "artifacts": {
                    logical_path: self._index_record(artifact)
                    for logical_path, artifact in sorted(bucket_artifacts.items())
                },
            }
            encoded = (json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")
            try:
                if index_path.is_file() and index_path.read_bytes() == encoded:
                    continue
            except OSError:
                pass
            temporary = index_path.with_name(f".{index_path.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
            try:
                with temporary.open("wb") as handle:
                    handle.write(encoded)
                    handle.flush()
                    os.fsync(handle.fileno())
                os.replace(temporary, index_path)
            finally:
                if temporary.exists():
                    temporary.unlink()
        for stale_path in sorted(index_dir.glob("*.json")):
            if stale_path not in expected_paths:
                stale_path.unlink()
        subdir_index_dir = self.archive_root / SUBDIR_INDEX_DIR_NAME
        if subdir_index_dir.is_symlink():
            raise ArchiveStoreError(f"Symlinked subdirectory index directories are not supported: {subdir_index_dir}")
        subdir_index_dir.mkdir(parents=True, exist_ok=True)
        expected_subdir_paths: Set[Path] = set()
        for subdir, subdir_artifacts in sorted(by_subdir.items()):
            subdir_path = subdir_index_dir / _subdir_index_name(subdir)
            expected_subdir_paths.add(subdir_path)
            payload = {
                "schema_version": INDEX_SCHEMA_VERSION,
                "subdir": subdir,
                "artifacts": {
                    logical_path: self._index_record(artifact)
                    for logical_path, artifact in sorted(subdir_artifacts.items())
                },
            }
            encoded = (json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")
            try:
                if subdir_path.is_file() and subdir_path.read_bytes() == encoded:
                    continue
            except OSError:
                pass
            temporary = subdir_path.with_name(f".{subdir_path.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
            try:
                with temporary.open("wb") as handle:
                    handle.write(encoded)
                    handle.flush()
                    os.fsync(handle.fileno())
                os.replace(temporary, subdir_path)
            finally:
                if temporary.exists():
                    temporary.unlink()
        for stale_path in sorted(subdir_index_dir.glob("*.json")):
            if stale_path not in expected_subdir_paths:
                stale_path.unlink()
        legacy_index_path = self.archive_root / INDEX_FILE
        if legacy_index_path.is_file() and not legacy_index_path.is_symlink():
            legacy_index_path.unlink()
        self._write_index_catalog(artifacts)
        _fsync_directory(index_dir)
        _fsync_directory(subdir_index_dir)
        _fsync_directory(self.archive_root)
        self._write_index_epoch()
        self._archived = dict(artifacts)
        self._subdir_archived.clear()
        self._index_loaded = True

    def _reset_cache(self) -> None:
        with self._cache_lock:
            for _, archive, _ in self._zip_reader_cache.values():
                archive.close()
            self._zip_reader_cache.clear()
            self._archives_loaded = False
            self._archived.clear()
            self._tombstones.clear()
            self._tombstones_loaded = False
            self._family_states.clear()
            self._states_loaded = False
            self._zip_manifests.clear()
            self._index_loaded = False
            self._subdir_archived.clear()
            self._index_catalog = None
            self._cache_epoch = None

    def _refresh_if_index_changed(self) -> None:
        self._assert_no_pending_index_update()
        epoch = self._read_index_epoch()
        with self._cache_lock:
            if self._cache_epoch is not None and epoch != self._cache_epoch:
                self._reset_cache()

    def _append_tombstone(self, logical_path: str, operation: str) -> None:
        subdir, name = logical_path.split("/", 1)
        _safe_logical_path(subdir, name)
        if operation not in {"delete", "undelete"}:
            raise ValueError(operation)
        self.archive_root.mkdir(parents=True, exist_ok=True)
        generation = self._next_generation()
        record = {
            "schema_version": ARCHIVE_SCHEMA_VERSION,
            "generation": generation,
            "operation": operation,
            "logical_path": logical_path,
        }
        tombstone_path = self.archive_root / TOMBSTONE_FILE
        if tombstone_path.is_symlink():
            raise ArchiveStoreError(f"Symlinked tombstone logs are not supported: {tombstone_path}")
        with tombstone_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        if tombstone_path.stat().st_size > TOMBSTONE_LOG_COMPACT_BYTES:
            self._reset_cache()
            self._load_tombstones()
            temporary = tombstone_path.with_name(f".{tombstone_path.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
            try:
                with temporary.open("w", encoding="utf-8") as handle:
                    for current_logical_path, (
                        current_generation,
                        current_operation,
                    ) in sorted(self._tombstones.items()):
                        compacted_record = {
                            "schema_version": ARCHIVE_SCHEMA_VERSION,
                            "generation": current_generation,
                            "operation": current_operation,
                            "logical_path": current_logical_path,
                        }
                        handle.write(
                            json.dumps(
                                compacted_record,
                                sort_keys=True,
                                separators=(",", ":"),
                            )
                            + "\n"
                        )
                    handle.flush()
                    os.fsync(handle.fileno())
                os.replace(temporary, tombstone_path)
                _fsync_directory(self.archive_root)
            finally:
                if temporary.exists():
                    temporary.unlink()
        self._write_index_epoch()
        self._reset_cache()

    def _store_family_state(self, family_id: str, status: str, run_token: str) -> None:
        if not family_id:
            raise ValueError("family_id is required")
        if status not in {"running", "complete", "failed"}:
            raise ValueError(f"Unsupported family state: {status}")
        state_path = self._state_bucket_path(family_id)
        states = self._read_state_bucket(state_path)
        previous = states.get(family_id)
        generation = 1 if previous is None else previous[0] + 1
        states[family_id] = (generation, status, run_token)
        self._write_state_bucket(state_path, states)
        self._family_states.clear()
        self._states_loaded = False

    def mark_family_state(self, family_id: str, status: str, run_token: str = "") -> bool:
        self._ensure_state_index()
        with (
            family_bucket_lock(
                self.archive_root,
                family_id,
                exclusive=False,
            ),
            state_bucket_lock(self.archive_root, family_id),
        ):
            state_path = self._state_bucket_path(family_id)
            states = self._read_state_bucket(state_path)
            if status != "running" and run_token:
                current = states.get(family_id)
                if current is not None and current[2] != run_token:
                    return False
            self._store_family_state(family_id, status, run_token)
        return True

    def _family_id_for_logical_path(self, logical_path: str) -> Optional[str]:
        subdir, name = logical_path.split("/", 1)
        _safe_logical_path(subdir, name)
        with producer_read_lock(self.archive_root):
            self._refresh_if_index_changed()
            artifact = self._load_subdir_artifacts(subdir).get(logical_path)
            if artifact is not None and artifact.family_id is not None:
                return artifact.family_id
            for legacy_subdir, legacy_name in _legacy_output_candidates(
                subdir,
                name,
            ):
                legacy_logical_path = _safe_logical_path(
                    legacy_subdir,
                    legacy_name,
                )
                artifact = self._load_subdir_artifacts(legacy_subdir).get(legacy_logical_path)
                if artifact is not None and artifact.family_id is not None:
                    return artifact.family_id
        return orthogroup_id_from_name(name)

    def _managed_family_id(
        self,
        logical_path: str,
        family_id: Optional[str],
        operation: str,
    ) -> str:
        inferred_family_id = self._family_id_for_logical_path(logical_path)
        if family_id not in {None, ""} and inferred_family_id is not None and family_id != inferred_family_id:
            raise ArchiveStoreError(
                f"The supplied family ID {family_id!r} differs from the "
                f"inferred family ID {inferred_family_id!r} for {logical_path}"
            )
        if family_id not in {None, ""} and inferred_family_id is None:
            _, name = logical_path.split("/", 1)
            if query_id_from_name(name, [str(family_id)]) != family_id:
                raise ArchiveStoreError(
                    f"The supplied family ID {family_id!r} does not match the logical filename {name!r}"
                )
        resolved_family_id = family_id or inferred_family_id
        if resolved_family_id is None:
            raise ArchiveStoreError(
                f"The gene-family ID could not be inferred for managed {operation}; supply --family-id explicitly"
            )
        return resolved_family_id

    def delete(
        self,
        logical_path: str,
        remove_live: bool = True,
        family_id: Optional[str] = None,
    ) -> None:
        subdir, name = logical_path.split("/", 1)
        _safe_logical_path(subdir, name)
        live_path = self.root / subdir / name
        equivalent_live_paths = [live_path]
        equivalent_live_paths.extend(
            self.root / legacy_subdir / legacy_name
            for legacy_subdir, legacy_name in _legacy_output_candidates(
                subdir,
                name,
            )
        )
        family_id = self._managed_family_id(logical_path, family_id, "deletion")
        family_context = family_bucket_lock(
            self.archive_root,
            family_id,
            exclusive=True,
        )
        with family_context as family_idle:
            if not family_idle:
                raise ArchiveStoreError("Failed to acquire the gene-family lock")
            with producer_quiescence_lock(self.archive_root) as readers_idle:
                if not readers_idle:
                    raise ArchiveStoreError("Failed to acquire the archive maintenance lock")
                with archive_lock(self.archive_root):
                    self._append_tombstone(logical_path, "delete")
                    if remove_live:
                        for equivalent_live_path in equivalent_live_paths:
                            if equivalent_live_path.is_file():
                                equivalent_live_path.unlink()

    def undelete(
        self,
        logical_path: str,
        *,
        _producer_locked: bool = False,
        _family_locked: bool = False,
        family_id: Optional[str] = None,
    ) -> None:
        if _family_locked:
            if family_id in {None, ""}:
                raise ValueError("family_id is required when the family lock is held")
        else:
            family_id = self._managed_family_id(logical_path, family_id, "undelete")
        assert family_id is not None
        family_context = (
            contextlib.nullcontext(True)
            if _family_locked
            else family_bucket_lock(self.archive_root, family_id, exclusive=True)
        )
        lock_context = contextlib.nullcontext(True) if _producer_locked else producer_quiescence_lock(self.archive_root)
        with family_context as family_idle:
            if not family_idle:
                raise ArchiveStoreError("Failed to acquire the gene-family lock")
            with lock_context as readers_idle:
                if not readers_idle:
                    raise ArchiveStoreError("Failed to acquire the archive maintenance lock")
                with archive_lock(self.archive_root):
                    self._append_tombstone(logical_path, "undelete")

    def restore(
        self,
        logical_path: str,
        overwrite: bool = False,
        family_id: Optional[str] = None,
    ) -> Path:
        family_id = self._managed_family_id(logical_path, family_id, "restore")
        family_context = family_bucket_lock(
            self.archive_root,
            family_id,
            exclusive=True,
        )
        with family_context as family_idle:
            if not family_idle:
                raise ArchiveStoreError("Failed to acquire the gene-family lock")
            with producer_quiescence_lock(self.archive_root) as readers_idle:
                if not readers_idle:
                    raise ArchiveStoreError("Failed to acquire the archive maintenance lock")
                self._reset_cache()
                self._load_archives()
                previous_tombstone = self._tombstones.get(logical_path)
                subdir, name = logical_path.split("/", 1)
                _safe_logical_path(subdir, name)
                physical_candidates = [(subdir, name)]
                physical_candidates.extend(_legacy_output_candidates(subdir, name))
                source_available = any(
                    (self.root / candidate_subdir / candidate_name).is_file()
                    or self._load_subdir_artifacts(candidate_subdir).get(
                        _safe_logical_path(candidate_subdir, candidate_name)
                    )
                    is not None
                    for candidate_subdir, candidate_name in physical_candidates
                )
                if not source_available:
                    raise FileNotFoundError(self.root / logical_path)
                self.undelete(
                    logical_path,
                    _producer_locked=True,
                    _family_locked=True,
                    family_id=family_id,
                )
                try:
                    return self.materialize(
                        subdir,
                        name,
                        overwrite=overwrite,
                        _producer_locked=True,
                    )
                except Exception:
                    if previous_tombstone is not None and previous_tombstone[1] == "delete":
                        with archive_lock(self.archive_root):
                            self._append_tombstone(logical_path, "delete")
                    raise

    def verify(
        self,
        progress_callback: Optional[Callable[..., None]] = None,
        *,
        deep: bool = True,
    ) -> List[Path]:
        with producer_read_lock(self.archive_root):
            return self._verify_unlocked(
                progress_callback=progress_callback,
                deep=deep,
            )

    def _verify_unlocked(
        self,
        progress_callback: Optional[Callable[..., None]] = None,
        *,
        deep: bool = True,
    ) -> List[Path]:
        verified: List[Path] = []
        if not self.archive_root.is_dir():
            return verified
        self._load_archives()
        catalog = self._load_index_catalog()
        if catalog is not None and self.family_filter is None:
            subdir_paths: Set[str] = set()
            for subdir in catalog["subdirs"]:
                subdir_artifacts = self._load_subdir_artifacts(str(subdir))
                expected_count = int(catalog["subdir_counts"][subdir])
                if len(subdir_artifacts) != expected_count:
                    raise ArchiveStoreError(
                        f"Subdirectory index count differs from its catalog for {subdir}: "
                        f"{len(subdir_artifacts)} != {expected_count}"
                    )
                for logical_path, artifact in subdir_artifacts.items():
                    if logical_path in subdir_paths:
                        raise ArchiveStoreError(
                            f"Subdirectory indexes contain a duplicate logical path: {logical_path}"
                        )
                    subdir_paths.add(logical_path)
                    indexed_artifact = self._archived.get(logical_path)
                    if indexed_artifact is not None and artifact != indexed_artifact:
                        raise ArchiveStoreError(f"Family and subdirectory archive metadata disagree: {logical_path}")
            if len(subdir_paths) != int(catalog["artifact_count"]):
                raise ArchiveStoreError("Subdirectory indexes do not contain the cataloged artifact count")
            family_paths = set(self._archived)
            if family_paths != subdir_paths:
                missing = sorted(family_paths - subdir_paths)
                extra = sorted(subdir_paths - family_paths)
                raise ArchiveStoreError(
                    f"Family and subdirectory archive indexes disagree: missing={missing[:10]}, extra={extra[:10]}"
                )
        manifest_modes: Dict[Path, str] = {}
        manifest_logical_paths: Set[str] = set()
        manifest_candidate_ranks: Dict[str, Tuple[int, int, str]] = {}
        manifest_candidate_records: Dict[
            str,
            Tuple[Path, Tuple[object, ...]],
        ] = {}
        archive_manifests: List[Tuple[Path, dict]] = []
        for zip_path in _physical_archive_paths(self.root):
            if zip_path.is_symlink() or zip_path.parent.is_symlink():
                raise ArchiveStoreError(f"Symlinked ZIP shards are not supported: {zip_path}")
            archive_manifests.append((zip_path, self._read_manifest(zip_path)))
        total_zip_files = len(archive_manifests)
        total_members = sum(len(manifest["members"]) for _, manifest in archive_manifests)
        total_bytes = (
            sum(
                int(manifest["members"][member_index]["size"])
                for _, manifest in archive_manifests
                for member_index in range(len(manifest["members"]))
            )
            if deep
            else 0
        )
        completed_members = 0
        completed_bytes = 0
        last_progress = time.monotonic()
        if progress_callback is not None:
            progress_callback(
                force=True,
                phase="verifying",
                verify_zip_files_completed=0,
                verify_zip_files_total=total_zip_files,
                verify_members_completed=0,
                verify_members_total=total_members,
                verify_bytes_completed=0,
                verify_bytes_total=total_bytes,
                current_zip="-",
                verify_mode="deep" if deep else "quick",
            )
        for zip_index, (zip_path, manifest) in enumerate(archive_manifests):
            with zipfile.ZipFile(zip_path, "r") as archive:
                for member in manifest["members"]:
                    member_name = str(member["member_name"])
                    if deep:
                        digest = hashlib.sha256()
                        try:
                            with archive.open(member_name, "r") as handle:
                                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                                    digest.update(chunk)
                                    completed_bytes += len(chunk)
                                    now = time.monotonic()
                                    if progress_callback is not None and now - last_progress >= 30.0:
                                        progress_callback(
                                            force=False,
                                            phase="verifying",
                                            verify_zip_files_completed=zip_index,
                                            verify_zip_files_total=total_zip_files,
                                            verify_members_completed=completed_members,
                                            verify_members_total=total_members,
                                            verify_bytes_completed=completed_bytes,
                                            verify_bytes_total=total_bytes,
                                            current_zip=zip_path.name,
                                            verify_mode="deep",
                                        )
                                        last_progress = now
                        except zipfile.BadZipFile as exc:
                            raise ArchiveStoreError(f"CRC verification failed in {zip_path}: {member_name}") from exc
                        if digest.hexdigest() != str(member.get("sha256", "")):
                            raise ArchiveStoreError(f"SHA256 verification failed in {zip_path}: {member_name}")
                    completed_members += 1
                    logical_path = str(member["logical_path"])
                    manifest_logical_paths.add(logical_path)
                    manifest_generation = int(member.get("generation", manifest["generation"]))
                    manifest_family_id = str(member["family_id"]) if member.get("family_id") not in {None, ""} else None
                    member_mtime_ns = member.get("mtime_ns")
                    source_signature = member.get("source_signature")
                    if member_mtime_ns is None and isinstance(source_signature, list) and len(source_signature) >= 4:
                        member_mtime_ns = source_signature[3]
                    manifest_mtime_ns = int(member_mtime_ns) if member_mtime_ns is not None else None
                    observed = (
                        member_name,
                        manifest_generation,
                        int(member["size"]),
                        int(member["crc"]),
                        str(member.get("sha256", "")),
                        manifest_mtime_ns,
                        manifest_family_id,
                    )
                    candidate_rank = (
                        manifest_generation,
                        int(manifest["generation"]),
                        str(zip_path),
                    )
                    if (
                        logical_path not in manifest_candidate_ranks
                        or candidate_rank > manifest_candidate_ranks[logical_path]
                    ):
                        manifest_candidate_ranks[logical_path] = candidate_rank
                        manifest_candidate_records[logical_path] = (
                            zip_path,
                            observed,
                        )
                    artifact = self._archived.get(logical_path)
                    if artifact is None or artifact.zip_path != zip_path:
                        continue
                    expected = (
                        artifact.member_name,
                        artifact.generation,
                        int(artifact.size or 0),
                        int(artifact.crc or 0),
                        artifact.sha256 or "",
                        artifact.mtime_ns,
                        artifact.family_id,
                    )
                    if observed != expected:
                        raise ArchiveStoreError(
                            f"Archive index differs from its manifest for {logical_path} in {zip_path}"
                        )
            manifest_modes[zip_path] = str(manifest["mode"])
            verified.append(zip_path)
            if progress_callback is not None:
                progress_callback(
                    force=True,
                    phase="verifying",
                    verify_zip_files_completed=zip_index + 1,
                    verify_zip_files_total=total_zip_files,
                    verify_members_completed=completed_members,
                    verify_members_total=total_members,
                    verify_bytes_completed=completed_bytes,
                    verify_bytes_total=total_bytes,
                    current_zip=zip_path.name,
                    verify_mode="deep" if deep else "quick",
                )
        for zip_path in {artifact.zip_path for artifact in self._archived.values() if artifact.zip_path is not None}:
            if zip_path not in manifest_modes:
                raise ArchiveStoreError(f"Archive index references a missing ZIP: {zip_path}")
        missing_manifest_members = set(self._archived) - manifest_logical_paths
        if missing_manifest_members:
            logical_path = sorted(missing_manifest_members)[0]
            artifact = self._archived[logical_path]
            raise ArchiveStoreError(
                f"Archive index references a member absent from its manifest: {logical_path} in {artifact.zip_path}"
            )
        if manifest_modes and not self._index_loaded:
            raise ArchiveStoreError(
                "ZIP shards exist but no authoritative archive index was found; "
                "run the repair command before reading this store"
            )
        if (
            self._index_loaded
            and self.family_filter is None
            and (unindexed := manifest_logical_paths - set(self._archived))
        ):
            raise ArchiveStoreError(
                "ZIP manifests contain logical members absent from the archive index: "
                + ", ".join(sorted(unindexed)[:10])
            )
        if self._index_loaded and self.family_filter is None:
            referenced_zip_paths = {
                artifact.zip_path for artifact in self._archived.values() if artifact.zip_path is not None
            }
            referenced_modes = {
                manifest_modes[zip_path] for zip_path in referenced_zip_paths if zip_path in manifest_modes
            }
            if len(referenced_modes) > 1:
                raise ArchiveStoreError(
                    "Referenced ZIP shards use mixed gene-family modes: " + ", ".join(sorted(referenced_modes))
                )
            metadata = _read_store_metadata(self.root)
            if metadata is not None and referenced_modes and referenced_modes != {str(metadata["mode"])}:
                raise ArchiveStoreError(
                    "Archive store metadata mode differs from ZIP manifests: "
                    f"metadata={metadata['mode']}, manifests=" + ",".join(sorted(referenced_modes))
                )
            orphaned = sorted(zip_path for zip_path in manifest_modes if zip_path not in referenced_zip_paths)
            if orphaned:
                raise ArchiveStoreError(
                    "Unreferenced ZIP shards were found: " + ", ".join(str(path) for path in orphaned[:10])
                )
            for logical_path, artifact in self._archived.items():
                candidate = manifest_candidate_records.get(logical_path)
                if candidate is None:
                    continue
                winner_path, winner_record = candidate
                indexed_record = (
                    artifact.member_name,
                    artifact.generation,
                    int(artifact.size or 0),
                    int(artifact.crc or 0),
                    artifact.sha256 or "",
                    artifact.mtime_ns,
                    artifact.family_id,
                )
                if artifact.zip_path != winner_path or indexed_record != winner_record:
                    raise ArchiveStoreError(
                        f"Archive index does not select the latest manifest generation for {logical_path}"
                    )
        return verified


@contextlib.contextmanager
def archive_lock(archive_root: Path, nonblocking: bool = False) -> Iterator[bool]:
    _validate_archive_root(archive_root)
    archive_root.mkdir(parents=True, exist_ok=True)
    with _bucket_lock(archive_root / LOCK_FILE, exclusive=True, nonblocking=nonblocking) as acquired:
        yield acquired


@contextlib.contextmanager
def producer_quiescence_lock(
    archive_root: Path,
    nonblocking: bool = False,
) -> Iterator[bool]:
    """Acquire the exclusive side of the logical-store maintenance lock.

    ZIP readers hold the shared side while opening or materializing members.
    Index-changing maintenance takes the exclusive side so referenced shards
    cannot be replaced or removed while a reader has them open.
    """

    _validate_archive_root(archive_root)
    archive_root.mkdir(parents=True, exist_ok=True)
    with _bucket_lock(archive_root / PRODUCER_LOCK_FILE, exclusive=True, nonblocking=nonblocking) as acquired:
        yield acquired


def default_completion_paths(family_id: str) -> List[str]:
    return [
        f"stat_branch/{family_id}_stat.branch.tsv",
        f"stat_tree/{family_id}_stat.tree.tsv",
        f"tree_plot/{family_id}_tree_plot.pdf",
    ]


def _completion_outputs_present_unchecked(
    store: GeneFamilyOutputStore,
    family_id: str,
) -> bool:
    for logical_path in default_completion_paths(family_id):
        subdir, name = logical_path.split("/", 1)
        artifact = store._artifact_unchecked(subdir, name)
        if artifact is None or (artifact.size or 0) <= 0:
            return False
    return True


def completion_outputs_present(store: GeneFamilyOutputStore, family_id: str) -> bool:
    with producer_read_lock(store.archive_root):
        store._refresh_if_index_changed()
        return _completion_outputs_present_unchecked(store, family_id)


def completed_family_ids(store: GeneFamilyOutputStore, family_ids: Iterable[str]) -> Set[str]:
    completed: Set[str] = set()
    legacy_candidates: List[str] = []
    for family_id in family_ids:
        state = store.family_state(family_id)
        if state == "complete":
            completed.add(family_id)
        elif state is None:
            legacy_candidates.append(family_id)
    store._refresh_if_index_changed()
    for family_id in legacy_candidates:
        if _completion_outputs_present_unchecked(store, family_id):
            # Legacy outputs created before the family-state log was introduced.
            completed.add(family_id)
    return completed


def _verify_zip_crc(
    path: Path,
    progress_callback: Optional[Callable[[int, int, int, int], None]] = None,
) -> None:
    """Read a completed ZIP once so ZipExtFile validates every member CRC."""

    try:
        with zipfile.ZipFile(path, "r") as archive:
            infos = archive.infolist()
            total_members = len(infos)
            total_bytes = sum(int(info.file_size) for info in infos)
            completed_members = 0
            completed_bytes = 0
            last_progress = time.monotonic()
            if progress_callback is not None:
                progress_callback(0, total_members, 0, total_bytes)
            for info in infos:
                with archive.open(info, "r") as handle:
                    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                        completed_bytes += len(chunk)
                        now = time.monotonic()
                        if progress_callback is not None and now - last_progress >= 30.0:
                            progress_callback(
                                completed_members,
                                total_members,
                                completed_bytes,
                                total_bytes,
                            )
                            last_progress = now
                completed_members += 1
            if progress_callback is not None:
                progress_callback(
                    completed_members,
                    total_members,
                    completed_bytes,
                    total_bytes,
                )
    except (OSError, zipfile.BadZipFile) as exc:
        raise ArchiveStoreError(f"CRC verification failed while reading {path}: {exc}") from exc


def _archive_chunk(
    root: Path,
    payload_root: Path,
    subdir: str,
    paths: Sequence[Path],
    mode: str,
    generation: int,
    family_from_name: Callable[[str], Optional[str]],
    *,
    compression: str = "adaptive",
    compression_level: int = 6,
    destination_path: Optional[Path] = None,
    byte_progress_callback: Optional[Callable[[int], None]] = None,
    verification_progress_callback: Optional[Callable[[int, int, int, int], None]] = None,
) -> Tuple[Path, List[Artifact], Dict[Path, ArchivedSourceSignature]]:
    shard_dir = Path(destination_path).resolve().parent if destination_path is not None else payload_root / subdir
    if shard_dir.is_symlink():
        raise ArchiveStoreError(f"Symlinked archive shard directories are not supported: {shard_dir}")
    shard_dir.mkdir(parents=True, exist_ok=True)
    final_path = (
        Path(destination_path).resolve()
        if destination_path is not None
        else shard_dir / f"{subdir}.part-{generation:06d}.zip"
    )
    if final_path.is_symlink():
        raise ArchiveStoreError(f"Symlinked ZIP destinations are not supported: {final_path}")
    if destination_path is not None and final_path.is_file() and not _zip_has_genegalleon_manifest(final_path):
        raise ArchiveStoreError(f"Refusing to replace an unrelated ZIP file: {final_path}")
    partial_path = shard_dir / (f".{final_path.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
    source_signatures = {path: _stat_signature(path) for path in paths}
    signatures: Dict[Path, ArchivedSourceSignature] = {}
    members: List[dict] = []
    last_progress = time.monotonic()
    try:
        with zipfile.ZipFile(
            partial_path,
            mode="w",
            compression=(zipfile.ZIP_STORED if compression == "store" else zipfile.ZIP_DEFLATED),
            compresslevel=(None if compression == "store" else compression_level),
            allowZip64=True,
        ) as archive:
            for path in paths:
                logical_path = _safe_logical_path(subdir, path.name)
                member_name = logical_path
                digest = hashlib.sha256()
                zip_info = zipfile.ZipInfo.from_file(path, arcname=member_name)
                zip_info.compress_type = _compression_for(path, compression)
                with (
                    path.open("rb") as source,
                    archive.open(
                        zip_info,
                        mode="w",
                        force_zip64=True,
                    ) as target,
                ):
                    for chunk in iter(lambda: source.read(1024 * 1024), b""):
                        digest.update(chunk)
                        target.write(chunk)
                        now = time.monotonic()
                        if byte_progress_callback is not None and now - last_progress >= 30.0:
                            try:
                                byte_progress_callback(partial_path.stat().st_size)
                            except FileNotFoundError:
                                pass
                            last_progress = now
                info = archive.getinfo(member_name)
                members.append(
                    {
                        "logical_path": logical_path,
                        "member_name": member_name,
                        "generation": generation,
                        "size": int(info.file_size),
                        "crc": int(info.CRC),
                        "sha256": digest.hexdigest(),
                        "mtime_ns": int(source_signatures[path][3]),
                        "family_id": family_from_name(path.name),
                        "source_signature": list(source_signatures[path]),
                    }
                )
                signatures[path] = (*source_signatures[path], digest.hexdigest())
            manifest = {
                "schema_version": ARCHIVE_SCHEMA_VERSION,
                "generation": generation,
                "mode": mode,
                "subdir": subdir,
                "compression": compression,
                "compression_level": compression_level,
                "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "members": members,
            }
            archive.writestr(
                MANIFEST_MEMBER,
                json.dumps(manifest, sort_keys=True, separators=(",", ":")),
                compress_type=zipfile.ZIP_DEFLATED,
            )
        _verify_zip_crc(partial_path, verification_progress_callback)
        with partial_path.open("rb") as handle:
            os.fsync(handle.fileno())
        os.replace(partial_path, final_path)
        _fsync_directory(shard_dir)
        if byte_progress_callback is not None:
            byte_progress_callback(final_path.stat().st_size)
    finally:
        if partial_path.exists():
            partial_path.unlink()
    artifacts = [
        Artifact(
            logical_path=str(member["logical_path"]),
            subdir=subdir,
            name=str(member["logical_path"]).split("/", 1)[1],
            generation=int(member["generation"]),
            zip_path=final_path,
            member_name=str(member["member_name"]),
            size=int(member["size"]),
            crc=int(member["crc"]),
            sha256=str(member["sha256"]),
            mtime_ns=int(member["mtime_ns"]),
            family_id=(str(member["family_id"]) if member.get("family_id") not in {None, ""} else None),
        )
        for member in members
    ]
    return final_path, artifacts, signatures


def _compact_artifact_chunk(
    payload_root: Path,
    subdir: str,
    artifacts: Sequence[Artifact],
    mode: str,
    archive_generation: int,
    *,
    compression: str = "adaptive",
    compression_level: int = 6,
    destination_path: Optional[Path] = None,
    byte_progress_callback: Optional[Callable[[int], None]] = None,
    verification_progress_callback: Optional[Callable[[int, int, int, int], None]] = None,
) -> Tuple[Path, List[Artifact]]:
    shard_dir = Path(destination_path).resolve().parent if destination_path is not None else payload_root / subdir
    if shard_dir.is_symlink():
        raise ArchiveStoreError(f"Symlinked archive shard directories are not supported: {shard_dir}")
    shard_dir.mkdir(parents=True, exist_ok=True)
    final_path = (
        Path(destination_path).resolve()
        if destination_path is not None
        else shard_dir / f"{subdir}.pack-{archive_generation:06d}.zip"
    )
    if final_path.is_symlink():
        raise ArchiveStoreError(f"Symlinked ZIP destinations are not supported: {final_path}")
    if destination_path is not None and final_path.is_file() and not _zip_has_genegalleon_manifest(final_path):
        raise ArchiveStoreError(f"Refusing to replace an unrelated ZIP file: {final_path}")
    partial_path = shard_dir / (f".{final_path.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
    members: List[dict] = []
    last_progress = time.monotonic()
    try:
        with (
            zipfile.ZipFile(
                partial_path,
                mode="w",
                compression=(zipfile.ZIP_STORED if compression == "store" else zipfile.ZIP_DEFLATED),
                compresslevel=(None if compression == "store" else compression_level),
                allowZip64=True,
            ) as destination,
            contextlib.ExitStack() as source_stack,
        ):
            source_archives: OrderedDict[Path, zipfile.ZipFile] = OrderedDict()
            for artifact in artifacts:
                if artifact.zip_path is None or artifact.member_name is None:
                    raise ArchiveStoreError(f"Cannot compact live artifact: {artifact.logical_path}")
                digest = hashlib.sha256()
                zip_info = zipfile.ZipInfo(artifact.logical_path)
                if artifact.mtime_ns is not None:
                    zip_info.date_time = time.localtime(artifact.mtime_ns / 1_000_000_000)[:6]
                zip_info.compress_type = _compression_for(Path(artifact.name), compression)
                if artifact.zip_path.is_symlink() or artifact.zip_path.parent.is_symlink():
                    raise ArchiveStoreError(f"Symlinked ZIP shards are not supported: {artifact.zip_path}")
                source_archive = source_archives.get(artifact.zip_path)
                if source_archive is None:
                    if len(source_archives) >= MAX_OPEN_COMPACTION_SOURCES:
                        _, evicted_archive = source_archives.popitem(last=False)
                        evicted_archive.close()
                    source_archive = source_stack.enter_context(zipfile.ZipFile(artifact.zip_path, "r"))
                    source_archives[artifact.zip_path] = source_archive
                else:
                    source_archives.move_to_end(artifact.zip_path)
                with source_archive.open(artifact.member_name, "r") as source:
                    with destination.open(zip_info, mode="w", force_zip64=True) as target:
                        for chunk in iter(lambda: source.read(1024 * 1024), b""):
                            digest.update(chunk)
                            target.write(chunk)
                            now = time.monotonic()
                            if byte_progress_callback is not None and now - last_progress >= 30.0:
                                try:
                                    byte_progress_callback(partial_path.stat().st_size)
                                except FileNotFoundError:
                                    pass
                                last_progress = now
                if artifact.sha256 is not None and digest.hexdigest() != artifact.sha256:
                    raise ArchiveStoreError(
                        f"SHA256 verification failed while compacting {artifact.zip_path}: {artifact.member_name}"
                    )
                info = destination.getinfo(artifact.logical_path)
                members.append(
                    {
                        "logical_path": artifact.logical_path,
                        "member_name": artifact.logical_path,
                        "generation": artifact.generation,
                        "size": int(info.file_size),
                        "crc": int(info.CRC),
                        "sha256": digest.hexdigest(),
                        "mtime_ns": artifact.mtime_ns,
                        "family_id": artifact.family_id,
                    }
                )
            manifest = {
                "schema_version": ARCHIVE_SCHEMA_VERSION,
                "generation": archive_generation,
                "mode": mode,
                "subdir": subdir,
                "compression": compression,
                "compression_level": compression_level,
                "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "members": members,
            }
            destination.writestr(
                MANIFEST_MEMBER,
                json.dumps(manifest, sort_keys=True, separators=(",", ":")),
                compress_type=zipfile.ZIP_DEFLATED,
            )
        _verify_zip_crc(partial_path, verification_progress_callback)
        with partial_path.open("rb") as handle:
            os.fsync(handle.fileno())
        os.replace(partial_path, final_path)
        _fsync_directory(shard_dir)
        if byte_progress_callback is not None:
            byte_progress_callback(final_path.stat().st_size)
    finally:
        if partial_path.exists():
            partial_path.unlink()
    compacted = [
        Artifact(
            logical_path=str(member["logical_path"]),
            subdir=subdir,
            name=str(member["logical_path"]).split("/", 1)[1],
            generation=int(member["generation"]),
            zip_path=final_path,
            member_name=str(member["member_name"]),
            size=int(member["size"]),
            crc=int(member["crc"]),
            sha256=str(member["sha256"]),
            mtime_ns=(int(member["mtime_ns"]) if member.get("mtime_ns") is not None else None),
            family_id=(str(member["family_id"]) if member.get("family_id") not in {None, ""} else None),
        )
        for member in members
    ]
    return final_path, compacted


def _remove_archived_sources(
    root: Path,
    signatures: Dict[Path, ArchivedSourceSignature],
) -> int:
    removed = 0
    touched_subdirs: Set[Path] = set()
    for path, signature in signatures.items():
        touched_subdirs.add(path.parent)
        metadata_signature = signature[:5]
        expected_sha256 = signature[5]
        if _signature_matches(path, metadata_signature):
            digest = hashlib.sha256()
            try:
                with path.open("rb") as handle:
                    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                        digest.update(chunk)
            except FileNotFoundError:
                continue
            if digest.hexdigest() != expected_sha256 or not _signature_matches(path, metadata_signature):
                continue
            path.unlink()
            removed += 1
    for live_dir in touched_subdirs:
        _fsync_directory(live_dir)
        try:
            if live_dir.is_dir() and not any(live_dir.iterdir()):
                live_dir.rmdir()
        except OSError:
            pass
    return removed


def _cleanup_partial_archives(archive_root: Path) -> None:
    partial_paths = list(archive_root.glob("*/.*.zip.partial*")) if archive_root.is_dir() else []
    if archive_root.name == ACTIVE_ARCHIVE_DIR_NAME:
        partial_paths.extend(archive_root.parent.glob(".*.zip.partial.*"))
    for partial_path in partial_paths:
        try:
            if partial_path.parent.is_symlink():
                raise ArchiveStoreError(f"Symlinked archive shard directories are not supported: {partial_path.parent}")
            if partial_path.is_file() and not partial_path.is_symlink():
                partial_path.unlink()
        except ArchiveStoreError:
            raise
        except OSError:
            pass


def _cleanup_partial_materializations(root: Path) -> None:
    if not root.is_dir():
        return
    for directory in sorted(root.iterdir()):
        if (
            not directory.is_dir()
            or directory.is_symlink()
            or directory.name.startswith(".")
            or directory.name in EXCLUDED_SUBDIRS
        ):
            continue
        for partial_path in sorted(directory.glob(".*.materialize.*")):
            if partial_path.is_file() and not partial_path.is_symlink():
                partial_path.unlink()
        _fsync_directory(directory)


def _balanced_archive_chunks(
    paths: Sequence[Path],
    max_files_per_shard: int,
    max_bytes_per_shard: int = 0,
) -> List[List[Path]]:
    if not paths:
        return []
    max_files = max(1, int(max_files_per_shard))
    max_bytes = max(0, int(max_bytes_per_shard))
    if max_bytes > 0:
        chunks: List[List[Path]] = []
        current: List[Path] = []
        current_bytes = 0
        for path in paths:
            path_bytes = int(os.stat(path, follow_symlinks=False).st_size)
            if current and (len(current) >= max_files or current_bytes + path_bytes > max_bytes):
                chunks.append(current)
                current = []
                current_bytes = 0
            current.append(path)
            current_bytes += path_bytes
        if current:
            chunks.append(current)
        return chunks
    chunk_count = (len(paths) + max_files - 1) // max_files
    chunk_size = (len(paths) + chunk_count - 1) // chunk_count
    return [list(paths[start : start + chunk_size]) for start in range(0, len(paths), chunk_size)]


def _balanced_artifact_chunks(
    artifacts: Sequence[Artifact],
    max_files_per_shard: int,
    max_bytes_per_shard: int = 0,
) -> List[List[Artifact]]:
    if not artifacts:
        return []
    max_files = max(1, int(max_files_per_shard))
    max_bytes = max(0, int(max_bytes_per_shard))
    chunks: List[List[Artifact]] = []
    current: List[Artifact] = []
    current_bytes = 0
    for artifact in artifacts:
        artifact_bytes = int(artifact.size or 0)
        if current and (len(current) >= max_files or (max_bytes > 0 and current_bytes + artifact_bytes > max_bytes)):
            chunks.append(current)
            current = []
            current_bytes = 0
        current.append(artifact)
        current_bytes += artifact_bytes
    if current:
        chunks.append(current)
    return chunks


def archive_completed_outputs(
    root: Path,
    mode: str,
    family_ids: Sequence[str],
    family_from_name: Callable[[str], Optional[str]],
    min_files: int = 1,
    max_files_per_shard: int = 5000,
    nonblocking: bool = False,
    include_incomplete: bool = False,
    compression: str = "adaptive",
    compression_level: int = 6,
    workers: int = 1,
    catalog_sources: Optional[Sequence[Path]] = None,
    catalog_family_ids: Optional[Sequence[str]] = None,
    preserve_existing_catalog: bool = False,
    progress_callback: Optional[Callable[..., None]] = None,
    direct_final: bool = False,
    direct_final_max_bytes: int = 0,
    max_bytes_per_shard: int = 0,
) -> List[Tuple[Path, int]]:
    root = Path(root).resolve()
    archive_root = _archive_state_root(root)
    payload_root = _archive_payload_root(root)
    created_paths: List[Path] = []
    removed_total = 0
    with (
        lock_available_family_ids(
            archive_root,
            family_ids,
            nonblocking=nonblocking,
        ) as lockable_family_ids,
        producer_quiescence_lock(
            archive_root,
            nonblocking=nonblocking,
        ) as producers_idle,
    ):
        if not producers_idle:
            return []
        with archive_lock(archive_root, nonblocking=nonblocking) as acquired:
            if not acquired:
                return []
            _cleanup_partial_archives(payload_root)
            store = GeneFamilyOutputStore(root)
            _assert_archive_mode(store, mode)
            _write_store_metadata(
                root,
                mode,
                family_ids if catalog_family_ids is None else catalog_family_ids,
                compression=compression,
                compression_level=compression_level,
                catalog_sources=catalog_sources,
                preserve_existing_catalog=preserve_existing_catalog,
            )
            valid_family_ids = set(lockable_family_ids)
            if not valid_family_ids:
                return []
            potential_by_subdir: Dict[str, List[Tuple[Path, str]]] = {}
            for subdir in sorted(
                entry.name
                for entry in root.iterdir()
                if entry.is_dir()
                and not entry.is_symlink()
                and not entry.name.startswith(".")
                and entry.name not in EXCLUDED_SUBDIRS
            ):
                potential = []
                for entry in sorted((root / subdir).iterdir()):
                    if not entry.is_file() or entry.is_symlink() or entry.name.startswith("."):
                        continue
                    family_id = family_from_name(entry.name)
                    if family_id is not None and family_id in valid_family_ids:
                        potential.append((entry, family_id))
                if len(potential) >= min_files:
                    potential_by_subdir[subdir] = potential
            if not potential_by_subdir:
                return []

            if store._load_index_catalog() is None and _physical_archive_paths(root):
                store._load_archives()
                if store._archived and not store._index_loaded:
                    store._write_index()
            live_family_ids = {family_id for potential in potential_by_subdir.values() for _, family_id in potential}
            eligible_family_ids = (
                live_family_ids if include_incomplete else completed_family_ids(store, live_family_ids)
            )
            if not eligible_family_ids:
                return []

            # An offline raw-to-ZIP conversion can write each previously raw
            # subdirectory directly to its user-facing ZIP.  Normal workflow
            # archiving keeps using immutable parts because it must coexist with
            # later live overrides and newly completed families.
            direct_specs: List[Tuple[str, List[Path], int]] = []
            if direct_final:
                for subdir, potential in sorted(potential_by_subdir.items()):
                    candidates = sorted(entry for entry, family_id in potential if family_id in eligible_family_ids)
                    if len(candidates) < min_files:
                        continue
                    if store._load_subdir_artifacts(subdir):
                        continue
                    if (
                        direct_final_max_bytes > 0
                        and sum(int(os.stat(path, follow_symlinks=False).st_size) for path in candidates)
                        > direct_final_max_bytes
                    ):
                        continue
                    direct_specs.append((subdir, candidates, store._next_generation()))

            direct_subdirs = {subdir for subdir, _, _ in direct_specs}
            direct_subdir_file_counts = {subdir: len(candidates) for subdir, candidates, _ in direct_specs}
            direct_completed = 0
            direct_progress_lock = threading.Lock()
            direct_zip_bytes: Dict[str, int] = {}

            def create_direct_final(
                spec: Tuple[str, List[Path], int],
            ) -> Tuple[
                str,
                Path,
                List[Artifact],
                Dict[Path, ArchivedSourceSignature],
            ]:
                subdir, candidates, generation = spec

                def report_direct_bytes(byte_count: int) -> None:
                    if progress_callback is None:
                        return
                    with direct_progress_lock:
                        direct_zip_bytes[subdir] = byte_count
                        progress_callback(
                            force=False,
                            phase="archiving-final",
                            subdir=subdir,
                            current_subdir=subdir,
                            current_zip=f"{subdir}.zip",
                            current_zip_bytes=byte_count,
                            subdir_files=len(candidates),
                            subdir_shards=1,
                            verify_members_completed=0,
                            verify_members_total=0,
                            verify_bytes_completed=0,
                            verify_bytes_total=0,
                            active_subdirs=",".join(sorted(direct_zip_bytes)),
                            active_zip_bytes=sum(direct_zip_bytes.values()),
                            finalized_subdirs=direct_completed,
                            total_subdirs=len(direct_specs),
                        )

                def report_direct_verification(
                    members_completed: int,
                    members_total: int,
                    bytes_completed: int,
                    bytes_total: int,
                ) -> None:
                    if progress_callback is None:
                        return
                    with direct_progress_lock:
                        progress_callback(
                            force=(members_completed == members_total),
                            phase="verifying-final-zip",
                            subdir=subdir,
                            current_subdir=subdir,
                            current_zip=f"{subdir}.zip",
                            current_zip_bytes=direct_zip_bytes.get(subdir, 0),
                            subdir_files=len(candidates),
                            subdir_shards=1,
                            verify_members_completed=members_completed,
                            verify_members_total=members_total,
                            verify_bytes_completed=bytes_completed,
                            verify_bytes_total=bytes_total,
                            active_subdirs=",".join(sorted(direct_zip_bytes)),
                            active_zip_bytes=sum(direct_zip_bytes.values()),
                            finalized_subdirs=direct_completed,
                            total_subdirs=len(direct_specs),
                        )

                report_direct_bytes(0)
                try:
                    zip_path, artifacts, signatures = _archive_chunk(
                        root,
                        payload_root,
                        subdir,
                        candidates,
                        mode,
                        generation,
                        family_from_name,
                        compression=compression,
                        compression_level=compression_level,
                        destination_path=_final_archive_path(root, subdir),
                        byte_progress_callback=report_direct_bytes,
                        verification_progress_callback=(
                            report_direct_verification if progress_callback is not None else None
                        ),
                    )
                    return subdir, zip_path, artifacts, signatures
                finally:
                    if progress_callback is not None:
                        with direct_progress_lock:
                            direct_zip_bytes.pop(subdir, None)
                            remaining_subdirs = sorted(direct_zip_bytes)
                            next_subdir = remaining_subdirs[0] if remaining_subdirs else "-"
                            progress_callback(
                                force=False,
                                phase="archiving-final",
                                subdir=next_subdir,
                                current_subdir=next_subdir,
                                current_zip=(f"{next_subdir}.zip" if remaining_subdirs else "-"),
                                current_zip_bytes=(direct_zip_bytes[next_subdir] if remaining_subdirs else 0),
                                subdir_files=(direct_subdir_file_counts[next_subdir] if remaining_subdirs else 0),
                                subdir_shards=(1 if remaining_subdirs else 0),
                                verify_members_completed=0,
                                verify_members_total=0,
                                verify_bytes_completed=0,
                                verify_bytes_total=0,
                                active_subdirs=(",".join(remaining_subdirs) or "-"),
                                active_zip_bytes=sum(direct_zip_bytes.values()),
                                finalized_subdirs=direct_completed,
                                total_subdirs=len(direct_specs),
                            )

            def commit_direct_final(
                result: Tuple[
                    str,
                    Path,
                    List[Artifact],
                    Dict[Path, ArchivedSourceSignature],
                ],
            ) -> None:
                nonlocal direct_completed, removed_total
                subdir, zip_path, artifacts, signatures = result
                final_subdir_artifacts = {artifact.logical_path: artifact for artifact in artifacts}
                store._merge_index_subdirs({subdir: final_subdir_artifacts})
                removed_total += _remove_archived_sources(root, signatures)
                referenced_paths = {
                    artifact.zip_path.resolve() for artifact in artifacts if artifact.zip_path is not None
                }
                _remove_unreferenced_subdir_archives(
                    root,
                    payload_root,
                    subdir,
                    referenced_paths,
                )
                created_paths.append(zip_path)
                direct_completed += 1
                with direct_progress_lock:
                    direct_zip_bytes.pop(subdir, None)
                    if progress_callback is not None:
                        progress_callback(
                            force=True,
                            phase="committed",
                            subdir=subdir,
                            current_subdir=subdir,
                            current_zip=zip_path.name,
                            subdir_files=len(artifacts),
                            subdir_shards=1,
                            verify_members_completed=0,
                            verify_members_total=0,
                            verify_bytes_completed=0,
                            verify_bytes_total=0,
                            finalized_subdirs=direct_completed,
                            total_subdirs=len(direct_specs),
                            current_zip_bytes=zip_path.stat().st_size,
                            active_subdirs=(",".join(sorted(direct_zip_bytes)) or "-"),
                            active_zip_bytes=sum(direct_zip_bytes.values()),
                            archived_live_files=removed_total,
                            created_zip_shards=len(created_paths),
                        )

            if direct_specs:
                if progress_callback is not None:
                    progress_callback(
                        force=True,
                        phase="archiving-final",
                        current_subdir="-",
                        current_zip="-",
                        current_zip_bytes=0,
                        subdir_files=0,
                        subdir_shards=0,
                        verify_members_completed=0,
                        verify_members_total=0,
                        verify_bytes_completed=0,
                        verify_bytes_total=0,
                        finalized_subdirs=0,
                        total_subdirs=len(direct_specs),
                        archived_live_files=removed_total,
                        created_zip_shards=len(created_paths),
                    )
                direct_workers = min(max(1, workers), len(direct_specs))
                try:
                    if direct_workers == 1:
                        for spec in direct_specs:
                            commit_direct_final(create_direct_final(spec))
                    else:
                        with concurrent.futures.ThreadPoolExecutor(max_workers=direct_workers) as executor:
                            futures = [executor.submit(create_direct_final, spec) for spec in direct_specs]
                            for future in concurrent.futures.as_completed(futures):
                                commit_direct_final(future.result())
                except BaseException:
                    _cleanup_partial_archives(payload_root)
                    raise

            for subdir, potential in sorted(potential_by_subdir.items()):
                if subdir in direct_subdirs:
                    continue
                candidates = sorted(entry for entry, family_id in potential if family_id in eligible_family_ids)
                if len(candidates) < min_files:
                    continue
                subdir_signatures: Dict[
                    Path,
                    Tuple[int, int, int, int, int],
                ] = {}
                subdir_artifacts = dict(store._load_subdir_artifacts(subdir))
                chunk_specs = []
                for chunk in _balanced_archive_chunks(
                    candidates,
                    max_files_per_shard,
                    max_bytes_per_shard,
                ):
                    generation = store._next_generation()
                    chunk_specs.append((generation, chunk))
                if progress_callback is not None:
                    progress_callback(
                        force=True,
                        phase="archiving",
                        subdir=subdir,
                        subdir_files=len(candidates),
                        subdir_shards=len(chunk_specs),
                        archived_live_files=removed_total,
                        created_zip_shards=len(created_paths),
                    )
                chunk_results: Dict[
                    int,
                    Tuple[
                        Path,
                        List[Artifact],
                        Dict[Path, ArchivedSourceSignature],
                    ],
                ] = {}
                futures: List[concurrent.futures.Future] = []
                try:
                    if min(max(1, workers), len(chunk_specs)) == 1:
                        for index, (generation, chunk) in enumerate(chunk_specs):
                            chunk_results[index] = _archive_chunk(
                                root,
                                payload_root,
                                subdir,
                                chunk,
                                mode,
                                generation,
                                family_from_name,
                                compression=compression,
                                compression_level=compression_level,
                            )
                    else:
                        with concurrent.futures.ThreadPoolExecutor(
                            max_workers=min(max(1, workers), len(chunk_specs))
                        ) as executor:
                            future_indexes = {}
                            for index, (generation, chunk) in enumerate(chunk_specs):
                                future = executor.submit(
                                    _archive_chunk,
                                    root,
                                    payload_root,
                                    subdir,
                                    chunk,
                                    mode,
                                    generation,
                                    family_from_name,
                                    compression=compression,
                                    compression_level=compression_level,
                                )
                                futures.append(future)
                                future_indexes[future] = index
                            for future in concurrent.futures.as_completed(future_indexes):
                                chunk_results[future_indexes[future]] = future.result()
                except BaseException:
                    for result in chunk_results.values():
                        try:
                            result[0].unlink()
                        except OSError:
                            pass
                    for future in futures:
                        if not future.done() or future.cancelled():
                            continue
                        try:
                            exc = future.exception()
                        except BaseException:
                            continue
                        if exc is None:
                            try:
                                future.result()[0].unlink()
                            except OSError:
                                pass
                    _cleanup_partial_archives(payload_root)
                    raise
                for index in sorted(chunk_results):
                    zip_path, artifacts, signatures = chunk_results[index]
                    created_paths.append(zip_path)
                    subdir_signatures.update(signatures)
                    for artifact in artifacts:
                        subdir_artifacts[artifact.logical_path] = artifact

                subdir_artifact_values = list(subdir_artifacts.values())
                referenced_shards = {
                    artifact.zip_path for artifact in subdir_artifact_values if artifact.zip_path is not None
                }
                ordered_artifacts = sorted(
                    subdir_artifact_values,
                    key=lambda artifact: artifact.logical_path,
                )
                target_artifact_chunks = _balanced_artifact_chunks(
                    ordered_artifacts,
                    max_files_per_shard,
                    max_bytes_per_shard,
                )
                logical_bytes = sum(int(artifact.size or 0) for artifact in ordered_artifacts)
                retain_named_parts = max_bytes_per_shard > 0 and logical_bytes > max_bytes_per_shard
                minimum_shard_count = max(1, len(target_artifact_chunks))
                if len(referenced_shards) <= minimum_shard_count + MAX_REFERENCED_SHARDS_PER_SUBDIR - 1:
                    final_subdir_artifacts = subdir_artifacts
                else:
                    compacted_artifacts: List[Artifact] = []
                    ordered = ordered_artifacts
                    final_archive = _final_archive_path(root, subdir)
                    keep_final_archive = final_archive.resolve() in {
                        path.resolve() for path in referenced_shards
                    } and not (retain_named_parts)
                    if keep_final_archive:
                        generation = store._next_generation()
                        zip_path, compacted = _compact_artifact_chunk(
                            payload_root,
                            subdir,
                            ordered,
                            mode,
                            generation,
                            compression=compression,
                            compression_level=compression_level,
                            destination_path=final_archive,
                        )
                        created_paths.append(zip_path)
                        compacted_artifacts.extend(compacted)
                    else:
                        for artifact_chunk in target_artifact_chunks:
                            generation = store._next_generation()
                            destination_path = (
                                payload_root / subdir / f"{subdir}.part-{generation:06d}.zip"
                                if retain_named_parts
                                else None
                            )
                            zip_path, compacted = _compact_artifact_chunk(
                                payload_root,
                                subdir,
                                artifact_chunk,
                                mode,
                                generation,
                                compression=compression,
                                compression_level=compression_level,
                                destination_path=destination_path,
                            )
                            created_paths.append(zip_path)
                            compacted_artifacts.extend(compacted)
                    final_subdir_artifacts = {artifact.logical_path: artifact for artifact in compacted_artifacts}

                # Commit and reclaim one logical subdirectory at a time. This bounds
                # peak disk usage to the old and new shards for one subdirectory,
                # instead of duplicating every changed subdirectory until one large
                # index commit finishes.
                store._merge_index_subdirs({subdir: final_subdir_artifacts})
                removed_total += _remove_archived_sources(root, subdir_signatures)
                referenced_paths = {
                    artifact.zip_path.resolve()
                    for artifact in final_subdir_artifacts.values()
                    if artifact.zip_path is not None
                }
                try:
                    _remove_unreferenced_subdir_archives(
                        root,
                        payload_root,
                        subdir,
                        referenced_paths,
                    )
                except OSError as exc:
                    raise ArchiveStoreError(f"Failed to remove obsolete ZIP payload for {subdir}: {exc}") from exc

                # Drop paths compacted away from the result list immediately; the
                # caller should only report shards that are still authoritative.
                created_paths = [path for path in created_paths if path.is_file()]
                if progress_callback is not None:
                    progress_callback(
                        force=True,
                        phase="committed",
                        subdir=subdir,
                        archived_live_files=removed_total,
                        created_zip_shards=len(created_paths),
                    )
            for shard_dir in sorted(payload_root.iterdir()) if payload_root.is_dir() else []:
                try:
                    if shard_dir.is_dir() and not any(shard_dir.iterdir()):
                        shard_dir.rmdir()
                except OSError:
                    pass
            existing_created_paths = [path for path in created_paths if path.is_file()]
            return [(path, removed_total if index == 0 else 0) for index, path in enumerate(existing_created_paths)]


def _visible_live_output_files(root: Path) -> Iterator[Tuple[str, Path]]:
    if not root.is_dir():
        raise FileNotFoundError(f"Gene-family output root was not found: {root}")
    for directory in sorted(root.iterdir()):
        if (
            not directory.is_dir()
            or directory.is_symlink()
            or directory.name.startswith(".")
            or directory.name in EXCLUDED_SUBDIRS
        ):
            continue
        for path in sorted(directory.iterdir()):
            if path.is_file() and not path.is_symlink() and not path.name.startswith("."):
                yield directory.name, path


def _unsupported_output_symlinks(root: Path) -> List[Path]:
    unsupported: List[Path] = []
    if not root.is_dir():
        return unsupported
    for directory in sorted(root.iterdir()):
        if directory.name.startswith(".") or directory.name in EXCLUDED_SUBDIRS:
            continue
        if directory.is_symlink():
            unsupported.append(directory)
            continue
        if not directory.is_dir():
            continue
        unsupported.extend(
            path for path in sorted(directory.iterdir()) if not path.name.startswith(".") and path.is_symlink()
        )
    return unsupported


def scan_live_output_ownership(
    root: Path,
    valid_family_ids: Set[str],
    family_from_name: Callable[[str], Optional[str]],
    *,
    workers: int = 1,
) -> Tuple[List[Path], List[Path]]:
    owned, unmatched, _, _, _, _, _, _ = _scan_live_output_inventory(
        root,
        valid_family_ids,
        family_from_name,
        workers=workers,
    )
    return owned, unmatched


def _scan_live_output_inventory(
    root: Path,
    valid_family_ids: Set[str],
    family_from_name: Callable[[str], Optional[str]],
    *,
    workers: int = 1,
) -> Tuple[
    List[Path],
    List[Path],
    int,
    List[Path],
    List[Path],
    int,
    Dict[str, Tuple[int, int]],
    int,
]:
    root = Path(root).resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Gene-family output root was not found: {root}")
    owned: List[Path] = []
    unmatched: List[Path] = []
    owned_bytes = 0
    unsupported_symlinks: List[Path] = []
    shared_raw: List[Path] = []
    shared_raw_bytes = 0
    owned_subdir_stats: Dict[str, Tuple[int, int]] = {}
    unmatched_bytes = 0
    with os.scandir(root) as root_entries:
        directories = sorted(root_entries, key=lambda entry: entry.name)
    scan_directories: List[Path] = []
    for directory_entry in directories:
        if directory_entry.name.startswith(".") or directory_entry.name in EXCLUDED_SUBDIRS:
            continue
        directory_path = Path(directory_entry.path)
        if directory_entry.is_symlink():
            unsupported_symlinks.append(directory_path)
            continue
        if not directory_entry.is_dir(follow_symlinks=False):
            continue
        scan_directories.append(directory_path)

    def scan_directory(
        directory_path: Path,
    ) -> Tuple[List[Path], List[Path], int, List[Path], List[Path], int, int]:
        directory_owned: List[Path] = []
        directory_unmatched: List[Path] = []
        directory_owned_bytes = 0
        directory_symlinks: List[Path] = []
        directory_shared_raw: List[Path] = []
        directory_shared_raw_bytes = 0
        directory_unmatched_bytes = 0
        with os.scandir(directory_path) as child_entries:
            entries = sorted(child_entries, key=lambda entry: entry.name)
        for entry in entries:
            if entry.name.startswith("."):
                continue
            path = Path(entry.path)
            if entry.is_symlink():
                directory_symlinks.append(path)
                continue
            try:
                stat_result = entry.stat(follow_symlinks=False)
            except FileNotFoundError:
                continue
            if not stat.S_ISREG(stat_result.st_mode):
                continue
            family_id = family_from_name(entry.name)
            if family_id is not None:
                if family_id in valid_family_ids:
                    directory_owned.append(path)
                    directory_owned_bytes += int(stat_result.st_size)
                else:
                    directory_unmatched.append(path)
                    directory_unmatched_bytes += int(stat_result.st_size)
            elif directory_path.name in SHARED_OUTPUT_SUBDIRS:
                directory_shared_raw.append(path)
                directory_shared_raw_bytes += int(stat_result.st_size)
            else:
                directory_unmatched.append(path)
                directory_unmatched_bytes += int(stat_result.st_size)
        return (
            directory_owned,
            directory_unmatched,
            directory_owned_bytes,
            directory_symlinks,
            directory_shared_raw,
            directory_shared_raw_bytes,
            directory_unmatched_bytes,
        )

    worker_count = max(1, min(int(workers), len(scan_directories) or 1))
    if worker_count == 1:
        results = map(scan_directory, scan_directories)
        executor = None
    else:
        executor = concurrent.futures.ThreadPoolExecutor(max_workers=worker_count)
        results = executor.map(scan_directory, scan_directories)
    try:
        for (
            directory_owned,
            directory_unmatched,
            directory_owned_bytes,
            directory_symlinks,
            directory_shared_raw,
            directory_shared_raw_bytes,
            directory_unmatched_bytes,
        ) in results:
            owned.extend(directory_owned)
            unmatched.extend(directory_unmatched)
            owned_bytes += directory_owned_bytes
            unsupported_symlinks.extend(directory_symlinks)
            shared_raw.extend(directory_shared_raw)
            shared_raw_bytes += directory_shared_raw_bytes
            unmatched_bytes += directory_unmatched_bytes
            if directory_owned:
                owned_subdir_stats[directory_owned[0].parent.name] = (
                    directory_owned_bytes,
                    len(directory_owned),
                )
    finally:
        if executor is not None:
            executor.shutdown(wait=True)
    return (
        owned,
        unmatched,
        owned_bytes,
        unsupported_symlinks,
        shared_raw,
        shared_raw_bytes,
        owned_subdir_stats,
        unmatched_bytes,
    )


def _read_storage_conversion_marker(root: Path) -> Optional[dict]:
    marker = _archive_state_root(root) / STORAGE_CONVERSION_FILE
    if marker.is_symlink():
        raise ArchiveStoreError(f"Symlinked storage conversion markers are not supported: {marker}")
    if not marker.is_file():
        return None
    try:
        with marker.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        if not isinstance(payload, dict):
            raise ArchiveStoreError(f"Storage conversion marker is not a JSON object: {marker}")
        if payload.get("target") not in {"zip", "raw"}:
            raise ArchiveStoreError(f"Storage conversion marker has an invalid target: {marker}")
        if payload.get("mode") not in {"query2family", "orthogroup"}:
            raise ArchiveStoreError(f"Storage conversion marker has an invalid mode: {marker}")
        return payload
    except ArchiveStoreError:
        raise
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        raise ArchiveStoreError(f"Failed to read storage conversion marker {marker}: {exc}") from exc


def _write_storage_conversion_marker(root: Path, mode: str, target: str) -> bool:
    archive_root = _archive_state_root(root)
    archive_root.mkdir(parents=True, exist_ok=True)
    existing = _read_storage_conversion_marker(root)
    if existing is not None:
        if existing.get("mode") != mode or existing.get("target") != target:
            raise ArchiveStoreError(
                "A different storage conversion is already pending: "
                f"mode={existing.get('mode')}, target={existing.get('target')}"
            )
        return True
    marker = archive_root / STORAGE_CONVERSION_FILE
    temporary = marker.with_name(f".{marker.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
    payload = {
        "schema_version": ARCHIVE_SCHEMA_VERSION,
        "mode": mode,
        "target": target,
        "pid": os.getpid(),
        "hostname": os.uname().nodename,
        "created_ns": time.time_ns(),
        "updated_ns": time.time_ns(),
        "phase": "starting",
        "completed_subdirs": [],
        "finalized_subdirs": 0,
        "total_subdirs": 0,
        "archived_live_files": 0,
        "created_zip_shards": 0,
        "materialized_files": 0,
    }
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, sort_keys=True, separators=(",", ":"))
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, marker)
        _fsync_directory(archive_root)
    finally:
        if temporary.exists():
            temporary.unlink()
    return False


def _update_storage_conversion_marker(root: Path, **updates: object) -> dict:
    root = Path(root).resolve()
    marker = _archive_state_root(root) / STORAGE_CONVERSION_FILE
    payload = _read_storage_conversion_marker(root)
    if payload is None:
        raise ArchiveStoreError(f"Storage conversion marker disappeared unexpectedly: {marker}")
    payload = dict(payload)
    payload.update(updates)
    payload["pid"] = os.getpid()
    payload["hostname"] = os.uname().nodename
    payload["updated_ns"] = time.time_ns()
    temporary = marker.with_name(f".{marker.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, sort_keys=True, separators=(",", ":"))
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, marker)
        _fsync_directory(marker.parent)
    finally:
        if temporary.exists():
            temporary.unlink()
    return payload


def _finish_storage_conversion(root: Path) -> None:
    marker = _archive_state_root(root) / STORAGE_CONVERSION_FILE
    if marker.is_file() and not marker.is_symlink():
        marker.unlink()
        _fsync_directory(marker.parent)


@contextlib.contextmanager
def storage_conversion_session(
    root: Path,
    mode: str,
    target: str,
    *,
    require_resume: bool = False,
) -> Iterator[bool]:
    archive_root = _archive_state_root(root)
    archive_root.mkdir(parents=True, exist_ok=True)
    with _bucket_lock(
        archive_root / STORAGE_CONVERSION_LOCK_FILE,
        exclusive=True,
        nonblocking=True,
    ) as acquired:
        if not acquired:
            raise ArchiveStoreError(f"Another storage conversion is active below {root}")
        existing = _read_storage_conversion_marker(root)
        if require_resume and existing is None:
            raise ArchiveStoreError("--resume was requested, but no interrupted storage conversion exists")
        resumed = _write_storage_conversion_marker(root, mode, target)
        yield resumed


def _assert_archive_mode(store: GeneFamilyOutputStore, mode: str) -> None:
    resolve_archive_mode(store.root, explicit_mode=mode, required=True)


def _zip_to_raw_requirements(
    root: Path,
    artifacts: Sequence[Artifact],
    *,
    block_size: Optional[int] = None,
) -> Tuple[int, int]:
    root = Path(root).resolve()
    if block_size is None:
        filesystem_stats = os.statvfs(root)
        block_size = int(filesystem_stats.f_frsize or filesystem_stats.f_bsize or 4096)
    block_size = max(1, int(block_size))
    artifact_subdirs = {root / artifact.subdir for artifact in artifacts}
    missing_subdirs = {subdir for subdir in artifact_subdirs if not subdir.is_dir()}
    allocated_file_bytes = sum(
        ((int(artifact.size or 0) + block_size - 1) // block_size) * block_size for artifact in artifacts
    )
    required_bytes = allocated_file_bytes + len(missing_subdirs) * block_size
    required_inodes = len(artifacts) + len(missing_subdirs) + (1 if artifacts else 0)
    return required_bytes, required_inodes


def _raw_to_zip_requirements(
    subdir_stats: Dict[str, Tuple[int, int]],
    *,
    workers: int = 1,
) -> Tuple[int, str, int, int]:
    """Return a conservative temporary-byte estimate for raw-to-ZIP conversion.

    Completed ZIPs replace their raw inputs, but up to ``workers`` ZIPs can be
    under construction while all corresponding sources still exist.  Deflate
    can expand incompressible input slightly, so include one percent plus
    per-member and per-archive overhead rather than assuming compression.
    """

    total_bytes = sum(logical_bytes for logical_bytes, _ in subdir_stats.values())
    total_files = sum(file_count for _, file_count in subdir_stats.values())

    estimates = []
    for subdir, (logical_bytes, file_count) in subdir_stats.items():
        upper_bound = logical_bytes + max(64 * 1024, (logical_bytes + 99) // 100) + file_count * 512 + 1024 * 1024
        estimates.append((upper_bound, subdir))
    estimates.sort(reverse=True)
    concurrent = estimates[: max(1, min(int(workers), len(estimates) or 1))]
    # Already committed ZIPs generally replace at least as many raw bytes.
    # Keep a small aggregate allowance for ZIP/filesystem overhead while later
    # subdirectories are still being written.
    aggregate_overhead = max(64 * 1024, (total_bytes + 99) // 100) + total_files * 512 + len(subdir_stats) * 1024 * 1024
    peak_bytes = sum(value for value, _ in concurrent) + aggregate_overhead
    peak_subdirs = ",".join(name for _, name in concurrent) or "-"
    largest_bytes = estimates[0][0] if estimates else 0
    return peak_bytes, peak_subdirs, largest_bytes, aggregate_overhead


def storage_conversion_summary(
    root: Path,
    valid_family_ids: Set[str],
    family_from_name: Callable[[str], Optional[str]],
    *,
    workers: int = 1,
    large_zip_warning_bytes: int = DEFAULT_LARGE_ZIP_WARNING_BYTES,
    max_final_zip_bytes: int = 0,
) -> Tuple[dict, List[Path]]:
    root = Path(root).resolve()
    (
        owned,
        unmatched,
        owned_bytes,
        unsupported_symlinks,
        shared_raw,
        shared_raw_bytes,
        owned_subdir_stats,
        unmatched_bytes,
    ) = _scan_live_output_inventory(
        root,
        valid_family_ids,
        family_from_name,
        workers=workers,
    )
    shard_paths = _physical_archive_paths(root)
    archived_only_count = 0
    archived_only_bytes = 0
    archived_only_artifacts: List[Artifact] = []
    if shard_paths:
        store = GeneFamilyOutputStore(root)
        live_logical_paths = {
            _safe_logical_path(path.parent.name, path.name) for path in (*owned, *unmatched, *shared_raw)
        }
        with producer_read_lock(store.archive_root):
            store._refresh_if_index_changed()
            store._load_archives()
            visible_archived = {
                logical_path: artifact
                for logical_path, artifact in store._archived.items()
                if not store._archived_artifact_is_deleted(artifact)
            }
            logical_count = len(live_logical_paths | set(visible_archived))
            archived_only_artifacts = [
                artifact
                for logical_path, artifact in visible_archived.items()
                if logical_path not in live_logical_paths
            ]
            archived_only_count = len(archived_only_artifacts)
            archived_only_bytes = sum(int(artifact.size or 0) for artifact in archived_only_artifacts)
            indexed_artifact_count = len(store._archived)
            archived_live_overrides = len(live_logical_paths & set(visible_archived))
    else:
        logical_count = len(owned) + len(unmatched) + len(shared_raw)
        indexed_artifact_count = 0
        archived_live_overrides = 0
    marker = _read_storage_conversion_marker(root)
    raw_allocated_bytes, raw_peak_new_files = _zip_to_raw_requirements(
        root,
        archived_only_artifacts,
    )
    (
        raw_zip_peak_new_bytes,
        raw_zip_peak_subdirs,
        raw_zip_largest_subdir_bytes,
        raw_zip_max_net_growth_bytes,
    ) = _raw_to_zip_requirements(owned_subdir_stats, workers=workers)
    raw_subdir_bytes: Dict[str, int] = {
        subdir: logical_bytes for subdir, (logical_bytes, _) in owned_subdir_stats.items()
    }
    for artifact in archived_only_artifacts:
        raw_subdir_bytes[artifact.subdir] = raw_subdir_bytes.get(artifact.subdir, 0) + int(artifact.size or 0)
    raw_subdir_files: Dict[str, int] = {subdir: file_count for subdir, (_, file_count) in owned_subdir_stats.items()}
    for artifact in archived_only_artifacts:
        raw_subdir_files[artifact.subdir] = raw_subdir_files.get(artifact.subdir, 0) + 1
    combined_upper_bounds = {
        subdir: (
            logical_bytes
            + max(64 * 1024, (logical_bytes + 99) // 100)
            + raw_subdir_files.get(subdir, 0) * 512
            + 1024 * 1024
        )
        for subdir, logical_bytes in raw_subdir_bytes.items()
    }
    if combined_upper_bounds:
        largest_combined_subdir, largest_combined_bytes = max(
            combined_upper_bounds.items(),
            key=lambda item: (item[1], item[0]),
        )
        if largest_combined_bytes > raw_zip_peak_new_bytes:
            raw_zip_peak_new_bytes = largest_combined_bytes
            raw_zip_peak_subdirs = largest_combined_subdir
        raw_zip_largest_subdir_bytes = max(
            raw_zip_largest_subdir_bytes,
            largest_combined_bytes,
        )
    else:
        raw_zip_peak_new_bytes = 0
        raw_zip_peak_subdirs = "-"
        raw_zip_largest_subdir_bytes = 0
    zip_bytes = sum(path.stat().st_size for path in shard_paths)
    combined_upper_total = sum(combined_upper_bounds.values())
    raw_zip_max_net_growth_bytes = max(
        0,
        combined_upper_total - owned_bytes - zip_bytes,
    )
    large_zip_subdirs = sorted(
        subdir
        for subdir, logical_bytes in raw_subdir_bytes.items()
        if large_zip_warning_bytes > 0 and logical_bytes > large_zip_warning_bytes
    )
    split_zip_subdirs = sorted(
        subdir
        for subdir, logical_bytes in raw_subdir_bytes.items()
        if max_final_zip_bytes > 0 and logical_bytes > max_final_zip_bytes
    )
    if marker is not None:
        storage = "conversion-in-progress"
    elif shard_paths and (owned or archived_live_overrides):
        storage = "mixed"
    elif shard_paths:
        storage = "zip"
    else:
        storage = "raw"
    conflicting_final_zips = _conflicting_final_zip_paths(
        root,
        {path.parent.name for path in owned},
    )
    shard_resolved_paths = {path.resolve() for path in shard_paths}
    metadata_files = (
        sum(
            1
            for path in _archive_state_root(root).rglob("*")
            if path.is_file() and path.resolve() not in shard_resolved_paths
        )
        if _archive_state_root(root).is_dir()
        else 0
    )
    metadata_files += sum(
        int(path.is_file())
        for path in (
            root / "README_GENE_FAMILY_OUTPUTS.txt",
            root / "ARCHIVE_STATUS.tsv",
        )
    )
    metadata_bytes = (
        sum(
            path.stat().st_size
            for path in _archive_state_root(root).rglob("*")
            if path.is_file() and path.resolve() not in shard_resolved_paths
        )
        if _archive_state_root(root).is_dir()
        else 0
    )
    metadata_bytes += sum(
        path.stat().st_size
        for path in (
            root / "README_GENE_FAMILY_OUTPUTS.txt",
            root / "ARCHIVE_STATUS.tsv",
        )
        if path.is_file()
    )
    summary = {
        "storage": storage,
        "mode": (str(metadata["mode"]) if (metadata := _read_store_metadata(root)) is not None else "unknown"),
        "logical_files": logical_count,
        "managed_logical_files": len(owned) + archived_only_count,
        "managed_logical_bytes": owned_bytes + archived_only_bytes,
        "owned_live_files": len(owned),
        "owned_live_bytes": owned_bytes,
        "unmatched_live_files": len(unmatched),
        "unmatched_live_bytes": unmatched_bytes,
        "shared_raw_files": len(shared_raw),
        "shared_raw_bytes": shared_raw_bytes,
        "indexed_artifacts": indexed_artifact_count,
        "archived_live_overrides": archived_live_overrides,
        "zip_files": len(shard_paths),
        "zip_bytes": zip_bytes,
        "raw_zip_peak_new_bytes": raw_zip_peak_new_bytes,
        "raw_zip_peak_subdirs": raw_zip_peak_subdirs,
        "raw_zip_largest_subdir_bytes": raw_zip_largest_subdir_bytes,
        "raw_zip_max_net_growth_bytes": raw_zip_max_net_growth_bytes,
        "large_zip_warning_bytes": int(large_zip_warning_bytes),
        "large_zip_subdir_count": len(large_zip_subdirs),
        "large_zip_subdirs": ",".join(large_zip_subdirs) or "-",
        "max_final_zip_bytes": int(max_final_zip_bytes),
        "split_zip_subdir_count": len(split_zip_subdirs),
        "split_zip_subdirs": ",".join(split_zip_subdirs) or "-",
        "raw_materialize_files": archived_only_count,
        "raw_materialize_bytes": archived_only_bytes,
        "raw_materialize_allocated_bytes": raw_allocated_bytes,
        "raw_peak_new_files": raw_peak_new_files,
        "unsupported_symlinks": len(unsupported_symlinks),
        "conflicting_final_zip_files": len(conflicting_final_zips),
        "metadata_files": metadata_files,
        "metadata_bytes": metadata_bytes,
        "physical_store_files": (
            len(owned)
            + len(unmatched)
            + len(shared_raw)
            + len(unsupported_symlinks)
            + len(shard_paths)
            + metadata_files
        ),
        "physical_managed_files": len(owned) + len(shard_paths) + metadata_files,
        "physical_managed_bytes": owned_bytes + zip_bytes + metadata_bytes,
        "physical_store_bytes": (owned_bytes + shared_raw_bytes + zip_bytes + metadata_bytes + unmatched_bytes),
    }
    if marker is not None:
        summary["conversion_target"] = str(marker["target"])
        summary["conversion_mode"] = str(marker["mode"])
        summary["conversion_phase"] = str(marker.get("phase", "unknown"))
        summary["conversion_updated_ns"] = int(marker.get("updated_ns", 0))
    return summary, unmatched


def storage_conversion_status(root: Path) -> dict:
    root = Path(root).resolve()
    archive_root = _archive_state_root(root)
    payload_root = _archive_payload_root(root)
    marker = _read_storage_conversion_marker(root)
    metadata = _read_store_metadata(root)
    partial_paths = sorted(payload_root.glob("*/.*.zip.partial*"))
    if payload_root.name == ACTIVE_ARCHIVE_DIR_NAME:
        partial_paths.extend(sorted(root.glob(".*.zip.partial.*")))
    zip_paths = _physical_archive_paths(root)
    lock_paths = sorted(archive_root.glob("**/*.lock")) if archive_root.is_dir() else []
    status = {
        "conversion": "in-progress" if marker is not None else "idle",
        "mode": (
            str(marker["mode"])
            if marker is not None
            else (str(metadata["mode"]) if metadata is not None else "unknown")
        ),
        "target": str(marker["target"]) if marker is not None else "none",
        "phase": str(marker.get("phase", "unknown")) if marker is not None else "none",
        "zip_files": len(zip_paths),
        "partial_zip_files": len(partial_paths),
        "archive_lock_files": len(lock_paths),
    }
    if metadata is not None:
        status.update(
            {
                "catalog_family_count": int(metadata.get("catalog_family_count", 0)),
                "catalog_family_ids_sha256": str(metadata.get("catalog_family_ids_sha256", "")),
                "compression": str(metadata.get("compression", "unknown")),
                "compression_level": metadata.get("compression_level", "unknown"),
                "family_lock_stripes": int(metadata.get("family_lock_stripes", FAMILY_LOCK_STRIPES)),
            }
        )
    if marker is not None:
        for key in (
            "created_ns",
            "updated_ns",
            "resumed",
            "owned_live_files",
            "owned_live_bytes",
            "archived_live_files",
            "created_zip_shards",
            "materialized_files",
            "materialized_bytes",
            "current_path",
            "current_subdir",
            "current_zip",
            "current_zip_bytes",
            "active_subdirs",
            "active_zip_bytes",
            "completed_subdirs",
            "finalized_subdirs",
            "total_subdirs",
            "repair_zip_files_completed",
            "repair_zip_files_total",
            "repair_zip_bytes_processed",
            "repair_zip_bytes_total",
            "repair_members_verified",
            "repair_bytes_verified",
            "verify_members_completed",
            "verify_members_total",
            "verify_bytes_completed",
            "verify_bytes_total",
            "verify_zip_files_completed",
            "verify_zip_files_total",
            "workers",
        ):
            if key in marker:
                value = marker[key]
                status[key] = ",".join(str(item) for item in value) if isinstance(value, list) else value
    return status


def optimize_archive_metadata(root: Path) -> dict:
    """Prune obsolete 256-way lock files after all project jobs are stopped."""

    root = Path(root).resolve()
    archive_root = _archive_state_root(root)
    if not archive_root.is_dir():
        return {"removed_legacy_lock_files": 0, "remaining_lock_files": 0}
    current_names = {f"{value:02x}.lock" for value in range(FAMILY_LOCK_STRIPES)}
    all_stripe_paths = []
    candidate_paths = []
    for directory_name in (FAMILY_LOCK_DIR_NAME, FAMILY_STATE_LOCK_DIR_NAME):
        directory = archive_root / directory_name
        if directory.is_symlink():
            raise ArchiveStoreError(f"Symlinked lock directories are not supported: {directory}")
        directory_paths = [
            path for path in sorted(directory.glob("*.lock")) if path.is_file() and not path.is_symlink()
        ]
        all_stripe_paths.extend(directory_paths)
        candidate_paths.extend(path for path in directory_paths if path.name not in current_names)
    removed = 0
    metadata_created = False
    with _bucket_lock(
        archive_root / STORAGE_CONVERSION_LOCK_FILE,
        exclusive=True,
        nonblocking=True,
    ) as conversion_idle:
        if not conversion_idle:
            raise ArchiveStoreError("Cannot optimize metadata during a storage conversion")
        if _read_storage_conversion_marker(root) is not None:
            raise ArchiveStoreError("Cannot optimize metadata while a conversion needs resuming")
        with producer_quiescence_lock(archive_root, nonblocking=True) as readers_idle:
            if not readers_idle:
                raise ArchiveStoreError("Cannot optimize metadata while archive readers are active")
            with archive_lock(archive_root, nonblocking=True) as archive_idle:
                if not archive_idle:
                    raise ArchiveStoreError("Cannot optimize metadata while archive maintenance is active")
                with contextlib.ExitStack() as stack:
                    acquired_paths: List[Path] = []
                    for path in all_stripe_paths:
                        acquired = stack.enter_context(_bucket_lock(path, exclusive=True, nonblocking=True))
                        if not acquired:
                            raise ArchiveStoreError(f"Cannot remove a lock that is still in use: {path}")
                        acquired_paths.append(path)
                    removable_paths = set(candidate_paths)
                    for path in acquired_paths:
                        if path not in removable_paths:
                            continue
                        try:
                            path.unlink()
                            removed += 1
                        except FileNotFoundError:
                            pass
                for directory_name in (FAMILY_LOCK_DIR_NAME, FAMILY_STATE_LOCK_DIR_NAME):
                    _fsync_directory(archive_root / directory_name)
                if _read_store_metadata(root) is None and _physical_archive_paths(root):
                    mode = resolve_archive_mode(root, required=True)
                    assert mode is not None
                    store = GeneFamilyOutputStore(root)
                    store._load_archives()
                    _write_store_metadata(
                        root,
                        mode,
                        (artifact.family_id for artifact in store._archived.values() if artifact.family_id is not None),
                    )
                    metadata_created = True
    remaining = sum(1 for _ in archive_root.glob("**/*.lock"))
    return {
        "removed_legacy_lock_files": removed,
        "remaining_lock_files": remaining,
        "family_lock_stripes": FAMILY_LOCK_STRIPES,
        "store_metadata_created": metadata_created,
    }


def convert_storage_to_zip(
    root: Path,
    mode: str,
    family_ids: Sequence[str],
    family_from_name: Callable[[str], Optional[str]],
    *,
    max_files_per_shard: int = 5000,
    strict_unmatched: bool = False,
    compression: str = "adaptive",
    compression_level: int = 6,
    workers: int = 1,
    catalog_sources: Optional[Sequence[Path]] = None,
    require_resume: bool = False,
    available_bytes: Optional[int] = None,
    large_zip_warning_bytes: int = DEFAULT_LARGE_ZIP_WARNING_BYTES,
    max_final_zip_bytes: int = 0,
    progress_callback: Optional[Callable[..., None]] = None,
) -> dict:
    root = Path(root).resolve()
    valid_family_ids = set(family_ids)
    summary, unmatched = storage_conversion_summary(
        root,
        valid_family_ids,
        family_from_name,
        workers=workers,
        large_zip_warning_bytes=large_zip_warning_bytes,
        max_final_zip_bytes=max_final_zip_bytes,
    )
    conflicting_final_zips = _catalog_conflicting_final_zip_paths(
        root,
        valid_family_ids,
        family_from_name,
    )
    if conflicting_final_zips:
        raise ArchiveStoreError(
            "Unrelated ZIP files would collide with GeneGalleon final archive "
            "destinations; move or rename them before converting: "
            + ", ".join(str(path) for path in conflicting_final_zips[:10])
        )
    if strict_unmatched and unmatched:
        raise ArchiveStoreError(
            "Unmatched live output files were found; rerun without --strict-unmatched "
            "to leave them in place: " + ", ".join(str(path) for path in unmatched[:10])
        )
    if summary["unsupported_symlinks"]:
        unsupported_symlinks = _unsupported_output_symlinks(root)
        raise ArchiveStoreError(
            "Symlinked output paths cannot be converted safely: "
            + ", ".join(str(path) for path in unsupported_symlinks[:10])
        )
    filesystem_free_bytes = int(shutil.disk_usage(root).free)
    effective_available_bytes = filesystem_free_bytes
    if available_bytes is not None:
        if int(available_bytes) < 0:
            raise ArchiveStoreError("available_bytes must be non-negative")
        effective_available_bytes = min(
            effective_available_bytes,
            int(available_bytes),
        )
    required_peak_bytes = int(summary["raw_zip_peak_new_bytes"])
    summary["filesystem_free_bytes"] = filesystem_free_bytes
    summary["effective_available_bytes"] = effective_available_bytes
    if required_peak_bytes > effective_available_bytes:
        raise ArchiveStoreError(
            "Insufficient temporary space for raw-to-ZIP conversion: "
            f"estimated_peak_new_bytes={required_peak_bytes}, "
            f"effective_available_bytes={effective_available_bytes}, "
            f"peak_subdirs={summary['raw_zip_peak_subdirs']}. "
            "Filesystem free space may exceed a user/project quota; pass the "
            "quota remaining as --available-bytes for a quota-aware preflight."
        )
    archived_removed = 0
    created_paths: Set[Path] = set()
    completed_subdirs: Set[str] = set()
    with storage_conversion_session(
        root,
        mode,
        "zip",
        require_resume=require_resume,
    ) as resumed:
        if resumed:
            _cleanup_partial_archives(_archive_payload_root(root))

            def repair_progress(*, force: bool = False, **fields: object) -> None:
                _update_storage_conversion_marker(root, resumed=True, **fields)
                if progress_callback is not None:
                    progress_callback(force=force, resumed=True, **fields)

            repair_progress(
                force=True,
                phase="repairing-index",
                repair_zip_files_completed=0,
                repair_zip_files_total=0,
                current_zip="-",
            )
            repair_archive_index(
                root,
                remove_orphans=True,
                progress_callback=repair_progress,
            )
            repair_progress(force=True, phase="preflight", current_zip="-")
        preflight_store = GeneFamilyOutputStore(root)
        preflight_store.verify()
        _assert_archive_mode(preflight_store, mode)
        _write_store_metadata(
            root,
            mode,
            family_ids,
            compression=compression,
            compression_level=compression_level,
            catalog_sources=catalog_sources,
        )
        existing_marker = _read_storage_conversion_marker(root) or {}
        completed_subdirs.update(str(value) for value in existing_marker.get("completed_subdirs", []))
        _update_storage_conversion_marker(
            root,
            phase="inventory-complete",
            resumed=resumed,
            owned_live_files=int(summary["owned_live_files"]),
            owned_live_bytes=int(summary["owned_live_bytes"]),
            compression=compression,
            compression_level=int(compression_level),
            workers=int(workers),
            estimated_peak_new_bytes=required_peak_bytes,
            effective_available_bytes=effective_available_bytes,
        )
        if progress_callback is not None:
            progress_callback(
                force=True,
                phase="inventory-complete",
                resumed=resumed,
                owned_live_files=int(summary["owned_live_files"]),
                owned_live_bytes=int(summary["owned_live_bytes"]),
            )
        for _ in range(3):
            removed_offset = archived_removed
            shard_offset = len(created_paths)

            def conversion_progress(
                *,
                force: bool = False,
                _removed_offset: int = removed_offset,
                _shard_offset: int = shard_offset,
                **fields: object,
            ) -> None:
                nonlocal completed_subdirs
                progress_fields = dict(fields)
                if "archived_live_files" in progress_fields:
                    progress_fields["archived_live_files"] = _removed_offset + int(
                        progress_fields["archived_live_files"]
                    )
                if "created_zip_shards" in progress_fields:
                    progress_fields["created_zip_shards"] = _shard_offset + int(progress_fields["created_zip_shards"])
                if progress_fields.get("phase") == "committed" and progress_fields.get("subdir"):
                    completed_subdirs.add(str(progress_fields["subdir"]))
                _update_storage_conversion_marker(
                    root,
                    completed_subdirs=sorted(completed_subdirs),
                    **progress_fields,
                )
                if progress_callback is not None:
                    progress_callback(force=force, **progress_fields)

            results = archive_completed_outputs(
                root=root,
                mode=mode,
                family_ids=family_ids,
                family_from_name=family_from_name,
                min_files=1,
                max_files_per_shard=max(1, max_files_per_shard),
                nonblocking=False,
                include_incomplete=True,
                compression=compression,
                compression_level=compression_level,
                workers=workers,
                catalog_sources=catalog_sources,
                progress_callback=conversion_progress,
                direct_final=True,
                direct_final_max_bytes=max_final_zip_bytes,
                max_bytes_per_shard=max_final_zip_bytes,
            )
            for zip_path, removed in results:
                created_paths.add(zip_path)
                archived_removed += removed
            remaining, _ = scan_live_output_ownership(
                root,
                valid_family_ids,
                family_from_name,
                workers=workers,
            )
            if not remaining:
                break
        else:
            raise ArchiveStoreError(
                "Owned live files kept changing during raw-to-ZIP conversion; "
                "stop older gene-family jobs and rerun the same command"
            )
        _update_storage_conversion_marker(
            root,
            phase="finalizing",
            finalized_subdirs=0,
            total_subdirs=0,
            current_subdir="-",
            current_zip="-",
            current_zip_bytes=0,
        )
        if progress_callback is not None:
            progress_callback(
                force=True,
                phase="finalizing",
                finalized_subdirs=0,
                total_subdirs=0,
                current_subdir="-",
                subdir="-",
                current_zip="-",
                current_zip_bytes=0,
            )
        finalized_paths = finalize_archives(
            root,
            mode,
            family_ids,
            nonblocking=False,
            compression=compression,
            compression_level=compression_level,
            max_final_zip_bytes=max_final_zip_bytes,
            progress_callback=conversion_progress,
        )
        created_paths.update(finalized_paths)
        if max_final_zip_bytes > 0:
            # A workspace may already contain a single final ZIP produced before
            # the byte limit was configured.  Recompact here so the limit also
            # applies when a ZIP-backed store has no live files left to archive.
            # Individual members larger than the limit necessarily remain in a
            # one-member shard.
            compacted_paths = compact_archives(
                root,
                mode,
                max_files_per_shard=max(1, max_files_per_shard),
                nonblocking=False,
                compression=compression,
                compression_level=compression_level,
                max_bytes_per_shard=max_final_zip_bytes,
            )
            created_paths.update(compacted_paths)
        _update_storage_conversion_marker(
            root,
            phase="verifying",
            current_subdir="-",
            current_zip="-",
            current_zip_bytes=0,
        )
        if progress_callback is not None:
            progress_callback(
                force=True,
                phase="verifying",
                current_subdir="-",
                subdir="-",
                current_zip="-",
                current_zip_bytes=0,
            )
        final_store = GeneFamilyOutputStore(root)
        # Every newly written archive was already read completely by
        # _verify_zip_crc before its atomic rename.  Recheck the final
        # manifests/indexes here without rereading every payload byte.
        final_store.verify(progress_callback=conversion_progress, deep=False)
        _assert_archive_mode(final_store, mode)
        _update_storage_conversion_marker(
            root,
            phase="complete",
            current_subdir="-",
            current_zip="-",
            current_zip_bytes=0,
            completed_subdirs=sorted(completed_subdirs),
            archived_live_files=archived_removed,
            created_zip_shards=len(created_paths),
        )
        if progress_callback is not None:
            progress_callback(
                force=True,
                phase="complete",
                current_subdir="-",
                subdir="-",
                current_zip="-",
                current_zip_bytes=0,
                archived_live_files=archived_removed,
                created_zip_shards=len(created_paths),
            )
        _finish_storage_conversion(root)
    final_summary, final_unmatched = storage_conversion_summary(
        root,
        valid_family_ids,
        family_from_name,
        workers=workers,
        large_zip_warning_bytes=large_zip_warning_bytes,
        max_final_zip_bytes=max_final_zip_bytes,
    )
    final_summary.update(
        {
            "archived_live_files": archived_removed,
            "created_zip_shards": len([path for path in created_paths if path.is_file()]),
            "unmatched_live_files": len(final_unmatched),
            "conversion_resumed": resumed,
            "compression": compression,
            "compression_level": compression_level,
            "workers": workers,
            "estimated_peak_new_bytes": required_peak_bytes,
            "filesystem_free_bytes": filesystem_free_bytes,
            "effective_available_bytes": effective_available_bytes,
            "large_zip_warning_bytes": int(large_zip_warning_bytes),
            "large_zip_subdir_count": int(summary["large_zip_subdir_count"]),
            "large_zip_subdirs": str(summary["large_zip_subdirs"]),
            "max_final_zip_bytes": int(max_final_zip_bytes),
            "split_zip_subdir_count": int(summary["split_zip_subdir_count"]),
            "split_zip_subdirs": str(summary["split_zip_subdirs"]),
        }
    )
    return final_summary


def convert_storage_to_raw(
    root: Path,
    mode: str,
    *,
    max_files_per_shard: int = 5000,
    pure_raw: bool = False,
    require_resume: bool = False,
    available_bytes: Optional[int] = None,
    progress_callback: Optional[Callable[..., None]] = None,
) -> dict:
    root = Path(root).resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Gene-family output root was not found: {root}")
    if available_bytes is not None and int(available_bytes) < 0:
        raise ArchiveStoreError("available_bytes must be non-negative")
    if pure_raw:
        for generated_status_path in (
            root / "README_GENE_FAMILY_OUTPUTS.txt",
            root / "ARCHIVE_STATUS.tsv",
        ):
            if generated_status_path.is_symlink():
                raise ArchiveStoreError(
                    f"Refusing pure-raw conversion with a symlinked generated status path: {generated_status_path}"
                )
    retired_archive_root = root / PURE_RAW_RETIRED_DIR_NAME
    if pure_raw and (retired_archive_root.exists() or retired_archive_root.is_symlink()):
        if retired_archive_root.is_symlink() or not retired_archive_root.is_dir():
            raise ArchiveStoreError(f"Invalid retired pure-raw metadata path: {retired_archive_root}")
        if any(retired_archive_root.glob("*/*.zip")):
            raise ArchiveStoreError(
                f"A retired pure-raw metadata directory unexpectedly contains ZIP payload: {retired_archive_root}"
            )
        shutil.rmtree(retired_archive_root)
        _fsync_directory(root)
    # Reject quota failures before creating a resumable conversion marker.  A
    # second capacity check under the maintenance locks below protects against
    # concurrent state changes between this read-only preflight and extraction.
    preliminary_store = GeneFamilyOutputStore(root)
    preliminary_store.verify(deep=False)
    _assert_archive_mode(preliminary_store, mode)
    preliminary_summary, _ = storage_conversion_summary(
        root,
        set(),
        lambda _name: None,
    )
    preliminary_required_bytes = int(preliminary_summary["raw_materialize_allocated_bytes"])
    preliminary_filesystem_free = int(shutil.disk_usage(root).free)
    preliminary_effective_available = preliminary_filesystem_free
    if available_bytes is not None:
        if int(available_bytes) < 0:
            raise ArchiveStoreError("available_bytes must be non-negative")
        preliminary_effective_available = min(
            preliminary_effective_available,
            int(available_bytes),
        )
    if preliminary_required_bytes > preliminary_effective_available:
        raise ArchiveStoreError(
            "Insufficient temporary space for ZIP-to-raw conversion: "
            f"allocated_required={preliminary_required_bytes}, "
            f"effective_available_bytes={preliminary_effective_available}. "
            "Pass quota remaining as --available-bytes for a quota-aware preflight."
        )
    restored = 0
    restored_bytes = 0
    with storage_conversion_session(
        root,
        mode,
        "raw",
        require_resume=require_resume,
    ) as resumed:
        _cleanup_partial_archives(_archive_payload_root(root))
        _cleanup_partial_materializations(root)
        preflight_store = GeneFamilyOutputStore(root)
        preflight_store.verify()
        _assert_archive_mode(preflight_store, mode)
        _update_storage_conversion_marker(
            root,
            phase="preflight",
            resumed=resumed,
        )
        if progress_callback is not None:
            progress_callback(force=True, phase="preflight", resumed=resumed)
        archive_root = _archive_state_root(root)
        with all_family_bucket_locks(archive_root), producer_quiescence_lock(archive_root) as readers_idle:
            if not readers_idle:
                raise ArchiveStoreError("Failed to acquire the archive maintenance lock")
            with archive_lock(archive_root):
                _cleanup_partial_materializations(root)
                store = GeneFamilyOutputStore(root)
                store._verify_unlocked()
                _assert_archive_mode(store, mode)
                pending: List[Artifact] = []
                for subdir in store._logical_subdirs_unlocked():
                    for name in store._physical_file_names_unlocked(subdir):
                        artifact = store._physical_artifact_unchecked(subdir, name)
                        if artifact is not None and not artifact.is_live:
                            pending.append(artifact)
                filesystem_stats = os.statvfs(root)
                block_size = int(filesystem_stats.f_frsize or filesystem_stats.f_bsize or 4096)
                required_bytes, required_inodes = _zip_to_raw_requirements(
                    root,
                    pending,
                    block_size=block_size,
                )
                filesystem_free_bytes = int(shutil.disk_usage(root).free)
                effective_available_bytes = filesystem_free_bytes
                if available_bytes is not None:
                    if int(available_bytes) < 0:
                        raise ArchiveStoreError("available_bytes must be non-negative")
                    effective_available_bytes = min(
                        effective_available_bytes,
                        int(available_bytes),
                    )
                if required_bytes > effective_available_bytes:
                    raise ArchiveStoreError(
                        "Insufficient filesystem free space for ZIP-to-raw conversion: "
                        f"allocated_required={required_bytes}, "
                        f"effective_available_bytes={effective_available_bytes}. "
                        "Pass quota remaining as --available-bytes for a "
                        "quota-aware preflight."
                    )
                total_inodes = int(filesystem_stats.f_files)
                available_inodes = int(filesystem_stats.f_favail)
                if (
                    total_inodes > 0
                    and available_inodes >= 0
                    and required_inodes > available_inodes
                ):
                    raise ArchiveStoreError(
                        "Insufficient filesystem inodes for ZIP-to-raw conversion: "
                        f"required={required_inodes}, available={available_inodes}. "
                        "User or project quota may provide fewer inodes than the "
                        "filesystem-wide available value"
                    )
                _update_storage_conversion_marker(
                    root,
                    phase="materializing",
                    pending_files=len(pending),
                    required_bytes=required_bytes,
                    required_inodes=required_inodes,
                    materialized_files=restored,
                )
                if progress_callback is not None:
                    progress_callback(
                        force=True,
                        phase="materializing",
                        pending_files=len(pending),
                        required_bytes=required_bytes,
                        required_inodes=required_inodes,
                        materialized_files=restored,
                    )
                marker_update_started = time.monotonic()
                for artifact in sorted(
                    pending,
                    key=lambda item: item.logical_path,
                ):
                    path = store.materialize(
                        artifact.subdir,
                        artifact.name,
                        _producer_locked=True,
                    )
                    if path.stat().st_size != int(artifact.size or 0):
                        raise ArchiveStoreError(f"Materialized size differs for {artifact.logical_path}")
                    if artifact.sha256 and _sha256_path(path) != artifact.sha256:
                        raise ArchiveStoreError(f"Materialized SHA256 differs for {artifact.logical_path}")
                    restored += 1
                    restored_bytes += int(artifact.size or 0)
                    if progress_callback is not None:
                        progress_callback(
                            phase="materializing",
                            materialized_files=restored,
                            materialized_bytes=restored_bytes,
                            current_path=artifact.logical_path,
                        )
                    if restored % 1000 == 0 or time.monotonic() - marker_update_started >= 30:
                        _update_storage_conversion_marker(
                            root,
                            phase="materializing",
                            materialized_files=restored,
                            materialized_bytes=restored_bytes,
                            current_path=artifact.logical_path,
                        )
                        marker_update_started = time.monotonic()

        # All non-deleted logical artifacts now have authoritative live files.
        # Purge can therefore remove ZIP payload, indexes, and tombstones without
        # changing the visible raw view. Family state and bounded lock metadata remain.
        _update_storage_conversion_marker(
            root,
            phase="removing-archives",
            materialized_files=restored,
            materialized_bytes=restored_bytes,
        )
        if progress_callback is not None:
            progress_callback(
                force=True,
                phase="removing-archives",
                materialized_files=restored,
                materialized_bytes=restored_bytes,
            )
        purge_archives(
            root,
            mode,
            max_files_per_shard=max(1, max_files_per_shard),
        )
        final_store = GeneFamilyOutputStore(root)
        final_store.verify()
        remaining_shards = _physical_archive_paths(root)
        if remaining_shards:
            raise ArchiveStoreError(
                "ZIP shards remained after raw conversion: " + ", ".join(str(path) for path in remaining_shards[:10])
            )
        if pure_raw:
            archive_root = _archive_state_root(root)
            if archive_root.is_symlink():
                raise ArchiveStoreError(f"Refusing to remove symlinked archive root: {archive_root}")
            if archive_root.is_dir():
                os.replace(archive_root, retired_archive_root)
                _fsync_directory(root)
                shutil.rmtree(retired_archive_root)
                _fsync_directory(root)
            for generated_status_path in (
                root / "README_GENE_FAMILY_OUTPUTS.txt",
                root / "ARCHIVE_STATUS.tsv",
            ):
                if generated_status_path.is_symlink():
                    raise ArchiveStoreError(
                        f"Refusing to remove symlinked generated status path: {generated_status_path}"
                    )
                generated_status_path.unlink(missing_ok=True)
            _fsync_directory(root)
        else:
            _update_storage_conversion_marker(
                root,
                phase="complete",
                materialized_files=restored,
                materialized_bytes=restored_bytes,
            )
            _finish_storage_conversion(root)
        if progress_callback is not None:
            progress_callback(
                force=True,
                phase="complete",
                materialized_files=restored,
                materialized_bytes=restored_bytes,
            )
    return {
        "storage": "raw",
        "materialized_files": restored,
        "materialized_bytes": restored_bytes,
        "required_peak_new_bytes": required_bytes,
        "filesystem_free_bytes": filesystem_free_bytes,
        "effective_available_bytes": effective_available_bytes,
        "pure_raw": pure_raw,
        "conversion_resumed": resumed,
    }


def compact_archives(
    root: Path,
    mode: str,
    max_files_per_shard: int = 5000,
    nonblocking: bool = False,
    compression: str = "adaptive",
    compression_level: int = 6,
    max_bytes_per_shard: int = 0,
) -> List[Path]:
    root = Path(root).resolve()
    archive_root = _archive_state_root(root)
    payload_root = _archive_payload_root(root)
    if not archive_root.is_dir():
        return []
    with producer_quiescence_lock(archive_root, nonblocking=nonblocking) as producers_idle:
        if not producers_idle:
            return []
        with archive_lock(archive_root, nonblocking=nonblocking) as acquired:
            if not acquired:
                return []
            _cleanup_partial_archives(payload_root)
            store = GeneFamilyOutputStore(root)
            _assert_archive_mode(store, mode)
            created: List[Path] = []
            for subdir in store._logical_subdirs_unlocked():
                artifacts = list(store._load_subdir_artifacts(subdir).values())
                referenced_shards = {artifact.zip_path for artifact in artifacts if artifact.zip_path is not None}
                logical_bytes = sum(int(artifact.size or 0) for artifact in artifacts)
                oversized_single_archive = (
                    len(referenced_shards) == 1 and max_bytes_per_shard > 0 and logical_bytes > max_bytes_per_shard
                )
                if len(referenced_shards) <= 1 and not oversized_single_archive:
                    continue
                ordered = sorted(artifacts, key=lambda artifact: artifact.logical_path)
                artifact_chunks = _balanced_artifact_chunks(
                    ordered,
                    max_files_per_shard,
                    max_bytes_per_shard,
                )
                compacted_all: List[Artifact] = []
                final_archive = _final_archive_path(root, subdir)
                keep_final_archive = final_archive.resolve() in {path.resolve() for path in referenced_shards} and not (
                    max_bytes_per_shard > 0 and logical_bytes > max_bytes_per_shard
                )
                retain_named_parts = max_bytes_per_shard > 0 and logical_bytes > max_bytes_per_shard
                if keep_final_archive:
                    generation = store._next_generation()
                    zip_path, compacted = _compact_artifact_chunk(
                        payload_root,
                        subdir,
                        ordered,
                        mode,
                        generation,
                        compression=compression,
                        compression_level=compression_level,
                        destination_path=final_archive,
                    )
                    created.append(zip_path)
                    compacted_all.extend(compacted)
                else:
                    for artifact_chunk in artifact_chunks:
                        generation = store._next_generation()
                        destination_path = (
                            payload_root / subdir / f"{subdir}.part-{generation:06d}.zip"
                            if retain_named_parts
                            else None
                        )
                        zip_path, compacted = _compact_artifact_chunk(
                            payload_root,
                            subdir,
                            artifact_chunk,
                            mode,
                            generation,
                            compression=compression,
                            compression_level=compression_level,
                            destination_path=destination_path,
                        )
                        created.append(zip_path)
                        compacted_all.extend(compacted)
                updated_subdir = {artifact.logical_path: artifact for artifact in compacted_all}
                store._merge_index_subdirs({subdir: updated_subdir})
                referenced_paths = {
                    artifact.zip_path.resolve() for artifact in updated_subdir.values() if artifact.zip_path is not None
                }
                _remove_unreferenced_subdir_archives(
                    root,
                    payload_root,
                    subdir,
                    referenced_paths,
                )
            return [path for path in created if path.is_file()]


def _write_archive_status(
    root: Path,
    store: Optional[GeneFamilyOutputStore] = None,
    *,
    nonblocking: bool = False,
) -> Optional[Path]:
    root = Path(root).resolve()
    store = store or GeneFamilyOutputStore(root)
    rows: List[Tuple[str, str, int, int, int, str]] = []
    with producer_read_lock(
        store.archive_root,
        nonblocking=nonblocking,
    ) as acquired:
        if not acquired:
            return None
        store._refresh_if_index_changed()
        for subdir in store._logical_subdirs_unlocked():
            artifacts = store._load_subdir_artifacts(subdir)
            live_dir = root / subdir
            live_files = (
                sum(
                    1
                    for path in live_dir.iterdir()
                    if path.is_file() and not path.is_symlink() and not path.name.startswith(".")
                )
                if live_dir.is_dir() and not live_dir.is_symlink()
                else 0
            )
            zip_paths = sorted(
                {artifact.zip_path.resolve() for artifact in artifacts.values() if artifact.zip_path is not None}
            )
            final_path = _final_archive_path(root, subdir).resolve()
            if zip_paths and set(zip_paths) == {final_path}:
                storage = "finalized"
            elif zip_paths and final_path in zip_paths:
                storage = "base+parts"
            elif zip_paths:
                storage = "parts"
            else:
                storage = "raw"
            if zip_paths and live_files:
                storage += "+live"
            location = ",".join(path.relative_to(root).as_posix() for path in zip_paths)
            rows.append((subdir, storage, live_files, len(artifacts), len(zip_paths), location))
    status_path = root / "ARCHIVE_STATUS.tsv"
    temporary = status_path.with_name(f".{status_path.name}.partial.{os.getpid()}.{uuid.uuid4().hex}")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            handle.write("logical_subdirectory\tstorage\tlive_files\tarchived_files\tzip_files\tlocation\n")
            for row in rows:
                handle.write("\t".join(str(value) for value in row) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, status_path)
        _fsync_directory(root)
    finally:
        if temporary.exists():
            temporary.unlink()
    return status_path


def _remove_unreferenced_subdir_archives(
    root: Path,
    payload_root: Path,
    subdir: str,
    referenced: Set[Path],
) -> None:
    candidates = list((payload_root / subdir).glob("*.zip"))
    final_path = _final_archive_path(root, subdir)
    if final_path.is_file() and _zip_has_genegalleon_manifest(final_path):
        candidates.append(final_path)
    for path in sorted(set(candidates)):
        if path.resolve() not in referenced:
            path.unlink()
    shard_dir = payload_root / subdir
    try:
        if shard_dir.is_dir() and not any(shard_dir.iterdir()):
            shard_dir.rmdir()
    except OSError:
        pass
    try:
        if payload_root.is_dir() and not any(payload_root.iterdir()):
            payload_root.rmdir()
    except OSError:
        pass


def _finalize_subdirs_locked(
    root: Path,
    mode: str,
    family_ids: Sequence[str],
    family_from_name: Callable[[str], Optional[str]],
    subdirs: Optional[Iterable[str]] = None,
    *,
    compression: str = "adaptive",
    compression_level: int = 6,
    max_final_zip_bytes: int = 0,
    progress_callback: Optional[Callable[..., None]] = None,
) -> List[Path]:
    root = Path(root).resolve()
    payload_root = _archive_payload_root(root)
    store = GeneFamilyOutputStore(root)
    _assert_archive_mode(store, mode)
    catalog = store._load_index_catalog()
    archived_subdirs = (
        sorted(str(value) for value in catalog["subdirs"])
        if catalog is not None
        else (store._load_archives() or sorted({artifact.subdir for artifact in store._archived.values()}))
    )
    selected = sorted(set(str(value) for value in subdirs)) if subdirs is not None else archived_subdirs
    valid_family_ids = set(family_ids)
    created: List[Path] = []
    finalized_count = 0
    total_subdirs = len(selected)
    for subdir in selected:
        _safe_logical_path(subdir, "__archive__")
        artifacts = list(store._load_subdir_artifacts(subdir).values())
        if not artifacts:
            continue
        logical_bytes = sum(int(artifact.size or 0) for artifact in artifacts)
        if max_final_zip_bytes > 0 and logical_bytes > max_final_zip_bytes:
            if progress_callback is not None:
                progress_callback(
                    force=True,
                    phase="retaining-parts",
                    subdir=subdir,
                    current_subdir=subdir,
                    current_zip="-",
                    finalized_subdirs=finalized_count,
                    total_subdirs=total_subdirs,
                    logical_bytes=logical_bytes,
                    max_final_zip_bytes=max_final_zip_bytes,
                )
            continue
        archived_logical_paths = {artifact.logical_path for artifact in artifacts}
        live_dir = root / subdir
        live_files = (
            [
                path
                for path in sorted(live_dir.iterdir())
                if path.is_file() and not path.is_symlink() and not path.name.startswith(".")
            ]
            if live_dir.is_dir() and not live_dir.is_symlink()
            else []
        )
        blocking_live_files = [
            path
            for path in live_files
            if _safe_logical_path(subdir, path.name) in archived_logical_paths
            or family_from_name(path.name) in valid_family_ids
        ]
        if blocking_live_files:
            raise ArchiveStoreError(
                f"Cannot finalize {subdir} while live output files remain: "
                + ", ".join(str(path) for path in blocking_live_files[:10])
            )
        final_path = _final_archive_path(root, subdir)
        if progress_callback is not None:
            progress_callback(
                force=True,
                phase="finalizing",
                subdir=subdir,
                current_subdir=subdir,
                current_zip=final_path.name,
                finalized_subdirs=finalized_count,
                total_subdirs=total_subdirs,
                verify_members_completed=0,
                verify_members_total=0,
                verify_bytes_completed=0,
                verify_bytes_total=0,
                current_zip_bytes=(final_path.stat().st_size if final_path.is_file() else 0),
            )
        referenced_before = {artifact.zip_path.resolve() for artifact in artifacts if artifact.zip_path is not None}
        if referenced_before == {final_path.resolve()}:
            _remove_unreferenced_subdir_archives(
                root,
                payload_root,
                subdir,
                referenced_before,
            )
            finalized_count += 1
            if progress_callback is not None:
                progress_callback(
                    force=True,
                    phase="finalizing",
                    subdir=subdir,
                    current_subdir=subdir,
                    current_zip=final_path.name,
                    finalized_subdirs=finalized_count,
                    total_subdirs=total_subdirs,
                    verify_members_completed=0,
                    verify_members_total=0,
                    verify_bytes_completed=0,
                    verify_bytes_total=0,
                    current_zip_bytes=final_path.stat().st_size,
                )
            continue
        generation = store._next_generation()

        def report_zip_bytes(
            byte_count: int,
            *,
            _subdir: str = subdir,
            _final_path: Path = final_path,
            _finalized_count: int = finalized_count,
        ) -> None:
            if progress_callback is not None:
                progress_callback(
                    force=False,
                    phase="finalizing",
                    subdir=_subdir,
                    current_subdir=_subdir,
                    current_zip=_final_path.name,
                    finalized_subdirs=_finalized_count,
                    total_subdirs=total_subdirs,
                    verify_members_completed=0,
                    verify_members_total=0,
                    verify_bytes_completed=0,
                    verify_bytes_total=0,
                    current_zip_bytes=byte_count,
                )

        def report_zip_verification(
            members_completed: int,
            members_total: int,
            bytes_completed: int,
            bytes_total: int,
            *,
            _subdir: str = subdir,
            _final_path: Path = final_path,
            _finalized_count: int = finalized_count,
        ) -> None:
            if progress_callback is not None:
                progress_callback(
                    force=(members_completed == members_total),
                    phase="verifying-final-zip",
                    subdir=_subdir,
                    current_subdir=_subdir,
                    current_zip=_final_path.name,
                    finalized_subdirs=_finalized_count,
                    total_subdirs=total_subdirs,
                    verify_members_completed=members_completed,
                    verify_members_total=members_total,
                    verify_bytes_completed=bytes_completed,
                    verify_bytes_total=bytes_total,
                    current_zip_bytes=(_final_path.stat().st_size if _final_path.is_file() else 0),
                )

        _, finalized = _compact_artifact_chunk(
            payload_root,
            subdir,
            sorted(artifacts, key=lambda artifact: artifact.logical_path),
            mode,
            generation,
            compression=compression,
            compression_level=compression_level,
            destination_path=final_path,
            byte_progress_callback=report_zip_bytes,
            verification_progress_callback=report_zip_verification,
        )
        updated = {artifact.logical_path: artifact for artifact in finalized}
        store._merge_index_subdirs({subdir: updated})
        referenced = {artifact.zip_path.resolve() for artifact in updated.values() if artifact.zip_path is not None}
        _remove_unreferenced_subdir_archives(
            root,
            payload_root,
            subdir,
            referenced,
        )
        created.append(final_path)
        finalized_count += 1
        if progress_callback is not None:
            progress_callback(
                force=True,
                phase="finalizing",
                subdir=subdir,
                current_subdir=subdir,
                current_zip=final_path.name,
                finalized_subdirs=finalized_count,
                total_subdirs=total_subdirs,
                verify_members_completed=0,
                verify_members_total=0,
                verify_bytes_completed=0,
                verify_bytes_total=0,
                current_zip_bytes=final_path.stat().st_size,
            )
    return created


def finalize_archives(
    root: Path,
    mode: str,
    family_ids: Sequence[str],
    *,
    subdirs: Optional[Iterable[str]] = None,
    nonblocking: bool = False,
    compression: str = "adaptive",
    compression_level: int = 6,
    max_final_zip_bytes: int = 0,
    progress_callback: Optional[Callable[..., None]] = None,
) -> List[Path]:
    root = Path(root).resolve()
    archive_root = _archive_state_root(root)
    if not archive_root.is_dir():
        return []
    finalized: List[Path] = []
    family_matchers = query_id_matchers(family_ids)
    family_from_name = (
        orthogroup_id_from_name if mode == "orthogroup" else lambda name: query_id_from_name(name, family_matchers)
    )
    with lock_available_family_ids(
        archive_root,
        family_ids,
        nonblocking=nonblocking,
    ) as lockable_family_ids:
        if set(lockable_family_ids) != set(family_ids):
            return []
        with producer_quiescence_lock(
            archive_root,
            nonblocking=nonblocking,
        ) as producers_idle:
            if not producers_idle:
                return []
            with archive_lock(archive_root, nonblocking=nonblocking) as acquired:
                if not acquired:
                    return []
                _cleanup_partial_archives(_archive_payload_root(root))
                finalized = _finalize_subdirs_locked(
                    root,
                    mode,
                    family_ids,
                    family_from_name,
                    subdirs,
                    compression=compression,
                    compression_level=compression_level,
                    max_final_zip_bytes=max_final_zip_bytes,
                    progress_callback=progress_callback,
                )
    _write_archive_status(root)
    return finalized


def migrate_archive_layout(root: Path) -> List[Tuple[Path, Path]]:
    """Move legacy mixed metadata/payload stores into the visible layout."""

    root = Path(root).resolve()
    legacy_root = root / LEGACY_ARCHIVE_DIR_NAME
    state_root = root / STORE_DIR_NAME
    payload_root = root / ACTIVE_ARCHIVE_DIR_NAME
    if legacy_root.is_symlink() or state_root.is_symlink() or payload_root.is_symlink():
        raise ArchiveStoreError("Symlinked archive layout roots are not supported")
    if legacy_root.exists() and state_root.exists():
        raise ArchiveStoreError(f"Both legacy and current archive metadata roots exist: {legacy_root}, {state_root}")
    migrated_legacy = legacy_root.is_dir()
    if migrated_legacy:
        legacy_store = GeneFamilyOutputStore(root)
        legacy_store.verify()
        with all_family_bucket_locks(legacy_root), producer_quiescence_lock(legacy_root) as readers_idle:
            if not readers_idle:
                raise ArchiveStoreError("Failed to acquire the archive maintenance lock")
            with archive_lock(legacy_root):
                os.replace(legacy_root, state_root)
                _fsync_directory(root)
    if not state_root.is_dir():
        return []

    embedded_zip_paths = [
        path
        for directory in state_root.iterdir()
        if directory.is_dir() and not directory.is_symlink()
        for path in directory.glob("*.zip")
    ]
    if not migrated_legacy and not embedded_zip_paths:
        _write_archive_readme(root)
        _write_archive_status(root)
        return []

    payload_root.mkdir(parents=True, exist_ok=True)
    moved: List[Tuple[Path, Path]] = []
    marker = state_root / INDEX_UPDATE_FILE
    if not marker.exists():
        marker.write_text(
            json.dumps(
                {
                    "schema_version": INDEX_SCHEMA_VERSION,
                    "pid": os.getpid(),
                    "hostname": os.uname().nodename,
                    "created_ns": time.time_ns(),
                    "layout_migration": True,
                },
                sort_keys=True,
                separators=(",", ":"),
            )
            + "\n",
            encoding="utf-8",
        )
        _fsync_directory(state_root)
    for source_dir in sorted(state_root.iterdir()):
        if not source_dir.is_dir() or source_dir.is_symlink():
            continue
        zip_paths = sorted(source_dir.glob("*.zip"))
        if not zip_paths:
            continue
        destination_dir = payload_root / source_dir.name
        if destination_dir.exists():
            if not destination_dir.is_dir() or destination_dir.is_symlink():
                raise ArchiveStoreError(f"Invalid visible archive directory: {destination_dir}")
        else:
            destination_dir.mkdir(parents=True)
        for source_path in zip_paths:
            with zipfile.ZipFile(source_path, "r") as archive:
                with archive.open(MANIFEST_MEMBER, "r") as handle:
                    manifest = json.load(handle)
            generation = int(manifest["generation"])
            kind = "pack" if source_path.name.startswith("pack-") else "part"
            destination_path = destination_dir / (f"{source_dir.name}.{kind}-{generation:06d}.zip")
            if destination_path.exists():
                raise ArchiveStoreError(f"Visible archive destination already exists: {destination_path}")
            os.replace(source_path, destination_path)
            moved.append((source_path, destination_path))
        try:
            source_dir.rmdir()
        except OSError:
            pass
        _fsync_directory(destination_dir)
    _fsync_directory(payload_root)
    repair_archive_index(root, remove_orphans=False)
    _write_archive_readme(root)
    _write_archive_status(root)
    return moved


def repair_archive_index(
    root: Path,
    *,
    remove_orphans: bool = False,
    progress_callback: Optional[Callable[..., None]] = None,
) -> List[Path]:
    root = Path(root).resolve()
    archive_root = _archive_state_root(root)
    payload_root = _archive_payload_root(root)
    if not archive_root.is_dir() and not _physical_archive_paths(root):
        return []
    with producer_quiescence_lock(archive_root) as readers_idle:
        if not readers_idle:
            raise ArchiveStoreError("Failed to acquire the archive maintenance lock")
        with archive_lock(archive_root):
            _cleanup_partial_archives(payload_root)
            store = GeneFamilyOutputStore(root)
            rebuilt: Dict[str, Artifact] = {}
            rebuilt_ranks: Dict[str, Tuple[int, int, str]] = {}
            mode_by_path: Dict[Path, str] = {}
            physical_paths = _physical_archive_paths(root)
            physical_sizes = {path: path.stat().st_size for path in physical_paths}
            processed_zip_bytes = 0
            verified_members = 0
            verified_bytes = 0
            last_progress = time.monotonic()
            if progress_callback is not None:
                progress_callback(
                    force=True,
                    phase="repairing-index",
                    repair_zip_files_completed=0,
                    repair_zip_files_total=len(physical_paths),
                    repair_zip_bytes_processed=0,
                    repair_zip_bytes_total=sum(physical_sizes.values()),
                    repair_members_verified=0,
                    repair_bytes_verified=0,
                    current_zip="-",
                )
            for zip_index, zip_path in enumerate(physical_paths):
                if zip_path.is_symlink() or zip_path.parent.is_symlink():
                    raise ArchiveStoreError(f"Symlinked ZIP shards are not supported: {zip_path}")
                manifest = store._read_manifest(zip_path)
                mode_by_path[zip_path.resolve()] = str(manifest["mode"])
                manifest_generation = int(manifest["generation"])
                with zipfile.ZipFile(zip_path, "r") as archive:
                    for member in manifest["members"]:
                        digest = hashlib.sha256()
                        member_name = str(member["member_name"])
                        with archive.open(member_name, "r") as handle:
                            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                                digest.update(chunk)
                                verified_bytes += len(chunk)
                                now = time.monotonic()
                                if progress_callback is not None and now - last_progress >= 30.0:
                                    progress_callback(
                                        force=False,
                                        phase="repairing-index",
                                        repair_zip_files_completed=zip_index,
                                        repair_zip_files_total=len(physical_paths),
                                        repair_zip_bytes_processed=processed_zip_bytes,
                                        repair_zip_bytes_total=sum(physical_sizes.values()),
                                        repair_members_verified=verified_members,
                                        repair_bytes_verified=verified_bytes,
                                        current_zip=zip_path.name,
                                    )
                                    last_progress = now
                        if digest.hexdigest() != str(member.get("sha256", "")):
                            raise ArchiveStoreError(f"SHA256 verification failed in {zip_path}: {member_name}")
                        verified_members += 1
                        logical_path = str(member["logical_path"])
                        subdir, name = logical_path.split("/", 1)
                        source_signature = member.get("source_signature")
                        mtime_ns = member.get("mtime_ns")
                        if mtime_ns is None and isinstance(source_signature, list) and len(source_signature) >= 4:
                            mtime_ns = source_signature[3]
                        artifact = Artifact(
                            logical_path=logical_path,
                            subdir=subdir,
                            name=name,
                            generation=int(member.get("generation", manifest["generation"])),
                            zip_path=zip_path,
                            member_name=member_name,
                            size=int(member["size"]),
                            crc=int(member["crc"]),
                            sha256=str(member.get("sha256", "")) or None,
                            mtime_ns=int(mtime_ns) if mtime_ns is not None else None,
                            family_id=(str(member["family_id"]) if member.get("family_id") not in {None, ""} else None),
                        )
                        previous = rebuilt.get(logical_path)
                        candidate_rank = (
                            artifact.generation,
                            manifest_generation,
                            str(zip_path),
                        )
                        if previous is None or candidate_rank > rebuilt_ranks[logical_path]:
                            rebuilt[logical_path] = artifact
                            rebuilt_ranks[logical_path] = candidate_rank
                processed_zip_bytes += physical_sizes[zip_path]
                if progress_callback is not None:
                    progress_callback(
                        force=True,
                        phase="repairing-index",
                        repair_zip_files_completed=zip_index + 1,
                        repair_zip_files_total=len(physical_paths),
                        repair_zip_bytes_processed=processed_zip_bytes,
                        repair_zip_bytes_total=sum(physical_sizes.values()),
                        repair_members_verified=verified_members,
                        repair_bytes_verified=verified_bytes,
                        current_zip=zip_path.name,
                    )
            referenced = {artifact.zip_path.resolve() for artifact in rebuilt.values() if artifact.zip_path is not None}
            referenced_modes = {mode_by_path[path] for path in referenced if path in mode_by_path}
            if len(referenced_modes) > 1:
                raise ArchiveStoreError(
                    "Cannot repair indexes whose referenced ZIP shards use mixed "
                    "gene-family modes: " + ", ".join(sorted(referenced_modes))
                )
            store._write_index(rebuilt, recover_pending=True)
            if referenced_modes and _read_store_metadata(root) is None:
                repaired_mode = next(iter(referenced_modes))
                _write_store_metadata(
                    root,
                    repaired_mode,
                    (artifact.family_id for artifact in rebuilt.values() if artifact.family_id is not None),
                )
            orphaned = [path for path in physical_paths if path.resolve() not in referenced]
            if remove_orphans:
                for path in orphaned:
                    path.unlink()
            return orphaned


def purge_archives(
    root: Path,
    mode: str,
    *,
    max_files_per_shard: int = 5000,
    valid_family_ids: Optional[Set[str]] = None,
    family_from_name: Optional[Callable[[str], Optional[str]]] = None,
    compression: str = "adaptive",
    compression_level: int = 6,
) -> List[Path]:
    root = Path(root).resolve()
    archive_root = _archive_state_root(root)
    payload_root = _archive_payload_root(root)
    if not archive_root.is_dir():
        return []
    # Finish one-time migration before taking the maintenance lock below.
    # This avoids recursively acquiring archive.lock and gives purge a single
    # bucketed state representation to prune.
    GeneFamilyOutputStore(root)._ensure_state_index()
    with all_family_bucket_locks(archive_root), producer_quiescence_lock(archive_root) as readers_idle:
        if not readers_idle:
            raise ArchiveStoreError("Failed to acquire the archive maintenance lock")
        with archive_lock(archive_root):
            store = GeneFamilyOutputStore(root)
            _assert_archive_mode(store, mode)
            store._load_tombstones()
            created: List[Path] = []
            archived_subdirs = list(store._logical_subdirs_unlocked())
            for subdir in archived_subdirs:
                archived = store._load_subdir_artifacts(subdir)
                retained: List[Artifact] = []
                for logical_path, artifact in archived.items():
                    if any(
                        tombstone is not None and tombstone[1] == "delete" and tombstone[0] >= artifact.generation
                        for equivalent_path in _equivalent_output_logical_paths(
                            artifact.subdir,
                            artifact.name,
                        )
                        if (tombstone := store._tombstones.get(equivalent_path)) is not None
                    ):
                        continue
                    family_id = artifact.family_id
                    if family_id is None and family_from_name is not None:
                        family_id = family_from_name(artifact.name)
                    if valid_family_ids is not None and family_id not in valid_family_ids:
                        continue
                    if (root / logical_path).is_file():
                        continue
                    retained.append(artifact)
                ordered = sorted(retained, key=lambda item: item.logical_path)
                remaining_subdir: Dict[str, Artifact] = {}
                for start in range(0, len(ordered), max(1, max_files_per_shard)):
                    generation = store._next_generation()
                    zip_path, compacted = _compact_artifact_chunk(
                        payload_root,
                        subdir,
                        ordered[start : start + max(1, max_files_per_shard)],
                        mode,
                        generation,
                        compression=compression,
                        compression_level=compression_level,
                    )
                    created.append(zip_path)
                    for artifact in compacted:
                        remaining_subdir[artifact.logical_path] = artifact
                store._merge_index_subdirs({subdir: remaining_subdir})
                referenced = {
                    artifact.zip_path.resolve()
                    for artifact in remaining_subdir.values()
                    if artifact.zip_path is not None
                }
                physical_subdir_paths = list((payload_root / subdir).glob("*.zip"))
                final_path = _final_archive_path(root, subdir)
                if final_path.is_file() and _zip_has_genegalleon_manifest(final_path):
                    physical_subdir_paths.append(final_path)
                for path in sorted(physical_subdir_paths):
                    if path.resolve() not in referenced:
                        path.unlink()
                _fsync_directory(payload_root / subdir)
            tombstone_path = archive_root / TOMBSTONE_FILE
            if tombstone_path.is_file() and not tombstone_path.is_symlink():
                tombstone_path.unlink()
            legacy_state_path = archive_root / FAMILY_STATE_FILE
            if legacy_state_path.is_file() and not legacy_state_path.is_symlink():
                legacy_state_path.unlink()
            if valid_family_ids is not None:
                state_dir = archive_root / FAMILY_STATE_DIR_NAME
                for state_path in sorted(state_dir.glob("*.json")):
                    bucket = state_path.stem
                    if re.fullmatch(r"[0-9a-f]{2}", bucket) is None:
                        continue
                    lock_path = archive_root / FAMILY_STATE_LOCK_DIR_NAME / f"{bucket}.lock"
                    with _bucket_lock(lock_path, exclusive=True):
                        states = store._read_state_bucket(state_path)
                        retained_states = {
                            family_id: state for family_id, state in states.items() if family_id in valid_family_ids
                        }
                        if retained_states:
                            store._write_state_bucket(state_path, retained_states)
                        elif state_path.is_file() and not state_path.is_symlink():
                            state_path.unlink()
            store._write_index_epoch()
            _fsync_directory(archive_root)
            store._reset_cache()
            return [path for path in created if path.is_file()]


def _sha256_path(path: Path) -> str:
    return cached_sha256_file(path)[0]


def cleanup_materialization_receipt(
    receipt_path: Path,
    *,
    nonblocking: bool = False,
    _family_locked: bool = False,
) -> List[Path]:
    raw_receipt_path = Path(receipt_path)
    if raw_receipt_path.is_symlink():
        raise ArchiveStoreError(f"Symlinked materialization receipts are not supported: {raw_receipt_path}")
    receipt_path = raw_receipt_path.resolve()
    if not receipt_path.is_file():
        return []
    records: List[dict] = []
    try:
        with receipt_path.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                record = json.loads(line)
                if not isinstance(record, dict):
                    raise ArchiveStoreError(
                        f"Materialization receipt record is not a JSON object at {receipt_path}:{line_number}"
                    )
                if int(record.get("schema_version", -1)) != ARCHIVE_SCHEMA_VERSION:
                    raise ArchiveStoreError(
                        f"Unsupported materialization receipt schema at {receipt_path}:{line_number}"
                    )
                records.append(record)
    except ArchiveStoreError:
        raise
    except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
        raise ArchiveStoreError(f"Failed to read materialization receipt {receipt_path}: {exc}") from exc
    if not records:
        receipt_path.unlink()
        return []
    try:
        roots = {Path(str(record["root"])).resolve() for record in records}
        family_ids = {str(record["family_id"]) for record in records}
        for record in records:
            logical_path = str(record["logical_path"])
            subdir, name = logical_path.split("/", 1)
            _safe_logical_path(subdir, name)
    except (KeyError, TypeError, ValueError) as exc:
        raise ArchiveStoreError(
            f"Materialization receipt has an invalid required field: {receipt_path}: {exc}"
        ) from exc
    if len(roots) != 1 or len(family_ids) != 1 or "" in family_ids:
        raise ArchiveStoreError(f"Materialization receipt mixes roots or families: {receipt_path}")
    root = next(iter(roots))
    family_id = next(iter(family_ids))
    expected_tmp_root = (root / "tmp").resolve()
    if (
        receipt_path.name != MATERIALIZATION_RECEIPT_NAME
        or receipt_path.parent.parent != expected_tmp_root
        or re.fullmatch(r"[0-9]+_.+", receipt_path.parent.name) is None
        or receipt_path.parent.name.split("_", 1)[1] != family_id
    ):
        raise ArchiveStoreError(f"Materialization receipt is outside its family task directory: {receipt_path}")

    @contextlib.contextmanager
    def lock_context() -> Iterator[bool]:
        if _family_locked:
            yield True
        else:
            with family_bucket_lock(
                _archive_state_root(root),
                family_id,
                exclusive=True,
                nonblocking=nonblocking,
            ) as acquired:
                yield acquired

    removed: List[Path] = []
    with lock_context() as acquired:
        if not acquired:
            return []
        touched_dirs: Set[Path] = set()
        for record in records:
            logical_path = str(record["logical_path"])
            subdir, name = logical_path.split("/", 1)
            _safe_logical_path(subdir, name)
            path = root / subdir / name
            if not path.is_file() or path.is_symlink():
                continue
            expected_size = int(record.get("size", -1))
            expected_sha256 = str(record.get("sha256", ""))
            if expected_size < 0 or re.fullmatch(r"[0-9a-f]{64}", expected_sha256) is None:
                continue
            if path.stat().st_size != expected_size:
                continue
            if _sha256_path(path) != expected_sha256:
                continue
            path.unlink()
            removed.append(path)
            touched_dirs.add(path.parent)
        receipt_path.unlink()
        for directory in sorted(touched_dirs, reverse=True):
            try:
                if directory.is_dir() and not any(directory.iterdir()):
                    directory.rmdir()
            except OSError:
                pass
    return removed


def cleanup_stale_tmp(
    root: Path,
    older_than_days: float,
    max_directories: int = 0,
    max_bytes: int = 0,
    max_files: int = 0,
    nonblocking: bool = False,
) -> List[Path]:
    if older_than_days < 0:
        raise ValueError("older_than_days must be non-negative")
    if max_directories < 0:
        raise ValueError("max_directories must be non-negative")
    if max_bytes < 0 or max_files < 0:
        raise ValueError("max_bytes and max_files must be non-negative")
    if older_than_days == 0 and max_directories == 0 and max_bytes == 0 and max_files == 0:
        return []
    root = Path(root).resolve()
    tmp_root = root / "tmp"
    if not tmp_root.is_dir() or tmp_root.is_symlink():
        return []
    cutoff_ns = time.time_ns() - int(older_than_days * 86400 * 1_000_000_000) if older_than_days > 0 else None
    removed: List[Path] = []
    archive_root = _archive_state_root(root)
    task_dirs = [
        candidate
        for candidate in sorted(tmp_root.iterdir())
        if candidate.is_dir() and not candidate.is_symlink() and re.fullmatch(r"[0-9]+_.+", candidate.name) is not None
    ]
    family_id_by_candidate = {candidate: candidate.name.split("_", 1)[1] for candidate in task_dirs}
    with lock_available_family_ids(
        archive_root,
        family_id_by_candidate.values(),
        nonblocking=nonblocking,
    ) as inactive_family_ids:
        candidates: List[Tuple[Path, int, int, int]] = []
        for candidate in task_dirs:
            if family_id_by_candidate[candidate] not in inactive_family_ids:
                continue
            try:
                newest_mtime_ns = candidate.stat().st_mtime_ns
                total_bytes = 0
                total_files = 0
                for descendant in candidate.rglob("*"):
                    if descendant.is_symlink():
                        continue
                    descendant_stat = descendant.stat()
                    newest_mtime_ns = max(newest_mtime_ns, descendant_stat.st_mtime_ns)
                    if descendant.is_file():
                        total_files += 1
                        total_bytes += int(descendant_stat.st_size)
            except FileNotFoundError:
                continue
            candidates.append((candidate, newest_mtime_ns, total_bytes, total_files))
        keep_within_limits: Set[Path] = set()
        kept_bytes = 0
        kept_files = 0
        for candidate, _, candidate_bytes, candidate_files in sorted(
            candidates,
            key=lambda item: (-item[1], item[0].name),
        ):
            if max_directories > 0 and len(keep_within_limits) >= max_directories:
                continue
            if max_bytes > 0 and kept_bytes + candidate_bytes > max_bytes:
                continue
            if max_files > 0 and kept_files + candidate_files > max_files:
                continue
            keep_within_limits.add(candidate)
            kept_bytes += candidate_bytes
            kept_files += candidate_files
        for candidate, newest_mtime_ns, _, _ in candidates:
            expired = cutoff_ns is not None and newest_mtime_ns < cutoff_ns
            limited = max_directories > 0 or max_bytes > 0 or max_files > 0
            excess = limited and candidate not in keep_within_limits
            if not expired and not excess:
                continue
            receipt_path = candidate / MATERIALIZATION_RECEIPT_NAME
            if receipt_path.is_file():
                cleanup_materialization_receipt(
                    receipt_path,
                    nonblocking=False,
                    _family_locked=True,
                )
            shutil.rmtree(candidate)
            removed.append(candidate)
        try:
            if tmp_root.is_dir() and not any(tmp_root.iterdir()):
                tmp_root.rmdir()
        except OSError:
            pass
    return removed


def _parse_subdirs(value: Optional[str]) -> Optional[Set[str]]:
    if value is None or not value.strip():
        return None
    return {item.strip() for item in value.split(",") if item.strip()}


def build_parser() -> argparse.ArgumentParser:
    from gene_family_output_cli import build_parser as build_cli_parser

    return build_cli_parser(
        description=__doc__ or "Gene-family output storage",
        default_large_zip_warning_bytes=DEFAULT_LARGE_ZIP_WARNING_BYTES,
    )


def _resolve_cli_mode(
    args: argparse.Namespace,
    root: Path,
    *,
    required: bool = True,
) -> Optional[str]:
    query_dir = getattr(args, "query_dir", None)
    genecount = getattr(args, "genecount", None)
    if query_dir is not None and genecount is not None:
        raise ValueError("--query-dir and --genecount cannot be used together")
    catalog_mode = "query2family" if query_dir is not None else ("orthogroup" if genecount is not None else None)
    explicit_mode = getattr(args, "mode", None)
    if explicit_mode is not None and catalog_mode is not None and explicit_mode != catalog_mode:
        raise ValueError(f"--mode {explicit_mode} conflicts with the supplied {catalog_mode} catalog")
    return resolve_archive_mode(
        root,
        explicit_mode=explicit_mode or catalog_mode,
        required=required,
    )


def _catalog_sources(args: argparse.Namespace) -> List[Path]:
    return [
        Path(path)
        for path in (
            getattr(args, "query_dir", None),
            getattr(args, "genecount", None),
            getattr(args, "family_id_file", None),
        )
        if path is not None
    ]


def _validate_zip_write_options(args: argparse.Namespace) -> None:
    if not 0 <= int(args.compression_level) <= 9:
        raise ValueError("--compression-level must be between 0 and 9")
    if not 1 <= int(args.workers) <= 4:
        raise ValueError("--workers must be between 1 and 4")


def run_cli(args: argparse.Namespace) -> int:
    if args.command == "cleanup-materialized":
        for removed_path in cleanup_materialization_receipt(
            args.receipt,
            nonblocking=args.nonblocking,
        ):
            print(f"removed-materialized\t{removed_path}")
        return 0

    root = args.root.resolve()
    if args.command == "conversion-status":
        status_payload = storage_conversion_status(root)
        if args.json:
            print(json.dumps(status_payload, indent=2, sort_keys=True))
        else:
            for key, value in status_payload.items():
                print(f"{key}\t{value}")
        return 0
    if args.command == "optimize-metadata":
        for key, value in optimize_archive_metadata(root).items():
            print(f"{key}\t{value}")
        return 0
    if args.command == "migrate-layout":
        for old_path, new_path in migrate_archive_layout(root):
            print(f"moved\t{old_path}\t{new_path}")
        return 0
    if args.command == "refresh-status":
        print(f"status\t{_write_archive_status(root)}")
        return 0
    if args.command == "archive-completed":
        _validate_zip_write_options(args)
        if args.progress_interval < 0:
            raise ValueError("--progress-interval must be non-negative")
        if args.max_final_zip_bytes < 0:
            raise ValueError("--max-final-zip-bytes must be non-negative")
        mode = _resolve_cli_mode(args, root)
        assert mode is not None
        family_ids, family_from_name = family_context(
            mode,
            query_dir=args.query_dir,
            genecount=args.genecount,
        )
        results = archive_completed_outputs(
            root=root,
            mode=mode,
            family_ids=family_ids,
            family_from_name=family_from_name,
            min_files=max(1, args.min_files),
            max_files_per_shard=max(1, args.max_files_per_shard),
            nonblocking=args.nonblocking,
            compression=args.compression,
            compression_level=args.compression_level,
            workers=args.workers,
            catalog_sources=_catalog_sources(args),
            catalog_family_ids=family_ids,
            max_bytes_per_shard=args.max_final_zip_bytes,
        )
        for zip_path, removed in results:
            print(f"archived\t{removed}\t{zip_path}")
        _write_archive_status(root, nonblocking=args.nonblocking)
        if mode == "orthogroup":
            completion_store = GeneFamilyOutputStore(root)
            if completed_family_ids(completion_store, family_ids) == set(family_ids):
                reporter = ProgressReporter(args.progress_interval)
                reporter.start()
                try:
                    finalized = finalize_archives(
                        root,
                        mode,
                        family_ids,
                        nonblocking=args.nonblocking,
                        compression=args.compression,
                        compression_level=args.compression_level,
                        max_final_zip_bytes=args.max_final_zip_bytes,
                        progress_callback=reporter.update,
                    )
                finally:
                    reporter.close()
                for zip_path in finalized:
                    print(f"finalized\t{zip_path}")
        return 0
    if args.command == "archive-family":
        _validate_zip_write_options(args)
        if args.max_final_zip_bytes < 0:
            raise ValueError("--max-final-zip-bytes must be non-negative")
        mode = _resolve_cli_mode(args, root)
        assert mode is not None
        if args.query_dir is not None or args.genecount is not None:
            family_ids, family_from_name = family_context(
                mode,
                query_dir=args.query_dir,
                genecount=args.genecount,
            )
        elif mode == "query2family":
            family_ids = [args.family_id]

            def family_from_name(name: str) -> Optional[str]:
                return query_id_from_name(name, [args.family_id])
        else:
            family_ids = [args.family_id]
            family_from_name = orthogroup_id_from_name
        if args.family_id not in set(family_ids):
            raise ValueError(f"Family ID is absent from the current input catalog: {args.family_id}")
        results = archive_completed_outputs(
            root=root,
            mode=mode,
            family_ids=[args.family_id],
            family_from_name=family_from_name,
            min_files=1,
            max_files_per_shard=max(1, args.max_files_per_shard),
            nonblocking=args.nonblocking,
            include_incomplete=True,
            compression=args.compression,
            compression_level=args.compression_level,
            workers=args.workers,
            catalog_sources=_catalog_sources(args),
            catalog_family_ids=family_ids,
            preserve_existing_catalog=(args.query_dir is None and args.genecount is None),
            max_bytes_per_shard=args.max_final_zip_bytes,
        )
        for zip_path, removed in results:
            print(f"archived-partial\t{removed}\t{zip_path}")
        _write_archive_status(root, nonblocking=args.nonblocking)
        return 0
    if args.command in {"convert-storage", "storage-status"}:
        if args.command == "convert-storage":
            _validate_zip_write_options(args)
            if args.progress_interval < 0:
                raise ValueError("--progress-interval must be non-negative")
            if args.available_bytes is not None and args.available_bytes < 0:
                raise ValueError("--available-bytes must be non-negative")
            if args.large_zip_warning_bytes < 0:
                raise ValueError("--large-zip-warning-bytes must be non-negative")
            if args.max_final_zip_bytes < 0:
                raise ValueError("--max-final-zip-bytes must be non-negative")
            if args.resume and args.dry_run:
                raise ValueError("--resume cannot be combined with --dry-run")
        target = "raw" if args.command == "convert-storage" and args.to == "files" else getattr(args, "to", None)
        if args.command == "convert-storage":
            if target == "zip" and args.pure_raw:
                raise ValueError("--pure-raw is valid only with --to raw or files")
            if target == "raw" and args.strict_unmatched:
                raise ValueError("--strict-unmatched is valid only with --to zip")
        mode = _resolve_cli_mode(
            args,
            root,
            required=(args.command == "convert-storage"),
        )
        needs_catalog = (
            target == "zip"
            or args.query_dir is not None
            or args.genecount is not None
            or args.family_id_file is not None
        )
        if needs_catalog:
            assert mode is not None
            family_ids, family_from_name = family_context_with_supplement(
                mode,
                query_dir=args.query_dir,
                genecount=args.genecount,
                family_id_file=args.family_id_file,
            )
        else:
            family_ids = []
            if mode is not None and _physical_archive_paths(root):
                archived_store = GeneFamilyOutputStore(root)
                archived_store._load_archives()
                family_ids = sorted(
                    {
                        artifact.family_id
                        for artifact in archived_store._archived.values()
                        if artifact.family_id is not None
                    }
                )
            if mode == "orthogroup":
                family_from_name = orthogroup_id_from_name
            elif mode == "query2family":
                archived_matchers = query_id_matchers(family_ids)

                def family_from_name(name: str) -> Optional[str]:
                    return query_id_from_name(name, archived_matchers)
            else:

                def family_from_name(name: str) -> Optional[str]:
                    return None

        if args.command == "storage-status" or args.dry_run:
            dry_run_reporter = None
            if args.command == "convert-storage" and args.dry_run:
                dry_run_reporter = ProgressReporter(args.progress_interval)
                dry_run_reporter.update(phase="inventory")
                dry_run_reporter.start()
            try:
                summary, unmatched = storage_conversion_summary(
                    root,
                    set(family_ids),
                    family_from_name,
                    workers=max(1, int(getattr(args, "workers", 1))),
                    large_zip_warning_bytes=int(
                        getattr(
                            args,
                            "large_zip_warning_bytes",
                            DEFAULT_LARGE_ZIP_WARNING_BYTES,
                        )
                    ),
                    max_final_zip_bytes=int(getattr(args, "max_final_zip_bytes", 0)),
                )
            finally:
                if dry_run_reporter is not None:
                    dry_run_reporter.update(force=True, phase="complete")
                    dry_run_reporter.close()
            if args.command == "convert-storage":
                summary["requested_target"] = target
                summary["pure_raw"] = bool(args.pure_raw)
                summary["mode"] = mode or "unknown"
                filesystem_free_bytes = int(shutil.disk_usage(root).free)
                effective_available_bytes = filesystem_free_bytes
                if args.available_bytes is not None:
                    effective_available_bytes = min(
                        effective_available_bytes,
                        int(args.available_bytes),
                    )
                required_peak_new_bytes = int(
                    summary["raw_zip_peak_new_bytes" if target == "zip" else "raw_materialize_allocated_bytes"]
                )
                summary["filesystem_free_bytes"] = filesystem_free_bytes
                summary["effective_available_bytes"] = effective_available_bytes
                summary["required_peak_new_bytes"] = required_peak_new_bytes
                summary["temporary_space_sufficient"] = int(required_peak_new_bytes <= effective_available_bytes)
            if getattr(args, "json", False):
                json_payload = dict(summary)
                json_payload["unmatched_paths"] = [str(path) for path in unmatched]
                if needs_catalog:
                    json_payload["conflicting_final_zip_paths"] = [
                        str(path)
                        for path in _catalog_conflicting_final_zip_paths(
                            root,
                            set(family_ids),
                            family_from_name,
                        )
                    ]
                print(json.dumps(json_payload, indent=2, sort_keys=True))
            else:
                for key, value in summary.items():
                    print(f"{key}\t{value}")
            if needs_catalog:
                if not getattr(args, "json", False):
                    for path in unmatched:
                        print(f"unmatched\t{path}")
                    for path in _catalog_conflicting_final_zip_paths(
                        root,
                        set(family_ids),
                        family_from_name,
                    ):
                        print(f"conflicting-final-zip\t{path}")
            if args.command == "storage-status" or args.dry_run:
                return 0
        if target == "zip":
            assert mode is not None
            reporter = ProgressReporter(args.progress_interval)
            reporter.start()
            try:
                summary = convert_storage_to_zip(
                    root,
                    mode,
                    family_ids,
                    family_from_name,
                    max_files_per_shard=max(1, args.max_files_per_shard),
                    strict_unmatched=args.strict_unmatched,
                    compression=args.compression,
                    compression_level=args.compression_level,
                    workers=args.workers,
                    catalog_sources=_catalog_sources(args),
                    require_resume=args.resume,
                    available_bytes=args.available_bytes,
                    large_zip_warning_bytes=args.large_zip_warning_bytes,
                    max_final_zip_bytes=args.max_final_zip_bytes,
                    progress_callback=reporter.update,
                )
            finally:
                reporter.close()
            _, unmatched_after = scan_live_output_ownership(
                root,
                set(family_ids),
                family_from_name,
                workers=max(1, args.workers),
            )
        else:
            assert mode is not None
            reporter = ProgressReporter(args.progress_interval)
            reporter.start()
            try:
                summary = convert_storage_to_raw(
                    root,
                    mode,
                    max_files_per_shard=max(1, args.max_files_per_shard),
                    pure_raw=args.pure_raw,
                    require_resume=args.resume,
                    available_bytes=args.available_bytes,
                    progress_callback=reporter.update,
                )
            finally:
                reporter.close()
        if args.json:
            json_payload = dict(summary)
            if target == "zip":
                json_payload["unmatched_paths"] = [str(path) for path in unmatched_after]
            print(json.dumps(json_payload, indent=2, sort_keys=True))
        else:
            for key, value in summary.items():
                print(f"{key}\t{value}")
            if target == "zip":
                for path in unmatched_after:
                    print(f"unmatched\t{path}")
        return 0
    if args.command == "cleanup-tmp":
        for removed_path in cleanup_stale_tmp(
            root,
            args.older_than_days,
            max_directories=args.max_directories,
            max_bytes=args.max_bytes,
            max_files=args.max_files,
            nonblocking=args.nonblocking,
        ):
            print(f"removed-tmp\t{removed_path}")
        return 0
    if args.command == "finalize":
        _validate_zip_write_options(args)
        if args.progress_interval < 0:
            raise ValueError("--progress-interval must be non-negative")
        if args.max_final_zip_bytes < 0:
            raise ValueError("--max-final-zip-bytes must be non-negative")
        mode = _resolve_cli_mode(args, root)
        assert mode is not None
        family_ids, _ = family_context(
            mode,
            query_dir=args.query_dir,
            genecount=args.genecount,
        )
        reporter = ProgressReporter(args.progress_interval)
        reporter.start()
        try:
            finalized = finalize_archives(
                root,
                mode,
                family_ids,
                subdirs=_parse_subdirs(args.subdirs),
                nonblocking=args.nonblocking,
                compression=args.compression,
                compression_level=args.compression_level,
                max_final_zip_bytes=args.max_final_zip_bytes,
                progress_callback=reporter.update,
            )
        finally:
            reporter.close()
        for zip_path in finalized:
            print(f"finalized\t{zip_path}")
        return 0
    if args.command == "repair":
        if args.progress_interval < 0:
            raise ValueError("--progress-interval must be non-negative")
        reporter = ProgressReporter(args.progress_interval)
        reporter.start()
        try:
            orphaned = repair_archive_index(
                root,
                remove_orphans=args.remove_orphans,
                progress_callback=reporter.update,
            )
        finally:
            reporter.close()
        for path in orphaned:
            action = "removed-orphan" if args.remove_orphans else "orphan"
            print(f"{action}\t{path}")
        return 0
    if args.command == "lock-path":
        print(family_lock_path(_archive_state_root(root), args.family_id))
        return 0

    family_filter = (
        args.family_id
        if args.command
        in {
            "materialize-family",
            "is-complete",
            "mark-running",
            "mark-complete",
            "mark-failed",
        }
        else None
    )
    store = GeneFamilyOutputStore(root, family_filter=family_filter)
    if args.command == "materialize-family":
        mode = _resolve_cli_mode(args, root)
        if mode == "query2family":
            if args.query_dir is not None:
                _, family_from_name = family_context(
                    "query2family",
                    query_dir=args.query_dir,
                )
            else:

                def family_from_name(name: str) -> Optional[str]:
                    return query_id_from_name(name, [args.family_id])
        else:
            family_from_name = orthogroup_id_from_name
        restored = store.materialize_family(
            args.family_id,
            family_from_name,
            destination_root=args.destination_root,
            subdirs=_parse_subdirs(args.subdirs),
            receipt_path=getattr(args, "receipt", None),
            run_token=getattr(args, "run_token", ""),
        )
        for path in restored:
            print(path)
        return 0
    if args.command == "materialize-families":
        mode = _resolve_cli_mode(args, root)
        with args.family_id_file.open("r", encoding="utf-8") as handle:
            selected_family_ids = {family_id for line in handle if (family_id := line.rstrip("\r\n"))}
        if mode == "query2family":
            matchers = query_id_matchers(selected_family_ids)

            def family_from_name(name: str) -> Optional[str]:
                return query_id_from_name(name, matchers)
        else:
            family_from_name = orthogroup_id_from_name
        for path in store.materialize_families(
            selected_family_ids,
            family_from_name,
            destination_root=args.destination_root,
            subdirs=_parse_subdirs(args.subdirs),
        ):
            print(path)
        return 0
    if args.command == "export-current":
        if args.destination_root.is_symlink():
            raise ValueError("--destination-root must not be a symlink")
        destination_root = args.destination_root.resolve()
        if destination_root == root or root in destination_root.parents:
            raise ValueError("--destination-root must be outside the logical store root")
        if destination_root.exists():
            if not destination_root.is_dir():
                raise ValueError("--destination-root must be a directory")
            if any(destination_root.iterdir()):
                raise ValueError("--destination-root must be absent or empty")
        for path in store.materialize_all(
            destination_root,
        ):
            print(path)
        return 0
    if args.command == "is-complete":
        state = store.family_state(args.family_id)
        if state == "complete":
            return 0
        if state is None and completion_outputs_present(store, args.family_id):
            return 0
        return 1
    if args.command == "delete":
        store.delete(
            args.path,
            remove_live=not args.keep_live,
            family_id=getattr(args, "family_id", None),
        )
        return 0
    if args.command == "undelete":
        store.undelete(
            args.path,
            family_id=getattr(args, "family_id", None),
        )
        return 0
    if args.command == "restore":
        print(
            store.restore(
                args.path,
                overwrite=args.overwrite,
                family_id=getattr(args, "family_id", None),
            )
        )
        return 0
    if args.command == "verify":
        if args.progress_interval < 0:
            raise ValueError("--progress-interval must be non-negative")
        reporter = ProgressReporter(args.progress_interval)
        reporter.start()
        try:
            verification_mode = "quick" if args.quick else "deep"
            verified_paths = store.verify(
                progress_callback=reporter.update,
                deep=(verification_mode == "deep"),
            )
        finally:
            reporter.close()
        if args.json:
            print(
                json.dumps(
                    {
                        "verification_mode": verification_mode,
                        "verified_zip_files": [str(path) for path in verified_paths],
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
        else:
            print(f"verification_mode\t{verification_mode}")
            for zip_path in verified_paths:
                print(f"verified\t{zip_path}")
        return 0
    if args.command == "list":
        subdirs = [args.subdir] if args.subdir else store.logical_subdirs()
        for subdir in subdirs:
            for artifact in store.artifacts(subdir):
                source = "live" if artifact.is_live else str(artifact.zip_path)
                print(f"{artifact.logical_path}\t{artifact.size or 0}\t{source}")
        return 0
    if args.command == "has-files":
        return (
            0
            if any(
                artifact.name.endswith(args.suffix) and (artifact.size or 0) > 0
                for artifact in store.artifacts(args.subdir)
            )
            else 1
        )
    if args.command == "status":
        mode = _resolve_cli_mode(args, root)
        family_ids, family_from_name = family_context(
            mode,
            query_dir=args.query_dir,
            genecount=args.genecount,
        )
        valid_ids = set(family_ids)
        for subdir in store.logical_subdirs():
            found_ids = {
                family_id
                for family_id in (family_from_name(name) for name in store.file_names(subdir))
                if family_id is not None and family_id in valid_ids
            }
            print(
                f"{subdir}\tpresent={len(found_ids)}\tmissing={len(valid_ids - found_ids)}\texpected={len(valid_ids)}"
            )
        return 0
    if args.command == "compact":
        _validate_zip_write_options(args)
        if args.max_final_zip_bytes < 0:
            raise ValueError("--max-final-zip-bytes must be non-negative")
        mode = _resolve_cli_mode(args, root)
        for zip_path in compact_archives(
            root,
            mode,
            max_files_per_shard=max(1, args.max_files_per_shard),
            nonblocking=args.nonblocking,
            compression=args.compression,
            compression_level=args.compression_level,
            max_bytes_per_shard=args.max_final_zip_bytes,
        ):
            print(f"compacted\t{zip_path}")
        _write_archive_status(root)
        return 0
    if args.command == "purge":
        _validate_zip_write_options(args)
        mode = _resolve_cli_mode(args, root)
        valid_family_ids: Optional[Set[str]] = None
        family_from_name: Optional[Callable[[str], Optional[str]]] = None
        if args.drop_unlisted:
            family_ids, family_from_name = family_context(
                mode,
                query_dir=args.query_dir,
                genecount=args.genecount,
            )
            valid_family_ids = set(family_ids)
        for zip_path in purge_archives(
            root,
            mode,
            max_files_per_shard=max(1, args.max_files_per_shard),
            valid_family_ids=valid_family_ids,
            family_from_name=family_from_name,
            compression=args.compression,
            compression_level=args.compression_level,
        ):
            print(f"purged\t{zip_path}")
        _write_archive_status(root)
        return 0
    if args.command in {"mark-running", "mark-complete", "mark-failed"}:
        store.mark_family_state(
            args.family_id,
            args.family_status,
            args.run_token,
        )
        return 0
    raise AssertionError(args.command)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    try:
        exit_code = run_cli(args)
    except (ArchiveStoreError, OSError, ValueError) as exc:
        print(f"gene_family_output_store: {exc}", file=sys.stderr)
        exit_code = 2
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
