#!/usr/bin/env python3
# coding: utf-8

import argparse
import hashlib
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import threading
import time
import uuid
import zipfile
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from functools import lru_cache
from importlib import metadata
from itertools import repeat
from multiprocessing import get_context
from pathlib import Path, PurePosixPath

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import csubst_site_wrapper as site_wrapper  # noqa: E402, I001


DEFAULT_MIN_SUPPORT = 5
DEFAULT_Q_COLUMN = "q_rate_enrichment_global"
DEFAULT_Q_THRESHOLD = 0.05
RUN_LOCK_HEARTBEAT_SECONDS = 5
RUN_LOCK_STALE_SECONDS = 300
BESTHIT_COLUMNS = [
    "besthit_0.05",
    "besthit_0.25",
    "besthit_0.5",
    "besthit_0.75",
    "besthit_0.95",
]
CANDIDATE_REQUIRED_COLUMNS = [
    "orthogroup",
    "trait",
    "state_change",
    "codon_site_alignment",
    "support_unit_count",
    "support_unit_ids",
    "support_branch_ids",
]
INTERNAL_COLUMNS = {
    "_analysis_key",
    "_cache_name",
    "_candidate_id",
    "_canonical_support_branch_ids",
    "_candidate_rank",
    "_source_summary_tsv",
    "_selection_min_support",
}
ARCHIVE_MANIFEST_COLUMNS = [
    "min_support",
    "candidate_count",
    "summary_tsv",
    "summary_sha256",
    "analysis_engine_signature",
    "csubst_version",
    "csubst_signature",
    "runtime_dependency_versions",
    "q_column",
    "q_threshold",
    "archive_zip",
    "status",
    "error",
]
CANDIDATE_MANIFEST_COLUMNS = [
    "candidate_rank",
    "candidate_id",
    "orthogroup",
    "trait",
    "codon_site_alignment",
    "state_change",
    "support_unit_count",
    "support_unit_ids",
    "support_branch_ids",
    "q_column",
    "q_value",
    "candidate_tsv",
    "focused_tree_site_pdf",
    "report_pdf",
    "csubst_sites_dir",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run CSUBST sites for significant CSUBST scan candidates once, then "
            "package self-contained candidate reports into descending min_support ZIP files."
        )
    )
    parser.add_argument(
        "--summary_prefix",
        metavar="PATH_PREFIX",
        required=True,
        help=(
            "Prefix before _min_support_<N>_summary.tsv, for example "
            "gene_summary/orthogroup/orthogroup_csubst_aa_change."
        ),
    )
    parser.add_argument("--dir_orthogroup", metavar="PATH", required=True)
    parser.add_argument("--file_trait", metavar="PATH", required=True)
    parser.add_argument("--out_dir", metavar="PATH", required=True)
    parser.add_argument("--min_support", metavar="INT", default=DEFAULT_MIN_SUPPORT, type=int)
    parser.add_argument("--q_column", metavar="COLUMN", default=DEFAULT_Q_COLUMN)
    parser.add_argument("--q_threshold", metavar="FLOAT", default=DEFAULT_Q_THRESHOLD, type=float)
    parser.add_argument(
        "--max_candidates",
        metavar="INT",
        default=0,
        type=int,
        help="Maximum candidates per min_support ZIP. Zero retains all significant candidates.",
    )
    parser.add_argument("--ncpu", metavar="INT", default=1, type=int)
    parser.add_argument(
        "--csubst_nonsyn_recode",
        metavar="STR",
        default="no",
        type=site_wrapper.normalize_csubst_nonsyn_recode,
    )
    parser.add_argument(
        "--pdb",
        choices=("none", "besthit"),
        default="none",
        help="Optional protein-structure search for each unique candidate analysis.",
    )
    return parser.parse_args()


def validate_args(args):
    if args.min_support < 2:
        raise ValueError("--min_support must be an integer >= 2.")
    if not np.isfinite(args.q_threshold) or not 0.0 <= args.q_threshold <= 1.0:
        raise ValueError("--q_threshold must be between 0 and 1.")
    if args.max_candidates < 0:
        raise ValueError("--max_candidates must be >= 0.")
    if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_.]*", str(args.q_column)) is None:
        raise ValueError(f"Unsafe q-value column name: {args.q_column}")
    args.ncpu = max(1, int(args.ncpu))
    dir_orthogroup = Path(args.dir_orthogroup).resolve()
    file_trait = Path(args.file_trait).resolve()
    if not dir_orthogroup.is_dir():
        raise FileNotFoundError(
            f"--dir_orthogroup directory was not found: {dir_orthogroup}"
        )
    if not file_trait.is_file():
        raise FileNotFoundError(f"--file_trait file was not found: {file_trait}")
    args.summary_prefix = str(Path(args.summary_prefix).resolve())
    args.dir_orthogroup = str(dir_orthogroup)
    args.file_trait = str(file_trait)
    args.out_dir = str(Path(args.out_dir).resolve())
    return args


def format_float_token(value):
    return format(float(value), ".12g")


def sanitize_token(value, fallback="value", max_length=80):
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_.-")
    if token == "":
        token = fallback
    return token[:max_length]


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@lru_cache(maxsize=1)
def csubst_version():
    try:
        return metadata.version("csubst")
    except metadata.PackageNotFoundError:
        return "unavailable"


@lru_cache(maxsize=1)
def csubst_signature():
    try:
        distribution = metadata.distribution("csubst")
    except metadata.PackageNotFoundError:
        return "unavailable"
    digest = hashlib.sha256()
    distribution_files = distribution.files or []
    source_files = sorted(
        path
        for path in distribution_files
        if path.parts
        and path.parts[0] == "csubst"
        and path.suffix in {".py", ".so"}
        and "__pycache__" not in path.parts
    )
    for relative_path in source_files:
        installed_path = Path(distribution.locate_file(relative_path))
        digest.update(str(relative_path).encode("utf-8"))
        if installed_path.is_file():
            digest.update(installed_path.read_bytes())
        else:
            digest.update(b"missing")
    return digest.hexdigest()


@lru_cache(maxsize=1)
def runtime_dependency_versions():
    versions = {}
    for package_name in ("numpy", "pandas", "pypdf", "reportlab"):
        try:
            versions[package_name] = metadata.version(package_name)
        except metadata.PackageNotFoundError:
            versions[package_name] = "unavailable"
    return json.dumps(versions, sort_keys=True, separators=(",", ":"))


@lru_cache(maxsize=1)
def analysis_engine_signature():
    digest = hashlib.sha256()
    engine_paths = [
        Path(__file__).resolve(),
        SCRIPT_DIR / "csubst_site_wrapper.py",
        SCRIPT_DIR / "stat_branch2tree_plot.r",
    ]
    treevis_dir = SCRIPT_DIR / "treevis"
    if treevis_dir.is_dir():
        engine_paths.extend(
            path
            for path in treevis_dir.rglob("*")
            if path.is_file()
            and (path.suffix.lower() == ".r" or path.name in {"DESCRIPTION", "NAMESPACE"})
        )
    for path in sorted(engine_paths):
        try:
            path_label = str(path.relative_to(SCRIPT_DIR))
        except ValueError:
            path_label = path.name
        digest.update(path_label.encode("utf-8"))
        if path.is_file():
            digest.update(path.read_bytes())
        else:
            digest.update(b"missing")
    digest.update(csubst_version().encode("utf-8"))
    digest.update(csubst_signature().encode("utf-8"))
    digest.update(runtime_dependency_versions().encode("utf-8"))
    return digest.hexdigest()


def discover_summary_tables(summary_prefix, minimum_support):
    prefix = Path(summary_prefix)
    pattern = re.compile(
        rf"^{re.escape(prefix.name)}_min_support_([0-9]+)_summary\.tsv$"
    )
    all_discovered = {}
    for candidate in prefix.parent.glob(f"{prefix.name}_min_support_*_summary.tsv"):
        match = pattern.match(candidate.name)
        if match is None:
            continue
        threshold = int(match.group(1))
        all_discovered[threshold] = candidate
    if not all_discovered:
        raise FileNotFoundError(
            f"No min_support summary TSVs were found for prefix: {prefix}"
        )
    discovered = {
        threshold: path
        for threshold, path in all_discovered.items()
        if threshold >= minimum_support
    }
    if not discovered and max(all_discovered) < minimum_support:
        return {}
    maximum_support = max(discovered)
    expected = set(range(minimum_support, maximum_support + 1))
    missing = sorted(expected.difference(discovered))
    if missing:
        raise FileNotFoundError(
            "The min_support summary series is incomplete. Missing threshold(s): "
            + ", ".join(str(value) for value in missing)
        )
    return {threshold: discovered[threshold] for threshold in sorted(discovered, reverse=True)}


def parse_branch_ids(value):
    if pd.isna(value):
        raise ValueError("support_branch_ids is blank.")
    tokens = [token.strip() for token in str(value).split(",") if token.strip()]
    if not tokens:
        raise ValueError("support_branch_ids is blank.")
    branch_ids = []
    for token in tokens:
        if re.fullmatch(r"[0-9]+", token) is None:
            raise ValueError(f"Invalid branch ID in support_branch_ids: {token}")
        branch_ids.append(int(token))
    return sorted(set(branch_ids))


def integer_value(value, column, minimum=None):
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if not np.isfinite(numeric) or float(numeric) != int(numeric):
        raise ValueError(f"{column} must contain finite integers: {value}")
    numeric = int(numeric)
    if minimum is not None and numeric < minimum:
        raise ValueError(f"{column} must be >= {minimum}: {value}")
    return numeric


def candidate_identity(row):
    identity_columns = [
        "orthogroup",
        "trait",
        "target_class",
        "scan_match",
        "codon_site_alignment",
        "from_state",
        "to_state",
        "state_change",
        "_canonical_support_branch_ids",
    ]
    return {
        column: None if column not in row or pd.isna(row[column]) else str(row[column])
        for column in identity_columns
    }


def assign_candidate_ids(frame, csubst_nonsyn_recode, pdb):
    candidate_ids = []
    analysis_keys = []
    cache_names = []
    for _, row in frame.iterrows():
        identity = candidate_identity(row)
        identity_text = json.dumps(identity, sort_keys=True, separators=(",", ":"))
        identity_digest = hashlib.sha256(identity_text.encode("utf-8")).hexdigest()[:16]
        candidate_id = "{}_site{}_{}_{}".format(
            sanitize_token(row["orthogroup"], fallback="orthogroup", max_length=40),
            int(row["codon_site_alignment"]),
            sanitize_token(row["state_change"], fallback="change", max_length=24),
            identity_digest,
        )
        analysis_payload = {
            "candidate": identity,
            "csubst_nonsyn_recode": csubst_nonsyn_recode,
            "pdb": pdb,
            "analysis_engine_signature": analysis_engine_signature(),
        }
        analysis_text = json.dumps(analysis_payload, sort_keys=True, separators=(",", ":"))
        analysis_key = hashlib.sha256(analysis_text.encode("utf-8")).hexdigest()
        candidate_ids.append(candidate_id)
        analysis_keys.append(analysis_key)
        cache_names.append(f"{candidate_id}_{analysis_key[:16]}")
    frame = frame.copy()
    frame["_candidate_id"] = candidate_ids
    frame["_analysis_key"] = analysis_keys
    frame["_cache_name"] = cache_names
    if frame["_candidate_id"].duplicated().any():
        duplicated = frame.loc[frame["_candidate_id"].duplicated(keep=False), "_candidate_id"]
        raise ValueError(
            "Candidate IDs are not unique within a min_support summary: "
            + ", ".join(duplicated.drop_duplicates().head(5).tolist())
        )
    return frame


def load_threshold_candidates(
    summary_path,
    minimum_support,
    q_column,
    q_threshold,
    max_candidates,
    csubst_nonsyn_recode,
    pdb,
):
    frame = pd.read_csv(summary_path, sep="\t", low_memory=False)
    required = [*CANDIDATE_REQUIRED_COLUMNS, q_column]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(
            f"{summary_path} is missing required candidate column(s): {', '.join(missing)}"
        )
    support = pd.to_numeric(frame["support_unit_count"], errors="coerce")
    invalid_support = ~np.isfinite(support.to_numpy(dtype=float)) | (support < minimum_support)
    if invalid_support.any():
        raise ValueError(
            f"{summary_path} contains {int(invalid_support.sum()):,} row(s) below min_support={minimum_support}."
        )
    qvalues = pd.to_numeric(frame[q_column], errors="coerce")
    keep = np.isfinite(qvalues.to_numpy(dtype=float)) & (qvalues <= q_threshold)
    selected = frame.loc[keep, :].copy()
    selected[q_column] = qvalues.loc[keep].astype(float)
    selected["support_unit_count"] = support.loc[keep].astype(int)
    selected["codon_site_alignment"] = [
        integer_value(value, "codon_site_alignment", minimum=1)
        for value in selected["codon_site_alignment"]
    ]
    canonical_branch_ids = []
    for row_index, value in selected["support_branch_ids"].items():
        try:
            branch_ids = parse_branch_ids(value)
        except ValueError as error:
            raise ValueError(f"{summary_path}, row {row_index + 2}: {error}") from error
        canonical_branch_ids.append(",".join(str(branch_id) for branch_id in branch_ids))
    selected["_canonical_support_branch_ids"] = canonical_branch_ids
    selected["orthogroup"] = selected["orthogroup"].astype("string").str.strip()
    selected["trait"] = selected["trait"].astype("string").str.strip()
    if selected["orthogroup"].isna().any() or selected["orthogroup"].eq("").any():
        raise ValueError(f"{summary_path} contains blank orthogroup values among selected candidates.")
    if selected["trait"].isna().any() or selected["trait"].eq("").any():
        raise ValueError(f"{summary_path} contains blank trait values among selected candidates.")
    sort_columns = [q_column]
    ascending = [True]
    if "p_rate_enrichment" in selected.columns:
        selected["p_rate_enrichment"] = pd.to_numeric(
            selected["p_rate_enrichment"], errors="coerce"
        )
        sort_columns.append("p_rate_enrichment")
        ascending.append(True)
    sort_columns.extend(["support_unit_count", "orthogroup", "codon_site_alignment", "state_change"])
    ascending.extend([False, True, True, True])
    selected = selected.sort_values(
        sort_columns,
        ascending=ascending,
        na_position="last",
        kind="mergesort",
    ).reset_index(drop=True)
    if max_candidates > 0:
        selected = selected.head(max_candidates).copy()
    selected = assign_candidate_ids(selected, csubst_nonsyn_recode, pdb)
    selected["_candidate_rank"] = np.arange(1, selected.shape[0] + 1, dtype=int)
    selected["_source_summary_tsv"] = Path(summary_path).name
    selected["_selection_min_support"] = int(minimum_support)
    return selected


def write_trait_color_tables(file_trait, traits, output_dir):
    trait_frame = pd.read_csv(file_trait, sep="\t", low_memory=False)
    if "species" not in trait_frame.columns:
        raise ValueError(f"Trait table is missing the species column: {file_trait}")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for trait in sorted(set(str(value) for value in traits)):
        if trait not in trait_frame.columns:
            raise ValueError(f"Trait column was not found in {file_trait}: {trait}")
        colors = pd.DataFrame(
            {
                "species": trait_frame["species"],
                "color": np.where(
                    trait_frame[trait].map(site_wrapper.is_foreground_trait_value),
                    "firebrick",
                    "black",
                ),
            }
        )
        trait_digest = hashlib.sha256(trait.encode("utf-8")).hexdigest()[:8]
        path = output_dir / (
            f"trait_{sanitize_token(trait, max_length=70)}_{trait_digest}.color.tsv"
        )
        colors.to_csv(path, sep="\t", index=False)
        paths[trait] = str(path.resolve())
    return paths


def safe_extract_zip(zip_path, destination):
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "r") as archive:
        for member in archive.infolist():
            target = (destination / member.filename).resolve()
            if os.path.commonpath([str(destination), str(target)]) != str(destination):
                raise ValueError(f"Unsafe path in ZIP archive {zip_path}: {member.filename}")
        archive.extractall(destination)


def candidate_cache_complete(cache_dir, record):
    cache_dir = Path(cache_dir)
    marker = cache_dir / "analysis.complete.tsv"
    focused = cache_dir / f"{record['_candidate_id']}.focused_tree_site.pdf"
    if not marker.is_file() or not focused.is_file() or focused.stat().st_size == 0:
        return False
    try:
        marker_frame = pd.read_csv(marker, sep="\t", dtype=str)
        if marker_frame.shape[0] != 1:
            return False
        if marker_frame.loc[0, "analysis_key"] != record["_analysis_key"]:
            return False
        artifacts = site_wrapper.resolve_site_artifacts(
            str(cache_dir), record["_canonical_support_branch_ids"]
        )
        for key in ("site_table_tsv", "site_summary_pdf"):
            path = artifacts[key]
            if path is None or not os.path.isfile(path) or os.path.getsize(path) == 0:
                return False
    except Exception:
        return False
    return True


def analyze_candidate(record, cache_root, effective_dir_orthogroup, trait_color_paths, nonsyn_recode, pdb):
    cache_root = Path(cache_root)
    cache_dir = cache_root / record["_cache_name"]
    if candidate_cache_complete(cache_dir, record):
        return {"candidate_id": record["_candidate_id"], "status": "cached", "error": ""}
    if cache_dir.exists():
        shutil.rmtree(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    previous_cwd = os.getcwd()
    try:
        os.chdir(cache_dir)
        iqtree_zip = site_wrapper.get_iqtree_anc_zip_path(
            effective_dir_orthogroup, record["orthogroup"]
        )
        if not os.path.isfile(iqtree_zip):
            raise FileNotFoundError(f"IQ-TREE ancestral-state ZIP was not found: {iqtree_zip}")
        safe_extract_zip(iqtree_zip, cache_dir)
        iqtree_anc_dir = site_wrapper.get_iqtree_anc_dir(
            str(cache_dir), record["orthogroup"]
        )
        iqtree_anc_rel_dir = os.path.basename(iqtree_anc_dir)
        command = site_wrapper.build_csubst_sites_command(
            iqtree_anc_rel_dir=iqtree_anc_rel_dir,
            iqtree_anc_dir=iqtree_anc_dir,
            branch_id_str=record["_canonical_support_branch_ids"],
            ncpu=1,
            csubst_nonsyn_recode=nonsyn_recode,
            pdb=pdb,
        )
        print("COMMAND: {}".format(site_wrapper.shell_join_command(command)), flush=True)
        subprocess.run(command, check=True)
        artifacts = site_wrapper.resolve_site_artifacts(
            str(cache_dir), record["_canonical_support_branch_ids"]
        )
        if artifacts["site_table_tsv"] is None:
            raise FileNotFoundError(f"CSUBST sites table was not found for {record['_candidate_id']}")
        site_table = pd.read_csv(artifacts["site_table_tsv"], sep="\t", low_memory=False)
        if "codon_site_alignment" not in site_table.columns:
            raise ValueError(
                f"CSUBST sites table lacks codon_site_alignment: {artifacts['site_table_tsv']}"
            )
        observed_sites = pd.to_numeric(site_table["codon_site_alignment"], errors="coerce")
        target_site = int(record["codon_site_alignment"])
        if not observed_sites.eq(target_site).any():
            raise ValueError(
                f"Candidate alignment site {target_site} was not found in CSUBST sites output "
                f"for {record['_candidate_id']}."
            )
        focused_name = f"{record['_candidate_id']}.focused_tree_site.pdf"
        site_wrapper.run_stat_branch2tree_plot(
            og=record["orthogroup"],
            branch_id_str=record["_canonical_support_branch_ids"],
            file_trait_color=trait_color_paths[record["trait"]],
            dir_out_og=str(cache_dir),
            dir_og=effective_dir_orthogroup,
            ncpu=1,
            csubst_nonsyn_recode=nonsyn_recode,
            convergent_sites=[target_site],
            file_tree_plot_out=str(cache_dir / focused_name),
        )
        shutil.rmtree(iqtree_anc_dir, ignore_errors=True)
        for temporary_name in (
            "tmp.csubst.sub_tensor.S.mmap",
            "tmp.csubst.sub_tensor.N.mmap",
            "Rplots.pdf",
        ):
            (cache_dir / temporary_name).unlink(missing_ok=True)
        marker = pd.DataFrame(
            [
                {
                    "analysis_key": record["_analysis_key"],
                    "candidate_id": record["_candidate_id"],
                    "orthogroup": record["orthogroup"],
                    "codon_site_alignment": target_site,
                    "support_branch_ids": record["_canonical_support_branch_ids"],
                    "csubst_nonsyn_recode": nonsyn_recode,
                    "pdb": pdb,
                    "analysis_engine_signature": analysis_engine_signature(),
                    "csubst_version": csubst_version(),
                    "csubst_signature": csubst_signature(),
                    "runtime_dependency_versions": runtime_dependency_versions(),
                    "command": site_wrapper.shell_join_command(command),
                }
            ]
        )
        marker.to_csv(cache_dir / "analysis.complete.tsv", sep="\t", index=False)
        return {"candidate_id": record["_candidate_id"], "status": "completed", "error": ""}
    except Exception as error:
        print(f"Candidate analysis failed for {record['_candidate_id']}: {error}", flush=True)
        return {
            "candidate_id": record["_candidate_id"],
            "status": "failed",
            "error": str(error),
        }
    finally:
        os.chdir(previous_cwd)


def analyze_orthogroup_batch(
    orthogroup,
    records,
    cache_root,
    dir_orthogroup,
    materialization_parent,
    trait_color_paths,
    nonsyn_recode,
    pdb,
):
    materialized_context = None
    effective_dir_orthogroup = dir_orthogroup
    try:
        if any(
            os.path.isdir(os.path.join(dir_orthogroup, name))
            for name in (".gg_store", ".gg_archives")
        ):
            materialized_context = site_wrapper.LockedMaterializationDirectory(
                materialization_parent, orthogroup
            )
            effective_dir_orthogroup = materialized_context.name
            site_wrapper.materialize_csubst_site_inputs(
                dir_og=dir_orthogroup,
                og=orthogroup,
                destination_root=effective_dir_orthogroup,
            )
        return [
            analyze_candidate(
                record=record,
                cache_root=cache_root,
                effective_dir_orthogroup=effective_dir_orthogroup,
                trait_color_paths=trait_color_paths,
                nonsyn_recode=nonsyn_recode,
                pdb=pdb,
            )
            for record in records
        ]
    finally:
        if materialized_context is not None:
            materialized_context.cleanup()


def ensure_candidate_analyses(
    candidates,
    cache_root,
    dir_orthogroup,
    materialization_parent,
    trait_color_paths,
    nonsyn_recode,
    pdb,
    ncpu,
):
    pending_records = []
    seen_analysis_keys = set()
    for record in candidates.to_dict(orient="records"):
        if record["_analysis_key"] in seen_analysis_keys:
            continue
        seen_analysis_keys.add(record["_analysis_key"])
        cache_dir = Path(cache_root) / record["_cache_name"]
        if not candidate_cache_complete(cache_dir, record):
            pending_records.append(record)
    if not pending_records:
        return []
    records_by_orthogroup = {}
    for record in pending_records:
        records_by_orthogroup.setdefault(record["orthogroup"], []).append(record)
    orthogroups = list(records_by_orthogroup)
    worker_count = min(max(1, int(ncpu)), len(orthogroups))
    failures = []
    # Use spawn so analysis workers do not inherit the parent process's run-lock
    # descriptor.  Inheriting that descriptor under fork can make a concurrent
    # invocation enter the shared work directory before the first run finishes.
    with ProcessPoolExecutor(
        max_workers=worker_count,
        mp_context=get_context("spawn"),
    ) as executor:
        results = executor.map(
            analyze_orthogroup_batch,
            orthogroups,
            [records_by_orthogroup[orthogroup] for orthogroup in orthogroups],
            repeat(str(cache_root)),
            repeat(str(dir_orthogroup)),
            repeat(str(materialization_parent)),
            repeat(trait_color_paths),
            repeat(nonsyn_recode),
            repeat(pdb),
        )
        for batch_results in results:
            for result in batch_results:
                if result["status"] == "failed":
                    failures.append(result)
    return failures


def printable_value(value):
    if pd.isna(value):
        return "NA"
    return str(value)


def candidate_annotation_text(row, q_column, q_threshold):
    lines = [
        "CSUBST scan candidate",
        "",
        f"Candidate ID: {row['_candidate_id']}",
        f"Orthogroup: {row['orthogroup']}",
        f"Trait: {row['trait']}",
        f"Alignment site: {int(row['codon_site_alignment'])}",
        f"State change: {printable_value(row['state_change'])}",
        f"Support unit count: {int(row['support_unit_count'])}",
        f"Support unit IDs: {printable_value(row['support_unit_ids'])}",
        f"Support branch IDs: {row['_canonical_support_branch_ids']}",
        "",
        "Selection",
        "",
        f"min_support: {int(row['_selection_min_support'])}",
        f"q column: {q_column}",
        f"q threshold: <= {format_float_token(q_threshold)}",
        f"q value: {printable_value(row[q_column])}",
        f"Analytical P value: {printable_value(row.get('p_rate_enrichment', np.nan))}",
        f"Source summary: {row['_source_summary_tsv']}",
        "",
        "Representative best hits",
        "",
    ]
    lines.extend(
        f"{column}: {printable_value(row.get(column, np.nan))}"
        for column in BESTHIT_COLUMNS
    )
    return "\n\n".join(lines)


def candidate_output_frame(row, q_column, q_threshold):
    source_columns = [column for column in row.index if column not in INTERNAL_COLUMNS]
    output = pd.DataFrame([row[source_columns].to_dict()])
    output.insert(0, "candidate_id", row["_candidate_id"])
    output.insert(0, "candidate_rank", int(row["_candidate_rank"]))
    output.insert(2, "selection_min_support", int(row["_selection_min_support"]))
    output.insert(3, "selection_q_column", q_column)
    output.insert(4, "selection_q_threshold", float(q_threshold))
    output["support_branch_ids"] = row["_canonical_support_branch_ids"]
    return output


def copy_candidate_cache(cache_dir, destination):
    shutil.copytree(cache_dir, destination)
    marker = destination / "analysis.complete.tsv"
    marker.unlink(missing_ok=True)
    make_csubst_manifests_portable(destination)


def make_csubst_manifests_portable(candidate_dir):
    candidate_dir = Path(candidate_dir).resolve()
    for manifest_path in candidate_dir.rglob("csubst.outputs.tsv"):
        manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
        if "output_file" not in manifest.columns:
            continue
        portable_paths = []
        file_exists = []
        file_sizes = []
        for output_file in manifest["output_file"].fillna(""):
            output_file = str(output_file).strip()
            output_path = (manifest_path.parent / output_file).resolve()
            is_internal = (
                output_file != ""
                and os.path.commonpath([str(candidate_dir), str(output_path)])
                == str(candidate_dir)
            )
            exists = is_internal and output_path.is_file()
            portable_paths.append(output_file if is_internal else "")
            file_exists.append("Y" if exists else "N")
            file_sizes.append(output_path.stat().st_size if exists else -1)
        manifest["output_path"] = portable_paths
        manifest["file_exists"] = file_exists
        manifest["file_size_bytes"] = file_sizes
        self_rows = manifest["output_file"].fillna("").map(
            lambda value: (manifest_path.parent / str(value)).resolve()
            == manifest_path.resolve()
        )
        manifest_text = ""
        for _ in range(10):
            manifest_text = manifest.to_csv(sep="\t", index=False, lineterminator="\n")
            manifest_size = len(manifest_text.encode("utf-8"))
            previous_sizes = pd.to_numeric(
                manifest.loc[self_rows, "file_size_bytes"], errors="coerce"
            )
            if previous_sizes.eq(manifest_size).all():
                break
            manifest.loc[self_rows, "file_size_bytes"] = manifest_size
        manifest_text = manifest.to_csv(sep="\t", index=False, lineterminator="\n")
        temporary_path = manifest_path.with_name(
            f".{manifest_path.name}.tmp-{os.getpid()}-{uuid.uuid4().hex}"
        )
        try:
            temporary_path.write_text(manifest_text, encoding="utf-8")
            os.replace(temporary_path, manifest_path)
        finally:
            temporary_path.unlink(missing_ok=True)


def package_candidate(row, package_root, cache_root, q_column, q_threshold):
    candidate_dir_name = f"candidate_{int(row['_candidate_rank']):04d}_{row['_candidate_id']}"
    candidate_dir = package_root / candidate_dir_name
    cache_dir = Path(cache_root) / row["_cache_name"]
    if not candidate_cache_complete(cache_dir, row.to_dict()):
        raise RuntimeError(f"Candidate cache is incomplete: {cache_dir}")
    copy_candidate_cache(cache_dir, candidate_dir)
    candidate_tsv = candidate_dir / "candidate.tsv"
    candidate_output_frame(row, q_column, q_threshold).to_csv(
        candidate_tsv, sep="\t", index=False
    )
    annotation_pdf = candidate_dir / "annotation.pdf"
    site_wrapper.create_pdf(
        candidate_annotation_text(row, q_column, q_threshold),
        str(annotation_pdf),
    )
    focused_pdf = candidate_dir / f"{row['_candidate_id']}.focused_tree_site.pdf"
    artifacts = site_wrapper.resolve_site_artifacts(
        str(candidate_dir), row["_canonical_support_branch_ids"]
    )
    report_pdf = candidate_dir / f"{row['_candidate_id']}.report.pdf"
    report_parts = [str(annotation_pdf), str(focused_pdf), artifacts["site_summary_pdf"]]
    if artifacts["pymol_summary_pdf"] is not None:
        report_parts.append(artifacts["pymol_summary_pdf"])
    site_wrapper.combine_pdfs(report_parts, str(report_pdf))
    annotation_pdf.unlink(missing_ok=True)
    return {
        "candidate_rank": int(row["_candidate_rank"]),
        "candidate_id": row["_candidate_id"],
        "orthogroup": row["orthogroup"],
        "trait": row["trait"],
        "codon_site_alignment": int(row["codon_site_alignment"]),
        "state_change": row["state_change"],
        "support_unit_count": int(row["support_unit_count"]),
        "support_unit_ids": row["support_unit_ids"],
        "support_branch_ids": row["_canonical_support_branch_ids"],
        "q_column": q_column,
        "q_value": row[q_column],
        "candidate_tsv": f"{candidate_dir_name}/{candidate_tsv.name}",
        "focused_tree_site_pdf": f"{candidate_dir_name}/{focused_pdf.name}",
        "report_pdf": f"{candidate_dir_name}/{report_pdf.name}",
        "csubst_sites_dir": f"{candidate_dir_name}/{site_wrapper.CSUBST_SITES_OUTDIR}",
    }


def write_package_readme(path, threshold, q_column, q_threshold, candidate_count, source_summary):
    text = "\n".join(
        [
            "CSUBST scan candidate sites",
            "",
            f"Source summary: {source_summary}",
            f"Selection: support_unit_count >= {threshold}",
            f"Selection: {q_column} <= {format_float_token(q_threshold)}",
            f"Selected candidates: {candidate_count}",
            "",
            "Each candidate directory contains the one-row source table, raw CSUBST sites outputs,",
            "a focused tree plot for the selected alignment site, and a combined PDF report.",
            "",
        ]
    )
    Path(path).write_text(text, encoding="utf-8")


def write_package_metadata(
    path,
    threshold,
    q_column,
    q_threshold,
    candidate_count,
    source_summary,
):
    pd.DataFrame(
        [
            {
                "min_support": int(threshold),
                "q_column": q_column,
                "q_threshold": float(q_threshold),
                "candidate_count": int(candidate_count),
                "source_summary_tsv": Path(source_summary).name,
                "source_summary_sha256": file_sha256(source_summary),
                "analysis_engine_signature": analysis_engine_signature(),
                "csubst_version": csubst_version(),
                "csubst_signature": csubst_signature(),
                "runtime_dependency_versions": runtime_dependency_versions(),
            }
        ]
    ).to_csv(path, sep="\t", index=False)


def archive_site_outputs_are_complete(archive, manifest_member):
    required_columns = {
        "output_kind",
        "output_file",
        "output_path",
        "file_exists",
        "file_size_bytes",
    }
    required_kinds = {"site_table_tsv", "site_summary_pdf", "output_manifest"}
    try:
        with archive.open(manifest_member) as handle:
            manifest = pd.read_csv(handle, sep="\t", dtype=str, keep_default_na=False)
    except Exception:
        return False
    if not required_columns.issubset(manifest.columns):
        return False
    manifest_parent = PurePosixPath(manifest_member).parent
    for output_kind in required_kinds:
        rows = manifest.loc[manifest["output_kind"] == output_kind, :]
        if rows.shape[0] != 1:
            return False
        row = rows.iloc[0]
        output_file = str(row["output_file"]).strip()
        output_path = PurePosixPath(output_file)
        if (
            row["file_exists"] != "Y"
            or output_file == ""
            or output_path.is_absolute()
            or ".." in output_path.parts
            or str(row["output_path"]).strip() != output_file
        ):
            return False
        archive_member = str(manifest_parent / output_path)
        try:
            member_info = archive.getinfo(archive_member)
            recorded_size = int(row["file_size_bytes"])
        except (KeyError, TypeError, ValueError):
            return False
        if member_info.file_size <= 0 or recorded_size != member_info.file_size:
            return False
    return True


def archive_matches_source(archive_path, source_summary):
    archive_path = Path(archive_path)
    if not archive_path.is_file() or not zipfile.is_zipfile(archive_path):
        return False
    metadata_member = f"{archive_path.stem}/package_metadata.tsv"
    candidate_manifest_member = f"{archive_path.stem}/candidate_manifest.tsv"
    try:
        with zipfile.ZipFile(archive_path, "r") as archive:
            if archive.testzip() is not None:
                return False
            member_names = set(archive.namelist())
            root_prefix = f"{archive_path.stem}/"
            if any(
                not name.startswith(root_prefix)
                or ".." in PurePosixPath(name).parts
                or PurePosixPath(name).is_absolute()
                for name in member_names
            ):
                return False
            required_top_level = {
                metadata_member,
                candidate_manifest_member,
                f"{archive_path.stem}/README.txt",
            }
            if not required_top_level.issubset(member_names):
                return False
            with archive.open(metadata_member) as handle:
                package_metadata = pd.read_csv(handle, sep="\t", dtype=str)
            with archive.open(candidate_manifest_member) as handle:
                candidate_manifest = pd.read_csv(handle, sep="\t", dtype=str)
            if candidate_manifest.columns.tolist() != CANDIDATE_MANIFEST_COLUMNS:
                return False
            if package_metadata.shape[0] != 1:
                return False
            if int(package_metadata.loc[0, "candidate_count"]) != candidate_manifest.shape[0]:
                return False
            for _, candidate in candidate_manifest.iterrows():
                for column in ("candidate_tsv", "focused_tree_site_pdf", "report_pdf"):
                    if root_prefix + candidate[column] not in member_names:
                        return False
                csubst_prefix = root_prefix + candidate["csubst_sites_dir"].rstrip("/") + "/"
                output_manifests = [
                    name
                    for name in member_names
                    if name.startswith(csubst_prefix) and name.endswith("/csubst.outputs.tsv")
                ]
                if len(output_manifests) != 1 or not archive_site_outputs_are_complete(
                    archive, output_manifests[0]
                ):
                    return False
        return (
            package_metadata.loc[0, "source_summary_sha256"] == file_sha256(source_summary)
            and package_metadata.loc[0, "analysis_engine_signature"]
            == analysis_engine_signature()
            and package_metadata.loc[0, "csubst_version"] == csubst_version()
            and package_metadata.loc[0, "csubst_signature"] == csubst_signature()
            and package_metadata.loc[0, "runtime_dependency_versions"]
            == runtime_dependency_versions()
        )
    except Exception:
        return False


def create_zip_atomic(package_root, archive_path):
    archive_path = Path(archive_path)
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_base = archive_path.parent / (
        f".{archive_path.stem}.tmp-{os.getpid()}-{uuid.uuid4().hex}"
    )
    temporary_zip = Path(str(temporary_base) + ".zip")
    try:
        shutil.make_archive(
            str(temporary_base),
            "zip",
            root_dir=str(package_root.parent),
            base_dir=package_root.name,
        )
        os.replace(temporary_zip, archive_path)
    finally:
        temporary_zip.unlink(missing_ok=True)


def package_threshold(
    candidates,
    threshold,
    source_summary,
    archive_path,
    packages_root,
    cache_root,
    q_column,
    q_threshold,
):
    package_root = Path(packages_root) / Path(archive_path).stem
    if package_root.exists():
        shutil.rmtree(package_root)
    package_root.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    try:
        for _, row in candidates.iterrows():
            manifest_rows.append(
                package_candidate(
                    row=row,
                    package_root=package_root,
                    cache_root=cache_root,
                    q_column=q_column,
                    q_threshold=q_threshold,
                )
            )
        pd.DataFrame(manifest_rows, columns=CANDIDATE_MANIFEST_COLUMNS).to_csv(
            package_root / "candidate_manifest.tsv", sep="\t", index=False
        )
        write_package_readme(
            package_root / "README.txt",
            threshold=threshold,
            q_column=q_column,
            q_threshold=q_threshold,
            candidate_count=candidates.shape[0],
            source_summary=Path(source_summary).name,
        )
        write_package_metadata(
            package_root / "package_metadata.tsv",
            threshold=threshold,
            q_column=q_column,
            q_threshold=q_threshold,
            candidate_count=candidates.shape[0],
            source_summary=source_summary,
        )
        create_zip_atomic(package_root, archive_path)
        if not archive_matches_source(archive_path, source_summary):
            raise RuntimeError(f"Candidate-site ZIP validation failed: {archive_path}")
    except Exception:
        raise
    else:
        shutil.rmtree(package_root)


def output_suffix(q_column, q_threshold, max_candidates, nonsyn_recode, pdb):
    suffix = f"{sanitize_token(q_column)}_le_{sanitize_token(format_float_token(q_threshold))}"
    if max_candidates > 0:
        suffix += f"_top{max_candidates}"
    suffix += site_wrapper.csubst_nonsyn_recode_output_suffix(nonsyn_recode)
    if pdb != "none":
        suffix += f"_pdb-{sanitize_token(pdb)}"
    return suffix


def archive_path_for_threshold(
    summary_prefix,
    out_dir,
    threshold,
    q_column,
    q_threshold,
    max_candidates,
    nonsyn_recode,
    pdb,
):
    suffix = output_suffix(
        q_column, q_threshold, max_candidates, nonsyn_recode, pdb
    )
    return Path(out_dir) / (
        f"{Path(summary_prefix).name}_candidate_sites_min_support_{threshold}_{suffix}.zip"
    )


def write_archive_manifest(rows, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    pd.DataFrame(rows, columns=ARCHIVE_MANIFEST_COLUMNS).to_csv(
        temporary, sep="\t", index=False
    )
    os.replace(temporary, path)


@contextmanager
def exclusive_run_lock(lock_path):
    lock_directory = Path(f"{lock_path}.d")
    lock_directory.parent.mkdir(parents=True, exist_ok=True)
    owner_path = lock_directory / "owner.json"
    owner_token = uuid.uuid4().hex
    owner = {
        "hostname": socket.gethostname(),
        "pid": os.getpid(),
        "token": owner_token,
        "created_unix": time.time(),
    }
    while True:
        try:
            lock_directory.mkdir(mode=0o700)
            owner_path.write_text(json.dumps(owner, sort_keys=True), encoding="utf-8")
            break
        except FileExistsError:
            if run_lock_is_stale(lock_directory):
                stale_directory = lock_directory.with_name(
                    f"{lock_directory.name}.stale-{uuid.uuid4().hex}"
                )
                try:
                    os.replace(lock_directory, stale_directory)
                except (FileNotFoundError, OSError):
                    continue
                shutil.rmtree(stale_directory, ignore_errors=True)
                continue
            time.sleep(0.2)

    heartbeat_stop = threading.Event()

    def update_heartbeat():
        while not heartbeat_stop.wait(RUN_LOCK_HEARTBEAT_SECONDS):
            try:
                os.utime(owner_path)
            except FileNotFoundError:
                return

    heartbeat_thread = threading.Thread(target=update_heartbeat, daemon=True)
    heartbeat_thread.start()
    try:
        yield
    finally:
        heartbeat_stop.set()
        heartbeat_thread.join(timeout=RUN_LOCK_HEARTBEAT_SECONDS + 1)
        try:
            current_owner = json.loads(owner_path.read_text(encoding="utf-8"))
        except (FileNotFoundError, json.JSONDecodeError, OSError):
            current_owner = {}
        if current_owner.get("token") == owner_token:
            shutil.rmtree(lock_directory, ignore_errors=True)


def run_lock_is_stale(lock_directory):
    lock_directory = Path(lock_directory)
    owner_path = lock_directory / "owner.json"
    try:
        owner = json.loads(owner_path.read_text(encoding="utf-8"))
        heartbeat_age = time.time() - owner_path.stat().st_mtime
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        try:
            heartbeat_age = time.time() - lock_directory.stat().st_mtime
        except FileNotFoundError:
            return False
        return heartbeat_age > RUN_LOCK_STALE_SECONDS

    if owner.get("hostname") == socket.gethostname():
        try:
            os.kill(int(owner["pid"]), 0)
        except (ProcessLookupError, ValueError, KeyError):
            return True
        except PermissionError:
            pass
    return heartbeat_age > RUN_LOCK_STALE_SECONDS


def run(args):
    args = validate_args(args)
    output_dir = Path(args.out_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    run_suffix = output_suffix(
        args.q_column,
        args.q_threshold,
        args.max_candidates,
        args.csubst_nonsyn_recode,
        args.pdb,
    )
    lock_path = output_dir / (
        f".{Path(args.summary_prefix).name}_candidate_sites_{run_suffix}.lock"
    )
    with exclusive_run_lock(lock_path):
        return run_locked(args, output_dir, run_suffix)


def run_locked(args, output_dir, run_suffix):
    summary_tables = discover_summary_tables(args.summary_prefix, args.min_support)
    manifest_path = output_dir / (
        f"{Path(args.summary_prefix).name}_candidate_sites_{run_suffix}_manifest.tsv"
    )
    work_root = output_dir / (
        f".{Path(args.summary_prefix).name}_candidate_sites_{run_suffix}.work"
    )
    cache_root = work_root / "cache"
    packages_root = work_root / "packages"
    materialization_parent = work_root / "materialized"
    trait_color_dir = work_root / "trait_colors"
    for path in (cache_root, packages_root, materialization_parent, trait_color_dir):
        path.mkdir(parents=True, exist_ok=True)

    threshold_candidates = {}
    all_traits = set()
    archive_rows = []
    for threshold, summary_path in summary_tables.items():
        candidates = load_threshold_candidates(
            summary_path=summary_path,
            minimum_support=threshold,
            q_column=args.q_column,
            q_threshold=args.q_threshold,
            max_candidates=args.max_candidates,
            csubst_nonsyn_recode=args.csubst_nonsyn_recode,
            pdb=args.pdb,
        )
        threshold_candidates[threshold] = candidates
        all_traits.update(candidates["trait"].dropna().astype(str).tolist())
        archive_path = archive_path_for_threshold(
            summary_prefix=args.summary_prefix,
            out_dir=output_dir,
            threshold=threshold,
            q_column=args.q_column,
            q_threshold=args.q_threshold,
            max_candidates=args.max_candidates,
            nonsyn_recode=args.csubst_nonsyn_recode,
            pdb=args.pdb,
        )
        archive_rows.append(
            {
                "min_support": threshold,
                "candidate_count": int(candidates.shape[0]),
                "summary_tsv": summary_path.name,
                "summary_sha256": file_sha256(summary_path),
                "analysis_engine_signature": analysis_engine_signature(),
                "csubst_version": csubst_version(),
                "csubst_signature": csubst_signature(),
                "runtime_dependency_versions": runtime_dependency_versions(),
                "q_column": args.q_column,
                "q_threshold": float(args.q_threshold),
                "archive_zip": archive_path.name,
                "status": "pending",
                "error": "",
            }
        )
    trait_color_paths = write_trait_color_tables(
        args.file_trait, all_traits, trait_color_dir
    ) if all_traits else {}
    write_archive_manifest(archive_rows, manifest_path)

    failures_detected = False
    for row_index, archive_row in enumerate(archive_rows):
        threshold = int(archive_row["min_support"])
        candidates = threshold_candidates[threshold]
        archive_path = output_dir / archive_row["archive_zip"]
        if archive_matches_source(archive_path, summary_tables[threshold]):
            archive_rows[row_index]["status"] = "existing"
            write_archive_manifest(archive_rows, manifest_path)
            print(f"min_support={threshold}: existing ZIP retained: {archive_path}", flush=True)
            continue
        if archive_path.exists():
            print(
                f"min_support={threshold}: existing ZIP is stale or incomplete and will be replaced: "
                f"{archive_path}",
                flush=True,
            )
        try:
            failures = ensure_candidate_analyses(
                candidates=candidates,
                cache_root=cache_root,
                dir_orthogroup=args.dir_orthogroup,
                materialization_parent=materialization_parent,
                trait_color_paths=trait_color_paths,
                nonsyn_recode=args.csubst_nonsyn_recode,
                pdb=args.pdb,
                ncpu=args.ncpu,
            )
        except Exception as error:
            failures_detected = True
            archive_rows[row_index]["status"] = "failed"
            archive_rows[row_index]["error"] = str(error)
            write_archive_manifest(archive_rows, manifest_path)
            break
        if failures:
            failures_detected = True
            archive_rows[row_index]["status"] = "failed"
            archive_rows[row_index]["error"] = "; ".join(
                f"{failure['candidate_id']}: {failure['error']}" for failure in failures
            )
            write_archive_manifest(archive_rows, manifest_path)
            break
        try:
            package_threshold(
                candidates=candidates,
                threshold=threshold,
                source_summary=summary_tables[threshold],
                archive_path=archive_path,
                packages_root=packages_root,
                cache_root=cache_root,
                q_column=args.q_column,
                q_threshold=args.q_threshold,
            )
        except Exception as error:
            failures_detected = True
            archive_rows[row_index]["status"] = "failed"
            archive_rows[row_index]["error"] = str(error)
            write_archive_manifest(archive_rows, manifest_path)
            break
        archive_rows[row_index]["status"] = "completed"
        write_archive_manifest(archive_rows, manifest_path)
        print(
            f"min_support={threshold}: {candidates.shape[0]:,} candidates packaged: {archive_path}",
            flush=True,
        )

    if failures_detected:
        raise RuntimeError(
            f"CSUBST scan candidate-site packaging failed. See: {manifest_path}"
        )
    if all(row["status"] in {"completed", "existing"} for row in archive_rows):
        shutil.rmtree(work_root, ignore_errors=True)
    return manifest_path


def main():
    args = parse_args()
    manifest_path = run(args)
    print(f"Candidate-site archive manifest written: {manifest_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
