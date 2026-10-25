#!/usr/bin/env python3
"""Measure warm trait-index reads with saved/current helper sources.

Run through workflow/tests/run_in_runtime.sh. Index preparation runs in separate
processes; a measured worker fails if its loader is called. Wall/CPU time and
Linux VmHWM are captured before hashing the complete returned DataFrame. The
fixture and index caches remain in --fixtures for inspection and reuse.
"""

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import pickle
import platform
import resource
import sqlite3
import statistics
import subprocess
import sys
import time
from contextlib import closing
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROFILES = ("narrow", "wide", "sparse", "small")
DATASETS = {
    "narrow": {"rows": 150000, "numeric_columns": 12, "duplicate_every": 257, "encoding": "text"},
    "wide": {"rows": 100000, "numeric_columns": 128, "duplicate_every": 257, "encoding": "blob"},
    "small": {"rows": 64, "numeric_columns": 8, "duplicate_every": 17, "encoding": "text"},
}
REQUIRE_WARM = False
LOADER_CALLS = 0


def sha256_file(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def dataset_for(profile):
    return "wide" if profile == "sparse" else profile


def identifier(row, duplicate_every):
    ordinal = row - 1 if row and row % duplicate_every == 0 else row
    return f"g{ordinal:07d}"


def selected_ids(profile):
    shape = DATASETS[dataset_for(profile)]
    if profile == "sparse":
        # Include a duplicated identity, as well as distributed sparse hits.
        ordinals = [0, 256] + [i * shape["rows"] // 19 for i in range(1, 19)]
    else:
        ordinals = range(shape["rows"])
    return {identifier(row, shape["duplicate_every"]) for row in ordinals}


def load_fixture_table(path):
    """Use one stable loader identity in preparation and measured processes."""
    global LOADER_CALLS
    if REQUIRE_WARM:
        raise AssertionError("Measured trait-index read attempted a cold build")
    LOADER_CALLS += 1
    import numpy as np
    import pandas as pd

    frame = pd.read_csv(path, sep="\t", converters={0: str})
    frame.index = pd.Index(np.arange(len(frame), dtype=np.int64) * 3 + 7, name="source_row")
    return frame


def load_indexer(support):
    path = support.resolve() / "trait_table_index.py"
    # A baseline snapshot may contain just this helper; unchanged dependencies
    # then come from the active checkout, identically for both revisions.
    sys.path[:0] = [str(support.resolve()), str(ROOT / "workflow/support")]
    spec = importlib.util.spec_from_file_location("trait_index_benchmark", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def manifest_path(fixtures):
    return fixtures / "fixtures.json"


def prepare_fixtures(args):
    import numpy as np

    args.fixtures.mkdir(parents=True, exist_ok=True)
    path = manifest_path(args.fixtures)
    manifest = json.loads(path.read_text()) if path.exists() else {"schema": 1, "datasets": {}}
    if manifest.get("schema") != 1:
        raise ValueError("Unsupported fixture manifest")
    for name in dict.fromkeys(dataset_for(profile) for profile in args.profiles):
        shape = DATASETS[name]
        source = args.fixtures / (name + ".tsv")
        saved = manifest["datasets"].get(name)
        if saved is not None:
            if saved["shape"] != shape or saved["source"] != str(source.resolve()):
                raise ValueError("Fixture parameters changed: " + name)
            if sha256_file(source) != saved["sha256"]:
                raise ValueError("Fixture content changed: " + name)
            continue
        if source.exists():
            raise FileExistsError("Refusing to overwrite an unmanifested fixture: " + str(source))
        staging = source.with_suffix(".tsv.partial")
        if staging.exists():
            raise FileExistsError("Interrupted fixture preparation remains: " + str(staging))
        rng = np.random.default_rng(20261025 + shape["numeric_columns"])
        with staging.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene_id"] + [f"trait{i:03d}" for i in range(shape["numeric_columns"])])
            for start in range(0, shape["rows"], 4096):
                values = rng.uniform(-100, 100, size=(min(4096, shape["rows"] - start), shape["numeric_columns"]))
                for offset, row in enumerate(values):
                    ordinal = start + offset
                    if ordinal % 97 == 0:
                        row[2] = np.nan
                    writer.writerow([identifier(ordinal, shape["duplicate_every"])]
                                     + [format(float(value), ".17g") for value in row])
        os.replace(staging, source)
        manifest["datasets"][name] = {"shape": shape, "source": str(source.resolve()),
                                      "bytes": source.stat().st_size, "sha256": sha256_file(source)}
        staging_manifest = path.with_suffix(".json.partial")
        staging_manifest.write_text(json.dumps(manifest, indent=2) + "\n")
        os.replace(staging_manifest, path)
    print(json.dumps(manifest))


def fixture_details(args):
    manifest = json.loads(manifest_path(args.fixtures).read_text())
    name = dataset_for(args.profile)
    return name, manifest["datasets"][name]


def prepare_index(args):
    module = load_indexer(args.support)
    name, details = fixture_details(args)
    source = Path(details["source"])
    result = module.read_subset(source, set(), load_fixture_table, args.cache)
    assert result.empty
    databases = list(args.cache.glob("*.sqlite3"))
    if len(databases) != 1:
        raise AssertionError("Expected exactly one index in the source-specific namespace")
    database = databases[0]
    with closing(sqlite3.connect(database.as_uri() + "?mode=ro", uri=True)) as connection:
        metadata = connection.execute("SELECT payload FROM metadata").fetchone()[0]
        rows, text_rows, blob_rows = connection.execute(
            "SELECT COUNT(*), SUM(typeof(payload) = 'text'), SUM(typeof(payload) = 'blob') FROM rows"
        ).fetchone()
    assert rows == DATASETS[name]["rows"]
    assert (text_rows if DATASETS[name]["encoding"] == "text" else blob_rows) == rows
    print(json.dumps({"database": str(database), "database_bytes": database.stat().st_size,
                      "database_sha256": sha256_file(database),
                      "metadata_sha256": hashlib.sha256(metadata.encode()).hexdigest(),
                      "rows": rows, "text_rows": text_rows, "blob_rows": blob_rows,
                      "loader_calls": LOADER_CALLS, "fixture_sha256": details["sha256"]}))


def peak_rss():
    status = Path("/proc/self/status")
    if status.is_file():
        for line in status.read_text().splitlines():
            if line.startswith("VmHWM:"):
                fields = line.split()
                if fields[2] != "kB":
                    raise ValueError("Unexpected VmHWM unit")
                return int(fields[1]) * 1024, "Linux /proc/self/status VmHWM"
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value * (1 if sys.platform == "darwin" else 1024), "resource.ru_maxrss"


def frame_digest(frame):
    """Hash every value in row order, including duplicate rows and exact floats."""
    metadata = {"columns": list(frame.columns), "dtypes": list(map(str, frame.dtypes)),
                "index_type": type(frame.index).__name__, "index_dtype": str(frame.index.dtype),
                "index_name": frame.index.name, "shape": list(frame.shape)}
    proof = hashlib.sha256(json.dumps(metadata, sort_keys=True).encode())
    for row in frame.itertuples(index=True, name=None):
        payload = pickle.dumps(row, protocol=4)
        proof.update(len(payload).to_bytes(8, "big"))
        proof.update(payload)
    return proof.hexdigest(), metadata


def worker(args):
    global REQUIRE_WARM
    module = load_indexer(args.support)
    _, details = fixture_details(args)
    identifiers = selected_ids(args.profile)
    REQUIRE_WARM = True
    started = time.perf_counter()
    started_cpu = time.process_time()
    frame = module.read_subset(Path(details["source"]), identifiers, load_fixture_table, args.cache)
    cpu_seconds = time.process_time() - started_cpu
    seconds = time.perf_counter() - started
    rss_bytes, rss_source = peak_rss()
    # Hashing/verification may allocate memory: capture the measurement first.
    digest, metadata = frame_digest(frame)
    if args.profile != "sparse":
        assert len(frame) == DATASETS[dataset_for(args.profile)]["rows"]
    assert LOADER_CALLS == 0
    print(json.dumps({"seconds": seconds, "cpu_seconds": cpu_seconds,
                      "peak_rss_bytes": rss_bytes, "rss_source": rss_source,
                      "digest": digest, "frame": metadata, "selected_ids": len(identifiers),
                      "duplicate_result_rows": len(frame) - frame.iloc[:, 0].nunique(),
                      "loader_calls": LOADER_CALLS, "fixture_sha256": details["sha256"]}))


def child(mode, fixtures, *extra):
    command = [sys.executable, str(Path(__file__).resolve()), "--mode", mode,
               "--fixtures", str(fixtures), *map(str, extra)]
    return json.loads(subprocess.check_output(command, text=True))


def compare(args):
    import numpy as np
    import pandas as pd

    if args.baseline_support is None or args.output is None:
        raise ValueError("--baseline-support and --output are required")
    if not 1 <= args.repeats <= 20:
        raise ValueError("--repeats must be 1..20")
    args.fixtures = args.fixtures.resolve()
    manifest = child("prepare-fixtures", args.fixtures, "--profiles", *args.profiles)
    supports = {"baseline": args.baseline_support.resolve()}
    if not args.baseline_only:
        supports["implementation"] = args.implementation_support.resolve()
    sources = {label: {"support": str(support), "sha256": sha256_file(support / "trait_table_index.py")}
               for label, support in supports.items()}
    indexes, caches = {}, {}
    for label, support in supports.items():
        indexes[label], caches[label] = {}, {}
        for name in dict.fromkeys(dataset_for(profile) for profile in args.profiles):
            cache = args.fixtures / "indexes" / (label + "-" + sources[label]["sha256"]) / name
            caches[label][name] = cache
            indexes[label][name] = child("prepare-index", args.fixtures, "--profile", name,
                                        "--support", support, "--cache", cache)
            assert indexes[label][name]["fixture_sha256"] == manifest["datasets"][name]["sha256"]
    for name in dict.fromkeys(dataset_for(profile) for profile in args.profiles):
        assert len({group[name]["metadata_sha256"] for group in indexes.values()}) == 1, name
        assert len({group[name]["database_sha256"] for group in indexes.values()}) == 1, name
    results = {"schema": 1, "python": sys.version, "platform": platform.platform(),
               "pandas": pd.__version__, "numpy": np.__version__,
               "benchmark_sha256": sha256_file(__file__), "sources": sources,
               "fixture_manifest": manifest, "indexes": indexes,
               "warmups": 1, "repeats": args.repeats, "baseline_only": args.baseline_only,
               "measurement": "Fresh process; warm content-verified index; time/RSS before result hashing; imports included in RSS",
               "samples": {}, "warmup_samples": {}}
    for profile in args.profiles:
        samples = {label: [] for label in supports}
        warmups = {}
        for trial in range(args.repeats + 1):
            labels = list(supports) if trial % 2 == 0 else list(reversed(supports))
            for label in labels:
                sample = child("worker", args.fixtures, "--profile", profile,
                               "--support", supports[label], "--cache", caches[label][dataset_for(profile)])
                print(json.dumps({"profile": profile, "source": label,
                                  "trial": "warmup" if trial == 0 else trial, **sample}), flush=True)
                if trial:
                    samples[label].append(sample)
                else:
                    warmups[label] = sample
        signatures = {sample["digest"] for group in samples.values() for sample in group}
        signatures.update(sample["digest"] for sample in warmups.values())
        assert len(signatures) == 1, "Returned DataFrames differ: " + profile
        results["samples"][profile] = samples
        results["warmup_samples"][profile] = warmups
    for name in dict.fromkeys(dataset_for(profile) for profile in args.profiles):
        assert sha256_file(manifest["datasets"][name]["source"]) == manifest["datasets"][name]["sha256"], name
    results["medians"] = {
        profile: {label: {metric: statistics.median(sample[metric] for sample in group)
                          for metric in ("seconds", "cpu_seconds", "peak_rss_bytes")}
                  for label, group in groups.items()}
        for profile, groups in results["samples"].items()
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results["medians"], indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-support", type=Path)
    parser.add_argument("--implementation-support", type=Path, default=ROOT / "workflow/support")
    parser.add_argument("--fixtures", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--profiles", nargs="+", choices=PROFILES, default=list(PROFILES))
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--baseline-only", action="store_true")
    parser.add_argument("--mode", choices=("compare", "prepare-fixtures", "prepare-index", "worker"),
                        default="compare", help=argparse.SUPPRESS)
    parser.add_argument("--profile", choices=PROFILES, help=argparse.SUPPRESS)
    parser.add_argument("--support", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--cache", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.mode in {"prepare-index", "worker"} and any(value is None for value in (args.profile, args.support, args.cache)):
        parser.error("Internal workers require --profile, --support, and --cache")
    {"compare": compare, "prepare-fixtures": prepare_fixtures,
     "prepare-index": prepare_index, "worker": worker}[args.mode](args)


if __name__ == "__main__":
    main()
