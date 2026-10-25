#!/usr/bin/env python3
"""Compare complete database BH-FDR updates with identical logical outputs.

Run through workflow/tests/run_in_runtime.sh with a saved --baseline and --output.
Fixtures are built in the parent, then copied into each fresh worker before its
timer starts. Operation timing excludes imports, copying and output hashing;
Linux VmHWM is sampled before hashing. The default large profiles have two
million rows. --bh-source applies one saved BH implementation to both database
helpers to isolate their SQL-result buffering; omitting it uses the BH module
beside each helper for a combined comparison. --baseline-only establishes the
pre-change measurements.
"""

import argparse
import hashlib
import importlib
import importlib.util
import json
import platform
import resource
import shutil
import sqlite3
import statistics
import struct
import subprocess
import sys
import tempfile
import time
from pathlib import Path

SUPPORT = Path(__file__).resolve().parents[1] / "support"
P_COLUMN = "p_rate_enrichment_asymptotic"
PROFILES = ("dense", "missing", "textual", "small")
DEPENDENCIES = ("csubst_scan_support", "gene_family_output_store", "pgls_multiplicity", "shared_namespace_lock")


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def source_proofs(args):
    labels = ("baseline",) if args.baseline_only else ("baseline", "implementation")
    paths = {label: getattr(args, label) for label in labels}
    paths.update({label + "_bh": selected_bh_source(args, label) for label in labels})
    paths.update({name: SUPPORT / (name + ".py") for name in DEPENDENCIES})
    return {name: {"path": str(path.resolve()), "sha256": file_sha256(path)}
            for name, path in paths.items()}


def selected_bh_source(args, label):
    return args.bh_source or getattr(args, label).resolve().parent / "csubst_scan_support.py"


def verify_sources(proofs):
    for name, proof in proofs.items():
        if file_sha256(proof["path"]) != proof["sha256"]:
            raise ValueError("Benchmark source changed: " + name)


def write_fixture(path, profile, rows):
    """Use bounded insert batches; include exact rowids outside float64 precision."""
    probability_type = "" if profile == "textual" else " REAL"
    with sqlite3.connect(path) as connection:
        connection.execute(
            'CREATE TABLE aa_change (fixture_id INTEGER PRIMARY KEY, '
            + P_COLUMN + probability_type + ', "rowid" TEXT, orthogroup TEXT, trait TEXT, scan_match TEXT)'
        )
        for start in range(0, rows, 8192):
            values = []
            for index in range(start, min(start + 8192, rows)):
                rowid = -23 if index == 0 else (2 ** 63 - 2 if index == rows - 1 else 2 ** 53 + index)
                probability = (index % 1001) / 1000
                if profile == "missing" and index % 4:
                    probability = None
                elif profile == "textual":
                    probability = (None, probability, str(probability), index % 2)[index % 4]
                values.append((rowid, probability, "repeated-label", f"OG{index % 997:07d}",
                               f"trait{index % 5}", ("any2spe", "spe2spe")[index % 2]))
            connection.executemany('INSERT INTO aa_change VALUES (?, ?, ?, ?, ?, ?)', values)


def hash_value(digest, value):
    """Include SQLite value types and float bits without decimal round trips."""
    if value is None:
        digest.update(b"N")
    elif isinstance(value, int):
        digest.update(b"I" + struct.pack(">q", value))
    elif isinstance(value, float):
        digest.update(b"F" + struct.pack(">d", value))
    else:
        encoded = value.encode("utf-8") if isinstance(value, str) else value
        digest.update((b"S" if isinstance(value, str) else b"B") + struct.pack(">Q", len(encoded)) + encoded)


def fingerprint(path):
    digest = hashlib.sha256()
    tables = {}
    with sqlite3.connect(path) as connection:
        schema = connection.execute(
            "SELECT type, name, tbl_name, sql FROM sqlite_master ORDER BY type, name"
        ).fetchall()
        for record in schema:
            for value in record:
                hash_value(digest, value)
        schema_sha256 = digest.hexdigest()
        for table in ("aa_change", "aa_change_fdr_metadata"):
            table_digest = hashlib.sha256()
            columns = connection.execute(f'PRAGMA table_info("{table}")').fetchall()
            for record in columns:
                for value in record:
                    hash_value(table_digest, value)
            identity = "_rowid_" if table == "aa_change" else "rowid"
            cursor = connection.execute(f'SELECT {identity}, * FROM "{table}" ORDER BY {identity}')
            count = 0
            while records := cursor.fetchmany(1024):
                for record in records:
                    table_digest.update(b"R")
                    for value in record:
                        hash_value(table_digest, value)
                count += len(records)
            tables[table] = {"rows": count, "sha256": table_digest.hexdigest()}
            digest.update(table.encode() + table_digest.digest())
    return {"sha256": digest.hexdigest(), "schema_sha256": schema_sha256, "tables": tables}


def peak_rss_bytes():
    status = Path("/proc/self/status")
    if status.is_file():
        for line in status.read_text().splitlines():
            if line.startswith("VmHWM:"):
                return int(line.split()[1]) * 1024, "Linux /proc/self/status VmHWM"
        raise RuntimeError("Linux process status lacks VmHWM")
    scale = 1 if sys.platform == "darwin" else 1024
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale, "resource.ru_maxrss"


def worker(args):
    # Preload current dependencies: the saved helper can reside beside old
    # modules, while the selected database source is the only changed variable.
    sys.path.insert(0, str(SUPPORT))
    for name in DEPENDENCIES:
        importlib.import_module(name)
    before_sha256 = file_sha256(args.source)
    module = load(args.source, "benchmark_database_fdr")
    bh_source = args.bh_source or args.source.resolve().parent / "csubst_scan_support.py"
    bh_sha256 = file_sha256(bh_source)
    module.calculate_bh_fdr = load(bh_source, "benchmark_shared_bh").calculate_bh_fdr
    with tempfile.TemporaryDirectory(prefix="gg-database-fdr-worker-") as temporary:
        database = Path(temporary) / "output.sqlite3"
        shutil.copyfile(args.fixture, database)
        if file_sha256(database) != args.fixture_sha256:
            raise ValueError("Copied input database differs from frozen fixture")
        engine = module.sqlalchemy.create_engine("sqlite:///" + str(database))
        started = time.perf_counter()
        started_cpu = time.process_time()
        try:
            module.add_analytical_aa_change_fdr(engine)
            seconds = time.perf_counter() - started
            cpu_seconds = time.process_time() - started_cpu
            peak, peak_method = peak_rss_bytes()
        finally:
            engine.dispose()
        proof = fingerprint(database)
    if file_sha256(args.source) != before_sha256 or file_sha256(bh_source) != bh_sha256:
        raise ValueError("Benchmark source changed during worker")
    print(json.dumps({"seconds": seconds, "cpu_seconds": cpu_seconds, "peak_rss_bytes": peak,
                      "peak_rss_method": peak_method, "source_sha256": before_sha256,
                      "bh_source_sha256": bh_sha256, "input_sha256": args.fixture_sha256,
                      "logical_output": proof}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--implementation", type=Path, default=SUPPORT / "generate_orthogroup_database.py")
    parser.add_argument("--bh-source", type=Path, help="Use this BH source for both helpers; default: beside each helper.")
    parser.add_argument("--baseline-only", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--rows", type=int, default=2000000)
    parser.add_argument("--small-rows", type=int, default=128)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--profiles", nargs="+", choices=PROFILES, default=list(PROFILES))
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--source", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--fixture", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--fixture-sha256", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    if args.baseline is None or args.output is None:
        parser.error("--baseline and --output are required")
    if args.rows < 3 or args.small_rows < 3 or not 1 <= args.repeats <= 20:
        parser.error("row counts must be >= 3 and repeats must be 1..20")
    import numpy
    import pandas
    import sqlalchemy

    labels = ("baseline",) if args.baseline_only else ("baseline", "implementation")
    proofs = source_proofs(args)
    results = {"python": sys.version, "platform": platform.platform(), "numpy": numpy.__version__,
               "pandas": pandas.__version__, "sqlalchemy": sqlalchemy.__version__, "sqlite": sqlite3.sqlite_version,
               "warmup": 1, "repeats": args.repeats, "sources": proofs, "inputs": {}, "samples": {},
               "cache": "warm OS caches; fresh process and input database copy per sample",
               "timer_scope": "complete add_analytical_aa_change_fdr; excludes imports, setup and hashing",
               "peak_rss_scope": "fresh worker imports/setup and operation; sampled before output hashing",
               "bh_policy": ("one explicit BH source for both helpers" if args.bh_source
                             else "BH source beside each helper; combined comparison")}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="gg-database-fdr-fixtures-") as temporary:
        for profile in args.profiles:
            fixture = Path(temporary) / (profile + ".sqlite3")
            rows = args.small_rows if profile == "small" else args.rows
            write_fixture(fixture, profile, rows)
            input_sha256 = file_sha256(fixture)
            results["inputs"][profile] = {"rows": rows, "bytes": fixture.stat().st_size, "sha256": input_sha256}
            samples = {label: [] for label in labels}
            expected = None
            for trial in range(args.repeats + 1):
                for label in labels if trial % 2 == 0 else tuple(reversed(labels)):
                    verify_sources(proofs)
                    command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--source",
                               str(getattr(args, label).resolve()), "--fixture", str(fixture),
                               "--fixture-sha256", input_sha256, "--bh-source",
                               str(selected_bh_source(args, label).resolve())]
                    sample = json.loads(subprocess.check_output(command, text=True))
                    if expected is not None and sample["logical_output"] != expected:
                        raise ValueError("Complete output differs for " + profile + "/" + label)
                    expected = sample["logical_output"]
                    print(profile, label, "warmup" if trial == 0 else trial, json.dumps(sample), flush=True)
                    if trial:
                        samples[label].append(sample)
            if file_sha256(fixture) != input_sha256:
                raise ValueError("Frozen fixture changed")
            results["samples"][profile] = samples
            results["medians"] = {
                name: {label: {metric: statistics.median(sample[metric] for sample in group)
                               for metric in ("seconds", "cpu_seconds", "peak_rss_bytes")}
                       for label, group in groups.items()}
                for name, groups in results["samples"].items()
            }
            verify_sources(proofs)
            args.output.write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
