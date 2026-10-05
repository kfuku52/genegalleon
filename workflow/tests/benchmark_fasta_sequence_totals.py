"""Compare exact FASTA-store byte totals in fresh processes in the GG runtime.

Before: python workflow/tests/benchmark_fasta_sequence_totals.py \
    --work-dir /tmp/gg-fasta-total --implementations legacy > before.json
After: repeat with --implementations legacy aggregate --expected-result before.json.
Input generation and cold/migration/warm ensure are measured separately from
total retrieval.  A cold index does not imply a flushed operating-system cache.
The optional Unicode control checks UTF-8 byte semantics; it is not biological
DNA.  This benchmark does not measure BLAST search or annotation accuracy.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import platform
import resource
import sqlite3
import statistics
import subprocess
import sys
import tempfile
import time
import zlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from workflow.support import fasta_sequence_store as store


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256")
    return digest.hexdigest()


def fingerprint(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def legacy_total(database: Path) -> int:
    """Match the selected TBLASTN core's existing decompressed-byte sum."""
    with sqlite3.connect(database.as_uri() + "?mode=ro", uri=True) as connection:
        return sum(len(zlib.decompress(row[0])) for row in connection.execute("SELECT sequence FROM sequences"))


def peak_rss() -> tuple[int, str]:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rss = int(value if sys.platform == "darwin" else value * 1024)
    status = Path("/proc/self/status")
    if status.is_file():
        for line in status.read_text().splitlines():
            if line.startswith("VmHWM:"):
                return int(line.split()[1]) * 1024, "proc_vm_hwm"
    return rss, "getrusage"


def index_state(directory: Path) -> dict:
    database, manifest = directory / "sequences.sqlite", directory / "manifest.json"
    if not database.exists():
        return {"exists": False}
    with sqlite3.connect(database.as_uri() + "?mode=ro", uri=True) as connection:
        metadata = dict(connection.execute("SELECT key, value FROM metadata"))
    return {"exists": True, "schema_version": metadata.get("schema_version"),
            "metadata_keys": sorted(metadata), "database_signature": store.file_signature(database),
            "manifest_sha256": sha256(manifest) if manifest.is_file() else None}


def worker(args: argparse.Namespace) -> dict:
    directory = args.work_dir.resolve()
    database, manifest = directory / "sequences.sqlite", directory / "manifest.json"
    before_state = index_state(directory) if args.worker_mode == "ensure" else None
    before_cpu, before_wall = time.process_time(), time.perf_counter()
    if args.worker_mode == "ensure":
        # Use the production atomic ensure path, including old-index migration.
        with contextlib.redirect_stdout(io.StringIO()):
            result = store.ensure(argparse.Namespace(database=database, manifest=manifest,
                source_list=directory / "sources.tsv", max_database_bytes=store.DEFAULT_MAX_DATABASE_BYTES,
                minimum_free_bytes=store.DEFAULT_MINIMUM_FREE_BYTES))
        if result:
            raise RuntimeError("FASTA store ensure failed")
        total = None
    elif args.implementation == "legacy":
        total = legacy_total(database)
    else:
        # Additive production API; there is deliberately no unvalidated fallback.
        if not hasattr(store, "total_sequence_bytes"):
            raise RuntimeError("aggregate requires fasta_sequence_store.total_sequence_bytes")
        total = store.total_sequence_bytes(database, manifest)
    wall, cpu = time.perf_counter() - before_wall, time.process_time() - before_cpu
    rss, rss_source = peak_rss()
    result = {"implementation": args.implementation, "operation": args.worker_mode,
              "wall_seconds": wall, "cpu_seconds": cpu, "peak_rss_bytes": rss,
              "peak_rss_source": rss_source, "total_sequence_bytes": total}
    if before_state is not None:
        after_state = index_state(directory)
        result.update(before_state=before_state, after_state=after_state,
                      state="cold_build" if not before_state["exists"] else
                      "unchanged" if before_state == after_state else "existing_index_update")
    result["process_cpu_seconds"] = time.process_time()
    return result


def run_worker(directory: Path, implementation: str, mode: str = "sum") -> dict:
    command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--worker-mode", mode,
               "--implementation", implementation, "--work-dir", str(directory)]
    started = time.perf_counter()
    output = subprocess.run(command, check=True, capture_output=True, text=True)
    result = json.loads(output.stdout)
    result["process_wall_seconds"] = time.perf_counter() - started
    return result


def prepare_inputs(directory: Path, args: argparse.Namespace) -> tuple[dict, float]:
    started = time.perf_counter()
    configuration = {"producer": "shake256-low2bit-v1", "species": args.species,
                     "records_per_species": args.records, "sequence_length": args.length,
                     "seed": args.seed, "unicode_control": args.unicode_control}
    metadata_path = directory / "benchmark_inputs.json"
    if metadata_path.exists():
        fixture = json.loads(metadata_path.read_text())
        if fixture["configuration"] != configuration:
            raise ValueError("Existing --work-dir was prepared with different workload parameters")
        for row in fixture["sources"]:
            if sha256(directory / row["filename"]) != row["sha256"]:
                raise RuntimeError("Saved benchmark FASTA source changed")
        return fixture, time.perf_counter() - started
    directory.mkdir(parents=True, exist_ok=True)
    if any(directory.iterdir()):
        raise ValueError("New --work-dir must be empty; refusing to replace existing data")
    sources, samples = [], []
    expected_sample = io.StringIO()
    alphabet = bytes(b"ACGT"[value % 4] for value in range(256))
    total_bytes = 0
    for species_index in range(args.species):
        species = f"Species_{species_index:03d}"
        path = directory / f"{species}.fa"
        with path.open("w", encoding="utf-8", newline="\n") as handle:
            for record in range(args.records):
                identifier = f"{species}_g{record:06d}"
                seed = f"{args.seed}:{species_index}:{record}".encode()
                sequence = hashlib.shake_256(seed).digest(args.length).translate(alphabet).decode("ascii")
                if args.unicode_control and species_index == 0 and record == 0:
                    sequence = "é" + sequence[1:]
                total_bytes += len(sequence.encode("utf-8"))
                store.write_record(handle, identifier, sequence)
                if record in {0, args.records - 1}:
                    samples.append(identifier)
                    store.write_record(expected_sample, identifier, sequence)
        sources.append({"species": species, "filename": path.name, "sha256": sha256(path),
                        "record_count": args.records, "file_bytes": path.stat().st_size})
    (directory / "sources.tsv").write_text("".join(
        f"{directory / row['filename']}\t{row['species']}\n" for row in sources), encoding="utf-8")
    (directory / "samples.txt").write_text("\n".join(samples) + "\n", encoding="utf-8")
    fixture = {"configuration": configuration, "sources": sources, "total_sequence_bytes": total_bytes,
               "sample_records": len(samples), "sample_fasta_sha256": hashlib.sha256(
                   expected_sample.getvalue().encode("utf-8")).hexdigest()}
    metadata_path.write_text(json.dumps(fixture, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return fixture, time.perf_counter() - started


def verify_outputs(directory: Path, fixture: dict) -> dict:
    for row in fixture["sources"]:
        if sha256(directory / row["filename"]) != row["sha256"]:
            raise RuntimeError("Benchmark FASTA source changed during measurement")
    manifest_path = directory / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    identity = [{key: row[key] for key in ("species", "sha256", "record_count")} for row in fixture["sources"]]
    if manifest["sources"] != identity:
        raise RuntimeError("FASTA-store manifest source identity differs from the producer")
    extracted = directory / "sample_extract.fa"
    with contextlib.redirect_stdout(io.StringIO()):
        status = store.extract(argparse.Namespace(database=directory / "sequences.sqlite",
            pattern_file=directory / "samples.txt", output=extracted, ignore_case=False,
            query_variants=False, prefix_species=False, require_all=True))
    sample_hash = sha256(extracted)
    if status or sample_hash != fixture["sample_fasta_sha256"]:
        raise RuntimeError("Sample extraction differs from the deterministic source FASTA")
    return {"input_sha256": fingerprint([fixture["configuration"], identity]),
            "source_identity": identity, "manifest_schema_version": manifest["schema_version"],
            "manifest_bytes_sha256": sha256(manifest_path), "sample_fasta_sha256": sample_hash,
            "total_sequence_bytes": fixture["total_sequence_bytes"]}


def benchmark(directory: Path, args: argparse.Namespace) -> dict:
    fixture, preparation_wall = prepare_inputs(directory, args)
    # Existing public manifest bytes must survive a private aggregate migration.
    manifest_path = directory / "manifest.json"
    saved_manifest_hash = sha256(manifest_path) if manifest_path.is_file() else None
    ensure_runs = [run_worker(directory, "legacy", "ensure"), run_worker(directory, "legacy", "ensure")]
    scientific_before = verify_outputs(directory, fixture)
    if saved_manifest_hash and scientific_before["manifest_bytes_sha256"] != saved_manifest_hash:
        raise RuntimeError("Public manifest bytes changed during existing-index ensure")
    database_hash_before = sha256(directory / "sequences.sqlite")
    measurements, warmups = [], []
    for implementation in args.implementations:
        for repetition in range(args.warmups + args.repeats):
            row = run_worker(directory, implementation)
            if row["total_sequence_bytes"] != fixture["total_sequence_bytes"]:
                raise RuntimeError(f"{implementation} byte total differs from the deterministic producer")
            row["repetition"] = repetition
            (warmups if repetition < args.warmups else measurements).append(row)
    scientific_after = verify_outputs(directory, fixture)
    if scientific_before != scientific_after or sha256(directory / "sequences.sqlite") != database_hash_before:
        raise RuntimeError("Read-only total measurements changed source, manifest, extraction or index bytes")
    if args.expected_result:
        expected = json.loads(args.expected_result.read_text())["scientific"]
        if expected != scientific_after:
            raise RuntimeError("Scientific identities differ from --expected-result")
    summaries = {}
    for implementation in args.implementations:
        rows = [row for row in measurements if row["implementation"] == implementation]
        summaries[implementation] = {"wall_seconds_median": statistics.median(row["wall_seconds"] for row in rows),
            "wall_seconds_min": min(row["wall_seconds"] for row in rows),
            "wall_seconds_max": max(row["wall_seconds"] for row in rows),
            "cpu_seconds_median": statistics.median(row["cpu_seconds"] for row in rows),
            "process_wall_seconds_median": statistics.median(row["process_wall_seconds"] for row in rows),
            "peak_rss_bytes_max": max(row["peak_rss_bytes"] for row in rows)}
    return {"schema": 1, "workload": fixture["configuration"], "scientific": scientific_after,
            "outputs_equivalent": True, "environment": {"python": sys.version, "platform": platform.platform(),
                "sqlite": sqlite3.sqlite_version, "zlib": zlib.ZLIB_RUNTIME_VERSION,
                "implementation_sha256": {path.name: sha256(path) for path in (Path(__file__), Path(store.__file__))}},
            "input_preparation_wall_seconds": preparation_wall, "store_preparation": ensure_runs,
            "database_bytes": (directory / "sequences.sqlite").stat().st_size,
            "database_sha256": database_hash_before, "repeats": args.repeats, "warmup_count": args.warmups,
            "summaries": summaries, "measurements": measurements, "warmups": warmups,
            "aggregate_speed_ratio": summaries["legacy"]["wall_seconds_median"] /
                max(1e-12, summaries["aggregate"]["wall_seconds_median"])
                if {"legacy", "aggregate"}.issubset(summaries) else None}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", type=int, default=9)
    parser.add_argument("--records", type=int, default=2000, help="FASTA records per species")
    parser.add_argument("--length", type=int, default=3000, help="Characters per sequence")
    parser.add_argument("--seed", type=int, default=918273)
    parser.add_argument("--unicode-control", action="store_true")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--implementations", nargs="+", choices=("legacy", "aggregate"), default=["legacy"])
    parser.add_argument("--work-dir", type=Path, help="Retain and reuse identical source FASTAs and index")
    parser.add_argument("--expected-result", type=Path, help="Require scientific identity with a saved before result")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--worker-mode", choices=("sum", "ensure"), default="sum", help=argparse.SUPPRESS)
    parser.add_argument("--implementation", choices=("legacy", "aggregate"), default="legacy", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if min(args.species, args.records, args.length, args.repeats) < 1 or args.warmups < 0:
        parser.error("species, records, length and repeats must be positive; warmups must be nonnegative")
    if args.worker:
        if not args.work_dir:
            parser.error("workers require --work-dir")
        result = worker(args)
    elif args.work_dir:
        result = benchmark(args.work_dir.resolve(), args)
    else:
        with tempfile.TemporaryDirectory(prefix="gg-fasta-total-benchmark-") as temporary:
            result = benchmark(Path(temporary), args)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
