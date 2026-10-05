"""Benchmark verification of a real published representative-input bundle.

Run in the GeneGalleon runtime. Fixture construction and plan/run publication
are excluded from verification measurements. The held-out three-species,
seven-locus fixture gains an unannotated N-only contig to enlarge each genome
without changing original sequence bytes, annotated coordinates or CDS. Reuse --work-dir for an exact
before/after bundle; --expected-report checks source/scientific equivalence.

Each sample uses a fresh Python process running the actual verify-inputs CLI
through runpy. CLI timing includes its imports and content verification; outer
subprocess timing additionally includes interpreter/instrumentation startup.
Filesystem caches are warm after preparation/fingerprinting; they are not
flushed. Linux RSS comes from this fresh process' own VmHWM, avoiding inherited
parent high-water marks. No integrity check is bypassed or monkeypatched.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import importlib
import importlib.metadata
import io
import json
import os
import platform
import resource
import runpy
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

TESTS = Path(__file__).resolve().parent
SUPPORT = TESTS.parent / "support"
SCRIPT = SUPPORT / "gene_model_refinement.py"
FIXTURE = TESTS / "test_gene_model_refinement_runtime.py"
MIB = 1024 ** 2


def sha256(path, limit=None):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        remaining = limit
        while remaining is None or remaining:
            block = handle.read(MIB if remaining is None else min(MIB, remaining))
            if not block:
                if remaining:
                    raise ValueError("Source genome prefix was truncated")
                break
            digest.update(block)
            if remaining is not None:
                remaining -= len(block)
    return digest.hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def table(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def cli(*command):
    result = subprocess.run([sys.executable, str(SCRIPT), *map(str, command)], text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError("Fixture publication failed: " + " ".join(map(str, command))
                           + "\n" + result.stdout[-8000:] + result.stderr[-8000:])


def pad_genome(genome, genome_bytes):
    """Append one valid 80-column contig, keeping the original bytes intact."""
    initial_bytes, initial_hash = genome.stat().st_size, sha256(genome)
    with genome.open("rb") as handle:
        handle.seek(-1, os.SEEK_END)
        if handle.read(1) != b"\n":
            raise ValueError("Fixture genome must end with a newline")
    header = b">benchmark_padding\n"
    remaining = genome_bytes - initial_bytes - len(header)
    if remaining < 2:
        raise ValueError("--genome-mib must allow at least one padding base")
    # A remainder of one byte cannot encode a nonempty final FASTA line.
    # Add a harmless two-byte description so the final sequence line is valid.
    if remaining % 81 == 1:
        header = b">benchmark_padding N\n"
        remaining -= 2
    full_lines, last_bytes = divmod(remaining, 81)
    bases = full_lines * 80 + (last_bytes - 1 if last_bytes else 0)
    with genome.open("ab") as handle:
        handle.write(header)
        chunk_lines = min(MIB // 81, full_lines)
        block = (b"N" * 80 + b"\n") * chunk_lines
        while full_lines:
            count = min(chunk_lines, full_lines)
            handle.write(block[:count * 81])
            full_lines -= count
        if last_bytes:
            handle.write(b"N" * (last_bytes - 1) + b"\n")
    if genome.stat().st_size != genome_bytes or sha256(genome, initial_bytes) != initial_hash:
        raise RuntimeError("Genome padding changed original sequence bytes")
    return {"original_genome_bytes": initial_bytes, "original_genome_sha256": initial_hash,
            "source_prefix_unchanged": True, "padding_contig": "benchmark_padding",
            "padding_bases": bases, "padding_fasta_bytes": genome_bytes - initial_bytes,
            "padding_line_width": 80, "padding_header": header.decode().rstrip("\n")}


def prepare(directory, genome_bytes):
    """Construct once, or reuse recorded immutable benchmark sources/bundle."""
    metadata_path = directory / "benchmark_workload.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text())
        if metadata.get("schema") != 1 or metadata["genome_bytes_per_species"] != genome_bytes:
            raise ValueError("Existing benchmark workload differs; use a new --work-dir")
        if not (directory / "run/effective/receipt.json").is_file():
            raise ValueError("Existing benchmark publication is incomplete; use a new --work-dir")
        if metadata["source_hashes"] != {relative: sha256(directory / relative) for relative in metadata["source_hashes"]}:
            raise ValueError("Existing benchmark source bytes changed")
        return metadata, True
    if directory.exists() and any(directory.iterdir()):
        raise ValueError("Benchmark --work-dir must be empty or contain its recorded workload")
    directory.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(TESTS))
    fixture = importlib.import_module("test_gene_model_refinement_runtime")
    source_dir = directory / "sources"
    source_dir.mkdir()
    inputs, edges, _rna, sources, truth = fixture.truth_fixture(source_dir, rna=False)
    prefix_evidence = []
    for row in sources:
        genome = Path(row["genome"])
        prefix_evidence.append({"species": row["species"], **pad_genome(genome, genome_bytes)})
    source_paths = [inputs, edges, *(Path(row[role]) for row in sources for role in ("cds", "gff", "genome"))]
    source_hashes = {str(path.relative_to(directory)): sha256(path) for path in source_paths}
    semantic_sources = [{"species": row["species"], "genetic_code": row["genetic_code"],
                         **{role + "_sha256": sha256(row[role]) for role in ("cds", "gff", "genome")}}
                        for row in sources]
    original_cds_gff = {relative: digest for relative, digest in source_hashes.items()
                        if relative.endswith((".cds.fa", ".gff3"))}
    root = directory / "run"
    cli("plan", "--inputs", inputs, "--edges", edges, "--output", root, "--mode", "off", "--policy", "longest")
    cli("run", "--output", root, "--cpus", 1)
    if source_hashes != {relative: sha256(directory / relative) for relative in source_hashes}:
        raise RuntimeError("Published workload modified source bytes")
    metadata = {"schema": 1, "species": len(sources), "loci_per_species": len(truth) // len(sources),
                "genome_bytes_per_species": genome_bytes, "total_genome_source_bytes": genome_bytes * len(sources),
                "parameters": {"mode": "off", "policy": "longest", "rna": False, "cpus": 1},
                "genome_padding_method": "new_unannotated_contig_uniform_80nt_lines",
                "fixture_sha256": sha256(FIXTURE), "original_genome_prefixes": prefix_evidence,
                "source_hashes": source_hashes, "source_cds_gff_sha256": original_cds_gff,
                "semantic_sources": semantic_sources, "edges_sha256": sha256(edges),
                "input_sha256": fingerprint({"sources": semantic_sources, "edges_sha256": sha256(edges)})}
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return metadata, False


def scientific_outputs(effective):
    """Include every bundle sequence/GFF/map, with path-independent fields."""
    manifest = table(effective / "inputs.tsv")
    roles = ("cds", "protein", "gff", "genome", "analysis_cds", "analysis_gff", "representative_map")
    normalized_manifest = [{key: row[key] for key in ("species", "genetic_code", *(role + "_sha256" for role in roles))}
                           for row in manifest]
    files = {str(path.relative_to(effective)): sha256(path) for path in sorted(effective.rglob("*"))
             if path.is_file() and (path.suffix in {".fa", ".fasta", ".gff", ".gff3"}
                                   or path.name in {"representative_map.tsv", "species_genetic_code.tsv", "changes.json"})}
    return {"sha256": fingerprint({"manifest": normalized_manifest, "files": files}),
            "manifest_sha256": sha256(effective / "inputs.tsv"),
            "normalized_manifest_sha256": fingerprint(normalized_manifest),
            "manifest": normalized_manifest, "files": files}


def peak_rss():
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    value = int(usage if sys.platform == "darwin" else usage * 1024)
    status = Path("/proc/self/status")
    if status.is_file():
        for line in status.read_text().splitlines():
            if line.startswith("VmHWM:"):
                return int(line.split()[1]) * 1024, "proc_vm_hwm"
    return value, "getrusage"


def worker(manifest, metrics):
    sys.path.insert(0, str(SUPPORT))
    command = ["verify-inputs", "--inputs", str(manifest), "--field", "layout"]
    sys.argv = [str(SCRIPT), *command]
    captured = io.StringIO()
    before_cpu, before_wall = time.process_time(), time.perf_counter()
    with contextlib.redirect_stdout(captured):
        runpy.run_path(str(SCRIPT), run_name="__main__")
    wall, cpu = time.perf_counter() - before_wall, time.process_time() - before_cpu
    maximum, source = peak_rss()
    Path(metrics).write_text(json.dumps({"command": [sys.executable, str(SCRIPT), *command],
                                        "layout": captured.getvalue().rstrip("\n"),
                                        "wall_seconds": wall, "cpu_seconds": cpu,
                                        "peak_rss_bytes": maximum, "peak_rss_source": source}, sort_keys=True) + "\n")


def measure(manifest, scratch, repetition):
    metrics = scratch / f"{repetition}.json"
    command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--manifest", str(manifest), "--metrics", str(metrics)]
    before = time.perf_counter()
    process = subprocess.run(command, text=True, capture_output=True)
    wall = time.perf_counter() - before
    if process.returncode:
        raise RuntimeError("Verification CLI failed:\n" + process.stdout[-8000:] + process.stderr[-8000:])
    result = json.loads(metrics.read_text())
    return {**result, "subprocess_wall_seconds": wall, "worker_command": command,
            "stderr": process.stderr, "repetition": repetition}


def summarize(rows):
    result = {key: {"median": statistics.median(row[key] for row in rows),
                    "min": min(row[key] for row in rows), "max": max(row[key] for row in rows)}
              for key in ("wall_seconds", "cpu_seconds", "subprocess_wall_seconds", "peak_rss_bytes")}
    return {"repeats": len(rows), **result}


def benchmark(args, directory):
    before_preparation = time.perf_counter()
    workload, reused = prepare(directory, round(args.genome_mib * MIB))
    preparation_seconds = time.perf_counter() - before_preparation
    effective, manifest = directory / "run/effective", directory / "run/effective/inputs.tsv"
    published_plan = json.loads((directory / "run/plan.json").read_text())
    tracked = {str(SUPPORT / name) for name in published_plan["request"]["implementation"]}
    tracked.update(map(str, (Path(__file__).resolve(), FIXTURE, TESTS / "test_gene_model_refinement.py")))
    implementation = {path: sha256(path) for path in sorted(tracked)}
    outputs = scientific_outputs(effective)
    measurements, warmups = [], []
    with tempfile.TemporaryDirectory(prefix=".verification-samples-", dir=directory) as temporary:
        for repetition in range(args.warmups + args.repeats):
            row = measure(manifest, Path(temporary), repetition)
            if len(row["layout"].split("\t")) != 6:
                raise RuntimeError("Unexpected representative layout schema")
            (warmups if repetition < args.warmups else measurements).append(row)
    layouts = {row["layout"] for row in [*warmups, *measurements]}
    if len(layouts) != 1:
        raise RuntimeError("Verification returned different layouts")
    layout = next(iter(layouts))
    relative_layout = [str(Path(path).relative_to(effective)) for path in layout.split("\t")]
    if outputs != scientific_outputs(effective):
        raise RuntimeError("Verification changed published scientific bytes or manifest")
    if workload["source_hashes"] != {relative: sha256(directory / relative) for relative in workload["source_hashes"]}:
        raise RuntimeError("Verification changed source bytes")
    if implementation != {path: sha256(path) for path in implementation}:
        raise RuntimeError("Implementation changed during benchmark")
    summary = summarize(measurements)
    report = {"schema": 1, "benchmark": "representative_input_verification", "work_dir": str(directory),
              "methodology": __doc__, "workload": workload, "repeats": args.repeats, "warmups": args.warmups,
              "preparation": {"wall_seconds": preparation_seconds, "excluded_from_measurements": True, "bundle_reused": reused},
              "environment": {"python": sys.version, "platform": platform.platform(),
                              "versions": {name: importlib.metadata.version(name) for name in ("biopython", "pysam", "pytest")},
                              "runtime": {"image_label": args.runtime_image_label or os.getenv("GG_CONTAINER_DOCKER_IMAGE"),
                                          "image_id": args.runtime_image_id, "kind": os.getenv("GG_TEST_RUNTIME")},
                              "implementation_sha256": implementation},
              "published_implementation_sha256": published_plan["request"]["implementation"],
              "main_helper_sha256": sha256(SCRIPT), "benchmark_helper_sha256": sha256(__file__),
              "current_fixture_sha256": sha256(FIXTURE), "layout": layout, "layout_relative_to_effective": relative_layout,
              "source_bytes_unchanged": True, "published_scientific_bytes_unchanged": True,
              "scientific_outputs": outputs, "summary": summary, "measurements": measurements,
              "warmup_measurements": warmups, "comparison": None}
    if args.expected_report:
        previous = json.loads(args.expected_report.read_text())
        if previous["workload"]["input_sha256"] != workload["input_sha256"]:
            raise RuntimeError("Before/after benchmark scientific source inputs differ")
        if previous["scientific_outputs"]["sha256"] != outputs["sha256"]:
            raise RuntimeError("Before/after benchmark scientific bundle outputs differ")
        if previous["layout_relative_to_effective"] != relative_layout:
            raise RuntimeError("Before/after verification layout schemas differ")
        same_directory = previous["work_dir"] == str(directory)
        if same_directory and previous["layout"] != layout:
            raise RuntimeError("Before/after verification changed the exact returned layout")
        comparable = all(previous["environment"][key] == report["environment"][key]
                         for key in ("python", "platform", "versions", "runtime"))
        comparable = comparable and previous["benchmark_helper_sha256"] == report["benchmark_helper_sha256"]
        report["comparison"] = {"previous_report": str(args.expected_report), "scientific_inputs_equivalent": True,
                                "scientific_outputs_equivalent": True, "relative_layout_equivalent": True,
                                "raw_layouts_equal": previous["layout"] == layout, "environment_method_comparable": comparable,
                                "cli_wall_speed_ratio": previous["summary"]["wall_seconds"]["median"] /
                                summary["wall_seconds"]["median"] if comparable else None}
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-mib", type=float, default=64, help="Approximate target FASTA bytes per species in MiB (default:64)")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--work-dir", type=Path, help="Retain/reuse exact prepared sources and published bundle for before/after runs")
    parser.add_argument("--output", type=Path, help="New JSON report path; existing reports are never overwritten")
    parser.add_argument("--expected-report", type=Path, help="Before report whose scientific inputs, bundle and layout must match")
    parser.add_argument("--runtime-image-label")
    parser.add_argument("--runtime-image-id", help="Image identity recorded by the host runtime freshness wrapper")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--manifest", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--metrics", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        if not args.manifest or not args.metrics:
            parser.error("Worker requires manifest and metrics")
        worker(args.manifest, args.metrics)
        return
    if not args.output or args.output.exists():
        parser.error("Provide a new --output report path")
    if not args.genome_mib > 0 or args.repeats < 1 or args.warmups < 0:
        parser.error("genome-mib/repeats must be positive and warmups nonnegative")
    if args.work_dir:
        report = benchmark(args, args.work_dir.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix="gg-representative-inputs-benchmark-") as temporary:
            report = benchmark(args, Path(temporary))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        json.dump(report, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"report": str(args.output.resolve()), "summary": report["summary"], "comparison": report["comparison"]}, sort_keys=True))


if __name__ == "__main__":
    main()
