"""Bounded, deterministic isoform-selection benchmark; run in the GG runtime.

Example: python workflow/tests/benchmark_gene_model_selection.py --families 32
Each implementation runs in a fresh child process.  The uncached implementation
is the behavioral oracle; equivalence is checked before reporting a speed ratio.
This measures generated sparse graphs, not real-genome annotation accuracy or
protein-to-genome prediction, and does not extrapolate their runtime.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import random
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from workflow.support.gene_model_selection import select_representatives
from workflow.support.gene_model_store import build_store, select_from_store


def workload(families: int, species_count: int, protein_length: int) -> tuple[list[dict], list[dict]]:
    if families < 1 or species_count < 3 or protein_length < 30:
        raise ValueError("Use at least one family, three species and 30 amino acids")
    species = [f"Species_{index:03d}" for index in range(species_count)]
    catalogs = [{"schema": 1, "species": name, "loci": []} for name in species]
    edges = []
    alphabet = "ACDEFGHIKLNPQRSTVWY"
    for family in range(families):
        randomizer = random.Random(8271 + family)
        core = "M" + "".join(randomizer.choices(alphabet, k=protein_length - 1))
        gene = f"g{family:06d}"
        for index, name in enumerate(species):
            tail = "".join(random.Random(family * 10000 + index).choices(alphabet, k=protein_length // 4))
            insertion = "".join(random.Random(family * 10000 + index + 100).choices(alphabet, k=15))
            middle = protein_length // 2
            proteins = [core + tail, core, core[:middle], core[:middle] + insertion + core[middle:]]
            candidates = [{"candidate_id": f"{name}_{gene}_t{transcript}",
                           "source_transcript_id": f"{name}_{gene}_t{transcript}",
                           "protein": protein, "cds": "ATG" * len(protein),
                           "blocks": [[family * 100000, family * 100000 + len(protein) * 3, 0]],
                           "quality": {"usable": True, "valid_orf": True}, "origin": "original"}
                          for transcript, protein in enumerate(proteins)]
            # A subset has independent complete-transcript evidence; the other
            # species retain genuinely ambiguous alternatives in this workload.
            if index % 3 == 0:
                candidates[1]["quality"]["full_length_supported"] = True
            catalogs[index]["loci"].append({"species": name, "gene_id": gene, "candidates": candidates})
        pairs = {tuple(sorted((index, (index + distance) % species_count))) for index in range(species_count)
                 for distance in [1, 2]}
        for a, b in sorted(pairs):
            edges.append({"species_a": species[a], "gene_a": gene,
                          "species_b": species[b], "gene_b": gene, "weight": 1.0,
                          "evidence": "generated_frozen_copy_correspondence"})
    return catalogs, edges


def fingerprint(value: dict | list) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def worker(args: argparse.Namespace) -> dict:
    if args.implementation == "store":
        if not args.store_input:
            raise ValueError("Store workers require --store-input")
        directory = Path(args.store_input)
        inputs = json.loads((directory / "benchmark_inputs.json").read_text())
        edges, input_hash = inputs["edges"], inputs["input_sha256"]
    else:
        catalogs, edges = workload(args.families, args.species, args.length)
        input_hash = fingerprint([catalogs, edges])
    before_cpu, before_wall = time.process_time(), time.perf_counter()
    if args.implementation == "store":
        result = select_from_store(directory / "loci.sqlite", edges)
    else:
        result = select_representatives(catalogs, edges, policy="longest" if args.implementation == "longest" else "conserved",
                                       cache_pair_scores=args.implementation != "reference")
    wall_seconds, cpu_seconds = time.perf_counter() - before_wall, time.process_time() - before_cpu
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # Linux task accounting can retain the pre-exec parent's peak after fork.
    # VmHWM belongs to this process' current memory map and excludes the driver's
    # separate input preparation.  getrusage is the fallback on other hosts.
    peak_rss_source = "getrusage"
    peak_rss_bytes = int(peak_rss if sys.platform == "darwin" else peak_rss * 1024)
    status = Path("/proc/self/status")
    if status.is_file():
        for line in status.read_text().splitlines():
            if line.startswith("VmHWM:"):
                peak_rss_bytes = int(line.split()[1]) * 1024
                peak_rss_source = "proc_vm_hwm"
                break
    scientific_result = {key: value for key, value in result.items() if key != "metrics"}
    return {"implementation": args.implementation, "input_sha256": input_hash,
            "output_sha256": fingerprint(scientific_result), "wall_seconds": wall_seconds,
            "cpu_seconds": cpu_seconds, "peak_rss_bytes": peak_rss_bytes, "peak_rss_source": peak_rss_source,
            "metrics": result["metrics"]}


def prepare_store(directory: Path, args: argparse.Namespace) -> dict:
    """Prepare disk inputs once; fresh workers measure selector RSS separately."""
    before_wall = time.perf_counter()
    catalogs, edges = workload(args.families, args.species, args.length)
    input_hash = fingerprint([catalogs, edges])
    catalog_dirs = []
    for catalog in catalogs:
        catalog_dir = directory / catalog["species"]
        catalog_dir.mkdir()
        metadata = {key: value for key, value in catalog.items() if key != "loci"}
        metadata["summary"] = {"loci": len(catalog["loci"]),
                               "candidates": sum(len(locus["candidates"]) for locus in catalog["loci"])}
        (catalog_dir / "catalog_metadata.json").write_text(json.dumps(metadata))
        with (catalog_dir / "loci.jsonl").open("w") as handle:
            for locus in catalog["loci"]:
                handle.write(json.dumps(locus, sort_keys=True, separators=(",", ":")) + "\n")
        catalog_dirs.append(catalog_dir)
    del catalogs
    (directory / "benchmark_inputs.json").write_text(json.dumps({"input_sha256": input_hash, "edges": edges}))
    receipt = build_store(catalog_dirs, directory / "loci.sqlite")
    return {"wall_seconds": time.perf_counter() - before_wall, "database_bytes": (directory / "loci.sqlite").stat().st_size,
            "loci": receipt["loci_count"], "candidates": receipt["candidate_count"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--families", type=int, default=32)
    parser.add_argument("--species", type=int, default=9)
    parser.add_argument("--length", type=int, default=300)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    implementations = ["reference", "cached", "store", "longest"]
    parser.add_argument("--implementation", choices=implementations, default="cached")
    parser.add_argument("--implementations", choices=implementations, nargs="+", default=["reference", "cached", "longest"])
    parser.add_argument("--store-input", help=argparse.SUPPRESS)
    parser.add_argument("--expected-output-sha256", help="Check conserved output against a previous same-input result")
    args = parser.parse_args()
    if args.repeats < 1 or args.warmups < 0:
        parser.error("repeats must be positive and warmups nonnegative")
    if args.worker:
        print(json.dumps(worker(args), sort_keys=True, allow_nan=False))
        return
    measurements, preparation = [], None
    with tempfile.TemporaryDirectory(prefix="gg-selection-benchmark-") as temporary:
        if "store" in args.implementations:
            preparation = prepare_store(Path(temporary), args)
        for implementation in args.implementations:
            for repetition in range(args.warmups + args.repeats):
                command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--implementation", implementation,
                           "--families", str(args.families), "--species", str(args.species), "--length", str(args.length)]
                if implementation == "store":
                    command += ["--store-input", temporary]
                output = subprocess.run(command, check=True, capture_output=True, text=True)
                measurement = json.loads(output.stdout)
                if repetition >= args.warmups:
                    measurements.append(measurement)
    fingerprints = {measurement["output_sha256"] for measurement in measurements
                    if measurement["implementation"] in {"reference", "cached", "store"}}
    if len(fingerprints) > 1 or (args.expected_output_sha256 and fingerprints != {args.expected_output_sha256}):
        raise RuntimeError("Conserved scientific outputs differ")
    summaries = {}
    for implementation in args.implementations:
        rows = [measurement for measurement in measurements if measurement["implementation"] == implementation]
        summaries[implementation] = {"wall_seconds_median": statistics.median(row["wall_seconds"] for row in rows),
                                     "wall_seconds_min": min(row["wall_seconds"] for row in rows),
                                     "wall_seconds_max": max(row["wall_seconds"] for row in rows),
                                     "cpu_seconds_median": statistics.median(row["cpu_seconds"] for row in rows),
                                     "peak_rss_bytes_max": max(row["peak_rss_bytes"] for row in rows),
                                     "metrics": rows[0]["metrics"]}
    print(json.dumps({"schema": 1, "workload": {"families": args.families, "species": args.species,
                                               "isoforms_per_locus": 4, "protein_length": args.length},
                      "environment": {"python": sys.version, "platform": platform.platform(),
                                      "implementation_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                                                                for path in (Path(__file__),
                                                                             Path(__file__).resolve().parents[1] / "support/gene_model_selection.py",
                                                                             Path(__file__).resolve().parents[1] / "support/gene_model_store.py")}},
                      "repeats": args.repeats, "warmups": args.warmups,
                      "input_sha256": measurements[0]["input_sha256"],
                      "reference_cached_outputs_equivalent": True if {"reference", "cached"}.issubset(summaries) else None,
                      "conserved_outputs_equivalent": True if len(set(args.implementations) - {"longest"}) > 1
                      or args.expected_output_sha256 else None,
                      "expected_output_sha256": args.expected_output_sha256,
                      "store_input_preparation": preparation,
                      "cache_speed_ratio": summaries["reference"]["wall_seconds_median"] /
                      max(1e-12, summaries["cached"]["wall_seconds_median"]) if {"reference", "cached"}.issubset(summaries) else None,
                      "summaries": summaries, "measurements": measurements}, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
