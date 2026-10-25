#!/usr/bin/env python3
"""Compare complete verified refinement review collection on synthetic publications.

Fixtures and source/output hashing are outside timers. Actual verify_inputs,
receipt hashes, locus-store build/read, and the complete collect call are used;
there are no proof or scientific-tool stubs. Imports are outside timers.
Each comparison uses one discarded warmup and alternating fresh worker processes.
Run in the GeneGalleon runtime:

python workflow/benchmarks/benchmark_refinement_review.py \\
    --baseline /mounted/baseline/plot_gene_model_refinement.py
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import platform
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow/support"
sys.path.insert(0, str(SUPPORT))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def write_table(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def receipt(directory, plan_hash, members):
    write_json(directory / "receipt.json", {"key": {"plan": plan_hash},
               "files": {str(path.relative_to(directory)): sha(path) for path in members}})


def fixture(root, species_count, genes):
    from gene_model_store import build_store

    root.mkdir(parents=True)
    write_json(root / "plan.json", {"benchmark": "synthetic-refinement-review", "schema": 1})
    plan_hash = sha(root / "plan.json")
    effective = root / "effective"
    effective.mkdir()
    repmap = effective / "representative_map.tsv"
    repmap.write_text("species\tgene_id\tcandidate_id\n")
    rows, changes, catalog_dirs = [], [], []
    names = [f"Species_s{index:03d}" for index in range(species_count)]
    for name in names:
        roles = {}
        for role, extension in [("cds", "fa"), ("protein", "fa"), ("gff", "gff3"), ("genome", "fa")]:
            path = effective / ("species_" + role) / (name + "." + extension)
            path.parent.mkdir(exist_ok=True)
            with path.open("w") as handle:
                if role == "genome":
                    handle.write(">chr1\n" + "ATGAAACCCTAA" * genes + "\n")
                elif role == "gff":
                    handle.write("##gff-version 3\n")
                    for index in range(genes):
                        handle.write(f"chr1\ts\tCDS\t{12 * index + 1}\t{12 * index + 12}\t.\t+\t0\tID={name}_t{index};Parent={name}_g{index}\n")
                else:
                    for index in range(genes):
                        handle.write(f">{name}_g{index}\n" + ("ATGAAACCCTAA" if role == "cds" else "MKP*") + "\n")
            roles[role] = str(path)
        roles["representative_map"] = str(repmap)
        rows.append({"species": name, "genetic_code": 1, **roles,
                     **{key + "_sha256": sha(value) for key, value in roles.items()}})
        directory = root / "catalog" / name
        directory.mkdir(parents=True)
        metadata = {"schema": 1, "species": name, "summary": {"loci": genes, "candidates": genes,
                    "fasta_mapping": {"mapped": genes, "unmapped": 0}}}
        write_json(directory / "catalog_metadata.json", metadata)
        models = []
        with (directory / "loci.jsonl").open("w") as handle:
            for index in range(genes):
                candidate = {"candidate_id": f"{name}_t{index}", "source_transcript_id": f"t{index}",
                             "origin": "original", "blocks": [[12 * index, 12 * index + 12, 0]],
                             "cds": "ATGAAACCCTAA", "protein": "MKP", "quality": {"usable": True}}
                locus = {"species": name, "gene_id": f"{name}_g{index}", "seqid": "chr1", "strand": "+",
                         "source_baseline_candidate_id": candidate["candidate_id"], "candidates": [candidate]}
                handle.write(json.dumps(locus, separators=(",", ":")) + "\n")
                if index % 50 == 0:
                    models.append({"gene_id": locus["gene_id"], "status": "accepted" if index % 100 == 0 else "proposal",
                                   "change_type": "model_revision", "problems": [] if index % 100 == 0 else ["weak_support"],
                                   "candidate": candidate, "donors": [], "rna_paths": [], "alignments": []})
        receipt(directory, plan_hash, [directory / "catalog_metadata.json", directory / "loci.jsonl"])
        catalog_dirs.append(directory)
        predictions = root / "predictions" / name
        write_json(predictions / "predictions.json", models)
        receipt(predictions, plan_hash, [predictions / "predictions.json"])
    # Interleave species and retain duplicate selection records. Summary counts
    # must see every record while the detailed selection uses the last duplicate.
    for index in range(genes):
        for name in names:
            changes.append({"species": name, "gene_id": f"{name}_g{index}", "candidate_id": f"{name}_t{index}",
                            "status": "conserved" if index % 41 == 0 else "retained",
                            "selected_origin": "original"})
    changes.append({**changes[0], "status": "sequence_repaired"})
    write_json(effective / "changes.json", changes)
    write_table(effective / "inputs.tsv", list(rows[0]), rows)
    write_table(effective / "effective_exclusions.tsv", ["species", "gene_id"], [])
    write_table(effective / "translation_admission.tsv", ["species", "gene_id", "status", "quality"],
                [{"species": name, "gene_id": f"{name}_g0", "status": "included",
                  "quality": json.dumps({"phase_inferred": True})} for name in names])
    receipt(effective, plan_hash, [path for path in effective.rglob("*") if path.is_file()])
    db = root / "catalog_index_final" / "loci.sqlite3"
    build_store(catalog_dirs, db)
    receipt(db.parent, plan_hash, [db])


def worker(args):
    spec = importlib.util.spec_from_file_location("measured_review", args.source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    start, cpu = time.perf_counter(), time.process_time()
    try:
        result = module.collect(args.fixture, max_loci=8, preferred_species="Species_s001")
        error = None
    except Exception as exc:
        result, error = None, [type(exc).__name__, str(exc)]
    elapsed, cpu = time.perf_counter() - start, time.process_time() - cpu
    peak = next(int(line.split()[1]) * 1024 for line in Path("/proc/self/status").read_text().splitlines()
                if line.startswith("VmHWM:"))
    raw = json.dumps(result, sort_keys=False, separators=(",", ":"), allow_nan=False).encode()
    args.payload.write_bytes(raw)
    print(json.dumps({"wall_s": elapsed, "cpu_s": cpu, "peak_rss_bytes": peak,
                      "output_sha256": hashlib.sha256(raw).hexdigest(), "error": error,
                      "species_count": len(result["species"]) if result else 0,
                      "details_count": len(result["details"]) if result else 0}))


def controller(args, work):
    shapes = {"small": (4, 32), "medium": (32, 500), "large": (128, 500)}
    samples, summary = [], {}
    for profile in args.profiles:
        root = work / profile
        fixture(root, *shapes[profile])
        for repeat in range(args.repeats + 1):
            pair = {}
            for label in (["baseline", "candidate"] if repeat % 2 == 0 else ["candidate", "baseline"]):
                payload = work / (label + ".json")
                cmd = [sys.executable, str(Path(__file__).resolve()), "--worker", "--baseline", str(args.baseline),
                       "--source", str(getattr(args, label)), "--fixture", str(root), "--payload", str(payload)]
                output = subprocess.run(cmd, check=True, capture_output=True, text=True,
                                        env=dict(os.environ, PYTHONHASHSEED="0", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1"))
                pair[label] = json.loads(output.stdout)
                assert pair[label]["error"] is None, pair[label]
                assert pair[label]["species_count"] == shapes[profile][0], pair[label]
                assert pair[label]["details_count"] == min(8, shapes[profile][0]), pair[label]
            assert pair["baseline"]["output_sha256"] == pair["candidate"]["output_sha256"], pair
            assert (work / "baseline.json").read_bytes() == (work / "candidate.json").read_bytes()
            if repeat:
                samples.append({"profile": profile, "pair": repeat, **pair})
        subset = [sample for sample in samples if sample["profile"] == profile]
        summary[profile] = {label: {metric: statistics.median(row[label][metric] for row in subset)
                                    for metric in ("wall_s", "cpu_s", "peak_rss_bytes")}
                            for label in ("baseline", "candidate")}
        summary[profile]["paired_speedups"] = [row["baseline"]["wall_s"] / row["candidate"]["wall_s"] for row in subset]
        print(profile, summary[profile], flush=True)
        # Corruption must still fail at the same real proof boundary.
        target = root / "predictions" / "Species_s000" / "predictions.json"
        original = target.read_bytes()
        target.write_bytes(original + b" ")
        errors = []
        for label in ("baseline", "candidate"):
            cmd = [sys.executable, str(Path(__file__).resolve()), "--worker", "--baseline", str(args.baseline),
                   "--source", str(getattr(args, label)), "--fixture", str(root),
                   "--payload", str(work / (label + "-error.json"))]
            output = subprocess.run(cmd, check=True, capture_output=True, text=True)
            errors.append(json.loads(output.stdout)["error"])
        assert errors[0] == errors[1] and errors[0][0] == "ValueError", errors
        target.write_bytes(original)
    result = {"environment": {"platform": platform.platform(), "python": sys.version}, "repeats": args.repeats,
              "source_sha256": {label: sha(getattr(args, label)) for label in ("baseline", "candidate")},
              "summaries": summary, "measurements": samples,
              "equivalence": "complete collect JSON bytes, insertion order, bounded gallery, actual tamper rejection"}
    args.output.write_text(json.dumps(result, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True, help="Saved baseline plot_gene_model_refinement.py.")
    parser.add_argument("--candidate", type=Path, default=SUPPORT / "plot_gene_model_refinement.py")
    parser.add_argument("--profiles", nargs="+", choices=["small", "medium", "large"], default=["small", "medium", "large"])
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--output", type=Path, default=Path("tmp/refinement-review-benchmark.json"))
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--source", type=Path)
    parser.add_argument("--fixture", type=Path)
    parser.add_argument("--payload", type=Path)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    if args.worker:
        worker(args)
        return
    args.baseline = args.baseline.resolve(strict=True)
    args.candidate = args.candidate.resolve(strict=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="gg-refinement-review-benchmark-") as temporary:
        controller(args, Path(temporary))


if __name__ == "__main__":
    main()
