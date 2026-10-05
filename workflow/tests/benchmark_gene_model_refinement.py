"""Repeatable complete-pipeline benchmark using held-out models and real miniprot.

Run in the GeneGalleon runtime, for example:
  python workflow/tests/benchmark_gene_model_refinement.py --repeats 3 --cpus 2

Fixture preparation and container startup are excluded. Fresh interpreters run
the actual CLI; one batched invocation is compared with separate stage CLIs.
Longest/off is a different-behavior resource baseline, never a speedup oracle.
Cold means a new output directory, not flushed operating-system disk caches.
Peak RSS is the maximum individual process RSS, not a sum of concurrent workers.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import importlib.metadata
import json
import os
import platform
import resource
import runpy
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

TESTS = Path(__file__).resolve().parent
SUPPORT = TESTS.parent / "support"
SCRIPT = SUPPORT / "gene_model_refinement.py"
STAGES = ("catalog", "correspondence", "select", "predict", "finalize")
POLICIES = {"longest_off": ("longest", "off"), "conserved_conservative": ("conserved", "conservative")}


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def table(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def peak_rss():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    own = int(value if sys.platform == "darwin" else value * 1024)
    status = Path("/proc/self/status")
    if status.exists():
        for line in status.read_text().splitlines():
            if line.startswith("VmHWM:"):
                own = int(line.split()[1]) * 1024
                break
    child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    child = int(child if sys.platform == "darwin" else child * 1024)
    return max(own, child), own, child


def worker(command, metrics_path):
    """Keep fixture imports/preparation outside each measured CLI process."""
    sys.path.insert(0, str(SUPPORT))
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    before_cpu, before_wall = time.process_time(), time.perf_counter()
    sys.argv = [str(SCRIPT), *command]
    runpy.run_path(str(SCRIPT), run_name="__main__")
    wall = time.perf_counter() - before_wall
    after_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = time.process_time() - before_cpu + after_children.ru_utime + after_children.ru_stime \
        - before_children.ru_utime - before_children.ru_stime
    maximum, own, children = peak_rss()
    Path(metrics_path).write_text(json.dumps({"cli_wall_seconds": wall, "cli_cpu_seconds": cpu,
                                             "peak_rss_bytes": maximum, "cli_peak_rss_bytes": own,
                                             "children_peak_rss_bytes": children}, sort_keys=True) + "\n")


def measure(command, scratch, index):
    """Wait for one fresh child, accounting for CLI and predictor descendants."""
    metric_path = scratch / f"{index}.metrics.json"
    argv = [sys.executable, str(Path(__file__).resolve()), "--worker-command", json.dumps(command),
            "--worker-metrics", str(metric_path)]
    with tempfile.TemporaryFile(mode="w+") as stdout, tempfile.TemporaryFile(mode="w+") as stderr:
        before = time.perf_counter()
        process = subprocess.Popen(argv, stdout=stdout, stderr=stderr, text=True)
        _pid, status, usage = os.wait4(process.pid, 0)
        process.returncode = os.waitstatus_to_exitcode(status)
        wall = time.perf_counter() - before
        stdout.seek(0)
        stderr.seek(0)
        if process.returncode:
            raise RuntimeError(f"CLI failed ({process.returncode}): {command}\n"
                               + stdout.read()[-8000:] + stderr.read()[-8000:])
    result = json.loads(metric_path.read_text())
    return {"command": [sys.executable, str(SCRIPT), *command], "wall_seconds": wall,
            "cpu_seconds": usage.ru_utime + usage.ru_stime, **result}


def scientific_outputs(root):
    """Compare exported bytes and path-independent scientific manifest fields."""
    effective = root / "effective"
    roles = ("cds", "protein", "gff", "genome", "representative_map")
    manifest = [{key: row[key] for key in ("species", "genetic_code", *(role + "_sha256" for role in roles))}
                for row in table(effective / "inputs.tsv")]
    files = {str(path.relative_to(effective)): sha256(path) for path in sorted(effective.rglob("*"))
             if path.is_file() and (path.suffix in {".fa", ".gff3"}
                                   or path.name in {"representative_map.tsv", "species_genetic_code.tsv", "changes.json"})}
    predictions = {}
    for path in sorted((root / "predictions").glob("*/predictions.json")):
        predictions[path.parent.name] = json.loads(path.read_text())
    value = {"manifest": manifest, "files": files, "predictions": predictions}
    return {"sha256": fingerprint(value), "manifest_sha256": fingerprint(manifest), "files": files,
            "prediction_sha256": fingerprint(predictions)}


def truth_counts(root, truth, original_records):
    sys.path.insert(0, str(SUPPORT))
    fasta_records = importlib.import_module("fasta_sequence_store").fasta_records
    target = "Species_target"
    selected = {identifier.removeprefix(target + "_"): sequence for identifier, _header, sequence in
                fasta_records(root / "effective/species_cds" / (target + ".fa"))}
    predictions = json.loads((root / "predictions" / target / "predictions.json").read_text())
    accepted = [row for row in predictions if row["status"] == "accepted"]
    negatives = ("intact", "pseudogene", "stop", "gap")
    choices = [row for row in table(root / "effective/representative_map.tsv") if row["species"] == target]
    predicted_ids = {row["candidate"]["candidate_id"] for row in accepted}
    exact_revisions = sum(row["change_type"] == "model_revision"
                          and row["candidate"]["cds"] == truth[target, row["gene_id"].removeprefix(target + "_")]["cds"]
                          and row["candidate"]["blocks"] == truth[target, row["gene_id"].removeprefix(target + "_")]["blocks"]
                          for row in accepted)
    alternative = [row for row in accepted if row["gene_id"] == target + "_alternate"]
    return {"heldout_exon_loci": 2,
            "heldout_exon_loci_reconstructed_exactly": sum(selected[kind] == truth[target, kind]["cds"]
                                                          for kind in ("repair", "minus")),
            "accepted_exact_model_revisions": exact_revisions,
            "accepted_isoform_additions": sum(row["change_type"] == "isoform_addition" for row in accepted),
            "accepted_exon_skipped_path_matches_truth": bool(alternative and alternative[0]["candidate"]["cds"]
                                                             == truth[target, "alternate"]["shorter_cds"]),
            "biological_negative_or_intact_loci": len(negatives),
            "negative_or_intact_loci_changed": sum(selected[kind] != original_records[kind + ".long"] for kind in negatives),
            "negative_or_intact_predictions_accepted": sum(row["gene_id"].removeprefix(target + "_") in negatives
                                                           for row in accepted),
            "selected_predicted_representatives": sum(row["candidate_id"] in predicted_ids for row in choices),
            "selected_exon_skipped_representative": selected["alternate"] == truth[target, "alternate"]["shorter_cds"]}


def summarize(rows):
    return {"repeats": len(rows), "wall_seconds_median": statistics.median(row["wall_seconds"] for row in rows),
            "wall_seconds_min": min(row["wall_seconds"] for row in rows),
            "wall_seconds_max": max(row["wall_seconds"] for row in rows),
            "cpu_seconds_median": statistics.median(row["cpu_seconds"] for row in rows),
            "peak_rss_bytes_max": max(row["peak_rss_bytes"] for row in rows)}


def benchmark(args):
    sys.path.insert(0, str(TESTS))
    fixture = importlib.import_module("test_gene_model_refinement_runtime")
    measurements, fingerprints, counts = [], {}, {}
    with tempfile.TemporaryDirectory(prefix="gg-refinement-benchmark-") as temporary:
        temporary = Path(temporary)
        data = temporary / "inputs"
        data.mkdir()
        inputs, edges, rna, sources, truth = fixture.truth_fixture(data, rna=not args.without_rna)
        source_files = [inputs, edges, *([rna] if rna else []),
                        *(Path(row[role]) for row in sources for role in ("cds", "gff", "genome"))]
        input_hashes = {path.name: sha256(path) for path in source_files}
        semantic_inputs = {"sources": [{"species": row["species"], "genetic_code": row["genetic_code"],
                                         **{role + "_sha256": sha256(row[role]) for role in ("cds", "gff", "genome")}}
                                        for row in sources], "edges_sha256": sha256(edges),
                           "rna_sha256": sha256(rna) if rna else None}
        original_records = {identifier: sequence for identifier, _header, sequence in
                            fixture.fasta_records(Path(next(row["cds"] for row in sources if row["species"] == "Species_target")))}
        measurement_index = 0
        for policy, (representative_policy, refinement_mode) in POLICIES.items():
            for execution in ("batched", "staged"):
                for repetition in range(args.repeats):
                    root = temporary / f"{policy}-{execution}-{repetition}"
                    plan = ["plan", "--inputs", str(inputs), "--edges", str(edges), "--output", str(root),
                            "--policy", representative_policy, "--mode", refinement_mode,
                            "--padding", "650", "--max-intron", "2000"]
                    if rna:
                        plan += ["--rna", str(rna)]
                    for cache in ("cold", "resume"):
                        commands = ([plan] if cache == "cold" else []) + [
                            [command, "--output", str(root), "--cpus", str(args.cpus)]
                            for command in (("run",) if execution == "batched" else STAGES)]
                        stages = []
                        for command in commands:
                            measurement_index += 1
                            stages.append({"stage": command[0], **measure(command, temporary, measurement_index)})
                        signature = scientific_outputs(root)
                        summary = truth_counts(root, truth, original_records)
                        if policy in fingerprints and signature["sha256"] != fingerprints[policy]:
                            raise RuntimeError(f"Scientific outputs differ across execution styles/restarts for {policy}")
                        fingerprints[policy] = signature["sha256"]
                        if summary["negative_or_intact_loci_changed"] or summary["negative_or_intact_predictions_accepted"]:
                            raise RuntimeError("Benchmark damaged a known intact model or biological negative")
                        if refinement_mode == "conservative" and (summary["heldout_exon_loci_reconstructed_exactly"] != 2
                                                                 or not summary["accepted_exon_skipped_path_matches_truth"]):
                            raise RuntimeError("Real miniprot benchmark did not reconstruct held-out truth")
                        counts[policy] = summary
                        measurements.append({"policy": policy, "execution": execution, "cache": cache,
                                             "repetition": repetition, "wall_seconds": sum(row["wall_seconds"] for row in stages),
                                             "cpu_seconds": sum(row["cpu_seconds"] for row in stages),
                                             "peak_rss_bytes": max(row["peak_rss_bytes"] for row in stages),
                                             "output_bytes": sum(path.stat().st_size for path in root.rglob("*") if path.is_file()),
                                             "scientific_outputs": signature, "truth_counts": summary, "stages": stages})
        if input_hashes != {path.name: sha256(path) for path in source_files}:
            raise RuntimeError("Benchmark modified source inputs")
    summaries = {}
    for policy in POLICIES:
        summaries[policy] = {}
        for execution in ("batched", "staged"):
            for cache in ("cold", "resume"):
                rows = [row for row in measurements if (row["policy"], row["execution"], row["cache"])
                        == (policy, execution, cache)]
                key = execution + "_" + cache
                summaries[policy][key] = summarize(rows)
                if execution == "staged":
                    summaries[policy][key]["stages"] = {
                        stage: summarize([item for row in rows for item in row["stages"] if item["stage"] == stage])
                        for stage in (("plan",) + STAGES if cache == "cold" else STAGES)}
    return {"schema": 1, "environment": {"python": sys.version, "platform": platform.platform(),
                                          "cpu_count": os.cpu_count(), "cpus_per_prediction": args.cpus,
                                          "miniprot": subprocess.check_output(["miniprot", "--version"], text=True).strip(),
                                          "dependencies": {name: importlib.metadata.version(name)
                                                           for name in ("biopython", "jcvi", "kfFractBias")},
                                          "implementation_sha256": {path.name: sha256(path) for path in
                                                                     (SCRIPT, SUPPORT / "gene_model_catalog.py",
                                                                      SUPPORT / "gene_model_selection.py", SUPPORT / "gene_model_store.py",
                                                                      Path(__file__))}},
            "workload": {"species": 3, "loci_per_species": 7, "protein_length": 280,
                         "rna_whole_path": not args.without_rna, "explicit_frozen_correspondence": True,
                         "synteny_discovery_measured": False, "fixture_sha256": sha256(Path(fixture.__file__)),
                         "input_sha256": fingerprint(semantic_inputs), "semantic_inputs": semantic_inputs,
                         "source_hashes": input_hashes},
            "measurement": {"repeats": args.repeats, "cold_definition": "new output directory; OS caches not flushed",
                            "resume_definition": "same completed output directory in a fresh interpreter",
                            "wall_includes": "fresh interpreter startup and all CLI commands; container startup excluded",
                            "cpu_includes": "CLI and completed predictor descendants",
                            "peak_rss_definition": "maximum individual CLI/child process RSS; not simultaneous aggregate",
                            "fixture_preparation_measured": False,
                            "policy_baseline_outputs_equivalent": False,
                            "same_policy_scientific_outputs_equivalent": True,
                            "staged_batched_comparison": "same scientific behavior; CLI startup/verification overhead differs"},
            "scientific_output_sha256": fingerprints, "truth_counts": counts,
            "summaries": summaries, "measurements": measurements}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--cpus", type=int, default=2)
    parser.add_argument("--without-rna", action="store_true")
    parser.add_argument("--output", type=Path, help="Save the complete report; stdout then contains compact summaries")
    parser.add_argument("--worker-command", help=argparse.SUPPRESS)
    parser.add_argument("--worker-metrics", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker_command:
        if not args.worker_metrics:
            parser.error("Worker requires a metrics path")
        worker(json.loads(args.worker_command), args.worker_metrics)
        return
    if args.repeats < 1 or args.cpus < 1:
        parser.error("repeats and cpus must be positive")
    if not shutil.which("miniprot"):
        parser.error("miniprot is required; run in the GeneGalleon runtime")
    result = benchmark(args)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
        result = {"report": str(args.output), "summaries": result["summaries"], "truth_counts": result["truth_counts"],
                  "scientific_output_sha256": result["scientific_output_sha256"], "measurement": result["measurement"]}
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
