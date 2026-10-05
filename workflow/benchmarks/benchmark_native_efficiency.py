#!/usr/bin/env python3
"""Compare fresh native verification and genome formatting in isolated fixtures."""
import argparse
import contextlib
import gzip
import hashlib
import io
import json
import os
import platform
import random
import resource
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def tree_rss(pid):
    try:
        status = Path(f"/proc/{pid}/status").read_text()
        rss = int(next(line.split()[1] for line in status.splitlines() if line.startswith("VmRSS:")))
        children = Path(f"/proc/{pid}/task/{pid}/children").read_text().split()
        return rss + sum(tree_rss(int(child)) for child in children)
    except (OSError, StopIteration, ValueError):
        return 0


def fingerprint(paths):
    result = hashlib.sha256()
    for path in paths:
        opener = gzip.open if path.name.endswith(".gz") else open
        with opener(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                result.update(chunk)
    return result.hexdigest()


def worker(args):
    sys.path.insert(0, str(args.repo / "workflow/support"))
    import input_generation_array_state as state
    root = args.fixture_root
    root.mkdir(parents=True)
    reads = []
    original = state.digest
    def digest(path):
        reads.append(Path(path).stat().st_size)
        return original(path)
    state.digest = digest
    os.environ["GG_TASK_CPUS"] = "1"
    if args.mode == "prepared":
        plan = root / "plan.json"
        state.atomic_json(plan, {"task_count": args.species, "tasks": [
            {"species_prefix": f"Species_{index}"} for index in range(args.species)]})
        state.atomic_json(str(plan) + ".settings.json", {})
        staging = Path(str(plan) + ".tasks")
        staging.mkdir()
        files = {}
        for index in range(1, args.species + 1):
            for suffix in (".json", ".resolved.tsv"):
                path = staging / f"{index}{suffix}"
                path.write_text(str(index) + "x" * 1536)
                files[str(path)] = original(path)
        state.atomic_json(str(plan) + ".prepared.json", {
            "plan_sha256": original(plan), "settings_sha256": original(str(plan) + ".settings.json"), "files": files})
        assert state.prepared(plan)  # The global prepare/finalize contract stays intact.
        reads.clear()
        stamp = time.perf_counter()
        import inspect
        options = {"task_index": 1} if "task_index" in inspect.signature(state.prepared).parameters else {}
        assert state.prepared(plan, **options)
        seconds = time.perf_counter() - stamp
        proof = original(plan)
    elif args.mode == "genome":
        from format_species_discovery import format_genome
        assert shutil.which("seqkit"), "Require real seqkit from the qualified runtime"
        path = root / "Example_species.genome.fa"
        fragment = "".join(random.Random(82431).choices("ACGT", k=4096))
        with path.open("w") as handle:
            handle.write(">chr1 original header\n")
            for _ in range(args.genome_mib * 256):
                handle.write(fragment + "\n")
        stamp = time.perf_counter()
        output = format_genome({"species_prefix": "Example_species", "genome_path": path}, root, True, False)
        seconds = time.perf_counter() - stamp
        proof = fingerprint([output["output_path"]])
    else:
        # Shared fixture utilities drive each variant's own core/support files,
        # so the baseline QC has a baseline implementation identity as well.
        sys.path.insert(0, str(ROOT / "workflow/tests"))
        import test_gg_input_generation_end_to_end as fixture
        fixture.REPO_ROOT = args.repo
        fixture.CORE_PATH = args.repo / "workflow/core/gg_input_generation_core.sh"
        inputs = fixture._write_direct_species_fixture(root)
        for path in inputs.glob("*/*.genome.fa"):
            with path.open("a") as handle:
                for _ in range(args.genome_mib * 1024):
                    handle.write("ACGT" * 256 + "\n")
        tools = fixture._install_fake_toolchain(root)
        source = root / "source"
        fixture._write_minimal_ete_taxonomy_db(source)
        fixture._write_runtime_busco_dataset(source)
        fixture._run_core(source, inputs, tools, "array_prepare")
        for index in ((1, 2) if args.mode == "final" else (1,)):
            fixture._run_core(source, inputs, tools, "array_worker", index)
        native = source / "output/input_generation"
        if args.mode == "final":
            from format_species_summary import read_species_summary_rows, write_species_summary_rows
            rows = {}
            for path in (native / "tmp/species_summary_shards").glob("*.tsv"):
                rows.update(read_species_summary_rows(path))
            summary = native / "summary.tsv"
            write_species_summary_rows(summary, rows)
            import input_validation_reuse as reuse
            import validate_cds_gff_mapping as mapping
            import validate_longest_cds_selection as ownership
            common = ["--species-cds-dir", str(native / "species_cds"), "--species-summary", str(summary),
                      "--reuse-validation-root", str(native), "--reuse-task-plan", str(native / "tmp/task_plan.json"),
                      "--format-contract-version", "24"]
            specific = ["--species-gff-dir", str(native / "species_gff"), "--species-genome-dir", str(native / "species_genome")]
            output_paths = [root / "mapping.json", root / "ownership.json"]
            reads.clear()
            reuse.implementation_identity.cache_clear()
            stdout = io.StringIO()
            stamp = time.perf_counter()
            with contextlib.redirect_stdout(stdout):
                if hasattr(reuse, "VerificationSession"):
                    sys.argv = ["final-qc", *common, *specific, "--mapping-stats-output", str(output_paths[0]),
                                "--ownership-stats-output", str(output_paths[1])]
                    assert reuse.main() == 0
                else:
                    for module, path, options in ((mapping, output_paths[0], specific), (ownership, output_paths[1], [])):
                        reuse.implementation_identity.cache_clear()
                        sys.argv = ["validator", *common, *options, "--stats-output", str(path)]
                        assert module.main() == 0
            seconds = time.perf_counter() - stamp
            assert stdout.getvalue().count("Reused verified") == 4
            payloads = [json.loads(path.read_text()) for path in output_paths]
            for payload in payloads:
                payload.pop("validation_implementation")  # Deliberately changed code identity, not scientific output.
            proof = hashlib.sha256(json.dumps(payloads, sort_keys=True).encode()).hexdigest()
        else:
            plan = native / "tmp/task_plan.json"
            target = root / "target"
            fixture._write_minimal_ete_taxonomy_db(target)
            fixture._write_runtime_busco_dataset(target, "embryophyta_odb12")
            common = dict(overwrite="0", busco_lineage="embryophyta_odb12", resume_from_task_plan=str(plan),
                          resume_from_task_plan_sha256=original(plan), resume_from_input_generation_root=str(native))
            env = fixture._core_env(target, inputs, tools, "array_prepare")
            env.update(common)
            subprocess.run(["bash", str(fixture.CORE_PATH)], env=env, capture_output=True, check=True)
            env = fixture._core_env(target, inputs, tools, "array_worker", 1)
            env.update(common)
            stamp = time.perf_counter()
            result = subprocess.run(["bash", str(fixture.CORE_PATH)], env=env, capture_output=True, text=True, check=True)
            seconds = time.perf_counter() - stamp
            assert "Reused verified input formatting" in result.stdout and "Reused verified CDS/GFF validation" in result.stdout
            generated = target / "output/input_generation"
            data_sha = fingerprint([path for role in ("species_cds", "species_gff", "species_genome",
                                    "species_cds_fx2tab", "species_cds_busco_full", "species_cds_busco_short")
                                    for path in sorted((generated / role).glob("*")) if path.is_file()])
            payloads = [json.loads((generated / f"tmp/task_stats_shards/1.{label}.json").read_text())
                        for label in ("mapping", "longest")]
            for payload in payloads:
                payload.pop("validation_implementation")
            proof = hashlib.sha256(json.dumps([data_sha, payloads], sort_keys=True).encode()).hexdigest()
    print(json.dumps({"seconds": seconds, "logical_sha256": proof, "sha256_calls": len(reads),
                      "sha256_bytes": sum(reads), "process_peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("genome", "prepared", "final", "worker"), required=True)
    parser.add_argument("--baseline-repo", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--genome-mib", type=int, default=64)
    parser.add_argument("--species", type=int, default=560)
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument("--repo", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--fixture-root", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if min(args.genome_mib, args.species, args.trials) < 1:
        parser.error("Positive workload sizes and trials required")
    if args.repo:
        worker(args)
        return
    if not args.baseline_repo or args.output.exists():
        parser.error("Require a baseline repository and a new output path")
    roots = {"baseline": args.baseline_repo.resolve(strict=True), "current": ROOT}
    samples = {key: [] for key in roots}
    with tempfile.TemporaryDirectory(prefix="gg-native-efficiency-") as temporary:
        fixture = Path(temporary) / "fixture"
        for trial in range(args.trials + 1):
            for name in list(roots)[::(-1 if trial % 2 else 1)]:
                shutil.rmtree(fixture, ignore_errors=True)
                command = [sys.executable, str(Path(__file__).resolve()), "--mode", args.mode, "--repo", str(roots[name]),
                           "--output", str(args.output), "--fixture-root", str(fixture), "--genome-mib", str(args.genome_mib),
                           "--species", str(args.species)]
                proc = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
                peak = 0
                while proc.poll() is None:
                    if args.mode == "genome":
                        peak = max(peak, tree_rss(proc.pid))
                    time.sleep(0.01)
                stdout = proc.communicate()[0]
                if proc.returncode:
                    raise RuntimeError(f"{name} worker failed: {proc.returncode}")
                sample = json.loads(stdout)
                if args.mode == "genome":
                    sample["sampled_tree_peak_rss_kib"] = peak
                print(name, trial, sample, flush=True)
                if trial:
                    samples[name].append(sample)
    if len({row["logical_sha256"] for rows in samples.values() for row in rows}) != 1:
        raise RuntimeError("Outputs differ")
    payload = {"mode": args.mode, "genome_mib": args.genome_mib, "species": args.species,
               "platform": platform.platform(), "python": sys.version, "warmups": 1, "trials": samples,
               "medians": {name: {key: statistics.median(row[key] for row in rows) for key in rows[0] if key != "logical_sha256"}
                           for name, rows in samples.items()}, "equivalent": True,
               "scope": "Warm filesystem, alternating fresh processes; synthetic fixtures only. Genome uses real seqkit, "
                        "one thread, process-tree RSS sampled at 10 ms. Final QC timing excludes fixture setup and module imports, "
                        "includes full scientific species_results, and recomputes baseline identity per validator. "
                        "Worker timing uses real core/validators but fake BUSCO/seqkit/Rscript and local taxonomy. "
                        "Worker SHA counters exclude its subprocesses; do not use them as whole-worker I/O. "
                        "Peak process RSS for final/worker includes setup; not isolated phase memory. Not an HPC/NAS benchmark."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload["medians"]))


if __name__ == "__main__":
    main()
