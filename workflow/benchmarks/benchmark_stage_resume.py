#!/usr/bin/env python3
"""Compare native checkpoint imports, including exact outputs and full-read I/O."""

import argparse
import contextlib
import hashlib
import io
import json
import platform
import resource
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def fixture(support, root, genome_mib):
    sys.path.insert(0, str(support))
    import input_generation_stage_resume as resume
    import run_input_generation_task as runner

    raw = root / "raw"
    raw.mkdir(exist_ok=True)
    roles = {key: raw / name for key, name in (("cds_path", "Example_species.cds.fa"),
             ("gff_path", "Example_species.gff"), ("genome_path", "Example_species.genome.fa"))}
    roles["cds_path"].write_text(">gene1\nATGAAA\n")
    roles["gff_path"].write_text("##gff-version 3\nchr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=gene1\n")
    with roles["genome_path"].open("wb") as handle:
        handle.write(b">chr1\n")
        chunk = (b"ATG" * 1024 + b"\n") * 341
        for _ in range(genome_mib):
            handle.write(chunk)
    plans = []
    for name in ("source", "target"):
        workspace = root / name
        output = workspace / "output/input_generation"
        plan = output / "tmp/task_plan.json"
        plan.parent.mkdir(parents=True, exist_ok=True)
        task = {"species_prefix": "Example_species", "species_key": "Example_species", "provider": "direct",
                "gene_grouping_mode": "rescue_overlap", "gff_repair_mode": "safe", "format_strict": True,
                "genetic_code": 1, **{k: str(v) for k, v in roles.items()}}
        task["input_sha256"] = resume.digest_paths(roles.values())
        resume.atomic_json(plan, {"version": 2, "task_count": 1, "tasks": [task]})
        settings = dict(provider="direct", gene_grouping_mode="rescue_overlap", gff_repair_mode="safe",
                        strict="1", run_validate_inputs="1", genetic_code="1", run_cds_fx2tab="0")
        settings.update({"species_" + role + "_dir": str(workspace / "input" / role)
                         for role in ("cds", "gff", "genome")})
        resume.atomic_json(str(plan) + ".settings.json", settings)
        resume.claim_workspace(plan, output, create=True)
        resume.atomic_json(str(plan) + ".prepared.json", {"plan_sha256": resume.digest(plan),
                           "settings_sha256": resume.digest(str(plan) + ".settings.json"), "files": {}})
        plans.append((plan, output, settings))
        if name != "source":
            continue
        meta = output / "tmp/task_meta_shards/1.json"
        sys.argv = ["runner", "--task-plan", str(plan), "--task-index", "1", "--describe-only",
                    "--task-meta-output", str(meta)]
        for key in ("species_cds_dir", "species_gff_dir", "species_genome_dir"):
            sys.argv.extend(["--" + key.replace("_", "-"), settings[key]])
        with contextlib.redirect_stdout(io.StringIO()):
            assert runner.main() == 0
        paths = resume.context(plan, 1, output, "format")[3]
        for label, source in (("cds", roles["cds_path"]), ("gff", roles["gff_path"]), ("genome", roles["genome_path"])):
            Path(paths[label]).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, paths[label])
        Path(paths["stats"]).parent.mkdir(parents=True, exist_ok=True)
        Path(paths["stats"]).write_text('{"fixture": true}\n')
        Path(paths["summary"]).parent.mkdir(parents=True, exist_ok=True)
        Path(paths["summary"]).write_text("species_prefix\tcds_output_path\nExample_species\t" + paths["cds"] + "\n")
        for name in ("1.longest.json", "1.mapping.json"):
            (output / "tmp/task_stats_shards" / name).write_text('{"fixture": true}\n')
        for stage in ("format", "validate"):
            resume.record(plan, 1, output, stage, "24")
        resume.atomic_json(str(plan) + ".completed/1.json", {"plan_sha256": resume.digest(plan), "task_index": 1,
                           "files": {item["path"]: item["sha256"] for item in resume.snapshot(paths).values()}})
    return resume, plans


def worker(args):
    resume, (source, target) = fixture(args.support_root, args.fixture_root, args.genome_mib)
    import input_generation_array_state as state
    original = state.digest
    reads = []
    original_copy = resume.copy_atomic
    copies = []

    def counted(path):
        reads.append(Path(path).stat().st_size)
        return original(path)

    state.digest = resume.digest = counted
    def counted_copy(source, destination, **kwargs):
        size = Path(source).stat().st_size
        copies.append(size)
        if kwargs.get("expected_sha256"):
            reads.append(size)  # The streaming copy also computes a SHA-256.
        return original_copy(source, destination, **kwargs)
    resume.copy_atomic = counted_copy
    started = time.perf_counter()
    with contextlib.redirect_stdout(io.StringIO()):
        resume.import_stages(argparse.Namespace(task_plan=target[0], root=target[1], source_plan=source[0],
                             source_root=source[1], source_plan_sha256=original(source[0]),
                             format_contract_version="24", task_index=None))
    elapsed = time.perf_counter() - started
    proof = hashlib.sha256()
    for path in sorted((args.fixture_root / "target").rglob("*")):
        if path.is_file() and ".namespace-v1" not in str(path) and not path.name.endswith(".lock"):
            proof.update(str(path.relative_to(args.fixture_root)).encode())
            with path.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    proof.update(chunk)
    print(json.dumps({"seconds": elapsed, "sha256_calls": len(reads), "sha256_bytes": sum(reads),
                      "output_sha256": proof.hexdigest(),
                      "copy_bytes": sum(copies),
                      "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                      * (1 if sys.platform == "darwin" else 1024),
                      "children_peak_rss_bytes": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
                      * (1 if sys.platform == "darwin" else 1024)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-root", type=Path, default=Path(__file__).resolve().parents[1] / "support")
    parser.add_argument("--baseline-support", type=Path)
    parser.add_argument("--genome-mib", type=int, default=256)
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--fixture-root", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not 1 <= args.genome_mib <= 4096 or not 1 <= args.trials <= 10:
        parser.error("Require 1..4096 genome MiB and 1..10 trials")
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error("--output is required")
    roots = {"current": args.support_root.resolve()}
    if args.baseline_support:
        roots["baseline"] = args.baseline_support.resolve()
    result = {"python": sys.version, "platform": platform.platform(), "genome_mib": args.genome_mib,
              "cache_mode": "warm; no filesystem cache eviction", "samples": {}}
    with tempfile.TemporaryDirectory(prefix="gg-resume-benchmark-") as directory:
        root = Path(directory)
        for trial in range(args.trials + 1):
            for label, support in list(roots.items())[::(-1 if trial % 2 else 1)]:
                for name in ("source", "target"):
                    shutil.rmtree(root / name, ignore_errors=True)
                command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--support-root", str(support),
                           "--fixture-root", str(root), "--genome-mib", str(args.genome_mib)]
                sample = json.loads(subprocess.check_output(command, text=True))
                print(label, "warmup" if trial == 0 else trial, sample, flush=True)
                if trial:
                    result["samples"].setdefault(label, []).append(sample)
    assert len({sample["output_sha256"] for samples in result["samples"].values() for sample in samples}) == 1
    result["medians"] = {label: {key: statistics.median(sample[key] for sample in samples)
                                for key in ("seconds", "sha256_calls", "sha256_bytes", "copy_bytes",
                                            "peak_rss_bytes", "children_peak_rss_bytes")}
                         for label, samples in result["samples"].items()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(result["medians"])


if __name__ == "__main__":
    main()
