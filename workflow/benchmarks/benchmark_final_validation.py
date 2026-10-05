#!/usr/bin/env python3
"""Compare complete native final QC with receipt-bound worker-QC reuse."""

import argparse
import contextlib
import hashlib
import io
import json
import os
import platform
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def prepare(root, genome_mib):
    sys.path.insert(0, str(ROOT / "workflow/tests"))
    import test_gg_input_generation_end_to_end as fixture

    inputs = fixture._write_direct_species_fixture(root)
    for genome in inputs.glob("*/*.genome.fa"):
        with genome.open("a") as handle:
            for _ in range(genome_mib * 1024):
                handle.write("ACGT" * 256 + "\n")
    workspace = root / "workspace"
    tools = fixture._install_fake_toolchain(root)
    fixture._write_minimal_ete_taxonomy_db(workspace)
    fixture._write_runtime_busco_dataset(workspace)
    fixture._run_core(workspace, inputs, tools, "array_prepare")
    for index in (1, 2):
        fixture._run_core(workspace, inputs, tools, "array_worker", index)
    native = workspace / "output/input_generation"
    sys.path.insert(0, str(ROOT / "workflow/support"))
    from format_species_summary import read_species_summary_rows, write_species_summary_rows
    rows = {}
    for path in sorted((native / "tmp/species_summary_shards").glob("*.tsv")):
        rows.update(read_species_summary_rows(path))
    write_species_summary_rows(native / "summary.tsv", rows)
    return native


def worker(args):
    native = args.worker_root
    # Match the native core's enabled per-process performance instrumentation.
    os.environ["GG_PERFORMANCE_DIR"] = str(native / "tmp/benchmark_performance" / str(os.getpid()))
    sys.path.insert(0, str(args.support))
    import input_generation_array_state as state
    import validate_cds_gff_mapping as mapping
    import validate_longest_cds_selection as longest
    from format_species_annotation import reference

    count = {"sha256_bytes": 0, "genome_reference_scans": 0}
    real_digest, real_index = state.digest, reference.genome_reference_index

    def digest(path):
        count["sha256_bytes"] += Path(path).stat().st_size
        return real_digest(path)

    def index(path):
        count["genome_reference_scans"] += 1
        return real_index(path)

    state.digest, reference.genome_reference_index = digest, index
    common = ["--species-cds-dir", str(native / "species_cds"), "--species-summary", str(native / "summary.tsv")]
    extra = []
    if args.reuse:
        import re
        core = (ROOT / "workflow/core/gg_input_generation_core.sh").read_text()
        contract = re.search(r'^format_contract_version="?(\d+)', core, re.MULTILINE)[1]
        extra = ["--reuse-validation-root", str(native), "--reuse-task-plan", str(native / "tmp/task_plan.json"),
                 "--format-contract-version", contract]
    output = io.StringIO()
    started = time.perf_counter()
    summaries = {}
    with contextlib.redirect_stdout(output):
        for name, module, specific in (
            ("mapping", mapping, ["--species-gff-dir", str(native / "species_gff"),
                                  "--species-genome-dir", str(native / "species_genome")]),
            ("ownership", longest, []),
        ):
            if args.reuse:
                # The core launches each validator in a separate process.
                from input_validation_reuse import implementation_identity
                implementation_identity.cache_clear()
            path = native / f"benchmark-{name}.json"
            sys.argv = [str(module.__file__), *common, *specific, *extra, "--stats-output", str(path)]
            if module.main():
                raise RuntimeError("Final validation failed")
            payload = json.loads(path.read_text())
            payload.pop("species_results", None)
            payload.pop("validation_options", None)
            payload.pop("validation_implementation", None)
            summaries[name] = payload
    result = {"seconds": time.perf_counter() - started,
              "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              "logical_sha256": hashlib.sha256(json.dumps(summaries, sort_keys=True).encode()).hexdigest(),
              "reused_species_stage_count": output.getvalue().count("Reused verified"), **count}
    print(json.dumps(result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-support", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--genome-mib", type=int, default=64)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--worker-root", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--support", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--reuse", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker_root:
        worker(args)
        return
    if not args.baseline_support or not args.output or args.output.exists():
        parser.error("--baseline-support and a new --output directory are required")
    if min(args.genome_mib, args.repeats) < 1:
        parser.error("Workload and repetitions must be positive")
    args.output = args.output.resolve()
    args.baseline_support = args.baseline_support.resolve(strict=True)
    args.output.mkdir(parents=True)
    native = prepare(args.output, args.genome_mib)
    variants = {"baseline": args.baseline_support, "current": ROOT / "workflow/support"}
    results = {name: [] for name in variants}
    for trial in range(args.repeats + 1):
        names = list(variants)
        if trial % 2:
            names.reverse()
        for name in names:
            command = [sys.executable, str(Path(__file__).resolve()), "--worker-root", str(native),
                       "--support", str(variants[name])]
            if name == "current":
                command.append("--reuse")
            result = json.loads(subprocess.check_output(command, text=True))
            if trial:
                results[name].append(result)
    if len({r["logical_sha256"] for trials in results.values() for r in trials}) != 1:
        raise RuntimeError("Final QC counts or source ownership differ")
    payload = {"species": 2, "genome_mib_each": args.genome_mib, "warmups": 1,
               "repeats": args.repeats, "equivalent": True, "trials": results,
               "python": sys.version, "platform": platform.platform(),
               "median_seconds": {name: statistics.median(r["seconds"] for r in trials)
                                  for name, trials in results.items()},
               "support_source_sha256": {
                   name: hashlib.sha256(json.dumps({str(path.relative_to(support)): hashlib.sha256(path.read_bytes()).hexdigest()
                                                   for path in sorted(support.rglob("*.py"))}, sort_keys=True).encode()).hexdigest()
                   for name, support in variants.items()},
               "scope": "Native worker-generated synthetic fixtures and proofs; real validators and GFF reader; "
                        "fake BUSCO/seqkit and local taxonomy only during fixture setup; "
                        "fresh process per trial; native performance instrumentation enabled; "
                        "implementation identity recomputed per validator; "
                        "fixture preparation/imports excluded; Linux RSS; "
                        "not a complete finalize or production NAS benchmark."}
    (args.output / "result.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload))


if __name__ == "__main__":
    main()
