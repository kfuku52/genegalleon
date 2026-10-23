#!/usr/bin/env python3
"""Compare resampling time, memory, results and consumed random state."""

import argparse
import hashlib
import importlib.util
import json
import platform
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path

SUPPORT = Path(__file__).resolve().parents[1] / "support"
WORKLOADS = {
    "mc-small": ("mc", 17, 100000),
    "mc-medium": ("mc", 1001, 100000),
    "mc-large": ("mc", 10001, 100000),
    "inference": ("inference", 1001, 2000),
    "bootstrap": ("bootstrap", 10001, 2000),
}


def worker(args):
    import numpy as np

    sys.path.insert(0, str(SUPPORT))
    spec = importlib.util.spec_from_file_location("resampling_under_test", args.source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    mode, blocks, replicates = WORKLOADS[args.worker]
    sums = np.linspace(-0.73, 1.11, blocks)
    counts = np.arange(blocks) % 5 + 1
    rng = np.random.default_rng(931)
    started, cpu_started = time.perf_counter(), time.process_time()
    if mode == "mc":
        result = module.sign_flip(sums, rng, replicates, 1)
    elif mode == "bootstrap":
        result = module.bootstrap_estimates(sums, counts, replicates, rng).tolist()
    else:
        result = module.inference(sums, list(range(blocks)), replicates, rng,
                                  permutation_replicates=100000, exact_max_states=1)
    seconds, cpu_seconds = time.perf_counter() - started, time.process_time() - cpu_started
    payload = json.dumps([result, rng.bit_generator.state], sort_keys=True).encode()
    return {
        "seconds": seconds, "cpu_seconds": cpu_seconds,
        "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss *
        (1 if sys.platform == "darwin" else 1024),
        "sha256": hashlib.sha256(payload).hexdigest(), "numpy": np.__version__,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--implementation", type=Path, default=SUPPORT / "subgenome_dominance.py")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--cases", nargs="+", choices=WORKLOADS, default=list(WORKLOADS))
    parser.add_argument("--worker", choices=WORKLOADS)
    parser.add_argument("--source", type=Path)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args)))
        return
    if not args.baseline or args.repeats < 1:
        parser.error("--baseline and a positive --repeats are required")
    sources = {"baseline": args.baseline.resolve(), "current": args.implementation.resolve()}
    report = {
        "environment": {"python": sys.version, "platform": platform.platform()},
        "sources": {name: {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                    for name, path in sources.items()},
        "warmups": 1, "repeats": args.repeats, "samples": [], "medians": {},
    }
    for case in args.cases:
        collected = {name: [] for name in sources}
        for trial in range(-1, args.repeats):
            pair = {}
            order = list(sources) if trial % 2 else list(reversed(sources))
            for name in order:
                command = [sys.executable, str(Path(__file__).resolve()), "--worker", case,
                           "--source", str(sources[name])]
                sample = json.loads(subprocess.check_output(command, text=True))
                sample.update(case=case, implementation=name, trial=trial, warmup=trial == -1)
                pair[name] = sample
                report["samples"].append(sample)
                if trial >= 0:
                    collected[name].append(sample)
                print(json.dumps(sample), flush=True)
            if pair["baseline"]["sha256"] != pair["current"]["sha256"]:
                raise RuntimeError(f"Output or random-state mismatch: {case}")
        report["medians"][case] = {
            name: {metric: statistics.median(row[metric] for row in samples)
                   for metric in ("seconds", "cpu_seconds", "peak_rss_bytes")}
            for name, samples in collected.items()
        }
        if args.output:
            args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["medians"], indent=2))


if __name__ == "__main__":
    main()
