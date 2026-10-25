#!/usr/bin/env python3
"""Compare BH helper memory in fresh processes with exact float-bit outputs.

Fixtures are generated and saved by the parent, outside measured workers. Each
sample includes helper import and input loading in peak RSS; timings cover only
calculate_bh_fdr. Peak RSS is captured before result hashing. Run without other
benchmarks, tests, or builds in the same GeneGalleon runtime.
"""

import argparse
import hashlib
import importlib.util
import json
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

SUPPORT = Path(__file__).resolve().parents[1] / "support"
PROFILES = ("dense", "sparse", "ties", "small")


def file_sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def bit_sha256(values):
    # A byte view avoids a full-size temporary that could distort worker RSS.
    return hashlib.sha256(memoryview(values.view("uint64")).cast("B")).hexdigest()


def peak_rss():
    if sys.platform.startswith("linux"):
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith("VmHWM:"):
                return int(line.split()[1]) * 1024, "/proc/self/status:VmHWM"
        raise RuntimeError("Linux process status does not report VmHWM")
    scale = 1 if sys.platform == "darwin" else 1024
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale, "getrusage:ru_maxrss"


def worker(args):
    import numpy as np

    sys.path.insert(0, str(SUPPORT))
    spec = importlib.util.spec_from_file_location("benchmark_bh_source", args.source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    values = np.load(args.input, allow_pickle=False)
    before = bit_sha256(values)
    started = time.perf_counter()
    started_cpu = time.process_time()
    result = module.calculate_bh_fdr(values)
    cpu_seconds = time.process_time() - started_cpu
    seconds = time.perf_counter() - started
    rss_bytes, rss_source = peak_rss()
    if result.shape != values.shape or result.dtype != np.dtype("float64"):
        raise AssertionError("BH output shape or dtype changed")
    if bit_sha256(values) != before:
        raise AssertionError("BH helper modified its input")
    print(json.dumps({
        "seconds": seconds, "cpu_seconds": cpu_seconds,
        "peak_rss_bytes": rss_bytes, "peak_rss_source": rss_source,
        "input_uint64_sha256": before, "q_uint64_sha256": bit_sha256(result),
        "rows": len(values), "input_unchanged": True,
    }))


def make_fixture(profile, rows):
    import numpy as np

    rng = np.random.default_rng(20261025)
    if profile == "small":
        return np.resize(np.array([
            0.0, -0.0, np.nextafter(0.0, 1.0), 1e-300, 0.01, 0.01, 0.04, 1.0, np.nan,
        ], dtype=np.float64), 4096)
    if profile == "ties":
        return rng.choice(np.array([
            0.0, -0.0, np.nextafter(0.0, 1.0), 1e-300, 0.0001, 0.01, 0.01, 0.2, 1.0, np.nan,
        ], dtype=np.float64), size=rows)
    values = rng.random(rows, dtype=np.float64)
    if profile == "sparse":
        for offset in range(1, 5):
            values[offset::5] = np.nan
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--implementation", type=Path, default=SUPPORT / "csubst_scan_support.py")
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--baseline-only", action="store_true",
                        help="Measure only --baseline (or --implementation when omitted).")
    parser.add_argument("--rows", type=int, default=5000000,
                        help="Rows in dense/sparse/ties fixtures; small always has 4096.")
    parser.add_argument("--profiles", nargs="+", choices=PROFILES, default=list(PROFILES))
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--source", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--input", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        if args.source is None or args.input is None:
            parser.error("Worker requires --source and --input")
        worker(args)
        return
    if args.output is None or (args.baseline is None and not args.baseline_only):
        parser.error("Require --output and either --baseline or --baseline-only")
    if args.rows < 1 or not 1 <= args.repeats <= 20:
        parser.error("Require --rows >= 1 and --repeats in 1..20")
    sources = {"baseline": (args.baseline or args.implementation).resolve()}
    if not args.baseline_only:
        sources["implementation"] = args.implementation.resolve()
    import numpy as np
    import pandas as pd

    report = {
        "python": sys.version, "platform": platform.platform(),
        "numpy": np.__version__, "pandas": pd.__version__,
        "command": sys.argv, "seed": 20261025, "warmup": 1, "repeats": args.repeats,
        "cache": "warm OS caches; fresh process per sample",
        "timing_scope": "calculate_bh_fdr only",
        "rss_scope": "worker including imports and input loading, captured before result hashing",
        "sources": {label: {"path": str(path), "sha256": file_sha256(path)}
                    for label, path in sources.items()},
        "inputs": {}, "warmups": {}, "samples": {}, "medians": {},
    }
    with tempfile.TemporaryDirectory(prefix="gg-bh-memory-") as temporary:
        root = Path(temporary)
        for profile in args.profiles:
            path = root / (profile + ".npy")
            values = make_fixture(profile, args.rows)
            np.save(path, values, allow_pickle=False)
            report["inputs"][profile] = {
                "rows": len(values), "finite_rows": int(np.isfinite(values).sum()),
                "dtype": str(values.dtype), "npy_bytes": path.stat().st_size,
                "npy_sha256": file_sha256(path), "uint64_sha256": bit_sha256(values),
            }
            del values
        for profile in args.profiles:
            samples = {label: [] for label in sources}
            warmups = {}
            proofs = set()
            for trial in range(args.repeats + 1):
                labels = list(sources)
                if trial % 2:
                    labels.reverse()
                for label in labels:
                    sample = json.loads(subprocess.check_output([
                        sys.executable, str(Path(__file__).resolve()), "--worker",
                        "--source", str(sources[label]), "--input", str(root / (profile + ".npy")),
                    ], text=True))
                    if sample["input_uint64_sha256"] != report["inputs"][profile]["uint64_sha256"]:
                        raise AssertionError(f"Worker input differs: {profile}/{label}")
                    proofs.add(sample["q_uint64_sha256"])
                    if len(proofs) != 1:
                        raise AssertionError(f"BH outputs differ in float bits: {profile}")
                    if trial:
                        samples[label].append(sample)
                    else:
                        warmups[label] = sample
                    print(profile, label, "warmup" if trial == 0 else trial, json.dumps(sample), flush=True)
            report["warmups"][profile] = warmups
            report["samples"][profile] = samples
            report["medians"][profile] = {
                label: {metric: statistics.median(sample[metric] for sample in rows)
                        for metric in ("seconds", "cpu_seconds", "peak_rss_bytes")}
                for label, rows in samples.items()
            }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["medians"], indent=2))


if __name__ == "__main__":
    main()
