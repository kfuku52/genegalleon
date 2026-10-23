#!/usr/bin/env python3
"""Compare BLAST coverage aggregation with an explicit saved baseline.

Use workflow/tests/run_in_runtime.sh. Each timed sample is a fresh process;
fixture construction, TSV hashing and DataFrame contract hashing are untimed.
Peak RSS covers the whole worker, including fixture and proof construction.
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
import time
from pathlib import Path

PROFILE_SHAPES = {"singleton": (20000, (1,)), "mixed": (12500, (1, 1, 1, 1, 4)), "multihit": (500, (40,))}


def load(path):
    spec = importlib.util.spec_from_file_location("coverage", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fixture(profile, scale):
    import numpy
    import pandas

    groups, sizes = PROFILE_SHAPES[profile]
    sizes = numpy.resize(sizes, groups * scale)
    group = numpy.repeat(numpy.arange(len(sizes)), sizes)
    row = numpy.arange(len(group))
    start = (row * 23) % 800 + 1
    end = start + (row % 97) + 5
    reverse = row % 3 == 0
    qlen = numpy.full(len(row), 1000.0)
    slen = numpy.full(len(row), 2000.0)
    qlen[row % 997 == 0] = 0
    qlen[row % 991 == 0] = -1
    qlen[row % 983 == 0] = numpy.nan
    slen[row % 977 == 0] = numpy.inf
    return pandas.DataFrame({
        "qacc": [f"q{value:07d}" for value in group],
        "sacc": [f"s{value % 251:04d}" for value in group],
        "pident": numpy.where(row % 17 == 0, numpy.nan, 85.0),
        "length": end - start + 1,
        "mismatch": row % 8,
        "gapopen": row % 3,
        "qstart": numpy.where(reverse, end, start),
        "qend": numpy.where(reverse, start, end),
        "sstart": numpy.where(reverse, end + 30, start + 30),
        "send": numpy.where(reverse, start + 30, end + 30),
        "evalue": numpy.where(row % 31 == 0, "invalid", "1e-20"),
        "bitscore": numpy.where(row % 29 == 0, numpy.nan, 100.0 + row % 200),
        "frames": numpy.where(reverse, "0/-1", "0/1"),
        "qlen": qlen,
        "slen": slen,
    })


def contract_fingerprint(frame):
    proof = hashlib.sha256()
    metadata = {
        "columns": frame.columns.tolist(), "dtypes": [str(value) for value in frame.dtypes],
        "index": frame.index.tolist(), "index_type": type(frame.index).__qualname__,
        "index_names": frame.index.names, "column_names": frame.columns.names,
    }
    proof.update(json.dumps(metadata, sort_keys=True).encode())
    for column in frame:
        values = frame[column].to_numpy()
        if values.dtype.hasobject:
            typed_values = [(type(value).__module__, type(value).__qualname__, repr(value)) for value in values]
            proof.update(json.dumps(typed_values).encode())
        else:
            proof.update(values.tobytes())
    return proof.hexdigest()


def worker(args):
    import numpy

    module = load(args.script)
    frame = fixture(args.profile, args.scale)
    with numpy.errstate(divide="ignore", invalid="ignore", over="ignore"):
        started = time.perf_counter()
        started_cpu = time.process_time()
        output = module.annotate_blast_coverage(frame)
        seconds = time.perf_counter() - started
        cpu_seconds = time.process_time() - started_cpu
    print(json.dumps({
        "seconds": seconds, "cpu_seconds": cpu_seconds, "input_rows": len(frame), "output_rows": len(output),
        "tsv_sha256": hashlib.sha256(output.to_csv(sep="\t", index=False).encode()).hexdigest(),
        "contract_sha256": contract_fingerprint(output),
        "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024),
    }))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--implementation", type=Path, default=Path(__file__).resolve().parents[1] / "support" / "annotate_blast_coverage.py")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--scale", type=int, default=1)
    parser.add_argument("--profiles", nargs="+", choices=PROFILE_SHAPES, default=list(PROFILE_SHAPES))
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--script", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--profile", choices=PROFILE_SHAPES, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    if args.baseline is None or args.output is None:
        parser.error("--baseline and --output are required")
    if not 1 <= args.repeats <= 20 or not 1 <= args.scale <= 20:
        parser.error("--repeats and --scale must each be 1..20")
    import numpy
    import pandas

    results = {
        "python": sys.version, "platform": platform.platform(), "pandas": pandas.__version__, "numpy": numpy.__version__,
        "warmup": 1, "repeats": args.repeats, "scale": args.scale,
        "cache": "warm OS caches; fresh process per sample", "profiles": args.profiles,
        "scripts": {label: {"path": str(getattr(args, label)), "sha256": hashlib.sha256(getattr(args, label).read_bytes()).hexdigest()}
                    for label in ("baseline", "implementation")},
        "samples": {},
    }
    for profile in args.profiles:
        samples = {"baseline": [], "implementation": []}
        for trial in range(args.repeats + 1):
            labels = ["baseline", "implementation"] if trial % 2 == 0 else ["implementation", "baseline"]
            for label in labels:
                command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--script", str(getattr(args, label).resolve()),
                           "--profile", profile, "--scale", str(args.scale)]
                sample = json.loads(subprocess.check_output(command, text=True))
                print(profile, label, "warmup" if trial == 0 else trial, sample, flush=True)
                if trial:
                    samples[label].append(sample)
        for metric in ("tsv_sha256", "contract_sha256"):
            assert len({sample[metric] for group in samples.values() for sample in group}) == 1, (profile, metric)
        results["samples"][profile] = samples
    results["medians"] = {
        profile: {label: {metric: statistics.median(sample[metric] for sample in samples)
                         for metric in ("seconds", "cpu_seconds", "peak_rss_bytes")}
                  for label, samples in groups.items()}
        for profile, groups in results["samples"].items()
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results["medians"], indent=2))


if __name__ == "__main__":
    main()
