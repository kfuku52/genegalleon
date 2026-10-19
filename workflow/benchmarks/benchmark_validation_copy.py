#!/usr/bin/env python3
"""Compare equivalent copies of real validation overlays from an owned JSONL.

Run through the runtime wrapper. Each trial uses a fresh process so maximum
RSS remains comparable; source parsing and output verification are not timed.
"""
import argparse
import copy
import hashlib
import json
import os
import platform
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from rescue_raw_validation import copy_validation_data  # noqa: E402


def signature(rows):
    digest = hashlib.sha256()
    for row in rows:
        # Include object field order as well as values.
        digest.update(json.dumps(row, ensure_ascii=False, allow_nan=False, separators=(",", ":")).encode())
        digest.update(b"\n")
    return digest.hexdigest()


def worker(args):
    digest, rows = hashlib.sha256(), []
    with args.input.open("rb") as stream:
        for line in stream:
            digest.update(line)
            rows.append(json.loads(line))
    if digest.hexdigest() != args.input_sha256:
        raise ValueError("Validation-copy benchmark input changed")
    expected = signature(rows)
    operation = copy.deepcopy if args.worker == "baseline" else copy_validation_data
    cpu, start = time.process_time(), time.perf_counter()
    copied = [operation(row) for row in rows]
    seconds, cpu_seconds = time.perf_counter() - start, time.process_time() - cpu
    scale = 1 / 1024 if sys.platform == "darwin" else 1
    maximum = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale
    actual = signature(copied)
    if actual != expected:
        raise RuntimeError("Validation copy changed data or field order")
    return dict(records=len(rows), seconds=seconds, cpu_seconds=cpu_seconds, peak_rss_kib=maximum,
                output_sha256=actual)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--worker", choices=("baseline", "candidate"))
    parser.add_argument("--input-sha256")
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args)))
        return
    if not args.output or args.repeats < 1:
        parser.error("--output and positive --repeats are required")
    digest = hashlib.sha256()
    with args.input.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    results = {"baseline": [], "candidate": []}
    for repeat in range(args.repeats + 1):
        for mode in (("baseline", "candidate") if repeat % 2 == 0 else ("candidate", "baseline")):
            command = [sys.executable, __file__, "--worker", mode, "--input", str(args.input),
                       "--input-sha256", digest.hexdigest()]
            row = json.loads(subprocess.check_output(command, text=True,
                                                    env={**os.environ, "PYTHONHASHSEED": "0"}))
            if repeat:
                results[mode].append(row)
    if len({(r["records"], r["output_sha256"]) for rows in results.values() for r in rows}) != 1:
        raise RuntimeError("Validation-copy trials disagree")
    report = dict(input=str(args.input.resolve()), input_sha256=digest.hexdigest(),
                  implementation_sha256=hashlib.sha256((Path(__file__).resolve().parents[1] / "support" /
                                                       "rescue_raw_validation.py").read_bytes()).hexdigest(),
                  python=sys.version, platform=platform.platform(), warmups=1, repeats=args.repeats,
                  equal=True, results=results,
                  summary={mode: {"median_seconds": statistics.median(r["seconds"] for r in rows),
                                  "median_CPU_seconds": statistics.median(r["cpu_seconds"] for r in rows),
                                  "median_peak_rss_kib": statistics.median(r["peak_rss_kib"] for r in rows)}
                           for mode, rows in results.items()})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["summary"]))


if __name__ == "__main__":
    main()
