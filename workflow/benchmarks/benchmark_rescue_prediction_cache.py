#!/usr/bin/env python3
"""Compare equivalent JSON reads on an owned, bounded rescue-result sample."""
import argparse
import hashlib
import json
import os
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from rescue_prediction_cache import stream_json_array  # noqa: E402


def file_hash(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def worker(mode, path):
    start = time.perf_counter()
    if mode == "baseline":
        with path.open() as handle:
            rows = json.load(handle)
    else:
        rows = stream_json_array(path)
    result, count = hashlib.sha256(), 0
    for row in rows:
        result.update(json.dumps(row, sort_keys=True, separators=(",", ":"), allow_nan=False).encode())
        count += 1
    return {"seconds": time.perf_counter() - start, "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            "records": count, "semantic_sha256": result.hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--records", type=int, default=20000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--worker", choices=("baseline", "stream"))
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args.worker, args.input)))
        return
    if args.records < 1 or args.repeats < 1 or not args.output:
        parser.error("positive --records/--repeats and --output are required")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="prediction-json-benchmark-", dir=args.output.parent) as scratch:
        sample = Path(scratch) / "models.json"
        count = 0
        with sample.open("w") as handle:
            handle.write("[")
            for row in stream_json_array(args.input):
                handle.write(("," if count else "") + json.dumps(row, sort_keys=True, allow_nan=False))
                count += 1
                if count == args.records:
                    break
            handle.write("]\n")
        results = {"baseline": [], "stream": []}
        for repetition in range(args.repeats + 1):
            for mode in ("baseline", "stream"):
                command = [sys.executable, __file__, "--worker", mode, "--input", str(sample)]
                value = json.loads(subprocess.check_output(command, text=True, env={**os.environ, "PYTHONHASHSEED": "0"}))
                if repetition:
                    results[mode].append(value)
        all_results = [row for values in results.values() for row in values]
        equivalent = len({(row["records"], row["semantic_sha256"]) for row in all_results}) == 1
        if not equivalent:
            raise RuntimeError("Prediction JSON readers produced different records")
        report = {"input": str(args.input.resolve()), "input_sha256": file_hash(args.input),
                  "sample_records": count, "sample_bytes": sample.stat().st_size, "sample_sha256": file_hash(sample),
                  "python": sys.version, "platform": platform.platform(), "repeats": args.repeats, "warmups": 1,
                  "results": results, "equivalent": equivalent,
                  "summary": {mode: {"median_seconds": statistics.median(row["seconds"] for row in values),
                                     "median_peak_rss_kib": statistics.median(row["peak_rss_kib"] for row in values)}
                              for mode, values in results.items()}}
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report["summary"], indent=2))


if __name__ == "__main__":
    main()
