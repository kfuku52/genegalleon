#!/usr/bin/env python3
"""Measure closed-log file reduction and packing cost with exact body/metadata checks.

This synthetic fixture compares raw storage with the new batch ZIP, not two
scientific implementations. It does not predict shared-filesystem throughput.
"""

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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow/support"))
import directory_output_store as store  # noqa: E402
import orthofinder_output_store as backend  # noqa: E402


def measure(args):
    with tempfile.TemporaryDirectory(prefix="gg-scheduler-logs-") as temporary:
        base = Path(temporary).resolve()
        root = base / "batch"
        logs = root / "scheduler_output"
        logs.mkdir(parents=True)
        expected = {}
        for task in range(args.tasks):
            for stream in ("out", "err"):
                seed = hashlib.sha256(f"{task}.{stream}".encode()).hexdigest().encode()
                body = (seed * ((args.bytes_per_log + 63) // 64))[:args.bytes_per_log]
                path = logs / f"gene.{task}.{stream}"
                path.write_bytes(body)
                path.chmod(0o640)
                os.utime(path, ns=(1700000000123456789, 1700000000123456789))
                expected[path.name] = (hashlib.sha256(body).hexdigest(), len(body), 0o640, path.stat().st_mtime_ns)
        def counts():
            entries = list(base.rglob("*"))
            return {"regular_files": sum(p.is_file() for p in entries),
                    "directories": sum(p.is_dir() for p in entries)}
        before = counts()
        began = time.perf_counter()
        store.convert(root, kind="scheduler_logs", writers_stopped=True)
        elapsed = time.perf_counter() - began
        after = counts()
        records = backend.archive_records(root, "scheduler_output", deep=True)["files"]
        observed = {row["path"]: (row["sha256"], row["size"], row["mode"], row["mtime_ns"]) for row in records}
        if observed != expected:
            raise RuntimeError("Stored scheduler log content or metadata differs")
        # Exercise the public logical reader without materializing all logs.
        with store.open_binary(root, "scheduler_output/gene.0.out", kind="scheduler_logs") as handle:
            if backend._hash(handle)[:2] != expected["gene.0.out"][:2]:
                raise RuntimeError("Logical log reader differs")
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024)
        return {"seconds": elapsed, "peak_rss_bytes": rss, "before": before, "after": after,
                "payload_files_before": len(expected), "payload_files_after": 1,
                "zip_bytes": (root / "scheduler_output.zip").stat().st_size,
                "logical_bytes": sum(row[1] for row in expected.values()), "equivalent": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tasks", type=int, default=1000)
    parser.add_argument("--bytes-per-log", type=int, default=4096)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if min(args.tasks, args.bytes_per_log, args.repeats) < 1:
        parser.error("workload sizes must be positive")
    if args.worker:
        print(json.dumps(measure(args)))
        return
    rows = []
    for _ in range(args.repeats + 1):
        command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--tasks", str(args.tasks),
                   "--bytes-per-log", str(args.bytes_per_log)]
        rows.append(json.loads(subprocess.check_output(command, text=True)))
    measured = rows[1:]
    print(json.dumps({"environment": {"python": sys.version, "platform": platform.platform()},
                      "workload": {"tasks": args.tasks, "logs_per_task": 2, "bytes_per_log": args.bytes_per_log,
                                   "discarded_warmups": 1, "repeats": args.repeats},
                      "samples": measured, "equivalent": all(row["equivalent"] for row in measured),
                      "median_seconds": statistics.median(row["seconds"] for row in measured),
                      "median_peak_rss_bytes": statistics.median(row["peak_rss_bytes"] for row in measured)}, indent=2))


if __name__ == "__main__":
    main()
