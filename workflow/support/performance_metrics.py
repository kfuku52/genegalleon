"""Optional process-local timing and I/O counters; never workflow authority."""

import argparse
import atexit
import contextlib
import json
import os
import resource
import socket
import sys
import threading
import time
from pathlib import Path

_started = time.perf_counter()
_counters = {}
_lock = threading.Lock()


def count(name, amount=1):
    if os.environ.get("GG_PERFORMANCE_DIR"):
        with _lock:
            _counters[name] = _counters.get(name, 0) + amount


def emit(phase, started, counters, status):
    scale = 1 if sys.platform == "darwin" else 1024
    record = {"schema_version": 1, "phase": phase, "pid": os.getpid(),
              "seconds": time.perf_counter() - started, "status": status,
              "process_peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale,
              "children_peak_rss_bytes": resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * scale,
              "counters": counters, "counters_are_inclusive": True}
    write_row(record)


def write_row(record):
    directory = os.environ.get("GG_PERFORMANCE_DIR")
    if not directory:
        return
    record = {**record, "recorded_at_epoch": time.time(), "host": socket.gethostname(),
              "job_id": os.environ.get("GG_JOB_ID", os.environ.get("SLURM_JOB_ID", os.environ.get("JOB_ID", ""))),
              "array_task_id": os.environ.get("GG_ARRAY_TASK_ID", os.environ.get("SLURM_ARRAY_TASK_ID", ""))}
    try:
        root = Path(directory)
        root.mkdir(parents=True, exist_ok=True)
        with (root / f"{os.getpid()}.jsonl").open("a") as handle:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    except OSError as exc:
        # Diagnostics cannot invalidate scientific outputs or certify success.
        print(f"Warning: performance telemetry unavailable: {exc}", file=sys.stderr)


@contextlib.contextmanager
def measure(phase):
    started = time.perf_counter()
    with _lock:
        before = dict(_counters)
    status = "failed"
    try:
        yield
        status = "ok"
    finally:
        with _lock:
            delta = {key: value - before.get(key, 0) for key, value in _counters.items()}
            emit(phase, started, delta, status)


def _finish():
    emit("process", _started, dict(_counters), "exit")


atexit.register(_finish)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("start", "elapsed"))
    parser.add_argument("--phase")
    parser.add_argument("--started", type=float)
    parser.add_argument("--status", choices=("ok", "failed"), default="ok")
    args = parser.parse_args()
    if args.action == "start":
        print(time.perf_counter())
    else:
        if args.phase is None or args.started is None:
            parser.error("Elapsed requires --phase and --started")
        write_row({"schema_version": 1, "phase": args.phase, "pid": os.getppid(),
                   "seconds": time.perf_counter() - args.started, "status": args.status, "wall_only": True})


if __name__ == "__main__":
    main()
