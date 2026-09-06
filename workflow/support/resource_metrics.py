#!/usr/bin/env python3
"""Measure a command without changing its output, status, or scientific settings.

Run inside the GeneGalleon container. Records are per attempt, never appended
concurrently to one shared file. ru_maxrss is the OS-reported per-process high
water mark, not the simultaneous sum of all parallel workers.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
from pathlib import Path
import signal
import subprocess
import sys
import time
import uuid


def read_events(path):
    """Treat every malformed observation as unavailable, never a job failure."""
    try:
        events = [json.loads(line) for line in path.read_text().splitlines()]
        for row in events:
            if (not isinstance(row, dict) or not isinstance(row.get("stage"), str)
                    or row.get("kind") not in {"start", "skip"}
                    or any(type(row.get(key)) not in (int, float)
                           or not math.isfinite(row[key]) or row[key] < 0 for key in ("at", "cpu"))):
                return []
        return events
    except (ValueError, OSError):
        return []


def run(command, directory, *, workflow, input_sha256, cpus, memory_gb,
        runtime_id, stage="workflow", server_id="local", expected_stages=()):
    if cpus < 1 or memory_gb <= 0 or not command:
        raise ValueError("positive CPU/memory budgets and a command are required")
    directory = Path(directory)
    writable = True
    try:
        directory.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        writable = False
        print(f"Resource metrics unavailable: {exc}", file=sys.stderr)
    start = time.monotonic()
    started_at = time.time()
    event_path = directory / ("events-" + uuid.uuid4().hex + ".jsonl")
    times_path = event_path.with_suffix(".times")
    env = dict(os.environ)
    if writable:
        env.update(GG_RESOURCE_EVENT_FILE=str(event_path), GG_RESOURCE_TIMES_FILE=str(times_path))
    else:
        env.pop("GG_RESOURCE_EVENT_FILE", None)
        env.pop("GG_RESOURCE_TIMES_FILE", None)
    child = subprocess.Popen(command, start_new_session=True, env=env)
    previous = {}

    def forward(signum, _frame):
        try:
            os.killpg(child.pid, signum)
        except ProcessLookupError:
            pass

    try:
        for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            previous[sig] = signal.signal(sig, forward)
        _, status, usage = os.wait4(child.pid, 0)
        child.returncode = os.waitstatus_to_exitcode(status)
    finally:
        for sig, handler in previous.items():
            signal.signal(sig, handler)
    wall = time.monotonic() - start
    record = {"schema_version": 1, "workflow": workflow, "stage": stage,
              "input_sha256": input_sha256, "runtime_id": runtime_id,
              "server_id": server_id, "cpus": cpus, "memory_gb": memory_gb,
              "started_at": started_at, "wall_seconds": wall,
              "cpu_seconds": usage.ru_utime + usage.ru_stime,
              "peak_process_rss_gb": usage.ru_maxrss / (1024**3 if sys.platform == "darwin" else 1024**2),
              "exit_code": child.returncode}
    record["cpu_efficiency"] = record["cpu_seconds"] / (wall * cpus)
    observations = read_events(event_path)
    events = [row for row in observations if row["kind"] == "start"]
    record["expected_stages"] = sorted(expected_stages)
    record["skipped_stages"] = [row["stage"] for row in observations if row["kind"] == "skip"]
    record["learning_eligible"] = bool(expected_stages) and (
        sorted(row["stage"] for row in events) == sorted(expected_stages)
        and not set(record["skipped_stages"]).intersection(expected_stages))
    boundaries = []
    for before, after in zip(events, events[1:]):
        boundaries.append({"stage": before["stage"],
                           "wall_seconds": max(0, after["at"] - before["at"]),
                           "cpu_seconds_reaped": max(0, after["cpu"] - before["cpu"])})
    if events:
        boundaries.append({"stage": events[-1]["stage"],
                           "wall_seconds": max(0, time.monotonic() - events[-1]["at"]),
                           "cpu_seconds_reaped": max(0, record["cpu_seconds"] - events[-1]["cpu"])})
    record["stage_boundaries"] = boundaries
    try:
        times_path.unlink(missing_ok=True)
    except OSError:
        pass
    target = directory / ("attempt-" + uuid.uuid4().hex + ".json")
    temporary = target.with_suffix(".tmp")
    try:
        with temporary.open("x") as stream:
            json.dump(record, stream, sort_keys=True, allow_nan=False)
            stream.write("\n")
        temporary.replace(target)
        if child.returncode == 0 and record["learning_eligible"] and input_sha256 and runtime_id != "unidentified" and server_id != "unidentified":
            key = {"workflow": workflow, "input_sha256": input_sha256,
                   "runtime_id": runtime_id, "server_id": server_id, "stage": stage}
            key_hash = hashlib.sha256(json.dumps(key, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
            index_dir = directory / "history" / key_hash
            index_dir.mkdir(parents=True, exist_ok=True)
            indexed = index_dir / f"cpus-{cpus}.json"
            pending = index_dir / (uuid.uuid4().hex + ".tmp")
            pending.write_text(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
            pending.replace(indexed)
    except OSError as exc:
        print(f"Resource metrics could not be saved: {exc}", file=sys.stderr)
    return child.returncode if child.returncode >= 0 else 128 - child.returncode


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--profile", type=Path)
    parser.add_argument("--workflow")
    parser.add_argument("--runtime-id", required=True)
    parser.add_argument("--server-id", required=True)
    parser.add_argument("--stage", default="workflow")
    parser.add_argument("--cpus", type=int, default=int(os.environ.get("GG_TASK_CPUS", "1")))
    parser.add_argument("--memory-gb", type=float, default=float(os.environ.get("GG_MEM_TOTAL_GB", "1")))
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    profile_path = args.profile or (Path(os.environ["GG_RESOURCE_PROFILE"]) if os.environ.get("GG_RESOURCE_PROFILE") else None)
    value = {"workflow": args.workflow, "input_sha256": None, "stages": []}
    if profile_path:
        try:
            candidate = json.loads(profile_path.read_text())
            if (not isinstance(candidate, dict)
                    or not isinstance(candidate.get("workflow"), str)
                    or (args.workflow and candidate["workflow"] != args.workflow)
                    or not isinstance(candidate.get("input_sha256"), str)
                    or re.fullmatch(r"[0-9a-f]{64}", candidate["input_sha256"]) is None
                    or not isinstance(candidate.get("stages"), list)
                    or any(not isinstance(s, str) or not s for s in candidate["stages"])):
                raise ValueError("invalid resource profile shape or workflow")
            value = candidate
        except (OSError, ValueError) as exc:
            if args.profile:
                raise
            print(f"Resource profile unavailable; recording unclassified timing: {exc}", file=sys.stderr)
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    return run(command, args.directory, workflow=value["workflow"], input_sha256=value["input_sha256"],
               cpus=args.cpus, memory_gb=args.memory_gb, runtime_id=args.runtime_id,
               server_id=args.server_id, stage=args.stage, expected_stages=value["stages"])


def event(stage, kind="start"):
    target = os.environ.get("GG_RESOURCE_EVENT_FILE")
    times_file = os.environ.get("GG_RESOURCE_TIMES_FILE")
    if not target or not times_file:
        return
    raw = Path(times_file).read_text()
    values = re.findall(r"(\d+)m([0-9.]+)s", raw)
    cpu = sum(int(minutes) * 60 + float(seconds) for minutes, seconds in values)
    with open(target, "a") as handle:
        handle.write(json.dumps({"stage": stage, "kind": kind, "at": time.monotonic(), "cpu": cpu}) + "\n")


if __name__ == "__main__":
    if sys.argv[1:2] == ["event"]:
        event(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "start")
    else:
        raise SystemExit(main())
