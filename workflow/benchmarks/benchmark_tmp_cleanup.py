#!/usr/bin/env python3
"""Concurrent stale-temp housekeeping with retained-content equivalence proof."""

import argparse
import hashlib
import json
import os
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def manifest(root):
    rows = []
    for p in sorted(root.rglob("*")):
        if p.is_file():
            rows.append((p.relative_to(root).as_posix(), hashlib.sha256(p.read_bytes()).hexdigest()))
    return hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest(), len(rows)


def worker(args):
    sys.path.insert(0, str(Path(args.source) / "workflow/support"))
    from gene_family_output_store import cleanup_stale_tmp, family_bucket_lock

    if args.writer_probe and not args.writer_family:
        original_scandir = os.scandir
        original_rglob = Path.rglob
        marker = Path(args.gate).with_name("scanning")

        def notify(path):
            if isinstance(path, int):
                return
            current = Path(path)
            if current.parent == Path(args.worker_root).resolve() / "tmp" and not marker.exists():
                try:
                    with marker.open("x"):
                        pass
                except FileExistsError:
                    pass
                else:
                    time.sleep(0.1)

        def notify_scan(path):
            notify(path)
            return original_scandir(path)

        def notify_rglob(path, *a, **kw):
            notify(path)
            return original_rglob(path, *a, **kw)

        os.scandir = notify_scan
        Path.rglob = notify_rglob
    print("ready", flush=True)
    while not Path(args.gate).exists():
        time.sleep(0.005)
    if args.writer_family:
        deadline = time.monotonic() + 10
        while not Path(args.gate).with_name("scanning").exists():
            if time.monotonic() > deadline:
                raise RuntimeError("no cleanup reached directory measurement")
            time.sleep(0.005)
        start = time.perf_counter()
        with family_bucket_lock(Path(args.worker_root) / ".gg_store", args.writer_family, exclusive=False) as acquired:
            assert acquired
        removed = []
    else:
        start = time.perf_counter()
        removed = cleanup_stale_tmp(Path(args.worker_root), older_than_days=7, nonblocking=True)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(
        json.dumps(
            dict(
                role="writer" if args.writer_family else "cleaner",
                seconds=time.perf_counter() - start,
                removed=len(removed),
                max_rss_bytes=rss if sys.platform == "darwin" else rss * 1024,
            )
        ),
        flush=True,
    )


def one(args, root):
    tree = root / "orthogroup" / "tmp"
    tree.mkdir(parents=True)
    old = time.time_ns() - 10 * 86400 * 10**9
    for i in range(args.families):
        d = tree / (str(i + 1) + "_OG" + str(i).zfill(7))
        d.mkdir()
        for j in range(args.files):
            p = d / ("file" + str(j) + ".txt")
            p.write_text("synthetic archived input " + str(j) + "\n")
            if i < args.families * 3 // 4:
                os.utime(p, ns=(old, old))
        if i < args.families * 3 // 4:
            os.utime(d, ns=(old, old))
    gate = root / "start"
    children = []
    for index in range(args.workers):
        command = [
            sys.executable,
            __file__,
            "--source",
            args.source,
            "--worker-root",
            str(tree.parent),
            "--gate",
            str(gate),
        ]
        if args.writer_probe:
            command += ["--writer-probe"]
            if index == args.workers - 1:
                command += ["--writer-family", "OG" + str(args.families - 1).zfill(7)]
        children.append(subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True))
    for p in children:
        if p.stdout.readline().strip() != "ready":
            raise RuntimeError(p.stderr.read())
    start = time.perf_counter()
    gate.touch()
    records = []
    for p in children:
        out, err = p.communicate(timeout=180)
        if p.returncode:
            for child in children:
                if child.poll() is None:
                    child.terminate()
                child.communicate(timeout=10)
            raise RuntimeError(err)
        records.append(json.loads(out))
    # Nonblocking housekeeping may defer an eligible family under contention.
    # Include one quiescent pass in the measured workload so both versions
    # converge to the same retention contract before comparing outputs.
    sys.path.insert(0, str(Path(args.source) / "workflow/support"))
    from gene_family_output_store import cleanup_stale_tmp

    removed = cleanup_stale_tmp(tree.parent, older_than_days=7)
    seconds = time.perf_counter() - start
    digest, files = manifest(tree)
    assert files == (args.families - args.families * 3 // 4) * args.files
    return dict(
        seconds=seconds,
        cleaner_seconds=[v["seconds"] for v in records if v["role"] == "cleaner"],
        writer_seconds=next((v["seconds"] for v in records if v["role"] == "writer"), None),
        max_worker_rss_bytes=max(v["max_rss_bytes"] for v in records),
        retained_sha256=digest,
        retained_files=files,
        removed_directories=len(removed) + sum(v["removed"] for v in records),
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source", required=True)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--families", type=int, default=256)
    p.add_argument("--files", type=int, default=16)
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--scratch")
    p.add_argument("--worker-root")
    p.add_argument("--gate")
    p.add_argument("--writer-probe", action="store_true")
    p.add_argument("--writer-family")
    a = p.parse_args()
    if a.worker_root:
        worker(a)
        return
    runs = []
    for i in range(a.repeats + 1):
        with tempfile.TemporaryDirectory(prefix="gg-tmp-cleanup-", dir=a.scratch) as d:
            record = one(a, Path(d))
        if i:
            runs.append(record)
    assert len({v["retained_sha256"] for v in runs}) == 1
    import statistics

    print(
        json.dumps(
            dict(
                workers=a.workers,
                families=a.families,
                files_per_family=a.files,
                warmups=1,
                runs=runs,
                median_seconds=statistics.median(v["seconds"] for v in runs),
                median_cleaner_seconds=statistics.median(x for v in runs for x in v["cleaner_seconds"]),
                median_writer_seconds=statistics.median(v["writer_seconds"] for v in runs) if a.writer_probe else None,
                median_peak_worker_rss_bytes=statistics.median(v["max_worker_rss_bytes"] for v in runs),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
