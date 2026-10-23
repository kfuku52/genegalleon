#!/usr/bin/env python3
"""Compare raw/ZIP logical reads and 1/2/4-worker lossless conversion on one filesystem."""
from __future__ import annotations

import argparse
import json
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/support"))
import artifact_provenance as provenance
import orthofinder_output_store as store


def trial(work_root, files, workers):
    with tempfile.TemporaryDirectory(prefix="orthofinder-bench-", dir=work_root) as temporary:
        root = Path(temporary) / "orthofinder"
        root.mkdir()
        (root / "Log.txt").write_text("OrthoFinder run completed\n")
        for name in store.MANAGED_DIRECTORIES:
            directory = root / name
            directory.mkdir(parents=True)
            for index in range(files):
                (directory / f"OG{index:07d}.fa").write_text(f">Species_gene{index}\n" + "ACGT--ACGT" * 16 + "\n")
        start = time.perf_counter()
        before = provenance.sha256_path(root)
        raw_read = time.perf_counter() - start
        start = time.perf_counter()
        store.convert(root, workers=workers)
        conversion = time.perf_counter() - start
        start = time.perf_counter()
        after = provenance.sha256_path(root)
        zip_read = time.perf_counter() - start
        if before != after:
            raise RuntimeError("Logical output differs")
        return dict(workers=workers, files_per_directory=files, raw_read_seconds=raw_read,
                    conversion_seconds=conversion, zip_read_seconds=zip_read,
                    peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 if sys.platform == "darwin" else 1),
                    logical_sha256=before[0], logical_bytes=before[1],
                    raw_files=files * 4 + 1, zip_files=sum(p.is_file() for p in root.rglob("*")))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-root", required=True, type=Path)
    parser.add_argument("--files-per-directory", type=int, default=256)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--trial", type=int, choices=(1, 2, 4))
    args = parser.parse_args()
    if args.files_per_directory < 1 or args.repeats < 1:
        parser.error("positive file and repetition counts are required")
    if args.trial:
        print(json.dumps(trial(args.work_root, args.files_per_directory, args.trial)))
        return
    rows = []
    for workers in (1, 2, 4):
        for repetition in range(args.repeats + 1):
            command = [sys.executable, "-B", __file__, "--work-root", str(args.work_root),
                       "--files-per-directory", str(args.files_per_directory), "--trial", str(workers)]
            completed = subprocess.run(command, check=True, capture_output=True, text=True)
            row = json.loads(completed.stdout)
            row.update(warmup=repetition == 0, repetition=repetition)
            rows.append(row)
            print(json.dumps(row), flush=True)
    if len({row["logical_sha256"] for row in rows}) != 1:
        raise RuntimeError("Different logical outputs across benchmark trials")
    summary = {str(workers): {key: statistics.median(row[key] for row in rows
                                                  if row["workers"] == workers and not row["warmup"])
                              for key in ("raw_read_seconds", "conversion_seconds", "zip_read_seconds", "peak_rss_kib")}
               for workers in (1, 2, 4)}
    print(json.dumps(dict(summary=summary, python=sys.version, platform=sys.platform,
                         filesystem_path=str(args.work_root), output_equivalence=True)), flush=True)


if __name__ == "__main__":
    main()
