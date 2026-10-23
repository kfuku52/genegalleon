#!/usr/bin/env python3
"""Count ZIP writer I/O and optionally model per-operation network latency.

This isolates output writing; it is not a measured NFS throughput claim. Both
writers use the same raw filesystem, workload, compression and reader validation.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import resource
import statistics
import subprocess
import sys
import tempfile
import time
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/support"))
from orthofinder_output_store import CHUNK, _StreamingWriter


class CountedRaw(io.RawIOBase):
    def __init__(self, path, latency):
        super().__init__()
        self.handle = io.FileIO(path, "x+")
        self.latency = latency
        self.counts = dict(write=0, seek=0)

    def writable(self):
        return True

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.handle.tell()

    def readinto(self, target):
        return self.handle.readinto(target)

    def _operation(self, kind):
        self.counts[kind] += 1
        if self.latency:
            time.sleep(self.latency)

    def write(self, data):
        self._operation("write")
        return self.handle.write(data)

    def seek(self, *args):
        self._operation("seek")
        return self.handle.seek(*args)

    def close(self):
        if not self.closed:
            super().close()
            self.handle.close()


def trial(writer, files, latency):
    with tempfile.TemporaryDirectory(prefix="orthofinder-writer-") as temporary:
        path = Path(temporary) / "outputs.zip"
        raw = CountedRaw(path, latency)
        digest = hashlib.sha256()
        started = time.perf_counter()
        with raw, (io.BufferedRandom(raw) if writer == "seekable" else io.BufferedWriter(raw, CHUNK)) as buffered:
            output = buffered if writer == "seekable" else _StreamingWriter(buffered)
            with zipfile.ZipFile(output, "w", allowZip64=True) as archive:
                for index in range(files):
                    name = f"Orthogroup_Sequences/OG{index:07d}.fa"
                    payload = (f">Species_gene{index}\n" + "ACGT--ACGT" * 16 + "\n").encode()
                    info = zipfile.ZipInfo(name, date_time=(2026, 1, 1, 0, 0, 0))
                    info.compress_type = zipfile.ZIP_DEFLATED
                    info._compresslevel = 6
                    with archive.open(info, "w", force_zip64=True) as member:
                        member.write(payload)
            buffered.flush()
        elapsed = time.perf_counter() - started
        with zipfile.ZipFile(path) as archive:
            for info in archive.infolist():
                payload = archive.read(info)  # CRC checked by the standard reader.
                digest.update(info.filename.encode())
                digest.update(hashlib.sha256(payload).digest())
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return dict(writer=writer, files=files, latency_ms=latency * 1000, seconds=elapsed,
                    peak_rss_kib=rss / 1024 if sys.platform == "darwin" else rss,
                    io=raw.counts, zip_bytes=path.stat().st_size, logical_sha256=digest.hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--files", type=int, default=8192)
    parser.add_argument("--latency-ms", type=float, default=0)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--trial", choices=("seekable", "streaming"))
    args = parser.parse_args()
    if args.files < 1 or args.repeats < 1 or args.latency_ms < 0:
        parser.error("positive counts and nonnegative latency required")
    if args.trial:
        print(json.dumps(trial(args.trial, args.files, args.latency_ms / 1000)))
        return
    rows = []
    for repetition in range(args.repeats + 1):
        for writer in ("seekable", "streaming"):
            command = [sys.executable, "-B", __file__, "--files", str(args.files),
                       "--latency-ms", str(args.latency_ms), "--trial", writer]
            row = json.loads(subprocess.run(command, check=True, capture_output=True, text=True).stdout)
            row.update(warmup=repetition == 0, repetition=repetition)
            rows.append(row)
            print(json.dumps(row), flush=True)
    if len({row["logical_sha256"] for row in rows}) != 1:
        raise RuntimeError("Different logical outputs")
    print(json.dumps(dict(summary={writer: {key: statistics.median(row[key] for row in rows
                           if row["writer"] == writer and not row["warmup"])
                           for key in ("seconds", "peak_rss_kib", "zip_bytes")}
                           for writer in ("seekable", "streaming")},
                          output_equivalence=True, python=sys.version, platform=sys.platform)))


if __name__ == "__main__":
    main()
