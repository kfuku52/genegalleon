#!/usr/bin/env python3
"""Measure the complete audit content/identity fence with bounded I/O tails.

Run this same runner in fresh processes against baseline/candidate support
roots. The optional delay models heterogeneous filesystem latency; it is not
an HPC speedup measurement. Every source is really read and rehashed.
"""
import argparse
import hashlib
import json
import platform
import resource
import sys
import tempfile
import threading
import time
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, required=True)
    parser.add_argument('--files', type=int, default=1024)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--latency', choices=['local', 'tails'], default='tails')
    args = parser.parse_args()
    if not 1 <= args.files <= 250000 or not 1 <= args.workers <= 64:
        parser.error('files/workers exceed the audit bounds')
    sys.path.insert(0, str(args.support_root.resolve()))
    from artifact_audit_runtime import AuditDigests, signature

    with tempfile.TemporaryDirectory(prefix='gg-audit-revalidation-') as tmp:
        root = Path(tmp)
        memo = AuditDigests()
        content_reads, guard_reads = [0] * args.files, [0] * args.files
        lock = threading.Lock()
        validating = False

        def delay(index):
            if args.latency == 'tails':
                time.sleep(0.02 if index % (args.workers * 8) == args.workers * 8 - 1 else 0.0005)

        for index in range(args.files):
            path = root / str(index)
            path.write_bytes((str(index).encode() + b'\n').ljust(8192, b'A'))

            def compute(path=path, index=index):
                if validating:
                    delay(index)
                raw = path.read_bytes()
                with lock:
                    content_reads[index] += 1
                return hashlib.sha256(raw).hexdigest(), len(raw)

            def guard(path=path, index=index):
                delay(index)
                with lock:
                    guard_reads[index] += 1
                return signature(path)

            memo.read(path, compute)
            memo.guard(str(index), signature(path), guard)

        validating = True
        started = time.perf_counter()
        memo.validate(workers=args.workers)
        elapsed = time.perf_counter() - started
        assert content_reads == [2] * args.files and guard_reads == [1] * args.files
        assert memo.metrics()['digest_revalidation_reads'] == args.files
        values = [(Path(key[0]).name, entry[1]) for key, entry in memo.entries.items()]
        fingerprint = hashlib.sha256(json.dumps(values, separators=(',', ':')).encode()).hexdigest()
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        print(json.dumps(dict(elapsed_seconds=elapsed, peak_rss_MiB=peak / (1024**2 if sys.platform == 'darwin' else 1024),
                              files=args.files, workers=args.workers, latency=args.latency,
                              complete_source_set_sha256=fingerprint, metrics=memo.metrics(),
                              each_content_read_twice=True, every_guard_checked=True,
                              python=platform.python_version(), platform=platform.platform()), sort_keys=True))


if __name__ == '__main__':
    main()
