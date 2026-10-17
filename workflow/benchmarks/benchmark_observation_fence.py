#!/usr/bin/env python3
"""Measure the complete final metadata fence on a private, stable fixture."""
from __future__ import annotations

import argparse
import hashlib
import json
import resource
import sys
import threading
import time
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', required=True, type=Path)
    parser.add_argument('--create', action='store_true')
    parser.add_argument('--paths', type=int, default=4000)
    parser.add_argument('--support-root', type=Path)
    parser.add_argument('--metadata-delay-ms', type=float, default=0)
    args = parser.parse_args()
    if not 1 <= args.paths <= 250000 or not 0 <= args.metadata_delay_ms <= 100:
        parser.error('paths must be 1..250000 and metadata-delay-ms must be 0..100')
    root = args.workspace
    if args.create:
        root.mkdir()
        for index in range(args.paths):
            (root / f'source-{index:07d}.tsv').write_bytes(b'unchanged\n')
        (root / '.gg-observation-benchmark.json').write_text(json.dumps({'paths': args.paths}))
        return
    if args.support_root is None:
        parser.error('--support-root is required for a measurement')
    marker = json.loads((root / '.gg-observation-benchmark.json').read_text())
    paths = [root / f'source-{index:07d}.tsv' for index in range(marker['paths'])]
    sys.path.insert(0, str(args.support_root.resolve(strict=True)))
    import workflow_observation as observer
    original = observer.path_signature
    final = False
    checked = {}
    lock = threading.Lock()
    def read(path):
        if final and args.metadata_delay_ms:
            time.sleep(args.metadata_delay_ms / 1000)
        signature = original(path)
        if final:
            with lock:
                assert path not in checked
                checked[path] = signature
        return signature
    observer.path_signature = read
    with observer.observe_files():
        for path in paths:
            observer.observe_path(path)
        expected = dict(observer._OBSERVED_PATHS.get())
        final = True
        began = time.perf_counter()
    elapsed = time.perf_counter() - began
    assert checked == expected and len(checked) == marker['paths']
    body = [(str(path), expected[path]) for path in paths]
    print(json.dumps(dict(seconds=elapsed, checked_paths=len(checked),
                          complete_signature_set_sha256=hashlib.sha256(json.dumps(body).encode()).hexdigest(),
                          every_fresh_signature_verified=True, metadata_delay_ms=args.metadata_delay_ms,
                          peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)))


if __name__ == '__main__':
    main()
