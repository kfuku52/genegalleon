#!/usr/bin/env python3
"""Measure large inventory identities beside small, unrelated audit guards."""
from __future__ import annotations

import argparse
import hashlib
import json
import resource
import sys
import time
from pathlib import Path


def read_io():
    source = Path('/proc/self/io')
    if not source.exists():
        return None
    return {key: int(value) for key, value in
            (line.split(':', 1) for line in source.read_text().splitlines())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', required=True, type=Path)
    parser.add_argument('--inventory-members', type=int, default=107841)
    parser.add_argument('--guards', type=int, default=1000)
    args = parser.parse_args()
    if not 1 <= args.inventory_members <= 150000 or not 1 <= args.guards <= 100000:
        parser.error('inventory-members must be 1..150000 and guards must be 1..100000')
    sys.path.insert(0, str(args.support_root.resolve(strict=True)))
    from artifact_audit_guards import GuardSpill, encode
    inventory = tuple(f'OG{index:07d}.summary_statistics.json' for index in range(args.inventory_members))
    spill = GuardSpill()
    try:
        spill.guard(('/fixture/output/orthogroup', 'collection_inventory', 'artifact_provenance'),
                    inventory, lambda: inventory)
        before = read_io()
        began = time.perf_counter()
        for index in range(args.guards):
            key = ('/fixture/output/orthogroup', f'artifact_provenance/OG{index:07d}.summary_statistics.json')
            identity = (str(index), index, None, 9223372036854775807)
            spill.guard(key, identity, lambda identity=identity: identity)
        elapsed = time.perf_counter() - began
        after = read_io()
        rows = list(spill.values())
        assert len(rows) == args.guards + 1
        assert all(identity == encode(reader()) for identity, reader in rows)
        # Reader order is an internal index; compare the complete identity set.
        proof = hashlib.sha256(b'\n'.join(sorted(identity for identity, _ in rows))).hexdigest()
        print(json.dumps(dict(seconds=elapsed, inventory_members=args.inventory_members,
                              distinct_guards=len(rows), identity_set_sha256=proof,
                              identities_and_readers_verified=True,
                              io_delta=({key: after[key] - before[key] for key in ['rchar', 'syscr']}
                                        if before is not None else None),
                              peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)))
    finally:
        spill.close()


if __name__ == '__main__':
    main()
