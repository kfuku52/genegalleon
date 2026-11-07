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
    parser.add_argument('--kind', choices=['guards', 'digests'], default='guards')
    parser.add_argument('--identity-bytes', type=int, default=0,
                        help='Extra bytes per identity for disk-capacity qualification')
    args = parser.parse_args()
    if not 1 <= args.inventory_members <= 150000 or not 1 <= args.guards <= 100000:
        parser.error('inventory-members must be 1..150000 and guards must be 1..100000')
    if not 0 <= args.identity_bytes <= 1048576:
        parser.error('identity-bytes must be 0..1048576')
    sys.path.insert(0, str(args.support_root.resolve(strict=True)))
    from artifact_audit_guards import DigestSpill, GuardSpill, encode
    inventory = tuple(f'OG{index:07d}.summary_statistics.json' for index in range(args.inventory_members))
    padding = 'X' * args.identity_bytes
    spill = GuardSpill() if args.kind == 'guards' else DigestSpill()
    try:
        spill.guard(('/fixture/output/orthogroup', 'collection_inventory', 'artifact_provenance'),
                    inventory, lambda: inventory)
        before = read_io()
        initial_pages = spill.database.execute('PRAGMA page_count').fetchone()[0]
        began = time.perf_counter()
        for index in range(args.guards):
            key = ('/fixture/output/orthogroup', f'artifact_provenance/OG{index:07d}.summary_statistics.json')
            identity = (str(index), index, None, 9223372036854775807)
            if padding:
                identity = identity, padding
            spill.guard(key, identity, lambda identity=identity: identity)
        # Flush dirty pages for a comparable disk footprint and write count.
        # This private database does not have to survive process termination.
        spill.database.commit()
        elapsed = time.perf_counter() - began
        after = read_io()
        verification_started = time.perf_counter()
        identities = []
        for identity, reader in spill.values():
            assert identity == encode(reader())
            identities.append(hashlib.sha256(identity).digest())
        assert len(identities) == args.guards + 1
        # Reader order is an internal index; compare the complete identity set.
        proof = hashlib.sha256(b''.join(sorted(identities))).hexdigest()
        verification_seconds = time.perf_counter() - verification_started
        pages = spill.database.execute('PRAGMA page_count').fetchone()[0]
        page_size = spill.database.execute('PRAGMA page_size').fetchone()[0]
        print(json.dumps(dict(seconds=elapsed, verification_seconds=verification_seconds,
                              inventory_members=args.inventory_members, kind=args.kind,
                              identity_bytes=args.identity_bytes,
                              database_bytes=pages * page_size, initial_database_bytes=initial_pages * page_size,
                              distinct_guards=len(identities), identity_set_sha256=proof,
                              identities_and_readers_verified=True,
                              io_delta=({key: after[key] - before[key] for key in ['rchar', 'syscr', 'wchar', 'syscw']}
                                        if before is not None else None),
                              peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)))
    finally:
        spill.close()


if __name__ == '__main__':
    main()
