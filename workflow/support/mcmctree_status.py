#!/usr/bin/env python3
"""Bind advisory MCMC status to tree bytes; warn on missing or stale evidence."""
import argparse
import hashlib
import json
from pathlib import Path
import sys


def tree_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_status(status, tree):
    try:
        data = json.loads(status.read_text())
        if data.get('schema_version') != 1 or data.get('tree_sha256') != tree_hash(tree):
            return {'state': 'legacy_unverified', 'reason': 'Missing or mismatched tree fingerprint'}
        return data
    except (OSError, ValueError, TypeError, AttributeError):
        return {'state': 'legacy_unverified', 'reason': 'Missing or unreadable diagnostics'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['bind', 'copy', 'warn'])
    parser.add_argument('--tree', type=Path, required=True)
    parser.add_argument('--status', type=Path, required=True)
    parser.add_argument('--target', type=Path)
    parser.add_argument('--cached', action='store_true', help='Refresh only a verified conversion of this source tree')
    args = parser.parse_args()
    if args.action == 'bind':
        data = json.loads(args.status.read_text())
        data['schema_version'] = 1
        data['tree_sha256'] = tree_hash(args.tree)
        args.status.write_text(json.dumps(data, indent=2) + '\n')
    else:
        data = read_status(args.status, args.tree)
        if args.action == 'copy':
            if args.target is None:
                parser.error('--target required for copy')
            target_status = Path(str(args.target) + '.convergence.json')
            if args.cached:
                prior = read_status(target_status, args.target)
                if prior.get('source_tree_sha256') != tree_hash(args.tree):
                    return
            data['schema_version'] = 1
            data['source_tree_sha256'] = tree_hash(args.tree)
            data['tree_sha256'] = tree_hash(args.target)
            pending = target_status.with_suffix('.pending')
            pending.write_text(json.dumps(data, indent=2) + '\n')
            pending.replace(target_status)
        elif data.get('state') != 'passed':
            print(f"WARNING: Species dating status={data.get('state', 'legacy_unverified')}. "
                  f"Downstream analyses continue with provisional ages. Diagnostics: {args.status}", file=sys.stderr)


if __name__ == '__main__':
    main()
