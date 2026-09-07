#!/usr/bin/env python3
"""Run inside GeneGalleon: compare archive-family with bounded queue collection.

Use --module-root to point at a baseline checkout's workflow/support directory.
Data stays on the container's filesystem unless --work-dir is supplied.
"""
import argparse
import concurrent.futures
import hashlib
import importlib
import json
import os
from pathlib import Path
import resource
import sys
import tempfile
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--implementation", choices=["per-family", "queued"], required=True)
parser.add_argument("--module-root", type=Path, default=Path(__file__).resolve().parents[1] / "support")
parser.add_argument("--families", type=int, default=10000)
parser.add_argument("--completed", type=int, default=200)
parser.add_argument("--work-dir", type=Path)
parser.add_argument("--queue-workers", type=int, default=1)
args = parser.parse_args()
if not 0 <= args.completed <= args.families or args.queue_workers < 1:
    parser.error("Require 0 <= completed <= families and positive queue-workers")
sys.path.insert(0, str(args.module_root))
m = importlib.import_module("gene_family_output_store")
metrics = {"path_iterdir_entries": 0, "zip_bytes_written": 0, "zip_writes": 0, "compactions": 0}
original_iterdir = Path.iterdir

def counted_iterdir(path):
    for entry in original_iterdir(path):
        metrics["path_iterdir_entries"] += 1
        yield entry

Path.iterdir = counted_iterdir
for function in ("_archive_chunk", "_compact_artifact_chunk"):
    original = getattr(m, function)
    def measured(*a, _original=original, _function=function, **kw):
        result = _original(*a, **kw)
        metrics["zip_bytes_written"] += result[0].stat().st_size
        metrics["zip_writes"] += 1
        metrics["compactions"] += int(_function == "_compact_artifact_chunk")
        return result
    setattr(m, function, measured)

with tempfile.TemporaryDirectory(dir=args.work_dir) as temporary:
    root = Path(temporary) / "orthogroup"
    expected = {}
    for i in range(args.families):
        family = f"OG{i:07d}"
        path = root / "mafft" / f"{family}_cds.aln.fa.gz"
        path.parent.mkdir(parents=True, exist_ok=True)
        data = hashlib.sha256(family.encode()).digest() * 32
        path.write_bytes(data)
        expected[path.name] = hashlib.sha256(data).hexdigest()
    started = time.monotonic()
    if args.implementation == "per-family":
        for i in range(args.completed):
            m.archive_completed_outputs(root, "orthogroup", [f"OG{i:07d}"], m.orthogroup_id_from_name,
                                        include_incomplete=True, preserve_existing_catalog=True)
        queued_seconds = None
    else:
        def enqueue(i):
            family = f"OG{i:07d}"
            path = root / "mafft" / f"{family}_cds.aln.fa.gz"
            inventory = m.family_inventory_path(root, family)
            inventory.mkdir(parents=True, exist_ok=True)
            (inventory / "worker.paths").write_bytes(os.fsencode(path) + b"\0")
            m.enqueue_family_archive(root, "orthogroup", family)
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.queue_workers) as executor:
            list(executor.map(enqueue, range(args.completed)))
        queued_seconds = time.monotonic() - started
        while m.drain_archive_queue(root, "orthogroup", m.orthogroup_id_from_name)["status"] == "committed":
            pass
    elapsed = time.monotonic() - started
    measured_metrics = dict(metrics)
    logical = m.GeneFamilyOutputStore(root)
    for name, digest in expected.items():
        with logical.open_binary("mafft", name) as handle:
            assert hashlib.sha256(handle.read()).hexdigest() == digest
    assert len(list((root / "mafft").glob("*"))) == args.families - args.completed
    print(json.dumps({"implementation": args.implementation, "families": args.families,
                      "completed": args.completed, "queue_workers": args.queue_workers, "elapsed_seconds": elapsed,
                      "enqueue_seconds": queued_seconds,
                      "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                      "output_sha256_verified": len(expected), **measured_metrics}, sort_keys=True))
