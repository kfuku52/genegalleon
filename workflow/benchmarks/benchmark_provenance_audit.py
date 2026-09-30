#!/usr/bin/env python3
"""Compare complete audits on a private reproducible raw/ZIP/mixed fixture."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import subprocess
import sys
import time
from pathlib import Path


def digest(data):
    return hashlib.sha256(data).hexdigest()


def create_fixture(root, families, steps, shared_mib, layout="raw", support=None):
    if root.exists():
        raise ValueError("fixture destination must not exist")
    store = root / "output/orthogroup"
    shared = root / "input/shared.tsv"
    shared.parent.mkdir(parents=True)
    data = (b"species\tgene\tannotation\n" * 50000)
    data = (data * (shared_mib * 1048576 // len(data) + 1))[:shared_mib * 1048576]
    shared.write_bytes(data)
    manifests = store / "artifact_provenance"
    tables = store / "stat_branch"
    manifests.mkdir(parents=True)
    tables.mkdir()
    shared_entry = {"scope": "workspace", "path": "input/shared.tsv", "artifact_type": "file",
                    "sha256": digest(data), "size_bytes": len(data), "label": "shared"}
    for index in range(families):
        family = f"OG{index:07d}"
        data = random.Random(index).randbytes(16384)
        name = family + ".tsv"
        (tables / name).write_bytes(data)
        entry = {"scope": "logical", "path": "stat_branch/" + name, "artifact_type": "file",
                 "sha256": digest(data), "size_bytes": len(data), "label": "table"}
        for step_index in range(steps):
            step = "summary_statistics" if step_index == 0 else f"step_{step_index:02d}"
            payload = {"schema_version": 1, "family_id": family, "step": step,
                       "inputs": [shared_entry] if step_index == 0 else [entry], "outputs": [entry]}
            (manifests / f"{family}.{step}.json").write_text(json.dumps(payload))
    if layout != "raw":
        if support is None:
            raise ValueError("ZIP fixture creation requires --support-root")
        catalog = root / "families.txt"
        catalog.write_text("".join(f"OG{index:07d}\n" for index in range(families)))
        subprocess.run([sys.executable, "-B", str(support / "gene_family_output_store.py"),
                        "convert-storage", "--root", str(store), "--mode", "orthogroup",
                        "--to", "zip", "--family-id-file", str(catalog), "--compression", "store",
                        "--max-files-per-shard", "5000", "--progress-interval", "0"], check=True)
        if layout == "mixed":
            tables.mkdir(exist_ok=True)
            for index in range(0, families, 2):
                (tables / f"OG{index:07d}.tsv").write_bytes(random.Random(index).randbytes(16384))
    marker = {"families": families, "steps": steps, "shared_mib": shared_mib, "layout": layout}
    (root / ".gg-audit-benchmark-fixture.json").write_text(json.dumps(marker))


def run(args):
    root = args.workspace.resolve(strict=True)
    marker = json.loads((root / ".gg-audit-benchmark-fixture.json").read_text())
    support = args.support_root.resolve(strict=True)
    report = root / f"benchmark-audit-{args.label}.tsv"
    code = r'''
import json,resource,sys,time
sys.path.insert(0,sys.argv[1])
import artifact_provenance as p
begin=time.perf_counter()
status=p.dispatch(sys.argv[2:])
print(json.dumps({'status':status,'seconds':time.perf_counter()-begin,
                 'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))
'''
    command = [sys.executable, "-B", "-c", code, str(support), "audit", "--logical-root",
               str(root / "output/orthogroup"), "--workspace-root", str(root), "--mode", "orthogroup",
               "--output-tsv", str(report)]
    if args.workers is not None:
        command.extend(["--workers", str(args.workers)])
    if args.no_progress:
        command.extend(["--progress-interval", "0"])
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
    began = time.perf_counter()
    completed = subprocess.run(command, capture_output=True, text=True, env=env, check=False)
    if completed.returncode:
        raise RuntimeError(completed.stderr + completed.stdout)
    result = json.loads(completed.stdout.splitlines()[-1])
    if result["status"] != 0:
        raise RuntimeError(completed.stderr + completed.stdout)
    result.update(label=args.label, wall_seconds=time.perf_counter()-began,
                  report_sha256=digest(report.read_bytes()), fixture=marker,
                  python=sys.version.split()[0], support_root=str(support),
                  command=command[4:], cpu_limit=os.environ.get("GG_BENCHMARK_CPU_LIMIT"),
                  cache_condition="run-scoped RAM cache empty; baseline persistent digest cache warm after warmup; OS cache not forcibly cleared")
    print(json.dumps(result, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True, type=Path)
    parser.add_argument("--create", action="store_true")
    parser.add_argument("--families", type=int, default=1000)
    parser.add_argument("--steps", type=int, default=15)
    parser.add_argument("--shared-mib", type=int, default=24)
    parser.add_argument("--layout", choices=["raw", "zip", "mixed"], default="raw")
    parser.add_argument("--support-root", type=Path)
    parser.add_argument("--label", default="trial")
    parser.add_argument("--workers", type=int)
    parser.add_argument("--no-progress", action="store_true")
    args = parser.parse_args()
    if args.create:
        create_fixture(args.workspace, args.families, args.steps, args.shared_mib, args.layout, args.support_root)
    else:
        if args.support_root is None:
            parser.error("--support-root is required for a measurement")
        run(args)


if __name__ == "__main__":
    main()
