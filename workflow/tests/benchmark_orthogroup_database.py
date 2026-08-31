"""Measure DB builds in a GeneGalleon runtime, with deterministic logical-output hashes.

Run --source against a saved before/after implementation and compare the JSON
outputs. The many-files case models 512 families with 256 branches; large-file
is a separate synthetic stress case with the same total rows in one input.
"""

import argparse
import hashlib
import importlib.util
import json
import os
import resource
import sqlite3
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def fingerprint(database):
    digest = hashlib.sha256()
    tables = {}
    with sqlite3.connect(database) as connection:
        schema = connection.execute(
            "SELECT name, sql FROM sqlite_master WHERE type='table' ORDER BY name"
        ).fetchall()
        for table, statement in schema:
            digest.update(statement.encode())
            quoted = '"' + table.replace('"', '""') + '"'
            columns = connection.execute(f"PRAGMA table_info({quoted})").fetchall()
            order = ",".join('"' + column[1].replace('"', '""') + '"' for column in columns)
            cursor = connection.execute(f"SELECT * FROM {quoted} ORDER BY {order}")
            count = 0
            while rows := cursor.fetchmany(1024):
                count += len(rows)
                for row in rows:
                    digest.update(json.dumps(row, ensure_ascii=False, separators=(",", ":")).encode())
                    digest.update(b"\n")
            tables[table] = count
    return digest.hexdigest(), tables


def worker(args):
    sys.path.insert(0, str(REPO_ROOT / "workflow/support"))
    spec = importlib.util.spec_from_file_location("benchmark_database", args.source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    root = args.worker.resolve()
    os.chdir(root)
    sys.argv = [str(args.source), "--overwrite", "1", "--dbpath", str(root / "result.sqlite3"),
                "--dir_stat_tree", str(root / "stat_tree"), "--dir_stat_branch", str(root / "stat_branch"),
                "--row_threshold", "4096", "--ncpu", "2"]
    started = time.monotonic()
    module.main()
    seconds = time.monotonic() - started
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    output_hash, tables = fingerprint(root / "result.sqlite3")
    (root / "measurement.json").write_text(json.dumps({
        "seconds": seconds, "peak_rss_kib_linux": peak_rss,
        "logical_sha256": output_hash, "table_rows": tables,
    }, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=REPO_ROOT / "workflow/support/generate_orthogroup_database.py")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker", type=Path)
    args = parser.parse_args()
    args.source = args.source.resolve()
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error("--output is required")
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    results = {"source": str(args.source), "workers": 2, "row_threshold": 4096, "numeric_metrics": 24, "cases": []}
    with tempfile.TemporaryDirectory(prefix="gg-db-benchmark-") as temporary:
        for case, families, branches in (("many-files", 512, 256), ("large-file", 1, 131072)):
            root = Path(temporary) / case
            (root / "stat_tree").mkdir(parents=True)
            (root / "stat_branch").mkdir()
            for family in range(families):
                name = f"OG{family:07d}"
                (root / "stat_tree" / f"{name}_stat.tree.tsv").write_text(
                    "num_branch\tnum_spe\tnum_dup\tnum_sp\n" + f"{branches}\t128\t0\t127\n"
                )
                with (root / "stat_branch" / f"{name}_stat.branch.tsv").open("w") as handle:
                    handle.write("branch_id\tnode_name\tnum_sp\tso_event\t" + "\t".join(f"metric_{n}" for n in range(24)) + "\n")
                    for branch in range(branches):
                        handle.write(f"{branch}\tnode{branch}\t2\tS\t" + "\t".join(str(family + branch + n) for n in range(24)) + "\n")
            with output.with_name(output.stem + f"-{case}.log").open("w") as log:
                subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker", str(root),
                                "--source", str(args.source)], check=True, stdout=log, stderr=subprocess.STDOUT)
            measurement = json.loads((root / "measurement.json").read_text())
            measurement.update(case=case, families=families, branches_per_family=branches)
            results["cases"].append(measurement)
            output.write_text(json.dumps(results, indent=2) + "\n")
            print(json.dumps(measurement), flush=True)


if __name__ == "__main__":
    main()
