#!/usr/bin/env python3
"""Compare exact terminal checks in fresh processes; use a fixed container."""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import importlib.util
import io
import json
import os
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from types import SimpleNamespace


class Store:
    def __init__(self, root):
        self.root = root

    def source_artifact(self, subdir, name):
        path = self.root / subdir
        return SimpleNamespace(logical_path=f"{subdir}/{name}", mtime_ns=path.stat().st_mtime_ns)

    def open_binary(self, subdir, name):
        return (self.root / subdir).open("rb")


def load_validator(support):
    helper = support / "terminal_artifact_validation.py"
    if helper.exists():
        spec = importlib.util.spec_from_file_location("terminal_benchmark", helper)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module.verify_terminal_profile, helper
    source = support / "workflow_api.py"
    tree = ast.parse(source.read_text())
    function = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                    and node.name == "verify_terminal_profile")
    namespace = {"csv": csv, "hashlib": hashlib, "io": io}
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(source), "exec"), namespace)
    return namespace["verify_terminal_profile"], source


def worker(args):
    validate, source = load_validator(args.support)
    start = time.perf_counter()
    proof = validate(Store(args.data), "family")
    elapsed = time.perf_counter() - start
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(json.dumps(dict(wall_seconds=elapsed, peak_rss_bytes=peak * (1 if sys.platform == "darwin" else 1024),
                          proof_sha256=hashlib.sha256(json.dumps(proof, sort_keys=True).encode()).hexdigest(),
                          state=proof["state"], source_sha256=hashlib.sha256(source.read_bytes()).hexdigest())))


def compare(args):
    with tempfile.TemporaryDirectory(prefix="gg-terminal-benchmark-") as directory:
        data = Path(directory)
        with (data / "stat_branch").open("w") as handle:
            handle.write("branch_id\tnode_name\tvalue\n")
            for index in range(args.rows):
                handle.write(f"{index}\tnode{index}\t{'x' * 96}\n")
        (data / "stat_tree").write_text("num_leaf\n100000\n")
        (data / "tree_plot").write_bytes(b"%PDF-1.7\n" + b"x" * 1024 * 1024 + b"\n%%EOF\n")
        paths = dict(baseline=args.baseline_support, candidate=args.candidate_support)
        samples = {key: [] for key in paths}
        for trial in range(args.trials + 1):
            order = list(paths) if trial % 2 == 0 else list(reversed(paths))
            for name in order:
                output = subprocess.check_output([sys.executable, "-B", str(Path(__file__).resolve()),
                                                  "--support", str(paths[name]), "--data", str(data)], text=True)
                result = json.loads(output)
                if trial:
                    samples[name].append(result)
        fingerprints = {row["proof_sha256"] for rows in samples.values() for row in rows}
        assert len(fingerprints) == 1 and all(row["state"] == "verified" for rows in samples.values() for row in rows)
        report = dict(python=sys.version, platform=sys.platform, cpu_count=os.cpu_count(), rows=args.rows,
                      input_bytes={p.name: p.stat().st_size for p in data.iterdir()}, warmups=1, trials=args.trials,
                      scope="native terminal structural validation only; warm filesystem cache; imports excluded from wall timer and included in RSS", equivalent_proofs=True,
                      samples=samples, medians={name: {key: statistics.median(row[key] for row in values)
                                                      for key in ("wall_seconds", "peak_rss_bytes")}
                                                for name, values in samples.items()})
        print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support", type=Path)
    parser.add_argument("--data", type=Path)
    parser.add_argument("--baseline-support", type=Path)
    parser.add_argument("--candidate-support", type=Path)
    parser.add_argument("--rows", type=int, default=100000)
    parser.add_argument("--trials", type=int, default=5)
    args = parser.parse_args()
    if args.support:
        worker(args)
    else:
        if not (args.baseline_support and args.candidate_support and 1 <= args.rows <= 1000000 and 1 <= args.trials <= 20):
            parser.error("provide both support paths and bounded positive rows/trials")
        compare(args)
