#!/usr/bin/env python3
"""Paired, fresh-worker comparison of real expression minima and full CLI.

Run through workflow/tests/run_in_runtime.sh. Setup, input generation, source
loading (kernel), and output fingerprinting are excluded from kernel timings.
CLI timings include source imports, parser/main, table I/O and replacements.
Workers report Linux VmHWM immediately after the measured operation.
"""

import argparse
import ast
import contextlib
import hashlib
import importlib.util
import io
import json
import os
import platform
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "tmp/performance-20261023-round4/expression"
CANDIDATE = ROOT / "workflow/support/orthogroup_statistics.py"
PROFILES = {
    "clade16": ("balanced", 16, 16, 0, 200),
    "small": ("balanced", 32, 32, 0, 100),
    "clade64": ("balanced", 64, 64, 0, 100),
    "clade128": ("balanced", 128, 128, 0, 50),
    "pectinate512": ("pectinate", 512, 512, 0, 1),
    "pectinate1024": ("pectinate", 1024, 1024, 0, 1),
    "constant1024": ("pectinate", 1024, 1024, 0, 1),
    "deep804": ("deep199", 804, 804, 0, 1),
    "balanced1024": ("balanced", 1024, 1024, 0, 1),
    "star1024": ("star", 1024, 1024, 0, 1),
    "missing1024": ("pectinate", 1024, 1024, 4, 10),
    "extra_catalog": ("balanced", 32, 2048, 0, 20),
}


def load_helper(path):
    import numpy

    parsed = ast.parse(path.read_text(), filename=str(path))
    main = next(node for node in parsed.body if isinstance(node, ast.FunctionDef) and node.name == "main")
    functions = [node for node in parsed.body
                 if isinstance(node, ast.FunctionDef) and node.name == "node_is_leaf"]
    functions.extend(node for node in main.body if isinstance(node, ast.FunctionDef) and node.name in {
        "get_node_prop", "get_node_label", "compute_min_expression_corr_by_node",
    })
    namespace = {"numpy": numpy}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), "exec"), namespace)
    return namespace["compute_min_expression_corr_by_node"]


def make_tree(shape, tips):
    import ete4

    leaves = []
    for index in range(tips):
        leaf = ete4.PhyloTree()
        leaf.name = f"Species_{index % 17}_gene{index:05}"
        leaf.dist = 1.0
        leaves.append(leaf)
    if shape == "star":
        tree = ete4.PhyloTree()
        for leaf in leaves:
            tree.add_child(leaf)
    elif shape == "deep199":
        tree = ete4.PhyloTree()
        spine = tree
        for leaf in leaves[:199]:
            spine.add_child(leaf)
            spine = spine.add_child()
        for leaf in leaves[199:]:
            spine.add_child(leaf)
    elif shape == "pectinate":
        tree = leaves[-1]
        for leaf in reversed(leaves[:-1]):
            parent = ete4.PhyloTree()
            parent.add_child(leaf)
            parent.add_child(tree)
            tree = parent
    else:
        level = leaves
        while len(level) > 1:
            following = []
            for start in range(0, len(level), 2):
                if start + 1 == len(level):
                    following.append(level[start])
                    continue
                parent = ete4.PhyloTree()
                parent.add_child(level[start])
                parent.add_child(level[start + 1])
                following.append(parent)
            level = following
        tree = level[0]
    for index, node in enumerate(tree.traverse()):
        node.add_prop("branch_id", index)
        if not node.is_leaf:
            node.name = f"node{index}"
            node.dist = 0.5
    return tree, leaves


def fixture(profile, need_correlation=True):
    import numpy as np

    shape, tips, catalog, missing, repetitions = PROFILES[profile]
    tree, leaves = make_tree(shape, tips)
    rng = np.random.default_rng(20261023)
    expression = rng.normal(size=(catalog, 32))
    if profile == "constant1024":
        # The deepest leaf has a constant profile, so every ancestor contains
        # real corrcoef NaNs. This bounds the reduction fallback overhead.
        expression[-1, :] = 1.0
    corr = np.corrcoef(expression) if need_correlation else None
    names = [leaf.name if not missing or index % missing else f"unused_missing{index}"
             for index, leaf in enumerate(leaves)]
    names.extend(f"unused_extra{index}" for index in range(tips, catalog))
    locations = {name: index for index, name in enumerate(names)}
    digest = hashlib.sha256(expression.tobytes() + repr((shape, list(locations.items()))).encode()).hexdigest()
    return tree, names, locations, corr, expression, repetitions, digest


def fingerprint(result):
    import numpy as np

    payload = [(type(key).__name__, repr(key), type(value).__name__,
                str(np.asarray(value).dtype), np.asarray(value).tobytes().hex())
               for key, value in result.items()]
    return hashlib.sha256(json.dumps(payload).encode()).hexdigest()


def peak_mib():
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("VmHWM:"):
            return int(line.split()[1]) / 1024
    raise RuntimeError("Linux VmHWM is unavailable")


def worker(args):
    source = Path(args.source).resolve()
    tree, names, locations, corr, expression, repetitions, input_digest = fixture(
        args.profile, need_correlation=args.scope == "kernel")
    if args.scope == "kernel":
        helper = load_helper(source)
        wall_start, cpu_start = time.perf_counter(), time.process_time()
        for _ in range(repetitions):
            result = helper(tree, locations, corr)
        wall = (time.perf_counter() - wall_start) / repetitions
        cpu = (time.process_time() - cpu_start) / repetitions
        rss = peak_mib()
        output_digest = fingerprint(result)
        output_count = len(result)
    else:
        # The real CLI is called in a fresh import with a fixed argv[0]. This
        # avoids path-dependent command banners while preserving all I/O.
        with tempfile.TemporaryDirectory(prefix="expression-cli-") as directory:
            directory = Path(directory)
            tree_path = directory / "rooted.nwk"
            tree_path.write_text(tree.write(parser=1, props=[], format_root_node=True) + "\n")
            expression_path = directory / "expression.tsv"
            with expression_path.open("w") as handle:
                handle.write("gene_id\t" + "\t".join(f"tissue{i}" for i in range(32)) + "\n")
                for index, name in enumerate(names):
                    handle.write(name + "\t" + "\t".join(map(str, expression[index])) + "\n")
            # This inner worker has imported only the stdlib before its timer;
            # numpy/ete4 used to prepare fixtures are not preloaded there.
            completed = subprocess.run([
                sys.executable, str(Path(__file__).resolve()), "--cli-worker", "--source", str(source),
                "--directory", str(directory),
            ], text=True, capture_output=True)
            if completed.returncode:
                print(completed.stderr, end="", file=sys.stderr)
                completed.check_returncode()
            measured = json.loads(completed.stdout)
            wall, cpu, rss = (measured[key] for key in ("wall_s", "cpu_s", "peak_mib"))
            output_digest = measured["output_sha256"]
            output_count = 2
    print(json.dumps({"source": str(source), "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                      "profile": args.profile, "scope": args.scope, "input_sha256": input_digest,
                      "wall_s": wall, "cpu_s": cpu, "peak_mib": rss,
                      "output_sha256": output_digest, "output_count": output_count,
                      "repetitions": repetitions if args.scope == "kernel" else 1}))


def cli_worker(args):
    directory, source = Path(args.directory).resolve(), Path(args.source).resolve()
    sys.path.insert(0, str(ROOT / "workflow/support"))
    sys.argv = ["orthogroup_statistics.py", "--rooted_tree", str(directory / "rooted.nwk"),
                "--expression", str(directory / "expression.tsv")]
    stdout, stderr = io.StringIO(), io.StringIO()
    os.chdir(directory)
    wall_start, cpu_start = time.perf_counter(), time.process_time()
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
        spec = importlib.util.spec_from_file_location("expression_perf_worker", source)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.main()
    wall = time.perf_counter() - wall_start
    cpu = time.process_time() - cpu_start
    rss = peak_mib()
    output = [(name, hashlib.sha256((directory / name).read_bytes()).hexdigest())
              for name in ("orthogroup.branch.tsv", "orthogroup.tree.tsv")]
    # Only canonicalize the per-worker fixture path in process text.
    output.extend((name, text.replace(str(directory), "<fixture>"))
                  for name, text in (("stdout", stdout.getvalue()), ("stderr", stderr.getvalue())))
    print(json.dumps({"wall_s": wall, "cpu_s": cpu, "peak_mib": rss,
                      "output_sha256": hashlib.sha256(json.dumps(output).encode()).hexdigest()}))


def controller(args):
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    worker_environment = os.environ.copy()
    worker_environment.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    profiles = args.profiles.split(",")
    records = []
    for scope in args.scopes.split(","):
        for profile in profiles:
            for pair in range(-1, args.pairs):
                order = ("before", "after") if pair % 2 == 0 else ("after", "before")
                pair_records = {}
                for variant in order:
                    path = args.baseline if variant == "before" else args.implementation
                    command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--source", str(path),
                               "--profile", profile, "--scope", scope]
                    completed = subprocess.run(command, text=True, capture_output=True, env=worker_environment)
                    if completed.returncode:
                        print(completed.stdout, end="", flush=True)
                        print(completed.stderr, end="", file=sys.stderr, flush=True)
                        completed.check_returncode()
                    record = json.loads(completed.stdout)
                    record.update(pair=pair, variant=variant)
                    records.append(record)
                    pair_records[variant] = record
                    print(json.dumps(record), flush=True)
                for key in ("input_sha256", "output_sha256", "output_count"):
                    assert pair_records["before"][key] == pair_records["after"][key], (scope, profile, pair, key)
                Path(args.output).write_text(json.dumps({"records": records}, indent=2) + "\n")
    summaries = []
    for scope in args.scopes.split(","):
        for profile in profiles:
            selected = [r for r in records if r["scope"] == scope and r["profile"] == profile and r["pair"] >= 0]
            summary = {"scope": scope, "profile": profile}
            for variant in ("before", "after"):
                values = [r for r in selected if r["variant"] == variant]
                summary[variant] = {key: statistics.median(r[key] for r in values)
                                    for key in ("wall_s", "cpu_s", "peak_mib")}
            ratios = []
            for pair in range(args.pairs):
                values = {r["variant"]: r for r in selected if r["pair"] == pair}
                ratios.append(values["before"]["wall_s"] / values["after"]["wall_s"])
            summary["paired_speedup_range"] = [min(ratios), max(ratios)]
            summaries.append(summary)
    import numpy

    Path(args.output).write_text(json.dumps({"environment": {
        "python": sys.version, "platform": platform.platform(), "machine": platform.machine(),
        "numpy": numpy.__version__, "worker_math_threads": 1,
    }, "records": records, "summaries": summaries}, indent=2) + "\n")
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--cli-worker", action="store_true")
    parser.add_argument("--directory")
    parser.add_argument("--source")
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--implementation", type=Path, default=CANDIDATE)
    parser.add_argument("--profile", choices=PROFILES)
    parser.add_argument("--scope", choices=("kernel", "cli"))
    parser.add_argument("--profiles", default=",".join(PROFILES))
    parser.add_argument("--scopes", default="kernel")
    parser.add_argument("--pairs", type=int, default=5)
    parser.add_argument("--output", default=str(WORK / "comparison.json"))
    options = parser.parse_args()
    if not options.worker and not options.cli_worker and options.baseline is None:
        parser.error("--baseline is required")
    if options.pairs < 1:
        parser.error("--pairs must be positive")
    if options.cli_worker:
        cli_worker(options)
    elif options.worker:
        worker(options)
    else:
        controller(options)
