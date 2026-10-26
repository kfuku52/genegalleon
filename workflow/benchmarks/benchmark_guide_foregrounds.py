#!/usr/bin/env python3
"""Differential fresh-process BUSCO projection and binary trait benchmarks.

Run via run_in_runtime.sh without concurrent benchmarks. Timers exclude fixture
creation/imports and complete output comparison. Projection executes the actual
source assignment extracted from build's AST, including the membership set.
"""
import argparse
import ast
import copy
import hashlib
import importlib.util
import json
import os
import pickle
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

SUPPORT = Path(__file__).resolve().parents[1] / "support"


def load(source):
    sys.path.insert(0, str(SUPPORT))
    spec = importlib.util.spec_from_file_location("measured", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def marker_projection(source):
    tree = ast.parse(Path(source).read_text())
    build = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "build")
    projection = next(node for node in ast.walk(build)
                      if isinstance(node, ast.Assign) and isinstance(node.value, ast.DictComp)
                      and any(isinstance(child, ast.Name) and child.id in {"markers", "marker_set"}
                              for child in ast.walk(node.value)))
    prerequisites = [node for node in ast.walk(build) if isinstance(node, ast.Assign)
                     and any(isinstance(target, ast.Name) and target.id == "marker_set"
                             for target in node.targets)]
    function = ast.FunctionDef(
        name="project", args=ast.arguments(posonlyargs=[], args=[ast.arg(arg="data"), ast.arg(arg="markers")],
                                          kwonlyargs=[], kw_defaults=[], defaults=[]),
        body=[*copy.deepcopy(prerequisites), copy.deepcopy(projection),
              ast.Return(value=ast.Name(id="data", ctx=ast.Load()))], decorator_list=[])
    namespace = {}
    exec(compile(ast.fix_missing_locations(ast.Module(body=[function], type_ignores=[])), str(source), "exec"), namespace)
    return namespace["project"]


def worker(args):
    with Path(args.fixture).open("rb") as handle:
        fixture = pickle.load(handle)
    if args.operation == "rows":
        module = load(args.source)
    elif args.operation == "markers":
        implementation = marker_projection(args.source)
    else:
        module = load(args.source)
        header, rows, tree = fixture
    result = error = None
    start_cpu, start_wall = time.process_time(), time.perf_counter()
    try:
        if args.operation == "rows":
            result = module.single_copy_rows(fixture)
        elif args.operation == "markers":
            result = implementation(*fixture)
        else:
            result = (module.resolve_binary_foregrounds(header, rows, tree), rows)
    except Exception as exc:
        error = (type(exc).__name__, str(exc))
    wall, cpu = time.perf_counter() - start_wall, time.process_time() - start_cpu
    rss = next(int(line.split()[1]) for line in Path("/proc/self/status").read_text().splitlines()
               if line.startswith("VmHWM:"))
    if args.operation == "rows" and result is not None:
        result = (list(result[0].items()), sorted(result[1]))
    if args.operation == "markers" and result is not None:
        result = [(species, list(records.items())) for species, records in result.items()]
    with Path(args.payload).open("wb") as handle:
        pickle.dump((result, error), handle, protocol=5)
    print(json.dumps({"wall_seconds": wall, "cpu_seconds": cpu, "rss_kib": rss}))


def fixtures(root):
    from Bio.Phylo.BaseTree import Clade, Tree
    profiles = []
    specs = (("rows-small", 32, 1, 0), ("rows-singletons", 10000, 1, 0),
             ("rows-duplicated", 10000, 8, 0), ("rows-wide", 10000, 8, 256))
    for name, markers, copies, width in specs:
        path = root / (name + ".tsv")
        with path.open("w") as handle:
            handle.write("# synthetic BUSCO full table\n")
            for copy_index in range(copies):
                for marker in range(markers):
                    if copy_index and marker % 4:
                        continue
                    status = "Missing" if marker % 4 == 1 else "Complete"
                    handle.write(f"B{marker}\t{status}\tg{marker}_{copy_index}\t100\t200\t" + "x" * width + "\n")
        profiles.append((name, "rows", path))
    for name, species, total, selected in (("markers-small", 3, 20, 10),
                                           ("markers-default", 64, 5000, 200),
                                           ("markers-large", 128, 10000, 1000)):
        data = {f"S{i}": {f"B{j}": f"M{j}" for j in range(total) if (j + i) % 9}
                for i in range(species)}
        markers = [f"B{j}" for j in range(0, total, max(1, total // selected))][:selected]
        profiles.append((name, "markers", (data, markers)))
    specs = (("traits-small", 8, 1, "balanced"), ("traits-many", 1024, 64, "balanced"),
             ("traits-ladder", 128, 32, "ladder"), ("traits-star", 1024, 32, "star"),
             ("traits-single", 1024, 1, "balanced"), ("traits-nonbinary", 128, 32, "balanced"))
    for name, tips, traits, topology in specs:
        leaves = [Clade(name=f"Plant_species{i:04d}") for i in range(tips)]
        level = leaves[:]
        if topology == "star":
            clade = Clade(clades=level)
        elif topology == "ladder":
            clade = level[0]
            for leaf in level[1:]:
                clade = Clade(clades=[clade, leaf])
        else:
            while len(level) > 1:
                level = [Clade(clades=level[i:i + 2]) for i in range(0, len(level), 2)]
            clade = level[0]
        header = ["species", *(f"trait{j}" for j in range(traits))]
        rows = [[leaf.name, *("manual" if name == "traits-nonbinary" else "NA" if (i + j) % 23 == 0
                              else str(int((i // (2 ** (j % 8))) % 3 == 0))
                              for j in range(traits))] for i, leaf in enumerate(leaves)]
        profiles.append((name, "traits", (header, rows, Tree(root=clade))))
    edges = {"rows-error-late": "A\tComplete\t\nB\tbad\tid\n",
             "rows-error-id": "bad/id\tMissing\nA\tComplete\t\n",
             "rows-error-sole": "A\tComplete\t\n",
             "rows-duplicate-invalid-sole": "A\tComplete\nA\tComplete\tid\n",
             "rows-empty": "# comment\n\n"}
    for name, content in edges.items():
        path = root / (name + ".tsv")
        path.write_text(content)
        profiles.append((name, "rows", path))
    return profiles


def main(args):
    root = Path(args.work_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    sources = {"baseline": Path(args.baseline).resolve()}
    if not args.baseline_only:
        sources["candidate"] = Path(args.candidate).resolve()
    report = {"platform": platform.platform(), "python": platform.python_version(),
              "sources": {key: {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                          for key, path in sources.items()}, "repeats": args.repeats, "profiles": {}}
    for name, operation, value in fixtures(root):
        if args.profiles and name not in args.profiles:
            continue
        fixture = root / (name + ".pkl")
        with fixture.open("wb") as handle:
            pickle.dump(value, handle, protocol=5)
        runs = {key: [] for key in sources}
        edge = any(part in name for part in ("error", "invalid", "empty"))
        for trial in range(1 if edge else args.repeats + 1):
            outputs = {}
            for key in list(sources)[::1 if trial % 2 == 0 else -1]:
                payload = root / (key + ".payload")
                command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--source", str(sources[key]),
                           "--fixture", str(fixture), "--payload", str(payload), "--operation", operation]
                completed = subprocess.run(command, check=True, capture_output=True, text=True,
                                           env={**os.environ, "PYTHONHASHSEED": "0", "OPENBLAS_NUM_THREADS": "1",
                                                "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
                measurement = json.loads(completed.stdout)
                with payload.open("rb") as handle:
                    outputs[key] = pickle.load(handle)
                if trial or edge:
                    runs[key].append(measurement)
            if len(outputs) == 2:
                assert outputs["baseline"] == outputs["candidate"], (name, outputs)
        summary = {key: {metric: statistics.median(run[metric] for run in values)
                         for metric in ("wall_seconds", "cpu_seconds", "rss_kib")}
                   for key, values in runs.items() if values}
        report["profiles"][name] = {"runs": runs, "median": summary, "complete_output_equal": len(sources) == 2}
        print(name, json.dumps(summary), flush=True)
        Path(args.output).write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline")
    parser.add_argument("--candidate")
    parser.add_argument("--baseline-only", action="store_true")
    parser.add_argument("--profiles", nargs="*")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--work-dir")
    parser.add_argument("--output")
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--source")
    parser.add_argument("--fixture")
    parser.add_argument("--payload")
    parser.add_argument("--operation", choices=("rows", "markers", "traits"))
    args = parser.parse_args()
    if args.worker:
        worker(args)
    else:
        main(args)
