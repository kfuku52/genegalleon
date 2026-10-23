#!/usr/bin/env python3
"""Compare promoter N-coordinate extraction in two statistics scripts.

Run in the GeneGalleon runtime with --baseline /saved/orthogroup_statistics.py.
The nested open_text/get_N_coordinate functions are compiled unchanged from each
script's AST; no fake FASTA parser or stand-in kernel is timed. Input creation,
function extraction, pandas import and output hashing are outside timing. Each
case has one warmup and alternating fresh processes; DataFrame schema, dtypes,
row order and every TSV byte must match. Synthetic workloads represent N runs,
gzip input, isolated alternating Ns, N-free promoters and a long-run stress.
This measures the annotation helper, not whole-workflow speed.
"""

import argparse
import ast
import gzip
import hashlib
import json
import platform
import resource
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pandas

PROFILES = ("runs", "runs-gzip", "alternating", "no-n", "long-run")


def extract(path):
    import re

    tree = ast.parse(path.read_text(), filename=str(path))
    main = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "main")
    functions = [node for node in main.body
                 if isinstance(node, ast.FunctionDef) and node.name in {"open_text", "get_N_coordinate"}]
    assert len(functions) == 2
    namespace = {"gzip": gzip, "pandas": pandas, "re": re}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), "exec"), namespace)
    return namespace["get_N_coordinate"]


def memory():
    # VmHWM resets for the fresh executable, unlike inherited ru_maxrss on
    # some launch paths. Keep both measures for an inspectable Linux result.
    values = {}
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith(("VmRSS:", "VmHWM:")):
            key, value, _ = line.split()
            values[key[:-1]] = int(value) * 1024
    assert "VmHWM" in values
    return values


def sha256_path(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def signature(table):
    return {
        "columns": table.columns.tolist(),
        "dtypes": list(map(str, table.dtypes)),
        "rows": len(table),
        "sha256": hashlib.sha256(table.to_csv(sep="\t", index=False).encode()).hexdigest(),
    }


def worker(args):
    function = extract(args.implementation)
    before = memory()
    start_wall, start_cpu = time.perf_counter(), time.process_time()
    table = function(str(args.input))
    wall, cpu = time.perf_counter() - start_wall, time.process_time() - start_cpu
    after = memory()
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    print(json.dumps({"wall_s": wall, "cpu_s": cpu, "peak_rss_bytes": after["VmHWM"],
                      "ru_maxrss_bytes": peak, "memory_before": before, "memory_after": after,
                      "signature": signature(table)}))


def fixtures(root, records, length, long_length):
    root.mkdir(parents=True, exist_ok=True)
    result = {}
    for profile in PROFILES:
        if profile == "runs-gzip":
            continue
        path = root / f"{profile}.fa"
        with path.open("w") as handle:
            count = 1 if profile == "long-run" else records
            size = long_length if profile == "long-run" else length
            for record in range(count):
                if profile == "runs":
                    block = "Nn" * 100 + "ACGT" * 150
                    sequence = (block * (size // len(block) + 1))[:size]
                    sequence = sequence[:-min(100, size)] + "n" * min(100, size)
                elif profile == "alternating":
                    sequence = ("NAnC" * (size // 4 + 1))[:size]
                elif profile == "no-n":
                    sequence = ("ACGT" * (size // 4 + 1))[:size]
                else:
                    sequence = ("Nn" * (size // 2 + 1))[:size]
                handle.write(f">node{record:07d}\n")
                for offset in range(0, len(sequence), 80):
                    handle.write(sequence[offset:offset + 80] + "\n")
        result[profile] = path
    compressed = root / "runs.fa.gz"
    with result["runs"].open("rb") as source, gzip.open(compressed, "wb", compresslevel=3) as target:
        shutil.copyfileobj(source, target, length=1024 * 1024)
    result["runs-gzip"] = compressed
    return result


def run(args):
    implementations = {"baseline": args.baseline.resolve()}
    if not args.baseline_only:
        implementations["candidate"] = args.implementation.resolve()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="gg-promoter-n-benchmark-") as temporary:
        inputs = fixtures(Path(temporary), args.records, args.length, args.long_length)
        rows = []
        for profile in args.profiles:
            for repeat in range(-1, args.repeats):
                order = list(implementations)
                if repeat % 2 == 0:
                    order.reverse()
                comparison = []
                for variant in order:
                    command = [sys.executable, str(Path(__file__).resolve()), "--worker",
                               "--implementation", str(implementations[variant]), "--input", str(inputs[profile])]
                    completed = subprocess.run(command, check=True, text=True, capture_output=True)
                    row = json.loads(completed.stdout)
                    row.update(profile=profile, repeat=repeat, variant=variant, warmup=repeat < 0)
                    rows.append(row)
                    comparison.append(row["signature"])
                assert all(value == comparison[0] for value in comparison), (profile, repeat, comparison)
        summary = {}
        for profile in args.profiles:
            summary[profile] = {}
            for variant in implementations:
                sample = [row for row in rows if row["profile"] == profile and row["variant"] == variant
                          and not row["warmup"]]
                summary[profile][variant] = {key: statistics.median(row[key] for row in sample)
                                             for key in ("wall_s", "cpu_s", "peak_rss_bytes")}
            if len(implementations) == 2:
                before, after = summary[profile]["baseline"], summary[profile]["candidate"]
                summary[profile]["change_percent"] = {
                    key: 100 * (after[key] / before[key] - 1) for key in before
                }
        report = {
            "environment": {"python": sys.version, "platform": platform.platform(),
                            "machine": platform.machine(), "pandas": pandas.__version__},
            "inputs": {profile: {"bytes": inputs[profile].stat().st_size, "sha256": sha256_path(inputs[profile])}
                       for profile in args.profiles},
            "implementations": {variant: {"path": str(path), "sha256": sha256_path(path)}
                                for variant, path in implementations.items()},
            "configuration": {"records": args.records, "length": args.length,
                              "long_length": args.long_length, "repeats": args.repeats},
            "summary": summary, "trials": rows,
        }
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps({"output": str(args.output), "summary": summary}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--baseline-only", action="store_true")
    parser.add_argument("--implementation", type=Path, default=Path("workflow/support/orthogroup_statistics.py"))
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--profiles", nargs="+", choices=PROFILES, default=list(PROFILES))
    parser.add_argument("--records", type=int, default=2000)
    parser.add_argument("--length", type=int, default=2000)
    parser.add_argument("--long-length", type=int, default=2000000)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--output", type=Path, default=Path("tmp/promoter-n-benchmark.json"))
    arguments = parser.parse_args()
    if arguments.worker:
        worker(arguments)
    else:
        if arguments.baseline is None:
            parser.error("--baseline is required")
        if min(arguments.records, arguments.length, arguments.long_length, arguments.repeats) < 1:
            parser.error("workload sizes and repeats must be positive")
        run(arguments)
