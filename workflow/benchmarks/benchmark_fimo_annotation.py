#!/usr/bin/env python3
"""Compare FIMO annotation loading and aggregation in orthogroup_statistics.main.

Fixtures, module imports and verification are outside timers. Each worker
executes the complete actual FIMO branch, including its strict/fallback parser,
not an independently rewritten reference. Complete typed frames and TSV bytes
must match. No scientific tools are stubbed. Run in the GeneGalleon runtime:

python workflow/benchmarks/benchmark_fimo_annotation.py \\
    --baseline /mounted/baseline/orthogroup_statistics.py
"""
import argparse
import ast
import contextlib
import hashlib
import importlib.util
import io
import json
import os
import pickle
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow/support"


def select_statement(source):
    main = next(node for node in ast.parse(source).body if isinstance(node, ast.FunctionDef) and node.name == "main")
    selected = [node for node in main.body if isinstance(node, ast.If) and ast.unparse(node.test) == "os.path.exists(params['fimo'])"]
    if len(selected) != 1:
        raise RuntimeError(f"Expected exactly one main FIMO branch; found {len(selected)}")
    return compile(ast.Module(body=selected, type_ignores=[]), "<actual-main-fimo-statement>", "exec")


def worker(args):
    sys.path.insert(0, str(SUPPORT))
    spec = importlib.util.spec_from_file_location("fimo_measured_implementation", args.source)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    statement = select_statement(Path(args.source).read_text())
    outputs = []
    env = dict(mod.__dict__, params={"fimo": args.fixture}, node_left_merge_tables=outputs)
    stdout = io.StringIO()
    error = None
    wall = time.perf_counter()
    cpu = time.process_time()
    try:
        with contextlib.redirect_stdout(stdout):
            exec(statement, env)
    except Exception as exc:
        error = [type(exc).__name__, str(exc)]
    cpu = time.process_time() - cpu
    wall = time.perf_counter() - wall
    hwm = next(int(line.split()[1]) for line in Path("/proc/self/status").read_text().splitlines() if line.startswith("VmHWM:"))
    # Serialize full frames for strict typed/order equality in the controller.
    with open(args.frame, "wb") as handle:
        pickle.dump(outputs, handle)
    digest = hashlib.sha256()
    for frame in outputs:
        digest.update(json.dumps([list(frame.columns), [str(dtype) for dtype in frame.dtypes], str(frame.index.dtype)]).encode())
        digest.update(frame.to_csv(sep="\t", index=True).encode())
    digest.update(json.dumps([stdout.getvalue(), error], ensure_ascii=False).encode())
    print(json.dumps({"wall_seconds": wall, "cpu_seconds": cpu, "peak_rss_kib": hwm,
                      "digest": digest.hexdigest(), "stdout": stdout.getvalue(), "error": error,
                      "rows": sum(len(frame) for frame in outputs)}))


def fixture(path, groups, hits):
    # Interleave every group so a grouping optimization must retain actual
    # first-occurrence order and nonadjacent annotation order.
    with path.open("w") as handle:
        handle.write("motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n")
        for hit in range(hits):
            for group in range(groups):
                name = ["NA", "NULL", "nan", "001"][group] if group < 4 else f"Species_example_gene{group}"
                alt = "." if hit % 7 == 0 else f"motif_alt{hit % 31}"
                handle.write(f"M{hit % 31}\t{alt}\t{name}\t{hit * 4 + 1}\t{hit * 4 + 6}\t{'+' if hit % 2 else '-'}\t9.3\t0.001\t{hit % 3 * 0.001}\tACGTAC\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, help="Saved baseline orthogroup_statistics.py.")
    parser.add_argument("--candidate", default=str(SUPPORT / "orthogroup_statistics.py"))
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--profiles", nargs="+", default=["small", "family", "singleton", "medium", "large", "dense", "legacy", "invalid", "empty"])
    parser.add_argument("--output", default="tmp/fimo-annotation-benchmark.json")
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--source")
    parser.add_argument("--fixture")
    parser.add_argument("--frame")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    if args.worker:
        worker(args)
        return
    args.baseline = str(Path(args.baseline).resolve(strict=True))
    args.candidate = str(Path(args.candidate).resolve(strict=True))
    import pandas
    work = Path(args.output).with_suffix("")
    work.mkdir(parents=True, exist_ok=True)
    profiles = {"small": (8, 8), "family": (512, 16), "singleton": (512, 1),
                "medium": (2000, 16), "large": (8000, 16), "dense": (100, 1000)}
    all_results = {}
    for profile in args.profiles:
        input_path = work / f"{profile}.tsv"
        if profile in profiles:
            fixture(input_path, *profiles[profile])
        elif profile == "legacy":
            input_path.write_text("#pattern name\tsequence name\tstart\tend\tstrand\tp-value\nM1\tNULL\t1\t8\t+\t0.002\nM2\tNA\t2\t9\t-\t0.005\nM3\tNULL\t8\t15\t.\t0.0001\n")
        elif profile == "invalid":
            input_path.write_text("motif_id\tsequence_name\tstart\tstop\tq-value\nM1\t001\twrong\t1\t.\nM2\tNA\t2\t8\t0.03\nM3\tNA\t5\t9\t0.02\n")
        elif profile == "empty":
            input_path.write_text("motif_id\tsequence_name\tstart\tstop\tq-value\n")
        else:
            raise ValueError(profile)
        measured = []
        for trial in range(-1, args.repeats):
            pair = {}
            frames = {}
            for label in (["baseline", "candidate"] if trial % 2 == 0 else ["candidate", "baseline"]):
                frame_path = work / f"{profile}-{label}.pickle"
                cmd = [sys.executable, str(Path(__file__).resolve()), "--worker", "--baseline", args.baseline,
                       "--source", getattr(args, label),
                       "--fixture", str(input_path), "--frame", str(frame_path)]
                result = subprocess.run(cmd, check=True, text=True, capture_output=True,
                                        env=dict(os.environ, PYTHONHASHSEED="0", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1"))
                pair[label] = json.loads(result.stdout)
                expected_rows = {"small": 8, "family": 512, "singleton": 512,
                                 "medium": 2000, "large": 8000, "dense": 100,
                                 "legacy": 2, "invalid": 1, "empty": 0}[profile]
                assert pair[label]["error"] is None, (profile, label, pair[label])
                assert pair[label]["rows"] == expected_rows, (profile, label, pair[label])
                with frame_path.open("rb") as handle:
                    frames[label] = pickle.load(handle)
            if pair["baseline"]["digest"] != pair["candidate"]["digest"]:
                raise AssertionError(f"Complete output/error mismatch in {profile}: {pair}")
            if len(frames["baseline"]) != len(frames["candidate"]):
                raise AssertionError(f"Frame count mismatch in {profile}")
            for before, after in zip(frames["baseline"], frames["candidate"], strict=True):
                pandas.testing.assert_frame_equal(before, after, check_exact=True)
            if trial >= 0:
                measured.append(pair)
        medians = {label: {key: statistics.median(pair[label][key] for pair in measured)
                           for key in ["wall_seconds", "cpu_seconds", "peak_rss_kib"]}
                   for label in ["baseline", "candidate"]}
        all_results[profile] = {"samples": measured, "medians": medians,
                                "paired_speedups": [p["baseline"]["wall_seconds"] / p["candidate"]["wall_seconds"] for p in measured]}
        print(json.dumps({profile: {"medians": medians, "paired_speedups": all_results[profile]["paired_speedups"]}}), flush=True)
    report = {"environment": {"platform": platform.platform(), "python": platform.python_version(), "pandas": pandas.__version__},
              "source_hashes": {label: hashlib.sha256(Path(getattr(args, label)).read_bytes()).hexdigest() for label in ["baseline", "candidate"]},
              "profiles": all_results}
    Path(args.output).write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
