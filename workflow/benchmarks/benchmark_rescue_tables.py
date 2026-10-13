#!/usr/bin/env python3
"""Compare the exact rescue TSV writer statements on a frozen model sample.

Run in the GeneGalleon runtime. Repetition is a serialization stress workload,
not additional biological candidates. Input loading, parsing and process startup
are outside writer wall time. Whole producer memory is not measured.
"""
import argparse
import ast
import csv
import hashlib
import json
import os
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def fence(path):
    st = Path(path).stat()
    return (st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns)

def writer(source):
    source = Path(source)
    paths = {name: source / "workflow/support" / name
             for name in ("rescue_gene_models.py", "pairwise_synteny.py")}
    identities = {name: fence(path) for name, path in paths.items()}
    content = {name: path.read_bytes() for name, path in paths.items()}
    proof = {name: {"path": str(path), "identity": identities[name],
                    "sha256": hashlib.sha256(content[name]).hexdigest()}
             for name, path in paths.items()}
    verify_sources(proof)
    tree = ast.parse(content["rescue_gene_models.py"].decode())
    outer = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "rescue")
    build = next(n for n in outer.body if isinstance(n, ast.FunctionDef) and n.name == "build")
    def output_statement(node, filename):
        return (isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
                and isinstance(node.value.func, ast.Name) and node.value.func.id == "write_tsv"
                and isinstance(node.value.args[0], ast.BinOp)
                and isinstance(node.value.args[0].right, ast.Constant)
                and node.value.args[0].right.value == filename)
    lo = next(i for i, n in enumerate(build.body) if output_statement(n, "quality_flags.tsv"))
    hi = next(i for i, n in enumerate(build.body) if output_statement(n, "audit.tsv"))
    body = build.body[lo:hi + 1]
    function = ast.parse("def tables(tmp, models, regions, validated):\n    pass\n").body[0]
    function.body = body
    pairwise = ast.parse(content["pairwise_synteny.py"].decode())
    tsv = next(n for n in pairwise.body if isinstance(n, ast.FunctionDef) and n.name == "write_tsv")
    module = ast.fix_missing_locations(ast.Module(body=[tsv, function], type_ignores=[]))
    scope = {"Path": Path, "csv": csv}
    exec(compile(module, str(source), "exec"), scope)
    verify_sources(proof)
    return scope["tables"], ast.dump(module, include_attributes=False), proof

def verify_sources(proof):
    for item in proof.values():
        path = Path(item["path"])
        if (fence(path) != tuple(item["identity"]) or sha(path) != item["sha256"]
                or fence(path) != tuple(item["identity"])):
            raise ValueError("Benchmark writer source changed: " + str(path))

def worker(args):
    sample = Path(args.models)
    initial = fence(sample)
    raw = sample.read_bytes()
    if hashlib.sha256(raw).hexdigest() != args.sample_sha256:
        raise ValueError("Frozen sample SHA differs")
    if fence(sample) != initial:
        raise ValueError("Frozen sample changed")
    records = json.loads(raw)
    if not isinstance(records, list) or not records:
        raise ValueError("Nonempty model sample required")
    models = records * args.repeat
    detected = {m["query"] for m in records}
    regions = [{"id": name, "query": name} for name in sorted(detected)]
    regions += [{"id": 'unmapped\t"query"\n', "query": "no_alignment_fixture"}]
    tables, ast_text, source_proof = writer(args.source)
    before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    tables(output, models, regions, models)
    elapsed = time.perf_counter() - start
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    verify_sources(source_proof)
    results = {name: {"bytes": (output / name).stat().st_size, "sha256": sha(output / name)}
               for name in ("quality_flags.tsv", "audit.tsv")}
    for name, expected in (("quality_flags.tsv", len(models)), ("audit.tsv", len(models) + 1)):
        with (output / name).open(newline="") as handle:
            rows = sum(1 for _ in csv.reader(handle, delimiter="\t"))
        if rows != expected + 1:
            raise ValueError("Wrong TSV row count: " + name)
    if fence(sample) != initial or sha(sample) != args.sample_sha256:
        raise ValueError("Frozen sample changed during benchmark")
    result = {"mode": args.label, "wall_seconds": elapsed, "rss_before_writer_kib": before,
              "peak_rss_kib": peak, "peak_above_before_kib": max(0, peak - before),
              "sample_records": len(records), "model_rows": len(models), "repeat": args.repeat,
              "outputs": results, "writer_ast_sha256": hashlib.sha256(ast_text.encode()).hexdigest(),
              "implementation_sha256": source_proof["rescue_gene_models.py"]["sha256"],
              "writer_sources": {name: item["sha256"] for name, item in source_proof.items()}}
    (output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--models", type=Path, required=True)
    parser.add_argument("--sample-sha256", required=True)
    parser.add_argument("--repeat", type=int, default=128)
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--label", default="legacy")
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.repeat < 1 or args.trials < 1:
        raise ValueError("Positive repeat/trial counts required")
    if args.worker:
        worker(args)
        return
    args.output.mkdir(parents=True, exist_ok=False)
    modes = [(args.label, args.source)]
    if args.candidate is not None:
        modes.append(("streamed", args.candidate))
    records = []
    for trial in range(args.trials + 1):
        order = modes if trial % 2 == 0 else list(reversed(modes))
        for label, source in order:
            out = args.output / f"{label}_{trial}"
            command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--source", str(source),
                       "--models", str(args.models), "--sample-sha256", args.sample_sha256,
                       "--repeat", str(args.repeat), "--label", label, "--output", str(out)]
            process = subprocess.run(command, capture_output=True, text=True)
            (args.output / f"{label}_{trial}.log").write_text(process.stdout + process.stderr)
            if process.returncode:
                raise RuntimeError(f"Benchmark worker {label}/{trial} exit {process.returncode}")
            record = json.loads((out / "result.json").read_text())
            record["trial"] = trial
            record["warmup"] = trial == 0
            records.append(record)
    for label, _ in modes:
        signatures = {(record["writer_ast_sha256"], json.dumps(record["writer_sources"], sort_keys=True))
                      for record in records if record["mode"] == label}
        if len(signatures) != 1:
            raise ValueError("Benchmark writer changed between trials: " + label)
    outputs = [record["outputs"] for record in records]
    if any(value != outputs[0] for value in outputs):
        raise ValueError("Rescue TSV output bytes differ")
    summary = {"status": "passed", "sample": str(args.models), "sample_sha256": args.sample_sha256,
               "warmups": 1, "measured_trials": args.trials, "repeat": args.repeat,
               "records": records, "byte_identical_outputs": outputs[0], "medians": {},
               "command": sys.argv, "python": sys.version, "affinity": sorted(os.sched_getaffinity(0)),
               "scope": "Exact production TSV writer statements; real frozen sample repeated for serialization stress; not whole producer, workflow, cold-cache or 500-species speed"}
    for label, _ in modes:
        measured = [r for r in records if r["mode"] == label and not r["warmup"]]
        summary["medians"][label] = {key: statistics.median(r[key] for r in measured)
            for key in ("wall_seconds", "rss_before_writer_kib", "peak_rss_kib", "peak_above_before_kib")}
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"status": summary["status"], "medians": summary["medians"],
                      "model_rows": records[0]["model_rows"], "byte_identical": True}))

if __name__ == "__main__":
    main()
