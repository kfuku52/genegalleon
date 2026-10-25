#!/usr/bin/env python3
"""Compare GFF species task planning plus complete small/multispecies CLI calls.

Planning profiles execute the real main prefix through source assignment, then
return the complete task list. The AST cutoff is the first ``frames = []`` after
planning; GFF parsing and output publication are excluded only for that scope.
Complete profiles execute the unmodified main and compare complete TSV bytes,
parsed frame types, logs, errors and input hashes. Imports/fixture generation and
verification are untimed. Each profile uses alternating fresh workers, a discarded
warmup pair and measured repeats. Linux VmHWM includes imports and retained input.
"""
import argparse
import ast
import contextlib
import gzip
import hashlib
import importlib.util
import io
import json
import os
import pickle
import platform
import re
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve()
while not (ROOT / "workflow/support").is_dir() and ROOT != ROOT.parent:
    ROOT = ROOT.parent
SUPPORT = ROOT / "workflow/support"
PROFILES = {
    "small": ("main", 4, 4, False),
    "multispecies": ("main", 64, 16, False),
    "gzip": ("main", 8, 8, True),
    "planning_128": ("plan", 128, 2000, False),
    "planning_512": ("plan", 512, 1000, False),
    "planning_one": ("plan", 1, 100000, False),
    "edge_qualifiers": ("main", 4, 3, False),
    "edge_duplicate": ("main", 2, 2, False),
    "edge_malformed": ("main", 2, 2, False),
    "edge_no_matches": ("main", 2, 2, False),
    "edge_threads": ("main", 8, 4, False),
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def make_fixture(root, name):
    scope, species_count, genes, compressed = PROFILES[name]
    root.mkdir()
    directory = root / "gff"
    directory.mkdir()
    labels = [f"Genus_species{i:04d}" for i in range(species_count)]
    if name == "edge_qualifiers":
        labels = ["Oldenlandia_corymbosa_var._corymbosa", "Oryza_sativa_Indica_Group",
                  "Plant_cf._alba", "Genus_sp._isolate"]
    fasta = root / ("queries.fa.gz" if compressed else "queries.fa")
    opener = gzip.open if compressed else open
    with opener(fasta, "wt") as handle:
        # Interleaving species must retain per-species FASTA order.
        for index in reversed(range(genes)):
            for label in reversed(labels):
                handle.write(f">{label}_g{index:06d} header α\nATGAAATAA\n")
        if name == "edge_duplicate":
            handle.write(f">{labels[0]}_g000000 repeated\nATG\n")
        if name == "edge_malformed":
            handle.write(">\nATG\n")
    for label in reversed(labels):
        with opener(directory / (label + (".gff.gz" if compressed else ".gff")), "wt") as handle:
            handle.write("##gff-version 3\n")
            # Planning inputs are real GFFs too, but deliberately do not repeat
            # genome annotation construction, which is outside this scope.
            for index in range(genes if scope == "main" else 1):
                gene = f"g{index:06d}"
                start = 1 + index * 30
                handle.write(f"chr1\tfixture\tgene\t{start}\t{start+8}\t.\t+\t.\tID={gene}\n")
                handle.write(f"chr1\tfixture\tmRNA\t{start}\t{start+8}\t.\t+\t.\tID={gene}.t;Parent={gene}\n")
                handle.write(f"chr1\tfixture\tCDS\t{start}\t{start+8}\t.\t+\t0\tID={gene}.c;Parent={gene}.t\n")
    (directory / ".hidden_species.gff").write_text("invalid\n")
    (directory / "unrelated.txt").write_text("invalid\n")
    (directory / "Absent_species.gff").write_text("invalid\n")
    if name == "edge_no_matches":
        fasta.write_text(">Missing_species_x\nATG\n")
    return fasta, directory


def normalized(text):
    text = re.sub(r"\d{4}-\d\d-\d\d \d\d:\d\d:\d\d(?:\.\d+)?", "<datetime>", text)
    return re.sub(r"completed in [\d,]+ secs:", "completed in <elapsed> secs:", text)


def worker(args):
    sys.path.insert(0, str(SUPPORT))
    spec = importlib.util.spec_from_file_location("candidate_gff_stats", args.source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if args.scope == "plan":
        tree = ast.parse(Path(args.source).read_text())
        main = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "main")
        cut = next(i for i, node in enumerate(main.body) if isinstance(node, ast.Assign)
                   and any(isinstance(target, ast.Name) and target.id == "frames" for target in node.targets))
        assert isinstance(main.body[cut].value, ast.List) and main.body[cut].value.elts == []
        main.body = main.body[:cut] + [ast.Return(value=ast.Name(id="tasks", ctx=ast.Load()))]
        prefix = ast.fix_missing_locations(ast.Module(body=[main], type_ignores=[]))
        exec(compile(prefix, args.source, "exec"), module.__dict__)
    output = Path(args.output_dir)
    output.mkdir()
    sys.argv = [args.source, "--mode", "gene_delim", "--dir_gff", args.gff_dir,
                "--seqfile", args.fasta, "--outfile", str(output / "stats.tsv"),
                "--ncpu", "2" if args.profile == "edge_threads" else "1"]
    if args.profile == "edge_no_matches":
        sys.argv.append("--require-matches")
    stdout, stderr = io.StringIO(), io.StringIO()
    error, result = None, None
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
        begin, cpu_begin = time.perf_counter(), time.process_time()
        try:
            result = module.main()
        except Exception as exc:
            error = (type(exc).__name__, str(exc))
        wall, cpu = time.perf_counter() - begin, time.process_time() - cpu_begin
    rss = next(int(line.split()[1]) for line in Path("/proc/self/status").read_text().splitlines()
               if line.startswith("VmHWM:"))
    payload = {"result": result, "error": error, "stdout": normalized(stdout.getvalue()),
               "stderr": normalized(stderr.getvalue())}
    path = output / "stats.tsv"
    if path.exists():
        import pandas
        payload["bytes"] = path.read_bytes()
        payload["frame"] = pandas.read_csv(path, sep="\t")
    with Path(args.payload).open("wb") as handle:
        pickle.dump(payload, handle, protocol=5)
    print(json.dumps({"wall_seconds": wall, "cpu_seconds": cpu, "peak_rss_kib": rss,
                      "source_sha256": sha256(args.source)}))


def main(args):
    import pandas
    report = {"python": platform.python_version(), "platform": platform.platform(),
              "pandas": pandas.__version__, "repeats": args.repeats,
              "harness_sha256": sha256(__file__),
              "sources": {key: {"path": str(getattr(args, key)), "sha256": sha256(getattr(args, key))}
                          for key in ("baseline", "candidate")}, "profiles": {}}
    with tempfile.TemporaryDirectory(prefix="gg-gff-tasks-benchmark-") as temporary:
        root = Path(temporary)
        for name in args.profiles or PROFILES:
            scope = PROFILES[name][0]
            fixture = root / name
            fasta, directory = make_fixture(fixture, name)
            before_inputs = {str(path.relative_to(fixture)): sha256(path) for path in fixture.rglob("*") if path.is_file()}
            runs = {"baseline": [], "candidate": []}
            for trial in range(1 if name.startswith("edge_") else args.repeats + 1):
                payloads = {}
                for key in ("baseline", "candidate") if trial % 2 == 0 else ("candidate", "baseline"):
                    output_dir = root / f"{name}-{trial}-{key}"
                    payload = root / f"{name}-{key}.pkl"
                    command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--source",
                               str(getattr(args, key)), "--profile", name, "--scope", scope, "--fasta", str(fasta),
                               "--gff-dir", str(directory), "--output-dir", str(output_dir), "--payload", str(payload)]
                    environment = os.environ.copy()
                    environment.update({"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                                        "PYTHONHASHSEED": "0"})
                    completed = subprocess.run(command, check=True, capture_output=True, text=True, env=environment)
                    measurement = json.loads(completed.stdout)
                    with payload.open("rb") as handle:
                        payloads[key] = pickle.load(handle)
                    if trial or name.startswith("edge_"):
                        runs[key].append(measurement)
                before, after = payloads["baseline"], payloads["candidate"]
                for key in ("result", "error", "bytes"):
                    assert before.get(key) == after.get(key), (name, key)
                # Worker output directories differ; no content inside the GFF
                # processor's logs depends on output filenames in these calls.
                for key in ("stdout", "stderr"):
                    if name != "edge_threads":  # thread scheduling can interleave progress prints
                        assert before[key] == after[key], (name, key, before[key], after[key])
                    else:
                        assert sorted(before[key].splitlines()) == sorted(after[key].splitlines()), (name, key)
                if "frame" in before:
                    pandas.testing.assert_frame_equal(before["frame"], after["frame"], check_exact=True)
                if name in {"edge_qualifiers", "edge_threads"}:
                    assert before["error"] is None, (name, before["error"])
                    assert len(before["frame"]) == PROFILES[name][1] * PROFILES[name][2], name
                if name in {"edge_duplicate", "edge_malformed", "edge_no_matches"}:
                    assert before["error"] is not None, name
                if not name.startswith("edge_"):
                    assert before["error"] is None, (name, before["error"])
                    assert len(before["result"]) == PROFILES[name][1] if scope == "plan" else len(before["frame"]) > 0
            after_inputs = {str(path.relative_to(fixture)): sha256(path) for path in fixture.rglob("*") if path.is_file()}
            assert before_inputs == after_inputs, name
            report["profiles"][name] = {"scope": scope, "species": PROFILES[name][1], "genes_per_species": PROFILES[name][2],
                                        "runs": runs, "input_hashes": before_inputs}
            print(name, {key: statistics.median(row["wall_seconds"] for row in rows) for key, rows in runs.items()}, flush=True)
    Path(args.output).write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--candidate", type=Path, default=SUPPORT / "gff2genestat.py")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--profiles", nargs="+", choices=PROFILES)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--source")
    parser.add_argument("--scope", choices=("plan", "main"))
    parser.add_argument("--profile")
    parser.add_argument("--fasta")
    parser.add_argument("--gff-dir")
    parser.add_argument("--output-dir")
    parser.add_argument("--payload")
    arguments = parser.parse_args()
    if arguments.worker:
        worker(arguments)
    else:
        if arguments.baseline is None or arguments.output is None:
            parser.error("--baseline and --output are required")
        main(arguments)
