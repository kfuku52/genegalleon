#!/usr/bin/env python3
"""Compare FASTA retention in complete GFF and cold/warm synteny-cache calls.

Fresh workers begin with stdlib imports only. CLI timing includes real source
imports/main. Cache timing includes importing the actual caller and its complete
ensure_species_gene_cache call, including a real GFF subprocess on cache misses.
Fixture generation, warm-cache preparation and output checking are untimed.
Parent and largest-child RSS are separate, never a process-tree memory sum.
"""

import argparse
import contextlib
import gzip
import hashlib
import importlib.util
import io
import json
import os
import platform
import re
import resource
import runpy
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow/support"
PROFILES = {"small-raw": (32, False), "normal-raw": (30000, False), "normal-gzip": (30000, True)}


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def make_fixture(root, genes, compressed):
    root.mkdir(parents=True)
    fasta = root / ("Species_a.fa.gz" if compressed else "Species_a.fa")
    with fasta.open("wb") as raw:
        with contextlib.ExitStack() as stack:
            stream = stack.enter_context(gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=1)) if compressed else raw
            for index in range(genes):
                stream.write(f">Species_a_gene{index:07d} synthetic CDS header α {index}\n".encode())
                for _ in range(25):
                    stream.write(b"ATG" * 20 + b"\n")
    gff = root / "Species_a.gff"
    with gff.open("w") as handle:
        handle.write("##gff-version 3\n")
        for index in range(genes):
            gene = f"gene{index:07d}"
            start, end = index * 2000 + 1, index * 2000 + 1500
            strand = "+" if index % 2 else "-"
            for kind, phase, attributes in (
                ("gene", ".", f"ID={gene}"),
                ("mRNA", ".", f"ID={gene}.t1;Parent={gene}"),
                ("CDS", "0", f"ID={gene}.c1;Parent={gene}.t1"),
            ):
                handle.write(f"chr1\tfixture\t{kind}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attributes}\n")
    return fasta


def cache_arguments(args):
    return dict(species_name="Species_a", species_cds_path=str(args.fasta), dir_sp_gff=str(args.fasta.parent),
                cache_dir=str(args.cache), lock_dir=str(args.cache / "locks"),
                gff2genestat_script=str(args.source), threads=1)


def import_cache_caller():
    sys.path.insert(0, str(SUPPORT))
    spec = importlib.util.spec_from_file_location("gff_fasta_cache_caller", SUPPORT / "synteny_neighbors.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def normalise_log(text, args):
    text = text.replace(str(args.fasta.parent), "<fixture>").replace(str(args.output_dir), "<output>")
    text = re.sub(r"\d{4}-\d\d-\d\d \d\d:\d\d:\d\d(?:\.\d+)?", "<datetime>", text)
    return re.sub(r"completed in [\d,]+ secs:", "completed in <elapsed> secs:", text)


def output_signature(path, args):
    # This import is after timing and RSS capture in each measured worker.
    import pandas as pd

    frame = pd.read_csv(path, sep="\t")
    frame.to_pickle(args.output_dir / "frame.pkl")
    return {"sha256": digest(path), "rows": len(frame), "columns": frame.columns.tolist(),
            "dtypes": list(map(str, frame.dtypes)), "index_dtype": str(frame.index.dtype),
            "index_name": frame.index.name}


def worker(args):
    args.source = args.source.resolve()
    if args.scope == "prepare-warm":
        module = import_cache_caller()
        module.ensure_species_gene_cache(**cache_arguments(args))
        return
    args.output_dir.mkdir(parents=True)
    os.chdir(args.output_dir)
    stdout, stderr = io.StringIO(), io.StringIO()
    child_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started, cpu_started = time.perf_counter(), time.process_time()
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
        if args.scope == "cli":
            sys.path.insert(0, str(SUPPORT))
            output = args.output_dir / "stats.tsv"
            sys.argv = [str(args.source), "--mode", "gene_delim", "--dir_gff", str(args.fasta.parent),
                        "--seqfile", str(args.fasta), "--feature", "CDS", "--multiple_hits", "longest",
                        "--phase-policy", "report", "--structure-policy", "report", "--ncpu", "1",
                        "--outfile", str(output)]
            runpy.run_path(str(args.source), run_name="__main__")
            reader_calls = 0
        else:
            module = import_cache_caller()
            reader_calls = 0
            actual_run = module.run_cmd

            def counted_run(command):
                nonlocal reader_calls
                reader_calls += 1
                return actual_run(command)

            module.run_cmd = counted_run
            output = Path(module.ensure_species_gene_cache(**cache_arguments(args)))
    wall, parent_cpu = time.perf_counter() - started, time.process_time() - cpu_started
    child_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    parent_peak = next(int(line.split()[1]) for line in Path("/proc/self/status").read_text().splitlines()
                       if line.startswith("VmHWM:"))
    child_peak = int(child_after.ru_maxrss)
    child_cpu = ((child_after.ru_utime + child_after.ru_stime)
                 - (child_before.ru_utime + child_before.ru_stime))
    assert reader_calls == {"cli": 0, "cache-cold": 1, "cache-warm": 0}[args.scope]
    signature = output_signature(output, args)
    receipt = None
    if args.scope != "cli":
        receipt = json.loads(Path(str(output) + ".provenance.json").read_text())
        expected = module.species_gene_cache_contract("Species_a", str(args.fasta),
                                                     str(args.fasta.parent / "Species_a.gff"),
                                                     gff2genestat_script=str(args.source))
        assert receipt == {**expected, "output_sha256": signature["sha256"]}
        assert receipt["inputs"]["gff_reader_sha256"] == digest(args.source)
        # The actual source hash was checked above. Every other receipt field
        # must match across variants, including all input and output hashes.
        receipt["inputs"]["gff_reader_sha256"] = "<verified implementation>"
    (args.output_dir / "result.json").write_text(json.dumps({
        "wall_s": wall, "parent_cpu_s": parent_cpu, "child_cpu_s": child_cpu,
        "total_cpu_s": parent_cpu + child_cpu,
        "parent_peak_rss_kib": parent_peak, "largest_child_rss_kib": child_peak,
        "max_parent_or_child_rss_kib": max(parent_peak, child_peak),
        "reader_subprocesses": reader_calls, "output": signature, "cache_receipt": receipt,
        "stdout": normalise_log(stdout.getvalue(), args), "stderr": normalise_log(stderr.getvalue(), args),
    }, indent=2) + "\n")


def controller(args):
    sources = {"before": args.before.resolve(), "after": args.after.resolve()}
    args.output = args.output.resolve()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    report = {"source_sha256": {name: digest(path) for name, path in sources.items()},
              "benchmark_sha256": digest(Path(__file__)), "python": sys.version,
              "platform": platform.platform(), "pairs": args.pairs, "warmup_pairs": 1,
              "sequence_nt": 1500, "profiles": {}}
    worker_environment = os.environ.copy()
    worker_environment["PYTHONPATH"] = str(SUPPORT) + (
        os.pathsep + worker_environment["PYTHONPATH"] if worker_environment.get("PYTHONPATH") else ""
    )
    with tempfile.TemporaryDirectory(prefix="gg-gff-fasta-benchmark-") as temporary:
        temporary = Path(temporary)
        for profile in args.profiles:
            genes, compressed = PROFILES[profile]
            fasta = make_fixture(temporary / profile, genes, compressed)
            inputs = {path.name: {"bytes": path.stat().st_size, "sha256": digest(path)}
                      for path in sorted(fasta.parent.iterdir())}
            profile_data = report["profiles"][profile] = {"genes": genes, "compressed": compressed,
                                                         "inputs": inputs, "scopes": {}}
            warm_caches = {name: temporary / f"warm-{profile}-{name}" for name in sources}
            if "cache-warm" in args.scopes:
                for name, source in sources.items():
                    subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker", "--source", str(source),
                                    "--scope", "prepare-warm", "--fasta", str(fasta),
                                    "--cache", str(warm_caches[name])], check=True, capture_output=True, text=True,
                                   env=worker_environment)
            for scope in args.scopes:
                pairs = profile_data["scopes"][scope] = []
                for pair in range(-1, args.pairs):
                    measured = {}
                    for name in (("before", "after") if pair % 2 == 0 else ("after", "before")):
                        root = temporary / f"{profile}-{scope}-{pair}-{name}"
                        cache = warm_caches[name] if scope == "cache-warm" else root / "cache"
                        command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--source", str(sources[name]),
                                   "--scope", scope, "--fasta", str(fasta), "--cache", str(cache), "--output-dir", str(root)]
                        completed = subprocess.run(command, capture_output=True, text=True, env=worker_environment)
                        if completed.returncode:
                            print(completed.stdout, end="")
                            print(completed.stderr, end="", file=sys.stderr)
                            completed.check_returncode()
                        measured[name] = json.loads((root / "result.json").read_text())
                    import pandas as pd

                    before = pd.read_pickle(temporary / f"{profile}-{scope}-{pair}-before/frame.pkl")
                    after = pd.read_pickle(temporary / f"{profile}-{scope}-{pair}-after/frame.pkl")
                    pd.testing.assert_frame_equal(before, after, check_exact=True)
                    for key in ("output", "cache_receipt", "stdout", "stderr", "reader_subprocesses"):
                        assert measured["before"][key] == measured["after"][key], (profile, scope, pair, key, measured)
                    if pair >= 0:
                        pairs.append(measured)
                    args.output.write_text(json.dumps(report, indent=2) + "\n")
                    print(json.dumps({"profile": profile, "scope": scope, "pair": pair, **measured}), flush=True)
            shutil.rmtree(fasta.parent)
    summaries = []
    for profile, values in report["profiles"].items():
        for scope, pairs in values["scopes"].items():
            medians = {variant: {metric: statistics.median(pair[variant][metric] for pair in pairs)
                                 for metric in ("wall_s", "parent_cpu_s", "child_cpu_s", "total_cpu_s",
                                                "parent_peak_rss_kib", "largest_child_rss_kib", "max_parent_or_child_rss_kib")}
                       for variant in sources}
            ratios = [pair["before"]["wall_s"] / pair["after"]["wall_s"] for pair in pairs]
            summaries.append({"profile": profile, "scope": scope, **medians, "paired_speedup_range": [min(ratios), max(ratios)]})
    report["summaries"] = summaries
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path)
    parser.add_argument("--after", type=Path, default=SUPPORT / "gff2genestat.py")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--profiles", nargs="+", choices=PROFILES, default=list(PROFILES))
    parser.add_argument("--scopes", nargs="+", choices=["cli", "cache-cold", "cache-warm"], default=["cli", "cache-cold", "cache-warm"])
    parser.add_argument("--pairs", type=int, default=5)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--source", type=Path)
    parser.add_argument("--scope", choices=["cli", "cache-cold", "cache-warm", "prepare-warm"])
    parser.add_argument("--fasta", type=Path)
    parser.add_argument("--cache", type=Path)
    parser.add_argument("--output-dir", type=Path)
    options = parser.parse_args()
    if not options.worker and (not options.before or not options.output or options.pairs < 1):
        parser.error("--before, --output, and positive --pairs are required")
    worker(options) if options.worker else controller(options)
