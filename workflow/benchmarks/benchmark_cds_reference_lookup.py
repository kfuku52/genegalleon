#!/usr/bin/env python3
"""Compare CDS reconstruction against many FASTA references in one runtime."""

import argparse
import hashlib
import importlib.util
import json
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

SUPPORT = Path(__file__).resolve().parents[1] / "support"


def worker(source, root, references, queries):
    sys.path.insert(0, str(SUPPORT))
    import pysam

    spec = importlib.util.spec_from_file_location("measured_normaliser", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    normaliser = module.CdsModelNormaliser(
        {"genome": str(root / "genome.fa"), "gff": str(root / "source.gff")}, root, "benchmark")
    digest = hashlib.sha256()
    started = time.perf_counter()
    try:
        for number in range(queries):
            contig = f"chr{number * 7919 % references:06d}"
            strand = "+" if number % 2 else "-"
            rows = [{"seqid": contig, "strand": strand, "start": start, "end": start + 9}
                    for start in (0, 16, 32, 48)]
            digest.update(normaliser.fetch(rows).encode("ascii"))
    finally:
        normaliser.close()
    return {"seconds": time.perf_counter() - started,
            "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            "output_sha256": digest.hexdigest(), "pysam": pysam.__version__}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-source", type=Path)
    parser.add_argument("--source", type=Path, default=SUPPORT / "cds_model_normalisation.py")
    parser.add_argument("--references", type=int, default=89579)
    parser.add_argument("--queries", type=int, default=1000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker-root", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if min(args.references, args.queries, args.repeats) < 1:
        parser.error("Workload sizes and repetitions must be positive")
    if args.worker_root:
        print(json.dumps(worker(args.source, args.worker_root, args.references, args.queries)))
        return
    if args.output is None or args.output.exists():
        parser.error("--output must name a new file")
    variants = {"current": args.source.resolve()}
    if args.baseline_source:
        variants = {"baseline": args.baseline_source.resolve(), **variants}
    results = {name: [] for name in variants}
    with tempfile.TemporaryDirectory(prefix="gg-reference-lookup-") as temporary:
        root = Path(temporary)
        with (root / "genome.fa").open("w") as handle:
            for number in range(args.references):
                handle.write(f">chr{number:06d}\n" + "ACGT" * 16 + "\n")
        (root / "source.gff").write_text("##gff-version 3\n")
        for trial in range(args.repeats + 1):
            names = list(variants)
            if trial % 2:
                names.reverse()
            for name in names:
                command = [sys.executable, str(Path(__file__).resolve()), "--source", str(variants[name]),
                           "--worker-root", str(root), "--references", str(args.references),
                           "--queries", str(args.queries)]
                result = json.loads(subprocess.check_output(command, text=True))
                if trial:
                    results[name].append(result)
    hashes = {item["output_sha256"] for trials in results.values() for item in trials}
    if len(hashes) != 1:
        raise RuntimeError("CDS reconstruction outputs differ")
    payload = {"references": args.references, "queries": args.queries, "blocks_per_query": 4,
               "warmups": 1, "repeats": args.repeats, "trials": results, "equivalent": True,
               "python": sys.version, "platform": platform.platform(),
               "source_sha256": {name: hashlib.sha256(path.read_bytes()).hexdigest()
                                 for name, path in variants.items()},
               "median_seconds": {name: statistics.median(x["seconds"] for x in trials)
                                  for name, trials in results.items()},
               "scope": "Synthetic fragmented genome; includes first FASTA index/open and CDS fetches; "
                        "fixture generation and GFF parsing excluded; fresh process per trial; Linux RSS."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload))


if __name__ == "__main__":
    main()
