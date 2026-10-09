#!/usr/bin/env python3
"""Compare bounded GFF/CDS extraction against a Git revision in one runtime."""
import argparse
import gzip
import hashlib
import json
import resource
import statistics
import subprocess
import sys
import tarfile
import tempfile
import time
from pathlib import Path


def worker(support, directory, mode):
    sys.path.insert(0, str(support))
    start = time.perf_counter()
    task = dict(provider="direct", species_key="Fixture_species", species_prefix="Fixture_species",
                gff_path=directory / "models.gff", genome_path=directory / "genome.fa")
    digest = hashlib.sha256()
    count = 0
    if mode == "extract":
        from format_species_annotation.genbank import derive_cds_records_from_gff_and_genome
        for header, sequence in derive_cds_records_from_gff_and_genome(task):
            digest.update((header + "\n" + sequence + "\n").encode())
            count += 1
    else:
        from format_species_discovery import format_genome
        output = directory / "formatted"
        output.mkdir(exist_ok=True)
        task["gff_path"] = None
        result = format_genome(task, output, True, False)
        count = result["written"]
        with gzip.open(result["output_path"], "rb") as handle:
            while block := handle.read(1024 * 1024):
                digest.update(block)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(json.dumps(dict(seconds=time.perf_counter() - start,
                          peak_rss_mib=rss / (1024 ** 2 if sys.platform == "darwin" else 1024),
                          sha256=digest.hexdigest(), records=count)))


def fixture(directory, mib, genes):
    # One long chromosome exercises the observed whole-record allocation;
    # variable case and line widths retain parser/output equivalence coverage.
    block = ("aCGtN" * 205)[:1024] + "\n"
    with (directory / "genome.fa").open("w") as handle:
        handle.write(">lcl|chr1 description\n")
        for _ in range(mib * 1024):
            handle.write(block)
    with (directory / "models.gff").open("w") as handle:
        handle.write("##gff-version 3\n")
        for i in range(genes):
            start = i * 400 + 1
            strand = "-" if i % 2 else "+"
            handle.write(f"chr1\tsrc\tgene\t{start}\t{start+299}\t.\t{strand}\t.\tID=g{i}\n")
            handle.write(f"chr1\tsrc\tmRNA\t{start}\t{start+299}\t.\t{strand}\t.\tID=t{i};Parent=g{i}\n")
            for a, b, phase in ((start, start + 98, 0), (start + 200, start + 299, 0)):
                handle.write(f"chr1\tsrc\tCDS\t{a}\t{b}\t.\t{strand}\t{phase}\tParent=t{i}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--revision", required=False)
    parser.add_argument("--mib", type=int, default=256)
    parser.add_argument("--genes", type=int, default=2000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--mode", choices=("extract", "format"), default="extract")
    parser.add_argument("--worker", nargs=2, metavar=("SUPPORT", "FIXTURE"))
    args = parser.parse_args()
    if args.worker:
        worker(*(Path(x) for x in args.worker), args.mode)
        return
    if min(args.mib, args.genes, args.repeats) < 1 or args.genes * 400 > args.mib * 1024 ** 2:
        parser.error("Positive sizes and gene coordinates within the fixture genome are required")
    repo = Path(__file__).resolve().parents[2]
    with tempfile.TemporaryDirectory(prefix="gg-genome-benchmark-") as temp:
        directory = Path(temp)
        fixture(directory, args.mib, args.genes)
        supports = {"after": repo / "workflow/support"}
        if args.revision:
            archive = directory / "before.tar"
            with archive.open("wb") as handle:
                subprocess.run(["git", "archive", args.revision, "workflow/support"], cwd=repo,
                               stdout=handle, check=True)
            before = directory / "before"
            before.mkdir()
            with tarfile.open(archive) as handle:
                handle.extractall(before, filter="data")
            supports = {"before": before / "workflow/support", **supports}
        results = {name: [] for name in supports}
        for trial in range(args.repeats + 1):
            for name, support in (list(supports.items()) if trial % 2 == 0 else list(supports.items())[::-1]):
                run = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker",
                                      str(support), str(directory), "--mode", args.mode], capture_output=True, text=True, check=True)
                value = json.loads(run.stdout)
                if trial:
                    results[name].append(value)
        proofs = {(row["sha256"], row["records"]) for rows in results.values() for row in rows}
        if len(proofs) != 1:
            raise ValueError("Extraction output differs between runs/revisions")
        print(json.dumps(dict(mode=args.mode, genome_mib=args.mib, genes=args.genes, warmups=1,
                              repeats=args.repeats, equivalent=True, python=sys.version,
                              results={name: dict(median_seconds=statistics.median(r["seconds"] for r in rows),
                                                  median_peak_rss_mib=statistics.median(r["peak_rss_mib"] for r in rows),
                                                  trials=rows) for name, rows in results.items()}), indent=2))


if __name__ == "__main__":
    main()
