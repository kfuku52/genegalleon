#!/usr/bin/env python3
"""Compare legacy/optimized searches on a bounded, read-only sample of real rescue evidence."""
import argparse
import gzip
import hashlib
import itertools
import json
import platform
import re
import resource
import shutil
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
import rescue_gene_models as rescue  # noqa: E402


def signature(paths):
    digest = hashlib.sha256()
    for path in paths:
        digest.update(path.read_bytes())
    return digest.hexdigest()


def measured(label, directory, operation, outputs):
    directory.mkdir()
    start = time.perf_counter()
    operation(directory)
    seconds = time.perf_counter() - start
    result = {"label": label, "seconds": seconds, "output_sha256": signature(outputs(directory)),
              "parent_peak_rss_KiB_linux": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    print(json.dumps(result), flush=True)
    return result


class FrozenEvidence:
    """Keep producer hashes immutable and account for every original interval."""
    def __init__(self, root, species):
        self.source = root / "rescued" / species
        self.hashes = {}
        self.plan = self.read_json(root / "plan.json")
        self.receipt = self.read_json(self.source / "receipt.json")
        key = self.receipt["key"]
        if key.get("plan") != self.hashes[str(root / "plan.json")] or key.get("species") != species:
            raise ValueError("Producer receipt does not belong to this plan/species")
        self.files = self.receipt["files"]

    def read_json(self, path):
        data = path.read_bytes()
        self.hashes[str(path)] = hashlib.sha256(data).hexdigest()
        return json.loads(data)

    def checked(self, path):
        relative = str(path.relative_to(self.source))
        expected = self.files.get(relative)
        if expected is None or not path.is_file() or rescue.digest(path) != expected:
            raise ValueError("Benchmark source changed or absent from receipt: " + str(path))
        self.hashes[str(path)] = expected
        return path

    def interval_ids(self):
        ids = set()
        for relative in self.files:
            if relative.startswith("intervals/"):
                parts = relative.split("/")
                if len(parts) < 3 or not re.fullmatch(r"[1-9][0-9]*", parts[1]):
                    raise ValueError("Malformed interval in producer receipt: " + relative)
                ids.add(int(parts[1]))
        if sorted(ids) != list(range(1, len(ids) + 1)):
            raise ValueError("Noncontiguous intervals in producer receipt")
        directory = self.source / "intervals"
        if not directory.is_dir() or {p.name for p in directory.iterdir()} != {str(i) for i in ids}:
            raise ValueError("Original interval directories differ from producer receipt")
        for i in sorted(ids):
            for name in ("region.fa", "queries.fa", "models.gff"):
                relative = f"intervals/{i}/{name}"
                if relative not in self.files or not (self.source / relative).is_file():
                    raise ValueError("Original interval file missing: " + relative)
        return sorted(ids)

    def recheck(self):
        self.interval_ids()
        for path, expected in self.hashes.items():
            if not Path(path).is_file() or rescue.digest(path) != expected:
                raise ValueError("Benchmark source changed: " + path)


def require_identical(actual, expected, message):
    if actual != expected:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", required=True, type=Path)
    parser.add_argument("--species", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpus", type=int, default=8)
    parser.add_argument("--interval-count", type=int, default=256)
    parser.add_argument("--query-count", type=int, default=4096)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--check-existing", action="store_true", help="Check optimized searches over ALL intervals/queries against the original receipt; no legacy rerun or warmup")
    args = parser.parse_args()
    if min(args.cpus, args.interval_count, args.query_count, args.repeats) < 1:
        parser.error("Counts must be positive")
    args.output = args.output.resolve()
    args.evidence = args.evidence.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    evidence = FrozenEvidence(args.evidence, args.species)
    plan, source, receipt = evidence.plan, evidence.source, evidence.receipt
    count = len(evidence.interval_ids())
    selected = count if args.check_existing else min(count, args.interval_count)
    chosen = sorted({1 + int((count - 1) * i / max(1, selected - 1)) for i in range(selected)})
    windows, proteins, dna = {}, {"donor": {}}, {}
    for i in chosen:
        directory = source / "intervals" / str(i)
        region, query = evidence.checked(directory / "region.fa"), evidence.checked(directory / "queries.fa")
        evidence.checked(directory / "models.gff")
        sequence = next(rescue.fasta_records(region))[2]
        dna[str(i)] = sequence
        regions = []
        for identifier, _, protein in rescue.fasta_records(query):
            proteins["donor"][identifier] = protein
            regions.append({"id": identifier, "donor": "donor", "query": identifier})
        windows[(str(i), 0, len(sequence))] = regions
    class Genome:
        def fetch(self, seqid, start, end):
            return dna[seqid][start:end]
    result = {"species": args.species, "cpus": args.cpus, "intervals": len(windows),
              "repeats": args.repeats, "warmups": 1, "platform": platform.platform(), "python": sys.version,
              "implementation_sha256": rescue.digest(rescue.__file__),
              "source_plan_sha256": evidence.hashes[str(args.evidence / "plan.json")],
              "source_receipt_sha256": evidence.hashes[str(source / "receipt.json")], "source_hashes": evidence.hashes,
              "original_intervals": count, "skipped": {}, "genome_queries": 0, "unique_queries": 0,
              "tools": rescue.identities(), "samples": {"interval_serial": [], "interval_parallel": [],
                                                        "genome_full": [], "genome_unique": []}}
    result["check_existing"] = args.check_existing
    def interval_outputs(directory):
        return [directory / "intervals" / str(i) / "models.gff" for i in range(1, len(windows) + 1)]
    if not windows:
        result["skipped"]["interval"] = "No interval search recorded by producer"
    for trial in range((1 if args.check_existing else args.repeats + 1) if windows else 0):
        pairs = []
        cases = [("interval_parallel", None)] if args.check_existing else [("interval_serial", 1), ("interval_parallel", None)]
        for label, workers in cases:
            sample = measured(label, args.output / (label + "_" + str(trial)),
                              lambda d, workers=workers: rescue.search_intervals(d, windows, proteins, Genome(),
                                  plan["request"]["sources"][args.species]["genetic_code"],
                                  plan["request"]["parameters"]["max_intron"], args.cpus, workers), interval_outputs)
            pairs.append(sample)
            if trial or args.check_existing:
                result["samples"][label].append(sample)
        if args.check_existing:
            directory = args.output / "interval_parallel_0"
            for index, original_index in enumerate(chosen, 1):
                expected = receipt["files"][f"intervals/{original_index}/models.gff"]
                require_identical(rescue.digest(directory / "intervals" / str(index) / "models.gff"), expected,
                                  f"Interval {original_index} differs")
        else:
            require_identical(pairs[0]["output_sha256"], pairs[1]["output_sha256"], "Interval outputs differ")
    fallback = [name in receipt["files"] for name in ("unresolved.fa", "genome.gff")]
    if any(fallback) and not all(fallback):
        raise ValueError("Incomplete fallback evidence in producer receipt")
    if not any(fallback):
        if any((source / name).exists() for name in ("unresolved.fa", "genome.gff")):
            raise ValueError("Fallback evidence absent from producer receipt")
        result["skipped"]["genome"] = "No fallback search recorded by producer"
    else:
        benchmark_genome(args, evidence, result)
    evidence.recheck()
    result["medians_seconds"] = {label: statistics.median(s["seconds"] for s in samples)
                                 for label, samples in result["samples"].items() if samples}
    if args.check_existing:
        result.update(repeats=1, warmups=0)
    result["outputs_identical"] = True
    (args.output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["medians_seconds"]), flush=True)


def benchmark_genome(args, evidence, result):
    plan, source, receipt = evidence.plan, evidence.source, evidence.receipt
    unresolved = evidence.checked(source / "unresolved.fa")
    evidence.checked(source / "genome.gff")
    queries = list(itertools.islice(rescue.fasta_records(unresolved), None if args.check_existing else args.query_count))
    if not queries:
        raise ValueError("Recorded fallback search has no queries")
    fixture = args.output / "genome_fixture"
    fixture.mkdir()
    original_genome = Path(plan["request"]["sources"][args.species]["genome"])
    genome_hash = plan["request"]["files"][str(original_genome)]
    if rescue.digest(original_genome) != genome_hash:
        raise ValueError("Frozen genome changed")
    evidence.hashes[str(original_genome)] = genome_hash
    genome = fixture / "genome.fa"
    opener = gzip.open if original_genome.suffix == ".gz" else open
    with opener(original_genome, "rb") as handle, genome.open("wb") as out:
        shutil.copyfileobj(handle, out)
    full = fixture / "full.fa"
    full.write_text("".join(f">{identifier}\n{sequence}\n" for identifier, _, sequence in queries))
    regions = [{"id": identifier, "donor": "donor", "query": identifier} for identifier, _, _ in queries]
    proteins = {"donor": {identifier: sequence for identifier, _, sequence in queries}}
    code = plan["request"]["sources"][args.species]["genetic_code"]
    index = fixture / "genome.mpi"
    rescue.run(["miniprot", "-T", code, "-t", args.cpus, "-d", index, genome], fixture, "index")
    result.update(genome_sha256=genome_hash, query_sample_sha256=rescue.digest(full), genome_queries=len(queries),
                  unique_queries=len({sequence for _, _, sequence in queries}))
    def genome_search(directory, dedup):
        if dedup:
            query = directory / "unique.fa"
            mapping = rescue.write_unique_queries(regions, proteins, query)
            raw = directory / "unique.gff"
        else:
            query, raw = full, directory / "models.gff"
        rescue.run(["miniprot", "-u", "-t", args.cpus, "-G", plan["request"]["parameters"]["max_intron"],
                    "--gff", index, query], directory, "map", raw)
        if dedup:
            rescue.expand_miniprot_queries(raw, directory / "models.gff", mapping)
    for trial in range(1 if args.check_existing else args.repeats + 1):
        pairs = []
        cases = [("genome_unique", True)] if args.check_existing else [("genome_full", False), ("genome_unique", True)]
        for label, dedup in cases:
            directory = args.output / (label + "_" + str(trial))
            sample = measured(label, directory, lambda d, dedup=dedup: genome_search(d, dedup),
                              lambda d: [d / "models.gff"])
            match = re.search(r"Peak RSS:\s*([\d.]+) GB", (directory / "logs/map.log").read_text())
            sample["miniprot_peak_rss_GB"] = float(match[1]) if match else None
            pairs.append(sample)
            if trial or args.check_existing:
                result["samples"][label].append(sample)
        if args.check_existing:
            require_identical(rescue.digest(args.output / "genome_unique_0/models.gff"), receipt["files"]["genome.gff"],
                              "Full genome evidence differs")
        else:
            require_identical(pairs[0]["output_sha256"], pairs[1]["output_sha256"], "Genome outputs differ")


if __name__ == "__main__":
    main()
