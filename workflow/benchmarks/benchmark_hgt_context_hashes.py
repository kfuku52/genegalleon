#!/usr/bin/env python3
"""Compare complete transfer-context summaries on exact synthetic reconciliations.

Use the GeneGalleon runtime wrapper. Fixture generation, imports, serialization,
and exact output comparison are outside the measured summarize() call. RSS is
the whole worker peak through summarize(), including fixtures and imports.
Workers are independent processes; warmups precede alternating paired samples.
Small control calls use batches to reduce timer and filesystem noise; timings
are per complete summarize() call, excluding output destruction between calls.
"""

import argparse
import csv
import hashlib
import json
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

PROFILES = {"single": (1, 0), "few_large": (4, 8192),
            "many": (512, 0), "many_large": (512, 8192)}
BATCHES = {"single": 100, "few_large": 3, "many": 3, "many_large": 3}


def leaf(name, species, transfer=False):
    event = f'<transferBack destinationSpecies="{species}"/>' if transfer else ""
    return (f'<clade><name>{name}</name><eventsRec>{event}'
            f'<leaf speciesLocation="{species}"/></eventsRec></clade>')


def fixture(root, profile):
    count, padding = PROFILES[profile]
    clades, branches, genes, stats = [], [], [], []
    for i in range(count):
        donor, recipient = f"donor_{i}", f"recipient_{i}"
        clades.append(f'<clade><name>event_{i}</name><eventsRec>'
                      '<branchingOut speciesLocation="D_species"/></eventsRec>'
                      + leaf(donor, "D_species") + leaf(recipient, "R_species", True) + '</clade>')
        branch = dict(orthogroup="OG0001", branch_id=str(i), node_name=f"event_{i}",
                      generax_transfer="Y@D_species@R_species", candidate_genes=f"{donor}; {recipient}")
        branches.append(branch)
        stats.append(dict(branch_id=str(i), gene_labels=branch["candidate_genes"],
                          support_generax_ufboot=str(80 + i % 21), num_leaf="2", so_event="S"))
        for gene, species in ((donor, "D_species"), (recipient, "R_species")):
            genes.append(dict(orthogroup="OG0001", gene_id=gene, gene_taxon=species,
                              host_scaffold_status="not_available"))
    clades.extend(leaf(f"background_{i}", "D_species" if i % 2 else "R_species") for i in range(padding))
    xml = ('<recPhylo><spTree><phylogeny><clade><name>root</name>'
           '<clade><name>D_species</name></clade><clade><name>R_species</name></clade>'
           '</clade></phylogeny></spTree><recGeneTree><phylogeny><clade><name>genes</name>'
           + ''.join(clades) + '</clade></phylogeny></recGeneTree></recPhylo>')
    (root / "generax_xml").mkdir()
    xml_path = root / "generax_xml/OG0001_generax.xml"
    xml_path.write_text(xml)
    (root / "stat_branch").mkdir()
    stat_path = root / "stat_branch/OG0001_stat.branch.tsv"
    with stat_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(stats[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(stats)
    return branches, genes, xml_path, stat_path


def worker(args):
    sys.path.insert(0, str(args.support_root))
    from gene_family_output_store import GeneFamilyOutputStore
    from summarize_hgt_transfer_context import summarize

    branches, genes, xml, stat = fixture(args.fixture_root, args.profile)
    wall_seconds = cpu_seconds = 0.0
    for _ in range(BATCHES[args.profile]):
        result = None
        started = time.perf_counter()
        started_cpu = time.process_time()
        result = summarize(branches, genes, GeneFamilyOutputStore(args.fixture_root))
        wall_seconds += time.perf_counter() - started
        cpu_seconds += time.process_time() - started_cpu
    wall_seconds /= BATCHES[args.profile]
    cpu_seconds /= BATCHES[args.profile]
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024)
    assert len(result[0]) == PROFILES[args.profile][0]
    assert len(result[1]) == 2 * len(result[0])
    assert all(event["mapping_status"] == "matched" for event in result[0])
    payload = json.dumps(result, sort_keys=True).encode()
    args.worker_output.write_bytes(payload)
    print(json.dumps(dict(wall_seconds=wall_seconds, cpu_seconds=cpu_seconds,
                         peak_rss_bytes=peak_rss, calls_per_sample=BATCHES[args.profile],
                         event_count=len(result[0]), link_count=len(result[1]),
                         xml_bytes=xml.stat().st_size, stat_bytes=stat.stat().st_size,
                         output_sha256=hashlib.sha256(payload).hexdigest())))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-support", type=Path)
    parser.add_argument("--candidate-support", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--profiles", nargs="+", choices=PROFILES, default=list(PROFILES))
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--support-root", type=Path)
    parser.add_argument("--fixture-root", type=Path)
    parser.add_argument("--worker-output", type=Path)
    parser.add_argument("--profile", choices=PROFILES)
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    if not all((args.baseline_support, args.candidate_support, args.output)):
        parser.error("--baseline-support, --candidate-support, --output are required")
    roots = {"baseline": args.baseline_support.resolve(), "candidate": args.candidate_support.resolve()}
    report = dict(python=sys.version, platform=platform.platform(), arguments=sys.argv,
                  implementations={name: dict(root=str(root), sha256=hashlib.sha256(
                      (root / "summarize_hgt_transfer_context.py").read_bytes()).hexdigest())
                                   for name, root in roots.items()}, profiles={})
    with tempfile.TemporaryDirectory(prefix="gg-transfer-context-benchmark-") as temporary:
        for profile in args.profiles:
            samples = {name: [] for name in roots}
            reference = None
            for trial in range(args.repeats + 1):
                order = list(roots) if trial % 2 == 0 else list(reversed(roots))
                for implementation in order:
                    fixture_root = Path(temporary) / f"{profile}-{implementation}-{trial}"
                    fixture_root.mkdir()
                    output = fixture_root / "output.json"
                    raw = subprocess.check_output([
                        sys.executable, str(Path(__file__).resolve()), "--worker",
                        "--support-root", str(roots[implementation]), "--fixture-root", str(fixture_root),
                        "--worker-output", str(output), "--profile", profile], text=True)
                    row = json.loads(raw)
                    payload = output.read_bytes()
                    if reference is None:
                        reference = payload
                    else:
                        assert payload == reference, (profile, implementation, trial)
                    if trial:
                        samples[implementation].append(row)
                    print(profile, implementation, "warmup" if trial == 0 else trial,
                          row["wall_seconds"], row["cpu_seconds"], row["peak_rss_bytes"], flush=True)
            report["profiles"][profile] = samples
            for implementation, rows in samples.items():
                report["profiles"][profile][implementation] = dict(samples=rows, medians={
                    metric: statistics.median(row[metric] for row in rows)
                    for metric in ("wall_seconds", "cpu_seconds", "peak_rss_bytes")})
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
