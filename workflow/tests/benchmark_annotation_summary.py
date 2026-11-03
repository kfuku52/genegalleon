"""Compare annotation import and multi-branch summaries on synthetic inputs.

The fixture contains no project data. Run each version in a fresh process, on
the same source DB/TSVs, and compare every output byte as well as peak RSS.
"""
from __future__ import annotations

import argparse
import cProfile
import csv
import hashlib
import json
import pstats
import resource
import sqlite3
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

SUPPORT = Path(__file__).resolve().parents[1] / "support"


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fixture(root, families, products, distinct_go=False):
    sys.path.insert(0, str(SUPPORT))
    import gene_annotation_database as annotations
    import reconciliation_database as reconciliation

    taxa = [f"species_{index:02d}" for index in range(24)]
    catalogue = [f"OG{index:07d}" for index in range(families)]
    model, family_contract = "a" * 64, reconciliation.digest(catalogue)
    branches = [{"branch_id": f"branch_{index}", "descendants": taxa[:index + 1]}
                for index in range(15)]
    dbpath = root / "source.db"
    with sqlite3.connect(dbpath) as db:
        db.executescript("""
            CREATE TABLE reconciliation_meta (key TEXT PRIMARY KEY,value TEXT NOT NULL);
            CREATE TABLE reconciliation_family (orthogroup TEXT PRIMARY KEY,status TEXT NOT NULL);
            CREATE TABLE reconciliation_gene (gene_id TEXT PRIMARY KEY,orthogroup TEXT,species TEXT);
            CREATE INDEX reconciliation_gene_family ON reconciliation_gene(orthogroup,species);
            CREATE TABLE reconciliation_species_branch (branch_id TEXT PRIMARY KEY,descendant_taxa TEXT,terminal INTEGER);
            CREATE TABLE reconciliation_event (orthogroup TEXT,event_type TEXT,branch_id TEXT);
            CREATE INDEX reconciliation_event_branch ON reconciliation_event(branch_id,orthogroup,event_type);
        """)
        db.executemany("INSERT INTO reconciliation_meta VALUES (?,?)", [
            ("schema", reconciliation.SCHEMA), ("family_catalog_sha256", family_contract),
            ("gene_model_contract_sha256", model)])
        db.executemany("INSERT INTO reconciliation_family VALUES (?,?)", [(family, "imported") for family in catalogue])
        descendants = {tuple([species]) for species in taxa} | {tuple(item["descendants"]) for item in branches}
        db.executemany("INSERT INTO reconciliation_species_branch VALUES (?,?,?)", [
            (reconciliation.branch_id(species), json.dumps(species), int(len(species) == 1))
            for species in sorted(descendants)])
        db.executemany("INSERT INTO reconciliation_gene VALUES (?,?,?)", (
            (f"{family}_{species}", family, species) for family in catalogue for species in taxa))
        db.executemany("INSERT INTO reconciliation_event VALUES (?,?,?)", (
            (family, "loss" if index % 3 else "duplication", reconciliation.branch_id(branch["descendants"]))
            for branch in branches for index, family in enumerate(catalogue) if index % 5 != 0))
    sources = []
    fields = ["gene_id", "orthogroup", "sprot_recname", *annotations.GO_COLUMNS.values(),
              *[f"annotation_{index:02d}" for index in range(32)]]
    for species_index, species in enumerate(taxa):
        path = root / (species + ".tsv")
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(fields)
            for index, family in enumerate(catalogue):
                product = f"Enzyme {species_index} with literal annotation, family {family}" if products else ""
                writer.writerow([f"{family}_{species}", family, product,
                                 ";".join(f"GO:{go:07d}" for go in range(1, 7)),
                                 ";".join(["P"] * 6),
                                 ";".join([f"synthetic GO term {index} {species}" if distinct_go else "synthetic GO term"] * 6),
                                 ";".join(["synthetic"] * 6),
                                 *[f"value {index % 17} column {column}" for column in range(32)]])
        sources.append(dict(species=species, path=str(path), sha256=sha(path)))
    (root / "annotations.json").write_text(json.dumps(dict(
        schema=annotations.SCHEMA, gene_model_contract_sha256=model,
        family_catalog_sha256=family_contract, source_database_sha256=sha(dbpath), sources=sources)))
    inputs = {}
    path = root / "branches.json"
    path.write_text(json.dumps(branches))
    inputs["branches"] = dict(path=str(path), sha256=sha(path))
    path = root / "fits.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["orthogroup", "status", "logFC", "p_value"])
        for index, family in enumerate(catalogue):
            writer.writerow([family, "success", index % 4, (index + 1) / (families + 1)])
    for view in ("transcript", "cds"):
        inputs["expression_" + view] = dict(path=str(path), sha256=sha(path))
    (root / "config.json").write_text(json.dumps(dict(
        schema="reconciliation_summary_v1", family_catalog_sha256=family_contract,
        gene_model_contract_sha256=model, inputs=inputs,
        expression_p_value_semantics="one_sided_trap_enrichment_positive_logfc_v1")))


def measure(args):
    sys.path[:0] = [str(args.support.resolve()), str(SUPPORT)]
    import gene_annotation_database as annotations
    import reconciliation_summary as summary

    root, target = args.worker, args.worker / args.result.stem
    target.mkdir()
    profiler = cProfile.Profile() if args.profile else None
    if profiler:
        profiler.enable()
    began = time.perf_counter()
    annotations.build_database(root / "source.db", root / "annotations.json", target / "annotated.db", target)
    annotation_seconds = time.perf_counter() - began
    began = time.perf_counter()
    summary.build_results(target / "annotated.db", root / "config.json", target / "summary")
    summary_seconds = time.perf_counter() - began
    if profiler:
        profiler.disable()
        pstats.Stats(profiler).sort_stats("cumulative").print_stats(25)
    result = dict(annotation_seconds=annotation_seconds, summary_seconds=summary_seconds,
                  seconds=annotation_seconds + summary_seconds,
                  peak_rss_kib_linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  output_sha256={str(path.relative_to(target)): sha(path) for path in
                                 [target / "annotated.db", *sorted((target / "summary").iterdir())]})
    args.result.write_text(json.dumps(result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-support", type=Path)
    parser.add_argument("--families", type=int, default=2048)
    parser.add_argument("--samples", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--empty-products", action="store_true")
    parser.add_argument("--distinct-go", action="store_true", help="Control workload without repeated GO tuples")
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--support", type=Path)
    parser.add_argument("--result", type=Path)
    args = parser.parse_args()
    if args.worker:
        measure(args)
        return
    if not args.output or not 1 <= args.families <= 16384 or not 1 <= args.samples <= 10 or not 0 <= args.warmups <= 3:
        parser.error("Specify output, families 1..16384, samples 1..10 and warmups 0..3")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    implementations = [("baseline", args.baseline_support)] if args.baseline_support else []
    implementations.append(("candidate", SUPPORT))
    result = dict(families=args.families, species=24, branches=15, products=not args.empty_products, distinct_go=args.distinct_go,
                  samples=args.samples, warmups=args.warmups, measurements=[])
    with tempfile.TemporaryDirectory(prefix="gg-annotation-summary-benchmark-") as temporary:
        root = Path(temporary)
        fixture(root, args.families, not args.empty_products, args.distinct_go)
        for iteration in range(-args.warmups, args.samples):
            for label, support in implementations:
                destination = root / f"{label}-{iteration}.json"
                argv = [sys.executable, str(Path(__file__).resolve()), "--worker", str(root),
                        "--support", str(support.resolve()), "--result", str(destination)]
                if args.profile:
                    argv.append("--profile")
                with args.output.with_name(f"{args.output.stem}-{label}-{iteration}.log").open("w") as handle:
                    subprocess.run(argv, check=True, stdout=handle, stderr=subprocess.STDOUT)
                row = dict(label=label, iteration=iteration, **json.loads(destination.read_bytes()))
                if iteration >= 0:
                    result["measurements"].append(row)
                print(json.dumps(row), flush=True)
    if len({json.dumps(row["output_sha256"], sort_keys=True) for row in result["measurements"]}) != 1:
        raise AssertionError("Before/after complete output bytes differ")
    result["medians"] = {label: {key: statistics.median(
        row[key] for row in result["measurements"] if row["label"] == label)
        for key in ("annotation_seconds", "summary_seconds", "seconds", "peak_rss_kib_linux")}
        for label, _ in implementations}
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
