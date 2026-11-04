"""Compare complete audited DB builds on one synthetic raw/ZIP/SQLite input.

Use --baseline-support against a saved support tree. Every measured build runs
in its own process, retains the final content fence, and compares the complete
logical DB, audit report and audit inventory. This is not a real-project ETA.
"""
from __future__ import annotations

import argparse
import cProfile
import hashlib
import importlib.util
import json
import pstats
import resource
import statistics
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

SUPPORT = Path(__file__).resolve().parents[1] / "support"


def fixture(root, support, families, packed, files_per_family=10, plot_arguments=False):
    sys.path.insert(0, str(support))
    import metadata_store
    from gene_family_output_store import archive_completed_outputs

    logical = root / "output/orthogroup"
    shared = root / "shared.tsv"
    shared.write_bytes(b"shared input\n" * 100)
    manifests = logical / "artifact_provenance"
    manifests.mkdir(parents=True)
    catalog = {f"OG{index:07d}" for index in range(families)}
    for family in sorted(catalog):
        contents = {
            "stat_tree/" + family + "_stat.tree.tsv": b"num_branch\tnum_spe\tnum_dup\tnum_sp\n1\t2\t0\t1\n",
            "stat_branch/" + family + "_stat.branch.tsv": b"branch_id\tnode_name\tnum_sp\tso_event\n0\tg1\t1\tS\n",
            **{f"fixture_{index}/{family}.txt": (f"{family} fixture {index}\n" * 32).encode()
               for index in range(files_per_family)},
        }
        outputs = []
        for relative, body in contents.items():
            path = logical / relative
            path.parent.mkdir(exist_ok=True)
            path.write_bytes(body)
            outputs.append(dict(label=relative.split("/")[0], scope="logical", path=relative,
                                artifact_type="file", sha256=hashlib.sha256(body).hexdigest(), size_bytes=len(body)))
        for step in range(4):
            payload = dict(schema_version=1, family_id=family, step=f"example_{step}",
                           inputs=[dict(label="shared", scope="workspace", path="shared.tsv", artifact_type="file",
                                        sha256=hashlib.sha256(shared.read_bytes()).hexdigest(), size_bytes=shared.stat().st_size)],
                           outputs=outputs)
            (manifests / f"{family}.example_{step}.json").write_text(json.dumps(payload))
        if plot_arguments:
            from gene_tree_plot_config import record
            snapshot = manifests / f'{family}.tree_plot.args.json'
            record(snapshot, logical, ['--panel1=tree'], 'taxonomic')
            payload = dict(schema_version=1, family_id=family, step='tree_plot', inputs=[], outputs=[
                dict(label='tree_plot_arguments', scope='logical', path='artifact_provenance/' + snapshot.name,
                     artifact_type='file', sha256=hashlib.sha256(snapshot.read_bytes()).hexdigest(),
                     size_bytes=snapshot.stat().st_size)])
            (manifests / f'{family}.tree_plot.json').write_text(json.dumps(payload))
    if packed:
        metadata_store.pack(manifests, writers_stopped=True, delete_originals=True)
        archive_completed_outputs(logical, "orthogroup", catalog, lambda name: name.split("_")[0].split(".")[0],
                                  include_incomplete=True)
    return logical


def measure(args):
    sys.path.insert(0, str(args.support))
    import zipfile

    import artifact_audit_runtime
    import gene_family_database_pipeline as pipeline
    import metadata_store
    from gene_family_output_store import GeneFamilyOutputStore

    spec = importlib.util.spec_from_file_location("db_benchmark", Path(__file__).with_name("benchmark_orthogroup_database.py"))
    fingerprint = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fingerprint)
    root, target = args.worker, args.result
    logical = root / "output/orthogroup"
    global_root = root / "output/.gg_global_artifacts"
    global_root.mkdir(exist_ok=True)
    database = target.with_suffix(".db")
    report = target.with_suffix(".tsv")
    manifest = target.with_suffix(".manifest.json")
    common = ["--workspace-root", str(root), "--logical-root", str(logical)]
    audit = [*common, "--output-tsv", str(report), "--mode", "orthogroup", "--workers", "4", "--progress-interval", "0"]
    db = ["--overwrite", "1", "--dbpath", str(database), "--dir_gene_family", str(logical),
          "--dir_stat_tree", str(logical / "stat_tree"), "--dir_stat_branch", str(logical / "stat_branch")]
    record = ["--workspace-root", str(root), "--logical-root", str(global_root), "--family-id", "orthogroup",
              "--step", "gene_family_database", "--manifest", str(manifest),
              "--input-gene-family-store", "families=" + str(logical), "--output", "database=" + str(database)]
    counts, lock = {"metadata_headers": 0, "zip_opens": 0, "index_epochs": 0}, threading.Lock()
    original_signature, original_zip = metadata_store._db_signature, zipfile.ZipFile.__init__
    original_epoch = GeneFamilyOutputStore._read_index_epoch
    def signature(*a, **kw):
        with lock:
            counts["metadata_headers"] += 1
        return original_signature(*a, **kw)
    def archive(*a, **kw):
        with lock:
            counts["zip_opens"] += 1
        return original_zip(*a, **kw)
    def epoch(*a, **kw):
        with lock:
            counts["index_epochs"] += 1
        return original_epoch(*a, **kw)
    metadata_store._db_signature, zipfile.ZipFile.__init__ = signature, archive
    GeneFamilyOutputStore._read_index_epoch = epoch
    phases = []
    original_phase = artifact_audit_runtime.AuditProgress.phase
    def phase(progress, name, total=None):
        phases.append((name, time.monotonic()))
        return original_phase(progress, name, total)
    artifact_audit_runtime.AuditProgress.phase = phase
    profile = cProfile.Profile() if args.profile else None
    if profile:
        profile.enable()
    began = time.monotonic()
    assert pipeline.run(audit, db, record) == 0
    seconds = time.monotonic() - began
    if profile:
        profile.disable()
        profile.dump_stats(str(target.with_suffix(".prof")))
        pstats.Stats(profile).sort_stats("cumulative").print_stats(20)
    phases.append(("finished", time.monotonic()))
    output_hash, tables = fingerprint.fingerprint(database)
    result = json.loads(report.with_suffix(".tsv.result.json").read_bytes())
    measured = dict(seconds=seconds, peak_rss_kib_linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                    logical_sha256=output_hash, table_rows=tables, counts=counts,
                    audit_sha256=hashlib.sha256(report.read_bytes()).hexdigest(),
                    inventory_sha256=result["inventory_sha256"], status_counts=result["status_counts"],
                    audit_metrics=result["metrics"],
                    phases=[dict(phase=name, seconds=phases[i+1][1]-stamp)
                            for i, (name, stamp) in enumerate(phases[:-1])])
    target.write_text(json.dumps(measured, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-support", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--families", type=int, default=128)
    parser.add_argument("--files-per-family", type=int, default=10, help="Unique extra artifacts, excluding two statistics")
    parser.add_argument("--plot-arguments", action="store_true", help="Include real renderer argument outputs and manifests")
    parser.add_argument("--samples", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--layout", choices=["raw", "packed"], default="packed")
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--support", type=Path, default=SUPPORT)
    parser.add_argument("--result", type=Path)
    args = parser.parse_args()
    if args.worker:
        measure(args)
        return
    if (not args.output or not 1 <= args.families <= 20000 or not 1 <= args.files_per_family <= 32
            or not 1 <= args.samples <= 10
            or not 0 <= args.warmups <= 3):
        parser.error("--output, families 1..20000, files-per-family 1..32, samples 1..10 and warmups 0..3 required")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    implementations = [("baseline", args.baseline_support)] if args.baseline_support else []
    implementations.append(("candidate", SUPPORT))
    result = dict(families=args.families, files_per_family=args.files_per_family, plot_arguments=args.plot_arguments,
                  layout=args.layout, workers=4, samples=args.samples,
                  warmups=args.warmups, measurements=[])
    with tempfile.TemporaryDirectory(prefix="gg-audited-db-benchmark-") as temporary:
        root = Path(temporary)
        fixture(root, SUPPORT, args.families, args.layout == "packed", args.files_per_family, args.plot_arguments)
        for iteration in range(-args.warmups, args.samples):
            for label, support in (implementations if iteration % 2 == 0 else implementations[::-1]):
                target = root / f"{label}-{iteration}.json"
                argv = [sys.executable, str(Path(__file__).resolve()), "--worker", str(root),
                        "--support", str(support.resolve()), "--result", str(target)]
                if args.profile:
                    argv.append("--profile")
                with args.output.with_name(f"{args.output.stem}-{label}-{iteration}.log").open("w") as handle:
                    subprocess.run(argv, check=True, cwd=root, stdout=handle, stderr=subprocess.STDOUT)
                row = dict(label=label, iteration=iteration, **json.loads(target.read_bytes()))
                if iteration >= 0:
                    result["measurements"].append(row)
                args.output.write_text(json.dumps(result, indent=2) + "\n")
                print(json.dumps(row), flush=True)
    for key in ("logical_sha256", "audit_sha256", "inventory_sha256"):
        if len({row[key] for row in result["measurements"]}) != 1:
            raise AssertionError("Before/after output mismatch: " + key)
    for key in ("unique_digest_reads", "digest_revalidation_reads", "digest_bytes_read", "source_guard_count"):
        if len({row["audit_metrics"][key] for row in result["measurements"]}) != 1:
            raise AssertionError("Before/after verification coverage mismatch: " + key)
    result["medians"] = {label: {key: statistics.median(row[key] for row in result["measurements"] if row["label"] == label)
                                  for key in ("seconds", "peak_rss_kib_linux")} for label, _ in implementations}
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
