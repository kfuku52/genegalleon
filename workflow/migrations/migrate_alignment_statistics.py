#!/usr/bin/env python3
"""One-time, offline migration of legacy alignment statistics (raw or ZIP)."""
import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
import artifact_provenance as provenance
import gene_family_output_store as storage
import orthogroup_output_summary
import query2family_output_summary

# Historical names belong only to this one-time migration, never live readers.
LEGACY = "amas"


def legacy_entries(store):
    entries = {}
    for stage in ("original", "cleaned"):
        for subdir, suffix in ((f"{LEGACY}_{stage}", f"_{LEGACY}.{stage}.tsv"),
                               (f"{LEGACY}.{stage}", f".{LEGACY}.{stage}.tsv")):
            for name in store.file_names(subdir):
                if not name.endswith(suffix):
                    raise ValueError(f"Unexpected legacy file: {subdir}/{name}")
                family = name[:-len(suffix)]
                if not family or any(c in family for c in "/\\\t\n"):
                    raise ValueError(f"Invalid family ID: {family!r}")
                entries.setdefault((family, stage), []).append(f"{subdir}/{name}")
    # A stopped cleanup may have deleted the TSV but not yet its manifest.
    for name in store.file_names("artifact_provenance"):
        for stage in ("original", "cleaned"):
            suffix = f".{LEGACY}_{stage}.json"
            if name.endswith(suffix):
                entries.setdefault((name[:-len(suffix)], stage), [])
    return entries


def read_overrides(path, root):
    result = {}
    if path:
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = (row["family_id"], row["stage"])
                if key in result:
                    raise ValueError(f"Duplicate mapping: {key}")
                source = Path(row["alignment"])
                result[key] = (source if source.is_absolute() else root / source, row["seq_type"])
    return result


def alignment_for(store, key, workspace, overrides):
    if key in overrides:
        return overrides[key]
    family, stage = key
    name = f"{family}.{LEGACY}_{stage}.json"
    try:
        with store.open_binary("artifact_provenance", name) as handle:
            manifest = json.load(handle)
    except FileNotFoundError as exc:
        raise ValueError(f"No provenance for {family}/{stage}; supply --alignments TSV. No input guessing is performed.") from exc
    if manifest.get("family_id") != family or manifest.get("step") != f"{LEGACY}_{stage}":
        raise ValueError(f"Unexpected provenance identity: {name}")
    inputs = [entry for entry in manifest["inputs"] if entry["label"] == "alignment"]
    if len(inputs) != 1:
        raise ValueError(f"No unique alignment in {name}")
    return provenance.resolve_reference(inputs[0], store.root, workspace), manifest["parameters"]["data_type"]


def copy_source(store, source, destination):
    try:
        relative = source.relative_to(store.root)
    except ValueError:
        relative = None
    if relative is not None and len(relative.parts) == 2:
        with store.open_binary(*relative.parts) as inp, destination.open("wb") as out:
            shutil.copyfileobj(inp, out)
    else:
        shutil.copyfile(source, destination)


def migrate(args):
    root = args.root.resolve()
    workspace = args.workspace_root.resolve()
    store = storage.GeneFamilyOutputStore(root)
    overrides = read_overrides(args.alignments, root)
    # The exclusive family gate blocks workflow producers for the entire migration.
    with storage._bucket_lock(storage.family_gate_path(store.archive_root), exclusive=True, nonblocking=True) as idle:
        if not idle:
            raise RuntimeError("A workflow is using this output store; stop it before migrating.")
        entries = legacy_entries(store)
        unknown = {
            (family, stage) for family, stage in set(overrides) - set(entries)
            if not (
                store.logical_exists(f"alignment_stats_{stage}/{family}_alignment_stats.{stage}.tsv")
                and store.logical_exists(f"artifact_provenance/{family}.alignment_stats_{stage}.json")
            )
        }
        if unknown:
            raise ValueError(f"Mappings without legacy output: {sorted(unknown)}")
        with tempfile.TemporaryDirectory(prefix="alignment-statistics-") as temporary:
            tmp = Path(temporary)
            prepared = []
            for index, key in enumerate(sorted(entries)):
                family, stage = key
                source, seq_type = alignment_for(store, key, workspace, overrides)
                if seq_type not in {"dna", "aa"}:
                    raise ValueError(f"Invalid sequence type for {key}: {seq_type}")
                job = tmp / str(index)
                job.mkdir()
                compressed = job / ("source.gz" if source.name.endswith(".gz") else "source.fa")
                copy_source(store, source, compressed)
                fasta = job / f"{family}.alignment_stats.{stage}.input.fasta"
                subprocess.run(["seqkit", "seq", "--threads", "1", str(compressed), "--out-file", str(fasta)], check=True)
                generated = job / "summary.tsv"
                subprocess.run(["cdskit", "stats", "--mode", "alignment", "--seq_type", seq_type,
                                "--seq_file", str(fasta), "--out_file", str(generated)], check=True)
                with generated.open() as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                if len(rows) != 1 or int(rows[0]["No_of_taxa"]) < 1:
                    raise ValueError(f"Invalid statistics output for {key}")
                prepared.append((key, source, seq_type, generated, compressed))
            # All alignments must validate before any existing artifact is changed.
            for (family, stage), source, seq_type, generated, snapshot in prepared:
                output = root / f"alignment_stats_{stage}" / f"{family}_alignment_stats.{stage}.tsv"
                output.parent.mkdir(parents=True, exist_ok=True)
                if output.is_symlink() or output.parent.is_symlink() or (root / "artifact_provenance").is_symlink():
                    raise ValueError(f"Symlinked output: {output}")
                with tempfile.NamedTemporaryFile(dir=output.parent, delete=False) as handle:
                    pending = Path(handle.name)
                    handle.write(generated.read_bytes())
                os.replace(pending, output)
                contract_args = provenance.build_parser().parse_args([
                    "record", "--manifest", str(root / "artifact_provenance" / f"{family}.alignment_stats_{stage}.json"),
                    "--step", f"alignment_stats_{stage}", "--family-id", family,
                    "--logical-root", str(root), "--workspace-root", str(workspace),
                    "--input", f"alignment={snapshot}", "--output", f"alignment_stats={output}",
                    "--parameter", f"data_type={seq_type}",
                    "--parameter", "statistics_engine=cdskit-stats-alignment",
                ])
                contract = provenance.build_contract(contract_args, include_diagnostics=True)
                contract["inputs"][0].update(provenance.path_reference(source, root, workspace))
                provenance.write_manifest_atomic(contract_args.manifest, contract)
            for (family, stage), old_paths in entries.items():
                for old in [*old_paths, f"artifact_provenance/{family}.{LEGACY}_{stage}.json"]:
                    store.delete(old, family_id=family, _family_locked=True)
            for stage in ("original", "cleaned"):
                for subdir in (f"{LEGACY}_{stage}", f"{LEGACY}.{stage}"):
                    directory = root / subdir
                    if directory.is_dir() and not directory.is_symlink() and not any(directory.iterdir()):
                        directory.rmdir()
    # Purge takes its own exclusive family gate and removes deleted ZIP members.
    storage.purge_archives(root, args.mode)
    current_store = storage.GeneFamilyOutputStore(root)
    for stage in ("original", "cleaned"):
        for subdir in (f"{LEGACY}_{stage}", f"{LEGACY}.{stage}"):
            directory = current_store.payload_root / subdir
            if directory.is_dir() and not directory.is_symlink() and not any(directory.iterdir()):
                directory.rmdir()
    summary = args.summary_out.resolve()
    summary.parent.mkdir(parents=True, exist_ok=True)
    if args.mode == "orthogroup":
        orthogroup_output_summary.run(SimpleNamespace(
            dir_og=str(root), genecount=str(args.genecount), out=str(summary), ncpu=1,
            updated_genecount_out=str(summary.with_name(summary.stem + ".genecount.alignment_stats.tsv")),
        ))
    else:
        query2family_output_summary.run(SimpleNamespace(
            dir_query2family=str(root), dir_query_gene=str(args.dir_query_gene), out=str(summary), ncpu=1,
        ))
    if args.mode == "orthogroup":
        old_aggregates = {
            summary.parent / f"orthogroup_genecount.{LEGACY}.tsv",
            summary.with_name(summary.stem + f".genecount.{LEGACY}.tsv"),
            args.genecount.with_name(args.genecount.stem + f".{LEGACY}.tsv"),
        }
        for old in old_aggregates:
            if old.is_file() and not old.is_symlink() and old.resolve() not in {summary, args.genecount.resolve()}:
                old.unlink()
    print(f"Migrated {len(entries)} statistics tables; summary: {summary}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--workspace-root", required=True, type=Path)
    parser.add_argument("--mode", required=True, choices=["orthogroup", "query2family"])
    parser.add_argument("--summary-out", required=True, type=Path)
    parser.add_argument("--genecount", type=Path)
    parser.add_argument("--dir-query-gene", type=Path)
    parser.add_argument("--alignments", type=Path, help="Optional TSV: family_id, stage, alignment, seq_type. Overrides old provenance.")
    args = parser.parse_args(argv)
    if args.mode == "orthogroup" and (args.genecount is None or not args.genecount.is_file()):
        parser.error("orthogroup mode requires an existing --genecount")
    if args.mode == "query2family" and (args.dir_query_gene is None or not args.dir_query_gene.is_dir()):
        parser.error("query2family mode requires an existing --dir-query-gene")
    migrate(args)


if __name__ == "__main__":
    main()
