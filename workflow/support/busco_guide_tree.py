#!/usr/bin/env python3
"""Preserve pre-rescue BUSCO proteins and infer a cached protein k-mer/NJ guide.

Distances compare the same single-copy BUSCO marker, with equal marker weights.
They are guide distances, not calibrated substitutions per site or divergence times.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import subprocess
import tempfile
import time
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path

import numpy as np
from Bio import Phylo

try:
    from busco_reference_quality import COMPARABLE_QUALITY, busco_quality, patristic_distances, safe_token
    from fasta_sequence_store import exclusive_lock, fasta_records, open_text
    from input_generation_array_state import FreshDigestBatch, atomic_json, digest, digest_paths
    from input_generation_stage_resume import copy_atomic, reject_output_overlap
    from species_labeling import extract_species_label
except ImportError:
    from .busco_reference_quality import COMPARABLE_QUALITY, busco_quality, patristic_distances, safe_token
    from .fasta_sequence_store import exclusive_lock, fasta_records, open_text
    from .input_generation_array_state import FreshDigestBatch, atomic_json, digest, digest_paths
    from .input_generation_stage_resume import copy_atomic, reject_output_overlap
    from .species_labeling import extract_species_label

SCHEMA = 1
# Earlier cache producers did not fence concurrent cache/tool changes before
# storing results. Keep those files, but never reuse them under this contract.
CACHE_SCHEMA = 2


def token(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


@contextmanager
def attempt_directory(output):
    path = Path(tempfile.mkdtemp(prefix=".guide-", dir=output))
    try:
        yield path
    except BaseException:
        path.rename(output / (".failed" + path.name))
        raise
    else:
        shutil.rmtree(path)


def cache_matches(path, metadata):
    try:
        return path.is_file() and json.loads(metadata.read_text()).get("sha256") == digest(path)
    except (OSError, ValueError, AttributeError):
        return False


def reject_cache_overlap(destinations, protected):
    # Atomic replacement protects other hard links; resolving names also
    # protects an input symlink whose target would otherwise be replaced.
    for path in destinations:
        if str(path.resolve()) in protected:
            raise ValueError("BUSCO cache output paths overlap input/tool files: " + str(path))


def single_copy_rows(path):
    groups = defaultdict(list)
    with open_text(Path(path)) as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                cells = line.rstrip("\n").split("\t")
                if len(cells) < 2 or cells[1] not in {"Complete", "Duplicated", "Fragmented", "Missing"}:
                    raise ValueError("Invalid BUSCO full-table row: " + str(path))
                groups[cells[0]].append(cells)
    if not groups:
        raise ValueError("Empty BUSCO table: " + str(path))
    result = {}
    for marker, rows in groups.items():
        safe_token(marker, "BUSCO marker")
        if len(rows) == 1 and rows[0][1] == "Complete":
            if len(rows[0]) < 3 or not rows[0][2]:
                raise ValueError("Single-copy BUSCO lacks sequence ID")
            result[marker] = rows[0][2]
    return result, set(groups)


def validate_protein_id(identifier, source_id, mode):
    # BUSCO HMMER reports the inner gene ID of MetaEuk's
    # reference|contig:start-end|strand target, or the base ID of a
    # contig:start-end|strand target in the prokaryotic transcriptome pipeline.
    # Preserve both IDs exactly;
    # checking this relation does not infer coordinates or translate CDS.
    parts = identifier.split("|")
    gene_id = "|".join(parts[1:-1]) if len(parts) > 2 else parts[0]
    if (mode in {"transcriptome", "proteins"} and identifier != source_id
            and not (mode == "transcriptome" and len(parts) > 1 and gene_id == source_id)):
        raise ValueError("BUSCO protein ID does not match its full-table sequence ID: " + source_id)


def preserve(args):
    """Save exact BUSCO AA, including MetaEuk's predicted sequence IDs."""
    paths = [args.full, args.short, args.input]
    boundary = FreshDigestBatch()
    before = boundary.read(paths)
    safe_token(args.species, "species")
    reject_output_overlap([args.output], paths)
    rows, universe = single_copy_rows(args.full)
    quality = busco_quality(args.short)
    if len(universe) != quality["markers"]:
        raise ValueError("BUSCO full/short marker counts differ")
    run_dir = args.run_dir if args.run_dir.is_dir() else args.run_dir.parent
    sequence_dir = run_dir / "busco_sequences/single_copy_busco_sequences"
    protein_paths = {}
    for marker in sorted(rows):
        paths = [sequence_dir / (marker + suffix) for suffix in (".faa", ".faa.gz")]
        matches = [p for p in paths if p.is_file()]
        if len(matches) != 1:
            raise ValueError(f"Expected one single-copy protein for {marker}: {sequence_dir}")
        protein_paths[marker] = matches[0]
    reject_output_overlap([args.output], list(protein_paths.values()))
    boundary.read(protein_paths.values())
    records = {}
    for marker, source_id in sorted(rows.items()):
        proteins = list(fasta_records(protein_paths[marker]))
        if len(proteins) != 1:
            raise ValueError("Single-copy protein file must contain one record: " + marker)
        identifier, header, sequence = proteins[0]
        sequence = sequence.upper().removesuffix("*")
        if not sequence or "*" in sequence or not sequence.isascii() or not sequence.isalpha():
            raise ValueError("Invalid BUSCO protein: " + marker)
        validate_protein_id(identifier, source_id, quality["mode"])
        records[marker] = {"sequence": sequence, "busco_sequence_id": source_id,
                           "protein_id": identifier, "protein_header": header}
    payload = {"schema": SCHEMA, "species": args.species, "quality": quality,
               "source_hashes": before, "records": records, "marker_ids": sorted(universe)}
    boundary.check()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=".single-copy-", dir=args.output.parent)
    os.close(fd)
    try:
        with open(name, "wb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as handle:
            handle.write(json.dumps(payload, sort_keys=True).encode())
        boundary.check()
        os.replace(name, args.output)
    finally:
        Path(name).unlink(missing_ok=True)


def read_archive(path, species, full, short, cds, source_hashes=None, expected_marker_ids=None):
    with gzip.open(path, "rt") as handle:
        data = json.load(handle)
    rows, universe = single_copy_rows(full)
    quality = busco_quality(short)
    if expected_marker_ids is not None and universe != expected_marker_ids:
        raise ValueError("Guide BUSCO marker universe differs: " + species)
    if (data.get("schema") != SCHEMA or data.get("species") != species or data.get("quality") != quality
            or set(data.get("records", {})) != set(rows) or set(data.get("marker_ids", [])) != universe
            or len(universe) != quality["markers"]):
        raise ValueError("BUSCO protein archive does not match tables: " + species)
    # Input paths can change when verified outputs are imported into another run.
    expected = data.get("source_hashes", {})
    current = source_hashes if source_hashes is not None else digest_paths([full, short, cds])
    if sorted(expected.values()) != sorted(current.values()):
        raise ValueError("BUSCO protein archive inputs changed: " + species)
    records = data["records"]
    for marker, record in records.items():
        sequence = record.get("sequence", "")
        if (record.get("busco_sequence_id") != rows[marker] or not sequence or "*" in sequence
                or not sequence.isascii() or not sequence.isalpha() or sequence != sequence.upper()):
            raise ValueError("Invalid archived BUSCO protein: " + species + "/" + marker)
        identifier = record.get("protein_id")
        header = record.get("protein_header")
        if (not isinstance(identifier, str) or not identifier or not isinstance(header, str)
                or not header.split() or header.split()[0] != identifier):
            raise ValueError("Invalid archived BUSCO protein ID/header: " + species + "/" + marker)
        validate_protein_id(identifier, rows[marker], quality["mode"])
    return records, quality


def execute(command, root, label):
    (root / (label + ".command.json")).write_text(json.dumps([str(x) for x in command]) + "\n")
    start = time.perf_counter()
    with (root / (label + ".log")).open("w") as log:
        subprocess.run([str(x) for x in command], stdout=log, stderr=log, check=True)
    return time.perf_counter() - start


def nj(matrix, names, root, label, rapidnj):
    aliases = {f"GG{i:07d}": n for i, n in enumerate(names)}
    phylip = root / (label + ".phylip")
    with phylip.open("w") as handle:
        handle.write(str(len(names)) + "\n")
        for i, alias in enumerate(aliases):
            handle.write(alias + " " + " ".join(format(x, ".17g") for x in matrix[i]) + "\n")
    output = root / (label + ".raw.nwk")
    seconds = execute([rapidnj, phylip, "-i", "pd", "-o", "t", "-n", "-c", "1", "-x", output], root, label)
    tree = Phylo.read(output, "newick")
    leaves = tree.get_terminals()
    if len(leaves) != len(names) or {c.name for c in leaves} != set(aliases):
        raise ValueError("RapidNJ output has missing/duplicated/unknown species")
    for clade in tree.find_clades():
        if clade is not tree.root and (clade.branch_length is None or not math.isfinite(clade.branch_length)
                                      or clade.branch_length < 0):
            raise ValueError("RapidNJ output has invalid branch lengths")
    for leaf in leaves:
        leaf.name = aliases[leaf.name]
    if not any((c.branch_length or 0) > 0 for c in tree.find_clades() if c is not tree.root):
        raise ValueError("No informative guide-tree branch lengths; provide an external tree")
    Phylo.write(tree, root / (label + ".nwk"), "newick", format_branch_length="%1.17g")
    return tree, seconds


def neighbors(tree, names, count, qualities):
    # Precompute distances once, avoiding repeated tree traversal in sorting.
    distances = patristic_distances(tree)
    return {n: sorted((m for m in names if m != n),
                      key=lambda m: (distances[n][m], -qualities[m]["complete_pct"], m))[:count]
            for n in names}


def build(args):
    if (args.markers < 20 or args.k < 3 or args.k > 10 or args.sketch_size < 16
            or args.sketch_size > 65536 or args.cpus < 1 or args.minimum_shared < 2
            or not 0 < args.occupancy <= 1 or args.nearest < 0):
        raise ValueError("Invalid guide-tree parameters")
    kernel, rapidnj = shutil.which("gg-kmer-distance"), shutil.which("rapidnj")
    if not kernel or not rapidnj:
        raise ValueError("BUSCO guide tree requires gg-kmer-distance and rapidnj in the GeneGalleon runtime")
    sources = {}
    for path in sorted(args.cds_dir.iterdir()):
        if not path.is_file() or not path.name.removesuffix(".gz").endswith((".fa", ".fas", ".fasta", ".fna")):
            continue
        name = extract_species_label(path.name)
        safe_token(name, "species")
        if name in sources:
            raise ValueError("Duplicate CDS species: " + name)
        sources[name] = {"cds": str(path.resolve()),
                         "full": str((args.full_dir / (name + ".busco.full.tsv")).resolve()),
                         "short": str((args.short_dir / (name + ".busco.short.txt")).resolve()),
                         "archive": str((args.full_dir / "single_copy" / (name + ".json.gz")).resolve())}
    names = sorted(sources)
    if len(names) < 3:
        raise ValueError("BUSCO guide tree needs at least three species")
    files = [p for source in sources.values() for p in source.values()]
    missing = [p for p in files if not Path(p).is_file()]
    if missing:
        raise ValueError("Missing pre-rescue BUSCO proteins/tables; regenerate BUSCO to preserve its exact AA: " + missing[0])
    boundary = FreshDigestBatch()
    hashes = boundary.read(files)
    support = {name: function.__globals__["__file__"] for name, function in
               (("reader", fasta_records), ("state", digest), ("selection", patristic_distances),
                ("publication", copy_atomic), ("species_labels", extract_species_label))}
    tools_boundary = FreshDigestBatch()
    tool_hashes = tools_boundary.read([kernel, rapidnj, __file__, *support.values()])
    protected = {str(Path(p).resolve()) for p in [*hashes, *tool_hashes]}
    request = {"schema": SCHEMA, "sources": sources, "files": hashes,
               "parameters": {k: getattr(args, k) for k in
                              ("markers", "k", "sketch_size", "occupancy", "minimum_shared", "nearest")},
               "tools": {"kernel": tool_hashes[kernel], "rapidnj": tool_hashes[rapidnj],
                         "implementation": tool_hashes[__file__],
                         "support": {name: tool_hashes[path] for name, path in support.items()}}}
    args.output.mkdir(parents=True, exist_ok=True)
    with exclusive_lock(args.output / ".guide.lock"):
        receipt = args.output / "receipt.json"
        if receipt.exists():
            old = json.loads(receipt.read_text())
            if old.get("request") != request:
                raise ValueError("Frozen BUSCO guide differs; use a new guide output directory")
            saved_files = old.get("files")
            if (not isinstance(saved_files, dict)
                    or not {"guide_tree.nwk", "stability.json", "markers.json"} <= set(saved_files)
                    or any(not p or Path(p).name != p or p in {".", ".."} for p in saved_files)
                    or not all((args.output / p).is_file() and digest(args.output / p) == h for p, h in saved_files.items())):
                raise ValueError("Frozen BUSCO guide outputs changed")
            boundary.check()
            tools_boundary.check()
            return old
        with attempt_directory(args.output) as root:
            start = time.perf_counter()
            data, qualities, counts = {}, {}, defaultdict(int)
            _, marker_universe = single_copy_rows(sources[names[0]]["full"])
            for name in names:
                source = sources[name]
                data[name], qualities[name] = read_archive(source["archive"], name, source["full"], source["short"], source["cds"],
                                                         {source[k]: hashes[source[k]] for k in ("full", "short", "cds")},
                                                         expected_marker_ids=marker_universe)
                data[name] = {marker: record["sequence"] for marker, record in data[name].items()}
                for marker in data[name]:
                    counts[marker] += 1
            if len({tuple(q[k] for k in COMPARABLE_QUALITY) for q in qualities.values()}) != 1:
                raise ValueError("Guide BUSCO lineage/version/mode/date/marker counts must match")
            minimum = math.ceil(args.occupancy * len(names))
            # High occupancy first; hash ordering avoids BUSCO-ID/order artefacts.
            markers = sorted((g for g, count in counts.items() if count >= minimum),
                             key=lambda g: (-counts[g], token(g)))[:args.markers]
            if len(markers) < args.minimum_shared:
                raise ValueError("Too few high-occupancy single-copy BUSCO markers")
            data = {n: {g: seq for g, seq in records.items() if g in markers} for n, records in data.items()}
            atomic_json(root / "markers.json", {"selected": markers, "occupancy": {g: counts[g] for g in markers}})
            args.cache.mkdir(parents=True, exist_ok=True)
            cache_boundary = FreshDigestBatch()
            sketch_paths, sketch_hashes, reused = [], [], 0
            sketch_seconds = 0.0
            for name in names:
                sequences = [data[name].get(g, "") for g in markers]
                key = token({"markers": markers, "sequences": sequences, "k": args.k,
                             "size": args.sketch_size, "kernel": request["tools"]["kernel"], "cache_schema": CACHE_SCHEMA})
                destination = args.cache / (key + ".bin")
                meta = args.cache / (key + ".json")
                reject_cache_overlap([destination, meta], protected)
                with exclusive_lock(args.cache / (key + ".lock")):
                    if cache_matches(destination, meta):
                        reused += 1
                    else:
                        input_path, output_path = root / "sketch.input", root / "sketch.bin"
                        input_path.write_text(str(len(markers)) + "\n" + "\n".join(sequences) + "\n")
                        sketch_seconds += execute([kernel, "sketch", input_path, output_path, args.k, args.sketch_size], root, "sketch-" + name)
                        tools_boundary.check()
                        copy_atomic(output_path, destination, expected_sha256=digest(output_path))
                        atomic_json(meta, {"key": key, "sha256": digest(destination)})
                    observed = digest(destination)
                    if observed != json.loads(meta.read_text())["sha256"]:
                        raise ValueError("Guide sketch changed after cache verification")
                    sketch_hashes.append(observed)
                sketch_paths.append(str(destination.resolve()))
            observed_sketches = cache_boundary.read(sketch_paths)
            if any(observed_sketches[path] != expected for path, expected in zip(sketch_paths, sketch_hashes, strict=True)):
                raise ValueError("Guide sketch changed after cache verification")
            (root / "sketches.txt").write_text("\n".join(sketch_paths) + "\n")
            matrix_key = token({"sketches": sketch_hashes, "kernel": request["tools"]["kernel"], "cache_schema": CACHE_SCHEMA})
            cached_matrix, matrix_meta = args.cache / (matrix_key + ".tsv"), args.cache / (matrix_key + ".matrix.json")
            reject_cache_overlap([cached_matrix, matrix_meta], protected)
            compare_seconds, compare_reused = 0.0, False
            with exclusive_lock(args.cache / (matrix_key + ".matrix.lock")):
                if cache_matches(cached_matrix, matrix_meta):
                    observed = cache_boundary.read([cached_matrix])[str(cached_matrix)]
                    if observed != json.loads(matrix_meta.read_text())["sha256"]:
                        raise ValueError("Guide distances changed after cache verification")
                    copy_atomic(cached_matrix, root / "pairwise.tsv", expected_sha256=observed)
                    compare_reused = True
                else:
                    compare_seconds = execute([kernel, "compare", root / "sketches.txt", root / "pairwise.tsv", args.cpus], root, "compare")
                    cache_boundary.check()
                    tools_boundary.check()
                    copy_atomic(root / "pairwise.tsv", cached_matrix, expected_sha256=digest(root / "pairwise.tsv"))
                    atomic_json(matrix_meta, {"key": matrix_key, "sha256": digest(cached_matrix)})
            cache_boundary.check()
            matrices = [np.zeros((len(names), len(names)), dtype=float) for _ in range(3)]
            saturated, shared_total, minimum_shared = 0, 0, len(markers)
            saturated_partners = [0] * len(names)
            with (root / "pairwise.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            if len(rows) != len(names)*(len(names)-1)//2:
                raise ValueError("Incomplete pairwise distance matrix")
            seen = set()
            for row in rows:
                i, j = int(row["i"]), int(row["j"])
                if not 0 <= i < j < len(names) or (i,j) in seen:
                    raise ValueError("Invalid pairwise matrix indices")
                seen.add((i,j))
                shared = int(row["shared"])
                if shared < args.minimum_shared or min(int(row["count_0"]), int(row["count_1"])) < max(1, args.minimum_shared//2):
                    raise ValueError(f"Too few shared usable BUSCO markers: {names[i]} / {names[j]}")
                minimum_shared = min(minimum_shared, shared)
                saturated += int(row["saturated"])
                shared_total += shared
                if int(row["saturated"]) == shared:
                    saturated_partners[i] += 1
                    saturated_partners[j] += 1
                for matrix, key in zip(matrices, ("distance", "panel_0", "panel_1"), strict=True):
                    value = float(row[key])
                    if not math.isfinite(value) or not 0 <= value <= 1:
                        raise ValueError("Invalid k-mer distance")
                    matrix[i,j] = matrix[j,i] = value
            if saturated == shared_total:
                raise ValueError("All usable BUSCO marker pairs are k-mer saturated; provide an external tree")
            isolated = [name for name, count in zip(names, saturated_partners, strict=True) if count == len(names)-1]
            if isolated:
                raise ValueError("All comparisons are k-mer saturated for species " + ", ".join(isolated)
                                 + "; provide an external tree")
            trees, nj_seconds = [], 0.0
            for matrix, label in zip(matrices, ("guide_tree", "panel_0", "panel_1"), strict=True):
                tree, seconds = nj(matrix, names, root, label, rapidnj)
                trees.append(tree)
                nj_seconds += seconds
            selected = [neighbors(tree, names, args.nearest, qualities) for tree in trees]
            stability = {n: {"nearest": selected[0][n], "panel_0": selected[1][n], "panel_1": selected[2][n],
                             "stable": set(selected[0][n]) == set(selected[1][n]) == set(selected[2][n])}
                         for n in names}
            atomic_json(root / "stability.json", stability)
            performance = {"total_seconds": time.perf_counter()-start, "sketch_seconds": sketch_seconds,
                           "compare_seconds": compare_seconds, "nj_seconds": nj_seconds,
                           "species": len(names), "markers": len(markers), "sketches_reused": reused,
                           "distances_reused": compare_reused,
                           "minimum_shared": minimum_shared, "saturated_marker_pairs": saturated,
                           "stable_species": sum(v["stable"] for v in stability.values()),
                           "cpus": args.cpus, "rooting": "unrooted", "negative_lengths": "RapidNJ -n adjustment"}
            atomic_json(root / "performance.json", performance)
            boundary.check()
            tools_boundary.check()
            cache_boundary.check()
            products = [p for p in root.iterdir() if p.is_file() and p.name not in {"sketch.input", "sketch.bin", "sketches.txt"}]
            result = {"request": request, "files": {p.name: digest(p) for p in products}, "performance": performance}
            reject_output_overlap([args.output / p.name for p in products] + [receipt], files)
            for path in products:
                copy_atomic(path, args.output / path.name, expected_sha256=result["files"][path.name])
            if any(digest(args.output / p) != h for p, h in result["files"].items()):
                raise ValueError("Guide outputs changed during publication")
            boundary.check()
            tools_boundary.check()
            cache_boundary.check()
            atomic_json(receipt, result)
            return result


def requires_single_copy(manifest, species):
    """Keep a sealed legacy table-only result valid without inventing AA data.

    The caller still verifies every declared input/output with provenance.
    Newly generated runs and existing three-output contracts require AA data.
    """
    if not manifest.exists():
        return True
    value = json.loads(manifest.read_text())
    outputs = value.get("outputs", [])
    roles = [item.get("label") for item in outputs]
    if (value.get("schema_version") != 1 or value.get("step") != "input_generation_species_busco"
            or value.get("family_id") != species
            or [item.get("label") for item in value.get("inputs", [])] != ["species_cds"]
            or len(set(roles)) != len(roles)
            or set(roles) not in ({"busco_full", "busco_short"},
                                  {"busco_full", "busco_short", "busco_single_copy"})):
        raise ValueError("Unsupported BUSCO output contract: " + species)
    return "busco_single_copy" in roles


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    contract = commands.add_parser("output-contract")
    contract.add_argument("--manifest", required=True, type=Path)
    contract.add_argument("--species", required=True)
    export = commands.add_parser("preserve")
    for key in ("run-dir", "full", "short", "input", "output"):
        export.add_argument("--" + key, required=True, type=Path)
    export.add_argument("--species", required=True)
    guide = commands.add_parser("build")
    for key in ("cds-dir", "full-dir", "short-dir", "output", "cache"):
        guide.add_argument("--" + key, required=True, type=Path)
    guide.add_argument("--markers", type=int, default=200)
    guide.add_argument("--k", type=int, default=5)
    guide.add_argument("--sketch-size", type=int, default=256)
    guide.add_argument("--occupancy", type=float, default=.8)
    guide.add_argument("--minimum-shared", type=int, default=50)
    guide.add_argument("--nearest", type=int, default=3)
    guide.add_argument("--cpus", type=int, default=1)
    return result


def main():
    args = parser().parse_args()
    if args.command == "preserve":
        preserve(args)
    elif args.command == "output-contract":
        print(int(requires_single_copy(args.manifest, args.species)))
    else:
        print(json.dumps(build(args)["performance"], sort_keys=True))


if __name__ == "__main__":
    main()
