#!/usr/bin/env python3
"""Inventory and apply reviewed calibrations without querying or reinterpreting TimeTree."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

FIELDS = (
    "topology_sha256", "node_id", "descendant_tips", "calibration", "decision",
    "source_type", "source_references", "interval_kind", "time_unit",
    "dependency_groups", "sequence_overlap", "calibration_overlap",
    "node_basis", "distribution_basis", "reviewer", "review_note",
)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def json_text(value) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def node_id(node) -> str:
    return digest(json_text(sorted(node.leaf_names())).encode())


def load_tree(path, *, require_binary=True):
    from nwkit.util import read_tree

    tree = read_tree(str(path), "auto", True, quiet=True)
    names = list(tree.leaf_names())
    if not names or any(not name for name in names) or len(set(names)) != len(names):
        raise ValueError("Calibration trees require unique nonempty tip labels.")
    if require_binary and any(len(node.children) != 2 for node in tree.traverse() if not node.is_leaf):
        raise ValueError("Calibration trees must be rooted and binary.")
    return tree


def topology_digest(tree) -> str:
    return digest(json_text(sorted(node_id(node) for node in tree.traverse())).encode())


def calibration_digest(tree) -> str:
    from nwkit.time_tree import parse_mcmctree_calibration

    labels = {}
    for node in tree.traverse():
        if not node.is_leaf:
            record = parse_mcmctree_calibration(node.name)
            if record:
                labels[node_id(node)] = {key: value for key, value in record.items() if key != "raw"}
    return digest(json_text(labels).encode())


def tree_text(tree) -> str:
    from nwkit.util import write_tree

    stream = io.StringIO()
    write_tree(tree, SimpleNamespace(outfile=stream), format=1, quiet=True, props=[])
    return stream.getvalue().strip() + "\n"


def table_text(rows) -> str:
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def atomic_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                     delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(text)
    try:
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def inventory(path, outdir, manifest=None):
    """Snapshot observed labels; never infer historical source/review information."""
    from nwkit.time_tree import parse_mcmctree_calibration

    path, outdir = Path(path), Path(outdir)
    tree = load_tree(path, require_binary=False)
    topology = topology_digest(tree)
    rows = []
    for node in tree.traverse():
        if node.is_leaf:
            continue
        record = parse_mcmctree_calibration(node.name)
        row = dict.fromkeys(FIELDS, "")
        row.update(topology_sha256=topology, node_id=node_id(node),
                   descendant_tips=json_text(sorted(node.leaf_names())),
                   calibration=record["raw"] if record else "", decision="unreviewed", source_type="unknown",
                   interval_kind="unknown", time_unit="Ma", dependency_groups="[]",
                   sequence_overlap="unknown", calibration_overlap="unknown")
        rows.append(row)
    # These immutable snapshots are separate from the legacy artifact contract.
    snapshot = outdir / digest(path.read_bytes())
    if not snapshot.exists():
        outdir.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=outdir) as work:
            staging = Path(work) / "snapshot"
            staging.mkdir()
            (staging / "candidates.tsv").write_text(table_text(rows), encoding="utf-8")
            shutil.copyfile(path, staging / "observed.nwk")
            (staging / "status.json").write_text(json_text({
                "schema_version": 1, "status": "unreviewed",
                "tree_sha256": digest(path.read_bytes()), "topology_sha256": topology,
                "observed_utc": datetime.now(timezone.utc).isoformat(),
                "candidate_count": len(rows),
                "source_provenance": "unknown; labels alone do not establish TimeTree sources",
                "historical_mcmc_diagnostics": "unknown",
            }) + "\n", encoding="utf-8")
            staging.rename(snapshot)
    current = {"schema_version": 1, "status": "unreviewed", "tree_sha256": digest(path.read_bytes()),
               "topology_sha256": topology, "calibration_sha256": calibration_digest(tree),
               "inventory": str(snapshot / "candidates.tsv"), "mcmc_diagnostics": "unknown"}
    reviewed_path = path.parent / "reviewed_calibrations.json"
    if manifest is not None and reviewed_path.is_file():
        reviewed = json.loads(reviewed_path.read_text(encoding="utf-8"))
        if (reviewed.get("topology_sha256") == topology
                and reviewed.get("calibration_sha256") == current["calibration_sha256"]
                and reviewed.get("manifest_sha256") == digest(Path(manifest).read_bytes())):
            current.update(status="reviewed_inputs", reviewed_record=str(reviewed_path),
                           scientific_validation="not_established_by_manifest")
    atomic_text(outdir / "current_status.json", json_text(current) + "\n")
    print(f"Calibration evidence status: {current['status']}; candidate template: {snapshot / 'candidates.tsv'}")
    return snapshot


def reviewed_rows(manifest, tree):
    from nwkit.time_tree import parse_mcmctree_calibration

    with Path(manifest).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(FIELDS):
            raise ValueError("Unexpected calibration manifest columns; start with inventory candidates.tsv.")
        rows = list(reader)
    topology = topology_digest(tree)
    nodes = {node_id(node): node for node in tree.traverse() if not node.is_leaf}
    seen, accepted = set(), []
    for row in rows:
        if None in row or any(value is None for value in row.values()):
            raise ValueError("Malformed calibration TSV row.")
        key = row["node_id"]
        if key in seen or key not in nodes or row["topology_sha256"] != topology:
            raise ValueError(f"Duplicate node, unknown node, or changed topology: {key}")
        seen.add(key)
        if json.loads(row["descendant_tips"]) != sorted(nodes[key].leaf_names()):
            raise ValueError(f"Descendant tips do not match node {key}.")
        if row["decision"] not in {"accept", "exclude", "unreviewed"}:
            raise ValueError(f"Invalid decision for {key}.")
        if row["decision"] != "accept":
            continue
        for field in ("source_references", "node_basis", "distribution_basis", "reviewer", "review_note"):
            if row[field].strip().lower() in {"", "unknown", "na", "n/a"}:
                raise ValueError(f"Accepted calibration requires {field}: {key}")
        if row["source_type"] not in {"primary", "secondary", "timetree"}:
            raise ValueError(f"Accepted calibration requires a known source_type: {key}")
        if row["interval_kind"] not in {"empirical_rule", "min_max", "study_posterior", "fossil_bounds", "user_defined"}:
            raise ValueError(f"Accepted calibration requires a known interval_kind: {key}")
        if row["time_unit"] != "Ma":
            raise ValueError("Reviewed calibration ages must use the public unit Ma.")
        groups = json.loads(row["dependency_groups"])
        if (not isinstance(groups, list) or not groups
                or any(not isinstance(group, str) or not group.strip() for group in groups)
                or len(set(groups)) != len(groups)):
            raise ValueError(f"Accepted calibration needs nonempty dependency_groups: {key}")
        for field in ("sequence_overlap", "calibration_overlap"):
            if row[field] not in {"none", "shared", "partial", "unknown"}:
                raise ValueError(f"Invalid {field}: {key}")
        record = parse_mcmctree_calibration(row["calibration"])
        if record is None or record["type"] not in {"bounded", "lower", "upper"}:
            raise ValueError("Reviewed mode supports explicit B/L/U distributions only.")
        if not row["calibration"].startswith(("B(", "L(", "U(")):
            raise ValueError("Use explicit B/L/U distributions with tail probabilities.")
        if record["type"] == "bounded":
            if record["lower"] >= record["upper"]:
                raise ValueError("B bounds must enclose a nonzero interval.")
            if record["lower_tail"] + record["upper_tail"] >= 1:
                raise ValueError("B tail probabilities must sum to less than one.")
        if "upper" in record and record["upper"] <= 0:
            raise ValueError("Upper age must be positive.")
        accepted.append(row)
    if not accepted:
        raise ValueError("Reviewed manifest accepts no calibrations.")
    return rows, accepted


def apply_manifest(tree_path, manifest, outfile, audit_out):
    inputs = {Path(tree_path).resolve(), Path(manifest).resolve()}
    if Path(outfile).resolve() in inputs or Path(audit_out).resolve() in inputs:
        raise ValueError("Outputs must not replace input files.")
    if Path(outfile).resolve() == Path(audit_out).resolve():
        raise ValueError("Tree and audit outputs must differ.")
    tree = load_tree(tree_path)
    rows, accepted = reviewed_rows(manifest, tree)
    labels = {row["node_id"]: row["calibration"] for row in accepted}
    for node in tree.traverse():
        if not node.is_leaf:
            node.name = labels.get(node_id(node), "")
    result = tree_text(tree)
    status = {
        "schema_version": 1, "status": "reviewed_inputs",
        "scientific_validation": "not_established_by_manifest",
        "mcmc_diagnostics": "not_run",
        "topology_sha256": topology_digest(tree),
        "calibration_sha256": calibration_digest(tree),
        "manifest_sha256": digest(Path(manifest).read_bytes()),
        "tree_sha256": digest(result.encode()), "accepted_count": len(accepted),
        "reviewed_utc": datetime.now(timezone.utc).isoformat(), "candidates": rows,
    }
    atomic_text(audit_out, json_text(status) + "\n")
    atomic_text(outfile, result)


def archive_existing(species_dir, provenance_dir=None):
    """Keep the previous dating inputs and outputs before an explicit rebuild."""
    species_dir = Path(species_dir)
    sources = [species_dir / name for name in
               ("constrained_tree", "mcmctree_parameter_estimation", "mcmctree_main")]
    sources = [path for path in sources if path.exists()]
    if not sources:
        return
    parent = species_dir / "calibration_history"
    parent.mkdir(exist_ok=True)
    destination = Path(tempfile.mkdtemp(prefix=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S-"), dir=parent))
    for source in sources:
        shutil.copytree(source, destination / source.name)
    if provenance_dir is not None:
        provenance = destination / "artifact_provenance"
        provenance.mkdir()
        for source in Path(provenance_dir).glob("species_tree.*.json"):
            shutil.copyfile(source, provenance / source.name)
    (destination / "status.json").write_text(json_text({
        "status": "archived_before_calibration_rebuild",
        "historical_validation": "not_inferred",
    }) + "\n", encoding="utf-8")
    print(f"Preserved previous dating artifacts: {destination}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    inv = sub.add_parser("inventory")
    inv.add_argument("--tree", required=True, type=Path)
    inv.add_argument("--outdir", required=True, type=Path)
    inv.add_argument("--manifest", type=Path, help="Currently selected reviewed manifest, if any")
    apply = sub.add_parser("apply")
    apply.add_argument("--tree", required=True, type=Path)
    apply.add_argument("--manifest", required=True, type=Path)
    apply.add_argument("--outfile", required=True, type=Path)
    apply.add_argument("--audit-out", required=True, type=Path)
    archive = sub.add_parser("archive")
    archive.add_argument("--species-dir", required=True, type=Path)
    archive.add_argument("--provenance-dir", type=Path)
    args = parser.parse_args()
    if args.command == "inventory":
        inventory(args.tree, args.outdir, args.manifest)
    elif args.command == "apply":
        apply_manifest(args.tree, args.manifest, args.outfile, args.audit_out)
    else:
        archive_existing(args.species_dir, args.provenance_dir)


if __name__ == "__main__":
    main()
