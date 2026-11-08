"""Method-separated, lossless reconciliation bundles and species-edge counts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
from collections import Counter
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from reconciliation_database import branch_id, canonical, digest, parse_family, reference_branches

SCHEMA = "species_reconciliation_v1"
COLUMNS = (
    "orthogroup",
    "source_id",
    "method",
    "species_branch_id",
    "status",
    "loss_count",
    "duplication_count",
    "speciation_count",
    "transfer_in_count",
    "transfer_out_count",
    "extant_gene_count",
    "extant_species_count",
)


def input_bytes(path):
    return Path(path).read_bytes() if path and Path(path).is_file() else None


def source_record(family, method, inputs, metadata):
    identity = dict(
        schema=SCHEMA,
        orthogroup=family,
        method=method,
        inputs={k: hashlib.sha256(v).hexdigest() for k, v in inputs.items() if v is not None},
        metadata=metadata,
    )
    return dict(identity, source_id=digest(identity), status="analysed", clades=[], events=[], genes=[])


def lca_records(source, table_data, gene_tree_data, species_tree_data, reference):
    from nwkit.clade_index import CladeIndex
    from nwkit.util import read_tree

    species = read_tree(species_tree_data.decode(), 1, True, quiet=True)
    gene_tree = read_tree(gene_tree_data.decode(), 1, True, quiet=True)
    species_index, gene_index = CladeIndex(species), CladeIndex(gene_tree)
    species_ids = {
        species_index.clade_id_for_node(n): branch_id(species_index.names_for_node(n)) for n in species.traverse()
    }
    expected_parents = {
        gene_index.clade_id_for_node(n): None if n.is_root else gene_index.clade_id_for_node(n.up)
        for n in gene_tree.traverse()
    }
    reader = csv.DictReader(io.StringIO(table_data.decode()), delimiter="\t")
    required = {
        "tree_id",
        "gene_clade_id",
        "parent_gene_clade_id",
        "event_type",
        "event_source",
        "event_status",
        "mapping_status",
        "species_event_id",
        "gene_name",
        "species_name",
        "implied_losses",
        "implied_loss_species_clade_ids",
    }
    if not required <= set(reader.fieldnames or []):
        raise ValueError("LCA table lacks loss locations; regenerate with the qualified NWKIT runtime")
    rows = list(reader)
    keys = [r["gene_clade_id"] for r in rows]
    if len(keys) != len(set(keys)) or set(keys) != set(expected_parents):
        raise ValueError("LCA table and rooted gene tree clades disagree")
    ordinal = {key: i for i, key in enumerate(keys)}
    family = source["orthogroup"]
    root_branch = None
    for row in rows:
        if (
            row["tree_id"] != family
            or row["event_source"] != "lca"
            or row["mapping_status"] != "mapped"
            or row["event_type"] not in {"leaf", "speciation", "duplication"}
            or row["event_status"]
            not in ({"not-applicable", "resolved"} if row["event_type"] == "leaf" else {"resolved"})
        ):
            raise ValueError("LCA events must be completely mapped and resolved")
        key = row["gene_clade_id"]
        parent = row["parent_gene_clade_id"] or None
        if parent != expected_parents[key]:
            raise ValueError("LCA gene parent differs from rooted input tree")
        mapped = species_ids.get(row["species_event_id"])
        if mapped is None:
            raise ValueError("LCA species mapping differs from complete reference tree")
        if parent is None:
            root_branch = mapped
        clade = ordinal[key]
        source["clades"].append((family, clade, ordinal[parent] if parent else None, row["gene_name"] or key))
        source["events"].append(
            (
                family,
                clade,
                0,
                row["event_type"],
                row["species_name"],
                mapped,
                canonical(dict(gene_clade_id=key, species_event_id=row["species_event_id"])),
            )
        )
        lost = json.loads(row["implied_loss_species_clade_ids"])
        if (
            not isinstance(lost, list)
            or any(k not in species_ids for k in lost)
            or int(row["implied_losses"]) != len(lost)
        ):
            raise ValueError("LCA loss locations and totals disagree")
        for i, location in enumerate(lost, 1):
            source["events"].append(
                (
                    family,
                    clade,
                    i,
                    "loss",
                    "",
                    species_ids[location],
                    canonical(dict(causal_gene_clade_id=key, implied_loss_species_clade_id=location)),
                )
            )
        if row["event_type"] == "leaf":
            if lost:
                raise ValueError("A gene leaf cannot imply a loss")
            source["genes"].append((row["gene_name"], family, row["species_name"], clade))
    if root_branch is None:
        raise ValueError("LCA gene root is missing")
    source["gene_root_species_branch_id"] = root_branch
    source["loss_semantics"] = "lca_implied_loss_below_gene_root"
    leaves = {str(n.name) for n in gene_tree.leaves()}
    if {g[0] for g in source["genes"]} != leaves:
        raise ValueError("LCA gene leaves disagree with the rooted tree")
    return source


def branch_rows(source, reference):
    taxa = {r[0]: set(json.loads(r[3])) for r in reference}
    counts = Counter((e[5], e[3]) for e in source["events"])
    genes = source["genes"]
    all_taxa = set().union(*taxa.values())
    if len({g[0] for g in genes}) != len(genes) or any(g[2] not in all_taxa for g in genes):
        raise ValueError("Invalid or duplicate reconciled genes")
    gene_counts = Counter(g[2] for g in genes)
    scope = taxa.get(source.get("gene_root_species_branch_id"))
    if source["status"] == "analysed" and source["method"] == "lca" and scope is None:
        raise ValueError("Missing LCA gene-root scope")
    transfer = source["method"] == "generax" and source["metadata"].get("reconciliation_model") == "UndatedDTL"
    if any(e[3] in {"branchingOut", "transferBack"} for e in source["events"]):
        transfer = True
    rows = []
    for key, descendants in taxa.items():
        status = source["status"]
        if status == "analysed" and source["method"] == "lca" and not descendants <= scope:
            status = "outside_gene_root_scope"
        numeric = status == "analysed"
        rows.append(
            dict(
                zip(
                    COLUMNS,
                    (
                        source["orthogroup"],
                        source["source_id"],
                        source["method"],
                        key,
                        status,
                        counts[key, "loss"] if numeric else None,
                        counts[key, "duplication"] if numeric else None,
                        counts[key, "speciation"] if numeric else None,
                        counts[key, "transferBack"] if numeric and transfer else None,
                        counts[key, "branchingOut"] if numeric and transfer else None,
                        sum(gene_counts[s] for s in descendants) if source["status"] == "analysed" else None,
                        sum(gene_counts[s] > 0 for s in descendants) if source["status"] == "analysed" else None,
                    ),
                    strict=True,
                )
            )
        )
    if sum(r["loss_count"] or 0 for r in rows) != sum(e[3] == "loss" for e in source["events"]):
        raise ValueError("Species-branch loss totals disagree with detailed events")
    return rows


def build_bundle(
    family,
    species_tree,
    gene_tree="",
    lca="",
    generax_xml="",
    generax_manifest="",
    generax_model="",
    generax_tree="",
    lca_manifest="",
):
    species_data = Path(species_tree).read_bytes()
    reference, topology = reference_branches(species_data)
    gene_data, lca_data, xml_data = map(input_bytes, (gene_tree, lca, generax_xml))
    generax_gene_data = input_bytes(generax_tree) if generax_tree else None
    sources = []
    for method in ("generax", "lca"):
        data = xml_data if method == "generax" else lca_data
        metadata = {"aggregation_contract": 1}
        if method == "generax":
            metadata["declared_reconciliation_model"] = generax_model or None
            metadata["reconciliation_model"] = None
            producer = input_bytes(generax_manifest)
            if producer is not None:
                record = json.loads(producer)
                if record.get("family_id") != family or record.get("step") != "generax":
                    raise ValueError("GeneRax producer identity mismatch")
                declared = [e for e in record.get("outputs", []) if e.get("label") == "generax_xml"]
                if data is not None and (
                    len(declared) != 1 or declared[0]["sha256"] != hashlib.sha256(data).hexdigest()
                ):
                    raise ValueError("GeneRax XML differs from its producer")
                metadata["producer"] = record
                actual_model = record.get("parameters", {}).get("reconciliation_model")
                if actual_model and generax_model and actual_model != generax_model:
                    raise ValueError("Declared GeneRax model differs from producer")
                metadata["reconciliation_model"] = actual_model
                if generax_gene_data is not None:
                    trees = [e for e in record.get("outputs", []) if e.get("label") == "generax_nhx"]
                    if len(trees) != 1 or trees[0]["sha256"] != hashlib.sha256(generax_gene_data).hexdigest():
                        raise ValueError("Reused GeneRax rooted tree differs from its producer")
        else:
            producer = input_bytes(lca_manifest)
            if producer is not None:
                record = json.loads(producer)
                if record.get("family_id") != family or record.get("step") != "reconciliation":
                    raise ValueError("LCA producer identity mismatch")
                declared = [e for e in record.get("outputs", []) if e.get("label") == "reconciliation"]
                if data is not None and (
                    len(declared) != 1 or declared[0]["sha256"] != hashlib.sha256(data).hexdigest()
                ):
                    raise ValueError("LCA table differs from its producer")
                for label, content in (("rooted_tree", gene_data), ("species_tree", species_data)):
                    original = [e for e in record.get("inputs", []) if e.get("label") == label]
                    if (
                        content is None
                        or len(original) != 1
                        or original[0].get("sha256") != hashlib.sha256(content).hexdigest()
                    ):
                        raise ValueError("LCA input differs from its producer: " + label)
                metadata["producer"] = record
        inputs = dict(
            species_tree=species_data,
            rooted_gene_tree=generax_gene_data if method == "generax" else gene_data,
            reconciliation=data,
        )
        source = source_record(family, method, inputs, metadata)
        if data is None:
            source["status"] = "not_available"
        elif method == "generax":
            source["clades"], source["events"], source["genes"] = parse_family(data, family, reference)
            source["loss_semantics"] = "reconstructed_lineage_extinction_not_localised_molecular_deletion"
        else:
            if gene_data is None:
                raise ValueError("LCA statistics require the exact rooted gene tree")
            lca_records(source, data, gene_data, species_data, reference)
        sources.append(source)
    return dict(
        schema=SCHEMA,
        orthogroup=family,
        species_tree_sha256=hashlib.sha256(species_data).hexdigest(),
        topology_sha256=topology,
        reference=reference,
        sources=sources,
    )


def validate_bundle(bundle, family, reference=None):
    if bundle["schema"] != SCHEMA or bundle["orthogroup"] != family:
        raise ValueError("Reconciliation bundle identity mismatch")
    if reference is not None and {r[0]: r[1:] for r in bundle["reference"]} != {r[0]: list(r[1:]) for r in reference}:
        # Node labels/sibling order may differ; rooted topology and descendant sets must not.
        actual = {(r[0], r[1], r[3], r[4]) for r in bundle["reference"]}
        expected = {(r[0], r[1], r[3], r[4]) for r in reference}
        if actual != expected:
            raise ValueError("Reconciliation bundle reference topology differs")
    seen = set()
    allowed = {r[0] for r in bundle["reference"]}
    for source in bundle["sources"]:
        identity = {k: source[k] for k in ("schema", "orthogroup", "method", "inputs", "metadata")}
        if (
            source["orthogroup"] != family
            or source["method"] not in {"lca", "generax"}
            or source["method"] in seen
            or digest(identity) != source["source_id"]
        ):
            raise ValueError("Reconciliation source contract mismatch")
        seen.add(source["method"])
        if source["status"] not in {"analysed", "not_available"}:
            raise ValueError("Unknown reconciliation source status")
        ids = {r[1] for r in source["clades"]}
        if len(ids) != len(source["clades"]) or any(
            r[0] != family or r[2] is not None and r[2] not in ids for r in source["clades"]
        ):
            raise ValueError("Reconciliation clade identity mismatch")
        event_keys = {(r[1], r[2]) for r in source["events"]}
        if len(event_keys) != len(source["events"]) or any(
            r[0] != family or r[1] not in ids or r[5] not in allowed for r in source["events"]
        ):
            raise ValueError("Reconciliation event identity mismatch")
        if any(g[1] != family or g[3] not in ids for g in source["genes"]):
            raise ValueError("Reconciled gene clade is missing")
        if source["status"] != "analysed" and any(source[k] for k in ("clades", "events", "genes")):
            raise ValueError("Unavailable reconciliation contains events")
        branch_rows(source, bundle["reference"])
    if seen != {"lca", "generax"}:
        raise ValueError("Bundle must record both methods, including unavailable status")
    return bundle


def create_tables(conn):
    conn.exec_driver_sql(
        "CREATE TABLE species_reconciliation_source (source_id TEXT PRIMARY KEY, orthogroup TEXT, method TEXT, status TEXT, metadata_json TEXT)"
    )
    conn.exec_driver_sql(
        "CREATE TABLE species_reconciliation_clade (source_id TEXT, orthogroup TEXT, clade_id INTEGER, parent_id INTEGER, label TEXT, PRIMARY KEY(source_id,clade_id))"
    )
    conn.exec_driver_sql(
        "CREATE TABLE species_reconciliation_event (source_id TEXT, orthogroup TEXT, clade_id INTEGER, event_order INTEGER, event_type TEXT, branch_id TEXT, attributes_json TEXT, PRIMARY KEY(source_id,clade_id,event_order))"
    )
    conn.exec_driver_sql(
        "CREATE TABLE species_reconciliation_gene (source_id TEXT, orthogroup TEXT, gene_id TEXT, species TEXT, clade_id INTEGER, PRIMARY KEY(source_id,gene_id))"
    )
    columns = ",".join('"' + c + '" ' + ("INTEGER" if c.endswith("_count") else "TEXT") for c in COLUMNS)
    conn.exec_driver_sql("CREATE TABLE stat_species_branch (" + columns + ",PRIMARY KEY(source_id,species_branch_id))")


def import_bundle(conn, bundle, stats_data):
    validate_bundle(bundle, bundle["orthogroup"])
    expected = [r for s in bundle["sources"] for r in branch_rows(s, bundle["reference"])]
    actual = list(csv.DictReader(io.StringIO(stats_data.decode()), delimiter="\t"))
    formatted = [{k: "" if v is None else str(v) for k, v in row.items()} for row in expected]
    if actual != formatted:
        raise ValueError("Species-branch TSV differs from its detailed reconciliation bundle")
    for source in bundle["sources"]:
        sid, family = source["source_id"], source["orthogroup"]
        meta = {k: v for k, v in source.items() if k not in {"clades", "events", "genes"}}
        conn.exec_driver_sql(
            "INSERT INTO species_reconciliation_source VALUES (?,?,?,?,?)",
            (sid, family, source["method"], source["status"], canonical(meta)),
        )
        inserts = [
            ("species_reconciliation_clade", [(sid, *r) for r in source["clades"]], 5),
            ("species_reconciliation_event", [(sid, r[0], r[1], r[2], r[3], r[5], r[6]) for r in source["events"]], 7),
            ("species_reconciliation_gene", [(sid, r[1], r[0], r[2], r[3]) for r in source["genes"]], 5),
        ]
        for table, rows, slots in inserts:
            for start in range(0, len(rows), 1000):
                conn.exec_driver_sql(
                    f"INSERT INTO {table} VALUES ({','.join('?' for _ in range(slots))})", rows[start : start + 1000]
                )
    conn.exec_driver_sql(
        "INSERT INTO stat_species_branch VALUES (" + ",".join("?" for _ in COLUMNS) + ")",
        [tuple(r[c] for c in COLUMNS) for r in expected],
    )


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--family-id", required=True)
    p.add_argument("--species-tree", required=True)
    for name in (
        "gene-tree",
        "lca",
        "lca-manifest",
        "generax-xml",
        "generax-manifest",
        "generax-model",
        "generax-tree",
    ):
        p.add_argument("--" + name, default="")
    p.add_argument("--bundle-out", required=True)
    p.add_argument("--stats-out", required=True)
    args = p.parse_args(argv)
    from nwkit.output_transaction import output_transaction
    from nwkit.util import validate_distinct_output_paths, validate_outputs_do_not_replace_inputs

    outputs = [("--bundle-out", args.bundle_out), ("--stats-out", args.stats_out)]
    validate_distinct_output_paths(outputs)
    validate_outputs_do_not_replace_inputs(
        [
            (n, getattr(args, n[2:].replace("-", "_")))
            for n in (
                "--species-tree",
                "--gene-tree",
                "--lca",
                "--lca-manifest",
                "--generax-xml",
                "--generax-manifest",
                "--generax-tree",
            )
        ],
        outputs,
    )
    bundle = build_bundle(
        args.family_id,
        args.species_tree,
        args.gene_tree,
        args.lca,
        args.generax_xml,
        args.generax_manifest,
        args.generax_model,
        args.generax_tree,
        args.lca_manifest,
    )
    validate_bundle(bundle, args.family_id)
    with output_transaction([args.bundle_out, args.stats_out], follow_symlinks=False) as stage:
        stage.write_text(
            args.bundle_out,
            lambda f: json.dump(bundle, f, sort_keys=True, ensure_ascii=False, allow_nan=False, separators=(",", ":")),
        )

        def write_stats(f):
            writer = csv.DictWriter(f, COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for source in bundle["sources"]:
                writer.writerows(branch_rows(source, bundle["reference"]))

        stage.write_text(args.stats_out, write_stats)


if __name__ == "__main__":
    main()
