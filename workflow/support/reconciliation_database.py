"""Normalised GeneRax recPhyloXML evidence, separate from NWKIT LCA statistics.

Loss records describe extinction of a reconstructed lineage. They do not identify
an extant gene, a molecular deletion, or the precise time of that extinction.
"""
from __future__ import annotations

import hashlib
import io
import json
import re
from pathlib import Path
from xml.etree import ElementTree as ET

from Bio import Phylo

SCHEMA = "generax_recphyloxml_v1"
EVENTS = frozenset(("leaf", "speciation", "duplication", "loss", "branchingOut", "transferBack"))


def canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"))


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def tag(element):
    return element.tag.rsplit("}", 1)[-1]


def children(element, name):
    return [child for child in element if tag(child) == name]


def one(element, name):
    found = children(element, name)
    if len(found) != 1:
        raise ValueError(f"Expected one {name}, found {len(found)}")
    return found[0]


def text_name(element):
    name = one(element, "name").text
    if not name or not name.strip():
        raise ValueError("Empty node name")
    return name.strip()


def branch_id(taxa):
    return "clade_" + digest(sorted(taxa))


def reference_branches(data):
    tree = Phylo.read(io.StringIO(data.decode("utf-8")), "newick")
    rows = []
    seen_tips, seen_clades = set(), set()

    def visit(node, parent=None):
        taxa = sorted(tip.name for tip in node.get_terminals())
        if any(not name for name in taxa) or len(taxa) != len(set(taxa)):
            raise ValueError("Species tree requires unique nonempty tip names")
        key = branch_id(taxa)
        if key in seen_clades:
            raise ValueError("Species tree contains unary or repeated clades")
        seen_clades.add(key)
        if not node.clades:
            seen_tips.add(node.name)
        rows.append((key, parent, node.name, canonical(taxa), int(not node.clades)))
        for child in node.clades:
            visit(child, key)

    visit(tree.root)
    topology = digest(sorted(row[0] for row in rows))
    return rows, topology


def xml_species_map(root, reference):
    node = one(one(one(root, "spTree"), "phylogeny"), "clade")
    mapping, sets = {}, set()

    def visit(clade):
        label = text_name(clade)
        if label in mapping:
            raise ValueError("Duplicate XML species label: " + label)
        # Reserve before descending, including ancestor/descendant duplicates.
        mapping[label] = None
        descendants = children(clade, "clade")
        taxa = [label] if not descendants else [tip for child in descendants for tip in visit(child)]
        if len(taxa) != len(set(taxa)):
            raise ValueError("Duplicate XML species tip")
        key = branch_id(taxa)
        if key in sets:
            raise ValueError("Repeated XML species clade")
        sets.add(key)
        mapping[label] = key
        return taxa

    visit(node)
    if sets != {row[0] for row in reference}:
        raise ValueError("XML species topology/taxa differ from the reference tree")
    return mapping


def parse_family(data, family, reference):
    root = ET.fromstring(data)
    if tag(root) != "recPhylo":
        raise ValueError("Expected recPhyloXML")
    mapping = xml_species_map(root, reference)
    tips = {row[0]: json.loads(row[3])[0] for row in reference if row[4]}
    gene_root = one(one(one(root, "recGeneTree"), "phylogeny"), "clade")
    clades, events, genes, seen_genes = [], [], [], set()

    pending = [(gene_root, None)]
    while pending:
        node, parent = pending.pop()
        ordinal = len(clades)
        label = text_name(node)
        clades.append((family, ordinal, parent, label))
        descendants = children(node, "clade")
        node_events = list(one(node, "eventsRec"))
        if not node_events:
            raise ValueError("Gene clade has no reconciliation event")
        leaf_events = []
        for order, event in enumerate(node_events):
            kind = tag(event)
            if kind not in EVENTS:
                raise ValueError("Unsupported reconciliation event: " + kind)
            attr = "destinationSpecies" if kind == "transferBack" else "speciesLocation"
            location = event.get(attr)
            if location not in mapping:
                raise ValueError(f"Unmapped {kind} species: {location}")
            events.append((family, ordinal, order, kind, location, mapping[location], canonical(event.attrib)))
            if kind == "leaf":
                leaf_events.append(mapping[location])
        if leaf_events:
            if descendants or len(leaf_events) != 1 or leaf_events[0] not in tips or label in seen_genes:
                raise ValueError("Invalid or duplicate extant gene leaf: " + label)
            seen_genes.add(label)
            genes.append((label, family, tips[leaf_events[0]], ordinal))
        elif not descendants and not any(tag(event) == "loss" for event in node_events):
            raise ValueError("Terminal gene clade is neither a leaf nor a loss")
        pending.extend((child, ordinal) for child in reversed(descendants))
    return clades, events, genes


def write_tables(engine, store, species_tree, method_manifest=""):
    """Stream one family at a time into the builder's private DB transaction."""
    if store is None:
        raise ValueError("Reconciliation import requires --dir_gene_family")
    tree_data = Path(species_tree).read_bytes()
    reference, topology = reference_branches(tree_data)
    manifest_data = Path(method_manifest).read_bytes() if method_manifest else b"{}"
    manifest = json.loads(manifest_data)
    if not isinstance(manifest, dict) or not isinstance(manifest.get("default", {}), dict):
        raise ValueError("Method manifest must contain a default object and family overrides")
    model_contract = manifest.get("gene_model_contract_sha256", "")
    if not isinstance(model_contract, str) or (model_contract and not re.fullmatch(r"[0-9a-f]{64}", model_contract)):
        raise ValueError("Invalid gene-model contract SHA-256")
    overrides = manifest.get("families", {})
    if not isinstance(overrides, dict) or any(not isinstance(value, dict) for value in overrides.values()):
        raise ValueError("Method family overrides must be objects")
    bundle_names = set(store.file_names("reconciliation_events"))
    stats_names = set(store.file_names("stat_species_branch"))
    with engine.begin() as conn:
        if conn.exec_driver_sql("SELECT name FROM sqlite_master WHERE name='reconciliation_meta'").first():
            raise ValueError("Reconciliation append is unsupported; rebuild the complete audited database")
        families = [row[0] for row in conn.exec_driver_sql("SELECT DISTINCT orthogroup FROM tree ORDER BY orthogroup")]
        if not families or not all(isinstance(family, str) and family for family in families):
            raise ValueError("Reconciliation requires a nonempty family catalogue")
        if set(overrides) - set(families):
            raise ValueError("Method overrides include families outside the database catalogue")
        conn.exec_driver_sql("CREATE TABLE reconciliation_meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
        conn.exec_driver_sql("CREATE TABLE reconciliation_species_branch (branch_id TEXT PRIMARY KEY, parent_id TEXT, reference_label TEXT, descendant_taxa TEXT NOT NULL, terminal INTEGER NOT NULL)")
        conn.exec_driver_sql("CREATE TABLE reconciliation_family (orthogroup TEXT PRIMARY KEY, status TEXT NOT NULL, xml_sha256 TEXT, method_json TEXT NOT NULL, generax_manifest_sha256 TEXT, generax_manifest_json TEXT)")
        conn.exec_driver_sql("CREATE TABLE reconciliation_clade (orthogroup TEXT, clade_id INTEGER, parent_id INTEGER, label TEXT, PRIMARY KEY(orthogroup,clade_id))")
        conn.exec_driver_sql("CREATE TABLE reconciliation_event (orthogroup TEXT, clade_id INTEGER, event_order INTEGER, event_type TEXT, source_species_label TEXT, branch_id TEXT, attributes_json TEXT, PRIMARY KEY(orthogroup,clade_id,event_order))")
        conn.exec_driver_sql("CREATE TABLE reconciliation_gene (gene_id TEXT PRIMARY KEY, orthogroup TEXT, species TEXT, clade_id INTEGER)")
        conn.exec_driver_sql("INSERT INTO reconciliation_species_branch VALUES (?,?,?,?,?)", reference)
        meta = {"schema": SCHEMA, "topology_sha256": topology,
                "species_tree_sha256": hashlib.sha256(tree_data).hexdigest(),
                "family_catalog_sha256": digest(families),
                "method_manifest_sha256": hashlib.sha256(manifest_data).hexdigest(),
                "gene_model_contract_sha256": manifest.get("gene_model_contract_sha256", ""),
                "loss_semantics": "reconstructed_lineage_extinction_not_localised_molecular_deletion"}
        conn.exec_driver_sql("INSERT INTO reconciliation_meta VALUES (?,?)", list(meta.items()))
        if bundle_names:
            from species_branch_statistics import create_tables
            create_tables(conn)
        available = set(store.file_names("generax_xml"))
        provenance_names = set(store.file_names("artifact_provenance"))
        for family in families:
            name = family + "_generax.xml"
            method = {**manifest.get("default", {}), **overrides.get(family, {})}
            provenance_name = family + ".generax.json"
            source, source_sha = None, None
            bundle = None
            bundle_name, stats_name = family + "_reconciliation.json", family + "_stat.species_branch.tsv"
            if bundle_name in bundle_names:
                from species_branch_statistics import import_bundle, validate_bundle
                if stats_name not in stats_names:
                    raise ValueError("Reconciliation bundle lacks species-branch statistics: " + family)
                with store.open_binary("reconciliation_events", bundle_name) as handle:
                    bundle_data = handle.read()
                with store.open_binary("stat_species_branch", stats_name) as handle:
                    stats_data = handle.read()
                summary_name = family + ".summary_statistics.json"
                if summary_name not in provenance_names:
                    raise ValueError("Species-branch bundle lacks producer provenance: " + family)
                with store.open_binary("artifact_provenance", summary_name) as handle:
                    summary_producer = json.load(handle)
                if summary_producer.get("family_id") != family or summary_producer.get("step") != "summary_statistics":
                    raise ValueError("Species-branch producer identity mismatch: " + family)
                for label, content in (("reconciliation_events", bundle_data), ("stat_species_branch", stats_data)):
                    declared = [e for e in summary_producer.get("outputs", []) if e.get("label") == label]
                    if len(declared) != 1 or declared[0].get("sha256") != hashlib.sha256(content).hexdigest():
                        raise ValueError("Species-branch output differs from its producer: " + family)
                bundle = validate_bundle(json.loads(bundle_data), family, reference)
                import_bundle(conn, bundle, stats_data)
            elif stats_name in stats_names:
                raise ValueError("Species-branch statistics lack their detailed bundle: " + family)
            if provenance_name in provenance_names:
                with store.open_binary("artifact_provenance", provenance_name) as handle:
                    source_data = handle.read()
                source = json.loads(source_data)
                if source.get("family_id") != family or source.get("step") != "generax":
                    raise ValueError("GeneRax manifest family/step mismatch: " + family)
                source_sha = hashlib.sha256(source_data).hexdigest()
            if name not in available:
                conn.exec_driver_sql("INSERT INTO reconciliation_family VALUES (?,?,?,?,?,?)",
                                     (family, "missing_xml", None, canonical(method), source_sha,
                                      canonical(source) if source else None))
                continue
            with store.open_binary("generax_xml", name) as handle:
                data = handle.read()
            xml_sha = hashlib.sha256(data).hexdigest()
            if bundle is None:
                clades, events, genes = parse_family(data, family, reference)
            else:
                normal = next(s for s in bundle["sources"] if s["method"] == "generax")
                if normal["status"] != "analysed" or normal["inputs"].get("reconciliation") != xml_sha:
                    raise ValueError("Normalised GeneRax source differs from XML: " + family)
                clades, events, genes = ([tuple(row) for row in normal[k]] for k in ("clades", "events", "genes"))
            if source is not None:
                declared = [entry for entry in source.get("outputs", []) if entry.get("label") == "generax_xml"]
                if len(declared) != 1 or declared[0].get("sha256") != xml_sha:
                    raise ValueError("GeneRax XML differs from its producer manifest: " + family)
            conn.exec_driver_sql("INSERT INTO reconciliation_family VALUES (?,?,?,?,?,?)",
                                 (family, "imported", xml_sha, canonical(method), source_sha,
                                  canonical(source) if source else None))
            for table, rows, slots in (("reconciliation_clade", clades, 4),
                                       ("reconciliation_event", events, 7), ("reconciliation_gene", genes, 4)):
                if rows:
                    conn.exec_driver_sql(f"INSERT INTO {table} VALUES ({','.join('?' for _ in range(slots))})", rows)
        conn.exec_driver_sql("CREATE INDEX reconciliation_event_branch ON reconciliation_event(branch_id,orthogroup,event_type)")
        conn.exec_driver_sql("CREATE INDEX reconciliation_gene_family ON reconciliation_gene(orthogroup,species)")
