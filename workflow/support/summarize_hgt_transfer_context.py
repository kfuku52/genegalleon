#!/usr/bin/env python3
"""Event-resolved, bilateral scaffold context from existing reconciliations.

No origin is inferred from best hits, no support is imputed, and no threshold
is applied. The reconciliation species tree is authoritative; present-day
descendants of an ancestral branch are explicitly identified as proxies.
"""

import argparse
import csv
import hashlib
import io
import math
import os
import re
import tempfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

from gene_family_output_store import GeneFamilyOutputStore, read_only_observation
from hgt_species_tree import read_species_tree
from scaffold_taxonomy import CONTEXT_COLUMNS, GENE_COLUMNS, METRICS, RANKS, species_key

IDENTITY_COLUMNS = ["event_id", "orthogroup", "branch_id", "node_name", "event_index",
                    "generax_transfer", "generax_donor_node", "generax_recipient_node"]
SIDE_COLUMNS = ["branch_type", "evidence_basis", "descendant_species", "context_status",
                "all_descendant_gene_count", "retained_gene_count", "excluded_gene_count",
                "mapped_gene_count", "scaffold_count", "retained_genes"]
EVENT_COLUMNS = [*IDENTITY_COLUMNS, "mapping_status", "mapping_reason", "xml_event_id",
                 "species_tree_mapping_status", "generax_xml_sha256", "stat_branch_sha256",
                 "reconciliation_species_tree_sha256", "support_generax_ufboot",
                 "support_status", "support_source",
                 *[f"{side}_{c}" for side in ("donor", "recipient")
                   for c in (*SIDE_COLUMNS, *CONTEXT_COLUMNS)]]
AUX_COLUMNS = ["intron_supported", "expression_measured", "synteny_support_score",
               "contamination_lca_taxid", "contamination_lca_sciname",
               "contamination_is_compatible_lineage"]
LINK_COLUMNS = [*IDENTITY_COLUMNS, "side", "gene_id", "gene_species", "lineage_status",
                "eligible_for_context", "context_reason", *GENE_COLUMNS, *AUX_COLUMNS]


def read_tsv(path):
    csv.field_size_limit(100_000_000)
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, columns):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows({c: row.get(c, "") for c in columns} for row in rows)


def labels(value):
    return frozenset(x.strip() for x in str(value).split(";") if x.strip())


def transfers(value):
    # Each token is an independent event, including malformed/empty tokens.
    return re.split(r"\s*[;,|]\s*", str(value).strip()) if str(value).strip() else [""]


def species_nodes(clade):
    nodes = {}
    for node in reversed(list(clade.iter("clade"))):
        children = node.findall("clade")
        name = species_key(node.findtext("name"))
        leaves = set().union(*(nodes[id(c)][1] for c in children)) if children else {name}
        if not children and not name:
            raise ValueError("Unnamed reconciliation species-tree terminal")
        nodes[id(node)] = (name, leaves, "internal" if children else "terminal")
    result = {}
    for name, leaves, kind in nodes.values():
        if name:
            if name in result:
                raise ValueError(f"Duplicate reconciliation species-tree name: {name}")
            result[name] = (frozenset(leaves), kind)
    return result


def parse_reconciliation(text):
    root = ET.fromstring(text)
    for element in root.iter():
        element.tag = element.tag.split("}")[-1]
    tree = root.find("recGeneTree/phylogeny/clade")
    if tree is None:
        raise ValueError("Missing recGeneTree in GeneRax XML")
    species = root.find("spTree/phylogeny/clade")
    nodes = species_nodes(species) if species is not None else {}
    species_sha = hashlib.sha256(ET.tostring(species)).hexdigest() if species is not None else ""
    leafsets, gene_species = {}, {}
    for clade in reversed(list(tree.iter("clade"))):
        children = clade.findall("clade")
        leaves = set().union(*(leafsets[c] for c in children)) if children else set()
        events = clade.find("eventsRec")
        for event in events if events is not None else []:
            if event.tag == "leaf":
                gene = clade.findtext("name")
                if not gene or gene in gene_species:
                    raise ValueError(f"Missing/duplicate reconciliation gene: {gene}")
                gene_species[gene] = species_key(event.get("speciesLocation"))
                leaves.add(gene)
        leafsets[clade] = frozenset(leaves)
    matches = defaultdict(list)
    for index, clade in enumerate(tree.iter("clade"), 1):
        events = clade.find("eventsRec")
        donors = [species_key(e.get("speciesLocation")) for e in events
                  if e.tag == "branchingOut"] if events is not None else []
        children = clade.findall("clade")
        for child_index, child in enumerate(children, 1):
            events = child.find("eventsRec")
            destinations = [species_key(e.get("destinationSpecies") or e.get("speciesLocation"))
                            for e in events if e.tag == "transferBack"] if events is not None else []
            for donor in donors:
                for recipient in destinations:
                    matches[(donor, recipient, leafsets[clade])].append({
                        "xml_event_id": f"{index}:{child_index}", "recipient": [child],
                        "donor": [other for other in children if other is not child],
                        "ambiguous": len(donors) != 1 or len(destinations) != 1,
                    })
    return nodes, gene_species, leafsets, matches, species_sha


def retained_leaves(clade, leafsets, allow_initial_transfer=False):
    leaves, pending = set(), [(clade, allow_initial_transfer)]
    while pending:
        node, allowed = pending.pop()
        events = node.find("eventsRec")
        if not allowed and events is not None and any(e.tag == "transferBack" for e in events):
            continue
        children = node.findall("clade")
        if children:
            pending.extend((child, False) for child in children)
        else:
            leaves.update(leafsets[node])
    return frozenset(leaves)


def external_nodes(path):
    if not path:
        return {}
    tree = read_species_tree(path)
    return {species_key(n.name): frozenset(species_key(t.name) for t in n.get_terminals())
            for n in tree.find_clades() if n.name}


def numeric(value):
    if value is None or str(value).strip().lower() in {"", "na", "nan", "none"}:
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Nonfinite scaffold measurement: {value}")
    return result


def validated_context(row):
    result = {c: row.get(c, "") for c in GENE_COLUMNS}
    for background in ("", "background_"):
        for rank in RANKS:
            prefix = f"host_scaffold_{background}{rank}_"
            values = [numeric(row.get(prefix + m)) for m in METRICS]
            if all(v is None for v in values):
                continue
            if any(v is None for v in values[:5]):
                raise ValueError("Partial scaffold count measurements")
            n, c, i, u, cds = values[:5]
            if any(v < 0 or not v.is_integer() for v in values[:5]) or c + i + u != n or cds > n:
                raise ValueError("Invalid scaffold count measurements")
            expected = [(c+i)/n if n else None, c/(c+i) if c+i else None, c/n if n else None]
            for observed, calculated in zip(values[5:], expected, strict=True):
                if (observed is None) != (calculated is None) or (
                        observed is not None and abs(observed - calculated) > 1e-9):
                    raise ValueError("Inconsistent scaffold fraction")
            result.update({prefix + m: int(v) if j < 5 else (v if v is not None else "")
                           for j, (m, v) in enumerate(zip(METRICS, values, strict=True))})
    return result


def summarize_side(event, side, match, nodes, gene_species, leafsets, genes):
    roots = match[side]
    all_genes = frozenset().union(*(leafsets[c] for c in roots))
    retained = frozenset().union(*(retained_leaves(c, leafsets, side == "recipient") for c in roots))
    node = event[f"generax_{side}_node"]
    allowed, kind = nodes[node]
    links, scaffolds, eligible = [], {}, []
    for gene in sorted(all_genes):
        species = gene_species[gene]
        lineage = ("transferred_out" if gene not in retained else
                   "species_unresolved" if not species else
                   "outside_species_branch" if species not in allowed else "retained")
        source = genes.get((event["orthogroup"], gene))
        row = {c: event[c] for c in IDENTITY_COLUMNS}
        row.update(side=side, gene_id=gene, gene_species=species, lineage_status=lineage,
                   eligible_for_context=lineage == "retained", context_reason="",
                   host_scaffold_status="gene_not_in_summary")
        if source is not None:
            if species_key(source.get("gene_taxon")) != species:
                row["host_scaffold_status"] = "gene_species_mismatch"
            else:
                row.update(validated_context(source))
            row.update({c: source.get(c, "") for c in AUX_COLUMNS})
        if lineage == "retained":
            eligible.append(gene)
            if row["host_scaffold_status"] == "measured" and row.get("host_scaffold_id"):
                key = (species, row["host_scaffold_id"])
                metrics = {c: row.get(c, "") for c in CONTEXT_COLUMNS}
                if key in scaffolds and scaffolds[key] != metrics:
                    raise ValueError(f"Conflicting scaffold context: {key}")
                scaffolds[key] = metrics
            else:
                row["context_reason"] = row["host_scaffold_status"]
        else:
            row["context_reason"] = lineage
        links.append(row)
    mapped = sum(r["eligible_for_context"] and r["host_scaffold_status"] == "measured"
                 and bool(r.get("host_scaffold_id")) for r in links)
    status = ("no_retained_extant_gene" if not eligible else "no_mapped_scaffold" if not scaffolds else
              "measured" if mapped == len(eligible) else "partial")
    summary = dict(branch_type=kind, evidence_basis="extant_descendant_proxy" if kind == "internal"
                   else "extant_terminal_genome", descendant_species="; ".join(sorted(allowed)),
                   context_status=status, all_descendant_gene_count=len(all_genes),
                   retained_gene_count=len(eligible), excluded_gene_count=len(all_genes)-len(eligible),
                   mapped_gene_count=mapped, scaffold_count=len(scaffolds), retained_genes="; ".join(eligible))
    for background in ("", "background_"):
        for rank in RANKS:
            prefix = f"host_scaffold_{background}{rank}_"
            values = [[numeric(s.get(prefix + m)) for m in METRICS[:5]] for s in scaffolds.values()]
            if not values or any(any(v is None for v in row) for row in values):
                continue
            totals = [int(sum(row[j] for row in values)) for j in range(5)]
            n, c, i, _, _ = totals
            totals.extend([(c+i)/n if n else "", c/(c+i) if c+i else "", c/n if n else ""])
            summary.update({prefix + m: v for m, v in zip(METRICS, totals, strict=True)})
    event.update({f"{side}_{c}": v for c, v in summary.items()})
    return links


def summarize(branches, gene_rows, store, species_tree=""):
    genes = {}
    for row in gene_rows:
        key = (row["orthogroup"], row["gene_id"])
        if key in genes:
            raise ValueError(f"Duplicate gene-summary key: {key}")
        genes[key] = row
    grouped = defaultdict(list)
    seen = set()
    for row in branches:
        key = (row["orthogroup"], row["branch_id"])
        if key in seen:
            raise ValueError(f"Duplicate branch-summary key: {key}")
        seen.add(key)
        grouped[row["orthogroup"]].append(row)
    expected_nodes = external_nodes(species_tree)
    events, links = [], []
    for family, rows in grouped.items():
        xml_text = stat_text = None
        try:
            with store.open_binary("generax_xml", family + "_generax.xml") as handle:
                xml_text = handle.read().decode("utf-8")
        except FileNotFoundError:
            pass
        try:
            with store.open_binary("stat_branch", family + "_stat.branch.tsv") as handle:
                stat_text = handle.read().decode("utf-8")
        except FileNotFoundError:
            pass
        parsed = parse_reconciliation(xml_text) if xml_text is not None else None
        stats = {}
        if stat_text is not None:
            for stat in csv.DictReader(io.StringIO(stat_text), delimiter="\t"):
                if stat["branch_id"] in stats:
                    raise ValueError(f"Duplicate stat_branch key in {family}")
                stats[stat["branch_id"]] = stat
        for branch in rows:
            branch_genes = labels(branch.get("candidate_genes", branch.get("gene_labels", "")))
            stat = stats.get(branch["branch_id"])
            for index, transfer in enumerate(transfers(branch.get("generax_transfer", "")), 1):
                event = {"event_id": f"{family}:{branch['branch_id']}:{index}", "orthogroup": family,
                         "branch_id": branch["branch_id"], "node_name": branch.get("node_name", ""),
                         "event_index": index, "generax_transfer": transfer,
                         "generax_donor_node": "", "generax_recipient_node": "",
                         "mapping_status": "unresolved", "mapping_reason": "",
                         "support_status": "missing_stat_branch", "support_source": "",
                         "generax_xml_sha256": hashlib.sha256(xml_text.encode()).hexdigest() if xml_text else "",
                         "stat_branch_sha256": hashlib.sha256(stat_text.encode()).hexdigest() if stat_text else ""}
                events.append(event)
                stat_matches = stat is not None and labels(stat.get("gene_labels", "")) == branch_genes
                if stat is not None:
                    event["support_status"] = "stat_branch_mismatch" if not stat_matches else "missing_ufboot"
                    if stat_matches:
                        support = numeric(stat.get("support_generax_ufboot"))
                        terminal = numeric(stat.get("num_leaf")) == 1 or stat.get("so_event") == "L"
                        event["support_status"] = "terminal_branch" if terminal else "missing_ufboot"
                        if not terminal and support is not None:
                            if not 0 <= support <= 100:
                                raise ValueError("UFBoot outside 0..100")
                            event.update(support_generax_ufboot=support, support_status="measured",
                                         support_source="stat_branch.support_generax_ufboot")
                parts = transfer.split("@")
                if len(parts) != 3 or parts[0] != "Y" or not all(parts[1:]):
                    event["mapping_reason"] = "invalid_transfer_annotation"
                    continue
                donor, recipient = map(species_key, parts[1:])
                event.update(generax_donor_node=donor, generax_recipient_node=recipient)
                if parsed is None:
                    event["mapping_reason"] = "missing_generax_xml"
                    continue
                nodes, gene_species, leafsets, matches, species_sha = parsed
                event["reconciliation_species_tree_sha256"] = species_sha
                if donor not in nodes or recipient not in nodes:
                    event["mapping_reason"] = "species_branch_unresolved"
                    continue
                compatible = not expected_nodes or all(expected_nodes.get(n) == leaves for n, (leaves, _) in nodes.items())
                event["species_tree_mapping_status"] = ("matched_external_tree" if expected_nodes and compatible else
                                                        "xml_species_tree_only" if compatible else "external_tree_mismatch")
                if not compatible:
                    event["mapping_reason"] = "external_tree_mismatch"
                    continue
                matched = matches.get((donor, recipient, branch_genes), [])
                if len(matched) != 1 or matched[0].get("ambiguous", False):
                    event["mapping_reason"] = "xml_event_missing" if not matched else "xml_event_ambiguous"
                    continue
                event.update(mapping_status="matched", mapping_reason="exact_transfer_and_branch_descendants",
                             xml_event_id=matched[0]["xml_event_id"])
                for side in ("donor", "recipient"):
                    links.extend(summarize_side(event, side, matched[0], nodes, gene_species, leafsets, genes))
    return events, links


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--branch_tsv", required=True)
    parser.add_argument("--gene_tsv", required=True)
    parser.add_argument("--dir_gene_family", required=True)
    parser.add_argument("--species_tree", default="")
    parser.add_argument("--event_out", required=True)
    parser.add_argument("--event_gene_out", required=True)
    args = parser.parse_args()
    store = GeneFamilyOutputStore(args.dir_gene_family)
    branches = read_tsv(args.branch_tsv)
    wanted = {(b["orthogroup"], g) for b in branches
              for g in labels(b.get("candidate_genes", b.get("gene_labels", "")))}
    csv.field_size_limit(100_000_000)
    with open(args.gene_tsv, newline="", encoding="utf-8") as handle:
        gene_rows = [{c: row.get(c, "") for c in ("orthogroup", "gene_id", "gene_taxon", *GENE_COLUMNS, *AUX_COLUMNS)}
                     for row in csv.DictReader(handle, delimiter="\t")
                     if (row["orthogroup"], row["gene_id"]) in wanted]
    branch_groups, gene_groups = defaultdict(list), defaultdict(list)
    for row in branches:
        branch_groups[row["orthogroup"]].append(row)
    for row in gene_rows:
        gene_groups[row["orthogroup"]].append(row)
    outputs = [(args.event_out, EVENT_COLUMNS), (args.event_gene_out, LINK_COLUMNS)]
    temporary = []
    try:
        for path, columns in outputs:
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            handle = tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                                  dir=Path(path).parent, prefix=".hgt_context_", delete=False)
            temporary.append(handle)
        writers = [csv.DictWriter(h, columns, delimiter="\t", extrasaction="raise")
                   for h, (_, columns) in zip(temporary, outputs, strict=True)]
        for writer in writers:
            writer.writeheader()
        # Bound event-gene expansion to one family, including on large projects.
        with read_only_observation(), store.read_snapshot():
            for family, rows in branch_groups.items():
                events, links = summarize(rows, gene_groups[family], store, args.species_tree)
                for writer, records, (_, columns) in zip(writers, (events, links), outputs, strict=True):
                    writer.writerows({c: row.get(c, "") for c in columns} for row in records)
        for handle, (path, _) in zip(temporary, outputs, strict=True):
            handle.close()
            os.replace(handle.name, path)
    finally:
        for handle in temporary:
            handle.close()
            Path(handle.name).unlink(missing_ok=True)


if __name__ == "__main__":
    main()
