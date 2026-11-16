#!/usr/bin/env python3
"""Read-only DB downstream: lineage events and independently qualified evidence.

One result database and three tabular/JSON exports, irrespective of family count.
All external inputs are SHA-bound in a configuration; relative paths resolve
against that configuration, never the caller's working directory.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import sqlite3
import tempfile
from collections import Counter, defaultdict
from contextlib import ExitStack
from pathlib import Path

from gene_annotation_database import SCHEMA as ANNOTATION_SCHEMA
from reconciliation_database import SCHEMA, branch_id, canonical
from shared_namespace_lock import namespace_lock


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for data in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(data)
    return result.hexdigest()


def identity(path):
    stat = Path(path).stat()
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


def load_config(path):
    path = Path(path).resolve(strict=True)
    guards = {path: identity(path)}
    config = json.loads(path.read_bytes())
    if identity(path) != guards[path]:
        raise ValueError("Configuration changed while reading")
    if config.get("schema") != "reconciliation_summary_v1":
        raise ValueError("Unsupported reconciliation summary configuration")
    unknown = set(config) - {"schema", "inputs", "family_catalog_sha256", "gene_model_contract_sha256",
                             "expression_p_value_semantics", "expression_fdr", "expression_min_logfc",
                             "proteome_species", "require_complete_reconciliation", "event_source"}
    if unknown:
        raise ValueError("Unknown summary configuration keys: " + ",".join(sorted(unknown)))
    if not isinstance(config.get("require_complete_reconciliation", False), bool):
        raise ValueError("require_complete_reconciliation must be a JSON boolean")
    if config.get("event_source", "generax") not in {"generax", "lca"}:
        raise ValueError("event_source must be generax or lca")
    inputs = {}
    for key, entry in config.get("inputs", {}).items():
        if key not in {"branches", "cafe", "expression_transcript", "expression_cds",
                       "proteome", "shared_peptides", "proteome_members", "gene_annotations"}:
            raise ValueError("Unknown evidence input: " + key)
        target = (path.parent / entry["path"]).resolve(strict=True)
        before = identity(target)
        if sha256(target) != entry["sha256"] or identity(target) != before:
            raise ValueError("Evidence SHA mismatch or changed input: " + key)
        if "\n" in str(target):
            raise ValueError("Newlines in input paths are unsupported")
        inputs[key] = target
        guards[target] = before
    if "branches" not in inputs:
        raise ValueError("A SHA-bound branch target JSON is required")
    if ("proteome" in inputs or "shared_peptides" in inputs) and "proteome_members" not in inputs:
        raise ValueError("Proteome evidence requires the complete searched-species gene membership")
    return config, inputs, guards


def read_tsv(path, required):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames) or not set(required) <= set(reader.fieldnames):
            raise ValueError(f"Invalid evidence header: {path}")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Malformed evidence row: " + str(path))
            yield row


def finite(value, name):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("Non-finite " + name)
    return number


def bh(values):
    """BH among tested candidates of ONE branch, including valid P=1."""
    ordered = sorted(values.items(), key=lambda pair: pair[1])
    result, previous = {}, 1.0
    for rank in range(len(ordered), 0, -1):
        key, value = ordered[rank - 1]
        if not 0 <= value <= 1:
            raise ValueError("P value outside [0,1]")
        previous = min(previous, value * len(ordered) / rank)
        result[key] = previous
    return result


def read_expression(path, families):
    result = {}
    for row in read_tsv(path, ("orthogroup", "status", "logFC", "p_value")):
        family = row["orthogroup"]
        if family not in families or family in result:
            raise ValueError("Expression catalogue mismatch or duplicate family: " + family)
        if row["status"] == "success":
            pvalue, logfc = finite(row["p_value"], "p_value"), finite(row["logFC"], "logFC")
            if not 0 <= pvalue <= 1:
                raise ValueError("P value outside [0,1]")
            result[family] = ("success", logfc, pvalue)
        else:
            result[family] = (row["status"], None, None)
    if set(result) != families:
        raise ValueError("Expression input must include all families, including untested rows")
    return result


def collect_evidence(inputs, config, db, families, gene_query=None):
    gene_query = gene_query or "SELECT gene_id,orthogroup,species FROM reconciliation_gene"
    gene_family = {}
    if "proteome" in inputs or "gene_annotations" in inputs:
        gene_family = dict(db.execute("SELECT gene_id,orthogroup FROM ("+gene_query+")"))
    expression = {view: read_expression(inputs["expression_" + view], families)
                  for view in ("transcript", "cds") if "expression_" + view in inputs}
    if expression and config.get("expression_p_value_semantics") != "one_sided_trap_enrichment_positive_logfc_v1":
        raise ValueError("Expression requires explicit one-sided trap-enrichment P-value semantics")
    members, proteins, shared, products = defaultdict(set), defaultdict(set), Counter(), defaultdict(set)
    audit = {"identified_proteome_genes_without_selected_family": 0}
    database_products = {}
    if db.execute("SELECT 1 FROM sqlite_master WHERE name='gene_annotation_meta'").fetchone():
        annotation_meta = dict(db.execute("SELECT key,value FROM gene_annotation_meta"))
        meta = dict(db.execute("SELECT key,value FROM reconciliation_meta"))
        if (annotation_meta.get("schema") != ANNOTATION_SCHEMA
                or annotation_meta.get("gene_model_contract_sha256") != meta["gene_model_contract_sha256"]
                or annotation_meta.get("family_catalog_sha256") != meta["family_catalog_sha256"]):
            raise ValueError("Embedded gene annotations disagree with the reconciliation contract")
        for gene, family, product in db.execute(
                "SELECT r.gene_id,r.orthogroup,a.product FROM ("+gene_query+") r "
                "JOIN gene_annotation_product a ON a.species=r.species AND a.gene_id=r.gene_id"):
            database_products[gene] = product or ""
            if product:
                products[family].add(product)
        audit["gene_annotation_database"] = annotation_meta
        audit["reconciliation_genes_with_embedded_annotations"] = len(database_products)
    if "proteome_members" in inputs:
        species = config.get("proteome_species")
        if not species:
            raise ValueError("Explicit proteome_species is required")
        expected = {gene: family for gene, family in db.execute(
            "SELECT gene_id,orthogroup FROM ("+gene_query+") WHERE species=?", (species,))}
        found = {}
        for row in read_tsv(inputs["proteome_members"], ("gene_id", "orthogroup")):
            gene, family = row["gene_id"], row["orthogroup"]
            if gene in found or expected.get(gene) != family:
                raise ValueError("Proteome gene membership differs from the reconciliation DB: " + gene)
            found[gene] = family
            members[family].add(gene)
        if found != expected:
            raise ValueError("Incomplete searched-species gene membership")
        if "proteome" in inputs:
            seen_proteins = set()
            for row in read_tsv(inputs["proteome"], ("gene_id", "orthogroup", "unique_unmodified_peptides_IL_equivalent", "best_peptide_q", "best_protein_q")):
                gene, family = row["gene_id"], row["orthogroup"]
                if not gene or gene in seen_proteins:
                    raise ValueError("Invalid/duplicate proteome gene: " + gene)
                seen_proteins.add(gene)
                peptide_count = int(row["unique_unmodified_peptides_IL_equivalent"])
                peptide_q, protein_q = finite(row["best_peptide_q"], "peptide_q"), finite(row["best_protein_q"], "protein_q")
                if peptide_count < 1 or not 0 <= peptide_q <= 0.01 or not 0 <= protein_q <= 0.01:
                    raise ValueError("Proteome row does not meet the admitted unique-peptide/q<=0.01 contract")
                if not family and gene not in gene_family:
                    audit["identified_proteome_genes_without_selected_family"] += 1
                    continue
                if found.get(gene) != family:
                    raise ValueError("Proteome gene membership differs from the reconciliation DB: " + gene)
                proteins[family].add(gene)
        if "shared_peptides" in inputs:
            seen = set()
            for row in read_tsv(inputs["shared_peptides"], ("orthogroup", "peptide", "MS_file", "spectrum", "proteins", "peptide_q", "protein_q")):
                family = row["orthogroup"]
                key = (row["MS_file"], row["spectrum"])
                if family not in members or key in seen:
                    raise ValueError("Invalid/duplicate shared family PSM")
                ids = row["proteins"].split(";")
                if len(ids) < 2 or len(ids) != len(set(ids)) or any(found.get(gene) != family for gene in ids):
                    raise ValueError("Shared peptide does not map exclusively to multiple searched-species genes in one family")
                seen.add(key)
                if any(not 0 <= finite(row[name], name) <= 0.01 for name in ("peptide_q", "protein_q")):
                    raise ValueError("Shared peptide q value exceeds the admitted contract")
                shared[family] += 1
    if "gene_annotations" in inputs:
        seen = set()
        for row in read_tsv(inputs["gene_annotations"], ("gene_id", "product")):
            gene = row["gene_id"]
            if gene not in gene_family or gene in seen:
                raise ValueError("Annotation gene differs from the DB or is duplicated: " + gene)
            seen.add(gene)
            if gene in database_products and row["product"] != database_products[gene]:
                raise ValueError("External product conflicts with the embedded gene annotation: " + gene)
            if row["product"]:
                products[gene_family[gene]].add(row["product"])
    return expression, members, proteins, shared, products, audit


def build_results(database, config_path, output_dir, record_argv=None):
    target = Path(output_dir).absolute()
    source = Path(database).resolve(strict=True)
    target.parent.mkdir(parents=True, exist_ok=True)
    with ExitStack() as stack:
        # The normal DB builder uses this namespace exclusively. Keep it stable
        # through both result and provenance publication, without locking other DBs.
        stack.enter_context(namespace_lock(source.with_name("." + source.name + ".build.guard"), exclusive=False))
        stack.enter_context(namespace_lock(target.with_name("." + target.name + ".build.guard"), exclusive=True))
        record = None
        if record_argv:
            import artifact_provenance as provenance
            record = provenance.build_parser().parse_args(["record", *record_argv])
            stack.enter_context(provenance.artifact_manifest_lock(record))
        return _build_results(database, config_path, output_dir, record)


def _build_results(database, config_path, output_dir, record=None):
    database, config_path, output_dir = Path(database).resolve(strict=True), Path(config_path).resolve(strict=True), Path(output_dir).absolute()
    config, inputs, before = load_config(config_path)
    if output_dir.is_symlink() or (output_dir.exists() and not output_dir.is_dir()):
        raise ValueError("Output must be a directory, not a symlink or file")
    if any(Path(str(database) + suffix).exists() for suffix in ("-wal", "-journal")):
        raise ValueError("Input database must be a closed, published snapshot")
    if output_dir == database.parent or output_dir in database.parents or output_dir in config_path.parents:
        raise ValueError("Output directory must be separate from the source database/configuration")
    for path in inputs.values():
        if output_dir == path or output_dir in path.parents:
            raise ValueError("Output directory contains an evidence input")
    before[database] = identity(database)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix="." + output_dir.name + ".", dir=output_dir.parent))
    try:
        with sqlite3.connect(database.as_uri() + "?mode=ro", uri=True) as db, sqlite3.connect(temporary / "reconciliation_summary.db") as result:
            db.execute("PRAGMA query_only=ON")
            db.execute("BEGIN")
            meta = dict(db.execute("SELECT key,value FROM reconciliation_meta"))
            if meta.get("schema") != SCHEMA or config.get("family_catalog_sha256") != meta["family_catalog_sha256"]:
                raise ValueError("Reconciliation schema or family catalogue contract mismatch")
            has_evidence = set(inputs) - {"branches", "cafe"}
            has_annotations = db.execute("SELECT 1 FROM sqlite_master WHERE name='gene_annotation_meta'").fetchone()
            if (has_evidence or has_annotations) and (not meta.get("gene_model_contract_sha256") or
                                 config.get("gene_model_contract_sha256") != meta["gene_model_contract_sha256"]):
                raise ValueError("Functional evidence requires the identical gene-model contract")
            families = dict(db.execute("SELECT orthogroup,status FROM reconciliation_family"))
            method = config.get("event_source", "generax")
            common = db.execute("SELECT 1 FROM sqlite_master WHERE name='stat_species_branch'").fetchone() is not None
            gene_query = "SELECT gene_id,orthogroup,species FROM reconciliation_gene"
            if common:
                source_rows = list(db.execute("SELECT orthogroup,status FROM species_reconciliation_source WHERE method=?", (method,)))
                if len(source_rows) != len(dict(source_rows)):
                    raise ValueError("Ambiguous reconciliation sources for selected method")
                source_status = dict(source_rows)
                families = {f: "imported" if source_status.get(f) == "analysed" else source_status.get(f, "missing_statistics") for f in families}
                gene_query = ("SELECT g.gene_id,g.orthogroup,g.species FROM species_reconciliation_gene g "
                              "JOIN species_reconciliation_source s ON g.source_id=s.source_id WHERE s.method='"+method+"'")
            elif method != "generax":
                raise ValueError("LCA branch statistics are missing; regenerate gg_gene_evolution statistics")
            if config.get("require_complete_reconciliation", False) and any(status != "imported" for status in families.values()):
                raise ValueError("Complete reconciliation was required, but selected-method results are missing")
            gene_counts = defaultdict(Counter)
            for family, species, count in db.execute("SELECT orthogroup,species,COUNT(*) FROM ("+gene_query+") GROUP BY orthogroup,species"):
                gene_counts[family][species] = count
            branches = json.loads(inputs["branches"].read_bytes())
            if not isinstance(branches, list) or not branches:
                raise ValueError("Target branches must be a nonempty list")
            db_branches = {key: json.loads(taxa) for key, taxa in db.execute("SELECT branch_id,descendant_taxa FROM reconciliation_species_branch")}
            target_keys, labels = {}, set()
            for branch in branches:
                label, taxa = branch["branch_id"], branch["descendants"]
                if not isinstance(label, str) or not label or not taxa or len(taxa) != len(set(taxa)):
                    raise ValueError("Invalid target branch")
                key = branch_id(taxa)
                if key not in db_branches or key in target_keys or label in labels:
                    raise ValueError("Unknown or duplicate target branch: " + label)
                labels.add(label)
                target_keys[key] = label
            cafe = {}
            if "cafe" in inputs:
                for row in read_tsv(inputs["cafe"], ("branch_id", "orthogroup", "copy_change", "significant_decrease")):
                    key = (row["branch_id"], row["orthogroup"])
                    status = row.get("status", "analysed")
                    if key[0] not in labels or key[1] not in families or key in cafe or status not in ("analysed", "not_analysed"):
                        raise ValueError("Invalid CAFE branch/family evidence")
                    if status == "not_analysed":
                        if row["copy_change"] != "" or row["significant_decrease"] != "":
                            raise ValueError("Unanalysed CAFE evidence must have blank change/significance")
                        cafe[key] = (None, None, status)
                        continue
                    if row["significant_decrease"] not in ("0", "1"):
                        raise ValueError("Invalid CAFE significance")
                    change = int(row["copy_change"])
                    significant = int(row["significant_decrease"])
                    if significant and change >= 0:
                        raise ValueError("CAFE significant decrease must have a negative copy change")
                    cafe[key] = (change, significant, status)
                if len(cafe) != len(labels) * len(families):
                    raise ValueError("CAFE comparison must include every target branch/family, not only significant hits")
            expression, members, proteins, shared, products, evidence_audit = collect_evidence(inputs, config, db, set(families), gene_query)
            evidence_audit["cafe_analysed_branch_family_pairs"] = sum(row[2] == "analysed" for row in cafe.values())
            evidence_audit["cafe_explicitly_unanalysed_branch_family_pairs"] = sum(row[2] == "not_analysed" for row in cafe.values())
            fdr = finite(config.get("expression_fdr", 0.05), "FDR")
            logfc = finite(config.get("expression_min_logfc", 1.0), "logFC cutoff")
            if not 0 < fdr <= 1:
                raise ValueError("Expression FDR must be in (0,1]")
            columns = {"branch_id": "TEXT", "orthogroup": "TEXT", "reconciliation_status": "TEXT",
                       "lineage_extinction_events": "INTEGER", "duplication_events": "INTEGER",
                       "extant_genes_inside": "INTEGER", "extant_genes_outside": "INTEGER",
                       "cafe_status": "TEXT", "cafe_copy_change": "INTEGER", "cafe_significant_decrease": "INTEGER", "cafe_generax_comparison": "TEXT",
                       "proteome_status": "TEXT", "unique_proteome_genes": "TEXT", "shared_family_psms": "INTEGER", "products": "TEXT"}
            for view in ("transcript", "cds"):
                columns.update({view + "_" + field: dtype for field, dtype in
                                (("status", "TEXT"), ("logFC", "REAL"), ("p_value", "REAL"), ("q_branch_loss", "REAL"), ("trap_enriched", "INTEGER"))})
            columns.update(loss_events="INTEGER", reconciliation_method="TEXT", cafe_reconciliation_comparison="TEXT")
            result.execute("CREATE TABLE branch_family (" + ",".join('"' + key + '" ' + dtype for key, dtype in columns.items()) + ", PRIMARY KEY(branch_id,orthogroup))")
            # These values are independent of the target branch. Serialize once
            # per family, without caching branch-local event counts or BH results.
            family_products = {family: canonical(sorted(products[family])) for family in families}
            family_proteins = {family: canonical(sorted(proteins[family])) for family in families}
            family_totals = {family: sum(gene_counts[family].values()) for family in families}
            insert = "INSERT INTO branch_family VALUES (" + ",".join("?" for _ in columns) + ")"
            with (temporary / "branch_family.tsv").open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(columns)
                batch = []
                for key, label in target_keys.items():
                    branch_taxa = db_branches[key]
                    branch_members = set(branch_taxa)
                    unique_taxa = len(branch_members) == len(branch_taxa)
                    events = defaultdict(Counter)
                    branch_status = {}
                    if common:
                        for family, status, losses, duplications in db.execute(
                                "SELECT orthogroup,status,loss_count,duplication_count FROM stat_species_branch WHERE method=? AND species_branch_id=?", (method, key)):
                            branch_status[family] = "imported" if status == "analysed" else status
                            if status == "analysed":
                                events[family]["loss"], events[family]["duplication"] = losses, duplications
                    else:
                        for family, kind, count in db.execute("SELECT orthogroup,event_type,COUNT(*) FROM reconciliation_event WHERE branch_id=? GROUP BY orthogroup,event_type", (key,)):
                            events[family][kind] = count
                    candidates = {family for family in families if events[family]["loss"] > 0}
                    qvalues = {view: bh({family: data[family][2] for family in candidates if data[family][0] == "success"})
                               for view, data in expression.items()}
                    for family, status in families.items():
                        if common:
                            status = branch_status.get(family, "missing_statistics")
                        valid = status == "imported"
                        counts = gene_counts[family]
                        if unique_taxa and len(counts) < len(branch_taxa):
                            inside = sum(count for taxon, count in counts.items() if taxon in branch_members)
                        else:
                            inside = sum(counts[taxon] for taxon in branch_taxa)
                        cafe_change, cafe_significant, cafe_status = cafe.get((label, family), (None, None, "not_analysed"))
                        row = dict.fromkeys(columns)
                        row.update(branch_id=label, orthogroup=family, reconciliation_status=status,
                                   lineage_extinction_events=events[family]["loss"] if valid and method == "generax" else None,
                                   duplication_events=events[family]["duplication"] if valid else None,
                                   extant_genes_inside=inside if valid else None,
                                   extant_genes_outside=family_totals[family] - inside if valid else None,
                                   cafe_status=cafe_status, cafe_copy_change=cafe_change,
                                   cafe_significant_decrease=cafe_significant,
                                   cafe_generax_comparison=("not_analysed" if not valid or cafe_status != "analysed" else
                                                           "both" if family in candidates and cafe_significant else
                                                           "generax_extinction_only" if family in candidates else
                                                           "cafe_decrease_only" if cafe_significant else "neither"),
                                   proteome_status=("not_analysed" if "proteome" not in inputs else
                                                    "not_assessable_no_gene" if not members[family] else
                                                    "unique_gene_support" if proteins[family] else
                                                    "shared_family_support" if shared[family] else "not_detected"),
                                   unique_proteome_genes=family_proteins[family],
                                   shared_family_psms=shared[family] if "shared_peptides" in inputs else None,
                                   products=family_products[family])
                        for view in ("transcript", "cds"):
                            fit_status, effect, pvalue = expression.get(view, {}).get(family, ("not_analysed", None, None))
                            qvalue = qvalues.get(view, {}).get(family)
                            row.update({view + "_status": fit_status, view + "_logFC": effect, view + "_p_value": pvalue,
                                        view + "_q_branch_loss": qvalue,
                                        view + "_trap_enriched": int(qvalue < fdr and effect > logfc) if qvalue is not None else None})
                        comparison = row["cafe_generax_comparison"]
                        row.update(loss_events=events[family]["loss"] if valid else None,
                                   reconciliation_method=method,
                                   cafe_reconciliation_comparison="reconciliation_loss_only" if comparison == "generax_extinction_only" else comparison)
                        if method != "generax":
                            row["cafe_generax_comparison"] = None
                        values = tuple(row.values())
                        writer.writerow(values)
                        batch.append(values)
                        if len(batch) >= 1000:
                            result.executemany(insert, batch)
                            batch.clear()
                if batch:
                    result.executemany(insert, batch)
            result.execute("CREATE TABLE branch_summary AS SELECT branch_id, COUNT(*) AS families, SUM(reconciliation_status='imported') AS imported_families, SUM(lineage_extinction_events>0) AS families_with_lineage_extinction, SUM(lineage_extinction_events) AS lineage_extinction_events, SUM(extant_genes_inside>0) AS families_with_extant_descendants, SUM(cafe_status='analysed') AS cafe_analysed_families, SUM(cafe_status='not_analysed') AS cafe_unanalysed_families, SUM(cafe_generax_comparison!='not_analysed') AS cafe_generax_comparable_families, SUM(cafe_significant_decrease) AS cafe_significant_decreases, SUM(cafe_generax_comparison='both') AS generax_extinction_and_cafe_decrease, SUM(cafe_generax_comparison='generax_extinction_only') AS generax_extinction_only, SUM(cafe_generax_comparison='cafe_decrease_only') AS cafe_decrease_only, SUM(transcript_trap_enriched) AS transcript_trap_enriched_loss_families, SUM(cds_trap_enriched) AS cds_trap_enriched_loss_families, SUM(loss_events>0) AS families_with_loss, SUM(loss_events) AS loss_events, SUM(cafe_reconciliation_comparison!='not_analysed') AS reconciliation_comparable_families, SUM(cafe_reconciliation_comparison='both') AS loss_and_cafe_decrease FROM branch_family GROUP BY branch_id")
            with (temporary / "branch_summary.tsv").open("w", newline="") as handle:
                cursor = result.execute("SELECT * FROM branch_summary ORDER BY branch_id")
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow([col[0] for col in cursor.description])
                writer.writerows(cursor)
            manifest = {"schema": "reconciliation_summary_v1", "database_sha256": sha256(database),
                        "config_sha256": sha256(config_path), "database_contract": meta,
                        "evidence_audit": evidence_audit,
                        "inputs": {key: {"path": str(path), "sha256": config["inputs"][key]["sha256"]} for key, path in inputs.items()},
                        "expression_correction": "BH_within_each_branch_among_successfully_tested_lineage_extinction_families",
                        "event_source": method,
                        "loss_semantics": "lca_implied_loss_below_gene_root" if method == "lca" else meta["loss_semantics"],
                        "interpretation": ["Lineage extinction is not a localised molecular gene deletion.",
                                           "Extant descendant presence is not an ancestral copy-number estimate.",
                                           "Expression supports an orthogroup, not a particular extinct copy.",
                                           "Unanalysed CAFE families are not nonsignificant or unchanged; comparisons use analysed families only.",
                                           "Shared peptides support a family only; nondetection is not biological absence."],
                        "outputs": {name: sha256(temporary / name) for name in ("branch_family.tsv", "branch_summary.tsv")}}
            if method == "lca":
                manifest["expression_correction"] = "BH_within_each_branch_among_successfully_tested_lca_loss_families"
                manifest["interpretation"][0] = "LCA implied losses are conditional on the rooted gene/species trees, not localised molecular deletions."
            result.execute("CREATE TABLE provenance (json TEXT NOT NULL)")
            result.execute("INSERT INTO provenance VALUES (?)", (canonical(manifest),))
            result.commit()
            (temporary / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        payload = None
        if record is not None:
            import artifact_provenance as provenance
            recorded_inputs = provenance.parse_path_pairs(record.input, "--input")
            recorded_outputs = provenance.parse_path_pairs(record.output, "--output")
            expected_inputs = {database, config_path, *inputs.values()}
            expected_outputs = {output_dir / name for name in ("reconciliation_summary.db", "branch_family.tsv", "branch_summary.tsv", "manifest.json")}
            if ({Path(path).resolve() for _, path in recorded_inputs} != expected_inputs or
                    {Path(path).absolute() for _, path in recorded_outputs} != expected_outputs or
                    record.optional_output or record.input_gene_family_store or record.input_gene_family_subdir or
                    record.input_gene_family_artifact or record.input_logical_directory or record.output_logical_directory or
                    Path(record.manifest).absolute() in expected_inputs or output_dir in Path(record.manifest).absolute().parents):
                raise ValueError("Downstream provenance must bind the exact input/output file set")
            payload = provenance.build_contract(record, include_diagnostics=True,
                output_sources={label: str(temporary / Path(path).name) for label, path in recorded_outputs})
        if (any(identity(path) != token for path, token in before.items()) or
                any(Path(str(database) + suffix).exists() for suffix in ("-wal", "-journal"))):
            raise ValueError("An input changed before downstream publication")
        backup = temporary.with_name(temporary.name + ".previous")
        if output_dir.exists():
            os.rename(output_dir, backup)
        try:
            os.rename(temporary, output_dir)
        except BaseException:
            if backup.exists():
                os.rename(backup, output_dir)
            raise
        if backup.exists():
            shutil.rmtree(backup)
        if payload is not None:
            provenance.write_manifest_atomic(record.manifest, payload)
            if os.environ.get("GG_OBSERVATION_ATTEMPT_DIR"):
                from workflow_observation import record_contract_result
                record_contract_result(record, 0)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dbpath")
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir")
    parser.add_argument("--list-inputs", action="store_true")
    parser.add_argument("--record", action="append", help="Optional exact provenance record argv token; use --record=TOKEN")
    args = parser.parse_args()
    if args.list_inputs:
        _, inputs, _ = load_config(args.config)
        for key, path in inputs.items():
            print(key + "=" + str(path))
    else:
        if not args.dbpath or not args.output_dir:
            parser.error("--dbpath and --output-dir are required")
        build_results(args.dbpath, args.config, args.output_dir, args.record)


if __name__ == "__main__":
    main()
