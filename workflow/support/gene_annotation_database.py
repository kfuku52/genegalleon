"""Add source-bound per-gene annotations to a copy of a published family DB.

The original tree/event database is never modified or reconstructed. Native TSV
columns remain directly queryable, including localisation scores and their
model, calibration and decision metadata. GO tuples have a normalised table.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sqlite3
import sys
import tempfile
from collections import Counter
from contextlib import ExitStack
from functools import lru_cache
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from cds_annotation_scope import SCHEMA as CDS_SCOPE_SCHEMA
from cds_annotation_scope import ids_digest, validate_qc_row
from reconciliation_database import SCHEMA as RECONCILIATION_SCHEMA
from shared_namespace_lock import namespace_lock

SCHEMA = "gene_annotation_database_v1"
GO_COLUMNS = {"ids": "go_ids", "aspects": "go_aspects", "terms": "go_terms", "evidence": "go_evidence"}
RESERVED = {"annotation_source_id", "source_orthogroup"}
GO_ID = re.compile(r"GO:[0-9]{7}")


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)


def identity(path):
    info = Path(path).stat()
    return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns


def quoted(name):
    return '"' + name.replace('"', '""') + '"'


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def load_manifest(path):
    path = Path(path).resolve(strict=True)
    before = identity(path)
    config = json.loads(path.read_bytes())
    if identity(path) != before:
        raise ValueError("Annotation manifest changed while reading")
    allowed = {"schema", "gene_model_contract_sha256", "family_catalog_sha256", "sources",
               "source_database_sha256", "require_complete_genes"}
    if not isinstance(config, dict) or config.get("schema") != SCHEMA or set(config) - allowed:
        raise ValueError("Unsupported annotation manifest schema or keys")
    for key in ("gene_model_contract_sha256", "family_catalog_sha256", "source_database_sha256"):
        value = config.get(key, "")
        if (key != "source_database_sha256" or value) and (not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value)):
            raise ValueError("Invalid annotation contract: " + key)
    if not isinstance(config.get("require_complete_genes", True), bool):
        raise ValueError("require_complete_genes must be a JSON boolean")
    sources, species = [], set()
    if not isinstance(config.get("sources"), list) or not config["sources"]:
        raise ValueError("A nonempty annotation source list is required")
    for entry in config["sources"]:
        if not isinstance(entry, dict) or set(entry) - {"species", "path", "sha256", "product_column", "go_columns", "cds_scope"}:
            raise ValueError("Unsupported annotation source keys")
        name = entry.get("species")
        if not isinstance(name, str) or not name.strip() or name != name.strip() or name in species:
            raise ValueError("Empty or duplicate annotation species")
        species.add(name)
        if not re.fullmatch(r"[0-9a-f]{64}", entry.get("sha256", "")):
            raise ValueError("Each annotation source requires its SHA-256")
        target = (path.parent / entry["path"]).resolve(strict=True)
        if "\n" in str(target):
            raise ValueError("Newlines in annotation paths are unsupported")
        item = {**entry, "path": target}
        sources.append(item)
    return config, sources, {path: before}


def check_cds_scope(source, ids, accepted, statuses):
    scope = source.get("cds_scope")
    if scope is None:
        return
    if (not isinstance(scope, dict) or scope.get("schema") != CDS_SCOPE_SCHEMA
            or any(not re.fullmatch(r"[0-9a-f]{64}", scope.get(key, "")) for key in
                   ("catalogue_sha256", "resolved_cds_sha256", "resolution_sha256"))
            or scope.get("catalogue_gene_ids_sha256") != ids_digest(ids)
            or scope.get("accepted_gene_ids_sha256") != ids_digest(accepted)
            or scope.get("input_count") != len(ids) or scope.get("accepted_count") != len(accepted)
            or scope.get("excluded_count") != statuses["excluded"]
            or len(ids) != statuses["accepted"] + statuses["excluded"]):
        raise ValueError("Annotation catalogue/CDS eligibility scope mismatch: " + source["species"])


def preflight_annotations(source_db, config, sources, guards):
    """Fail scope/coverage errors before copying a potentially large family DB."""
    with sqlite3.connect(source_db.as_uri() + "?mode=ro", uri=True) as database:
        if database.execute("SELECT 1 FROM sqlite_master WHERE name='gene_annotation'").fetchone():
            raise ValueError("Use the original family DB, not an already annotated derivative")
        if not database.execute("SELECT 1 FROM sqlite_master WHERE name='reconciliation_meta'").fetchone():
            raise ValueError("Annotation import requires a model-bound reconciliation DB; enable gene_family_database_reconciliation")
        meta = dict(database.execute("SELECT key,value FROM reconciliation_meta"))
        if (meta.get("schema") != RECONCILIATION_SCHEMA
                or meta.get("gene_model_contract_sha256") != config["gene_model_contract_sha256"]
                or meta.get("family_catalog_sha256") != config["family_catalog_sha256"]):
            raise ValueError("Annotation/reconciliation model or family contract mismatch")
        expected = {(species, gene): family for gene, species, family in database.execute(
            "SELECT gene_id,species,orthogroup FROM reconciliation_gene")}
        expected_species = {json.loads(row[0])[0] for row in database.execute(
            "SELECT descendant_taxa FROM reconciliation_species_branch WHERE terminal=1")}
    if config.get("require_complete_genes", True) and {v["species"] for v in sources} != expected_species:
        raise ValueError("Complete annotation requires the exact DB species set")
    covered = set()
    for source in sources:
        guards[source["path"]] = identity(source["path"])
        header, product, go = source_header(source)
        ids, accepted, statuses = set(), set(), Counter()
        for row in source_rows(source, header):
            gene = row["gene_id"]
            if not gene or gene != gene.strip() or gene in ids or ("species" in row and row["species"] != source["species"]):
                raise ValueError("Empty/ambiguous/duplicate gene ID or conflicting annotation species")
            ids.add(gene)
            key = source["species"], gene
            if key in expected:
                if row["orthogroup"] != expected[key]:
                    raise ValueError("Annotation family disagrees with exact reconciliation membership: " + gene)
                covered.add(key)
            status = validate_qc_row(row, product, go)
            statuses[status] += 1
            if status == "accepted":
                accepted.add(gene)
        check_cds_scope(source, ids, accepted, statuses)
    missing = Counter(species for species, _ in expected.keys() - covered)
    if config.get("require_complete_genes", True) and missing:
        raise ValueError(f"Annotations omit {sum(missing.values())} exact reconciliation genes; by species: {canonical(dict(missing))}")


def source_header(source):
    with source["path"].open(encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"), None)
    if (not header or not {"gene_id", "orthogroup"} <= set(header)
            or len({name.casefold() for name in header}) != len(header)
            or any(not name or "\x00" in name for name in header)
            or RESERVED & {name.casefold() for name in header}
            or any(name.casefold() in {"gene_id", "species", "orthogroup"} and name != name.casefold() for name in header)):
        raise ValueError("Invalid or reserved annotation columns")
    product = source.get("product_column", "sprot_recname")
    if product not in header or product in {"species", "gene_id", "orthogroup"}:
        raise ValueError("Declared annotation product column is absent: " + product)
    go = source.get("go_columns", GO_COLUMNS if "go_ids" in header else {})
    if not isinstance(go, dict) or (go and (set(go) != set(GO_COLUMNS)
               or len(set(go.values())) != 4 or not set(go.values()) <= set(header))):
        raise ValueError("GO IDs, aspects, terms and evidence require four aligned columns")
    return header, product, go


def source_rows(source, header):
    """Hash exact input bytes during the CSV pass, including quoted newlines."""
    csv.field_size_limit(max(csv.field_size_limit(), 64 * 1024 * 1024))
    checksum = hashlib.sha256()
    with source["path"].open("rb") as handle:
        def lines():
            first = True
            for line in handle:
                checksum.update(line)
                yield line.decode("utf-8-sig" if first else "utf-8")
                first = False
        reader = csv.DictReader(lines(), delimiter="\t")
        if reader.fieldnames != header:
            raise ValueError("Annotation header changed before import")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Malformed annotation TSV row")
            yield row
    if checksum.hexdigest() != source["sha256"]:
        raise ValueError("Annotation source SHA mismatch: " + str(source["path"]))


def go_tuples(row, columns):
    if not columns:
        return []
    return list(_parse_go_values(tuple(row[columns[key]] for key in GO_COLUMNS)))


def _parse_go_values(values):
    if not any(values):
        return ()
    parts = [[item.strip() for item in value.split(";")] for value in values]
    if len({len(value) for value in parts}) != 1:
        raise ValueError("GO annotation tuple cardinality mismatch")
    result = []
    for index, (go_id, aspect, term, evidence) in enumerate(zip(*parts, strict=True)):
        if not GO_ID.fullmatch(go_id) or aspect not in {"P", "F", "C"}:
            raise ValueError("Invalid GO ID/aspect; no annotation tuples are guessed")
        result.append((index, go_id, aspect, term, evidence))
    return tuple(result)


def import_annotations(database, config, sources, guards):
    if not database.execute("SELECT 1 FROM sqlite_master WHERE name='reconciliation_meta'").fetchone():
        raise ValueError("Annotation import requires a model-bound reconciliation DB; enable gene_family_database_reconciliation")
    meta = dict(database.execute("SELECT key,value FROM reconciliation_meta"))
    if (meta.get("schema") != RECONCILIATION_SCHEMA
            or meta.get("gene_model_contract_sha256") != config["gene_model_contract_sha256"]
            or meta.get("family_catalog_sha256") != config["family_catalog_sha256"]):
        raise ValueError("Annotation/reconciliation model or family contract mismatch")
    if database.execute("SELECT 1 FROM sqlite_master WHERE name='gene_annotation'").fetchone():
        raise ValueError("Use the original family DB, not an already annotated derivative")
    families = {row[0] for row in database.execute("SELECT orthogroup FROM reconciliation_family")}
    expected = {(species, gene): family for gene, species, family in database.execute(
        "SELECT gene_id,species,orthogroup FROM reconciliation_gene")}
    expected_species = {json.loads(row[0])[0] for row in database.execute(
        "SELECT descendant_taxa FROM reconciliation_species_branch WHERE terminal=1")}
    if config.get("require_complete_genes", True) and {item["species"] for item in sources} != expected_species:
        raise ValueError("Complete annotation requires the exact DB species set")
    fields, headers = {}, []
    for source in sources:
        before = identity(source["path"])
        if guards.setdefault(source["path"], before) != before:
            raise ValueError("Annotation source changed after preflight")
        header, product, go = source_header(source)
        headers.append((header, product, go))
        for name in header:
            prior = fields.setdefault(name.casefold(), name)
            if prior != name:
                raise ValueError("Case-ambiguous annotation columns across sources")
    columns = sorted(name for name in fields.values() if name not in {"species", "gene_id", "orthogroup"})
    if len(columns) > 1000:
        raise ValueError("Annotation column count exceeds the bounded wide-table schema")
    database.execute("CREATE TABLE gene_annotation_source (annotation_source_id INTEGER PRIMARY KEY, species TEXT UNIQUE NOT NULL, path TEXT NOT NULL, sha256 TEXT NOT NULL, columns_json TEXT NOT NULL, product_column TEXT NOT NULL, go_columns_json TEXT NOT NULL, cds_scope_json TEXT NOT NULL)")
    database.execute("CREATE TABLE gene_annotation (species TEXT NOT NULL, gene_id TEXT NOT NULL, orthogroup TEXT, source_orthogroup TEXT NOT NULL, annotation_source_id INTEGER NOT NULL, " +
                     ",".join(quoted(name) + " TEXT" for name in columns) + ", PRIMARY KEY(species,gene_id)) WITHOUT ROWID")
    database.execute("CREATE TABLE gene_annotation_go (species TEXT NOT NULL, gene_id TEXT NOT NULL, ordinal INTEGER NOT NULL, go_id TEXT NOT NULL, aspect TEXT NOT NULL, term TEXT NOT NULL, evidence TEXT NOT NULL, annotation_source_id INTEGER NOT NULL, PRIMARY KEY(species,gene_id,ordinal)) WITHOUT ROWID")
    slots = ",".join("?" for _ in range(5 + len(columns)))
    # Only immutable parsed tuples are cached, within this import. Large fields
    # bypass the cache; both entry count and admitted key size remain bounded.
    parse_go = lru_cache(maxsize=256)(_parse_go_values)
    covered, total, terms, source_counts, qc_counts = set(), 0, 0, {}, Counter()
    for source_id, (source, (header, product, go)) in enumerate(zip(sources, headers, strict=True), 1):
        species = source["species"]
        database.execute("INSERT INTO gene_annotation_source VALUES (?,?,?,?,?,?,?,?)",
                         (source_id, species, str(source["path"]), source["sha256"], canonical(header), product, canonical(go),
                          canonical(source.get("cds_scope"))))
        batch, go_batch, count, batch_bytes = [], [], 0, 0
        source_ids, accepted_ids, source_statuses = set(), set(), Counter()
        for row in source_rows(source, header):
            gene, family = row["gene_id"], row["orthogroup"]
            if not gene or gene != gene.strip() or ("species" in row and row["species"] != species):
                raise ValueError("Empty/ambiguous gene ID or conflicting annotation species")
            key = species, gene
            source_ids.add(gene)
            status = validate_qc_row(row, product, go)
            source_statuses[status] += 1
            qc_counts[status] += 1
            if status == "accepted":
                accepted_ids.add(gene)
            if key in expected:
                if expected[key] != family:
                    raise ValueError("Annotation family disagrees with exact reconciliation membership: " + gene)
                covered.add(key)
            batch.append((species, gene, family if family in families else None, family, source_id,
                          *map(row.get, columns)))
            batch_bytes += sum(map(len, row.values()))
            values = tuple(row[go[key]] for key in GO_COLUMNS) if go else ()
            entries = (parse_go(values) if sum(map(len, values)) <= 4096 else _parse_go_values(values))
            go_batch.extend((species, gene, *entry, source_id) for entry in entries)
            terms += len(entries)
            count += 1
            if len(batch) >= 1000 or batch_bytes >= 8 * 1024 * 1024:
                database.executemany("INSERT INTO gene_annotation VALUES (" + slots + ")", batch)
                database.executemany("INSERT INTO gene_annotation_go VALUES (?,?,?,?,?,?,?,?)", go_batch)
                batch.clear()
                go_batch.clear()
                batch_bytes = 0
        if batch:
            database.executemany("INSERT INTO gene_annotation VALUES (" + slots + ")", batch)
            database.executemany("INSERT INTO gene_annotation_go VALUES (?,?,?,?,?,?,?,?)", go_batch)
        source_counts[species] = count
        check_cds_scope(source, source_ids, accepted_ids, source_statuses)
        total += count
    missing = len(expected.keys() - covered)
    if config.get("require_complete_genes", True) and missing:
        raise ValueError(f"Annotations omit {missing} exact reconciliation genes")
    database.execute("CREATE INDEX gene_annotation_family ON gene_annotation(orthogroup,species)")
    database.execute("CREATE INDEX gene_annotation_go_id ON gene_annotation_go(go_id,species,gene_id)")
    cases = " ".join("WHEN " + str(index) + " THEN " + quoted(product)
                     for index, (_, product, _) in enumerate(headers, 1))
    database.execute("CREATE VIEW gene_annotation_product AS SELECT species,gene_id,orthogroup,CASE annotation_source_id " + cases + " END AS product,annotation_source_id FROM gene_annotation")
    audit = {"genes": total, "GO_tuples": terms, "genes_by_species": source_counts,
             "reconciliation_genes": len(expected), "reconciliation_genes_without_annotation": missing,
             "cds_qc_status_counts": dict(qc_counts)}
    database.execute("CREATE TABLE gene_annotation_meta (key TEXT PRIMARY KEY,value TEXT NOT NULL)")
    database.executemany("INSERT INTO gene_annotation_meta VALUES (?,?)", [
        ("schema", SCHEMA), ("gene_model_contract_sha256", config["gene_model_contract_sha256"]),
        ("family_catalog_sha256", config["family_catalog_sha256"]), ("audit", canonical(audit)),
        ("localization_semantics", "Preserved prediction scores, labels and model/decision/calibration metadata; not experimental localisation."),
        ("annotation_semantics", "Extant gene annotations, not identified extinct-copy products; source_orthogroup preserves membership outside the selected catalogue.")])
    return audit


def build_database(source_db, manifest, output_db, scratch=None, record_argv=None):
    source, manifest, target = Path(source_db).resolve(strict=True), Path(manifest).resolve(strict=True), Path(output_db).absolute()
    config, sources, guards = load_manifest(manifest)
    inputs = {source, manifest, *(item["path"] for item in sources)}
    if target.resolve() in inputs or target.is_symlink() or (target.exists() and not target.is_file()):
        raise ValueError("Annotation output must be a separate database file")
    target.parent.mkdir(parents=True, exist_ok=True)
    with ExitStack() as stack:
        stack.enter_context(namespace_lock(source.with_name("." + source.name + ".build.guard"), exclusive=False))
        stack.enter_context(namespace_lock(target.with_name("." + target.name + ".build.guard"), exclusive=True))
        record = None
        if record_argv:
            import artifact_provenance as provenance
            record = provenance.build_parser().parse_args(["record", *record_argv])
            if ({Path(path).resolve() for _, path in provenance.parse_path_pairs(record.input, "--input")} != inputs
                    or {Path(path).absolute() for _, path in provenance.parse_path_pairs(record.output, "--output")} != {target}
                    or record.optional_output or record.input_gene_family_store or record.input_gene_family_subdir
                    or record.input_gene_family_artifact or record.input_logical_directory or record.output_logical_directory
                    or Path(record.manifest).resolve() in inputs | {target}):
                raise ValueError("Annotation provenance must bind exactly the source DB, manifest, TSVs and sole output DB")
            stack.enter_context(provenance.artifact_manifest_lock(record))
        if any(Path(str(source) + suffix).exists() for suffix in ("-wal", "-journal")):
            raise ValueError("Source database must be a closed published snapshot")
        guards[source] = identity(source)
        preflight_annotations(source, config, sources, guards)
        temp_root = Path(stack.enter_context(tempfile.TemporaryDirectory(prefix="gg-gene-annotations-", dir=scratch)))
        private = temp_root / "database.db"
        checksum = hashlib.sha256()
        with source.open("rb") as reader, private.open("xb") as writer:
            for chunk in iter(lambda: reader.read(1024 * 1024), b""):
                checksum.update(chunk)
                writer.write(chunk)
        source_sha = checksum.hexdigest()
        if config.get("source_database_sha256") and source_sha != config["source_database_sha256"]:
            raise ValueError("Source database SHA differs from the admitted snapshot")
        with sqlite3.connect(private) as database:
            database.execute("BEGIN")
            audit = import_annotations(database, config, sources, guards)
            database.executemany("INSERT INTO gene_annotation_meta VALUES (?,?)", [
                ("source_database_sha256", source_sha), ("manifest_sha256", sha256(manifest))])
            database.commit()
        fd, stage_name = tempfile.mkstemp(prefix="." + target.name + ".", dir=target.parent)
        staged = Path(stage_name)
        try:
            with os.fdopen(fd, "wb") as writer, private.open("rb") as reader:
                shutil.copyfileobj(reader, writer, 1024 * 1024)
                writer.flush()
                os.fsync(writer.fileno())
                os.fchmod(writer.fileno(), 0o400)
            payload = None
            if record is not None:
                payload = provenance.build_contract(record, include_diagnostics=True,
                    output_sources={label: str(staged) for label, _ in provenance.parse_path_pairs(record.output, "--output")})
            if (any(identity(path) != before for path, before in guards.items())
                    or any(Path(str(source) + suffix).exists() for suffix in ("-wal", "-journal"))):
                raise ValueError("An annotation input changed before publication")
            os.replace(staged, target)
            if payload is not None:
                provenance.write_manifest_atomic(record.manifest, payload)
                if os.environ.get("GG_OBSERVATION_ATTEMPT_DIR"):
                    from workflow_observation import record_contract_result
                    record_contract_result(record, 0)
        finally:
            if staged.exists():
                staged.unlink()
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-db")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-db")
    parser.add_argument("--scratch")
    parser.add_argument("--list-inputs", action="store_true")
    parser.add_argument("--record", action="append", help="Use --record=TOKEN for native provenance arguments")
    args = parser.parse_args()
    if args.list_inputs:
        _, sources, _ = load_manifest(args.manifest)
        for index, source in enumerate(sources):
            print("annotation_" + str(index) + "=" + str(source["path"]))
    else:
        if not args.source_db or not args.output_db:
            parser.error("--source-db and --output-db are required")
        print(canonical(build_database(args.source_db, args.manifest, args.output_db, args.scratch, args.record)))


if __name__ == "__main__":
    main()
