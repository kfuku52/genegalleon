"""Resolve ANNEVO profiles from explicit, genome-bound NCBI taxonomy evidence."""
from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
import sys
from pathlib import Path

from annevo_runtime_download import atomic_json, canonical_bytes, file_sha

# ANNEVO uses NCBI clades, not similarly named traditional classifications.
PROFILES = {"Mammalia": 40674, "Aves": 8782, "Actinopteri": 186623,
            "Insecta": 50557, "Fungi": 4751, "Magnoliopsida": 3398}
OVERLAP_PROFILES = {"Mammalia", "Actinopteri"}


def stable_signature(path):
    info = path.stat()
    return [info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns]


def require_quiescent_taxonomy(path):
    if any(Path(str(path) + suffix).exists() for suffix in ("-wal", "-journal")):
        raise ValueError("Taxonomy database has an active transaction; use a quiescent snapshot")


def taxonomy_sha(path):
    path = Path(path)
    require_quiescent_taxonomy(path)
    before = stable_signature(path)
    digest = file_sha(path)
    require_quiescent_taxonomy(path)
    if stable_signature(path) != before:
        raise ValueError("Taxonomy database changed while hashing")
    return digest


def protect_outputs(outputs, inputs):
    protected = {Path(path).expanduser().resolve() for path in inputs if path}
    if any(Path(path).expanduser().resolve() in protected for path in outputs):
        raise ValueError("ANNEVO outputs must not overwrite inputs")


def write_unchanged(path, value):
    path = Path(path)
    if path.is_symlink():
        raise ValueError("Redirected ANNEVO selection receipt")
    if not path.exists() or path.read_bytes() != canonical_bytes(value):
        atomic_json(path, value)


def evidence_rows(path, species, genome_sha, *, metadata=False):
    path = Path(path).expanduser().absolute()
    before = stable_signature(path)
    digest = file_sha(path)
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        required = {"species", "taxid"} | ({"genome_sha256"} if metadata else set())
        if not required <= set(fields) or len(set(fields)) != len(fields):
            raise ValueError("TaxID evidence needs unique columns: " + ", ".join(sorted(required)))
        rows = []
        seen = set()
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError("Malformed TaxID evidence row")
            row = {k: v.strip() for k, v in row.items()}
            if not row["species"] or row["species"] in seen:
                raise ValueError("Empty or duplicate species in TaxID evidence: " + row["species"])
            seen.add(row["species"])
            if not re.fullmatch(r"[1-9][0-9]*", row["taxid"]):
                raise ValueError("Invalid TaxID in evidence for " + row["species"])
            if row["species"] == species:
                if metadata and row["genome_sha256"] != genome_sha:
                    raise ValueError("Genome metadata checksum differs from the selected genome")
                rows.append(row)
    if stable_signature(path) != before:
        raise ValueError("TaxID evidence changed while reading")
    return rows, {"path": str(path), "sha256": digest,
                  "kind": "genome_metadata" if metadata else "explicit_taxid_map"}


class Taxonomy:
    def __init__(self, path):
        # Read the exact main-file snapshot whose bytes are hashed, never a WAL
        # that appeared after the initial sidecar check. External main-file or
        # sidecar changes still invalidate selection at the closing fence.
        self.conn = sqlite3.connect(Path(path).resolve().as_uri() + "?mode=ro&immutable=1", uri=True)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute("BEGIN")
        self.tables = {r[0] for r in self.conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}

    def node(self, taxid):
        row = self.conn.execute("SELECT taxid,spname,rank,track FROM species WHERE taxid=?", (taxid,)).fetchone()
        return dict(row) if row else None

    def resolve(self, taxid):
        visited = []
        while self.node(taxid) is None:
            if taxid in visited or "merged" not in self.tables:
                raise ValueError("Unresolved or deleted NCBI TaxID: " + str(taxid))
            visited.append(taxid)
            merged = self.conn.execute("SELECT taxid_new FROM merged WHERE taxid_old=?", (taxid,)).fetchone()
            if merged is None:
                raise ValueError("Unresolved or deleted NCBI TaxID: " + str(taxid))
            taxid = int(merged[0])
        return taxid, visited

    def lineage(self, taxid):
        terminal = self.node(taxid)
        ids = [int(v) for v in terminal["track"].split(",") if v]
        if not ids or ids[0] != taxid or ids[-1] != 1 or len(set(ids)) != len(ids):
            raise ValueError("Invalid NCBI ancestor track")
        nodes = [self.node(v) for v in ids]
        for index, node in enumerate(nodes):
            if node is None or [int(v) for v in node["track"].split(",") if v] != ids[index:]:
                raise ValueError("Incomplete or contradictory NCBI ancestor tracks")
        # A clade-level TaxID alone does not identify the genome's species.
        if not any(node["rank"] == "species" for node in nodes):
            raise ValueError("TaxID evidence does not identify a species or an infraspecific taxon")
        return list(reversed(nodes))

    def check_name(self, name, lineage):
        if not name:
            return
        hits = {r[0] for r in self.conn.execute("SELECT taxid FROM species WHERE spname=? COLLATE NOCASE", (name,))}
        if "synonym" in self.tables:
            hits |= {r[0] for r in self.conn.execute("SELECT taxid FROM synonym WHERE spname=? COLLATE NOCASE", (name,))}
        if len(hits) != 1 or not hits <= {n["taxid"] for n in lineage}:
            raise ValueError("Ambiguous or contradictory scientific_name in genome metadata: " + name)


def select(args):
    # Protect both publication paths, including the error receipt, before any
    # exception handler can write to a caller-supplied input path.
    protect_outputs([args.output, str(args.output) + ".failed.json"],
                    [args.genome, args.taxid_map, args.genome_metadata, args.taxonomy_db, __file__])
    genome = Path(args.genome).absolute()
    record = {"schema_version": 1, "kind": "genegalleon-annevo-model-selection",
              "species": args.species, "requested_profile": args.model_profile,
              "genome_path": str(genome), "genome_sha256": file_sha(genome),
              "selector_sha256": file_sha(Path(__file__)), "status": "unresolved", "evidence": []}
    try:
        if args.model_profile != "auto":
            if args.model_profile not in PROFILES:
                raise ValueError("Unsupported explicit ANNEVO model profile")
            record.update(model_profile=args.model_profile, reason="explicit_model_profile", taxid=None,
                          lineage=[], taxonomy=None)
        else:
            candidates = []
            for path, metadata in ((args.taxid_map, False), (args.genome_metadata, True)):
                if path:
                    rows, evidence = evidence_rows(path, args.species, record["genome_sha256"], metadata=metadata)
                    record["evidence"].append(evidence)
                    candidates += rows
            if not candidates:
                raise ValueError("TaxID unresolved: supply an explicit species/TaxID map or hash-bound genome metadata; species labels are not guessed")
            db = Path(args.taxonomy_db).expanduser().absolute()
            if not db.is_file() or not db.stat().st_size:
                raise ValueError("NCBI taxonomy cache is missing; prepare it online or provide annevo_taxonomy_db")
            before = stable_signature(db)
            record["taxonomy"] = {"path": str(db), "sha256": taxonomy_sha(db), "format": "ETE NCBI SQLite"}
            taxonomy = Taxonomy(db)
            try:
                resolved = [taxonomy.resolve(int(row["taxid"])) for row in candidates]
                ids = {taxid for taxid, _ in resolved}
                record["input_taxids"] = [row["taxid"] for row in candidates]
                if len(ids) != 1:
                    raise ValueError("Conflicting NCBI TaxIDs in genome metadata and explicit mapping")
                taxid = ids.pop()
                lineage = taxonomy.lineage(taxid)
                record.update(taxid=taxid, merged_taxids=sorted({v for _, old in resolved for v in old}),
                              scientific_name=lineage[-1]["spname"], lineage=lineage)
                for row in candidates:
                    taxonomy.check_name(row.get("scientific_name", ""), lineage)
                matches = [profile for profile, anchor in PROFILES.items() if anchor in {n["taxid"] for n in lineage}]
                if len(matches) != 1:
                    raise ValueError("No supported ANNEVO model in NCBI ancestors" if not matches else "Conflicting ANNEVO clades in NCBI ancestors")
                profile = matches[0]
                anchor = next(n for n in lineage if n["taxid"] == PROFILES[profile])
                if anchor["spname"] != profile:
                    raise ValueError("NCBI model anchor name contradicts its TaxID")
                record.update(model_profile=profile, anchor_taxid=PROFILES[profile], reason="ncbi_ancestor_taxid")
            finally:
                taxonomy.conn.close()
            require_quiescent_taxonomy(db)
            if stable_signature(db) != before:
                raise ValueError("Taxonomy database changed during model selection")
        record.update(status="selected", overlap_pred=record["model_profile"] in OVERLAP_PROFILES)
        write_unchanged(args.output, record)
        return record
    except (OSError, ValueError, sqlite3.Error) as exc:
        record.update(status="error", reason=str(exc))
        write_unchanged(str(args.output) + ".failed.json", record)
        raise ValueError(str(exc)) from exc


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    for name in ("species", "genome", "output"):
        result.add_argument("--" + name, required=True)
    result.add_argument("--model-profile", default="auto")
    for name in ("taxonomy-db", "taxid-map", "genome-metadata"):
        result.add_argument("--" + name, default="")
    return result


if __name__ == "__main__":
    try:
        print(json.dumps(select(parser().parse_args()), sort_keys=True))
    except (OSError, ValueError, RuntimeError) as exc:
        print("ANNEVO model selection: " + str(exc), file=sys.stderr)
        raise SystemExit(1) from None
