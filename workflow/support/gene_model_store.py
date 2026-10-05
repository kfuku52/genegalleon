"""Streaming locus index and bounded-component conserved isoform selection.

Only locus identifiers and frozen correspondence edges are retained globally.
CDS/protein strings are loaded from SQLite for one connected component at a
time.  The store is a disposable derived artifact; source catalogs stay intact.
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import sqlite3
import tempfile
from collections import defaultdict
from pathlib import Path

if __package__:
    from .gene_model_selection import locus_can_vote, prepare_correspondence, select_representatives
else:
    from gene_model_selection import locus_can_vote, prepare_correspondence, select_representatives

SCHEMA_VERSION = 1


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _identity(path):
    stat = path.stat()
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


def build_store(catalog_dirs, sqlite_path):
    """Atomically import ``catalog_metadata.json`` and streaming ``loci.jsonl``.

    Reject source changes during import, duplicate species/loci/transcripts, and
    incompatible metadata.  Source-receipt and tool-identity cache policy belongs
    to the enclosing refinement stage, not this storage helper.
    """
    destination = Path(sqlite_path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=".locus-store-", suffix=".sqlite", dir=destination.parent)
    os.close(descriptor)
    temporary = Path(temporary_name)
    database = sqlite3.connect(temporary)
    sources, species_seen, locus_count, candidate_count, largest = [], set(), 0, 0, 0
    try:
        database.executescript("""
            PRAGMA journal_mode=DELETE;
            PRAGMA foreign_keys=ON;
            -- Keep bulky JSON outside the primary-key btree. WITHOUT ROWID
            -- makes foreign-key probes read overflowing record payloads even
            -- when they only need a species/locus key. Rowid tables retain
            -- separate compact unique indices and the same checked relations.
            CREATE TABLE catalogs (species TEXT NOT NULL PRIMARY KEY, json TEXT NOT NULL);
            CREATE TABLE loci (
                species TEXT NOT NULL, gene_id TEXT NOT NULL, seqid TEXT NOT NULL, strand TEXT NOT NULL,
                json TEXT NOT NULL, PRIMARY KEY(species, gene_id),
                FOREIGN KEY(species) REFERENCES catalogs(species)
            );
            CREATE TABLE candidate_owners (
                species TEXT NOT NULL, candidate_id TEXT NOT NULL, gene_id TEXT NOT NULL,
                PRIMARY KEY(species, candidate_id), FOREIGN KEY(species,gene_id) REFERENCES loci(species,gene_id)
            ) WITHOUT ROWID;
        """)
        with database:
            for directory in sorted({Path(path).expanduser().resolve() for path in catalog_dirs}):
                metadata_path, loci_path = directory / "catalog_metadata.json", directory / "loci.jsonl"
                before = {path: _identity(path) for path in (metadata_path, loci_path)}
                metadata_bytes = metadata_path.read_bytes()
                metadata = json.loads(metadata_bytes)
                if metadata.get("schema_version", metadata.get("schema", 1)) != SCHEMA_VERSION:
                    raise ValueError("Unsupported gene model catalog schema")
                if "loci" in metadata:
                    raise ValueError("Catalog metadata must omit loci; use loci.jsonl")
                species = str(metadata.get("species", ""))
                if not species or species in species_seen:
                    raise ValueError("Missing or duplicate catalog species: " + species)
                species_seen.add(species)
                database.execute("INSERT INTO catalogs(species,json) VALUES (?,?)", (species, _json(metadata)))
                content_hash = hashlib.sha256()
                species_loci = species_candidates = 0
                with loci_path.open("rb") as handle:
                    for line_number, line in enumerate(handle, start=1):
                        content_hash.update(line)
                        if not line.strip():
                            raise ValueError(f"Empty locus JSONL record: {loci_path}:{line_number}")
                        locus = json.loads(line)
                        if str(locus.get("species", species)) != species:
                            raise ValueError("Locus species disagrees with catalog: " + str(locus.get("species")))
                        locus["species"] = species
                        gene_id, seqid, strand = str(locus.get("gene_id", "")), str(locus.get("seqid", "")), str(locus.get("strand", ""))
                        records = locus.get("candidates", [])
                        if not gene_id or not records or not isinstance(records, list):
                            raise ValueError("Locus requires gene ID and coding candidates")
                        serialized = _json(locus)
                        try:
                            database.execute("INSERT INTO loci(species,gene_id,seqid,strand,json) VALUES (?,?,?,?,?)",
                                             (species, gene_id, seqid, strand, serialized))
                            for candidate in records:
                                identifier = str(candidate.get("candidate_id", ""))
                                if not identifier:
                                    raise ValueError("Candidate ID is empty")
                                if str(candidate.get("species", species)) != species or str(candidate.get("gene_id", gene_id)) != gene_id:
                                    raise ValueError("Candidate ownership disagrees with containing locus")
                                database.execute("INSERT INTO candidate_owners(species,candidate_id,gene_id) VALUES (?,?,?)",
                                                 (species, identifier, gene_id))
                        except sqlite3.IntegrityError as error:
                            raise ValueError(f"Duplicate locus or candidate ownership in {species}/{gene_id}") from error
                        largest = max(largest, len(serialized.encode()))
                        locus_count += 1
                        candidate_count += len(records)
                        species_loci += 1
                        species_candidates += len(records)
                for key, observed in (("loci", species_loci), ("candidates", species_candidates)):
                    expected = metadata.get("summary", {}).get(key)
                    if expected is not None and int(expected) != observed:
                        raise ValueError(f"Catalog {species} summary {key} count disagrees with loci.jsonl")
                if any(_identity(path) != identity for path, identity in before.items()):
                    raise OSError("Catalog source changed during locus-store import")
                sources.append({"species": species, "metadata_path": str(metadata_path), "loci_path": str(loci_path),
                                "metadata_sha256": hashlib.sha256(metadata_bytes).hexdigest(),
                                "loci_sha256": content_hash.hexdigest()})
            database.execute("PRAGMA user_version=" + str(SCHEMA_VERSION))
        database.close()
        os.replace(temporary, destination)
    except BaseException:
        database.close()
        temporary.unlink(missing_ok=True)
        raise
    return {"schema": SCHEMA_VERSION, "database": str(destination), "species_count": len(species_seen),
            "loci_count": locus_count, "candidate_count": candidate_count,
            "max_locus_json_bytes": largest, "sources": sources}


@contextlib.contextmanager
def _connection(db):
    owned = not isinstance(db, sqlite3.Connection)
    connection = sqlite3.connect(Path(db).expanduser().resolve().as_uri() + "?mode=ro", uri=True) if owned else db
    try:
        if owned:
            # Keep one consistent read snapshot and one read lock for this
            # scope, rather than renewing NFS locks for every locus query.
            connection.execute("BEGIN")
        if connection.execute("PRAGMA user_version").fetchone()[0] != SCHEMA_VERSION:
            raise ValueError("Unsupported gene model store schema")
        yield connection
    finally:
        if owned:
            connection.close()


def load_locus(db, species, gene_id):
    """Load one owning gene; absent loci raise KeyError rather than fallback."""
    with _connection(db) as connection:
        row = connection.execute("SELECT json FROM loci WHERE species=? AND gene_id=?", (species, gene_id)).fetchone()
        if row is None:
            raise KeyError((species, gene_id))
        return json.loads(row[0])


def iter_locus_keys(db, species=None):
    """Yield stable ``(species,gene_id)`` keys without loading coding sequences."""
    with _connection(db) as connection:
        cursor = connection.execute("SELECT species,gene_id FROM loci ORDER BY species,gene_id") if species is None else connection.execute(
            "SELECT species,gene_id FROM loci WHERE species=? ORDER BY gene_id", (species,))
        for row in cursor:
            yield row[0], row[1]


def iter_loci(db, species=None):
    """Yield one decoded locus at a time in stable species/gene order."""
    with _connection(db) as connection:
        cursor = connection.execute("SELECT json FROM loci ORDER BY species,gene_id") if species is None else connection.execute(
            "SELECT json FROM loci WHERE species=? ORDER BY gene_id", (species,))
        for row in cursor:
            yield json.loads(row[0])


def select_from_store(db, edges, *, max_component_loci=5000, **selectionkwargs):
    """Apply the selection engine by graph component without loading all CDSs.

    Oversized components abstain with a source-representative record for every
    locus; their sequences are decoded individually.  Isolated genes similarly
    use one-record batches.  Result decisions/scores/audits remain in memory but
    contain no CDS or protein strings.  The graph metadata cost is linear in
    loci and edges, and is reported separately from loaded component size.
    """
    if max_component_loci < 1:
        raise ValueError("max_component_loci must be positive")
    with _connection(db) as connection:
        parent = {node: node for node in iter_locus_keys(connection)}
        nonvoting, coordinate_ambiguity = set(), set()
        locus = None
        # Admission is streamed before component formation.  An excluded
        # coding locus must not bridge otherwise independent components or
        # push a valid component over its size cap.  Keep only key/boolean
        # metadata from this pass, never the decoded sequences.
        for locus in iter_loci(connection):
            node = str(locus["species"]), str(locus["gene_id"])
            if not locus_can_vote(locus, selectionkwargs.get("candidate_limit", 32)):
                nonvoting.add(node)
            if locus.get("ambiguous_coordinates", False):
                coordinate_ambiguity.add(node)
        del locus

        def find(node):
            while parent[node] != node:
                parent[node] = parent[parent[node]]
                node = parent[node]
            return node

        edge_rows = list(edges)
        prepared = prepare_correspondence(parent, edge_rows, coordinate_ambiguity=coordinate_ambiguity,
                                          nonvoting=nonvoting)
        voting_edges = [{"species_a": a[0], "gene_a": a[1], "species_b": b[0], "gene_b": b[1], "weight": weight}
                        for (a, b), weight in sorted(prepared["edge_by_pair"].items())]
        for edge in voting_edges:
            a, b = (edge["species_a"], edge["gene_a"]), (edge["species_b"], edge["gene_b"])
            root_a, root_b = find(a), find(b)
            if root_a != root_b:
                parent[max(root_a, root_b)] = min(root_a, root_b)
        components, component_edges = defaultdict(list), defaultdict(list)
        for node in sorted(parent):
            components[find(node)].append(node)
        for edge in voting_edges:
            component_edges[find((str(edge["species_a"]), str(edge["gene_a"])))].append(edge)
        result = None
        total_metrics = defaultdict(int)
        peak_loaded_loci = limited_components = 0

        def merge(part):
            nonlocal result
            if result is None:
                result = {key: value for key, value in part.items() if key not in {"metrics", "selections", "scores", "audit", "omitted_edges"}}
                result.update(selections=[], scores=[], audit=[], omitted_edges=[])
            for key in ("selections", "scores", "audit", "omitted_edges"):
                result[key].extend(part.get(key, []))
            for key, value in part["metrics"].items():
                if key in {"pair_cache_limit", "pair_cache_peak_entries"}:
                    total_metrics[key] = max(total_metrics[key], value)
                else:
                    total_metrics[key] += value

        for root, keys in sorted(components.items()):
            if len(keys) == 1:
                continue
            if len(keys) > max_component_loci:
                limited_components += 1
                for species, gene_id in keys:
                    locus = load_locus(connection, species, gene_id)
                    part = select_representatives([{"schema": 1, "species": species, "loci": [locus]}], [],
                                                  _copy_ambiguity={(species, gene_id)} & prepared["ambiguity"],
                                                  **selectionkwargs)
                    for selection in part["selections"]:
                        selection.update(status="insufficient_evidence", reason="component_limit_exceeded")
                    merge(part)
                    del locus
                peak_loaded_loci = max(peak_loaded_loci, 1)
                continue
            catalogs = defaultdict(list)
            for species, gene_id in keys:
                catalogs[species].append(load_locus(connection, species, gene_id))
            peak_loaded_loci = max(peak_loaded_loci, len(keys))
            part = select_representatives([{"schema": 1, "species": species, "loci": loci}
                                           for species, loci in sorted(catalogs.items())], component_edges[root],
                                          _copy_ambiguity={node for node in keys if node in prepared["ambiguity"]},
                                          **selectionkwargs)
            merge(part)
            del catalogs
        # Most loci have no voting neighbors. Stream their records through one
        # cursor instead of doing one indexed query and NFS read-lock cycle per
        # gene. Each locus still uses the unchanged single-locus selection path.
        for locus in iter_loci(connection):
            node = str(locus["species"]), str(locus["gene_id"])
            if len(components[find(node)]) != 1:
                continue
            part = select_representatives([{"schema": 1, "species": node[0], "loci": [locus]}], [],
                                          _copy_ambiguity={node} & prepared["ambiguity"],
                                          **selectionkwargs)
            merge(part)
            peak_loaded_loci = max(peak_loaded_loci, 1)
        if result is None:
            result = select_representatives([], [], **selectionkwargs)
            total_metrics.update(result["metrics"])
        result["selections"].sort(key=lambda row: (row["species"], row["gene_id"]))
        result["omitted_edges"] = prepared["omitted_edges"] if selectionkwargs.get("policy", "conserved") == "conserved" else []
        result["scores"].sort(key=lambda row: (row["species"], row["gene_id"], row["candidate_id"]))
        result["audit"].sort(key=lambda row: row["loci"])
        result["omitted_edges"].sort(key=lambda row: (row["a"], row["b"]))
        total_metrics["components"] = len(components)
        total_metrics.update(store_loci=len(parent), graph_edges=len(edge_rows),
                             max_component_loci=max_component_loci, peak_loaded_loci=peak_loaded_loci,
                             limited_components=limited_components)
        result["metrics"] = dict(total_metrics)
        return result
