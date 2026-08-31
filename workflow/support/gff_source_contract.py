"""Resolve duplicate species annotations from the exact indexed FASTA source.

Only an unambiguous conventional CDS/annotation basename pair is accepted.
This is not a best-prefix search: multiple releases or unknown naming schemes
require an explicit input correction, rather than an arbitrary annotation.
"""

import json
import re
import sqlite3
import zlib
from pathlib import Path


def paired_gff_name(source_path, candidates):
    name = Path(source_path).name
    stem = re.split(r"\.(?:cds|pep|protein)(?=\.|_|$)", name, maxsplit=1, flags=re.IGNORECASE)[0]
    if stem == name:
        raise ValueError("Cannot determine an exact GFF source pair for {}".format(name))
    pattern = re.compile(re.escape(stem) + r"(?:\.\d+)?(?:\.gene(?:_exons)?)?\.(?:gff3?|gtf)(?:\.gz)?$")
    matches = sorted(candidate for candidate in candidates if pattern.fullmatch(candidate))
    if len(matches) != 1:
        raise ValueError("Ambiguous GFF source pair for {}: {}".format(name, ", ".join(matches) or "none"))
    return matches[0]


def source_bound_gff_names(database, sequences, candidates):
    """Return gene -> GFF only after exact identifier and sequence agreement.

    candidates maps each requested FASTA identifier to same-species GFF names.
    The existing sequence store is read-only and must still describe its inputs.
    """
    selected = {}
    connection = sqlite3.connect(Path(database).resolve().as_uri() + "?mode=ro", uri=True)
    try:
        connection.execute("PRAGMA query_only=ON")
        connection.execute("BEGIN")
        metadata = dict(connection.execute("SELECT key,value FROM metadata"))
        if metadata.get("schema_version") != "2":
            raise ValueError("Unsupported FASTA sequence-store source contract")
        sources = json.loads(metadata["source_signatures"])
        checked_sources = set()
        for identifier, sequence in sequences.items():
            rows = connection.execute(
                "SELECT sequence,source_order FROM sequences WHERE identifier=?", (identifier,)
            ).fetchall()
            if len(rows) != 1 or zlib.decompress(rows[0][0]).decode().upper() != sequence.upper():
                raise ValueError("GFF source FASTA identity/sequence mismatch for {}".format(identifier))
            order = rows[0][1]
            if not isinstance(order, int) or not 0 <= order < len(sources):
                raise ValueError("Invalid indexed FASTA source for {}".format(identifier))
            source = sources[order]
            if order not in checked_sources:
                path = Path(source["path"])
                current = path.stat()
                fields = {"device": current.st_dev, "inode": current.st_ino, "size": current.st_size,
                          "mtime_ns": current.st_mtime_ns, "ctime_ns": current.st_ctime_ns}
                if path.is_symlink() or any(source.get(key) != value for key, value in fields.items()):
                    raise ValueError("Indexed GFF source FASTA changed: {}".format(path.name))
                checked_sources.add(order)
            selected[identifier] = paired_gff_name(source["path"], candidates[identifier])
    finally:
        connection.close()
    return selected
