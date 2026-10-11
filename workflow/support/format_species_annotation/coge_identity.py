"""Prove CoGe orphan CDS ownership and split feature identity from author RNA aliases."""

from collections import defaultdict
from functools import lru_cache
from pathlib import Path

from format_species_writers import open_text

from .common import parse_gff_attributes


def coge_identity_evidence(path):
    source = Path(path).resolve()
    stat = source.stat()
    identity = (stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
    result = _evidence(str(source), identity)
    stat = source.stat()
    if identity != (stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns):
        raise OSError("CoGe annotation changed during identity qualification")
    return result


@lru_cache(maxsize=4)
def _evidence(path, identity):
    genes, rnas, blocks = {}, {}, defaultdict(list)
    with open_text(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) != 9 or parts[1].lower() != "coge":
                continue
            attrs = parse_gff_attributes(parts[8])
            ids = attrs.get("ID", ())
            row = (parts[0], parts[6], int(parts[3]), int(parts[4]), attrs)
            kind = parts[2].lower()
            if kind == "gene" or kind.endswith("rna") or "transcript" in kind:
                if len(ids) != 1:
                    continue
                target = genes if kind == "gene" else rnas
                if ids[0] in target and target[ids[0]] != row:
                    raise ValueError("Conflicting CoGe gene/RNA identity: " + ids[0])
                target[ids[0]] = row
            elif kind == "cds" and not attrs.get("Parent") and attrs.get("coge_fid"):
                names = attrs.get("Name", ())
                if len(names) == 1:
                    blocks[names[0]].append(row)
    alias_rnas = defaultdict(set)
    for identifier, row in rnas.items():
        for alias in (identifier,) + row[4].get("Alias", ()) + row[4].get("Name", ()):
            alias_rnas[alias].add(identifier)
    owners, splits = {}, set()
    for transcript, rows in blocks.items():
        candidates = set(alias_rnas.get(transcript, ()))
        explicit = []
        for row in rows:
            declared = set(row[4].get("Alias", ())) & rnas.keys()
            candidates.update(declared)
            explicit.append(declared)
        parents = {parent for rna in candidates for parent in rnas[rna][4].get("Parent", ())}
        if not parents:
            continue
        if len(parents) != 1:
            raise ValueError("Ambiguous CoGe CDS gene aliases: " + transcript)
        gene = next(iter(parents))
        if gene not in genes:
            continue
        anchor = genes[gene]
        if any(row[:2] != anchor[:2] or not 1 <= anchor[2] <= row[2] <= row[3] <= anchor[3]
               for row in rows):
            raise ValueError("CoGe CDS gene alias disagrees with coordinates: " + transcript)
        owners[transcript] = gene
        # Internal CoGe feature IDs can label pieces of one transcript. A
        # shared exact RNA ID in every CDS Alias, not a common name or gene,
        # is required before joining different feature IDs.
        shared = set.intersection(*explicit) if explicit else set()
        if len(shared) != 1:
            continue
        rna = rnas[next(iter(shared))]
        if any(row[:2] != rna[:2] or not rna[2] <= row[2] <= row[3] <= rna[3] for row in rows):
            continue
        positions = sorted((row[2], row[3]) for row in rows)
        metadata = {tuple(sorted((key, values) for key, values in row[4].items()
                                 if key not in {"ID", "coge_fid"})) for row in rows}
        if len(metadata) == 1 and all(left[1] < right[0] for left, right in zip(positions, positions[1:], strict=False)):
            splits.add(transcript)
    return {"gene_owners": owners, "split_transcripts": frozenset(splits)}
