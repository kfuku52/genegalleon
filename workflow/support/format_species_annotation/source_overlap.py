"""Prove source CDS overlaps without inventing a biological exception."""

import hashlib
from collections import defaultdict
from pathlib import Path

from gff_feature_structure import ordered_feature_blocks

from .common import (
    first_token,
    iter_fasta_records,
    parse_gff_attributes,
    reverse_complement,
)
from .grouping import extract_cds_header_alias_tiers
from .organelle import iter_non_organelle_gff_lines


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_overlap_key(attributes):
    attrs = parse_gff_attributes(attributes)
    ids, parents = attrs.get("ID", ()), attrs.get("Parent", ())
    if len(ids) != 1 or len(parents) != 1 or not ids[0] or not parents[0]:
        return None
    return ids[0], parents[0]


def mark_source_overlap(attributes, confirmed):
    # This reserved marker is produced by this audit, never trusted from an
    # unaudited input file. Keep other attributes and their encoding intact.
    if not confirmed and not any(field.partition("=")[0].strip() == "gg_source_overlap"
                                 for field in attributes.split(";")):
        return attributes
    fields = [field for field in attributes.split(";")
              if field.partition("=")[0].strip() != "gg_source_overlap" and field]
    if confirmed:
        fields.append("gg_source_overlap=confirmed")
    return ";".join(fields)


def audit_source_overlaps(gff_path, task):
    groups = defaultdict(list)
    for line in iter_non_organelle_gff_lines(gff_path):
        parts = line.rstrip("\n\r").split("\t")
        if len(parts) != 9 or parts[2].lower() != "cds":
            continue
        key = source_overlap_key(parts[8])
        if key is not None:
            groups[key].append(parts)
    candidates = {}
    for key, rows in groups.items():
        if any(parse_gff_attributes(row[8]).get(name) for row in rows
               for name in ("exception", "pseudo")):
            continue  # Declared biological exceptions retain their own contract.
        try:
            blocks = [(row[0], row[6], int(row[3]), int(row[4])) for row in rows]
        except ValueError:
            continue  # Preserve malformed rows for the normal strict validator.
        if (len({block[:2] for block in blocks}) != 1
                or any(strand not in ("+", "-") or start < 1 or end < start
                       for _seqid, strand, start, end in blocks)):
            continue  # Ordered fragments/trans-splicing have separate rules.
        ordered = ordered_feature_blocks(blocks, key[0], allow_overlap=True)
        genomic = sorted(ordered, key=lambda block: (block[2], block[3]))
        if any(right[2] <= left[3] for left, right in zip(genomic, genomic[1:], strict=False)):
            candidates[key] = (rows, ordered)
    audit = {"version": 1, "confirmed": [], "unconfirmed": []}
    if not candidates:
        return set(), audit
    source_cds, genome = task.get("cds_path"), task.get("genome_path")
    if source_cds is None or genome is None:
        audit["unconfirmed"] = [{"cds_id": key[0], "parent": key[1],
                                  "reason": "source CDS and genome are required"} for key in candidates]
        return set(), audit
    paths = {"gff": Path(gff_path), "cds": Path(source_cds), "genome": Path(genome)}
    audit["inputs"] = {name: {"path": str(path.resolve()), "sha256": sha256_file(path)}
                       for name, path in paths.items()}
    required_seqids = {block[0] for _rows, blocks in candidates.values() for block in blocks}
    sequences = {}
    for header, sequence in iter_fasta_records(genome):
        identifier = first_token(header)
        if identifier in required_seqids:
            if identifier in sequences:
                raise ValueError("Duplicate source genome ID: " + identifier)
            sequences[identifier] = sequence.upper()
    source_records = defaultdict(list)
    for header, sequence in iter_fasta_records(source_cds):
        aliases = set().union(*extract_cds_header_alias_tiers(task, header))
        digest = hashlib.sha256(sequence.upper().encode("ascii")).hexdigest()
        source_records[digest].append((first_token(header), aliases))
    confirmed = set()
    for key, (rows, blocks) in candidates.items():
        identities = set(key)
        for row in rows:
            attrs = parse_gff_attributes(row[8])
            for name in ("protein_id", "transcript_id", "locus_tag", "gene_id"):
                identities.update(attrs.get(name, ()))
        record = {"cds_id": key[0], "parent": key[1], "blocks": [list(block) for block in blocks]}
        try:
            chunks = []
            for seqid, strand, start, end in blocks:
                sequence = sequences[seqid]
                if end > len(sequence):
                    raise ValueError("GFF coordinate exceeds source genome")
                chunk = sequence[start - 1:end]
                chunks.append(reverse_complement(chunk) if strand == "-" else chunk)
            joined = "".join(chunks)
            # The first phase trims a partial 5' codon, never every exon.
            first = next(row for row in rows if (row[0], row[6], int(row[3]), int(row[4])) == blocks[0])
            if first[7] not in ("0", "1", "2"):
                raise ValueError("Source overlap requires an explicit first CDS phase")
            joined = joined[int(first[7]):]
            digest = hashlib.sha256(joined.encode("ascii")).hexdigest()
            matches = sorted(identifier for identifier, aliases in source_records.get(digest, ())
                             if aliases & identities)
            if not matches:
                raise ValueError("Complete source CDS sequence and feature identity do not match")
            record.update(source_cds_ids=matches, cds_sha256=digest, length=len(joined), first_phase=int(first[7]))
            audit["confirmed"].append(record)
            confirmed.add(key)
        except (KeyError, ValueError) as exc:
            record["reason"] = str(exc)
            audit["unconfirmed"].append(record)
    for name, path in paths.items():
        if sha256_file(path) != audit["inputs"][name]["sha256"]:
            raise ValueError("Source changed during overlap audit: " + name)
    return confirmed, audit
