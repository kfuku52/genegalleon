#!/usr/bin/env python3
"""Canonical nucleotide-reader coordinates for immutable selected bundles."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import os
import re
import sys
import tempfile
from pathlib import Path
from urllib.parse import quote

if __package__:
    from . import gff2genestat as genestat
    from .fasta_sequence_store import exclusive_lock, fasta_records
    from .representative_selection import RepresentativeMap, file_digest, load_effective_inputs
else:
    import gff2genestat as genestat
    from fasta_sequence_store import exclusive_lock, fasta_records
    from representative_selection import RepresentativeMap, file_digest, load_effective_inputs

GFF_COLUMNS = ("sequence", "source", "feature", "start", "end", "score", "strand", "phase", "attributes")
HELPERS = ("fractionation_selected_inputs.py", "gff2genestat.py", "representative_selection.py",
           "gff_feature_structure.py", "species_labeling.py", "fasta_sequence_store.py",
           "gff_source_contract.py", "content_digest_cache.py", "gff_attribute_syntax.py",
           "format_species_writers.py", "format_species_common.py", "format_species_constants.py",
           "format_species_provider_config.py", "format_species_taxonomy.py")


def implementation():
    support = Path(__file__).resolve().parent
    paths = [support / name for name in HELPERS]
    paths.extend(sorted((support / "format_species_annotation").glob("*.py")))
    return {str(path.relative_to(support)): file_digest(path) for path in paths}


def runtime_identity():
    return {"python": sys.version, "numpy": genestat.numpy.__version__, "pandas": genestat.pandas.__version__}


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _identifier(value):
    if not value or value.startswith("#") or ";" in value or "||" in value or any(c.isspace() or ord(c) < 32 or ord(c) == 127 for c in value):
        raise ValueError("Invalid qualified fractionation gene identifier: " + repr(value))
    return value


def _new_id(kind, identity, used):
    stem = "ggfractionation_" + kind + "_" + hashlib.sha256(_json(identity).encode()).hexdigest()[:24]
    candidate, suffix = stem, 0
    while candidate in used:
        suffix += 1
        candidate = stem + "_" + str(suffix)
    used.add(candidate)
    return candidate


def _phase(value):
    if str(value).strip() in {"", ".", "nan"}:
        return "."
    numeric = float(value)
    if numeric not in {0, 1, 2}:
        raise ValueError("Invalid source CDS phase")
    return str(int(numeric))


def canonical_gff(species, source, choices):
    """Use genestat's explicit map selection; never choose longest anew."""
    lengths = {}
    for identifier, _header, sequence in fasta_records(Path(source["cds"])):
        _identifier(identifier)
        if identifier in lengths or not sequence or re.search(r"[^ACGTURYSWKMBDHVNacgturyswkmbdhvn]", sequence):
            raise ValueError("Invalid selected nucleotide FASTA: " + identifier)
        lengths[identifier] = len(sequence)
    if not lengths:
        raise ValueError("Selected fractionation CDS is empty: " + species)
    frame = genestat.read_gff_table(source["gff"])
    if frame.shape[1] != 9:
        raise ValueError("Selected fractionation GFF requires nine columns")
    frame.columns = GFF_COLUMNS
    with contextlib.redirect_stdout(sys.stderr):
        selected = genestat.extract_by_ids(frame, list(lengths), "CDS", "longest",
                                          representative_map=choices, species=species)
    if set(selected["gene_id"]) != set(lengths):
        raise ValueError("Selected fractionation GFF/CDS membership differs")
    lines, audit, used = ["##gff-version 3\n"], [], set(lengths)
    for gene_id, group in selected.groupby("gene_id", sort=True):
        choice = choices.choice(gene_id, species)
        source_transcript = choice["source_transcript_id"]
        selected_ids = set(group["selected_transcript"])
        if selected_ids != {source_transcript}:
            # Parentless CDS can use its locus as the reader's model key.
            # Admit this representation only when every source feature itself
            # has the exact requested identity, never a different alias.
            source_ids = {genestat.parse_attribute_fields(value)[0] for value in group["attributes"]}
            if source_ids != {source_transcript}:
                raise ValueError("Fractionation transcript identity differs from representative map: " + gene_id)
        blocks, mode = genestat.transcript_blocks(group, gene_id)
        axes = {(seqid, strand) for seqid, strand, _start, _end in blocks}
        if len(axes) != 1 or next(iter(axes))[1] not in {"+", "-"}:
            raise ValueError("Fractionation coordinates require one source contig and strand")
        if sum(end - start + 1 for _seqid, _strand, start, end in blocks) != lengths[gene_id]:
            raise ValueError("Fractionation selected CDS length differs from its exact source path: " + gene_id)
        seqid, strand = next(iter(axes))
        start, end = min(block[2] for block in blocks), max(block[3] for block in blocks)
        transcript = _new_id("transcript", [species, gene_id, source_transcript], used)
        provenance = "gg_source_transcript_id=" + quote(source_transcript, safe="._-:")
        gene_attribute = quote(gene_id, safe="._-:")
        lines.append(f"{seqid}\tGeneGalleon\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene_attribute};{provenance}\n")
        lines.append(f"{seqid}\tGeneGalleon\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={transcript};Parent={gene_attribute};{provenance}\n")
        parts = sorted({(str(row.sequence), str(row.strand), int(row.start), int(row.end), _phase(row.phase))
                        for row in group.itertuples(index=False)}, key=lambda part: (part[2], part[3], part[4]))
        for index, (block_seqid, block_strand, block_start, block_end, phase) in enumerate(parts):
            child = _new_id("cds", [species, gene_id, index, block_start, block_end, phase], used)
            lines.append(f"{block_seqid}\tGeneGalleon\tCDS\t{block_start}\t{block_end}\t.\t{block_strand}\t{phase}\tID={child};Parent={transcript};{provenance}\n")
        audit.append({"gene_id": gene_id, "candidate_id": choice["candidate_id"],
                      "source_transcript_id": source_transcript, "reader_transcript_id": transcript,
                      "seqid": seqid, "strand": strand, "cds_length": lengths[gene_id],
                      "source_blocks_1based": [list(block) for block in blocks], "source_splice_mode": mode})
    return "".join(lines), audit


def prepare_selected_inputs(inputs, species_names, cache_root):
    manifest = Path(inputs).resolve()
    rows = load_effective_inputs(manifest)
    manifest_hash, helpers, runtime = file_digest(manifest), implementation(), runtime_identity()
    output = []
    for species in species_names:
        if species not in rows:
            raise ValueError("Selected input bundle lacks fractionation species: " + species)
        source = rows[species]
        identities = {role: {"path": source[role], "sha256": source[role + "_sha256"]}
                      for role in ("cds", "gff", "representative_map")}
        expected = {"schema": 1, "view": "fractionation_nucleotide_coordinates", "species": species,
                    "manifest": {"path": str(manifest), "sha256": manifest_hash},
                    "inputs": identities, "genetic_code": int(source["genetic_code"]),
                    "implementation": helpers, "runtime": runtime}
        key = hashlib.sha256(_json(expected).encode()).hexdigest()
        parent = Path(cache_root).resolve() / "representative" / manifest_hash / "fractionation"
        parent.mkdir(parents=True, exist_ok=True)
        directory = parent / key
        receipt_path, gff_path = directory / "receipt.json", directory / "reader.gff3"
        with exclusive_lock(parent / (key + ".lock")):
            if directory.exists():
                try:
                    receipt = json.loads(receipt_path.read_text())
                    if receipt["identity"] != expected or receipt["key"] != key or file_digest(gff_path) != receipt["gff_sha256"]:
                        raise ValueError("Changed fractionation selected reader view")
                except (OSError, KeyError, json.JSONDecodeError) as exc:
                    raise ValueError("Incomplete fractionation selected reader view") from exc
            else:
                text, audit = canonical_gff(species, source, RepresentativeMap(source["representative_map"]))
                if (file_digest(manifest) != manifest_hash or implementation() != helpers or runtime_identity() != runtime or
                        any(file_digest(item["path"]) != item["sha256"] for item in identities.values())):
                    raise ValueError("Fractionation selected sources changed during adaptation")
                with tempfile.TemporaryDirectory(prefix=".fractionation-view-", dir=parent) as temporary:
                    stage = Path(temporary)
                    (stage / "reader.gff3").write_text(text)
                    receipt = {"key": key, "identity": expected, "gff_sha256": file_digest(stage / "reader.gff3"),
                               "mapping": audit, "reader_feature": "gene", "reader_attribute": "ID"}
                    (stage / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
                    os.replace(stage, directory)
        output.append({"species": species, "cds": source["cds"], "gff": str(gff_path), "receipt": str(receipt_path)})
    if (file_digest(manifest) != manifest_hash or implementation() != helpers or runtime_identity() != runtime or
            any(file_digest(rows[species][role]) != rows[species][role + "_sha256"]
                for species in species_names for role in ("cds", "gff", "representative_map"))):
        raise ValueError("Fractionation selected identity changed before delivery")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", required=True, type=Path)
    parser.add_argument("--species", required=True, action="append")
    parser.add_argument("--cache-root", required=True, type=Path)
    args = parser.parse_args()
    for row in prepare_selected_inputs(args.inputs, args.species, args.cache_root):
        print("\t".join(row[key] for key in ("cds", "gff", "receipt")))


if __name__ == "__main__":
    main()
