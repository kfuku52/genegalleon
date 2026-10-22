#!/usr/bin/env python3
"""Apply explicit UTR coordinate corrections with independent sequence proof.

This opt-in curation does not change CDS, exon, RNA or gene models. It is not
part of automatic safe repair: unsupported proposals fail before publication.
"""

import argparse
import gzip
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path

SUPPORT = Path(__file__).resolve().parent
if str(SUPPORT) not in sys.path:
    sys.path.insert(0, str(SUPPORT))

from Bio.Data import CodonTable
from Bio.Seq import Seq
from fasta_sequence_store import fasta_records
from format_species_annotation.common import reverse_complement
from format_species_annotation.genome_intervals import genome_intervals
from format_species_annotation.locus_identity import rows
from format_species_annotation.source_overlap import sha256_file
from gff_feature_structure import ordered_annotated_blocks, ordered_feature_blocks
from input_validation_reuse import implementation_identity

UTR_TYPES = {"five_prime_utr": "five_prime", "three_prime_utr": "three_prime"}
CHANGE_KEYS = {"line", "transcript_id", "feature_id", "old_start", "old_end", "new_start", "new_end"}


def coordinate(row):
    return row[0], row[6], int(row[3]), int(row[4])


def selected_records(path, identifiers):
    found = {}
    for identifier, _header, sequence in fasta_records(Path(path)):
        if identifier in identifiers:
            if identifier in found:
                raise ValueError("Duplicate sequence identity: " + identifier)
            found[identifier] = sequence.upper()
    if set(found) != identifiers:
        raise ValueError("Missing exact transcript IDs in " + str(path))
    return found


def transcript_projection(exons, start, end):
    """Project a half-open spliced interval to closed genomic blocks."""
    result, offset = [], 0
    for seqid, strand, first, last in exons:
        length = last - first + 1
        a, b = max(start, offset), min(end, offset + length)
        if a < b:
            if strand == "+":
                result.append((seqid, strand, first + a - offset, first + b - offset - 1))
            else:
                result.append((seqid, strand, last - (b - offset) + 1, last - (a - offset)))
        offset += length
    return result


def coding_interval(exons, coding):
    """Require one contiguous coding interval in the explicit spliced RNA."""
    spans = []
    for seqid, strand, first, last in coding:
        offset, matches = 0, []
        for exon_seqid, exon_strand, exon_first, exon_last in exons:
            if (seqid, strand) == (exon_seqid, exon_strand) and exon_first <= first <= last <= exon_last:
                start = offset + (first - exon_first if strand == "+" else exon_last - last)
                matches.append((start, start + last - first + 1))
            offset += exon_last - exon_first + 1
        if len(matches) != 1:
            raise ValueError("CDS is not uniquely contained in an explicit exon")
        spans.extend(matches)
    if any(left[1] != right[0] for left, right in zip(spans, spans[1:], strict=False)):
        raise ValueError("CDS is not contiguous in explicit exons")
    return spans[0][0], spans[-1][1]


def model_proof(transcript, model, genome, supplied_cds, supplied_rna, genetic_code):
    anchors = model["rna"]
    if len(anchors) != 1 or anchors[0][2] not in {"mRNA", "transcript", "RNA"}:
        raise ValueError("One explicit RNA model is required: " + transcript)
    rna = coordinate(anchors[0])
    if not model["exon"] or not model["CDS"]:
        raise ValueError("Explicit exon and CDS blocks are required: " + transcript)
    exons = ordered_feature_blocks([coordinate(r) for r in model["exon"]], transcript)
    coding, mode = ordered_annotated_blocks(
        [(*coordinate(r), r[8]) for r in model["CDS"]], transcript)
    if mode != "cis" or any(block[:2] != rna[:2] for block in exons + coding):
        raise ValueError("An ordinary single-axis transcript is required: " + transcript)
    if (min(b[2] for b in exons), max(b[3] for b in exons)) != rna[2:]:
        raise ValueError("RNA bounds differ from explicit exons: " + transcript)
    start, end = coding_interval(exons, coding)

    def sequence(blocks):
        chunks = []
        for seqid, strand, first, last in blocks:
            chunk = genome.fetch(seqid, first - 1, last)
            chunks.append(chunk if strand == "+" else reverse_complement(chunk))
        return "".join(chunks)

    if sequence(exons) != supplied_rna:
        raise ValueError("Published RNA differs from genomic exons: " + transcript)
    if sequence(coding) != supplied_cds:
        raise ValueError("Supplied CDS differs from genomic CDS: " + transcript)
    table = CodonTable.unambiguous_dna_by_id[genetic_code]
    if (set(supplied_cds) - set("ACGT") or len(supplied_cds) % 3
            or supplied_cds[:3] not in table.start_codons
            or supplied_cds[-3:] not in table.stop_codons
            or "*" in str(Seq(supplied_cds[:-3]).translate(table=genetic_code))):
        raise ValueError("A complete unambiguous coding frame is required: " + transcript)
    offset, phases = 0, {}
    for block in coding:
        phases[block] = str((-offset) % 3)
        offset += block[3] - block[2] + 1
    if any(r[7] != phases[coordinate(r)] for r in model["CDS"]):
        raise ValueError("CDS phases are inconsistent: " + transcript)
    if supplied_rna.count(supplied_cds) != 1:
        raise ValueError("CDS must occur exactly once in the published RNA: " + transcript)
    if supplied_rna[start:end] != supplied_cds or transcript_projection(exons, start, end) != coding:
        raise ValueError("Published RNA coding interval differs from GFF: " + transcript)
    expected = {"five_prime": transcript_projection(exons, 0, start),
                "three_prime": transcript_projection(exons, end, len(supplied_rna))}
    return expected, {"transcript_id": transcript, "rna_length": len(supplied_rna),
                      "cds_length": len(supplied_cds), "cds_rna_offset": start,
                      "rna_sha256": hashlib.sha256(supplied_rna.encode()).hexdigest(),
                      "cds_sha256": hashlib.sha256(supplied_cds.encode()).hexdigest()}


def curate(gff, genome, cds, transcripts, proposal, outfile, audit, genetic_code=1):
    inputs = dict(gff=Path(gff), genome=Path(genome), cds=Path(cds),
                  transcripts=Path(transcripts), proposal=Path(proposal))
    outputs = [Path(outfile), Path(audit)]
    protected = {p.resolve() for p in inputs.values()}
    if (len({p.resolve() for p in outputs}) != 2
            or any(p.exists() or p.is_symlink() or p.resolve() in protected for p in outputs)):
        raise ValueError("Curation requires two distinct new output paths")
    before = {key: sha256_file(path) for key, path in inputs.items()}
    request = json.loads(inputs["proposal"].read_text())
    if (set(request) != {"schema_version", "gff_sha256", "changes"}
            or request["schema_version"] != 1 or request["gff_sha256"] != before["gff"]
            or not isinstance(request["changes"], list) or not 1 <= len(request["changes"]) <= 10000):
        raise ValueError("Proposal does not bind the exact source GFF")
    changes = {}
    for item in request["changes"]:
        if (not isinstance(item, dict) or set(item) != CHANGE_KEYS or any(type(item[key]) is not int or item[key] < 1
                for key in ("line", "old_start", "old_end", "new_start", "new_end"))
                or item["line"] in changes
                or any(not isinstance(item[key], str) or not item[key] for key in ("transcript_id", "feature_id"))
                or item["new_end"] < item["new_start"]):
            raise ValueError("Malformed or repeated UTR correction")
        changes[item["line"]] = item
    identifiers = {item["transcript_id"] for item in changes.values()}
    watched = identifiers | {item["feature_id"] for item in changes.values()}
    counts = dict.fromkeys(watched, 0)
    models = {name: {"rna": [], "exon": [], "CDS": [], "utr": []} for name in identifiers}
    selected = {}
    for number, fields, attrs in rows(inputs["gff"]):
        ids, parents = attrs.get("ID", ()), attrs.get("Parent", ())
        for value in set(ids) & watched:
            counts[value] += 1
        for name in set(ids) & identifiers:
            models[name]["rna"].append(fields)
        for name in set(parents) & identifiers:
            if len(parents) != 1:
                raise ValueError("Shared transcript children cannot be curated: " + name)
            kind = fields[2] if fields[2] in {"exon", "CDS"} else "utr" if fields[2].lower() in UTR_TYPES else ""
            if kind:
                models[name][kind].append(fields)
        if number in changes:
            item = changes[number]
            if (ids != (item["feature_id"],) or parents != (item["transcript_id"],)
                    or fields[2].lower() not in UTR_TYPES
                    or (int(fields[3]), int(fields[4])) != (item["old_start"], item["old_end"])):
                raise ValueError("Proposed source UTR row differs at line " + str(number))
            selected[number] = fields
    if set(selected) != set(changes) or any(count != 1 for count in counts.values()):
        raise ValueError("Curation requires unique explicit RNA and UTR identities")
    source_cds, source_rna = selected_records(inputs["cds"], identifiers), selected_records(inputs["transcripts"], identifiers)
    intervals = [coordinate(r) for model in models.values() for r in model["exon"] + model["CDS"]]
    proofs = {}
    replacements = {id(selected[n]): changes[n] for n in changes}
    with genome_intervals(inputs["genome"], [(s, a - 1, b) for s, _strand, a, b in intervals]) as reader:
        for name, model in models.items():
            expected, proofs[name] = model_proof(name, model, reader, source_cds[name], source_rna[name], genetic_code)
            occupied = set()
            for row in model["utr"]:
                seqid, strand, first, last = coordinate(row)
                update = replacements.get(id(row))
                if update:
                    first, last = update["new_start"], update["new_end"]
                    if (seqid, strand, update["old_start"], update["old_end"]) in expected[UTR_TYPES[row[2].lower()]]:
                        raise ValueError("Proposed UTR already matches the sequence-supported model")
                block = seqid, strand, first, last
                if block not in expected[UTR_TYPES[row[2].lower()]] or block in occupied:
                    raise ValueError("UTR coordinates lack unique sequence support: " + name)
                occupied.add(block)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=outputs[0].parent, delete=False) as raw:
            temporary = Path(raw.name)
            writer = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) if str(outputs[0]).endswith(".gz") else raw
            opener = gzip.open if str(inputs["gff"]).endswith(".gz") else open
            with opener(inputs["gff"], "rt", newline="") as source:
                for number, line in enumerate(source, 1):
                    if number in changes:
                        fields = line.rstrip("\r\n").split("\t")
                        fields[3:5] = [str(changes[number][key]) for key in ("new_start", "new_end")]
                        line = "\t".join(fields) + ("\r\n" if line.endswith("\r\n") else "\n")
                    writer.write(line.encode())
            if writer is not raw:
                writer.close()
        after = {key: sha256_file(path) for key, path in inputs.items()}
        if after != before:
            raise ValueError("Curation input changed during proof/publication")
        result = {"schema_version": 1, "status": "sequence_verified_utr_curation", "genetic_code": genetic_code,
                  "implementation_sha256": implementation_identity(),
                  "inputs": {key: {"path": str(path.resolve()), "sha256": before[key]} for key, path in inputs.items()},
                  "changes": list(changes.values()), "models": proofs, "cds_file_unchanged": True,
                  "cds_exon_rna_gene_features_unchanged": True, "output_sha256": sha256_file(temporary)}
        with outputs[1].open("x") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
        os.link(temporary, outputs[0])
        return result
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("gff", "genome", "cds", "transcripts", "proposal", "outfile", "audit"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--genetic-code", type=int, default=1)
    args = vars(parser.parse_args())
    try:
        result = curate(**args)
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    print(json.dumps(result))


if __name__ == "__main__":
    main()
