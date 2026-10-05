"""Preserve every annotated coding isoform with immutable source provenance.

Coordinates are zero based, half open, and listed in transcript order. GFF
phases describe bases before the next complete codon; internal CDS blocks are
concatenated without deleting phase bases. No source sequence is padded,
masked, corrected, or edited by this module.
"""

from __future__ import annotations

import contextlib
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

from Bio.Data import CodonTable
from Bio.Seq import Seq
from cds_model_normalisation import CdsModelNormaliser, phase_conventions, protein
from fasta_sequence_store import fasta_records
from format_species_annotation.common import parse_gff_attributes, sanitize_identifier
from format_species_annotation.grouping import (
    build_gff_cds_grouping_index,
    extract_cds_header_alias_tiers,
    resolve_cds_header_gff_gene,
)
from format_species_annotation.organelle import fasta_header_is_organelle
from gff_attribute_syntax import validate_gff

SCHEMA_VERSION = 1


def source_signature(path):
    """Hash a stable source; reject mutation while reading it."""
    path = Path(path).expanduser().resolve()
    before = path.stat()
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    after = path.stat()
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (
            after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns):
        raise OSError("Source changed while cataloging: " + str(path))
    return {"path": str(path), "size": after.st_size, "sha256": digest}


def _stat_identity(path):
    stat = Path(path).stat()
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


def _formatted_id(species, source_id):
    value = sanitize_identifier(source_id)
    return value if value.startswith(species + "_") else sanitize_identifier(species + "_" + value)


@contextlib.contextmanager
def indexed_genome(path):
    """Index a scratch link/copy, keeping genome FASTA and its directory intact."""
    import pysam

    with tempfile.TemporaryDirectory(prefix="gg-isoform-genome-") as scratch:
        staged = Path(scratch) / "genome.fa"
        path = Path(path).resolve()
        if path.name.lower().endswith(".gz"):
            with gzip.open(path, "rb") as source, staged.open("wb") as destination:
                shutil.copyfileobj(source, destination)
        else:
            staged.symlink_to(path)
        indexing = subprocess.run([sys.executable, "-c", "import pysam,sys; pysam.faidx(sys.argv[1])", str(staged)],
                                  capture_output=True, text=True)
        if indexing.returncode or indexing.stderr.strip():
            raise ValueError("FASTA index warning or failure: " + indexing.stderr.strip())
        with pysam.FastaFile(str(staged)) as genome:
            yield genome


def reconstruct_sequence(blocks, seqid, strand, genome):
    """Fetch biological CDS without per-exon phase trimming."""
    if strand not in {"+", "-"} or seqid not in genome.references:
        raise ValueError("Unknown CDS contig or strand: " + seqid)
    pieces = []
    for start, end, _phase in blocks:
        if start < 0 or end <= start or end > genome.get_reference_length(seqid):
            raise ValueError("CDS outside genome bounds: " + seqid)
        piece = genome.fetch(seqid, start, end).upper()
        pieces.append(str(Seq(piece).reverse_complement()) if strand == "-" else piece)
    return "".join(pieces)


def validate_candidate(candidate, genetic_code=1):
    """Return independent translation/structure quality for a coding candidate.

    Partial source models remain partial. Only their known first-phase offset
    and incomplete final codon are ignored when translating; CDS itself stays
    intact. ``usable`` is comparative eligibility, not a claim of a complete
    or experimentally confirmed ORF.
    """
    code = int(genetic_code)
    if code not in CodonTable.unambiguous_dna_by_id:
        raise ValueError("Unsupported genetic code: " + str(code))
    sequence = str(candidate.get("cds", "")).upper()
    blocks = candidate.get("blocks", [])
    previous = candidate.get("quality", {})
    structure_problem = previous.get("structure_problem", "")
    ascending = sorted(blocks)
    if not structure_problem:
        if candidate.get("strand", "+") not in {"+", "-"}:
            structure_problem = "invalid_strand"
        elif any(block[0] < 0 or block[1] <= block[0] for block in blocks):
            structure_problem = "invalid_cds_coordinates"
        elif any(a[1] > b[0] for a, b in zip(ascending, ascending[1:], strict=False)):
            structure_problem = "overlapping_cds"
        elif sum(block[1] - block[0] for block in blocks) != len(sequence):
            structure_problem = "cds_length_coordinate_mismatch"
    phases = [str(block[2]) for block in blocks]
    first = int(phases[0]) if phases and phases[0] in {"0", "1", "2"} else 0
    unknown = any(phase not in {"0", "1", "2"} for phase in phases)
    allowed_first = {first} if phases and phases[0] in {"0", "1", "2"} else {0, 1, 2}
    cumulative = 0
    for block in blocks:
        if str(block[2]) in {"0", "1", "2"}:
            allowed_first.intersection_update({(int(block[2]) + cumulative) % 3})
        cumulative += block[1] - block[0]
    phase_conflict = not blocks or not allowed_first
    phase_offset = first if not phase_conflict else 0
    coding = sequence[phase_offset:]
    codon_table = CodonTable.unambiguous_dna_by_id[code]
    dual_codons = set(codon_table.stop_codons) & set(codon_table.forward_table)
    uncertain_codons = (sorted({coding[position:position + 3] for position in range(0, len(coding) - 2, 3)
                               if coding[position:position + 3] in dual_codons}) if dual_codons else [])
    # Existing GG translators use context-free codon meaning, rather than
    # complete-CDS initiation/termination semantics. A dual stop/amino-acid
    # codon cannot establish an intact ORF under that contract, even at the end.
    translation_uncertain = bool(previous.get('translation_uncertain') or uncertain_codons)
    incomplete_codon = bool(len(coding) % 3)
    invalid_base = any(base not in "ACGTRYSWKMBDHVN" for base in sequence)
    translated = (str(Seq(coding[:len(coding) // 3 * 3]).translate(table=code))
                  if len(coding) >= 3 and not invalid_base else "")
    internal_stop = "*" in translated.removesuffix("*")
    starts = codon_table.start_codons
    stops = codon_table.stop_codons
    has_start = bool(coding) and coding[:3] in starts
    has_stop = len(coding) >= 3 and not incomplete_codon and coding[-3:] in stops
    partial = bool(previous.get("annotated_partial") or first or incomplete_codon or not has_start or not has_stop)
    ambiguous = any(base not in "ACGT" for base in sequence)
    sequence_mismatch = bool(previous.get("sequence_mismatch", False))
    exception = previous.get("annotated_exception", "")
    usable = bool(translated) and not any((phase_conflict, unknown, internal_stop, ambiguous, invalid_base,
                                         sequence_mismatch, structure_problem, exception, translation_uncertain))
    quality = {**previous, "structure_problem": structure_problem,
               "valid_orf": usable and not partial, "internal_stop": internal_stop,
               "partial": partial, "ambiguous": ambiguous, "invalid_base": invalid_base,
               "phase_conflict": phase_conflict, "sequence_mismatch": sequence_mismatch,
               "phase_unknown": unknown or bool(previous.get("phase_unknown")),
               "phase_unresolved": unknown, "phase_inferred": bool(previous.get("phase_inferred")),
               "has_start": has_start, "has_stop": has_stop, "incomplete_codon": incomplete_codon,
               "translation_offset": phase_offset, "genetic_code": code, "usable": usable,
               "translation_uncertain": translation_uncertain,
               "translation_uncertain_codons": uncertain_codons,
               "translation_convention": "context_free_codons_without_terminal_definite_stop",
               "annotated_pseudogene": exception == "annotated_pseudogene",
               "translation_exception": exception == "annotated_translation_exception",
               "sequence_exception": exception == "annotated_sequence_exception"}
    candidate["protein"] = translated.removesuffix("*")
    return quality


def _source_matches(sequence, candidate):
    return bool(_source_convention(sequence, candidate))


def _source_convention(sequence, candidate):
    """Recognise exact supplied/genomic identity and bounded stop conventions."""
    biological = candidate["cds"]
    offset = candidate["quality"]["translation_offset"]
    partial = candidate["quality"]["partial"]
    if sequence == biological:
        return "exact_genomic_cds"
    if any(sequence == biological + "N" * count for count in (1, 2)):
        return "terminal_codon_padding"
    code = candidate["quality"]["genetic_code"]
    terminal_stop = (len(biological) >= 3 and not len(biological) % 3
                     and biological[-3:] in CodonTable.unambiguous_dna_by_id[code].stop_codons)
    if terminal_stop and sequence == biological[:-3]:
        return "omitted_terminal_stop"
    if terminal_stop and sequence == biological[:-3] + "NNN":
        return "masked_terminal_stop"
    if offset and partial:
        coding = biological[offset:]
        coding = coding[:len(coding) // 3 * 3]
        if sequence == coding or any(sequence == coding + "N" * count for count in (1, 2)):
            return "annotated_partial_frame"
    return ""


def _infer_missing_phases(candidate, genetic_code):
    """Infer continuation phases only from a source-bound complete genomic ORF."""
    quality = candidate["quality"]
    if (not quality["phase_unresolved"] or not quality.get("source_sequence_agreement")
            or quality["partial"] or any(quality.get(key) for key in
                                        ("phase_conflict", "internal_stop", "ambiguous", "invalid_base",
                                         "sequence_mismatch", "structure_problem", "annotated_exception",
                                         "translation_uncertain"))):
        return
    cumulative = 0
    for block in candidate["blocks"]:
        expected = (3 - cumulative % 3) % 3
        if block[2] in {0, 1, 2} and block[2] != expected:
            return
        cumulative += block[1] - block[0]
    candidate["source_blocks"] = [list(block) for block in candidate["blocks"]]
    cumulative = 0
    for block in candidate["blocks"]:
        block[2] = (3 - cumulative % 3) % 3
        cumulative += block[1] - block[0]
    quality["phase_inferred"] = True
    quality["phase_inference_evidence"] = "complete_genomic_cds_and_bound_source"
    candidate["quality"] = validate_candidate(candidate, genetic_code)
    candidate["junctions"] = _junctions(candidate["blocks"], candidate["strand"])


def _corrected_cds_length(sequence, candidate, global_convention):
    """Infer only formatter corrections with the same bound-model evidence.

    Valid supplied translations are unchanged by CdsModelNormaliser. Partial
    frame trimming requires its unambiguous GFF3 convention. UTR mismatches or
    uncertain phase conventions intentionally have no inferred corrected length.
    """
    quality = candidate["quality"]
    code = quality["genetic_code"]
    valid_bases = not (set(sequence) - set("ACGTRYSWKMBDHVN"))
    if quality.get("annotated_exception") or (valid_bases and protein(sequence, code)):
        return len(sequence)
    if not _source_matches(sequence, candidate) or quality["phase_conflict"]:
        return None
    blocks = candidate["blocks"]
    conventions = phase_conventions([{"start": a, "end": b, "phase": str(p)} for a, b, p in blocks])
    convention = next(iter(conventions)) if len(conventions) == 1 else (
        global_convention if global_convention in conventions else None)
    first = blocks[0][2]
    if first and convention == "gff3" and quality["usable"]:
        return (len(candidate["cds"]) - first) // 3 * 3
    if not first and sequence == candidate["cds"] and quality["usable"]:
        return len(sequence)
    return None


def _attribute_parser(text):
    parsed = parse_gff_attributes(text)
    for key in ("ID", "gene_id", "transcript_id"):
        if len(set(parsed.get(key, ()))) > 1:
            raise ValueError("Conflicting structural GFF attribute: " + key)
    return {key: ",".join(dict.fromkeys(values)) for key, values in parsed.items()}


class _CatalogNormaliser(CdsModelNormaliser):
    """Use typed source graph edges before percent decoding can erase commas.

    The shared normaliser still provides phase/exception interpretation. This
    catalog adapter keeps each parsed Parent identity intact and supplies GTF's
    implicit transcript-to-gene edges without altering source attributes.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.models.clear()
        self.exons.clear()
        self.nodes.clear()
        self.parent_ids = defaultdict(set)
        for row in self.features:
            parsed = parse_gff_attributes(row["source_line"].split("\t", 8)[8])
            tid, gid = row["attributes"].get("transcript_id"), row["attributes"].get("gene_id")
            identifier = row["attributes"].get("ID")
            if not identifier and row["feature"] in {"mRNA", "transcript", "pseudogenic_transcript"}:
                identifier = tid
            if not identifier and row["feature"].lower().endswith("gene"):
                identifier = gid
            if identifier:
                self.nodes[identifier].append(row)
                self.parent_ids[identifier].update(parsed.get("Parent", ()))
            if tid and gid:
                self.parent_ids[tid].add(gid)
            if row["feature"] in {"CDS", "exon"}:
                owners = parsed.get("Parent", parsed.get("transcript_id", parsed.get("ID", ())))
                for owner in owners:
                    (self.models if row["feature"] == "CDS" else self.exons)[owner].append(row)
        for blocks in self.models.values():
            blocks.sort(key=lambda row: row["start"], reverse=blocks[0]["strand"] == "-")

    def ancestry(self, identifiers):
        found, pending = set(), list(identifiers)
        while pending:
            identifier = pending.pop()
            if identifier not in found:
                found.add(identifier)
                pending.extend(self.parent_ids.get(identifier, ()))
        return found


def _junctions(blocks, strand):
    cumulative = 0
    junctions = []
    for index, block in enumerate(blocks[:-1]):
        cumulative += block[1] - block[0]
        following = blocks[index + 1]
        donor, acceptor = (block[1], following[0]) if strand == "+" else (block[0], following[1])
        junctions.append([donor, acceptor, (cumulative - int(blocks[0][2])) % 3])
    return junctions


def build_catalog(species, cds_path, gff_path, genome_path, genetic_code=1):
    """Create a complete original-annotation catalog before representative choice."""
    species = str(species)
    if not species or sanitize_identifier(species) != species:
        raise ValueError("Species must be a nonempty formatted species identifier")
    paths = {"cds": Path(cds_path), "gff": Path(gff_path), "genome": Path(genome_path)}
    frozen = {name: _stat_identity(path) for name, path in paths.items()}
    sources = {name: source_signature(path) for name, path in paths.items()}
    validate_gff(paths["gff"])
    task = {"provider": "direct", "species_prefix": species, "species_key": species,
            "cds_path": paths["cds"], "gff_path": paths["gff"], "genome_path": paths["genome"],
            "gene_grouping_mode": "strict"}
    grouping = build_gff_cds_grouping_index(task)
    normaliser = _CatalogNormaliser({"species": species, "gff": str(paths["gff"]),
                                    "genome": str(paths["genome"]), "genetic_code": int(genetic_code)},
                                   tempfile.gettempdir(), "catalog", attribute_parser=_attribute_parser)
    # Formatter phase evidence comes only from nuclear annotations. Retain
    # organelle provenance here without letting their separate codes/phases vote.
    normaliser.votes = Counter(convention for rows in normaliser.models.values()
                               if not any(row["seqid"] in grouping["organelle_seqids"] for row in rows)
                               and all(row["source"] != "genegalleon_rescue" for row in rows)
                               for conventions in [phase_conventions(rows)] if len(conventions) == 1
                               for convention in conventions)
    normaliser.global_convention = next(iter(normaliser.votes)) if len(normaliser.votes) == 1 else None
    loci, candidates, transcript_aliases, excluded_loci = {}, {}, defaultdict(set), {}
    genes_by_id, transcripts_by_id = {}, {}
    try:
        with indexed_genome(paths["genome"]) as genome:
            for transcript, rows in sorted(normaliser.models.items()):
                organelle_rows = [row for row in rows if row["seqid"] in grouping["organelle_seqids"]]
                if organelle_rows:
                    if len(organelle_rows) != len(rows):
                        raise ValueError("Coding transcript mixes nuclear and organelle contigs: " + transcript)
                    owners = normaliser.ancestry([transcript])
                    source_genes = {row["attributes"]["ID"] for owner in owners
                                    for row in normaliser.nodes.get(owner, [])
                                    if row["feature"].lower().endswith("gene") and row["attributes"].get("ID")}
                    gene = next(iter(source_genes)) if len(source_genes) == 1 else transcript
                    gene_id, candidate_id = _formatted_id(species, gene), _formatted_id(species, transcript)
                    if candidate_id in transcripts_by_id and transcripts_by_id[candidate_id] != transcript:
                        raise ValueError("Distinct GFF transcripts collide after identifier sanitization: " + candidate_id)
                    transcripts_by_id[candidate_id] = transcript
                    blocks = [[row["start"], row["end"], int(row["phase"]) if row["phase"] in {"0", "1", "2"} else -1]
                              for row in rows]
                    seqid, strand = rows[0]["seqid"], rows[0]["strand"]
                    try:
                        sequence = reconstruct_sequence(blocks, seqid, strand, genome)
                        problem = ""
                    except ValueError as error:
                        sequence, problem = "", str(error)
                    candidate = {"candidate_id": candidate_id, "species": species, "gene_id": gene_id,
                                 "source_transcript_id": transcript, "source_gene_id": gene,
                                 "seqid": seqid, "strand": strand, "blocks": blocks, "cds": sequence,
                                 "protein": "", "origin": "original", "junctions": [],
                                 "quality": {"usable": False, "valid_orf": False, "genetic_code": None,
                                             "exclusion_reason": "organelle_annotation", "structure_problem": problem},
                                 "coding_key": hashlib.sha256(json.dumps([sources["genome"]["sha256"], seqid, strand,
                                                                          blocks, sequence]).encode()).hexdigest()}
                    excluded = excluded_loci.setdefault(gene_id, {"species": species, "gene_id": gene_id,
                                                                  "source_gene_id": gene, "seqid": seqid, "strand": strand,
                                                                  "exclusion_reason": "organelle_annotation", "candidates": []})
                    excluded["candidates"].append(candidate)
                    continue
                if transcript in grouping["ambiguous_transcript_gene_tokens"]:
                    raise ValueError("Ambiguous GFF locus ownership: " + transcript)
                gene = grouping["transcript_gene_tokens"].get(transcript, "")
                if not gene:
                    raise ValueError("GFF coding transcript has no unambiguous locus: " + transcript)
                gene_id, candidate_id = _formatted_id(species, gene), _formatted_id(species, transcript)
                if gene_id in genes_by_id and genes_by_id[gene_id] != gene:
                    raise ValueError("Distinct GFF genes collide after identifier sanitization: " + gene_id)
                if candidate_id in transcripts_by_id and transcripts_by_id[candidate_id] != transcript:
                    raise ValueError("Distinct GFF transcripts collide after identifier sanitization: " + candidate_id)
                genes_by_id[gene_id], transcripts_by_id[candidate_id] = gene, transcript
                coordinates = {(row["seqid"], row["strand"]) for row in rows}
                seqid, strand = rows[0]["seqid"], rows[0]["strand"]
                blocks = [[row["start"], row["end"], int(row["phase"]) if row["phase"] in {"0", "1", "2"} else -1]
                          for row in rows]
                ascending = sorted(blocks)
                structure = ("mixed_contig_or_strand" if len(coordinates) != 1 else
                             "overlapping_cds" if any(a[1] > b[0] for a, b in zip(ascending, ascending[1:], strict=False))
                             else "")
                owners = normaliser.ancestry([transcript])
                owned_rows = [*rows, *(row for owner in owners for row in normaliser.nodes.get(owner, []))]
                source_gene_ids = {row["attributes"]["ID"] for row in owned_rows
                                   if row["feature"].lower().endswith("gene") and row["attributes"].get("ID")}
                if len(source_gene_ids) > 1:
                    raise ValueError("Coding transcript has multiple source GFF gene owners: " + transcript)
                source_gene = next(iter(source_gene_ids), gene)
                if not structure and any(row["seqid"] != seqid or
                                         (row["strand"] in {"+", "-"} and row["strand"] != strand)
                                         for row in owned_rows):
                    structure = "ancestor_contig_or_strand_mismatch"
                partial = any(row["attributes"].get("partial", "").lower() in {"true", "1", "yes"}
                              or any(name in row["attributes"] for name in ("start_range", "end_range"))
                              for row in owned_rows)
                exception = normaliser.exceptional(transcript) or (
                    "annotated_pseudogene" if any(row["feature"].lower().startswith("pseudogenic_")
                                                  for row in owned_rows) else "")
                quality = {"structure_problem": structure, "annotated_partial": partial,
                           "annotated_exception": exception,
                           "sequence_mismatch": False}
                try:
                    sequence = reconstruct_sequence(blocks, seqid, strand, genome) if not structure else ""
                except ValueError as error:
                    sequence, quality["structure_problem"] = "", str(error)
                candidate = {"candidate_id": candidate_id, "species": species, "gene_id": gene_id,
                             "source_transcript_id": transcript, "source_gene_id": source_gene, "gene_token": gene,
                             "seqid": seqid, "strand": strand, "cds": sequence, "protein": "",
                             "blocks": blocks, "junctions": [], "quality": quality, "origin": "original",
                             "source_fasta_ids": [], "source_cds_sha256": [], "source_cds": []}
                candidate["quality"] = validate_candidate(candidate, genetic_code)
                if not candidate["quality"]["phase_conflict"] and not candidate["quality"]["phase_unresolved"]:
                    candidate["junctions"] = _junctions(blocks, strand)
                coding = [sources["genome"]["sha256"], seqid, strand, blocks, hashlib.sha256(sequence.encode()).hexdigest()]
                candidate["coding_key"] = hashlib.sha256(json.dumps(coding, separators=(",", ":")).encode()).hexdigest()
                locus = loci.setdefault(gene_id, {"species": species, "gene_id": gene_id,
                                                  "source_gene_id": source_gene, "gene_token": gene,
                                                  "seqid": seqid, "strand": strand,
                                                  "candidates": [], "source_baseline_candidate_id": ""})
                if (locus["seqid"], locus["strand"]) != (seqid, strand):
                    locus["ambiguous_coordinates"] = True
                locus["candidates"].append(candidate)
                candidates[candidate_id] = candidate
                for alias in (transcript, candidate_id):
                    transcript_aliases[alias].add(candidate_id)
                for row in [*rows, *normaliser.nodes.get(transcript, [])]:
                    for name in ("protein_id", "transcript_id", "orig_protein_id", "orig_transcript_id", "Alias", "Name", "Accession"):
                        for alias in filter(None, row["attributes"].get(name, "").split(",")):
                            transcript_aliases[alias].add(candidate_id)
    finally:
        normaliser.close()
    mapping, seen_fasta, bound_by_locus = [], set(), defaultdict(list)
    for identifier, header, raw in fasta_records(paths["cds"]):
        if identifier in seen_fasta:
            raise ValueError("Duplicate source FASTA identifier: " + identifier)
        seen_fasta.add(identifier)
        sequence = "".join(raw.split()).upper()
        primary, _gene_aliases = extract_cds_header_alias_tiers(task, header)
        if (fasta_header_is_organelle(header, grouping["organelle_seqids"], grouping["organelle_aliases"])
                or set(primary).intersection(grouping["organelle_aliases"])):
            mapping.append({"source_fasta_id": identifier,
                            "source_cds_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                            "mapping_status": "excluded_organelle", "exclusion_reason": "organelle_annotation",
                            "candidate_ids": [], "possible_candidate_ids": [], "sequence_agreement": None})
            continue
        explicit = set().union(*(transcript_aliases.get(alias, set()) for alias in primary))
        match = resolve_cds_header_gff_gene(task, header, grouping_index=grouping)
        possible = explicit
        if not possible and match["status"] == "mapped":
            locus = loci.get(_formatted_id(species, match["gene_token"]))
            possible = {candidate["candidate_id"] for candidate in locus["candidates"]} if locus else set()
        exact = {identifier for identifier in possible if _source_matches(sequence, candidates[identifier])}
        chosen = explicit if len(explicit) == 1 else exact
        status = "mapped" if len(chosen) == 1 else "ambiguous" if chosen or possible else "unmapped"
        bound = sorted(chosen) if status == "mapped" else []
        row = {"source_fasta_id": identifier, "source_cds_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
               "mapping_status": status, "candidate_ids": bound, "possible_candidate_ids": sorted(possible),
               "sequence_agreement": bool(bound and bound[0] in exact),
               "source_convention": _source_convention(sequence, candidates[bound[0]]) if bound else ""}
        mapping.append(row)
        if bound:
            candidate = candidates[bound[0]]
            candidate["source_fasta_ids"].append(identifier)
            candidate["source_cds_sha256"].append(row["source_cds_sha256"])
            # Keep supplied sequence/header evidence separate from genomic
            # reconstruction, including unresolved conflicts and source case.
            candidate["source_cds"].append({"source_fasta_id": identifier, "header": header,
                                             "cds": "".join(raw.split()),
                                             "sha256": hashlib.sha256("".join(raw.split()).encode()).hexdigest(),
                                             "normalized_sha256": row["source_cds_sha256"],
                                             "sequence_agreement": row["sequence_agreement"],
                                             "source_convention": row["source_convention"]})
            candidate["quality"]["source_sequence_agreement"] = row["sequence_agreement"]
            candidate.setdefault("source_conventions", []).append(row["source_convention"])
            corrected_length = _corrected_cds_length(sequence, candidate, normaliser.global_convention)
            if corrected_length is not None:
                candidate["corrected_cds_length"] = corrected_length
            if not row["sequence_agreement"]:
                candidate["quality"]["sequence_mismatch"] = True
                candidate["quality"] = validate_candidate(candidate, genetic_code)
            else:
                _infer_missing_phases(candidate, genetic_code)
            bound_by_locus[candidate["gene_id"]].append(bound[0])
    for gene_id, bound in bound_by_locus.items():
        if len(bound) == 1:
            loci[gene_id]["source_baseline_candidate_id"] = bound[0]
    if any(_stat_identity(path) != frozen[name] for name, path in paths.items()):
        raise OSError("Sources changed during coding isoform catalog construction")
    summary = {"loci": len(loci), "candidates": len(candidates),
               "usable_candidates": sum(candidate["quality"]["usable"] for candidate in candidates.values()),
               "coding_paths": len({candidate["coding_key"] for candidate in candidates.values()}),
               "fasta_mapping": dict(Counter(row["mapping_status"] for row in mapping))}
    if excluded_loci:
        summary.update(excluded_loci=len(excluded_loci),
                       excluded_candidates=sum(len(locus["candidates"]) for locus in excluded_loci.values()))
    return {"schema_version": SCHEMA_VERSION, "species": species, "genetic_code": int(genetic_code),
            "sources": sources, "loci": [loci[key] for key in sorted(loci)], "fasta_mapping": mapping,
            "excluded_loci": [excluded_loci[key] for key in sorted(excluded_loci)],
            "summary": summary}


def write_catalog(catalog, output_dir):
    """Publish catalog and complete transcript-keyed CDS/protein FASTAs atomically."""
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    paths = {"catalog": directory / "catalog.json", "cds": directory / "all_isoforms.cds.fa",
             "protein": directory / "all_isoforms.protein.fa", "loci": directory / "loci.jsonl",
             "metadata": directory / "catalog_metadata.json", "excluded": directory / "excluded_loci.jsonl"}
    with tempfile.TemporaryDirectory(prefix=".catalog-", dir=directory) as scratch:
        temporary = {name: Path(scratch) / path.name for name, path in paths.items()}
        temporary["catalog"].write_text(json.dumps(catalog, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        metadata = {key: value for key, value in catalog.items() if key not in {"loci", "excluded_loci"}}
        excluded_loci = catalog.get("excluded_loci", [])
        metadata["excluded_loci_summary"] = {"loci": len(excluded_loci),
                                             "candidates": sum(len(locus["candidates"]) for locus in excluded_loci),
                                             "reasons": dict(Counter(locus["exclusion_reason"] for locus in excluded_loci))}
        temporary["metadata"].write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        with temporary["loci"].open("w", encoding="utf-8") as handle:
            for locus in catalog["loci"]:
                handle.write(json.dumps(locus, sort_keys=True, separators=(",", ":")) + "\n")
        with temporary["excluded"].open("w", encoding="utf-8") as handle:
            for locus in excluded_loci:
                handle.write(json.dumps(locus, sort_keys=True, separators=(",", ":")) + "\n")
        with temporary["cds"].open("w", encoding="utf-8") as cds, temporary["protein"].open("w", encoding="utf-8") as protein:
            for locus in catalog["loci"]:
                for candidate in locus["candidates"]:
                    cds.write(">" + candidate["candidate_id"] + "\n" + candidate["cds"] + "\n")
                    if candidate["quality"]["usable"]:
                        protein.write(">" + candidate["candidate_id"] + "\n" + candidate["protein"] + "\n")
        for name, destination in paths.items():
            os.replace(temporary[name], destination)
    return {name: str(path.resolve()) for name, path in paths.items()}
