"""Annotation implementation: genbank."""

import re
import sys
from collections import defaultdict
from urllib.parse import quote

from format_species_writers import apply_common_replacements, open_text

try:
    from Bio import SeqIO
except Exception:  # pragma: no cover - runtime without biopython
    SeqIO = None



from .common import (
    build_gff_genome_seqid_map,
    choose_first_gff_attribute,
    collapse_transcript_suffix,
    first_token,
    load_genome_sequences,
    parse_gff_attributes,
    resolve_feature_gene_token,
    reverse_complement,
    sanitize_identifier,
    transcript_feature_gene_token,
)
from .grouping import (
    build_rescued_gene_tokens_for_transcripts,
)


def iter_genome_records_from_gbff(path):
    for record in iter_genbank_records(path):
        record_id = extract_seqrecord_id(record)
        sequence = re.sub(r"\s+", "", str(record.seq or "")).upper()
        if sequence == "":
            continue
        yield record_id, sequence


def require_biopython_for_gbff(path):
    if SeqIO is None:
        raise RuntimeError("Biopython is required to process GenBank/GBFF/EMBL inputs: {}".format(path))


def insdc_flatfile_format(path):
    name = str(getattr(path, "name", path) or "").lower()
    if name.endswith(".gz"):
        name = name[: -len(".gz")]
    return "embl" if name.endswith(".embl") else "genbank"


def iter_genbank_records(path):
    require_biopython_for_gbff(path)
    with open_text(path, "rt") as handle:
        for record in SeqIO.parse(handle, insdc_flatfile_format(path)):
            yield record


def choose_first_feature_qualifier(feature, keys):
    qualifiers = getattr(feature, "qualifiers", {}) or {}
    for key in keys:
        values = qualifiers.get(key, ())
        if isinstance(values, str):
            values = (values,)
        for value in values:
            text = str(value or "").strip()
            if text != "":
                return text
    return ""


def extract_genbank_gene_id_from_dbxref(feature):
    qualifiers = getattr(feature, "qualifiers", {}) or {}
    values = qualifiers.get("db_xref", ())
    if isinstance(values, str):
        values = (values,)
    for raw in values:
        text = str(raw or "").strip()
        if text.startswith("GeneID:"):
            candidate = text.split(":", 1)[1].strip()
            if candidate != "":
                return "GeneID{}".format(candidate)
    return ""


def normalize_genbank_feature_location(feature):
    location = getattr(feature, "location", None)
    if location is None:
        return []
    if hasattr(location, "parts") and len(getattr(location, "parts", ())) > 0:
        parts = list(location.parts)
    else:
        parts = [location]
    normalized = []
    for part in parts:
        try:
            start = int(part.start) + 1
            end = int(part.end)
        except Exception:
            continue
        if start > end:
            start, end = end, start
        normalized.append((start, end))
    return normalized


def get_genbank_feature_strand(feature):
    strand = getattr(getattr(feature, "location", None), "strand", None)
    return "-" if strand == -1 else "+"


def infer_genbank_gene_token(feature, fallback_prefix, fallback_index):
    gene_token = choose_first_feature_qualifier(
        feature,
        ("locus_tag", "gene", "gene_synonym", "Name", "ID"),
    )
    if gene_token == "":
        gene_token = extract_genbank_gene_id_from_dbxref(feature)
    if gene_token == "":
        gene_token = "{}_gene{}".format(fallback_prefix, fallback_index)
    return sanitize_identifier(gene_token)


def infer_genbank_transcript_token(feature, gene_token, fallback_index):
    transcript_token = choose_first_feature_qualifier(
        feature,
        ("protein_id", "transcript_id", "mrna", "rna_id"),
    )
    if transcript_token == "":
        transcript_token = "{}.t{}".format(gene_token, fallback_index)
    return sanitize_identifier(transcript_token)


def parse_genbank_codon_start(feature):
    raw = choose_first_feature_qualifier(feature, ("codon_start",))
    if raw == "":
        return 1
    try:
        value = int(raw)
    except Exception:
        return 1
    if value not in (1, 2, 3):
        return 1
    return value


def extract_seqrecord_id(record):
    candidate = first_token(apply_common_replacements(str(getattr(record, "id", "") or "").strip()))
    if candidate != "":
        return candidate
    candidate = first_token(apply_common_replacements(str(getattr(record, "name", "") or "").strip()))
    if candidate != "":
        return candidate
    return "unnamed_record"


def format_gff_attributes(attributes):
    fields = []
    for key, value in attributes:
        text = str(value or "").strip()
        if text == "":
            continue
        fields.append("{}={}".format(key, quote(text, safe="._:-|")))
    return ";".join(fields) if len(fields) > 0 else "."


def compute_gff_phases(parts, strand, codon_start):
    if len(parts) == 0:
        return []
    transcript_order = list(parts)
    if strand == "-":
        transcript_order = list(reversed(transcript_order))
    phases_by_part = {}
    effective_bases = 0
    for index, part in enumerate(transcript_order):
        length = part["end"] - part["start"] + 1
        phase = codon_start - 1 if index == 0 else (3 - (effective_bases % 3)) % 3
        phases_by_part[(part["start"], part["end"])] = str(phase)
        effective_bases += max(0, length - phase)
    return [phases_by_part[(part["start"], part["end"])] for part in parts]


def iter_gbff_coding_entries(path):
    gene_entries = {}
    transcripts_by_gene = defaultdict(list)
    fallback_gene_index = 0
    transcript_counts_by_gene = defaultdict(int)

    for record in iter_genbank_records(path):
        seqid = extract_seqrecord_id(record)
        for feature in getattr(record, "features", ()):
            feature_type = str(getattr(feature, "type", "") or "").lower()
            parts = normalize_genbank_feature_location(feature)
            if len(parts) == 0:
                continue
            strand = get_genbank_feature_strand(feature)
            if feature_type == "gene":
                fallback_gene_index += 1
                gene_token = infer_genbank_gene_token(feature, seqid, fallback_gene_index)
                gene_key = (seqid, gene_token)
                gene_name = choose_first_feature_qualifier(feature, ("gene", "locus_tag"))
                start = min(start for start, _end in parts)
                end = max(end for _start, end in parts)
                existing = gene_entries.get(gene_key)
                if existing is None:
                    gene_entries[gene_key] = {
                        "seqid": seqid,
                        "gene_token": gene_token,
                        "gene_name": gene_name,
                        "start": start,
                        "end": end,
                        "strand": strand,
                    }
                else:
                    existing["start"] = min(existing["start"], start)
                    existing["end"] = max(existing["end"], end)
                    if existing["gene_name"] == "" and gene_name != "":
                        existing["gene_name"] = gene_name
                continue
            if feature_type != "cds":
                continue
            fallback_gene_index += 1
            gene_token = infer_genbank_gene_token(feature, seqid, fallback_gene_index)
            gene_key = (seqid, gene_token)
            gene_name = choose_first_feature_qualifier(feature, ("gene", "locus_tag"))
            transcript_counts_by_gene[gene_key] += 1
            transcript_token = infer_genbank_transcript_token(
                feature,
                gene_token,
                transcript_counts_by_gene[gene_key],
            )
            sequence = str(feature.extract(record.seq)).upper()
            codon_start = parse_genbank_codon_start(feature)
            if codon_start > 1:
                sequence = sequence[codon_start - 1 :]
            if sequence == "":
                continue
            start = min(start for start, _end in parts)
            end = max(end for _start, end in parts)
            existing = gene_entries.get(gene_key)
            if existing is None:
                gene_entries[gene_key] = {
                    "seqid": seqid,
                    "gene_token": gene_token,
                    "gene_name": gene_name,
                    "start": start,
                    "end": end,
                    "strand": strand,
                }
            else:
                existing["start"] = min(existing["start"], start)
                existing["end"] = max(existing["end"], end)
                if existing["gene_name"] == "" and gene_name != "":
                    existing["gene_name"] = gene_name
            transcript_entry = {
                "seqid": seqid,
                "gene_key": gene_key,
                "gene_token": gene_token,
                "gene_name": gene_name,
                "transcript_token": transcript_token,
                "start": start,
                "end": end,
                "strand": strand,
                "parts": [{"start": part_start, "end": part_end} for part_start, part_end in parts],
                "codon_start": codon_start,
                "sequence": sequence,
            }
            transcripts_by_gene[gene_key].append(transcript_entry)

    for gene_key, gene_entry in gene_entries.items():
        transcripts = transcripts_by_gene.get(gene_key, [])
        if len(transcripts) == 0:
            continue
        yield (
            gene_entry,
            sorted(
                transcripts,
                key=lambda item: (item["start"], item["end"], item["transcript_token"]),
            ),
        )


def derive_cds_records_from_gff_and_genome(task):
    gff_path = task.get("gff_path")
    genome_path = task.get("genome_path")
    if gff_path is None or genome_path is None:
        raise ValueError("GFF+genome inputs are required to derive CDS for {}".format(task.get("species_key", "")))

    genome_sequences = load_genome_sequences(genome_path)
    feature_records = {}
    cds_features_by_transcript = defaultdict(list)
    utr_features_by_transcript = defaultdict(list)
    gene_cache = {}

    with open_text(gff_path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if line == "":
                continue
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, feature_type, start_text, end_text, _score, strand, phase_text, attr_text = parts[:9]
            feature_type_lower = str(feature_type or "").strip().lower()
            attrs = parse_gff_attributes(attr_text)
            feature_id = choose_first_gff_attribute(attrs, ("ID", "transcript_id", "protein_id", "Name"))
            parents = attrs.get("Parent", ())
            if feature_id != "":
                feature_records[feature_id] = {
                    "feature_type": feature_type_lower,
                    "parents": parents,
                    "attrs": attrs,
                }
            if feature_type_lower in ("five_prime_utr", "three_prime_utr"):
                try:
                    start = int(start_text)
                    end = int(end_text)
                except Exception:
                    continue
                if start > end:
                    start, end = end, start
                for parent_id in parents:
                    if str(parent_id).strip() == "":
                        continue
                    utr_features_by_transcript[parent_id].append(
                        {
                            "seqid": str(seqid or "").strip(),
                            "start": start,
                            "end": end,
                            "strand": str(strand or "").strip() or "+",
                        }
                    )
                continue
            if feature_type_lower != "cds":
                continue
            try:
                start = int(start_text)
                end = int(end_text)
            except Exception:
                continue
            if start > end:
                start, end = end, start
            transcript_ids = [value for value in parents if str(value).strip() != ""]
            if len(transcript_ids) == 0:
                fallback_id = feature_id
                if fallback_id == "":
                    fallback_id = choose_first_gff_attribute(attrs, ("transcript_id", "protein_id", "Name"))
                if fallback_id == "":
                    fallback_id = "{}:{}-{}".format(seqid, start, end)
                transcript_ids = [fallback_id]
            explicit_gene = choose_first_gff_attribute(
                attrs,
                ("gene", "gene_id", "locus_tag", "geneName", "Parent_Accession", "Accession"),
            )
            for transcript_id in transcript_ids:
                gene_token = explicit_gene
                if gene_token == "":
                    gene_token = resolve_feature_gene_token(
                        transcript_id,
                        feature_records,
                        task["provider"],
                        gene_cache,
                        set(),
                    )
                if gene_token == "":
                    gene_token = collapse_transcript_suffix(task["provider"], transcript_id)
                cds_features_by_transcript[transcript_id].append(
                    {
                        "seqid": str(seqid or "").strip(),
                        "start": start,
                        "end": end,
                        "strand": str(strand or "").strip() or "+",
                        "gene_token": gene_token,
                    }
                )

    rescued_gene_tokens = build_rescued_gene_tokens_for_transcripts(task, cds_features_by_transcript)
    required_gff_seqids = {
        feature["seqid"]
        for features in cds_features_by_transcript.values()
        for feature in features
        if str(feature.get("seqid", "") or "").strip() != ""
    }
    try:
        genome_seqid_map, missing_seqids = build_gff_genome_seqid_map(
            genome_sequences,
            required_gff_seqids,
        )
    except ValueError as exc:
        raise ValueError(
            "Genome FASTA for {} has ambiguous sequence-ID aliases required by {}: {}".format(
                task.get("species_key", ""),
                gff_path,
                exc,
            )
        ) from exc
    if len(missing_seqids) > 0:
        raise ValueError(
            "Genome FASTA for {} is missing sequence(s) {} required by {}".format(
                task.get("species_key", ""),
                ",".join(missing_seqids[:5]),
                gff_path,
            )
        )
    for transcript_id in sorted(cds_features_by_transcript.keys()):
        features = cds_features_by_transcript[transcript_id]
        if len(features) == 0:
            continue
        strands = sorted({feature["strand"] for feature in features})
        if len(strands) > 1:
            sys.stderr.write(
                "Warning: skipping transcript '{}' in {} because CDS features use mixed strands: {}\n".format(
                    transcript_id,
                    gff_path,
                    ",".join(strands),
                )
            )
            continue
        utr_features = utr_features_by_transcript.get(transcript_id, ())
        strand = strands[0]
        transcript_gene_token = rescued_gene_tokens.get(transcript_id, "") or transcript_feature_gene_token(features)
        trimmed_features = []
        for feature in features:
            segments = [(feature["start"], feature["end"])]
            for utr in utr_features:
                if (
                    utr["seqid"] != feature["seqid"]
                    or utr["strand"] != feature["strand"]
                    or utr["end"] < feature["start"]
                    or feature["end"] < utr["start"]
                ):
                    continue
                next_segments = []
                for seg_start, seg_end in segments:
                    if utr["end"] < seg_start or seg_end < utr["start"]:
                        next_segments.append((seg_start, seg_end))
                        continue
                    if seg_start < utr["start"]:
                        next_segments.append((seg_start, utr["start"] - 1))
                    if utr["end"] < seg_end:
                        next_segments.append((utr["end"] + 1, seg_end))
                segments = [(seg_start, seg_end) for seg_start, seg_end in next_segments if seg_start <= seg_end]
                if len(segments) == 0:
                    break
            for seg_start, seg_end in segments:
                trimmed_features.append(
                    {
                        "seqid": feature["seqid"],
                        "start": seg_start,
                        "end": seg_end,
                        "strand": feature["strand"],
                        "gene_token": transcript_gene_token,
                    }
                )
        if len(trimmed_features) == 0:
            continue
        trimmed_strands = sorted({feature["strand"] for feature in trimmed_features})
        if len(trimmed_strands) > 1:
            sys.stderr.write(
                "Warning: skipping transcript '{}' in {} because CDS features use mixed strands: {}\n".format(
                    transcript_id,
                    gff_path,
                    ",".join(trimmed_strands),
                )
            )
            continue
        ordered = sorted(trimmed_features, key=lambda item: item["start"], reverse=(strand == "-"))
        pieces = []
        gene_token = ""
        for feature in ordered:
            if feature["strand"] != strand:
                raise ValueError("Mixed-strand CDS features for transcript '{}' in {}".format(transcript_id, gff_path))
            seqid = feature["seqid"]
            fasta_seqid = genome_seqid_map.get(seqid, seqid)
            genome_seq = genome_sequences.get(fasta_seqid, "")
            if genome_seq == "":
                raise ValueError(
                    "Genome FASTA for {} is missing sequence '{}' required by {}".format(
                        task.get("species_key", ""),
                        seqid,
                        gff_path,
                    )
                )
            start_idx = feature["start"] - 1
            end_idx = feature["end"]
            piece = genome_seq[start_idx:end_idx]
            # GFF3 phase annotates codon frame continuity; it does not indicate
            # bases to trim from the CDS feature when reconstructing the spliced
            # CDS sequence. We therefore concatenate the full CDS spans in
            # transcript order and leave any frame repair to the final padding
            # step used elsewhere in the pipeline.
            if strand == "-":
                piece = reverse_complement(piece)
            if piece != "":
                pieces.append(piece)
            if gene_token == "":
                gene_token = feature["gene_token"]
        sequence = "".join(pieces)
        if sequence == "":
            continue
        header = str(transcript_id or "").strip()
        if gene_token != "":
            header = "{} [gene={}]".format(header, gene_token)
        yield header, sequence


def derive_cds_records_from_gbff(task):
    gbff_path = task.get("gbff_path")
    if gbff_path is None:
        raise ValueError("GBFF input is required to derive CDS for {}".format(task.get("species_key", "")))
    for _gene_entry, transcripts in iter_gbff_coding_entries(gbff_path):
        for transcript in transcripts:
            header = transcript["transcript_token"]
            if transcript["gene_token"] != "":
                header = "{} [gene={}]".format(header, transcript["gene_token"])
            yield header, transcript["sequence"]


def iter_gff_lines_from_gbff(task):
    gbff_path = task.get("gbff_path")
    if gbff_path is None:
        raise ValueError("GBFF input is required to derive GFF for {}".format(task.get("species_key", "")))
    yield "##gff-version 3\n"
    for gene_entry, transcripts in sorted(
        iter_gbff_coding_entries(gbff_path),
        key=lambda item: (item[0]["seqid"], item[0]["start"], item[0]["gene_token"]),
    ):
        seqid = gene_entry["seqid"]
        strand = gene_entry["strand"]
        gene_id = "gene:{}".format(gene_entry["gene_token"])
        gene_attrs = [("ID", gene_id), ("Name", gene_entry["gene_name"] or gene_entry["gene_token"])]
        yield (
            "\t".join(
                [
                    seqid,
                    "genbank_derived",
                    "gene",
                    str(gene_entry["start"]),
                    str(gene_entry["end"]),
                    ".",
                    strand,
                    ".",
                    format_gff_attributes(gene_attrs),
                ]
            )
            + "\n"
        )
        for transcript in transcripts:
            transcript_id = "transcript:{}".format(transcript["transcript_token"])
            transcript_attrs = [
                ("ID", transcript_id),
                ("Parent", gene_id),
                ("Name", transcript["transcript_token"]),
            ]
            yield (
                "\t".join(
                    [
                        seqid,
                        "genbank_derived",
                        "mRNA",
                        str(transcript["start"]),
                        str(transcript["end"]),
                        ".",
                        transcript["strand"],
                        ".",
                        format_gff_attributes(transcript_attrs),
                    ]
                )
                + "\n"
            )
            ordered_parts = sorted(transcript["parts"], key=lambda item: item["start"])
            phases = compute_gff_phases(ordered_parts, transcript["strand"], transcript["codon_start"])
            for index, part in enumerate(ordered_parts, start=1):
                cds_attrs = [
                    ("ID", "cds:{}.{}".format(transcript["transcript_token"], index)),
                    ("Parent", transcript_id),
                ]
                yield (
                    "\t".join(
                        [
                            seqid,
                            "genbank_derived",
                            "CDS",
                            str(part["start"]),
                            str(part["end"]),
                            ".",
                            transcript["strand"],
                            phases[index - 1],
                            format_gff_attributes(cds_attrs),
                        ]
                    )
                    + "\n"
                )
