"""Annotation implementation: common."""

import io
import re
import tarfile
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

from format_species_provider_config import ENSEMBL_LIKE_PROVIDERS, ORYZA_MINUTA_PROVIDER
from format_species_writers import apply_common_replacements, open_text

try:
    from Bio import SeqIO
except Exception:  # pragma: no cover - runtime without biopython
    SeqIO = None

from format_species_common import (
    is_fasta_filename,
)
from format_species_constants import (
    ENSEMBL_GENE_ID_PATTERN,
    FASTA_ARCHIVE_EXTENSIONS,
    FASTA_EXTENSIONS,
    GENE_GROUPING_MODES,
    GFF_EXTENSIONS,
    INVALID_ID_CHARS,
    SPECIES_KEY_FILENAME_SUFFIXES,
    TAXONOMIC_INFRASPECIFIC_RANK_ALIASES,
    TAXONOMIC_PROXIMITY_QUALIFIERS,
)
from format_species_taxonomy import (
    is_hybrid_connector_token,
    species_prefix_token_count,
    taxonomic_token_key,
)


def task_has_usable_source_bundle(cds_path, gff_path, gbff_path, genome_path):
    return cds_path is not None or gbff_path is not None or (gff_path is not None and genome_path is not None)


def gff_seqids(path):
    seqids = set()
    with open_text(path, "rt") as handle:
        for raw_line in handle:
            if raw_line.startswith("#") or raw_line.strip() == "":
                continue
            parts = raw_line.rstrip("\n\r").split("\t")
            if len(parts) < 9:
                continue
            seqid = str(parts[0] or "").strip()
            if seqid != "":
                seqids.add(seqid)
    return seqids


def fasta_headers_overlap_gff_seqids(fasta_path, gff_path):
    seqids = gff_seqids(gff_path)
    if len(seqids) == 0:
        return False
    for header, _sequence in iter_fasta_records(fasta_path):
        if first_token(header) in seqids:
            return True
    return False


def fasta_looks_like_gff_genome(fasta_path, gff_paths):
    for gff_path in gff_paths:
        try:
            if fasta_headers_overlap_gff_seqids(fasta_path, gff_path):
                return True
        except Exception:
            continue
    return False


def first_token(text):
    parts = text.split()
    if len(parts) == 0:
        return ""
    return parts[0]


def first_fasta_header_id(path):
    with open_text(path, "rt") as handle:
        for raw_line in handle:
            if not raw_line.startswith(">"):
                continue
            return first_token(raw_line[1:].strip())
    return ""


def gff_has_seqid(path, expected_seqid):
    target = str(expected_seqid or "").strip()
    if target == "":
        return False
    with open_text(path, "rt") as handle:
        for raw_line in handle:
            if raw_line.startswith("#") or raw_line.strip() == "":
                continue
            parts = raw_line.rstrip("\n\r").split("\t", 1)
            if len(parts) >= 1 and str(parts[0] or "").strip() == target:
                return True
    return False


def fasta_gff_lcl_prefix_mismatch(fasta_path, gff_path):
    fasta_seqid = first_fasta_header_id(fasta_path)
    if not fasta_seqid.startswith("lcl|"):
        return None
    stripped_seqid = fasta_seqid[len("lcl|") :]
    if stripped_seqid == "" or not gff_has_seqid(gff_path, stripped_seqid):
        return None
    return fasta_seqid, stripped_seqid


def extract_header_tag_value(header, tag):
    pattern = r"\[{}=([^\]]+)\]".format(re.escape(tag))
    match = re.search(pattern, header)
    if match is not None:
        return match.group(1).strip()
    plain_pattern = r"(?:^|\s){}=([^\s]+)".format(re.escape(tag))
    match = re.search(plain_pattern, header)
    if match is not None:
        return match.group(1).strip()
    return ""


def extract_ncbi_gene_id_from_header(header):
    db_xref = extract_header_tag_value(header, "db_xref")
    if db_xref == "":
        return ""
    for token in db_xref.split(","):
        text = token.strip()
        if text.startswith("GeneID:"):
            return text[len("GeneID:") :]
    return ""


def extract_ncbi_ensembl_gene_id_from_header(header):
    db_xref = extract_header_tag_value(header, "db_xref")
    if db_xref == "":
        return ""
    for token in db_xref.split(","):
        text = token.strip()
        if not text.startswith("Ensembl:"):
            continue
        candidate = text[len("Ensembl:") :].strip()
        if ENSEMBL_GENE_ID_PATTERN.match(candidate):
            return candidate
    return ""


def strip_species_key_filename_suffixes(value):
    name = Path(str(value or "")).name
    lower = name.lower()
    for suffix in sorted(SPECIES_KEY_FILENAME_SUFFIXES, key=len, reverse=True):
        if lower.endswith(suffix):
            return name[: len(name) - len(suffix)]
    return name


def species_prefix_filename_token(token):
    text = str(token or "").strip()
    if is_hybrid_connector_token(text):
        return "x"
    key = taxonomic_token_key(text)
    if key == "sp" or key in TAXONOMIC_PROXIMITY_QUALIFIERS or key in TAXONOMIC_INFRASPECIFIC_RANK_ALIASES.values():
        return text
    return text.split(".", 1)[0]


def species_prefix_from_value(value):
    name = strip_species_key_filename_suffixes(value)
    tokens = [token for token in name.split("_") if token != ""]
    count = species_prefix_token_count(tokens)
    if count == 0:
        return ""
    prefix_tokens = [species_prefix_filename_token(token) for token in tokens[:count]]
    prefix_tokens = [token for token in prefix_tokens if token != ""]
    return "_".join(prefix_tokens)


def normalize_output_basename(source_name, species_prefix):
    if source_name.startswith(species_prefix + "_"):
        return source_name
    dotted_prefix = species_prefix + "."
    if source_name.startswith(dotted_prefix):
        return species_prefix + "_" + source_name[len(dotted_prefix) :]
    return species_prefix + "_" + source_name


def strip_suffix_case_insensitive(text, suffixes):
    lower = text.lower()
    for suffix in sorted(suffixes, key=len, reverse=True):
        if lower.endswith(suffix):
            return text[: len(text) - len(suffix)]
    return text


def normalize_cds_output_basename(source_name, species_prefix):
    normalized = normalize_output_basename(source_name, species_prefix)
    if normalized.lower().endswith(".cds.gz"):
        return normalized[: -len(".gz")] + ".fa.gz"
    stem = strip_suffix_case_insensitive(
        normalized,
        tuple(ext for ext in FASTA_EXTENSIONS if ext not in (".cds", ".cds.gz")),
    )
    return stem + ".fa.gz"


def normalize_genome_output_basename(source_name, species_prefix):
    normalized = normalize_output_basename(source_name, species_prefix)
    stem = strip_suffix_case_insensitive(normalized, FASTA_EXTENSIONS)
    return stem + ".fa.gz"


def normalize_gff_output_basename(source_name, species_prefix):
    normalized = normalize_output_basename(source_name, species_prefix)
    stem = strip_suffix_case_insensitive(normalized, GFF_EXTENSIONS)
    return stem + ".gff.gz"


def sanitize_identifier(identifier):
    out = identifier.replace("−", "-")
    out = apply_common_replacements(out)
    out = INVALID_ID_CHARS.sub("_", out)
    out = re.sub(r"_+", "_", out)
    out = out.strip("_")
    if out == "":
        out = "unnamed"
    return out


def pad_to_codon_length(seq):
    remainder = len(seq) % 3
    if remainder == 0:
        return seq
    return seq + ("N" * (3 - remainder))


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def normalize_gff_attribute_value(raw_value):
    value = str(raw_value or "").strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in ("'", '"'):
        value = value[1:-1]
    else:
        if value.startswith(("'", '"')):
            value = value[1:]
        if value.endswith(("'", '"')):
            value = value[:-1]
    return unquote(value.strip())


def parse_gff_attributes(attr_text):
    attrs = defaultdict(list)
    for raw_field in str(attr_text or "").split(";"):
        field = raw_field.strip()
        if field == "":
            continue
        equal_pos = field.find("=")
        space_pos = field.find(" ")
        if equal_pos != -1 and (space_pos == -1 or equal_pos < space_pos):
            key = field[:equal_pos].strip()
            raw_value = field[equal_pos + 1 :]
        elif space_pos > 0:
            key = field[:space_pos].strip()
            raw_value = field[space_pos + 1 :]
        else:
            continue
        if key == "":
            continue
        for raw_item in str(raw_value or "").split(","):
            value = normalize_gff_attribute_value(raw_item)
            if value != "":
                attrs[key].append(value)
    return {key: tuple(values) for key, values in attrs.items()}


def choose_first_gff_attribute(attrs, keys):
    for key in keys:
        values = attrs.get(key, ())
        if len(values) > 0 and str(values[0]).strip() != "":
            return str(values[0]).strip()
    return ""


def load_genome_sequences(path):
    sequences = {}
    for header, sequence in iter_fasta_records(path):
        full_name = apply_common_replacements(str(header or "").strip())
        token = first_token(full_name)
        seq = re.sub(r"\s+", "", str(sequence or "")).upper()
        if token != "" and token not in sequences:
            sequences[token] = seq
        if full_name != "" and full_name not in sequences:
            sequences[full_name] = seq
    return sequences


def resolve_feature_gene_token(feature_id, feature_records, provider, cache, active):
    feature_text = str(feature_id or "").strip()
    if feature_text == "":
        return ""
    if feature_text in cache:
        return cache[feature_text]
    if feature_text in active:
        collapsed = collapse_transcript_suffix(provider, feature_text)
        cache[feature_text] = collapsed
        return collapsed

    record = feature_records.get(feature_text)
    if record is None:
        collapsed = collapse_transcript_suffix(provider, feature_text)
        cache[feature_text] = collapsed if collapsed != "" else feature_text
        return cache[feature_text]

    active.add(feature_text)
    attrs = record.get("attrs", {})
    direct = choose_first_gff_attribute(
        attrs,
        ("gene", "gene_id", "locus_tag", "geneName", "Accession", "Parent_Accession"),
    )
    if direct == "":
        if record.get("feature_type") == "gene":
            direct = choose_first_gff_attribute(attrs, ("Name", "ID"))
        else:
            for parent_id in record.get("parents", ()):
                direct = resolve_feature_gene_token(parent_id, feature_records, provider, cache, active)
                if direct != "":
                    break
    if direct == "":
        direct = choose_first_gff_attribute(attrs, ("Name", "ID"))
    if direct == "":
        direct = collapse_transcript_suffix(provider, feature_text)
    active.remove(feature_text)
    cache[feature_text] = str(direct or "").strip()
    return cache[feature_text]


def build_coge_gff_gene_id_map(gff_path):
    feature_parents = {}
    direct_gene_tokens = {}
    alias_keys = ("ID", "Name", "Alias", "CDS", "mRNA", "transcript_id", "protein_id")
    with open_text(gff_path, "rt") as handle:
        for raw_line in handle:
            if raw_line.startswith("#") or raw_line.strip() == "":
                continue
            parts = raw_line.rstrip("\n\r").split("\t")
            if len(parts) < 9:
                continue
            feature_type = str(parts[2] or "").strip().lower()
            if feature_type not in ("gene", "mrna"):
                continue
            attrs = parse_gff_attributes(parts[8])
            feature_id = choose_first_gff_attribute(attrs, ("ID", "transcript_id", "protein_id", "Name"))
            if feature_id == "":
                continue
            feature_parents[feature_id] = tuple(attrs.get("Parent", ()))
            direct = choose_first_gff_attribute(
                attrs,
                ("gene", "gene_id", "locus_tag", "geneName", "Accession", "Parent_Accession"),
            )
            if direct == "" and feature_type == "gene":
                direct = choose_first_gff_attribute(attrs, ("Name", "ID"))
            if direct != "":
                direct_gene_tokens[feature_id] = direct

    cache = {}

    def resolve_gene_token(feature_id, active):
        feature_text = str(feature_id or "").strip()
        if feature_text == "":
            return ""
        if feature_text in cache:
            return cache[feature_text]
        if feature_text in active:
            return collapse_transcript_suffix("coge", feature_text)
        direct = direct_gene_tokens.get(feature_text, "")
        if direct != "":
            cache[feature_text] = direct
            return direct
        active.add(feature_text)
        for parent_id in feature_parents.get(feature_text, ()):
            direct = resolve_gene_token(parent_id, active)
            if direct != "":
                break
        active.remove(feature_text)
        if direct == "":
            direct = collapse_transcript_suffix("coge", feature_text)
        cache[feature_text] = direct
        return direct

    tokens_by_alias = defaultdict(set)
    with open_text(gff_path, "rt") as handle:
        for raw_line in handle:
            if raw_line.startswith("#") or raw_line.strip() == "":
                continue
            parts = raw_line.rstrip("\n\r").split("\t")
            if len(parts) < 9:
                continue
            feature_type = str(parts[2] or "").strip().lower()
            if feature_type not in ("gene", "mrna", "cds"):
                continue
            attrs = parse_gff_attributes(parts[8])
            feature_id = choose_first_gff_attribute(attrs, ("ID", "transcript_id", "protein_id", "Name"))
            direct = choose_first_gff_attribute(
                attrs,
                ("gene", "gene_id", "locus_tag", "geneName", "Accession", "Parent_Accession"),
            )
            if direct == "" and feature_type == "gene":
                direct = choose_first_gff_attribute(attrs, ("Name", "ID"))
            if direct == "" and feature_id in feature_parents:
                direct = resolve_gene_token(feature_id, set())
            if direct == "":
                for parent_id in attrs.get("Parent", ()):
                    direct = resolve_gene_token(parent_id, set())
                    if direct != "":
                        break
            if direct == "" and feature_id != "":
                direct = collapse_transcript_suffix("coge", feature_id)
            if direct == "":
                continue
            if feature_id != "":
                tokens_by_alias[feature_id].add(direct)
            for key in alias_keys:
                for value in attrs.get(key, ()):
                    alias = str(value or "").strip()
                    if alias != "":
                        tokens_by_alias[alias].add(direct)

    return {
        alias: next(iter(gene_tokens))
        for alias, gene_tokens in tokens_by_alias.items()
        if len(gene_tokens) == 1
    }


def merge_coordinate_intervals(intervals):
    ordered = sorted((int(start), int(end)) for start, end in intervals if int(start) <= int(end))
    merged = []
    for start, end in ordered:
        if len(merged) == 0 or start > merged[-1][1] + 1:
            merged.append([start, end])
            continue
        if end > merged[-1][1]:
            merged[-1][1] = end
    return [(start, end) for start, end in merged]


def compute_interval_overlap_bases(left_intervals, right_intervals):
    overlap = 0
    left_index = 0
    right_index = 0
    while left_index < len(left_intervals) and right_index < len(right_intervals):
        left_start, left_end = left_intervals[left_index]
        right_start, right_end = right_intervals[right_index]
        start = max(left_start, right_start)
        end = min(left_end, right_end)
        if start <= end:
            overlap += end - start + 1
        if left_end <= right_end:
            left_index += 1
        else:
            right_index += 1
    return overlap


def transcript_feature_gene_token(features):
    counts = defaultdict(int)
    for feature in features:
        token = str(feature.get("gene_token", "") or "").strip()
        if token != "":
            counts[token] += 1
    if len(counts) == 0:
        return ""
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]


def _iter_fasta_records_from_handle(handle):
    header = None
    chunks = []
    for raw_line in handle:
        line = raw_line.rstrip("\n\r")
        if line == "":
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header = line[1:]
            chunks = []
            continue
        chunks.append(re.sub(r"\s+", "", line))
    if header is not None:
        yield header, "".join(chunks)


def iter_fasta_records(path):
    path_lower = path.name.lower()
    if any(path_lower.endswith(suffix) for suffix in FASTA_ARCHIVE_EXTENSIONS):
        with tarfile.open(path, "r:*") as archive:
            members = [member for member in archive.getmembers() if member.isfile()]
            fasta_members = [member for member in members if is_fasta_filename(Path(member.name).name)]
            if len(fasta_members) == 0 and len(members) == 1:
                fasta_members = members
            for member in fasta_members:
                extracted = archive.extractfile(member)
                if extracted is None:
                    continue
                with io.TextIOWrapper(extracted, encoding="utf-8") as handle:
                    yield from _iter_fasta_records_from_handle(handle)
        return

    with open_text(path, "rt") as handle:
        yield from _iter_fasta_records_from_handle(handle)


def count_fasta_records(path):
    count = 0
    first_name = ""
    for header, _ in iter_fasta_records(path):
        count += 1
        if first_name == "":
            first_name = first_token(header)
    return count, first_name


def extract_ensembl_id(header):
    match = re.search(r"(?:^|\s)gene:([^\s]+)", header)
    if match:
        return match.group(1)
    return first_token(header)


def extract_phycocosm_id(header):
    token = first_token(header)
    if token.startswith("jgi|"):
        parts = [part for part in token.split("|") if part != ""]
        if len(parts) >= 2:
            return parts[-1]
    if "|" in token:
        return token.split("|")[-1]
    return token


def extract_phytozome_id(header):
    return first_token(header)


def extract_coge_id(header):
    fields = [str(field or "").strip() for field in str(header or "").split("||")]
    if len(fields) >= 5 and fields[4] != "":
        return fields[4]
    return first_token(header)


def extract_gwh_id(header):
    for tag in ("Gene", "OriGeneID"):
        candidate = extract_header_tag_value(header, tag)
        if candidate != "":
            return candidate
    return first_token(header)


def extract_provider_id(provider, header):
    if provider in ENSEMBL_LIKE_PROVIDERS:
        return extract_ensembl_id(header)
    if provider == "phycocosm":
        return extract_phycocosm_id(header)
    if provider == "coge":
        return extract_coge_id(header)
    if provider == "gwh":
        return extract_gwh_id(header)
    if provider in ("ncbi", "refseq", "genbank", "plantgarden", "plantaedb"):
        return first_token(header)
    return extract_phytozome_id(header)


def extract_provider_transcript_id(provider, header):
    token = first_token(header)
    if provider in ENSEMBL_LIKE_PROVIDERS:
        token = re.sub(r"^(?:transcript|mrna|rna):", "", token, flags=re.IGNORECASE)
        if token != "":
            return token
    extracted = extract_provider_id(provider, header)
    if extracted != "":
        return extracted
    return token


def collapse_transcript_suffix(provider, identifier):
    text = identifier
    if provider == "phytozome":
        text = re.sub(r"[.][0-9]+$", "", text)
    if provider in (
        "phycocosm",
        "phytozome",
        "ensembl",
        "ensemblplants",
        "ensemblmetazoa",
        "ensemblprotists",
        "ncbi",
        "refseq",
        "genbank",
        "coge",
        "cngb",
        "gwh",
        "flybase",
        "wormbase",
        "vectorbase",
        "fernbase",
        "insectbase",
        ORYZA_MINUTA_PROVIDER,
        "direct",
        "local",
    ):
        text = re.sub(r"[._-]t[0-9]+$", "", text, flags=re.IGNORECASE)
        text = re.sub(r"[._-](?:transcript|mrna|rna|isoform)[._-]?[0-9]+$", "", text, flags=re.IGNORECASE)
        text = re.sub(r"[._-]amt[0-9]+$", "", text, flags=re.IGNORECASE)
    return text


def normalize_gene_grouping_mode(value):
    text = str(value or "").strip().lower()
    if text == "":
        return "strict"
    if text not in GENE_GROUPING_MODES:
        raise ValueError(
            "Unknown gene_grouping_mode: {} (allowed: {})".format(
                text,
                ",".join(GENE_GROUPING_MODES),
            )
        )
    return text


def gene_grouping_mode_for_task(task):
    return normalize_gene_grouping_mode(task.get("gene_grouping_mode", "strict"))
