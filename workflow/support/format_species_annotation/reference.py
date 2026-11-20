"""Keep GFF reference names consistent with the exported genome FASTA."""

import io
import re
import tarfile
from pathlib import Path

from format_species_common import is_fasta_filename
from format_species_constants import FASTA_ARCHIVE_EXTENSIONS
from format_species_writers import apply_common_replacements, open_text
from genome_fasta_repair import fasta_blocks

from .common import build_gff_genome_seqid_map, first_token, parse_gff_attributes
from .organelle import gff_organelle_seqids, iter_non_organelle_gff_lines


class GenomeReferenceIndex(dict):
    """Header aliases and lengths, without retaining genome sequences."""

    def __init__(self):
        super().__init__()
        self.canonical_ids = {}


def genome_header_aliases(header):
    normalized = apply_common_replacements(header.strip())
    token = first_token(normalized)
    aliases = {first_token(header), token, normalized}
    # CoGe whole-genome sequence exports wrap the chromosome name in gi|.
    # Do not interpret multi-field NCBI GI/accession headers as this wrapper.
    wrapped = re.fullmatch(r"gi\|([^|\s]+)", token)
    if wrapped:
        aliases.add(wrapped[1])
    original = re.search(r"(?:^|\s)OriSeqID=([^\s;]+)", header)
    declared = None
    if original:
        aliases.update((original[1], apply_common_replacements(original[1])))
        match = re.search(r"(?:^|\s)Len=(\d+)(?:\s|$)", header)
        declared = int(match[1]) if match else None
    return token, aliases, declared


def add_genome_reference(index, seen, header, length):
    token, aliases, declared = genome_header_aliases(header)
    if not token or token in seen:
        raise ValueError("Empty or duplicate genome FASTA ID: '{}'".format(token))
    seen.add(token)
    if declared is not None and declared != length:
        raise ValueError("FASTA OriSeqID length disagrees with sequence for '{}'".format(token))
    for alias in aliases:
        previous = index.canonical_ids.get(alias)
        if previous is not None and previous != token:
            raise ValueError("Ambiguous genome FASTA alias: '{}'".format(alias))
        index[alias] = length
        index.canonical_ids[alias] = token


_ASCII_WHITESPACE = dict.fromkeys(code for code in range(128) if chr(code).isspace())


def _strip_sequence_whitespace(text):
    if text.isascii():
        return text.translate(_ASCII_WHITESPACE)
    return "".join(text.split())


def genome_fragments(path, chunk_bytes=1024 * 1024, *, strict=True, repairs=None):
    """Yield bounded headers/sequence fragments, including unwrapped FASTA."""
    for handle in genome_text_handles(Path(path)):
        have_header = False
        edits = []
        for data, headers, _ in fasta_blocks(handle, edits, chunk_size=chunk_bytes):
            cursor = 0
            for start, stop in headers:
                if text := _strip_sequence_whitespace(data[cursor:start]):
                    if have_header:
                        yield "sequence", text
                    elif strict:
                        raise ValueError("Genome sequence precedes its FASTA header")
                yield "header", data[start + 1:stop]
                have_header = True
                cursor = stop + 1
            if text := _strip_sequence_whitespace(data[cursor:]):
                if have_header:
                    yield "sequence", text
                elif strict:
                    raise ValueError("Genome sequence precedes its FASTA header")
        if repairs is not None and edits:
            stream = getattr(handle, "fasta_member", path)
            repairs.extend({**edit, "stream": str(stream)} for edit in edits)


def genome_text_handles(path):
    if any(path.name.lower().endswith(suffix) for suffix in FASTA_ARCHIVE_EXTENSIONS):
        with tarfile.open(path, "r:*") as archive:
            members = [member for member in archive.getmembers() if member.isfile()]
            selected = [member for member in members if is_fasta_filename(Path(member.name).name)]
            if not selected and len(members) == 1:
                selected = members
            for member in selected:
                extracted = archive.extractfile(member)
                if extracted is not None:
                    with io.TextIOWrapper(extracted, encoding="utf-8") as handle:
                        handle.fasta_member = member.name
                        yield handle
    else:
        with open_text(path, "rt") as handle:
            yield handle


def genome_reference_index(path):
    path = Path(path)
    fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
    before = tuple(getattr(path.stat(), field) for field in fields)
    index = GenomeReferenceIndex()
    seen = set()
    header, length = None, 0

    def finish():
        if header is None:
            return
        add_genome_reference(index, seen, header, length)

    for kind, text in genome_fragments(path):
        if kind == "header":
            finish()
            header, length = text.strip(), 0
        elif text:
            if header is None:
                raise ValueError("Genome sequence precedes its FASTA header")
            length += len(text)
    finish()
    if before != tuple(getattr(path.stat(), field) for field in fields):
        raise OSError("Genome FASTA changed while reading reference names: " + str(path))
    if not seen:
        raise ValueError("Genome FASTA contains no records")
    return index


def gff_reference_mapping(gff_path, genome_path, *, reference_index=None):
    if genome_path is None:
        return None
    seqids, declared = set(), set()
    excluded = gff_organelle_seqids(gff_path)
    for line in iter_non_organelle_gff_lines(gff_path):
        parts = line.rstrip("\r\n").split("\t")
        if not line.startswith("#") and len(parts) >= 9:
            seqids.add(parts[0])
        elif line.startswith("##sequence-region "):
            directive = line.split()
            if len(directive) == 4 and directive[1] not in excluded:
                declared.add(directive[1])
    index = genome_reference_index(genome_path) if reference_index is None else reference_index
    mapping, missing = build_gff_genome_seqid_map(index, seqids)
    if missing:
        raise ValueError("GFF references absent from genome FASTA: " + ", ".join(missing[:20]))
    # Earlier exports removed organelle features/genome records while leaving
    # their sequence-region directives. An unused declaration is not a coding
    # coordinate; keep it only when the exported genome contains that reference.
    declared_mapping, _unused_missing = build_gff_genome_seqid_map(index, declared - seqids)
    mapping.update(declared_mapping)
    return mapping


def normalize_gff_reference_lines(lines, mapping):
    for line in lines:
        if mapping is None:
            yield line
            continue
        stripped = line.rstrip("\r\n")
        newline = line[len(stripped):]
        if stripped.startswith("##sequence-region "):
            parts = stripped.split()
            if len(parts) == 4 and parts[1] not in mapping:
                continue
            if len(parts) == 4:
                parts[1] = mapping[parts[1]]
                line = " ".join(parts) + newline
        elif not stripped.startswith("#"):
            parts = stripped.split("\t")
            if len(parts) >= 9 and parts[0] in mapping:
                parts[0] = mapping[parts[0]]
                line = "\t".join(parts) + newline
        yield line


def validate_gff_genome_references(gff_path, genome_path):
    index = genome_reference_index(genome_path)
    canonical = {name: index[name] for name in set(index.canonical_ids.values())}
    checked = 0
    gene_axes = {}
    rna_parents = []
    for line in iter_non_organelle_gff_lines(gff_path):
        if line.startswith("##sequence-region "):
            parts = line.split()
            if len(parts) == 4 and (parts[1] not in canonical or not 1 <= int(parts[2]) <= int(parts[3]) <= canonical[parts[1]]):
                raise ValueError("Formatted GFF sequence-region is incompatible with genome FASTA: " + line.strip())
        if line.startswith("#"):
            continue
        parts = line.rstrip("\r\n").split("\t")
        if len(parts) < 9:
            continue
        if parts[0] not in canonical:
            raise ValueError("Formatted GFF reference '{}' is absent from genome FASTA IDs".format(parts[0]))
        start, end = int(parts[3]), int(parts[4])
        if not 1 <= start <= end <= canonical[parts[0]]:
            raise ValueError("GFF coordinates exceed genome sequence '{}': {}-{}".format(parts[0], start, end))
        feature_type = parts[2].lower()
        if feature_type in ("gene", "mrna", "rna", "transcript"):
            keys = ("ID",) if feature_type == "gene" else ("Parent",)
            attrs = parse_gff_attributes(parts[8], keys=keys)
            if feature_type == "gene":
                for gene_id in attrs.get("ID", ()):
                    gene_axes.setdefault(gene_id, set()).add((parts[0], parts[6]))
            else:
                for parent in attrs.get("Parent", ()):
                    rna_parents.append((parent, parts[0], parts[6]))
        checked += 1
    for parent, seqid, strand in rna_parents:
        axes = gene_axes.get(parent)
        # A missing gene row can still carry an authoritative declared Parent.
        # Validate existing rows without inventing a replacement gene identity.
        if axes and not any(
            gene_seqid == seqid and (gene_strand not in ("+", "-") or strand not in ("+", "-") or gene_strand == strand)
            for gene_seqid, gene_strand in axes
        ):
            raise ValueError("Formatted GFF gene/RNA reference or strand disagreement for Parent '{}'".format(parent))
    return checked
