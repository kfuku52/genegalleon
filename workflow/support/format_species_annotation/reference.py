"""Keep GFF reference names consistent with the exported genome FASTA."""

import io
import re
import tarfile
from pathlib import Path

from format_species_common import is_fasta_filename
from format_species_constants import FASTA_ARCHIVE_EXTENSIONS
from format_species_writers import apply_common_replacements, open_text

from .common import build_gff_genome_seqid_map, first_token
from .organelle import gff_organelle_seqids, iter_non_organelle_gff_lines


class GenomeReferenceIndex(dict):
    """Header aliases and lengths, without retaining genome sequences."""

    def __init__(self):
        super().__init__()
        self.canonical_ids = {}


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
        normalized = apply_common_replacements(header)
        token = first_token(normalized)
        if not token or token in seen:
            raise ValueError("Empty or duplicate genome FASTA ID: '{}'".format(token))
        seen.add(token)
        aliases = {first_token(header), token, normalized}
        original = re.search(r"(?:^|\s)OriSeqID=([^\s;]+)", header)
        if original:
            aliases.update((original[1], apply_common_replacements(original[1])))
            declared = re.search(r"(?:^|\s)Len=(\d+)(?:\s|$)", header)
            if declared and int(declared[1]) != length:
                raise ValueError("FASTA OriSeqID length disagrees with sequence for '{}'".format(token))
        for alias in aliases:
            previous = index.canonical_ids.get(alias)
            if previous is not None and previous != token:
                raise ValueError("Ambiguous genome FASTA alias: '{}'".format(alias))
            index[alias] = length
            index.canonical_ids[alias] = token

    for handle in genome_text_handles(path):
        line_start = True
        for line in iter(lambda handle=handle: handle.readline(1024 * 1024), ""):
            if line_start and line.startswith(">"):
                if not line.endswith("\n") and len(line) == 1024 * 1024:
                    raise ValueError("Genome FASTA header exceeds 1 MiB")
                finish()
                header, length = line[1:].strip(), 0
            elif line.strip():
                if header is None:
                    raise ValueError("Genome sequence precedes its FASTA header")
                length += len(re.sub(r"\s+", "", line))
            line_start = line.endswith("\n")
        finish()
        header, length = None, 0
    if before != tuple(getattr(path.stat(), field) for field in fields):
        raise OSError("Genome FASTA changed while reading reference names: " + str(path))
    if not seen:
        raise ValueError("Genome FASTA contains no records")
    return index


def gff_reference_mapping(gff_path, genome_path):
    if genome_path is None:
        return None
    seqids = set()
    excluded = gff_organelle_seqids(gff_path)
    for line in iter_non_organelle_gff_lines(gff_path):
        parts = line.rstrip("\r\n").split("\t")
        if not line.startswith("#") and len(parts) >= 9:
            seqids.add(parts[0])
        elif line.startswith("##sequence-region "):
            directive = line.split()
            if len(directive) == 4 and directive[1] not in excluded:
                seqids.add(directive[1])
    mapping, missing = build_gff_genome_seqid_map(genome_reference_index(genome_path), seqids)
    if missing:
        raise ValueError("GFF references absent from genome FASTA: " + ", ".join(missing[:20]))
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
        checked += 1
    return checked
