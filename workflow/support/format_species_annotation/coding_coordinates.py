"""Recover coding intervals only from unique, exact published CDS evidence.

Some publishers distribute a gene/exon-only GFF and one CDS per gene. Exons
may include UTRs, so they cannot simply be relabelled as CDS. Require a complete
one-to-one identity match and a unique CDS occurrence in each spliced model.
The shared view is private; all original source files remain untouched.
"""
import hashlib
import tempfile
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from urllib.parse import quote

from format_species_writers import open_text
from gff_attribute_syntax import validated_lines

from .common import first_token, iter_fasta_records, parse_gff_attributes
from .genome_intervals import genome_intervals, stat_identity

VERSION = 1
_scratch = tempfile.TemporaryDirectory(prefix="gg-coding-coordinates-")
_complement = str.maketrans("ACGTRYMKBDHVN", "TGCAYRKMVHDBN")


def unique_exact_offset(spliced, coding):
    offset = spliced.find(coding)
    if offset < 0 or spliced.find(coding, offset + 1) >= 0:
        return None
    if set(spliced + coding) <= set("ACGT"):
        return offset
    if set(spliced + coding) - set("ACGTN"):
        return None
    # N is unknown, not a discriminating base. Shift-And checks every possible
    # alignment under N compatibility as well as the literal exact occurrence.
    masks = dict.fromkeys("ACGTN", 0)
    for i, base in enumerate(coding):
        for symbol in ("ACGTN" if base == "N" else base + "N"):
            masks[symbol] |= 1 << i
    state, terminal, matches = 0, 1 << (len(coding) - 1), 0
    for base in spliced:
        state = ((state << 1) | 1) & masks[base]
        if state & terminal:
            matches += 1
            if matches > 1:
                return None
    return offset if matches == 1 else None


def signature(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def coding_view(task, source):
    if (str(task.get("gff_repair_mode", "safe")).strip().lower() == "off"
            or not task.get("cds_path") or not task.get("genome_path")):
        return Path(source), {"version": VERSION, "status": "not_applicable"}
    paths = tuple(str(Path(path).resolve()) for path in (source, task["cds_path"], task["genome_path"]))
    identities = tuple(stat_identity(path) for path in paths)
    result = _cached_view(paths, identities)
    if tuple(stat_identity(path) for path in paths) != identities:
        raise OSError("Sources changed during coding-coordinate reconstruction")
    return result


@lru_cache(maxsize=4)
def _cached_view(paths, identities):
    try:
        return reconstruct(*paths)
    finally:
        if tuple(stat_identity(path) for path in paths) != identities:
            raise OSError("Sources changed during coding-coordinate reconstruction")


def reconstruct(gff, cds, genome):
    # Ordinary annotations retain their existing path, including their native
    # source-syntax repairs. Do not apply this adapter to mixed/CDS/RNA GFFs.
    kinds = set()
    with open_text(gff, "rt") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) != 9 or parts[2].lower() not in {"gene", "exon"}:
                return Path(gff), {"version": VERSION, "status": "not_applicable"}
            kinds.add(parts[2].lower())
    if kinds != {"gene", "exon"}:
        return Path(gff), {"version": VERSION, "status": "not_applicable"}
    genes, exons, identifiers = {}, defaultdict(list), set()
    with open_text(gff, "rt") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) != 9:
                return Path(gff), {"version": VERSION, "status": "not_applicable"}
            kind = parts[2].lower()
            if kind not in {"gene", "exon"}:
                return Path(gff), {"version": VERSION, "status": "not_applicable"}
            attrs = parse_gff_attributes(parts[8])
            ids = attrs.get("ID", ())
            if len(ids) > 1 or (ids and ids[0] in identifiers):
                raise ValueError("Duplicate coding-coordinate source feature ID")
            identifiers.update(ids)
            start, end = int(parts[3]), int(parts[4])
            if start < 1 or end < start or parts[6] not in {"+", "-"}:
                raise ValueError("Invalid coding-coordinate source axis")
            if kind == "gene":
                if len(ids) != 1 or attrs.get("Parent"):
                    raise ValueError("Coding-coordinate reconstruction requires explicit gene IDs")
                genes[ids[0]] = parts
            else:
                parents = attrs.get("Parent", ())
                if len(parents) != 1:
                    raise ValueError("Coding-coordinate reconstruction requires one exon parent")
                exons[parents[0]].append(parts)
    sequences = {}
    for header, sequence in iter_fasta_records(Path(cds)):
        name, sequence = first_token(header), sequence.upper()
        if name in sequences or not sequence or set(sequence) - set("ACGTN"):
            raise ValueError("Coding-coordinate reconstruction requires unique published DNA CDS")
        sequences[name] = sequence
    if set(genes) != set(exons) or set(genes) != set(sequences):
        raise ValueError("Coding-coordinate gene/exon/published-CDS IDs do not match completely")
    for gene, rows in exons.items():
        parent = genes[gene]
        rows.sort(key=lambda row: int(row[3]))
        for i, row in enumerate(rows):
            if ((row[0], row[6]) != (parent[0], parent[6])
                    or int(row[3]) < int(parent[3]) or int(row[4]) > int(parent[4])
                    or (i and int(rows[i - 1][4]) >= int(row[3]))):
                raise ValueError("Conflicting or overlapping coding-coordinate exons: " + gene)
        if parent[6] == "-":
            rows.reverse()
    intervals = ((row[0], int(row[3]) - 1, int(row[4])) for rows in exons.values() for row in rows)
    segments, offsets, transcript_ids = {}, {}, {}
    with genome_intervals(genome, intervals, scratch_dir=_scratch.name) as reader:
        for gene, rows in exons.items():
            ref = genes[gene][0]
            # This reconstruction requires exact publisher reference identities.
            if ref not in reader.regions:
                raise ValueError("Missing exact coding-coordinate genome reference: " + ref)
            pieces = [reader.fetch(ref, int(row[3]) - 1, int(row[4])) for row in rows]
            if genes[gene][6] == "-":
                pieces = [piece.translate(_complement)[::-1] for piece in pieces]
            spliced, coding = "".join(pieces), sequences[gene]
            offset = unique_exact_offset(spliced, coding)
            if offset is None:
                raise ValueError("Published CDS does not have one exact spliced-exon occurrence: " + gene)
            transcript = gene + ".gg_published_cds"
            if transcript in identifiers:
                raise ValueError("Reconstructed transcript ID collides with source: " + transcript)
            transcript_ids[gene], offsets[gene] = transcript, offset
            cursor, emitted = 0, 0
            segments[gene] = []
            for row, piece in zip(rows, pieces, strict=True):
                left, right = max(cursor, offset), min(cursor + len(piece), offset + len(coding))
                if left < right:
                    a, b = left - cursor, right - cursor
                    if row[6] == "+":
                        start, end = int(row[3]) + a, int(row[3]) + b - 1
                    else:
                        start, end = int(row[4]) - b + 1, int(row[4]) - a
                    segments[gene].append((start, end, (3 - emitted % 3) % 3))
                    emitted += right - left
                cursor += len(piece)
            if emitted != len(coding):
                raise ValueError("Incomplete published coding-coordinate coverage: " + gene)
    hashes = {role: signature(path) for role, path in zip(("gff", "cds", "genome"), (gff, cds, genome), strict=True)}
    digest = hashlib.sha256((str(VERSION) + str(sorted(hashes.items()))).encode()).hexdigest()
    output = Path(_scratch.name) / (digest + ".gff")
    with output.open("w") as dest:
        with open_text(gff, "rt") as source:
            for line in validated_lines(source, gff):
                if line.startswith("#") or not line.strip():
                    dest.write(line)
                    continue
                parts = line.rstrip("\r\n").split("\t")
                attrs = parse_gff_attributes(parts[8])
                if parts[2].lower() == "exon":
                    # Keep the original exon/UTR bounds under the explicit RNA.
                    fields = parts[8].split(";")
                    parts[8] = ";".join("Parent=" + quote(transcript_ids[attrs["Parent"][0]], safe="._-:")
                                        if field.strip().startswith("Parent=") else field for field in fields)
                    dest.write("\t".join(parts) + "\n")
                    continue
                dest.write(line if line.endswith("\n") else line + "\n")
                gene = attrs["ID"][0]
                transcript = quote(transcript_ids[gene], safe="._-:")
                basis = ";annotation_basis=unique_exact_published_CDS_on_spliced_exons"
                rna = parts.copy()
                rna[1:3], rna[7] = ["GeneGalleon", "mRNA"], "."
                rna[8] = ("ID=" + transcript + ";Parent=" + quote(gene, safe="._-:") + basis
                          + ";spliced_cds_offset=" + str(offsets[gene]))
                dest.write("\t".join(rna) + "\n")
                for start, end, phase in segments[gene]:
                    coding_row = rna.copy()
                    coding_row[2:5], coding_row[7] = ["CDS", str(start), str(end)], str(phase)
                    coding_row[8] = "Parent=" + transcript + basis
                    dest.write("\t".join(coding_row) + "\n")
    return output, {"version": VERSION, "status": "reconstructed", "source_sha256": hashes,
                    "gene_models": len(genes), "published_cds_records": len(sequences),
                    "cds_models_with_unknown_bases": sum("N" in sequence for sequence in sequences.values()),
                    "cds_segments": sum(map(len, segments.values())), "output_sha256": signature(output),
                    "evidence": "one unique exact published CDS per complete spliced-exon model"}
