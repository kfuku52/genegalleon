"""Recover author model identities from supported Btu and Lavan exports.

Bambusa tulda exports reuse locus tags and point CDS at unrelated genes. The
original Btu gene/transcript IDs survive in every model's Note. Reconstruct
only this explicit convention, retaining coordinates, phases and protein IDs.
All consumers use the same cached annotation view; source fingerprints continue
to refer to the untouched input file.
"""

import hashlib
import re
import tempfile
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from urllib.parse import quote

from format_species_writers import open_text

from .common import parse_gff_attributes

SOURCE_IDENTITY_VERSION = 1
AUTHOR_ID = re.compile(r"(?:^|[;~\s])ID:(Btu[A-Za-z0-9.]*[.]g\d+)(t\d+)(?:[.]CDS)?(?:$|[;~\s])")
LAVAN_MODEL = re.compile(r"^(Lavan[.](?:\d+G\d+|S\d+))[.]\d+$")
_scratch = tempfile.TemporaryDirectory(prefix="gg-source-identity-")


def author_identity(attrs):
    matches = {m.groups() for note in attrs.get("Note", ()) for m in AUTHOR_ID.finditer(note)}
    if len(matches) > 1:
        raise ValueError("Conflicting original Btu identities in GFF Note")
    if not matches:
        return None
    gene, suffix = next(iter(matches))
    return gene, gene + suffix


def attributes_text(attrs):
    return ";".join(key + "=" + ",".join(quote(str(v), safe=":._-|@") for v in values)
                    for key, values in attrs.items() if values)


def source_annotation_path(path, *, repair_locus_ids=False):
    source = Path(path).resolve()
    stat = source.stat()
    view = _annotation_view(str(source), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
    if repair_locus_ids:
        return _locus_view(view)[0]
    return view


def task_annotation_path(task):
    return source_annotation_path(task["gff_path"],
        repair_locus_ids=str(task.get("gff_repair_mode", "safe")).strip().lower() != "off")


def locus_identity_audit(task):
    if not task.get("gff_path") or str(task.get("gff_repair_mode", "safe")).strip().lower() == "off":
        return {"version": 1, "status": "off", "mappings": [], "problems": []}
    return _locus_view(source_annotation_path(task["gff_path"]))[1]


def _locus_view(source):
    source = Path(source)
    stat = source.stat()
    return _cached_locus_view(str(source), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


@lru_cache(maxsize=4)
def _cached_locus_view(source, _size, _mtime, _ctime):
    from .locus_identity import normalise_locus_ids
    try:
        return normalise_locus_ids(source, _scratch.name)
    finally:
        stat = Path(source).stat()
        if (stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns) != (_size, _mtime, _ctime):
            raise OSError("GFF changed during locus identity repair: " + str(source))


@lru_cache(maxsize=4)
def _annotation_view(source, _size, _mtime, _ctime):
    genes, transcripts = {}, {}
    feature_counts = defaultdict(int)
    missing_identity = 0
    lavandula_models = {}
    other_rna_models = 0
    conflicting_lavan = False
    active = False
    # A full first pass ensures mixed conventions and inconsistent models fail
    # before any repaired annotation is exposed to another consumer.
    with open_text(source, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("##genegalleon-original-id-normalization version={}".format(SOURCE_IDENTITY_VERSION)):
                return Path(source)
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue
            kind = parts[2].lower()
            feature_counts[kind] += 1
            if kind not in ("mrna", "cds", "exon"):
                continue
            attrs = parse_gff_attributes(parts[8])
            if kind == "mrna":
                identifier = next(iter(attrs.get("ID", ())), "")
                match = LAVAN_MODEL.fullmatch(identifier)
                if match:
                    if attrs.get("Parent") or identifier in lavandula_models:
                        conflicting_lavan = True
                    lavandula_models[identifier] = dict(gene=match[1], axis=(parts[0], parts[6]),
                                                      start=int(parts[3]), end=int(parts[4]), source=parts[1])
                else:
                    other_rna_models += 1
            identity = author_identity(attrs)
            if identity is None:
                missing_identity += 1
                continue
            active = True
            gene, transcript = identity
            axis = parts[0], parts[6]
            start, end = int(parts[3]), int(parts[4])
            if start > end or axis[1] not in ("+", "-"):
                raise ValueError("Invalid original Btu model coordinates: " + transcript)
            key = gene, axis
            model = genes.setdefault(key, dict(start=start, end=end, source=parts[1]))
            model["start"], model["end"] = min(model["start"], start), max(model["end"], end)
            key = transcript, axis
            model = transcripts.setdefault(key, dict(rna=None, cds=[], exon=[], pseudo_cds=True))
            if kind == "mrna":
                span = start, end
                if model["rna"] is not None and model["rna"] != span:
                    raise ValueError("Conflicting original Btu RNA model: " + transcript)
                model["rna"] = span
            else:
                model[kind].append((start, end))
                if kind == "cds" and "true" not in attrs.get("pseudo", ()):
                    model["pseudo_cds"] = False
    if not active:
        if lavandula_models and not feature_counts["gene"]:
            if other_rna_models or conflicting_lavan:
                raise ValueError("Mixed original Lavan RNA identity conventions")
            return _lavandula_view(source, lavandula_models)
        return Path(source)
    if missing_identity:
        raise ValueError("Mixed or missing original Btu model identities in GFF")
    for (transcript, _axis), model in transcripts.items():
        rna = model["rna"]
        if rna is None and model["cds"] and model["pseudo_cds"] and not model["exon"]:
            # The export also contains an explicitly marked pseudogene CDS
            # without RNA on another contig. Retain its coding coordinates as
            # a pseudo transcript; do not borrow the RNA from the other locus.
            model["pseudo_transcript"] = True
            continue
        if rna is None or any(start < rna[0] or end > rna[1] for start, end in model["cds"] + model["exon"]):
            raise ValueError("Original Btu CDS/exon is outside its RNA model: " + transcript)

    axes_by_id = defaultdict(set)
    for identifier, axis in list(genes) + list(transcripts):
        axes_by_id[identifier].add(axis)

    def canonical(identifier, axis):
        if len(axes_by_id[identifier]) == 1:
            return identifier
        # Original author IDs can recur on another contig or strand. They are
        # different physical loci and cannot share a selected CDS identity.
        return identifier + ".locus" + hashlib.sha256((axis[0] + "|" + axis[1]).encode()).hexdigest()[:12]

    digest = hashlib.sha256(Path(source).read_bytes()).hexdigest()
    output = Path(_scratch.name) / (digest + ".gff")
    with output.open("w") as dest:
        dest.write("##gff-version 3\n")
        dest.write("##genegalleon-original-id-normalization version={} source_sha256={} "
                   "replaced_gene_features={} author_gene_features={}\n".format(
                       SOURCE_IDENTITY_VERSION, digest, feature_counts["gene"], len(genes)))
        for (gene, axis), model in sorted(genes.items()):
            identifier = canonical(gene, axis)
            attrs = dict(ID=[identifier], locus_tag=[identifier], orig_author_gene_id=[gene])
            dest.write("\t".join((axis[0], model["source"], "gene", str(model["start"]),
                                  str(model["end"]), ".", axis[1], ".", attributes_text(attrs))) + "\n")
        for (transcript, axis), model in sorted(transcripts.items()):
            if not model.get("pseudo_transcript"):
                continue
            gene = re.sub(r"t\d+$", "", transcript)
            attrs = dict(ID=[canonical(transcript, axis)], Parent=[canonical(gene, axis)],
                         pseudo=["true"], orig_author_transcript_id=[transcript],
                         annotation_basis=["original_Note_and_pseudo_CDS"])
            dest.write("\t".join((axis[0], "GeneGalleon", "transcript",
                                  str(min(start for start, _end in model["cds"])),
                                  str(max(end for _start, end in model["cds"])),
                                  ".", axis[1], ".", attributes_text(attrs))) + "\n")
        with open_text(source, "rt", errors="replace") as handle:
            for line in handle:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("##gff-version"):
                    continue
                parts = line.rstrip("\r\n").split("\t")
                if line.startswith("#") or len(parts) < 9:
                    dest.write(line)
                    continue
                kind = parts[2].lower()
                if kind == "gene":
                    continue
                attrs = parse_gff_attributes(parts[8])
                if kind not in ("mrna", "cds", "exon"):
                    if attrs.get("Parent"):
                        raise ValueError("Unsupported parent relationship in original Btu annotation")
                    dest.write(line)
                    continue
                gene, transcript = author_identity(attrs)
                axis = parts[0], parts[6]
                gene_id, rna_id = canonical(gene, axis), canonical(transcript, axis)
                for key in ("ID", "Parent", "locus_tag"):
                    if attrs.get(key):
                        attrs["orig_export_" + key.lower()] = attrs[key]
                attrs["locus_tag"] = [gene_id]
                attrs["orig_author_transcript_id"] = [transcript]
                if kind == "mrna":
                    attrs["ID"], attrs["Parent"], attrs["Name"] = [rna_id], [gene_id], [rna_id]
                elif kind == "cds":
                    attrs["Parent"] = [rna_id]
                    # Keep CDS/protein accession aliases, including anonymous
                    # accessions used with exact NCBI FASTA locations.
                else:
                    attrs["Parent"] = [rna_id]
                    attrs["ID"] = ["exon-{}-{}-{}".format(rna_id, parts[3], parts[4])]
                parts[8] = attributes_text(attrs)
                dest.write("\t".join(parts) + "\n")
    return output


def _lavandula_view(source, models):
    """Restore the documented Lavan author gene number, including disjoint CDS.

    This export omitted gene features, while all RNA IDs explicitly encode the
    author's chromosome/unplaced gene number followed by the isoform number.
    These IDs are stronger evidence than coordinate overlap. Generic numeric
    suffix inference elsewhere still requires a connected coding locus.
    """
    genes = {}
    axes = defaultdict(set)
    for model in models.values():
        gene, axis = model["gene"], model["axis"]
        axes[gene].add(axis)
        record = genes.setdefault((gene, axis), dict(model))
        record["start"] = min(record["start"], model["start"])
        record["end"] = max(record["end"], model["end"])

    def canonical(gene, axis):
        if len(axes[gene]) == 1:
            return gene
        return gene + ".locus" + hashlib.sha256((axis[0] + "|" + axis[1]).encode()).hexdigest()[:12]

    digest = hashlib.sha256(Path(source).read_bytes()).hexdigest()
    output = Path(_scratch.name) / (digest + ".gff")
    with output.open("w") as dest:
        dest.write("##gff-version 3\n")
        dest.write("##genegalleon-original-id-normalization version={} convention=Lavan "
                   "source_sha256={} author_gene_features={}\n".format(SOURCE_IDENTITY_VERSION, digest, len(genes)))
        for (gene, axis), model in sorted(genes.items()):
            identifier = canonical(gene, axis)
            attrs = dict(ID=[identifier], locus_tag=[identifier], orig_author_gene_id=[gene])
            dest.write("\t".join((axis[0], model["source"], "gene", str(model["start"]),
                                  str(model["end"]), ".", axis[1], ".", attributes_text(attrs))) + "\n")
        with open_text(source, "rt", errors="replace") as handle:
            for line in handle:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("##gff-version"):
                    continue
                parts = line.rstrip("\r\n").split("\t")
                if line.startswith("#") or len(parts) < 9:
                    dest.write(line)
                    continue
                attrs = parse_gff_attributes(parts[8])
                kind = parts[2].lower()
                if kind == "mrna":
                    identifier = attrs["ID"][0]
                    model = models[identifier]
                    attrs["Parent"] = [canonical(model["gene"], model["axis"])]
                    attrs["orig_author_gene_id"] = [model["gene"]]
                else:
                    for parent in attrs.get("Parent", ()):
                        if parent not in models:
                            raise ValueError("Unknown original Lavan RNA parent: " + parent)
                        model = models[parent]
                        if ((parts[0], parts[6]) != model["axis"] or int(parts[3]) < model["start"]
                                or int(parts[4]) > model["end"]):
                            raise ValueError("Original Lavan feature is outside its RNA: " + parent)
                parts[8] = attributes_text(attrs)
                dest.write("\t".join(parts) + "\n")
    return output


def iter_source_gff_lines(path):
    with open_text(source_annotation_path(path), "rt", errors="replace") as handle:
        yield from handle
