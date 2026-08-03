"""Annotation implementation: tasks."""

import re

try:
    from Bio import SeqIO
except Exception:  # pragma: no cover - runtime without biopython
    SeqIO = None

from format_species_common import (
    is_fasta_filename,
    is_gbff_filename,
    is_gff_filename,
    is_probable_cds_filename,
    is_probable_genome_filename,
    pick_single_file,
)
from format_species_constants import (
    FASTA_EXTENSIONS,
    GENBANK_EXTENSIONS,
    GFF_EXTENSIONS,
)
from format_species_taxonomy import (
    invalid_species_key_error,
    normalize_species_key_for_runtime,
)

from .common import (
    collapse_transcript_suffix,
    extract_header_tag_value,
    extract_ncbi_ensembl_gene_id_from_header,
    extract_ncbi_gene_id_from_header,
    extract_provider_id,
    fasta_looks_like_gff_genome,
    iter_fasta_records,
    normalize_output_basename,
    sanitize_identifier,
    species_prefix_from_value,
    strip_suffix_case_insensitive,
    task_has_usable_source_bundle,
)
from .genbank import (
    derive_cds_records_from_gbff,
    derive_cds_records_from_gff_and_genome,
)
from .grouping import resolve_cds_header_gff_gene


def task_missing_annotation_label(cds_path, gff_path, gbff_path, genome_path):
    if task_has_usable_source_bundle(cds_path, gff_path, gbff_path, genome_path):
        return ""
    if cds_path is None and gff_path is None and gbff_path is None and genome_path is None:
        return "CDS-or-GBFF-or-(GFF+genome)"
    if gff_path is not None and genome_path is None:
        return "CDS-or-GBFF-or-genome"
    if genome_path is not None and gff_path is None and gbff_path is None:
        return "CDS-or-GBFF-or-GFF"
    if gbff_path is None:
        return "CDS-or-GBFF-or-(GFF+genome)"
    return ""


def describe_task_cds_input(task):
    cds_path = task.get("cds_path")
    if cds_path is not None:
        return cds_path
    gff_path = task.get("gff_path")
    gbff_path = task.get("gbff_path")
    genome_path = task.get("genome_path")
    if gff_path is not None and genome_path is not None:
        return "{} + {} (derived CDS)".format(gff_path, genome_path)
    if gbff_path is not None and genome_path is not None:
        return "{} + {} (derived CDS)".format(gbff_path, genome_path)
    if gbff_path is not None:
        return "{} (derived CDS)".format(gbff_path)
    return ""


def describe_task_gff_input(task):
    gff_path = task.get("gff_path")
    if gff_path is not None:
        return gff_path
    gbff_path = task.get("gbff_path")
    if gbff_path is not None:
        return "{} (derived GFF)".format(gbff_path)
    return ""


def describe_task_genome_input(task):
    genome_path = task.get("genome_path")
    if genome_path is not None:
        return genome_path
    gbff_path = task.get("gbff_path")
    if gbff_path is not None:
        return "{} (derived genome)".format(gbff_path)
    return ""


def build_derived_cds_output_basename(task):
    gff_path = task.get("gff_path")
    species_prefix = task["species_prefix"]
    if gff_path is not None:
        normalized = normalize_output_basename(gff_path.name, species_prefix)
        stem = strip_suffix_case_insensitive(normalized, GFF_EXTENSIONS)
        return stem + ".derived.cds.fa.gz"
    gbff_path = task.get("gbff_path")
    if gbff_path is not None:
        normalized = normalize_output_basename(gbff_path.name, species_prefix)
        stem = strip_suffix_case_insensitive(normalized, GENBANK_EXTENSIONS)
        return stem + ".derived.cds.fa.gz"
    genome_path = task.get("genome_path")
    if genome_path is not None:
        normalized = normalize_output_basename(genome_path.name, species_prefix)
        stem = strip_suffix_case_insensitive(normalized, FASTA_EXTENSIONS)
        return stem + ".derived.cds.fa.gz"
    return species_prefix + ".derived.cds.fa.gz"


def build_derived_gff_output_basename(task):
    gbff_path = task.get("gbff_path")
    species_prefix = task["species_prefix"]
    if gbff_path is not None:
        normalized = normalize_output_basename(gbff_path.name, species_prefix)
        stem = strip_suffix_case_insensitive(normalized, GENBANK_EXTENSIONS)
        return stem + ".derived.gff.gz"
    return species_prefix + ".derived.gff.gz"


def build_derived_genome_output_basename(task):
    gbff_path = task.get("gbff_path")
    species_prefix = task["species_prefix"]
    if gbff_path is not None:
        normalized = normalize_output_basename(gbff_path.name, species_prefix)
        stem = strip_suffix_case_insensitive(normalized, GENBANK_EXTENSIONS)
        return stem + ".derived.genome.fa.gz"
    return species_prefix + ".derived.genome.fa.gz"


def iter_task_cds_records(task):
    cds_path = task.get("cds_path")
    if cds_path is not None:
        yield from iter_fasta_records(cds_path)
        return
    gff_path = task.get("gff_path")
    genome_path = task.get("genome_path")
    gbff_path = task.get("gbff_path")
    if gff_path is not None and genome_path is not None:
        yield from derive_cds_records_from_gff_and_genome(task)
        return
    if gbff_path is not None:
        yield from derive_cds_records_from_gbff(task)
        return
    raise ValueError("No CDS source is available for {}".format(task.get("species_key", "")))


def discover_generic_species_dir_tasks(provider, input_dir, allowed_species_keys=None):
    warnings = []
    errors = []
    tasks = []

    for species_dir in sorted(input_dir.iterdir()):
        if not species_dir.is_dir():
            continue
        species_key = normalize_species_key_for_runtime(species_dir.name)
        if allowed_species_keys is not None and species_key not in allowed_species_keys:
            continue
        invalid_species_key = invalid_species_key_error(provider, species_key)
        if invalid_species_key != "":
            errors.append(invalid_species_key)
            continue
        species_prefix = species_prefix_from_value(species_key)
        if species_prefix == "":
            warnings.append("[{}] skipped '{}': unable to parse species prefix.".format(provider, species_key))
            continue

        files = [path for path in species_dir.iterdir() if path.is_file()]
        gff_matches = [path for path in files if is_gff_filename(path.name)]
        gbff_matches = [path for path in files if is_gbff_filename(path.name)]
        cds_matches = [path for path in files if is_probable_cds_filename(provider, path.name)]
        genome_matches = [
            path for path in files if path not in cds_matches and is_probable_genome_filename(provider, path.name)
        ]
        if len(cds_matches) == 0:
            unclassified_fasta = [path for path in files if is_fasta_filename(path.name) and path not in genome_matches]
            if (
                provider == "direct"
                and len(genome_matches) == 0
                and len(unclassified_fasta) == 1
                and len(gff_matches) > 0
                and fasta_looks_like_gff_genome(unclassified_fasta[0], gff_matches)
            ):
                genome_matches = unclassified_fasta
                unclassified_fasta = []
            cds_matches = [path for path in unclassified_fasta]

        cds_path = pick_single_file(cds_matches, provider, species_key, "CDS", warnings)
        gff_path = pick_single_file(gff_matches, provider, species_key, "GFF", warnings)
        gbff_path = pick_single_file(gbff_matches, provider, species_key, "GBFF", warnings)
        genome_path = pick_single_file(genome_matches, provider, species_key, "genome", warnings)

        missing_label = task_missing_annotation_label(cds_path, gff_path, gbff_path, genome_path)
        if missing_label != "":
            errors.append("[{}] {}: missing {}".format(provider, species_key, missing_label))
            continue
        tasks.append(
            {
                "provider": provider,
                "species_key": species_key,
                "species_prefix": species_prefix,
                "cds_path": cds_path,
                "gff_path": gff_path,
                "gbff_path": gbff_path,
                "genome_path": genome_path,
            }
        )
    return tasks, warnings, errors


def build_header_gene_aggregate_id(task, header, transcript_id):
    provider = task["provider"]
    species_prefix = task["species_prefix"]
    if provider == "coge":
        extracted = extract_provider_id(provider, header)
        gene_token = task.get("_provider_gene_id_map", {}).get(extracted, "")
        if gene_token != "":
            prefixed = "{}_{}".format(species_prefix, sanitize_identifier(gene_token))
            return sanitize_identifier(prefixed)
    if provider in ("direct", "local"):
        ensembl_gene_id = extract_ncbi_ensembl_gene_id_from_header(header)
        gene_symbol = extract_header_tag_value(header, "gene")
        gene_id = extract_ncbi_gene_id_from_header(header)
        locus_tag = extract_header_tag_value(header, "locus_tag")
        protein_id = extract_header_tag_value(header, "protein_id")
        if ensembl_gene_id != "":
            gene_token = ensembl_gene_id
        elif locus_tag != "":
            gene_token = locus_tag
        elif gene_id != "":
            gene_token = "GeneID{}".format(gene_id)
        elif gene_symbol != "":
            gene_token = gene_symbol
        elif protein_id != "":
            gene_token = protein_id
        else:
            gene_token = ""
        if gene_token != "":
            prefixed = "{}_{}".format(species_prefix, sanitize_identifier(gene_token))
            return sanitize_identifier(prefixed)
    if provider in ("ncbi", "refseq", "genbank", "plantaedb"):
        ensembl_gene_id = extract_ncbi_ensembl_gene_id_from_header(header)
        gene_symbol = extract_header_tag_value(header, "gene")
        gene_id = extract_ncbi_gene_id_from_header(header)
        locus_tag = extract_header_tag_value(header, "locus_tag")
        protein_id = extract_header_tag_value(header, "protein_id")
        if ensembl_gene_id != "":
            gene_token = ensembl_gene_id
        elif locus_tag != "":
            gene_token = locus_tag
        elif gene_id != "":
            gene_token = "GeneID{}".format(gene_id)
        elif gene_symbol != "":
            gene_token = gene_symbol
        elif protein_id != "":
            gene_token = protein_id
        else:
            gene_token = extract_provider_id(provider, header)
    else:
        gene_symbol = extract_header_tag_value(header, "gene")
        extracted = extract_provider_id(provider, header)
        collapsed = collapse_transcript_suffix(provider, extracted)
        if gene_symbol != "":
            if (
                collapsed != ""
                and collapsed != gene_symbol
                and re.search(r"(?:^|[._-]){}$".format(re.escape(gene_symbol)), collapsed, re.IGNORECASE)
            ):
                gene_token = collapsed
            else:
                gene_token = gene_symbol
        else:
            gene_token = collapsed if collapsed != "" else extracted
    prefixed = "{}_{}".format(species_prefix, sanitize_identifier(gene_token))
    return sanitize_identifier(prefixed)


def resolve_gene_aggregate_id(task, header, transcript_id):
    gff_match = resolve_cds_header_gff_gene(task, header)
    if gff_match["status"] == "mapped":
        prefixed = "{}_{}".format(task["species_prefix"], sanitize_identifier(gff_match["gene_token"]))
        return sanitize_identifier(prefixed), gff_match
    return build_header_gene_aggregate_id(task, header, transcript_id), gff_match


def build_gene_aggregate_id(task, header, transcript_id):
    gene_id, _gff_match = resolve_gene_aggregate_id(task, header, transcript_id)
    return gene_id
