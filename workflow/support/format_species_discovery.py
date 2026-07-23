"""Local input discovery and formatted output writers."""

import re
import time
from collections import defaultdict
from pathlib import Path

from format_species_annotations import (
    audit_matches_inputs,
    build_coge_gff_gene_id_map,
    build_derived_cds_output_basename,
    build_derived_genome_output_basename,
    build_derived_gff_output_basename,
    build_gene_aggregate_id,
    count_fasta_records,
    describe_task_cds_input,
    discover_generic_species_dir_tasks,
    extract_provider_transcript_id,
    first_token,
    gff_repair_audit_path,
    gff_repair_mode_for_task,
    iter_fasta_records,
    iter_genome_records_from_gbff,
    iter_gff_lines_from_gbff,
    iter_task_cds_records,
    normalize_cds_output_basename,
    normalize_genome_output_basename,
    normalize_gff_output_basename,
    pad_to_codon_length,
    read_json,
    repair_result_fields,
    sanitize_identifier,
    species_prefix_from_value,
    task_missing_annotation_label,
    write_repaired_gff,
)
from format_species_common import (
    is_fasta_filename,
    is_gbff_filename,
    is_gff_filename,
    is_probable_genome_filename,
    pick_single_file,
)
from format_species_constants import ENSEMBL_CDS_PATTERN
from format_species_provider_config import ENSEMBL_LIKE_PROVIDERS, ORYZA_MINUTA_PROVIDER
from format_species_taxonomy import invalid_species_key_error, normalize_species_key_for_runtime
from format_species_writers import (
    apply_common_replacements,
    write_fasta_records_gzip,
    write_gff_gzip,
    write_gff_lines_gzip,
)


def discover_ensembl_like_tasks(input_dir, provider, allowed_species_keys=None):
    warnings = []
    errors = []
    tasks = []

    cds_by_species = defaultdict(list)
    gff_by_species = defaultdict(list)
    gbff_by_species = defaultdict(list)
    genome_by_species = defaultdict(list)

    for path in sorted(input_dir.iterdir()):
        if not path.is_file():
            continue
        name = path.name
        species_key = normalize_species_key_for_runtime(species_prefix_from_value(name))
        if species_key == "":
            continue
        if is_fasta_filename(name) and ENSEMBL_CDS_PATTERN.search(name):
            cds_by_species[species_key].append(path)
            continue
        if is_gff_filename(name):
            gff_by_species[species_key].append(path)
            continue
        if is_gbff_filename(name):
            gbff_by_species[species_key].append(path)
            continue
        if is_probable_genome_filename(provider, name):
            genome_by_species[species_key].append(path)

    for species_key in sorted(
        set(cds_by_species.keys())
        | set(gff_by_species.keys())
        | set(gbff_by_species.keys())
        | set(genome_by_species.keys())
    ):
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

        cds_path = pick_single_file(
            cds_by_species.get(species_key, []),
            provider,
            species_key,
            "CDS",
            warnings,
        )
        gff_path = pick_single_file(
            gff_by_species.get(species_key, []),
            provider,
            species_key,
            "GFF",
            warnings,
        )
        gbff_path = pick_single_file(
            gbff_by_species.get(species_key, []),
            provider,
            species_key,
            "GBFF",
            warnings,
        )
        genome_path = pick_single_file(
            genome_by_species.get(species_key, []),
            provider,
            species_key,
            "genome",
            warnings,
        )
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


def discover_phycocosm_tasks(input_dir, allowed_species_keys=None):
    warnings = []
    errors = []
    tasks = []

    for species_dir in sorted(input_dir.iterdir()):
        if not species_dir.is_dir():
            continue
        species_key = normalize_species_key_for_runtime(species_dir.name)
        if allowed_species_keys is not None and species_key not in allowed_species_keys:
            continue
        invalid_species_key = invalid_species_key_error("phycocosm", species_key)
        if invalid_species_key != "":
            errors.append(invalid_species_key)
            continue
        species_prefix = species_prefix_from_value(species_key)
        if species_prefix == "":
            warnings.append("[phycocosm] skipped '{}': unable to parse species prefix.".format(species_key))
            continue

        files = [path for path in species_dir.iterdir() if path.is_file()]
        cds_matches = [path for path in files if "fasta" in path.name.lower() and is_fasta_filename(path.name)]
        gff_matches = [path for path in files if "gff" in path.name.lower() and is_gff_filename(path.name)]
        gbff_matches = [path for path in files if is_gbff_filename(path.name)]
        genome_matches = [path for path in files if is_probable_genome_filename("phycocosm", path.name)]

        cds_path = pick_single_file(cds_matches, "phycocosm", species_key, "CDS", warnings)
        gff_path = pick_single_file(gff_matches, "phycocosm", species_key, "GFF", warnings)
        gbff_path = pick_single_file(gbff_matches, "phycocosm", species_key, "GBFF", warnings)
        genome_path = pick_single_file(genome_matches, "phycocosm", species_key, "genome", warnings)

        missing_label = task_missing_annotation_label(cds_path, gff_path, gbff_path, genome_path)
        if missing_label != "":
            errors.append("[phycocosm] {}: missing {}".format(species_key, missing_label))
            continue
        tasks.append(
            {
                "provider": "phycocosm",
                "species_key": species_key,
                "species_prefix": species_prefix,
                "cds_path": cds_path,
                "gff_path": gff_path,
                "gbff_path": gbff_path,
                "genome_path": genome_path,
            }
        )

    return tasks, warnings, errors


def discover_phytozome_tasks(input_dir, allowed_species_keys=None):
    warnings = []
    errors = []
    tasks = []

    for species_dir in sorted(input_dir.iterdir()):
        if not species_dir.is_dir():
            continue
        species_key = normalize_species_key_for_runtime(species_dir.name)
        if allowed_species_keys is not None and species_key not in allowed_species_keys:
            continue
        invalid_species_key = invalid_species_key_error("phytozome", species_key)
        if invalid_species_key != "":
            errors.append(invalid_species_key)
            continue
        species_prefix = species_prefix_from_value(species_key)
        if species_prefix == "":
            warnings.append("[phytozome] skipped '{}': unable to parse species prefix.".format(species_key))
            continue

        files = [path for path in species_dir.iterdir() if path.is_file()]
        cds_matches = [
            path
            for path in files
            if ("cds_" in path.name.lower() or ".cds." in path.name.lower()) and is_fasta_filename(path.name)
        ]
        gff_matches = [path for path in files if "gene.gff3" in path.name.lower() and is_gff_filename(path.name)]
        gbff_matches = [path for path in files if is_gbff_filename(path.name)]
        genome_matches = [path for path in files if is_probable_genome_filename("phytozome", path.name)]

        cds_path = pick_single_file(cds_matches, "phytozome", species_key, "CDS", warnings)
        gff_path = pick_single_file(gff_matches, "phytozome", species_key, "GFF", warnings)
        gbff_path = pick_single_file(gbff_matches, "phytozome", species_key, "GBFF", warnings)
        genome_path = pick_single_file(genome_matches, "phytozome", species_key, "genome", warnings)

        missing_label = task_missing_annotation_label(cds_path, gff_path, gbff_path, genome_path)
        if missing_label != "":
            errors.append("[phytozome] {}: missing {}".format(species_key, missing_label))
            continue
        tasks.append(
            {
                "provider": "phytozome",
                "species_key": species_key,
                "species_prefix": species_prefix,
                "cds_path": cds_path,
                "gff_path": gff_path,
                "gbff_path": gbff_path,
                "genome_path": genome_path,
            }
        )

    return tasks, warnings, errors


def discover_ncbi_like_tasks(input_dir, provider, allowed_species_keys=None):
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
        cds_matches = [
            path for path in files if "cds_from_genomic" in path.name.lower() and is_fasta_filename(path.name)
        ]
        gff_matches = [path for path in files if "genomic.gff" in path.name.lower() and is_gff_filename(path.name)]
        gbff_matches = [path for path in files if "genomic.gbff" in path.name.lower() and is_gbff_filename(path.name)]
        genome_matches = [path for path in files if is_probable_genome_filename(provider, path.name)]

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


def discover_tasks(provider, input_dir, allowed_species_keys=None):
    if provider in ENSEMBL_LIKE_PROVIDERS:
        return discover_ensembl_like_tasks(input_dir, provider, allowed_species_keys)
    if provider == "phycocosm":
        return discover_phycocosm_tasks(input_dir, allowed_species_keys)
    if provider == "phytozome":
        return discover_phytozome_tasks(input_dir, allowed_species_keys)
    if provider in ("ncbi", "refseq", "genbank", "plantgarden", "plantaedb"):
        return discover_ncbi_like_tasks(input_dir, provider, allowed_species_keys)
    if provider == "coge":
        return discover_generic_species_dir_tasks(provider, input_dir, allowed_species_keys)
    if provider == "cngb":
        return discover_generic_species_dir_tasks(provider, input_dir, allowed_species_keys)
    if provider in (
        "ddbj",
        "gwh",
        "citrusgenomedb",
        "figshare",
        "plantgarden",
        "flybase",
        "wormbase",
        "vectorbase",
        "fernbase",
        "veupathdb",
        "dictybase",
        "insectbase",
        ORYZA_MINUTA_PROVIDER,
        "direct",
    ):
        return discover_generic_species_dir_tasks(provider, input_dir, allowed_species_keys)
    if provider == "local":
        return discover_generic_species_dir_tasks(provider, input_dir, allowed_species_keys)
    raise ValueError("Unknown provider: {}".format(provider))


def build_formatted_cds_id(task, header):
    extracted = extract_provider_transcript_id(task["provider"], header)
    sanitized = sanitize_identifier(extracted)
    prefixed = "{}_{}".format(task["species_prefix"], sanitized)
    return sanitize_identifier(prefixed)


def prepare_cds_identifier_task(task):
    if task.get("provider") != "coge" or task.get("gff_path") is None:
        return task
    prepared = dict(task)
    prepared["_provider_gene_id_map"] = build_coge_gff_gene_id_map(task["gff_path"])
    return prepared


def format_cds(task, output_dir, overwrite, dry_run):
    cds_input = task.get("cds_path")
    if cds_input is not None:
        output_name = normalize_cds_output_basename(cds_input.name, task["species_prefix"])
    else:
        output_name = build_derived_cds_output_basename(task)
    output_path = output_dir / output_name

    if output_path.exists() and output_path.stat().st_size > 0 and not overwrite:
        before_count = 0
        for _header, _sequence in iter_task_cds_records(task):
            before_count += 1
        after_count, first_existing = count_fasta_records(output_path)
        return {
            "status": "skip",
            "output_path": output_path,
            "input_path": describe_task_cds_input(task),
            "written": after_count,
            "duplicates": max(0, before_count - after_count),
            "before_count": before_count,
            "after_count": after_count,
            "first_sequence_name": first_existing,
        }

    before_count = 0
    aggregated_away = 0
    first_sequence_name = ""
    records_by_gene = {}
    cds_task = prepare_cds_identifier_task(task)
    for header, sequence in iter_task_cds_records(task):
        before_count += 1
        transcript_id = build_formatted_cds_id(cds_task, header)
        gene_id = build_gene_aggregate_id(cds_task, header, transcript_id)
        seq = re.sub(r"\s+", "", sequence).upper()
        # Keep codon-frame-safe length (equivalent role of `cdskit pad` in shell pipelines).
        seq = pad_to_codon_length(seq)

        previous = records_by_gene.get(gene_id)
        if previous is None:
            records_by_gene[gene_id] = {
                "id": gene_id,
                "sequence": seq,
                "transcript_id": transcript_id,
            }
            continue

        previous_seq = previous["sequence"]
        previous_transcript_id = previous["transcript_id"]
        # Keep one representative CDS per gene; prefer longer CDS and then lexicographically stable tie-breaker.
        if len(seq) > len(previous_seq) or (len(seq) == len(previous_seq) and transcript_id < previous_transcript_id):
            records_by_gene[gene_id] = {
                "id": gene_id,
                "sequence": seq,
                "transcript_id": transcript_id,
            }

    ordered_ids = sorted(records_by_gene.keys())
    after_count = len(ordered_ids)
    aggregated_away = max(0, before_count - after_count)
    if len(ordered_ids) > 0:
        first_sequence_name = ordered_ids[0]

    if task.get("cds_path") is None and after_count == 0:
        if output_path.exists():
            output_path.unlink()
        return {
            "status": "empty",
            "output_path": None,
            "input_path": describe_task_cds_input(task),
            "written": 0,
            "duplicates": aggregated_away,
            "before_count": before_count,
            "after_count": after_count,
            "first_sequence_name": first_sequence_name,
        }

    if not dry_run:
        write_fasta_records_gzip(
            output_path,
            ((records_by_gene[gene_id]["id"], records_by_gene[gene_id]["sequence"]) for gene_id in ordered_ids),
        )

    status = "dry-run" if dry_run else "write"
    return {
        "status": status,
        "output_path": output_path,
        "input_path": describe_task_cds_input(task),
        "written": after_count,
        "duplicates": aggregated_away,
        "before_count": before_count,
        "after_count": after_count,
        "first_sequence_name": first_sequence_name,
    }


def format_genome(task, output_dir, overwrite, dry_run):
    genome_path = task.get("genome_path")
    gbff_path = task.get("gbff_path")
    if genome_path is None and gbff_path is None:
        return {"status": "missing", "output_path": None, "written": 0}

    if genome_path is not None:
        output_name = normalize_genome_output_basename(genome_path.name, task["species_prefix"])
    else:
        output_name = build_derived_genome_output_basename(task)
    output_path = output_dir / output_name
    if output_path.exists() and output_path.stat().st_size > 0 and not overwrite:
        return {"status": "skip", "output_path": output_path, "written": 0}
    if dry_run:
        return {"status": "dry-run", "output_path": output_path, "written": 0}

    written = 0

    def iter_genome_output_records():
        nonlocal written
        if genome_path is not None:
            for header, sequence in iter_fasta_records(genome_path):
                record_id = first_token(apply_common_replacements(header))
                if record_id == "":
                    record_id = "unnamed"
                seq = re.sub(r"\s+", "", sequence).upper()
                written += 1
                yield record_id, seq
            return
        for record_id, sequence in iter_genome_records_from_gbff(gbff_path):
            written += 1
            yield record_id, sequence

    write_fasta_records_gzip(output_path, iter_genome_output_records())
    return {"status": "write", "output_path": output_path, "written": written}


def format_gff(task, output_dir, overwrite, dry_run, formatted_cds_path=None):
    gff_path = task.get("gff_path")
    gbff_path = task.get("gbff_path")
    repair_mode = gff_repair_mode_for_task(task)
    if gff_path is not None:
        output_name = normalize_gff_output_basename(gff_path.name, task["species_prefix"])
    elif gbff_path is not None:
        output_name = build_derived_gff_output_basename(task)
    else:
        return {
            "status": "missing",
            "output_path": None,
            "lines": 0,
            "repair_mode": repair_mode,
            "repair_status": "not_applicable",
            "repair_audit_path": None,
            "repair_gene_ids": 0,
            "repair_references": 0,
            "repair_ambiguous": 0,
            "repair_collisions": 0,
        }
    output_path = output_dir / output_name
    if output_path.exists() and output_path.stat().st_size > 0 and not overwrite:
        if gff_path is not None and formatted_cds_path is not None and Path(formatted_cds_path).exists():
            audit_path = gff_repair_audit_path(output_path)
            audit = read_json(audit_path)
            if audit_matches_inputs(
                audit,
                repair_mode,
                gff_path,
                formatted_cds_path,
                output_path,
            ):
                result = {"status": "skip", "output_path": output_path, "lines": 0}
                result.update(repair_result_fields(audit, output_path))
                return result
            if audit is None and repair_mode == "off":
                result = {"status": "skip", "output_path": output_path, "lines": 0}
                result.update(repair_result_fields(None, output_path))
                result["repair_mode"] = repair_mode
                result["repair_status"] = "legacy_untracked"
                return result
        elif gff_path is None or formatted_cds_path is None:
            result = {"status": "skip", "output_path": output_path, "lines": 0}
            result.update(repair_result_fields(None, output_path))
            result["repair_mode"] = repair_mode
            result["repair_status"] = "not_applied"
            return result
    if dry_run:
        result = {"status": "dry-run", "output_path": output_path, "lines": 0}
        result.update(repair_result_fields(None, output_path))
        result["repair_mode"] = repair_mode
        result["repair_status"] = "planned"
        return result

    feature_count = 0
    if gff_path is not None:
        if formatted_cds_path is not None and Path(formatted_cds_path).exists():
            audit = write_repaired_gff(
                gff_path=gff_path,
                cds_path=formatted_cds_path,
                output_path=output_path,
                species_prefix=task["species_prefix"],
                mode=repair_mode,
            )
            line_count = int(audit.get("line_count", 0) or 0)
            repair_fields = repair_result_fields(audit, output_path)
        else:
            line_count = write_gff_gzip(gff_path, output_path)
            repair_fields = repair_result_fields(None, output_path)
            repair_fields["repair_mode"] = repair_mode
            repair_fields["repair_status"] = "not_applied"
    else:
        line_count, feature_count = write_gff_lines_gzip(output_path, iter_gff_lines_from_gbff(task))
        repair_fields = repair_result_fields(None, output_path)
        repair_fields["repair_mode"] = repair_mode
        repair_fields["repair_status"] = "not_applicable"
        if feature_count == 0:
            if output_path.exists():
                output_path.unlink()
            result = {"status": "empty", "output_path": None, "lines": line_count}
            result.update(repair_fields)
            return result
    result = {"status": "write", "output_path": output_path, "lines": line_count}
    result.update(repair_fields)
    return result


def utc_now_iso():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
