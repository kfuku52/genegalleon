"""Local input discovery and formatted output writers."""

import csv
import hashlib
import json
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
    build_gff_cds_grouping_index,
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
    resolve_gene_aggregate_id,
    sanitize_identifier,
    species_prefix_from_value,
    task_missing_annotation_label,
    write_repaired_gff,
)
from format_species_common import (
    is_fasta_filename,
    is_gbff_filename,
    is_gff_filename,
    is_probable_cds_filename,
    is_probable_genome_filename,
    pick_single_file,
)
from format_species_constants import (
    ENSEMBL_CDS_PATTERN,
    KNOWN_ALLOWED_MISSING_CDS_IDS,
    gff_mapping_fallback_is_tolerable,
)
from format_species_provider_config import ENSEMBL_LIKE_PROVIDERS, ORYZA_MINUTA_PROVIDER
from format_species_taxonomy import invalid_species_key_error, normalize_species_key_for_runtime
from format_species_writers import (
    apply_common_replacements,
    write_fasta_records_gzip,
    write_gff_gzip,
    write_gff_lines_gzip,
)

CDS_GFF_GROUPING_AUDIT_VERSION = 7


def cds_gff_grouping_audit_paths(output_path):
    return (
        Path(str(output_path) + ".gff-grouping.json"),
        Path(str(output_path) + ".gff-grouping.tsv"),
    )


def cds_gff_source_signature(path):
    resolved = Path(path).expanduser().resolve()
    before_stat = resolved.stat()
    digest = hashlib.sha256()
    with open(resolved, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    after_stat = resolved.stat()
    if (
        int(before_stat.st_size) != int(after_stat.st_size)
        or int(before_stat.st_mtime_ns) != int(after_stat.st_mtime_ns)
    ):
        raise OSError("File changed while computing its signature: {}".format(resolved))
    return {
        "path": str(resolved),
        "size": int(after_stat.st_size),
        "sha256": digest.hexdigest(),
    }


def cds_gff_artifact_signature(path):
    return cds_gff_source_signature(path)


def cds_gff_signatures_equal(recorded, current):
    """Compare content contracts while ignoring legacy mtime diagnostics."""

    if not isinstance(recorded, dict) or not isinstance(current, dict):
        return False
    keys = ("path", "size", "sha256")
    return all(recorded.get(key) == current.get(key) for key in keys)


def cds_gff_grouping_audit_matches(audit, task, output_path, strict_mode):
    if not isinstance(audit, dict):
        return False
    if int(audit.get("version", 0) or 0) != CDS_GFF_GROUPING_AUDIT_VERSION:
        return False
    if str(audit.get("gene_grouping_mode", "") or "") != str(task.get("gene_grouping_mode", "strict") or "strict"):
        return False
    if bool(audit.get("format_strict", False)) != bool(strict_mode):
        return False
    if task.get("cds_path") is None or task.get("gff_path") is None:
        return False
    try:
        if not cds_gff_signatures_equal(
            audit.get("cds_input"), cds_gff_source_signature(task["cds_path"])
        ):
            return False
        if not cds_gff_signatures_equal(
            audit.get("gff_input"), cds_gff_source_signature(task["gff_path"])
        ):
            return False
        if not cds_gff_signatures_equal(
            audit.get("output_fingerprint"), cds_gff_artifact_signature(output_path)
        ):
            return False
    except (FileNotFoundError, OSError):
        return False
    if str(audit.get("output_path", "") or "") != str(Path(output_path).expanduser().resolve()):
        return False
    _json_path, records_path = cds_gff_grouping_audit_paths(output_path)
    try:
        return cds_gff_signatures_equal(
            audit.get("records_audit_fingerprint"), cds_gff_artifact_signature(records_path)
        )
    except (FileNotFoundError, OSError):
        return False


def cds_gff_result_fields(audit=None):
    payload = audit if isinstance(audit, dict) else {}
    stats = payload.get("stats") if isinstance(payload.get("stats"), dict) else {}
    return {
        "grouping_source": str(payload.get("grouping_source", "header") or "header"),
        "gff_grouping_audit_path": str(payload.get("records_audit_path", "") or ""),
        "gff_records_mapped": int(stats.get("mapped", 0) or 0),
        "gff_records_unmapped": int(stats.get("unmapped", 0) or 0),
        "gff_records_ambiguous": int(stats.get("ambiguous", 0) or 0),
        "gff_unexpected_mapping_records": int(stats.get("unexpected_mapping_records", 0) or 0),
        "gff_mapping_fallback_tolerated": int(stats.get("mapping_fallback_tolerated", 0) or 0),
        "gff_coordinate_rescued_transcripts": int(stats.get("coordinate_rescued_transcripts", 0) or 0),
        "gff_coordinate_rescued_groups": int(stats.get("coordinate_rescued_groups", 0) or 0),
    }


def write_cds_gff_grouping_audit(task, output_path, audit_rows, payload, strict_mode):
    json_path, records_path = cds_gff_grouping_audit_paths(output_path)
    json_tmp = Path(str(json_path) + ".tmp.{}.{}".format(time.time_ns(), len(audit_rows)))
    records_tmp = Path(str(records_path) + ".tmp.{}.{}".format(time.time_ns(), len(audit_rows)))
    fieldnames = (
        "record_index",
        "raw_cds_id",
        "formatted_transcript_id",
        "mapping_status",
        "matched_aliases",
        "candidate_gene_tokens",
        "selected_gene_id",
        "raw_sequence_length",
        "sequence_length",
        "selected_longest",
    )
    try:
        with open(records_tmp, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for row in audit_rows:
                writer.writerow({key: row.get(key, "") for key in fieldnames})
        records_tmp.replace(records_path)
        final_payload = dict(payload)
        final_payload.update(
            {
                "version": CDS_GFF_GROUPING_AUDIT_VERSION,
                "provider": task["provider"],
                "species_key": task["species_key"],
                "species_prefix": task["species_prefix"],
                "gene_grouping_mode": task.get("gene_grouping_mode", "strict"),
                "format_strict": bool(strict_mode),
                "cds_input": cds_gff_source_signature(task["cds_path"]),
                "gff_input": cds_gff_source_signature(task["gff_path"]),
                "output_path": str(Path(output_path).expanduser().resolve()),
                "output_fingerprint": cds_gff_artifact_signature(output_path),
                "records_audit_path": str(records_path.resolve()),
                "records_audit_fingerprint": cds_gff_artifact_signature(records_path),
            }
        )
        with open(json_tmp, "wt", encoding="utf-8") as handle:
            json.dump(final_payload, handle, ensure_ascii=True, indent=2, sort_keys=True)
        json_tmp.replace(json_path)
        return final_payload
    finally:
        for tmp_path in (json_tmp, records_tmp):
            try:
                tmp_path.unlink()
            except FileNotFoundError:
                pass


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
        cds_matches = [
            path
            for path in files
            if is_probable_cds_filename("phycocosm", path.name)
            and not is_probable_genome_filename("phycocosm", path.name)
        ]
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

        all_matches = cds_matches + gff_matches + gbff_matches + genome_matches
        accessions = sorted(
            {
                match.group(1).upper()
                for path in all_matches
                for match in [re.search(r"(GC[AF]_[0-9]+\.[0-9]+)", path.name, re.IGNORECASE)]
                if match is not None
            }
        )
        selected_accession = ""
        if len(accessions) > 0:
            candidate_bundles = []
            for accession in accessions:
                bundle_cds = [path for path in cds_matches if accession.lower() in path.name.lower()]
                bundle_gff = [path for path in gff_matches if accession.lower() in path.name.lower()]
                bundle_gbff = [path for path in gbff_matches if accession.lower() in path.name.lower()]
                bundle_genome = [path for path in genome_matches if accession.lower() in path.name.lower()]
                coherent = len(bundle_gbff) > 0 or len(bundle_gff) > 0 and (
                    len(bundle_cds) > 0 or len(bundle_genome) > 0
                )
                if coherent:
                    candidate_bundles.append(
                        (
                            -(len(bundle_cds) > 0) - (len(bundle_gff) > 0) - (len(bundle_genome) > 0),
                            accession,
                            bundle_cds,
                            bundle_gff,
                            bundle_gbff,
                            bundle_genome,
                        )
                    )
            if len(candidate_bundles) == 0 and len(accessions) > 1:
                errors.append(
                    "[{}] {}: annotation files span incompatible assembly accessions: {}".format(
                        provider, species_key, ",".join(accessions)
                    )
                )
                continue
            if len(candidate_bundles) > 0:
                _score, selected_accession, cds_matches, gff_matches, gbff_matches, genome_matches = sorted(
                    candidate_bundles, key=lambda item: (item[0], item[1])
                )[0]
                if len(accessions) > 1:
                    warnings.append(
                        "[{}] {}: multiple assembly bundles found. Using coherent bundle '{}'.".format(
                            provider, species_key, selected_accession
                        )
                    )

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
                "source_bundle_id": selected_accession,
                "gff_auto_selected_from_multiple": len(gff_matches) > 1,
                "gff_selection_candidates": tuple(sorted(path.name for path in gff_matches)),
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
    if provider in ("ncbi", "refseq", "genbank", "plantaedb"):
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
    prepared = dict(task)
    if task.get("cds_path") is not None and task.get("gff_path") is not None:
        prepared["_gff_cds_grouping_index"] = build_gff_cds_grouping_index(task)
    if task.get("provider") == "coge" and task.get("gff_path") is not None:
        prepared["_provider_gene_id_map"] = build_coge_gff_gene_id_map(task["gff_path"])
    return prepared


def format_cds(task, output_dir, overwrite, dry_run, strict=None, reuse_existing=False):
    cds_input = task.get("cds_path")
    if cds_input is not None:
        output_name = normalize_cds_output_basename(cds_input.name, task["species_prefix"])
    else:
        output_name = build_derived_cds_output_basename(task)
    output_path = output_dir / output_name

    use_gff_grouping = task.get("cds_path") is not None and task.get("gff_path") is not None
    strict_mode = bool(task.get("format_strict", False)) if strict is None else bool(strict)
    audit_json_path, _audit_records_path = cds_gff_grouping_audit_paths(output_path)
    existing_audit = read_json(audit_json_path) if use_gff_grouping else None
    # gg-cache-guard: audited - the grouping audit covers inputs/parameters; the outer formatter contract handles rebuild policy.
    if (
        output_path.exists()
        and output_path.stat().st_size > 0
        and not overwrite
        and (
            reuse_existing
            or
            not use_gff_grouping
            or cds_gff_grouping_audit_matches(existing_audit, task, output_path, strict_mode)
        )
    ):
        if use_gff_grouping:
            result = {
                "status": "skip",
                "output_path": output_path,
                "input_path": describe_task_cds_input(task),
                "written": int(existing_audit.get("after_count", 0) or 0),
                "duplicates": int(existing_audit.get("duplicates", 0) or 0),
                "before_count": int(existing_audit.get("before_count", 0) or 0),
                "after_count": int(existing_audit.get("after_count", 0) or 0),
                "first_sequence_name": str(existing_audit.get("first_sequence_name", "") or ""),
            }
            result.update(cds_gff_result_fields(existing_audit))
            return result
        before_count = 0
        for _header, _sequence in iter_task_cds_records(task):
            before_count += 1
        after_count, first_existing = count_fasta_records(output_path)
        result = {
            "status": "skip",
            "output_path": output_path,
            "input_path": describe_task_cds_input(task),
            "written": after_count,
            "duplicates": max(0, before_count - after_count),
            "before_count": before_count,
            "after_count": after_count,
            "first_sequence_name": first_existing,
        }
        result.update(cds_gff_result_fields())
        return result

    before_count = 0
    aggregated_away = 0
    first_sequence_name = ""
    records_by_gene = {}
    cds_task = prepare_cds_identifier_task(task)
    grouping_index = cds_task.get("_gff_cds_grouping_index")
    mapping_counts = {"mapped": 0, "unmapped": 0, "ambiguous": 0}
    raw_gff_tokens_by_gene_id = defaultdict(set)
    audit_rows = []
    for header, sequence in iter_task_cds_records(task):
        before_count += 1
        transcript_id = build_formatted_cds_id(cds_task, header)
        gene_id, gff_match = resolve_gene_aggregate_id(cds_task, header, transcript_id)
        mapping_status = str(gff_match.get("status", "not_applicable") or "not_applicable")
        if mapping_status in mapping_counts:
            mapping_counts[mapping_status] += 1
        if mapping_status == "mapped":
            raw_gene_token = str(gff_match.get("gene_token", "") or "").strip()
            if raw_gene_token != "":
                raw_gff_tokens_by_gene_id[gene_id].add(raw_gene_token)
        raw_seq = re.sub(r"\s+", "", sequence).upper()
        # Keep codon-frame-safe length (equivalent role of `cdskit pad` in shell pipelines).
        seq = pad_to_codon_length(raw_seq)
        record_index = before_count
        audit_rows.append(
            {
                "record_index": record_index,
                "raw_cds_id": first_token(header),
                "formatted_transcript_id": transcript_id,
                "mapping_status": mapping_status,
                "matched_aliases": ",".join(gff_match.get("matched_aliases", ())),
                "candidate_gene_tokens": ",".join(gff_match.get("candidate_gene_tokens", ())),
                "selected_gene_id": gene_id,
                "raw_sequence_length": len(raw_seq),
                "sequence_length": len(seq),
                "selected_longest": 0,
            }
        )

        previous = records_by_gene.get(gene_id)
        if previous is None:
            records_by_gene[gene_id] = {
                "id": gene_id,
                "sequence": seq,
                "raw_sequence_length": len(raw_seq),
                "transcript_id": transcript_id,
                "record_index": record_index,
            }
            continue

        previous_raw_sequence_length = previous["raw_sequence_length"]
        previous_transcript_id = previous["transcript_id"]
        # Compare biological input lengths before codon padding; padded ties can otherwise hide a longer CDS.
        if len(raw_seq) > previous_raw_sequence_length or (
            len(raw_seq) == previous_raw_sequence_length and transcript_id < previous_transcript_id
        ):
            records_by_gene[gene_id] = {
                "id": gene_id,
                "sequence": seq,
                "raw_sequence_length": len(raw_seq),
                "transcript_id": transcript_id,
                "record_index": record_index,
            }

    ordered_ids = sorted(records_by_gene.keys())
    after_count = len(ordered_ids)
    aggregated_away = max(0, before_count - after_count)
    if len(ordered_ids) > 0:
        first_sequence_name = ordered_ids[0]
    selected_record_indexes = {records_by_gene[gene_id]["record_index"] for gene_id in ordered_ids}
    for row in audit_rows:
        row["selected_longest"] = int(row["record_index"] in selected_record_indexes)

    sanitized_gene_id_collisions = {
        gene_id: tuple(sorted(raw_tokens))
        for gene_id, raw_tokens in raw_gff_tokens_by_gene_id.items()
        if len(raw_tokens) > 1
    }
    if len(sanitized_gene_id_collisions) > 0:
        examples = []
        for gene_id in sorted(sanitized_gene_id_collisions)[:5]:
            examples.append("{} <- {}".format(gene_id, ",".join(sanitized_gene_id_collisions[gene_id])))
        raise ValueError(
            "GFF-backed CDS grouping for {} has distinct gene IDs that collide after identifier "
            "sanitization: {}".format(task.get("species_prefix", ""), "; ".join(examples))
        )

    allowed_missing_ids = KNOWN_ALLOWED_MISSING_CDS_IDS.get(task.get("species_prefix", ""), frozenset())
    unexpected_unmapped_ids = sorted(
        {
            row["selected_gene_id"]
            for row in audit_rows
            if row["mapping_status"] == "unmapped" and row["selected_gene_id"] not in allowed_missing_ids
        }
    )
    unexpected_mapping_count = len(unexpected_unmapped_ids) + mapping_counts["ambiguous"]
    mapping_fallback_tolerated = (
        grouping_index is not None
        and unexpected_mapping_count > 0
        and not strict_mode
        and gff_mapping_fallback_is_tolerable(
            before_count,
            len(unexpected_unmapped_ids),
            mapping_counts["ambiguous"],
        )
    )
    if grouping_index is not None and unexpected_mapping_count > 0 and not mapping_fallback_tolerated:
        candidate_names = ",".join(task.get("gff_selection_candidates", ())) or Path(task["gff_path"]).name
        raise ValueError(
            "GFF-backed CDS grouping for {} failed: unexpected_unmapped={} ambiguous={} "
            "GFF='{}' candidates={} sample={}".format(
                task.get("species_prefix", ""),
                len(unexpected_unmapped_ids),
                mapping_counts["ambiguous"],
                Path(task["gff_path"]).name,
                candidate_names,
                ",".join(unexpected_unmapped_ids[:5]),
            )
        )

    if task.get("cds_path") is None and after_count == 0:
        if output_path.exists():
            output_path.unlink()
        result = {
            "status": "empty",
            "output_path": None,
            "input_path": describe_task_cds_input(task),
            "written": 0,
            "duplicates": aggregated_away,
            "before_count": before_count,
            "after_count": after_count,
            "first_sequence_name": first_sequence_name,
        }
        result.update(cds_gff_result_fields())
        return result

    if not dry_run:
        write_fasta_records_gzip(
            output_path,
            ((records_by_gene[gene_id]["id"], records_by_gene[gene_id]["sequence"]) for gene_id in ordered_ids),
        )

    if grouping_index is None:
        grouping_source = "header"
    elif mapping_counts["unmapped"] > 0 or mapping_counts["ambiguous"] > 0:
        grouping_source = "gff_with_allowed_missing"
    else:
        grouping_source = "gff"
    audit_payload = {
        "grouping_source": grouping_source,
        "before_count": before_count,
        "after_count": after_count,
        "duplicates": aggregated_away,
        "first_sequence_name": first_sequence_name,
        "stats": {
            "mapped": mapping_counts["mapped"],
            "unmapped": mapping_counts["unmapped"],
            "ambiguous": mapping_counts["ambiguous"],
            "unexpected_mapping_records": unexpected_mapping_count,
            "mapping_fallback_tolerated": int(mapping_fallback_tolerated),
            "gff_transcripts_total": int((grouping_index or {}).get("transcripts_total", 0) or 0),
            "coordinate_rescued_transcripts": int(
                (grouping_index or {}).get("coordinate_rescued_transcripts", 0) or 0
            ),
            "coordinate_rescued_groups": int((grouping_index or {}).get("coordinate_rescued_groups", 0) or 0),
        },
    }
    if grouping_index is not None and not dry_run:
        audit_payload = write_cds_gff_grouping_audit(
            task,
            output_path,
            audit_rows,
            audit_payload,
            strict_mode,
        )

    status = "dry-run" if dry_run else "write"
    result = {
        "status": status,
        "output_path": output_path,
        "input_path": describe_task_cds_input(task),
        "written": after_count,
        "duplicates": aggregated_away,
        "before_count": before_count,
        "after_count": after_count,
        "first_sequence_name": first_sequence_name,
    }
    result.update(cds_gff_result_fields(audit_payload))
    return result


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
    # gg-cache-guard: audited - the outer formatter contract deletes this species output on rebuild.
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


def format_gff(
    task,
    output_dir,
    overwrite,
    dry_run,
    formatted_cds_path=None,
    reuse_existing=False,
):
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
    # gg-cache-guard: audited - reuse is explicit or input/output hashes are checked below; outer provenance handles rebuild policy.
    if output_path.exists() and output_path.stat().st_size > 0 and not overwrite:
        if reuse_existing:
            result = {"status": "skip", "output_path": output_path, "lines": 0}
            result.update(repair_result_fields(None, output_path))
            result["repair_mode"] = repair_mode
            result["repair_status"] = "reused_by_policy"
            return result
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
