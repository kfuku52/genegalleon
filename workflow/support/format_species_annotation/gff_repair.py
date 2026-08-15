"""Conservative repair of formatted GFF identifiers against final CDS IDs."""

import json
import os
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote

from format_species_writers import (
    apply_common_replacements,
    inspect_invalid_utf8,
    open_text,
    write_gff_lines_gzip,
)

from .common import (
    first_token,
    iter_fasta_records,
    normalize_gff_attribute_value,
    parse_gff_attributes,
    sanitize_identifier,
)

GFF_REPAIR_VERSION = 2
GFF_REPAIR_MODES = ("off", "safe", "strict")
GENE_ALIAS_KEYS = ("Name", "Alias", "gene", "gene_id", "locus_tag", "geneName", "ID")
GENE_REFERENCE_KEYS = frozenset(("Parent", "Derives_from", "gene", "gene_id"))


def normalize_gff_repair_mode(value):
    text = str(value or "").strip().lower()
    if text == "":
        return "safe"
    if text not in GFF_REPAIR_MODES:
        raise ValueError(
            "Unknown gff_repair_mode: {} (allowed: {})".format(
                text,
                ",".join(GFF_REPAIR_MODES),
            )
        )
    return text


def gff_repair_mode_for_task(task):
    return normalize_gff_repair_mode(task.get("gff_repair_mode", "safe"))


def gff_repair_audit_path(output_path):
    return Path(str(output_path) + ".repair.json")


def file_fingerprint(path):
    target = Path(path).expanduser().resolve()
    stat = target.stat()
    return {
        "path": str(target),
        "size": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def read_json(path):
    try:
        with open(path, "rt", encoding="utf-8") as handle:
            payload = json.load(handle)
    except (FileNotFoundError, OSError, ValueError, TypeError):
        return None
    return payload if isinstance(payload, dict) else None


def write_json_atomic(path, payload):
    path = Path(path)
    tmp_path = Path(str(path) + ".tmp.{}".format(os.getpid()))
    try:
        with open(tmp_path, "wt", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=True, indent=2, sort_keys=True)
            handle.write("\n")
        tmp_path.replace(path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def audit_matches_inputs(audit, mode, source_path, cds_path, output_path):
    if not isinstance(audit, dict):
        return False
    if int(audit.get("repair_version", 0) or 0) != GFF_REPAIR_VERSION:
        return False
    if str(audit.get("mode") or "") != mode:
        return False
    for key, path in (
        ("source_fingerprint", source_path),
        ("cds_fingerprint", cds_path),
        ("output_fingerprint", output_path),
    ):
        try:
            current = file_fingerprint(path)
        except (FileNotFoundError, OSError):
            return False
        if audit.get(key) != current:
            return False
    return True


def read_formatted_cds_gene_ids(cds_path, species_prefix):
    prefix = "{}_".format(species_prefix)
    out = set()
    for header, _sequence in iter_fasta_records(Path(cds_path)):
        record_id = first_token(header)
        if record_id.startswith(prefix):
            gene_id = record_id[len(prefix) :]
        else:
            gene_id = record_id
        if gene_id != "":
            out.add(gene_id)
    return out


def iter_gff_feature_rows(gff_path):
    with open_text(Path(gff_path), "rt", errors="replace") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = apply_common_replacements(raw_line.rstrip("\n\r"))
            if line == "" or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            yield line_number, parts[2].strip().lower(), parse_gff_attributes(parts[8])


def choose_gene_id_repairs(gff_path, cds_gene_ids):
    gene_attrs_by_id = defaultdict(lambda: defaultdict(set))
    id_feature_types = defaultdict(set)
    missing_gene_id_lines = []

    for line_number, feature_type, attrs in iter_gff_feature_rows(gff_path):
        feature_ids = tuple(str(value or "").strip() for value in attrs.get("ID", ()))
        for feature_id in feature_ids:
            if feature_id != "":
                id_feature_types[feature_id].add(feature_type)
        if feature_type != "gene":
            continue
        raw_gene_id = feature_ids[0] if len(feature_ids) > 0 else ""
        if raw_gene_id == "":
            missing_gene_id_lines.append(line_number)
            continue
        for key in GENE_ALIAS_KEYS:
            for raw_value in attrs.get(key, ()):
                value = str(raw_value or "").strip()
                if value != "":
                    gene_attrs_by_id[raw_gene_id][key].add(value)

    proposed = {}
    evidence_by_id = {}
    ambiguous = []
    for raw_gene_id in sorted(gene_attrs_by_id):
        hits = defaultdict(list)
        for key in GENE_ALIAS_KEYS:
            for value in sorted(gene_attrs_by_id[raw_gene_id].get(key, ())):
                if value in cds_gene_ids:
                    hits[value].append({"attribute": key, "value": value, "match": "exact"})
                normalized = sanitize_identifier(value)
                if normalized in cds_gene_ids and normalized != value:
                    hits[normalized].append(
                        {"attribute": key, "value": value, "match": "sanitized"}
                    )
        if len(hits) == 0:
            continue
        if len(hits) > 1:
            ambiguous.append(
                {
                    "source_id": raw_gene_id,
                    "candidate_cds_ids": sorted(hits),
                }
            )
            continue
        target_id = next(iter(hits))
        proposed[raw_gene_id] = target_id
        evidence_by_id[raw_gene_id] = sorted(
            hits[target_id],
            key=lambda item: (
                0 if item["match"] == "exact" else 1,
                GENE_ALIAS_KEYS.index(item["attribute"]),
                item["value"],
            ),
        )[0]

    target_sources = defaultdict(list)
    for source_id, target_id in proposed.items():
        target_sources[target_id].append(source_id)

    collisions = []
    blocked = set()
    for target_id, source_ids in sorted(target_sources.items()):
        if len(source_ids) > 1:
            collisions.append(
                {
                    "target_id": target_id,
                    "source_ids": sorted(source_ids),
                    "reason": "multiple_gene_features_map_to_one_cds_id",
                }
            )
            blocked.update(source_ids)
            continue
        source_id = source_ids[0]
        if target_id != source_id and target_id in id_feature_types:
            collisions.append(
                {
                    "target_id": target_id,
                    "source_ids": [source_id],
                    "existing_feature_types": sorted(id_feature_types[target_id]),
                    "reason": "target_id_already_exists",
                }
            )
            blocked.add(source_id)
        if any(feature_type != "gene" for feature_type in id_feature_types.get(source_id, ())):
            collisions.append(
                {
                    "target_id": target_id,
                    "source_ids": [source_id],
                    "existing_feature_types": sorted(id_feature_types[source_id]),
                    "reason": "source_id_is_shared_with_non_gene_feature",
                }
            )
            blocked.add(source_id)

    accepted = {
        source_id: target_id
        for source_id, target_id in proposed.items()
        if source_id not in blocked and source_id != target_id
    }
    repairs = []
    for source_id, target_id in sorted(accepted.items()):
        evidence = evidence_by_id[source_id]
        repairs.append(
            {
                "source_id": source_id,
                "target_id": target_id,
                "attribute": evidence["attribute"],
                "source_value": evidence["value"],
                "reason": "{}_{}_to_cds_id".format(evidence["match"], evidence["attribute"]),
            }
        )

    return {
        "id_mapping": accepted,
        "repairs": repairs,
        "ambiguous": ambiguous,
        "collisions": collisions,
        "missing_gene_id_lines": missing_gene_id_lines,
        "gene_features": len(gene_attrs_by_id),
        "cds_gene_ids": len(cds_gene_ids),
    }


def split_attribute_field(field):
    equal_pos = field.find("=")
    space_pos = field.find(" ")
    if equal_pos != -1 and (space_pos == -1 or equal_pos < space_pos):
        return field[:equal_pos].strip(), field[equal_pos + 1 :], "="
    return "", "", ""


def normalize_raw_attribute_value(value):
    text = str(value or "").strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in ("'", '"'):
        text = text[1:-1]
    return normalize_gff_attribute_value(text)


def rewrite_attribute_value(raw_value, id_mapping, multiple=False):
    values = raw_value.split(",") if multiple else [raw_value]
    changed = 0
    rewritten = []
    for raw_item in values:
        normalized = normalize_raw_attribute_value(raw_item)
        replacement = id_mapping.get(normalized)
        if replacement is None:
            rewritten.append(raw_item)
            continue
        rewritten.append(replacement)
        changed += 1
    return ",".join(rewritten), changed


def rewrite_gff_attributes(attr_text, feature_type, id_mapping):
    bare_fields = [field.strip() for field in str(attr_text or "").split(";") if field.strip() != ""]
    if (
        len(bare_fields) == 1
        and bare_fields[0] != "."
        and "=" not in bare_fields[0]
        and not re.search(r"\s", bare_fields[0])
    ):
        token = normalize_gff_attribute_value(bare_fields[0])
        encoded = quote(token, safe="._:-|")
        if feature_type == "gene":
            return "ID={};gene_id={}".format(encoded, encoded), 0, 0, 1
        return "Parent={};gene_id={}".format(encoded, encoded), 0, 0, 1

    fields = str(attr_text).split(";")
    changed = 0
    reference_changes = 0
    rewritten = []
    for field in fields:
        key, raw_value, separator = split_attribute_field(field)
        if separator == "":
            rewritten.append(field)
            continue
        should_rewrite = key in GENE_REFERENCE_KEYS or (feature_type == "gene" and key == "ID")
        if not should_rewrite:
            rewritten.append(field)
            continue
        next_value, value_changes = rewrite_attribute_value(
            raw_value,
            id_mapping,
            multiple=key in ("Parent", "Derives_from"),
        )
        if value_changes == 0:
            rewritten.append(field)
            continue
        rewritten.append("{}{}{}".format(key, separator, next_value))
        changed += value_changes
        if key in GENE_REFERENCE_KEYS:
            reference_changes += value_changes
    return ";".join(rewritten), changed, reference_changes, 0


def iter_repaired_gff_lines(gff_path, id_mapping, counters):
    with open_text(Path(gff_path), "rt", errors="replace") as handle:
        for raw_line in handle:
            line = apply_common_replacements(raw_line)
            stripped = line.rstrip("\n\r")
            newline = line[len(stripped) :]
            if stripped == "" or stripped.startswith("#"):
                yield line
                continue
            parts = stripped.split("\t")
            if len(parts) < 9:
                yield line
                continue
            feature_type = parts[2].strip().lower()
            attributes, value_changes, reference_changes, normalized_bare = rewrite_gff_attributes(
                parts[8],
                feature_type,
                id_mapping,
            )
            if value_changes > 0 or normalized_bare > 0:
                parts[8] = attributes
                counters["changed_lines"] += 1
                counters["changed_values"] += value_changes
                counters["changed_references"] += reference_changes
                counters["normalized_bare_attribute_lines"] += normalized_bare
                line = "\t".join(parts) + newline
            yield line


def write_repaired_gff(gff_path, cds_path, output_path, species_prefix, mode):
    mode = normalize_gff_repair_mode(mode)
    source_fingerprint = file_fingerprint(gff_path)
    cds_fingerprint = file_fingerprint(cds_path)
    encoding_audit = inspect_invalid_utf8(gff_path)
    cds_gene_ids = read_formatted_cds_gene_ids(cds_path, species_prefix)
    plan = choose_gene_id_repairs(gff_path, cds_gene_ids) if mode != "off" else {
        "id_mapping": {},
        "repairs": [],
        "ambiguous": [],
        "collisions": [],
        "missing_gene_id_lines": [],
        "gene_features": 0,
        "cds_gene_ids": len(cds_gene_ids),
    }
    if mode == "strict" and (len(plan["ambiguous"]) > 0 or len(plan["collisions"]) > 0):
        raise ValueError(
            "GFF repair is ambiguous for {}: ambiguous={}, collisions={}".format(
                gff_path,
                len(plan["ambiguous"]),
                len(plan["collisions"]),
            )
        )

    counters = {
        "changed_lines": 0,
        "changed_values": 0,
        "changed_references": 0,
        "normalized_bare_attribute_lines": 0,
    }
    line_count, _feature_count = write_gff_lines_gzip(
        Path(output_path),
        iter_repaired_gff_lines(gff_path, plan["id_mapping"], counters),
    )
    status = (
        "repaired"
        if (
            len(plan["id_mapping"]) > 0
            or counters["normalized_bare_attribute_lines"] > 0
            or encoding_audit["invalid_utf8_bytes"] > 0
        )
        else "unchanged"
    )
    if mode == "off":
        status = "off"
    audit = {
        "repair_version": GFF_REPAIR_VERSION,
        "mode": mode,
        "status": status,
        "species_prefix": species_prefix,
        "source_fingerprint": source_fingerprint,
        "cds_fingerprint": cds_fingerprint,
        "output_fingerprint": file_fingerprint(output_path),
        "line_count": line_count,
        "cds_gene_ids": plan["cds_gene_ids"],
        "gene_features": plan["gene_features"],
        "renamed_gene_ids": len(plan["id_mapping"]),
        "changed_lines": counters["changed_lines"],
        "changed_values": counters["changed_values"],
        "changed_references": counters["changed_references"],
        "normalized_bare_attribute_lines": counters["normalized_bare_attribute_lines"],
        "ambiguous_count": len(plan["ambiguous"]),
        "collision_count": len(plan["collisions"]),
        "missing_gene_id_count": len(plan["missing_gene_id_lines"]),
        "repairs": plan["repairs"],
        "ambiguous": plan["ambiguous"],
        "collisions": plan["collisions"],
        "missing_gene_id_lines": plan["missing_gene_id_lines"],
    }
    audit.update(encoding_audit)
    audit_path = gff_repair_audit_path(output_path)
    write_json_atomic(audit_path, audit)
    return audit


def repair_result_fields(audit, output_path):
    if not isinstance(audit, dict):
        return {
            "repair_mode": "",
            "repair_status": "unavailable",
            "repair_audit_path": None,
            "repair_gene_ids": 0,
            "repair_references": 0,
            "repair_ambiguous": 0,
            "repair_collisions": 0,
            "normalized_bare_attribute_lines": 0,
            "invalid_utf8_bytes": 0,
            "invalid_utf8_sequences": 0,
            "invalid_utf8_line_count": 0,
            "invalid_utf8_lines": [],
        }
    return {
        "repair_mode": str(audit.get("mode") or ""),
        "repair_status": str(audit.get("status") or "unavailable"),
        "repair_audit_path": gff_repair_audit_path(output_path),
        "repair_gene_ids": int(audit.get("renamed_gene_ids", 0) or 0),
        "repair_references": int(audit.get("changed_references", 0) or 0),
        "repair_ambiguous": int(audit.get("ambiguous_count", 0) or 0),
        "repair_collisions": int(audit.get("collision_count", 0) or 0),
        "normalized_bare_attribute_lines": int(audit.get("normalized_bare_attribute_lines", 0) or 0),
        "invalid_utf8_bytes": int(audit.get("invalid_utf8_bytes", 0) or 0),
        "invalid_utf8_sequences": int(audit.get("invalid_utf8_sequences", 0) or 0),
        "invalid_utf8_line_count": int(audit.get("invalid_utf8_line_count", 0) or 0),
        "invalid_utf8_lines": list(audit.get("invalid_utf8_lines", ()) or ()),
    }
