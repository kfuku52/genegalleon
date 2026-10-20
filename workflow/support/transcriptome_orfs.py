"""Explicit transcript FASTA to CDS conversion through CDSKIT's native ORF caller."""

import hashlib
import importlib.metadata
import json
import re
from pathlib import Path


def input_type(value):
    kind = str(value.get("cds_input_type") or "cds").strip().lower()
    if kind not in {"cds", "transcriptome"}:
        raise ValueError("cds_input_type must be cds or transcriptome")
    return kind


def validate_manifest_row(row):
    if input_type(row) == "transcriptome":
        onekp_sample = (str(row.get("provider", "")).strip().lower() == "onekp"
                        and re.fullmatch(r"[A-Za-z]{4}", str(row.get("id", "")).strip()))
        # The native 1KP resolver verifies the official sample identity and
        # selects its transcript assembly before any nucleotide download.
        if not (row.get("cds_url") or row.get("local_cds_path") or onekp_sample):
            raise ValueError("Transcriptome input requires an explicit sequence URL/path")
        if any(row.get(prefix + role + suffix) for role in ("gff", "gbff", "genome")
               for prefix, suffix in (("", "_url"), ("", "_path"), ("local_", "_path"))) or row.get("paired_curation"):
            raise ValueError("Transcript ORFs cannot be paired with genomic annotation roles")


def apply_input_type(task, row):
    validate_manifest_row(row)
    if input_type(row) == "transcriptome":
        task["cds_input_type"] = "transcriptome"
    return task


def verify_input_type(task, planned):
    if input_type(task) != input_type(planned):
        raise ValueError("Staged nucleotide input type differs from the frozen plan")


def extract_orf(task, identifier, sequence):
    from cdskit.longestcds import choose_best_candidate

    if any(task.get(role + "_path") for role in ("gff", "gbff", "genome")):
        raise ValueError("Transcript ORFs cannot be paired with genomic annotation roles")
    invalid = set(sequence.upper()) - set("ACGTRYSWKMBDHVNU")
    if invalid:
        raise ValueError("Invalid transcript nucleotides for " + identifier + ": " + repr(sorted(invalid)))
    # This is the same complete-first policy as native cdskit longestorf in
    # gg_transcriptome_generation; do not substitute a project-local ORF caller.
    candidate = choose_best_candidate(sequence.upper().replace("U", "T"), int(task.get("genetic_code", 1)))
    audit = task.setdefault("_transcriptome_orf_audit", [])
    row = dict(source_id=identifier, transcript_length=len(sequence))
    if candidate is None:
        audit.append(dict(row, status="no_orf"))
        return None
    audit.append(dict(row, status="extracted", strand=candidate.strand,
                      start=candidate.start_1based, end=candidate.end_1based,
                      category=candidate.category, cds_length=len(candidate.output_seq)))
    return candidate.output_seq


def signature(path):
    path = Path(path).resolve()
    before = path.stat()
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    after = path.stat()
    fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
    if any(getattr(before, key) != getattr(after, key) for key in fields):
        raise ValueError("ORF input/output changed while hashing: " + str(path))
    return dict(path=str(path), size=after.st_size, sha256=h.hexdigest())


def contract(task):
    from cdskit import longestcds

    return dict(version=1, cds_input_type=input_type(task), selection="complete-first",
                gene_grouping="declared-trinity-gene-v1",
                genetic_code=int(task.get("genetic_code", 1)), source=signature(task["cds_path"]),
                implementation=signature(__file__)["sha256"],
                cdskit_version=importlib.metadata.version("cdskit"),
                cdskit_implementation=signature(longestcds.__file__)["sha256"])


def audit_path(output):
    return Path(str(output) + ".orf-extraction.json")


def verify_existing(task, output):
    path = audit_path(output)
    if input_type(task) == "cds":
        if path.exists():
            raise ValueError("Existing CDS was derived from transcripts; input type changed")
        return
    try:
        value = json.loads(path.read_text())
        if value["contract"] == contract(task) and value["output"] == signature(output):
            return
    except (OSError, ValueError, KeyError, TypeError):
        pass
    raise ValueError("Missing or stale transcript ORF audit; use a new workspace or explicit overwrite")


def write_audit(task, output, writer):
    if input_type(task) == "transcriptome":
        writer(audit_path(output), dict(contract=contract(task), output=signature(output),
                                      coordinate_system="source transcript; not genomic",
                                      records=task.get("_transcriptome_orf_audit", [])))
