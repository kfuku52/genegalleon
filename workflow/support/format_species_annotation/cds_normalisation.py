"""Normalise evidenced CDS/GFF inconsistencies before representative selection."""
import contextlib
import hashlib
import json
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import quote

from Bio.Data import CodonTable
from cds_model_normalisation import CdsModelNormaliser

from .common import first_token, parse_gff_attributes
from .genome_intervals import genome_intervals
from .grouping import extract_cds_header_alias_tiers, resolve_cds_header_gff_gene
from .organelle import iter_non_organelle_gff_lines
from .reference import gff_reference_mapping
from .source_identity import task_annotation_path

VERSION = 2


def audit_path(cds_path):
    return Path(str(cds_path) + ".cds-normalisation.json")


def signature(path):
    path = Path(path).resolve()
    before = path.stat()
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise OSError("Source changed during CDS normalisation: " + str(path))
    return {"path": str(path), "size": after.st_size, "sha256": digest}


def contract(task):
    return {"version": VERSION, "genetic_code": int(task.get("genetic_code", 1)),
            "gff_repair_mode": str(task.get("gff_repair_mode", "safe")),
            "implementation": signature(__file__)["sha256"],
            "evidence_implementation": signature(__import__(CdsModelNormaliser.__module__, fromlist=[""]).__file__)["sha256"],
            "inputs": {key: signature(task[key]) for key in ("cds_path", "gff_path", "genome_path", "gbff_path") if task.get(key)}}


def applicable(task):
    return task.get("gff_path") is not None and task.get("genome_path") is not None


def current_audit(task, cds_path):
    if not applicable(task):
        return None
    try:
        payload = json.loads(audit_path(cds_path).read_text())
        if isinstance(payload, dict) and payload.get("contract") == contract(task) and payload.get("output") == signature(cds_path):
            return payload
    except (FileNotFoundError, OSError, ValueError, TypeError):
        pass
    return None


def paired_audit(task, cds_path):
    """Do not publish the source GFF against an unverified corrected CDS."""
    payload = current_audit(task, cds_path)
    if not applicable(task):
        return payload
    required = task.get("_cds_normalisation_required", False)
    try:
        grouping = json.loads(Path(str(cds_path) + ".gff-grouping.json").read_text())
        required = required or grouping.get("cds_normalisation_required", False)
    except (FileNotFoundError, OSError, ValueError, TypeError, AttributeError):
        pass
    if payload is None and (required or audit_path(cds_path).exists()):
        raise ValueError("Missing or stale CDS normalisation audit; reformat CDS before paired GFF: " + str(cds_path))
    return payload


def cds_changes(normaliser, key, offset, trailing, reason):
    blocks = normaliser.models[key]
    ranges = [[row["start"], row["end"]] for row in blocks]
    for index, row in enumerate(blocks):
        trim = min(offset, ranges[index][1] - ranges[index][0])
        ranges[index][1 if row["strand"] == "-" else 0] += -trim if row["strand"] == "-" else trim
        offset -= trim
    for index in reversed(range(len(blocks))):
        row = blocks[index]
        trim = min(trailing, ranges[index][1] - ranges[index][0])
        ranges[index][0 if row["strand"] == "-" else 1] += trim if row["strand"] == "-" else -trim
        trailing -= trim
    if offset or trailing:
        raise ValueError("Partial CDS trimming exceeds annotated spans: " + key)
    updates, cumulative = {}, 0
    for row, (start, end) in zip(blocks, ranges, strict=True):
        phase = (3 - cumulative % 3) % 3
        update = {"start": start, "end": end, "phase": str(phase), "reason": reason,
                  "transcript": key, "original_start": row["start"], "original_end": row["end"], "original_phase": row["phase"]}
        updates[row["source_line"]] = update
        cumulative += end - start
    return updates


def iter_normalised_cds_records(task, state=None):
    """Yield (header, effective CDS, decision), retaining unsupported originals."""
    from .tasks import iter_task_cds_records

    state = state if state is not None else {}
    if not applicable(task) or state.get("dry_run"):
        for header, sequence in iter_task_cds_records(task):
            yield header, sequence, {"status": "not_evaluated", "reason": "missing_pair_or_dry_run",
                                     "source_length": len(sequence)}
        return
    code = int(task.get("genetic_code", 1))
    if code not in CodonTable.unambiguous_dna_by_id:
        raise ValueError("Unknown CDS genetic code: " + str(code))
    frozen = contract(task)
    reference_cache = []
    def reference_mapping():
        if not reference_cache:
            reference_cache.append(gff_reference_mapping(task["gff_path"], task["genome_path"],
                                                        reference_index=readers[0].index))
        return reference_cache[0]
    genome_scope = contextlib.ExitStack()
    readers = []
    def genome_context():
        if not readers:
            regions = ((row["seqid"], row["start"], row["end"]) for row in normaliser.features
                       if row["feature"] in {"CDS", "exon"})
            readers.append(genome_scope.enter_context(genome_intervals(
                task["genome_path"], regions, scratch_dir=task.get("_normalisation_scratch", tempfile.gettempdir()))))
        return contextlib.nullcontext(readers[0])
    def attribute_parser(text):
        return {key: ",".join(values) for key, values in parse_gff_attributes(text).items()}
    def gff_lines():
        for line in iter_non_organelle_gff_lines(task_annotation_path(task)):
            if line.startswith("#") or len(line.rstrip("\r\n").split("\t")) == 9:
                yield line
    normaliser = CdsModelNormaliser(
        {"species": task["species_prefix"], "gff": str(task["gff_path"]), "genome": str(task["genome_path"]), "genetic_code": code},
        Path(task.get("_normalisation_scratch", tempfile.gettempdir())), "format", genome_context=genome_context,
        reference_mapping=reference_mapping, attribute_parser=attribute_parser, gff_lines=gff_lines, coding_only=True)
    index = task.get("_gff_cds_grouping_index")
    if index is None:
        from .grouping import build_gff_cds_grouping_index
        index = build_gff_cds_grouping_index(task)
    by_gene = defaultdict(set)
    for transcript, gene in index["transcript_gene_tokens"].items():
        if transcript in normaliser.models:
            by_gene[gene].add(transcript)
    updates = {}
    try:
        for number, (header, sequence) in enumerate(iter_task_cds_records({**task, "_genome_interval_context": genome_context}), 1):
            raw = "".join(sequence.split()).upper()
            match = resolve_cds_header_gff_gene(task, header, grouping_index=index)
            keys = by_gene.get(match["gene_token"], set()) if match["status"] == "mapped" else set()
            primary, _gene_aliases = extract_cds_header_alias_tiers(task, header)
            bound = {key for key in keys if key in primary or any(
                value in primary for row in normaliser.models[key] for name in ("protein_id", "transcript_id")
                for value in row["attributes"].get(name, "").split(","))}
            if bound:
                keys = bound
            identifier = str(number)
            normaliser.normalise(identifier, raw, keys)
            decision = normaliser.rows[-1]
            decision.update(raw_cds_id=first_token(header), gene_token=match["gene_token"], source_length=len(raw))
            correction = normaliser.corrected.get(identifier)
            if correction is not None:
                raw = correction["sequence"]
                for evidence in correction["evidence"]:
                    if evidence["reason"] == "annotated_partial_frame":
                        changes = cds_changes(normaliser, evidence["transcript"], evidence["offset"],
                                              evidence["trailing_bases_omitted"], evidence["reason"])
                        for source_line, update in changes.items():
                            if source_line in updates and updates[source_line] != update:
                                raise ValueError("Conflicting CDS/GFF normalisations: " + evidence["transcript"])
                            updates[source_line] = update
            else:
                if decision["status"] == "excluded":
                    decision["status"] = "retained_unresolved"
            yield header, raw, decision
        if contract(task) != frozen:
            raise OSError("Inputs changed during CDS normalisation")
        state.update(contract=frozen, counts=dict(Counter(row["status"] for row in normaliser.rows)),
                     reasons=dict(Counter(row["reason"] for row in normaliser.rows)),
                     records=[row for row in normaliser.rows if row["status"] != "unchanged"],
                     phase_convention_votes=dict(normaliser.votes), gff_updates=updates)
    finally:
        try:
            normaliser.close()
        finally:
            genome_scope.close()


def write_audit(task, cds_path, state):
    if not applicable(task) or state.get("dry_run"):
        return None
    from .gff_repair import write_json_atomic
    payload = {**state, "output": signature(cds_path), "original_sources_modified": False,
               "counting_unit": "input CDS record before longest-isoform selection"}
    write_json_atomic(audit_path(cds_path), payload)
    return payload


def normalise_gff_lines(lines, updates):
    for line in lines:
        update = updates.get(line.rstrip("\r\n"))
        if update is None:
            yield line
            continue
        if update["start"] == update["end"]:
            continue
        fields = line.rstrip("\r\n").split("\t")
        fields[3], fields[4], fields[7] = str(update["start"] + 1), str(update["end"]), update["phase"]
        fields[8] += ";gg_cds_normalisation=" + quote(update["reason"], safe="")
        yield "\t".join(fields) + "\n"
