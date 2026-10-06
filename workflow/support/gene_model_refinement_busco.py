#!/usr/bin/env python3
"""Evaluate paired representative DNA CDS sets with one frozen BUSCO contract."""

from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from busco_quality_metadata import parse_short_summary
from dated_tree_presentation import STATUS, STATUS_COLOURS, STATUS_LABELS
from gene_model_refinement import verify_inputs
from input_generation_array_state import FreshDigestBatch, atomic_json, digest
from species_labeling import extract_species_label

RESCUE_SUPPORT = ("nearest_only", "balanced_only", "both")
RESCUE_SUPPORT_COLOURS = ("#5275b5", "#9a6ab2", "#37956f")
RESCUE_SUPPORT_LABELS = ("Nearest relatives only", "Phylogenetically balanced only", "Both reference groups")
RESCUE_SELF_COLOUR = "#b59b5b"
RESCUE_SELF_LABEL = "Self-species only"
PATH_SUPPORT = (*RESCUE_SUPPORT, "self_only", "other_interspecies")
PATH_SUPPORT_COLOURS = (*RESCUE_SUPPORT_COLOURS, RESCUE_SELF_COLOUR, "#8a929b")
PATH_SUPPORT_LABELS = (*RESCUE_SUPPORT_LABELS, RESCUE_SELF_LABEL, "Other interspecies only")
# Keep the legacy count fields above readable; the current view uses all three
# support types, including self support in a mixed donor set.
SUPPORT_GROUPS = ("self_only", "relative_only", "phylogenetic_only", "multiple")
SUPPORT_GROUP_COLOURS = (RESCUE_SELF_COLOUR, *RESCUE_SUPPORT_COLOURS)
SUPPORT_GROUP_LABELS = ("S only (self)", "R only (relatives)", "P only (phylogenetic)", "Multiple (at least two of S/R/P)")
PATH_SUPPORT_GROUPS = (*SUPPORT_GROUPS, "other_interspecies")
PATH_SUPPORT_GROUP_COLOURS = (*SUPPORT_GROUP_COLOURS, "#8a929b")
PATH_SUPPORT_GROUP_LABELS = (*SUPPORT_GROUP_LABELS, "Other interspecies only")
REPEAT_GROUPS = ("te_hit", "other_repeat_hit", "no_repeat_hit", "not_assessed")
REPEAT_GROUP_COLOURS = ("#b45158", "#d99843", "#477b80", "#d6dbe1")
REPEAT_GROUP_LABELS = ("TE overlap", "Other / unclassified repeat overlap", "Assessed: no repeat overlap", "Not assessed")


def repeat_group(record):
    """Classify CDS overlap without inferring whether a locus is a functional gene."""
    if record.get("status") == "not_provided":
        if any(record.get(key) is not None for key in ("masked_fraction", "te_annotated_fraction", "classes")):
            raise ValueError("Unavailable repeat evidence must have null measurements")
        return "not_assessed"
    masked, te = record.get("masked_fraction"), record.get("te_annotated_fraction")
    if (record.get("status") != "available" or type(masked) not in (int, float) or type(te) not in (int, float)
            or not 0 <= te <= masked <= 1 or not isinstance(record.get("classes"), list)
            or any(not isinstance(c, str) or not c for c in record["classes"])
            or bool(masked) != bool(record["classes"])):
        raise ValueError("Invalid repeat overlap measurements")
    return "te_hit" if te > 0 else "other_repeat_hit" if masked > 0 else "no_repeat_hit"


def collect_rescue_repeat_evidence(changes, evidence_dir=None):
    """Join immutable evidence audits to the exact rescued loci in Before."""
    from rescue_model_evidence import read_json_snapshot

    if evidence_dir is not None and not Path(evidence_dir).is_dir():
        raise ValueError("Rescue evidence directory does not exist")
    updates = {}
    snapshots = {}
    for name, value in changes["species"].items():
        if value["refinement_status"] == "not_analysed":
            updates[name] = (None, None)
            continue
        total = value["prior_rescued_loci"]
        counts = dict.fromkeys(REPEAT_GROUPS, 0)
        counts["not_assessed"] = total
        directory = Path(evidence_dir) / name if evidence_dir is not None else None
        if directory is None or not directory.exists():
            updates[name] = (counts, None)
            continue
        source = changes.get("evidence", {}).get(name, {})
        loci = source.get("rescued_loci_support")
        if not isinstance(loci, dict) or len(loci) != total:
            raise ValueError("Repeat audit needs complete rescued-locus identities: " + name)
        receipt, receipt_hash = read_json_snapshot(directory / "receipt.json")
        records, records_hash = read_json_snapshot(directory / "evidence.json")
        key = receipt.get("key", {})
        inputs = key.get("inputs", {})
        selection = changes.get("rescue_reference_selection", {})
        expected = {f"/rescued/{name}/models.json": source.get("rescue_models_sha256"),
                    f"/rescued/{name}/receipt.json": source.get("rescue_receipt_sha256")}
        if (key.get("schema") != 1 or key.get("species") != name or not isinstance(inputs, dict)
                or receipt.get("files", {}).get("evidence.json") != records_hash
                or any(wanted is None or [h for p, h in inputs.items() if p.endswith(suffix)] != [wanted]
                       for suffix, wanted in expected.items())
                or selection.get("plan_sha256") not in [h for p, h in inputs.items() if p.endswith("/plan.json")]):
            raise ValueError("Repeat audit belongs to different rescue models: " + name)
        snapshots[str(directory / "receipt.json")] = receipt_hash
        snapshots[str(directory / "evidence.json")] = records_hash
        by_id = {}
        if not isinstance(records, list):
            raise ValueError("Repeat audit records must be a list")
        for record in records:
            identifier = record["model_id"]
            if identifier in by_id or record.get("rescue_status") != "accepted":
                raise ValueError("Duplicate or unaccepted repeat audit model: " + identifier)
            by_id[identifier] = record
        if set(by_id) != {r["source_model_id"] for r in loci.values()}:
            raise ValueError("Repeat audit model membership differs from rescued loci: " + name)
        counts = dict.fromkeys(REPEAT_GROUPS, 0)
        per_locus = {}
        for gene, model in loci.items():
            record = by_id[model["source_model_id"]]["repeat"]
            category = repeat_group(record)
            counts[category] += 1
            per_locus[gene] = {"source_model_id": model["source_model_id"], "category": category, **record}
        updates[name] = (counts, {"audit_receipt_sha256": receipt_hash, "audit_records_sha256": records_hash,
                                 "audit_directory": str(directory.resolve()), "loci": per_locus})
    boundary = FreshDigestBatch()
    if boundary.read(snapshots) != snapshots:
        raise ValueError("Repeat audit changed while loading")
    boundary.check()
    for name, (counts, evidence) in updates.items():
        changes["species"][name]["rescue_repeat_groups"] = counts
        if evidence is not None:
            changes.setdefault("evidence", {}).setdefault(name, {})["rescued_loci_repeat"] = evidence
        else:
            changes.get("evidence", {}).get(name, {}).pop("rescued_loci_repeat", None)
    changes["repeat_group_classification"] = (
        "CDS overlap with the supplied repeat annotation: any TE-labelled overlap takes priority; otherwise any "
        "other/unclassified repeat overlap; assessed no overlap; or not assessed. Each rescued locus counts once. "
        "Overlap does not establish TE origin or lack of gene function; no overlap does not establish a true gene."
    )
    return changes


def group_support_evidence(evidence, species, selection, *, allow_other=False):
    """Count each locus/path once by S/R/P membership, not donor-species count."""
    nearest = set(selection["nearest_references"][species]) - {species}
    balanced = set(selection["common_references"]) - {species}
    groups = PATH_SUPPORT_GROUPS if allow_other else SUPPORT_GROUPS
    counts = dict.fromkeys(groups, 0)
    for identifier, record in evidence.items():
        donors = record.get("supporting_donors")
        if (not isinstance(donors, list) or not donors
                or any(not isinstance(d, str) or not d for d in donors)):
            raise ValueError("Missing supporting donor evidence: " + identifier)
        donors = set(donors)
        if not allow_other and donors - nearest - balanced - {species}:
            raise ValueError("Rescue donor differs from the frozen reference selection: " + identifier)
        flags = (species in donors, bool(donors & nearest), bool(donors & balanced))
        number = sum(flags)
        category = "multiple" if number >= 2 else SUPPORT_GROUPS[flags.index(True)] if number else "other_interspecies"
        if category not in counts:
            raise ValueError("Rescue lacks S/R/P supporting donors: " + identifier)
        counts[category] += 1
    return counts


def regroup_model_support(changes):
    """Add the S/R/P view from per-model evidence without changing legacy fields."""
    selection = changes.get("rescue_reference_selection")
    if selection is None:
        return changes
    updates = {}
    for species, value in changes["species"].items():
        updates[species] = {}
        if value["refinement_status"] == "not_analysed":
            if any(value.get(k) is not None for k in ("rescue_support_groups", "accepted_path_support_groups")):
                raise ValueError("Unanalysed support groups must be unavailable")
            continue
        for legacy, output, evidence_key, total in (
            ("rescue_support_counts", "rescue_support_groups", "rescued_loci_support", value["prior_rescued_loci"]),
            ("accepted_path_support_counts", "accepted_path_support_groups", "accepted_paths_support",
             value["accepted_repair_paths"] + value["accepted_isoform_paths"]),
        ):
            if value.get(legacy) is None:
                if value.get(output) is not None:
                    raise ValueError("Support groups lack their source counts: " + species)
                continue
            evidence = changes.get("evidence", {}).get(species, {}).get(evidence_key)
            if not isinstance(evidence, dict) or len(evidence) != total:
                raise ValueError("Support regrouping needs complete per-model evidence: " + species)
            counts = group_support_evidence(evidence, species, selection, allow_other=legacy == "accepted_path_support_counts")
            if value.get(output) is not None and value[output] != counts:
                raise ValueError("Support groups disagree with per-model evidence: " + species)
            updates[species][output] = counts
    for species, value in updates.items():
        changes["species"][species].update(value)
    changes["support_group_classification"] = (
        "S = self-species homology; R = nearest relatives; P = phylogenetically balanced references. "
        "Only means exactly one of S/R/P; multiple means at least two of these support types, not two donor species. "
        "A donor in both frozen R/P lists supplies both types. Other interspecies only means no S/R/P support; "
        "additional unselected donors remain in per-model evidence. Target RNA is separate from S."
    )
    return changes


def classify_accepted_path_support(models, species, selection, allowed_species):
    """Partition accepted paths by recorded homology donors, keeping RNA separate."""
    nearest = set(selection["nearest_references"][species]) - {species}
    balanced = set(selection["common_references"]) - {species}
    counts, evidence = dict.fromkeys(PATH_SUPPORT, 0), {}
    for model in models:
        if model.get("status") != "accepted":
            continue
        candidate = model["candidate"]
        identifier = candidate["candidate_id"]
        donors = model.get("donors")
        if (not isinstance(donors, list) or not donors
                or any(not isinstance(d, str) or d not in allowed_species for d in donors)):
            raise ValueError("Accepted coding path lacks valid supporting donors: " + identifier)
        donors = set(donors)
        # Target species is bound by the prediction receipt, not a candidate field.
        if (candidate["source_transcript_id"] != identifier or identifier in evidence
                or model["change_type"] not in {"model_revision", "isoform_addition"}
                or donors != set(candidate["support"]["donors"])
                or donors != {a["donor_species"] for a in model["alignments"]}):
            raise ValueError("Accepted coding-path identity/support records disagree: " + identifier)
        near, common = bool(donors & nearest), bool(donors & balanced)
        category = ("both" if near and common else "nearest_only" if near else "balanced_only" if common
                    else "self_only" if donors == {species} else "other_interspecies")
        counts[category] += 1
        evidence[identifier] = {"gene_id": model["gene_id"], "change_type": model["change_type"],
                                "supporting_donors": sorted(donors), "category": category,
                                "other_supporting_donors": sorted(donors - nearest - balanced - {species})}
    return counts, evidence


def collect_accepted_path_support(root, changes):
    """Extend a saved model summary from verified predictions without rereading rescue models."""
    from plot_gene_model_refinement import verified_json

    plan_hash = digest(root / "plan.json")
    request = json.loads((root / "plan.json").read_text())["request"]
    if (changes["plan_sha256"] != plan_hash
            or changes["effective_receipt_sha256"] != digest(root / "effective/receipt.json")
            or {s for s, c in changes["species"].items() if c["refinement_status"] == "analysed"}
            != set(request["sources"])):
        raise ValueError("Saved model counts belong to another refinement publication")
    selection = changes["rescue_reference_selection"]
    updates = {}
    for species, value in changes["species"].items():
        if value["refinement_status"] == "not_analysed":
            continue
        directory = root / "predictions" / species
        receipt_hash = digest(directory / "receipt.json")
        if receipt_hash != changes["evidence"][species]["prediction_receipt_sha256"]:
            raise ValueError("Saved prediction receipt changed: " + species)
        models = verified_json(directory, "predictions.json", plan_hash)
        counts, evidence = classify_accepted_path_support(models, species, selection, set(request["sources"]))
        for key, kind in (("accepted_repair_paths", "model_revision"), ("accepted_isoform_paths", "isoform_addition")):
            if value[key] != sum(e["change_type"] == kind for e in evidence.values()):
                raise ValueError("Saved accepted coding-path counts changed: " + species)
        if digest(directory / "receipt.json") != receipt_hash:
            raise ValueError("Prediction receipt changed while reading support: " + species)
        updates[species] = (counts, evidence)
    if digest(root / "plan.json") != plan_hash:
        raise ValueError("Refinement plan changed while collecting path support")
    for species, (counts, evidence) in updates.items():
        changes["species"][species]["accepted_path_support_counts"] = counts
        changes["evidence"][species]["accepted_paths_support"] = evidence
    return regroup_model_support(changes)


def classify_rescue_support(models, species, rescued, plan):
    """Count each accepted gene once using all consolidated supporting donors."""
    nearest = set(plan["nearest_references"][species]) - {species}
    balanced = set(plan["common_references"]) - {species}
    allowed = set(plan["donors"][species])
    if allowed != nearest | balanced:
        raise ValueError("Rescue donor groups differ from the frozen donor selection")
    counts = dict.fromkeys(RESCUE_SUPPORT, 0)
    evidence = {}
    aliases = {}
    for gene in rescued:
        for identifier in rescued[gene] if isinstance(rescued, dict) else [gene]:
            if identifier in aliases and aliases[identifier] != gene:
                raise ValueError("Rescue model maps to multiple source genes: " + identifier)
            aliases[identifier] = gene
    self_comparisons = {job["id"] for job in plan.get("synteny_jobs", [])
                        if job["a"] == job["b"] == species and job.get("kind") == "self"}
    for model in models:
        if model.get("status") != "accepted":
            continue
        identifier = model["model_id"]
        gene = aliases.get(identifier)
        if gene is None:
            raise ValueError("Accepted rescue gene IDs differ from the source annotation: " + species)
        if gene in evidence:
            raise ValueError("Duplicate accepted rescue gene: " + identifier)
        support = model.get("support")
        if not isinstance(support, list) or not support:
            raise ValueError("Accepted rescue gene lacks supporting donors: " + identifier)
        donors = set()
        for entry in support:
            donor = entry.get("donor")
            if (entry.get("target") != species
                    or (donor not in allowed and not (donor == species and entry.get("comparison") in self_comparisons))):
                raise ValueError("Rescue support donor/target differs from the frozen plan")
            donors.add(donor)
        near, common = bool(donors & nearest), bool(donors & balanced)
        category = "both" if near and common else "nearest_only" if near else "balanced_only" if common else "self_only"
        if category != "self_only":
            counts[category] += 1
        evidence[gene] = {"source_model_id": identifier, "supporting_donors": sorted(donors), "category": category}
    if set(evidence) != set(rescued):
        raise ValueError("Accepted rescue gene IDs differ from the source annotation: " + species)
    return counts, evidence


def input_pairs(root, cds_dir=None):
    """Native sources are frozen; extra CDS-only species are explicit passthroughs."""
    plan = json.loads((root / "plan.json").read_text())
    effective = {r["species"]: r for r in verify_inputs(root / "effective" / "inputs.tsv")}
    if set(effective) != set(plan["request"]["sources"]):
        raise ValueError("Refinement species membership differs from the frozen plan")
    pairs = {}
    for name, source in plan["request"]["sources"].items():
        before = Path(source["fasta"])
        if digest(before) != plan["request"]["files"][str(before)]:
            raise ValueError("Refinement source CDS changed: " + str(before))
        pairs[name] = {"species": name, "before": str(before), "after": effective[name]["cds"],
                       "refinement_status": "analysed", "reason": ""}
    seen = set()
    if cds_dir:
        for path in sorted(cds_dir.iterdir()):
            if path.name.startswith(".") or not path.is_file() or not re.search(r"\.(fa|fasta|fna)(\.gz)?$", path.name):
                continue
            name = extract_species_label(path.name, strip_extension=True)
            if not name:
                raise ValueError("Cannot identify CDS species: " + path.name)
            if name in seen:
                raise ValueError("Ambiguous CDS files for species: " + name)
            seen.add(name)
            if name not in pairs:
                pairs[name] = {"species": name, "before": str(path.resolve()), "after": str(path.resolve()),
                               "refinement_status": "not_analysed",
                               "reason": "No matching genome/GFF in the refinement plan; CDS retained unchanged"}
    return [pairs[n] for n in sorted(pairs)]


def read_result(path):
    text = path.read_text()
    row = parse_short_summary(text, path)
    fields = {
        "complete": r"Complete BUSCOs \(C\)", "single": r"Complete and single-copy BUSCOs \(S\)",
        "duplicated": r"Complete and duplicated BUSCOs \(D\)", "fragmented": r"Fragmented BUSCOs \(F\)",
        "missing": r"Missing BUSCOs \(M\)", "total": r"Total BUSCO groups searched",
    }
    for name, label in fields.items():
        hits = re.findall(r"^\s*(\d+)\s+" + label, text, re.M)
        if len(hits) != 1:
            raise ValueError("Missing/ambiguous BUSCO count: " + name)
        row[name] = int(hits[0])
    if (row["total"] <= 0 or row["complete"] != row["single"] + row["duplicated"]
            or row["complete"] + row["fragmented"] + row["missing"] != row["total"]):
        raise ValueError("Inconsistent BUSCO counts")
    match = re.search(r"Creation date:\s*([^,)]+)", text)
    row["lineage_creation_date"] = match.group(1).strip() if match else ""
    row["dependencies"] = dict(re.findall(r"^\s*(hmmsearch|metaeuk):\s*(\S+)", text, re.M))
    for key in ("lineage", "mode", "busco_version", "lineage_creation_date"):
        if not row[key]:
            raise ValueError("Missing BUSCO comparability metadata: " + key)
    return row


def paired_result(pair, before, after):
    keys = ("lineage", "lineage_creation_date", "total", "mode", "busco_version", "dependencies")
    if any(before[k] != after[k] for k in keys):
        raise ValueError("Noncomparable BUSCO pair: " + pair["species"])
    return dict(pair, before_result=before, after_result=after,
                delta_complete=after["complete"] - before["complete"],
                delta_complete_pp=100 * (after["complete"] - before["complete"]) / before["total"],
                delta_duplicated=after["duplicated"] - before["duplicated"])


def run_one(pair, phase, report, contract, cpus):
    source = Path(pair[phase])
    source_hash = digest(source)
    key = {"contract": contract, "source_sha256": source_hash}
    directory = report / "runs" / pair["species"] / phase
    directory.mkdir(parents=True, exist_ok=True)
    summary, receipt = directory / "summary.txt", directory / "receipt.json"
    full_table = directory / "full_table.tsv"
    with (directory / ".lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if receipt.exists():
            saved = json.loads(receipt.read_text())
            if (saved["key"] != key or digest(summary) != saved["summary_sha256"]
                    or not full_table.is_file() or digest(full_table) != saved.get("full_table_sha256")):
                raise ValueError("BUSCO cache changed; use a new report directory")
            return read_result(summary)
        # Identical bytes and an identical complete contract permit reuse, including excluded species.
        previous = report / "runs" / pair["species"] / "before"
        if phase == "after" and (previous / "receipt.json").exists():
            saved = json.loads((previous / "receipt.json").read_text())
            if (saved["key"] == key and digest(previous / "summary.txt") == saved["summary_sha256"]
                    and (previous / "full_table.tsv").is_file()
                    and digest(previous / "full_table.tsv") == saved.get("full_table_sha256")):
                shutil.copyfile(previous / "summary.txt", summary)
                shutil.copyfile(previous / "full_table.tsv", full_table)
                atomic_json(receipt, {"key": key, "summary_sha256": digest(summary),
                                      "full_table_sha256": digest(full_table), "reused_identical_before": True})
                return read_result(summary)
        # Preserve small QC tables and logs; discard large temporary predictor files.
        with tempfile.TemporaryDirectory(prefix="gg-refinement-busco-") as scratch:
            scratch = Path(scratch)
            opener = gzip.open if source.suffix == ".gz" else open
            with opener(source, "rb") as src, (scratch / "input.fasta").open("wb") as dst:
                shutil.copyfileobj(src, dst)
            support = Path(__file__).resolve().parent
            command = ["bash", "-c", 'source "$1"; gg_run_busco_with_metaeuk_modified_fas_compat "${@:2}"',
                       "gg-busco", str(support / "gg_busco.sh"), "--in", str(scratch / "input.fasta"),
                       "--mode", "transcriptome", "--out", "busco", "--out_path", str(scratch),
                       "--cpu", str(cpus), "--evalue", "1e-03", "--limit", "20", "--offline",
                       "--lineage_dataset", contract["lineage_path"], "--download_path", contract["download_path"]]
            with (directory / "busco.log").open("w") as log:
                subprocess.run(command, cwd=scratch, stdout=log, stderr=subprocess.STDOUT, check=True)
            summaries = sorted((scratch / "busco").glob("short_summary.specific.*.txt"))
            if len(summaries) != 1:
                raise ValueError("Expected one BUSCO specific summary")
            shutil.copyfile(summaries[0], summary)
            tables = list((scratch / "busco").glob("run_*/full_table.tsv")) + list((scratch / "busco").glob("run_*/full_table.tsv.gz"))
            if len(tables) != 1:
                raise ValueError("Expected one BUSCO full table")
            raw_table = tables[0]
            if raw_table.suffix != ".gz":
                shutil.copyfile(raw_table, full_table)
            else:
                with gzip.open(raw_table, "rb") as src, full_table.open("wb") as dst:
                    shutil.copyfileobj(src, dst)
        if digest(source) != source_hash:
            raise OSError("BUSCO input mutated during evaluation")
        result = read_result(summary)
        atomic_json(receipt, {"key": key, "summary_sha256": digest(summary),
                              "full_table_sha256": digest(full_table), "reused_identical_before": False})
        return result


def collect_model_changes(root, pairs, rescue_output=None):
    """Bind gene rescue and accepted coding-path counts to the same publication."""
    from format_species_annotation.common import parse_gff_attributes
    from plot_gene_model_refinement import verified_json

    plan_hash = digest(root / "plan.json")
    request = json.loads((root / "plan.json").read_text())["request"]
    analysed = {r["species"] for r in pairs if r["refinement_status"] == "analysed"}
    if analysed != set(request["sources"]):
        raise ValueError("Model-count species differ from the refinement plan")
    stats, evidence = {}, {}
    batch = FreshDigestBatch()
    anchor = request.get("rescue_output")
    if rescue_output is not None and anchor and Path(rescue_output).resolve() != Path(anchor).resolve():
        raise ValueError("Rescue support directory differs from the refinement plan")
    rescue_root = Path(rescue_output or anchor) if rescue_output or anchor else None
    # Refinement can reuse synteny anchors before any missing-gene augmentation.
    if rescue_output is None and rescue_root is not None and not (rescue_root / "augmented/receipt.json").is_file():
        rescue_root = None
    rescue_plan = augmented = None
    if rescue_root is not None:
        from rescue_model_evidence import read_json_snapshot
        rescue_plan, rescue_hash = read_json_snapshot(rescue_root / "plan.json")
        augmented, augmented_hash = read_json_snapshot(rescue_root / "augmented/receipt.json")
        if augmented["key"]["plan"] != rescue_hash:
            raise ValueError("Rescue augmentation belongs to a different plan")
        expected = request["files"].get(str(rescue_root / "plan.json"))
        if expected is not None and expected != rescue_hash:
            raise ValueError("Frozen rescue plan changed")
        hashes = batch.read([rescue_root / "plan.json", rescue_root / "augmented/receipt.json"])
        if (hashes[str(rescue_root / "plan.json")] != rescue_hash
                or hashes[str(rescue_root / "augmented/receipt.json")] != augmented_hash):
            raise ValueError("Rescue metadata changed while loading support")
    for pair in pairs:
        name = pair["species"]
        if pair["refinement_status"] == "not_analysed":
            stats[name] = {"refinement_status": "not_analysed", "prior_rescued_loci": None,
                           "accepted_repair_paths": None, "accepted_isoform_paths": None}
            continue
        source = request["sources"][name]
        if (Path(pair["before"]).resolve() != Path(source["fasta"]).resolve()
                or Path(pair["after"]).resolve() != (root / "effective/species_cds" / (name + ".fa")).resolve()):
            raise ValueError("Model counts belong to different BUSCO inputs: " + name)
        metadata = verified_json(root / "catalog" / name, "catalog_metadata.json", plan_hash)
        models = verified_json(root / "predictions" / name, "predictions.json", plan_hash)
        gff = Path(source["gff"])
        source_hash = batch.read([gff])[str(gff)]
        if source_hash != metadata["sources"]["gff"]["sha256"] or source_hash != request["files"][str(gff)]:
            raise ValueError("Rescue source annotation changed: " + name)
        rescued, transcripts = {}, []
        opener = gzip.open if gff.suffix == ".gz" else open
        with opener(gff, "rt") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) == 9 and fields[1:3] == ["genegalleon_rescue", "gene"]:
                    identifiers = parse_gff_attributes(fields[8]).get("ID", [])
                    if len(identifiers) != 1:
                        raise ValueError("Rescued gene lacks a unique ID: " + name)
                    rescued.setdefault(identifiers[0], set()).add(identifiers[0])
                elif len(fields) == 9 and fields[1] == "genegalleon_rescue" and fields[2] in {"mRNA", "transcript"}:
                    attr = parse_gff_attributes(fields[8])
                    if len(attr.get("ID", [])) != 1 or not attr.get("Parent"):
                        raise ValueError("Rescued transcript lacks an ID or gene parent: " + name)
                    transcripts.append((attr["ID"][0], attr["Parent"]))
        for identifier, parents in transcripts:
            for parent in parents:
                if parent in rescued:
                    rescued[parent].add(identifier)
        accepted = [m for m in models if m["status"] == "accepted"]
        stats[name] = {
            "refinement_status": "analysed", "prior_rescued_loci": len(rescued),
            "accepted_repair_paths": sum(m["change_type"] == "model_revision" for m in accepted),
            "accepted_isoform_paths": sum(m["change_type"] == "isoform_addition" for m in accepted),
        }
        evidence[name] = {"source_gff_sha256": source_hash,
                          "catalog_receipt_sha256": digest(root / "catalog" / name / "receipt.json"),
                          "prediction_receipt_sha256": digest(root / "predictions" / name / "receipt.json")}
        if rescue_root is not None:
            path_counts, path_support = classify_accepted_path_support(models, name, rescue_plan, set(request["sources"]))
            stats[name]["accepted_path_support_counts"] = path_counts
            evidence[name]["accepted_paths_support"] = path_support
            from rescue_model_evidence import json_array, read_json_snapshot
            worker = rescue_root / "rescued" / name
            receipt, receipt_hash = read_json_snapshot(worker / "receipt.json")
            expected = request["files"].get(str(worker / "receipt.json"))
            if (receipt["key"]["plan"] != rescue_hash or receipt["key"]["species"] != name
                    or augmented["key"]["rescue_receipts"].get(name) != receipt_hash
                    or (expected is not None and expected != receipt_hash)
                    or augmented["files"].get("species_gff/" + name + ".rescue.gff3") != source_hash):
                raise ValueError("Rescue support publication differs from the source annotation: " + name)
            hashes = batch.read([worker / "receipt.json", worker / "models.json"])
            if (hashes[str(worker / "receipt.json")] != receipt_hash
                    or hashes[str(worker / "models.json")] != receipt["files"]["models.json"]):
                raise ValueError("Frozen rescue models changed: " + name)
            counts, support = classify_rescue_support(json_array(worker / "models.json"), name, rescued, rescue_plan)
            stats[name]["rescue_support_counts"] = counts
            stats[name]["rescue_self_only_loci"] = sum(s["category"] == "self_only" for s in support.values())
            evidence[name].update(rescue_receipt_sha256=receipt_hash,
                                  rescue_models_sha256=receipt["files"]["models.json"], rescued_loci_support=support)
    batch.check()
    if digest(root / "plan.json") != plan_hash:
        raise ValueError("Refinement plan changed while collecting model counts")
    result = {"schema": 1, "plan_sha256": plan_hash,
            "effective_receipt_sha256": digest(root / "effective/receipt.json"),
            "species": stats, "evidence": evidence,
            "units": "Prior missing-gene rescue counts unique source gene loci already present before refinement; repairs and additional isoforms count accepted coding paths, not unique loci."}
    if rescue_root is not None:
        result["rescue_reference_selection"] = {
            "rescue_output": str(rescue_root), "plan_sha256": rescue_hash, "augmented_receipt_sha256": augmented_hash,
            "nearest_references": rescue_plan["nearest_references"], "common_references": rescue_plan["common_references"],
            "classification": "Supporting donors are deduplicated across all consolidated support records. A donor in both frozen reference lists supports both groups; both does not require two distinct donor species. Self-species-only loci are recorded separately and excluded from the three interspecies support groups.",
        }
    return regroup_model_support(result)


def plot_comparison(rows, output, model_changes=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.ticker import MaxNLocator

    stacked_rescue = stacked_paths = grouped_support = False
    if model_changes is not None:
        stats = model_changes["species"]
        stacked_rescue = any(v.get("rescue_support_counts") is not None for v in stats.values())
        stacked_paths = any(v.get("accepted_path_support_counts") is not None for v in stats.values())
        grouped_support = any(v.get("rescue_support_groups") is not None or v.get("accepted_path_support_groups") is not None
                              for v in stats.values())
        if len(rows) != len(stats) or set(stats) != {r["species"] for r in rows}:
            raise ValueError("Model-count and BUSCO species membership differ")
        for row in rows:
            value = stats[row["species"]]
            if value["refinement_status"] != row["refinement_status"]:
                raise ValueError("Model-count analysis status differs from BUSCO")
            for key in ("prior_rescued_loci", "accepted_repair_paths", "accepted_isoform_paths"):
                count = value[key]
                if row["refinement_status"] == "not_analysed":
                    if count is not None:
                        raise ValueError("Unanalysed model counts must be unavailable, not zero")
                elif type(count) is not int or count < 0:
                    raise ValueError("Model counts must be nonnegative integers")
            support = value.get("rescue_support_counts")
            repeats = value.get("rescue_repeat_groups")
            if row["refinement_status"] == "not_analysed":
                if repeats is not None:
                    raise ValueError("Unanalysed repeat counts must be unavailable")
            elif repeats is not None:
                if (not isinstance(repeats, dict) or set(repeats) != set(REPEAT_GROUPS)
                        or any(type(c) is not int or c < 0 for c in repeats.values())
                        or sum(repeats.values()) != value["prior_rescued_loci"]):
                    raise ValueError("Repeat groups must sum to the rescued gene count")
            if row["refinement_status"] == "not_analysed":
                if support is not None or value.get("rescue_self_only_loci") is not None:
                    raise ValueError("Unanalysed rescue support counts must be unavailable")
            elif stacked_rescue:
                self_only = value.get("rescue_self_only_loci", 0)
                if (not isinstance(support, dict) or set(support) != set(RESCUE_SUPPORT)
                        or any(type(c) is not int or c < 0 for c in support.values())
                        or type(self_only) is not int or self_only < 0
                        or sum(support.values()) + self_only != value["prior_rescued_loci"]):
                    raise ValueError("Rescue support categories plus self-only loci must sum to the rescued gene count")
            path_support = value.get("accepted_path_support_counts")
            if row["refinement_status"] == "not_analysed":
                if path_support is not None:
                    raise ValueError("Unanalysed coding-path support counts must be unavailable")
            elif stacked_paths:
                if (not isinstance(path_support, dict) or set(path_support) != set(PATH_SUPPORT)
                        or any(type(c) is not int or c < 0 for c in path_support.values())
                        or sum(path_support.values()) != value["accepted_repair_paths"] + value["accepted_isoform_paths"]):
                    raise ValueError("Coding-path support categories must sum to the accepted path count")
            for key, categories, stacked, total in (
                ("rescue_support_groups", SUPPORT_GROUPS, stacked_rescue, value["prior_rescued_loci"]),
                ("accepted_path_support_groups", PATH_SUPPORT_GROUPS, stacked_paths,
                 None if row["refinement_status"] == "not_analysed" else value["accepted_repair_paths"] + value["accepted_isoform_paths"]),
            ):
                groups = value.get(key)
                if row["refinement_status"] == "not_analysed" or not stacked:
                    if groups is not None:
                        raise ValueError("Support groups must be unavailable without analysed source counts")
                elif grouped_support:
                    if (not isinstance(groups, dict) or set(groups) != set(categories)
                            or any(type(c) is not int or c < 0 for c in groups.values()) or sum(groups.values()) != total):
                        raise ValueError("S/R/P support groups must sum to the source model count")
    extra = model_changes is not None
    support_legend = stacked_rescue or stacked_paths
    margin_left = .20 if extra else .24
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11 if extra else 10, "svg.fonttype": "none"})
    figure_height = max(8, .52 * len(rows) + 5) if extra else max(6, .43 * len(rows) + 2)
    fig, axes = plt.subplots(1, 5 if extra else 3, figsize=(25 if extra else 19, figure_height),
                             gridspec_kw={"width_ratios": [1, 1, .65, .75, .9] if extra else [1, 1, .65],
                                          "wspace": .16 if extra else .12}, sharey=True)
    colors = STATUS_COLOURS
    for i, row in enumerate(rows):
        if row["refinement_status"] == "not_analysed":
            for ax in axes:
                ax.axhspan(i - .48, i + .48, color="#eef0f3", zorder=0)
        for j, phase in enumerate(("before", "after")):
            result, left = row[phase + "_result"], 0
            for category, color in zip(STATUS, colors, strict=True):
                value = 100 * result[category] / result["total"]
                axes[j].barh(i, value, left=left, color=color, height=.7)
                left += value
            axes[j].text(101, i, f'{100 * result["complete"] / result["total"]:.2f}', va="center", fontsize=9)
        delta = row["delta_complete_pp"]
        axes[2].barh(i, delta, color="#2a8d75" if delta >= 0 else "#b85250", height=.65)
        axes[2].text(.98, i, f'{row["delta_complete"]:+d} ({delta:+.2f} pp)',
                     transform=axes[2].get_yaxis_transform(), ha="right", va="center", fontsize=9)
        if extra:
            value = stats[row["species"]]
            if row["refinement_status"] == "not_analysed":
                for ax in axes[3:]:
                    ax.text(.03, i, "Not analysed", transform=ax.get_yaxis_transform(), va="center", color="#657585", fontsize=9)
            else:
                rescued = value["prior_rescued_loci"]
                repair, isoform = value["accepted_repair_paths"], value["accepted_isoform_paths"]
                rescue_y = i - .20
                if stacked_rescue:
                    offset = 0
                    support = value["rescue_support_groups"] if grouped_support else {
                        **value["rescue_support_counts"], "self_only": value.get("rescue_self_only_loci", 0)}
                    categories = SUPPORT_GROUPS if grouped_support else (*RESCUE_SUPPORT, "self_only")
                    group_colors = SUPPORT_GROUP_COLOURS if grouped_support else (*RESCUE_SUPPORT_COLOURS, RESCUE_SELF_COLOUR)
                    for category, color in zip(categories, group_colors, strict=True):
                        count = support[category]
                        axes[3].barh(rescue_y, count, left=offset, color=color, height=.32)
                        offset += count
                else:
                    axes[3].barh(rescue_y, rescued, color="#5275b5", height=.32)
                axes[3].annotate(str(rescued), (rescued, rescue_y),
                                 xytext=(4, 0), textcoords="offset points", va="center", fontsize=9)
                repeats = value.get("rescue_repeat_groups") or {**dict.fromkeys(REPEAT_GROUPS, 0), "not_assessed": rescued}
                offset = 0
                for category, color in zip(REPEAT_GROUPS, REPEAT_GROUP_COLOURS, strict=True):
                    axes[3].barh(i + .20, repeats[category], left=offset, color=color, height=.32,
                                 hatch="///" if category == "not_assessed" else None)
                    offset += repeats[category]
                axes[3].annotate(str(rescued), (rescued, i + .20), xytext=(4, 0),
                                 textcoords="offset points", va="center", fontsize=9)
                type_y = i - .20 if stacked_paths else i
                height = .32 if stacked_paths else .7
                axes[4].barh(type_y, repair, color="#187d97", height=height)
                axes[4].barh(type_y, isoform, left=repair, color="#d38b21", height=height)
                axes[4].annotate(f"{repair} / {isoform}", (repair + isoform, type_y), xytext=(4, 0),
                                 textcoords="offset points", va="center", fontsize=9)
                if stacked_paths:
                    offset = 0
                    support = value["accepted_path_support_groups"] if grouped_support else value["accepted_path_support_counts"]
                    categories = PATH_SUPPORT_GROUPS if grouped_support else PATH_SUPPORT
                    group_colors = PATH_SUPPORT_GROUP_COLOURS if grouped_support else PATH_SUPPORT_COLOURS
                    for category, color in zip(categories, group_colors, strict=True):
                        count = support[category]
                        axes[4].barh(i + .20, count, left=offset, color=color, height=.32)
                        offset += count
                    axes[4].annotate(str(offset), (offset, i + .20), xytext=(4, 0),
                                     textcoords="offset points", va="center", fontsize=9)
    labels = [r["species"].replace("_", " ") + ("  [not analysed]" if r["refinement_status"] == "not_analysed" else "") for r in rows]
    axes[0].set_yticks(range(len(rows)), labels)
    axes[0].invert_yaxis()
    titles = ["Before refinement", "After refinement", "Change in complete BUSCOs"]
    if extra:
        titles += ["Missing-gene rescue", "Accepted coding paths"]
    for ax, title in zip(axes, titles, strict=True):
        ax.set_title(title, fontweight="bold", pad=16, fontsize=11 if extra else 12)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.grid(axis="x", alpha=.15)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0)
    for ax in axes[:2]:
        ax.set_xlim(0, 116)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_xlabel("BUSCO groups (%)     C (%) at right")
    limit = max(1, max(abs(r["delta_complete_pp"]) for r in rows) * 1.2)
    axes[2].set_xlim(-limit, limit * 2.5)
    axes[2].axvline(0, color="#647383", linewidth=.7)
    axes[2].set_xlabel("Percentage points (pp)")
    if extra:
        rescued_max = max((v["prior_rescued_loci"] or 0) for v in stats.values())
        paths_max = max((v["accepted_repair_paths"] or 0) + (v["accepted_isoform_paths"] or 0) for v in stats.values())
        axes[3].set_xlim(0, max(1, rescued_max) * 1.25)
        axes[4].set_xlim(0, max(1, paths_max) * 1.65)
        for ax in axes[3:]:
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
        axes[3].set_xlabel("Upper: support; lower: repeats\nPreviously rescued gene loci")
        axes[4].set_xlabel("Upper: repair / isoform; lower: support" if stacked_paths else "Accepted paths; labels: repair / isoform")
    # Reserve footer space in inches so legends/notes stay separated for both
    # small cohorts and whole-dataset figures.
    fig.subplots_adjust(left=margin_left, right=.98, top=.90, bottom=4.6 / figure_height if extra else .14)
    fig.suptitle("Representative CDS completeness and gene-model improvement" if extra else
                 "Representative CDS completeness before and after refinement", x=margin_left, ha="left", y=.98, fontsize=17, fontweight="bold")
    identity = rows[0]["before_result"]
    fig.text(margin_left, .94, f'BUSCO {identity["busco_version"]}; {identity["lineage"]} ({identity["lineage_creation_date"]}); '
             f'n = {identity["total"]}; transcriptome mode; one representative per locus', fontsize=10)
    fig.legend([Patch(facecolor=c) for c in colors], STATUS_LABELS,
               loc="lower left", bbox_to_anchor=(margin_left, 3.25 / figure_height if extra else .065), ncol=4, frameon=False)
    if support_legend:
        if grouped_support:
            group_colors = PATH_SUPPORT_GROUP_COLOURS if stacked_paths else SUPPORT_GROUP_COLOURS
            group_labels = PATH_SUPPORT_GROUP_LABELS if stacked_paths else SUPPORT_GROUP_LABELS
        else:
            group_colors = PATH_SUPPORT_COLOURS if stacked_paths else (*RESCUE_SUPPORT_COLOURS, RESCUE_SELF_COLOUR)
            group_labels = PATH_SUPPORT_LABELS if stacked_paths else (*RESCUE_SUPPORT_LABELS, RESCUE_SELF_LABEL)
        fig.legend([Patch(facecolor=c) for c in group_colors], group_labels,
                   loc="lower left", bbox_to_anchor=(margin_left, 2.4 / figure_height), ncol=3, frameon=False,
                   title="Supporting donor groups (upper rescue and lower coding-path bars)")
        fig.legend([Patch(facecolor=c) for c in ("#187d97", "#d38b21")],
                   ["Repair coding paths", "Additional isoform paths"],
                   loc="lower left", bbox_to_anchor=(.70, 3.25 / figure_height), ncol=2, frameon=False)
    elif extra:
        fig.legend([Patch(facecolor=c) for c in ("#5275b5", "#187d97", "#d38b21")],
                   ["Previously rescued gene loci", "Repair coding paths", "Additional isoform paths"],
                   loc="lower left", bbox_to_anchor=(margin_left, 2.4 / figure_height), ncol=3, frameon=False)
    if extra:
        fig.legend([Patch(facecolor=c, hatch="///" if k == "not_assessed" else None)
                    for k, c in zip(REPEAT_GROUPS, REPEAT_GROUP_COLOURS, strict=True)], REPEAT_GROUP_LABELS,
                   loc="lower left", bbox_to_anchor=(margin_left, 1.65 / figure_height), ncol=2, frameon=False,
                   title="Repeat annotation (lower rescue bars; any CDS overlap)")
    note = "Grey rows: excluded from structural refinement; unchanged CDS are still evaluated by BUSCO.\n"
    note += "Before = refinement source CDS (including earlier rescued genes); after = selected DNA CDS, not all isoforms."
    if extra:
        note += "\nRescue counts are gene loci already in Before; repair / isoform counts are accepted paths and may share a locus."
    if support_legend:
        if grouped_support:
            note += "\nS = self-species homology; R = nearest relatives; P = phylogenetically balanced references. Target RNA is separate from S."
            note += "\nMultiple = at least two support types among S/R/P; a donor in both frozen R/P lists supplies both types."
        else:
            note += "\nBoth = support from both frozen reference groups; a donor belonging to both lists also qualifies."
            note += "\nSelf-species only = no interspecies support; mixed self/interspecies support uses the interspecies group."
    if stacked_paths:
        note += ("\nOnly = exactly one of S/R/P; other-only = no S/R/P donor. Additional unselected donors remain in the evidence."
                 if grouped_support else
                 "\nLower bars count each accepted path once by donor-group membership; other-only = no selected-group donor. Target RNA is separate.")
    elif stacked_rescue:
        note += " Labels: total loci."
    if extra:
        note += "\nRepeat overlap is advisory, not proof of TE origin; no hit does not establish a true gene. Missing annotation = not assessed."
    fig.text(margin_left, .25 / figure_height if extra else .025, note, fontsize=10)
    for suffix in ("png", "svg"):
        fig.savefig(output / ("busco_comparison." + suffix), dpi=180, facecolor="white")
    plt.close(fig)


def render_existing(report, root=None, rescue_output=None, rescue_evidence_dir=None):
    """Redraw a historical evaluation without executing its predictor again."""
    value = json.loads((report / "busco_comparison.json").read_text())
    frozen = json.loads((report / "contract.json").read_text())
    pairs = {p["species"]: p for p in frozen["pairs"]}
    rows = value["species"]
    if (value["contract"] != frozen["contract"] or len(pairs) != len(frozen["pairs"])
            or len(rows) != len(pairs) or {r["species"] for r in rows} != set(pairs)):
        raise ValueError("Comparison contract or species membership changed")
    verified_tables = 0
    for row in rows:
        pair = pairs[row["species"]]
        for phase in ("before", "after"):
            directory = report / "runs" / pair["species"] / phase
            saved = json.loads((directory / "receipt.json").read_text())
            if (saved["key"]["contract"] != frozen["contract"]
                    or saved["key"]["source_sha256"] != digest(pair[phase])
                    or saved["summary_sha256"] != digest(directory / "summary.txt")
                    or read_result(directory / "summary.txt") != row[phase + "_result"]):
                raise ValueError("Comparison input or score changed")
            # Summary-only historical publications remain renderable as such.
            if "full_table_sha256" in saved:
                if digest(directory / "full_table.tsv") != saved["full_table_sha256"]:
                    raise ValueError("Comparison full table changed")
                verified_tables += 1
        if paired_result(pair, row["before_result"], row["after_result"]) != row:
            raise ValueError("Comparison delta changed")
    changes = collect_model_changes(root, rows, rescue_output) if root is not None else None
    if changes is not None:
        collect_rescue_repeat_evidence(changes, rescue_evidence_dir)
        atomic_json(report / "model_change_summary.json", changes)
    plot_comparison(rows, report, changes)
    atomic_json(report / "rendering_provenance.json", {
        "comparison_sha256": digest(report / "busco_comparison.json"),
        "evaluation_contract": frozen["contract"], "renderer_sha256": digest(Path(__file__)),
        "verified_full_tables": verified_tables, "predictor_executed": False,
        "model_change_summary_sha256": digest(report / "model_change_summary.json") if changes is not None else None,
    })
    return rows


def evaluate(pairs, report, lineage, download_path, cpus=4, jobs=1, model_changes=None):
    if not pairs:
        raise ValueError("Empty BUSCO comparison")
    report.mkdir(parents=True, exist_ok=True)
    files = sorted(p for p in lineage.rglob("*") if p.is_file() and not p.name.startswith("."))
    if not (lineage / "dataset.cfg").is_file() or not files:
        raise ValueError("A complete local BUSCO lineage is required")
    batch = FreshDigestBatch()
    hashes = batch.read(files)
    batch.read([pair[phase] for pair in pairs for phase in ("before", "after")])
    lineage_hash = hashlib.sha256(json.dumps({str(p.relative_to(lineage)): hashes[str(p)] for p in files}, sort_keys=True).encode()).hexdigest()
    version = subprocess.check_output(["busco", "--version"], text=True).strip()
    support = Path(__file__).resolve().parent
    contract = {"schema": 1, "busco_version": version, "lineage_sha256": lineage_hash,
                "lineage_path": str(lineage), "download_path": str(download_path),
                "mode": "transcriptome", "evalue": "1e-03", "limit": 20,
                "cpus_per_job": cpus,
                "tool_sha256": {tool: digest(shutil.which(tool)) for tool in ("busco", "metaeuk", "hmmsearch")},
                "implementation": digest(Path(__file__)), "wrapper_sha256": digest(support / "gg_busco.sh"),
                "hmmsearch_wrapper_sha256": digest(support / "gg_wrapper_bin/hmmsearch")}
    atomic_json(report / "contract.json", {"contract": contract, "pairs": pairs}, immutable=True)
    def one(pair):
        before = run_one(pair, "before", report, contract, cpus)
        after = run_one(pair, "after", report, contract, cpus)
        row = paired_result(pair, before, after)
        print(json.dumps({"species": pair["species"], "delta_complete": row["delta_complete"]}), flush=True)
        return row
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        rows = list(pool.map(one, pairs))
    batch.check()
    atomic_json(report / "busco_comparison.json", {"contract": contract, "species": rows})
    with (report / "busco_comparison.tsv").open("w", newline="") as handle:
        columns = ("species", "refinement_status", "reason", "before_complete", "after_complete", "delta_complete", "delta_complete_pp", "delta_duplicated")
        writer = csv.DictWriter(handle, columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in columns if k not in ("before_complete", "after_complete")}
                            | {"before_complete": row["before_result"]["complete"], "after_complete": row["after_result"]["complete"]})
    if model_changes is not None:
        atomic_json(report / "model_change_summary.json", model_changes)
    plot_comparison(rows, report, model_changes)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Refinement publication; optional with --plot-only to include model-change counts")
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--cds-dir", type=Path, help="All dataset species, including unrefined CDS-only species")
    parser.add_argument("--rescue-output", type=Path,
                        help="Original completed rescue publication for imported inputs; otherwise inferred from the refinement plan")
    parser.add_argument("--rescue-evidence-dir", type=Path,
                        help="Separate evidence audits at DIR/SPECIES/{receipt,evidence}.json; absent species are not assessed")
    parser.add_argument("--lineage", type=Path, help="Frozen local lineage directory")
    parser.add_argument("--download-path", type=Path)
    parser.add_argument("--plot-only", action="store_true", help="Validate and redraw an existing comparison using its original evaluation contract")
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--jobs", type=int, default=1, help="Total CPU budget = jobs times cpus")
    args = parser.parse_args()
    if (args.rescue_output or args.rescue_evidence_dir) and not args.output:
        parser.error("Rescue support/evidence requires --output to bind it to the source annotation")
    if args.plot_only:
        if args.lineage or args.download_path or args.cds_dir:
            parser.error("--plot-only uses the saved evaluation; do not supply new inputs or lineage settings")
        render_existing(args.report.resolve(), args.output.resolve() if args.output else None, args.rescue_output, args.rescue_evidence_dir)
        return
    if not args.output or not args.lineage or not args.download_path:
        parser.error("--output, --lineage and --download-path are required for evaluation")
    root, report = args.output.resolve(), args.report.resolve()
    if args.cpus < 1 or args.jobs < 1:
        parser.error("CPU and job counts must be positive")
    if root == report or root in report.parents or report in root.parents:
        parser.error("Report must be separate from the immutable refinement tree")
    pairs = input_pairs(root, args.cds_dir)
    changes = collect_model_changes(root, pairs, args.rescue_output)
    collect_rescue_repeat_evidence(changes, args.rescue_evidence_dir)
    evaluate(pairs, report, args.lineage.resolve(), args.download_path.resolve(), args.cpus, args.jobs, model_changes=changes)


if __name__ == "__main__":
    main()
