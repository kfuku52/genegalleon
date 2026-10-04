#!/usr/bin/env python3
"""Audit and explicitly curate an already formatted CDS/GFF/genome pair.

This is a source-curation step before native input generation, not another
isoform selector. Only SHA-bound, individually approved records can be removed
or accepted as source exceptions. The reference genome is never rewritten.
"""

import argparse
import contextlib
import gzip
import hashlib
import io
import json
import re
import shutil
import tempfile
from collections import defaultdict
from pathlib import Path

from format_species_annotation.common import (
    build_gff_genome_seqid_map,
    first_token,
    iter_fasta_records,
    merge_coordinate_intervals,
    parse_gff_attributes,
)
from format_species_annotation.reference import genome_reference_index, validate_gff_genome_references
from format_species_writers import open_text

CONTRACT_VERSION = 1
STAT_FIELDS = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")


def fingerprint(path):
    path = Path(path).resolve(strict=True)
    before = tuple(getattr(path.stat(), key) for key in STAT_FIELDS)
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    if before != tuple(getattr(path.stat(), key) for key in STAT_FIELDS):
        raise OSError("Input changed during curation: " + str(path))
    return dict(path=str(path), sha256=digest, size=before[2], stat=before)


def gff_rows(path):
    with open_text(path, "rt") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                raise ValueError("Use a separate reference FASTA; embedded GFF FASTA is unsupported")
            if line.startswith("#") or not line.strip():
                yield line, None, None
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) != 9 or not 1 <= int(parts[3]) <= int(parts[4]):
                raise ValueError("Invalid GFF feature at line " + str(number))
            yield line, parts, parse_gff_attributes(parts[8])


def inspect_pair(species, cds, gff, genome):
    cds, gff, genome = (Path(path) for path in (cds, gff, genome))
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", species):
        raise ValueError("Invalid species prefix")
    inputs = {key: fingerprint(value) for key, value in dict(cds=cds, gff=gff, genome=genome).items()}
    lengths = {}
    for header, sequence in iter_fasta_records(cds):
        identifier = first_token(header)
        if not identifier.startswith(species + "_") or identifier in lengths or not sequence:
            raise ValueError("Expected unique, nonempty formatted CDS IDs: " + identifier)
        lengths[identifier] = len(sequence)
    if not lengths:
        raise ValueError("CDS FASTA contains no records")
    genes, parents, axes, coding, seqids = set(), defaultdict(set), defaultdict(set), [], set()
    for line, parts, attrs in gff_rows(gff):
        if parts is None:
            if line.startswith("##sequence-region "):
                fields = line.split()
                if len(fields) != 4:
                    raise ValueError("Malformed GFF sequence-region")
                seqids.add(fields[1])
            continue
        seqids.add(parts[0])
        for node in attrs.get("ID", ()):
            axes[node].add(parts[0])
            parents[node].update(attrs.get("Parent", ()))
            if parts[2].lower() == "gene":
                genes.add(node)
        if parts[2].lower() == "cds":
            if not attrs.get("Parent"):
                raise ValueError("CDS feature lacks explicit Parent")
            coding.append((parts, attrs["Parent"]))
    index = genome_reference_index(genome)
    mapping, missing = build_gff_genome_seqid_map(index, seqids)
    missing = set(missing)
    # Removing a reference must not break an ID or Parent that also belongs to
    # a retained reference. Do not guess how to split a mixed model.
    for node, references in axes.items():
        if references & missing and references - missing:
            raise ValueError("GFF ID spans present and absent references: " + node)
    for _line, parts, attrs in gff_rows(gff):
        if parts is None:
            continue
        for parent in attrs.get("Parent", ()):
            if parent in axes and any((parts[0] in missing) != (ref in missing) for ref in axes[parent]):
                raise ValueError("GFF Parent crosses the reference exclusion boundary: " + parent)

    complete = set()
    for node in list(parents):
        active, pending = set(), [(node, False)]
        while pending:
            current, finish = pending.pop()
            if finish:
                active.remove(current)
                complete.add(current)
            elif current not in complete:
                if current in active:
                    raise ValueError("Cyclic GFF Parent: " + current)
                active.add(current)
                pending.append((current, True))
                pending.extend((parent, False) for parent in parents.get(current, ()))

    def owners(node, visited=frozenset()):
        if node in visited:
            raise ValueError("Cyclic GFF Parent: " + node)
        if node in genes or node in lengths or (not parents[node] and species + "_" + node in lengths):
            return {node}
        result = set()
        for parent in parents[node]:
            result.update(owners(parent, visited | {node}))
        return result

    refs, spans, owner_sources = defaultdict(set), defaultdict(lambda: defaultdict(set)), {}
    for parts, pids in coding:
        for parent in pids:
            for owner in owners(parent):
                candidates = {owner, species + "_" + owner} & lengths.keys()
                if len(candidates) > 1:
                    raise ValueError("Ambiguous species-prefix mapping: " + owner)
                if not candidates:
                    continue
                identifier = next(iter(candidates))
                if owner_sources.setdefault(identifier, owner) != owner:
                    raise ValueError("Distinct GFF gene owners collide: " + identifier)
                refs[identifier].add(parts[0])
                spans[identifier][(parent, parts[0], parts[6])].add((int(parts[3]) - 1, int(parts[4])))
    affected = sorted(identifier for identifier, references in refs.items() if references & missing)
    for identifier in affected:
        if refs[identifier] - missing:
            raise ValueError("CDS owner spans present and absent references: " + identifier)
    span_lengths = {identifier: sorted({sum(end - start for start, end in merge_coordinate_intervals(intervals))
                                       for intervals in models.values()}) for identifier, models in spans.items()}
    report = dict(contract_version=CONTRACT_VERSION, species=species, inputs=inputs,
                  cds_records=len(lengths), missing_genome_references=sorted(missing),
                  cds_on_missing_references=affected, cds_without_gff_counterpart=sorted(lengths.keys() - refs.keys()),
                  counting_unit="Already selected, formatted CDS record; explicit GFF Parent ownership",
                  does_not_certify_genome_CDS_sequence_identity=True)
    verify_inputs(inputs)
    return report, lengths, span_lengths, mapping


def verify_inputs(inputs):
    for value in inputs.values():
        if tuple(getattr(Path(value["path"]).stat(), key) for key in STAT_FIELDS) != tuple(value["stat"]):
            raise OSError("Input changed during curation: " + value["path"])


def approved_decisions(path, report, lengths, span_lengths):
    policy = json.loads(Path(path).read_text())
    if policy.get("schema_version") != 1 or policy.get("species") != report["species"] or not policy.get("decision_basis", "").strip():
        raise ValueError("A species-specific schema-1 decision with decision_basis is required")
    expected = {key: value["sha256"] for key, value in report["inputs"].items()}
    if policy.get("input_sha256") != expected:
        raise ValueError("Decision manifest does not match exact CDS/GFF/genome hashes")
    excluded, retained = set(), {}
    for row in policy.get("records", ()):
        identifier, action, reason = row.get("cds_id"), row.get("action"), row.get("reason")
        if identifier not in lengths or identifier in excluded or identifier in retained:
            raise ValueError("Missing or duplicate approved CDS ID: " + str(identifier))
        if action == "exclude" and reason == "missing_genome_reference":
            excluded.add(identifier)
        elif action == "retain_and_flag" and reason == "missing_gff_counterpart":
            if identifier not in report["cds_without_gff_counterpart"]:
                raise ValueError("Approved missing-GFF exception is no longer observed: " + identifier)
            retained[identifier] = row
        elif action == "retain_and_flag" and reason == "coding_span_conflict":
            observed = span_lengths.get(identifier, ())
            if (row.get("cds_length") != lengths[identifier] or row.get("gff_coding_span_length") not in observed
                    or lengths[identifier] in observed):
                raise ValueError("Approved coding-span conflict does not match observed evidence: " + identifier)
            retained[identifier] = row
        else:
            raise ValueError("Unsupported explicit CDS decision: " + str(identifier))
    if excluded != set(report["cds_on_missing_references"]):
        raise ValueError("Approved exclusion IDs differ from all observed CDS on absent references")
    if report["missing_genome_references"] and policy.get("exclude_missing_reference_annotations") is not True:
        raise ValueError("Removing annotations on absent references requires explicit approval")
    if set(report["cds_without_gff_counterpart"]) != {key for key, row in retained.items() if row["reason"] == "missing_gff_counterpart"}:
        raise ValueError("Unapproved CDS without GFF counterpart")
    if len(excluded) == len(lengths):
        raise ValueError("Curation would remove every CDS record")
    return policy, excluded, retained


@contextlib.contextmanager
def gzip_text(path):
    with Path(path).open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                yield text


def curate_pair(species, cds, gff, genome, decision_manifest, output_dir):
    """Publish a new local source pair only after all checks succeed."""
    policy_fingerprint = fingerprint(decision_manifest)
    report, lengths, span_lengths, mapping = inspect_pair(species, cds, gff, genome)
    policy, excluded, retained = approved_decisions(decision_manifest, report, lengths, span_lengths)
    target = Path(output_dir).resolve()
    if target.exists():
        raise FileExistsError("Curation destination already exists: " + str(target))
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".gg-paired-curation-", dir=target.parent) as temporary:
        root = Path(temporary)
        cds_output, gff_output = root / (species + "_cds.fa.gz"), root / (species + "_annotations.gff.gz")
        # Copy decoded lines literally; retained sequences, descriptions, case,
        # wrapping and IDs must not change during source curation.
        keep, removed_features = False, 0
        with open_text(cds, "rt") as source, gzip_text(cds_output) as dest:
            for line in source:
                if line.startswith(">"):
                    keep = first_token(line[1:]) not in excluded
                if keep:
                    dest.write(line)
        with gzip_text(gff_output) as dest:
            for line, parts, _attrs in gff_rows(gff):
                if parts is not None:
                    if parts[0] not in mapping:
                        removed_features += 1
                        continue
                    if mapping[parts[0]] != parts[0]:
                        parts[0] = mapping[parts[0]]
                        line = "\t".join(parts) + "\n"
                elif line.startswith("##sequence-region "):
                    fields = line.split()
                    if len(fields) != 4:
                        raise ValueError("Malformed GFF sequence-region")
                    if fields[1] not in mapping:
                        continue
                    if mapping[fields[1]] != fields[1]:
                        fields[1] = mapping[fields[1]]
                        line = " ".join(fields) + "\n"
                dest.write(line)
        checked = validate_gff_genome_references(gff_output, genome)
        if not checked:
            raise ValueError("Curation would remove every GFF feature")
        observed = {first_token(header): len(sequence) for header, sequence in iter_fasta_records(cds_output)}
        if observed != {key: value for key, value in lengths.items() if key not in excluded}:
            raise ValueError("Curated CDS identity/count verification failed")
        verify_inputs(report["inputs"])
        verify_inputs({"policy": policy_fingerprint})
        report.update(decision_manifest=policy_fingerprint, decision_basis=policy["decision_basis"],
                      excluded_cds_ids=sorted(excluded), retained_source_exceptions=list(retained.values()),
                      remaining_cds_records=len(observed), excluded_gff_features=removed_features,
                      checked_gff_features=checked, original_sources_modified=False,
                      cds_output=fingerprint(cds_output), gff_output=fingerprint(gff_output))
        for key in ("cds_output", "gff_output"):
            report[key]["path"] = str(target / Path(report[key]["path"]).name)
            report[key].pop("stat")
        (root / "curation.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        # Source hashes and curation receipt stay together. Native input
        # generation can use these outputs through its existing local manifest.
        target.mkdir()  # Exclusive creation; never merge with another writer.
        try:
            for path in sorted(root.iterdir(), key=lambda path: path.name == "curation.json"):
                path.rename(target / path.name)
        except BaseException:
            shutil.rmtree(target)
            raise
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("audit", "curate"))
    for key in ("species", "cds", "gff", "genome"):
        parser.add_argument("--" + key, required=True)
    parser.add_argument("--decision-manifest")
    parser.add_argument("--output-dir")
    parser.add_argument("--report", help="Audit JSON destination (audit mode only).")
    args = parser.parse_args()
    if args.mode == "curate" and (not args.decision_manifest or not args.output_dir):
        parser.error("curate requires --decision-manifest and --output-dir")
    if args.mode == "audit" and (args.decision_manifest or args.output_dir):
        parser.error("Use curate to apply explicit decisions")
    try:
        if args.mode == "curate":
            report = curate_pair(args.species, args.cds, args.gff, args.genome, args.decision_manifest, args.output_dir)
        else:
            report = inspect_pair(args.species, args.cds, args.gff, args.genome)[0]
            if args.report:
                with Path(args.report).open("x") as handle:
                    json.dump(report, handle, indent=2, sort_keys=True)
                    handle.write("\n")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.error(str(exc))
    print(json.dumps(dict(species=args.species, cds_records=report["cds_records"],
                          missing_genome_references=len(report["missing_genome_references"]),
                          cds_on_missing_references=len(report["cds_on_missing_references"]),
                          cds_without_gff_counterpart=len(report["cds_without_gff_counterpart"]),
                          remaining_cds_records=report.get("remaining_cds_records"))))


if __name__ == "__main__":
    main()
