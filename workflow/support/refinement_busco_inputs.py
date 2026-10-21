"""Prepare immutable BUSCO protein views from frozen refinement admission.

No catalogue or genomic reconstruction runs here. Native source records use
their published annotated coding paths; the selected view copies its admitted
protein bytes. CDS-only species have an explicit, audited frame-zero policy.
"""
from __future__ import annotations

import csv
import hashlib
import json
import shutil
import sqlite3
import tempfile
from pathlib import Path

from Bio.Data import CodonTable
from Bio.Seq import Seq

try:
    from fasta_sequence_store import exclusive_lock, fasta_records
    from input_generation_array_state import FreshDigestBatch, atomic_json, digest, safe_component
except ImportError:
    from .fasta_sequence_store import exclusive_lock, fasta_records
    from .input_generation_array_state import FreshDigestBatch, atomic_json, digest, safe_component

ADMITTED_PROTEIN_POLICY = {
    "schema": 1,
    "native": "frozen_source_mapping_and_existing_exportable_usable_admission",
    "native_partial": "retain_published_annotated_phase_and_complete_codons",
    "native_ambiguous_identity": "published_same_coding_path_baseline_or_exclude",
    "after": "exact_effective_protein_and_translation_admission",
    "cds_only": {"frame": 0, "genetic_code": "explicit_species_code",
                 "terminal_stop": "remove_one", "partial_terminal_codon": "exclude",
                 "internal_stop": "exclude", "invalid_base": "exclude",
                 "ambiguous_codon": "IUPAC_translation_or_X",
                 "partial_start_or_stop": "retain_without_inventing_residues",
                 "context_dependent_genetic_code": "reject"},
    "identifier": "source_fasta_id_for_source_phases_effective_gene_id_after",
    "duplicate_identifier": "reject",
    "exclusion": "record_every_excluded_source_identifier_and_reason",
}


class _Snapshots:
    def __init__(self):
        self.batch = FreshDigestBatch()
        self.hashes = {}

    def bind(self, path, expected=None):
        path = Path(path).resolve()
        sha = self.batch.read([path])[str(path)]
        if expected is not None and sha != expected:
            raise ValueError("Frozen protein input changed: " + str(path))
        self.hashes[str(path)] = sha
        return sha

    def json(self, path, expected=None):
        sha = self.bind(path, expected)
        raw = Path(path).read_bytes()
        if hashlib.sha256(raw).hexdigest() != sha:
            raise ValueError("Protein input changed while reading: " + str(path))
        self.batch.check()
        return json.loads(raw)

    def publication(self, directory, plan_sha, members):
        directory = Path(directory)
        receipt_path = directory / "receipt.json"
        receipt = self.json(receipt_path)
        if (not isinstance(receipt, dict) or receipt.get("key", {}).get("plan") != plan_sha
                or not isinstance(receipt.get("files"), dict)):
            raise ValueError("Protein evidence belongs to another refinement plan: " + str(directory))
        paths = [receipt_path]
        for member in members:
            expected = receipt["files"].get(member)
            if not isinstance(expected, str):
                raise ValueError("Protein evidence lacks published member: " + member)
            path = directory / member
            self.bind(path, expected)
            paths.append(path)
        return receipt, paths


def _records(path):
    seen = set()
    for identifier, header, dna in fasta_records(Path(path)):
        if identifier in seen:
            raise ValueError("Duplicate protein-source identifier: " + identifier)
        seen.add(identifier)
        yield identifier, header, dna.upper()


def _quality_reasons(quality):
    reasons = [key for key in ("internal_stop", "ambiguous", "invalid_base", "phase_conflict",
               "phase_unknown", "phase_unresolved", "sequence_mismatch", "translation_uncertain",
               "annotated_pseudogene", "translation_exception", "sequence_exception") if quality.get(key)]
    for key in ("structure_problem", "annotated_exception"):
        if quality.get(key):
            reasons.append(key + ":" + str(quality[key]))
    return reasons or ["native_translation_not_admitted"]


def _native_catalogue(root, species, source, plan, frozen, *, load_candidates=True):
    directory = root / "catalog" / species
    _receipt, paths = frozen.publication(directory, frozen.hashes[str(root / "plan.json")],
                                         ["catalog_metadata.json", "loci.jsonl"])
    metadata = frozen.json(directory / "catalog_metadata.json")
    if metadata.get("species") != species or metadata.get("genetic_code") != source["genetic_code"]:
        raise ValueError("Native catalogue species/code differs: " + species)
    for role, catalogue_role in (("fasta", "cds"), ("gff", "gff"), ("genome", "genome")):
        path = Path(source[role]).resolve()
        expected = plan["request"]["files"][str(path)]
        frozen.bind(path, expected)
        if metadata.get("sources", {}).get(catalogue_role, {}).get("sha256") != expected:
            raise ValueError("Native catalogue source differs: " + species + "/" + role)
        paths.append(path)
    mappings = {}
    wanted = set()
    for row in metadata["fasta_mapping"]:
        identifier = row["source_fasta_id"]
        if identifier in mappings:
            raise ValueError("Duplicate native source mapping: " + identifier)
        mappings[identifier] = row
        wanted.update(row.get("candidate_ids", []))
        wanted.update(row.get("sequence_agreeing_candidate_ids", []))
    candidates = {}
    if not load_candidates:
        return mappings, candidates, paths
    with (directory / "loci.jsonl").open() as handle:
        for line in handle:
            locus = json.loads(line)
            for candidate in locus["candidates"]:
                identifier = candidate["candidate_id"]
                if identifier not in wanted:
                    continue
                if identifier in candidates:
                    raise ValueError("Duplicate native candidate identity: " + identifier)
                candidates[identifier] = (candidate, {
                    key: locus.get(key) for key in ("gene_id", "ambiguous_coordinates",
                        "source_baseline_candidate_id", "source_baseline_coding_candidate_id")})
    frozen.batch.check()
    if wanted - candidates.keys():
        raise ValueError("Native source mapping references absent candidates: " + species)
    return mappings, candidates, paths


def _native_record(identifier, dna, mappings, candidates):
    mapping = mappings.get(identifier)
    if mapping is None:
        return None, ["source_identifier_absent_from_native_mapping"], {}
    binding = {"mapping": mapping}
    if hashlib.sha256(dna.encode()).hexdigest() != mapping.get("source_cds_sha256"):
        return None, ["source_sequence_not_bound_to_native_catalogue"], binding
    chosen = mapping.get("candidate_ids", []) if mapping.get("mapping_status") == "mapped" else []
    if len(chosen) != 1:
        agreeing = mapping.get("sequence_agreeing_candidate_ids", [])
        if agreeing:
            paths = [candidates[key] for key in agreeing]
            baselines = {locus.get("source_baseline_coding_candidate_id") for _candidate, locus in paths}
            if (len({(c["gene_id"], c["coding_key"]) for c, _locus in paths}) == 1
                    and len(baselines) == 1 and next(iter(baselines)) in agreeing):
                chosen = [next(iter(baselines))]
        if len(chosen) != 1:
            return None, [mapping.get("exclusion_reason") or "native_mapping_" + mapping.get("mapping_status", "unmapped")], binding
    candidate, locus = candidates[chosen[0]]
    quality = candidate["quality"]
    binding.update(candidate_id=chosen[0], gene_id=candidate["gene_id"], quality=quality,
                   source_convention=mapping.get("source_convention", ""))
    exportable = (bool(candidate["cds"]) and bool(candidate["blocks"])
                  and not locus.get("ambiguous_coordinates") and not quality.get("sequence_mismatch")
                  and not quality.get("structure_problem"))
    admitted = exportable and quality.get("usable", False) and bool(candidate["protein"]) and "*" not in candidate["protein"]
    if not admitted:
        reasons = ["unreconstructible_or_ambiguous_coordinates"] if locus.get("ambiguous_coordinates") else _quality_reasons(quality)
        return None, reasons, binding
    return candidate["protein"], [], binding


def _frame_zero_record(dna, code):
    reasons = []
    if not dna:
        reasons.append("empty_sequence")
    if len(dna) % 3:
        reasons.append("partial_terminal_codon")
    if set(dna) - set("ACGTRYSWKMBDHVN"):
        reasons.append("invalid_base")
    protein = str(Seq(dna).translate(table=code)) if not reasons else ""
    if "*" in protein.removesuffix("*"):
        reasons.append("internal_stop")
    protein = protein.removesuffix("*")
    if not protein and not reasons:
        reasons.append("empty_protein")
    return (None if reasons else protein), reasons, {"basis": "declared_frame_zero"}


def _source_view(source, output, code, native=None):
    records, excluded, included = [], [], 0
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, _header, dna in _records(source):
            protein, reasons, binding = (_native_record(identifier, dna, *native) if native is not None
                                         else _frame_zero_record(dna, code))
            row = {"source_fasta_id": identifier, "source_cds_sha256": hashlib.sha256(dna.encode()).hexdigest(),
                   "status": "included" if protein is not None else "excluded", **binding}
            if protein is None:
                row["reasons"] = reasons
                excluded.append({"source_fasta_id": identifier, "identifier": identifier,
                                 "reason": reasons[0], "reasons": reasons})
            else:
                handle.write(">" + identifier + "\n" + protein + "\n")
                included += 1
            records.append(row)
    return records, excluded, included, {}


def _admission_index(effective, destination, species, frozen):
    """One global TSV pass; use disk storage rather than cohort-sized rows."""
    connection = sqlite3.connect(destination)
    try:
        connection.execute('CREATE TABLE admission (species TEXT, gene_id TEXT, record TEXT, PRIMARY KEY(species,gene_id))')
        with connection, (effective / "translation_admission.tsv").open(newline="") as handle:
            for record in csv.DictReader(handle, delimiter="\t"):
                if record["species"] not in species:
                    continue
                if record["status"] not in {"included", "excluded"}:
                    raise ValueError("Invalid effective translation admission status")
                try:
                    connection.execute('INSERT INTO admission VALUES (?,?,?)',
                                       (record["species"], record["gene_id"], json.dumps(record)))
                except sqlite3.IntegrityError as error:
                    raise ValueError("Duplicate effective translation admission") from error
        frozen.batch.check()
    finally:
        connection.close()
    return destination


def _selected_view(species, source, output, effective, row, code, frozen, admission_database):
    proteins = dict((identifier, protein) for identifier, _header, protein in _records(row["protein"]))
    admission = {}
    connection = sqlite3.connect(Path(admission_database).resolve().as_uri() + '?mode=ro', uri=True)
    try:
        for (raw,) in connection.execute('SELECT record FROM admission WHERE species=? ORDER BY rowid', (species,)):
            record = json.loads(raw)
            identifier = record["gene_id"]
            if identifier in admission or record["status"] not in {"included", "excluded"}:
                raise ValueError("Invalid effective translation admission: " + species)
            quality = json.loads(record["quality"])
            if quality.get("genetic_code") != code:
                raise ValueError("Effective translation code differs: " + species)
            admission[identifier] = {**record, "quality": quality}
    finally:
        connection.close()
    if set(proteins) != {identifier for identifier, record in admission.items() if record["status"] == "included"}:
        raise ValueError("Effective protein membership differs from translation admission: " + species)
    if any(not protein or "*" in protein for protein in proteins.values()):
        raise ValueError("Effective admitted protein is empty or contains a stop: " + species)
    records, excluded, seen = [], [], set()
    for identifier, _header, dna in _records(source):
        seen.add(identifier)
        if identifier not in admission:
            raise ValueError("Effective CDS lacks translation admission: " + identifier)
        record = admission[identifier]
        records.append({"source_fasta_id": identifier, "candidate_id": record["candidate_id"],
                        "status": record["status"], "quality": record["quality"],
                        "source_cds_sha256": hashlib.sha256(dna.encode()).hexdigest()})
        if record["status"] == "excluded":
            reasons = _quality_reasons(record["quality"])
            excluded.append({"source_fasta_id": identifier, "identifier": identifier,
                             "reason": reasons[0], "reasons": reasons})
    if set(proteins) - seen:
        raise ValueError("Effective protein lacks its source CDS: " + species)
    shutil.copyfile(row["protein"], output)
    frozen.batch.check()
    return records, excluded, len(proteins), {
        "annotation_excluded_records": [record for identifier, record in admission.items() if identifier not in seen],
        "effective_protein_sha256": row["protein_sha256"],
        "translation_admission_sha256": frozen.hashes[str(effective / "translation_admission.tsv")],
    }


def prepare_admitted_protein_inputs(pairs, root, destination, genetic_codes):
    """Prepare/verify one immutable publication and return phase descriptors."""
    root, destination = Path(root).resolve(), Path(destination).resolve()
    pairs = list(pairs)
    names = [pair["species"] for pair in pairs]
    if not pairs or len(set(names)) != len(names) or not all(safe_component(name) for name in names):
        raise ValueError("Protein preparation needs unique safe species identifiers")
    frozen = _Snapshots()
    plan_path = root / "plan.json"
    plan = frozen.json(plan_path)
    plan_sha = frozen.hashes[str(plan_path)]
    sources = plan["request"]["sources"]
    effective = root / "effective"
    analysed = [name for name in names if name in sources]
    effective_rows = {}
    if analysed:
        effective_members = ["inputs.tsv", "translation_admission.tsv"]
        effective_members += ["species_protein/" + name + ".fa" for name in analysed]
        _receipt, _paths = frozen.publication(effective, plan_sha, effective_members)
        with (effective / "inputs.tsv").open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["species"] in effective_rows:
                    raise ValueError("Duplicate effective species")
                effective_rows[row["species"]] = row
    frozen.bind(Path(__file__))
    protected = [effective, *(root / "catalog" / name for name in analysed)]
    if (destination == root or root.is_relative_to(destination)
            or any(destination == path or destination.is_relative_to(path) for path in protected)):
        raise ValueError("Protein outputs overlap a frozen refinement publication")
    work = []
    for pair in pairs:
        name = pair["species"]
        code = genetic_codes.get(name)
        if type(code) is not int or code not in CodonTable.unambiguous_dna_by_id:
            raise ValueError("Explicit supported genetic code required: " + name)
        native = None
        species_paths = [plan_path, Path(__file__).resolve()]
        if name in sources:
            if pair.get("refinement_status") != "analysed" or code != sources[name]["genetic_code"]:
                raise ValueError("Native protein scope/code differs: " + name)
            paths = _native_catalogue(root, name, sources[name], plan, frozen, load_candidates=False)[2]
            native = True
            species_paths += paths
        elif pair.get("refinement_status") != "not_analysed":
            raise ValueError("Protein species lacks native annotation: " + name)
        elif set(CodonTable.unambiguous_dna_by_id[code].stop_codons) & CodonTable.unambiguous_dna_by_id[code].forward_table.keys():
            raise ValueError("Context-dependent CDS-only genetic code: " + name)
        for phase in (("pre_rescue", "before", "after") if "pre_rescue" in pair else ("before", "after")):
            source_path = Path(pair[phase]).resolve()
            if source_path == destination or source_path.is_relative_to(destination):
                raise ValueError("Protein destination overlaps source CDS")
            source_sha = frozen.bind(source_path)
            paths = [*species_paths, source_path]
            selected = None
            if native is not None and phase == "before" and source_path != Path(sources[name]["fasta"]).resolve():
                raise ValueError("Before protein source differs from native catalogue: " + name)
            if native is not None and phase == "after":
                row = effective_rows[name]
                if (source_path != Path(row["cds"]).resolve() or source_sha != row["cds_sha256"]
                        or int(row["genetic_code"]) != code
                        or Path(row["protein"]).resolve() != effective / "species_protein" / (name + ".fa")):
                    raise ValueError("After protein inputs differ from effective publication: " + name)
                frozen.bind(row["protein"], row["protein_sha256"])
                paths += [effective / member for member in ("receipt.json", "inputs.tsv", "translation_admission.tsv")]
                paths.append(Path(row["protein"]))
                selected = row
            work.append((name, phase, source_path, source_sha, code, native, selected, paths))
    frozen.batch.check()
    key = {"schema": 1, "policy": ADMITTED_PROTEIN_POLICY, "inputs": dict(frozen.hashes),
           "pairs": pairs, "genetic_codes": {name: genetic_codes[name] for name in names}}
    destination.parent.mkdir(parents=True, exist_ok=True)
    with exclusive_lock(destination.parent / ("." + destination.name + ".lock")):
        if destination.exists():
            receipt = frozen.json(destination / "receipt.json")
            if receipt.get("key") != key:
                raise ValueError("Admitted protein publication has different frozen inputs; use a new destination")
            required = {"manifest.json"}
            for name, phase, *_rest in work:
                required.update({name + "/" + phase + suffix for suffix in (".faa", ".admission.json")})
            if set(receipt.get("files", {})) != required:
                raise ValueError("Admitted protein publication has incomplete member scope")
            for member, sha in receipt["files"].items():
                if Path(member).is_absolute() or ".." in Path(member).parts:
                    raise ValueError("Unsafe protein publication member")
                frozen.bind(destination / member, sha)
            manifest = frozen.json(destination / "manifest.json", receipt["files"].get("manifest.json"))
            if set(manifest) != set(names):
                raise ValueError("Admitted protein manifest has different species scope")
            for name in names:
                if set(manifest[name]) != {phase for species, phase, *_rest in work if species == name}:
                    raise ValueError("Admitted protein manifest has different phase scope")
            for name, phase, source_path, source_sha, code, *_rest in work:
                descriptor = manifest[name][phase]
                protein_path = destination / name / (phase + ".faa")
                audit_path = destination / name / (phase + ".admission.json")
                audit = frozen.json(audit_path)
                if (descriptor.get("protein_path") != str(protein_path) or descriptor.get("audit_path") != str(audit_path)
                        or descriptor.get("protein_sha256") != frozen.hashes[str(protein_path)]
                        or descriptor.get("audit_sha256") != frozen.hashes[str(audit_path)]
                        or descriptor.get("source_sha256") != source_sha or descriptor.get("policy") != ADMITTED_PROTEIN_POLICY
                        or audit.get("species") != name or audit.get("phase") != phase
                        or audit.get("genetic_code") != code or audit.get("source_sha256") != source_sha
                        or audit.get("protein_sha256") != descriptor["protein_sha256"]
                        or audit.get("policy") != ADMITTED_PROTEIN_POLICY):
                    raise ValueError("Admitted protein descriptor differs from its frozen audit")
            frozen.batch.check()
            return manifest
        temporary = Path(tempfile.mkdtemp(prefix=".admitted-proteins-", dir=destination.parent))
        try:
            manifest = {}
            admission_database = (_admission_index(effective, temporary / '.translation_admission.sqlite3',
                                                   set(analysed), frozen) if analysed else None)
            current_name, native_records = None, None
            for name, phase, source_path, source_sha, code, native, selected, paths in work:
                if name != current_name:
                    current_name = name
                    # Release the preceding species before constructing the
                    # next one, including during the RHS catalogue allocation.
                    native_records = None
                    native_records = (_native_catalogue(root, name, sources[name], plan, frozen)[:2]
                                      if native else None)
                folder = temporary / name
                folder.mkdir(exist_ok=True)
                protein, audit_path = folder / (phase + ".faa"), folder / (phase + ".admission.json")
                records, excluded, included, extras = (
                    _selected_view(name, source_path, protein, effective, selected, code, frozen, admission_database)
                    if selected is not None
                    else _source_view(source_path, protein, code, native_records))
                audit = {"schema": 1, "status": "translated", "policy": ADMITTED_PROTEIN_POLICY,
                         "species": name, "phase": phase, "genetic_code": code,
                         "source_path": str(source_path), "source_sha256": source_sha,
                         "protein_sha256": digest(protein),
                         "counts": {"records": len(records), "included": included, "excluded": len(excluded)},
                         "excluded_records": excluded, "records": records,
                         "inputs": {str(Path(path).resolve()): frozen.hashes[str(Path(path).resolve())] for path in paths},
                         "implementation_sha256": frozen.hashes[str(Path(__file__).resolve())],
                         "native_implementation": plan["request"].get("implementation") if native else None, **extras}
                atomic_json(audit_path, audit)
                manifest.setdefault(name, {})[phase] = {
                    "protein_path": str(destination / name / protein.name), "audit_path": str(destination / name / audit_path.name),
                    "protein_sha256": audit["protein_sha256"], "audit_sha256": digest(audit_path),
                    "source_sha256": source_sha, "policy": ADMITTED_PROTEIN_POLICY}
            if admission_database is not None:
                admission_database.unlink()
            atomic_json(temporary / "manifest.json", manifest)
            files = {str(path.relative_to(temporary)): digest(path) for path in temporary.rglob("*") if path.is_file()}
            atomic_json(temporary / "receipt.json", {"key": key, "files": files})
            frozen.batch.check()
            temporary.rename(destination)
            return manifest
        finally:
            if temporary.exists():
                shutil.rmtree(temporary)
