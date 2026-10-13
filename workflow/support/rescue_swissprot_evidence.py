#!/usr/bin/env python3
"""Advisory Swiss-Prot homology for published missing-gene rescue candidates.

Only published rescue CDSs are searched. No model, admission decision or BUSCO
score is changed. Other-protein support means absence of a TE annotation in the
matched Swiss-Prot entry, not proof that the candidate is a functional host gene.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import inspect
import json
import math
import re
import shutil
import sqlite3
import subprocess
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path

from Bio.Data import CodonTable
from Bio.Seq import Seq

try:
    from fasta_sequence_store import fasta_records, open_text
    from input_generation_array_state import FreshDigestBatch, atomic_json, digest, safe_component
    from rescue_gene_models import exclusive_lock, stage
    from rescue_model_evidence import attributes, read_json_snapshot
except ImportError:
    from .fasta_sequence_store import fasta_records, open_text
    from .input_generation_array_state import FreshDigestBatch, atomic_json, digest, safe_component
    from .rescue_gene_models import exclusive_lock, stage
    from .rescue_model_evidence import attributes, read_json_snapshot

GROUPS = ("te_only", "other_only", "both", "no_informative_hit", "not_assessed")
COLOURS = ("#b45158", "#477b80", "#9a6ab2", "#c7cdd4", "#edf0f3")
LABELS = ("TE-related support only", "Other protein support only", "TE-related + other support",
          "No informative support", "Not assessed")
FIELDS = ("query,target,pident,alnlen,evalue,bits,qlen,tlen,qstart,qend,tstart,tend,qaln,taln")
DEFAULTS = {"evalue": 1e-5, "query_coverage": .5, "target_coverage": .5,
            "minimum_alignment": 50, "score_fraction": .9, "max_hits": 50, "sensitivity": 7.5,
            "search_evalue": 1e-5}
NO_SUPPORT_REASONS = ("no_returned_hits", "weak_hit", "short_hit", "partial_hit", "annotation_unknown")
TE_PATTERN = re.compile(r"\b(?:transpos(?:ase|on|able)|retrotranspos\w*|retroposon)\b", re.I)
UNINFORMATIVE = re.compile(r"\b(?:uncharacteri[sz]ed|hypothetical|putative protein|unknown function)\b", re.I)
METHOD = ("Swiss-Prot protein homology; qualifying alignments satisfy the recorded E-value, paired query/target "
          "coverage and length thresholds, and are within the recorded score fraction of the best qualifying hit. "
          "TE support uses explicit TE protein names, the Transposable element keyword or transposase activity. "
          "TE silencing/regulation GO terms and the Transposition process keyword alone do not imply TE origin. Other support uses informative "
          "entries lacking these annotations. Each gene locus counts once across its published coding sequences. "
          "Heuristic support categories are advisory, not calibrated TE probabilities or functional-gene calls. "
          "Partial TE matches remain in the hit table; no informative support does not exclude TE origin.")


def encoded_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def protein(sequence, code):
    sequence = sequence.upper()
    table = CodonTable.unambiguous_dna_by_id[code]
    if (len(sequence) < 6 or len(sequence) % 3 or set(sequence) - set("ACGT")
            or sequence[:3] not in table.start_codons or sequence[-3:] not in table.stop_codons):
        raise ValueError("Published rescue candidate lacks an intact code-compatible CDS")
    if set(table.stop_codons) & set(table.forward_table):
        return None  # Context-dependent translation is not inferred by this audit.
    result = str(Seq(sequence).translate(table=code))
    if not result.endswith("*") or "*" in result[:-1]:
        raise ValueError("Published rescue candidate contains an internal stop")
    return result[:-1]


def candidates(root, boundary):
    """Read small finalized publications, not genome-wide prediction arrays."""
    root = Path(root).resolve()
    plan, plan_hash = read_json_snapshot(root / "plan.json")
    receipt, receipt_hash = read_json_snapshot(root / "augmented/receipt.json")
    expected = {str(root / "plan.json"): plan_hash, str(root / "augmented/receipt.json"): receipt_hash}
    if receipt.get("key", {}).get("plan") != plan_hash:
        raise ValueError("Augmented rescue receipt belongs to another plan")
    augmented = root / "augmented"
    index = augmented / "inputs.tsv"
    expected[str(index)] = receipt["files"]["inputs.tsv"]
    if boundary.read(expected) != expected:
        raise ValueError("Frozen rescue publication changed")
    with index.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    seen, records, bindings = set(), [], {}
    for row in rows:
        name = row["species"]
        if not safe_component(name) or name in seen or name not in plan["request"]["sources"]:
            raise ValueError("Invalid/duplicate published rescue species")
        seen.add(name)
        paths = {key: Path(row[key]).resolve() for key in ("cds", "gff")}
        for path in paths.values():
            relative = str(path.relative_to(augmented))
            expected[str(path)] = receipt["files"][relative]
        if boundary.read(paths.values()) != {str(p): expected[str(p)] for p in paths.values()}:
            raise ValueError("Frozen rescue publication changed")
        genes, owners = set(), {}
        transcripts = []
        with open_text(paths["gff"]) as handle:
            for line in handle:
                if line.strip() == "##FASTA":
                    break
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip().split("\t")
                if len(f) != 9:
                    raise ValueError("Invalid published rescue GFF")
                if f[1] != "genegalleon_rescue" or f[2] not in {"gene", "mRNA", "transcript"}:
                    continue
                a = attributes(f[8])
                identifier = a.get("ID", "")
                if not identifier or identifier in owners:
                    raise ValueError("Duplicate/missing rescued GFF identity")
                owners[identifier] = identifier if f[2] == "gene" else None
                if f[2] == "gene":
                    genes.add(identifier)
                else:
                    transcripts.append((identifier, a.get("Parent", "")))
        for identifier, parent in transcripts:
            if parent not in genes:
                raise ValueError("Rescued transcript has no unique rescued gene parent")
            owners[identifier] = parent
        if len(genes) != int(row["rescued_models"]):
            raise ValueError("Rescued gene count differs from publication")
        code = plan["request"]["sources"][name]["genetic_code"]
        found, fasta_ids = set(), set()
        for identifier, _, sequence in fasta_records(paths["cds"]):
            if identifier in fasta_ids:
                raise ValueError("Duplicate CDS identity in published rescue FASTA")
            fasta_ids.add(identifier)
            if identifier not in owners:
                continue
            gene = owners[identifier]
            found.add(gene)
            pep = protein(sequence, code)
            records.append({"species": name, "gene_id": gene, "cds_id": identifier, "genetic_code": code,
                            "cds_sha256": hashlib.sha256(sequence.upper().encode()).hexdigest(),
                            "protein_sha256": hashlib.sha256(pep.encode()).hexdigest() if pep else None,
                            "protein": pep})
        if found != genes:
            raise ValueError("Published CDS is missing a rescued gene")
        bindings[name] = {"source_gff_sha256": expected[str(paths["gff"])],
                          "source_cds_sha256": expected[str(paths["cds"])], "gene_ids": sorted(genes)}
    if set(receipt["key"]["rescue_receipts"]) != seen:
        raise ValueError("Augmented rescue species membership differs from receipt")
    if boundary.read(expected) != expected:
        raise ValueError("Frozen rescue publication changed")
    return records, bindings, expected, plan_hash, receipt_hash


def accession(identifier):
    if identifier.startswith("sp|"):
        return identifier.split("|")[1]
    return identifier


def reference_annotations(prefix, metadata, targets):
    wanted = {accession(t) for t in targets}
    descriptions, annotations = {}, {}
    for identifier, header, sequence in fasta_records(Path(str(prefix) + ".pep")):
        if accession(identifier) in wanted:
            if accession(identifier) in descriptions:
                raise ValueError("Duplicate Swiss-Prot FASTA accession")
            descriptions[accession(identifier)] = (header.split(None, 1)[1].split(" OS=", 1)[0], sequence)
    with open_text(metadata) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"accession", "keywords", "go_terms"} <= set(reader.fieldnames or []):
            raise ValueError("Swiss-Prot metadata lacks required annotation columns")
        for row in reader:
            if row["accession"] in wanted:
                if row["accession"] in annotations:
                    raise ValueError("Duplicate Swiss-Prot metadata accession")
                annotations[row["accession"]] = row
    if wanted != set(descriptions) or wanted != set(annotations):
        raise ValueError("Swiss-Prot sequence/metadata is incomplete or from a different release")
    result = {}
    for identifier in wanted:
        name, sequence = descriptions[identifier]
        row = annotations[identifier]
        keywords = {k.strip().lower() for k in row["keywords"].split(";")}
        go_terms = {term.strip().lower() for term in row["go_terms"].split(";")}
        named_te = bool(TE_PATTERN.search(name)) and not re.search(r"\b(?:silencing|suppression|suppressor|regulator|repressor|restriction)\b", name, re.I)
        basis = (["protein_name"] if named_te else []) + (["keyword:Transposable element"] if "transposable element" in keywords else [])
        if "transposase activity" in go_terms:
            basis.append("GO:transposase activity")
        te = bool(basis)
        informative = not UNINFORMATIVE.search(name) or bool(row["go_terms"] or row.get("ec_numbers"))
        result[identifier] = {"accession": identifier, "protein_name": name, "keywords": row["keywords"],
                              "go_terms": row["go_terms"], "organism": row.get("organism", ""),
                              "taxid": row.get("taxid", ""), "target_length": len(sequence), "sequence": sequence,
                              "te_annotation_basis": basis,
                              "annotation_group": "te_related" if te else "other" if informative else "uninformative"}
    return result


def read_raw_hits(path, queries):
    rows, targets = [], set()
    with path.open() as handle:
        for line in handle:
            f = line.rstrip().split("\t")
            if len(f) != 14 or f[0] not in queries:
                raise ValueError("Invalid Swiss-Prot alignment row/query")
            numeric = list(map(float, f[2:6]))
            aln, qlen, tlen, qs, qe, ts, te = int(f[3]), *map(int, f[6:12])
            if (any(not math.isfinite(v) for v in numeric) or not 0 <= numeric[0] <= 100
                    or numeric[2] < 0 or numeric[3] < 0 or qlen != len(queries[f[0]])
                    or not 1 <= qs <= qe <= qlen or not 1 <= ts <= te <= tlen
                    or aln != len(f[12]) or aln != len(f[13])
                    or f[12].replace("-", "") != queries[f[0]][qs - 1:qe]
                    or len(f[13].replace("-", "")) != te - ts + 1):
                raise ValueError("Invalid Swiss-Prot alignment measurements")
            paired = sum(q != "-" and t != "-" for q, t in zip(f[12], f[13], strict=True))
            rows.append({"query": f[0], "target": f[1], "identity_pct": numeric[0], "alignment_length": aln,
                         "evalue": numeric[2], "bits": numeric[3], "query_length": qlen,
                         "target_length": tlen, "query_start": qs, "query_end": qe,
                         "target_start": ts, "target_end": te,
                         "paired_residues": paired,
                         "query_coverage": paired / qlen, "target_coverage": paired / tlen,
                         "target_aligned_sequence": f[13].replace("-", "")})
            targets.add(f[1])
    result = {q: [] for q in queries}
    for row in rows:
        result[row["query"]].append(row)
    return result


def annotate_hits(raw, annotations):
    seen = set()
    result = {q: [] for q in raw}
    for original in (h for rows in raw.values() for h in rows):
        row = dict(original)
        annotation = annotations[accession(row["target"])]
        if (row["target_length"] != annotation["target_length"]
                or row.pop("target_aligned_sequence") != annotation["sequence"][row["target_start"] - 1:row["target_end"]]
                or (row["query"], row["target"]) in seen):
            raise ValueError("Duplicate Swiss-Prot hit or target length/DB mismatch")
        seen.add((row["query"], row["target"]))
        result[row["query"]].append({**row, **{k: v for k, v in annotation.items() if k != "sequence"}})
    return {q: sorted(hits, key=lambda h: (-h["bits"], h["evalue"], h["accession"])) for q, hits in result.items()}


def read_hits(path, queries, prefix, metadata):
    raw = read_raw_hits(path, queries)
    targets = {h["target"] for rows in raw.values() for h in rows}
    annotations = reference_annotations(prefix, metadata, targets) if targets else {}
    return annotate_hits(raw, annotations)


def diagnostics(hits, parameters):
    """Disjoint failure reasons and an independent, unranked partial-TE flag."""
    def length(h):
        return h.get("paired_residues", h["alignment_length"])
    significant = [h for h in hits if h["evalue"] <= parameters["evalue"]]
    covered = [h for h in significant if h["query_coverage"] >= parameters["query_coverage"]
               and h["target_coverage"] >= parameters["target_coverage"]]
    complete = [h for h in covered if length(h) >= parameters["minimum_alignment"]]
    category, _ = classify(hits, parameters)
    reason = ""
    if category == "no_informative_hit":
        reason = ("no_returned_hits" if not hits else "weak_hit" if not significant
                  else "short_hit" if covered and not complete
                  else "partial_hit" if not covered else "annotation_unknown")
    partial = sorted({h["accession"] for h in significant if h["annotation_group"] == "te_related"
                      and length(h) >= parameters["minimum_alignment"]
                      and h["query_coverage"] >= parameters["query_coverage"]
                      and h["target_coverage"] < parameters["target_coverage"]})
    return {"no_support_reason": reason, "partial_te_accessions": partial,
            "partial_te_homology": bool(partial)}


def cached_annotations(prefix, metadata, targets, cache, boundary):
    signature = {"files": boundary.read([Path(str(prefix) + ".pep"), metadata]),
                 "implementation": encoded_hash([inspect.getsource(reference_annotations),
                                                  TE_PATTERN.pattern, UNINFORMATIVE.pattern])}
    directory = cache / "annotations" / encoded_hash(signature)
    directory.mkdir(parents=True, exist_ok=True)
    wanted = {accession(t) for t in targets}
    result = {}
    with exclusive_lock(directory / "cache.lock"), sqlite3.connect(directory / "annotations.sqlite3") as connection:
        connection.execute("CREATE TABLE IF NOT EXISTS annotations (accession TEXT PRIMARY KEY, payload TEXT, sha256 TEXT)")
        # One read cursor holds the database lock once. Per-accession SELECTs
        # otherwise cause thousands of network-filesystem lock round trips.
        connection.execute("PRAGMA temp_store=MEMORY")
        connection.execute("CREATE TEMP TABLE wanted (accession TEXT PRIMARY KEY)")
        connection.executemany("INSERT INTO wanted VALUES (?)", ((a,) for a in sorted(wanted)))
        for identifier, payload, sha in connection.execute(
                "SELECT a.accession,a.payload,a.sha256 FROM wanted w CROSS JOIN annotations a ON a.accession=w.accession"):
            if hashlib.sha256(payload.encode()).hexdigest() == sha:
                value = json.loads(payload)
                if isinstance(value, dict) and value.get("accession") == identifier:
                    result[identifier] = value
        missing = wanted - set(result)
        if missing:
            fresh = reference_annotations(prefix, metadata, missing)
            boundary.check()
            for identifier, value in fresh.items():
                payload = json.dumps(value, sort_keys=True)
                connection.execute("INSERT OR REPLACE INTO annotations VALUES (?,?,?)",
                                   (identifier, payload, hashlib.sha256(payload.encode()).hexdigest()))
            result.update(fresh)
        boundary.check()
    return result, signature


def classify(hits, parameters):
    qualifying = [h for h in hits if h["evalue"] <= parameters["evalue"]
                  and h["query_coverage"] >= parameters["query_coverage"]
                  and h["target_coverage"] >= parameters["target_coverage"]
                  and h.get("paired_residues", h["alignment_length"]) >= parameters["minimum_alignment"]]
    best = max((h["bits"] for h in qualifying), default=0)
    supported = [h for h in qualifying if h["bits"] >= best * parameters["score_fraction"]]
    groups = {h["annotation_group"] for h in supported}
    category = ("both" if {"te_related", "other"} <= groups else "te_only" if "te_related" in groups
                else "other_only" if "other" in groups else "no_informative_hit")
    return category, [h["accession"] for h in supported if h["annotation_group"] != "uninformative"]


def search(queries, prefix, metadata, cache, parameters, tool, boundary, cpus, memory_gb, scratch):
    db = Path(str(prefix) + ".mmseqs")
    paths = sorted(p for p in db.parent.glob(db.name + "*") if p.is_file() and not p.name.endswith(".ready"))
    if not db.is_file() or not Path(str(db) + ".dbtype").is_file():
        raise ValueError("Swiss-Prot MMseqs2 database is not ready")
    signature = {"schema": 2, "files": boundary.read([*paths, Path(str(prefix) + ".pep")]),
                 "parameters": {k: parameters[k] for k in ("search_evalue", "max_hits", "sensitivity")},
                 "tool": tool, "implementation": encoded_hash(inspect.getsource(read_raw_hits)),
                 "binary": boundary.read([Path(shutil.which("mmseqs"))])}
    key = encoded_hash(signature)
    directory = cache / key
    directory.mkdir(parents=True, exist_ok=True)
    with exclusive_lock(directory / "cache.lock"):
        with sqlite3.connect(directory / "queries.sqlite3") as connection:
            connection.execute("CREATE TABLE IF NOT EXISTS hits (query TEXT PRIMARY KEY, payload TEXT, sha256 TEXT)")
            connection.execute("PRAGMA temp_store=MEMORY")
            connection.execute("CREATE TEMP TABLE wanted (query TEXT PRIMARY KEY)")
            connection.executemany("INSERT INTO wanted VALUES (?)", ((q,) for q in sorted(queries)))
            result = {}
            for query, payload, sha in connection.execute(
                    "SELECT h.query,h.payload,h.sha256 FROM wanted w CROSS JOIN hits h ON h.query=w.query"):
                if hashlib.sha256(payload.encode()).hexdigest() == sha:
                    value = json.loads(payload)
                    if isinstance(value, list) and all(h.get("query") == query for h in value):
                        result[query] = value
            missing = {q: pep for q, pep in queries.items() if q not in result}
            if missing:
                with tempfile.TemporaryDirectory(prefix="gg-swissprot-", dir=scratch) as temporary:
                    tmp = Path(temporary)
                    (tmp / "query.fa").write_text("".join(f">{q}\n{pep}\n" for q, pep in sorted(missing.items())))
                    commands = [
                        ["mmseqs", "createdb", str(tmp / "query.fa"), str(tmp / "queryDB")],
                        ["mmseqs", "search", str(tmp / "queryDB"), str(db), str(tmp / "resultDB"), str(tmp / "tmp"),
                         "--threads", str(cpus), "--split-memory-limit", f"{memory_gb}G", "--max-seqs", str(parameters["max_hits"]),
                         "-e", str(parameters["search_evalue"]), "-s", str(parameters["sensitivity"]), "-a", "1"],
                        ["mmseqs", "convertalis", str(tmp / "queryDB"), str(db), str(tmp / "resultDB"), str(tmp / "hits.tsv"),
                         "--threads", str(cpus), "--format-output", FIELDS]]
                    for command in commands:
                        subprocess.run(command, check=True)
                    fresh = read_raw_hits(tmp / "hits.tsv", missing)
                    boundary.check()
                    for query, hits in fresh.items():
                        payload = json.dumps(hits, sort_keys=True)
                        connection.execute("INSERT OR REPLACE INTO hits VALUES (?,?,?)",
                                           (query, payload, hashlib.sha256(payload.encode()).hexdigest()))
                    result.update(fresh)
            boundary.check()
    targets = {h["target"] for rows in result.values() for h in rows}
    annotations, annotation_signature = cached_annotations(prefix, metadata, targets, cache, boundary)
    return annotate_hits(result, annotations), {"alignment": signature, "annotation": annotation_signature}, len(missing)


def audit(args):
    started = time.monotonic()
    root, output = args.rescue_output.resolve(), args.output.resolve()
    if root == output or root in output.parents or output in root.parents:
        raise ValueError("Swiss-Prot audit must be separate from the frozen rescue directory")
    cache = args.cache.resolve()
    if (cache == output or cache in output.parents or output in cache.parents
            or cache == root or root in cache.parents):
        raise ValueError("Swiss-Prot cache must be separate from the audit and frozen rescue directories")
    boundary = FreshDigestBatch()
    records, bindings, inputs, plan_hash, receipt_hash = candidates(root, boundary)
    parameters = dict(DEFAULTS)
    for name in parameters:
        parameters[name] = getattr(args, name, parameters[name])
    if getattr(args, "search_evalue", None) is None:
        # Preserve callers that only configured the support cutoff before the
        # raw search bound became a separate option. A tighter support cutoff
        # keeps the default raw search cache; a looser one needs broader hits.
        parameters["search_evalue"] = max(DEFAULTS["search_evalue"], parameters["evalue"])
    if (not 0 < parameters["evalue"] <= parameters["search_evalue"] <= 1 or not all(0 < parameters[n] <= 1 for n in ("query_coverage", "target_coverage", "score_fraction"))
            or parameters["minimum_alignment"] < 1 or parameters["max_hits"] < 2
            or not 1 <= parameters["sensitivity"] <= 7.5 or args.cpus < 1 or args.memory_gb < 1
            or any(not math.isfinite(v) for v in parameters.values())):
        raise ValueError("Invalid Swiss-Prot search/support parameters")
    queries = {r["protein_sha256"]: r["protein"] for r in records if r["protein"] is not None}
    version = subprocess.check_output(["mmseqs", "version"], text=True).strip()
    hits, signature, searched = search(queries, args.db_prefix.resolve(), args.metadata.resolve(), args.cache.resolve(),
                                        parameters, version, boundary, args.cpus, args.memory_gb, args.scratch)
    key = {"schema": 1, "rescue_plan_sha256": plan_hash, "augmented_receipt_sha256": receipt_hash,
           "inputs": inputs, "species": bindings, "search": signature, "parameters": parameters,
           "implementation": boundary.read([Path(__file__)]), "candidates_sha256": encoded_hash(records)}
    if any(output == Path(p) or output in Path(p).parents for p in boundary.paths):
        raise ValueError("Audit destination would overwrite an input or database")
    def build(tmp):
        per_gene = defaultdict(list)
        for record in records:
            peptide = record["protein_sha256"]
            category, accessions = classify(hits[peptide], parameters) if peptide else ("not_assessed", [])
            per_gene[(record["species"], record["gene_id"])].append(
                {k: v for k, v in record.items() if k != "protein"}
                | {"category": category, "support_accessions": accessions, "hits": hits.get(peptide, []),
                   **(diagnostics(hits[peptide], parameters) if peptide else
                      {"no_support_reason": "", "partial_te_accessions": [], "partial_te_homology": False}),
                   "reason": "translation_uncertain" if peptide is None else ""})
        species = {name: {"counts": dict.fromkeys(GROUPS, 0), "loci": {},
                          "no_support_reasons": dict.fromkeys(NO_SUPPORT_REASONS, 0),
                          "partial_te_homology_loci": 0} for name in bindings}
        for (name, gene), coding in sorted(per_gene.items()):
            groups = {r["category"] for r in coding}
            te, other = bool(groups & {"te_only", "both"}), bool(groups & {"other_only", "both"})
            category = ("both" if te and other else "te_only" if te else "other_only" if other
                        else "not_assessed" if "not_assessed" in groups else "no_informative_hit")
            species[name]["counts"][category] += 1
            flags = sorted({a for c in coding for a in c["partial_te_accessions"]})
            reasons = {c["no_support_reason"] for c in coding} - {""}
            # The most informative available coding sequence determines a
            # single reason; no-hit isoforms cannot hide a partial/short hit.
            reason = next((r for r in reversed(NO_SUPPORT_REASONS) if r in reasons), "") if category == "no_informative_hit" else ""
            if reason:
                species[name]["no_support_reasons"][reason] += 1
            species[name]["partial_te_homology_loci"] += bool(flags)
            species[name]["loci"][gene] = {"category": category, "coding_sequences": coding,
                                          "no_support_reason": reason, "partial_te_accessions": flags,
                                          "partial_te_homology": bool(flags)}
        atomic_json(tmp / "evidence.json", {"schema": 1, "method": METHOD, "parameters": parameters, "species": species})
        with (tmp / "loci.tsv").open("w") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["species", "gene_id", "category", "cds_ids", "support_accessions",
                             "partial_te_homology", "partial_te_accessions", "no_support_reason"])
            for name, data in sorted(species.items()):
                for gene, record in sorted(data["loci"].items()):
                    writer.writerow([name, gene, record["category"], ";".join(c["cds_id"] for c in record["coding_sequences"]),
                                     ";".join(sorted({a for c in record["coding_sequences"] for a in c["support_accessions"]})),
                                     record["partial_te_homology"], ";".join(record["partial_te_accessions"]),
                                     record["no_support_reason"]])
        with (tmp / "hits.tsv").open("w") as handle:
            columns = ["species", "gene_id", "cds_id", "category", "query", "accession", "protein_name", "annotation_group",
                       "evalue", "bits", "identity_pct", "alignment_length", "paired_residues", "query_coverage", "target_coverage",
                       "query_start", "query_end", "target_start", "target_end", "keywords", "go_terms", "te_annotation_basis", "organism"]
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for name, data in sorted(species.items()):
                for gene, record in sorted(data["loci"].items()):
                    for coding in record["coding_sequences"]:
                        for hit in coding["hits"]:
                            writer.writerow({**hit, "species": name, "gene_id": gene, "cds_id": coding["cds_id"], "category": coding["category"],
                                             "te_annotation_basis": ";".join(hit["te_annotation_basis"])})
        atomic_json(tmp / "summary.json", {"schema": 1, "species": {n: v["counts"] for n, v in species.items()},
                                           "coding_sequences": len(records), "unique_proteins": len(queries),
                                           "method": METHOD, "sequence_changes": 0})
    destination = stage(output.parent, output.name, key, build, guard=boundary.check)
    atomic_json(output.parent / (output.name + ".execution.json"), {"elapsed_seconds": time.monotonic() - started,
                "searched_unique_proteins": searched, "cached_unique_proteins": len(queries) - searched,
                "coding_sequences": len(records), "cpus": args.cpus, "memory_gb": args.memory_gb,
                "receipt_sha256": digest(destination / "receipt.json")})
    return destination


def collect(changes, directory):
    """Bind plotted categories to exactly the rescued loci in Before."""
    if directory is None:
        return changes
    directory = Path(directory)
    receipt, receipt_hash = read_json_snapshot(directory / "receipt.json")
    evidence, evidence_hash = read_json_snapshot(directory / "evidence.json")
    key = receipt["key"]
    selection = changes.get("rescue_reference_selection", {})
    if (key.get("schema") != 1 or key.get("rescue_plan_sha256") != selection.get("plan_sha256")
            or key.get("augmented_receipt_sha256") != selection.get("augmented_receipt_sha256")
            or receipt["files"].get("evidence.json") != evidence_hash or evidence.get("schema") != 1
            or evidence.get("parameters") != key.get("parameters")
            or set(key["species"]) != set(evidence["species"])):
        raise ValueError("Swiss-Prot audit belongs to different rescue inputs")
    updates, diagnostic_updates = {}, {}
    for name, value in changes["species"].items():
        if value["refinement_status"] == "not_analysed":
            if name in evidence["species"]:
                raise ValueError("Swiss-Prot audit includes an unanalysed species")
            updates[name] = None
            continue
        loci = changes.get("evidence", {}).get(name, {}).get("rescued_loci_support")
        data = evidence["species"].get(name)
        if (not isinstance(loci, dict) or data is None or set(loci) != set(data["loci"])
                or set(loci) != set(key["species"][name]["gene_ids"])
                or key["species"][name]["source_gff_sha256"] != changes["evidence"][name]["source_gff_sha256"]):
            raise ValueError("Swiss-Prot audit rescued-locus identities/source differ: " + name)
        counts = Counter(r["category"] for r in data["loci"].values())
        expected = {g: counts[g] for g in GROUPS}
        if set(counts) - set(GROUPS) or expected != data["counts"] or sum(counts.values()) != value["prior_rescued_loci"]:
            raise ValueError("Swiss-Prot categories do not sum to rescued loci")
        updates[name] = expected
        diagnostic_updates[name] = validate_diagnostics(data, expected)
    if set(evidence["species"]) != {n for n, v in changes["species"].items() if v["refinement_status"] != "not_analysed"}:
        raise ValueError("Swiss-Prot audit species membership differs")
    if digest(directory / "receipt.json") != receipt_hash or digest(directory / "evidence.json") != evidence_hash:
        raise ValueError("Swiss-Prot audit changed while loading")
    for name, counts in updates.items():
        changes["species"][name]["rescue_swissprot_groups"] = counts
        for field in ("rescue_partial_te_groups", "rescue_no_support_reasons"):
            changes["species"][name][field] = diagnostic_updates.get(name, {}).get(field)
    changes["swissprot_evidence"] = {"directory": str(directory.resolve()), "receipt_sha256": receipt_hash,
                                     "evidence_sha256": evidence_hash, "method": evidence["method"], "parameters": evidence["parameters"]}
    return changes


def validate_diagnostics(data, counts):
    if "no_support_reasons" not in data:
        return {}
    reasons = Counter(r.get("no_support_reason", "") for r in data["loci"].values()
                      if r["category"] == "no_informative_hit")
    expected = {r: reasons[r] for r in NO_SUPPORT_REASONS}
    if (set(reasons) - set(NO_SUPPORT_REASONS) or expected != data["no_support_reasons"]
            or sum(expected.values()) != counts["no_informative_hit"]
            or any(type(r.get("partial_te_homology")) is not bool
                   or not isinstance(r.get("partial_te_accessions"), list)
                   or r["partial_te_homology"] != bool(r["partial_te_accessions"])
                   for r in data["loci"].values())
            or sum(r["partial_te_homology"] for r in data["loci"].values()) != data["partial_te_homology_loci"]):
        raise ValueError("Swiss-Prot diagnostic counts differ from rescued loci")
    groups = Counter("primary_te_support" if r["category"] in {"te_only", "both"} else
                     "not_assessed" if r["category"] == "not_assessed" else
                     "partial_te_only" if r["partial_te_homology"] else "no_te_support"
                     for r in data["loci"].values())
    return {"rescue_partial_te_groups": {g: groups[g] for g in
                                       ("primary_te_support", "partial_te_only", "no_te_support", "not_assessed")},
            "rescue_no_support_reasons": expected}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("rescue-output", "output", "db-prefix", "metadata", "cache"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--memory-gb", type=int, default=8)
    parser.add_argument("--scratch", type=Path, help="Local temporary storage; final evidence remains in --output")
    for name, value in DEFAULTS.items():
        parser.add_argument("--" + name.replace("_", "-"), type=type(value),
                            default=None if name == "search_evalue" else value)
    audit(parser.parse_args())


if __name__ == "__main__":
    main()
