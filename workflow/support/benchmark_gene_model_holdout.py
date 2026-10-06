#!/usr/bin/env python3
"""Controlled coding-path holdouts from real reference CDS/GFF/genome triples.

This evaluates interval prediction and admission with exact-source proteins.
It does not estimate cohort-wide discovery precision or independent homology
support. Reference annotations are the benchmark truth, not experimental proof.
"""
import argparse
import csv
import hashlib
import inspect
import json
import re
import shutil
from collections import defaultdict
from pathlib import Path

from Bio.Seq import Seq
from fasta_sequence_store import fasta_records, open_text
from gene_model_catalog import indexed_genome, reconstruct_sequence, validate_candidate
from gene_model_species_profiles import parameters_for, read_profiles
from input_generation_array_state import FreshDigestBatch, atomic_json, digest
from rescue_gene_models import attributes, consolidate, search_intervals, validate_model

DEFAULTS = {"minimum_identity": .5, "minimum_coverage": .95, "max_intron": 20000,
            "max_interval": 200000, "padding": 300}


def references(gff, cds, genome, code):
    """Associate protein IDs only after exact CDS/assembly sequence agreement."""
    sequences = {}
    for _identifier, header, sequence in fasta_records(cds):
        match = re.search(r"\[protein_id=([^\]]+)\]", header)
        if match and match[1] not in sequences:
            sequences[match[1]] = sequence.upper()
    blocks, metadata = defaultdict(list), {}
    with open_text(gff) as handle:
        for line in handle:
            if line.strip() == "##FASTA":
                break
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attr = attributes(fields[8])
            protein = attr.get("protein_id")
            if (protein not in sequences or attr.get("transl_table", "1") != str(code)
                    or attr.get("exception") or attr.get("transl_except") or attr.get("pseudo") == "true"):
                continue
            blocks[protein].append([int(fields[3]) - 1, int(fields[4]), int(fields[7])])
            metadata[protein] = fields[0], fields[6], attr.get("Parent", protein)
    candidates = []
    for protein, rows in sorted(blocks.items()):
        seqid, strand, transcript = metadata[protein]
        rows.sort(reverse=strand == "-")
        if seqid not in genome.references:
            continue
        sequence = reconstruct_sequence(rows, seqid, strand, genome)
        if sequence != sequences[protein]:
            continue
        candidate = {"seqid": seqid, "strand": strand, "blocks": rows, "cds": sequence,
                     "origin": "original"}
        quality = validate_candidate(candidate, code)
        if not quality["valid_orf"] or not 100 <= len(sequence) // 3 <= 600:
            continue
        candidates.append({**candidate, "protein_id": protein, "transcript": transcript,
                           "protein": str(Seq(sequence).translate(table=code)).removesuffix("*")})
    return candidates


def benchmark(args):
    args.output.mkdir(parents=True, exist_ok=True)
    with args.inputs.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    names = [r["species"] for r in rows]
    if not rows or len(names) != len(set(names)):
        raise ValueError("Benchmark needs unique species inputs")
    profiles = read_profiles(args.species_profiles, names)
    boundary = FreshDigestBatch()
    paths = [args.inputs, Path(__file__)] + [Path(r[k]) for r in rows for k in ("cds", "gff", "genome")]
    paths += [Path(inspect.getmodule(f).__file__) for f in
              (indexed_genome, parameters_for, validate_model, atomic_json, fasta_records)]
    paths.append(Path(shutil.which("miniprot")))
    if args.species_profiles:
        paths.append(args.species_profiles)
    identities = boundary.read(paths)
    result = {}
    for row in rows:
        name, code = row["species"], int(row.get("genetic_code") or 1)
        params = parameters_for({"parameters": DEFAULTS, "species_profiles": profiles}, name)
        with indexed_genome(Path(row["genome"])) as genome:
            truth = references(Path(row["gff"]), Path(row["cds"]), genome, code)
            # Deterministic split: equal single/multiple-exon strata where
            # available, without selecting examples by predictor success.
            strata = [[c for c in truth if (len(c["blocks"]) > 1) == multiple] for multiple in (False, True)]
            selected = []
            for group in strata:
                group.sort(key=lambda c: hashlib.sha256((str(args.seed) + c["protein_id"]).encode()).hexdigest())
                selected.extend(group[:args.loci // 2])
            if len(selected) < args.loci:
                chosen = {c["protein_id"] for c in selected}
                selected.extend(c for c in truth if c["protein_id"] not in chosen)
                selected = selected[:args.loci]
            if len(selected) != args.loci:
                raise ValueError("Insufficient eligible truth loci: " + name)
            details = []
            for i, candidate in enumerate(selected):
                start = max(0, min(b[0] for b in candidate["blocks"]) - params["padding"])
                end = min(genome.get_reference_length(candidate["seqid"]),
                          max(b[1] for b in candidate["blocks"]) + params["padding"])
                if end - start > params["max_interval"]:
                    details.append({"protein_id": candidate["protein_id"], "status": "interval_bound", "accepted": 0})
                    continue
                window = genome.fetch(candidate["seqid"], start, end)
                for condition in ("held_out", "assembly_gap_control"):
                    directory = args.output / name / f"{i:03d}" / condition
                    directory.mkdir(parents=True, exist_ok=True)
                    sequence = window
                    if condition == "assembly_gap_control":
                        # Remove most of every coding exon. No intact truth
                        # coding path remains in this deliberately damaged locus.
                        chars = list(window)
                        for a, b, _phase in candidate["blocks"]:
                            for position in range(a + 3, b - 3):
                                chars[position - start] = "N"
                        sequence = "".join(chars)
                    fasta = directory / "interval.fa"
                    fasta.write_text(">interval\n" + sequence + "\n")
                    region = {"id": "q", "seqid": "interval", "start": 0, "end": len(sequence),
                              "query": candidate["protein_id"], "donor": name}
                    with indexed_genome(fasta) as local:
                        models = search_intervals(directory, {("interval", 0, len(sequence)): [region]},
                                                  {name: {candidate["protein_id"]: candidate["protein"]}},
                                                  local, code, params["max_intron"], args.cpus)
                        checked = []
                        for model in models:
                            model["seqid"] = "interval"
                            model["evidence"] = region
                            checked.append(validate_model(model, local, code, params))
                        published = consolidate(checked, [], name)
                    accepted = [m for m in published if m["status"] == "accepted"]
                    expected = [[a - start, b - start, p] for a, b, p in candidate["blocks"]]
                    exact = sum(m["strand"] == candidate["strand"] and m["cds"] == expected
                                and m["sequence"] == candidate["cds"] for m in accepted)
                    details.append({"protein_id": candidate["protein_id"], "condition": condition,
                                    "exons": len(expected), "strand": candidate["strand"],
                                    "accepted": len(accepted), "exact_truth_paths": exact,
                                    "false_or_inexact_paths": len(accepted) - exact,
                                    "problems": sorted({p for m in published for p in m["problems"]})})
            positive = [r for r in details if r.get("condition") == "held_out"]
            negative = [r for r in details if r.get("condition") == "assembly_gap_control"]
            result[name] = {"eligible_truth_loci": len(truth), "held_out_loci": args.loci,
                            "exact_recovered": sum(r["exact_truth_paths"] > 0 for r in positive),
                            "missed_truth_loci": args.loci - sum(r["exact_truth_paths"] > 0 for r in positive),
                            "inexact_accepted_paths": sum(r["false_or_inexact_paths"] for r in positive),
                            "gap_control_false_additions": sum(r["accepted"] for r in negative),
                            "parameters": params, "details": details}
        print(json.dumps({name: {k: v for k, v in result[name].items() if k != "details"}}), flush=True)
    boundary.check()
    atomic_json(args.output / "benchmark.json", {"schema": 1, "inputs": identities, "seed": args.seed,
                "species": result, "scope": __doc__, "miniprot_sha256": digest(shutil.which("miniprot"))})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True, help="species,cds,gff,genome,genetic_code TSV")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--species-profiles", type=Path)
    parser.add_argument("--loci", type=int, default=20)
    parser.add_argument("--seed", type=int, default=169)
    parser.add_argument("--cpus", type=int, default=2)
    args = parser.parse_args()
    args.inputs = args.inputs.resolve()
    args.output = args.output.resolve()
    if args.species_profiles:
        args.species_profiles = args.species_profiles.resolve()
    if args.loci < 2 or args.cpus < 1:
        parser.error("At least two loci and positive CPUs are required")
    benchmark(args)


if __name__ == "__main__":
    main()
