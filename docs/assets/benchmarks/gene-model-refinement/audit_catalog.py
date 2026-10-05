#!/usr/bin/env python3
"""Read-only admission audit of the existing seven-species curated inputs.

Run through workflow/tests/run_in_runtime.sh from the repository root. No
annotation/CDS/genome is regenerated, and scratch FASTA indices are temporary.
The elapsed time is an observation from a non-isolated audit, not a benchmark.
"""

import argparse
import hashlib
import importlib.metadata
import json
import platform
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def signature(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--runtime-image", required=True)
    args = parser.parse_args()
    support = Path.cwd() / "workflow" / "support"
    sys.path.insert(0, str(support))
    from gene_model_catalog import build_catalog

    species = ("Amborella_trichopoda", "Arabidopsis_thaliana", "Cephalotus_follicularis",
               "Dionaea_muscipula", "Nepenthes_gracilis", "Oryza_sativa", "Spinacia_oleracea")
    before, rows, aggregate = {}, [], Counter()
    started = time.perf_counter()
    for name in species:
        paths = {}
        for role in ("cds", "gff", "genome"):
            choices = sorted((args.input / ("species_" + role)).glob(name + "*"))
            if len(choices) != 1:
                raise ValueError("Expected one existing curated source for " + name + ":" + role)
            paths[role] = choices[0].resolve()
            before[str(paths[role])] = signature(paths[role])
        began = time.perf_counter()
        catalog = build_catalog(name, paths["cds"], paths["gff"], paths["genome"], genetic_code=1)
        candidates = [candidate for locus in catalog["loci"] for candidate in locus["candidates"]]
        mapping = Counter(row["mapping_status"] for row in catalog["fasta_mapping"])
        conventions = Counter(row.get("source_convention", "") for row in catalog["fasta_mapping"])
        withheld = [candidate for candidate in candidates if not candidate["quality"]["usable"]]
        withheld_reasons = Counter(key for candidate in withheld for key in
                                  ("internal_stop", "ambiguous", "invalid_base", "phase_conflict",
                                   "phase_unresolved", "sequence_mismatch", "structure_problem", "annotated_exception")
                                  if candidate["quality"].get(key))
        disagreements = []
        for candidate in candidates:
            for source in candidate.get("source_cds", []):
                if source["sequence_agreement"]:
                    continue
                supplied, genomic = source["cds"].upper(), candidate["cds"]
                differences = [(index, a, b) for index, (a, b) in enumerate(zip(supplied, genomic, strict=False)) if a != b]
                same_length = len(supplied) == len(genomic)
                masked_only = same_length and all(a == "N" for _, a, _ in differences)
                disagreements.append({"source_fasta_id": source["source_fasta_id"],
                                      "candidate_id": candidate["candidate_id"],
                                      "source_cds_sha256": source["sha256"],
                                      "genomic_cds_sha256": hashlib.sha256(genomic.encode()).hexdigest(),
                                      "source_length": len(supplied), "genomic_length": len(genomic),
                                      "same_length_source_masking_only": masked_only,
                                      "first_ungapped_differences_start0": differences[:8],
                                      "quality": {key: candidate["quality"].get(key) for key in
                                                  ("usable", "internal_stop", "ambiguous", "partial",
                                                   "phase_unresolved", "sequence_mismatch")}})
        counts = {"loci": len(catalog["loci"]), "candidates": len(candidates),
                  "usable_candidates": len(candidates) - len(withheld), "withheld_candidates": len(withheld),
                  "source_baseline_loci": sum(bool(locus["source_baseline_candidate_id"]) for locus in catalog["loci"]),
                  "source_cds_records": sum(mapping.values()), "uniquely_mapped_source_cds": mapping["mapped"],
                  "inferred_phase_candidates": sum(candidate["quality"].get("phase_inferred", False) for candidate in candidates),
                  "masked_terminal_stop_source_records": conventions["masked_terminal_stop"],
                  "sequence_mismatch_source_records": len(disagreements),
                  "same_length_source_masking_only_records": sum(row["same_length_source_masking_only"] for row in disagreements),
                  "source_genomic_length_disagreement_records": sum(row["source_length"] != row["genomic_length"] for row in disagreements),
                  "excluded_organelle_loci": len(catalog.get("excluded_loci", []))}
        aggregate.update(counts)
        rows.append({"species": name, "genetic_code": 1, "counts": counts,
                     "sources": catalog["sources"], "mapping_status_counts": dict(mapping),
                     "source_convention_counts": dict(conventions), "withheld_reason_counts_nonexclusive": dict(withheld_reasons),
                     "sequence_mismatches": disagreements, "observed_audit_seconds": time.perf_counter() - began})
    after = {path: signature(path) for path in before}
    if after != before:
        raise ValueError("Curated source bytes changed during admission audit")
    tracked = [support / "gene_model_catalog.py", support / "cds_model_normalisation.py",
               support / "fasta_sequence_store.py", support / "gff_attribute_syntax.py"]
    tracked += sorted((support / "format_species_annotation").glob("*.py"))
    report = {"schema": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "read-only catalog admission; no predictions, source regeneration, or isolated performance claims",
              "runtime": {"docker_image": args.runtime_image, "python": platform.python_version(),
                          "packages": {package: importlib.metadata.version(package) for package in ("biopython", "pysam")}},
              "script_sha256": signature(__file__),
              "implementation_sha256": {str(path.relative_to(support)): signature(path) for path in tracked},
              "aggregate": dict(aggregate), "source_files_checked": len(before),
              "all_source_bytes_unchanged": before == after,
              "source_sha256_before": before, "source_sha256_after": after,
              "species": rows, "observed_total_audit_seconds": time.perf_counter() - started}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"report": str(args.output.resolve()), "aggregate": report["aggregate"],
                      "source_files_checked": len(before), "all_source_bytes_unchanged": before == after}))


if __name__ == "__main__":
    main()
