"""Candidate-only Swiss-Prot audit, real search, cache and immutable joins."""
import copy
import csv
import json
import random
import sqlite3
import subprocess
from argparse import Namespace
from pathlib import Path

import pytest
from Bio.Data import CodonTable

from workflow.support import rescue_swissprot_evidence as swiss
from workflow.support.input_generation_array_state import FreshDigestBatch, atomic_json, digest


def fixture(tmp_path):
    rng = random.Random(171)
    peptides = ["M" + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(99)) for _ in range(5)]
    codons = {aa: codon for codon, aa in CodonTable.unambiguous_dna_by_id[1].forward_table.items()}
    root = tmp_path / "rescue"
    (root / "augmented/species_cds").mkdir(parents=True)
    (root / "augmented/species_gff").mkdir()
    sources, files, rows = {}, {}, []
    for name in ("Animal_one", "Fungus_two"):
        sources[name] = {"genetic_code": 1}
        cds = root / "augmented/species_cds" / (name + ".fa")
        gff = root / "augmented/species_gff" / (name + ".gff3")
        cds.write_text(">original\nATGAAATAA\n" + "".join(
            f">{name}_g{i}\n" + "".join(codons[aa] for aa in p) + "TAA\n" for i, p in enumerate(peptides)))
        gff.write_text("##gff-version 3\n" + "".join(
            f"chr\tgenegalleon_rescue\tgene\t{i*400+1}\t{i*400+303}\t.\t+\t.\tID={name}_g{i}\n"
            f"chr\tgenegalleon_rescue\tmRNA\t{i*400+1}\t{i*400+303}\t.\t+\t.\tID={name}_g{i}.t1;Parent={name}_g{i}\n"
            for i in range(len(peptides))))
        for path in (cds, gff):
            files[str(path.relative_to(root / "augmented"))] = digest(path)
        rows.append({"species": name, "rescued_models": len(peptides), "cds": str(cds), "gff": str(gff)})
    with (root / "augmented/inputs.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    files["inputs.tsv"] = digest(root / "augmented/inputs.tsv")
    atomic_json(root / "plan.json", {"request": {"sources": sources}})
    atomic_json(root / "augmented/receipt.json", {"key": {"plan": digest(root / "plan.json"),
                                                           "rescue_receipts": dict.fromkeys(sources, "frozen")}, "files": files})
    prefix = tmp_path / "sprot"
    reference = [("TE1", "Mariner transposase", peptides[0], "Transposable element", "DNA transposition"),
                 ("HOST1", "DNA helicase", peptides[1], "DNA-binding", "DNA helicase activity"),
                 ("TE2", "Retrotransposon protein", peptides[2], "", ""),
                 ("HOST2", "ATPase", peptides[2], "", "ATPase activity"),
                 ("UNKNOWN", "Uncharacterized protein", peptides[4], "", "")]
    Path(str(prefix) + ".pep").write_text("".join(f">sp|{a}|{a}_REF {n} OS=Test organism\n{p}\n"
                                                  for a, n, p, _, _ in reference))
    meta = Path(str(prefix) + ".meta.tsv")
    with meta.open("w") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["accession", "keywords", "go_terms"])
        writer.writerows((a, k, go) for a, _, _, k, go in reference)
    args = Namespace(rescue_output=root, output=tmp_path / "audit", db_prefix=prefix, metadata=meta,
                     cache=tmp_path / "cache", scratch=tmp_path, cpus=2, memory_gb=2, **swiss.DEFAULTS)
    return args, peptides


def hit(group="te_related", bits=100, **values):
    return {"annotation_group": group, "accession": group, "evalue": 1e-20, "bits": bits,
            "query_coverage": .8, "target_coverage": .8, "alignment_length": 100, **values}


@pytest.mark.parametrize(("hits", "category"), [
    ([hit()], "te_only"), ([hit("other")], "other_only"),
    ([hit(), hit("other", 95)], "both"), ([hit(), hit("other", 70)], "te_only"),
    ([hit(query_coverage=.1)], "no_informative_hit"),
    ([hit(target_coverage=.1)], "no_informative_hit"),
    ([hit(alignment_length=20)], "no_informative_hit"),
    ([hit("uninformative")], "no_informative_hit"), ([], "no_informative_hit"),
])
def test_support_requires_coverage_and_competitive_scores(hits, category):
    assert swiss.classify(hits, swiss.DEFAULTS)[0] == category


@pytest.mark.parametrize(("hits", "reason"), [
    ([], "no_returned_hits"), ([hit(evalue=.1)], "weak_hit"),
    ([hit(alignment_length=20)], "short_hit"),
    ([hit(target_coverage=.1)], "partial_hit"),
    ([hit("uninformative")], "annotation_unknown"),
])
def test_missing_support_reasons_are_disjoint(hits, reason):
    assert swiss.diagnostics(hits, swiss.DEFAULTS)["no_support_reason"] == reason


def test_partial_te_flag_is_independent_of_primary_host_and_competing_scores():
    hits = [hit("other", 1000), hit(bits=100, target_coverage=.1)]
    assert swiss.classify(hits, swiss.DEFAULTS)[0] == "other_only"
    assert swiss.diagnostics(hits, swiss.DEFAULTS)["partial_te_accessions"] == ["te_related"]


@pytest.mark.parametrize("implicit", ["absent", "none"])
def test_legacy_loose_support_cutoff_expands_implicit_search_bound(tmp_path, implicit):
    args, _ = fixture(tmp_path)
    subprocess.run(["mmseqs", "createdb", str(args.db_prefix) + ".pep", str(args.db_prefix) + ".mmseqs"], check=True)
    args.evalue, args.search_evalue = 1e-3, None
    if implicit == "absent":
        del args.search_evalue
    swiss.audit(args)
    evidence = json.loads((args.output / "evidence.json").read_text())
    assert evidence["parameters"]["evalue"] == evidence["parameters"]["search_evalue"] == 1e-3
    # An explicitly chosen insufficient raw bound must still fail, rather than
    # silently making support claims from a truncated search.
    args.search_evalue = 1e-5
    with pytest.raises(ValueError, match="search/support parameters"):
        swiss.audit(args)


def test_support_threshold_and_metadata_changes_reuse_raw_search(tmp_path):
    args, _ = fixture(tmp_path)
    subprocess.run(["mmseqs", "createdb", str(args.db_prefix) + ".pep", str(args.db_prefix) + ".mmseqs"], check=True)
    swiss.audit(args)
    args.minimum_alignment = 150
    swiss.audit(args)
    assert json.loads((args.output.parent / "audit.execution.json").read_text())["searched_unique_proteins"] == 0
    assert all(data["counts"]["no_informative_hit"] == 5
               for data in json.loads((args.output / "evidence.json").read_text())["species"].values())
    args.minimum_alignment = 50
    args.metadata.write_text(args.metadata.read_text().replace("HOST1\tDNA-binding", "HOST1\tTransposable element"))
    swiss.audit(args)
    assert json.loads((args.output.parent / "audit.execution.json").read_text())["searched_unique_proteins"] == 0
    assert all(data["counts"]["te_only"] == 2
               for data in json.loads((args.output / "evidence.json").read_text())["species"].values())


def test_candidates_ignore_original_genes_deduplicate_only_for_search_and_bind_inputs(tmp_path):
    args, peptides = fixture(tmp_path)
    records, bindings, _, _, _ = swiss.candidates(args.rescue_output, FreshDigestBatch())
    assert len(records) == 10 and len({r["protein_sha256"] for r in records}) == 5
    assert {r["protein"] for r in records} == set(peptides)
    assert all(len(b["gene_ids"]) == 5 for b in bindings.values())
    cds = next((args.rescue_output / "augmented/species_cds").iterdir())
    cds.write_text(cds.read_text().replace("ATGAAATAA", "ATGCCCTAA"))
    with pytest.raises(ValueError, match="publication changed"):
        swiss.candidates(args.rescue_output, FreshDigestBatch())


def test_translation_uses_species_code_and_preserves_uncertain_context():
    assert swiss.protein("ATGTGATAA", 4) == "MW"
    assert swiss.protein("TTGAAATAA", 1) == "LK"
    assert swiss.protein("ATGAAATGA", 27) is None
    with pytest.raises(ValueError, match="internal stop"):
        swiss.protein("ATGTGATAA", 1)


def test_metadata_distinguishes_explicit_te_annotations_and_uncharacterized_protein(tmp_path):
    args, _ = fixture(tmp_path)
    annotations = swiss.reference_annotations(args.db_prefix, args.metadata, {"sp|TE1|TE1_REF", "UNKNOWN"})
    assert annotations["TE1"]["annotation_group"] == "te_related"
    assert annotations["UNKNOWN"]["annotation_group"] == "uninformative"
    args.metadata.write_text("accession\tkeywords\tgo_terms\n")
    with pytest.raises(ValueError, match="incomplete"):
        swiss.reference_annotations(args.db_prefix, args.metadata, {"TE1"})


@pytest.mark.parametrize("go_term", ["transposable element silencing by siRNA-mediated DNA methylation",
                                    "transposable element silencing by piRNA-mediated heterochromatin formation",
                                    "negative regulation of DNA transposition", "DNA transposition"])
def test_host_te_silencing_and_process_annotations_do_not_mark_te_origin(tmp_path, go_term):
    args, _ = fixture(tmp_path)
    args.metadata.write_text(args.metadata.read_text().replace("DNA-binding\tDNA helicase activity", "Transposition\t" + go_term))
    annotations = swiss.reference_annotations(args.db_prefix, args.metadata, {"HOST1"})
    assert annotations["HOST1"]["annotation_group"] == "other"
    assert annotations["HOST1"]["te_annotation_basis"] == []


def test_gapped_alignment_cannot_inflate_paired_coverage(tmp_path):
    args, peptides = fixture(tmp_path)
    query = peptides[0]
    # Both spans cover the complete proteins, but only the first ten residues pair.
    qaln, taln = query + "-" * 90, query[:10] + "-" * 90 + query[10:]
    path = tmp_path / "hits.tsv"
    path.write_text("\t".join(map(str, ["q", "TE1", 100, 190, 1e-20, 100, 100, 100, 1, 100, 1, 100, qaln, taln])) + "\n")
    rows = swiss.read_hits(path, {"q": query}, args.db_prefix, args.metadata)["q"]
    assert rows[0]["query_coverage"] == rows[0]["target_coverage"] == .1
    assert swiss.classify(rows, swiss.DEFAULTS)[0] == "no_informative_hit"


def test_real_mmseqs_search_cache_and_verified_gene_counts(tmp_path):
    args, _ = fixture(tmp_path)
    subprocess.run(["mmseqs", "createdb", str(args.db_prefix) + ".pep", str(args.db_prefix) + ".mmseqs"], check=True)
    original = {str(p): digest(p) for p in args.rescue_output.rglob("*") if p.is_file()}
    swiss.audit(args)
    evidence = json.loads((args.output / "evidence.json").read_text())
    for data in evidence["species"].values():
        assert data["counts"] == dict(zip(swiss.GROUPS, [1, 1, 1, 2, 0], strict=True))
    assert json.loads((args.output / "summary.json").read_text())["unique_proteins"] == 5
    first = digest(args.output / "evidence.json")
    swiss.audit(args)
    assert digest(args.output / "evidence.json") == first
    execution = json.loads((args.output.parent / "audit.execution.json").read_text())
    assert execution["searched_unique_proteins"] == 0 and execution["cached_unique_proteins"] == 5
    cache_db = next(args.cache.glob("*/queries.sqlite3"))
    with sqlite3.connect(cache_db) as connection:
        connection.execute("UPDATE hits SET sha256='corrupted' WHERE query=(SELECT query FROM hits LIMIT 1)")
    swiss.audit(args)
    execution = json.loads((args.output.parent / "audit.execution.json").read_text())
    assert execution["searched_unique_proteins"] == 1 and execution["cached_unique_proteins"] == 4
    assert digest(args.output / "evidence.json") == first
    assert {str(p): digest(p) for p in args.rescue_output.rglob("*") if p.is_file()} == original
    receipt = json.loads((args.output / "receipt.json").read_text())
    changes = {"rescue_reference_selection": {"plan_sha256": receipt["key"]["rescue_plan_sha256"],
               "augmented_receipt_sha256": receipt["key"]["augmented_receipt_sha256"]}, "species": {}, "evidence": {}}
    for name, data in evidence["species"].items():
        changes["species"][name] = {"refinement_status": "analysed", "prior_rescued_loci": 5}
        changes["evidence"][name] = {"source_gff_sha256": receipt["key"]["species"][name]["source_gff_sha256"],
                                     "rescued_loci_support": dict.fromkeys(data["loci"], {})}
    changes["species"]["Excluded_species"] = {"refinement_status": "not_analysed", "prior_rescued_loci": None}
    swiss.collect(changes, args.output)
    assert changes["species"]["Excluded_species"]["rescue_swissprot_groups"] is None
    broken = copy.deepcopy(changes)
    broken["evidence"]["Animal_one"]["source_gff_sha256"] = "changed"
    with pytest.raises(ValueError, match="identities/source differ"):
        swiss.collect(broken, args.output)
    evidence["species"]["Animal_one"]["counts"]["te_only"] = 99
    atomic_json(args.output / "evidence.json", evidence)
    with pytest.raises(ValueError, match="different rescue inputs"):
        swiss.collect(changes, args.output)


@pytest.mark.parametrize("enabled", [False, True])
def test_refinement_finisher_reuses_annotation_and_forwards_verified_audit(tmp_path, enabled):
    import os
    core = Path(__file__).resolve().parents[1] / "core/gg_input_generation_core.sh"
    function = core.read_text().split("finish_gene_model_refinement() {", 1)[1].split("\n}\n", 1)[0]
    anchor = tmp_path / "anchor"
    (anchor / "augmented").mkdir(parents=True)
    (anchor / "augmented/receipt.json").write_text("{}")
    log = tmp_path / "commands"
    script = '''set -euo pipefail
python() { printf '%s\\n' "$*" >> "$COMMAND_LOG"; }
annotate_gene_model_rescue_swissprot() { printf 'audit %s\\n' "$1" >> "$COMMAND_LOG"; }
ensure_shared_busco_lineage_ready() { busco_lineage_resolved=embryophyta_odb12; }
ensure_busco_download_path() { printf '%s\\n' /db; }
gg_memory_parallel_job_cap() { printf '%s\\n' 2; }
finish_gene_model_refinement() {''' + function + "\n}\nfinish_gene_model_refinement\n"
    env = dict(os.environ, COMMAND_LOG=str(log), run_species_busco="1", gg_support_dir="/support",
               run_gene_model_rescue_swissprot=str(int(enabled)), gene_model_rescue_dir=str(anchor),
               gene_model_refinement_rescue_dir="", gene_model_refinement_inputs="", gene_model_rescue_swissprot_dir="",
               gene_model_refinement_dir="/refinement", species_cds_dir="/all-cds", GG_TASK_CPUS="8",
               species_busco_parallel_jobs="auto", task_plan_output="/plan", gg_workspace_dir="/workspace",
               GG_MEM_TOOL_GB="64", species_busco_memory_gb_per_job="16", busco_lineage_resolved="")
    result = subprocess.run(["bash", "-c", script], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    commands = log.read_text().splitlines()
    assert any(c.startswith("audit ") for c in commands) == enabled
    plot = next(c for c in commands if "gene_model_refinement_busco.py" in c)
    assert ("--rescue-swissprot-dir " + str(anchor) + ".swissprot" in plot) == enabled
