"""Rescue admission regressions for reused and newly predicted evidence.

Discovery and external prediction are controlled on tiny genomes. Actual genomic
QC, terminal completion, donor/locus admission, annotation ownership, biological
ambiguity, serialization and normal atomic publication remain real. No old
implementation snapshot is imported or copied into repository tests.
"""
import copy
import hashlib
import json
from collections import Counter
from pathlib import Path

import pytest

from workflow.support import rescue_gene_models as rescue
from workflow.support import rescue_model_store


def make_inputs(directory):
    directory.mkdir()
    dna = {
        "chr_intact": "ATGAAATAA",
        "chr_partial": "ATGAAACCC",
        "chr_owned": "ATGCCCTAA",
        "chr_ambiguous": "ATGAAAATGCCCTAA",
        "chr_empty": "CCCCCCCCC",
    }
    genome = directory / "genome.fa"
    genome.write_text("".join(f">{contig}\n{sequence}\n" for contig, sequence in dna.items()))
    gff = directory / "source.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr_owned\tsource\tgene\t1\t9\t.\t+\t.\tID=original_owner\n"
        "chr_owned\tsource\tmRNA\t1\t9\t.\t+\t.\tID=original_tx;Parent=original_owner\n"
        "chr_owned\tsource\tCDS\t1\t9\t.\t+\t0\tParent=original_tx\n")
    cds = directory / "cds.fa"
    cds.write_text(">Target_species_original_owner\nATGCCCTAA\n")
    name = "Target_species"
    proteins = {
        name: {"self": "MP"},
        "Relative_species": {"good": "MK", "partial": "MKP", "owned": "MP",
                             "long": "MKMP", "short": "MP", "empty": "WW"},
        "Balanced_species": {"good": "MK"},
    }
    specs = [
        ("good_relative", "Relative_species", "good", "chr_intact", 0, 9, True),
        ("good_balanced", "Balanced_species", "good", "chr_intact", 0, 9, True),
        ("partial", "Relative_species", "partial", "chr_partial", 0, 9, True),
        ("owned", "Relative_species", "owned", "chr_owned", 0, 9, False),
        ("ambiguous_long", "Relative_species", "long", "chr_ambiguous", 0, 15, False),
        ("ambiguous_short", "Relative_species", "short", "chr_ambiguous", 6, 15, False),
        ("empty", "Relative_species", "empty", "chr_empty", 0, 9, True),
    ]
    regions, raw = [], {}
    for identifier, donor, query, contig, start, end, genome_only in specs:
        region = {"id": identifier, "target": name, "donor": donor, "query": query,
                  "seqid": contig, "start": 0, "end": len(dna[contig]),
                  "expected_start": 0, "expected_end": len(dna[contig]),
                  "expected_strand": "+", "comparison": "toy_frozen_comparison"}
        if genome_only:
            region.update(genome_only=True, placement="unanchored", orthology="unassigned",
                          expected_copy="unassigned", nomination={"reason": "no_target_match"})
        regions.append(region)
        if identifier == "empty":
            continue
        length = len(proteins[donor][query])
        identity = .8 if identifier.startswith("ambiguous") else .9
        raw[identifier] = {
            "query": identifier, "seqid": contig, "strand": "+", "cds": [[start, end, 0]],
            "frameshift": False, "coverage": 1.0, "identity": identity,
            "query_start": 0, "query_end": length, "query_length": length,
            "query_span_coverage": 1.0, "id": "raw_" + identifier,
            "paf_query": identifier, "paf_seqid": contig, "paf_strand": "+",
            "paf": f"##PAF\t{identifier}\t{length}\t0\t{length}\t+\t{contig}\t{len(dna[contig])}\t{start}\t{end}\t{length}\t{length}\t60\tcg:Z:{length}M",
            "evidence": copy.deepcopy(region), "search": "genome_fallback"}
    source = {"genome": str(genome), "gff": str(gff), "fasta": str(cds), "genetic_code": 1}
    hashes = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in (genome, gff, cds)}
    return {"name": name, "source": source, "files": hashes, "dna": dna,
            "proteins": proteins, "regions": regions, "raw": raw}


def execute(module, root, inputs, monkeypatch, *, cache, gemoma, fallback,
            storage="legacy", corrupt_cache=False, invalid_raw=False):
    root.mkdir()
    name = inputs["name"]
    for donor, proteins in inputs["proteins"].items():
        prepared = root / "prepared" / donor
        prepared.mkdir(parents=True)
        (prepared / "genes.pep").write_text(
            "".join(f">{query}\n{protein}\n" for query, protein in proteins.items()))
        (prepared / "positions.json").write_text("[]\n")
        (prepared / "receipt.json").write_text(json.dumps({"toy_prepared_donor": donor}) + "\n")
    request = {
        "sources": {name: copy.deepcopy(inputs["source"])},
        "files": copy.deepcopy(inputs["files"]),
        "gemoma_jar": "toy.jar" if gemoma else None,
        "output_storage": {"format": storage, "retain_search_inputs": False},
        "parameters": {"minimum_coverage": .95, "minimum_identity": .5,
                       "max_intron": 100, "unanchored_min_species": 2,
                       "genome_fallback": int(fallback), "terminal_max_extension": 30},
    }
    if cache != "none":
        request["prediction_cache"] = {"root": "verified_toy_cache"}
    plan = {"species": [name], "donors": {name: ["Relative_species", "Balanced_species"]},
            "nearest_references": {name: ["Relative_species"]}, "request": request}
    calls = {"cache_check": 0, "verify_sources": 0, "external": [], "gemoma_regions": []}
    regions = copy.deepcopy(inputs["regions"])
    rows = copy.deepcopy(inputs["raw"])
    covered = (set(rows) | {"empty"} if cache.startswith("complete") else
               {"good_relative", "partial", "empty"} if cache == "partial" else set())
    local_covered = (covered if cache == "complete" else set())
    if invalid_raw:
        rows["good_relative"]["cds"][0][1] = len(inputs["dna"]["chr_intact"]) + 3

    class Cache:
        search_contract = {"toy": "exact_complete_coverage_with_empty_result"}
        local_search_compatible = True
        genome_search_compatible = True

        def candidate_ids(self):
            return {r["id"] for r in regions
                    if not r.get("genome_only") and r["id"] in local_covered}

        def genome_candidate_ids(self):
            return set(covered)

        def iter_models(self):
            self.check()
            for region in regions:
                if region["id"] in covered and region["id"] in rows:
                    yield copy.deepcopy(rows[region["id"]])
            self.check()

        def check(self):
            calls["cache_check"] += 1
            if corrupt_cache and calls["cache_check"] == 3:
                raise ValueError("Verified toy prediction cache generation changed")

    def sources_guard(current, names, kinds):
        assert current is plan and names == [name] and kinds == ["genome", "fasta", "gff"]
        calls["verify_sources"] += 1
        for path, expected in inputs["files"].items():
            assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == expected

    def interval_predictions(tmp, windows, *args, **kwargs):
        local = []
        for _window, queries in windows.items():
            for region in queries:
                if region["id"] in rows:
                    item = copy.deepcopy(rows[region["id"]])
                    item.pop("evidence")
                    item.pop("search")
                    local.append(item)
        return local

    def external_run(command, tmp, label, stdout=None):
        calls["external"].append(label)
        if label == "miniprot_index":
            (tmp / "genome.mpi").write_bytes(b"owned toy index")
            return
        assert label == "miniprot_genome" and stdout is not None
        text = ["##gff-version 3\n"]
        for number, (identifier, _header, _protein) in enumerate(
                module.fasta_records(tmp / "unresolved.unique.fa"), 1):
            if identifier not in rows:
                # miniprot -u reports every unmapped query as an embedded PAF.
                text.append(f"##PAF\t{identifier}\t{len(_protein)}\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0\n")
                continue
            row = rows[identifier]
            model_id = f"MP{number:06d}"
            text.append(row["paf"] + "\n")
            start, end, phase = row["cds"][0]
            text.append(f"{row['seqid']}\tminiprot\tmRNA\t{start+1}\t{end}\t.\t+\t.\tID={model_id};Target={identifier} 1 {row['query_length']};Identity={row['identity']}\n")
            text.append(f"{row['seqid']}\tminiprot\tCDS\t{start+1}\t{end}\t.\t+\t{phase}\tParent={model_id}\n")
        Path(stdout).write_text("".join(text))

    def gemoma_spy(tmp, producer, current, source, selected, genome, validated, cpus):
        assert current is plan and cpus == 1 and source["genetic_code"] == 1
        calls["gemoma_regions"].append([r["id"] for r in selected])
        # Genuine unresolved selection enters this existing optional branch.
        # Its independent genomic terminal validation is exercised by the
        # unchanged public test_rescue_gemoma_completion.py suite.

    key = {"toy_immutable_source_sha": inputs["files"], "species": name}
    with monkeypatch.context() as patch:
        patch.delenv("GG_GENOME_INDEX_CACHE", raising=False)
        patch.setattr(module, "prepared", lambda *_: None)
        patch.setattr(module, "rescue_key", lambda *_: copy.deepcopy(key))
        patch.setattr(module, "candidates", lambda *_: copy.deepcopy(regions))
        patch.setattr(module, "nominate_genome_only_candidates", lambda *_: ([], {"nominated": 0}))
        patch.setattr(module, "verify_sources", sources_guard)
        patch.setattr(module, "verify_prediction_cache", lambda *_: Cache())
        patch.setattr(module, "search_intervals", interval_predictions)
        patch.setattr(module, "run", external_run)
        patch.setattr(module, "refine_gemoma", gemoma_spy)
        patch.setattr(module, "plan_digest", lambda *_: "b" * 64)
        output = module.rescue(root, plan, name, 1)
    assert calls["verify_sources"] == 2
    assert (output / "receipt.json").is_file()
    if storage == "compact":
        models = list(rescue_model_store.iter_models(output))
        partials = list(rescue_model_store.iter_partial_models(output))
        revisions = list(rescue_model_store.iter_revision_models(output))
    else:
        models = json.loads((output / "models.json").read_text())
        partials = json.loads((output / "partial_models.json").read_text())
        revisions = json.loads((output / "revision_candidates.json").read_text())
    outputs = {
        "models": models, "partials": partials, "revisions": revisions,
        "placement": json.loads((output / "placement_audit.json").read_text()),
        "candidates": json.loads((output / "candidates.json").read_text()),
        "reuse": json.loads((output / "prediction_reuse.json").read_text()),
        "memo": json.loads((output / "raw_validation_memo.json").read_text()),
        "audit": (output / "audit.tsv").read_bytes(),
        "quality_flags": (output / "quality_flags.tsv").read_bytes(),
        "search_inputs": json.loads((output / "search_inputs.json").read_text()),
        "status_counts": dict(Counter(m["status"] for m in models)),
    }
    return outputs, calls


@pytest.mark.parametrize("cache", ["none", "partial", "complete", "complete_genome_only"])
@pytest.mark.parametrize("gemoma", [False, True])
@pytest.mark.parametrize("fallback", [False, True])
@pytest.mark.parametrize("storage", ["legacy", "compact"])
def test_rescue_admission_preserves_intact_partial_owned_ambiguous_and_empty_evidence(
        tmp_path, monkeypatch, cache, gemoma, fallback, storage):
    inputs = make_inputs(tmp_path / "inputs")
    before = {p: Path(p).read_bytes() for p in inputs["files"]}
    results, calls = execute(rescue, tmp_path / "output", inputs, monkeypatch,
                             cache=cache, gemoma=gemoma, fallback=fallback, storage=storage)
    # Explicit oracle for these supplied raw/cached/search alignments, including
    # ordering and duplicate evidence. The same oracle applies to both formats.
    expected_order = {
        ("none", False): ["owned", "ambiguous_long", "ambiguous_short"],
        ("none", True): ["owned", "ambiguous_long", "ambiguous_short",
                        "good_relative", "good_balanced", "partial",
                        "owned", "ambiguous_long", "ambiguous_short"],
        ("partial", False): ["good_relative", "partial", "owned", "ambiguous_long", "ambiguous_short"],
        ("partial", True): ["good_relative", "partial", "owned", "ambiguous_long", "ambiguous_short",
                           "good_balanced", "owned", "ambiguous_long", "ambiguous_short"],
        ("complete", False): ["good_relative", "good_balanced", "partial", "owned", "ambiguous_long", "ambiguous_short"],
        ("complete", True): ["good_relative", "good_balanced", "partial", "owned", "ambiguous_long", "ambiguous_short"],
        ("complete_genome_only", False): ["good_relative", "good_balanced", "partial", "owned", "ambiguous_long", "ambiguous_short",
                                        "owned", "ambiguous_long", "ambiguous_short"],
        ("complete_genome_only", True): ["good_relative", "good_balanced", "partial", "owned", "ambiguous_long", "ambiguous_short",
                                       "owned", "ambiguous_long", "ambiguous_short"],
    }[cache, fallback]
    models = results["models"]
    assert [m["query"] for m in models] == expected_order
    accepted = [m for m in models if m["status"] == "accepted"]
    independently_supported = cache in ("complete", "complete_genome_only") or fallback
    assert len(accepted) == int(independently_supported)
    expected_status = ({"accepted": 1, "duplicate_support": 1, "unresolved": len(models) - 2}
                       if independently_supported else {"unresolved": len(models)})
    assert results["status_counts"] == expected_status
    if accepted:
        model = accepted[0]
        assert model["seqid"] == "chr_intact" and model["strand"] == "+"
        assert model["sequence"] == "ATGAAATAA" and model["cds"] == [[0, 9, 0]]
        assert model["problems"] == []
        assert {e["donor"] for e in model["support"]} == {"Relative_species", "Balanced_species"}
        assert model["placement_evidence"]["orthology"] == "unassigned"
        assert model["placement_evidence"]["expected_copy"] == "unassigned"
        assert model["quality_evidence"]["translation_initiation"] == "not_established"
        assert model["quality_evidence"]["native_terminal_completeness"] == "not_established"
    elif cache == "partial":
        model = next(m for m in models if m["query"] == "good_relative")
        assert model["sequence"] == "ATGAAATAA" and model["status"] == "unresolved"
        assert model["placement_evidence"]["reason"] == "insufficient_unique_donor_support"
        assert "unanchored_genome_search" in model["problems"]
    owned = [m for m in models if m["query"] == "owned"]
    assert owned and all(m["status"] == "unresolved" for m in owned)
    assert all(m["revision_owner_ids"] == ["original_owner"] for m in owned)
    assert all("overlap_existing_annotation" in m["problems"] for m in owned)
    assert len(results["revisions"]) == 1
    assert results["revisions"][0]["revision_owner_ids"] == ["original_owner"]
    ambiguous = [m for m in models if m["query"].startswith("ambiguous") and
                 m["seqid"] == "chr_ambiguous"]
    assert {tuple(map(tuple, m["cds"])) for m in ambiguous} == {((0, 15, 0),), ((6, 15, 0),)}
    assert all(m["status"] == "unresolved" and
               "ambiguous_coding_path_representative" in m["problems"] for m in ambiguous)
    # The same MP donor protein also maps to the original annotated locus.
    # Genome-query deduplication restores that placement to its short alias.
    # It is outside the alias's anchor and cannot become an accepted path or
    # contribute support to an existing-model revision.
    off_locus = [m for m in models if m["query"].startswith("ambiguous") and
                 m["seqid"] != "chr_ambiguous"]
    assert len(off_locus) == int(fallback and cache in ("none", "partial"))
    for model in off_locus:
        assert model["query"] == "ambiguous_short" and model["seqid"] == "chr_owned"
        assert model["cds"] == [[0, 9, 0]] and model["sequence"] == "ATGCCCTAA"
        assert model["status"] == "unresolved" and model["revision_owner_ids"] == ["original_owner"]
        assert set(model["problems"]) == {"outside_expected_synteny_interval", "overlap_existing_annotation"}
    assert all(row["id"] == "owned" for revision in results["revisions"] for row in revision["support"])
    partials = results["partials"]
    assert len(partials) == int(cache != "none" or fallback)
    if partials:
        model = partials[0]
        assert model["query"] == "partial" and model["sequence"] == "ATGAAACCC"
        assert model["status"] == "unresolved" and "missing_stop" in model["problems"]
        assert model["partial_evidence"]["partial"] is True
        assert model["partial_evidence"]["representative_eligible"] is False
    empty = [line.split("\t") for line in results["audit"].decode().splitlines()
             if line.startswith("empty\t")]
    assert empty == [["empty", "empty", "unresolved", "no_alignment", "", "", ""]]
    assert "empty" not in {m["query"] for m in models}
    assert len(results["quality_flags"].decode().splitlines()) == len(models) + 1
    expected_gemoma = ([["owned", "ambiguous_long", "ambiguous_short"]]
                       if gemoma and cache in ("none", "partial") else [])
    assert calls["gemoma_regions"] == expected_gemoma
    assert {p: Path(p).read_bytes() for p in inputs["files"]} == before


def test_complete_cache_still_fails_closed_if_prediction_generation_changes(tmp_path, monkeypatch):
    inputs = make_inputs(tmp_path / "inputs")
    with pytest.raises(ValueError, match="prediction cache generation changed"):
        execute(rescue, tmp_path / "output", inputs, monkeypatch,
                cache="complete", gemoma=False, fallback=True, corrupt_cache=True)
    assert not (tmp_path / "output/rescued" / inputs["name"] / "receipt.json").exists()


def test_complete_cache_still_rejects_out_of_genome_raw_prediction(tmp_path, monkeypatch):
    inputs = make_inputs(tmp_path / "inputs")
    with pytest.raises(ValueError, match="Predicted CDS outside genome"):
        execute(rescue, tmp_path / "output", inputs, monkeypatch,
                cache="complete", gemoma=False, fallback=True, invalid_raw=True)
    assert not (tmp_path / "output/rescued" / inputs["name"] / "receipt.json").exists()
