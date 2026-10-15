"""Source-bound convention decoding preserves biological CDS and strict gates."""
import copy
import hashlib
import json
from pathlib import Path

import pytest
from Bio.Seq import Seq

from workflow.tests.test_gene_model_catalog import MODULE, digest, fixture, row
from workflow.tests.test_gene_model_refinement import refinement


def model_fixture(tmp_path, *, strand="+", offset=0, tail="", source_mode="full",
                  source="IGDBv2", phases=None, coding="ATGAAACCCTAA", two_transcripts=False):
    """Two blocks whose nonzero continuation distinguishes the conventions."""
    biological = "TC"[:offset] + coding + tail
    spliced_genome = biological[:4] + "NNNNNN" + biological[4:]
    end = len(spliced_genome)
    annotations = row("gene", 1, end, "g", strand=strand)
    source_phases = phases if phases is not None else [(-offset) % 3, ((-offset) + 4) % 3]
    blocks = []
    for tx in ("t", "other") if two_transcripts else ("t",):
        annotations += row("mRNA", 1, end, tx, parent="g", strand=strand)
        for index, (a, b) in enumerate([(1, 4), (11, end)]):
            if strand == "-":
                a, b = end + 1 - b, end + 1 - a
            annotations += row("CDS", a, b, tx + str(index), parent=tx,
                               phase=str(source_phases[index]), strand=strand)
            if tx == "t":
                blocks.append([a - 1, b, source_phases[index]])
    annotations = annotations.replace("\ts\t", "\t" + source + "\t")
    supplied = biological
    identifier = "t"
    if source_mode == "trimmed":
        supplied = biological[offset:len(biological) - len(tail)] if tail else biological[offset:]
    elif source_mode == "mismatch":
        supplied = biological.replace("AAA", "CCC")
    elif source_mode == "unbound":
        identifier = "unknown"
    elif source_mode == "gene":
        identifier = "g"
    supplied = supplied[0].lower() + supplied[1:]
    genomic = spliced_genome if strand == "+" else str(Seq(spliced_genome).reverse_complement())
    paths = fixture(tmp_path, ">" + identifier + " raw source header\n" + supplied + "\n", annotations, genomic)
    return paths, biological, blocks, supplied


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("offset", [0, 1, 2])
@pytest.mark.parametrize("tail", ["", "C"])
def test_bound_complementary_models_preserve_all_internal_bases_and_partials(tmp_path, strand, offset, tail):
    paths, biological, blocks, supplied = model_fixture(tmp_path, strand=strand, offset=offset, tail=tail)
    before = [digest(path) for path in paths]
    catalog = MODULE.build_catalog("Species_one", *paths)
    candidate = catalog["loci"][0]["candidates"][0]
    assert [digest(path) for path in paths] == before
    assert candidate["cds"] == biological and candidate["source_blocks"] == blocks
    assert candidate["blocks"] == [[a, b, (-phase) % 3] for a, b, phase in blocks]
    assert candidate["protein"] == "MKP"
    quality = candidate["quality"]
    assert quality["usable"] and not quality["phase_conflict"]
    assert quality["translation_offset"] == offset
    assert quality["partial"] is bool(offset or tail)
    assert quality["valid_orf"] is (not bool(offset or tail))
    assert quality["incomplete_codon"] is bool(tail)
    evidence = candidate["source_phase_interpretation"]
    assert evidence["binding"] == "source_transcript"
    assert evidence["source_convention"] == "complementary" and evidence["canonical_convention"] == "gff3"
    assert evidence["source_first_phase"] == (-offset) % 3 and evidence["translation_offset"] == offset
    assert evidence["source_cds_sha256"] == hashlib.sha256(supplied.upper().encode()).hexdigest()
    assert evidence["gff_sha256"] == catalog["sources"]["gff"]["sha256"]
    assert evidence["genome_sha256"] == catalog["sources"]["genome"]["sha256"]
    assert candidate["source_cds"][0]["cds"] == supplied
    assert candidate["source_cds"][0]["sequence_agreement"]
    assert catalog["fasta_mapping"][0]["sequence_agreement"]
    frozen = copy.deepcopy(candidate)
    candidate["quality"] = MODULE.validate_candidate(candidate)
    assert candidate == frozen
    again = MODULE.build_catalog("Species_one", *paths)["loci"][0]["candidates"][0]
    assert again == frozen
    encoded = [catalog["sources"]["genome"]["sha256"], candidate["seqid"], candidate["strand"],
               candidate["blocks"], hashlib.sha256(biological.encode()).hexdigest()]
    assert candidate["coding_key"] == hashlib.sha256(json.dumps(encoded, separators=(",", ":")).encode()).hexdigest()
    written = MODULE.write_catalog(catalog, tmp_path / "catalogue")
    archived = json.loads(Path(written["catalog"]).read_text())["loci"][0]["candidates"][0]
    assert archived == frozen


@pytest.mark.parametrize("offset", [1, 2])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_existing_formatter_partial_cds_proves_the_bound_complementary_offset(tmp_path, offset, strand):
    paths, biological, _blocks, supplied = model_fixture(
        tmp_path, offset=offset, strand=strand, tail="C", source_mode="trimmed")
    value = MODULE.build_catalog("Species_one", *paths)
    candidate = value["loci"][0]["candidates"][0]
    assert candidate["cds"] == biological and candidate["protein"] == "MKP"
    assert candidate["quality"]["usable"] and candidate["quality"]["partial"]
    assert candidate["source_phase_interpretation"]["sequence_convention"] == "annotated_partial_frame"
    assert candidate["source_cds"][0]["cds"] == supplied
    assert candidate["corrected_cds_length"] == 12


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("offset", [0, 1, 2])
def test_shared_analysis_derivative_reproduces_interpreted_partial_peptide_and_gff(tmp_path, strand, offset):
    paths, biological, _blocks, _supplied = model_fixture(tmp_path, strand=strand, offset=offset, tail="C")
    before = [digest(path) for path in paths]
    candidate = MODULE.build_catalog("Species_one", *paths)["loci"][0]["candidates"][0]
    frozen = copy.deepcopy(candidate)
    analysis = refinement.analysis_coding_candidate(candidate)
    assert candidate == frozen
    assert analysis["cds"] == biological[offset:-1] == "ATGAAACCCTAA"
    assert analysis["protein"] == str(Seq(analysis["cds"]).translate()).removesuffix("*") == "MKP"
    assert analysis["analysis"] == {"removed_first_bases": offset, "removed_final_bases": 1,
                                    "source_cds_length": len(biological), "analysis_cds_length": 12}
    assert analysis["source_blocks"] == frozen["source_blocks"]
    assert analysis["source_phase_interpretation"] == frozen["source_phase_interpretation"]
    selected = refinement.selected_gff_rows(paths[1].read_text(), {candidate["source_transcript_id"]})
    derived = refinement.analysis_gff_rows(selected, [analysis])
    exported = []
    for line in derived.splitlines():
        fields = line.split("\t")
        if len(fields) == 9 and fields[2] == "CDS":
            exported.append([int(fields[3]) - 1, int(fields[4]), int(fields[7])])
    exported.sort(reverse=strand == "-")
    assert exported == analysis["blocks"]
    assert MODULE.phase_conventions([{"start": a, "end": b, "phase": str(p)} for a, b, p in exported]) >= {"gff3"}
    with MODULE.indexed_genome(paths[2]) as genome:
        assert MODULE.reconstruct_sequence(exported, "chr1", strand, genome) == analysis["cds"]
    assert MODULE.validate_candidate(candidate) == frozen["quality"] and candidate == frozen
    assert [digest(path) for path in paths] == before


@pytest.mark.parametrize("source_mode", ["mismatch", "unbound"])
def test_complementary_adoption_requires_existing_source_identity_and_sequence_agreement(tmp_path, source_mode):
    paths, biological, blocks, _supplied = model_fixture(tmp_path, source_mode=source_mode)
    value = MODULE.build_catalog("Species_one", *paths)
    candidate = value["loci"][0]["candidates"][0]
    assert candidate["blocks"] == blocks and candidate["cds"] == biological
    assert candidate["quality"]["phase_conflict"] and not candidate["quality"]["usable"]
    assert "source_blocks" not in candidate and "source_phase_interpretation" not in candidate
    if source_mode == "mismatch":
        assert candidate["quality"]["sequence_mismatch"]
    else:
        assert value["fasta_mapping"][0]["mapping_status"] == "unmapped"


@pytest.mark.parametrize("coding,pseudogene,flag", [("ATGTGACCCTAA", False, "internal_stop"),
                                                  ("ATGNNNCCCTAA", False, "ambiguous"),
                                                  ("ATGAAACCCTAA", True, "annotated_pseudogene")])
def test_interpreting_source_phases_does_not_clear_independent_quality_failures(tmp_path, coding, pseudogene, flag):
    paths, biological, blocks, _supplied = model_fixture(tmp_path, coding=coding)
    if pseudogene:
        paths[1].write_text(paths[1].read_text().replace("ID=g\n", "ID=g;pseudo=true\n"))
    candidate = MODULE.build_catalog("Species_one", *paths)["loci"][0]["candidates"][0]
    assert candidate["cds"] == biological and candidate["source_blocks"] == blocks
    assert "source_phase_interpretation" in candidate and not candidate["quality"]["phase_conflict"]
    assert candidate["quality"][flag] and not candidate["quality"]["usable"]


def test_generic_and_predicted_validation_does_not_accept_complementary_frames():
    for origin in ("original", "predicted"):
        candidate = {"origin": origin, "strand": "+", "cds": "ATGAAATAA", "blocks": [[0, 4, 0], [10, 15, 1]]}
        quality = MODULE.validate_candidate(candidate)
        assert candidate["protein"] == "MK" and quality["phase_conflict"] and not quality["usable"]
        assert candidate["blocks"] == [[0, 4, 0], [10, 15, 1]]


@pytest.mark.parametrize("generated_source", ["GeneGalleon", "genegalleon_refinement", "genegalleon_rescue"])
def test_generated_annotation_frames_stay_strict_even_with_source_agreement(tmp_path, generated_source):
    paths, biological, blocks, _supplied = model_fixture(tmp_path, source=generated_source)
    candidate = MODULE.build_catalog("Species_one", *paths)["loci"][0]["candidates"][0]
    assert candidate["blocks"] == blocks and candidate["cds"] == biological
    assert candidate["quality"]["source_sequence_agreement"]
    assert candidate["quality"]["phase_conflict"] and not candidate["quality"]["usable"]
    assert "source_phase_interpretation" not in candidate


def test_each_model_in_mixed_native_and_generated_gff_keeps_its_own_convention(tmp_path):
    biological = "ATGAAATAA"
    genome = "ATGANNNNNNAATAA" * 3
    annotations, records = "", ""
    for index, (tx, source, phase) in enumerate([("native", "IGDBv2", "1"),
                                                ("rescue", "genegalleon_rescue", "2"),
                                                ("generated_bad", "GeneGalleon", "1")]):
        shift = index * 15
        part = row("gene", shift + 1, shift + 15, "g" + str(index))
        part += row("mRNA", shift + 1, shift + 15, tx, parent="g" + str(index))
        part += row("CDS", shift + 1, shift + 4, tx + "a", parent=tx, phase="0")
        part += row("CDS", shift + 11, shift + 15, tx + "b", parent=tx, phase=phase)
        annotations += part.replace("\ts\t", "\t" + source + "\t")
        records += ">" + tx + "\n" + biological + "\n"
    value = MODULE.build_catalog("Species_one", *fixture(tmp_path, records, annotations, genome))
    candidates = {c["source_transcript_id"]: c for locus in value["loci"] for c in locus["candidates"]}
    assert candidates["native"]["quality"]["valid_orf"]
    assert "source_phase_interpretation" in candidates["native"]
    assert candidates["rescue"]["quality"]["valid_orf"] and "source_blocks" not in candidates["rescue"]
    assert candidates["generated_bad"]["quality"]["phase_conflict"] and not candidates["generated_bad"]["quality"]["usable"]
    assert "source_blocks" not in candidates["generated_bad"]


def test_valid_gff3_is_unchanged_including_nonzero_partial_offset(tmp_path):
    paths, biological, _blocks, _supplied = model_fixture(tmp_path, offset=1, phases=[1, 0])
    candidate = MODULE.build_catalog("Species_one", *paths)["loci"][0]["candidates"][0]
    assert candidate["quality"]["usable"] and candidate["protein"] == "MKP"
    assert candidate["cds"] == biological and candidate["quality"]["translation_offset"] == 1
    assert "source_blocks" not in candidate and "source_phase_interpretation" not in candidate


def test_ambiguous_single_block_convention_cannot_choose_the_stop_free_reading_frame(tmp_path):
    biological = "CCATGAATGAATAA"
    annotations = row("gene", 1, len(biological), "g") + row("mRNA", 1, len(biological), "t", parent="g")
    annotations += row("CDS", 1, len(biological), "c", parent="t", phase="1")
    candidate = MODULE.build_catalog("Species_one", *fixture(tmp_path, ">t\n" + biological + "\n", annotations, biological))["loci"][0]["candidates"][0]
    assert candidate["quality"]["internal_stop"] and not candidate["quality"]["usable"]
    assert candidate["quality"]["translation_offset"] == 1 and candidate["cds"] == biological
    assert str(Seq(biological[2:]).translate()).removesuffix("*") == "MNE"
    assert "source_phase_interpretation" not in candidate and "source_blocks" not in candidate


@pytest.mark.parametrize("phases", [[0, 0], [0, "."], [".", 1]])
def test_incompatible_or_unknown_partial_phase_does_not_authorise_complementary_adoption(tmp_path, phases):
    paths, _biological, blocks, _supplied = model_fixture(tmp_path, phases=phases, tail="C")
    candidate = MODULE.build_catalog("Species_one", *paths)["loci"][0]["candidates"][0]
    assert not candidate["quality"]["usable"] and "source_phase_interpretation" not in candidate
    assert "source_blocks" not in candidate
    assert [[a, b] for a, b, _p in candidate["blocks"]] == [[a, b] for a, b, _p in blocks]


def test_one_unambiguous_bound_coding_path_preserves_transcript_identity_ambiguity(tmp_path):
    paths, _biological, blocks, _supplied = model_fixture(tmp_path, source_mode="gene", two_transcripts=True)
    value = MODULE.build_catalog("Species_one", *paths)
    locus = value["loci"][0]
    assert value["fasta_mapping"][0]["mapping_status"] == "ambiguous"
    assert locus["source_baseline_candidate_id"] == ""
    assert locus["source_baseline_coding_candidate_id"] in {c["candidate_id"] for c in locus["candidates"]}
    assert len({c["coding_key"] for c in locus["candidates"]}) == 1
    for candidate in locus["candidates"]:
        assert candidate["quality"]["valid_orf"] and candidate["source_blocks"] == blocks
        assert candidate["source_phase_interpretation"]["binding"] == "unique_source_coding_path"
        assert candidate["source_fasta_ids"] == []
        assert MODULE.validate_candidate(copy.deepcopy(candidate)) == candidate["quality"]


@pytest.mark.parametrize("mutation", ["cds", "blocks", "source_blocks", "source_cds", "source_case", "raw_hash",
                                     "origin", "sources", "offset", "coding_key", "transcript", "gene", "deleted_proof"])
def test_repeated_validation_rejects_changed_interpretation_evidence(tmp_path, mutation):
    paths, _biological, _blocks, _supplied = model_fixture(tmp_path)
    candidate = MODULE.build_catalog("Species_one", *paths)["loci"][0]["candidates"][0]
    if mutation == "cds":
        candidate["cds"] = candidate["cds"].replace("AAA", "CCC")
    elif mutation == "blocks":
        candidate["blocks"][1][2] = 1
    elif mutation == "source_blocks":
        candidate["source_blocks"][0][0] += 1
    elif mutation == "source_cds":
        candidate["source_cds"][0]["cds"] = "ATGCCC" + candidate["source_cds"][0]["cds"][6:]
    elif mutation == "source_case":
        candidate["source_cds"][0]["cds"] = candidate["source_cds"][0]["cds"].upper()
    elif mutation == "raw_hash":
        candidate["source_cds"][0]["sha256"] = "0" * 64
    elif mutation == "origin":
        candidate["origin"] = "predicted"
    elif mutation == "sources":
        candidate["source_annotation_sources"].append("GeneGalleon")
    elif mutation == "offset":
        candidate["source_phase_interpretation"]["translation_offset"] = 1
    elif mutation == "transcript":
        candidate["source_transcript_id"] = "another_transcript"
    elif mutation == "gene":
        candidate["source_gene_id"] = "another_gene"
    elif mutation == "deleted_proof":
        del candidate["source_phase_interpretation"]
    else:
        candidate["coding_key"] = "0" * 64
    quality = MODULE.validate_candidate(candidate)
    assert not quality["usable"] and quality["structure_problem"] == "source_phase_interpretation_mismatch"
