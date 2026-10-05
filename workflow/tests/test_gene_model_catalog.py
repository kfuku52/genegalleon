"""Complete isoform preservation and independent CDS admission contracts."""

import gzip
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import pytest
from Bio.Seq import Seq

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
SPEC = importlib.util.spec_from_file_location("gene_model_catalog", SUPPORT / "gene_model_catalog.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def fixture(tmp_path, records, annotations, genome_sequence, *, compressed=False):
    suffix = ".gz" if compressed else ""
    paths = [tmp_path / (name + suffix) for name in ("source.cds.fa", "source.gff", "source.genome.fa")]
    data = [records, "##gff-version 3\n" + annotations, ">chr1\n" + genome_sequence + "\n"]
    for path, content in zip(paths, data, strict=True):
        if compressed:
            with gzip.open(path, "wt") as handle:
                handle.write(content)
        else:
            path.write_text(content)
    return paths


def row(feature, start, end, identifier, *, parent="", phase=".", strand="+"):
    attrs = "ID=" + identifier if identifier else ""
    if parent:
        attrs += (";" if attrs else "") + "Parent=" + parent
    return f"chr1\ts\t{feature}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}\n"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.mark.parametrize("compressed", [False, True])
def test_duplicate_genome_contig_names_cannot_choose_one_conflicting_assembly_record(tmp_path, compressed):
    annotations = row("gene", 1, 12, "g") + row("mRNA", 1, 12, "t", parent="g")
    annotations += row("CDS", 1, 12, "c", parent="t", phase="0")
    paths = fixture(tmp_path, ">t\nATGAAACCCTAA\n", annotations,
                    "ATGAAACCCTAA\n>chr1\nATGCCCCCCTAA", compressed=compressed)
    before = [digest(path) for path in paths]
    with pytest.raises(ValueError, match="duplicate sequence"):
        MODULE.build_catalog("Species_one", *paths)
    assert [digest(path) for path in paths] == before
    assert not any(path.name.endswith((".fai", ".gzi")) for path in tmp_path.iterdir())


@pytest.mark.parametrize("ancestor", ["gene", "mRNA"])
@pytest.mark.parametrize("axis", ["seqid", "strand"])
def test_coding_path_cannot_disagree_with_explicit_ancestor_axis(tmp_path, ancestor, axis):
    annotations = row("gene", 1, 12, "g") + row("mRNA", 1, 12, "t", parent="g")
    annotations += row("CDS", 1, 12, "c", parent="t", phase="0")
    lines = []
    for line in annotations.splitlines():
        fields = line.split("\t")
        if fields[2] == ancestor:
            fields[0 if axis == "seqid" else 6] = "chr2" if axis == "seqid" else "-"
        lines.append("\t".join(fields) + "\n")
    paths = fixture(tmp_path, ">t\nATGAAACCCTAA\n", "".join(lines), "ATGAAACCCTAA\n>chr2\nATGAAACCCTAA")
    catalog = MODULE.build_catalog("Species_one", *paths)
    candidate = catalog["loci"][0]["candidates"][0]
    assert candidate["quality"]["structure_problem"] == "ancestor_contig_or_strand_mismatch"
    assert not candidate["quality"]["usable"]
    assert not candidate["quality"]["valid_orf"]
    assert candidate["source_cds"][0]["cds"] == "ATGAAACCCTAA"
    written = MODULE.write_catalog(catalog, tmp_path / "catalog")
    assert Path(written["protein"]).read_text() == ""


@pytest.mark.parametrize("compressed", [False, True])
def test_all_transcript_identities_survive_identical_coding_paths_and_longest_source(tmp_path, compressed):
    genome = list("N" * 60)
    genome[0:4], genome[10:15], genome[20:28] = "ATGA", "AATAA", "AACCTTAA"
    annotations = row("gene", 1, 28, "g1")
    for transcript, second_end, phase in (("t1", 15, "2"), ("t2", 28, "2"), ("t3", 15, "2")):
        annotations += row("mRNA", 1, second_end, transcript, parent="g1")
        annotations += row("CDS", 1, 4, transcript + "-a", parent=transcript, phase="0")
        annotations += row("CDS", 11 if transcript != "t2" else 21, second_end,
                           transcript + "-b", parent=transcript, phase=phase)
    paths = fixture(tmp_path, ">g1\nATGAAACCTTAA\n", annotations, "".join(genome), compressed=compressed)
    before = [digest(path) for path in paths]
    catalog = MODULE.build_catalog("Species_one", *paths)
    assert catalog["summary"] == {"loci": 1, "candidates": 3, "usable_candidates": 3,
                                  "coding_paths": 2, "fasta_mapping": {"mapped": 1}}
    locus = catalog["loci"][0]
    assert locus["gene_id"] == "Species_one_g1"
    assert locus["source_baseline_candidate_id"] == "Species_one_t2"
    candidates = {candidate["source_transcript_id"]: candidate for candidate in locus["candidates"]}
    assert candidates["t1"]["coding_key"] == candidates["t3"]["coding_key"]
    assert candidates["t1"]["cds"] == "ATGAAATAA"
    assert candidates["t1"]["protein"] == "MK"
    assert candidates["t1"]["blocks"] == [[0, 4, 0], [10, 15, 2]]
    assert candidates["t1"]["junctions"] == [[4, 10, 1]]
    assert candidates["t2"]["corrected_cds_length"] == 12
    written = MODULE.write_catalog(catalog, tmp_path / "output")
    assert Path(written["cds"]).read_text().count(">") == 3
    assert Path(written["protein"]).read_text().count(">") == 3
    assert json.loads(Path(written["catalog"]).read_text()) == catalog
    streamed = [json.loads(line) for line in Path(written["loci"]).read_text().splitlines()]
    metadata = json.loads(Path(written["metadata"]).read_text())
    assert "loci" not in metadata
    excluded = [json.loads(line) for line in Path(written["excluded"]).read_text().splitlines()]
    assert {**{key: value for key, value in metadata.items() if key != "excluded_loci_summary"},
            "loci": streamed, "excluded_loci": excluded} == catalog
    assert [digest(path) for path in paths] == before
    assert not any(path.name.endswith((".fai", ".gzi")) for path in tmp_path.iterdir())


def test_minus_strand_phase_and_partial_source_have_no_padding(tmp_path):
    genome = list("N" * 80)
    genome[30:35], genome[40:44] = str(Seq("AATAA").reverse_complement()), str(Seq("ATGA").reverse_complement())
    genome[60:70] = "AATGAAATAA"
    annotations = row("gene", 31, 44, "minus", strand="-")
    annotations += row("mRNA", 31, 44, "mt", parent="minus", strand="-")
    annotations += row("CDS", 31, 35, "m-a", parent="mt", phase="2", strand="-")
    annotations += row("CDS", 41, 44, "m-b", parent="mt", phase="0", strand="-")
    annotations += row("gene", 61, 70, "partial") + row("mRNA", 61, 70, "pt", parent="partial")
    annotations += row("CDS", 61, 70, "p-c", parent="pt", phase="1")
    paths = fixture(tmp_path, ">mt\nATGAAATAA\n>pt\nATGAAATAA\n", annotations, "".join(genome))
    catalog = MODULE.build_catalog("Species_one", *paths)
    candidates = {candidate["source_transcript_id"]: candidate for locus in catalog["loci"]
                  for candidate in locus["candidates"]}
    assert candidates["mt"]["cds"] == "ATGAAATAA"
    assert candidates["mt"]["blocks"] == [[40, 44, 0], [30, 35, 2]]
    assert candidates["mt"]["junctions"] == [[40, 35, 1]]
    assert candidates["mt"]["quality"]["valid_orf"]
    assert candidates["pt"]["cds"] == "AATGAAATAA"
    assert candidates["pt"]["protein"] == "MK"
    assert candidates["pt"]["quality"]["partial"]
    assert candidates["pt"]["quality"]["usable"]
    assert not candidates["pt"]["quality"]["valid_orf"]
    assert candidates["pt"]["corrected_cds_length"] == 9
    assert all(row["sequence_agreement"] for row in catalog["fasta_mapping"])


@pytest.mark.parametrize("code,protein,usable", [(1, "M*", False), (4, "MW", True)])
def test_genetic_code_preserves_real_internal_stop(tmp_path, code, protein, usable):
    annotations = row("gene", 1, 9, "g") + row("mRNA", 1, 9, "t", parent="g")
    annotations += row("CDS", 1, 9, "c", parent="t", phase="0")
    paths = fixture(tmp_path, ">t\nATGTGATAA\n", annotations, "ATGTGATAA")
    candidate = MODULE.build_catalog("Species_one", *paths, genetic_code=code)["loci"][0]["candidates"][0]
    assert candidate["cds"] == "ATGTGATAA"
    assert candidate["protein"] == protein
    assert candidate["quality"]["usable"] is usable
    assert candidate["quality"]["internal_stop"] is (code == 1)


@pytest.mark.parametrize("code,dual", [(27, "TGA"), (28, "TAA"), (31, "TAA")])
@pytest.mark.parametrize("internal", [False, True])
@pytest.mark.parametrize("phase", ["0", "."])
def test_context_dependent_dual_codons_withhold_coding_admission_and_phase_inference(tmp_path, code, dual, internal, phase):
    sequence = "ATG" + (dual if internal else "AAA") + dual
    annotations = row("gene", 1, 9, "g") + row("mRNA", 1, 9, "t", parent="g")
    annotations += row("CDS", 1, 9, "c", parent="t", phase=phase)
    paths = fixture(tmp_path, ">t\n" + sequence + "\n", annotations, sequence)
    before = [digest(path) for path in paths]
    catalog = MODULE.build_catalog("Species_one", *paths, genetic_code=code)
    candidate = catalog["loci"][0]["candidates"][0]
    assert candidate["cds"] == sequence
    assert candidate["source_cds"][0]["cds"] == sequence
    quality = candidate["quality"]
    assert quality["translation_uncertain"]
    assert quality["translation_uncertain_codons"] == [dual]
    assert not quality["usable"] and not quality["valid_orf"]
    assert not quality["phase_inferred"]
    assert quality["source_sequence_agreement"]
    assert candidate["blocks"][0][2] == (0 if phase == "0" else -1)
    written = MODULE.write_catalog(catalog, tmp_path / "catalog")
    assert Path(written["protein"]).read_text() == ""
    assert sequence in Path(written["cds"]).read_text()
    assert [digest(path) for path in paths] == before


@pytest.mark.parametrize("code,sequence,expected", [(1, "TTGAAATAA", "LK"), (11, "GTGAAATAA", "VK"),
                                                   (4, "ATGTGATAA", "MW"), (27, "ATGAAACCC", "MKP")])
def test_non_dual_codons_keep_existing_context_free_translation_contract(code, sequence, expected):
    candidate = {"strand": "+", "blocks": [[0, len(sequence), 0]], "cds": sequence}
    quality = MODULE.validate_candidate(candidate, code)
    assert quality["usable"]
    assert not quality["translation_uncertain"]
    assert quality["translation_uncertain_codons"] == []
    assert quality["translation_convention"] == "context_free_codons_without_terminal_definite_stop"
    assert candidate["protein"] == expected


def test_gene_only_fasta_with_identical_paths_does_not_invent_identity(tmp_path):
    annotations = row("gene", 1, 9, "g")
    for transcript in ("t1", "t2"):
        annotations += row("mRNA", 1, 9, transcript, parent="g")
        annotations += row("CDS", 1, 9, "c-" + transcript, parent=transcript, phase="0")
    paths = fixture(tmp_path, ">g\nATGAAATAA\n", annotations, "ATGAAATAA")
    catalog = MODULE.build_catalog("Species_one", *paths)
    assert catalog["fasta_mapping"][0]["mapping_status"] == "ambiguous"
    assert catalog["fasta_mapping"][0]["candidate_ids"] == []
    assert catalog["loci"][0]["source_baseline_candidate_id"] == ""
    assert all(candidate["quality"]["usable"] for candidate in catalog["loci"][0]["candidates"])


@pytest.mark.parametrize("attribute", ["Name", "Accession"])
def test_transcript_specific_provider_alias_binds_identical_coding_paths_exactly(tmp_path, attribute):
    annotations = row("gene", 1, 9, "g")
    for transcript in ("t1", "t2"):
        annotations += row("mRNA", 1, 9, transcript, parent="g").rstrip("\n") + ";" + attribute + "=RNA" + transcript + "\n"
        annotations += row("CDS", 1, 9, "c-" + transcript, parent=transcript, phase="0")
    paths = fixture(tmp_path, ">Species_one_RNAt2\nATGAAATAA\n", annotations, "ATGAAATAA")
    catalog = MODULE.build_catalog("Species_one", *paths)
    assert catalog["fasta_mapping"][0]["mapping_status"] == "mapped"
    assert catalog["fasta_mapping"][0]["candidate_ids"] == ["Species_one_t2"]
    assert catalog["loci"][0]["source_baseline_candidate_id"] == "Species_one_t2"


def test_explicit_identity_mismatch_cannot_be_relabelled_as_other_isoform(tmp_path):
    annotations = row("gene", 1, 21, "g")
    for transcript, start, end in (("t1", 1, 9), ("t2", 13, 21)):
        annotations += row("mRNA", start, end, transcript, parent="g")
        annotations += row("CDS", start, end, "c-" + transcript, parent=transcript, phase="0")
    paths = fixture(tmp_path, ">t1 original conflicting source\naTgCcCtAa\n", annotations, "ATGAAATAANNNATGCCCTAA")
    catalog = MODULE.build_catalog("Species_one", *paths)
    candidates = catalog["loci"][0]["candidates"]
    assert candidates[0]["source_transcript_id"] == "t1"
    assert candidates[0]["quality"]["sequence_mismatch"]
    assert not candidates[0]["quality"]["usable"]
    assert candidates[0]["cds"] == "ATGAAATAA"
    assert candidates[0]["source_cds"] == [{"source_fasta_id": "t1", "header": "t1 original conflicting source",
                                             "cds": "aTgCcCtAa", "sha256": hashlib.sha256(b"aTgCcCtAa").hexdigest(),
                                             "normalized_sha256": hashlib.sha256(b"ATGCCCTAA").hexdigest(),
                                             "sequence_agreement": False, "source_convention": ""}]
    archived = MODULE.write_catalog(catalog, tmp_path / "audit")
    assert json.loads(Path(archived["catalog"]).read_text())["loci"][0]["candidates"][0]["source_cds"] == candidates[0]["source_cds"]
    assert candidates[1]["quality"]["usable"]
    assert catalog["fasta_mapping"][0]["candidate_ids"] == ["Species_one_t1"]


def test_escaped_comma_in_source_parent_is_one_exact_transcript_identity(tmp_path):
    annotations = row("gene", 1, 9, "g") + row("mRNA", 1, 9, "t%2C1", parent="g")
    annotations += row("CDS", 1, 9, "c", parent="t%2C1", phase="0")
    paths = fixture(tmp_path, ">g\nATGAAATAA\n", annotations, "ATGAAATAA")
    catalog = MODULE.build_catalog("Species_one", *paths)
    locus = catalog["loci"][0]
    assert len(locus["candidates"]) == 1
    candidate = locus["candidates"][0]
    assert candidate["source_transcript_id"] == "t,1"
    assert candidate["source_gene_id"] == "g"
    assert candidate["candidate_id"] == "Species_one_t,1"
    assert candidate["quality"]["valid_orf"]
    assert locus["source_baseline_candidate_id"] == candidate["candidate_id"]


def test_source_attribute_parser_distinguishes_gtf_scalars_from_gff3_lists_and_escapes():
    assert MODULE.parse_gff_attributes('gene_id "g"; transcript_id "t%2C1"; Note "literal,semicolon;percent%25";') == {
        "gene_id": ("g",), "transcript_id": ("t%2C1",), "Note": ("literal,semicolon;percent%25",)}
    assert MODULE.parse_gff_attributes('ID=t%2C1;Parent=g1,g%2C2;Alias=t%252C1') == {
        "ID": ("t,1",), "Parent": ("g1", "g,2"), "Alias": ("t%2C1",)}


def test_phase_conflict_is_retained_without_fixing_bases(tmp_path):
    annotations = row("gene", 1, 15, "g") + row("mRNA", 1, 15, "t", parent="g")
    annotations += row("CDS", 1, 4, "c1", parent="t", phase="0")
    annotations += row("CDS", 11, 15, "c2", parent="t", phase="1")
    paths = fixture(tmp_path, ">t\nATGAAATAA\n", annotations, "ATGANNNNNNAATAA")
    candidate = MODULE.build_catalog("Species_one", *paths)["loci"][0]["candidates"][0]
    assert candidate["cds"] == "ATGAAATAA"
    assert candidate["quality"]["phase_conflict"]
    assert not candidate["quality"]["usable"]


@pytest.mark.parametrize("kind", ["gene", "transcript"])
def test_identifier_sanitization_collision_is_rejected(tmp_path, kind):
    annotations = ""
    for number, identifier in enumerate(("a:b", "a_b")):
        gene, transcript = (identifier, f"t{number}") if kind == "gene" else ("g", identifier)
        if kind == "gene" or number == 0:
            annotations += row("gene", 1, 9, gene)
        annotations += row("mRNA", 1, 9, transcript, parent=gene)
        annotations += row("CDS", 1, 9, "c" + str(number), parent=transcript, phase="0")
    paths = fixture(tmp_path, ">unmapped\nATGAAATAA\n", annotations, "ATGAAATAA")
    with pytest.raises(ValueError, match="collide"):
        MODULE.build_catalog("Species_one", *paths)


def test_two_gene_parents_are_not_arbitrarily_collapsed(tmp_path):
    annotations = row("gene", 1, 9, "g1") + row("gene", 1, 9, "g2")
    annotations += row("mRNA", 1, 9, "t", parent="g1,g2")
    annotations += row("CDS", 1, 9, "c", parent="t", phase="0")
    paths = fixture(tmp_path, ">t\nATGAAATAA\n", annotations, "ATGAAATAA")
    with pytest.raises(ValueError, match="Ambiguous|unambiguous"):
        MODULE.build_catalog("Species_one", *paths)


def test_prefixed_effective_ids_are_not_prefixed_again(tmp_path):
    annotations = row("gene", 1, 9, "Species_one_g")
    annotations += row("mRNA", 1, 9, "Species_one_t", parent="Species_one_g")
    annotations += row("CDS", 1, 9, "c", parent="Species_one_t", phase="0")
    paths = fixture(tmp_path, ">Species_one_g\nATGAAATAA\n", annotations, "ATGAAATAA")
    locus = MODULE.build_catalog("Species_one", *paths)["loci"][0]
    assert locus["gene_id"] == "Species_one_g"
    assert locus["candidates"][0]["candidate_id"] == "Species_one_t"
    assert locus["candidates"][0]["source_transcript_id"] == "Species_one_t"


def test_provider_gene_token_and_actual_gff_parent_identity_are_preserved_separately(tmp_path):
    annotations = row("gene", 1, 9, "g").replace("ID=g", "ID=g;Dbxref=GeneID:123")
    annotations += row("mRNA", 1, 9, "t", parent="g")
    annotations += row("CDS", 1, 9, "c", parent="t", phase="0")
    paths = fixture(tmp_path, ">Species_one_GeneID123\nATGAAATAA\n", annotations, "ATGAAATAA")
    locus = MODULE.build_catalog("Species_one", *paths)["loci"][0]
    assert locus["gene_id"] == "Species_one_GeneID123"
    assert locus["source_gene_id"] == "g"
    assert locus["candidates"][0]["source_gene_id"] == "g"
    assert locus["candidates"][0]["gene_token"] == "GeneID123"


def test_uncertain_partial_phase_convention_has_no_invented_formatter_corrected_length(tmp_path):
    annotations = row("gene", 1, 10, "g") + row("mRNA", 1, 10, "t", parent="g")
    annotations += row("CDS", 1, 10, "c", parent="t", phase="1")
    paths = fixture(tmp_path, ">t\nAATGAAATAA\n", annotations, "AATGAAATAA")
    candidate = MODULE.build_catalog("Species_one", *paths)["loci"][0]["candidates"][0]
    assert candidate["cds"] == "AATGAAATAA"
    assert "corrected_cds_length" not in candidate


@pytest.mark.parametrize("organelle_cds_id", ["", "cp_t", "Species_one_cp_g"])
def test_organelle_models_remain_untranslated_provenance_and_never_reenter_nuclear_candidates(tmp_path, organelle_cds_id):
    annotations = row("gene", 1, 9, "g") + row("mRNA", 1, 9, "t", parent="g")
    annotations += row("CDS", 1, 9, "c", parent="t", phase="0")
    annotations += "cp\ts\tregion\t1\t9\t.\t+\t.\tID=cp_region;genome=chloroplast\n"
    annotations += (row("gene", 1, 9, "cp_g") + row("mRNA", 1, 9, "cp_t", parent="cp_g")
                    + row("CDS", 1, 9, "cp_c", parent="cp_t", phase="0")).replace("chr1\t", "cp\t")
    records = ">t\nATGAAATAA\n" + (">" + organelle_cds_id + "\nATGTGATAA\n" if organelle_cds_id else "")
    paths = fixture(tmp_path, records, annotations, "ATGAAATAA")
    paths[2].write_text(paths[2].read_text() + ">cp\nATGTGATAA\n")
    before = [digest(path) for path in paths]
    catalog = MODULE.build_catalog("Species_one", *paths)
    assert [locus["gene_id"] for locus in catalog["loci"]] == ["Species_one_g"]
    assert catalog["summary"]["candidates"] == 1
    assert catalog["summary"]["excluded_candidates"] == 1
    excluded = catalog["excluded_loci"][0]
    assert excluded["source_gene_id"] == "cp_g"
    assert excluded["exclusion_reason"] == "organelle_annotation"
    organelle = excluded["candidates"][0]
    assert organelle["cds"] == "ATGTGATAA"
    assert organelle["protein"] == ""
    assert not organelle["quality"]["usable"]
    assert organelle["quality"]["genetic_code"] is None
    if organelle_cds_id:
        assert catalog["fasta_mapping"][-1]["mapping_status"] == "excluded_organelle"
        assert catalog["fasta_mapping"][-1]["exclusion_reason"] == "organelle_annotation"
    written = MODULE.write_catalog(catalog, tmp_path / "output")
    assert "cp_t" not in Path(written["cds"]).read_text()
    assert "cp_t" not in Path(written["loci"]).read_text()
    assert "excluded_loci" not in json.loads(Path(written["metadata"]).read_text())
    assert "ATGTGATAA" not in Path(written["metadata"]).read_text()
    assert "cp_t" in Path(written["excluded"]).read_text()
    assert [digest(path) for path in paths] == before


@pytest.mark.parametrize("supplied,convention", [("ATGAAATGA", "exact_genomic_cds"),
                                                ("ATGAAA", "omitted_terminal_stop"),
                                                ("ATGAAANNN", "masked_terminal_stop")])
def test_real_genomic_terminal_stop_conventions_are_audited_without_source_edits(tmp_path, supplied, convention):
    annotations = row("gene", 1, 9, "g") + row("mRNA", 1, 9, "t", parent="g")
    annotations += row("CDS", 1, 9, "c", parent="t", phase="0")
    paths = fixture(tmp_path, ">t\n" + supplied + "\n", annotations, "ATGAAATGA")
    before = [digest(path) for path in paths]
    catalog = MODULE.build_catalog("Species_one", *paths)
    candidate = catalog["loci"][0]["candidates"][0]
    assert candidate["cds"] == "ATGAAATGA"
    assert candidate["protein"] == "MK"
    assert candidate["quality"]["usable"]
    assert candidate["quality"]["valid_orf"]
    assert catalog["fasta_mapping"][0]["source_convention"] == convention
    assert candidate["source_conventions"] == [convention]
    assert [digest(path) for path in paths] == before


@pytest.mark.parametrize("genomic,supplied,code", [("ATGAAATGA", "ATGCCCNNN", 1),
                                                ("ATGAAATGA", "ATGAAANNN", 4),
                                                ("ATGNNNTGA", "ATGNNNNNN", 1),
                                                ("ATGTGATGA", "ATGTGANNN", 1)])
def test_stop_conventions_do_not_repair_internal_changes_gaps_or_nonstop_codons(tmp_path, genomic, supplied, code):
    annotations = row("gene", 1, 9, "g") + row("mRNA", 1, 9, "t", parent="g")
    annotations += row("CDS", 1, 9, "c", parent="t", phase="0")
    paths = fixture(tmp_path, ">t\n" + supplied + "\n", annotations, genomic)
    candidate = MODULE.build_catalog("Species_one", *paths, genetic_code=code)["loci"][0]["candidates"][0]
    assert candidate["cds"] == genomic
    assert not candidate["quality"]["usable"]


@pytest.mark.parametrize("partial,unknown_first", [(False, False), (False, True), (True, True)])
def test_unknown_phases_are_inferred_only_for_bound_complete_intact_genomic_orf(tmp_path, partial, unknown_first):
    annotations = row("gene", 1, 15, "g") + row("mRNA", 1, 15, "t", parent="g")
    annotations += row("CDS", 1, 4, "c1", parent="t", phase="." if unknown_first else "0")
    annotations += row("CDS", 11, 14 if partial else 15, "c2", parent="t", phase=".")
    sequence = "ATGAAATA" if partial else "ATGAAATAA"
    paths = fixture(tmp_path, ">t\n" + sequence + "\n", annotations, "ATGANNNNNNAATAA")
    candidate = MODULE.build_catalog("Species_one", *paths)["loci"][0]["candidates"][0]
    assert candidate["quality"]["phase_unknown"]
    assert not candidate["quality"]["phase_conflict"]
    assert candidate["quality"]["phase_inferred"] is (not partial)
    assert candidate["quality"]["usable"] is (not partial)
    if not partial:
        assert candidate["blocks"] == [[0, 4, 0], [10, 15, 2]]
        assert candidate["source_blocks"] == [[0, 4, -1 if unknown_first else 0], [10, 15, -1]]
        assert candidate["quality"]["phase_inference_evidence"] == "complete_genomic_cds_and_bound_source"
    else:
        assert candidate["quality"]["phase_unresolved"]
        assert "source_blocks" not in candidate


@pytest.mark.parametrize("sequence,blocks,reason", [
    ("ATGZZZTAA", [[0, 9, 0]], "invalid_base"),
    ("ATGNNNTAA", [[0, 9, 0]], "ambiguous"),
    ("ATGAAATAA", [[0, 6, 0], [5, 8, 0]], "structure_problem"),
    ("ATGAAATAA", [[0, 8, 0]], "structure_problem"),
])
def test_candidate_validation_withholds_invalid_sequences_and_geometry(sequence, blocks, reason):
    candidate = {"cds": sequence, "blocks": blocks, "strand": "+"}
    quality = MODULE.validate_candidate(candidate)
    assert quality[reason]
    assert not quality["usable"]


@pytest.mark.parametrize("exception,flag", [
    ("annotated_pseudogene", "annotated_pseudogene"),
    ("annotated_translation_exception", "translation_exception"),
    ("annotated_sequence_exception", "sequence_exception"),
])
def test_annotation_exceptions_are_explicit_protected_quality_flags(exception, flag):
    candidate = {"cds": "ATGAAATAA", "blocks": [[0, 9, 0]],
                 "quality": {"annotated_exception": exception}}
    quality = MODULE.validate_candidate(candidate)
    assert quality[flag]
    assert not quality["usable"]
