"""Advisory model evidence, optional tracks, legacy receipts and data preservation."""
import copy
import json
from argparse import Namespace
from pathlib import Path

import pysam
import pytest
from Bio.Seq import Seq

from workflow.support import rescue_model_evidence as evidence
from workflow.support.input_generation_array_state import digest
from workflow.support.rescue_model_quality import model_quality


@pytest.mark.parametrize("start", ["ATG", "TTG", "CTG", "AAA"])
def test_start_table_does_not_establish_translation(start):
    model = {"sequence": start + "AAATAA", "query_start": 0, "query_end": 99, "query_length": 100,
             "evidence": {"donor": "Plant_one", "expected_strand": "+"}, "strand": "+"}
    quality = model_quality(model)
    assert ("alternative_start" in quality["flags"]) == (start in {"TTG", "CTG"})
    assert quality["translation_initiation"] == "not_established"
    assert quality["native_terminal_completeness"] == "not_established"
    assert "donor_c_terminus_unaligned" in quality["flags"]
    assert quality["terminal_alignment"] == {"n_aligned": True, "c_aligned": False, "query_length": 100}


def test_same_species_isoforms_are_one_donor_and_inversion_is_advisory():
    model = {"sequence": "ATGAAATAA", "strand": "-", "support": [
        {"donor": "Plant_one", "query": "iso1", "expected_strand": "+"},
        {"donor": "Plant_one", "query": "iso2", "expected_strand": "-"}]}
    quality = model_quality(model)
    assert quality["donor_species"] == ["Plant_one"]
    assert "single_donor_species" in quality["flags"]
    assert "expected_strand_conflict" in quality["flags"]
    assert quality["terminal_alignment"]["n_aligned"] is None


@pytest.mark.parametrize("begin,end,flags", [(0, 100, set()), (5, 100, {"donor_n_terminus_unaligned"}),
                                          (0, 95, {"donor_c_terminus_unaligned"}),
                                          (2, 97, {"donor_n_terminus_unaligned", "donor_c_terminus_unaligned"})])
def test_95pct_coverage_is_not_native_terminal_completeness(begin, end, flags):
    quality = model_quality({"sequence": "ATGAAATAA", "coverage": .95, "query_start": begin,
                             "query_end": end, "query_length": 100})
    if end - begin > 95:
        flags = flags | {"donor_internal_unaligned_query"}
    assert set(quality["flags"]) == flags
    assert quality["native_terminal_completeness"] == "not_established"


def test_end_residues_aligned_do_not_hide_a_large_internal_deletion():
    from Bio.Align import PairwiseAligner
    donor, fragment = "M" + "W" * 100 + "K", "MK"
    aligned = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                             open_gap_score=-5, extend_gap_score=-1).align(donor, fragment)[0]
    paired = int(sum(b - a for a, b in aligned.aligned[0]))
    quality = model_quality({"query_start": int(aligned.aligned[0][0][0]),
                             "query_end": int(aligned.aligned[0][-1][1]), "query_length": len(donor),
                             "coverage": paired / len(donor)})
    assert quality["terminal_alignment"]["n_aligned"] and quality["terminal_alignment"]["c_aligned"]
    assert quality["query_alignment"]["aligned_query_fraction"] == pytest.approx(2 / 102)
    assert quality["query_alignment"]["internal_unaligned_query_fraction"] == pytest.approx(100 / 102)
    assert "donor_internal_unaligned_query" in quality["flags"]
    assert quality["native_terminal_completeness"] == "not_established"


@pytest.mark.parametrize("text", ["", "{}", "[{}", "[{},]", "[{}]x", "[null]", "[{} {}]"])
def test_stream_rejects_incomplete_or_malformed_models(tmp_path, text):
    path = tmp_path / "models.json"
    path.write_text(text)
    with pytest.raises(ValueError):
        list(evidence.json_array(path, chunk_size=3))


def test_stream_handles_arbitrary_boundaries(tmp_path):
    rows = [{"message": 'with " quotes and UTF-8 植物', "x": [1, 2, 3]}, {}]
    path = tmp_path / "models.json"
    path.write_text(json.dumps(rows, ensure_ascii=False))
    assert list(evidence.json_array(path, chunk_size=1)) == rows
    path.write_text("[ ]")
    assert list(evidence.json_array(path, chunk_size=1)) == []


def test_junctions_do_not_invent_a_single_rna_transcript():
    model = {"seqid": "chr", "strand": "+", "cds": [[10, 20, 0], [30, 40, 2], [50, 60, 1]]}
    tracks = {("chr", "+"): [(0, 40, [(0, 20), (30, 40)], "left", "leaf"),
                                (30, 70, [(30, 40), (50, 70)], "right", "leaf")]}
    result = evidence.rna_evidence(model, {"rna_junctions": [{}], "rna_transcripts": [{}]},
                                   {("chr", "+"): {(20, 30), (40, 50)}}, evidence.IntervalIndex(tracks))
    assert result["junctions_supported"] == 2
    assert result["exon_chain_transcripts"] == []
    assert result["exon_chain_independence_groups"] == []


def test_technical_runs_and_unknown_strand_not_independent_confirmation():
    model = {"seqid": "chr", "strand": "+", "cds": [[10, 20, 0], [30, 40, 2]]}
    tracks = {("chr", "+"): [(0, 50, [(0, 20), (30, 50)], str(n), "one_individual_leaf") for n in range(3)],
              ("chr", "."): [(0, 50, [(0, 20), (30, 50)], "unknown", "one_individual_leaf")]}
    result = evidence.rna_evidence(model, {"rna_transcripts": [{}]}, {}, evidence.IntervalIndex(tracks))
    assert len(result["exon_chain_transcripts"]) == 3
    assert result["exon_chain_independence_groups"] == ["one_individual_leaf"]
    assert result["unstranded_exon_chain_transcripts"] == ["unknown"]
    assert result["junctions_supported"] is None
    assert result["translation_start_confirmed"] is False


def test_repeat_overlap_union_does_not_double_count_or_reject_gene():
    model = {"seqid": "chr", "cds": [[0, 100, 0]]}
    tracks = evidence.IntervalIndex({"chr": [(0, 60, "LINE/L1"), (30, 100, "LINE/L1"), (0, 100, "Simple_repeat")]})
    result = evidence.repeat_evidence(model, {"repeats": [{}]}, tracks)
    assert result["masked_fraction"] == result["te_annotated_fraction"] == 1
    assert result["classes"] == ["LINE/L1", "Simple_repeat"]
    assert "decision" not in result


def test_literal_identifiers_in_rna_are_not_split_or_reinterpreted():
    assert evidence.transcript_parents("Parent=iso%2C1,iso%3B2") == ["iso,1", "iso;2"]
    assert evidence.transcript_parents('transcript_id "iso,1=a";') == ["iso,1=a"]


def test_escaped_gtf_identifier_does_not_merge_distinct_rna_paths():
    assert evidence.transcript_parents(r'transcript_id "iso\"1;=a"; gene_id "g";') == ['iso"1;=a']
    assert evidence.transcript_parents(r'transcript_id "iso\\1";') == [r"iso\1"]


@pytest.mark.parametrize("text", ['Parent=t1;Parent=t2', 'transcript_id "t1"; transcript_id "t2";'])
def test_ambiguous_transcript_identity_fails(text):
    with pytest.raises(ValueError, match="Duplicate"):
        evidence.transcript_parents(text)


def test_escaped_transcripts_cannot_fabricate_a_complete_spliced_path(tmp_path):
    track = tmp_path / "rna.gtf"
    track.write_text('chr\ts\texon\t1\t20\t.\t+\t.\ttranscript_id "iso\\\"1";\n'
                     'chr\ts\texon\t31\t50\t.\t+\t.\ttranscript_id "iso\\\"2";\n')
    spec = {"rna_transcripts": [{"path": str(track), "independence_group": "leaf"}]}
    junctions, transcripts, _ = evidence.load_tracks(spec, {"chr": 50})
    result = evidence.rna_evidence({"seqid": "chr", "strand": "+", "cds": [[10, 20, 0], [30, 40, 2]]},
                                   spec, junctions, transcripts)
    assert result["exon_chain_transcripts"] == []


def fixture(tmp_path, strand="+"):
    root = tmp_path / "rescue"
    worker = root / "rescued" / "Plant_example"
    worker.mkdir(parents=True)
    seq = "TTGAAATAA"
    genome = tmp_path / "source.fa"
    genome.write_text(">chr\n" + (str(Seq(seq).reverse_complement()) if strand == "-" else seq) + "\n")
    model = {"model_id": "Plant_example_ggrescue_a", "seqid": "chr", "strand": strand, "cds": [[0, 9, 0]],
             "sequence": seq, "status": "accepted", "query_start": 0, "query_end": 2, "query_length": 2,
             "support": [{"donor": "Plant_one", "expected_strand": strand}]}
    plan = {"request": {"sources": {"Plant_example": {"genome": str(genome), "genetic_code": 1}},
                        "files": {str(genome): digest(genome)}}}
    (root / "plan.json").write_text(json.dumps(plan))
    (worker / "models.json").write_text(json.dumps([model, {**model, "status": "duplicate_support"}]))
    (worker / "receipt.json").write_text(json.dumps({"key": {"plan": digest(root / "plan.json")},
                                                   "files": {"models.json": digest(worker / "models.json")}}))
    args = Namespace(rescue_output=root, species="Plant_example", output=tmp_path / "audit", evidence_manifest=None)
    return args, genome


@pytest.mark.parametrize("strand", ["+", "-"])
def test_audit_legacy_models_missing_inputs_and_preserved_sources(tmp_path, strand):
    args, genome = fixture(tmp_path, strand)
    before = {p: p.read_bytes() for p in args.rescue_output.rglob("*") if p.is_file()}
    summary = evidence.audit(args)
    assert summary["accepted_models"] == 1
    assert summary["sequence_changes"] == 0
    row = json.loads((args.output / "evidence.json").read_text())[0]
    assert all(row[k]["status"] == "not_provided" for k in ("rna", "repeat", "dna"))
    assert row["decision"] == "unchanged" and "alternative_start" in row["review_flags"]
    assert before == {p: p.read_bytes() for p in before}
    assert not Path(str(genome) + ".fai").exists()
    assert evidence.audit(args) == summary
    (args.output / "evidence.json").write_text("corrupt")
    assert evidence.audit(args) == summary


def test_plan_changed_after_read_cannot_label_old_parsed_values_with_new_hash(tmp_path, monkeypatch):
    args, _ = fixture(tmp_path)
    path = args.rescue_output / "plan.json"
    old = path.read_bytes()
    original = json.loads
    changed = []
    def changing_loads(value, *positional, **kwargs):
        parsed = original(value, *positional, **kwargs)
        if not changed and (value.encode() if isinstance(value, str) else value) == old:
            changed.append(True)
            new = copy.deepcopy(parsed)
            new["request"]["sources"][args.species]["genetic_code"] = 2
            path.write_text(json.dumps(new))
            receipt = args.rescue_output / "rescued" / args.species / "receipt.json"
            saved = original(receipt.read_text())
            saved["key"]["plan"] = digest(path)
            receipt.write_text(json.dumps(saved))
        return parsed
    monkeypatch.setattr(json, "loads", changing_loads)
    with pytest.raises(ValueError):
        evidence.audit(args)
    assert changed and not (args.output / "receipt.json").exists()


@pytest.mark.parametrize("kind", ["models", "genome"])
def test_frozen_file_changed_after_initial_verification_is_not_rebound(tmp_path, monkeypatch, kind):
    args, genome = fixture(tmp_path)
    path = args.rescue_output / "rescued" / args.species / "models.json" if kind == "models" else genome
    original = evidence.digest
    changed = []
    def changing_digest(value):
        result = original(value)
        if Path(value) == path and not changed:
            changed.append(True)
            path.write_text("[]" if kind == "models" else ">chr\nCTGAAATAA\n")
        return result
    monkeypatch.setattr(evidence, "digest", changing_digest)
    with pytest.raises(ValueError):
        evidence.audit(args)
    assert changed and not (args.output / "receipt.json").exists()


def test_manifest_changed_after_parse_cannot_hide_supplied_rna(tmp_path, monkeypatch):
    args, genome = fixture(tmp_path)
    track = tmp_path / "rna.gtf"
    track.write_text('chr\ts\texon\t1\t9\t.\t+\t.\ttranscript_id "t";\n')
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps({"schema_version": 1, "species": {args.species: {
        "reference_genome_sha256": digest(genome)}}}))
    args.evidence_manifest = path
    old = path.read_bytes()
    original = json.loads
    changed = []
    def changing_loads(value, *positional, **kwargs):
        parsed = original(value, *positional, **kwargs)
        if not changed and (value.encode() if isinstance(value, str) else value) == old:
            changed.append(True)
            new = copy.deepcopy(parsed)
            new["species"][args.species]["rna_transcripts"] = [{"path": str(track), "format": "exon_gff_gtf",
                                                              "independence_group": "leaf"}]
            path.write_text(json.dumps(new))
        return parsed
    monkeypatch.setattr(json, "loads", changing_loads)
    with pytest.raises(ValueError, match="Input changed while loading"):
        evidence.audit(args)
    assert changed and not (args.output / "receipt.json").exists()


@pytest.mark.parametrize("strand", ["+", "-"])
def test_dna_split_initiator_follows_spliced_transcription_order(tmp_path, strand):
    coding, genomic = "TTGAAATAA", "TGTAGTGAAATAA"
    blocks = [[0, 1, 0], [5, 13, 2]]
    if strand == "-":
        genomic = str(Seq(genomic).reverse_complement())
        blocks = [[12, 13, 0], [0, 8, 2]]
    genome_path, bam_path = tmp_path / "genome.fa", tmp_path / "reads.bam"
    genome_path.write_text(">chr\n" + genomic + "\n")
    pysam.faidx(str(genome_path))
    with pysam.AlignmentFile(str(bam_path), "wb", header={"SQ": [{"SN": "chr", "LN": 13}]}) as bam:
        for number, (mapq, baseq) in enumerate([(60, "I"), (255, "I"), (10, "I"), (60, "!")]):
            read = pysam.AlignedSegment()
            read.query_name, read.query_sequence = str(number), genomic
            read.flag, read.reference_id, read.reference_start = 0, 0, 0
            read.mapping_quality, read.cigar = mapq, [(0, 13)]
            read.query_qualities = pysam.qualitystring_to_array(baseq * 13)
            bam.write(read)
    pysam.index(str(bam_path))
    model = {"seqid": "chr", "strand": strand, "cds": blocks, "sequence": coding}
    spec = {"dna": {"min_mapq": 20, "min_baseq": 20}}
    with pysam.FastaFile(str(genome_path)) as genome, pysam.AlignmentFile(str(bam_path), "rb") as bam:
        evidence.check_model(model, genome)
        result = evidence.dna_evidence(model, spec, bam, genome)
    assert result["reference_identity"] == "not_verified"
    assert result["hq_depth_min"] == result["hq_depth_median"] == 1
    assert [row["genomic_position_1based"] for row in result["start_base_support"]] == ([1, 6, 7] if strand == "+" else [13, 8, 7])
    assert [row["matching"] for row in result["start_base_support"]] == [1, 1, 1]


@pytest.mark.parametrize("fault", ["model", "receipt", "genome", "inside", "ancestor"])
def test_reject_changed_sources_and_dangerous_output(tmp_path, fault):
    args, genome = fixture(tmp_path)
    worker = args.rescue_output / "rescued" / args.species
    if fault == "model":
        (worker / "models.json").write_text("[]")
    elif fault == "receipt":
        receipt = json.loads((worker / "receipt.json").read_text())
        receipt["key"]["plan"] = "foreign"
        (worker / "receipt.json").write_text(json.dumps(receipt))
    elif fault == "genome":
        genome.write_text(">chr\nATGAAATAA\n")
    elif fault == "inside":
        args.output = worker / "audit"
    else:
        args.output = tmp_path
    with pytest.raises(ValueError):
        evidence.audit(args)


def test_identical_coding_sequences_at_distinct_loci_are_preserved(tmp_path):
    args, genome = fixture(tmp_path)
    genome.write_text(">chr\nTTGAAATAACCCCTTGAAATAA\n")
    worker = args.rescue_output / "rescued" / args.species
    model = json.loads((worker / "models.json").read_text())[0]
    other = {**model, "model_id": "Plant_example_ggrescue_b", "cds": [[13, 22, 0]]}
    (worker / "models.json").write_text(json.dumps([model, other]))
    plan_path = args.rescue_output / "plan.json"
    plan = json.loads(plan_path.read_text())
    plan["request"]["files"][str(genome)] = digest(genome)
    plan_path.write_text(json.dumps(plan))
    (worker / "receipt.json").write_text(json.dumps({"key": {"plan": digest(plan_path)},
                                                   "files": {"models.json": digest(worker / "models.json")}}))
    assert evidence.audit(args)["accepted_models"] == 2
    rows = json.loads((args.output / "evidence.json").read_text())
    assert len({row["model_id"] for row in rows}) == 2
    assert len({row["cds_sha256"] for row in rows}) == 1


@pytest.mark.parametrize("strand", ["+", "-"])
def test_real_bam_start_bases_repeat_and_rna_are_advisory(tmp_path, strand):
    args, genome = fixture(tmp_path, strand)
    bam_path = tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header={"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr", "LN": 9}]}) as bam:
        for number, flag in enumerate([0, 0, 256, 2048, 1024, 512, 0]):
            read = pysam.AlignedSegment()
            read.query_name = str(number)
            read.query_sequence = genome.read_text().splitlines()[1]
            read.flag, read.reference_id, read.reference_start = flag, 0, 0
            read.mapping_quality, read.cigar = 60, [(0, 9)]
            if number == 6:
                read.mapping_quality = 255  # unavailable, rather than very high quality
            read.query_qualities = pysam.qualitystring_to_array("I" * 9)
            bam.write(read)
    pysam.index(str(bam_path))
    hints = tmp_path / "hints.gff"
    hints.write_text("chr\tProtHint\tintron\t4\t6\t1\t+\t.\tsrc=P\n")
    transcripts = tmp_path / "rna.gtf"
    transcripts.write_text(f'chr\tStringTie\texon\t1\t9\t.\t{strand}\t.\ttranscript_id "t1";\n')
    repeats = tmp_path / "repeats.out"
    repeats.write_text("100 0 0 0 chr 1 9 (0) + repeat LINE/L1 1 9 (0) 1\n")
    manifest = {"schema_version": 1, "species": {args.species: {"reference_genome_sha256": digest(genome),
        "rna_junctions": [{"path": str(hints), "format": "braker_hints", "independence_group": "leaf"}],
        "rna_transcripts": [{"path": str(transcripts), "format": "exon_gff_gtf", "independence_group": "leaf"}],
        "repeats": [{"path": str(repeats), "format": "repeatmasker_out"}],
        "dna": {"path": str(bam_path), "index": str(bam_path) + ".bai", "format": "bam", "min_mapq": 20, "min_baseq": 20}}}}
    args.evidence_manifest = tmp_path / "manifest.json"
    args.evidence_manifest.write_text(json.dumps(manifest))
    evidence.audit(args)
    row = json.loads((args.output / "evidence.json").read_text())[0]
    assert row["dna"]["hq_depth_min"] == row["dna"]["hq_depth_median"] == 2
    assert [base["matching"] for base in row["dna"]["start_base_support"]] == [2, 2, 2]
    assert row["rna"]["exon_chain_transcripts"] == ["0:t1"]
    assert row["repeat"]["te_annotated_fraction"] == 1
    assert row["decision"] == "unchanged"
    assert "te_overlap_ge_50pct" in row["review_flags"]
    bad = copy.deepcopy(manifest)
    bad["species"][args.species]["reference_genome_sha256"] = "foreign"
    args.evidence_manifest.write_text(json.dumps(bad))
    with pytest.raises(ValueError, match="reference differs"):
        evidence.audit(args)


@pytest.mark.parametrize("fault", [None, "no_reference", "changed_bases", "changed_length",
                                  "missing_target", "wrong_original_header", "duplicate_fasta", "blank_header", "empty_record"])
def test_pre_filter_bam_requires_verified_original_reference(tmp_path, fault):
    args, genome = fixture(tmp_path)
    original = tmp_path / "original.fa"
    original.write_text(">chr\n" + ("CTGAAATAA" if fault == "changed_bases" else "TTGAAATAA")
                        + "\n>excluded\n" + ("ACGTA" if fault == "wrong_original_header" else "ACGT") + "\n")
    if fault == "duplicate_fasta":
        with original.open("a") as handle:
            handle.write(">excluded\nACGT\n")
    if fault in {"blank_header", "empty_record"}:
        with original.open("a") as handle:
            handle.write(">\n" if fault == "blank_header" else ">empty\n")
    refs = [{"SN": "chr", "LN": 10 if fault == "changed_length" else 9}, {"SN": "excluded", "LN": 4}]
    if fault == "missing_target":
        refs = refs[1:]
    bam_path = tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header={"HD": {"SO": "coordinate"}, "SQ": refs}):
        pass
    pysam.index(str(bam_path))
    dna = {"path": str(bam_path), "index": str(bam_path) + ".bai", "format": "bam", "min_mapq": 20, "min_baseq": 20}
    if fault != "no_reference":
        dna["reference_genome"] = str(original)
    args.evidence_manifest = tmp_path / "manifest.json"
    args.evidence_manifest.write_text(json.dumps({"schema_version": 1, "species": {args.species: {
        "reference_genome_sha256": digest(genome), "dna": dna}}}))
    if fault:
        with pytest.raises(ValueError):
            evidence.audit(args)
    else:
        summary = evidence.audit(args)
        assert summary["dna_reference"]["extra_bam_contigs"] == 1
        assert summary["dna_reference"]["reference_genome_sha256"] == digest(original)
        assert summary["accepted_models"] == 1 and summary["sequence_changes"] == 0
        row = json.loads((args.output / "evidence.json").read_text())[0]
        assert row["dna"]["hq_depth_min"] == 0
        assert row["decision"] == "unchanged"
        assert not Path(str(original) + ".fai").exists()


def test_runtime_cli_help_has_no_writes(tmp_path):
    from workflow.tests.test_support_script_help_smoke import test_support_script_help_smoke
    test_support_script_help_smoke("rescue_model_evidence.py", tmp_path)
