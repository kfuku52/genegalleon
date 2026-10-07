"""GeMoMa uses the same genomic terminal repair and bounded donor QC.

Only the external Java process is simulated.  Genomic reconstruction, protein
alignment, strict CDS validation and terminal completion are real.
"""
import copy
from pathlib import Path

import pytest
from Bio.Seq import Seq

from workflow.support import rescue_gene_models as rescue


class Genome:
    def __init__(self, sequence):
        self.sequence = sequence

    def get_reference_length(self, seqid):
        assert seqid == "chr1"
        return len(self.sequence)

    def fetch(self, seqid, start, end):
        assert seqid == "chr1" and 0 <= start <= end <= len(self.sequence)
        return self.sequence[start:end]


def predict(tmp_path, monkeypatch, coding, protein, interval, *, strand="+", parameters=None):
    donor, target, query = "Species_donor", "Species_target", "Species_donor_t1"
    prepared = tmp_path / "prepared" / donor
    prepared.mkdir(parents=True)
    rescue.write_tsv(prepared / "genes.id_map.tsv", ("original_id", "jcvi_id", "locus_id", "status"),
                     [(query, query, "donor_gene", "selected")])
    (prepared / "genes.anchor_admission.json").write_text('{"records": []}')
    (prepared / "genes.pep").write_text(f">{query}\n{protein}\n")
    gff, jar, java = [tmp_path / name for name in ("donor.gff3", "test.jar", "java")]
    gff.write_text(f"chr1\tsource\tmRNA\t1\t{len(coding)}\t.\t+\t.\tID={query};Parent=donor_gene\n")
    jar.write_text("simulated GeMoMa jar")
    java.write_text("simulated Java executable")
    plan = {"request": {"gemoma_jar": str(jar), "gemoma_java": str(java),
                        "files": {str(path): rescue.digest(path) for path in (jar, java)},
                        "parameters": {"minimum_coverage": .6, "minimum_identity": .5, "max_intron": 100,
                                       "terminal_max_extension": 30, "terminal_max_unaligned_c_overhang": 2,
                                       **(parameters or {})},
                        "sources": {donor: {"gff": str(gff), "genome": str(tmp_path / "donor.fa")}}}}
    begin, end = interval
    if strand == "-":
        begin, end = len(coding) - end, len(coding) - begin
    region = {"target": target, "donor": donor, "query": query, "id": "query1", "seqid": "chr1",
              "expected_strand": strand, "start": 0, "end": len(coding),
              "expected_start": 0, "expected_end": len(coding)}
    genome = Genome(coding if strand == "+" else str(Seq(coding).reverse_complement()))
    before = (genome.sequence, copy.deepcopy(region), gff.read_bytes(), (prepared / "genes.pep").read_bytes())
    monkeypatch.setattr(rescue, "verify_sources", lambda *_: None)

    def fake_java(command, directory, label):
        assert label == region["id"] and "GeMoMaPipeline" in command
        output = Path(next(str(item).split("=", 1)[1] for item in command if str(item).startswith("outdir=")))
        output.mkdir()
        (output / "final_annotation.gff").write_text(
            f"chr1\tGeMoMa\tmRNA\t{begin + 1}\t{end}\t.\t{strand}\t.\tID=prediction\n"
            f"chr1\tGeMoMa\tCDS\t{begin + 1}\t{end}\t.\t{strand}\t0\tParent=prediction\n")

    monkeypatch.setattr(rescue, "run", fake_java)
    rows = []
    rescue.refine_gemoma(tmp_path, tmp_path, plan, {"species": target, "genetic_code": 1},
                         [region], genome, rows, 1)
    assert before == (genome.sequence, region, gff.read_bytes(), (prepared / "genes.pep").read_bytes())
    assert len(rows) == 1
    return rows[0]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_gemoma_repairs_both_termini_with_real_source_and_short_c_overhang(tmp_path, monkeypatch, strand):
    model = predict(tmp_path, monkeypatch, "ATGAAACCCGGGTAA", "MKP", (3, 9), strand=strand)
    assert model["sequence"] == "ATGAAACCCGGGTAA" and model["cds"] == [[0, 15, 0]]
    assert model["problems"] == []
    assert model["terminal_completion"]["status"] == "completed"
    assert model["partial_evidence"]["completed_from_partial"] and not model["partial_evidence"]["partial"]
    assert model["coverage"] == model["identity"] == 1.0
    assert model["query_start"] == 0 and model["query_end"] == model["query_length"] == 3
    raw = model["raw_prediction"]
    assert raw["evidence"] == model["evidence"]
    assert raw["cds"] == ([[3, 9, 0]] if strand == "+" else [[6, 12, 0]])
    assert raw["coverage"] == pytest.approx(2 / 3)
    assert raw["query_start"] == 1 and raw["query_end"] == 3
    assert model["terminal_completion"]["attempts"][0]["alignment"]["unaligned_c_overhang_residues"] == 1


def test_gemoma_keeps_true_partial_as_ineligible_evidence(tmp_path, monkeypatch):
    model = predict(tmp_path, monkeypatch, "CCCAAACCCTAA", "MKP", (3, 9))
    assert model["sequence"] == "AAACCCTAA" and model["cds"] == [[3, 12, 0]]
    assert "missing_start" in model["problems"]
    assert model["terminal_completion"]["status"] == "rejected"
    assert model["partial_evidence"]["partial"] and not model["partial_evidence"]["representative_eligible"]
    assert model["raw_prediction"]["cds"] == [[3, 9, 0]]


def test_gemoma_never_uses_terminal_extension_to_conceal_internal_stop(tmp_path, monkeypatch):
    model = predict(tmp_path, monkeypatch, "ATGAAATAACCCTAA", "MKP", (3, 12))
    assert model["sequence"] == "AAATAACCCTAA" and model["cds"] == [[3, 15, 0]]
    assert {"missing_start", "internal_stop"} <= set(model["problems"])
    assert model["terminal_completion"]["reasons"] == ["non_terminal_failure_preserved"]
    assert model["partial_evidence"]["partial"]


def test_gemoma_donor_alignment_budget_is_audited_before_quadratic_dp(tmp_path, monkeypatch):
    import Bio.Align

    def forbidden_alignment(*args, **kwargs):
        pytest.fail("A >25M-cell donor alignment must be withheld before allocating the aligner")

    monkeypatch.setattr(Bio.Align, "PairwiseAligner", forbidden_alignment)
    coding = "ATG" + "AAA" * 5000 + "TAA"
    model = predict(tmp_path, monkeypatch, coding, "M" + "K" * 5000, (3, len(coding) - 3))
    assert "donor_alignment_budget_exceeded" in model["problems"]
    assert model["partial_evidence"]["partial"] and not model["partial_evidence"]["representative_eligible"]
    assert model["sequence"] == coding[3:]
    assert "raw_prediction" in model
