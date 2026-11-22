import csv
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from cdskit.homology import parse_hit, write_selection_report

SUPPORT = Path(__file__).resolve().parents[1] / "support"


@pytest.mark.parametrize(("baseline_target", "baseline_end", "baseline_score", "baseline_length", "retains"), [
    ("reference", 98, 386, 211, True),
    ("reference", 98, 254, 211, True),
    ("reference", 62, 254, 211, False),
    ("reference", 59, 386, 211, False),
    ("reference", 98, 253, 211, False),
    ("other", 98, 386, 211, False),
    ("reference", 98, 386, 110, False),
])
def test_representative_guard_uses_dominating_hit_to_same_target(
    tmp_path, monkeypatch, baseline_target, baseline_end, baseline_score, baseline_length, retains,
):
    monkeypatch.syspath_prepend(str(SUPPORT))
    from transcriptome_homology import retain_supported_longer_orfs

    baseline = SeqRecord(Seq("AAA" * baseline_length), id="long", description="long original ORF")
    selected = SeqRecord(Seq("CCC" * 110), id="short")
    unchanged = SeqRecord(Seq("GCC" * 120), id="unchanged")
    candidates = tmp_path / "candidates.fa"
    output, report, raw = (tmp_path / name for name in ("selected.fa", "selection.tsv", "diamond.tsv"))
    SeqIO.write((baseline, selected, unchanged), candidates, "fasta")
    SeqIO.write((selected, unchanged), output, "fasta")
    selected_hit = "orf2\treference\t1e-12\t254\t10\t100\t110\t1\t62\t100"
    write_selection_report(report, [
        ("group", baseline, selected, "supported_alternative", None, parse_hit(selected_hit.split("\t"))[1]),
        ("other_group", unchanged, unchanged, "single_candidate", None, None),
    ])
    raw.write_text(
        f"orf1\t{baseline_target}\t1e-20\t{baseline_score}\t1\t50\t{baseline_length}\t1\t{baseline_end}\t100\n"
        f"{selected_hit}\n"
    )
    retain_supported_longer_orfs(raw, [("long",), ("short",)], candidates, output, report)

    expected = baseline if retains else selected
    assert [(record.id, record.seq) for record in SeqIO.parse(output, "fasta")] == [
        (expected.id, expected.seq), (unchanged.id, unchanged.seq)
    ]
    with report.open(newline="") as stream:
        row, other = csv.DictReader(stream, delimiter="\t")
    assert row["selected_id"] == expected.id
    assert row["reason"] == ("retain_dominant_baseline_hit" if retains else "supported_alternative")
    assert other["selected_id"] == unchanged.id
    if retains:
        assert row["selected_description"] == baseline.description
        assert float(row["baseline_target_coverage"]) == float(row["selected_target_coverage"]) == .98
        assert float(row["selected_query_coverage"]) < .60
