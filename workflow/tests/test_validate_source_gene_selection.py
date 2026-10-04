import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
import validate_longest_cds_selection as longest
from validate_source_gene_selection import SourceGeneSelection


@pytest.fixture
def annotation(tmp_path):
    path = tmp_path / "source.gff"
    # Children intentionally precede their parents; ancestry must be order-free.
    path.write_text(
        "chr1\ts\tmRNA\t1\t20\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr1\ts\tmRNA\t1\t20\t.\t+\t.\tID=t2;Parent=g1\n"
        "chr1\ts\tgene\t1\t20\t.\t+\t.\tID=g1\n"
        "chr1\ts\tgene\t30\t50\t.\t+\t.\tID=g2\n"
        "chr1\ts\tmRNA\t30\t50\t.\t+\t.\tID=t3;Parent=g2\n")
    return path


def test_independent_source_check_rejects_retained_isoforms(annotation):
    check = SourceGeneSelection(annotation)
    check.observe("t1", "Species_one_gene1", True)
    check.observe("t2", "Species_one_transcript2", True)
    with pytest.raises(ValueError, match="retained_isoform_gene_groups=1"):
        check.validate()


def test_independent_source_check_rejects_distinct_gene_collapse(annotation):
    check = SourceGeneSelection(annotation)
    check.observe("t1", "Species_one_collapsed", True)
    check.observe("t3", "Species_one_collapsed", False)
    with pytest.raises(ValueError, match="distinct_source_gene_merges=1"):
        check.validate()


def test_independent_source_check_only_counts_retained_representatives(annotation):
    check = SourceGeneSelection(annotation)
    check.observe("t1", "one", True)
    check.observe("t2", "unresolved", False)
    check.observe("no_source_identity", "unresolved", True)
    stats = check.validate()
    assert stats["source_gene_resolved_records"] == 2
    assert stats["source_gene_unresolved_records"] == 1
    assert stats["source_gene_check"] == "partial_explicit_parent"
    assert not stats["source_gene_complete"]
    assert stats["source_gene_unresolved_selected_records"] == 1


def test_independent_source_check_accepts_correct_one_gene_selection(annotation):
    check = SourceGeneSelection(annotation)
    check.observe("t1", "one", True)
    check.observe("t2", "one", True)
    check.observe("t3", "two", True)
    assert check.validate()["source_gene_retained_isoform_groups"] == 0


def test_native_longest_validation_does_not_trust_its_own_wrong_grouping(annotation, monkeypatch):
    task = dict(provider="direct", species_key="Species_one", species_prefix="Species_one",
                gff_path=annotation, cds_path=annotation)
    monkeypatch.setattr(longest.formatter, "prepare_cds_identifier_task", lambda task: task)
    monkeypatch.setattr(longest.formatter, "iter_normalised_cds_records", lambda task: iter([("t1", "ATGAAA", {}), ("t2", "ATGAAATTT", {})]))
    monkeypatch.setattr(longest.formatter, "build_formatted_cds_id", lambda task, header: header)
    monkeypatch.setattr(longest.formatter, "build_gene_aggregate_id", lambda task, header, tid: "Species_one_" + tid)
    with pytest.raises(ValueError, match="retained_isoform_gene_groups=1"):
        longest.collect_expected_longest_records(task)


def test_source_check_does_not_guess_missing_parent_from_suffix(tmp_path):
    path = tmp_path / "source.gff"
    path.write_text("chr1\ts\tCDS\t1\t9\t.\t+\t0\tID=g1.2\n")
    check = SourceGeneSelection(path)
    check.observe("g1.2", "gene1", True)
    assert check.validate()["source_gene_unresolved_records"] == 1
    assert check.validate()["source_gene_check"] == "unresolved"
    assert not check.validate()["source_gene_complete"]


def test_missing_gene_features_preserve_explicit_parent_ownership(tmp_path):
    path = tmp_path / "source.gff"
    path.write_text(
        "chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=gene%3Ag1\n"
        "chr1\ts\tmRNA\t30\t38\t.\t+\t.\tID=t2;Parent=gene%3Ag1\n"
        "chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=t3;Parent=g2\n")
    check = SourceGeneSelection(path)
    check.observe("t1", "g1", True)
    check.observe("t2", "g1", False)
    check.observe("t3", "g2", True)
    stats = check.validate()
    assert stats["source_gene_complete"]
    assert stats["source_gene_resolved_records"] == 3
    split = SourceGeneSelection(path)
    split.observe("t1", "first", True)
    split.observe("t2", "second", True)
    with pytest.raises(ValueError, match="retained_isoform_gene_groups=1"):
        split.validate()
    merge = SourceGeneSelection(path)
    merge.observe("t1", "merged", True)
    merge.observe("t3", "merged", False)
    with pytest.raises(ValueError, match="distinct_source_gene_merges=1"):
        merge.validate()


@pytest.mark.parametrize("axis", [("chr2", "+"), ("chr1", "-")])
def test_missing_parent_conflicting_axes_fail_without_guessing(tmp_path, axis):
    path = tmp_path / "source.gff"
    path.write_text("chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n"
                    f"{axis[0]}\ts\tmRNA\t1\t9\t.\t{axis[1]}\t.\tID=t2;Parent=g1\n")
    with pytest.raises(ValueError, match="Conflicting axes"):
        SourceGeneSelection(path)
