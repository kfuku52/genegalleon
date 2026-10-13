"""Synthetic counterparts of measured placement holds, without genomic inputs."""
import copy

from workflow.support.rescue_additional_candidates import reassess_unanchored_models
from workflow.support.rescue_coding_paths import OwnershipIndex
from workflow.support.rescue_gene_models import consolidate


def model(donor, query, start=0, length=99, *, seqid="s", strand="+", cds=None, problems=None):
    blocks = cds or [[start, start + length, 0]]
    bases = sum(end - start for start, end, *_ in blocks)
    return {"seqid": seqid, "strand": strand, "cds": blocks,
            "sequence": "ATG" + "AAA" * (bases // 3 - 2) + "TAA",
            "coverage": .99, "identity": .8,
            "problems": ["unanchored_genome_search"] if problems is None else problems,
            "support": [], "query": donor + "_" + query,
            "evidence": {"target": "T", "donor": donor, "query": query, "genome_only": True}}


def owner(gene, start, end, *, strand="+", cds=None, seqid="s"):
    return {"gene_id": gene, "seqid": seqid, "strand": strand,
            "start": start, "end": end, "cds": cds or [[start, end]]}


def test_shared_025679_009647_locus_keeps_gene_discovery_separate_from_copy_assignment():
    rows = [model("Trip", "025679", length=630), model("Trip", "009647", length=630),
            model("Dion", "q", length=630),
            model("Trip", "025679", 2000, 633), model("Trip", "009647", 2000, 633),
            model("Dion", "q", 2000, 633)]
    before = copy.deepcopy(rows)
    existing = [owner("original_raw_gene", 2000, 2633)]
    diagnostic = reassess_unanchored_models(rows, ownership=OwnershipIndex(existing))
    assert all(not row["problems"] for row in rows[:3])
    assert all(row["problems"] for row in rows[3:])
    assert diagnostic["annotation_aware_uniqueness"]["newly_supported_missing_annotation_records"] == 3
    for row, original in zip(rows, before, strict=True):
        assert (row["coverage"], row["identity"], row["support"], row["evidence"]) == (
            original["coverage"], original["identity"], original["support"], original["evidence"])
        evidence = row["placement_evidence"]
        assert evidence["orthology"] == evidence["expected_copy"] == "unassigned"
        assert evidence["original_all_loci_assessment"]["reason"] == "insufficient_unique_donor_support"
    admitted = consolidate(rows, existing, "T", ("Trip", "Dion"))
    assert len([row for row in admitted if row["status"] == "accepted"]) == 1
    assert all(row.get("revision_owner_ids", []) in ([], ["original_raw_gene"]) for row in admitted)


def test_019127_compatible_paths_keep_separate_support_after_annotated_copy_exclusion():
    # Shifted exon geometry reproduces the measured 263 shared in-frame bases;
    # sequences are synthetic, not copied from a research genome.
    short = [[1810, 1891, 0], [1641, 1678, 0], [704, 886, 2]]
    long = [[1800, 1891, 0], [704, 922, 2]]
    rows = [model("Trip", "019127", strand="-", cds=short),
            model("Dion", "q", strand="-", cds=long),
            model("Trip", "019127", 3000, 303, strand="-"),
            model("Dion", "q", 3000, 303, strand="-")]
    existing = [owner("existing_copy", 3000, 3303, strand="-")]
    reassess_unanchored_models(rows, ownership=OwnershipIndex(existing))
    assert not rows[0]["problems"] and not rows[1]["problems"]
    resolved = consolidate(rows, existing, "T", ("Trip", "Dion"))
    primary = next(row for row in resolved if row["status"] == "accepted")
    assert len(primary["sequence"]) == 309
    assert primary["locus_support"]["independent_donor_species"] == ["Dion", "Trip"]
    assert [row["donor"] for row in primary["support"]] == ["Dion"]
    alternative = primary["alternative_coding_paths"][0]
    assert [row["donor"] for row in alternative["support"]] == ["Trip"]
    assert primary["path_selection"]["representative_status"] == "ambiguous"


def test_010425_true_unannotated_alternatives_still_hold_even_after_owned_paths_excluded():
    rows = [model("Trip", "010425", length=1518), model("Nep", "short", length=555),
            model("Nep", "short", 3000, 555), model("Nep", "short", 5000, 555),
            model("Nep", "short", 8000, 555)]
    reassess_unanchored_models(rows, ownership=OwnershipIndex([owner("existing", 8000, 8555)]))
    assert all(row["problems"] for row in rows)
    assert rows[0]["placement_evidence"]["independent_donor_species"] == ["Trip"]
    assert rows[0]["placement_evidence"]["reason"] == "insufficient_unique_donor_support"


def test_owned_existing_model_revision_scoring_and_raw_gene_id_are_preserved():
    rows = [model("D", "q"), model("E", "q")]
    legacy = copy.deepcopy(rows)
    reassess_unanchored_models(legacy)
    existing = [owner("raw_gene%2C1", 0, 99)]
    reassess_unanchored_models(rows, ownership=OwnershipIndex(existing))
    assert [row["problems"] for row in rows] == [row["problems"] for row in legacy]
    assert rows[0]["placement_evidence"]["uniqueness_scope"] == "all_coding_loci_for_existing_model_revision"
    assert rows[0]["placement_evidence"]["original_annotation_owner_ids"] == ["raw_gene%2C1"]
    resolved = consolidate(rows, existing, "T")
    assert not any(row["status"] == "accepted" for row in resolved)
    assert all(row["revision_owner_ids"] == ["raw_gene%2C1"] for row in resolved)
    assert all(row["problems"] == ["overlap_existing_annotation"] for row in resolved)


def test_unknown_strand_annotation_remains_conservatively_owned():
    rows = [model("D", "q"), model("E", "q"), model("D", "q", 200), model("E", "q", 200)]
    reassess_unanchored_models(rows, ownership=OwnershipIndex([owner("unknown", 200, 299, strand=".")]))
    assert not rows[0]["problems"] and not rows[1]["problems"]
    assert rows[2]["placement_evidence"]["original_annotation_owner_ids"] == ["unknown"]


def test_opposite_strand_annotation_does_not_hide_a_real_unannotated_competitor():
    rows = [model("D", "q"), model("E", "q"), model("D", "q", 200), model("E", "q", 200)]
    reassess_unanchored_models(rows, ownership=OwnershipIndex([owner("opposite", 200, 299, strand="-")]))
    assert all(row["problems"] for row in rows)
    assert all(row["placement_evidence"]["uniqueness_scope"] == "unannotated_coding_loci" for row in rows)


def test_intron_nested_coding_locus_is_not_misclassified_as_existing_ownership():
    rows = [model("D", "q", 300), model("E", "q", 300), model("D", "q", 1200), model("E", "q", 1200)]
    existing = [owner("spanning", 0, 1000, cds=[[0, 99], [900, 1000]])]
    reassess_unanchored_models(rows, ownership=OwnershipIndex(existing))
    assert all(row["problems"] for row in rows)
    assert all(row["placement_evidence"]["uniqueness_scope"] == "unannotated_coding_loci" for row in rows)


def test_hard_qc_failure_neither_supplies_support_nor_blocks_unique_placement():
    rows = [model("D", "q"), model("E", "q"),
            model("D", "q", 200, problems=["frameshift", "unanchored_genome_search"])]
    reassess_unanchored_models(rows, ownership=OwnershipIndex([]))
    assert not rows[0]["problems"] and not rows[1]["problems"]
    assert rows[2]["problems"] == ["frameshift", "unanchored_genome_search"]
    assert "placement_evidence" not in rows[2]


def test_one_species_paralogs_and_self_never_meet_two_species_threshold():
    rows = [model("D", "q1"), model("D", "q2"), model("T", "self")]
    reassess_unanchored_models(rows, ownership=OwnershipIndex([]))
    assert all(row["problems"] for row in rows)
    assert rows[0]["placement_evidence"]["independent_donor_species"] == ["D"]


def test_incompatible_coding_bridge_remains_held():
    rows = [model("D", "q", 0), model("E", "q", 90), model("F", "q", 180)]
    reassess_unanchored_models(rows, ownership=OwnershipIndex([]))
    assert all(row["problems"] for row in rows)
    assert rows[0]["placement_evidence"]["reason"] == "incompatible_unanchored_paths"


def test_ownership_queries_are_bounded_by_distinct_coding_paths_not_donor_aliases():
    class CountingOwnership:
        calls = 0
        def overlapping(self, _path):
            self.calls += 1
            return []
    ownership = CountingOwnership()
    rows = [model(donor, f"q{i}") for donor in ("D", "E") for i in range(100)]
    reassess_unanchored_models(rows, ownership=ownership)
    assert ownership.calls == 1
    assert all(not row["problems"] for row in rows)
