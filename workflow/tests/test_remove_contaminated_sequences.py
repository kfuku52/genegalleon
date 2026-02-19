import numpy

from workflow.script.remove_contaminated_sequences import resolve_rank_taxid_from_lineage


class FakeNcbi:
    def __init__(self, rank_map):
        self.rank_map = rank_map

    def get_rank(self, lineage_taxids):
        return {taxid: self.rank_map.get(taxid, "no rank") for taxid in lineage_taxids}


def test_resolve_rank_taxid_from_lineage_missing_returns_minus_one():
    ncbi = FakeNcbi({})
    assert resolve_rank_taxid_from_lineage(ncbi, "class", "nan") == -1
    assert resolve_rank_taxid_from_lineage(ncbi, "class", "") == -1


def test_resolve_rank_taxid_from_lineage_returns_rank_taxid_when_present():
    ncbi = FakeNcbi({2: "superkingdom", 2157: "class"})
    assert resolve_rank_taxid_from_lineage(ncbi, "class", "2;2157;9999") == 2157


def test_resolve_rank_taxid_from_lineage_returns_zero_when_rank_absent():
    ncbi = FakeNcbi({2: "superkingdom", 1224: "phylum"})
    assert resolve_rank_taxid_from_lineage(ncbi, "class", "2;1224;1236") == 0


def test_aligned_taxid_fallback_logic_matches_expected():
    lineage_to_rank_taxid = {"2;10;70": 70, "2;10;71": 0, "nan": -1}
    lineage = numpy.array(["2;10;70", "2;10;71", "nan"])
    lca = numpy.array([1000, 1001, 1002])
    rank_taxid = numpy.array([lineage_to_rank_taxid[x] for x in lineage], dtype=int)
    aligned = numpy.where(rank_taxid == -1, 0, numpy.where(rank_taxid == 0, lca, rank_taxid))
    assert aligned.tolist() == [70, 1001, 0]
