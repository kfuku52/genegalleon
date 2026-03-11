import numpy
import pandas

from workflow.support.remove_contaminated_sequences import (
    build_taxid_rank_lookup,
    resolve_lineage_rank_index,
    resolve_rank_taxid_from_lineage,
)


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


def test_resolve_rank_taxid_from_lineage_treats_domain_and_superkingdom_as_aliases():
    ncbi = FakeNcbi({2: "domain", 33090: "kingdom"})
    assert resolve_rank_taxid_from_lineage(ncbi, "superkingdom", "1;2;33090") == 2


def test_build_taxid_rank_lookup_and_bulk_rank_resolution():
    ncbi = FakeNcbi({2: "superkingdom", 1224: "phylum", 28211: "class"})
    taxid_to_rank = build_taxid_rank_lookup(ncbi, ["2;1224;28211", "2;1224"])
    assert resolve_rank_taxid_from_lineage(ncbi, "class", "2;1224;28211", taxid_to_rank=taxid_to_rank) == 28211


def test_resolve_lineage_rank_index_treats_domain_and_superkingdom_as_aliases():
    df_lineage = pandas.DataFrame(
        {
            "rank": ["no rank", "domain", "kingdom", "species"],
            "name": ["root", "Eukaryota", "Viridiplantae", "Arabidopsis thaliana"],
            "taxid": [1, 2759, 33090, 3702],
        }
    )
    assert resolve_lineage_rank_index(df_lineage, "superkingdom") == 1


def test_resolve_lineage_rank_index_raises_clear_error_when_rank_absent():
    df_lineage = pandas.DataFrame(
        {
            "rank": ["no rank", "domain", "kingdom", "species"],
            "name": ["root", "Eukaryota", "Viridiplantae", "Arabidopsis thaliana"],
            "taxid": [1, 2759, 33090, 3702],
        }
    )
    try:
        resolve_lineage_rank_index(df_lineage, "class")
    except ValueError as exc:
        assert "Requested rank 'class' was not found" in str(exc)
    else:
        raise AssertionError("Expected ValueError for missing rank")


def test_aligned_taxid_fallback_logic_matches_expected():
    lineage_to_rank_taxid = {"2;10;70": 70, "2;10;71": 0, "nan": -1}
    lineage = numpy.array(["2;10;70", "2;10;71", "nan"])
    lca = numpy.array([1000, 1001, 1002])
    rank_taxid = numpy.array([lineage_to_rank_taxid[x] for x in lineage], dtype=int)
    aligned = numpy.where(rank_taxid == -1, 0, numpy.where(rank_taxid == 0, lca, rank_taxid))
    assert aligned.tolist() == [70, 1001, 0]
