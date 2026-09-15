import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))

import annotate_hgt_tree_plot as annotator  # noqa: E402
import score_hgt_candidates as scorer  # noqa: E402


class Taxonomy:
    def get_name_translator(self, names):
        ids = {"Host species": [3], "Other species": [4], "Ambiguous species": [3, 4]}
        return {name: ids[name] for name in names if name in ids}

    def get_lineage(self, taxid):
        return {3: [1, 2, 10, 11, 3], 4: [1, 4]}.get(taxid, [])

    def get_rank(self, ids):
        ranks = {1: "no rank", 2: "domain", 3: "species", 4: "species", 10: "clade", 11: "clade"}
        return {taxid: ranks[taxid] for taxid in ids}

    def get_taxid_translator(self, ids):
        names = {1: "root", 2: "Eukaryota", 3: "Host species", 4: "Other species", 10: "Outer", 11: "Inner"}
        return {taxid: names[taxid] for taxid in ids}


def resolver():
    result = scorer.TaxonomyResolver("")
    result.enabled = True
    result.ncbi = Taxonomy()
    return result


def test_unresolved_domain_is_not_a_mismatch():
    r = resolver()
    assert pd.isna(r.compare("Host species", 4)["same_superkingdom"])
    assert r.compare("Host species", 3)["same_superkingdom"] == 1
    result = scorer.besthit_support_from_leaf_rows(pd.DataFrame([
        dict(node_name="g", taxon="Host species", organism="Other species", taxid_y=4),
    ]), r)
    assert result["gene_count"] == 1
    assert pd.isna(result["same_superkingdom_fraction"])


def test_ambiguous_taxonomy_name_is_not_arbitrarily_resolved():
    r = resolver()
    assert r.resolve_name_taxid("Ambiguous species") == 0
    assert r.resolve_name_taxid("Ambiguous species") == 0  # cached outcome
    assert scorer.resolve_taxonomy_annotation("Ambiguous species", "", r)["domain"] == ""


def test_complete_lineage_keeps_nested_clades_and_cache():
    r = resolver()
    for _ in range(2):
        annotation = scorer.resolve_taxonomy_annotation("Host species", "", r)
        assert annotation["taxonomy"] == "domain:Eukaryota; clade:Outer; clade:Inner; species:Host species"
        assert annotation["species"] == "Host species"


@pytest.mark.parametrize("missing", [None, float("nan"), pd.NA, ""])
def test_missing_besthits_are_not_counted(missing):
    result = scorer.besthit_support_from_leaf_rows(pd.DataFrame([
        dict(node_name="g", taxon="Host species", organism=missing, sprot_best=missing, taxid_y=missing),
    ]), resolver())
    assert result["gene_count"] == 0
    assert result["per_gene"] == {}


def test_unmeasured_intron_remains_missing_in_tree_evidence():
    branch = pd.Series(dict(orthogroup="OG1", branch_id=3, node_name="n", gene_labels="g", generax_event="H"))
    leaves = pd.DataFrame([dict(node_name="g", so_event="L", taxon="Host species", num_intron=pd.NA)])
    _, records = scorer.summarize_candidate_branch(branch, leaves, pd.DataFrame(), [], resolver())
    genes = scorer.aggregate_gene_records(pd.DataFrame(records))
    assert pd.isna(genes.iloc[0].intron_supported)
    annotated = annotator.annotate_leaf_rows(leaves, genes, "OG1")
    assert pd.isna(annotated.iloc[0].hgt_Intron)
