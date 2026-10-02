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


@pytest.mark.parametrize("has_taxon", [True, False])
@pytest.mark.parametrize("dtype", [None, object])
def test_candidate_taxa_keep_first_leaf_and_requested_gene_order(has_taxon, dtype):
    branch = pd.Series(dict(orthogroup="OG1", branch_id=3, node_name="n",
                            gene_labels="missing; NA; g2; g1; g2", generax_event="H"))
    leaves = pd.DataFrame([
        dict(node_name="g1", taxon="Host species"),
        dict(node_name="g1", taxon="Other species"),
        dict(node_name="g2", taxon=pd.NA),
        dict(node_name="NA", taxon="Other species"),
    ], dtype=dtype)
    missing_taxon = str(leaves.loc[2, "taxon"])
    if not has_taxon:
        leaves = leaves.drop(columns="taxon")
    summary, records = scorer.summarize_candidate_branch(branch, leaves, pd.DataFrame(), [], resolver())
    assert summary["candidate_gene_count"] == 5
    assert summary["matched_leaf_count"] == 3
    assert [row["gene_id"] for row in records] == ["missing", "NA", "g2", "g1", "g2"]
    assert [row["gene_taxon"] for row in records] == (
        ["", "Other species", missing_taxon, "Host species", missing_taxon] if has_taxon else [""] * 5
    )
    assert records[3]["recipient_domain"] == ("Eukaryota" if has_taxon else "")


@pytest.mark.parametrize("branch_type", [int, float, str])
@pytest.mark.parametrize("categorical", [False, True])
def test_gene_aggregation_keeps_ties_branch_types_and_optional_defaults(branch_type, categorical):
    rows = pd.DataFrame([
        dict(orthogroup="OG2", gene_id="NA", candidate_branch_id=branch_type(10), gene_taxon="ten"),
        dict(orthogroup="OG2", gene_id="NA", candidate_branch_id=branch_type(2), gene_taxon="first two"),
        dict(orthogroup="OG2", gene_id="NA", candidate_branch_id=branch_type(2), gene_taxon="later two"),
        dict(orthogroup="OG1", gene_id="g", candidate_branch_id=branch_type(7), gene_taxon="one"),
        dict(orthogroup=None, gene_id="excluded", candidate_branch_id=branch_type(1), gene_taxon="missing group"),
    ], index=[6, 3, 5, 7, 1])
    if categorical:
        rows["orthogroup"] = pd.Categorical(rows.orthogroup, categories=["unused", "OG2", "OG1"])
    output = scorer.aggregate_gene_records(rows)
    assert output.columns.tolist() == scorer.GENE_OUTPUT_COLUMNS
    assert output.orthogroup.tolist() == (["OG2", "OG1"] if categorical else ["OG1", "OG2"])
    observed = output.set_index(["orthogroup", "gene_id"]).loc[("OG2", "NA")]
    assert observed.gene_taxon == ("ten" if branch_type is str else "first two")
    assert observed.candidate_branch_count == 2
    assert observed.candidate_branch_ids == ("10; 2" if branch_type is str else "2.0; 10.0" if branch_type is float else "2; 10")
    assert observed.besthit_accession == ""
    assert pd.isna(observed.intron_supported)
    assert not observed.expression_measured


def test_gene_aggregation_does_not_accept_missing_branch_identity():
    rows = pd.DataFrame([dict(orthogroup="OG1", gene_id="g", candidate_branch_id=pd.NA)])
    with pytest.raises(TypeError):
        scorer.aggregate_gene_records(rows)
