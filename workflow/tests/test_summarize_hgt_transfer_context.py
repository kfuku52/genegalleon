import csv
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))

from gene_family_output_store import (  # noqa: E402
    GeneFamilyOutputStore,
    convert_storage_to_zip,
    orthogroup_id_from_name,
)
from summarize_hgt_transfer_context import (  # noqa: E402
    EVENT_COLUMNS,
    LINK_COLUMNS,
    parse_reconciliation,
    summarize,
    write_tsv,
)


def leaf(gene, species):
    return f'<clade><name>{gene}</name><eventsRec><leaf speciesLocation="{species}"/></eventsRec></clade>'


def fixture(tmp_path, support="90", num_leaf="5", transfer="Y@D_species@007"):
    later = ('<clade><name>later</name><eventsRec><branchingOut speciesLocation="R_species"/></eventsRec>'
             + leaf("r1", "R_species")
             + '<clade><name>t1</name><eventsRec><transferBack destinationSpecies="T_species"/>'
               '<leaf speciesLocation="T_species"/></eventsRec></clade></clade>')
    xml = ('<recPhylo><spTree><phylogeny><clade><name>root</name>'
           '<clade><name>D_species</name></clade><clade><name>007</name>'
           '<clade><name>R_species</name></clade><clade><name>S_species</name></clade></clade>'
           '<clade><name>T_species</name></clade></clade></phylogeny></spTree>'
           '<recGeneTree><phylogeny><clade><name>parent</name>'
           '<eventsRec><branchingOut speciesLocation="D_species"/></eventsRec>'
           '<clade><name>copies</name><eventsRec><duplication speciesLocation="D_species"/></eventsRec>'
           + leaf("d1", "D_species") + leaf("d2", "D_species") + '</clade>'
           '<clade><name>recipient</name><eventsRec><transferBack destinationSpecies="007"/>'
           '<speciation speciesLocation="007"/></eventsRec>' + later + leaf("s1", "S_species")
           + '</clade></clade></phylogeny></recGeneTree></recPhylo>')
    root = tmp_path / "families"
    (root / "generax_xml").mkdir(parents=True)
    (root / "generax_xml/OG0001_generax.xml").write_text(xml)
    (root / "stat_branch").mkdir()
    branch = dict(orthogroup="OG0001", branch_id="9", node_name="parent", generax_transfer=transfer,
                  candidate_genes="d1; d2; r1; s1; t1")
    stat = dict(branch_id="9", gene_labels=branch["candidate_genes"], support_generax_ufboot=support,
                num_leaf=num_leaf, so_event="S")
    write_tsv(root / "stat_branch/OG0001_stat.branch.tsv", [stat], list(stat))
    genes = []
    for gene, species, scaffold in [("d1", "D_species", "donor_scaf"), ("d2", "D_species", "donor_scaf"),
                                    ("r1", "R_species", "r_scaf"), ("s1", "S_species", "s_scaf"),
                                    ("t1", "T_species", "t_scaf")]:
        row = dict(orthogroup="OG0001", gene_id=gene, gene_taxon=species, host_scaffold_id=scaffold,
                   host_scaffold_status="measured", host_scaffold_count_unit="gff_locus",
                   host_scaffold_locus_id=gene, intron_supported="", synteny_support_score="")
        prefix = "host_scaffold_background_class_"
        row.update({prefix + k: v for k, v in dict(total_count=20, compatible_count=18,
                    incompatible_count=1, unresolved_count=1, cds_id_count=0,
                    classified_fraction=.95, compatible_fraction=18/19, compatible_all_fraction=.9).items()})
        genes.append(row)
    return root, branch, genes


def test_exact_event_roles_subsequent_transfer_and_scaffold_copy_deduplication(tmp_path):
    root, branch, genes = fixture(tmp_path)
    events, links = summarize([branch], genes, GeneFamilyOutputStore(root))
    event = events[0]
    assert event["event_id"] == "OG0001:9:1"
    assert event["mapping_status"] == "matched"
    assert event["support_generax_ufboot"] == 90
    assert event["donor_retained_gene_count"] == 2
    assert event["donor_scaffold_count"] == 1
    assert event["donor_host_scaffold_background_class_total_count"] == 20
    assert event["recipient_retained_genes"] == "r1; s1"
    assert event["recipient_host_scaffold_background_class_total_count"] == 40
    assert event["recipient_evidence_basis"] == "extant_descendant_proxy"
    assert event["donor_evidence_basis"] == "extant_terminal_genome"
    moved = next(row for row in links if row["gene_id"] == "t1")
    assert moved["lineage_status"] == "transferred_out"
    assert not moved["eligible_for_context"]
    assert next(row for row in links if row["gene_id"] == "d1")["intron_supported"] == ""


@pytest.mark.parametrize("support,num_leaf,status", [("", "5", "missing_ufboot"),
                         ("100", "1", "terminal_branch"), ("89.9", "5", "measured")])
def test_missing_terminal_and_inclusive_support_are_not_imputed(tmp_path, support, num_leaf, status):
    root, branch, genes = fixture(tmp_path, support=support, num_leaf=num_leaf)
    events, _ = summarize([branch], genes, GeneFamilyOutputStore(root))
    assert events[0]["support_status"] == status
    if status != "measured":
        assert "support_generax_ufboot" not in events[0]
    else:
        assert events[0]["support_generax_ufboot"] == 89.9


def test_each_transfer_token_evaluated_independently_without_besthit_proxy(tmp_path):
    root, branch, genes = fixture(tmp_path, transfer="Y@D_species@007; Y@T_species@007")
    for gene in genes:
        gene["donor_class"] = "Insecta"
    events, _ = summarize([branch], genes, GeneFamilyOutputStore(root))
    assert [e["event_index"] for e in events] == [1, 2]
    assert events[0]["mapping_status"] == "matched"
    assert events[1]["mapping_reason"] == "xml_event_missing"
    assert "donor_retained_gene_count" not in events[1]


def test_wrong_family_branch_join_does_not_reuse_support(tmp_path):
    root, branch, genes = fixture(tmp_path)
    branch["branch_id"] = "10"
    events, _ = summarize([branch], genes, GeneFamilyOutputStore(root))
    assert events[0]["support_status"] == "missing_stat_branch"
    assert "support_generax_ufboot" not in events[0]


def test_missing_gene_context_keeps_partial_measurement(tmp_path):
    root, branch, genes = fixture(tmp_path)
    events, links = summarize([branch], genes[1:], GeneFamilyOutputStore(root))
    assert events[0]["donor_context_status"] == "partial"
    assert events[0]["donor_retained_gene_count"] == 2
    assert events[0]["donor_mapped_gene_count"] == 1
    missing = next(r for r in links if r["gene_id"] == "d1")
    assert missing["host_scaffold_status"] == "gene_not_in_summary"
    assert "intron_supported" not in missing


def test_missing_xml_keeps_unresolved_event_and_blank_counts(tmp_path):
    root, branch, genes = fixture(tmp_path)
    (root / "generax_xml/OG0001_generax.xml").unlink()
    events, links = summarize([branch], genes, GeneFamilyOutputStore(root))
    assert events[0]["mapping_reason"] == "missing_generax_xml"
    assert "donor_scaffold_count" not in events[0]
    assert links == []


def test_species_tree_internal_labels_and_external_mismatch(tmp_path):
    root, branch, genes = fixture(tmp_path)
    tree = tmp_path / "species.nwk"
    tree.write_text("(D_species,(R_species,S_species)007,T_species)root;")
    events, _ = summarize([branch], genes, GeneFamilyOutputStore(root), str(tree))
    assert events[0]["species_tree_mapping_status"] == "matched_external_tree"
    tree.write_text("(D_species,(R_species,T_species)007,S_species)root;")
    events, links = summarize([branch], genes, GeneFamilyOutputStore(root), str(tree))
    assert events[0]["mapping_reason"] == "external_tree_mismatch"
    assert links == []


def test_corrupt_scaffold_metrics_and_duplicate_keys_fail(tmp_path):
    root, branch, genes = fixture(tmp_path)
    with pytest.raises(ValueError, match="Duplicate gene-summary"):
        summarize([branch], genes + genes[:1], GeneFamilyOutputStore(root))
    genes[0]["host_scaffold_background_class_compatible_count"] = 21
    with pytest.raises(ValueError, match="Invalid scaffold count"):
        summarize([branch], genes, GeneFamilyOutputStore(root))


def test_duplicate_species_labels_are_not_arbitrarily_resolved(tmp_path):
    root, _, _ = fixture(tmp_path)
    xml = (root / "generax_xml/OG0001_generax.xml").read_text()
    xml = xml.replace("<name>S_species</name>", "<name>R_species</name>")
    with pytest.raises(ValueError, match="Duplicate reconciliation species-tree"):
        parse_reconciliation(xml)


def test_cli_output_schema_including_empty_results(tmp_path):
    event_path, link_path = tmp_path / "events.tsv", tmp_path / "links.tsv"
    write_tsv(event_path, [], EVENT_COLUMNS)
    write_tsv(link_path, [], LINK_COLUMNS)
    assert next(csv.reader(event_path.open(), delimiter="\t")) == EVENT_COLUMNS
    assert next(csv.reader(link_path.open(), delimiter="\t")) == LINK_COLUMNS


def test_zip_backed_reconciliation_and_stat_branch_match_raw_outputs(tmp_path):
    root, branch, genes = fixture(tmp_path)
    expected = summarize([branch], genes, GeneFamilyOutputStore(root))
    convert_storage_to_zip(root, "orthogroup", ["OG0001"], orthogroup_id_from_name)
    assert not (root / "generax_xml/OG0001_generax.xml").exists()
    assert summarize([branch], genes, GeneFamilyOutputStore(root)) == expected


def test_ambiguous_xml_origin_is_held_instead_of_cross_joining(tmp_path):
    root, branch, genes = fixture(tmp_path)
    path = root / "generax_xml/OG0001_generax.xml"
    path.write_text(path.read_text().replace('<branchingOut speciesLocation="D_species"/>',
                    '<branchingOut speciesLocation="D_species"/><branchingOut speciesLocation="T_species"/>', 1))
    events, links = summarize([branch], genes, GeneFamilyOutputStore(root))
    assert events[0]["mapping_reason"] == "xml_event_ambiguous"
    assert links == []


def test_species_mismatch_holds_scaffold_but_preserves_source_gene_observations(tmp_path):
    root, branch, genes = fixture(tmp_path)
    genes[0].update(gene_taxon="Other_species", intron_supported="False")
    _, links = summarize([branch], genes, GeneFamilyOutputStore(root))
    row = next(r for r in links if r["gene_id"] == "d1")
    assert row["host_scaffold_status"] == "gene_species_mismatch"
    assert "host_scaffold_id" not in row
    assert row["intron_supported"] == "False"
