import argparse
import csv

import pytest

from workflow.support import pairwise_synteny as synteny


def write_inputs(tmp_path, fasta=">Species_name_t1\nMPEPTIDE\n>Species_name_t2\nMPEPTIDEAAA\n"):
    source = tmp_path / "Species_name.protein.fa"
    source.write_text(fasta)
    gff = tmp_path / "Species_name.gff3"
    gff.write_text("##gff-version 3\nchr1\ttest\tgene\t1\t100\t.\t+\t.\tID=g1\n"
                   "chr1\ttest\tmRNA\t1\t40\t.\t+\t.\tID=t1;Parent=g1\n"
                   "chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=t2;Parent=g1\n")
    return {"species": "Species_name", "fasta": str(source), "gff": str(gff), "mode": "protein",
            "feature": "", "attribute": "", "genetic_code": None}


def test_prepare_selects_one_isoform_and_preserves_id_mapping(tmp_path):
    source = write_inputs(tmp_path)
    genes, metadata = synteny.prepare_genome(source, tmp_path, "target", 1)
    assert [(g.gene_id, g.start, g.end) for g in genes] == [("Species_name_t2", 0, 100)]
    assert metadata["collapsed_isoform_count"] == 1
    assert (tmp_path / "target.pep").read_text() == ">Species_name_t2\nMPEPTIDEAAA\n"
    with (tmp_path / "target.id_map.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["status"] for row in rows] == ["isoform_excluded", "selected"]
    assert all(row["locus_id"] == "g1" for row in rows)


@pytest.mark.parametrize("fasta,message", [
    (">t1\nMPEPTIDE\n>t1\nMPEPTIDE\n", "Duplicate"),
    (">Species_name_t1\nMPEPTIDE\n>t1\nMPEPTIDE\n", "Ambiguous"),
    (">unmapped\nMPEPTIDE\n", "GFF"),
    (">t1\nMP*EPTIDE\n", "internal stop"),
])
def test_invalid_fasta_annotation_mapping_fails(tmp_path, fasta, message):
    source = write_inputs(tmp_path, fasta)
    with pytest.raises(ValueError, match=message):
        synteny.prepare_genome(source, tmp_path, "target", 1)


def test_cds_translation_uses_selected_genetic_code(tmp_path):
    source = write_inputs(tmp_path, ">t1\nATGTAGTGA\n")
    source.update(mode="cds", genetic_code=6)
    genes, _ = synteny.prepare_genome(source, tmp_path, "query", 1)
    assert genes[0].gene_id == "Species_name_t1"
    assert (tmp_path / "query.pep").read_text() == ">Species_name_t1\nMQ\n"


def test_unknown_display_chromosome_is_rejected(tmp_path):
    source = write_inputs(tmp_path)
    genes, _ = synteny.prepare_genome(source, tmp_path, "target", 1)
    with pytest.raises(ValueError, match="annotated chromosomes"):
        synteny.selected_seqids("Chr1", genes)


def test_source_discovery_rejects_multiple_annotation_releases(tmp_path):
    directory = tmp_path / "input/species_gff"
    directory.mkdir(parents=True)
    for release in (1, 2):
        (directory / f"Species_name.release{release}.gff3").write_text("##gff-version 3\n")
    with pytest.raises(ValueError, match="exactly one"):
        synteny.source_file(tmp_path, "", "species_gff", "Species_name", (".gff3",))


def test_plan_rejects_duplicate_pair_ids_before_running_tools(tmp_path):
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text("analysis_id\ttarget_species\tquery_species\na\tTarget_species\tQuery_species\na\tTarget_species\tQuery_species\n")
    args = argparse.Namespace(workspace=tmp_path, pairs=pairs, sequence_mode="auto", genetic_code=1,
                              cscore=0.7, min_anchors=4, distance=20, minimum_mapping_fraction=1, formats="pdf")
    with pytest.raises(ValueError, match="unique"):
        synteny.build_plan(args)


def test_plot_settings_do_not_change_analysis_contract(tmp_path):
    plan = {"workspace": str(tmp_path), "parameters": {}, "tools": {"jcvi": "example"}, "formats": ["pdf"],
            "pairs": [{"analysis_id": "pair", "target_species": "Target_species", "query_species": "Query_species",
                       "target": {"fasta": "/input/t.fa", "gff": "/input/t.gff"},
                       "query": {"fasta": "/input/q.fa", "gff": "/input/q.gff"},
                       "target_seqids": "", "query_seqids": ""}]}
    before_analysis = synteny.contract_args(plan, "analysis")
    before_plots = synteny.contract_args(plan, "plots")
    plan["formats"] = ["png"]
    plan["pairs"][0]["query_seqids"] = "chr2,chr1"
    assert synteny.contract_args(plan, "analysis") == before_analysis
    assert synteny.contract_args(plan, "plots") != before_plots


def test_empty_anchors_are_not_a_completed_analysis(tmp_path):
    (tmp_path / "target.query.lifted.anchors").write_text("###\n")
    with pytest.raises(ValueError, match="No syntenic blocks"):
        synteny.summarize_anchors(tmp_path, ((), ()))


def test_changed_source_is_rejected_before_publication(tmp_path):
    source = tmp_path / "source.fa"
    source.write_text(">gene\nMPEPTIDE\n")
    plan = {"input_hashes": {str(source): synteny.digest(source)}}
    synteny.verify_inputs(plan)
    source.write_text(">gene\nMPEPTIDEX\n")
    with pytest.raises(ValueError, match="Input changed"):
        synteny.verify_inputs(plan)
