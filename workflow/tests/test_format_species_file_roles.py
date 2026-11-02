"""Roles and species identities must not depend on incidental substrings."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from build_download_manifest import discover_species_dir_based  # noqa: E402
from format_species_annotation.tasks import discover_generic_species_dir_tasks  # noqa: E402
from format_species_common import (  # noqa: E402
    is_probable_cds_filename,
    is_probable_genome_filename,
    provider_candidate_sort_key,
)
from format_species_download.local import resolve_local_manifest_row  # noqa: E402
from format_species_providers.common import (  # noqa: E402
    is_probable_cds_url,
    is_probable_genome_url,
    select_best_url_for_label,
)
from format_species_providers.web import (  # noqa: E402
    select_best_plantgarden_candidate,
    select_figshare_file_record,
)


@pytest.mark.parametrize("species", ["Cucurbita_pepo", "Euphorbia_peplus", "Peperomia_pellucida"])
@pytest.mark.parametrize("role,cds,genome", [
    ("cds", True, False), ("transcripts", True, False), ("genome", False, True),
    ("proteins", False, False), ("peptides", False, False), ("ncrna", False, False), ("rna", False, False),
    ("rna_from_genomic", False, False), ("translated_cds", False, False),
])
def test_species_substrings_do_not_override_actual_roles(species, role, cds, genome):
    name = species + "." + role + ".fa.gz"
    assert is_probable_cds_filename("local", name) is cds
    assert is_probable_genome_filename("local", name) is genome


@pytest.mark.parametrize("provider,name", [
    ("local", "Crocus_chrysanthus.fa"),
    ("local", "Chasmanthe_floribunda.fa"),
    ("figshare", "Haptanthus_hazlettii.fa"),
    ("local", "Plant_one.genome.proteins.fa"),
    ("local", "Crocus_chrysanthus.transcripts.fa"),
])
def test_genome_markers_are_complete_tokens_and_exclude_non_genomes(provider, name):
    assert not is_probable_genome_filename(provider, name)


@pytest.mark.parametrize("provider,name,expected", [
    ("ncbi", "GCF_000001735.4_TAIR10.1_genomic.fna.gz", "genome"),
    ("ncbi", "GCF_000001735.4_TAIR10.1_cds_from_genomic.fna.gz", "cds"),
    ("ensemblplants", "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz", "genome"),
    ("local", "Plant_one.ASM9206v1.fa", "genome"),
    ("local", "Plant_one.Chr1.fa", "genome"),
    ("local", "Plant_one.Chrall.fa", "genome"),
    ("local", "Plant_one.hap1.fa", "genome"),
    ("figshare", "Plant_one.haplotype2.fa", "genome"),
    ("fernbase", "Azfi_v1.1.fa", "genome"),
    ("fernbase", "Azfi_v1.1.cds.gz", "cds"),
    ("phycocosm", "MicrYARC1_GeneCatalog_CDS.fasta", "cds"),
    ("plantgarden", "APO_r1.1.pmol.fasta.gz", "genome"),
    ("veupathdb", "Plant_one_AnnotatedCDSs.fasta", "cds"),
])
def test_provider_native_names_remain_supported(provider, name, expected):
    assert is_probable_cds_filename(provider, name) is (expected == "cds")
    assert is_probable_genome_filename(provider, name) is (expected == "genome")


def test_known_species_prefix_is_excluded_before_role_checks_and_ranking():
    name = "Plant_transcript.genome.fa"
    assert not is_probable_cds_filename("local", name, species_key="Plant_transcript")
    assert is_probable_genome_filename("local", name, species_key="Plant_transcript")
    key = provider_candidate_sort_key("citrusgenomedb", "GFF", "Plant_go.gene.gff3", species_key="Plant_go")
    assert key[0] == 0
    assert not is_probable_genome_filename("local", "Plant_transcriptome.fa", species_key="Plant_transcript")


def test_candidate_ranking_does_not_find_pep_or_go_inside_species():
    assert provider_candidate_sort_key("citrusgenomedb", "GFF", "Gossypium_hirsutum.gene.gff3")[0] == 0
    for provider in ("local", "figshare", "citrusgenomedb"):
        assert provider_candidate_sort_key(provider, "CDS", "Cucurbita_pepo.cds.fa")[:-1] == (0, 0, 0)
    assert provider_candidate_sort_key("citrusgenomedb", "GFF", "Plant_one.go.gff3")[0] == 1


@pytest.mark.parametrize("directory", ["genome", "transcriptomes", "cds", "protein", "annotation"])
def test_url_directories_do_not_determine_roles(directory):
    url = "https://example.invalid/" + directory + "/Plant_one.genome.fa?download=1"
    assert not is_probable_cds_url("local", url)
    assert is_probable_genome_url("local", url)
    assert not is_probable_genome_url("local", "https://example.invalid/" + directory + "/Plant_one.fa")
    assert is_probable_cds_url("local", "https://example.invalid/" + directory + "/Plant_one%2Ecds.fa")


def test_citrus_documented_assembly_directory_only_accepts_non_annotation_fasta():
    base = "https://example.invalid/assembly/"
    assert is_probable_genome_url("citrusgenomedb", base + "opaque.fa")
    for name in ("Plant_one.cds.fa", "Plant_one.proteins.fa", "Plant_one.transcripts.fa"):
        assert not is_probable_genome_url("citrusgenomedb", base + name)


def test_encoded_filename_roles_are_also_used_when_ranking_urls():
    urls = ["https://example.invalid/a.transcripts.fa", "https://example.invalid/z%2Ecds.fa"]
    assert select_best_url_for_label("local", "CDS", urls) == urls[1]


@pytest.mark.parametrize("species", ["Cucurbita_pepo", "Euphorbia_peplus", "Peperomia_pellucida"])
def test_plantgarden_uses_common_roles_without_species_substring_exclusion(species):
    base = "https://example.invalid/"
    links = [base + species + suffix for suffix in
             (".genome.fa", ".cds.fa", ".transcripts.fa", ".proteins.fa", ".repeat.gff3",
              ".EDTA.gff3", ".a.longest.gff3", ".gene.gff3")]
    assert select_best_plantgarden_candidate(links, "GENOME") == base + species + ".genome.fa"
    assert select_best_plantgarden_candidate(links, "CDS") == base + species + ".cds.fa"
    assert select_best_plantgarden_candidate(links, "GFF") == base + species + ".gene.gff3"


@pytest.mark.parametrize("role", ["repeat", "repeatmask", "repeatmasker", "EDTA", "TE"])
def test_repeat_annotations_never_become_plantgarden_gene_gff(role):
    url = "https://example.invalid/Plant_one." + role + ".gff3"
    assert select_best_plantgarden_candidate([url], "GFF") == ""


@pytest.mark.parametrize("names", [
    ["Drosera_capensis.genome.fa"],
    ["Drosera_capensis.genome.fa", "Dionaea_muscipula.genome.fa"],
    ["Drosera_regiae.genome.fa"],
    ["Other_regia.genome.fa"],
    ["opaque.genome.fa"],
    ["Drosera_regia.v1.genome.fa", "Drosera_regia.v2.genome.fa"],
])
def test_figshare_rejects_missing_partial_and_ambiguous_species_matches(names):
    with pytest.raises(ValueError, match="set genome_filename"):
        select_figshare_file_record([{"name": n} for n in names], "GENOME", "", "Drosera_regia", "fixture")


@pytest.mark.parametrize("name", ["Drosera_regia.genome.fa", "Drosera-regia.genome.fa",
                                  "Drosera.regia.genome.fa", "DroseraRegia.genome.fa"])
def test_figshare_selects_unique_complete_species_match(name):
    target = {"name": name}
    files = [{"name": "Drosera_capensis.genome.fa"}, target]
    assert select_figshare_file_record(files, "GENOME", "", "Drosera_regia", "fixture") is target


def test_figshare_explicit_filename_retains_authoritative_selection():
    target = {"name": "opaque.fasta"}
    assert select_figshare_file_record([target], "GENOME", "opaque.fasta", "Drosera_regia", "fixture") is target


def test_figshare_duplicate_explicit_filenames_require_a_url():
    files = [{"name": "genome.fa", "id": 1}, {"name": "genome.fa", "id": 2}]
    with pytest.raises(ValueError, match="duplicate filename"):
        select_figshare_file_record(files, "GENOME", "genome.fa", "Drosera_regia", "fixture")


def test_discovery_manifest_and_local_resolution_agree_and_exclude_protein_fallback(tmp_path):
    species = "Crocus_chrysanthus"
    directory = tmp_path / species
    directory.mkdir()
    genome = directory / (species + ".genome.fa")
    genome.write_text(">chr1\nATGAAATAA\n")
    (directory / (species + ".gene.gff3")).write_text("chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=g1\n")
    for role in ("proteins", "ncrna", "peptides", "translated_cds"):
        (directory / (species + "." + role + ".fa")).write_text(">other\nMEVK\n")
    tasks, _warnings, errors = discover_generic_species_dir_tasks("local", tmp_path)
    assert not errors and tasks[0]["cds_path"] is None and tasks[0]["genome_path"] == genome
    rows, _warnings, errors = discover_species_dir_based("local", tmp_path)
    assert not errors and rows[0]["cds_filename"] == "" and rows[0]["genome_filename"] == genome.name
    row = resolve_local_manifest_row("local", str(directory), species, {}, tmp_path, [], 1)
    assert row["cds_filename"] == "" and row["genome_filename"] == genome.name


def test_transcript_fasta_is_never_reused_as_genome_in_generated_manifest(tmp_path):
    directory = tmp_path / "Crocus_chrysanthus"
    directory.mkdir()
    transcript = directory / "Crocus_chrysanthus.transcripts.fa"
    transcript.write_text(">t1\nATGAAATAA\n")
    rows, _warnings, errors = discover_species_dir_based("local", tmp_path)
    assert not errors and rows[0]["cds_filename"] == transcript.name and rows[0]["genome_filename"] == ""


def test_protein_only_directory_does_not_silently_become_cds(tmp_path):
    directory = tmp_path / "Plant_one"
    directory.mkdir()
    (directory / "Plant_one.proteins.fa").write_text(">p1\nMEVK\n")
    tasks, _warnings, errors = discover_generic_species_dir_tasks("local", tmp_path)
    assert not tasks and errors
    rows, _warnings, errors = discover_species_dir_based("local", tmp_path)
    assert not rows and errors
    with pytest.raises(ValueError, match="requires CDS"):
        resolve_local_manifest_row("local", str(directory), "Plant_one", {}, tmp_path, [], 1)
