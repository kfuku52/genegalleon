from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"


def load_module(filename: str, module_name: str):
    script_path = SUPPORT_DIR / filename
    spec = spec_from_file_location(module_name, script_path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_get_trait_matrix_preserves_qualified_species_labels():
    mod = load_module("get_trait_matrix.py", "get_trait_matrix_module")

    assert mod.get_species_name("Dictyostelium_cf_discoideum_gene1.1") == "Dictyostelium_cf_discoideum"
    assert mod.get_gene_name("Dictyostelium_cf_discoideum_gene1.1") == "gene1.1"
    assert mod.trait_filename_to_species_name("Dictyostelium_cf_discoideum_trait.tsv") == "Dictyostelium_cf_discoideum"


def test_synteny_neighbors_guess_species_name_preserves_qualified_species_labels():
    mod = load_module("synteny_neighbors.py", "synteny_neighbors_module")

    assert mod.guess_species_name("Bacillus_subtilis_subsp_subtilis_geneA") == "Bacillus_subtilis_subsp_subtilis"
    assert mod.guess_species_name("Amoeba_sp_JDSRuffled_geneB") == "Amoeba_sp_JDSRuffled"


def test_get_promoter_fasta_extracts_qualified_species_label_from_gene_id():
    mod = load_module("get_promoter_fasta.py", "get_promoter_fasta_module")

    assert mod.extract_species_label("Dictyostelium_cf_discoideum_gene1") == "Dictyostelium_cf_discoideum"
    assert mod.extract_species_label("Amoeba_sp_JDSRuffled_gene2") == "Amoeba_sp_JDSRuffled"


def test_parse_grampa_matches_prefix_labeled_qualified_species_gene_names():
    mod = load_module("parse_grampa.py", "parse_grampa_module")
    species_names = ("Dictyostelium_cf_discoideum", "Arabidopsis_thaliana")
    _ordered, species_set, species_suffixes = mod.build_species_matcher(species_names)

    gt_id, species_gene_map = mod.summarize_gene_tree(
        (0, "(Dictyostelium_cf_discoideum_gene1,Arabidopsis_thaliana_gene2);"),
        species_names=species_names,
        species_set=species_set,
        species_suffixes=species_suffixes,
    )

    assert gt_id == "GT-1"
    assert species_gene_map["Dictyostelium_cf_discoideum"] == "gene1"
    assert species_gene_map["Arabidopsis_thaliana"] == "gene2"


def test_species_labeling_builds_qualified_labels_from_scientific_text():
    mod = load_module("species_labeling.py", "species_labeling_module")

    assert mod.species_label_from_taxonomic_text("Dictyostelium cf. discoideum") == "Dictyostelium_cf_discoideum"
    assert mod.species_label_from_taxonomic_text("Bacillus subtilis subsp. subtilis") == "Bacillus_subtilis_subsp_subtilis"
    assert mod.species_label_from_taxonomic_text("Amoeba sp. JDS-Ruffled") == "Amoeba_sp_JDSRuffled"
    assert mod.species_label_from_taxonomic_text("Amoeba sp.") == "Amoeba_sp_unknown"
    assert mod.species_label_from_taxonomic_text("Solanum lycopersicum cultivar Heinz 1706") == "Solanum_lycopersicum_cultivar_Heinz1706"
    assert mod.species_label_from_taxonomic_text("Escherichia coli serovar O157") == "Escherichia_coli_serovar_O157"
