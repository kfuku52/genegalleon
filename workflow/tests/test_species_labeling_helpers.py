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


def test_synteny_neighbors_find_species_file_uses_exact_species_label(tmp_path):
    mod = load_module("synteny_neighbors.py", "synteny_neighbors_module")

    (tmp_path / "Species_a_subsp_x.fa").write_text(">x\nATG\n", encoding="utf-8")
    assert mod.find_species_file(str(tmp_path), "Species_a", mod.FASTA_EXTENSIONS) == ""

    target = tmp_path / "Species_a.fa"
    target.write_text(">a\nATG\n", encoding="utf-8")
    assert mod.find_species_file(str(tmp_path), "Species_a", mod.FASTA_EXTENSIONS) == str(target)


def test_synteny_neighbors_load_gene_info_replaces_numeric_chromosome_dtype(tmp_path):
    mod = load_module("synteny_neighbors.py", "synteny_neighbors_dtype_module")
    cache = tmp_path / "species.gff_info.tsv"
    cache.write_text(
        "gene_id\tchromosome\tstart\tend\tstrand\n"
        "gene_2\t2\t20\t29\t-\n"
        "gene_1\t1\t10\t19\t+\n",
        encoding="utf-8",
    )

    observed = mod.load_gene_info(str(cache))

    assert observed["chromosome"].tolist() == ["1", "2"]
    assert observed["gene_id"].tolist() == ["gene_1", "gene_2"]


def test_synteny_neighbors_load_gene_info_preserves_both_reversed_interval_bounds(tmp_path):
    mod = load_module("synteny_neighbors.py", "synteny_neighbors_reversed_interval_module")
    cache = tmp_path / "species.gff_info.tsv"
    cache.write_text(
        "gene_id\tchromosome\tstart\tend\tstrand\n"
        "gene_reversed\tchr1\t20\t10\t-\n",
        encoding="utf-8",
    )

    observed = mod.load_gene_info(str(cache))

    assert observed.loc[0, "start"] == 10
    assert observed.loc[0, "end"] == 20


def test_synteny_species_gene_cache_tracks_input_and_output_content(tmp_path):
    mod = load_module("synteny_neighbors.py", "synteny_neighbors_cache_module")
    cds = tmp_path / "Species_a.fa"
    gff = tmp_path / "Species_a.gff"
    output = tmp_path / "Species_a.gff_info.tsv"
    manifest = tmp_path / "Species_a.gff_info.tsv.provenance.json"
    cds.write_text(">a\nATG\n", encoding="utf-8")
    gff.write_text("chr1\t.\tCDS\t1\t3\t.\t+\t.\tID=a\n", encoding="utf-8")
    output.write_text("gene_id\nA\n", encoding="utf-8")
    contract = mod.species_gene_cache_contract("Species_a", str(cds), str(gff))

    assert not mod.species_gene_cache_is_current(str(output), str(manifest), contract)
    recorded = dict(contract, output_sha256=mod.sha256_file(str(output)))
    mod.write_species_gene_cache_manifest(str(manifest), recorded)
    assert manifest.is_file()
    assert mod.species_gene_cache_is_current(str(output), str(manifest), contract)

    output.write_text("gene_id\nB\n", encoding="utf-8")
    assert not mod.species_gene_cache_is_current(str(output), str(manifest), contract)
    output.write_text("gene_id\nA\n", encoding="utf-8")
    gff.write_text("chr1\t.\tCDS\t2\t4\t.\t+\t.\tID=a\n", encoding="utf-8")
    changed_contract = mod.species_gene_cache_contract("Species_a", str(cds), str(gff))
    assert not mod.species_gene_cache_is_current(str(output), str(manifest), changed_contract)


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

    assert mod.species_label_from_taxonomic_text("Dictyostelium cf. discoideum") == "Dictyostelium_cf._discoideum"
    assert mod.species_label_from_taxonomic_text("Bacillus subtilis subsp. subtilis") == "Bacillus_subtilis_subsp._subtilis"
    assert mod.species_label_from_taxonomic_text("Amoeba sp. JDS-Ruffled") == "Amoeba_sp._JDS-Ruffled"
    assert mod.species_label_from_taxonomic_text("Amoeba sp.") == "Amoeba_sp."
    assert mod.species_label_from_taxonomic_text("Asimitellaria furusei var. furusei") == "Asimitellaria_furusei_var._furusei"
    assert mod.species_label_from_taxonomic_text("Solanum lycopersicum cultivar Heinz 1706") == "Solanum_lycopersicum_cultivar_Heinz1706"
    assert mod.species_label_from_taxonomic_text("Escherichia coli serovar O157") == "Escherichia_coli_serovar_O157"
    assert (
        mod.species_label_from_taxonomic_text("Cenchrus americanus x Cenchrus purpureus")
        == "Cenchrus_americanus_x_Cenchrus_purpureus"
    )
    assert (
        mod.species_label_from_taxonomic_text("Cenchrus americanus \u00d7 Cenchrus purpureus")
        == "Cenchrus_americanus_x_Cenchrus_purpureus"
    )
    assert mod.extract_species_label("Citrus_x_limon_gene1") == "Citrus_x_limon"
    assert mod.extract_species_label("Citrus_\u00d7_limon_gene1") == "Citrus_x_limon"


def test_species_labeling_extracts_dotted_rank_labels_from_filenames():
    mod = load_module("species_labeling.py", "species_labeling_module")

    assert mod.extract_species_label("Arisaema_sp._aooni_longestCDS.fa.gz", strip_extension=True) == "Arisaema_sp._aooni"
    assert mod.extract_species_label("Amphizonella_sp_longestCDS.fa.gz", strip_extension=True) == "Amphizonella_sp"
    assert mod.extract_species_label("Cunea_sp_longestCDS_contamination_removal.fa.gz", strip_extension=True) == "Cunea_sp"
    assert mod.extract_species_label("Vannella_sp_longestCDS.transcript.fa.gz", strip_extension=True) == "Vannella_sp"
    assert mod.extract_species_label("Vexillifera_sp_longestCDS.busco.full.tsv", strip_extension=True) == "Vexillifera_sp"
    assert mod.extract_species_label("Homo_sapiens.isoform.fa.gz", strip_extension=True) == "Homo_sapiens"
    assert (
        mod.extract_species_label("Dictyostelium_cf_discoideum.longestCDS.fa.gz", strip_extension=True)
        == "Dictyostelium_cf_discoideum"
    )
    assert (
        mod.extract_species_label("Asimitellaria_furusei_var._subramosa_busco.full.tsv", strip_extension=True)
        == "Asimitellaria_furusei_var._subramosa"
    )
    assert (
        mod.extract_species_label("Asimitellaria_furusei_var._subramosa.fa.busco.full.tsv", strip_extension=True)
        == "Asimitellaria_furusei_var._subramosa"
    )
    assert (
        mod.extract_species_label("Asimitellaria_furusei_var._furusei.cds.fa.gz", strip_extension=True)
        == "Asimitellaria_furusei_var._furusei"
    )
    assert (
        mod.extract_species_label("homo_sapiens.GRCh38.cds.all.fa.gz", strip_extension=True)
        == "homo_sapiens"
    )
    assert (
        mod.extract_species_label("Arisaema_sp._aooni.dna.primary_assembly.fa.gz", strip_extension=True)
        == "Arisaema_sp._aooni"
    )
    assert (
        mod.extract_species_label("Cenchrus_americanus_x_Cenchrus_purpureus_metadata.tsv", strip_extension=True)
        == "Cenchrus_americanus_x_Cenchrus_purpureus"
    )
    assert mod.scientific_name_from_label("Asimitellaria_furusei_var._furusei") == "Asimitellaria furusei var. furusei"
    assert (
        mod.scientific_name_from_label("Cenchrus_americanus_x_Cenchrus_purpureus")
        == "Cenchrus americanus x Cenchrus purpureus"
    )
    assert (
        mod.base_species_label("Cenchrus_americanus_x_Cenchrus_purpureus")
        == "Cenchrus_americanus_x_Cenchrus_purpureus"
    )


def test_species_labeling_matches_species_label_exactly_after_suffix_stripping():
    mod = load_module("species_labeling.py", "species_labeling_module")

    assert mod.matches_species_label("Species_a.fa", "Species_a", strip_extension=True)
    assert not mod.matches_species_label("Species_a_subsp_x.fa", "Species_a", strip_extension=True)
