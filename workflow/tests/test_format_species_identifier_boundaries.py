"""Adversarial identifiers must not select another assembly or discard nuclear data."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from build_download_manifest import (  # noqa: E402
    discover_species_dir_based,
    infer_coge_genome_id_from_files,
    infer_ddbj_wgs_id_from_files,
    infer_gwh_accession_from_files,
)
from format_species_accessions import gwh_accession_matches, select_ncbi_bundle  # noqa: E402
from format_species_annotation.organelle import is_organelle_annotation  # noqa: E402
from format_species_discovery import discover_ncbi_like_tasks  # noqa: E402
from format_species_provider_urls import extract_ncbi_accession_from_source_id  # noqa: E402
from format_species_providers import gwh, web  # noqa: E402
from format_species_providers.catalogs import extract_plantaedb_ncbi_accession  # noqa: E402


@pytest.mark.parametrize("with_cds", [False, True])
def test_ncbi_versions_are_kept_separate_in_runtime_and_generated_manifest(tmp_path, with_cds):
    directory = tmp_path / "Plant_one"
    directory.mkdir()
    for name in ["GCA_000001.1_old_genomic.gff", "GCA_000001.10_new_genomic.fna"] + (
            ["GCA_000001.1_old_cds_from_genomic.fna"] if with_cds else []):
        (directory / name).write_text(">chr1\nATG\n" if name.endswith("fna") else
                                      "chr1\ts\tCDS\t1\t3\t.\t+\t0\tID=t1\n")
    tasks, _warnings, errors = discover_ncbi_like_tasks(tmp_path, "ncbi")
    rows, _warnings, manifest_errors = discover_species_dir_based("ncbi", tmp_path)
    if with_cds:
        assert not errors and not manifest_errors
        assert tasks[0]["source_bundle_id"] == "GCA_000001.1"
        assert tasks[0]["genome_path"] is None and rows[0]["genome_filename"] == ""
    else:
        assert not tasks and not rows
        assert "incompatible assembly" in errors[0] and "incompatible assembly" in manifest_errors[0]


def test_ncbi_complete_versions_never_share_members():
    a, b = "GCA_000001.1", "GCA_000001.10"
    groups = tuple([Path(a + suffix), Path(b + suffix)] for suffix in
                   ("_cds_from_genomic.fna", "_genomic.gff", "_genomic.gbff", "_genomic.fna"))
    selected, chosen, ids = select_ncbi_bundle(groups)
    assert selected == a and ids == (a, b)
    assert all(len(group) == 1 and group[0].name.startswith(a + "_") for group in chosen)


@pytest.mark.parametrize("text,expected", [
    ("GCA_000001.1", "GCA_000001.1"),
    ("https://example.invalid/GCA_000001.10_genomic.fna", "GCA_000001.10"),
    ("NotGCA_000001.1", ""), ("GCA_000001.1wrong", ""),
    ("GCA_000001.1.10", ""), ("GCA_000001.1_GCA_000002.1", ""),
])
def test_ncbi_extraction_requires_one_complete_accession(text, expected):
    assert extract_ncbi_accession_from_source_id(text) == expected


@pytest.mark.parametrize("species", ["Oryza_sativa", "Pinus_taeda", "Plant_one"])
def test_ddbj_does_not_turn_species_and_filename_words_into_wgs_ids(species):
    assert infer_ddbj_wgs_id_from_files(species, [Path(species + ".genome.fa")], []) == ""


@pytest.mark.parametrize("species,filename,expected", [
    ("ABCD", "genome.fa", "ABCD00000000"),
    ("abcd", "genome.fa", "ABCD00000000"),
    ("Plant_one", "Plant_one.ABCD01000001.1.genome.fa", "ABCD00000000"),
    ("Plant_one", "GHIJKL020000001.genome.fa", "GHIJKL000000000"),
])
def test_ddbj_explicit_prefix_and_full_accessions_still_work(species, filename, expected):
    assert infer_ddbj_wgs_id_from_files(species, [Path(filename)], []) == expected


def test_gwh_inference_and_matches_require_token_and_version_boundaries():
    assert infer_gwh_accession_from_files("Plant_one", [Path("notGWHFAKE.genome.fa")], []) == ""
    assert not gwh_accession_matches("Plant_GWHABCD00000000.10", "GWHABCD00000000.1")
    assert gwh_accession_matches("Plant_GWHABCD00000000.1", "GWHABCD00000000.1")
    assert gwh_accession_matches("Plant_GWHABCD00000000.10", "GWHABCD00000000")
    assert not gwh_accession_matches("GWHABCD00000000.1 GWHABCD00000000.10", "GWHABCD00000000.1")
    with pytest.raises(ValueError, match="multiple assembly"):
        gwh.extract_gwh_accession_candidate("GWHABCD00000000.1_GWHABCD00000000.10")


def test_gwh_folder_lookup_rejects_a_prefix_version_match(monkeypatch):
    root = "https://example.invalid/gwh/"
    monkeypatch.setattr(gwh, "resolve_gwh_download_base_url", lambda: root)
    monkeypatch.setattr(gwh, "fetch_text_with_headers", lambda url, timeout, headers:
                        '<a href="Plants/">Plants</a>' if url == root else
                        '<a href="Plant_GWHABCD00000000.10/">other version</a>')
    with pytest.raises(ValueError, match="not found"):
        gwh.resolve_gwh_folder_url_from_id("GWHABCD00000000.1", 1, {})


def test_gwh_show_page_lookup_rejects_stem_only_matches(monkeypatch):
    monkeypatch.setattr(gwh, "fetch_json_with_headers", lambda *args, **kwargs:
                        {"data": [{"url": "https://example.invalid/show", "attrs": {}}]})
    monkeypatch.setattr(gwh, "resolve_gwh_bundle_from_show_page", lambda *args, **kwargs:
                        {"genome_url": "https://example.invalid/GWHABCD00000000.10.genome.fa"})
    with pytest.raises(ValueError, match="did not match"):
        gwh.resolve_gwh_show_page_url_from_id("GWHABCD00000000.1", 1, {})


def veupath_record(identifier, species="Plant one", genome=True):
    return {"attributes": {"primary_key": identifier, "species": species,
            "URLGenomeFasta": f"https://example.invalid/Current_Release/{identifier}/{identifier}_Genome.fasta" if genome else "",
            "URLgff": f"https://example.invalid/Current_Release/{identifier}/{identifier}.gff"}}


def test_veupathdb_requires_exact_source_identity_not_partial_or_species_fallback(monkeypatch):
    records = [veupath_record("PlantV10")]
    monkeypatch.setattr(web, "fetch_veupathdb_records", lambda timeout, headers: records)
    with pytest.raises(ValueError, match="not found"):
        web.resolve_veupathdb_download_urls_from_id("PlantV1", "Plant_one", 1, {})
    records.append(veupath_record("PlantV1"))
    result = web.resolve_veupathdb_download_urls_from_id("PlantV1", "Plant_one", 1, {})
    assert "/PlantV1/" in result["genome_url"]
    with pytest.raises(ValueError, match="not found"):
        web.resolve_veupathdb_download_urls_from_id("Plant.V1", "", 1, {})
    with pytest.raises(ValueError, match="multiple assemblies"):
        web.resolve_veupathdb_download_urls_from_id("Plant_one", "", 1, {})


def test_veupathdb_does_not_replace_an_exact_incomplete_record_with_another_version(monkeypatch):
    monkeypatch.setattr(web, "fetch_veupathdb_records", lambda timeout, headers:
                        [veupath_record("PlantV10"), veupath_record("PlantV1", genome=False)])
    with pytest.raises(ValueError, match="GFF/genome"):
        web.resolve_veupathdb_download_urls_from_id("PlantV1", "", 1, {})


def test_veupathdb_record_cannot_mix_assembly_versions(monkeypatch):
    record = veupath_record("PlantV1")
    record["attributes"]["URLgff"] = record["attributes"]["URLgff"].replace("PlantV1", "PlantV10")
    monkeypatch.setattr(web, "fetch_veupathdb_records", lambda timeout, headers: [record])
    with pytest.raises(ValueError, match="mixes genome and GFF"):
        web.resolve_veupathdb_download_urls_from_id("PlantV1", "", 1, {})


@pytest.mark.parametrize("infer,names", [
    (infer_coge_genome_id_from_files, ["Plant_one.gid1.cds.fa", "Plant_one.gid10.gff"]),
    (infer_gwh_accession_from_files, ["GWHABCD00000000.1.cds.fa", "GWHABCD00000000.10.gff"]),
    (infer_ddbj_wgs_id_from_files, ["ABCD01000001.genome.fa", "EFGH01000001.cds.fa"]),
])
def test_manifest_source_id_inference_rejects_conflicting_sources(infer, names):
    with pytest.raises(ValueError, match="multiple"):
        infer("Plant_one", [Path(name) for name in names], [])


def test_ncbi_links_need_real_hostnames_not_matching_url_text():
    page = '<a href="https://ncbi.nlm.nih.gov.example.invalid/GCA_000001.1">spoof</a>'
    assert extract_plantaedb_ncbi_accession("https://example.invalid/page", page) == ""
    page += '<a href="https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000002.1">valid</a>'
    assert extract_plantaedb_ncbi_accession("https://example.invalid/page", page) == "GCA_000002.1"


def test_ncbi_visible_text_fallback_excludes_hidden_script_identifiers():
    page = '<script>const unrelated="GCA_000001.1";</script><p>Assembly: GCA_000002.1</p>'
    assert extract_plantaedb_ncbi_accession("https://example.invalid/page", page) == "GCA_000002.1"


def test_ncbi_custom_service_authority_is_exact_including_port(monkeypatch):
    monkeypatch.setenv("GG_NCBI_DATASETS_BASE_URL", "http://127.0.0.1:12345/datasets")
    page = '<a href="http://127.0.0.1:12346/datasets/genome/GCA_000001.1">wrong port</a>'
    assert extract_plantaedb_ncbi_accession("https://example.invalid/page", page) == ""
    page += '<a href="http://127.0.0.1:12345/datasets/genome/GCA_000002.1">configured</a>'
    assert extract_plantaedb_ncbi_accession("https://example.invalid/page", page) == "GCA_000002.1"


@pytest.mark.parametrize("value", ["non-mitochondrial", "nuclear chloroplast insertion", "chloroplast-like nuclear",
                                  "notmitochondrial", "nuclear"])
def test_nuclear_or_negated_genome_labels_do_not_trigger_organelle_exclusion(value):
    assert not is_organelle_annotation("region", {"genome": (value,)})


@pytest.mark.parametrize("value", ["mitochondrion", "mitochondrial DNA", "chloroplast", "chloroplast-genome",
                                  "plastid", "apicoplast", "mitochondria"])
def test_positive_organelle_evidence_remains_supported(value):
    assert is_organelle_annotation("region", {"genome": (value.lower(),)})
