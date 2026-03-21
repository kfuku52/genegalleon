import importlib.util
from pathlib import Path
import zipfile
import json
import threading
import subprocess
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

import pandas


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "generate_species_trait.py"


def load_script_module():
    spec = importlib.util.spec_from_file_location("generate_species_trait_module", SCRIPT_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_script(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_generate_species_trait_from_manifest_bulk_source(tmp_path):
    manifest = tmp_path / "input" / "input_generation" / "download_plan.tsv"
    write_text(
        manifest,
        (
            "provider\tid\tspecies_key\n"
            "ncbi\tGCF_000001405.40 (Homo sapiens)\tHomo_sapiens\n"
            "ncbi\tGCF_000001635.27 (Mus musculus)\t\n"
        ),
    )

    trait_source = tmp_path / "downloads" / "trait_datasets" / "austraits" / "austraits.tsv"
    write_text(
        trait_source,
        (
            "species\tbody_mass_g\n"
            "Homo_sapiens\t70\n"
            "Homo_sapiens\t72\n"
            "Mus_musculus\t22\n"
            "Pan_troglodytes\t45\n"
        ),
    )

    trait_plan = tmp_path / "input" / "input_generation" / "trait_plan.tsv"
    write_text(
        trait_plan,
        (
            "database\tsource_column\toutput_trait\tvalue_type\taggregation\n"
            "austraits\tbody_mass_g\tbody_mass_g\tnumeric\tmedian\n"
        ),
    )

    db_sources = tmp_path / "input" / "input_generation" / "trait_database_sources.tsv"
    write_text(
        db_sources,
        (
            "database\tacquisition_mode\turi\tspecies_column\tdelimiter\n"
            f"austraits\tbulk\t{trait_source}\tspecies\ttsv\n"
        ),
    )

    output = tmp_path / "input" / "species_trait" / "species_trait.tsv"
    completed = run_script(
        "--download-manifest",
        str(manifest),
        "--trait-plan",
        str(trait_plan),
        "--database-sources",
        str(db_sources),
        "--downloads-dir",
        str(tmp_path / "downloads" / "trait_datasets"),
        "--output",
        str(output),
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout
    assert output.exists()

    df = pandas.read_csv(output, sep="\t")
    assert list(df.columns) == ["species", "body_mass_g"]
    assert df.shape[0] == 2
    mass_by_species = dict(zip(df["species"], df["body_mass_g"]))
    assert mass_by_species["Homo_sapiens"] == 71
    assert mass_by_species["Mus_musculus"] == 22


def test_generate_species_trait_from_species_cds_source(tmp_path):
    species_cds = tmp_path / "output" / "input_generation" / "species_cds"
    species_cds.mkdir(parents=True, exist_ok=True)
    (species_cds / "Canis_lupus.fa").write_text(">a\nATG\n", encoding="utf-8")
    (species_cds / "Felis_catus.fna.gz").write_text(">b\nATG\n", encoding="utf-8")

    trait_source = tmp_path / "downloads" / "trait_datasets" / "animaltraits" / "animaltraits.tsv"
    write_text(
        trait_source,
        (
            "scientific_name\tlifespan_years\n"
            "Canis lupus\t13\n"
            "Felis catus\t15\n"
            "Panthera leo\t14\n"
        ),
    )

    trait_plan = tmp_path / "input" / "input_generation" / "trait_plan.tsv"
    write_text(
        trait_plan,
        (
            "database\tsource_column\toutput_trait\tvalue_type\taggregation\n"
            "animaltraits\tlifespan_years\tlifespan_years\tnumeric\tmedian\n"
        ),
    )

    db_sources = tmp_path / "input" / "input_generation" / "trait_database_sources.tsv"
    write_text(
        db_sources,
        (
            "database\tacquisition_mode\turi\tspecies_column\tdelimiter\n"
            f"animaltraits\tbulk\t{trait_source}\tscientific_name\ttsv\n"
        ),
    )

    output = tmp_path / "species_trait.tsv"
    completed = run_script(
        "--species-source",
        "species_cds",
        "--species-cds-dir",
        str(species_cds),
        "--trait-plan",
        str(trait_plan),
        "--database-sources",
        str(db_sources),
        "--downloads-dir",
        str(tmp_path / "downloads" / "trait_datasets"),
        "--output",
        str(output),
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout
    df = pandas.read_csv(output, sep="\t")
    assert set(df["species"].tolist()) == {"Canis_lupus", "Felis_catus"}
    value_map = dict(zip(df["species"], df["lifespan_years"]))
    assert value_map["Canis_lupus"] == 13
    assert value_map["Felis_catus"] == 15


def test_generate_species_trait_species_api_mode_downloads_target_species_only(tmp_path):
    manifest = tmp_path / "input" / "input_generation" / "download_plan.tsv"
    write_text(
        manifest,
        (
            "provider\tid\tspecies_key\n"
            "local\tid1\tHomo_sapiens\n"
            "local\tid2\tMus_musculus\n"
        ),
    )

    api_root = tmp_path / "api"
    write_text(
        api_root / "Homo_sapiens.tsv",
        "species\tis_carnivorous\nHomo_sapiens\t0\n",
    )
    write_text(
        api_root / "Mus_musculus.tsv",
        "species\tis_carnivorous\nMus_musculus\t0\n",
    )

    trait_plan = tmp_path / "input" / "input_generation" / "trait_plan.tsv"
    write_text(
        trait_plan,
        (
            "database\tsource_column\toutput_trait\tvalue_type\taggregation\n"
            "gift\tis_carnivorous\tis_carnivorous\tbinary\tany\n"
        ),
    )

    db_sources = tmp_path / "input" / "input_generation" / "trait_database_sources.tsv"
    write_text(
        db_sources,
        (
            "database\tacquisition_mode\turi\tspecies_column\tresponse_format\n"
            f"gift\tspecies_api\tfile://{api_root}/{{species}}.tsv\tspecies\ttsv\n"
        ),
    )

    output = tmp_path / "input" / "species_trait" / "species_trait.tsv"
    completed = run_script(
        "--download-manifest",
        str(manifest),
        "--trait-plan",
        str(trait_plan),
        "--database-sources",
        str(db_sources),
        "--downloads-dir",
        str(tmp_path / "downloads" / "trait_datasets"),
        "--output",
        str(output),
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout
    df = pandas.read_csv(output, sep="\t")
    assert df.shape[0] == 2
    assert df["is_carnivorous"].tolist() == [0, 0]


def test_generate_species_trait_normalizes_taxonomic_qualifiers_for_matching():
    mod = load_script_module()

    assert mod.normalize_species_name("Dictyostelium cf. discoideum") == "Dictyostelium_cf_discoideum"
    assert mod.normalize_species_name("Bacillus subtilis subsp. subtilis") == "Bacillus_subtilis_subsp_subtilis"
    assert mod.normalize_species_name("Amoeba sp. JDS-Ruffled") == "Amoeba_sp_JDSRuffled"
    assert mod.split_genus_epithet("Dictyostelium_cf_discoideum") == ("Dictyostelium", "discoideum")
    assert mod.split_genus_epithet("Amoeba_sp_JDSRuffled") is None


def test_generate_species_trait_maps_base_species_rows_to_unique_qualified_target(tmp_path):
    species_cds = tmp_path / "output" / "input_generation" / "species_cds"
    species_cds.mkdir(parents=True, exist_ok=True)
    (species_cds / "Dictyostelium_cf_discoideum.fa").write_text(">a\nATG\n", encoding="utf-8")

    trait_source = tmp_path / "downloads" / "trait_datasets" / "animaltraits" / "animaltraits.tsv"
    write_text(
        trait_source,
        (
            "scientific_name\tlifespan_years\n"
            "Dictyostelium discoideum\t7\n"
        ),
    )

    trait_plan = tmp_path / "input" / "input_generation" / "trait_plan.tsv"
    write_text(
        trait_plan,
        (
            "database\tsource_column\toutput_trait\tvalue_type\taggregation\n"
            "animaltraits\tlifespan_years\tlifespan_years\tnumeric\tmedian\n"
        ),
    )

    db_sources = tmp_path / "input" / "input_generation" / "trait_database_sources.tsv"
    write_text(
        db_sources,
        (
            "database\tacquisition_mode\turi\tspecies_column\tdelimiter\n"
            f"animaltraits\tbulk\t{trait_source}\tscientific_name\ttsv\n"
        ),
    )

    output = tmp_path / "species_trait.tsv"
    completed = run_script(
        "--species-source",
        "species_cds",
        "--species-cds-dir",
        str(species_cds),
        "--trait-plan",
        str(trait_plan),
        "--database-sources",
        str(db_sources),
        "--downloads-dir",
        str(tmp_path / "downloads" / "trait_datasets"),
        "--output",
        str(output),
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout

    df = pandas.read_csv(output, sep="\t")
    assert df.shape[0] == 1
    assert df.loc[0, "species"] == "Dictyostelium_cf_discoideum"
    assert df.loc[0, "lifespan_years"] == 7


def test_generate_species_trait_does_not_duplicate_base_rows_when_exact_target_exists(tmp_path):
    species_cds = tmp_path / "output" / "input_generation" / "species_cds"
    species_cds.mkdir(parents=True, exist_ok=True)
    (species_cds / "Dictyostelium_discoideum.fa").write_text(">a\nATG\n", encoding="utf-8")
    (species_cds / "Dictyostelium_cf_discoideum.fa").write_text(">b\nATG\n", encoding="utf-8")

    trait_source = tmp_path / "downloads" / "trait_datasets" / "animaltraits" / "animaltraits.tsv"
    write_text(
        trait_source,
        (
            "scientific_name\tlifespan_years\n"
            "Dictyostelium discoideum\t7\n"
        ),
    )

    trait_plan = tmp_path / "input" / "input_generation" / "trait_plan.tsv"
    write_text(
        trait_plan,
        (
            "database\tsource_column\toutput_trait\tvalue_type\taggregation\n"
            "animaltraits\tlifespan_years\tlifespan_years\tnumeric\tmedian\n"
        ),
    )

    db_sources = tmp_path / "input" / "input_generation" / "trait_database_sources.tsv"
    write_text(
        db_sources,
        (
            "database\tacquisition_mode\turi\tspecies_column\tdelimiter\n"
            f"animaltraits\tbulk\t{trait_source}\tscientific_name\ttsv\n"
        ),
    )

    output = tmp_path / "species_trait.tsv"
    completed = run_script(
        "--species-source",
        "species_cds",
        "--species-cds-dir",
        str(species_cds),
        "--trait-plan",
        str(trait_plan),
        "--database-sources",
        str(db_sources),
        "--downloads-dir",
        str(tmp_path / "downloads" / "trait_datasets"),
        "--output",
        str(output),
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout

    df = pandas.read_csv(output, sep="\t")
    value_map = dict(zip(df["species"], df["lifespan_years"]))
    assert value_map["Dictyostelium_discoideum"] == 7
    assert pandas.isna(value_map["Dictyostelium_cf_discoideum"])


def test_generate_species_trait_strict_fails_when_source_column_missing(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    write_text(
        manifest,
        (
            "provider\tid\tspecies_key\n"
            "local\tid\tArabidopsis_thaliana\n"
        ),
    )
    trait_source = tmp_path / "traits.tsv"
    write_text(
        trait_source,
        (
            "species\tleaf_n\n"
            "Arabidopsis_thaliana\t2.1\n"
        ),
    )
    trait_plan = tmp_path / "trait_plan.tsv"
    write_text(
        trait_plan,
        (
            "database\tsource_column\toutput_trait\n"
            "austraits\tmissing_column\tleaf_N_percent\n"
        ),
    )
    db_sources = tmp_path / "trait_database_sources.tsv"
    write_text(
        db_sources,
        (
            "database\tacquisition_mode\turi\tspecies_column\n"
            f"austraits\tbulk\t{trait_source}\tspecies\n"
        ),
    )

    completed = run_script(
        "--download-manifest",
        str(manifest),
        "--trait-plan",
        str(trait_plan),
        "--database-sources",
        str(db_sources),
        "--downloads-dir",
        str(tmp_path / "downloads" / "trait_datasets"),
        "--strict",
    )
    assert completed.returncode != 0
    assert "missing source column" in completed.stdout


def test_generate_species_trait_supports_long_format_with_trait_key_and_zip_source(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    write_text(
        manifest,
        (
            "provider\tid\tspecies_key\n"
            "local\ta\tEucalyptus_regnans\n"
            "local\tb\tAcacia_dealbata\n"
        ),
    )

    zip_path = tmp_path / "austraits.zip"
    zip_member = "austraits-7.0.0/traits.csv"
    with zipfile.ZipFile(zip_path, mode="w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(
            zip_member,
            (
                "taxon_name,trait_name,value\n"
                "Eucalyptus regnans,leaf_N_per_dry_mass,1.1\n"
                "Eucalyptus regnans,genome_size,1.7\n"
                "Acacia dealbata,leaf_N_per_dry_mass,2.2\n"
                "Acacia dealbata,genome_size,0.9\n"
            ),
        )

    trait_plan = tmp_path / "trait_plan.tsv"
    write_text(
        trait_plan,
        (
            "database\tsource_column\toutput_trait\tvalue_type\taggregation\ttrait_key\ttrait_key_column\n"
            "austraits\tvalue\tleaf_N_per_dry_mass\tnumeric\tmedian\tleaf_N_per_dry_mass\ttrait_name\n"
            "austraits\tvalue\tgenome_size\tnumeric\tmedian\tgenome_size\ttrait_name\n"
        ),
    )

    db_sources = tmp_path / "trait_database_sources.tsv"
    write_text(
        db_sources,
        (
            "database\tacquisition_mode\turi\tspecies_column\tdelimiter\tarchive_member\ttrait_key_column\n"
            f"austraits\tbulk\t{zip_path}\ttaxon_name\tcsv\t{zip_member}\ttrait_name\n"
        ),
    )

    output = tmp_path / "species_trait.tsv"
    completed = run_script(
        "--download-manifest",
        str(manifest),
        "--trait-plan",
        str(trait_plan),
        "--database-sources",
        str(db_sources),
        "--downloads-dir",
        str(tmp_path / "downloads" / "trait_datasets"),
        "--output",
        str(output),
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout
    df = pandas.read_csv(output, sep="\t")
    assert set(df.columns) == {"species", "leaf_N_per_dry_mass", "genome_size"}
    out = dict(zip(df["species"], zip(df["leaf_N_per_dry_mass"], df["genome_size"])))
    assert out["Eucalyptus_regnans"] == (1.1, 1.7)
    assert out["Acacia_dealbata"] == (2.2, 0.9)


def test_generate_species_trait_bulk_uri_accepts_multiple_sources(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    write_text(
        manifest,
        (
            "provider\tid\tspecies_key\n"
            "local\ta\tCanis_lupus\n"
            "local\tb\tFelis_catus\n"
        ),
    )
    source1 = tmp_path / "elton_1.tsv"
    source2 = tmp_path / "elton_2.tsv"
    write_text(source1, "Scientific\tBodyMass-Value\nCanis lupus\t30\n")
    write_text(source2, "Scientific\tBodyMass-Value\nFelis catus\t4\n")

    trait_plan = tmp_path / "trait_plan.tsv"
    write_text(
        trait_plan,
        "database\tsource_column\toutput_trait\tvalue_type\taggregation\neltontraits\tBodyMass-Value\tbody_mass_g\tnumeric\tmedian\n",
    )
    db_sources = tmp_path / "trait_database_sources.tsv"
    write_text(
        db_sources,
        (
            "database\tacquisition_mode\turi\tspecies_column\tdelimiter\n"
            f"eltontraits\tbulk\t{source1},{source2}\tScientific\ttab\n"
        ),
    )
    output = tmp_path / "species_trait.tsv"
    completed = run_script(
        "--download-manifest",
        str(manifest),
        "--trait-plan",
        str(trait_plan),
        "--database-sources",
        str(db_sources),
        "--downloads-dir",
        str(tmp_path / "downloads" / "trait_datasets"),
        "--output",
        str(output),
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout
    df = pandas.read_csv(output, sep="\t")
    vals = dict(zip(df["species"], df["body_mass_g"]))
    assert vals["Canis_lupus"] == 30
    assert vals["Felis_catus"] == 4


def test_generate_species_trait_gift_api_mode_uses_species_lookup_and_trait_pages(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    write_text(
        manifest,
        (
            "provider\tid\tspecies_key\n"
            "local\ta\tHomo_sapiens\n"
            "local\tb\tMus_musculus\n"
        ),
    )
    trait_plan = tmp_path / "trait_plan.tsv"
    write_text(
        trait_plan,
        (
            "database\tsource_column\toutput_trait\tvalue_type\taggregation\ttrait_key\ttrait_key_column\n"
            "gift\ttrait_value\twoodiness\tcategorical\tfirst\t1.1.1\ttrait_ID\n"
        ),
    )

    request_log = []

    class GiftHandler(BaseHTTPRequestHandler):
        def do_GET(self):  # noqa: N802
            parsed = urlparse(self.path)
            query_map = parse_qs(parsed.query)
            query_name = query_map.get("query", [""])[0]
            request_log.append(
                {
                    "path": parsed.path,
                    "query": query_name,
                    "genus": query_map.get("genus", [""])[0],
                    "epithet": query_map.get("epithet", [""])[0],
                    "traitid": query_map.get("traitid", [""])[0],
                    "startat": query_map.get("startat", [""])[0],
                }
            )

            if parsed.path != "/api/index1.0.php":
                payload = []
            elif query_name == "names_matched_unique":
                genus = query_map.get("genus", [""])[0]
                epithet = query_map.get("epithet", [""])[0]
                if (genus, epithet) == ("Homo", "sapiens"):
                    payload = [{"work_ID": "101", "work_species": "Homo sapiens", "accepted": "1"}]
                elif (genus, epithet) == ("Mus", "musculus"):
                    payload = [{"work_ID": "202", "work_species": "Mus musculus", "accepted": "1"}]
                else:
                    payload = []
            elif query_name == "traits":
                trait_id = query_map.get("traitid", [""])[0]
                start_at = int(query_map.get("startat", ["0"])[0] or "0")
                if trait_id != "1.1.1":
                    payload = []
                elif start_at == 0:
                    payload = [
                        {"work_ID": "101", "trait_value": "woody", "agreement": "1"},
                        {"work_ID": "999", "trait_value": "woody", "agreement": "1"},
                    ]
                elif start_at == 2:
                    payload = [
                        {"work_ID": "202", "trait_value": "non-woody", "agreement": "1"},
                    ]
                else:
                    payload = []
            else:
                payload = []

            body = json.dumps(payload).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, fmt, *args):  # noqa: A003
            return

    server = ThreadingHTTPServer(("127.0.0.1", 0), GiftHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    completed = None
    try:
        db_sources = tmp_path / "trait_database_sources.tsv"
        write_text(
            db_sources,
            (
                "database\tacquisition_mode\turi\tspecies_column\tgift_version\tgift_page_size\n"
                f"gift\tgift_api\thttp://127.0.0.1:{server.server_port}/api/\tspecies\t1.0\t2\n"
            ),
        )
        output = tmp_path / "species_trait.tsv"
        completed = run_script(
            "--download-manifest",
            str(manifest),
            "--trait-plan",
            str(trait_plan),
            "--database-sources",
            str(db_sources),
            "--downloads-dir",
            str(tmp_path / "downloads" / "trait_datasets"),
            "--output",
            str(output),
        )
    finally:
        server.shutdown()
        thread.join(timeout=3)
        server.server_close()

    assert completed is not None
    assert completed.returncode == 0, completed.stderr + completed.stdout
    df = pandas.read_csv(output, sep="\t")
    values = dict(zip(df["species"], df["woodiness"]))
    assert values["Homo_sapiens"] == "woody"
    assert values["Mus_musculus"] == "non-woody"
    lookup_calls = [row for row in request_log if row["query"] == "names_matched_unique"]
    assert {(row["genus"], row["epithet"]) for row in lookup_calls} == {("Homo", "sapiens"), ("Mus", "musculus")}
    trait_calls = [row for row in request_log if row["query"] == "traits"]
    assert any(row["startat"] == "0" and row["traitid"] == "1.1.1" for row in trait_calls)
    assert any(row["startat"] == "2" and row["traitid"] == "1.1.1" for row in trait_calls)


def test_generate_species_trait_non_strict_warns_and_continues_on_source_load_failure(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    write_text(
        manifest,
        (
            "provider\tid\tspecies_key\n"
            "local\ta\tHomo_sapiens\n"
        ),
    )
    trait_plan = tmp_path / "trait_plan.tsv"
    write_text(
        trait_plan,
        (
            "database\tsource_column\toutput_trait\tvalue_type\taggregation\ttrait_key\ttrait_key_column\n"
            "gift\ttrait_value\twoodiness\tcategorical\tmode\t1.1.1\ttrait_ID\n"
        ),
    )
    db_sources = tmp_path / "trait_database_sources.tsv"
    write_text(
        db_sources,
        (
            "database\tacquisition_mode\turi\tspecies_column\tgift_version\n"
            "gift\tgift_api\thttp://127.0.0.1:9/api/\tspecies\t1.0\n"
        ),
    )
    output = tmp_path / "species_trait.tsv"

    completed = run_script(
        "--download-manifest",
        str(manifest),
        "--trait-plan",
        str(trait_plan),
        "--database-sources",
        str(db_sources),
        "--downloads-dir",
        str(tmp_path / "downloads" / "trait_datasets"),
        "--output",
        str(output),
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout
    assert "WARNING: [gift] failed to load source:" in completed.stdout
    df = pandas.read_csv(output, sep="\t")
    assert list(df.columns) == ["species"]
    assert df["species"].tolist() == ["Homo_sapiens"]


def test_generate_species_trait_gift_api_resolves_trait_name_via_traits_meta(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    write_text(
        manifest,
        (
            "provider\tid\tspecies_key\n"
            "local\ta\tHomo_sapiens\n"
        ),
    )
    trait_plan = tmp_path / "trait_plan.tsv"
    write_text(
        trait_plan,
        (
            "database\tsource_column\toutput_trait\tvalue_type\taggregation\ttrait_key\ttrait_key_column\n"
            "gift\ttrait_value\twoodiness\tcategorical\tfirst\tWoodiness_1\ttrait_ID\n"
        ),
    )
    request_log = []

    class GiftHandler(BaseHTTPRequestHandler):
        def do_GET(self):  # noqa: N802
            parsed = urlparse(self.path)
            query_map = parse_qs(parsed.query)
            query_name = query_map.get("query", [""])[0]
            request_log.append(
                {
                    "query": query_name,
                    "traitid": query_map.get("traitid", [""])[0],
                }
            )
            if parsed.path != "/api/index1.0.php":
                payload = []
            elif query_name == "names_matched_unique":
                payload = [{"work_ID": "101", "work_species": "Homo sapiens", "accepted": "1"}]
            elif query_name == "traits_meta":
                payload = [
                    {"Lvl3": "1.1.1", "Trait1": "Woodiness", "Trait2": "Woodiness_1", "count": "200"},
                    {"Lvl3": "1.1.2", "Trait1": "Woodiness", "Trait2": "Woodiness_2", "count": "100"},
                ]
            elif query_name == "traits":
                if query_map.get("traitid", [""])[0] == "1.1.1":
                    payload = [{"work_ID": "101", "trait_value": "woody", "agreement": "1"}]
                else:
                    payload = []
            else:
                payload = []
            body = json.dumps(payload).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, fmt, *args):  # noqa: A003
            return

    server = ThreadingHTTPServer(("127.0.0.1", 0), GiftHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    completed = None
    try:
        db_sources = tmp_path / "trait_database_sources.tsv"
        write_text(
            db_sources,
            (
                "database\tacquisition_mode\turi\tspecies_column\tgift_version\n"
                f"gift\tgift_api\thttp://127.0.0.1:{server.server_port}/api/\tspecies\t1.0\n"
            ),
        )
        output = tmp_path / "species_trait.tsv"
        completed = run_script(
            "--download-manifest",
            str(manifest),
            "--trait-plan",
            str(trait_plan),
            "--database-sources",
            str(db_sources),
            "--downloads-dir",
            str(tmp_path / "downloads" / "trait_datasets"),
            "--output",
            str(output),
        )
    finally:
        server.shutdown()
        thread.join(timeout=3)
        server.server_close()

    assert completed is not None
    assert completed.returncode == 0, completed.stderr + completed.stdout
    df = pandas.read_csv(output, sep="\t")
    assert dict(zip(df["species"], df["woodiness"]))["Homo_sapiens"] == "woody"
    assert any(row["query"] == "traits_meta" for row in request_log)
    assert any(row["query"] == "traits" and row["traitid"] == "1.1.1" for row in request_log)


def test_generate_species_trait_print_gift_traits_supports_search_and_limit():
    class GiftHandler(BaseHTTPRequestHandler):
        def do_GET(self):  # noqa: N802
            parsed = urlparse(self.path)
            query_map = parse_qs(parsed.query)
            query_name = query_map.get("query", [""])[0]
            if parsed.path == "/api/index.php" and query_name == "versions":
                payload = [{"version": "1.0"}]
            elif parsed.path == "/api/index1.0.php" and query_name == "traits_meta":
                payload = [
                    {"Lvl3": "1.1.1", "Trait2": "Woodiness_1", "Trait1": "Woodiness", "type": "categorical", "Units": "woody, non-woody", "count": "249169"},
                    {"Lvl3": "3.2.3", "Trait2": "Seed_mass_mean", "Trait1": "Seed mass", "type": "numeric", "Units": "g", "count": "24034"},
                ]
            else:
                payload = []
            body = json.dumps(payload).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, fmt, *args):  # noqa: A003
            return

    server = ThreadingHTTPServer(("127.0.0.1", 0), GiftHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    completed = None
    try:
        completed = run_script(
            "--print-gift-traits",
            "--gift-api",
            f"http://127.0.0.1:{server.server_port}/api/",
            "--gift-version",
            "latest",
            "--gift-trait-search",
            "wood",
            "--gift-trait-limit",
            "1",
        )
    finally:
        server.shutdown()
        thread.join(timeout=3)
        server.server_close()

    assert completed is not None
    assert completed.returncode == 0, completed.stderr + completed.stdout
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    assert lines[0] == "trait_id\ttrait_name\ttrait_group\tvalue_type\tunits\tcount"
    assert len(lines) == 2
    assert "1.1.1\tWoodiness_1\tWoodiness\tcategorical" in lines[1]


def test_generate_species_trait_print_supported_databases():
    completed = run_script("--print-supported-databases")
    assert completed.returncode == 0, completed.stderr
    assert "database\tacquisition_mode\tretrieval_scope\tnotes" in completed.stdout
    assert "eltontraits\tbulk\tall_species_snapshot" in completed.stdout
    assert "gift\tgift_api\ttrait_subset_api" in completed.stdout


def test_aggregate_trait_column_mode_prefers_sorted_value_on_tie():
    module = load_script_module()
    plan_row = module.TraitPlanRow(
        database="gift",
        source_column="trait_value",
        output_trait="woodiness",
        value_type="categorical",
        aggregation="mode",
        positive_values=set(),
        trait_key="1.1.1",
        trait_key_column="trait_ID",
    )
    db_df = pandas.DataFrame(
        {
            "__species_norm": ["Homo_sapiens", "Homo_sapiens", "Homo_sapiens", "Mus_musculus"],
            "trait_value": ["woody", "non-woody", "woody", "non-woody"],
        }
    )
    aggregated = module.aggregate_trait_column(db_df=db_df, plan_row=plan_row)
    assert aggregated.to_dict() == {
        "Homo_sapiens": "woody",
        "Mus_musculus": "non-woody",
    }


def test_aggregate_trait_column_binary_positive_values():
    module = load_script_module()
    plan_row = module.TraitPlanRow(
        database="gift",
        source_column="is_carnivorous",
        output_trait="is_carnivorous",
        value_type="binary",
        aggregation="any",
        positive_values={"yes", "present"},
        trait_key="",
        trait_key_column="",
    )
    db_df = pandas.DataFrame(
        {
            "__species_norm": ["Homo_sapiens", "Homo_sapiens", "Mus_musculus", "Mus_musculus"],
            "is_carnivorous": ["", "yes", "absent", "0"],
        }
    )
    aggregated = module.aggregate_trait_column(db_df=db_df, plan_row=plan_row)
    assert aggregated.to_dict() == {
        "Homo_sapiens": 1,
        "Mus_musculus": 0,
    }
