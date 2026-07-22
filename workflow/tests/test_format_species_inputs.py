import csv
import gzip
import io
import json
import os
import shutil
import sqlite3
import ssl
import subprocess
import sys
import tarfile

import pytest
from http.client import RemoteDisconnected
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from urllib.error import URLError

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "format_species_inputs.py"
VALIDATE_MAPPING_SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "validate_cds_gff_mapping.py"
SMALL_DATASET_ROOT = Path(__file__).resolve().parent / "data" / "small_gfe_dataset"


def load_module():
    spec = spec_from_file_location("format_species_inputs_module", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_script(*args, env=None):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), *args],
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )


def run_validate_mapping_script(*args):
    return subprocess.run(
        [sys.executable, str(VALIDATE_MAPPING_SCRIPT_PATH), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def test_cli_default_outputs_honor_isolated_output_root(monkeypatch, tmp_path):
    module = load_module()
    output_root = tmp_path / "isolated-output"
    monkeypatch.setenv("GG_INPUT_GENERATION_OUTPUT_ROOT", str(output_root))

    args = module.build_arg_parser().parse_args(["--provider", "local"])

    assert Path(args.species_cds_dir) == output_root / "species_cds"
    assert Path(args.species_gff_dir) == output_root / "species_gff"
    assert Path(args.species_genome_dir) == output_root / "species_genome"
    assert Path(args.species_summary_output) == output_root / "gg_input_generation_species.tsv"


def test_provider_resolvers_share_one_registry_interface():
    module = load_module()

    assert "coge" in module.PROVIDER_ID_RESOLVERS
    assert "figshare" in module.PROVIDER_ID_RESOLVERS
    assert all(hasattr(resolver, "resolve") for resolver in module.PROVIDER_ID_RESOLVERS.values())


def test_transient_network_error_treats_ssl_eof_as_retryable():
    module = load_module()

    error = URLError(ssl.SSLEOFError(8, "UNEXPECTED_EOF_WHILE_READING"))

    assert module.is_transient_network_error(error)


def test_cds_extension_is_treated_as_fasta():
    module = load_module()

    assert module.is_fasta_filename("GZX_Primary.cds")
    assert module.is_fasta_filename("GZX_Primary.cds.gz")
    assert not module.is_fasta_filename("GZX_Primary.gff")
    assert (
        module.normalize_cds_output_basename("GZX_Primary.cds", "Fakus_species")
        == "Fakus_species_GZX_Primary.cds.fa.gz"
    )
    assert (
        module.normalize_cds_output_basename("GZX_Primary.cds.gz", "Fakus_species")
        == "Fakus_species_GZX_Primary.cds.fa.gz"
    )


def test_coge_cds_header_uses_delimited_feature_id():
    module = load_module()
    headers = (
        "Gilia yorkii||Gy1||47683||49946||GY000001-RA||-1||CDS||3560718418||1",
        "Gilia yorkii||Gy1||92931||101228||GY000002-RA||-1||CDS||3560718421||2",
        (
            "Sarracenia purpurea ||CM117427.1_RagTag||149705||152539||"
            "evm.model.CM117427.1_RagTag.1||-1||CDS||4559620697||1"
        ),
    )

    assert [module.extract_provider_id("coge", header) for header in headers] == [
        "GY000001-RA",
        "GY000002-RA",
        "evm.model.CM117427.1_RagTag.1",
    ]
    gilia_task = {"provider": "coge", "species_prefix": "Gilia_yorkii"}
    assert module.build_formatted_cds_id(gilia_task, headers[0]) == "Gilia_yorkii_GY000001-RA"
    assert module.build_formatted_cds_id(gilia_task, headers[1]) == "Gilia_yorkii_GY000002-RA"


def test_coge_coding_export_requires_cds_features(tmp_path):
    module = load_module()
    gff_path = tmp_path / "gene_only.gff"
    gff_path.write_text("chr1\tCoGe\tgene\t1\t9\t.\t+\t.\tID=g1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="coding-only GFF export with no CDS records"):
        module.validate_coge_export_gff_file(gff_path, gid="69349")


def test_format_coge_delimited_cds_headers_map_to_gff_genes(tmp_path):
    input_dir = tmp_path / "CoGe" / "species_wise_original"
    species_dir = input_dir / "Gilia_yorkii"
    species_dir.mkdir(parents=True)
    (species_dir / "Gilia_yorkii.coge.gid62042.cds.fa").write_text(
        (
            ">Gilia yorkii||Gy1||1||9||GY000001-RA||1||CDS||101||1\nATGAAATTT\n"
            ">Gilia yorkii||Gy1||20||28||GY000002-RA||1||CDS||102||2\nATGCCCTTT\n"
        ),
        encoding="utf-8",
    )
    sarracenia_dir = input_dir / "Sarracenia_purpurea"
    sarracenia_dir.mkdir()
    (sarracenia_dir / "Sarracenia_purpurea.coge.gid69349.cds.fa").write_text(
        (
            ">Sarracenia purpurea ||CM117427.1_RagTag||1||9||"
            "evm.model.CM117427.1_RagTag.1||1||CDS||201||1\nATGAAATTT\n"
        ),
        encoding="utf-8",
    )
    (sarracenia_dir / "Sarracenia_purpurea.gid69349.gff").write_text(
        "\n".join(
            [
                "##gff-version\t3",
                (
                    "CM117427.1_RagTag\tCoGe\tgene\t1\t9\t.\t+\t.\t"
                    "ID=EVM%20prediction%20CM117427.1_RagTag.1;"
                    "Name=EVM%20prediction%20CM117427.1_RagTag.1;"
                    "gene=EVM%20prediction%20CM117427.1_RagTag.1"
                ),
                (
                    "CM117427.1_RagTag\tCoGe\tmRNA\t1\t9\t.\t+\t.\t"
                    "ID=evm.model.CM117427.1_RagTag.1;"
                    "Parent=EVM%20prediction%20CM117427.1_RagTag.1"
                ),
                (
                    "CM117427.1_RagTag\tCoGe\tCDS\t1\t9\t.\t+\t0\t"
                    "ID=evm.model.CM117427.1_RagTag.1.cds;"
                    "Parent=evm.model.CM117427.1_RagTag.1;"
                    "Alias=evm.model.CM117427.1_RagTag.1"
                ),
                "",
            ]
        ),
        encoding="utf-8",
    )
    (species_dir / "Gilia_yorkii.gid62042.gff").write_text(
        "\n".join(
            [
                "##gff-version\t3",
                "Gy1\tCoGe\tgene\t1\t9\t.\t+\t.\tID=GY000001;Name=GY000001;gene=GY000001",
                "Gy1\tCoGe\tmRNA\t1\t9\t.\t+\t.\tID=GY000001.mRNA1;Parent=GY000001;Alias=GY000001-RA",
                "Gy1\tCoGe\tCDS\t1\t9\t.\t+\t0\tID=GY000001-RA;Parent=GY000001.mRNA1;CDS=GY000001-RA",
                "Gy1\tCoGe\tgene\t20\t28\t.\t+\t.\tID=GY000002;Name=GY000002;gene=GY000002",
                "Gy1\tCoGe\tmRNA\t20\t28\t.\t+\t.\tID=GY000002.mRNA1;Parent=GY000002;Alias=GY000002-RA",
                "Gy1\tCoGe\tCDS\t20\t28\t.\t+\t0\tID=GY000002-RA;Parent=GY000002.mRNA1;CDS=GY000002-RA",
                "",
            ]
        ),
        encoding="utf-8",
    )

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    completed = run_script(
        "--provider",
        "coge",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(tmp_path / "species_genome"),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Gilia_yorkii_coge.gid62042.cds.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [">Gilia_yorkii_GY000001", ">Gilia_yorkii_GY000002"]
    sarracenia_cds = out_cds / "Sarracenia_purpurea_coge.gid69349.cds.fa.gz"
    with gzip.open(sarracenia_cds, "rt", encoding="utf-8") as handle:
        sarracenia_headers = [line.strip() for line in handle if line.startswith(">")]
    assert sarracenia_headers == [">Sarracenia_purpurea_EVM_prediction_CM117427.1_RagTag.1"]

    mapping = run_validate_mapping_script(
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
    )
    assert mapping.returncode == 0, mapping.stderr + "\n" + mapping.stdout
    assert "[Gilia_yorkii] CDS-to-GFF mapping OK: 2/2 IDs" in mapping.stdout
    assert "[Sarracenia_purpurea] CDS-to-GFF mapping OK: 1/1 IDs" in mapping.stdout


def test_local_gff_genome_lcl_prefix_is_resolved_when_deriving_cds(tmp_path):
    input_dir = tmp_path / "Local" / "species_wise_original"
    species_dir = input_dir / "Gilia_yorkii"
    species_dir.mkdir(parents=True)
    (species_dir / "Gilia_yorkii.genome.fa").write_text(
        ">lcl|Gy1\nATGAAATTT\n>lcl|Gy2\nATGCCCTTT\n",
        encoding="utf-8",
    )
    (species_dir / "Gilia_yorkii.gff").write_text(
        "Gy1\tCoGe\tgene\t1\t9\t.\t+\t.\tID=GY000001\n"
        "Gy1\tCoGe\tmRNA\t1\t9\t.\t+\t.\tID=GY000001.mRNA1;Parent=GY000001\n"
        "Gy1\tCoGe\tCDS\t1\t9\t.\t+\t0\tID=GY000001-RA;Parent=GY000001.mRNA1\n"
        "Gy2\tCoGe\tgene\t1\t9\t.\t+\t.\tID=GY000002\n"
        "Gy2\tCoGe\tmRNA\t1\t9\t.\t+\t.\tID=GY000002.mRNA1;Parent=GY000002\n"
        "Gy2\tCoGe\tCDS\t1\t9\t.\t+\t0\tID=GY000002-RA;Parent=GY000002.mRNA1\n",
        encoding="utf-8",
    )

    out_cds = tmp_path / "species_cds"
    completed = run_script(
        "--provider",
        "local",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(tmp_path / "species_gff"),
        "--species-genome-dir",
        str(tmp_path / "species_genome"),
    )

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    formatted_cds = next(out_cds.glob("*.fa.gz"))
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert ">Gilia_yorkii_GY000001" in text
    assert ">Gilia_yorkii_GY000002" in text
    assert "ATGAAATTT" in text
    assert "ATGCCCTTT" in text


def test_gff_genome_seqid_map_rejects_many_to_one_lcl_alias():
    module = load_module()

    with pytest.raises(ValueError, match="both resolve to FASTA sequence ID 'lcl\\|Gy1'"):
        module.build_gff_genome_seqid_map(
            {"lcl|Gy1": "ATGAAATTT"},
            {"Gy1", "lcl|Gy1"},
        )


def test_gff_genome_seqid_map_reports_unresolved_seqids():
    module = load_module()

    seqid_map, missing_seqids = module.build_gff_genome_seqid_map(
        {"lcl|Gy1": "ATGAAATTT", "Gy2": "ATGCCCTTT"},
        {"Gy1", "Gy2", "Gy3"},
    )

    assert seqid_map == {"Gy1": "lcl|Gy1", "Gy2": "Gy2"}
    assert missing_seqids == ("Gy3",)


def test_direct_ncbi_like_cds_header_uses_locus_tag_for_gene_grouping():
    module = load_module()
    task = {
        "provider": "direct",
        "species_prefix": "Vanilla_planifolia",
    }
    header = (
        "lcl|CM028150.1_cds_KAG0495310.1_1 "
        "[locus_tag=HPP92_000001] [protein_id=KAG0495310.1] [db_xref=NCBI_GP:KAG0495310.1]"
    )

    gene_id = module.build_gene_aggregate_id(
        task,
        header,
        "Vanilla_planifolia_lcl_CM028150.1_cds_KAG0495310.1_1",
    )

    assert gene_id == "Vanilla_planifolia_HPP92_000001"


def test_direct_species_discovery_treats_lone_chr_fasta_as_genome_when_gff_exists(tmp_path):
    module = load_module()
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Santalum_album"
    species_dir.mkdir(parents=True)
    (species_dir / "tanxiang.FINAL.chr.fa").write_text(">lcl|Chr01\nATGAAATTT\n", encoding="utf-8")
    (species_dir / "tanxiang.FINAL.chr_modified.gff").write_text(
        "Chr01\tGnomon\tmRNA\t1\t9\t.\t+\t.\tID=SA1G00001\n"
        "Chr01\tGnomon\tCDS\t1\t9\t.\t+\t0\tID=SA1G00001.cds1;Parent=SA1G00001\n",
        encoding="utf-8",
    )

    tasks, warnings, errors = module.discover_tasks("direct", input_dir)

    assert errors == []
    assert warnings == []
    assert len(tasks) == 1
    assert tasks[0]["cds_path"] is None
    assert tasks[0]["genome_path"].name == "tanxiang.FINAL.chr.fa"
    assert tasks[0]["gff_path"].name == "tanxiang.FINAL.chr_modified.gff"


def test_direct_species_discovery_keeps_lone_nonmatching_fasta_as_cds_when_gff_exists(tmp_path):
    module = load_module()
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Fakus_species"
    species_dir.mkdir(parents=True)
    (species_dir / "sequences.fa").write_text(">gene1\nATGAAA\n", encoding="utf-8")
    (species_dir / "annotation.gff3").write_text(
        "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cds1;Parent=gene1\n",
        encoding="utf-8",
    )

    tasks, warnings, errors = module.discover_tasks("direct", input_dir)

    assert errors == []
    assert warnings == []
    assert len(tasks) == 1
    assert tasks[0]["cds_path"].name == "sequences.fa"
    assert tasks[0]["genome_path"] is None
    assert tasks[0]["gff_path"].name == "annotation.gff3"


class FakeTextPipe:
    def __init__(self):
        self.parts = []

    def write(self, text):
        self.parts.append(text)
        return len(text)

    def close(self):
        return None

    def getvalue(self):
        return "".join(self.parts)


class FakeBinaryResponse:
    def __init__(self, payload):
        self._payload = payload

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def read(self):
        return self._payload


def write_test_taxonomy_fixture(tmp_path):
    db_path = tmp_path / "taxa.sqlite"
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute(
        "CREATE TABLE species (taxid INT PRIMARY KEY, parent INT, spname VARCHAR(50) COLLATE NOCASE, common VARCHAR(50) COLLATE NOCASE, rank VARCHAR(50), track TEXT)"
    )
    cur.execute("CREATE TABLE synonym (taxid INT,spname VARCHAR(50) COLLATE NOCASE, PRIMARY KEY (spname, taxid))")
    species_rows = [
        (1, 1, "root", "", "no rank", "1"),
        (131567, 1, "cellular organisms", "", "cellular root", "131567,1"),
        (2759, 131567, "Eukaryota", "", "domain", "2759,131567,1"),
        (33090, 2759, "Viridiplantae", "", "kingdom", "33090,2759,131567,1"),
        (242159, 33090, "Ostreococcus lucimarinus", "", "species", "242159,33090,2759,131567,1"),
        (5786, 2759, "Dictyostelium discoideum", "", "species", "5786,2759,131567,1"),
        (5911, 2759, "Tetrahymena thermophila", "", "species", "5911,2759,131567,1"),
    ]
    cur.executemany(
        "INSERT INTO species (taxid, parent, spname, common, rank, track) VALUES (?, ?, ?, ?, ?, ?)", species_rows
    )
    cur.execute("INSERT INTO synonym (taxid, spname) VALUES (?, ?)", (242159, "Ostreococcus_lucimarinus"))
    conn.commit()
    conn.close()

    def nodes_line(taxid, parent, rank, gc_id, mito_gc_id):
        return "{}\t|\t{}\t|\t{}\t|\t\t|\t0\t|\t0\t|\t{}\t|\t0\t|\t{}\t|\t0\t|\t0\t|\t0\t|\t\t|\n".format(
            taxid, parent, rank, gc_id, mito_gc_id
        )

    gencode_text = (
        "1\t|\tSGC0\t|\tStandard\t|\t\t|\t\t|\n"
        "4\t|\tSGC4\t|\tMold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma\t|\t\t|\t\t|\n"
        "6\t|\tSGC6\t|\tCiliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear\t|\t\t|\t\t|\n"
        "11\t|\tSGC11\t|\tBacterial, Archaeal and Plant Plastid\t|\t\t|\t\t|\n"
    )
    nodes_text = "".join(
        [
            nodes_line(1, 1, "no rank", 1, 0),
            nodes_line(131567, 1, "cellular root", 1, 0),
            nodes_line(2759, 131567, "domain", 1, 0),
            nodes_line(33090, 2759, "kingdom", 1, 0),
            nodes_line(242159, 33090, "species", 1, 1),
            nodes_line(5786, 2759, "species", 1, 4),
            nodes_line(5911, 2759, "species", 6, 4),
        ]
    )

    taxdump_path = tmp_path / "taxdump.tar.gz"
    with tarfile.open(taxdump_path, "w:gz") as archive:
        for name, text in {
            "gencode.dmp": gencode_text,
            "nodes.dmp": nodes_text,
            "readme.txt": "test fixture\n",
        }.items():
            payload = text.encode("utf-8")
            info = tarfile.TarInfo(name=name)
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))

    return db_path, taxdump_path


def test_format_species_inputs_with_small_fixture_all_providers(tmp_path):
    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    stats_json = tmp_path / "stats.json"

    completed = run_script(
        "--provider",
        "all",
        "--input-dir",
        str(SMALL_DATASET_ROOT),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--stats-output",
        str(stats_json),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    ensembl_cds = out_cds / "Ostreococcus_lucimarinus_ASM9206v1.cds.all.fa.gz"
    phycocosm_cds = out_cds / "Microglena_spYARC_MicrYARC1_GeneCatalog_CDS_20220803.fa.gz"
    phytozome_cds = out_cds / "Hydrocotyle_leucocephala_HleucocephalaHAP1_768_v2.1.cds_primaryTranscriptOnly.fa.gz"
    assert ensembl_cds.exists()
    assert phycocosm_cds.exists()
    assert phytozome_cds.exists()

    with gzip.open(ensembl_cds, "rt", encoding="utf-8") as handle:
        ensembl_text = handle.read()
    assert ensembl_text.count(">Ostreococcus_lucimarinus_OSTLU_25062") == 1
    assert ">Ostreococcus_lucimarinus_OSTLU_99999" in ensembl_text
    assert "ATGAAN" in ensembl_text

    with gzip.open(phycocosm_cds, "rt", encoding="utf-8") as handle:
        phycocosm_text = handle.read()
    assert ">Microglena_spYARC_mRNA.MigICE15955" in phycocosm_text
    assert ">Microglena_spYARC_mRNA.MigICE_00468" in phycocosm_text

    with gzip.open(phytozome_cds, "rt", encoding="utf-8") as handle:
        phytozome_text = handle.read()
    assert ">Hydrocotyle_leucocephala_HyleuH1.06G006800" in phytozome_text
    assert "ATGATGAN" in phytozome_text

    ensembl_gff = out_gff / "Ostreococcus_lucimarinus_ASM9206v1.56.gff.gz"
    phytozome_gff = out_gff / "Hydrocotyle_leucocephala_HleucocephalaHAP1_768_v2.1.gene.gff.gz"
    phytozome_exons = out_gff / "Hydrocotyle_leucocephala_HleucocephalaHAP1_768_v2.1.gene_exons.gff.gz"
    assert ensembl_gff.exists()
    assert phytozome_gff.exists()
    assert not phytozome_exons.exists()

    with gzip.open(ensembl_gff, "rt", encoding="utf-8") as handle:
        gff_text = handle.read()
    assert "evm.model." not in gff_text
    assert "Oropetium_20150105_" not in gff_text

    with gzip.open(phytozome_gff, "rt", encoding="utf-8") as handle:
        phytozome_gff_text = handle.read()
    assert "evm_27.model." not in phytozome_gff_text
    stats = json.loads(stats_json.read_text(encoding="utf-8"))
    assert stats["species_processed"] == 3
    assert stats["num_species_cds_files"] == 3
    assert stats["num_species_gff_files"] == 3
    assert stats["cds_sequences_before"] >= stats["cds_sequences_after"]
    assert stats["cds_first_sequence_name"] != ""


def test_species_taxonomy_metadata_resolver_supports_nonstandard_nuclear_codes(tmp_path):
    mod = load_module()
    db_path, taxdump_path = write_test_taxonomy_fixture(tmp_path)
    resolver = mod.SpeciesTaxonomyMetadataResolver(str(db_path), str(taxdump_path))

    tetrahymena = resolver.resolve("Tetrahymena thermophila")
    assert tetrahymena["taxid"] == "5911"
    assert tetrahymena["nuclear_genetic_code_id"] == "6"
    assert tetrahymena["nuclear_genetic_code_name"] == "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear"
    assert tetrahymena["mitochondrial_genetic_code_id"] == "4"
    assert tetrahymena["plastid_genetic_code_id"] == ""


def test_parse_species_key_candidate_preserves_taxonomic_qualifiers():
    mod = load_module()

    assert mod.parse_species_key_candidate("Dictyostelium cf. discoideum") == "Dictyostelium_cf_discoideum"
    assert mod.parse_species_key_candidate("Bacillus subtilis subsp. subtilis") == "Bacillus_subtilis_subsp_subtilis"
    assert mod.parse_species_key_candidate("Amoeba sp. JDS-Ruffled") == "Amoeba_sp_JDSRuffled"
    assert mod.parse_species_key_candidate("Amoeba sp.") == "Amoeba_sp_unknown"
    assert (
        mod.parse_species_key_candidate("Solanum lycopersicum cultivar Heinz 1706")
        == "Solanum_lycopersicum_cultivar_Heinz1706"
    )
    assert mod.parse_species_key_candidate("Escherichia coli serovar O157") == "Escherichia_coli_serovar_O157"
    assert mod.parse_species_key_candidate("Citrus x limon") == "Citrus_x_limon"
    assert (
        mod.parse_species_key_candidate("Cenchrus americanus x Cenchrus purpureus")
        == "Cenchrus_americanus_x_Cenchrus_purpureus"
    )
    assert (
        mod.parse_species_key_candidate("Cenchrus americanus \u00d7 Cenchrus purpureus")
        == "Cenchrus_americanus_x_Cenchrus_purpureus"
    )


def test_source_id_candidates_add_ensemblplants_species_and_accession_variants():
    mod = load_module()

    candidates = mod.source_id_candidates(
        "ensemblplants",
        "GCA_910589775.1",
        "Avena_eriantha",
    )
    assert "Avena_eriantha" in candidates
    assert "Avena_eriantha_gca910589775v1cm" in candidates

    candidates = mod.source_id_candidates(
        "ensemblplants",
        "GCA_000695525.1",
        "Brassica_oleracea_var._oleracea",
    )
    assert "Brassica_oleracea" in candidates

    candidates = mod.source_id_candidates(
        "ensemblplants",
        "GCA_001952365.2",
        "Oryza_sativa_aus_subgroup",
    )
    assert "Oryza_sativa_aus" in candidates


def test_source_id_candidates_add_ensemblmetazoa_accession_suffix_variants():
    mod = load_module()

    candidates = mod.source_id_candidates(
        "ensemblmetazoa",
        "GCA_003676215.3",
        "Rhopalosiphum_maidis",
    )

    assert "Rhopalosiphum_maidis" in candidates
    assert "Rhopalosiphum_maidis_gca003676215v3" in candidates


def test_expand_ensemblgenomes_id_candidates_prefers_matching_directory_suffix(monkeypatch):
    mod = load_module()
    provider_module = sys.modules[mod.expand_ensemblgenomes_id_candidates.__module__]

    monkeypatch.setattr(
        provider_module,
        "fetch_ensemblgenomes_dir_ids",
        lambda provider, timeout, headers: ["rhopalosiphum_maidis_gca003676215v3"],
    )

    candidates = mod.expand_ensemblgenomes_id_candidates(
        "ensemblmetazoa",
        ["Rhopalosiphum_maidis"],
        timeout=1,
        headers={},
    )

    assert candidates[0] == "rhopalosiphum_maidis_gca003676215v3"
    assert "Rhopalosiphum_maidis" in candidates


def test_infer_provider_from_ensemblgenomes_prefix_and_url():
    mod = load_module()

    assert mod.infer_provider_from_id("ensemblmetazoa:anopheles_gambiae") == "ensemblmetazoa"
    assert mod.infer_provider_from_id("ensemblprotists:phytophthora_parasitica") == "ensemblprotists"
    assert (
        mod.infer_provider_from_id(
            "https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/metazoa/fasta/anopheles_gambiae/cds/"
        )
        == "ensemblmetazoa"
    )
    assert (
        mod.infer_provider_from_id(
            "https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/protists/gff3/phytophthora_parasitica/"
        )
        == "ensemblprotists"
    )


def test_discover_ncbi_like_tasks_normalizes_bare_sp_species_key(tmp_path):
    mod = load_module()
    input_dir = tmp_path / "NCBI_Genome" / "species_wise_original"
    species_dir = input_dir / "Amoeba_sp"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "GCA_000000001.1_demo_cds_from_genomic.fna.gz").write_text("", encoding="utf-8")
    (species_dir / "GCA_000000001.1_demo_genomic.gff.gz").write_text("", encoding="utf-8")
    (species_dir / "GCA_000000001.1_demo_genomic.fna.gz").write_text("", encoding="utf-8")

    tasks, warnings, errors = mod.discover_ncbi_like_tasks(input_dir, "ncbi")

    assert warnings == []
    assert errors == []
    assert len(tasks) == 1
    assert tasks[0]["species_key"] == "Amoeba_sp_unknown"
    assert tasks[0]["species_prefix"] == "Amoeba_sp_unknown"


def test_species_prefix_keeps_hybrid_marker_with_following_epithet():
    mod = load_module()
    assert mod.species_prefix_from_value("Citrus_x_limon") == "Citrus_x_limon"
    assert mod.species_prefix_from_value("Petunia_x_hybrida") == "Petunia_x_hybrida"
    assert (
        mod.species_prefix_from_value("Cenchrus_americanus_x_Cenchrus_purpureus.cds.fa.gz")
        == "Cenchrus_americanus_x_Cenchrus_purpureus"
    )


def test_species_prefix_preserves_dotted_taxonomic_qualifiers():
    mod = load_module()

    assert (
        mod.species_prefix_from_value("Asimitellaria_furusei_var._furusei_demo.fa.gz")
        == "Asimitellaria_furusei_var._furusei"
    )
    assert (
        mod.species_prefix_from_value("Asimitellaria_furusei_var._subramosa.cds.fa.gz")
        == "Asimitellaria_furusei_var._subramosa"
    )
    assert mod.species_prefix_from_value("Arisaema_sp._aooni_demo.fa.gz") == "Arisaema_sp._aooni"
    assert (
        mod.species_prefix_from_value("Bacillus_subtilis_subsp._subtilis_demo.gff.gz")
        == "Bacillus_subtilis_subsp._subtilis"
    )
    assert mod.species_prefix_from_value("homo_sapiens.GRCh38.cds.all.fa.gz") == "homo_sapiens"


def test_discover_ensembl_like_tasks_preserves_dotted_taxonomic_qualifiers(tmp_path):
    mod = load_module()
    input_dir = tmp_path / "Ensembl" / "original_files"
    input_dir.mkdir(parents=True)
    for species_key in (
        "Asimitellaria_furusei_var._furusei",
        "Asimitellaria_furusei_var._subramosa",
    ):
        (input_dir / f"{species_key}.ASM.cds.all.fa.gz").write_text("", encoding="utf-8")
        (input_dir / f"{species_key}.ASM.1.gff3.gz").write_text("", encoding="utf-8")
        (input_dir / f"{species_key}.ASM.dna.primary_assembly.fa.gz").write_text("", encoding="utf-8")

    tasks, warnings, errors = mod.discover_ensembl_like_tasks(input_dir, "ensembl")

    assert warnings == []
    assert errors == []
    assert [task["species_key"] for task in tasks] == [
        "Asimitellaria_furusei_var._furusei",
        "Asimitellaria_furusei_var._subramosa",
    ]
    assert [task["species_prefix"] for task in tasks] == [
        "Asimitellaria_furusei_var._furusei",
        "Asimitellaria_furusei_var._subramosa",
    ]


def test_discover_ncbi_like_tasks_rejects_incomplete_taxonomic_qualifier_species_key(tmp_path):
    mod = load_module()
    input_dir = tmp_path / "NCBI_Genome" / "species_wise_original"
    species_dir = input_dir / "Dictyostelium_cf"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "GCA_054859205.1_ASM5485920v1_cds_from_genomic.fna.gz").write_text("", encoding="utf-8")
    (species_dir / "GCA_054859205.1_ASM5485920v1_genomic.gff.gz").write_text("", encoding="utf-8")
    (species_dir / "GCA_054859205.1_ASM5485920v1_genomic.fna.gz").write_text("", encoding="utf-8")

    tasks, warnings, errors = mod.discover_ncbi_like_tasks(input_dir, "ncbi")

    assert tasks == []
    assert warnings == []
    assert len(errors) == 1
    assert "requires a following epithet" in errors[0]
    assert "Dictyostelium_cf" in errors[0]


def test_species_taxonomy_metadata_resolver_falls_back_from_qualified_name_to_base_species(tmp_path):
    mod = load_module()
    db_path, taxdump_path = write_test_taxonomy_fixture(tmp_path)
    resolver = mod.SpeciesTaxonomyMetadataResolver(str(db_path), str(taxdump_path))

    dicty = resolver.resolve("Dictyostelium_cf_discoideum")
    assert dicty["taxid"] == "5786"
    assert dicty["nuclear_genetic_code_id"] == "1"
    assert dicty["mitochondrial_genetic_code_id"] == "4"


def test_manifest_declared_providers_preserves_manifest_order_and_skips_unlisted_rows():
    mod = load_module()

    providers = mod.manifest_declared_providers(
        [
            {"provider": "ncbi", "id": "GCA_000000001.1"},
            {"provider": "direct", "id": "sample_direct"},
            {"provider": "NCBI", "id": "GCA_000000002.1"},
            {"provider": "oryza_minuta", "id": "gramene_tetraploids"},
            {"provider": "", "id": "missing_provider"},
            {"provider": "unsupported", "id": "unsupported_provider"},
            {"provider": "local", "id": "sample_local"},
        ],
        provider_filter="all",
    )

    assert providers == ["ncbi", "direct", "oryza_minuta", "local"]


def test_species_summary_includes_taxid_and_genetic_codes_when_taxonomy_cache_is_available(tmp_path):
    db_path, taxdump_path = write_test_taxonomy_fixture(tmp_path)
    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    species_summary = tmp_path / "gg_input_generation_species.tsv"
    env = dict(os.environ)
    env["GG_TAXONOMY_DBFILE"] = str(db_path)
    env["GG_TAXONOMY_TAXDUMPFILE"] = str(taxdump_path)

    completed = run_script(
        "--provider",
        "ensemblplants",
        "--input-dir",
        str(SMALL_DATASET_ROOT / "20230216_EnsemblPlants" / "original_files"),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
        env=env,
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert row["species_prefix"] == "Ostreococcus_lucimarinus"
    assert row["taxid"] == "242159"
    assert row["nuclear_genetic_code_id"] == "1"
    assert row["nuclear_genetic_code_name"] == "Standard"
    assert row["mitochondrial_genetic_code_id"] == "1"
    assert row["mitochondrial_genetic_code_name"] == "Standard"
    assert row["plastid_genetic_code_id"] == "11"
    assert row["plastid_genetic_code_name"] == "Bacterial, Archaeal and Plant Plastid"


def test_format_species_inputs_strict_mode_accepts_cds_only_inputs(tmp_path):
    dataset_copy = tmp_path / "small_gfe_dataset_copy"
    shutil.copytree(SMALL_DATASET_ROOT, dataset_copy)
    missing_gff = (
        dataset_copy / "20230216_EnsemblPlants" / "original_files" / "Ostreococcus_lucimarinus.ASM9206v1.56.gff3"
    )
    missing_gff.unlink()

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    species_summary = tmp_path / "gg_input_generation_species.tsv"
    completed = run_script(
        "--provider",
        "ensemblplants",
        "--input-dir",
        str(dataset_copy / "20230216_EnsemblPlants" / "original_files"),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-summary-output",
        str(species_summary),
        "--strict",
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert any(path.name.endswith(".fa.gz") for path in out_cds.iterdir())
    assert not any(out_gff.iterdir())
    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["gff_status"] == "missing"
    assert rows[0]["genome_status"] == "missing"


def test_format_species_inputs_derives_cds_from_gff_and_genome_when_cds_is_missing(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Arabidopsis_thaliana"
    species_dir.mkdir(parents=True, exist_ok=True)
    gff_path = species_dir / "Arabidopsis_thaliana.annotation.gff3"
    genome_path = species_dir / "Arabidopsis_thaliana.genome.fa"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                "chr1\tsrc\tCDS\t7\t9\t.\t+\t0\tID=cds2;Parent=gene1.t1",
                "chr2\tsrc\tgene\t1\t9\t.\t-\t.\tID=gene2",
                "chr2\tsrc\tmRNA\t1\t9\t.\t-\t.\tID=gene2.t1;Parent=gene2",
                "chr2\tsrc\tCDS\t1\t3\t.\t-\t0\tID=cds3;Parent=gene2.t1",
                "chr2\tsrc\tCDS\t7\t9\t.\t-\t0\tID=cds4;Parent=gene2.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(
        ">chr1\nATGAAATTT\n>chr2\nTTTAAACAT\n",
        encoding="utf-8",
    )

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    species_summary = tmp_path / "gg_input_generation_species.tsv"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Arabidopsis_thaliana_annotation.derived.cds.fa.gz"
    assert formatted_cds.exists()
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert text.count(">Arabidopsis_thaliana_gene1") == 1
    assert text.count(">Arabidopsis_thaliana_gene2") == 1
    assert "ATGTTT" in text
    assert "ATGAAA" in text

    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert str(gff_path) in row["cds_input_path"]
    assert str(genome_path) in row["cds_input_path"]
    assert "derived CDS" in row["cds_input_path"]
    assert row["gff_input_path"] == str(gff_path)
    assert row["genome_input_path"] == str(genome_path)


def test_format_species_inputs_skips_mixed_strand_transcript_when_deriving_cds(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Ophrys_sphegodes"
    species_dir.mkdir(parents=True, exist_ok=True)
    gff_path = species_dir / "Ophrys_sphegodes.annotation.gff3"
    genome_path = species_dir / "Ophrys_sphegodes.genome.fa"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tID=gene1",
                "chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                "chr1\tsrc\tCDS\t10\t12\t.\t-\t0\tID=cds2;Parent=gene1.t1",
                "chr1\tsrc\tgene\t20\t25\t.\t+\t.\tID=gene2",
                "chr1\tsrc\tmRNA\t20\t25\t.\t+\t.\tID=gene2.t1;Parent=gene2",
                "chr1\tsrc\tCDS\t20\t25\t.\t+\t0\tID=cds3;Parent=gene2.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(">chr1\nATGAAATTTCCCGGGAAATGAAATAA\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert "skipping transcript 'gene1.t1'" in completed.stderr

    formatted_cds = out_cds / "Ophrys_sphegodes_annotation.derived.cds.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert ">Ophrys_sphegodes_gene1" not in text
    assert ">Ophrys_sphegodes_gene2" in text


def test_format_species_inputs_skips_mixed_strand_transcript_before_utr_trimming(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Ophrys_sphegodes"
    species_dir.mkdir(parents=True, exist_ok=True)
    gff_path = species_dir / "Ophrys_sphegodes.annotation.gff3"
    genome_path = species_dir / "Ophrys_sphegodes.genome.fa"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tID=gene1",
                "chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                "chr1\tsrc\tfive_prime_UTR\t1\t3\t.\t+\t.\tParent=gene1.t1",
                "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                "chr1\tsrc\tCDS\t10\t12\t.\t-\t0\tID=cds2;Parent=gene1.t1",
                "chr1\tsrc\tgene\t20\t25\t.\t+\t.\tID=gene2",
                "chr1\tsrc\tmRNA\t20\t25\t.\t+\t.\tID=gene2.t1;Parent=gene2",
                "chr1\tsrc\tCDS\t20\t25\t.\t+\t0\tID=cds3;Parent=gene2.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(">chr1\nATGAAATTTCCCGGGAAATGAAATAA\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert "skipping transcript 'gene1.t1'" in completed.stderr

    formatted_cds = out_cds / "Ophrys_sphegodes_annotation.derived.cds.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert ">Ophrys_sphegodes_gene1" not in text
    assert ">Ophrys_sphegodes_gene2" in text


def test_format_species_inputs_treats_softmasked_direct_fasta_as_genome(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Arabidopsis_thaliana"
    species_dir.mkdir(parents=True, exist_ok=True)
    gff_path = species_dir / "Arabidopsis_thaliana.genes.gff3"
    genome_path = species_dir / "Arabidopsis_thaliana.softmasked.fa"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(">chr1\nATGAAATTT\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Arabidopsis_thaliana_genes.derived.cds.fa.gz"
    formatted_genome = out_genome / "Arabidopsis_thaliana_softmasked.fa.gz"
    assert formatted_cds.exists()
    assert formatted_genome.exists()
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert ">Arabidopsis_thaliana_gene1" in text
    assert "ATGAAATTT" in text


def test_format_species_inputs_derives_cds_without_trimming_nonzero_phase(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Arabidopsis_thaliana"
    species_dir.mkdir(parents=True, exist_ok=True)
    gff_path = species_dir / "Arabidopsis_thaliana.annotation.gff3"
    genome_path = species_dir / "Arabidopsis_thaliana.genome.fa"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t15\t.\t+\t.\tID=gene1",
                "chr1\tsrc\tmRNA\t1\t15\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                "chr1\tsrc\tCDS\t10\t15\t.\t+\t2\tID=cds2;Parent=gene1.t1",
                "chr2\tsrc\tgene\t1\t15\t.\t-\t.\tID=gene2",
                "chr2\tsrc\tmRNA\t1\t15\t.\t-\t.\tID=gene2.t1;Parent=gene2",
                "chr2\tsrc\tCDS\t1\t6\t.\t-\t1\tID=cds3;Parent=gene2.t1",
                "chr2\tsrc\tCDS\t10\t15\t.\t-\t0\tID=cds4;Parent=gene2.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(
        ">chr1\nATGAAACCCGGGTTT\n>chr2\nAAACCCGGGTTTCAT\n",
        encoding="utf-8",
    )

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Arabidopsis_thaliana_annotation.derived.cds.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert ">Arabidopsis_thaliana_gene1" in text
    assert ">Arabidopsis_thaliana_gene2" in text
    assert "ATGAAAGGGTTT" in text
    assert "ATGAAATTT" not in text


def test_build_formatted_cds_id_ensembl_uses_transcript_token_for_tie_breaks():
    module = load_module()
    task = {
        "provider": "ensemblplants",
        "species_prefix": "Arabidopsis_thaliana",
    }
    raw_header = "AT1G01520.5 cds chromosome:TAIR10:1:190478:192436:1 gene:AT1G01520 gene_symbol:ASG4"
    derived_header = "transcript:AT1G01520.1 gene=AT1G01520"

    assert module.build_formatted_cds_id(task, raw_header) == "Arabidopsis_thaliana_AT1G01520.5"
    assert module.build_formatted_cds_id(task, derived_header) == "Arabidopsis_thaliana_AT1G01520.1"
    assert (
        module.build_gene_aggregate_id(task, raw_header, "Arabidopsis_thaliana_AT1G01520.5")
        == "Arabidopsis_thaliana_AT1G01520"
    )
    assert (
        module.build_gene_aggregate_id(task, derived_header, "Arabidopsis_thaliana_AT1G01520.1")
        == "Arabidopsis_thaliana_AT1G01520"
    )


def test_format_species_inputs_derives_cds_excluding_overlapping_utrs(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Arabidopsis_thaliana"
    species_dir.mkdir(parents=True, exist_ok=True)
    gff_path = species_dir / "Arabidopsis_thaliana.annotation.gff3"
    genome_path = species_dir / "Arabidopsis_thaliana.genome.fa"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t15\t.\t+\t.\tID=gene1",
                "chr1\tsrc\tmRNA\t1\t15\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                "chr1\tsrc\tfive_prime_UTR\t1\t6\t.\t+\t.\tParent=gene1.t1",
                "chr1\tsrc\tCDS\t1\t15\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                "chr2\tsrc\tgene\t1\t15\t.\t-\t.\tID=gene2",
                "chr2\tsrc\tmRNA\t1\t15\t.\t-\t.\tID=gene2.t1;Parent=gene2",
                "chr2\tsrc\tfive_prime_UTR\t10\t15\t.\t-\t.\tParent=gene2.t1",
                "chr2\tsrc\tCDS\t1\t15\t.\t-\t0\tID=cds2;Parent=gene2.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(
        ">chr1\nATGAAACCCGGGTTT\n>chr2\nAAACCCGGGTTTCAT\n",
        encoding="utf-8",
    )

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Arabidopsis_thaliana_annotation.derived.cds.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert ">Arabidopsis_thaliana_gene1" in text
    assert ">Arabidopsis_thaliana_gene2" in text
    assert "CCCGGGTTT" in text
    assert "ATGAAACCCGGGTTT" not in text


def test_format_species_inputs_rescue_overlap_merges_misassigned_gene_ids(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Arabidopsis_thaliana"
    species_dir.mkdir(parents=True, exist_ok=True)
    gff_path = species_dir / "Arabidopsis_thaliana.annotation.gff3"
    genome_path = species_dir / "Arabidopsis_thaliana.genome.fa"
    species_summary = tmp_path / "gg_input_generation_species.tsv"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=badGeneA",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=locusX.t1;Parent=badGeneA",
                "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cds1;Parent=locusX.t1",
                "chr1\tsrc\tCDS\t13\t18\t.\t+\t0\tID=cds2;Parent=locusX.t1",
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=badGeneB",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=locusX.t2;Parent=badGeneB",
                "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=cds3;Parent=locusX.t2",
                "chr1\tsrc\tCDS\t13\t18\t.\t+\t0\tID=cds4;Parent=locusX.t2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(">chr1\nATGAAACCCGGGTTTCCC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Arabidopsis_thaliana_annotation.derived.cds.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [">Arabidopsis_thaliana_locusX"]

    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["gene_grouping_mode"] == "rescue_overlap"
    assert rows[0]["aggregated_cds_removed"] == "1"


def test_format_species_inputs_uses_locus_tag_for_genbank_style_ncbi_cds(tmp_path):
    input_dir = tmp_path / "NCBI_Genome" / "species_wise_original"
    species_dir = input_dir / "Dictyostelium_cf_discoideum"
    species_dir.mkdir(parents=True, exist_ok=True)
    cds_path = species_dir / "GCA_054859205.1_ASM5485920v1_cds_from_genomic.fna.gz"
    gff_path = species_dir / "GCA_054859205.1_ASM5485920v1_genomic.gff.gz"
    genome_path = species_dir / "GCA_054859205.1_ASM5485920v1_genomic.fna.gz"
    with gzip.open(cds_path, "wt", encoding="utf-8") as handle:
        handle.write(
            (
                ">lcl|JBTAPH010000036.1_cds_KAM9986187.1_1 [locus_tag=ACTFIY_010592] [protein=hypothetical protein] [protein_id=KAM9986187.1] [location=complement(join(9..158,229..966))] [gbkey=CDS]\n"
                "ATGTCTACCACTGTTAACAATAATGATGCCTCTAGTAGTAGTAGCTCTGCCTCTAATAACGATGAATCCTTTGATTTAAGAATGAAATCAATGGAGGATCAAATCAATAACCTTTCATTAGCCTTTACCAGATTCATGAAAGAACCTATGTTCTCTTCTAATACCAAATCACGTAGCCAACCTTCTCATGATAACTCTGACACTGAGAATGAACAAAGTGATGACGAATCAAGTAACAAT\n"
                ">lcl|JBTAPH010000036.1_cds_KAM9986188.1_2 [locus_tag=ACTFIY_010593] [protein=hypothetical protein] [protein_id=KAM9986188.1] [location=complement(join(2114..2192,2246..2381,2478..2691))] [gbkey=CDS]\n"
                "ATGGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATC\n"
            )
        )
    with gzip.open(gff_path, "wt", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "##gff-version 3",
                    "JBTAPH010000036.1\tGenbank\tgene\t9\t966\t.\t-\t.\tID=gene-ACTFIY_010592;Name=ACTFIY_010592;gbkey=Gene;gene_biotype=protein_coding;locus_tag=ACTFIY_010592",
                    "JBTAPH010000036.1\tGenbank\tmRNA\t9\t966\t.\t-\t.\tID=rna-mrna.DD_M4_00007442-RA:cds;Parent=gene-ACTFIY_010592;gbkey=mRNA;locus_tag=ACTFIY_010592;orig_protein_id=gnl|WGS:JBTAPH|DD_M4_00007442-RA:cds;orig_transcript_id=gnl|WGS:JBTAPH|mrna.DD_M4_00007442-RA:cds;product=hypothetical protein",
                    "JBTAPH010000036.1\tGenbank\tCDS\t229\t966\t.\t-\t0\tID=cds-KAM9986187.1;Parent=rna-mrna.DD_M4_00007442-RA:cds;Dbxref=NCBI_GP:KAM9986187.1;Name=KAM9986187.1;gbkey=CDS;locus_tag=ACTFIY_010592;orig_transcript_id=gnl|WGS:JBTAPH|mrna.DD_M4_00007442-RA:cds;product=hypothetical protein;protein_id=KAM9986187.1",
                    "JBTAPH010000036.1\tGenbank\tCDS\t9\t158\t.\t-\t0\tID=cds-KAM9986187.1;Parent=rna-mrna.DD_M4_00007442-RA:cds;Dbxref=NCBI_GP:KAM9986187.1;Name=KAM9986187.1;gbkey=CDS;locus_tag=ACTFIY_010592;orig_transcript_id=gnl|WGS:JBTAPH|mrna.DD_M4_00007442-RA:cds;product=hypothetical protein;protein_id=KAM9986187.1",
                    "JBTAPH010000036.1\tGenbank\tgene\t2114\t2691\t.\t-\t.\tID=gene-ACTFIY_010593;Name=ACTFIY_010593;gbkey=Gene;gene_biotype=protein_coding;locus_tag=ACTFIY_010593",
                    "JBTAPH010000036.1\tGenbank\tmRNA\t2114\t2691\t.\t-\t.\tID=rna-mrna.DD_M4_00007443-RA:cds;Parent=gene-ACTFIY_010593;gbkey=mRNA;locus_tag=ACTFIY_010593;orig_protein_id=gnl|WGS:JBTAPH|DD_M4_00007443-RA:cds;orig_transcript_id=gnl|WGS:JBTAPH|mrna.DD_M4_00007443-RA:cds;product=hypothetical protein",
                    "JBTAPH010000036.1\tGenbank\tCDS\t2478\t2691\t.\t-\t0\tID=cds-KAM9986188.1;Parent=rna-mrna.DD_M4_00007443-RA:cds;Dbxref=NCBI_GP:KAM9986188.1;Name=KAM9986188.1;gbkey=CDS;locus_tag=ACTFIY_010593;orig_transcript_id=gnl|WGS:JBTAPH|mrna.DD_M4_00007443-RA:cds;product=hypothetical protein;protein_id=KAM9986188.1",
                    "JBTAPH010000036.1\tGenbank\tCDS\t2246\t2381\t.\t-\t0\tID=cds-KAM9986188.1;Parent=rna-mrna.DD_M4_00007443-RA:cds;Dbxref=NCBI_GP:KAM9986188.1;Name=KAM9986188.1;gbkey=CDS;locus_tag=ACTFIY_010593;orig_transcript_id=gnl|WGS:JBTAPH|mrna.DD_M4_00007443-RA:cds;product=hypothetical protein;protein_id=KAM9986188.1",
                    "JBTAPH010000036.1\tGenbank\tCDS\t2114\t2192\t.\t-\t0\tID=cds-KAM9986188.1;Parent=rna-mrna.DD_M4_00007443-RA:cds;Dbxref=NCBI_GP:KAM9986188.1;Name=KAM9986188.1;gbkey=CDS;locus_tag=ACTFIY_010593;orig_transcript_id=gnl|WGS:JBTAPH|mrna.DD_M4_00007443-RA:cds;product=hypothetical protein;protein_id=KAM9986188.1",
                    "",
                ]
            )
        )
    with gzip.open(genome_path, "wt", encoding="utf-8") as handle:
        handle.write(">JBTAPH010000036.1\n" + ("A" * 3000) + "\n")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "ncbi",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Dictyostelium_cf_discoideum_GCA_054859205.1_ASM5485920v1_cds_from_genomic.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [
        ">Dictyostelium_cf_discoideum_ACTFIY_010592",
        ">Dictyostelium_cf_discoideum_ACTFIY_010593",
    ]

    mapping = run_validate_mapping_script(
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
    )
    assert mapping.returncode == 0, mapping.stderr + "\n" + mapping.stdout
    assert "[Dictyostelium_cf_discoideum] CDS-to-GFF mapping OK: 2/2 IDs" in mapping.stdout


def test_format_species_inputs_uses_protein_id_when_ncbi_header_lacks_gene_tags(tmp_path):
    input_dir = tmp_path / "NCBI_Genome" / "species_wise_original"
    species_dir = input_dir / "Dictyostelium_firmibasis"
    species_dir.mkdir(parents=True, exist_ok=True)
    cds_path = species_dir / "GCA_036169595.1_ASM3616959v1_cds_from_genomic.fna.gz"
    gff_path = species_dir / "GCA_036169595.1_ASM3616959v1_genomic.gff.gz"
    genome_path = species_dir / "GCA_036169595.1_ASM3616959v1_genomic.fna.gz"
    with gzip.open(cds_path, "wt", encoding="utf-8") as handle:
        handle.write(
            (
                ">lcl|CM069765.1_cds_KAK5581746.1_1 [protein=hypothetical protein] [protein_id=KAK5581746.1] [location=join(5022..5093,5192..5424)] [gbkey=CDS]\n"
                "ATGCAAACAAATACATTTAGCAATGTACCTGGCTCACTTAATATTGAAGACCTATTAAATAAAATAGAAACTGTAGTATT\n"
                ">lcl|CM069765.1_cds_KAK5581747.1_2 [protein=hypothetical protein] [protein_id=KAK5581747.1] [location=8575..9000] [gbkey=CDS]\n"
                "ATGGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAA\n"
            )
        )
    with gzip.open(gff_path, "wt", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "##gff-version 3",
                    "CM069765.1\tGenbank\tgene\t5010\t5424\t.\t+\t.\tID=gene-RB653_003324;Name=RB653_003324;gbkey=Gene;gene_biotype=protein_coding;locus_tag=RB653_003324",
                    "CM069765.1\tGenbank\tmRNA\t5010\t5424\t.\t+\t.\tID=rna-DFI_00002-RA;Parent=gene-RB653_003324;gbkey=mRNA;locus_tag=RB653_003324;orig_protein_id=gnl|WGS:JAVFKY|DFI_00002-RA:cds;orig_transcript_id=gnl|WGS:JAVFKY|DFI_00002-RA;product=hypothetical protein",
                    "CM069765.1\tGenbank\tCDS\t5022\t5093\t.\t+\t0\tID=cds-KAK5581746.1;Parent=rna-DFI_00002-RA;Dbxref=NCBI_GP:KAK5581746.1;Name=KAK5581746.1;gbkey=CDS;locus_tag=RB653_003324;orig_transcript_id=gnl|WGS:JAVFKY|DFI_00002-RA;product=hypothetical protein;protein_id=KAK5581746.1",
                    "CM069765.1\tGenbank\tCDS\t5192\t5424\t.\t+\t0\tID=cds-KAK5581746.1;Parent=rna-DFI_00002-RA;Dbxref=NCBI_GP:KAK5581746.1;Name=KAK5581746.1;gbkey=CDS;locus_tag=RB653_003324;orig_transcript_id=gnl|WGS:JAVFKY|DFI_00002-RA;product=hypothetical protein;protein_id=KAK5581746.1",
                    "CM069765.1\tGenbank\tgene\t8575\t9000\t.\t+\t.\tID=gene-RB653_003325;Name=RB653_003325;gbkey=Gene;gene_biotype=protein_coding;locus_tag=RB653_003325",
                    "CM069765.1\tGenbank\tmRNA\t8575\t9000\t.\t+\t.\tID=rna-DFI_00003-RA;Parent=gene-RB653_003325;gbkey=mRNA;locus_tag=RB653_003325;orig_protein_id=gnl|WGS:JAVFKY|DFI_00003-RA:cds;orig_transcript_id=gnl|WGS:JAVFKY|DFI_00003-RA;product=hypothetical protein",
                    "CM069765.1\tGenbank\tCDS\t8575\t9000\t.\t+\t0\tID=cds-KAK5581747.1;Parent=rna-DFI_00003-RA;Dbxref=NCBI_GP:KAK5581747.1;Name=KAK5581747.1;gbkey=CDS;locus_tag=RB653_003325;orig_transcript_id=gnl|WGS:JAVFKY|DFI_00003-RA;product=hypothetical protein;protein_id=KAK5581747.1",
                    "",
                ]
            )
        )
    with gzip.open(genome_path, "wt", encoding="utf-8") as handle:
        handle.write(">CM069765.1\n" + ("A" * 12000) + "\n")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "ncbi",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Dictyostelium_firmibasis_GCA_036169595.1_ASM3616959v1_cds_from_genomic.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [
        ">Dictyostelium_firmibasis_KAK5581746.1",
        ">Dictyostelium_firmibasis_KAK5581747.1",
    ]

    mapping = run_validate_mapping_script(
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
    )
    assert mapping.returncode == 0, mapping.stderr + "\n" + mapping.stdout
    assert "[Dictyostelium_firmibasis] CDS-to-GFF mapping OK: 2/2 IDs" in mapping.stdout


def test_format_species_inputs_strict_gene_grouping_keeps_misassigned_gene_ids_separate(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Arabidopsis_thaliana"
    species_dir.mkdir(parents=True, exist_ok=True)
    gff_path = species_dir / "Arabidopsis_thaliana.annotation.gff3"
    genome_path = species_dir / "Arabidopsis_thaliana.genome.fa"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=badGeneA",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=locusX.t1;Parent=badGeneA",
                "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cds1;Parent=locusX.t1",
                "chr1\tsrc\tCDS\t13\t18\t.\t+\t0\tID=cds2;Parent=locusX.t1",
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=badGeneB",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=locusX.t2;Parent=badGeneB",
                "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=cds3;Parent=locusX.t2",
                "chr1\tsrc\tCDS\t13\t18\t.\t+\t0\tID=cds4;Parent=locusX.t2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(">chr1\nATGAAACCCGGGTTTCCC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--gene-grouping-mode",
        "strict",
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Arabidopsis_thaliana_annotation.derived.cds.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [
        ">Arabidopsis_thaliana_badGeneA",
        ">Arabidopsis_thaliana_badGeneB",
    ]


def test_format_species_inputs_fernbase_prefers_primary_annotation_files(tmp_path):
    input_dir = tmp_path / "FernBase" / "species_wise_original"
    species_dir = input_dir / "Azolla_filiculoides"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "Azolla_filiculoides.CDS.lowconfidence_v1.1.fasta").write_text(
        ">Azfi_g1.t1\nATG\n", encoding="utf-8"
    )
    (species_dir / "Azolla_filiculoides.CDS.highconfidence_v1.1.fasta").write_text(
        ">Azfi_g1.t1\nATGAA\n", encoding="utf-8"
    )
    (species_dir / "Azolla_filiculoides.transcript.highconfidence_v1.1.fasta").write_text(
        ">Azfi_g1.t1\nATGAAA\n", encoding="utf-8"
    )
    (species_dir / "Azolla_filiculoides.gene_models.lowconfidence_v1.1.gff").write_text(
        "chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=Azfi_g1\n",
        encoding="utf-8",
    )
    (species_dir / "Azolla_filiculoides.gene_models.highconfidence_v1.1.gff").write_text(
        "chr1\tsrc\tgene\t1\t5\t.\t+\t.\tID=Azfi_g1\n",
        encoding="utf-8",
    )
    (species_dir / "Azolla_filiculoides.genome_v1.2.fasta").write_text(">chr1\nATGCATGC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "fernbase",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Azolla_filiculoides_CDS.highconfidence_v1.1.fa.gz"
    formatted_gff = out_gff / "Azolla_filiculoides_gene_models.highconfidence_v1.1.gff.gz"
    formatted_genome = out_genome / "Azolla_filiculoides_genome_v1.2.fa.gz"
    assert formatted_cds.exists()
    assert formatted_gff.exists()
    assert formatted_genome.exists()

    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert ">Azolla_filiculoides_Azfi_g1" in cds_text
    assert "ATGAAN" in cds_text
    assert "lowconfidence" not in formatted_gff.name


def test_format_species_inputs_fernbase_accepts_markerless_genome_fasta(tmp_path):
    input_dir = tmp_path / "FernBase" / "species_wise_original"
    species_dir = input_dir / "Ceratopteris_richardii"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "Crichardii_676_v2.0_cds.fa").write_text(">Crich_g1.t1\nATGAA\n", encoding="utf-8")
    (species_dir / "Crichardii_676_v2.1.gene.gff3").write_text(
        "chr1\tsrc\tgene\t1\t5\t.\t+\t.\tID=Crich_g1\n",
        encoding="utf-8",
    )
    (species_dir / "Crichardii_676_v2.0.fa").write_text(">chr1\nATGCATGC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "fernbase",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_genome = out_genome / "Ceratopteris_richardii_Crichardii_676_v2.0.fa.gz"
    assert formatted_genome.exists()


def test_format_species_inputs_fernbase_uses_plain_gene_tag_for_aggregation(tmp_path):
    input_dir = tmp_path / "FernBase" / "species_wise_original"
    species_dir = input_dir / "Pteridium_aquilinum"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "pt_aq.cds.fa").write_text(
        (
            ">pteridium_mrna-253 gene=pteridium_gene-252 seq_id=scaff1 type=cds\n"
            "ATGAA\n"
            ">pteridium_mrna-303 gene=pteridium_gene-295 seq_id=scaff1 type=cds\n"
            "ATGAAA\n"
            ">pteridium_mrna-308 gene=pteridium_gene-295 seq_id=scaff1 type=cds\n"
            "ATGAAATAG\n"
        ),
        encoding="utf-8",
    )
    (species_dir / "pt_aq.gff3").write_text(
        (
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=pteridium_gene-252\n"
            "chr1\tsrc\tgene\t20\t40\t.\t+\t.\tID=pteridium_gene-295\n"
        ),
        encoding="utf-8",
    )
    (species_dir / "pt_aq_final_genome.fasta").write_text(">chr1\nATGCATGCATGC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    species_summary = tmp_path / "gg_input_generation_species.tsv"
    completed = run_script(
        "--provider",
        "fernbase",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Pteridium_aquilinum_pt_aq.cds.fa.gz"
    assert formatted_cds.exists()
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [
        ">Pteridium_aquilinum_pteridium_gene-252",
        ">Pteridium_aquilinum_pteridium_gene-295",
    ]

    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["species_prefix"] == "Pteridium_aquilinum"
    assert rows[0]["aggregated_cds_removed"] == "1"
    assert rows[0]["cds_sequences_before"] == "3"
    assert rows[0]["cds_sequences_after"] == "2"


def test_write_fasta_records_gzip_prefers_seqkit(monkeypatch, tmp_path):
    mod = load_module()
    output_path = tmp_path / "species.fa.gz"
    calls = {}

    class FakeSeqkitProcess:
        def __init__(self, command, **kwargs):
            calls["command"] = command
            calls["kwargs"] = kwargs
            self.command = command
            self._stdin_pipe = FakeTextPipe()
            self.stdin = self._stdin_pipe
            self.stderr = io.StringIO("")

        def wait(self):
            output_arg = self.command[self.command.index("-o") + 1]
            with gzip.open(output_arg, "wt", encoding="utf-8") as handle:
                handle.write(self._stdin_pipe.getvalue())
            return 0

        def kill(self):
            return None

    monkeypatch.setenv("GG_TASK_CPUS", "3")
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/{}".format(name))
    monkeypatch.setattr(subprocess, "Popen", FakeSeqkitProcess)

    mod.write_fasta_records_gzip(output_path, [("seq1", "ATG"), ("seq2", "ATGA")])

    with gzip.open(output_path, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert ">seq1" in text
    assert ">seq2" in text
    assert calls["command"][:4] == ["/usr/bin/seqkit", "seq", "--threads", "3"]
    assert calls["command"][-1] == "-"


def test_write_gff_gzip_prefers_pigz(monkeypatch, tmp_path):
    mod = load_module()
    input_path = tmp_path / "input.gff3"
    input_path.write_text("chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=evm.model.g1\n", encoding="utf-8")
    output_path = tmp_path / "output.gff.gz"
    calls = {}

    class FakePigzProcess:
        def __init__(self, command, **kwargs):
            calls["command"] = command
            calls["kwargs"] = kwargs
            self.stdout = kwargs["stdout"]
            self._stdin_pipe = FakeTextPipe()
            self.stdin = self._stdin_pipe
            self.stderr = io.StringIO("")

        def wait(self):
            with gzip.GzipFile(fileobj=self.stdout, mode="wb") as handle:
                handle.write(self._stdin_pipe.getvalue().encode("utf-8"))
            return 0

        def kill(self):
            return None

    monkeypatch.setenv("GG_TASK_CPUS", "4")
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/{}".format(name))
    monkeypatch.setattr(subprocess, "Popen", FakePigzProcess)

    line_count = mod.write_gff_gzip(input_path, output_path)

    with gzip.open(output_path, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert line_count == 1
    assert "evm.model." not in text
    assert calls["command"] == ["/usr/bin/pigz", "-p", "4", "-c"]


def test_resolve_provider_download_limits_keeps_fernbase_and_insectbase_default_caps_at_two(monkeypatch):
    mod = load_module()
    monkeypatch.delenv("GG_INPUT_MAX_CONCURRENT_DOWNLOADS_FERNBASE", raising=False)
    monkeypatch.delenv("GG_INPUT_MAX_CONCURRENT_DOWNLOADS_INSECTBASE", raising=False)
    limits = mod.resolve_provider_download_limits(8)
    assert limits["fernbase"] == 2
    assert limits["insectbase"] == 2


def test_resolve_ncbi_download_urls_from_id_retries_transient_remote_disconnect(monkeypatch):
    mod = load_module()
    calls = {"esearch": 0, "esummary": 0}

    def fake_urlopen(request, timeout):
        url = request.full_url
        if "esearch.fcgi" in url:
            calls["esearch"] += 1
            if calls["esearch"] == 1:
                raise RemoteDisconnected("Remote end closed connection without response")
            payload = {
                "header": {"type": "esearch", "version": "0.3"},
                "esearchresult": {"idlist": ["12345"]},
            }
            return FakeBinaryResponse(json.dumps(payload).encode("utf-8"))
        if "esummary.fcgi" in url:
            calls["esummary"] += 1
            payload = {
                "header": {"type": "esummary", "version": "0.3"},
                "result": {
                    "uids": ["12345"],
                    "12345": {
                        "ftppath_genbank": (
                            "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/036/169/595/GCA_036169595.1_ASM3616959v1"
                        ),
                        "organism": "Dictyostelium firmibasis",
                        "speciesname": "Dictyostelium firmibasis",
                    },
                },
            }
            return FakeBinaryResponse(json.dumps(payload).encode("utf-8"))
        raise AssertionError(url)

    provider_module = sys.modules[mod.resolve_ncbi_download_urls_from_id.__module__]
    monkeypatch.setattr(provider_module, "urlopen", fake_urlopen)
    monkeypatch.setattr(provider_module, "throttle_ncbi_eutils_request", lambda: None)
    monkeypatch.setattr(provider_module.time, "sleep", lambda _seconds: None)

    resolved = mod.resolve_ncbi_download_urls_from_id("GCA_036169595.1", timeout=1.0)
    assert calls["esearch"] == 2
    assert calls["esummary"] == 1
    assert resolved["species_key"] == "Dictyostelium_firmibasis"
    assert resolved["cds_filename"] == "GCA_036169595.1_ASM3616959v1_cds_from_genomic.fna.gz"
    assert resolved["ncbi_source_db"] == "genbank"


def test_iter_fasta_records_reads_tar_bz2_archive(tmp_path):
    mod = load_module()
    archive_path = tmp_path / "example.genome.fa.tar.bz2"
    payload = ">chr1 description\nATGC\nATGC\n".encode("utf-8")
    with tarfile.open(archive_path, "w:bz2") as archive:
        info = tarfile.TarInfo(name="nested/example.genome.fa")
        info.size = len(payload)
        archive.addfile(info, io.BytesIO(payload))
        note = tarfile.TarInfo(name="README.txt")
        note_payload = b"fixture\n"
        note.size = len(note_payload)
        archive.addfile(note, io.BytesIO(note_payload))

    records = list(mod.iter_fasta_records(archive_path))
    assert records == [("chr1 description", "ATGCATGC")]


def test_format_species_inputs_fernbase_prefers_namespaced_transcript_gene_id_when_gene_tag_is_short(tmp_path):
    input_dir = tmp_path / "FernBase" / "species_wise_original"
    species_dir = input_dir / "Salvinia_molesta"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "sg2.cds").write_text(
        (">sg2.g13251.t1 gene=g13251\nATGAAATAG\n>sg2.g9.t1 gene=g9\nATGAAATAA\n"),
        encoding="utf-8",
    )
    (species_dir / "sg2.gff3").write_text(
        ("Chr_1\tgmst\tgene\t1\t9\t.\t+\t.\tID=sg2.g13251\nChr_1\tgmst\tgene\t20\t28\t.\t+\t.\tID=sg2.g9\n"),
        encoding="utf-8",
    )
    (species_dir / "sg2_genome.fasta").write_text(">Chr_1\nATGCATGCATGC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "fernbase",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Salvinia_molesta_sg2.cds.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [
        ">Salvinia_molesta_sg2.g13251",
        ">Salvinia_molesta_sg2.g9",
    ]


def test_format_species_inputs_fernbase_strips_amt_suffix_when_gff_uses_base_gene_id(tmp_path):
    input_dir = tmp_path / "FernBase" / "species_wise_original"
    species_dir = input_dir / "Azolla_filiculoides"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "Azolla_filiculoides.CDS.highconfidence_v1.1.fasta").write_text(
        (">Azfi_s0034.g025227.AMT2\nATGAAATAG\n>Azfi_s0093.g043301.AMT2\nATGAAATAA\n"),
        encoding="utf-8",
    )
    (species_dir / "Azolla_filiculoides.gene_models.highconfidence_v1.1.gff").write_text(
        (
            "SCAF_1\tAUGUSTUS\tgene\t1\t9\t.\t+\t.\tID=Azfi_s0034.g025227\n"
            "SCAF_1\tAUGUSTUS\tgene\t20\t28\t.\t+\t.\tID=Azfi_s0093.g043301\n"
        ),
        encoding="utf-8",
    )
    (species_dir / "Azolla_filiculoides.genome_v1.2.fasta").write_text(">SCAF_1\nATGCATGCATGC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "fernbase",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Azolla_filiculoides_CDS.highconfidence_v1.1.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [
        ">Azolla_filiculoides_Azfi_s0034.g025227",
        ">Azolla_filiculoides_Azfi_s0093.g043301",
    ]


def test_species_summary_is_incremental_and_persistent_across_runs(tmp_path):
    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    species_summary = tmp_path / "gg_input_generation_species.tsv"

    first = run_script(
        "--provider",
        "ensemblplants",
        "--input-dir",
        str(SMALL_DATASET_ROOT / "20230216_EnsemblPlants" / "original_files"),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert first.returncode == 0, first.stderr + "\n" + first.stdout
    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows_first = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows_first) == 1
    assert rows_first[0]["provider"] == "ensemblplants"
    assert rows_first[0]["species_prefix"] == "Ostreococcus_lucimarinus"
    assert int(rows_first[0]["cds_sequences_before"]) >= int(rows_first[0]["cds_sequences_after"])
    assert rows_first[0]["cds_first_sequence_name"] != ""

    second = run_script(
        "--provider",
        "phycocosm",
        "--input-dir",
        str(SMALL_DATASET_ROOT / "PhycoCosm" / "species_wise_original"),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert second.returncode == 0, second.stderr + "\n" + second.stdout
    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows_second = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows_second) == 2
    by_key = {(row["provider"], row["species_prefix"]): row for row in rows_second}
    assert ("ensemblplants", "Ostreococcus_lucimarinus") in by_key
    assert ("phycocosm", "Microglena_spYARC") in by_key
    assert by_key[("phycocosm", "Microglena_spYARC")]["cds_first_sequence_name"] != ""


def test_format_species_inputs_derives_from_gbff_and_genome_when_gff_and_cds_are_missing(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Arabidopsis_thaliana"
    species_dir.mkdir(parents=True, exist_ok=True)
    gbff_path = species_dir / "Arabidopsis_thaliana.genomic.gbff"
    genome_path = species_dir / "Arabidopsis_thaliana.genome.fa"
    gbff_path.write_text(
        "\n".join(
            [
                "LOCUS       chr1               9 bp    DNA     linear   PLN 01-JAN-2000",
                "DEFINITION  test.",
                "ACCESSION   chr1",
                "VERSION     chr1",
                "FEATURES             Location/Qualifiers",
                "     gene            1..9",
                '                     /locus_tag="gene1"',
                '                     /gene="gene1"',
                "     CDS             join(1..3,7..9)",
                '                     /locus_tag="gene1"',
                '                     /gene="gene1"',
                '                     /protein_id="gene1.t1"',
                "ORIGIN",
                "        1 atgaaattt",
                "//",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(">chr1\nATGAAATTT\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    species_summary = tmp_path / "gg_input_generation_species.tsv"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Arabidopsis_thaliana_genomic.derived.cds.fa.gz"
    formatted_gff = out_gff / "Arabidopsis_thaliana_genomic.derived.gff.gz"
    formatted_genome = out_genome / "Arabidopsis_thaliana_genome.fa.gz"
    assert formatted_cds.exists()
    assert formatted_gff.exists()
    assert formatted_genome.exists()

    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert cds_text.count(">Arabidopsis_thaliana_gene1") == 1
    assert "ATGTTT" in cds_text

    with gzip.open(formatted_gff, "rt", encoding="utf-8") as handle:
        gff_text = handle.read()
    assert "##gff-version 3" in gff_text
    assert "\tgene\t" in gff_text
    assert "\tmRNA\t" in gff_text
    assert "\tCDS\t" in gff_text

    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert str(gbff_path) in row["cds_input_path"]
    assert "derived CDS" in row["cds_input_path"]
    assert str(gbff_path) in row["gff_input_path"]
    assert "derived GFF" in row["gff_input_path"]


def test_format_species_inputs_does_not_write_empty_gbff_derived_outputs(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Fakus_emptyus"
    species_dir.mkdir(parents=True, exist_ok=True)
    gbff_path = species_dir / "Fakus_emptyus.genomic.gbff"
    genome_path = species_dir / "Fakus_emptyus.genome.fa"
    gbff_path.write_text(
        "\n".join(
            [
                "LOCUS       chr1               9 bp    DNA     linear   PLN 01-JAN-2000",
                "DEFINITION  test.",
                "ACCESSION   chr1",
                "VERSION     chr1",
                "FEATURES             Location/Qualifiers",
                "     gene            1..9",
                '                     /locus_tag="gene1"',
                '                     /gene="gene1"',
                "ORIGIN",
                "        1 atgaaattt",
                "//",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(">chr1\nATGAAATTT\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    species_summary = tmp_path / "gg_input_generation_species.tsv"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert completed.returncode == 2
    assert "derived CDS contained no records" in completed.stderr
    assert "derived GFF contained no feature rows" in completed.stderr
    assert "CDS=Fakus_emptyus.genome.fa (derived CDS) (empty, NA" in completed.stdout
    assert "GFF=Fakus_emptyus.genomic.gbff (derived GFF) (empty, lines=1)" in completed.stdout
    assert not (out_cds / "Fakus_emptyus_genomic.derived.cds.fa.gz").exists()
    assert not (out_gff / "Fakus_emptyus_genomic.derived.gff.gz").exists()
    assert (out_genome / "Fakus_emptyus_genome.fa.gz").exists()
    assert species_summary.exists()
    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        assert len(list(csv.DictReader(handle, delimiter="\t"))) == 0


def test_ensembl_like_gff_selection_prefers_full_annotation_over_partial_files():
    module = load_module()
    candidates = [
        "https://ftp.example/Saccharomyces_cerevisiae.R64-1-1.115.abinitio.gff3.gz",
        "https://ftp.example/Saccharomyces_cerevisiae.R64-1-1.115.chromosome.I.gff3.gz",
        "https://ftp.example/Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz",
    ]
    assert module.select_best_url_for_label("ensembl", "GFF", candidates).endswith(".115.gff3.gz")

    candidates = [
        "https://ftp.example/Penaeus_japonicus_gca017312705v1.Mj_TUMSAT_v1.0.62.chr.gff3.gz",
        "https://ftp.example/Penaeus_japonicus_gca017312705v1.Mj_TUMSAT_v1.0.62.primary_assembly.NC_007010.1.gff3.gz",
        "https://ftp.example/Penaeus_japonicus_gca017312705v1.Mj_TUMSAT_v1.0.62.gff3.gz",
    ]
    assert module.select_best_url_for_label("ensemblmetazoa", "GFF", candidates).endswith(".62.gff3.gz")


def test_ensembl_like_explicit_partial_gff_url_is_replaced_with_full_annotation(monkeypatch):
    module = load_module()
    provider_module = sys.modules[module.resolve_preferred_ensembl_like_gff_url.__module__]

    partial_url = "https://ftp.example/current_gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.abinitio.gff3.gz"
    full_url = "https://ftp.example/current_gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz"

    def fake_resolve_urls_from_index_url(provider, index_url, timeout, headers):
        assert provider == "ensembl"
        assert index_url == "https://ftp.example/current_gff3/saccharomyces_cerevisiae/"
        return {"gff_url": full_url}

    monkeypatch.setattr(provider_module, "resolve_urls_from_index_url", fake_resolve_urls_from_index_url)

    assert module.resolve_preferred_ensembl_like_gff_url("ensembl", partial_url, 1, {}) == full_url


def test_remove_stale_ensembl_like_partial_gff_outputs_keeps_full_annotation(tmp_path):
    module = load_module()
    gff_dir = tmp_path / "species_gff"
    gff_dir.mkdir()
    full = gff_dir / "Saccharomyces_cerevisiae_R64-1-1.115.gff.gz"
    abinitio = gff_dir / "Saccharomyces_cerevisiae_R64-1-1.115.abinitio.gff.gz"
    chromosome = gff_dir / "Saccharomyces_cerevisiae_R64-1-1.115.chromosome.I.gff.gz"
    unrelated = gff_dir / "Drosophila_melanogaster_BDGP6.54.115.abinitio.gff.gz"
    longer_species = gff_dir / "Saccharomyces_cerevisiae_subsp_x_R64-1-1.115.abinitio.gff.gz"
    for path in (full, abinitio, chromosome, unrelated, longer_species):
        path.write_bytes(b"dummy")
    abinitio_audit = Path(str(abinitio) + ".repair.json")
    chromosome_audit = Path(str(chromosome) + ".repair.json")
    abinitio_audit.write_text("{}", encoding="utf-8")
    chromosome_audit.write_text("{}", encoding="utf-8")

    removed = module.remove_stale_ensembl_like_partial_gff_outputs(
        "ensembl",
        "Saccharomyces_cerevisiae",
        gff_dir,
        full,
    )

    assert sorted(removed) == sorted([abinitio.name, chromosome.name])
    assert full.exists()
    assert unrelated.exists()
    assert longer_species.exists()
    assert not abinitio.exists()
    assert not chromosome.exists()
    assert not abinitio_audit.exists()
    assert not chromosome_audit.exists()
