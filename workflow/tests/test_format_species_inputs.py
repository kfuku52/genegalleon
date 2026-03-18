from importlib.util import module_from_spec, spec_from_file_location
import io
from pathlib import Path
import csv
import gzip
import json
import os
import sqlite3
import shutil
import subprocess
import sys
import tarfile


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "format_species_inputs.py"
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


def write_test_taxonomy_fixture(tmp_path):
    db_path = tmp_path / "taxa.sqlite"
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("CREATE TABLE species (taxid INT PRIMARY KEY, parent INT, spname VARCHAR(50) COLLATE NOCASE, common VARCHAR(50) COLLATE NOCASE, rank VARCHAR(50), track TEXT)")
    cur.execute("CREATE TABLE synonym (taxid INT,spname VARCHAR(50) COLLATE NOCASE, PRIMARY KEY (spname, taxid))")
    species_rows = [
        (1, 1, "root", "", "no rank", "1"),
        (131567, 1, "cellular organisms", "", "cellular root", "131567,1"),
        (2759, 131567, "Eukaryota", "", "domain", "2759,131567,1"),
        (33090, 2759, "Viridiplantae", "", "kingdom", "33090,2759,131567,1"),
        (242159, 33090, "Ostreococcus lucimarinus", "", "species", "242159,33090,2759,131567,1"),
        (5911, 2759, "Tetrahymena thermophila", "", "species", "5911,2759,131567,1"),
    ]
    cur.executemany("INSERT INTO species (taxid, parent, spname, common, rank, track) VALUES (?, ?, ?, ?, ?, ?)", species_rows)
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


def test_format_species_inputs_strict_mode_fails_on_missing_pair(tmp_path):
    dataset_copy = tmp_path / "small_gfe_dataset_copy"
    shutil.copytree(SMALL_DATASET_ROOT, dataset_copy)
    missing_gff = dataset_copy / "20230216_EnsemblPlants" / "original_files" / "Ostreococcus_lucimarinus.ASM9206v1.56.gff3"
    missing_gff.unlink()

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    completed = run_script(
        "--provider",
        "ensemblplants",
        "--input-dir",
        str(dataset_copy / "20230216_EnsemblPlants" / "original_files"),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--strict",
    )
    assert completed.returncode == 1, completed.stderr + "\n" + completed.stdout
    assert "missing GFF" in completed.stderr


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
    raw_header = (
        "AT1G01520.5 cds chromosome:TAIR10:1:190478:192436:1 "
        "gene:AT1G01520 gene_symbol:ASG4"
    )
    derived_header = "transcript:AT1G01520.1 gene=AT1G01520"

    assert module.build_formatted_cds_id(task, raw_header) == "Arabidopsis_thaliana_AT1G01520.5"
    assert module.build_formatted_cds_id(task, derived_header) == "Arabidopsis_thaliana_AT1G01520.1"
    assert module.build_gene_aggregate_id(task, raw_header, "Arabidopsis_thaliana_AT1G01520.5") == "Arabidopsis_thaliana_AT1G01520"
    assert module.build_gene_aggregate_id(task, derived_header, "Arabidopsis_thaliana_AT1G01520.1") == "Arabidopsis_thaliana_AT1G01520"


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


def test_format_species_inputs_fernbase_prefers_primary_annotation_files(tmp_path):
    input_dir = tmp_path / "FernBase" / "species_wise_original"
    species_dir = input_dir / "Azolla_filiculoides"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "Azolla_filiculoides.CDS.lowconfidence_v1.1.fasta").write_text(">Azfi_g1.t1\nATG\n", encoding="utf-8")
    (species_dir / "Azolla_filiculoides.CDS.highconfidence_v1.1.fasta").write_text(">Azfi_g1.t1\nATGAA\n", encoding="utf-8")
    (species_dir / "Azolla_filiculoides.transcript.highconfidence_v1.1.fasta").write_text(">Azfi_g1.t1\nATGAAA\n", encoding="utf-8")
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
    monkeypatch.setattr(mod.shutil, "which", lambda name: "/usr/bin/{}".format(name))
    monkeypatch.setattr(mod.subprocess, "Popen", FakeSeqkitProcess)

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
    monkeypatch.setattr(mod.shutil, "which", lambda name: "/usr/bin/{}".format(name))
    monkeypatch.setattr(mod.subprocess, "Popen", FakePigzProcess)

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
        (
            ">sg2.g13251.t1 gene=g13251\n"
            "ATGAAATAG\n"
            ">sg2.g9.t1 gene=g9\n"
            "ATGAAATAA\n"
        ),
        encoding="utf-8",
    )
    (species_dir / "sg2.gff3").write_text(
        (
            "Chr_1\tgmst\tgene\t1\t9\t.\t+\t.\tID=sg2.g13251\n"
            "Chr_1\tgmst\tgene\t20\t28\t.\t+\t.\tID=sg2.g9\n"
        ),
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
        (
            ">Azfi_s0034.g025227.AMT2\n"
            "ATGAAATAG\n"
            ">Azfi_s0093.g043301.AMT2\n"
            "ATGAAATAA\n"
        ),
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
