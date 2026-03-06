from pathlib import Path
import csv
import gzip
import json
import shutil
import subprocess
import sys


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "format_species_inputs.py"
SMALL_DATASET_ROOT = Path(__file__).resolve().parent / "data" / "small_gfe_dataset"


def run_script(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), *args],
        capture_output=True,
        text=True,
        check=False,
    )


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

    ensembl_gff = out_gff / "Ostreococcus_lucimarinus_ASM9206v1.56.gff3"
    phytozome_gff = out_gff / "Hydrocotyle_leucocephala_HleucocephalaHAP1_768_v2.1.gene.gff3"
    phytozome_exons = out_gff / "Hydrocotyle_leucocephala_HleucocephalaHAP1_768_v2.1.gene_exons.gff3"
    assert ensembl_gff.exists()
    assert phytozome_gff.exists()
    assert not phytozome_exons.exists()

    gff_text = ensembl_gff.read_text(encoding="utf-8")
    assert "evm.model." not in gff_text
    assert "Oropetium_20150105_" not in gff_text

    phytozome_gff_text = phytozome_gff.read_text(encoding="utf-8")
    assert "evm_27.model." not in phytozome_gff_text
    stats = json.loads(stats_json.read_text(encoding="utf-8"))
    assert stats["species_processed"] == 3
    assert stats["num_species_cds_files"] == 3
    assert stats["num_species_gff_files"] == 3
    assert stats["cds_sequences_before"] >= stats["cds_sequences_after"]
    assert stats["cds_first_sequence_name"] != ""


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
    formatted_gff = out_gff / "Azolla_filiculoides_gene_models.highconfidence_v1.1.gff"
    formatted_genome = out_genome / "Azolla_filiculoides_genome_v1.2.fa.gz"
    assert formatted_cds.exists()
    assert formatted_gff.exists()
    assert formatted_genome.exists()

    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert ">Azolla_filiculoides_Azfi_g1" in cds_text
    assert "ATGAAN" in cds_text
    assert "lowconfidence" not in formatted_gff.name


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
