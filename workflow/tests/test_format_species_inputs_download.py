from pathlib import Path
import csv
import subprocess
import sys


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "format_species_inputs.py"


def run_script(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def make_manifest(path, rows):
    with open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["provider", "species_key", "cds_url", "gff_url", "cds_filename", "gff_filename"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def to_file_url(path):
    return path.resolve().as_uri()


def test_download_manifest_then_format_ensembl(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()

    cds_source = source_dir / "ensembl_cds.fa"
    gff_source = source_dir / "ensembl_gene.gff3"
    cds_source.write_text(
        (
            ">x cds chromosome:chr1:1:5:1 gene:OSTLU_25062\n"
            "ATGAA\n"
        ),
        encoding="utf-8",
    )
    gff_source.write_text(
        "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=evm.model.OSTLU_25062\n",
        encoding="utf-8",
    )

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "ensemblplants",
                "species_key": "Ostreococcus_lucimarinus",
                "cds_url": to_file_url(cds_source),
                "gff_url": to_file_url(gff_source),
                "cds_filename": "Ostreococcus_lucimarinus.ASM9206v1.cds.all.fa",
                "gff_filename": "Ostreococcus_lucimarinus.ASM9206v1.56.gff3",
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    out_cds = tmp_path / "out_cds"
    out_gff = tmp_path / "out_gff"
    completed = run_script(
        "--provider",
        "ensemblplants",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(download_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    raw_cds = download_dir / "20230216_EnsemblPlants" / "original_files" / "Ostreococcus_lucimarinus.ASM9206v1.cds.all.fa"
    raw_gff = download_dir / "20230216_EnsemblPlants" / "original_files" / "Ostreococcus_lucimarinus.ASM9206v1.56.gff3"
    assert raw_cds.exists()
    assert raw_gff.exists()

    formatted_cds = out_cds / "Ostreococcus_lucimarinus_ASM9206v1.cds.all.fa"
    formatted_gff = out_gff / "Ostreococcus_lucimarinus_ASM9206v1.56.gff3"
    assert formatted_cds.exists()
    assert formatted_gff.exists()

    cds_text = formatted_cds.read_text(encoding="utf-8")
    assert ">Ostreococcus_lucimarinus_OSTLU_25062" in cds_text
    assert "ATGAAN" in cds_text
    assert "evm.model." not in formatted_gff.read_text(encoding="utf-8")


def test_download_only_writes_raw_layout(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    cds_source = source_dir / "phyco_cds.fasta"
    gff_source = source_dir / "phyco_gene.gff3"
    cds_source.write_text(">jgi|X|1|mRNA.A\nATG\n", encoding="utf-8")
    gff_source.write_text("scaf\tsrc\tgene\t1\t3\t.\t+\t.\tID=A\n", encoding="utf-8")

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "phycocosm",
                "species_key": "Microglena_spYARC_MicrYARC1",
                "cds_url": to_file_url(cds_source),
                "gff_url": to_file_url(gff_source),
                "cds_filename": "MicrYARC1_GeneCatalog_CDS_20220803.fasta",
                "gff_filename": "MicrYARC1_GeneCatalog_genes_20220803.gff3",
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    completed = run_script(
        "--provider",
        "phycocosm",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(download_dir),
        "--download-only",
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    raw_dir = download_dir / "PhycoCosm" / "species_wise_original" / "Microglena_spYARC_MicrYARC1"
    assert (raw_dir / "MicrYARC1_GeneCatalog_CDS_20220803.fasta").exists()
    assert (raw_dir / "MicrYARC1_GeneCatalog_genes_20220803.gff3").exists()


def test_download_manifest_missing_required_columns_fails(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    cds_source = source_dir / "cds.fa"
    cds_source.write_text(">a\nATG\n", encoding="utf-8")

    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "provider\tspecies_key\tcds_url\n"
        "phycocosm\tMicroglena_spYARC_MicrYARC1\t{}\n".format(to_file_url(cds_source)),
        encoding="utf-8",
    )

    completed = run_script(
        "--provider",
        "phycocosm",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(tmp_path / "download_cache"),
        "--download-only",
    )
    assert completed.returncode == 1, completed.stderr + "\n" + completed.stdout
    assert "missing required columns" in completed.stderr
