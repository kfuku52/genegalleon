from pathlib import Path
import csv
import gzip
import io
import json
import os
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import subprocess
import sys
import threading
from urllib.parse import parse_qs, urlparse
import zipfile

from openpyxl import Workbook


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "format_species_inputs.py"
SMALL_DATASET_ROOT = Path(__file__).resolve().parent / "data" / "small_gfe_dataset"


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
            fieldnames=[
                "provider",
                "id",
                "species_key",
                "cds_url",
                "gff_url",
                "genome_url",
                "cds_filename",
                "gff_filename",
                "genome_filename",
                "cds_url_template",
                "gff_url_template",
                "genome_url_template",
                "local_cds_path",
                "local_gff_path",
                "local_genome_path",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            payload = dict(row)
            if str(payload.get("id", "") or "").strip() == "":
                payload["id"] = str(payload.get("species_key", "") or "").strip()
            writer.writerow(payload)


def make_manifest_xlsx(path, headers, rows):
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "download_plan"
    sheet.append(list(headers))
    for row in rows:
        sheet.append([row.get(key, "") for key in headers])
    workbook.save(path)
    workbook.close()


def to_file_url(path):
    return path.resolve().as_uri()


class _NcbiFixtureHandler(SimpleHTTPRequestHandler):
    def __init__(self, *args, root_dir=None, **kwargs):
        self._root_dir = root_dir
        super().__init__(*args, directory=str(root_dir), **kwargs)

    def _send_json(self, payload):
        body = json.dumps(payload).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path.startswith("/eutils/esearch.fcgi"):
            self._send_json(
                {
                    "header": {"type": "esearch", "version": "0.3"},
                    "esearchresult": {"idlist": ["12345"]},
                }
            )
            return
        if self.path.startswith("/eutils/esummary.fcgi"):
            self._send_json(
                {
                    "header": {"type": "esummary", "version": "0.3"},
                    "result": {
                        "uids": ["12345"],
                        "12345": {
                            "ftppath_refseq": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14",
                            "organism": "Homo sapiens (human)",
                            "speciesname": "Homo sapiens",
                        },
                    },
                }
            )
            return
        super().do_GET()


def test_download_manifest_then_format_ensembl(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()

    cds_source = source_dir / "ensembl_cds.fa"
    gff_source = source_dir / "ensembl_gene.gff3"
    genome_source = source_dir / "ensembl_genome.fa"
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
    genome_source.write_text(
        ">chr1 chromosome\nATGCGT\n",
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
                "genome_url": to_file_url(genome_source),
                "cds_filename": "Ostreococcus_lucimarinus.ASM9206v1.cds.all.fa",
                "gff_filename": "Ostreococcus_lucimarinus.ASM9206v1.56.gff3",
                "genome_filename": "Ostreococcus_lucimarinus.ASM9206v1.dna.toplevel.fa",
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    out_cds = tmp_path / "out_cds"
    out_gff = tmp_path / "out_gff"
    out_genome = tmp_path / "out_genome"
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
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    raw_cds = download_dir / "20230216_EnsemblPlants" / "original_files" / "Ostreococcus_lucimarinus.ASM9206v1.cds.all.fa"
    raw_gff = download_dir / "20230216_EnsemblPlants" / "original_files" / "Ostreococcus_lucimarinus.ASM9206v1.56.gff3"
    raw_genome = download_dir / "20230216_EnsemblPlants" / "original_files" / "Ostreococcus_lucimarinus.ASM9206v1.dna.toplevel.fa"
    assert raw_cds.exists()
    assert raw_gff.exists()
    assert raw_genome.exists()

    formatted_cds = out_cds / "Ostreococcus_lucimarinus_ASM9206v1.cds.all.fa.gz"
    formatted_gff = out_gff / "Ostreococcus_lucimarinus_ASM9206v1.56.gff3"
    formatted_genome = out_genome / "Ostreococcus_lucimarinus_ASM9206v1.dna.toplevel.fa.gz"
    assert formatted_cds.exists()
    assert formatted_gff.exists()
    assert formatted_genome.exists()

    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert ">Ostreococcus_lucimarinus_OSTLU_25062" in cds_text
    assert "ATGAAN" in cds_text
    assert "evm.model." not in formatted_gff.read_text(encoding="utf-8")
    with gzip.open(formatted_genome, "rt", encoding="utf-8") as handle:
        genome_text = handle.read()
    assert ">chr1" in genome_text


def test_download_manifest_xlsx_reads_provider_and_id_columns(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    cds_source = source_dir / "ensembl_cds.fa"
    gff_source = source_dir / "ensembl_gene.gff3"
    cds_source.write_text(">x cds gene:OSTLU_25062\nATGAA\n", encoding="utf-8")
    gff_source.write_text("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=OSTLU_25062\n", encoding="utf-8")

    manifest = tmp_path / "manifest.xlsx"
    make_manifest_xlsx(
        manifest,
        headers=("provider", "id", "species_key", "cds_url", "gff_url"),
        rows=[
            {
                "provider": "ensemblplants",
                "id": "Ostreococcus_lucimarinus",
                "species_key": "Ostreococcus_lucimarinus",
                "cds_url": to_file_url(cds_source),
                "gff_url": to_file_url(gff_source),
            }
        ],
    )

    completed = run_script(
        "--provider",
        "ensemblplants",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(tmp_path / "download_cache"),
        "--download-only",
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout


def test_download_manifest_xlsx_requires_provider_id_first_columns(tmp_path):
    manifest = tmp_path / "manifest_bad.xlsx"
    make_manifest_xlsx(
        manifest,
        headers=("id", "provider", "species_key"),
        rows=[{"id": "foo", "provider": "ensembl", "species_key": "foo"}],
    )
    completed = run_script(
        "--provider",
        "ensembl",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(tmp_path / "download_cache"),
        "--download-only",
    )
    assert completed.returncode == 1
    assert "first two columns" in completed.stderr


def test_download_manifest_provider_local_uses_phytozome_local_dir(tmp_path):
    local_species_dir = (
        SMALL_DATASET_ROOT
        / "Phytozome"
        / "species_wise_original"
        / "Hydrocotyle_leucocephala_HAP1v2.1"
    )
    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "local",
                "id": str(local_species_dir),
                "species_key": "",
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    out_cds = tmp_path / "out_cds"
    out_gff = tmp_path / "out_gff"
    out_genome = tmp_path / "out_genome"
    completed = run_script(
        "--provider",
        "local",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(download_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    raw_dir = download_dir / "Local" / "species_wise_original" / "Hydrocotyle_leucocephala_HAP1v2.1"
    assert (raw_dir / "HleucocephalaHAP1_768_v2.1.cds_primaryTranscriptOnly.fa").exists()
    assert (raw_dir / "HleucocephalaHAP1_768_v2.1.gene.gff3").exists()

    formatted_cds = out_cds / "Hydrocotyle_leucocephala_HleucocephalaHAP1_768_v2.1.cds_primaryTranscriptOnly.fa.gz"
    formatted_gff = out_gff / "Hydrocotyle_leucocephala_HleucocephalaHAP1_768_v2.1.gene.gff3"
    assert formatted_cds.exists()
    assert formatted_gff.exists()
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert ">Hydrocotyle_leucocephala_HyleuH1.06G006800" in text


def test_download_manifest_provider_local_accepts_labeled_id_choice(tmp_path):
    local_species_dir = (
        SMALL_DATASET_ROOT
        / "Phytozome"
        / "species_wise_original"
        / "Hydrocotyle_leucocephala_HAP1v2.1"
    )
    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "local",
                "id": "{} (Hydrocotyle leucocephala)".format(local_species_dir),
                "species_key": "",
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    completed = run_script(
        "--provider",
        "local",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(download_dir),
        "--download-only",
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert "provider=local id looked like a labeled choice" in completed.stderr

    raw_dir = download_dir / "Local" / "species_wise_original" / "Hydrocotyle_leucocephala_HAP1v2.1"
    assert (raw_dir / "HleucocephalaHAP1_768_v2.1.cds_primaryTranscriptOnly.fa").exists()
    assert (raw_dir / "HleucocephalaHAP1_768_v2.1.gene.gff3").exists()


def test_download_manifest_ncbi_accepts_explicit_cds_gff_without_accession_id(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    cds_source = source_dir / "custom_cds.fa"
    gff_source = source_dir / "custom_gene.gff3"
    cds_source.write_text(">tx1\nATGAA\n", encoding="utf-8")
    gff_source.write_text("chr1\tsrc\tgene\t1\t5\t.\t+\t.\tID=gene1\n", encoding="utf-8")

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "ncbi",
                "id": "not_an_accession",
                "species_key": "Homo_sapiens",
                "cds_url": to_file_url(cds_source),
                "gff_url": to_file_url(gff_source),
                "genome_url": "",
                "cds_filename": "",
                "gff_filename": "",
                "genome_filename": "",
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    completed = run_script(
        "--provider",
        "ncbi",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(download_dir),
        "--download-only",
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert "failed to resolve id 'not_an_accession'" not in completed.stderr

    raw_dir = download_dir / "NCBI_Genome" / "species_wise_original" / "Homo_sapiens"
    assert (raw_dir / "custom_cds.fa").exists()
    assert (raw_dir / "custom_gene.gff3").exists()


def test_download_manifest_rejects_phycocosm_provider(tmp_path):
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
    assert completed.returncode != 0
    assert "does not support provider 'phycocosm'" in completed.stderr


def test_download_manifest_rejects_phytozome_provider(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    cds_source = source_dir / "phyto_cds.fa"
    gff_source = source_dir / "phyto_gene.gff3"
    cds_source.write_text(">tx1\nATG\n", encoding="utf-8")
    gff_source.write_text("scaf\tsrc\tgene\t1\t3\t.\t+\t.\tID=g1\n", encoding="utf-8")

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "phytozome",
                "species_key": "Hydrocotyle_leucocephala_HleucocephalaHAP1_768_v2.1",
                "cds_url": to_file_url(cds_source),
                "gff_url": to_file_url(gff_source),
                "cds_filename": "Hydrocotyle.cds.fa",
                "gff_filename": "Hydrocotyle.gff3",
            }
        ],
    )

    completed = run_script(
        "--provider",
        "phytozome",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(tmp_path / "download_cache"),
        "--download-only",
    )
    assert completed.returncode != 0
    assert "does not support provider 'phytozome'" in completed.stderr


def test_download_manifest_writes_resolved_manifest_tsv(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    (source_dir / "homo_sapiens.cds.fa").write_text(">t1\nATG\n", encoding="utf-8")
    (source_dir / "homo_sapiens.gff3").write_text(
        "scaf\tsrc\tgene\t1\t3\t.\t+\t.\tID=g1\n",
        encoding="utf-8",
    )

    manifest = tmp_path / "manifest.tsv"
    source_uri_root = source_dir.resolve().as_uri()
    make_manifest(
        manifest,
        [
            {
                "provider": "ensembl",
                "id": "homo_sapiens",
                "species_key": "Homo_sapiens",
                "cds_url_template": source_uri_root + "/{id}.cds.fa",
                "gff_url_template": source_uri_root + "/{id}.gff3",
                "cds_filename": "Homo_sapiens.cds.fa",
                "gff_filename": "Homo_sapiens.gff3",
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    resolved_manifest = tmp_path / "resolved_download_plan.tsv"
    completed = run_script(
        "--provider",
        "ensembl",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(download_dir),
        "--resolved-manifest-output",
        str(resolved_manifest),
        "--download-only",
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert resolved_manifest.exists()

    with open(resolved_manifest, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert row["provider"] == "ensembl"
    assert row["id"] == "homo_sapiens"
    assert row["species_key"] == "Homo_sapiens"
    assert row["cds_url"] == source_uri_root + "/homo_sapiens.cds.fa"
    assert row["gff_url"] == source_uri_root + "/homo_sapiens.gff3"
    assert row["cds_filename"] == "Homo_sapiens.cds.fa"
    assert row["gff_filename"] == "Homo_sapiens.gff3"


def test_download_manifest_missing_required_columns_fails(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    cds_source = source_dir / "cds.fa"
    cds_source.write_text(">a\nATG\n", encoding="utf-8")

    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "species_key\tcds_url\n"
        "Microglena_spYARC_MicrYARC1\t{}\n".format(to_file_url(cds_source)),
        encoding="utf-8",
    )

    completed = run_script(
        "--provider",
        "all",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(tmp_path / "download_cache"),
        "--download-only",
    )
    assert completed.returncode == 1, completed.stderr + "\n" + completed.stdout
    assert "required columns provider,id" in completed.stderr


def test_download_manifest_requires_provider_value(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    cds_source = source_dir / "ncbi_cds.fna.gz"
    gff_source = source_dir / "ncbi_genomic.gff.gz"
    genome_source = source_dir / "ncbi_genomic.fna.gz"
    with gzip.open(cds_source, "wt", encoding="utf-8") as handle:
        handle.write(">lcl|NC_000001.11_cds_NP_000001.1_1 [gene=ABC1] [db_xref=GeneID:111]\nATGAA\n")
    with gzip.open(gff_source, "wt", encoding="utf-8") as handle:
        handle.write("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1\n")
    with gzip.open(genome_source, "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGCATGC\n")

    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        (
            "provider\tid\tspecies_key\tcds_url\tgff_url\tgenome_url\tcds_filename\tgff_filename\tgenome_filename\n"
            "\tGCF_000001405.40\tHomo_sapiens\t{cds}\t{gff}\t{genome}\tGCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz\tGCF_000001405.40_GRCh38.p14_genomic.gff.gz\tGCF_000001405.40_GRCh38.p14_genomic.fna.gz\n"
        ).format(
            cds=to_file_url(cds_source),
            gff=to_file_url(gff_source),
            genome=to_file_url(genome_source),
        ),
        encoding="utf-8",
    )

    download_dir = tmp_path / "download_cache"
    completed = run_script(
        "--provider",
        "all",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(download_dir),
        "--download-only",
    )
    assert completed.returncode == 1, completed.stderr + "\n" + completed.stdout
    assert "provider is required" in completed.stderr


def test_download_manifest_requires_id_value(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    cds_source = source_dir / "ncbi_cds.fna.gz"
    gff_source = source_dir / "ncbi_genomic.gff.gz"
    with gzip.open(cds_source, "wt", encoding="utf-8") as handle:
        handle.write(">lcl|NC_000001.11_cds_NP_000001.1_1 [gene=ABC1] [db_xref=GeneID:111]\nATGAA\n")
    with gzip.open(gff_source, "wt", encoding="utf-8") as handle:
        handle.write("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1\n")

    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        (
            "provider\tid\tspecies_key\tcds_url\tgff_url\n"
            "ncbi\t\tHomo_sapiens\t{cds}\t{gff}\n"
        ).format(cds=to_file_url(cds_source), gff=to_file_url(gff_source)),
        encoding="utf-8",
    )
    completed = run_script(
        "--provider",
        "ncbi",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(tmp_path / "download_cache"),
        "--download-only",
    )
    assert completed.returncode == 1, completed.stderr + "\n" + completed.stdout
    assert "id is required" in completed.stderr


def test_download_manifest_ensembl_provider_from_id_prefix(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    species_id = "homo_sapiens"
    cds_source = source_dir / (species_id + ".cds.fa")
    gff_source = source_dir / (species_id + ".gene.gff3")
    genome_source = source_dir / (species_id + ".genome.fa")
    cds_source.write_text(">ENSG000001_t1\nATGAA\n", encoding="utf-8")
    gff_source.write_text("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=ENSG000001\n", encoding="utf-8")
    genome_source.write_text(">chr1\nATGCATGC\n", encoding="utf-8")

    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        (
            "provider\tid\tspecies_key\tcds_url\tgff_url\tgenome_url\n"
            "ensembl\tensembl:homo_sapiens (Homo sapiens)\tHomo_sapiens\t\t\t\n"
        ),
        encoding="utf-8",
    )

    env = dict(os.environ)
    env["GG_ENSEMBL_CDS_URL_TEMPLATE"] = source_dir.resolve().as_uri() + "/{id}.cds.fa"
    env["GG_ENSEMBL_GFF_URL_TEMPLATE"] = source_dir.resolve().as_uri() + "/{id}.gene.gff3"
    env["GG_ENSEMBL_GENOME_URL_TEMPLATE"] = source_dir.resolve().as_uri() + "/{id}.genome.fa"

    download_dir = tmp_path / "download_cache"
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--provider",
            "all",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(download_dir),
            "--download-only",
        ],
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    raw_dir = download_dir / "Ensembl" / "original_files"
    assert (raw_dir / "Homo_sapiens.homo_sapiens.cds.fa").exists()
    assert (raw_dir / "Homo_sapiens.homo_sapiens.gene.gff3").exists()
    assert (raw_dir / "Homo_sapiens.homo_sapiens.genome.fa").exists()


def test_download_manifest_supports_coge_and_cngb_with_id_inference(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    species_key = "Arabidopsis_thaliana"

    provider_layout = {
        "coge": "CoGe",
        "cngb": "CNGB",
    }
    provider_ids = {
        "coge": "24739",
        "cngb": "cngb:arabidopsis_thaliana_v1",
    }

    for provider in ("coge", "cngb"):
        cds_source = source_dir / "{}_cds.fa".format(provider)
        gff_source = source_dir / "{}_gene.gff3".format(provider)
        genome_source = source_dir / "{}_genome.fa".format(provider)
        cds_source.write_text(
            (
                ">AT1G01010_t1\n"
                "ATGAA\n"
                ">AT1G01010_t2\n"
                "ATGAAATTT\n"
            ),
            encoding="utf-8",
        )
        gff_source.write_text(
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=AT1G01010\n",
            encoding="utf-8",
        )
        genome_source.write_text(
            ">chr1\nATGCATGC\n",
            encoding="utf-8",
        )

        manifest = tmp_path / "manifest_{}.tsv".format(provider)
        make_manifest(
            manifest,
            [
                {
                    "provider": provider,
                    "id": provider_ids[provider],
                    "species_key": species_key,
                    "cds_url": to_file_url(cds_source),
                    "gff_url": to_file_url(gff_source),
                    "genome_url": to_file_url(genome_source),
                    "cds_filename": species_key + ".cds.fa",
                    "gff_filename": species_key + ".gene.gff3",
                    "genome_filename": species_key + ".genome.fa",
                }
            ],
        )

        download_dir = tmp_path / "download_cache" / provider
        out_cds = tmp_path / "out_cds" / provider
        out_gff = tmp_path / "out_gff" / provider
        out_genome = tmp_path / "out_genome" / provider
        completed = run_script(
            "--provider",
            "all",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(download_dir),
            "--species-cds-dir",
            str(out_cds),
            "--species-gff-dir",
            str(out_gff),
            "--species-genome-dir",
            str(out_genome),
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

        raw_dir = download_dir / provider_layout[provider] / "species_wise_original" / species_key
        assert (raw_dir / (species_key + ".cds.fa")).exists()
        assert (raw_dir / (species_key + ".gene.gff3")).exists()
        assert (raw_dir / (species_key + ".genome.fa")).exists()

        formatted_cds = out_cds / (species_key + "_cds.fa.gz")
        formatted_gff = out_gff / (species_key + "_gene.gff3")
        formatted_genome = out_genome / (species_key + "_genome.fa.gz")
        assert formatted_cds.exists()
        assert formatted_gff.exists()
        assert formatted_genome.exists()

        with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
            cds_text = handle.read()
        assert cds_text.count(">Arabidopsis_thaliana_AT1G01010") == 1
        assert "ATGAAATTT" in cds_text


def test_download_manifest_resolves_coge_urls_from_id_without_templates(tmp_path):
    class _CoGeFixtureHandler(SimpleHTTPRequestHandler):
        def _send_json(self, payload):
            body = json.dumps(payload).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def _send_text(self, payload):
            body = payload.encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "text/plain; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            parsed = urlparse(self.path)
            if parsed.path == "/coge/api/v1/genomes/search/arabidopsis_thaliana_v1/":
                self._send_json(
                    {
                        "genomes": [
                            {
                                "id": 24739,
                                "name": "Arabidopsis thaliana v10, 50x masked",
                                "certified": False,
                                "deleted": False,
                                "organism": {"name": "Arabidopsis thaliana"},
                            }
                        ]
                    }
                )
                return
            if parsed.path == "/coge/GenomeInfo.pl":
                self._send_json(
                    {
                        "file": "Arabidopsis_thaliana.gid24739.gff",
                        "files": [
                            "http://127.0.0.1:{}/files/Arabidopsis_thaliana.gid24739.gff".format(
                                self.server.server_port
                            )
                        ],
                    }
                )
                return
            if parsed.path == "/files/Arabidopsis_thaliana.gid24739.gff":
                self._send_text("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=AT1G01010\n")
                return
            if parsed.path == "/coge/get_seqs_for_feattype_for_genome.pl":
                self._send_text(">AT1G01010.1\nATGAAATTT\n")
                return
            if parsed.path == "/coge/api/v1/genomes/24739/sequence":
                self._send_text(">chr1\nATGCATGC\n")
                return
            self.send_response(404)
            self.end_headers()

    server = ThreadingHTTPServer(("127.0.0.1", 0), _CoGeFixtureHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "coge",
                    "id": "24739 (Arabidopsis thaliana)",
                    "species_key": "Arabidopsis_thaliana",
                    "cds_url": "",
                    "gff_url": "",
                    "genome_url": "",
                    "cds_filename": "",
                    "gff_filename": "",
                    "genome_filename": "",
                }
            ],
        )

        download_dir = tmp_path / "download_cache"
        env = dict(os.environ)
        env["GG_COGE_API_BASE_URL"] = "http://127.0.0.1:{}/coge/api/v1".format(server.server_port)
        env["GG_COGE_WEB_BASE_URL"] = "http://127.0.0.1:{}/coge".format(server.server_port)
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--provider",
                "coge",
                "--download-manifest",
                str(manifest),
                "--download-dir",
                str(download_dir),
                "--download-only",
            ],
            capture_output=True,
            text=True,
            check=False,
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

        raw_dir = download_dir / "CoGe" / "species_wise_original" / "Arabidopsis_thaliana"
        assert (raw_dir / "Arabidopsis_thaliana.coge.gid24739.cds.fa").exists()
        assert (raw_dir / "Arabidopsis_thaliana.gid24739.gff").exists()
        assert (raw_dir / "Arabidopsis_thaliana.coge.gid24739.genome.fa").exists()
    finally:
        server.shutdown()
        server.server_close()


def test_download_manifest_coge_requires_numeric_genome_id(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "coge",
                "id": "arabidopsis_thaliana_v1",
                "species_key": "Arabidopsis_thaliana",
                "cds_url": "",
                "gff_url": "",
                "genome_url": "",
                "cds_filename": "",
                "gff_filename": "",
                "genome_filename": "",
            }
        ],
    )

    completed = run_script(
        "--provider",
        "coge",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(tmp_path / "download_cache"),
        "--download-only",
    )
    assert completed.returncode == 1
    assert "requires numeric genome_id (gid)" in completed.stderr


def test_download_manifest_resolves_cngb_id_via_cnsa_then_ncbi(tmp_path):
    ftp_root = tmp_path / "ftp_root"
    ftp_dir = ftp_root / "genomes" / "all" / "GCF" / "000" / "001" / "405" / "GCF_000001405.40_GRCh38.p14"
    ftp_dir.mkdir(parents=True)

    with gzip.open(ftp_dir / "GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz", "wt", encoding="utf-8") as handle:
        handle.write(">NC_000001.11_cds_NP_000001.1_1 [gene=ABC1] [db_xref=GeneID:111]\nATGAA\n")
    with gzip.open(ftp_dir / "GCF_000001405.40_GRCh38.p14_genomic.gff.gz", "wt", encoding="utf-8") as handle:
        handle.write("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1\n")
    with gzip.open(ftp_dir / "GCF_000001405.40_GRCh38.p14_genomic.fna.gz", "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGCATGC\n")

    class _CngbNcbiFixtureHandler(SimpleHTTPRequestHandler):
        def __init__(self, *args, root_dir=None, **kwargs):
            self._root_dir = root_dir
            super().__init__(*args, directory=str(root_dir), **kwargs)

        def _send_json(self, payload):
            body = json.dumps(payload).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            parsed = urlparse(self.path)
            if parsed.path.startswith("/cngb/cnsa/ajax/assembly/public_view/"):
                query = parse_qs(parsed.query)
                if (query.get("q") or [""])[0] == "CNA0012345":
                    self._send_json(
                        {
                            "code": 0,
                            "data": {
                                "summary_data": {
                                    "accession_id": "CNA0012345",
                                    "refseq_assembly_accession": "GCF_000001405.40",
                                    "organism": {"name": "Homo sapiens"},
                                }
                            },
                        }
                    )
                    return
                self._send_json({"code": 2, "error": {"code": 2001, "content": "invalid parameters"}})
                return
            if parsed.path.startswith("/eutils/esearch.fcgi"):
                self._send_json(
                    {
                        "header": {"type": "esearch", "version": "0.3"},
                        "esearchresult": {"idlist": ["12345"]},
                    }
                )
                return
            if parsed.path.startswith("/eutils/esummary.fcgi"):
                self._send_json(
                    {
                        "header": {"type": "esummary", "version": "0.3"},
                        "result": {
                            "uids": ["12345"],
                            "12345": {
                                "ftppath_refseq": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14",
                                "organism": "Homo sapiens (human)",
                                "speciesname": "Homo sapiens",
                            },
                        },
                    }
                )
                return
            super().do_GET()

    handler = lambda *args, **kwargs: _CngbNcbiFixtureHandler(*args, root_dir=ftp_root, **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "cngb",
                    "id": "cngb:CNA0012345",
                    "species_key": "",
                    "cds_url": "",
                    "gff_url": "",
                    "genome_url": "",
                    "cds_filename": "",
                    "gff_filename": "",
                    "genome_filename": "",
                }
            ],
        )
        download_dir = tmp_path / "download_cache"
        env = dict(os.environ)
        env["GG_CNGB_CNSA_BASE_URL"] = "http://127.0.0.1:{}/cngb/cnsa/ajax".format(server.server_port)
        env["GG_NCBI_EUTILS_BASE_URL"] = "http://127.0.0.1:{}/eutils".format(server.server_port)
        env["GG_NCBI_FTP_BASE_URL"] = "http://127.0.0.1:{}".format(server.server_port)

        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--provider",
                "cngb",
                "--download-manifest",
                str(manifest),
                "--download-dir",
                str(download_dir),
                "--download-only",
            ],
            capture_output=True,
            text=True,
            check=False,
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

        raw_dir = download_dir / "CNGB" / "species_wise_original" / "Homo_sapiens"
        assert (raw_dir / "GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz").exists()
        assert (raw_dir / "GCF_000001405.40_GRCh38.p14_genomic.gff.gz").exists()
        assert (raw_dir / "GCF_000001405.40_GRCh38.p14_genomic.fna.gz").exists()
    finally:
        server.shutdown()
        server.server_close()


def test_download_manifest_recovers_stale_lock_file(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    cds_source = source_dir / "coge_cds.fa"
    gff_source = source_dir / "coge_gene.gff3"
    cds_source.write_text(">AT1G01010_t1\nATG\n", encoding="utf-8")
    gff_source.write_text("chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=AT1G01010\n", encoding="utf-8")

    manifest = tmp_path / "manifest.tsv"
    species_key = "Arabidopsis_thaliana"
    cds_name = "Arabidopsis_thaliana.cds.fa"
    gff_name = "Arabidopsis_thaliana.gene.gff3"
    make_manifest(
        manifest,
        [
            {
                "provider": "coge",
                "id": "24739",
                "species_key": species_key,
                "cds_url": to_file_url(cds_source),
                "gff_url": to_file_url(gff_source),
                "cds_filename": cds_name,
                "gff_filename": gff_name,
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    raw_dir = download_dir / "CoGe" / "species_wise_original" / species_key
    raw_dir.mkdir(parents=True)
    stale_lock = raw_dir / (cds_name + ".lock")
    stale_lock.write_text("999999\t0\tstale\n", encoding="utf-8")

    completed = run_script(
        "--provider",
        "coge",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(download_dir),
        "--download-only",
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert (raw_dir / cds_name).exists()
    assert (raw_dir / gff_name).exists()
    assert "recovered stale lock" in completed.stderr.lower()


def test_download_manifest_resolves_urls_from_id_templates_for_non_ncbi(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    source_id_token = "24739"
    cds_source = source_dir / (source_id_token + ".cds.fasta")
    gff_source = source_dir / (source_id_token + ".genes.gff3")
    cds_source.write_text(">AT1G01010_t1\nATG\n", encoding="utf-8")
    gff_source.write_text("chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=AT1G01010\n", encoding="utf-8")

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "coge",
                "id": source_id_token,
                "species_key": "Arabidopsis_thaliana",
                "cds_url": "",
                "gff_url": "",
                "cds_filename": "Arabidopsis_thaliana.cds.fa",
                "gff_filename": "Arabidopsis_thaliana.gene.gff3",
                "cds_url_template": source_dir.resolve().as_uri() + "/{id}.cds.fasta",
                "gff_url_template": source_dir.resolve().as_uri() + "/{id}.genes.gff3",
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    completed = run_script(
        "--provider",
        "coge",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(download_dir),
        "--download-only",
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    raw_dir = download_dir / "CoGe" / "species_wise_original" / "Arabidopsis_thaliana"
    assert (raw_dir / "Arabidopsis_thaliana.cds.fa").exists()
    assert (raw_dir / "Arabidopsis_thaliana.gene.gff3").exists()


def test_download_manifest_resolves_urls_from_id_for_all_non_ncbi_providers_via_env_templates(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()

    rows = []
    provider_species = {
        "ensembl": "Homo_sapiens",
        "ensemblplants": "Ostreococcus_lucimarinus",
        "coge": "Arabidopsis_thaliana",
        "cngb": "Arabidopsis_thaliana",
        "flybase": "Drosophila_melanogaster",
        "wormbase": "Caenorhabditis_elegans",
        "vectorbase": "Anopheles_gambiae",
    }
    provider_ids = {
        "ensembl": "homo_sapiens",
        "ensemblplants": "Ostreococcus_lucimarinus",
        "coge": "24739",
        "cngb": "cngb:arabidopsis_thaliana_v1",
        "flybase": "dmel_r6.61",
        "wormbase": "celegans_prjna13758_ws290",
        "vectorbase": "anopheles_gambiae_pest",
    }

    for provider, source_id in provider_ids.items():
        token = source_id.split(":", 1)[-1] if ":" in source_id else source_id
        cds_source = source_dir / (token + ".cds.fa")
        gff_source = source_dir / (token + ".gene.gff3")
        genome_source = source_dir / (token + ".genome.fa")
        cds_source.write_text(
            ">geneA.t1\nATGAA\n>geneA.t2\nATGAAATTT\n",
            encoding="utf-8",
        )
        gff_source.write_text("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=geneA\n", encoding="utf-8")
        genome_source.write_text(">chr1\nATGCATGC\n", encoding="utf-8")
        rows.append(
            {
                "provider": provider,
                "id": source_id,
                "species_key": provider_species[provider],
                "cds_url": "",
                "gff_url": "",
                "genome_url": "",
                "cds_filename": "",
                "gff_filename": "",
                "genome_filename": "",
            }
        )

    manifest = tmp_path / "manifest.tsv"
    make_manifest(manifest, rows)

    env = dict(os.environ)
    provider_to_env = {
        "ensembl": "GG_ENSEMBL",
        "ensemblplants": "GG_ENSEMBLPLANTS",
        "coge": "GG_COGE",
        "cngb": "GG_CNGB",
        "flybase": "GG_FLYBASE",
        "wormbase": "GG_WORMBASE",
        "vectorbase": "GG_VECTORBASE",
    }
    for provider, env_prefix in provider_to_env.items():
        env[env_prefix + "_CDS_URL_TEMPLATE"] = source_dir.resolve().as_uri() + "/{id}.cds.fa"
        env[env_prefix + "_GFF_URL_TEMPLATE"] = source_dir.resolve().as_uri() + "/{id}.gene.gff3"
        env[env_prefix + "_GENOME_URL_TEMPLATE"] = source_dir.resolve().as_uri() + "/{id}.genome.fa"

    download_dir = tmp_path / "download_cache"
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--provider",
            "all",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(download_dir),
            "--download-only",
        ],
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    expected_roots = {
        "ensembl": download_dir / "Ensembl" / "original_files",
        "ensemblplants": download_dir / "20230216_EnsemblPlants" / "original_files",
        "coge": download_dir / "CoGe" / "species_wise_original" / provider_species["coge"],
        "cngb": download_dir / "CNGB" / "species_wise_original" / provider_species["cngb"],
        "flybase": download_dir / "FlyBase" / "species_wise_original" / provider_species["flybase"],
        "wormbase": download_dir / "WormBase" / "species_wise_original" / provider_species["wormbase"],
        "vectorbase": download_dir / "VectorBase" / "species_wise_original" / provider_species["vectorbase"],
    }
    for provider, root in expected_roots.items():
        assert root.exists(), provider
        files = sorted([p.name for p in root.iterdir() if p.is_file()])
        assert any(name.endswith(".cds.fa") for name in files), provider
        assert any(name.endswith(".gene.gff3") for name in files), provider
        assert any(name.endswith(".genome.fa") for name in files), provider


def test_download_manifest_ncbi_id_only_auto_resolve(tmp_path):
    ftp_root = tmp_path / "ftp_root"
    ftp_dir = ftp_root / "genomes" / "all" / "GCF" / "000" / "001" / "405" / "GCF_000001405.40_GRCh38.p14"
    ftp_dir.mkdir(parents=True)

    cds_content = ">NC_000001.11_cds_NP_000001.1_1\nATGAA\n"
    gff_content = "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1\n"
    with gzip.open(ftp_dir / "GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz", "wt", encoding="utf-8") as handle:
        handle.write(cds_content)
    with gzip.open(ftp_dir / "GCF_000001405.40_GRCh38.p14_genomic.gff.gz", "wt", encoding="utf-8") as handle:
        handle.write(gff_content)
    with gzip.open(ftp_dir / "GCF_000001405.40_GRCh38.p14_genomic.fna.gz", "wt", encoding="utf-8") as handle:
        handle.write(">NC_000001.11 chromosome 1\nATGCATGC\n")

    handler = lambda *args, **kwargs: _NcbiFixtureHandler(*args, root_dir=ftp_root, **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "ncbi",
                    "id": "GCF_000001405.40 (Homo sapiens)",
                    "species_key": "",
                    "cds_url": "",
                    "gff_url": "",
                    "cds_filename": "",
                    "gff_filename": "",
                }
            ],
        )

        download_dir = tmp_path / "download_cache"
        out_cds = tmp_path / "out_cds"
        out_gff = tmp_path / "out_gff"
        out_genome = tmp_path / "out_genome"
        env = dict(os.environ)
        env["GG_NCBI_EUTILS_BASE_URL"] = "http://127.0.0.1:{}/eutils".format(server.server_port)
        env["GG_NCBI_FTP_BASE_URL"] = "http://127.0.0.1:{}".format(server.server_port)
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--provider",
                "ncbi",
                "--download-manifest",
                str(manifest),
                "--download-dir",
                str(download_dir),
                "--species-cds-dir",
                str(out_cds),
                "--species-gff-dir",
                str(out_gff),
                "--species-genome-dir",
                str(out_genome),
            ],
            capture_output=True,
            text=True,
            check=False,
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

        raw_dir = download_dir / "NCBI_Genome" / "species_wise_original" / "Homo_sapiens"
        assert (raw_dir / "GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz").exists()
        assert (raw_dir / "GCF_000001405.40_GRCh38.p14_genomic.gff.gz").exists()
        assert (raw_dir / "GCF_000001405.40_GRCh38.p14_genomic.fna.gz").exists()

        formatted_cds = out_cds / "Homo_sapiens_GCF_000001405.40_GRCh38.p14_cds_from_genomic.fa.gz"
        formatted_gff = out_gff / "Homo_sapiens_GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
        formatted_genome = out_genome / "Homo_sapiens_GCF_000001405.40_GRCh38.p14_genomic.fa.gz"
        assert formatted_cds.exists()
        assert formatted_gff.exists()
        assert formatted_genome.exists()
        with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
            cds_text = handle.read()
        assert ">Homo_sapiens_NC_000001.11_cds_NP_000001.1_1" in cds_text
    finally:
        server.shutdown()
        server.server_close()


def test_download_manifest_refseq_and_genbank_id_only_auto_resolve(tmp_path):
    ftp_root = tmp_path / "ftp_root"
    refseq_dir = ftp_root / "genomes" / "all" / "GCF" / "000" / "001" / "405" / "GCF_000001405.40_GRCh38.p14"
    genbank_dir = ftp_root / "genomes" / "all" / "GCA" / "000" / "001" / "405" / "GCA_000001405.29_GRCh38.p14"
    refseq_dir.mkdir(parents=True)
    genbank_dir.mkdir(parents=True)

    with gzip.open(refseq_dir / "GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz", "wt", encoding="utf-8") as handle:
        handle.write(">NC_000001.11_cds_NP_refseq_1.1_1\nATGAA\n")
    with gzip.open(refseq_dir / "GCF_000001405.40_GRCh38.p14_genomic.gff.gz", "wt", encoding="utf-8") as handle:
        handle.write("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene_refseq\n")
    with gzip.open(refseq_dir / "GCF_000001405.40_GRCh38.p14_genomic.fna.gz", "wt", encoding="utf-8") as handle:
        handle.write(">NC_000001.11 chromosome 1 refseq\nATGCATGC\n")

    with gzip.open(genbank_dir / "GCA_000001405.29_GRCh38.p14_cds_from_genomic.fna.gz", "wt", encoding="utf-8") as handle:
        handle.write(">NC_000001.11_cds_NP_genbank_1.1_1\nATGAA\n")
    with gzip.open(genbank_dir / "GCA_000001405.29_GRCh38.p14_genomic.gff.gz", "wt", encoding="utf-8") as handle:
        handle.write("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene_genbank\n")
    with gzip.open(genbank_dir / "GCA_000001405.29_GRCh38.p14_genomic.fna.gz", "wt", encoding="utf-8") as handle:
        handle.write(">NC_000001.11 chromosome 1 genbank\nATGCATGC\n")

    class _NcbiRefseqGenbankFixtureHandler(SimpleHTTPRequestHandler):
        def __init__(self, *args, root_dir=None, **kwargs):
            self._root_dir = root_dir
            super().__init__(*args, directory=str(root_dir), **kwargs)

        def _send_json(self, payload):
            body = json.dumps(payload).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            if self.path.startswith("/eutils/esearch.fcgi"):
                self._send_json(
                    {
                        "header": {"type": "esearch", "version": "0.3"},
                        "esearchresult": {"idlist": ["12345"]},
                    }
                )
                return
            if self.path.startswith("/eutils/esummary.fcgi"):
                self._send_json(
                    {
                        "header": {"type": "esummary", "version": "0.3"},
                        "result": {
                            "uids": ["12345"],
                            "12345": {
                                "ftppath_refseq": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14",
                                "ftppath_genbank": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14",
                                "organism": "Homo sapiens (human)",
                                "speciesname": "Homo sapiens",
                            },
                        },
                    }
                )
                return
            super().do_GET()

    handler = lambda *args, **kwargs: _NcbiRefseqGenbankFixtureHandler(*args, root_dir=ftp_root, **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        env = dict(os.environ)
        env["GG_NCBI_EUTILS_BASE_URL"] = "http://127.0.0.1:{}/eutils".format(server.server_port)
        env["GG_NCBI_FTP_BASE_URL"] = "http://127.0.0.1:{}".format(server.server_port)

        for provider, expected_dir, expected_stem in (
            ("refseq", "NCBI_RefSeq", "GCF_000001405.40_GRCh38.p14"),
            ("genbank", "NCBI_GenBank", "GCA_000001405.29_GRCh38.p14"),
        ):
            manifest = tmp_path / ("manifest_{}.tsv".format(provider))
            make_manifest(
                manifest,
                [
                    {
                        "provider": provider,
                        "id": "GCF_000001405.40 (Homo sapiens)",
                        "species_key": "",
                        "cds_url": "",
                        "gff_url": "",
                        "genome_url": "",
                        "cds_filename": "",
                        "gff_filename": "",
                        "genome_filename": "",
                    }
                ],
            )

            download_dir = tmp_path / ("download_cache_{}".format(provider))
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    "--provider",
                    provider,
                    "--download-manifest",
                    str(manifest),
                    "--download-dir",
                    str(download_dir),
                    "--download-only",
                ],
                capture_output=True,
                text=True,
                check=False,
                env=env,
            )
            assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

            raw_dir = download_dir / expected_dir / "species_wise_original" / "Homo_sapiens"
            assert (raw_dir / (expected_stem + "_cds_from_genomic.fna.gz")).exists()
            assert (raw_dir / (expected_stem + "_genomic.gff.gz")).exists()
            assert (raw_dir / (expected_stem + "_genomic.fna.gz")).exists()
    finally:
        server.shutdown()
        server.server_close()


def test_download_manifest_ncbi_dataset_url_id_auto_resolve(tmp_path):
    ftp_root = tmp_path / "ftp_root"
    ftp_dir = ftp_root / "genomes" / "all" / "GCF" / "000" / "001" / "405" / "GCF_000001405.40_GRCh38.p14"
    ftp_dir.mkdir(parents=True)

    cds_content = ">NC_000001.11_cds_NP_000001.1_1\nATGAA\n"
    gff_content = "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1\n"
    with gzip.open(ftp_dir / "GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz", "wt", encoding="utf-8") as handle:
        handle.write(cds_content)
    with gzip.open(ftp_dir / "GCF_000001405.40_GRCh38.p14_genomic.gff.gz", "wt", encoding="utf-8") as handle:
        handle.write(gff_content)
    with gzip.open(ftp_dir / "GCF_000001405.40_GRCh38.p14_genomic.fna.gz", "wt", encoding="utf-8") as handle:
        handle.write(">NC_000001.11 chromosome 1\nATGCATGC\n")

    handler = lambda *args, **kwargs: _NcbiFixtureHandler(*args, root_dir=ftp_root, **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        manifest.write_text(
            "provider\tid\tspecies_key\nncbi\thttps://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/\tHomo_sapiens\n",
            encoding="utf-8",
        )

        download_dir = tmp_path / "download_cache"
        env = dict(os.environ)
        env["GG_NCBI_EUTILS_BASE_URL"] = "http://127.0.0.1:{}/eutils".format(server.server_port)
        env["GG_NCBI_FTP_BASE_URL"] = "http://127.0.0.1:{}".format(server.server_port)
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--provider",
                "all",
                "--download-manifest",
                str(manifest),
                "--download-dir",
                str(download_dir),
                "--download-only",
            ],
            capture_output=True,
            text=True,
            check=False,
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

        raw_dir = download_dir / "NCBI_Genome" / "species_wise_original" / "Homo_sapiens"
        assert (raw_dir / "GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz").exists()
        assert (raw_dir / "GCF_000001405.40_GRCh38.p14_genomic.gff.gz").exists()
        assert (raw_dir / "GCF_000001405.40_GRCh38.p14_genomic.fna.gz").exists()
    finally:
        server.shutdown()
        server.server_close()


def test_download_manifest_ncbi_falls_back_to_datasets_api_when_ftp_files_missing(tmp_path):
    ftp_root = tmp_path / "ftp_root"
    ftp_root.mkdir()

    cds_content = ">lcl|ctg1_cds_XP_1.1_1 [gene=catA] [db_xref=GeneID:765915]\nATGAA\n"
    gff_content = "ctg1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1\n"
    genome_content = ">ctg1\nATGCATGC\n"

    class _NcbiFallbackHandler(SimpleHTTPRequestHandler):
        def __init__(self, *args, root_dir=None, **kwargs):
            self._root_dir = root_dir
            super().__init__(*args, directory=str(root_dir), **kwargs)

        def _send_json(self, payload):
            body = json.dumps(payload).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def _send_zip(self, files):
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, mode="w", compression=zipfile.ZIP_DEFLATED) as archive:
                archive.writestr("README.md", "fixture\n")
                archive.writestr("ncbi_dataset/data/assembly_data_report.jsonl", "{}\n")
                archive.writestr("ncbi_dataset/data/dataset_catalog.json", "{}\n")
                for name, text in files.items():
                    archive.writestr(name, text)
            payload = buf.getvalue()
            self.send_response(200)
            self.send_header("Content-Type", "application/zip")
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)

        def do_GET(self):
            if self.path.startswith("/eutils/esearch.fcgi"):
                self._send_json(
                    {
                        "header": {"type": "esearch", "version": "0.3"},
                        "esearchresult": {"idlist": ["765915"]},
                    }
                )
                return
            if self.path.startswith("/eutils/esummary.fcgi"):
                self._send_json(
                    {
                        "header": {"type": "esummary", "version": "0.3"},
                        "result": {
                            "uids": ["765915"],
                            "765915": {
                                "ftppath_refseq": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/102/555/GCA_002102555.1_Catan2",
                                "organism": "Catenaria anguillulae PL171",
                                "speciesname": "Catenaria anguillulae",
                            },
                        },
                    }
                )
                return
            if self.path.startswith("/datasets/genome/accession/GCA_002102555.1/download"):
                query = parse_qs(urlparse(self.path).query)
                include_type = query.get("include_annotation_type", [""])[0]
                if include_type == "CDS_FASTA":
                    files = {"ncbi_dataset/data/GCA_002102555.1/cds_from_genomic.fna": cds_content}
                elif include_type == "GENOME_GFF":
                    files = {"ncbi_dataset/data/GCA_002102555.1/genomic.gff": gff_content}
                elif include_type == "GENOME_FASTA":
                    files = {"ncbi_dataset/data/GCA_002102555.1/GCA_002102555.1_Catan2_genomic.fna": genome_content}
                else:
                    self.send_error(400, "invalid include_annotation_type")
                    return
                self._send_zip(files)
                return
            super().do_GET()

    handler = lambda *args, **kwargs: _NcbiFallbackHandler(*args, root_dir=ftp_root, **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "ncbi",
                    "id": "GCA_002102555.1",
                    "species_key": "",
                    "cds_url": "",
                    "gff_url": "",
                    "genome_url": "",
                    "cds_filename": "",
                    "gff_filename": "",
                    "genome_filename": "",
                }
            ],
        )

        download_dir = tmp_path / "download_cache"
        env = dict(os.environ)
        env["GG_NCBI_EUTILS_BASE_URL"] = "http://127.0.0.1:{}/eutils".format(server.server_port)
        env["GG_NCBI_FTP_BASE_URL"] = "http://127.0.0.1:{}".format(server.server_port)
        env["GG_NCBI_DATASETS_BASE_URL"] = "http://127.0.0.1:{}/datasets".format(server.server_port)

        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--provider",
                "ncbi",
                "--download-manifest",
                str(manifest),
                "--download-dir",
                str(download_dir),
                "--download-only",
            ],
            capture_output=True,
            text=True,
            check=False,
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

        raw_dir = download_dir / "NCBI_Genome" / "species_wise_original" / "Catenaria_anguillulae"
        cds_path = raw_dir / "GCA_002102555.1_Catan2_cds_from_genomic.fna.gz"
        gff_path = raw_dir / "GCA_002102555.1_Catan2_genomic.gff.gz"
        genome_path = raw_dir / "GCA_002102555.1_Catan2_genomic.fna.gz"
        assert cds_path.exists()
        assert gff_path.exists()
        assert genome_path.exists()

        with gzip.open(cds_path, "rt", encoding="utf-8") as handle:
            assert "GeneID:765915" in handle.read()
        with gzip.open(gff_path, "rt", encoding="utf-8") as handle:
            assert "ID=gene1" in handle.read()
        with gzip.open(genome_path, "rt", encoding="utf-8") as handle:
            assert ">ctg1" in handle.read()
        assert "fallback via NCBI Datasets API" in completed.stderr
    finally:
        server.shutdown()
        server.server_close()


def test_download_manifest_ncbi_gene_level_aggregate_keeps_one_longest_per_gene(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    cds_source = source_dir / "ncbi_cds.fna.gz"
    gff_source = source_dir / "ncbi_genomic.gff.gz"
    genome_source = source_dir / "ncbi_genomic.fna.gz"
    with gzip.open(cds_source, "wt", encoding="utf-8") as handle:
        handle.write(
            ">lcl|NC_000001.11_cds_XP_1.1_1 [gene=ABC1] [db_xref=GeneID:111]\n"
            "ATGAA\n"
            ">lcl|NC_000001.11_cds_XP_1.2_2 [gene=ABC1] [db_xref=GeneID:111]\n"
            "ATGAAATTT\n"
        )
    with gzip.open(gff_source, "wt", encoding="utf-8") as handle:
        handle.write("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1\n")
    with gzip.open(genome_source, "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGCATGC\n")

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "ncbi",
                "species_key": "Homo_sapiens",
                "id": "GCF_000001405.40",
                "cds_url": to_file_url(cds_source),
                "gff_url": to_file_url(gff_source),
                "genome_url": to_file_url(genome_source),
                "cds_filename": "GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz",
                "gff_filename": "GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
                "genome_filename": "GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
            }
        ],
    )

    out_cds = tmp_path / "out_cds"
    out_gff = tmp_path / "out_gff"
    out_genome = tmp_path / "out_genome"
    completed = run_script(
        "--provider",
        "ncbi",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(tmp_path / "download_cache"),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Homo_sapiens_GCF_000001405.40_GRCh38.p14_cds_from_genomic.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert text.count(">Homo_sapiens_GeneID111") == 1
    assert "ATGAAATTT" in text


def test_download_manifest_ncbi_gene_level_aggregate_prefers_ensembl_gene_id(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    cds_source = source_dir / "ncbi_cds.fna.gz"
    gff_source = source_dir / "ncbi_genomic.gff.gz"
    genome_source = source_dir / "ncbi_genomic.fna.gz"
    with gzip.open(cds_source, "wt", encoding="utf-8") as handle:
        handle.write(
            ">lcl|NC_000001.11_cds_XP_1.1_1 [gene=ABC1] [db_xref=GeneID:111,Ensembl:ENSG00000111111]\n"
            "ATGAA\n"
            ">lcl|NC_000001.11_cds_XP_1.2_2 [gene=ABC1] [db_xref=GeneID:111,Ensembl:ENSG00000111111]\n"
            "ATGAAATTT\n"
        )
    with gzip.open(gff_source, "wt", encoding="utf-8") as handle:
        handle.write("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1\n")
    with gzip.open(genome_source, "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGCATGC\n")

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "ncbi",
                "species_key": "Homo_sapiens",
                "id": "GCF_000001405.40",
                "cds_url": to_file_url(cds_source),
                "gff_url": to_file_url(gff_source),
                "genome_url": to_file_url(genome_source),
                "cds_filename": "GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz",
                "gff_filename": "GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
                "genome_filename": "GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
            }
        ],
    )

    out_cds = tmp_path / "out_cds"
    out_gff = tmp_path / "out_gff"
    out_genome = tmp_path / "out_genome"
    completed = run_script(
        "--provider",
        "ncbi",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(tmp_path / "download_cache"),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Homo_sapiens_GCF_000001405.40_GRCh38.p14_cds_from_genomic.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert text.count(">Homo_sapiens_ENSG00000111111") == 1
    assert ">Homo_sapiens_GeneID111" not in text
    assert "ATGAAATTT" in text
