from pathlib import Path
import bz2
import csv
import gzip
import io
import json
import os
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import socket
import subprocess
import sys
import tarfile
import threading
import time
from urllib.parse import parse_qs, urlparse
import zipfile

from openpyxl import Workbook


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "format_species_inputs.py"
SMALL_DATASET_ROOT = Path(__file__).resolve().parent / "data" / "small_gfe_dataset"


def run_script(*args, env=None):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), *args],
        capture_output=True,
        text=True,
        check=False,
        env=env,
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
                "gbff_url",
                "genome_url",
                "cds_archive_member",
                "gff_archive_member",
                "gbff_archive_member",
                "genome_archive_member",
                "cds_filename",
                "gff_filename",
                "gbff_filename",
                "genome_filename",
                "cds_url_template",
                "gff_url_template",
                "genome_url_template",
                "local_cds_path",
                "local_gff_path",
                "local_gbff_path",
                "local_genome_path",
                "fernbase_confidence_mode",
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


def _current_boot_id():
    boot_id_path = Path("/proc/sys/kernel/random/boot_id")
    try:
        if boot_id_path.is_file():
            return boot_id_path.read_text(encoding="utf-8").strip()
    except OSError:
        pass
    completed = subprocess.run(
        ["sysctl", "-n", "kern.bootsessionuuid"],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        return ""
    return completed.stdout.strip()


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


class _VEuPathDbFixtureHandler(SimpleHTTPRequestHandler):
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
        if self.path.startswith("/veupathdb/service/record-types/organism/searches/GenomeDataTypes/reports/standard"):
            base = "http://127.0.0.1:{}".format(self.server.server_port)
            self._send_json(
                {
                    "records": [
                        {
                            "attributes": {
                                "primary_key": "Entamoeba nuttalli P19",
                                "species": "Entamoeba nuttalli",
                                "project_id": "AmoebaDB",
                                "URLGenomeFasta": (
                                    base
                                    + "/common/downloads/Current_Release/EnuttalliP19/fasta/data/"
                                    "AmoebaDB-68_EnuttalliP19_Genome.fasta"
                                ),
                                "URLgff": (
                                    base
                                    + "/common/downloads/Current_Release/EnuttalliP19/gff/data/"
                                    "AmoebaDB-68_EnuttalliP19.gff"
                                ),
                                "URLproteinFasta": (
                                    base
                                    + "/common/downloads/Current_Release/EnuttalliP19/fasta/data/"
                                    "AmoebaDB-68_EnuttalliP19_AnnotatedProteins.fasta"
                                ),
                            }
                        }
                    ]
                }
            )
            return
        super().do_GET()


class _InsectBaseFixtureHandler(SimpleHTTPRequestHandler):
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
        if self.path.startswith("/api/genome/genomes/IBG_00001/"):
            self._send_json(
                {
                    "ibg_id": "IBG_00001",
                    "species": "Abrostola tripartita",
                }
            )
            return
        super().do_GET()


def write_tar_bz2(path, members):
    with tarfile.open(path, "w:bz2") as archive:
        for name, content in members.items():
            payload = str(content).encode("utf-8")
            info = tarfile.TarInfo(name=name)
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))


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
    formatted_gff = out_gff / "Ostreococcus_lucimarinus_ASM9206v1.56.gff.gz"
    formatted_genome = out_genome / "Ostreococcus_lucimarinus_ASM9206v1.dna.toplevel.fa.gz"
    assert formatted_cds.exists()
    assert formatted_gff.exists()
    assert formatted_genome.exists()

    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert ">Ostreococcus_lucimarinus_OSTLU_25062" in cds_text
    assert "ATGAAN" in cds_text
    with gzip.open(formatted_gff, "rt", encoding="utf-8") as handle:
        formatted_gff_text = handle.read()
    assert "evm.model." not in formatted_gff_text
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
    formatted_gff = out_gff / "Hydrocotyle_leucocephala_HleucocephalaHAP1_768_v2.1.gene.gff.gz"
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
        formatted_gff = out_gff / (species_key + "_gene.gff.gz")
        formatted_genome = out_genome / (species_key + "_genome.fa.gz")
        assert formatted_cds.exists()
        assert formatted_gff.exists()
        assert formatted_genome.exists()

        with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
            cds_text = handle.read()
        assert cds_text.count(">Arabidopsis_thaliana_AT1G01010") == 1
        assert "ATGAAATTT" in cds_text


def test_download_manifest_supports_direct_with_explicit_urls(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    species_key = "Medicago_sativa"
    cds_source = source_dir / "medicago_sativa.cds.fa"
    gff_source = source_dir / "medicago_sativa.gene.gff3"
    genome_source = source_dir / "medicago_sativa.genome.fa"
    cds_source.write_text(">MsG1.t1\nATGAA\n>MsG1.t2\nATGAAATTT\n", encoding="utf-8")
    gff_source.write_text("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=MsG1\n", encoding="utf-8")
    genome_source.write_text(">chr1\nATGCATGC\n", encoding="utf-8")

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "direct",
                "id": "medicago_sativa_direct",
                "species_key": species_key,
                "cds_url": to_file_url(cds_source),
                "gff_url": to_file_url(gff_source),
                "genome_url": to_file_url(genome_source),
                "cds_filename": species_key + ".direct.cds.fa",
                "gff_filename": species_key + ".direct.gene.gff3",
                "genome_filename": species_key + ".direct.genome.fa",
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    out_cds = tmp_path / "out_cds"
    out_gff = tmp_path / "out_gff"
    out_genome = tmp_path / "out_genome"
    completed = run_script(
        "--provider",
        "direct",
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

    raw_dir = download_dir / "Direct" / "species_wise_original" / species_key
    assert (raw_dir / (species_key + ".direct.cds.fa")).exists()
    assert (raw_dir / (species_key + ".direct.gene.gff3")).exists()
    assert (raw_dir / (species_key + ".direct.genome.fa")).exists()

    formatted_cds = out_cds / (species_key + "_direct.cds.fa.gz")
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert cds_text.count(">Medicago_sativa_MsG1") == 1


def test_download_manifest_all_provider_only_scans_providers_declared_in_manifest_xlsx(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    species_key = "Medicago_sativa"
    cds_source = source_dir / "medicago_sativa.cds.fa"
    gff_source = source_dir / "medicago_sativa.gene.gff3"
    genome_source = source_dir / "medicago_sativa.genome.fa"
    cds_source.write_text(">MsG1.t1\nATGAA\n>MsG1.t2\nATGAAATTT\n", encoding="utf-8")
    gff_source.write_text("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=MsG1\n", encoding="utf-8")
    genome_source.write_text(">chr1\nATGCATGC\n", encoding="utf-8")

    manifest = tmp_path / "manifest.xlsx"
    headers = [
        "provider",
        "id",
        "species_key",
        "cds_url",
        "gff_url",
        "genome_url",
        "cds_filename",
        "gff_filename",
        "genome_filename",
    ]
    make_manifest_xlsx(
        manifest,
        headers,
        [
            {
                "provider": "direct",
                "id": "medicago_sativa_direct",
                "species_key": species_key,
                "cds_url": to_file_url(cds_source),
                "gff_url": to_file_url(gff_source),
                "genome_url": to_file_url(genome_source),
                "cds_filename": species_key + ".direct.cds.fa",
                "gff_filename": species_key + ".direct.gene.gff3",
                "genome_filename": species_key + ".direct.genome.fa",
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    out_cds = tmp_path / "out_cds"
    out_gff = tmp_path / "out_gff"
    out_genome = tmp_path / "out_genome"
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
    assert "input directory not found" not in completed.stderr

    raw_dir = download_dir / "Direct" / "species_wise_original" / species_key
    assert (raw_dir / (species_key + ".direct.cds.fa")).exists()
    assert (raw_dir / (species_key + ".direct.gene.gff3")).exists()
    assert (raw_dir / (species_key + ".direct.genome.fa")).exists()


def test_download_manifest_supports_direct_with_cds_only(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    species_key = "Croton_tiglium"
    cds_source = source_dir / "VVPY-SOAPdenovo-Trans-assembly.fa.gz"
    with gzip.open(cds_source, "wt", encoding="utf-8") as handle:
        handle.write(">ctg1.t1\nATGAA\n>ctg1.t2\nATGAAATTT\n")

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "direct",
                "id": "VVPY-Croton_tiglium",
                "species_key": species_key,
                "cds_url": to_file_url(cds_source),
                "cds_filename": species_key + ".cds.fa.gz",
                "gff_filename": "",
                "genome_filename": "",
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    out_cds = tmp_path / "out_cds"
    out_gff = tmp_path / "out_gff"
    out_genome = tmp_path / "out_genome"
    species_summary = tmp_path / "gg_input_generation_species.tsv"
    completed = run_script(
        "--provider",
        "direct",
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
        "--species-summary-output",
        str(species_summary),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    raw_dir = download_dir / "Direct" / "species_wise_original" / species_key
    assert (raw_dir / (species_key + ".cds.fa.gz")).exists()

    formatted_cds = out_cds / (species_key + "_cds.fa.gz")
    assert formatted_cds.exists()
    assert not any(out_gff.iterdir())
    assert not any(out_genome.iterdir())

    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["gff_status"] == "missing"
    assert rows[0]["genome_status"] == "missing"


def test_download_manifest_supports_direct_archive_members(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    species_key = "Lens_culinaris"
    archive_path = source_dir / "Lcu.2RBY.zip"
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(
            "Lcu.2RBY/Lcu.2RBY.cds.fasta",
            ">LcG1.t1\nATGAA\n>LcG1.t2\nATGAAATTT\n",
        )
        archive.writestr(
            "Lcu.2RBY/Lcu.2RBY.hc_genes.gff3",
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=LcG1\n",
        )
        archive.writestr(
            "Lcu.2RBY/Lcu.2RBY.fasta",
            ">chr1\nATGCATGC\n",
        )

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "direct",
                "id": "lens_culinaris_direct_archive",
                "species_key": species_key,
                "cds_url": to_file_url(archive_path),
                "gff_url": to_file_url(archive_path),
                "genome_url": to_file_url(archive_path),
                "cds_archive_member": "Lcu.2RBY/Lcu.2RBY.cds.fasta",
                "gff_archive_member": "Lcu.2RBY/Lcu.2RBY.hc_genes.gff3",
                "genome_archive_member": "Lcu.2RBY/Lcu.2RBY.fasta",
                "cds_filename": "",
                "gff_filename": "",
                "genome_filename": "",
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    out_cds = tmp_path / "out_cds"
    out_gff = tmp_path / "out_gff"
    out_genome = tmp_path / "out_genome"
    completed = run_script(
        "--provider",
        "direct",
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

    raw_dir = download_dir / "Direct" / "species_wise_original" / species_key
    assert (raw_dir / "Lcu.2RBY.cds.fasta").exists()
    assert (raw_dir / "Lcu.2RBY.hc_genes.gff3").exists()
    assert (raw_dir / "Lcu.2RBY.fasta").exists()

    formatted_files = list(out_cds.glob("*.cds.fa.gz"))
    assert len(formatted_files) == 1
    with gzip.open(formatted_files[0], "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert cds_text.count(">Lens_culinaris_LcG1") == 1


def test_download_manifest_derives_cds_when_manifest_has_only_gff_and_genome(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    species_key = "Arabidopsis_thaliana"
    gff_source = source_dir / "arabidopsis_thaliana.gene.gff3"
    genome_source = source_dir / "arabidopsis_thaliana.genome.fa"
    gff_source.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                "chr1\tsrc\tCDS\t7\t9\t.\t+\t0\tID=cds2;Parent=gene1.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_source.write_text(">chr1\nATGAAATTT\n", encoding="utf-8")

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "direct",
                "id": "arabidopsis_thaliana_direct",
                "species_key": species_key,
                "cds_url": "",
                "gff_url": to_file_url(gff_source),
                "genome_url": to_file_url(genome_source),
                "cds_filename": "",
                "gff_filename": species_key + ".direct.gene.gff3",
                "genome_filename": species_key + ".direct.genome.fa",
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    out_cds = tmp_path / "out_cds"
    out_gff = tmp_path / "out_gff"
    out_genome = tmp_path / "out_genome"
    completed = run_script(
        "--provider",
        "direct",
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

    raw_dir = download_dir / "Direct" / "species_wise_original" / species_key
    assert not any(path.name.endswith(".cds.fa") for path in raw_dir.iterdir())
    assert (raw_dir / (species_key + ".direct.gene.gff3")).exists()
    assert (raw_dir / (species_key + ".direct.genome.fa")).exists()

    formatted_cds = out_cds / (species_key + "_direct.gene.derived.cds.fa.gz")
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert cds_text.count(">Arabidopsis_thaliana_gene1") == 1
    assert "ATGTTT" in cds_text


def test_download_manifest_derives_cds_when_direct_urls_need_explicit_target_filenames(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    species_key = "Populus_tremula_x_Populus_alba"
    gff_source = source_dir / "sPta717_v1.1.gene.gff3.gz"
    genome_source = source_dir / "sPta717_v1.1.fa.gz"
    with gzip.open(gff_source, "wt", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1",
                    "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                    "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                    "chr1\tsrc\tCDS\t7\t9\t.\t+\t0\tID=cds2;Parent=gene1.t1",
                    "",
                ]
            )
        )
    with gzip.open(genome_source, "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGAAATTT\n")

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "direct",
                "id": "Populus_tremula_x_Populus_alba_AspenDB_sPta717_v1.1",
                "species_key": species_key,
                "cds_url": "",
                "gff_url": to_file_url(gff_source),
                "genome_url": to_file_url(genome_source),
                "cds_filename": "",
                "gff_filename": species_key + ".AspenDB.gene.gff3.gz",
                "genome_filename": species_key + ".AspenDB.genome.fa.gz",
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    out_cds = tmp_path / "out_cds"
    out_gff = tmp_path / "out_gff"
    out_genome = tmp_path / "out_genome"
    completed = run_script(
        "--provider",
        "direct",
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

    raw_dir = download_dir / "Direct" / "species_wise_original" / species_key
    assert (raw_dir / (species_key + ".AspenDB.gene.gff3.gz")).exists()
    assert (raw_dir / (species_key + ".AspenDB.genome.fa.gz")).exists()
    assert not (raw_dir / "sPta717_v1.1.fa.gz").exists()

    matches = list(out_cds.glob("*.derived.cds.fa.gz"))
    assert len(matches) == 1
    formatted_cds = matches[0]
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert cds_text.count(">Populus_tremula_gene1") == 1
    assert "ATGTTT" in cds_text


def test_download_manifest_derives_from_gbff_and_genome(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    species_key = "Moringa_oleifera"
    gbff_source = source_dir / "moringa_oleifera.genomic.gbff"
    genome_source = source_dir / "moringa_oleifera.genome.fa"
    gbff_source.write_text(
        "\n".join(
            [
                "LOCUS       chr1               9 bp    DNA     linear   PLN 01-JAN-2000",
                "DEFINITION  test.",
                "ACCESSION   chr1",
                "VERSION     chr1",
                "FEATURES             Location/Qualifiers",
                "     gene            1..9",
                "                     /locus_tag=\"MoG1\"",
                "                     /gene=\"MoG1\"",
                "     CDS             join(1..3,7..9)",
                "                     /locus_tag=\"MoG1\"",
                "                     /gene=\"MoG1\"",
                "                     /protein_id=\"MoG1.t1\"",
                "ORIGIN",
                "        1 atgaaattt",
                "//",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_source.write_text(">chr1\nATGAAATTT\n", encoding="utf-8")

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "direct",
                "id": "moringa_oleifera_direct",
                "species_key": species_key,
                "cds_url": "",
                "gff_url": "",
                "gbff_url": to_file_url(gbff_source),
                "genome_url": to_file_url(genome_source),
                "cds_filename": "",
                "gff_filename": "",
                "gbff_filename": species_key + ".direct.genomic.gbff",
                "genome_filename": species_key + ".direct.genome.fa",
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    out_cds = tmp_path / "out_cds"
    out_gff = tmp_path / "out_gff"
    out_genome = tmp_path / "out_genome"
    completed = run_script(
        "--provider",
        "direct",
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

    raw_dir = download_dir / "Direct" / "species_wise_original" / species_key
    assert (raw_dir / (species_key + ".direct.genomic.gbff")).exists()
    assert (raw_dir / (species_key + ".direct.genome.fa")).exists()

    formatted_cds = out_cds / (species_key + "_direct.genomic.derived.cds.fa.gz")
    formatted_gff = out_gff / (species_key + "_direct.genomic.derived.gff.gz")
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert cds_text.count(">Moringa_oleifera_MoG1") == 1
    assert "ATGTTT" in cds_text
    with gzip.open(formatted_gff, "rt", encoding="utf-8") as handle:
        gff_text = handle.read()
    assert "\tgene\t" in gff_text
    assert "\tCDS\t" in gff_text


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


def test_download_manifest_resolves_plantaedb_page_to_ncbi_bundle(tmp_path):
    server_root = tmp_path / "server_root"
    plantae_page = server_root / "taxa" / "species" / "erigeron-breviscapus"
    plantae_page.parent.mkdir(parents=True, exist_ok=True)

    ftp_dir = (
        server_root
        / "genomes"
        / "all"
        / "GCA"
        / "999"
        / "999"
        / "999"
        / "GCA_999999999.1_ASM999999v1"
    )
    ftp_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(ftp_dir / "GCA_999999999.1_ASM999999v1_cds_from_genomic.fna.gz", "wt", encoding="utf-8") as handle:
        handle.write(">gene1.t1\nATGAAATTT\n")
    with gzip.open(ftp_dir / "GCA_999999999.1_ASM999999v1_genomic.gff.gz", "wt", encoding="utf-8") as handle:
        handle.write("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1\n")
    with gzip.open(ftp_dir / "GCA_999999999.1_ASM999999v1_genomic.fna.gz", "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGAAATTT\n")

    class _PlantaeDbFixtureHandler(SimpleHTTPRequestHandler):
        def __init__(self, *args, directory=None, **kwargs):
            super().__init__(*args, directory=str(server_root), **kwargs)

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
                        "esearchresult": {"idlist": ["99999"]},
                    }
                )
                return
            if self.path.startswith("/eutils/esummary.fcgi"):
                self._send_json(
                    {
                        "header": {"type": "esummary", "version": "0.3"},
                        "result": {
                            "uids": ["99999"],
                            "99999": {
                                "ftppath_genbank": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/999/999/999/GCA_999999999.1_ASM999999v1",
                                "organism": "Erigeron breviscapus",
                                "speciesname": "Erigeron breviscapus",
                            },
                        },
                    }
                )
                return
            super().do_GET()

        def log_message(self, *args):
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), _PlantaeDbFixtureHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        plantae_page.write_text(
            (
                "<html><body>"
                '<a href="http://127.0.0.1:{port}/data-hub/genome/GCA_999999999.1/">NCBI genome</a>'
                "</body></html>"
            ).format(port=server.server_port),
            encoding="utf-8",
        )
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "plantaedb",
                    "id": "http://127.0.0.1:{}/taxa/species/erigeron-breviscapus".format(server.server_port),
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

        env = dict(os.environ)
        env["GG_NCBI_EUTILS_BASE_URL"] = "http://127.0.0.1:{}/eutils".format(server.server_port)
        env["GG_NCBI_FTP_BASE_URL"] = "http://127.0.0.1:{}".format(server.server_port)
        env["GG_PLANTAEDB_WEB_BASE_URL"] = "http://127.0.0.1:{}".format(server.server_port)

        download_dir = tmp_path / "download_cache"
        completed = run_script(
            "--provider",
            "plantaedb",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(download_dir),
            "--download-only",
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    finally:
        server.shutdown()
        thread.join(timeout=5)

    raw_dir = download_dir / "PlantaeDB" / "species_wise_original" / "Erigeron_breviscapus"
    assert (raw_dir / "GCA_999999999.1_ASM999999v1_cds_from_genomic.fna.gz").exists()
    assert (raw_dir / "GCA_999999999.1_ASM999999v1_genomic.gff.gz").exists()
    assert (raw_dir / "GCA_999999999.1_ASM999999v1_genomic.fna.gz").exists()


def test_download_manifest_plantaedb_formats_after_download(tmp_path):
    server_root = tmp_path / "server_root"
    plantae_page = server_root / "taxa" / "species" / "berberis-thunbergii"
    plantae_page.parent.mkdir(parents=True, exist_ok=True)

    ftp_dir = (
        server_root
        / "genomes"
        / "all"
        / "GCA"
        / "999"
        / "999"
        / "998"
        / "GCA_999999998.1_ASM999998v1"
    )
    ftp_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(ftp_dir / "GCA_999999998.1_ASM999998v1_genomic.gbff.gz", "wt", encoding="utf-8") as handle:
        handle.write(
            """LOCUS       chr1                      9 bp    DNA     linear   PLN 01-JAN-2000\n"""
            """DEFINITION  synthetic test record.\n"""
            """ACCESSION   chr1\n"""
            """VERSION     chr1\n"""
            """FEATURES             Location/Qualifiers\n"""
            """     source          1..9\n"""
            """                     /organism=\"Berberis thunbergii\"\n"""
            """     gene            1..9\n"""
            """                     /gene=\"gene1\"\n"""
            """                     /locus_tag=\"LOC1\"\n"""
            """     mRNA            1..9\n"""
            """                     /gene=\"gene1\"\n"""
            """                     /locus_tag=\"LOC1\"\n"""
            """                     /transcript_id=\"rna-LOC1\"\n"""
            """     CDS             1..9\n"""
            """                     /gene=\"gene1\"\n"""
            """                     /locus_tag=\"LOC1\"\n"""
            """                     /protein_id=\"prot-LOC1\"\n"""
            """                     /codon_start=1\n"""
            """                     /translation=\"MKF\"\n"""
            """ORIGIN\n"""
            """        1 atgaaattt\n"""
            """//\n"""
        )
    with gzip.open(ftp_dir / "GCA_999999998.1_ASM999998v1_genomic.fna.gz", "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGAAATTT\n")

    class _PlantaeDbFormatFixtureHandler(SimpleHTTPRequestHandler):
        def __init__(self, *args, directory=None, **kwargs):
            super().__init__(*args, directory=str(server_root), **kwargs)

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
                        "esearchresult": {"idlist": ["99998"]},
                    }
                )
                return
            if self.path.startswith("/eutils/esummary.fcgi"):
                self._send_json(
                    {
                        "header": {"type": "esummary", "version": "0.3"},
                        "result": {
                            "uids": ["99998"],
                            "99998": {
                                "ftppath_genbank": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/999/999/998/GCA_999999998.1_ASM999998v1",
                                "organism": "Berberis thunbergii",
                                "speciesname": "Berberis thunbergii",
                            },
                        },
                    }
                )
                return
            super().do_GET()

        def log_message(self, *args):
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), _PlantaeDbFormatFixtureHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        plantae_page.write_text(
            (
                "<html><body>"
                '<a href="http://127.0.0.1:{port}/data-hub/genome/GCA_999999998.1/">NCBI genome</a>'
                "</body></html>"
            ).format(port=server.server_port),
            encoding="utf-8",
        )
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "plantaedb",
                    "id": "http://127.0.0.1:{}/taxa/species/berberis-thunbergii".format(server.server_port),
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

        env = dict(os.environ)
        env["GG_NCBI_EUTILS_BASE_URL"] = "http://127.0.0.1:{}/eutils".format(server.server_port)
        env["GG_NCBI_FTP_BASE_URL"] = "http://127.0.0.1:{}".format(server.server_port)
        env["GG_PLANTAEDB_WEB_BASE_URL"] = "http://127.0.0.1:{}".format(server.server_port)

        download_dir = tmp_path / "download_cache"
        cds_dir = tmp_path / "out_cds"
        gff_dir = tmp_path / "out_gff"
        genome_dir = tmp_path / "out_genome"
        summary_path = tmp_path / "species_summary.tsv"
        stats_path = tmp_path / "stats.json"
        completed = run_script(
            "--provider",
            "plantaedb",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(download_dir),
            "--species-cds-dir",
            str(cds_dir),
            "--species-gff-dir",
            str(gff_dir),
            "--species-genome-dir",
            str(genome_dir),
            "--species-summary-output",
            str(summary_path),
            "--stats-output",
            str(stats_path),
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    finally:
        server.shutdown()
        thread.join(timeout=5)

    assert sorted(path.name for path in cds_dir.glob("*.fa.gz")) == [
        "Berberis_thunbergii_GCA_999999998.1_ASM999998v1_genomic.derived.cds.fa.gz"
    ]
    assert sorted(path.name for path in gff_dir.glob("*.gz")) == [
        "Berberis_thunbergii_GCA_999999998.1_ASM999998v1_genomic.derived.gff.gz"
    ]
    assert sorted(path.name for path in genome_dir.glob("*.fa.gz")) == [
        "Berberis_thunbergii_GCA_999999998.1_ASM999998v1_genomic.fa.gz"
    ]
    summary_text = summary_path.read_text(encoding="utf-8")
    assert "Berberis_thunbergii" in summary_text
    stats = json.loads(stats_path.read_text(encoding="utf-8"))
    assert stats["cds_sequences_before"] == 1
    assert stats["cds_sequences_after"] == 1
    assert stats["aggregated_cds_removed"] == 0


def test_download_manifest_resolves_citrusgenomedb_organism_page_to_public_bundle(tmp_path):
    server_root = tmp_path / "server_root"
    organism_page = server_root / "organism" / "5799"
    published_analysis = server_root / "Analysis" / "2530647"
    unpublished_analysis = server_root / "Analysis" / "9999999"
    organism_page.parent.mkdir(parents=True, exist_ok=True)
    published_analysis.parent.mkdir(parents=True, exist_ok=True)

    cds_path = server_root / "citrus_downloads" / "Citrus_australasica" / "Cau_AZM_v1.0" / "genes" / "AZM.v1.0.CDS.fa.gz"
    gff_path = server_root / "citrus_downloads" / "Citrus_australasica" / "Cau_AZM_v1.0" / "genes" / "AZM.v1.0.gene.model.gff3.gz"
    genome_path = server_root / "citrus_downloads" / "Citrus_australasica" / "Cau_AZM_v1.0" / "assembly" / "AZM.v1.0.genome.fa.gz"
    cds_path.parent.mkdir(parents=True, exist_ok=True)
    genome_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(cds_path, "wt", encoding="utf-8") as handle:
        handle.write(">AZM_gene1.t1\nATGAAATTT\n")
    with gzip.open(gff_path, "wt", encoding="utf-8") as handle:
        handle.write(
            "##gff-version 3\n"
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=AZM_gene1\n"
            "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=AZM_gene1.t1;Parent=AZM_gene1\n"
            "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds_AZM_gene1.t1;Parent=AZM_gene1.t1\n"
        )
    with gzip.open(genome_path, "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGAAATTT\n")

    organism_page.write_text(
        (
            "<html><head><title>Citrus australasica | Citrus Genome Database</title></head><body>"
            '<a href="/Analysis/9999999">draft assembly</a>'
            '<a href="/Analysis/2530647">Citrus australasica cv. AZM genome v1.0</a>'
            "</body></html>"
        ),
        encoding="utf-8",
    )
    unpublished_analysis.write_text(
        (
            "<html><head><title>Draft citrus analysis</title></head><body>"
            "Bulk download of the assembly files will become available once the data is published"
            "</body></html>"
        ),
        encoding="utf-8",
    )
    published_analysis.write_text(
        (
            "<html><head><title>Citrus australasica cv. AZM genome v1.0 | Citrus Genome Database</title></head><body>"
            '<a href="/citrus_downloads/Citrus_australasica/Cau_AZM_v1.0/assembly/AZM.v1.0.genome.fa.gz">genome</a>'
            '<a href="/citrus_downloads/Citrus_australasica/Cau_AZM_v1.0/genes/AZM.v1.0.CDS.fa.gz">CDS</a>'
            '<a href="/citrus_downloads/Citrus_australasica/Cau_AZM_v1.0/genes/AZM.v1.0.gene.model.gff3.gz">GFF3</a>'
            "</body></html>"
        ),
        encoding="utf-8",
    )

    handler = lambda *args, **kwargs: SimpleHTTPRequestHandler(*args, directory=str(server_root), **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "citrusgenomedb",
                    "id": "http://127.0.0.1:{}/organism/5799".format(server.server_port),
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

        env = dict(os.environ)
        env["GG_CITRUSGENOMEDB_WEB_BASE_URL"] = "http://127.0.0.1:{}".format(server.server_port)
        download_dir = tmp_path / "download_cache"
        completed = run_script(
            "--provider",
            "citrusgenomedb",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(download_dir),
            "--download-only",
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    finally:
        server.shutdown()
        server.server_close()

    raw_dir = download_dir / "CitrusGenomeDB" / "species_wise_original" / "Citrus_australasica"
    assert (raw_dir / "AZM.v1.0.CDS.fa.gz").exists()
    assert (raw_dir / "AZM.v1.0.gene.model.gff3.gz").exists()
    assert (raw_dir / "AZM.v1.0.genome.fa.gz").exists()


def test_download_manifest_citrusgenomedb_analysis_page_derives_cds_from_gff_and_genome(tmp_path):
    server_root = tmp_path / "server_root"
    analysis_page = server_root / "Analysis" / "3267147"
    analysis_page.parent.mkdir(parents=True, exist_ok=True)

    gff_path = server_root / "citrus_downloads" / "Citrus_limon" / "Cl_Xiangshui_v1" / "genes" / "GWHCBFU00000000.cgd_gene.gff.gz"
    genome_path = server_root / "citrus_downloads" / "Citrus_limon" / "Cl_Xiangshui_v1" / "assembly" / "GWHCBFU00000000.genome.fasta_NewID.fasta.gz"
    protein_path = server_root / "citrus_downloads" / "Citrus_limon" / "Cl_Xiangshui_v1" / "genes" / "GWHCBFU00000000.Protein_editIDFinial.faa.gz"
    rna_path = server_root / "citrus_downloads" / "Citrus_limon" / "Cl_Xiangshui_v1" / "genes" / "GWHCBFU00000000.RNA_editIDFinial.faa.gz"
    gff_path.parent.mkdir(parents=True, exist_ok=True)
    genome_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(gff_path, "wt", encoding="utf-8") as handle:
        handle.write(
            "#OriSeqID=CTG_101\tAccession=GWHCBFU00000054\n"
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=limon_gene1\n"
            "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=limon_gene1.t1;Parent=limon_gene1\n"
            "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=cds_limon_gene1.t1a;Parent=limon_gene1.t1\n"
            "chr1\tsrc\tCDS\t7\t9\t.\t+\t0\tID=cds_limon_gene1.t1b;Parent=limon_gene1.t1\n"
        )
    with gzip.open(genome_path, "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGAAATTT\n")
    with gzip.open(protein_path, "wt", encoding="utf-8") as handle:
        handle.write(">protein1\nMKF\n")
    with gzip.open(rna_path, "wt", encoding="utf-8") as handle:
        handle.write(">rna1\nMKF\n")

    analysis_page.write_text(
        (
            "<html><head><title>Citrus limon genome v1.0 | Citrus Genome Database</title></head><body>"
            '<a href="/citrus_downloads/Citrus_limon/Cl_Xiangshui_v1/assembly/GWHCBFU00000000.genome.fasta_NewID.fasta.gz">genome</a>'
            '<a href="/citrus_downloads/Citrus_limon/Cl_Xiangshui_v1/genes/GWHCBFU00000000.cgd_gene.gff.gz">GFF3</a>'
            '<a href="/citrus_downloads/Citrus_limon/Cl_Xiangshui_v1/genes/GWHCBFU00000000.Protein_editIDFinial.faa.gz">protein</a>'
            '<a href="/citrus_downloads/Citrus_limon/Cl_Xiangshui_v1/genes/GWHCBFU00000000.RNA_editIDFinial.faa.gz">RNA</a>'
            "</body></html>"
        ),
        encoding="utf-8",
    )

    handler = lambda *args, **kwargs: SimpleHTTPRequestHandler(*args, directory=str(server_root), **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "citrusgenomedb",
                    "id": "http://127.0.0.1:{}/Analysis/3267147".format(server.server_port),
                    "species_key": "Citrus_x_limon",
                    "cds_url": "",
                    "gff_url": "",
                    "genome_url": "",
                    "cds_filename": "",
                    "gff_filename": "",
                    "genome_filename": "",
                }
            ],
        )

        env = dict(os.environ)
        env["GG_CITRUSGENOMEDB_WEB_BASE_URL"] = "http://127.0.0.1:{}".format(server.server_port)
        download_dir = tmp_path / "download_cache"
        out_cds = tmp_path / "out_cds"
        out_gff = tmp_path / "out_gff"
        out_genome = tmp_path / "out_genome"
        completed = run_script(
            "--provider",
            "citrusgenomedb",
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
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    finally:
        server.shutdown()
        server.server_close()

    raw_dir = download_dir / "CitrusGenomeDB" / "species_wise_original" / "Citrus_x_limon"
    assert not any("cds" in path.name.lower() for path in raw_dir.iterdir())
    assert (raw_dir / "GWHCBFU00000000.cgd_gene.gff.gz").exists()
    assert (raw_dir / "GWHCBFU00000000.genome.fasta_NewID.fasta.gz").exists()

    formatted_cds_candidates = sorted(out_cds.glob("Citrus_x_limon*.fa.gz"))
    assert len(formatted_cds_candidates) == 1
    with gzip.open(formatted_cds_candidates[0], "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert ">Citrus_x_limon_limon_gene1" in cds_text
    assert "ATGTTT" in cds_text


def test_download_manifest_citrusgenomedb_analysis_page_prefers_cds_over_cdna(tmp_path):
    server_root = tmp_path / "server_root"
    analysis_page = server_root / "Analysis" / "6981406"
    analysis_page.parent.mkdir(parents=True, exist_ok=True)

    cds_path = server_root / "citrus_downloads" / "Citrus_aurantium" / "Ca_ZGSC_v1.0" / "genes" / "CGD_ZGSC-M.CDS.fa.gz"
    cdna_path = server_root / "citrus_downloads" / "Citrus_aurantium" / "Ca_ZGSC_v1.0" / "genes" / "CGD_ZGSC-M.cDNA.fa.gz"
    gff_path = server_root / "citrus_downloads" / "Citrus_aurantium" / "Ca_ZGSC_v1.0" / "genes" / "CGD_ZGSC-M.gene.model.gff3.gz"
    genome_path = server_root / "citrus_downloads" / "Citrus_aurantium" / "Ca_ZGSC_v1.0" / "assembly" / "ZGSC-M.genome.fa.gz"
    cds_path.parent.mkdir(parents=True, exist_ok=True)
    genome_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(cds_path, "wt", encoding="utf-8") as handle:
        handle.write(">citrus_gene1.t1\nATGAAATTT\n")
    with gzip.open(cdna_path, "wt", encoding="utf-8") as handle:
        handle.write(">citrus_gene1.t1\nATGAAATTTAAA\n")
    with gzip.open(gff_path, "wt", encoding="utf-8") as handle:
        handle.write(
            "##gff-version 3\n"
            "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tID=citrus_gene1\n"
            "chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=citrus_gene1.t1;Parent=citrus_gene1\n"
            "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds_citrus_gene1.t1;Parent=citrus_gene1.t1\n"
        )
    with gzip.open(genome_path, "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGAAATTTAAA\n")

    analysis_page.write_text(
        (
            "<html><head><title>Citrus aurantium cv. ZGSC genome v1.0 | Citrus Genome Database</title></head><body>"
            '<a href="/citrus_downloads/Citrus_aurantium/Ca_ZGSC_v1.0/assembly/ZGSC-M.genome.fa.gz">genome</a>'
            '<a href="/citrus_downloads/Citrus_aurantium/Ca_ZGSC_v1.0/genes/CGD_ZGSC-M.gene.model.gff3.gz">GFF3</a>'
            '<a href="/citrus_downloads/Citrus_aurantium/Ca_ZGSC_v1.0/genes/CGD_ZGSC-M.cDNA.fa.gz">cDNA</a>'
            '<a href="/citrus_downloads/Citrus_aurantium/Ca_ZGSC_v1.0/genes/CGD_ZGSC-M.CDS.fa.gz">CDS</a>'
            "</body></html>"
        ),
        encoding="utf-8",
    )

    handler = lambda *args, **kwargs: SimpleHTTPRequestHandler(*args, directory=str(server_root), **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "citrusgenomedb",
                    "id": "http://127.0.0.1:{}/Analysis/6981406".format(server.server_port),
                    "species_key": "Citrus_x_aurantium",
                    "cds_url": "",
                    "gff_url": "",
                    "genome_url": "",
                    "cds_filename": "",
                    "gff_filename": "",
                    "genome_filename": "",
                }
            ],
        )
        env = dict(os.environ)
        env["GG_CITRUSGENOMEDB_WEB_BASE_URL"] = "http://127.0.0.1:{}".format(server.server_port)
        resolved_manifest = tmp_path / "resolved.tsv"
        completed = run_script(
            "--provider",
            "all",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(tmp_path / "download_cache"),
            "--resolved-manifest-output",
            str(resolved_manifest),
            "--download-only",
            "--dry-run",
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)

    rows = list(csv.DictReader(resolved_manifest.open(encoding="utf-8"), delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert row["cds_url"].endswith("/genes/CGD_ZGSC-M.CDS.fa.gz")


def test_download_manifest_citrusgenomedb_organism_page_prefers_assembly_bundle_over_gene_fasta(tmp_path):
    server_root = tmp_path / "server_root"
    organism_page = server_root / "organism" / "5753"
    analysis_page = server_root / "Analysis" / "3267147"
    organism_page.parent.mkdir(parents=True, exist_ok=True)
    analysis_page.parent.mkdir(parents=True, exist_ok=True)

    gff_path = server_root / "citrus_downloads" / "Citrus_limon" / "Cl_Xiangshui_v1" / "genes" / "GWHCBFU00000000.cgd_gene.gff.gz"
    genome_path = server_root / "citrus_downloads" / "Citrus_limon" / "Cl_Xiangshui_v1" / "assembly" / "GWHCBFU00000000.genome.fasta_NewID.fasta.gz"
    gff_path.parent.mkdir(parents=True, exist_ok=True)
    genome_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(gff_path, "wt", encoding="utf-8") as handle:
        handle.write(
            "##gff-version 3\n"
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=limon_gene1\n"
            "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=limon_gene1.t1;Parent=limon_gene1\n"
            "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds_limon_gene1.t1;Parent=limon_gene1.t1\n"
        )
    with gzip.open(genome_path, "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGAAATTT\n")

    organism_page.write_text(
        (
            "<html><head><title>Citrus limon | Citrus Genome Database</title></head><body>"
            '<a href="/citrus_downloads/Citrus_limon/C.limon_EMF-UC_v1-Primary_genome/annotation/Climon_v1_primary-mRNA.fa">mRNA</a>'
            '<a href="/citrus_downloads/Citrus_limon/C.limon_EMF-UC_v1-Primary_genome/annotation/Climon_v1_primary-annotation.gff">GFF</a>'
            '<a href="/citrus_downloads/Citrus_limon/C.limon_EMF-UC_v1-Primary_genome/annotation/Climon_v1_primary-genes.fa">genes</a>'
            '<a href="/Analysis/3267147">Citrus limon genome v1.0</a>'
            "</body></html>"
        ),
        encoding="utf-8",
    )
    analysis_page.write_text(
        (
            "<html><head><title>Citrus limon genome v1.0 | Citrus Genome Database</title></head><body>"
            '<a href="/citrus_downloads/Citrus_limon/Cl_Xiangshui_v1/assembly/GWHCBFU00000000.genome.fasta_NewID.fasta.gz">genome</a>'
            '<a href="/citrus_downloads/Citrus_limon/Cl_Xiangshui_v1/genes/GWHCBFU00000000.cgd_gene.gff.gz">GFF3</a>'
            "</body></html>"
        ),
        encoding="utf-8",
    )

    handler = lambda *args, **kwargs: SimpleHTTPRequestHandler(*args, directory=str(server_root), **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "citrusgenomedb",
                    "id": "http://127.0.0.1:{}/organism/5753".format(server.server_port),
                    "species_key": "Citrus_x_limon",
                    "cds_url": "",
                    "gff_url": "",
                    "genome_url": "",
                    "cds_filename": "",
                    "gff_filename": "",
                    "genome_filename": "",
                }
            ],
        )

        env = dict(os.environ)
        env["GG_CITRUSGENOMEDB_WEB_BASE_URL"] = "http://127.0.0.1:{}".format(server.server_port)
        resolved_manifest = tmp_path / "resolved.tsv"
        completed = run_script(
            "--provider",
            "all",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(tmp_path / "download_cache"),
            "--resolved-manifest-output",
            str(resolved_manifest),
            "--download-only",
            "--dry-run",
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    finally:
        server.shutdown()
        server.server_close()

    rows = list(csv.DictReader(resolved_manifest.open(encoding="utf-8"), delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert row["species_key"] == "Citrus_x_limon"
    assert row["gff_url"].endswith("/genes/GWHCBFU00000000.cgd_gene.gff.gz")
    assert row["genome_url"].endswith("/assembly/GWHCBFU00000000.genome.fasta_NewID.fasta.gz")
    assert row["cds_url"] == ""


def test_download_manifest_citrusgenomedb_analysis_page_follows_repository_index(tmp_path):
    server_root = tmp_path / "server_root"
    analysis_page = server_root / "analysis" / "189"
    repository_index = server_root / "citrus_downloads" / "Citrus_medica" / "C.medica_Hzau_v1_genome" / "index.html"
    annotation_index = server_root / "citrus_downloads" / "Citrus_medica" / "C.medica_Hzau_v1_genome" / "annotation" / "index.html"
    assembly_index = server_root / "citrus_downloads" / "Citrus_medica" / "C.medica_Hzau_v1_genome" / "assembly" / "index.html"
    analysis_page.parent.mkdir(parents=True, exist_ok=True)
    repository_index.parent.mkdir(parents=True, exist_ok=True)
    annotation_index.parent.mkdir(parents=True, exist_ok=True)
    assembly_index.parent.mkdir(parents=True, exist_ok=True)

    cds_path = server_root / "citrus_downloads" / "Citrus_medica" / "C.medica_Hzau_v1_genome" / "annotation" / "C.medica_Hzau_v1.cds.fa.gz"
    gff_path = server_root / "citrus_downloads" / "Citrus_medica" / "C.medica_Hzau_v1_genome" / "annotation" / "C.medica_Hzau_v1.gff3.gz"
    genome_path = server_root / "citrus_downloads" / "Citrus_medica" / "C.medica_Hzau_v1_genome" / "assembly" / "C.medica_Hzau_v1.genome.fa.gz"
    with gzip.open(cds_path, "wt", encoding="utf-8") as handle:
        handle.write(">medica_gene1.t1\nATGAAATTT\n")
    with gzip.open(gff_path, "wt", encoding="utf-8") as handle:
        handle.write(
            "##gff-version 3\n"
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=medica_gene1\n"
            "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=medica_gene1.t1;Parent=medica_gene1\n"
            "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds_medica_gene1.t1;Parent=medica_gene1.t1\n"
        )
    with gzip.open(genome_path, "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGAAATTT\n")

    analysis_page.write_text(
        (
            "<html><head><title>Citrus medica genome v1.0 | Citrus Genome Database</title></head><body>"
            '<a href="/citrus_downloads/Citrus_medica/C.medica_Hzau_v1_genome/">CGD data repository</a>'
            "</body></html>"
        ),
        encoding="utf-8",
    )
    repository_index.write_text(
        (
            "<html><body>"
            '<a href="/citrus_downloads/Citrus_medica/C.medica_Hzau_v1_genome/annotation/">annotation/</a>'
            '<a href="/citrus_downloads/Citrus_medica/C.medica_Hzau_v1_genome/assembly/">assembly/</a>'
            "</body></html>"
        ),
        encoding="utf-8",
    )
    annotation_index.write_text(
        (
            "<html><body>"
            '<a href="/citrus_downloads/Citrus_medica/C.medica_Hzau_v1_genome/annotation/C.medica_Hzau_v1.cds.fa.gz">CDS</a>'
            '<a href="/citrus_downloads/Citrus_medica/C.medica_Hzau_v1_genome/annotation/C.medica_Hzau_v1.gff3.gz">GFF3</a>'
            "</body></html>"
        ),
        encoding="utf-8",
    )
    assembly_index.write_text(
        (
            "<html><body>"
            '<a href="/citrus_downloads/Citrus_medica/C.medica_Hzau_v1_genome/assembly/C.medica_Hzau_v1.genome.fa.gz">genome</a>'
            "</body></html>"
        ),
        encoding="utf-8",
    )

    handler = lambda *args, **kwargs: SimpleHTTPRequestHandler(*args, directory=str(server_root), **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "citrusgenomedb",
                    "id": "http://127.0.0.1:{}/analysis/189".format(server.server_port),
                    "species_key": "Citrus_medica",
                    "cds_url": "",
                    "gff_url": "",
                    "genome_url": "",
                    "cds_filename": "",
                    "gff_filename": "",
                    "genome_filename": "",
                }
            ],
        )

        env = dict(os.environ)
        env["GG_CITRUSGENOMEDB_WEB_BASE_URL"] = "http://127.0.0.1:{}".format(server.server_port)
        resolved_manifest = tmp_path / "resolved.tsv"
        completed = run_script(
            "--provider",
            "all",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(tmp_path / "download_cache"),
            "--resolved-manifest-output",
            str(resolved_manifest),
            "--download-only",
            "--dry-run",
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    finally:
        server.shutdown()
        server.server_close()

    rows = list(csv.DictReader(resolved_manifest.open(encoding="utf-8"), delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert row["species_key"] == "Citrus_medica"
    assert row["cds_url"].endswith("/annotation/C.medica_Hzau_v1.cds.fa.gz")
    assert row["gff_url"].endswith("/annotation/C.medica_Hzau_v1.gff3.gz")
    assert row["genome_url"].endswith("/assembly/C.medica_Hzau_v1.genome.fa.gz")


def test_download_manifest_figshare_article_resolves_requested_filenames(tmp_path):
    server_root = tmp_path / "server_root"
    api_dir = server_root / "v2" / "articles"
    api_dir.mkdir(parents=True, exist_ok=True)

    article_payload = {
        "id": 28759280,
        "title": "Construction of the super pan-genome for the genus <i>Actinidia</i>",
        "files": [
            {
                "id": 53524346,
                "name": "MW_GeneModels.gff3",
                "download_url": "http://127.0.0.1:9/files/53524346",
            },
            {
                "id": 53524460,
                "name": "MW_chr.fasta",
                "download_url": "http://127.0.0.1:9/files/53524460",
            },
            {
                "id": 58581835,
                "name": "MW.CDS.fasta",
                "download_url": "http://127.0.0.1:9/files/58581835",
            },
            {
                "id": 58581799,
                "name": "KY.CDS.fasta",
                "download_url": "http://127.0.0.1:9/files/58581799",
            },
        ],
    }
    (api_dir / "28759280").write_text(json.dumps(article_payload), encoding="utf-8")

    handler = lambda *args, **kwargs: SimpleHTTPRequestHandler(*args, directory=str(server_root), **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "figshare",
                    "id": "https://figshare.com/articles/dataset/example_bundle/28759280",
                    "species_key": "Actinidia_deliciosa",
                    "cds_filename": "MW.CDS.fasta",
                    "gff_filename": "MW_GeneModels.gff3",
                    "genome_filename": "MW_chr.fasta",
                }
            ],
        )

        env = dict(os.environ)
        env["GG_FIGSHARE_API_BASE_URL"] = "http://127.0.0.1:{}/v2".format(server.server_port)
        resolved_manifest = tmp_path / "resolved.tsv"
        completed = run_script(
            "--provider",
            "all",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(tmp_path / "download_cache"),
            "--resolved-manifest-output",
            str(resolved_manifest),
            "--download-only",
            "--dry-run",
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    finally:
        server.shutdown()
        server.server_close()

    rows = list(csv.DictReader(resolved_manifest.open(encoding="utf-8"), delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert row["species_key"] == "Actinidia_deliciosa"
    assert row["cds_url"] == "http://127.0.0.1:9/files/58581835"
    assert row["gff_url"] == "http://127.0.0.1:9/files/53524346"
    assert row["genome_url"] == "http://127.0.0.1:9/files/53524460"
    assert row["cds_filename"] == "MW.CDS.fasta"
    assert row["gff_filename"] == "MW_GeneModels.gff3"
    assert row["genome_filename"] == "MW_chr.fasta"


def test_download_manifest_figshare_article_supports_archive_members(tmp_path):
    server_root = tmp_path / "server_root"
    api_dir = server_root / "v2" / "articles"
    files_dir = server_root / "files"
    api_dir.mkdir(parents=True, exist_ok=True)
    files_dir.mkdir(parents=True, exist_ok=True)

    bundle_path = files_dir / "YouCha.annotation.tar.gz"
    with tarfile.open(bundle_path, "w:gz") as archive:
        for name, content in {
            "YouCha.annotation/A/camellia_meiocarpa.gene.cds.fa": ">gene1\nATGAAATTT\n",
            "YouCha.annotation/A/camellia_meiocarpa.gene.gff": "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1\n",
        }.items():
            payload = content.encode("utf-8")
            info = tarfile.TarInfo(name=name)
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))

    genome_path = files_dir / "youcha.Changed.A.fasta.gz"
    with gzip.open(genome_path, "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGAAATTT\n")

    article_payload = {
        "id": 26926918,
        "title": "Haplotype-resolved genome assembly of the tetraploid Youcha tree Camellia meiocarpa Hu",
        "files": [
            {
                "id": 48978784,
                "name": "YouCha.annotation.tar.gz",
                "download_url": "http://127.0.0.1:0/files/YouCha.annotation.tar.gz",
            },
            {
                "id": 48976912,
                "name": "youcha.Changed.A.fasta.gz",
                "download_url": "http://127.0.0.1:0/files/youcha.Changed.A.fasta.gz",
            },
        ],
    }

    class _FigshareFixtureHandler(SimpleHTTPRequestHandler):
        def __init__(self, *args, root_dir=None, **kwargs):
            self._root_dir = root_dir
            super().__init__(*args, directory=str(root_dir), **kwargs)

        def do_GET(self):
            if self.path == "/v2/articles/26926918":
                payload = json.dumps(article_payload).replace(":0/", ":{}".format(self.server.server_port) + "/").encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(payload)))
                self.end_headers()
                self.wfile.write(payload)
                return
            super().do_GET()

    handler = lambda *args, **kwargs: _FigshareFixtureHandler(*args, root_dir=server_root, **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "figshare",
                    "id": "https://figshare.com/articles/dataset/camellia_meiocarpa_bundle/26926918",
                    "species_key": "Camellia_meiocarpa",
                    "cds_filename": "YouCha.annotation.tar.gz",
                    "gff_filename": "YouCha.annotation.tar.gz",
                    "genome_filename": "youcha.Changed.A.fasta.gz",
                    "cds_archive_member": "YouCha.annotation/A/camellia_meiocarpa.gene.cds.fa",
                    "gff_archive_member": "YouCha.annotation/A/camellia_meiocarpa.gene.gff",
                }
            ],
        )

        env = dict(os.environ)
        env["GG_FIGSHARE_API_BASE_URL"] = "http://127.0.0.1:{}/v2".format(server.server_port)
        resolved_manifest = tmp_path / "resolved.tsv"
        completed = run_script(
            "--provider",
            "all",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(tmp_path / "download_cache"),
            "--resolved-manifest-output",
            str(resolved_manifest),
            "--download-only",
            "--dry-run",
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    finally:
        server.shutdown()
        server.server_close()

    rows = list(csv.DictReader(resolved_manifest.open(encoding="utf-8"), delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert row["species_key"] == "Camellia_meiocarpa"
    assert row["cds_url"].endswith("/files/YouCha.annotation.tar.gz")
    assert row["gff_url"].endswith("/files/YouCha.annotation.tar.gz")
    assert row["genome_url"].endswith("/files/youcha.Changed.A.fasta.gz")
    assert row["cds_archive_member"] == "YouCha.annotation/A/camellia_meiocarpa.gene.cds.fa"
    assert row["gff_archive_member"] == "YouCha.annotation/A/camellia_meiocarpa.gene.gff"
    assert row["cds_filename"] == "camellia_meiocarpa.gene.cds.fa"
    assert row["gff_filename"] == "camellia_meiocarpa.gene.gff"
    assert row["genome_filename"] == "youcha.Changed.A.fasta.gz"


def test_download_manifest_resolves_plantgarden_assembly_page_to_public_bundle(tmp_path):
    server_root = tmp_path / "server_root"
    assembly_page = server_root / "en" / "list" / "t64480" / "genome" / "t64480.G001"
    download_dir = server_root / "en" / "download" / "t64480" / "t64480.G001"
    assembly_page.parent.mkdir(parents=True, exist_ok=True)
    download_dir.mkdir(parents=True, exist_ok=True)

    with gzip.open(download_dir / "APO1.1.cds.fasta.gz", "wt", encoding="utf-8") as handle:
        handle.write(">Apo1.1ch29g20170.1\nATGAAATTT\n")
    with gzip.open(download_dir / "APO1.1.genes.gff.gz", "wt", encoding="utf-8") as handle:
        handle.write(
            "##gff-version 3\n"
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=Apo1.1ch29g20170\n"
            "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=Apo1.1ch29g20170.1;Parent=Apo1.1ch29g20170\n"
            "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds_Apo1.1ch29g20170.1;Parent=Apo1.1ch29g20170.1\n"
        )
    with gzip.open(download_dir / "APO_r1.1.pmol.fasta.gz", "wt", encoding="utf-8") as handle:
        handle.write(">APO1.1ch01\nATGAAATTT\n")

    assembly_page.write_text(
        (
            "<html><body>"
            "<div id='k_contents' "
            'data-pgtag-genome_assembly_id="t64480.G001" '
            'data-pgtag-species_id="t64480" '
            'data-pgtag-species_name="Actinidia polygama" '
            'data-pgtag-assembly_version="APO_r1.1"></div>'
            '<div class="_download" onClick="k_vexDownloadModal(\'en\',\'t64480.G001\',\'gid\')"></div>'
            "</body></html>"
        ),
        encoding="utf-8",
    )

    handler = lambda *args, **kwargs: SimpleHTTPRequestHandler(*args, directory=str(server_root), **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "plantgarden",
                    "id": "http://127.0.0.1:{}/en/list/t64480/genome/t64480.G001".format(server.server_port),
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

        env = dict(os.environ)
        env["GG_PLANTGARDEN_WEB_BASE_URL"] = "http://127.0.0.1:{}".format(server.server_port)
        resolved_manifest = tmp_path / "resolved.tsv"
        completed = run_script(
            "--provider",
            "all",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(tmp_path / "download_cache"),
            "--resolved-manifest-output",
            str(resolved_manifest),
            "--download-only",
            "--dry-run",
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    finally:
        server.shutdown()
        server.server_close()

    rows = list(csv.DictReader(resolved_manifest.open(encoding="utf-8"), delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert row["species_key"] == "Actinidia_polygama"
    assert row["cds_url"].endswith("/en/download/t64480/t64480.G001/APO1.1.cds.fasta.gz")
    assert row["gff_url"].endswith("/en/download/t64480/t64480.G001/APO1.1.genes.gff.gz")
    assert row["genome_url"].endswith("/en/download/t64480/t64480.G001/APO_r1.1.pmol.fasta.gz")


def test_download_manifest_plantgarden_genome_list_page_falls_back_to_transcripts_when_cds_is_absent(tmp_path):
    server_root = tmp_path / "server_root"
    genome_list_dir = server_root / "en" / "list" / "t385388" / "genome"
    genome_list_page = genome_list_dir / "index.html"
    assembly_page = genome_list_dir / "t385388.G001"
    download_dir = server_root / "en" / "download" / "t385388" / "t385388.G001"
    genome_list_dir.mkdir(parents=True, exist_ok=True)
    assembly_page.parent.mkdir(parents=True, exist_ok=True)
    download_dir.mkdir(parents=True, exist_ok=True)

    with gzip.open(download_dir / "CON_genome_assembly_v1.0_final.fa.gz", "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGAAATTT\n")
    with gzip.open(download_dir / "Chrall.genes.gff.gz", "wt", encoding="utf-8") as handle:
        handle.write(
            "##gff-version 3\n"
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=CON_gene1\n"
            "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=CON_gene1.t1;Parent=CON_gene1\n"
            "chr1\tsrc\texon\t1\t9\t.\t+\t.\tID=exon_CON_gene1.t1;Parent=CON_gene1.t1\n"
        )
    with gzip.open(download_dir / "Chrall.transcripts.fasta.gz", "wt", encoding="utf-8") as handle:
        handle.write(">CON_gene1.t1\nATGAAATTT\n")
    with gzip.open(download_dir / "Chrall.proteins.fasta.gz", "wt", encoding="utf-8") as handle:
        handle.write(">CON_gene1.p1\nMKF\n")

    genome_list_page.write_text(
        (
            "<html><body>"
            '<a href="/en/list/t385388/genome/t385388.G001">CON_genome_v1.0</a>'
            "</body></html>"
        ),
        encoding="utf-8",
    )
    assembly_page.write_text(
        (
            "<html><body>"
            "<div id='k_contents' "
            'data-pgtag-genome_assembly_id="t385388.G001" '
            'data-pgtag-species_id="t385388" '
            'data-pgtag-species_name="Camellia oleifera" '
            'data-pgtag-sub_rank="var." '
            'data-pgtag-sub_name="Nanyongensis" '
            'data-pgtag-assembly_version="CON_genome_v1.0"></div>'
            "</body></html>"
        ),
        encoding="utf-8",
    )

    handler = lambda *args, **kwargs: SimpleHTTPRequestHandler(*args, directory=str(server_root), **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "plantgarden",
                    "id": "http://127.0.0.1:{}/en/list/t385388/genome".format(server.server_port),
                    "species_key": "Camellia_oleifera",
                    "cds_url": "",
                    "gff_url": "",
                    "genome_url": "",
                    "cds_filename": "",
                    "gff_filename": "",
                    "genome_filename": "",
                }
            ],
        )

        env = dict(os.environ)
        env["GG_PLANTGARDEN_WEB_BASE_URL"] = "http://127.0.0.1:{}".format(server.server_port)
        resolved_manifest = tmp_path / "resolved.tsv"
        completed = run_script(
            "--provider",
            "all",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(tmp_path / "download_cache"),
            "--resolved-manifest-output",
            str(resolved_manifest),
            "--download-only",
            "--dry-run",
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    finally:
        server.shutdown()
        server.server_close()

    rows = list(csv.DictReader(resolved_manifest.open(encoding="utf-8"), delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert row["species_key"] == "Camellia_oleifera"
    assert row["cds_url"].endswith("/en/download/t385388/t385388.G001/Chrall.transcripts.fasta.gz")
    assert row["gff_url"].endswith("/en/download/t385388/t385388.G001/Chrall.genes.gff.gz")
    assert row["genome_url"].endswith("/en/download/t385388/t385388.G001/CON_genome_assembly_v1.0_final.fa.gz")


def test_download_manifest_jgi_credentials_enable_protected_direct_download(tmp_path):
    server_root = tmp_path / "server_root"
    protected_dir = server_root / "protected"
    protected_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(protected_dir / "test.cds.fa.gz", "wt", encoding="utf-8") as handle:
        handle.write(">gene1\nATGAAATTT\n")

    expected_login = "user@example.org"
    expected_password = "secret-password"
    expected_token = "csrf-token-123"

    class _JgiFixtureHandler(SimpleHTTPRequestHandler):
        def __init__(self, *args, directory=None, **kwargs):
            super().__init__(*args, directory=str(server_root), **kwargs)

        def do_GET(self):
            if self.path == "/signon":
                body = (
                    '<html><body><form action="/signon/create" method="post">'
                    '<input type="hidden" name="utf8" value="&#x2713;" />'
                    f'<input type="hidden" name="authenticity_token" value="{expected_token}" />'
                    '<input type="text" name="login" />'
                    '<input type="password" name="password" />'
                    '<input type="submit" name="commit" value="Sign In" />'
                    "</form></body></html>"
                ).encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Content-Length", str(len(body)))
                self.send_header("Set-Cookie", "anon=1; Path=/")
                self.end_headers()
                self.wfile.write(body)
                return
            if self.path == "/welcome":
                body = b"welcome"
                self.send_response(200)
                self.send_header("Content-Type", "text/plain")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return
            if self.path == "/protected/test.cds.fa.gz":
                cookie = self.headers.get("Cookie", "")
                if "jgi_session=ok" not in cookie:
                    self.send_error(403, "Forbidden")
                    return
            super().do_GET()

        def do_POST(self):
            if self.path != "/signon/create":
                self.send_error(404)
                return
            length = int(self.headers.get("Content-Length", "0"))
            payload = self.rfile.read(length).decode("utf-8")
            parsed = parse_qs(payload)
            if (
                parsed.get("login", [""])[0] != expected_login
                or parsed.get("password", [""])[0] != expected_password
                or parsed.get("authenticity_token", [""])[0] != expected_token
            ):
                self.send_error(403, "Invalid credentials")
                return
            self.send_response(302)
            self.send_header("Location", "/welcome")
            self.send_header("Set-Cookie", "jgi_session=ok; Path=/")
            self.end_headers()

        def log_message(self, *args):
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), _JgiFixtureHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "direct",
                    "id": "jgi_protected_test",
                    "species_key": "Jgi_cookie_test",
                    "cds_url": "http://127.0.0.1:{}/protected/test.cds.fa.gz".format(server.server_port),
                    "gff_url": "",
                    "genome_url": "",
                    "cds_filename": "test.cds.fa.gz",
                    "gff_filename": "",
                    "genome_filename": "",
                }
            ],
        )
        download_dir = tmp_path / "download_cache"
        env = dict(os.environ)
        env["GG_JGI_SIGNON_BASE_URL"] = "http://127.0.0.1:{}".format(server.server_port)
        env["GG_TEST_JGI_LOGIN"] = expected_login
        env["GG_TEST_JGI_PASSWORD"] = expected_password
        completed = run_script(
            "--provider",
            "direct",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(download_dir),
            "--download-only",
            "--jgi-login-env",
            "GG_TEST_JGI_LOGIN",
            "--jgi-password-env",
            "GG_TEST_JGI_PASSWORD",
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    finally:
        server.shutdown()
        thread.join(timeout=5)

    raw_dir = download_dir / "Direct" / "species_wise_original" / "Jgi_cookie_test"
    assert (raw_dir / "test.cds.fa.gz").exists()


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


def test_download_manifest_resolves_gwh_id_via_public_index(tmp_path):
    server_root = tmp_path / "server_root"
    folder = server_root / "gwh" / "Plants" / "Medicago_sativa_Zhongmu3_GWHIGRM00000000.1"
    folder.mkdir(parents=True, exist_ok=True)

    with gzip.open(folder / "GWHIGRM00000000.1.CDS.fasta.gz", "wt", encoding="utf-8") as handle:
        handle.write(
            ">GWHTIGRM000001.1 Protein=GWHPIGRM000001.1 Gene=GWHGIGRM000001.1 OriGeneID=MsaZM3G010000001\nATGAAATTT\n"
        )
    with gzip.open(folder / "GWHIGRM00000000.1.gff.gz", "wt", encoding="utf-8") as handle:
        handle.write("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1\n")
    with gzip.open(folder / "GWHIGRM00000000.1.genome.fasta.gz", "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGCATGC\n")

    handler = lambda *args, **kwargs: SimpleHTTPRequestHandler(*args, directory=str(server_root), **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "gwh",
                    "id": "GWHIGRM00000000.1",
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
        out_cds = tmp_path / "out_cds"
        env = dict(os.environ)
        env["GG_GWH_DOWNLOAD_BASE_URL"] = "http://127.0.0.1:{}/gwh".format(server.server_port)
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--provider",
                "gwh",
                "--download-manifest",
                str(manifest),
                "--download-dir",
                str(download_dir),
                "--species-cds-dir",
                str(out_cds),
                "--species-gff-dir",
                str(tmp_path / "out_gff"),
                "--species-genome-dir",
                str(tmp_path / "out_genome"),
            ],
            capture_output=True,
            text=True,
            check=False,
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

        raw_dir = download_dir / "GWH" / "species_wise_original" / "Medicago_sativa"
        assert (raw_dir / "GWHIGRM00000000.1.CDS.fasta.gz").exists()
        assert (raw_dir / "GWHIGRM00000000.1.gff.gz").exists()
        assert (raw_dir / "GWHIGRM00000000.1.genome.fasta.gz").exists()

        formatted_cds = out_cds / "Medicago_sativa_GWHIGRM00000000.1.CDS.fa.gz"
        with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
            cds_text = handle.read()
        assert ">Medicago_sativa_GWHGIGRM000001.1" in cds_text
    finally:
        server.shutdown()
        server.server_close()


def test_download_manifest_resolves_gwh_id_via_show_page_fallback(tmp_path):
    server_root = tmp_path / "server_root"
    download_folder = server_root / "downloads" / "Aegle_marmelos_GWHBKHK01000000"
    download_folder.mkdir(parents=True, exist_ok=True)

    with gzip.open(download_folder / "GWHBKHK01000000.CDS.fasta.gz", "wt", encoding="utf-8") as handle:
        handle.write(">gene1\nATGAAATTT\n")
    with gzip.open(download_folder / "GWHBKHK01000000.gff.gz", "wt", encoding="utf-8") as handle:
        handle.write("chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1\n")
    with gzip.open(download_folder / "GWHBKHK01000000.genome.fasta.gz", "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGCATGC\n")

    class _GwhShowHandler(SimpleHTTPRequestHandler):
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
            if parsed.path == "/gwh/gwhSearch/api":
                query = parse_qs(parsed.query)
                term = (query.get("term") or [""])[0]
                payload = {"totalHits": 0, "data": []}
                if term == "GWHBKHK01000000":
                    payload = {
                        "totalHits": 1,
                        "data": [
                            {
                                "title": "Aegle marmelos",
                                "url": "http://127.0.0.1:{}/gwh/Assembly/26289/show".format(self.server.server_port),
                                "attrs": {"source": "Direct submission", "has_annotation": "Yes"},
                            }
                        ],
                    }
                self._send_json(payload)
                return
            if parsed.path == "/gwh/Assembly/26289/show":
                body = (
                    "<html><body>"
                    "<div>Aegle marmelos</div>"
                    "<div>Scientific Name Aegle marmelos Common Names bael</div>"
                    '<a href="http://127.0.0.1:{}/downloads/Aegle_marmelos_GWHBKHK01000000/'
                    'GWHBKHK01000000.genome.fasta.gz">DNA</a>'
                    '<a href="http://127.0.0.1:{}/downloads/Aegle_marmelos_GWHBKHK01000000/'
                    'GWHBKHK01000000.gff.gz">GFF</a>'
                    '<a href="http://127.0.0.1:{}/downloads/Aegle_marmelos_GWHBKHK01000000/'
                    'GWHBKHK01000000.CDS.fasta.gz">CDS</a>'
                    "</body></html>"
                ).format(self.server.server_port, self.server.server_port, self.server.server_port)
                encoded = body.encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Content-Length", str(len(encoded)))
                self.end_headers()
                self.wfile.write(encoded)
                return
            super().do_GET()

    handler = lambda *args, **kwargs: _GwhShowHandler(*args, root_dir=server_root, **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "gwh",
                    "id": "GWHBKHK01000000",
                    "species_key": "Aegle_marmelos",
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
        out_cds = tmp_path / "out_cds"
        env = dict(os.environ)
        env["GG_GWH_DOWNLOAD_BASE_URL"] = "http://127.0.0.1:{}/missing_index".format(server.server_port)
        env["GG_GWH_WEB_BASE_URL"] = "http://127.0.0.1:{}/gwh".format(server.server_port)
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--provider",
                "gwh",
                "--download-manifest",
                str(manifest),
                "--download-dir",
                str(download_dir),
                "--species-cds-dir",
                str(out_cds),
                "--species-gff-dir",
                str(tmp_path / "out_gff"),
                "--species-genome-dir",
                str(tmp_path / "out_genome"),
            ],
            capture_output=True,
            text=True,
            check=False,
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

        raw_dir = download_dir / "GWH" / "species_wise_original" / "Aegle_marmelos"
        assert (raw_dir / "GWHBKHK01000000.CDS.fasta.gz").exists()
        assert (raw_dir / "GWHBKHK01000000.gff.gz").exists()
        assert (raw_dir / "GWHBKHK01000000.genome.fasta.gz").exists()
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def test_download_manifest_resolves_ddbj_bioproject_to_public_wgs_gbff(tmp_path):
    server_root = tmp_path / "server_root"
    gbff_dir = server_root / "public" / "ddbj_database" / "wgs" / "BA" / "AH"
    gbff_dir.mkdir(parents=True, exist_ok=True)

    gbff_text = (
        "LOCUS       BAAHMP010000001         9 bp    DNA     linear   PLN 08-JUL-2025\n"
        "DEFINITION  Triphyophyllum peltatum test sequence.\n"
        "ACCESSION   BAAHMP010000001 BAAHMP010000000\n"
        "VERSION     BAAHMP010000001.1\n"
        "DBLINK      BioProject:PRJDB15739\n"
        "SOURCE      Triphyophyllum peltatum\n"
        "  ORGANISM  Triphyophyllum peltatum\n"
        "            Eukaryota; Viridiplantae.\n"
        "FEATURES             Location/Qualifiers\n"
        "     source          1..9\n"
        "                     /mol_type=\"genomic DNA\"\n"
        "                     /organism=\"Triphyophyllum peltatum\"\n"
        "     CDS             1..9\n"
        "                     /codon_start=1\n"
        "                     /locus_tag=\"Tripe_000001\"\n"
        "                     /product=\"hypothetical protein\"\n"
        "                     /protein_id=\"GAB0000001.1\"\n"
        "                     /translation=\"MK\"\n"
        "ORIGIN\n"
        "        1 atgaaataa\n"
        "//\n"
    )
    with gzip.open(gbff_dir / "BAAHMP.gz", "wt", encoding="utf-8") as handle:
        handle.write(gbff_text)

    class _DdbjFixtureHandler(SimpleHTTPRequestHandler):
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
            if parsed.path == "/search/api/entries/bioproject/PRJDB15739":
                self._send_json(
                    {
                        "identifier": "PRJDB15739",
                        "organism": {"identifier": "63090", "name": "Triphyophyllum peltatum"},
                        "dbXrefs": [
                            {
                                "identifier": "BAAHMP000000000",
                                "type": "insdc-master",
                                "url": "https://ddbj.nig.ac.jp/search/entry/ddbj/BAAHMP000000000",
                            }
                        ],
                    }
                )
                return
            super().do_GET()

    handler = lambda *args, **kwargs: _DdbjFixtureHandler(*args, root_dir=server_root, **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "ddbj",
                    "id": "PRJDB15739",
                    "species_key": "",
                    "cds_url": "",
                    "gff_url": "",
                    "gbff_url": "",
                    "genome_url": "",
                    "cds_filename": "",
                    "gff_filename": "",
                    "gbff_filename": "",
                    "genome_filename": "",
                }
            ],
        )

        download_dir = tmp_path / "download_cache"
        out_cds = tmp_path / "out_cds"
        out_gff = tmp_path / "out_gff"
        out_genome = tmp_path / "out_genome"
        env = dict(os.environ)
        env["GG_DDBJ_SEARCH_API_BASE_URL"] = "http://127.0.0.1:{}/search/api".format(server.server_port)
        env["GG_DDBJ_PUBLIC_WGS_BASE_URL"] = "http://127.0.0.1:{}/public/ddbj_database/wgs".format(
            server.server_port
        )

        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--provider",
                "ddbj",
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

        raw_dir = download_dir / "DDBJ" / "species_wise_original" / "Triphyophyllum_peltatum"
        assert (raw_dir / "BAAHMP.gbff.gz").exists()

        cds_outputs = list(out_cds.glob("*.gz"))
        gff_outputs = list(out_gff.glob("*.gz"))
        genome_outputs = list(out_genome.glob("*.gz"))
        assert len(cds_outputs) == 1
        assert len(gff_outputs) == 1
        assert len(genome_outputs) == 1

        with gzip.open(cds_outputs[0], "rt", encoding="utf-8") as handle:
            cds_text = handle.read()
        assert "Triphyophyllum_peltatum_" in cds_text

        with gzip.open(gff_outputs[0], "rt", encoding="utf-8") as handle:
            gff_text = handle.read()
        assert "gene:Tripe_000001" in gff_text

        with gzip.open(genome_outputs[0], "rt", encoding="utf-8") as handle:
            genome_text = handle.read()
        assert ">BAAHMP010000001" in genome_text
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


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
    stale_lock.write_text(
        json.dumps(
            {
                "format": "shared-lock-v2",
                "pid": 999999999,
                "hostname": socket.gethostname(),
                "boot_id": _current_boot_id(),
                "created_at": time.time() - 3600.0,
            }
        )
        + "\n",
        encoding="utf-8",
    )

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


def test_download_manifest_does_not_reclaim_fresh_foreign_lock_file(tmp_path):
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
    foreign_lock = raw_dir / (cds_name + ".lock")
    foreign_lock.write_text(
        json.dumps(
            {
                "format": "shared-lock-v2",
                "pid": 999999999,
                "hostname": "foreign-node",
                "boot_id": "foreign-boot",
                "created_at": time.time(),
            }
        )
        + "\n",
        encoding="utf-8",
    )

    env = os.environ.copy()
    env["GG_DOWNLOAD_LOCK_ACQUIRE_TIMEOUT_SECONDS"] = "1"
    env["GG_DOWNLOAD_LOCK_POLL_SECONDS"] = "0.1"
    completed = run_script(
        "--provider",
        "coge",
        "--download-manifest",
        str(manifest),
        "--download-dir",
        str(download_dir),
        "--download-only",
        env=env,
    )
    assert completed.returncode == 1, completed.stderr + "\n" + completed.stdout
    assert not (raw_dir / cds_name).exists()
    assert "waiting for shared lock" in completed.stderr.lower()
    assert "timed out waiting for shared lock" in completed.stderr.lower()


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
        "fernbase": "Azolla_filiculoides",
        "veupathdb": "Entamoeba_nuttalli",
        "dictybase": "Dictyostelium_discoideum",
    }
    provider_ids = {
        "ensembl": "homo_sapiens",
        "ensemblplants": "Ostreococcus_lucimarinus",
        "coge": "24739",
        "cngb": "cngb:arabidopsis_thaliana_v1",
        "flybase": "dmel_r6.61",
        "wormbase": "celegans_prjna13758_ws290",
        "vectorbase": "anopheles_gambiae_pest",
        "fernbase": "Azolla_filiculoides",
        "veupathdb": "veupathdb:EnuttalliP19",
        "dictybase": "dictybase:Dictyostelium_discoideum",
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
        "fernbase": "GG_FERNBASE",
        "veupathdb": "GG_VEUPATHDB",
        "dictybase": "GG_DICTYBASE",
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
        "fernbase": download_dir / "FernBase" / "species_wise_original" / provider_species["fernbase"],
        "veupathdb": download_dir / "VEuPathDB" / "species_wise_original" / provider_species["veupathdb"],
        "dictybase": download_dir / "dictyBase" / "species_wise_original" / provider_species["dictybase"],
    }
    for provider, root in expected_roots.items():
        assert root.exists(), provider
        files = sorted([p.name for p in root.iterdir() if p.is_file()])
        assert any(name.endswith(".cds.fa") for name in files), provider
        assert any(name.endswith(".gene.gff3") for name in files), provider
        assert any(name.endswith(".genome.fa") for name in files), provider


def test_download_manifest_fernbase_provider_follows_latest_version_subdir(tmp_path):
    server_root = tmp_path / "server_root"
    v1_dir = server_root / "ftp" / "Azolla_filiculoides" / "Azolla_asm_v1.0"
    v2_dir = server_root / "ftp" / "Azolla_filiculoides" / "Azolla_asm_v1.1"
    v1_dir.mkdir(parents=True, exist_ok=True)
    v2_dir.mkdir(parents=True, exist_ok=True)

    (v1_dir / "Azolla_filiculoides.CDS.highconfidence_v1.0.fasta").write_text(">old.t1\nATG\n", encoding="utf-8")
    (v1_dir / "Azolla_filiculoides.gene_models.highconfidence_v1.0.gff").write_text(
        "chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=old\n",
        encoding="utf-8",
    )
    (v1_dir / "Azolla_filiculoides.genome_v1.0.fasta").write_text(">chr1\nATGC\n", encoding="utf-8")

    (v2_dir / "Azolla_filiculoides.CDS.lowconfidence_v1.1.fasta").write_text(">low.t1\nATG\n", encoding="utf-8")
    (v2_dir / "Azolla_filiculoides.CDS.highconfidence_v1.1.fasta").write_text(">new.t1\nATGAA\n", encoding="utf-8")
    (v2_dir / "Azolla_filiculoides.gene_models.lowconfidence_v1.1.gff").write_text(
        "chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=low\n",
        encoding="utf-8",
    )
    (v2_dir / "Azolla_filiculoides.gene_models.highconfidence_v1.1.gff").write_text(
        "chr1\tsrc\tgene\t1\t5\t.\t+\t.\tID=new\n",
        encoding="utf-8",
    )
    (v2_dir / "Azolla_filiculoides.genome_v1.2.fasta").write_text(">chr1\nATGCATGC\n", encoding="utf-8")

    handler = lambda *args, **kwargs: SimpleHTTPRequestHandler(*args, directory=str(server_root), **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "fernbase",
                    "id": "Azolla_filiculoides",
                    "species_key": "Azolla_filiculoides",
                    "cds_url": "",
                    "gff_url": "",
                    "genome_url": "",
                }
            ],
        )

        env = dict(os.environ)
        env["GG_FERNBASE_ID_URL_TEMPLATE"] = "http://127.0.0.1:{}/ftp/{{id}}/".format(server.server_port)
        download_dir = tmp_path / "download_cache"
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--provider",
                "fernbase",
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
    finally:
        server.shutdown()
        thread.join(timeout=5)

    raw_dir = download_dir / "FernBase" / "species_wise_original" / "Azolla_filiculoides"
    assert (raw_dir / "Azolla_filiculoides.CDS.highconfidence_v1.1.fasta").exists()
    assert (raw_dir / "Azolla_filiculoides.gene_models.highconfidence_v1.1.gff").exists()
    assert (raw_dir / "Azolla_filiculoides.genome_v1.2.fasta").exists()
    assert not (raw_dir / "Azolla_filiculoides.CDS.lowconfidence_v1.1.fasta").exists()


def test_download_manifest_fernbase_combined_mode_merges_non_overlapping_low_confidence_genes(tmp_path):
    server_root = tmp_path / "server_root"
    v1_dir = server_root / "ftp" / "Azolla_filiculoides" / "Azolla_asm_v1.1"
    v1_dir.mkdir(parents=True, exist_ok=True)

    (v1_dir / "Azolla_filiculoides.CDS.highconfidence_v1.1.fasta").write_text(
        ">Azfi_high1.t1 gene=Azfi_high1\nATGAAATTT\n",
        encoding="utf-8",
    )
    (v1_dir / "Azolla_filiculoides.CDS.lowconfidence_v1.1.fasta").write_text(
        ">Azfi_low_overlap.t1 gene=Azfi_low_overlap\nATGAAATTT\n>Azfi_low_keep.t1 gene=Azfi_low_keep\nATGAAATTT\n",
        encoding="utf-8",
    )
    (v1_dir / "Azolla_filiculoides.gene_models.highconfidence_v1.1.gff").write_text(
        "\n".join(
            [
                "##gff-version 3",
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=Azfi_high1",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=Azfi_high1.t1;Parent=Azfi_high1",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=Azfi_high1.t1.cds;Parent=Azfi_high1.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    (v1_dir / "Azolla_filiculoides.gene_models.lowconfidence_v1.1.gff").write_text(
        "\n".join(
            [
                "##gff-version 3",
                "chr1\tsrc\tgene\t5\t12\t.\t-\t.\tID=Azfi_low_overlap",
                "chr1\tsrc\tmRNA\t5\t12\t.\t-\t.\tID=Azfi_low_overlap.t1;Parent=Azfi_low_overlap",
                "chr1\tsrc\tCDS\t5\t12\t.\t-\t0\tID=Azfi_low_overlap.t1.cds;Parent=Azfi_low_overlap.t1",
                "chr1\tsrc\tgene\t20\t28\t.\t+\t.\tID=Azfi_low_keep",
                "chr1\tsrc\tmRNA\t20\t28\t.\t+\t.\tID=Azfi_low_keep.t1;Parent=Azfi_low_keep",
                "chr1\tsrc\tCDS\t20\t28\t.\t+\t0\tID=Azfi_low_keep.t1.cds;Parent=Azfi_low_keep.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    (v1_dir / "Azolla_filiculoides.genome_v1.2.fasta").write_text(
        ">chr1\nATGAAATTTCCCAAAGGGATGAAATTTCCCAAAGGG\n",
        encoding="utf-8",
    )

    handler = lambda *args, **kwargs: SimpleHTTPRequestHandler(*args, directory=str(server_root), **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "fernbase",
                    "id": "Azolla_filiculoides",
                    "species_key": "Azolla_filiculoides",
                    "fernbase_confidence_mode": "high-low combined",
                }
            ],
        )

        env = dict(os.environ)
        env["GG_FERNBASE_ID_URL_TEMPLATE"] = "http://127.0.0.1:{}/ftp/{{id}}/".format(server.server_port)
        download_dir = tmp_path / "download_cache"
        out_cds = tmp_path / "species_cds"
        out_gff = tmp_path / "species_gff"
        out_genome = tmp_path / "species_genome"
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--provider",
                "fernbase",
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
    finally:
        server.shutdown()
        thread.join(timeout=5)

    raw_dir = download_dir / "FernBase" / "species_wise_original" / "Azolla_filiculoides"
    combined_raw_cds = raw_dir / "Azolla_filiculoides.CDS.highlowcombined_v1.1.fasta"
    combined_raw_gff = raw_dir / "Azolla_filiculoides.gene_models.highlowcombined_v1.1.gff"
    assert (raw_dir / "Azolla_filiculoides.CDS.highconfidence_v1.1.fasta").exists()
    assert (raw_dir / "Azolla_filiculoides.CDS.lowconfidence_v1.1.fasta").exists()
    assert combined_raw_cds.exists()
    assert combined_raw_gff.exists()

    formatted_cds = out_cds / "Azolla_filiculoides_CDS.highlowcombined_v1.1.fa.gz"
    formatted_gff = out_gff / "Azolla_filiculoides_gene_models.highlowcombined_v1.1.gff.gz"
    assert formatted_cds.exists()
    assert formatted_gff.exists()

    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert ">Azolla_filiculoides_Azfi_high1" in cds_text
    assert ">Azolla_filiculoides_Azfi_low_keep" in cds_text
    assert "Azfi_low_overlap" not in cds_text

    with gzip.open(formatted_gff, "rt", encoding="utf-8") as handle:
        gff_text = handle.read()
    assert "ID=Azfi_high1" in gff_text
    assert "ID=Azfi_low_keep" in gff_text
    assert "ID=Azfi_low_overlap" not in gff_text


def test_download_manifest_fernbase_provider_accepts_markerless_top_level_genome_fasta(tmp_path):
    server_root = tmp_path / "server_root"
    species_dir = server_root / "ftp" / "Ceratopteris_richardii"
    species_dir.mkdir(parents=True, exist_ok=True)

    (species_dir / "Crichardii_676_v2.0_cds.fa").write_text(">Crich_g1.t1\nATGAA\n", encoding="utf-8")
    (species_dir / "Crichardii_676_v2.1.gene.gff3").write_text(
        "chr1\tsrc\tgene\t1\t5\t.\t+\t.\tID=Crich_g1\n",
        encoding="utf-8",
    )
    (species_dir / "Crichardii_676_v2.0.fa").write_text(">chr1\nATGCATGC\n", encoding="utf-8")

    handler = lambda *args, **kwargs: SimpleHTTPRequestHandler(*args, directory=str(server_root), **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "fernbase",
                    "id": "Ceratopteris_richardii",
                    "species_key": "Ceratopteris_richardii",
                    "cds_url": "",
                    "gff_url": "",
                    "genome_url": "",
                }
            ],
        )

        env = dict(os.environ)
        env["GG_FERNBASE_ID_URL_TEMPLATE"] = "http://127.0.0.1:{}/ftp/{{id}}/".format(server.server_port)
        download_dir = tmp_path / "download_cache"
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--provider",
                "fernbase",
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
    finally:
        server.shutdown()
        thread.join(timeout=5)

    raw_dir = download_dir / "FernBase" / "species_wise_original" / "Ceratopteris_richardii"
    assert (raw_dir / "Crichardii_676_v2.0_cds.fa").exists()
    assert (raw_dir / "Crichardii_676_v2.1.gene.gff3").exists()
    assert (raw_dir / "Crichardii_676_v2.0.fa").exists()


def test_download_manifest_veupathdb_provider_resolves_from_service(tmp_path):
    server_root = tmp_path / "server_root"
    data_dir = server_root / "common" / "downloads" / "Current_Release" / "EnuttalliP19"
    fasta_dir = data_dir / "fasta" / "data"
    gff_dir = data_dir / "gff" / "data"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    gff_dir.mkdir(parents=True, exist_ok=True)

    (fasta_dir / "AmoebaDB-68_EnuttalliP19_AnnotatedCDSs.fasta").write_text(">gene1.t1\nATGAA\n", encoding="utf-8")
    (fasta_dir / "AmoebaDB-68_EnuttalliP19_Genome.fasta").write_text(">chr1\nATGCATGC\n", encoding="utf-8")
    (fasta_dir / "AmoebaDB-68_EnuttalliP19_AnnotatedProteins.fasta").write_text(">gene1.p1\nMKK\n", encoding="utf-8")
    (gff_dir / "AmoebaDB-68_EnuttalliP19.gff").write_text(
        "chr1\tsrc\tgene\t1\t5\t.\t+\t.\tID=gene1\n",
        encoding="utf-8",
    )

    handler = lambda *args, **kwargs: _VEuPathDbFixtureHandler(*args, root_dir=server_root, **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "veupathdb",
                    "id": "EnuttalliP19",
                    "species_key": "",
                    "cds_url": "",
                    "gff_url": "",
                    "genome_url": "",
                }
            ],
        )

        env = dict(os.environ)
        env["GG_VEUPATHDB_SERVICE_BASE_URL"] = "http://127.0.0.1:{}/veupathdb/service".format(server.server_port)
        download_dir = tmp_path / "download_cache"
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--provider",
                "veupathdb",
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
    finally:
        server.shutdown()
        thread.join(timeout=5)

    raw_dir = download_dir / "VEuPathDB" / "species_wise_original" / "Entamoeba_nuttalli"
    assert (raw_dir / "Entamoeba_nuttalli.veupathdb.EnuttalliP19.cds.fa").exists()
    assert (raw_dir / "Entamoeba_nuttalli.veupathdb.EnuttalliP19.gene.gff3").exists()
    assert (raw_dir / "Entamoeba_nuttalli.veupathdb.EnuttalliP19.genome.fa").exists()


def test_download_manifest_insectbase_provider_resolves_from_api_and_formats_archived_genome(tmp_path):
    server_root = tmp_path / "server_root"
    species_token = "Abrostola_tripartita"
    data_dir = server_root / "data" / "genome" / species_token
    data_dir.mkdir(parents=True, exist_ok=True)

    (data_dir / (species_token + ".cds.fa")).write_text(
        ">gene1.t1\nATGAA\n>gene1.t2\nATGAAATGA\n",
        encoding="utf-8",
    )
    (data_dir / (species_token + ".gff3")).write_text(
        "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1\n",
        encoding="utf-8",
    )
    write_tar_bz2(
        data_dir / (species_token + ".genome.fa.tar.bz2"),
        {
            "nested/" + species_token + ".genome.fa": ">chr1 insectbase\nATGCATGC\n",
            "README.txt": "fixture\n",
        },
    )

    handler = lambda *args, **kwargs: _InsectBaseFixtureHandler(*args, root_dir=server_root, **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "insectbase",
                    "id": "IBG_00001",
                    "species_key": "",
                    "cds_url": "",
                    "gff_url": "",
                    "genome_url": "",
                }
            ],
        )

        download_dir = tmp_path / "download_cache"
        out_cds = tmp_path / "out_cds"
        out_gff = tmp_path / "out_gff"
        out_genome = tmp_path / "out_genome"
        env = dict(os.environ)
        env["GG_INSECTBASE_API_BASE_URL"] = "http://127.0.0.1:{}/api/genome".format(server.server_port)
        completed = run_script(
            "--provider",
            "insectbase",
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
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    finally:
        server.shutdown()
        thread.join(timeout=5)

    raw_dir = download_dir / "InsectBase" / "species_wise_original" / species_token
    assert (raw_dir / (species_token + ".cds.fa")).exists()
    assert (raw_dir / (species_token + ".gff3")).exists()
    archive_path = raw_dir / (species_token + ".genome.fa.tar.bz2")
    assert archive_path.exists()

    formatted_cds = list(out_cds.glob("Abrostola_tripartita*.fa.gz"))
    formatted_gff = list(out_gff.glob("Abrostola_tripartita*.gff.gz"))
    formatted_genome = list(out_genome.glob("Abrostola_tripartita*.fa.gz"))
    assert len(formatted_cds) == 1
    assert len(formatted_gff) == 1
    assert len(formatted_genome) == 1

    with gzip.open(formatted_cds[0], "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert ">Abrostola_tripartita_gene1" in cds_text
    with gzip.open(formatted_genome[0], "rt", encoding="utf-8") as handle:
        genome_text = handle.read()
    assert ">chr1" in genome_text
    assert "ATGCATGC" in genome_text
    with bz2.open(archive_path, "rb") as handle:
        assert handle.read(4) != b""


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
                elif include_type == "GENOME_GBFF":
                    files = {
                        "ncbi_dataset/data/GCA_002102555.1/genomic.gbff": "\n".join(
                            [
                                "LOCUS       ctg1               8 bp    DNA     linear   PLN 01-JAN-2000",
                                "DEFINITION  test.",
                                "ACCESSION   ctg1",
                                "VERSION     ctg1",
                                "FEATURES             Location/Qualifiers",
                                "     gene            1..8",
                                "                     /locus_tag=\"catA\"",
                                "                     /gene=\"catA\"",
                                "     CDS             1..8",
                                "                     /locus_tag=\"catA\"",
                                "                     /gene=\"catA\"",
                                "                     /protein_id=\"catA.t1\"",
                                "ORIGIN",
                                "        1 atgcatgc",
                                "//",
                                "",
                            ]
                        )
                    }
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
        gbff_path = raw_dir / "GCA_002102555.1_Catan2_genomic.gbff.gz"
        genome_path = raw_dir / "GCA_002102555.1_Catan2_genomic.fna.gz"
        assert cds_path.exists()
        assert gff_path.exists()
        assert gbff_path.exists()
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


def test_download_manifest_ncbi_uses_gbff_when_cds_and_gff_are_missing(tmp_path):
    ftp_root = tmp_path / "ftp_root"
    ftp_root.mkdir()

    gbff_content = "\n".join(
        [
            "LOCUS       ctg1               9 bp    DNA     linear   PLN 01-JAN-2000",
            "DEFINITION  test.",
            "ACCESSION   ctg1",
            "VERSION     ctg1",
            "FEATURES             Location/Qualifiers",
            "     gene            1..9",
            "                     /locus_tag=\"gene1\"",
            "                     /gene=\"gene1\"",
            "     CDS             join(1..3,7..9)",
            "                     /locus_tag=\"gene1\"",
            "                     /gene=\"gene1\"",
            "                     /protein_id=\"gene1.t1\"",
            "ORIGIN",
            "        1 atgaaattt",
            "//",
            "",
        ]
    )
    genome_content = ">ctg1\nATGAAATTT\n"

    class _NcbiGbffFallbackHandler(SimpleHTTPRequestHandler):
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
                        "esearchresult": {"idlist": ["40086895"]},
                    }
                )
                return
            if self.path.startswith("/eutils/esummary.fcgi"):
                self._send_json(
                    {
                        "header": {"type": "esummary", "version": "0.3"},
                        "result": {
                            "uids": ["40086895"],
                            "40086895": {
                                "ftppath_genbank": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/040/086/895/GCA_040086895.1_UAZ_Capgrnd_2",
                                "organism": "Capsella grandiflora",
                                "speciesname": "Capsella grandiflora",
                            },
                        },
                    }
                )
                return
            if self.path.startswith("/datasets/genome/accession/GCA_040086895.1/download"):
                query = parse_qs(urlparse(self.path).query)
                include_type = query.get("include_annotation_type", [""])[0]
                if include_type == "GENOME_GBFF":
                    files = {"ncbi_dataset/data/GCA_040086895.1/genomic.gbff": gbff_content}
                elif include_type == "GENOME_FASTA":
                    files = {"ncbi_dataset/data/GCA_040086895.1/GCA_040086895.1_UAZ_Capgrnd_2_genomic.fna": genome_content}
                else:
                    self.send_error(404, "not available")
                    return
                self._send_zip(files)
                return
            super().do_GET()

    handler = lambda *args, **kwargs: _NcbiGbffFallbackHandler(*args, root_dir=ftp_root, **kwargs)
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
                    "id": "GCA_040086895.1",
                    "species_key": "",
                    "cds_url": "",
                    "gff_url": "",
                    "gbff_url": "",
                    "genome_url": "",
                    "cds_filename": "",
                    "gff_filename": "",
                    "gbff_filename": "",
                    "genome_filename": "",
                }
            ],
        )

        download_dir = tmp_path / "download_cache"
        out_cds = tmp_path / "out_cds"
        out_gff = tmp_path / "out_gff"
        out_genome = tmp_path / "out_genome"
        species_summary = tmp_path / "gg_input_generation_species.tsv"
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
                "--species-cds-dir",
                str(out_cds),
                "--species-gff-dir",
                str(out_gff),
                "--species-genome-dir",
                str(out_genome),
                "--species-summary-output",
                str(species_summary),
            ],
            capture_output=True,
            text=True,
            check=False,
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

        raw_dir = download_dir / "NCBI_Genome" / "species_wise_original" / "Capsella_grandiflora"
        assert (raw_dir / "GCA_040086895.1_UAZ_Capgrnd_2_genomic.gbff.gz").exists()
        assert (raw_dir / "GCA_040086895.1_UAZ_Capgrnd_2_genomic.fna.gz").exists()

        formatted_cds = out_cds / "Capsella_grandiflora_GCA_040086895.1_UAZ_Capgrnd_2_genomic.derived.cds.fa.gz"
        formatted_gff = out_gff / "Capsella_grandiflora_GCA_040086895.1_UAZ_Capgrnd_2_genomic.derived.gff.gz"
        formatted_genome = out_genome / "Capsella_grandiflora_GCA_040086895.1_UAZ_Capgrnd_2_genomic.fa.gz"
        assert formatted_cds.exists()
        assert formatted_gff.exists()
        assert formatted_genome.exists()

        with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
            assert ">Capsella_grandiflora_gene1" in handle.read()
        with gzip.open(formatted_gff, "rt", encoding="utf-8") as handle:
            assert "\tCDS\t" in handle.read()
        with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == 1
        assert "derived CDS" in rows[0]["cds_input_path"]
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


def test_download_manifest_oryza_minuta_merges_public_gramene_subgenomes(tmp_path):
    server_root = tmp_path / "server_root"
    for token, seqid, gene_id, seq in (
        ("oryza_minutabb", "chrBB", "bb_gene1", "ATGAAATTT"),
        ("oryza_minutacc", "chrCC", "cc_gene1", "ATGCCCTTT"),
    ):
        cds_path = server_root / "oryza" / "tetraploids" / "fasta" / token / "cds" / "{}.{}.cds.all.fa.gz".format(
            "Oryza_minutabb" if token.endswith("bb") else "Oryza_minutacc",
            token,
        )
        gff_path = server_root / "oryza" / "tetraploids" / "gff3" / token / "{}.{}.gff3.gz".format(
            "Oryza_minutabb" if token.endswith("bb") else "Oryza_minutacc",
            token,
        )
        genome_path = server_root / "oryza" / "tetraploids" / "fasta" / token / "dna" / "{}.{}.dna.toplevel.fa.gz".format(
            "Oryza_minutabb" if token.endswith("bb") else "Oryza_minutacc",
            token,
        )
        cds_path.parent.mkdir(parents=True, exist_ok=True)
        gff_path.parent.mkdir(parents=True, exist_ok=True)
        genome_path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(cds_path, "wt", encoding="utf-8") as handle:
            handle.write(">{}.t1\n{}\n".format(gene_id, seq))
        with gzip.open(gff_path, "wt", encoding="utf-8") as handle:
            handle.write(
                "##gff-version 3\n"
                "{seqid}\tsrc\tgene\t1\t9\t.\t+\t.\tID={gene}\n"
                "{seqid}\tsrc\tmRNA\t1\t9\t.\t+\t.\tID={gene}.t1;Parent={gene}\n"
                "{seqid}\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds_{gene}.t1;Parent={gene}.t1\n".format(
                    seqid=seqid,
                    gene=gene_id,
                )
            )
        with gzip.open(genome_path, "wt", encoding="utf-8") as handle:
            handle.write(">{}\n{}\n".format(seqid, seq))

    handler = lambda *args, **kwargs: SimpleHTTPRequestHandler(*args, directory=str(server_root), **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "oryza_minuta",
                    "id": "gramene_tetraploids",
                    "species_key": "",
                }
            ],
        )

        download_dir = tmp_path / "download_cache"
        out_cds = tmp_path / "out_cds"
        out_gff = tmp_path / "out_gff"
        out_genome = tmp_path / "out_genome"
        resolved_manifest = tmp_path / "resolved.tsv"
        env = dict(os.environ)
        env["GG_ORYZA_MINUTA_GRAMENE_BASE_URL"] = "http://127.0.0.1:{}/oryza/tetraploids".format(server.server_port)

        completed = run_script(
            "--provider",
            "oryza_minuta",
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
            "--resolved-manifest-output",
            str(resolved_manifest),
            env=env,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

        raw_dir = download_dir / "OryzaMinuta" / "species_wise_original" / "Oryza_minuta"
        assert (raw_dir / "Oryza_minuta.gramene_tetraploids.merged.cds.all.fa.gz").exists()
        assert (raw_dir / "Oryza_minuta.gramene_tetraploids.merged.gff3.gz").exists()
        assert (raw_dir / "Oryza_minuta.gramene_tetraploids.merged.dna.toplevel.fa.gz").exists()

        with gzip.open(raw_dir / "Oryza_minuta.gramene_tetraploids.merged.cds.all.fa.gz", "rt", encoding="utf-8") as handle:
            merged_cds = handle.read()
        assert ">bb_gene1.t1" in merged_cds
        assert ">cc_gene1.t1" in merged_cds

        formatted_cds_candidates = sorted(out_cds.glob("*.gz"))
        formatted_gff_candidates = sorted(out_gff.glob("*.gz"))
        formatted_genome_candidates = sorted(out_genome.glob("*.gz"))
        assert len(formatted_cds_candidates) == 1
        assert len(formatted_gff_candidates) == 1
        assert len(formatted_genome_candidates) == 1

        with gzip.open(formatted_cds_candidates[0], "rt", encoding="utf-8") as handle:
            formatted_cds = handle.read()
        assert ">Oryza_minuta_bb_gene1" in formatted_cds
        assert ">Oryza_minuta_cc_gene1" in formatted_cds

        resolved_rows = list(csv.DictReader(resolved_manifest.open("rt", encoding="utf-8"), delimiter="\t"))
        assert len(resolved_rows) == 1
        assert resolved_rows[0]["provider"] == "oryza_minuta"
        assert resolved_rows[0]["id"] == "gramene_tetraploids"
        assert resolved_rows[0]["species_key"] == "Oryza_minuta"
        assert resolved_rows[0]["cds_filename"] == "Oryza_minuta.gramene_tetraploids.merged.cds.all.fa.gz"
        assert resolved_rows[0]["gff_filename"] == "Oryza_minuta.gramene_tetraploids.merged.gff3.gz"
        assert resolved_rows[0]["genome_filename"] == "Oryza_minuta.gramene_tetraploids.merged.dna.toplevel.fa.gz"
    finally:
        server.shutdown()
        server.server_close()
