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
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


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
        "species_key\tcds_url\n"
        "Microglena_spYARC_MicrYARC1\t{}\n".format(to_file_url(cds_source)),
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
    assert "at least one of provider/id columns" in completed.stderr


def test_download_manifest_infers_provider_from_id_when_provider_column_missing(tmp_path):
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
            "id\tspecies_key\tcds_url\tgff_url\tgenome_url\tcds_filename\tgff_filename\tgenome_filename\n"
            "GCF_000001405.40\tHomo_sapiens\t{cds}\t{gff}\t{genome}\tGCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz\tGCF_000001405.40_GRCh38.p14_genomic.gff.gz\tGCF_000001405.40_GRCh38.p14_genomic.fna.gz\n"
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
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    raw_dir = download_dir / "NCBI_Genome" / "species_wise_original" / "Homo_sapiens"
    assert (raw_dir / "GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz").exists()
    assert (raw_dir / "GCF_000001405.40_GRCh38.p14_genomic.gff.gz").exists()
    assert (raw_dir / "GCF_000001405.40_GRCh38.p14_genomic.fna.gz").exists()
    assert "provider was inferred as 'ncbi'" in completed.stderr


def test_download_manifest_supports_coge_and_cngb_with_id_inference(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    species_key = "Arabidopsis_thaliana"

    provider_layout = {
        "coge": "CoGe",
        "cngb": "CNGB",
    }
    provider_ids = {
        "coge": "coge:arabidopsis_thaliana_v1",
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
                    "provider": "",
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
        assert "provider was inferred as '{}'".format(provider) in completed.stderr


def test_download_manifest_recovers_stale_lock_file(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    cds_source = source_dir / "phyco_cds.fasta"
    gff_source = source_dir / "phyco_gene.gff3"
    cds_source.write_text(">jgi|X|1|mRNA.A\nATG\n", encoding="utf-8")
    gff_source.write_text("scaf\tsrc\tgene\t1\t3\t.\t+\t.\tID=A\n", encoding="utf-8")

    manifest = tmp_path / "manifest.tsv"
    species_key = "Microglena_spYARC_MicrYARC1"
    cds_name = "MicrYARC1_GeneCatalog_CDS_20220803.fasta"
    gff_name = "MicrYARC1_GeneCatalog_genes_20220803.gff3"
    make_manifest(
        manifest,
        [
            {
                "provider": "phycocosm",
                "species_key": species_key,
                "cds_url": to_file_url(cds_source),
                "gff_url": to_file_url(gff_source),
                "cds_filename": cds_name,
                "gff_filename": gff_name,
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    raw_dir = download_dir / "PhycoCosm" / "species_wise_original" / species_key
    raw_dir.mkdir(parents=True)
    stale_lock = raw_dir / (cds_name + ".lock")
    stale_lock.write_text("999999\t0\tstale\n", encoding="utf-8")

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
    assert (raw_dir / cds_name).exists()
    assert (raw_dir / gff_name).exists()
    assert "recovered stale lock" in completed.stderr.lower()


def test_download_manifest_resolves_urls_from_id_templates_for_non_ncbi(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    species_id = "Microglena_spYARC_MicrYARC1"
    cds_source = source_dir / (species_id + ".cds.fasta")
    gff_source = source_dir / (species_id + ".genes.gff3")
    cds_source.write_text(">jgi|X|1|mRNA.A\nATG\n", encoding="utf-8")
    gff_source.write_text("scaf\tsrc\tgene\t1\t3\t.\t+\t.\tID=A\n", encoding="utf-8")

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "phycocosm",
                "id": species_id,
                "species_key": "",
                "cds_url": "",
                "gff_url": "",
                "cds_filename": "MicrYARC1_GeneCatalog_CDS_20220803.fasta",
                "gff_filename": "MicrYARC1_GeneCatalog_genes_20220803.gff3",
                "cds_url_template": source_dir.resolve().as_uri() + "/{id}.cds.fasta",
                "gff_url_template": source_dir.resolve().as_uri() + "/{id}.genes.gff3",
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
    raw_dir = download_dir / "PhycoCosm" / "species_wise_original" / species_id
    assert (raw_dir / "MicrYARC1_GeneCatalog_CDS_20220803.fasta").exists()
    assert (raw_dir / "MicrYARC1_GeneCatalog_genes_20220803.gff3").exists()


def test_download_manifest_resolves_urls_from_id_for_all_non_ncbi_providers_via_env_templates(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()

    rows = []
    provider_species = {
        "ensemblplants": "Ostreococcus_lucimarinus",
        "phycocosm": "Microglena_spYARC_MicrYARC1",
        "phytozome": "Hydrocotyle_leucocephala_HleucocephalaHAP1_768_v2.1",
        "coge": "Arabidopsis_thaliana",
        "cngb": "Arabidopsis_thaliana",
    }
    provider_ids = {
        "ensemblplants": "Ostreococcus_lucimarinus",
        "phycocosm": "Microglena_spYARC_MicrYARC1",
        "phytozome": "Hydrocotyle_leucocephala_HleucocephalaHAP1_768_v2.1",
        "coge": "coge:arabidopsis_thaliana_v1",
        "cngb": "cngb:arabidopsis_thaliana_v1",
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
        "ensemblplants": "GG_ENSEMBLPLANTS",
        "phycocosm": "GG_PHYCOCOSM",
        "phytozome": "GG_PHYTOZOME",
        "coge": "GG_COGE",
        "cngb": "GG_CNGB",
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
        "ensemblplants": download_dir / "20230216_EnsemblPlants" / "original_files",
        "phycocosm": download_dir / "PhycoCosm" / "species_wise_original" / provider_species["phycocosm"],
        "phytozome": download_dir / "Phytozome" / "species_wise_original" / provider_species["phytozome"],
        "coge": download_dir / "CoGe" / "species_wise_original" / provider_species["coge"],
        "cngb": download_dir / "CNGB" / "species_wise_original" / provider_species["cngb"],
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
                    "id": "GCF_000001405.40",
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
            "id\tspecies_key\nhttps://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/\tHomo_sapiens\n",
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
        assert "provider was inferred as 'ncbi'" in completed.stderr
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
