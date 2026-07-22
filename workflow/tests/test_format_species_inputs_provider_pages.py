# ruff: noqa: E501,E731,F403,F405

import shutil
import struct
import zlib

import pytest
from format_species_download_helpers import *


def _build_rar4_store_archive(member_name, payload):
    def block(body):
        return struct.pack("<H", zlib.crc32(body) & 0xFFFF) + body

    name = member_name.encode("utf-8")
    main_header = struct.pack("<BHHHI", 0x73, 0, 13, 0, 0)
    file_header = struct.pack(
        "<BHHIIBIIBBHI",
        0x74,
        0x8000,
        32 + len(name),
        len(payload),
        len(payload),
        3,
        zlib.crc32(payload) & 0xFFFFFFFF,
        0,
        20,
        0x30,
        len(name),
        0o100644,
    ) + name
    end_header = struct.pack("<BHH", 0x7B, 0, 7)
    return b"Rar!\x1a\x07\x00" + block(main_header) + block(file_header) + payload + block(end_header)


def test_download_manifest_resolves_coge_urls_from_id_without_templates(tmp_path):
    gff_queries = []

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
                gff_queries.append(parse_qs(parsed.query))
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
                self._send_text(
                    "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=AT1G01010\n"
                    "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=AT1G01010.1;Parent=AT1G01010\n"
                    "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=AT1G01010.1;Parent=AT1G01010.1\n"
                )
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
            timeout=TEST_COMMAND_TIMEOUT,
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

        raw_dir = download_dir / "CoGe" / "species_wise_original" / "Arabidopsis_thaliana"
        assert (raw_dir / "Arabidopsis_thaliana.coge.gid24739.cds.fa").exists()
        assert (raw_dir / "Arabidopsis_thaliana.gid24739.gff").exists()
        assert (raw_dir / "Arabidopsis_thaliana.coge.gid24739.genome.fa").exists()
        assert gff_queries
        assert gff_queries[0]["cds"] == ["1"]
        assert gff_queries[0]["annos"] == ["0"]
    finally:
        server.shutdown()
        server.server_close()


def test_download_manifest_rejects_header_only_coge_gff_before_bundle_download(tmp_path):
    requested_paths = []

    class _HeaderOnlyCoGeFixtureHandler(SimpleHTTPRequestHandler):
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
            requested_paths.append(parsed.path)
            if parsed.path == "/coge/GenomeInfo.pl" and parse_qs(parsed.query).get("fname") == ["get_gff"]:
                self._send_json(
                    {
                        "file": "Sarracenia_purpurea.gid69349.gff",
                        "files": [
                            "http://127.0.0.1:{}/files/Sarracenia_purpurea.gid69349.gff".format(
                                self.server.server_port
                            )
                        ],
                    }
                )
                return
            if parsed.path == "/files/Sarracenia_purpurea.gid69349.gff":
                self._send_text("##gff-version\t3\n##Generated by CoGe\n")
                return
            self.send_response(404)
            self.end_headers()

    server = ThreadingHTTPServer(("127.0.0.1", 0), _HeaderOnlyCoGeFixtureHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [{"provider": "coge", "id": "69349", "species_key": "Sarracenia_purpurea"}],
        )
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
                str(tmp_path / "download_cache"),
                "--download-only",
            ],
            capture_output=True,
            text=True,
            check=False,
            env=env,
            timeout=TEST_COMMAND_TIMEOUT,
        )
        assert completed.returncode == 1
        assert "header-only GFF export for gid 69349" in completed.stderr
        assert "/coge/get_seqs_for_feattype_for_genome.pl" not in requested_paths
        assert "/coge/api/v1/genomes/69349/sequence" not in requested_paths
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
        env["GG_NCBI_DATASETS_BASE_URL"] = "http://127.0.0.1:{}/datasets".format(server.server_port)
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
            timeout=TEST_COMMAND_TIMEOUT,
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
        env["GG_NCBI_DATASETS_BASE_URL"] = "http://127.0.0.1:{}/datasets".format(server.server_port)
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
            timeout=TEST_COMMAND_TIMEOUT,
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
            timeout=TEST_COMMAND_TIMEOUT,
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
            timeout=TEST_COMMAND_TIMEOUT,
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
            timeout=TEST_COMMAND_TIMEOUT,
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
            timeout=TEST_COMMAND_TIMEOUT,
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
            timeout=TEST_COMMAND_TIMEOUT,
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
            timeout=TEST_COMMAND_TIMEOUT,
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
        download_cache = tmp_path / "download_cache"
        out_cds = tmp_path / "out_cds"
        out_gff = tmp_path / "out_gff"
        out_genome = tmp_path / "out_genome"
        completed = run_script(
            "--provider",
            "all",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(download_cache),
            "--resolved-manifest-output",
            str(resolved_manifest),
            "--species-cds-dir",
            str(out_cds),
            "--species-gff-dir",
            str(out_gff),
            "--species-genome-dir",
            str(out_genome),
            env=env,
            timeout=TEST_COMMAND_TIMEOUT,
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
    raw_dir = download_cache / "Figshare" / "species_wise_original" / "Camellia_meiocarpa"
    assert (raw_dir / row["cds_filename"]).is_file()
    assert (raw_dir / row["gff_filename"]).is_file()
    assert (raw_dir / row["genome_filename"]).is_file()
    assert len(list(out_cds.glob("*.cds.fa.gz"))) == 1
    assert len(list(out_gff.glob("*.gff.gz"))) == 1


def test_download_manifest_figshare_extracts_rar_member_from_extensionless_url(tmp_path):
    if shutil.which("bsdtar") is None:
        pytest.skip("bsdtar is supplied by the GeneGalleon container runtime")

    server_root = tmp_path / "server_root"
    api_dir = server_root / "v2" / "articles"
    files_dir = server_root / "files"
    api_dir.mkdir(parents=True, exist_ok=True)
    files_dir.mkdir(parents=True, exist_ok=True)

    member_name = "data/test.txt"
    member_payload = b"original fixture\r\n"
    rar_payload = _build_rar4_store_archive(member_name, member_payload)
    (files_dir / "45433333").write_bytes(rar_payload)

    article_payload = {
        "id": 25533064,
        "title": "Chromosome-level genome assembly and annotation of Reaumuria soongarica",
        "files": [
            {
                "id": 45433333,
                "name": "Chromosome-level genome assembly and annotation of Reaumuria soongarica.rar",
                "download_url": "http://127.0.0.1:0/files/45433333",
            }
        ],
    }

    class _FigshareRarFixtureHandler(SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=str(server_root), **kwargs)

        def do_GET(self):
            if self.path == "/v2/articles/25533064":
                payload = json.dumps(article_payload).replace(":0/", ":{}".format(self.server.server_port) + "/").encode("utf-8")
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(payload)))
                self.end_headers()
                self.wfile.write(payload)
                return
            super().do_GET()

    server = ThreadingHTTPServer(("127.0.0.1", 0), _FigshareRarFixtureHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "figshare",
                    "id": "https://springernature.figshare.com/articles/dataset/25533064",
                    "species_key": "Reaumuria_soongarica",
                    "cds_filename": article_payload["files"][0]["name"],
                    "cds_archive_member": member_name,
                }
            ],
        )
        env = dict(os.environ)
        env["GG_FIGSHARE_API_BASE_URL"] = "http://127.0.0.1:{}/v2".format(server.server_port)
        completed = run_script(
            "--provider",
            "figshare",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(tmp_path / "download_cache"),
            "--download-only",
            env=env,
            timeout=TEST_COMMAND_TIMEOUT,
        )
    finally:
        server.shutdown()
        server.server_close()

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    extracted = (
        tmp_path
        / "download_cache"
        / "Figshare"
        / "species_wise_original"
        / "Reaumuria_soongarica"
        / "test.txt"
    )
    assert extracted.read_bytes() == member_payload


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
            timeout=TEST_COMMAND_TIMEOUT,
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
            timeout=TEST_COMMAND_TIMEOUT,
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
            timeout=TEST_COMMAND_TIMEOUT,
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
        env["GG_NCBI_DATASETS_BASE_URL"] = "http://127.0.0.1:{}/datasets".format(server.server_port)

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
            timeout=TEST_COMMAND_TIMEOUT,
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
            timeout=TEST_COMMAND_TIMEOUT,
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
            timeout=TEST_COMMAND_TIMEOUT,
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
