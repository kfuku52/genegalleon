# ruff: noqa: E501,E731,F403,F405

from format_species_download_helpers import *


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
            timeout=TEST_COMMAND_TIMEOUT,
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
    assert "Download diagnostics:" in completed.stdout
    assert "stale_locks recovered=1" in completed.stdout


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
        timeout=TEST_COMMAND_TIMEOUT,
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
        "ensemblmetazoa": "Anopheles_gambiae",
        "ensemblprotists": "Phytophthora_parasitica",
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
        "ensemblmetazoa": "anopheles_gambiae",
        "ensemblprotists": "phytophthora_parasitica",
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
        "ensemblmetazoa": "GG_ENSEMBLMETAZOA",
        "ensemblprotists": "GG_ENSEMBLPROTISTS",
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
        timeout=TEST_COMMAND_TIMEOUT,
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    expected_roots = {
        "ensembl": download_dir / "Ensembl" / "original_files",
        "ensemblplants": download_dir / "20230216_EnsemblPlants" / "original_files",
        "ensemblmetazoa": download_dir / "EnsemblMetazoa" / "original_files",
        "ensemblprotists": download_dir / "EnsemblProtists" / "original_files",
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
            timeout=TEST_COMMAND_TIMEOUT,
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
            timeout=TEST_COMMAND_TIMEOUT,
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
            timeout=TEST_COMMAND_TIMEOUT,
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
            timeout=TEST_COMMAND_TIMEOUT,
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
            timeout=TEST_COMMAND_TIMEOUT,
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
    gff_content = (
        "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene-gene1;locus_tag=gene1\n"
        "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=rna-1;Parent=gene-gene1;locus_tag=gene1\n"
        "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds-1;Parent=rna-1;protein_id=NP_000001.1;locus_tag=gene1\n"
    )
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
            ],
            capture_output=True,
            text=True,
            check=False,
            env=env,
            timeout=TEST_COMMAND_TIMEOUT,
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
        assert ">Homo_sapiens_gene1" in cds_text
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
        env["GG_NCBI_DATASETS_BASE_URL"] = "http://127.0.0.1:{}/datasets".format(server.server_port)

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
                timeout=TEST_COMMAND_TIMEOUT,
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
        env["GG_NCBI_DATASETS_BASE_URL"] = "http://127.0.0.1:{}/datasets".format(server.server_port)
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
            timeout=TEST_COMMAND_TIMEOUT,
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
            timeout=TEST_COMMAND_TIMEOUT,
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
            timeout=TEST_COMMAND_TIMEOUT,
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
        handle.write(
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene-111;Dbxref=GeneID:111;gene=ABC1\n"
            "chr1\tsrc\tmRNA\t1\t6\t.\t+\t.\tID=rna-XP_1.1;Parent=gene-111\n"
            "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cds-XP_1.1;Parent=rna-XP_1.1;Dbxref=GeneID:111;gene=ABC1;protein_id=XP_1.1\n"
            "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=rna-XP_1.2;Parent=gene-111\n"
            "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds-XP_1.2;Parent=rna-XP_1.2;Dbxref=GeneID:111;gene=ABC1;protein_id=XP_1.2\n"
        )
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
    with open(str(formatted_cds) + ".gff-grouping.json", "rt", encoding="utf-8") as handle:
        audit = json.load(handle)
    assert audit["grouping_source"] == "gff"
    assert audit["stats"]["mapped"] == 2
    assert audit["stats"]["unmapped"] == 0


def test_download_manifest_ncbi_gene_level_aggregate_prefers_locus_tag_over_geneid(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    cds_source = source_dir / "ncbi_cds.fna.gz"
    gff_source = source_dir / "ncbi_genomic.gff.gz"
    genome_source = source_dir / "ncbi_genomic.fna.gz"
    with gzip.open(cds_source, "wt", encoding="utf-8") as handle:
        handle.write(
            ">lcl|NC_003070.9_cds_NP_171609.1_1 [gene=NAC001] [locus_tag=AT1G01010] [db_xref=Araport:AT1G01010,TAIR:AT1G01010,GeneID:839580] [protein_id=NP_171609.1]\n"
            "ATGAA\n"
            ">lcl|NC_003070.9_cds_NP_001000000.1_2 [gene=NAC001] [locus_tag=AT1G01010] [db_xref=Araport:AT1G01010,GeneID:839580] [protein_id=NP_001000000.1]\n"
            "ATGAAATTT\n"
        )
    with gzip.open(gff_source, "wt", encoding="utf-8") as handle:
        handle.write(
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene-AT1G01010;Dbxref=Araport:AT1G01010,TAIR:AT1G01010,GeneID:839580;gene=NAC001;locus_tag=AT1G01010\n"
        )
    with gzip.open(genome_source, "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGCATGC\n")

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "ncbi",
                "species_key": "Arabidopsis_thaliana",
                "id": "GCF_000001735.4",
                "cds_url": to_file_url(cds_source),
                "gff_url": to_file_url(gff_source),
                "genome_url": to_file_url(genome_source),
                "cds_filename": "GCF_000001735.4_TAIR10.1_cds_from_genomic.fna.gz",
                "gff_filename": "GCF_000001735.4_TAIR10.1_genomic.gff.gz",
                "genome_filename": "GCF_000001735.4_TAIR10.1_genomic.fna.gz",
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

    formatted_cds = out_cds / "Arabidopsis_thaliana_GCF_000001735.4_TAIR10.1_cds_from_genomic.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert text.count(">Arabidopsis_thaliana_AT1G01010") == 1
    assert ">Arabidopsis_thaliana_GeneID839580" not in text
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
        handle.write(
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene-ENSG00000111111;Dbxref=Ensembl:ENSG00000111111\n"
            "chr1\tsrc\tmRNA\t1\t5\t.\t+\t.\tID=rna-XP_1.1;Parent=gene-ENSG00000111111\n"
            "chr1\tsrc\tCDS\t1\t5\t.\t+\t0\tID=cds-XP_1.1;Parent=rna-XP_1.1;protein_id=XP_1.1\n"
            "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=rna-XP_1.2;Parent=gene-ENSG00000111111\n"
            "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds-XP_1.2;Parent=rna-XP_1.2;protein_id=XP_1.2\n"
        )
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
            timeout=TEST_COMMAND_TIMEOUT,
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
