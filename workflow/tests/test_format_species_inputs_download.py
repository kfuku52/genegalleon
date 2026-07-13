# ruff: noqa: E501,F403,F405

from format_species_download_helpers import *


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
        timeout=TEST_COMMAND_TIMEOUT,
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


def test_download_manifest_xlsx_direct_catalog_sheet_fills_runtime_urls(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    species_key = "Catalog_species"
    gff_source = source_dir / "catalog_species.gene.gff3"
    genome_source = source_dir / "catalog_species.genome.fa"
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

    manifest = tmp_path / "manifest.xlsx"
    manifest_headers = [
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
    catalog_headers = [
        "choice",
        "id",
        "species",
        "species_key",
        "cds_url",
        "gff_url",
        "gbff_url",
        "genome_url",
        "cds_filename",
        "gff_filename",
        "gbff_filename",
        "genome_filename",
    ]
    make_manifest_xlsx_with_direct_catalog(
        manifest,
        manifest_headers,
        [{"provider": "direct", "id": "catalog_direct (Catalog species)"}],
        catalog_headers,
        [
            {
                "choice": "catalog_direct (Catalog species)",
                "id": "catalog_direct",
                "species": "Catalog species",
                "species_key": species_key,
                "gff_url": to_file_url(gff_source),
                "genome_url": to_file_url(genome_source),
                "gff_filename": species_key + ".catalog.gene.gff3",
                "genome_filename": species_key + ".catalog.genome.fa",
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
    assert (raw_dir / (species_key + ".catalog.gene.gff3")).exists()
    assert (raw_dir / (species_key + ".catalog.genome.fa")).exists()

    formatted_cds = out_cds / (species_key + "_catalog.gene.derived.cds.fa.gz")
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert cds_text.count(">Catalog_species_gene1") == 1
    assert "ATGTTT" in cds_text


def test_download_manifest_direct_recognizes_haplotype_genome_filename(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    species_key = "Fagopyrum_esculentum"
    cds_source = source_dir / "Fagopyrum_esculentum_hap1.cds.long.fa.gz"
    gff_source = source_dir / "Fagopyrum_esculentum_hap1.gff3.gz"
    genome_source = source_dir / "Fagopyrum_esculentum_hap1.fasta.gz"
    with gzip.open(cds_source, "wt", encoding="utf-8") as handle:
        handle.write(">gene1.t1\nATGTTTAAA\n")
    with gzip.open(gff_source, "wt", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1",
                    "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                    "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                    "",
                ]
            )
        )
    with gzip.open(genome_source, "wt", encoding="utf-8") as handle:
        handle.write(">chr1\nATGTTTAAA\n")

    manifest = tmp_path / "manifest.tsv"
    make_manifest(
        manifest,
        [
            {
                "provider": "direct",
                "id": "Fagopyrum_esculentum_Figshare_hap1",
                "species_key": species_key,
                "cds_url": to_file_url(cds_source),
                "gff_url": to_file_url(gff_source),
                "genome_url": to_file_url(genome_source),
                "cds_filename": cds_source.name,
                "gff_filename": gff_source.name,
                "genome_filename": genome_source.name,
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

    formatted_genome = out_genome / "Fagopyrum_esculentum_hap1.fa.gz"
    assert formatted_genome.exists()
    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["genome_status"] == "write"
    assert rows[0]["genome_output_path"] == str(formatted_genome)


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


def test_download_manifest_redownloads_corrupt_gzip_cache(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    species_key = "Croton_tiglium"
    cds_source = source_dir / "Croton_tiglium.cds.fa.gz"
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
            }
        ],
    )

    download_dir = tmp_path / "download_cache"
    raw_dir = download_dir / "Direct" / "species_wise_original" / species_key
    raw_dir.mkdir(parents=True)
    cached_cds = raw_dir / (species_key + ".cds.fa.gz")
    cached_cds.write_bytes(b"\x1f\x8b\x08\x00partial")

    out_cds = tmp_path / "out_cds"
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
        str(tmp_path / "out_gff"),
        "--species-genome-dir",
        str(tmp_path / "out_genome"),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert "found corrupt gzip cache" in completed.stderr
    assert "Download diagnostics:" in completed.stdout
    assert "corrupt_cache_recoveries=1" in completed.stdout
    assert list(raw_dir.glob("*.corrupt.*"))
    with gzip.open(cached_cds, "rt", encoding="utf-8") as handle:
        assert "ATGAAATTT" in handle.read()

    formatted_cds = out_cds / (species_key + "_cds.fa.gz")
    assert formatted_cds.exists()


def test_download_manifest_retries_transient_http_errors(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    species_key = "Croton_tiglium"
    cds_source = source_dir / "Croton_tiglium.cds.fa.gz"
    with gzip.open(cds_source, "wt", encoding="utf-8") as handle:
        handle.write(">ctg1.t1\nATGAA\n")

    class FlakyHandler(SimpleHTTPRequestHandler):
        attempts = 0

        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=str(source_dir), **kwargs)

        def log_message(self, format, *args):
            return

        def do_GET(self):
            if self.path.endswith("/Croton_tiglium.cds.fa.gz"):
                type(self).attempts += 1
                if type(self).attempts == 1:
                    self.send_response(502)
                    self.end_headers()
                    self.wfile.write(b"temporary upstream error")
                    return
            super().do_GET()

    server = ThreadingHTTPServer(("127.0.0.1", 0), FlakyHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "direct",
                    "id": "VVPY-Croton_tiglium",
                    "species_key": species_key,
                    "cds_url": "http://127.0.0.1:{}/Croton_tiglium.cds.fa.gz".format(server.server_port),
                    "cds_filename": species_key + ".cds.fa.gz",
                }
            ],
        )

        env = dict(os.environ)
        env["GG_DOWNLOAD_RETRY_BASE_SECONDS"] = "0"
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_PATH),
                "--provider",
                "direct",
                "--download-manifest",
                str(manifest),
                "--download-dir",
                str(tmp_path / "download_cache"),
                "--download-only",
            ],
            cwd=str(source_dir),
            capture_output=True,
            text=True,
            check=False,
            env=env,
            timeout=TEST_COMMAND_TIMEOUT,
        )
    finally:
        server.shutdown()
        thread.join(timeout=5)

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert "download attempt 1/4 failed transiently; retrying" in completed.stderr
    assert "Download diagnostics:" in completed.stdout
    assert "retries transient=1,corrupt_gzip=0" in completed.stdout
    assert FlakyHandler.attempts == 2


def test_download_manifest_uses_default_user_agent_for_direct_downloads(tmp_path):
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    species_key = "Croton_tiglium"
    cds_source = source_dir / "Croton_tiglium.cds.fa.gz"
    with gzip.open(cds_source, "wt", encoding="utf-8") as handle:
        handle.write(">ctg1.t1\nATGAA\n")

    class UserAgentHandler(SimpleHTTPRequestHandler):
        seen_user_agent = ""

        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=str(source_dir), **kwargs)

        def log_message(self, format, *args):
            return

        def do_GET(self):
            if self.path.endswith("/Croton_tiglium.cds.fa.gz"):
                type(self).seen_user_agent = self.headers.get("User-Agent", "")
                if type(self).seen_user_agent != "genegalleon-input-generation":
                    self.send_response(403)
                    self.end_headers()
                    self.wfile.write(b"missing expected user agent")
                    return
            super().do_GET()

    server = ThreadingHTTPServer(("127.0.0.1", 0), UserAgentHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        manifest = tmp_path / "manifest.tsv"
        make_manifest(
            manifest,
            [
                {
                    "provider": "direct",
                    "id": "VVPY-Croton_tiglium",
                    "species_key": species_key,
                    "cds_url": "http://127.0.0.1:{}/Croton_tiglium.cds.fa.gz".format(server.server_port),
                    "cds_filename": species_key + ".cds.fa.gz",
                }
            ],
        )

        completed = run_script(
            "--provider",
            "direct",
            "--download-manifest",
            str(manifest),
            "--download-dir",
            str(tmp_path / "download_cache"),
            "--download-only",
        )
    finally:
        server.shutdown()
        thread.join(timeout=5)

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert UserAgentHandler.seen_user_agent == "genegalleon-input-generation"


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
    assert cds_text.count(">Populus_tremula_x_Populus_alba_gene1") == 1
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
