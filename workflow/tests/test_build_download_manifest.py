from pathlib import Path
import csv
import json
import shutil
import subprocess
import sys
import zipfile

from openpyxl import load_workbook


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "build_download_manifest.py"
SMALL_DATASET_ROOT = Path(__file__).resolve().parent / "data" / "small_gfe_dataset"


def run_script(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def read_manifest(path):
    if path.suffix.lower() == ".xlsx":
        workbook = load_workbook(path, read_only=True, data_only=True)
        try:
            sheet = workbook.active
            rows = list(sheet.iter_rows(values_only=True))
        finally:
            workbook.close()
        if len(rows) == 0:
            return []
        headers = [str(value or "").strip() for value in rows[0]]
        out = []
        for values in rows[1:]:
            if values is None:
                continue
            record = {}
            nonempty = False
            for idx, key in enumerate(headers):
                if key == "":
                    continue
                value = values[idx] if idx < len(values) else ""
                text = str(value or "").strip()
                if text != "":
                    nonempty = True
                record[key] = text
            if nonempty:
                out.append(record)
        return out
    with open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_list_column_values(sheet, column):
    values = []
    row_idx = 1
    while True:
        value = sheet.cell(row=row_idx, column=column).value
        if value is None or str(value).strip() == "":
            break
        values.append(str(value).strip())
        row_idx += 1
    return values


def read_comment_vml(path):
    with zipfile.ZipFile(path) as archive:
        for name in archive.namelist():
            if name.endswith(".vml"):
                return archive.read(name).decode("utf-8")
    return ""


def test_build_download_manifest_from_small_fixture_all(tmp_path):
    out = tmp_path / "download_manifest.xlsx"
    completed = run_script(
        "--provider",
        "all",
        "--input-dir",
        str(SMALL_DATASET_ROOT),
        "--output",
        str(out),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert out.exists()

    rows = read_manifest(out)
    assert len(rows) == 1
    providers = sorted([row["provider"] for row in rows])
    assert providers == ["ensemblplants"]

    for row in rows:
        assert row["cds_url"].startswith("file://")
        assert row["gff_url"].startswith("file://")
        assert "genome_url" in row
        assert row["cds_filename"] != ""
        assert row["gff_filename"] != ""
        assert "genome_filename" in row


def test_build_download_manifest_strict_fails_on_missing_pair(tmp_path):
    dataset_copy = tmp_path / "small_dataset_copy"
    shutil.copytree(SMALL_DATASET_ROOT, dataset_copy)
    missing_gff = dataset_copy / "20230216_EnsemblPlants" / "original_files" / "Ostreococcus_lucimarinus.ASM9206v1.56.gff3"
    missing_gff.unlink()

    out = tmp_path / "download_manifest.xlsx"
    completed = run_script(
        "--provider",
        "ensemblplants",
        "--input-dir",
        str(dataset_copy),
        "--output",
        str(out),
        "--strict",
    )
    assert completed.returncode == 1, completed.stderr + "\n" + completed.stdout
    rows = read_manifest(out)
    assert len(rows) == 0


def test_build_download_manifest_ncbi_like_provider_from_species_dir_fixture(tmp_path):
    input_root = tmp_path / "dataset"
    fixtures = (
        ("NCBI_Genome", "GCF_000001405.40", "GCF_000001405.40_GRCh38.p14"),
        ("NCBI_RefSeq", "GCF_000001405.40", "GCF_000001405.40_GRCh38.p14"),
        ("NCBI_GenBank", "GCA_000001405.29", "GCA_000001405.29_GRCh38.p14"),
    )
    for provider_dir, accession, stem in fixtures:
        species_dir = input_root / provider_dir / "species_wise_original" / accession
        species_dir.mkdir(parents=True, exist_ok=True)
        (species_dir / (stem + "_cds_from_genomic.fna.gz")).write_bytes(b"dummy")
        (species_dir / (stem + "_genomic.gff.gz")).write_bytes(b"dummy")
        (species_dir / (stem + "_genomic.fna.gz")).write_bytes(b"dummy")

    out = tmp_path / "download_manifest_ncbi.xlsx"
    completed = run_script(
        "--provider",
        "ncbi",
        "--input-dir",
        str(input_root),
        "--output",
        str(out),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    rows = read_manifest(out)
    assert len(rows) == 2
    assert [row["provider"] for row in rows] == ["ncbi", "ncbi"]
    assert sorted([row["id"] for row in rows]) == ["GCA_000001405.29", "GCF_000001405.40"]
    for row in rows:
        assert row["species_key"] in ("GCA_000001405.29", "GCF_000001405.40")
        assert row["genome_url"].startswith("file://")
        assert row["genome_filename"].endswith("_genomic.fna.gz")


def test_build_download_manifest_ensembl_provider_from_flat_fixture(tmp_path):
    input_root = tmp_path / "dataset"
    ensembl_dir = input_root / "Ensembl" / "original_files"
    ensembl_dir.mkdir(parents=True)
    (ensembl_dir / "homo_sapiens.GRCh38.cds.all.fa.gz").write_bytes(b"dummy")
    (ensembl_dir / "homo_sapiens.GRCh38.112.gff3.gz").write_bytes(b"dummy")
    (ensembl_dir / "homo_sapiens.GRCh38.dna.primary_assembly.fa.gz").write_bytes(b"dummy")

    out = tmp_path / "download_manifest_ensembl.xlsx"
    completed = run_script(
        "--provider",
        "ensembl",
        "--input-dir",
        str(input_root),
        "--output",
        str(out),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    rows = read_manifest(out)
    assert len(rows) == 1
    assert rows[0]["provider"] == "ensembl"
    assert rows[0]["id"] == "homo_sapiens"
    assert rows[0]["species_key"] == "homo_sapiens"
    assert rows[0]["cds_filename"] == "homo_sapiens.GRCh38.cds.all.fa.gz"
    assert rows[0]["gff_filename"] == "homo_sapiens.GRCh38.112.gff3.gz"
    assert rows[0]["genome_filename"] == "homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"


def test_build_download_manifest_coge_and_cngb_provider_from_species_dir_fixture(tmp_path):
    input_root = tmp_path / "dataset"
    species_key = "Arabidopsis_thaliana"
    provider_dir_name = {"coge": "CoGe", "cngb": "CNGB"}

    for provider in ("coge", "cngb"):
        species_dir = input_root / provider_dir_name[provider] / "species_wise_original" / species_key
        species_dir.mkdir(parents=True, exist_ok=True)
        if provider == "coge":
            cds_name = "Arabidopsis_thaliana.coge.gid24739.cds.fa"
            gff_name = "Arabidopsis_thaliana.gid24739.gff3"
            genome_name = "Arabidopsis_thaliana.coge.gid24739.genome.fa"
        else:
            cds_name = species_key + ".cds.fa"
            gff_name = species_key + ".gene.gff3"
            genome_name = species_key + ".genome.fa"
        (species_dir / cds_name).write_text(">AT1G01010_t1\nATGAA\n", encoding="utf-8")
        (species_dir / gff_name).write_text(
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=AT1G01010\n",
            encoding="utf-8",
        )
        (species_dir / genome_name).write_text(">chr1\nATGCATGC\n", encoding="utf-8")

        out = tmp_path / ("download_manifest_{}.xlsx".format(provider))
        completed = run_script(
            "--provider",
            provider,
            "--input-dir",
            str(input_root),
            "--output",
            str(out),
        )
        assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
        rows = read_manifest(out)
        assert len(rows) == 1
        assert rows[0]["provider"] == provider
        expected_id = "24739" if provider == "coge" else species_key
        assert rows[0]["id"] == expected_id
        assert rows[0]["species_key"] == species_key
        assert rows[0]["cds_filename"] == cds_name
        assert rows[0]["gff_filename"] == gff_name
        assert rows[0]["genome_filename"] == genome_name


def test_build_download_manifest_fernbase_provider_prefers_primary_annotation_files(tmp_path):
    input_root = tmp_path / "dataset"
    species_dir = input_root / "FernBase" / "species_wise_original" / "Azolla_filiculoides"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "Azolla_filiculoides.CDS.lowconfidence_v1.1.fasta").write_text(">g1.t1\nATG\n", encoding="utf-8")
    (species_dir / "Azolla_filiculoides.CDS.highconfidence_v1.1.fasta").write_text(">g1.t1\nATGAA\n", encoding="utf-8")
    (species_dir / "Azolla_filiculoides.transcript.highconfidence_v1.1.fasta").write_text(">g1.t1\nATGAAA\n", encoding="utf-8")
    (species_dir / "Azolla_filiculoides.gene_models.lowconfidence_v1.1.gff").write_text(
        "chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=g1\n",
        encoding="utf-8",
    )
    (species_dir / "Azolla_filiculoides.gene_models.highconfidence_v1.1.gff").write_text(
        "chr1\tsrc\tgene\t1\t5\t.\t+\t.\tID=g1\n",
        encoding="utf-8",
    )
    (species_dir / "Azolla_filiculoides.genome_v1.2.fasta").write_text(">chr1\nATGCATGC\n", encoding="utf-8")

    out = tmp_path / "download_manifest_fernbase.xlsx"
    completed = run_script(
        "--provider",
        "fernbase",
        "--input-dir",
        str(input_root),
        "--output",
        str(out),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    rows = read_manifest(out)
    assert len(rows) == 1
    assert rows[0]["provider"] == "fernbase"
    assert rows[0]["id"] == "Azolla_filiculoides"
    assert rows[0]["species_key"] == "Azolla_filiculoides"
    assert rows[0]["cds_filename"] == "Azolla_filiculoides.CDS.highconfidence_v1.1.fasta"
    assert rows[0]["gff_filename"] == "Azolla_filiculoides.gene_models.highconfidence_v1.1.gff"
    assert rows[0]["genome_filename"] == "Azolla_filiculoides.genome_v1.2.fasta"


def test_build_download_manifest_xlsx_has_provider_and_id_dropdowns(tmp_path):
    out = tmp_path / "download_manifest.xlsx"
    completed = run_script(
        "--provider",
        "all",
        "--input-dir",
        str(SMALL_DATASET_ROOT),
        "--output",
        str(out),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    vml = read_comment_vml(out)
    assert "width:144px;height:79px" not in vml
    assert "visibility:hidden" in vml

    workbook = load_workbook(out)
    try:
        sheet = workbook["download_plan"]
        headers = [cell.value for cell in sheet[1]]
        assert headers[0] == "provider"
        assert headers[1] == "id"
        header_comments = {str(cell.value): (cell.comment.text if cell.comment else "") for cell in sheet[1]}
        assert set(header_comments.keys()) == set(headers)
        assert all(header_comments[name] != "" for name in headers)
        assert "Allowed values:" in header_comments["provider"]
        assert "Provider-specific identifier." in header_comments["id"]
        assert "Supported placeholders: {id}, {species_key}, {provider}." in header_comments["cds_url_template"]
        assert "provider=local" in header_comments["local_cds_path"]

        validations = list(sheet.data_validations.dataValidation)
        provider_validation = next(v for v in validations if str(v.sqref) == "A2:A5000")
        id_validation = next(v for v in validations if str(v.sqref) == "B2:B5000")
        assert "_lists!$A$1:$A$" in str(provider_validation.formula1)
        assert 'INDIRECT("id_opts_"&$A2)' in str(id_validation.formula1)

        list_sheet = workbook["_lists"]
        provider_values = [list_sheet.cell(row=i, column=1).value for i in range(1, 11)]
        assert provider_values == [
            "ensembl",
            "ensemblplants",
            "ncbi",
            "coge",
            "cngb",
            "flybase",
            "wormbase",
            "vectorbase",
            "fernbase",
            "local",
        ]
        assert "id_opts_ensembl" in workbook.defined_names
        assert "id_opts_ncbi" in workbook.defined_names
        assert "id_opts_coge" in workbook.defined_names
        assert "id_opts_fernbase" in workbook.defined_names
        assert "id_opts_local" in workbook.defined_names
    finally:
        workbook.close()


def test_build_download_manifest_xlsx_id_lists_are_provider_specific(tmp_path):
    input_root = tmp_path / "dataset"

    coge_species_gid = (
        ("Arabidopsis_thaliana", "24739"),
        ("Oryza_sativa", "24740"),
    )
    for species, gid in coge_species_gid:
        species_dir = input_root / "CoGe" / "species_wise_original" / species
        species_dir.mkdir(parents=True, exist_ok=True)
        (species_dir / (species + ".coge.gid{}.cds.fa".format(gid))).write_text(">g1.t1\nATG\n", encoding="utf-8")
        (species_dir / (species + ".gid{}.gff3".format(gid))).write_text(
            "chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=g1\n",
            encoding="utf-8",
        )

    local_species = ("Hydrocotyle_leucocephala_HAP1v2.1", "Marchantia_polymorpha_v6")
    for species in local_species:
        species_dir = input_root / "Local" / "species_wise_original" / species
        species_dir.mkdir(parents=True, exist_ok=True)
        (species_dir / (species + ".cds.fa")).write_text(">g1.t1\nATG\n", encoding="utf-8")
        (species_dir / (species + ".gene.gff3")).write_text("chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=g1\n", encoding="utf-8")

    out = tmp_path / "download_manifest_provider_lists.xlsx"
    completed = run_script(
        "--provider",
        "all",
        "--input-dir",
        str(input_root),
        "--output",
        str(out),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    workbook = load_workbook(out)
    try:
        list_sheet = workbook["_lists"]
        coge_values = read_list_column_values(list_sheet, 5)
        cngb_values = read_list_column_values(list_sheet, 6)
        fernbase_values = read_list_column_values(list_sheet, 10)
        local_values = read_list_column_values(list_sheet, 11)
        ncbi_values = read_list_column_values(list_sheet, 4)

        assert coge_values == [
            "24739 (Arabidopsis thaliana)",
            "29177 (Oryza sativa)",
            "42091 (Zea mays)",
            "78085 (Drosophila melanogaster)",
            "72417 (Caenorhabditis elegans)",
        ]
        assert cngb_values == [
            "CNA0012345 (Homo sapiens)",
            "GCF_000001405.40 (Homo sapiens)",
            "GCA_000001635.9 (Mus musculus)",
            "GCF_049306965.1 (Danio rerio)",
            "GCA_000001215.4 (Drosophila melanogaster)",
        ]
        assert fernbase_values == [
            "Azolla_filiculoides (Azolla filiculoides)",
            "Salvinia_cucullata_v2 (Salvinia cucullata v2)",
        ]
        assert local_values == [
            "Hydrocotyle_leucocephala_HAP1v2.1",
            "Marchantia_polymorpha_v6",
        ]
        assert ncbi_values == [
            "GCF_000001405.40 (Homo sapiens)",
            "GCA_000001635.9 (Mus musculus)",
            "GCF_049306965.1 (Danio rerio)",
            "GCA_000001215.4 (Drosophila melanogaster)",
            "GCF_000002985.6 (Caenorhabditis elegans)",
        ]
    finally:
        workbook.close()


def test_build_download_manifest_xlsx_prefers_snapshot_for_full_providers(tmp_path):
    out = tmp_path / "download_manifest.xlsx"
    snapshot = tmp_path / "id_options_snapshot.json"
    snapshot.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "providers": {
                    "ensembl": [
                        {"id": "homo_sapiens", "species": "Homo sapiens"},
                        {"id": "pan_troglodytes", "species": "Pan troglodytes"},
                    ],
                    "ensemblplants": [
                        {"id": "arabidopsis_thaliana", "species": "Arabidopsis thaliana"},
                    ],
                    "flybase": [
                        {"id": "dmel_r6.66", "species": "Drosophila melanogaster"},
                    ],
                    "wormbase": [
                        {"id": "caenorhabditis_elegans_prjna13758", "species": "Caenorhabditis elegans"},
                    ],
                    "vectorbase": [
                        {"id": "AgambiaePEST", "species": "Anopheles gambiae"},
                    ],
                    "fernbase": [
                        {"id": "Ceratopteris_richardii", "species": "Ceratopteris richardii"},
                    ],
                    "local": [
                        {"id": "/data/local_species_1", "species": "Local species 1"},
                        {"id": "/data/local_species_2", "species": "Local species 2"},
                    ],
                    "ncbi": [
                        {"id": "FAKE_NCBI_ID", "species": "Fake species"},
                    ],
                    "coge": [
                        {"id": "99999", "species": "Fake coge species"},
                    ],
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    completed = run_script(
        "--provider",
        "all",
        "--input-dir",
        str(SMALL_DATASET_ROOT),
        "--output",
        str(out),
        "--id-options-snapshot",
        str(snapshot),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    workbook = load_workbook(out)
    try:
        list_sheet = workbook["_lists"]
        ensembl_values = read_list_column_values(list_sheet, 2)
        ensemblplants_values = read_list_column_values(list_sheet, 3)
        ncbi_values = read_list_column_values(list_sheet, 4)
        coge_values = read_list_column_values(list_sheet, 5)
        cngb_values = read_list_column_values(list_sheet, 6)
        flybase_values = read_list_column_values(list_sheet, 7)
        wormbase_values = read_list_column_values(list_sheet, 8)
        vectorbase_values = read_list_column_values(list_sheet, 9)
        fernbase_values = read_list_column_values(list_sheet, 10)
        local_values = read_list_column_values(list_sheet, 11)

        assert ensembl_values == [
            "homo_sapiens (Homo sapiens)",
            "pan_troglodytes (Pan troglodytes)",
        ]
        assert ensemblplants_values == [
            "arabidopsis_thaliana (Arabidopsis thaliana)",
        ]
        assert "FAKE_NCBI_ID (Fake species)" not in ncbi_values
        assert ncbi_values == [
            "GCF_000001405.40 (Homo sapiens)",
            "GCA_000001635.9 (Mus musculus)",
            "GCF_049306965.1 (Danio rerio)",
            "GCA_000001215.4 (Drosophila melanogaster)",
            "GCF_000002985.6 (Caenorhabditis elegans)",
        ]
        assert coge_values == [
            "24739 (Arabidopsis thaliana)",
            "29177 (Oryza sativa)",
            "42091 (Zea mays)",
            "78085 (Drosophila melanogaster)",
            "72417 (Caenorhabditis elegans)",
        ]
        assert cngb_values == [
            "CNA0012345 (Homo sapiens)",
            "GCF_000001405.40 (Homo sapiens)",
            "GCA_000001635.9 (Mus musculus)",
            "GCF_049306965.1 (Danio rerio)",
            "GCA_000001215.4 (Drosophila melanogaster)",
        ]
        assert flybase_values == ["dmel_r6.66 (Drosophila melanogaster)"]
        assert wormbase_values == ["caenorhabditis_elegans_prjna13758 (Caenorhabditis elegans)"]
        assert vectorbase_values == ["AgambiaePEST (Anopheles gambiae)"]
        assert fernbase_values == ["Ceratopteris_richardii (Ceratopteris richardii)"]
        assert local_values == ["/data/local_species_1", "/data/local_species_2"]
    finally:
        workbook.close()


def test_build_download_manifest_local_provider_from_phytozome_like_fixture(tmp_path):
    input_root = tmp_path / "dataset"
    local_species_dir = input_root / "Local" / "species_wise_original" / "Hydrocotyle_leucocephala_HAP1v2.1"
    local_species_dir.mkdir(parents=True)

    source_dir = SMALL_DATASET_ROOT / "Phytozome" / "species_wise_original" / "Hydrocotyle_leucocephala_HAP1v2.1"
    for name in (
        "HleucocephalaHAP1_768_v2.1.cds_primaryTranscriptOnly.fa",
        "HleucocephalaHAP1_768_v2.1.gene.gff3",
    ):
        shutil.copy2(source_dir / name, local_species_dir / name)

    out = tmp_path / "download_manifest_local.xlsx"
    completed = run_script(
        "--provider",
        "local",
        "--input-dir",
        str(input_root),
        "--output",
        str(out),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    rows = read_manifest(out)
    assert len(rows) == 1
    assert rows[0]["provider"] == "local"
    assert rows[0]["id"] == "Hydrocotyle_leucocephala_HAP1v2.1"
    assert rows[0]["cds_filename"].endswith(".cds_primaryTranscriptOnly.fa")
    assert rows[0]["gff_filename"].endswith(".gene.gff3")


def test_build_download_manifest_all_rows_keep_local_last(tmp_path):
    input_root = tmp_path / "dataset"

    ncbi_species_dir = input_root / "NCBI_Genome" / "species_wise_original" / "GCF_000001405.40"
    ncbi_species_dir.mkdir(parents=True, exist_ok=True)
    (ncbi_species_dir / "GCF_000001405.40_GRCh38.p14_cds_from_genomic.fna.gz").write_bytes(b"dummy")
    (ncbi_species_dir / "GCF_000001405.40_GRCh38.p14_genomic.gff.gz").write_bytes(b"dummy")
    (ncbi_species_dir / "GCF_000001405.40_GRCh38.p14_genomic.fna.gz").write_bytes(b"dummy")

    local_species_dir = input_root / "Local" / "species_wise_original" / "Hydrocotyle_leucocephala_HAP1v2.1"
    local_species_dir.mkdir(parents=True, exist_ok=True)
    (local_species_dir / "HleucocephalaHAP1_768_v2.1.cds_primaryTranscriptOnly.fa").write_text(">g1\nATG\n", encoding="utf-8")
    (local_species_dir / "HleucocephalaHAP1_768_v2.1.gene.gff3").write_text(
        "chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=g1\n",
        encoding="utf-8",
    )

    out = tmp_path / "download_manifest.tsv"
    completed = run_script(
        "--provider",
        "all",
        "--input-dir",
        str(input_root),
        "--output",
        str(out),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    rows = read_manifest(out)
    assert [row["provider"] for row in rows] == ["ncbi", "local"]
