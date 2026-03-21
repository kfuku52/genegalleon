from pathlib import Path
import gzip
import json
import subprocess
import sys


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "validate_cds_gff_mapping.py"


def run_script(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def write_gzip_text(path, text):
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def load_module():
    import importlib.util

    spec = importlib.util.spec_from_file_location("validate_cds_gff_mapping_module", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_species_prefix_from_name_preserves_taxonomic_qualifiers():
    mod = load_module()

    assert mod.species_prefix_from_name("Dictyostelium_cf_discoideum_demo.fa.gz") == "Dictyostelium_cf_discoideum"
    assert mod.species_prefix_from_name("Bacillus_subtilis_subsp_subtilis_demo.gff.gz") == "Bacillus_subtilis_subsp_subtilis"
    assert mod.species_prefix_from_name("Amoeba_sp_JDSRuffled_demo.fa.gz") == "Amoeba_sp_JDSRuffled"


def test_validate_cds_gff_mapping_passes_on_matching_ids(tmp_path):
    cds_dir = tmp_path / "species_cds"
    gff_dir = tmp_path / "species_gff"
    cds_dir.mkdir()
    gff_dir.mkdir()

    cds_path = cds_dir / "Arabidopsis_thaliana_demo.fa.gz"
    gff_path = gff_dir / "Arabidopsis_thaliana_demo.gff.gz"
    write_gzip_text(
        cds_path,
        (
            ">Arabidopsis_thaliana_gene1\nATGAAA\n"
            ">Arabidopsis_thaliana_gene2\nATGCCC\n"
        ),
    )
    write_gzip_text(
        gff_path,
        (
            "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=gene1.CDS1;Parent=gene1;\n"
            "chr1\tsrc\tCDS\t11\t16\t.\t+\t0\tID=gene2.CDS1;Parent=gene2;\n"
        ),
    )

    completed = run_script(
        "--species-cds-dir",
        str(cds_dir),
        "--species-gff-dir",
        str(gff_dir),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert "[Arabidopsis_thaliana] CDS-to-GFF mapping OK: 2/2 IDs" in completed.stdout
    assert "species_checked=1, species_passed=1" in completed.stdout


def test_validate_cds_gff_mapping_fails_on_missing_ids(tmp_path):
    cds_dir = tmp_path / "species_cds"
    gff_dir = tmp_path / "species_gff"
    cds_dir.mkdir()
    gff_dir.mkdir()

    cds_path = cds_dir / "Arabidopsis_thaliana_demo.fa.gz"
    gff_path = gff_dir / "Arabidopsis_thaliana_demo.gff.gz"
    write_gzip_text(
        cds_path,
        (
            ">Arabidopsis_thaliana_gene1\nATGAAA\n"
            ">Arabidopsis_thaliana_gene2\nATGCCC\n"
        ),
    )
    write_gzip_text(
        gff_path,
        "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=gene1.CDS1;Parent=gene1;\n",
    )

    completed = run_script(
        "--species-cds-dir",
        str(cds_dir),
        "--species-gff-dir",
        str(gff_dir),
        "--missing-limit",
        "5",
    )
    assert completed.returncode == 1
    assert "missing=1 sample=Arabidopsis_thaliana_gene2" in completed.stderr


def test_validate_cds_gff_mapping_allows_known_azolla_filiculoides_orphan_cds_ids(tmp_path):
    cds_dir = tmp_path / "species_cds"
    gff_dir = tmp_path / "species_gff"
    cds_dir.mkdir()
    gff_dir.mkdir()

    write_gzip_text(
        cds_dir / "Azolla_filiculoides_demo.fa.gz",
        (
            ">Azolla_filiculoides_Azfi_s0001.g000001\nATGAAA\n"
            ">Azolla_filiculoides_Azfi_s0034.g025227\nATGCCC\n"
            ">Azolla_filiculoides_Azfi_s0093.g043301\nATGCCC\n"
        ),
    )
    write_gzip_text(
        gff_dir / "Azolla_filiculoides_demo.gff.gz",
        "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=Azfi_s0001.g000001.t1.CDS1;Parent=Azfi_s0001.g000001.t1;\n",
    )

    completed = run_script(
        "--species-cds-dir",
        str(cds_dir),
        "--species-gff-dir",
        str(gff_dir),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert "[Azolla_filiculoides] CDS-to-GFF mapping OK: 1/3 IDs (allowed_missing=2)" in completed.stdout


def test_validate_cds_gff_mapping_accepts_nthreads_and_reports_it(tmp_path):
    cds_dir = tmp_path / "species_cds"
    gff_dir = tmp_path / "species_gff"
    stats_path = tmp_path / "stats.json"
    cds_dir.mkdir()
    gff_dir.mkdir()

    write_gzip_text(
        cds_dir / "Arabidopsis_thaliana_demo.fa.gz",
        ">Arabidopsis_thaliana_gene1\nATGAAA\n",
    )
    write_gzip_text(
        gff_dir / "Arabidopsis_thaliana_demo.gff.gz",
        "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=gene1.CDS1;Parent=gene1;\n",
    )
    write_gzip_text(
        cds_dir / "Oryza_sativa_demo.fa.gz",
        ">Oryza_sativa_gene1\nATGAAA\n",
    )
    write_gzip_text(
        gff_dir / "Oryza_sativa_demo.gff.gz",
        "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=gene1.CDS1;Parent=gene1;\n",
    )

    completed = run_script(
        "--species-cds-dir",
        str(cds_dir),
        "--species-gff-dir",
        str(gff_dir),
        "--nthreads",
        "2",
        "--stats-output",
        str(stats_path),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    stats = json.loads(stats_path.read_text(encoding="utf-8"))
    assert stats["nthreads"] == 2
    assert stats["species_checked"] == 2
    assert stats["species_passed"] == 2


def test_validate_cds_gff_mapping_keeps_legacy_ncpu_alias(tmp_path):
    cds_dir = tmp_path / "species_cds"
    gff_dir = tmp_path / "species_gff"
    stats_path = tmp_path / "stats.json"
    cds_dir.mkdir()
    gff_dir.mkdir()

    write_gzip_text(
        cds_dir / "Arabidopsis_thaliana_demo.fa.gz",
        ">Arabidopsis_thaliana_gene1\nATGAAA\n",
    )
    write_gzip_text(
        gff_dir / "Arabidopsis_thaliana_demo.gff.gz",
        "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=gene1.CDS1;Parent=gene1;\n",
    )

    completed = run_script(
        "--species-cds-dir",
        str(cds_dir),
        "--species-gff-dir",
        str(gff_dir),
        "--ncpu",
        "3",
        "--stats-output",
        str(stats_path),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    stats = json.loads(stats_path.read_text(encoding="utf-8"))
    assert stats["nthreads"] == 3
