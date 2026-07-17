import gzip
import json
import subprocess
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
FORMAT_SCRIPT = SUPPORT_DIR / "format_species_inputs.py"
VALIDATE_SCRIPT = SUPPORT_DIR / "validate_cds_gff_mapping.py"


def load_format_module():
    spec = spec_from_file_location("format_species_inputs_gff_repair_test", FORMAT_SCRIPT)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_gzip_text(path, text):
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def read_gzip_text(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return handle.read()


def reaumuria_fixture_text():
    return (
        "##gff-version 3\n"
        "Chr06\ttransdecoder\tgene\t1\t9\t.\t+\t.\tID=MSTRG.10121_Chr06_+;\n"
        "Chr06\ttransdecoder\tmRNA\t1\t9\t.\t+\t.\tID=MSTRG.10121.1.p1;Parent=MSTRG.10121_Chr06_+;\n"
        "Chr06\ttransdecoder\tCDS\t1\t9\t.\t+\t0\tID=cds.MSTRG.10121.1.p1;Parent=MSTRG.10121.1.p1;\n"
        "UnChr100\tGAF\tgene\t20\t28\t.\t-\t.\tID=gmmUnChr100G000010;Name=gene_0;\n"
        "UnChr100\tGeMoMa\tmRNA\t20\t28\t.\t-\t.\tID=gmmUnChr100G000010.1;Name=reference_transcript;Parent=gene_0;\n"
        "UnChr100\tGeMoMa\tCDS\t20\t28\t.\t-\t0\tParent=gmmUnChr100G000010.1;\n"
    )


def test_safe_repair_aligns_transdecoder_and_gemoma_gene_ids(tmp_path):
    mod = load_format_module()
    source_gff = tmp_path / "source.gff"
    source_gff.write_text(reaumuria_fixture_text(), encoding="utf-8")
    cds_dir = tmp_path / "species_cds"
    gff_dir = tmp_path / "species_gff"
    cds_dir.mkdir()
    gff_dir.mkdir()
    cds_path = cds_dir / "Reaumuria_soongarica_demo.fa.gz"
    output_gff = gff_dir / "Reaumuria_soongarica_demo.gff.gz"
    write_gzip_text(
        cds_path,
        (
            ">Reaumuria_soongarica_MSTRG.10121_Chr06\nATGAAATTT\n"
            ">Reaumuria_soongarica_gene_0\nATGCCCTAA\n"
        ),
    )

    audit = mod.write_repaired_gff(
        source_gff,
        cds_path,
        output_gff,
        "Reaumuria_soongarica",
        "safe",
    )

    repaired = read_gzip_text(output_gff)
    assert "ID=MSTRG.10121_Chr06;" in repaired
    assert "Parent=MSTRG.10121_Chr06;" in repaired
    assert "ID=MSTRG.10121_Chr06_+;" not in repaired
    assert "ID=gene_0;Name=gene_0;" in repaired
    assert "ID=gmmUnChr100G000010;Name=gene_0;" not in repaired
    assert audit["renamed_gene_ids"] == 2
    assert audit["changed_references"] == 1
    assert [item["target_id"] for item in audit["repairs"]] == [
        "MSTRG.10121_Chr06",
        "gene_0",
    ]

    completed = subprocess.run(
        [
            sys.executable,
            str(VALIDATE_SCRIPT),
            "--species-cds-dir",
            str(cds_dir),
            "--species-gff-dir",
            str(gff_dir),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert "CDS-to-GFF mapping OK: 2/2 IDs" in completed.stdout


def test_safe_repair_does_not_guess_ambiguous_gene_alias(tmp_path):
    mod = load_format_module()
    source_gff = tmp_path / "ambiguous.gff"
    source_gff.write_text(
        "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=raw1;Name=gene_a;Alias=gene_b\n",
        encoding="utf-8",
    )
    cds_path = tmp_path / "Species_name.fa.gz"
    output_gff = tmp_path / "Species_name.gff.gz"
    write_gzip_text(
        cds_path,
        ">Species_name_gene_a\nATG\n>Species_name_gene_b\nATG\n",
    )

    audit = mod.write_repaired_gff(source_gff, cds_path, output_gff, "Species_name", "safe")

    assert audit["renamed_gene_ids"] == 0
    assert audit["ambiguous_count"] == 1
    assert "ID=raw1;" in read_gzip_text(output_gff)
    with pytest.raises(ValueError, match="GFF repair is ambiguous"):
        mod.write_repaired_gff(source_gff, cds_path, output_gff, "Species_name", "strict")


def test_safe_repair_rejects_target_id_collision(tmp_path):
    mod = load_format_module()
    source_gff = tmp_path / "collision.gff"
    source_gff.write_text(
        (
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=raw_gene;Name=gene_a\n"
            "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=gene_a;Parent=raw_gene\n"
        ),
        encoding="utf-8",
    )
    cds_path = tmp_path / "Species_name.fa.gz"
    output_gff = tmp_path / "Species_name.gff.gz"
    write_gzip_text(cds_path, ">Species_name_gene_a\nATG\n")

    audit = mod.write_repaired_gff(source_gff, cds_path, output_gff, "Species_name", "safe")

    assert audit["renamed_gene_ids"] == 0
    assert audit["collision_count"] == 1
    assert "ID=raw_gene;Name=gene_a" in read_gzip_text(output_gff)


def test_format_gff_tracks_repair_cache_and_refreshes_legacy_output(tmp_path):
    mod = load_format_module()
    source_gff = tmp_path / "broken.gff"
    source_gff.write_text(reaumuria_fixture_text(), encoding="utf-8")
    cds_path = tmp_path / "Reaumuria_soongarica.fa.gz"
    write_gzip_text(
        cds_path,
        (
            ">Reaumuria_soongarica_MSTRG.10121_Chr06\nATGAAATTT\n"
            ">Reaumuria_soongarica_gene_0\nATGCCCTAA\n"
        ),
    )
    output_dir = tmp_path / "species_gff"
    output_dir.mkdir()
    task = {
        "provider": "figshare",
        "species_key": "Reaumuria_soongarica",
        "species_prefix": "Reaumuria_soongarica",
        "gff_path": source_gff,
        "gbff_path": None,
        "gff_repair_mode": "safe",
    }

    first = mod.format_gff(task, output_dir, False, False, formatted_cds_path=cds_path)
    second = mod.format_gff(task, output_dir, False, False, formatted_cds_path=cds_path)

    assert first["status"] == "write"
    assert first["repair_status"] == "repaired"
    assert second["status"] == "skip"
    audit_path = Path(str(first["output_path"]) + ".repair.json")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    assert audit["repair_version"] == mod.GFF_REPAIR_VERSION

    audit_path.unlink()
    write_gzip_text(first["output_path"], reaumuria_fixture_text())
    refreshed = mod.format_gff(task, output_dir, False, False, formatted_cds_path=cds_path)
    assert refreshed["status"] == "write"
    assert "ID=MSTRG.10121_Chr06;" in read_gzip_text(first["output_path"])


def test_off_mode_preserves_source_gene_ids(tmp_path):
    mod = load_format_module()
    source_gff = tmp_path / "source.gff"
    source_gff.write_text(reaumuria_fixture_text(), encoding="utf-8")
    cds_path = tmp_path / "Reaumuria_soongarica.fa.gz"
    output_gff = tmp_path / "Reaumuria_soongarica.gff.gz"
    write_gzip_text(cds_path, ">Reaumuria_soongarica_MSTRG.10121_Chr06\nATG\n")

    audit = mod.write_repaired_gff(
        source_gff,
        cds_path,
        output_gff,
        "Reaumuria_soongarica",
        "off",
    )

    assert audit["status"] == "off"
    assert audit["renamed_gene_ids"] == 0
    assert "ID=MSTRG.10121_Chr06_+;" in read_gzip_text(output_gff)
