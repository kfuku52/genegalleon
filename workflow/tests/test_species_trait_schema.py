import json
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = [pytest.mark.runtime, pytest.mark.integration]

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
from species_trait_schema import schema_path, schema_payload, select_traits  # noqa: E402


def test_declared_categories_are_excluded_even_when_numeric(tmp_path):
    table = tmp_path / "traits.tsv"
    table.write_text('species\tlatitude\thabit\tcode\na\t-35\t"a\tb"\t1\nb\t42\therb\t2\n')
    schema_path(table).write_bytes(schema_payload(table.read_bytes(),
                                  {"latitude": "numeric", "habit": "text", "code": "categorical"}))
    report = select_traits(table)
    assert [row["status"] for row in report] == ["selected", "excluded", "excluded"]
    assert [row["reason"] for row in report] == ["", "declared_text", "declared_categorical"]
    with pytest.raises(ValueError, match="explicitly encode"):
        select_traits(table, "code")
    table.write_text(table.read_text().replace("-35", "-36"))
    with pytest.raises(ValueError, match="does not match"):
        select_traits(table)


@pytest.mark.parametrize("value", ["bad", "Inf", "-Inf", "TRUE"])
def test_numeric_corruption_is_never_auto_excluded(tmp_path, value):
    table = tmp_path / "traits.tsv"
    table.write_text(f"species\theight\na\t{value}\n")
    for declared in (False, True):
        if declared:
            schema_path(table).write_bytes(schema_payload(table.read_bytes(), {"height": "numeric"}))
        with pytest.raises(ValueError, match="invalid numeric|non-finite"):
            select_traits(table)


@pytest.mark.parametrize("text", ["species\tx\tx\na\t1\t2\n", "species\tx\na\t1\t2\n",
                                      "species\t x\na\t1\n"])
def test_malformed_table_rejected(tmp_path, text):
    table = tmp_path / "traits.tsv"
    table.write_text(text)
    with pytest.raises(ValueError):
        select_traits(table)


def test_legacy_numeric_table_and_aggregated_binary_keep_numeric_values(tmp_path):
    table = tmp_path / "traits.tsv"
    table.write_text("species\tlatitude\tfraction\na\t-35\t0.5\nb\tNA\t1\n")
    assert all(row["status"] == "selected" for row in select_traits(table))
    schema_path(table).write_bytes(schema_payload(table.read_bytes(), {"latitude": "numeric", "fraction": "binary"}))
    assert all(row["status"] == "selected" for row in select_traits(table))


def test_generator_schema_and_r_analysis_connection(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("provider\tid\tspecies_key\n" + "".join(f"local\t{i}\tGenus_sp{i}\n" for i in range(1, 7)))
    source = tmp_path / "source.tsv"
    source.write_text("species\tlatitude\thabit\n" + "".join(f"Genus_sp{i}\t{-40+i*10}\ttree\n" for i in range(1, 7)))
    plan = tmp_path / "plan.tsv"
    plan.write_text("database\tsource_column\toutput_trait\tvalue_type\taggregation\n"
                    "austraits\tlatitude\tlatitude\tnumeric\tmean\naustraits\thabit\thabit\ttext\tunique\n")
    sources = tmp_path / "sources.tsv"
    sources.write_text(f"database\tacquisition_mode\turi\tspecies_column\tdelimiter\naustraits\tbulk\t{source}\tspecies\ttsv\n")
    table = tmp_path / "traits.tsv"
    cmd = [sys.executable, str(SUPPORT / "generate_species_trait.py"), "--download-manifest", str(manifest),
           "--trait-plan", str(plan), "--database-sources", str(sources), "--output", str(table),
           "--downloads-dir", str(tmp_path / "downloads"), "--species-source", "download_manifest", "--strict"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads(schema_path(table).read_text())["traits"] == {"latitude": "numeric", "habit": "text"}
    tree = tmp_path / "tree.nwk"
    tree.write_text("((Genus_sp1:1,Genus_sp2:1):1,(Genus_sp3:1,Genus_sp4:1):1,(Genus_sp5:1,Genus_sp6:1):1);")
    counts = tmp_path / "counts.tsv"
    counts.write_text("Orthogroup\t" + "\t".join(f"Genus_sp{i}" for i in range(1, 7)) + "\nOG1\t1\t3\t2\t4\t3\t5\n")
    output = tmp_path / "results"
    cmd = ["Rscript", str(SUPPORT / "orthogroup_copy_number_trait_pgls.r"), f"--file_trait={table}",
           f"--file_orthogroup_copy_number={counts}", f"--file_sptree={tree}", f"--outdir={output}"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "habit\ttext\texcluded\tdeclared_text" in (output / "trait_selection.tsv").read_text()
    assert "latitude" in (output / "orthogroup_copy_number_trait_pgls.tsv").read_text()
    assert "habit" not in (output / "orthogroup_copy_number_trait_pgls.tsv").read_text()
    # No partially updated result bundle when the sidecar no longer matches.
    before = {p.name: p.read_bytes() for p in output.iterdir()}
    table.write_text(table.read_text().replace("-30", "-31"))
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode != 0
    assert {p.name: p.read_bytes() for p in output.iterdir()} == before


@pytest.mark.parametrize("alias", ["direct", "hardlink", "symlink"])
def test_selection_report_cannot_replace_input(tmp_path, alias):
    table = tmp_path / "traits.tsv"
    payload = b"species\theight\na\t10\n"
    table.write_bytes(payload)
    target = table if alias == "direct" else tmp_path / "report.tsv"
    if alias == "hardlink":
        target.hardlink_to(table)
    elif alias == "symlink":
        target.symlink_to(table)
    result = subprocess.run([sys.executable, str(SUPPORT / "species_trait_schema.py"),
                             "--table", str(table), "--report", str(target)], capture_output=True, text=True)
    assert result.returncode != 0
    assert table.read_bytes() == payload
