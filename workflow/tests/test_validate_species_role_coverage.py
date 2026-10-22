import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from validate_species_role_coverage import validate


def fixture(tmp_path):
    dirs = {role: tmp_path / role for role in ("cds", "gff", "genome")}
    for path in dirs.values():
        path.mkdir()
    rows = []
    for species in ("Alpha_one", "Beta_two"):
        cds = dirs["cds"] / (species + ".fa")
        cds.write_text(f">{species}_gene1\nATGAAATAA\n")
        rows.append(dict(species_prefix=species, cds_output_path=str(cds),
                         gff_output_path="", genome_output_path="", gff_input_path=""))
    gff = dirs["gff"] / "Beta_two.gff"
    gff.write_text("chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=gene1\n")
    rows[1].update(gff_input_path="source.gff", gff_output_path=str(gff))
    return dirs, rows


def summary(tmp_path, rows):
    path = tmp_path / "species.tsv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_optional_roles_are_checked_per_species(tmp_path):
    dirs, rows = fixture(tmp_path)
    assert validate(summary(tmp_path, rows), dirs, ["cds"]) == dict(species=2, cds=2, gff=1, genome=0)


@pytest.mark.parametrize("change", ["missing_declared", "missing_required", "undeclared", "hidden_declaration",
                                    "misplaced", "duplicate_file", "duplicate_species", "unpaired"])
def test_mixed_roles_never_hide_missing_or_mismatched_outputs(tmp_path, change):
    dirs, rows = fixture(tmp_path)
    required = ["cds"]
    if change == "missing_declared":
        Path(rows[1]["gff_output_path"]).unlink()
    elif change == "missing_required":
        required.append("gff")
    elif change == "undeclared":
        (dirs["gff"] / "Orphan_extra.gff").write_text("orphan\n")
    elif change == "hidden_declaration":
        Path(rows[1]["gff_output_path"]).unlink()
        rows[1]["gff_output_path"] = ""
    elif change == "misplaced":
        wrong = tmp_path / "Beta_two.gff"
        wrong.write_text(Path(rows[1]["gff_output_path"]).read_text())
        rows[1]["gff_output_path"] = str(wrong)
    elif change == "duplicate_file":
        rows[1]["cds_output_path"] = rows[0]["cds_output_path"]
    elif change == "duplicate_species":
        rows[1]["species_prefix"] = rows[0]["species_prefix"]
    elif change == "unpaired":
        Path(rows[1]["cds_output_path"]).unlink()
        rows[1]["cds_output_path"] = ""
        required = []
    with pytest.raises(ValueError):
        validate(summary(tmp_path, rows), dirs, required)
