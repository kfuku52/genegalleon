import csv
import shlex
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path

import pytest

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))

import amalgkit_metadata_accessions as ama  # noqa: E402


def write_tsv(path, fieldnames, rows):
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_missing_accessions_preserves_requested_order(tmp_path, capsys):
    metadata = tmp_path / "metadata.tsv"
    accessions = tmp_path / "accessions.txt"
    write_tsv(metadata, ["run", "lib_source"], [{"run": "SRR1", "lib_source": "TRANSCRIPTOMIC"}])
    accessions.write_text("SRR1\nSRR2\nSRR2\nSRR3\n", encoding="utf-8")

    ama.print_missing_accessions(metadata, accessions)

    assert capsys.readouterr().out.splitlines() == ["SRR2", "SRR3"]


def test_extract_transcriptomic_rows_filters_requested_runs(tmp_path):
    metadata = tmp_path / "metadata.tsv"
    accessions = tmp_path / "accessions.txt"
    output = tmp_path / "relaxed.tsv"
    write_tsv(
        metadata,
        ["run", "lib_source", "lib_strategy"],
        [
            {"run": "SRR1", "lib_source": "TRANSCRIPTOMIC", "lib_strategy": "RNA-Seq"},
            {"run": "SRR2", "lib_source": "GENOMIC", "lib_strategy": "WGS"},
            {"run": "SRR3", "lib_source": "", "lib_strategy": "CLONE"},
        ],
    )
    accessions.write_text("SRR1\nSRR2\nSRR3\n", encoding="utf-8")

    ama.extract_transcriptomic_rows(metadata, accessions, output)

    assert [row["run"] for row in read_tsv(output)] == ["SRR1", "SRR3"]


def test_merge_metadata_tables_by_run_keeps_primary_rows_first(tmp_path):
    primary = tmp_path / "primary.tsv"
    extra = tmp_path / "extra.tsv"
    output = tmp_path / "merged.tsv"
    write_tsv(primary, ["run", "species"], [{"run": "SRR1", "species": "A"}])
    write_tsv(
        extra,
        ["run", "species", "extra"],
        [
            {"run": "SRR1", "species": "A2", "extra": "x"},
            {"run": "SRR2", "species": "B", "extra": "y"},
        ],
    )

    ama.merge_metadata_tables_by_run(primary, extra, output)

    rows = read_tsv(output)
    assert [row["run"] for row in rows] == ["SRR1", "SRR2"]
    assert rows[0]["species"] == "A"
    assert rows[1]["extra"] == "y"


def test_species_filter_matches_the_scientific_name_column_only(tmp_path):
    source, output = tmp_path / "source.tsv", tmp_path / "selected.tsv"
    fields = ["run", "scientific_name", "description"]
    write_tsv(source, fields, [
        {"run": "SRR1", "scientific_name": "Plant one", "description": "quoted\tfield\nand line"},
        {"run": "SRR2", "scientific_name": "Plant oneformis", "description": ""},
        {"run": "SRR3", "scientific_name": "Other plant", "description": "comparison with Plant one"},
        {"run": "SRR4", "scientific_name": "Plant_one", "description": ""},
    ])
    ama.extract_species_rows(source, "Plant_one", output)
    assert [r["run"] for r in read_tsv(output)] == ["SRR1", "SRR4"]
    assert read_tsv(output)[0]["description"] == "quoted\tfield\nand line"
    ama.extract_species_rows(source, "Missing plant", output)
    assert not read_tsv(output) and output.read_text().splitlines() == ["\t".join(fields)]


@pytest.mark.parametrize("header", ["run\tdescription", "run\tscientific_name\tscientific_name"])
def test_missing_or_ambiguous_species_column_preserves_previous_output(tmp_path, header):
    source, output = tmp_path / "source.tsv", tmp_path / "selected.tsv"
    source.write_text(header + "\n")
    output.write_text("previous verified output\n")
    with pytest.raises(SystemExit, match="required 'scientific_name'"):
        ama.extract_species_rows(source, "Plant one", output)
    assert output.read_text() == "previous verified output\n"


def test_streamed_extraction_supports_in_place_output_and_preserves_mode(tmp_path):
    source, accessions = tmp_path / "source.tsv", tmp_path / "runs.txt"
    write_tsv(source, ["run", "lib_source"], [
        {"run": "SRR1", "lib_source": "TRANSCRIPTOMIC"}, {"run": "SRR2", "lib_source": "GENOMIC"}])
    source.chmod(0o640)
    accessions.write_text("SRR1\nSRR2\n")
    ama.extract_transcriptomic_rows(source, accessions, source)
    assert [r["run"] for r in read_tsv(source)] == ["SRR1"]
    assert source.stat().st_mode & 0o777 == 0o640
    assert not list(tmp_path.glob("*.tmp"))


def test_streamed_extraction_does_not_publish_partial_output_on_bad_selected_row(tmp_path):
    source, output, accessions = tmp_path / "source.tsv", tmp_path / "out.tsv", tmp_path / "runs.txt"
    source.write_text("run\tlib_source\nSRR1\tTRANSCRIPTOMIC\nSRR2\tTRANSCRIPTOMIC\textra\n")
    accessions.write_text("SRR1\nSRR2\n")
    output.write_text("previous output\n")
    with pytest.raises(ValueError, match="fields not in fieldnames"):
        ama.extract_transcriptomic_rows(source, accessions, output)
    assert output.read_text() == "previous output\n"
    assert not list(tmp_path.glob("*.tmp"))


def test_streamed_merge_supports_primary_as_destination_with_schema_union(tmp_path):
    primary, extra = tmp_path / "primary.tsv", tmp_path / "extra.tsv"
    write_tsv(primary, ["run", "value"], [{"run": "SRR1", "value": "original"}])
    write_tsv(extra, ["run", "added"], [{"run": "SRR1", "added": "duplicate"},
                                     {"run": "SRR2", "added": "new"}])
    ama.merge_metadata_tables_by_run(primary, extra, primary)
    assert read_tsv(primary) == [{"run": "SRR1", "value": "original", "added": ""},
                                 {"run": "SRR2", "value": "", "added": "new"}]


def test_streamed_merge_detects_a_schema_change_before_publication(tmp_path, monkeypatch):
    source, extra, output = tmp_path / "source.tsv", tmp_path / "extra.tsv", tmp_path / "out.tsv"
    write_tsv(source, ["run", "value"], [{"run": "SRR1", "value": "old"}])
    output.write_text("previous output\n")
    original = ama.metadata_writer

    @contextmanager
    def changed_source(*args, **kwargs):
        write_tsv(source, ["run", "changed"], [{"run": "SRR2", "changed": "new"}])
        with original(*args, **kwargs) as writer:
            yield writer

    monkeypatch.setattr(ama, "metadata_writer", changed_source)
    with pytest.raises(ValueError, match="header changed"):
        ama.merge_metadata_tables_by_run(source, extra, output)
    assert output.read_text() == "previous output\n"


def test_sample_flags_are_updated_by_column_without_rewriting_description_text(tmp_path):
    source, output = tmp_path / "source.tsv", tmp_path / "out.tsv"
    fields = ["run", "description", "is_qualified", "scientific_name", "exclusion", "is_sampled"]
    description = "literal\t\t\tno\tdata"
    write_tsv(source, fields, [
        dict(zip(fields, ["SRR1", description, "", "Plant one", "no", ""], strict=True)),
        dict(zip(fields, ["SRR2", description, "", "Plant one", "excluded", ""], strict=True)),
        dict(zip(fields, ["SRR3", description, "no", "Plant one", "no", ""], strict=True)),
    ])
    ama.extract_species_rows(source, "Plant_one", output, qualify_unmarked=True)
    rows = read_tsv(output)
    assert rows[0]["is_sampled"] == rows[0]["is_qualified"] == "yes"
    assert rows[1]["is_sampled"] == rows[1]["is_qualified"] == ""
    assert rows[2]["is_sampled"] == "" and rows[2]["is_qualified"] == "no"
    assert all(row["description"] == description for row in rows)


def test_real_shell_species_filter_keeps_only_exact_species_rows(tmp_path):
    from shell_static_helpers import function_body
    support = Path(ama.__file__).parent
    core = support.parent / "core/gg_transcriptome_generation_core.sh"
    body = function_body(core.read_text(), "extract_sraid_metadata_rows_for_species")
    source, output = tmp_path / "source.tsv", tmp_path / "out.tsv"
    write_tsv(source, ["run", "scientific_name", "description"], [
        {"run": "SRR1", "scientific_name": "Plant one", "description": "literal\t\t\tno\tdata"},
        {"run": "SRR2", "scientific_name": "Other plant", "description": "Plant one"},
        {"run": "SRR3", "scientific_name": "Plant oneformis", "description": ""}])
    script = ('set -euo pipefail\ngg_support_dir=' + shlex.quote(str(support)) +
              '\n' + body + '\nextract_sraid_metadata_rows_for_species "$1" "$2" "$3"\n')
    result = subprocess.run(["bash", "-c", script, "test", str(source), "Plant one", str(output)],
                            text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    assert [r["run"] for r in read_tsv(output)] == ["SRR1"]
    assert read_tsv(output)[0]["description"] == "literal\t\t\tno\tdata"
