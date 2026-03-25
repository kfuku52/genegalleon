import csv
import re
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
CORE_SCRIPT = REPO_ROOT / "workflow" / "core" / "gg_transcriptome_generation_core.sh"


def _function_body(text: str, function_name: str) -> str:
    pattern = re.compile(rf"^\s*{re.escape(function_name)}\(\)\s*\{{", re.MULTILINE)
    match = pattern.search(text)
    if match is None:
        raise AssertionError(f"Function not found: {function_name}")
    start = match.start()
    next_match = re.search(r"^\s*[A-Za-z_][A-Za-z0-9_]*\(\)\s*\{", text[match.end():], re.MULTILINE)
    if next_match is None:
        return text[start:]
    return text[start:match.end() + next_match.start()]


def _embedded_python(function_name: str) -> str:
    text = CORE_SCRIPT.read_text(encoding="utf-8")
    body = _function_body(text, function_name)
    match = re.search(r"<<'PY'\n(.*?)\nPY", body, re.DOTALL)
    if match is None:
        raise AssertionError(f"Embedded Python not found in function: {function_name}")
    return match.group(1)


def test_merge_metadata_tables_by_run_appends_new_columns_from_relaxed_metadata(tmp_path):
    primary_path = tmp_path / "metadata.tsv"
    extra_path = tmp_path / "metadata.relaxed.tsv"
    output_path = tmp_path / "metadata.merged.tsv"

    primary_path.write_text(
        "\n".join(
            [
                "run\tscientific_name\tlib_layout",
                "SRR1\tSpecies one\tsingle",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    extra_path.write_text(
        "\n".join(
            [
                "run\tscientific_name\tlib_layout\tsample_attribute_tissue\tdev_stage",
                "SRR1\tSpecies one\tsingle\tleaf\tadult",
                "SRR2\tSpecies two\tpaired\troot\tjuvenile",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    completed = subprocess.run(
        [sys.executable, "-", str(primary_path), str(extra_path), str(output_path)],
        input=_embedded_python("merge_metadata_tables_by_run"),
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    with output_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    assert reader.fieldnames == [
        "run",
        "scientific_name",
        "lib_layout",
        "sample_attribute_tissue",
        "dev_stage",
    ]
    assert [row["run"] for row in rows] == ["SRR1", "SRR2"]
    assert rows[0]["scientific_name"] == "Species one"
    assert rows[0]["sample_attribute_tissue"] == ""
    assert rows[0]["dev_stage"] == ""
    assert rows[1]["scientific_name"] == "Species two"
    assert rows[1]["sample_attribute_tissue"] == "root"
    assert rows[1]["dev_stage"] == "juvenile"
