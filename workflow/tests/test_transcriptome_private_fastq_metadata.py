import csv
import re
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
CORE_SCRIPT = REPO_ROOT / "workflow" / "core" / "gg_transcriptome_generation_core.sh"
SUPPORT_DIR = REPO_ROOT / "workflow" / "support"


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


def _run_repair_helper(metadata_path: Path, species_label: str, output_path: Path):
    return subprocess.run(
        [sys.executable, "-", str(metadata_path), species_label, str(SUPPORT_DIR), str(output_path)],
        input=_embedded_python("repair_private_fastq_metadata_scientific_names"),
        text=True,
        capture_output=True,
        check=False,
    )


def _read_rows(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def test_repair_private_fastq_metadata_fills_placeholder_from_species_label(tmp_path: Path):
    metadata_path = tmp_path / "metadata_private_fastq.tsv"
    output_path = tmp_path / "metadata_private_fastq.fixed.tsv"
    metadata_path.write_text(
        "\n".join(
            [
                "run\tscientific_name\tlib_layout",
                "Cs1\tPlease add in format: Genus species\tpaired",
                "Cs2\tPlease add in format: Genus species\tpaired",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    completed = _run_repair_helper(metadata_path, "Chamaegastrodia_sikokiana", output_path)

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    rows = _read_rows(output_path)
    assert [row["scientific_name"] for row in rows] == [
        "Chamaegastrodia sikokiana",
        "Chamaegastrodia sikokiana",
    ]


def test_repair_private_fastq_metadata_canonicalizes_species_label_suffixes(tmp_path: Path):
    metadata_path = tmp_path / "metadata_private_fastq.tsv"
    output_path = tmp_path / "metadata_private_fastq.fixed.tsv"
    metadata_path.write_text(
        "\n".join(
            [
                "run\tscientific_name\tlib_layout",
                "Bu39-1\tPlease add in format: Genus species\tpaired",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    completed = _run_repair_helper(metadata_path, "Burmannia_cryptopetala-previous", output_path)

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    rows = _read_rows(output_path)
    assert rows[0]["scientific_name"] == "Burmannia cryptopetala"


def test_repair_private_fastq_metadata_preserves_existing_scientific_names(tmp_path: Path):
    metadata_path = tmp_path / "metadata_private_fastq.tsv"
    output_path = tmp_path / "metadata_private_fastq.fixed.tsv"
    metadata_path.write_text(
        "\n".join(
            [
                "run\tscientific_name\tlib_layout",
                "P1\tPetrosavia sakuraii\tpaired",
                "P2\tPlease add in format: Genus species\tpaired",
                "P3\t\tpaired",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    completed = _run_repair_helper(metadata_path, "Petrosavia_sakuraii-previous", output_path)

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    rows = _read_rows(output_path)
    assert [row["scientific_name"] for row in rows] == [
        "Petrosavia sakuraii",
        "Petrosavia sakuraii",
        "Petrosavia sakuraii",
    ]
