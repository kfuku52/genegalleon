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


def test_stage_quant_reference_fasta_aliases_uses_metadata_scientific_name_columns(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    reference_path = tmp_path / "Species_one_longestCDS.fa.gz"
    output_dir = tmp_path / "fasta"

    metadata_path.write_text(
        "\n".join(
            [
                "run\tscientific_name\tscientific_name_original",
                "SRR1\tSpecies one strain X\tSpecies one/original",
                "SRR2\tSpecies one strain X\tSpecies one/original",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    reference_path.write_text(">seq\nAAAA\n", encoding="utf-8")

    completed = subprocess.run(
        [sys.executable, "-", str(metadata_path), str(reference_path), str(output_dir), "Species_one"],
        input=_embedded_python("stage_quant_reference_fasta_aliases"),
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    canonical_alias = output_dir / "Species_one_for_kallisto_index.fasta"
    scientific_name_alias = output_dir / "Species_one_strain_X_for_kallisto_index.fasta"
    original_name_alias = output_dir / "Species_one_original_for_kallisto_index.fasta"

    for alias_path in (canonical_alias, scientific_name_alias, original_name_alias):
        assert alias_path.is_symlink(), f"Expected symlink alias to exist: {alias_path}"
        assert alias_path.resolve() == reference_path.resolve()


def test_stage_quant_reference_fasta_aliases_still_stages_canonical_prefix_without_metadata_columns(tmp_path):
    metadata_path = tmp_path / "metadata.tsv"
    reference_path = tmp_path / "Species_one_longestCDS.fa.gz"
    output_dir = tmp_path / "fasta"

    metadata_path.write_text(
        "\n".join(
            [
                "run\tlib_layout",
                "SRR1\tpaired",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    reference_path.write_text(">seq\nAAAA\n", encoding="utf-8")

    completed = subprocess.run(
        [sys.executable, "-", str(metadata_path), str(reference_path), str(output_dir), "Species_one"],
        input=_embedded_python("stage_quant_reference_fasta_aliases"),
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    canonical_alias = output_dir / "Species_one_for_kallisto_index.fasta"
    assert canonical_alias.is_symlink()
    assert canonical_alias.resolve() == reference_path.resolve()
