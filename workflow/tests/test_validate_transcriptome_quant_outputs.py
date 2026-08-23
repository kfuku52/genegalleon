from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow" / "support" / "validate_transcriptome_quant_outputs.py"
CORE = REPO_ROOT / "workflow" / "core" / "gg_transcriptome_generation_core.sh"
ABUNDANCE_HEADER = "target_id\tlength\teff_length\test_counts\ttpm\n"


def write_metadata(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("run", "scientific_name", "exclusion", "is_sampled"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_run(quant_root: Path, run_id: str, *, abundance: str | None = None) -> None:
    run_dir = quant_root / run_id
    run_dir.mkdir(parents=True)
    (run_dir / f"{run_id}_abundance.tsv").write_text(
        abundance or ABUNDANCE_HEADER + "tx1\t100\t80\t4\t2.5\n",
        encoding="utf-8",
    )
    (run_dir / f"{run_id}_run_info.json").write_text(
        json.dumps({"p_pseudoaligned": 92.5}) + "\n",
        encoding="utf-8",
    )


def run_validator(metadata: Path, quant_root: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--metadata",
            str(metadata),
            "--quant-root",
            str(quant_root),
        ],
        text=True,
        capture_output=True,
        check=False,
    )


def test_empty_quant_directory_is_not_a_complete_legacy_artifact(tmp_path: Path):
    metadata = tmp_path / "metadata.tsv"
    quant_root = tmp_path / "quant"
    quant_root.mkdir()
    write_metadata(
        metadata,
        [{"run": "SRR1", "scientific_name": "Test species", "exclusion": "no", "is_sampled": "yes"}],
    )

    completed = run_validator(metadata, quant_root)

    assert completed.returncode == 1
    assert "quant output directory is missing for run: SRR1" in completed.stderr


def test_every_selected_run_requires_content_valid_abundance_and_run_info(tmp_path: Path):
    metadata = tmp_path / "metadata.tsv"
    quant_root = tmp_path / "quant"
    quant_root.mkdir()
    write_metadata(
        metadata,
        [
            {"run": "SRR1", "scientific_name": "Test species", "exclusion": "no", "is_sampled": "yes"},
            {"run": "SRR2", "scientific_name": "Test species", "exclusion": "no", "is_sampled": "yes"},
        ],
    )
    write_run(quant_root, "SRR1")

    missing = run_validator(metadata, quant_root)
    assert missing.returncode == 1
    assert "quant output directory is missing for run: SRR2" in missing.stderr

    write_run(
        quant_root,
        "SRR2",
        abundance=ABUNDANCE_HEADER + "tx1\t100\t80\t-1\t2.5\n",
    )
    invalid = run_validator(metadata, quant_root)
    assert invalid.returncode == 1
    assert "abundance table for SRR2 has invalid est_counts" in invalid.stderr


def test_complete_selected_run_set_passes_and_ignores_unselected_rows(tmp_path: Path):
    metadata = tmp_path / "metadata.tsv"
    quant_root = tmp_path / "quant"
    quant_root.mkdir()
    write_metadata(
        metadata,
        [
            {"run": "SRR1", "scientific_name": "Test species", "exclusion": "no", "is_sampled": "yes"},
            {"run": "SRR2", "scientific_name": "Test species", "exclusion": "excluded", "is_sampled": "no"},
        ],
    )
    write_run(quant_root, "SRR1")

    completed = run_validator(metadata, quant_root)

    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout) == {"run_count": 1, "runs": ["SRR1"], "status": "valid"}


def test_transcriptome_quant_stage_validates_before_reuse_and_before_publish():
    text = CORE.read_text(encoding="utf-8")
    validator = '"${gg_support_dir}/validate_transcriptome_quant_outputs.py"'

    assert text.count(validator) == 2
    assert text.index(validator) < text.index(
        "gg_artifact_prepare_stage quant_needs_update run_amalgkit_quant"
    )
    assert text.index(validator, text.index(validator) + 1) < text.index(
        'mv_out_replace_dir "./quant" "${dir_amalgkit_quant}/${sp_ub}"'
    )
    assert "quant_outputs=(./quant/*)" not in text
