from __future__ import annotations

import csv
import json
import subprocess
import sys
import zipfile
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


def test_reordered_abundance_columns_keep_quoted_ids_and_check_short_rows(tmp_path: Path):
    metadata = tmp_path / "metadata.tsv"
    quant_root = tmp_path / "quant"
    write_metadata(metadata, [{"run": "SRR1", "scientific_name": "Test species"}])
    header = "tpm\ttarget_id\test_counts\teff_length\tlength\n\n"
    row = '2.5\t"tx\t1"\t4\t80'
    write_run(quant_root, "SRR1", abundance=header + row + "\t100\n")
    completed = run_validator(metadata, quant_root)
    assert completed.returncode == 0, completed.stderr

    (quant_root / "SRR1" / "SRR1_abundance.tsv").write_text(header + row + "\n")
    truncated = run_validator(metadata, quant_root)
    assert truncated.returncode == 1
    assert "abundance table for SRR1 has invalid length" in truncated.stderr


def test_transcriptome_quant_stage_validates_before_reuse_and_before_publish():
    text = CORE.read_text(encoding="utf-8")
    validator = '"${gg_support_dir}/validate_transcriptome_quant_outputs.py"'

    assert text.count(validator) == 2
    assert text.index(validator) < text.index(
        "gg_artifact_prepare_stage quant_needs_update run_amalgkit_quant"
    )
    assert text.index(validator, text.index(validator) + 1) < text.index(
        'pack --raw "${quant_species_output_dir}" --source "./quant"'
    )
    assert "quant_outputs=(./quant/*)" not in text


def test_zip_migration_preserves_identity_and_only_stages_selected_merge_inputs(tmp_path: Path, monkeypatch):
    monkeypatch.syspath_prepend(str(SCRIPT.parent))
    from archive_transcriptome_quant import pack_quant, stage_merge
    from artifact_provenance import sha256_path
    from validate_transcriptome_quant_outputs import validate_quant_outputs

    raw = tmp_path / "output" / "transcriptome_quantification" / "amalgkit_quant" / "Test_species"
    metadata = tmp_path / "metadata.tsv"
    write_metadata(metadata, [
        {"run": "R1", "scientific_name": "Test species"},
        {"run": "R1-old", "scientific_name": "Test species"},
    ])
    for run in ("R1", "R1-old", "excluded"):
        write_run(raw, run)
        (raw / run / f"{run}_abundance.h5").write_bytes(b"binary-hdf5-" + run.encode())
    (raw / "transcript_gene.tsv").write_text("tx1\tg1\n")
    expected = sha256_path(raw)
    manifest = tmp_path / "Test_species.quant.json"
    payload = {
        "schema_version": 1, "step": "transcriptome_quant", "family_id": raw.name,
        "outputs": [{
            "artifact_type": "directory", "label": "quant_directory", "path": str(raw),
            "sha256": expected[0], "size_bytes": expected[1],
        }],
    }
    manifest.write_text(json.dumps(payload))
    archive = pack_quant(raw, raw, manifest)
    assert not raw.exists()
    assert sha256_path(raw) == expected
    assert validate_quant_outputs(metadata, raw, require_transcript_map=True)["run_count"] == 2

    destination = tmp_path / "quant"
    stage_merge(raw, metadata, destination)
    assert {path.relative_to(destination).as_posix() for path in destination.rglob("*") if path.is_file()} == {
        f"{run}/{run}_{suffix}"
        for run in ("R1", "R1-old") for suffix in ("abundance.tsv", "run_info.json")
    }
    with zipfile.ZipFile(archive) as reader:
        for path in destination.rglob("*"):
            if path.is_file():
                assert path.read_bytes() == reader.read(raw.name + "/" + path.relative_to(destination).as_posix())
        assert len(reader.namelist()) == 10


def test_pack_refuses_unverified_migration_and_validator_rejects_unsafe_zip(tmp_path: Path, monkeypatch):
    monkeypatch.syspath_prepend(str(SCRIPT.parent))
    import pytest
    from archive_transcriptome_quant import pack_quant

    raw = tmp_path / "output" / "transcriptome_quantification" / "amalgkit_quant" / "Test_species"
    write_run(raw, "R1")
    manifest = tmp_path / "Test_species.quant.json"
    manifest.write_text(json.dumps({
        "schema_version": 1, "step": "transcriptome_quant", "family_id": raw.name,
        "outputs": [{
            "artifact_type": "directory", "label": "quant_directory", "path": str(raw),
            "sha256": "0" * 64, "size_bytes": 0,
        }],
    }))
    archive = raw.with_suffix(".zip")
    with pytest.raises(ValueError, match="differs from recorded"):
        pack_quant(raw, raw, manifest)
    assert raw.is_dir()
    assert not archive.exists()

    raw.rename(tmp_path / "retained")
    metadata = tmp_path / "selected.tsv"
    write_metadata(metadata, [{"run": "R1", "scientific_name": "Test species"}])
    with zipfile.ZipFile(archive, "w") as writer:
        writer.writestr(f"{raw.name}/../outside", "bad")
    result = run_validator(metadata, raw)
    assert result.returncode == 1
    assert "Unsafe or unexpected ZIP member" in result.stderr
