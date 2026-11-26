import subprocess
import sys
from pathlib import Path

import pytest

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))

import transcript_quantification_qc as qc  # noqa: E402


@pytest.fixture
def merge_dir(tmp_path):
    root = tmp_path / "amalgkit_merge"
    species = root / "Acorus_calamus_var._americanus"
    species.mkdir(parents=True)
    (species / f"{species.name}_metadata.tsv").write_text(
        "run\tscientific_name\tsample_group\tlib_layout\texclusion\tmapping_rate\t"
        "fastp_duplication_rate\tfastp_insert_size_peak\ttotal_spots\ttotal_bases\n"
        "SRR1\tAcorus calamus var. americanus\tflower\tpaired\tno\t91\t24\t200\t1000000\t130000000\n"
        "SRR2\tAcorus calamus var. americanus\tleaf\tsingle\tno\t90\t\t\t1500000\t140000000\n"
        "SRR3\tAcorus calamus var. americanus\troot\tpaired\treated\t92\t22\t220\t1600000\t150000000\n",
        encoding="utf-8",
    )
    (species / f"{species.name}_tpm.tsv").write_text(
        "target_id\tSRR2\tSRR3\tSRR1\n"
        "tx1\t1000000\t0\t700000\n"
        "tx2\t0\t500000\t0\n"
        "tx3\t0\t500000\t300000\n",
        encoding="utf-8",
    )
    return root


def test_collect_and_reuse_without_opening_tpm_or_unrelated_species(merge_dir, tmp_path):
    pytest.importorskip("pandas")
    pytest.importorskip("matplotlib")
    unrelated = merge_dir / "Other_species"
    unrelated.mkdir()
    (unrelated / "Other_species_metadata.tsv").write_text(
        "run\tscientific_name\nSRR4\tOther species\n", encoding="utf-8",
    )
    retained = tmp_path / "retained.tsv"
    retained.write_text("run\nSRR1\n", encoding="utf-8")
    summary = tmp_path / "summary.tsv"
    source_plot, reused_plot = tmp_path / "from_merge.svg", tmp_path / "from_summary.svg"
    command = [sys.executable, str(SUPPORT_DIR / "transcript_quantification_qc.py"), "--retained-runs", str(retained)]
    subprocess.run(
        [*command, "--merge-dir", str(merge_dir), "--write-summary", str(summary), "--output", str(source_plot)],
        check=True,
    )
    rows = qc.read_table(summary, qc.SUMMARY_COLUMNS)
    assert {row["run"]: float(row["mean_tpm"]) for row in rows} == {"SRR2": 1e6, "SRR1": 5e5}
    for path in merge_dir.glob("*/*_tpm.tsv"):
        path.unlink()
    subprocess.run([*command, "--summary", str(summary), "--output", str(reused_plot)], check=True)
    for plot in (source_plot, reused_plot):
        svg = plot.read_text()
        assert 'width="518.4pt"' in svg
        assert "Mean TPM uses transcripts with TPM &gt; 0 in 2 original runs." in svg
        assert "Other columns use 1 retained runs." in svg


def test_axes_jitter_and_taxonomic_rank(merge_dir):
    pytest.importorskip("pandas")
    pytest.importorskip("matplotlib")
    originals = qc.collect_merge_qc(merge_dir)
    fig = qc.draw_qc(originals, originals)
    assert fig.get_figwidth() == pytest.approx(7.2)
    mapping, duplication, insert = fig.axes[1:4]
    assert mapping.get_xlim() == duplication.get_xlim() == (0, 100)
    assert insert.get_xlim() == (0, 300)
    assert not any(axis.patches for axis in fig.axes[1:6])
    points = mapping.collections[-1].get_offsets()
    assert points[:, 0].tolist() == [90, 91]
    assert points[0, 1] - points[1, 1] > 2
    labels = {text.get_text().strip(): text.get_fontstyle() for text in fig.texts}
    assert labels["Acorus calamus"] == labels["americanus"] == "italic"
    assert labels["var."] == "normal"
    assert lines_without_fliers(fig.axes[-1])


def lines_without_fliers(axis):
    return bool(axis.patches) and all(line.get_marker() in ("None", "", None) for line in axis.lines)


def test_tpm_column_order_and_missing_values(merge_dir):
    pytest.importorskip("pandas")
    path = next(merge_dir.glob("*/*_tpm.tsv"))
    assert qc.mean_detected_tpm(path, ["SRR1", "SRR2"]) == {"SRR1": 5e5, "SRR2": 1e6}
    path.write_text("target_id\tSRR1\ntranscript\tNA\n", encoding="utf-8")
    with pytest.raises(ValueError, match="finite"):
        qc.mean_detected_tpm(path, ["SRR1"])
