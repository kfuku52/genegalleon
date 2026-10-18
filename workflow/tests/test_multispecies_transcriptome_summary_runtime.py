import csv
import subprocess
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
SPECIES = "Osmunda_regalis"


def run_summary(tmp_path, busco_rows):
    inputs = {name: tmp_path / name for name in ("assembly", "metadata", "merge", "busco_isoform", "busco_cds")}
    for path in inputs.values():
        path.mkdir()
    source = inputs["busco_cds"] / f"{SPECIES}_busco.full.tsv"
    source.write_text("# Busco id\tStatus\tSequence\tScore\tLength\n" + busco_rows)
    (inputs["merge"] / f"{SPECIES}_tpm.tsv").write_text(
        "target_id\tERR2040894\n" + f"{SPECIES}_g2\t7.5\n{SPECIES}_g1\t2.5\n"
    )
    (inputs["assembly"] / f"{SPECIES}_assembly_stat.tsv").write_text(
        "file\tformat\ttype\tnum_seqs\tsum_len\tN50\n"
        f"{SPECIES}_isoform.fa.gz\tFASTA\tDNA\t3\t600\t224\n"
        f"{SPECIES}_longestCDS.fa.gz\tFASTA\tDNA\t3\t300\t126\n"
    )
    output = tmp_path / "summary"
    output.mkdir()
    # A previous expression table must not survive a new empty BUSCO result.
    (output / "expression.tsv").write_text("stale expression data\n")
    result = subprocess.run(
        [
            "Rscript", str(SUPPORT / "multispecies_transcriptome_summary.r"),
            f"--dir_assembly_stat={inputs['assembly']}",
            f"--dir_amalgkit_metadata={inputs['metadata']}",
            f"--dir_amalgkit_merge={inputs['merge']}",
            f"--dir_busco_isoform={inputs['busco_isoform']}",
            f"--dir_busco_longest_cds={inputs['busco_cds']}",
            f"--dir_myscript={SUPPORT}",
        ],
        cwd=output, text=True, capture_output=True, timeout=90,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert source.read_text().endswith(busco_rows)
    with (output / "assembly_stat_summary.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 2
    assert {row["N50"] for row in rows} == {"224", "126"}
    assert (output / "assembly_stat_summary.pdf").stat().st_size > 0
    assert "Ending multispecies_transcriptome_summary.r" in result.stdout
    return output, result


@pytest.mark.parametrize("busco_rows", ["88at2759\tMissing\n844at2759\tMissing\n", ""])
def test_empty_busco_preserves_assembly_summary_without_fabricated_expression(tmp_path, busco_rows):
    output, result = run_summary(tmp_path, busco_rows)
    assert "no detected BUSCO genes" in result.stdout
    with (output / "busco_table.tsv").open() as handle:
        assert list(csv.DictReader(handle, delimiter="\t")) == []
    assert not (output / "expression.tsv").exists()
    assert not (output / "expression.imputed.tsv").exists()
    assert not list(output.glob("multisp_*.pdf"))


@pytest.mark.parametrize("second_status", ["Complete", "Fragmented"])
def test_detected_busco_keeps_expression_values_and_busco_order(tmp_path, second_status):
    output, result = run_summary(
        tmp_path,
        f"88at2759\tComplete\t{SPECIES}_g1\t50\t100\n"
        f"844at2759\t{second_status}\t{SPECIES}_g2\t40\t90\n",
    )
    assert "no detected BUSCO genes" not in result.stdout
    with (output / "expression.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["busco_id"]: float(row["ERR2040894"]) for row in rows} == {
        "88at2759": 2.5, "844at2759": 7.5,
    }
    assert not list(output.glob("multisp_*.pdf"))
