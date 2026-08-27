import csv
import subprocess
import sys
from pathlib import Path

from shell_static_helpers import REPO_ROOT

SELECT_CORE_SPECIES = REPO_ROOT / "workflow" / "support" / "select_orthofinder_core_species.py"
SOURCE_REVISIONS = Path("/opt/pg/logs/source_revisions.tsv")
EXPECTED_SOURCES = {
    "amalgkit",
    "cdskit",
    "csubst",
    "nwkit",
    "BUSCO",
    "paml",
    "kfl1ou",
    "kfFractBias",
    "kftools",
    "rkftools",
    "RADTE",
}


def _write_protein(path: Path, count: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(f">{path.stem}_gene{index}\nMPEP\n" for index in range(1, count + 1)),
        encoding="utf-8",
    )


def _write_busco(path: Path, complete_pct: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"C:{complete_pct}%[S:{complete_pct}%,D:0.0%],F:0.0%,M:0.0%,n:100\n",
        encoding="utf-8",
    )


def test_runtime_records_exact_revisions_for_every_moving_source():
    assert SOURCE_REVISIONS.is_file(), f"Missing container source manifest: {SOURCE_REVISIONS}"
    with SOURCE_REVISIONS.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    revisions = {row["source"]: row["revision"] for row in rows}
    assert set(revisions) == EXPECTED_SOURCES
    assert all(len(revision) == 40 for revision in revisions.values())
    assert all(set(revision) <= set("0123456789abcdef") for revision in revisions.values())


def test_genegalleon_core_species_selection_executes_real_nwkit_sample(tmp_path: Path):
    for species, protein_count, busco_pct in (
        ("Alpha_species", 3, 95.0),
        ("Beta_species", 2, 90.0),
        ("Gamma_species", 1, 85.0),
    ):
        _write_protein(tmp_path / "protein" / f"{species}.fa", protein_count)
        _write_busco(tmp_path / "busco" / f"{species}.busco.short.txt", busco_pct)

    tree = tmp_path / "species.nwk"
    tree.write_text("(Alpha_species:1,Beta_species:1,Gamma_species:1);\n", encoding="utf-8")
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    command = [
        sys.executable,
        str(SELECT_CORE_SPECIES),
        "--protein-dir",
        str(tmp_path / "protein"),
        "--busco-short-dir",
        str(tmp_path / "busco"),
        "--species-tree",
        str(tree),
        "--max-core-species",
        "2",
        "--filters",
        "busco_complete_pct:ge:80,num_seq:le:100000",
        "--rank",
        "num_seq:asc,busco_complete_pct:desc",
        "--method",
        "max-pd",
        "--candidates-table",
        str(output_dir / "candidates.tsv"),
        "--selected-table",
        str(output_dir / "selected.tsv"),
        "--selected-list",
        str(output_dir / "selected.txt"),
        "--core-tree",
        str(output_dir / "core.nwk"),
    ]

    completed = subprocess.run(command, capture_output=True, text=True, check=False)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    selected = (output_dir / "selected.txt").read_text(encoding="utf-8").splitlines()
    assert len(selected) == 2
    assert set(selected) <= {"Alpha_species.fa", "Beta_species.fa", "Gamma_species.fa"}
    report_header = (output_dir / "selected.tsv").read_text(encoding="utf-8").splitlines()[0]
    assert report_header.startswith("sample_order\tleaf_name\t")
    assert (output_dir / "core.nwk").is_file()
