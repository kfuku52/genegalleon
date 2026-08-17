import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "support" / "get_trait_matrix.py"


def run_get_trait_matrix(tmp_path: Path, trait_text: str | None):
    trait_dir = tmp_path / "traits"
    trait_dir.mkdir()
    if trait_text is not None:
        (trait_dir / "Species_a.tsv").write_text(trait_text, encoding="utf-8")

    seqfile = tmp_path / "genes.fa"
    seqfile.write_text(">Species_a_gene1\nATG\n", encoding="utf-8")
    outfile = tmp_path / "expression.tsv"
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--dir_trait",
            str(trait_dir),
            "--seqfile",
            str(seqfile),
            "--outfile",
            str(outfile),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    return completed, outfile


def test_get_trait_matrix_omits_output_without_usable_traits(tmp_path):
    for case_name, trait_text in (
        ("no_files", None),
        ("identifier_only", "gene_id\ngene1\n"),
    ):
        case_dir = tmp_path / case_name
        case_dir.mkdir()
        completed, outfile = run_get_trait_matrix(case_dir, trait_text)

        assert completed.returncode == 0, completed.stderr
        assert not outfile.exists()
        assert "No usable trait matrix was generated" in completed.stdout


def test_get_trait_matrix_writes_matched_trait_values(tmp_path):
    completed, outfile = run_get_trait_matrix(
        tmp_path,
        "gene_id\troot_1\troot_2\ngene1\t1.0\t2.0\n",
    )

    assert completed.returncode == 0, completed.stderr
    assert outfile.read_text(encoding="utf-8") == (
        "gene_id\troot_1\troot_2\n"
        "Species_a_gene1\t1.0\t2.0\n"
    )
