import subprocess
import sys
from pathlib import Path

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT_DIR))

from score_hgt_candidates import BRANCH_OUTPUT_COLUMNS, GENE_OUTPUT_COLUMNS, ORTHOGROUP_OUTPUT_COLUMNS  # noqa: E402

SCRIPT_PATH = SUPPORT_DIR / "write_hgt_output_readme.py"


def write_header(path: Path, columns):
    path.write_text("\t".join(columns) + "\n", encoding="utf-8")


def test_hgt_output_readme_documents_all_hgt_table_columns(tmp_path: Path):
    branch_tsv = tmp_path / "hgt_branch_candidates.tsv"
    gene_tsv = tmp_path / "hgt_gene_candidates.tsv"
    orthogroup_tsv = tmp_path / "hgt_orthogroup_summary.tsv"
    readme_md = tmp_path / "README.md"
    write_header(branch_tsv, BRANCH_OUTPUT_COLUMNS)
    write_header(gene_tsv, GENE_OUTPUT_COLUMNS)
    write_header(orthogroup_tsv, ORTHOGROUP_OUTPUT_COLUMNS)

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--output",
            str(readme_md),
            "--branch_tsv",
            str(branch_tsv),
            "--gene_tsv",
            str(gene_tsv),
            "--orthogroup_tsv",
            str(orthogroup_tsv),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    readme_text = readme_md.read_text(encoding="utf-8")
    assert "# GeneGalleon HGT output tables" in readme_text
    assert "空欄は「陰性」ではなく「未測定・比較不能」" in readme_text
    assert "未定義" not in readme_text
    for column in BRANCH_OUTPUT_COLUMNS + GENE_OUTPUT_COLUMNS + ORTHOGROUP_OUTPUT_COLUMNS:
        assert f"| `{column}` |" in readme_text
