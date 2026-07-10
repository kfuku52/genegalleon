import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow" / "support" / "prepare_orthogroup_copy_number.py"
GENOME_EVOLUTION_CORE = REPO_ROOT / "workflow" / "core" / "gg_genome_evolution_core.sh"


def _run_prepare(tmp_path: Path, genecount_text: str, tree_text: str, *extra_args: str) -> subprocess.CompletedProcess[str]:
    genecount = tmp_path / "Orthogroups.GeneCount.selected.tsv"
    tree = tmp_path / "dated_species_tree.nwk"
    outdir = tmp_path / "orthogroup_copy_number"
    genecount.write_text(genecount_text, encoding="utf-8")
    tree.write_text(tree_text, encoding="utf-8")
    env = os.environ.copy()
    mpl_config = tmp_path / "mplconfig"
    mpl_config.mkdir(parents=True, exist_ok=True)
    env["MPLCONFIGDIR"] = str(mpl_config)
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--genecount",
            str(genecount),
            "--dated_species_tree",
            str(tree),
            "--output_dir",
            str(outdir),
            *extra_args,
        ],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )


def test_prepare_orthogroup_copy_number_filters_and_writes_matrix(tmp_path):
    proc = _run_prepare(
        tmp_path,
        "\n".join(
            [
                "besthit_0.95\tOrthogroup\tsp1\tsp2\tsp3\tsp4",
                "hit1\tOG1\t1\t2\t3\t4",
                "hit2\tOG_TOO_WIDE\t1\t1\t1\t99",
                "",
            ]
        ),
        "((sp1:1,sp2:1):1,(sp3:1,sp4:1):1);\n",
        "--max_size_differential",
        "10",
    )
    assert proc.returncode == 0, proc.stderr

    outdir = tmp_path / "orthogroup_copy_number"
    matrix = outdir / "orthogroup_copy_number.tsv"
    removed = outdir / "removed_orthogroups.tsv"
    assert matrix.exists()
    assert removed.exists()
    assert "OG1" in matrix.read_text(encoding="utf-8")
    assert "OG_TOO_WIDE" not in matrix.read_text(encoding="utf-8")
    assert "OG_TOO_WIDE" in removed.read_text(encoding="utf-8")


def test_prepare_orthogroup_copy_number_reports_missing_tree_species(tmp_path):
    proc = _run_prepare(
        tmp_path,
        "\n".join(
            [
                "besthit_0.95\tOrthogroup\tsp1\tsp2\tsp3",
                "hit1\tOG1\t1\t2\t3",
                "",
            ]
        ),
        "((sp1:1,sp2:1):1,(sp3:1,sp4:1):1);\n",
    )
    assert proc.returncode != 0
    assert "missing species column(s)" in proc.stderr
    assert "sp4" in proc.stderr


def test_prepare_orthogroup_copy_number_reports_nonnumeric_counts(tmp_path):
    proc = _run_prepare(
        tmp_path,
        "\n".join(
            [
                "besthit_0.95\tOrthogroup\tsp1\tsp2\tsp3\tsp4",
                "hit1\tOG1\t1\tbad\t3\t4",
                "",
            ]
        ),
        "((sp1:1,sp2:1):1,(sp3:1,sp4:1):1);\n",
    )
    assert proc.returncode != 0
    assert "non-numeric copy numbers" in proc.stderr
    assert "sp2" in proc.stderr


def test_genome_evolution_wires_trait_pgls_without_requiring_cafe():
    core = GENOME_EVOLUTION_CORE.read_text(encoding="utf-8")
    prep_condition = (
        'if [[ ! -s "${file_orthogroup_copy_number}" && '
        '( ${run_cafe} -eq 1 || ${run_orthogroup_copy_number_trait_pgls} -eq 1 ) ]]'
    )
    trait_input_check = (
        'disable_if_no_input_file "run_orthogroup_copy_number_trait_pgls" '
        '"${file_orthogroup_genecount_selected}" "${file_dated_species_tree}" "${file_trait}"'
    )
    assert prep_condition in core
    assert 'task="Orthogroup copy-number matrix preparation"' in core
    assert trait_input_check in core
    assert 'cafe5 \\' in core
    assert '--infile "${file_orthogroup_copy_number}"' in core
