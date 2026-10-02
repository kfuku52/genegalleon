import csv
import os
import subprocess
import sys
from pathlib import Path

import pytest
from shell_static_helpers import read_text

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


@pytest.mark.parametrize('identifiers', [['001', '1'], ['NA', 'NULL', 'nan', '001', '1']])
def test_copy_number_preparation_retains_literal_orthogroup_identifiers(tmp_path, identifiers):
    source = 'besthit_0.95\tOrthogroup\tsp1\tsp2\tsp3\tsp4\n'
    source += ''.join(f'hit\t{identifier}\t1\t2\t3\t4\n' for identifier in identifiers)
    proc = _run_prepare(tmp_path, source, '((sp1:1,sp2:1):1,(sp3:1,sp4:1):1);\n')
    assert proc.returncode == 0, proc.stderr
    output = tmp_path / 'orthogroup_copy_number/orthogroup_copy_number.tsv'
    with output.open(newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert [row['Orthogroup'] for row in rows] == identifiers
    assert all([row[column] for column in ['sp1', 'sp2', 'sp3', 'sp4']] == ['1', '2', '3', '4'] for row in rows)


@pytest.mark.parametrize('compressed', [False, True])
def test_copy_number_identifier_reader_keeps_numeric_types_and_metadata_missingness(tmp_path, compressed):
    import gzip

    import pandas

    from workflow.support.prepare_orthogroup_copy_number import load_gene_count_table

    path = tmp_path / ('counts.tsv.gz' if compressed else 'counts.tsv')
    source = 'besthit_0.95\tOrthogroup\tsp1\tsp2\nNA\t001\t1\t2.5\nhit\t1\t3\t4.5\n'
    if compressed:
        with gzip.open(path, 'wt') as handle:
            handle.write(source)
    else:
        path.write_text(source)
    original = path.read_bytes()
    frame = load_gene_count_table(str(path))
    assert frame.Orthogroup.tolist() == ['001', '1']
    assert frame.sp1.dtype == 'int64' and frame.sp2.dtype == 'float64'
    assert frame.sp1.tolist() == [1, 3] and frame.sp2.tolist() == [2.5, 4.5]
    assert pandas.isna(frame.at[0, 'besthit_0.95']) and frame.at[1, 'besthit_0.95'] == 'hit'
    assert path.read_bytes() == original


@pytest.mark.parametrize('empty_id', ['', ' '])
def test_copy_number_identifier_reader_still_rejects_empty_ids(tmp_path, empty_id):
    from workflow.support.prepare_orthogroup_copy_number import load_gene_count_table

    path = tmp_path / 'counts.tsv'
    path.write_text(f'besthit_0.95\tOrthogroup\tsp1\nhit\t{empty_id}\t1\n')
    with pytest.raises(SystemExit, match='contains empty Orthogroup IDs'):
        load_gene_count_table(str(path))


def test_copy_number_identifier_reader_still_rejects_actual_duplicate_ids(tmp_path):
    from workflow.support.prepare_orthogroup_copy_number import load_gene_count_table

    path = tmp_path / 'counts.tsv'
    path.write_text('besthit_0.95\tOrthogroup\tsp1\nhit\t001\t1\nhit\t001\t2\n')
    with pytest.raises(SystemExit, match='contains duplicate Orthogroup IDs: 001'):
        load_gene_count_table(str(path))


@pytest.mark.parametrize('count, message', [('NA', 'non-numeric copy numbers'), ('-1', 'negative copy numbers')])
def test_copy_number_identifier_fix_keeps_copy_number_validation(tmp_path, count, message):
    from workflow.support.prepare_orthogroup_copy_number import get_species_count_table, load_gene_count_table

    path = tmp_path / 'counts.tsv'
    path.write_text(f'besthit_0.95\tOrthogroup\tsp1\nhit\tNA\t{count}\n')
    frame = load_gene_count_table(str(path))
    assert frame.Orthogroup.tolist() == ['NA']
    with pytest.raises(SystemExit, match=message):
        get_species_count_table(frame, ['sp1'])


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


def test_duplicate_species_diagnostics_precede_missing_count_columns(tmp_path):
    proc = _run_prepare(tmp_path, 'besthit_0.95\tOrthogroup\tother\nhit\tOG1\t1\n',
                        '(z:1,a:1,b:1,z:1,a:1);\n')
    assert proc.returncode != 0
    assert proc.stderr.strip() == 'ERROR: Dated species tree has duplicate leaf label(s): a, z'


def test_genome_evolution_wires_trait_pgls_without_requiring_cafe():
    core = read_text(GENOME_EVOLUTION_CORE)
    prep_condition = (
        'if [[ ${copy_number_needs_update} -eq 1 && '
        '${run_orthogroup_copy_number_stage} -eq 1 ]]'
    )
    trait_input_check = (
        'disable_if_no_input_file "run_orthogroup_copy_number_trait_pgls" '
        '"${file_orthogroup_genecount_selected}" "${file_dated_species_tree}" "${file_trait}"'
    )
    assert prep_condition in core
    assert 'task="Orthogroup copy-number matrix preparation"' in core
    assert 'gg_artifact_prepare_stage copy_number_needs_update run_orthogroup_copy_number_stage' in core
    assert '--parameter "max_size_differential=${orthogroup_copy_number_max_size_differential}"' in core
    assert trait_input_check in core
    assert 'cafe5 \\' in core
    assert '--infile "${file_orthogroup_copy_number}"' in core
