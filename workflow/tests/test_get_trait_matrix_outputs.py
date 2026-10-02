import subprocess
import sys
from pathlib import Path

import pytest

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


@pytest.mark.parametrize("ncpu", [1, 2])
@pytest.mark.parametrize("qualified", [False, True])
def test_get_trait_matrix_keeps_species_with_colliding_gene_ids(tmp_path, ncpu, qualified):
    traits = tmp_path / "traits"
    traits.mkdir()
    for species, value in (("Species_a", 11), ("Species_b", 99)):
        gene = species + "_gene1" if qualified else "gene1"
        # A fully qualified ID from another species must not match this file.
        foreign = "Species_b_gene2\t123\n" if species == "Species_a" else ""
        (traits / (species + ".tsv")).write_text(f"gene_id\troot\n{gene}\t{value}\n{foreign}")
    fasta = tmp_path / "genes.fa"
    fasta.write_text(">Species_a_gene1\nATG\n>Species_b_gene1\nATG\n>Species_b_gene2\nATG\n")
    output = tmp_path / "expression.tsv"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--dir_trait", str(traits), "--seqfile", str(fasta),
         "--outfile", str(output), "--ncpu", str(ncpu)], capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert output.read_text() == "gene_id\troot\nSpecies_a_gene1\t11\nSpecies_b_gene1\t99\n"


@pytest.mark.parametrize('ncpu', [1, 2])
@pytest.mark.parametrize('qualified', [False, True])
def test_trait_identifiers_keep_literal_na_words_and_leading_zeroes(tmp_path, ncpu, qualified):
    import pandas as pd

    genes = ['NA', '001', 'NULL', 'nan', '000']
    traits = tmp_path / 'traits'
    traits.mkdir()
    expected_ids, expected_values = [], []
    for species, offset in [('Species_a', 0), ('Species_b', 100)]:
        lines = ['gene_id\tmeasurement\tmissing_value']
        for i, gene in enumerate(genes):
            name = species + '_' + gene
            key = name if qualified else gene
            lines.append(f'{key}\t{offset + i + 1}\tNA')
            expected_ids.append(name)
            expected_values.append(offset + i + 1)
        lines.append('unmatched\t999\tNULL')
        (traits / (species + '.tsv')).write_text('\n'.join(lines) + '\n')
    fasta = tmp_path / 'genes.fa'
    fasta.write_text(''.join(f'>{name}\nATG\n' for name in expected_ids))
    output = tmp_path / 'traits.tsv'
    result = subprocess.run([sys.executable, str(SCRIPT), '--dir_trait', str(traits),
        '--seqfile', str(fasta), '--outfile', str(output), '--ncpu', str(ncpu)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    actual = pd.read_csv(output, sep='\t')
    expected = pd.DataFrame({'gene_id': expected_ids, 'measurement': expected_values,
                             'missing_value': [float('nan')] * len(expected_ids)})
    pd.testing.assert_frame_equal(actual, expected)


@pytest.mark.parametrize('inferred_index', [False, True])
def test_numeric_gene_identifiers_keep_zero_padding_and_trait_types(tmp_path, inferred_index):
    from workflow.support.get_trait_matrix import process_trait_file

    path = tmp_path / 'traits.tsv'
    header = 'measurement\tmissing\n' if inferred_index else 'gene_id\tmeasurement\tmissing\n'
    path.write_text(header + '001\t1\tNA\n000\t2\tNULL\n0\t3\tnan\n')
    mapping = {gene: 'Species_a_' + gene for gene in ['001', '000', '0']}
    result = process_trait_file(path, set(mapping), mapping)
    assert result.gene_id.tolist() == list(mapping.values())
    assert result.measurement.tolist() == [1, 2, 3]
    assert str(result.measurement.dtype) == 'int64'
    assert result['missing'].isna().all()
