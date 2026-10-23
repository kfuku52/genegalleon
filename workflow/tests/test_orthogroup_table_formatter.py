from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import pandas
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "orthogroup_table_formatter.py"


def load_module():
    spec = spec_from_file_location("orthogroup_table_formatter", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_hog2og_long_cells_do_not_require_fixed_width_numpy_strings(tmp_path, monkeypatch):
    mod = load_module()
    input_path = tmp_path / "N1.tsv"
    out_dir = tmp_path / "Orthogroups"
    unrelated_tmp = tmp_path / "tmp.user.tsv"
    huge_cell = ", ".join(f"gene{i}" for i in range(6000))
    unrelated_tmp.write_text("keep\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    pandas.DataFrame(
        [
            {
                "OG": "OG0000001",
                "HOG": "N1.HOG0000001",
                "Gene Tree Parent Clade": "N1",
                "spA": huge_cell,
                "spB": "beta1, beta2",
            },
            {
                "OG": "OG0000002",
                "HOG": "N1.HOG0000002",
                "Gene Tree Parent Clade": "N1",
                "spA": "",
                "spB": "beta3",
            },
        ]
    ).to_csv(input_path, sep="\t", index=False)

    original_to_numpy = pandas.DataFrame.to_numpy

    def guarded_to_numpy(self, *args, **kwargs):
        dtype = kwargs.get("dtype")
        if dtype is None and args:
            dtype = args[0]
        if dtype is str:
            raise AssertionError("DataFrame.to_numpy(dtype=str) should not be used here")
        return original_to_numpy(self, *args, **kwargs)

    monkeypatch.setattr(pandas.DataFrame, "to_numpy", guarded_to_numpy)

    mod.run(
        SimpleNamespace(
            file_orthogroup_table=str(input_path),
            mode="hog2og",
            dir_out=str(out_dir),
        )
    )

    orthogroups = pandas.read_csv(out_dir / "Orthogroups.tsv", sep="\t")
    gene_counts = pandas.read_csv(out_dir / "Orthogroups.GeneCount.tsv", sep="\t")

    assert orthogroups["Orthogroup"].tolist() == ["HOG0000001", "HOG0000002"]
    assert gene_counts.loc[0, "spA"] == 6000
    assert gene_counts.loc[0, "spB"] == 2
    assert gene_counts.loc[0, "Total"] == 6002
    assert gene_counts.loc[1, "spA"] == 0
    assert gene_counts.loc[1, "spB"] == 1
    assert gene_counts.loc[1, "Total"] == 1
    assert unrelated_tmp.read_text(encoding="utf-8") == "keep\n"


def test_hog2og_preserves_numeric_gene_ids_and_removes_only_leading_prefix(tmp_path):
    mod = load_module()
    input_path = tmp_path / "Ntsv1.tsv"
    out_dir = tmp_path / "Orthogroups"
    pandas.DataFrame(
        [
            {
                "OG": "OG0000001",
                "HOG": "Ntsv1.HOGNtsv1.0000001",
                "Gene Tree Parent Clade": "Ntsv1",
                "numeric_species": 101,
            },
            {
                "OG": "OG0000002",
                "HOG": "HOGNtsv1.0000002",
                "Gene Tree Parent Clade": "Ntsv1",
                "numeric_species": "NA",
            },
        ]
    ).to_csv(input_path, sep="\t", index=False)

    mod.run(
        SimpleNamespace(
            file_orthogroup_table=str(input_path),
            mode="hog2og",
            dir_out=str(out_dir),
        )
    )

    orthogroups = pandas.read_csv(
        out_dir / "Orthogroups.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    gene_counts = pandas.read_csv(out_dir / "Orthogroups.GeneCount.tsv", sep="\t")
    assert orthogroups["Orthogroup"].tolist() == [
        "HOGNtsv1.0000001",
        "HOGNtsv1.0000002",
    ]
    assert orthogroups["numeric_species"].tolist() == ["101", "NA"]
    assert gene_counts["numeric_species"].tolist() == [1, 1]


@pytest.mark.parametrize("second_id", ["N0.HOG0000004", "HOG0000004"])
def test_duplicate_hog_identity_refuses_publication_and_preserves_outputs(tmp_path, second_id):
    mod = load_module()
    source = tmp_path / "N0.tsv"
    source.write_text(
        "HOG\tOG\tGene Tree Parent Clade\tspA\n"
        "N0.HOG0000004\tOG0000002\tn116\tgene_a\n"
        f"{second_id}\tOG0000003\tn333\tgene_b\n",
        encoding="utf-8",
    )
    original = source.read_bytes()
    out = tmp_path / "published"
    out.mkdir()
    previous = {name: b"previous verified output\n" for name in
                ["Orthogroups.tsv", "Orthogroups.GeneCount.tsv", "README.txt"]}
    for name, data in previous.items():
        (out / name).write_bytes(data)
    with pytest.raises(ValueError, match="duplicate family identities: HOG0000004"):
        mod.run(SimpleNamespace(file_orthogroup_table=str(source), mode="hog2og", dir_out=str(out)))
    assert source.read_bytes() == original
    assert {name: (out / name).read_bytes() for name in previous} == previous


def test_empty_hog_identity_does_not_create_output_directory(tmp_path):
    mod = load_module()
    source = tmp_path / "N0.tsv"
    source.write_text("HOG\tOG\tGene Tree Parent Clade\tspA\n\tOG1\tn1\tgene_a\n", encoding="utf-8")
    out = tmp_path / "new_output"
    with pytest.raises(ValueError, match="empty family identity"):
        mod.run(SimpleNamespace(file_orthogroup_table=str(source), mode="hog2og", dir_out=str(out)))
    assert not out.exists()


@pytest.mark.parametrize("duplicate_index", [False, True])
def test_gene_counts_preserve_cell_semantics_index_and_column_order(duplicate_index):
    mod = load_module()
    cells = [
        None,
        float("nan"),
        pandas.NA,
        "",
        "NA",
        " ",
        ",",
        ",gene,,",
        101,
        0,
        True,
        -2.5,
        "gene1, gene2",
        "gene1\ngene2",
    ]
    index = pandas.Index([i % 3 for i in range(len(cells))], name="row") if duplicate_index else None
    families = [f"OG{i}" for i in range(len(cells))]
    source = pandas.DataFrame({"speciesA": cells, "Orthogroup": families, "speciesB": "gene"}, index=index)
    original = source.copy(deep=True)
    counts = [0, 0, 0, 0, 1, 1, 2, 4, 1, 1, 1, 1, 2, 1]
    expected = pandas.DataFrame(
        {"Orthogroup": families, "speciesB": 1, "speciesA": counts, "Total": [value + 1 for value in counts]},
        index=index,
    )

    actual = mod.build_gene_count_table(source, ["speciesB", "speciesA"])

    pandas.testing.assert_frame_equal(actual, expected)
    pandas.testing.assert_frame_equal(source, original)


@pytest.mark.parametrize("dtype", ["object", "string"])
def test_gene_counts_keep_integer_columns_for_empty_tables(dtype):
    mod = load_module()
    source = pandas.DataFrame({"Orthogroup": pandas.Series(dtype=dtype), "speciesA": pandas.Series(dtype=dtype)})
    expected = pandas.DataFrame(
        {"Orthogroup": source["Orthogroup"], "speciesA": pandas.Series(dtype=int), "Total": pandas.Series(dtype=int)}
    )
    pandas.testing.assert_frame_equal(mod.build_gene_count_table(source, ["speciesA"]), expected)


def test_gene_counts_without_species_keep_family_index():
    mod = load_module()
    source = pandas.DataFrame({"Orthogroup": ["OG1", "OG2"]}, index=[9, 9])
    expected = source.assign(Total=0)
    pandas.testing.assert_frame_equal(mod.build_gene_count_table(source, []), expected)


def test_sp2og_preserves_formatted_genes_and_integer_counts(tmp_path):
    mod = load_module()
    source = tmp_path / "sonicparanoid.tsv"
    source.write_text("group_id\tspA.fa\tspB.faa\n1\tg*1,g2*\t\n12\tNA\t007\n", encoding="utf-8")
    out = tmp_path / "output"
    mod.run(SimpleNamespace(file_orthogroup_table=str(source), mode="sp2og", dir_out=str(out)))

    assert (out / "Orthogroups.tsv").read_text() == ("Orthogroup\tspA\tspB\nSP0000001\tg1, g2\t\nSP0000012\tNA\t007\n")
    assert (out / "Orthogroups.GeneCount.tsv").read_text() == (
        "Orthogroup\tspA\tspB\tTotal\nSP0000001\t2\t0\t2\nSP0000012\t1\t1\t2\n"
    )
