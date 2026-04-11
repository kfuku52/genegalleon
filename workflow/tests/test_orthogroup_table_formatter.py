from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import pandas


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
    huge_cell = ", ".join(f"gene{i}" for i in range(6000))

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
