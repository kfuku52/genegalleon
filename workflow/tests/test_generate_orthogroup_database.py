from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "script" / "generate_orthogroup_database.py"


def load_module():
    spec = spec_from_file_location("generate_orthogroup_database", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_read_header_columns_returns_tsv_header_columns(tmp_path):
    mod = load_module()
    infile = tmp_path / "sample.tsv"
    infile.write_text("col1\tcol2\tcol3\n1\t2\t3\n", encoding="utf-8")

    assert mod.read_header_columns(str(infile)) == ["col1", "col2", "col3"]


def test_process_files_uses_single_read_csv_call(tmp_path, monkeypatch):
    mod = load_module()
    infile = tmp_path / "OG0001.test.tsv"
    pandas.DataFrame({"a": [1, 2], "b": [3, 4]}).to_csv(infile, sep="\t", index=False)

    call_count = {"n": 0}
    original_read_csv = mod.pd.read_csv

    def spy_read_csv(*args, **kwargs):
        call_count["n"] += 1
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(mod.pd, "read_csv", spy_read_csv)
    out = mod.process_files(str(infile), ["orthogroup", "a", "b"])

    assert call_count["n"] == 1
    assert out.columns.tolist() == ["orthogroup", "a", "b"]
    assert out["orthogroup"].iloc[0] == "OG0001"


def test_process_files_returns_empty_dataframe_when_columns_missing(tmp_path):
    mod = load_module()
    infile = tmp_path / "OG0002.test.tsv"
    pandas.DataFrame({"a": [1], "b": [2]}).to_csv(infile, sep="\t", index=False)

    out = mod.process_files(str(infile), ["orthogroup", "a", "b", "c"])
    assert out.empty


def test_parse_cutoff_stat_parses_valid_tokens_and_ignores_invalid():
    mod = load_module()
    parsed = mod.parse_cutoff_stat("OCNany2spe,0.8|badtoken|,1|ABC,notfloat|XYZ,1.5")
    assert parsed == [("OCNany2spe", 0.8), ("XYZ", 1.5)]


def test_apply_cutoff_accepts_preparsed_cutoff_list():
    mod = load_module()
    df = pandas.DataFrame(
        [
            {"A": "0.9", "B": "0.1"},
            {"A": "0.7", "B": "0.9"},
            {"A": "0.95", "B": "0.95"},
        ]
    )
    out = mod.apply_cutoff(df, [("A", 0.8), ("B", 0.5)])
    assert out.shape[0] == 1
    assert out.iloc[0]["A"] == "0.95"
