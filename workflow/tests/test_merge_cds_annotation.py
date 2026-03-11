from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "merge_cds_annotation.py"


def load_module():
    spec = spec_from_file_location("merge_cds_annotation", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_load_expression_sets_gene_id_index(tmp_path):
    mod = load_module()
    infile = tmp_path / "expression.tsv"
    infile.write_text(
        "Unnamed: 0\tleaf\troot\n"
        "geneA\t1\t2\n"
        "geneB\t3\t4\n",
        encoding="utf-8",
    )

    out = mod.load_expression(str(infile))

    assert out.index.name == "gene_id"
    assert "gene_id" not in out.columns
    assert out.loc["geneA", "leaf"] == 1
    assert out.loc["geneB", "root"] == 4


def test_load_busco_aggregates_by_gene_id_and_indexes(tmp_path):
    mod = load_module()
    infile = tmp_path / "busco.tsv"
    infile.write_text(
        "# BUSCO v5\n"
        "BUSCO1\tComplete\tgeneA:cds1\t100\t200\turl1\tdesc1\n"
        "BUSCO2\tDuplicated\tgeneA:cds2\t150\t210\turl2\tdesc2\n"
        "BUSCO3\tFragmented\tgeneB\t90\t190\turl3\tdesc3\n",
        encoding="utf-8",
    )

    out = mod.load_busco(str(infile))

    assert out.index.name == "gene_id"
    assert "gene_id" not in out.columns
    assert out.loc["geneA", "busco_id"] == "BUSCO1; BUSCO2"
    assert out.loc["geneA", "busco_status"] == "Complete; Duplicated"
    assert out.loc["geneB", "busco_sequence"] == "geneB"


def test_join_if_available_accepts_preindexed_tables():
    mod = load_module()
    df = pandas.DataFrame(index=pandas.Index(["geneA", "geneB"], name="gene_id"))
    tmp = pandas.DataFrame(
        {"annotation": ["hitA", "hitB"]},
        index=pandas.Index(["geneA", "geneC"], name="gene_id"),
    )

    out = mod.join_if_available(df, tmp)

    assert out.index.tolist() == ["geneA", "geneB"]
    assert out.loc["geneA", "annotation"] == "hitA"
    assert pandas.isna(out.loc["geneB", "annotation"])
