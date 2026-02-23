import pandas

from workflow.support.gff2genestat import add_intron_info


def test_add_intron_info_uses_semicolon_delimiter():
    gff = pandas.DataFrame(
        {
            "gene_id": ["gene1", "gene1", "gene1"],
            "start": [1, 21, 41],
            "end": [10, 30, 50],
        }
    )
    df_all = pandas.DataFrame(columns=["gene_id", "feature_size", "num_intron", "intron_positions"])

    out = add_intron_info(gff=gff, df_all=df_all, id_col="gene_id")
    row = out.loc[out["gene_id"] == "gene1"].iloc[0]

    assert row["num_intron"] == 2
    assert row["intron_positions"] == "10;20"
    assert "," not in row["intron_positions"]


def test_add_intron_info_keeps_empty_when_no_intron():
    gff = pandas.DataFrame(
        {
            "gene_id": ["gene2"],
            "start": [100],
            "end": [150],
        }
    )
    df_all = pandas.DataFrame(columns=["gene_id", "feature_size", "num_intron", "intron_positions"])

    out = add_intron_info(gff=gff, df_all=df_all, id_col="gene_id")
    row = out.loc[out["gene_id"] == "gene2"].iloc[0]

    assert row["num_intron"] == 0
    assert row["intron_positions"] == ""
