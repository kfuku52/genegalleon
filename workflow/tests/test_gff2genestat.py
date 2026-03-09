import pandas

from workflow.support.gff2genestat import add_id_column
from workflow.support.gff2genestat import add_intron_info
from workflow.support.gff2genestat import extract_by_ids
from workflow.support.gff2genestat import summarize_gene_features


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


def test_add_id_column_assigns_full_seq_name_from_parent_gene_id():
    gff = pandas.DataFrame(
        {
            "attributes": [
                "ID=gene1.CDS1;Parent=gene1;",
                "ID=gene2.CDS1;Parent=gene2;",
            ]
        }
    )
    seq_names = pandas.Series(
        [
            "Arabidopsis_thaliana_gene1",
            "Arabidopsis_thaliana_gene2",
        ]
    )

    out = add_id_column(gff=gff, seq_names=seq_names)

    assert out.loc[0, "gene_id"] == "Arabidopsis_thaliana_gene1"
    assert out.loc[1, "gene_id"] == "Arabidopsis_thaliana_gene2"


def test_add_id_column_accepts_hyphen_to_underscore_gene_name_variant():
    gff = pandas.DataFrame(
        {
            "attributes": [
                "ID=gene_1.CDS1;Parent=gene_1;",
            ]
        }
    )
    seq_names = pandas.Series(
        [
            "Arabidopsis_thaliana_gene-1",
        ]
    )

    out = add_id_column(gff=gff, seq_names=seq_names)

    assert out.loc[0, "gene_id"] == "Arabidopsis_thaliana_gene-1"


def test_add_id_column_accepts_ncbi_dbxref_geneid_namespace():
    gff = pandas.DataFrame(
        {
            "attributes": [
                "ID=cds-XP_024525354.1;Parent=rna-XM_024669586.1;Dbxref=GeneID:112340394,Genbank:XP_024525354.1;Name=XP_024525354.1;gene=LOC112340394;protein_id=XP_024525354.1;",
            ]
        }
    )
    seq_names = pandas.Series(
        [
            "Selaginella_moellendorffii_GeneID112340394",
        ]
    )

    out = add_id_column(gff=gff, seq_names=seq_names)

    assert out.loc[0, "gene_id"] == "Selaginella_moellendorffii_GeneID112340394"


def test_extract_by_ids_matches_cds_parent_prefix_without_regex_scan():
    gff = pandas.DataFrame(
        {
            "feature": ["CDS"],
            "attributes": ["ID=gene1.1.model.CDS1;Parent=gene1.1.model;"],
        }
    )
    seq_names = pandas.Series(["Arabidopsis_thaliana_gene1.1"])

    out = extract_by_ids(gff=gff, seq_names=seq_names, feature="CDS", multiple_hits="longest")

    assert out.shape[0] == 1
    assert out.iloc[0]["gene_id"] == "Arabidopsis_thaliana_gene1.1"


def test_extract_by_ids_resolves_ncbi_cds_via_dbxref_geneid():
    gff = pandas.DataFrame(
        {
            "feature": ["CDS"],
            "attributes": [
                "ID=cds-XP_024525354.1;Parent=rna-XM_024669586.1;Dbxref=GeneID:112340394,Genbank:XP_024525354.1;Name=XP_024525354.1;gene=LOC112340394;protein_id=XP_024525354.1;",
            ],
        }
    )
    seq_names = pandas.Series(["Selaginella_moellendorffii_GeneID112340394"])

    out = extract_by_ids(gff=gff, seq_names=seq_names, feature="CDS", multiple_hits="longest")

    assert out.shape[0] == 1
    assert out.iloc[0]["gene_id"] == "Selaginella_moellendorffii_GeneID112340394"


def test_extract_by_ids_resolves_cds_via_parent_feature_match():
    gff = pandas.DataFrame(
        {
            "feature": ["mRNA", "CDS"],
            "attributes": [
                "ID=mrna1;Name=gene1.1;Parent=gene1;",
                "ID=cds1;Parent=mrna1;",
            ],
        }
    )
    seq_names = pandas.Series(["Arabidopsis_thaliana_gene1.1"])

    out = extract_by_ids(gff=gff, seq_names=seq_names, feature="CDS", multiple_hits="longest")

    assert out.shape[0] == 1
    assert out.iloc[0]["gene_id"] == "Arabidopsis_thaliana_gene1.1"


def test_extract_by_ids_resolves_exact_parent_chain_via_gene_feature():
    gff = pandas.DataFrame(
        {
            "feature": ["gene", "mRNA", "CDS"],
            "attributes": [
                "ID=gene1;",
                "ID=mrna1;Parent=gene1;",
                "ID=cds1;Parent=mrna1;",
            ],
        }
    )
    seq_names = pandas.Series(["Arabidopsis_thaliana_gene1"])

    out = extract_by_ids(gff=gff, seq_names=seq_names, feature="CDS", multiple_hits="longest")

    assert out.shape[0] == 1
    assert out.iloc[0]["gene_id"] == "Arabidopsis_thaliana_gene1"


def test_extract_by_ids_resolves_namespaced_parent_identifier():
    gff = pandas.DataFrame(
        {
            "feature": ["CDS"],
            "attributes": [
                "ID=transcript:AT1G08450.1.Araport11.447.CDS.1;Parent=transcript:AT1G08450.1.Araport11.447;",
            ],
        }
    )
    seq_names = pandas.Series(["Arabidopsis_thaliana_AT1G08450.1"])

    out = extract_by_ids(gff=gff, seq_names=seq_names, feature="CDS", multiple_hits="longest")

    assert out.shape[0] == 1
    assert out.iloc[0]["gene_id"] == "Arabidopsis_thaliana_AT1G08450.1"


def test_extract_by_ids_accepts_gtf_style_transcript_id():
    gff = pandas.DataFrame(
        {
            "feature": ["CDS"],
            "attributes": [
                'gene_id "AT1G08450"; transcript_id "AT1G08450.1";',
            ],
        }
    )
    seq_names = pandas.Series(["Arabidopsis_thaliana_AT1G08450.1"])

    out = extract_by_ids(gff=gff, seq_names=seq_names, feature="CDS", multiple_hits="longest")

    assert out.shape[0] == 1
    assert out.iloc[0]["gene_id"] == "Arabidopsis_thaliana_AT1G08450.1"


def test_summarize_gene_features_handles_interleaved_gene_rows():
    gff = pandas.DataFrame(
        {
            "gene_id": ["gene1", "gene2", "gene1"],
            "sequence": ["chr1", "chr2", "chr1"],
            "strand": ["+", "-", "+"],
            "start": [1, 100, 21],
            "end": [10, 120, 30],
        }
    )
    out_cols = ["gene_id", "feature_size", "num_intron", "intron_positions", "chromosome", "start", "end", "strand"]

    out = summarize_gene_features(gff=gff, out_cols=out_cols)

    row1 = out.loc[out["gene_id"] == "gene1"].iloc[0]
    row2 = out.loc[out["gene_id"] == "gene2"].iloc[0]
    assert row1["feature_size"] == 20
    assert row1["num_intron"] == 1
    assert row1["intron_positions"] == "10"
    assert row1["chromosome"] == "chr1"
    assert row1["start"] == 1
    assert row1["end"] == 30
    assert row2["feature_size"] == 21
    assert row2["num_intron"] == 0
    assert row2["intron_positions"] == ""
