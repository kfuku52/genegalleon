import subprocess
import sys
from pathlib import Path

import pandas
import pytest

from workflow.support.gff2genestat import (
    add_id_column,
    add_intron_info,
    extract_by_ids,
    process_single_gff,
    summarize_gene_features,
    trim_species_prefix,
)

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "gff2genestat.py"
OUT_COLS = ["gene_id", "feature_size", "num_intron", "intron_positions", "chromosome", "start", "end", "strand", "feature_blocks", "feature_type"]


def test_longest_selects_one_transcript_before_summarizing():
    gff = pandas.DataFrame(
        [
            ["chr1", "gene", 1, 500, "+", "ID=gene1;"],
            ["chr1", "mRNA", 1, 500, "+", "ID=long;Parent=gene1;"],
            ["chr1", "mRNA", 1, 60, "+", "ID=short;Parent=gene1;"],
            ["chr1", "CDS", 401, 500, "+", "ID=long.c2;Parent=long;"],
            ["chr1", "CDS", 1, 200, "+", "ID=long.c1;Parent=long;"],
            ["chr1", "CDS", 1, 60, "+", "ID=short.c1;Parent=short;"],
            ["chr1", "CDS", 1, 200, "+", "ID=long.duplicate;Parent=long;"],
        ], columns=["sequence", "feature", "start", "end", "strand", "attributes"]
    )
    for ordered in (gff, gff.iloc[::-1]):
        selected = extract_by_ids(ordered, pandas.Series(["Species_a_gene1"]), "CDS", "longest")
        out = summarize_gene_features(selected, OUT_COLS)
        assert len(out) == 1
        assert out.iloc[0]["feature_size"] == 300
        assert out.iloc[0]["num_intron"] == 1
        assert out.iloc[0]["intron_positions"] == "200"
        assert out.iloc[0]["feature_blocks"] == "1-200;401-500"
        assert out.iloc[0]["feature_type"] == "CDS"


def test_longest_gtf_transcripts_do_not_merge_isoforms():
    gff = pandas.DataFrame(
        [
            ["chr1", "CDS", 1, 90, "+", 'gene_id "gene1"; transcript_id "t1";'],
            ["chr1", "CDS", 1, 30, "+", 'gene_id "gene1"; transcript_id "t2";'],
        ], columns=["sequence", "feature", "start", "end", "strand", "attributes"]
    )
    out = extract_by_ids(gff, pandas.Series(["Species_a_gene1"]), "CDS", "longest")
    assert out["end"].tolist() == [90]


def test_longest_breaks_equal_length_ties_by_transcript_id(capsys):
    gff = pandas.DataFrame(
        [
            ["chr1", "CDS", 101, 190, "+", "Parent=gene1.t2;"],
            ["chr1", "CDS", 1, 90, "+", "Parent=gene1.t1;"],
        ], columns=["sequence", "feature", "start", "end", "strand", "attributes"]
    )
    selected = extract_by_ids(gff, pandas.Series(["Species_a_gene1"]), "CDS", "longest")

    assert selected[["start", "end"]].values.tolist() == [[1, 90]]
    assert (
        "Equal-length longest transcripts for Species_a_gene1; "
        "selected gene1.t1 deterministically from gene1.t1,gene1.t2."
    ) in capsys.readouterr().err


def test_reverse_strand_cds_uses_transcript_order_and_genomic_bounds():
    gff = pandas.DataFrame(
        {"gene_id": ["gene1"] * 3, "sequence": ["chr1"] * 3,
         "strand": ["-"] * 3, "start": [1, 201, 101], "end": [10, 230, 120]}
    )
    out = summarize_gene_features(gff, OUT_COLS).iloc[0]
    assert out["feature_size"] == 60
    assert out["intron_positions"] == "30;50"
    assert out["feature_blocks"] == "201-230;101-120;1-10"
    assert (out["start"], out["end"]) == (1, 230)


@pytest.mark.parametrize("sequence,strand", [("chr2", "+"), ("chr1", "-")])
def test_summary_rejects_mixed_coordinate_systems(sequence, strand):
    gff = pandas.DataFrame(
        {"gene_id": ["gene1", "gene1"], "sequence": ["chr1", sequence],
         "strand": ["+", strand], "start": [1, 101], "end": [10, 120]}
    )
    with pytest.raises(ValueError, match="coordinate"):
        summarize_gene_features(gff, OUT_COLS)


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


def test_trim_species_prefix_preserves_dotted_taxonomic_qualifier_boundaries():
    assert trim_species_prefix("Asimitellaria_furusei_var._furusei_gene1") == "gene1"
    assert trim_species_prefix("Asimitellaria_furusei_var._subramosa_gene1") == "gene1"
    assert trim_species_prefix("Arisaema_sp._aooni_gene1") == "gene1"
    assert trim_species_prefix("Bacillus_subtilis_subsp._subtilis_gene1") == "gene1"


def test_main_does_not_assign_short_species_gff_to_longer_species_seqfile(tmp_path: Path):
    gff_dir = tmp_path / "gff"
    gff_dir.mkdir()
    seqfile = tmp_path / "Species_a_subsp_x.fa"
    outfile = tmp_path / "out.tsv"
    seqfile.write_text(">Species_a_subsp_x_gene1\nATGAAA\n", encoding="utf-8")
    (gff_dir / "Species_a.gff").write_text(
        "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=gene1.CDS1;Parent=gene1;\n",
        encoding="utf-8",
    )

    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--dir_gff",
            str(gff_dir),
            "--seqfile",
            str(seqfile),
            "--outfile",
            str(outfile),
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stdout + completed.stderr
    out = pandas.read_csv(outfile, sep="\t")
    assert out.empty


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


def test_process_single_gff_skips_file_with_fewer_than_nine_columns(tmp_path, capsys):
    gff_dir = tmp_path / "gff"
    gff_dir.mkdir()
    (gff_dir / "Arabidopsis_thaliana.bad.gff3").write_text("chr1\tsrc\tCDS\n", encoding="utf-8")
    out_cols = ["gene_id", "feature_size", "num_intron", "intron_positions", "chromosome", "start", "end", "strand"]
    gff_cols = ["sequence", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"]

    out = process_single_gff(
        gff_file="Arabidopsis_thaliana.bad.gff3",
        dir_gff=str(gff_dir),
        seq_sp_values=["Arabidopsis_thaliana_gene1"],
        feature="CDS",
        multiple_hits="longest",
        gff_cols=gff_cols,
        out_cols=out_cols,
    )

    captured = capsys.readouterr()
    assert out.empty
    assert "Skipping malformed GFF with fewer than 9 columns" in captured.err


def test_utr_annotation_follows_selected_transcript_and_preserves_phase():
    from workflow.support.gff2genestat import attach_transcript_structure
    gff = pandas.DataFrame([
        ["chr1", "CDS", 101, 160, "-", "1", "Parent=gene1.long"],
        ["chr1", "CDS", 21, 50, "-", "1", "Parent=gene1.long"],
        ["chr1", "CDS", 101, 130, "-", "0", "Parent=gene1.short"],
        ["chr1", "five_prime_UTR", 161, 180, "-", ".", "Parent=gene1.long"],
        ["chr1", "three_prime_UTR", 1, 20, "-", ".", "Parent=gene1.long"],
        ["chr1", "five_prime_UTR", 131, 250, "-", ".", "Parent=gene1.short"],
    ], columns=["sequence", "feature", "start", "end", "strand", "phase", "attributes"])
    selected = extract_by_ids(gff, pandas.Series(["Species_a_gene1"]), "CDS", "longest")
    annotated = attach_transcript_structure(selected, gff)
    cols = OUT_COLS + ["gff_transcript_id", "utr_blocks", "cds_first_phase"]
    row = summarize_gene_features(annotated, cols).iloc[0]
    assert row.gff_transcript_id == "gene1.long"
    assert row.utr_blocks == "161-180;1-20"
    assert row.cds_first_phase == 1
    # A transcript interval alone does not establish UTR sequence.
    no_utr = attach_transcript_structure(selected, gff[gff.feature == "CDS"])
    assert no_utr.utr_blocks.eq("").all()


@pytest.mark.parametrize('strand,starts,ends', [('+',[1,4,20],[3,6,22]), ('-',[20,17,1],[22,19,3])])
def test_adjacent_cds_fragments_are_not_introns(strand, starts, ends):
    gff = pandas.DataFrame(dict(gene_id=['g']*3,sequence=['chr1']*3,strand=[strand]*3,start=starts,end=ends))
    row = summarize_gene_features(gff,OUT_COLS).iloc[0]
    assert row.feature_size == 9
    assert row.num_intron == 1
    assert row.intron_positions == '6'
    legacy = add_intron_info(gff, pandas.DataFrame(columns=OUT_COLS)).iloc[0]
    assert legacy.num_intron == 1 and legacy.intron_positions == '6'


def test_all_cds_phases_are_validated_and_missing_initial_phase_is_inferred():
    from workflow.support.gff2genestat import attach_transcript_structure
    cds = pandas.DataFrame(dict(gene_id=['g','g'],selected_transcript=['t','t'],
        sequence=['chr1','chr1'],feature=['CDS','CDS'],start=[1,10],end=[4,14],
        strand=['+','+'],phase=[float("nan"),2.0],attributes=['Parent=t','Parent=t']))
    annotated = attach_transcript_structure(cds,cds)
    assert annotated.cds_first_phase.eq(0).all()
    cds['phase'] = ['0','1']
    with pytest.raises(ValueError,match='Conflicting CDS phases'):
        attach_transcript_structure(cds,cds)
    cds['phase'] = ['.','.']
    assert attach_transcript_structure(cds,cds).cds_first_phase.isna().all()


def test_cds_length_validation_rejects_wrong_transcript_and_preserves_missing():
    from workflow.support.gff2genestat import validate_cds_lengths
    traits=pandas.DataFrame(dict(gene_id=['a'],feature_size=[6]))
    records=[('a','a','ATGNNN'),('missing','missing','ATG')]
    validate_cds_lengths(traits,records)
    traits.loc[0,'feature_size']=3
    with pytest.raises(ValueError,match='GFF=3, CDS=6'):
        validate_cds_lengths(traits,records)
    with pytest.raises(ValueError,match='ungapped nucleotide'):
        validate_cds_lengths(traits,[('a','a','MKE')])
