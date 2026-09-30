import gzip
import subprocess
import sys
from pathlib import Path

import pandas
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "reformat_uniprot_diamond.py"


@pytest.mark.parametrize('with_query_fasta', [False, True])
@pytest.mark.parametrize('compressed', [False, True])
def test_diamond_annotations_preserve_literal_query_and_subject_ids(tmp_path, with_query_fasta, compressed):
    genes = ['NA', 'NULL', 'nan', '001']
    subjects = ['001', 'nan', 'NULL', 'NA']

    def write(name, text):
        path = tmp_path / (name + ('.gz' if compressed else ''))
        if compressed:
            with gzip.open(path, 'wt') as handle:
                handle.write(text)
        else:
            path.write_text(text)
        return path

    diamond = write('diamond.tsv', ''.join(f'{gene}\t{subject}\t90\t30\t1e-10\t200\t30\n'
                                           for gene, subject in zip(genes, subjects, strict=True))
                    + 'NA\tother\t95\t25\tNA\t190\t30\n')
    fasta = write('queries.fa', ''.join(f'>{gene}\nMKT\n' for gene in genes))
    uniprot = write('subjects.fa', ''.join(f'>{subject} desc_{subject}\nMKT\n' for subject in subjects))
    metadata = write('metadata.tsv', 'accession\tgene_name_primary\n' + ''.join(
        f'{subject}\tname_{subject}\n' for subject in subjects))
    output = tmp_path / 'annotations.tsv'
    command = [sys.executable, str(SCRIPT_PATH), '--diamond_tsv', str(diamond),
               '--uniprot_fasta', str(uniprot), '--uniprot_meta_tsv', str(metadata), '--outfile', str(output)]
    if with_query_fasta:
        command += ['--query_fasta', str(fasta)]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stderr == ''
    actual = pandas.read_csv(output, sep='\t', dtype=str, keep_default_na=False)
    expected_genes = genes if with_query_fasta else sorted(genes)
    assert actual.gene_id.tolist() == expected_genes
    mapping = dict(zip(genes, subjects, strict=True))
    for row in actual.itertuples():
        subject = mapping[row.gene_id]
        assert row.sprot_best == row.sprot_alias == subject
        assert row.sprot_recname == 'desc_' + subject
        assert row.gene_name_primary == 'name_' + subject
        assert row.sprot_coverage == '100' and row.sprot_identity == '90'
        assert row.sprot_evalue == '1.000e-10'


def test_diamond_reader_keeps_empty_ids_and_numeric_na_handling(tmp_path):
    from workflow.support.reformat_uniprot_diamond import load_best_diamond_hits

    path = tmp_path / 'missing.tsv'
    path.write_text('q1\t\tNA\tNULL\tnan\t\tNA\n'
                    '\tP1\t90\t30\t1e-10\t200\t30\n')
    result = load_best_diamond_hits(path)
    first = result.loc[result.qseqid == 'q1'].iloc[0]
    assert pandas.isna(first.sseqid)
    assert result.qseqid.isna().sum() == 1
    assert first[['pident', 'length', 'evalue', 'bitscore', 'qlen', 'coverage']].isna().all()


def test_reformat_uniprot_diamond_adds_metadata_columns_except_rnammer(tmp_path):
    query_fasta = tmp_path / "query.fa"
    query_fasta.write_text(">gene1\nMKT\n>gene2\nAAA\n", encoding="utf-8")

    diamond_tsv = tmp_path / "diamond.tsv"
    diamond_tsv.write_text("gene1\tP12345\t87.5\t80\t1e-50\t120\t100\n", encoding="utf-8")

    uniprot_fasta = tmp_path / "uniprot.pep"
    uniprot_fasta.write_text(">P12345 Protein kinase\nMKT\n", encoding="utf-8")

    uniprot_meta = tmp_path / "uniprot.meta.tsv"
    uniprot_meta.write_text(
        (
            "accession\tsignal_start\tsignal_end\ttransmem_aa\ttransmem_count\ttransmem_regions\tkegg_gene\t"
            "kegg_orthology\tgo_ids\tgo_aspects\tgo_terms\tgo_evidence\tgene_name_primary\tgene_name_synonyms\t"
            "ec_numbers\tsubcellular_location\tkeywords\tinterpro_ids\tpfam_ids\treactome_ids\torganism\ttaxid\n"
            "P12345\t1\t20\t44\t2\t45-67;70-90\thsa:1234\tK00844\tGO:0005524;GO:0004672\tF;P\t"
            "ATP binding; protein kinase activity\tIEA\tAKT1\tPKB\t2.7.11.1\tCytoplasm\tKinase;Transferase\t"
            "IPR000719\tPF00069\tR-HSA-123456\tHomo sapiens\t9606\n"
        ),
        encoding="utf-8",
    )

    outfile = tmp_path / "uniprot.annotation.tsv"
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--diamond_tsv",
            str(diamond_tsv),
            "--query_fasta",
            str(query_fasta),
            "--uniprot_fasta",
            str(uniprot_fasta),
            "--uniprot_meta_tsv",
            str(uniprot_meta),
            "--outfile",
            str(outfile),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr

    out = pandas.read_csv(outfile, sep="\t")
    assert out.columns.tolist() == [
        "gene_id",
        "sprot_best",
        "sprot_alias",
        "sprot_coverage",
        "sprot_identity",
        "sprot_evalue",
        "sprot_recname",
        "signal_start",
        "signal_end",
        "transmem_aa",
        "transmem_count",
        "transmem_regions",
        "kegg_gene",
        "kegg_orthology",
        "go_ids",
        "go_aspects",
        "go_terms",
        "go_evidence",
        "gene_name_primary",
        "gene_name_synonyms",
        "ec_numbers",
        "subcellular_location",
        "keywords",
        "interpro_ids",
        "pfam_ids",
        "reactome_ids",
        "organism",
        "taxid",
    ]
    assert "rnammer" not in out.columns
    assert "signalp_score" not in out.columns
    assert "tmhmm_predhel" not in out.columns

    row_hit = out.loc[out["gene_id"] == "gene1", :].iloc[0]
    assert row_hit["sprot_best"] == "P12345"
    assert row_hit["sprot_alias"] == "P12345"
    assert row_hit["sprot_coverage"] == 80
    assert row_hit["sprot_identity"] == 87.5
    assert float(row_hit["sprot_evalue"]) == 1e-50
    assert row_hit["sprot_recname"] == "Protein kinase"
    assert row_hit["signal_start"] == 1
    assert row_hit["signal_end"] == 20
    assert row_hit["transmem_aa"] == 44
    assert row_hit["transmem_count"] == 2
    assert row_hit["transmem_regions"] == "45-67;70-90"
    assert row_hit["kegg_gene"] == "hsa:1234"
    assert row_hit["kegg_orthology"] == "K00844"
    assert row_hit["go_ids"] == "GO:0005524;GO:0004672"
    assert row_hit["go_aspects"] == "F;P"
    assert row_hit["go_terms"] == "ATP binding; protein kinase activity"
    assert row_hit["go_evidence"] == "IEA"
    assert row_hit["gene_name_primary"] == "AKT1"
    assert row_hit["gene_name_synonyms"] == "PKB"
    assert row_hit["ec_numbers"] == "2.7.11.1"
    assert row_hit["subcellular_location"] == "Cytoplasm"
    assert row_hit["keywords"] == "Kinase;Transferase"
    assert row_hit["interpro_ids"] == "IPR000719"
    assert row_hit["pfam_ids"] == "PF00069"
    assert row_hit["reactome_ids"] == "R-HSA-123456"
    assert row_hit["organism"] == "Homo sapiens"
    assert row_hit["taxid"] == 9606

    row_no_hit = out.loc[out["gene_id"] == "gene2", :].iloc[0]
    assert pandas.isna(row_no_hit["sprot_best"])
    assert pandas.isna(row_no_hit["kegg_gene"])


def test_reformat_uniprot_diamond_supports_gzipped_inputs(tmp_path):
    query_fasta_gz = tmp_path / "query.fa.gz"
    with gzip.open(query_fasta_gz, "wt", encoding="utf-8") as handle:
        handle.write(">geneA\nMKT\n")

    uniprot_fasta_gz = tmp_path / "uniprot.pep.gz"
    with gzip.open(uniprot_fasta_gz, "wt", encoding="utf-8") as handle:
        handle.write(">P99999 Example protein\nMKT\n")

    uniprot_meta_gz = tmp_path / "uniprot.meta.tsv.gz"
    with gzip.open(uniprot_meta_gz, "wt", encoding="utf-8") as handle:
        handle.write(
            (
                "accession\tsignal_start\tsignal_end\ttransmem_aa\ttransmem_count\ttransmem_regions\tkegg_gene\t"
                "kegg_orthology\tgo_ids\tgo_aspects\tgo_terms\tgo_evidence\tgene_name_primary\tgene_name_synonyms\t"
                "ec_numbers\tsubcellular_location\tkeywords\tinterpro_ids\tpfam_ids\treactome_ids\torganism\ttaxid\n"
                "P99999\t2\t18\t0\t\t\tath:AT1G01010\tK00001\tGO:0003674\tF\tmolecular_function\tIEA\tGENEA\t\t"
                "1.1.1.1\tPlastid\tPlastid\tIPR000001\tPF00001\tR-ATH-1\tArabidopsis thaliana\t3702\n"
            )
        )

    diamond_tsv = tmp_path / "diamond.tsv"
    diamond_tsv.write_text("geneA\tP99999\t90\t30\t2e-10\t60\t30\n", encoding="utf-8")

    outfile = tmp_path / "out.tsv"
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT_PATH),
            "--diamond_tsv",
            str(diamond_tsv),
            "--query_fasta",
            str(query_fasta_gz),
            "--uniprot_fasta",
            str(uniprot_fasta_gz),
            "--uniprot_meta_tsv",
            str(uniprot_meta_gz),
            "--outfile",
            str(outfile),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr

    out = pandas.read_csv(outfile, sep="\t")
    assert out.shape[0] == 1
    assert out.loc[0, "gene_id"] == "geneA"
    assert out.loc[0, "sprot_best"] == "P99999"
    assert out.loc[0, "sprot_recname"] == "Example protein"
    assert out.loc[0, "kegg_gene"] == "ath:AT1G01010"
    assert out.loc[0, "go_ids"] == "GO:0003674"
    assert out.loc[0, "organism"] == "Arabidopsis thaliana"
