import subprocess
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "merge_cds_annotation.py"


def load_module():
    spec = spec_from_file_location("merge_cds_annotation", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize('reader', ['uniprot', 'cdskit_localize', 'gff_info', 'fx2tab', 'expression', 'mmseqs', 'busco'])
@pytest.mark.parametrize('identifiers', [['001', '010'], ['NA', 'nan', 'None'],
                                      ['gene#1', 'gene:001', 'gene:1.2', 'gene:cds1']])
def test_annotation_loaders_preserve_lexical_identifiers(tmp_path, reader, identifiers):
    mod = load_module()
    path = tmp_path / 'annotation.tsv'
    if reader == 'busco':
        text = ''.join(f'B{i}\tComplete\t{gene}\t100\t200\tu\td\n' for i, gene in enumerate(identifiers))
    elif reader == 'mmseqs':
        text = ''.join(f'{gene}\t1\tspecies\tname\t1\t1\t1\t0.5\t1\n' for gene in identifiers)
    else:
        id_column = {'cdskit_localize': 'seq_id', 'fx2tab': '#id', 'expression': 'Identifier'}.get(reader, 'gene_id')
        text = f'{id_column}\tmetric\n' + ''.join(f'{gene}\t7\n' for gene in identifiers)
    path.write_text(text)
    loaded = getattr(mod, 'load_' + reader)(str(path))
    assert loaded.index.tolist() == identifiers
    assert loaded.index.name == 'gene_id'


def test_missing_busco_sequence_does_not_become_a_nan_gene(tmp_path):
    mod = load_module()
    path = tmp_path / 'busco.tsv'
    path.write_text('missing\tMissing\t\t\t\t\t\nreal\tComplete\tnan\t7\t9\tu\td\n')
    result = mod.load_busco(str(path))
    assert result.index.tolist() == ['nan']
    assert result.loc['nan', 'busco_id'] == 'real'
    path.write_text('missing\tMissing\t\t\t\t\t\n')
    assert mod.load_busco(str(path)) is None


def test_busco_annotation_reader_strips_only_terminal_coordinate_ranges(tmp_path):
    mod = load_module()
    path = tmp_path / "busco.tsv"
    path.write_text("# comment\nb1\tComplete\tgene:001:1-9\t100\t3\tu\tcontains#literal\n")
    result = mod.load_busco(str(path))
    assert result.index.tolist() == ["gene:001"]
    assert result.loc["gene:001", "busco_description"] == "contains#literal"


def test_busco_missing_metadata_keeps_legacy_nan_text(tmp_path):
    mod = load_module()
    path = tmp_path / 'busco.tsv'
    path.write_text('B1\tComplete\tgeneA\t100\t200\t\t\nB2\tDuplicated\tgeneA\t\t210\turl\tdesc\n')
    result = mod.load_busco(str(path))
    assert result.loc['geneA', 'busco_score'] == '100.0; nan'
    assert result.loc['geneA', 'busco_description'] == 'nan; desc'


def test_expression_renames_only_the_identifier_header(tmp_path):
    mod = load_module()
    path = tmp_path / 'expression.tsv'
    path.write_text('Identifier\tIdentifier_score\n001\t7\n')
    result = mod.load_expression(str(path))
    assert result.columns.tolist() == ['Identifier_score']
    assert result.loc['001', 'Identifier_score'] == 7


def test_fx2tab_renames_only_the_exact_standard_headers(tmp_path):
    mod = load_module()
    path = tmp_path / 'fx2tab.tsv'
    path.write_text('#id\tlength\tlength_ratio\t#identity\n001\t99\t0.5\t7\n')
    result = mod.load_fx2tab(str(path))
    assert result.columns.tolist() == ['cds_length', 'length_ratio', '#identity']
    assert result.loc['001', 'cds_length'] == 99


def test_text_identifier_parsing_retains_missing_metadata(tmp_path):
    mod = load_module()
    path = tmp_path / 'uniprot.tsv'
    path.write_text('gene_id\tvalue\nNA\tNA\n001\t7\n')
    result = mod.load_uniprot(str(path))
    assert pandas.isna(result.loc['NA', 'value'])
    assert result.loc['001', 'value'] == 7


@pytest.mark.parametrize('identifiers', [['001', '010'], ['NA', 'nan', 'None']])
def test_orthogroup_map_preserves_literal_member_and_group_ids(tmp_path, identifiers):
    mod = load_module()
    path = tmp_path / 'orthogroups.tsv'
    path.write_text('Orthogroup\tSpecies\n' + ''.join(f'{gene}\t{gene}\n' for gene in identifiers)
                    + 'empty\t\n')
    result = mod.load_orthogroup_map(str(path), 'Species')
    assert result.index.tolist() == identifiers
    assert result.tolist() == identifiers


@pytest.mark.parametrize('rows,identifiers,groups', [
    ('\t , 001, NA, , 001\nOG2\tNA, beta\nOG3\t beta, ,gamma\nempty\t , , \n',
     ['001', 'NA', 'beta', 'gamma'], [float('nan'), float('nan'), 'OG2', 'OG3']),
    ('OG1\t , , \nOG2\t\u2003, \n', [], []),
])
def test_orthogroup_map_keeps_first_missing_group_and_empty_mapping_schema(tmp_path, rows, identifiers, groups):
    mod = load_module()
    path = tmp_path / 'orthogroups.tsv'
    path.write_text('Orthogroup\tSpecies\n' + rows)
    # The parser's group dtype remains part of the Series schema, including
    # when every nonempty source cell contains only blank member tokens.
    parsed = pandas.read_csv(path, sep='\t', converters={
        'Orthogroup': mod.identifier_text, 'Species': mod.identifier_text,
    })
    expected = pandas.Series(groups, dtype=parsed['Orthogroup'].dtype, name='Orthogroup',
                             index=pandas.Index(identifiers, dtype=str, name='gene_id'))

    result = mod.load_orthogroup_map(str(path), 'Species')

    pandas.testing.assert_series_equal(result, expected, check_exact=True)


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
        "BUSCO1\tComplete\tgeneA:1-600\t100\t200\turl1\tdesc1\n"
        "BUSCO2\tDuplicated\tgeneA:4-633\t150\t210\turl2\tdesc2\n"
        "BUSCO3\tFragmented\tgeneB\t90\t190\turl3\tdesc3\n",
        encoding="utf-8",
    )

    out = mod.load_busco(str(infile))

    assert out.index.name == "gene_id"
    assert "gene_id" not in out.columns
    assert out.loc["geneA", "busco_id"] == "BUSCO1; BUSCO2"
    assert out.loc["geneA", "busco_status"] == "Complete; Duplicated"
    assert out.loc["geneB", "busco_sequence"] == "geneB"


@pytest.mark.parametrize('genes', [
    ['geneC', '001', 'NA', 'nan'],
    ['geneC', '001', 'geneC', 'NA', '001', 'nan'],
    ['geneC', '001', 'geneC', '001'],
])
def test_busco_annotation_preserves_singletons_duplicates_order_and_dtypes(tmp_path, genes):
    """Single records and interleaved duplicate groups share the same schema."""
    mod = load_module()
    path = tmp_path / 'busco.tsv'
    rows = [[f'B{i}', 'Complete', f'{gene}:{i * 3 + 1}-{i * 3 + 600}', '100', '200', f'url{i}', f'desc{i}']
            for i, gene in enumerate(genes)]
    path.write_text(''.join('\t'.join(row) + '\n' for row in rows), encoding='utf-8')
    grouped = {}
    for gene, row in zip(genes, rows, strict=True):
        grouped.setdefault(gene, []).append(row)
    expected_rows = [[gene, *('; '.join(row[column] for row in records) for column in range(7))]
                     for gene, records in grouped.items()]
    expected = pandas.DataFrame(expected_rows, columns=[
        'gene_id', 'busco_id', 'busco_status', 'busco_sequence', 'busco_score', 'busco_length',
        'busco_orthodb_url', 'busco_description',
    ]).astype(str).set_index('gene_id')

    pandas.testing.assert_frame_equal(mod.load_busco(str(path)), expected)


def test_busco_annotation_preserves_quoted_metadata_and_missing_numeric_spelling(tmp_path):
    mod = load_module()
    path = tmp_path / 'busco.tsv'
    path.write_text('B1\tComplete\tgeneC:1-600\t100\t200\tu\t"説明 α\tfield\ncontinued"\n'
                    'B2\tFragmented\tgeneA:4-300\t\t\t\t\n'
                    'B3\tDuplicated\tgeneC:7-609\t101\t201\tNA\tnan\n', encoding='utf-8')
    result = mod.load_busco(str(path))
    assert result.index.tolist() == ['geneC', 'geneA']
    assert result.loc['geneC', 'busco_description'] == '説明 α\tfield\ncontinued; nan'
    assert result.loc['geneC', 'busco_score'] == '100.0; 101.0'
    assert result.loc['geneA', 'busco_score'] == 'nan'
    assert result.loc['geneA', 'busco_length'] == 'nan'
    assert result.loc['geneA', 'busco_orthodb_url'] == 'nan'


@pytest.mark.parametrize('invalid_tail', [
    b'B2\tComplete\tgeneB\t100\t200\tu\t\xff\n',
    b'B2\tComplete\tgeneB\t100\t200\tu\td\textra\n',
    b'B2\tComplete\tgeneB\t100\t200\tu\t"unterminated\n',
])
def test_merge_cli_reads_complete_busco_before_replacing_output(tmp_path, invalid_tail):
    fasta = tmp_path / 'cds.fa'
    fasta.write_text('>geneA\nATG\n')
    busco = tmp_path / 'busco.tsv'
    busco.write_bytes(b'B1\tComplete\tgeneA\t100\t200\tu\td\n' + invalid_tail)
    output = tmp_path / 'annotation.tsv'
    output.write_bytes(b'previous complete annotation\n')
    result = subprocess.run([sys.executable, str(SCRIPT_PATH), '--cds_fasta', str(fasta),
                             '--busco_tsv', str(busco), '--out_tsv', str(output)],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert output.read_bytes() == b'previous complete annotation\n'


def test_load_cdskit_localize_prefixes_columns_and_indexes(tmp_path):
    mod = load_module()
    infile = tmp_path / "cdskit_localize.tsv"
    infile.write_text(
        "seq_id\tpredicted_class\tp_noTP\tp_SP\tp_mTP\tp_cTP\tp_lTP\tp_peroxisome\tperox_signal_type\n"
        "geneA\tSP\t0.1\t0.7\t0.1\t0.05\t0.05\t0.02\t-\n"
        "geneB\tperoxisome\t0.2\t0.1\t0.1\t0.1\t0.1\t0.4\tPTS1\n",
        encoding="utf-8",
    )

    out = mod.load_cdskit_localize(str(infile))

    assert out.index.name == "gene_id"
    assert "seq_id" not in out.columns
    assert "predicted_class" not in out.columns
    assert out.loc["geneA", "cdskit_localize_predicted_class"] == "SP"
    assert out.loc["geneB", "cdskit_localize_p_peroxisome"] == 0.4
    assert out.loc["geneB", "cdskit_localize_perox_signal_type"] == "PTS1"


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


@pytest.mark.parametrize('ncpu', [1, 2])
def test_merge_cli_preserves_identifiers_and_exact_metric_names(tmp_path, ncpu):
    ids = ['nan', '001', 'NA', 'None']
    fasta = tmp_path / 'cds.fa'
    fasta.write_text(''.join(f'>{gene}\nACGT\n' for gene in ids))
    orthogroups = tmp_path / 'orthogroups.tsv'
    orthogroups.write_text('Orthogroup\tSpecies\n' + ''.join(f'{gene}\t{gene}\n' for gene in ids))
    uniprot = tmp_path / 'uniprot.tsv'
    uniprot.write_text('gene_id\tuniprot_label\n' + ''.join(f'{gene}\tlabel{i}\n' for i, gene in enumerate(ids[::-1])))
    expression = tmp_path / 'expression.tsv'
    expression.write_text('Identifier\tIdentifier_score\n' + ''.join(f'{gene}\t7\n' for gene in ids))
    fx2tab = tmp_path / 'fx2tab.tsv'
    fx2tab.write_text('#id\tlength\tlength_ratio\n' + ''.join(f'{gene}\t99\t0.5\n' for gene in ids))
    busco = tmp_path / 'busco.tsv'
    busco.write_text(''.join(f'B{i}\tComplete\t{gene}\t100\t200\t\t\n' for i, gene in enumerate(ids)))
    output = tmp_path / 'result.tsv'
    subprocess.run([sys.executable, str(SCRIPT_PATH), '--cds_fasta', str(fasta), '--uniprot_tsv', str(uniprot),
                    '--orthogroup_tsv', str(orthogroups), '--scientific_name', 'Species',
                    '--expression_tsv', str(expression), '--fx2tab', str(fx2tab), '--busco_tsv', str(busco),
                    '--out_tsv', str(output), '--ncpu', str(ncpu)], check=True, capture_output=True, text=True)
    result = pandas.read_csv(output, sep='\t', keep_default_na=False, dtype={'gene_id': str})
    assert result.gene_id.tolist() == ids
    assert result.orthogroup.tolist() == ids
    assert result.uniprot_label.tolist() == ['label3', 'label2', 'label1', 'label0']
    assert result['Identifier_score'].tolist() == [7] * 4
    assert result.cds_length.tolist() == [99] * 4
    assert result.length_ratio.tolist() == [0.5] * 4
    assert result.busco_id.tolist() == ['B0', 'B1', 'B2', 'B3']
    assert result.busco_description.tolist() == ['nan'] * 4
