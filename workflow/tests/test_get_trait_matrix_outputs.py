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


@pytest.mark.parametrize('header', ['gene_id\tvalue\tmissing\tlabel\n', 'value\tmissing\tlabel\n'])
def test_shared_index_preserves_types_duplicate_order_and_literal_ids(tmp_path, header):
    import pandas as pd

    from workflow.support.get_trait_matrix import process_trait_file

    path = tmp_path / 'traits.tsv'
    path.write_text(header + '001\t1.2345678901234567\tNA\ttext\nNA\t-2\tNULL\t001\n'
                    '001\t3\tnan\tNA\nNULL\t4\tNA\tother\n')
    mapping = {'001': 'Species_a_001', 'NA': 'Species_a_NA', 'NULL': 'Species_a_NULL'}
    for ids in [set(mapping), {'001'}, {'absent'}, {'NA', 'NULL'}]:
        expected = process_trait_file(path, ids, mapping)
        for _ in range(2):
            actual = process_trait_file(path, ids, mapping, tmp_path / 'cache')
            pd.testing.assert_frame_equal(actual, expected)
    path.write_text(header + '001\t9\tNA\tnew\n')
    pd.testing.assert_frame_equal(process_trait_file(path, set(mapping), mapping, tmp_path / 'cache'),
                                  process_trait_file(path, set(mapping), mapping))


def test_shared_index_build_once_and_rebuild_corrupt_index(tmp_path):
    import pandas as pd

    from workflow.support.trait_table_index import read_subset

    path = tmp_path / 'traits.tsv'
    path.write_text('gene_id\tvalue\ng\t1\ng\t2\n')
    calls = []
    def loader(path):
        calls.append(path)
        return pd.read_csv(path, sep='\t')
    cache = tmp_path / 'cache'
    for _ in range(2):
        result = read_subset(path, {'g'}, loader, cache)
        assert result.value.tolist() == [1, 2]
    assert len(calls) == 1
    next(cache.glob('*.sqlite3')).write_bytes(b'corrupt index')
    assert read_subset(path, {'g'}, loader, cache).value.tolist() == [1, 2]
    assert len(calls) == 2


@pytest.mark.parametrize('change', ['rewrite', 'replace'])
def test_shared_index_rejects_source_change_during_parse(tmp_path, change):
    import pandas as pd

    from workflow.support.trait_table_index import read_subset

    path = tmp_path / 'traits.tsv'
    path.write_text('gene_id\tvalue\ng\t1\n')
    def loader(path):
        result = pd.read_csv(path, sep='\t')
        if change == 'rewrite':
            path.write_text('gene_id\tvalue\ng\t2\n')
        else:
            replacement = tmp_path / 'replacement'
            replacement.write_bytes(path.read_bytes())
            replacement.replace(path)
        return result
    with pytest.raises(OSError, match='changed'):
        read_subset(path, {'g'}, loader, tmp_path / 'cache')


def test_concurrent_family_queries_publish_one_complete_species_index(tmp_path):
    from concurrent.futures import ThreadPoolExecutor

    import pandas as pd

    from workflow.support.get_trait_matrix import process_trait_file

    path = tmp_path / 'traits.tsv'
    path.write_text('gene_id\tvalue\n' + ''.join(f'g{i}\t{i}\n' for i in range(100)))
    mapping = {f'g{i}': f'Species_a_g{i}' for i in range(100)}
    def query(i):
        selected = {f'g{i}', f'g{i+1}'}
        return process_trait_file(path, selected, mapping, tmp_path / 'cache')
    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(query, range(8)))
    for i, actual in enumerate(results):
        pd.testing.assert_frame_equal(actual, process_trait_file(path, {f'g{i}', f'g{i+1}'}, mapping))
    assert len(list((tmp_path / 'cache').glob('*.sqlite3'))) == 1


@pytest.mark.parametrize('compressed', [False, True])
def test_indexed_wide_rows_preserve_unicode_unsigned_booleans_and_float_values(tmp_path, compressed):
    import pandas as pd

    from workflow.support.get_trait_matrix import process_trait_file

    path = tmp_path / ('traits.tsv.gz' if compressed else 'traits.tsv')
    pd.DataFrame({'gene_id': ['001', 'NA', '001'], 'text': ['λ'*2000, 'NULL', 'quoted\twith\nnewline'],
                  'unsigned': [2**64-2, 2**64-3, 2**64-4], 'flag': [True, False, True],
                  'value': [float('nan'), float('inf'), -0.0]}).to_csv(path, sep='\t', index=False)
    mapping = {'001': 'Species_a_001', 'NA': 'Species_a_NA'}
    expected = process_trait_file(path, set(mapping), mapping)
    for _ in range(2):
        actual = process_trait_file(path, set(mapping), mapping, tmp_path / 'cache')
        pd.testing.assert_frame_equal(actual, expected)


@pytest.mark.parametrize('compressed', [False, True])
def test_index_binding_preserves_filename_based_compression_errors(tmp_path, compressed):
    import gzip

    from workflow.support.get_trait_matrix import process_trait_file

    source = tmp_path / ('traits.tsv.gz' if compressed else 'traits.tsv')
    content = b'gene_id\tvalue\ng\t1\n'
    source.write_bytes(gzip.compress(content) if compressed else content)
    alias = tmp_path / ('alias.tsv' if compressed else 'alias.tsv.gz')
    alias.symlink_to(source)
    mapping = {'g': 'Species_a_g'}
    process_trait_file(source, {'g'}, mapping, tmp_path / 'cache')
    with pytest.raises((UnicodeError, OSError)) as expected:
        process_trait_file(alias, {'g'}, mapping)
    with pytest.raises(type(expected.value)):
        process_trait_file(alias, {'g'}, mapping, tmp_path / 'cache')


@pytest.mark.parametrize('target', ['source', 'index'])
@pytest.mark.parametrize('change', ['rewrite', 'replace'])
def test_shared_index_rejects_mutation_while_decoding_cached_rows(tmp_path, monkeypatch, target, change):
    import pandas as pd

    from workflow.support import trait_table_index as indexer

    source = tmp_path / 'traits.tsv'
    source.write_text('gene_id\tvalue\ng0\t1\ng1\t2\ng0\t3\n')
    cache = tmp_path / 'cache'

    def loader(path):
        return pd.read_csv(path, sep='\t')

    indexer.read_subset(source, {'g0'}, loader, cache)
    database = next(cache.glob('*.sqlite3'))
    mutated = []
    original_loads = indexer.json.loads

    def loads_then_mutate(payload, *args, **kwargs):
        value = original_loads(payload, *args, **kwargs)
        if not mutated and isinstance(value, list) and value and value[0] == 'g0':
            path = source if target == 'source' else database
            if change == 'rewrite':
                # Preserve valid bytes so the final proof, rather than a parser
                # or SQLite error, is responsible for rejecting the mutation.
                with path.open('r+b') as handle:
                    first = handle.read(1)
                    handle.seek(0)
                    handle.write(first)
            else:
                replacement = path.with_name(path.name + '.replacement')
                replacement.write_bytes(path.read_bytes())
                replacement.replace(path)
            mutated.append(path)
        return value

    monkeypatch.setattr(indexer.json, 'loads', loads_then_mutate)
    with pytest.raises(OSError, match='changed'):
        indexer.read_subset(source, {'g0'}, loader, cache, source_proof=indexer.FreshDigestBatch())
    assert mutated == [source if target == 'source' else database]


@pytest.mark.parametrize('index_kind', ['integer', 'string'])
def test_shared_index_keeps_extension_types_with_nonunique_named_index(tmp_path, index_kind):
    import pandas as pd

    from workflow.support.trait_table_index import read_subset

    source = tmp_path / 'traits.tsv'
    source.write_text('gene_id\textension_float\tlabel\tmissing\n'
                      '001\t1.25\tλ\t\nNA\t-2.5\tNULL\t3.5\n'
                      '001\t0\t日本語\t\nother\t4\tkept\t5\n')

    def loader(path):
        frame = pd.read_csv(path, sep='\t', keep_default_na=False, dtype={'gene_id': str})
        # These extension values are JSON-representable; missing traits use the
        # existing ordinary float/NaN contract of the index format.
        frame['extension_float'] = frame['extension_float'].astype('Float64')
        frame['label'] = frame['label'].astype('string')
        frame['missing'] = pd.to_numeric(frame['missing'], errors='coerce')
        labels = [9, 2, 9, 4] if index_kind == 'integer' else ['row9', 'row2', 'row9', 'row4']
        frame.index = pd.Index(labels, name='original_row')
        return frame

    cache = tmp_path / 'cache'
    for identifiers in ({'001'}, {'NA'}, {'001', 'NA'}, {'absent'}, set()):
        expected = read_subset(source, identifiers, loader)
        for _ in range(2):
            actual = read_subset(source, identifiers, loader, cache)
            pd.testing.assert_frame_equal(actual, expected)
