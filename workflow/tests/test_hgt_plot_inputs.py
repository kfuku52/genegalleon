import hashlib
import os
import sqlite3
import subprocess
import sys
from pathlib import Path

import pandas
import pytest

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / 'workflow/support'
sys.path.insert(0, str(SUPPORT))
import hgt_plot_inputs as index  # noqa: E402


def candidates(tmp_path):
    branch = tmp_path / 'branch.tsv'
    gene = tmp_path / 'gene.tsv'
    branch.write_text('orthogroup\tbranch_id\tcandidate_gene_count\toptional\n'
                      'OG1\t1\t2\t001\nOG2\t3\t1\tlate-text\n')
    gene.write_text('orthogroup\tgene_id\tgene_taxon\tcandidate_branch_count\tintron_supported\n'
                    'OG1\ta\tZ species\t1\tTrue\nOG2\tb\tNA\t2\t\n'
                    'OG1\tc\tA species\t3\tFalse\n')
    return branch, gene


def test_family_frames_preserve_whole_table_types_rows_and_prefixes(tmp_path):
    branch, gene = candidates(tmp_path)
    db = tmp_path / 'index.db'
    index.prepare(db, branch, gene)
    for family, prefix in [('OG1', 'A_species_'), ('OG2', 'NA_'), ('missing', '')]:
        actual = index.read_family(db, branch, gene, family)
        for frame, path in zip(actual[:2], [branch, gene], strict=True):
            expected = pandas.read_csv(path, sep='\t', low_memory=False)
            expected = expected.loc[expected.orthogroup.astype(str) == family]
            pandas.testing.assert_frame_equal(frame, expected)
        assert actual[2] == prefix


def test_unchanged_index_reuse_does_not_parse_candidate_tables(tmp_path, monkeypatch):
    branch, gene = candidates(tmp_path)
    db = tmp_path / 'index.db'
    index.prepare(db, branch, gene)
    before = db.read_bytes()
    def reject(*args, **kwargs):
        raise AssertionError('whole candidate TSV reparsed during warm reuse')
    monkeypatch.setattr(index.pandas, 'read_csv', reject)
    index.prepare(db, branch, gene)
    assert db.read_bytes() == before
    assert len(index.read_family(db, branch, gene, 'OG1')[0]) == 1


def test_same_size_source_replacement_rejects_old_index_then_rebuilds(tmp_path):
    branch, gene = candidates(tmp_path)
    db = tmp_path / 'index.db'
    index.prepare(db, branch, gene)
    before = gene.stat()
    gene.write_text(gene.read_text().replace('Z species', 'B species'))
    os.utime(gene, ns=(before.st_atime_ns, before.st_mtime_ns))
    with pytest.raises(ValueError, match='changed since index'):
        index.read_family(db, branch, gene, 'OG1')
    index.prepare(db, branch, gene)
    assert 'B species' in index.read_family(db, branch, gene, 'OG1')[1].gene_taxon.tolist()


def test_corrupt_spool_payload_is_rejected_before_deserialization(tmp_path):
    branch, gene = candidates(tmp_path)
    db = tmp_path / 'index.db'
    index.prepare(db, branch, gene)
    with sqlite3.connect(db) as connection:
        connection.execute("UPDATE frames SET payload=? WHERE family='OG1'", (b'not-json',))
    with pytest.raises(ValueError, match='Corrupt private'):
        index.read_family(db, branch, gene, 'OG1')


def test_invalid_input_cannot_replace_an_existing_index(tmp_path):
    branch, gene = candidates(tmp_path)
    db = tmp_path / 'index.db'
    index.prepare(db, branch, gene)
    before = hashlib.sha256(db.read_bytes()).hexdigest()
    branch.write_text('wrong\tcolumns\n1\t2\n')
    with pytest.raises(ValueError, match='orthogroup column'):
        index.prepare(db, branch, gene)
    assert hashlib.sha256(db.read_bytes()).hexdigest() == before
    link = tmp_path / 'unsafe.tsv'
    link.symlink_to(gene)
    with pytest.raises(ValueError, match='owned regular'):
        index.prepare(db, branch, link)


def test_indexed_annotation_matches_original_cli_bytes(tmp_path):
    branch, gene = candidates(tmp_path)
    db = tmp_path / 'index.db'
    index.prepare(db, branch, gene)
    source = tmp_path / 'stat.tsv'
    source.write_text('branch_id\tnode_name\tso_event\n1\ta\tL\n2\tc\tL\n')
    args = [sys.executable, str(SUPPORT / 'annotate_hgt_tree_plot.py'), '--stat_branch', str(source),
            '--branch_tsv', str(branch), '--gene_tsv', str(gene), '--orthogroup', 'OG1']
    original, indexed = tmp_path / 'original.tsv', tmp_path / 'indexed.tsv'
    subprocess.run([*args, '--outfile', str(original)], check=True, capture_output=True)
    prefix = tmp_path / 'prefix.txt'
    subprocess.run([*args, '--candidate_index', str(db), '--ortholog_prefix_out', str(prefix),
                    '--outfile', str(indexed)], check=True, capture_output=True)
    assert indexed.read_bytes() == original.read_bytes()
    assert prefix.read_text() == 'A_species_\n'


def shell_functions():
    source = (ROOT / 'workflow/core/gg_hgt_core.sh').read_text()
    return source[source.index('hgt_plot_worker_limit() {'):source.index('\nif [[ -s "${file_hgt_branch}"', source.index('hgt_plot_worker_limit() {'))]


@pytest.mark.parametrize('cpus,memory,requested,expected', [(4, 60, 'auto', 4), (16, 60, 'auto', 7),
    (8, 16, 'auto', 2), (4, 60, '2', 2), (1, 4, 'auto', 1)])
def test_worker_limit_respects_cpus_and_memory(cpus, memory, requested, expected):
    result = subprocess.run(['bash', '-c', shell_functions() + '\nhgt_plot_worker_limit'],
        env={**os.environ, 'GG_TASK_CPUS': str(cpus), 'GG_MEM_TOOL_GB': str(memory),
             'hgt_plot_workers': requested}, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == str(expected)


@pytest.mark.parametrize('requested', ['0', '-1', 'five', '5', '01'])
def test_worker_limit_rejects_invalid_or_overallocated_settings(requested):
    result = subprocess.run(['bash', '-c', shell_functions() + '\nhgt_plot_worker_limit'],
        env={**os.environ, 'GG_TASK_CPUS': '4', 'GG_MEM_TOOL_GB': '60', 'hgt_plot_workers': requested},
        capture_output=True, text=True)
    assert result.returncode != 0


def test_parallel_pool_reaps_failures_and_does_not_serialize_fast_workers(tmp_path):
    source = (ROOT / 'workflow/core/gg_hgt_core.sh').read_text()
    start = source.index('  hgt_plot_pids=()\n')
    end = source.index('\nfi\n\nhgt_cleanup_materialization_run', start)
    pool = source[start:end]
    script = shell_functions() + '''
hgt_plot_family() { sleep 0.15; echo "$1" >> "$RESULT"; [[ "$1" != bad ]]; }
hgt_plot_limit=3
hgt_orthogroups=(one bad three four five)
''' + pool
    result = subprocess.run(['bash', '-c', script], env={**os.environ, 'RESULT': str(tmp_path / 'done')},
                            capture_output=True, text=True, timeout=10)
    assert result.returncode != 0
    done = (tmp_path / 'done').read_text().splitlines()
    assert {'one', 'bad', 'three'}.issubset(done)
    assert len(done) == len(set(done))
    assert 'five' not in done


def test_typed_index_roundtrip_retains_nonfinite_signs_and_object_boolean_nulls():
    import numpy
    frame = pandas.DataFrame({'floating': [numpy.inf, -numpy.inf, -0.0, numpy.nan],
                              'mixed': pandas.Series([True, False, numpy.nan, True], dtype=object),
                              'integer': [2**60, 2**60 + 1, 2**60 + 2, 2**60 + 3]})
    actual = index.decode_frame(index.encode_frame(frame))
    pandas.testing.assert_frame_equal(actual, frame)
    assert numpy.signbit(actual.floating.iloc[2])
