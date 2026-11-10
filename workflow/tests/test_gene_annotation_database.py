import csv
import json
import os
import shlex
import sqlite3
import subprocess
from pathlib import Path

import pytest

from workflow.support import gene_annotation_database as annotation
from workflow.support.reconciliation_summary import build_results
from workflow.tests.shell_static_helpers import function_body
from workflow.tests.test_reconciliation_database import build_fixture, config_fixture


def write_tsv(path, fields, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


@pytest.fixture
def cohort(tmp_path):
    db, *_ = build_fixture(tmp_path)
    with sqlite3.connect(db) as conn:
        meta = dict(conn.execute('SELECT key,value FROM reconciliation_meta'))
        genes = list(conn.execute('SELECT gene_id,orthogroup,species FROM reconciliation_gene'))
    fields = ['gene_id', 'orthogroup', 'sprot_recname', *annotation.GO_COLUMNS.values(),
              'cdskit_localize_p_extracellular', 'cdskit_localize_score_kind',
              'cdskit_localize_calibration_status']
    sources, by_species = [], {}
    for species in ('A', 'B', 'C'):
        rows = []
        for gene, family, taxon in genes:
            if taxon == species:
                rows.append(dict(gene_id=gene, orthogroup=family, sprot_recname='Extant enzyme\nsubunit',
                                 go_ids='GO:0000001;GO:0000002', go_aspects='P;F',
                                 go_terms='process; function', go_evidence='IDA;IEA',
                                 cdskit_localize_p_extracellular='0.123456789',
                                 cdskit_localize_score_kind='sigmoid_score',
                                 cdskit_localize_calibration_status='unverified'))
        if species == 'A':
            rows.append(dict.fromkeys(fields, '') | {'gene_id': 'outside_selected', 'orthogroup': 'OG9999'})
        path = tmp_path / (species + '.tsv')
        write_tsv(path, fields, rows)
        sources.append({'species': species, 'path': path.name, 'sha256': annotation.sha256(path)})
        by_species[species] = rows
    manifest = tmp_path / 'annotations.json'
    config = {'schema': annotation.SCHEMA, 'gene_model_contract_sha256': meta['gene_model_contract_sha256'],
              'family_catalog_sha256': meta['family_catalog_sha256'], 'source_database_sha256': annotation.sha256(db),
              'sources': sources}
    manifest.write_text(json.dumps(config))
    return {'db': db, 'manifest': manifest, 'config': config, 'rows': by_species,
            'fields': fields, 'target': tmp_path / 'annotated.db'}


def update_source(cohort, species, rows=None, fields=None):
    item = next(x for x in cohort['config']['sources'] if x['species'] == species)
    path = cohort['manifest'].parent / item['path']
    write_tsv(path, fields or cohort['fields'], rows if rows is not None else cohort['rows'][species])
    item['sha256'] = annotation.sha256(path)
    cohort['manifest'].write_text(json.dumps(cohort['config']))


def build(cohort, **kwargs):
    return annotation.build_database(cohort['db'], cohort['manifest'], cohort['target'], **kwargs)


def test_raw_catalogue_resolution_merge_GO_and_database_keep_excluded_genes(cohort, tmp_path):
    import sys

    from workflow.support import cds_resolution
    from workflow.support.cds_annotation_scope import QC_COLUMNS
    from workflow.support.reformat_uniprot_diamond import refresh_go_annotations
    helper = Path(__file__).resolve().parents[1] / 'support/merge_cds_annotation.py'
    metadata = tmp_path / 'uniprot-meta.tsv'
    metadata.write_text('accession\tgo_ids\tgo_aspects\tgo_terms\tgo_evidence\nP12345\tGO:0000001\tP\tprocess\tIDA\n')
    totals = {'accepted': 0, 'excluded': 0}
    for species in ('A', 'B', 'C'):
        work = tmp_path / species
        work.mkdir()
        rows = cohort['rows'][species]
        raw = work / 'catalogue.fa'
        sequences = {row['gene_id']: 'ATGAAATAA' for row in rows}
        if species == 'A':
            # Keep a selected reconciliation gene even though its supplied CDS fails QC.
            sequences[rows[0]['gene_id']] = 'ATGNNNTAA'
            sequences.update(stop='TAATAAATAAATAAATGAAA', frame='ATGTAATAA')
        raw.write_text(''.join(f'>{gene}\n{seq}\n' for gene, seq in sequences.items()))
        resolved, _, report = cds_resolution.resolve(raw, None, None, work / 'resolved')
        decisions = json.loads(report.read_text())['decisions']
        accepted = {v['gene_id'] for v in decisions if v['selected_source'] is not None}
        saved = work / 'saved.tsv'
        fields = ['gene_id', 'orthogroup', 'sprot_recname', 'sprot_best', 'exact_score', *annotation.GO_COLUMNS.values()]
        saved_rows = []
        for row in rows:
            if row['gene_id'] in accepted:
                saved_rows.append(dict.fromkeys(fields, '') | dict(gene_id=row['gene_id'], orthogroup=row['orthogroup'],
                    sprot_recname='literal product' if row['orthogroup'] else '',
                    sprot_best='P12345' if row['orthogroup'] else '', exact_score='001.2300'))
        write_tsv(saved, fields, saved_rows)
        og = work / 'orthogroups.tsv'
        og.write_text('Orthogroup\t' + species + '\n' + ''.join(f'{row["orthogroup"]}\t{row["gene_id"]}\n' for row in rows if row['orthogroup']))
        merged, scope = work / 'merged.tsv', work / 'scope.json'
        subprocess.run([sys.executable, str(helper), '--scientific_name', species, '--cds_fasta', str(raw),
            '--resolved_cds_fasta', str(resolved), '--cds_resolution_json', str(report), '--scope_json', str(scope),
            '--existing_annotation_tsv', str(saved), '--orthogroup_tsv', str(og), '--out_tsv', str(merged)], check=True)
        refresh_go_annotations([str(merged)], str(metadata), str(work / 'GO'))
        result = work / 'GO/merged.tsv'
        with result.open() as handle:
            actual = list(csv.DictReader(handle, delimiter='\t'))
        assert {r['gene_id'] for r in actual} == set(sequences)
        for row in actual:
            assert set(QC_COLUMNS) <= row.keys()
            totals[row['cds_qc_status']] += 1
            if row['cds_qc_status'] == 'excluded':
                assert row['sprot_recname'] == row['go_ids'] == row['exact_score'] == ''
                assert row['functional_annotation_status'] == 'not_assessed_cds_qc_excluded'
                assert json.loads(row['cds_qc_exclusion_reasons'])
            else:
                assert row['exact_score'] == '001.2300'
                assert row['sprot_recname'] == ('literal product' if row['orthogroup'] else '')
        source = next(v for v in cohort['config']['sources'] if v['species'] == species)
        source.update(path=str(result), sha256=annotation.sha256(result), cds_scope=json.loads(scope.read_text()))
    cohort['manifest'].write_text(json.dumps(cohort['config']))
    audit = build(cohort)
    assert audit['genes'] == 10 and audit['reconciliation_genes_without_annotation'] == 0
    assert audit['cds_qc_status_counts'] == totals
    with sqlite3.connect(cohort['target']) as db:
        excluded = db.execute("SELECT gene_id,cds_qc_exclusion_reasons FROM gene_annotation WHERE cds_qc_status='excluded'").fetchall()
        assert len(excluded) == 3
        assert any('internal_ambiguous_bases' in reason for _, reason in excluded)
        assert any('premature_stop' in reason for _, reason in excluded)
        assert any('ambiguous_frame' in reason for _, reason in excluded)
        assert db.execute('SELECT count(*) FROM gene_annotation_source WHERE cds_scope_json!=\'null\'').fetchone() == (3,)


def test_missing_scope_fails_before_large_database_copy(cohort, monkeypatch):
    cohort['rows']['A'].pop(0)
    update_source(cohort, 'A')
    def unexpected_copy(*args, **kwargs):
        raise AssertionError('Large DB staging must not start for missing annotation coverage')
    monkeypatch.setattr(annotation.tempfile, 'TemporaryDirectory', unexpected_copy)
    with pytest.raises(ValueError, match='Annotations omit 1.*by species'):
        build(cohort)


@pytest.mark.parametrize('fault', ['functional_hit', 'no_reason', 'accepted_candidate', 'wrong_scope'])
def test_excluded_qc_evidence_is_strictly_validated(cohort, fault):
    from workflow.support.cds_annotation_scope import QC_COLUMNS
    fields = [*cohort['fields'], *QC_COLUMNS]
    rows = []
    for original in cohort['rows']['A']:
        row = dict(original)
        row.update(dict(zip(QC_COLUMNS, ('excluded', 'no_valid_cds', '{"supplied":"premature_stop"}',
            '{"supplied":{"accepted":false,"reason":"premature_stop"}}', '', 'not_assessed_cds_qc_excluded'), strict=True)))
        for key in fields:
            if key not in {'gene_id', 'orthogroup', *QC_COLUMNS}:
                row[key] = ''
        rows.append(row)
    if fault == 'functional_hit':
        rows[0]['sprot_recname'] = 'invented product'
    elif fault == 'no_reason':
        rows[0]['cds_qc_exclusion_reasons'] = '{}'
    elif fault == 'accepted_candidate':
        rows[0]['cds_qc_candidates_json'] = '{"supplied":{"accepted":true,"reason":"accepted"}}'
    else:
        cohort['config']['sources'][0]['cds_scope'] = {'schema': 'annotation_cds_scope_v1', 'input_count': 0}
    update_source(cohort, 'A', rows, fields)
    with pytest.raises(ValueError):
        build(cohort)
    assert not cohort['target'].exists()


def test_copies_all_tables_and_keeps_per_gene_GO_localization_and_source(cohort):
    before = annotation.sha256(cohort['db'])
    audit = build(cohort)
    assert audit['genes'] == 8 and audit['GO_tuples'] == 14
    assert audit['reconciliation_genes_without_annotation'] == 0
    assert annotation.sha256(cohort['db']) == before
    with sqlite3.connect(cohort['db']) as original, sqlite3.connect(cohort['target']) as derived:
        for (table,) in original.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall():
            assert original.execute('SELECT * FROM ' + table).fetchall() == derived.execute('SELECT * FROM ' + table).fetchall()
        assert derived.execute("SELECT orthogroup,source_orthogroup FROM gene_annotation WHERE gene_id='outside_selected'").fetchone() == (None, 'OG9999')
        assert derived.execute("SELECT product FROM gene_annotation_product WHERE gene_id='OG0001_A'").fetchone() == ('Extant enzyme\nsubunit',)
        assert derived.execute("SELECT go_id,aspect,term,evidence FROM gene_annotation_go WHERE gene_id='OG0001_A' ORDER BY ordinal").fetchall() == [
            ('GO:0000001', 'P', 'process', 'IDA'), ('GO:0000002', 'F', 'function', 'IEA')]
        assert derived.execute("SELECT cdskit_localize_p_extracellular,cdskit_localize_score_kind,cdskit_localize_calibration_status FROM gene_annotation WHERE gene_id='OG0001_A'").fetchone() == ('0.123456789', 'sigmoid_score', 'unverified')
        assert derived.execute('SELECT species,sha256 FROM gene_annotation_source ORDER BY species').fetchall() == [
            (x['species'], x['sha256']) for x in cohort['config']['sources']]
    assert cohort['target'].stat().st_mode & 0o222 == 0


def test_embedded_products_used_by_existing_reconciliation_downstream(cohort, tmp_path):
    build(cohort)
    config, _ = config_fixture(tmp_path, cohort['target'])
    before = annotation.sha256(cohort['target'])
    build_results(cohort['target'], config, tmp_path / 'summary')
    with sqlite3.connect(tmp_path / 'summary/reconciliation_summary.db') as db:
        row = db.execute("SELECT products FROM branch_family WHERE branch_id='loss_B' AND orthogroup='OG0001'").fetchone()
        assert json.loads(row[0]) == ['Extant enzyme\nsubunit']
    assert annotation.sha256(cohort['target']) == before


@pytest.mark.parametrize('key', ['gene_model_contract_sha256', 'family_catalog_sha256', 'source_database_sha256'])
def test_incompatible_contract_preserves_previous_output(cohort, key):
    cohort['target'].write_bytes(b'previous publication')
    cohort['config'][key] = '0' * 64
    cohort['manifest'].write_text(json.dumps(cohort['config']))
    with pytest.raises(ValueError):
        build(cohort)
    assert cohort['target'].read_bytes() == b'previous publication'


@pytest.mark.parametrize('fault', ['duplicate_gene', 'wrong_family', 'missing_gene', 'blank_gene',
                                  'GO_count', 'GO_id', 'GO_aspect', 'source_sha'])
def test_invalid_annotations_fail_before_publication(cohort, fault):
    rows = cohort['rows']['A']
    if fault == 'duplicate_gene':
        rows.append(rows[0].copy())
    elif fault == 'wrong_family':
        rows[0]['orthogroup'] = 'OG0002' if rows[0]['orthogroup'] != 'OG0002' else 'OG0001'
    elif fault == 'missing_gene':
        rows.pop(0)
    elif fault == 'blank_gene':
        rows[0]['gene_id'] = ''
    elif fault == 'GO_count':
        rows[0]['go_evidence'] = 'IDA'
    elif fault == 'GO_id':
        rows[0]['go_ids'] = 'bad;GO:0000002'
    elif fault == 'GO_aspect':
        rows[0]['go_aspects'] = 'UNKNOWN;F'
    update_source(cohort, 'A')
    if fault == 'source_sha':
        cohort['config']['sources'][0]['sha256'] = '0' * 64
        cohort['manifest'].write_text(json.dumps(cohort['config']))
    with pytest.raises((ValueError, sqlite3.IntegrityError)):
        build(cohort)
    assert not cohort['target'].exists()


def test_explicit_partial_annotation_records_missing_genes(cohort):
    cohort['config']['require_complete_genes'] = False
    cohort['rows']['A'].pop(0)
    update_source(cohort, 'A')
    audit = build(cohort)
    assert audit['reconciliation_genes_without_annotation'] == 1


def test_species_specific_missing_columns_remain_SQL_NULL(cohort):
    fields = [x for x in cohort['fields'] if not x.startswith('cdskit_')]
    rows = [{key: value for key, value in row.items() if key in fields} for row in cohort['rows']['B']]
    update_source(cohort, 'B', rows, fields)
    build(cohort)
    with sqlite3.connect(cohort['target']) as db:
        assert db.execute("SELECT cdskit_localize_score_kind FROM gene_annotation WHERE species='B'").fetchone() == (None,)


def test_same_literal_gene_id_in_different_species_is_not_collapsed(cohort):
    for species in ('A', 'B'):
        cohort['rows'][species].append(dict.fromkeys(cohort['fields'], '') | {'gene_id': 'same', 'orthogroup': ''})
        update_source(cohort, species)
    build(cohort)
    with sqlite3.connect(cohort['target']) as db:
        assert db.execute("SELECT species FROM gene_annotation WHERE gene_id='same' ORDER BY species").fetchall() == [('A',), ('B',)]


def test_reconciled_gene_id_in_another_species_keeps_its_own_annotation(cohort):
    gene = cohort['rows']['A'][0]['gene_id']
    cohort['rows']['B'].append(dict.fromkeys(cohort['fields'], '') | {'gene_id': gene, 'orthogroup': 'OG9999'})
    update_source(cohort, 'B')
    build(cohort)
    with sqlite3.connect(cohort['target']) as db:
        assert db.execute(
            "SELECT orthogroup,source_orthogroup FROM gene_annotation WHERE species='B' AND gene_id=?", (gene,)
        ).fetchone() == (None, 'OG9999')


def test_cross_method_gene_conflict_fails_before_annotation_publication(cohort):
    with sqlite3.connect(cohort['db']) as db:
        gene, species = db.execute("SELECT gene_id,species FROM reconciliation_gene LIMIT 1").fetchone()
        db.execute("CREATE TABLE species_reconciliation_gene (gene_id TEXT, species TEXT, orthogroup TEXT)")
        db.execute("INSERT INTO species_reconciliation_gene VALUES (?,?,?)", (gene, species, 'OG9999'))
    cohort['config']['source_database_sha256'] = annotation.sha256(cohort['db'])
    cohort['manifest'].write_text(json.dumps(cohort['config']))
    with pytest.raises(ValueError, match="Reconciliation methods disagree on gene identity"):
        build(cohort)
    assert not cohort['target'].exists()


@pytest.mark.parametrize('suffix', ['-wal', '-journal'])
def test_active_database_not_used(cohort, suffix):
    (cohort['db'].parent / (cohort['db'].name + suffix)).write_bytes(b'active')
    with pytest.raises(ValueError, match='closed published'):
        build(cohort)


def test_no_in_place_mutation_or_nested_annotation(cohort):
    before = annotation.sha256(cohort['db'])
    with pytest.raises(ValueError, match='separate database'):
        annotation.build_database(cohort['db'], cohort['manifest'], cohort['db'])
    assert annotation.sha256(cohort['db']) == before
    build(cohort)
    cohort['config'].pop('source_database_sha256')
    cohort['manifest'].write_text(json.dumps(cohort['config']))
    with pytest.raises(ValueError, match='original family DB'):
        annotation.build_database(cohort['target'], cohort['manifest'], cohort['target'].with_name('nested.db'))


def test_input_change_during_import_preserves_output(cohort, monkeypatch):
    original = annotation.import_annotations
    def changed(*args):
        result = original(*args)
        cohort['manifest'].write_text(cohort['manifest'].read_text() + ' ')
        return result
    monkeypatch.setattr(annotation, 'import_annotations', changed)
    cohort['target'].write_bytes(b'previous')
    with pytest.raises(ValueError, match='input changed'):
        build(cohort)
    assert cohort['target'].read_bytes() == b'previous'


def test_conflicting_external_product_rejected(cohort, tmp_path):
    build(cohort)
    config_path, config = config_fixture(tmp_path, cohort['target'])
    product = tmp_path / 'products.tsv'
    write_tsv(product, ['gene_id', 'product'], [{'gene_id': 'OG0001_A', 'product': 'Conflicting product'}])
    config['inputs']['gene_annotations'] = {'path': product.name, 'sha256': annotation.sha256(product)}
    config_path.write_text(json.dumps(config))
    with pytest.raises(ValueError, match='conflicts'):
        build_results(cohort['target'], config_path, tmp_path / 'summary')


@pytest.mark.parametrize('field', ['Species', 'annotation_source_id', 'source_orthogroup'])
def test_metadata_column_collisions_fail(cohort, field):
    rows = [dict(row, **{field: 'value'}) for row in cohort['rows']['A']]
    update_source(cohort, 'A', rows, cohort['fields'] + [field])
    with pytest.raises(ValueError, match='reserved annotation columns'):
        build(cohort)


def test_long_native_structure_field_is_preserved(cohort):
    rows = [dict(row, feature_block_sequences='A' * 200000) for row in cohort['rows']['A']]
    update_source(cohort, 'A', rows, cohort['fields'] + ['feature_block_sequences'])
    build(cohort)
    with sqlite3.connect(cohort['target']) as db:
        assert db.execute("SELECT length(feature_block_sequences) FROM gene_annotation WHERE gene_id='OG0001_A'").fetchone() == (200000,)


def test_repeated_go_parsing_is_scoped_to_one_import(cohort, monkeypatch):
    original, seen = annotation._parse_go_values, []

    def counted(values):
        seen.append(values)
        return original(values)

    monkeypatch.setattr(annotation, '_parse_go_values', counted)
    build(cohort)
    assert len(seen) == 2
    cohort['target'] = cohort['target'].with_name('second.db')
    build(cohort)
    assert len(seen) == 4  # No cross-command retention of parsed input.


def test_large_go_keys_bypass_cache_without_truncating_terms(cohort, monkeypatch):
    for species, rows in cohort['rows'].items():
        for row in rows:
            if row['go_ids']:
                row['go_terms'] = ('x' * 3000) + ';' + ('y' * 3000)
        update_source(cohort, species)
    original, seen = annotation._parse_go_values, []

    def counted(values):
        seen.append(values)
        return original(values)

    monkeypatch.setattr(annotation, '_parse_go_values', counted)
    audit = build(cohort)
    assert len(seen) == audit['genes']
    with sqlite3.connect(cohort['target']) as db:
        assert db.execute('SELECT max(length(term)) FROM gene_annotation_go').fetchone() == (3000,)


@pytest.mark.parametrize('field,value', [('go_aspects', 'UNKNOWN;F'), ('go_evidence', 'IDA')])
def test_valid_cached_ids_do_not_hide_changed_invalid_tuple(cohort, field, value):
    cohort['rows']['C'][-1][field] = value
    update_source(cohort, 'C')
    with pytest.raises(ValueError):
        build(cohort)
    assert not cohort['target'].exists()


def test_go_tuple_result_is_mutable_without_affecting_later_calls():
    row = dict(go_ids='GO:0000001', go_aspects='P', go_terms='literal', go_evidence='IDA')
    first = annotation.go_tuples(row, dict(reversed(list(annotation.GO_COLUMNS.items()))))
    first.clear()
    assert annotation.go_tuples(row, annotation.GO_COLUMNS) == [(0, 'GO:0000001', 'P', 'literal', 'IDA')]


def test_native_annotation_launch_under_nounset_uses_runtime_scratch(cohort, tmp_path):
    workflow = Path(__file__).resolve().parents[1]
    function = function_body((workflow / 'core/gg_gene_summary_core.sh').read_text(),
                             'run_gene_annotation_database_for_source')
    family, summary, scratch = (tmp_path / name for name in ('native-family', 'summary', 'runtime-scratch'))
    family.mkdir()
    (summary / 'artifact_provenance').mkdir(parents=True)
    scratch.mkdir(mode=0o700)
    values = dict(run_gene_annotation_database='1', dir_gene_family=str(family),
                  summary_output_dir=str(summary), gene_annotation_manifest=str(cohort['manifest']),
                  file_gene_family_db=str(cohort['db']), gg_support_dir=str(workflow / 'support'),
                  gene_family_source='orthogroup', gg_workspace_output_dir=str(tmp_path / 'output'),
                  gg_workspace_dir=str(tmp_path))
    script = 'set -euo pipefail\nunset gg_tmp_dir\n'
    script += '\n'.join(key + '=' + shlex.quote(value) for key, value in values.items()) + '\n'
    script += 'gg_artifact_prepare_stage() { printf -v "$1" "%s" 1; }\n'
    script += function + '\nrun_gene_annotation_database_for_source\n'
    original_sha = annotation.sha256(cohort['db'])
    result = subprocess.run(['bash', '-c', script], env={**os.environ, 'TMPDIR': str(scratch)},
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    target = family / 'gg_orthogroup.annotated.db'
    assert target.exists() and annotation.sha256(cohort['db']) == original_sha
    with sqlite3.connect(target) as db:
        meta = dict(db.execute('SELECT key,value FROM gene_annotation_meta'))
        assert json.loads(meta['audit'])['reconciliation_genes_without_annotation'] == 0
    assert (summary / 'artifact_provenance/orthogroup.gene_annotation_database.json').is_file()


def test_native_provenance_requires_exact_sources_and_records_derivative(cohort, tmp_path):
    record_path = tmp_path / 'provenance.json'
    argv = ['--manifest', str(record_path), '--workspace-root', str(tmp_path),
            '--logical-root', str(tmp_path / 'global'), '--step', 'gene_annotation_database',
            '--family-id', 'orthogroup', '--input', 'database=' + str(cohort['db']),
            '--input', 'manifest=' + str(cohort['manifest']), '--output', 'database=' + str(cohort['target'])]
    with pytest.raises(ValueError, match='exactly'):
        build(cohort, record_argv=argv)
    assert not cohort['target'].exists() and not record_path.exists()
    for index, source in enumerate(cohort['config']['sources']):
        argv += ['--input', 'annotation_' + str(index) + '=' + str(tmp_path / source['path'])]
    build(cohort, record_argv=argv)
    record = json.loads(record_path.read_text())
    assert len(record['inputs']) == 5 and len(record['outputs']) == 1
    assert record['outputs'][0]['sha256'] == annotation.sha256(cohort['target'])


def test_unbound_legacy_database_does_not_claim_model_compatibility(cohort):
    with sqlite3.connect(cohort['db']) as db:
        db.execute('DROP TABLE reconciliation_meta')
    cohort['config']['source_database_sha256'] = annotation.sha256(cohort['db'])
    cohort['manifest'].write_text(json.dumps(cohort['config']))
    with pytest.raises(ValueError, match='model-bound reconciliation DB'):
        build(cohort)
    assert not cohort['target'].exists()
