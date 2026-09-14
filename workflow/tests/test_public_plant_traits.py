import io
import json
import sys
from pathlib import Path
from urllib.error import HTTPError

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
import public_plant_traits as p


class Fixture:
    def __init__(self, tables=None, payloads=None):
        self.tables = tables or {}
        self.payloads = iter(payloads or [])
        self.decisions = []
    def csv(self, url, **kwargs):
        return self.tables[int(url.rsplit('/', 1)[-1])].copy()
    def get(self, url, **kwargs):
        return json.dumps(next(self.payloads)).encode()
    def json(self, url, **kwargs):
        return next(self.payloads)


def test_brot_does_not_pool_representations_or_promote_subspecies():
    table = pd.DataFrame([
        ['1', 'Quercus robur', 'RespFire', 'yes', '[2]', 'boolean', 'paper'],
        ['2', 'Quercus robur', 'RespFire', '80', '%', 'quantitative', 'paper'],
        ['3', 'Quercus robur subsp. robur', 'RespFire', 'yes', '[2]', 'boolean', 'paper'],
    ], columns=['ID', 'Taxon', 'Trait', 'Data', 'Units', 'DataType', 'SourceID'])
    rows = p.brot(Fixture({11194784: table}), {}, {'Quercus robur': 'Quercus_robur'})
    assert len(rows) == 2
    assert {row[1] for row in rows} == {'RespFire|boolean|[2]', 'RespFire|quantitative|%'}
    assert json.loads(rows[0][-1])['SourceID'] == 'paper'


def cpt_tables():
    tables = {file_id: pd.DataFrame({'SAMPLE ID': ['1'], 'x': ['2']}) for file_id in p.CPT_FILES.values()}
    tables[p.CPT_FILES['species']] = pd.DataFrame({'SAMPLE ID': ['1'], 'SPECIES ID': ['SP1'], 'ACCEPTED GENUS': ['Allium'], 'ACCEPTED SPECIES': ['senescens']})
    tables[p.CPT_FILES['taxonomy']] = tables[p.CPT_FILES['species']].drop(columns='SAMPLE ID')
    tables[p.CPT_FILES['pathway']] = pd.DataFrame({'SPECIES ID': ['SP1'], 'Photo_Path': ['C3']})
    return tables


def test_cpt_joins_samples_and_species_and_excludes_flagged():
    tables = cpt_tables()
    tables[p.CPT_FILES['chemical']]['flagged'] = 'unreliable'
    rows = p.cpt(Fixture(tables), {}, {'Allium senescens': 'Allium_senescens'})
    assert 'pathway:Photo_Path' in {r[1] for r in rows}
    assert not any(r[1].startswith('chemical:') for r in rows)
    assert all(r[0] == 'Allium_senescens' for r in rows)


@pytest.mark.parametrize('problem', ['duplicate', 'unmapped'])
def test_cpt_rejects_bad_joins(problem):
    tables = cpt_tables()
    if problem == 'duplicate':
        tables[p.CPT_FILES['species']] = pd.concat([tables[p.CPT_FILES['species']]] * 2)
    else:
        tables[p.CPT_FILES['chemical']]['SAMPLE ID'] = 'missing'
    with pytest.raises(ValueError):
        p.cpt(Fixture(tables), {}, {'Allium senescens': 'Allium_senescens'})


def test_bien_validates_species_and_preserves_unit_and_citation():
    record = dict(scrubbed_species_binomial='Quercus robur', trait_name='height', trait_value='3', unit='m', id=1, source_citation='paper', access='public')
    rows = p.bien(Fixture(payloads=[{'data': [record]}]), {'response_format': 'json'}, {'Quercus robur': 'Quercus_robur'})
    assert rows[0][1:4] == ['height|m', '3', 'm']
    assert rows[0][5] == 'paper'
    with pytest.raises(ValueError, match='unexpected species'):
        p.bien(Fixture(payloads=[{'data': [record]}]), {'response_format': 'json'}, {'Pinus pinaster': 'Pinus_pinaster'})


def test_bien_service_failure_is_not_missing_data(tmp_path):
    def fail(request, **kwargs):
        raise HTTPError(request.full_url, 404, 'not found', {}, io.BytesIO(b'<html>Cannot GET route</html>'))
    with pytest.raises(ValueError, match='endpoint unavailable'):
        p.load_public_traits('bien', {}, ['Quercus_robur'], tmp_path, 1, False, opener=fail)
    report = json.loads(next((tmp_path / 'bien/runs').glob('*.json')).read_text())
    assert report['status'] == 'failed'
    assert not list((tmp_path / 'bien/runs').glob('*.tsv'))


def test_bien_documented_no_data_response(tmp_path):
    def empty(request, **kwargs):
        data = b'{"message":"No data found for species Quercus robur"}'
        raise HTTPError(request.full_url, 404, 'not found', {}, io.BytesIO(data))
    table = p.load_public_traits('bien', {}, ['Quercus_robur'], tmp_path, 1, False, opener=empty)
    assert table.empty
    report = json.loads(next((tmp_path / 'bien/runs').glob('*.json')).read_text())
    assert report['status'] == 'complete' and len(report['requests']) == 1


def test_algae_does_not_inherit_and_preserves_measurement_context():
    hit = dict(scientificname='Ulva lactuca', rank='Species', status='accepted', AphiaID=1)
    attribute = dict(AphiaID=1, measurementTypeID=15, measurementType='Body size', measurementValue='20', AphiaID_Inherited=1, source_id=2, children=[{'measurementType': 'Unit', 'measurementValue': 'cm'}])
    inherited = dict(attribute, AphiaID_Inherited=9)
    rows = p.algaetraits(Fixture(payloads=[[hit], [attribute, inherited]]), {}, {'Ulva lactuca': 'Ulva_lactuca'})
    assert len(rows) == 1
    assert json.loads(rows[0][2])['qualifiers'][0]['measurementValue'] == 'cm'


def test_algae_ambiguous_name_is_not_selected():
    hit = dict(scientificname='Ulva lactuca', rank='Species', status='accepted', AphiaID=1)
    assert p.algaetraits(Fixture(payloads=[[hit, dict(hit, AphiaID=2)]]), {}, {'Ulva lactuca': 'Ulva_lactuca'}) == []


def test_dry_run_has_no_download_or_write(tmp_path):
    def forbidden(*args, **kwargs):
        raise AssertionError('network access')
    assert p.load_public_traits('brot', {}, ['Quercus_robur'], tmp_path, 1, True, opener=forbidden) is None
    assert not list(tmp_path.iterdir())


def test_generator_dispatches_and_publishes_schema(tmp_path):
    from test_generate_species_trait import run_script
    source = tmp_path / 'brot.csv'
    source.write_text('ID,Taxon,Trait,Data,Units,DataType,SourceID\n1,Quercus robur,Height,4,m,quantitative,paper\n')
    # Exercise anonymous HTTP acquisition through the complete generator CLI.
    import threading
    from functools import partial
    from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
    server = ThreadingHTTPServer(('127.0.0.1', 0), partial(SimpleHTTPRequestHandler, directory=str(tmp_path)))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        manifest = tmp_path / 'manifest.tsv'
        manifest.write_text('provider\tid\tspecies_key\nlocal\tx\tQuercus_robur\nlocal\ty\tPinus_pinaster\n')
        sources = tmp_path / 'sources.tsv'
        sources.write_text('database\tacquisition_mode\turi\tspecies_column\ttrait_key_column\nbrot\tpublic_plant_traits\thttp://127.0.0.1:%s/brot.csv\tspecies\ttrait_key\n' % server.server_port)
        plan = tmp_path / 'plan.tsv'
        plan.write_text('database\tsource_column\toutput_trait\tvalue_type\taggregation\ttrait_key\ttrait_key_column\nbrot\tvalue\theight_m\tnumeric\tmedian\tHeight|quantitative|m\ttrait_key\n')
        output = tmp_path / 'traits.tsv'
        result = run_script('--download-manifest', str(manifest), '--trait-plan', str(plan), '--database-sources', str(sources), '--downloads-dir', str(tmp_path / 'downloads'), '--output', str(output), '--strict')
        assert result.returncode == 0, result.stdout + result.stderr
        table = pd.read_csv(output, sep='\t').set_index('species')
        assert table.loc['Quercus_robur', 'height_m'] == 4
        assert pd.isna(table.loc['Pinus_pinaster', 'height_m'])
        assert Path(str(output) + '.schema.json').exists()
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


def test_cpt_pathway_uses_taxonomy_not_sample_subset():
    tables = cpt_tables()
    tables[p.CPT_FILES['taxonomy']] = pd.DataFrame({'SPECIES ID': ['SP2'], 'ACCEPTED GENUS': ['Pinus'], 'ACCEPTED SPECIES': ['pinaster']})
    tables[p.CPT_FILES['pathway']] = pd.DataFrame({'SPECIES ID': ['SP2'], 'Photo_Path': ['C3']})
    rows = p.cpt(Fixture(tables), {}, {'Pinus pinaster': 'Pinus_pinaster'})
    assert len(rows) == 1 and rows[0][2] == 'C3'


def test_hydraulic_csv_encoding(tmp_path):
    data = 'SAMPLE ID,Ψtlp\n1,-1.2\n'.encode('gb18030')
    acquisition = p.Acquisition(tmp_path, 1, lambda *args, **kwargs: io.BytesIO(data))
    assert list(acquisition.csv('https://example.org/hydraulic', encoding='gb18030')) == ['SAMPLE ID', 'Ψtlp']


def test_missing_source_tokens_not_text_labels(tmp_path):
    data = b'ID,Taxon,Trait,Data,Units,DataType,SourceID\n1,Quercus robur,GrowthForm,NA,[13],categorical,paper\n'
    table = p.load_public_traits('brot', {}, ['Quercus_robur'], tmp_path, 1, False, opener=lambda *args, **kwargs: io.BytesIO(data))
    assert table.empty
    report = json.loads(next((tmp_path / 'brot/runs').glob('*.json')).read_text())
    assert report['missing_value_rows'] == 1


def test_bien_csv_live_contract_and_nonpublic_exclusion(tmp_path):
    data = b'trait_id,scrubbed_species_binomial,trait_name,trait_value,unit,access,url_source\n1,Quercus robur,height,4,m,public,paper\n2,Quercus robur,height,7,m,private,paper2\n'
    table = p.load_public_traits('bien', {}, ['Quercus_robur'], tmp_path, 1, False, opener=lambda *args, **kwargs: io.BytesIO(data))
    assert table.value.tolist() == ['4']
    assert table.source_reference.tolist() == ['paper']
    report = json.loads(next((tmp_path / 'bien/runs').glob('*.json')).read_text())
    assert report['species_decisions'][0]['status'] == 'excluded_nonpublic'
