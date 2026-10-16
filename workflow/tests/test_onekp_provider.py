import sys
from pathlib import Path
from unittest import mock

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
import format_species_inputs as fsi
from format_species_providers import onekp


@pytest.fixture(autouse=True)
def empty_index_cache():
    onekp.sample_index.cache_clear()
    yield
    onekp.sample_index.cache_clear()


def records(**updates):
    row = dict(code='BMJR', species='Adiantum raddianum', is_mislabeled=False,
               assembly_url='https://ftp.cngb.org/pub/SciRAID/onekp/assemblies/BMJR-Adiantum_raddianum/')
    row.update(updates)
    return {'data': {'samples': [row], 'total': 1}}


def test_code_resolves_published_cds_and_preserves_curated_species():
    with mock.patch.object(onekp, 'fetch_json_with_headers', return_value=records()) as fetch:
        result = onekp.resolve_onekp_download_urls_from_id('bmjr', 'Curated_species', 10, {})
        assert result['cds_url'].endswith('/BMJR-translated-nucleotides.fa.gz')
        assert result['species_key'] == 'Curated_species'
        assert result['onekp_reported_species'] == 'Adiantum raddianum'
        onekp.resolve_onekp_download_urls_from_id('BMJR', '', 10, {})
        fetch.assert_called_once()


@pytest.mark.parametrize('update', [{'is_mislabeled': True}, {'is_mislabeled': None}, {'assembly_url': 'https://other.example/file'}, {'code': 'OTHER'}])
def test_reject_unverified_sample_identity_and_location(update):
    with mock.patch.object(onekp, 'fetch_json_with_headers', return_value=records(**update)):
        with pytest.raises(ValueError):
            onekp.resolve_onekp_download_urls_from_id('BMJR', '', 10, {})


def test_incomplete_sample_index_is_not_cached():
    data = records()
    data['data']['total'] = 2
    with mock.patch.object(onekp, 'fetch_json_with_headers', return_value=data):
        with pytest.raises(ValueError, match='incomplete'):
            onekp.resolve_onekp_download_urls_from_id('BMJR', '', 10, {})


def test_provider_discovers_only_coding_nucleotides(tmp_path):
    species = tmp_path / 'Adiantum_raddianum'
    species.mkdir()
    (species/'BMJR.cds.fa.gz').write_bytes(b'fixture')
    (species/'BMJR-translated-protein.fa.gz').write_bytes(b'protein')
    tasks, _, errors = fsi.discover_tasks('onekp', tmp_path)
    assert not errors
    assert tasks[0]['cds_path'].name == 'BMJR.cds.fa.gz'
