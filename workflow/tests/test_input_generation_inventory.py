import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
from input_generation_inventory import inspect


def plan(tmp, name, species, done, lineage='embryophyta_odb12', finalized=False):
    path = tmp / name / 'output/input_generation/tmp/task_plan.json'
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps({'task_count': len(species), 'tasks': [{'species_prefix': s} for s in species]}))
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    Path(str(path) + '.settings.json').write_text(json.dumps({'busco_lineage': lineage}))
    receipts = Path(str(path) + '.completed')
    receipts.mkdir()
    for i, sp in enumerate(species, 1):
        if sp in done:
            (receipts / f'{i}.json').write_text(json.dumps({'plan_sha256': sha, 'task_index': i,
                                                         'species_prefix': sp, 'files': {'raw': 'a'*64}}))
    if finalized:
        (path.parent.parent / 'gg_input_generation_runs.tsv').write_text(
            'input_generation_mode\texit_code\tstage_multispecies_summary_status\narray_finalize\t0\tok\n')
    return {'path': str(path), 'sha256': sha}


def test_updates_are_deduplicated_and_subset_finalize_is_not_whole_completion(tmp_path):
    original = plan(tmp_path, 'original', ['A', 'B', 'C'], {'A'})
    update = plan(tmp_path, 'update', ['A', 'B', 'D'], {'A', 'B', 'D'}, finalized=True)
    old = plan(tmp_path, 'old', ['C'], {'C'}, lineage='eukaryota_odb12')
    request = {'cohort': [original, update], 'sources': [original, update, old], 'lineage': 'embryophyta_odb12'}
    result = inspect(request)
    assert (result['expected'], result['recorded'], result['remaining']) == (4, 3, 1)
    assert not result['finalize_recorded']
    assert result['workflow_verified'] is False
    assert result['checksum_verified'] is False
    receipt = Path(original['path'] + '.completed') / '2.json'
    receipt.write_text(json.dumps({'plan_sha256': 'b'*64, 'task_index': 2, 'species_prefix': 'B', 'files': {'raw': 'a'*64}}))
    # A corrupted superseded receipt does not erase a valid updated publication.
    assert inspect(request)['invalid'] == 0
    Path(original['path']).write_text('{}')
    with pytest.raises(ValueError, match='plan changed'):
        inspect(request)


def test_all_species_records_still_do_not_certify_checksums(tmp_path):
    source = plan(tmp_path, 'all', ['A', 'B'], {'A', 'B'}, finalized=True)
    result = inspect({'cohort': [source], 'sources': [source], 'lineage': 'embryophyta_odb12'})
    assert result['remaining'] == 0 and result['finalize_recorded']
    assert not result['workflow_verified'] and not result['checksum_verified']


def test_public_api_inventory_is_readonly_and_bounded(tmp_path, monkeypatch):
    source = plan(tmp_path, 'all', ['A'], {'A'})
    request = {'cohort': [source], 'sources': [source], 'lineage': 'embryophyta_odb12'}
    before = {str(p): (p.stat().st_mtime_ns, p.read_bytes()) for p in tmp_path.rglob('*') if p.is_file()}
    api = Path(__file__).resolve().parents[1] / 'support/workflow_api.py'
    result = subprocess.run([sys.executable, str(api), 'input-inventory', '--request-json', json.dumps(request)],
                            capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    payload = json.loads(result.stdout)
    assert payload['read_only'] and not payload['execution_authorized']
    assert payload['inventory']['recorded'] == 1
    assert before == {str(p): (p.stat().st_mtime_ns, p.read_bytes()) for p in tmp_path.rglob('*') if p.is_file()}
    monkeypatch.setattr('input_generation_inventory.MAX_BYTES', 1)
    with pytest.raises(ValueError, match='exceeds'):
        inspect(request)
    raw = Path(source['path'])
    target = raw.with_suffix('.original')
    raw.rename(target)
    raw.symlink_to(target)
    monkeypatch.setattr('input_generation_inventory.MAX_BYTES', 8 * 1024 * 1024)
    with pytest.raises(OSError):
        inspect(request)
