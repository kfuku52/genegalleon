"""Fresh-byte equivalence, bounded pair fences, and immutable-generation refusal."""
import copy
import hashlib
import json
import os
import shutil
import sys
from collections import Counter
from importlib import import_module
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / 'support'
sys.path.insert(0, str(SUPPORT))
rescue = import_module('rescue_gene_models')
snapshot_module = import_module('rescue_prepared_snapshot')
PreparedSnapshot = snapshot_module.PreparedSnapshot
state = import_module(rescue.digest.__module__)


def fixture(tmp_path, count=3):
    root = tmp_path / 'run'
    root.mkdir()
    sources, hashes = {}, {}
    names = [f'Species_{number}' for number in range(count)]
    for number, name in enumerate(names):
        sources[name] = {}
        for key in ('fasta', 'gff', 'genome'):
            path = tmp_path / (name + '.' + key)
            path.write_text(f'{key}:{number}\n' + 'A' * (8192 if key == 'genome' else 16))
            sources[name][key] = str(path)
            hashes[str(path)] = rescue.digest(path)
    jobs = [dict(a=names[0], b=name, kind='pair', index=index, id=f'comparison_{index:06d}')
            for index, name in enumerate(names[1:], 1)]
    owners = ('kfFractBias', 'jcvi', 'biopython', 'numpy', 'natsort', 'more-itertools', 'python',
              'diamond', 'diamond_sha256', 'lastal', 'lastal_sha256', 'lastdb', 'lastdb_sha256')
    plan = {'request': {'sources': sources, 'files': hashes,
                        'parameters': dict(cscore=.7, min_anchors=3, distance=20, diagonal_bound=4),
                        'tools': {**dict.fromkeys(owners, 'fixture'), 'source_hashes': {'jcvi.fake': 'fixture'}}},
            'synteny_jobs': jobs}
    rescue.atomic_json(root / 'plan.json', plan)
    plan_hash = rescue.digest(root / 'plan.json')
    for number, name in enumerate(names):
        directory = root / 'prepared' / name
        directory.mkdir(parents=True)
        for filename, text in [('genes.bed', f'chr1\t0\t12\tg{number}\n'),
                               ('genes.pep', f'>g{number}\nMKP\n'),
                               ('genes.id_map.tsv', 'original_id\tjcvi_id\tlocus_id\tstatus\ng\tg\tg\tselected\n')]:
            (directory / filename).write_text(text)
        rescue.atomic_json(directory / 'receipt.json', {'key': {'plan': plan_hash, 'species': name},
                           'files': rescue.hash_outputs(directory)})
    for job in jobs:
        directory = root / 'synteny' / job['id']
        directory.mkdir(parents=True)
        rescue.atomic_json(directory / 'blocks.json', [])
        rescue.atomic_json(directory / 'receipt.json', {'key': rescue.comparison_key(root, job),
                                                       'files': rescue.hash_outputs(directory)})
    return root, plan, names


def metrics(monkeypatch):
    values = Counter()
    monkeypatch.setattr(state, 'count', lambda name, amount=1: values.update({name: amount}))
    return values


def test_reused_pairs_keep_exact_keys_and_hash_each_frozen_input_once(tmp_path, monkeypatch):
    root, plan, names = fixture(tmp_path)
    keys = [(rescue.comparison_key(root, job), rescue.comparison_cache_key(root, plan, job))
            for job in plan['synteny_jobs']]
    measured = metrics(monkeypatch)
    for job in plan['synteny_jobs']:
        rescue.synteny(root, plan, job['index'], 1)
    baseline = dict(measured)
    measured.clear()
    expected_paths = [root / 'plan.json'] + [Path(path) for path in plan['request']['files']]
    expected_paths += [path for name in names for path in (root / 'prepared' / name).iterdir()]
    expected_paths += [root / 'synteny' / job['id'] / 'blocks.json' for job in plan['synteny_jobs']]
    with PreparedSnapshot(root, plan, rescue) as snapshot:
        for job, expected in zip(plan['synteny_jobs'], keys, strict=True):
            assert rescue.comparison_key(root, job, snapshot) == expected[0]
            assert rescue.comparison_cache_key(root, plan, job, snapshot) == expected[1]
            rescue.synteny(root, plan, job['index'], 1, prepared_snapshot=snapshot)
    assert measured['sha256_reads'] == len(expected_paths)
    assert measured['sha256_bytes'] == sum(path.stat().st_size for path in expected_paths)
    assert measured['sha256_bytes'] < baseline['sha256_bytes']
    first = dict(measured)
    measured.clear()
    # A new command must re-read content, even without any stat changes.
    with PreparedSnapshot(root, plan, rescue) as snapshot:
        for job in plan['synteny_jobs']:
            rescue.synteny(root, plan, job['index'], 1, prepared_snapshot=snapshot)
    assert dict(measured) == first


@pytest.mark.parametrize('target', ['source', 'member', 'receipt', 'plan'])
@pytest.mark.parametrize('mutation', ['rewrite', 'replace', 'remove', 'permission', 'symlink'])
def test_observed_generation_changes_fail_without_refresh(tmp_path, monkeypatch, target, mutation):
    root, plan, names = fixture(tmp_path)
    paths = {'source': Path(plan['request']['sources'][names[0]]['genome']),
             'member': root / 'prepared' / names[0] / 'genes.pep',
             'receipt': root / 'prepared' / names[0] / 'receipt.json', 'plan': root / 'plan.json'}
    snapshot = PreparedSnapshot(root, plan, rescue)
    snapshot.prepared(names[0])
    path = paths[target]
    old = path.stat()
    data = path.read_bytes()
    if mutation == 'rewrite':
        path.write_bytes(bytes([data[0] ^ 1]) + data[1:])
        os.utime(path, ns=(old.st_atime_ns, old.st_mtime_ns))
    elif mutation == 'replace':
        replacement = path.with_name(path.name + '.replacement')
        replacement.write_bytes(data)
        os.utime(replacement, ns=(old.st_atime_ns, old.st_mtime_ns))
        replacement.replace(path)
    elif mutation == 'remove':
        path.unlink()
    elif mutation == 'permission':
        path.chmod(0o600 if old.st_mode & 0o077 else 0o644)
    else:
        replacement = path.with_name(path.name + '.target')
        replacement.write_bytes(data)
        path.unlink()
        path.symlink_to(replacement)
    monkeypatch.setattr(rescue, 'prepared', lambda *_args, **_kwargs: pytest.fail('Changed proof must not rebuild'))
    with pytest.raises((OSError, ValueError)):
        snapshot.prepared(names[0])
    snapshot.close()
    with pytest.raises(ValueError, match='closed'):
        snapshot.plan_digest()


def test_symlink_replaced_with_same_target_is_not_a_new_proof(tmp_path):
    root, plan, names = fixture(tmp_path)
    path = Path(plan['request']['sources'][names[0]]['genome'])
    target = path.with_suffix('.real')
    path.rename(target)
    path.symlink_to(target)
    snapshot = PreparedSnapshot(root, plan, rescue)
    snapshot.prepared(names[0])
    path.unlink()
    path.symlink_to(target)
    with pytest.raises(OSError):
        snapshot.check([names[0]])
    snapshot.close()


@pytest.mark.parametrize('damage', ['missing', 'corrupt_member', 'malformed_receipt', 'wrong_key', 'journal'])
def test_first_unverified_publication_keeps_normal_recovery(tmp_path, monkeypatch, damage):
    root, plan, names = fixture(tmp_path)
    name = names[0]
    directory = root / 'prepared' / name
    if damage == 'missing':
        shutil.rmtree(directory)
    elif damage == 'corrupt_member':
        (directory / 'genes.pep').write_text('broken')
    elif damage == 'malformed_receipt':
        (directory / 'receipt.json').write_text('{}')
    elif damage == 'wrong_key':
        value = json.loads((directory / 'receipt.json').read_text())
        value['key']['plan'] = 'stale'
        rescue.atomic_json(directory / 'receipt.json', value)
    else:
        token = hashlib.sha256(('prepared/' + name).encode()).hexdigest()[:12]
        backup = directory.parent / ('.previous-' + token + '-fixture')
        directory.rename(backup)
        lock = root / '.locks' / ('prepared__' + name + '.lock')
        lock.parent.mkdir()
        rescue.atomic_json(lock.with_suffix('.publish.json'),
                           {'destination': name, 'backup': backup.name,
                            'temporary': '.working-' + token + '-fixture',
                            'key': {'plan': rescue.digest(root / 'plan.json'), 'species': name}})
    calls = []

    def rebuild(_source, temporary, _prefix, _score):
        calls.append('build')
        (temporary / 'genes.pep').write_text('>fixed\nMK\n')
        (temporary / 'genes.bed').write_text('chr1\t0\t9\tfixed\n')
        return [], {}

    monkeypatch.setattr(rescue, 'prepare_rescue_genome', rebuild)
    snapshot = PreparedSnapshot(root, plan, rescue)
    assert snapshot.prepared(name) == directory
    assert rescue.verified(directory, {'plan': snapshot.plan_hash, 'species': name})
    assert bool(calls) == (damage != 'journal')
    (directory / 'genes.pep').write_text('changed after fresh proof')
    with pytest.raises(OSError):
        snapshot.prepared(name)
    snapshot.close()


def test_source_corruption_cannot_be_repaired_as_a_prepared_output(tmp_path, monkeypatch):
    root, plan, names = fixture(tmp_path)
    Path(plan['request']['sources'][names[0]]['genome']).write_text('changed')
    monkeypatch.setattr(rescue, 'prepared', lambda *_args, **_kwargs: pytest.fail('No producer for changed input'))
    snapshot = PreparedSnapshot(root, plan, rescue)
    with pytest.raises(ValueError, match='Frozen rescue input changed'):
        snapshot.prepared(names[0])
    snapshot.close()


def test_stale_or_mutated_in_memory_plan_cannot_stamp_new_keys(tmp_path):
    root, plan, names = fixture(tmp_path)
    stale = copy.deepcopy(plan)
    stale['request']['sources'][names[0]]['genetic_code'] = 4
    with pytest.raises(ValueError, match='Frozen rescue plan changed'):
        PreparedSnapshot(root, stale, rescue)
    snapshot = PreparedSnapshot(root, plan, rescue)
    plan['synteny_jobs'][0]['b'] = names[-1]
    with pytest.raises(ValueError, match='Frozen rescue comparison changed'):
        rescue.comparison_key(root, plan['synteny_jobs'][0], snapshot)
    snapshot.close()


def test_pair_fences_do_not_visit_all_previously_used_species(tmp_path, monkeypatch):
    root, plan, names = fixture(tmp_path, count=50)
    snapshot = PreparedSnapshot(root, plan, rescue)
    for name in names:
        snapshot.prepared(name)
    for name in names[2:]:
        monkeypatch.setattr(snapshot.species[name]['proof'], 'check',
                            lambda: pytest.fail('Pair check scanned an unrelated species'))
        monkeypatch.setattr(snapshot.sources[name], 'check',
                            lambda: pytest.fail('Pair check scanned an unrelated source'))
    rescue.comparison_key(root, plan['synteny_jobs'][0], snapshot)
    rescue.comparison_cache_key(root, plan, plan['synteny_jobs'][0], snapshot)
    snapshot.close()


def test_exception_discards_command_proof_and_does_not_hide_failure(tmp_path):
    root, plan, names = fixture(tmp_path)
    snapshot = PreparedSnapshot(root, plan, rescue)
    with pytest.raises(RuntimeError, match='builder failed'):
        with snapshot:
            snapshot.prepared(names[0])
            raise RuntimeError('builder failed')
    assert snapshot.closed and snapshot.species == {} and snapshot.sources == {}
    with PreparedSnapshot(root, plan, rescue) as next_command:
        assert next_command.prepared(names[0]).is_dir()


def test_source_only_proof_is_fenced_at_exit(tmp_path):
    root, plan, names = fixture(tmp_path)
    snapshot = PreparedSnapshot(root, plan, rescue)
    with pytest.raises(OSError):
        with snapshot:
            snapshot.verify_source(names[0])
            assert names[0] not in snapshot.species
            Path(plan['request']['sources'][names[0]]['genome']).write_text('changed')
    assert snapshot.closed


@pytest.mark.parametrize('mutation', ['add_member', 'remove_member', 'replace_parent'])
def test_prepared_directory_and_members_keep_the_same_generation(tmp_path, mutation):
    root, plan, names = fixture(tmp_path)
    snapshot = PreparedSnapshot(root, plan, rescue)
    directory = snapshot.prepared(names[0])
    if mutation == 'add_member':
        (directory / 'unreceipted.txt').write_text('new')
    elif mutation == 'remove_member':
        (directory / 'genes.id_map.tsv').unlink()
    else:
        # Replacing the parent with an equivalent symlink cannot authorize the
        # old hashes, even though all member inodes/content remain unchanged.
        parent = directory.parent
        renamed = parent.with_name('prepared_original')
        parent.rename(renamed)
        parent.symlink_to(renamed, target_is_directory=True)
    with pytest.raises(OSError):
        snapshot.check([names[0]])
    snapshot.close()


def test_mutation_during_first_hash_is_not_hidden_by_initial_repair(tmp_path, monkeypatch):
    root, plan, names = fixture(tmp_path)
    path = root / 'prepared' / names[0] / 'genes.pep'
    original = state.digest

    def changing_digest(target):
        result = original(target)
        if Path(target) == path:
            path.write_text('changed after content read')
        return result

    monkeypatch.setattr(state, 'digest', changing_digest)
    monkeypatch.setattr(rescue, 'prepared', lambda *_args, **_kwargs: pytest.fail('Do not hide a hashing race'))
    snapshot = PreparedSnapshot(root, plan, rescue)
    with pytest.raises(snapshot_module.PreparationChanged):
        snapshot.prepared(names[0])
    snapshot.close()


def test_unreceipted_required_file_has_no_unhashed_fallback(tmp_path):
    root, plan, names = fixture(tmp_path)
    path = root / 'prepared' / names[0] / 'receipt.json'
    receipt = json.loads(path.read_text())
    del receipt['files']['genes.pep']
    rescue.atomic_json(path, receipt)
    snapshot = PreparedSnapshot(root, plan, rescue)
    with pytest.raises(ValueError, match='absent from receipt'):
        rescue.comparison_cache_key(root, plan, plan['synteny_jobs'][0], snapshot)
    snapshot.close()


@pytest.mark.parametrize('publication', ['comparison_cache', 'synteny'])
def test_source_mutated_during_output_hashing_refuses_intermediate_publication(tmp_path, monkeypatch, publication):
    root, plan, names = fixture(tmp_path)
    job = plan['synteny_jobs'][0]
    destination = root / 'synteny' / job['id']
    shutil.rmtree(destination)
    cache = tmp_path / 'cache'
    original = rescue.hash_outputs

    def changing_hashes(directory, **kwargs):
        result = original(directory, **kwargs)
        parent = cache / 'comparisons' if publication == 'comparison_cache' else root / 'synteny'
        if directory.parent == parent:
            (root / 'prepared' / names[0] / 'genes.pep').write_text('changed while hashing stage outputs')
        return result

    monkeypatch.setattr(rescue, 'hash_outputs', changing_hashes)
    monkeypatch.setattr(rescue, 'build_comparison',
                        lambda directory, *_args: rescue.atomic_json(directory / 'blocks.json', []))
    # Compute algorithm-bound keys after installing the fake builder too.
    cache_key = rescue.comparison_cache_key(root, plan, job)
    cache_id = hashlib.sha256(json.dumps(cache_key, sort_keys=True).encode()).hexdigest()
    target = cache / 'comparisons' / cache_id if publication == 'comparison_cache' else destination
    snapshot = PreparedSnapshot(root, plan, rescue)
    with pytest.raises(OSError):
        rescue.synteny(root, plan, job['index'], 1, cache, prepared_snapshot=snapshot)
    assert not (target / 'receipt.json').exists()
    assert target.with_name(target.name + '.failed').is_dir()
    assert not (destination / 'receipt.json').exists()
    snapshot.close()


def test_rebuilt_preparation_refuses_source_mutation_during_output_hashing(tmp_path, monkeypatch):
    root, plan, names = fixture(tmp_path)
    name = names[0]
    directory = root / 'prepared' / name
    shutil.rmtree(directory)
    original = rescue.hash_outputs

    def rebuild(_source, temporary, _prefix, _score):
        (temporary / 'genes.pep').write_text('>fixed\nMK\n')
        return [], {}

    def changing_hashes(path, **kwargs):
        result = original(path, **kwargs)
        Path(plan['request']['sources'][name]['genome']).write_text('changed during preparation output hash')
        return result

    monkeypatch.setattr(rescue, 'prepare_rescue_genome', rebuild)
    monkeypatch.setattr(rescue, 'hash_outputs', changing_hashes)
    snapshot = PreparedSnapshot(root, plan, rescue)
    with pytest.raises(OSError):
        snapshot.prepared(name)
    assert not (directory / 'receipt.json').exists()
    assert directory.with_name(name + '.failed').is_dir()
    snapshot.close()


def test_receipt_becoming_malformed_after_digest_is_not_initial_corruption(tmp_path, monkeypatch):
    root, plan, names = fixture(tmp_path)
    receipt = root / 'prepared' / names[0] / 'receipt.json'
    original = Path.read_bytes

    def changing_read(path):
        if path == receipt:
            path.write_bytes(b'{invalid JSON')
        return original(path)

    monkeypatch.setattr(Path, 'read_bytes', changing_read)
    monkeypatch.setattr(rescue, 'prepared', lambda *_args, **_kwargs: pytest.fail('Do not repair an observed receipt race'))
    snapshot = PreparedSnapshot(root, plan, rescue)
    with pytest.raises(snapshot_module.PreparationChanged):
        snapshot.prepared(names[0])
    snapshot.close()
