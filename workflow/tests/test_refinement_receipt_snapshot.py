"""Fresh receipt proofs, bounded reuse, publication races and lazy resume reads."""
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
refinement = import_module('gene_model_refinement')
rescue = refinement.rescue
snapshots = import_module('refinement_receipt_snapshot')
state = import_module(rescue.digest.__module__)


def publication(tmp_path, name='Species_target', size=1024 * 1024):
    directory = tmp_path / 'rescued' / name
    directory.mkdir(parents=True)
    (directory / 'intervals').mkdir()
    (directory / 'intervals/regions.fa').write_text('>region\n' + 'ACGT' * (size // 4) + '\n')
    rescue.atomic_json(directory / 'models.json', [{'model_id': name, 'sequence': 'ATGAAATAA'}])
    key = {'plan': 'frozen-fixture-plan', 'species': name}
    rescue.atomic_json(directory / 'receipt.json', {'key': key, 'files': rescue.hash_outputs(directory)})
    return directory, key, {str(directory / 'receipt.json'): rescue.digest(directory / 'receipt.json')}


def metrics(monkeypatch):
    values = Counter()
    monkeypatch.setattr(state, 'count', lambda name, amount=1: values.update({name: amount}))
    return values


def test_same_receipt_and_output_bytes_with_one_fresh_read_per_command(tmp_path, monkeypatch):
    directory, key, receipts = publication(tmp_path)
    before = {str(p.relative_to(directory)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in directory.rglob('*') if p.is_file()}
    measured = metrics(monkeypatch)
    for _ in range(4):
        assert rescue.digest(directory / 'receipt.json') == next(iter(receipts.values()))
        assert rescue.verified(directory, key)
    baseline = dict(measured)
    measured.clear()
    with snapshots.ReceiptSnapshot(rescue.verified) as snapshot:
        for _ in range(4):
            snapshot.verify(receipts)
            snapshot.check(receipts)
    optimized = dict(measured)
    assert optimized['sha256_reads'] == len(before)
    assert optimized['sha256_bytes'] == sum(p.stat().st_size for p in directory.rglob('*') if p.is_file())
    assert baseline == {name: value * 4 for name, value in optimized.items()}
    measured.clear()
    with snapshots.ReceiptSnapshot(rescue.verified) as fresh_command:
        fresh_command.verify(receipts)
    assert dict(measured) == optimized
    assert before == {str(p.relative_to(directory)): hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in directory.rglob('*') if p.is_file()}


@pytest.mark.parametrize('target', ['member', 'receipt'])
@pytest.mark.parametrize('mutation', ['rewrite', 'replace', 'remove', 'permission', 'symlink'])
def test_verified_generation_cannot_refresh(tmp_path, target, mutation):
    directory, _key, receipts = publication(tmp_path)
    calls = []

    def verifier(*args):
        calls.append(args)
        return rescue.verified(*args)

    snapshot = snapshots.ReceiptSnapshot(verifier)
    snapshot.verify(receipts)
    path = directory / ('models.json' if target == 'member' else 'receipt.json')
    info, raw = path.stat(), path.read_bytes()
    if mutation == 'rewrite':
        path.write_bytes(bytes([raw[0] ^ 1]) + raw[1:])
        os.utime(path, ns=(info.st_atime_ns, info.st_mtime_ns))
    elif mutation == 'replace':
        replacement = path.with_suffix('.replacement')
        replacement.write_bytes(raw)
        os.utime(replacement, ns=(info.st_atime_ns, info.st_mtime_ns))
        replacement.replace(path)
    elif mutation == 'remove':
        path.unlink()
    elif mutation == 'permission':
        path.chmod(0o600 if info.st_mode & 0o077 else 0o644)
    else:
        replacement = path.with_suffix('.target')
        replacement.write_bytes(raw)
        path.unlink()
        path.symlink_to(replacement)
    with pytest.raises(snapshots.ReceiptChanged):
        snapshot.verify(receipts)
    with pytest.raises(snapshots.ReceiptChanged, match='previously changed'):
        snapshot.verify(receipts)
    assert len(calls) == 1
    snapshot.close()
    with pytest.raises(ValueError, match='closed'):
        snapshot.verify(receipts)


@pytest.mark.parametrize('mutation', ['add_nested', 'replace_directory', 'parent_symlink'])
def test_publication_and_ancestor_identity_changes_are_fenced(tmp_path, mutation):
    directory, _key, receipts = publication(tmp_path)
    snapshot = snapshots.ReceiptSnapshot(rescue.verified)
    snapshot.verify(receipts)
    if mutation == 'add_nested':
        (directory / 'intervals/unreceipted.txt').write_text('added')
    elif mutation == 'replace_directory':
        old = directory.with_name(directory.name + '.previous')
        directory.rename(old)
        shutil.copytree(old, directory)
    else:
        parent = directory.parent
        target = parent.with_name('same_contents')
        parent.rename(target)
        parent.symlink_to(target, target_is_directory=True)
    with pytest.raises(snapshots.ReceiptChanged):
        snapshot.check(receipts)
    snapshot.close()


def test_per_dependency_reuse_and_all_used_exit_fence(tmp_path, monkeypatch):
    a, _key, first = publication(tmp_path, 'Species_first')
    b, _key, second = publication(tmp_path, 'Species_second')
    seen = []
    original = snapshots._identity

    def observed(path, **kwargs):
        seen.append(path)
        return original(path, **kwargs)

    with snapshots.ReceiptSnapshot(rescue.verified) as snapshot:
        snapshot.verify(first | second)
        monkeypatch.setattr(snapshots, '_identity', observed)
        snapshot.verify(first)
        assert any(path.is_relative_to(a) for path in seen)
        assert not any(path.is_relative_to(b) for path in seen)
    seen.clear()
    with pytest.raises(snapshots.ReceiptChanged):
        with snapshots.ReceiptSnapshot(rescue.verified) as snapshot:
            snapshot.verify(first | second)
            (b / 'models.json').write_text('changed only at command exit')
    assert snapshot.closed


def test_receipt_changed_to_malformed_json_between_hash_and_parse_is_not_verified(tmp_path, monkeypatch):
    directory, _key, receipts = publication(tmp_path)
    original = Path.read_bytes
    path = directory / 'receipt.json'

    def changed_read(self):
        if self == path:
            self.write_text('{malformed replacement')
        return original(self)

    monkeypatch.setattr(Path, 'read_bytes', changed_read)
    snapshot = snapshots.ReceiptSnapshot(lambda *_args: pytest.fail('Changed receipt must not reach verifier'))
    with pytest.raises(snapshots.ReceiptChanged):
        snapshot.verify(receipts)
    assert snapshot.entries == {}
    snapshot.close()


def test_member_change_during_initial_full_hash_cannot_be_bound_as_verified(tmp_path, monkeypatch):
    directory, _key, receipts = publication(tmp_path)
    original = rescue.digest
    member = directory / 'models.json'

    def changed_hash(path):
        result = original(path)
        if Path(path) == member:
            member.write_text('changed after its old content hash')
        return result

    monkeypatch.setattr(rescue, 'digest', changed_hash)
    snapshot = snapshots.ReceiptSnapshot(rescue.verified)
    with pytest.raises(snapshots.ReceiptChanged):
        snapshot.verify(receipts)
    assert snapshot.entries == {}
    snapshot.close()


def test_observed_receipt_hash_failure_poisons_command_without_retry(tmp_path, monkeypatch):
    directory, _key, receipts = publication(tmp_path)
    original = snapshots.digest
    path = directory / 'receipt.json'
    raw = path.read_bytes()

    def removed_receipt(filename):
        if Path(filename) == path:
            path.unlink()
        return original(filename)

    monkeypatch.setattr(snapshots, 'digest', removed_receipt)
    snapshot = snapshots.ReceiptSnapshot(lambda *_args: pytest.fail('Removed receipt reached verifier'))
    with pytest.raises(snapshots.ReceiptChanged, match='verification failed'):
        snapshot.verify(receipts)
    path.write_bytes(raw)
    monkeypatch.setattr(snapshots, 'digest', original)
    with pytest.raises(snapshots.ReceiptChanged, match='previously changed'):
        snapshot.verify(receipts)
    snapshot.close()


@pytest.mark.parametrize('name', ['', '.', '..', '../outside.fa', '/outside.fa', 'receipt.json'])
def test_malformed_receipt_paths_are_rejected_before_reading_members(tmp_path, name):
    directory, key, _receipts = publication(tmp_path)
    rescue.atomic_json(directory / 'receipt.json', {'key': key, 'files': {name: 'bad'}})
    receipts = {str(directory / 'receipt.json'): rescue.digest(directory / 'receipt.json')}
    with snapshots.ReceiptSnapshot(lambda *_args: pytest.fail('Unsafe member reached verifier')) as snapshot:
        with pytest.raises(ValueError, match='Unsafe'):
            snapshot.verify(receipts)


def stage_fixture(tmp_path, monkeypatch):
    directory, _key, receipts = publication(tmp_path)
    root = tmp_path / 'refinement'
    root.mkdir()
    rescue.atomic_json(root / 'plan.json', {'request': {'revision_candidates': {
        'Species_target': str(directory / 'revision_candidates.json')}}})
    monkeypatch.setattr(refinement, 'load', lambda *_args: None)
    key = {'dependencies': {'rescue_workers': {'Species_target': next(iter(receipts.values()))}}}
    return root, directory, key


def test_output_repair_unchanged_but_cached_dependency_change_refuses_publication(tmp_path, monkeypatch):
    root, dependency, key = stage_fixture(tmp_path, monkeypatch)
    calls = []

    def builder(tmp):
        calls.append(tmp)
        rescue.atomic_json(tmp / 'predictions.json', [{'gene_id': 'existing', 'status': 'accepted'}])

    with refinement.invocation_context():
        output = refinement.stage(root, 'prediction_fixture', key, builder)
    original = (output / 'predictions.json').read_bytes()
    (output / 'predictions.json').write_text('initial damaged output')
    with pytest.raises(snapshots.ReceiptChanged):
        with refinement.invocation_context():
            repaired = refinement.stage(root, 'prediction_fixture', key, builder)
            assert (repaired / 'predictions.json').read_bytes() == original
            (dependency / 'models.json').write_text('changed after successful proof')
            refinement.stage(root, 'second_output', key, builder)
    assert len(calls) == 2
    assert not (root / 'second_output').exists()


def test_mutation_after_output_hash_refuses_publication_and_keeps_diagnostics(tmp_path, monkeypatch):
    root, dependency, key = stage_fixture(tmp_path, monkeypatch)
    original = rescue.hash_outputs

    def changing_hashes(directory, **kwargs):
        result = original(directory, **kwargs)
        (dependency / 'models.json').write_text('changed in hash-to-rename gap')
        return result

    monkeypatch.setattr(rescue, 'hash_outputs', changing_hashes)
    with pytest.raises(snapshots.ReceiptChanged):
        with refinement.invocation_context():
            refinement.stage(root, 'unpublished', key, lambda tmp: (tmp / 'result.json').write_text('result'))
    assert not (root / 'unpublished').exists()
    assert (root / 'unpublished.failed/result.json').read_text() == 'result'
    assert refinement._INVOCATION_CACHE is None
    assert refinement._RECEIPT_SNAPSHOT is None


def test_context_cleanup_keeps_library_calls_ordinary(tmp_path, monkeypatch):
    root, _dependency, key = stage_fixture(tmp_path, monkeypatch)
    original = rescue.verified
    calls = []

    def counted(directory, current_key, **kwargs):
        calls.append(directory)
        return original(directory, current_key, **kwargs)

    monkeypatch.setattr(rescue, 'verified', counted)
    with pytest.raises(RuntimeError, match='abort'):
        with refinement.invocation_context():
            refinement.stage(root, 'first', key, lambda tmp: (tmp / 'result').write_text('same'))
            raise RuntimeError('abort')
    assert refinement._INVOCATION_CACHE is None and refinement._RECEIPT_SNAPSHOT is None
    calls.clear()
    refinement.stage(root, 'first', key, lambda _tmp: pytest.fail('Completed output must reuse'))
    refinement.stage(root, 'first', key, lambda _tmp: pytest.fail('Completed output must reuse'))
    assert sum(Path(path).parent.name == 'rescued' for path in calls) == 2


def test_verified_resume_does_not_decode_unused_catalog_loci_or_edges(tmp_path, monkeypatch):
    fixtures = import_module('test_gene_model_refinement')
    inputs, edges, _rows = fixtures.tiny_inputs(tmp_path)
    root = tmp_path / 'run'
    value = refinement.plan(root, inputs=inputs, edges=edges, mode='off')
    with refinement.invocation_context():
        refinement.select(root, value)
        refinement.predict_species(root, value, 'Species_target')
    expected = (root / 'predictions/Species_target/predictions.json').read_bytes()
    original = Path.read_text

    def unused_read(path, *args, **kwargs):
        if path.name in {'catalog_metadata.json', 'edges.json'}:
            pytest.fail('Verified resume must not deserialize unused large data')
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, 'read_text', unused_read)
    monkeypatch.setattr(refinement, 'iter_loci', lambda *_args: pytest.fail('Verified resume must not load loci'))
    with refinement.invocation_context():
        refinement.select(root, value)
        refinement.predict_species(root, value, 'Species_target')
    assert (root / 'predictions/Species_target/predictions.json').read_bytes() == expected


def test_whole_refinement_completed_resume_keeps_every_output_and_receipt_byte(tmp_path):
    fixtures = import_module('test_gene_model_refinement')
    inputs, edges, _rows = fixtures.tiny_inputs(tmp_path)
    root = tmp_path / 'run'
    value = refinement.plan(root, inputs=inputs, edges=edges, mode='off')
    refinement.finalize(root, value)
    expected = {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in root.rglob('*') if path.is_file()}
    with refinement.invocation_context():
        refinement.finalize(root, value)
    assert expected == {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                        for path in root.rglob('*') if path.is_file()}
    refinement.verify_inputs(root / 'effective/inputs.tsv')


def correspondence_fixture(tmp_path, edges, *, legacy=False, projection=None):
    directory = tmp_path / 'correspondence'
    directory.mkdir()
    rescue.atomic_json(directory / 'edges.json', edges)
    if not legacy:
        rescue.atomic_json(directory / 'donor_species.json', projection if projection is not None
                           else refinement.correspondence_donor_species(edges))
    rescue.atomic_json(directory / 'receipt.json', {'key': {'plan': 'fixture'},
                                                   'files': rescue.hash_outputs(directory)})
    return directory


@pytest.mark.parametrize('legacy', [False, True])
@pytest.mark.parametrize('edges', [[],
    [dict(species_a='A', species_b='B', ambiguous=True)],
    [dict(species_a='A', species_b='B', ambiguous=False)] * 2,
    [dict(species_a='A', species_b='B', ambiguous=False),
     dict(species_a='C', species_b='A', ambiguous=False),
     dict(species_a='A', species_b='D', ambiguous=True)]])
def test_donor_projection_and_legacy_scan_keep_exact_existing_source_guards(tmp_path, legacy, edges):
    directory = correspondence_fixture(tmp_path, edges, legacy=legacy)
    names = ['A', 'B', 'C', 'D', 'target_only']
    with refinement.invocation_context():
        for name in names:
            expected = {name}
            for edge in edges:
                if not edge['ambiguous'] and name in {edge['species_a'], edge['species_b']}:
                    expected.update((edge['species_a'], edge['species_b']))
            assert refinement.prediction_donor_names(directory, name, names) == sorted(expected)


def test_producer_projects_only_final_ambiguous_flags(tmp_path):
    fixtures = import_module('test_gene_model_refinement')
    inputs, edge_table, _rows = fixtures.tiny_inputs(tmp_path)
    rows = refinement.read_table(edge_table)
    # A conflicting target copy renders both donor links ambiguous.
    target = next(row for row in refinement.read_table(inputs) if row['species'] == 'Species_target')
    gff = Path(target['gff'])
    gff.write_text(gff.read_text() + 'chr1\ts\tgene\t1\t12\t.\t+\t.\tID=g2\n'
                   'chr1\ts\tmRNA\t1\t12\t.\t+\t.\tID=t_extra;Parent=g2\n'
                   'chr1\ts\tCDS\t1\t12\t.\t+\t0\tID=c_extra;Parent=t_extra\n')
    fasta = Path(target['cds'])
    fasta.write_text(fasta.read_text() + '>g2\nATGAAACCCTAA\n')
    conflicting = {**rows[0], 'gene_a': 'Species_target_g2'}
    fixtures.write_tsv(edge_table, list(rows[0]), rows + [conflicting])
    root = tmp_path / 'run'
    value = refinement.plan(root, inputs=inputs, edges=edge_table, mode='off')
    directory = refinement.correspondence(root, value)
    final = json.loads((directory / 'edges.json').read_text())
    assert any(edge['ambiguous'] for edge in final)
    assert json.loads((directory / 'donor_species.json').read_text()) == refinement.correspondence_donor_species(final)
    assert 'donor_species.json' in json.loads((directory / 'receipt.json').read_text())['files']


def test_unreceipted_projection_is_not_used_and_legacy_still_scans_edges(tmp_path):
    edges = [dict(species_a='A', species_b='B', ambiguous=False)]
    directory = correspondence_fixture(tmp_path, edges, legacy=True)
    rescue.atomic_json(directory / 'donor_species.json', {'A': ['C']})
    with refinement.invocation_context():
        assert refinement.prediction_donor_names(directory, 'A', ['A', 'B', 'C']) == ['A', 'B']


@pytest.mark.parametrize('projection', [[], {'A': ['C']}, {'A': ['B', 'B']}, {'A': [False]}, {'C': []}])
def test_invalid_receipted_projection_is_not_a_legacy_fallback(tmp_path, projection):
    directory = correspondence_fixture(tmp_path, [], projection=projection)
    with refinement.invocation_context():
        with pytest.raises(ValueError, match='Invalid correspondence donor'):
            refinement.prediction_donor_names(directory, 'A', ['A', 'B'])


def test_changed_projection_refuses_cached_generation_and_publication(tmp_path):
    directory = correspondence_fixture(tmp_path, [dict(species_a='A', species_b='B', ambiguous=False)])
    with pytest.raises(snapshots.ReceiptChanged):
        with refinement.invocation_context():
            assert refinement.prediction_donor_names(directory, 'A', ['A', 'B']) == ['A', 'B']
            rescue.atomic_json(directory / 'donor_species.json', {})
            refinement.prediction_donor_names(directory, 'A', ['A', 'B'])
