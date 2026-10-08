"""Immutable input generations, exact scopes, aliases and publication races."""
import copy
import hashlib
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
snapshots = import_module('refinement_input_snapshot')
state = import_module('input_generation_array_state')
fixtures = import_module('test_gene_model_refinement')


def proof_fixture(tmp_path):
    root = tmp_path / 'run'
    root.mkdir()
    sources, files = {}, {}
    for name in ('Species_a', 'Species_b'):
        sources[name] = {}
        for field in ('fasta', 'gff', 'genome'):
            path = tmp_path / (name + '.' + field)
            path.write_text(name + field + '\n')
            sources[name][field] = str(path)
            files[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    other = tmp_path / 'global_revision.json'
    other.write_text('{"records": []}\n')
    files[str(other)] = hashlib.sha256(other.read_bytes()).hexdigest()
    value = {'request': {'sources': sources, 'files': files, 'parameters': {'policy': 'conserved'}},
             'species': list(sources)}
    state.atomic_json(root / 'plan.json', value)
    return root, value, other


def select_files(value, names=None):
    request = value['request']
    paths = {s[k] for s in request['sources'].values() for k in ('fasta', 'gff', 'genome')}
    files = request['files'] if names is None else {p: sha for p, sha in request['files'].items() if p not in paths}
    if names is not None:
        files.update({request['sources'][name][k]: request['files'][request['sources'][name][k]]
                      for name in names for k in ('fasta', 'gff', 'genome')})
    return files


def bind(snapshot, root, value, files=None):
    snapshot.verify_plan(root, value)
    snapshot.verify_files(value['request']['files'] if files is None else files)


def mutate(path, kind):
    info, raw = path.stat(), path.read_bytes()
    if kind == 'rewrite':
        path.write_bytes(bytes([raw[0] ^ 1]) + raw[1:])
        os.utime(path, ns=(info.st_atime_ns, info.st_mtime_ns))
    elif kind == 'replace':
        replacement = path.with_suffix('.replacement')
        replacement.write_bytes(raw)
        os.utime(replacement, ns=(info.st_atime_ns, info.st_mtime_ns))
        replacement.replace(path)
    elif kind == 'remove':
        path.unlink()
    elif kind == 'permission':
        path.chmod(0o600 if info.st_mode & 0o077 else 0o644)
    elif kind == 'symlink':
        replacement = path.with_suffix('.target')
        replacement.write_bytes(raw)
        path.unlink()
        path.symlink_to(replacement)
    return raw


@pytest.mark.parametrize('target', ['input', 'plan'])
@pytest.mark.parametrize('kind', ['rewrite', 'replace', 'remove', 'permission', 'symlink'])
def test_generation_changes_are_sticky_even_after_restoration(tmp_path, target, kind):
    root, value, other = proof_fixture(tmp_path)
    snapshot = snapshots.RefinementInputSnapshot()
    bind(snapshot, root, value)
    path = other if target == 'input' else root / 'plan.json'
    raw = mutate(path, kind)
    with pytest.raises(snapshots.InputChanged):
        bind(snapshot, root, value)
    if path.is_symlink():
        path.unlink()
    path.write_bytes(raw)
    with pytest.raises(snapshots.InputChanged, match='previously changed'):
        bind(snapshot, root, value)
    snapshot.close()
    with pytest.raises(ValueError, match='closed'):
        bind(snapshot, root, value)


@pytest.mark.parametrize('kind', ['parent_permission', 'parent_replacement', 'parent_symlink'])
def test_ancestor_generation_and_permission_changes_are_rejected(tmp_path, kind):
    container = tmp_path / 'inputs'
    container.mkdir()
    root, value, _other = proof_fixture(container)
    snapshot = snapshots.RefinementInputSnapshot()
    bind(snapshot, root, value)
    if kind == 'parent_permission':
        container.chmod(0o700)
    else:
        old = container.with_name('previous_inputs')
        container.rename(old)
        if kind == 'parent_replacement':
            shutil.copytree(old, container)
        else:
            container.symlink_to(old, target_is_directory=True)
    with pytest.raises(snapshots.InputChanged):
        snapshot.check_all()
    snapshot.close()


def test_unrelated_sibling_generation_is_permitted(tmp_path):
    root, value, _other = proof_fixture(tmp_path)
    with snapshots.RefinementInputSnapshot() as snapshot:
        bind(snapshot, root, value)
        (tmp_path / 'new_output_sibling').mkdir()
        (tmp_path / 'new_output_sibling/receipt.json').write_text('{}\n')
        bind(snapshot, root, value)


@pytest.mark.parametrize('changed_field', ['parameters', 'species', 'files'])
def test_plan_semantics_cannot_be_absorbed_while_builder_keeps_old_value(tmp_path, changed_field):
    root, value, _other = proof_fixture(tmp_path)
    snapshot = snapshots.RefinementInputSnapshot()
    bind(snapshot, root, value)
    replacement = copy.deepcopy(value)
    if changed_field == 'parameters':
        replacement['request']['parameters']['policy'] = 'longest'
    elif changed_field == 'species':
        replacement['species'].reverse()
    else:
        replacement['request']['files'].pop(next(iter(replacement['request']['files'])))
    state.atomic_json(root / 'plan.json', replacement)
    with pytest.raises(snapshots.InputChanged):
        bind(snapshot, root, replacement)
    snapshot.close()


def test_frozen_plan_deepcopy_cannot_follow_in_memory_mutation(tmp_path):
    root, value, _other = proof_fixture(tmp_path)
    snapshot = snapshots.RefinementInputSnapshot()
    bind(snapshot, root, value)
    value['request']['parameters']['policy'] = 'longest'
    with pytest.raises(snapshots.InputChanged):
        snapshot.check_all()
    snapshot.close()


def test_midhash_change_is_not_bound_as_a_completed_proof(tmp_path, monkeypatch):
    root, value, other = proof_fixture(tmp_path)
    original = state.digest

    def changed(path):
        result = original(path)
        if Path(path) == other:
            mutate(other, 'rewrite')
        return result

    monkeypatch.setattr(state, 'digest', changed)
    snapshot = snapshots.RefinementInputSnapshot()
    with pytest.raises(snapshots.InputChanged):
        bind(snapshot, root, value)
    with pytest.raises(snapshots.InputChanged, match='previously changed'):
        snapshot.check_all()
    snapshot.close()


@pytest.mark.parametrize('link', ['hardlink', 'symlink'])
def test_shared_physical_aliases_read_bytes_once_with_independent_path_fences(tmp_path, monkeypatch, link):
    root, value, other = proof_fixture(tmp_path)
    target = Path(value['request']['sources']['Species_a']['fasta'])
    alias = Path(value['request']['sources']['Species_b']['fasta'])
    alias.unlink()
    if link == 'hardlink':
        os.link(target, alias)
    else:
        alias.symlink_to(target)
    value['request']['files'][str(alias)] = value['request']['files'][str(target)]
    state.atomic_json(root / 'plan.json', value)
    original, reads = state.digest, Counter()

    def counted(path):
        info = Path(path).stat()
        reads[info.st_dev, info.st_ino] += 1
        return original(path)

    monkeypatch.setattr(state, 'digest', counted)
    snapshot = snapshots.RefinementInputSnapshot()
    bind(snapshot, root, value, select_files(value, ['Species_a']))
    bind(snapshot, root, value, select_files(value, ['Species_b']))
    info = target.stat()
    assert reads[info.st_dev, info.st_ino] == 1
    mutate(alias, 'replace')
    with pytest.raises(snapshots.InputChanged):
        bind(snapshot, root, value, select_files(value, ['Species_b']))
    assert other.is_file()
    snapshot.close()


def test_alias_reuse_cannot_ignore_changed_original_path(tmp_path):
    root, value, _other = proof_fixture(tmp_path)
    target = Path(value['request']['sources']['Species_a']['fasta'])
    alias = Path(value['request']['sources']['Species_b']['fasta'])
    alias.unlink()
    alias.symlink_to(target)
    value['request']['files'][str(alias)] = value['request']['files'][str(target)]
    state.atomic_json(root / 'plan.json', value)
    snapshot = snapshots.RefinementInputSnapshot()
    bind(snapshot, root, value, select_files(value, ['Species_a']))
    original = target.with_suffix('.original')
    target.rename(original)
    target.symlink_to(original)
    with pytest.raises(snapshots.InputChanged):
        bind(snapshot, root, value, select_files(value, ['Species_b']))
    snapshot.close()


def test_new_alias_cannot_absorb_identical_byte_replacement_with_new_inode(tmp_path):
    root, value, _other = proof_fixture(tmp_path)
    target = Path(value['request']['sources']['Species_a']['fasta'])
    alias = Path(value['request']['sources']['Species_b']['fasta'])
    alias.unlink()
    alias.symlink_to(target)
    value['request']['files'][str(alias)] = value['request']['files'][str(target)]
    state.atomic_json(root / 'plan.json', value)
    snapshot = snapshots.RefinementInputSnapshot()
    bind(snapshot, root, value, select_files(value, ['Species_a']))
    mutate(target, 'replace')
    with pytest.raises(snapshots.InputChanged):
        bind(snapshot, root, value, select_files(value, ['Species_b']))
    snapshot.close()


def test_independent_files_with_identical_content_do_not_share_a_stat_proof(tmp_path, monkeypatch):
    root, value, _other = proof_fixture(tmp_path)
    target = Path(value['request']['sources']['Species_a']['fasta'])
    independent = Path(value['request']['sources']['Species_b']['fasta'])
    independent.write_bytes(target.read_bytes())
    value['request']['files'][str(independent)] = value['request']['files'][str(target)]
    state.atomic_json(root / 'plan.json', value)
    original, reads = state.digest, Counter()

    def counted(path):
        reads[Path(path)] += 1
        return original(path)

    monkeypatch.setattr(state, 'digest', counted)
    with snapshots.RefinementInputSnapshot() as snapshot:
        bind(snapshot, root, value, select_files(value, ['Species_a']))
        bind(snapshot, root, value, select_files(value, ['Species_b']))
    assert reads[target] == reads[independent] == 1


def test_many_global_inputs_keep_generation_checks_bounded_per_guard(tmp_path, monkeypatch):
    root, value, _other = proof_fixture(tmp_path)
    for index in range(40):
        path = tmp_path / ('revision_' + str(index) + '.json')
        path.write_text('{"records": []}\n')
        value['request']['files'][str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    state.atomic_json(root / 'plan.json', value)
    snapshot = snapshots.RefinementInputSnapshot()
    bind(snapshot, root, value, select_files(value, []))
    global_proof = snapshot.groups[None]
    original, checks = global_proof.check, []

    def checked():
        checks.append(True)
        return original()

    monkeypatch.setattr(global_proof, 'check', checked)
    snapshot.verify_files(select_files(value, []))
    # Full owner and pathname checks must stay independent of the number of
    # members, rather than scanning all members once per member in a guard.
    assert len(checks) <= 3
    mutate(tmp_path / 'revision_39.json', 'replace')
    with pytest.raises(snapshots.InputChanged):
        snapshot.verify_files(select_files(value, []))
    snapshot.close()


def real_plan(tmp_path):
    inputs, edges, rows = fixtures.tiny_inputs(tmp_path)
    root = tmp_path / 'refinement'
    value = refinement.plan(root, inputs=inputs, edges=edges, mode='off')
    return root, value, rows


def request_metrics(monkeypatch, value):
    targets = {Path(p).resolve() for p in value['request']['files']}
    measured = Counter()
    original = state.digest

    def counted(path):
        if Path(path).resolve() in targets:
            measured['reads'] += 1
            measured['bytes'] += Path(path).stat().st_size
        return original(path)

    monkeypatch.setattr(state, 'digest', counted)
    return measured


def test_four_load_guards_reduce_request_bytes_and_every_new_cli_is_fresh(tmp_path, monkeypatch):
    root, value, _rows = real_plan(tmp_path)
    measured = request_metrics(monkeypatch, value)
    for _ in range(4):
        refinement.load(root)
    ordinary = dict(measured)
    measured.clear()
    with refinement.invocation_context():
        for _ in range(4):
            assert refinement.load(root) == value
    local = dict(measured)
    assert ordinary == {key: amount * 4 for key, amount in local.items()}
    assert local['bytes'] == sum(Path(p).stat().st_size for p in value['request']['files'])
    measured.clear()
    with refinement.invocation_context():
        refinement.load(root)
    assert dict(measured) == local
    measured.clear()
    refinement.load(root)
    assert dict(measured) == local
    assert refinement._INPUT_SNAPSHOT is None


@pytest.mark.parametrize('names', [[], ['Species_target'], ['Species_target', 'Species_donor1']])
def test_names_selection_preserves_global_files_and_omits_other_species(tmp_path, monkeypatch, names):
    root, value, _rows = real_plan(tmp_path)
    expected = select_files(value, names)
    omitted = Path(value['request']['sources']['Species_donor2']['genome'])
    omitted.write_text('not a valid frozen input, but outside this scope\n')
    measured = request_metrics(monkeypatch, value)
    with refinement.invocation_context():
        for _ in range(2):
            refinement.load(root, iter(names))
    assert measured['reads'] == len(expected)
    assert measured['bytes'] == sum(Path(p).stat().st_size for p in expected)
    with pytest.raises(ValueError, match='Frozen refinement input changed'):
        refinement.load(root)


def test_runtime_identity_checks_stay_ordinary_every_load(tmp_path, monkeypatch):
    root, value, _rows = real_plan(tmp_path)
    tool = tmp_path / 'miniprot'
    tool.write_bytes(b'fixture tool identity\n')
    value['request']['miniprot'] = {'path': str(tool), 'sha256': refinement.digest(tool)}
    state.atomic_json(root / 'plan.json', value)
    measured = Counter()
    for name in ('implementation', 'dependency_identities'):
        original = getattr(refinement, name)

        def counted(original=original, name=name):
            measured[name] += 1
            return original()

        monkeypatch.setattr(refinement, name, counted)
    original_digest = refinement.digest

    def counted_digest(path):
        if Path(path) == tool:
            measured['tool_reads'] += 1
        return original_digest(path)

    monkeypatch.setattr(refinement, 'digest', counted_digest)
    with refinement.invocation_context():
        for _ in range(4):
            refinement.load(root)
    assert measured == {'implementation': 4, 'dependency_identities': 4, 'tool_reads': 4}


@pytest.mark.parametrize('change', ['input', 'plan'])
def test_post_output_hash_change_blocks_publication_and_preserves_old_result(tmp_path, monkeypatch, change):
    root, value, _rows = real_plan(tmp_path)
    name = 'Species_target'
    refinement.stage(root, 'result', {'generation': 'old'}, lambda tmp: (tmp / 'result').write_text('old'), [name])
    before = {p.name: p.read_bytes() for p in (root / 'result').iterdir()}
    original = refinement.rescue.hash_outputs

    def changed_hashes(*args, **kwargs):
        output = original(*args, **kwargs)
        if change == 'input':
            mutate(Path(value['request']['sources'][name]['gff']), 'rewrite')
        else:
            changed = copy.deepcopy(value)
            changed['request']['parameters']['policy'] = 'longest'
            state.atomic_json(root / 'plan.json', changed)
        return output

    monkeypatch.setattr(refinement.rescue, 'hash_outputs', changed_hashes)
    with pytest.raises(snapshots.InputChanged):
        with refinement.invocation_context():
            refinement.stage(root, 'result', {'generation': 'new'}, lambda tmp: (tmp / 'result').write_text('new'), [name])
    assert before == {p.name: p.read_bytes() for p in (root / 'result').iterdir()}
    assert (root / 'result.failed/result').read_text() == 'new'
    assert refinement._INPUT_SNAPSHOT is None


def test_successful_exit_checks_every_used_scope_and_closes(tmp_path):
    root, value, _rows = real_plan(tmp_path)
    with pytest.raises(snapshots.InputChanged):
        with refinement.invocation_context():
            refinement.load(root)
            snapshot = refinement._INPUT_SNAPSHOT
            mutate(Path(value['request']['sources']['Species_donor2']['genome']), 'rewrite')
    assert snapshot.closed and refinement._INPUT_SNAPSHOT is None


def test_malformed_plan_read_cannot_recover_inside_command(tmp_path):
    root, _value, _rows = real_plan(tmp_path)
    raw = (root / 'plan.json').read_bytes()
    with pytest.raises(snapshots.InputChanged, match='previously changed'):
        with refinement.invocation_context():
            refinement.load(root)
            (root / 'plan.json').write_text('{malformed')
            with pytest.raises(snapshots.InputChanged):
                refinement.load(root)
            (root / 'plan.json').write_bytes(raw)
            refinement.load(root)
    assert refinement._INPUT_SNAPSHOT is None


def test_nested_context_restores_outer_proof_without_cross_root_reuse(tmp_path):
    first = tmp_path / 'first'
    second = tmp_path / 'second'
    first.mkdir()
    second.mkdir()
    a, _value, _rows = real_plan(first)
    b, _value, _rows = real_plan(second)
    with refinement.invocation_context():
        refinement.load(a)
        outer = refinement._INPUT_SNAPSHOT
        with refinement.invocation_context():
            refinement.load(b)
        assert refinement._INPUT_SNAPSHOT is outer
        refinement.load(a)
    assert outer.closed and refinement._INPUT_SNAPSHOT is None


def test_completed_refinement_real_cli_resume_preserves_all_output_and_receipt_bytes(tmp_path, monkeypatch):
    root, value, _rows = real_plan(tmp_path)
    refinement.finalize(root, value)
    before = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in root.rglob('*') if p.is_file()}
    monkeypatch.setattr(sys, 'argv', ['gene_model_refinement.py', 'run', '--output', str(root), '--cpus', '1'])
    refinement.main()
    assert before == {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in root.rglob('*') if p.is_file()}
    refinement.verify_inputs(root / 'effective/inputs.tsv')
