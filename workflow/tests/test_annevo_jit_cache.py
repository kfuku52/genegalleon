import importlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'support'))
jit=importlib.import_module('annevo_jit_cache')


@pytest.fixture
def runtime(tmp_path):
    source=tmp_path/'source'
    source.mkdir()
    return dict(generation=str(source),runtime_digest='a'*64,source_commit='b'*40)


def test_cache_reuse_and_runtime_source_cpu_invalidation(tmp_path,runtime,monkeypatch):
    monkeypatch.setattr(jit,'cpu_identity',lambda:{'features':'CPU one'})
    root=tmp_path/'jit'
    first,key=jit.prepare(root,runtime,[])
    assert jit.prepare(root,runtime,[])==(first,key)
    assert json.loads((first/'cache_identity.json').read_text())['runtime_digest']==runtime['runtime_digest']
    assert jit.prepare(root,{**runtime,'runtime_digest':'c'*64},[])[0]!=first
    assert jit.prepare(root,{**runtime,'source_commit':'d'*40},[])[0]!=first
    monkeypatch.setattr(jit,'cpu_identity',lambda:{'features':'CPU two'})
    assert jit.prepare(root,runtime,[])[0]!=first


@pytest.mark.parametrize('problem',['redirected_root','writable_root','redirected_member','writable_member','changed_identity'])
def test_unsafe_or_changed_caches_stop_before_native_loading(tmp_path,runtime,problem):
    root=tmp_path/'jit'
    first,_=jit.prepare(root,runtime,[])
    if problem=='redirected_root':
        redirected=tmp_path/'redirected'
        redirected.symlink_to(root,target_is_directory=True)
        root=redirected
    elif problem=='writable_root':
        root.chmod(0o777)
    elif problem=='redirected_member':
        (first/'unsafe.nbc').symlink_to(tmp_path/'foreign')
    elif problem=='writable_member':
        member=first/'unsafe.nbc'
        member.write_bytes(b'compiled kernel')
        member.chmod(0o666)
    else:
        identity=first/'cache_identity.json'
        identity.chmod(0o600)
        identity.write_text('{}')
    with pytest.raises(ValueError,match='cache'):
        jit.prepare(root,runtime,[])


def test_cache_cannot_contain_inputs_or_modify_installed_runtime(tmp_path,runtime):
    protected=tmp_path/'input.fa'
    protected.write_text('>seq\nATG\n')
    with pytest.raises(ValueError,match='separate'):
        jit.prepare(tmp_path,runtime,[protected])
    with pytest.raises(ValueError,match='separate'):
        jit.prepare(Path(runtime['generation'])/'jit',runtime,[])
    assert not (Path(runtime['generation'])/'jit').exists()
