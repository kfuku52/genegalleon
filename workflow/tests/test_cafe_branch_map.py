import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / 'support/cafe_branch_map.py'
spec = importlib.util.spec_from_file_location('branch_map', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

OLD = '((Tri<0>_4:1,Anc<1>_3:1)<3>_4:1,Out<2>_1:2)<4>_1;'
NEW = '(((Tri<0>_4:1,Dio<1>_4:1)<4>_4:1,Anc<2>_3:2)<5>_4:1,Out<3>_1:3)<6>_1;'


def paths(tmp_path, old=OLD, new=NEW):
    a, b = tmp_path / 'old.tre', tmp_path / 'new.tre'
    a.write_text('#NEXUS\nBEGIN TREES;\n TREE F1 = ' + old + '\nEND;\n')
    b.write_text('#NEXUS\nBEGIN TREES;\n TREE F2 = ' + new + '\nEND;\n')
    return a, b


def test_added_species_are_explicit_on_unique_shared_clade(tmp_path):
    a, b = paths(tmp_path)
    result = module.map_branch(a, b, '<3>')
    assert result['eligible'] and result['new_target'] == '<5>'
    assert result['state'] == 'unique_clade_with_added_species'
    assert result['candidates'][0]['descendants'] == ['Anc', 'Dio', 'Tri']
    assert result['candidates'][0]['added_descendants'] == ['Dio']


def test_added_sister_species_does_not_silently_choose_old_edge(tmp_path):
    a, b = paths(tmp_path)
    result = module.map_branch(a, b, 'Tri<0>')
    assert not result['eligible'] and result['new_target'] is None
    assert result['state'] == 'ambiguous_edge_split'
    assert {row['native_key'] for row in result['candidates']} == {'Tri<0>', '<4>'}


def test_unchanged_terminal_and_internal_clades(tmp_path):
    a, b = paths(tmp_path, new=OLD)
    for old, expected in [('<3>', '<3>'), ('<1>', 'Anc<1>')]:
        result = module.map_branch(a, b, old)
        assert result['eligible'] and result['new_target'] == expected
        assert result['state'] == 'exact_clade'


@pytest.mark.parametrize('target', ['<4>', '<40>', 'Anc'])
def test_root_unknown_and_partial_identifiers_are_rejected(tmp_path, target):
    a, b = paths(tmp_path)
    with pytest.raises(ValueError, match='exactly one non-root'):
        module.map_branch(a, b, target)


def test_missing_target_taxon_and_changed_topology_do_not_invent_mapping(tmp_path):
    a, b = paths(tmp_path, new='(Tri<0>_1:1,Out<2>_1:1)<3>_1;')
    result = module.map_branch(a, b, '<3>')
    assert not result['eligible'] and result['state'] == 'missing_target_species'
    a, b = paths(tmp_path, new='((Tri<0>_1:1,Out<2>_1:1)<3>_1:1,Anc<1>_1:2)<4>_1;')
    result = module.map_branch(a, b, '<3>')
    assert not result['eligible'] and result['state'] == 'unmapped_clade'


@pytest.mark.parametrize('new', [
    '(Tri<0>_1:1,Tri<1>_1:1)<2>_1;',
    '(Tri<0>_1:1,Anc<0>_1:1)<2>_1;',
    '(Tri:1,Anc<1>_1:1)<2>_1;'])
def test_invalid_species_and_branch_names_stop(tmp_path, new):
    a, b = paths(tmp_path, new=new)
    with pytest.raises(ValueError):
        module.map_branch(a, b, '<3>')


def test_multiline_nexus_uses_first_tree_and_preserves_full_species_keys(tmp_path):
    name = 'Tri_qualified_v1'
    a, b = paths(tmp_path, old=OLD.replace('Tri', name), new=OLD.replace('Tri', name).replace('_3:', '*_3:'))
    b.write_text(b.read_text().replace('Anc<1>*_3:', '\nAnc<1>*_3:') + 'TREE ignored = malformed;\n')
    result = module.map_branch(a, b, '<3>')
    assert result['eligible'] and name in result['old_target']['descendants']


def test_cli_binds_table_headers_and_never_overwrites_mapping(tmp_path, capsys):
    a, b = paths(tmp_path)
    table = tmp_path / 'change.tab'
    table.write_text('FamilyID\t<5>\nF1\t-1\n')
    out = tmp_path / 'mapping.json'
    argv = ['--old-tree', str(a), '--new-tree', str(b), '--old-target', '<3>',
            '--new-branch-table', str(table), '--output', str(out)]
    assert module.main(argv) == 0
    assert 'new_branch_table_sha256' in capsys.readouterr().out
    with pytest.raises(FileExistsError):
        module.main(argv)
    table.write_text('FamilyID\t<9>\nF1\t-1\n')
    with pytest.raises(ValueError, match='native branch table'):
        module.main(argv[:-2])


def test_cli_ambiguous_mapping_has_nonzero_exit(tmp_path, capsys):
    a, b = paths(tmp_path)
    assert module.main(['--old-tree', str(a), '--new-tree', str(b), '--old-target', '<0>']) == 2
    assert 'ambiguous_edge_split' in capsys.readouterr().out
