"""CLI contracts for relation strings using hand-known branch ranks and orders."""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = next(parent for parent in Path(__file__).resolve().parents
            if (parent / 'workflow/support/orthogroup_statistics.py').is_file())
SCRIPT = ROOT / 'workflow/support/orthogroup_statistics.py'


@pytest.mark.parametrize('newick, expected', [
    ('((a:1,b:1)x:1,(c:1,(d:1,e:1)z:1)y:1)root;', {
        'a': ('0', '2,8', ''), 'b': ('1', '2,8', ''), 'x': ('2', '8', '0,1'),
        'c': ('3', '7,8', ''), 'd': ('4', '6,7,8', ''), 'e': ('5', '6,7,8', ''),
        'z': ('6', '7,8', '4,5'), 'y': ('7', '8', '3,6,4,5'),
        'root': ('8', '', '2,7,0,1,3,6,4,5'),
    }),
    ('((b:1,a:1)x:1,((e:1,d:1)z:1,c:1)y:1)root;', {
        'a': ('0', '2,8', ''), 'b': ('1', '2,8', ''), 'x': ('2', '8', '1,0'),
        'c': ('3', '7,8', ''), 'd': ('4', '6,7,8', ''), 'e': ('5', '6,7,8', ''),
        'z': ('6', '7,8', '5,4'), 'y': ('7', '8', '6,3,5,4'),
        'root': ('8', '', '2,7,1,0,6,3,5,4'),
    }),
    ('(NA:1,001:1,NULL:1,nan:1)root;', {
        '001': ('0', '4', ''), 'NA': ('1', '4', ''), 'NULL': ('2', '4', ''),
        'nan': ('3', '4', ''), 'root': ('4', '', '1,0,2,3'),
    }),
    ('a:1;', {'a': ('0', '', '')}),
    ('((a:1,b:2,c:3)x:4,d:5)root;', {
        'a': ('0', '3,5', ''), 'b': ('1', '3,5', ''), 'c': ('2', '3,5', ''),
        'x': ('3', '5', '0,1,2'), 'd': ('4', '5', ''), 'root': ('5', '', '3,4,0,1,2'),
    }),
    ('((a:1)x:2,b:3)root;', {
        'a': ('1', '0,3', ''), 'x': ('0', '3', '1'), 'b': ('2', '3', ''),
        'root': ('3', '', '0,2,1'),
    }),
])
def test_cli_preserves_nearest_ancestor_and_level_order_descendant_strings(tmp_path, newick, expected):
    tree = tmp_path / 'rooted.nwk'
    tree.write_text(newick + '\n')
    original = tree.read_bytes()
    # The subprocess loads actual container dependencies, without package stubs.
    completed = subprocess.run([sys.executable, str(SCRIPT), '--rooted_tree', str(tree)],
                               cwd=tmp_path, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    result = pd.read_csv(tmp_path / 'orthogroup.branch.tsv', sep='\t', keep_default_na=False, dtype=str)
    observed = {row.node_name: (row.branch_id, row.ancestors, row.descendants)
                for row in result.itertuples(index=False)}
    assert observed == expected
    assert len(result) == len(expected)
    for row in result.itertuples(index=False):
        siblings = result.loc[(result.parent == row.parent) & (result.branch_id != row.branch_id)]
        expected_sister = siblings.iloc[0].branch_id if row.parent != '-999' and len(siblings) == 1 else '-999'
        assert row.sister == expected_sister
    summary = pd.read_csv(tmp_path / 'orthogroup.tree.tsv', sep='\t')
    assert len(summary) == 1
    assert summary.loc[0, 'num_branch'] == len(expected)
    assert tree.read_bytes() == original
