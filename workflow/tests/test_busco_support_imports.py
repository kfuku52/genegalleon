"""BUSCO preservation must not import the optional synteny toolchain."""
import os
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.parametrize("mode", ["package", "script"])
def test_busco_support_does_not_load_rescue_or_synteny_dependencies(mode, tmp_path):
    root = Path(__file__).resolve().parents[2]
    code = r'''
import importlib
import importlib.abc
from pathlib import Path
import runpy
import sys

class RejectSynteny(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[-1] in {'rescue_gene_models', 'pairwise_synteny'} or fullname.split('.')[0] in {'kffractbias', 'jcvi'}:
            raise ImportError('BUSCO imported an unrelated synteny dependency: ' + fullname)

sys.meta_path.insert(0, RejectSynteny())
if sys.argv[1] == 'package':
    module = importlib.import_module('workflow.support.busco_guide_tree')
    assert module.safe_token('BUSCO_1', 'marker') == 'BUSCO_1'
    assert module.COMPARABLE_QUALITY == ('lineage', 'version', 'mode', 'lineage_date', 'markers')
else:
    script = Path(sys.argv[2]) / 'workflow/support/busco_guide_tree.py'
    sys.path.insert(0, str(script.parent))
    sys.argv = [str(script), '--help']
    runpy.run_path(str(script), run_name='__main__')
'''
    run_cwd = tmp_path / "run_cwd"
    run_cwd.mkdir()
    result = subprocess.run([sys.executable, "-B", "-c", code, mode, str(root)],
                            cwd=run_cwd if mode == "script" else root,
                            capture_output=True, text=True,
                            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    assert result.returncode == 0, result.stdout + result.stderr
    if mode == "script":
        assert "usage:" in result.stdout.lower()
        assert list(run_cwd.iterdir()) == []
