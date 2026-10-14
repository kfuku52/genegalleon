"""The comparison retains coordinates, scores, identifiers and phases."""
import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location('pipeline_benchmark', Path(__file__).resolve().parents[1]/'benchmarks/benchmark_annevo_pipeline.py')
bm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bm)


@pytest.mark.parametrize('change', ['phase', 'score', 'coordinate', 'identifier', 'metadata_order'])
def test_gff_equivalence_only_ignores_comments_and_feature_order(tmp_path, change):
    a, b = tmp_path/'a.gff3', tmp_path/'b.gff3'
    lines = ['seq\tANNEVO\tgene\t1\t12\t0.95\t+\t.\tID=g1',
             'seq\tANNEVO\tCDS\t1\t12\t0.95\t+\t0\tID=c1;Parent=g1']
    a.write_text('\n'.join(lines)+'\n')
    replacements = {'phase': ('+\t0', '+\t1'), 'score': ('0.95', '0.94'),
                    'coordinate': ('1\t12', '2\t12'), 'identifier': ('ID=c1', 'ID=c2')}
    text = '\n'.join(reversed(lines))+'\n'
    if change in replacements:
        text = text.replace(*replacements[change])
    b.write_text('# native metadata\n'+text)
    assert (bm.feature_digest(a)==bm.feature_digest(b)) == (change=='metadata_order')
