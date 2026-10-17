import importlib
import sys
from pathlib import Path

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'support'))
eq=importlib.import_module('annevo_prediction_equivalence')


def test_comments_and_order_do_not_change_any_annotation_field(tmp_path):
    a,b=tmp_path/'a',tmp_path/'b'
    lines=['seq\tANNEVO\tgene\t1\t9\t0.9\t+\t.\tID=g',
           'seq\tANNEVO\tCDS\t1\t9\t0.9\t+\t0\tID=c;Parent=g']
    a.write_text('# sequence 1\n'+'\n'.join(lines)+'\n')
    b.write_text('# sequence 2\n'+'\n'.join(reversed(lines))+'\n')
    assert eq.matching_features(a,b)
    b.write_text(b.read_text().replace('+\t0','+\t1'))
    with pytest.raises(ValueError,match='phase'):
        eq.matching_features(a,b)


def test_only_execution_artifact_fields_are_excluded_from_decision_comparison():
    value={'receipt_sha256':'execution','GFF_sha256':'comments','GFF_bytes':10,
           'model_sha256':'model','runtime_digest':'runtime','blocks':[[1,9,0]],
           'chosen':True,'alignment':{'identity':.99}}
    result=eq.decision_fields(value)
    assert result=={k:v for k,v in value.items() if k not in eq.ARTIFACT_FIELDS}
    assert result['blocks']==[[1,9,0]] and result['model_sha256']=='model'
