"""Paired native routes must retain actual decisions, not just aggregate counts."""
import importlib.util
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'benchmarks'))
spec = importlib.util.spec_from_file_location('genome_mask_benchmark', Path(__file__).resolve().parents[1]/'benchmarks/benchmark_annevo_genome_models.py')
bm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bm)


def fixture(root):
    locations = []
    for scope in ('unmasked_control', 'gene_42', 'exon_42'):
        roles = ('gene', 'exon') if scope == 'unmasked_control' else ('',)
        for condition in ('with_ANNEVO', 'with_compared_ANNEVO'):
            for role in roles:
                d = root/scope/condition/role
                d.mkdir(parents=True)
                case = dict(candidate_id='c1', blocks=[[0, 9, 0]], chosen=True,
                            evidence={'receipt_sha256': condition, 'GFF_sha256': condition,
                                      'GFF_bytes': len(condition), 'model_sha256':'same', 'alignment':{'identity':.99}})
                (d/'accepted_cases.json').write_text(json.dumps([case]))
                if role == 'exon' or scope.startswith('exon_'):
                    (d/'selected_predictions.gff3').write_text('seq\tGG\tCDS\t1\t9\t.\t+\t0\tID=c1\n')
                    (d/'selected_predictions.cds.fa').write_text('>c1\nATGAAATAA\n')
                    (d/'selection.json').write_text(json.dumps({'candidate_id':'c1','score':.99}))
                locations.append(d)
    return locations


def test_independent_execution_receipt_is_the_only_ignored_case_metadata(tmp_path):
    fixture(tmp_path)
    result = bm.compare_admissions(tmp_path, [42])
    assert len(result)==4
    assert all(r['accepted_models']==1 and r['accepted_cases_identical'] for r in result)


@pytest.mark.parametrize('change', ['phase', 'alignment', 'chosen_model', 'missing_CDS'])
def test_equal_counts_do_not_hide_changed_decisions_or_incomplete_artifacts(tmp_path, change):
    fixture(tmp_path)
    d=tmp_path/'exon_42/with_compared_ANNEVO'
    if change in ('phase', 'alignment'):
        cases=json.loads((d/'accepted_cases.json').read_text())
        if change=='phase':
            cases[0]['blocks'][0][2]=1
        else:
            cases[0]['evidence']['alignment']['identity']=.98
        (d/'accepted_cases.json').write_text(json.dumps(cases))
    elif change=='chosen_model':
        (d/'selection.json').write_text(json.dumps({'candidate_id':'c2','score':.99}))
    else:
        (d/'selected_predictions.cds.fa').unlink()
    with pytest.raises(ValueError, match='decisions or selected artifacts changed'):
        bm.compare_admissions(tmp_path,[42])
