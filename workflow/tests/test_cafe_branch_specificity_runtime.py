"""Small real-CAFE integration, including the R GO consumer and cache replay."""
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[2]


def table(path, fields, rows):
    with path.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(fields)
        writer.writerows(rows)


def test_native_internal_target_rate_mapping(tmp_path):
    if not shutil.which('cafe5'):
        pytest.skip('GeneGalleon CAFE runtime is required')
    from workflow.support.cafe_branch_specificity import NativeCafe, read_cafe_tree, lr_statistic
    asr = tmp_path / 'Gamma_asr.tre'
    asr.write_text('''#nexus
BEGIN TREES;
 TREE family = ((A<1>_9:1,B<2>_9:1)<5>_9:1,C<3>_3:2)<4>_3;
END;
''')
    engine = NativeCafe(read_cafe_tree(asr, '<5>'), tmp_path / 'native',
                        'cafe5', 5, 1000, 1, 180)
    # Shared high counts in the sister tips should fit a fast internal branch,
    # verifying native lambda ordering for an ancestral (not terminal) target.
    null = engine.fit([9, 9, 3], False)
    alternative = engine.fit([9, 9, 3], True)
    assert lr_statistic(null, alternative) > 0
    assert alternative.lambdas[1] > alternative.lambdas[0]


def test_native_cafe_lrt_go_both_directions_and_resume(tmp_path):
    if not shutil.which('cafe5') or not shutil.which('Rscript'):
        pytest.skip('GeneGalleon CAFE/R runtime is required')
    prefix = tmp_path / 'Gamma'
    (tmp_path / 'Gamma_asr.tre').write_text('''#nexus
BEGIN TREES;
 TREE gain = ((A<1>_5:1,B<2>_3:1)<5>_4:1,C<3>_4:2)<4>_4;
END;
''')
    labels = ['A<1>', 'B<2>', 'C<3>', '<4>', '<5>']
    table(tmp_path / 'Gamma_count.tab', ['FamilyID', *labels],
          [['gain', 5, 3, 4, 4, 4], ['loss', 3, 5, 4, 4, 4]])
    table(tmp_path / 'Gamma_change.tab', ['FamilyID', *labels],
          [['gain', 1, -1, 0, 0, 0], ['loss', -1, 1, 0, 0, 0]])
    table(tmp_path / 'Gamma_branch_probabilities.tab', ['#FamilyID', *labels],
          [['gain', .01, .01, .9, 'N/A', .9], ['loss', .01, .01, .9, 'N/A', .9]])
    table(tmp_path / 'ids.tsv', ['Orthogroup', 'ref'], [['gain', 'g1'], ['loss', 'g2']])
    table(tmp_path / 'ref.annotation.tsv', ['gene_id', 'go_ids', 'go_aspects', 'go_terms'],
          [['g1', 'GO:G', 'BP', 'gain'], ['g2', 'GO:L', 'BP', 'loss']])
    output = tmp_path / 'output'
    command = ['Rscript', str(ROOT / 'workflow/support/cafe_go_enrichment.r'),
               str(prefix) + '_change.tab', str(prefix) + '_branch_probabilities.tab',
               str(tmp_path / 'ids.tsv'), str(tmp_path / 'ref.annotation.tsv'),
               str(output), 'A<1>', 'both', 'BP', 'cafe_lrt', '.05', '1', '8', '1000', '1']
    env = {**os.environ, 'PYTHONDONTWRITEBYTECODE': '1'}
    # Up to 64 native optimizations plus native reconstruction/simulation.
    # The wall-clock budget must also cover the converged boundary fits.
    first = subprocess.run(command, text=True, capture_output=True, timeout=1800, env=env)
    assert first.returncode == 0, first.stdout + first.stderr
    with (output / 'family_specificity.tsv').open() as handle:
        families = list(csv.DictReader(handle, delimiter='\t'))
    assert {r['direction'] for r in families} == {'increase', 'decrease'}
    assert all(r['status'] == 'tested' for r in families)
    # With one bootstrap draw, P is >=1/2; never pretend this smoke has power.
    assert all(float(r['p_value']) >= .5 and r['selected'] == 'FALSE' for r in families)
    with (output / 'enrichment_significant_both_A<1>_all_go.tsv').open() as handle:
        go = list(csv.DictReader(handle, delimiter='\t'))
    assert {(r['go_ids'], r['direction']) for r in go} == {('GO:G', 'increase'), ('GO:L', 'decrease')}
    assert all(float(r['p_value_adjusted']) == 1 for r in go)
    metadata = json.loads((output / 'native_cafe/metadata.json').read_text())
    assert metadata['status'] == 'complete'
    assert metadata['null_model'] == 'one_lambda_per_family'
    logs = {p: p.stat().st_mtime_ns for p in (output / 'native_cafe/runs').glob('*/*/attempt_*/cafe.log')}
    assert logs
    saved_results = (output / 'family_specificity.tsv').read_text()
    replay = subprocess.run(command, text=True, capture_output=True, timeout=60, env=env)
    assert replay.returncode == 0, replay.stdout + replay.stderr
    assert saved_results == (output / 'family_specificity.tsv').read_text()
    assert logs == {p: p.stat().st_mtime_ns for p in (output / 'native_cafe/runs').glob('*/*/attempt_*/cafe.log')}
