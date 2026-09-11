"""Consume actual standard CAFE output without a custom optimizer or CAFE patch."""
import csv
import hashlib
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]


def table(path, fields, rows):
    with path.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(fields)
        writer.writerows(rows)


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


@pytest.mark.parametrize("model", ["Base", "Gamma"])
def test_standard_cafe_outputs_support_gain_loss_flags_and_replay(tmp_path, model):
    if not shutil.which('cafe5') or not shutil.which('Rscript'):
        pytest.skip('GeneGalleon CAFE/R runtime is required')
    (tmp_path / 'tree.nwk').write_text('((A:1,B:1):1,(C:1,D:1):1);\n')
    table(tmp_path / 'families.tsv', ['Desc', 'FamilyID', 'A', 'B', 'C', 'D'], [
        ['(null)', 'gain', 20, 2, 2, 2], ['(null)', 'loss', 1, 20, 20, 20],
        ['(null)', 'stable_small', 2, 2, 2, 2], ['(null)', 'stable_large', 20, 20, 20, 20]])
    native = tmp_path / 'native'
    # Fixed lambda is an ordinary CAFE option used to generate this fixture.
    # The GO analysis below consumes saved outputs and never calls CAFE itself.
    command = ['cafe5', '-i', str(tmp_path / 'families.tsv'), '-t', str(tmp_path / 'tree.nwk'),
               '-l', '0.02', '-c', '1', '-o', str(native)]
    if model == 'Gamma':
        command.extend(['-k', '2', '-a', '1'])
    result = subprocess.run(command, text=True, capture_output=True, timeout=300)
    (tmp_path / 'cafe.log').write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    changes = read(native / f'{model}_change.tab')
    target = next(k for k in changes[0] if k.startswith('A<'))
    table(tmp_path / 'ids.tsv', ['Orthogroup', 'ref'], [
        ['gain', 'g1, g1b'], ['loss', 'g2'], ['stable_small', 'g3'], ['stable_large', 'g4']])
    table(tmp_path / 'ref.annotation.tsv', ['gene_id', 'go_ids', 'go_aspects', 'go_terms'], [
        ['g1', 'GO:G', 'BP', 'gain'], ['g1b', 'GO:G', 'BP', 'gain'],
        ['g2', 'GO:L', 'BP', 'loss'], ['g3', 'GO:U', 'BP', 'other'], ['g4', 'GO:U', 'BP', 'other']])
    hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in native.iterdir() if p.is_file()}
    output = tmp_path / 'go'
    command = ['Rscript', str(ROOT / 'workflow/support/cafe_go_enrichment.r'),
        str(native / f'{model}_change.tab'), str(native / f'{model}_branch_probabilities.tab'),
        str(tmp_path / 'ids.tsv'), str(tmp_path / 'ref.annotation.tsv'), str(output),
        target, 'both', 'BP', 'cafe_branch_flags', '.05']
    result = subprocess.run(command, text=True, capture_output=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    families = read(output / 'family_branch_flags.tsv')
    expected_selected = {'gain', 'loss'} if model == 'Base' else {'gain'}
    assert {r['FamilyID'] for r in families if r['selected'] == 'TRUE'} == expected_selected
    assert len(families) == 4
    assert {r['direction'] for r in families if r['selected'] == 'TRUE'} == (
        {'increase', 'decrease'} if model == 'Base' else {'increase'})
    if model == 'Gamma':
        # The mixture reconstruction also flags other branches for this loss.
        # It must be rejected even though the target loss is strongly flagged.
        loss = next(r for r in families if r['FamilyID'] == 'loss')
        assert loss['selection_status'] == 'other_branches_flagged'
        assert int(loss['n_other_flagged']) > 0
    go = read(output / f'enrichment_significant_both_{target}_all_go.tsv')
    assert {(r['go_ids'], r['direction']) for r in go} == {('GO:G', 'increase'), ('GO:L', 'decrease')}
    assert next(r for r in go if r['go_ids'] == 'GO:G')['n_selected_in_go'] == '1'
    loss_go = next(r for r in go if r['go_ids'] == 'GO:L')
    assert loss_go['n_selected_in_go'] == ('1' if model == 'Base' else '0')
    if model == 'Gamma':
        assert float(loss_go['p_value']) == 1.0
    metadata = read(output / 'branch_flags_metadata.tsv')[0]
    assert metadata['interpretation'] == 'exploratory_target_restricted_native_flags_not_a_rate_contrast'
    saved = {p: p.read_bytes() for p in output.iterdir() if p.is_file()}
    replay = subprocess.run(command, text=True, capture_output=True, timeout=60)
    assert replay.returncode == 0, replay.stdout + replay.stderr
    assert saved == {p: p.read_bytes() for p in output.iterdir() if p.is_file()}
    assert hashes == {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in hashes}

    # An identically shaped report from a different output prefix must not be
    # silently combined with these changes. Failed reruns invalidate summaries.
    mismatched = tmp_path / 'Other_branch_probabilities.tab'
    shutil.copyfile(native / f'{model}_branch_probabilities.tab', mismatched)
    command[3] = str(mismatched)
    rejected = subprocess.run(command, text=True, capture_output=True, timeout=60)
    assert rejected.returncode != 0
    assert 'matching native *_branch_probabilities.tab' in rejected.stderr
    assert not (output / 'family_branch_flags.tsv').exists()
    assert not (output / 'branch_flags_metadata.tsv').exists()
    assert not (output / f'enrichment_significant_both_{target}_all_go.tsv').exists()
    assert not (output / f'enrichment_significant_both_{target}_significant_go.tsv').exists()
