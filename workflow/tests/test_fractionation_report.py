import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'workflow/support'))
from fractionation_report import paired_contrasts, run  # noqa: E402


def loci():
    return [dict(analysis_id=a, target_seqid='chr1', target_gene=f'g{i}', target_rank=str(i + 1),
                 common_callable='1', retained_copies=str(c), expected_copies=str(k))
            for a, c, k in [('a', 2, 2), ('b', 2, 5)] for i in range(12)]


def test_paired_contrasts_use_disjoint_reference_blocks_and_depth():
    rows = paired_contrasts(loci(), block_size=3, replicates=100, permutations=100)
    r = rows[0]
    assert r['effect'] == pytest.approx(.6)
    assert r['n_loci'] == 12 and r['n_blocks'] == 4
    assert r['ci_low'] == pytest.approx(.6) and r['ci_high'] == pytest.approx(.6)
    assert r['p_value'] == .125 and r['q_value'] == .125


def test_shared_missingness_and_schema():
    rows = loci()
    for r in rows:
        if int(r['target_rank']) <= 3:
            r['common_callable'] = '0'
    assert paired_contrasts(rows, block_size=3, replicates=100, permutations=100)[0]['n_loci'] == 9
    rows[-1]['common_callable'] = '0'
    with pytest.raises(ValueError, match='mask'):
        paired_contrasts(rows, replicates=100, permutations=100)


def test_report_config_rejects_unknown_keys_before_outputs(tmp_path):
    config = tmp_path / 'report.json'
    config.write_text(json.dumps({'manifest': 'missing.tsv', 'unknown': True}))
    with pytest.raises(ValueError, match='supported'):
        run(config, tmp_path / 'out')
    assert not (tmp_path / 'out').exists()


def test_report_core_calls_native_dependency_without_pair_table(tmp_path):
    # Reuse the owning dependency's public synthetic report fixture.
    import kffractbias.comparison  # noqa: F401
    reference = tmp_path / 'ref.bed'
    query = tmp_path / 'query.bed'
    anchors = tmp_path / 'pairs.anchors'
    mapping = tmp_path / 'mapping.tsv'
    reference.write_text(''.join(f'chr1\t{i*100}\t{i*100+90}\tt{i}\n' for i in range(12)))
    query.write_text(''.join(f'c{c}\t{i*100}\t{i*100+90}\tq{c}_{i}\n' for c in 'AB' for i in range(12)))
    anchors.write_text(''.join('###\n' + ''.join(f't{i}\tq{c}_{i}\t10\n' for i in range(12)) for c in 'AB'))
    mapping.write_text('gene_id\tgroup_id\tcopy\tsegment_id\n' + ''.join(f'q{c}_{i}\tG1\t{c}\tc{c}\n' for c in 'AB' for i in range(12)))
    manifest = tmp_path / 'manifest.tsv'
    manifest.write_text('analysis_id\tspecies\texpected_copies\ttarget_bed\tquery_bed\tsynteny_file\tmapping_file\n' + ''.join(f'{a}\t{a}_species\t2\t{reference}\t{query}\t{anchors}\t{mapping}\n' for a in ['a', 'b']))
    config = tmp_path / 'config.json'
    config.write_text(json.dumps(dict(manifest=str(manifest), window_size=3, minimum_anchors=2, trim=0,
                                     statistics=dict(block_size=3, replicates=100, permutations=100))))
    workspace = tmp_path / 'workspace'
    env = dict(os.environ, GG_COMMON_TMP_ROOT='workspace', gg_workspace_dir=str(workspace),
               GG_ARRAY_TASK_ID='1', GG_JOB_ID='test', GG_TASK_CPUS='1', GG_MEM_PER_CPU_GB='8',
               kffractbias_report_config=str(config))
    p = subprocess.run(['bash', str(ROOT / 'workflow/core/gg_fractionation_bias_core.sh')], env=env,
                       capture_output=True, text=True, timeout=90)
    assert p.returncode == 0, p.stdout + p.stderr
    result = workspace / 'output/kffractbias_reports'
    summary = json.loads((result / 'comparison.report.json').read_text())
    assert summary['common_callable_loci'] == 12
    contrasts = json.loads((result / 'comparison.contrasts.json').read_text())
    assert contrasts['contrasts'][0]['effect'] == 0
