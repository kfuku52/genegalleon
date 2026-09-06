"""Resource contracts across planning, runtime, and optional telemetry."""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'support'))
from resource_metrics import run  # noqa: E402 - support directory is added above


@pytest.mark.parametrize('slots,per_slot,expected', [(32, '1056M', 33), (32, '128M', 4), (3, '1366M', 4)])
def test_gridengine_exact_total_budget(slots, per_slot, expected):
    code = f'''source "{ROOT}/support/gg_util.sh"
qstat() {{ echo 'hard_resource_list: s_vmem={per_slot}'; }}
NSLOTS={slots}
JOB_ID=123
gg_normalize_scheduler_env >/dev/null
echo "$GG_MEM_TOTAL_GB"
'''
    env = {k: v for k, v in os.environ.items() if not k.startswith(('GG_', 'MEM_', 'SLURM_', 'PBS_'))}
    result = subprocess.run(['bash', '-c', code], env=env, text=True, capture_output=True, check=True)
    assert int(result.stdout.strip()) == expected


def test_notung_workers_share_the_memory_divisor():
    core = (ROOT / 'core/gg_genome_evolution_core.sh').read_text()
    for line in core.splitlines():
        if 'wait_until_jobn_le' in line:
            assert 'GG_TASK_CPUS' not in line
    start = core.index('GG_GENOME_PARALLEL_JOBS=${GG_TASK_CPUS}')
    block = core[start:core.index('iqtree_full_mem_args=', start)]
    result = subprocess.run(['bash', '-c', f'''source "{ROOT}/support/gg_util.sh"
GG_TASK_CPUS=32; GG_MEM_TOOL_GB=28; genome_parallel_jobs=2; genome_parallel_memory_gb_per_job=2
{block}
echo "budget=$((GG_GENOME_PARALLEL_JOBS * memory_notung))"
'''], text=True, capture_output=True, check=True)
    assert 'budget=28' in result.stdout


@pytest.mark.parametrize('profile', [{}, [], {'workflow': 'gg_genome_evolution'}, None])
def test_optional_bad_profile_preserves_command_status(tmp_path, profile):
    path = tmp_path / 'profile.json'
    path.write_text(json.dumps(profile))
    result = subprocess.run([sys.executable, str(ROOT / 'support/resource_metrics.py'),
        '--directory', str(tmp_path / 'metrics'), '--workflow', 'gg_genome_evolution',
        '--runtime-id', 'test', '--server-id', 'test', '--', sys.executable, '-c',
        "print('science-started'); raise SystemExit(7)"],
        env={**os.environ, 'GG_RESOURCE_PROFILE': str(path)}, text=True, capture_output=True)
    assert result.returncode == 7
    assert 'science-started' in result.stdout


@pytest.mark.parametrize('event', ['{}', '[]', '{"stage":"x","kind":"start","at":"bad","cpu":0}'])
def test_bad_event_does_not_fail_successful_command(tmp_path, event):
    command = [sys.executable, '-c', f"import os; open(os.environ['GG_RESOURCE_EVENT_FILE'],'w').write({event!r}+'\\n')"]
    assert run(command, tmp_path, workflow='gg_genome_evolution', input_sha256='a'*64,
               cpus=4, memory_gb=32, runtime_id='test', server_id='test', expected_stages=['x']) == 0
    record = json.loads(next(tmp_path.glob('attempt-*.json')).read_text())
    assert record['learning_eligible'] is False
    assert not list(tmp_path.glob('history/*/*.json'))
