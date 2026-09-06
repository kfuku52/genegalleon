"""Exercise the stdin execution bridge and its real telemetry subprocess."""
import json
import shlex
import subprocess
from pathlib import Path


def test_container_profile_override_does_not_leak_between_runs():
    util = Path(__file__).resolve().parents[1] / "support/gg_util.sh"
    code = f'''source {shlex.quote(str(util))}
export GG_RESOURCE_PROFILE=/workspace/input/current.json
export SINGULARITYENV_GG_RESOURCE_PROFILE=/old/project/profile.json
export APPTAINERENV_GG_RESOURCE_PROFILE=/old/project/profile.json
unset_singularity_envs
[[ -z "${{SINGULARITYENV_GG_RESOURCE_PROFILE:-}}" && -z "${{APPTAINERENV_GG_RESOURCE_PROFILE:-}}" ]]
[[ "$GG_RESOURCE_PROFILE" == /workspace/input/current.json ]]
'''
    assert subprocess.run(["bash", "-e", "-c", code]).returncode == 0


def test_bridge_preserves_failure_output_and_records_stage_boundaries(tmp_path):
    workflow = Path(__file__).resolve().parents[1]
    util = workflow / "support/gg_util.sh"
    helper = workflow / "support/resource_metrics.py"
    output = tmp_path / "metrics"
    core = tmp_path / "gg_gene_evolution_core.sh"
    core.write_text(f'''source {shlex.quote(str(util))}
gg_support_dir={shlex.quote(str(workflow / 'support'))}
export GG_RESOURCE_OWNER_PID=$BASHPID
gg_step_start first
python -c "sum(i*i for i in range(1000000))"
gg_step_start second
echo scientific-output
exit 7
''')
    # A local adapter maps the same /script and /workspace mounts as the real
    # container adapter. Python and Bash execution and signal/status handling
    # are real; scheduler calls are absent.
    adapter = tmp_path / "adapter.py"
    adapter.write_text('''import os,sys
args=sys.argv[3:]
args=[sys.argv[1] if x=="/script/support/resource_metrics.py" else sys.argv[2] if x=="/workspace/output/resource_metrics" else x for x in args]
os.execvp(args[0],args)
''')
    code = f'''source {shlex.quote(str(util))}
container_adapter() {{ shift 2; python {shlex.quote(str(adapter))} {shlex.quote(str(helper))} {shlex.quote(str(output))} "$@"; }}
singularity_command=(container_adapter exec)
GG_TASK_CPUS=4
GG_MEM_TOTAL_GB=8
gg_run_container_shell_script image.sif {shlex.quote(str(core))}
'''
    result = subprocess.run(["bash", "-c", code], text=True, capture_output=True)
    assert result.returncode == 7, result.stderr
    assert "scientific-output" in result.stdout
    records = list(output.glob("attempt-*.json"))
    assert len(records) == 1
    record = json.loads(records[0].read_text())
    assert record["cpus"] == 4
    assert record["exit_code"] == 7
    assert [stage["stage"] for stage in record["stage_boundaries"]] == ["first", "second"]
    assert record["stage_boundaries"][0]["cpu_seconds_reaped"] > 0
