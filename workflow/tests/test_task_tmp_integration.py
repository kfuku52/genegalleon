"""Full workflow checks with isolated external scratch and deterministic tools."""
import sys
from pathlib import Path

SUPPORT = Path(__file__).resolve().parents[1] / 'support'
RUNNER = SUPPORT / 'task_tmp.py'


def test_input_array_workers_with_separate_scratch_roots(tmp_path, monkeypatch):
    # Exercise the full prepare/parallel workers/finalize core flow with the
    # existing deterministic toolchain fixture. Each job has isolated scratch,
    # as on separate compute nodes; shared plans and shards must still work.
    import shlex
    import test_gg_input_generation_end_to_end as integration

    original_core = integration.CORE_PATH
    original_env = integration._core_env
    wrapper = tmp_path / 'external-core.sh'
    wrapper.write_text(
        'exec ' + shlex.join([sys.executable, str(RUNNER), '--workflow',
                             'gg_input_generation', '--', 'bash', str(original_core)]) + '\n'
    )

    def external_env(*args, **kwargs):
        env = original_env(*args, **kwargs)
        scratch = tmp_path / ('scratch-' + env['input_generation_mode'] + '-' + env.get('GG_ARRAY_TASK_ID', '1'))
        scratch.mkdir(exist_ok=True)
        env.update(GG_COMMON_TMP_ROOT=str(scratch), GG_TMP_MOUNT=str(scratch),
                   GG_TMP_WORKSPACE_ID=env['gg_workspace_dir'])
        return env

    monkeypatch.setattr(integration, 'CORE_PATH', wrapper)
    monkeypatch.setattr(integration, '_core_env', external_env)
    integration.test_gg_input_generation_array_mode_end_to_end_with_parallel_workers(tmp_path)
    assert not list(tmp_path.glob('scratch-*/genegalleon-*/*/*/run-*'))



def test_genome_protein_staging_uses_selected_scratch(tmp_path, monkeypatch):
    import test_genome_evolution_protein_mode as integration

    scratch = tmp_path / 'scratch'
    scratch.mkdir()
    monkeypatch.setenv('GG_COMMON_TMP_ROOT', str(scratch))
    monkeypatch.setenv('GG_TMP_TASK_ROOT', str(scratch / 'work'))
    integration.test_genome_evolution_protein_mode_prefers_species_protein_inputs(tmp_path)
    inputs = (tmp_path / 'capture/input_files.txt').read_text().splitlines()
    assert inputs
    assert all(Path(p).is_relative_to(scratch) for p in inputs)
