import re
import subprocess

from shell_static_helpers import (
    REPO_ROOT,
    WORKFLOW_DIR,
    entrypoint_scheduler_header,
    function_body,
    read_text,
)


def test_common_scratch_defaults_to_node_local_tmp():
    common_params = read_text(WORKFLOW_DIR / "gg_common_params.sh")
    assert ': "${GG_COMMON_TMP_ROOT:=/tmp}"' in common_params
    assert "GG_COMMON_TMP_ROOT=workspace" in read_text(
        REPO_ROOT / "docs" / "temporary-storage.md"
    )
    assert "default is\nnode-local `/tmp`" in read_text(
        REPO_ROOT / "docs" / "execution-model.md"
    )


def test_site_runtime_exec_command_uses_output_parameter():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    site_runtime_path = WORKFLOW_DIR / "support" / "gg_site_runtime.sh"
    util_text = read_text(util_path)
    site_text = read_text(site_runtime_path)
    site_body = function_body(site_text, "gg_site_container_shell_command")
    helper_body = function_body(site_text, "gg_set_command_array")
    assert "local out_var=${2:-}" in site_body
    assert 'gg_set_command_array "${out_var}" "${runtime_bin}" exec || return 1' in site_body
    assert 'gg_set_command_array "${out_var}" "${runtime_bin}" exec --contain || return 1' in site_body
    assert 'printf -v "${out_var}" \'%s\' "${command_text}"' not in site_body
    assert 'echo "${runtime_bin} shell"' not in site_body
    assert 'echo "${runtime_bin} shell --contain"' not in site_body
    assert 'echo "${runtime_bin} exec"' not in site_body
    assert 'echo "${runtime_bin} exec --contain"' not in site_body
    assert 'site profile = nig" >&2' not in site_body
    assert 'site profile = nhr-fau" >&2' not in site_body
    assert 'site profile = default" >&2' not in site_body
    assert 'eval "${out_var}=()"' in helper_body
    assert 'eval "${out_var}+=( ${quoted_arg} )"' in helper_body
    assert 'gg_site_container_shell_command "${runtime_bin}" singularity_command' in util_text
    assert 'singularity_command="$(gg_site_container_shell_command "${runtime_bin}")"' not in util_text
    assert 'singularity_command=( "${runtime_bin}" exec )' in util_text


def test_set_singularityenv_does_not_dump_singularityenv_values():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = read_text(util_path)
    assert 'set | grep "^SINGULARITY"' not in text
    assert "forwarded_container_env_vars" in text


def test_all_entrypoints_call_set_singularity_command():
    entrypoints = sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    assert entrypoints, "No entrypoint scripts were found."
    for script in entrypoints:
        text = read_text(script)
        assert "set_singularity_command" in text or "gg_entrypoint_prepare_container_runtime" in text, (
            f"Missing container runtime preparation call: {script}"
        )


def test_gg_trigger_versions_dump_is_runtime_agnostic():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = read_text(util_path)
    body = function_body(text, "gg_trigger_versions_dump")
    assert 'export SINGULARITYENV_GG_VERSION="${gg_version}"' in body
    assert 'export APPTAINERENV_GG_VERSION="${gg_version}"' in body
    assert 'command -v "${container_runtime_bin}"' in body
    assert 'container_runtime_bin="$(gg_container_shell_command_runtime_bin || true)"' in body
    assert 'gg_run_container_shell_script "${gg_container_image_path}" "${versions_script}"' in body
    assert '"${container_runtime_bin}" inspect "${gg_container_image_path}"' in body
    assert '"${container_runtime_bin}" version || {' in body
    assert 'singularity inspect "${gg_container_image_path}"' not in body
    assert "singularity version || {" not in body


def test_entrypoint_activate_container_runtime_prints_version_summary():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = read_text(util_path)
    body = function_body(text, "gg_entrypoint_activate_container_runtime")
    assert "gg_entrypoint_print_version_summary" in body


def test_container_build_metadata_includes_repo_version_label():
    dockerfile = read_text(REPO_ROOT / "container" / "Dockerfile")
    buildx = read_text(REPO_ROOT / "container" / "buildx.sh")
    local_build = read_text(REPO_ROOT / "container" / "apptainer_local_build.sh")
    definition_template = read_text(REPO_ROOT / "container" / "apptainer_local_build.def.template")

    assert 'org.opencontainers.image.version="${GG_VERSION}"' in dockerfile
    assert '--build-arg GG_VERSION="${gg_version}"' in buildx
    assert '--build-arg SECURITY_REFRESH_EPOCH="${SECURITY_REFRESH_EPOCH}"' in buildx
    assert "s|@@GG_VERSION@@|" in local_build
    assert "s|@@SECURITY_REFRESH_EPOCH@@|" in local_build
    assert "org.opencontainers.image.version @@GG_VERSION@@" in definition_template
    assert "io.genegalleon.security-refresh-epoch @@SECURITY_REFRESH_EPOCH@@" in definition_template


def test_run_container_shell_script_uses_exec_with_bash_stdin_bridge():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = read_text(util_path)
    body = function_body(text, "gg_run_container_shell_script")
    assert "subcommand=$(gg_container_shell_command_subcommand || true)" in body
    assert '"${singularity_command[@]}" "${image_path}" "${shell_argv[@]}" < "${script_path}"' in body
    assert '${singularity_command} "${image_path}" "${shell_argv[@]}" < "${script_path}"' in body
    assert '"${singularity_command[@]}" "${image_path}" < "${script_path}"' in body


def test_entrypoints_stream_core_scripts_via_container_shell_helper():
    entrypoints = sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    assert entrypoints, "No entrypoint scripts were found."
    for script in entrypoints:
        text = read_text(script)
        if "gg_core_dir=" not in text:
            continue
        assert 'gg_run_container_shell_script "${gg_container_image_path}"' in text
        assert '${singularity_command} "${gg_container_image_path}" <' not in text


def test_progress_summary_entrypoint_uses_auto_forwarding_and_normalized_nslots():
    entrypoint = WORKFLOW_DIR / "gg_progress_summary_entrypoint.sh"
    text = read_text(entrypoint)

    assert "forward_config_vars_to_container_env()" not in text
    assert 'gg_entrypoint_name="gg_progress_summary_entrypoint.sh"' in text
    assert 'forward_config_vars_to_container_env "${gg_entrypoint_name}"' in text
    assert "unset -f forward_config_vars_to_container_env" not in text
    assert "for exported_name in mode_transcriptome_assembly ncpu_progress_summary; do" not in text
    assert 'ncpu_progress_summary="${ncpu_progress_summary:-${GG_TASK_CPUS:-1}}"' not in text

    idx_variable_sgenizer = text.index("gg_entrypoint_prepare_container_runtime")
    idx_ncpu_default = text.index(': "${ncpu_progress_summary:=${GG_TASK_CPUS:-1}}"')
    assert idx_ncpu_default > idx_variable_sgenizer


def test_input_generation_entrypoint_forwards_env_driven_overrides():
    entrypoint = WORKFLOW_DIR / "gg_input_generation_entrypoint.sh"
    text = read_text(entrypoint)

    assert 'gg_apply_registered_env_overrides "${gg_entrypoint_name}"' in text
    assert "gg_apply_named_env_overrides \\" not in text
    assert "for gg_input_var_name in ${!GG_INPUT_@}; do" not in text
    assert 'gg_forward_env_vars_with_prefix_to_container_env "GG_INPUT_MAX_CONCURRENT_DOWNLOADS_"' in text
    assert 'gg_forward_env_vars_with_prefix_to_container_env "GG_DOWNLOAD_"' in text
    assert 'GG_INPUT_MAX_CONCURRENT_DOWNLOADS_CNGB:=1' in text
    assert 'export "SINGULARITYENV_${gg_input_var_name}=${!gg_input_var_name}"' not in text
    assert 'export "APPTAINERENV_${gg_input_var_name}=${!gg_input_var_name}"' not in text


def test_input_generation_entrypoint_is_array_ready():
    entrypoint = WORKFLOW_DIR / "gg_input_generation_entrypoint.sh"
    text = read_text(entrypoint)

    assert "#SBATCH -a 1" in text
    assert "#$ -t 1" in text
    assert "#PBS -J 1" in text
    assert "gg_input_generation_entrypoint.sh_%A_%a.out" in text
    assert "gg_input_generation_entrypoint.sh_%A_%a.err" in text
    assert "task_plan.json" in text


def test_input_generation_trait_profile_preset_is_wired():
    entrypoint = WORKFLOW_DIR / "gg_input_generation_entrypoint.sh"
    core = WORKFLOW_DIR / "core" / "gg_input_generation_core.sh"
    entry_text = read_text(entrypoint)
    core_text = read_text(core)

    assert 'trait_profile="none"' in entry_text
    assert "run_cds_fx2tab=1" in entry_text
    assert "run_species_busco=1" in entry_text
    assert "run_multispecies_summary=1" in entry_text
    assert "gg_apply_registered_env_overrides" in entry_text
    assert "trait_profile" in read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")
    assert "GG_INPUT_" not in core_text
    assert "apply_env_override()" not in core_text
    assert 'case "${trait_profile}" in' in core_text
    assert "gift_starter" in core_text
    assert "gbif_distribution" in core_text


def test_registered_config_vars_are_consumed_by_core_or_shared_runtime():
    registry = WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh"
    shared_consumers = {
        "artifact_stale_policy": read_text(
            WORKFLOW_DIR / "support" / "gg_util" / "06_workspace_validation.sh"
        ),
    }
    for variable, consumer in shared_consumers.items():
        assert variable in consumer

    for entrypoint in sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh")):
        core = WORKFLOW_DIR / "core" / entrypoint.name.replace(
            "_entrypoint.sh",
            "_core.sh",
        )
        assert core.is_file(), f"Missing core script for {entrypoint.name}"
        completed = subprocess.run(
            [
                "bash",
                "-lc",
                'source "$1"; gg_print_entrypoint_config_vars "$2"',
                "bash",
                str(registry),
                entrypoint.name,
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        assert completed.returncode == 0, completed.stderr
        core_text = read_text(core)
        unused = [
            variable
            for variable in completed.stdout.splitlines()
            if variable and variable not in core_text and variable not in shared_consumers
        ]
        assert not unused, (
            f"Registered variables are not consumed by {core.name} or an explicit "
            f"shared runtime consumer: {', '.join(unused)}"
        )


def test_entrypoints_use_active_scheduler_directives_in_header_template():
    entrypoints = sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    assert entrypoints, "No entrypoint scripts were found."
    for script in entrypoints:
        header = entrypoint_scheduler_header(script)
        assert "##PBS" not in header, f"Use active #PBS directives in {script}"
        assert "##SBATCH -N" not in header, f"Drop legacy commented node-count example from {script}"
        assert "##SBATCH -n" not in header, f"Drop legacy commented task-count example from {script}"
        assert "#PBS -S /bin/bash" in header, f"Missing PBS shell directive in {script}"
        assert "#PBS -V" in header, f"Missing PBS environment export directive in {script}"
        assert "#SBATCH --ignore-pbs" in header, f"Missing Slurm PBS-ignore guard in {script}"


def test_entrypoints_use_shared_slurm_partition_fallbacks():
    entrypoints = sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    assert entrypoints, "No entrypoint scripts were found."
    for script in entrypoints:
        header = entrypoint_scheduler_header(script)
        assert "#SBATCH -p epyc,rome,medium" in header, f"Missing shared Slurm partitions in {script}"


def test_entrypoint_scheduler_directives_are_left_aligned():
    entrypoints = sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    assert entrypoints, "No entrypoint scripts were found."
    bad_lines = []
    for script in entrypoints:
        for lineno, line in enumerate(entrypoint_scheduler_header(script).splitlines(), start=1):
            if re.match(r"^ (?:#SBATCH|#PBS|#\$)", line):
                bad_lines.append(f"{script}:{lineno}: {line}")
    assert not bad_lines, "Left-align scheduler directives in entrypoint headers:\n" + "\n".join(bad_lines)
