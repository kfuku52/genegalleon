import re

from shell_static_helpers import (
    REPO_ROOT,
    WORKFLOW_DIR,
    entrypoint_scheduler_header,
    function_body,
    read_text,
)


def test_set_singularity_command_supports_apptainer_fallback():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = read_text(util_path)
    assert "gg_detect_container_runtime_binary()" in text
    assert "command -v singularity" in text
    assert "command -v apptainer" in text
    assert "Neither singularity nor apptainer was found on PATH." in text


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
    assert "s|@@GG_VERSION@@|" in local_build
    assert "org.opencontainers.image.version @@GG_VERSION@@" in definition_template


def test_container_defaults_install_program_sources_from_moving_branches():
    dockerfile = read_text(REPO_ROOT / "container" / "Dockerfile")
    buildx = read_text(REPO_ROOT / "container" / "buildx.sh")
    local_build = read_text(REPO_ROOT / "container" / "apptainer_local_build.sh")
    ensure_latest = read_text(REPO_ROOT / "container" / "scripts" / "ensure_kfuku52_tools_latest.sh")
    capability_matrix = read_text(REPO_ROOT / "container" / "CAPABILITY_MATRIX.md")

    assert 'ARG KFU52_AMALGKIT_AUTO_SELECT_REF="0"' in dockerfile
    source_branches = read_text(REPO_ROOT / "container" / "source_branches.env")
    assert 'ARG KFU52_AMALGKIT_REPO_REF=""' in dockerfile
    assert "KFU52_AMALGKIT_AUTO_SELECT_REF=${KFU52_AMALGKIT_AUTO_SELECT_REF:-0}" in buildx
    assert "KFU52_AMALGKIT_REPO_REF=${KFU52_AMALGKIT_REPO_REF:-${GG_SOURCE_AMALGKIT_REPO_REF}}" in buildx
    assert "KFU52_AMALGKIT_AUTO_SELECT_REF=${KFU52_AMALGKIT_AUTO_SELECT_REF:-0}" in local_build
    assert "KFU52_AMALGKIT_REPO_REF=${KFU52_AMALGKIT_REPO_REF:-${GG_SOURCE_AMALGKIT_REPO_REF}}" in local_build
    assert "amalgkit_auto_select_ref=${KFU52_AMALGKIT_AUTO_SELECT_REF:-0}" in ensure_latest
    assert "amalgkit_repo_ref_override=${KFU52_AMALGKIT_REPO_REF-master}" in ensure_latest
    assert "Git-sourced programs follow the moving branches" in capability_matrix
    assert "auto-selects newer commit among `master`, `kfdevel`, and `devel`" not in capability_matrix
    assert 'ARG KFU52_CSUBST_REPO_REF=""' in dockerfile
    assert "KFU52_CSUBST_REPO_REF=${KFU52_CSUBST_REPO_REF:-${GG_SOURCE_CSUBST_REPO_REF}}" in buildx
    assert "KFU52_CSUBST_REPO_SHA=${KFU52_CSUBST_REPO_SHA:-}" in buildx
    assert "KFU52_CSUBST_REPO_REF=${KFU52_CSUBST_REPO_REF:-${GG_SOURCE_CSUBST_REPO_REF}}" in local_build
    assert "KFU52_CSUBST_REPO_SHA=${KFU52_CSUBST_REPO_SHA:-}" in local_build
    assert "csubst_repo_ref=${KFU52_CSUBST_REPO_REF:-master}" in ensure_latest
    assert "GG_SOURCE_NWKIT_REPO_REF=master" in source_branches
    sha_overrides = (
        "KFU52_AMALGKIT_REPO_SHA",
        "KFU52_CDSKIT_REPO_SHA",
        "KFU52_CSUBST_REPO_SHA",
        "KFU52_NWKIT_REPO_SHA",
        "BUSCO_REPO_SHA",
        "PAML_REPO_SHA",
        "KFL1OU_REPO_SHA",
        "KFFRACTBIAS_REPO_SHA",
        "KFTOOLS_REPO_SHA",
        "RKFTOOLS_REPO_SHA",
        "RADTE_REPO_SHA",
    )
    for sha_var in sha_overrides:
        expected = f"{sha_var}=${{{sha_var}:-}}"
        assert expected in buildx
        assert expected in local_build


def test_run_container_shell_script_uses_exec_with_bash_stdin_bridge():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = read_text(util_path)
    body = function_body(text, "gg_run_container_shell_script")
    assert "subcommand=$(gg_container_shell_command_subcommand || true)" in body
    assert '"${singularity_command[@]}" "${image_path}" bash -s -- < "${script_path}"' in body
    assert '${singularity_command} "${image_path}" bash -s -- < "${script_path}"' in body
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


def test_gg_util_has_common_forward_config_export_helpers():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = read_text(util_path)
    assert "gg_export_var_to_container_env_if_set()" in text
    assert "gg_apply_named_env_overrides()" in text
    assert "gg_apply_registered_env_overrides()" in text
    assert "gg_entrypoint_env_override_prefix()" in text
    assert "gg_forward_env_vars_with_prefix_to_container_env()" in text
    assert "gg_print_named_config_summary()" in text
    assert "gg_print_registered_config_summary()" in text
    assert "gg_require_versions_dump()" in text
    assert "gg_prepare_entrypoint_runtime_snapshot()" in text
    assert "gg_entrypoint_runtime_snapshot_dir()" in text
    assert "gg_resolve_physical_path()" in text
    body = function_body(text, "forward_config_vars_to_container_env")
    assert "gg_print_entrypoint_config_vars" in body
    assert 'gg_export_var_to_container_env_if_set "gg_debug_mode"' in body
    body = function_body(text, "gg_apply_registered_env_overrides")
    assert "gg_print_entrypoint_config_vars" in body
    assert "gg_entrypoint_env_override_prefix" in body


def test_entrypoint_config_var_registry_covers_all_entrypoints():
    registry_text = read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")
    entrypoints = sorted(path.name for path in WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    assert entrypoints
    for entrypoint_name in entrypoints:
        assert f"{entrypoint_name})" in registry_text


def test_entrypoints_define_scheduler_sections_in_canonical_order():
    entrypoints = sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    assert entrypoints, "No entrypoint scripts were found."
    section_tokens = [
        "# SLURM",
        "## UGE",
        "## PBS",
        '# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):',
    ]
    for script in entrypoints:
        header = entrypoint_scheduler_header(script)
        positions = []
        for token in section_tokens:
            assert token in header, f"Missing scheduler header section {token!r} in {script}"
            positions.append(header.index(token))
        assert positions == sorted(positions), f"Scheduler header section order drifted in {script}"


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


def test_gg_util_uses_explicit_conda_shell_initialization():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = read_text(util_path)
    body = function_body(text, "gg_initialize_conda_shell")
    assert "/home/.bashrc" not in text
    assert "micromamba shell hook --shell bash" in body
    assert "conda() {" in text
    assert 'micromamba "$@"' in text
    assert "source /opt/conda/etc/profile.d/conda.sh" in text


def test_gg_add_container_bind_mount_skips_duplicate_destinations():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = read_text(util_path)
    body = function_body(text, "gg_add_container_bind_mount")
    assert "gg_container_mount_destination()" in text
    assert "gg_container_bind_destination_exists()" in text
    assert 'if gg_container_bind_destination_exists "${mount_spec}"; then' in body
    assert 'GG_CONTAINER_BIND_MOUNTS=$(gg_csv_prepend "${mount_spec}" "${GG_CONTAINER_BIND_MOUNTS:-}")' in body
    assert "gg_sync_container_bind_envs" in body
