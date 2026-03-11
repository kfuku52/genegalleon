from pathlib import Path
import re


REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_DIR = REPO_ROOT / "workflow"
CORE_DIR = WORKFLOW_DIR / "core"
CONTAINER_SCRIPTS_DIR = REPO_ROOT / "container" / "scripts"


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _function_body(text: str, function_name: str) -> str:
    pattern = re.compile(rf"^\s*{re.escape(function_name)}\(\)\s*\{{", re.MULTILINE)
    match = pattern.search(text)
    if match is None:
        raise AssertionError(f"Function not found: {function_name}")
    start = match.start()
    next_match = re.search(r"^\s*[A-Za-z_][A-Za-z0-9_]*\(\)\s*\{", text[match.end():], re.MULTILINE)
    if next_match is None:
        return text[start:]
    return text[start:match.end() + next_match.start()]


def _workflow_shell_scripts():
    return sorted(WORKFLOW_DIR.rglob("*.sh"))


def _container_shell_scripts():
    return sorted(CONTAINER_SCRIPTS_DIR.rglob("*.sh"))


def _core_and_entrypoint_scripts():
    core_scripts = sorted(CORE_DIR.glob("*.sh"))
    entrypoint_scripts = sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    return core_scripts + entrypoint_scripts


def _strict_mode_header(script: Path) -> str:
    max_lines = 60 if script.name.endswith("_entrypoint.sh") else 30
    return "\n".join(_read_text(script).splitlines()[:max_lines])


def _entrypoint_scheduler_header(script: Path) -> str:
    text = _read_text(script)
    marker = "set -euo pipefail"
    assert marker in text, f"Missing strict mode marker in {script}"
    return text.split(marker, 1)[0]


def _entrypoint_modify_block_assignments(script: Path):
    in_block = False
    for lineno, line in enumerate(_read_text(script).splitlines(), start=1):
        if "### Start: Modify this block to tailor your analysis ###" in line:
            in_block = True
            continue
        if "### End: Modify this block to tailor your analysis ###" in line:
            in_block = False
            continue
        if not in_block:
            continue
        stripped = line.strip()
        if re.match(r"^[A-Za-z_][A-Za-z0-9_]*=", stripped):
            yield lineno, stripped


def _common_param_assignments(script: Path):
    for lineno, line in enumerate(_read_text(script).splitlines(), start=1):
        stripped = line.strip()
        if re.match(r'^: "\$\{[A-Za-z0-9_]+:=', stripped):
            yield lineno, stripped


def _set_e_scripts():
    scripts = []
    for script in _workflow_shell_scripts():
        text = _read_text(script)
        header = _strict_mode_header(script)
        if ("set -eo pipefail" in header) or ("set -euo pipefail" in header):
            scripts.append((script, text))
    return scripts


def _is_within_double_quotes(line: str, index: int) -> bool:
    in_double = False
    escaped = False
    i = 0
    while i < index and i < len(line):
        ch = line[i]
        if escaped:
            escaped = False
            i += 1
            continue
        if ch == "\\":
            escaped = True
            i += 1
            continue
        if ch == '"':
            in_double = not in_double
        i += 1
    return in_double


def _unquoted_brace_expansions(line: str):
    # Remove simple command substitutions so inner quotes do not confuse
    # the line-level double-quote tracker.
    cmd_sub = re.compile(r"\$\([^()]*\)")
    normalized = line
    while True:
        next_line = cmd_sub.sub("CMD_SUB", normalized)
        if next_line == normalized:
            break
        normalized = next_line

    out = []
    for m in re.finditer(r"\$\{[^}]+\}", normalized):
        if not _is_within_double_quotes(normalized, m.start()):
            out.append(m.group(0))
    return out


def test_core_and_entrypoint_scripts_set_pipefail():
    scripts = _core_and_entrypoint_scripts()
    assert scripts, "No core/entrypoint scripts were found."
    for script in scripts:
        header = _strict_mode_header(script)
        has_pipefail = ("set -eo pipefail" in header) or ("set -euo pipefail" in header)
        assert has_pipefail, f"Missing pipefail guard in script header: {script}"


def test_core_and_entrypoint_scripts_use_strict_euo_pipefail():
    scripts = _core_and_entrypoint_scripts()
    assert scripts, "No core/entrypoint scripts were found."
    for script in scripts:
        header = _strict_mode_header(script)
        assert "set -euo pipefail" in header, f"Use strict mode (set -euo pipefail): {script}"


def test_non_library_workflow_shell_scripts_use_strict_euo_pipefail():
    allowed_non_strict = {
        WORKFLOW_DIR / "gg_common_params.sh",
        WORKFLOW_DIR / "support" / "gg_util.sh",
    }
    scripts = _workflow_shell_scripts()
    assert scripts, "No workflow shell scripts were found."
    for script in scripts:
        if script in allowed_non_strict:
            continue
        header = _strict_mode_header(script)
        assert "set -euo pipefail" in header, f"Use strict mode (set -euo pipefail): {script}"


def test_support_directory_has_no_numbered_duplicate_scripts():
    duplicates = sorted((WORKFLOW_DIR / "support").glob("* 2.*"))
    assert not duplicates, f"Remove accidental duplicate support scripts: {duplicates}"


def test_gg_versions_lists_support_dir_without_undefined_alias():
    text = _read_text(WORKFLOW_DIR / "support" / "gg_versions.sh")
    assert "dir_myscript" not in text
    assert 'ls -la "${gg_support_dir}"' in text


def test_gg_versions_uses_shared_core_bootstrap_runtime():
    text = _read_text(WORKFLOW_DIR / "support" / "gg_versions.sh")
    assert 'source "${gg_core_bootstrap}"' in text
    assert 'gg_bootstrap_core_runtime "${BASH_SOURCE[0]:-$0}" "" 1 1' in text
    assert 'gg_workspace_dir="/workspace"' not in text
    assert 'source "${gg_support_dir}/gg_util.sh"' not in text


def test_entrypoint_bootstrap_sets_python_pycacheprefix_outside_repo():
    text = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_bootstrap.sh")
    body = _function_body(text, "gg_configure_python_pycacheprefix")
    assert 'default_pycache_prefix="${TMPDIR:-/tmp}/genegalleon_pycache"' in body
    assert 'mkdir -p -- "${default_pycache_prefix}" 2>/dev/null || true' in body
    assert 'export PYTHONPYCACHEPREFIX="${default_pycache_prefix}"' in body
    init_body = _function_body(text, "gg_entrypoint_initialize")
    assert "gg_configure_python_pycacheprefix" in init_body


def test_core_bootstrap_sets_python_pycacheprefix_under_tmp():
    text = _read_text(WORKFLOW_DIR / "support" / "gg_core_bootstrap.sh")
    body = _function_body(text, "gg_configure_python_pycacheprefix_from_core")
    assert 'default_pycache_prefix="/tmp/genegalleon_pycache"' in body
    assert 'mkdir -p -- "${default_pycache_prefix}" 2>/dev/null || true' in body
    assert 'export PYTHONPYCACHEPREFIX="${default_pycache_prefix}"' in body
    runtime_body = _function_body(text, "gg_bootstrap_core_runtime")
    assert "gg_configure_python_pycacheprefix_from_core" in runtime_body


def test_busco_download_lock_is_per_lineage_shared_artifact_lock():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    assert 'lock_file="${runtime_busco_db}/lineages/.busco_${busco_lineage}.download.lock"' in text
    assert 'lock_file="${runtime_busco_db}/locks/busco_downloads.lock"' not in text


def test_no_find_exec_rm_rf_in_workflow_shell_scripts():
    for script in _workflow_shell_scripts():
        text = _read_text(script)
        assert "-exec rm -rf" not in text, f"Dangerous find+rm pattern found: {script}"


def test_no_bare_rm_glob_in_set_e_core_scripts():
    bare_rm_glob = re.compile(
        r"^[ \t]*rm[ \t]+(?!-f\b|-rf\b|-fr\b|--\b|-[^\n]*f\b).*[*?]",
        re.MULTILINE,
    )
    for script, text in _set_e_scripts():
        assert bare_rm_glob.search(text) is None, f"Bare rm glob found in {script}"


def test_no_plain_rm_r_in_set_e_scripts():
    bare_rm_r = re.compile(r"^[ \t]*rm[ \t]+-r(?!f)\b", re.MULTILINE)
    for script, text in _set_e_scripts():
        assert bare_rm_r.search(text) is None, f"Use rm -rf -- for recursive deletion: {script}"


def test_no_rm_rf_glob_in_set_e_scripts():
    rm_rf_glob = re.compile(r"^[ \t]*rm[ \t]+-[^\n]*r[^\n]*f[^\n]*[ \t]+[^\n]*[*?]", re.MULTILINE)
    for script, text in _set_e_scripts():
        assert rm_rf_glob.search(text) is None, f"Use nullglob+array guard instead of rm -rf glob: {script}"


def test_workflow_and_container_scripts_use_rm_rf_with_double_dash():
    rm_rf_without_dd = re.compile(r"\brm[ \t]+-rf(?![ \t]+--)")
    scripts = _workflow_shell_scripts() + _container_shell_scripts()
    for script in scripts:
        text = _read_text(script)
        assert rm_rf_without_dd.search(text) is None, f"Use rm -rf -- for option-safe recursive delete: {script}"


def test_workflow_and_container_scripts_do_not_use_for_in_seq_command_substitution():
    pattern = re.compile(r"for[ \t]+[A-Za-z_][A-Za-z0-9_]*[ \t]+in[ \t]+\$\(\s*seq\b")
    scripts = _workflow_shell_scripts() + _container_shell_scripts()
    for script in scripts:
        text = _read_text(script)
        assert pattern.search(text) is None, (
            f"Use arithmetic for-loops instead of `for ... in $(seq ...)`: {script}"
        )


def test_workflow_and_container_scripts_do_not_use_for_in_command_substitution():
    pattern = re.compile(r"for[ \t]+[A-Za-z_][A-Za-z0-9_]*[ \t]+in[ \t]+\$\(")
    scripts = _workflow_shell_scripts() + _container_shell_scripts()
    for script in scripts:
        text = _read_text(script)
        assert pattern.search(text) is None, (
            f"Use arrays/mapfile instead of `for ... in $(...)` to avoid word-splitting bugs: {script}"
        )


def test_workflow_and_container_scripts_avoid_nonportable_awk_match_array_capture():
    pattern = re.compile(r"match\([^)]*,[^)]*,[ \t]*[A-Za-z_][A-Za-z0-9_]*\)")
    scripts = _workflow_shell_scripts() + _container_shell_scripts()
    for script in scripts:
        text = _read_text(script)
        assert pattern.search(text) is None, (
            f"Avoid awk match(..., ..., array) for portability: {script}"
        )


def test_no_unquoted_cd_var_in_set_e_scripts():
    unquoted_cd_var = re.compile(
        r"^[ \t]*cd[ \t]+(\$[A-Za-z_][A-Za-z0-9_]*|\$\{[A-Za-z_][A-Za-z0-9_]*\})\b",
        re.MULTILINE,
    )
    for script, text in _set_e_scripts():
        assert unquoted_cd_var.search(text) is None, f"Unquoted variable in cd command: {script}"


def test_disable_if_no_input_file_uses_quoted_expansions():
    for script in sorted(CORE_DIR.glob("*.sh")):
        for lineno, line in enumerate(_read_text(script).splitlines(), start=1):
            if "disable_if_no_input_file" not in line:
                continue
            bad = _unquoted_brace_expansions(line)
            assert not bad, f"Unquoted variable expansion in disable_if_no_input_file: {script}:{lineno}: {line}"


def test_cp_out_and_mv_out_use_quoted_expansions():
    for script in sorted(CORE_DIR.glob("*.sh")):
        for lineno, line in enumerate(_read_text(script).splitlines(), start=1):
            stripped = line.lstrip()
            if not (stripped.startswith("cp_out ") or stripped.startswith("mv_out ")):
                continue
            bad = _unquoted_brace_expansions(line)
            assert not bad, f"Unquoted variable expansion in cp_out/mv_out: {script}:{lineno}: {line}"


def test_safe_directory_clear_helper_is_used_for_mcmctree_dirs():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    assert "clear_directory_contents_safe()" in text
    assert 'clear_directory_contents_safe "$(dirname "${file_iq2mc_ctl}")"' in text
    assert 'clear_directory_contents_safe "${dir_mcmctree2}"' in text


def test_support_python_shebangs_use_python3():
    support_dir = WORKFLOW_DIR / "support"
    scripts = sorted(support_dir.glob("*.py"))
    assert scripts, "No support Python scripts were found."
    for script in scripts:
        first_line = _read_text(script).splitlines()[0]
        if first_line.startswith("#!"):
            assert re.match(r"^#! ?/usr/bin/env python3$", first_line), f"Use python3 shebang: {script}"


def test_progress_summary_entrypoint_runs_core_script_in_container():
    entrypoint = WORKFLOW_DIR / "gg_progress_summary_entrypoint.sh"
    text = _read_text(entrypoint)
    assert 'gg_run_container_shell_script "${gg_container_image_path}" "${gg_core_dir}/gg_progress_summary_core.sh"' in text
    assert "orthogroup_output_summary.py" not in text
    assert "transcriptome_assembly_output_summary.py" not in text


def test_no_known_unquoted_query_gene_file_expansions():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "head --bytes 1 ${file_query_gene}",
        "seqkit stats --tabular ${file_query_gene}",
        "seqkit translate --allow-unknown-codon --transl-table ${genetic_code} --threads ${GG_TASK_CPUS} ${file_query_gene}",
        "cp_out ${file_query_gene} ${dir_output_active}/query_gene/$(basename \"${file_query_gene}\")",
    ]
    for token in banned_tokens:
        assert token not in text, f"Found unquoted file_query_gene expansion: {token}"


def test_no_known_unquoted_file_sp_trait_expansions():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "binarize_species_trait ${file_sp_trait} species_trait_binary.tsv",
        "if [[ -s ${file_sp_trait} ]]; then",
        "if head -n1 ${file_sp_trait} | grep -q ' '; then",
        "sed '2,$ s/\\t/_.*\\t/' ${file_sp_trait} > \"foreground.tsv\"",
        "\t${file_sp_trait}",
    ]
    for token in banned_tokens:
        assert token not in text, f"Found unquoted file_sp_trait expansion: {token}"


def test_no_known_unquoted_array_appends_for_output_file_lists():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "codeml_output_files+=(${codeml_output_file})",
        "missing_files+=(${codeml_output_file})",
        "relax_output_files+=(${relax_output_file})",
        "missing_files+=(${relax_output_file})",
    ]
    for token in banned_tokens:
        assert token not in text, f"Found unquoted array append: {token}"


def test_no_unquoted_dir_script_invocations_in_core_scripts():
    banned_tokens = [
        "python ${gg_support_dir}/",
        "Rscript ${gg_support_dir}/",
        "bash ${gg_support_dir}/",
    ]
    for script in sorted(CORE_DIR.glob("*.sh")):
        text = _read_text(script)
        for token in banned_tokens:
            assert token not in text, f"Found unquoted gg_support_dir invocation in {script}: {token}"


def test_no_unquoted_long_option_value_expansions_in_core_scripts():
    option_value = re.compile(r"--[A-Za-z0-9_.-]+=\$\{[^}]+\}")
    for script in sorted(CORE_DIR.glob("*.sh")):
        for lineno, line in enumerate(_read_text(script).splitlines(), start=1):
            stripped = line.lstrip()
            if stripped.startswith("#"):
                continue
            for match in option_value.finditer(line):
                assert _is_within_double_quotes(line, match.start()), (
                    f"Unquoted --key=${{var}} expansion in {script}:{lineno}: {line}"
                )


def test_no_unquoted_path_like_option_expansions_in_core_scripts():
    brace_pat = re.compile(r"--[A-Za-z0-9_.-]+\s+\$\{(?:file_|dir_|sp_)[^}]*\}")
    plain_pat = re.compile(r"--[A-Za-z0-9_.-]+\s+\$(?:file_|dir_|sp_)[A-Za-z0-9_]*")
    for script in sorted(CORE_DIR.glob("*.sh")):
        for lineno, line in enumerate(_read_text(script).splitlines(), start=1):
            stripped = line.lstrip()
            if stripped.startswith("#"):
                continue
            assert brace_pat.search(line) is None, (
                f"Unquoted path-like option expansion in {script}:{lineno}: {line}"
            )
            assert plain_pat.search(line) is None, (
                f"Unquoted path-like option expansion in {script}:{lineno}: {line}"
            )


def test_no_unquoted_infile_option_expansions_in_core_scripts():
    pattern = re.compile(r"--infile\s+\$\{[^}]+\}")
    for script in sorted(CORE_DIR.glob("*.sh")):
        for lineno, line in enumerate(_read_text(script).splitlines(), start=1):
            stripped = line.lstrip()
            if stripped.startswith("#"):
                continue
            assert pattern.search(line) is None, (
                f"Unquoted --infile expansion in {script}:{lineno}: {line}"
            )


def test_no_unquoted_outfile_option_expansions_in_core_scripts():
    pattern = re.compile(r"--outfile\s+\$\{[^}]+\}")
    for script in sorted(CORE_DIR.glob("*.sh")):
        for lineno, line in enumerate(_read_text(script).splitlines(), start=1):
            stripped = line.lstrip()
            if stripped.startswith("#"):
                continue
            assert pattern.search(line) is None, (
                f"Unquoted --outfile expansion in {script}:{lineno}: {line}"
            )


def test_no_eval_in_gg_test_shell_commands():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    assert 'eval "${command_text}"' not in text


def test_workflow_scripts_avoid_which_subshell_usage():
    pattern = re.compile(r"\$\(\s*which\s+")
    for script in _workflow_shell_scripts():
        text = _read_text(script)
        assert pattern.search(text) is None, f"Use command -v instead of which in subshell: {script}"


def test_print_softmasked_percentage_handles_zero_length_input_safely():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "print_softmasked_percentage")
    assert "tr -d '\\n'" in body
    assert 'if [[ "${num_total_bp}" -eq 0 ]]; then' in body
    assert 'echo "0.0% masked (0/0 bp)"' in body
    assert 'python -c' in body
    assert ' ${num_masked_bp} ${num_total_bp}' not in body


def test_is_output_older_than_inputs_uses_compgen_variable_listing():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "is_output_older_than_inputs")
    assert 'compgen -A variable | grep -E -- "${input_file_variable_regex}"' in body
    assert 'set | grep "${input_file_variable_regex}"' not in body


def test_set_singularity_command_supports_apptainer_fallback():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    assert "gg_detect_container_runtime_binary()" in text
    assert "command -v singularity" in text
    assert "command -v apptainer" in text
    assert 'Neither singularity nor apptainer was found on PATH.' in text


def test_site_runtime_exec_command_uses_output_parameter():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    site_runtime_path = WORKFLOW_DIR / "support" / "gg_site_runtime.sh"
    util_text = _read_text(util_path)
    site_text = _read_text(site_runtime_path)
    site_body = _function_body(site_text, "gg_site_container_shell_command")
    helper_body = _function_body(site_text, "gg_set_command_array")
    assert 'local out_var=${2:-}' in site_body
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
    text = _read_text(util_path)
    assert 'set | grep "^SINGULARITY"' not in text
    assert "forwarded_container_env_vars" in text


def test_all_entrypoints_call_set_singularity_command():
    entrypoints = sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    assert entrypoints, "No entrypoint scripts were found."
    for script in entrypoints:
        text = _read_text(script)
        assert (
            "set_singularity_command" in text
            or "gg_entrypoint_prepare_container_runtime" in text
        ), f"Missing container runtime preparation call: {script}"


def test_gg_trigger_versions_dump_is_runtime_agnostic():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "gg_trigger_versions_dump")
    assert 'export SINGULARITYENV_GG_VERSION="${gg_version}"' in body
    assert 'export APPTAINERENV_GG_VERSION="${gg_version}"' in body
    assert 'command -v "${container_runtime_bin}"' in body
    assert 'container_runtime_bin="$(gg_container_shell_command_runtime_bin || true)"' in body
    assert 'gg_run_container_shell_script "${gg_container_image_path}" "${versions_script}"' in body
    assert '"${container_runtime_bin}" inspect "${gg_container_image_path}"' in body
    assert '"${container_runtime_bin}" version || {' in body
    assert 'singularity inspect "${gg_container_image_path}"' not in body
    assert 'singularity version || {' not in body


def test_entrypoint_activate_container_runtime_prints_version_summary():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "gg_entrypoint_activate_container_runtime")
    assert "gg_entrypoint_print_version_summary" in body


def test_container_build_metadata_includes_repo_version_label():
    dockerfile = _read_text(REPO_ROOT / "container" / "Dockerfile")
    buildx = _read_text(REPO_ROOT / "container" / "buildx.sh")
    local_build = _read_text(REPO_ROOT / "container" / "apptainer_local_build.sh")
    definition_template = _read_text(REPO_ROOT / "container" / "apptainer_local_build.def.template")

    assert 'org.opencontainers.image.version="${GG_VERSION}"' in dockerfile
    assert '--build-arg GG_VERSION="${gg_version}"' in buildx
    assert 's|@@GG_VERSION@@|' in local_build
    assert "org.opencontainers.image.version @@GG_VERSION@@" in definition_template


def test_run_container_shell_script_uses_exec_with_bash_stdin_bridge():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "gg_run_container_shell_script")
    assert 'subcommand=$(gg_container_shell_command_subcommand || true)' in body
    assert '"${singularity_command[@]}" "${image_path}" bash -s -- < "${script_path}"' in body
    assert '${singularity_command} "${image_path}" bash -s -- < "${script_path}"' in body
    assert '"${singularity_command[@]}" "${image_path}" < "${script_path}"' in body


def test_entrypoints_stream_core_scripts_via_container_shell_helper():
    entrypoints = sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    assert entrypoints, "No entrypoint scripts were found."
    for script in entrypoints:
        text = _read_text(script)
        if 'gg_core_dir=' not in text:
            continue
        assert 'gg_run_container_shell_script "${gg_container_image_path}"' in text
        assert '${singularity_command} "${gg_container_image_path}" <' not in text


def test_progress_summary_entrypoint_uses_auto_forwarding_and_normalized_nslots():
    entrypoint = WORKFLOW_DIR / "gg_progress_summary_entrypoint.sh"
    text = _read_text(entrypoint)

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
    text = _read_text(entrypoint)

    assert "gg_apply_named_env_overrides \\" in text
    assert "provider GG_INPUT_PROVIDER" in text
    assert "trait_profile GG_INPUT_TRAIT_PROFILE" in text
    assert "for gg_input_var_name in ${!GG_INPUT_@}; do" not in text
    assert 'gg_forward_env_vars_with_prefix_to_container_env "GG_INPUT_MAX_CONCURRENT_DOWNLOADS_"' in text
    assert 'export "SINGULARITYENV_${gg_input_var_name}=${!gg_input_var_name}"' not in text
    assert 'export "APPTAINERENV_${gg_input_var_name}=${!gg_input_var_name}"' not in text


def test_input_generation_trait_profile_preset_is_wired():
    entrypoint = WORKFLOW_DIR / "gg_input_generation_entrypoint.sh"
    core = WORKFLOW_DIR / "core" / "gg_input_generation_core.sh"
    entry_text = _read_text(entrypoint)
    core_text = _read_text(core)

    assert 'trait_profile="none"' in entry_text
    assert "trait_profile GG_INPUT_TRAIT_PROFILE" in entry_text
    assert "GG_INPUT_" not in core_text
    assert "apply_env_override()" not in core_text
    assert 'case "${trait_profile}" in' in core_text
    assert "gift_starter" in core_text


def test_gg_util_has_common_forward_config_export_helpers():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    assert "gg_export_var_to_container_env_if_set()" in text
    assert "gg_apply_named_env_overrides()" in text
    assert "gg_forward_env_vars_with_prefix_to_container_env()" in text
    assert "gg_print_named_config_summary()" in text
    assert "gg_print_registered_config_summary()" in text
    assert "gg_print_entrypoint_config_summary()" in text
    assert "gg_require_versions_dump()" in text
    assert "gg_resolve_physical_path()" in text
    body = _function_body(text, "forward_config_vars_to_container_env")
    assert "gg_print_entrypoint_config_vars" in body
    assert 'gg_export_var_to_container_env_if_set "gg_debug_mode"' in body


def test_entrypoint_config_var_registry_covers_all_entrypoints():
    registry_text = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")
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
        header = _entrypoint_scheduler_header(script)
        positions = []
        for token in section_tokens:
            assert token in header, f"Missing scheduler header section {token!r} in {script}"
            positions.append(header.index(token))
        assert positions == sorted(positions), f"Scheduler header section order drifted in {script}"


def test_entrypoints_use_active_scheduler_directives_in_header_template():
    entrypoints = sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    assert entrypoints, "No entrypoint scripts were found."
    for script in entrypoints:
        header = _entrypoint_scheduler_header(script)
        assert "##PBS" not in header, f"Use active #PBS directives in {script}"
        assert "##SBATCH -N" not in header, f"Drop legacy commented node-count example from {script}"
        assert "##SBATCH -n" not in header, f"Drop legacy commented task-count example from {script}"
        assert "#PBS -S /bin/bash" in header, f"Missing PBS shell directive in {script}"
        assert "#PBS -V" in header, f"Missing PBS environment export directive in {script}"
        assert "#SBATCH --ignore-pbs" in header, f"Missing Slurm PBS-ignore guard in {script}"


def test_entrypoint_scheduler_directives_are_left_aligned():
    entrypoints = sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    assert entrypoints, "No entrypoint scripts were found."
    bad_lines = []
    for script in entrypoints:
        for lineno, line in enumerate(_entrypoint_scheduler_header(script).splitlines(), start=1):
            if re.match(r"^ (?:#SBATCH|#PBS|#\$)", line):
                bad_lines.append(f"{script}:{lineno}: {line}")
    assert not bad_lines, "Left-align scheduler directives in entrypoint headers:\n" + "\n".join(bad_lines)


def test_gg_util_uses_explicit_conda_shell_initialization():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "gg_initialize_conda_shell")
    assert "/home/.bashrc" not in text
    assert "micromamba shell hook --shell bash" in body
    assert 'conda() {' in text
    assert 'micromamba "$@"' in text
    assert "source /opt/conda/etc/profile.d/conda.sh" in text


def test_gg_add_container_bind_mount_skips_duplicate_destinations():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "gg_add_container_bind_mount")
    assert "gg_container_mount_destination()" in text
    assert "gg_container_bind_destination_exists()" in text
    assert 'if gg_container_bind_destination_exists "${mount_spec}"; then' in body
    assert 'GG_CONTAINER_BIND_MOUNTS=$(gg_csv_prepend "${mount_spec}" "${GG_CONTAINER_BIND_MOUNTS:-}")' in body
    assert "gg_sync_container_bind_envs" in body


def test_ensure_latest_jaspar_file_uses_set_e_safe_assignments():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "_prepare_latest_jaspar_file_locked")

    assert "if resolved_filename=$(_jaspar_find_latest_meme_filename_remote); then" in body
    assert "if resolved_filename=$(_jaspar_find_latest_meme_filename_local" in body
    assert "if ! resolved_path=$(_ensure_jaspar_file_named" in body

    unsafe_remote = re.compile(
        r"^[ \t]*resolved_filename=\$\(_jaspar_find_latest_meme_filename_remote\)\s*$",
        re.MULTILINE,
    )
    unsafe_local = re.compile(
        r"^[ \t]*resolved_filename=\$\(_jaspar_find_latest_meme_filename_local[^\n]*\)\s*$",
        re.MULTILINE,
    )
    unsafe_ensure = re.compile(
        r"^[ \t]*resolved_path=\$\(_ensure_jaspar_file_named[^\n]*\)\s*$",
        re.MULTILINE,
    )

    assert unsafe_remote.search(body) is None
    assert unsafe_local.search(body) is None
    assert unsafe_ensure.search(body) is None


def test_gene_convergence_entrypoint_forwards_plot_runtime_envs():
    entrypoint = WORKFLOW_DIR / "gg_gene_convergence_entrypoint.sh"
    text = _read_text(entrypoint)

    assert 'gg_export_var_to_container_env_if_set "PYMOL_HEADLESS"' in text
    assert 'gg_export_var_to_container_env_if_set "QT_QPA_PLATFORM"' in text
    assert 'export "SINGULARITYENV_PYMOL_HEADLESS=${PYMOL_HEADLESS}"' not in text
    assert 'export "APPTAINERENV_PYMOL_HEADLESS=${PYMOL_HEADLESS}"' not in text
    assert 'export "SINGULARITYENV_QT_QPA_PLATFORM=${QT_QPA_PLATFORM}"' not in text
    assert 'export "APPTAINERENV_QT_QPA_PLATFORM=${QT_QPA_PLATFORM}"' not in text


def test_gene_evolution_core_uses_run_fimo_not_legacy_run_meme_toggle():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert "run_meme=0" not in text
    assert "${run_meme}" not in text
    assert "run_fimo=0" in text
    assert "${run_fimo}" in text


def test_gene_evolution_scripts_do_not_use_removed_run_pgls_gene_tree_toggle():
    core_script = CORE_DIR / "gg_gene_evolution_core.sh"
    entrypoint = WORKFLOW_DIR / "gg_gene_evolution_entrypoint.sh"
    core_text = _read_text(core_script)
    entrypoint_text = _read_text(entrypoint)
    assert "run_pgls_gene_tree" not in core_text
    assert "run_pgls_gene_tree" not in entrypoint_text


def test_gene_evolution_core_passes_gg_task_cpus_to_kfl1ou():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert 'CPU_PER_HOST=' not in text
    assert 'cpu_pick="${GG_TASK_CPUS}"' not in text
    assert '--nslots="${GG_TASK_CPUS}"' in text
    assert 'taskset' not in text
    assert 'cpu_id=$(python -c' not in text
    assert '"${l1ou_cmd[@]}"' in text


def test_gene_evolution_core_uses_kfl1ou_wrapper_with_supported_args_only():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    block_start = text.index('task="l1ou"')
    block_end = text.index('mv_out fit_ind.RData', block_start)
    l1ou_block = text[block_start:block_end]

    assert 'detect_OU_shift_l1ou.r' not in l1ou_block
    assert 'Rscript "${gg_support_dir}/detect_OU_shift_kfl1ou.r"' in l1ou_block
    assert '--require_internal_node_labels="${require_internal_node_labels:-1}"' not in l1ou_block
    assert '--clade_collapse_similarity_method="${clade_collapse_similarity_method}"' not in l1ou_block
    assert '--clade_collapse_similarity_threshold="${clade_collapse_similarity_threshold}"' not in l1ou_block
    assert '--ceil_negative=0' not in l1ou_block
    assert '--replicate_sep="_"' in l1ou_block


def test_detect_ou_shift_kfl1ou_enables_measurement_error_by_default():
    script = WORKFLOW_DIR / "support" / "detect_OU_shift_kfl1ou.r"
    text = _read_text(script)
    assert 'measurement_error = TRUE' in text
    assert 'input_error = input_error_fit' in text


def test_gene_evolution_core_uses_explicit_ne_and_grouped_logic_for_tree_pruning_gate():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert "if [[ ! ${run_tree_pruning} -eq 1 && ${run_l1ou} -eq 1 ]]; then" not in text
    assert "if [[ ${run_tree_pruning} -ne 1 && ${run_l1ou} -eq 1 ]]; then" in text


def test_genome_evolution_core_uses_ne_for_busco_count_mismatch_checks():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "if [[ ! ${num_busco_ids} -eq ${num_singlecopy_fasta} && ${run_individual_get_fasta} -eq 1 ]]; then",
        "if [[ ! ${num_busco_ids} -eq ${num_mafft_fasta} && ${run_individual_mafft} -eq 1 ]]; then",
        "if [[ ! ${num_busco_ids} -eq ${num_trimal_fasta} && ${run_individual_trimal} -eq 1 ]]; then",
        "if [[ ! ${num_busco_ids} -eq ${num_iqtree_pep} && ${run_individual_iqtree_pep} -eq 1 ]]; then",
        "if [[ ! ${num_busco_ids} -eq ${num_iqtree_dna} && ${run_individual_iqtree_dna} -eq 1 ]]; then",
    ]
    expected_tokens = [
        "if [[ ${num_busco_ids} -ne ${num_singlecopy_fasta} && ${run_individual_get_fasta} -eq 1 ]]; then",
        "if [[ ${num_busco_ids} -ne ${num_mafft_fasta} && ${run_individual_mafft} -eq 1 ]]; then",
        "if [[ ${num_busco_ids} -ne ${num_trimal_fasta} && ${run_individual_trimal} -eq 1 ]]; then",
        "if [[ ${num_busco_ids} -ne ${num_iqtree_pep} && ${run_individual_iqtree_pep} -eq 1 ]]; then",
        "if [[ ${num_busco_ids} -ne ${num_iqtree_dna} && ${run_individual_iqtree_dna} -eq 1 ]]; then",
    ]
    for token in banned_tokens:
        assert token not in text
    for token in expected_tokens:
        assert token in text


def test_set_singularityenv_forwards_gg_common_variables():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "set_singularityenv")

    assert 'resolved_workspace_dir=$(gg_resolve_physical_path "${gg_workspace_dir}")' in body
    assert 'resolved_workflow_dir=$(gg_resolve_physical_path "${gg_workflow_dir}")' in body
    assert 'resolved_container_image_path=$(gg_resolve_physical_path "${gg_container_image_path}")' in body
    assert 'gg_workspace_dir="${resolved_workspace_dir}"' in body
    assert 'gg_workflow_dir="${resolved_workflow_dir}"' in body
    assert 'gg_container_image_path="${resolved_container_image_path}"' in body
    assert 'resolved_workspace_layout=$(gg_resolve_workspace_layout "${gg_workspace_dir}")' in body
    assert "export SINGULARITYENV_PYTHONPYCACHEPREFIX=/tmp/genegalleon_pycache" in body
    assert "export APPTAINERENV_PYTHONPYCACHEPREFIX=/tmp/genegalleon_pycache" in body
    assert "for gg_common_var_name in ${!GG_COMMON_@}; do" in body
    assert 'export "SINGULARITYENV_${gg_common_var_name}=${!gg_common_var_name}"' in body
    assert 'export "APPTAINERENV_${gg_common_var_name}=${!gg_common_var_name}"' in body


def test_cp_out_and_mv_out_use_option_safe_copy_and_move():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    cp_body = _function_body(text, "cp_out")
    mv_body = _function_body(text, "mv_out")
    assert 'cp -- "$@"' in cp_body
    assert 'cp "$@"' not in cp_body
    assert 'mv -- "$@"' in mv_body
    assert 'mv "$@"' not in mv_body


def test_cp_out_and_mv_out_prepare_destination_dir_for_multisource_or_trailing_slash():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    cp_body = _function_body(text, "cp_out")
    mv_body = _function_body(text, "mv_out")
    expected_guard = 'if [[ $# -gt 2 || "${dest}" == */ ]]; then'
    assert expected_guard in cp_body
    assert expected_guard in mv_body
    assert 'ensure_dir "${dest%/}"' in cp_body
    assert 'ensure_dir "${dest%/}"' in mv_body
    assert 'ensure_parent_dir "${dest}"' in cp_body
    assert 'ensure_parent_dir "${dest}"' in mv_body


def test_no_if_bracket_dollar_question_nonzero_checks():
    pattern = re.compile(r"if[ \t]+\[\[[^\n]*\$\?[^\n]*\]\]", re.MULTILINE)
    for script in _workflow_shell_scripts():
        text = _read_text(script)
        assert pattern.search(text) is None, f"Use direct command checks instead of $? inside [[ ... ]]: {script}"


def test_no_if_not_command_then_status_assignment_from_dollar_question():
    # Pattern that loses the original command exit code under `if ! cmd; then`.
    pattern = re.compile(
        r"if[ \t]+![^\n]*\n(?:[ \t]*#.*\n)*[ \t]*[A-Za-z_][A-Za-z0-9_]*[ \t]*=[ \t]*\$\?",
        re.MULTILINE,
    )
    for script in _workflow_shell_scripts():
        text = _read_text(script)
        assert pattern.search(text) is None, (
            f"Unsafe status capture after `if ! ...`; use set +e or branch-specific capture: {script}"
        )


def test_container_scripts_do_not_use_pipefail_fragile_grep_q_pipelines():
    pattern = re.compile(r"\|\s*grep\s+-q(?:\s|$)")
    for script in _container_shell_scripts():
        text = _read_text(script)
        assert pattern.search(text) is None, (
            f"Use awk/string matching instead of pipefail-fragile grep -q pipeline: {script}"
        )


def test_core_scripts_do_not_use_head_dash_dash_bytes():
    for script in sorted(CORE_DIR.glob("*.sh")):
        text = _read_text(script)
        assert "head --bytes" not in text, f"Use portable `head -c` instead of `head --bytes`: {script}"


def test_core_scripts_do_not_use_double_bracket_without_space_after_open():
    pattern = re.compile(r"\[\[\(")
    for script in sorted(CORE_DIR.glob("*.sh")):
        text = _read_text(script)
        assert pattern.search(text) is None, f"Use `[[ (` spacing for readability/stability: {script}"


def test_core_scripts_do_not_use_echo_dash_e():
    pattern = re.compile(r"^[ \t]*echo[ \t]+-e\b", re.MULTILINE)
    for script in sorted(CORE_DIR.glob("*.sh")):
        text = _read_text(script)
        assert pattern.search(text) is None, f"Use printf instead of `echo -e`: {script}"


def test_benchmark_r_scripts_do_not_hardcode_tmp_epgls_lib():
    support_dir = WORKFLOW_DIR / "support"
    benchmark_scripts = sorted(support_dir.glob("benchmark_*.R"))
    if not benchmark_scripts:
        return
    for script in benchmark_scripts:
        text = _read_text(script)
        assert '"/tmp/epgls_lib"' not in text, f"Hardcoded /tmp path found: {script}"


def test_gg_trigger_versions_dump_preserves_failure_exit_code_tracking():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "gg_trigger_versions_dump")
    assert "local versions_exit_code=0" in body
    assert "local block_exit_code=0" in body
    assert "versions_exit_code=$?" not in body
    assert "block_exit_code=$?" in body
    assert "if [[ ${block_exit_code} -ne 0 && ${versions_exit_code} -eq 0 ]]; then" in body


def test_gg_require_versions_dump_preserves_exit_code_under_set_e():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "gg_require_versions_dump")
    assert "local versions_exit_code=0" in body
    assert "gg_trigger_versions_dump" in body
    assert "versions_exit_code=$?" in body
    assert 'echo "gg_require_versions_dump: gg_versions trigger failed for ${trigger_name}." >&2' in body


def test_gg_trigger_versions_dump_reuses_one_log_per_container():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "gg_trigger_versions_dump")
    assert 'log_file="${versions_dir}/container.${container_key_hash}.versions.log"' in body
    assert 'if [[ -s "${log_file}" ]]; then' in body
    assert 'gg_trigger_versions_dump: skipped existing ${log_file}' in body
    assert 'timestamp=$(date' not in body


def test_gg_trigger_versions_dump_key_tracks_versions_script_and_repo_version():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "gg_trigger_versions_dump")
    assert 'container_key_seed="gg_container_image_path=${gg_container_image_path};runtime=${container_runtime_bin};gg_version=${gg_version}"' in body
    assert 'versions_script_hash=$(sha256sum "${versions_script}" | awk \'{print $1}\')' in body
    assert 'container_key_seed="${container_key_seed};versions_script_sha256=${versions_script_hash}"' in body


def test_gg_trigger_versions_dump_does_not_print_full_log_contents():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "gg_trigger_versions_dump")
    assert 'cat "${log_file}"' not in body
    assert 'cat "${failed_log_file}"' not in body


def test_entrypoints_require_versions_dump_success():
    entrypoints = sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    assert entrypoints, "No entrypoint scripts were found."
    for script in entrypoints:
        text = _read_text(script)
        assert 'Warning: gg_versions trigger failed.' not in text, f"Legacy warning-only versions dump handling remains: {script}"
        assert 'gg_require_versions_dump "${gg_entrypoint_name}"' in text, f"Entry point must fail on versions dump error: {script}"


def test_gg_prepare_cds_fasta_stream_pipes_to_cdskit_pad():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "gg_prepare_cds_fasta_stream")
    assert "| cdskit pad --codontable" in body
    assert "| cdskit pad" in body


def test_debug_tree_annotation_paths_do_not_use_legacy_gg_pipeline_repo():
    annotation_summary = _read_text(WORKFLOW_DIR / "support" / "annotation_summary.r")
    stat_branch2tree = _read_text(WORKFLOW_DIR / "support" / "stat_branch2tree_plot.r")
    legacy_path_token = "/repos/gg_pipeline/gg_pipeline/script/tree_annotation"
    assert legacy_path_token not in annotation_summary
    assert legacy_path_token not in stat_branch2tree


def test_entrypoints_define_fixed_tmp_cleanup_flags():
    expected_tokens = {
        "gg_genome_annotation_entrypoint.sh": "delete_tmp_dir=1",
        "gg_transcriptome_generation_entrypoint.sh": "delete_tmp_dir=1",
        "gg_gene_evolution_entrypoint.sh": "delete_tmp_dir=1",
        "gg_gene_evolution_entrypoint.sh#preexisting": "delete_preexisting_tmp_dir=1",
        "gg_genome_evolution_entrypoint.sh": "delete_tmp_dir=1",
    }

    for key, token in expected_tokens.items():
        script_name = key.split("#")[0]
        text = _read_text(WORKFLOW_DIR / script_name)
        assert token in text, f"Missing fixed tmp cleanup flag in {script_name}: {token}"


def test_entrypoints_forward_cleanup_flags_defined_outside_config_block():
    expected_tokens = {
        "gg_genome_annotation_entrypoint.sh": 'forward_config_vars_to_container_env "${gg_entrypoint_name}" "delete_tmp_dir"',
        "gg_transcriptome_generation_entrypoint.sh": 'forward_config_vars_to_container_env "${gg_entrypoint_name}" "delete_tmp_dir"',
        "gg_gene_evolution_entrypoint.sh": 'forward_config_vars_to_container_env "${gg_entrypoint_name}" "delete_tmp_dir" "delete_preexisting_tmp_dir"',
        "gg_gene_evolution_entrypoint.sh#preexisting": '"delete_preexisting_tmp_dir"',
    }
    for key, token in expected_tokens.items():
        script_name = key.split("#")[0]
        text = _read_text(WORKFLOW_DIR / script_name)
        assert token in text, f"Missing cleanup var forwarding in {script_name}: {token}"


def test_entrypoints_print_config_summary_near_startup():
    expected_tokens = {
        "gg_gene_convergence_entrypoint.sh": 'gg_entrypoint_print_config_summary_if_available "${gg_entrypoint_name}"',
        "gg_gene_database_entrypoint.sh": 'gg_entrypoint_print_config_summary_if_available "${gg_entrypoint_name}"',
        "gg_gene_evolution_entrypoint.sh": 'gg_entrypoint_print_config_summary_if_available "${gg_entrypoint_name}" "delete_tmp_dir" "delete_preexisting_tmp_dir"',
        "gg_genome_annotation_entrypoint.sh": 'gg_entrypoint_print_config_summary_if_available "${gg_entrypoint_name}" "delete_tmp_dir"',
        "gg_genome_evolution_entrypoint.sh": 'gg_entrypoint_print_config_summary_if_available "${gg_entrypoint_name}"',
        "gg_input_generation_entrypoint.sh": 'gg_entrypoint_print_config_summary_if_available "${gg_entrypoint_name}"',
        "gg_progress_summary_entrypoint.sh": 'gg_entrypoint_print_config_summary_if_available "${gg_entrypoint_name}"',
        "gg_transcriptome_generation_entrypoint.sh": 'gg_entrypoint_print_config_summary_if_available "${gg_entrypoint_name}" "delete_tmp_dir"',
    }
    for script_name, token in expected_tokens.items():
        text = _read_text(WORKFLOW_DIR / script_name)
        assert token in text, f"Missing config summary print in {script_name}: {token}"


def test_entrypoint_bootstrap_provides_config_summary_fallback_helper():
    text = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_bootstrap.sh")
    assert "gg_entrypoint_print_config_summary_if_available() {" in text
    assert 'declare -F gg_print_entrypoint_config_summary >/dev/null 2>&1' in text
    assert 'Config summary helper is unavailable in gg_util.sh; skipping entrypoint config summary for ${entrypoint_name}.' in text


def test_gene_convergence_entrypoint_does_not_define_unused_delete_tmp_dir():
    text = _read_text(WORKFLOW_DIR / "gg_gene_convergence_entrypoint.sh")
    assert "delete_tmp_dir=" not in text


def test_gene_evolution_entrypoint_allows_debug_runner_mode_overrides():
    text = _read_text(WORKFLOW_DIR / "gg_gene_evolution_entrypoint.sh")
    assert 'mode_gene_evolution="${mode_gene_evolution:-query2family}"' in text
    assert 'run_hyphy_relax="${run_hyphy_relax:-0}"' in text
    assert 'run_hyphy_relax_reversed="${run_hyphy_relax_reversed:-0}"' in text
    assert "mode_orthogroup=0 # Analyze OrthoFinder orthogroups" not in text
    assert "mode_query2family=1 # Analyze all homologs of input genelist" not in text


def test_genome_evolution_core_prints_effective_config_summary():
    text = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")
    assert "print_effective_genome_evolution_config_summary()" in text
    assert 'gg_print_registered_config_summary \\' in text
    assert '"effective config summary (gg_genome_evolution_core.sh)" \\' in text
    assert "print_effective_genome_evolution_config_summary" in text.split("root_species_tree()", 1)[0]


def test_entrypoints_with_exit_if_running_call_duplicate_guard():
    scripts = [
        "gg_gene_convergence_entrypoint.sh",
        "gg_genome_annotation_entrypoint.sh",
        "gg_gene_evolution_entrypoint.sh",
    ]
    for script_name in scripts:
        text = _read_text(WORKFLOW_DIR / script_name)
        assert "exit_if_running=" in text
        assert (
            "exit_if_running_qstat" in text
            or "gg_entrypoint_prepare_container_runtime 1" in text
        ), f"Missing duplicate-job guard call in {script_name}"


def test_download_pfam_helper_guards_output_dir_before_recursive_delete():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "_download_pfam_le_to_dir")
    assert 'if [[ -z "${output_dir}" || "${output_dir}" == "/" ]]; then' in body
    assert 'rm -rf -- "${output_dir}.tmp"' in body
    assert 'rm -rf -- "${output_dir}"' in body


def test_download_helpers_use_set_e_safe_command_guards():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)

    uniprot_body = _function_body(text, "_download_uniprot_sprot_to_prefix")
    assert 'if ! curl -fsSL "${uniprot_url}" | gzip -dc > "${pep_tmp}"; then' in uniprot_body
    assert 'if ! diamond makedb --in "${pep_tmp}" --db "${dmnd_tmp_prefix}"; then' in uniprot_body

    pfam_body = _function_body(text, "_download_pfam_le_to_dir")
    assert 'if ! curl -fsSL "${url}" -o "${archive_path}"; then' in pfam_body
    assert 'if ! tar -xzf "${archive_path}" -C "${tmp_dir}"; then' in pfam_body


def test_download_lock_helper_uses_shared_lock_metadata_and_heartbeat():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "gg_array_download_once")
    assert 'if gg_artifact_ready "${artifact_path}"; then' in body
    assert 'if ! gg_shared_lock_acquire "${lock_file}" "${description}"; then' in body
    assert 'gg_shared_lock_start_heartbeat "${lock_file}"' in body
    assert 'heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}' in body
    assert 'gg_shared_lock_stop_heartbeat "${heartbeat_pid}"' in body
    assert 'gg_shared_lock_release "${lock_file}"' in body
    assert '.dlock' not in body


def test_shared_lock_helpers_encode_owner_metadata_and_safe_heartbeat_reclaim_rules():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    assert '"format": "shared-lock-v2"' in text
    assert "os.O_WRONLY | os.O_CREAT | os.O_EXCL" in text
    assert 'touch -c -- "${lock_file}" 2>/dev/null || true' in text
    assert 'stale_reason="same_host_same_boot_dead_pid"' in text
    assert 'stale_reason="heartbeat_timeout"' in text
    assert 'waiting for shared lock: ${description} (${owner_summary})' in text
    assert 'timed out waiting for shared lock: ${description} (${owner_summary})' in text


def test_download_entrypoints_use_shared_lock_helper_for_busco_and_ete():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    busco_body = _function_body(text, "ensure_busco_download_path")
    ete_body = _function_body(text, "ensure_ete_taxonomy_db")
    assert 'runtime_ready_marker="${runtime_busco_lineage}/.download.ready"' in busco_body
    assert 'gg_array_download_once "${lock_file}" "${runtime_ready_marker}"' in busco_body
    assert 'gg_array_download_once "${lock_file}" "${db_file}" "ETE taxonomy DB"' in ete_body


def test_gene_evolution_core_uses_shared_lock_helpers_for_db_builds_and_shared_copies():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert "command -v flock" not in text
    assert " flock " not in text
    assert 'db_lock_file="${sp_cds_blastdb}.tblastn.build.lock"' in text
    assert 'db_lock_file="${sp_cds_blastdb}.diamond.build.lock"' in text
    assert 'gg_shared_lock_acquire "${db_lock_file}" "TBLASTN database build (${sp})"' in text
    assert 'gg_shared_lock_acquire "${db_lock_file}" "DIAMOND database build (${sp})"' in text
    assert 'gg_shared_lock_acquire "${lock_file}" "GeneRax species tree copy"' in text
    assert 'gg_shared_lock_acquire "${lock_file}" "parameter artifact copy (${file_to})"' in text


def test_genome_annotation_core_uses_shared_lock_for_species_cds_validation_stamp():
    script = CORE_DIR / "gg_genome_annotation_core.sh"
    text = _read_text(script)
    assert "command -v flock" not in text
    assert "species_cds_validation_lock_dir" not in text
    assert 'gg_shared_lock_acquire "${species_cds_validation_lock}" "species CDS validation stamp"' in text
    assert 'gg_shared_lock_start_heartbeat "${species_cds_validation_lock}"' in text
    assert 'heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}' in text


def test_ete_taxonomy_helper_uses_explicit_shared_taxdump_path():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    set_env_body = _function_body(text, "gg_set_taxonomy_cache_env")
    locked_body = _function_body(text, "_ensure_ete_taxonomy_db_locked")
    ensure_body = _function_body(text, "ensure_ete_taxonomy_db")
    assert 'workspace_taxonomy_taxdumpfile()' in text
    assert 'ensure_dir "${dir_taxonomy}/ete4"' in set_env_body
    assert 'export GG_TAXONOMY_TAXDUMPFILE="${dir_taxonomy}/taxdump.tar.gz"' in set_env_body
    assert 'os.makedirs(os.path.join(cache_dir, "ete4"), exist_ok=True)' in locked_body
    assert 'GG_TAXONOMY_TAXDUMPFILE="${taxdump_file}"' in text
    assert 'NCBITaxa(dbfile=db_file, taxdump_file=ensure_taxdump_file(), update=True)' in locked_body
    assert 'ensure_dir "${dir_taxonomy}/ete4"' in ensure_body
    assert 'taxdump_file=$(workspace_taxonomy_taxdumpfile "${gg_workspace_dir}")' in ensure_body
    assert '_ensure_ete_taxonomy_db_locked "${db_file}" "${taxdump_file}"' in ensure_body


def test_latest_jaspar_lock_uses_shared_lock_and_marker_resolution():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "ensure_latest_jaspar_file")
    assert 'if resolved_path=$(_resolve_latest_jaspar_path_from_marker "${latest_marker}" "${sys_dir}" "${runtime_dir}"); then' in body
    assert 'if ! gg_shared_lock_acquire "${lock_file}" "latest JASPAR motif file"; then' in body
    assert 'gg_shared_lock_start_heartbeat "${lock_file}"' in body
    assert 'heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}' in body
    assert 'gg_shared_lock_stop_heartbeat "${heartbeat_pid}"' in body
    assert 'gg_shared_lock_release "${lock_file}"' in body


def test_shared_lock_heartbeat_is_not_started_via_command_substitution():
    disallowed = 'heartbeat_pid=$(gg_shared_lock_start_heartbeat'
    for script in _workflow_shell_scripts():
      assert disallowed not in _read_text(script), (
          f"Do not start shared-lock heartbeat via command substitution: {script}"
      )


def test_pfam_helpers_use_only_new_runtime_layout_and_function_name():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    gene_core = CORE_DIR / "gg_gene_evolution_core.sh"
    util_text = _read_text(util_path)
    gene_text = _read_text(gene_core)

    assert "legacy_runtime_dir" not in util_text
    assert "downloads/Pfam_LE" not in util_text
    assert "ensure_pfam_domain_db()" not in util_text
    assert "ensure_pfam_domain_db" not in gene_text
    assert 'ensure_pfam_le_db "${gg_workspace_dir}"' in gene_text


def test_get_total_fastq_len_uses_bash3_compatible_read_loop_and_excludes_hidden_files():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "get_total_fastq_len")
    assert "mapfile -d '' -t files" not in body
    assert "while IFS= read -r -d '' f; do" in body
    assert 'find "${input_dir}" -type f ! -name \'.*\' -name "${regex}" -print0' in body


def test_gg_util_avoids_mapfile_for_host_bash_compatibility():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    assert "mapfile" not in text


def test_busco_dataset_download_merges_staged_directory_contents():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "_download_busco_lineage_to_runtime")
    assert 'gg_merge_directory_contents "busco_downloads" "${runtime_busco_db}"' in body
    assert "-exec mv -f -- {}" not in body


def test_shared_busco_summary_stage_normalizes_and_checks_species_set_before_collect():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    gate = "if [[ ${run_species_get_busco_summary} -ne 1 ]]; then"
    prepare = "prepare_species_tree_input_dir"
    normalize = 'normalize_busco_table_naming "${dir_species_busco_full}" "${dir_species_busco_short}"'
    check = 'if ! is_species_set_identical "${species_tree_input_dir}" "${dir_species_busco_full}"; then'
    collect = 'python "${gg_support_dir}/collect_common_BUSCO_genes.py" \\'
    assert "run_shared_busco_summary_stage() {" in text
    assert gate in text
    assert prepare in text
    assert normalize in text
    assert check in text
    assert collect in text
    gate_idx = text.index(gate)
    prepare_idx = text.index(prepare, gate_idx)
    normalize_idx = text.index(normalize, prepare_idx)
    check_idx = text.index(check, normalize_idx)
    collect_idx = text.index(collect, check_idx)
    assert gate_idx < prepare_idx < normalize_idx < check_idx < collect_idx


def test_genome_busco_summary_syncs_from_shared_summary_and_gates_busco_getfasta():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    sync_fn = "sync_genome_busco_summary_table_from_shared() {"
    cmp_stmt = 'cmp -s "${file_species_busco_summary_table}" "${file_genome_busco_summary_table}"'
    copy_stmt = 'cp_out "${file_species_busco_summary_table}" "${file_genome_busco_summary_table}"'
    sync_call = "sync_genome_busco_summary_table_from_shared || true"
    gate = 'disable_if_no_input_file "run_busco_getfasta" "${file_genome_busco_summary_table}"'
    assert sync_fn in text
    assert cmp_stmt in text
    assert copy_stmt in text
    assert sync_call in text
    assert gate in text
    assert text.index(sync_call) < text.index(gate)


def test_shared_species_busco_stage_runs_for_species_or_genome_busco_flags():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    stage_fn = "run_shared_species_busco_stage() {"
    gate = "if [[ ${run_species_busco} -ne 1 ]]; then"
    species_call = "run_shared_species_busco_stage"
    assert stage_fn in text
    assert gate in text
    assert text.count(species_call) >= 2


def test_busco_getfasta_step_is_gated_by_summary_table_presence():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    gate = 'disable_if_no_input_file "run_busco_getfasta" "${file_genome_busco_summary_table}"'
    step = "if [[ ${run_busco_getfasta} -eq 1 ]]; then"
    assert gate in text
    assert step in text
    assert text.index(gate) < text.index(step)


def test_genome_evolution_core_uses_safe_busco_summary_count_helper():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    assert "get_busco_summary_gene_count()" in text
    assert 'num_busco_ids=$(get_busco_summary_gene_count "${file_species_busco_summary_table}")' in text


def test_remove_empty_subdirs_calls_use_quoted_variable_arguments():
    pattern = re.compile(
        r"^[ \t]*remove_empty_subdirs[ \t]+(\$\{[A-Za-z_][A-Za-z0-9_]*\}|\$[A-Za-z_][A-Za-z0-9_]*)[ \t]*$",
        re.MULTILINE,
    )
    for script in sorted(CORE_DIR.glob("*.sh")):
        text = _read_text(script)
        assert pattern.search(text) is None, (
            f"Unquoted variable passed to remove_empty_subdirs in {script}"
        )


def test_codeml_step_guards_against_empty_trait_columns_before_paste():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    guard = "if [[ ${#colname_array[@]} -le 1 ]]; then"
    message = "No trait columns were detected in ${file_sp_trait}. Skipping ${task}."
    paste_cmd = "paste -d$'\\t' \"${header_files[@]}\" > combined_header.tsv"
    assert guard in text
    assert message in text
    assert paste_cmd in text
    assert text.index(guard) < text.index(paste_cmd)


def test_codeml_omega_array_access_uses_nounset_safe_defaults():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert "codeml_out_background_omega=${codeml_out_omegas[0]:-}" in text
    assert "codeml_out_foreground_omega=${codeml_out_omegas[1]:-}" in text
    assert "codeml_out_background_omega=${codeml_out_omegas[0]}" not in text
    assert "codeml_out_foreground_omega=${codeml_out_omegas[1]}" not in text


def test_genome_annotation_core_uses_clean_safe_mkdir_for_task_tmp_dirs():
    script = CORE_DIR / "gg_genome_annotation_core.sh"
    text = _read_text(script)
    assert "mkdir input_gff" not in text
    assert "if [[ -e input_gff ]]; then" in text
    assert "mkdir -p input_gff" in text
    assert "mkdir ./tmp" not in text
    assert "if [[ -e ./tmp ]]; then" in text
    assert "mkdir -p ./tmp" in text


def test_gene_evolution_core_cleans_stale_packaging_dirs_before_mkdir():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert 'if [[ -e "${og_id}.mapdnds_parameter" ]]; then' in text
    assert 'mkdir -p "${og_id}.mapdnds_parameter"' in text
    assert 'if [[ -e "${og_id}.iqtree.anc" ]]; then' in text
    assert 'mkdir -p "${og_id}.iqtree.anc"' in text


def test_genome_evolution_core_parses_mcmctree_constraints_with_validation_helper():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    assert "parse_mcmctree_constraint_record()" in text
    assert "mcmctree_divergence_time_constraints=()" in text
    assert 'parse_mcmctree_constraint_record "${mcmctree_divergence_time_constraints[0]}" mcmctree_params' in text
    assert 'for mdtc in "${mcmctree_divergence_time_constraints[@]}"; do' in text
    assert "mcmctree_params=(${mcmctree_divergence_time_constraints//,/ })" not in text
    assert "mcmctree_params=(${mdtc//,/ })" not in text


def test_genome_evolution_core_uses_rerun_safe_mkdir_for_orthogroup_grampa_tmp_input():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    assert "mkdir ./tmp.orthogroup_grampa_indir" not in text
    assert "mkdir -p ./tmp.orthogroup_grampa_indir" in text


def test_genome_evolution_core_quotes_grampa_output_and_cafe_option_values():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)

    assert 'busco_grampa "${dir_busco_rooted_nwk_dna}" "$(dirname "${file_busco_grampa_dna}")" ${file_busco_grampa_dna}' not in text
    assert 'busco_grampa "${dir_busco_rooted_nwk_pep}" "$(dirname "${file_busco_grampa_pep}")" ${file_busco_grampa_pep}' not in text
    assert 'busco_grampa "./tmp.orthogroup_grampa_indir" "$(dirname "${file_orthogroup_grampa}")" ${file_orthogroup_grampa}' not in text
    assert '--genecount ${file_orthogroup_genecount_selected}' not in text
    assert '--dated_species_tree ${file_dated_species_tree}' not in text
    assert "--max_size_differential ${max_size_differential_cafe}" not in text
    assert "--tree ${file_dated_species_tree}" not in text
    assert "--n_gamma_cats ${n_gamma_cats_cafe}" not in text
    assert "--cores ${GG_TASK_CPUS}" not in text
    assert "--output_prefix ${dir_cafe_output}" not in text

    assert 'busco_grampa "${dir_busco_rooted_nwk_dna}" "$(dirname "${file_busco_grampa_dna}")" "${file_busco_grampa_dna}"' in text
    assert 'busco_grampa "${dir_busco_rooted_nwk_pep}" "$(dirname "${file_busco_grampa_pep}")" "${file_busco_grampa_pep}"' in text
    assert 'busco_grampa "./tmp.orthogroup_grampa_indir" "$(dirname "${file_orthogroup_grampa}")" "${file_orthogroup_grampa}"' in text
    assert '--genecount "${file_orthogroup_genecount_selected}"' in text
    assert '--dated_species_tree "${file_dated_species_tree}"' in text
    assert '--max_size_differential "${max_size_differential_cafe}"' in text
    assert '--tree "${file_dated_species_tree}"' in text
    assert '--n_gamma_cats "${n_gamma_cats_cafe}"' in text
    assert '--cores "${GG_TASK_CPUS}"' in text
    assert '--output_prefix "${dir_cafe_output}"' in text


def test_genome_evolution_core_builds_grampa_arguments_with_array():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    assert "h1_param=\"-h1 " not in text
    assert "grampa_args=(" in text
    assert 'grampa_args+=(-h1 "${grampa_h1_normalized}")' in text
    assert 'grampa.py "${grampa_args[@]}"' in text


def test_genome_evolution_core_uses_option_safe_grep_for_orthogroup_id_removal():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    assert 'grep -v -Fx -- "${og_id}"' in text
    assert 'grep -v -Fx "${og_id}"' not in text


def test_core_scripts_use_mkdir_p_for_known_task_directories():
    genome_annotation = _read_text(CORE_DIR / "gg_genome_annotation_core.sh")
    transcriptome = _read_text(CORE_DIR / "gg_transcriptome_generation_core.sh")

    assert 'mkdir "${sp_ub}.genomescope"' not in genome_annotation
    assert 'mkdir -p "${sp_ub}.genomescope"' in genome_annotation
    assert 'mkdir "${sp_ub}.jcvi_dotplot"' not in genome_annotation
    assert 'mkdir -p "${sp_ub}.jcvi_dotplot"' in genome_annotation
    assert 'mkdir "${selected_fastq_dir}"' not in transcriptome
    assert 'mkdir -p "${selected_fastq_dir}"' in transcriptome
    assert 'mkdir "${assembly_input_fastq_dir}"' not in transcriptome
    assert 'mkdir -p "${assembly_input_fastq_dir}"' in transcriptome


def test_all_entrypoints_create_workspace_dir_before_cd():
    entrypoints = sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    assert entrypoints, "No entrypoint scripts were found."
    util_text = _read_text(WORKFLOW_DIR / "support" / "gg_util.sh")
    helper_body = _function_body(util_text, "gg_entrypoint_enter_workspace")
    assert 'mkdir -p "${gg_workspace_dir}"' in helper_body
    assert 'cd "${gg_workspace_dir}"' in helper_body
    assert helper_body.index('mkdir -p "${gg_workspace_dir}"') < helper_body.index('cd "${gg_workspace_dir}"')
    for script in entrypoints:
        text = _read_text(script)
        mkdir_token = 'mkdir -p "${gg_workspace_dir}"'
        cd_token = 'cd "${gg_workspace_dir}"'
        if "gg_entrypoint_enter_workspace" in text:
            continue
        assert mkdir_token in text, f"Missing workspace mkdir guard in {script}"
        assert cd_token in text, f"Missing workspace cd in {script}"
        assert text.index(mkdir_token) < text.index(cd_token), f"Workspace mkdir must come before cd in {script}"


def test_genome_evolution_core_quotes_orthogroup_iq2mc_and_busco_summary_options():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)

    banned_tokens = [
        "--dir_orthofinder_og ${dir_orthofinder_filtered}",
        "--dir_species_protein ${dir_sp_protein}",
        "--min_gene_num ${min_num_gene}",
        "--max_gene_num ${max_num_gene}",
        "--min_species_num ${min_num_species}",
        "--min_percent_species_coverage ${min_percent_species_coverage}",
        "--ncpu ${GG_TASK_CPUS}; then",
        "--orthofinder_hog_genecount ${dir_orthofinder_hog2og}/Orthogroups.GeneCount.tsv",
        "normalize_iq2mc_constraint_tree ${file_constrained_tree}",
        "if ! ${iq2mc_binary} \\",
        "--busco_outdir ${dir_species_busco_full}",
    ]
    for token in banned_tokens:
        assert token not in text, f"Found unquoted token: {token}"

    expected_tokens = [
        '--dir_orthofinder_og "${dir_orthofinder_filtered}"',
        '--dir_species_protein "${dir_sp_protein}"',
        '--min_gene_num "${min_num_gene}"',
        '--max_gene_num "${max_num_gene}"',
        '--min_species_num "${min_num_species}"',
        '--min_percent_species_coverage "${min_percent_species_coverage}"',
        '--ncpu "${GG_TASK_CPUS}"; then',
        '--orthofinder_hog_genecount "${dir_orthofinder_hog2og}/Orthogroups.GeneCount.tsv"',
        'normalize_iq2mc_constraint_tree "${file_constrained_tree}"',
        'if ! "${iq2mc_binary}" \\',
        '--busco_outdir "${dir_species_busco_full}"',
    ]
    for token in expected_tokens:
        assert token in text, f"Missing quoted token: {token}"

    iq2mc_start = text.index('if ! "${iq2mc_binary}" \\')
    iq2mc_end = text.index("--prefix iq2mc; then", iq2mc_start) + len("--prefix iq2mc; then")
    iq2mc_block = text[iq2mc_start:iq2mc_end]
    assert "-m ${nucleotide_model}" not in iq2mc_block
    assert "-te ${file_constrained_tree}" not in iq2mc_block
    assert "--mcmc-bds ${mcmc_birth_death_sampling}" not in iq2mc_block
    assert "--mcmc-clock ${mcmc_clock_model}" not in iq2mc_block
    assert "--mcmc-iter ${mcmc_burnin},${mcmc_sampfreq},${mcmc_nsample}" not in iq2mc_block
    assert "-T ${GG_TASK_CPUS}" not in iq2mc_block
    assert '-m "${nucleotide_model}"' in iq2mc_block
    assert '-te "${file_constrained_tree}"' in iq2mc_block
    assert '--mcmc-bds "${mcmc_birth_death_sampling}"' in iq2mc_block
    assert '--mcmc-clock "${mcmc_clock_model}"' in iq2mc_block
    assert '--mcmc-iter "${mcmc_burnin},${mcmc_sampfreq},${mcmc_nsample}"' in iq2mc_block
    assert '-T "${GG_TASK_CPUS}"' in iq2mc_block

    busco_summary_start = text.index('python "${gg_support_dir}/collect_common_BUSCO_genes.py" \\')
    busco_summary_end = text.index('--outfile "tmp.busco_summary_table.tsv"', busco_summary_start) + len('--outfile "tmp.busco_summary_table.tsv"')
    busco_summary_block = text[busco_summary_start:busco_summary_end]
    assert "--busco_outdir ${dir_species_busco_full}" not in busco_summary_block
    assert "--ncpu ${GG_TASK_CPUS}" not in busco_summary_block
    assert '--busco_outdir "${dir_species_busco_full}"' in busco_summary_block
    assert '--ncpu "${GG_TASK_CPUS}"' in busco_summary_block


def test_transcriptome_core_quotes_known_path_sensitive_options_and_symlinks():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)

    banned_tokens = [
        "--fastq_dir ${dir_species_fastq}",
        "--out_dir ${dir_tmp}",
        "--download_dir ${dir_amalgkit_download_dir}",
        "--metadata ${file_amalgkit_metadata}",
        "--rrna_filter ${amalgkit_rrna_filter}",
        "--contam_filter ${amalgkit_contam_filter}",
        "--contam_filter_rank ${contamination_removal_rank_for_amalgkit}",
        "--contam_filter_db ${dir_mmseqs2_db}/UniRef90_DB",
        "ln -s ${dir_amalgkit_getfastq_sp} \"./getfastq\"",
        "--fasta_file ${file_longestcds}",
        "--mmseqs2taxonomy_tsv ${file_longestcds_mmseqs2taxonomy}",
        "--fx2tab_tsv ${file_longestcds_fx2tab}",
        "--species_name ${sp_ub}",
        "--rank ${contamination_removal_rank_for_remove_contaminated_sequences}",
        "seqkit seq --threads ${GG_TASK_CPUS} ${file_isoform} --out-file \"busco_infile_cdna.fa\"",
        "seqkit seq --threads ${GG_TASK_CPUS} ${file_longestcds} --out-file \"busco_infile_cds.fa\"",
        "seqkit seq --threads ${GG_TASK_CPUS} ${file_longestcds_contamination_removal_fasta} --out-file \"busco_infile_cds.fa\"",
        "--lineage_dataset ${dir_busco_lineage}",
        "--download_path ${dir_busco_db}",
        "if [[ -e ${file_kallisto_reference_fasta} ]]; then",
        "ln -s ${file_kallisto_reference_fasta} ${file_reference_fasta_link}",
        "ln -s ${dir_amalgkit_quant}/${sp_ub} ./quant",
        'grep -e "${sp_space}" "./metadata/metadata.tsv"',
        "d.loc[:,'scientific_name']='${sp_ub}'",
        "mv_out ./metadata_private_fastq.tsv ./metadata.tsv",
    ]
    for token in banned_tokens:
        assert token not in text, f"Found unquoted transcriptome token: {token}"

    expected_tokens = [
        '--fastq_dir "${dir_species_fastq}"',
        '--out_dir "${dir_tmp}"',
        '--download_dir "${dir_amalgkit_download_dir}"',
        '--metadata "${file_amalgkit_metadata}"',
        '--rrna_filter "${amalgkit_rrna_filter}"',
        '--contam_filter "${amalgkit_contam_filter}"',
        '--contam_filter_rank "${contamination_removal_rank_for_amalgkit}"',
        '--contam_filter_db "${dir_mmseqs2_db}/UniRef90_DB"',
        'ln -s "${dir_amalgkit_getfastq_sp}" "./getfastq"',
        '--fasta_file "${file_longestcds}"',
        '--mmseqs2taxonomy_tsv "${file_longestcds_mmseqs2taxonomy}"',
        '--fx2tab_tsv "${file_longestcds_fx2tab}"',
        '--species_name "${contamination_removal_target_taxon:-${sp_ub}}"',
        '--rank "${contamination_removal_rank_for_remove_contaminated_sequences}"',
        'seqkit seq --threads "${GG_TASK_CPUS}" "${file_isoform}" --out-file "busco_infile_cdna.fa"',
        'seqkit seq --threads "${GG_TASK_CPUS}" "${file_longestcds}" --out-file "busco_infile_cds.fa"',
        'seqkit seq --threads "${GG_TASK_CPUS}" "${file_longestcds_contamination_removal_fasta}" --out-file "busco_infile_cds.fa"',
        '--lineage_dataset "${dir_busco_lineage}"',
        '--download_path "${dir_busco_db}"',
        'if [[ -e "${file_kallisto_reference_fasta}" ]]; then',
        'ln -s "${file_kallisto_reference_fasta}" "${file_reference_fasta_link}"',
        'ln -s "${dir_amalgkit_quant}/${sp_ub}" "./quant"',
        'grep -F -- "${sp_space}" "./metadata/metadata.tsv"',
        'mv_out "./metadata_private_fastq.tsv" "./metadata.tsv"',
    ]
    for token in expected_tokens:
        assert token in text, f"Missing quoted transcriptome token: {token}"


def test_transcriptome_core_sraid_metadata_filter_handles_zero_match_explicitly():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    assert 'grep -F -- "${sp_space}" "./metadata/metadata.tsv" || true' in text
    assert 'if [[ $(wc -l < "./metadata.tsv") -le 1 ]]; then' in text
    assert "No metadata rows matched species" in text


def test_transcriptome_core_requires_taxid_for_contam_filter():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    assert 'effective_amalgkit_contam_filter=' not in text
    assert 'Continuing with effective_amalgkit_contam_filter=no.' not in text
    assert 'amalgkit_contam_filter=yes requires a taxid column in metadata: ${file_amalgkit_metadata}. Exiting.' in text


def test_transcriptome_core_passes_download_dir_to_amalgkit_integrate():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    integrate_start = text.index("    amalgkit integrate \\")
    integrate_end = text.index('    mv_out "./metadata_private_fastq.tsv" "./metadata.tsv"', integrate_start)
    integrate_block = text[integrate_start:integrate_end]
    assert '--download_dir "${dir_amalgkit_download_dir}"' in integrate_block


def test_transcriptome_core_passes_shared_mmseqs_db_to_amalgkit_getfastq():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    getfastq_start = text.index("  if amalgkit getfastq \\")
    getfastq_end = text.index('    echo "amalgkit getfastq safely finished."', getfastq_start)
    getfastq_block = text[getfastq_start:getfastq_end]
    assert 'dir_mmseqs2_db="${gg_workspace_downloads_dir}/mmseqs2"' in text
    assert '--contam_filter_db "${dir_mmseqs2_db}/UniRef90_DB"' in getfastq_block


def test_transcriptome_core_invalidates_stale_cached_query_tables_on_species_prefix_change():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    body = _function_body(text, "invalidate_cached_query_table_if_prefix_mismatch")

    assert "first_query=$(awk -F '\\t' -v skip=\"${header_lines}\" 'NR > skip && $1 != \"\" { print $1; exit }' \"${table_file}\")" in body
    assert 'if [[ "${first_query}" != ${expected_prefix}* ]]; then' in body
    assert 'stale_file="${table_file}.stale.$(date +%Y%m%d%H%M%S)"' in body
    assert 'mv -f -- "${table_file}" "${stale_file}"' in body
    assert 'Archived stale file to: ${stale_file}' in body
    assert 'invalidate_cached_query_table_if_prefix_mismatch "${file_longestcds_fx2tab}" "${sp_ub}_" "${task}" 1' in text
    assert 'invalidate_cached_query_table_if_prefix_mismatch "${file_longestcds_mmseqs2taxonomy}" "${sp_ub}_" "${task}" 0' in text


def test_common_params_do_not_define_contamination_removal_rank():
    text = _read_text(WORKFLOW_DIR / "gg_common_params.sh")
    assert "GG_COMMON_CONTAMINATION_REMOVAL_RANK" not in text


def test_common_busco_lineage_defaults_to_auto():
    text = _read_text(WORKFLOW_DIR / "gg_common_params.sh")
    assert ': "${GG_COMMON_BUSCO_LINEAGE:=auto}"' in text


def test_common_params_define_reference_species_auto_only_once():
    text = _read_text(WORKFLOW_DIR / "gg_common_params.sh")
    assert ': "${GG_COMMON_REFERENCE_SPECIES:=auto}"' in text
    assert "GG_COMMON_ANNOTATION_SPECIES" not in text
    assert "GG_COMMON_ANNOTATION_REPRESENTATIVE_SPECIES" not in text
    assert "GG_COMMON_MCMCTREE_DIVERGENCE_TIME_CONSTRAINTS_STR" not in text
    assert "GG_COMMON_TREEVIS_CLADE_ORTHOLOG_PREFIX" not in text


def test_common_params_do_not_define_genome_evolution_specific_grampa_or_go_target():
    text = _read_text(WORKFLOW_DIR / "gg_common_params.sh")
    assert "GG_COMMON_OUTGROUP_LABELS" not in text
    assert "GG_COMMON_GRAMPA_H1" not in text
    assert "GG_COMMON_TARGET_BRANCH_GO" not in text


def test_core_scripts_resolve_busco_lineage_through_shared_helper():
    genome_annotation = _read_text(CORE_DIR / "gg_genome_annotation_core.sh")
    transcriptome = _read_text(CORE_DIR / "gg_transcriptome_generation_core.sh")
    genome_evolution = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")

    assert 'busco_lineage="${busco_lineage:-${GG_COMMON_BUSCO_LINEAGE:-auto}}"' in genome_annotation
    assert 'busco_lineage="${busco_lineage:-${GG_COMMON_BUSCO_LINEAGE:-auto}}"' in transcriptome
    assert 'busco_lineage="${busco_lineage:-${GG_COMMON_BUSCO_LINEAGE:-auto}}"' in genome_evolution
    assert 'gg_resolve_busco_lineage "${gg_workspace_dir}" "${busco_lineage}" "${sp_ub}"' in genome_annotation
    assert 'gg_resolve_busco_lineage "${gg_workspace_dir}" "${busco_lineage}" "${sp_ub}"' in transcriptome
    assert 'gg_resolve_busco_lineage "${gg_workspace_dir}" "${busco_lineage}" "$@"' in genome_evolution


def test_core_scripts_resolve_reference_species_through_shared_helper():
    gene_evolution = _read_text(CORE_DIR / "gg_gene_evolution_core.sh")
    genome_evolution = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")

    assert 'annotation_species="${annotation_species:-${GG_COMMON_REFERENCE_SPECIES:-auto}}"' in gene_evolution
    assert 'annotation_species="${annotation_species:-${GG_COMMON_REFERENCE_SPECIES:-auto}}"' in genome_evolution
    assert "GG_COMMON_ANNOTATION_SPECIES" not in gene_evolution
    assert "GG_COMMON_ANNOTATION_SPECIES" not in genome_evolution


def test_genome_evolution_uses_local_optional_grampa_and_go_target_parameters():
    entrypoint = _read_text(WORKFLOW_DIR / "gg_genome_evolution_entrypoint.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")
    core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")

    assert 'grampa_h1="" # Optional GRAMPA H1 hypothesis. Leave empty to skip GRAMPA steps. Example: "2" or "x,y,z".' in entrypoint
    assert 'target_branch_go="" # Optional GO-enrichment target branch. Leave empty to skip GO enrichment. Example: "<1>" or "Arabidopsis_thaliana".' in entrypoint
    assert "GG_COMMON_GRAMPA_H1" not in core
    assert "GG_COMMON_TARGET_BRANCH_GO" not in core
    assert 'grampa_h1="${grampa_h1:-}"' in core
    assert 'target_branch_go="${target_branch_go:-}"' in core
    assert 'Disabling GRAMPA tasks because grampa_h1 is empty. Set grampa_h1 in gg_genome_evolution_entrypoint.sh to enable them.' in core
    assert 'Disabling run_go_enrichment because target_branch_go is empty. Set target_branch_go in gg_genome_evolution_entrypoint.sh to enable it.' in core
    assert ': "${grampa_h1:?' not in core
    assert ': "${target_branch_go:?' not in core
    assert "grampa_h1" in config_vars
    assert "target_branch_go" in config_vars


def test_genome_evolution_uses_local_species_tree_rooting_parameter():
    entrypoint = _read_text(WORKFLOW_DIR / "gg_genome_evolution_entrypoint.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")
    core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")
    common = _read_text(WORKFLOW_DIR / "gg_common_params.sh")

    assert 'species_tree_rooting="taxonomy" # taxonomy[,ncbi[,opentree,timetree...]] | outgroup,GENUS_SPECIES[,GENUS_SPECIES...] | midpoint | mad | mv' in entrypoint
    assert "GG_COMMON_OUTGROUP_LABELS" not in common
    assert 'species_tree_rooting="${species_tree_rooting:-taxonomy}"' in core
    assert 'parse_species_tree_rooting "${species_tree_rooting}" species_tree_rooting_method species_tree_rooting_value' in core
    assert 'species_tree_rooting must be one of "outgroup,GENUS_SPECIES[,GENUS_SPECIES...]", "midpoint", "mad", "mv", or "taxonomy[,ncbi[,opentree,timetree...]]".' in core
    assert 'nwkit_root_args=(--method "${root_method}" --infile "${infile}" --outfile "${outfile}")' in core
    assert 'nwkit_root_args+=(--outgroup "${root_value}")' in core
    assert 'nwkit_root_args+=(--download_dir "${dir_nwkit_download_dir}")' in core
    assert "species_tree_rooting" in config_vars


def test_genome_evolution_supports_protein_input_mode_and_species_code_overrides():
    entrypoint = _read_text(WORKFLOW_DIR / "gg_genome_evolution_entrypoint.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")
    core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")

    assert 'input_sequence_mode="cds" # {cds,protein}; protein mode uses species_protein inputs or per-species CDS->protein translation with optional species_genetic_code/species_genetic_code.tsv overrides.' in entrypoint
    assert "input_sequence_mode" in config_vars
    assert 'input_sequence_mode="${input_sequence_mode:-cds}"' in core
    assert 'species_genetic_code_table_path() {' in core
    assert 'echo "${gg_workspace_input_dir}/species_genetic_code/species_genetic_code.tsv"' in core
    assert 'echo "${gg_workspace_input_dir}/species_protein"' in core
    assert 'prepare_species_genetic_code_table() {' in core
    assert 'lookup_species_genetic_code() {' in core
    assert 'check_species_protein_dir "${dir_sp_protein_input}"' in core
    assert 'species_genetic_code.tsv is ignored because species_protein inputs are provided' in core
    assert 'if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then' in core
    assert 'Ignoring species_protein inputs in cds mode: ${dir_sp_protein_input}' in core
    assert 'cds mode always generates temporary species_protein FASTA files from species_cds.' in core
    assert 'run_cds_translation must be 1 when species proteins need to be generated from species_cds.' in core
    assert 'Translation started: ${cds} (genetic_code=${species_code})' in core
    assert 'refresh_dir_for_shared_protein_input_signature "${dir_genome_evolution}" "genome_evolution" "${shared_protein_input_signature}"' in core
    assert 'mapfile -t annotation_species_candidates < <(gg_species_names_from_fasta_dir "${dir_sp_protein_input}")' in core
    protein_candidates_index = core.index('if [[ ${#annotation_species_candidates[@]} -eq 0 && "${input_sequence_mode}" == "protein" ]]; then')
    cds_candidates_index = core.index('if [[ ${#annotation_species_candidates[@]} -eq 0 ]]; then\n  mapfile -t annotation_species_candidates < <(gg_species_names_from_fasta_dir "${dir_sp_cds}")\nfi')
    assert protein_candidates_index < cds_candidates_index


def test_genome_evolution_protein_mode_disables_incompatible_dna_and_busco_steps():
    core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")

    assert 'if [[ "${input_sequence_mode}" == "protein" ]]; then' in core
    assert 'protein mode does not support undated_species_tree=${undated_species_tree}.' in core
    assert "Use iqtree_pep or astral_pep instead." in core
    assert "Disabling DNA-only species-tree steps in protein mode" in core
    assert "Disabling CDS-only dating steps in protein mode" in core
    assert "Disabling DNA-only BUSCO genome-evolution steps in protein mode" in core
    assert 'dir_species_busco_full="${gg_workspace_output_dir}/species_protein_busco_full"' in core
    assert 'dir_species_busco_short="${gg_workspace_output_dir}/species_protein_busco_short"' in core
    assert 'prepare_species_tree_input_dir' in core
    assert '--mode "${species_tree_busco_mode}"' in core
    assert 'outfile2="${dir_busco_fasta}/${busco_id}${genome_busco_fasta_suffix}"' in core
    assert 'outfile=${dir_busco_mafft}/${infile_base}${genome_busco_aln_suffix}' in core
    assert 'outfile="${dir_busco_trimal}/${infile_base}${genome_busco_trimal_suffix}"' in core


def test_annotation_and_transcriptome_use_local_contamination_removal_rank_parameter():
    annotation_entrypoint = _read_text(WORKFLOW_DIR / "gg_genome_annotation_entrypoint.sh")
    transcriptome_entrypoint = _read_text(WORKFLOW_DIR / "gg_transcriptome_generation_entrypoint.sh")
    annotation_core = _read_text(CORE_DIR / "gg_genome_annotation_core.sh")
    transcriptome_core = _read_text(CORE_DIR / "gg_transcriptome_generation_core.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")
    common = _read_text(WORKFLOW_DIR / "gg_common_params.sh")

    assert 'contamination_removal_rank="domain" # Taxonomic rank for contamination removal. Canonical value is domain; GeneGalleon normalizes tool-specific synonyms automatically.' in annotation_entrypoint
    assert 'contamination_removal_rank="domain" # Taxonomic rank for contamination removal. Canonical value is domain; GeneGalleon normalizes tool-specific synonyms automatically.' in transcriptome_entrypoint
    assert 'contamination_removal_target_taxon="${contamination_removal_target_taxon:-}" # Optional NCBI taxon name used as the lineage anchor for contamination removal (for example, Eukaryota when the sample species name is unknown).' in annotation_entrypoint
    assert 'contamination_removal_target_taxon="${contamination_removal_target_taxon:-}" # Optional NCBI taxon name used as the lineage anchor for contamination removal (for example, Eukaryota when the sample species name is unknown).' in transcriptome_entrypoint
    assert "GG_COMMON_CONTAMINATION_REMOVAL_RANK" not in common
    assert 'contamination_removal_rank="${contamination_removal_rank:-domain}"' in annotation_core
    assert 'contamination_removal_rank="${contamination_removal_rank:-domain}"' in transcriptome_core
    assert 'contamination_removal_target_taxon="${contamination_removal_target_taxon:-}"' in annotation_core
    assert 'contamination_removal_target_taxon="${contamination_removal_target_taxon:-}"' in transcriptome_core
    assert '--species_name "${contamination_removal_target_taxon:-${sp_ub}}"' in annotation_core
    assert '--species_name "${contamination_removal_target_taxon:-${sp_ub}}"' in transcriptome_core
    assert "GG_COMMON_CONTAMINATION_REMOVAL_RANK" not in annotation_core
    assert "GG_COMMON_CONTAMINATION_REMOVAL_RANK" not in transcriptome_core
    assert "contamination_removal_rank" in config_vars
    assert "contamination_removal_target_taxon" in config_vars


def test_transcriptome_entrypoint_uses_descriptive_busco_flag_names():
    entrypoint = _read_text(WORKFLOW_DIR / "gg_transcriptome_generation_entrypoint.sh")
    core = _read_text(CORE_DIR / "gg_transcriptome_generation_core.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")

    assert 'run_busco_isoforms=1 # BUSCO for transcriptome assembly with isoforms.' in entrypoint
    assert 'run_busco_longest_cds=1 # BUSCO for longest CDS.' in entrypoint
    assert 'run_busco_contamination_removed_longest_cds=0 # BUSCO for contamination-removed longest CDS.' in entrypoint
    assert 'disable_if_no_input_file "run_busco_isoforms" "${file_isoform}"' in core
    assert 'disable_if_no_input_file "run_busco_longest_cds" "${file_longestcds}"' in core
    assert 'disable_if_no_input_file "run_busco_contamination_removed_longest_cds" "${file_longestcds_contamination_removal_fasta}"' in core
    assert "run_busco1" not in entrypoint
    assert "run_busco2" not in entrypoint
    assert "run_busco3" not in entrypoint
    assert "run_busco1" not in core
    assert "run_busco2" not in core
    assert "run_busco3" not in core
    assert "run_busco1" not in config_vars
    assert "run_busco2" not in config_vars
    assert "run_busco3" not in config_vars


def test_transcriptome_wrapper_uses_amalgkit_default_filter_order():
    entrypoint = _read_text(WORKFLOW_DIR / "gg_transcriptome_generation_entrypoint.sh")
    core = _read_text(CORE_DIR / "gg_transcriptome_generation_core.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")

    assert "amalgkit_filter_order=" not in entrypoint
    assert "--filter_order" not in core
    assert "amalgkit_filter_order" not in config_vars


def test_entrypoint_modify_block_parameters_have_inline_comments():
    scripts = sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    assert scripts, "No entrypoint scripts were found."
    missing = []
    for script in scripts:
        for lineno, line in _entrypoint_modify_block_assignments(script):
            if "#" not in line:
                missing.append(f"{script}:{lineno}: {line}")
    assert not missing, "Add inline comments to parameter assignments:\n" + "\n".join(missing)


def test_common_parameters_have_inline_comments():
    script = WORKFLOW_DIR / "gg_common_params.sh"
    missing = [
        f"{script}:{lineno}: {line}"
        for lineno, line in _common_param_assignments(script)
        if "#" not in line
    ]
    assert not missing, "Add inline comments to common parameters:\n" + "\n".join(missing)


def test_genome_annotation_core_quotes_known_path_sensitive_options():
    script = CORE_DIR / "gg_genome_annotation_core.sh"
    text = _read_text(script)

    banned_tokens = [
        "--lineage_dataset ${dir_busco_lineage}",
        "--download_path ${dir_busco_db}",
        "seqkit seq --threads ${GG_TASK_CPUS} ${file_sp_genome} > \"busco_genome_input.fa\"",
        "mmseqs createdb ${file_sp_cds} queryDB",
        "--fasta_file ${file_sp_cds}",
        "--mmseqs2taxonomy_tsv ${file_sp_cds_mmseqs2taxonomy}",
        "--fx2tab_tsv ${file_sp_cds_fx2tab}",
        "--species_name ${sp_ub}",
        "--cds_fasta ${file_sp_cds}",
        "--uniprot_tsv ${file_sp_uniprot_annotation}",
        "--busco_tsv ${file_sp_cds_busco_full}",
        "--expression_tsv ${file_sp_expression}",
        "--gff_info ${file_sp_gff_info}",
        "--orthogroup_tsv ${file_orthogroup}",
        "--mmseqs2taxonomy_tsv ${file_sp_cds_mmseqs2taxonomy}",
        "--fx2tab ${file_sp_cds_fx2tab}",
        "--fx2tab ${file_sp_genome_fx2tab}",
        "mmseqs createdb ${file_sp_genome} queryDB",
        "--fasta_file ${file_sp_genome}",
        "--mmseqs2taxonomy_tsv ${file_sp_genome_mmseqs2taxonomy}",
        "--fx2tab_tsv ${file_sp_genome_fx2tab}",
        'grep -e "^${scaffold}[[:space:]]" tmp.species1.bed',
    ]
    for token in banned_tokens:
        assert token not in text, f"Found unquoted genome-annotation token: {token}"

    expected_tokens = [
        '--lineage_dataset "${dir_busco_lineage}"',
        '--download_path "${dir_busco_db}"',
        'seqkit seq --threads "${GG_TASK_CPUS}" "${file_sp_genome}" > "busco_genome_input.fa"',
        'mmseqs createdb "${file_sp_cds}" queryDB',
        '--fasta_file "${file_sp_cds}"',
        '--mmseqs2taxonomy_tsv "${file_sp_cds_mmseqs2taxonomy}"',
        '--fx2tab_tsv "${file_sp_cds_fx2tab}"',
        '--species_name "${contamination_removal_target_taxon:-${sp_ub}}"',
        '--cds_fasta "${file_sp_cds}"',
        '--uniprot_tsv "${file_sp_uniprot_annotation}"',
        '--busco_tsv "${file_sp_cds_busco_full}"',
        '--expression_tsv "${file_sp_expression}"',
        '--gff_info "${file_sp_gff_info}"',
        '--orthogroup_tsv "${file_orthogroup}"',
        '--mmseqs2taxonomy_tsv "${file_sp_cds_mmseqs2taxonomy}"',
        '--fx2tab "${file_sp_cds_fx2tab}"',
        '--fx2tab "${file_sp_genome_fx2tab}"',
        'mmseqs createdb "${file_sp_genome}" queryDB',
        '--fasta_file "${file_sp_genome}"',
        '--mmseqs2taxonomy_tsv "${file_sp_genome_mmseqs2taxonomy}"',
        '--fx2tab_tsv "${file_sp_genome_fx2tab}"',
        'awk -v scaffold="${scaffold}" \'$1 == scaffold\' tmp.species1.bed',
    ]
    for token in expected_tokens:
        assert token in text, f"Missing quoted genome-annotation token: {token}"


def test_genome_evolution_core_quotes_busco_lineage_and_trimal_input_paths():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)

    banned_tokens = [
        "--lineage_dataset ${dir_busco_lineage}",
        "--download_path ${dir_busco_db}",
        "seqkit seq --remove-gaps --threads 1 ${dir_single_copy_mafft}/${infile}",
        "seqkit translate --transl-table ${genetic_code} --threads 1 ${dir_single_copy_mafft}/${infile}",
        "seqkit seq --remove-gaps --threads 1 ${dir_busco_mafft}/${infile}",
        "seqkit translate --transl-table ${genetic_code} --threads 1 ${dir_busco_mafft}/${infile}",
        'trimal -in ${dir_single_copy_mafft}/${infile}',
        'trimal -in ${dir_busco_mafft}/${infile}',
    ]
    for token in banned_tokens:
        assert token not in text, f"Found unquoted genome-evolution token: {token}"

    expected_tokens = [
        '--lineage_dataset "${dir_busco_lineage}"',
        '--download_path "${dir_busco_db}"',
        'seqkit seq --remove-gaps --threads 1 "${dir_single_copy_mafft}/${infile}"',
        'seqkit translate --transl-table "${genetic_code}" --threads 1 "${dir_single_copy_mafft}/${infile}"',
        'seqkit seq --remove-gaps --threads 1 "${dir_busco_mafft}/${infile}"',
        'seqkit translate --transl-table "${genetic_code}" --threads 1 "${dir_busco_mafft}/${infile}"',
        'seqkit seq --threads 1 "${dir_single_copy_mafft}/${infile}" --out-file "tmp.${infile_base}.pep.aln.fasta"',
        'seqkit seq --threads 1 "${dir_busco_mafft}/${infile}" --out-file "tmp.${infile_base}.pep.aln.fasta"',
        '-in "tmp.${infile_base}.pep.aln.fasta"',
    ]
    for token in expected_tokens:
        assert token in text, f"Missing quoted genome-evolution token: {token}"


def test_gene_evolution_core_quotes_trait_promoter_and_summary_path_options():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)

    banned_tokens = [
        "--dir_trait ${dir_sp_expression}",
        "--dir_genome ${dir_sp_genome}",
        "--geneinfo_tsv ${file_og_gff_info}",
        "--species_tree ${species_tree_pruned}",
        "--unrooted_tree ${file_og_unrooted_tree_analysis}",
        "--rooted_tree ${file_og_rooted_tree_analysis}",
        "--rooting_log ${file_og_rooted_log}",
        "--dated_tree ${file_og_dated_tree_analysis}",
        "--dated_log ${file_og_dated_tree_log}",
        "--generax_nhx ${generax2orthogroup_statistics}",
        "--hyphy_dnds_json ${file_og_hyphy_dnds}",
        "--hyphy_relax_json ${file_og_hyphy_relax}",
        "--hyphy_relax_reversed_json ${file_og_hyphy_relax_reversed}",
        "--l1ou_tree ${file_og_l1ou_fit_tree}",
        "--l1ou_regime ${file_og_l1ou_fit_regime}",
        "--l1ou_leaf ${file_og_l1ou_fit_leaf}",
        "--expression ${file_og_expression}",
        "--mapdnds_tree_dn ${file_og_mapdnds_dn}",
        "--mapdnds_tree_ds ${file_og_mapdnds_ds}",
        "--codeml_tsv ${file_og_codeml_two_ratio}",
        "--character_gff ${file_og_gff_info}",
        "--fimo ${file_og_fimo}",
        "--scm_intron ${file_og_scm_intron_summary}",
        "--csubst_b ${file_og_csubst_b}",
        "--gene_pgls_stats ${file_og_gene_pgls}",
        "--species_pgls_stats ${file_og_species_pgls}",
        "--rpsblast ${file_og_rpsblast}",
        "--uniprot ${file_og_uniprot_annotation}",
        "--clade_ortholog_prefix ${treevis_clade_ortholog_prefix}",
        "run_hyphy_relax_for_all_traits 1 ${file_og_hyphy_relax}",
        "run_hyphy_relax_for_all_traits 0 ${file_og_hyphy_relax_reversed}",
        "get_hyphy_genetic_code ${genetic_code}",
    ]
    for token in banned_tokens:
        assert token not in text, f"Found unquoted gene-evolution token: {token}"

    expected_tokens = [
        '--dir_trait "${dir_sp_expression}"',
        '--dir_genome "${dir_sp_genome}"',
        '--geneinfo_tsv "${file_og_gff_info}"',
        '--species_tree "${species_tree_pruned}"',
        '--unrooted_tree "${file_og_unrooted_tree_analysis}"',
        '--rooted_tree "${file_og_rooted_tree_analysis}"',
        '--rooting_log "${file_og_rooted_log}"',
        '--dated_tree "${file_og_dated_tree_analysis}"',
        '--dated_log "${file_og_dated_tree_log}"',
        '--generax_nhx "${generax2orthogroup_statistics}"',
        '--hyphy_dnds_json "${file_og_hyphy_dnds}"',
        '--hyphy_relax_json "${file_og_hyphy_relax}"',
        '--hyphy_relax_reversed_json "${file_og_hyphy_relax_reversed}"',
        '--l1ou_tree "${file_og_l1ou_fit_tree}"',
        '--l1ou_regime "${file_og_l1ou_fit_regime}"',
        '--l1ou_leaf "${file_og_l1ou_fit_leaf}"',
        '--expression "${file_og_expression}"',
        '--mapdnds_tree_dn "${file_og_mapdnds_dn}"',
        '--mapdnds_tree_ds "${file_og_mapdnds_ds}"',
        '--codeml_tsv "${file_og_codeml_two_ratio}"',
        '--character_gff "${file_og_gff_info}"',
        '--fimo "${file_og_fimo}"',
        '--scm_intron "${file_og_scm_intron_summary}"',
        '--csubst_b "${file_og_csubst_b}"',
        '--gene_pgls_stats "${file_og_gene_pgls}"',
        '--species_pgls_stats "${file_og_species_pgls}"',
        '--rpsblast "${file_og_rpsblast}"',
        '--uniprot "${file_og_uniprot_annotation}"',
        '--clade_ortholog_prefix "${treevis_clade_ortholog_prefix}"',
        'run_hyphy_relax_for_all_traits 1 "${file_og_hyphy_relax}"',
        'run_hyphy_relax_for_all_traits 0 "${file_og_hyphy_relax_reversed}"',
        'get_hyphy_genetic_code "${genetic_code}"',
    ]
    for token in expected_tokens:
        assert token in text, f"Missing quoted gene-evolution token: {token}"


def test_gene_evolution_core_uses_exact_header_match_for_missing_query_gene_detection():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert 'if ! grep -q -e "^>${gene_name}" "${og_id}.query.cds.2.fasta"; then' not in text
    assert "if ! awk -v gene=\"${gene_name}\" '" in text
    assert "sub(/[[:space:]].*$/, \"\", header)" in text
    assert "if (header == gene)" in text


def test_gene_evolution_core_filters_empty_translated_records_before_diamond_makedb():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert "filter_translated_fasta_for_diamond() {" in text
    assert 'gsub(/\\*/, "", $0)' in text
    assert 'if ($0 != "") {' in text
    assert 'printf("Dropped %d translated protein records with empty sequence after stop-codon removal.\\n", dropped) > "/dev/stderr"' in text
    assert text.count('filter_translated_fasta_for_diamond \\') >= 2


def test_transcriptome_core_quotes_mmseqs_createdb_input_path():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    assert "mmseqs createdb ${file_longestcds} queryDB" not in text
    assert 'mmseqs createdb "${file_longestcds}" queryDB' in text


def test_gene_evolution_core_quotes_orthogroup_lookup_and_makeblastdb_args():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert 'og_id=$(python -c "import sys,pandas; df=pandas.read_csv(sys.argv[1],sep=\'\\t\',header=0); print(df.loc[int(sys.argv[2]),:].iloc[0])" ${file_orthogroup_genecount_selected} ${ind})' not in text
    assert 'og_id=$(python -c "import sys,pandas; df=pandas.read_csv(sys.argv[1],sep=\'\\t\',header=0); print(df.loc[int(sys.argv[2]),:].iloc[0])" "${file_orthogroup_genecount_selected}" "${ind}")' not in text
    assert "df=pandas.read_csv(sys.argv[1],sep='\\t',header=0); print(df.loc[int(sys.argv[2]),:].iloc[0])" not in text
    assert 'og_id=$(awk -F\'\\t\' -v row="${GG_ARRAY_TASK_ID}" \'NR == (row + 1) { print $1; exit }\' "${file_orthogroup_genecount_selected}")' in text
    assert 'makeblastdb -dbtype nucl -title ${sp_cds} -out ${sp_cds_blastdb}' not in text
    assert 'makeblastdb -dbtype nucl -in ${sp_cds} -out ${sp_cds_blastdb}' not in text
    assert 'makeblastdb -dbtype nucl -title "${sp_cds}" -out "${sp_cds_blastdb}"' in text
    assert 'makeblastdb -dbtype nucl -in "${sp_cds}" -out "${sp_cds_blastdb}"' in text


def test_transcriptome_core_avoids_direct_mv_out_glob_for_getfastq_and_quant():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    assert 'mv_out "${dir_tmp}"/getfastq/* "${dir_amalgkit_getfastq_sp}"' not in text
    assert "mv_out ./quant/* \"${dir_amalgkit_quant}/${sp_ub}\"" not in text
    assert 'getfastq_outputs=("${dir_tmp}"/getfastq/*)' in text
    assert 'mv_out "${getfastq_outputs[@]}" "${dir_amalgkit_getfastq_sp}"' in text
    assert "quant_outputs=(./quant/*)" in text
    assert 'mv_out "${quant_outputs[@]}" "${dir_amalgkit_quant}/${sp_ub}"' in text


def test_no_cp_out_or_mv_out_glob_arguments_in_core_scripts():
    pattern = re.compile(r"^[ \t]*(cp_out|mv_out)\b[^\n]*\*", re.MULTILINE)
    for script in sorted(CORE_DIR.glob("*.sh")):
        text = _read_text(script)
        assert pattern.search(text) is None, (
            f"Use nullglob+array guard instead of cp_out/mv_out glob in {script}"
        )


def test_gene_evolution_core_quotes_notung_zip_and_summary_presence_checks():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert "if [[ -s ${file_og_notung_reconcil} ]]; then" not in text
    assert "unzip -qf ${file_og_notung_reconcil}" not in text
    assert "! -s ${file_og_stat_branch}" not in text
    assert "! -s ${file_og_stat_tree}" not in text
    assert 'if [[ -s "${file_og_notung_reconcil}" ]]; then' in text
    assert 'unzip -qf "${file_og_notung_reconcil}"' in text
    assert '! -s "${file_og_stat_branch}"' in text
    assert '! -s "${file_og_stat_tree}"' in text


def test_genome_evolution_core_quotes_notung_unzip_and_rooting_temp_paths():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "unzip -q ${infile}",
        "2>&1 | tee ${busco_id}.root.txt",
        "if [[ -s ${busco_id}.root.nwk ]]; then",
        "run_mafft ${input_alignment_file} &",
        "run_trimal ${input_alignment_file} &",
        "--seqfile tmp.${infile_base}.cds.fasta",
        "--aa_aln tmp.${infile_base}.pep.aln.fasta",
        "--outfile tmp.${infile_base}.cds.aln.fasta",
        "if [[ -s tmp.${infile_base}.cds.aln.fasta ]]; then",
        "seqkit seq --threads 1 tmp.${infile_base}.cds.aln.fasta --out-file",
        "rm -f -- tmp.${infile_base}*",
        "if [[ -s tmp.${infile_base}.trimal.fasta ]]; then",
        "seqkit seq --threads 1 tmp.${infile_base}.trimal.fasta --out-file",
        "rm -f -- tmp.${infile_base}.*",
    ]
    for token in banned_tokens:
        assert token not in text, f"Found unquoted genome-evolution temp/rooting token: {token}"

    expected_tokens = [
        'unzip -q "${infile}"',
        '2>&1 | tee "${busco_id}.root.txt"',
        'if [[ -s "${busco_id}.root.nwk" ]]; then',
        'run_mafft "${input_alignment_file}" &',
        'run_trimal "${input_alignment_file}" &',
        '--seqfile "tmp.${infile_base}.cds.fasta"',
        '--aa_aln "tmp.${infile_base}.pep.aln.fasta"',
        '--outfile "tmp.${infile_base}.cds.aln.fasta"',
        'if [[ -s "tmp.${infile_base}.cds.aln.fasta" ]]; then',
        'seqkit seq --threads 1 "tmp.${infile_base}.cds.aln.fasta" --out-file',
        'rm -f -- "tmp.${infile_base}"*',
        'if [[ -s "tmp.${infile_base}.trimal.fasta" ]]; then',
        'seqkit seq --threads 1 "tmp.${infile_base}.trimal.fasta" --out-file',
        'rm -f -- "tmp.${infile_base}."*',
    ]
    for token in expected_tokens:
        assert token in text, f"Missing quoted genome-evolution temp/rooting token: {token}"


def test_gene_and_genome_core_quote_model_and_codon_options():
    targets = [
        CORE_DIR / "gg_gene_evolution_core.sh",
        CORE_DIR / "gg_genome_evolution_core.sh",
    ]
    patterns = [
        re.compile(r"--transl-table\s+\$\{[^}]+\}"),
        re.compile(r"--codontable\s+\$\{[^}]+\}"),
        re.compile(r"(^|\s)-m\s+\$\{[^}]+\}", re.MULTILINE),
        re.compile(r"--prefix\s+tmp\.\$\{[^}]+\}"),
    ]
    for script in targets:
        text = _read_text(script)
        for pattern in patterns:
            assert pattern.search(text) is None, (
                f"Found unquoted model/codon option in {script}: {pattern.pattern}"
            )


def test_gene_evolution_core_quotes_trimal_tmp_paths():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert "-out ${og_id}.cds.trimal.tmp1.fasta" not in text
    assert "--seqfile ${og_id}.cds.trimal.tmp1.fasta" not in text
    assert '-out "${og_id}.cds.trimal.tmp1.fasta"' in text
    assert '--seqfile "${og_id}.cds.trimal.tmp1.fasta"' in text


def test_genome_evolution_core_quotes_parallel_function_call_args():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "run_iqtree_pep ${input_alignment_file} &",
        "run_iqtree_dna ${input_alignment_file} &",
        "busco_iqtree_dna ${input_alignment_file} ${dir_busco_trimal} ${dir_busco_iqtree_dna} &",
        "busco_iqtree_pep ${input_alignment_file} ${dir_busco_trimal} ${dir_busco_iqtree_pep} &",
        "busco_notung ${infile} ${dir_busco_iqtree_dna} \"${dir_busco_notung_dna}\" &",
        "busco_notung ${infile} ${dir_busco_iqtree_pep} \"${dir_busco_notung_pep}\" &",
        "busco_species_tree_assisted_gene_tree_rooting ${infile} \"${dir_busco_notung_dna}\" ${dir_busco_iqtree_dna} \"${dir_busco_rooted_txt_dna}\" \"${dir_busco_rooted_nwk_dna}\" &",
        "busco_species_tree_assisted_gene_tree_rooting ${infile} \"${dir_busco_notung_pep}\" ${dir_busco_iqtree_pep} \"${dir_busco_rooted_txt_pep}\" \"${dir_busco_rooted_nwk_pep}\" &",
    ]
    for token in banned_tokens:
        assert token not in text, f"Found unquoted parallel call args: {token}"

    expected_tokens = [
        'run_iqtree_pep "${input_alignment_file}" &',
        'run_iqtree_dna "${input_alignment_file}" &',
        'busco_iqtree_dna "${input_alignment_file}" "${dir_busco_trimal}" "${dir_busco_iqtree_dna}" &',
        'busco_iqtree_pep "${input_alignment_file}" "${dir_busco_trimal}" "${dir_busco_iqtree_pep}" &',
        'busco_notung "${infile}" "${dir_busco_iqtree_dna}" "${dir_busco_notung_dna}" &',
        'busco_notung "${infile}" "${dir_busco_iqtree_pep}" "${dir_busco_notung_pep}" &',
        'busco_species_tree_assisted_gene_tree_rooting "${infile}" "${dir_busco_notung_dna}" "${dir_busco_iqtree_dna}" "${dir_busco_rooted_txt_dna}" "${dir_busco_rooted_nwk_dna}" &',
        'busco_species_tree_assisted_gene_tree_rooting "${infile}" "${dir_busco_notung_pep}" "${dir_busco_iqtree_pep}" "${dir_busco_rooted_txt_pep}" "${dir_busco_rooted_nwk_pep}" &',
    ]
    for token in expected_tokens:
        assert token in text, f"Missing quoted parallel call args token: {token}"


def test_genome_evolution_core_uses_array_for_optional_orthofinder_species_tree_arg():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "param_species_tree=''",
        "param_species_tree=\"-s ${species_tree}\"",
        "if [[ -n \"${param_species_tree}\" ]]; then",
        "${param_species_tree}; then",
    ]
    for token in banned_tokens:
        assert token not in text, f"Found fragile optional species-tree arg token: {token}"

    expected_tokens = [
        "param_species_tree=()",
        'param_species_tree=(-s "${species_tree}")',
        "if [[ ${#param_species_tree[@]} -gt 0 ]]; then",
        '"${param_species_tree[@]}"; then',
    ]
    for token in expected_tokens:
        assert token in text, f"Missing robust optional species-tree arg token: {token}"


def test_gene_evolution_core_quotes_notung_and_mapdnds_args():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "--prefix ${og_id} \\",
        "if [[ ! -s ${species_tree_pruned} ]]; then",
        "java -jar -Xmx${memory_notung}g ${notung_jar} \\",
        "-s ${species_tree_pruned} \\",
        "-g ${file_og_unrooted_tree_analysis} \\",
        "-g ${og_id}.root.nwk \\",
        "--seqtype CODON${genetic_code} \\",
        "--prefix ${og_id}.iqtree2mapdNdS \\",
        "--iqtree ${og_id}.iqtree2mapdNdS.iqtree \\",
        "--log ${og_id}.iqtree2mapdNdS.log \\",
        "--state ${og_id}.iqtree2mapdNdS.state \\",
        "--treefile ${og_id}.iqtree2mapdNdS.treefile \\",
        "--genetic_code ${genetic_code}",
    ]
    for token in banned_tokens:
        assert token not in text, f"Found unquoted notung/mapdNdS token: {token}"

    expected_tokens = [
        '--prefix "${og_id}" \\',
        'if [[ ! -s "${species_tree_pruned}" ]]; then',
        'java -jar -Xmx${memory_notung}g "${notung_jar}" \\',
        '-s "${species_tree_pruned}" \\',
        '-g "${file_og_unrooted_tree_analysis}" \\',
        '-g "${og_id}.root.nwk" \\',
        '--seqtype "CODON${genetic_code}" \\',
        '--prefix "${og_id}.iqtree2mapdNdS" \\',
        '--iqtree "${og_id}.iqtree2mapdNdS.iqtree" \\',
        '--log "${og_id}.iqtree2mapdNdS.log" \\',
        '--state "${og_id}.iqtree2mapdNdS.state" \\',
        '--treefile "${og_id}.iqtree2mapdNdS.treefile" \\',
        '--genetic_code "${genetic_code}"',
    ]
    for token in expected_tokens:
        assert token in text, f"Missing quoted notung/mapdNdS token: {token}"


def test_genome_evolution_core_quotes_orthofinder_cleanup_calls():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    assert "orthofinder_output_directory_cleanup ${dir_orthofinder} ${GG_TASK_CPUS}" not in text
    assert "orthofinder_output_directory_cleanup ${dir_orthofinder}/core ${GG_TASK_CPUS}" not in text
    assert 'orthofinder_output_directory_cleanup "${dir_orthofinder}" "${GG_TASK_CPUS}"' in text
    assert 'orthofinder_output_directory_cleanup "${dir_orthofinder}/core" "${GG_TASK_CPUS}"' in text


def test_no_line_start_option_uses_unquoted_variable_in_core_scripts():
    long_option = re.compile(r"^[ \t]*--[A-Za-z0-9_-]+\s+\$\{[^}]+\}", re.MULTILINE)
    short_option = re.compile(r"^[ \t]*-[A-Za-z]\s+\$\{[^}]+\}", re.MULTILINE)
    single_dash_long_option = re.compile(r"^[ \t]*-[A-Za-z0-9][A-Za-z0-9_-]+\s+\$\{[^}]+\}", re.MULTILINE)
    for script in sorted(CORE_DIR.glob("*.sh")):
        text = _read_text(script)
        assert long_option.search(text) is None, (
            f"Found unquoted variable in long option argument in {script}"
        )
        assert short_option.search(text) is None, (
            f"Found unquoted variable in short option argument in {script}"
        )
        assert single_dash_long_option.search(text) is None, (
            f"Found unquoted variable in single-dash long option argument in {script}"
        )


def test_wait_for_background_jobs_helper_exists():
    util_script = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_script)
    assert re.search(r"^\s*wait_for_background_jobs\(\)\s*\{", text, re.MULTILINE)
    body = _function_body(text, "wait_for_background_jobs")
    assert "mapfile -t pids" not in body
    assert "done < <(jobs -pr)" in body


def test_gene_and_genome_core_do_not_use_plain_wait():
    targets = [
        CORE_DIR / "gg_gene_evolution_core.sh",
        CORE_DIR / "gg_genome_evolution_core.sh",
    ]
    plain_wait = re.compile(r"^\s*wait\s*$", re.MULTILINE)
    for script in targets:
        text = _read_text(script)
        assert plain_wait.search(text) is None, f"Found plain wait in {script}"
        assert "wait_for_background_jobs" in text, f"Missing wait helper usage in {script}"


def test_gene_evolution_core_uses_nullglob_array_for_query_hits_merge():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert "for query_hits_tmp_file in ${query_hits_tmp_dir}/*.hits.fasta; do" not in text
    assert 'query_hits_tmp_files=("${query_hits_tmp_dir}"/*.hits.fasta)' in text
    assert 'for query_hits_tmp_file in "${query_hits_tmp_files[@]}"; do' in text


def test_genome_annotation_core_skips_identity_rename_in_jcvi_output_loop():
    script = CORE_DIR / "gg_genome_annotation_core.sh"
    text = _read_text(script)
    assert 'mv_out "${file}" "${file/species1.species2/${sp_ub}}"' not in text
    assert 'renamed_file="${file/species1.species2/${sp_ub}}"' in text
    assert 'if [[ "${renamed_file}" != "${file}" ]]; then' in text


def test_no_for_seq_command_substitution_in_core_scripts():
    pattern = re.compile(r"^\s*for\s+[A-Za-z_][A-Za-z0-9_]*\s+in\s+\$\(seq\b", re.MULTILINE)
    for script in sorted(CORE_DIR.glob("*.sh")):
        text = _read_text(script)
        assert pattern.search(text) is None, (
            f"Use arithmetic for-loop instead of for-in $(seq ...) in {script}"
        )


def test_transcriptome_core_uses_array_args_for_trinity_and_rnaspades_inputs():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    assert "trinity_input=" not in text
    assert "${trinity_input}" not in text
    assert "--single \"${in_single}\"" in text
    assert "--left \"${in_left}\"" in text
    assert "--right \"${in_right}\"" in text

    assert "rnaspades_input=$(for i in" not in text
    assert "rnaspades_input_args=()" in text
    assert '"${rnaspades_input_args[@]}"' in text


def test_transcriptome_core_uses_array_for_assembly_stat_input_files():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    assert "input_files=${file_isoform}" not in text
    assert 'input_files="${input_files} ${file_longestcds}"' not in text
    assert 'input_files="${input_files} ${file_longestcds_contamination_removal_fasta}"' not in text
    assert "${input_files}" not in text
    assert 'input_files=("${file_isoform}")' in text
    assert 'input_files+=("${file_longestcds}")' in text
    assert 'input_files+=("${file_longestcds_contamination_removal_fasta}")' in text
    assert '"${input_files[@]}"' in text


def test_transcriptome_core_quotes_get_total_fastq_len_dir_argument():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    banned_tokens = [
        'get_total_fastq_len ${selected_fastq_dir} "*.amalgkit.fastq.gz"',
        'get_total_fastq_len ${selected_fastq_dir} "*_1.amalgkit.fastq.gz"',
        'get_total_fastq_len ${selected_fastq_dir} "*_2.amalgkit.fastq.gz"',
        'get_total_fastq_len ${assembly_input_fastq_dir} "*.amalgkit.fastq.gz"',
    ]
    for token in banned_tokens:
        assert token not in text, f"Found unquoted get_total_fastq_len dir arg token: {token}"

    expected_tokens = [
        'get_total_fastq_len "${selected_fastq_dir}" "*.amalgkit.fastq.gz"',
        'get_total_fastq_len "${selected_fastq_dir}" "*_1.amalgkit.fastq.gz"',
        'get_total_fastq_len "${selected_fastq_dir}" "*_2.amalgkit.fastq.gz"',
        'get_total_fastq_len "${assembly_input_fastq_dir}" "*.amalgkit.fastq.gz"',
    ]
    for token in expected_tokens:
        assert token in text, f"Missing quoted get_total_fastq_len dir arg token: {token}"


def test_transcriptome_core_guards_non_positive_assembly_resources():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    expected_tokens = [
        "if [[ ${assembly_cpus} -lt 1 ]]; then",
        "assembly_cpus=1",
        "if [[ ${assembly_mem_gb} -lt 1 ]]; then",
        "assembly_mem_gb=1",
    ]
    for token in expected_tokens:
        assert token in text, f"Missing assembly resource guard token: {token}"


def test_genome_evolution_core_uses_array_args_for_nwkit_mcmctree_constraints():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        'bound_params="${bound_params} --lower_bound ${mcmctree_params[2]}"',
        'bound_params="${bound_params} --upper_bound ${mcmctree_params[3]}"',
        'left_right="--left_species ${mcmctree_params[0]} --right_species ${mcmctree_params[1]}"',
        'tree_string=$(printf \'%s\\n\' "${tree_string}" | nwkit mcmctree ${left_right} ${bound_params})',
        'echo -e "${tree_string}" > "tmp.constrained.tree.nwk"',
        'nwkit mcmctree \\n    --infile "${file_undated_species_tree}"',
    ]
    for token in banned_tokens:
        assert token not in text, f"Found string-concatenated mcmctree args token: {token}"

    expected_tokens = [
        'dir_nwkit_download_dir="${gg_workspace_downloads_dir}/nwkit_downloads"',
        'ensure_dir "${dir_nwkit_download_dir}"',
        'nwkit mcmctree \\',
        '--download_dir "${dir_nwkit_download_dir}"',
        "nwkit_args=(",
        '--download_dir "${dir_nwkit_download_dir}"',
        '--left_species "${mcmctree_params[0]}"',
        '--right_species "${mcmctree_params[1]}"',
        'nwkit_args+=(--lower_bound "${mcmctree_params[2]}")',
        'nwkit_args+=(--upper_bound "${mcmctree_params[3]}")',
        'tree_string=$(printf \'%s\\n\' "${tree_string}" | nwkit mcmctree "${nwkit_args[@]}")',
        'printf \'%s\\n\' "${tree_string}" > "tmp.constrained.tree.nwk"',
    ]
    for token in expected_tokens:
        assert token in text, f"Missing array-based mcmctree args token: {token}"


def test_genome_evolution_core_initializes_concat_iqtree_optional_args_as_arrays():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "bootstrap_params=''",
        "iqtree_mem_arg=''",
        "iqtree_mem_arg=\"-mem ${GG_MEM_TOTAL_GB}G\"",
        "${iqtree_mem_arg} \\",
        "${bootstrap_params}; then",
    ]
    for token in banned_tokens:
        assert token not in text, f"Found fragile concat-IQ-TREE args token: {token}"

    expected_tokens = [
        "bootstrap_params=(--ufboot 1000 --bnni)",
        "bootstrap_params=()",
        "iqtree_mem_args=()",
        'iqtree_mem_args=(-mem "${GG_MEM_TOTAL_GB}G")',
        '"${iqtree_mem_args[@]}" \\',
        '"${bootstrap_params[@]}"; then',
    ]
    for token in expected_tokens:
        assert token in text, f"Missing array-based concat-IQ-TREE args token: {token}"


def test_gene_evolution_core_quotes_nwkit_subtree_leaves_argument():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert "--leaves ${comma_separated_genes}" not in text
    assert '--leaves "${comma_separated_genes}"' in text
    assert "run_nwkit_subtree ${subtree_infile}" not in text
    assert 'run_nwkit_subtree "${subtree_infile}"' in text


def test_gene_evolution_core_uses_array_optional_args_for_iqtree_and_csubst():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        'other_iqtree_params="--ufboot 1000 --bnni"',
        'other_iqtree_params="${other_iqtree_params} --fast"',
        "${other_iqtree_params}",
        'foreground_params="--foreground foreground.tsv --fg_format 2"',
        'foreground_params=""',
        "${foreground_params}",
    ]
    for token in banned_tokens:
        assert token not in text, f"Found string-based optional args token: {token}"

    expected_tokens = [
        "other_iqtree_params=()",
        "other_iqtree_params=(--ufboot 1000 --bnni)",
        "other_iqtree_params+=(--fast)",
        '"${other_iqtree_params[@]}"',
        "foreground_params=(--foreground foreground.tsv --fg_format 2)",
        "foreground_params=()",
        '"${foreground_params[@]}"',
    ]
    for token in expected_tokens:
        assert token in text, f"Missing array-based optional args token: {token}"


def test_gene_evolution_core_disables_initial_ufboot_when_fast_mode_is_enabled():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert 'if [[ ${num_seq} -gt ${iqtree_fast_mode_gt} ]]; then' in text
    assert 'if [[ ${use_ufboot} -eq 1 ]]; then' in text
    assert 'Disabling IQ-TREE UFBOOT because fast mode is enabled for large alignments (${num_seq} > ${iqtree_fast_mode_gt}).' in text
    assert 'other_iqtree_params=()' in text
    assert 'file_tree="${og_id}.treefile"' in text
    assert 'other_iqtree_params+=(--fast)' in text


def test_gene_evolution_core_keeps_generax_ufboot_task_free_of_fast_flag():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert 'Skipping IQ-TREE --fast in UFBOOT-on-GeneRax mode because the options are incompatible.' in text

    ufboot_block_start = text.index('task="IQ-TREE UFBOOT on GeneRax topology"')
    ufboot_block_end = text.index('build_iqtree_mem_args', ufboot_block_start)
    ufboot_block = text[ufboot_block_start:ufboot_block_end]
    assert 'other_iqtree_params+=( --fast )' not in ufboot_block


def test_no_pipe_to_grep_q_in_core_and_support_scripts():
    scripts = sorted(CORE_DIR.glob("*.sh")) + sorted((WORKFLOW_DIR / "support").glob("*.sh"))
    pattern = re.compile(r"\|\s*(?:z?grep)\s+-q|\|\s*grep\s+-Fq|\|\s*grep\s+-Fxq")
    for script in scripts:
        for line_no, line in enumerate(_read_text(script).splitlines(), start=1):
            stripped = line.lstrip()
            if stripped.startswith("#"):
                continue
            assert pattern.search(line) is None, (
                f"Use non-pipeline grep checks under pipefail in {script}:{line_no}: {line}"
            )


def test_no_pipe_to_head_n1_in_core_and_support_scripts():
    scripts = sorted(CORE_DIR.glob("*.sh")) + sorted((WORKFLOW_DIR / "support").glob("*.sh"))
    pattern = re.compile(r"\|\s*head\s+-n\s*1|\|\s*head\s+-1")
    for script in scripts:
        for line_no, line in enumerate(_read_text(script).splitlines(), start=1):
            stripped = line.lstrip()
            if stripped.startswith("#"):
                continue
            assert pattern.search(line) is None, (
                f"Avoid pipe-to-head under pipefail in {script}:{line_no}: {line}"
            )


def test_no_pipe_to_awk_exit_in_core_and_support_scripts():
    scripts = sorted(CORE_DIR.glob("*.sh")) + sorted((WORKFLOW_DIR / "support").glob("*.sh"))
    pattern = re.compile(r"\|\s*awk\b[^\n]*\bexit\b")
    for script in scripts:
        for line_no, line in enumerate(_read_text(script).splitlines(), start=1):
            stripped = line.lstrip()
            if stripped.startswith("#"):
                continue
            assert pattern.search(line) is None, (
                f"Avoid pipe-to-awk-exit under pipefail in {script}:{line_no}: {line}"
            )


def test_no_standalone_unquoted_file_or_dir_variable_argument_lines_in_core_scripts():
    pattern = re.compile(r"^[ \t]*\$\{(?:file|dir)_[^}]+\}[ \t]*\\?$", re.MULTILINE)
    for script in sorted(CORE_DIR.glob("*.sh")):
        text = _read_text(script)
        assert pattern.search(text) is None, (
            f"Found standalone unquoted file/dir variable argument line in {script}"
        )


def test_core_scripts_quote_is_output_older_than_inputs_target_path():
    targets = [
        CORE_DIR / "gg_gene_evolution_core.sh",
        CORE_DIR / "gg_genome_annotation_core.sh",
        CORE_DIR / "gg_transcriptome_generation_core.sh",
    ]
    for script in targets:
        text = _read_text(script)
        assert re.search(r'is_output_older_than_inputs\s+"\^file_[^"]+"\s+\$\{[^}]+\}', text) is None, (
            f"Found unquoted is_output_older_than_inputs path arg in {script}"
        )


def test_gene_evolution_core_uses_exact_qacc_match_for_rpsblast_hit_presence():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert 'if ! grep -F -q -- "${gene}" "${og_id}.rpsblast.tmp.tsv"; then' not in text
    assert "if ! awk -F '\\t' -v gene=\"${gene}\" '$1 == gene {found=1; exit} END {exit(found ? 0 : 1)}' \"${og_id}.rpsblast.tmp.tsv\"; then" in text
    assert "qlen=$(seqkit fx2tab --length ungapped_translated_cds.fas | awk -F '\\t' -v gene=\"${gene}\" '$1 == gene {print $NF}')" in text


def test_gene_evolution_core_tree_pruning_uses_awk_id_match_instead_of_python_split_filter():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "python -c \"import sys,re; keys = open(sys.argv[1]).read().split('\\n'); entries = open(sys.argv[2]).read().split('>'); [ sys.stdout.write('>'+e) for e in entries if (re.sub('\\n.*','',e) in keys)&(len(e)!=0) ]\"",
    ]
    for token in banned_tokens:
        assert token not in text
    assert "sub(/[[:space:]].*$/, \"\", id)" in text
    assert "keep_seq = (id in keep)" in text
    assert "' target_genes.txt \"${og_id}.untrimmed.input.fasta\" > \"${og_id}.untrimmed.pruned.tmp.fasta\"" in text
    assert "' target_genes.txt \"${og_id}.trimmed.input.fasta\" > \"${og_id}.trimmed.pruned.tmp.fasta\"" in text


def test_gg_util_python_package_check_uses_importlib_with_argv():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "gg_test_python_packages")
    assert 'python -c "import ${pypackage}"' not in body
    assert "importlib.import_module(sys.argv[1])" in body


def test_orthogroup_statistics_flushes_terminal_contiguous_n_run():
    script = WORKFLOW_DIR / "support" / "orthogroup_statistics.py"
    text = _read_text(script)
    assert "# Flush the final contiguous N-run when it reaches the sequence end." in text
    assert "new_slice = N_slices[len(N_slices)-1]+':'+N_extension_until" in text
    assert "re.finditer('[Nn]', seqs[i])" in text
    assert "if len(seqs) == 0 or all((s.strip() == '') for s in seqs):" in text


def test_orthogroup_statistics_handles_ete4_root_support_assertion_on_reroot():
    script = WORKFLOW_DIR / "support" / "orthogroup_statistics.py"
    text = _read_text(script)
    assert "root has branch property: support" in text
    assert "clear_root_branch_property_compat(tree, 'support')" in text


def test_orthogroup_statistics_skips_invalid_regime2tree_summary_instead_of_aborting():
    script = WORKFLOW_DIR / "support" / "orthogroup_statistics.py"
    text = _read_text(script)
    assert "tree_tmp = regime2tree(params[method+'_regime'])" in text
    assert "except ValueError as exc:" in text
    assert "Skipping {} regime summary due to invalid regime parameters" in text


def test_support_python_scalar_conditions_use_logical_and_not_bitwise_and():
    orthogroup_stats = _read_text(WORKFLOW_DIR / "support" / "orthogroup_statistics.py")
    csubst_wrapper = _read_text(WORKFLOW_DIR / "support" / "csubst_site_wrapper.py")

    orthogroup_banned = [
        '(os.path.exists(params["unrooted_tree"]))&(os.path.exists(params["rooted_tree"]))',
        '(os.path.exists(params["dated_tree"]))&(os.path.exists(params["species_tree"]))',
        '(os.path.exists(params["l1ou_regime"]))&(os.path.exists(params["l1ou_leaf"]))',
        '(os.path.exists(params["expression"]))&(os.path.exists(params["rooted_tree"]))',
        "(params['clade_ortholog_prefix']!='')&(os.path.exists(params[\"rooted_tree\"]))",
    ]
    for token in orthogroup_banned:
        assert token not in orthogroup_stats
    assert 'os.path.exists(params["unrooted_tree"]) and os.path.exists(params["rooted_tree"])' in orthogroup_stats
    assert 'os.path.exists(params["dated_tree"]) and os.path.exists(params["species_tree"])' in orthogroup_stats

    csubst_banned = [
        "(og_indices.shape[0] > args.max_per_og)&(args.max_per_og > 0)",
        "(not os.path.exists(dir_out)) & (not os.path.exists(out_zip))",
        "os.path.exists(dir_out) & (not os.path.exists(out_zip))",
    ]
    for token in csubst_banned:
        assert token not in csubst_wrapper
    assert "(og_indices.shape[0] > args.max_per_og) and (args.max_per_og > 0)" in csubst_wrapper
    assert "(not os.path.exists(dir_out)) and (not os.path.exists(out_zip))" in csubst_wrapper
    assert "os.path.exists(dir_out) and (not os.path.exists(out_zip))" in csubst_wrapper


def test_gene_evolution_core_uses_csubst_search_namespace():
    core = _read_text(CORE_DIR / "gg_gene_evolution_core.sh")
    assert "csubst analyze \\" not in core
    assert 'csubst search \\' in core
    assert 'csubst_search_dir="csubst_search"' in core
    assert '"${csubst_search_dir}/csubst_cb_stats.tsv"' in core
    for redundant_flag in [
        '--infile_type "iqtree"',
        '--iqtree_redo "no"',
        '--mg_sister "no"',
        '--exclude_sister_pair "yes"',
        '--ml_anc "no"',
        '--b "yes"',
        '--s "no"',
        '--cs "no"',
        '--cb "yes"',
        '--bs "no"',
        '--cbs "no"',
        '--asrv "each"',
        '--calibrate_longtail "yes"',
        '--outdir "${csubst_search_dir}"',
    ]:
        assert redundant_flag not in core


def test_csubst_site_wrapper_omits_redundant_sites_defaults():
    wrapper = _read_text(WORKFLOW_DIR / "support" / "csubst_site_wrapper.py")
    assert "cmd = ['csubst', 'sites']" in wrapper
    assert "cmd += ['--ml_anc', 'no']" not in wrapper
    assert "cmd += ['--mafft_exe', 'mafft']" not in wrapper


def test_is_fastq_requiring_downstream_analysis_done_quotes_path_checks():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "is_fastq_requiring_downstream_analysis_done")
    assert '-s ${file_isoform}' not in body
    assert '-s ${file_amalgkit_merge_count}' not in body
    assert '-s "${file_isoform}"' in body
    assert '-s "${file_amalgkit_merge_count}"' in body
    assert 'echo "${out}"' in body


def test_genome_annotation_core_has_jcvi_scaffold_fallback_for_small_genomes():
    script = CORE_DIR / "gg_genome_annotation_core.sh"
    text = _read_text(script)
    assert "No scaffolds >= ${minimum_scaffold_size} bp were found. Falling back to top 20 longest scaffolds." in text
    assert "sort_values('length',ascending=False)" in text
    assert "No scaffolds were selected from ${file_sp_genome_fx2tab}. Exiting." in text


def test_gg_util_direct_cp_mv_calls_use_option_separator():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    banned_tokens = [
        'mv "${pep_tmp}" "${output_prefix}.pep"',
        'mv "${dmnd_tmp_prefix}.dmnd" "${output_prefix}.dmnd"',
        'cp "${f}" "${staged_dir}/"',
        'mv "${staged_dir}" "${output_dir}.tmp"',
        'mv "${output_dir}.tmp" "${output_dir}"',
        'mv "${tmp_file}" "${output_file}"',
        'mv "${latest_marker}.tmp" "${latest_marker}"',
    ]
    for token in banned_tokens:
        assert token not in text
    expected_tokens = [
        'mv -- "${pep_tmp}" "${output_prefix}.pep"',
        'mv -- "${dmnd_tmp_prefix}.dmnd" "${output_prefix}.dmnd"',
        'cp -- "${f}" "${staged_dir}/"',
        'mv -- "${staged_dir}" "${output_dir}.tmp"',
        'mv -- "${output_dir}.tmp" "${output_dir}"',
        'mv -- "${tmp_file}" "${output_file}"',
        'mv -- "${latest_marker}.tmp" "${latest_marker}"',
    ]
    for token in expected_tokens:
        assert token in text


def test_support_shell_scripts_direct_cp_mv_require_option_separator():
    scripts = sorted((WORKFLOW_DIR / "support").glob("*.sh"))
    assert scripts, "No support shell scripts were found."
    pattern = re.compile(r"^[ \t]*(?:if[ \t]+)?(cp|mv)[ \t]+(?!--)", re.MULTILINE)
    for script in scripts:
        for line_no, line in enumerate(_read_text(script).splitlines(), start=1):
            stripped = line.lstrip()
            if stripped.startswith("#"):
                continue
            assert pattern.search(line) is None, (
                f"Use `cp --`/`mv --` in support shell scripts: {script}:{line_no}: {line}"
            )


def test_core_scripts_do_not_use_echo_triple_quote_blocks():
    pattern = re.compile(r'echo\s+"""')
    for script in sorted(CORE_DIR.glob("*.sh")):
        text = _read_text(script)
        assert pattern.search(text) is None, f"Use printf or heredoc instead of echo triple quotes: {script}"


def test_genome_annotation_core_uses_z_check_for_optional_expression_and_gff_paths():
    script = CORE_DIR / "gg_genome_annotation_core.sh"
    text = _read_text(script)
    assert "if [[ ! -n ${file_sp_expression} ]]; then" not in text
    assert "if [[ ! -n ${file_sp_gff} ]]; then" not in text
    assert 'if [[ -z "${file_sp_expression}" ]]; then' in text
    assert 'if [[ -z "${file_sp_gff}" ]]; then' in text


def test_gene_evolution_core_quotes_path_in_deactivation_messages():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "echo '${run_get_gff_info} is deactivated. Empty input:' ${dir_sp_gff}",
        "echo '${run_scm_intron} is deactivated. Empty input:' ${dir_sp_gff}",
        "echo '${dir_sp_expression} is not empty. Continued:' ${dir_sp_expression}",
        "echo '${dir_sp_expression} is empty:' ${dir_sp_expression}",
        "echo '${dir_sp_genome} is not empty. Continued:' ${dir_sp_genome}",
        "echo '${dir_sp_genome} is empty:' ${dir_sp_genome}",
    ]
    for token in banned_tokens:
        assert token not in text
    expected_tokens = [
        'echo "\\${run_get_gff_info} is deactivated. Empty input: ${dir_sp_gff}"',
        'echo "\\${run_scm_intron} is deactivated. Empty input: ${dir_sp_gff}"',
        'echo "\\${dir_sp_expression} is not empty. Continued: ${dir_sp_expression}"',
        'echo "\\${dir_sp_expression} is empty: ${dir_sp_expression}"',
        'echo "\\${dir_sp_genome} is not empty. Continued: ${dir_sp_genome}"',
        'echo "\\${dir_sp_genome} is empty: ${dir_sp_genome}"',
    ]
    for token in expected_tokens:
        assert token in text


def test_genome_evolution_core_quotes_shuf_count_for_orthofinder_core_sampling():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    assert 'shuf -n ${max_orthofinder_core_species}' not in text
    assert 'shuf -n "${max_orthofinder_core_species}"' in text


def test_transcriptome_gene_and_genome_evolution_core_guard_tmp_delete_against_root_path():
    transcriptome = _read_text(CORE_DIR / "gg_transcriptome_generation_core.sh")
    gene_evolution = _read_text(CORE_DIR / "gg_gene_evolution_core.sh")
    genome_evolution = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")

    assert 'if [[ -n "${dir_tmp:-}" && "${dir_tmp}" != "/" ]]; then' in transcriptome
    assert "Refusing to delete unsafe tmp directory" in transcriptome

    assert 'if [[ -n "${dir_tmp:-}" && -d "${dir_tmp}" && "${dir_tmp}" != "/" ]]; then' in gene_evolution
    assert "Refusing to delete unsafe tmp directory" in gene_evolution

    assert 'if [[ -d "${dir_tmp}" && -n "${dir_tmp:-}" && "${dir_tmp}" != "/" ]]; then' in genome_evolution
    assert "Refusing to delete unsafe tmp directory" in genome_evolution
    assert 'echo "$(date): end"' in genome_evolution


def test_transcriptome_core_busco_summary_loop_guards_missing_dir_before_find():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    assert 'if [[ -z "$(find "${dir_busco}" -mindepth 1 -print -quit)" ]]; then' not in text
    assert 'if [[ ! -d "${dir_busco}" || -z "$(find "${dir_busco}" -mindepth 1 -print -quit 2> /dev/null)" ]]; then' in text


def test_gene_evolution_core_quotes_key_s_checks_in_downstream_tasks():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "if [[ -s ${file_og_expression} && ${run_l1ou} -eq 1 ]]; then",
        "if [[ ! -s ${file_og_hyphy_relax_reversed} && ${run_hyphy_relax_reversed} -eq 1 ]]; then",
        "if [[ ! -s ${file_og_scm_intron_summary} && ${run_scm_intron} -eq 1 ]]; then",
        "if [[ ( ! -s ${file_og_l1ou_fit_rdata} || ! -s ${file_og_l1ou_fit_tree} || ! -s ${file_og_l1ou_fit_regime} || ! -s ${file_og_l1ou_fit_leaf} ) && ${run_l1ou} -eq 1 ]]; then",
        "if ( [[ ${summary_flag} -eq 1 || ! -s ${file_og_tree_plot} ]] ) && [[ ${run_tree_plot} -eq 1 ]]; then",
        "if [[ -s ${file_og_stat_branch} && -s ${file_og_stat_tree} && -s ${file_og_tree_plot} && ${gg_debug_mode:-0} -eq 0 ]]; then",
    ]
    for token in banned_tokens:
        assert token not in text

    expected_tokens = [
        'if [[ -s "${file_og_expression}" && ${run_l1ou} -eq 1 ]]; then',
        'if [[ ! -s "${file_og_hyphy_relax_reversed}" && ${run_hyphy_relax_reversed} -eq 1 ]]; then',
        'if [[ ! -s "${file_og_scm_intron_summary}" && ${run_scm_intron} -eq 1 ]]; then',
        'if [[ (! -s "${file_og_l1ou_fit_rdata}" || ! -s "${file_og_l1ou_fit_tree}" || ! -s "${file_og_l1ou_fit_regime}" || ! -s "${file_og_l1ou_fit_leaf}") && ${run_l1ou} -eq 1 ]]; then',
        'if ([[ ${summary_flag} -eq 1 || ! -s "${file_og_tree_plot}" ]]) && [[ ${run_tree_plot} -eq 1 ]]; then',
        'if [[ -s "${file_og_stat_branch}" && -s "${file_og_stat_tree}" && -s "${file_og_tree_plot}" && ${gg_debug_mode:-0} -eq 0 ]]; then',
    ]
    for token in expected_tokens:
        assert token in text

    assert "run_phylogeneticem" not in text


def test_gene_evolution_core_guards_array_task_id_range_before_input_indexing():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)

    assert 'mapfile -t files < <(find "${dir_genelist}" -mindepth 1 -maxdepth 1 -type f ! -name \'.*\' | sort)' in text
    assert 'files=( "${dir_genelist}"/* )' not in text
    assert 'if [[ ! "${GG_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${GG_ARRAY_TASK_ID} -lt 1 ]]; then' in text
    assert 'num_orthogroups=$(awk \'END { print (NR > 0 ? NR - 1 : 0) }\'' in text
    assert 'if [[ ${GG_ARRAY_TASK_ID} -gt ${num_orthogroups} ]]; then' in text
    assert "df=pandas.read_csv(sys.argv[1],sep='\\t',header=0); print(df.loc[int(sys.argv[2]),:].iloc[0])" not in text
    assert 'og_id=$(awk -F\'\\t\' -v row="${GG_ARRAY_TASK_ID}" \'NR == (row + 1) { print $1; exit }\'' in text

    idx_guard = 'if [[ ${idx} -ge ${#files[@]} ]]; then'
    idx_use = 'file_query_gene="${files[${idx}]}"'
    assert idx_guard in text
    assert idx_use in text
    assert text.index(idx_guard) < text.index(idx_use)


def test_transcriptome_core_guards_array_task_id_range_before_array_indexing():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)

    assert 'dir_input_fastq="${dir_transcriptome_assembly_input}/input_fastq"' not in text
    assert 'dir_input_sra_list="${dir_transcriptome_assembly_input}/input_sra_list"' not in text
    assert "Backward compatibility for legacy input layout" not in text
    assert 'mapfile -t fastq_mode_dirs < <(find "${dir_input_fastq}" -mindepth 1 -maxdepth 1 -type d ! -name \'.*\' | sort)' in text
    assert 'mapfile -t files_fastq < <(find "${dir_species_fastq}" -maxdepth 1 -type f ! -name \'.*\' \\( -name "*.fq" -o -name "*.fastq" -o -name "*.fq.gz" -o -name "*.fastq.gz" \\) | sort)' in text
    assert 'mapfile -t sra_mode_files < <(find "${dir_input_sra_list}" -mindepth 1 -maxdepth 1 -type f ! -name \'.*\' | sort)' in text
    assert 'mapfile -t metadata_mode_files < <(find "${dir_input_amalgkit_metadata}" -mindepth 1 -maxdepth 1 -type f ! -name \'.*\' | sort)' in text
    assert 'dirs=( "${fastq_mode_dirs[@]}" )' in text
    assert 'files=( "${sra_mode_files[@]}" )' in text
    assert 'files=( "${metadata_mode_files[@]}" )' in text
    assert re.search(r'^[ \t]*dirs=\( "\$\{dir_input_fastq\}"/\* \)', text, re.MULTILINE) is None
    assert re.search(r'^[ \t]*files=\( "\$\{dir_input_sra_list\}"/\* \)', text, re.MULTILINE) is None
    assert re.search(r'^[ \t]*files=\( "\$\{dir_input_amalgkit_metadata\}"/\* \)', text, re.MULTILINE) is None
    assert 'files_fastq=( "${dir_species_fastq}"/* )' not in text
    assert 'if [[ ! "${GG_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${GG_ARRAY_TASK_ID} -lt 1 ]]; then' in text
    assert 'if [[ ${#fastq_mode_dirs[@]} -eq 0 ]]; then' in text
    assert 'if [[ ${#sra_mode_files[@]} -eq 0 ]]; then' in text
    assert 'if [[ ${#metadata_mode_files[@]} -eq 0 ]]; then' in text

    fastq_guard = 'if [[ ${id} -ge ${#dirs[@]} ]]; then'
    fastq_use = 'dir_species_fastq="${dirs[${id}]}"'
    assert fastq_guard in text
    assert fastq_use in text
    assert text.index(fastq_guard) < text.index(fastq_use)

    sraid_guard = 'if [[ ${id} -ge ${#files[@]} ]]; then'
    sraid_use = 'file_input_sra_list="${files[${id}]}"'
    assert sraid_guard in text
    assert sraid_use in text
    assert text.index(sraid_guard) < text.index(sraid_use)

    metadata_guard = 'if [[ ${id} -ge ${#files[@]} ]]; then'
    metadata_use = 'file_metadata="${files[${id}]}"'
    assert metadata_use in text
    assert text.rindex(metadata_guard) < text.index(metadata_use)


def test_transcriptome_core_stores_generated_metadata_under_output_dir():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)

    assert 'dir_input_amalgkit_metadata="${gg_workspace_input_dir}/amalgkit_metadata"' in text
    assert 'dir_generated_amalgkit_metadata="${dir_transcriptome_assembly_output}/amalgkit_metadata"' in text
    assert 'file_input_amalgkit_metadata="${dir_input_amalgkit_metadata}/${sp_ub}_metadata.tsv"' in text
    assert 'file_generated_amalgkit_metadata="${dir_generated_amalgkit_metadata}/${sp_ub}_metadata.tsv"' in text
    assert 'if [[ "${selected_transcriptome_mode}" == "metadata" ]]; then' in text
    assert 'file_amalgkit_metadata="${file_input_amalgkit_metadata}"' in text
    assert 'file_amalgkit_metadata="${file_generated_amalgkit_metadata}"' in text
    assert 'file_amalgkit_metadata="${dir_amalgkit_metadata}/${sp_ub}_metadata.tsv"' not in text


def test_transcriptome_core_guards_getfastq_outputs_before_assembly_and_quant():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    dir_guard = 'if [[ ! -d "${dir_amalgkit_getfastq_sp}" ]]; then'
    fastq_guard = 'if [[ -z "$(find "${dir_amalgkit_getfastq_sp}" -type f -name "*.amalgkit.fastq.gz" -print -quit 2> /dev/null)" ]]; then'
    assert text.count(dir_guard) >= 2
    assert text.count(fastq_guard) >= 2
    assert "amalgkit getfastq output directory not found" in text
    assert "No amalgkit getfastq FASTQ files were found in" in text


def test_transcriptome_core_uses_type_f_for_fastq_find_patterns():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    expected_tokens = [
        'find "${dir_amalgkit_getfastq_sp}" -type f -name "*_2.amalgkit.fastq.gz"',
        'find "${dir_amalgkit_getfastq_sp}" -type f -name "*.amalgkit.fastq.gz"',
        'find "${dir_amalgkit_getfastq_sp}" -type f -name "*_1.amalgkit.fastq.gz"',
        'find "${selected_fastq_dir}" -type f -name "*.amalgkit.fastq.gz"',
        'find "${selected_fastq_dir}" -type f -name "*_1.amalgkit.fastq.gz"',
        'find "${selected_fastq_dir}" -type f -name "*_2.amalgkit.fastq.gz"',
        'find "${assembly_input_fastq_dir}" -type f -name "*.amalgkit.fastq.gz"',
        'find "${assembly_input_fastq_dir}" -type f -name "*_1.amalgkit.fastq.gz"',
        'find "${assembly_input_fastq_dir}" -type f -name "*_2.amalgkit.fastq.gz"',
    ]
    for token in expected_tokens:
        assert token in text

    banned_tokens = [
        'find "${dir_amalgkit_getfastq_sp}" -name "*_2.amalgkit.fastq.gz"',
        'find "${dir_amalgkit_getfastq_sp}" -name "*.amalgkit.fastq.gz"',
        'find "${dir_amalgkit_getfastq_sp}" -name "*_1.amalgkit.fastq.gz"',
        'find "${selected_fastq_dir}" -name "*.amalgkit.fastq.gz"',
        'find "${selected_fastq_dir}" -name "*_1.amalgkit.fastq.gz"',
        'find "${selected_fastq_dir}" -name "*_2.amalgkit.fastq.gz"',
        'find "${assembly_input_fastq_dir}" -name "*.amalgkit.fastq.gz"',
        'find "${assembly_input_fastq_dir}" -name "*_1.amalgkit.fastq.gz"',
        'find "${assembly_input_fastq_dir}" -name "*_2.amalgkit.fastq.gz"',
    ]
    for token in banned_tokens:
        assert token not in text


def test_gene_evolution_core_escapes_embedded_quotes_in_seqtype_error_message():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert 'echo "Unsupported sequence type \'${seqtype}\' in \'${file_query_gene}\'. Only "DNA" or "Protein" are allowed. Exiting."' not in text
    assert 'echo "Unsupported sequence type \'${seqtype}\' in \'${file_query_gene}\'. Only \\"DNA\\" or \\"Protein\\" are allowed. Exiting."' in text


def test_genome_annotation_core_guards_array_task_id_before_task_index_math():
    script = CORE_DIR / "gg_genome_annotation_core.sh"
    text = _read_text(script)
    assert 'if [[ ! -d "${dir_sp_cds}" ]]; then' in text
    assert 'echo "Input directory not found: ${dir_sp_cds}. Exiting."' in text
    assert 'find "${dir_sp_cds}" -maxdepth 1 -type f ! -name \'.*\'' in text
    assert 'find "${dir_sp_dnaseq}/${sp_ub}" -type f ! -name \'.*\'' in text
    guard = 'if [[ ! "${GG_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${GG_ARRAY_TASK_ID} -lt 1 ]]; then'
    task_index = "task_index=$((GG_ARRAY_TASK_ID - 1))"
    assert guard in text
    assert task_index in text
    assert text.index(guard) < text.index(task_index)


def test_genome_annotation_core_multispecies_summary_requires_real_summary_inputs():
    script = CORE_DIR / "gg_genome_annotation_core.sh"
    text = _read_text(script)
    assert "summary_inputs_available=0" in text
    assert 'No multispecies summary inputs are available yet. Skipping summary generation.' in text
    assert 'if is_output_older_than_inputs "^file_sp_" "${file_multispecies_summary}"; then' not in text
    assert 'if is_output_older_than_inputs "^(dir_summary_|file_summary_)" "${file_multispecies_summary}"; then' in text
    assert 'dir_summary_species_annotation="${gg_workspace_output_dir}/species_cds_annotation"' in text
    assert 'file_summary_species_tree_dated="${dir_summary_species_tree}/dated_species_tree.nwk"' in text


def test_genome_evolution_core_excludes_hidden_files_when_listing_species_proteins():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    assert 'find "${dir_sp_protein}" -maxdepth 1 -type f ! -name \'.*\'' in text


def test_input_generation_core_excludes_hidden_files_when_listing_species_inputs():
    script = CORE_DIR / "gg_input_generation_core.sh"
    text = _read_text(script)
    expected_tokens = [
        'find "${species_cds_dir}" -maxdepth 1 -type f ! -name \'.*\'',
        'find "${species_gff_dir}" -maxdepth 1 -type f ! -name \'.*\'',
    ]
    for token in expected_tokens:
        assert text.count(token) >= 2


def test_input_generation_core_runs_cds_gff_mapping_validation():
    script = CORE_DIR / "gg_input_generation_core.sh"
    text = _read_text(script)
    assert 'validate_cds_gff_mapping.py' in text
    assert '--species-cds-dir "${species_cds_dir}"' in text
    assert '--species-gff-dir "${species_gff_dir}"' in text
    assert '--nthreads "${GG_TASK_CPUS:-1}"' in text


def test_input_generation_core_populates_species_summary_taxonomy_metadata_nonfatally():
    script = CORE_DIR / "gg_input_generation_core.sh"
    text = _read_text(script)
    assert 'if ! ensure_ete_taxonomy_db "${gg_workspace_dir}"; then' in text
    assert "species_summary taxonomy metadata" in text
    assert "Continuing without taxid/genetic code annotation." in text


def test_mmseqs_uniref90_download_retries_and_reports_disk_context():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "_download_mmseqs_uniref90_db")
    assert 'local output_db="${db_dir}/${uniref_db}_DB"' in body
    assert 'mmseqs databases "${uniref_db}" "${output_db}" "${db_dir}" --threads "${nthreads}"' in body
    assert "for attempt in 1 2 3; do" in body
    assert 'Preparing MMseqs2 UniRef90 taxonomy DB in: ${db_dir} (attempt ${attempt}/${max_attempts})' in body
    assert 'MMseqs2 UniRef90 taxonomy DB preparation failed in: ${db_dir} (attempt ${attempt}/${max_attempts})' in body
    assert 'df -h "${db_dir}" >&2 || true' in body
    assert "sleep 5" in body


def test_genome_evolution_core_does_not_include_legacy_output_migration_code():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "legacy_busco_files=()",
        "legacy_singlecopy_files=()",
        "legacy_mafft_files=()",
        "legacy_trimal_files=()",
        "migrate_legacy_fasta_outputs()",
        "migrate_legacy_fasta_outputs \"${dir_busco_fasta}\"",
        "migrate_legacy_fasta_outputs \"${dir_busco_mafft}\"",
        "migrate_legacy_fasta_outputs \"${dir_busco_trimal}\"",
        "-name \"*_busco.full.tsv\"",
        "-name \"*_busco.short.txt\"",
    ]
    for token in banned_tokens:
        assert token not in text


def test_genome_evolution_core_does_not_generate_unused_partitions_txt():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "partitions.txt",
        "partition_start=",
        "partition_end=",
        "Failed to read alignment length from",
    ]
    for token in banned_tokens:
        assert token not in text
