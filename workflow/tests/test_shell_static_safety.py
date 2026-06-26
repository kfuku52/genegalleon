# ruff: noqa: E501

import re
from pathlib import Path

from shell_static_helpers import (
    CONTAINER_SCRIPTS_DIR,
    CORE_DIR,
    GITHUB_WORKFLOWS_DIR,
    WORKFLOW_DIR,
)
from shell_static_helpers import (
    container_shell_scripts as _container_shell_scripts,
)
from shell_static_helpers import (
    core_and_entrypoint_scripts as _core_and_entrypoint_scripts,
)
from shell_static_helpers import (
    function_body as _function_body,
)
from shell_static_helpers import (
    read_text as _read_text,
)
from shell_static_helpers import (
    strict_mode_header as _strict_mode_header,
)
from shell_static_helpers import (
    workflow_shell_scripts as _workflow_shell_scripts,
)


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
        WORKFLOW_DIR / "support" / "gg_shared_lock.sh",
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
    assert 'default_pycache_prefix="${TMPDIR:-/tmp}/genegalleon_pycache"' in body
    assert 'mkdir -p -- "${default_pycache_prefix}" 2>/dev/null || true' in body
    assert 'export PYTHONPYCACHEPREFIX="${default_pycache_prefix}"' in body
    runtime_body = _function_body(text, "gg_bootstrap_core_runtime")
    assert "gg_configure_python_pycacheprefix_from_core" in runtime_body


def test_conda_activation_temporarily_disables_nounset():
    text = _read_text(WORKFLOW_DIR / "support" / "gg_util.sh")
    activate_body = _function_body(text, "gg_activate_conda_env")
    deactivate_body = _function_body(text, "gg_deactivate_conda_env")
    assert "if [[ $- == *u* ]]; then" in activate_body
    assert "set +u" in activate_body
    assert 'conda activate "${conda_env}"' in activate_body
    assert "set -u" in activate_body
    assert "if [[ $- == *u* ]]; then" in deactivate_body
    assert "set +u" in deactivate_body
    assert "conda deactivate >/dev/null 2>&1 || true" in deactivate_body
    assert "set -u" in deactivate_body


def test_prepare_cmd_runtime_makes_ulimit_best_effort():
    text = _read_text(WORKFLOW_DIR / "support" / "gg_util.sh")
    body = _function_body(text, "gg_prepare_cmd_runtime")
    assert "ulimit -s unlimited 2>/dev/null || true" in body


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
        assert pattern.search(text) is None, f"Use arithmetic for-loops instead of `for ... in $(seq ...)`: {script}"


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
        assert pattern.search(text) is None, f"Avoid awk match(..., ..., array) for portability: {script}"


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
    assert (
        'gg_runtime_core_script="$(gg_prepare_entrypoint_runtime_snapshot "${gg_entrypoint_name}" "${gg_core_dir}/gg_progress_summary_core.sh")"'
        in text
    )
    assert 'gg_run_container_shell_script "${gg_container_image_path}" "${gg_runtime_core_script}"' in text
    assert "orthogroup_output_summary.py" not in text
    assert "transcriptome_assembly_output_summary.py" not in text


def test_no_known_unquoted_query_gene_file_expansions():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "head --bytes 1 ${file_query_gene}",
        "seqkit stats --tabular ${file_query_gene}",
        "seqkit translate --allow-unknown-codon --transl-table ${genetic_code} --threads ${GG_TASK_CPUS} ${file_query_gene}",
        'cp_out ${file_query_gene} ${dir_output_active}/query_gene/$(basename "${file_query_gene}")',
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
        'sed \'2,$ s/\\t/_.*\\t/\' "${file_sp_trait}" > "foreground.tsv"',
        "sed '2,$ s/\\t/_.*\\t/' species_trait_binary.tsv > foreground.tsv",
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
            assert brace_pat.search(line) is None, f"Unquoted path-like option expansion in {script}:{lineno}: {line}"
            assert plain_pat.search(line) is None, f"Unquoted path-like option expansion in {script}:{lineno}: {line}"


def test_no_unquoted_infile_option_expansions_in_core_scripts():
    pattern = re.compile(r"--infile\s+\$\{[^}]+\}")
    for script in sorted(CORE_DIR.glob("*.sh")):
        for lineno, line in enumerate(_read_text(script).splitlines(), start=1):
            stripped = line.lstrip()
            if stripped.startswith("#"):
                continue
            assert pattern.search(line) is None, f"Unquoted --infile expansion in {script}:{lineno}: {line}"


def test_no_unquoted_outfile_option_expansions_in_core_scripts():
    pattern = re.compile(r"--outfile\s+\$\{[^}]+\}")
    for script in sorted(CORE_DIR.glob("*.sh")):
        for lineno, line in enumerate(_read_text(script).splitlines(), start=1):
            stripped = line.lstrip()
            if stripped.startswith("#"):
                continue
            assert pattern.search(line) is None, f"Unquoted --outfile expansion in {script}:{lineno}: {line}"


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
    assert "python -c" in body
    assert " ${num_masked_bp} ${num_total_bp}" not in body


def test_is_output_older_than_inputs_uses_compgen_variable_listing():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "is_output_older_than_inputs")
    assert 'compgen -A variable | grep -E -- "${input_file_variable_regex}"' in body
    assert 'set | grep "${input_file_variable_regex}"' not in body


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


def test_convergent_sites_entrypoint_forwards_plot_runtime_envs():
    entrypoint = WORKFLOW_DIR / "gg_convergent_sites_entrypoint.sh"
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
    assert "CPU_PER_HOST=" not in text
    assert 'cpu_pick="${GG_TASK_CPUS}"' not in text
    assert '--nslots="${GG_TASK_CPUS}"' in text
    assert "taskset" not in text
    assert "cpu_id=$(python -c" not in text
    assert '"${l1ou_cmd[@]}"' in text


def test_gene_evolution_core_uses_kfl1ou_wrapper_with_supported_args_only():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    block_start = text.index('task="l1ou"')
    block_end = text.index("mv_out fit_ind.RData", block_start)
    l1ou_block = text[block_start:block_end]

    assert "detect_OU_shift_l1ou.r" not in l1ou_block
    assert 'Rscript "${gg_support_dir}/detect_OU_shift_kfl1ou.r"' in l1ou_block
    assert '--require_internal_node_labels="${require_internal_node_labels:-1}"' not in l1ou_block
    assert '--clade_collapse_similarity_method="${clade_collapse_similarity_method}"' not in l1ou_block
    assert '--clade_collapse_similarity_threshold="${clade_collapse_similarity_threshold}"' not in l1ou_block
    assert "--ceil_negative=0" not in l1ou_block
    assert '--replicate_sep="_"' in l1ou_block


def test_detect_ou_shift_kfl1ou_enables_measurement_error_by_default():
    script = WORKFLOW_DIR / "support" / "detect_OU_shift_kfl1ou.r"
    text = _read_text(script)
    assert "measurement_error = TRUE" in text
    assert "input_error = input_error_fit" in text


def test_gene_evolution_core_uses_explicit_ne_and_grouped_logic_for_tree_pruning_gate():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert "if [[ ! ${run_tree_pruning} -eq 1 && ${run_l1ou} -eq 1 ]]; then" not in text
    assert "if [[ ${run_tree_pruning} -ne 1 && ${run_l1ou} -eq 1 ]]; then" in text


def test_genome_evolution_core_uses_ne_for_busco_count_mismatch_checks():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "if [[ ! ${num_busco_ids} -eq ${num_singlecopy_fasta} && ${run_extract_species_tree_fasta} -eq 1 ]]; then",
        "if [[ ! ${num_busco_ids} -eq ${num_mafft_fasta} && ${run_individual_mafft} -eq 1 ]]; then",
        "if [[ ! ${num_busco_ids} -eq ${num_trimal_fasta} && ${run_individual_trimal} -eq 1 ]]; then",
        "if [[ ! ${num_busco_ids} -eq ${num_iqtree_pep} && ${run_individual_iqtree_pep} -eq 1 ]]; then",
        "if [[ ! ${num_busco_ids} -eq ${num_iqtree_dna} && ${run_individual_iqtree_dna} -eq 1 ]]; then",
    ]
    expected_tokens = [
        "if [[ ${num_busco_ids} -ne ${num_singlecopy_fasta} && ${run_extract_species_tree_fasta} -eq 1 ]]; then",
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
    assert "export SINGULARITYENV_PYTHONNOUSERSITE=1" in body
    assert "export APPTAINERENV_PYTHONNOUSERSITE=1" in body
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


def test_mv_out_replace_dir_replaces_existing_destination_tree_safely():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "mv_out_replace_dir")
    assert 'rm -rf -- "${dest_dir}"' in body
    assert 'ensure_parent_dir "${dest_dir}"' in body
    assert 'mv -- "${staged_dir}" "${dest_dir}"' in body
    assert 'mv "${staged_dir}" "${dest_dir}"' not in body


def test_resolve_rnaspades_transcript_fasta_checks_supported_output_names():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "resolve_rnaspades_transcript_fasta")
    assert '"transcripts.fasta"' in body
    assert '"soft_filtered_transcripts.fasta"' in body
    assert '"hard_filtered_transcripts.fasta"' in body


def test_capture_busco_repro_artifacts_recreates_dir_and_copies_artifacts_safely():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "capture_busco_repro_artifacts")
    assert 'recreate_dir "${repro_dir}"' in body
    assert 'gg_copy_file_portable "${input_fasta}" "${repro_dir}/"' in body
    assert 'gg_copy_dir_portable "${busco_tmp_dir}" "${repro_dir}/busco_tmp"' in body
    assert 'gg_copy_file_portable "${stderr_log}" "${repro_dir}/busco.stderr.log"' in body


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
    assert "gg_trigger_versions_dump: skipped existing ${log_file}" in body
    assert "timestamp=$(date" not in body


def test_gg_trigger_versions_dump_key_tracks_versions_script_and_repo_version():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "gg_trigger_versions_dump")
    assert (
        'container_key_seed="gg_container_image_path=${gg_container_image_path};runtime=${container_runtime_bin};gg_version=${gg_version}"'
        in body
    )
    assert "versions_script_hash=$(sha256sum \"${versions_script}\" | awk '{print $1}')" in body
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
        assert "Warning: gg_versions trigger failed." not in text, (
            f"Legacy warning-only versions dump handling remains: {script}"
        )
        assert 'gg_require_versions_dump "${gg_entrypoint_name}"' in text, (
            f"Entry point must fail on versions dump error: {script}"
        )


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


def test_entrypoints_apply_registered_env_overrides_before_forwarding_config():
    expected_tokens = {
        "gg_genome_annotation_entrypoint.sh": 'gg_apply_registered_env_overrides "${gg_entrypoint_name}" "delete_tmp_dir"',
        "gg_transcriptome_generation_entrypoint.sh": 'gg_apply_registered_env_overrides "${gg_entrypoint_name}" "delete_tmp_dir"',
        "gg_gene_evolution_entrypoint.sh": 'gg_apply_registered_env_overrides "${gg_entrypoint_name}" "delete_tmp_dir" "delete_preexisting_tmp_dir"',
    }
    entrypoints = sorted(WORKFLOW_DIR.glob("gg_*_entrypoint.sh"))
    assert entrypoints
    for script in entrypoints:
        text = _read_text(script)
        if script.name in expected_tokens:
            token = expected_tokens[script.name]
        else:
            token = 'gg_apply_registered_env_overrides "${gg_entrypoint_name}"'
        assert token in text, f"Missing registered env override call in {script.name}: {token}"


def test_convergent_sites_entrypoint_does_not_define_unused_delete_tmp_dir():
    text = _read_text(WORKFLOW_DIR / "gg_convergent_sites_entrypoint.sh")
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
    assert "gg_print_registered_config_summary \\" in text
    assert '"effective config summary (gg_genome_evolution_core.sh)" \\' in text
    assert "print_effective_genome_evolution_config_summary" in text.split("root_species_tree()", 1)[0]


def test_entrypoints_with_exit_if_running_call_duplicate_guard():
    scripts = [
        "gg_convergent_sites_entrypoint.sh",
        "gg_genome_annotation_entrypoint.sh",
        "gg_gene_evolution_entrypoint.sh",
    ]
    for script_name in scripts:
        text = _read_text(WORKFLOW_DIR / script_name)
        assert "exit_if_running=" in text
        assert "exit_if_running_qstat" in text or "gg_entrypoint_prepare_container_runtime 1" in text, (
            f"Missing duplicate-job guard call in {script_name}"
        )


def test_remove_empty_subdirs_calls_use_quoted_variable_arguments():
    pattern = re.compile(
        r"^[ \t]*remove_empty_subdirs[ \t]+(\$\{[A-Za-z_][A-Za-z0-9_]*\}|\$[A-Za-z_][A-Za-z0-9_]*)[ \t]*$",
        re.MULTILINE,
    )
    for script in sorted(CORE_DIR.glob("*.sh")):
        text = _read_text(script)
        assert pattern.search(text) is None, f"Unquoted variable passed to remove_empty_subdirs in {script}"


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


def test_genome_evolution_core_runs_mcmctree_time_scaling_in_scratch():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)

    expected_tokens = [
        'mcmctree_time_scale.py"',
        'iq2mc_work_dir=$(mktemp -d "${dir_tmp}/tmp.iq2mc.work.XXXXXX")',
        'scale_mcmctree_calibrations_file "${file_constrained_tree}" "${iq2mc_scaled_constraint_tree}" "${mcmctree_time_scale_factor}" "down"',
        '-te "${iq2mc_scaled_constraint_tree}"',
        'scale_mcmctree_calibrations_file "${iq2mc_work_dir}/iq2mc.rooted.nwk" "${file_iq2mc_rooted_tree}" "${mcmctree_time_scale_factor}" "up"',
        'mcmctree_work_dir=$(mktemp -d "${dir_tmp}/tmp.mcmctree.work.XXXXXX")',
        'extract_scaled_mcmctree_figtree "${mcmctree_work_dir}/$(basename "${file_mcmctree_raw_output}")" "${file_mcmctree_figtree_tre}" "${mcmctree_time_scale_factor}"',
        "Raw scaled MCMCTree output is not retained by default.",
    ]
    for token in expected_tokens:
        assert token in text, f"Missing MCMCTree scaling token: {token}"


def test_genome_evolution_core_uses_rerun_safe_mkdir_for_orthogroup_grampa_tmp_input():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    assert "mkdir ./tmp.orthogroup_grampa_indir" not in text
    assert "mkdir -p ./tmp.orthogroup_grampa_indir" in text
    assert '[[ "${file_name}" == "${og_id}"* ]]' not in text
    assert 'gg_orthogroup_file_matches_id "${file_name}" "${og_id}"' in text


def test_genome_evolution_core_quotes_grampa_output_and_cafe_option_values():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)

    assert (
        'busco_grampa "${dir_busco_rooted_nwk_dna}" "$(dirname "${file_busco_grampa_dna}")" ${file_busco_grampa_dna}'
        not in text
    )
    assert (
        'busco_grampa "${dir_busco_rooted_nwk_pep}" "$(dirname "${file_busco_grampa_pep}")" ${file_busco_grampa_pep}'
        not in text
    )
    assert (
        'busco_grampa "./tmp.orthogroup_grampa_indir" "$(dirname "${file_orthogroup_grampa}")" ${file_orthogroup_grampa}'
        not in text
    )
    assert "--genecount ${file_orthogroup_genecount_selected}" not in text
    assert "--dated_species_tree ${file_dated_species_tree}" not in text
    assert "--max_size_differential ${max_size_differential_cafe}" not in text
    assert "--tree ${file_dated_species_tree}" not in text
    assert "--n_gamma_cats ${n_gamma_cats_cafe}" not in text
    assert "--cores ${GG_TASK_CPUS}" not in text
    assert "--output_prefix ${dir_cafe_output}" not in text
    assert "--file_cafe_input=${dir_cafe_orthogroup_selection}/cafe_input.tsv" not in text
    assert "--file_sptree=${file_dated_species_tree}" not in text
    assert "--file_trait=${file_trait}" not in text

    assert (
        'busco_grampa "${dir_busco_rooted_nwk_dna}" "$(dirname "${file_busco_grampa_dna}")" "${file_busco_grampa_dna}"'
        in text
    )
    assert (
        'busco_grampa "${dir_busco_rooted_nwk_pep}" "$(dirname "${file_busco_grampa_pep}")" "${file_busco_grampa_pep}"'
        in text
    )
    assert (
        'busco_grampa "./tmp.orthogroup_grampa_indir" "$(dirname "${file_orthogroup_grampa}")" "${file_orthogroup_grampa}"'
        in text
    )
    assert '--genecount "${file_orthogroup_genecount_selected}"' in text
    assert '--dated_species_tree "${file_dated_species_tree}"' in text
    assert '--max_size_differential "${max_size_differential_cafe}"' in text
    assert '--tree "${file_dated_species_tree}"' in text
    assert '--n_gamma_cats "${n_gamma_cats_cafe}"' in text
    assert '--cores "${GG_TASK_CPUS}"' in text
    assert '--output_prefix "${dir_cafe_output}"' in text
    assert '--file_cafe_input="${dir_cafe_orthogroup_selection}/cafe_input.tsv"' in text
    assert '--file_sptree="${file_dated_species_tree}"' in text
    assert '--file_trait="${file_trait}"' in text


def test_genome_evolution_core_builds_grampa_arguments_with_array():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    assert 'h1_param="-h1 ' not in text
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
    assert '-te "${iq2mc_scaled_constraint_tree}"' in iq2mc_block
    assert '--mcmc-bds "${mcmc_birth_death_sampling}"' in iq2mc_block
    assert '--mcmc-clock "${mcmc_clock_model}"' in iq2mc_block
    assert '--mcmc-iter "${mcmc_burnin},${mcmc_sampfreq},${mcmc_nsample}"' in iq2mc_block
    assert '-T "${GG_TASK_CPUS}"' in iq2mc_block

    busco_summary_start = text.index('python "${gg_support_dir}/collect_common_BUSCO_genes.py" \\')
    busco_summary_end = text.index('--outfile "tmp.busco_summary_table.tsv"', busco_summary_start) + len(
        '--outfile "tmp.busco_summary_table.tsv"'
    )
    busco_summary_block = text[busco_summary_start:busco_summary_end]
    assert "--busco_outdir ${dir_species_busco_full}" not in busco_summary_block
    assert "--ncpu ${GG_TASK_CPUS}" not in busco_summary_block
    assert '--busco_outdir "${dir_species_busco_full}"' in busco_summary_block
    assert '--ncpu "${GG_TASK_CPUS}"' in busco_summary_block


def test_genome_evolution_core_does_not_fallback_to_nonroot_hog_table():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    start = text.index('hog_table="${dir_orthofinder}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"')
    end = text.index('task="OMArk analysis of species-wise protein input files"', start)
    hog_block = text[start:end]

    assert "Falling back" not in hog_block
    assert 'hog_table="${hog_candidates[0]}"' not in hog_block
    assert "Required root-level HOG table was not found" in hog_block
    assert "will not be selected automatically" in hog_block
    assert 'orthogroup_table="OG"' in hog_block


def test_genome_evolution_core_treats_v31_orthogroups_as_root_hog_equivalent():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)

    assert "detect_orthofinder_version() {" in text
    assert "orthofinder_supports_root_hog_equivalent() {" in text
    assert "copy_root_hog_equivalent_from_orthogroups() {" in text
    assert 'if orthofinder_supports_root_hog_equivalent "${orthofinder_version}"; then' in text
    assert "OrthoFinder version >= 3.1" in text
    assert "Orthogroups/Orthogroups.tsv is treated as the root-level HOG equivalent" in text
    assert 'cp_out "${source_og}" "${target_dir}/Orthogroups.tsv"' in text
    assert 'cp_out "${source_genecount}" "${target_dir}/Orthogroups.GeneCount.tsv"' in text


def test_common_params_do_not_define_contamination_removal_rank():
    text = _read_text(WORKFLOW_DIR / "gg_common_params.sh")
    assert "GG_COMMON_CONTAMINATION_REMOVAL_RANK" not in text


def test_common_busco_lineage_defaults_to_eukaryota_odb12():
    text = _read_text(WORKFLOW_DIR / "gg_common_params.sh")
    assert ': "${GG_COMMON_BUSCO_LINEAGE:=eukaryota_odb12}"' in text


def test_common_input_sequence_mode_defaults_to_cds():
    text = _read_text(WORKFLOW_DIR / "gg_common_params.sh")
    assert ': "${GG_COMMON_INPUT_SEQUENCE_MODE:=cds}"' in text


def test_common_csubst_nonsyn_recode_defaults_to_no():
    text = _read_text(WORKFLOW_DIR / "gg_common_params.sh")
    assert ': "${GG_COMMON_CSUBST_NONSYN_RECODE:=no}"' in text


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

    assert (
        'grampa_h1="" # Optional GRAMPA H1 hypothesis. Leave empty to skip GRAMPA steps. Example: "2" or "x,y,z".'
        in entrypoint
    )
    assert (
        'target_branch_go="" # Optional GO-enrichment target branch. Leave empty to skip GO enrichment. Example: "<1>" or "Arabidopsis_thaliana".'
        in entrypoint
    )
    assert "GG_COMMON_GRAMPA_H1" not in core
    assert "GG_COMMON_TARGET_BRANCH_GO" not in core
    assert 'grampa_h1="${grampa_h1:-}"' in core
    assert 'target_branch_go="${target_branch_go:-}"' in core
    assert (
        "Disabling GRAMPA tasks because grampa_h1 is empty. Set grampa_h1 in gg_genome_evolution_entrypoint.sh to enable them."
        in core
    )
    assert (
        "Disabling run_go_enrichment because target_branch_go is empty. Set target_branch_go in gg_genome_evolution_entrypoint.sh to enable it."
        in core
    )
    assert ': "${grampa_h1:?' not in core
    assert ': "${target_branch_go:?' not in core
    assert "grampa_h1" in config_vars
    assert "target_branch_go" in config_vars


def test_cafe_go_enrichment_resolves_target_branch_by_exact_column_name():
    text = _read_text(WORKFLOW_DIR / "support" / "cafe_go_enrichment.r")

    assert "target_branch_id <- grep(target_branch" not in text
    assert "resolve_target_branch_id <- function(target_branch, branch_columns)" in text
    assert "exact_matches <- branch_columns[branch_columns == target_branch]" in text
    assert "grep(target_branch, branch_columns, value = TRUE, fixed = TRUE)" in text


def test_genome_evolution_exposes_cafe_trait_pgls_parameters():
    entrypoint = _read_text(WORKFLOW_DIR / "gg_genome_evolution_entrypoint.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")
    core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")

    expected_entrypoint_tokens = [
        "run_cafe_trait_pgls=0 # Test associations between CAFE-selected orthogroup copy numbers and species traits with species-tree PGLS.",
        'cafe_trait="all" # Trait column name(s) in species_trait.tsv to test against CAFE-selected copy numbers, or "all".',
        "cafe_trait_min_species=4 # Minimum number of tree-matched species required for each CAFE copy-number trait PGLS fit.",
        'cafe_trait_family_ids="" # Optional comma/space-separated CAFE family IDs to test; empty means use max_families.',
        'cafe_trait_family_file="" # Optional file listing CAFE family IDs to test.',
        'cafe_trait_max_families="all" # Maximum CAFE families tested: all|auto|0 for unlimited, or a non-negative integer.',
        'file_trait="auto" # Species trait table path for CAFE copy-number trait PGLS, or auto for workspace/input/species_trait/species_trait.tsv.',
    ]
    for token in expected_entrypoint_tokens:
        assert token in entrypoint

    for token in [
        "run_cafe_trait_pgls",
        "cafe_trait",
        "cafe_trait_min_species",
        "cafe_trait_family_ids",
        "cafe_trait_family_file",
        "cafe_trait_max_families",
        "cafe_trait_p_adjust_method",
        "cafe_trait_alpha",
        "cafe_trait_plot_top_n",
        "file_trait",
    ]:
        assert token in config_vars
        assert token in core

    assert 'file_trait="${gg_workspace_input_dir}/species_trait/species_trait.tsv"' in core
    assert 'task="CAFE copy-number trait PGLS"' in core
    assert 'Rscript "${gg_support_dir}/cafe_trait_pgls.r"' in core


def test_genome_evolution_uses_local_species_tree_rooting_parameter():
    entrypoint = _read_text(WORKFLOW_DIR / "gg_genome_evolution_entrypoint.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")
    core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")
    common = _read_text(WORKFLOW_DIR / "gg_common_params.sh")

    assert (
        'species_tree_rooting="taxonomy" # taxonomy[,ncbi[,opentree,timetree...]] | outgroup,GENUS_SPECIES[,GENUS_SPECIES...] | midpoint | mad | mv; selects how species trees are rooted before dating, using taxonomy providers, explicit outgroups, or topology/branch-length rooting methods.'
        in entrypoint
    )
    assert "GG_COMMON_OUTGROUP_LABELS" not in common
    assert 'species_tree_rooting="${species_tree_rooting:-taxonomy}"' in core
    assert (
        'parse_species_tree_rooting "${species_tree_rooting}" species_tree_rooting_method species_tree_rooting_value'
        in core
    )
    assert (
        'species_tree_rooting must be one of "outgroup,GENUS_SPECIES[,GENUS_SPECIES...]", "midpoint", "mad", "mv", or "taxonomy[,ncbi[,opentree,timetree...]]".'
        in core
    )
    assert 'nwkit_root_args=(--method "${root_method}" --infile "${infile}" --outfile "${outfile}")' in core
    assert 'nwkit_root_args+=(--outgroup "${root_value}")' in core
    assert 'nwkit_root_args+=(--download_dir "${dir_nwkit_download_dir}")' in core
    assert "species_tree_rooting" in config_vars


def test_genome_evolution_supports_protein_input_mode_and_species_code_overrides():
    entrypoint = _read_text(WORKFLOW_DIR / "gg_genome_evolution_entrypoint.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")
    core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")
    util = _read_text(WORKFLOW_DIR / "support" / "gg_util.sh")

    assert (
        'input_sequence_mode="${input_sequence_mode:-${GG_COMMON_INPUT_SEQUENCE_MODE:-cds}}" # {cds,protein}; protein mode uses species_protein inputs or per-species CDS->protein translation with optional species_genetic_code/species_genetic_code.tsv overrides.'
        in entrypoint
    )
    assert "input_sequence_mode" in config_vars
    assert 'input_sequence_mode="${input_sequence_mode:-${GG_COMMON_INPUT_SEQUENCE_MODE:-cds}}"' in core
    assert 'input_sequence_mode=$(gg_normalize_input_sequence_mode "${input_sequence_mode}")' in core
    assert "gg_normalize_input_sequence_mode() {" in util
    assert "species_genetic_code_table_path() {" in core
    assert 'echo "${gg_workspace_input_dir}/species_genetic_code/species_genetic_code.tsv"' in core
    assert 'echo "${gg_workspace_input_dir}/species_protein"' in core
    assert "prepare_species_genetic_code_table() {" in core
    assert "lookup_species_genetic_code() {" in core
    assert 'check_species_protein_dir "${dir_sp_protein_input}"' in core
    assert "species_genetic_code.tsv is ignored because species_protein inputs are provided" in core
    assert 'if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then' in core
    assert "Ignoring species_protein inputs in cds mode: ${dir_sp_protein_input}" in core
    assert "cds mode always generates temporary species_protein FASTA files from species_cds." in core
    assert "run_cds_translation must be 1 when species proteins need to be generated from species_cds." in core
    assert "Translation started: ${cds} (genetic_code=${species_code})" in core
    assert 'translated_file="${sp_ub}.fa.gz"' in core
    assert 'Copying protein FASTA: $(basename "${protein_path}") -> ${translated_file}' in core
    assert "stage_species_protein_fasta() {" in core
    assert "prepare_species_protein_orthofinder_dir() {" in core
    assert 'prepare_species_protein_orthofinder_dir "${dir_sp_protein}" "${dir_sp_protein_orthofinder}"' in core
    assert (
        'refresh_dir_for_shared_protein_input_signature "${dir_genome_evolution}" "genome_evolution" "${shared_protein_input_signature}"'
        in core
    )
    assert (
        'mapfile -t annotation_species_candidates < <(gg_species_names_from_fasta_dir "${dir_sp_protein_input}")'
        in core
    )
    protein_candidates_index = core.index(
        'if [[ ${#annotation_species_candidates[@]} -eq 0 && "${input_sequence_mode}" == "protein" ]]; then'
    )
    cds_candidates_index = core.index(
        'if [[ ${#annotation_species_candidates[@]} -eq 0 ]]; then\n  mapfile -t annotation_species_candidates < <(gg_species_names_from_fasta_dir "${dir_sp_cds}")\nfi'
    )
    assert protein_candidates_index < cds_candidates_index


def test_gene_evolution_uses_shared_input_mode_and_limits_protein_mode_to_supported_stages():
    entrypoint = _read_text(WORKFLOW_DIR / "gg_gene_evolution_entrypoint.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")
    core = _read_text(CORE_DIR / "gg_gene_evolution_core.sh")
    util = _read_text(WORKFLOW_DIR / "support" / "gg_util.sh")

    gene_block_start = config_vars.index("gg_gene_evolution_entrypoint.sh)")
    gene_block_end = config_vars.index("gg_hgt_entrypoint.sh)")
    gene_block = config_vars[gene_block_start:gene_block_end]

    assert (
        'input_sequence_mode="${input_sequence_mode:-${GG_COMMON_INPUT_SEQUENCE_MODE:-cds}}" # {cds,protein}; protein mode is partial and deactivates CDS-only analyses.'
        in entrypoint
    )
    assert "\ninput_sequence_mode\n" in gene_block
    assert 'input_sequence_mode="${input_sequence_mode:-${GG_COMMON_INPUT_SEQUENCE_MODE:-cds}}"' in core
    assert 'input_sequence_mode=$(gg_normalize_input_sequence_mode "${input_sequence_mode}")' in core
    assert "disable_flag_with_reason() {" in util
    assert "gg_species_genetic_code_table_path() {" in util
    assert "gg_species_protein_input_dir_path() {" in util
    assert "gg_prepare_species_genetic_code_table() {" in util
    assert "gg_lookup_species_genetic_code() {" in util
    assert "apply_gene_evolution_input_sequence_mode() {" in core
    assert "gene_evolution_model_is_aa() {" in core
    assert "assert_gene_evolution_aa_model_for_protein_mode() {" in core
    assert 'if [[ "${input_sequence_mode}" != "protein" ]]; then' in core
    assert 'if [[ "${mode_gene_evolution}" == "query2family" ]]; then' in core
    assert "query2family currently requires species_cds-backed search and CDS extraction." in core
    assert 'dir_sp_protein_input="$(gg_species_protein_input_dir_path "${gg_workspace_input_dir}")"' in core
    assert 'file_species_genetic_code="$(gg_species_genetic_code_table_path "${gg_workspace_input_dir}")"' in core
    assert 'file_og_pep_fasta="${dir_output_active}/protein_fasta/${og_id}_pep.fa.gz"' in core
    assert (
        'if [[ "${input_sequence_mode}" == "protein" ]]; then\n  file_og_primary_fasta="${file_og_pep_fasta}"' in core
    )
    assert 'if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then' in core
    assert 'check_species_protein_dir "${dir_sp_protein_input}"' in core
    assert (
        'gg_prepare_species_genetic_code_table "${dir_sp_cds}" "${genetic_code}" "${file_species_genetic_code_resolved}" "${file_species_genetic_code}"'
        in core
    )
    assert (
        'translate_orthogroup_cds_to_protein_fasta "${file_og_cds_fasta}" "${file_og_pep_fasta}" "${file_species_genetic_code_resolved}"'
        in core
    )
    assert 'assert_gene_evolution_aa_model_for_protein_mode "${task}"' in core
    assert 'disable_if_no_input_file "run_collect_gff_info" "${file_og_primary_fasta}"' in core
    assert (
        'seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_primary_fasta}" --out-file "${og_id}.gff2genestat_input.fasta"'
        in core
    )
    assert 'disable_if_no_input_file "run_extract_promoter_fasta" "${file_og_gff_info}"' in core
    assert 'disable_if_no_input_file "run_fimo" "${file_og_promoter_fasta}" "${jaspar_path}"' in core
    assert (
        'disable_flag_with_reason "run_mapdnds_parameter_estimation" "input_sequence_mode=protein: mapdNdS parameter estimation requires codon alignments."'
        in core
    )
    assert 'disable_flag_with_reason "run_collect_gff_info"' not in core
    assert 'disable_flag_with_reason "run_scm_intron"' not in core
    assert 'disable_flag_with_reason "run_extract_promoter_fasta"' not in core
    assert 'disable_flag_with_reason "run_fimo"' not in core
    assert 'disable_flag_with_reason "treevis_synteny"' not in core
    assert 'synteny_source_dir="${dir_sp_cds}"' in core
    assert (
        'if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then\n    synteny_source_dir="${dir_sp_protein_input}"'
        in core
    )
    assert '--input_sequence_mode "${synteny_sequence_mode}"' in core
    assert (
        'apply_gene_evolution_input_sequence_mode\nif [[ "${mode_gene_evolution}" == "query2family" && ${run_query_blast} -eq 1 ]]; then'
        in core
    )


def test_gene_evolution_supports_auto_query_blast_evalue_by_query_length():
    entrypoint = _read_text(WORKFLOW_DIR / "gg_gene_evolution_entrypoint.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")
    core = _read_text(CORE_DIR / "gg_gene_evolution_core.sh")

    gene_block_start = config_vars.index("gg_gene_evolution_entrypoint.sh)")
    gene_block_end = config_vars.index("gg_hgt_entrypoint.sh)")
    gene_block = config_vars[gene_block_start:gene_block_end]

    assert 'query_blast_evalue="auto"' in entrypoint
    assert 'query_blast_auto_evalue_maxlen_cutoffs="40:1000,80:100,150:10,300:1,inf:0.01"' in entrypoint
    assert "\nquery_blast_auto_evalue_maxlen_cutoffs\n" in gene_block
    assert 'query_blast_evalue="${query_blast_evalue:-auto}"' in core
    assert "resolve_query_blast_evalue() {" in core
    assert "max_len <= cutoff + 0" in core
    assert (
        'resolve_query_blast_evalue "${query_blast_evalue}" "${query_aa_local}" "${query_blast_auto_evalue_maxlen_cutoffs}"'
        in core
    )
    assert '-evalue "${effective_query_blast_evalue}"' in core
    assert '--evalue "${effective_query_blast_evalue}"' in core


def test_genome_evolution_places_run_cds_translation_before_dependent_run_flags():
    entrypoint = _read_text(WORKFLOW_DIR / "gg_genome_evolution_entrypoint.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")

    entrypoint_translation_index = entrypoint.index("run_cds_translation=1")
    entrypoint_species_busco_index = entrypoint.index("run_species_busco=1")
    entrypoint_species_omark_index = entrypoint.index("run_species_omark=0")
    entrypoint_orthofinder_index = entrypoint.index("run_orthofinder=1")
    entrypoint_busco_getfasta_index = entrypoint.index("run_busco_dupaware_extract_fasta=")

    assert entrypoint_translation_index < entrypoint_species_busco_index
    assert entrypoint_translation_index < entrypoint_species_omark_index
    assert entrypoint_translation_index < entrypoint_orthofinder_index
    assert entrypoint_translation_index < entrypoint_busco_getfasta_index

    config_translation_index = config_vars.index("run_cds_translation")
    config_species_busco_index = config_vars.index("run_species_busco")
    config_species_omark_index = config_vars.index("run_species_omark")
    config_species_omark_summary_index = config_vars.index("run_build_species_omark_summary")
    config_orthofinder_index = config_vars.index("run_orthofinder")
    config_og_selection_index = config_vars.index("run_og_selection")
    config_busco_getfasta_index = config_vars.index("run_busco_dupaware_extract_fasta")

    assert config_translation_index < config_species_busco_index
    assert config_translation_index < config_species_omark_index
    assert config_translation_index < config_orthofinder_index
    assert config_translation_index < config_busco_getfasta_index
    assert (
        config_orthofinder_index
        < config_species_omark_index
        < config_species_omark_summary_index
        < config_og_selection_index
    )


def test_genome_evolution_exposes_omark_controls_and_summary_stage():
    entrypoint = _read_text(WORKFLOW_DIR / "gg_genome_evolution_entrypoint.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")
    core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")

    assert (
        "run_species_omark=0 # Run OMArk proteome quality assessment for each species using shared protein inputs."
        in entrypoint
    )
    assert (
        "run_build_species_omark_summary=1 # Build the shared OMArk summary table for species-wise proteome quality assessment."
        in entrypoint
    )
    assert (
        'omark_db_path="auto" # Path to the OMArk OMAmer LUCA.h5 database, or "auto" to download it under workspace/downloads/omark/.'
        in entrypoint
    )
    assert "run_species_omark" in config_vars
    assert "run_build_species_omark_summary" in config_vars
    assert "omark_db_path" in config_vars
    assert 'omark_db_path="${omark_db_path:-auto}"' in core
    assert 'dir_species_omamer="${dir_genome_evolution}/omamer_search"' in core
    assert 'dir_species_omark="${dir_genome_evolution}/omark"' in core
    assert 'file_species_omark_summary_table="${dir_genome_evolution}/omark_summary_table/omark_summary.tsv"' in core
    assert 'ensure_omark_database "${gg_workspace_dir}" "${omark_db_path}"' in core
    assert "omamer search \\" in core
    assert "omark \\" in core
    assert 'python "${gg_support_dir}/summarize_omark.py" \\' in core


def test_genome_evolution_entrypoint_groups_omark_with_orthofinder_stage():
    entrypoint = _read_text(WORKFLOW_DIR / "gg_genome_evolution_entrypoint.sh")

    species_tree_heading_index = entrypoint.index("# Species-tree workflow flags")
    species_tree_fasta_index = entrypoint.index("run_extract_species_tree_fasta=1")
    orthogroup_heading_index = entrypoint.index("# Orthogroup and species-protein QC workflow flags")
    orthofinder_index = entrypoint.index("run_orthofinder=1")
    omark_index = entrypoint.index("run_species_omark=0")
    omark_summary_index = entrypoint.index("run_build_species_omark_summary=1")
    og_selection_index = entrypoint.index("run_og_selection=1")

    assert species_tree_heading_index < species_tree_fasta_index < orthogroup_heading_index
    assert orthogroup_heading_index < orthofinder_index
    assert orthofinder_index < omark_index < omark_summary_index < og_selection_index


def test_genome_evolution_core_defaults_shared_protein_flags_for_legacy_launchers():
    core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")

    assert 'run_cds_translation="${run_cds_translation:-1}"' in core
    assert 'run_species_omark="${run_species_omark:-0}"' in core
    assert 'run_build_species_busco_summary="${run_build_species_busco_summary:-1}"' in core
    assert 'run_build_species_omark_summary="${run_build_species_omark_summary:-1}"' in core
    assert 'run_extract_species_tree_fasta="${run_extract_species_tree_fasta:-1}"' in core


def test_common_params_and_rooting_helpers_expose_shared_species_label_parser():
    common = _read_text(WORKFLOW_DIR / "gg_common_params.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")
    genome_core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")
    gene_core = _read_text(CORE_DIR / "gg_gene_evolution_core.sh")
    rooting_helper = _read_text(WORKFLOW_DIR / "support" / "species_tree_guided_gene_tree_rooting.r")
    treevis_core = _read_text(WORKFLOW_DIR / "support" / "treevis" / "R" / "00_core.R")

    assert ': "${GG_COMMON_SPECIES_LABEL_PARSER:=taxonomic}"' in common
    assert ': "${GG_COMMON_SPECIES_LABEL_REGEX:=}"' in common
    assert ': "${GG_COMMON_SPECIES_LABEL_MAP_TSV:=}"' in common
    assert "species_label_parser" in config_vars
    assert config_vars.count("species_label_regex") >= 2
    assert config_vars.count("species_label_map_tsv") >= 2
    assert 'species_label_parser="${species_label_parser:-${GG_COMMON_SPECIES_LABEL_PARSER:-taxonomic}}"' in genome_core
    assert 'species_label_regex="${species_label_regex:-${GG_COMMON_SPECIES_LABEL_REGEX:-}}"' in genome_core
    assert 'species_label_map_tsv="${species_label_map_tsv:-${GG_COMMON_SPECIES_LABEL_MAP_TSV:-}}"' in genome_core
    assert '"--species_parser=${species_label_parser}" \\' in genome_core
    assert "species_parser = args[['species_parser']]" in rooting_helper
    assert "species_parser = 'taxonomic'" in rooting_helper
    assert "get_species_overlap_score(phy=rt, dc_cutoff=0, species_parser=species_parser)" in rooting_helper
    assert 'species_label_parser="${species_label_parser:-${GG_COMMON_SPECIES_LABEL_PARSER:-taxonomic}}"' in gene_core
    assert 'species_label_regex="${species_label_regex:-${GG_COMMON_SPECIES_LABEL_REGEX:-}}"' in gene_core
    assert 'species_label_map_tsv="${species_label_map_tsv:-${GG_COMMON_SPECIES_LABEL_MAP_TSV:-}}"' in gene_core
    assert 'TREEVIS_SPECIES_PARSER="${species_label_parser}" \\' in gene_core
    assert "treevis_species_parser = Sys.getenv('TREEVIS_SPECIES_PARSER', unset='taxonomic')" in treevis_core
    assert "get_species_name(spp, species_parser=treevis_species_parser)" in treevis_core


def test_gene_evolution_core_passes_species_label_parser_options_to_radte():
    core = _read_text(CORE_DIR / "gg_gene_evolution_core.sh")

    assert 'radte_args+=("--species-parser=${species_label_parser}")' in core
    assert 'if [[ -n "${species_label_regex}" ]]; then' in core
    assert 'radte_args+=("--species-regex=${species_label_regex}")' in core
    assert 'if [[ -n "${species_label_map_tsv}" ]]; then' in core
    assert 'radte_args+=("--species-map-tsv=${species_label_map_tsv}")' in core


def test_nwkit_call_sites_receive_species_label_parser_options():
    gene_core = _read_text(CORE_DIR / "gg_gene_evolution_core.sh")
    genome_core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")

    assert 'nwkit_root_args+=(--species-parser "${species_label_parser}")' in gene_core
    assert 'nwkit_root_args+=(--species-regex "${species_label_regex}")' in gene_core
    assert 'nwkit_root_args+=(--species-map-tsv "${species_label_map_tsv}")' in gene_core
    assert 'nwkit_root_args+=(--species-parser "${species_label_parser}")' in genome_core
    assert 'nwkit_root_args+=(--species-regex "${species_label_regex}")' in genome_core
    assert 'nwkit_root_args+=(--species-map-tsv "${species_label_map_tsv}")' in genome_core
    assert '--species-parser "${species_label_parser}"' in genome_core
    assert 'nwkit_args+=(--species-regex "${species_label_regex}")' in genome_core
    assert 'nwkit_args+=(--species-map-tsv "${species_label_map_tsv}")' in genome_core


def test_genome_evolution_reuse_check_precedes_busco_lineage_resolution():
    core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")

    missing_outputs_index = core.index("if [[ ${missing_busco_outputs} -ne 1 ]]; then")
    reuse_return_index = core.index("    return 0\n  fi", missing_outputs_index)
    resolve_lineage_index = core.index('if ! resolve_busco_lineage_for_species_set "${input_species_set[@]}"; then')

    assert reuse_return_index < resolve_lineage_index


def test_genome_evolution_core_uses_canonical_duplicate_aware_busco_flag_names():
    core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")

    assert 'run_busco_dupaware_extract_fasta="${run_busco_dupaware_extract_fasta:-0}"' in core
    assert 'run_busco_dupaware_grampa_pep="${run_busco_dupaware_grampa_pep:-0}"' in core
    assert "run_busco_dupaware_extract_fasta" in config_vars


def test_genome_evolution_exposes_single_copy_ortholog_decay_plot():
    entrypoint = _read_text(WORKFLOW_DIR / "gg_genome_evolution_entrypoint.sh")
    core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")

    assert "run_single_copy_ortholog_decay_plot=1" in entrypoint
    assert "orthogroup_decay_replicates=1000" in entrypoint
    assert 'orthogroup_decay_species_counts="auto"' in entrypoint
    assert "orthogroup_decay_seed=1" in entrypoint
    assert "run_single_copy_ortholog_decay_plot" in config_vars
    assert "orthogroup_decay_replicates" in config_vars
    assert "orthogroup_decay_species_counts" in config_vars
    assert "orthogroup_decay_seed" in config_vars
    assert 'run_single_copy_ortholog_decay_plot="${run_single_copy_ortholog_decay_plot:-1}"' in core
    assert 'task="Single-copy ortholog decay plot"' in core
    assert 'python "${gg_support_dir}/single_copy_ortholog_decay_plot.py" \\' in core
    assert '--orthogroup-genecount "${orthogroup_decay_genecount}"' in core
    assert '--replicates "${orthogroup_decay_replicates}"' in core
    assert '--species-counts "${orthogroup_decay_species_counts}"' in core
    assert '--seed "${orthogroup_decay_seed}"' in core


def test_genome_evolution_runs_omark_after_orthofinder_and_before_og_selection():
    core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")

    orthofinder_index = core.index(
        'task="OrthoFinder"\nif [[ ! -s "${file_orthofinder_done_marker}" && ${run_orthofinder} -eq 1 ]]; then'
    )
    omark_index = core.index(
        'task="OMArk analysis of species-wise protein input files"\nrun_shared_species_omark_stage'
    )
    omark_summary_index = core.index('task="Summarizing OMArk species quality results"\nrun_shared_omark_summary_stage')
    og_selection_index = core.index('task="Selecting orthogroups based on gene and species numbers"')

    assert orthofinder_index < omark_index
    assert omark_index < omark_summary_index
    assert omark_summary_index < og_selection_index


def test_genome_evolution_protein_mode_allows_species_tree_plotting_with_pep_trees_only():
    core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")

    assert (
        'if [[ "${input_sequence_mode}" == "protein" ]]; then\n    disable_if_no_input_file "run_plot_species_trees" "${file_concat_iqtree_pep_root}" "${file_astral_tree_pep}"'
        in core
    )
    assert (
        'disable_if_no_input_file "run_plot_species_trees" "${file_concat_iqtree_dna_root}" "${file_concat_iqtree_pep_root}" "${file_astral_tree_dna}" "${file_astral_tree_pep}"'
        in core
    )


def test_plot_species_trees_r_filters_missing_tree_inputs():
    script = _read_text(WORKFLOW_DIR / "support" / "plot_species_trees.r")

    assert "has_nonempty_file <- function(file_path)" in script
    assert "tree_specs = Filter(function(spec) has_nonempty_file(spec[['nwk']]), tree_specs)" in script
    assert "if (length(tree_specs) == 0) {" in script
    assert "plot_ncol = min(2, length(plots))" in script
    assert "out = plot_list(gglist = plots, ncol = plot_ncol)" in script


def test_genome_evolution_protein_mode_disables_incompatible_dna_and_busco_steps():
    core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")

    assert 'if [[ "${input_sequence_mode}" == "protein" ]]; then' in core
    assert "protein mode does not support undated_species_tree=${undated_species_tree}." in core
    assert "Use iqtree_pep or astral_pep instead." in core
    assert "Disabling DNA-only species-tree steps in protein mode" in core
    assert "Disabling CDS-only dating steps in protein mode" in core
    assert "Disabling DNA-only duplicate-aware BUSCO genome-evolution steps in protein mode" in core
    assert 'dir_species_busco_full="${gg_workspace_output_dir}/species_protein_busco_full"' in core
    assert 'dir_species_busco_short="${gg_workspace_output_dir}/species_protein_busco_short"' in core
    assert "prepare_species_tree_input_dir" in core
    assert '--mode "${species_tree_busco_mode}"' in core
    assert 'outfile2="${dir_busco_fasta}/${busco_id}${genome_busco_fasta_suffix}"' in core
    assert "outfile=${dir_busco_mafft}/${infile_base}${genome_busco_aln_suffix}" in core
    assert 'outfile="${dir_busco_trimal}/${infile_base}${genome_busco_trimal_suffix}"' in core


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
    missing = [f"{script}:{lineno}: {line}" for lineno, line in _common_param_assignments(script) if "#" not in line]
    assert not missing, "Add inline comments to common parameters:\n" + "\n".join(missing)


def test_genome_annotation_core_quotes_known_path_sensitive_options():
    script = CORE_DIR / "gg_genome_annotation_core.sh"
    text = _read_text(script)

    banned_tokens = [
        "--lineage_dataset ${dir_busco_lineage}",
        "--download_path ${dir_busco_db}",
        'seqkit seq --threads ${GG_TASK_CPUS} ${file_sp_genome} > "busco_genome_input.fa"',
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
        "awk -v scaffold=\"${scaffold}\" '$1 == scaffold' tmp.species1.bed",
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
        "trimal -in ${dir_single_copy_mafft}/${infile}",
        "trimal -in ${dir_busco_mafft}/${infile}",
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
    assert 'if ! awk -v gene="${gene_name}" \'' in text
    assert 'sub(/[[:space:]].*$/, "", header)' in text
    assert "if (header == gene)" in text


def test_gene_evolution_core_filters_empty_translated_records_before_diamond_makedb():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert "filter_translated_fasta_for_diamond() {" in text
    assert 'gsub(/\\*/, "", $0)' in text
    assert 'if ($0 != "") {' in text
    assert (
        'printf("Dropped %d translated protein records with empty sequence after stop-codon removal.\\n", dropped) > "/dev/stderr"'
        in text
    )
    assert text.count("filter_translated_fasta_for_diamond \\") >= 2


def test_gene_evolution_core_quotes_orthogroup_lookup_and_makeblastdb_args():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert (
        "og_id=$(python -c \"import sys,pandas; df=pandas.read_csv(sys.argv[1],sep='\\t',header=0); print(df.loc[int(sys.argv[2]),:].iloc[0])\" ${file_orthogroup_genecount_selected} ${ind})"
        not in text
    )
    assert (
        'og_id=$(python -c "import sys,pandas; df=pandas.read_csv(sys.argv[1],sep=\'\\t\',header=0); print(df.loc[int(sys.argv[2]),:].iloc[0])" "${file_orthogroup_genecount_selected}" "${ind}")'
        not in text
    )
    assert "df=pandas.read_csv(sys.argv[1],sep='\\t',header=0); print(df.loc[int(sys.argv[2]),:].iloc[0])" not in text
    assert (
        "og_id=$(awk -F'\\t' -v row=\"${GG_ARRAY_TASK_ID}\" 'NR == (row + 1) { print $1; exit }' \"${file_orthogroup_genecount_selected}\")"
        in text
    )
    assert "makeblastdb -dbtype nucl -title ${sp_cds} -out ${sp_cds_blastdb}" not in text
    assert "makeblastdb -dbtype nucl -in ${sp_cds} -out ${sp_cds_blastdb}" not in text
    assert 'makeblastdb -dbtype nucl -title "${sp_cds}" -out "${sp_cds_blastdb}"' in text
    assert 'makeblastdb -dbtype nucl -in "${sp_cds}" -out "${sp_cds_blastdb}"' in text


def test_no_cp_out_or_mv_out_glob_arguments_in_core_scripts():
    pattern = re.compile(r"^[ \t]*(cp_out|mv_out)\b[^\n]*\*", re.MULTILINE)
    for script in sorted(CORE_DIR.glob("*.sh")):
        text = _read_text(script)
        assert pattern.search(text) is None, f"Use nullglob+array guard instead of cp_out/mv_out glob in {script}"


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
            assert pattern.search(text) is None, f"Found unquoted model/codon option in {script}: {pattern.pattern}"


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
        'busco_notung ${infile} ${dir_busco_iqtree_dna} "${dir_busco_notung_dna}" &',
        'busco_notung ${infile} ${dir_busco_iqtree_pep} "${dir_busco_notung_pep}" &',
        'busco_species_tree_assisted_gene_tree_rooting ${infile} "${dir_busco_notung_dna}" ${dir_busco_iqtree_dna} "${dir_busco_rooted_txt_dna}" "${dir_busco_rooted_nwk_dna}" &',
        'busco_species_tree_assisted_gene_tree_rooting ${infile} "${dir_busco_notung_pep}" ${dir_busco_iqtree_pep} "${dir_busco_rooted_txt_pep}" "${dir_busco_rooted_nwk_pep}" &',
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
        'param_species_tree="-s ${species_tree}"',
        'if [[ -n "${param_species_tree}" ]]; then',
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


def test_genome_evolution_core_requires_requested_species_tree_before_orthofinder():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)

    assert "species_tree_summary_generation_requested() {" in text
    assert "species_tree_requested_for_orthofinder=0" in text
    assert "if species_tree_summary_generation_requested; then" in text
    assert 'refresh_species_tree_for_shared_protein_input_signature "${shared_protein_input_signature}"' in text
    assert "Refusing to run OrthoFinder without a species tree." in text
    assert "Species-tree generation was requested, but no summary tree is available." in text
    assert "Running OrthoFinder without species tree constraints." not in text
    assert (
        "Refusing to run OrthoFinder without species tree constraints because a species tree was found but does not match the current OrthoFinder species set."
        in text
    )


def test_genome_evolution_core_only_uses_orthofinder_core_tree_when_species_tree_is_available():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)

    assert "orthofinder_core_species_tree_args=()" in text
    species_tree_arg = 'select_orthofinder_core_args+=(--species-tree "${species_tree}")'
    core_tree_check = 'if [[ ! -s "${file_orthofinder_core_species_tree}" ]]; then'
    core_arg = 'orthofinder_core_species_tree_args=(-s "${file_orthofinder_core_species_tree}")'
    orthofinder_call_arg = '"${orthofinder_core_species_tree_args[@]}"; then'
    assert "nwkit prune --invert_match yes" not in text
    assert species_tree_arg in text
    assert core_tree_check in text
    assert core_arg in text
    assert orthofinder_call_arg in text
    assert (
        text.index(species_tree_arg)
        < text.index(core_tree_check)
        < text.index(core_arg)
        < text.index(orthofinder_call_arg)
    )


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
        assert long_option.search(text) is None, f"Found unquoted variable in long option argument in {script}"
        assert short_option.search(text) is None, f"Found unquoted variable in short option argument in {script}"
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
        assert pattern.search(text) is None, f"Use arithmetic for-loop instead of for-in $(seq ...) in {script}"


def test_busco_support_script_uses_shared_hmmsearch_compat_wrapper():
    script = WORKFLOW_DIR / "support" / "gg_busco.sh"
    text = _read_text(script)
    assert "gg_busco_hmmsearch_wrapper_path() {" in text
    assert "gg_run_busco_with_metaeuk_modified_fas_compat() {" in text
    assert "gg_busco_stderr_matches_known_metaeuk_modified_fas_bug() {" in text
    assert 'GG_REAL_HMMSEARCH="${real_hmmsearch}" \\' in text
    assert "GG_BUSCO_METAEUK_MODIFIED_FAS_COMPAT=1 \\" in text


def test_busco_call_sites_use_shared_hmmsearch_compat_helper():
    transcriptome = _read_text(CORE_DIR / "gg_transcriptome_generation_core.sh")
    annotation = _read_text(CORE_DIR / "gg_genome_annotation_core.sh")
    input_generation = _read_text(CORE_DIR / "gg_input_generation_core.sh")
    genome_evolution = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")
    assert "gg_run_busco_with_metaeuk_modified_fas_compat \\" in transcriptome
    assert "gg_run_busco_with_metaeuk_modified_fas_compat \\" in annotation
    assert "gg_run_busco_with_metaeuk_modified_fas_compat \\" in input_generation
    assert "gg_run_busco_with_metaeuk_modified_fas_compat \\" in genome_evolution


def test_genome_evolution_core_uses_array_args_for_nwkit_mcmctree_constraints():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        'bound_params="${bound_params} --lower_bound ${mcmctree_params[2]}"',
        'bound_params="${bound_params} --upper_bound ${mcmctree_params[3]}"',
        'left_right="--left_species ${mcmctree_params[0]} --right_species ${mcmctree_params[1]}"',
        "tree_string=$(printf '%s\\n' \"${tree_string}\" | nwkit mcmctree ${left_right} ${bound_params})",
        'echo -e "${tree_string}" > "tmp.constrained.tree.nwk"',
        'nwkit mcmctree \\n    --infile "${file_undated_species_tree}"',
    ]
    for token in banned_tokens:
        assert token not in text, f"Found string-concatenated mcmctree args token: {token}"

    expected_tokens = [
        'dir_nwkit_download_dir="${gg_workspace_downloads_dir}/nwkit_downloads"',
        'ensure_dir "${dir_nwkit_download_dir}"',
        '--download_dir "${dir_nwkit_download_dir}"',
        "nwkit_args=(",
        '--download_dir "${dir_nwkit_download_dir}"',
        '--left_species "${mcmctree_params[0]}"',
        '--right_species "${mcmctree_params[1]}"',
        'nwkit_args+=(--lower_bound "${mcmctree_params[2]}")',
        'nwkit_args+=(--upper_bound "${mcmctree_params[3]}")',
        'echo "nwkit mcmctree params: ${nwkit_args[*]}"',
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
        'iqtree_mem_arg="-mem ${GG_MEM_TOTAL_GB}G"',
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


def test_gene_evolution_core_builds_species_foreground_regexes_with_label_boundaries():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)

    assert "write_species_trait_foreground_regex_table()" in text
    assert "from species_labeling import extract_species_label" in text
    assert "RANK_OR_QUALIFIER_TOKENS" in text
    assert 'return r"^{}_(?!(?:{})(?:\\.|_)).*"' in text
    assert "write_species_trait_foreground_regex_table species_trait_binary.tsv foreground.tsv" in text
    assert 'write_species_trait_foreground_regex_table "${file_sp_trait}" "foreground.tsv"' in text


def test_gene_evolution_core_anchors_query_ids_in_maxalign_keep_regex():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)

    assert 'patterns = [f"(?i:{re.escape(gene)}.*)" for gene in normalized_ids]' not in text
    assert 'patterns = [f"(?i:^{re.escape(gene)}$)" for gene in normalized_ids]' in text


def test_gene_evolution_core_disables_initial_ufboot_when_fast_mode_is_enabled():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert "if [[ ${num_seq} -gt ${iqtree_fast_mode_gt} ]]; then" in text
    assert "if [[ ${use_ufboot} -eq 1 ]]; then" in text
    assert (
        "Disabling IQ-TREE UFBOOT because fast mode is enabled for large alignments (${num_seq} > ${iqtree_fast_mode_gt})."
        in text
    )
    assert "other_iqtree_params=()" in text
    assert 'file_tree="${og_id}.treefile"' in text
    assert "other_iqtree_params+=(--fast)" in text


def test_gene_evolution_core_keeps_generax_ufboot_task_free_of_fast_flag():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert "Skipping IQ-TREE --fast in UFBOOT-on-GeneRax mode because the options are incompatible." in text

    ufboot_block_start = text.index('task="IQ-TREE UFBOOT on GeneRax topology"')
    ufboot_block_end = text.index("build_iqtree_mem_args", ufboot_block_start)
    ufboot_block = text[ufboot_block_start:ufboot_block_end]
    assert "other_iqtree_params+=( --fast )" not in ufboot_block


def test_gene_evolution_core_drops_all_branch_lengths_from_generax_constraint_tree():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    ufboot_block_start = text.index('task="IQ-TREE UFBOOT on GeneRax topology"')
    ufboot_block_end = text.index("build_iqtree_mem_args", ufboot_block_start)
    ufboot_block = text[ufboot_block_start:ufboot_block_end]

    assert (
        'nwkit drop --target all --length yes --outformat 9 --outfile "${og_id}.generax_ufboot.constraint.nwk"'
        in ufboot_block
    )
    assert (
        'nwkit drop --target root --length yes --outfile "${og_id}.generax_ufboot.constraint.nwk"' not in ufboot_block
    )


def test_gene_evolution_core_uses_container_safe_generax_mpi_launcher():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    generax_block_start = text.index('task="GeneRax"')
    generax_block_end = text.index('task="IQ-TREE UFBOOT on GeneRax topology"', generax_block_start)
    generax_block = text[generax_block_start:generax_block_end]

    assert 'mpiexec_args=(mpiexec -oversubscribe -np "${GG_TASK_CPUS}")' in generax_block
    assert (
        "mpi_env_args=(env OMPI_MCA_plm=isolated OMPI_MCA_plm_rsh_agent=/bin/false OMPI_MCA_btl=^openib)"
        in generax_block
    )
    assert "running_under_scheduler" not in generax_block


def test_generax_container_smoke_test_uses_runtime_mpi_launcher():
    script = CONTAINER_SCRIPTS_DIR / "ensure_generax_stable.sh"
    body = _function_body(_read_text(script), "run_smoke_test")

    assert "mpi_env_args=(env OMPI_MCA_plm=isolated OMPI_MCA_plm_rsh_agent=/bin/false OMPI_MCA_btl=^openib)" in body
    assert "mpiexec_args=(mpiexec --allow-run-as-root -oversubscribe -np 1)" in body
    assert '"${mpi_env_args[@]}" "${mpiexec_args[@]}"' in body
    assert "env OMPI_MCA_plm=isolated \\" not in body


def test_container_ghcr_builds_arm64_on_native_runner_without_qemu():
    workflow = _read_text(GITHUB_WORKFLOWS_DIR / "container-ghcr.yml")
    build_start = workflow.index("  build-and-push:")
    publish_start = workflow.index("  publish-manifest:", build_start)
    build_block = workflow[build_start:publish_start]

    assert "runs-on: ${{ matrix.runner }}" in build_block
    assert "- platform: linux/amd64\n            runner: ubuntu-latest" in build_block
    assert "- platform: linux/arm64\n            runner: ubuntu-24.04-arm" in build_block
    assert "docker/setup-qemu-action" not in build_block


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
            assert pattern.search(line) is None, f"Avoid pipe-to-head under pipefail in {script}:{line_no}: {line}"


def test_no_pipe_to_awk_exit_in_core_and_support_scripts():
    scripts = sorted(CORE_DIR.glob("*.sh")) + sorted((WORKFLOW_DIR / "support").glob("*.sh"))
    pattern = re.compile(r"\|\s*awk\b[^\n]*\bexit\b")
    for script in scripts:
        for line_no, line in enumerate(_read_text(script).splitlines(), start=1):
            stripped = line.lstrip()
            if stripped.startswith("#"):
                continue
            assert pattern.search(line) is None, f"Avoid pipe-to-awk-exit under pipefail in {script}:{line_no}: {line}"


def test_no_standalone_unquoted_file_or_dir_variable_argument_lines_in_core_scripts():
    pattern = re.compile(r"^[ \t]*\$\{(?:file|dir)_[^}]+\}[ \t]*\\?$", re.MULTILINE)
    for script in sorted(CORE_DIR.glob("*.sh")):
        text = _read_text(script)
        assert pattern.search(text) is None, f"Found standalone unquoted file/dir variable argument line in {script}"


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
    assert (
        "if ! awk -F '\\t' -v gene=\"${gene}\" '$1 == gene {found=1; exit} END {exit(found ? 0 : 1)}' \"${og_id}.rpsblast.tmp.tsv\"; then"
        in text
    )
    assert (
        "qlen=$(seqkit fx2tab --length ungapped_translated_cds.fas | awk -F '\\t' -v gene=\"${gene}\" '$1 == gene {print $NF}')"
        in text
    )


def test_gene_evolution_core_tree_pruning_uses_awk_id_match_instead_of_python_split_filter():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "python -c \"import sys,re; keys = open(sys.argv[1]).read().split('\\n'); entries = open(sys.argv[2]).read().split('>'); [ sys.stdout.write('>'+e) for e in entries if (re.sub('\\n.*','',e) in keys)&(len(e)!=0) ]\"",
    ]
    for token in banned_tokens:
        assert token not in text
    assert 'sub(/[[:space:]].*$/, "", id)' in text
    assert "keep_seq = (id in keep)" in text
    assert '\' target_genes.txt "${og_id}.untrimmed.input.fasta" > "${og_id}.untrimmed.pruned.tmp.fasta"' in text
    assert '\' target_genes.txt "${og_id}.trimmed.input.fasta" > "${og_id}.trimmed.pruned.tmp.fasta"' in text


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
    assert "root has a distance" in text
    assert "clear_root_branch_property_compat(tree, 'dist')" in text


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
    assert "csubst search \\" in core
    assert 'csubst_search_dir="csubst_search"' in core
    assert 'csubst_nonsyn_recode=$(echo "${csubst_nonsyn_recode:-${GG_COMMON_CSUBST_NONSYN_RECODE:-no}}" | tr' in core
    assert '--nonsyn_recode "${csubst_nonsyn_recode}"' in core
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
    assert "if recode != 'no':" in wrapper
    assert "cmd += ['--nonsyn_recode', recode]" in wrapper


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
        "echo '${run_collect_gff_info} is deactivated. Empty input:' ${dir_sp_gff}",
        "echo '${run_scm_intron} is deactivated. Empty input:' ${dir_sp_gff}",
        "echo '${dir_sp_expression} is not empty. Continued:' ${dir_sp_expression}",
        "echo '${dir_sp_expression} is empty:' ${dir_sp_expression}",
        "echo '${dir_sp_genome} is not empty. Continued:' ${dir_sp_genome}",
        "echo '${dir_sp_genome} is empty:' ${dir_sp_genome}",
    ]
    for token in banned_tokens:
        assert token not in text
    expected_tokens = [
        'echo "\\${run_collect_gff_info} is deactivated. Empty input: ${dir_sp_gff}"',
        'echo "\\${run_scm_intron} is deactivated. Empty input: ${dir_sp_gff}"',
        'echo "\\${dir_sp_expression} is not empty. Continued: ${dir_sp_expression}"',
        'echo "\\${dir_sp_expression} is empty: ${dir_sp_expression}"',
        'echo "\\${dir_sp_genome} is not empty. Continued: ${dir_sp_genome}"',
        'echo "\\${dir_sp_genome} is empty: ${dir_sp_genome}"',
    ]
    for token in expected_tokens:
        assert token in text


def test_genome_evolution_core_uses_filtered_orthofinder_core_sampling():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    assert "shuf -n" not in text
    assert "select_orthofinder_core_species.py" in text
    assert 'orthofinder_core_filters="${orthofinder_core_filters:-busco_complete_pct:ge:80,num_seq:le:100000}"' in text
    assert '--filters "${orthofinder_core_filters}"' in text
    assert '--rank "${orthofinder_core_rank}"' in text
    assert '--method "${orthofinder_core_method}"' in text


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

    assert "mapfile -t files < <(find \"${dir_genelist}\" -mindepth 1 -maxdepth 1 -type f ! -name '.*' | sort)" in text
    assert 'files=( "${dir_genelist}"/* )' not in text
    assert 'if [[ ! "${GG_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${GG_ARRAY_TASK_ID} -lt 1 ]]; then' in text
    assert "num_orthogroups=$(awk 'END { print (NR > 0 ? NR - 1 : 0) }'" in text
    assert "if [[ ${GG_ARRAY_TASK_ID} -gt ${num_orthogroups} ]]; then" in text
    assert "df=pandas.read_csv(sys.argv[1],sep='\\t',header=0); print(df.loc[int(sys.argv[2]),:].iloc[0])" not in text
    assert "og_id=$(awk -F'\\t' -v row=\"${GG_ARRAY_TASK_ID}\" 'NR == (row + 1) { print $1; exit }'" in text

    idx_guard = "if [[ ${idx} -ge ${#files[@]} ]]; then"
    idx_use = 'file_query_gene="${files[${idx}]}"'
    assert idx_guard in text
    assert idx_use in text
    assert text.index(idx_guard) < text.index(idx_use)


def test_gene_evolution_core_escapes_embedded_quotes_in_seqtype_error_message():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert (
        'echo "Unsupported sequence type \'${seqtype}\' in \'${file_query_gene}\'. Only "DNA" or "Protein" are allowed. Exiting."'
        not in text
    )
    assert (
        'echo "Unsupported sequence type \'${seqtype}\' in \'${file_query_gene}\'. Only \\"DNA\\" or \\"Protein\\" are allowed. Exiting."'
        in text
    )


def test_genome_annotation_core_guards_array_task_id_before_task_index_math():
    script = CORE_DIR / "gg_genome_annotation_core.sh"
    text = _read_text(script)
    assert 'if [[ ! -d "${dir_sp_cds}" ]]; then' in text
    assert 'echo "Input directory not found: ${dir_sp_cds}. Exiting."' in text
    assert "find \"${dir_sp_cds}\" -maxdepth 1 -type f ! -name '.*'" in text
    assert "find \"${dir_sp_dnaseq}/${sp_ub}\" -type f ! -name '.*'" in text
    guard = 'if [[ ! "${GG_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${GG_ARRAY_TASK_ID} -lt 1 ]]; then'
    task_index = "task_index=$((GG_ARRAY_TASK_ID - 1))"
    assert guard in text
    assert task_index in text
    assert text.index(guard) < text.index(task_index)


def test_genome_annotation_core_multispecies_summary_requires_real_summary_inputs():
    script = CORE_DIR / "gg_genome_annotation_core.sh"
    text = _read_text(script)
    assert "summary_inputs_available=0" in text
    assert "No multispecies summary inputs are available yet. Skipping summary generation." in text
    assert 'if is_output_older_than_inputs "^file_sp_" "${file_multispecies_summary}"; then' not in text
    assert 'if is_output_older_than_inputs "^(dir_summary_|file_summary_)" "${file_multispecies_summary}"; then' in text
    assert 'dir_summary_species_annotation="${gg_workspace_output_dir}/species_cds_annotation"' in text
    assert 'file_summary_species_tree_dated="${dir_summary_species_tree}/dated_species_tree.nwk"' in text


def test_genome_evolution_core_excludes_hidden_files_when_listing_species_proteins():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    assert "find \"${dir_sp_protein_orthofinder}\" -maxdepth 1 -type f ! -name '.*'" in text


def test_input_generation_core_excludes_hidden_files_when_listing_species_inputs():
    script = CORE_DIR / "gg_input_generation_core.sh"
    text = _read_text(script)
    assert "find \"${search_dir}\" -maxdepth 1 -type f ! -name '.*'" in text
    expected_tokens = [
        'list_nonhidden_matching_files "${species_cds_dir}"',
        'list_nonhidden_matching_files "${species_gff_dir}"',
    ]
    for token in expected_tokens:
        assert text.count(token) >= 2


def test_input_generation_core_runs_cds_gff_mapping_validation():
    script = CORE_DIR / "gg_input_generation_core.sh"
    text = _read_text(script)
    assert "validate_cds_gff_mapping.py" in text
    assert '--species-cds-dir "${species_cds_dir}"' in text
    assert '--species-gff-dir "${species_gff_dir}"' in text
    assert '--nthreads "${GG_TASK_CPUS:-1}"' in text


def test_input_generation_core_populates_species_summary_taxonomy_metadata_nonfatally():
    script = CORE_DIR / "gg_input_generation_core.sh"
    text = _read_text(script)
    assert 'if ! ensure_ete_taxonomy_db "${gg_workspace_dir}"; then' in text
    assert "species_summary taxonomy metadata" in text
    assert "Continuing without taxid/genetic code annotation." in text


def test_input_generation_core_wires_array_modes_and_busco_outputs_under_input_generation():
    script = CORE_DIR / "gg_input_generation_core.sh"
    text = _read_text(script)
    assert 'input_generation_mode="${input_generation_mode:-single}"' in text
    assert "single|array_prepare|array_worker|array_finalize" in text
    assert 'run_species_busco="${run_species_busco:-1}"' in text
    assert 'run_cds_fx2tab="${run_cds_fx2tab:-1}"' in text
    assert 'run_multispecies_summary="${run_multispecies_summary:-1}"' in text
    assert 'species_cds_fx2tab_dir="${input_generation_root}/species_cds_fx2tab"' in text
    assert 'species_busco_full_dir="${input_generation_root}/species_cds_busco_full"' in text
    assert 'species_busco_short_dir="${input_generation_root}/species_cds_busco_short"' in text
    assert 'file_multispecies_summary="${input_generation_root}/annotation_summary/annotation_summary.tsv"' in text
    assert 'cmd=(Rscript "${gg_support_dir}/annotation_summary.r")' in text
    assert 'cmd+=(--dir_species_cds_fx2tab="${species_cds_fx2tab_dir}")' in text
    assert '--dir_species_tree="${gg_workspace_output_dir}/species_tree"' not in text
    assert 'python "${gg_support_dir}/plan_input_generation_tasks.py"' in text
    assert 'python "${gg_support_dir}/run_input_generation_task.py"' in text
    assert 'python "${gg_support_dir}/merge_input_generation_shards.py"' in text


def test_mmseqs_uniref90_download_retries_and_reports_disk_context():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "_download_mmseqs_uniref90_db")
    assert 'local output_db="${db_dir}/${uniref_db}_DB"' in body
    assert 'mmseqs databases "${uniref_db}" "${output_db}" "${db_dir}" --threads "${nthreads}"' in body
    assert "for attempt in 1 2 3; do" in body
    assert "Preparing MMseqs2 UniRef90 taxonomy DB in: ${db_dir} (attempt ${attempt}/${max_attempts})" in body
    assert "MMseqs2 UniRef90 taxonomy DB preparation failed in: ${db_dir} (attempt ${attempt}/${max_attempts})" in body
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
        'migrate_legacy_fasta_outputs "${dir_busco_fasta}"',
        'migrate_legacy_fasta_outputs "${dir_busco_mafft}"',
        'migrate_legacy_fasta_outputs "${dir_busco_trimal}"',
        '-name "*_busco.full.tsv"',
        '-name "*_busco.short.txt"',
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
