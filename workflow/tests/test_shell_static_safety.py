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


def _set_e_scripts():
    scripts = []
    for script in _workflow_shell_scripts():
        text = _read_text(script)
        header = "\n".join(text.splitlines()[:30])
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
        text = _read_text(script)
        header = "\n".join(text.splitlines()[:20])
        has_pipefail = ("set -eo pipefail" in header) or ("set -euo pipefail" in header)
        assert has_pipefail, f"Missing pipefail guard in script header: {script}"


def test_core_and_entrypoint_scripts_use_strict_euo_pipefail():
    scripts = _core_and_entrypoint_scripts()
    assert scripts, "No core/entrypoint scripts were found."
    for script in scripts:
        text = _read_text(script)
        header = "\n".join(text.splitlines()[:20])
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
        text = _read_text(script)
        header = "\n".join(text.splitlines()[:30])
        assert "set -euo pipefail" in header, f"Use strict mode (set -euo pipefail): {script}"


def test_busco_download_lock_is_global():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    assert 'lock_file="${runtime_busco_db}/locks/busco_downloads.lock"' in text
    assert '${runtime_busco_db}/${busco_lineage}.lock' not in text


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
    assert '< "${dir_script}/core/gg_progress_summary_core.sh"' in text
    assert "orthogroup_output_summary.py" not in text
    assert "transcriptome_assembly_output_summary.py" not in text


def test_no_known_unquoted_query_gene_file_expansions():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "head --bytes 1 ${file_query_gene}",
        "seqkit stats --tabular ${file_query_gene}",
        "seqkit translate --allow-unknown-codon --transl-table ${genetic_code} --threads ${NSLOTS} ${file_query_gene}",
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
        "python ${dir_script}/",
        "Rscript ${dir_script}/",
        "bash ${dir_script}/",
    ]
    for script in sorted(CORE_DIR.glob("*.sh")):
        text = _read_text(script)
        for token in banned_tokens:
            assert token not in text, f"Found unquoted dir_script invocation in {script}: {token}"


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
        assert "set_singularity_command" in text, f"Missing set_singularity_command call: {script}"


def test_gg_trigger_versions_dump_is_runtime_agnostic():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "gg_trigger_versions_dump")
    assert 'export SINGULARITYENV_GG_VERSION="${gg_version}"' in body
    assert 'export APPTAINERENV_GG_VERSION="${gg_version}"' in body
    assert 'command -v "${container_runtime_bin}"' in body
    assert '"${container_runtime_bin}" inspect "${gg_image}"' in body
    assert '"${container_runtime_bin}" version || {' in body
    assert 'singularity inspect "${gg_image}"' not in body
    assert 'singularity version || {' not in body


def test_progress_summary_entrypoint_uses_auto_forwarding_and_normalized_nslots():
    entrypoint = WORKFLOW_DIR / "gg_progress_summary_entrypoint.sh"
    text = _read_text(entrypoint)

    assert "forward_config_vars_to_container_env()" in text
    assert 'forward_config_vars_to_container_env "${BASH_SOURCE[0]}"' in text
    assert "unset -f forward_config_vars_to_container_env" in text
    assert "for exported_name in mode_transcriptome_assembly ncpu_progress_summary; do" not in text
    assert 'ncpu_progress_summary="${ncpu_progress_summary:-${NSLOTS:-1}}"' not in text

    idx_variable_sgenizer = text.index("variable_SGEnizer")
    idx_ncpu_default = text.index(': "${ncpu_progress_summary:=${NSLOTS:-1}}"')
    assert idx_ncpu_default > idx_variable_sgenizer


def test_input_generation_entrypoint_forwards_env_driven_overrides():
    entrypoint = WORKFLOW_DIR / "gg_input_generation_entrypoint.sh"
    text = _read_text(entrypoint)

    assert "for gg_input_var_name in ${!GG_INPUT_@}; do" in text
    assert 'export "SINGULARITYENV_${gg_input_var_name}=${!gg_input_var_name}"' in text
    assert 'export "APPTAINERENV_${gg_input_var_name}=${!gg_input_var_name}"' in text
    assert 'export "SINGULARITYENV_GG_INPUTPREP_CONFIG=${GG_INPUTPREP_CONFIG}"' in text
    assert 'export "APPTAINERENV_GG_INPUTPREP_CONFIG=${GG_INPUTPREP_CONFIG}"' in text
    assert 'export "SINGULARITYENV_GG_DATASET_ROOT=${GG_DATASET_ROOT}"' in text
    assert 'export "APPTAINERENV_GG_DATASET_ROOT=${GG_DATASET_ROOT}"' in text


def test_ensure_latest_jaspar_file_uses_set_e_safe_assignments():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "ensure_latest_jaspar_file")

    assert "if resolved_filename=$(_jaspar_find_latest_meme_filename_remote); then" in body
    assert "if resolved_filename=$(_jaspar_find_latest_meme_filename_local" in body
    assert "if resolved_path=$(_ensure_jaspar_file_named" in body

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

    assert 'export "SINGULARITYENV_PYMOL_HEADLESS=${PYMOL_HEADLESS}"' in text
    assert 'export "APPTAINERENV_PYMOL_HEADLESS=${PYMOL_HEADLESS}"' in text
    assert 'export "SINGULARITYENV_QT_QPA_PLATFORM=${QT_QPA_PLATFORM}"' in text
    assert 'export "APPTAINERENV_QT_QPA_PLATFORM=${QT_QPA_PLATFORM}"' in text


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


def test_gene_evolution_core_clamps_l1ou_cpu_selection_to_available_cores():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert 'CPU_PER_HOST=$(grep -c processor /proc/cpuinfo)' not in text
    assert 'cpu_pick="${NSLOTS}"' in text
    assert 'if [[ "${cpu_pick}" -gt "${CPU_PER_HOST}" ]]; then' in text
    assert 'cpu_id=$(python -c' in text
    assert '"${cpu_pick}" "${CPU_PER_HOST}"' in text
    assert "if command -v taskset >/dev/null 2>&1 && [[ -n \"${cpu_id}\" ]]; then" in text
    assert 'taskset -c "${cpu_id}" "${l1ou_cmd[@]}"' in text
    assert '"${l1ou_cmd[@]}"' in text


def test_gene_evolution_core_uses_explicit_ne_and_grouped_logic_for_tree_pruning_gate():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    assert "if [[ ! ${run_tree_pruning} -eq 1 && ${run_phylogeneticem} -eq 1 || ${run_l1ou} -eq 1 ]]; then" not in text
    assert "if [[ ${run_tree_pruning} -ne 1 && ( ${run_phylogeneticem} -eq 1 || ${run_l1ou} -eq 1 ) ]]; then" in text


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


def test_gg_trigger_versions_dump_reuses_one_log_per_container():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "gg_trigger_versions_dump")
    assert 'log_file="${versions_dir}/container.${container_key_hash}.versions.log"' in body
    assert 'if [[ -s "${log_file}" ]]; then' in body
    assert 'gg_trigger_versions_dump: skipped existing ${log_file}' in body
    assert 'timestamp=$(date' not in body


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


def test_entrypoints_allow_env_override_for_tmp_cleanup_flags():
    expected_tokens = {
        "gg_genome_annotation_entrypoint.sh": 'delete_tmp_dir="${delete_tmp_dir:-1}"',
        "gg_transcriptome_generation_entrypoint.sh": 'delete_tmp_dir="${delete_tmp_dir:-1}"',
        "gg_gene_evolution_entrypoint.sh": 'delete_tmp_dir="${delete_tmp_dir:-1}"',
        "gg_gene_evolution_entrypoint.sh#preexisting": 'delete_preexisting_tmp_dir="${delete_preexisting_tmp_dir:-1}"',
        "gg_genome_evolution_entrypoint.sh": 'delete_tmp_dir="${delete_tmp_dir:-1}"',
    }

    for key, token in expected_tokens.items():
        script_name = key.split("#")[0]
        text = _read_text(WORKFLOW_DIR / script_name)
        assert token in text, f"Missing env-overridable tmp cleanup flag in {script_name}: {token}"


def test_entrypoints_forward_cleanup_flags_defined_outside_config_block():
    expected_tokens = {
        "gg_genome_annotation_entrypoint.sh": 'export "SINGULARITYENV_delete_tmp_dir=${delete_tmp_dir}"',
        "gg_transcriptome_generation_entrypoint.sh": 'export "SINGULARITYENV_delete_tmp_dir=${delete_tmp_dir}"',
        "gg_gene_evolution_entrypoint.sh": 'export "SINGULARITYENV_delete_tmp_dir=${delete_tmp_dir}"',
        "gg_gene_evolution_entrypoint.sh#preexisting": 'export "SINGULARITYENV_delete_preexisting_tmp_dir=${delete_preexisting_tmp_dir}"',
    }
    for key, token in expected_tokens.items():
        script_name = key.split("#")[0]
        text = _read_text(WORKFLOW_DIR / script_name)
        assert token in text, f"Missing cleanup var forwarding in {script_name}: {token}"


def test_gene_convergence_entrypoint_does_not_define_unused_delete_tmp_dir():
    text = _read_text(WORKFLOW_DIR / "gg_gene_convergence_entrypoint.sh")
    assert "delete_tmp_dir=" not in text


def test_entrypoints_with_exit_if_running_call_duplicate_guard():
    scripts = [
        "gg_gene_convergence_entrypoint.sh",
        "gg_genome_annotation_entrypoint.sh",
        "gg_gene_evolution_entrypoint.sh",
    ]
    for script_name in scripts:
        text = _read_text(WORKFLOW_DIR / script_name)
        assert "exit_if_running=" in text
        assert "exit_if_running_qstat" in text, f"Missing duplicate-job guard call in {script_name}"


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


def test_pfam_helpers_use_only_new_runtime_layout_and_function_name():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    gene_core = CORE_DIR / "gg_gene_evolution_core.sh"
    util_text = _read_text(util_path)
    gene_text = _read_text(gene_core)

    assert "legacy_runtime_dir" not in util_text
    assert "downloads/Pfam_LE" not in util_text
    assert "ensure_pfam_domain_db()" not in util_text
    assert "ensure_pfam_domain_db" not in gene_text
    assert 'ensure_pfam_le_db "${dir_pg}"' in gene_text


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


def test_busco_dataset_move_uses_mv_option_separator():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "ensure_busco_download_path")
    assert "-exec mv -f -- {}" in body
    assert "-exec mv -f {}" not in body


def test_species_busco_set_check_is_gated_by_run_flag():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    gate = "if [[ ${run_species_get_busco_summary} -eq 1 ]]; then"
    normalize = 'normalize_busco_table_naming "${dir_species_busco_full}" "${dir_species_busco_short}"'
    check = 'if ! is_species_set_identical "${dir_sp_cds}" "${dir_species_busco_full}"; then'
    assert gate in text
    assert normalize in text
    assert check in text
    assert text.index(gate) < text.index(normalize) < text.index(check)


def test_genome_busco_summary_normalizes_legacy_table_names_before_collect_common():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    gate = "if [[ ${run_genome_get_busco_summary} -eq 1 ]]; then"
    normalize = 'normalize_busco_table_naming "${dir_species_busco_full}" "${dir_species_busco_short}"'
    collect = 'python "${dir_script}/collect_common_BUSCO_genes.py" \\'
    assert gate in text
    assert normalize in text
    assert collect in text
    gate_idx = text.index(gate)
    normalize_idx = text.index(normalize, gate_idx)
    collect_idx = text.index(collect, normalize_idx)
    assert gate_idx < normalize_idx < collect_idx


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
    assert "--cores ${NSLOTS}" not in text
    assert "--output_prefix ${dir_cafe_output}" not in text

    assert 'busco_grampa "${dir_busco_rooted_nwk_dna}" "$(dirname "${file_busco_grampa_dna}")" "${file_busco_grampa_dna}"' in text
    assert 'busco_grampa "${dir_busco_rooted_nwk_pep}" "$(dirname "${file_busco_grampa_pep}")" "${file_busco_grampa_pep}"' in text
    assert 'busco_grampa "./tmp.orthogroup_grampa_indir" "$(dirname "${file_orthogroup_grampa}")" "${file_orthogroup_grampa}"' in text
    assert '--genecount "${file_orthogroup_genecount_selected}"' in text
    assert '--dated_species_tree "${file_dated_species_tree}"' in text
    assert '--max_size_differential "${max_size_differential_cafe}"' in text
    assert '--tree "${file_dated_species_tree}"' in text
    assert '--n_gamma_cats "${n_gamma_cats_cafe}"' in text
    assert '--cores "${NSLOTS}"' in text
    assert '--output_prefix "${dir_cafe_output}"' in text


def test_genome_evolution_core_builds_grampa_arguments_with_array():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    assert "h1_param=\"-h1 " not in text
    assert "grampa_args=(" in text
    assert "grampa_args+=( -h1 \"${grampa_h1_normalized}\" )" in text
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
    for script in entrypoints:
        text = _read_text(script)
        mkdir_token = 'mkdir -p "${dir_pg}"'
        cd_token = 'cd "${dir_pg}"'
        assert mkdir_token in text, f"Missing workspace mkdir guard in {script}"
        assert cd_token in text, f"Missing workspace cd in {script}"
        assert text.index(mkdir_token) < text.index(cd_token), (
            f"Workspace mkdir must come before cd in {script}"
        )


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
        "--ncpu ${NSLOTS}; then",
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
        '--ncpu "${NSLOTS}"; then',
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
    assert "-T ${NSLOTS}" not in iq2mc_block
    assert '-m "${nucleotide_model}"' in iq2mc_block
    assert '-te "${file_constrained_tree}"' in iq2mc_block
    assert '--mcmc-bds "${mcmc_birth_death_sampling}"' in iq2mc_block
    assert '--mcmc-clock "${mcmc_clock_model}"' in iq2mc_block
    assert '--mcmc-iter "${mcmc_burnin},${mcmc_sampfreq},${mcmc_nsample}"' in iq2mc_block
    assert '-T "${NSLOTS}"' in iq2mc_block

    busco_summary_start = text.index('python "${dir_script}/collect_common_BUSCO_genes.py" \\')
    busco_summary_end = text.index('--outfile "tmp.busco_summary_table.tsv"', busco_summary_start) + len('--outfile "tmp.busco_summary_table.tsv"')
    busco_summary_block = text[busco_summary_start:busco_summary_end]
    assert "--busco_outdir ${dir_species_busco_full}" not in busco_summary_block
    assert "--ncpu ${NSLOTS}" not in busco_summary_block
    assert '--busco_outdir "${dir_species_busco_full}"' in busco_summary_block
    assert '--ncpu "${NSLOTS}"' in busco_summary_block


def test_transcriptome_core_quotes_known_path_sensitive_options_and_symlinks():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)

    banned_tokens = [
        "--fastq_dir ${dir_species_fastq}",
        "--out_dir ${dir_tmp}",
        "--metadata ${file_amalgkit_metadata}",
        "ln -s ${dir_amalgkit_getfastq_sp} \"./getfastq\"",
        "ln -s ${silva_rrna_ref} ${file_reference_fasta_link}",
        "kallisto index --make-unique -i \"./index/${sp_ub}.idx\" ${file_reference_fasta_link}",
        "--fasta_file ${file_longestcds}",
        "--mmseqs2taxonomy_tsv ${file_longestcds_mmseqs2taxonomy}",
        "--fx2tab_tsv ${file_longestcds_fx2tab}",
        "--species_name ${sp_ub}",
        "--rank ${contamination_removal_rank}",
        "seqkit seq --threads ${NSLOTS} ${file_isoform} --out-file \"busco_infile_cdna.fa\"",
        "seqkit seq --threads ${NSLOTS} ${file_longestcds} --out-file \"busco_infile_cds.fa\"",
        "seqkit seq --threads ${NSLOTS} ${file_longestcds_contamination_removal_fasta} --out-file \"busco_infile_cds.fa\"",
        "--lineage_dataset ${dir_busco_lineage}",
        "--download_path ${dir_busco_db}",
        "if [[ -e ${file_kallisto_reference_fasta} ]]; then",
        "ln -s ${file_kallisto_reference_fasta} ${file_reference_fasta_link}",
        "ln -s ${dir_amalgkit_quant}/${sp_ub} ./quant",
        'grep -e "${sp_space}" "./metadata/metadata.tsv"',
        "d.loc[:,'scientific_name']='${sp_ub}'",
    ]
    for token in banned_tokens:
        assert token not in text, f"Found unquoted transcriptome token: {token}"

    expected_tokens = [
        '--fastq_dir "${dir_species_fastq}"',
        '--out_dir "${dir_tmp}"',
        '--metadata "${file_amalgkit_metadata}"',
        'ln -s "${dir_amalgkit_getfastq_sp}" "./getfastq"',
        'ln -s "${silva_rrna_ref}" "${file_reference_fasta_link}"',
        'kallisto index --make-unique -i "./index/${sp_ub}.idx" "${file_reference_fasta_link}"',
        '--fasta_file "${file_longestcds}"',
        '--mmseqs2taxonomy_tsv "${file_longestcds_mmseqs2taxonomy}"',
        '--fx2tab_tsv "${file_longestcds_fx2tab}"',
        '--species_name "${sp_ub}"',
        '--rank "${contamination_removal_rank}"',
        'seqkit seq --threads "${NSLOTS}" "${file_isoform}" --out-file "busco_infile_cdna.fa"',
        'seqkit seq --threads "${NSLOTS}" "${file_longestcds}" --out-file "busco_infile_cds.fa"',
        'seqkit seq --threads "${NSLOTS}" "${file_longestcds_contamination_removal_fasta}" --out-file "busco_infile_cds.fa"',
        '--lineage_dataset "${dir_busco_lineage}"',
        '--download_path "${dir_busco_db}"',
        'if [[ -e "${file_kallisto_reference_fasta}" ]]; then',
        'ln -s "${file_kallisto_reference_fasta}" "${file_reference_fasta_link}"',
        'ln -s "${dir_amalgkit_quant}/${sp_ub}" "./quant"',
        'grep -F -- "${sp_space}" "./metadata/metadata.tsv"',
        "d.loc[:,'scientific_name']=sys.argv[1]",
    ]
    for token in expected_tokens:
        assert token in text, f"Missing quoted transcriptome token: {token}"


def test_transcriptome_core_sraid_metadata_filter_handles_zero_match_explicitly():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    assert 'grep -F -- "${sp_space}" "./metadata/metadata.tsv" || true;' in text
    assert 'if [[ $(wc -l < "./metadata.tsv") -le 1 ]]; then' in text
    assert "No metadata rows matched species" in text


def test_genome_annotation_core_quotes_known_path_sensitive_options():
    script = CORE_DIR / "gg_genome_annotation_core.sh"
    text = _read_text(script)

    banned_tokens = [
        "--lineage_dataset ${dir_busco_lineage}",
        "--download_path ${dir_busco_db}",
        "seqkit seq --threads ${NSLOTS} ${file_sp_genome} > \"busco_genome_input.fa\"",
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
        'seqkit seq --threads "${NSLOTS}" "${file_sp_genome}" > "busco_genome_input.fa"',
        'mmseqs createdb "${file_sp_cds}" queryDB',
        '--fasta_file "${file_sp_cds}"',
        '--mmseqs2taxonomy_tsv "${file_sp_cds_mmseqs2taxonomy}"',
        '--fx2tab_tsv "${file_sp_cds_fx2tab}"',
        '--species_name "${sp_ub}"',
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
        "--phylogeneticem_tree ${file_og_pem_tree}",
        "--phylogeneticem_regime ${file_og_pem_regime}",
        "--phylogeneticem_leaf ${file_og_pem_leaf}",
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
        '--phylogeneticem_tree "${file_og_pem_tree}"',
        '--phylogeneticem_regime "${file_og_pem_regime}"',
        '--phylogeneticem_leaf "${file_og_pem_leaf}"',
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
    assert 'og_id=$(awk -F\'\\t\' -v row="${SGE_TASK_ID}" \'NR == (row + 1) { print $1; exit }\' "${file_orthogroup_genecount_selected}")' in text
    assert 'makeblastdb -dbtype nucl -title ${sp_cds} -out ${sp_cds_blastdb}' not in text
    assert 'makeblastdb -dbtype nucl -in ${sp_cds} -out ${sp_cds_blastdb}' not in text
    assert 'makeblastdb -dbtype nucl -title "${sp_cds}" -out "${sp_cds_blastdb}"' in text
    assert 'makeblastdb -dbtype nucl -in "${sp_cds}" -out "${sp_cds_blastdb}"' in text


def test_transcriptome_core_avoids_direct_mv_out_glob_for_getfastq_and_quant():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    assert 'mv_out "${dir_tmp}"/getfastq/* "${dir_amalgkit_getfastq_sp}"' not in text
    assert "mv_out ./quant/* \"${dir_amalgkit_quant}/${sp_ub}\"" not in text
    assert 'getfastq_outputs=( "${dir_tmp}"/getfastq/* )' in text
    assert 'mv_out "${getfastq_outputs[@]}" "${dir_amalgkit_getfastq_sp}"' in text
    assert "quant_outputs=( ./quant/* )" in text
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
        'param_species_tree=( -s "${species_tree}" )',
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
    assert "orthofinder_output_directory_cleanup ${dir_orthofinder} ${NSLOTS}" not in text
    assert "orthofinder_output_directory_cleanup ${dir_orthofinder}/core ${NSLOTS}" not in text
    assert 'orthofinder_output_directory_cleanup "${dir_orthofinder}" "${NSLOTS}"' in text
    assert 'orthofinder_output_directory_cleanup "${dir_orthofinder}/core" "${NSLOTS}"' in text


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
    assert 'query_hits_tmp_files=( "${query_hits_tmp_dir}"/*.hits.fasta )' in text
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
    assert 'input_files=( "${file_isoform}" )' in text
    assert 'input_files+=( "${file_longestcds}" )' in text
    assert 'input_files+=( "${file_longestcds_contamination_removal_fasta}" )' in text
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
        "if [[ ${nslots_assembly} -lt 1 ]]; then",
        "nslots_assembly=1",
        "if [[ ${memory_assembly} -lt 1 ]]; then",
        "memory_assembly=1",
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
    ]
    for token in banned_tokens:
        assert token not in text, f"Found string-concatenated mcmctree args token: {token}"

    expected_tokens = [
        "nwkit_args=(",
        '--left_species "${mcmctree_params[0]}"',
        '--right_species "${mcmctree_params[1]}"',
        'nwkit_args+=( --lower_bound "${mcmctree_params[2]}" )',
        'nwkit_args+=( --upper_bound "${mcmctree_params[3]}" )',
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
        "iqtree_mem_arg=\"-mem ${MEM_PER_HOST}G\"",
        "${iqtree_mem_arg} \\",
        "${bootstrap_params}; then",
    ]
    for token in banned_tokens:
        assert token not in text, f"Found fragile concat-IQ-TREE args token: {token}"

    expected_tokens = [
        "bootstrap_params=( --ufboot 1000 --bnni )",
        "bootstrap_params=()",
        "iqtree_mem_args=()",
        'iqtree_mem_args=( -mem "${MEM_PER_HOST}G" )',
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
        "other_iqtree_params=( --ufboot 1000 --bnni )",
        "other_iqtree_params+=( --fast )",
        '"${other_iqtree_params[@]}"',
        "foreground_params=( --foreground foreground.tsv --fg_format 2 )",
        "foreground_params=()",
        '"${foreground_params[@]}"',
    ]
    for token in expected_tokens:
        assert token in text, f"Missing array-based optional args token: {token}"


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
        '(os.path.exists(params["phylogeneticem_regime"]))&(os.path.exists(params["phylogeneticem_leaf"]))',
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


def test_is_fastq_requiring_downstream_analysis_done_quotes_path_checks():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "is_fastq_requiring_downstream_analysis_done")
    assert '-s ${file_isoform}' not in body
    assert '-s ${file_rRNA_contamination_report}' not in body
    assert '-s ${file_amalgkit_merge_count}' not in body
    assert '-s "${file_isoform}"' in body
    assert '-s "${file_rRNA_contamination_report}"' in body
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
    assert 'if [[ ! -d "${dir_busco}" || -z "$(find "${dir_busco}" -mindepth 1 -print -quit 2>/dev/null)" ]]; then' in text


def test_gene_evolution_core_quotes_key_s_checks_in_downstream_tasks():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)
    banned_tokens = [
        "if [[ -s ${file_og_expression} && ( ${run_l1ou} -eq 1 || ${run_phylogeneticem} -eq 1 ) ]]; then",
        "if [[ ! -s ${file_og_hyphy_relax_reversed} && ${run_hyphy_relax_reversed} -eq 1 ]]; then",
        "if [[ ! -s ${file_og_scm_intron_summary} && ${run_scm_intron} -eq 1 ]]; then",
        "if [[ ${run_phylogeneticem} -eq 1 && ( ! -s ${file_og_pem_rdata} || ! -s ${file_og_pem_tree} ||  ! -s ${file_og_pem_regime} || ! -s ${file_og_pem_leaf} ) && ${num_gene_after_maxalign:-0} -gt 3 ]]; then",
        "if [[ ( ! -s ${file_og_l1ou_fit_rdata} || ! -s ${file_og_l1ou_fit_tree} || ! -s ${file_og_l1ou_fit_regime} || ! -s ${file_og_l1ou_fit_leaf} ) && ${run_l1ou} -eq 1 ]]; then",
        "if ( [[ ${summary_flag} -eq 1 || ! -s ${file_og_tree_plot} ]] ) && [[ ${run_tree_plot} -eq 1 ]]; then",
        "if [[ -s ${file_og_stat_branch} && -s ${file_og_stat_tree} && -s ${file_og_tree_plot} && ${gg_debug_mode:-0} -eq 0 ]]; then",
    ]
    for token in banned_tokens:
        assert token not in text

    expected_tokens = [
        'if [[ -s "${file_og_expression}" && ( ${run_l1ou} -eq 1 || ${run_phylogeneticem} -eq 1 ) ]]; then',
        'if [[ ! -s "${file_og_hyphy_relax_reversed}" && ${run_hyphy_relax_reversed} -eq 1 ]]; then',
        'if [[ ! -s "${file_og_scm_intron_summary}" && ${run_scm_intron} -eq 1 ]]; then',
        'if [[ ${run_phylogeneticem} -eq 1 && ( ! -s "${file_og_pem_rdata}" || ! -s "${file_og_pem_tree}" || ! -s "${file_og_pem_regime}" || ! -s "${file_og_pem_leaf}" ) && ${num_gene_after_maxalign:-0} -gt 3 ]]; then',
        'if [[ ( ! -s "${file_og_l1ou_fit_rdata}" || ! -s "${file_og_l1ou_fit_tree}" || ! -s "${file_og_l1ou_fit_regime}" || ! -s "${file_og_l1ou_fit_leaf}" ) && ${run_l1ou} -eq 1 ]]; then',
        'if ( [[ ${summary_flag} -eq 1 || ! -s "${file_og_tree_plot}" ]] ) && [[ ${run_tree_plot} -eq 1 ]]; then',
        'if [[ -s "${file_og_stat_branch}" && -s "${file_og_stat_tree}" && -s "${file_og_tree_plot}" && ${gg_debug_mode:-0} -eq 0 ]]; then',
    ]
    for token in expected_tokens:
        assert token in text


def test_gene_evolution_core_guards_sge_task_id_range_before_input_indexing():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = _read_text(script)

    assert 'mapfile -t files < <(find "${dir_genelist}" -mindepth 1 -maxdepth 1 -type f ! -name \'.*\' | sort)' in text
    assert 'files=( "${dir_genelist}"/* )' not in text
    assert 'if [[ ! "${SGE_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${SGE_TASK_ID} -lt 1 ]]; then' in text
    assert 'num_orthogroups=$(awk \'END { print (NR > 0 ? NR - 1 : 0) }\'' in text
    assert 'if [[ ${SGE_TASK_ID} -gt ${num_orthogroups} ]]; then' in text
    assert "df=pandas.read_csv(sys.argv[1],sep='\\t',header=0); print(df.loc[int(sys.argv[2]),:].iloc[0])" not in text
    assert 'og_id=$(awk -F\'\\t\' -v row="${SGE_TASK_ID}" \'NR == (row + 1) { print $1; exit }\'' in text

    idx_guard = 'if [[ ${idx} -ge ${#files[@]} ]]; then'
    idx_use = 'file_query_gene="${files[${idx}]}"'
    assert idx_guard in text
    assert idx_use in text
    assert text.index(idx_guard) < text.index(idx_use)


def test_transcriptome_core_guards_sge_task_id_range_before_array_indexing():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)

    assert 'dir_input_fastq="${dir_transcriptome_assembly_input}/input_fastq"' not in text
    assert 'dir_input_sra_list="${dir_transcriptome_assembly_input}/input_sra_list"' not in text
    assert "Backward compatibility for legacy input layout" not in text
    assert 'mapfile -t fastq_mode_dirs < <(find "${dir_input_fastq}" -mindepth 1 -maxdepth 1 -type d ! -name \'.*\' | sort)' in text
    assert 'mapfile -t files_fastq < <(find "${dir_species_fastq}" -maxdepth 1 -type f ! -name \'.*\' \\( -name "*.fq" -o -name "*.fastq" -o -name "*.fq.gz" -o -name "*.fastq.gz" \\) | sort)' in text
    assert 'mapfile -t dirs < <(find "${dir_input_fastq}" -mindepth 1 -maxdepth 1 -type d ! -name \'.*\' | sort)' in text
    assert 'mapfile -t sra_mode_files < <(find "${dir_input_sra_list}" -mindepth 1 -maxdepth 1 -type f ! -name \'.*\' | sort)' in text
    assert 'mapfile -t metadata_mode_files < <(find "${dir_amalgkit_metadata}" -mindepth 1 -maxdepth 1 -type f ! -name \'.*\' | sort)' in text
    assert 'mapfile -t files < <(find "${dir_input_sra_list}" -mindepth 1 -maxdepth 1 -type f ! -name \'.*\' | sort)' in text
    assert 'mapfile -t files < <(find "${dir_amalgkit_metadata}" -mindepth 1 -maxdepth 1 -type f ! -name \'.*\' | sort)' in text
    assert re.search(r'^[ \t]*dirs=\( "\$\{dir_input_fastq\}"/\* \)', text, re.MULTILINE) is None
    assert re.search(r'^[ \t]*files=\( "\$\{dir_input_sra_list\}"/\* \)', text, re.MULTILINE) is None
    assert re.search(r'^[ \t]*files=\( "\$\{dir_amalgkit_metadata\}"/\* \)', text, re.MULTILINE) is None
    assert 'files_fastq=( "${dir_species_fastq}"/* )' not in text
    assert 'if [[ ! "${SGE_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${SGE_TASK_ID} -lt 1 ]]; then' in text
    assert 'if [[ ${#dirs[@]} -eq 0 ]]; then' in text

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


def test_transcriptome_core_guards_getfastq_outputs_before_rrna_assembly_and_quant():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    dir_guard = 'if [[ ! -d "${dir_amalgkit_getfastq_sp}" ]]; then'
    fastq_guard = 'if [[ -z "$(find "${dir_amalgkit_getfastq_sp}" -type f -name "*.amalgkit.fastq.gz" -print -quit 2>/dev/null)" ]]; then'
    assert text.count(dir_guard) >= 3
    assert text.count(fastq_guard) >= 3
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


def test_genome_annotation_core_guards_sge_task_id_before_task_index_math():
    script = CORE_DIR / "gg_genome_annotation_core.sh"
    text = _read_text(script)
    assert 'if [[ ! -d "${dir_sp_cds}" ]]; then' in text
    assert 'echo "Input directory not found: ${dir_sp_cds}. Exiting."' in text
    assert 'find "${dir_sp_cds}" -maxdepth 1 -type f ! -name \'.*\'' in text
    assert 'find "${dir_sp_dnaseq}/${sp_ub}" -type f ! -name \'.*\'' in text
    guard = 'if [[ ! "${SGE_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${SGE_TASK_ID} -lt 1 ]]; then'
    task_index = "task_index=$((SGE_TASK_ID-1))"
    assert guard in text
    assert task_index in text
    assert text.index(guard) < text.index(task_index)


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
