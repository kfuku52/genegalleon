# ruff: noqa: E501

import subprocess

from shell_static_helpers import (
    CORE_DIR,
    REPO_ROOT,
    WORKFLOW_DIR,
)
from shell_static_helpers import (
    function_body as _function_body,
)
from shell_static_helpers import (
    read_text as _read_text,
)
from shell_static_helpers import (
    workflow_shell_scripts as _workflow_shell_scripts,
)


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


def test_nonconda_download_helpers_use_archive_files_and_wget_fallback():
    script_path = REPO_ROOT / "container" / "scripts" / "install_nonconda_fallbacks.sh"
    text = _read_text(script_path)
    download_script_path = REPO_ROOT / "container" / "scripts" / "download_url.sh"
    download_script = _read_text(download_script_path)

    download_body = _function_body(text, "download_url_to_file")
    assert 'bash "${download_url_script}" "${url}" "${dest}"' in download_body

    assert "max_attempts=${DOWNLOAD_URL_MAX_ATTEMPTS:-8}" in download_script
    assert "retry_delay_sec=${DOWNLOAD_URL_RETRY_DELAY_SEC:-5}" in download_script
    assert "connect_timeout_sec=${DOWNLOAD_URL_CONNECT_TIMEOUT_SEC:-30}" in download_script
    assert "max_time_sec=${DOWNLOAD_URL_MAX_TIME_SEC:-600}" in download_script
    assert '--retry "${max_attempts}" \\' in download_script
    assert "--retry-all-errors \\" in download_script
    assert '--retry-delay "${retry_delay_sec}" \\' in download_script
    assert '--connect-timeout "${connect_timeout_sec}" \\' in download_script
    assert '--max-time "${max_time_sec}" \\' in download_script
    assert '--tries="${max_attempts}" \\' in download_script
    assert '--waitretry="${retry_delay_sec}" \\' in download_script
    assert '-O "${output_path}" \\' in download_script
    assert 'rm -f -- "${output_path}"' in download_script

    tag_body = _function_body(text, "download_github_tag_tarball")
    assert "archive_path=$(mktemp)" in tag_body
    assert 'if ! download_checked_url_to_file "${url}" "${archive_path}" "${expected_sha256}"; then' in tag_body
    assert 'if ! tar -xzf "${archive_path}" -C "${dest}" --strip-components=1; then' in tag_body

    cafe_body = _function_body(text, "install_cafe5")
    assert 'archive_path="${workdir}/CAFE5-5.1.0.tar.gz"' in cafe_body
    assert "if ! download_checked_url_to_file \\" in cafe_body
    assert 'install_r_cran_packages "${r_env_name}" Rphylopars Rtsne' in text
    assert "pkgs <- c('Rphylopars','Rtsne'" in text
    assert 'if ! tar -xzf "${archive_path}" -C "${workdir}"; then' in cafe_body


def test_download_lock_helper_uses_shared_lock_metadata_and_heartbeat():
    lock_path = WORKFLOW_DIR / "support" / "gg_shared_lock.sh"
    text = _read_text(lock_path)
    body = _function_body(text, "gg_array_download_once")
    assert 'if gg_artifact_ready "${artifact_path}"; then' in body
    assert 'if ! gg_shared_lock_acquire "${lock_file}" "${description}"; then' in body
    assert 'gg_shared_lock_start_heartbeat "${lock_file}"' in body
    assert "heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}" in body
    assert '"$@" >&2' in body
    assert 'gg_shared_lock_stop_heartbeat "${heartbeat_pid}"' in body
    assert 'gg_shared_lock_release "${lock_file}"' in body
    assert ".dlock" not in body


def test_download_lock_helper_redirects_artifact_stdout(tmp_path):
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    artifact_path = tmp_path / "artifact.ready"
    lock_path = tmp_path / "artifact.lock"
    script = f"""
set -euo pipefail
source "{util_path}"
GG_ARRAY_TASK_ID=1
GG_LOCK_HEARTBEAT_SECONDS=1
fake_builder() {{
  echo noisy-tool-stdout
  gg_write_ready_marker "{artifact_path}"
}}
gg_array_download_once "{lock_path}" "{artifact_path}" "fake artifact" fake_builder
"""

    result = subprocess.run(
        ["bash", "-c", script],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        timeout=10,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout == ""
    assert "noisy-tool-stdout" in result.stderr


def test_uniprot_mmseqs_helper_stdout_is_only_db_prefix(tmp_path):
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    workspace = tmp_path / "workspace"
    runtime_prefix = workspace / "downloads" / "uniprot_sprot" / "uniprot_sprot"
    runtime_prefix.parent.mkdir(parents=True)
    (runtime_prefix.with_suffix(".pep")).write_text(">sp|P1\nMA\n", encoding="utf-8")
    (runtime_prefix.with_suffix(".meta.tsv.gz")).write_text("metadata\n", encoding="utf-8")
    script = f"""
set -euo pipefail
source "{util_path}"
GG_ARRAY_TASK_ID=1
GG_LOCK_HEARTBEAT_SECONDS=1
mmseqs() {{
  if [[ "$1" == "createdb" ]]; then
    echo "createdb $2 $3"
    printf 'db\\n' > "$3"
    printf '0\\n' > "$3.dbtype"
    return 0
  fi
  return 1
}}
ensure_uniprot_sprot_mmseqs_db "{workspace}"
"""

    result = subprocess.run(
        ["bash", "-c", script],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        timeout=10,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout == str(runtime_prefix) + "\n"
    assert "createdb" in result.stderr


def test_uniprot_blast_helper_stdout_is_only_db_prefix(tmp_path):
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    workspace = tmp_path / "workspace"
    runtime_prefix = workspace / "downloads" / "uniprot_sprot" / "uniprot_sprot"
    runtime_prefix.parent.mkdir(parents=True)
    (runtime_prefix.with_suffix(".pep")).write_text(">sp|P1\nMA\n", encoding="utf-8")
    (runtime_prefix.with_suffix(".meta.tsv.gz")).write_text("metadata\n", encoding="utf-8")
    script = f"""
set -euo pipefail
source "{util_path}"
GG_ARRAY_TASK_ID=1
GG_LOCK_HEARTBEAT_SECONDS=1
makeblastdb() {{
  local out_prefix=""
  while [[ $# -gt 0 ]]; do
    if [[ "$1" == "-out" ]]; then
      out_prefix="$2"
      shift 2
    else
      shift
    fi
  done
  echo "makeblastdb $out_prefix"
  printf 'pin\\n' > "$out_prefix.pin"
  printf 'phr\\n' > "$out_prefix.phr"
  printf 'psq\\n' > "$out_prefix.psq"
}}
ensure_uniprot_sprot_blast_db "{workspace}"
"""

    result = subprocess.run(
        ["bash", "-c", script],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        timeout=10,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout == str(runtime_prefix) + "\n"
    assert "makeblastdb" in result.stderr


def test_uniprot_db_prefix_validator_rejects_contaminated_values(tmp_path):
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    script = f"""
set -euo pipefail
source "{util_path}"
validate_uniprot_sprot_db_prefix $'createdb /tmp/uniprot.pep\\n/tmp/uniprot.mmseqs' mmseqs2
"""

    result = subprocess.run(
        ["bash", "-c", script],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        timeout=10,
        check=False,
    )

    assert result.returncode != 0
    assert "DB prefix is malformed" in result.stderr


def test_uniprot_db_prefix_validator_allows_spaces_in_paths(tmp_path):
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    db_prefix = tmp_path / "uniprot sprot"
    (tmp_path / "uniprot sprot.mmseqs").write_text("db\n", encoding="utf-8")
    (tmp_path / "uniprot sprot.mmseqs.dbtype").write_text("0\n", encoding="utf-8")
    script = f"""
set -euo pipefail
source "{util_path}"
validate_uniprot_sprot_db_prefix "{db_prefix}" mmseqs2
"""

    result = subprocess.run(
        ["bash", "-c", script],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        timeout=10,
        check=False,
    )

    assert result.returncode == 0, result.stderr


def test_uniprot_db_prefixes_are_validated_after_capture():
    for script in [
        CORE_DIR / "gg_genome_annotation_core.sh",
        CORE_DIR / "gg_gene_evolution_core.sh",
        CORE_DIR / "gg_genome_evolution_core.sh",
    ]:
        text = _read_text(script)
        assert 'validate_uniprot_sprot_db_prefix "${uniprot_db_prefix}" "blastp"' in text
        assert 'validate_uniprot_sprot_db_prefix "${uniprot_db_prefix}" "mmseqs2"' in text


def test_download_entrypoints_use_shared_lock_helper_for_busco_and_ete():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    busco_body = _function_body(text, "ensure_busco_download_path")
    ete_body = _function_body(text, "ensure_ete_taxonomy_db")
    assert 'runtime_ready_marker="${runtime_busco_lineage}/.download.ready"' in busco_body
    assert 'gg_array_download_once "${lock_file}" "${runtime_ready_marker}"' in busco_body
    assert 'gg_array_download_once "${lock_file}" "${db_file}" "ETE taxonomy DB"' in ete_body


def test_ete_taxonomy_helper_uses_explicit_shared_taxdump_path():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    set_env_body = _function_body(text, "gg_set_taxonomy_cache_env")
    locked_body = _function_body(text, "_ensure_ete_taxonomy_db_locked")
    ensure_body = _function_body(text, "ensure_ete_taxonomy_db")
    assert "workspace_taxonomy_taxdumpfile()" in text
    assert 'ensure_dir "${dir_taxonomy}/ete4"' in set_env_body
    assert 'export GG_TAXONOMY_TAXDUMPFILE="${dir_taxonomy}/taxdump.tar.gz"' in set_env_body
    assert 'os.makedirs(os.path.join(cache_dir, "ete4"), exist_ok=True)' in locked_body
    assert 'GG_TAXONOMY_TAXDUMPFILE="${taxdump_file}"' in text
    assert "NCBITaxa(dbfile=db_file, taxdump_file=ensure_taxdump_file(), update=True)" in locked_body
    assert 'ensure_dir "${dir_taxonomy}/ete4"' in ensure_body
    assert 'taxdump_file=$(workspace_taxonomy_taxdumpfile "${gg_workspace_dir}")' in ensure_body
    assert '_ensure_ete_taxonomy_db_locked "${db_file}" "${taxdump_file}"' in ensure_body


def test_latest_jaspar_lock_uses_shared_lock_and_marker_resolution():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "ensure_latest_jaspar_file")
    assert (
        'if resolved_path=$(_resolve_latest_jaspar_path_from_marker "${latest_marker}" "${sys_dir}" "${runtime_dir}"); then'
        in body
    )
    assert 'if ! gg_shared_lock_acquire "${lock_file}" "latest JASPAR motif file"; then' in body
    assert 'gg_shared_lock_start_heartbeat "${lock_file}"' in body
    assert "heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}" in body
    assert 'gg_shared_lock_stop_heartbeat "${heartbeat_pid}"' in body
    assert 'gg_shared_lock_release "${lock_file}"' in body


def test_shared_lock_heartbeat_is_not_started_via_command_substitution():
    disallowed = "heartbeat_pid=$(gg_shared_lock_start_heartbeat"
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
    gate = "if [[ ${run_build_species_busco_summary} -ne 1 ]]; then"
    source_dir = "source_species_input_dir=$(effective_species_input_source_dir_path)"
    normalize = 'normalize_busco_table_naming "${dir_species_busco_full}" "${dir_species_busco_short}"'
    check = 'if ! is_species_set_identical "${source_species_input_dir}" "${dir_species_busco_full}"; then'
    collect = 'python "${gg_support_dir}/collect_common_BUSCO_genes.py" \\'
    assert "run_shared_busco_summary_stage() {" in text
    assert gate in text
    assert source_dir in text
    assert normalize in text
    assert check in text
    assert collect in text
    gate_idx = text.index(gate)
    source_dir_idx = text.index(source_dir, gate_idx)
    normalize_idx = text.index(normalize, source_dir_idx)
    check_idx = text.index(check, normalize_idx)
    collect_idx = text.index(collect, check_idx)
    assert gate_idx < source_dir_idx < normalize_idx < check_idx < collect_idx


def test_genome_busco_summary_syncs_from_shared_summary_and_gates_busco_getfasta():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    sync_fn = "sync_genome_busco_summary_table_from_shared() {"
    provenance = 'gg_artifact_prepare_stage sync_needs_update run_sync "${sync_provenance_args[@]}"'
    copy_stmt = 'cp_out "${file_species_busco_summary_table}" "${file_genome_busco_summary_table}"'
    sync_call = "sync_genome_busco_summary_table_from_shared || return $?"
    gate = 'disable_if_no_input_file "run_busco_dupaware_extract_fasta" "${file_genome_busco_summary_table}"'
    assert sync_fn in text
    assert provenance in text
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


def test_shared_species_busco_stage_matches_existing_outputs_by_species_name():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    body = _function_body(text, "run_shared_species_busco_stage")

    assert 'sp_ub=$(gg_species_name_from_path_or_dot "${seq_file}")' in body
    assert '-name "*busco.full.tsv"' in body
    assert '-name "*busco.short.txt"' in body
    assert 'busco_species=$(gg_species_name_from_path_or_dot "${busco_base}")' in body
    assert 'busco_output_exists_for_species "${dir_species_busco_full}" "${sp_ub}" "*busco.full.tsv"' in body
    assert 'busco_output_exists_for_species "${dir_species_busco_short}" "${sp_ub}" "*busco.short.txt"' in body


def test_busco_getfasta_step_is_gated_by_summary_table_presence():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    gate = 'disable_if_no_input_file "run_busco_dupaware_extract_fasta" "${file_genome_busco_summary_table}"'
    step = "if [[ ${busco_extract_needs_update} -eq 1 && ${run_busco_dupaware_extract_fasta} -eq 1 ]]; then"
    assert gate in text
    assert step in text
    assert text.index(gate) < text.index(step)


def test_busco_getfasta_step_defines_and_uses_its_duplicate_aware_helper():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    step = "if [[ ${busco_extract_needs_update} -eq 1 && ${run_busco_dupaware_extract_fasta} -eq 1 ]]; then"
    helper = "generate_genome_dupaware_busco_fasta() {"
    invoke = 'generate_genome_dupaware_busco_fasta "${busco_idx}" &'
    step_idx = text.index(step)
    helper_idx = text.index(helper, step_idx)
    invoke_idx = text.index(invoke, helper_idx)
    assert helper_idx < invoke_idx


def test_genome_evolution_core_uses_safe_busco_summary_count_helper():
    script = CORE_DIR / "gg_genome_evolution_core.sh"
    text = _read_text(script)
    assert 'num_busco_ids=$(get_busco_summary_gene_count "${file_species_busco_summary_table}")' in text
