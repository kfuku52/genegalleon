import gzip
import shlex
import subprocess
import time
from pathlib import Path

import pytest

GG_UTIL_PATH = Path(__file__).resolve().parents[1] / "support" / "gg_util.sh"
GG_ENTRYPOINT_CONFIG_VARS_PATH = Path(__file__).resolve().parents[1] / "support" / "gg_entrypoint_config_vars.sh"
GG_ENTRYPOINT_BOOTSTRAP_PATH = Path(__file__).resolve().parents[1] / "support" / "gg_entrypoint_bootstrap.sh"


def run_bash(cmd: str, cwd: Path):
    return subprocess.run(
        ["bash", "-lc", cmd],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )


def test_forward_config_vars_trims_registry_whitespace_for_genome_evolution_entrypoint(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "astral_min_tips=7; "
        "run_species_busco=1; "
        "forward_config_vars_to_container_env gg_genome_evolution_entrypoint.sh; "
        'printf "astral=%s\\nrun_species_busco=%s\\n" '
        '"${SINGULARITYENV_astral_min_tips:-}" "${SINGULARITYENV_run_species_busco:-}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        "astral=7",
        "run_species_busco=1",
    ]


def test_forward_config_vars_includes_species_tree_zip_options(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "species_tree_output_storage=files; "
        "species_tree_zip_compression=store; "
        "species_tree_zip_compression_level=0; "
        "forward_config_vars_to_container_env gg_genome_evolution_entrypoint.sh; "
        'printf "storage=%s\\ncompression=%s\\nlevel=%s\\n" '
        '"${SINGULARITYENV_species_tree_output_storage:-}" '
        '"${SINGULARITYENV_species_tree_zip_compression:-}" '
        '"${SINGULARITYENV_species_tree_zip_compression_level:-}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        "storage=files",
        "compression=store",
        "level=0",
    ]


def test_forward_config_vars_includes_gene_evolution_csubst_nonsyn_recode(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "csubst_nonsyn_recode=dayhoff6; "
        "forward_config_vars_to_container_env gg_gene_evolution_entrypoint.sh; "
        'printf "csubst_nonsyn_recode=%s\\n" '
        '"${SINGULARITYENV_csubst_nonsyn_recode:-}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "csubst_nonsyn_recode=dayhoff6"


def test_forward_config_vars_includes_binary_foreground_resolver_option(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "csubst_resolve_binary_foreground=yes; "
        "forward_config_vars_to_container_env gg_gene_evolution_entrypoint.sh; "
        'printf "resolver=%s\\n" '
        '"${SINGULARITYENV_csubst_resolve_binary_foreground:-}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "resolver=yes"


def test_forward_config_vars_includes_gene_evolution_csubst_scan_options(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "run_csubst_scan=1; "
        "csubst_scan_match=all; "
        "csubst_scan_n_permutations=25; "
        "forward_config_vars_to_container_env gg_gene_evolution_entrypoint.sh; "
        'printf "run=%s\\nmatch=%s\\nperm=%s\\n" '
        '"${SINGULARITYENV_run_csubst_scan:-}" '
        '"${SINGULARITYENV_csubst_scan_match:-}" '
        '"${SINGULARITYENV_csubst_scan_n_permutations:-}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        "run=1",
        "match=all",
        "perm=25",
    ]


def test_forward_config_vars_includes_gene_summary_csubst_site_nonsyn_recode(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "csubst_site_nonsyn_recode=dayhoff9; "
        "forward_config_vars_to_container_env gg_gene_summary_entrypoint.sh; "
        'printf "csubst_site_nonsyn_recode=%s\\n" '
        '"${SINGULARITYENV_csubst_site_nonsyn_recode:-}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "csubst_site_nonsyn_recode=dayhoff9"


def test_gene_summary_run_csubst_scan_aa_change_summary_env_override_and_forwarding(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"source {shlex.quote(str(GG_ENTRYPOINT_CONFIG_VARS_PATH))}; "
        "run_csubst_scan_aa_change_summary=0; "
        "GG_GENE_SUMMARY_RUN_CSUBST_SCAN_AA_CHANGE_SUMMARY=1; "
        "gg_apply_registered_env_overrides gg_gene_summary_entrypoint.sh; "
        "forward_config_vars_to_container_env gg_gene_summary_entrypoint.sh; "
        'printf "local=%s\\nforwarded=%s\\n" '
        '"${run_csubst_scan_aa_change_summary}" "${SINGULARITYENV_run_csubst_scan_aa_change_summary:-}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        "local=1",
        "forwarded=1",
    ]


def test_gene_summary_candidate_site_options_are_forwarded(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "run_csubst_scan_candidate_sites=1; "
        "csubst_scan_candidate_sites_min_support=5; "
        "csubst_scan_candidate_sites_q_threshold=0.01; "
        "forward_config_vars_to_container_env gg_gene_summary_entrypoint.sh; "
        'printf "run=%s\\nmin=%s\\nq=%s\\n" '
        '"${SINGULARITYENV_run_csubst_scan_candidate_sites:-}" '
        '"${SINGULARITYENV_csubst_scan_candidate_sites_min_support:-}" '
        '"${SINGULARITYENV_csubst_scan_candidate_sites_q_threshold:-}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        "run=1",
        "min=5",
        "q=0.01",
    ]


def test_gene_summary_candidate_site_scoped_env_override_is_applied(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"source {shlex.quote(str(GG_ENTRYPOINT_CONFIG_VARS_PATH))}; "
        "run_csubst_scan_candidate_sites=0; "
        "csubst_scan_candidate_sites_min_support=5; "
        "GG_GENE_SUMMARY_RUN_CSUBST_SCAN_CANDIDATE_SITES=1; "
        "GG_GENE_SUMMARY_CSUBST_SCAN_CANDIDATE_SITES_MIN_SUPPORT=7; "
        "gg_apply_registered_env_overrides gg_gene_summary_entrypoint.sh; "
        "forward_config_vars_to_container_env gg_gene_summary_entrypoint.sh; "
        'printf "run=%s\\nmin=%s\\nforwarded=%s\\n" '
        '"${run_csubst_scan_candidate_sites}" '
        '"${csubst_scan_candidate_sites_min_support}" '
        '"${SINGULARITYENV_run_csubst_scan_candidate_sites:-}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        "run=1",
        "min=7",
        "forwarded=1",
    ]


def test_export_var_to_container_env_ignores_invalid_variable_name(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "gg_export_var_to_container_env_if_set ' invalid variable'; "
        'printf "%s\\n" ok'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "ok"
    assert "ignoring invalid variable name" in completed.stderr


def test_scheduler_defaults_derives_tool_memory_budget(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "GG_TASK_CPUS=4; "
        "GG_MEM_PER_CPU_GB=8; "
        "ensure_gg_scheduler_defaults >/dev/null; "
        'printf "total=%s\\nreserve=%s\\ntool=%s\\n" '
        '"${GG_MEM_TOTAL_GB}" "${GG_MEM_TOOL_RESERVE_GB}" "${GG_MEM_TOOL_GB}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.splitlines() == [
        "total=32",
        "reserve=4",
        "tool=28",
    ]


def test_scheduler_defaults_clamps_manual_tool_memory_to_reserved_budget(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "GG_TASK_CPUS=4; "
        "GG_MEM_PER_CPU_GB=8; "
        "GG_MEM_TOTAL_GB=32; "
        "GG_MEM_TOOL_RESERVE_GB=4; "
        "GG_MEM_TOOL_GB=32; "
        "ensure_gg_scheduler_defaults >/dev/null; "
        'printf "tool=%s\\n" "${GG_MEM_TOOL_GB}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "tool=28"


def test_scheduler_defaults_uses_constant_like_reserve_for_large_jobs(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "GG_TASK_CPUS=8; "
        "GG_MEM_PER_CPU_GB=32; "
        "ensure_gg_scheduler_defaults >/dev/null; "
        'printf "total=%s\\nreserve=%s\\ntool=%s\\n" '
        '"${GG_MEM_TOTAL_GB}" "${GG_MEM_TOOL_RESERVE_GB}" "${GG_MEM_TOOL_GB}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.splitlines() == [
        "total=256",
        "reserve=4",
        "tool=252",
    ]


def test_apply_registered_env_overrides_uses_entrypoint_prefixes_and_empty_values(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"source {shlex.quote(str(GG_ENTRYPOINT_CONFIG_VARS_PATH))}; "
        "run_assembly=1; "
        "delete_tmp_dir=1; "
        "amalgkit_sra_strategy_query=default; "
        "GG_TRANSCRIPTOME_RUN_ASSEMBLY=0; "
        "GG_TRANSCRIPTOME_DELETE_TMP_DIR=0; "
        "GG_TRANSCRIPTOME_AMALGKIT_SRA_STRATEGY_QUERY=; "
        "gg_apply_registered_env_overrides gg_transcriptome_generation_entrypoint.sh delete_tmp_dir; "
        'printf "run_assembly=%s\\ndelete_tmp_dir=%s\\nstrategy=%s\\n" '
        '"${run_assembly}" "${delete_tmp_dir}" "${amalgkit_sra_strategy_query}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.splitlines() == [
        "run_assembly=0",
        "delete_tmp_dir=0",
        "strategy=",
    ]


def test_apply_registered_env_overrides_keeps_gg_input_prefix(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"source {shlex.quote(str(GG_ENTRYPOINT_CONFIG_VARS_PATH))}; "
        "provider=NCBI; "
        "trait_profile=none; "
        "GG_INPUT_PROVIDER=direct; "
        "GG_INPUT_TRAIT_PROFILE=gift_starter; "
        "gg_apply_registered_env_overrides gg_input_generation_entrypoint.sh; "
        'printf "provider=%s\\ntrait_profile=%s\\n" "${provider}" "${trait_profile}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.splitlines() == [
        "provider=direct",
        "trait_profile=gift_starter",
    ]


def test_prepare_entrypoint_runtime_snapshot_copies_core_script_to_job_task_dir(tmp_path):
    workflow_dir = tmp_path / "workflow"
    workspace_dir = tmp_path / "workspace"
    workflow_dir.mkdir()
    (workspace_dir / "output").mkdir(parents=True)
    core_script = workflow_dir / "gg_test_core.sh"
    entrypoint_script = workflow_dir / "gg_test_entrypoint.sh"
    core_script.write_text("echo core\\n", encoding="utf-8")
    entrypoint_script.write_text("echo entrypoint\\n", encoding="utf-8")
    (tmp_path / "VERSION").write_text("0.0.0\\n", encoding="utf-8")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"gg_workflow_dir={shlex.quote(str(workflow_dir))}; "
        f"gg_workspace_dir={shlex.quote(str(workspace_dir))}; "
        f"gg_workspace_output_dir={shlex.quote(str(workspace_dir / 'output'))}; "
        "GG_JOB_ID='job:42'; "
        "GG_ARRAY_TASK_ID='task/7'; "
        f"snapshot=$(gg_prepare_entrypoint_runtime_snapshot gg_test_entrypoint.sh {shlex.quote(str(core_script))}); "
        'snapshot_dir="$(dirname "${snapshot}")"; '
        'printf "snapshot=%s\\n" "${snapshot}"; '
        'printf "core=%s\\n" "$([[ -s "${snapshot}" ]] && echo 1 || echo 0)"; '
        'printf "entrypoint=%s\\n" "$([[ -s "${snapshot_dir}/gg_test_entrypoint.sh" ]] && echo 1 || echo 0)"; '
        'printf "version=%s\\n" "$([[ -s "${snapshot_dir}/VERSION" ]] && echo 1 || echo 0)"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    lines = completed.stdout.splitlines()
    assert lines[0].endswith("/runtime/gg_test_entrypoint/job_42_task_7/gg_test_core.sh")
    assert lines[1:] == [
        "core=1",
        "entrypoint=1",
        "version=1",
    ]


def test_workspace_pfam_le_dir_is_under_downloads_dedicated_folder(tmp_path):
    project_dir = tmp_path / "project"
    project_dir.mkdir()
    command = f"source {shlex.quote(str(GG_UTIL_PATH))}; workspace_pfam_le_dir {shlex.quote(str(project_dir))}"

    completed = run_bash(command, cwd=tmp_path)
    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(project_dir / "downloads" / "pfam" / "Pfam_LE")


def test_ensure_uniprot_sprot_metadata_tsv_builds_runtime_meta_from_runtime_dat(tmp_path):
    workspace_dir = tmp_path / "workspace"
    runtime_prefix = workspace_dir / "downloads" / "uniprot_sprot" / "uniprot_sprot"
    runtime_dat = runtime_prefix.with_suffix(".dat.gz")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"workspace_dir={shlex.quote(str(workspace_dir))}; "
        f"runtime_prefix={shlex.quote(str(runtime_prefix))}; "
        f"runtime_dat={shlex.quote(str(runtime_dat))}; "
        'mkdir -p "$(dirname "${runtime_prefix}")" "$(dirname "${workspace_dir}/downloads/locks/.x")"; '
        "printf 'ID   TEST1_HUMAN              Reviewed;         100 AA.\\nAC   P12345;\\nOS   Homo sapiens (Human).\\nOX   NCBI_TaxID=9606;\\n//\\n' "
        '| gzip -c > "${runtime_dat}"; '
        'meta_path=$(ensure_uniprot_sprot_metadata_tsv "${workspace_dir}" "${runtime_prefix}"); '
        'printf "meta=%s\\n" "${meta_path}"; '
        'test -s "${meta_path}"; printf "size=%s\\n" "$?"; '
        "python -c \"import gzip,sys; t=gzip.open(sys.argv[1],'rt',encoding='utf-8').read(); "
        "print('has_taxid=' + ('1' if 'taxid' in t else '0')); "
        "print('has_accession=' + ('1' if 'P12345' in t else '0'))\" \"${meta_path}\""
    )

    completed = run_bash(command, cwd=tmp_path)
    assert completed.returncode == 0, completed.stderr
    lines = completed.stdout.strip().splitlines()
    assert lines[0] == f"meta={runtime_prefix}.meta.tsv.gz"
    assert lines[1] == "size=0"
    assert lines[2] == "has_taxid=1"
    assert lines[3] == "has_accession=1"


def test_ensure_uniprot_sprot_metadata_tsv_falls_back_to_runtime_path_for_sys_prefix(tmp_path):
    workspace_dir = tmp_path / "workspace"
    runtime_prefix = workspace_dir / "downloads" / "uniprot_sprot" / "uniprot_sprot"
    runtime_dat = runtime_prefix.with_suffix(".dat.gz")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"workspace_dir={shlex.quote(str(workspace_dir))}; "
        f"runtime_prefix={shlex.quote(str(runtime_prefix))}; "
        f"runtime_dat={shlex.quote(str(runtime_dat))}; "
        'mkdir -p "$(dirname "${runtime_prefix}")" "$(dirname "${workspace_dir}/downloads/locks/.x")"; '
        "printf 'ID   TEST2_HUMAN              Reviewed;         100 AA.\\nAC   Q99999;\\nOS   Homo sapiens (Human).\\nOX   NCBI_TaxID=9606;\\n//\\n' "
        '| gzip -c > "${runtime_dat}"; '
        'meta_path=$(ensure_uniprot_sprot_metadata_tsv "${workspace_dir}" "/usr/local/db/uniprot_sprot"); '
        'printf "meta=%s\\n" "${meta_path}"; '
        'test -s "${meta_path}"; printf "size=%s\\n" "$?"'
    )

    completed = run_bash(command, cwd=tmp_path)
    assert completed.returncode == 0, completed.stderr
    lines = completed.stdout.strip().splitlines()
    assert lines[0] == f"meta={runtime_prefix}.meta.tsv.gz"
    assert lines[1] == "size=0"


def test_gg_array_download_once_accepts_nonempty_ready_marker(tmp_path):
    lock_file = tmp_path / "locks" / "artifact.lock"
    marker_file = tmp_path / "runtime" / ".artifact.ready"
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"marker_file={shlex.quote(str(marker_file))}; "
        "SECONDS=0; "
        f'gg_array_download_once {shlex.quote(str(lock_file))} "$marker_file" "marker artifact" '
        'gg_write_ready_marker "$marker_file"; '
        "status=$?; "
        "heartbeat_elapsed=$SECONDS; "
        'printf "%s\\n%s\\n" "$status" "$heartbeat_elapsed"; '
        'wc -c < "$marker_file"; '
        'cat "$marker_file"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    lines = completed.stdout.strip().splitlines()
    assert lines[0] == "0"
    assert int(lines[1]) < 5, f"heartbeat shutdown took {lines[1]}s"
    assert int(lines[2]) > 0
    assert lines[3] == "ready"


def test_download_busco_lineage_to_runtime_merges_into_existing_runtime_db(tmp_path):
    runtime_db = tmp_path / "workspace" / "downloads" / "busco_downloads"
    existing_lineage = runtime_db / "lineages" / "existing_odb12"
    existing_lineage.mkdir(parents=True)
    (existing_lineage / "dataset.cfg").write_text("existing\n", encoding="utf-8")
    placement_dir = runtime_db / "placement_files"
    placement_dir.mkdir(parents=True)
    (placement_dir / "mapping.txt").write_text("mapping\n", encoding="utf-8")
    ready_marker = runtime_db / "lineages" / "eukaryota_odb12" / ".download.ready"
    runtime_lineage = runtime_db / "lineages" / "eukaryota_odb12"

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "busco() { "
        "mkdir -p busco_downloads/lineages/eukaryota_odb12/info; "
        "mkdir -p busco_downloads/placement_files; "
        "printf 'new\\n' > busco_downloads/lineages/eukaryota_odb12/dataset.cfg; "
        "printf 'placement\\n' > busco_downloads/placement_files/new_mapping.txt; "
        "printf 'versions\\n' > busco_downloads/file_versions.tsv; "
        "}; "
        f'_download_busco_lineage_to_runtime "eukaryota_odb12" '
        f"{shlex.quote(str(runtime_db))} "
        f"{shlex.quote(str(runtime_lineage))} "
        f"{shlex.quote(str(ready_marker))}; "
        "status=$?; "
        'printf "%s\\n" "$status"; '
        f'test -s {shlex.quote(str(runtime_lineage / "dataset.cfg"))}; printf "new=%s\\n" "$?"; '
        f'test -s {shlex.quote(str(existing_lineage / "dataset.cfg"))}; printf "existing=%s\\n" "$?"; '
        f'test -s {shlex.quote(str(placement_dir / "mapping.txt"))}; printf "mapping=%s\\n" "$?"; '
        f'test -s {shlex.quote(str(placement_dir / "new_mapping.txt"))}; printf "new_mapping=%s\\n" "$?"; '
        f'test -s {shlex.quote(str(ready_marker))}; printf "ready=%s\\n" "$?"; '
        'test ! -e busco_downloads; printf "staging_removed=%s\\n" "$?"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        "0",
        "new=0",
        "existing=0",
        "mapping=0",
        "new_mapping=0",
        "ready=0",
        "staging_removed=0",
    ]


def test_contamination_rank_normalizes_superkingdom_for_amalgkit(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; gg_normalize_contamination_removal_rank_for_amalgkit superkingdom"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "domain"


def test_contamination_rank_normalizes_domain_for_remove_contaminated_sequences(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "gg_normalize_contamination_removal_rank_for_remove_contaminated_sequences domain"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "superkingdom"


def test_entrypoint_bootstrap_resolves_workflow_dir_from_sourced_bootstrap_when_entrypoint_is_spooled(tmp_path):
    unrelated_submit_dir = tmp_path / "submit"
    unrelated_submit_dir.mkdir()
    expected_workflow_dir = GG_ENTRYPOINT_BOOTSTRAP_PATH.parents[1]

    command = (
        f"source {shlex.quote(str(GG_ENTRYPOINT_BOOTSTRAP_PATH))}; "
        f"export SLURM_SUBMIT_DIR={shlex.quote(str(unrelated_submit_dir))}; "
        f"export PBS_O_WORKDIR={shlex.quote(str(unrelated_submit_dir))}; "
        f"cd {shlex.quote(str(unrelated_submit_dir))}; "
        "gg_set_workflow_dir /var/spool/slurm/d/job15582515/slurm_script; "
        'printf "%s\\n" "${gg_workflow_dir}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(expected_workflow_dir)


def test_all_entrypoint_locators_try_later_project_directory_for_spooled_scripts(tmp_path):
    workflow_dir = GG_ENTRYPOINT_BOOTSTRAP_PATH.parents[1]
    repo_root = workflow_dir.parent
    unrelated_submit_dir = tmp_path / "submit"
    unrelated_submit_dir.mkdir()

    for entrypoint in sorted(workflow_dir.glob("gg_*_entrypoint.sh")):
        text = entrypoint.read_text(encoding="utf-8")
        start = text.index("# Resolve workflow paths for local and scheduler-spooled execution.")
        end = text.index("gg_entrypoint_name=", start)
        spooled_script = tmp_path / f"{entrypoint.stem}-slurm_script"
        spooled_script.write_text(
            "#!/usr/bin/env bash\n"
            "set -euo pipefail\n" + text[start:end] + 'printf "resolved=%s\\n" "${gg_workflow_dir}"\n',
            encoding="utf-8",
        )
        command = (
            f"export SLURM_SUBMIT_DIR={shlex.quote(str(unrelated_submit_dir))}; "
            f"export PBS_O_WORKDIR={shlex.quote(str(repo_root))}; "
            f"cd {shlex.quote(str(repo_root))}; "
            f"bash {shlex.quote(str(spooled_script))}"
        )

        completed = run_bash(command, cwd=tmp_path)

        assert completed.returncode == 0, f"{entrypoint.name}: {completed.stderr}"
        assert completed.stdout.strip() == f"resolved={workflow_dir}"


def test_add_container_bind_mount_uses_only_singularity_bindpath_for_singularity_runtime(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "singularity_command='singularity exec'; "
        "gg_add_container_bind_mount '/host/workspace:/workspace'; "
        "gg_add_container_bind_mount '/host/workflow:/script'; "
        'printf "GG=%s\\nSB=%s\\nSBP=%s\\nAB=%s\\nABP=%s\\n" '
        '"${GG_CONTAINER_BIND_MOUNTS:-}" "${SINGULARITY_BIND:-}" "${SINGULARITY_BINDPATH:-}" '
        '"${APPTAINER_BIND:-}" "${APPTAINER_BINDPATH:-}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        "GG=/host/workflow:/script,/host/workspace:/workspace",
        "SB=",
        "SBP=/host/workflow:/script,/host/workspace:/workspace",
        "AB=",
        "ABP=",
    ]


def test_add_container_bind_mount_uses_only_apptainer_bindpath_for_apptainer_runtime(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "singularity_command='apptainer exec'; "
        "gg_add_container_bind_mount '/host/workspace:/workspace'; "
        "gg_add_container_bind_mount '/host/workflow:/script'; "
        'printf "GG=%s\\nSB=%s\\nSBP=%s\\nAB=%s\\nABP=%s\\n" '
        '"${GG_CONTAINER_BIND_MOUNTS:-}" "${SINGULARITY_BIND:-}" "${SINGULARITY_BINDPATH:-}" '
        '"${APPTAINER_BIND:-}" "${APPTAINER_BINDPATH:-}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        "GG=/host/workflow:/script,/host/workspace:/workspace",
        "SB=",
        "SBP=",
        "AB=",
        "ABP=/host/workflow:/script,/host/workspace:/workspace",
    ]


def test_resolve_annotation_species_prefers_known_model_species(tmp_path):
    species_dir = tmp_path / "species_cds"
    species_dir.mkdir()
    (species_dir / "Cephalotus_follicularis.fa").write_text(">a\nATG\n")
    (species_dir / "Arabidopsis_thaliana.fa").write_text(">a\nATG\n")
    (species_dir / "Oryza_sativa.fa").write_text(">a\nATG\n")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"set -- $(gg_species_names_from_fasta_dir {shlex.quote(str(species_dir))}); "
        'gg_resolve_annotation_species auto "$@"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "Arabidopsis_thaliana"


def test_resolve_annotation_species_uses_first_available_when_no_model_species_exists(tmp_path):
    species_dir = tmp_path / "species_cds"
    species_dir.mkdir()
    (species_dir / "Cephalotus_follicularis.fa").write_text(">a\nATG\n")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"set -- $(gg_species_names_from_fasta_dir {shlex.quote(str(species_dir))}); "
        'gg_resolve_annotation_species auto "$@"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "Cephalotus_follicularis"


def test_resolve_annotation_species_prefers_cross_clade_model_species(tmp_path):
    species_dir = tmp_path / "species_cds"
    species_dir.mkdir()
    (species_dir / "Cephalotus_follicularis.fa").write_text(">a\nATG\n")
    (species_dir / "Danio_rerio.fa").write_text(">a\nATG\n")
    (species_dir / "Escherichia_coli.fa").write_text(">a\nATG\n")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"set -- $(gg_species_names_from_fasta_dir {shlex.quote(str(species_dir))}); "
        'gg_resolve_annotation_species auto "$@"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "Danio_rerio"


def test_resolve_annotation_species_normalizes_legacy_trailing_underscore(tmp_path):
    command = f'source {shlex.quote(str(GG_UTIL_PATH))}; gg_resolve_annotation_species "Arabidopsis_thaliana_"'

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "Arabidopsis_thaliana"


def test_species_name_from_path_or_dot_preserves_taxonomic_qualifiers(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        'printf "%s\\n%s\\n%s\\n%s\\n%s\\n%s\\n%s\\n%s\\n%s\\n%s\\n%s\\n" '
        '"$(gg_species_name_from_path_or_dot "Dictyostelium_cf_discoideum_GCA_054859205.1.fa.gz")" '
        '"$(gg_species_name_from_path_or_dot "Bacillus_subtilis_subsp._subtilis_demo.fa.gz")" '
        '"$(gg_species_name_from_path_or_dot "Amoeba_sp._JDS-Ruffled.tsv")" '
        '"$(gg_species_name_from_path_or_dot "Arisaema_sp._aooni_longestCDS.fa.gz")" '
        '"$(gg_species_name_from_path_or_dot "Asimitellaria_furusei_var._subramosa_busco.full.tsv")" '
        '"$(gg_species_name_from_path_or_dot "Asimitellaria_furusei_var._furusei.cds.fa.gz")" '
        '"$(gg_species_name_from_path_or_dot "homo_sapiens.GRCh38.cds.all.fa.gz")" '
        '"$(gg_species_name_from_path_or_dot "Amphizonella_sp_longestCDS.fa.gz")" '
        '"$(gg_species_name_from_path_or_dot "Cunea_sp_longestCDS_contamination_removal.fa.gz")" '
        '"$(gg_species_name_from_path_or_dot "Vannella_sp_longestCDS.transcript.fa.gz")" '
        '"$(gg_species_name_from_path_or_dot "Vexillifera_sp_longestCDS.busco.full.tsv")"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        "Dictyostelium_cf_discoideum",
        "Bacillus_subtilis_subsp._subtilis",
        "Amoeba_sp._JDS-Ruffled",
        "Arisaema_sp._aooni",
        "Asimitellaria_furusei_var._subramosa",
        "Asimitellaria_furusei_var._furusei",
        "homo_sapiens",
        "Amphizonella_sp",
        "Cunea_sp",
        "Vannella_sp",
        "Vexillifera_sp",
    ]


def test_check_species_sequences_accept_transcriptome_longest_cds_for_genus_sp(tmp_path):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    seqkit = bin_dir / "seqkit"
    seqkit.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "[[ ${1:-} == seq ]]\n"
        "shift\n"
        "[[ ${1:-} == --name ]] && shift\n"
        "if [[ ${1:-} == --threads ]]; then shift 2; fi\n"
        'gzip -cd -- ${1:?} | awk \'/^>/ {sub(/^>/, ""); sub(/[[:space:]].*$/, ""); print}\'\n',
        encoding="utf-8",
    )
    seqkit.chmod(0o755)

    cases = (
        ("species_cds", "check_species_cds_dir", "ATG", "All per-species CDS files are valid."),
        (
            "species_protein",
            "check_species_protein_dir",
            "MPEP",
            "All per-species protein files are valid.",
        ),
    )
    for directory_name, validator, sequence, success_message in cases:
        species_dir = tmp_path / directory_name
        species_dir.mkdir()
        fasta_path = species_dir / "Amphizonella_sp_longestCDS.fa.gz"
        with gzip.open(fasta_path, "wt", encoding="utf-8") as handle:
            handle.write(f">Amphizonella_sp_g0\n{sequence}\n>Amphizonella_sp_g1\n{sequence}\n")

        command = (
            f"export PATH={shlex.quote(str(bin_dir))}:$PATH; "
            f"source {shlex.quote(str(GG_UTIL_PATH))}; "
            "GG_TASK_CPUS=1; "
            f"{validator} {shlex.quote(str(species_dir))}"
        )
        completed = run_bash(command, cwd=tmp_path)

        assert completed.returncode == 0, completed.stderr + completed.stdout
        assert success_message in completed.stdout


def test_fasta_relabel_headers_to_species_preserves_taxonomic_qualifiers(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "printf '>Dictyostelium_cf_discoideum_gene1\\nATG\\n>Bacillus_subtilis_subsp_subtilis_gene2\\nATG\\n>Asimitellaria_furusei_var._subramosa_gene3\\nATG\\n>Arisaema_sp._aooni_gene4\\nATG\\n' "
        "| gg_fasta_relabel_headers_to_species"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        ">Dictyostelium_cf_discoideum",
        "ATG",
        ">Bacillus_subtilis_subsp_subtilis",
        "ATG",
        ">Asimitellaria_furusei_var._subramosa",
        "ATG",
        ">Arisaema_sp._aooni",
        "ATG",
    ]


def test_resolve_busco_lineage_from_lineages_prefers_deepest_mapped_taxon(tmp_path):
    mapping_dir = tmp_path / "busco_mappings"
    mapping_dir.mkdir()
    (mapping_dir / "mapping_taxids-busco_dataset_name.eukaryota_odb12.test.txt").write_text(
        "2759\teukaryota_odb12\n33090\tviridiplantae_odb12\n3193\tembryophyta_odb12\n3744\trosales_odb12\n"
    )

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"gg_resolve_busco_lineage_from_lineages auto {shlex.quote(str(mapping_dir))} "
        '"1,131567,2759,33090,3193,3744"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "rosales_odb12"


def test_resolve_busco_lineage_from_lineages_uses_deepest_common_taxon(tmp_path):
    mapping_dir = tmp_path / "busco_mappings"
    mapping_dir.mkdir()
    (mapping_dir / "mapping_taxids-busco_dataset_name.eukaryota_odb12.test.txt").write_text(
        "2759\teukaryota_odb12\n"
        "33090\tviridiplantae_odb12\n"
        "3193\tembryophyta_odb12\n"
        "3700\tbrassicales_odb12\n"
        "4530\tpoales_odb12\n"
    )

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"gg_resolve_busco_lineage_from_lineages auto {shlex.quote(str(mapping_dir))} "
        '"1,131567,2759,33090,3193,3700" '
        '"1,131567,2759,33090,3193,4530"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "embryophyta_odb12"


def test_resolve_busco_lineage_from_lineages_passes_through_explicit_value(tmp_path):
    mapping_dir = tmp_path / "busco_mappings"
    mapping_dir.mkdir()

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"gg_resolve_busco_lineage_from_lineages metazoa_odb13 {shlex.quote(str(mapping_dir))} "
        '"1,131567,2759,33208"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "metazoa_odb13"


def test_resolve_busco_lineage_from_lineages_prefers_latest_common_odb_version(tmp_path):
    mapping_dir = tmp_path / "busco_mappings"
    mapping_dir.mkdir()
    (mapping_dir / "mapping_taxids-busco_dataset_name.archaea_odb13.test.txt").write_text("2157\tarchaea_odb13\n")
    (mapping_dir / "mapping_taxids-busco_dataset_name.bacteria_odb13.test.txt").write_text("2\tbacteria_odb13\n")
    (mapping_dir / "mapping_taxids-busco_dataset_name.eukaryota_odb12.test.txt").write_text(
        "2759\teukaryota_odb12\n33090\tviridiplantae_odb12\n3193\tembryophyta_odb12\n"
    )
    (mapping_dir / "mapping_taxids-busco_dataset_name.eukaryota_odb13.test.txt").write_text(
        "2759\teukaryota_odb13\n33090\tviridiplantae_odb13\n3193\tembryophyta_odb13\n3744\trosales_odb13\n"
    )

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"gg_resolve_busco_lineage_from_lineages auto {shlex.quote(str(mapping_dir))} "
        '"1,131567,2759,33090,3193,3744"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "rosales_odb13"


def test_finalize_auto_busco_lineage_name_appends_requested_odb_suffix(tmp_path):
    command = f"source {shlex.quote(str(GG_UTIL_PATH))}; gg_finalize_auto_busco_lineage_name brassicales 13"

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "brassicales_odb13"


def test_finalize_auto_busco_lineage_name_preserves_existing_suffix(tmp_path):
    command = f"source {shlex.quote(str(GG_UTIL_PATH))}; gg_finalize_auto_busco_lineage_name embryophyta_odb12"

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "embryophyta_odb12"


def test_cdskit_localize_organism_group_infers_from_busco_lineage(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "printf '%s\\n%s\\n%s\\n' "
        '"$(gg_cdskit_localize_organism_group_from_busco_lineage embryophyta_odb12)" '
        '"$(gg_cdskit_localize_organism_group_from_busco_lineage metazoa_odb12)" '
        '"$(gg_cdskit_localize_organism_group_from_busco_lineage eukaryota_odb12)"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == ["plant", "non_plant", "unknown"]


def test_cdskit_localize_organism_group_resolver_uses_explicit_group_or_lineage(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "printf '%s\\n%s\\n' "
        '"$(gg_resolve_cdskit_localize_organism_group nonplant '
        f'{shlex.quote(str(tmp_path))} auto)" '
        '"$(gg_resolve_cdskit_localize_organism_group auto '
        f'{shlex.quote(str(tmp_path))} viridiplantae_odb12)"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == ["non_plant", "plant"]


def test_workspace_layout_defaults_to_split_for_empty_workspace(tmp_path):
    project_dir = tmp_path / "project"
    project_dir.mkdir()
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f'gg_prepare_cmd_runtime {shlex.quote(str(project_dir))} "" 0 0; '
        "printf '%s\\n%s\\n%s\\n%s\\n' "
        '"${gg_workspace_layout_resolved}" '
        '"${gg_workspace_input_dir}" '
        '"${gg_workspace_output_dir}" '
        '"${gg_workspace_downloads_dir}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    lines = completed.stdout.strip().splitlines()
    assert lines == [
        "split",
        str(project_dir / "input"),
        str(project_dir / "output"),
        str(project_dir / "downloads"),
    ]
    assert (project_dir / "input").is_dir()
    assert (project_dir / "output").is_dir()
    assert (project_dir / "downloads").is_dir()


def test_workspace_layout_ignores_root_level_entries_and_stays_split(tmp_path):
    project_dir = tmp_path / "project"
    (project_dir / "species_cds").mkdir(parents=True)
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "printf '%s\\n%s\\n%s\\n' "
        f'"$(gg_resolve_workspace_layout {shlex.quote(str(project_dir))})" '
        f'"$(workspace_input_root {shlex.quote(str(project_dir))})" '
        f'"$(workspace_output_root {shlex.quote(str(project_dir))})"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        "split",
        str(project_dir / "input"),
        str(project_dir / "output"),
    ]


def test_workspace_layout_no_longer_honors_legacy_override(tmp_path):
    project_dir = tmp_path / "project"
    project_dir.mkdir()
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f'gg_prepare_cmd_runtime {shlex.quote(str(project_dir))} "" 0 0; '
        "printf '%s\\n%s\\n%s\\n%s\\n' "
        '"${gg_workspace_layout_resolved}" '
        '"${gg_workspace_input_dir}" '
        '"${gg_workspace_output_dir}" '
        '"${gg_workspace_downloads_dir}"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        "split",
        str(project_dir / "input"),
        str(project_dir / "output"),
        str(project_dir / "downloads"),
    ]
    assert (project_dir / "input").is_dir()
    assert (project_dir / "output").is_dir()
    assert (project_dir / "downloads").is_dir()


def test_ensure_pfam_le_db_uses_new_workspace_layout_without_migrating_legacy_dir(tmp_path):
    project_dir = tmp_path / "project"
    legacy_dir = project_dir / "downloads" / "Pfam_LE"
    new_dir = project_dir / "downloads" / "pfam" / "Pfam_LE"
    legacy_dir.mkdir(parents=True)
    new_dir.mkdir(parents=True)
    (legacy_dir / "Pfam.pal").write_text("dummy\n")
    (new_dir / "Pfam.pal").write_text("dummy\n")

    command = f"source {shlex.quote(str(GG_UTIL_PATH))}; ensure_pfam_le_db {shlex.quote(str(project_dir))}"
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(new_dir)
    assert (new_dir / "Pfam.pal").is_file()
    assert legacy_dir.is_dir()


def test_ensure_pfam_le_db_backfills_nonempty_ready_marker(tmp_path):
    project_dir = tmp_path / "project"
    new_dir = project_dir / "downloads" / "pfam" / "Pfam_LE"
    ready_file = new_dir / ".pfam_le.ready"
    new_dir.mkdir(parents=True)
    (new_dir / "Pfam.pal").write_text("dummy\n")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"runtime_dir=$(ensure_pfam_le_db {shlex.quote(str(project_dir))}); "
        "status=$?; "
        'printf "%s\\n%s\\n" "$status" "$runtime_dir"; '
        f"wc -c < {shlex.quote(str(ready_file))}; "
        f"cat {shlex.quote(str(ready_file))}"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    lines = completed.stdout.strip().splitlines()
    assert lines[0] == "0"
    assert lines[1] == str(new_dir)
    assert int(lines[2]) > 0
    assert lines[3] == "ready"


def test_mv_out_accepts_pipe_input(tmp_path):
    outfile = tmp_path / "nested" / "out.txt"
    command = f"source {shlex.quote(str(GG_UTIL_PATH))}; printf 'hello\\n' | mv_out {shlex.quote(str(outfile))}"

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert outfile.read_text() == "hello\n"


def test_cp_out_accepts_pipe_input(tmp_path):
    outfile = tmp_path / "nested" / "out.txt"
    command = f"source {shlex.quote(str(GG_UTIL_PATH))}; printf 'world\\n' | cp_out {shlex.quote(str(outfile))}"

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert outfile.read_text() == "world\n"


def test_cp_out_single_argument_without_pipe_fails(tmp_path):
    outfile = tmp_path / "out.txt"
    command = f"source {shlex.quote(str(GG_UTIL_PATH))}; cp_out {shlex.quote(str(outfile))}"

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode != 0
    assert "at least 2 arguments are required unless stdin is piped" in completed.stdout


def test_cp_out_creates_destination_dir_when_target_has_trailing_slash(tmp_path):
    src = tmp_path / "src.txt"
    src.write_text("abc\n")
    dest_dir = tmp_path / "nested" / "dest_dir"
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; cp_out {shlex.quote(str(src))} {shlex.quote(str(dest_dir) + '/')}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    copied = dest_dir / src.name
    assert copied.read_text() == "abc\n"


def test_cp_out_preserves_existing_destination_when_copy_fails(tmp_path):
    src = tmp_path / "src.txt"
    destination = tmp_path / "out" / "result.txt"
    src.write_text("new\n")
    destination.parent.mkdir()
    destination.write_text("old\n")
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        'cp() { local destination="${!#}"; printf "partial\\n" > "${destination}"; return 9; }; '
        f"cp_out {shlex.quote(str(src))} {shlex.quote(str(destination))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode != 0
    assert destination.read_text() == "old\n"
    assert not list(destination.parent.glob("*.gg-copy.*"))


def test_mv_out_creates_destination_dir_for_multi_source_move(tmp_path):
    src1 = tmp_path / "a.txt"
    src2 = tmp_path / "b.txt"
    src1.write_text("a\n")
    src2.write_text("b\n")
    dest_dir = tmp_path / "nested" / "mv_dest"
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"mv_out {shlex.quote(str(src1))} {shlex.quote(str(src2))} {shlex.quote(str(dest_dir))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert not src1.exists()
    assert not src2.exists()
    assert (dest_dir / "a.txt").read_text() == "a\n"
    assert (dest_dir / "b.txt").read_text() == "b\n"


def test_mv_out_bundle_publishes_all_pairs_after_complete_staging(tmp_path):
    first_source = tmp_path / "stage" / "first.txt"
    second_source = tmp_path / "stage" / "second.txt"
    first_destination = tmp_path / "out-a" / "first.txt"
    second_destination = tmp_path / "out-b" / "second.txt"
    first_source.parent.mkdir()
    first_destination.parent.mkdir()
    second_destination.parent.mkdir()
    first_source.write_text("new-first\n")
    second_source.write_text("new-second\n")
    first_destination.write_text("old-first\n")
    second_destination.write_text("old-second\n")
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "mv_out_bundle "
        f"{shlex.quote(str(first_source))} {shlex.quote(str(first_destination))} "
        f"{shlex.quote(str(second_source))} {shlex.quote(str(second_destination))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert first_destination.read_text() == "new-first\n"
    assert second_destination.read_text() == "new-second\n"
    assert not first_source.exists()
    assert not second_source.exists()
    assert not list(tmp_path.rglob("*.gg-stage.*"))
    assert not list(tmp_path.rglob("*.gg-backup.*"))


def test_mv_out_bundle_rolls_back_every_pair_after_publish_failure(tmp_path):
    first_source = tmp_path / "stage" / "first.txt"
    second_source = tmp_path / "stage" / "second.txt"
    first_destination = tmp_path / "out-a" / "first.txt"
    second_destination = tmp_path / "out-b" / "second.txt"
    first_source.parent.mkdir()
    first_destination.parent.mkdir()
    second_destination.parent.mkdir()
    first_source.write_text("new-first\n")
    second_source.write_text("new-second\n")
    first_destination.write_text("old-first\n")
    second_destination.write_text("old-second\n")
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "GG_TEST_MV_COUNT=0; "
        "mv() { GG_TEST_MV_COUNT=$((GG_TEST_MV_COUNT + 1)); "
        'if [[ ${GG_TEST_MV_COUNT} -eq 6 ]]; then return 1; fi; command mv "$@"; }; '
        "mv_out_bundle "
        f"{shlex.quote(str(first_source))} {shlex.quote(str(first_destination))} "
        f"{shlex.quote(str(second_source))} {shlex.quote(str(second_destination))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode != 0
    assert first_source.read_text() == "new-first\n"
    assert second_source.read_text() == "new-second\n"
    assert first_destination.read_text() == "old-first\n"
    assert second_destination.read_text() == "old-second\n"
    assert not list(tmp_path.rglob("*.gg-stage.*"))
    assert not list(tmp_path.rglob("*.gg-backup.*"))


@pytest.mark.parametrize("failed_move", range(1, 7))
def test_mv_out_bundle_recovers_when_mv_changes_state_then_fails(tmp_path, failed_move):
    first_source = tmp_path / "stage" / "first.txt"
    second_source = tmp_path / "stage" / "second.txt"
    first_destination = tmp_path / "out-a" / "first.txt"
    second_destination = tmp_path / "out-b" / "second.txt"
    first_source.parent.mkdir()
    first_destination.parent.mkdir()
    second_destination.parent.mkdir()
    first_source.write_text("new-first\n")
    second_source.write_text("new-second\n")
    first_destination.write_text("old-first\n")
    second_destination.write_text("old-second\n")
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "GG_TEST_MV_COUNT=0; "
        "mv() { GG_TEST_MV_COUNT=$((GG_TEST_MV_COUNT + 1)); "
        'command mv "$@" || return $?; '
        f"if [[ ${{GG_TEST_MV_COUNT}} -eq {failed_move} ]]; then return 1; fi; }}; "
        "mv_out_bundle "
        f"{shlex.quote(str(first_source))} {shlex.quote(str(first_destination))} "
        f"{shlex.quote(str(second_source))} {shlex.quote(str(second_destination))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode != 0
    assert first_source.read_text() == "new-first\n"
    assert second_source.read_text() == "new-second\n"
    assert first_destination.read_text() == "old-first\n"
    assert second_destination.read_text() == "old-second\n"
    assert not list(tmp_path.rglob("*.gg-stage.*"))
    assert not list(tmp_path.rglob("*.gg-backup.*"))


def test_mv_out_bundle_removes_partial_cross_filesystem_stage(tmp_path):
    source = tmp_path / "source.txt"
    destination = tmp_path / "out" / "result.txt"
    source.write_text("new\n")
    destination.parent.mkdir()
    destination.write_text("old\n")
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        'mv() { command cp -- "$1" "$2"; return 1; }; '
        "mv_out_bundle "
        f"{shlex.quote(str(source))} {shlex.quote(str(destination))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode != 0
    assert source.read_text() == "new\n"
    assert destination.read_text() == "old\n"
    assert not list(tmp_path.rglob("*.gg-stage.*"))
    assert not list(tmp_path.rglob("*.gg-backup.*"))


@pytest.mark.parametrize("interrupted_move", [1, 3, 5, 6])
def test_mv_out_bundle_recovers_from_signal_after_completed_move(tmp_path, interrupted_move):
    bashpid = subprocess.run(
        ["bash", "-lc", 'printf "%s" "${BASHPID:-}"'],
        capture_output=True,
        text=True,
        check=False,
    )
    if not bashpid.stdout:
        pytest.skip("BASHPID is required to target the mv_out_bundle subshell")
    first_source = tmp_path / "stage" / "first.txt"
    second_source = tmp_path / "stage" / "second.txt"
    first_destination = tmp_path / "out-a" / "first.txt"
    second_destination = tmp_path / "out-b" / "second.txt"
    first_source.parent.mkdir()
    first_destination.parent.mkdir()
    second_destination.parent.mkdir()
    first_source.write_text("new-first\n")
    second_source.write_text("new-second\n")
    first_destination.write_text("old-first\n")
    second_destination.write_text("old-second\n")
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "GG_TEST_MV_COUNT=0; "
        "mv() { GG_TEST_MV_COUNT=$((GG_TEST_MV_COUNT + 1)); "
        'command mv "$@" || return $?; '
        f"if [[ ${{GG_TEST_MV_COUNT}} -eq {interrupted_move} ]]; then "
        'kill -TERM "${BASHPID}"; fi; }; '
        "mv_out_bundle "
        f"{shlex.quote(str(first_source))} {shlex.quote(str(first_destination))} "
        f"{shlex.quote(str(second_source))} {shlex.quote(str(second_destination))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 130
    assert first_source.read_text() == "new-first\n"
    assert second_source.read_text() == "new-second\n"
    assert first_destination.read_text() == "old-first\n"
    assert second_destination.read_text() == "old-second\n"
    assert not list(tmp_path.rglob("*.gg-stage.*"))
    assert not list(tmp_path.rglob("*.gg-backup.*"))


def test_mv_out_bundle_rejects_canonical_destination_aliases(tmp_path):
    first_source = tmp_path / "first.txt"
    second_source = tmp_path / "second.txt"
    first_source.write_text("first\n")
    second_source.write_text("second\n")
    destination = tmp_path / "out" / "result.txt"
    alias = tmp_path / "out" / "." / "result.txt"
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "mv_out_bundle "
        f"{shlex.quote(str(first_source))} {shlex.quote(str(destination))} "
        f"{shlex.quote(str(second_source))} {shlex.quote(str(alias))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode != 0
    assert "duplicate destination" in completed.stdout
    assert first_source.read_text() == "first\n"
    assert second_source.read_text() == "second\n"
    assert not destination.exists()


def test_mv_out_bundle_rejects_hardlinked_destination_aliases(tmp_path):
    first_source = tmp_path / "first.txt"
    second_source = tmp_path / "second.txt"
    first_source.write_text("first\n")
    second_source.write_text("second\n")
    first_destination = tmp_path / "first-destination.txt"
    second_destination = tmp_path / "second-destination.txt"
    first_destination.write_text("old\n")
    second_destination.hardlink_to(first_destination)
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "mv_out_bundle "
        f"{shlex.quote(str(first_source))} {shlex.quote(str(first_destination))} "
        f"{shlex.quote(str(second_source))} {shlex.quote(str(second_destination))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode != 0
    assert "duplicate destination" in completed.stdout
    assert first_source.read_text() == "first\n"
    assert second_source.read_text() == "second\n"


def test_mv_out_bundle_serializes_concurrent_publishers(tmp_path):
    destination = tmp_path / "out" / "result.txt"
    destination.parent.mkdir()
    destination.write_text("old\n")
    first_source = tmp_path / "first.txt"
    second_source = tmp_path / "second.txt"
    first_source.write_text("first\n")
    second_source.write_text("second\n")
    acquired_marker = tmp_path / "first-entered.txt"
    release_marker = tmp_path / "release-first.txt"
    first_command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "GG_LOCK_POLL_SECONDS=1; GG_LOCK_ACQUIRE_TIMEOUT_SECONDS=10; "
        "mkdir() { command mkdir \"$@\"; local status=$?; "
        f"if [[ $* == *gg-bundle.lock* && ${{status}} -eq 0 ]]; then touch {shlex.quote(str(acquired_marker))}; "
        f"while [[ ! -e {shlex.quote(str(release_marker))} ]]; do sleep 0.05; done; fi; "
        "return ${status}; }; "
        f"mv_out_bundle {shlex.quote(str(first_source))} {shlex.quote(str(destination))}"
    )
    second_command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "GG_LOCK_POLL_SECONDS=1; GG_LOCK_ACQUIRE_TIMEOUT_SECONDS=10; "
        f"mv_out_bundle {shlex.quote(str(second_source))} {shlex.quote(str(destination))}"
    )
    first = subprocess.Popen(
        ["bash", "-lc", first_command],
        cwd=tmp_path,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    deadline = time.monotonic() + 20
    while not acquired_marker.exists() and time.monotonic() < deadline:
        time.sleep(0.05)
    if not acquired_marker.exists():
        release_marker.touch()
        first_stdout, first_stderr = first.communicate(timeout=10)
        pytest.fail("first publisher did not acquire lock: " + first_stderr + first_stdout)
    second = subprocess.Popen(
        ["bash", "-lc", second_command],
        cwd=tmp_path,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    time.sleep(0.2)
    assert second.poll() is None
    release_marker.touch()
    first_stdout, first_stderr = first.communicate(timeout=10)
    second_stdout, second_stderr = second.communicate(timeout=10)

    assert first.returncode == 0, first_stderr + first_stdout
    assert second.returncode == 0, second_stderr + second_stdout
    assert destination.read_text() == "second\n"
    assert not first_source.exists()
    assert not second_source.exists()
    assert not list(destination.parent.glob("*.gg-bundle.lock"))
    assert not list(tmp_path.rglob("*.gg-stage.*"))
    assert not list(tmp_path.rglob("*.gg-backup.*"))


def test_mv_out_replace_dir_replaces_existing_nonempty_directory(tmp_path):
    staged_dir = tmp_path / "staged" / "SRR000001"
    dest_dir = tmp_path / "runtime" / "SRR000001"
    (staged_dir / "nested").mkdir(parents=True)
    (staged_dir / "nested" / "new.tsv").write_text("new\n")
    (dest_dir / "nested").mkdir(parents=True)
    (dest_dir / "nested" / "old.tsv").write_text("old\n")
    (dest_dir / "stale.txt").write_text("stale\n")
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"mv_out_replace_dir {shlex.quote(str(staged_dir))} {shlex.quote(str(dest_dir))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert not staged_dir.exists()
    assert not (dest_dir / "stale.txt").exists()
    assert not (dest_dir / "nested" / "old.tsv").exists()
    assert (dest_dir / "nested" / "new.tsv").read_text() == "new\n"


def test_mv_out_replace_dir_restores_existing_directory_when_publish_fails(tmp_path):
    staged_dir = tmp_path / "staged" / "RUN1"
    dest_dir = tmp_path / "runtime" / "RUN1"
    staged_dir.mkdir(parents=True)
    dest_dir.mkdir(parents=True)
    (staged_dir / "new.txt").write_text("new\n")
    (dest_dir / "old.txt").write_text("old\n")
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "GG_TEST_MV_COUNT=0; "
        "mv() { GG_TEST_MV_COUNT=$((GG_TEST_MV_COUNT + 1)); "
        'if [[ ${GG_TEST_MV_COUNT} -eq 3 ]]; then return 9; fi; command mv "$@"; }; '
        f"mv_out_replace_dir {shlex.quote(str(staged_dir))} {shlex.quote(str(dest_dir))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode != 0
    assert (dest_dir / "old.txt").read_text() == "old\n"
    assert not (dest_dir / "new.txt").exists()
    assert (staged_dir / "new.txt").read_text() == "new\n"
    assert not list(tmp_path.rglob("*.gg-stage.*"))
    assert not list(tmp_path.rglob("*.gg-backup.*"))


def test_resolve_rnaspades_transcript_fasta_prefers_primary_output(tmp_path):
    output_dir = tmp_path / "rnaspades_output"
    output_dir.mkdir()
    (output_dir / "transcripts.fasta").write_text(">primary\nAAAA\n")
    (output_dir / "soft_filtered_transcripts.fasta").write_text(">soft\nCCCC\n")
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; resolve_rnaspades_transcript_fasta {shlex.quote(str(output_dir))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(output_dir / "transcripts.fasta")


def test_resolve_rnaspades_transcript_fasta_falls_back_to_soft_then_hard(tmp_path):
    output_dir = tmp_path / "rnaspades_output"
    output_dir.mkdir()
    (output_dir / "soft_filtered_transcripts.fasta").write_text(">soft\nCCCC\n")
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; resolve_rnaspades_transcript_fasta {shlex.quote(str(output_dir))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(output_dir / "soft_filtered_transcripts.fasta")

    (output_dir / "soft_filtered_transcripts.fasta").unlink()
    (output_dir / "hard_filtered_transcripts.fasta").write_text(">hard\nGGGG\n")
    completed = run_bash(command, cwd=tmp_path)
    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(output_dir / "hard_filtered_transcripts.fasta")


def test_capture_busco_repro_artifacts_saves_input_workdir_and_stderr(tmp_path):
    repro_dir = tmp_path / "repro"
    input_fasta = tmp_path / "busco_infile_cdna.fa"
    busco_tmp = tmp_path / "busco_tmp" / "run_eukaryota_odb12"
    stderr_log = tmp_path / "busco.stderr.log"
    input_fasta.write_text(">seq1\nAAAA\n")
    busco_tmp.mkdir(parents=True)
    (busco_tmp / "marker.txt").write_text("marker\n")
    stderr_log.write_text("Failed to open sequence file\n")
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"capture_busco_repro_artifacts "
        f"{shlex.quote(str(repro_dir))} "
        f"{shlex.quote(str(input_fasta))} "
        f"{shlex.quote(str(tmp_path / 'busco_tmp'))} "
        f"eukaryota_odb12 cdna_isoforms "
        f"{shlex.quote(str(stderr_log))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert (repro_dir / "busco_infile_cdna.fa").read_text() == ">seq1\nAAAA\n"
    assert (repro_dir / "busco_tmp" / "run_eukaryota_odb12" / "marker.txt").read_text() == "marker\n"
    assert (repro_dir / "busco.stderr.log").read_text() == "Failed to open sequence file\n"
    info = (repro_dir / "capture_info.tsv").read_text()
    assert "stage_key\tcdna_isoforms\n" in info
    assert "lineage\teukaryota_odb12\n" in info


def test_ensure_jaspar_file_latest_prefers_highest_local_release(tmp_path):
    project_dir = tmp_path / "project"
    jaspar_dir = project_dir / "downloads" / "jaspar"
    jaspar_dir.mkdir(parents=True)
    old_file = jaspar_dir / "JASPAR2022_CORE_plants_non-redundant_pfms_meme.txt"
    new_file = jaspar_dir / "JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt"
    old_file.write_text("old\n")
    new_file.write_text("new\n")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"GG_JASPAR_SKIP_REMOTE_LOOKUP=1; "
        f"ensure_jaspar_file {shlex.quote(str(project_dir))} latest"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(new_file)
    marker = jaspar_dir / "latest_core_plants_non-redundant_pfms_meme.filename"
    assert marker.read_text().strip() == new_file.name


def test_ensure_jaspar_file_latest_uses_cached_marker_without_remote_lookup(tmp_path):
    project_dir = tmp_path / "project"
    jaspar_dir = project_dir / "downloads" / "jaspar"
    jaspar_dir.mkdir(parents=True)
    filename = "JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt"
    target_file = jaspar_dir / filename
    marker = jaspar_dir / "latest_core_plants_non-redundant_pfms_meme.filename"
    target_file.write_text("cached\n")
    marker.write_text(f"{filename}\n")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"GG_JASPAR_SKIP_REMOTE_LOOKUP=1; "
        f"ensure_jaspar_file {shlex.quote(str(project_dir))} auto"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(target_file)


def test_recreate_dir_rejects_root_path(tmp_path):
    command = f"source {shlex.quote(str(GG_UTIL_PATH))}; recreate_dir /"
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode != 0
    assert "Refusing to recreate unsafe directory path: /" in completed.stdout


def test_is_species_set_identical_ignores_hidden_files(tmp_path):
    dir1 = tmp_path / "d1"
    dir2 = tmp_path / "d2"
    dir1.mkdir()
    dir2.mkdir()
    (dir1 / "Arabidopsis_thaliana.fa").write_text(">a\nATG\n")
    (dir2 / "Arabidopsis_thaliana.fa").write_text(">a\nATG\n")
    (dir1 / ".DS_Store").write_text("x\n")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"if is_species_set_identical {shlex.quote(str(dir1))} {shlex.quote(str(dir2))}; then echo OK; else echo NG; fi"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert "OK" in completed.stdout
    assert "NG" not in completed.stdout


def test_gg_species_name_from_path_strips_terminal_extensions(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        'printf "%s\\n" '
        '"$(gg_species_name_from_path Acanthamoeba_castellanii.fa)" '
        '"$(gg_species_name_from_path Arabidopsis_thaliana.fa.gz)" '
        '"$(gg_species_name_from_path Arabidopsis_thaliana_Athaliana_447_Araport11.cds_primaryTranscriptOnly.fa)"'
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip().splitlines() == [
        "Acanthamoeba_castellanii",
        "Arabidopsis_thaliana",
        "Arabidopsis_thaliana",
    ]


def test_is_species_set_identical_ignores_duplicate_busco_name_variants(tmp_path):
    dir1 = tmp_path / "inputs"
    dir2 = tmp_path / "busco"
    dir1.mkdir()
    dir2.mkdir()
    (dir1 / "Acanthamoeba_castellanii.fa").write_text(">a\nATG\n")
    (dir2 / "Acanthamoeba_castellanii_busco.full.tsv").write_text("BUSCO1\tComplete\n")
    (dir2 / "Acanthamoeba_castellanii.fa.busco.full.tsv").write_text("BUSCO1\tComplete\n")

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"if is_species_set_identical {shlex.quote(str(dir1))} {shlex.quote(str(dir2))}; then echo OK; else echo NG; fi"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert "OK" in completed.stdout
    assert "NG" not in completed.stdout


def test_is_species_set_identical_returns_nonzero_when_directory_is_missing(tmp_path):
    dir1 = tmp_path / "d1"
    dir2 = tmp_path / "d2"
    dir1.mkdir()

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"is_species_set_identical {shlex.quote(str(dir1))} {shlex.quote(str(dir2))}"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode != 0
    assert "Directory not found for species-set comparison" in completed.stdout


def test_check_if_species_files_unique_ignores_hidden_files(tmp_path):
    species_dir = tmp_path / "species"
    species_dir.mkdir()
    (species_dir / "Arabidopsis_thaliana.fa").write_text(">a\nATG\n")
    (species_dir / ".DS_Store").write_text("x\n")

    command = f"source {shlex.quote(str(GG_UTIL_PATH))}; check_if_species_files_unique {shlex.quote(str(species_dir))}"
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert "Species files are unique in" in completed.stdout


def test_gg_find_fasta_files_excludes_hidden_files(tmp_path):
    fasta_dir = tmp_path / "fasta"
    fasta_dir.mkdir()
    visible = fasta_dir / "Arabidopsis_thaliana.fa"
    hidden = fasta_dir / ".Arabidopsis_thaliana.fa"
    visible.write_text(">a\nATG\n")
    hidden.write_text(">h\nATG\n")

    command = f"source {shlex.quote(str(GG_UTIL_PATH))}; gg_find_fasta_files {shlex.quote(str(fasta_dir))} 1"
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    out_lines = [line for line in completed.stdout.splitlines() if line.strip()]
    assert str(visible) in out_lines
    assert str(hidden) not in out_lines


def test_gg_find_helpers_follow_symlinked_search_root(tmp_path):
    fasta_dir = tmp_path / "fasta"
    fasta_dir.mkdir()
    visible = fasta_dir / "Arabidopsis_thaliana.fa"
    visible.write_text(">a\nATG\n")
    linked_dir = tmp_path / "linked_fasta"
    linked_dir.symlink_to(fasta_dir, target_is_directory=True)

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"gg_find_fasta_files {shlex.quote(str(linked_dir))} 1; "
        f"gg_find_file_basenames {shlex.quote(str(linked_dir))} '*.fa' 1"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    out_lines = [line for line in completed.stdout.splitlines() if line.strip()]
    assert str(linked_dir / visible.name) in out_lines
    assert visible.name in out_lines


def test_genome_annotation_species_cds_contract_accepts_symlinked_search_root(tmp_path):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    seqkit = bin_dir / "seqkit"
    seqkit.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "[[ ${1:-} == seq ]]\n"
        "shift\n"
        "[[ ${1:-} == --name ]] && shift\n"
        "if [[ ${1:-} == --threads ]]; then shift 2; fi\n"
        'awk \'/^>/ {sub(/^>/, ""); sub(/[[:space:]].*$/, ""); print}\' ${1:?}\n',
        encoding="utf-8",
    )
    seqkit.chmod(0o755)

    species_dir = tmp_path / "species_cds"
    species_dir.mkdir()
    (species_dir / "Arabidopsis_thaliana.fa").write_text(">Arabidopsis_thaliana_gene1\nATG\n")
    (species_dir / "Oryza_sativa.fa").write_text(">Oryza_sativa_gene1\nATG\n")
    linked_dir = tmp_path / "linked_species_cds"
    linked_dir.symlink_to(species_dir, target_is_directory=True)

    find_expression = (
        f"find -H {shlex.quote(str(linked_dir))} -maxdepth 1 -type f ! -name '.*' "
        "\\( -name '*.fa' -o -name '*.fa.gz' -o -name '*.fas' -o -name '*.fas.gz' "
        "-o -name '*.fasta' -o -name '*.fasta.gz' -o -name '*.fna' -o -name '*.fna.gz' \\) | sort"
    )
    command = (
        f"export PATH={shlex.quote(str(bin_dir))}:$PATH; "
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "GG_TASK_CPUS=1; "
        f"{find_expression}; "
        f"check_species_cds_dir {shlex.quote(str(linked_dir) + '/.')}; "
        f"check_if_species_files_unique {shlex.quote(str(linked_dir) + '/.')}"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr + completed.stdout
    assert str(linked_dir / "Arabidopsis_thaliana.fa") in completed.stdout
    assert str(linked_dir / "Oryza_sativa.fa") in completed.stdout
    assert "All per-species CDS files are valid." in completed.stdout
    assert "Species files are unique in" in completed.stdout
