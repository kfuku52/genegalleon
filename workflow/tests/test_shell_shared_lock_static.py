from shell_static_helpers import CORE_DIR, WORKFLOW_DIR, function_body, read_text


def test_shared_lock_helpers_delegate_metadata_to_shared_python_module():
    util_text = read_text(WORKFLOW_DIR / "support" / "gg_util.sh")
    helper_text = read_text(WORKFLOW_DIR / "support" / "shared_lock.py")

    assert 'SHARED_LOCK_FORMAT = "shared-lock-v2"' in helper_text
    assert "os.O_CREAT | os.O_EXCL | os.O_WRONLY" in helper_text
    assert "def stale_lock_reason(" in helper_text
    assert "same_host_same_boot_dead_pid" in helper_text
    assert "heartbeat_timeout" in helper_text

    assert "gg_shared_lock_helper_script()" in util_text
    assert '"${helper_script}" read-metadata "${lock_file}"' in util_text
    assert '"${helper_script}" owner-summary "${lock_file}"' in util_text
    assert '"${helper_script}" try-create "${lock_file}" --pid "${owner_pid}"' in util_text
    assert '"${helper_script}" reclaim-if-stale "${lock_file}" --stale-seconds "${stale_seconds}"' in util_text
    assert 'touch -c -- "${lock_file}" 2>/dev/null || true' in util_text
    assert "waiting for shared lock: ${description} (${owner_summary})" in util_text
    assert "timed out waiting for shared lock: ${description} (${owner_summary})" in util_text


def test_shared_semaphore_helpers_reuse_shared_lock_primitives():
    util_text = read_text(WORKFLOW_DIR / "support" / "gg_util.sh")
    acquire_body = function_body(util_text, "gg_shared_semaphore_acquire")

    assert 'slot_lock="${semaphore_dir}/slot.${slot_idx}.lock"' in acquire_body
    assert 'if gg_shared_lock_try_create "${slot_lock}"; then' in acquire_body
    assert (
        'if gg_shared_lock_reclaim_if_stale "${slot_lock}" "${description} slot ${slot_idx}/${max_slots}"; then'
    ) in acquire_body
    assert 'GG_SHARED_SEMAPHORE_SLOT_LOCK_FILE="${slot_lock}"' in acquire_body
    assert (
        "waiting for shared semaphore slot: ${description} (max_concurrent=${max_slots}; lock_dir=${semaphore_dir})"
    ) in acquire_body
    assert (
        "timed out waiting for shared semaphore slot: ${description} "
        "(max_concurrent=${max_slots}; lock_dir=${semaphore_dir})"
    ) in acquire_body
    assert 'gg_shared_lock_start_heartbeat "${slot_lock}"' in util_text
    assert "trap 'cleanup_shared_semaphore; restore_shared_semaphore_traps' EXIT" in util_text
    assert "trap 'shared_semaphore_signal_handler INT' INT" in util_text
    assert "trap 'shared_semaphore_signal_handler TERM' TERM" in util_text
    assert 'gg_shared_lock_stop_heartbeat "${heartbeat_pid}"' in util_text
    assert 'gg_shared_semaphore_release "${slot_lock}"' in util_text


def test_gene_evolution_core_uses_shared_lock_helpers_for_db_builds_and_shared_copies():
    script = CORE_DIR / "gg_gene_evolution_core.sh"
    text = read_text(script)
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
    text = read_text(script)
    assert "command -v flock" not in text
    assert "species_cds_validation_lock_dir" not in text
    assert 'gg_shared_lock_acquire "${species_cds_validation_lock}" "species CDS validation stamp"' in text
    assert 'gg_shared_lock_start_heartbeat "${species_cds_validation_lock}"' in text
    assert "heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}" in text
