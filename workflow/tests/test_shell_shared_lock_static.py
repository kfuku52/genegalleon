from shell_static_helpers import CORE_DIR, read_text


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
