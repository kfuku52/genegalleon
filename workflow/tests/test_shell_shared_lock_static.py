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
