from pathlib import Path
import shlex
import subprocess


GG_UTIL_PATH = Path(__file__).resolve().parents[1] / "support" / "gg_util.sh"


def run_bash(cmd: str, cwd: Path):
    return subprocess.run(
        ["bash", "-lc", cmd],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )


def test_workspace_pfam_le_dir_is_under_downloads_dedicated_folder(tmp_path):
    project_dir = tmp_path / "project"
    project_dir.mkdir()
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"workspace_pfam_le_dir {shlex.quote(str(project_dir))}"
    )

    completed = run_bash(command, cwd=tmp_path)
    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == str(project_dir / "downloads" / "pfam" / "Pfam_LE")


def test_gg_array_download_once_accepts_nonempty_ready_marker(tmp_path):
    lock_file = tmp_path / "locks" / "artifact.lock"
    marker_file = tmp_path / "runtime" / ".artifact.ready"
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"marker_file={shlex.quote(str(marker_file))}; "
        f'gg_array_download_once {shlex.quote(str(lock_file))} "$marker_file" "marker artifact" '
        'gg_write_ready_marker "$marker_file"; '
        'status=$?; '
        'printf "%s\\n" "$status"; '
        'wc -c < "$marker_file"; '
        'cat "$marker_file"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    lines = completed.stdout.strip().splitlines()
    assert lines[0] == "0"
    assert int(lines[1]) > 0
    assert lines[2] == "ready"


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
        f'{shlex.quote(str(runtime_db))} '
        f'{shlex.quote(str(runtime_lineage))} '
        f'{shlex.quote(str(ready_marker))}; '
        'status=$?; '
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
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "gg_normalize_contamination_removal_rank_for_amalgkit superkingdom"
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
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        'gg_resolve_annotation_species "Arabidopsis_thaliana_"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "Arabidopsis_thaliana"


def test_resolve_busco_lineage_from_lineages_prefers_deepest_mapped_taxon(tmp_path):
    mapping_dir = tmp_path / "busco_mappings"
    mapping_dir.mkdir()
    (
        mapping_dir / "mapping_taxids-busco_dataset_name.eukaryota_odb12.test.txt"
    ).write_text(
        "2759\teukaryota_odb12\n"
        "33090\tviridiplantae_odb12\n"
        "3193\tembryophyta_odb12\n"
        "3744\trosales_odb12\n"
    )

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f'gg_resolve_busco_lineage_from_lineages auto {shlex.quote(str(mapping_dir))} '
        '"1,131567,2759,33090,3193,3744"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "rosales_odb12"


def test_resolve_busco_lineage_from_lineages_uses_deepest_common_taxon(tmp_path):
    mapping_dir = tmp_path / "busco_mappings"
    mapping_dir.mkdir()
    (
        mapping_dir / "mapping_taxids-busco_dataset_name.eukaryota_odb12.test.txt"
    ).write_text(
        "2759\teukaryota_odb12\n"
        "33090\tviridiplantae_odb12\n"
        "3193\tembryophyta_odb12\n"
        "3700\tbrassicales_odb12\n"
        "4530\tpoales_odb12\n"
    )

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f'gg_resolve_busco_lineage_from_lineages auto {shlex.quote(str(mapping_dir))} '
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
        f'gg_resolve_busco_lineage_from_lineages metazoa_odb13 {shlex.quote(str(mapping_dir))} '
        '"1,131567,2759,33208"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "metazoa_odb13"


def test_resolve_busco_lineage_from_lineages_prefers_latest_common_odb_version(tmp_path):
    mapping_dir = tmp_path / "busco_mappings"
    mapping_dir.mkdir()
    (
        mapping_dir / "mapping_taxids-busco_dataset_name.archaea_odb13.test.txt"
    ).write_text("2157\tarchaea_odb13\n")
    (
        mapping_dir / "mapping_taxids-busco_dataset_name.bacteria_odb13.test.txt"
    ).write_text("2\tbacteria_odb13\n")
    (
        mapping_dir / "mapping_taxids-busco_dataset_name.eukaryota_odb12.test.txt"
    ).write_text(
        "2759\teukaryota_odb12\n"
        "33090\tviridiplantae_odb12\n"
        "3193\tembryophyta_odb12\n"
    )
    (
        mapping_dir / "mapping_taxids-busco_dataset_name.eukaryota_odb13.test.txt"
    ).write_text(
        "2759\teukaryota_odb13\n"
        "33090\tviridiplantae_odb13\n"
        "3193\tembryophyta_odb13\n"
        "3744\trosales_odb13\n"
    )

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f'gg_resolve_busco_lineage_from_lineages auto {shlex.quote(str(mapping_dir))} '
        '"1,131567,2759,33090,3193,3744"'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "rosales_odb13"


def test_finalize_auto_busco_lineage_name_appends_requested_odb_suffix(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        'gg_finalize_auto_busco_lineage_name brassicales 13'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "brassicales_odb13"


def test_finalize_auto_busco_lineage_name_preserves_existing_suffix(tmp_path):
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        'gg_finalize_auto_busco_lineage_name embryophyta_odb12'
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "embryophyta_odb12"


def test_workspace_layout_defaults_to_split_for_empty_workspace(tmp_path):
    project_dir = tmp_path / "project"
    project_dir.mkdir()
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"gg_prepare_cmd_runtime {shlex.quote(str(project_dir))} \"\" 0 0; "
        "printf '%s\\n%s\\n%s\\n%s\\n' "
        "\"${gg_workspace_layout_resolved}\" "
        "\"${gg_workspace_input_dir}\" "
        "\"${gg_workspace_output_dir}\" "
        "\"${gg_workspace_downloads_dir}\""
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
        f"\"$(gg_resolve_workspace_layout {shlex.quote(str(project_dir))})\" "
        f"\"$(workspace_input_root {shlex.quote(str(project_dir))})\" "
        f"\"$(workspace_output_root {shlex.quote(str(project_dir))})\""
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
        f"gg_prepare_cmd_runtime {shlex.quote(str(project_dir))} \"\" 0 0; "
        "printf '%s\\n%s\\n%s\\n%s\\n' "
        "\"${gg_workspace_layout_resolved}\" "
        "\"${gg_workspace_input_dir}\" "
        "\"${gg_workspace_output_dir}\" "
        "\"${gg_workspace_downloads_dir}\""
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

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"ensure_pfam_le_db {shlex.quote(str(project_dir))}"
    )
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
        'status=$?; '
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
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"printf 'hello\\n' | mv_out {shlex.quote(str(outfile))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert outfile.read_text() == "hello\n"


def test_cp_out_accepts_pipe_input(tmp_path):
    outfile = tmp_path / "nested" / "out.txt"
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"printf 'world\\n' | cp_out {shlex.quote(str(outfile))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert outfile.read_text() == "world\n"


def test_cp_out_single_argument_without_pipe_fails(tmp_path):
    outfile = tmp_path / "out.txt"
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"cp_out {shlex.quote(str(outfile))}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode != 0
    assert "at least 2 arguments are required unless stdin is piped" in completed.stdout


def test_cp_out_creates_destination_dir_when_target_has_trailing_slash(tmp_path):
    src = tmp_path / "src.txt"
    src.write_text("abc\n")
    dest_dir = tmp_path / "nested" / "dest_dir"
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"cp_out {shlex.quote(str(src))} {shlex.quote(str(dest_dir) + '/')}"
    )

    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    copied = dest_dir / src.name
    assert copied.read_text() == "abc\n"


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
    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        "recreate_dir /"
    )
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

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"check_if_species_files_unique {shlex.quote(str(species_dir))}"
    )
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

    command = (
        f"source {shlex.quote(str(GG_UTIL_PATH))}; "
        f"gg_find_fasta_files {shlex.quote(str(fasta_dir))} 1"
    )
    completed = run_bash(command, cwd=tmp_path)

    assert completed.returncode == 0, completed.stderr
    out_lines = [line for line in completed.stdout.splitlines() if line.strip()]
    assert str(visible) in out_lines
    assert str(hidden) not in out_lines
