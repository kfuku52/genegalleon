# ruff: noqa: E501

import re

from shell_static_helpers import CORE_DIR, WORKFLOW_DIR
from shell_static_helpers import function_body as _function_body
from shell_static_helpers import read_text as _read_text


def test_transcriptome_core_quotes_known_path_sensitive_options_and_symlinks():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)

    banned_tokens = [
        "--fastq_dir ${dir_species_fastq}",
        "--out_dir ${dir_tmp}",
        "--download_dir ${dir_amalgkit_download_dir}",
        "--download_lock_dir ${dir_amalgkit_download_lock_dir}",
        "--metadata ${file_amalgkit_metadata}",
        "--rrna_filter ${amalgkit_rrna_filter}",
        "--contam_filter ${amalgkit_contam_filter}",
        "--contam_filter_rank ${contamination_removal_rank_for_amalgkit}",
        "--contam_filter_db ${dir_mmseqs2_db}/UniRef90_DB",
        'ln -s ${dir_amalgkit_getfastq_sp} "./getfastq"',
        "--fasta_file ${file_longestcds}",
        "--mmseqs2taxonomy_tsv ${file_longestcds_mmseqs2taxonomy}",
        "--fx2tab_tsv ${file_longestcds_fx2tab}",
        "--species_name ${sp_ub}",
        "--rank ${contamination_removal_rank_for_remove_contaminated_sequences}",
        'seqkit seq --threads ${GG_TASK_CPUS} ${file_isoform} --out-file "busco_infile_cdna.fa"',
        'seqkit seq --threads ${GG_TASK_CPUS} ${file_longestcds} --out-file "busco_infile_cds.fa"',
        'seqkit seq --threads ${GG_TASK_CPUS} ${file_longestcds_contamination_removal_fasta} --out-file "busco_infile_cds.fa"',
        "--lineage_dataset ${dir_busco_lineage}",
        "--download_path ${dir_busco_db}",
        "if [[ -e ${file_kallisto_reference_fasta} ]]; then",
        "ln -s ${file_kallisto_reference_fasta} ${file_reference_fasta_link}",
        "stage_quant_reference_fasta_aliases ${file_amalgkit_metadata} ${file_kallisto_reference_fasta} ./fasta ${sp_ub}",
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
        '--download_lock_dir "${dir_amalgkit_download_lock_dir}"',
        '--metadata "${file_amalgkit_metadata}"',
        '--rrna_filter "${rrna_filter_value}"',
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
        "stage_quant_reference_fasta_aliases \\",
        '        "${file_amalgkit_metadata}" \\',
        '        "${file_kallisto_reference_fasta}" \\',
        '        "./fasta" \\',
        '        "${sp_ub}" \\',
        '        "${gg_support_dir}"',
        'ln -s "${dir_amalgkit_quant}/${sp_ub}" "./quant"',
        "merge_output_prefix=$(resolve_amalgkit_merge_output_prefix \\",
        '    "${file_amalgkit_metadata}" \\',
        '    "./merge" \\',
        '    "${sp_ub}")',
        'mv_out "${merge_output_dir}/${merge_output_prefix}_eff_length.tsv" "${file_amalgkit_merge_efflen}"',
        'mv_out "${merge_output_dir}/${merge_output_prefix}_est_counts.tsv" "${file_amalgkit_merge_count}"',
        'mv_out "${merge_output_dir}/${merge_output_prefix}_tpm.tsv" "${file_amalgkit_merge_tpm}"',
        'grep -F -- "${species_name}" "${metadata_source}"',
        'mv_out "./metadata_private_fastq.tsv" "./metadata.tsv"',
    ]
    for token in expected_tokens:
        assert token in text, f"Missing quoted transcriptome token: {token}"


def test_transcriptome_core_sraid_metadata_filter_handles_zero_match_explicitly():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    extract_body = _function_body(text, "extract_sraid_metadata_rows_for_species")
    assert 'grep -F -- "${species_name}" "${metadata_source}" || true' in extract_body
    assert 'if ! metadata_table_has_data_rows "./metadata.tsv"; then' in text
    assert "No metadata rows matched species" in text


def test_transcriptome_core_sraid_metadata_search_accepts_non_illumina_short_reads():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    assert (
        'amalgkit_sra_strategy_query="${amalgkit_sra_strategy_query:-\\"RNA-seq\\"[Strategy] OR \\"EST\\"[Strategy] OR \\"CLONE\\"[Strategy]}"'
        in text
    )
    assert 'search_string="${search_string} AND (${amalgkit_sra_strategy_query})"' in text
    assert '\\"Illumina\\"[Platform]' not in text
    assert '\\"CLONE\\"[Strategy]' in text


def test_transcriptome_core_relaxes_sraid_strategy_filter_for_missing_explicit_accessions():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    helper = _read_text(WORKFLOW_DIR / "support" / "amalgkit_metadata_accessions.py")
    build_body = _function_body(text, "build_entrez_or_search_string_from_file")
    missing_body = _function_body(text, "extract_requested_accessions_missing_from_metadata")
    relaxed_body = _function_body(text, "extract_transcriptomic_rows_for_requested_accessions")
    merge_body = _function_body(text, "merge_metadata_tables_by_run")

    assert "printf '(%s)\\n' \"${joined_terms}\"" in build_body
    assert (
        'mapfile -t missing_requested_sra_ids < <(extract_requested_accessions_missing_from_metadata "./metadata.tsv" "${file_input_sra_list}")'
        in text
    )
    assert "Retrying the missing accessions without the Entrez strategy filter" in text
    assert 'python "${gg_support_dir}/amalgkit_metadata_accessions.py" \\' in missing_body
    assert "missing \\" in missing_body
    assert 'python "${gg_support_dir}/amalgkit_metadata_accessions.py" \\' in relaxed_body
    assert "extract-transcriptomic \\" in relaxed_body
    assert 'lib_source == "transcriptomic"' in helper
    assert 'lib_strategy in {"rna-seq", "est", "clone"}' in helper
    assert (
        'extract_transcriptomic_rows_for_requested_accessions "./metadata/metadata.tsv" "./metadata_missing_accessions.txt" "./metadata.relaxed.tsv"'
        in text
    )
    assert 'merge_metadata_tables_by_run "./metadata.tsv" "./metadata.relaxed.tsv" "./metadata.merged.tsv"' in text
    assert "Relaxed accession-driven metadata fallback retained" in text
    assert 'python "${gg_support_dir}/amalgkit_metadata_accessions.py" \\' in merge_body
    assert "merge-by-run \\" in merge_body
    assert "Could not determine metadata header" in helper
    assert 'if "run" not in fieldnames:' in helper


def test_transcriptome_core_detects_long_read_platforms_from_metadata():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    body = _function_body(text, "detect_transcriptome_read_technology_from_metadata")
    configure_body = _function_body(text, "configure_transcriptome_runtime_from_detected_metadata")

    assert (
        'file_amalgkit_read_technology="${dir_transcriptome_assembly_output}/amalgkit_read_technology/${sp_ub}_read_technology.tsv"'
        in text
    )
    assert 'file_amalgkit_read_technology_summary_sh="${dir_tmp}/metadata/read_technology.summary.sh"' in text
    assert 'python "${gg_support_dir}/detect_amalgkit_read_technology.py" \\' in body
    assert '--metadata "${metadata_tsv}"' in body
    assert '--classification-out "${classification_tsv}"' in body
    assert '--summary-sh "${summary_sh}"' in body
    assert (
        'detect_transcriptome_read_technology_from_metadata "${file_amalgkit_metadata}" "${file_amalgkit_read_technology}" "${file_amalgkit_read_technology_summary_sh}"'
        in text
    )
    assert 'source "${summary_sh}"' in body
    assert 'effective_assembly_method="rna-bloom2"' in configure_body
    assert 'effective_assembly_method="rnaspades"' in configure_body
    assert "Mixed PacBio and ONT long-read runs were detected in metadata" in configure_body
    assert "Mixed ONT cDNA and direct-RNA runs were detected in metadata" in configure_body
    assert (
        "run_amalgkit_quant remains enabled; amalgkit quant will use quant_backend=${amalgkit_quant_backend}."
        in configure_body
    )
    assert (
        "run_amalgkit_merge remains enabled because amalgkit merge accepts normalized abundance tables from long-read quant."
        in configure_body
    )
    assert "amalgkit quant cannot auto-resolve oarfish sequencing technology for these runs." in configure_body


def test_transcriptome_core_can_recover_public_original_fastqs_after_getfastq_failure():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    body = _function_body(text, "download_public_original_fastqs_for_metadata")

    assert (
        'download_public_original_fastqs_for_metadata "${file_amalgkit_metadata}" "${dir_amalgkit_getfastq_sp}"' in text
    )
    assert (
        "amalgkit getfastq did not safely finish. Attempting fallback download of public original FASTQ files." in text
    )
    assert "Fallback download of public original FASTQ files succeeded." in text
    assert "Fallback direct FASTQ recovery also failed. Exiting." in text
    assert 'xml_url = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc={}".format(' in body
    assert 'if node.attrib.get("semantic_name") != "fastq":' in body
    assert 'if node.attrib.get("supertype") != "Original":' in body
    assert 'if payload[:2] == b"\\x1f\\x8b":' in body
    assert 'dest = run_dir / "{}_{}.amalgkit.fastq.gz".format(run, idx)' in body


def test_transcriptome_core_detects_fatal_getfastq_logs_and_retries_without_rrna():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    cleanup_body = _function_body(text, "cleanup_partial_getfastq_outputs")
    detect_body = _function_body(text, "amalgkit_getfastq_log_has_fatal_message")
    attempt_body = _function_body(text, "run_amalgkit_getfastq_attempt")

    assert 'rm -rf -- "${dir_tmp}/getfastq"' in cleanup_body
    assert 'rm -rf -- "${dir_amalgkit_getfastq_sp}"' in cleanup_body
    assert "grep -Eq '^ERROR: '" in detect_body
    assert "Detected fatal message in amalgkit getfastq log despite a zero exit code" in attempt_body
    assert '--download_lock_dir "${dir_amalgkit_download_lock_dir}"' in attempt_body
    assert '--ncbi_download_max_concurrency "${amalgkit_ncbi_download_max_concurrency}"' in attempt_body
    assert '--aws_download_max_concurrency "${amalgkit_aws_download_max_concurrency}"' in attempt_body
    assert '--gcp_download_max_concurrency "${amalgkit_gcp_download_max_concurrency}"' in attempt_body
    assert 'run_amalgkit_getfastq_attempt "no" "retry_rrna_filter_no"' in text
    assert "Exiting without fallback download so partial outputs do not reach downstream steps." in text


def test_transcriptome_entrypoint_exposes_auto_assembly_and_metadata_detection():
    entrypoint = _read_text(WORKFLOW_DIR / "gg_transcriptome_generation_entrypoint.sh")
    core = _read_text(CORE_DIR / "gg_transcriptome_generation_core.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")

    assert (
        'amalgkit_sra_strategy_query="${amalgkit_sra_strategy_query:-\\"RNA-seq\\"[Strategy] OR \\"EST\\"[Strategy] OR \\"CLONE\\"[Strategy]}" # Entrez strategy clause appended in mode_transcriptome_assembly=sraid; include CLONE so capillary/Sanger cDNA libraries are eligible. Explicit-accession fallback automatically retries without this clause when transcriptomic runs are missed. Set empty to disable strategy filtering.'
        in entrypoint
    )
    assert (
        'amalgkit_sra_strategy_query="${amalgkit_sra_strategy_query:-\\"RNA-seq\\"[Strategy] OR \\"EST\\"[Strategy] OR \\"CLONE\\"[Strategy]}"'
        in core
    )
    assert "amalgkit_sra_strategy_query" in config_vars
    assert (
        'assembly_method="auto" # {auto,Trinity,rnaSPAdes,RNA-Bloom2}; auto picks rnaSPAdes for short-read metadata and RNA-Bloom2 for detected PacBio/ONT metadata.'
        in entrypoint
    )
    assert (
        "requested_assembly_method=$(printf '%s' \"${assembly_method:-auto}\" | tr '[:upper:]' '[:lower:]' | tr '_' '-')"
        in core
    )
    assert "assembly_method" in config_vars
    assert (
        'amalgkit_quant_backend="${amalgkit_quant_backend:-auto}" # {auto,kallisto,oarfish}; auto selects kallisto for short-read runs and oarfish for long-read runs.'
        in entrypoint
    )
    assert (
        'amalgkit_oarfish_seq_tech="${amalgkit_oarfish_seq_tech:-auto}" # {auto,ont-cdna,ont-drna,pac-bio,pac-bio-hifi}; auto infers long-read subtype from metadata for oarfish.'
        in entrypoint
    )
    assert (
        'amalgkit_oarfish_options="${amalgkit_oarfish_options:-}" # Optional extra shell-style option string forwarded to amalgkit quant --oarfish_options.'
        in entrypoint
    )
    assert 'amalgkit_quant_backend="${amalgkit_quant_backend:-auto}"' in core
    assert 'amalgkit_oarfish_seq_tech="${amalgkit_oarfish_seq_tech:-auto}"' in core
    assert 'amalgkit_oarfish_options="${amalgkit_oarfish_options:-}"' in core
    assert "amalgkit_quant_backend" in config_vars
    assert "amalgkit_oarfish_seq_tech" in config_vars
    assert "amalgkit_oarfish_options" in config_vars
    assert "amalgkit_long_read_instrument_pattern" not in entrypoint
    assert "amalgkit_long_read_instrument_pattern" not in core
    assert "amalgkit_long_read_instrument_pattern" not in config_vars


def test_transcriptome_core_long_read_branch_uses_rnabloom2_and_corset():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)

    assert "if [[ \"${effective_assembly_method}\" == 'rna-bloom2' ]]; then" in text
    assert "rnabloom_input_args=(-long)" in text
    assert "rnabloom_extra_args+=(-lrpb)" in text
    assert "rnabloom_extra_args+=(-stranded)" in text
    assert "rnabloom.transcripts.fa" in text
    assert "task='Corset clustering of long-read transcripts'" in text
    assert 'corset_minimap2_preset="map-pb"' in text
    assert 'corset_minimap2_preset="map-ont"' in text
    assert "corset \\" in text
    assert (
        'file_corset_clusters="${dir_transcriptome_assembly_output}/corset_clusters/${sp_ub}_corset.clusters.tsv"'
        in text
    )
    assert 'python "${gg_support_dir}/rename_rnabloom_transcripts.py" \\' in text
    assert 'aggregate_expression="\\-i[0-9].*"' in text


def test_transcriptome_entrypoint_exposes_amalgkit_ncbi_concurrency_limits():
    entrypoint = _read_text(WORKFLOW_DIR / "gg_transcriptome_generation_entrypoint.sh")
    core = _read_text(CORE_DIR / "gg_transcriptome_generation_core.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")

    assert (
        'amalgkit_ncbi_metadata_max_concurrency="${amalgkit_ncbi_metadata_max_concurrency:-20}" # Maximum concurrent NCBI Entrez metadata requests across array tasks. Forwarded to amalgkit metadata/getfastq --ncbi_metadata_max_concurrency. Set 0 or auto to disable throttling.'
        in entrypoint
    )
    assert (
        'amalgkit_ncbi_download_max_concurrency="${amalgkit_ncbi_download_max_concurrency:-20}" # Maximum concurrent NCBI cloud-object downloads across array tasks. Forwarded to amalgkit getfastq --ncbi_download_max_concurrency. Set 0 or auto to disable throttling.'
        in entrypoint
    )
    assert (
        'amalgkit_aws_download_max_concurrency="${amalgkit_aws_download_max_concurrency:-20}" # Maximum concurrent AWS cloud-object downloads across array tasks. Forwarded to amalgkit getfastq --aws_download_max_concurrency. Set 0 or auto to disable throttling.'
        in entrypoint
    )
    assert (
        'amalgkit_gcp_download_max_concurrency="${amalgkit_gcp_download_max_concurrency:-20}" # Maximum concurrent GCP cloud-object downloads across array tasks. Forwarded to amalgkit getfastq --gcp_download_max_concurrency. Set 0 or auto to disable throttling.'
        in entrypoint
    )
    assert 'amalgkit_ncbi_metadata_max_concurrency="${amalgkit_ncbi_metadata_max_concurrency:-20}"' in core
    assert 'amalgkit_ncbi_download_max_concurrency="${amalgkit_ncbi_download_max_concurrency:-20}"' in core
    assert 'amalgkit_aws_download_max_concurrency="${amalgkit_aws_download_max_concurrency:-20}"' in core
    assert 'amalgkit_gcp_download_max_concurrency="${amalgkit_gcp_download_max_concurrency:-20}"' in core
    assert "amalgkit_ncbi_metadata_max_concurrency" in config_vars
    assert "amalgkit_ncbi_download_max_concurrency" in config_vars
    assert "amalgkit_aws_download_max_concurrency" in config_vars
    assert "amalgkit_gcp_download_max_concurrency" in config_vars
    assert "amalgkit_metadata_max_concurrent_jobs" not in config_vars
    assert "amalgkit_getfastq_max_concurrent_jobs" not in config_vars


def test_transcriptome_core_requires_taxid_for_contam_filter():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    assert "effective_amalgkit_contam_filter=" not in text
    assert "Continuing with effective_amalgkit_contam_filter=no." not in text
    assert "amalgkit_contam_filter=yes requires a taxid column in metadata: ${file_amalgkit_metadata}. Exiting." in text


def test_transcriptome_core_passes_download_dir_to_amalgkit_integrate():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    integrate_start = text.index("    amalgkit integrate \\")
    integrate_end = text.index('    mv_out "./metadata_private_fastq.tsv" "./metadata.tsv"', integrate_start)
    integrate_block = text[integrate_start:integrate_end]
    assert '--download_dir "${dir_amalgkit_download_dir}"' in integrate_block
    assert "repair_private_fastq_metadata_scientific_names \\" in integrate_block
    assert '      "./metadata_private_fastq.tsv" \\' in integrate_block
    assert '      "${sp_ub}" \\' in integrate_block
    assert '      "${gg_support_dir}"' in integrate_block


def test_transcriptome_core_passes_shared_mmseqs_db_to_amalgkit_getfastq():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    getfastq_block = _function_body(text, "run_amalgkit_getfastq_attempt")
    assert 'dir_mmseqs2_db="${gg_workspace_downloads_dir}/mmseqs2"' in text
    assert '--contam_filter_db "${dir_mmseqs2_db}/UniRef90_DB"' in getfastq_block


def test_transcriptome_core_delegates_ncbi_parallelism_to_amalgkit():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    normalize_body = _function_body(text, "normalize_amalgkit_download_limit_value")
    require_body = _function_body(text, "require_amalgkit_supported_options")
    metadata_body = _function_body(text, "run_amalgkit_metadata_query")
    getfastq_body = _function_body(text, "run_amalgkit_getfastq_attempt")

    assert "gg_run_with_shared_semaphore" not in text
    assert "dir_amalgkit_metadata_semaphore=" not in text
    assert "dir_amalgkit_getfastq_semaphore=" not in text
    assert 'dir_amalgkit_download_lock_dir="${dir_amalgkit_download_dir}/locks"' in text
    assert (
        'amalgkit_ncbi_metadata_max_concurrency="$(normalize_amalgkit_download_limit_value "${amalgkit_ncbi_metadata_max_concurrency}" "amalgkit_ncbi_metadata_max_concurrency")"'
        in text
    )
    assert (
        'amalgkit_ncbi_download_max_concurrency="$(normalize_amalgkit_download_limit_value "${amalgkit_ncbi_download_max_concurrency}" "amalgkit_ncbi_download_max_concurrency")"'
        in text
    )
    assert (
        'amalgkit_aws_download_max_concurrency="$(normalize_amalgkit_download_limit_value "${amalgkit_aws_download_max_concurrency}" "amalgkit_aws_download_max_concurrency")"'
        in text
    )
    assert (
        'amalgkit_gcp_download_max_concurrency="$(normalize_amalgkit_download_limit_value "${amalgkit_gcp_download_max_concurrency}" "amalgkit_gcp_download_max_concurrency")"'
        in text
    )
    assert "if (( 10#${normalized_value} == 0 )); then" in normalize_body
    assert "printf '0\\n'" in normalize_body
    assert "Use a non-negative integer or auto" in normalize_body
    assert (
        "Installed amalgkit ${subcommand} does not support ${option_name}. Update amalgkit to a build with shared download throttling. Exiting."
        in require_body
    )
    assert (
        'require_amalgkit_supported_options "metadata" "--download_lock_dir" "--ncbi_metadata_max_concurrency"'
        in metadata_body
    )
    assert '--download_lock_dir "${dir_amalgkit_download_lock_dir}"' in metadata_body
    assert '--ncbi_metadata_max_concurrency "${amalgkit_ncbi_metadata_max_concurrency}"' in metadata_body
    assert 'require_amalgkit_supported_options "getfastq" \\' in getfastq_body
    assert '--ncbi_download_max_concurrency "${amalgkit_ncbi_download_max_concurrency}"' in getfastq_body
    assert '--aws_download_max_concurrency "${amalgkit_aws_download_max_concurrency}"' in getfastq_body
    assert '--gcp_download_max_concurrency "${amalgkit_gcp_download_max_concurrency}"' in getfastq_body


def test_transcriptome_core_invalidates_stale_cached_query_tables_on_species_prefix_change():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    body = _function_body(text, "invalidate_cached_query_table_if_prefix_mismatch")

    assert (
        'first_query=$(awk -F \'\\t\' -v skip="${header_lines}" \'NR > skip && $1 != "" { print $1; exit }\' "${table_file}")'
        in body
    )
    assert "expected_species=${expected_prefix%_}" in body
    assert 'first_query_species=$(gg_species_name_from_path_or_dot "${first_query}")' in body
    assert 'if [[ "${first_query_species}" != "${expected_species}" ]]; then' in body
    assert 'stale_file="${table_file}.stale.$(date +%Y%m%d%H%M%S)"' in body
    assert 'mv -f -- "${table_file}" "${stale_file}"' in body
    assert "Archived stale file to: ${stale_file}" in body
    assert (
        'invalidate_cached_query_table_if_prefix_mismatch "${file_longestcds_fx2tab}" "${sp_ub}_" "${task}" 1' in text
    )
    assert (
        'invalidate_cached_query_table_if_prefix_mismatch "${file_longestcds_mmseqs2taxonomy}" "${sp_ub}_" "${task}" 0'
        in text
    )


def test_annotation_and_transcriptome_use_local_contamination_removal_rank_parameter():
    annotation_entrypoint = _read_text(WORKFLOW_DIR / "gg_genome_annotation_entrypoint.sh")
    transcriptome_entrypoint = _read_text(WORKFLOW_DIR / "gg_transcriptome_generation_entrypoint.sh")
    annotation_core = _read_text(CORE_DIR / "gg_genome_annotation_core.sh")
    transcriptome_core = _read_text(CORE_DIR / "gg_transcriptome_generation_core.sh")
    config_vars = _read_text(WORKFLOW_DIR / "support" / "gg_entrypoint_config_vars.sh")
    common = _read_text(WORKFLOW_DIR / "gg_common_params.sh")

    assert (
        'contamination_removal_rank="domain" # Taxonomic rank for contamination removal. Canonical value is domain; GeneGalleon normalizes tool-specific synonyms automatically.'
        in annotation_entrypoint
    )
    assert (
        'contamination_removal_rank="domain" # Taxonomic rank for contamination removal. Canonical value is domain; GeneGalleon normalizes tool-specific synonyms automatically.'
        in transcriptome_entrypoint
    )
    assert (
        'contamination_removal_target_taxon="${contamination_removal_target_taxon:-}" # Optional NCBI taxon name used as the lineage anchor for contamination removal (for example, Eukaryota when the sample species name is unknown).'
        in annotation_entrypoint
    )
    assert (
        'contamination_removal_target_taxon="${contamination_removal_target_taxon:-}" # Optional NCBI taxon name used as the lineage anchor for contamination removal (for example, Eukaryota when the sample species name is unknown).'
        in transcriptome_entrypoint
    )
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

    assert "run_busco_isoforms=1 # BUSCO for transcriptome assembly with isoforms." in entrypoint
    assert "run_busco_longest_cds=1 # BUSCO for longest CDS." in entrypoint
    assert "run_busco_contamination_removed_longest_cds=0 # BUSCO for contamination-removed longest CDS." in entrypoint
    assert 'disable_if_no_input_file "run_busco_isoforms" "${file_isoform}"' in core
    assert 'disable_if_no_input_file "run_busco_longest_cds" "${file_longestcds}"' in core
    assert (
        'disable_if_no_input_file "run_busco_contamination_removed_longest_cds" "${file_longestcds_contamination_removal_fasta}"'
        in core
    )
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


def test_transcriptome_core_quotes_mmseqs_createdb_input_path():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    assert "mmseqs createdb ${file_longestcds} queryDB" not in text
    assert 'mmseqs createdb "${file_longestcds}" queryDB' in text


def test_transcriptome_core_uses_rerun_safe_directory_replacement_for_staged_outputs():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    assert 'mv_out "${dir_tmp}"/getfastq/* "${dir_amalgkit_getfastq_sp}"' not in text
    assert 'mv_out ./quant/* "${dir_amalgkit_quant}/${sp_ub}"' not in text
    assert 'mv_out "./merge/${sp_ub}" "$(dirname "$(dirname "${file_amalgkit_merge_tpm}")")"' not in text
    assert 'mv_out_replace_dir "./merge/${sp_ub}" "$(dirname "${file_amalgkit_merge_tpm}")"' not in text
    assert 'getfastq_outputs=("${dir_tmp}"/getfastq/*)' in text
    assert 'mv_out_replace_dir "${dir_tmp}/getfastq" "${dir_amalgkit_getfastq_sp}"' in text
    assert "quant_outputs=(./quant/*)" in text
    assert 'mv_out_replace_dir "./quant" "${dir_amalgkit_quant}/${sp_ub}"' in text
    assert "resolve_amalgkit_merge_output_prefix" in text
    assert 'mv_out "${merge_output_dir}/${merge_output_prefix}_eff_length.tsv" "${file_amalgkit_merge_efflen}"' in text
    assert 'mv_out "${merge_output_dir}/${merge_output_prefix}_est_counts.tsv" "${file_amalgkit_merge_count}"' in text
    assert 'mv_out "${merge_output_dir}/${merge_output_prefix}_tpm.tsv" "${file_amalgkit_merge_tpm}"' in text


def test_transcriptome_core_uses_array_args_for_trinity_and_rnaspades_inputs():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    assert "trinity_input=" not in text
    assert "${trinity_input}" not in text
    assert '--single "${in_single}"' in text
    assert '--left "${in_left}"' in text
    assert '--right "${in_right}"' in text

    assert "rnaspades_input=$(for i in" not in text
    assert "rnaspades_input_args=()" in text
    assert '"${rnaspades_input_args[@]}"' in text
    assert 'OMP_NUM_THREADS="${assembly_cpus}" \\' in text
    assert 'OMP_THREAD_LIMIT="${assembly_cpus}" \\' in text
    assert 'rnaspades_transcript_fasta=$(resolve_rnaspades_transcript_fasta "${dir_tmp}/rnaspades_output")' in text
    assert 'echo "Using rnaSPAdes transcript fasta: ${rnaspades_transcript_fasta}"' in text
    assert "Checked: transcripts.fasta, soft_filtered_transcripts.fasta, hard_filtered_transcripts.fasta" in text


def test_transcriptome_core_filters_invalid_paired_fastq_before_assembly():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    count_body = _function_body(text, "fastq_num_seqs_from_file")
    filter_body = _function_body(text, "filter_valid_paired_fastq_files")

    assert 'seqkit stats --tabular "${fastq_path}"' in count_body
    assert '$i == "num_seqs"' in count_body
    assert 'right_file="${left_file%_1.amalgkit.fastq.gz}_2.amalgkit.fastq.gz"' in filter_body
    assert 'expected_left="${right_file%_2.amalgkit.fastq.gz}_1.amalgkit.fastq.gz"' in filter_body
    assert "read_count_mismatch" in filter_body
    assert 'files_left=("${valid_left[@]}")' in filter_body
    assert 'files_right=("${valid_right[@]}")' in filter_body
    assert 'if ! filter_valid_paired_fastq_files "${dir_tmp}/paired_fastq_validation.tsv"; then' in text


def test_transcriptome_core_captures_busco_repro_artifacts_on_failure_paths():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)
    assert "capture_busco_failure_context() {" in text
    assert "cleanup_busco_stage_temp_artifacts() {" in text
    assert "run_busco_with_capture() {" in text
    assert "capture_busco_repro_artifacts \\" in text
    assert '2> >(tee "${stderr_log}" >&2)' in text
    assert "if ! gg_run_busco_with_metaeuk_modified_fas_compat \\" in text
    assert 'if gg_busco_stderr_matches_known_metaeuk_modified_fas_bug "${stderr_log}"; then' in text
    assert "return 10" in text
    assert "Skipping BUSCO outputs for longest CDS because BUSCO hit the known MetaEuk transcriptome bug." in text
    assert 'run_busco_with_capture "cdna_isoforms" "busco_infile_cdna.fa"' in text
    assert 'run_busco_with_capture "longest_cds" "busco_infile_cds.fa"' in text
    assert 'run_busco_with_capture "contamination_removed_longest_cds" "busco_infile_cds.fa"' in text
    assert 'capture_busco_failure_context "cdna_isoforms" "busco_infile_cdna.fa" "./busco_tmp.stderr.log"' in text


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


def test_is_fastq_requiring_downstream_analysis_done_quotes_path_checks():
    util_path = WORKFLOW_DIR / "support" / "gg_util.sh"
    text = _read_text(util_path)
    body = _function_body(text, "is_fastq_requiring_downstream_analysis_done")
    assert "-s ${file_isoform}" not in body
    assert "-s ${file_amalgkit_merge_count}" not in body
    assert '-s "${file_isoform}"' in body
    assert '-s "${file_amalgkit_merge_count}"' in body
    assert 'echo "${out}"' in body


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
    assert (
        'if [[ ! -d "${dir_busco}" || -z "$(find "${dir_busco}" -mindepth 1 -print -quit 2> /dev/null)" ]]; then'
        in text
    )


def test_transcriptome_core_guards_array_task_id_range_before_array_indexing():
    script = CORE_DIR / "gg_transcriptome_generation_core.sh"
    text = _read_text(script)

    assert 'dir_input_fastq="${dir_transcriptome_assembly_input}/input_fastq"' not in text
    assert 'dir_input_sra_list="${dir_transcriptome_assembly_input}/input_sra_list"' not in text
    assert "Backward compatibility for legacy input layout" not in text
    assert (
        "mapfile -t fastq_mode_dirs < <(find \"${dir_input_fastq}\" -mindepth 1 -maxdepth 1 -type d ! -name '.*' | sort)"
        in text
    )
    assert (
        'mapfile -t files_fastq < <(find "${dir_species_fastq}" -maxdepth 1 -type f ! -name \'.*\' \\( -name "*.fq" -o -name "*.fastq" -o -name "*.fq.gz" -o -name "*.fastq.gz" \\) | sort)'
        in text
    )
    assert (
        "mapfile -t sra_mode_files < <(find \"${dir_input_sra_list}\" -mindepth 1 -maxdepth 1 -type f ! -name '.*' | sort)"
        in text
    )
    assert (
        "mapfile -t metadata_mode_files < <(find \"${dir_input_amalgkit_metadata}\" -mindepth 1 -maxdepth 1 -type f ! -name '.*' | sort)"
        in text
    )
    assert 'dirs=( "${fastq_mode_dirs[@]}" )' in text
    assert 'files=( "${sra_mode_files[@]}" )' in text
    assert 'files=( "${metadata_mode_files[@]}" )' in text
    assert re.search(r'^[ \t]*dirs=\( "\$\{dir_input_fastq\}"/\* \)', text, re.MULTILINE) is None
    assert re.search(r'^[ \t]*files=\( "\$\{dir_input_sra_list\}"/\* \)', text, re.MULTILINE) is None
    assert re.search(r'^[ \t]*files=\( "\$\{dir_input_amalgkit_metadata\}"/\* \)', text, re.MULTILINE) is None
    assert 'files_fastq=( "${dir_species_fastq}"/* )' not in text
    assert 'if [[ ! "${GG_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${GG_ARRAY_TASK_ID} -lt 1 ]]; then' in text
    assert "if [[ ${#fastq_mode_dirs[@]} -eq 0 ]]; then" in text
    assert "if [[ ${#sra_mode_files[@]} -eq 0 ]]; then" in text
    assert "if [[ ${#metadata_mode_files[@]} -eq 0 ]]; then" in text

    fastq_guard = "if [[ ${id} -ge ${#dirs[@]} ]]; then"
    fastq_use = 'dir_species_fastq="${dirs[${id}]}"'
    assert fastq_guard in text
    assert fastq_use in text
    assert text.index(fastq_guard) < text.index(fastq_use)

    sraid_guard = "if [[ ${id} -ge ${#files[@]} ]]; then"
    sraid_use = 'file_input_sra_list="${files[${id}]}"'
    assert sraid_guard in text
    assert sraid_use in text
    assert text.index(sraid_guard) < text.index(sraid_use)

    metadata_guard = "if [[ ${id} -ge ${#files[@]} ]]; then"
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
