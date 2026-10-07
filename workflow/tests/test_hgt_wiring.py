from pathlib import Path

from shell_static_helpers import read_text

REPO_ROOT = Path(__file__).resolve().parents[2]
GENE_EVOLUTION_ENTRYPOINT = REPO_ROOT / "workflow" / "gg_gene_evolution_entrypoint.sh"
GENE_EVOLUTION_CORE = REPO_ROOT / "workflow" / "core" / "gg_gene_evolution_core.sh"
GENE_SUMMARY_ENTRYPOINT = REPO_ROOT / "workflow" / "gg_gene_summary_entrypoint.sh"
GENE_SUMMARY_CORE = REPO_ROOT / "workflow" / "core" / "gg_gene_summary_core.sh"
GENOME_ANNOTATION_CORE = REPO_ROOT / "workflow" / "core" / "gg_genome_annotation_core.sh"
HGT_CORE = REPO_ROOT / "workflow" / "core" / "gg_hgt_core.sh"


def test_category1_focus_is_default_and_preserves_project_cohort_forwarding():
    entry = GENE_SUMMARY_ENTRYPOINT.read_text()
    core = read_text(GENE_SUMMARY_CORE)
    hgt = read_text(HGT_CORE)
    assert 'run_hgt_trait_focus="${run_hgt_trait_focus:-1}"' in entry
    assert 'run_hgt_focus="${run_hgt_trait_focus}"' in core
    assert 'hgt_focus_event_tsv="${hgt_summary_focus_event_tsv:-auto}"' in core
    assert 'hgt_focus_event_gene_tsv="${hgt_summary_focus_event_gene_tsv:-auto}"' in core
    assert 'python "${gg_support_dir}/focus_hgt_traits.py"' in hgt
    assert '--event_tsv "${hgt_focus_events}" --event_gene_tsv "${hgt_focus_links}"' in hgt
    assert '--output "result_bundle=${dir_hgt_trait_focus}"' in hgt
    assert '--parameter "plots=${run_hgt_plot}"' in hgt
    assert '--gene_family_root "${dir_orthogroup}"' in hgt
    assert '--input-gene-family-subdir "gene_tree_stats=${dir_orthogroup}::stat_branch"' in hgt
    assert '--input "gene_tree_focus_helper=${gg_support_dir}/focus_hgt_gene_trees.py"' in hgt
    assert '--gff_info_root "${gg_workspace_output_dir}/species_gff_info"' in hgt
    assert '--filter_audit_tsv "${hgt_focus_filter_audit_tsv}"' in hgt
    assert 'hgt_summary_focus_context_annotations_tsv="${hgt_summary_focus_context_annotations_tsv:-}"' in entry
    assert 'hgt_focus_context_annotations_tsv="${hgt_summary_focus_context_annotations_tsv:-}"' in core
    assert '--context_annotations_tsv "${hgt_focus_context_annotations_tsv}"' in hgt
    assert '--input "context_annotations=${hgt_focus_context_annotations_tsv}"' in hgt
    assert '--input "gene_context_annotations=${gg_support_dir}/focus_hgt_context_annotations.py"' in hgt
    assert '--mmseqs2_taxonomy_dir "${gg_workspace_output_dir}/species_cds_mmseqs2taxonomy"' in hgt
    assert '--scaffold_taxonomy_dir "${gg_workspace_output_dir}/species_scaffold_taxonomy"' in hgt
    assert '"context_mmseqs2_classifications" "${gg_workspace_output_dir}/species_cds_mmseqs2taxonomy"' in hgt
    assert '"context_gene_host_labels" "${gg_workspace_output_dir}/species_scaffold_taxonomy"' in hgt
    assert 'hgt_summary_focus_context_annotations_tsv' in (REPO_ROOT/'workflow/support/gg_entrypoint_config_vars.sh').read_text()


def test_gene_evolution_records_the_same_arguments_it_renders_for_focused_replay():
    text = read_text(GENE_EVOLUTION_CORE)
    assert 'tree_plot_render_args=(' in text
    assert '--species-parser "${species_label_parser}" -- "${tree_plot_render_args[@]}"' in text
    assert 'Rscript "${gg_support_dir}/stat_branch2tree_plot.r" "${tree_plot_render_args[@]}"' in text
    assert '--output "tree_plot_arguments=${file_og_tree_plot_args}"' in text


def test_shared_query_pfam_filter_defaults_and_provenance_apply_without_plotting():
    entry, summary, hgt = [read_text(path) for path in (GENE_SUMMARY_ENTRYPOINT, GENE_SUMMARY_CORE, HGT_CORE)]
    forwarded = (REPO_ROOT/'workflow/support/gg_entrypoint_config_vars.sh').read_text()
    for suffix, default in [('require_shared_pfam', '1'), ('allow_both_no_pfam', '0'), ('min_shared_pfam_coverage', '0.5')]:
        assert f'hgt_summary_focus_{suffix}="${{hgt_summary_focus_{suffix}:-{default}}}"' in entry
        assert f'hgt_focus_{suffix}="${{hgt_summary_focus_{suffix}:-{default}}}"' in summary
        assert f'--{suffix} "${{hgt_focus_{suffix}}}"' in hgt
        assert f'--parameter "{suffix}=${{hgt_focus_{suffix}}}"' in hgt
        assert f'hgt_summary_focus_{suffix}' in forwarded
    pfam = hgt.index('--input-gene-family-subdir "pfam_query_hits=')
    plotting = hgt.index('if [[ ${run_hgt_plot} -eq 1 ]]; then', pfam)
    assert pfam < plotting
    assert '--input "pfam_filter=${gg_support_dir}/focus_hgt_pfam.py"' in hgt


def test_species_direction_configuration_and_taxonomy_provenance_apply_without_plotting():
    entry, summary, hgt = [read_text(path) for path in (GENE_SUMMARY_ENTRYPOINT, GENE_SUMMARY_CORE, HGT_CORE)]
    forwarded = (REPO_ROOT/'workflow/support/gg_entrypoint_config_vars.sh').read_text()
    for suffix, default in [('direction_filter', 'any'), ('species_taxonomy', 'auto')]:
        assert f'hgt_summary_focus_{suffix}="${{hgt_summary_focus_{suffix}:-{default}}}"' in entry
        assert f'hgt_focus_{suffix}="${{hgt_summary_focus_{suffix}:-{default}}}"' in summary
        assert f'hgt_summary_focus_{suffix}' in forwarded
    assert '--direction_filter "${hgt_focus_direction_filter}"' in hgt
    assert '--species_taxonomy "${hgt_focus_taxonomy_path}"' in hgt
    assert '--input "direction_species_taxonomy=${hgt_focus_taxonomy_path}"' in hgt
    assert '--parameter "direction_filter=${hgt_focus_direction_filter}"' in hgt
    taxonomy = hgt.index('--input "direction_species_taxonomy=')
    assert taxonomy < hgt.index('if [[ ${run_hgt_plot} -eq 1 ]]; then', taxonomy)


def test_gene_evolution_core_passes_uniprot_metadata_and_synteny_to_summary():
    text = read_text(GENE_EVOLUTION_CORE)
    assert '--uniprot_meta_tsv "${uniprot_meta_tsv}"' in text
    assert '--synteny "${file_og_synteny}"' in text
    assert 'summary_untrimmed_fasta="${og_id}.summary.untrimmed.fasta"' in text
    assert (
        'seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_untrimmed_aln_analysis}" '
        '--out-file "${summary_untrimmed_fasta}"' in text
    )
    assert '--untrimmed_aln "${summary_untrimmed_fasta}"' in text
    assert '--trimmed_aln "${summary_trimmed_fasta}"' in text
    assert 'synteny_source_dir="${dir_sp_cds}"' in text
    assert '--input_sequence_mode "${synteny_sequence_mode}"' in text
    assert 'if [[ ${gene_evolution_plot_only} -ne 1 ]] && [[ ${treevis_synteny} -eq 1 || ${treevis_synteny_similarity} -eq 1 ]] && { [[ ${run_summary} -eq 1 ]] || [[ ${run_tree_plot} -eq 1 ]]; }; then' in text


def test_gene_evolution_hgt_profile_is_wired_as_a_preset():
    entry_text = GENE_EVOLUTION_ENTRYPOINT.read_text(encoding="utf-8")
    core_text = read_text(GENE_EVOLUTION_CORE)

    assert 'gene_evolution_profile="${gene_evolution_profile:-default}"' in entry_text
    assert 'input_sequence_mode="${input_sequence_mode:-${GG_COMMON_INPUT_SEQUENCE_MODE:-cds}}"' in entry_text
    assert 'apply_gene_evolution_profile()' in core_text
    assert 'apply_gene_evolution_input_sequence_mode()' in core_text
    assert 'gene_evolution_profile=$(echo "${gene_evolution_profile:-default}"' in core_text
    assert 'input_sequence_mode=$(gg_normalize_input_sequence_mode "${input_sequence_mode}")' in core_text
    assert 'mode_gene_evolution="orthogroup"' in core_text
    assert 'set_profile_default_override run_generax "0" "1"' in core_text
    assert 'set_profile_default_override generax_rec_model "UndatedDL" "UndatedDTL"' in core_text


def test_genome_annotation_core_passes_uniprot_metadata_to_reformatter():
    text = GENOME_ANNOTATION_CORE.read_text(encoding="utf-8")
    assert '--uniprot_meta_tsv "${uniprot_meta_tsv}"' in text


def test_hgt_core_uses_optional_direct_contamination_input_directory():
    entrypoint_text = GENE_SUMMARY_ENTRYPOINT.read_text(encoding="utf-8")
    summary_core_text = GENE_SUMMARY_CORE.read_text(encoding="utf-8")
    core_text = HGT_CORE.read_text(encoding="utf-8")

    assert 'run_hgt_candidate_summary="${run_hgt_candidate_summary:-0}"' in entrypoint_text
    assert 'run_hgt_summary_plots="${run_hgt_summary_plots:-0}"' in entrypoint_text
    assert 'hgt_summary_contamination_dir="${hgt_summary_contamination_dir:-}"' in entrypoint_text
    assert 'hgt_summary_taxonomy_flow_rank="${hgt_summary_taxonomy_flow_rank:-phylum}"' in entrypoint_text
    assert 'hgt_summary_species_tree="${hgt_summary_species_tree:-auto}"' in entrypoint_text
    assert 'hgt_summary_species_trait="${hgt_summary_species_trait:-auto}"' in entrypoint_text
    assert '--species_trait "${hgt_species_trait_path}"' in core_text
    assert 'species_trait.tsv' in core_text
    assert 'hgt_summary_transfer_tree_max_edges="${hgt_summary_transfer_tree_max_edges:-200}"' in entrypoint_text
    assert 'hgt_summary_transfer_arrow_alpha="${hgt_summary_transfer_arrow_alpha:-0.55}"' in entrypoint_text
    assert 'hgt_summary_tree_width_mm="${hgt_summary_tree_width_mm:-60}"' in entrypoint_text
    assert "hgt_min_branch_score" not in entrypoint_text
    assert 'bash "${gg_core_dir}/gg_hgt_core.sh"' in summary_core_text
    assert 'run_hgt_eval="${run_hgt_candidate_summary}"' in summary_core_text
    assert 'run_hgt_plot="${run_hgt_summary_plots}"' in summary_core_text
    assert 'hgt_contamination_dir="${hgt_summary_contamination_dir:-}"' in summary_core_text
    assert 'hgt_species_tree="${hgt_summary_species_tree:-auto}"' in summary_core_text
    assert 'hgt_transfer_tree_max_edges="${hgt_summary_transfer_tree_max_edges:-200}"' in summary_core_text
    assert 'hgt_transfer_arrow_alpha="${hgt_summary_transfer_arrow_alpha:-0.55}"' in summary_core_text
    assert 'run_hgt_plot="${run_hgt_plot:-1}"' in core_text
    assert 'hgt_tree_width_mm="${hgt_tree_width_mm:-60}"' in core_text
    assert 'hgt_species_tree="${hgt_species_tree:-auto}"' in core_text
    assert 'hgt_transfer_tree_max_edges="${hgt_transfer_tree_max_edges:-200}"' in core_text
    assert 'hgt_transfer_arrow_alpha="${hgt_transfer_arrow_alpha:-0.55}"' in core_text
    assert 'hgt_contamination_dir="${hgt_contamination_dir:-}"' in core_text
    assert 'default_hgt_contamination_dir="${gg_workspace_output_dir}/species_cds_contamination_removal_tsv"' in core_text
    assert 'file_hgt_readme="${dir_hgt}/README.md"' in core_text
    assert 'file_hgt_transfer_tree_pdf="${dir_hgt_plot}/hgt_transfer_tree.pdf"' in core_text
    assert 'file_hgt_transfer_edges="${dir_hgt_plot}/hgt_transfer_edges.tsv"' in core_text
    assert '--input "hgt_candidate_scorer=${gg_support_dir}/score_hgt_candidates.py"' in core_text
    assert '--parameter "schema_version=2"' in core_text
    assert 'if [[ -n "${hgt_contamination_dir}" ]]; then' in core_text
    assert '--dir_contamination_tsv "${contamination_arg}"' in core_text
    assert "--min_branch_score" not in core_text
    assert 'python "${gg_support_dir}/plot_hgt_summary.py"' in core_text
    assert '--transfer_tree_pdf "${file_hgt_transfer_tree_pdf}"' in core_text
    assert '--transfer_edges_tsv "${file_hgt_transfer_edges}"' in core_text
    assert '--species_tree "${hgt_species_tree_path}"' in core_text
    assert '--transfer_tree_max_edges "${hgt_transfer_tree_max_edges}"' in core_text
    assert core_text.count('--transfer_arrow_alpha "${hgt_transfer_arrow_alpha}"') == 2
    assert core_text.count('--parameter "transfer_arrow_alpha=${hgt_transfer_arrow_alpha}"') == 2
    assert 'python "${gg_support_dir}/write_hgt_output_readme.py"' in core_text
    assert '--output "${file_hgt_readme}"' in core_text
    assert 'python "${gg_support_dir}/annotate_hgt_tree_plot.py"' in core_text
    assert "mapfile -t hgt_orthogroups" not in core_text
    assert 'while IFS= read -r og_id; do' in core_text
    assert '--panel4="heatmap,no,colrel,_,hgt_,HGT evidence (column max=1)"' in core_text
    assert '--panel15="meme,${file_og_meme}"' in core_text
    assert '--panel8="categorical,besthit_lca_rank_display,Hit LCA,-"' in core_text
    assert '--panel9="signal_peptide"' in core_text
    assert '--panel10="transmembrane_domain"' in core_text
