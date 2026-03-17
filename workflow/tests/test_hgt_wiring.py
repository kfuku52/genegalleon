from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
GENE_EVOLUTION_ENTRYPOINT = REPO_ROOT / "workflow" / "gg_gene_evolution_entrypoint.sh"
GENE_EVOLUTION_CORE = REPO_ROOT / "workflow" / "core" / "gg_gene_evolution_core.sh"
GENOME_ANNOTATION_CORE = REPO_ROOT / "workflow" / "core" / "gg_genome_annotation_core.sh"
HGT_ENTRYPOINT = REPO_ROOT / "workflow" / "gg_hgt_entrypoint.sh"
HGT_CORE = REPO_ROOT / "workflow" / "core" / "gg_hgt_core.sh"


def test_gene_evolution_core_passes_uniprot_metadata_and_synteny_to_summary():
    text = GENE_EVOLUTION_CORE.read_text(encoding="utf-8")
    assert '--uniprot_meta_tsv "${uniprot_meta_tsv}"' in text
    assert '--synteny "${file_og_synteny}"' in text
    assert 'summary_unaligned_fasta="${og_id}.summary.unaligned.fasta"' in text
    assert 'seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_primary_fasta}" --out-file "${summary_unaligned_fasta}"' in text
    assert 'synteny_source_dir="${dir_sp_cds}"' in text
    assert '--input_sequence_mode "${synteny_sequence_mode}"' in text
    assert 'if [[ ${treevis_synteny} -eq 1 ]] && { [[ ${run_summary} -eq 1 ]] || [[ ${run_tree_plot} -eq 1 ]]; }; then' in text


def test_gene_evolution_hgt_profile_is_wired_as_a_preset():
    entry_text = GENE_EVOLUTION_ENTRYPOINT.read_text(encoding="utf-8")
    core_text = GENE_EVOLUTION_CORE.read_text(encoding="utf-8")

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
    entrypoint_text = HGT_ENTRYPOINT.read_text(encoding="utf-8")
    core_text = HGT_CORE.read_text(encoding="utf-8")

    assert 'hgt_contamination_dir=""' in entrypoint_text
    assert 'hgt_contamination_dir="${hgt_contamination_dir:-}"' in core_text
    assert 'default_hgt_contamination_dir="${gg_workspace_output_dir}/species_cds_contamination_removal_tsv"' in core_text
    assert 'if [[ -n "${hgt_contamination_dir}" ]]; then' in core_text
    assert '--dir_contamination_tsv "${contamination_arg}"' in core_text
