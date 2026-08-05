#!/usr/bin/env bash
set -euo pipefail

gg_core_bootstrap="/script/support/gg_core_bootstrap.sh"
if [[ ! -s "${gg_core_bootstrap}" ]]; then
  gg_core_bootstrap="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../support/gg_core_bootstrap.sh"
fi
# shellcheck disable=SC1090
source "${gg_core_bootstrap}"
unset gg_core_bootstrap

### Start: Job-supplied configuration ###
# Configuration variables are provided by gg_gene_summary_entrypoint.sh.
### End: Job-supplied configuration ###

gg_bootstrap_core_runtime "${BASH_SOURCE[0]:-$0}" "base" 1 1

hgt_gene_family_mode="${hgt_gene_family_mode:-orthogroup}"

run_hgt_eval="${run_hgt_eval:-1}"
run_hgt_plot="${run_hgt_plot:-1}"
hgt_use_taxonomy_db="${hgt_use_taxonomy_db:-1}"
hgt_contamination_dir="${hgt_contamination_dir:-}"
hgt_taxonomy_flow_rank="${hgt_taxonomy_flow_rank:-phylum}"
hgt_taxonomy_flow_max_categories="${hgt_taxonomy_flow_max_categories:-12}"
hgt_tree_plot_width="${hgt_tree_plot_width:-24}"
hgt_promoter_bp="${hgt_promoter_bp:-2000}"
hgt_fimo_qvalue="${hgt_fimo_qvalue:-0.05}"

dir_orthogroup="${hgt_gene_family_dir:-${gg_workspace_output_dir}/orthogroup}"
file_orthogroup_db="${hgt_db_path:-${dir_orthogroup}/gg_orthogroup.db}"
dir_hgt="${hgt_output_dir:-${gg_workspace_output_dir}/hgt}"
default_hgt_contamination_dir="${gg_workspace_output_dir}/species_cds_contamination_removal_tsv"
file_hgt_branch="${dir_hgt}/hgt_branch_candidates.tsv"
file_hgt_gene="${dir_hgt}/hgt_gene_candidates.tsv"
file_hgt_orthogroup="${dir_hgt}/hgt_orthogroup_summary.tsv"
dir_hgt_plot="${dir_hgt}/plots"
dir_hgt_tree_plot="${dir_hgt}/tree_plot"
dir_hgt_tree_input="${dir_hgt}/tree_plot_input"
dir_hgt_tmp="${dir_hgt}/tmp"
file_hgt_overview_pdf="${dir_hgt_plot}/hgt_branch_overview.pdf"
file_hgt_taxonomy_flow_pdf="${dir_hgt_plot}/hgt_taxonomy_flow.pdf"

enable_all_run_flags_for_debug_mode

if [[ "${run_hgt_eval}" != "0" && "${run_hgt_eval}" != "1" ]]; then
  echo "Invalid binary flag value: run_hgt_eval=${run_hgt_eval} (expected 0 or 1)"
  exit 1
fi
if [[ "${run_hgt_plot}" != "0" && "${run_hgt_plot}" != "1" ]]; then
  echo "Invalid binary flag value: run_hgt_plot=${run_hgt_plot} (expected 0 or 1)"
  exit 1
fi
if [[ "${hgt_use_taxonomy_db}" != "0" && "${hgt_use_taxonomy_db}" != "1" ]]; then
  echo "Invalid binary flag value: hgt_use_taxonomy_db=${hgt_use_taxonomy_db} (expected 0 or 1)"
  exit 1
fi
if ! [[ "${hgt_taxonomy_flow_max_categories}" =~ ^[0-9]+$ ]]; then
  echo "Invalid hgt_taxonomy_flow_max_categories: ${hgt_taxonomy_flow_max_categories}"
  exit 1
fi
if ! [[ "${hgt_tree_plot_width}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "Invalid hgt_tree_plot_width: ${hgt_tree_plot_width}"
  exit 1
fi
if ! [[ "${hgt_promoter_bp}" =~ ^[0-9]+$ ]]; then
  echo "Invalid hgt_promoter_bp: ${hgt_promoter_bp}"
  exit 1
fi
if ! [[ "${hgt_fimo_qvalue}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "Invalid hgt_fimo_qvalue: ${hgt_fimo_qvalue}"
  exit 1
fi
mkdir -p "${dir_hgt}"

hgt_majority_ortholog_prefix() {
  local orthogroup_id=$1
  local gene_tsv=$2
  if [[ ! -s "${gene_tsv}" ]]; then
    return 0
  fi
  awk -F $'\t' -v og="${orthogroup_id}" '
      NR==1 {
        for (i=1; i<=NF; i++) {
          if ($i=="orthogroup") og_col=i
          if ($i=="gene_taxon") tax_col=i
        }
        next
      }
      (og_col>0 && tax_col>0 && $og_col==og && $tax_col!="") {
        counts[$tax_col]++
      }
      END {
        max_count=0
        best=""
        for (taxon in counts) {
          if (counts[taxon] > max_count || (counts[taxon] == max_count && taxon < best)) {
            max_count=counts[taxon]
            best=taxon
          }
        }
        gsub(/ /, "_", best)
        if (best != "") {
          print best "_"
        }
      }
    ' "${gene_tsv}"
}

hgt_count_tree_tips() {
  local stat_branch=$1
  awk -F $'\t' '
      NR==1 {
        for (i=1; i<=NF; i++) {
          if ($i=="so_event") so_col=i
          if ($i=="num_leaf") leaf_col=i
        }
        next
      }
      ((so_col>0 && $so_col=="L") || (leaf_col>0 && $leaf_col=="1")) {n++}
      END {print n+0}
    ' "${stat_branch}"
}

hgt_select_alignment_input() {
  local num_tip=$1
  shift
  local candidate=""
  local candidate_n=0
  for candidate in "$@"; do
    candidate_n=$(gg_count_fasta_records "${candidate}")
    if [[ ${candidate_n} -ge ${num_tip} ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  printf '%s\n' "$1"
}

hgt_materialization_root="${dir_hgt_tmp}/materialized"
hgt_materialization_global_lock="${dir_hgt_tmp}/materialized.cleanup.lock"
hgt_materialization_run_dir=""
hgt_materialization_lock_held=0
hgt_materialization_locking_available=0
if command -v flock >/dev/null 2>&1; then
  hgt_materialization_locking_available=1
fi

hgt_cleanup_materialization_run() {
  local target="${hgt_materialization_run_dir:-}"
  if [[ -n "${target}" && ${hgt_materialization_locking_available} -eq 1 ]]; then
    mkdir -p -- "${dir_hgt_tmp}"
    exec 198> "${hgt_materialization_global_lock}"
    flock -x 198
  fi
  if [[ ${hgt_materialization_lock_held:-0} -eq 1 ]]; then
    flock -u 197 2>/dev/null || true
    exec 197>&-
    hgt_materialization_lock_held=0
  fi
  if [[ -z "${target}" ]]; then
    return 0
  fi
  case "${target}" in
    "${hgt_materialization_root}/"*)
      rm -rf -- "${target}"
      rmdir -- "${hgt_materialization_root}" 2>/dev/null || true
      ;;
    *)
      echo "Refusing to remove unexpected HGT materialization run: ${target}" >&2
      return 1
      ;;
  esac
  hgt_materialization_run_dir=""
  if [[ ${hgt_materialization_locking_available} -eq 1 ]]; then
    flock -u 198 2>/dev/null || true
    exec 198>&-
  fi
}

hgt_cleanup_stale_materialization_runs() {
  local candidate=""
  if [[ ! -d "${hgt_materialization_root}" || -L "${hgt_materialization_root}" ]]; then
    return 0
  fi
  for candidate in "${hgt_materialization_root}"/*; do
    [[ -d "${candidate}" && ! -L "${candidate}" ]] || continue
    exec 196> "${candidate}/.run.lock"
    if flock -n -x 196; then
      case "${candidate}" in
        "${hgt_materialization_root}/"*)
          rm -rf -- "${candidate}"
          ;;
      esac
      flock -u 196 2>/dev/null || true
    fi
    exec 196>&-
  done
  rmdir -- "${hgt_materialization_root}" 2>/dev/null || true
}

hgt_prepare_materialization_run() {
  mkdir -p -- "${dir_hgt_tmp}"
  if [[ ${hgt_materialization_locking_available} -eq 0 ]]; then
    mkdir -p "${hgt_materialization_root}"
    hgt_materialization_run_dir="${hgt_materialization_root}/${GG_JOB_ID:-local}_${GG_ARRAY_TASK_ID:-1}_$$"
    mkdir -p "${hgt_materialization_run_dir}"
    return 0
  fi
  exec 198> "${hgt_materialization_global_lock}"
  flock -x 198
  hgt_cleanup_stale_materialization_runs
  mkdir -p "${hgt_materialization_root}"
  hgt_materialization_run_dir="${hgt_materialization_root}/${GG_JOB_ID:-local}_${GG_ARRAY_TASK_ID:-1}_$$"
  mkdir -p "${hgt_materialization_run_dir}"
  exec 197> "${hgt_materialization_run_dir}/.run.lock"
  flock -x 197
  hgt_materialization_lock_held=1
  flock -u 198
  exec 198>&-
}

hgt_cleanup_family_materialization() {
  local target="${1:-}"
  [[ -n "${target}" ]] || return 0
  case "${target}" in
    "${hgt_materialization_run_dir}/"*)
      rm -rf -- "${target}"
      ;;
    *)
      echo "Refusing to remove unexpected HGT family materialization: ${target}" >&2
      return 1
      ;;
  esac
}

trap hgt_cleanup_materialization_run EXIT

if [[ ${run_hgt_eval} -eq 1 ]]; then
  if [[ ! -s "${file_orthogroup_db}" ]]; then
    echo "Skipping HGT evaluation because the orthogroup database was not found: ${file_orthogroup_db}"
    exit 0
  fi

  hgt_taxonomy_dbfile=""
  if [[ ${hgt_use_taxonomy_db} -eq 1 ]]; then
    if ensure_ete_taxonomy_db "${gg_workspace_dir}"; then
      hgt_taxonomy_dbfile=$(workspace_taxonomy_dbfile "${gg_workspace_dir}")
    else
      echo "Warning: Failed to prepare the ETE taxonomy DB. Continuing with best-hit name heuristics only." >&2
    fi
  fi

  contamination_arg=""
  if [[ -n "${hgt_contamination_dir}" ]]; then
    if [[ -d "${hgt_contamination_dir}" ]]; then
      contamination_arg="${hgt_contamination_dir}"
    else
      echo "Warning: HGT contamination directory was provided but not found. Skipping contamination input: ${hgt_contamination_dir}" >&2
    fi
  elif [[ -d "${default_hgt_contamination_dir}" ]]; then
    contamination_arg="${default_hgt_contamination_dir}"
  fi

  python "${gg_support_dir}/score_hgt_candidates.py" \
    --dbpath "${file_orthogroup_db}" \
    --branch_out "${file_hgt_branch}" \
    --gene_out "${file_hgt_gene}" \
    --orthogroup_out "${file_hgt_orthogroup}" \
    --dir_contamination_tsv "${contamination_arg}" \
    --taxonomy_dbfile "${hgt_taxonomy_dbfile}"
fi

if [[ ${run_hgt_plot} -eq 1 ]]; then
  if [[ ! -s "${file_hgt_branch}" || ! -s "${file_hgt_gene}" ]]; then
    echo "Skipping HGT plotting because candidate tables were not found: ${file_hgt_branch}, ${file_hgt_gene}"
    echo "$(date): Exiting Singularity environment"
    exit 0
  fi

  mkdir -p "${dir_hgt_plot}" "${dir_hgt_tree_plot}" "${dir_hgt_tree_input}" "${dir_hgt_tmp}"
  hgt_prepare_materialization_run

  hgt_taxonomy_dbfile=""
  if [[ ${hgt_use_taxonomy_db} -eq 1 ]] && ensure_ete_taxonomy_db "${gg_workspace_dir}"; then
    hgt_taxonomy_dbfile=$(workspace_taxonomy_dbfile "${gg_workspace_dir}")
  fi

  python "${gg_support_dir}/plot_hgt_summary.py" \
    --branch_tsv "${file_hgt_branch}" \
    --gene_tsv "${file_hgt_gene}" \
    --overview_pdf "${file_hgt_overview_pdf}" \
    --taxonomy_flow_pdf "${file_hgt_taxonomy_flow_pdf}" \
    --taxonomy_dbfile "${hgt_taxonomy_dbfile}" \
    --flow_rank "${hgt_taxonomy_flow_rank}" \
    --flow_max_categories "${hgt_taxonomy_flow_max_categories}"

  if ! Rscript -e "if (!requireNamespace('ggimage', quietly=TRUE)) quit(status=1)" > /dev/null 2>&1; then
    echo "ggimage package is unavailable. Skipping HGT tree plots."
    echo "$(date): Exiting Singularity environment"
    exit 0
  fi

  hgt_orthogroups=()
  while IFS= read -r og_id; do
    [[ -z "${og_id}" ]] && continue
    hgt_orthogroups+=("${og_id}")
  done < <(
    awk -F $'\t' '
        NR==1 {
          for (i=1; i<=NF; i++) {
            if ($i=="orthogroup") col=i
          }
          next
        }
        (col>0 && $col!="") {seen[$col]=1}
        END {
          for (k in seen) print k
        }
      ' "${file_hgt_branch}" | sort
  )

  for og_id in "${hgt_orthogroups[@]}"; do
    [[ -z "${og_id}" ]] && continue
    dir_og_input_root="${dir_orthogroup}"
    dir_hgt_materialized="${hgt_materialization_run_dir}/${og_id}"
    hgt_current_materialized=""
    if [[ -d "${dir_orthogroup}/.gg_store" || -d "${dir_orthogroup}/.gg_archives" ]]; then
      hgt_current_materialized="${dir_hgt_materialized}"
      materialize_args=(
        materialize-family
        --root "${dir_orthogroup}"
        --mode "${hgt_gene_family_mode}"
        --family-id "${og_id}"
        --destination-root "${dir_hgt_materialized}"
        --subdirs "stat_branch,synteny,rpsblast,clipkit,orthogroup_extraction_fasta,maxalign,mafft,protein_fasta,cds_fasta,dated_tree,fimo,meme"
      )
      if [[ "${hgt_gene_family_mode}" == "query2family" ]]; then
        materialize_args+=(--query-dir "${gg_workspace_input_dir}/query_gene")
      fi
      if ! python "${gg_support_dir}/gene_family_output_store.py" "${materialize_args[@]}"; then
        echo "Skipping HGT tree plot for ${og_id}: failed to materialize ZIP-backed inputs." >&2
        hgt_cleanup_family_materialization "${hgt_current_materialized}"
        hgt_current_materialized=""
        continue
      fi
      dir_og_input_root="${dir_hgt_materialized}"
    fi
    file_og_stat_branch="${dir_og_input_root}/stat_branch/${og_id}_stat.branch.tsv"
    file_og_synteny="${dir_og_input_root}/synteny/${og_id}_synteny.tsv"
    file_og_rpsblast="${dir_og_input_root}/rpsblast/${og_id}_rpsblast.tsv"
    file_og_clipkit="${dir_og_input_root}/clipkit/${og_id}_cds.clipkit.fa.gz"
    file_og_orthogroup_extraction_fasta="${dir_og_input_root}/orthogroup_extraction_fasta/${og_id}_orthogroup_extraction.fa.gz"
    file_og_maxalign="${dir_og_input_root}/maxalign/${og_id}_cds.maxalign.fa.gz"
    file_og_mafft="${dir_og_input_root}/mafft/${og_id}_cds.aln.fa.gz"
    file_og_pep_fasta="${dir_og_input_root}/protein_fasta/${og_id}_pep.fa.gz"
    file_og_cds_fasta="${dir_og_input_root}/cds_fasta/${og_id}_cds.fa.gz"
    file_og_dated_tree="${dir_og_input_root}/dated_tree/${og_id}_dated.nwk"
    file_og_fimo="${dir_og_input_root}/fimo/${og_id}_fimo.tsv"
    file_og_meme="${dir_og_input_root}/meme/${og_id}_meme.xml"
    file_hgt_stat_branch="${dir_hgt_tree_input}/${og_id}_hgt_stat.branch.tsv"
    file_hgt_tree_plot="${dir_hgt_tree_plot}/${og_id}_hgt_tree_plot.pdf"
    ortholog_prefix=$(hgt_majority_ortholog_prefix "${og_id}" "${file_hgt_gene}")
    if [[ -z "${ortholog_prefix}" ]]; then
      ortholog_prefix="HGT_UNRESOLVED_"
    fi

    if [[ ! -s "${file_og_stat_branch}" ]]; then
      echo "Skipping HGT tree plot for ${og_id}: stat_branch not found (${file_og_stat_branch})"
      hgt_cleanup_family_materialization "${hgt_current_materialized}"
      hgt_current_materialized=""
      continue
    fi

    python "${gg_support_dir}/annotate_hgt_tree_plot.py" \
      --stat_branch "${file_og_stat_branch}" \
      --branch_tsv "${file_hgt_branch}" \
      --gene_tsv "${file_hgt_gene}" \
      --orthogroup "${og_id}" \
      --outfile "${file_hgt_stat_branch}"

    num_tip_treeplot=$(hgt_count_tree_tips "${file_hgt_stat_branch}")
    panel_trimmed_aln=$(hgt_select_alignment_input "${num_tip_treeplot}" \
      "${file_og_clipkit}" \
      "${file_og_orthogroup_extraction_fasta}" \
      "${file_og_maxalign}" \
      "${file_og_mafft}" \
      "${file_og_pep_fasta}" \
      "${file_og_cds_fasta}")
    panel_untrimmed_aln=$(hgt_select_alignment_input "${num_tip_treeplot}" \
      "${file_og_orthogroup_extraction_fasta}" \
      "${file_og_mafft}" \
      "${file_og_pep_fasta}" \
      "${file_og_cds_fasta}")

    (
      cd "${dir_hgt_tmp}"
      rm -f -- stat_branch2tree_plot.pdf
      Rscript "${gg_support_dir}/stat_branch2tree_plot.r" \
        --stat_branch="${file_hgt_stat_branch}" \
        --max_delta_intron_present="-0.5" \
        --width="${hgt_tree_plot_width}" \
        --rel_widths="tree,2.2,heatmap,0.65,pointplot,0.55,cluster_membership,0.55,synteny,0.9,tiplabel,0.45,categorical,0.8,signal_peptide,0.12,tm,0.12,intron,0.12,domain,1.0,alignment,1.2,meme,0.95,ortholog,0.8" \
        --panel1="tree,bl_rooted,support_unrooted,species,L" \
        --panel2="heatmap,no,abs,_,expression_,Expression" \
        --panel3="pointplot,no,rel,_,expression_" \
        --panel4="heatmap,no,abs,_,hgt_,HGT evidence" \
        --panel5="cluster_membership,100000" \
        --panel6="synteny,${file_og_synteny},5" \
        --panel7="tiplabel" \
        --panel8="categorical,besthit_lca_rank_display,Hit LCA,-" \
        --panel9="signal_peptide" \
        --panel10="transmembrane_domain" \
        --panel11="intron_number" \
        --panel12="domain,${file_og_rpsblast}" \
        --panel13="alignment,${panel_trimmed_aln},${panel_untrimmed_aln}" \
        --panel14="fimo,${hgt_promoter_bp},${hgt_fimo_qvalue}" \
        --panel15="meme,${file_og_meme}" \
        --panel16="ortholog,${ortholog_prefix},${file_og_dated_tree}" \
        --show_branch_id="yes" \
        --event_method="generax" \
        --species_color_table="PLACEHOLDER" \
        --pie_chart_value_transformation="identity" \
        --long_branch_display="auto" \
        --long_branch_ref_quantile="0.95" \
        --long_branch_detect_ratio="5" \
        --long_branch_cap_ratio="2.5" \
        --long_branch_tail_shrink="0.02" \
        --long_branch_max_fraction="0.1"
      if [[ -s "stat_branch2tree_plot.pdf" ]]; then
        mv_out "stat_branch2tree_plot.pdf" "${file_hgt_tree_plot}"
      else
        echo "Warning: HGT tree plot was not generated for ${og_id}."
      fi
    )
    hgt_cleanup_family_materialization "${hgt_current_materialized}"
    hgt_current_materialized=""
  done
fi

hgt_cleanup_materialization_run
trap - EXIT
echo "$(date): Exiting Singularity environment"
