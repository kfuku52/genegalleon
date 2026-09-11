#!/usr/bin/env bash

# Scheduler header notes:
# - Keep sections ordered as SLURM -> UGE -> PBS across entrypoints.
# - Update job name, CPU count, memory, walltime, log paths, and array size together.
# - UGE resources default to SHIROKANE AGE; SLURM/PBS site-specific lines remain examples.

# SLURM
# Common parameters: job name, cores per task, memory per core, walltime, log files, and working directory.
#SBATCH -J gg_gene_evolution
#SBATCH -c 4
#SBATCH --mem-per-cpu=8G
#SBATCH -t 7-00:00:00
#SBATCH --output=gg_gene_evolution_entrypoint.sh_%A_%a.out
#SBATCH --error=gg_gene_evolution_entrypoint.sh_%A_%a.err
#SBATCH --chdir=.
#SBATCH --ignore-pbs
# Array example for array-aware entrypoints.
#SBATCH -a 1
# Site-specific partition example.
#SBATCH -p epyc,rome,medium
# Optional notifications and single-node examples.
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>

## UGE
# SHIROKANE AGE defaults: shell, working directory, slot count, memory per slot, and ljob.
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -l s_vmem=8G
#$ -l ljob
# Array example for array-aware entrypoints.
#$ -t 1

## PBS
# Common parameters: shell, CPU count, total memory, and exported environment.
#PBS -S /bin/bash
#PBS -l ncpus=4
#PBS -l mem=32G
# Array example for array-aware entrypoints.
#PBS -J 1
# Site-specific queue example.
#PBS -q small
#PBS -V

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# N = Number of rows (excluding the header) in workspace/output/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv for mode_gene_evolution=orthogroup
# N = Number of files in workspace/input/query_gene for mode_gene_evolution=query2family

set -euo pipefail

echo "$(date): Starting"

# Resolve workflow paths for local and scheduler-spooled execution.
gg_bootstrap_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
gg_bootstrap_checked_bases=""
if [[ -n "${KFAUTO_RUNTIME_RELEASE_ROOT:-}" ]]; then
  gg_bootstrap_bases=("${KFAUTO_RUNTIME_RELEASE_ROOT}")
else
  gg_bootstrap_bases=(
    "${SLURM_SUBMIT_DIR:-}"
    "${PBS_O_WORKDIR:-}"
    "${PWD:-}"
    "${gg_bootstrap_script_dir}"
  )
fi
for gg_bootstrap_base in "${gg_bootstrap_bases[@]}"
do
  [[ -n "${gg_bootstrap_base}" ]] || continue
  case ":${gg_bootstrap_checked_bases}:" in
    *":${gg_bootstrap_base}:"*) continue ;;
  esac
  gg_bootstrap_checked_bases="${gg_bootstrap_checked_bases:+${gg_bootstrap_checked_bases}:}${gg_bootstrap_base}"
  for bootstrap_path in \
    "${gg_bootstrap_base}/support/gg_entrypoint_bootstrap.sh" \
    "${gg_bootstrap_base}/workflow/support/gg_entrypoint_bootstrap.sh"
  do
    if [[ -s "${bootstrap_path}" ]]; then
      # shellcheck disable=SC1090
      source "${bootstrap_path}"
      break
    fi
  done
  if declare -F gg_entrypoint_initialize >/dev/null 2>&1; then
    break
  fi
done
unset gg_bootstrap_base gg_bootstrap_bases gg_bootstrap_checked_bases gg_bootstrap_script_dir bootstrap_path
if ! declare -F gg_entrypoint_initialize >/dev/null 2>&1; then
  echo "Failed to locate gg_entrypoint_bootstrap.sh from BASH_SOURCE[0]=${BASH_SOURCE[0]}" >&2
  exit 1
fi
if ! gg_entrypoint_initialize "${BASH_SOURCE[0]}" 1 "gg_gene_evolution"; then
  exit 1
fi
gg_entrypoint_name="gg_gene_evolution_entrypoint.sh"

### Start: Modify this block to tailor your analysis ###

# Mode
mode_gene_evolution="${mode_gene_evolution:-query2family}" # query2family|orthogroup; query2family starts from query gene IDs or FASTA and retrieves homolog families, while orthogroup analyzes existing OrthoFinder orthogroups directly.
gene_evolution_profile="${gene_evolution_profile:-default}" # default|hgt; hgt switches to orthogroup mode and enables HGT-oriented defaults for UniProt annotation, GeneRax DTL reconciliation, GFF/expression/intron evidence, and tree plotting.
input_sequence_mode="${input_sequence_mode:-${GG_COMMON_INPUT_SEQUENCE_MODE:-cds}}" # {cds,protein}; protein mode is partial and deactivates CDS-only analyses.
gene_family_output_storage="${gene_family_output_storage:-${GG_COMMON_GENE_FAMILY_OUTPUT_STORAGE:-zip}}" # zip|files|raw; raw aliases files, while zip queues family artifacts for progress-summary collection.
gene_family_zip_min_batch_files="${gene_family_zip_min_batch_files:-${GG_COMMON_GENE_FAMILY_ZIP_MIN_BATCH_FILES:-100}}" # Deprecated compatibility setting; array tasks queue their own output inventories; collectors control batching.
gene_family_zip_compression="${gene_family_zip_compression:-${GG_COMMON_GENE_FAMILY_ZIP_COMPRESSION:-adaptive}}" # adaptive|deflate|store.
gene_family_zip_compression_level="${gene_family_zip_compression_level:-${GG_COMMON_GENE_FAMILY_ZIP_COMPRESSION_LEVEL:-6}}" # Deflate level 0-9.
gene_family_zip_workers="${gene_family_zip_workers:-${GG_COMMON_GENE_FAMILY_ZIP_WORKERS:-1}}" # Concurrent ZIP writers, capped at 4.
gene_family_final_zip_max_bytes="${gene_family_final_zip_max_bytes:-${GG_COMMON_GENE_FAMILY_FINAL_ZIP_MAX_BYTES:-0}}" # 0 permits one final ZIP of any size; a positive byte limit retains human-readable part ZIPs for larger subdirectories.
gene_family_tmp_retention_days="${gene_family_tmp_retention_days:-${GG_COMMON_GENE_FAMILY_TMP_RETENTION_DAYS:-7}}" # Failed task directories below the active output root's tmp directory are removed after this many days when their family lock is idle; 0 disables the age limit.
gene_family_tmp_max_dirs="${gene_family_tmp_max_dirs:-${GG_COMMON_GENE_FAMILY_TMP_MAX_DIRS:-100}}" # Maximum failed task directories retained below the active output root; oldest excess directories are removed and 0 disables the count limit.
gene_family_tmp_max_bytes="${gene_family_tmp_max_bytes:-${GG_COMMON_GENE_FAMILY_TMP_MAX_BYTES:-107374182400}}" # Maximum bytes retained across inactive failed task directories; 0 disables the byte limit.
gene_family_tmp_max_files="${gene_family_tmp_max_files:-${GG_COMMON_GENE_FAMILY_TMP_MAX_FILES:-100000}}" # Maximum files retained across inactive failed task directories; 0 disables the file-count limit.

# Query2family workflow flags
run_extract_query_fasta=1 # Activated if mode_gene_evolution=query2family. Generate amino acid fasta file for query BLAST.
run_query_blast=1 # Activated if mode_gene_evolution=query2family; search query amino-acid sequences against formatted species CDS/protein databases to retrieve homolog candidates.
run_extract_primary_fasta=1 # Generate in-frame CDS fasta file.
run_rps_blast=1 # RPS-BLAST protein domain search.
run_uniprot_annotation=0 # Annotation against UniProt Swiss-Prot.
run_cdskit_localize=1 # Predict targeting-peptide and peroxisome localization signals with cdskit localize.

# Alignment and tree workflow flags
run_mafft=1 # In-frame nucleotide alignment using MAFFT.
run_alignment_stats_original=1 # Alignment statistics before MaxAlign and TrimAl using cdskit stats.
run_maxalign=0 # Remove anomalous sequences by cdskit maxalign.
run_trimal=0 # Remove less-alignable codon sites.
run_clipkit=1 # Remove less-alignable codon sites.
run_alignment_stats_cleaned=1 # Alignment statistics after MaxAlign and TrimAl using cdskit stats.
run_iqtree=1 # Maximum-likelihood phylogenetic reconstruction.
run_tree_root=1 # Root gene tree using tree_rooting_method.

# Reconciliation and dating workflow flags
run_orthogroup_extraction=0 # Optional query2family refinement; extract from the pre-GeneRax rooted homolog tree, then run GeneRax on the extracted tree and FASTA when enabled.
run_generax=0 # GeneRax off by default for local/smoke environments without MPI setup.
reconciliation_duplication_cost=1.5 # Duplication weight for species-tree-assisted rooting.
reconciliation_loss_cost=1 # Implied-loss weight for species-tree-assisted rooting.
run_reconciliation=0 # Run NWKIT LCA reconciliation for RADTE and event statistics.
run_tree_dating=0 # Species-tree-guided divergence time estimation with NWKIT; native or IQ-TREE sequence engine.

# Trait and promoter workflow flags
run_generate_expression_matrix=0 # Generate trait matrix of gene expression level.
run_collect_gff_info=0 # Generate trait matrix on introns and gene positions.
run_extract_promoter_fasta=0 # Generate a promoter fasta file from reference genomes.
run_fimo=0 # Identify promoter motif sequences.
run_tree_pruning=0 # If 1, discard genes without expression data.

# Workflow control flags
check_pruned=0 # Delete downstream outputs if inconsistent to run_tree_pruning.

# dN/dS workflow flags
run_mapdnds_parameter_estimation=0 # Parameter estimation for mapdNdS.
run_mapdnds=0 # Stochastic substitution mapping to estimate dN/dS by mapdNdS.
run_codeml_two_ratio=0 # Run codeml two-ratio model for dN/dS.
run_hyphy_dnds=0 # ML dN/dS estimation by HyPhy FitMG94.bf.
run_hyphy_relax="${run_hyphy_relax:-0}" # Run HyPhy RELAX.
run_hyphy_relax_reversed="${run_hyphy_relax_reversed:-0}" # Run HyPhy RELAX with reversed foreground/background.

# Comparative-analysis workflow flags
run_asr_intron=0 # NWKIT ancestral-state probabilities of intron presence/absence.
run_native_ou=0 # OU modeling of gene expression using NWKIT; set 1 to enable.
run_expression_trait_pgls=0 # Unified expression ~ species-trait analysis with selectable RSC and species-tree PGLS methods.
run_iqtree_anc=0 # Ancestral state reconstruction required for CSUBST.
run_csubst=0 # Protein convergence analysis with CSUBST.
run_csubst_scan=0 # Direct recurrent amino-acid/state-change scan with CSUBST.
run_summary=1 # Generate summary tables.
run_tree_plot=1 # Tree visualization pdf.

# Query2family parameters
query_blast_method="diamond" # diamond|tblastn; diamond translates species CDS and runs protein BLASTP-style searches, while tblastn searches protein queries against nucleotide CDS databases.
query_blast_evalue="auto" # BLAST E-value threshold. Set "auto" to choose by query max AA length.
query_blast_auto_evalue_maxlen_cutoffs="40:1000,80:100,150:10,300:1,inf:0.01" # Auto E-value table: max_query_aa_length:evalue.
query_blast_coverage="0.25" # BLAST coverage threshold.
max_num_gene_blast_hit_retrieval=5000 # Maximum number of genes to retrieve.
retain_query_in_maxalign=1 # Keep query sequences during MaxAlign filtering in query2family mode so seed genes are not removed as alignment outliers.

# Annotation parameters
uniprot_annotation_method="mmseqs2" # blastp|mmseqs2 for UniProt Swiss-Prot annotation search engine.
cdskit_localize_model="${cdskit_localize_model:-targeting5-perox-deeploc21-et-v1}" # cdskit localize model path or alias; default includes the peroxisome head.
cdskit_localize_organism_group="${cdskit_localize_organism_group:-auto}" # auto|unknown|plant|non_plant; auto infers from GG_COMMON_BUSCO_LINEAGE/busco_lineage.
cdskit_localize_include_features=0 # Include internal cdskit localize sequence features in the output TSV.
cdskit_localize_no_model_download=0 # Set 1 to require the localize model to already exist in the cdskit model cache.

# Phylogeny reconstruction and reconciliation parameters
iqtree_fast_mode_gt=2000 # Sequence-count threshold above which IQ-TREE runs with --fast and disables UFBOOT for large alignments.
tree_rooting_method="${tree_rooting_method:-mad}" # mad|reconciliation|midpoint|md; reconciliation uses the pruned species tree with NWKIT, and md maps to NWKIT method "mv".
generax_model="GTR+G4" # GeneRax substitution model.
generax_rec_model="UndatedDL" # "UndatedDTL" or "UndatedDL"; GeneRax reconciliation model, with DL modeling duplication/loss and DTL also allowing transfer events for HGT-oriented analyses.
radte_max_age=1000 # Upper limit of estimated divergence time in MY.
radte_sequence_engine="${radte_sequence_engine:-native}" # native|iqtree; NWKIT handles dates for both.
radte_iqtree_interface="${radte_iqtree_interface:-auto}" # auto|cli|library; optional externally built library worker.
radte_iqtree_worker="${radte_iqtree_worker:-}" # Optional external worker executable.
radte_iqtree_model="${radte_iqtree_model:-}" # Optional complete IQ-TREE model, e.g. GY+F3X4+R4.
radte_substitution_model="${radte_substitution_model:-auto}" # auto: GY94 for CDS, LG for protein; also ecmk07, ecmrest, hky, gtr.
radte_codon_frequencies="${radte_codon_frequencies:-}" # Empty: GY94 f3x4, ECM published frequencies; f selects alignment codon frequencies.
radte_kappa="${radte_kappa:-}" # Empty estimates HKY/GY94 kappa; a number fixes it.
radte_omega="${radte_omega:-}" # Empty estimates GY94 omega; a number fixes it.
radte_gamma_shape="${radte_gamma_shape:-}" # Empty estimates gamma shape.
radte_gamma_categories="${radte_gamma_categories:-4}" # Number of site-rate gamma categories.
radte_inference="${radte_inference:-auto}" # auto|marginal|joint-map.
radte_likelihood="${radte_likelihood:-auto}" # auto|exact|quadratic.
radte_uncertainty="${radte_uncertainty:-profile}" # profile|laplace|bootstrap|none; conditional on fixed species ages.
radte_interval_level="${radte_interval_level:-0.95}" # Nominal conditional interval level; broad coverage is not established.
radte_rate_sd="${radte_rate_sd:-}" # Empty estimates relaxed-clock SD; 0 selects a strict clock.
radte_maxiter="${radte_maxiter:-1000}" # Maximum optimizer iterations.
radte_seed="${radte_seed:-1}" # Random seed.
radte_species_intervals_tsv="${radte_species_intervals_tsv:-}" # Optional external species intervals for display only; workspace-relative or absolute TSV.


# species_expression data (value in input files)
exp_value_type="log2p1" # Expression scale used in species_expression input tables.
pgls_methods="rsc" # rsc,species-nwkit, or all; selected methods share the same prepared expression/trait inputs.
species_expression_aggregation="sum" # sum|mean|max|all; paralogs are combined within each biological sample on the linear expression scale before species-tree PGLS.
species_paralog_missing="error" # error|ignore for incomplete paralog measurements within a species/sample.
species_paralog_sampling_covariance="" # Optional workspace-relative or absolute TSV of response,gene_name_1,gene_name_2,sampling_covariance (and optional tree_id) for known-SE species PGLS aggregation.

# Reconciled speciation contrast (RSC) PGLS
rsc_responses="all" # all or comma-separated expression response names after removing the replicate suffix (for example root,leaf).
rsc_predictors="all" # Ordinary traits automatically; explicitly name GBIF observed columns with metadata, or comma-separated other columns.
rsc_predictor_mode="separate" # separate|joint; separate fits one predictor at a time and is the recommended screening mode.
rsc_event_source="auto" # auto|nhx|lca|species-overlap; auto uses NHX only when the selected reconciliation tree contains GeneRax D annotations.
rsc_speciation_coverage="complete" # complete|any; complete excludes partially sampled daughter clades.
rsc_event_weighting="event" # event|contrast; event gives each species-tree event equal total weight across repeated paralog contrasts.
rsc_model="hierarchical" # hierarchical|replicate-reml|cluster-hc1.
rsc_gene_branch_length="original" # original|unit; original requires a dated gene tree.
rsc_gene_evolution_model="brownian" # brownian|lambda|ou|kappa|delta|eb|acdc|independent.
rsc_gene_evolution_parameter="auto" # auto or a fixed numeric shape parameter where the selected model has one.
rsc_species_branch_length="original" # original|unit; original requires the dated species tree.
rsc_species_evolution_model="brownian" # Evolutionary covariance model for species-trait contrasts.
rsc_species_evolution_parameter="auto" # auto or a fixed numeric shape parameter.
rsc_inference="wald" # wald|parametric-bootstrap (continuous expression responses).
rsc_bootstrap_replicates=1000 # Simulations used only for parametric-bootstrap inference.
rsc_seed=1 # Reproducible NWKIT simulation seed.
rsc_confidence_level="0.95" # Two-sided confidence interval level.
rsc_reml="yes" # yes|no; REML for Gaussian variance components.
rsc_min_species_events=2 # Families below this number of eligible independent species events are recorded as not_estimable.
rsc_unmatched="error" # error|warn|ignore for unmatched tree/table labels.

# RSC response/predictor replicate and factor metadata
rsc_replicate_separator="_" # Numeric suffix delimiter in wide expression columns, for example root_1 and root_2.
rsc_expression_sample_metadata="" # Optional workspace-relative or absolute TSV mapping expression column,response,biological_id and optional technical_id,batch,SE columns.
rsc_within_variance="pooled" # pooled|leaf|known-se for biological expression replicates.
rsc_technical_aggregation="error" # error|mean for technical expression replicates.
rsc_predictor_biological_id="" # Optional species-trait table column identifying biological predictor replicates.
rsc_predictor_technical_id="" # Optional technical-replicate column nested within predictor biological IDs.
rsc_predictor_batch="" # Optional predictor batch column.
rsc_predictor_within_variance="pooled" # pooled|leaf|known-se for predictor replicates.
rsc_predictor_technical_aggregation="error" # error|mean for technical predictor replicates.
rsc_predictor_standard_error_columns="" # Comma-separated known-SE columns corresponding to rsc_predictors.
rsc_predictor_sample_size_columns="" # Optional comma-separated sample-size columns corresponding to known predictor SEs.
rsc_categorical_predictors="" # Numeric-coded unordered predictor columns; string predictors are detected automatically.
rsc_ordered_predictors="" # TRAIT=LOW|MIDDLE|HIGH entries, comma-separated.
rsc_factor_reference="" # TRAIT=REFERENCE entries for unordered factors, comma-separated.
rsc_factor_coding="treatment" # treatment|sum.
rsc_categorical_replicate_policy="latent" # error|latent; latent propagates uncertainty when categorical predictor replicates disagree within a species.

# RSC hierarchical and delayed-origin diagnostics
rsc_event_random_effect="auto" # auto|yes|no.
rsc_lineage_random_slope="auto" # auto|yes|no.
rsc_lineage_inference="none" # none|likelihood-ratio|parametric-bootstrap.
rsc_lineage_leave_one_out="no" # yes|no.
rsc_categorical_origin_diagnostics="none" # none|stochastic-map; applied only to categorical-predictor analyses.
rsc_origin_map_replicates=200 # Stochastic maps per categorical predictor.
rsc_origin_map_threads=1 # Parallel workers for stochastic mapping.
rsc_origin_min_posterior="0.5" # Minimum transition posterior used for origin leave-one-out.
rsc_origin_leave_one_out="no" # yes|no.
rsc_allow_large_dense="no" # yes attempts dense fits beyond NWKIT's normal memory guard.

# Promoter cis-element analysis
promoter_bp=2000 # Upstream promoter length in bp extracted from reference genomes for FIMO motif scans.
fimo_qvalue="0.05" # False discovery rate threshold for FIMO motif search
jaspar_file="latest" # "latest"/"auto" or explicit JASPAR filename in ${dir_jaspardb}

# Ornstein-Uhlenbeck modeling of gene expression evolution
native_ou_criterion="AICc" # AICc|AIC|BIC|pBIC|bootstrap; path requires AIC or AICc.
native_ou_max_shifts="auto" # Tree- and budget-aware cap; integer overrides remain available.
native_ou_calibration_replicates=199 # Complete searches per test when criterion=bootstrap.
native_ou_calibration_level="0.05" # Sequential plug-in bootstrap test level.
native_ou_bootstrap=0 # Complete selection repetitions for stability frequencies.
native_ou_seed=1 # Calibration random stream.
native_ou_bootstrap_seed=2 # Independent stability-bootstrap random stream.
native_ou_convergence=1 # Shared regimes and nested returns; requires auto/exhaustive/lasso.
native_ou_root_model="OUfixedRoot" # OUfixedRoot|OUrandomRoot.
native_ou_trait_covariance="diagonal" # diagonal|full; full estimates evolutionary covariance across traits, shared across regimes.
native_ou_alpha_model="trait-specific" # trait-specific|shared; use shared for a justified common time scale or sensitivity analysis.
native_ou_estimate_measurement_error="yes" # yes|no; estimate extra observation variance beyond sampling SE.
native_ou_search_strategy="auto" # auto|exhaustive|lasso|native-path; path requires convergence=0.
native_ou_replicate_separator="_" # Final suffix delimiter for replicate groups; empty keeps columns separate.

# CSUBST options
csubst_max_arity=10 # Maximum foreground arity considered by CSUBST.
csubst_exhaustive_until=1 # Exhaustively enumerate foreground combinations up to this arity.
csubst_cutoff_stat="OCNany2spe,2.0|omegaCany2spe,5.0" # CSUBST branch-statistic cutoffs used to retain combinations.
csubst_max_combination=10000 # Maximum number of CSUBST combinations retained after filtering.
csubst_fg_exclude_wg="no" # Exclude whole-genome duplication branches from CSUBST foregrounds.
csubst_fg_stem_only="yes" # Restrict CSUBST foreground candidates to stem branches only.
csubst_resolve_binary_foreground="${csubst_resolve_binary_foreground:-no}" # yes|no; assign distinct lineage IDs to disconnected foreground clades in 0/1 trait columns before CSUBST.
csubst_nonsyn_recode="${csubst_nonsyn_recode:-${GG_COMMON_CSUBST_NONSYN_RECODE:-no}}" # no|3di20|dayhoff6|sr6|kgb6|sr4|dayhoff9|dayhoff12|dayhoff15|dayhoff18|srchisq6|kgbauto6; optional amino-acid recoding scheme used for CSUBST nonsynonymous convergence tests.
csubst_scan_unit_mode="${csubst_scan_unit_mode:-clade}" # lineage|stem|clade definition of independent foreground support units; clade follows the current CSUBST default.
csubst_scan_match="${csubst_scan_match:-any2spe}" # CSUBST scan recurrent substitution pattern classes: any2spe by default, or comma-separated classes/all.
csubst_scan_min_event_pp="${csubst_scan_min_event_pp:-0.5}" # Posterior probability threshold for candidate discovery and support calls.
csubst_scan_min_support="${csubst_scan_min_support:-2}" # Minimum foreground-unit support; "1" means one unit, fractional values such as "0.5" are proportions, and "1.0" means 100%.
csubst_scan_rate_event_mode="${csubst_scan_rate_event_mode:-posterior_sum}" # posterior_sum|called; event mass used for scan rate tests.
csubst_scan_rate_length="${csubst_scan_rate_length:-n_rescaled}" # raw|sn_rescaled|n_rescaled branch-length scale for scan rate tests.
csubst_scan_rate_exposure="${csubst_scan_rate_exposure:-q_weighted}" # q_weighted|state_aware|raw_branch_length exposure model for scan rate tests.
csubst_scan_other_scope="${csubst_scan_other_scope:-all}" # all|sister control branch set for scan rate tests.
csubst_scan_site_plot="${csubst_scan_site_plot:-yes}" # Generate csubst scan tree + detected-site summary plot.
csubst_scan_tree_site_plot_format="${csubst_scan_tree_site_plot_format:-pdf}" # pdf|png|svg format for csubst scan site plot.
csubst_scan_tree_site_plot_max_sites="${csubst_scan_tree_site_plot_max_sites:-30}" # Maximum detected sites shown in the csubst scan site plot.

# Intron and chromosomal character evolution
intron_gain_rate="0.0001" # Fixed absent-to-present rate for intron ASR, per dated-tree branch-length unit.
retrotransposition_rate="0.001" # Fixed present-to-absent rate for intron ASR, per dated-tree branch-length unit.

# Tree visualization
treevis_event_method="species_overlap" # "auto", "generax", or "species_overlap"; source for duplication/transfer/loss event labels in treevis plots.
treevis_clade_ortholog=1 # Prefix clade-ortholog labels in treevis plots with the resolved annotation species when available.
treevis_support_value="auto" # "auto", "support_generax_ufboot", "support_unrooted", "dup_conf_score", "no"; auto prefers GeneRax-topology UFBoot when available.
treevis_branch_length="bl_rooted" # "bl_dated", "bl_rooted", "mapdnds_omega"; branch-length metric used for treevis tree geometry.
treevis_branch_color="ou_native_regime" # "species", "no", or *_regime; branch color source for treevis, including species colors or inferred regime assignments.
treevis_retrotransposition_delta_intron="-0.5" # Delta-intron cutoff used to flag retrotransposition candidates in plots.
treevis_heatmap_transform="no" # "no", "log2", "log10p1", "log2p1"; transform applied to numeric matrix values before rendering treevis heatmaps.
treevis_pie_chart_value_transformation="identity" # identity|delog2|delog2p1|delog10|delog10p1; transform expression-like values before pie-chart rendering.
treevis_max_intergenic_dist=100000 # Maximum distance between genes in bp.
treevis_synteny=1 # Generate and display local synteny evidence around orthogroup genes when genome/GFF inputs are available.
treevis_sequence_similarity=1 # Display pairwise amino-acid identity from the analysis alignment when available.
treevis_sequence_similarity_width_mm=15 # Physical width of the sequence identity triangle.
treevis_cis_similarity=1 # Display motif-set similarity when FIMO and promoter FASTA are available.
treevis_cis_similarity_width_mm=15 # Physical width of the promoter cis-element triangle.
treevis_synteny_similarity=1 # Display pairwise synteny similarity when neighborhood input is available.
treevis_synteny_similarity_width_mm=15 # Physical width of the triangular similarity column.
synteny_search_window=20 # Number of neighboring genes on each side included in the similarity search.
treevis_synteny_window=5 # Number of neighboring genes on each side displayed in the treevis synteny panel.
treevis_query_marker=1 # In query2family mode, display query best-hit markers as a treevis categorical panel.
treevis_long_branch_display="auto" # "auto" or "no"; auto detects unusually long branches and compresses their displayed length using the long-branch thresholds below.
treevis_long_branch_ref_quantile="0.95" # Reference branch-length quantile used for long-branch detection.
treevis_long_branch_detect_ratio=5 # Ratio above the reference used to flag long branches.
treevis_long_branch_cap_ratio="2.5" # Maximum displayed long-branch length as a multiple of the reference.
treevis_long_branch_tail_shrink="0.02" # Fraction of the shrunken long-branch tail kept after capping.
treevis_long_branch_max_fraction="0.1" # Maximum fraction of plot width allocated to a capped long branch.

### End: Modify this block to tailor your analysis ###

# Misc
exit_if_running=0 # Exit without main analysis if the same GG_ARRAY_TASK_ID is already running.
delete_tmp_dir=1 # After this run, delete tmp directory created for each job. Set 0 when debugging.
delete_preexisting_tmp_dir=1 # Before starting this job, delete tmp directory created by previous run.

source "${gg_support_dir}/gg_util.sh" # loading utility functions
# Reject retired SCM overrides instead of silently leaving ASR disabled.
if [[ -n "${run_scm_intron+x}" || -n "${GG_GENE_EVOLUTION_RUN_SCM_INTRON+x}" ]]; then
  echo "run_scm_intron is retired; use run_asr_intron (GG_GENE_EVOLUTION_RUN_ASR_INTRON for environment overrides)." >&2
  exit 1
fi
# Forward config variables (including external overrides) into container environment.
gg_apply_registered_env_overrides "${gg_entrypoint_name}" "delete_tmp_dir" "delete_preexisting_tmp_dir"
forward_config_vars_to_container_env "${gg_entrypoint_name}" "delete_tmp_dir" "delete_preexisting_tmp_dir"
if ! gg_entrypoint_prepare_container_runtime 1; then
  exit 1
fi
gg_entrypoint_activate_container_runtime

gg_entrypoint_enter_workspace
set +e
gg_runtime_core_script="$(gg_prepare_entrypoint_runtime_snapshot "${gg_entrypoint_name}" "${gg_core_dir}/gg_gene_evolution_core.sh")"
gg_run_container_shell_script "${gg_container_image_path}" "${gg_runtime_core_script}"
cmd_exit_code=$?
set -e
if [[ ${cmd_exit_code} -eq 8 ]]; then
  echo "Output files were detected. No more run of gg_gene_evolution_core.sh is necessary."
elif [[ ${cmd_exit_code} -ne 0 ]]; then
  echo "gg_gene_evolution_core.sh failed with code ${cmd_exit_code}."
  exit "${cmd_exit_code}"
fi
gg_require_versions_dump "${gg_entrypoint_name}"
echo "$(date): Ending"
