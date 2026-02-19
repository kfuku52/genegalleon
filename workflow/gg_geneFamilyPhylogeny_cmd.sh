#!/usr/bin/env bash

### Start: Modify this block to tailor your analysis ###

# Mode
# Values should be boolean (0 or 1).
# Only 1 mode should be activated.
mode_orthogroup=0 # Analyze OrthoFinder orthogroups
mode_query2family=1 # Analyze all homologs of input genelist in ${dir_pg}/input/query2family_input

# Workflow
# Values should be boolean (0 or 1).
run_get_query_fasta=1 # Activated if mode_query2family=1. Generate amino acid fasta file for query BLAST.
run_query_blast=1 # Activated if mode_query2family=1.
run_get_fasta=1 # Generate in-frame CDS fasta file.
run_rps_blast=${run_rps_blast:-1} # RPS-BLAST protein domain search. The protein domains will appear in tree_plot.pdf.
run_uniprot_annotation=0 # DIAMOND-based annotation against UniProt Swiss-Prot.
run_mafft=1 # In-frame nucleotide alignment using MAFFT.
run_amas_original=1 # Alignment statistics before MaxAlign and TrimAl using AMAS.
run_maxalign=0 # Remove anomalous sequences by cdskit maxalign. This will change the workflow. Fix to 0 or 1 within each project.
run_trimal=0 # Remove less-alignable codon sites. This will change the workflow. Fix to 0 or 1 within each project. If both run_trimal and run_clipkit are set to 1, run_clipkit is prioritized.
run_clipkit=1 # Remove less-alignable codon sites. This will change the workflow. Fix to 0 or 1 within each project. If both run_trimal and run_clipkit are set to 1, run_clipkit is prioritized.
run_amas_cleaned=1 # Alignment statistics after MaxAlign and TrimAl using AMAS.
run_iqtree=1 # Maximum-likelihood phylogenetic reconstruction.
run_tree_root=${run_tree_root:-0} # Root gene tree using tree_rooting_method.
tree_rooting_method="${tree_rooting_method:-mad}" # notung|midpoint|mad|md; md is mapped to nwkit method "mv".
run_orthogroup_extraction=0 # Activated if mode_query2family=1. Extracting orthogroup containing query genes from the IQ-TREE tree.
run_generax=${run_generax:-1} # Phylogeny reconciliation including gene tree rooting. This step improves tree topology.
run_notung_reconcil=0 # Run NOTUNG for RADTE. Alternatively, set run_generax=1 to use GeneRax reconciliation in RADTE.
run_tree_dating=0 # Species-tree-guided divergence time estimation with RADTE.
run_get_expression_matrix=0 # Generate trait matrix of gene expression level.
run_get_gff_info=0 # Generate trait matrix on introns and gene positions. Need gff input.
run_get_promoter_fasta=0 # Generate a promoter fasta file from reference genomes.
run_fimo=0 # Identify promoter motif sequences.
run_tree_pruning=0 # If set to 1, discard genes without expression data. When you change this option in the middle of analysis, make sure to set check_pruned=1 to delete downstream outputs generated previously. Use it with run_get_expression_matrix=1
check_pruned=0 # Delete downstream outputs if it's inconsistent to the run_tree_pruning setting.
run_mapdnds_parameter_estimation=0 # Parameter estimation for mapdNdS.
run_mapdnds=0 # Stochastic substitution mapping to estimate dN/dS by mapdNdS. https://academic.oup.com/mbe/article/35/3/734/4705835
run_codeml_two_ratio=0 # Run codeml's two ratio model for dN/dS. workspace/species_trait/species_trait.tsv is necessary to specify foreground taxa. Generax reconciliation is a prerequisite.
run_hyphy_dnds=0 # ML dN/dS estimation by HyPhy FitMG94.bf. https://github.com/veg/hyphy-analyses/tree/master/FitMG94
run_hyphy_relax=0 # Run HyPhy RELAX https://stevenweaver.github.io/hyphy-site/methods/selection-methods/
run_hyphy_relax_reversed=0 # Run HyPhy RELAX with reversed foreground and background branches.
run_scm_intron=0 # Stochastic character mapping of intron traits. Need gff3 input.
run_phylogeneticem=0 # OU modeling of gene expression by PhylogeneticEM. Use l1ou if no special reason.
run_l1ou=0 # OU modeling of gene expression by l1ou.
run_pgls_gene_tree=0 # Run PGLS analysis using per-species traits and per-gene expression with gene tree.
run_pgls_species_tree=0 # Run PGLS analysis using per-species traits and per-species expression with species tree.
run_iqtree_anc=0 # Ancestral state reconstruction required for CSUBST
run_csubst=0 # Run protein convergence analysis with CSUBST.
run_summary=1 # Generate summary tables. If any of intermediate files in previous steps are updated, the tables will also be updated automatically.
run_tree_plot=1 # Tree visualization pdf.

# Translation table
genetic_code=1 # Integer. See here https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes

# Sequence selection in mode_query2family=1
query_blast_method="diamond" # diamond|tblastn. Query search backend used when mode_query2family=1. In diamond mode, species CDS are translated to proteins and searched by diamond blastp.
query_blast_evalue="0.01" # BLAST E-value threshold. Used if mode_query2family=1. A very relaxed threshold (=1) will work fine in most cases if run_maxalign=1.
query_blast_coverage="0.25" # BLAST coverage threshold. e.g, if set to 0.5, BLAST hits that cover less than than 50% of query sequences will be removed.
max_num_gene_blast_hit_retrieval="5000" # INTEGER. Maximum number of genes to retrieve from BLAST hits. If the number of BLAST hits is larger than this value, genes with top query coverage (qjointcov) will be sampled.
retain_query_in_maxalign=1 # BOOL. Query genes will not be discarded by MaxAlign. Activated only when mode_query2family=1.

# Phylogeny reconstraction and reconciliation
iqtree_fast_mode_gt="2000" # INTEGER. If the number of genes in a tree is larger than this value, IQ-TREE will run in fast mode. See here http://www.iqtree.org/doc/Command-Reference
generax_model="GTR+G4" # Nucleotide or amino acid substitution model in GeneRax. If set, the same model will be forced in preceeding IQ-TREE tree reconstruction.GTR+G4 for nucleotides and LG+G4 for amino acids are recommended. See here for available models: https://github.com/amkozlov/raxml-ng/wiki/Input-data#evolutionary-model
generax_rec_model="UndatedDL" # "UndatedDTL" or "UndatedDL" See here https://github.com/BenoitMorel/GeneRax/wiki/GeneRax
radte_max_age="1000" # Upper limit of estimated divergence time in million years (MY).

# species_expression data (value in input files)
exp_value_type="log2p1"

# Promoter cis-element analysis
promoter_bp="2000" # Promoter length in bp
fimo_qvalue="0.05" # False discovery rate threshold for FIMO motif search
jaspar_file="latest" # "latest"/"auto" or explicit JASPAR filename in ${dir_jaspardb}

# Ornstein-Uhlenbeck modeling of gene expression evolution
clade_collapse_similarity_method="pearson" # Correlation method to calculate within-clade expression divergence. Clades with extremely similar expression values slow l1ou analysis. Such clades will be collapsed.
clade_collapse_similarity_threshold=0.99 # Correlation threshold for clade collapsing.
l1ou_criterion="AICc" # "pBIC", "mBIC", "BIC", or "AICc". See l1ou documentation.
l1ou_nbootstrap=0 # INTEGER. Number of bootstrap replicates in l1ou regime shift detection. See l1ou documentation.
l1ou_use_fit_file=1 # BOOL. Use intermediate file if any in l1ou.
l1ou_alpha_upper="PhylogeneticEM" # "PhylogeneticEM" calculate alpha upper limit in the way PhylogenenticEM does. This is more robust than the l1ou default, which is sensitive to extremely short branches.
l1ou_convergence=1 # BOOL. Estimate convergent expression regime shifts in l1ou.
large_tree_num_gene=1000 # INTEGER. Define large trees by the number of genes. Evoke fast calculation using this variable.
large_tree_max_nshift=10 # INTEGER. Upper limit of expression shifts in large trees.
phylogeneticem_use_fit_file=1 # BOOL. Use intermediate file if any in PhylogeneticEM.

# CSUBST options
csubst_max_arity=10
csubst_exhaustive_until=1
csubst_cutoff_stat="OCNany2spe,2.0|omegaCany2spe,5.0"
csubst_max_combination=10000
csubst_fg_exclude_wg="no"
csubst_fg_stem_only="yes"

# Substitution model in CSUBST and mapdNdS
if [[ ${genetic_code} -eq 1 ]]; then
  codon_model="ECMK07+F+R4" # Standard genetic code
else
  codon_model="GY+F+R4" # Non-standard genetic code
fi

# Intron and chromosomal character evolution
intron_gain_rate=0.0001 # FLOAT. Assumed intron gain rate in event/MY. For justification , see https://www.biorxiv.org/content/10.1101/409888v2
retrotransposition_rate=0.001 # FLOAT. Assumed retrotransposition rate in event/MY. For justification , see https://www.biorxiv.org/content/10.1101/409888v2

# Tree visualization
treevis_event_method="species_overlap" # Method to annotate branching events. "auto", "generax", or "species_overlap". "generax" is preferred if "auto".
treevis_clade_ortholog=1 # BOOL. Show clade orthologs. Useful to identify Arabidopsis ortholog to browse its function.
treevis_clade_ortholog_prefix="Arabidopsis_thaliana_" # Gene ID prefix for clade orthologs. Suppressed when treevis_clade_ortholog=0
treevis_support_value="support_unrooted" # Branch supports to show. "support_unrooted" shows IQ-TREE supports. "dup_conf_score" for duplication confidence score. "no" to suppress.
treevis_branch_length="bl_rooted" # Branch length. "bl_dated" for divergence time; "bl_rooted" for substitution/site; "mapdnds_omega" for dN/dS estimated by mapdNdS.
treevis_branch_color="l1ou_regime" # Branch colors. "species" for species-wise colors, "no" for black, and "*_regime" column in *.stat.branch.tsv, e.g. "l1ou_regime" for regime-wise colors
treevis_retrotransposition_delta_intron="-0.5" # Maximum delta_intron_present to define retrotransposition. Suppressed if intron scm is not available.
treevis_heatmap_transform="no" # "no", "log2", "log10p1", "log2", "log10p1". Log transformation for a better visualization in expression heatmap. log2p1 means log2(x+1). Don't activate if the expression value is already log-transformed.
treevis_pie_chart_value_transformation="identity" # identity|delog2|delog2p1|delog10|delog10p1
treevis_max_intergenic_dist="100000" # Maximum distance between genes in bp to call the cluster membership.
treevis_long_branch_display="auto" # "auto" compresses extreme long branches only for plotting; "no" disables.
treevis_long_branch_ref_quantile="0.95" # Baseline quantile for long-branch auto detection.
treevis_long_branch_detect_ratio="5" # Branches > (quantile * ratio) are treated as outlier-long.
treevis_long_branch_cap_ratio="2.5" # Display cap ratio applied before tail shrink.
treevis_long_branch_tail_shrink="0.02" # Shrink rate for the over-cap portion of outlier-long branches.
treevis_long_branch_max_fraction="0.1" # Skip compression when too many branches are detected as long.

### End: Modify this block to tailor your analysis ###

### ----------------------------------------------------------------------- ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

dir_pg="/workspace"
dir_myscript="/script/script"
source "${dir_myscript}/gg_util.sh" # Load utility functions
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "base" 1 1

build_iqtree_mem_args() {
  IQTREE_MEM_ARGS=()
  if [[ -n "${MEM_PER_HOST:-}" ]]; then
    IQTREE_MEM_ARGS=(--mem "${MEM_PER_HOST}G")
  fi
}

# Setting modes
if enable_all_run_flags_for_debug_mode "gg debug mode: All run_* variables are forced to set 1, except for too-time-consuming tasks."; then
  #run_orthogroup_extraction=0; echo "gg debug mode: run_orthogroup_extraction=${run_orthogroup_extraction}"
  run_codeml_two_ratio=0; echo "gg debug mode: run_codeml_two_ratio=${run_codeml_two_ratio}"
  csubst_cutoff_stat="OCNany2spe,0|omegaCany2spe,1"; echo "gg debug mode: csubst_cutoff_stat=${csubst_cutoff_stat}"
fi
query_blast_method=$(echo "${query_blast_method}" | tr '[:upper:]' '[:lower:]')
tree_rooting_method=$(echo "${tree_rooting_method}" | tr '[:upper:]' '[:lower:]')
if [[ "${tree_rooting_method}" != "notung" && "${tree_rooting_method}" != "midpoint" && "${tree_rooting_method}" != "mad" && "${tree_rooting_method}" != "md" ]]; then
  echo "Invalid tree_rooting_method: ${tree_rooting_method}"
  echo "tree_rooting_method must be one of notung, midpoint, mad, md. Exiting."
  exit 1
fi
if [[ ${mode_query2family} -eq 1 && ${run_query_blast} -eq 1 ]]; then
  if [[ "${query_blast_method}" != "tblastn" && "${query_blast_method}" != "diamond" ]]; then
	  echo "Invalid query_blast_method: ${query_blast_method}"
	  echo "query_blast_method must be either "tblastn" or "diamond". Exiting."
	  exit 1
  fi
fi
if [[ ${mode_orthogroup} -eq 1 && ${mode_query2family} -eq 1 ]]; then
	echo 'Only one of ${mode_orthogroup} and ${mode_query2family} should be set to 1. Exiting.'
	exit 0
elif [[ ${mode_orthogroup} -eq 0 && ${mode_query2family} -eq 0 ]]; then
	echo 'Either ${mode_orthogroup} or ${mode_query2family} should be set to 1. Exiting.'
	exit 0
elif [[ ${mode_orthogroup} -eq 1 ]]; then
	dir_og="${dir_pg_output}/orthogroup"
	file_genecount="${dir_pg_output}/orthofinder/Orthogroups/Orthogroups.GeneCount.selected.tsv"
	ind=$((${SGE_TASK_ID}-1))
	og_id=$(python -c "import sys,pandas; df=pandas.read_csv(sys.argv[1],sep='\t',header=0); print(df.loc[int(sys.argv[2]),:].iloc[0])" ${file_genecount} ${ind})
	echo "OrthoGroup ID: ${og_id}"
	run_get_query_fasta=0
	run_query_blast=0
  run_orthogroup_extraction=0
elif [[ ${mode_query2family} -eq 1 ]]; then
	dir_og="${dir_pg_output}/query2family"
	dir_genelist="${dir_pg_input}/query2family_input"
	if [[ ! -d "${dir_genelist}" ]]; then
	  echo "Input directory does not exist: ${dir_genelist}"
	  exit 1
	fi
	shopt -s nullglob
	files=( "${dir_genelist}"/* )
	shopt -u nullglob
	if [[ ${#files[@]} -eq 0 ]]; then
	  echo "Input directory is empty: ${dir_genelist}"
	  exit 1
	fi
	idx=$((SGE_TASK_ID-1))
	file_query2family_input="${files[${idx}]}"
	og_id="$(basename "${file_query2family_input}")"
	echo "mode_query2family: ${#files[@]} input genelist file(s) were detected in ${dir_genelist}/"
	echo "mode_query2family: Input genelist file = ${file_query2family_input}"
	echo "output file prefix: ${og_id}"
	if [[ ! -f "${file_query2family_input}" ]]; then
	  echo "Input genelist file not found, probably due to out-of-range SGE_TASK_ID=${SGE_TASK_ID}: ${file_query2family_input}"
	  exit 1
	fi
fi
if [[ -z "$og_id" ]]; then
    echo "og_id is empty. Exiting."
    exit 1
fi

dir_sp_genome="${dir_pg_input}/species_genome"
dir_sp_gff="${dir_pg_input}/species_gff"
dir_sp_expression="${dir_pg_input}/species_expression"
dir_sp_cds="${dir_pg_input}/species_cds"
dir_sp_blastdb="${dir_pg_output}/species_cds_blastdb"
dir_sp_trait="${dir_pg_input}/species_trait"
file_sp_trait="${dir_sp_trait}/species_trait.tsv"
dir_hyphy_batch_files="/opt/conda/envs/base/lib/hyphy/TemplateBatchFiles"
file_og="${dir_pg_output}/orthofinder/Orthogroups/Orthogroups.selected.tsv"
if [[ -s "${dir_pg_output}/species_tree/dated_species_tree.nwk" ]]; then
	file_species_tree_base="dated_species_tree"
elif [[ -s "${dir_pg_output}/species_tree/undated_species_tree.nwk" ]]; then
	file_species_tree_base="undated_species_tree"
else
	file_species_tree_base="tmp_species_tree"
fi
species_tree="${dir_pg_output}/species_tree/${file_species_tree_base}.nwk"
species_tree_generax="${dir_pg_output}/species_tree/${file_species_tree_base}.generax.nwk" # generated later
species_tree_pruned="${dir_pg_output}/species_tree/${file_species_tree_base}.pruned.nwk"
species_list_pruned="${dir_pg_output}/species_tree/${file_species_tree_base}.pruned.txt"
species_map_pruned="${dir_pg_output}/species_tree/${file_species_tree_base}.pruned.smap"
notung_jar="/usr/local/bin/Notung.jar"
dir_rpsblastdb="/usr/local/db/Pfam_LE"
dir_jaspardb="/usr/local/db/jaspar"

# Directory PATHs
# Directories for temporary files
dir_tmp_main="${dir_og}/tmp"
dir_tmp="${dir_tmp_main}/${SGE_TASK_ID}_${og_id}" #_${RANDOM}
# Alignment and gene tree preparation and others
dir_og_query2family_input="${dir_og}/query2family_input"
dir_og_query_aa_fasta="${dir_og}/query_aa.fasta"
dir_og_query_blast="${dir_og}/query_blast"
dir_og_cds_fasta="${dir_og}/cds.fasta"
dir_og_rpsblast="${dir_og}/rpsblast"
dir_og_uniprot_annotation="${dir_og}/uniprot_annotation"
dir_og_mafft="${dir_og}/mafft"
dir_og_maxalign="${dir_og}/maxalign"
dir_og_trimal="${dir_og}/trimal"
dir_og_clipkit="${dir_og}/clipkit"
dir_og_clipkit_log="${dir_og}/clipkit.log"
dir_og_iqtree_tree="${dir_og}/iqtree.tree"
dir_og_orthogroup_extraction_nwk="${dir_og}/orthogroup_extraction_nwk"
dir_og_orthogroup_extraction_fasta="${dir_og}/orthogroup_extraction_fasta"
dir_og_generax_nhx="${dir_og}/generax.tree"
dir_og_generax_nwk="${dir_og}/generax.nwk"
dir_og_generax_xml="${dir_og}/generax.xml"
dir_og_rooted_tree="${dir_og}/rooted_tree"
dir_og_rooted_log="${dir_og}/rooted_tree.log"
dir_og_notung_reconcil="${dir_og}/notung.reconcil"
dir_og_dated_tree="${dir_og}/dated_tree"
dir_og_dated_tree_log="${dir_og}/dated_tree.log"
dir_og_mapdnds_parameter="${dir_og}/mapdNdS.parameter"
dir_og_mapdnds_dn="${dir_og}/mapdNdS.dN.tree"
dir_og_mapdnds_ds="${dir_og}/mapdNdS.dS.tree"
dir_og_codeml_two_ratio="${dir_og}/codeml.two_ratio"
dir_og_hyphy_dnds="${dir_og}/hyphy.dnds"
dir_og_hyphy_relax="${dir_og}/hyphy.relax"
dir_og_hyphy_relax_reversed="${dir_og}/hyphy.relax.reversed"
dir_og_expression="${dir_og}/character.expression"
dir_og_gff_info="${dir_og}/character.gff"
dir_og_scm_intron_summary="${dir_og}/scm.intron.summary"
dir_og_scm_intron_plot="${dir_og}/scm.intron.plot"
# Cis-regulatory motif
dir_og_promoter_fasta="${dir_og}/promoter.fasta"
dir_og_meme="${dir_og}/meme"
dir_og_fimo="${dir_og}/fimo"
dir_og_fimo_collapsed="${dir_og}/fimo.collapsed"
# OU expression modeling
dir_og_pem_rdata="${dir_og}/PhylogeneticEM.rdata"
dir_og_pem_tree="${dir_og}/PhylogeneticEM.tree"
dir_og_pem_regime="${dir_og}/PhylogeneticEM.regime"
dir_og_pem_leaf="${dir_og}/PhylogeneticEM.leaf"
dir_og_pem_plot="${dir_og}/PhylogeneticEM.plot"
dir_og_l1ou_fit_rdata="${dir_og}/l1ou.fit.rdata"
dir_og_l1ou_fit_conv_rdata="${dir_og}/l1ou.fit.conv.rdata"
dir_og_l1ou_fit_tree="${dir_og}/l1ou.fit.tree"
dir_og_l1ou_fit_regime="${dir_og}/l1ou.fit.regime"
dir_og_l1ou_fit_leaf="${dir_og}/l1ou.fit.leaf"
dir_og_l1ou_fit_plot="${dir_og}/l1ou.fit.plot"
# Protein convergence analysis
dir_og_iqtree_anc="${dir_og}/iqtree.anc"
dir_og_csubst_b="${dir_og}/csubst.b"
dir_og_csubst_cb_2="${dir_og}/csubst.cb_2"
dir_og_csubst_cb_stats="${dir_og}/csubst.cb_stats"
# PGLS output
dir_og_gene_pgls="${dir_og}/pgls.gene_tree"
dir_og_gene_pgls_plot="${dir_og}/pgls.gene_tree.plot"
dir_og_species_pgls="${dir_og}/pgls.species_tree"
dir_og_species_pgls_plot="${dir_og}/pgls.species_tree.plot"
# Summary
dir_og_amas_original="${dir_og}/amas.original"
dir_og_amas_cleaned="${dir_og}/amas.cleaned"
dir_og_stat_branch="${dir_og}/stat.branch"
dir_og_stat_tree="${dir_og}/stat.tree"
dir_og_tree_plot="${dir_og}/tree_plot"
dir_og_parameters="${dir_og}/parameters"
# Pruned datasets
dir_og_untrimmed_aln_pruned="${dir_og}/pruned.untrimmed_alignment"
dir_og_trimmed_aln_pruned="${dir_og}/pruned.trimmed_alignment"
dir_og_unrooted_tree_pruned="${dir_og}/pruned.unrooted_tree"
dir_og_rooted_tree_pruned="${dir_og}/pruned.rooted_tree"
dir_og_dated_tree_pruned="${dir_og}/pruned.dated_tree"

# File PATHs
# Alignment and gene tree preparation and others
file_og_query_aa_fasta="${dir_og_query_aa_fasta}/${og_id}_query.aa.fa.gz"
file_og_query_blast="${dir_og_query_blast}/${og_id}_query_blast.tsv"
file_og_cds_fasta="${dir_og_cds_fasta}/${og_id}_cds.fa.gz"
file_og_rpsblast="${dir_og_rpsblast}/${og_id}_rpsblast.tsv"
file_og_uniprot_annotation="${dir_og_uniprot_annotation}/${og_id}_uniprot.tsv"
file_og_mafft="${dir_og_mafft}/${og_id}_cds.aln.fa.gz"
file_og_maxalign="${dir_og_maxalign}/${og_id}_cds.maxalign.fa.gz"
file_og_trimal="${dir_og_trimal}/${og_id}_cds.trimal.fa.gz"
file_og_clipkit="${dir_og_clipkit}/${og_id}_cds.clipkit.fa.gz"
file_og_clipkit_log="${dir_og_clipkit_log}/${og_id}_cds.clipkit.log"
file_og_iqtree_tree="${dir_og_iqtree_tree}/${og_id}_iqtree.nwk"
file_og_orthogroup_extraction_nwk="${dir_og_orthogroup_extraction_nwk}/${og_id}_orthogroup_extraction.nwk"
file_og_orthogroup_extraction_fasta="${dir_og_orthogroup_extraction_fasta}/${og_id}_orthogroup_extraction.fa.gz"
file_og_generax_nhx="${dir_og_generax_nhx}/${og_id}_generax.nhx"
file_og_generax_nwk="${dir_og_generax_nwk}/${og_id}_generax.nwk"
file_og_generax_xml="${dir_og_generax_xml}/${og_id}_generax.xml"
file_og_rooted_tree="${dir_og_rooted_tree}/${og_id}_root.nwk"
file_og_rooted_log="${dir_og_rooted_log}/${og_id}_root.txt"
file_og_notung_reconcil="${dir_og_notung_reconcil}/${og_id}_notung.reconcil.zip"
file_og_dated_tree="${dir_og_dated_tree}/${og_id}_dated.nwk"
file_og_dated_tree_log="${dir_og_dated_tree_log}/${og_id}_dated.log.txt"
file_og_mapdnds_parameter="${dir_og_mapdnds_parameter}/${og_id}_parameter.zip"
file_og_mapdnds_dn="${dir_og_mapdnds_dn}/${og_id}_mapdNdS.dN.nwk"
file_og_mapdnds_ds="${dir_og_mapdnds_ds}/${og_id}_mapdNdS.dS.nwk"
file_og_codeml_two_ratio="${dir_og_codeml_two_ratio}/${og_id}_codeml.two_ratio.tsv"
file_og_hyphy_dnds="${dir_og_hyphy_dnds}/${og_id}_hyphy.dnds.json"
file_og_hyphy_relax="${dir_og_hyphy_relax}/${og_id}_hyphy.relax.json"
file_og_hyphy_relax_reversed="${dir_og_hyphy_relax_reversed}/${og_id}_hyphy.relax.reversed.json"
file_og_expression="${dir_og_expression}/${og_id}_expression.tsv"
file_og_gff_info="${dir_og_gff_info}/${og_id}_gff.tsv"
file_og_scm_intron_summary="${dir_og_scm_intron_summary}/${og_id}_scm.intron.tsv"
file_og_scm_intron_plot="${dir_og_scm_intron_plot}/${og_id}_scm.intron.pdf"
# Cis-regulatory motif
file_og_promoter_fasta="${dir_og_promoter_fasta}/${og_id}_promoter.fa.gz"
file_og_meme="${dir_og_meme}/${og_id}_meme.xml"
file_og_fimo="${dir_og_fimo}/${og_id}_fimo.tsv"
file_og_fimo_collapsed="${dir_og_fimo_collapsed}/${og_id}_fimo.collapsed.tsv"
# OU expression modeling
file_og_pem_rdata="${dir_og_pem_rdata}/${og_id}_PhylogeneticEM.RData"
file_og_pem_tree="${dir_og_pem_tree}/${og_id}_PhylogeneticEM.tree.tsv"
file_og_pem_regime="${dir_og_pem_regime}/${og_id}_PhylogeneticEM.regime.tsv"
file_og_pem_leaf="${dir_og_pem_leaf}/${og_id}_PhylogeneticEM.leaf.tsv"
file_og_pem_plot="${dir_og_pem_plot}/${og_id}_PhylogeneticEM.pdf"
file_og_l1ou_fit_rdata="${dir_og_l1ou_fit_rdata}/${og_id}_l1ou.RData"
file_og_l1ou_fit_conv_rdata="${dir_og_l1ou_fit_conv_rdata}/${og_id}_l1ou.conv.RData"
file_og_l1ou_fit_tree="${dir_og_l1ou_fit_tree}/${og_id}_l1ou.tree.tsv"
file_og_l1ou_fit_regime="${dir_og_l1ou_fit_regime}/${og_id}_l1ou.regime.tsv"
file_og_l1ou_fit_leaf="${dir_og_l1ou_fit_leaf}/${og_id}_l1ou.leaf.tsv"
file_og_l1ou_fit_plot="${dir_og_l1ou_fit_plot}/${og_id}_l1ou.pdf"
# Protein convergence analysis
file_og_iqtree_anc="${dir_og_iqtree_anc}/${og_id}_iqtree.anc.zip"
file_og_csubst_b="${dir_og_csubst_b}/${og_id}_csubst_b.tsv"
file_og_csubst_cb_2="${dir_og_csubst_cb_2}/${og_id}_csubst_cb_2.tsv"
file_og_csubst_cb_stats="${dir_og_csubst_cb_stats}/${og_id}_csubst_cb_stats.tsv"
if [[ ${csubst_max_arity} -gt 2 ]]; then
  for i in $(seq 3 ${csubst_max_arity}); do
    declare dir_og_csubst_cb_${i}="${dir_og}/csubst.cb_${i}"
    my_csubst_dir=dir_og_csubst_cb_${i}
    declare file_og_csubst_cb_${i}="${!my_csubst_dir}/${og_id}.csubst_cb_${i}.tsv"
  done
fi
# PGLS output
file_og_gene_pgls="${dir_og_gene_pgls}/${og_id}_gene_PGLS.tsv"
file_og_gene_pgls_plot="${dir_og_gene_pgls_plot}/${og_id}_gene_PGLS.barplot.pdf"
file_og_species_pgls="${dir_og_species_pgls}/${og_id}_species_PGLS.tsv"
file_og_species_pgls_plot="${dir_og_species_pgls_plot}/${og_id}_species_PGLS.barplot.pdf"
# Summary
file_og_stat_branch="${dir_og_stat_branch}/${og_id}_stat.branch.tsv"
file_og_stat_tree="${dir_og_stat_tree}/${og_id}_stat.tree.tsv"
file_og_amas_original="${dir_og_amas_original}/${og_id}_amas.original.tsv"
file_og_amas_cleaned="${dir_og_amas_cleaned}/${og_id}_amas.cleaned.tsv"
file_og_tree_plot="${dir_og_tree_plot}/${og_id}_tree_plot.pdf"
# Pruned datasets
file_og_untrimmed_aln_pruned="${dir_og_untrimmed_aln_pruned}/${og_id}_cds.untrimmed.pruned.fa.gz"
file_og_trimmed_aln_pruned="${dir_og_trimmed_aln_pruned}/${og_id}_cds.trimmed.pruned.fa.gz"
file_og_unrooted_tree_pruned="${dir_og_unrooted_tree_pruned}/${og_id}_unrooted.pruned.nwk"
file_og_rooted_tree_pruned="${dir_og_rooted_tree_pruned}/${og_id}_rooted.pruned.nwk"
file_og_dated_tree_pruned="${dir_og_dated_tree_pruned}/${og_id}_dated.pruned.nwk"

# Define intermediate files for downstream analysis.
# These variables are updated via helper functions to keep routing logic explicit.
set_analysis_file() {
  local slot=$1
  local path=$2
  case "${slot}" in
    untrimmed_aln) file_og_untrimmed_aln_analysis=${path} ;;
    trimmed_aln) file_og_trimmed_aln_analysis=${path} ;;
    unrooted_tree) file_og_unrooted_tree_analysis=${path} ;;
    rooted_tree) file_og_rooted_tree_analysis=${path} ;;
    dated_tree) file_og_dated_tree_analysis=${path} ;;
    *)
      echo "Error: Unknown analysis file slot: ${slot}"
      exit 1
      ;;
  esac
}

set_default_analysis_files() {
  set_analysis_file untrimmed_aln "${file_og_mafft}"
  set_analysis_file trimmed_aln "${file_og_mafft}"
  set_analysis_file unrooted_tree "${file_og_iqtree_tree}"
  set_analysis_file rooted_tree "${file_og_rooted_tree}"
  set_analysis_file dated_tree "${file_og_dated_tree}"
  if [[ ${run_generax} -eq 1 ]]; then
    if [[ -s "${species_tree_pruned}" ]]; then
      set_analysis_file rooted_tree "${file_og_generax_nhx}"
    else
      echo "run_generax is deactivated: missing species tree for GeneRax (${species_tree_pruned})"
      run_generax=0
    fi
  fi
}

switch_alignment_analysis_source() {
  local infile=$1
  set_analysis_file untrimmed_aln "${infile}"
  set_analysis_file trimmed_aln "${infile}"
}

set_default_analysis_files

memory_notung=$((MEM_PER_HOST/4))

echo "Checking parameter conflicts..."
if [[ ${run_trimal} -eq 1 && ${run_clipkit} -eq 1 ]]; then
  echo 'Both of ${run_trimal} and ${run_clipkit} are set to 1.'
	echo '${run_trimal} is deactivated. ${run_clipkit} is still active.'
	run_trimal=0
fi
if [[ ! -d "${dir_sp_gff}" ]] || [[ -z "$(find "${dir_sp_gff}" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
  if [[ ${run_get_gff_info} -eq 1 ]]; then
    echo '${run_get_gff_info} is deactivated. Empty input:' ${dir_sp_gff}
    run_get_gff_info=0
  fi
  if [[ ${run_scm_intron} -eq 1 ]]; then
    echo '${run_scm_intron} is deactivated. Empty input:' ${dir_sp_gff}
    run_scm_intron=0
  fi
fi
if [[ -d "${dir_sp_expression}" ]] && [[ -n "$(find "${dir_sp_expression}" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
  echo '${dir_sp_expression} is not empty. Continued:' ${dir_sp_expression}
else
  echo '${dir_sp_expression} is empty:' ${dir_sp_expression}
  echo '${run_get_expression_matrix}, ${run_tree_pruning}, and other options are deactivated.'
  run_tree_pruning=0
  run_get_expression_matrix=0
  run_phylogeneticem=0
  run_l1ou=0
fi
if [[ -d "${dir_sp_genome}" ]] && [[ -n "$(find "${dir_sp_genome}" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
  echo '${dir_sp_genome} is not empty. Continued:' ${dir_sp_genome}
else
  echo '${dir_sp_genome} is empty:' ${dir_sp_genome}
  echo '${run_get_promoter_fasta} and ${run_meme} are deactivated.'
  run_get_promoter_fasta=0
  run_meme=0
fi

echo "Checking preexisting tmp directory."
if [[ -e ${dir_tmp} && ${delete_preexisting_tmp_dir} -eq 1 ]]; then
	echo "$(date): Deleting preexisting ${dir_tmp}"
	rm -r ${dir_tmp_main}/${SGE_TASK_ID}_*
fi
if [[ ! -e ${dir_tmp} ]]; then
	echo "Making ${dir_tmp}"
	mkdir -p ${dir_tmp}
fi
cd ${dir_tmp}
echo "Working at: $(pwd)"

task="Species tree formatting check"
if [[ ! -s ${species_tree_pruned} || ! -s ${species_map_pruned} ]]; then
  echo "$(date): Warning: ${task}: formatted files were not found."
  echo "Missing: ${species_tree_pruned}"
  echo "Missing: ${species_map_pruned}"
  echo "Run gg_speciesTree to generate ${file_species_tree_base}.pruned.nwk/.pruned.smap before tree-dependent steps."
else
  gg_step_skip "${task}"
fi

task="Query fasta generation"
if [[ ! -s ${file_og_query_aa_fasta} && ${run_get_query_fasta} -eq 1 ]]; then
    gg_step_start "${task}"
    if [[ "$(head --bytes 1 ${file_query2family_input})" == ">" ]]; then
        seqtype=$(seqkit stats --tabular ${file_query2family_input} | awk 'NR>1 {print $3}')
        if [[ ${seqtype} == "DNA" ]]; then
            echo "DNA sequences were detected. The file will be treated as in-frame CDS sequences, translated into amino acids, and used as a ${query_blast_method} query: ${file_query2family_input}"
            seqkit translate --allow-unknown-codon --transl-table ${genetic_code} --threads ${NSLOTS} ${file_query2family_input} > "${og_id}.query.aa.tmp.fasta"
            seqkit seq --threads "${NSLOTS}" "${og_id}.query.aa.tmp.fasta" --out-file "${og_id}.query.aa.out.fa.gz"
            mv_out "${og_id}.query.aa.out.fa.gz" "${file_og_query_aa_fasta}"
            rm -f "${og_id}.query.aa.tmp.fasta"
        elif [[ ${seqtype} == "Protein" ]]; then
            echo "Amino acid sequences were detected. The file will be used as a ${query_blast_method} query: ${file_query2family_input}"
            seqkit seq --threads "${NSLOTS}" "${file_query2family_input}" --out-file "${og_id}.query.aa.out.fa.gz"
            mv_out "${og_id}.query.aa.out.fa.gz" "${file_og_query_aa_fasta}"
        else
            echo "Unsupported sequence type '${seqtype}' in '${file_query2family_input}'. Only "DNA" or "Protein" are allowed. Exiting."
            exit 1
        fi
    else
        echo "Gene IDs were detected. Extracting in-frame CDS sequences from species_cds: ${file_query2family_input}"
        cp_out ${file_query2family_input} ${dir_og_query2family_input}/$(basename "${file_query2family_input}")
        mapfile -t genes < <(sed -e '/^[[:space:]]*$/d' "${file_query2family_input}")
        mapfile -t cds_files < <(gg_find_fasta_files "${dir_sp_cds}" 1)
        if [[ -e pattern.txt ]]; then
            rm pattern.txt
        fi
        touch pattern.txt
        for gene in "${genes[@]}"; do
            echo "${gene}" >> pattern.txt
            if [[ "${gene}" == *"−"* ]]; then
                echo "Query sequence name contains minus sign. Searching the sequence name with hyphen as well: ${gene}"
                echo "${gene//−/-}" >> pattern.txt # Replace minus signs ("−") with hyphens ("-") and add to pattern.txt
            fi
        done
        if [[ -e ${og_id}.query.cds.fasta ]]; then
            rm ${og_id}.query.cds.fasta
        fi
        if [[ -e ${og_id}.query.cds.2.fasta ]]; then
            rm ${og_id}.query.cds.2.fasta
        fi
        touch ${og_id}.query.cds.fasta
        query_hits_tmp_dir="./tmp.query_hits"
        if [[ -d ${query_hits_tmp_dir} ]]; then
            rm -r ${query_hits_tmp_dir}
        fi
        mkdir -p ${query_hits_tmp_dir}
        for file_cds in "${cds_files[@]}"; do
            wait_until_jobn_le ${NSLOTS}
            (
                sp_ub=$(gg_species_name_from_path "${file_cds}")
                query_hits_tmp_file="${query_hits_tmp_dir}/$(basename "${file_cds}").hits.fasta"
                seqkit grep --threads ${NSLOTS} --ignore-case --pattern-file <(awk -v sp="${sp_ub}" '{print $0; print sp "_" $0}' pattern.txt) ${file_cds} \
                | sed -e "s/^>${sp_ub}_/>/" -e "s/^>${sp_ub}-/>/" -e "s/^>${sp_ub}[[:space:]]/>/" -e "s/^>${sp_ub}\./>/" -e "s/^>/>${sp_ub}_/" \
                > ${query_hits_tmp_file}
            ) &
        done
        wait
        for query_hits_tmp_file in ${query_hits_tmp_dir}/*.hits.fasta; do
            if [[ -e "${query_hits_tmp_file}" ]]; then
                cat "${query_hits_tmp_file}" >> ${og_id}.query.cds.fasta
            fi
        done
        rm -r ${query_hits_tmp_dir}
        gg_prepare_cds_fasta_stream "${NSLOTS}" "${genetic_code}" < "${og_id}.query.cds.fasta" \
        | sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
        > ${og_id}.query.cds.2.fasta
        num_query=${#genes[@]}
        num_result=$(grep -c -e "^>" "${og_id}.query.cds.2.fasta")
        echo "Number of gene names in query: ${num_query}"
        echo "Number of gene names in extracted fasta: ${num_result}"
        if [[ ${num_query} -ne ${num_result} ]]; then
            echo "Some gene names were not found in species_cds."
            for gene_name in "${genes[@]}"; do
                if ! grep -q -e "^>${gene_name}" "${og_id}.query.cds.2.fasta"; then
                    echo "Query gene not found in species_cds: ${gene_name}"
                fi
            done
            echo "Exiting."
            exit 1
        fi
        if [[ -s ${og_id}.query.cds.2.fasta ]]; then
            echo "Translating in-frame CDS sequences to amino acid sequences: ${og_id}.query.cds.2.fasta"
            seqkit translate --allow-unknown-codon --transl-table ${genetic_code} --threads ${NSLOTS} ${og_id}.query.cds.2.fasta > "${og_id}.query.aa.tmp.fasta"
            seqkit seq --threads "${NSLOTS}" "${og_id}.query.aa.tmp.fasta" --out-file "${og_id}.query.aa.out.fa.gz"
            mv_out "${og_id}.query.aa.out.fa.gz" "${file_og_query_aa_fasta}"
            rm -f "${og_id}.query.aa.tmp.fasta"
            rm ${og_id}.query.cds.2.fasta
        fi
    fi
else
    gg_step_skip "${task}"
fi

task="In-frame query BLAST (${query_blast_method})"
if [[ ! -s ${file_og_query_blast} && ${run_query_blast} -eq 1 && ${mode_query2family} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ ${query_blast_method} == "tblastn" ]]; then
    if ! type makeblastdb >/dev/null 2>&1; then
      echo "makeblastdb was not found but query_blast_method=tblastn. Exiting."
      exit 1
    fi
    if ! type tblastn >/dev/null 2>&1; then
      echo "tblastn was not found but query_blast_method=tblastn. Exiting."
      exit 1
    fi
  elif [[ ${query_blast_method} == "diamond" ]]; then
    if ! type diamond >/dev/null 2>&1; then
      echo "diamond was not found but query_blast_method=diamond. Exiting."
      exit 1
    fi
    echo "DIAMOND mode selected: species CDS will be translated to proteins because diamond makedb/blastp use protein reference databases."
  fi

  export BLASTDB_LMDB_MAP_SIZE=100000000
  check_species_cds ${dir_pg}
  check_if_species_files_unique ${dir_sp_cds}

  if [[ -e ${og_id}.blastQuery.fasta ]]; then
    rm ${og_id}.blastQuery.fasta
  fi
  touch ${og_id}.blastQuery.fasta

  db_files=()
  ensure_dir "${dir_sp_blastdb}"
  cds_files=()
  mapfile -t cds_files < <(gg_find_fasta_files "${dir_sp_cds}" 1)
  cds_spp=()
  for cds_file in "${cds_files[@]}"; do
    cds_spp+=( "$(gg_species_name_from_path "${cds_file}")" )
  done
  mapfile -t cds_spp < <(printf "%s\n" "${cds_spp[@]}" | sort -u)
  for sp in "${cds_spp[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    echo "sp: ${sp}"
    sp_cds_candidates=()
    for cds_file in "${cds_files[@]}"; do
      if [[ "$(gg_species_name_from_path "${cds_file}")" == "${sp}" ]]; then
        sp_cds_candidates+=( "${cds_file}" )
      fi
    done
    if [[ ${#sp_cds_candidates[@]} -eq 0 ]]; then
      echo "No CDS file was found for species: ${sp}. Skipping."
      continue
    fi
    mapfile -t sp_cds_candidates < <(printf "%s\n" "${sp_cds_candidates[@]}" | sort)
    sp_cds=${sp_cds_candidates[0]}
    sp_cds_blastdb="${dir_sp_blastdb}/$(basename ${sp_cds})"
    db_files+=("${sp_cds_blastdb}")
    if [[ ${query_blast_method} == "tblastn" ]]; then
      echo "makeblastdb input CDS file: ${sp_cds}"
      echo "makeblastdb output database file: ${sp_cds_blastdb}"
      if [[ ! -e ${sp_cds_blastdb}.nhr || ! -e ${sp_cds_blastdb}.nin || ! -e ${sp_cds_blastdb}.nsq || ! -e ${sp_cds_blastdb}.ndb ]]; then
        db_lock_file="${sp_cds_blastdb}.tblastn.lock"
        (
          flock 9
          if [[ ! -e ${sp_cds_blastdb}.nhr || ! -e ${sp_cds_blastdb}.nin || ! -e ${sp_cds_blastdb}.nsq || ! -e ${sp_cds_blastdb}.ndb ]]; then
            if zgrep -e "^>" "${sp_cds}" | grep -q -e "[[:blank:]]"; then
              echo "Space is detected. Please remove all annotation info after spaces in sequence names. Exiting: ${sp_cds}"
              exit 1
            fi
            if zgrep -e "^>" "${sp_cds}" | grep -q -e "|"; then
              echo "Bar (|) is detected. Bars in sequence names will be replaced with underlines (_): ${sp_cds}"
            fi
            echo "Generating BLAST database: ${sp_cds}" | tee >(cat >&2)
            if [[ ${sp_cds} =~ ".gz" ]]; then
              seqkit seq --threads "${NSLOTS}" "${sp_cds}" | makeblastdb -dbtype nucl -title ${sp_cds} -out ${sp_cds_blastdb}
            else
              makeblastdb -dbtype nucl -in ${sp_cds} -out ${sp_cds_blastdb}
            fi
          fi
        ) 9>"${db_lock_file}"
      fi
    elif [[ ${query_blast_method} == "diamond" ]]; then
      sp_cds_diamond_fasta="${sp_cds_blastdb}.diamond.fasta"
      echo "diamond input CDS file: ${sp_cds}"
      echo "diamond translated protein file: ${sp_cds_diamond_fasta}"
      echo "diamond database file: ${sp_cds_blastdb}.dmnd"
      if [[ ! -e ${sp_cds_blastdb}.dmnd ]]; then
        db_lock_file="${sp_cds_blastdb}.diamond.lock"
        (
          flock 9
          if [[ ! -e ${sp_cds_blastdb}.dmnd ]]; then
            if zgrep -e "^>" "${sp_cds}" | grep -q -e "[[:blank:]]"; then
              echo "Space is detected. Please remove all annotation info after spaces in sequence names. Exiting: ${sp_cds}"
              exit 1
            fi
            if zgrep -e "^>" "${sp_cds}" | grep -q -e "|"; then
              echo "Bar (|) is detected. Bars in sequence names will be replaced with underlines (_): ${sp_cds}"
            fi
            echo "Generating DIAMOND database: ${sp_cds}" | tee >(cat >&2)
            if [[ ${sp_cds} =~ ".gz" ]]; then
              seqkit seq --remove-gaps --threads "${NSLOTS}" "${sp_cds}" \
              | seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${NSLOTS}" \
              | sed '/^>/! s/\*//g' \
              > "${sp_cds_diamond_fasta}"
            else
              seqkit seq --remove-gaps --threads "${NSLOTS}" "${sp_cds}" \
              | seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${NSLOTS}" \
              | sed '/^>/! s/\*//g' \
              > "${sp_cds_diamond_fasta}"
            fi
            if [[ $(head -c 1 "${sp_cds_diamond_fasta}") != '>' ]]; then
              sed -e "1d" "${sp_cds_diamond_fasta}" > "${sp_cds_diamond_fasta}.tmp"
              mv_out "${sp_cds_diamond_fasta}.tmp" "${sp_cds_diamond_fasta}"
            fi
            if [[ ! -s "${sp_cds_diamond_fasta}" ]]; then
              echo "Translated FASTA for DIAMOND is empty: ${sp_cds_diamond_fasta}. Exiting."
              exit 1
            fi
            diamond makedb --in "${sp_cds_diamond_fasta}" --db "${sp_cds_blastdb}"
            if [[ $? -ne 0 ]]; then
              echo "diamond makedb failed for ${sp_cds}. Exiting."
              exit 1
            fi
            rm -f "${sp_cds_diamond_fasta}"
          fi
        ) 9>"${db_lock_file}"
      fi
    fi
  done
  wait
  echo "db_files: ${db_files[*]}"
  query_aa_local="${og_id}.query.aa.tmp.for_blast.fasta"
  seqkit seq --threads "${NSLOTS}" "${file_og_query_aa_fasta}" --out-file "${query_aa_local}"

  outfmt="qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore frames qlen slen"
  if [[ ${query_blast_method} == "tblastn" ]]; then
    db_files_str=$(printf " %s" "${db_files[@]}")
    db_files_str="${db_files_str# }"
    echo "Running tblastn."
    tblastn \
    -query "${query_aa_local}" \
    -db "${db_files_str}" \
    -out blast_out.tsv \
    -db_gencode ${genetic_code} \
    -evalue ${query_blast_evalue} \
    -max_target_seqs 50000 \
    -outfmt "6 ${outfmt}" \
    -num_threads ${NSLOTS}
    if [[ $? -ne 0 ]]; then
      echo "tblastn failed. Exiting."
      exit 1
    fi
  elif [[ ${query_blast_method} == "diamond" ]]; then
    echo "Running diamond blastp."
    rm -f blast_out.tsv
    touch blast_out.tsv
    for db_file in "${db_files[@]}"; do
      if [[ ! -e ${db_file}.dmnd ]]; then
        echo "DIAMOND database file is missing: ${db_file}.dmnd. Exiting."
        exit 1
      fi
      tmp_diamond_out="$(basename ${db_file}).diamond.out.tsv"
      diamond blastp \
      --query "${query_aa_local}" \
      --db "${db_file}" \
      --out "${tmp_diamond_out}" \
      --evalue "${query_blast_evalue}" \
      --max-target-seqs 50000 \
      --threads "${NSLOTS}" \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
      if [[ $? -ne 0 ]]; then
        echo "diamond blastp failed for database: ${db_file}. Exiting."
        exit 1
      fi
      if [[ -s ${tmp_diamond_out} ]]; then
        awk -F '\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,"0/1",$13,$14}' "${tmp_diamond_out}" >> blast_out.tsv
      fi
      rm -f "${tmp_diamond_out}"
    done
  fi
  rm -f "${query_aa_local}"

  python ${dir_myscript}/annotate_blast_coverage.py \
  --in blast_out.tsv \
  --ncpu ${NSLOTS} \
  --outfmt-columns "${outfmt}" \
  --frame-filter "0/1" \
  --out blast_out_inframe.tmp3.tsv

  if [[ -s blast_out_inframe.tmp3.tsv ]]; then
    mv_out blast_out_inframe.tmp3.tsv ${file_og_query_blast}
  else
    echo "No query BLAST hits were detected after in-frame filtering. Exiting."
    exit 1
  fi
else
	gg_step_skip "${task}"
fi

task="Fasta generation"
if [[ ! -s ${file_og_cds_fasta} && ${run_get_fasta} -eq 1 ]]; then
	gg_step_start "${task}"

	if [[ ${mode_orthogroup} -eq 1 ]]; then
    genes=()
	    read -r -a genes <<< "$(awk -v og="${og_id}" '$1==og {$1=""; sub(/^[[:space:]]+/, "", $0); gsub(",", "", $0); gsub(/\t/, " ", $0); sub(/[[:space:]]*$/, "", $0); gsub(/\047|"/, "", $0); print; exit}' "${file_og}")"
	elif [[ ${mode_query2family} -eq 1 ]]; then
	  python ${dir_myscript}/extract_gene_id_from_blast_table.py \
	  --infile ${file_og_query_blast} \
	  --outfile gene_id_list.txt \
	  --min_query_blast_coverage ${query_blast_coverage} \
	  --max_num_gene_blast_hit_retrieval ${max_num_gene_blast_hit_retrieval}
		mapfile -t genes < gene_id_list.txt
	fi
  cds_files=()
  mapfile -t cds_files < <(gg_find_fasta_files "${dir_sp_cds}" 1)
  echo "Number of CDS files in ${dir_sp_cds}: ${#cds_files[@]}"
  spp=()
  for file_cds in "${cds_files[@]}"; do
    sp_ub=$(gg_species_name_from_path "${file_cds}")
    spp+=("${sp_ub}")
  done
	if [[ -e pattern.txt ]]; then
		rm pattern.txt
	fi
	touch pattern.txt
	for gene in "${genes[@]}"; do
    echo "${gene}" >> pattern.txt
	done
	if [[ -e ${og_id}.cds.fasta ]]; then
		rm ${og_id}.cds.fasta
	fi
	touch ${og_id}.cds.fasta
  for file_cds in "${cds_files[@]}"; do
    sp_ub=$(gg_species_name_from_path "${file_cds}")
    seqkit grep --threads ${NSLOTS} --pattern-file pattern.txt "${file_cds}" \
    >> ${og_id}.cds.fasta
  done

  seqkit replace --pattern "X" --replacement "N" --by-seq --ignore-case --threads ${NSLOTS} "${og_id}.cds.fasta" \
  | seqkit replace --pattern " .*" --replacement "" --ignore-case --threads ${NSLOTS} \
  | seqkit replace --pattern "\+" --replacement "_" --ignore-case --threads ${NSLOTS} \
  | cdskit pad --codontable ${genetic_code} \
  | sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
  > ${og_id}.cds.2.fasta

  num_gene=${#genes[@]}
  fasta_genes=()
  mapfile -t fasta_genes < <(awk '/^>/ {sub(/^>/, "", $0); print}' "${og_id}.cds.2.fasta")
  num_seq=${#fasta_genes[@]}
  echo "Number of genes in the orthogroup or BLAST hit: ${num_gene}"
  echo "Number of sequences in the fasta: ${num_seq}"
  if [[ ${num_gene} -eq ${num_seq} ]]; then
      echo "Number of genes and sequences matched. Fasta generation completed!"
      seqkit seq --threads "${NSLOTS}" "${og_id}.cds.2.fasta" --out-file "${og_id}.cds.out.fa.gz"
      mv_out "${og_id}.cds.out.fa.gz" "${file_og_cds_fasta}"
  else
      echo "Number of genes and sequences did not match."
      echo "Genes in the orthogroup or BLAST hit:"
      echo -e "${genes[@]}"
      echo ""
      echo "Genes in the generated FASTA:"
			printf '%s\n' "${fasta_genes[@]}" | sort | tr '\n' ' '
      echo ""
      echo "There may be duplicated or missing sequences."
      echo "If you have recently replaced species_cds files, please make sure to remove species_cds_blastdb before rerunning."
      echo "Exiting."
      exit 1
  fi
else
	gg_step_skip "${task}"
fi

task="Protein RPS-BLAST"
disable_if_no_input_file "run_rps_blast" ${file_og_cds_fasta}
if [[ ! -s ${file_og_rpsblast} && ${run_rps_blast} -eq 1 ]]; then
    gg_step_start "${task}"
    dir_rpsblastdb=$(ensure_pfam_domain_db "${dir_pg}")
    if [[ $? -ne 0 ]]; then
        echo "Failed to prepare Pfam_LE DB. Exiting."
        exit 1
    fi
    if [[ -s "${dir_rpsblastdb}/Pfam.pal" ]]; then
        db_rpsblast="${dir_rpsblastdb}/Pfam"
    else
        rps_db_candidates=()
        mapfile -t rps_db_candidates < <(find "${dir_rpsblastdb}" -maxdepth 1 -type f -name "*.loo" | sort)
        if [[ ${#rps_db_candidates[@]} -eq 0 ]]; then
            echo "No RPS-BLAST DB index (*.loo) was found in: ${dir_rpsblastdb}. Exiting."
            exit 1
        fi
        if [[ ${#rps_db_candidates[@]} -gt 1 ]]; then
            echo "Multiple RPS-BLAST DB indices were found. Using the first sorted entry: ${rps_db_candidates[0]}"
        fi
        db_rpsblast="${rps_db_candidates[0]%.loo}"
    fi
    echo "db_rpsblast: ${db_rpsblast}"

    if [[ -e ${og_id}.rpsblast.tmp.tsv ]]; then
        rm ${og_id}.rpsblast.tmp.tsv
    fi

    seqkit seq --remove-gaps --threads ${NSLOTS} ${file_og_cds_fasta} \
    | seqkit translate --allow-unknown-codon --transl-table ${genetic_code} --threads ${NSLOTS} \
    > ungapped_translated_cds.fas

    if [[ $(head -c 1 ungapped_translated_cds.fas) != '>' ]]; then
        sed -e "1d" ungapped_translated_cds.fas > ungapped_translated_cds2.fas
        mv_out ungapped_translated_cds2.fas ungapped_translated_cds.fas
    fi

    outfmt="qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle"

    rpsblast \
    -query ungapped_translated_cds.fas \
    -db ${db_rpsblast} \
    -out ${og_id}.rpsblast.tmp.tsv \
    -evalue 0.01 \
    -outfmt "6 ${outfmt}" \
    -num_threads ${NSLOTS}

    if [[ $? -ne 0 ]]; then
        echo "RPS-BLAST failed. Exiting."
        exit 1
    fi

    genes=()
    mapfile -t genes < <(awk '/^>/ {sub(/^>/, "", $0); sub(/^[[:space:]]*/, "", $0); sub(/[[:space:]].*$/, "", $0); print}' ungapped_translated_cds.fas)
    for gene in "${genes[@]}"; do
        if ! grep -F -q -- "${gene}" "${og_id}.rpsblast.tmp.tsv"; then
            echo "${gene}: no hit in RPS-BLAST. Appending qlen to output tsv."
            qlen=$(seqkit grep --by-name --pattern "${gene}" ungapped_translated_cds.fas | seqkit fx2tab --length | awk 'NR==1 {print $NF; exit}')
            echo -e "${gene}\t\t\t\t\t\t\t\t\t\t\t\t${qlen}\t\t" >> ${og_id}.rpsblast.tmp.tsv
        else
            echo "${gene}: RPS-BLAST hit found."
        fi
    done
    {
      printf '%s\n' "${outfmt}" | tr ' ' '\t'
      cat "${og_id}.rpsblast.tmp.tsv"
    } > ${og_id}.rpsblast.tsv

    cp_out ${og_id}.rpsblast.tsv ${file_og_rpsblast}
else
	gg_step_skip "${task}"
fi

task="Gene trait extraction from gff files"
disable_if_no_input_file "run_get_gff_info" ${file_og_cds_fasta}
if [[ ! -s ${file_og_gff_info} && ${run_get_gff_info} -eq 1 ]]; then
  gg_step_start "${task}"
  if [[ -e gff2genestat.tsv ]]; then
    rm gff2genestat.tsv
  fi
  seqkit seq --threads "${NSLOTS}" "${file_og_cds_fasta}" --out-file "${og_id}.gff2genestat_input.fasta"

  python ${dir_myscript}/gff2genestat.py \
  --dir_gff ${dir_sp_gff} \
  --feature "CDS" \
  --multiple_hits "longest" \
  --seqfile "${og_id}.gff2genestat_input.fasta" \
  --ncpu ${NSLOTS} \
  --outfile gff2genestat.tsv
  rm -f "${og_id}.gff2genestat_input.fasta"

  if [[ -s gff2genestat.tsv ]]; then
    mv_out gff2genestat.tsv ${file_og_gff_info}
  fi
else
	gg_step_skip "${task}"
fi

task="UniProt annotation by DIAMOND"
disable_if_no_input_file "run_uniprot_annotation" ${file_og_cds_fasta}
if [[ ! -s ${file_og_uniprot_annotation} && ${run_uniprot_annotation} -eq 1 ]]; then
    gg_step_start "${task}"
    uniprot_db_prefix=$(ensure_uniprot_sprot_db "${dir_pg}")
    if [[ $? -ne 0 ]]; then
      echo "Failed to prepare UniProt Swiss-Prot DB. Exiting."
      exit 1
    fi

    seqkit seq --remove-gaps --only-id --threads ${NSLOTS} ${file_og_cds_fasta} \
    | seqkit translate --allow-unknown-codon --transl-table ${genetic_code} --threads ${NSLOTS} \
    > uniprot.query.pep.fas

    diamond blastp \
    --query uniprot.query.pep.fas \
    --threads ${NSLOTS} \
    --db ${uniprot_db_prefix} \
    --very-sensitive \
    --out uniprot.diamond.tsv \
    --outfmt 6 qseqid sseqid pident length evalue bitscore qlen \
    --max-target-seqs 1 \
    --evalue 1e-2

    python ${dir_myscript}/reformat_uniprot_diamond.py \
    --diamond_tsv uniprot.diamond.tsv \
    --query_fasta uniprot.query.pep.fas \
    --uniprot_fasta ${uniprot_db_prefix}.pep \
    --outfile uniprot.annotation.tsv

    cp_out uniprot.annotation.tsv ${file_og_uniprot_annotation}
else
		gg_step_skip "${task}"
fi

task="In-frame mafft alignment"
disable_if_no_input_file "run_mafft" ${file_og_cds_fasta}
if [[ ! -s ${file_og_mafft} && ${run_mafft} -eq 1 ]]; then
  gg_step_start "${task}"

  seqkit seq --threads "${NSLOTS}" "${file_og_cds_fasta}" --out-file tmp.cds.input.fasta
  cdskit mask --seqfile tmp.cds.input.fasta --codontable ${genetic_code} --outfile tmp.cds.fasta

	seqkit translate \
	--allow-unknown-codon \
	--transl-table ${genetic_code} \
	--threads ${NSLOTS} \
	tmp.cds.fasta \
	> tmp.pep.fasta

	mafft \
	--auto \
	--amino \
	--thread ${NSLOTS} \
	--quiet \
	tmp.pep.fasta \
	> tmp.pep.aln.fasta

	cdskit backalign \
	--seqfile tmp.cds.fasta \
	--aa_aln tmp.pep.aln.fasta \
	--codontable ${genetic_code} \
	--outfile ${og_id}.cds.aln.fasta

  seqkit seq --threads "${NSLOTS}" "${og_id}.cds.aln.fasta" --out-file "${og_id}.cds.aln.out.fa.gz"
  mv_out "${og_id}.cds.aln.out.fa.gz" "${file_og_mafft}"
  rm -f tmp.cds.input.fasta
else
	gg_step_skip "${task}"
fi

task="AMAS for original alignment"
disable_if_no_input_file "run_amas_original" ${file_og_untrimmed_aln_analysis}
if [[ ! -s ${file_og_amas_original} && ${run_amas_original} -eq 1 ]]; then
	gg_step_start "${task}"
  seqkit seq --threads "${NSLOTS}" "${file_og_untrimmed_aln_analysis}" --out-file "${og_id}.amas.original.input.fasta"

	AMAS.py summary \
	--in-format fasta \
	--data-type dna \
	--in-files "${og_id}.amas.original.input.fasta"

	mv_out summary.txt ${file_og_amas_original}
  rm -f "${og_id}.amas.original.input.fasta"
else
	gg_step_skip "${task}"
fi

task="MaxAlign"
disable_if_no_input_file "run_maxalign" ${file_og_untrimmed_aln_analysis}
if [[ ! -s ${file_og_maxalign} && ${run_maxalign} -eq 1 ]]; then
	gg_step_start "${task}"

	seqkit seq --threads "${NSLOTS}" "${file_og_untrimmed_aln_analysis}" --out-file "${og_id}.cds.aln.fasta"
	maxalign_keep_regex=""
	if [[ ${mode_query2family} -eq 1 && ${retain_query_in_maxalign} -eq 0 ]]; then
		echo "Query sequence(s) is NOT necessarily retained in MaxAlign."
	elif [[ ${mode_query2family} -eq 1 && ${retain_query_in_maxalign} -eq 1 ]]; then
		echo "Query sequence(s) is retained in MaxAlign."
		maxalign_keep_regex=$(python - "${file_query2family_input}" <<'PY'
import re
import sys

infile = sys.argv[1]
gene_ids = []
with open(infile, 'r', encoding='utf-8', errors='replace') as handle:
    first_char = handle.read(1)
    handle.seek(0)
    if first_char == '>':
        for line in handle:
            if not line.startswith('>'):
                continue
            gene = line[1:].strip().split()[0]
            if gene:
                gene_ids.append(gene)
    else:
        for line in handle:
            gene = line.strip()
            if gene:
                gene_ids.append(gene)

normalized_ids = []
seen = set()
for gene in gene_ids:
    for candidate in (gene, gene.replace('−', '-')):
        if candidate and candidate not in seen:
            seen.add(candidate)
            normalized_ids.append(candidate)

patterns = [f"(?i:{re.escape(gene)}.*)" for gene in normalized_ids]
print(','.join(patterns))
PY
)
		if [[ -z "${maxalign_keep_regex}" ]]; then
			echo "Warning: No query IDs were parsed for MaxAlign --keep. Running without keep constraints."
		fi
	else
		maxalign_keep_regex=""
	fi

	maxalign_cmd=( \
		cdskit maxalign \
		--seqfile "${og_id}.cds.aln.fasta" \
		--outfile "${og_id}.maxalign.output.fasta" \
	)
	if [[ -n "${maxalign_keep_regex}" ]]; then
		maxalign_cmd+=(--keep "${maxalign_keep_regex}")
	fi
	"${maxalign_cmd[@]}"

	echo Number of sequences before MaxAlign: $(gg_count_fasta_records "${og_id}.cds.aln.fasta")
	echo Number of sequences after MaxAlign: $(gg_count_fasta_records "${og_id}.maxalign.output.fasta")

	seqkit seq --threads "${NSLOTS}" "${og_id}.maxalign.output.fasta" --out-file "${og_id}.maxalign.out.fa.gz"
	mv_out "${og_id}.maxalign.out.fa.gz" "${file_og_maxalign}"
	rm -f "${og_id}.maxalign.output.fasta"
else
	gg_step_skip "${task}"
fi
if [[ ${run_maxalign} -eq 1 ]]; then
    switch_alignment_analysis_source "${file_og_maxalign}"
fi

task="TrimAl"
disable_if_no_input_file "run_trimal" ${file_og_untrimmed_aln_analysis}
if [[ ! -s ${file_og_trimal} && ${run_trimal} -eq 1 ]]; then
	gg_step_start "${task}"

	seqkit translate --allow-unknown-codon --transl-table ${genetic_code} --threads ${NSLOTS} ${file_og_untrimmed_aln_analysis} \
	| sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
	> untrimmed.pep.fasta

	seqkit seq --remove-gaps --threads ${NSLOTS} \
	${file_og_untrimmed_aln_analysis} \
	> untrimmed.cds.degap.fasta

	trimal \
	-in untrimmed.pep.fasta \
	-backtrans untrimmed.cds.degap.fasta \
	-out ${og_id}.cds.trimal.tmp1.fasta \
	-ignorestopcodon \
	-automated1

  cdskit rmseq \
  --seqfile ${og_id}.cds.trimal.tmp1.fasta \
  --problematic_percent 100 \
  | cdskit hammer \
  --seqfile "-" \
  --codontable ${genetic_code} \
  --nail 4 \
  --outfile "${og_id}.cds.trimal.tmp2.fasta"

  if [[ -s "${og_id}.cds.trimal.tmp2.fasta" ]]; then
    echo "Copying. Output file detected for the task: ${task}"
    seqkit seq --threads "${NSLOTS}" "${og_id}.cds.trimal.tmp2.fasta" --out-file "${og_id}.cds.trimal.out.fa.gz"
    mv_out "${og_id}.cds.trimal.out.fa.gz" "${file_og_trimal}"
  fi
else
	gg_step_skip "${task}"
fi
if [[ ${run_trimal} -eq 1 ]]; then
    set_analysis_file trimmed_aln "${file_og_trimal}"
fi

task="ClipKIT"
disable_if_no_input_file "run_clipkit" ${file_og_untrimmed_aln_analysis}
if [[ ! -s ${file_og_clipkit} && ${run_clipkit} -eq 1 ]]; then
	gg_step_start "${task}"

  seqkit seq --threads "${NSLOTS}" "${file_og_untrimmed_aln_analysis}" --out-file "${og_id}.cds.clipkit.input.fasta"

	clipkit \
	"${og_id}.cds.clipkit.input.fasta" \
	--mode smart-gap \
	--sequence_type nt \
	--codon \
	--input_file_format "fasta" \
	--output_file_format "fasta" \
	--output "${og_id}.cds.clipkit.tmp.fasta" \
	--log

	cdskit hammer \
  --codontable ${genetic_code} \
  --nail 4 \
  --seqfile "${og_id}.cds.clipkit.tmp.fasta" \
  | cdskit rmseq \
  --problematic_percent 100 \
  --outfile "${og_id}.cds.clipkit.hammer.fasta"

  if [[ -s "${og_id}.cds.clipkit.hammer.fasta" ]]; then
    echo "Copying. Output file detected for the task: ${task}"
    seqkit seq --threads "${NSLOTS}" "${og_id}.cds.clipkit.hammer.fasta" --out-file "${og_id}.cds.clipkit.out.fa.gz"
    mv_out "${og_id}.cds.clipkit.out.fa.gz" "${file_og_clipkit}"
    cp_out "${og_id}.cds.clipkit.tmp.fasta.log" ${file_og_clipkit_log}
  fi
  rm -f "${og_id}.cds.clipkit.input.fasta"
else
	gg_step_skip "${task}"
fi
if [[ ${run_clipkit} -eq 1 ]]; then
    set_analysis_file trimmed_aln "${file_og_clipkit}"
fi

task="AMAS for cleaned alignment"
disable_if_no_input_file "run_amas_cleaned" ${file_og_trimmed_aln_analysis}
if [[ ! -s ${file_og_amas_cleaned} && ${run_amas_cleaned} -eq 1 ]]; then
	gg_step_start "${task}"
  seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.amas.cleaned.input.fasta"

	AMAS.py summary \
	--in-format fasta \
	--data-type dna \
	--in-files "${og_id}.amas.cleaned.input.fasta"

	mv_out summary.txt ${file_og_amas_cleaned}
  rm -f "${og_id}.amas.cleaned.input.fasta"
else
	gg_step_skip "${task}"
fi
if [[ ${run_maxalign} -eq 1 ]]; then
	# This code block should be placed immediately after "AMAS for cleaned alignment".
	# orthogroup_summary.tsv will not include necessary info otherwise.
    num_gene_before_maxalign=$(gg_count_fasta_records "${file_og_mafft}")
    num_gene_after_maxalign=$(gg_count_fasta_records "${file_og_maxalign}")
    echo "Number of genes before MaxAlign: ${num_gene_before_maxalign}"
    echo "Number of genes after MaxAlign: ${num_gene_after_maxalign}"
    if [[ ${num_gene_after_maxalign} -lt 3 ]]; then
        echo "This is not sufficient for tree-based analysis (<3). Exitng."
        exit 1
    fi
fi

task="IQ-TREE"
disable_if_no_input_file "run_iqtree" ${file_og_trimmed_aln_analysis}
if [[ ! -s ${file_og_iqtree_tree} && ${run_iqtree} -eq 1 ]]; then
	gg_step_start "${task}"
	num_seq=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
	if [[ ${num_seq} -ge 4 && ${num_seq} -le ${iqtree_fast_mode_gt} ]]; then
		other_iqtree_params="--ufboot 1000 --bnni"
		file_tree="${og_id}.contree"
	else
		other_iqtree_params=""
		file_tree="${og_id}.treefile"
	fi
	if [[ ${num_seq} -gt ${iqtree_fast_mode_gt} ]]; then
	  other_iqtree_params="${other_iqtree_params} --fast"
	fi

		if [[ ${run_generax} -eq 1 ]]; then
			base_model=${generax_model%%+*}
		aa_models=( Blosum62 cpREV Dayhoff DCMut DEN FLU HIVb HIVw JTT JTT-DCMut LG mtART mtMAM mtREV mtZOA PMB rtREV stmtREV VT WAG LG4M LG4X PROTGTR )
		if printf -- '%s\n' "${aa_models[@]}" | grep -Fxq "${base_model}"; then
      is_aa_model=1
		else
      is_aa_model=0
    fi
		if [[ ${is_aa_model} -eq 1 ]]; then
			echo "Specified substitution model was interpreted as an amino acid model (base model = ${base_model})."
			seqkit translate --allow-unknown-codon --transl-table ${genetic_code} ${file_og_trimmed_aln_analysis} \
			| sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
			> iqtree_input.fa
		else
			echo "Specified substitution model was interpreted as a nucleotide model (base model = ${base_model})."
			seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file iqtree_input.fa
		fi
	else
		seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file iqtree_input.fa
	fi
	echo "IQ-TREE starting..."

	iqtree_model_string=${generax_model}
	build_iqtree_mem_args
	
	iqtree \
	-s iqtree_input.fa \
	-m ${iqtree_model_string} \
	-T AUTO \
	--threads-max ${NSLOTS} \
	--prefix ${og_id} \
	"${IQTREE_MEM_ARGS[@]}" \
	--seed 12345 \
	--redo \
	${other_iqtree_params}

	cp_out ${file_tree} ${file_og_iqtree_tree}
else
	gg_step_skip "${task}"
fi

task="Gene tree rooting"
disable_if_no_input_file "run_tree_root" ${file_og_unrooted_tree_analysis}
if [[ ( ! -s ${file_og_rooted_tree} || ! -s ${file_og_rooted_log} ) && ${run_tree_root} -eq 1 ]]; then
	gg_step_start "${task}"

  if [[ "${tree_rooting_method}" == "notung" ]]; then
    if [[ ! -s ${species_tree_pruned} ]]; then
      echo "tree_rooting_method=notung requires species tree: ${species_tree_pruned}"
      exit 1
    fi
    if [[ -e "./${og_id}.notung.root" ]]; then
      rm -r "./${og_id}.notung.root"
    fi

    echo "memory_notung: ${memory_notung}"
    java -jar -Xmx${memory_notung}g ${notung_jar} \
    -s ${species_tree_pruned} \
    -g ${file_og_unrooted_tree_analysis} \
    --root \
    --infertransfers "false" \
    --treeoutput newick \
    --log \
    --treestats \
    --events \
    --parsable \
    --speciestag prefix \
    --allopt \
    --maxtrees 1000 \
    --nolosses \
    --outputdir "./${og_id}.notung.root"

    rooted_candidates=()
    mapfile -t rooted_candidates < <(find "./${og_id}.notung.root" -maxdepth 1 -type f -name "${og_id}.iqtree.nwk.rooting.*" | sort -V)
    selected_rooted_tree=""
    for candidate in "${rooted_candidates[@]}"; do
      if [[ "${candidate}" =~ \.rooting\.[0-9]+$ ]]; then
        selected_rooted_tree="${candidate}"
        break
      fi
    done
    if [[ -z "${selected_rooted_tree}" ]]; then
      echo "NOTUNG did not generate rooted-tree candidates in ./${og_id}.notung.root"
      exit 1
    fi

    nwkit label --target intnode --force yes --infile "${selected_rooted_tree}" --outfile "${og_id}.root.tmp.nwk"
    mv_out "${og_id}.root.tmp.nwk" ${file_og_rooted_tree}
    {
      echo "tree_rooting_method=notung"
      echo "selected_rooting=${selected_rooted_tree}"
    } > "${og_id}.root.txt"
    mv_out "${og_id}.root.txt" ${file_og_rooted_log}
  else
    nwkit_root_method="${tree_rooting_method}"
    if [[ "${nwkit_root_method}" == "md" ]]; then
      nwkit_root_method="mv"
    fi
    nwkit root --method "${nwkit_root_method}" --infile ${file_og_unrooted_tree_analysis} \
    | nwkit label --target intnode --force yes --outfile "${og_id}.root.tmp.nwk"
    mv_out "${og_id}.root.tmp.nwk" ${file_og_rooted_tree}
    {
      echo "tree_rooting_method=${tree_rooting_method}"
      echo "nwkit_method=${nwkit_root_method}"
    } > "${og_id}.root.txt"
    mv_out "${og_id}.root.txt" ${file_og_rooted_log}
  fi
else
	gg_step_skip "${task}"
fi

task="Orthogroup extraction with NWKIT"
run_orthogroup_extraction_original=${run_orthogroup_extraction} # This variable may be disabled by disable_if_no_input_file but the original value is necessary to properly update file_og_*_analysis
disable_if_no_input_file "run_orthogroup_extraction" ${file_query2family_input} ${file_og_trimmed_aln_analysis} ${file_og_rooted_tree_analysis}
if [[ ( ! -s ${file_og_orthogroup_extraction_nwk} || ! -s ${file_og_orthogroup_extraction_fasta} ) && ${run_orthogroup_extraction} -eq 1 ]]; then
	gg_step_start "${task}"

  if [[ "$(head --bytes 1 ${file_query2family_input})" == ">" ]]; then
    echo "Fasta format was detected. Running run_orthogroup_extraction but gene names in the input fasta may not be compatible with this task."
    comma_separated_genes=$(awk '/^>/ {sub(/^>/, "", $0); sub(/[[:space:]].*$/, "", $0); gsub(/−/, "-", $0); print}' "${file_query2family_input}" | paste -sd, -)
	  else
	    echo "Gene IDs were detected."
	    cp_out ${file_query2family_input} ${dir_og_query2family_input}/$(basename "${file_query2family_input}")
		    comma_separated_genes=$(tr '\n' ',' < "${file_query2family_input}" | sed -e 's/,$//' | tr '−' '-')
	  fi
  echo "Seed genes for orthogroup extraction: ${comma_separated_genes}"

  run_nwkit_subtree () {
    local infile=$1
    echo "Running nwkit subtree for ${infile}"
  	local info_txt=$(nwkit subtree --infile ${infile} --leaves ${comma_separated_genes} --orthogroup "yes" --dup_conf_score_threshold 0 2> /dev/null | nwkit info 2> /dev/null)
	  	local num_leaf
	  	num_leaf=$(awk -F': *' '/Number of leaves/ {print $2; exit}' <<< "${info_txt}")
    echo -e "${num_leaf}\t${infile}" >> tmp_num_leaf.tsv
  }

  subtree_infiles=()
  if [[ "${tree_rooting_method}" == "notung" && -d "./${og_id}.notung.root" ]]; then
    mapfile -t subtree_infiles < <(find "./${og_id}.notung.root" -maxdepth 1 -type f | grep -E "${og_id}.iqtree.nwk.rooting.[0-9]+$" | sort -V)
  fi
  if [[ ${#subtree_infiles[@]} -eq 0 ]]; then
    if [[ -s "${file_og_rooted_tree_analysis}" ]]; then
      subtree_infiles=( "${file_og_rooted_tree_analysis}" )
    else
      echo "No rooted tree is available for orthogroup extraction."
      exit 1
    fi
  fi

  echo -e "num_leaf\tfile" > tmp_num_leaf.tsv
  for subtree_infile in "${subtree_infiles[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    run_nwkit_subtree ${subtree_infile}
  done

  min_leaf_file=$(awk 'NR>1 {print $1, $2}' tmp_num_leaf.tsv | sort -n | awk 'NR==1 {print $2}')
  min_leaf_num=$(awk 'NR>1 {print $1, $2}' tmp_num_leaf.tsv | sort -n | awk 'NR==1 {print $1}')
  max_leaf_file=$(awk 'NR>1 {print $1, $2}' tmp_num_leaf.tsv | sort -nr | awk 'NR==1 {print $2}')
  max_leaf_num=$(awk 'NR>1 {print $1, $2}' tmp_num_leaf.tsv | sort -nr | awk 'NR==1 {print $1}')
  echo "Minimum number of orthogroup subtree leaves after checking all rooting positions: ${min_leaf_num} in ${min_leaf_file} (will be used for orthogroup extraction)"
  echo "Maximum number of orthogroup subtree leaves after checking all rooting positions: ${max_leaf_num} in ${max_leaf_file} (shown just as a reference)"

	  nwkit subtree --infile ${min_leaf_file} --leaves ${comma_separated_genes} --orthogroup "yes" --dup_conf_score_threshold 0 \
	  | nwkit drop --target "intnode" --name "yes" --outfile "${og_id}.orthogroup_extraction.tmp.nwk"
	  mv_out "${og_id}.orthogroup_extraction.tmp.nwk" ${file_og_orthogroup_extraction_nwk}

  seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file tmp.trimmed.input.fasta
  nwkit intersection \
  --infile ${file_og_orthogroup_extraction_nwk} \
  --outfile /dev/null \
  --seqin tmp.trimmed.input.fasta \
	--seqout tmp.fasta \
	--match "complete"
	rm -f tmp.trimmed.input.fasta

	  cdskit hammer --nail 4 -s tmp.fasta -o "${og_id}.orthogroup_extraction.tmp.fasta"
	  seqkit seq --threads "${NSLOTS}" "${og_id}.orthogroup_extraction.tmp.fasta" --out-file "${og_id}.orthogroup_extraction.out.fa.gz"
	  mv_out "${og_id}.orthogroup_extraction.out.fa.gz" "${file_og_orthogroup_extraction_fasta}"
	  rm -f "${og_id}.orthogroup_extraction.tmp.fasta"
else
	gg_step_skip "${task}"
fi
if [[ ${run_orthogroup_extraction_original} -eq 1 ]]; then
  set_analysis_file unrooted_tree "${file_og_orthogroup_extraction_nwk}"
  set_analysis_file trimmed_aln "${file_og_orthogroup_extraction_fasta}"
fi

task="GeneRax"
disable_if_no_input_file "run_generax" ${file_og_trimmed_aln_analysis} ${file_og_unrooted_tree_analysis} ${species_tree_pruned}
if [[ ! -s ${file_og_generax_nhx} && ${run_generax} -eq 1 ]]; then
	gg_step_start "${task}"

	base_model=${generax_model%%+*}
	aa_models=( Blosum62 cpREV Dayhoff DCMut DEN FLU HIVb HIVw JTT JTT-DCMut LG mtART mtMAM mtREV mtZOA PMB rtREV stmtREV VT WAG LG4M LG4X PROTGTR )
	if printf -- '%s\n' "${aa_models[@]}" | grep -Fxq "${base_model}"; then
    is_aa_model=1
	else
    is_aa_model=0
  fi
	if [[ ${is_aa_model} -eq 1 ]]; then
		echo "Specified substitution model was interpreted as an amino acid model (base model = ${base_model})."
		seqkit translate --allow-unknown-codon --transl-table ${genetic_code} ${file_og_trimmed_aln_analysis} \
		| sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
		> generax_input_alignment.fas
	else
		echo "Specified substitution model was interpreted as a nucleotide model (base model = ${base_model})."
		seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file generax_input_alignment.fas
	fi

  nwkit drop --target intnode --support yes --name yes \
	--infile ${file_og_unrooted_tree_analysis} \
	--outfile generax_input_gene_tree.nwk

	#avoid multifurcating tree
	  R -q -e "library(ape); t=read.tree(\"generax_input_gene_tree.nwk\"); t=multi2di(t,random=FALSE); write.tree(t, \"generax_input_gene_tree_bi.nwk\")"

	generate_generax_mapfile () {
		# https://github.com/BenoitMorel/GeneRax/wiki/Gene-to-species-mapping
			my_aln_file=$1
				awk '/^>/ {sub(/^>/, "", $0); print}' "${my_aln_file}" > tmp.gene_names.txt
			while IFS= read -r gene_name; do gg_species_name_from_path "${gene_name}"; done < tmp.gene_names.txt > tmp.species_names.txt
			paste tmp.gene_names.txt tmp.species_names.txt > generax_map.txt
		rm tmp.gene_names.txt tmp.species_names.txt
	}
	generate_generax_mapfile generax_input_alignment.fas

	echo """
	[FAMILIES]
	- family_1
	starting_gene_tree = generax_input_gene_tree_bi.nwk
	alignment = generax_input_alignment.fas
	mapping = generax_map.txt
	subst_model = ${generax_model}
	""" | sed -e "s/^[[:space:]]*//" | grep -v "^$" > generax_families.txt

  #${host_mpiexec_path} -np ${NSLOTS} generax \
  mpiexec_args=(mpiexec -oversubscribe -np ${NSLOTS})
  mpi_env_args=()
  running_under_scheduler=0
  if [[ -n "${SLURM_JOB_ID:-}" || -n "${PBS_JOBID:-}" || -n "${PE_HOSTFILE:-}" || -n "${LSB_JOBID:-}" ]]; then
    running_under_scheduler=1
  fi
  if [[ ${running_under_scheduler} -eq 0 ]]; then
    # Local/container runs may not have ssh/rsh and often probe unavailable OpenIB transports.
    mpi_env_args=(env OMPI_MCA_plm=isolated OMPI_MCA_plm_rsh_agent=/bin/false OMPI_MCA_btl=^openib)
  fi
  if [[ "$(id -u)" -eq 0 ]]; then
    mpiexec_args+=(--allow-run-as-root)
  fi
	"${mpi_env_args[@]}" "${mpiexec_args[@]}" generax \
	--species-tree "${species_tree_pruned}" \
	--families generax_families.txt \
	--strategy "SPR" \
	--rec-model ${generax_rec_model} \
	--prefix "generax_${og_id}" \
	--per-family-rates \
  --skip-family-filtering \
	--mad-rooting \
	--seed 12345

  echo "GeneRax exit code = $?"

	generax_out_sptree="./generax_${og_id}/species_trees/starting_species_tree.newick" # generax v2.0
  if [[ -s ${generax_out_sptree} ]]; then
    lock_file="${species_tree_generax}.lock"
	    if command -v flock >/dev/null 2>&1; then
	      exec 9> "${lock_file}"
	      flock 9
	      if [[ ! -s ${species_tree_generax} ]]; then
	        echo "copying GeneRax output species tree (first writer only)."
	        cp_out ${generax_out_sptree} ${species_tree_generax}
	      fi
	      flock -u 9
	      exec 9>&-
	    else
	      echo "Error: flock command is required but not available."
	      exit 1
	    fi
	  elif [[ ! -s ${species_tree_generax} ]]; then
	    echo "GeneRax species tree file was not found yet: ${generax_out_sptree}"
	  fi
	echo "copying GeneRax output gene tree."
	reconciled_base="./generax_${og_id}/reconciliations/family_1_reconciliated"
	reconciled_xml="${reconciled_base}.xml"
	reconciled_nhx="${reconciled_base}.nhx"
		if [[ -e ${reconciled_nhx} ]]; then
			echo "GeneRax outfile was found. Copying."
			nwkit nhx2nwk --infile ${reconciled_nhx} --outfile "${og_id}.generax.tmp.nwk"
			mv_out "${og_id}.generax.tmp.nwk" ${file_og_generax_nwk}
			cp_out ${reconciled_xml} ${file_og_generax_xml}
			cp_out ${reconciled_nhx} ${file_og_generax_nhx}
	else
		echo "GeneRax outfile was not found. Exiting."
		exit 1
	fi
else
	gg_step_skip "${task}"
fi

task="NOTUNG reconciliation"
disable_if_no_input_file "run_notung_reconcil" ${file_og_rooted_tree} ${species_tree_pruned}
if [[ ! -s ${file_og_notung_reconcil} && ${run_notung_reconcil} -eq 1 ]]; then
	gg_step_start "${task}"

	echo "memory_notung: ${memory_notung}"

	if [[ -s ./${og_id}.root.nwk ]]; then
		rm ${og_id}.root.nwk
	fi
	if [[ -e ./${og_id}.notung.reconcil ]]; then
		rm -r ${og_id}.notung.reconcil
	fi

  nwkit drop --target intnode --support yes --name yes \
	--infile ${file_og_rooted_tree} \
	--outfile ${og_id}.root.nwk

	java -jar -Xmx${memory_notung}g ${notung_jar} \
	-s ${species_tree_pruned} \
	-g ${og_id}.root.nwk \
	--reconcile \
	--infertransfers "false" \
	--treeoutput newick \
	--log \
	--treestats \
	--events \
	--parsable \
	--speciestag prefix \
	--maxtrees 1 \
	--nolosses \
	--outputdir ./${og_id}.notung.reconcil

	if [[ -s ${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled.*parsable.txt || -s ${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled.parsable.txt ]]; then
		zip -rq ${og_id}.notung.reconcil.zip ${og_id}.notung.reconcil
		cp_out ${og_id}.notung.reconcil.zip ${file_og_notung_reconcil}
	fi
else
	gg_step_skip "${task}"
fi

task="Species-tree-guided divergence time estimation"
disable_if_no_input_file "run_tree_dating" ${species_tree_pruned} ${file_og_unrooted_tree_analysis}
if [[ ( ! -s ${file_og_dated_tree} || ! -s ${file_og_dated_tree_log} ) && ${run_tree_dating} -eq 1 ]]; then
	gg_step_start "${task}"

	if [[ ${run_generax} -eq 1 ]]; then
		radte_input="--species_tree=${species_tree_generax}"
		radte_input="${radte_input} --generax_nhx=${file_og_generax_nhx}"
	else
		if [[ -e ./${og_id}.notung.reconcil ]]; then
			rm -r ./${og_id}.notung.reconcil
		fi
		cp_out ${dir_og_notung_reconcil}/${og_id}.notung.reconcil.zip .
		unzip -q ${og_id}.notung.reconcil.zip
		if [[ -s ./${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled.0 ]]; then
			cp_out ./${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled.0 ./${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled
			cp_out ./${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled.0.parsable.txt ./${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled.parsable.txt
		fi
		radte_input="--species_tree=${species_tree_pruned}"
		radte_input="${radte_input} --gene_tree=./${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled"
		radte_input="${radte_input} --notung_parsable=./${og_id}.notung.reconcil/${og_id}.root.nwk.reconciled.parsable.txt"
	fi

	Rscript ${dir_myscript}/radte.r \
	${radte_input} \
	--max_age=${radte_max_age} \
	--chronos_lambda=1 \
	--chronos_model=discrete \
	--pad_short_edge=0.001 \
	2>&1 | tee radte.log

	constrained_node=$(awk -F': *' '/^Calibrated nodes:/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' radte.log)
	echo ${constrained_node} > ${og_id}.dated.log.txt

	if grep -q ":-" radte_gene_tree_output.nwk; then
      contain_negative_bl=1
    else
      contain_negative_bl=0
    fi
	if [[ ${contain_negative_bl} -eq 1 ]]; then
			echo "Dated tree has negative branch length. Deleting output files depending on the tree file."
			for key in l1ou pem scm dated stat tree_plot; do
				files=()
				mapfile -t files < <(compgen -A variable "file_og_${key}")
				for f in "${files[@]}"; do
					target_file="${!f}"
					if [[ -e "${target_file}" ]]; then
						echo "deleting: ${target_file}"
						rm "${target_file}"
					fi
				done
			done
		else
			echo "Dated tree has no negative branch length. Continue."
			cp_out radte_calibrated_nodes.txt ${file_og_dated_tree_log}
			cp_out radte_gene_tree_output.nwk ${file_og_dated_tree}
		fi
else
	gg_step_skip "${task}"
fi

task="Expression matrix preparation"
disable_if_no_input_file "run_get_expression_matrix" ${file_og_trimmed_aln_analysis}
if [[ ! -s ${file_og_expression} && ${run_get_expression_matrix} -eq 1 ]]; then
	gg_step_start "${task}"
  seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.trait_matrix_input.fasta"

	python ${dir_myscript}/get_trait_matrix.py \
	--dir_trait ${dir_sp_expression} \
	--seqfile "${og_id}.trait_matrix_input.fasta" \
	--ncpu ${NSLOTS} \
	--outfile expression_matrix.tsv
  rm -f "${og_id}.trait_matrix_input.fasta"
  if [[ -s expression_matrix.tsv ]]; then
  	mv_out expression_matrix.tsv ${file_og_expression}
  fi
else
	gg_step_skip "${task}"
fi

task="Promoter fasta generation"
disable_if_no_input_file "run_get_promoter_fasta" ${file_og_gff_info}
if [[ ! -s ${file_og_promoter_fasta} && ${run_get_promoter_fasta} -eq 1 ]]; then
	gg_step_start "${task}"

  python ${dir_myscript}/get_promoter_fasta.py \
  --dir_genome ${dir_sp_genome} \
  --geneinfo_tsv ${file_og_gff_info} \
  --seqkit_exe "seqkit" \
  --outfile "${og_id}.promoter.tmp.fasta" \
  --promoter_bp ${promoter_bp} \
  --ncpu ${NSLOTS}
  if [[ -s "${og_id}.promoter.tmp.fasta" ]]; then
    seqkit seq --threads "${NSLOTS}" "${og_id}.promoter.tmp.fasta" --out-file "${og_id}.promoter.out.fa.gz"
    mv_out "${og_id}.promoter.out.fa.gz" "${file_og_promoter_fasta}"
    rm -f "${og_id}.promoter.tmp.fasta"
  fi
else
	gg_step_skip "${task}"
fi

task="fimo"
jaspar_path=""
if [[ ${run_fimo} -eq 1 ]]; then
  jaspar_path=$(ensure_jaspar_file "${dir_pg}" "${jaspar_file}")
  if [[ $? -ne 0 || -z "${jaspar_path}" ]]; then
    echo "Failed to prepare JASPAR motif file (${jaspar_file}). Exiting."
    exit 1
  fi
fi
disable_if_no_input_file "run_fimo" "${file_og_promoter_fasta}" "${jaspar_path}"
if [[ ! -s ${file_og_fimo} && ${run_fimo} -eq 1 ]]; then
  gg_step_start "${task}"
  seqkit seq --threads "${NSLOTS}" "${file_og_promoter_fasta}" --out-file "${og_id}.fimo.input.fasta"

  fimo \
  --oc "fimo_out" \
  ${jaspar_path} \
  "${og_id}.fimo.input.fasta"
  rm -f "${og_id}.fimo.input.fasta"

  if [[ -s "./fimo_out/fimo.tsv" ]]; then
    mv_out "./fimo_out/fimo.tsv" ${file_og_fimo}
    rm -r "./fimo_out"
  fi
else
	gg_step_skip "${task}"
fi

task="Tree pruning"
disable_if_no_input_file "run_tree_pruning" ${file_og_expression} ${file_og_untrimmed_aln_analysis} ${file_og_trimmed_aln_analysis} ${file_og_unrooted_tree_analysis} ${file_og_rooted_tree_analysis}
if [[ ( ! -s ${file_og_untrimmed_aln_pruned} || ! -s ${file_og_trimmed_aln_pruned} || ! -s ${file_og_unrooted_tree_pruned} || ! -s ${file_og_rooted_tree_pruned} ) ]]; then
  is_all_outputs_exist=0
else
  if [[ ${run_tree_dating} -eq 1 && ! -s ${file_og_dated_tree_pruned} ]]; then
    is_all_outputs_exist=0
  else
    is_all_outputs_exist=1
  fi
fi
if [[ ${is_all_outputs_exist} -eq 0 && ${run_tree_pruning} -eq 1 ]]; then
	gg_step_start "${task}"

	cut -f 1 "${file_og_expression}" | tail -n +2 > target_genes.txt

  if [[ -s ${file_og_untrimmed_aln_analysis} ]]; then
    seqkit seq --threads "${NSLOTS}" "${file_og_untrimmed_aln_analysis}" --out-file "${og_id}.untrimmed.input.fasta"
    python -c "import sys,re; keys = open(sys.argv[1]).read().split('\n'); entries = open(sys.argv[2]).read().split('>'); [ sys.stdout.write('>'+e) for e in entries if (re.sub('\n.*','',e) in keys)&(len(e)!=0) ]" \
    target_genes.txt "${og_id}.untrimmed.input.fasta" > "${og_id}.untrimmed.pruned.tmp.fasta"
    rm -f "${og_id}.untrimmed.input.fasta"
    seqkit seq --threads "${NSLOTS}" "${og_id}.untrimmed.pruned.tmp.fasta" --out-file "${og_id}.untrimmed.pruned.out.fa.gz"
    mv_out "${og_id}.untrimmed.pruned.out.fa.gz" "${file_og_untrimmed_aln_pruned}"
    rm -f "${og_id}.untrimmed.pruned.tmp.fasta"
  fi

  if [[ -s ${file_og_trimmed_aln_analysis} ]]; then
    seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.trimmed.input.fasta"
    python -c "import sys,re; keys = open(sys.argv[1]).read().split('\n'); entries = open(sys.argv[2]).read().split('>'); [ sys.stdout.write('>'+e) for e in entries if (re.sub('\n.*','',e) in keys)&(len(e)!=0) ]" \
    target_genes.txt "${og_id}.trimmed.input.fasta" > "${og_id}.trimmed.pruned.tmp.fasta"
    rm -f "${og_id}.trimmed.input.fasta"
    seqkit seq --threads "${NSLOTS}" "${og_id}.trimmed.pruned.tmp.fasta" --out-file "${og_id}.trimmed.pruned.out.fa.gz"
    mv_out "${og_id}.trimmed.pruned.out.fa.gz" "${file_og_trimmed_aln_pruned}"
    rm -f "${og_id}.trimmed.pruned.tmp.fasta"
  fi

  mapfile -t prune_genes < <(sed -e '/^[[:space:]]*$/d' target_genes.txt)
  prune_pattern=""
  if [[ ${#prune_genes[@]} -gt 0 ]]; then
    prune_pattern=$(
      printf '%s\n' "${prune_genes[@]}" \
      | sed -e 's/[][(){}.^$+*?|\\-]/\\&/g' \
      | paste -sd'|' -
    )
  fi

  if [[ -s ${file_og_unrooted_tree_analysis} ]]; then
    if [[ ${#prune_genes[@]} -eq 0 ]]; then
      cat "${file_og_unrooted_tree_analysis}"
    else
      nwkit prune \
      --infile "${file_og_unrooted_tree_analysis}" \
      --pattern "^(${prune_pattern})$" \
      --invert_match yes
    fi \
    | nwkit drop --target root --length yes \
    > "${og_id}.unrooted.pruned.tmp.nwk"
    mv_out "${og_id}.unrooted.pruned.tmp.nwk" ${file_og_unrooted_tree_pruned}
  fi

  if [[ -s ${file_og_rooted_tree_analysis} ]]; then
    if [[ ${#prune_genes[@]} -eq 0 ]]; then
      cat "${file_og_rooted_tree_analysis}"
    else
      nwkit prune \
      --infile "${file_og_rooted_tree_analysis}" \
      --pattern "^(${prune_pattern})$" \
      --invert_match yes
    fi \
    | nwkit drop --target root --length yes \
    > "${og_id}.rooted.pruned.tmp.nwk"
    mv_out "${og_id}.rooted.pruned.tmp.nwk" ${file_og_rooted_tree_pruned}
  fi

  if [[ -s ${file_og_dated_tree_analysis} ]]; then
    if [[ ${#prune_genes[@]} -eq 0 ]]; then
      cat "${file_og_dated_tree_analysis}"
    else
      nwkit prune \
      --infile "${file_og_dated_tree_analysis}" \
      --pattern "^(${prune_pattern})$" \
      --invert_match yes
    fi \
    | nwkit drop --target root --length yes \
    > "${og_id}.dated.pruned.tmp.nwk"
    mv_out "${og_id}.dated.pruned.tmp.nwk" ${file_og_dated_tree_pruned}
  fi
else
	gg_step_skip "${task}"
fi
if [[ ${run_tree_pruning} -eq 1 ]]; then
	num_gene_before_pruning=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
	num_gene_after_pruning=$(gg_count_fasta_records "${file_og_trimmed_aln_pruned}")
    echo "Number of genes before pruning: ${num_gene_before_pruning}"
    echo "Number of genes after pruning: ${num_gene_after_pruning}"
	if [[ ${num_gene_after_pruning} -lt 3 ]]; then
		echo 'This is not sufficient for tree-based analysis (<3) . Exitng.'
		exit 0
	fi
		set_analysis_file untrimmed_aln "${file_og_untrimmed_aln_pruned}"
		set_analysis_file trimmed_aln "${file_og_trimmed_aln_pruned}"
		set_analysis_file unrooted_tree "${file_og_unrooted_tree_pruned}"
		set_analysis_file rooted_tree "${file_og_rooted_tree_pruned}"
		set_analysis_file dated_tree "${file_og_dated_tree_pruned}"
fi
if [[ -s ${file_og_expression} && ( ${run_l1ou} -eq 1 || ${run_phylogeneticem} -eq 1 ) ]]; then
  # This block should be run after tree pruning.
  num_gene_trait=$(( $(wc -l < "${file_og_expression}") - 1 )) # -1 for header
  num_gene_tree=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
  if [[ ${num_gene_trait} -eq ${num_gene_tree} ]]; then
      echo "num_gene_trait (${num_gene_trait}) and num_gene_tree (${num_gene_tree}) matched."
  else
      echo "num_gene_trait (${num_gene_trait}) and num_gene_tree (${num_gene_tree}) did not match."
      if [[ ! ${run_tree_pruning} -eq 1 && ${run_phylogeneticem} -eq 1 || ${run_l1ou} -eq 1 ]]; then
        echo "Set run_tree_pruning=1 to run phylogenetic comparative analysis. Exiting."
        exit 1
    fi
  fi
fi

task="Checking whether downstream output files are correctly unpruned/pruned."
if [[ ${check_pruned} -eq 1 ]]; then
    gg_step_start "${task}"
    file_to_remove=(
        "${file_og_mapdnds_dn}"
        "${file_og_mapdnds_ds}"
        "${file_og_hyphy_dnds}"
        "${file_og_gff_info}"
        "${file_og_scm_intron_summary}"
        "${file_og_scm_intron_plot}"
        "${file_og_pem_rdata}"
        "${file_og_pem_tree}"
        "${file_og_pem_regime}"
        "${file_og_pem_leaf}"
        "${file_og_pem_plot}"
        "${file_og_l1ou_fit_rdata}"
        "${file_og_l1ou_fit_conv_rdata}"
        "${file_og_l1ou_fit_tree}"
        "${file_og_l1ou_fit_regime}"
        "${file_og_l1ou_fit_leaf}"
        "${file_og_l1ou_fit_plot}"
        "${file_og_iqtree_anc}"
        "${file_og_csubst_b}"
        "${file_og_csubst_cb_stats}"
        "${file_og_gene_pgls}"
        "${file_og_gene_pgls_plot}"
        "${file_og_species_pgls}"
        "${file_og_species_pgls_plot}"
    )
    for ((i=2; i<=csubst_max_arity; i++)); do
        varname="file_og_csubst_cb_${i}"
        if [[ -n "${!varname:-}" ]]; then
            file_to_remove+=( "${!varname}" )
        fi
    done
    remove_flag=0
    if [[ ${run_tree_pruning} -eq 0 ]]; then
        if [[ -s ${file_og_untrimmed_aln_pruned} ]]; then
            remove_flag=1
        fi
    fi
    if [[ ${run_tree_pruning} -eq 1 ]]; then
        if [[ -s ${file_og_stat_tree} ]]; then
            num_gene_pruned=$(gg_count_fasta_records "${file_og_trimmed_aln_pruned}")
            num_stat_branch_row=$(wc -l < "${file_og_stat_branch}")
            if [[ $((${num_gene_pruned}*2)) -lt ${num_stat_branch_row} ]]; then
                remove_flag=1
            fi
        fi
    fi
    if [[ ${remove_flag} -eq 1 ]]; then
        echo "Downstream output files are inconsistent with the run_tree_pruning setting."
        for file in "${file_to_remove[@]}"; do
            echo "Not found: ${file}"
            if [[ -e ${file} ]]; then
                echo "Deleting: ${file}"
                rm ${file}
            fi
        done
        echo "Completed the deletion of inconsistent files."
    else
        echo "Downstream output files are consistent with the run_tree_pruning setting."
    fi
else
	gg_step_skip "${task}"
fi

task="Parameter estimation for mapdNdS"
disable_if_no_input_file "run_mapdnds_parameter_estimation" ${file_og_rooted_tree_analysis} ${file_og_trimmed_aln_analysis}
if [[ ! -s ${file_og_mapdnds_parameter} && ${run_mapdnds_parameter_estimation} -eq 1 ]]; then
	gg_step_start "${task}"

  nwkit drop --target intnode --support yes --name yes \
	--infile ${file_og_rooted_tree_analysis} \
	--outfile mapdnds_input.nwk

	seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file ./mapdnds_input.fasta

  # F3X4+G4 shouldn not be changed otherwise iqtree2mapnh.py has to be updated.
	build_iqtree_mem_args
	iqtree \
	-s mapdnds_input.fasta \
	-m "GY+F3X4+G4" \
	-te mapdnds_input.nwk \
	-T AUTO \
	--threads-max ${NSLOTS} \
	--seqtype CODON${genetic_code} \
	--prefix ${og_id}.iqtree2mapdNdS \
	"${IQTREE_MEM_ARGS[@]}" \
	--ancestral \
	--seed 12345 \
	--redo

	python ${dir_myscript}/iqtree2mapnh.py \
	--iqtree ${og_id}.iqtree2mapdNdS.iqtree \
	--log ${og_id}.iqtree2mapdNdS.log \
	--state ${og_id}.iqtree2mapdNdS.state \
	--alignment mapdnds_input.fasta \
	--treefile ${og_id}.iqtree2mapdNdS.treefile \
	--rooted_tree mapdnds_input.nwk \
	--genetic_code ${genetic_code}

	if [[ -s "iqtree2mapnh.params" && -s "iqtree2mapnh.nwk" ]]; then
	  echo "iqtree2mapnh was successfully completed."
	  mkdir ${og_id}.mapdnds.parameter
    mv_out "iqtree2mapnh.params" ./${og_id}.mapdnds.parameter
    mv_out "iqtree2mapnh.nwk" ./${og_id}.mapdnds.parameter
    zip -r ${og_id}.mapdnds.parameter.zip ${og_id}.mapdnds.parameter
	    mv_out ${og_id}.mapdnds.parameter.zip ${file_og_mapdnds_parameter}
  else
    echo "iqtree2mapnh.params was not generated."
  fi
else
	gg_step_skip "${task}"
fi

task="mapdNdS main run"
disable_if_no_input_file "run_mapdnds" ${file_og_mapdnds_parameter} ${file_og_trimmed_aln_analysis}
if [[ ( ! -s ${file_og_mapdnds_dn} || ! -s ${file_og_mapdnds_ds} ) && ${run_mapdnds} -eq 1 ]]; then
	gg_step_start "${task}"

	unzip -o ${file_og_mapdnds_parameter}
	cd ${dir_tmp}/${og_id}.mapdNdS.parameter
	seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file ./mapdnds_input.fasta

	mapnh \
	SEQ=mapdnds_input.fasta \
	TREE=iqtree2mapnh.nwk \
	OUT=${og_id} \
	param=iqtree2mapnh.params \
	2>&1 | tee mapnh.log.txt

  if [[ -s ${og_id}.dN.dnd && -s ${og_id}.dS.dnd ]]; then
    echo "mapnh successfully generated dN and dS trees."
    mv_out ${og_id}.dN.dnd ${file_og_mapdnds_dn}
  	mv_out ${og_id}.dS.dnd ${file_og_mapdnds_ds}
  else
    echo "mapnh failed to generate dN and dS trees."
  fi
  cd ${dir_tmp}
else
	gg_step_skip "${task}"
fi

task="CodeML two-ratio model"
disable_if_no_input_file "run_codeml_two_ratio" ${file_og_rooted_tree_analysis} ${file_og_trimmed_aln_analysis} ${file_sp_trait}
if [[ ! -s ${file_og_codeml_two_ratio} && ${run_codeml_two_ratio} -eq 1 ]]; then
	gg_step_start "${task}"

  binarize_species_trait ${file_sp_trait} species_trait_binary.tsv
  sed '2,$ s/\t/_.*\t/' species_trait_binary.tsv > foreground.tsv
  IFS=$'\t' read -r -a colname_array < foreground.tsv

  for ((i=1; i<${#colname_array[@]}; i++)); do
    trait="${colname_array[$i]}"
    echo "Processing trait: ${trait}"
    awk -F'\t' -v trait_col="$((i+1))" 'NR>1 && $trait_col == 1 { print $1 }' foreground.tsv > "foreground_${trait}.txt"

	    target_spnode=$(paste -sd'|' "foreground_${trait}.txt") # Bar-separated list of target species nodes
    echo "Regular expression for CodeML foreground node search: ${target_spnode}"

    nwkit drop \
    --infile ${file_og_rooted_tree_analysis} \
    --target "intnode" \
    --support "yes" \
    --name "yes" \
    | nwkit mark \
    --pattern "${target_spnode}" \
    --insert_txt "#1" \
    --insert_sep "" \
    --target "mrca" \
    --target_only_clade "yes" \
    --outfile "codeml_input_${trait}.nwk"
    echo "CodeML input tree: $(< "codeml_input_${trait}.nwk")"

    bash ${dir_myscript}/shorten_fasta_newick_names.sh \
    "${file_og_trimmed_aln_analysis}" "codeml_input2_${trait}.fasta" "codeml_input_${trait}.nwk" "codeml_input2_${trait}.nwk" 90

    grep -q "#1:0;$" "codeml_input_${trait}.nwk"; exit_code1=$?
    grep -q "#1" "codeml_input_${trait}.nwk"; exit_code2=$?
    flag_unanalyzable=0
    if [[ ${exit_code1} -eq 0 || ${exit_code2} -eq 1 ]]; then
      flag_unanalyzable=1
    fi
    if [[ ${flag_unanalyzable} -eq 1 ]]; then
      echo "Target species tree node (${target_spnode}) is the gene tree root node. Generating an empty output file."
      codeml_out_treelength="NA"
      codeml_out_treelength_dn="NA"
      codeml_out_treelength_ds="NA"
      codeml_out_kappa="NA"
      codeml_out_background_omega="NA"
      codeml_out_foreground_omega="NA"
      codeml_out_time="NA"
    else
      python -c 'from pathlib import Path; import sys; template = Path(sys.argv[1]).read_text(encoding="utf-8"); rendered = template.replace("__SEQFILE__", sys.argv[2]).replace("__TREEFILE__", sys.argv[3]).replace("__ICODE__", sys.argv[4]); Path(sys.argv[5]).write_text(rendered, encoding="utf-8")' \
        "${dir_myscript}/heredoc/codeml_two_ratio.ctl.template" \
        "codeml_input2_${trait}.fasta" \
        "codeml_input2_${trait}.nwk" \
        "$((genetic_code - 1))" \
        "my_codeml_${trait}.ctl"

      codeml "my_codeml_${trait}.ctl"
	      codeml_out_treelength=$(awk '/^tree length =/ {sub(/^tree length =[[:space:]]*/, "", $0); print; exit}' mlc)
	      codeml_out_treelength_dn=$(awk '/^tree length for dN:/ {sub(/^tree length for dN:[[:space:]]*/, "", $0); print; exit}' mlc)
	      codeml_out_treelength_ds=$(awk '/^tree length for dS:/ {sub(/^tree length for dS:[[:space:]]*/, "", $0); print; exit}' mlc)
	      codeml_out_kappa=$(awk '/^kappa \(ts\/tv\) =/ {sub(/^kappa \(ts\/tv\) =[[:space:]]*/, "", $0); print; exit}' mlc)
	      read -r -a codeml_out_omegas <<< "$(awk '/^w \(dN\/dS\) for branches:/ {sub(/^w \(dN\/dS\) for branches:[[:space:]]*/, "", $0); print; exit}' mlc)"
	      codeml_out_background_omega=${codeml_out_omegas[0]}
	      codeml_out_foreground_omega=${codeml_out_omegas[1]}
	      codeml_out_time=$(awk '/^Time used:/ {sub(/^Time used:[[:space:]]*/, "", $0); print; exit}' mlc)
    fi
    if [[ -n ${codeml_out_background_omega} && -n ${codeml_out_foreground_omega} ]]; then
      echo "The task '${task}' has completed successfully for trait '${trait}'."
      echo -e "tree_length_${trait}\ttree_length_dn_${trait}\ttree_length_ds_${trait}\tkappa_${trait}\tbackground_omega_${trait}\tforeground_omega_${trait}\tcodeml_time_${trait}" > "file_og_codeml_two_ratio_${trait}.tsv"
      echo -e "${codeml_out_treelength}\t${codeml_out_treelength_dn}\t${codeml_out_treelength_ds}\t${codeml_out_kappa}\t${codeml_out_background_omega}\t${codeml_out_foreground_omega}\t${codeml_out_time}" >> "file_og_codeml_two_ratio_${trait}.tsv"
    else
      echo "The task '${task}' failed for trait '${trait}'."
    fi
  done

  # Combine all codeml output files
  missing_files=()
  codeml_output_files=()
  for ((i=1; i<${#colname_array[@]}; i++)); do
    trait="${colname_array[$i]}"
    codeml_output_file="file_og_codeml_two_ratio_${trait}.tsv"
    if [[ -s ${codeml_output_file} ]]; then
      codeml_output_files+=(${codeml_output_file})
    else
      missing_files+=(${codeml_output_file})
    fi
  done
  if [[ ${#missing_files[@]} -gt 0 ]]; then
    echo "The following codeml output files are missing:"
    for f in "${missing_files[@]}"; do
      echo "${f}"
    done
    echo "The task has failed: ${task}"
  else
    echo "All codeml two-ratio model output files are generated. Combining them into a single tsv file: ${file_og_codeml_two_ratio}"
    header_files=()
    data_files=()
    for f in "${codeml_output_files[@]}"; do
      base=$(basename "$f")
      head -n1 "$f" > "header_${base}"
      tail -n1 "$f" > "data_${base}"
      header_files+=("header_${base}")
      data_files+=("data_${base}")
    done
	    paste -d$'\t' "${header_files[@]}" > combined_header.tsv
	    paste -d$'\t' "${data_files[@]}" > combined_data.tsv
	    cat combined_header.tsv combined_data.tsv > "${og_id}.codeml.two_ratio.tmp.tsv"
	    mv_out "${og_id}.codeml.two_ratio.tmp.tsv" ${file_og_codeml_two_ratio}
	    echo "The task has completed successfully: ${task}"
	  fi
else
	gg_step_skip "${task}"
fi

task="HyPhy dN-dS estimation"
disable_if_no_input_file "run_hyphy_dnds" ${file_og_rooted_tree_analysis} ${file_og_trimmed_aln_analysis}
if [[ ! -s ${file_og_hyphy_dnds} && ${run_hyphy_dnds} -eq 1 ]]; then
	gg_step_start "${task}"

	nwkit drop --target intnode --support yes --name yes \
	--infile ${file_og_rooted_tree_analysis} \
	--outfile "hyphy_input.nwk"

	seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file "hyphy_input.fasta"

	hyphy_genetic_code=$(get_hyphy_genetic_code ${genetic_code})


  hyphy "${dir_myscript}/hyphy/FitMG94.bf" \
  --alignment "hyphy_input.fasta" \
  --tree "hyphy_input.nwk" \
  --code "${hyphy_genetic_code}" \
  --frequencies "CF3x4" \
  --type "local" \
  --lrt "No" \
  --rooted "Yes" \
  --CPU ${NSLOTS}

  # --lrt "Yes" took too long time for some genes. 20 sec vs 10 min in a small tree.


	mv_out "hyphy_input.fasta.FITTER.json" "${file_og_hyphy_dnds}"
else
	gg_step_skip "${task}"
fi


# HyPhy RELAX
run_hyphy_relax_for_all_traits() {
  local foreground="$1"
  local out_json="$2"

  binarize_species_trait ${file_sp_trait} species_trait_binary.tsv
  sed '2,$ s/\t/_.*\t/' species_trait_binary.tsv > foreground.tsv
  IFS=$'\t' read -r -a colname_array < foreground.tsv

  local reversed_mark=""
  if [[ "${foreground}" == "1" ]]; then
    echo "Running HyPhy RELAX for foreground=1"
  elif [[ "${foreground}" == "0" ]]; then
    echo "Running HyPhy RELAX for foreground=0 (reversed)"
    reversed_mark="_reversed"
  else
    echo "Error: foreground must be 1 or 0"
    echo "Exiting."
    exit 1
  fi

  for ((i=1; i<${#colname_array[@]}; i++)); do
    trait="${colname_array[$i]}"
    echo "Processing trait: ${trait}"
    awk -F'\t' -v trait_col="$((i+1))" -v foreground="${foreground}" 'NR>1 && $trait_col == foreground { print $1 }' foreground.tsv > "foreground_${trait}${reversed_mark}.txt"

	    fg_regex=$(paste -sd'|' "foreground_${trait}${reversed_mark}.txt")
    echo "Foreground node search pattern: ${fg_regex}"
    nwkit drop --target intnode --support yes --name yes --infile ${file_og_rooted_tree_analysis} \
    | nwkit mark --pattern "${fg_regex}" --target "clade" --target_only_clade "yes" --insert_txt "{Foreground}" --outformat 1 --outfile "hyphy_input_${trait}${reversed_mark}.nwk"

    if grep -q "{Foreground}" "hyphy_input_${trait}${reversed_mark}.nwk"; then
      seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file "hyphy_input_${trait}${reversed_mark}.fasta"
      hyphy_genetic_code=$(get_hyphy_genetic_code ${genetic_code})
      hyphy relax \
      --alignment "hyphy_input_${trait}${reversed_mark}.fasta" \
      --tree "hyphy_input_${trait}${reversed_mark}.nwk" \
      --code "${hyphy_genetic_code}" \
      --mode "Classic mode" \
      --test "Foreground" \
      --muliple-hits "Yes" \
      --srv "Yes" \
      --rooted "Yes" \
      --CPU ${NSLOTS}
    else
      echo "No foreground lineage is included in this gene tree. Generating an empty JSON file."
      echo "{}" > "hyphy_input_${trait}${reversed_mark}.fasta.RELAX.json"
    fi
  done

  # Combine all hyphy relax output files
  missing_files=()
  relax_output_files=()
  for ((i=1; i<${#colname_array[@]}; i++)); do
    trait="${colname_array[$i]}"
    relax_output_file="hyphy_input_${trait}${reversed_mark}.fasta.RELAX.json"
    if [[ -s ${relax_output_file} ]]; then
      relax_output_files+=(${relax_output_file})
    else
      missing_files+=(${relax_output_file})
    fi
  done
  if [[ ${#missing_files[@]} -gt 0 ]]; then
    echo "The following HyPhy RELAX output files are missing:"
    for f in "${missing_files[@]}"; do
      echo "${f}"
    done
  else
    echo "All HyPhy RELAX output files are generated. Combining them into a single JSON file: ${out_json}"
    echo "{}" > "combined_relax_output${reversed_mark}.json"
    for ((i=0; i<${#relax_output_files[@]}; i++)); do
      file=${relax_output_files[$i]}
      trait=${colname_array[$((i+1))]}
      jq --arg key "${trait}${reversed_mark}" --slurpfile value "${file}" '. + {($key): $value[0]}' "combined_relax_output${reversed_mark}.json" > "tmp${reversed_mark}.json"
      mv_out "tmp${reversed_mark}.json" "combined_relax_output${reversed_mark}.json"
    done
	    jq . "combined_relax_output${reversed_mark}.json" > "combined_relax_output${reversed_mark}.tmp.json"
	    mv_out "combined_relax_output${reversed_mark}.tmp.json" ${out_json}
	  fi
}

task="HyPhy RELAX"
disable_if_no_input_file "run_hyphy_relax" ${file_og_rooted_tree_analysis} ${file_og_trimmed_aln_analysis} ${file_sp_trait}
if [[ ! -s ${file_og_hyphy_relax} && ${run_hyphy_relax} -eq 1 ]]; then
	gg_step_start "${task}"
  run_hyphy_relax_for_all_traits 1 ${file_og_hyphy_relax}
  if [[ -s ${file_og_hyphy_relax} ]]; then
    echo "The task has completed successfully: ${task}"
  else
    echo "The task has failed: ${task}"
  fi
else
	gg_step_skip "${task}"
fi

task="HyPhy RELAX with reversed foreground/background"
disable_if_no_input_file "run_hyphy_relax_reversed" ${file_og_rooted_tree_analysis} ${file_og_trimmed_aln_analysis} ${file_sp_trait}
if [[ ! -s ${file_og_hyphy_relax_reversed} && ${run_hyphy_relax_reversed} -eq 1 ]]; then
  gg_step_start "${task}"
  run_hyphy_relax_for_all_traits 0 ${file_og_hyphy_relax_reversed}
  if [[ -s ${file_og_hyphy_relax_reversed} ]]; then
    echo "The task has completed successfully: ${task}"
  else
    echo "The task has failed: ${task}"
  fi
else
  gg_step_skip "${task}"
fi

task="Stochastic character mapping of intron evolution"
disable_if_no_input_file "run_scm_intron" ${file_og_gff_info} ${file_og_dated_tree_analysis}
if [[ ! -s ${file_og_scm_intron_summary} && ${run_scm_intron} -eq 1 ]]; then
	gg_step_start "${task}"


  Rscript ${dir_myscript}/scm_intron_evolution.r \
  --tree_file=${file_og_dated_tree_analysis} \
  --trait_file=${file_og_gff_info} \
  --intron_gain_rate=${intron_gain_rate} \
  --retrotransposition_rate=${retrotransposition_rate} \
  --nrep=1000 \
  --nslots=${NSLOTS}

  cp_out intron_evolution_summary.tsv ${file_og_scm_intron_summary}
  if [[ -e intron_evolution_plot.pdf ]]; then
    cp_out intron_evolution_plot.pdf ${file_og_scm_intron_plot}
  fi
else
	gg_step_skip "${task}"
fi

task="PhylogeneticEM"
disable_if_no_input_file "run_phylogeneticem" ${file_og_expression} ${file_og_dated_tree_analysis}
if [[ ( ! -s ${file_og_pem_rdata} || ! -s ${file_og_pem_tree} ||  ! -s ${file_og_pem_regime} || ! -s ${file_og_pem_leaf} ) && ${num_gene_after_maxalign} -gt 3 && ${run_phylogeneticem} -eq 1 ]]; then
	gg_step_start "${task}"

	pem_fit_file=''
	if [[ ${phylogeneticem_use_fit_file} -eq 1 && -s ${file_og_pem_rdata} ]]; then
		pem_fit_file=${file_og_pem_rdata}
	fi

	Rscript ${dir_myscript}/detect_OU_shift_PhylogeneticEM.r \
	--tree_file=${file_og_dated_tree_analysis} \
	--trait_file=${file_og_expression} \
	--nslots=${NSLOTS} \
	--fit_file=${pem_fit_file} \
	--clade_collapse_similarity_method=${clade_collapse_similarity_method} \
	--clade_collapse_similarity_threshold=${clade_collapse_similarity_threshold} \
	--ceil_negative=0 \
  --replicate_sep="_"

	mv_out PhylogeneticEM.tree.tsv ${file_og_pem_tree}
	mv_out PhylogeneticEM.regime.tsv ${file_og_pem_regime}
	mv_out PhylogeneticEM.leaf.tsv ${file_og_pem_leaf}
	mv_out PhylogeneticEM.plot.pdf ${file_og_pem_plot}
	mv_out PhylogeneticEM.out.RData ${file_og_pem_rdata}

else
	gg_step_skip "${task}"
fi

task="l1ou"
disable_if_no_input_file "run_l1ou" ${file_og_trimmed_aln_analysis} ${file_og_expression} ${file_og_dated_tree_analysis}
if [[ ( ! -s ${file_og_l1ou_fit_rdata} || ! -s ${file_og_l1ou_fit_tree} || ! -s ${file_og_l1ou_fit_regime} || ! -s ${file_og_l1ou_fit_leaf} ) && ${run_l1ou} -eq 1 ]]; then
	gg_step_start "${task}"

  num_gene=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
  if [[ ${num_gene} -ge ${large_tree_num_gene} ]]; then
    max_nshift=${large_tree_max_nshift}
  else
    max_nshift=0
  fi

  # taskset is needed because l1ou is not good at handling multiple CPUs.
  if command -v nproc >/dev/null 2>&1; then
    CPU_PER_HOST=$(nproc)
  else
	  CPU_PER_HOST=$(grep -c processor /proc/cpuinfo)
  fi
	echo "CPU_PER_HOST: ${CPU_PER_HOST}"
	cpu_id=$(python -c 'import sys; from numpy import random; a = random.choice(range(int(sys.argv[2])), int(sys.argv[1]), replace=False); print(",".join([str(b) for b in a]))' ${NSLOTS} ${CPU_PER_HOST})
	echo "CPU IDs for l1ou: ${cpu_id}"


	fit_ind_file=''
	if [[ ${l1ou_use_fit_file} -eq 1 && -s ${file_og_l1ou_fit_rdata} ]]; then
		fit_ind_file=${file_og_l1ou_fit_rdata}
	fi

	taskset -c ${cpu_id} Rscript ${dir_myscript}/detect_OU_shift_l1ou.r \
	--max_nshift=${max_nshift} \
	--tree_file=${file_og_dated_tree_analysis} \
	--trait_file=${file_og_expression} \
	--nslots=${NSLOTS} \
	--clade_collapse_similarity_method=${clade_collapse_similarity_method} \
	--clade_collapse_similarity_threshold=${clade_collapse_similarity_threshold} \
	--ceil_negative=0 \
	--criterion=${l1ou_criterion} \
	--nbootstrap=${l1ou_nbootstrap} \
	--fit_ind_file=${fit_ind_file} \
	--fit_conv_file='' \
	--alpha_upper=${l1ou_alpha_upper} \
	--detect_convergence=${l1ou_convergence} \
	--replicate_sep="_"

	mv_out fit_ind.RData ${file_og_l1ou_fit_rdata}
	mv_out l1ou_tree.tsv ${file_og_l1ou_fit_tree}
	mv_out l1ou_regime.tsv ${file_og_l1ou_fit_regime}
	mv_out l1ou_leaf.tsv ${file_og_l1ou_fit_leaf}
	mv_out l1ou_plot.pdf ${file_og_l1ou_fit_plot}
	if [[ ${l1ou_nbootstrap} -gt 0 ]]; then
		cp_out l1ou_bootstrap.tsv ${dir_l1ou_bootstrap}/${l1ou_bootstrap}
	fi
	if [[ ${l1ou_convergence} -eq 1 ]]; then
		cp_out fit_conv.RData ${file_og_l1ou_fit_conv_rdata}
	fi

else
	gg_step_skip "${task}"
fi

task="Gene tree PGLS analysis"
disable_if_no_input_file "run_pgls_gene_tree" ${file_sp_trait} ${file_og_expression} ${file_og_dated_tree_analysis}
if [[ ! -s ${file_og_gene_pgls} && ${run_pgls_gene_tree} -eq 1 ]]; then
	gg_step_start "${task}"

	Rscript ${dir_myscript}/gene_tree_pgls.r \
	--prefix=${og_id} \
	--file_tree=${file_og_dated_tree_analysis} \
	--file_exp=${file_og_expression} \
	--file_trait=${file_sp_trait} \
	--replicate_sep="_" \
	--merge_replicates="yes" \
	2>&1 | tee pgls.log


	mv_out gene_PGLS.tsv ${file_og_gene_pgls}
	mv_out gene_PGLS.barplot.pdf ${file_og_gene_pgls_plot}
else
	gg_step_skip "${task}"
fi

task="Species tree PGLS analysis"
disable_if_no_input_file "run_pgls_species_tree" ${file_sp_trait} ${species_tree_pruned} ${file_og_expression}
if [[ ! -s ${file_og_species_pgls} && ${run_pgls_species_tree} -eq 1 ]]; then
	gg_step_start "${task}"

	Rscript ${dir_myscript}/species_tree_pgls.r \
	--file_sptree=${species_tree_pruned} \
	--file_exp=${file_og_expression} \
	--file_trait=${file_sp_trait} \
	--replicate_sep="_" \
	--exp_value_type=${exp_value_type} \
	--merge_replicates="yes" \
	2>&1 | tee pgls.log


	mv_out species_tree_PGLS.tsv ${file_og_species_pgls}
	mv_out species_tree_PGLS.barplot.pdf ${file_og_species_pgls_plot}
else
	gg_step_skip "${task}"
fi

task="IQ-TREE ancestral codon sequence reconstruction for CSUBST"
disable_if_no_input_file "run_iqtree_anc" ${file_og_trimmed_aln_analysis} ${file_og_rooted_tree_analysis}
if [[ ! -s ${file_og_iqtree_anc} && ${run_iqtree_anc} -eq 1 ]]; then
	gg_step_start "${task}"

	rm -rf csubst.*
  seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file "tmp.csubst.fasta"

  nwkit drop --target intnode --support yes --name no --infile ${file_og_rooted_tree_analysis} \
	| nwkit intersection --seqin tmp.csubst.fasta --seqout csubst.fasta --outfile csubst.nwk

  rm tmp.csubst.fasta

	build_iqtree_mem_args
	iqtree \
	-s csubst.fasta \
	-te csubst.nwk \
	-m ${codon_model} \
	-T AUTO \
	--threads-max ${NSLOTS} \
	--seqtype CODON${genetic_code} \
	--prefix csubst \
	--ancestral \
	--rate \
	"${IQTREE_MEM_ARGS[@]}" \
	--seed 12345 \
	--redo

	if [[ -s csubst.rate && -s csubst.state && -s csubst.treefile ]]; then
		mkdir ${og_id}.iqtree.anc
		mv_out csubst.* ${og_id}.iqtree.anc
		zip -rq ${og_id}.iqtree.anc.zip ${og_id}.iqtree.anc
		mv_out ${og_id}.iqtree.anc.zip ${file_og_iqtree_anc}
		rm -r ${og_id}.iqtree.anc
	fi
else
	gg_step_skip "${task}"
fi

task="CSUBST"
disable_if_no_input_file "run_csubst" ${file_og_iqtree_anc}
if [[ ( ! -s ${file_og_csubst_b} || ! -s ${file_og_csubst_cb_stats} ) && ${run_csubst} -eq 1 ]]; then
	gg_step_start "${task}"

	if [[ -s ${file_sp_trait} ]]; then
		echo "CSUBST foreground specification file: ${file_sp_trait}"
    if head -n1 ${file_sp_trait} | grep -q ' '; then
      echo "Column names should not contain spaces: ${file_sp_trait}"
      echo "Exiting."
      exit 1
    fi
		sed '2,$ s/\t/_.*\t/' ${file_sp_trait} > "foreground.tsv"
		foreground_params="--foreground foreground.tsv --fg_format 2"
	else
		echo 'Foreground specification file was not found. CSUBST will run without it.'
		foreground_params=""
	fi
  rm -rf csubst.*
  rm -rf ${og_id}.*
  unzip -q ${file_og_iqtree_anc}
  csubst_input_base="./${og_id}.iqtree.anc/csubst"

  csubst analyze \
  --genetic_code ${genetic_code} \
  --infile_type "iqtree" \
  --alignment_file "${csubst_input_base}.fasta" \
  --rooted_tree_file "${csubst_input_base}.nwk" \
  --iqtree_treefile "${csubst_input_base}.treefile" \
  --iqtree_state "${csubst_input_base}.state" \
  --iqtree_rate "${csubst_input_base}.rate" \
  --iqtree_iqtree "${csubst_input_base}.iqtree" \
  --iqtree_log "${csubst_input_base}.log" \
  --iqtree_model ${codon_model} \
  --iqtree_redo "no" \
  --max_arity ${csubst_max_arity} \
  --exhaustive_until ${csubst_exhaustive_until} \
  --cutoff_stat "${csubst_cutoff_stat}" \
  --max_combination ${csubst_max_combination} \
  --fg_exclude_wg ${csubst_fg_exclude_wg} \
  --fg_stem_only ${csubst_fg_stem_only} \
  --mg_sister "no" \
  --exclude_sister_pair "yes" \
  --ml_anc "no" \
  --b "yes" \
  --s "no" \
  --cs "no" \
  --cb "yes" \
  --bs "no" \
  --cbs "no" \
  --calc_quantile "no" \
  --omegaC_method "submodel" \
  --asrv "each" \
  --threads ${NSLOTS} \
  --calibrate_longtail "yes" \
  --float_type 32 \
  ${foreground_params}

	  if [[ -s csubst_cb_stats.tsv ]]; then
	    echo "CSUBST was successful."
	    mv_out csubst_b.tsv ${file_og_csubst_b}
	    mv_out csubst_cb_stats.tsv ${file_og_csubst_cb_stats}
	    if [[ ${csubst_max_arity} -gt 2 ]]; then
	      for i in $(seq 2 ${csubst_max_arity}); do
	        if [[ -e "csubst_cb_${i}.tsv" ]]; then
	          my_csubst_file=file_og_csubst_cb_${i}
	          mv_out csubst_cb_${i}.tsv ${!my_csubst_file}
	        fi
	      done
	    fi
  else
    echo "CSUBST failed."
  fi
else
	gg_step_skip "${task}"
fi

task="summary statistics"
is_output_older_than_inputs "^file_og_" ${file_og_tree_plot}; summary_flag=$?
disable_if_no_input_file "run_summary" ${file_og_rooted_tree_analysis}
if [[ ( ${summary_flag} -eq 1 || ! -s ${file_og_stat_branch} || ! -s ${file_og_stat_tree} ) && ${run_summary} -eq 1 ]]; then
	gg_step_start "${task}"

	if [[ -s ${file_og_notung_reconcil} ]]; then
		unzip -qf ${file_og_notung_reconcil}
	fi
  notung_root_log_for_summary="PLACEHOLDER"
  if [[ -d "./${og_id}.notung.root" ]]; then
    notung_log_candidates=()
    mapfile -t notung_log_candidates < <(find "./${og_id}.notung.root" -maxdepth 1 -type f -name "*.ntglog" | sort)
    if [[ ${#notung_log_candidates[@]} -gt 0 ]]; then
      notung_root_log_for_summary="${notung_log_candidates[0]}"
    fi
  fi
  notung_reconcil_stats_for_summary="PLACEHOLDER"
  if [[ -d "./${og_id}.notung.reconcil" ]]; then
    reconcil_stats_candidates=()
    mapfile -t reconcil_stats_candidates < <(find "./${og_id}.notung.reconcil" -maxdepth 1 -type f -name "*.reconciled*.parsable.txt" | sort)
    if [[ ${#reconcil_stats_candidates[@]} -gt 0 ]]; then
      notung_reconcil_stats_for_summary="${reconcil_stats_candidates[0]}"
    fi
  fi
  if [[ ${run_tree_pruning} -eq 1 ]]; then
    generax2orthogroup_statistics="PLACEHOLDER" # generax nhx should be pruned to get used here.
  else
    generax2orthogroup_statistics=${file_og_generax_nhx}
  fi
  seqkit seq --threads "${NSLOTS}" "${file_og_cds_fasta}" --out-file "${og_id}.summary.unaligned.fasta"
  seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.summary.trimmed.fasta"
  seqkit seq --threads "${NSLOTS}" "${file_og_promoter_fasta}" --out-file "${og_id}.summary.promoter.fasta"

	python ${dir_myscript}/orthogroup_statistics.py \
	--species_tree ${species_tree_pruned} \
	--unaligned_aln "${og_id}.summary.unaligned.fasta" \
	--trimal_aln "${og_id}.summary.trimmed.fasta" \
	--unrooted_tree ${file_og_unrooted_tree_analysis} \
	--rooted_tree ${file_og_rooted_tree_analysis} \
	--rooting_log ${file_og_rooted_log} \
	--notung_root_log "${notung_root_log_for_summary}" \
	--notung_reconcil_stats "${notung_reconcil_stats_for_summary}" \
	--dated_tree ${file_og_dated_tree_analysis} \
	--dated_log ${file_og_dated_tree_log} \
	--generax_nhx ${generax2orthogroup_statistics} \
	--hyphy_dnds_json ${file_og_hyphy_dnds} \
	--hyphy_relax_json ${file_og_hyphy_relax} \
	--hyphy_relax_reversed_json ${file_og_hyphy_relax_reversed} \
	--l1ou_tree ${file_og_l1ou_fit_tree} \
	--l1ou_regime ${file_og_l1ou_fit_regime} \
	--l1ou_leaf ${file_og_l1ou_fit_leaf} \
	--phylogeneticem_tree ${file_og_pem_tree} \
	--phylogeneticem_regime ${file_og_pem_regime} \
	--phylogeneticem_leaf ${file_og_pem_leaf} \
	--expression ${file_og_expression} \
	--mapdnds_tree_dn ${file_og_mapdnds_dn} \
	--mapdnds_tree_ds ${file_og_mapdnds_ds} \
	--codeml_tsv ${file_og_codeml_two_ratio} \
	--character_gff ${file_og_gff_info} \
	--fimo ${file_og_fimo} \
	--promoter_fasta "${og_id}.summary.promoter.fasta" \
	--scm_intron ${file_og_scm_intron_summary} \
	--csubst_b ${file_og_csubst_b} \
	--gene_pgls_stats ${file_og_gene_pgls} \
	--species_pgls_stats ${file_og_species_pgls} \
	--rpsblast ${file_og_rpsblast} \
	--uniprot ${file_og_uniprot_annotation} \
	--ncpu ${NSLOTS} \
	--clade_ortholog_prefix ${treevis_clade_ortholog_prefix}
  rm -f "${og_id}.summary.unaligned.fasta" "${og_id}.summary.trimmed.fasta" "${og_id}.summary.promoter.fasta"

	#--csubst_cb_stats ${file_og_csubst_cb_stats} \ # Does not support --arity 3 or larger

	cp_out orthogroup.branch.tsv ${file_og_stat_branch}
	cp_out orthogroup.tree.tsv ${file_og_stat_tree}

else
	gg_step_skip "${task}"
fi

task="stat_branch2tree_plot"
disable_if_no_input_file "run_tree_plot" ${file_og_stat_branch} ${file_og_stat_tree}
if ( [[ ${summary_flag} -eq 1 || ! -s ${file_og_tree_plot} ]] ) && [[ ${run_tree_plot} -eq 1 ]]; then
    gg_step_start "${task}"

	    if [[ ${treevis_clade_ortholog} -eq 1 ]]; then
	        ortholog_prefix=${treevis_clade_ortholog_prefix}
	    else
	        ortholog_prefix=""
	    fi
    cb_path=${file_og_csubst_cb_2/cb_2/cb_ARITY}

    Rscript ${dir_myscript}/stat_branch2tree_plot.r \
    --stat_branch=${file_og_stat_branch} \
    --tree_annotation_dir="${dir_myscript}/tree_annotation" \
    --max_delta_intron_present=${treevis_retrotransposition_delta_intron} \
    --width="7.2" \
    --rel_widths="" \
    --panel1="tree,${treevis_branch_length},${treevis_support_value},${treevis_branch_color},L" \
    --panel2="heatmap,${treevis_heatmap_transform},abs,_,expression_" \
    --panel3="pointplot,no,rel,_,expression_" \
    --panel4="cluster_membership,${treevis_max_intergenic_dist}" \
	    --panel5="tiplabel" \
	    --panel6="signal_peptide" \
	    --panel7="transmembrane_domain" \
	    --panel8="intron_number" \
    --panel9="domain,${file_og_rpsblast}" \
    --panel10="alignment,${file_og_trimmed_aln_analysis},${file_og_untrimmed_aln_analysis}" \
    --panel11="fimo,${promoter_bp},${fimo_qvalue}" \
    --panel12="meme,${file_og_meme}" \
    --panel13="ortholog,${ortholog_prefix},${file_og_dated_tree}" \
    --show_branch_id="yes" \
    --event_method=${treevis_event_method} \
    --species_color_table="PLACEHOLDER" \
    --pie_chart_value_transformation=${treevis_pie_chart_value_transformation} \
    --long_branch_display=${treevis_long_branch_display} \
    --long_branch_ref_quantile=${treevis_long_branch_ref_quantile} \
    --long_branch_detect_ratio=${treevis_long_branch_detect_ratio} \
    --long_branch_cap_ratio=${treevis_long_branch_cap_ratio} \
    --long_branch_tail_shrink=${treevis_long_branch_tail_shrink} \
    --long_branch_max_fraction=${treevis_long_branch_max_fraction} \
    --protein_convergence="100,100,yes,3-${csubst_max_arity},${cb_path},${csubst_cutoff_stat}"

	    if [[ -e "df_fimo.tsv" ]]; then
	    	mv_out "df_fimo.tsv" "${file_og_fimo_collapsed}"
	    fi
	    mv_out stat_branch2tree_plot.pdf ${file_og_tree_plot}
else
	gg_step_skip "${task}"
fi

# Copy parameter files and codes to ${dir_og_parameters} for record
file_params=(
	${file_sp_trait}
	${species_tree}
	${species_tree_pruned}
)
	for file_from in "${file_params[@]}"; do
		file_to="${dir_og_parameters}/$(basename "${file_from}")"
		if [[ ! -e "${file_from}" ]]; then
			continue
		fi
		lock_file="${file_to}.lock"
			if command -v flock >/dev/null 2>&1; then
				exec 7> "${lock_file}"
				flock 7
				filesize_from=$(stat -c%s "${file_from}")
				filesize_to=0
			if [[ -s ${file_to} ]]; then
				filesize_to=$(stat -c%s "${file_to}")
		fi
		if [[ ${filesize_from} -ne ${filesize_to} ]]; then
			echo "Storing important files for record: ${file_to}"
				cp_out "${file_from}" "${file_to}"
			fi
				flock -u 7
				exec 7>&-
			else
				echo "Error: flock command is required but not available."
				exit 1
			fi
		done

cd ${dir_pg}
remove_empty_subdirs ${dir_og}

if [[ -s ${file_og_stat_branch} && -s ${file_og_stat_tree} && -s ${file_og_tree_plot} && ${gg_debug_mode:-0} -eq 0 ]]; then
    echo "Output files detected."
    echo "${file_og_stat_branch}"
    echo "${file_og_stat_tree}"
    echo "${file_og_tree_plot}"
    if [[ ${delete_tmp_dir} -eq 1 ]]; then
        echo "Deleting ${dir_tmp}"
        rm -r ${dir_tmp}
    else
        echo "Leaving ${dir_tmp}"
    fi
    echo "$(date): Exiting Singularity environment"
    exit 8
elif [[ -s ${file_og_stat_branch} && -s ${file_og_stat_tree} && -s ${file_og_tree_plot} && ${gg_debug_mode:-0} -eq 1 ]]; then
    echo "Output files detected & debug mode. Leaving ${dir_tmp}"
else
    echo "Output files not found. Leaving ${dir_tmp}"
fi

###################
echo "$(date): Exiting Singularity environment"
