# %%
cli_args = commandArgs(trailingOnly = TRUE)
mode = ifelse('--debug' %in% cli_args, 'debug', 'batch')
cli_args = cli_args[cli_args != '--debug']

suppressPackageStartupMessages(library(ape, quietly = TRUE))
suppressPackageStartupMessages(library(aplot, quietly = TRUE))
suppressPackageStartupMessages(library(cowplot, quietly = TRUE))
suppressPackageStartupMessages(library(ggimage, quietly = TRUE))
suppressPackageStartupMessages(library(ggmsa, quietly = TRUE))
suppressPackageStartupMessages(library(ggplot2, quietly = TRUE))
suppressPackageStartupMessages(library(ggrepel, quietly = TRUE))
suppressPackageStartupMessages(library(ggtree, quietly = TRUE))
suppressPackageStartupMessages(library(grid, quietly = TRUE))
suppressPackageStartupMessages(library(viridis, quietly = TRUE))
suppressPackageStartupMessages(library(xml2, quietly = TRUE))
suppressPackageStartupMessages(library(svglite, quietly = TRUE))
suppressPackageStartupMessages(library(rkftools, quietly = TRUE))
options(stringsAsFactors = FALSE)


# %%
if (mode == "debug") {

  #dir_pg = '/Users/kf/Dropbox/repos/gg_pipeline/workspace/orthogroup'
  #og_id = 'SP0000270'

  dir_pg = '/Users/kf/Library/CloudStorage/GoogleDrive-kenji.fukushima@nig.ac.jp/My Drive/psnl/collaborators/Jun_Kitano/20240725_teleost/workspace_debug/orthogroup'
  og_id = 'SP0015715'

  tree_annotation_dir = '/Users/kf/Library/CloudStorage/GoogleDrive-kenji.fukushima@nig.ac.jp/My Drive/psnl/repos/gg_pipeline/gg_pipeline/script/tree_annotation'
  dir_pg = ifelse(endsWith(dir_pg, '/'), dir_pg, paste0(dir_pg, '/'))
  setwd(dir_pg)
  args = c()
  args = c(args, paste0('--stat_branch=', dir_pg, 'stat.branch/', og_id, '.stat.branch.tsv'))
  args = c(args, paste0('--max_delta_intron_present=', '-0.5')) # Threshold to define retrotransposition. Need intron SCM output to activate.
  args = c(args, paste0('--width=', '7.2')) # Width size in inches (7.2 for two-columns figures in many journal).
  args = c(args, paste0('--rel_widths=', 'heatmap,0.3')) # Manually set column relative widths. The value of 1.0 roughly corresponds to the original total widths.
  path_aln_candidates = c(
    paste0(dir_pg, 'clipkit/', og_id, '.cds.clipkit.fa.gz'),
    paste0(dir_pg, 'clipkit/', og_id, '.cds.clipkit.fasta'),
    paste0(dir_pg, 'clipkit/', og_id, '.cds.clipkit.fa')
  )
  path_aln = path_aln_candidates[file.exists(path_aln_candidates)][1]
  if (is.na(path_aln)) {
    stop(sprintf('Alignment file was not found in debug mode. Checked: %s', paste(path_aln_candidates, collapse = ', ')))
  }
  path_csubst_input_fasta = path_aln # '/Users/kf/Library/CloudStorage/GoogleDrive-kenji.fukushima@nig.ac.jp/My Drive/psnl/collaborators/Jun_Kitano/20240725_teleost/gg_pipeline/csubst_site_Freshwater_K10_minOCN1.8_minomegaC3.0_minOCNCoD0.0/SP0015715_14_39_82_85_95_107_117_136_141_216/SP0015715.iqtree.anc/csubst.fasta'
  path_rpsblast = paste0(dir_pg, 'rpsblast/', og_id, '.rpsblast.tsv')
  path_ortho_nwk = paste0(dir_pg, 'dated_tree/', og_id, '.dated.nwk')
  #path_ortho_nwk = ''
  args = c(args, paste0('--panel1=', 'tree,bl_rooted,support_unrooted,species,L'))
  args = c(args, paste0('--panel2=', 'tiplabel'))
  args = c(args, paste0('--panel3=', 'heatmap,log2p1,abs,_,expression_'))
  args = c(args, paste0('--panel4=', 'pointplot,no,rel,_,expression_'))
  args = c(args, paste0('--panel5=', 'cluster_membership,100000'))
  args = c(args, paste0('--panel6=', 'signal_peptide'))
  args = c(args, paste0('--panel7=', 'transmembrane_domain'))
  args = c(args, paste0('--panel8=', 'intron_number'))
  args = c(args, paste0('--panel9=', 'domain,', path_rpsblast))
  args = c(args, paste0('--panel10=', 'alignment,', path_aln))
  args = c(args, paste0('--panel11=', 'fimo,3000,0.05'))
  args = c(args, paste0('--panel12=', 'amino_acid_site,1,460:461:464:468:475:477:479,', path_csubst_input_fasta))
  #args = c(args, paste0('--panel13=', 'ortholog,Arabidopsis_thaliana_,', path_ortho_nwk))
  args = c(args, paste0('--panel13=', 'ortholog,Danio_rerio_,', path_ortho_nwk))
  args = c(args, paste0('--protein_convergence=', '0,0,yes,3-10,', paste0(dir_pg, 'csubst.cb_ARITY/', og_id, '.csubst_cb_ARITY.tsv'), ',OCNany2spe,2.0|omegaCany2spe,3.0'))
  #args = c(args, paste0('--protein_convergence=', '100,100,yes,3-10,', paste0(dir_pg, 'csubst.cb_ARITY/', og_id, '.csubst_cb_ARITY.tsv'), ',OCNany2spe,2.0|omegaCany2spe,3.0'))
  args = c(args, paste0('--branch_combination=', '14,39,82,85,95,107,117,136,141,216'))
  args = c(args, paste0('--show_branch_id=', 'yes'))
  args = c(args, paste0('--event_method=', 'species_overlap'))
  args = c(args, paste0('--pie_chart_value_transformation=', 'delog2p1'))
  args = c(args, paste0('--species_color_table=', './species_color.tsv'))
  args = c(args, paste0('--tree_annotation_dir=', tree_annotation_dir))
  options(repr.matrix.max.cols = 1000)
} else if (mode == "batch") {
  args = cli_args
}

cat('arguments:\n')
args = rkftools::get_parsed_args(args, print = TRUE)

Rfiles = list.files(file.path(args[['tree_annotation_dir']], 'R'))
for (Rfile in Rfiles) {
  source(file.path(args[['tree_annotation_dir']], 'R', Rfile))
}

# --panelINT: order-sensitive comma-delimited parameters for plotting panels (INT=1, leftmost). e.g., "--panel1 tree,FOO --panel2 heatmap,BAR --panel3 tiplabel"

# tree: phylogenetic tree. format = tree,BRANCH_LENGTH,NODE_LABEL,BRANCH_COLOR,ORIENTATION
# BRANCH_LENGTH: the column name for branch length in the --stat_branch file. e.g., "bl_dated", "bl_rooted", "mapdnds_omega"
# NODE_LABEL: the column name for node labels in the --stat_branch file. e.g., "support_unrooted", "no" to suppress.
# BRANCH_COLOR: the column name for branch colors. e.g., "species", "mapdnds_omega", "l1ou_regime", "no"
# ORIENTATION: tree orientation. "L" for left-to-right time flow, or "R" in reverse.

# heatmap: trait heatmap. format = heatmap,TRANSFORM,GENEWISE_SCALE,REPLICATE_SEPARATOR,INFILE
# TRANSFORM: value transformation. "no", "log2", "log2p1", "log10", "log10p1"
# GENEWISE_SCALE: whether trait values are scaled in each gene. "abs" for no-scaling and "rel" for scaling
# REPLICATE_SEPARATOR: separator between the experiment name and the replication number in column names. e.g. "_"
# TRAIT_PREFIX: Prefix string specifying the trait columns in the --stat_branch file.

# pointplot: trait point plot. format = pointplot,TRANSFORM,GENEWISE_SCALE,REPLICATE_SEPARATOR,INFILE
# TRANSFORM: value transformation. "no", "log2", "log2p1", "log10", "log10p1"
# GENEWISE_SCALE: whether trait values are scaled in each gene. "abs" for no-scaling and "rel" for scaling
# REPLICATE_SEPARATOR: separator between the experiment name and the replication number in column names. e.g. "_"
# TRAIT_PREFIX: Prefix string specifying the trait columns in the --stat_branch file.

# cluster_membership: Gene cluster membership plot. format = cluster_membership,MAXDIST
# MAXDIST: Maximum distance between genes in bp to call the cluster membership.

# signal_peptide: Stacked bar plot for TargetP results.

# transmembrane_domain: Number of transmembrane domains predicted by TMHMM.

# intron_number: Number of introns in CDS regions originally obtained from gff file.

# tiplabel: gene names. format = tiplabel

# domain: protein domains. format = domain,INFILE
# INFILE: rpsblast output tsv. When running rpsblast, please specify -outfmt "6 std qlen slen stitle"

# alignment: nucleotide or amino acid alignment. format = alignment,INFILE
# INFILE: alignment fasta file. This can be trimmed or untrimmed alignment.

# fimo: promoter cis-element plot. format = fimo,PROMOTER_BP,QVALUE
# PROMOTER_BP: the size of examined promoters
# QVALUE: q-value cutoff for predicted cis-regulatory motifs

# amino_acid_site: Amino acid site plot. format = amino_acid_site,GENETIC_CODE,SITES,ALIGNMENT_PATH
# GENETIC_CODE: genetic code number. 1 for standard, 4 for mold, 11 for the echinoderm mitochondrial code, etc.
# SITES: colon-delimited amino acid site numbers to plot.
# ALIGNMENT_PATH: path to the alignment file. This can be trimmed or untrimmed alignment.

# ortholog: clade orthologs. format = ortholog,ORTHOLOG_PREFIX,INFILE
# ORTHOLOG_PREFIX: gene label prefix to identify species to show orthologs. This string is removed from the clade labels. e.g., "Arabidopsis_thaliana_", "Oryza_sativa_"
# INFILE: "no" or a newick file path. A larger tree can be specified. Number of genes can be different, but the topology can be mergeable to the tree generated from the stat_branch file.

# --protein_convergence: visualization of csubst outputs. format = TOP_PERCENT_TO_SHOW,MAXIMUM_TO_SHOW,FOREGROUND_ONLY,ARITY_RANGE,CB_PATH,CUTOFF_STAT
# TOP_PERCENT_TO_SHOW: how many % of top convergent branch combinations to show.
# MAXIMUM_TO_SHOW: Further limit the number of branch combinations to visualize if TOP_PERCENT_TO_SHOW resulted in greater branch combinations.
# FOREGROUND_ONLY: 'yes' or 'no'. If 'yes', limit the visualization to the foreground taxa.
# ARITY: Minimum-Maximum arity (K) to visualize. e.g. 3-10.
# CB_PATH: PATH to the csubst_cb_?.tsv. If arity is specified as "ARITY", this script automatically finds the cb files for higher-order convergence.
# CUTOFF_STAT: Minimum convergence statistics to visualize. See csubst analyze --cutoff_stat.

# --show_branch_id: 'yes' or 'no'. Whether to show branch_id (a.k.a., 'branch_id') on depicted branches.

# --event_method: Method to annotate branching events. 'auto', 'generax', or 'species_overlap'. 'generax' is preferred if 'auto'.

# --tree_annotation_dir: Directory PATH for the associated local R package

# --pie_chart_value_transformation: 'identity|delog2|delog2p1|delog10|delog10p1'. 


# %%
cat(as.character(Sys.time()), 'Loading infiles.\n')
b = read.table(args[['stat_branch']], header = TRUE, sep = '\t', quote = '', stringsAsFactors = FALSE, fill = TRUE, check.name = FALSE)
b = enhance_branch_table(b, args, event_method = args[['event_method']])

# %%
cat(as.character(Sys.time()), 'Starting tree plot.\n')
args[['font_size']] = 6
args[['font_size_factor']] = 0.352777778 # https://stackoverflow.com/questions/17311917/ggplot2-the-unit-of-size
args[['scale_thickness']] = min(b[['branch_thickness']])
args[['scale_label']] = 'Million years ago'
args[['nodelabel_color']] = 'gray50'
args[['nodepoint_size']] = 0.8
args[['trait_point_size']] = 1
args[['margins']] = c(0.03, 0, 0.2, 0) # cm
args[['trait_axis_label']] = 'Expression'
args[['node_colors']] = c('D' = 'red', 'R' = 'orange', 'S' = 'blue', 'H' = 'darkturquoise')
cat('node_colors:', names(args[['node_colors']]), '=', args[['node_colors']], '\n')

tree_flag = 0
g = list()
for (col in unlist(args[grep("^panel", names(args))])) {
  if (grepl('^tree', col)) {
    dist_col = strsplit(col, ',')[[1]][2]
    nodelabel_col = strsplit(col, ',')[[1]][3]
    branch_color = strsplit(col, ',')[[1]][4]
    orientation = strsplit(col, ',')[[1]][5]
    if (endsWith(branch_color, '_regime')) {
      col_pattern = paste0(sub('_regime$', '', branch_color), '_mu_')
      num_trait = sum(startsWith(colnames(b)[colnames(b) != paste0(col_pattern, 'complementarity')], col_pattern))
      args[['trait_colors']] = get_trait_colors(num_trait, method = 'continuous')
    }
    g = add_tree_column(g, args, b, dist_col, nodelabel_col, branch_color, orientation, 
                        args[['pie_chart_value_transformation']], args[['species_color_table']])
    if (tree_flag == 0) {
      tree = table2phylo(df = b, name_col = 'node_name', dist_col = dist_col)
      if (dist_col == 'bl_dated') {
        try_out = try(rkftools::force_ultrametric(tree, stop_if_larger_change = 0.1), silent = TRUE)
        if (class(try_out) != "try-error") {
          tree = try_out
        }
      }
      names(g)[grepl('^tree,', names(g))] = 'tree'
      tree_flag = 1
    }
    if ('protein_convergence' %in% names(args)) {
      param_split = strsplit(args[['protein_convergence']], ',')[[1]]
      top_percent_to_show = as.numeric(param_split[1])
      max_num_to_show = as.numeric(param_split[2])
      is_target_only = ifelse(param_split[3] == 'yes', TRUE, FALSE)
      arity_range = param_split[4]
      cb_path = param_split[5]
      cutoff_stat = paste(param_split[6:length(param_split)], collapse = ',')
      g[['tree']] = overlay_convergence(g[['tree']], top_percent_to_show, max_num_to_show, is_target_only, 
                                        arity_range, cb_path, cutoff_stat, node_label='highest_arity')
    }
    if ('branch_combination' %in% names(args)) {
      branch_combination = as.numeric(strsplit(args[['branch_combination']], ',')[[1]])
      g[['tree']] = overlay_branch_combination(g[['tree']], branch_combination)
    }
  } else if (grepl('^heatmap', col)) {
    transform = strsplit(col, ',')[[1]][2]
    scale = strsplit(col, ',')[[1]][3]
    replicate_sep = strsplit(col, ',')[[1]][4]
    trait_prefix = strsplit(col, ',')[[1]][5]
    df_trait = get_df_trait(b, transform, scale, trait_prefix)
    df_trait = merge_replicates(trait_table = df_trait, replicate_sep = replicate_sep)
    args[['trait_colors']] = get_trait_colors(ncol(df_trait), method = 'continuous')
    g = add_heatmap_column(g, args, df_trait)
  } else if (grepl('^pointplot', col)) {
    transform = strsplit(col, ',')[[1]][2]
    scale = strsplit(col, ',')[[1]][3]
    replicate_sep = strsplit(col, ',')[[1]][4]
    trait_prefix = strsplit(col, ',')[[1]][5]
    df_trait = get_df_trait(b, transform, scale, trait_prefix)
    num_color = ncol(merge_replicates(trait_table = df_trait, replicate_sep = replicate_sep))
    args[['trait_colors']] = get_trait_colors(num_color, method = 'continuous')
    g = add_pointplot_column(g, args, df_trait, replicate_sep)
  } else if (grepl('^tiplabel', col)) {
    g = add_tiplabel_column(g, args)
  } else if (grepl('^signal_peptide', col)) {
    g = add_signal_peptide_column(g, args)
  } else if (grepl('^transmembrane_domain', col)) {
    g = add_integer_column(g = g, args = args, gname = 'tm', col = 'tmhmm_predhel', xlab = 'TM')
  } else if (grepl('^intron_number', col)) {
    g = add_integer_column(g = g, args = args, gname = 'intron', col = 'num_intron', xlab = 'I')
  } else if (grepl('^cluster_membership', col)) {
    max_bp_membership = as.numeric(strsplit(col, ',')[[1]][2])
    g = add_cluster_membership_column(g = g, args = args, gname = 'cluster_membership', max_bp_membership = max_bp_membership)
  } else if (grepl('^domain', col)) {
    path_rpsblast = strsplit(col, ',')[[1]][2]
    if (file.exists(path_rpsblast)) {
      df_rpsblast = read.table(path_rpsblast, sep = '\t', header = TRUE, stringsAsFactors = FALSE, comment.char = '', quote = '', check.name = FALSE)
    } else {
      df_rpsblast = NULL
    }
    g = add_protein_domain_column(g, args, df_rpsblast)
  } else if (grepl('^alignment', col)) {
    path_seqs = strsplit(col, ',')[[1]][2]
    if (file.exists(path_seqs)) {
      seqs = ape::read.FASTA(path_seqs, type = 'DNA')
    } else {
      seqs = NULL
    }
    g = add_alignment_column(g, args, seqs)
  } else if (grepl('^fimo', col)) {
    xmax = as.integer(strsplit(col, ',')[[1]][2])
    qvalue = as.numeric(strsplit(col, ',')[[1]][3])
    g = add_fimo_column(g, args, qname = 'fimo', xmax = xmax, qvalue, multiple_connection = 'align_from_tss')
  } else if (grepl('^amino_acid_site', col)) {
      genetic_code = as.integer(strsplit(col, ',')[[1]][2])
      selected_amino_acid_sites = as.integer(strsplit(strsplit(col, ',')[[1]][3], ':')[[1]])
      cds_aln_file = strsplit(col, ',')[[1]][4]
      if (file.exists(cds_aln_file)) {
          cds_aln = ape::read.FASTA(cds_aln_file, type = 'DNA')
          aa_aln = ape::trans(x = cds_aln, code = genetic_code, codonstart = 1)
          tidy_aln = ggmsa::tidy_msa(aa_aln)
          g = add_amino_acid_site_column(g, args, tidy_aln, selected_amino_acid_sites)
      } else {
          cat('CDS alignment file not found. Skipping.\n')
      }
  } else if (grepl('^ortholog', col)) {
    ortholog_prefix = strsplit(col, ',')[[1]][2]
    path_ortho_nwk = strsplit(col, ',')[[1]][3]
    g = add_ortholog_column(g, args, tree, ortholog_prefix, path_ortho_nwk)
  } else {
    cat('Unknown --plot specification. Skipping:', col, '\n')
  }
}

rel_widths = get_rel_widths(g, args[['rel_widths']])
base_width = as.numeric(args[['width']])
cp = cowplot::plot_grid(plotlist = g, nrow = 1, align = 'h', axis = 'bt', rel_widths = rel_widths)
if (mode == 'debug') { plot(cp) }
cat('Writing the plot pdf and svg.\n')
height = length(tree[['tip.label']]) / 10
if (mode == 'debug') {
  extensions = c('.pdf', '.svg')
} else if (mode == 'batch') {
  extensions = c('.pdf')
}
for (extension in extensions) {
  cowplot::save_plot(
    filename = paste0('stat_branch2tree_plot', extension),
    plot = cp,
    nrow = 1,
    base_height = max(3, height),
    base_width = base_width,
    units = 'in',
    dpi = 300,
    limitsize = FALSE,
    bg = 'transparent'
  )
}
cat(as.character(Sys.time()), 'stat_branch2tree_plot completed successfully. Exiting.\n')

# %%
