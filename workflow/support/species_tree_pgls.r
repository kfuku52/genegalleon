#!/usr/bin/env Rscript

cli_args = commandArgs(trailingOnly = TRUE)
args = cli_args

library(ape)
library(nlme)
library(Rphylopars)
library(ggplot2)
library(rkftools)
script_arg = grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir = if (length(script_arg) > 0) dirname(normalizePath(sub('^--file=', '', script_arg[1]))) else getwd()
source(file.path(script_dir, 'pgls_common.R'))

cat('arguments:\n')
args = rkftools::get_parsed_args(args, print = TRUE)


sum_per_species = function(exp, value_type) {
  my_fun = function(x) { paste(strsplit(x, '_')[[1]][1:2], collapse = '_') }
  exp_sci_name = sapply(exp[['gene_id']], my_fun)
  exp[, 'gene_id'] = exp_sci_name
  nongeneid_cols = colnames(exp)[(colnames(exp) != 'gene_id')]
  if (value_type == 'identity') {
    agg_fun = sum
  } else if (value_type == 'log') {
    agg_fun = function(x) { sum(exp(x)) }
    return_fun = function(x) { log(x) }
  } else if (value_type == 'log2') {
    agg_fun = function(x) { sum(2**x) }
    return_fun = function(x) { log2(x) }
  } else if (value_type == 'log2p1') {
    agg_fun = function(x) { sum(2**x) - length(x) }
    return_fun = function(x) { log2(x+1) }
  }
  exp = aggregate(exp[, nongeneid_cols, drop = FALSE], by = list(exp[['gene_id']]), agg_fun)
  colnames(exp)[1] = 'gene_id'
  for (col in colnames(exp)[2:ncol(exp)]) {
      exp[[col]] = return_fun(exp[[col]])
  }
  return(exp)
}

generate_dummy_df_stat = function(tree, exp, trait, output_cols) {
  fg_nums = rkftools::count_foreground_lineage(tree, trait)
  trait_cols = colnames(trait[2:length(trait)])
  exp_groups = unique(sub('_[0-9]+$', '', colnames(exp)[2:length(colnames(exp))]))
  df_stat_template = data.frame(matrix(NA, nrow = length(exp_groups), ncol = length(output_cols)))
  colnames(df_stat_template) = output_cols
  df_stat_template[['variable']] = exp_groups
  df_stat = data.frame()
  for (trait_value in trait_cols) {
    df_stat_tmp = df_stat_template
    df_stat_tmp[['trait']] = trait_value
    df_stat_tmp[['num_foreground_lineage']] = fg_nums[[trait_value]]
    df_stat = rbind(df_stat, df_stat_tmp)
  }
  return(df_stat)
}

save_pgls_plot = function(df_stat, font_size = 8, pval_line = 0.05) {
  df_stat[, 'trait'] = factor(df_stat[, 'trait'], levels = trait_cols)
  g = list()
  g[['pcc']] = ggplot(df_stat, aes(x = trait, y = PCC, fill = variable)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    ylim(-1, 1) +
    ylab('Pearson\'s correlation coeffient') +
    xlab('') +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = c(0.98, 0.98),
      legend.justification = c(0.98, 0.98),
      legend.direction = "horizontal",
      legend.box = "horizontal",
    )
  g[['r2']] = ggplot(df_stat, aes(x = trait, y = R2adj, fill = variable)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    ylim(0, 1) +
    ylab('R2 (PGLS)') +
    xlab('') +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = c(0.98, 0.98),
      legend.justification = c(0.98, 0.98),
      legend.direction = "horizontal",
      legend.box = "horizontal",
    )
  g[['pval']] = ggplot(df_stat, aes(x = trait, y = -log10(pval), fill = variable)) +
    geom_hline(yintercept = -log10(pval_line), linetype = "dashed", color = "blue", linewidth = 0.25) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    ylab('-log10(P value) (PGLS)') +
    xlab('') +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = c(0.98, 0.98),
      legend.justification = c(0.98, 0.98),
      legend.direction = "horizontal",
      legend.box = "horizontal",
    )
  for (gname in names(g)) {
    g[[gname]] = g[[gname]] +
      xlab('') +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        text = element_text(size = font_size, color = 'black'),
        axis.text = element_text(size = font_size, color = 'black'),
        strip.text = element_text(size = font_size, color = 'black'),
        axis.title = element_text(size = font_size, color = 'black'),
        legend.text = element_text(size = font_size, color = 'black'),
        legend.title = element_text(size = font_size, face = 'bold', hjust = 0.5),
        legend.background = element_rect(fill = alpha('white', 0.4)),
        legend.box.just = 'bottom',
      )
  }
  cp = cowplot::plot_grid(plotlist = g, ncol = 1, axis = 'bt', align = 'vh')
  extensions = c('.pdf', '.svg')
  for (extension in extensions) {
    cowplot::save_plot(
      filename = 'species_tree_PGLS.barplot.pdf',
      plot = cp,
      nrow = 1,
      base_height = 9.6,
      base_width = 7.2,
      units = 'in',
      dpi = 300,
      limitsize = FALSE,
      bg = 'transparent'
    )
  }
  return(cp)
}

tree = ape::read.tree(args[['file_sptree']])
trait = read.table(args[['file_trait']], header = TRUE, sep = '\t')
trait[, 'species'] = sub(' ', '_', trait[, 'species'])
trait_cols = colnames(trait[2:length(trait)])
output_cols = c("R2", "R2adj", "sigma", "Fstat", "pval", "logLik", "AIC", "BIC", "PCC", "p.adj", "trait", "variable", "fit_mode", 'num_foreground_lineage')

exp = read.table(args[['file_exp']], header = TRUE, sep = '\t')
expression_spp = unique(leaf2species(leaf_names = exp[['gene_id']], use_underbar = TRUE))
remove_tips = tree[['tip.label']][!tree[['tip.label']] %in% expression_spp]
tree = ape::drop.tip(tree, remove_tips)
if (length(expression_spp) == 1) {
  cat('Only one species is available and thus PGLS with the species tree cannot be performed. Generating a dummy output and exiting.\n')
  df_stat = generate_dummy_df_stat(tree, exp, trait, output_cols)
  file_name = 'species_tree_PGLS.tsv'
  write.table(df_stat, file_name, row.names = FALSE, sep = '\t', quote = FALSE)
  q()
}
if (args[['merge_replicates']] == 'yes') {
  cat('Expression replicates will be merged. Within-species variation will not be taken into account.\n')
  rownames(exp) = exp[, 1]
  exp = exp[, 2:ncol(exp), drop=FALSE]
  exp = rkftools::merge_replicates(exp, args[['replicate_sep']])
  exp = cbind(data.frame(gene_id = rownames(exp), stringsAsFactors = FALSE), exp)
  rownames(exp) = NULL
} else {
  cat('Expression replicates will not be merged. Within-species variation will be taken into account.\n')
}

num_gene = nrow(exp)
exp = sum_per_species(exp, value_type = args[['exp_value_type']])
num_sp = nrow(exp)
exp = sort_exp(exp, tree)
colnames(exp)[1] = 'species'
is_retained_cols = apply(exp, 2, function(x) { sum(!is.na(x)) > 1 })
removed_expression_cols = colnames(exp)[!is_retained_cols]
expression_cols = colnames(exp)[is_retained_cols]
expression_cols = expression_cols[2:length(expression_cols)]
exp = exp[, is_retained_cols]
expression_bases = get_expression_bases(exp[, 2:ncol(exp), drop=FALSE], args[['replicate_sep']])
exp_mean = add_expression_mean_cols(exp, expression_bases)
df_trait_exp = merge(exp_mean, trait, all.x = TRUE)
cat('Expression data removed due to insufficuent data points:', paste0(removed_expression_cols), '\n')
cat('Number of genes in the input expression table:', num_gene, '\n')
cat('Number of species in the tree:', length(tree[['tip.label']]), '\n')
cat('Number of species in the processed expression table:', nrow(exp), '\n')

head(df_trait_exp)

df_stat = run_phylopars_regression(
  df_trait_exp = df_trait_exp,
  tree = tree,
  trait_cols = trait_cols,
  expression_bases = expression_bases,
  output_cols = output_cols,
  include_foreground_lineage = TRUE,
  use_phenocov = (args[['merge_replicates']] != 'yes')
)
df_stat = finalize_pgls_stats(df_stat)
file_name = 'species_tree_PGLS.tsv'
write.table(df_stat, file_name, row.names = FALSE, sep = '\t', quote = FALSE)
cp = save_pgls_plot(df_stat, font_size = 8)

cp
head(df_stat)
