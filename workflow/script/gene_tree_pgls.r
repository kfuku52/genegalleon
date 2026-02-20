#!/usr/bin/env Rscript

cli_args = commandArgs(trailingOnly=TRUE)
run_mode = ifelse('--debug' %in% cli_args, 'debug', 'batch')
cli_args = cli_args[cli_args != '--debug']

library(ape)
library(nlme)
library(ggplot2)
library(Rphylopars)
library(rkftools)
script_arg = grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir = if (length(script_arg) > 0) dirname(normalizePath(sub('^--file=', '', script_arg[1]))) else getwd()
source(file.path(script_dir, 'pgls_common.R'))

if (run_mode=='debug') {
    dir_pg = '/Users/kf/Dropbox/collaborators/Yudai_Okuyama/20240606_Asarum_PGLS/workspace_debug'
    og_id = 'HOG0000272'
    args = c()
    args = c(args, paste0('--file_tree=', file.path(dir_pg, 'orthogroup', 'dated_tree', paste0(og_id, '.dated.nwk'))))
    args = c(args, paste0('--file_exp=', file.path(dir_pg, 'orthogroup', 'character.expression', paste0(og_id, '.expression.tsv'))))
    args = c(args, paste0('--file_trait=', file.path(dir_pg, 'species_trait', 'species_trait.tsv')))
    args = c(args, paste0('--replicate_sep=', '_'))
    args = c(args, paste0('--merge_replicates=', 'yes'))
    setwd(dir_pg)
} else {
    args = cli_args
}

cat('arguments:\n')
args = rkftools::get_parsed_args(args, print=TRUE)

#read input files
tree = ape::read.tree(args[['file_tree']])
trait = read.table(args[['file_trait']], header=TRUE, sep='\t')
trait[,'species'] = sub(' ', '_', trait[,'species'])
exp = read.table(args[['file_exp']], header=TRUE, sep='\t')
if (args[['merge_replicates']]=='yes') {
    cat('Expression replicates will be merged. Within-species variation will not be taken into account.\n')
    rownames(exp) = exp[,1]
    exp = exp[,2:ncol(exp), drop=FALSE]
    exp = merge_replicates(exp, args[['replicate_sep']])
    exp = cbind(data.frame(gene_id=rownames(exp), stringsAsFactors=FALSE), exp)
    rownames(exp) = NULL
} else {
    cat('Expression replicates will not be merged. Within-species variation will be taken into account.\n')
}
exp = sort_exp(exp, tree)
expression_bases = get_expression_bases(exp, args[['replicate_sep']])
output_cols = c("R2", "R2adj", "sigma", "Fstat", "pval", "logLik", "AIC", "BIC", "PCC", "p.adj", "trait", "variable", "fit_mode")

#extract trait names
trait_cols = colnames(trait[2:length(trait)])

#calculate mean expression
exp = add_expression_mean_cols(exp, expression_bases)

#prepare species column from the gene_id column
exp[,'species'] = exp[,'gene_id']
exp[,'species'] = sub('_', '|', exp[,'species'])
exp[,'species'] = sub('_.*', '', exp[,'species'])
exp[,'species'] = sub('\\|', '_', exp[,'species'])
exp = sort_exp(exp, tree)
#merge expression and trait matrices and rename the first column to species
df_trait_exp = merge(exp, trait, all.x=TRUE)
df_trait_exp = df_trait_exp[-1]
names(df_trait_exp)[1]='species'

head(df_trait_exp)

cat('Number of genes in the tree:', length(tree[['tip.label']]), '\n')
cat('Number of genes in the table:', nrow(exp), '\n')
#cat('Do names match?', geiger::name.check(tree, exp), '\n')

df_stat = run_phylopars_regression(
    df_trait_exp = df_trait_exp,
    tree = tree,
    trait_cols = trait_cols,
    expression_bases = expression_bases,
    output_cols = output_cols,
    include_foreground_lineage = FALSE,
    verbose_working = TRUE,
    use_phenocov = (args[['merge_replicates']] != 'yes')
)

#output1
df_stat = finalize_pgls_stats(df_stat)
file_name = 'gene_tree_PGLS.tsv'
write.table(df_stat, file_name, row.names=FALSE, sep='\t')
head(df_stat)

#bar plot pval
df_stat[,'trait'] = factor(df_stat[,'trait'], levels=trait_cols)
g = ggplot(df_stat, aes(x=trait, y=-log(pval), fill=variable)) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_hline(yintercept=-log(0.05), linetype="dashed", color = "red", linewidth=0.5)

file_name = 'gene_PGLS.barplot.pdf'
ggsave(file_name, g)
g
