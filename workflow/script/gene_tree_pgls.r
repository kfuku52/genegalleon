#!/usr/bin/env Rscript

cli_args = commandArgs(trailingOnly=TRUE)
run_mode = ifelse('--debug' %in% cli_args, 'debug', 'batch')
cli_args = cli_args[cli_args != '--debug']

library(ape)
library(nlme)
library(ggplot2)
library(Rphylopars)
library(rkftools)

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

run_phylopars = function(df_trait_exp, tree, trait_cols, expression_bases) {
    nrow = length(trait_cols)*length(expression_bases)
    stats = c("R2", "R2adj", "sigma", "Fstat", "pval", "logLik", "AIC", "BIC", "PCC", "p.adj", "trait", "variable")
    df_stat = data.frame(matrix(vector(), nrow, length(stats)))
    colnames(df_stat) = stats
    i = 1
    for (trait_col in trait_cols){
        for (expression_base in expression_bases) {
            cat('Working with', trait_col, 'vs', expression_base, '\n')
            expression_cols = colnames(df_trait_exp)[startsWith(colnames(df_trait_exp), expression_base)]
            expression_cols = expression_cols[(expression_cols==expression_base)|(grepl('.*[0-9]$', expression_cols))]      
            explanatory_variables = paste(expression_cols, collapse=' + ')
            explained_variable = trait_col
            formula_string = paste(explained_variable, '~', explanatory_variables)
            cat('PGLS formula:', formula_string, '\n')
                        
            out_phylopars = try(phylopars.lm(
                as.formula(formula_string), 
                trait_data=df_trait_exp, 
                tree=tree, 
                model = "BM", 
                pheno_error=TRUE,
                phylo_correlated=TRUE,
                pheno_correlated=TRUE
            ))
            

            if (class(out_phylopars) != "try-error") {
                for (stat in stats[1:6]){
                    df_stat[i,stat] = out_phylopars[stat]
                }
                df_stat[i,'AIC'] = AIC(out_phylopars)
                df_stat[i,'BIC'] = BIC(out_phylopars)
                PCC = cor(df_trait_exp[paste('mean_', expression_base, sep = "")], df_trait_exp[trait_col], method='pearson', use="complete.obs")
                df_stat[i,'PCC'] = PCC
                df_stat[i,'trait'] = trait_col
                df_stat[i,'variable'] = expression_base
            } else {
                cat('error during phylopars process:', formula_string, '\n')
            }
            i=i+1
        }
    }
    cat('adjusted pvalues are calculated:', '\n')
    cat('method = fdr', '\n')
    cat('length =', length(df_stat[,'pval']), '\n')
    df_stat[,'p.adj'] = p.adjust(df_stat[,'pval'], length(df_stat[,'pval']), method = 'fdr')
    return(df_stat)
}

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

#extract trait names
trait_cols = colnames(trait[2:length(trait)])

#calculate mean expression
for (col in expression_bases) {
    is_col = grepl(col, colnames(exp))
    mean_col = paste('mean_', col, sep='')
    if (sum(is_col)>1) {
        exp[,mean_col] = apply(exp[,is_col], 1, function(x){mean(x, na.rm=TRUE)})
    } else {
        exp[,mean_col] = exp[,is_col]
    }
}

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

df_stat = run_phylopars(df_trait_exp, tree, trait_cols, expression_bases)

#output1
df_stat = df_stat[order(df_stat[,'variable']),]
df_stat = df_stat[apply(df_stat, 1, function(x) !all(is.na(x))), ] # Remove NA-only rows
rownames(df_stat) = NULL
file_name = 'gene_PGLS.tsv'
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
