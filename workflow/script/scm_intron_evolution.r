cat('Starting scm_intron_evolution.\n')
cli_args = commandArgs(trailingOnly=TRUE)
mode = ifelse('--debug' %in% cli_args, 'debug', 'batch')
cli_args = cli_args[cli_args != '--debug']

suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phytools))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(rkftools))

cat(paste('ape version:', packageVersion('ape'), '\n'))
cat(paste('phytools version:', packageVersion('phytools'), '\n'))
cat(paste('rkftools version:', packageVersion('rkftools'), '\n'))

if (mode=="debug") {    
    #dir_pg = '/Users/kef74yk/Dropbox (Personal)/repos/gg_pipeline/workspace/query2family'
    dir_pg = '/Volumes/kfT7/Dropbox/data/Nepenthes_gracilis/20221011_pg_pipeline/workspace/query2family'
    og_id = 'HKT1'
    #setwd(dir_pg)
    args = c()
    args = c(args, paste0('--tree_file=', file.path(dir_pg, 'dated_tree', paste0(og_id, '.dated.nwk'))))
    args = c(args, paste0('--trait_file=', file.path(dir_pg, 'character.gff', paste0(og_id, '.gff.tsv'))))
    #args = c(args, paste0('--tree_file=/Users/kef74yk/Downloads/test_scm_rbind/OG0001897.dated.pruned.nwk'))
    #args = c(args, paste0('--trait_file=/Users/kef74yk/Downloads/test_scm_rbind/OG0001897.intron.tsv'))
    
    args = c(args, paste0('--intron_gain_rate=', '0.0001'))
    args = c(args, paste0('--retrotransposition_rate=', '0.001'))
    args = c(args, paste0('--nrep=', '100'))
    args = c(args, paste0('--nslots=', '4'))
} else if (mode=="batch") {
    args = cli_args
}

cat('arguments:\n')
args = rkftools::get_parsed_args(args, print=TRUE)

tree = read.tree(args[['tree_file']])
#tree = read.newick(args[['tree_file']]) # https://github.com/liamrevell/phytools/issues/118
tree = rkftools::pad_short_edges(tree, threshold=1e-3, external_only=FALSE)
tree = rkftools::force_ultrametric(tree, stop_if_larger_change=0.01)

trait = read.table(args[['trait_file']], header=TRUE, sep="\t", stringsAsFactors=FALSE)
is_intree = (trait[['gene_id']] %in% tree[['tip.label']])
trait = trait[is_intree,]
rownames(trait) = trait[['gene_id']]
trait = trait[tree[['tip.label']],]
rownames(trait) = tree[['tip.label']]
trait[,'gene_id'] = NULL

if (mode=='debug') {
    plot(tree)
    head(trait)
}

tree2 = tree
cols = c('intron_present', 'intron_absent')
num_introns = trait[,'num_intron']
tree2[['tip.label']] = paste(num_introns, tree2[['tip.label']])
names(num_introns) = tree2[['tip.label']]
x = ifelse(num_introns>0, 'intron_present', 'intron_absent')
names(x) = names(num_introns)
if (all(is.na(x))) {
    intron_xmatrix = matrix(0, nrow=length(x), ncol=length(cols))
    colnames(intron_xmatrix) = cols
    rownames(intron_xmatrix) = names(num_introns)
} else {
    intron_xmatrix = to.matrix(x,sort(unique(x)))
    for (chr in cols) {
        if (! chr %in% colnames(intron_xmatrix)) {
            xdf = data.frame(intron_xmatrix)
            xdf[chr] = rep(0, nrow(intron_xmatrix))
            intron_xmatrix = as.matrix(xdf)
        }
    }        
}
is_zero = (apply(intron_xmatrix, 1, sum)==0)
intron_xmatrix[is_zero,] = 1/ncol(intron_xmatrix) # give initial values for missing data

prepare_tree_and_xmatrix = function(tree, trait) {
    tree2 = tree
    cols = c('intron_present', 'intron_absent')
    num_introns = trait[,'num_intron']
    tree2[['tip.label']] = paste(num_introns, tree2[['tip.label']])
    names(num_introns) = tree2[['tip.label']]
    x = ifelse(num_introns>0, 'intron_present', 'intron_absent')
    names(x) = names(num_introns)
    if (all(is.na(x))) {
        intron_xmatrix = matrix(0, nrow=length(x), ncol=length(cols))
        colnames(intron_xmatrix) = cols
        rownames(intron_xmatrix) = names(num_introns)
    } else {
        intron_xmatrix = to.matrix(x,sort(unique(x)))
        for (chr in cols) {
            if (! chr %in% colnames(intron_xmatrix)) {
                xdf = data.frame(intron_xmatrix)
                xdf[chr] = rep(0, nrow(intron_xmatrix))
                intron_xmatrix = as.matrix(xdf)
            }
        }        
    }
    is_zero = (apply(intron_xmatrix, 1, sum)==0)
    intron_xmatrix[is_zero,] = 1/ncol(intron_xmatrix) # give initial values for missing data
    return(list(intron_xmatrix=intron_xmatrix, tree2=tree2))
}

my_time = proc.time()
out = prepare_tree_and_xmatrix(tree, trait)
tree2 = out[['tree2']]
intron_xmatrix = out[['intron_xmatrix']]
out_cols = c('white', 'black')
names(out_cols) = c('intron_absent', 'intron_present')

Q = matrix(
    c(
        -args[['intron_gain_rate']], 
        args[['intron_gain_rate']], 
        args[['retrotransposition_rate']], 
        -args[['retrotransposition_rate']]
    )
    , 2, 2, byrow=TRUE)

colnames(Q) = c('intron_absent', 'intron_present')
rownames(Q) = c('intron_absent', 'intron_present')
nsim_core = ceiling(args[['nrep']]/args[['nslots']])

each_simmap = function(seed, tree, x, nsim, Q) {
    set.seed(seed)
    out = phytools::make.simmap(tree=tree, x=x, Q=Q, nsim=nsim)
    return(out)
}

parallel_simmap = mclapply(
    X=1:args[['nslots']], 
    FUN=each_simmap,
    tree=tree2, 
    x=intron_xmatrix, 
    nsim=nsim_core,
    Q=Q, 
    mc.cores=args[['nslots']]
)
intron_mtrees = do.call(c, parallel_simmap)
if (!("multiSimmap" %in% class(intron_mtrees))) {
    print(intron_mtrees)
    class(intron_mtrees) = c("multiSimmap",class(intron_mtrees))
}
#intron_mtrees = make.simmap(tree2, intron_xmatrix, model='ARD', nsim=args[['nrep']], pi='equal', Q=Q) #Q='empirical', use.empirical=TRUE)
intron_summary = summary(intron_mtrees, plot=FALSE)
cat('Number of simulated scenarios:', length(intron_mtrees), '\n')
cat("Time elapsed for the stochastic character mapping: intron\n")
print(proc.time()-my_time)

get_summary_table = function(tree, trait, intron_summary) {
    cols = c('leaf', 'num_intron')
    nrow = length(tree[['node.label']]) + length(tree[['tip.label']])
    df_out = data.frame(matrix(NA, nrow, length(cols)))
    colnames(df_out) = cols
    node_labels = c(tree[['node.label']], tree[['tip.label']])
    df_out[,'leaf'] = node_labels
    df_out[1:nrow(trait),'num_intron'] = trait[['num_intron']]
    num_intron_present = sum(trait[['num_intron']]>0, na.rm=TRUE)
    num_gene = nrow(trait) - sum(is.na(trait[['num_intron']]))
    is_all_intron_present = (num_intron_present==num_gene)&(num_gene!=0)
    is_all_intron_absent = (num_intron_present==0)&(num_gene!=0)
    if (is_all_intron_present) {
        cat('All genes with input info are in the intron_present state.\n')
        df_out[,'intron_present'] = 1
        df_out[,'intron_absent'] = 0
        return(df_out)
    } else if (is_all_intron_absent) {
        cat('All genes with input info are in the intron_absent state.\n')
        df_out[,'intron_present'] = 0
        df_out[,'intron_absent'] = 1
        return(df_out)
    } else {
        cat('The gene set contains both intron-containing and intron-less genes.\n')
        cat('Preparing the output table from simmap results.\n')
        df_out[,'intron_present'] = intron_summary[['ace']][,'intron_present']
        df_out[,'intron_absent'] = intron_summary[['ace']][,'intron_absent']
        return(df_out)
    }
}

df_out = get_summary_table(tree, trait, intron_summary)
write.table(df_out, file="intron_evolution_summary.tsv", sep="\t", quote=FALSE, row.names=FALSE)
if (mode=='debug') {
    head(df_out)
}

fsize=0.35
my_plot = function() {
    plot(intron_summary, colors=out_cols, fsize=fsize, ftype="reg")
    add.simmap.legend(colors=out_cols, prompt=FALSE, x=0.9*par()[['usr']][1], y=5,fsize=fsize)
}

if (any(apply(intron_xmatrix, 2, sum)==nrow(intron_xmatrix))) {
    cat("Tree plotting was skipped because there is no change in the character state.\n")
} else {
    cat("Plotting the phylogenetic tree.\n")
    pdf("intron_evolution_plot.pdf", width=7.2, height=length(tree[['tip.label']])/8)
    my_plot()
    dev.off()
    if (mode=='debug') {
        my_plot()
    }
}
cat('Ending scm_intron_evolution.\n')
