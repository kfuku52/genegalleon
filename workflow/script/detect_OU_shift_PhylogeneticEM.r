cli_args = commandArgs(trailingOnly=TRUE)
mode = ifelse('--debug' %in% cli_args, 'debug', 'batch')
cli_args = cli_args[cli_args != '--debug']

library(PhylogeneticEM, quietly=TRUE)
library(Rphylopars, quietly=TRUE)
library(rkftools, quietly=TRUE)
options(stringsAsFactors=FALSE)

if (mode=="debug") {
    #og = '1_YABBY'
    og = 'SP0014252' # 'SP0014357'
    dir_pg = '/Users/kef74yk/Dropbox/repos/gg_pipeline/workspace/orthogroup/'
    setwd(dir_pg)
    args = c()
    args = c(args, paste0('--tree_file=', dir_pg, 'dated_tree/', og, '.dated.nwk'))
    args = c(args, paste0('--trait_file=', dir_pg, 'character.expression/', og, '.expression.tsv'))
    args = c(args, paste0('--nslots=', 8))
    args = c(args, paste0('--fit_file=', dir_pg, og, '.PhylogeneticEM.collapsed.RData'))
    args = c(args, paste0('--clade_collapse_similarity_method=', 'pearson'))
    args = c(args, paste0('--clade_collapse_similarity_threshold=', 0.99))
    args = c(args, paste0('--ceil_negative=', 0))
    args = c(args, paste0('--replicate_sep=', '_')) # Separator for replicate IDs. "" for no replicate.
} else if (mode=="batch") {
    args = cli_args
}

cat('arguments:\n')
args = rkftools::get_parsed_args(args, print=TRUE)

my_time = proc.time()
tree = read.tree(args[['tree_file']])
cat('Number of leaves:', length(tree$tip.label), '\n')
tree_original = tree
trait_table = read.table(args[['trait_file']], header=TRUE, row.names=1, sep="\t")
trait_table = merge_replicates(trait_table, args[['replicate_sep']])

if ((args[['clade_collapse_similarity_method']]=='complementarity')|(args[['ceil_negative']])) {
    if (args[['clade_collapse_similarity_method']]=='complementarity') {
        cat('clade_collapse_similarity_method is set to complementary: forced to ceil negative values.\n')
    }
    cat('Minimum trait value =', min(trait_table), ', Negative values are ceiled to 0.\n')
    if (min(trait_table) < 0) {
        trait_table[trait_table < 0] = 0
    }
}

tree = pad_short_edges(tree, threshold=1e-3, external_only=FALSE)
tree = force_ultrametric(tree, stop_if_larger_change=0.01)

original_trait_table = trait_table
original_traits = colnames(original_trait_table)

out = remove_invariant_traits(trait_table, small_dif=0.001)
trait_table = out[['trait_table']]
removed_traits = out[['removed_traits']]
if (ncol(trait_table)==0) {
    cat('All traits are invariant. Generating placeholder outputs without running l1ou.\n')
    leaf_table = get_placeholder_leaf(tree, original_trait_table)
    regime_table = get_placeholder_regime(tree, original_trait_table)
    tree_table = get_placeholder_tree(tree, original_trait_table)
    write.table(tree_table, file="l1ou_tree.tsv", sep="\t", quote=FALSE, row.names=FALSE) 
    write.table(regime_table, file="l1ou_regime.tsv", sep="\t", quote=FALSE, row.names=FALSE) 
    write.table(leaf_table, file="l1ou_leaf.tsv", sep="\t", quote=FALSE, row.names=FALSE)
    res = NA
    save(res, file='PhylogeneticEM.RData')
    cat('Placeholder outputs were generated. Exiting.\n')
    quit(save='no', status=0)
    stop('No available trait. This line may be printed in debugging, but should not appear in the batch mode.\n')
}

### Phylogenetic imputation
trait_table = phylogenetic_imputation(tree, trait_table)

if (args[['clade_collapse_similarity_threshold']]) {
    cat('Start cluster collapsing by similarity thresholding (cutoff =', args[['clade_collapse_similarity_threshold']], ')\n')
    high_comple_clades = get_high_similarity_clades(tree, trait_table, args[['clade_collapse_similarity_method']], 
                            args[['clade_collapse_similarity_threshold']], verbose=FALSE, num_test=0)
    out = collapse_clades(tree, trait_table, high_comple_clades)
    tree = out[['tree']]
    trait_table = out[['trait']]
} else {
    cat('Skipped cluster collapsing.', '\n')
}

### pad_short_branch() practically let this chunk obsoleted
tree_bifurcation = tree
tree = collapse_short_branches(tree, tol=1e-8)

trait_table = trait_table[tree$tip.label,]
trait_matrix = t(as.matrix(trait_table))
trait_matrix = trait_matrix[,tree$tip.label]
if (mode=="debug") {
    trait_matrix
    plot(tree, show.node.label=TRUE)
}

cat("Time elapsed for data preparation:\n")
print(proc.time()-my_time)

my_time = proc.time()

flag = FALSE
if ((!is.null(args[['fit_file']]))&(!is.na(args[['fit_file']]))) {
    if (file.exists(args[['fit_file']])) {
        if (file.info(args[['fit_file']])$size) {
            load(args[['fit_file']])
            if(class(res)!="try-error") {
                cat(args[['fit_file']], 'was found. Loading the shift configurations.\n')
                flag = TRUE
            } else {
                cat(args[['fit_file']], 'is malformed. Printing...\n')
                print(res)
                rm(res)
            }
        } else {
            cat(args[['fit_file']], 'has size zero.\n')
        }
    }
}

if (!flag) {
    res = PhyloEM(
        phylo=tree,
        Y_data=trait_matrix,
        process="scOU",
        #K_max=min(max(10, floor(sqrt(length(tree$tip.label)))), length(tree$tip.label)/2),
        K_max=max(10, floor(sqrt(length(tree$tip.label)))),
        methods.segmentation='lasso',
        independent=FALSE,
        random.root=TRUE,
        stationary.root=TRUE,
        alpha_grid=TRUE,
        light_result=TRUE,
        parallel_alpha=TRUE,
        Ncores=args[['nslots']],
        tol_tree=1e-8
    )
}
save(res, file='PhylogeneticEM.out.RData')

cat('Time elapsed for PhyloEM:\n')
print(proc.time()-my_time)

my_time = proc.time()

tree_table = get_tree_table(res, mode='PhylogeneticEM')
regime_table = get_regime_table(res, mode='PhylogeneticEM')
leaf_table = get_leaf_table(res, mode='PhylogeneticEM')

if (contains_polytomy(res[['phylo']])) {
    df_node_match = multi2bi_node_number_transfer(multifurcated_tree=res[['phylo']], bifurcated_tree=tree_bifurcation)
    for (i in 1:nrow(df_node_match)) {
        mtree_node_name = get_node_name_by_num(res[['phylo']], df_node_match[i,'mtree_node'])
        btree_node_name = get_node_name_by_num(tree_original, df_node_match[i,'btree_node'])
        if (mtree_node_name%in%regime_table$node_name) {
            cat('regime_table: polytomy node [', mtree_node_name, '] will be replaced with original binary node [', btree_node_name, ']\n')
            regime_table[regime_table$node_name==mtree_node_name,'node_name'] = btree_node_name
        } else {
            cat('regime_table: polytomy node [', mtree_node_name, '] was not present in regime_table (corresponding binary node [', btree_node_name, '])\n')
        }
    }
}

if (args[['clade_collapse_similarity_threshold']]) {
    cat('Starting the recovery of collapsed nodes.\n')
    node_num_mapping = map_node_num(tree_original, tree_bifurcation, out[['collapse_leaf_names']], verbose=FALSE)
    tree_table = tree_table_collapse2original(tree_table, tree_original)
    regime_table = regime_table_collapse2original(regime_table, tree_original, tree, node_num_mapping)
    leaf_table = leaf_table_collapse2original(leaf_table, tree_original, tree, node_num_mapping)
    leaf_table = restore_imputed_leaves(leaf_table, original_trait_table)
}

cat('Time elapsed for result processing:\n')
print(proc.time()-my_time)

write.table(tree_table, 'PhylogeneticEM.tree.tsv', row.names=FALSE, sep='\t', quote=FALSE)
write.table(regime_table, 'PhylogeneticEM.regime.tsv', row.names=FALSE, sep='\t', quote=FALSE)
write.table(leaf_table, 'PhylogeneticEM.leaf.tsv', row.names=FALSE, sep='\t', quote=FALSE)

pdf("PhylogeneticEM.plot.pdf", height=length(tree$tip.label)/10+1, width=nrow(trait_matrix)+5)
plot(res, show.tip.label=TRUE, label_cex=0.5)
graphics.off()
cat('PhylogeneticEM completed!\n')

if (mode=='debug') {
    plot(res, show.tip.label=TRUE, label_cex=0.2)
}
