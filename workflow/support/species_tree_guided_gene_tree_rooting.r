cli_args = commandArgs(trailingOnly = TRUE)
args = cli_args
program_start = proc.time()

library(ape)
library(phytools)
library(rkftools)

options(expressions=20000)

cat('arguments:\n')
args = rkftools::get_parsed_args(args, print = TRUE)
args[['ncpu']] = as.integer(args[['ncpu']])

get_root_position_dependent_species_overlap_scores = function(phy, nslots) {
    num_parallel = max(1L, min(as.integer(nslots), nrow(phy$edge)))
    cluster = parallel::makeCluster(num_parallel, type='PSOCK', outfile='')
    on.exit(parallel::stopCluster(cluster), add=TRUE)
    parallel::clusterExport(cluster, varlist = c('phy', 'reroot', 'get_species_overlap_score'), envir = environment())
    so_score = parallel::parLapply(cluster, seq_len(nrow(phy$edge)), function(i) {
        rt = reroot(tree=phy, node.number=phy$edge[i,2])
        get_species_overlap_score(phy=rt, dc_cutoff=0)
    })
    species_overlap_scores = unlist(so_score, use.names = FALSE)
    return(species_overlap_scores)
}

input_tree = read.tree(args[['in_tree']])
unrooted_tree = unroot(input_tree)
midpoint_tree = phytools::midpoint.root(unrooted_tree)
cat('Number of leaves in the input tree:', length(unrooted_tree[['tip.label']]), '\n')

cat('Starting MAD rooting using', args[['ncpu']], 'CPUs.\n')
start = proc.time()
unrooted_newick = write.tree(unrooted_tree)
res = try(MAD_parallel(unrooted_newick=unrooted_newick, output_mode='custom', ncpu=args[['ncpu']]))
if(class(res)=="try-error") {
    cat(res, '\n')
    cat('MAD cannot be completed correctly. Proceeding.\n')
    res = vector(mode='list', 7)
}
mad_tree = read.tree(text=res[[1]])
end = proc.time()
cat('Elapsed time for MAD rooting:', (end - start)[3], 'sec\n')

cat('Starting species overlap search.\n')
start = proc.time()
check_species_overlap_score = FALSE
if (check_species_overlap_score) {
    species_overlap_scores = get_root_position_dependent_species_overlap_scores(phy=unrooted_tree, nslots=nslots)
    cat('Number of root positions with the minimum species overlap score:', sum(species_overlap_scores==min(species_overlap_scores)), '/', length(species_overlap_scores), '\n')
    ind_min_so = c(1:length(species_overlap_scores))[species_overlap_scores==min(species_overlap_scores)]
    cat('Root positions with the minimum species overlap score:', ind_min_so, '\n')
    cat('Species overlap score: minimum:', min(species_overlap_scores), '\n')
    cat('Species overlap score: MAD:', species_overlap_scores[res[[4]]], '\n')
    midpoint_so_score = get_species_overlap_score(midpoint_tree)
    cat('Species overlap score: midpoint:', midpoint_so_score, '\n')
}
if (is.null(res[[5]])) {
    ind_top_ten_mad = NA
} else {
    ind_top_ten_mad = order(res[[5]])[1:10]
    ind_rho_peak = c(1:length(res[[7]]))[(res[[7]]!=0)&(res[[7]]!=1)]
    cat('Root positions with rho peak:', ind_rho_peak, '\n')
    cat('Top 10 MAD positions:', ind_top_ten_mad, '\n')
}

if (! is.null(mad_tree)) {
    if (is_same_root(mad_tree, midpoint_tree)) {
        cat('MAD and midpoint rootings are consistent.\n')
    } else {
        cat('MAD and midpoint rootings are inconsistent.\n')
    }
}
end = proc.time()
cat('Elapsed time for species overlap search:', (end - start)[3], 'sec\n')

cat('Starting the analysis of NOTUNG root positions.\n')
start = proc.time()
zip_members = unzip(args[['notung_root_zip']], list = TRUE)$Name
zip_members = zip_members[nzchar(zip_members)]
dir_candidates = unique(sub('/.*$', '', zip_members))
for (dir_candidate in dir_candidates) {
    candidate_path = file.path(getwd(), dir_candidate)
    if (dir.exists(candidate_path)) {
        unlink(candidate_path, recursive = TRUE)
    }
}
unzip(args[['notung_root_zip']])
dir_notung = NULL
for (dir_candidate in dir_candidates) {
    candidate_path = file.path(getwd(), dir_candidate)
    if (dir.exists(candidate_path)) {
        dir_notung = candidate_path
        break
    }
}
if (is.null(dir_notung)) {
    dir_notung_legacy = file.path(getwd(), sub('.zip$', '', basename(args[['notung_root_zip']])))
    if (dir.exists(dir_notung_legacy)) {
        dir_notung = dir_notung_legacy
    }
}
if (is.null(dir_notung)) {
    nwk_files = character(0)
} else {
    files = list.files(dir_notung)
    nwk_files = files[grep('[0-9]$', files)]
}
notung_roots = vector(mode='numeric', length(nwk_files))
cat('Number of rooted NOTUNG trees:', length(nwk_files), '\n')

is_mid_compatible_with_notung = FALSE
is_mad_compatible_with_notung = FALSE
cl = parallel::makeCluster(args[['ncpu']], type='PSOCK')
if (length(nwk_files) > 0) {
    parallel::clusterExport(cl, varlist = c('dir_notung', 'nwk_files', 'midpoint_tree', 'mad_tree', 'is_same_root'), envir = environment())
    results = parallel::parLapply(cl, seq_along(nwk_files), function(i) {
        cat('Processing', i, 'th NOTUNG tree.\n')
        notung_tree = ape::read.tree(file.path(dir_notung, nwk_files[i]))
        is_mid_compatible = is_same_root(midpoint_tree, notung_tree)
        is_mad_compatible = if (!is.null(mad_tree)) is_same_root(mad_tree, notung_tree) else FALSE
        c(is_mid_compatible, is_mad_compatible)
    })
    results = do.call(rbind, results)
} else {
    results = matrix(FALSE, nrow = 0, ncol = 2)
}
parallel::stopCluster(cl)
is_mid_compatible_with_notung = if (nrow(results) > 0) any(results[, 1]) else FALSE
is_mad_compatible_with_notung = if (nrow(results) > 0) any(results[, 2]) else FALSE
cat("is_mid_compatible_with_notung:", is_mid_compatible_with_notung, "\n")
cat("is_mad_compatible_with_notung:", is_mad_compatible_with_notung, "\n")

if (is_mad_compatible_with_notung) {
    cat('The MAD root position is in NOTUNG root positions. Returning the MAD tree.\n')
    out_tree = mad_tree
} else if (is_mid_compatible_with_notung) {
    cat('The MAD root position is not compatible with NOTUNG results.\n')
    cat('The midpoint root position is found to be compatible with NOTUNG resunts. Returning the midpoint tree.\n')
    out_tree = midpoint_tree
} else if (length(nwk_files) > 0) {
    cat('Neither MAD nor midpoint tree is compatible with NOTUNG results. Returning the first NOTUNG tree.\n')
    out_tree = read.tree(file.path(dir_notung, nwk_files[1]))
} else if (!is.null(mad_tree)) {
    cat('No rooted NOTUNG trees were found. Returning the MAD tree.\n')
    out_tree = mad_tree
} else {
    cat('No rooted NOTUNG trees were found. Returning the midpoint tree.\n')
    out_tree = midpoint_tree
}
cat('Elapsed time for the analysis of NOTUNG root positions:', (proc.time() - start)[3], 'sec\n')

cat('Writing output tree.\n')
out_tree[['node.label']] = NULL
write.tree(out_tree, args[['out_tree']])
cat('Elapsed time for reconciliation-assisted gene tree rooting:', (proc.time() - program_start)[3], 'sec\n')
cat('Reconciliation-assisted gene tree rooting completed.\n')
