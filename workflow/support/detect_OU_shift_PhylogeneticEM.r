cli_args = commandArgs(trailingOnly=TRUE)
args = cli_args

library(PhylogeneticEM, quietly=TRUE)
library(Rphylopars, quietly=TRUE)
library(rkftools, quietly=TRUE)
options(stringsAsFactors=FALSE)

cat('arguments:\n')
args = rkftools::get_parsed_args(args, print=TRUE)

safe_condition_message = function(obj) {
    cond = attr(obj, "condition")
    if (inherits(cond, "condition")) {
        return(conditionMessage(cond))
    }
    return(as.character(obj))
}

parse_bool_flag = function(x, default = 0L) {
    if (is.null(x) || (length(x) == 0) || is.na(x) || (x == "")) {
        return(as.integer(default))
    }
    if (is.logical(x)) {
        return(as.integer(x[[1]]))
    }
    sval = tolower(as.character(x[[1]]))
    if (sval %in% c("1", "true", "t", "yes", "y")) {
        return(1L)
    }
    if (sval %in% c("0", "false", "f", "no", "n")) {
        return(0L)
    }
    num = suppressWarnings(as.numeric(sval))
    if (is.finite(num)) {
        return(as.integer(num != 0))
    }
    warning("Unrecognized boolean flag value: ", sval, ". Falling back to default=", default)
    return(as.integer(default))
}

validate_internal_node_labels = function(tree, context = "tree") {
    num_node = tree$Nnode
    labels = rep(NA_character_, num_node)
    if (!is.null(tree$node.label) && (length(tree$node.label) > 0)) {
        raw_labels = as.character(tree$node.label)
        copy_n = min(length(raw_labels), num_node)
        labels[seq_len(copy_n)] = raw_labels[seq_len(copy_n)]
    }

    issues = c()
    missing_idx = which(is.na(labels) | (labels == ""))
    if (length(missing_idx) > 0) {
        issues = c(issues, sprintf("%d internal node label(s) are missing/empty", length(missing_idx)))
    }

    non_missing = labels[!(is.na(labels) | (labels == ""))]
    if (anyDuplicated(non_missing)) {
        dup = unique(non_missing[duplicated(non_missing)])
        issues = c(issues, paste0("duplicate internal node labels: ", paste(dup, collapse = ", ")))
    }

    overlap = intersect(unique(non_missing), unique(tree$tip.label))
    if (length(overlap) > 0) {
        issues = c(issues, paste0("internal labels overlap with tip labels: ", paste(overlap, collapse = ", ")))
    }

    if (length(issues) > 0) {
        stop(
            paste0(
                "`--require_internal_node_labels=1` failed for ", context, ": ",
                paste(issues, collapse = "; "),
                "."
            )
        )
    }
}

safe_force_ultrametric = function(tree, tolerances = c(0.01, 0.05, 0.1, 1, Inf)) {
    sanitize_tree_edges = function(x, min_len = 1e-8) {
        bad = which(!is.finite(x$edge.length) | (x$edge.length <= 0))
        if (length(bad) > 0) {
            cat("Replacing", length(bad), "non-positive edge lengths with", min_len, "\n")
            x$edge.length[bad] = min_len
        }
        return(x)
    }
    finalize_ultrametric = function(x) {
        x = sanitize_tree_edges(x)
        if (!isTRUE(ape::is.ultrametric(x))) {
            retry = try(force_ultrametric(x, stop_if_larger_change = Inf), silent = TRUE)
            if (!inherits(retry, "try-error")) {
                x = sanitize_tree_edges(retry)
            } else {
                cat("Final ultrametric recheck failed:", safe_condition_message(retry), "\n")
            }
        }
        if (!isTRUE(ape::is.ultrametric(x))) {
            chronos_try = try(ape::chronos(x, lambda = 1, quiet = TRUE), silent = TRUE)
            if (!inherits(chronos_try, "try-error")) {
                x = sanitize_tree_edges(chronos_try)
                cat("Applied ape::chronos fallback to enforce ultrametricity.\n")
            } else {
                cat("ape::chronos fallback failed:", safe_condition_message(chronos_try), "\n")
            }
        }
        if (!isTRUE(ape::is.ultrametric(x))) {
            cat("Tree remained non-ultrametric after fallback correction.\n")
        }
        return(x)
    }
    out_tree = tree
    for (tol in tolerances) {
        trial = try(force_ultrametric(out_tree, stop_if_larger_change = tol), silent = TRUE)
        if (!inherits(trial, "try-error")) {
            if (!identical(tol, tolerances[[1]])) {
                cat("force_ultrametric succeeded with relaxed threshold:", tol, "\n")
            }
            return(finalize_ultrametric(trial))
        }
        cat("force_ultrametric failed at threshold", tol, ":", safe_condition_message(trial), "\n")
    }
    cat("force_ultrametric failed for all thresholds; using padded tree as-is.\n")
    return(finalize_ultrametric(out_tree))
}

ensure_internal_node_names = function(tree, prefix = "n") {
    num_tip = length(tree$tip.label)
    num_node = tree$Nnode
    node_labels = tree$node.label
    if (is.null(node_labels)) {
        node_labels = rep("", num_node)
    }
    if (length(node_labels) < num_node) {
        node_labels = c(node_labels, rep("", num_node - length(node_labels)))
    }
    node_labels = as.character(node_labels)
    for (i in seq_len(num_node)) {
        if (is.na(node_labels[[i]]) || node_labels[[i]] == "") {
            node_labels[[i]] = paste0(prefix, num_tip + i)
        }
    }
    tree$node.label = make.unique(node_labels, sep = "_")
    return(tree)
}

stable_jitter_matrix = function(mat, scale_factor = 1e-6, seed = 1L) {
    out = mat
    if (length(out) == 0) {
        return(out)
    }
    row_sd = apply(out, 1, sd, na.rm = TRUE)
    row_sd[!is.finite(row_sd)] = 0
    noise_sd = pmax(row_sd * scale_factor, 1e-12)
    set.seed(seed)
    noise = matrix(rnorm(length(out)), nrow = nrow(out), ncol = ncol(out))
    noise = noise * noise_sd
    rownames(noise) = rownames(out)
    colnames(noise) = colnames(out)
    out = out + noise
    return(out)
}

compress_trait_matrix = function(trait_matrix, max_components) {
    if (nrow(trait_matrix) <= max_components) {
        return(trait_matrix)
    }
    x = t(trait_matrix)
    x_centered = scale(x, center = TRUE, scale = FALSE)
    sv = svd(x_centered)
    num_comp = min(max_components, ncol(sv$u), length(sv$d))
    if (num_comp < 1) {
        return(trait_matrix)
    }
    scores = sv$u[, seq_len(num_comp), drop = FALSE] %*% diag(sv$d[seq_len(num_comp)])
    out = t(scores)
    rownames(out) = paste0("pc", seq_len(nrow(out)))
    colnames(out) = colnames(trait_matrix)
    return(out)
}

run_phyloem_with_retries = function(tree, trait_matrix, nslots) {
    n_tip = length(tree$tip.label)
    if (n_tip < 15) {
        k_base_target = as.integer(max(4, floor(sqrt(n_tip))))
    } else {
        k_base_target = as.integer(max(10, floor(sqrt(n_tip))))
    }
    k_base = max(2L, min(k_base_target, as.integer(n_tip - 1)))
    k_small = max(2L, min(as.integer(max(3, floor(sqrt(n_tip)))), as.integer(n_tip - 1)))
    k_tiny = max(2L, min(3L, as.integer(n_tip - 1)))
    p_reduced = max(2L, min(as.integer(n_tip - 2), as.integer(floor(sqrt(n_tip)))))
    mat_jitter = stable_jitter_matrix(trait_matrix, scale_factor = 1e-6, seed = 1L)
    mat_reduced = compress_trait_matrix(mat_jitter, max_components = p_reduced)

    attempts = list(
        list(
            label = "baseline_lasso",
            Y_data = trait_matrix,
            K_max = k_base,
            method = "lasso",
            random_root = TRUE,
            parallel_alpha = FALSE,
            alpha_grid = TRUE
        ),
        list(
            label = "jitter_lasso_smallK",
            Y_data = mat_jitter,
            K_max = k_small,
            method = "lasso",
            random_root = TRUE,
            parallel_alpha = FALSE,
            alpha_grid = TRUE
        ),
        list(
            label = "reduced_lasso_smallK",
            Y_data = mat_reduced,
            K_max = k_small,
            method = "lasso",
            random_root = TRUE,
            parallel_alpha = FALSE,
            alpha_grid = TRUE
        ),
        list(
            label = "reduced_same_shifts_smallK",
            Y_data = mat_reduced,
            K_max = k_small,
            method = "same_shifts",
            random_root = TRUE,
            parallel_alpha = FALSE,
            alpha_grid = TRUE
        ),
        list(
            label = "reduced_best_single_move_tinyK_fixedRoot",
            Y_data = mat_reduced,
            K_max = k_tiny,
            method = "best_single_move",
            random_root = FALSE,
            parallel_alpha = FALSE,
            alpha_grid = TRUE
        )
    )

    last_try_error = NULL
    for (attempt in attempts) {
        cat("PhyloEM attempt:", attempt$label, "| p=", nrow(attempt$Y_data), "n=", ncol(attempt$Y_data), "K_max=", attempt$K_max, "method=", attempt$method, "\n")
        res_try = try(
            PhyloEM(
                phylo = tree,
                Y_data = attempt$Y_data,
                process = "scOU",
                K_max = attempt$K_max,
                methods.segmentation = attempt$method,
                independent = FALSE,
                random.root = attempt$random_root,
                stationary.root = TRUE,
                alpha_grid = attempt$alpha_grid,
                light_result = TRUE,
                parallel_alpha = attempt$parallel_alpha,
                Ncores = ifelse(attempt$parallel_alpha, nslots, 1),
                tol_tree = 1e-8
            ),
            silent = TRUE
        )
        if (!inherits(res_try, "try-error")) {
            cat("PhyloEM succeeded with attempt:", attempt$label, "\n")
            return(list(res = res_try, trait_matrix = attempt$Y_data, label = attempt$label))
        }
        last_try_error = res_try
        cat("PhyloEM failed for attempt", attempt$label, ":", safe_condition_message(res_try), "\n")
    }
    return(list(res = last_try_error, trait_matrix = trait_matrix, label = "all_failed"))
}

safely_get_phyloem_tables = function(res, tree, trait_table) {
    tree_try = try(get_tree_table(res, mode = 'PhylogeneticEM'), silent = TRUE)
    regime_try = try(get_regime_table(res, mode = 'PhylogeneticEM'), silent = TRUE)
    leaf_try = try(get_leaf_table(res, mode = 'PhylogeneticEM'), silent = TRUE)

    if (inherits(tree_try, "try-error")) {
        cat("get_tree_table(PhylogeneticEM) failed:", safe_condition_message(tree_try), "\n")
        tree_table = get_placeholder_tree(tree, trait_table)
    } else {
        tree_table = tree_try
    }
    if (inherits(regime_try, "try-error")) {
        cat("get_regime_table(PhylogeneticEM) failed:", safe_condition_message(regime_try), "\n")
        regime_table = get_placeholder_regime(tree, trait_table)
    } else {
        regime_table = regime_try
    }
    if (inherits(leaf_try, "try-error")) {
        cat("get_leaf_table(PhylogeneticEM) failed:", safe_condition_message(leaf_try), "\n")
        leaf_table = get_placeholder_leaf(tree, trait_table)
    } else {
        leaf_table = leaf_try
    }

    regime_set = sort(unique(suppressWarnings(as.integer(regime_table$regime))))
    regime_set = regime_set[!is.na(regime_set)]
    leaf_set = sort(unique(suppressWarnings(as.integer(leaf_table$regime))))
    leaf_set = leaf_set[!is.na(leaf_set)]
    if (!setequal(c(0L, regime_set), leaf_set)) {
        cat("Regime set mismatch between regime and leaf tables. Writing placeholder OU tables.\n")
        regime_table = get_placeholder_regime(tree, trait_table)
        leaf_table = get_placeholder_leaf(tree, trait_table)
    }
    return(list(tree_table = tree_table, regime_table = regime_table, leaf_table = leaf_table))
}

write_safe_tree_plot = function(file, tree, height, width, tip_cex = 0.5) {
    pdf(file, height = max(2, height), width = max(6, width))
    par(mar = c(2, 2, 2, 2))
    plot_try = try(plot(tree, show.tip.label = TRUE, cex = tip_cex), silent = TRUE)
    if (inherits(plot_try, "try-error")) {
        plot.new()
        text(0.5, 0.5, labels = "placeholder")
    }
    graphics.off()
}

my_time = proc.time()
tree = read.tree(args[['tree_file']])
require_internal_node_labels = parse_bool_flag(args[['require_internal_node_labels']], default = 0L)
cat('require_internal_node_labels =', require_internal_node_labels, '\n')
if (require_internal_node_labels == 1L) {
    validate_internal_node_labels(tree, context = args[['tree_file']])
}
tree = ensure_internal_node_names(tree)
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
tree = safe_force_ultrametric(tree, tolerances = c(0.01, 0.05, 0.1, 1, Inf))
tree = ensure_internal_node_names(tree)

original_trait_table = trait_table
original_traits = colnames(original_trait_table)

write_placeholder_outputs_and_quit = function(tree, trait_table, trait_matrix = NULL, reason = "") {
    if (nzchar(reason)) {
        cat(reason, "\n")
    }
    leaf_table = get_placeholder_leaf(tree, trait_table)
    regime_table = get_placeholder_regime(tree, trait_table)
    tree_table = get_placeholder_tree(tree, trait_table)
    write.table(tree_table, file="PhylogeneticEM.tree.tsv", sep="\t", quote=FALSE, row.names=FALSE) 
    write.table(regime_table, file="PhylogeneticEM.regime.tsv", sep="\t", quote=FALSE, row.names=FALSE) 
    write.table(leaf_table, file="PhylogeneticEM.leaf.tsv", sep="\t", quote=FALSE, row.names=FALSE)
    plot_width = 8
    if (!is.null(trait_matrix)) {
        plot_width = max(8, nrow(trait_matrix) + 5)
    }
    write_safe_tree_plot("PhylogeneticEM.plot.pdf", tree, height = length(tree$tip.label)/10+1, width = plot_width, tip_cex = 0.5)
    res = NA
    save(res, file='PhylogeneticEM.out.RData')
    cat('Placeholder outputs were generated. Exiting.\n')
    quit(save='no', status=0)
}

out = remove_invariant_traits(trait_table, small_dif=0.001)
trait_table = out[['trait_table']]
removed_traits = out[['removed_traits']]
if (ncol(trait_table)==0) {
    write_placeholder_outputs_and_quit(
        tree = tree,
        trait_table = original_trait_table,
        reason = 'All traits are invariant. Generating placeholder outputs without running PhylogeneticEM.'
    )
    stop('No available trait. This line may be printed in debugging, but should not appear in the batch mode.\n')
}

impute_traits_with_fallback = function(tree, trait_table) {
    trait_table2 = data.frame(
        species = as.character(rownames(trait_table)),
        trait_table,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    trait_table2 = sort_exp(trait_table2, tree, col = "species")
    for (j in 2:ncol(trait_table2)) {
        trait_table2[[j]] = suppressWarnings(as.numeric(trait_table2[[j]]))
    }
    num_missing = sum(is.na(trait_table2[, -1, drop = FALSE]))
    num_all = nrow(trait_table2) * (ncol(trait_table2) - 1)
    cat("Phylogenetic imputation with phylopars:", num_missing, "/", num_all, "traits will be imputed.\n")

    should_try_correlated = function(tbl, phy) {
        n_tip = length(phy$tip.label)
        n_trait = ncol(tbl) - 1
        if (n_trait <= 1) {
            return(TRUE)
        }
        if (n_trait >= (n_tip - 1L)) {
            cat("Skipping correlated-imputation attempts: num_traits (", n_trait, ") >= num_tips - 1 (", n_tip - 1L, ").\n", sep = "")
            return(FALSE)
        }
        miss = !is.na(as.matrix(tbl[, -1, drop = FALSE]))
        overlap = t(miss) %*% miss
        overlap_lower = overlap[lower.tri(overlap, diag = FALSE)]
        if (length(overlap_lower) > 0) {
            frac_sparse = mean(overlap_lower < 3L)
            if (is.finite(frac_sparse) && (frac_sparse > 0.25)) {
                cat("Skipping correlated-imputation attempts: sparse trait overlap (fraction pairwise overlap <3 =", sprintf("%.3f", frac_sparse), ").\n")
                return(FALSE)
            }
        }
        x = as.matrix(tbl[, -1, drop = FALSE])
        for (j in seq_len(ncol(x))) {
            vals = suppressWarnings(as.numeric(x[, j]))
            fill = mean(vals, na.rm = TRUE)
            if (!is.finite(fill)) {
                fill = 0
            }
            vals[is.na(vals)] = fill
            x[, j] = vals
        }
        x = scale(x, center = TRUE, scale = TRUE)
        x[!is.finite(x)] = 0
        d = try(svd(x, nu = 0, nv = 0)$d, silent = TRUE)
        if (!inherits(d, "try-error")) {
            d = d[is.finite(d) & (d > 0)]
            if (length(d) >= 2) {
                kappa_est = max(d) / min(d)
                if (is.finite(kappa_est) && (kappa_est > 1e8)) {
                    cat("Skipping correlated-imputation attempts: ill-conditioned trait matrix (kappa~", format(kappa_est, scientific = TRUE), ").\n", sep = "")
                    return(FALSE)
                }
            }
        }
        return(TRUE)
    }

    correlated_ok = should_try_correlated(trait_table2, tree)
    attempts = list()
    if (correlated_ok) {
        attempts[[length(attempts) + 1L]] = list(label = "phylo_correlated=TRUE, pheno_correlated=TRUE", phylo = TRUE, pheno = TRUE)
        attempts[[length(attempts) + 1L]] = list(label = "phylo_correlated=TRUE, pheno_correlated=FALSE", phylo = TRUE, pheno = FALSE)
    } else {
        cat("Proceeding directly with uncorrelated phylogenetic-imputation fallback.\n")
    }
    attempts[[length(attempts) + 1L]] = list(label = "phylo_correlated=FALSE, pheno_correlated=FALSE", phylo = FALSE, pheno = FALSE)
    for (attempt in attempts) {
        cat("Phylogenetic imputation attempt:", attempt$label, "\n")
        rp_out = try(
            rkftools:::.rphylopars_phylopars(
                tree = tree,
                trait_data = trait_table2,
                phylo_correlated = attempt$phylo,
                pheno_correlated = attempt$pheno
            ),
            silent = TRUE
        )
        if (!inherits(rp_out, "try-error")) {
            imputed_matrix = data.frame(rp_out[["anc_recon"]], check.names = FALSE)
            imputed_matrix = imputed_matrix[tree[["tip.label"]], , drop = FALSE]
            cat("Phylogenetic imputation succeeded with", attempt$label, "\n")
            return(imputed_matrix)
        }
        cat("Phylogenetic imputation failed:", attempt$label, "::", conditionMessage(attr(rp_out, "condition")), "\n")
    }

    cat("All phylogenetic-imputation attempts failed. Falling back to column-mean imputation.\n")
    imputed_matrix = trait_table[tree$tip.label, , drop = FALSE]
    for (j in 1:ncol(imputed_matrix)) {
        vals = suppressWarnings(as.numeric(imputed_matrix[[j]]))
        fill = mean(vals, na.rm = TRUE)
        if (!is.finite(fill)) {
            fill = 0
        }
        vals[is.na(vals)] = fill
        imputed_matrix[[j]] = vals
    }
    return(imputed_matrix)
}

### Phylogenetic imputation
trait_table = impute_traits_with_fallback(tree, trait_table)

if (args[['clade_collapse_similarity_threshold']]) {
    cat('Start cluster collapsing by similarity thresholding (cutoff =', args[['clade_collapse_similarity_threshold']], ')\n')
    high_comple_clades = get_high_similarity_clades(tree, trait_table, args[['clade_collapse_similarity_method']], 
                            args[['clade_collapse_similarity_threshold']], verbose=FALSE, num_test=0)
    out = collapse_clades(tree, trait_table, high_comple_clades)
    tree = out[['tree']]
    tree = ensure_internal_node_names(tree)
    trait_table = out[['trait']]
} else {
    cat('Skipped cluster collapsing.', '\n')
}

### pad_short_branch() practically let this chunk obsoleted
tree_bifurcation = tree
tree = collapse_short_branches(tree, tol=1e-8)
tree = ensure_internal_node_names(tree)

trait_table = trait_table[tree$tip.label,]
trait_matrix = t(as.matrix(trait_table))
trait_matrix = trait_matrix[,tree$tip.label]
if (nrow(trait_matrix) >= (length(tree$tip.label) - 1)) {
    write_placeholder_outputs_and_quit(
        tree = tree,
        trait_table = trait_table,
        trait_matrix = trait_matrix,
        reason = sprintf(
            'Too many traits for stable PhyloEM optimization (num_traits=%d, num_tips=%d). Falling back to placeholder outputs.',
            nrow(trait_matrix), length(tree$tip.label)
        )
    )
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
    retry_out = run_phyloem_with_retries(tree = tree, trait_matrix = trait_matrix, nslots = args[['nslots']])
    res = retry_out[['res']]
    trait_matrix = retry_out[['trait_matrix']]
    if (inherits(res, "try-error")) {
        write_placeholder_outputs_and_quit(
            tree = tree,
            trait_table = trait_table,
            trait_matrix = trait_matrix,
            reason = paste("PhyloEM failed:", safe_condition_message(res), "Generating placeholder outputs and exiting successfully.")
        )
    }
}
save(res, file='PhylogeneticEM.out.RData')

cat('Time elapsed for PhyloEM:\n')
print(proc.time()-my_time)

my_time = proc.time()

tab_out = safely_get_phyloem_tables(res, tree, trait_table)
tree_table = tab_out[['tree_table']]
regime_table = tab_out[['regime_table']]
leaf_table = tab_out[['leaf_table']]

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

pdf("PhylogeneticEM.plot.pdf", height=max(2, length(tree$tip.label)/10+1), width=max(8, nrow(trait_matrix)+5))
par(mar=c(2,2,2,2))
plot_try = try(plot(res, show.tip.label=TRUE, label_cex=0.5), silent = TRUE)
if (inherits(plot_try, "try-error")) {
    cat("plot(PhylogeneticEM result) failed:", safe_condition_message(plot_try), "\n")
    try(plot(tree, show.tip.label=TRUE, cex=0.5), silent = TRUE)
}
graphics.off()
cat('PhylogeneticEM completed!\n')
