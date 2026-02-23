cli_args = commandArgs(trailingOnly=TRUE)
args = cli_args

library(l1ou, quietly=TRUE)
library(rkftools, quietly=TRUE)
options(stringsAsFactors=FALSE)

# one bootstrap replicate needs approx. 1/2 time of the original fit estimation at nCores=1.
# OUfixedRoot is preferrable, https://github.com/khabbazian/l1ou/pull/5
# fixed.alpha =FALSE might take too long time.

cat('arguments:\n')
args = rkftools::get_parsed_args(args, print=TRUE)

# this function was copied from the package PhylogeneticEM
find_grid_alpha <- function(phy, alpha = NULL,
                            nbr_alpha = 10,
                            factor_up_alpha = 2,
                            factor_down_alpha = 3,
                            quantile_low_distance = 0.0001,
                            log_transform = TRUE, ...){
  if (!is.null(alpha)) return(alpha)
  dtips <- cophenetic(phy)
  d_min <- quantile(dtips[dtips > 0], quantile_low_distance)
  h_tree <- node.depth.edgelength(phy)[1]
  alpha_min <- 1 / (factor_down_alpha * h_tree)
  alpha_max <- factor_up_alpha / (2 * d_min)
  alpha_max_machine <- log(.Machine$double.xmax^0.98)/(2*h_tree)
  if (alpha_max > alpha_max_machine){
    warning("The chosen alpha_max was above the machine precision. Taking alpha_max as the largest possible on this machine.")
    alpha_max <- alpha_max_machine
  }
  if (log_transform){
    return(c(0, exp(seq(log(alpha_min), log(alpha_max), length.out = nbr_alpha))))
  } else {
    return(c(0, seq(alpha_min, alpha_max, length.out = nbr_alpha)))
  }
}

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

node_num_to_name = function(tree, node_num) {
    if (!is.finite(node_num)) {
        return(NA_character_)
    }
    num_tip = length(tree$tip.label)
    if ((node_num >= 1) && (node_num <= num_tip)) {
        return(as.character(tree$tip.label[[node_num]]))
    }
    node_index = node_num - num_tip
    if (is.null(tree$node.label)) {
        return(paste0("n", node_num))
    }
    if ((node_index >= 1) && (node_index <= length(tree$node.label))) {
        out = as.character(tree$node.label[[node_index]])
        if (!is.na(out) && (out != "")) {
            return(out)
        }
    }
    return(paste0("n", node_num))
}

extract_shift_param_values = function(param_obj, shift_idx, trait_names) {
    out = rep(NA_real_, length(trait_names))
    names(out) = trait_names
    if (is.null(param_obj) || (length(trait_names) == 0)) {
        return(out)
    }
    mat = try(as.matrix(param_obj), silent = TRUE)
    if (inherits(mat, "try-error")) {
        mat = NULL
    }
    vals = NULL
    if (!is.null(mat)) {
        if ((nrow(mat) == length(trait_names)) && (ncol(mat) >= shift_idx)) {
            vals = mat[, shift_idx]
        } else if ((ncol(mat) == length(trait_names)) && (nrow(mat) >= shift_idx)) {
            vals = mat[shift_idx, ]
        }
    }
    if (is.null(vals) && (length(trait_names) == 1) && (length(param_obj) >= shift_idx)) {
        vals = unname(param_obj[[shift_idx]])
    }
    if (is.null(vals)) {
        return(out)
    }
    vals = suppressWarnings(as.numeric(vals))
    if (length(vals) == length(trait_names)) {
        out = vals
    }
    return(out)
}

build_l1ou_regime_table_fallback = function(fit_obj, tree, trait_table, leaf_table) {
    regime_table = get_placeholder_regime(tree, trait_table)
    shift_cfg = fit_obj[['shift.configuration']]
    if (is.null(shift_cfg) || (length(shift_cfg) == 0)) {
        return(regime_table)
    }
    shift_nodes = suppressWarnings(as.integer(unname(shift_cfg)))
    shift_regimes = suppressWarnings(as.integer(names(shift_cfg)))
    if (all(is.na(shift_regimes))) {
        shift_regimes = seq_along(shift_nodes)
    }
    if (nrow(leaf_table) > 0 && ("regime" %in% colnames(leaf_table))) {
        valid_regimes = sort(unique(suppressWarnings(as.integer(leaf_table$regime))))
        valid_regimes = valid_regimes[!is.na(valid_regimes)]
        if (length(valid_regimes) > 0) {
            for (i in seq_along(shift_regimes)) {
                if (is.na(shift_regimes[[i]]) || (!(shift_regimes[[i]] %in% valid_regimes))) {
                    shift_regimes[[i]] = valid_regimes[[which.max(valid_regimes)]]
                }
            }
        }
    }
    trait_names = colnames(trait_table)
    shift_rows = list()
    for (i in seq_along(shift_nodes)) {
        node_name = node_num_to_name(tree, shift_nodes[[i]])
        if (is.na(node_name) || (node_name == "")) {
            next
        }
        val_row = as.list(rep(NA_real_, length(trait_names)))
        names(val_row) = trait_names
        val_row = c(
            list(regime = shift_regimes[[i]], node_name = node_name, param = "shift_value"),
            extract_shift_param_values(fit_obj[['shift.values']], i, trait_names)
        )
        mean_row = c(
            list(regime = shift_regimes[[i]], node_name = node_name, param = "shift_mean"),
            extract_shift_param_values(fit_obj[['shift.means']], i, trait_names)
        )
        shift_rows[[length(shift_rows) + 1]] = as.data.frame(val_row, check.names = FALSE)
        shift_rows[[length(shift_rows) + 1]] = as.data.frame(mean_row, check.names = FALSE)
    }
    if (length(shift_rows) == 0) {
        return(regime_table)
    }
    shift_df = do.call(rbind, shift_rows)
    for (col in colnames(regime_table)) {
        if (!(col %in% colnames(shift_df))) {
            shift_df[[col]] = NA
        }
    }
    shift_df = shift_df[, colnames(regime_table), drop = FALSE]
    out = rbind(shift_df, regime_table)
    rownames(out) = NULL
    return(out)
}

safely_get_l1ou_tables = function(fit_obj, tree, trait_table) {
    tree_try = try(get_tree_table(fit_obj, mode = 'l1ou'), silent = TRUE)
    regime_try = try(get_regime_table(fit_obj, mode = 'l1ou'), silent = TRUE)
    leaf_try = try(get_leaf_table(fit_obj, mode = 'l1ou'), silent = TRUE)

    if (inherits(tree_try, "try-error")) {
        cat("get_tree_table(l1ou) failed:", safe_condition_message(tree_try), "\n")
        tree_table = get_placeholder_tree(tree, trait_table)
    } else {
        tree_table = tree_try
    }
    if (inherits(leaf_try, "try-error")) {
        cat("get_leaf_table(l1ou) failed:", safe_condition_message(leaf_try), "\n")
        leaf_table = get_placeholder_leaf(tree, trait_table)
    } else {
        leaf_table = leaf_try
    }
    if (inherits(regime_try, "try-error")) {
        cat("get_regime_table(l1ou) failed:", safe_condition_message(regime_try), "\n")
        regime_table = build_l1ou_regime_table_fallback(fit_obj, tree, trait_table, leaf_table)
    } else {
        regime_table = regime_try
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

my_time = proc.time()
tree = read.tree(args[['tree_file']])
require_internal_node_labels = parse_bool_flag(args[['require_internal_node_labels']], default = 0L)
cat('require_internal_node_labels =', require_internal_node_labels, '\n')
if (require_internal_node_labels == 1L) {
    validate_internal_node_labels(tree, context = args[['tree_file']])
}
tree = ensure_internal_node_names(tree)
if (!isTRUE(ape::is.binary(tree))) {
    stop("Input tree is not binary. Resolve tree polytomies upstream before running detect_OU_shift_l1ou.r.")
}
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

write_l1ou_placeholder_outputs = function(tree, trait_table, fit_ind_obj = NA, fit_conv_obj = NA, reason = "", quit_after = FALSE) {
    if (nzchar(reason)) {
        cat(reason, "\n")
    }
    leaf_table = get_placeholder_leaf(tree, trait_table)
    regime_table = get_placeholder_regime(tree, trait_table)
    tree_table = get_placeholder_tree(tree, trait_table)
    write.table(tree_table, file="l1ou_tree.tsv", sep="\t", quote=FALSE, row.names=FALSE) 
    write.table(regime_table, file="l1ou_regime.tsv", sep="\t", quote=FALSE, row.names=FALSE) 
    write.table(leaf_table, file="l1ou_leaf.tsv", sep="\t", quote=FALSE, row.names=FALSE)
    pdf("l1ou_plot.pdf", height=length(tree$tip.label)/5+1, width=8)
    plot(tree, show.tip.label=TRUE, cex=0.5)
    graphics.off()
    fit_ind = fit_ind_obj
    fit_conv = fit_conv_obj
    save(fit_ind, file='fit_ind.RData')
    save(fit_conv, file='fit_conv.RData')
    if (quit_after) {
        cat('Placeholder outputs were generated. Exiting.\n')
        quit(save='no', status=0)
    }
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

out = remove_invariant_traits(trait_table, small_dif=1e-3)
trait_table = out[['trait_table']]
removed_traits = out[['removed_traits']]
if (ncol(trait_table)==0) {
    write_l1ou_placeholder_outputs(
        tree = tree,
        trait_table = original_trait_table,
        reason = 'All traits are invariant. Generating placeholder outputs without running l1ou.',
        quit_after = TRUE
    )
    stop('No available trait. This line may be printed in debugging, but should not appear in the batch mode.\n')
}

### Phylogenetic imputation
trait_table = impute_traits_with_fallback(tree, trait_table)

if (args[['clade_collapse_similarity_threshold']]) {
    cat('Start cluster collapsing by similarity thresholding ( cutoff =', args[['clade_collapse_similarity_threshold']], ')\n')
    high_comple_clades = get_high_similarity_clades(tree, trait_table, args[['clade_collapse_similarity_method']], 
                            args[['clade_collapse_similarity_threshold']], verbose=FALSE, num_test=0)
    out = collapse_clades(tree, trait_table, high_comple_clades)
    if (nrow(out[['trait']])<4) {
        cat('Only', nrow(out[['trait']]), 'genes remained. Clade collapsing is skipped.\n')
    } else {
        tree = out[['tree']]
        tree = ensure_internal_node_names(tree)
        trait_table = out[['trait']]
    }
} else {
    cat('Skipped cluster collapsing.', '\n')
}

trait_matrix = as.matrix(trait_table)
if (ncol(trait_matrix) >= (length(tree$tip.label) - 1)) {
    write_l1ou_placeholder_outputs(
        tree = tree,
        trait_table = trait_table,
        reason = sprintf(
            'Too many traits for stable l1ou optimization (num_traits=%d, num_tips=%d). Falling back to placeholder outputs.',
            ncol(trait_matrix), length(tree$tip.label)
        ),
        quit_after = TRUE
    )
}

cat("Time elapsed for data preparation:\n")
print(proc.time()-my_time)

my_time = proc.time()
fit_ind=NULL
fit_conv=NULL
adj_data=NULL
adj_data = l1ou::adjust_data(tree=tree, Y=trait_matrix, normalize=FALSE)

alpha.lower = NA
if (args[['alpha_upper']]=='Inf') {
    alpha.upper = Inf
} else if (args[['alpha_upper']]=='l1ou') {
    alpha.upper = l1ou:::alpha_upper_bound(adj_data$tree)
} else if (args[['alpha_upper']]=='PhylogeneticEM') {
    alphas = find_grid_alpha(phy=tree, nbr_alpha=10)
    alpha.upper = alphas[11]
    # alpha.lower = alphas[1] # l1ou fails if alpha.lower is set to 0.
} else {
    alpha.upper = as.numeric(args[['alpha_upper']])
}
cat('default alpha upper bound in l1ou =',  l1ou:::alpha_upper_bound(adj_data[['tree']]), '\n')
cat('the alpha upper bound actually used =', alpha.upper, '\n')

args[['max_nshift']] = as.integer(args[['max_nshift']])
num_leaf = length(tree[['tip.label']])
if (args[['max_nshift']]==0) {
    max_nshift = as.integer(max(min(100, num_leaf/2), floor(num_leaf^0.5)))
} else {
    max_nshift = args[['max_nshift']]
}
cat('maximum number of shifts =', max_nshift, '\n')

ind_flag = FALSE
if ((!is.null(args[['fit_ind_file']]))&(!is.na(args[['fit_ind_file']]))) {
    if (file.exists(args[['fit_ind_file']])) {
        if (file.info(args[['fit_ind_file']])$size) {
            load(args[['fit_ind_file']])
            if(class(fit_ind)!="try-error") {
                cat(args[['fit_ind_file']], 'was found. Loading the shift configurations.\n')
                ind_flag = TRUE
            } else {
                cat(args[['fit_ind_file']], 'is malformed. Printing...\n')
                print(fit_ind)
                rm(fit_ind)
            }
        } else {
            cat(args[['fit_ind_file']], 'has size zero.\n')
        }
    }
}
if (!ind_flag) {
    cat('fit_ind.RData was not found. Starting shift detection.\n')
    fit_ind = try(l1ou::estimate_shift_configuration(
        tree=adj_data[['tree']], 
        Y=adj_data[['Y']], 
        max.nShifts=max_nshift,
        criterion=args[['criterion']], 
        root.model="OUfixedRoot", 
        nCores=args[['nslots']], 
        rescale=FALSE, 
        alpha.upper=alpha.upper, 
        alpha.lower=alpha.lower,
        quietly=FALSE))
    save(fit_ind, file='fit_ind.RData')
    cat("Time elapsed for shift detection:\n")
}
print(proc.time()-my_time)

conv_flag = FALSE
if ((!is.null(args[['fit_conv_file']]))&(!is.na(args[['fit_conv_file']]))) {
    if (file.exists(args[['fit_conv_file']])) {
        if (file.info(args[['fit_conv_file']])[['size']]) {
            load(args[['fit_conv_file']])
            if(class(fit_conv)!="try-error") {
                cat(args[['fit_conv_file']], 'was found. Loading the shift configurations.\n')
                conv_flag = TRUE
            } else {
                cat(args[['fit_conv_file']], 'is malformed. Printing...\n')
                print(fit_conv)
                rm(fit_conv)
            }
        } else {
            cat(args[['fit_conv_file']], 'has size zero.\n')
        }
    }
}
if ((!conv_flag)&(args[['detect_convergence']])) {
    cat('Started convergence analysis.\n')
    if (class(fit_ind)=="l1ou") {
        if (fit_ind[['nShifts']] < 2) {
            fit_conv = fit_ind
        } else {
            fit_conv = try(l1ou::estimate_convergent_regimes(
                fit_ind, 
                criterion=args[['criterion']], 
                method="backward", 
                fixed.alpha=TRUE, 
                nCores=args[['nslots']]))
        }
    } else {
        warning("estimate_shift_configuration() failed.")
        warning(attr(fit_ind, "condition"))
    }
    save(fit_conv, file='fit_conv.RData')
    cat("Time elapsed for shift detection + convergence detection:\n")
    print(fit_conv)
    print(proc.time()-my_time)
} else {
    cat('Skipped convergence analysis.\n')
    fit_conv = fit_ind
}


if (class(fit_conv)=="l1ou") {
    if (args[['nbootstrap']] > 0) {
        cat("Start bootstrapping,", args[['nbootstrap']], "replicates,", args[['nslots']], "cores.\n")
        bootstrap_result = l1ou::l1ou_bootstrap_support(
            fit_conv, 
            nItrs=args[['nbootstrap']], 
            multicore=(args[['nslots']]>1), 
            nCores=args[['nslots']])
        bootstrap_table = get_bootstrap_table(fit_conv, bootstrap_result, mode='l1ou')
        write.table(bootstrap_table, file="l1ou_bootstrap.tsv", sep="\t", quote=FALSE, row.names=FALSE)
        cat("Time elapsed for shift detection + convergence detection + bootstrapping:\n")
        print(proc.time()-my_time)
    } 
    cat('Start preparing output tables.\n')
    tab_out = safely_get_l1ou_tables(fit_conv, tree, trait_table)
    tree_table = tab_out[['tree_table']]
    regime_table = tab_out[['regime_table']]
    leaf_table = tab_out[['leaf_table']]
    for (mt in removed_traits) {
        regime_table[mt] = NA
        leaf_table[mt] = 0 
    }
    for (mt in original_traits[!(original_traits %in% colnames(regime_table))]) {
        regime_table[mt] = NA
    }
    for (mt in original_traits[!(original_traits %in% colnames(leaf_table))]) {
        leaf_table[mt] = 0
    }
    regime_table = regime_table[,c(colnames(regime_table)[!colnames(regime_table) %in% original_traits], original_traits)]
    leaf_table = leaf_table[,c(colnames(leaf_table)[!colnames(leaf_table) %in% original_traits], original_traits)]
    
    if (args[['clade_collapse_similarity_threshold']]) {
        cat('Starting the recovery of collapsed nodes.\n')
        node_num_mapping = map_node_num(tree_original, tree, out[['collapse_leaf_names']], verbose=FALSE)
        tree_table = tree_table_collapse2original(tree_table, tree_original)
        regime_table = regime_table_collapse2original(regime_table, tree_original, tree, node_num_mapping)
        leaf_table = leaf_table_collapse2original(leaf_table, tree_original, tree, node_num_mapping)
    }
    
    write.table(tree_table, file="l1ou_tree.tsv", sep="\t", quote=FALSE, row.names=FALSE) 
    write.table(regime_table, file="l1ou_regime.tsv", sep="\t", quote=FALSE, row.names=FALSE) 
    write.table(leaf_table, file="l1ou_leaf.tsv", sep="\t", quote=FALSE, row.names=FALSE)
    pdf("l1ou_plot.pdf", height=length(tree$tip.label)/5+1, width=ncol(adj_data$Y)/2+5)
    par(mar=c(2,2,2,2))
    try(plot(fit_conv, edge.shift.ann=FALSE))
    graphics.off()
} else {
    warning("estimate_convergent_regimes() failed. Writing placeholder outputs.")
    if (inherits(fit_conv, "try-error")) {
        warning(attr(fit_conv, "condition"))
    }
    write_l1ou_placeholder_outputs(
        tree = tree,
        trait_table = trait_table,
        fit_ind_obj = fit_ind,
        fit_conv_obj = fit_conv,
        reason = 'l1ou fit did not converge to a valid object. Placeholder outputs were written.',
        quit_after = FALSE
    )
}

cat('l1ou completed!\n')
