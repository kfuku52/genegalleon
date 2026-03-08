cli_args = commandArgs(trailingOnly = TRUE)
args = cli_args

library(kfl1ou, quietly = TRUE)
options(stringsAsFactors = FALSE)

load_rkftools_helpers = function() {
    if (requireNamespace("rkftools", quietly = TRUE)) {
        library(rkftools, quietly = TRUE)
        return(invisible("package"))
    }

    script_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (!length(script_arg)) {
        stop("Package 'rkftools' is not installed, and the current script path could not be resolved for fallback loading.")
    }

    script_path = sub("^--file=", "", script_arg[[1]])
    script_path = gsub("~\\+~", " ", script_path, fixed = FALSE)
    script_path = normalizePath(script_path, mustWork = TRUE)
    repos_dir = normalizePath(file.path(dirname(script_path), "..", "..", ".."), mustWork = FALSE)
    rkftools_dir = file.path(repos_dir, "rkftools")
    source_files = c(
        file.path(rkftools_dir, "R", "kfutil.R"),
        file.path(rkftools_dir, "R", "kfphylo.R"),
        file.path(rkftools_dir, "R", "kfpcm.R")
    )

    missing_files = source_files[!file.exists(source_files)]
    if (length(missing_files)) {
        stop(
            "Package 'rkftools' is not installed, and fallback source files were not found: ",
            paste(missing_files, collapse = ", ")
        )
    }

    for (source_file in source_files) {
        source(source_file, local = .GlobalEnv)
    }
    cat("Package 'rkftools' is not installed; sourced helper functions from ", rkftools_dir, ".\n", sep = "")
    invisible("source")
}

load_rkftools_helpers()

cat("arguments:\n")
args = get_parsed_args(args, print = TRUE)

validate_args = function(args, allowed, required) {
    unknown = setdiff(names(args), allowed)
    if (length(unknown)) {
        stop(
            "Unknown argument(s): ",
            paste(unknown, collapse = ", "),
            ". Supported arguments are: ",
            paste(sort(allowed), collapse = ", "),
            "."
        )
    }

    missing_required = required[vapply(required, function(x) {
        is.null(args[[x]]) || (length(args[[x]]) == 0L) || is.na(args[[x]]) || (args[[x]] == "")
    }, logical(1))]
    if (length(missing_required)) {
        stop("Missing required argument(s): ", paste(missing_required, collapse = ", "), ".")
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
    if (is.null(x) || (length(x) == 0L) || is.na(x) || (x == "")) {
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

parse_integer_arg = function(x, default = NA_integer_) {
    if (is.null(x) || (length(x) == 0L) || is.na(x) || (x == "")) {
        return(as.integer(default))
    }
    out = suppressWarnings(as.integer(x[[1]]))
    if (is.na(out)) {
        return(as.integer(default))
    }
    return(out)
}

parse_string_arg = function(x, default = NULL) {
    if (is.null(x) || (length(x) == 0L) || is.na(x) || (x == "")) {
        return(default)
    }
    return(as.character(x[[1]]))
}

normalize_replicate_sep = function(x, default = "_") {
    out = parse_string_arg(x, default = default)
    if (is.null(out) || is.na(out)) {
        return(default)
    }
    return(out)
}

merge_replicates_with_input_error = function(raw_trait_table, replicate_sep = "_") {
    trait_table = as.data.frame(raw_trait_table, check.names = FALSE, stringsAsFactors = FALSE)
    if (is.null(rownames(trait_table))) {
        stop("trait_table must have row names before replicate processing.")
    }
    for (j in seq_len(ncol(trait_table))) {
        trait_table[[j]] = suppressWarnings(as.numeric(trait_table[[j]]))
    }

    merged_trait_table = merge_replicates(trait_table, replicate_sep)
    input_error = matrix(
        0,
        nrow = nrow(merged_trait_table),
        ncol = ncol(merged_trait_table),
        dimnames = list(rownames(merged_trait_table), colnames(merged_trait_table))
    )
    if (ncol(trait_table) == 0L) {
        return(list(trait_table = merged_trait_table, input_error = input_error))
    }

    original_cols = colnames(trait_table)
    for (trait_name in colnames(merged_trait_table)) {
        is_col = (original_cols == trait_name)
        if (!(isTRUE(replicate_sep == "") || is.na(replicate_sep))) {
            is_col = is_col | startsWith(original_cols, paste0(trait_name, replicate_sep))
        }
        values = as.matrix(trait_table[, is_col, drop = FALSE])
        if (!length(values)) {
            next
        }
        se2 = apply(values, 1, function(v) {
            vv = suppressWarnings(as.numeric(v))
            vv = vv[is.finite(vv)]
            n_obs = length(vv)
            if (n_obs <= 1L) {
                return(0)
            }
            trait_var = suppressWarnings(stats::var(vv))
            if (!is.finite(trait_var) || (trait_var < 0)) {
                return(0)
            }
            return(trait_var / n_obs)
        })
        se2[!is.finite(se2)] = 0
        se2[se2 < 0] = 0
        input_error[, trait_name] = se2
    }

    num_groups_with_replicates = sum(vapply(
        colnames(merged_trait_table),
        function(trait_name) {
            is_col = (original_cols == trait_name)
            if (!(isTRUE(replicate_sep == "") || is.na(replicate_sep))) {
                is_col = is_col | startsWith(original_cols, paste0(trait_name, replicate_sep))
            }
            sum(is_col) > 1L
        },
        logical(1)
    ))
    num_traits_with_nonzero_error = sum(colSums(input_error > 0, na.rm = TRUE) > 0L)
    cat(
        "Replicate preprocessing:",
        ncol(trait_table), "input columns ->",
        ncol(merged_trait_table), "merged traits.",
        num_groups_with_replicates, "trait groups had replicates and",
        num_traits_with_nonzero_error, "traits received non-zero tip-specific input_error.\n"
    )
    return(list(trait_table = merged_trait_table, input_error = input_error))
}

align_input_error_to_tree = function(tree, input_error, trait_names = NULL, context = "input_error") {
    if (is.null(input_error)) {
        return(NULL)
    }
    aligned = as.data.frame(input_error, check.names = FALSE, stringsAsFactors = FALSE)
    aligned = align_trait_table_to_tree(tree, aligned, context = context)
    if (!is.null(trait_names)) {
        missing_traits = setdiff(trait_names, colnames(aligned))
        extra_traits = setdiff(colnames(aligned), trait_names)
        if (length(missing_traits) || length(extra_traits)) {
            stop(
                context, " column names must match the trait names. Missing: ",
                ifelse(length(missing_traits), paste(missing_traits, collapse = ", "), "none"),
                ". Extra: ",
                ifelse(length(extra_traits), paste(extra_traits, collapse = ", "), "none"),
                "."
            )
        }
        aligned = aligned[, trait_names, drop = FALSE]
    }
    aligned = as.matrix(aligned)
    aligned[!is.finite(aligned)] = 0
    aligned[aligned < 0] = 0
    return(aligned)
}

fit_matches_current_data = function(fit_obj, tree, trait_matrix) {
    if (!inherits(fit_obj, "l1ou")) {
        return(FALSE)
    }
    fit_tree = fit_obj[["tree"]]
    fit_Y = fit_obj[["Y"]]
    if (is.null(fit_tree) || is.null(fit_Y) || !inherits(fit_tree, "phylo")) {
        return(FALSE)
    }
    fit_Y = as.matrix(fit_Y)
    trait_matrix = as.matrix(trait_matrix)
    if (!identical(fit_tree$tip.label, tree$tip.label)) {
        return(FALSE)
    }
    if (!identical(dim(fit_Y), dim(trait_matrix))) {
        return(FALSE)
    }
    if (!identical(rownames(fit_Y), rownames(trait_matrix))) {
        return(FALSE)
    }
    if (!identical(colnames(fit_Y), colnames(trait_matrix))) {
        return(FALSE)
    }
    return(TRUE)
}

align_trait_table_to_tree = function(tree, trait_table, context = "trait_table") {
    if (is.null(rownames(trait_table))) {
        stop(context, " must have row names matching tree tip labels.")
    }
    missing_rows = setdiff(tree$tip.label, rownames(trait_table))
    extra_rows = setdiff(rownames(trait_table), tree$tip.label)
    if (length(missing_rows) || length(extra_rows)) {
        stop(
            context, " row names must match tree tip labels. Missing: ",
            ifelse(length(missing_rows), paste(missing_rows, collapse = ", "), "none"),
            ". Extra: ",
            ifelse(length(extra_rows), paste(extra_rows, collapse = ", "), "none"),
            "."
        )
    }
    trait_table[tree$tip.label, , drop = FALSE]
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
    if (is.null(param_obj) || (length(trait_names) == 0L)) {
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
    if (is.null(vals) && (length(trait_names) == 1L) && (length(param_obj) >= shift_idx)) {
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
    shift_cfg = fit_obj[["shift.configuration"]]
    if (is.null(shift_cfg) || (length(shift_cfg) == 0L)) {
        return(regime_table)
    }

    shift_edges = suppressWarnings(as.integer(unname(shift_cfg)))
    shift_regimes = suppressWarnings(as.integer(names(shift_cfg)))
    if (all(is.na(shift_regimes))) {
        shift_regimes = seq_along(shift_edges)
    }

    if (nrow(leaf_table) > 0L && ("regime" %in% colnames(leaf_table))) {
        valid_regimes = sort(unique(suppressWarnings(as.integer(leaf_table$regime))))
        valid_regimes = valid_regimes[!is.na(valid_regimes)]
        if (length(valid_regimes) > 0L) {
            for (i in seq_along(shift_regimes)) {
                if (is.na(shift_regimes[[i]]) || (!(shift_regimes[[i]] %in% valid_regimes))) {
                    shift_regimes[[i]] = valid_regimes[[which.max(valid_regimes)]]
                }
            }
        }
    }

    trait_names = colnames(trait_table)
    shift_rows = list()
    for (i in seq_along(shift_edges)) {
        edge_idx = shift_edges[[i]]
        if (!is.finite(edge_idx) || (edge_idx < 1L) || (edge_idx > nrow(tree$edge))) {
            next
        }
        node_name = node_num_to_name(tree, tree$edge[edge_idx, 2])
        if (is.na(node_name) || (node_name == "")) {
            next
        }
        val_row = c(
            list(regime = shift_regimes[[i]], node_name = node_name, param = "shift_value"),
            as.list(extract_shift_param_values(fit_obj[["shift.values"]], i, trait_names))
        )
        mean_row = c(
            list(regime = shift_regimes[[i]], node_name = node_name, param = "shift_mean"),
            as.list(extract_shift_param_values(fit_obj[["shift.means"]], i, trait_names))
        )
        shift_rows[[length(shift_rows) + 1L]] = as.data.frame(val_row, check.names = FALSE)
        shift_rows[[length(shift_rows) + 1L]] = as.data.frame(mean_row, check.names = FALSE)
    }

    if (length(shift_rows) == 0L) {
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
    tree_try = try(get_tree_table(fit_obj, mode = "l1ou"), silent = TRUE)
    regime_try = try(get_regime_table(fit_obj, mode = "l1ou"), silent = TRUE)
    leaf_try = try(get_leaf_table(fit_obj, mode = "l1ou"), silent = TRUE)

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

    return(list(
        tree_table = tree_table,
        regime_table = regime_table,
        leaf_table = leaf_table
    ))
}

save_named_object = function(name, value, file) {
    tmp = new.env(parent = emptyenv())
    assign(name, value, envir = tmp)
    save(list = name, file = file, envir = tmp)
}

write_placeholder_outputs = function(tree, trait_table,
                                     fit_ind_obj = NA,
                                     fit_conv_obj = NA,
                                     reason = "",
                                     quit_after = FALSE) {
    if (nzchar(reason)) {
        cat(reason, "\n")
    }
    kfl1ou::write_l1ou_placeholder_outputs(
        tree = tree,
        Y = trait_table,
        dir = ".",
        prefix = "l1ou",
        reason = reason
    )
    save_named_object("fit_ind", fit_ind_obj, "fit_ind.RData")
    save_named_object("fit_conv", fit_conv_obj, "fit_conv.RData")
    if (quit_after) {
        cat("Placeholder outputs were generated. Exiting.\n")
        quit(save = "no", status = 0)
    }
}

load_cached_fit = function(path, object_name) {
    if (is.null(path) || is.na(path) || (path == "")) {
        return(NULL)
    }
    if (!file.exists(path)) {
        return(NULL)
    }
    if (!isTRUE(file.info(path)$size > 0)) {
        cat(path, "has size zero.\n")
        return(NULL)
    }
    tmp = new.env(parent = emptyenv())
    loaded = try(load(path, envir = tmp), silent = TRUE)
    if (inherits(loaded, "try-error")) {
        cat(path, "could not be loaded:", safe_condition_message(loaded), "\n")
        return(NULL)
    }
    if (object_name %in% loaded) {
        return(tmp[[object_name]])
    }
    if (length(loaded) == 1L) {
        return(tmp[[loaded[[1]]]])
    }
    cat(path, "does not contain", object_name, ".\n")
    return(NULL)
}

resolve_alpha_upper = function(arg, tree) {
    default_upper = kfl1ou:::alpha_upper_bound(tree)
    sval = parse_string_arg(arg, default = "kfl1ou")
    if (is.null(sval) || (sval %in% c("kfl1ou", "l1ou", "auto"))) {
        return(list(value = default_upper, default = default_upper))
    }
    if (sval == "Inf") {
        return(list(value = Inf, default = default_upper))
    }
    numeric_upper = suppressWarnings(as.numeric(sval))
    if (is.finite(numeric_upper)) {
        return(list(value = numeric_upper, default = default_upper))
    }
    warning("Unrecognized alpha_upper value: ", sval, ". Falling back to the kfl1ou default upper bound.")
    return(list(value = default_upper, default = default_upper))
}

resolve_max_nshift = function(arg, tree) {
    auto_limit = floor(nrow(tree$edge) / 2)
    sval = parse_string_arg(arg, default = NULL)
    if (is.null(sval) || (sval %in% c("0", "auto", "AUTO", "null", "NULL", "NA", "na"))) {
        return(list(value = NULL, reported = auto_limit))
    }
    numeric_limit = suppressWarnings(as.integer(sval))
    if (is.na(numeric_limit)) {
        warning("Unrecognized max_nshift value: ", sval, ". Falling back to auto.")
        return(list(value = NULL, reported = auto_limit))
    }
    return(list(value = numeric_limit, reported = numeric_limit))
}

restore_fit_for_output = function(model, original.tree) {
    if (!inherits(model, "l1ou")) {
        return(model)
    }
    if (length(model$tree$tip.label) == length(original.tree$tip.label) &&
        identical(sort(model$tree$tip.label), sort(original.tree$tip.label))) {
        return(model)
    }
    cat("Restoring fitted results to the repaired full tree.\n")
    return(kfl1ou::restore_original_tree_fit(model, original.tree))
}

restore_bootstrap_result_to_tree = function(model, bootstrap_result, original.tree) {
    if (!inherits(model, "l1ou")) {
        return(bootstrap_result)
    }
    if (length(model$tree$tip.label) == length(original.tree$tip.label) &&
        identical(sort(model$tree$tip.label), sort(original.tree$tip.label))) {
        return(bootstrap_result)
    }

    mapping = kfl1ou:::map_pruned_edges_to_original(
        pruned.tree = reorder(model$tree, "postorder"),
        original.tree = reorder(original.tree, "postorder"),
        representative = "tipward"
    )

    detection_rate = suppressWarnings(as.numeric(bootstrap_result$detection.rate))
    if (length(detection_rate) != length(mapping$representative.edge)) {
        stop(
            "Length mismatch while restoring bootstrap results: expected ",
            length(mapping$representative.edge), ", got ", length(detection_rate), "."
        )
    }

    restored = bootstrap_result
    restored$detection.rate = rep(0, length(mapping$original.to.pruned.edge))
    valid = !is.na(mapping$representative.edge)
    restored$detection.rate[mapping$representative.edge[valid]] = detection_rate[valid]

    if (!is.null(restored$all.shifts)) {
        restored$all.shifts = lapply(restored$all.shifts, function(shift_cfg) {
            shift_idx = suppressWarnings(as.integer(shift_cfg))
            shift_idx = shift_idx[!is.na(shift_idx)]
            mapping$representative.edge[shift_idx]
        })
    }

    return(restored)
}

tree_file = parse_string_arg(args[["tree_file"]], default = NULL)
trait_file = parse_string_arg(args[["trait_file"]], default = NULL)
validate_args(
    args = args,
    allowed = c(
        "alpha_upper",
        "criterion",
        "detect_convergence",
        "fit_conv_file",
        "fit_ind_file",
        "max_nshift",
        "nbootstrap",
        "nslots",
        "replicate_sep",
        "trait_file",
        "tree_file"
    ),
    required = c("tree_file", "trait_file")
)

tree_input = read.tree(tree_file)
if (!isTRUE(ape::is.binary(tree_input))) {
    stop("Input tree is not binary. Resolve tree polytomies upstream before running detect_OU_shift_kfl1ou.r.")
}

replicate_sep = normalize_replicate_sep(args[["replicate_sep"]], default = "_")
cat("replicate_sep =", replicate_sep, "\n")

raw_trait_table = read.table(trait_file, header = TRUE, row.names = 1, sep = "\t")
replicate_prep = merge_replicates_with_input_error(raw_trait_table, replicate_sep = replicate_sep)
trait_table = align_trait_table_to_tree(tree_input, replicate_prep$trait_table, context = "trait_table")
input_error = align_input_error_to_tree(
    tree_input,
    replicate_prep$input_error,
    trait_names = colnames(trait_table),
    context = "input_error"
)

prep_time = proc.time()
full_data = kfl1ou::adjust_data(
    tree = tree_input,
    Y = as.matrix(trait_table),
    normalize = FALSE,
    quietly = FALSE,
    repair.tree = TRUE,
    min.edge.length = 1e-3,
    ultrametric.tolerance = c(0.01, 0.05, 0.1, 1, Inf),
    drop.all.missing = FALSE,
    drop.invariant = FALSE
)

tree_full = full_data$tree
original_trait_table = as.data.frame(full_data$Y, check.names = FALSE, stringsAsFactors = FALSE)
original_traits = colnames(original_trait_table)
input_error_full = align_input_error_to_tree(
    tree_full,
    input_error,
    trait_names = original_traits,
    context = "input_error_full"
)

cat("Time elapsed for tree repair/alignment:\n")
print(proc.time() - prep_time)

adj_data = try(
    kfl1ou::adjust_data(
    tree = tree_full,
    Y = as.matrix(original_trait_table),
    normalize = FALSE,
    quietly = FALSE,
    repair.tree = FALSE,
    drop.all.missing = TRUE,
    drop.invariant = TRUE,
    invariant.tolerance = 1e-3
    ),
    silent = TRUE
)

if (inherits(adj_data, "try-error")) {
    adj_message = safe_condition_message(adj_data)
    if (grepl("all traits are invariant or entirely missing after adjustment", adj_message, fixed = TRUE)) {
        write_placeholder_outputs(
            tree = tree_full,
            trait_table = original_trait_table,
            reason = "All traits are invariant after kfl1ou preprocessing. Generating placeholder outputs without running kfl1ou.",
            quit_after = TRUE
        )
    }
    stop(adj_message)
}

removed_traits = adj_data$removed.traits
if (!length(removed_traits)) {
    removed_traits = character(0)
}

trait_matrix = as.matrix(adj_data$Y)
input_error_fit = align_input_error_to_tree(
    adj_data$tree,
    input_error_full,
    trait_names = colnames(trait_matrix),
    context = "input_error_fit"
)
if (!any(input_error_fit > 0, na.rm = TRUE)) {
    input_error_fit = NULL
    cat("No non-zero tip-specific input_error was derived from replicate variation. Shared measurement_error will still be estimated.\n")
} else {
    cat(
        "Using tip-specific input_error for",
        sum(colSums(input_error_fit > 0, na.rm = TRUE) > 0L),
        "trait(s) while also estimating shared measurement_error.\n"
    )
}
if (ncol(trait_matrix) >= (length(adj_data$tree$tip.label) - 1L)) {
    write_placeholder_outputs(
        tree = tree_full,
        trait_table = original_trait_table,
        reason = sprintf(
            "Too many traits for stable kfl1ou optimization (num_traits=%d, num_tips=%d). Falling back to placeholder outputs.",
            ncol(trait_matrix), length(adj_data$tree$tip.label)
        ),
        quit_after = TRUE
    )
}

alpha_bounds = resolve_alpha_upper(args[["alpha_upper"]], adj_data$tree)
cat("default alpha upper bound in kfl1ou =", alpha_bounds$default, "\n")
cat("the alpha upper bound actually used =", alpha_bounds$value, "\n")

max_shift_info = resolve_max_nshift(args[["max_nshift"]], adj_data$tree)
cat("maximum number of shifts =", max_shift_info$reported, "\n")

nslots = parse_integer_arg(args[["nslots"]], default = 1L)
if (is.na(nslots) || (nslots < 1L)) {
    nslots = 1L
}
nbootstrap = parse_integer_arg(args[["nbootstrap"]], default = 0L)
if (is.na(nbootstrap) || (nbootstrap < 0L)) {
    nbootstrap = 0L
}
detect_convergence = (parse_bool_flag(args[["detect_convergence"]], default = 1L) == 1L)
criterion = parse_string_arg(args[["criterion"]], default = "pBIC")

fit_ind = NULL
fit_conv = NULL
my_time = proc.time()

fit_ind_file = parse_string_arg(args[["fit_ind_file"]], default = NULL)
fit_conv_file = parse_string_arg(args[["fit_conv_file"]], default = NULL)

loaded_fit_ind = load_cached_fit(fit_ind_file, "fit_ind")
ind_flag = FALSE
if (!is.null(loaded_fit_ind)) {
    if (inherits(loaded_fit_ind, "try-error")) {
        cat(fit_ind_file, "is malformed. Printing...\n")
        print(loaded_fit_ind)
    } else if (inherits(loaded_fit_ind, "l1ou")) {
        if (fit_matches_current_data(loaded_fit_ind, adj_data$tree, trait_matrix)) {
            fit_ind = loaded_fit_ind
            cat(fit_ind_file, "was found. Loading the shift configurations.\n")
            ind_flag = TRUE
        } else {
            cat(fit_ind_file, "was found but does not match the current tree/trait matrix. Recomputing.\n")
        }
    }
}

if (!ind_flag) {
    cat("fit_ind.RData was not found. Starting shift detection.\n")
    fit_ind = try(
        kfl1ou::estimate_shift_configuration(
            tree = adj_data$tree,
            Y = trait_matrix,
            max.nShifts = max_shift_info$value,
            criterion = criterion,
            root.model = "OUfixedRoot",
            nCores = nslots,
            rescale = FALSE,
            alpha.upper = alpha_bounds$value,
            alpha.lower = NA,
            measurement_error = TRUE,
            input_error = input_error_fit,
            quietly = FALSE
        ),
        silent = TRUE
    )
    save_named_object("fit_ind", fit_ind, "fit_ind.RData")
    cat("Time elapsed for shift detection:\n")
}
print(proc.time() - my_time)

loaded_fit_conv = load_cached_fit(fit_conv_file, "fit_conv")
conv_flag = FALSE
if (!is.null(loaded_fit_conv)) {
    if (inherits(loaded_fit_conv, "try-error")) {
        cat(fit_conv_file, "is malformed. Printing...\n")
        print(loaded_fit_conv)
    } else if (inherits(loaded_fit_conv, "l1ou")) {
        if (fit_matches_current_data(loaded_fit_conv, adj_data$tree, trait_matrix)) {
            fit_conv = loaded_fit_conv
            cat(fit_conv_file, "was found. Loading the shift configurations.\n")
            conv_flag = TRUE
        } else {
            cat(fit_conv_file, "was found but does not match the current tree/trait matrix. Recomputing.\n")
        }
    }
}

if ((!conv_flag) && detect_convergence) {
    cat("Started convergence analysis.\n")
    if (inherits(fit_ind, "l1ou")) {
        if (fit_ind[["nShifts"]] < 2L) {
            fit_conv = fit_ind
        } else {
            fit_conv = try(
                kfl1ou::estimate_convergent_regimes(
                    fit_ind,
                    criterion = criterion,
                    method = "backward",
                    fixed.alpha = TRUE,
                    nCores = nslots
                ),
                silent = TRUE
            )
        }
    } else {
        warning("estimate_shift_configuration() failed.")
        if (inherits(fit_ind, "try-error")) {
            warning(attr(fit_ind, "condition"))
        }
    }
    save_named_object("fit_conv", fit_conv, "fit_conv.RData")
    cat("Time elapsed for shift detection + convergence detection:\n")
    print(fit_conv)
    print(proc.time() - my_time)
} else if (!conv_flag) {
    cat("Skipped convergence analysis.\n")
    fit_conv = fit_ind
} else {
    cat("Loaded fit_conv from cache. Skipping convergence analysis.\n")
}

if (inherits(fit_conv, "l1ou")) {
    fit_output = restore_fit_for_output(fit_conv, tree_full)

    if (nbootstrap > 0L) {
        cat("Start bootstrapping,", nbootstrap, "replicates,", nslots, "cores.\n")
        bootstrap_result = kfl1ou::l1ou_bootstrap_support(
            fit_conv,
            nItrs = nbootstrap,
            multicore = (nslots > 1L),
            nCores = nslots
        )
        bootstrap_output = restore_bootstrap_result_to_tree(fit_conv, bootstrap_result, tree_full)
        bootstrap_table = get_bootstrap_table(fit_output, bootstrap_output, mode = "l1ou")
        write.table(bootstrap_table, file = "l1ou_bootstrap.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
        cat("Time elapsed for shift detection + convergence detection + bootstrapping:\n")
        print(proc.time() - my_time)
    }

    cat("Start preparing output tables.\n")
    tab_out = safely_get_l1ou_tables(fit_output, fit_output$tree, original_trait_table)
    tree_table = tab_out[["tree_table"]]
    regime_table = tab_out[["regime_table"]]
    leaf_table = tab_out[["leaf_table"]]

    for (mt in removed_traits) {
        regime_table[[mt]] = NA
        leaf_table[[mt]] = 0
    }
    for (mt in original_traits[!(original_traits %in% colnames(regime_table))]) {
        regime_table[[mt]] = NA
    }
    for (mt in original_traits[!(original_traits %in% colnames(leaf_table))]) {
        leaf_table[[mt]] = 0
    }

    regime_table = regime_table[, c(colnames(regime_table)[!colnames(regime_table) %in% original_traits], original_traits), drop = FALSE]
    leaf_table = leaf_table[, c(colnames(leaf_table)[!colnames(leaf_table) %in% original_traits], original_traits), drop = FALSE]

    write.table(tree_table, file = "l1ou_tree.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(regime_table, file = "l1ou_regime.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(leaf_table, file = "l1ou_leaf.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

    grDevices::pdf("l1ou_plot.pdf",
                   height = length(fit_output$tree$tip.label) / 5 + 1,
                   width = max(8, ncol(adj_data$Y) / 2 + 5))
    par(mar = c(2, 2, 2, 2))
    try(plot(fit_output, edge.shift.ann = FALSE), silent = TRUE)
    grDevices::dev.off()
} else {
    warning("estimate_convergent_regimes() failed. Writing placeholder outputs.")
    if (inherits(fit_conv, "try-error")) {
        warning(attr(fit_conv, "condition"))
    }
    write_placeholder_outputs(
        tree = tree_full,
        trait_table = original_trait_table,
        fit_ind_obj = fit_ind,
        fit_conv_obj = fit_conv,
        reason = "kfl1ou fit did not converge to a valid object. Placeholder outputs were written.",
        quit_after = FALSE
    )
}

cat("kfl1ou completed!\n")
