gg_color_hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

add_motif_colors = function(df_fimo, min_count=2) {
    motif_keys = as.character(df_fimo[['motif_altid']])
    motif_tab = sort(table(motif_keys), decreasing=TRUE)
    df_count = data.frame(
        altid = names(motif_tab),
        count = as.integer(motif_tab),
        stringsAsFactors = FALSE
    )
    if (nrow(df_count) == 0) {
        df_fimo[,'count'] = rep(NA_integer_, nrow(df_fimo))
        df_fimo[,'color_id'] = factor(rep(NA_integer_, nrow(df_fimo)))
        return(df_fimo)
    }
    df_count[,'color_id'] = seq_len(nrow(df_count))
    df_count[(df_count[['count']]<min_count),'color_id'] = 0
    color_ids = unname(df_count[['color_id']][match(motif_keys, df_count[['altid']])])
    df_fimo[,'count'] = unname(df_count[['count']][match(motif_keys, df_count[['altid']])])
    df_fimo[,'color_id'] = factor(color_ids, levels=sort(unique(df_count[,'color_id'])))
    return(df_fimo)
}

add_motif_counts = function(df_fimo) {
    motif_keys = as.character(df_fimo[['motif_altid']])
    motif_tab = table(motif_keys)
    df_fimo[,'count'] = as.integer(unname(motif_tab[motif_keys]))
    return(df_fimo)
}

get_df_polygon = function(df_fimo, min_count, max_count, multiple_connection='align_from_tss', alpha_max=0.5, ncpu=1, verbose=TRUE) {
    empty_polygon_df = function() {
        data.frame(
            x = numeric(0),
            y = numeric(0),
            polygon_id = factor(),
            motif_altid = factor(),
            alpha_value = numeric(0)
        )
    }
    n = nrow(df_fimo)
    if (n == 0L) {
        if (verbose) {
            cat('Number of motifs with displayed connections: 0 \n')
        }
        return(empty_polygon_df())
    }
    motif_raw = df_fimo[['motif_altid']]
    motif_vals = as.character(motif_raw)
    valid_motif = !is.na(motif_raw) & nzchar(motif_vals)
    motif_tab = table(motif_vals[valid_motif])
    keep_altids = names(motif_tab)[(motif_tab >= min_count) & (motif_tab <= max_count)]
    if (verbose) {
        cat('Number of motifs with displayed connections:', length(keep_altids), '\n')
    }
    if (length(keep_altids) == 0) {
        return(empty_polygon_df())
    }

    y_num = suppressWarnings(as.numeric(df_fimo[['y']]))
    start_num = suppressWarnings(as.numeric(df_fimo[['start']]))
    end_num = suppressWarnings(as.numeric(df_fimo[['end']]))
    ymin_num = suppressWarnings(as.numeric(df_fimo[['ymin']]))
    ymax_num = suppressWarnings(as.numeric(df_fimo[['ymax']]))
    keep_map = match(motif_vals, keep_altids, nomatch = 0L)
    valid_row = valid_motif & is.finite(y_num) & is.finite(start_num) & is.finite(end_num) &
        is.finite(ymin_num) & is.finite(ymax_num) & (keep_map > 0L)
    if (!any(valid_row)) {
        return(empty_polygon_df())
    }
    num_max = length(unique(y_num[valid_row]))
    if (num_max == 0L) {
        return(empty_polygon_df())
    }
    ncpu = suppressWarnings(as.integer(ncpu))
    if (!is.finite(ncpu) || is.na(ncpu) || (ncpu < 1L)) {
        ncpu = 1L
    }

    idx_keep = which(valid_row)
    idx_by_altid = split(idx_keep, motif_vals[idx_keep], drop = TRUE)
    build_altid_polygon_parts = function(altid) {
        idx_altid = idx_by_altid[[altid]]
        if (is.null(idx_altid) || (length(idx_altid) < 2L)) {
            return(NULL)
        }
        ord = order(y_num[idx_altid], start_num[idx_altid], end_num[idx_altid], method = 'radix')
        idx_altid = idx_altid[ord]
        n_altid = length(idx_altid)
        if (n_altid < 2L) {
            return(NULL)
        }
        yvals = y_num[idx_altid]
        starts = start_num[idx_altid]
        ends = end_num[idx_altid]
        ymins = ymin_num[idx_altid]
        ymaxs = ymax_num[idx_altid]
        y_change = c(TRUE, yvals[-1L] != yvals[-n_altid])
        group_starts = which(y_change)
        n_group = length(group_starts)
        if (n_group < 2L) {
            return(NULL)
        }
        group_ends = c(group_starts[-1L] - 1L, n_altid)
        alpha_value = (n_group / num_max) * alpha_max
        x_parts_local = list()
        y_parts_local = list()
        id_parts_local = list()
        motif_parts_local = list()
        alpha_parts_local = list()
        part_idx_local = 1L
        counter_local = 0L
        append_local = function(x1, x2, x3, x4, y1, y2, y3, y4) {
            n_poly = length(x1)
            if (n_poly == 0L) {
                return()
            }
            total_len = n_poly * 4L
            pos = seq.int(1L, total_len, by = 4L)
            x_vec = numeric(total_len)
            y_vec = numeric(total_len)
            x_vec[pos] = x1
            x_vec[pos + 1L] = x2
            x_vec[pos + 2L] = x3
            x_vec[pos + 3L] = x4
            y_vec[pos] = y1
            y_vec[pos + 1L] = y2
            y_vec[pos + 2L] = y3
            y_vec[pos + 3L] = y4
            poly_ids = counter_local + seq_len(n_poly) - 1L
            x_parts_local[[part_idx_local]] <<- x_vec
            y_parts_local[[part_idx_local]] <<- y_vec
            id_parts_local[[part_idx_local]] <<- rep(poly_ids, each = 4L)
            motif_parts_local[[part_idx_local]] <<- rep(altid, total_len)
            alpha_parts_local[[part_idx_local]] <<- rep(alpha_value, total_len)
            part_idx_local <<- part_idx_local + 1L
            counter_local <<- counter_local + n_poly
        }
        for (k in seq_len(n_group - 1L)) {
            current_ind = seq.int(group_starts[k], group_ends[k])
            next_ind = seq.int(group_starts[k + 1L], group_ends[k + 1L])
            n_current = length(current_ind)
            n_next = length(next_ind)
            if ((n_current == 0L) || (n_next == 0L)) {
                next
            }
            if ((n_current == 1L) || (n_next == 1L)) {
                if (multiple_connection == 'align_from_tss') {
                    current_ind = tail(current_ind, 1L)
                    next_ind = tail(next_ind, 1L)
                    n_current = length(current_ind)
                    n_next = length(next_ind)
                }
                if (n_current == 1L) {
                    ci_vec = rep(current_ind, n_next)
                    ni_vec = next_ind
                } else if (n_next == 1L) {
                    ci_vec = current_ind
                    ni_vec = rep(next_ind, n_current)
                } else {
                    ci_vec = rep(current_ind, each = n_next)
                    ni_vec = rep(next_ind, times = n_current)
                }
                if (length(ci_vec) == 0L) {
                    next
                }
                append_local(
                    starts[ci_vec], ends[ci_vec], ends[ni_vec], starts[ni_vec],
                    ymaxs[ci_vec], ymaxs[ci_vec], ymins[ni_vec], ymins[ni_vec]
                )
            } else {
                if (multiple_connection == 'align_from_tss') {
                    min_num_ind = min(n_current, n_next)
                    if (min_num_ind == 0L) {
                        next
                    }
                    current_ind2 = tail(current_ind, min_num_ind)
                    next_ind2 = tail(next_ind, min_num_ind)
                    append_local(
                        starts[current_ind2], ends[current_ind2], ends[next_ind2], starts[next_ind2],
                        ymaxs[current_ind2], ymaxs[current_ind2], ymins[next_ind2], ymins[next_ind2]
                    )
                } else if (multiple_connection == 'all') {
                    current_mid = (starts[current_ind] + ends[current_ind]) / 2
                    next_mid = (starts[next_ind] + ends[next_ind]) / 2
                    node_x = mean(c(current_mid, next_mid))
                    current_size = ends[current_ind] - starts[current_ind]
                    next_size = ends[next_ind] - starts[next_ind]
                    node_width = mean(c(current_size, next_size))
                    node_y = mean(c(mean(yvals[current_ind]), mean(yvals[next_ind])))
                    if (is.finite(node_y) && (node_y %% 1 == 0)) {
                        node_y = node_y + 0.5
                    }
                    node_x_start = node_x - (node_width / 2)
                    node_x_end = node_x + (node_width / 2)
                    append_local(
                        starts[current_ind], ends[current_ind], rep(node_x_end, n_current), rep(node_x_start, n_current),
                        ymaxs[current_ind], ymaxs[current_ind], rep(node_y, n_current), rep(node_y, n_current)
                    )
                    append_local(
                        starts[next_ind], ends[next_ind], rep(node_x_end, n_next), rep(node_x_start, n_next),
                        ymins[next_ind], ymins[next_ind], rep(node_y, n_next), rep(node_y, n_next)
                    )
                }
            }
        }
        if (part_idx_local == 1L) {
            return(NULL)
        }
        return(list(
            x = unlist(x_parts_local, use.names = FALSE),
            y = unlist(y_parts_local, use.names = FALSE),
            id = unlist(id_parts_local, use.names = FALSE),
            motif = unlist(motif_parts_local, use.names = FALSE),
            alpha = unlist(alpha_parts_local, use.names = FALSE),
            n_poly = as.integer(counter_local)
        ))
    }
    use_parallel = (ncpu > 1L) && (length(keep_altids) > 1L) && (.Platform$OS.type != 'windows')
    if (use_parallel) {
        motif_results = parallel::mclapply(keep_altids, build_altid_polygon_parts, mc.cores = ncpu)
        keep_result = !vapply(motif_results, is.null, logical(1))
        motif_results = motif_results[keep_result]
        if (length(motif_results) == 0L) {
            return(empty_polygon_df())
        }
        seg_len = vapply(motif_results, function(x) length(x[['x']]), integer(1))
        total_len = sum(seg_len)
        x_all = numeric(total_len)
        y_all = numeric(total_len)
        id_all = integer(total_len)
        motif_all = character(total_len)
        alpha_all = numeric(total_len)
        write_pos = 1L
        poly_offset = 0L
        for (ri in seq_along(motif_results)) {
            seg_n = seg_len[ri]
            idx = write_pos:(write_pos + seg_n - 1L)
            res = motif_results[[ri]]
            x_all[idx] = res[['x']]
            y_all[idx] = res[['y']]
            id_all[idx] = as.integer(res[['id']]) + poly_offset
            motif_all[idx] = res[['motif']]
            alpha_all[idx] = res[['alpha']]
            poly_offset = poly_offset + as.integer(res[['n_poly']])
            write_pos = write_pos + seg_n
        }
        df_polygon = data.frame(
            x = x_all,
            y = y_all,
            polygon_id = factor(id_all),
            motif_altid = factor(motif_all),
            alpha_value = alpha_all,
            stringsAsFactors = FALSE
        )
        return(df_polygon)
    }

    x_parts = list()
    y_parts = list()
    id_parts = list()
    motif_parts = list()
    alpha_parts = list()
    part_idx = 1L
    counter = 0L

    append_polygon_batch = function(x1, x2, x3, x4, y1, y2, y3, y4, altid, alpha_value) {
        n_poly = length(x1)
        if (n_poly == 0L) {
            return()
        }
        total_len = n_poly * 4L
        pos = seq.int(1L, total_len, by = 4L)
        x_vec = numeric(total_len)
        y_vec = numeric(total_len)
        x_vec[pos] = x1
        x_vec[pos + 1L] = x2
        x_vec[pos + 2L] = x3
        x_vec[pos + 3L] = x4
        y_vec[pos] = y1
        y_vec[pos + 1L] = y2
        y_vec[pos + 2L] = y3
        y_vec[pos + 3L] = y4
        poly_ids = counter + seq_len(n_poly) - 1L
        x_parts[[part_idx]] <<- x_vec
        y_parts[[part_idx]] <<- y_vec
        id_parts[[part_idx]] <<- rep(poly_ids, each = 4L)
        motif_parts[[part_idx]] <<- rep(altid, total_len)
        alpha_parts[[part_idx]] <<- rep(alpha_value, total_len)
        part_idx <<- part_idx + 1L
        counter <<- counter + n_poly
    }

    for (altid in keep_altids) {
        idx_altid = idx_by_altid[[altid]]
        if (is.null(idx_altid) || (length(idx_altid) < 2L)) {
            next
        }
        ord = order(y_num[idx_altid], start_num[idx_altid], end_num[idx_altid], method = 'radix')
        idx_altid = idx_altid[ord]
        n_altid = length(idx_altid)
        if (n_altid < 2L) {
            next
        }
        yvals = y_num[idx_altid]
        starts = start_num[idx_altid]
        ends = end_num[idx_altid]
        ymins = ymin_num[idx_altid]
        ymaxs = ymax_num[idx_altid]
        y_change = c(TRUE, yvals[-1L] != yvals[-n_altid])
        group_starts = which(y_change)
        n_group = length(group_starts)
        if (n_group < 2L) {
            next
        }
        group_ends = c(group_starts[-1L] - 1L, n_altid)
        alpha_value = (n_group / num_max) * alpha_max

        for (k in seq_len(n_group - 1L)) {
            current_ind = seq.int(group_starts[k], group_ends[k])
            next_ind = seq.int(group_starts[k + 1L], group_ends[k + 1L])
            n_current = length(current_ind)
            n_next = length(next_ind)
            if ((n_current == 0L) || (n_next == 0L)) {
                next
            }
            if ((n_current == 1L) || (n_next == 1L)) {
                if (multiple_connection == 'align_from_tss') {
                    current_ind = tail(current_ind, 1L)
                    next_ind = tail(next_ind, 1L)
                    n_current = length(current_ind)
                    n_next = length(next_ind)
                }
                if (n_current == 1L) {
                    ci_vec = rep(current_ind, n_next)
                    ni_vec = next_ind
                } else if (n_next == 1L) {
                    ci_vec = current_ind
                    ni_vec = rep(next_ind, n_current)
                } else {
                    ci_vec = rep(current_ind, each = n_next)
                    ni_vec = rep(next_ind, times = n_current)
                }
                if (length(ci_vec) == 0L) {
                    next
                }
                append_polygon_batch(
                    starts[ci_vec], ends[ci_vec], ends[ni_vec], starts[ni_vec],
                    ymaxs[ci_vec], ymaxs[ci_vec], ymins[ni_vec], ymins[ni_vec],
                    altid, alpha_value
                )
            } else {
                if (multiple_connection == 'align_from_tss') {
                    min_num_ind = min(n_current, n_next)
                    if (min_num_ind == 0L) {
                        next
                    }
                    current_ind2 = tail(current_ind, min_num_ind)
                    next_ind2 = tail(next_ind, min_num_ind)
                    append_polygon_batch(
                        starts[current_ind2], ends[current_ind2], ends[next_ind2], starts[next_ind2],
                        ymaxs[current_ind2], ymaxs[current_ind2], ymins[next_ind2], ymins[next_ind2],
                        altid, alpha_value
                    )
                } else if (multiple_connection == 'all') {
                    current_mid = (starts[current_ind] + ends[current_ind]) / 2
                    next_mid = (starts[next_ind] + ends[next_ind]) / 2
                    node_x = mean(c(current_mid, next_mid))
                    current_size = ends[current_ind] - starts[current_ind]
                    next_size = ends[next_ind] - starts[next_ind]
                    node_width = mean(c(current_size, next_size))
                    node_y = mean(c(mean(yvals[current_ind]), mean(yvals[next_ind])))
                    if (is.finite(node_y) && (node_y %% 1 == 0)) {
                        node_y = node_y + 0.5
                    }
                    node_x_start = node_x - (node_width / 2)
                    node_x_end = node_x + (node_width / 2)
                    append_polygon_batch(
                        starts[current_ind], ends[current_ind], rep(node_x_end, n_current), rep(node_x_start, n_current),
                        ymaxs[current_ind], ymaxs[current_ind], rep(node_y, n_current), rep(node_y, n_current),
                        altid, alpha_value
                    )
                    append_polygon_batch(
                        starts[next_ind], ends[next_ind], rep(node_x_end, n_next), rep(node_x_start, n_next),
                        ymins[next_ind], ymins[next_ind], rep(node_y, n_next), rep(node_y, n_next),
                        altid, alpha_value
                    )
                }
            }
        }
    }

    if (part_idx == 1L) {
        return(empty_polygon_df())
    }
    df_polygon = data.frame(
        x = unlist(x_parts, use.names = FALSE),
        y = unlist(y_parts, use.names = FALSE),
        polygon_id = factor(unlist(id_parts, use.names = FALSE)),
        motif_altid = factor(unlist(motif_parts, use.names = FALSE)),
        alpha_value = unlist(alpha_parts, use.names = FALSE),
        stringsAsFactors = FALSE
    )
    return(df_polygon)
}

merge_overlap_group = function(df_fimo, merge_level='TF', ncpu=1, verbose=TRUE) {
    if (merge_level=='TF') {
        altid_col = 'motif_altid'
    } else if (merge_level=='TFfamily') {
        altid_col = 'merged_altid'
        motif_raw = df_fimo[['motif_altid']]
        motif_chr = as.character(motif_raw)
        uniq_motif = unique(motif_chr)
        uniq_family = sub('([a-zA-Z])[-_]{0,1}[0-9]{1,3}$', '\\1', uniq_motif)
        mapped_family = uniq_family[match(motif_chr, uniq_motif)]
        mapped_family[is.na(motif_raw)] = NA_character_
        df_fimo[,altid_col] = mapped_family
    }
    df_fimo = add_complete_overlap_groups(df_fimo, altid_col, ncpu = ncpu, verbose = verbose)
    n = nrow(df_fimo)
    if (n == 0) {
        return(df_fimo[0, , drop = FALSE])
    }
    labels_chr = as.character(df_fimo[['label']])
    altid_chr = as.character(df_fimo[[altid_col]])
    starts = suppressWarnings(as.numeric(df_fimo[['start']]))
    ends = suppressWarnings(as.numeric(df_fimo[['end']]))
    swap_idx = is.finite(starts) & is.finite(ends) & (starts > ends)
    if (any(swap_idx)) {
        tmp_swap = starts[swap_idx]
        starts[swap_idx] = ends[swap_idx]
        ends[swap_idx] = tmp_swap
    }
    sort_order = order(labels_chr, altid_chr, starts, ends, method='radix', na.last=TRUE)
    labels_chr = labels_chr[sort_order]
    altid_chr = altid_chr[sort_order]
    starts = starts[sort_order]
    ends = ends[sort_order]
    y_vals = suppressWarnings(as.numeric(df_fimo[['y']]))[sort_order]
    strand_vals = as.character(df_fimo[['strand']])
    strand_vals = strand_vals[sort_order]
    finite_start = is.finite(starts)
    finite_end = is.finite(ends)

    label_levels = unique(labels_chr)
    altid_levels = unique(altid_chr)
    label_id = match(labels_chr, label_levels)
    altid_id = match(altid_chr, altid_levels)
    label_key = label_id
    altid_key = altid_id
    label_key[is.na(label_key)] = 0L
    altid_key[is.na(altid_key)] = 0L
    same_key = c(FALSE, (label_key[-1] == label_key[-n]) & (altid_key[-1] == altid_key[-n]))
    key_new = c(TRUE, !same_key[-1])
    key_starts = which(key_new)
    key_ends = c(key_starts[-1] - 1L, n)
    end_scan = ends
    end_scan[!finite_end] = -Inf
    prev_max_end = rep(-Inf, n)
    for (k in seq_len(length(key_starts))) {
        s = key_starts[k]
        e = key_ends[k]
        if (s < e) {
            prev_max_end[(s + 1L):e] = cummax(end_scan[s:(e - 1L)])
        }
    }
    new_group = key_new | !finite_start | (starts > prev_max_end)
    num_groups = sum(new_group)

    y_clean = y_vals
    y_clean[!is.finite(y_clean)] = Inf
    start_clean = starts
    start_clean[!finite_start] = Inf
    end_clean = ends
    end_clean[!finite_end] = -Inf
    strand_mask = integer(n)
    strand_not_na = !is.na(strand_vals)
    strand_mask[strand_not_na & (strand_vals == '+')] = 1L
    strand_mask[strand_not_na & (strand_vals == '-')] = 2L
    strand_mask[strand_not_na & ((strand_vals == '+-') | (strand_vals == '-+'))] = 3L
    # Fast bitwise-union lookup for masks in {0,1,2,3}.
    mask_or = matrix(c(
        0L, 1L, 2L, 3L,
        1L, 1L, 3L, 3L,
        2L, 3L, 2L, 3L,
        3L, 3L, 3L, 3L
    ), nrow = 4L, byrow = TRUE)
    out_label_id = rep(NA_integer_, num_groups)
    out_altid_id = rep(NA_integer_, num_groups)
    out_group = seq_len(num_groups)
    out_y = rep(NA_real_, num_groups)
    out_start = rep(NA_real_, num_groups)
    out_end = rep(NA_real_, num_groups)
    out_strand = rep(NA_character_, num_groups)
    current_group = 1L
    current_label_id = label_id[1L]
    current_altid_id = altid_id[1L]
    current_y_min = y_clean[1L]
    current_start_min = start_clean[1L]
    current_end_max = end_clean[1L]
    current_strand_mask = strand_mask[1L]
    current_first_raw_idx = if (strand_not_na[1L]) 1L else 0L

    if (n > 1L) {
        for (i in 2:n) {
            if (new_group[i]) {
                out_label_id[current_group] = current_label_id
                out_altid_id[current_group] = current_altid_id
                if (is.finite(current_y_min)) {
                    out_y[current_group] = current_y_min
                }
                if (is.finite(current_start_min)) {
                    out_start[current_group] = current_start_min
                }
                if (is.finite(current_end_max)) {
                    out_end[current_group] = current_end_max
                }
                if (current_strand_mask == 3L) {
                    out_strand[current_group] = '+-'
                } else if (current_strand_mask == 1L) {
                    out_strand[current_group] = '+'
                } else if (current_strand_mask == 2L) {
                    out_strand[current_group] = '-'
                } else if (current_first_raw_idx > 0L) {
                    out_strand[current_group] = strand_vals[current_first_raw_idx]
                }
                current_group = current_group + 1L
                current_label_id = label_id[i]
                current_altid_id = altid_id[i]
                current_y_min = y_clean[i]
                current_start_min = start_clean[i]
                current_end_max = end_clean[i]
                current_strand_mask = strand_mask[i]
                current_first_raw_idx = if (strand_not_na[i]) i else 0L
            } else {
                if (y_clean[i] < current_y_min) {
                    current_y_min = y_clean[i]
                }
                if (start_clean[i] < current_start_min) {
                    current_start_min = start_clean[i]
                }
                if (end_clean[i] > current_end_max) {
                    current_end_max = end_clean[i]
                }
                current_strand_mask = mask_or[current_strand_mask + 1L, strand_mask[i] + 1L]
                if ((current_first_raw_idx == 0L) && strand_not_na[i]) {
                    current_first_raw_idx = i
                }
            }
        }
    }
    out_label_id[current_group] = current_label_id
    out_altid_id[current_group] = current_altid_id
    if (is.finite(current_y_min)) {
        out_y[current_group] = current_y_min
    }
    if (is.finite(current_start_min)) {
        out_start[current_group] = current_start_min
    }
    if (is.finite(current_end_max)) {
        out_end[current_group] = current_end_max
    }
    if (current_strand_mask == 3L) {
        out_strand[current_group] = '+-'
    } else if (current_strand_mask == 1L) {
        out_strand[current_group] = '+'
    } else if (current_strand_mask == 2L) {
        out_strand[current_group] = '-'
    } else if (current_first_raw_idx > 0L) {
        out_strand[current_group] = strand_vals[current_first_raw_idx]
    }
    out_label = label_levels[out_label_id]
    out_altid = altid_levels[out_altid_id]
    df_out = data.frame(
        label = out_label,
        motif_altid = out_altid,
        overlap_group = out_group,
        y = out_y,
        start = out_start,
        end = out_end,
        strand = out_strand,
        stringsAsFactors = FALSE
    )
    rownames(df_out) = NULL
    return(df_out)
}

add_complete_overlap_groups = function(df_fimo, altid_col, ncpu=1, verbose=TRUE) {
    all_intervals_overlap = function(iv_a, iv_b) {
        n_a = length(iv_a[['start']])
        n_b = length(iv_b[['start']])
        if (n_a == 0L || n_b == 0L) {
            return(FALSE)
        }
        idx = findInterval(iv_a[['end']], iv_b[['start']])
        if (any(idx == 0L)) {
            return(FALSE)
        }
        return(all(iv_b[['prefix_end_max']][idx] >= iv_a[['start']]))
    }
    altid_raw = df_fimo[[altid_col]]
    altid_vals = as.character(altid_raw)
    valid_altid = !is.na(altid_raw) & nzchar(altid_vals)
    df_fimo[,'complete_overlap_group'] = rep(NA_integer_, nrow(df_fimo))
    if (sum(valid_altid) < 2L) {
        return(df_fimo)
    }
    by_altid_index = split(which(valid_altid), altid_vals[valid_altid], drop = TRUE)
    altids = names(by_altid_index)
    if (length(altids) < 2L) {
        return(df_fimo)
    }
    label_vals = as.character(df_fimo[['label']])
    label_codes = as.integer(factor(label_vals))
    starts_num_all = suppressWarnings(as.numeric(df_fimo[['start']]))
    ends_num_all = suppressWarnings(as.numeric(df_fimo[['end']]))
    genes_by_altid = lapply(by_altid_index, function(idx) {
        codes = unique(label_codes[idx])
        codes = codes[!is.na(codes)]
        if (length(codes) <= 1L) {
            return(codes)
        }
        sort.int(codes, method = 'quick')
    })
    gene_signature = vapply(genes_by_altid, function(g) paste(g, collapse = ","), character(1))
    candidate_groups = split(names(gene_signature), gene_signature, drop = TRUE)
    candidate_groups = candidate_groups[vapply(candidate_groups, length, integer(1)) > 1L]
    if (length(candidate_groups) == 0) {
        return(df_fimo)
    }
    candidate_altids = unique(unlist(candidate_groups, use.names=FALSE))
    intervals_by_altid_gene = setNames(vector('list', length(candidate_altids)), candidate_altids)
    for (altid in candidate_altids) {
        idx_altid = by_altid_index[[altid]]
        label_index = split(idx_altid, label_codes[idx_altid], drop = TRUE)
        genes = genes_by_altid[[altid]]
        n_genes = length(genes)
        if (n_genes == 0L) {
            intervals_by_altid_gene[[altid]] = list()
            next
        }
        gene_intervals = vector('list', n_genes)
        for (gi in seq_len(n_genes)) {
            idx = label_index[[as.character(genes[gi])]]
            if (is.null(idx) || (length(idx) == 0L)) {
                gene_intervals[[gi]] = NULL
                next
            }
            starts = starts_num_all[idx]
            ends = ends_num_all[idx]
            is_valid = is.finite(starts) & is.finite(ends)
            if (!any(is_valid)) {
                gene_intervals[[gi]] = NULL
                next
            }
            starts = starts[is_valid]
            ends = ends[is_valid]
            start_sorted = pmin(starts, ends)
            end_sorted = pmax(starts, ends)
            ord = order(start_sorted, end_sorted, method = 'radix')
            start_sorted = start_sorted[ord]
            end_sorted = end_sorted[ord]
            gene_intervals[[gi]] = list(
                start = start_sorted,
                end = end_sorted,
                prefix_end_max = cummax(end_sorted),
                min_start = start_sorted[1L],
                max_end = end_sorted[length(end_sorted)]
            )
        }
        intervals_by_altid_gene[[altid]] = gene_intervals
    }

    ncpu = suppressWarnings(as.integer(ncpu))
    if (!is.finite(ncpu) || is.na(ncpu) || (ncpu < 1L)) {
        ncpu = 1L
    }
    get_components_in_group = function(grp) {
        n_grp = length(grp)
        if (n_grp < 2L) {
            return(list())
        }
        parent = seq_len(n_grp)
        uf_find_local = function(x) {
            while (parent[x] != x) {
                parent[x] <<- parent[parent[x]]
                x <- parent[x]
            }
            x
        }
        uf_union_local = function(a, b) {
            ra = uf_find_local(a)
            rb = uf_find_local(b)
            if (ra != rb) {
                parent[rb] <<- ra
            }
        }
        for (ii in seq_len(n_grp - 1L)) {
            iv1 = intervals_by_altid_gene[[grp[ii]]]
            n_gene = length(iv1)
            for (jj in seq.int(ii + 1L, n_grp)) {
                if (uf_find_local(ii) == uf_find_local(jj)) {
                    next
                }
                iv2 = intervals_by_altid_gene[[grp[jj]]]
                if ((n_gene == 0L) || (length(iv2) != n_gene)) {
                    next
                }
                is_overlap = TRUE
                for (gi in seq_len(n_gene)) {
                    g1 = iv1[[gi]]
                    g2 = iv2[[gi]]
                    if (is.null(g1) || is.null(g2)) {
                        is_overlap = FALSE
                        break
                    }
                    if ((g1[['min_start']] > g2[['max_end']]) || (g2[['min_start']] > g1[['max_end']])) {
                        is_overlap = FALSE
                        break
                    }
                    if (!all_intervals_overlap(g1, g2) || !all_intervals_overlap(g2, g1)) {
                        is_overlap = FALSE
                        break
                    }
                }
                if (is_overlap) {
                    uf_union_local(ii, jj)
                }
            }
        }
        roots = vapply(seq_len(n_grp), uf_find_local, integer(1))
        comps = split(seq_len(n_grp), roots)
        comps = comps[vapply(comps, length, integer(1)) > 1L]
        if (length(comps) == 0L) {
            return(list())
        }
        lapply(comps, function(idx) grp[idx])
    }
    use_parallel = (ncpu > 1L) && (length(candidate_groups) > 1L) && (.Platform$OS.type != 'windows')
    if (use_parallel) {
        component_lists = parallel::mclapply(candidate_groups, get_components_in_group, mc.cores = ncpu)
    } else {
        component_lists = lapply(candidate_groups, get_components_in_group)
    }
    complete_overlap_groups = unlist(component_lists, recursive = FALSE, use.names = FALSE)

    n_altid = length(altids)
    altid_to_idx = setNames(seq_len(n_altid), altids)
    if (length(complete_overlap_groups)>0) {
        if (verbose) {
            for (i in seq_len(length(complete_overlap_groups))) {
                members = complete_overlap_groups[[i]]
                preview_n = min(6L, length(members))
                preview = paste(members[seq_len(preview_n)], collapse=' ')
                if (length(members) > preview_n) {
                    preview = paste0(preview, ' ... (', length(members), ' motifs)')
                }
                cat('Detected universal overlap (motifs will be merged):', preview, '\n')
            }
        }
        mapped_altid = rep(NA_character_, n_altid)
        mapped_group = rep(NA_integer_, n_altid)
        counter = 1L
        for (i in seq_len(length(complete_overlap_groups))) {
            member_altids = complete_overlap_groups[[i]]
            member_idx = as.integer(altid_to_idx[member_altids])
            my_order = order(nchar(member_altids))
            new_altid = paste(member_altids[my_order], collapse=',')
            if (nchar(new_altid)>20) {
                new_altid = paste0(substr(new_altid, 1, 20), '...')
            }
            mapped_altid[member_idx] = new_altid
            mapped_group[member_idx] = counter
            counter = counter + 1
        }
        idx = as.integer(altid_to_idx[altid_vals])
        is_mapped = valid_altid & !is.na(idx) & !is.na(mapped_group[idx])
        if (any(is_mapped)) {
            df_fimo[is_mapped, 'complete_overlap_group'] = mapped_group[idx[is_mapped]]
            altid_vals[is_mapped] = mapped_altid[idx[is_mapped]]
            df_fimo[[altid_col]] = altid_vals
        }
    }
    return(df_fimo)
}

get_df_N = function(df_tip) {
    df_tip = data.frame(df_tip, stringsAsFactors=FALSE)
    cols = c('x_start', 'x_end', 'y_start', 'y_end')
    x_start_parts = list()
    x_end_parts = list()
    y_start_parts = list()
    y_end_parts = list()
    part_idx = 1L
    for (i in seq_len(nrow(df_tip))) {
        N_slice = df_tip[i,'promoter_N']
        if (is.na(N_slice)) {
            next
        }
        if (N_slice=='') {
            next
        }
        y_val = suppressWarnings(as.numeric(df_tip[i,'y']))
        if (!is.finite(y_val)) {
            next
        }
        slice_vec = strsplit(as.character(N_slice), ',', fixed = TRUE)[[1]]
        if (length(slice_vec) == 0L) {
            next
        }
        for (j in seq_along(slice_vec)) {
            token = trimws(slice_vec[j])
            if (token == '') {
                next
            }
            if (!grepl(':', token, fixed = TRUE)) {
                x = suppressWarnings(as.numeric(token))
                if (!is.finite(x)) {
                    next
                }
                x_start_parts[[part_idx]] = x - 0.5
                x_end_parts[[part_idx]] = x + 0.5
                y_start_parts[[part_idx]] = y_val
                y_end_parts[[part_idx]] = y_val
                part_idx = part_idx + 1L
            } else {
                start_stop = strsplit(token, ':', fixed = TRUE)[[1]]
                if (length(start_stop) < 2L) {
                    next
                }
                x1 = suppressWarnings(as.numeric(start_stop[1]))
                x2 = suppressWarnings(as.numeric(start_stop[2]))
                if (!is.finite(x1) || !is.finite(x2)) {
                    next
                }
                lo = min(x1, x2)
                hi = max(x1, x2)
                x_start_parts[[part_idx]] = lo - 0.5
                x_end_parts[[part_idx]] = hi + 0.5
                y_start_parts[[part_idx]] = y_val
                y_end_parts[[part_idx]] = y_val
                part_idx = part_idx + 1L
            }
        }
    }
    if (part_idx == 1L) {
        df_N = data.frame(matrix(numeric(0), nrow = 0, ncol = length(cols)))
    } else {
        df_N = data.frame(
            x_start = as.numeric(unlist(x_start_parts, use.names = FALSE)),
            x_end = as.numeric(unlist(x_end_parts, use.names = FALSE)),
            y_start = as.numeric(unlist(y_start_parts, use.names = FALSE)),
            y_end = as.numeric(unlist(y_end_parts, use.names = FALSE)),
            stringsAsFactors = FALSE
        )
    }
    colnames(df_N) = cols
    return(df_N)
}

add_fimo_column = function(g, args, qname, xmax=2000, qvalue=0.01, multiple_connection='align_from_tss', ncpu=1) {
    cat(as.character(Sys.time()), 'Adding fimo column.\n')
    if (!'fimo_start'%in%colnames(g[['tree']][['data']])) {
        cat('Scaffold positions are not available. FIMO plot will not be displayed.\n')
        return(g)
    }
    ncpu = suppressWarnings(as.integer(ncpu))
    if (!is.finite(ncpu) || is.na(ncpu) || (ncpu < 1L)) {
        ncpu = 1L
    }
    df_tip = get_df_tip(g[['tree']])
    df_tip = data.frame(df_tip, stringsAsFactors=FALSE)
    to_char_no_na = function(x) {
        v = as.character(x)
        v[is.na(v)] = ''
        v
    }
    split_column_tokens = function(v) {
        out = strsplit(v, ';', fixed = TRUE)
        if (length(out) == 0L) {
            return(out)
        }
        empty_idx = (v == '')
        if (any(empty_idx)) {
            out[empty_idx] = replicate(sum(empty_idx), character(0), simplify = FALSE)
        }
        has_space = grepl(' ', v, fixed = TRUE)
        if (any(has_space)) {
            out[has_space] = lapply(out[has_space], trimws)
            empty_idx2 = has_space & (v == '')
            if (any(empty_idx2)) {
                out[empty_idx2] = replicate(sum(empty_idx2), character(0), simplify = FALSE)
            }
        }
        out
    }
    fimo_start_list = split_column_tokens(to_char_no_na(df_tip[['fimo_start']]))
    fimo_end_list = split_column_tokens(to_char_no_na(df_tip[['fimo_end']]))
    fimo_strand_list = split_column_tokens(to_char_no_na(df_tip[['fimo_strand']]))
    fimo_altid_list = split_column_tokens(to_char_no_na(df_tip[['fimo_alt_id']]))
    fimo_qvalue_list = split_column_tokens(to_char_no_na(df_tip[['fimo_qvalue']]))
    labels_vec = as.character(df_tip[['label']])
    y_vec = suppressWarnings(as.numeric(df_tip[['y']]))
    parse_tip_fimo = function(i) {
        raw_start = fimo_start_list[[i]]
        if (length(raw_start) == 0L) {
            return(list(core = 0L, pass = 0L))
        }
        fimo_starts = suppressWarnings(as.integer(raw_start))
        fimo_ends = suppressWarnings(as.integer(fimo_end_list[[i]]))
        fimo_strands = fimo_strand_list[[i]]
        fimo_altids = fimo_altid_list[[i]]
        fimo_qvalues = suppressWarnings(as.numeric(fimo_qvalue_list[[i]]))
        n = min(length(fimo_starts), length(fimo_ends), length(fimo_strands), length(fimo_altids), length(fimo_qvalues))
        if (n == 0L) {
            return(list(core = 0L, pass = 0L))
        }
        if (length(fimo_starts) != n) {
            fimo_starts = fimo_starts[seq_len(n)]
        }
        if (length(fimo_ends) != n) {
            fimo_ends = fimo_ends[seq_len(n)]
        }
        if (length(fimo_strands) != n) {
            fimo_strands = fimo_strands[seq_len(n)]
        }
        if (length(fimo_altids) != n) {
            fimo_altids = fimo_altids[seq_len(n)]
        }
        if (length(fimo_qvalues) != n) {
            fimo_qvalues = fimo_qvalues[seq_len(n)]
        }
        keep_core = is.finite(fimo_starts) & is.finite(fimo_ends) & is.finite(fimo_qvalues) &
            !is.na(fimo_altids) & (fimo_altids != '')
        n_core = sum(keep_core)
        if (n_core == 0L) {
            return(list(core = 0L, pass = 0L))
        }
        keep_q = keep_core & (fimo_qvalues <= qvalue)
        n_pass = sum(keep_q)
        if (n_pass == 0L) {
            return(list(core = n_core, pass = 0L))
        }
        return(list(
            core = n_core,
            pass = n_pass,
            row_idx = rep.int(i, n_pass),
            start = fimo_starts[keep_q],
            end = fimo_ends[keep_q],
            strand = fimo_strands[keep_q],
            qvalue = fimo_qvalues[keep_q],
            motif_altid = fimo_altids[keep_q]
        ))
    }
    index_vec = seq_len(nrow(df_tip))
    use_parallel = (ncpu > 1L) && (length(index_vec) > 1L) && (.Platform$OS.type != 'windows')
    if (use_parallel) {
        parsed_parts = parallel::mclapply(index_vec, parse_tip_fimo, mc.cores = ncpu)
    } else {
        if ((ncpu > 1L) && (.Platform$OS.type == 'windows')) {
            cat('ncpu > 1 was requested, but mclapply is unsupported on Windows. Falling back to ncpu=1.\n')
        }
        parsed_parts = vector('list', length(index_vec))
        for (i in index_vec) {
            parsed_parts[[i]] = parse_tip_fimo(i)
        }
    }
    core_counts = vapply(parsed_parts, function(x) x[['core']], integer(1))
    pass_counts = vapply(parsed_parts, function(x) x[['pass']], integer(1))
    num_core_records = sum(core_counts)
    num_pass_records = sum(pass_counts)
    if (num_core_records == 0L) {
        cat('No FIMO records are available. FIMO plot will not be displayed.\n')
        return(g)
    }
    if (num_pass_records == 0L) {
        cat('No FIMO records passed the qvalue cutoff. FIMO plot will not be displayed.\n')
        return(g)
    }
    row_idx_out = integer(num_pass_records)
    start_out = numeric(num_pass_records)
    end_out = numeric(num_pass_records)
    strand_out = character(num_pass_records)
    qvalue_out = numeric(num_pass_records)
    motif_out = character(num_pass_records)
    write_pos = 1L
    for (i in seq_along(parsed_parts)) {
        n_pass = pass_counts[i]
        if (n_pass == 0L) {
            next
        }
        idx = write_pos:(write_pos + n_pass - 1L)
        part = parsed_parts[[i]]
        row_idx_out[idx] = part[['row_idx']]
        start_out[idx] = part[['start']]
        end_out[idx] = part[['end']]
        strand_out[idx] = part[['strand']]
        qvalue_out[idx] = part[['qvalue']]
        motif_out[idx] = part[['motif_altid']]
        write_pos = write_pos + n_pass
    }
    df_fimo = data.frame(
        y = y_vec[row_idx_out],
        label = labels_vec[row_idx_out],
        start = start_out,
        end = end_out,
        strand = strand_out,
        qvalue = qvalue_out,
        motif_altid = motif_out,
        stringsAsFactors = FALSE
    )
    df_fimo = merge_overlap_group(df_fimo, merge_level='TFfamily', ncpu=ncpu, verbose=FALSE)
    strand_chr = as.character(df_fimo[['strand']])
    delta_up = rep(0, nrow(df_fimo))
    delta_down = rep(0, nrow(df_fimo))
    hjust_vals = rep(0.5, nrow(df_fimo))
    is_plus = (strand_chr == '+')
    is_minus = (strand_chr == '-')
    is_both = (strand_chr == '+-') | (strand_chr == '-+')
    delta_up[is_plus] = 0.3
    delta_down[is_minus] = 0.3
    delta_up[is_both] = 0.15
    delta_down[is_both] = 0.15
    hjust_vals[is_plus] = -0.05
    hjust_vals[is_minus] = 1.05
    df_fimo[,'ymin'] = df_fimo[['y']] - delta_down
    df_fimo[,'ymax'] = df_fimo[['y']] + delta_up
    df_fimo[,'hjust'] = hjust_vals
    df_fimo[,'start'] = xmax - df_fimo[['start']]
    df_fimo[,'end'] = xmax - df_fimo[['end']]
    df_fimo[,'xmid'] = (df_fimo[,'start']+df_fimo[,'end'])/2
    df_fimo = add_motif_counts(df_fimo)
    df_fimo = df_fimo[order(df_fimo[['y']], df_fimo[['start']]),]
    df_polygon = get_df_polygon(df_fimo, min_count=2, max_count=Inf, multiple_connection=multiple_connection, ncpu=ncpu, verbose=FALSE)
    legend_texts = unique(df_fimo[['motif_altid']])
    unique_altid = unique(df_fimo[(df_fimo[['count']]==1),'motif_altid'])
    my_palette = rep('gray50', length(legend_texts))
    num_colored = length(legend_texts) - length(unique_altid)
    if (num_colored > 0L) {
        my_palette[!legend_texts%in%unique_altid] = gg_color_hue(num_colored)
    }
    df_fimo[,'motif_altid'] = factor(df_fimo[,'motif_altid'], levels=legend_texts)
    df_N = get_df_N(df_tip)
    df_N[,'x_start'] = xmax - df_N[['x_start']]
    df_N[,'x_end'] = xmax - df_N[['x_end']]
    df_backbone = df_tip[(df_tip[['promoter_available']]=='Y'),c('label','y')]
    for (power in 0:10) {
        if ((10**(power+1))>xmax) {
            break
        }
    }
    xbreaks = (0:(xmax%/%(10**power))) * (10**power)
    xlabels = xbreaks - xmax
    polygon_layer = NULL
    alpha_scale = NULL
    if (nrow(df_polygon) > 0L) {
        alpha_vals = df_polygon[['alpha_value']]
        alpha_vals = alpha_vals[is.finite(alpha_vals)]
        alpha_range = c(0, 0.5)
        if (length(alpha_vals) > 0L) {
            alpha_range = range(alpha_vals)
            if (alpha_range[1] == alpha_range[2]) {
                alpha_range = c(max(0, alpha_range[1] - 1e-6), alpha_range[2] + 1e-6)
            }
        }
        polygon_layer = geom_polygon(data=df_polygon, mapping=aes(x=x, y=y, fill=motif_altid, group=polygon_id, alpha=alpha_value))
        alpha_scale = scale_alpha_continuous(range=alpha_range)
    }
    g[[qname]] = ggplot(data=df_fimo) +
        geom_blank(data=df_tip, aes(y=label)) + 
        geom_segment(data=df_backbone, mapping=aes(y=y, yend=y), x=1, xend=xmax, linewidth=0.1, color='black') +
        geom_segment(data=df_N, mapping=aes(x=x_start, xend=x_end, y=y_start, yend=y_end), linewidth=0.1, color='gray70') +
        geom_rect(mapping=aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=motif_altid), alpha=0.5) +
        geom_text(mapping=aes(x=xmid, y=y, label=motif_altid, color=motif_altid, hjust=hjust), size=0.2, angle=90, alpha=0.1) +
        #geom_text_repel(mapping=aes(x=xmid, y=y, label=motif_altid), size=0.15, angle=90, max.overlaps=1000, segment.size=0.1) +
        polygon_layer +
        scale_x_continuous(breaks=xbreaks, labels=xlabels, limits=c(-0.5, xmax+0.5)) +
        scale_fill_manual(values=my_palette) +
        scale_color_manual(values=my_palette) +
        alpha_scale +
        #aplot::ylim2(gg=g[['tree']]) +
        xlab('Motif position') +
        theme_minimal(base_size=args[['font_size']]) +
        guides(
            fill=guide_legend(title="", nrow=6, byrow=FALSE),
            color='none',
            alpha='none'
        ) +
        coord_cartesian(clip = "off") +
        theme(
            axis.title.y=element_blank(),
            axis.title.x=element_text(size=args[['font_size']]),
            axis.line=element_blank(),
            axis.ticks=element_blank(),
            axis.text.y=element_blank(),
            axis.text.x=element_text(size=args[['font_size']], color='black'),
            panel.grid.major.y=element_blank(),
            panel.grid.minor.y=element_blank(),
            panel.grid.minor.x=element_blank(),
            legend.position="bottom",
            legend.justification = "right",
            legend.title=element_text(size=args[['font_size']]),
            legend.text=element_text(size=args[['font_size']]),
            legend.box.just='center',
            legend.key.size=unit(0.4, 'lines'),
            rect=element_rect(fill="transparent"),
            plot.margin=unit(args[['margins']], "cm")
        )
    return(g)
}

extract_alignment_sites <- function(tidy_aln, selected_amino_acid_sites) {
  # Filter for only the selected positions
  tidy_aln <- tidy_aln[tidy_aln[["position"]] %in% selected_amino_acid_sites, ]
  
  # Keep the original position
  tidy_aln[["original_position"]] <- tidy_aln[["position"]]
  
  # Re-index positions using match() rather than a loop
  unique_positions <- unique(tidy_aln[["position"]])
  tidy_aln[["position"]] <- match(tidy_aln[["position"]], unique_positions)
  
  rownames(tidy_aln) <- NULL
  tidy_aln
}

add_amino_acid_site_column <- function(g, args, tidy_aln, selected_amino_acid_sites) {
  cat(as.character(Sys.time()), "Adding amino acid site column.\n")
  df_tip <- get_df_tip(g[["tree"]])[,c('label','branch_id')]
  tidy_site <- extract_alignment_sites(tidy_aln, selected_amino_acid_sites)
  tidy_site <- merge(tidy_site, df_tip, by.x = "name", by.y = "label", all.x = TRUE)
  tidy_site[["name"]] <- factor(tidy_site[["name"]], levels = df_tip[["label"]])
  qname <- "amino_acid_site"
    g[[qname]] <- ggplot() + 
    geom_blank(data = df_tip, aes(y = label)) + 
    ggmsa::geom_msa(data = tidy_site, color = "LETTER") +
    scale_x_continuous(
        breaks = 1:length(selected_amino_acid_sites),
        labels = selected_amino_acid_sites
    ) +
    xlab("Amino acid position (aa)") +
    theme(
        axis.title.y        = element_blank(),
        axis.title.x        = element_text(size = args[["font_size"]]),
        axis.text.x         = element_text(angle = 90, hjust = 1, vjust = 0.5, color = 'black', size = args[["font_size"]]),
        axis.ticks.x        = element_blank(),
        axis.line.x         = element_blank(),
        axis.text.y         = element_blank(),
        axis.ticks.y        = element_blank(),
        axis.line.y         = element_blank(), 
        panel.grid.major.y  = element_blank(),
        rect                = element_blank(),
        plot.margin         = unit(args[["margins"]] / 4, "cm")
    )
  
  g
}
