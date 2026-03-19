enhance_branch_table = function(b, args, event_method=c('auto','generax','species_overlap')) {
    if ('mapdnds_omega' %in% colnames(b)) {
        b[is.na(b[['mapdnds_omega']]),'mapdnds_omega'] = 0
    }
    b[['branch_category']] = 'S'
    if ((event_method=='generax')|((event_method=='auto')&('generax_event' %in% colnames(b)))) {
        cat('generax_event_parent will be used to define branch events.\n')
        b[['is_parent_dup']] = (b[['generax_event_parent']]=='D')
        b[['is_parent_transfer']] = (b[['generax_event_parent']]=='H')
    } else if (event_method=='species_overlap') {
        cat('so_event_parent will be used to define branch events.\n')
        b[['is_parent_dup']] = (b[['so_event_parent']]=='D')
        b[['is_parent_transfer']] = FALSE
    } else {
        cat('Non-supported --event_method specified. Exiting.\n')
        quit()
    }
    if ('delta_intron_present' %in% colnames(b)) {
        cat('delta_intron_present will be used to define R branches.\n')
        cat('Nodes will be categorized into: S, D, R\n')
        b = attach_neighbor_stats(df=b, neighbor='sister', columns='delta_intron_present')
        b[['is_retrotransposition']] = (b[['delta_intron_present']] <= args[['max_delta_intron_present']])
        b[['is_retrotransposition']][is.na(b[['is_retrotransposition']])] = FALSE
        b[['is_sister_retrotransposition']] = (b[['sister_delta_intron_present']] <= args[['max_delta_intron_present']])
        b[['is_sister_retrotransposition']][is.na(b[['is_sister_retrotransposition']])] = FALSE
        b[['is_lower_delta_intron_present']] = (b[['delta_intron_present']] <= b[['sister_delta_intron_present']])    
        b[(b[['is_parent_dup']]),'branch_category'] = 'D'
        b[(b[['is_parent_dup']])&(b[['is_retrotransposition']])&(!b[['is_sister_retrotransposition']]),'branch_category'] = 'R'
    } else {
        cat('Nodes will be categorized into: S, D', '\n')
        b[(b[['is_parent_dup']]),'branch_category'] = 'D'
    }
    if (sum(b[['is_parent_transfer']])>0) {
        b[(b[['is_parent_transfer']]),'branch_category'] = 'H'
    }
    
    b = attach_neighbor_stats(df=b, neighbor='child1', columns='branch_category')
    b = attach_neighbor_stats(df=b, neighbor='child2', columns='branch_category')
    b[is.na(b[['child1_branch_category']]),'child1_branch_category'] = ''
    b[is.na(b[['child2_branch_category']]),'child2_branch_category'] = ''
    b[['node_category']] = ''
    child_cat = b[,c('child1_branch_category','child2_branch_category'), drop = FALSE]
    num_children_S = rowSums(child_cat == 'S')
    num_children_D = rowSums(child_cat == 'D')
    num_children_R = rowSums(child_cat == 'R')
    num_children_H = rowSums(child_cat == 'H')
    b[(num_children_S>=1),'node_category'] = 'S'
    b[(num_children_D>=1),'node_category'] = 'D'
    b[((num_children_D==1)&(num_children_R==1)),'node_category'] = 'R'
    b[(num_children_H>=1),'node_category'] = 'H'
    cat('detected node categories:', unique(b[['node_category']][b[['node_category']]!='NaN']), '\n')

    # Support values
    support_cols = colnames(b)[grepl('^support_', colnames(b))]
    root_nlabel = max(b[['branch_id']])
    subroot_nlabels = b[(b[['parent']]==root_nlabel),'branch_id']
    for (col in support_cols) {
        root_support = b[(b[['parent']]==-999),col]
        if (!is.na(root_support)) {
            cat('Warning:', col, ': support value at root node is not NA (', root_support, '). Masking.\n')
            cat('Non-NA root support is no problem if this is a pruned tree. If not, check the data carefully.\n')
            b[(b[['parent']]==-999),col] = NA
        }
        is_reconciled_branch = ((is.na(b[[col]]))&(b[['parent']]!=-999))
        b[is_reconciled_branch,col] = '-'
        # Suppress the support value of one of the subroot branches
        subroot_idx = which(b[['branch_id']] %in% subroot_nlabels)
        if (length(subroot_idx) > 0) {
            b[subroot_idx[length(subroot_idx)], col] = NA
        }
    }

    # dNdS
    omega_method = 'mapdnds'
    omega_column = paste0(omega_method, '_omega')
    if (omega_column %in% colnames(b)) {
        b = attach_neighbor_stats(df=b, neighbor='sister', columns=omega_column)
        b[['is_higher_dnds']] = (b[[omega_column]]>=b[[paste0('sister_', omega_column)]])
    }

    # csubst
    if ('csubst_is_foreground' %in% colnames(b)) {
        b[,'csubst_is_foreground'] = ifelse(b[['csubst_is_foreground']] %in% c('Y','yes'), TRUE, FALSE)
        b[,'csubst_is_marginal'] = ifelse(b[['csubst_is_marginal']] %in% c('Y','yes'), TRUE, FALSE)
    }

    # Branch width
    b[['branch_thickness']] = 0.25
    #b[(b[['is_higher_dnds']]&b[['is_parent_dup']]),'branch_thickness'] = 1

    #img_info = get_img_info(b, args[['dir_phylopic']])

    # sorting
    b = b[order(b[['branch_id']]),]
    rownames(b) = NULL
    return(b)
}

add_gene_cluster_membership = function(df_tip, max_bp_membership) {
    df_tip[,'cluster_membership'] = ''
    df_tip[,'species'] = ''
    is_annotated_leaf = (df_tip[['so_event']]=='L')
    is_annotated_leaf = is_annotated_leaf & (!is.na(df_tip[['start']]))
    is_annotated_leaf = is_annotated_leaf & (!is.na(df_tip[['end']]))
    spp = unique(df_tip[is_annotated_leaf, 'taxon', drop = FALSE][['taxon']])
    for (sp in spp) {
        sp_ub = sub(' ', '_', sp)
        gene_membership_counter = 0
        is_sp = ((df_tip[['so_event']]=='L')&(df_tip[['taxon']]==sp))
        df_tip[is_sp,'species'] = sp_ub
        chromosome_names = unique(df_tip[is_sp, 'chromosome', drop = FALSE][['chromosome']])
        chromosome_names = chromosome_names[chromosome_names!='']
        for (chromosome_name in chromosome_names) {
            is_chromosome = ((is_sp)&(df_tip[['chromosome']]==chromosome_name))
            is_chromosome[is.na(is_chromosome)] = FALSE
            n_chr = sum(is_chromosome)
            if (n_chr == 0) {
                next
            }
            b_cr = data.frame(df_tip[is_chromosome,], stringsAsFactors=FALSE)
            b_cr_start_new = apply(b_cr[,c('start','end')], 1, min)
            b_cr_end_new = apply(b_cr[,c('start','end')], 1, max)
            b_cr[,'start'] = b_cr_start_new
            b_cr[,'end'] = b_cr_end_new
            b_cr = b_cr[order(b_cr[['start']]), , drop=FALSE]
            rownames(b_cr) = NULL

            # Start a new cluster at each chromosome, then split by intergenic gaps.
            gene_membership_counter = gene_membership_counter + 1
            if (nrow(b_cr) == 1) {
                cluster_ids = gene_membership_counter
            } else {
                intergenic = b_cr[['start']][-1] - b_cr[['end']][-nrow(b_cr)]
                is_new_cluster = c(TRUE, intergenic > max_bp_membership)
                cluster_ids = (gene_membership_counter - 1L) + cumsum(is_new_cluster)
                gene_membership_counter = max(cluster_ids)
            }
            b_cr[['cluster_membership']] = paste0(sp_ub, '_', cluster_ids)
            df_tip[match(as.character(b_cr[['label']]), as.character(df_tip[['label']])), 'cluster_membership'] = b_cr[['cluster_membership']]
        }
    }
    return(df_tip)
}

add_cluster_membership_column = function(g, args, gname, max_bp_membership){
    cat(as.character(Sys.time()), 'Adding gene cluster membership column.\n')
    if (!'start'%in%colnames(g[['tree']][['data']])) {
        cat('Scaffold positions are not available. cluster_membership plot will not be displayed.\n')
        return(g)
    }   
    
    df_tip = get_df_tip(g[['tree']])
    df_tip = add_gene_cluster_membership(df_tip=df_tip, max_bp_membership=max_bp_membership)
    font_size = args[['font_size']]

    duplicated_membership = duplicated(df_tip[['cluster_membership']])
    duplicated_membership = duplicated_membership & (!is.na(df_tip[['cluster_membership']]))
    duplicated_membership = duplicated_membership & (df_tip[['cluster_membership']] != '')
    multimember_clusters = unique(df_tip[['cluster_membership']][duplicated_membership])
    is_multimember = (df_tip[['cluster_membership']]%in%multimember_clusters)
    is_na = (is.na(df_tip[['cluster_membership']]))|(df_tip[['cluster_membership']]=='')
    df_single = df_tip[(!is_multimember)&(!is_na),]
    df_multi = df_tip[(is_multimember)&(!is_na),]

    g[[gname]] = ggplot(mapping=aes(x=species, y=y)) +
        geom_blank(data=df_tip, aes(y=label)) + 
        geom_point(data=df_single, color='gray90', alpha=1) +
        geom_line(data=df_multi, mapping=aes(group=cluster_membership, color=cluster_membership), alpha=0.5) +
        geom_point(data=df_multi, mapping=aes(color=cluster_membership), alpha=1) +
        #aplot::ylim2(gg=g[['tree']]) +
        #ylim(c(1, max(df_tip[['y']]))) +
        xlab(paste0('Gene cluster\nmembership\n(max dist =\n', max_bp_membership, ' bp)')) +
        theme_void() +
        guides(
            fill=guide_legend(title=NULL, nrow=6, byrow=FALSE)
        ) +
        coord_cartesian(clip = "off") +
        theme(
            axis.text=element_blank(),
            axis.title.x=element_text(size=args[['font_size']]),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.x=element_blank(),
            legend.title=element_blank(),
            legend.text=element_blank(),
            legend.position="none",
            legend.key.size=unit(0.4, 'lines'),
            legend.box.just='center',
            plot.margin=unit(args[['margins']], "cm")
        )
  return(g)
}

add_synteny_column = function(g, args, gname, path_synteny, synteny_window = 5) {
    cat(as.character(Sys.time()), 'Adding synteny column.\n')
    if (is.null(path_synteny) || (!nzchar(path_synteny)) || (!file.exists(path_synteny))) {
        cat('Synteny table not found. Synteny column will not be added.\n')
        return(g)
    }
    synteny_window = suppressWarnings(as.integer(synteny_window))
    if (is.na(synteny_window) || (!is.finite(synteny_window)) || (synteny_window < 1)) {
        synteny_window = 5
    }
    df_tip = get_df_tip(g[['tree']])
    if (nrow(df_tip) < 2) {
        cat('The tree has fewer than two tips. Synteny column will not be added.\n')
        return(g)
    }
    df_syn = tryCatch(
        read.table(path_synteny, sep='\t', header=TRUE, stringsAsFactors=FALSE, quote='', comment.char='', check.names=FALSE),
        error=function(e) NULL
    )
    if (is.null(df_syn) || (nrow(df_syn) == 0)) {
        cat('Synteny table is empty. Synteny column will not be added.\n')
        return(g)
    }
    required_cols = c('node_name', 'offset', 'group_id')
    if (!all(required_cols %in% colnames(df_syn))) {
        cat('Synteny table lacks required columns (node_name, offset, group_id). Synteny column will not be added.\n')
        return(g)
    }
    if (!('direction' %in% colnames(df_syn))) {
        df_syn[['direction']] = ifelse(suppressWarnings(as.integer(df_syn[['offset']])) < 0, 'upstream', 'downstream')
    }
    df_syn[['offset']] = suppressWarnings(as.integer(df_syn[['offset']]))
    df_syn[['group_id']] = as.character(df_syn[['group_id']])
    df_syn = df_syn[
        is.finite(df_syn[['offset']]) &
        (!is.na(df_syn[['group_id']])) &
        (df_syn[['group_id']] != '') &
        (df_syn[['offset']] != 0),
        ,
        drop=FALSE
    ]
    if (nrow(df_syn) == 0) {
        cat('No valid synteny entries were found. Synteny column will not be added.\n')
        return(g)
    }
    df_syn = df_syn[abs(df_syn[['offset']]) <= synteny_window, , drop=FALSE]
    if (nrow(df_syn) == 0) {
        cat('No synteny entries are within the requested window. Synteny column will not be added.\n')
        return(g)
    }
    df_syn = merge(
        df_syn,
        df_tip[, c('label', 'y')],
        by.x='node_name',
        by.y='label',
        all.x=TRUE,
        sort=FALSE
    )
    df_syn = df_syn[is.finite(df_syn[['y']]), , drop=FALSE]
    if (nrow(df_syn) == 0) {
        cat('No synteny entries matched tree tip labels. Synteny column will not be added.\n')
        return(g)
    }
    df_syn[['direction']] = ifelse(df_syn[['offset']] < 0, 'upstream', 'downstream')
    df_syn[['x']] = synteny_window + 1 + df_syn[['offset']]

    # Keep neighbor occupancy before group filtering so blank edge positions stay empty.
    df_back = unique(df_syn[, c('x', 'y'), drop=FALSE])
    df_back = df_back[order(df_back[['y']], df_back[['x']], method='radix'), , drop=FALSE]

    df_group_tip = unique(df_syn[, c('group_id', 'node_name'), drop=FALSE])
    df_group_n = aggregate(node_name ~ group_id, data=df_group_tip, FUN=function(v) length(unique(v)))
    keep_groups = df_group_n[['group_id']][df_group_n[['node_name']] >= 2]
    if (length(keep_groups) == 0) {
        cat('No multi-tip synteny groups were detected. Synteny column will not be added.\n')
        return(g)
    }
    df_syn = df_syn[df_syn[['group_id']] %in% keep_groups, , drop=FALSE]
    if (nrow(df_syn) == 0) {
        cat('Synteny groups were filtered out. Synteny column will not be added.\n')
        return(g)
    }

    pick_nearest_offset = function(df) {
        if (nrow(df) <= 1) {
            return(df[1, , drop=FALSE])
        }
        if ('neighbor_gene' %in% colnames(df)) {
            df = df[order(abs(df[['offset']]), df[['neighbor_gene']], method='radix'), , drop=FALSE]
        } else {
            df = df[order(abs(df[['offset']]), method='radix'), , drop=FALSE]
        }
        return(df[1, , drop=FALSE])
    }

    # Keep one marker per (tip, direction, group), nearest to the focal gene.
    key = paste(df_syn[['node_name']], df_syn[['direction']], df_syn[['group_id']], sep='__')
    df_syn_split = split(df_syn, key, drop=TRUE)
    df_syn = do.call(
        rbind,
        lapply(df_syn_split, pick_nearest_offset)
    )
    rownames(df_syn) = NULL

    link_parts = list()
    link_id = 1L
    for (direction_name in c('upstream', 'downstream')) {
        df_dir = df_syn[df_syn[['direction']] == direction_name, , drop=FALSE]
        if (nrow(df_dir) == 0) {
            next
        }
        gids = unique(df_dir[['group_id']])
        for (gid in gids) {
            df_grp = df_dir[df_dir[['group_id']] == gid, , drop=FALSE]
            if (nrow(df_grp) < 2) {
                next
            }
            df_grp = df_grp[order(df_grp[['y']], method='radix'), , drop=FALSE]
            for (i in seq_len(nrow(df_grp) - 1)) {
                x_a = suppressWarnings(as.numeric(df_grp[i, 'x']))
                x_b = suppressWarnings(as.numeric(df_grp[i + 1, 'x']))
                y_a = suppressWarnings(as.numeric(df_grp[i, 'y']))
                y_b = suppressWarnings(as.numeric(df_grp[i + 1, 'y']))
                if (!(is.finite(x_a) && is.finite(x_b) && is.finite(y_a) && is.finite(y_b))) {
                    next
                }
                link_parts[[link_id]] = data.frame(
                    x = x_a,
                    y = y_a,
                    xend = x_b,
                    yend = y_b,
                    direction = direction_name,
                    link_id = paste0('synteny_link_', link_id),
                    group_id = gid,
                    stringsAsFactors=FALSE
                )
                link_id = link_id + 1L
            }
        }
    }
    if (length(link_parts) > 0) {
        df_link = do.call(rbind, link_parts)
    } else {
        df_link = data.frame(
            x = numeric(0),
            y = numeric(0),
            xend = numeric(0),
            yend = numeric(0),
            link_id = character(0),
            group_id = character(0),
            stringsAsFactors=FALSE
        )
    }

    df_point = df_syn[, c('x', 'y', 'group_id', 'node_name', 'offset', 'direction'), drop=FALSE]
    df_point = df_point[order(df_point[['y']], df_point[['x']], method='radix'), , drop=FALSE]
    group_levels = sort(unique(c(as.character(df_point[['group_id']]), as.character(df_link[['group_id']]))))
    point_palette = gg_color_hue(length(group_levels))
    names(point_palette) = group_levels
    df_point[['group_id']] = factor(df_point[['group_id']], levels=group_levels)
    if (nrow(df_link) > 0) {
        df_link[['group_id']] = factor(df_link[['group_id']], levels=group_levels)
    }

    x_center = synteny_window + 1
    x_max = synteny_window * 2 + 1
    x_tick_df = data.frame(
        x = c(1, synteny_window, synteny_window + 2, x_max),
        label = c(paste0('-', synteny_window), '-1', '+1', paste0('+', synteny_window)),
        stringsAsFactors=FALSE
    )
    x_tick_df = x_tick_df[!duplicated(x_tick_df[['x']]), , drop=FALSE]

    df_center = data.frame(
        x = rep(x_center, nrow(df_tip)),
        y = df_tip[['y']],
        stringsAsFactors=FALSE
    )

    link_layer = NULL
    if (nrow(df_link) > 0) {
        link_layer = geom_segment(
            data=df_link,
            mapping=aes(x=x, y=y, xend=xend, yend=yend, color=group_id, group=link_id),
            linewidth=1.8,
            lineend='round',
            alpha=0.5
        )
    }
    node_size = 1.5

    g[[gname]] = ggplot(data=df_point) +
        geom_blank(data=df_tip, aes(y=label)) +
        geom_segment(data=df_tip, mapping=aes(y=y, yend=y), x=1, xend=x_max, linewidth=0.25, color='gray90') +
        link_layer +
        geom_point(data=df_back, mapping=aes(x=x, y=y), color='gray90', size=node_size, alpha=1) +
        geom_point(data=df_center, mapping=aes(x=x, y=y), color='black', size=node_size, alpha=1) +
        geom_point(mapping=aes(x=x, y=y, color=group_id), size=node_size, alpha=1) +
        scale_x_continuous(
            breaks=x_tick_df[['x']],
            labels=x_tick_df[['label']],
            limits=c(0.5, x_max + 0.5)
        ) +
        scale_color_manual(values=point_palette) +
        xlab('Neighboring genes') +
        theme_minimal(base_size=args[['font_size']]) +
        guides(
            color='none'
        ) +
        coord_cartesian(clip='off') +
        theme(
            axis.title.y=element_blank(),
            axis.title.x=element_text(size=args[['font_size']]),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, color='black', size=args[['font_size']]),
            panel.grid.major.y=element_blank(),
            panel.grid.minor.y=element_blank(),
            panel.grid.minor.x=element_blank(),
            panel.grid.major.x=element_blank(),
            plot.margin=unit(args[['margins']], 'cm')
        )
    return(g)
}

