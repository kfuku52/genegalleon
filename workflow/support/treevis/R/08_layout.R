# Limit only graphical keys, preserving the complete fill scale and plotted hits.
# Occurrence counts are attached by the domain/FIMO producers before geometric
# splitting or polygon expansion, so repeated drawing does not inflate ranks.
treevis_fit_fill_legend = function(p, width_mm, name) {
    counts = attr(p, 'treevis_legend_counts')
    if (is.null(counts) || length(counts) == 0) return(p)
    keys = names(counts)[order(-as.numeric(counts), names(counts), method='radix')]
    candidate = function(n) {
        scale = p$scales$get_scales('fill')$clone()
        scale$breaks = keys[seq_len(n)]
        result = suppressMessages(p + scale) + ggplot2::guides(
            fill=if (n == 0) 'none' else ggplot2::guide_legend(
                title=NULL, nrow=min(6L, n), byrow=FALSE)) +
            ggplot2::theme(legend.margin=ggplot2::margin(0, 0, 0, 0),
                legend.text=ggplot2::element_text(margin=ggplot2::margin(l=1, r=0)))
        result
    }
    # With six rows, adding entries to this ranked prefix never reduces width.
    low = 0L
    high = length(keys)
    result = candidate(0L)
    while (low < high) {
        n = ceiling((low + high + 1L) / 2)
        trial = candidate(n)
        legend = cowplot::get_legend(trial)
        measured = if (is.null(legend)) 0 else grid::convertWidth(
            grid::grobWidth(legend), 'mm', valueOnly=TRUE)
        if (measured <= width_mm) {
            low = n
            result = trial
        } else {
            high = n - 1L
        }
    }
    attr(result, 'treevis_legend_visible') = keys[seq_len(low)]
    cat(name, 'legend:', low, 'of', length(keys), 'entries fit; omitted',
        length(keys)-low, 'lower-frequency keys.\n')
    result
}

treevis_keep_edge_ticks_inside = function(grob) {
    if (inherits(grob, 'text') && all(grob$rot == 0) && length(grob$x) > 0) {
        units = grid::unitType(grob$x)
        if (all(units %in% c('npc', 'native'))) {
            x = as.numeric(grob$x)
            grob$hjust = ifelse(x < 0.15, 0, ifelse(x > 0.85, 1, 0.5))
        }
    }
    for (i in seq_along(grob$children)) {
        grob$children[[i]] = treevis_keep_edge_ticks_inside(grob$children[[i]])
    }
    for (i in seq_along(grob$grobs)) {
        grob$grobs[[i]] = treevis_keep_edge_ticks_inside(grob$grobs[[i]])
    }
    grob
}

# Physical widths are resolved on a PDF device, before combining the columns.
# Overrides are lower bounds for the data panel, never fractions of the page.
treevis_layout_mm = function(g, panel_widths_mm = NULL, height_mm = NULL) {
    overrides = numeric()
    if (!is.null(panel_widths_mm) && !is.na(panel_widths_mm) && nzchar(panel_widths_mm)) {
        entries = strsplit(panel_widths_mm, ",", fixed=TRUE)[[1]]
        for (entry in entries) {
            pair = strsplit(entry, ":", fixed=TRUE)[[1]]
            value = if (length(pair) == 2) suppressWarnings(as.numeric(pair[2])) else NA_real_
            if (!is.finite(value) || value <= 0 || !nzchar(pair[1]))
                stop("Invalid --panel_widths_mm entry: ", entry, "; expected prefix:mm (> 0).")
            overrides[pair[1]] = value
        }
        # Optional panels may be absent for a particular family.
    }
    minimum_width = function(name, default) {
        value = default
        for (key in names(overrides)) {
            if (startsWith(name, key)) value = max(value, overrides[[key]])
        }
        value
    }
    domain_width = minimum_width('domain', 22.5)
    if ('domain' %in% names(g) && !is.null(attr(g[['domain']], 'treevis_width_mm'))) {
        domain_width = max(domain_width, attr(g[['domain']], 'treevis_width_mm'))
    }
    annotation_titles = c(domain='Amino acid\nposition (aa)',
        alignment='Alignment\nposition (nt)', fimo='Promoter motif\nposition (kb)')
    for (name in intersect(names(annotation_titles), names(g))) {
        g[[name]] = g[[name]] + ggplot2::labs(x=annotation_titles[[name]]) +
            ggplot2::guides(x=ggplot2::guide_axis(check.overlap=TRUE))
    }
    for (name in intersect(c('domain', 'fimo'), names(g))) {
        g[[name]] = treevis_fit_fill_legend(g[[name]],
            minimum_width(name, domain_width) - 2, name)
    }
    # Thin annotation columns use vertical, single-line axis titles.
    for (name in names(g)) {
        if (!grepl('^(localization$|signal_peptide$|peroxisome$|tm$|intron$|categorical,|pairwise_tip_suffix$)', name)) next
        label = g[[name]]$labels$x
        if (is.character(label)) g[[name]] = g[[name]] + ggplot2::labs(x=gsub('\n', ' ', label, fixed=TRUE))
        g[[name]] = g[[name]] + ggplot2::theme(
            axis.title.x=ggplot2::element_text(angle=90, hjust=1, vjust=0.5))
    }
    built_plots = lapply(g, ggplot2::ggplot_build)
    grobs = lapply(built_plots, ggplot2::ggplot_gtable)
    grobs = cowplot::align_plots(plotlist=grobs, align='h', axis='bt')
    names(grobs) = names(g)
    if (is.null(height_mm)) {
        n_tip = max(vapply(g, function(p) length(unique(p$data$label)), integer(1)))
        height_mm = (max(3, n_tip / 10) + if ('synteny' %in% names(g)) 1.4 else 0) * 25.4
    }
    vertical_margin = vapply(grobs, function(gt) {
        panel = gt$layout[gt$layout$name == 'panel', , drop=FALSE]
        outside = setdiff(seq_along(gt$heights), seq.int(panel$t[1], panel$b[1]))
        sum(grid::convertHeight(gt$heights[outside], 'mm', valueOnly=TRUE))
    }, numeric(1))
    panel_height_mm = height_mm - max(vertical_margin)

    widths = panel_widths = setNames(numeric(length(g)), names(g))
    for (name in names(g)) {
        p = g[[name]]
        built = built_plots[[name]]
        gt = grobs[[name]]
        if (name %in% c('domain', 'alignment', 'fimo')) {
            # Endpoint tick text must not spill into the next compact column.
            for (i in which(gt$layout$name == 'axis-b')) {
                gt$grobs[[i]] = treevis_keep_edge_ticks_inside(gt$grobs[[i]])
            }
        }
        panel_cells = gt$layout[grepl('^panel($|-)', gt$layout$name), , drop=FALSE]
        if (nrow(panel_cells) != 1) stop('Expected one data panel in column: ', name)
        panel_cols = seq.int(panel_cells$l[1], panel_cells$r[1])
        needed = if (grepl('^tree', name)) {
            60
        } else if (name %in% c('domain', 'fimo')) {
            domain_width
        } else if (name %in% c('alignment', 'synteny')) {
            22.5
        } else if (name == 'gene_structure') {
            23
        } else if (name == 'intron_sites') {
            45
        } else if (name == 'pointplot') {
            30
        } else if (grepl('^(localization$|signal_peptide$|peroxisome$|tm$|intron$|categorical,)', name)) {
            8
        } else {
            20
        }
        square_bar_height = attr(p, 'treevis_square_bar_height')
        if (!is.null(square_bar_height)) {
            if (panel_height_mm <= 0) stop('Figure height leaves no space for localization bars.')
            ranges = built$layout$panel_params[[1]]
            item_width = attr(p, 'treevis_square_item_width')
            if (is.null(item_width)) item_width = 1
            needed = panel_height_mm * square_bar_height / diff(ranges$y.range) * diff(ranges$x.range) / item_width
        }
        fixed = attr(p, 'treevis_width_mm')
        if (!is.null(fixed)) needed = fixed
        if (grepl('^heatmap($|,)', name)) {
            needed = max(8, length(unique(p$data$group)) * 4)
        }
        if (name == 'alignment') {
            xs = unlist(lapply(built$data, function(d) d$x))
            needed = max(needed, length(unique(xs[is.finite(xs)])) * 0.15)
        }
        if (grepl('^(amino_acid_site|site_state)', name)) {
            xs = unlist(lapply(built$data, function(d) d$x))
            needed = max(8, length(unique(xs[is.finite(xs)])) * 3)
        }
        # Text inside the data panel is not included in gtable's fixed widths.
        if (grepl('^(tiplabel$|text,|ortholog,|pairwise_tip_suffix$|tm$|intron$)', name)) {
            text_width = 0
            for (i in seq_along(p$layers)) {
                if (!inherits(p$layers[[i]]$geom, 'GeomText')) next
                d = built$data[[i]]
                for (j in seq_len(nrow(d))) {
                    if (is.na(d$label[j])) next
                    label = grid::textGrob(d$label[j], gp=grid::gpar(
                        fontsize=d$size[j] * ggplot2::.pt,
                        fontfamily=d$family[j], fontface=d$fontface[j],
                        lineheight=d$lineheight[j]))
                    text_width = max(text_width, grid::convertWidth(grid::grobWidth(label), 'mm', valueOnly=TRUE))
                }
            }
            needed = max(if (grepl('^(tiplabel$|text,|ortholog,)', name)) 2 else needed,
                         text_width * 1.1 + 2)
        }
        for (key in names(overrides)) {
            if (startsWith(name, key)) needed = max(needed, overrides[[key]])
        }
        outside = setdiff(seq_along(gt$widths), panel_cols)
        outside_mm = if (length(outside)) sum(grid::convertWidth(gt$widths[outside], 'mm', valueOnly=TRUE)) else 0
        data_width = needed
        # Bottom/top titles and legends can span the panel or the whole column.
        # Reserve their intrinsic width rather than letting them spill into neighbors.
        for (i in seq_along(gt$grobs)) {
            if (!grepl('^(guide-box|xlab|title|subtitle|caption)', gt$layout$name[i])) next
            child = gt$grobs[[i]]
            child_mm = grid::convertWidth(grid::grobWidth(child), 'mm', valueOnly=TRUE)
            if (!is.null(square_bar_height) && inherits(child, 'titleGrob')) {
                # titleGrob itself has a null width; measure its text children.
                child_mm = max(c(child_mm, vapply(child$children, function(text) {
                    grid::convertWidth(grid::grobWidth(text), 'mm', valueOnly=TRUE)
                }, numeric(1))))
            }
            span = seq.int(gt$layout$l[i], gt$layout$r[i])
            extra = intersect(span, outside)
            extra_mm = if (length(extra)) sum(grid::convertWidth(gt$widths[extra], 'mm', valueOnly=TRUE)) else 0
            needed = max(needed, child_mm - extra_mm + 2)
        }
        if (!is.null(square_bar_height)) {
            # Keep each probability bar square; let the legend/title use a
            # separate gutter instead of stretching the data panel.
            gutter = needed - data_width
            edges = c(1L, length(gt$widths))
            gt$widths[edges] = gt$widths[edges] + grid::unit(gutter / 2, 'mm')
            outside_mm = outside_mm + gutter
            labels = grepl('^(guide-box|xlab|title|subtitle|caption)', gt$layout$name)
            gt$layout$l[labels] = 1L
            gt$layout$r[labels] = length(gt$widths)
            needed = data_width
        }
        gt$widths[panel_cols] = grid::unit(rep(needed / length(panel_cols), length(panel_cols)), 'mm')
        panel_widths[name] = needed
        widths[name] = needed + outside_mm
        grobs[[name]] = gt
    }
    list(plots=grobs, widths_mm=widths, panel_widths_mm=panel_widths,
         width_mm=sum(widths),
         legend_entries=lapply(g, function(p) attr(p, 'treevis_legend_visible')))
}
