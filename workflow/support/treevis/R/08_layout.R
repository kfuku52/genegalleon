# Physical widths are resolved on a PDF device, before combining the columns.
# Overrides are lower bounds for the data panel, never fractions of the page.
treevis_layout_mm = function(g, panel_widths_mm = NULL) {
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
    grobs = vector('list', length(g))
    names(grobs) = names(g)
    widths = panel_widths = setNames(numeric(length(g)), names(g))
    for (name in names(g)) {
        p = g[[name]]
        built = ggplot2::ggplot_build(p)
        gt = ggplot2::ggplot_gtable(built)
        panel_cells = gt$layout[grepl('^panel($|-)', gt$layout$name), , drop=FALSE]
        if (nrow(panel_cells) != 1) stop('Expected one data panel in column: ', name)
        panel_cols = seq.int(panel_cells$l[1], panel_cells$r[1])
        needed = if (grepl('^tree', name)) {
            60
        } else if (grepl('^(domain|alignment|fimo|synteny|gene_structure|intron_sites)$', name)) {
            45
        } else if (name == 'pointplot') {
            30
        } else if (grepl('^(signal_peptide$|peroxisome$|tm$|intron$|categorical,)', name)) {
            8
        } else {
            20
        }
        fixed = attr(p, 'treevis_width_mm')
        if (!is.null(fixed)) needed = fixed
        if (grepl('^heatmap($|,)', name)) {
            needed = max(8, length(unique(p$data$group)) * 4)
        }
        if (name == 'alignment') {
            xs = unlist(lapply(built$data, function(d) d$x))
            needed = max(needed, length(unique(xs[is.finite(xs)])) * 0.3)
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
        # Bottom/top titles and legends can span the panel or the whole column.
        # Reserve their intrinsic width rather than letting them spill into neighbors.
        for (i in seq_along(gt$grobs)) {
            if (!grepl('^(guide-box|xlab|title|subtitle|caption)', gt$layout$name[i])) next
            child = gt$grobs[[i]]
            child_mm = grid::convertWidth(grid::grobWidth(child), 'mm', valueOnly=TRUE)
            span = seq.int(gt$layout$l[i], gt$layout$r[i])
            extra = intersect(span, outside)
            extra_mm = if (length(extra)) sum(grid::convertWidth(gt$widths[extra], 'mm', valueOnly=TRUE)) else 0
            needed = max(needed, child_mm - extra_mm + 2)
        }
        gt$widths[panel_cols] = grid::unit(rep(needed / length(panel_cols), length(panel_cols)), 'mm')
        panel_widths[name] = needed
        widths[name] = needed + outside_mm
        grobs[[name]] = gt
    }
    list(plots=grobs, widths_mm=widths, panel_widths_mm=panel_widths,
         width_mm=sum(widths))
}
