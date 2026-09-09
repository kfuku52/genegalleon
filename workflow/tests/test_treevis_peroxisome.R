suppressPackageStartupMessages(library(genegalleon.treevis))
suppressPackageStartupMessages(library(ggplot2))
tree <- list(data=data.frame(isTip=c(TRUE, TRUE, TRUE, TRUE, FALSE),
    label=c('zero', 'half', 'one', 'missing', 'internal'), y=c(4, 2, 1, 3, 5),
    cdskit_localize_p_peroxisome=c(0, 0.5, 1, NA, 0.9)))
args <- list(font_size=6, margins=c(0, 0, 0, 0))
g <- add_peroxisome_column(list(tree=tree), args)
p <- g$peroxisome
stopifnot(identical(as.character(p$data$label), c('one', 'half', 'missing', 'zero')))
built <- ggplot_build(p)
# Full-scale background and actual probability are independent, never normalized.
stopifnot(nrow(built$data[[2]]) == 3L,
          all(built$data[[2]]$x == 1),
          identical(as.numeric(built$data[[3]]$x), c(1, 0.5, 0)),
          identical(p$scales$get_scales('x')$limits, c(0, 1)),
          is.na(p$data$probability[3]))
missing <- tree
missing$data$cdskit_localize_p_peroxisome <- NULL
stopifnot(is.null(add_peroxisome_column(list(tree=missing), args)$peroxisome))
for (invalid in c(-0.1, 1.1, Inf)) {
    bad <- tree
    bad$data$cdskit_localize_p_peroxisome[1] <- invalid
    stopifnot(inherits(try(add_peroxisome_column(list(tree=bad), args), silent=TRUE), 'try-error'))
}
cat('Peroxisome panel tests passed.\n')

# Physical bars stay square even when legends require a wider column.
for (n_tip in c(10L, 60L)) {
    data <- data.frame(isTip=TRUE, label=paste0('gene', seq_len(n_tip)),
        y=seq_len(n_tip), cdskit_localize_p_peroxisome=0.5)
    for (key in c('noTP', 'SP', 'mTP', 'cTP', 'lTP')) data[[paste0('cdskit_localize_p_', key)]] <- 0.2
    plots <- add_signal_peptide_column(list(tree=list(data=data)), args)
    plots <- add_peroxisome_column(plots, args)
    plots$tree <- NULL
    height_mm <- max(3, n_tip / 10) * 25.4
    pdf(NULL)
    layout <- treevis_layout_mm(plots, height_mm=height_mm)
    for (name in names(plots)) {
        gt <- layout$plots[[name]]
        cell <- gt$layout[gt$layout$name == 'panel', ]
        outside <- setdiff(seq_along(gt$heights), seq.int(cell$t, cell$b))
        panel_height <- height_mm - sum(grid::convertHeight(gt$heights[outside], 'mm', valueOnly=TRUE))
        ranges <- ggplot_build(plots[[name]])$layout$panel_params[[1]]
        bar_height <- panel_height * 0.8 / diff(ranges$y.range)
        bar_width <- layout$panel_widths_mm[[name]] / diff(ranges$x.range)
        stopifnot(abs(bar_width - bar_height) < 1e-8,
                  layout$widths_mm[[name]] >= layout$panel_widths_mm[[name]])
    }
    dev.off()
}
cat('Square localization bar layout tests passed.\n')

# A single localization column contains two independent adjacent squares.
input <- tree
for (key in c('noTP','SP','mTP','cTP','lTP')) input$data[[paste0('cdskit_localize_p_',key)]] <- 0.2
combined <- add_localization_column(list(tree=input),args)
stopifnot(identical(names(combined),c('tree','localization')))
b <- ggplot_build(combined$localization)
rectangles <- b$data[[3]]
left <- rectangles[rectangles$xmin < 1,]
right <- rectangles[rectangles$xmin > 1,]
stopifnot(all(abs(tapply(left$xmax-left$xmin,left$ymin,sum)-1)<1e-10),
          identical(as.numeric(right$xmax-right$xmin),c(1,0.5,0)),
          all(abs(right$xmin-1.4)<1e-10))
# The peroxisome prediction never rescales the five-way targeting stack.
input$data$cdskit_localize_p_peroxisome <- NA_real_
stopifnot(!is.null(add_localization_column(list(tree=input),args)$localization))
input$data$cdskit_localize_p_peroxisome <- NULL
stopifnot(!is.null(add_localization_column(list(tree=input),args)$localization))
only_perox <- add_localization_column(list(tree=tree),args)
stopifnot(!is.null(only_perox$localization))
for (key in c('noTP','SP','mTP','cTP','lTP')) input$data[[paste0('cdskit_localize_p_',key)]] <- NULL
stopifnot(is.null(add_localization_column(list(tree=input),args)$localization))
cat('Combined localization probabilities and missing-data tests passed.\n')
