suppressPackageStartupMessages(library(genegalleon.treevis))
suppressPackageStartupMessages(library(ggplot2))
tree <- ape::read.tree(text='((g1:1,g2:1):1,g3:2);')
g <- list(tree=ggtree::ggtree(tree))
args <- list(font_size=8, margins=c(0.1, 0.1, 0.1, 0.1))
infile <- tempfile(fileext='.tsv')

# Hidden families must not consume palette colors, even if shared outside the
# display window or present on tips absent from the plotted tree.
visible <- data.frame(node_name=c('g1','g2','g1','g2'),
                      offset=c(-1,-1,1,1), group_id=c('B','B','D','D'),
                      evalue_cutoff=1e-10)
hidden <- data.frame(node_name=c('g1','g2','absent1','absent2','g1','g1','g2'),
                     offset=c(-20,-20,1,1,2,3,20),
                     group_id=c('A','A','C','C','E','F','F'),
                     evalue_cutoff=1e-10)
palette_for <- function(rows) {
    write.table(rows, infile, sep='\t', quote=FALSE, row.names=FALSE)
    p <- add_synteny_column(g, args, 'synteny', infile, 5)$synteny
    stopifnot(setequal(as.character(p$data$group_id), c('B','D')),
              nrow(p$data) == 4)
    built <- ggplot_build(p)
    scale <- built$plot$scales$get_scales('colour')
    stopifnot(setequal(scale$get_limits(), c('B','D')))
    scale$map(c('B','D'))
}
stopifnot(identical(palette_for(visible), palette_for(rbind(visible, hidden))))
cat('Hidden synteny families do not consume palette colors.\n')

unlink(infile)
