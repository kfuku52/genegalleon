suppressPackageStartupMessages(library(genegalleon.treevis))
suppressPackageStartupMessages(library(ggplot2))
# Grob measurement must not create Rplots.pdf in the working directory.
pdf(NULL)

outdir <- commandArgs(trailingOnly=TRUE)
if (!length(outdir)) outdir <- tempfile('synteny-legend-')
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
outdir <- normalizePath(outdir)
tree <- ape::read.tree(text='((g1:1,g2:1):1,g3:2);')
g <- list(tree=ggtree::ggtree(tree))
args <- list(font_size=8, margins=c(0.1, 0.1, 0.1, 0.1))
syn <- data.frame(
    node_name=rep(c('g1', 'g2', 'g3'), each=5),
    offset=c(-3,-2,-1,1,2, -3,-2,-1,1,3, -3,-2,-1,2,3),
    group_id=c('A','U1','B','C','U2', 'A','U3','B','C','U4', 'U5','A','B','C','U6'),
    evalue_cutoff=1e-10)
infile <- file.path(outdir, 'synteny.tsv')
write.table(syn, infile, sep='\t', quote=FALSE, row.names=FALSE)
plot <- add_synteny_column(g, args, 'synteny', infile, 3)$synteny
stopifnot(plot$theme$legend.position == 'bottom',
          is.null(plot$scales$get_scales('shape')$name),
          'Same gene family\n(E-value <= 1e-10)' %in% plot$scales$get_scales('shape')$breaks)
built <- ggplot_build(plot)
# Semantic legend entries never add plotted points or alter neighbor occupancy.
stopifnot(nrow(built$data[[4]]) == nrow(syn), nrow(built$data[[5]]) == 3,
          nrow(built$data[[6]]) == 9)
# A panel containing only adjacent links must retain their opacity, without
# rescaling it to the range of distances present in this particular tree.
stopifnot(all(built$data[[3]]$alpha > 0.5),
          all(built$data[[4]]$alpha == 1), all(built$data[[5]]$alpha == 1),
          all(built$data[[6]]$alpha == 1))
legend <- cowplot::get_legend(plot)
stopifnot(inherits(legend, 'gtable'))
combined <- cowplot::plot_grid(g$tree + ggtree::geom_tiplab() + xlim(0, 2.25), plot,
                              nrow=1, align='h', axis='bt', rel_widths=c(1.5,1))
ggsave(file.path(outdir, 'neighboring-genes-legend.pdf'), combined, width=6, height=4.5)
ggsave(file.path(outdir, 'neighboring-genes-legend.png'), combined, width=6, height=4.5, dpi=180, bg='white')
# Legacy inputs must not imply a cutoff, and conflicting metadata must fail.
syn$evalue_cutoff <- NULL
write.table(syn, infile, sep='\t', quote=FALSE, row.names=FALSE)
legacy <- add_synteny_column(g, args, 'synteny', infile, 3)$synteny
stopifnot(any(grepl('unavailable', legacy$scales$get_scales('shape')$breaks)))
syn$evalue_cutoff <- c(0.01, rep(1e-10, nrow(syn)-1))
write.table(syn, infile, sep='\t', quote=FALSE, row.names=FALSE)
stopifnot(inherits(try(add_synteny_column(g, args, 'synteny', infile, 3), silent=TRUE), 'try-error'))
syn$evalue_cutoff <- 1e-10
write.table(syn, infile, sep='\t', quote=FALSE, row.names=FALSE)
cat('Synteny legend tests passed. Example:', outdir, '\n')

# Cutoff metadata must be valid and retain meaningful precision in its label.
for (cutoff in c(0, 1e-100, 1.23456789e-10)) {
    variant <- syn
    variant$evalue_cutoff <- cutoff
    write.table(variant, infile, sep='\t', quote=FALSE, row.names=FALSE)
    panel <- add_synteny_column(g,args,'synteny',infile,3)$synteny
    expected <- paste0('Same gene family\n(E-value <= ',
                       format(cutoff, scientific=TRUE, trim=TRUE, digits=15), ')')
    stopifnot(expected %in% panel$scales$get_scales('shape')$breaks)
}
for (cutoff in c(-1, NA_real_, Inf, NaN)) {
    variant <- syn
    variant$evalue_cutoff <- cutoff
    write.table(variant, infile, sep='\t', quote=FALSE, row.names=FALSE)
    stopifnot(inherits(try(add_synteny_column(g,args,'synteny',infile,3),silent=TRUE),'try-error'))
}

# Missing, empty, unmatched and unshared inputs keep the panel absent.
stopifnot(identical(add_synteny_column(g,args,'synteny',paste0(infile,'.absent'),3),g))
for (kind in c('empty','unmatched','unshared','outside')) {
    variant <- syn
    if (kind == 'empty') variant <- variant[FALSE,]
    if (kind == 'unmatched') variant$node_name <- 'absent'
    if (kind == 'unshared') variant$group_id <- seq_len(nrow(variant))
    if (kind == 'outside') variant$offset <- 20
    write.table(variant, infile, sep='\t', quote=FALSE, row.names=FALSE)
    stopifnot(identical(add_synteny_column(g,args,'synteny',infile,3),g))
}

# A repeated family on one side retains its nearest member; its other
# recorded position remains occupied, and links never cross focal genes.
variant <- rbind(syn,data.frame(node_name='g1',offset=-2,group_id='B',evalue_cutoff=1e-10))
write.table(variant,infile,sep='\t',quote=FALSE,row.names=FALSE)
panel <- add_synteny_column(g,args,'synteny',infile,3)$synteny
stopifnot(nrow(panel$data[panel$data$node_name=='g1' & panel$data$group_id=='B',]) == 1,
          panel$data$offset[panel$data$node_name=='g1' & panel$data$group_id=='B'] == -1)
links <- ggplot_build(panel)$data[[3]]
stopifnot(all((links$x-4)*(links$xend-4) > 0))

# With shared groups only on opposing sides, an empty link layer must render.
variant <- data.frame(node_name=c('g1','g2'),offset=c(-1,1),group_id='A',evalue_cutoff=1e-10)
write.table(variant,infile,sep='\t',quote=FALSE,row.names=FALSE)
panel <- add_synteny_column(g,args,'synteny',infile,1)$synteny
stopifnot(inherits(ggplotGrob(panel),'gtable'), nrow(panel$data)==2)
write.table(syn,infile,sep='\t',quote=FALSE,row.names=FALSE)

# Opacity follows displayed vertical separation, including within one family.
# Horizontal shifts must not affect opacity, and faint links render underneath.
distance_tree <- ape::read.tree(text=paste0('(', paste0('t', 1:13, ':1', collapse=','), ');'))
distance_g <- list(tree=ggtree::ggtree(distance_tree))
tip_rows <- get_df_tip(distance_g$tree)
distance_syn <- data.frame(
    node_name=as.character(tip_rows$label[c(1,2,6,1,3,11,4,5)]),
    offset=c(-1,-1,-2,1,2,3,-3,-1),
    group_id=c('A','A','A','B','B','B','C','C'), evalue_cutoff=1e-10)
distance_infile <- file.path(outdir, 'synteny-distance.tsv')
write.table(distance_syn,distance_infile,sep='\t',quote=FALSE,row.names=FALSE)
distance_panel <- add_synteny_column(distance_g,args,'synteny',distance_infile,3)$synteny
distance_built <- ggplot_build(distance_panel)
distance_links <- distance_built$data[[3]]
gaps <- abs(distance_links$yend-distance_links$y)
opacity <- vapply(sort(unique(gaps)), function(gap) unique(distance_links$alpha[gaps == gap]), numeric(1))
stopifnot(identical(sort(unique(gaps)), c(1,2,4,8)),
          all(diff(opacity) < 0), opacity[1] > 0.5, tail(opacity,1) < 0.2,
          all(is.finite(distance_links$alpha)), all(distance_links$alpha > 0),
          all(distance_links$alpha <= 1), all(diff(distance_links$alpha) >= 0),
          all(distance_built$data[[6]]$alpha == 1),
          nrow(distance_panel$data) == nrow(distance_syn))
# Batched marker/link assembly keeps the same geometry with shuffled rows,
# duplicate records and farther copies of a family on the same side.
copies <- distance_syn[abs(distance_syn$offset) < 3,]
copies$offset <- copies$offset + sign(copies$offset)
variant <- rbind(distance_syn, distance_syn, copies)
variant$neighbor_gene <- paste0('neighbor',seq_len(nrow(variant)))
variant <- variant[rev(seq_len(nrow(variant))),]
write.table(variant,distance_infile,sep='\t',quote=FALSE,row.names=FALSE)
batched_panel <- add_synteny_column(distance_g,args,'synteny',distance_infile,3)$synteny
stopifnot(identical(distance_panel$data,batched_panel$data),
          identical(distance_links,ggplot_build(batched_panel)$data[[3]]))
write.table(distance_syn,distance_infile,sep='\t',quote=FALSE,row.names=FALSE)
distance_example <- cowplot::plot_grid(
    distance_g$tree + ggtree::geom_tiplab(), distance_panel,
    nrow=1, align='h', axis='bt', rel_widths=c(1.5,1))
ggsave(file.path(outdir,'neighboring-genes-distance.png'), distance_example,
       width=6,height=6,dpi=180,bg='white')
cat('Synteny link opacity decreases with displayed row distance.\n')

# Exercise the actual workflow plot driver, including its reserved legend height.
test_file <- sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE)[1])
repo <- normalizePath(file.path(dirname(test_file), '..', '..'))
branch <- data.frame(
    branch_id=c(4,2,0,1,3), parent=c(-999,4,2,2,4),
    sister=c(-999,3,1,0,2), child1=c(2,0,NA,NA,NA), child2=c(3,1,NA,NA,NA),
    node_name=c('root','n4','g1','g2','g3'), bl_rooted=c(0,1,1,1,2),
    so_event=c('S','S','L','L','L'), so_event_parent='S')
write.table(branch, file.path(outdir, 'stat_branch.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
oldwd <- setwd(outdir)
result <- system2('Rscript', c(file.path(repo, 'workflow/support/stat_branch2tree_plot.r'),
    '--stat_branch=stat_branch.tsv', '--panel_widths_mm=tree:60',
    '--panel1=tree,bl_rooted,no,no,L', paste0('--panel2=synteny,', infile, ',3'),
    '--show_branch_id=no', '--event_method=species_overlap',
    '--species_color_table=PLACEHOLDER', '--pie_chart_value_transformation=identity',
    '--max_delta_intron_present=-0.5', '--long_branch_display=no'), stdout=TRUE, stderr=TRUE)
setwd(oldwd)
if (!is.null(attr(result, 'status'))) stop(paste(result, collapse='\n'))
stopifnot(file.info(file.path(outdir, 'stat_branch2tree_plot.pdf'))$size > 1000)
cat('Workflow tree plot with synteny legend passed.\n')


# Keep the reusable example input consistent with the rendered example.
write.table(syn,infile,sep='\t',quote=FALSE,row.names=FALSE)
dev.off()
