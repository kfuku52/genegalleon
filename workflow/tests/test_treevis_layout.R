suppressPackageStartupMessages(library(genegalleon.treevis))
suppressPackageStartupMessages(library(ggplot2))
output = Sys.getenv('TREEVIS_LAYOUT_OUTPUT', unset=tempfile('treevis-layout-'))
dir.create(output, showWarnings=FALSE, recursive=TRUE)
pdf(file.path(output, 'measure.pdf'))
base = ggplot(data.frame(x=1:3,y=1:3),aes(x,y)) + geom_point() + theme_classic()
g = list(tree=base, domain=base)
a = treevis_layout_mm(g)
b = treevis_layout_mm(c(g,list(extra=base)))
stopifnot(identical(a$widths_mm,b$widths_mm[names(g)]),
          b$width_mm > a$width_mm,
          a$panel_widths_mm['tree'] >= 60,
          a$panel_widths_mm['domain'] == 22.5)
# Overrides may increase minima, but may never squeeze a column.
c = treevis_layout_mm(g, 'tree:1,domain:90')
stopifnot(c$panel_widths_mm['tree'] == a$panel_widths_mm['tree'],
          c$panel_widths_mm['domain'] >= 90)
for (bad in c('tree:0','tree:-1','tree:NaN','tree:Inf','tree','tree:1:2')) {
    stopifnot(inherits(try(treevis_layout_mm(g,bad),silent=TRUE),'try-error'))
}
text_plot = function(label, size=3) ggplot(data.frame(x=0,y=1,label=label),aes(x,y,label=label)) +
    geom_text(size=size,hjust=0) + scale_x_continuous(limits=c(0,1),expand=expansion(mult=0)) + theme_void()
w = function(p) treevis_layout_mm(list(tiplabel=p))$panel_widths_mm[[1]]
stopifnot(w(text_plot('WWWW')) > w(text_plot('iiii')),
          w(text_plot('label',6)) > w(text_plot('label',3)))
heat = function(n) ggplot(expand.grid(group=letters[seq_len(n)],y=1:3),aes(group,y)) + geom_tile() + theme_void()
stopifnot(treevis_layout_mm(list(heatmap=heat(10)))$panel_widths_mm >= 40,
          treevis_layout_mm(list(heatmap=heat(10)))$width_mm > treevis_layout_mm(list(heatmap=heat(2)))$width_mm)
# Check actual fixed panel units and the physical PDF page size.
for (name in names(a$plots)) {
    gt = a$plots[[name]]
    panel = gt$layout[gt$layout$name == 'panel',]
    actual = grid::convertWidth(sum(gt$widths[panel$l:panel$r]),'mm',valueOnly=TRUE)
    stopifnot(abs(actual-a$panel_widths_mm[name]) < 1e-8)
}
dev.off()
plot = cowplot::plot_grid(plotlist=a$plots,nrow=1,align='h',axis='bt',rel_widths=a$widths_mm)
file = file.path(output,'layout.pdf')
ggsave(file,plot,width=a$width_mm,height=80,units='mm',limitsize=FALSE)
read_width = "import re,sys; data=open(sys.argv[1],'rb').read(); m=re.search(rb'/MediaBox\\s*\\[\\s*0\\s+0\\s+([0-9.]+)',data); print(float(m.group(1)))"
points = as.numeric(system2('python', c('-c',shQuote(read_width),shQuote(file)),stdout=TRUE))
# R's PDF MediaBox rounds to whole points.
stopifnot(abs(points * 25.4/72 - a$width_mm) < 25.4/72)
cat('Physical layout tests passed. Output:',file,'\n')

# Compact annotation columns and independent legend filtering.
pdf(NULL)
keys <- sprintf('motif%02d', 1:30)
annotation <- ggplot(data.frame(x=seq_along(keys), y=1, key=factor(keys, levels=keys)),
    aes(x,y,fill=key)) + geom_tile() +
    scale_fill_manual(values=setNames(grDevices::rainbow(length(keys)),keys)) +
    theme_minimal(base_size=6) + theme(legend.position='bottom')
attr(annotation, 'treevis_legend_counts') <- setNames(30:1, keys)
fit_legend <- getFromNamespace('treevis_fit_fill_legend', 'genegalleon.treevis')
fitted <- fit_legend(annotation,20.5,'test')
visible <- attr(fitted,'treevis_legend_visible')
stopifnot(length(visible)>0, length(visible)<length(keys),
          identical(visible, keys[seq_along(visible)]),
          identical(ggplot_build(annotation)$data[[1]]$fill,
                    ggplot_build(fitted)$data[[1]]$fill),
          nrow(ggplot_build(annotation)$data[[1]]) == nrow(ggplot_build(fitted)$data[[1]]),
          grid::convertWidth(grid::grobWidth(cowplot::get_legend(fitted)), 'mm', valueOnly=TRUE) <= 20.5)
compact <- treevis_layout_mm(list(domain=annotation, fimo=annotation, alignment=base, synteny=base))
stopifnot(all(compact$panel_widths_mm == 22.5))
wider <- treevis_layout_mm(list(domain=annotation, fimo=annotation), 'domain:50')
stopifnot(all(wider$panel_widths_mm == 50),
          length(wider$legend_entries$fimo) >= length(compact$legend_entries$fimo))
# A query tile's 0.9 x 0.9 geometry must remain square in physical units.
input <- list(tree=list(data=data.frame(isTip=TRUE,label=paste0('g',1:10),y=1:10,
                                       query_marker=rep(c('Best hit','-'),5))))
q <- add_categorical_column(input, list(font_size=6,margins=rep(0,4)),
    'categorical,query_marker,Query','query_marker','Query')
q$tree <- NULL
ql <- treevis_layout_mm(q, height_mm=100)
gt <- ql$plots[[1]]
cell <- gt$layout[gt$layout$name == 'panel',]
outside <- setdiff(seq_along(gt$heights),seq.int(cell$t,cell$b))
y_mm <- 100 - sum(grid::convertHeight(gt$heights[outside],'mm',valueOnly=TRUE))
ranges <- ggplot_build(q[[1]])$layout$panel_params[[1]]
stopifnot(abs(ql$panel_widths_mm[[1]]*0.9/diff(ranges$x.range) - y_mm*0.9/diff(ranges$y.range)) < 1e-8)
dev.off()
cat('Compact panels, query squares, and frequency-limited legends passed.\n')
