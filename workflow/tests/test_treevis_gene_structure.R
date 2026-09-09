suppressPackageStartupMessages(library(genegalleon.treevis))
suppressPackageStartupMessages(library(ggplot2))
tips = data.frame(label=c('plus','minus','zero','missing'), y=1:4,
    feature_blocks=c('1-30;1031-1090','1061-1090;1-60','1-90',NA),
    feature_type='CDS', strand=c('+','-','+','+'), num_intron=c(1,1,0,NA),
    intron_feature_size=c(90,90,90,NA),
    tiplab_color=c('red','blue','green','black'))
make_data = genegalleon.treevis:::treevis_gene_structure_data
d = make_data(tips, 'linear')
stopifnot(nrow(d$boxes)==5, nrow(d$introns)==2,
          all(d$introns$length_bp==1000),
          identical(d$boxes$start,c(0,1030,0,1030,0)),
          identical(d$boxes$colour,c('red','red','blue','blue','green')))
compressed = make_data(tips, 'compressed')
stopifnot(all(compressed$introns$end-compressed$introns$start < 1000),
          all(compressed$boxes$end-compressed$boxes$start == d$boxes$end-d$boxes$start))
bad = tips; bad$feature_blocks[1]='1-30;20-90'
stopifnot(inherits(try(make_data(bad),silent=TRUE),'try-error'))
bad = tips; bad$num_intron[1]=2
stopifnot(inherits(try(make_data(bad),silent=TRUE),'try-error'))
tree = ggtree::ggtree(ape::read.tree(text='((plus:1,minus:1):1,zero:2,missing:2);'))
ii = match(tree$data$label,tips$label)
for (key in setdiff(names(tips), c('label','y'))) tree$data[[key]] = tips[[key]][ii]
g = list(tree=tree)
args = list(font_size=6, margins=rep(0.03,4))
p = add_gene_structure_column(g,args)$gene_structure
stopifnot(attr(p,'treevis_width_mm')==23,
          all(ggplot_build(p)$data[[1]]$fill == p$data$colour))
missing = g; missing$tree$data$feature_blocks=NA_character_
stopifnot(identical(add_gene_structure_column(missing,args),missing))
f = tempfile(fileext='.tsv')
write.table(data.frame(node_name=c('plus','minus'),offset=c(1,1),group_id='x'),
            f,sep='\t',quote=FALSE,row.names=FALSE)
suffix = add_synteny_similarity_column(g,args,f)$pairwise_tip_suffix
expected = get_df_tip(tree)$tiplab_color
stopifnot(identical(ggplot_build(suffix)$data[[1]]$colour, expected))
cat('Gene structure coordinates, missing data, compression and shared tip colors passed.\n')

utr_tips = tips
utr_tips$utr_blocks = c('1091-1100','1091-1100',NA,NA)
utr_data = make_data(utr_tips,'linear')
stopifnot(sum(utr_data$boxes$feature=='UTR')==2,
          all(utr_data$boxes$fill[utr_data$boxes$feature=='UTR'] != utr_data$boxes$colour[utr_data$boxes$feature=='UTR']),
          nrow(utr_data$introns)==2,
          utr_data$boxes$feature[utr_data$boxes$label=='minus'][1]=='UTR')

site_tips = data.frame(label=c('a','b','phase','gap','zero','bad_length','missing'),y=1:7,
    num_intron=c(1,1,1,1,0,1,NA), intron_positions=c('3','3','3','3','','3',NA),
    intron_feature_size=c(12,12,12,9,12,9,NA),cds_first_phase=c(0,0,1,0,0,0,NA),
    tiplab_color='red')
seqs = setNames(rep('AAACCCGGGTTT',7),site_tips$label)
seqs['gap']='AAA---GGGTTT'
sites = genegalleon.treevis:::treevis_intron_site_data(site_tips, seqs)
stopifnot(nrow(sites$sites)==2,
    sites$events$site_id[sites$events$node_name=='a']==sites$events$site_id[sites$events$node_name=='b'],
    sites$events$site_id[sites$events$node_name=='a']!=sites$events$site_id[sites$events$node_name=='phase'],
    is.na(sites$events$site_id[sites$events$node_name=='gap']),
    sites$events$status[sites$events$node_name=='gap']=='alignment_gap',
    'CDS_length_mismatch' %in% sites$diagnostics$reason,
    'missing_GFF' %in% sites$diagnostics$reason,
    sites$cells$state[sites$cells$node_name=='zero' & sites$cells$site_id=='I001']=='absent')
cat('CDS/UTR encoding, phase-specific correspondence and ambiguous junctions passed.\n')
unknown = site_tips
unknown$cds_first_phase[unknown$label=='b'] = NA
unknown_sites = genegalleon.treevis:::treevis_intron_site_data(unknown,seqs)
stopifnot(unknown_sites$events$status[unknown_sites$events$node_name=='b']=='position_only',
    all(unknown_sites$cells$state[unknown_sites$cells$node_name=='b']=='position_only'))

# Inline IDs attach by CDS offset rather than genomic/UTR intron order.
inline_introns = data.frame(label=c('a','a','b','gap'),cds_offset=c(0,3,3,3),
    y=c(1,1,2,4),start=c(0,50,60,70),end=c(10,60,70,80))
inline = genegalleon.treevis:::treevis_structure_intron_labels(inline_introns,unknown_sites$events)
stopifnot(nrow(inline)==3, !any(inline$cds_offset==0),
    inline$display_id[inline$label=='a']=='1',
    inline$display_id[inline$label=='gap']=='?',
    inline$display_id[inline$label=='b']=='?') # Unknown phase matches two phase groups.

# Correspondence is embedded in the structure plot, without adding another panel.
inline_tree = g$tree
inline_tree$data$cds_first_phase = 0
inline_tree$data$intron_positions = c('30','30','',NA)[match(inline_tree$data$label,tips$label)]
aln = tempfile(fileext='.fa')
writeLines(as.vector(rbind(paste0('>',tips$label),rep(paste(rep('A',90),collapse=''),4))),aln)
embedded = add_gene_structure_column(list(tree=inline_tree),args,path_alignment=aln)
stopifnot(identical(names(embedded),c('tree','gene_structure')))
labels = attr(embedded$gene_structure,'treevis_intron_labels')
stopifnot(nrow(labels)==2,all(labels$display_id=='1'),
    attr(embedded$gene_structure,'treevis_width_mm')==23,
    !any(vapply(embedded$gene_structure$layers, function(layer) inherits(layer$geom,'GeomText'), logical(1))),
    nrow(attr(embedded$gene_structure,'treevis_intron_connections'))>0)

# Bands connect actual introns in occurrence order. Unknown phase interrupts
# exact-match bands and changes the adjoining edges to weak bands.
band_labels = data.frame(site_id=c('I001','I001','I001','I001','I002'),
    status=c('mapped','mapped','position_only','mapped','mapped'),
    alignment_left=c(100,100,100,100,200),phase=c(0,0,NA,0,0),
    x=c(10,20,30,40,50),text_x=c(100,200,300,400,500),
    start=c(9,19,29,39,49),end=c(11,21,31,41,51),y=c(1,2,3,4,5))
bands = genegalleon.treevis:::treevis_intron_connection_polygons(band_labels,0.5)
stopifnot(nrow(bands)==12, length(unique(bands$connection))==3,
    bands$kind[1]=='exact_phase', all(bands$kind[5:12]=='position_candidate'),
    identical(bands$x[1:4],c(9.5,10.5,20.5,19.5)),
    identical(bands$y[5:8],c(2,2,3,3)),
    !any(bands$from_index==2 & bands$to_index==4))


near_labels = data.frame(site_id=c('I001',NA,'I002','I003','I004'),
    status=c('mapped','position_only','mapped','mapped','mapped'),
    alignment_left=c(100,101,102,106,106),phase=c(0,NA,2,0,1),
    x=1:5,start=(1:5)-0.5,end=(1:5)+0.5,y=1:5)
near_bands = genegalleon.treevis:::treevis_intron_connection_polygons(near_labels,0.3)
stopifnot(nrow(near_bands)==12,all(near_bands$kind=='position_candidate'),
    all(near_bands$opacity==0.05),!any(near_bands$from_index==3 & near_bands$to_index==4))
# Equidistant alternatives do not arbitrarily select a nearby intron.
tied = rbind(near_labels[1:2,],transform(near_labels[2,],alignment_left=99,x=3,start=2.5,end=3.5))
stopifnot(nrow(genegalleon.treevis:::treevis_intron_connection_polygons(tied,0.3))==0)

# A complete intron-free family and an all-missing family are valid empty results.
empty_tips = data.frame(label='a',y=1,num_intron=0,intron_feature_size=6,
    cds_first_phase=0,intron_positions='')
empty_sites = genegalleon.treevis:::treevis_intron_site_data(empty_tips,c(a='AAAAAA'))
stopifnot(nrow(empty_sites$events)==0,nrow(empty_sites$sites)==0,nrow(empty_sites$diagnostics)==0)
empty_tips$num_intron=NA
empty_sites = genegalleon.treevis:::treevis_intron_site_data(empty_tips,c(a='AAAAAA'))
stopifnot(nrow(empty_sites$events)==0,nrow(empty_sites$diagnostics)==1)
adjacent = tips[3,,drop=FALSE]; adjacent$feature_blocks='1-30;31-90'
stopifnot(nrow(make_data(adjacent)$introns)==0,nrow(make_data(adjacent)$boxes)==2)

# A wrong GFF transcript must not survive as a plausible zero/short structure.
bad_tree=inline_tree
bad_tree$data$intron_feature_size[bad_tree$data$label=='zero']=3
bad_tree$data$feature_blocks[bad_tree$data$label=='zero']='1-3'
stopifnot(inherits(try(add_gene_structure_column(list(tree=bad_tree),args,path_alignment=aln),silent=TRUE),'try-error'))
