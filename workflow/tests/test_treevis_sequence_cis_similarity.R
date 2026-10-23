suppressPackageStartupMessages(library(genegalleon.treevis))
suppressPackageStartupMessages(library(ggplot2))
similarity_source = Sys.getenv('GG_TREEVIS_SIMILARITY_SOURCE', '')
if (nzchar(similarity_source)) {
    implementation = new.env(parent=asNamespace('genegalleon.treevis'))
    sys.source(similarity_source, envir=implementation)
    add_sequence_similarity_column = implementation$add_sequence_similarity_column
}
sequence_add = add_sequence_similarity_column
g=list(tree=ggtree::ggtree(ape::read.tree(text='((a:1,b:1):1,c:2,d:2);')))
args=list(font_size=6)
tips=as.character(get_df_tip(g$tree)$label)
pairs=utils::combn(tips,2)
value=function(p,a,b) {
    k=which(apply(pairs,2,function(x)setequal(x,c(a,b))))
    p$data$similarity[4*(k-1)+1]
}
aln=tempfile(fileext='.fa.gz')
con=gzfile(aln,'wt'); writeLines(c('>a description','GctAAA---NNN','>b','GCCGAAACGNNN','>c','NNN---------'),con);close(con)
p=add_sequence_similarity_column(g,args,aln)$sequence_similarity
stopifnot(value(p,'a','b')==.5, is.na(value(p,'a','c')),is.na(value(p,'a','d')),
          attr(p,'treevis_width_mm')==15,grepl('amino acid',p$labels$x))
con=gzfile(aln,'wt'); writeLines(c('>a','ACDX-','>b','ACEAA','>c','XXXXX'),con);close(con)
p=add_sequence_similarity_column(g,args,aln,'protein')$sequence_similarity
stopifnot(abs(value(p,'a','b')-2/3)<1e-10,is.na(value(p,'a','c')))
con=gzfile(aln,'wt');writeLines(c('>a','GCN','>b','GCT'),con);close(con)
p=add_sequence_similarity_column(g,args,aln)$sequence_similarity
stopifnot(is.na(value(p,'a','b')))
con=gzfile(aln,'wt');writeLines(c('>a','TGA','>b','TGG'),con);close(con)
p=add_sequence_similarity_column(g,args,aln,genetic_code=1)$sequence_similarity
stopifnot(is.na(value(p,'a','b')))
p=add_sequence_similarity_column(g,args,aln,genetic_code=2)$sequence_similarity
stopifnot(value(p,'a','b')==1)
con=gzfile(aln,'wt');writeLines(c('>a','ACGT','>b','ACGT'),con);close(con)
stopifnot(inherits(try(add_sequence_similarity_column(g,args,aln),silent=TRUE),'try-error'))
con=gzfile(aln,'wt');writeLines(c('>a','ACG','>b','ACGT'),con);close(con)
stopifnot(inherits(try(add_sequence_similarity_column(g,args,aln),silent=TRUE),'try-error'))
con=gzfile(aln,'wt');writeLines(c('>a first','ACG','>a duplicate','ACG'),con);close(con)
stopifnot(inherits(try(add_sequence_similarity_column(g,args,aln),silent=TRUE),'try-error'))
stopifnot(identical(add_sequence_similarity_column(g,args,paste0(aln,'.missing')),g))

prom=tempfile(fileext='.fa.gz');con=gzfile(prom,'wt')
writeLines(c('>a description','ACGT','>b','ACGT','>c','ACGT'),con);close(con)
fimo=tempfile(fileext='.tsv')
hits=data.frame(motif_id=c('M1','M1','M2','M1','M3','M4'),
    sequence_name=c('a','a','a','b','b','c'),qvalue=c(.001,.002,.01,.001,.001,.1))
names(hits)[3]='q-value'
write.table(hits,fimo,sep='\t',quote=FALSE,row.names=FALSE)
p=add_cis_similarity_column(g,args,fimo,prom)$cis_similarity
stopifnot(abs(value(p,'a','b')-1/3)<1e-10,value(p,'a','c')==0,is.na(value(p,'a','d')))
p2=add_cis_similarity_column(g,args,fimo,prom,qvalue=.001)$cis_similarity
stopifnot(value(p2,'a','b')==.5)
write.table(hits[FALSE,],fimo,sep='\t',quote=FALSE,row.names=FALSE)
p=add_cis_similarity_column(g,args,fimo,prom)$cis_similarity
stopifnot(all(is.na(p$data$similarity)))
# Legacy commented header and trailing comments; never use p-values as q-values.
writeLines(c('#pattern name\tsequence name\tq-value','M1\ta\t0.001','M1\tb\t0.001','# command'),fimo)
p=add_cis_similarity_column(g,args,fimo,prom)$cis_similarity
stopifnot(value(p,'a','b')==1)
writeLines(c('motif_id\tsequence_name\tp-value','M1\ta\t0.001'),fimo)
stopifnot(inherits(try(add_cis_similarity_column(g,args,fimo,prom),silent=TRUE),'try-error'))
stopifnot(identical(add_cis_similarity_column(g,args,fimo,paste0(prom,'.missing')),g))
unlink(c(aln,prom,fimo))
cat('Sequence identity and promoter cis similarity tests passed.\n')

# The first available triangle owns one shared suffix column; later triangles
# never repeat it, and suffixes follow tree y coordinates, not FASTA order.
sg=list(tree=ggtree::ggtree(ape::read.tree(text='((gene320:1,geneABC:1):1,gene0.1:2);')))
fa=tempfile(fileext='.fa')
writeLines(c('>gene0.1','GCTAAA','>geneABC','GCCGAA','>gene320','GCCAAA'),fa)
sg=add_synteny_similarity_column(sg,args,paste0(fa,'.missing'))
stopifnot(!('pairwise_tip_suffix' %in% names(sg)))
sg=add_sequence_similarity_column(sg,args,fa)
stopifnot(identical(names(sg),c('tree','pairwise_tip_suffix','sequence_similarity')),
          attr(sg$pairwise_tip_suffix,'treevis_width_mm')==8,
          setequal(sg$pairwise_tip_suffix$data$suffix,c('...320','...ABC','...0.1')))
st=get_df_tip(sg$tree)
stopifnot(identical(sg$pairwise_tip_suffix$data$y, st$y))
sy=tempfile(fileext='.tsv')
write.table(data.frame(node_name=c('gene320','geneABC'),offset=1,group_id='G'),sy,
            sep='\t',quote=FALSE,row.names=FALSE)
sg=add_synteny_similarity_column(sg,args,sy)
stopifnot(sum(names(sg)=='pairwise_tip_suffix')==1,
          match('pairwise_tip_suffix',names(sg))+1==match('sequence_similarity',names(sg)))
unlink(c(fa,sy))
cat('Shared triangle suffix labels and missing-first-panel tests passed.\n')

# Identifier coercion must never merge biologically distinct motifs.
pr=tempfile(); f=tempfile()
writeLines(c('>a','ACGT','>b','ACGT','>c','ACGT'),pr)
writeLines(c('motif_id\tsequence_name\tq-value','001\ta\t0.001','1\tb\t0.001','NA\tc\t0.001'),f)
p=add_cis_similarity_column(g,args,f,pr)$cis_similarity
stopifnot(value(p,'a','b')==0,value(p,'a','c')==0)
for (bad in c('NA','NaN','Inf','-0.1','1.1','invalid','')) {
    writeLines(c('motif_id\tsequence_name\tq-value',paste0('M1\ta\t',bad),'M1\tb\t0.001'),f)
    stopifnot(inherits(try(add_cis_similarity_column(g,args,f,pr),silent=TRUE),'try-error'))
}
writeLines(c('not a fasta alignment'),f)
stopifnot(inherits(try(add_sequence_similarity_column(g,args,f),silent=TRUE),'try-error'))
unlink(c(f,pr))
cat('Exact motif IDs and invalid q-value/FASTA checks passed.\n')

# A position-wise oracle classifies each residue independently. Compare every
# panel row and built drawing layer, including the zero/NA distinction, rather
# than checking only a few displayed pairs.
pairwise_panel = getFromNamespace('treevis_pairwise_column','genegalleon.treevis')
read_similarity_fasta = getFromNamespace('treevis_similarity_fasta','genegalleon.treevis')
oracle_sequence_values = function(tips, sequences) {
    pairs = utils::combn(as.character(tips$label),2)
    vapply(seq_len(ncol(pairs)),function(k) {
        a = sequences[[pairs[1,k]]]; b = sequences[[pairs[2,k]]]
        if (is.null(a) || is.null(b)) return(NA_real_)
        a = strsplit(a,'',fixed=TRUE)[[1]]
        b = strsplit(b,'',fixed=TRUE)[[1]]
        valid = which(grepl('^[ACDEFGHIKLMNPQRSTVWY]$',a) &
                      grepl('^[ACDEFGHIKLMNPQRSTVWY]$',b))
        if (!length(valid)) return(NA_real_)
        mean(a[valid] == b[valid])
    },numeric(1))
}
render_similarity = function(plot) {
    path = tempfile(fileext='.png')
    grDevices::png(path,width=700,height=500)
    grob = ggplotGrob(plot)
    bbox_mm = c(grid::convertWidth(sum(grob$widths),'mm',valueOnly=TRUE),
                grid::convertHeight(sum(grob$heights),'mm',valueOnly=TRUE))
    print(plot)
    grDevices::dev.off()
    bytes = readBin(path,'raw',n=file.info(path)$size)
    unlink(path)
    list(bytes=bytes,bbox_mm=bbox_mm)
}
check_sequence_panel = function(sequences, render=FALSE) {
    path = tempfile(fileext='.fa')
    writeLines(as.vector(rbind(paste0('>',names(sequences),' descriptions'),unlist(sequences))),path)
    actual = sequence_add(g,args,path,'protein')
    tips = get_df_tip(g$tree); tips = tips[order(tips$y),,drop=FALSE]
    expected = pairwise_panel(g,args,tips,
        oracle_sequence_values(tips,as.list(read_similarity_fasta(path))),
        'sequence_similarity','Sequence\nidentity\n(amino acid)',15)
    stopifnot(identical(names(actual),names(expected)))
    for (name in c('pairwise_tip_suffix','sequence_similarity')) {
        stopifnot(identical(actual[[name]]$data,expected[[name]]$data),
            identical(actual[[name]]$labels,expected[[name]]$labels),
            identical(attr(actual[[name]],'treevis_width_mm'),attr(expected[[name]],'treevis_width_mm')),
            identical(ggplot_build(actual[[name]])$data,ggplot_build(expected[[name]])$data))
        if (render) stopifnot(identical(render_similarity(actual[[name]]),render_similarity(expected[[name]])))
    }
    unlink(path)
}
check_sequence_panel(list(c='XXXXXXXXXXXXXXXXXXXX',b='ACDEFGHIKLMNPQRSTVWY',
                          a='ACDEFGHIKLMNPQRSTVWY'),render=TRUE)
check_sequence_panel(list(a='ACDEFGHIKLMNPQRSTVWYXBZJOU*-?é',
                          c=strrep('-',30),
                          b='YWVTSRQPNMLKIHGFEDCAXXXXXXX??X'),render=TRUE)
check_sequence_panel(list(c='acdx-',a='ACDX-',b='ACEAA'))
set.seed(3714)
for (trial in seq_len(12)) {
    alphabet = strsplit('ACDEFGHIKLMNPQRSTVWYXBZJOU*-?','',fixed=TRUE)[[1]]
    sequences = setNames(lapply(seq_len(4),function(i) paste0(sample(alphabet,83,replace=TRUE),collapse='')),
                         sample(c('a','b','c','d')))
    check_sequence_panel(sequences)
}
cat('Complete sequence panel data, drawing layers and raster equivalence passed.\n')

# Check exact failures as well as successful values. These inputs fail before
# translation or drawing, so no external tool message is part of the contract.
check_sequence_error = function(lines,message,mode='protein') {
    path = tempfile(fileext='.fa')
    writeLines(lines,path)
    failure = tryCatch(sequence_add(g,args,path,mode),error=function(error) error)
    unlink(path)
    stopifnot(inherits(failure,'simpleError'),identical(conditionMessage(failure),message))
}
check_sequence_error(c('>a','ACD','>b','ACD'),'Sequence similarity mode must be cds or protein.','unknown')
check_sequence_error('ACDEFG','Nonempty FASTA input has no record headers.')
check_sequence_error(c('ACD','>a','ACD','>b','ACD'),'FASTA sequence appeared before its header.')
check_sequence_error(c('>a first','ACD','>a second','ACD'),'FASTA IDs must be nonempty and unique.')
check_sequence_error(c('>a','','>b','ACD'),
    'Sequence similarity requires a nonempty, equal-length alignment; raw unaligned FASTA is not valid.')
check_sequence_error(c('>a','ACD','>b','ACDE'),
    'Sequence similarity requires a nonempty, equal-length alignment; raw unaligned FASTA is not valid.')
check_sequence_error(c('>a','ACGT','>b','ACGT'),
    'Amino-acid identity requires a codon alignment with length divisible by three.','cds')
cat('Exact sequence similarity validation failures passed.\n')

# Read the installed FIMO's real output, including its header and footer.
run_dir=tempfile('cis-fimo-integration-');dir.create(run_dir)
motif=file.path(run_dir,'motif.meme');prom=file.path(run_dir,'prom.fa')
writeLines(c('MEME version 4','','ALPHABET= ACGT','','strands: + -','',
 'Background letter frequencies','A 0.25 C 0.25 G 0.25 T 0.25','',
 'MOTIF M1','letter-probability matrix: alength= 4 w= 8 nsites= 20 E= 0',
 rep(c('1 0 0 0','0 1 0 0','0 0 1 0','0 0 0 1'),2)),motif)
writeLines(c('>a','ACGTACGTACGTACGTACGTACGT','>b','ACGTACGTACGTACGTACGTACGT',
             '>c','AAAAAAAAAAAAAAAAAAAAAAAA'),prom)
status=system2('fimo',c('--oc',shQuote(file.path(run_dir,'fimo')),shQuote(motif),shQuote(prom)),
               stdout=file.path(run_dir,'stdout'),stderr=file.path(run_dir,'stderr'))
stopifnot(status==0)
p=add_cis_similarity_column(g,args,file.path(run_dir,'fimo','fimo.tsv'),prom)$cis_similarity
stopifnot(!is.null(p),value(p,'a','b')==1,value(p,'a','c')==0)
unlink(run_dir,recursive=TRUE)
cat('Real FIMO output integration passed.\n')
