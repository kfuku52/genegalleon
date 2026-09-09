suppressPackageStartupMessages(library(genegalleon.treevis))
suppressPackageStartupMessages(library(ggplot2))
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
