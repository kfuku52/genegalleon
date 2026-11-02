# Fixed-universe event sensitivity: analytical boundaries and real public CLI.
e <- new.env()
for (x in parse('workflow/support/cafe_go_enrichment.r')) {
  if (is.call(x) && as.character(x[[1]]) %in% c('<-', '=') &&
      is.call(x[[3]]) && identical(x[[3]][[1]], as.name('function'))) eval(x, e)
}
expect_error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error=identity), 'error'))

# Equality is untestable under the documented strict P < alpha convention.
t <- e$tarone_testability(c(.001, .025), .05)
stopifnot(t$factor==2, t$n_testable==1, identical(t$testable,c(TRUE,FALSE)), t$threshold==.025)
none <- e$tarone_testability(c(1,1),.05)
all <- e$tarone_testability(rep(0,3),.05)
stopifnot(none$factor==1, none$n_testable==0, none$union_bound==0,
          all$factor==3, all$n_testable==3, all$union_bound<=.05)
for (bad in list(NA_real_, NaN, -.01, 1.01, numeric())) expect_error(e$tarone_testability(bad))
expect_error(e$tarone_testability(.01,0))
expect_error(e$tarone_testability(.01,1))

# Enumerate every null label assignment. Fixed overlapping GO sets need not be
# independent; the Tarone union bound still cannot exceed alpha for these valid
# conditional tests. Tarone contains the ordinary Bonferroni rejection set.
N <- 8L; n <- 2L
sets <- list(c(1,2), c(1,2,3), c(4,5), 4L, c(6,7,8))
K <- lengths(sets)
min_p <- dhyper(pmin(n,K),K,N-K,n)
for (alpha in c(.05,.1,.2,.4)) {
  screen <- e$tarone_testability(min_p,alpha)
  labels <- combn(N,n,simplify=FALSE)
  rejected <- vapply(labels,function(label) {
    hits <- vapply(sets,function(s) sum(label %in% s),integer(1))
    p <- phyper(hits-1L,K,N-K,n,lower.tail=FALSE)
    tarone <- screen$testable & p<screen$threshold
    bonferroni <- p<alpha/length(sets)
    stopifnot(all(!bonferroni | tarone))
    any(tarone)
  },logical(1))
  stopifnot(mean(rejected)<=alpha+1e-12, screen$union_bound<=alpha+1e-12)
}

# The annotation background contains a GO with no events, a GO with only other
# branch events, and a universally present GO. All remain in every branch.
population <- data.frame(FamilyID=paste0('F',1:100),branch=c(rep('A',10),rep('B',90)),
                         is_target=c(rep(TRUE,10),rep(FALSE,90)),stringsAsFactors=FALSE)
linked <- population[c(1:4,11:16),];linked$go_ids<-'GO:T';linked$go_aspects<-'BP';linked$go_terms<-'target'
common <- population;common$go_ids<-'GO:C';common$go_aspects<-'CC';common$go_terms<-'common'
other <- population[31,];other$go_ids<-'GO:O';other$go_aspects<-'MF';other$go_terms<-'other'
events <- rbind(linked,common,other)
background <- events[,c('FamilyID','go_ids','go_aspects','go_terms')]
background <- rbind(background,data.frame(FamilyID='no_event',go_ids='GO:ZERO',go_aspects='BP',go_terms='no events'))
universe <- e$event_go_universe(background,c(population$FamilyID,'no_event'))
legacy <- e$summarise_go_enrichment(events,10,90,.05)
fixed <- e$summarise_event_all_go(events,10,90,universe,.05)
stopifnot(nrow(fixed$all)==4,nrow(legacy$all)==2,
          all(fixed$all$p_value[fixed$all$go_ids %in% c('GO:O','GO:ZERO')]==1),
          fixed$all$minimum_attainable_p[fixed$all$go_ids=='GO:ZERO']==1,
          fixed$all$n_change_out_target_in_go[fixed$all$go_ids=='GO:O']==1,
          fixed$all$minimum_attainable_p[fixed$all$go_ids=='GO:C']==1,
          identical(fixed$all$p_value_adjusted,p.adjust(fixed$all$p_value,'BH')))
idx <- match(legacy$all$go_ids,fixed$all$go_ids)
stopifnot(identical(legacy$all$p_value,fixed$all$p_value[idx]),
          all(fixed$all$minimum_attainable_p<=fixed$all$p_value+1e-12))
for (i in seq_len(nrow(fixed$all))) {
  row <- fixed$all[i,]
  support <- row$n_change_in_target_in_go+row$n_change_out_target_in_go
  attainable <- phyper(seq.int(max(0,10-(100-support)),min(10,support))-1,
                      support,100-support,10,lower.tail=FALSE)
  stopifnot(abs(min(attainable)-row$minimum_attainable_p)<1e-12)
}
zero <- e$summarise_event_all_go(events[0,],0,0,universe,.05)
stopifnot(nrow(zero$all)==4,all(zero$all$p_value==1),all(zero$all$minimum_attainable_p==1),
          zero$screen$n_testable==0,nrow(zero$significant)==0,nrow(zero$tarone_significant)==0)
expect_error(e$event_go_universe(rbind(background,transform(background[1,],go_terms='conflict')),population$FamilyID))
expect_error(e$summarise_event_all_go(events,1,1,universe))

root <- normalizePath('.')
tmp <- tempfile('go-tarone-');dir.create(tmp)
tryCatch({
  put <- function(d,name)e$write_tsv_base(d,file.path(tmp,name))
  families <- c(population$FamilyID,'no_event','not_in_cafe')
  changes <- data.frame(FamilyID=families[-length(families)],A=c(rep(1,10),rep(0,91)),B=c(rep(0,10),rep(1,90),0))
  put(changes,'change.tsv');put(transform(changes,A=.001,B=.001),'prob.tsv')
  put(data.frame(Orthogroup=families,ref=paste0('g',seq_along(families))),'ids.tsv')
  annotation <- rbind(background,data.frame(FamilyID='not_in_cafe',go_ids='GO:EXCLUDED',go_aspects='BP',go_terms='excluded'))
  annotation$gene_id <- paste0('g',match(annotation$FamilyID,families))
  put(annotation[,c('gene_id','go_ids','go_aspects','go_terms')],'ref_annotation.tsv')
  cli <- function(argv,success=TRUE) {
    out <- suppressWarnings(system2(file.path(R.home('bin'),'Rscript'),
      shQuote(c(file.path(root,'workflow/support/cafe_go_enrichment.r'),argv)),stdout=TRUE,stderr=TRUE))
    status <- attr(out,'status');if(is.null(status))status<-0
    stopifnot(if(success)status==0 else status!=0)
    out
  }
  freeze <- c('--freeze-event-go-universe',file.path(tmp,c('change.tsv','ids.tsv','ref_annotation.tsv','universe.tsv')),'BP,MF,CC','ref')
  cli(freeze)
  frozen <- e$read_tsv_base(file.path(tmp,'universe.tsv'))
  stopifnot(nrow(frozen)==4,!('GO:EXCLUDED' %in% frozen$go_ids))
  cli(freeze,success=FALSE) # Freeze cannot overwrite an existing catalog.
  put(transform(changes,A=0,B=0),'different_changes.tsv')
  independent <- freeze;independent[2]<-file.path(tmp,'different_changes.tsv');independent[5]<-file.path(tmp,'universe2.tsv')
  cli(independent)
  stopifnot(identical(readLines(file.path(tmp,'universe.tsv')),readLines(file.path(tmp,'universe2.tsv'))))
  args <- c(file.path(tmp,c('change.tsv','prob.tsv','ids.tsv','ref_annotation.tsv')),
            file.path(tmp,'branch_A'),'A','increase','BP,MF,CC','event_all_go','.05','ref',file.path(tmp,'universe.tsv'))
  cli(args)
  actual <- e$read_tsv_base(file.path(tmp,'branch_A/enrichment_significant_increase_A_all_go.tsv'))
  stopifnot(setequal(actual$go_ids,frozen$go_ids),max(abs(actual$p_value_adjusted-p.adjust(actual$p_value,'BH')))<1e-12)
  bargs<-args;bargs[5]<-file.path(tmp,'branch_B');bargs[6]<-'B';cli(bargs)
  b <- e$read_tsv_base(file.path(tmp,'branch_B/enrichment_significant_increase_B_all_go.tsv'))
  stopifnot(setequal(actual$go_ids,b$go_ids),actual$n_change_in_target_in_go[actual$go_ids=='GO:O']==0,
            b$n_change_in_target_in_go[b$go_ids=='GO:O']==1)
  meta <- e$read_tsv_base(file.path(tmp,'branch_A/event_all_go_metadata.tsv'))
  stopifnot(meta$n_background_go==4,meta$n_tarone_testable<=meta$tarone_factor,
            meta$tarone_union_bound<=.05+1e-12)
  frozen$go_terms[1]<-'tampered';put(frozen,'bad_universe.tsv')
  bad<-args;bad[12]<-file.path(tmp,'bad_universe.tsv');cli(bad,success=FALSE)
  stopifnot(!file.exists(file.path(tmp,'branch_A/enrichment_significant_increase_A_all_go.tsv')),
            file.exists(file.path(tmp,'universe.tsv')))
},finally=unlink(tmp,recursive=TRUE))
cat('Fixed event GO universe, Tarone testability and CLI regressions passed\n')
