# Run from the repository root in a GeneGalleon container.
e <- new.env()
for (x in parse('workflow/support/cafe_go_enrichment.r')) {
  if (is.call(x) && as.character(x[[1]]) %in% c('<-', '=') &&
      is.call(x[[3]]) && identical(x[[3]][[1]], as.name('function'))) eval(x, e)
}
expect_error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), 'error'))

# Existing family-wide P values are adjusted before branch/sign screening.
ids <- c('gain','loss','broad','mixed','flat','background','unreported')
changes <- data.frame(FamilyID=ids,check.names=FALSE)
changes[['A<1>']] <- c(3,-3,3,-3,0,1,1)
changes[['B<2>']] <- c(0,0,4,4,0,0,0)
changes[['C<3>']] <- 0L; changes[['<4>']] <- 0L; changes[['<5>']] <- 0L
prob <- changes[1:6,]; prob[-1] <- .9
prob[['A<1>']] <- .01; prob[['B<2>']][3:4] <- .01; prob[['<4>']] <- NA_real_
fam <- data.frame(FamilyID=ids,pvalue=c(.001,.002,.001,.001,.001,.9,.8))
native <- e$read_cafe_branch_families(changes,prob,fam,ids,'A<1>','<4>')
inc <- e$select_cafe_branch_families(native,.05,'increase')
dec <- e$select_cafe_branch_families(native,.05,'decrease')
both <- e$select_cafe_branch_families(native,.05,'both')
stopifnot(identical(inc$FamilyID[inc$selected],'gain'),identical(dec$FamilyID[dec$selected],'loss'),
  setequal(both$FamilyID[both$selected],c('gain','loss')),
  identical(both$family_p_value_adjusted,p.adjust(fam$pvalue,'BH')),
  identical(inc$family_p_value_adjusted,dec$family_p_value_adjusted),
  both$selection_status[both$FamilyID=='mixed']=='other_branches_flagged',
  is.na(both$n_other_flagged[both$FamilyID=='unreported']),!both$selected[7])
bad <- native; bad$family_p_value[7] <- .0001
expect_error(e$select_cafe_branch_families(bad,.05,'both'))
bad <- prob; bad[['B<2>']][1] <- NA_real_
expect_error(e$read_cafe_branch_families(changes,bad,fam,ids,'A<1>','<4>'))
bad <- prob; bad[['<4>']] <- 0
expect_error(e$read_cafe_branch_families(changes,bad,fam,ids,'A<1>','<4>'))
bad <- fam; bad$pvalue[1] <- NaN
expect_error(e$read_cafe_branch_families(changes,prob,bad,ids,'A<1>','<4>'))
expect_error(e$read_cafe_branch_families(changes,rbind(prob,prob[1,]),fam,ids,'A<1>','<4>'))
expect_error(e$read_cafe_branch_families(changes,prob,fam[-1,],ids,'A<1>','<4>'))
expect_error(e$read_cafe_branch_families(changes,prob,fam,ids,'<4>','<4>'))
none_reported <- e$read_cafe_branch_families(changes,prob[0,],transform(fam,pvalue=1),ids,'A<1>','<4>')
stopifnot(!any(e$select_cafe_branch_families(none_reported,.05,'both')$selected))
asr <- tempfile(); writeLines(c('#nexus','BEGIN TREES;',
  ' TREE family = ((A<1>_9:1,B<2>_3:1)<5>_4:1,C<3>_4:2)<4>_4;','END;'),asr)
stopifnot(e$read_cafe_root_branch(asr)=='<4>'); unlink(asr)
# One row per family; rejected candidate GO terms retain P=1.
ann <- data.frame(FamilyID=ids,go_ids=paste0('GO:',seq_along(ids)),go_aspects='BP',go_terms=ids)
candidate <- unique(ann[,c('go_ids','go_aspects','go_terms')])
g <- e$summarise_family_go(inc,rbind(ann,ann),candidate)
stopifnot(setequal(g$all$go_ids,candidate$go_ids),g$all$n_selected_in_go[g$all$go_ids=='GO:1']==1,
  all(g$all$p_value[g$all$go_ids!='GO:1']==1))
none <- inc; none$selected <- FALSE
g0 <- e$summarise_family_go(none,ann,candidate)
stopifnot(nrow(g0$all)==7,all(g0$all$p_value==1),nrow(g0$significant)==0,
  nrow(e$summarise_family_go(inc,ann,candidate[0,])$all)==0)

# Legacy default remains unchanged, including its original GO filtering rule.
events <- data.frame(FamilyID = paste0('F', 1:100), is_target = c(rep(TRUE, 10), rep(FALSE, 90)))
rows <- events[c(1:4, 11:16), ]; rows$go_ids <- 'GO:T'; rows$go_aspects <- 'BP'; rows$go_terms <- 'target'
common <- events; common$go_ids <- 'GO:C'; common$go_aspects <- 'BP'; common$go_terms <- 'common'
rows <- rbind(rows, common)
for (i in 1:8) {
  z <- events[17+i, ]; z$go_ids <- paste0('GO:O', i); z$go_aspects <- 'BP'; z$go_terms <- paste0('other', i)
  rows <- rbind(rows, z)
}
legacy <- e$summarise_go_enrichment(rows, 10, 90, .05)
stopifnot(nrow(legacy$all) == 2, abs(legacy$all$p_value_adjusted[legacy$all$go_ids == 'GO:T'] - .01644975288516) < 1e-12)

cat('CAFE GO native-output selection and legacy numerical regressions passed\n')

# Eight-argument legacy CLI and explicit event mode must remain equivalent,
# without any of the new native-comparison input files being present.
root <- normalizePath('.')
tmp <- tempfile('cafe-go-default-'); dir.create(tmp)
tryCatch({
  put <- function(d, name) e$write_tsv_base(d, file.path(tmp, name))
  put(data.frame(FamilyID = events$FamilyID, A = as.integer(events$is_target),
                 B = as.integer(!events$is_target)), 'change.tsv')
  put(data.frame(FamilyID = events$FamilyID, A = .001, B = .001), 'prob.tsv')
  put(data.frame(Orthogroup = events$FamilyID, ref = paste0('g', 1:100)), 'ids.tsv')
  annotation <- rows[, c('go_ids', 'go_aspects', 'go_terms')]
  annotation$gene_id <- paste0('g', match(rows$FamilyID, events$FamilyID))
  put(annotation, 'ref.annotation.tsv')
  run <- function(out, extra = character()) {
    argv <- c(file.path(root, 'workflow/support/cafe_go_enrichment.r'),
              file.path(tmp, c('change.tsv', 'prob.tsv', 'ids.tsv', 'ref.annotation.tsv')),
              file.path(tmp, out), 'A', 'increase', 'BP', extra)
    result <- suppressWarnings(system2(file.path(R.home('bin'), 'Rscript'), shQuote(argv), stdout = TRUE, stderr = TRUE))
    if (!is.null(attr(result, 'status'))) stop(paste(result, collapse = '\n'))
  }
  run('default'); run('explicit', 'event')
  for (name in list.files(file.path(tmp, 'default'))) {
    stopifnot(identical(readLines(file.path(tmp, 'default', name)), readLines(file.path(tmp, 'explicit', name))))
  }
  actual <- e$read_tsv_base(file.path(tmp, 'default', 'enrichment_significant_increase_A_all_go.tsv'))
  stopifnot(nrow(actual) == 2, abs(actual$p_value_adjusted[actual$go_ids == 'GO:T'] - .01644975288516) < 1e-12)
}, finally = unlink(tmp, recursive = TRUE))
cat('CAFE GO default CLI equivalence passed\n')

# A bad new input must not leave old optional GO summaries published.
tmp <- tempfile('cafe-go-failure-'); dir.create(tmp)
tryCatch({
  names <- c('family_branch_flags.tsv', 'branch_flags_metadata.tsv',
    'enrichment_significant_both_A<1>_all_go.tsv',
    'enrichment_significant_both_A<1>_significant_go.tsv')
  for (name in names) writeLines('old success', file.path(tmp, name))
  dir.create(file.path(tmp, 'source_evidence'))
  writeLines('retained native evidence', file.path(tmp, 'source_evidence', 'cache.txt'))
  argv <- c(file.path(root, 'workflow/support/cafe_go_enrichment.r'),
    rep(file.path(tmp, 'missing.tsv'), 4), tmp, 'A<1>', 'both', 'BP', 'cafe_branch_flags')
  result <- suppressWarnings(system2(file.path(R.home('bin'), 'Rscript'), shQuote(argv), stdout=TRUE, stderr=TRUE))
  stopifnot(!is.null(attr(result, 'status')), !any(file.exists(file.path(tmp, names))),
    file.exists(file.path(tmp, 'source_evidence', 'cache.txt')))
}, finally=unlink(tmp, recursive=TRUE))
cat('CAFE GO failed-input publication regression passed\n')

# Exercise the full CLI when CAFE reports no significant families, and when an
# internal target has both gains and losses.
tmp <- tempfile('cafe-go-boundary-'); dir.create(tmp)
tryCatch({
  put <- function(d, name) e$write_tsv_base(d, file.path(tmp, name))
  native_changes <- changes
  native_changes[['<4>']] <- 0L
  put(native_changes, 'Base_change.tab')
  writeLines(c('#nexus', 'BEGIN TREES;',
    ' TREE family = ((A<1>_9:1,B<2>_3:1)<5>_4:1,C<3>_4:2)<4>_4;', 'END;'), file.path(tmp, 'Base_asr.tre'))
  put(data.frame(Orthogroup=ids,ref=paste0('g',seq_along(ids))), 'ids.tsv')
  put(data.frame(gene_id=paste0('g',seq_along(ids)),go_ids=paste0('GO:',seq_along(ids)),
    go_aspects='BP',go_terms=ids), 'ref.annotation.tsv')
  run <- function(out, target='A<1>') {
    argv <- c(file.path(root,'workflow/support/cafe_go_enrichment.r'),
      file.path(tmp,c('Base_change.tab','Base_branch_probabilities.tab','ids.tsv','ref.annotation.tsv')),
      file.path(tmp,out),target,'both','BP','cafe_branch_flags','.05')
    result <- suppressWarnings(system2(file.path(R.home('bin'),'Rscript'),shQuote(argv),stdout=TRUE,stderr=TRUE))
    if (!is.null(attr(result,'status'))) stop(paste(result,collapse='\n'))
    file.path(tmp,out)
  }
  put(prob[0,], 'Base_branch_probabilities.tab')
  put(transform(fam,pvalue=1), 'Base_family_results.txt')
  out <- run('none')
  screened <- e$read_tsv_base(file.path(out,'family_branch_flags.tsv'))
  stopifnot(nrow(screened)==length(ids),!any(screened$selected),!any(screened$branch_reported),
    nrow(e$read_tsv_base(file.path(out,'enrichment_significant_both_A<1>_all_go.tsv')))==0)
  native_changes[['<5>']][1:2] <- c(3,-3)
  native_changes[['A<1>']][1:2] <- 0L
  put(native_changes,'Base_change.tab')
  native_prob <- prob; native_prob[['<5>']][1:2] <- .001
  # Write native root N/A, not the output-table NA spelling.
  write.table(native_prob,file.path(tmp,'Base_branch_probabilities.tab'),sep='\t',quote=FALSE,row.names=FALSE,na='N/A')
  put(fam,'Base_family_results.txt')
  out <- run('internal','<5>')
  screened <- e$read_tsv_base(file.path(out,'family_branch_flags.tsv'))
  stopifnot(setequal(screened$FamilyID[screened$selected],c('gain','loss')))
},finally=unlink(tmp,recursive=TRUE))
cat('CAFE GO empty-report and internal-target CLI regressions passed\n')
