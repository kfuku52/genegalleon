suppressPackageStartupMessages(library(genegalleon.treevis))
suppressPackageStartupMessages(library(ggplot2))
pdf(NULL)

# This independent point-membership oracle deliberately does not use endpoint
# searching or a sweep. Keep the public plot constructor and its private helper
# independently testable without exporting a new API.
pointwise_keys <- function(is_atgc, mapped_untrim_nt, df_seq_rps,
                           domain_rank, non_domain_fill, key_sep) {
  active <- vector('list', length(is_atgc))
  if (nrow(df_seq_rps)) {
    qlen <- suppressWarnings(as.integer(df_seq_rps[1, 'qlen']))
    if (!is.na(qlen) && qlen > 0) {
      ord <- order(ifelse(df_seq_rps$sacc %in% names(domain_rank),
                          domain_rank[df_seq_rps$sacc], Inf),
                   suppressWarnings(as.numeric(df_seq_rps$qstart)),
                   method = 'radix', na.last = TRUE)
      for (k in ord) {
        first <- suppressWarnings(as.integer(df_seq_rps[k, 'qstart']))
        last <- suppressWarnings(as.integer(df_seq_rps[k, 'qend']))
        if (!is.finite(first) || !is.finite(last)) next
        first <- max(1L, (first - 1L) * 3L + 1L)
        last <- min(qlen * 3L, last * 3L)
        if (first > last) next
        hit <- which(is_atgc & !is.na(mapped_untrim_nt) &
                       mapped_untrim_nt >= first & mapped_untrim_nt <= last)
        label <- as.character(df_seq_rps[k, 'sacc'])
        for (h in hit) {
          if (!(label %in% active[[h]])) active[[h]] <- c(active[[h]], label)
        }
      }
    }
  }
  out <- rep(NA_character_, length(is_atgc))
  out[which(is_atgc)] <- non_domain_fill
  for (h in which(lengths(active) > 0L)) {
    labels <- unique(as.character(active[[h]]))
    ord <- order(ifelse(labels %in% names(domain_rank), domain_rank[labels], Inf),
                 labels, method = 'radix', na.last = TRUE)
    out[h] <- paste(labels[ord], collapse = key_sep)
  }
  out
}

source_path <- Sys.getenv('GG_TREEVIS_ALIGNMENT_SOURCE', '')
if (nzchar(source_path)) {
  implementation <- new.env(parent = asNamespace('genegalleon.treevis'))
  sys.source(source_path, implementation)
  alignment_add <- implementation$add_alignment_column
  alignment_keys <- implementation$.alignment_domain_keys
} else {
  alignment_add <- add_alignment_column
  alignment_keys <- getFromNamespace('.alignment_domain_keys', 'genegalleon.treevis')
}
oracle_env <- new.env(parent = environment(alignment_add))
oracle_env$.alignment_domain_keys <- pointwise_keys
oracle_add <- alignment_add
environment(oracle_add) <- oracle_env
separator <- ':::DOMAINSEP:::'
ranks <- c(Z = 1L, A = 2L, B = 3L)
hits <- function(labels, start, end, qlen = 8) {
  data.frame(sacc = labels, qstart = start, qend = end, qlen = qlen,
             stringsAsFactors = FALSE)
}
check_keys <- function(mask, mapping, domains, rank = ranks) {
  expected <- pointwise_keys(mask, mapping, domains, rank, '__non_domain__', separator)
  actual <- alignment_keys(mask, mapping, domains, rank, '__non_domain__', separator)
  stopifnot(identical(actual, expected))
  invisible(actual)
}

# Inclusive codon endpoints, gaps, overlapping copies of one domain, unknown
# labels, reversed/clipped hits, fractional coordinates and first-row qlen.
mask <- c(TRUE, TRUE, FALSE, rep(TRUE, 19))
mapping <- rep(NA_integer_, length(mask))
mapping[which(mask)] <- c(1:6, 10:24)
edge <- hits(c('A', 'Z', 'A', 'B', 'unknown', NA, '', paste0('A', separator, 'B')),
             c(1, 2, 3, -2, 5.8, 4, 7, 6), c(3, 4, 6, 1, 9.9, 5, 8, 7),
             c(8, rep(100, 7)))
check_keys(mask, mapping, edge)
check_keys(mask, rep(NA_integer_, length(mask)), edge)
check_keys(rep(FALSE, 8), rep(NA_integer_, 8), edge)
check_keys(logical(0), integer(0), edge)
check_keys(mask, mapping, hits(c('A', 'A', 'A'), c(1, 2, 7), c(3, 5, 8)))
check_keys(mask, mapping, hits(c('A', 'Z', 'B', 'B'),
                             c(NA, Inf, 10, 4), c(1, 2, 2, -1)))
check_keys(mask, mapping, hits('A', 1, 8, qlen = NA_real_))
check_keys(mask, mapping, hits('A', 1, 8, qlen = 0))
check_keys(mask, mapping, hits(character(0), numeric(0), numeric(0), numeric(0)))
stopifnot(identical(check_keys(rep(TRUE, 9), 1:9, hits(c('A', 'B'), c(1, 2), c(2, 3))),
                    c(rep('A', 3), rep(paste0('A', separator, 'B'), 3), rep('B', 3))))
for (fn in list(pointwise_keys, alignment_keys)) {
  overflow <- suppressWarnings(try(fn(TRUE, 1L, hits('A', 1000000000, 1000000001),
                                     ranks, '__non_domain__', separator), silent = TRUE))
  stopifnot(inherits(overflow, 'try-error'))
}

# Randomized point membership is especially useful for coincident starts/ends,
# duplicate labels and deleted untrimmed bases. Fixed seed keeps CI bounded.
set.seed(1584)
for (i in seq_len(180)) {
  size <- sample(1:80, 1)
  mask <- sample(c(TRUE, FALSE), size, replace = TRUE)
  mapping <- rep(NA_integer_, size)
  mapping[which(mask)] <- cumsum(sample(1:3, sum(mask), replace = TRUE))
  if (i %% 13 == 0) mapping[] <- NA_integer_
  nhit <- sample(0:24, 1)
  domains <- hits(sample(c('A', 'Z', 'B', 'unknown', '', NA_character_), nhit, replace = TRUE),
                  sample(-3:35, nhit, replace = TRUE), sample(-3:35, nhit, replace = TRUE),
                  rep(sample(1:30, 1), nhit))
  check_keys(mask, mapping, domains)
}

g <- list(tree = ggtree::ggtree(ape::read.tree(text = '((a:1,b:1):1,c:2,d:2);')))
args <- list(font_size = 6, margins = rep(0, 4))
seqs <- list(a = c(rep('aa', 3), '04', rep('aa', 6)), b = rep('04', 10),
             c = c(rep('cc', 9), '04'))
untrim <- list(a = rep('aa', 9), b = rep('04', 12), c = rep('aa', 9))
rps <- data.frame(qacc = c('a', 'a', 'a', 'a', 'a', 'c', 'absent'),
                  sacc = paste0('PF', 1:7),
                  stitle = c('entry,A', 'entry,B', 'entry,A', 'entry,Z', 'entry,',
                             'entry,C', 'entry,unused'),
                  qstart = c(1, 2, 1, NA, 5, 1, 1), qend = c(2, 3, 1, 3, 5, 3, 3),
                  qlen = 3, slen = 3, stringsAsFactors = FALSE)
check_plot <- function(sequences, untrimmed, domain_hits) {
  actual <- alignment_add(g, args, sequences, domain_hits, untrimmed)
  expected <- oracle_add(g, args, sequences, domain_hits, untrimmed)
  stopifnot(identical(names(actual), names(expected)))
  if (!is.null(actual$alignment)) {
    stopifnot(identical(actual$alignment$labels, expected$alignment$labels),
              identical(lapply(actual$alignment$layers, function(layer) layer$data),
                        lapply(expected$alignment$layers, function(layer) layer$data)),
              identical(ggplot_build(actual$alignment)$data,
                        ggplot_build(expected$alignment)$data))
    rect <- actual$alignment$layers[[2]]$data
    stopifnot(is.integer(rect$split_index), is.integer(rect$split_total),
              is.double(rect$xmin), is.double(rect$xmax))
  }
  invisible(actual)
}
check_plot(seqs, untrim, rps)
check_plot(seqs, NULL, rps)
check_plot(seqs, untrim, NULL)
stopifnot(identical(check_plot(NULL, untrim, rps), g))
check_plot(list(b = rep('04', 10)), untrim, rps)
check_plot(list(a = character(0), b = rep('04', 10)), untrim, rps)

# A hand-calculated rectangle fixture also checks the zero-based endpoints and
# vertical split ordering, beyond agreement with the independent key oracle.
p <- check_plot(seqs, untrim, rps[1:2, ])$alignment
rect <- p$layers[[2]]$data
a <- rect[as.character(rect$label) == 'a', ]
stopifnot(identical(a$xmin, c(0, 4, 4, 7)), identical(a$xmax, c(2, 6, 6, 9)),
          identical(as.character(a$fill), c('A', 'A', 'B', 'B')),
          identical(a$split_index, c(1L, 1L, 2L, 1L)),
          identical(a$split_total, c(1L, 2L, 2L, 1L)),
          identical(a$ymax[2], a$ymin[3]),
          identical(levels(rect$label), as.character(get_df_tip(g$tree)$label)))
cat('Alignment interval, point-membership and complete ggplot layer equivalence tests passed.\n')
