#!/usr/bin/env Rscript
# Compare the real sequence-similarity calculation and panel construction in
# alternating fresh R processes. Input generation and output hashing are untimed.
arguments = commandArgs(TRUE)
option = function(name, default=NULL) {
    where = match(name, arguments)
    if (is.na(where)) return(default)
    if (where == length(arguments)) stop(paste('Missing value for', name))
    arguments[[where+1L]]
}
flag = function(name) name %in% arguments
script = normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[[1]]))
profiles = list(
    small=list(tips=7L, sites=300L, ambiguous=FALSE, kernel_calls=30L, panel_calls=3L),
    dense64=list(tips=64L, sites=300L, ambiguous=FALSE, kernel_calls=3L, panel_calls=1L),
    dense256=list(tips=256L, sites=1000L, ambiguous=FALSE, kernel_calls=1L, panel_calls=1L),
    ambiguous256=list(tips=256L, sites=1000L, ambiguous=TRUE, kernel_calls=1L, panel_calls=1L),
    dense512=list(tips=512L, sites=300L, ambiguous=FALSE, kernel_calls=1L, panel_calls=1L))

if (flag('--worker')) {
    suppressPackageStartupMessages(library(genegalleon.treevis))
    suppressPackageStartupMessages(library(ggplot2))
    profile_name = option('--profile')
    profile = profiles[[profile_name]]
    if (is.null(profile)) stop('Unknown profile')
    mode = option('--mode')
    if (!mode %in% c('kernel','panel')) stop('Mode must be kernel or panel')
    source_path = normalizePath(option('--source'))
    source_sha256 = digest::digest(file=source_path,algo='sha256')
    environment = new.env(parent=asNamespace('genegalleon.treevis'))
    sys.source(source_path, envir=environment)
    if (mode == 'kernel') {
        # Run the actual producer through its final boundary, replacing only the
        # plot constructor with a data capture. All FASTA/calculation work stays.
        environment$treevis_pairwise_column = function(g,args,tips,similarities,...) {
            list(tips=as.character(tips[['label']]), y=tips[['y']], values=similarities)
        }
    }
    ids = sprintf('Species_gene%04d', seq_len(profile$tips))
    tree = ape::stree(profile$tips, type='star', tip.label=ids)
    tree$edge.length = rep(1, nrow(tree$edge))
    g = list(tree=ggtree::ggtree(tree))
    alphabet = strsplit('ACDEFGHIKLMNPQRSTVWY','',fixed=TRUE)[[1]]
    sequences = vapply(seq_len(profile$tips),function(i) {
        positions = seq_len(profile$sites)
        sequence = alphabet[(positions*7L+i*11L+(positions*i %% 23L)) %% length(alphabet)+1L]
        if (profile$ambiguous) {
            bad = (positions+i) %% 5L == 0L
            sequence[bad] = c('X','B','Z','J','O','U','*','-','?')[(positions[bad]+i) %% 9L+1L]
        }
        paste0(sequence,collapse='')
    }, character(1))
    alignment = tempfile(fileext='.fa')
    writeLines(as.vector(rbind(paste0('>',rev(ids),' description'),rev(sequences))),alignment)
    input_sha256 = digest::digest(file=alignment,algo='sha256')
    calls = profile[[paste0(mode,'_calls')]]
    invisible(gc())
    started = proc.time()
    for (iteration in seq_len(calls)) {
        result = environment$add_sequence_similarity_column(g,list(font_size=6),alignment,'protein')
    }
    timing = proc.time() - started
    seconds = timing[['elapsed']]
    cpu_seconds = sum(timing[c('user.self','sys.self')])
    status = readLines('/proc/self/status')
    peak_kib = as.numeric(sub('^VmHWM:[[:space:]]*([0-9]+).*','\\1',grep('^VmHWM:',status,value=TRUE)))
    # Retain types, row/column order, NA values, panel membership and widths.
    canonical = if (mode == 'kernel') result else list(
        names=names(result),
        panels=lapply(result[setdiff(names(result),'tree')],function(panel) list(
            data=panel$data, width_mm=attr(panel,'treevis_width_mm'), labels=as.list(panel$labels))))
    fingerprint = digest::digest(canonical, algo='sha256', serialize=TRUE)
    unlink(alignment)
    stopifnot(identical(source_sha256,digest::digest(file=source_path,algo='sha256')))
    cat(jsonlite::toJSON(list(source=source_path,source_sha256=source_sha256,input_sha256=input_sha256,
        profile=profile_name, mode=mode,
        tips=profile$tips, sites=profile$sites, ambiguous=profile$ambiguous, calls=calls,
        seconds=seconds,cpu_seconds=cpu_seconds,peak_rss_kib_linux=peak_kib,output_sha256=fingerprint,
        R=R.version.string,platform=R.version$platform,
        package_versions=as.list(vapply(c('genegalleon.treevis','ggplot2','ggtree','ape','digest','jsonlite'),
            function(package) as.character(utils::packageVersion(package)),character(1)))),
        auto_unbox=TRUE,digits=17),'\n')
    quit(status=0L)
}

output = option('--output')
if (is.null(output)) stop('--output is required')
baseline = option('--baseline-source', option('--source'))
candidate = option('--candidate-source')
if (is.null(baseline)) stop('--baseline-source or --source is required')
sources = list(baseline=normalizePath(baseline))
if (!is.null(candidate)) sources$candidate = normalizePath(candidate)
source_sha256 = lapply(sources,function(path) digest::digest(file=path,algo='sha256'))
repeats = as.integer(option('--repeats','5'))
if (!is.finite(repeats) || repeats < 1L || repeats > 10L) stop('Require 1..10 repeats')
selected = strsplit(option('--profiles','small,dense64,dense256,ambiguous256'),',',fixed=TRUE)[[1]]
if (any(!selected %in% names(profiles))) stop('Unknown profile')
samples = list()
for (profile in selected) {
    modes = if (profile == 'dense512') 'kernel' else c('kernel','panel')
    for (mode in modes) {
        expected_fingerprint = NULL
        for (trial in 0:repeats) {
            order = names(sources)
            if (trial %% 2L) order = rev(order)
            fingerprints = character()
            for (label in order) {
                command = c('--vanilla',shQuote(script),'--worker','--source',shQuote(sources[[label]]),
                    '--profile',profile,'--mode',mode)
                raw = system2(file.path(R.home('bin'),'Rscript'),command,stdout=TRUE)
                if (!is.null(attr(raw,'status'))) stop('Worker failed')
                sample = jsonlite::fromJSON(tail(raw[nzchar(raw)],1L))
                if (!identical(sample$source_sha256,source_sha256[[label]])) stop('Source changed during benchmark')
                if (is.null(expected_fingerprint)) expected_fingerprint = sample$output_sha256
                if (!identical(expected_fingerprint,sample$output_sha256)) stop('Output fingerprints differ across trials')
                fingerprints[[label]] = sample$output_sha256
                sample$label = label
                sample$trial = trial
                if (trial) samples[[length(samples)+1L]] = sample
                cat(profile,mode,if (trial) trial else 'warmup',label,
                    sample$seconds,'s',sample$peak_rss_kib_linux,'KiB\n')
            }
            if (length(unique(fingerprints)) != 1L) stop('Output fingerprints differ')
        }
    }
}
summary = list(sources=sources,source_sha256=source_sha256,repeats=repeats,
    timing_scope='Producer call including FASTA decoding; kernel captures final values, panel creates ggplot objects without rendering.',
    profiles=profiles[selected],samples=samples)
dir.create(dirname(output), recursive=TRUE, showWarnings=FALSE)
writeLines(jsonlite::toJSON(summary,auto_unbox=TRUE,pretty=TRUE,digits=17),output)
