# Shared alignment boundary + coding phase denotes a positional correspondence
# candidate, not proof of common ancestry. Gap-spanning boundaries stay unresolved.
treevis_intron_site_data = function(tips, seqs) {
    seqs = seqs[names(seqs) %in% as.character(tips[['label']])]
    if (!length(seqs)) return(NULL)
    if (any(nchar(seqs) == 0) || length(unique(nchar(seqs))) != 1) stop('Intron sites require an equal-length untrimmed CDS alignment')
    if (any(grepl('[^ACGTRYSWKMBDHVN?.-]', seqs))) stop('Intron sites require nucleotide sequences')
    events = list()
    rows = list()
    diagnostics = list()
    for (i in seq_len(nrow(tips))) {
        id = as.character(tips[['label']][i])
        reason = NA_character_
        if (!(id %in% names(seqs))) reason = 'missing_alignment'
        count = tips[['num_intron']][i]
        length_cds = tips[['intron_feature_size']][i]
        phase = tips[['cds_first_phase']][i]
        if (is.na(count) || is.na(length_cds)) reason = 'missing_GFF'
        if (is.na(reason)) {
            chars = strsplit(seqs[[id]], '', fixed=TRUE)[[1]]
            nongap = which(!chars %in% c('-', '.'))
            if (length(nongap) != length_cds) reason = 'CDS_length_mismatch'
        }
        if (!is.na(reason)) {
            diagnostics[[length(diagnostics)+1]] = data.frame(node_name=id, reason=reason)
            next
        }
        if (!is.finite(count) || count < 0 || count != floor(count) || (!is.na(phase) && !phase %in% 0:2)) stop('Invalid intron count or CDS phase')
        value = as.character(tips[['intron_positions']][i])
        offsets = if (is.na(value) || !nzchar(value)) numeric() else suppressWarnings(as.numeric(strsplit(value, '[;,]')[[1]]))
        if (length(offsets) != count || any(!is.finite(offsets)) || any(offsets != floor(offsets)) ||
            any(offsets < 1 | offsets >= length_cds) || anyDuplicated(offsets) || is.unsorted(offsets)) stop('Invalid intron offsets for ', id)
        rows[[id]] = list(chars=chars, cumulative=cumsum(!chars %in% c('-', '.')), phase=phase, y=tips[['y']][i],
            colour=if ('tiplab_color' %in% names(tips)) as.character(tips[['tiplab_color']][i]) else 'black')
        for (k in seq_along(offsets)) {
            offset = offsets[k]
            left = nongap[offset]; right = nongap[offset+1]
            status = if (right != left + 1) 'alignment_gap' else if (!all(chars[c(left,right)] %in% c('A','C','G','T'))) 'ambiguous_bases' else if (is.na(phase)) 'position_only' else 'mapped'
            events[[length(events)+1]] = data.frame(node_name=id, intron_index=k, cds_offset=offset,
                alignment_left=left, alignment_right=right, phase=(offset-phase) %% 3, status=status)
        }
    }
    events = if (length(events)) do.call(rbind, events) else data.frame(node_name=character(), intron_index=integer(),
        cds_offset=numeric(), alignment_left=integer(), alignment_right=integer(), phase=numeric(), status=character())
    keys = unique(events[events$status == 'mapped',c('alignment_left','phase'),drop=FALSE])
    unknown = unique(events[events$status == 'position_only',c('alignment_left','phase'),drop=FALSE])
    keys = rbind(keys, unknown[!unknown$alignment_left %in% keys$alignment_left,,drop=FALSE])
    keys = keys[order(keys$alignment_left,keys$phase),,drop=FALSE]
    keys$site_id = sprintf('I%03d',seq_len(nrow(keys)))
    key = function(left,phase) paste(left,phase,sep=':')
    events$site_id = keys$site_id[match(key(events$alignment_left,events$phase),key(keys$alignment_left,keys$phase))]
    events$site_id[events$status != 'mapped'] = NA_character_
    events$candidate_site_ids = rep(NA_character_, nrow(events))
    for (k in which(events$status == 'position_only')) {
        candidates = keys$site_id[keys$alignment_left == events$alignment_left[k]]
        events$candidate_site_ids[k] = paste(candidates,collapse=';')
        if (length(candidates) == 1) events$site_id[k] = candidates
    }
    cells = list()
    for (id in names(rows)) {
        row = rows[[id]]
        for (j in seq_len(nrow(keys))) {
            b = keys$alignment_left[j]
            available = all(row$chars[c(b,b+1)] %in% c('A','C','G','T')) &&
                (is.na(row$phase) || is.na(keys$phase[j]) || (row$cumulative[b]-row$phase) %% 3 == keys$phase[j])
            present = any(events$node_name == id & events$status == 'mapped' & !is.na(events$site_id) & events$site_id == keys$site_id[j])
            possible = any(events$node_name == id & events$status == 'position_only' & events$alignment_left == b)
            cells[[length(cells)+1]] = data.frame(node_name=id,y=row$y,x=j,site_id=keys$site_id[j],
                state=if (present) 'present' else if (possible) 'position_only' else if (available) 'absent' else 'unresolved',colour=row$colour)
        }
    }
    list(sites=keys, events=events, cells=if (length(cells)) do.call(rbind,cells) else data.frame(),
        diagnostics=if (length(diagnostics)) do.call(rbind,diagnostics) else data.frame(node_name=character(),reason=character()))
}
