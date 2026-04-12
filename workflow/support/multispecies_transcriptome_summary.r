# %%
cat('Starting multispecies_transcriptome_summary.r\n')
cli_args = commandArgs(trailingOnly=TRUE)
args = cli_args
cat(paste0('Working at: ', getwd(), '\n'))

library(cowplot, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(missMDA, quietly=TRUE)
library(svglite, quietly=TRUE)
library(rkftools, quietly=TRUE)
#library(readxl, quietly=TRUE)
has_rtsne = requireNamespace("Rtsne", quietly=TRUE)

options(stringsAsFactors=FALSE)

cat('arguments:\n')
args = rkftools::get_parsed_args(args, print=TRUE)

font_size = 8
font_size_factor = 0.352777778 # https://stackoverflow.com/questions/17311917/ggplot2-the-unit-of-size
dir_assembly_stat = args[['dir_assembly_stat']]
dir_amalgkit_metadata = args[['dir_amalgkit_metadata']]
dir_amalgkit_merge = args[['dir_amalgkit_merge']]
dir_busco_isoform = args[['dir_busco_isoform']]
dir_busco_longest_cds = args[['dir_busco_longest_cds']]

normalize_busco_species_colname = function(x) {
    x = sub('\\.tsv$', '', x)
    x = sub('([._]busco([._]full)?)$', '', x, perl=TRUE)
    x = sub('([._]full)$', '', x, perl=TRUE)
    x
}

coerce_numeric_data_frame = function(df_in) {
    df_out = as.data.frame(lapply(df_in, function(x) {
        if (is.numeric(x)) {
            return(as.numeric(x))
        }
        suppressWarnings(as.numeric(as.character(x)))
    }), check.names=FALSE, stringsAsFactors=FALSE)
    rownames(df_out) = rownames(df_in)
    df_out
}

is_zero_variance_vector = function(x) {
    x = x[!is.na(x)]
    if (length(x) <= 1) {
        return(FALSE)
    }
    stats::var(x) == 0
}

filter_expression_for_dimensional_reduction = function(df_in, min_non_missing = 2) {
    df_out = coerce_numeric_data_frame(df_in)
    min_non_missing = max(2L, as.integer(min_non_missing))
    mat = as.matrix(df_out)
    inf_mask = is.infinite(mat)
    num_infinite = sum(inf_mask)
    if (num_infinite > 0) {
        cat(sprintf('Converting %d Inf/-Inf value(s) to NA before dimensional reduction.\n', num_infinite))
        mat[inf_mask] = NA_real_
        df_out = as.data.frame(mat, check.names=FALSE, stringsAsFactors=FALSE)
        rownames(df_out) = rownames(df_in)
    }

    repeat {
        if (nrow(df_out) == 0 || ncol(df_out) == 0) {
            break
        }
        row_all_na = apply(df_out, 1, function(x) all(is.na(x)))
        col_all_na = vapply(df_out, function(x) all(is.na(x)), logical(1))
        row_zero_var = apply(df_out, 1, is_zero_variance_vector)
        col_zero_var = vapply(df_out, is_zero_variance_vector, logical(1))
        row_too_sparse = apply(df_out, 1, function(x) sum(!is.na(x)) < min_non_missing)
        col_too_sparse = vapply(df_out, function(x) sum(!is.na(x)) < min_non_missing, logical(1))
        if (!(any(row_all_na) || any(col_all_na) || any(row_zero_var) || any(col_zero_var) || any(row_too_sparse) || any(col_too_sparse))) {
            break
        }
        if (any(row_all_na)) {
            cat(sprintf('Removing %d all-NA gene row(s) before dimensional reduction.\n', sum(row_all_na)))
        }
        if (any(col_all_na)) {
            cat(sprintf('Removing %d all-NA species column(s) before dimensional reduction.\n', sum(col_all_na)))
        }
        if (any(row_zero_var)) {
            cat(sprintf('Removing %d zero-variance gene row(s) before dimensional reduction.\n', sum(row_zero_var)))
        }
        if (any(col_zero_var)) {
            cat(sprintf('Removing %d zero-variance species column(s) before dimensional reduction.\n', sum(col_zero_var)))
        }
        if (any(row_too_sparse)) {
            cat(sprintf('Removing %d gene row(s) with <%d observed species before dimensional reduction.\n', sum(row_too_sparse), min_non_missing))
        }
        if (any(col_too_sparse)) {
            cat(sprintf('Removing %d species column(s) with <%d observed genes before dimensional reduction.\n', sum(col_too_sparse), min_non_missing))
        }
        if (any(row_all_na) || any(row_zero_var) || any(row_too_sparse)) {
            df_out = df_out[!(row_all_na | row_zero_var | row_too_sparse), , drop=FALSE]
        }
        if (nrow(df_out) == 0 || ncol(df_out) == 0) {
            break
        }
        if (any(col_all_na) || any(col_zero_var) || any(col_too_sparse)) {
            df_out = df_out[, !(col_all_na | col_zero_var | col_too_sparse), drop=FALSE]
        }
    }

    df_out
}

run_expression_imputation = function(df_in) {
    tryCatch({
        nb = missMDA::estim_ncpPCA(df_in, method.cv="Kfold", nbsim=100, threshold=1e-2, verbose=TRUE)
        cat(paste0('Number of components for expression data imputations: ', nb[['ncp']], '\n'))
        missMDA::imputePCA(df_in, ncp=nb[['ncp']])
    }, error=function(e) {
        cat(sprintf(
            'Skipping transcriptome dimensionality reduction because missMDA failed on the filtered expression matrix: %s\n',
            conditionMessage(e)
        ))
        NULL
    })
}

# %%
cat('Starting transcriptome dimensionality reduction plot.\n')

cat('Generating busco_table.tsv\n')
my_command = paste('python', shQuote(file.path(args[['dir_myscript']], 'collect_common_BUSCO_genes.py')))
my_command = paste(my_command, '--busco_outdir', shQuote(dir_busco_longest_cds))
my_command = paste(my_command, '--outfile', 'busco_table.tsv')
cat(paste0('Command: ', my_command, '\n'))
if (file.exists('busco_table.tsv')) {
    file.remove('busco_table.tsv')
}
busco_status = system(my_command)
if (busco_status != 0) {
    stop(sprintf('collect_common_BUSCO_genes.py failed with exit code %s', busco_status))
}
df_busco = read.table('busco_table.tsv', header=TRUE, sep='\t', quote='')
if (ncol(df_busco) >= 4) {
    original_busco_cols = colnames(df_busco)[4:ncol(df_busco)]
    normalized_busco_cols = vapply(original_busco_cols, normalize_busco_species_colname, character(1))
    normalized_busco_cols = make.unique(normalized_busco_cols, sep='__dup')
    for (i in seq_along(original_busco_cols)) {
        if (original_busco_cols[i] != normalized_busco_cols[i]) {
            cat(sprintf('Normalized BUSCO column: %s -> %s\n', original_busco_cols[i], normalized_busco_cols[i]))
        }
    }
    colnames(df_busco)[4:ncol(df_busco)] = normalized_busco_cols
}

cat('Generating expression.tsv\n')
if (file.exists('expression.tsv')) {
    file.remove('expression.tsv')
}
tpm_files = list.files(path=dir_amalgkit_merge, pattern="_tpm\\.tsv$", recursive=TRUE, full.names=TRUE)
if (length(tpm_files)==0) {
    cat(paste0('Skipping. No tsv file found in: ', dir_amalgkit_merge, '\n'))
} else {
    tpm_list = list()
    cat(paste('Processing', length(tpm_files), 'TPM files.\n'))
    for (i in 1:length(tpm_files)) {
        file_path = tpm_files[i]
        cat('Processing file', i, 'of', length(tpm_files), 'TPM files:', file_path, '\n')
        sp_ub = tpm_files[i]
        sp_ub = sub('^.*/', '', sp_ub)
        sp_ub = sub('_tpm.tsv$', '', sp_ub)
        #sp_ub = sub('_', 'PLACEHOLDER', sp_ub)
        #sp_ub = sub('_.*', '', sp_ub)
        #sp_ub = sub('PLACEHOLDER', '_', sp_ub)
        cat('Species unique base:', sp_ub, '\n')
        if (!(sp_ub %in% colnames(df_busco))) {
            cat('Skipping TPM file because BUSCO column is missing for species:', sp_ub, '\n')
            next
        }
        tpm = read.table(file_path, header=TRUE, sep='\t')
        tpm[,'target_id'] = sub('-i[0-9]+$', '', tpm[['target_id']])
        colnames(tpm) = sub('target_id', sp_ub, colnames(tpm))
        retained_cols = colnames(tpm)[2:length(colnames(tpm))]
        tmp_busco = data.frame(df_busco[,sp_ub])
        colnames(tmp_busco) = sp_ub
        tmp_busco[,'row_index'] = 1:nrow(tmp_busco)
        tpm = merge(tmp_busco, tpm, by=sp_ub, all.x=TRUE, all.y=FALSE, sort=FALSE)
        tpm = tpm[order(tpm[['row_index']]),]
        tpm = data.frame(tpm[,retained_cols])
        colnames(tpm) = retained_cols
        tpm_list[[length(tpm_list) + 1]] = tpm
    }
    cat('\n')
    if (length(tpm_list) == 0) {
        cat('Skipping. No TPM files matched BUSCO species columns.\n')
    } else {
        df_expression = do.call(cbind, tpm_list)
        df_expression = cbind(df_busco[,c('busco_id','orthodb_url', 'description')], df_expression)
        write.table(df_expression, 'expression.tsv', row.names=FALSE, sep='\t', quote=FALSE)
    }
}
min_species = 4
min_gene = 4

if (file.exists('expression.tsv')) {
    df_exp = read.table('expression.tsv', header=TRUE, sep='\t', quote='')
    df_exp = df_exp[,-c(1,2,3)]
    is_gene_with_ge_min_species = apply(df_exp, 1, function(x){((length(x)-sum(is.na(x)))<min_species)})
    num_enough_gene = nrow(df_exp) - sum(is_gene_with_ge_min_species)
    is_no_enough_data = (num_enough_gene <= min_gene)
    cat('Number of all genes:', nrow(df_exp), '\n')
    cat(paste0('Number of genes with >=', min_species, ' non-missing species: ', num_enough_gene, '\n'))
    if (is_no_enough_data) {
        cat('Expecting at least', min_gene, 'analyzable genes. Skipping.\n')
    } else {
        df_exp_filtered = filter_expression_for_dimensional_reduction(df_exp, min_non_missing=min_species)
        cat(sprintf(
            'Expression matrix retained %d gene row(s) and %d species column(s) after filtering.\n',
            nrow(df_exp_filtered), ncol(df_exp_filtered)
        ))
        if (nrow(df_exp_filtered) < 2 || ncol(df_exp_filtered) < 2) {
            cat(sprintf(
                'Skipping dimensional reduction because too few rows/columns remain after filtering: %d gene row(s), %d species column(s).\n',
                nrow(df_exp_filtered), ncol(df_exp_filtered)
            ))
        } else {
            cat('Generating expression.imputed.tsv\n')
            print(paste('Starting the expression level imputation:', Sys.time()))
            res_comp = run_expression_imputation(df_exp_filtered)
            if (!is.null(res_comp)) {
                imp = res_comp[['completeObs']]
                df_imp = as.data.frame(imp, check.names=FALSE, stringsAsFactors=FALSE)
                write.table(df_imp, 'expression.imputed.tsv', row.names=FALSE, sep='\t', quote=FALSE)
                print(paste('Ending the expression level imputation:', Sys.time()))

                for (method in c('pca', 'tsne', 'mds')) {
                cat(paste0('Starting dimensionality reduction: ', method, '\n'))
                input_data = df_imp
                if (method=='pca') {
                    if (nrow(input_data) < 2 || ncol(input_data) < 3) {
                        cat(sprintf(
                            'Skipping PCA plot because at least 2 gene rows and 3 species columns are required after filtering; observed %d gene row(s) and %d species column(s).\n',
                            nrow(input_data), ncol(input_data)
                        ))
                        next
                    }
                    pca_res = prcomp(t(input_data), scale.=FALSE)
                    if (ncol(pca_res[['x']]) < 2) {
                        cat(sprintf('Skipping PCA plot because only %d principal component(s) are available.\n', ncol(pca_res[['x']])))
                        next
                    }
                    df = data.frame(pca_res[['x']][,c(1,2)])
                    xlabel = paste0("PC 1 (", round(summary(pca_res)$importance[2,1]*100, digits=1), "%)")
                    ylabel = paste0("PC 2 (", round(summary(pca_res)$importance[2,2]*100, digits=1), "%)")
                } else if (method=='mds') {
                    if (ncol(input_data) < 3) {
                        cat(sprintf(
                            'Skipping MDS plot because at least 3 species columns are required after filtering; observed %d.\n',
                            ncol(input_data)
                        ))
                        next
                    }
                    corr_matrix = stats::cor(input_data, method = 'pearson', use = 'pairwise.complete.obs')
                    corr_matrix[is.na(corr_matrix)] = 0
                    diag(corr_matrix) = 1
                    tc_dist_dist = as.dist(1 - corr_matrix + 1e-09)
                    mds = stats::cmdscale(tc_dist_dist, k = 2)
                    df = data.frame(mds[, c(1,2)])
                    rownames(df) = colnames(input_data)
                    xlabel = 'MDS 1'
                    ylabel = 'MDS 2'
                } else if (method=='tsne') {
                    if (!has_rtsne) {
                        cat('Skipping tSNE plot because Rtsne is unavailable.\n')
                        next
                    }
                    perplexity = min(30, floor((ncol(input_data) - 1) / 3))
                    if (perplexity < 1) {
                        cat(sprintf(
                            'Skipping tSNE plot because at least 4 species columns are required after filtering; observed %d.\n',
                            ncol(input_data)
                        ))
                        next
                    }
                    out_tsne = Rtsne::Rtsne(as.matrix(t(input_data)), theta=0, check_duplicates=FALSE, 
                                            verbose=FALSE, dims=2, perplexity=perplexity)
                    df = data.frame(out_tsne[['Y']])
                    rownames(df) = colnames(input_data)
                    xlabel = 'tSNE 1'
                    ylabel = 'tSNE 2'
                }
                colnames(df) = c('x','y')
                df[,'label'] = rownames(df)

                g = ggplot() +
                    geom_text(data=df, mapping=aes(x=x, y=y, label=label), size=font_size*font_size_factor) +
                    xlab(xlabel) +
                    ylab(ylabel) +
                    theme_bw() +
                    theme(
                        axis.text=element_text(size=font_size),
                        axis.title=element_text(size=font_size),
                        panel.grid.major.y=element_blank(),
                        panel.grid.major.x=element_blank(),
                        panel.grid.minor.y=element_blank(),
                        panel.grid.minor.x=element_blank(),
                        #axis.title.y=element_blank(), 
                        #axis.text.y=element_blank(), 
                        #axis.ticks.y=element_blank(),
                        legend.title=element_blank(),
                        legend.text=element_text(size=font_size),
                        legend.position="inside",
                        legend.position.inside=c(0.1,0.9),
                        #legend.key.size=unit(0.4, 'lines'), 
                        #legend.box.just=0.5,
                        rect=element_rect(fill="transparent"),
                        plot.margin=unit(rep(0.1, 4), "cm")
                    )
                for (ext in c('pdf','svg')) {
                    file_name = paste0("multisp_", method, '.', ext)
                    ggsave(file_name, width = 7.2, height = 7.2)
                }
                g
                }
            }
        }
    }
}
cat('Ending transcriptome dimensionality reduction plot.\n')

# %%
cat('Starting assembly stat summary.\n')
file_base = 'assembly_stat_summary'

stat_files = list.files(dir_assembly_stat)
stat_files = stat_files[endsWith(stat_files, '.tsv')]
if (length(stat_files)==0) {
    cat(paste0('Skipping. No tsv file found in: ', dir_assembly_stat, '\n'))
} else {
    cat(paste0(length(stat_files), ' tsv file found in: ', dir_assembly_stat, '\n'))
    df_list = list()
    for (i in 1:length(stat_files)) {
        file_path = file.path(dir_assembly_stat, stat_files[i])
        df_list[[i]] = read.table(file_path, header=TRUE, sep='\t')
    }
    df = do.call(rbind, df_list)
    original_cols = colnames(df)
    df[,'species'] = sub('\\..*', '', sub('.*/', '', df[['file']]))
    df[,'species'] = sub('_longestCDS', '', df[['species']])
    df[,'species'] = sub('_', ' ', df[['species']])
    df[,'assembly_type'] = ''
    df[grepl('isoform', df[['file']]),'assembly_type'] = 'isoform'
    df[grepl('longestCDS', df[['file']]),'assembly_type'] = 'longestCDS'
    ordered_cols = c('species', 'assembly_type', original_cols[3:length(original_cols)], original_cols[1:2])
    df = df[,ordered_cols]
    df = df[order(df[['assembly_type']], df[['species']]),]
    write.table(df, paste0(file_base, '.tsv'), row.names=FALSE, sep='\t', quote=FALSE)

    g = ggplot(data=df) +
        geom_point(aes(x=log10(num_seqs+1), y=species), alpha=1) +
        xlim(0, ceiling(max(log10(df[['num_seqs']]+1)))) +
        ylab('') +
        xlab('log10(number of sequences + 1)') +
        facet_grid(. ~ assembly_type) +
        theme_bw() +
        theme(
            #panel.grid.major.x=element_blank(),
            panel.grid.minor.x=element_blank(),
            #panel.grid.major.y=element_blank(),
            #panel.grid.minor.y=element_blank(),
            text = element_text(size=font_size, color='black'),
            axis.text=element_text(size=font_size, color='black'),
            strip.text = element_text(size=font_size, color='black'),
            axis.title=element_text(size=font_size, color='black'),
            legend.text=element_text(size=font_size, color='black'),
            legend.title=element_text(size=font_size, face='bold', hjust=0.5),
            legend.background = element_rect(fill=alpha('white', 0.4)),
            legend.box.just='bottom',
        )
    extensions = c('.pdf','.svg')
    for (extension in extensions) {
        cowplot::save_plot(
            filename=paste0(file_base, extension), 
            plot=g, 
            nrow=1,
            base_height=max(1.5, nrow(df)/12),
            base_width=7.2,
            units='in', 
            dpi=300, 
            limitsize=FALSE,
            bg = 'transparent'
        )
    }
    cat('\n')    
    
    g
}


# %%
cat('Skipping legacy SNV summary: feature has been removed from the pipeline.\n')

# %%
cat('Ending multispecies_transcriptome_summary.r\n')

# %%
