#!/usr/bin/env python
# coding: utf-8

import argparse
import datetime
import io
import os
import shlex
import subprocess
import sys
import time

try:
    from distutils.util import strtobool
except ModuleNotFoundError:
    def strtobool(val):
        val = str(val).strip().lower()
        if val in ('y', 'yes', 't', 'true', 'on', '1'):
            return 1
        if val in ('n', 'no', 'f', 'false', 'off', '0'):
            return 0
        raise ValueError("invalid truth value {!r}".format(val))

import numpy
import pandas


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--remove_unannotated', metavar='yes|no', default='yes', type=strtobool, help='')
    parser.add_argument('--dir_orthofinder_og', metavar='PATH', default='', type=str, help='', required=True)
    parser.add_argument('--dir_species_protein', metavar='PATH', default='', type=str, help='', required=True)
    parser.add_argument('--path_diamond_db', metavar='PATH', default='', type=str, help='', required=True)
    parser.add_argument('--evalue', metavar='FLOAT', default='1e-2', type=str, help='E-value cutoff for DIAMOND BLAST search')
    parser.add_argument('--gene_size_quantiles', metavar='COMMA_SEPARATED_FLOAT', default='0.05,0.25,0.5,0.75,0.95', type=str, help='')
    parser.add_argument('--min_gene_num', metavar='INT', default=4, type=int, help='')
    parser.add_argument('--max_gene_num', metavar='INT', default=1000, type=int, help='')
    parser.add_argument('--min_species_num', metavar='INT', default=2, type=int, help='')
    parser.add_argument('--min_percent_species_coverage', metavar='FLOAT', default=50, type=int, help='')
    parser.add_argument('--ncpu', metavar='INT', default=4, type=int, help='Number of CPU threads.')
    parser.add_argument('--threads', dest='ncpu', metavar='INT', type=int, help=argparse.SUPPRESS)
    return parser


def get_pyplot():
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.rcParams['font.size'] = 8
    matplotlib.rcParams['font.family'] = 'Helvetica'
    matplotlib.rcParams['svg.fonttype'] = 'none'  # none, path, or svgfont
    return plt


def run_command(command, stdout_file=None, append=False, tool_name='tool', capture_output=False):
    print('Command:', ' '.join(shlex.quote(str(x)) for x in command))
    stdout_handle = None
    stdout_stream = None
    if stdout_file is not None:
        mode = 'ab' if append else 'wb'
        stdout_handle = open(stdout_file, mode)
    elif capture_output:
        stdout_stream = subprocess.PIPE
    try:
        out = subprocess.run(
            command,
            stdout=stdout_handle if stdout_handle is not None else stdout_stream,
            stderr=subprocess.PIPE,
            check=False,
        )
    finally:
        if stdout_handle is not None:
            stdout_handle.close()
    if out.returncode != 0:
        sys.stderr.write("{} did not finish safely.\n".format(tool_name))
        sys.stderr.write('{} stderr:\n'.format(tool_name))
        sys.stderr.write(out.stderr.decode('utf8'))
        sys.stderr.write('\n')
        sys.exit(1)
    return out.stdout if capture_output else None


def get_species_protein_files(dir_species_protein):
    species_protein_files = os.listdir(dir_species_protein)
    fasta_exts = ('.fa', '.fas', '.fasta', '.fna', '.fa.gz', '.fas.gz', '.fasta.gz', '.fna.gz')
    return [f for f in species_protein_files if f.endswith(fasta_exts)]


def get_concatenated_fx2tab(args):
    species_protein_files = get_species_protein_files(args.dir_species_protein)
    if len(species_protein_files) == 0:
        return pandas.DataFrame(columns=['#id', 'length'])
    species_protein_paths = [os.path.join(args.dir_species_protein, f) for f in species_protein_files]
    command_fx2tab = [
        'seqkit', 'fx2tab', '--threads', str(args.ncpu), '--length', '--name', '--header-line',
        '--only-id',
    ]
    command_fx2tab.extend(species_protein_paths)
    fx2tab_stdout = run_command(command_fx2tab, tool_name='seqkit', capture_output=True)
    if fx2tab_stdout is None or len(fx2tab_stdout) == 0:
        return pandas.DataFrame(columns=['#id', 'length'])
    return pandas.read_csv(io.StringIO(fx2tab_stdout.decode('utf8')), sep='\t', header=0, low_memory=False)


def get_df_gc_original(df_og_original, df_gc_original, fx2tab, args):
    quantiles = [float(val) for val in args.gene_size_quantiles.split(',')]
    quantile_cols = {quantile: 'geneid_' + str(quantile) for quantile in quantiles}
    for quantile in quantiles:
        df_gc_original[quantile_cols[quantile]] = ''

    gene2length = (
        fx2tab.loc[:, ['#id', 'length']]
        .drop_duplicates(subset=['#id'], keep='first')
        .set_index('#id')['length']
    )
    species_cols = [col for col in df_og_original.columns if col != 'Orthogroup']

    print('Identifying genes with specified quantile lengths in {:,} orthogroups'.format(df_og_original.shape[0]))
    long_frames = []
    for species_rank, col in enumerate(species_cols):
        values = df_og_original[col]
        mask = values.notna()
        if not mask.any():
            continue
        tmp = pandas.DataFrame(
            {
                'Orthogroup': df_og_original.loc[mask, 'Orthogroup'].to_numpy(),
                'genes': values.loc[mask].astype(str).to_numpy(),
                'species_rank': species_rank,
            }
        )
        long_frames.append(tmp)

    if len(long_frames) == 0:
        return df_gc_original

    long_df = pandas.concat(long_frames, ignore_index=True)
    long_df.loc[:, 'genes'] = long_df['genes'].str.strip()
    long_df = long_df.loc[long_df['genes'] != '', :]
    if long_df.empty:
        return df_gc_original

    long_df.loc[:, 'cell_idx'] = numpy.arange(long_df.shape[0], dtype=numpy.int64)
    long_df.loc[:, 'gene_id'] = long_df['genes'].str.split(',')
    long_df = long_df.explode('gene_id', ignore_index=True)
    long_df.loc[:, 'gene_id'] = long_df['gene_id'].astype(str).str.strip()
    long_df = long_df.loc[long_df['gene_id'] != '', :]
    if long_df.empty:
        return df_gc_original

    long_df.loc[:, 'gene_pos'] = long_df.groupby('cell_idx', sort=False).cumcount()
    long_df = long_df.sort_values(['Orthogroup', 'species_rank', 'gene_pos'], kind='mergesort')
    long_df.loc[:, 'length'] = long_df['gene_id'].map(gene2length)

    total_counts = long_df.groupby('Orthogroup', sort=False).size()
    matched_mask = long_df['length'].notna()
    matched_counts = long_df.loc[matched_mask, :].groupby('Orthogroup', sort=False).size()
    matched_counts = matched_counts.reindex(total_counts.index).fillna(0).astype(int)
    mismatched = total_counts.index[total_counts.to_numpy() != matched_counts.to_numpy()]
    for orthogroup in mismatched:
        num_total = int(total_counts.loc[orthogroup])
        num_matched = int(matched_counts.loc[orthogroup])
        txt = 'Only {:,} out of {:,} genes were matched between fx2tab and Orthogroups.tsv: {}\n'
        sys.stderr.write(txt.format(num_matched, num_total, orthogroup))

    matched = long_df.loc[matched_mask, ['Orthogroup', 'gene_id', 'length']]
    if matched.empty:
        return df_gc_original

    quantile_values = (
        matched.groupby('Orthogroup', sort=False)['length']
        .quantile(q=quantiles, interpolation='nearest')
        .unstack(level=-1)
        .reindex(columns=quantiles)
    )
    selected = pandas.DataFrame(index=quantile_values.index)
    matched_og = matched['Orthogroup']
    matched_length = matched['length'].to_numpy(copy=False)
    for quantile in quantiles:
        target_len = matched_og.map(quantile_values[quantile]).to_numpy()
        selected_rows = matched.loc[matched_length == target_len, ['Orthogroup', 'gene_id']]
        selected_gene = selected_rows.drop_duplicates(subset='Orthogroup', keep='first').set_index('Orthogroup')['gene_id']
        selected[quantile_cols[quantile]] = selected.index.to_series().map(selected_gene)
    for quantile in quantiles:
        col = quantile_cols[quantile]
        df_gc_original.loc[:, col] = df_gc_original['Orthogroup'].map(selected[col]).fillna(df_gc_original[col])
    return df_gc_original


def _attach_besthits(df_gc_original_quartile, df_diamond, cols):
    out = df_gc_original_quartile.copy()
    if df_diamond.empty:
        for col in cols:
            out[col.replace('geneid_', 'besthit_')] = ''
        return out
    besthit_map = (
        df_diamond.loc[:, ['qseqid', 'stitle']]
        .drop_duplicates(subset=['qseqid'], keep='first')
        .set_index('qseqid')['stitle']
    )
    for col in cols:
        out[col.replace('geneid_', 'besthit_')] = out[col].map(besthit_map)
    return out


def diamond_annotation(df_gc_original_quartile, args):
    tmp_query_fasta = 'tmp.query.fasta'
    tmp_query_list = 'tmp.query.txt'
    tmp_diamond_out = 'tmp.diamond.out.tsv'
    if os.path.exists(tmp_query_fasta):
        os.remove(tmp_query_fasta)
    with open(tmp_query_fasta, 'w'):
        pass
    cols = ['geneid_' + str(val) for val in args.gene_size_quantiles.split(',')]
    raw_query_ids = pandas.unique(df_gc_original_quartile.loc[:, cols].to_numpy().ravel())
    query_ids = [str(x) for x in raw_query_ids if pandas.notna(x) and str(x) != '']
    numpy.savetxt(tmp_query_list, query_ids, fmt='%s')
    species_protein_files = get_species_protein_files(args.dir_species_protein)
    if (not os.path.exists(tmp_diamond_out)) and (len(query_ids) > 0):
        print('Generating {}'.format(tmp_diamond_out))
        for species_protein_file in species_protein_files:
            species_protein_path = os.path.join(args.dir_species_protein, species_protein_file)
            command_seqkit = [
                'seqkit', 'grep', '--threads', str(args.ncpu), '--pattern-file', tmp_query_list, species_protein_path
            ]
            run_command(command_seqkit, stdout_file=tmp_query_fasta, append=True, tool_name='seqkit')
        command_diamond = [
            'diamond', 'blastp', '--query', tmp_query_fasta, '--threads', str(args.ncpu),
            '--db', args.path_diamond_db, '--outfmt', '6', 'qseqid', 'stitle', 'evalue', 'bitscore',
            '--max-target-seqs', '1', '--evalue', args.evalue, '--out', tmp_diamond_out
        ]
        run_command(command_diamond, tool_name='diamond')
    if not os.path.exists(tmp_diamond_out):
        return _attach_besthits(df_gc_original_quartile, pandas.DataFrame(columns=['qseqid', 'stitle']), cols)
    print('Loading {}'.format(tmp_diamond_out))
    try:
        df_diamond = pandas.read_csv(tmp_diamond_out, sep='\t', header=None, low_memory=False)
    except pandas.errors.EmptyDataError:
        df_diamond = pandas.DataFrame(columns=['qseqid', 'stitle'])
    if not df_diamond.empty:
        if df_diamond.shape[1] < 2:
            df_diamond = pandas.DataFrame(columns=['qseqid', 'stitle'])
        else:
            df_diamond = df_diamond.iloc[:, :2].copy()
            df_diamond.columns = ['qseqid', 'stitle']
    return _attach_besthits(df_gc_original_quartile, df_diamond, cols)


def prepare_annotation(args):
    file_genecount_annot = os.path.join(args.dir_orthofinder_og, 'Orthogroups.GeneCount.annotated.tsv')
    if os.path.exists(file_genecount_annot):
        print('Skipping annotation. Annotated file detected: {}'.format(file_genecount_annot))
        return file_genecount_annot

    print('Started functional annotation of one representative gene per orthogroup.')
    fx2tab = get_concatenated_fx2tab(args)
    if fx2tab.empty:
        sys.stderr.write('No fx2tab records were generated from species protein files. Exiting.\n')
        sys.exit(1)

    quartile_path = 'tmp.Orthogroups.GeneCount.quartile_genes.tsv'
    if os.path.exists(quartile_path):
        print('Loading: {}'.format(quartile_path))
        df_gc_original_quartile = pandas.read_csv(quartile_path, header=0, sep='\t', low_memory=False)
    else:
        print('Generating {}'.format(quartile_path))
        file_og_in = os.path.join(args.dir_orthofinder_og, 'Orthogroups.tsv')
        df_og_original = pandas.read_csv(file_og_in, header=0, sep='\t', low_memory=False)
        file_genecount_in = os.path.join(args.dir_orthofinder_og, 'Orthogroups.GeneCount.tsv')
        df_gc_original = pandas.read_csv(file_genecount_in, header=0, sep='\t', low_memory=False)
        df_gc_original_quartile = get_df_gc_original(df_og_original, df_gc_original, fx2tab, args)
        df_gc_original_quartile.to_csv(quartile_path, sep='\t', index=False)
    df_gc_original_annotated = diamond_annotation(df_gc_original_quartile, args)
    df_gc_original_annotated.to_csv(file_genecount_annot, sep='\t', index=False)
    return file_genecount_annot


def select_orthogroups(args, file_genecount_in):
    file_genecount_out = os.path.join(args.dir_orthofinder_og, 'Orthogroups.GeneCount.selected.tsv')
    df_gc_original = pandas.read_csv(file_genecount_in, header=0, sep='\t', low_memory=False)
    print('Number of orthogroups in Orthogroups.GeneCount.tsv: {:,}'.format(df_gc_original.shape[0]))
    is_big_enough = (df_gc_original['Total'] >= args.min_gene_num)
    is_small_enough = (df_gc_original['Total'] <= args.max_gene_num)
    print('Number of excluded orthogroups with less than {:,} genes: {:,}'.format(args.min_gene_num, (~is_big_enough).sum()))
    print('Number of excluded orthogroups with more than {:,} genes: {:,}'.format(args.max_gene_num, (~is_small_enough).sum()))
    if args.remove_unannotated:
        cols = df_gc_original.columns.str.startswith('besthit_')
        is_annotated = (~((df_gc_original.loc[:, cols] == '') | (df_gc_original.loc[:, cols].isnull()))).any(axis=1)
        print('Number of excluded orthogroups with no UniProt DIAMOND BLASTP hit: {:,}'.format((~is_annotated).sum()))
    else:
        print('Annotation info was not used for the orthogroup selection.')
        is_annotated = True

    ex_cols = ['Orthogroup', 'Total']
    ex_cols += df_gc_original.columns[df_gc_original.columns.str.startswith('geneid_')].tolist()
    ex_cols += df_gc_original.columns[df_gc_original.columns.str.startswith('besthit_')].tolist()
    sp_cols = [col for col in df_gc_original.columns if col not in ex_cols]

    num_sp = (df_gc_original.loc[:, sp_cols] >= 1).sum(axis=1)
    percent_sp = num_sp / len(sp_cols) * 100
    is_enough_num_species = (num_sp >= args.min_species_num)
    is_enough_percent_species = (percent_sp >= args.min_percent_species_coverage)
    print('Number of excluded orthogroups with less than {:,} species: {:,}'.format(args.min_species_num, (~is_enough_num_species).sum()))
    print('Number of excluded orthogroups with less than {:,}% of species: {:,}'.format(args.min_percent_species_coverage, (~is_enough_percent_species).sum()))

    is_selected = (is_big_enough & is_small_enough & is_annotated & is_enough_num_species & is_enough_percent_species)
    df_gc = df_gc_original.loc[is_selected, :].reset_index(drop=True)
    print('Number of orthogroups in Orthogroups.GeneCount.selected.tsv: {:,}'.format(df_gc.shape[0]))
    df_gc.to_csv(file_genecount_out, sep='\t', index=False)

    file_og_in = os.path.join(args.dir_orthofinder_og, 'Orthogroups.tsv')
    file_og_out = os.path.join(args.dir_orthofinder_og, 'Orthogroups.selected.tsv')
    df_og_original = pandas.read_csv(file_og_in, header=0, sep='\t', low_memory=False)
    print('Number of orthogroups in Orthogroups.tsv: {:,}'.format(df_og_original.shape[0]))
    df_og = df_og_original.loc[(df_og_original['Orthogroup'].isin(df_gc['Orthogroup'])), :]
    df_og = df_og.reset_index(drop=True)
    print('Number of orthogroups in Orthogroups.selected.tsv: {:,}'.format(df_og.shape[0]))
    df_og.to_csv(file_og_out, sep='\t', index=False)

    return df_gc_original, df_gc, sp_cols


def print_selection_stats(args, df_gc_original, df_gc, sp_cols):
    print('')
    for sp_col in sp_cols:
        num_sp_gene_before = df_gc_original.loc[:, sp_col].sum()
        num_sp_gene_after = df_gc.loc[:, sp_col].sum()
        txt = 'Number of genes before and after selection: {:,} and {:,} in {}'
        print(txt.format(num_sp_gene_before, num_sp_gene_after, sp_col))
    print('')
    print('Minimum gene number cutoff: {:,}'.format(args.min_gene_num))
    print('Minimum species number cutoff: {:,}'.format(args.min_species_num))
    print('Minimum % species_coverage cutoff: {:,}'.format(args.min_percent_species_coverage))
    print('Number of selected orthogroups: {:,}'.format(df_gc.shape[0]))
    num_single_sp = (df_gc.loc[:, sp_cols] == 1).sum(axis=1)
    num_single_og = (num_single_sp == len(sp_cols)).sum()
    print('Number of selected strictly single-copy orthogroups: {:,}'.format(num_single_og))
    num_sp = (df_gc.loc[:, sp_cols] >= 1).sum(axis=1)
    num_nomissing_og = (num_sp == len(sp_cols)).sum()
    print('Number of selected orthogroups without missing species: {:,}'.format(num_nomissing_og))
    txt = 'In gg_gene_evolution_job.sh, set the array job range of 1-{} to examine all selected orthogroups.'
    print(txt.format(df_gc.shape[0]))


def plot_gene_number(args, df_gc_original, df_gc):
    before_sorted = df_gc_original.sort_values(by='Total', ascending=False).reset_index(drop=True)
    after_sorted = df_gc.sort_values(by='Total', ascending=False).reset_index(drop=True)
    plt = get_pyplot()
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(3.2, 3.2), sharey=False, sharex=False)
    ax = axes
    ax.plot(before_sorted.index + 1, before_sorted['Total'], color='gray', label='Before selection')
    ax.plot(after_sorted.index + 1, after_sorted['Total'], color='black', label='After selection')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Orthogroup')
    ax.set_ylabel('Number of genes')
    ax.legend()
    fig.tight_layout()
    outbase = os.path.join(args.dir_orthofinder_og, 'gene_number')
    for ext in ['svg', 'pdf']:
        outpath = os.path.join(outbase + '.' + ext)
        fig.savefig(outpath, format=ext)


def remove_tmp_files():
    for tmp_file in [f for f in os.listdir(os.getcwd()) if f.startswith('tmp.')]:
        print('Removing: {}'.format(tmp_file))
        os.remove(tmp_file)


def run(args):
    file_genecount_in = prepare_annotation(args)
    df_gc_original, df_gc, sp_cols = select_orthogroups(args, file_genecount_in)
    print_selection_stats(args, df_gc_original, df_gc, sp_cols)
    plot_gene_number(args, df_gc_original, df_gc)
    remove_tmp_files()


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    args.ncpu = max(1, int(args.ncpu))
    start = time.time()
    print('Starting {} at {}'.format(sys.argv[0], datetime.datetime.now()))
    run(args)
    print('Ending {} at {}. Elapsed time: {:,} sec'.format(sys.argv[0], datetime.datetime.now(), int(time.time() - start)))


if __name__ == '__main__':
    main()
