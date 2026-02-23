#!/usr/bin/env python3

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime
from functools import lru_cache
import gzip
import numpy
import os
import pandas
import re
import sys
import time

pandas.options.mode.chained_assignment = None

def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', metavar='intron_num|gene_delim', default='', type=str, help='',
                        choices=['intron_num','gene_delim'])
    parser.add_argument('--dir_gff', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--seqfile', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--outfile', metavar='PATH', default='gff2genestat.tsv', type=str, help='')
    parser.add_argument('--feature', metavar='STR', default='CDS', type=str, help='')
    parser.add_argument('--multiple_hits', metavar='STR', default='longest', type=str, help='')
    parser.add_argument('--ncpu', metavar='INT', default=1, type=int, help='Number of worker threads.')
    return parser

def read_fasta_seqname(file_path):
    raw_seqnames = []
    if file_path.endswith('.gz'):
        with gzip.open(file_path,'r') as f:
            for line in f:
                line = line.decode("utf-8")
                if ">" in line:
                    raw_seqnames.append(line)
    else:
        with open(file_path,'r') as f:
            for line in f:
                if ">" in line:
                    raw_seqnames.append(line)
    seqnames = pandas.Series(raw_seqnames)
    seqnames = seqnames.str.replace('^>','',regex=True)
    seqnames = seqnames.str.replace('\n$','',regex=True)
    return seqnames

def get_gene_names(seq_names):
    gene_names = seq_names.str.replace('_', 'PLACEHOLDERTEXT', 2)
    gene_names = gene_names.str.replace(r'.*PLACEHOLDERTEXT', '', regex=True)
    return gene_names

@lru_cache(maxsize=8192)
def _get_search_regex(search_ids_tuple):
    regexes = [
        '[^a-zA-Z0-9]' + sid + '(?:$|[^a-zA-Z_0-9])'
        for sid in search_ids_tuple
    ]
    return '|'.join(regexes)


def get_index_bools(gff, search_ids):
    search_ids_tuple = tuple(
        dict.fromkeys(
            sid for sid in search_ids.fillna('').astype(str).tolist() if sid != ''
        )
    )
    if len(search_ids_tuple) == 0:
        return numpy.zeros(gff.shape[0], dtype=bool)
    regex = _get_search_regex(search_ids_tuple)
    conditions = gff.loc[:,'attributes'].str.contains(regex, regex=True, na=False).values
    return conditions

def get_feature_len(gff):
    len_per_line = gff.loc[:,'end'] - gff.loc[:,'start'] + 1
    len_feature = len_per_line.sum()
    return len_feature

def select_id_from_multiple_hits(gff_feat, feat2_ids, multiple_hits):
    if (feat2_ids.shape[0]==1):
        return feat2_ids[0]
    if (multiple_hits=='first'):
        print('First feature was selected: {}'.format(feat2_ids[0]))
        return feat2_ids[0]
    elif (multiple_hits=='longest'):
        longest_id = ''
        longest_size = 0
        for feat2_id in feat2_ids:
            is_feat2 = get_index_bools(gff=gff_feat, search_ids=pandas.Series([feat2_id,]))
            gff_feat2_id = gff_feat.loc[is_feat2,:]
            len_feature = get_feature_len(gff=gff_feat2_id)
            if (len_feature > longest_size):
                longest_size = len_feature
                longest_id = feat2_id
        print('Longest feature was selected: {}'.format(longest_id))
        return longest_id

def extract_by_ids(gff, seq_names, feature, multiple_hits):
    print('Extracting gene IDs: {}'.format(datetime.datetime.now()), flush=True)
    gff_feat = gff.loc[(gff.loc[:,'feature']==feature),:]
    gene_names = get_gene_names(seq_names)
    ub_gene_names = gene_names.str.replace('-','_', regex=False)
    gene_name_variants = pandas.concat([gene_names, ub_gene_names], ignore_index=True)
    search_ids = pandas.concat([seq_names, gene_name_variants], ignore_index=True)
    conditions = get_index_bools(gff=gff_feat, search_ids=search_ids)
    if (conditions.sum()==0): # This block was introduced to parse Ensembl's gff3 file.
        print('Gene IDs were not found in the feature \"{}\". Checking other features.'.format(feature), flush=True)
        search_ids_list = search_ids.tolist()
        for gene_name in gene_name_variants:
            gene_search_ids = pandas.Series([sid for sid in search_ids_list if gene_name in sid])
            for feature2 in ['mRNA','J_gene_segment','V_gene_segment','C_gene_segment','D_gene_segment','gene']:
                gff_feat2 = gff.loc[(gff.loc[:,'feature']==feature2),:]
                conditions = get_index_bools(gff=gff_feat2, search_ids=gene_search_ids)
                if (conditions.sum()==0):
                    print('Gene ID {} was not found in the feature \"{}\".'.format(gene_name, feature2))
                else:
                    feat2_entries = gff_feat2.loc[conditions,:]
                    feat2_attrs = feat2_entries.loc[:,'attributes'].drop_duplicates()
                    txt = 'Found {} gff entries for {} in the feature \"{}\".'
                    print(txt.format(feat2_attrs.shape[0], gene_name, feature2))
                    feat2_ids = feat2_attrs.replace('^ID=','',regex=True).replace(';.*','',regex=True).values
                    feat2_id = select_id_from_multiple_hits(gff_feat, feat2_ids, multiple_hits)
                    feat2_attr = [ attr for attr in feat2_attrs.values if feat2_id in attr ]
                    is_feat2 = get_index_bools(gff=gff_feat, search_ids=pandas.Series([feat2_id,]))
                    gff_feat.loc[is_feat2,'attributes'] = gff_feat.loc[is_feat2,'attributes'] + feat2_attr
                    break
        conditions = get_index_bools(gff=gff_feat, search_ids=search_ids)
    if (conditions.sum()==0):
        print('No match was found.')
    out = gff_feat.loc[conditions,:]
    return out

def add_id_column(gff, seq_names, new_col='gene_id'):
    print('Adding gene id column: {}'.format(datetime.datetime.now()), flush=True)
    gff.loc[:,new_col] = ''
    gene_names = get_gene_names(seq_names)
    ub_gene_names = gene_names.str.replace('-','_',regex=False)
    for i in range(seq_names.shape[0]):
        search_ids = pandas.Series([seq_names.iloc[i], gene_names.iloc[i], ub_gene_names.iloc[i] ])
        conditions = get_index_bools(gff, search_ids)
        if conditions.sum():
            gff.loc[conditions,new_col] = seq_names.iloc[i]
    return gff

def add_intron_info(gff, df_all, id_col='gene_id'):
    print('Adding intron information: {}'.format(datetime.datetime.now()), flush=True)
    seq_names = gff.loc[:,id_col].unique()
    df_intron = pandas.DataFrame(None, index=numpy.arange(len(seq_names)), columns=df_all.columns)
    for j,seq_name in enumerate(seq_names):
        gff_gene = gff.loc[(gff.loc[:,id_col]==seq_name),:].reset_index()
        feature_size = ((gff_gene.loc[:,'end']+1).values - gff_gene.loc[:,'start'].values).sum()
        num_intron = gff_gene.shape[0]-1
        intron_positions = list()
        current_offset = 0
        max_intron_pos = 0
        for i in numpy.arange(gff_gene.shape[0]-1):
            feature_block_size = gff_gene.at[i,'end'] - gff_gene.at[i,'start'] + 1
            intron_pos = current_offset + feature_block_size
            max_intron_pos = intron_pos if (intron_pos > max_intron_pos) else max_intron_pos
            intron_positions.append(str(intron_pos))
            current_offset += feature_block_size
        str_intron_positions = ';'.join(intron_positions)
        df_intron.at[j,'gene_id'] = seq_name
        df_intron.at[j,'feature_size'] = feature_size
        df_intron.at[j,'num_intron'] = num_intron
        df_intron.at[j,'intron_positions'] = str_intron_positions
        if (max_intron_pos > feature_size):
            txt = 'Intron position cannot be greater than feature size: {}, feature_size={}, max intron position = {}'
            raise Exception(txt.format(seq_name, feature_size, max_intron_pos))
    df_all = pandas.concat([df_all, df_intron], axis='index', ignore_index=True)
    return df_all

def add_gene_delim(gff, df_all, id_col='gene_id'):
    print('Adding gene location: {}'.format(datetime.datetime.now()), flush=True)
    seq_names = gff.loc[:,id_col].unique()
    for seq_name in seq_names:
        is_seqname = (df_all.loc[:,id_col]==seq_name)
        gff_gene = gff.loc[(gff.loc[:,id_col]==seq_name),:].reset_index(drop=True)
        df_all.loc[is_seqname,'chromosome'] = gff_gene.at[0,'sequence']
        df_all.loc[is_seqname,'strand'] = gff_gene.at[0,'strand']
        df_all.loc[is_seqname,'start'] = gff_gene.at[0,'start']
        df_all.loc[is_seqname,'end'] = gff_gene.loc[:,'end'].tail(1).values[0]
    return df_all

def filename2sciname(file_name):
    tmp = file_name.replace('_','|',1)
    tmp = re.sub(r'[-_.,].*', '', tmp)
    sci_name = tmp.replace('|', '_')
    return sci_name


def process_single_gff(gff_file, dir_gff, seq_sp_values, feature, multiple_hits, gff_cols, out_cols):
    print("{}: Started processing: {}".format(datetime.datetime.now(), gff_file), flush=True)
    gff_path = os.path.join(dir_gff, gff_file)
    if os.stat(gff_path).st_size == 0:
        sys.stderr.write('Empty file: {}\n'.format(gff_path))
        return pandas.DataFrame(columns=out_cols)
    gff = pandas.read_csv(gff_path, sep='\t', header=None, comment='#', low_memory=False, quoting=3)
    gff.columns = gff_cols
    seq_sp = pandas.Series(seq_sp_values)
    gff_id = extract_by_ids(gff=gff, seq_names=seq_sp, feature=feature, multiple_hits=multiple_hits)
    if gff_id.shape[0] == 0:
        return pandas.DataFrame(columns=out_cols)
    gff_id = add_id_column(gff=gff_id, seq_names=seq_sp)
    df_tmp = pandas.DataFrame(columns=out_cols)
    df_tmp = add_intron_info(gff=gff_id, df_all=df_tmp)
    df_tmp = add_gene_delim(gff=gff_id, df_all=df_tmp)
    return df_tmp

def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    args.ncpu = max(1, int(args.ncpu))
    start_time = time.time()
    print('gff2genestat.py started: {}'.format(datetime.datetime.now()))

    out_cols = ['gene_id', 'feature_size', 'num_intron', 'intron_positions', 'chromosome', 'start', 'end','strand']
    gff_cols = ['sequence','source','feature','start','end','score','strand','phase','attributes']
    seq_names = read_fasta_seqname(file_path=args.seqfile)
    gff_files = [gff for gff in os.listdir(args.dir_gff) if not gff.startswith('.')]
    gff_files = [gff for gff in gff_files if gff.endswith(('.gff', '.gtf', '.gff3', '.gff.gz', '.gtf.gz', '.gff3.gz'))]
    tasks = []
    for gff_file in gff_files:
        sp_ub = filename2sciname(file_name=gff_file)
        seq_sp = seq_names[seq_names.str.startswith(sp_ub)]
        if seq_sp.shape[0] == 0:
            continue
        tasks.append((gff_file, seq_sp.tolist()))

    frames = []
    if args.ncpu == 1 or len(tasks) <= 1:
        for gff_file, seq_sp_values in tasks:
            df_tmp = process_single_gff(
                gff_file=gff_file,
                dir_gff=args.dir_gff,
                seq_sp_values=seq_sp_values,
                feature=args.feature,
                multiple_hits=args.multiple_hits,
                gff_cols=gff_cols,
                out_cols=out_cols,
            )
            if df_tmp.shape[0] > 0:
                frames.append(df_tmp)
    else:
        max_workers = min(args.ncpu, len(tasks))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    process_single_gff,
                    gff_file,
                    args.dir_gff,
                    seq_sp_values,
                    args.feature,
                    args.multiple_hits,
                    gff_cols,
                    out_cols,
                ): (idx, gff_file)
                for idx, (gff_file, seq_sp_values) in enumerate(tasks)
            }
            ordered_results = {}
            for future in as_completed(futures):
                idx, gff_file = futures[future]
                try:
                    df_tmp = future.result()
                    ordered_results[idx] = df_tmp
                except Exception as exc:
                    raise RuntimeError('Failed processing {}: {}'.format(gff_file, exc)) from exc
            for idx in sorted(ordered_results.keys()):
                df_tmp = ordered_results[idx]
                if df_tmp.shape[0] > 0:
                    frames.append(df_tmp)

    if len(frames) > 0:
        df_all = pandas.concat(frames, axis='index', ignore_index=True)
    else:
        df_all = pandas.DataFrame(columns=out_cols)

    df_all = df_all.drop_duplicates()
    num_input = len(seq_names)
    num_output = df_all.shape[0]
    print('Number of input genes: {}'.format(num_input), flush=True)
    print('Number of output entries: {}'.format(num_output), flush=True)
    if (num_input != num_output):
        df_all['gene_id'] = df_all['gene_id'].fillna('').astype(str)
        for seq_name in seq_names.sort_values():
            gene_name = get_gene_names(pandas.Series([seq_name])).iloc[0]
            is_found = df_all.loc[:, 'gene_id'].str.contains(gene_name).any()
            if is_found:
                continue
            print('Missing in the output file: {}'.format(seq_name), flush=True)
    df_all.to_csv(args.outfile, sep='\t', header=True, index=False)
    print('gff2genestat.py completed in {:,} secs: {}'.format(int(time.time() - start_time), datetime.datetime.now()))


if __name__ == '__main__':
    main()
