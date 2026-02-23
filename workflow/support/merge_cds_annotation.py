#!/usr/bin/env python3
# coding: utf-8

import argparse
import datetime
import gzip
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas


def read_fasta_ids(path):
    ids = []
    opener = gzip.open if path.endswith('.gz') else open
    mode = 'rt' if path.endswith('.gz') else 'r'
    with opener(path, mode) as handle:
        for line in handle:
            if line.startswith('>'):
                ids.append(line[1:].strip().split()[0])
    return ids


def load_orthogroup_map(path, species_col):
    if not os.path.exists(path):
        print('File not found: {}'.format(path), flush=True)
        return None
    print('{}: Processing: {}'.format(datetime.datetime.now(), path), flush=True)

    header = pandas.read_csv(path, sep='\t', nrows=0)
    if ('Orthogroup' not in header.columns) or (species_col not in header.columns):
        print('Required columns were not found in orthogroup file: Orthogroup, {}'.format(species_col), flush=True)
        return None

    tmp = pandas.read_csv(path, sep='\t', usecols=['Orthogroup', species_col], low_memory=False)
    tmp = tmp.dropna(subset=[species_col]).copy()
    if tmp.empty:
        return None

    tmp['gene_id'] = tmp[species_col].astype(str).str.split(',')
    tmp = tmp.explode('gene_id')
    tmp['gene_id'] = tmp['gene_id'].str.strip()
    tmp = tmp.loc[tmp['gene_id'] != '', ['gene_id', 'Orthogroup']]
    tmp = tmp.drop_duplicates(subset=['gene_id'], keep='first')
    return tmp.set_index('gene_id')['Orthogroup']


def load_uniprot(path):
    if not os.path.exists(path):
        print('File not found: {}'.format(path), flush=True)
        return None
    print('{}: Processing: {}'.format(datetime.datetime.now(), path), flush=True)
    return pandas.read_csv(path, sep='\t', header=0, index_col=None, low_memory=False)


def load_busco(path):
    if not os.path.exists(path):
        print('File not found: {}'.format(path), flush=True)
        return None
    print('{}: Processing: {}'.format(datetime.datetime.now(), path), flush=True)
    colnames = ['busco_id', 'busco_status', 'busco_sequence', 'busco_score', 'busco_length', 'busco_orthodb_url', 'busco_description']
    tmp = pandas.read_csv(path, sep='\t', header=None, index_col=None, comment='#', low_memory=False, names=colnames)
    tmp['gene_id'] = tmp['busco_sequence'].astype(str).str.replace(':.*', '', regex=True)
    tmp = tmp.loc[tmp['gene_id'].notna() & (tmp['gene_id'] != ''), :].copy()
    if tmp.empty:
        return None
    agg_cols = [c for c in tmp.columns if c != 'gene_id']
    tmp_agg = tmp.loc[:, ['gene_id'] + agg_cols].astype(str)
    tmp_agg = tmp_agg.groupby('gene_id', as_index=False, sort=False)[agg_cols].agg('; '.join)
    return tmp_agg


def load_fx2tab(path):
    if not os.path.exists(path):
        print('File not found: {}'.format(path), flush=True)
        return None
    print('{}: Processing: {}'.format(datetime.datetime.now(), path), flush=True)
    tmp = pandas.read_csv(path, sep='\t', header=0, index_col=None, low_memory=False)
    tmp.columns = tmp.columns.str.replace('length', 'cds_length')
    tmp.columns = tmp.columns.str.replace('#id', 'gene_id')
    return tmp


def load_gff_info(path):
    if not os.path.exists(path):
        print('File not found: {}'.format(path), flush=True)
        return None
    print('{}: Processing: {}'.format(datetime.datetime.now(), path), flush=True)
    return pandas.read_csv(path, sep='\t', header=0, index_col=None, low_memory=False)


def load_expression(path):
    if not os.path.exists(path):
        print('File not found: {}'.format(path), flush=True)
        return None
    print('{}: Processing: {}'.format(datetime.datetime.now(), path), flush=True)
    tmp = pandas.read_csv(path, sep='\t', header=0, index_col=None, low_memory=False)
    if 'Unnamed: 0' in tmp.columns:
        tmp.columns = tmp.columns.str.replace('Unnamed: 0', 'gene_id')
    else:
        tmp.columns = tmp.columns.str.replace(tmp.columns[0], 'gene_id')
    return tmp


def load_mmseqs(path):
    if not os.path.exists(path):
        print('File not found: {}'.format(path), flush=True)
        return None
    print('{}: Processing: {}'.format(datetime.datetime.now(), path), flush=True)
    tmp = pandas.read_csv(path, sep='\t', header=None, index_col=None, low_memory=False)
    tmp.columns = [
        'gene_id', 'lca_taxid', 'lca_rank', 'lca_sciname',
        'num_assigned_protein_fragment', 'num_labeled_protein_fragment',
        'num_taxid_supporting_protein_fragment', 'fraction_evalue_support',
        'lineage_taxids',
    ]  # https://github.com/soedinglab/MMseqs2/wiki#taxonomy-format
    return tmp


def join_if_available(df, tmp):
    if tmp is None:
        return df
    if ('gene_id' not in tmp.columns):
        return df
    tmp_indexed = tmp.drop_duplicates(subset=['gene_id'], keep='first').set_index('gene_id')
    return df.join(tmp_indexed, how='left')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scientific_name', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--cds_fasta', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--orthogroup_tsv', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--uniprot_tsv', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--busco_tsv', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--expression_tsv', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--gff_info', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--mmseqs2taxonomy_tsv', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--fx2tab', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--out_tsv', metavar='PATH', default='cds_annotation.tsv', type=str, help='')
    parser.add_argument('--ncpu', metavar='INT', default=1, type=int, help='Number of worker threads.')
    args = parser.parse_args()
    args.ncpu = max(1, int(args.ncpu))

    seqnames = read_fasta_ids(args.cds_fasta)
    df = pandas.DataFrame(index=pandas.Index(seqnames, name='gene_id'))

    orthogroup_map = load_orthogroup_map(args.orthogroup_tsv, args.scientific_name)
    if orthogroup_map is not None:
        df.loc[:, 'orthogroup'] = df.index.to_series().map(orthogroup_map).fillna('')

    tasks = {
        'uniprot': (load_uniprot, args.uniprot_tsv),
        'busco': (load_busco, args.busco_tsv),
        'fx2tab': (load_fx2tab, args.fx2tab),
        'gff_info': (load_gff_info, args.gff_info),
        'expression': (load_expression, args.expression_tsv),
        'mmseqs': (load_mmseqs, args.mmseqs2taxonomy_tsv),
    }

    loaded = {k: None for k in tasks}
    if args.ncpu == 1:
        for key, (func, path) in tasks.items():
            loaded[key] = func(path)
    else:
        max_workers = min(args.ncpu, len(tasks))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(func, path): key for key, (func, path) in tasks.items()}
            for future in as_completed(futures):
                key = futures[future]
                loaded[key] = future.result()

    for key in ['uniprot', 'busco', 'fx2tab', 'gff_info', 'expression', 'mmseqs']:
        df = join_if_available(df, loaded[key])

    df = df.reset_index()
    df.to_csv(args.out_tsv, sep='\t', index=False)
    print('merge_cds_annotation.py completed!', flush=True)


if __name__ == '__main__':
    main()
