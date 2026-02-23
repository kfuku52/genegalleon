#!/usr/bin/env python3

import argparse
import gzip
import os

import pandas

OUTPUT_COLUMNS = [
    'gene_id',
    'sprot_best',
    'sprot_alias',
    'sprot_coverage',
    'sprot_identity',
    'sprot_evalue',
    'sprot_recname',
    'signal_start',
    'signal_end',
    'transmem_aa',
    'transmem_count',
    'transmem_regions',
    'kegg_gene',
    'kegg_orthology',
    'go_ids',
    'go_aspects',
    'go_terms',
    'go_evidence',
    'gene_name_primary',
    'gene_name_synonyms',
    'ec_numbers',
    'subcellular_location',
    'keywords',
    'interpro_ids',
    'pfam_ids',
    'reactome_ids',
    'organism',
    'taxid',
]

METADATA_COLUMNS = [
    'signal_start',
    'signal_end',
    'transmem_aa',
    'transmem_count',
    'transmem_regions',
    'kegg_gene',
    'kegg_orthology',
    'go_ids',
    'go_aspects',
    'go_terms',
    'go_evidence',
    'gene_name_primary',
    'gene_name_synonyms',
    'ec_numbers',
    'subcellular_location',
    'keywords',
    'interpro_ids',
    'pfam_ids',
    'reactome_ids',
    'organism',
    'taxid',
]

DIAMOND_COLUMNS = [
    'qseqid',
    'sseqid',
    'pident',
    'length',
    'evalue',
    'bitscore',
    'qlen',
]


def open_text(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def read_fasta_ids(path):
    ids = []
    if (path == '') or (not os.path.exists(path)):
        return ids
    with open_text(path) as fin:
        for line in fin:
            if line.startswith('>'):
                ids.append(line[1:].strip().split()[0])
    return ids


def read_uniprot_descriptions(path):
    desc = {}
    if (path == '') or (not os.path.exists(path)):
        return desc
    with open_text(path) as fin:
        for line in fin:
            if line.startswith('>'):
                header = line[1:].strip()
                if header == '':
                    continue
                parts = header.split(None, 1)
                acc = parts[0]
                text = parts[1] if len(parts) > 1 else ''
                desc[acc] = text
    return desc


def read_uniprot_metadata(path):
    if (path == '') or (not os.path.exists(path)):
        return pandas.DataFrame(columns=['accession'] + METADATA_COLUMNS)
    df = pandas.read_csv(path, sep='\t', dtype=str, low_memory=False, keep_default_na=False)
    if 'accession' not in df.columns:
        if 'sprot_best' in df.columns:
            df = df.rename(columns={'sprot_best': 'accession'})
        else:
            return pandas.DataFrame(columns=['accession'] + METADATA_COLUMNS)
    keep_cols = ['accession']
    for col in METADATA_COLUMNS:
        if col not in df.columns:
            df[col] = ''
        keep_cols.append(col)
    df = df.loc[:, keep_cols]
    for col in keep_cols:
        df[col] = df[col].fillna('')
    df = df.drop_duplicates(subset=['accession'], keep='first')
    return df


def fmt_float(value, precision=2):
    if pandas.isna(value):
        return ''
    txt = f"{float(value):.{precision}f}"
    return txt.rstrip('0').rstrip('.')


def fmt_evalue(value):
    if pandas.isna(value):
        return ''
    return f"{float(value):.3e}"


def init_output_frame(gene_ids):
    out = pandas.DataFrame({'gene_id': gene_ids})
    for col in OUTPUT_COLUMNS:
        if col != 'gene_id':
            out[col] = ''
    return out


def load_best_diamond_hits(path):
    if (path == '') or (not os.path.exists(path)) or (os.path.getsize(path) == 0):
        return pandas.DataFrame(columns=DIAMOND_COLUMNS + ['coverage'])

    df = pandas.read_csv(path, sep='\t', header=None, names=DIAMOND_COLUMNS, dtype=str, low_memory=False)
    if df.shape[0] == 0:
        return pandas.DataFrame(columns=DIAMOND_COLUMNS + ['coverage'])

    for col in ['pident', 'length', 'evalue', 'bitscore', 'qlen']:
        df[col] = pandas.to_numeric(df[col], errors='coerce')

    df = df.sort_values(
        by=['qseqid', 'bitscore', 'evalue', 'pident', 'length'],
        ascending=[True, False, True, False, False],
        na_position='last',
    )
    df = df.drop_duplicates(subset=['qseqid'], keep='first')
    df['coverage'] = (df['length'] / df['qlen']) * 100.0

    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--diamond_tsv', metavar='PATH', default='', type=str)
    parser.add_argument('--query_fasta', metavar='PATH', default='', type=str)
    parser.add_argument('--uniprot_fasta', metavar='PATH', default='', type=str)
    parser.add_argument('--uniprot_meta_tsv', metavar='PATH', default='', type=str)
    parser.add_argument('--outfile', metavar='PATH', required=True, type=str)
    args = parser.parse_args()

    query_ids = read_fasta_ids(args.query_fasta)
    hits = load_best_diamond_hits(args.diamond_tsv)

    if len(query_ids) == 0:
        query_ids = hits['qseqid'].dropna().drop_duplicates().tolist()

    out = init_output_frame(query_ids)

    if hits.shape[0] > 0:
        tmp = hits.loc[:, ['qseqid', 'sseqid', 'coverage', 'pident', 'evalue']].copy()
        tmp = tmp.rename(columns={'qseqid': 'gene_id', 'sseqid': 'sprot_best'})

        tmp['sprot_alias'] = tmp['sprot_best']
        tmp['sprot_coverage'] = tmp['coverage'].map(lambda x: fmt_float(x, precision=2))
        tmp['sprot_identity'] = tmp['pident'].map(lambda x: fmt_float(x, precision=2))
        tmp['sprot_evalue'] = tmp['evalue'].map(fmt_evalue)

        tmp = tmp.loc[:, ['gene_id', 'sprot_best', 'sprot_alias', 'sprot_coverage', 'sprot_identity', 'sprot_evalue']]
        out = pandas.merge(out, tmp, how='left', on='gene_id', suffixes=('', '_new'))

        for col in ['sprot_best', 'sprot_alias', 'sprot_coverage', 'sprot_identity', 'sprot_evalue']:
            new_col = f'{col}_new'
            if new_col in out.columns:
                out[col] = out[new_col].fillna(out[col])
                out = out.drop(columns=[new_col])

    uniprot_desc = read_uniprot_descriptions(args.uniprot_fasta)
    if len(uniprot_desc) > 0:
        out['sprot_recname'] = out['sprot_best'].map(uniprot_desc).fillna('')

    metadata = read_uniprot_metadata(args.uniprot_meta_tsv)
    if metadata.shape[0] > 0:
        tmp_meta = metadata.rename(columns={'accession': 'sprot_best'})
        out = pandas.merge(out, tmp_meta, how='left', on='sprot_best', suffixes=('', '_meta'))
        for col in METADATA_COLUMNS:
            new_col = f'{col}_meta'
            if new_col in out.columns:
                out[col] = out[new_col].fillna(out[col])
                out = out.drop(columns=[new_col])

    for col in OUTPUT_COLUMNS:
        if col not in out.columns:
            out[col] = ''
        out[col] = out[col].fillna('')

    out = out.loc[:, OUTPUT_COLUMNS]
    out.to_csv(args.outfile, sep='\t', index=False)


if __name__ == '__main__':
    main()
