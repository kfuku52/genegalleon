#!/usr/bin/env python3

import argparse
import datetime
import gzip
import os
import re
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
from Bio import SeqIO
from ete4 import NCBITaxa

pandas.options.mode.chained_assignment = None
DOMAIN_RANK_ALIASES = frozenset({'domain', 'superkingdom'})
MMSEQS_TAXONOMY_COLUMNS = ['query', 'lca_taxid', 'lca_sciname', 'lineage_taxids']


def rank_aliases(rank):
    rank = str(rank).strip().lower()
    if rank in DOMAIN_RANK_ALIASES:
        return DOMAIN_RANK_ALIASES
    return frozenset({rank})


def rank_matches(rank_name, requested_rank):
    return str(rank_name).strip().lower() in rank_aliases(requested_rank)


def resolve_lineage_rank_index(df_lineage, requested_rank):
    rank_mask = df_lineage['rank'].astype(str).map(lambda rank_name: rank_matches(rank_name, requested_rank))
    matching_indexes = df_lineage.index[rank_mask]
    if len(matching_indexes) == 0:
        available_ranks = ', '.join(df_lineage['rank'].astype(str).tolist())
        raise ValueError(
            'Requested rank {!r} was not found in the species lineage. Available ranks: {}'.format(
                requested_rank, available_ranks
            )
        )
    return int(matching_indexes[0])


def ncbi_taxa(dbfile=None):
    if dbfile:
        return NCBITaxa(dbfile=dbfile)
    return NCBITaxa()


def parse_lineage_taxids(lineage_taxid_str):
    lineage_taxid_str = str(lineage_taxid_str)
    if lineage_taxid_str in ('', 'nan'):
        return []
    return [int(x) for x in lineage_taxid_str.split(';') if x not in ('', 'nan')]


def build_taxid_rank_lookup(ncbi, lineage_taxid_strs):
    unique_taxids = sorted({taxid for lineage_taxid_str in lineage_taxid_strs for taxid in parse_lineage_taxids(lineage_taxid_str)})
    if len(unique_taxids) == 0:
        return {}
    rank_map = ncbi.get_rank(unique_taxids)
    return {int(taxid): str(rank_name) for taxid, rank_name in rank_map.items()}


def resolve_rank_taxid_from_lineage(ncbi, rank, lineage_taxid_str, taxid_to_rank=None):
    lineage_taxids = parse_lineage_taxids(lineage_taxid_str)
    if len(lineage_taxids) == 0:
        return -1

    if taxid_to_rank is None:
        lineage_ranks = ncbi.get_rank(lineage_taxids)
    else:
        lineage_ranks = taxid_to_rank
    for taxid in lineage_taxids:
        rank_name = lineage_ranks.get(taxid, '')
        if rank_matches(rank_name, rank):
            return int(taxid)
    return 0


def main():
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_file', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--mmseqs2taxonomy_tsv', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--fx2tab_tsv', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--species_name', metavar='GENUS_SPECIES', default='', type=str, help='')
    parser.add_argument(
        '--rank',
        metavar='domain|superkingdom|kingdom|phylum|subphylum|class|order|family|genus|species',
        default='class',
        type=str,
        help='Stringency for taxnomic rank incompatibility check. If "genus", different genera are considered incompatible.',
    )
    parser.add_argument('--rename_seq', metavar='yes|no', default='yes', type=strtobool, help='')
    parser.add_argument('--rename_prefix', metavar='STR', default='scaffold', type=str, help='')
    parser.add_argument('--verbose', metavar='yes|no', default='no', type=strtobool, help='')
    parser.add_argument('--ncpu', metavar='INT', default=1, type=int, help='Number of CPU processes.')
    args = parser.parse_args()
    args.ncpu = max(1, int(args.ncpu))
    print('{} started: {}'.format(sys.argv[0], datetime.datetime.now()))

    df = pandas.read_csv(
        args.mmseqs2taxonomy_tsv,
        header=None,
        sep='\t',
        low_memory=False,
        usecols=[0, 1, 3, 8],
        names=MMSEQS_TAXONOMY_COLUMNS,
        dtype={'query': str, 'lca_sciname': str, 'lineage_taxids': str},
    )  # https://github.com/soedinglab/MMseqs2/wiki#taxonomy-format
    df['lca_taxid'] = pandas.to_numeric(df['lca_taxid'], errors='coerce').fillna(0).astype('int64', copy=False)

    taxonomy_dbfile = os.environ.get('GG_TAXONOMY_DBFILE', '').strip()
    if taxonomy_dbfile:
        print('Using NCBI taxonomy DB file: {}'.format(taxonomy_dbfile))
        ncbi = ncbi_taxa(dbfile=taxonomy_dbfile)
    else:
        ncbi = ncbi_taxa()

    fx = pandas.read_csv(
        args.fx2tab_tsv,
        header=0,
        sep='\t',
        low_memory=False,
        usecols=['#id', 'length', 'GC'],
        dtype={'#id': str, 'length': 'int64', 'GC': 'float64'},
    ).rename(columns={'#id': 'query'})
    df = pandas.merge(df, fx, on='query')
    df.sort_values(by='length', ascending=False, inplace=True, ignore_index=True)

    lineage_series = df['lineage_taxids'].astype(str)
    lineage_tasks = pandas.unique(lineage_series).tolist()
    print('Processing {} unique lineage patterns for rank conversion: {}'.format(len(lineage_tasks), args.rank))

    taxid_to_rank = build_taxid_rank_lookup(ncbi, lineage_tasks)
    lineage_to_rank_taxid = {}
    for lineage_taxid_str in lineage_tasks:
        lineage_to_rank_taxid[lineage_taxid_str] = resolve_rank_taxid_from_lineage(
            ncbi,
            args.rank,
            lineage_taxid_str,
            taxid_to_rank=taxid_to_rank,
        )

    rank_taxid = lineage_series.map(lineage_to_rank_taxid).fillna(-1).to_numpy(dtype=int, copy=False)
    lca_taxid = df['lca_taxid'].to_numpy(dtype=int, copy=False)
    df['aligned_taxid'] = numpy.where(rank_taxid == -1, 0, numpy.where(rank_taxid == 0, lca_taxid, rank_taxid))

    sci_name = args.species_name.replace('_', ' ')
    print('Scientific name for the taxid search: {}'.format(sci_name))
    name2taxid = ncbi.get_name_translator([sci_name])
    if len(name2taxid) == 0:
        print('Species name not found in the NCBI Taxonomy database: {}'.format(sci_name))
        print('Trying to retrieve genus-level taxonomic information.')
        sci_name = re.sub(' .*', '', sci_name)
        name2taxid = ncbi.get_name_translator([sci_name])
    if len(name2taxid) == 0:
        print('Exiting. Genus name not found in the NCBI Taxonomy database: {}'.format(sci_name))
        sys.exit(1)

    lineage_taxids = ncbi.get_lineage(name2taxid[sci_name][0])
    lineage_ranks = ncbi.get_rank(lineage_taxids)
    lineage_names = ncbi.get_taxid_translator(lineage_taxids)
    df_lineage = pandas.DataFrame({'rank': lineage_ranks, 'name': lineage_names})
    df_lineage['taxid'] = df_lineage.index
    df_lineage = df_lineage.loc[lineage_taxids, :].reset_index(drop=True)
    tmp = [lr + ':' + ln for lr, ln in zip(lineage_ranks.values(), lineage_names.values())]
    print('NCBI Taxonomy lineage of {}: {}'.format(sci_name, ', '.join(tmp)))

    try:
        rank_index = resolve_lineage_rank_index(df_lineage, args.rank)
    except ValueError as exc:
        print('Exiting. {}'.format(exc))
        sys.exit(1)
    selected_lineage_taxids = df_lineage.loc[:rank_index, 'taxid'].tolist()
    selected_lineage_names = df_lineage.loc[:rank_index, 'name'].tolist()
    print('Lineages for taxonomic consistency check: {}'.format(', '.join(selected_lineage_names)))

    is_lineage_seq = df['aligned_taxid'].isin(selected_lineage_taxids)
    is_unclassified = (df['aligned_taxid'] == 0)
    is_compatible_seq = (is_lineage_seq | is_unclassified)
    df['is_compatible_lineage'] = is_compatible_seq

    num_lineage_seq = is_lineage_seq.sum()
    num_nonlineage_seq = df.shape[0] - num_lineage_seq
    num_unclassified = is_unclassified.sum()
    bp_lineage_seq = df.loc[is_lineage_seq, 'length'].sum()
    bp_compatible_seq = df.loc[is_compatible_seq, 'length'].sum()
    bp_nonlineage_seq = df.loc[(~is_lineage_seq) & (~is_unclassified), 'length'].sum()
    bp_unclassified = df.loc[is_unclassified, 'length'].sum()
    gc_lineage_seq = df.loc[is_lineage_seq, 'GC'].mean()
    gc_nonlineage_seq = df.loc[(~is_lineage_seq) & (~is_unclassified), 'GC'].mean()
    gc_unclassified = df.loc[is_unclassified, 'GC'].mean()

    print(
        'Seqeunces compatible with the lineage of {}: {:,} sequences totaling {:,} bp with the mean GC of {:.1f}%'.format(
            sci_name, num_lineage_seq, bp_lineage_seq, gc_lineage_seq
        )
    )
    print(
        'Seqeunces non-compatible with the lineage of {}: {:,} sequences totaling {:,} bp with the mean GC of {:.1f}%'.format(
            sci_name, num_nonlineage_seq, bp_nonlineage_seq, gc_nonlineage_seq
        )
    )
    print(
        'Unclassified sequences (taxid=0): {:,} sequenced totaling {:,} bp with the mean GC of {:.1f}%. These sequences will be retained.'.format(
            num_unclassified, bp_unclassified, gc_unclassified
        )
    )

    nonlineage_summary = (
        df.loc[~is_compatible_seq, ['lca_taxid', 'lca_sciname', 'length', 'GC']]
        .groupby(['lca_taxid', 'lca_sciname'], sort=False)
        .agg(num_seq=('lca_taxid', 'size'), bp_seq=('length', 'sum'), gc_seq=('GC', 'mean'))
        .reset_index()
    )
    for row in nonlineage_summary.itertuples(index=False):
        txt = '{:,} sequences totaling {:,} bp with the mean GC of {:.1f}% were identified as being from {} (taxid={}). These sequences will be removed.'
        print(txt.format(row.num_seq, row.bp_seq, row.gc_seq, row.lca_sciname, row.lca_taxid))

    sorted_compatible_seqnames = df.loc[is_compatible_seq, 'query'].tolist()
    if args.rename_seq:
        print('Sequences will be renamed. The longest sequence name is {}1.'.format(args.rename_prefix))
        new_seq_names = [args.rename_prefix + str(i + 1) for i in range(len(sorted_compatible_seqnames))]
        df['new_seq_name'] = ''
        df.loc[is_compatible_seq, 'new_seq_name'] = new_seq_names
        query_to_length = dict(zip(df['query'], df['length']))
    df.to_csv('lineage_compatibility.tsv', sep='\t', index=False)

    compatible_id_set = set(sorted_compatible_seqnames)
    compatible_records = {}
    if args.fasta_file.endswith('gz') or args.fasta_file.endswith('gzipped'):
        handler = SeqIO.parse(gzip.open(args.fasta_file, 'rt'), 'fasta')
    else:
        handler = SeqIO.parse(args.fasta_file, 'fasta')

    for record in handler:
        if record.id in compatible_id_set:
            compatible_records[record.id] = record

    newseq_index = 0
    print('Summary:')
    print('Input: {:,} sequences totaling {:,} bp'.format(df.shape[0], df['length'].sum()))
    print(
        'Output: {:,} sequences totaling {:,} bp'.format(
            len(sorted_compatible_seqnames),
            bp_compatible_seq,
        )
    )
    with open('clean_sequences.fa', 'w') as output_handle:
        for scs in sorted_compatible_seqnames:
            record = compatible_records[scs]
            record.name = ''
            record.description = ''
            if args.rename_seq:
                old_seq_name = record.id
                new_seq_name = new_seq_names[newseq_index]
                seq_bp = query_to_length[old_seq_name]
                record.id = new_seq_name
                if args.verbose:
                    print('Renaming: {} -> {} ({:,} bp)'.format(old_seq_name, new_seq_name, seq_bp))
                newseq_index += 1
            SeqIO.write(record, output_handle, 'fasta')

    print('{} ended: {} ({:,} sec)'.format(sys.argv[0], datetime.datetime.now(), int(time.time() - start_time)))


if __name__ == '__main__':
    main()
