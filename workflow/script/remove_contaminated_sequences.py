#! /usr/bin/env python

import argparse
import datetime
import gzip
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from distutils.util import strtobool

import numpy
import pandas
from Bio import SeqIO
from ete4 import NCBITaxa

pandas.options.mode.chained_assignment = None

_WORKER_NCBI = None
_WORKER_RANK = None


def ncbi_taxa(dbfile=None):
    if dbfile:
        return NCBITaxa(dbfile=dbfile)
    return NCBITaxa()


def _init_worker(rank, taxonomy_dbfile):
    global _WORKER_NCBI, _WORKER_RANK
    _WORKER_RANK = rank
    if taxonomy_dbfile:
        _WORKER_NCBI = ncbi_taxa(dbfile=taxonomy_dbfile)
    else:
        _WORKER_NCBI = ncbi_taxa()


def resolve_rank_taxid_from_lineage(ncbi, rank, lineage_taxid_str):
    lineage_taxid_str = str(lineage_taxid_str)
    if lineage_taxid_str in ('', 'nan'):
        return -1

    lineage_taxids = [int(x) for x in lineage_taxid_str.split(';') if x not in ('', 'nan')]
    if len(lineage_taxids) == 0:
        return -1

    lineage_ranks = ncbi.get_rank(lineage_taxids)
    for taxid, rank_name in lineage_ranks.items():
        if rank_name == rank:
            return int(taxid)
    return 0


def _resolve_rank_taxid_worker(lineage_taxid_str):
    rank_taxid = resolve_rank_taxid_from_lineage(_WORKER_NCBI, _WORKER_RANK, lineage_taxid_str)
    return lineage_taxid_str, rank_taxid


def main():
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_file', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--mmseqs2taxonomy_tsv', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--fx2tab_tsv', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--species_name', metavar='GENUS_SPECIES', default='', type=str, help='')
    parser.add_argument(
        '--rank',
        metavar='superkingdom|kingdom|phylum|subphylum|class|order|family|genus|species',
        default='class',
        type=str,
        help='Stringency for taxnomic rank incompatibility check. If "genus", different genera are considered incompatible.',
    )
    parser.add_argument('--rename_seq', metavar='yes|no', default='yes', type=strtobool, help='')
    parser.add_argument('--rename_prefix', metavar='STR', default='scaffold', type=str, help='')
    parser.add_argument('--verbose', metavar='yes|no', default='no', type=strtobool, help='')
    parser.add_argument('--ncpu', metavar='INT', default=1, type=int, help='Number of CPU processes.')
    # Backward-compatible alias.
    parser.add_argument('--threads', dest='ncpu', metavar='INT', type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.ncpu = max(1, int(args.ncpu))
    print('{} started: {}'.format(sys.argv[0], datetime.datetime.now()))

    df = pandas.read_csv(args.mmseqs2taxonomy_tsv, header=None, sep='\t', low_memory=False)
    df.columns = [
        'query', 'lca_taxid', 'lca_rank', 'lca_sciname',
        'num_assigned_protein_fragment', 'num_labeled_protein_fragment',
        'num_taxid_supporting_protein_fragment', 'fraction_evalue_support',
        'lineage_taxids',
    ]  # https://github.com/soedinglab/MMseqs2/wiki#taxonomy-format
    df['lca_taxid'] = pandas.to_numeric(df['lca_taxid'], errors='coerce').fillna(0).astype(int)

    taxonomy_dbfile = os.environ.get('GG_TAXONOMY_DBFILE', '').strip()
    if taxonomy_dbfile:
        print('Using NCBI taxonomy DB file: {}'.format(taxonomy_dbfile))
        ncbi = ncbi_taxa(dbfile=taxonomy_dbfile)
    else:
        ncbi = ncbi_taxa()

    fx = pandas.read_csv(args.fx2tab_tsv, header=0, sep='\t', low_memory=False)
    fx.columns = fx.columns.str.replace('#id', 'query')
    df = pandas.merge(df, fx, on='query')
    df = df.sort_values(by='length', ascending=False).reset_index(drop=True)

    unique_tax = df.loc[:, ['lca_taxid', 'lineage_taxids']].drop_duplicates().reset_index(drop=True)
    lineage_tasks = unique_tax['lineage_taxids'].astype(str).drop_duplicates().tolist()
    print('Processing {} unique lineage patterns for rank conversion: {}'.format(len(lineage_tasks), args.rank))

    lineage_to_rank_taxid = {}
    if args.ncpu > 1 and len(lineage_tasks) > 1:
        with ProcessPoolExecutor(
            max_workers=args.ncpu,
            initializer=_init_worker,
            initargs=(args.rank, taxonomy_dbfile),
        ) as executor:
            for lineage_taxid_str, rank_taxid in executor.map(_resolve_rank_taxid_worker, lineage_tasks):
                lineage_to_rank_taxid[str(lineage_taxid_str)] = int(rank_taxid)
    else:
        for lineage_taxid_str in lineage_tasks:
            lineage_to_rank_taxid[lineage_taxid_str] = resolve_rank_taxid_from_lineage(ncbi, args.rank, lineage_taxid_str)

    lineage_series = df['lineage_taxids'].astype(str)
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

    rank_index = df_lineage.index[df_lineage['rank'] == args.rank][0]
    selected_lineage_taxids = df_lineage.loc[:rank_index, 'taxid'].tolist()
    selected_lineage_names = df_lineage.loc[:rank_index, 'name'].tolist()
    print('Lineages for taxonomic consistency check: {}'.format(', '.join(selected_lineage_names)))

    df['is_compatible_lineage'] = False
    is_lineage_seq = df['aligned_taxid'].isin(selected_lineage_taxids)
    is_unclassified = (df['aligned_taxid'] == 0)
    df.loc[is_lineage_seq, 'is_compatible_lineage'] = True
    df.loc[is_unclassified, 'is_compatible_lineage'] = True

    num_lineage_seq = is_lineage_seq.sum()
    num_nonlineage_seq = df.shape[0] - num_lineage_seq
    num_unclassified = is_unclassified.sum()
    bp_lineage_seq = df.loc[is_lineage_seq, 'length'].sum()
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

    nonlineage_taxids = df.loc[~df['is_compatible_lineage'], 'lca_taxid'].unique()
    for nonlineage_taxid in nonlineage_taxids:
        is_nonlineage = (df['lca_taxid'] == nonlineage_taxid)
        num_seq = is_nonlineage.sum()
        bp_seq = df.loc[is_nonlineage, 'length'].sum()
        gc_seq = df.loc[is_nonlineage, 'GC'].mean()
        nonlineage_sciname = df.loc[is_nonlineage, 'lca_sciname'].values[0]
        txt = '{:,} sequences totaling {:,} bp with the mean GC of {:.1f}% were identified as being from {} (taxid={}). These sequences will be removed.'
        print(txt.format(num_seq, bp_seq, gc_seq, nonlineage_sciname, nonlineage_taxid))

    sorted_compatible_seqnames = df.loc[(is_lineage_seq | is_unclassified), 'query'].tolist()
    if args.rename_seq:
        print('Sequences will be renamed. The longest sequence name is {}1.'.format(args.rename_prefix))
        new_seq_names = [args.rename_prefix + str(i + 1) for i in range(len(sorted_compatible_seqnames))]
        df['new_seq_name'] = ''
        df.loc[(is_lineage_seq | is_unclassified), 'new_seq_name'] = new_seq_names
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

    sorted_compatible_record_list = []
    newseq_index = 0
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
        sorted_compatible_record_list.append(record)

    print('Summary:')
    print('Input: {:,} sequences totaling {:,} bp'.format(df.shape[0], df['length'].sum()))
    print(
        'Output: {:,} sequences totaling {:,} bp'.format(
            len(sorted_compatible_record_list),
            df.loc[(df['is_compatible_lineage']), 'length'].sum(),
        )
    )
    with open('clean_sequences.fa', 'w') as output_handle:
        SeqIO.write(sorted_compatible_record_list, output_handle, 'fasta')

    print('{} ended: {} ({:,} sec)'.format(sys.argv[0], datetime.datetime.now(), int(time.time() - start_time)))


if __name__ == '__main__':
    main()
