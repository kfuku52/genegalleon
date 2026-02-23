#!/usr/bin/env python3
# coding: utf-8

import argparse
import os
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd


def read_busco_table(path_to_table, col_names):
    table = pd.read_table(
        path_to_table,
        sep='\t',
        header=None,
        comment='#',
        names=col_names,
    )
    table.loc[:, 'sequence'] = table.loc[:, 'sequence'].str.replace(':[-\\.0-9]*$', '', regex=True)
    for col in ['sequence', 'orthodb_url', 'description']:
        table[col] = table[col].fillna('-').astype('string')
    return table


def process_species_table(path_to_table, species_infile, col_names):
    tmp_table = read_busco_table(path_to_table, col_names)

    tmp_metadata = tmp_table.loc[:, ['busco_id', 'orthodb_url', 'description']]
    tmp_metadata = tmp_metadata.drop_duplicates(keep='first', inplace=False, ignore_index=True)

    tmp_sequence = (
        tmp_table.loc[:, ['busco_id', 'sequence']]
        .groupby('busco_id', as_index=False)
        .agg(sequence=('sequence', ','.join))
        .rename(columns={'sequence': species_infile})
    )
    return species_infile, tmp_sequence, tmp_metadata


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--busco_outdir', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--outfile', metavar='PATH', default='busco_summary.tsv', type=str, help='')
    parser.add_argument('--ncpu', metavar='INT', default=1, type=int, help='Number of worker threads.')
    return parser


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    args.ncpu = max(1, int(args.ncpu))
    print("Starting collect_common_BUSCO_gene.py")

    col_names = ['busco_id', 'status', 'sequence', 'score', 'length', 'orthodb_url', 'description']
    species_infiles = sorted(
        [
            f
            for f in os.listdir(path=args.busco_outdir)
            if (not f.startswith('.'))
            and f.endswith('.tsv')
            and os.path.isfile(os.path.join(args.busco_outdir, f))
        ]
    )

    sequence_tables = {}
    metadata_tables = {}

    tasks = []
    for species_infile in species_infiles:
        path_to_table = os.path.join(args.busco_outdir, species_infile)
        if not os.path.exists(path_to_table):
            warnings.warn(f'full_table.tsv does not exist. Skipping: {species_infile}')
            continue
        tasks.append((species_infile, path_to_table))

    if len(tasks) == 0:
        pd.DataFrame(columns=['busco_id', 'orthodb_url', 'description']).to_csv(args.outfile, sep='\t', index=None)
        print("Ending collect_common_BUSCO_gene.py")
        raise SystemExit(0)

    if args.ncpu == 1:
        for species_infile, path_to_table in tasks:
            print(f'Working on {species_infile}', flush=True)
            species, tmp_sequence, tmp_metadata = process_species_table(path_to_table, species_infile, col_names)
            sequence_tables[species] = tmp_sequence
            metadata_tables[species] = tmp_metadata
    else:
        max_workers = min(args.ncpu, len(tasks))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(process_species_table, path_to_table, species_infile, col_names): species_infile
                for species_infile, path_to_table in tasks
            }
            for future in as_completed(futures):
                species_infile = futures[future]
                print(f'Working on {species_infile}', flush=True)
                species, tmp_sequence, tmp_metadata = future.result()
                sequence_tables[species] = tmp_sequence
                metadata_tables[species] = tmp_metadata

    if len(sequence_tables) == 0:
        pd.DataFrame(columns=['busco_id', 'orthodb_url', 'description']).to_csv(args.outfile, sep='\t', index=None)
        print("Ending collect_common_BUSCO_gene.py")
        raise SystemExit(0)

    ordered_species = [s for s in species_infiles if s in sequence_tables]
    merged_sequence = (
        pd.concat(
            [sequence_tables[species].set_index('busco_id') for species in ordered_species],
            axis=1,
            sort=False,
        )
        .reset_index()
    )

    metadata_frames = [metadata_tables[species] for species in ordered_species if species in metadata_tables]
    if len(metadata_frames) > 0:
        merged_metadata = pd.concat(metadata_frames, ignore_index=True)
        merged_metadata.loc[:, ['orthodb_url', 'description']] = (
            merged_metadata.loc[:, ['orthodb_url', 'description']].replace('-', pd.NA)
        )
        merged_metadata = (
            merged_metadata
            .groupby('busco_id', as_index=False, sort=False)[['orthodb_url', 'description']]
            .first()
            .fillna('-')
        )
    else:
        merged_metadata = pd.DataFrame(columns=['busco_id', 'orthodb_url', 'description'])

    merged_table = merged_metadata.merge(merged_sequence, on='busco_id', how='outer')
    species_columns = [c for c in species_infiles if c in merged_table.columns]

    for species_col in species_infiles:
        if species_col not in merged_table.columns:
            merged_table[species_col] = '-'
    if len(species_columns) > 0:
        merged_table.loc[:, species_infiles] = merged_table.loc[:, species_infiles].fillna('-').astype('string')
        keep_rows = (merged_table.loc[:, species_infiles] != '-').any(axis=1)
        dropped_rows = int((~keep_rows).sum())
        if dropped_rows > 0:
            print(f'Dropping {dropped_rows} BUSCO rows with no genes in any species.', flush=True)
        merged_table = merged_table.loc[keep_rows, :]

    for col in ['orthodb_url', 'description']:
        if col not in merged_table.columns:
            merged_table[col] = '-'
    merged_table.loc[:, ['orthodb_url', 'description']] = merged_table.loc[:, ['orthodb_url', 'description']].fillna('-').astype('string')

    column_order = ['busco_id', 'orthodb_url', 'description'] + species_infiles
    merged_table = merged_table.loc[:, column_order]
    merged_table = merged_table.sort_values('busco_id', kind='mergesort').reset_index(drop=True)

    merged_table.to_csv(args.outfile, sep='\t', index=None)
    print("Ending collect_common_BUSCO_gene.py")


if __name__ == '__main__':
    main()
