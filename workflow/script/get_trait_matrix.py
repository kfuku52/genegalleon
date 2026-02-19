#! /usr/bin/env python

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import os
import pandas

pandas.options.mode.chained_assignment = None

def read_fasta_seqname(file_path):
    seqnames = []
    if file_path.endswith('.gz'):
        with gzip.open(file_path, 'rt') as f:
            for line in f:
                if line.startswith('>'):
                    seqnames.append(line[1:].strip())
    else:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    seqnames.append(line[1:].strip())
    return seqnames


def get_gene_name(seq_name):
    parts = seq_name.split('_', 2)
    if len(parts) >= 3:
        return parts[2]
    return seq_name


def get_species_name(seq_name):
    parts = seq_name.split('_', 2)
    if len(parts) >= 2:
        return parts[0] + '_' + parts[1]
    return seq_name


def trait_filename_to_species_name(trait_file):
    name = trait_file
    if '.' in name:
        name = name.split('.', 1)[0]
    parts = name.split('_', 2)
    if len(parts) >= 2:
        return parts[0] + '_' + parts[1]
    return name


def process_trait_file(trait_path, search_ids, id_map):
    trait = pandas.read_csv(trait_path, sep='\t', header=0, comment='#')
    if trait.shape[0] == 0:
        return trait
    if trait.index.values[0] != 0:
        trait = trait.reset_index(drop=False)
    trait = trait.rename(columns={trait.columns.to_list()[0]: 'gene_id'})
    trait['gene_id'] = trait['gene_id'].astype(str)
    trait2 = trait.loc[trait['gene_id'].isin(search_ids), :].copy()
    if trait2.shape[0] == 0:
        return trait2
    trait2.loc[:, 'gene_id'] = trait2['gene_id'].map(id_map).fillna(trait2['gene_id'])
    return trait2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir_trait', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--seqfile', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--outfile', metavar='PATH', default='trait.tsv', type=str, help='')
    parser.add_argument('--ncpu', metavar='INT', default=1, type=int, help='Number of worker threads.')
    # Backward-compatible alias.
    parser.add_argument('--threads', dest='ncpu', metavar='INT', type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.ncpu = max(1, int(args.ncpu))
    print('get_trait_matrix.py started.')

    seq_names = read_fasta_seqname(file_path=args.seqfile)
    gene_names = [get_gene_name(sn) for sn in seq_names]
    gene_species_uniq = sorted(set([get_species_name(sn) for sn in seq_names]))
    search_ids = set(seq_names + gene_names)

    id_map = {}
    for seq_name, gene_name in zip(seq_names, gene_names):
        id_map[seq_name] = seq_name
        id_map[gene_name] = seq_name

    trait_files = sorted([tf for tf in os.listdir(args.dir_trait) if not tf.startswith('.')])
    tasks = []
    for trait_file in trait_files:
        species_name = trait_filename_to_species_name(trait_file)
        print("Started processing {}: species name = {}".format(trait_file, species_name))
        if species_name not in gene_species_uniq:
            print('Trait file for {} not found in {}. Skipping'.format(species_name, args.dir_trait))
            continue
        trait_path = os.path.join(args.dir_trait, trait_file)
        tasks.append((trait_file, trait_path))

    frames = []
    if args.ncpu == 1 or len(tasks) <= 1:
        for trait_file, trait_path in tasks:
            trait2 = process_trait_file(trait_path=trait_path, search_ids=search_ids, id_map=id_map)
            if trait2.shape[0] > 0:
                frames.append(trait2)
            print('Finished processing {}: Number of identified genes = {:,}'.format(trait_file, trait2.shape[0]))
    else:
        max_workers = min(args.ncpu, len(tasks))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(process_trait_file, trait_path, search_ids, id_map): (idx, trait_file)
                for idx, (trait_file, trait_path) in enumerate(tasks)
            }
            ordered_results = {}
            for future in as_completed(futures):
                idx, trait_file = futures[future]
                trait2 = future.result()
                ordered_results[idx] = trait2
                print('Finished processing {}: Number of identified genes = {:,}'.format(trait_file, trait2.shape[0]))
            for idx in sorted(ordered_results.keys()):
                trait2 = ordered_results[idx]
                if trait2.shape[0] > 0:
                    frames.append(trait2)

    if len(frames) > 0:
        df_trait_all = pandas.concat(frames, ignore_index=True)
    else:
        df_trait_all = pandas.DataFrame()
    print('Number of input genes: {:,}'.format(len(seq_names)))
    print('Number of output genes with traits: {:,}'.format(df_trait_all.shape[0]))
    df_trait_all.to_csv(args.outfile, sep='\t', header=True, index=False)
    print('get_trait_matrix.py done!')


if __name__ == '__main__':
    main()
