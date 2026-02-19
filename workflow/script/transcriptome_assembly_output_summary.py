#!/usr/bin/env python
# coding: utf-8

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import time

import numpy
import pandas


def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir_transcriptome_assembly', metavar='PATH', type=str, required=True, help='')
    parser.add_argument(
        '--mode',
        metavar='MODE',
        type=str,
        choices=['auto', 'sraid', 'fastq', 'metadata'],
        required=True,
        help='',
    )
    parser.add_argument('--out', metavar='PATH', type=str, required=True, help='')
    parser.add_argument('--ncpu', metavar='INT', default=1, type=int, help='Number of worker threads.')
    # Backward-compatible alias.
    parser.add_argument('--threads', dest='ncpu', metavar='INT', type=int, help=argparse.SUPPRESS)
    return parser


def sorted_entries(path):
    return sorted(os.listdir(path))


def collect_species_from_sra_list(base_dir):
    input_dir = os.path.join(base_dir, 'input_sra_list')
    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f'Directory not found: {input_dir}')
    return [filename.split('.')[0] for filename in sorted_entries(input_dir) if filename.endswith('.txt')]


def collect_species_from_fastq(base_dir):
    input_dir = os.path.join(base_dir, 'input_fastq')
    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f'Directory not found: {input_dir}')
    out = []
    for dirname in sorted_entries(input_dir):
        dir_path = os.path.join(input_dir, dirname)
        if os.path.isdir(dir_path):
            out.append(dirname)
    return out


def collect_species_from_metadata(base_dir):
    input_dir = os.path.join(base_dir, 'amalgkit_metadata')
    if not os.path.isdir(input_dir):
        raise FileNotFoundError(f'Directory not found: {input_dir}')
    return [filename.split('.')[0] for filename in sorted_entries(input_dir) if filename.endswith('.metadata.tsv')]


def collect_input_species(base_dir, mode):
    if mode == 'sraid':
        return collect_species_from_sra_list(base_dir)
    if mode == 'fastq':
        return collect_species_from_fastq(base_dir)
    if mode == 'metadata':
        return collect_species_from_metadata(base_dir)

    # auto mode: pick the first source with non-empty entries
    for collector in [collect_species_from_sra_list, collect_species_from_fastq, collect_species_from_metadata]:
        try:
            out = collector(base_dir)
        except FileNotFoundError:
            continue
        if out:
            return out
    return []


def collect_species_ids_for_subdir(base_dir, subdir):
    subdir_path = os.path.realpath(os.path.join(base_dir, subdir))
    files = sorted_entries(subdir_path)
    if subdir == 'amalgkit_getfastq':
        species_ids_in_files = set()
        for f in files:
            if f.endswith('.safely_removed.txt'):
                species_id = f.split('.')[0]
                species_ids_in_files.add(species_id)
            else:
                path = os.path.join(subdir_path, f)
                if os.path.isdir(path) and os.listdir(path):
                    species_ids_in_files.add(f)
        return species_ids_in_files
    return {f'_'.join(f.split('.', 1)[0].split('_')[:2]) for f in files}


def run(args):
    print('args:', vars(args))
    start = time.time()

    base_dir = args.dir_transcriptome_assembly
    input_species_list = collect_input_species(base_dir, args.mode)

    df = pandas.DataFrame(index=input_species_list)
    df.index.name = 'species'
    df['SGE_TASK_ID'] = numpy.arange(1, df.shape[0] + 1)

    subdirs = sorted_entries(base_dir)
    subdirs = [sd for sd in subdirs if os.path.isdir(os.path.join(base_dir, sd))]
    df = pandas.concat([df, pandas.DataFrame(data=0, index=df.index, columns=subdirs, dtype=int)], axis=1)
    species_index_set = set(df.index)
    candidate_subdirs = [sd for sd in subdirs if sd not in ['multispecies_summary', 'tmp']]
    species_ids_per_subdir = {}
    if args.ncpu > 1 and len(candidate_subdirs) > 1:
        max_workers = min(args.ncpu, len(candidate_subdirs))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_map = {
                executor.submit(collect_species_ids_for_subdir, base_dir, subdir): subdir
                for subdir in candidate_subdirs
            }
            for future in as_completed(future_map):
                subdir = future_map[future]
                species_ids_per_subdir[subdir] = future.result()
    else:
        for subdir in candidate_subdirs:
            species_ids_per_subdir[subdir] = collect_species_ids_for_subdir(base_dir, subdir)

    for subdir in subdirs:
        if subdir in ['multispecies_summary', 'tmp']:
            continue
        species_ids_in_files = species_ids_per_subdir.get(subdir, set())
        species_ids_with_file = list(species_ids_in_files & species_index_set)
        df.loc[species_ids_with_file, subdir] = 1
        num_missing = (df.loc[:, subdir] == 0).sum()
        txt = 'Subdirectory {}: {:,} / {:,} files are missing.'
        print(txt.format(subdir, num_missing, df.shape[0]))

    df = pandas.concat([df, pandas.DataFrame(data=0, index=df.index, columns=['safely_removed.txt'], dtype=int)], axis=1)
    flagdir_path = os.path.realpath(os.path.join(base_dir, 'amalgkit_getfastq'))
    if os.path.isdir(flagdir_path):
        flag_files = [f for f in sorted_entries(flagdir_path) if f.endswith('.safely_removed.txt')]
        flag_species_ids = {f.split('.')[0] for f in flag_files}
        flag_species_ids_with_file = list(flag_species_ids & set(df.index))
        df.loc[flag_species_ids_with_file, 'safely_removed.txt'] = 1
        num_missing = (df['safely_removed.txt'] == 0).sum()
        print(
            'amalgkit_getfastq/*.safely_removed.txt: {:,} / {:,} files are missing.'.format(
                num_missing, df.shape[0]
            )
        )

    incomplete_ids = df.loc[df['safely_removed.txt'] == 0, 'SGE_TASK_ID'].astype(int).tolist()
    print('Incomplete job IDs (based on amalgkit_getfastq/*.safely_removed.txt):', ','.join(str(x) for x in sorted(incomplete_ids)))

    print('Writing output file:', args.out, flush=True)
    df.to_csv(args.out, index=True, sep='\t')
    print('Done. Elapsed time: {:,} sec'.format(int(time.time() - start)))


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    args.ncpu = max(1, int(args.ncpu))
    run(args)


if __name__ == '__main__':
    main()
