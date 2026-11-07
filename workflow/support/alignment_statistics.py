"""Read and append family alignment statistics from live or archived outputs."""

import os
from concurrent.futures import ThreadPoolExecutor

import numpy
import pandas
from parallel_io import bounded_results


def _alignment_stats_columns():
    return [
        "No_of_taxa",
        "Alignment_length",
        "Total_matrix_cells",
        "Undetermined_characters",
        "Missing_percent",
        "No_variable_sites",
        "Parsimony_informative_sites",
        "GC_content",
    ]


def _visible_entries(path):
    return [entry for entry in os.listdir(path) if not entry.startswith(".")]


def _read_alignment_stats_file(file_path, og_id, alignment_stats_cols):
    tmp = pandas.read_csv(file_path, sep="\t", header=0, usecols=alignment_stats_cols, nrows=1, low_memory=False)
    # usecols filters columns but retains their order in the input header.
    return og_id, tmp.loc[:, alignment_stats_cols].iloc[0].to_list()


def _read_alignment_stats_store_file(store, subdir, file_name, og_id, alignment_stats_cols):
    with store.open_binary(subdir, file_name) as handle:
        tmp = pandas.read_csv(handle, sep="\t", header=0, usecols=alignment_stats_cols, nrows=1, low_memory=False)
    return og_id, tmp.loc[:, alignment_stats_cols].iloc[0].to_list()


def append_alignment_stats(df, dir_alignment_stats, extension, ncpu, extract_id, store=None, logical_subdir=None, *, columns_before_check=False):
    alignment_stats_cols = _alignment_stats_columns()
    alignment_stats_new_cols = [f"{col}_{extension}" for col in alignment_stats_cols]

    def initialize_columns():
        for column in alignment_stats_new_cols:
            df[column] = numpy.nan

    if columns_before_check:
        initialize_columns()
    if store is None and not os.path.isdir(dir_alignment_stats):
        print(f"{extension}: {dir_alignment_stats} was not found. Skipping.", flush=True)
        return df, False
    if store is not None and logical_subdir not in store.logical_subdirs():
        print(f"{extension}: logical subdirectory {logical_subdir} was not found. Skipping.", flush=True)
        return df, False
    if not columns_before_check:
        initialize_columns()

    files = sorted(_visible_entries(dir_alignment_stats)) if store is None else store.file_names(logical_subdir)
    queued = []
    seen_ogs = set()
    valid_ogs = set(df.index.astype(str))
    for file in files:
        file_path = os.path.join(dir_alignment_stats, file)
        if store is None and not os.path.isfile(file_path):
            continue
        og_id = extract_id(file)
        if og_id is None:
            continue
        if og_id in seen_ogs:
            continue
        if og_id not in valid_ogs:
            continue
        seen_ogs.add(og_id)
        queued.append((file, og_id))

    counter = 0
    result_rows = []
    if ncpu > 1 and len(queued) > 1:
        max_workers = min(ncpu, len(queued))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            if store is None:
                reader = _read_alignment_stats_file
                arguments = ((os.path.join(dir_alignment_stats, file), og_id, alignment_stats_cols)
                             for file, og_id in queued)
            else:
                reader = _read_alignment_stats_store_file
                arguments = ((store, logical_subdir, file, og_id, alignment_stats_cols)
                             for file, og_id in queued)
            for og_id, values in bounded_results(executor, reader, arguments, max_workers * 8):
                result_rows.append((og_id, values))
                counter += 1
    else:
        for file, og_id in queued:
            if store is None:
                file_path = os.path.join(dir_alignment_stats, file)
                og_id, values = _read_alignment_stats_file(file_path, og_id, alignment_stats_cols)
            else:
                og_id, values = _read_alignment_stats_store_file(
                    store,
                    logical_subdir,
                    file,
                    og_id,
                    alignment_stats_cols,
                )
            result_rows.append((og_id, values))
            counter += 1

    if result_rows:
        result_df = pandas.DataFrame.from_records(
            [values for _, values in result_rows],
            index=[og_id for og_id, _ in result_rows],
            columns=alignment_stats_new_cols,
        )
        df.loc[result_df.index, alignment_stats_new_cols] = result_df

    print(f"{extension}: {counter} alignment-statistics results were appended.", flush=True)
    return df, True
