#!/usr/bin/env python
# coding: utf-8

import argparse
import datetime
import glob
import os
import re
import time
import gc
from concurrent.futures import ThreadPoolExecutor, as_completed
import math
import logging
from tqdm import tqdm
import pandas as pd
import sqlalchemy
from tenacity import retry, wait_fixed, stop_after_attempt, retry_if_exception_type

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler("generate_orthogroup_database.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Set pandas display options
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 1000)

MAX_SQL_VARIABLES = 999


def read_header_columns(file_path):
    """
    Read the first line of a TSV file and return column names.
    """
    with open(file_path, "r", encoding="utf-8", errors="replace") as file_handle:
        header_line = file_handle.readline().strip()
    if not header_line:
        return []
    return header_line.split("\t")

def optimize_sqlite(engine):
    """
    Optimize SQLite settings and set a busy timeout.
    """
    with engine.begin() as conn:
        conn.execute(sqlalchemy.text("PRAGMA journal_mode = OFF;"))
        conn.execute(sqlalchemy.text("PRAGMA synchronous = OFF;"))
        conn.execute(sqlalchemy.text("PRAGMA cache_size = 100000;"))  # Adjust as needed
        conn.execute(sqlalchemy.text("PRAGMA temp_store = MEMORY;"))
        conn.execute(sqlalchemy.text("PRAGMA locking_mode = EXCLUSIVE;"))
        conn.execute(sqlalchemy.text("PRAGMA busy_timeout = 30000;"))
    logger.info("SQLite PRAGMA settings optimized for performance and busy timeout set.")

def calculate_chunksize(num_columns, max_sql_vars=MAX_SQL_VARIABLES):
    return max(1, math.floor(max_sql_vars / num_columns))

def initialize_buffers(infiles, columns):
    return {stat: [] for stat in infiles.keys()}

def create_indexes(engine, tables):
    """
    Create indexes after all data has been inserted and committed.
    """
    with engine.begin() as conn:
        for table in tables:
            try:
                index_name = f"idx_orthogroup_{table}"
                conn.execute(sqlalchemy.text(f"CREATE INDEX IF NOT EXISTS {index_name} ON {table} (orthogroup);"))
                logger.info(f"Created index '{index_name}' on table '{table}'.")
            except Exception as e:
                logger.error(f"Failed to create index on table '{table}': {e}")

def validate_directories(required_dirs, db_path):
    for dir_path in required_dirs:
        if not os.path.isdir(dir_path):
            logger.error(f"Directory does not exist: {dir_path}")
            exit(1)
    db_dir = os.path.dirname(db_path)
    if db_dir and not os.path.exists(db_dir):
        try:
            os.makedirs(db_dir)
            logger.info(f"Created database directory: {db_dir}")
        except Exception as e:
            logger.error(f"Failed to create database directory '{db_dir}': {e}")
            exit(1)

@retry(
    wait=wait_fixed(2),
    stop=stop_after_attempt(3),
    retry=retry_if_exception_type(Exception),
    reraise=True
)
def process_files(file_path, columns_to_read, available_cols_set=None):
    """
    Read a TSV file, add the 'orthogroup' column, ensure any missing columns become NaN,
    and return the trimmed DataFrame with the columns we want to keep (including 'orthogroup').
    If any required columns are missing, return an empty DataFrame.
    """
    try:
        # Derive orthogroup name from the file name.
        og = os.path.splitext(os.path.basename(file_path))[0].split('.')[0]

        # Define the desired columns to read (exclude 'orthogroup' because it isn’t in the file)
        desired_cols = [col for col in columns_to_read if col != 'orthogroup']

        # Use caller-provided header columns when available to avoid a second parse.
        if available_cols_set is None:
            available_cols_set = set(read_header_columns(file_path))
        else:
            available_cols_set = set(available_cols_set)

        # Filter desired columns to only those that exist in the file.
        filtered_cols = [c for c in desired_cols if c in available_cols_set]

        # Check if any desired columns are missing from the file
        missing_cols = [c for c in desired_cols if c not in available_cols_set]
        if missing_cols:
            # Convert the set to a string and truncate if it is too long
            msg = str(missing_cols)
            if len(msg) > 300:
                msg = msg[:300] + "..."
            logger.warning(
                f"Missing columns {msg} in file '{file_path}'. Skipping this file."
            )
            return pd.DataFrame()  # Return empty DataFrame

        # Read only the available columns using the filtered list.
        df = pd.read_csv(
            file_path, 
            sep="\t", 
            header=0, 
            usecols=filtered_cols, 
            low_memory=True
        )

        # --- Clean null characters ---
        # Apply cleaning only on object-type columns so that numeric columns remain unaffected.
        #object_cols = df.select_dtypes(include=[object]).columns
        #for col in object_cols:
        #    df[col] = df[col].str.replace('\x00', '')

        # Insert the 'orthogroup' column derived from the file name.
        df["orthogroup"] = og

        # Reorder the columns to ensure 'orthogroup' is first.
        # Use columns_to_read which includes 'orthogroup' and all expected columns
        df = df[columns_to_read]

        return df

    except Exception as e:
        logger.error(f"Error processing file {file_path}: {e}")
        raise

def apply_cutoff(df, cutoff_stat_str):
    try:
        cutoff_stats = [s.strip().replace('\'', '').replace('\"', '') for s in cutoff_stat_str.split('|')]
        for stat in cutoff_stats:
            stat_name, stat_value = stat.split(',')
            stat_value = float(stat_value)
            if stat_name in df.columns:
                df = df[df[stat_name].astype(float).fillna(0) >= stat_value]
        return df
    except Exception as e:
        logger.error(f"Error applying cutoff: {e}")
        return df

def main():
    parser = argparse.ArgumentParser(description="Optimize performance for database population script.")
    parser.add_argument('--overwrite', metavar='bool', default=0, type=int, help='Overwrite existing database if set to 1.')
    parser.add_argument('--dbpath', metavar='PATH', default='', type=str, help='Path to the SQLite database.')
    parser.add_argument('--dir_stat_tree', metavar='PATH', default='', type=str, help='Directory for stat_tree files.')
    parser.add_argument('--dir_stat_branch', metavar='PATH', default='', type=str, help='Directory for stat_branch files.')
    parser.add_argument('--dir_csubst_cb_prefix', metavar='PATH', default='', type=str, help='Prefix path for csubst_cb directories.')
    parser.add_argument('--row_threshold', metavar='INT', default=10000, type=int, help='Number of rows to accumulate before inserting into SQL.')
    parser.add_argument('--cb_categories', metavar='CAT1,CAT2,...', default='any2any,any2spe', type=str, help='CSUBST cb stat categories to incorporate.')
    parser.add_argument('--cutoff_stat', metavar='STAT1,VALUE1|STAT2,VALUE2|...', default='OCNany2spe,0.8', type=str, help='Cutoff statistics for filtering.')
    parser.add_argument('--ncpu', dest='max_workers', metavar='INT', default=4, type=int, help='Number of worker threads.')
    # Backward-compatible aliases.
    parser.add_argument('--threads', dest='max_workers', metavar='INT', type=int, help=argparse.SUPPRESS)
    parser.add_argument('--max_workers', dest='max_workers', metavar='INT', type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    logger.info("Starting the orthogroup database generation script.")

    params = vars(args)
    params['max_workers'] = max(1, int(params['max_workers']))
    db_path = params['dbpath']

    cb_categories = [cat.strip() for cat in args.cb_categories.split(',')]
    all_cb_categories = ['any2any','any2spe','spe2any','spe2spe','any2dif','dif2any','spe2dif','dif2spe','dif2dif']
    cb_remove_categories = list(set(all_cb_categories).difference(set(cb_categories)))

    required_dirs = [
        params['dir_stat_tree'],
        params['dir_stat_branch'],
    ]
    validate_directories(required_dirs, db_path)

    if params['overwrite'] and os.path.exists(db_path):
        try:
            os.remove(db_path)
            logger.info(f"Existing database '{db_path}' removed due to overwrite flag.")
        except Exception as e:
            logger.error(f"Failed to remove existing database '{db_path}': {e}")
            exit(1)

    engine = sqlalchemy.create_engine(
        f"sqlite:///{db_path}",
        poolclass=sqlalchemy.pool.QueuePool,
        echo=False,
        future=True,
        connect_args={"timeout": 30},
        pool_size=params['max_workers'],
        max_overflow=0
    )
    optimize_sqlite(engine)

    # We'll not connect yet; we will connect within transactions below.

    # Gather input directories
    infiles = {}
    columns = {}
    num_columns = {}
    header_columns_by_file = {}
    header_columns_set_by_file = {}
    indirs = {
        'tree': params['dir_stat_tree'],
        'branch': params['dir_stat_branch'],
    }
    # Identify csubst directories
    cb_dirs = [d for d in glob.glob(params['dir_csubst_cb_prefix'] + '*') if not d.endswith('csubst.cb_stats')]
    for cb_dir in cb_dirs:
        if os.path.isdir(cb_dir) and len(os.listdir(cb_dir)) > 0:
            logger.info(f"CSUBST higher-order convergence directory detected: {cb_dir}")
            arity = re.sub('.*_', '', cb_dir)
            table_name = f'cb{arity}'
            indirs[table_name] = cb_dir
    logger.info(f"Input directories to be appended: {', '.join(indirs.values())}")

    # Check column names of all input files
    for stat in indirs.keys():
        dir_path = indirs[stat]
        if not os.path.exists(dir_path):
            logger.warning(f"Directory does not exist. Skipping: {dir_path}")
            continue
        infiles[stat] = [f for f in os.listdir(dir_path) if not f.startswith('.')]
        logger.info(f"Number of infiles for '{stat}': {len(infiles[stat])}")
        num_columns[stat] = []
        header_columns_by_file[stat] = {}
        header_columns_set_by_file[stat] = {}
        column_names_set = set()
        for infile in infiles[stat]:
            file_path = os.path.join(dir_path, infile)
            try:
                infile_columns = read_header_columns(file_path)
                infile_column_names_set = set(infile_columns)
                header_columns_by_file[stat][infile] = infile_columns
                header_columns_set_by_file[stat][infile] = infile_column_names_set
                column_names_set.update(infile_column_names_set)
                num_columns[stat].append(len(infile_column_names_set))
            except Exception as e:
                logger.error(f"Error reading header from {file_path}: {e}")
                header_columns_by_file[stat][infile] = []
                header_columns_set_by_file[stat][infile] = set()
                num_columns[stat].append(0)

        # Sort files so that the one with the greatest number of columns is first
        sorted_pairs = sorted(
            zip(infiles[stat], num_columns[stat]),
            key=lambda item: item[1],
            reverse=True
        )
        infiles[stat] = [infile for infile, _ in sorted_pairs]
        if not infiles[stat]:
            logger.warning(f"No valid files found for '{stat}'. Skipping.")
            continue

        max_columns_file = os.path.join(dir_path, infiles[stat][0])
        try:
            max_columns = header_columns_by_file[stat].get(infiles[stat][0], [])
            logger.info(f"Max columns file: {max_columns_file}")
            logger.info(f"Number of all columns in input for '{stat}': {len(column_names_set)}")
            logger.info(f"Max number of columns in input tables for '{stat}': {len(max_columns)}")
            
            # Initialize columns list
            if stat.startswith('cb'):
                filtered_columns = [
                    col for col in max_columns
                    if not any(remove_cat in col for remove_cat in cb_remove_categories)
                ]
                additional_columns = list(column_names_set - set(max_columns))
                additional_filtered = [
                    col for col in additional_columns 
                    if not any(remove_cat in col for remove_cat in cb_remove_categories)
                ]
                columns[stat] = ['orthogroup'] + filtered_columns + additional_filtered
            else:
                columns[stat] = ['orthogroup'] \
                                + max_columns \
                                + list(column_names_set - set(max_columns))
            max_col_len = 300 # Upper limit to detect malformed column names like '\x00\x00\x00\x00...'
            columns[stat] = [col for col in columns[stat] if (len(col) <= max_col_len)]
            logger.info(f"Number of all columns for '{stat}': {len(columns[stat])}")
            preview_cols = columns[stat][:20]
            suffix = ' ...' if len(columns[stat]) > 20 else ''
            logger.info(f"First columns for '{stat}': {preview_cols}{suffix}")
        except Exception as e:
            logger.error(f"Error reading max columns file '{max_columns_file}': {e}")
            columns[stat] = ['orthogroup']
            logger.warning(f"Falling back to minimal columns for '{stat}'.")

    logger.info(f"{datetime.datetime.today()}: Started adding infiles to the database.")

    buffers = initialize_buffers(infiles, columns)
    buffer_row_counts = {stat: 0 for stat in infiles.keys()}
    processed_files = {stat: 0 for stat in infiles.keys()}
    total_files = {stat: len(files) for stat, files in infiles.items()}

    chunksizes = {}
    for stat, cols in columns.items():
        num_cols = len(cols)
        chunksizes[stat] = calculate_chunksize(num_cols)
    for stat, cs in chunksizes.items():
        if cs < 1:
            chunksizes[stat] = 1

    # Process files concurrently
    futures = {}
    with ThreadPoolExecutor(max_workers=params['max_workers']) as executor:
        for stat, files in infiles.items():
            for infile in files:
                file_path = os.path.join(indirs[stat], infile)
                if os.path.getsize(file_path) == 0:
                    logger.warning(f"Skipping empty file: {file_path}")
                    processed_files[stat] += 1
                    continue
                future = executor.submit(
                    process_files,
                    file_path,
                    columns[stat],
                    header_columns_set_by_file[stat].get(infile),
                )
                futures[future] = (stat, file_path)

        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing files"):
            stat, file_path = futures[future]
            try:
                df = future.result()
                # If it's a csubst table, apply cutoff
                if stat.startswith('cb'):
                    df = apply_cutoff(df, params['cutoff_stat'])
                if not df.empty:
                    buffers[stat].append(df)
                    buffer_row_counts[stat] += len(df)
                processed_files[stat] += 1
                # Check if buffer exceeds threshold and insert into DB
                if buffer_row_counts[stat] >= params['row_threshold']:
                    full_df = pd.concat(buffers[stat], ignore_index=True)
                    if not full_df.empty:
                        with engine.begin() as conn:
                            full_df.to_sql(
                                name=stat,
                                con=conn,
                                if_exists="append",
                                index=False,
                                dtype=None,
                                chunksize=chunksizes[stat],
                                method='multi'
                            )
                        remaining = total_files[stat] - processed_files[stat]
                        logger.info(f"{datetime.datetime.today()}: {stat}: Inserted {buffer_row_counts[stat]} rows. Files done: {processed_files[stat]}, remaining: {remaining}")
                    buffers[stat] = []
                    buffer_row_counts[stat] = 0
                    gc.collect()
            except Exception as e:
                logger.error(f"Error processing file {file_path}: {e}")

    # Insert any remaining rows for each table
    for stat, buffer_list in buffers.items():
        if buffer_row_counts[stat] > 0 and buffer_list:
            full_df = pd.concat(buffer_list, ignore_index=True)
            if not full_df.empty:
                with engine.begin() as conn:
                    full_df.to_sql(
                        name=stat,
                        con=conn,
                        if_exists="append",
                        index=False,
                        dtype=None,
                        chunksize=chunksizes[stat],
                        method='multi'
                    )
                logger.info(f"{datetime.datetime.today()}: Inserted remaining {buffer_row_counts[stat]} rows for '{stat}'.")
            else:
                logger.info(f"No rows to insert for '{stat}' (buffer empty).")
            buffers[stat] = []
            buffer_row_counts[stat] = 0
            gc.collect()

    logger.info(f"{datetime.datetime.today()}: Completed adding infiles to the database.")

    # Retrieve table info
    with engine.begin() as conn:
        try:
            tables = pd.read_sql_query(sql=sqlalchemy.text("SELECT name FROM sqlite_master WHERE type='table'"), con=conn)['name'].values
            logger.info(f"Existing tables before indexing: {tables}")
        except Exception as e:
            logger.error(f"Failed to retrieve tables after insertion: {e}")
            tables = []

    # Create indexes on the new tables
    create_indexes(engine, tables)

    # Show column info
    with engine.begin() as conn:
        for table in tables:
            try:
                columns_info = pd.read_sql_query(sql=sqlalchemy.text(f"PRAGMA TABLE_INFO({table})"), con=conn)
                logger.info(f"\n{table}\n{columns_info[['name', 'type']].to_string(index=False)}")
            except Exception as e:
                logger.error(f"Failed to retrieve columns from '{table}': {e}")

    # Dispose engine to close all connections
    engine.dispose()
    logger.info("All database operations completed and engine disposed.")

if __name__ == "__main__":
    start_time = time.time()
    main()
    elapsed_time = time.time() - start_time
    logger.info(f"Elapsed time: {elapsed_time:,.1f} [sec]")
    logger.info(f"{datetime.datetime.today()}: Database generation completed!")
