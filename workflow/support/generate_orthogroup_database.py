#!/usr/bin/env python3
# coding: utf-8

import argparse
import datetime
import gc
import glob
import logging
import math
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **_kwargs):
        return iterable
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from gene_family_output_store import GeneFamilyOutputStore

try:
    import sqlalchemy
except ImportError:
    sqlalchemy = None
try:
    from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_fixed
except ImportError:
    def retry(*_args, **_kwargs):
        def _decorator(func):
            return func
        return _decorator

    def wait_fixed(_seconds):
        return None

    def stop_after_attempt(_attempts):
        return None

    def retry_if_exception_type(_exception_type):
        return None

logger = logging.getLogger(__name__)


def configure_logging(log_file_path="generate_orthogroup_database.log"):
    handlers = [logging.StreamHandler()]
    if log_file_path:
        handlers.insert(0, logging.FileHandler(log_file_path))
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=handlers,
        force=True,
    )

# Set pandas display options
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 1000)

MAX_SQL_VARIABLES = 999
AA_CHANGE_TABLE = "aa_change"
AA_CHANGE_UNIT_TABLE = "aa_change_unit"
AA_CHANGE_FDR_PVALUE_COLUMNS = {
    "p_rate_enrichment": "q_rate_enrichment_global",
    "p_rate_enrichment_empirical": "q_rate_enrichment_empirical_global",
    "p_rate_enrichment_empirical_maxT": "q_rate_enrichment_empirical_maxT_global",
}
CSUBST_SCAN_BASELINE_COLUMNS = {
    AA_CHANGE_TABLE: {
        "site_rate_categorized",
        "q_rate_enrichment_empirical",
        "q_rate_enrichment_empirical_by_trait",
        "q_rate_enrichment_empirical_by_trait_match",
    },
    AA_CHANGE_UNIT_TABLE: {
        "fg_clade_branch_ids",
    },
}


def require_sqlalchemy():
    if sqlalchemy is None:
        raise ImportError("sqlalchemy is required to build the orthogroup database.")


def read_header_columns(file_path, store=None, logical_subdir=None, logical_name=None):
    """
    Read the first line of a TSV file and return column names.
    """
    if store is None:
        with open(file_path, "r", encoding="utf-8", errors="replace") as file_handle:
            header_line = file_handle.readline().strip()
    else:
        with store.open_binary(logical_subdir, logical_name) as file_handle:
            header_line = file_handle.readline().decode("utf-8", errors="replace").strip()
    if not header_line:
        return []
    return header_line.split("\t")


def visible_entries(path):
    return [entry for entry in os.listdir(path) if not entry.startswith(".")]


def visible_files(path):
    return [entry for entry in visible_entries(path) if os.path.isfile(os.path.join(path, entry))]


def has_visible_entries(path):
    return len(visible_entries(path)) > 0


def discover_csubst_cb_dirs(prefix):
    if not prefix:
        return []
    return [
        path
        for path in glob.glob(prefix + '*')
        if not path.endswith('csubst_cb_stats')
    ]


def validate_csubst_scan_schemas(scan_dirs, store=None):
    """
    Require current semantic marker columns while allowing optional columns to vary.

    The marker columns distinguish the current scan semantics from legacy
    outputs. Additional columns are discovered dynamically and need not be the
    same across every per-family TSV.
    """
    problems = []
    for table_name, scan_dir in scan_dirs:
        logical_subdir = os.path.basename(os.path.normpath(scan_dir)) if scan_dir else ""
        if not scan_dir:
            continue
        if store is None:
            if not os.path.isdir(scan_dir):
                continue
            files = visible_files(scan_dir)
        else:
            if logical_subdir not in store.logical_subdirs():
                continue
            files = store.file_names(logical_subdir)
        if not files:
            continue

        required_columns = CSUBST_SCAN_BASELINE_COLUMNS[table_name]
        for infile in files:
            file_path = os.path.join(scan_dir, infile)
            header_columns = read_header_columns(
                file_path,
                store=store,
                logical_subdir=logical_subdir,
                logical_name=infile,
            )
            duplicate_columns = sorted({col for col in header_columns if header_columns.count(col) > 1})
            if duplicate_columns:
                problems.append(
                    f"{file_path}: duplicate columns: {', '.join(duplicate_columns)}"
                )
                continue

            header_set = frozenset(header_columns)
            missing_columns = sorted(required_columns.difference(header_set))
            if missing_columns:
                problems.append(
                    f"{file_path}: missing current-schema columns: {', '.join(missing_columns)}"
                )

    if problems:
        details = "\n  - ".join(problems)
        raise ValueError(
            "Unsupported legacy CSUBST scan TSV schema. GeneGalleon requires the current "
            "semantic marker columns but accepts any number of additional columns. "
            "Regenerate legacy CSUBST scan outputs with the current "
            f"csubst version before rebuilding the database.\n  - {details}"
        )


def optimize_sqlite(engine):
    """
    Optimize SQLite settings and set a busy timeout.
    """
    require_sqlalchemy()
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
    require_sqlalchemy()
    with engine.begin() as conn:
        for table in tables:
            try:
                index_name = f"idx_orthogroup_{table}"
                conn.execute(sqlalchemy.text(f"CREATE INDEX IF NOT EXISTS {index_name} ON {table} (orthogroup);"))
                logger.info(f"Created index '{index_name}' on table '{table}'.")
            except Exception as e:
                logger.error(f"Failed to create index on table '{table}': {e}")


def quote_sql_identifier(identifier):
    if not re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", str(identifier)):
        raise ValueError(f"Unsafe SQL identifier: {identifier}")
    return f'"{identifier}"'


def table_exists(conn, table_name):
    query = sqlalchemy.text("SELECT 1 FROM sqlite_master WHERE type='table' AND name=:name")
    return conn.execute(query, {"name": table_name}).fetchone() is not None


def table_columns(conn, table_name):
    table_sql = quote_sql_identifier(table_name)
    info = pd.read_sql_query(sql=sqlalchemy.text(f"PRAGMA TABLE_INFO({table_sql})"), con=conn)
    if info.empty:
        return []
    return info["name"].tolist()


def calculate_bh_fdr(pvalues):
    pvalues = pd.to_numeric(pd.Series(pvalues), errors="coerce").to_numpy(dtype=float)
    qvalues = np.full(shape=pvalues.shape, fill_value=np.nan, dtype=float)
    finite = np.isfinite(pvalues)
    if not finite.any():
        return qvalues
    finite_index = np.flatnonzero(finite)
    finite_p = np.clip(pvalues[finite], 0.0, 1.0)
    order = np.argsort(finite_p, kind="mergesort")
    ranked = finite_p[order]
    ranks = np.arange(1, ranked.shape[0] + 1, dtype=float)
    ranked_q = ranked * ranked.shape[0] / ranks
    ranked_q = np.minimum.accumulate(ranked_q[::-1])[::-1]
    ranked_q = np.clip(ranked_q, 0.0, 1.0)
    qvalues[finite_index[order]] = ranked_q
    return qvalues


def add_global_aa_change_fdr_columns(engine, table_name=AA_CHANGE_TABLE):
    require_sqlalchemy()
    with engine.begin() as conn:
        if not table_exists(conn, table_name):
            logger.info(f"Skipping global FDR calculation because table '{table_name}' does not exist.")
            return []
        columns = table_columns(conn, table_name)
        pvalue_columns = [col for col in AA_CHANGE_FDR_PVALUE_COLUMNS if col in columns]
        if not pvalue_columns:
            logger.info(f"Skipping global FDR calculation because '{table_name}' has no recognized P-value columns.")
            return []

        table_sql = quote_sql_identifier(table_name)
        for p_col in pvalue_columns:
            q_col = AA_CHANGE_FDR_PVALUE_COLUMNS[p_col]
            if q_col not in columns:
                conn.execute(sqlalchemy.text(f"ALTER TABLE {table_sql} ADD COLUMN {quote_sql_identifier(q_col)} REAL"))
                logger.info(f"Added global FDR column '{q_col}' to table '{table_name}'.")

        select_cols = ["rowid AS _rowid"] + [quote_sql_identifier(col) for col in pvalue_columns]
        df = pd.read_sql_query(
            sql=sqlalchemy.text(f"SELECT {', '.join(select_cols)} FROM {table_sql}"),
            con=conn,
        )
        if df.empty:
            logger.info(f"Table '{table_name}' is empty; global FDR columns were added without row updates.")
            return [AA_CHANGE_FDR_PVALUE_COLUMNS[col] for col in pvalue_columns]

        update_df = pd.DataFrame({"_rowid": df["_rowid"].astype(int)})
        for p_col in pvalue_columns:
            q_col = AA_CHANGE_FDR_PVALUE_COLUMNS[p_col]
            update_df[q_col] = calculate_bh_fdr(df[p_col])

        temp_table = "_tmp_aa_change_global_fdr"
        update_df.to_sql(temp_table, con=conn, if_exists="replace", index=False, chunksize=calculate_chunksize(update_df.shape[1]), method="multi")
        conn.execute(sqlalchemy.text(f"CREATE INDEX IF NOT EXISTS idx_{temp_table}_rowid ON {quote_sql_identifier(temp_table)} (_rowid);"))
        temp_sql = quote_sql_identifier(temp_table)
        for q_col in update_df.columns:
            if q_col == "_rowid":
                continue
            q_col_sql = quote_sql_identifier(q_col)
            conn.execute(sqlalchemy.text(
                f"UPDATE {table_sql} "
                f"SET {q_col_sql} = (SELECT {q_col_sql} FROM {temp_sql} WHERE {temp_sql}._rowid = {table_sql}.rowid) "
                f"WHERE rowid IN (SELECT _rowid FROM {temp_sql});"
            ))
        conn.execute(sqlalchemy.text(f"DROP TABLE {temp_sql};"))
        logger.info(
            f"Calculated global BH FDR for {update_df.shape[0]:,} '{table_name}' rows using columns: {', '.join(pvalue_columns)}"
        )
        return [AA_CHANGE_FDR_PVALUE_COLUMNS[col] for col in pvalue_columns]

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
def process_files(
    file_path,
    columns_to_read,
    available_cols_set=None,
    fill_missing_columns=False,
    store=None,
    logical_subdir=None,
    logical_name=None,
):
    """
    Read a TSV file, add the 'orthogroup' column, ensure any missing columns become NaN,
    and return the trimmed DataFrame with the columns we want to keep (including 'orthogroup').
    If required columns are missing, raise an error unless fill_missing_columns
    is enabled. The latter is used for CSUBST scan tables, whose optional output
    columns may vary between files and csubst releases.
    """
    try:
        # Derive the gene-family ID from the file name. This script is used for
        # both orthogroup and query2family outputs, so strip GeneGalleon suffixes
        # before falling back to the legacy dot-delimited convention.
        og = gene_family_id_from_path(file_path)

        # Define the desired columns to read (exclude 'orthogroup' because it isn’t in the file)
        desired_cols = [col for col in columns_to_read if col != 'orthogroup']

        # Use caller-provided header columns when available to avoid a second parse.
        if available_cols_set is None:
            available_cols_set = set(
                read_header_columns(
                    file_path,
                    store=store,
                    logical_subdir=logical_subdir,
                    logical_name=logical_name,
                )
            )
        else:
            available_cols_set = set(available_cols_set)

        # Filter desired columns to only those that exist in the file.
        filtered_cols = [c for c in desired_cols if c in available_cols_set]

        # Check if any desired columns are missing from the file
        missing_cols = [c for c in desired_cols if c not in available_cols_set]
        if missing_cols and not fill_missing_columns:
            preview = ", ".join(missing_cols[:20])
            if len(missing_cols) > 20:
                preview += f", ... ({len(missing_cols)} total)"
            raise ValueError(
                f"Missing required columns in '{file_path}': {preview}"
            )

        # Read only the available columns using the filtered list.
        if store is None:
            df = pd.read_csv(
                file_path,
                sep="\t",
                header=0,
                usecols=filtered_cols,
                low_memory=True,
            )
        else:
            with store.open_binary(logical_subdir, logical_name) as file_handle:
                df = pd.read_csv(
                    file_handle,
                    sep="\t",
                    header=0,
                    usecols=filtered_cols,
                    low_memory=True,
                )

        if fill_missing_columns:
            for col in missing_cols:
                df[col] = np.nan

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

    except Exception:
        logger.exception("Error processing file %s", file_path)
        raise


def gene_family_id_from_path(file_path):
    stem = os.path.splitext(os.path.basename(file_path))[0]
    for suffix in (
        '_stat.branch',
        '_stat.tree',
        '_csubst_scan_units',
        '_csubst_scan',
        '_csubst_cb_stats',
        '.stat.branch',
        '.stat.tree',
        '.csubst_scan_units',
        '.csubst_scan',
        '.csubst_cb_stats',
    ):
        if stem.endswith(suffix):
            return stem[:-len(suffix)]
    for marker in (
        '.csubst_cb_',
        '_csubst_cb_',
    ):
        if marker in stem:
            return stem.split(marker, 1)[0]
    return stem.split('.')[0]

def parse_cutoff_stat(cutoff_stat):
    parsed = []
    if cutoff_stat is None:
        return parsed

    if isinstance(cutoff_stat, (list, tuple)):
        tokens = cutoff_stat
    else:
        tokens = [s.strip() for s in str(cutoff_stat).split('|')]

    for token in tokens:
        if isinstance(token, (list, tuple)):
            if len(token) != 2:
                continue
            stat_name = str(token[0]).strip().replace("'", "").replace('"', "")
            stat_value_raw = token[1]
        else:
            token_str = str(token).strip().replace("'", "").replace('"', "")
            if not token_str or ',' not in token_str:
                continue
            stat_name, stat_value_raw = token_str.split(',', 1)
            stat_name = stat_name.strip()
            stat_value_raw = stat_value_raw.strip()

        if not stat_name:
            continue
        try:
            stat_value = float(stat_value_raw)
        except (TypeError, ValueError):
            continue
        parsed.append((stat_name, stat_value))

    return parsed


def apply_cutoff(df, cutoff_stat):
    try:
        cutoff_stats = parse_cutoff_stat(cutoff_stat)
        for stat_name, stat_value in cutoff_stats:
            if stat_name in df.columns:
                values = pd.to_numeric(df[stat_name], errors='coerce').fillna(0)
                df = df[values >= stat_value]
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
    parser.add_argument('--dir_csubst_aa_change', metavar='PATH', default='', type=str, help='Directory for csubst scan candidate state-change files.')
    parser.add_argument('--dir_csubst_aa_change_unit', metavar='PATH', default='', type=str, help='Directory for csubst scan foreground-unit files.')
    parser.add_argument('--dir_gene_family', metavar='PATH', default='', type=str, help='Optional query2family/orthogroup root used to read live and ZIP-backed logical subdirectories.')
    parser.add_argument('--row_threshold', metavar='INT', default=10000, type=int, help='Number of rows to accumulate before inserting into SQL.')
    parser.add_argument('--cb_categories', metavar='CAT1,CAT2,...', default='any2any,any2spe', type=str, help='CSUBST cb stat categories to incorporate.')
    parser.add_argument('--cutoff_stat', metavar='STAT1,VALUE1|STAT2,VALUE2|...', default='OCNany2spe,0.8', type=str, help='Cutoff statistics for filtering.')
    parser.add_argument('--ncpu', dest='max_workers', metavar='INT', default=4, type=int, help='Number of worker threads.')
    args = parser.parse_args()
    configure_logging()
    require_sqlalchemy()
    logger.info("Starting the orthogroup database generation script.")

    params = vars(args)
    params['max_workers'] = max(1, int(params['max_workers']))
    db_path = params['dbpath']
    output_store = (
        GeneFamilyOutputStore(params['dir_gene_family'])
        if params['dir_gene_family']
        else None
    )

    cb_categories = [cat.strip() for cat in args.cb_categories.split(',')]
    all_cb_categories = ['any2any','any2spe','spe2any','spe2spe','any2dif','dif2any','spe2dif','dif2spe','dif2dif']
    cb_remove_categories = list(set(all_cb_categories).difference(set(cb_categories)))

    required_dirs = [
        params['dir_stat_tree'],
        params['dir_stat_branch'],
    ]
    if output_store is None:
        validate_directories(required_dirs, db_path)
    else:
        logical_subdirs = set(output_store.logical_subdirs())
        for required_dir in required_dirs:
            required_subdir = os.path.basename(os.path.normpath(required_dir))
            if required_subdir not in logical_subdirs:
                raise FileNotFoundError(
                    f"Logical input subdirectory does not exist: {required_subdir} "
                    f"under {params['dir_gene_family']}"
                )
        db_dir = os.path.dirname(db_path)
        if db_dir:
            os.makedirs(db_dir, exist_ok=True)

    scan_dirs = [
        (AA_CHANGE_TABLE, params['dir_csubst_aa_change']),
        (AA_CHANGE_UNIT_TABLE, params['dir_csubst_aa_change_unit']),
    ]
    try:
        validate_csubst_scan_schemas(scan_dirs, store=output_store)
    except ValueError as exc:
        logger.error(str(exc))
        raise SystemExit(1) from exc

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
    if output_store is None:
        cb_dirs = discover_csubst_cb_dirs(params['dir_csubst_cb_prefix'])
    else:
        cb_prefix = os.path.basename(os.path.normpath(params['dir_csubst_cb_prefix']))
        cb_dirs = [
            os.path.join(params['dir_gene_family'], subdir)
            for subdir in output_store.logical_subdirs()
            if subdir.startswith(cb_prefix) and subdir != "csubst_cb_stats"
        ]
    for cb_dir in cb_dirs:
        cb_subdir = os.path.basename(os.path.normpath(cb_dir))
        cb_has_entries = (
            os.path.isdir(cb_dir) and has_visible_entries(cb_dir)
            if output_store is None
            else bool(output_store.file_names(cb_subdir))
        )
        if cb_has_entries:
            logger.info(f"CSUBST higher-order convergence directory detected: {cb_dir}")
            arity = re.sub('.*_', '', cb_dir)
            table_name = f'cb{arity}'
            indirs[table_name] = cb_dir
    for table_name, scan_dir in scan_dirs:
        scan_subdir = os.path.basename(os.path.normpath(scan_dir)) if scan_dir else ""
        scan_has_entries = (
            bool(scan_dir) and os.path.isdir(scan_dir) and has_visible_entries(scan_dir)
            if output_store is None
            else bool(scan_dir) and bool(output_store.file_names(scan_subdir))
        )
        if scan_has_entries:
            logger.info(f"CSUBST scan directory detected for table '{table_name}': {scan_dir}")
            indirs[table_name] = scan_dir
    logger.info(f"Input directories to be appended: {', '.join(indirs.values())}")

    # Check column names of all input files
    for stat in indirs.keys():
        dir_path = indirs[stat]
        logical_subdir = os.path.basename(os.path.normpath(dir_path))
        if output_store is None and not os.path.exists(dir_path):
            logger.warning(f"Directory does not exist. Skipping: {dir_path}")
            continue
        infiles[stat] = (
            visible_files(dir_path)
            if output_store is None
            else output_store.file_names(logical_subdir)
        )
        logger.info(f"Number of infiles for '{stat}': {len(infiles[stat])}")
        num_columns[stat] = []
        header_columns_by_file[stat] = {}
        header_columns_set_by_file[stat] = {}
        column_names_set = set()
        for infile in infiles[stat]:
            file_path = os.path.join(dir_path, infile)
            try:
                infile_columns = read_header_columns(
                    file_path,
                    store=output_store,
                    logical_subdir=logical_subdir,
                    logical_name=infile,
                )
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

    # Process files concurrently. Any unreadable or malformed input makes the
    # complete database build fail; a partial database must never be reported
    # as a successful result.
    futures = {}
    processing_errors = []
    with ThreadPoolExecutor(max_workers=params['max_workers']) as executor:
        for stat, files in infiles.items():
            for infile in files:
                file_path = os.path.join(indirs[stat], infile)
                logical_subdir = os.path.basename(os.path.normpath(indirs[stat]))
                artifact = (
                    None
                    if output_store is None
                    else output_store.artifact(logical_subdir, infile)
                )
                file_size = (
                    os.path.getsize(file_path)
                    if output_store is None
                    else int(artifact.size or 0) if artifact is not None else 0
                )
                if file_size == 0:
                    processing_errors.append((file_path, "input file is empty"))
                    continue
                future = executor.submit(
                    process_files,
                    file_path,
                    columns[stat],
                    header_columns_set_by_file[stat].get(infile),
                    stat in CSUBST_SCAN_BASELINE_COLUMNS,
                    output_store,
                    logical_subdir,
                    infile,
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
                processing_errors.append((file_path, str(e)))

    if processing_errors:
        engine.dispose()
        details = "\n  - ".join(
            f"{file_path}: {message}" for file_path, message in processing_errors
        )
        raise RuntimeError(
            "Orthogroup database generation failed; no partial build will be "
            f"reported as successful.\n  - {details}"
        )

    incomplete_tables = {
        stat: (processed_files[stat], total_files[stat])
        for stat in total_files
        if processed_files[stat] != total_files[stat]
    }
    if incomplete_tables:
        engine.dispose()
        details = ", ".join(
            f"{stat}={processed}/{total}"
            for stat, (processed, total) in sorted(incomplete_tables.items())
        )
        raise RuntimeError(f"Input-file accounting mismatch: {details}")

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

    add_global_aa_change_fdr_columns(engine)

    with engine.begin() as conn:
        try:
            tables = pd.read_sql_query(sql=sqlalchemy.text("SELECT name FROM sqlite_master WHERE type='table'"), con=conn)['name'].values
            logger.info(f"Existing tables after global FDR calculation: {tables}")
        except Exception as e:
            logger.error(f"Failed to retrieve tables after FDR calculation: {e}")
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
