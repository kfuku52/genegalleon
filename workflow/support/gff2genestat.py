#!/usr/bin/env python3

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime
import gzip
import numpy
import os
import pandas
import re
import sys
import time

pandas.options.mode.chained_assignment = None

MATCH_ATTRIBUTE_KEYS = {
    '',
    'ID',
    'Parent',
    'Name',
    'gene',
    'geneName',
    'gene_id',
    'locus_tag',
    'protein_id',
    'transcript',
    'transcript_id',
}

FALLBACK_FEATURES = {
    'mRNA',
    'transcript',
    'J_gene_segment',
    'V_gene_segment',
    'C_gene_segment',
    'D_gene_segment',
    'gene',
}

NAMESPACE_SUFFIX_DELIMITERS = {
    ':',
    '|',
}

ATTRIBUTE_KEY_ALIASES = {
    'gene_name': 'geneName',
}

def build_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', metavar='intron_num|gene_delim', default='', type=str, help='',
                        choices=['intron_num','gene_delim'])
    parser.add_argument('--dir_gff', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--seqfile', metavar='PATH', default='', type=str, help='')
    parser.add_argument('--outfile', metavar='PATH', default='gff2genestat.tsv', type=str, help='')
    parser.add_argument('--feature', metavar='STR', default='CDS', type=str, help='')
    parser.add_argument('--multiple_hits', metavar='STR', default='longest', type=str, help='')
    parser.add_argument('--ncpu', metavar='INT', default=1, type=int, help='Number of worker threads.')
    return parser

def read_fasta_seqname(file_path):
    seqnames = []
    opener = gzip.open if file_path.endswith('.gz') else open
    with opener(file_path, 'rt', encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('>'):
                seqnames.append(line[1:].rstrip('\n'))
    return pandas.Series(seqnames)

def get_gene_names(seq_names):
    return seq_names.astype(str).str.split('_', n=2).str[-1]


def is_identifier_char(char):
    return char.isalnum() or char == '_'


def is_match_better(candidate, current):
    if candidate is None:
        return False
    if current is None:
        return True
    candidate_key = (candidate[1], -candidate[2], -candidate[3])
    current_key = (current[1], -current[2], -current[3])
    return candidate_key > current_key


def choose_better_match(current, candidate):
    if is_match_better(candidate, current):
        return candidate
    return current


def build_search_term_lookup(seq_names):
    if isinstance(seq_names, pandas.Series):
        seq_name_values = seq_names.fillna('').astype(str).tolist()
    else:
        seq_name_values = [str(seq_name) for seq_name in seq_names]
    lookup = {}
    min_len = None
    max_len = 0
    for order, seq_name in enumerate(seq_name_values):
        gene_name = seq_name.split('_', 2)[-1] if seq_name != '' else ''
        ub_gene_name = gene_name.replace('-', '_')
        term_specs = (
            (seq_name, 2),
            (gene_name, 1),
            (ub_gene_name, 0),
        )
        for term, priority in term_specs:
            term = str(term).strip()
            if term == '':
                continue
            current = lookup.get(term)
            candidate = (seq_name, priority, order)
            if current is None or priority > current[1] or (priority == current[1] and order < current[2]):
                lookup[term] = candidate
            term_len = len(term)
            min_len = term_len if min_len is None else min(min_len, term_len)
            max_len = max(max_len, term_len)
    if min_len is None:
        min_len = 0
    return lookup, min_len, max_len


def normalize_attribute_value(raw_value):
    value = raw_value if isinstance(raw_value, str) else str(raw_value)
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {'"', "'"}:
        return value[1:-1]
    if value[:1] in {'"', "'"}:
        value = value[1:]
    if value[-1:] in {'"', "'"}:
        value = value[:-1]
    return value


def unique_values(values):
    if len(values) <= 1:
        return tuple(values)
    seen = set()
    out = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return tuple(out)


def match_value_to_gene_id(value, lookup, min_len, max_len, value_cache):
    if value in value_cache:
        return value_cache[value]
    value = str(value).strip()
    if value == '':
        value_cache[value] = None
        return None

    exact_hit = lookup.get(value)
    if exact_hit is not None:
        result = (exact_hit[0], len(value), exact_hit[1], exact_hit[2])
        value_cache[value] = result
        return result

    delimiters = [idx for idx, char in enumerate(value) if not is_identifier_char(char)]
    best = None
    best_len = 0
    seen_segments = set()
    segments = [value]
    segments.extend(
        value[idx + 1:]
        for idx in delimiters
        if value[idx] in NAMESPACE_SUFFIX_DELIMITERS
    )

    for segment in segments:
        if segment == '' or segment in seen_segments:
            continue
        seen_segments.add(segment)
        segment_len = len(segment)
        if segment_len < min_len:
            continue
        if best_len > 0 and segment_len <= best_len:
            continue
        segment_delims = [idx for idx, char in enumerate(segment) if not is_identifier_char(char)]
        candidate_ends = [segment_len]
        candidate_ends.extend(idx for idx in reversed(segment_delims) if idx > 0)
        for end in candidate_ends:
            candidate_len = end
            if candidate_len < min_len:
                break
            if max_len > 0 and candidate_len > max_len:
                continue
            if best_len > 0 and candidate_len <= best_len:
                break
            candidate = segment[:end]
            hit = lookup.get(candidate)
            if hit is None:
                continue
            match = (hit[0], candidate_len, hit[1], hit[2])
            best = choose_better_match(best, match)
            if best is not None:
                best_len = best[1]
            break
    value_cache[value] = best
    return best


def parse_attribute_summary(attr_text, lookup, min_len, max_len, value_cache, skip_id_match_when_parent=False):
    attr_id = ''
    parents = []
    values = []
    for raw_field in str(attr_text).split(';'):
        field = raw_field.strip()
        if field == '':
            continue
        equal_pos = field.find('=')
        space_pos = field.find(' ')
        if equal_pos != -1 and (space_pos == -1 or equal_pos < space_pos):
            key = field[:equal_pos]
            raw_value = field[equal_pos + 1:]
        elif space_pos > 0:
            key = field[:space_pos]
            raw_value = field[space_pos + 1:]
        else:
            continue
        key = ATTRIBUTE_KEY_ALIASES.get(key, key)
        if key not in MATCH_ATTRIBUTE_KEYS:
            continue
        if key == 'Parent':
            for raw_parent in raw_value.split(','):
                value = normalize_attribute_value(raw_parent)
                if value != '':
                    parents.append(value)
                    values.append((key, value))
            continue
        value = normalize_attribute_value(raw_value)
        if value == '':
            continue
        if key == 'ID' and attr_id == '':
            attr_id = value
        values.append((key, value))

    best = None
    has_parents = len(parents) > 0
    for key, value in values:
        if key not in MATCH_ATTRIBUTE_KEYS:
            continue
        if skip_id_match_when_parent and has_parents and key == 'ID':
            continue
        best = choose_better_match(
            best,
            match_value_to_gene_id(
                value=value,
                lookup=lookup,
                min_len=min_len,
                max_len=max_len,
                value_cache=value_cache,
            ),
        )
    return attr_id, tuple(parents), best


def parse_attribute_fields(attr_text):
    attr_id = ''
    parents = []
    match_values = []
    for raw_field in str(attr_text).split(';'):
        field = raw_field.strip()
        if field == '':
            continue
        equal_pos = field.find('=')
        space_pos = field.find(' ')
        if equal_pos != -1 and (space_pos == -1 or equal_pos < space_pos):
            key = field[:equal_pos]
            raw_value = field[equal_pos + 1:]
        elif space_pos > 0:
            key = field[:space_pos]
            raw_value = field[space_pos + 1:]
        else:
            continue
        key = ATTRIBUTE_KEY_ALIASES.get(key, key)
        if key not in MATCH_ATTRIBUTE_KEYS:
            continue
        if key == 'Parent':
            for raw_parent in raw_value.split(','):
                value = normalize_attribute_value(raw_parent)
                if value != '':
                    parents.append(value)
            continue
        value = normalize_attribute_value(raw_value)
        if value == '':
            continue
        if key == 'ID' and attr_id == '':
            attr_id = value
        match_values.append(value)
    return attr_id, unique_values(parents), unique_values(match_values)


def match_values_to_gene_id(values, lookup, min_len, max_len, value_cache):
    best = None
    for value in values:
        best = choose_better_match(
            best,
            match_value_to_gene_id(
                value=value,
                lookup=lookup,
                min_len=min_len,
                max_len=max_len,
                value_cache=value_cache,
            ),
        )
    return best


def merge_id_info(existing, parents, match, values=()):
    if existing is None:
        return {
            'parents': parents,
            'match': match,
            'values': values,
        }
    merged_parents = existing['parents'] if len(parents) == 0 else unique_values(existing['parents'] + parents)
    merged_match = choose_better_match(existing['match'], match)
    existing_values = existing.get('values', ())
    merged_values = existing_values if len(values) == 0 else unique_values(existing_values + values)
    return {'parents': merged_parents, 'match': merged_match, 'values': merged_values}


def resolve_feature_match(feature_id, id_info, resolved_cache, active_stack, lookup, min_len, max_len, value_cache):
    if feature_id == '':
        return None
    cached = resolved_cache.get(feature_id)
    if cached is not None or feature_id in resolved_cache:
        return cached
    if feature_id in active_stack:
        resolved_cache[feature_id] = None
        return None
    info = id_info.get(feature_id)
    if info is None:
        resolved_cache[feature_id] = None
        return None
    active_stack.add(feature_id)
    best = info['match']
    for parent_id in info['parents']:
        parent_match = resolve_feature_match(
            parent_id,
            id_info,
            resolved_cache,
            active_stack,
            lookup,
            min_len,
            max_len,
            value_cache,
        )
        best = choose_better_match(best, parent_match)
    if best is None and len(info.get('values', ())) > 0:
        best = match_values_to_gene_id(
            values=info['values'],
            lookup=lookup,
            min_len=min_len,
            max_len=max_len,
            value_cache=value_cache,
        )
    active_stack.remove(feature_id)
    resolved_cache[feature_id] = best
    return best

def extract_by_ids(gff, seq_names, feature, multiple_hits):
    print('Extracting gene IDs: {}'.format(datetime.datetime.now()), flush=True)
    lookup, min_len, max_len = build_search_term_lookup(seq_names)
    gff_feat = gff.loc[(gff.loc[:,'feature']==feature),:].copy()
    if gff_feat.shape[0] == 0:
        return gff_feat.assign(gene_id='')

    value_cache = {}
    fallback_mask = gff.loc[:, 'feature'].isin(FALLBACK_FEATURES)
    id_info = {}
    resolved_cache = {}
    if fallback_mask.any():
        for attr_text in gff.loc[fallback_mask, 'attributes'].fillna('').astype(str).tolist():
            attr_id, parents, match_values = parse_attribute_fields(attr_text)
            if attr_id == '' and len(parents) == 0 and len(match_values) == 0:
                continue
            if attr_id == '':
                continue
            match = None
            if len(parents) == 0 and len(match_values) > 0:
                match = match_values_to_gene_id(
                    values=match_values,
                    lookup=lookup,
                    min_len=min_len,
                    max_len=max_len,
                    value_cache=value_cache,
                )
            existing = id_info.get(attr_id)
            if existing is None:
                id_info[attr_id] = {'parents': parents, 'match': match, 'values': match_values}
            else:
                id_info[attr_id] = merge_id_info(existing, parents, match, match_values)

    gene_ids = []
    for attr_text in gff_feat.loc[:, 'attributes'].fillna('').astype(str).tolist():
        _, parents, match_values = parse_attribute_fields(attr_text)
        best = None
        for parent_id in parents:
            parent_match = resolve_feature_match(
                parent_id,
                id_info,
                resolved_cache,
                set(),
                lookup,
                min_len,
                max_len,
                value_cache,
            )
            best = choose_better_match(best, parent_match)
        if best is None and len(parents) > 0:
            best = match_values_to_gene_id(
                values=parents,
                lookup=lookup,
                min_len=min_len,
                max_len=max_len,
                value_cache=value_cache,
            )
        if best is None and len(match_values) > 0:
            best = match_values_to_gene_id(
                values=match_values,
                lookup=lookup,
                min_len=min_len,
                max_len=max_len,
                value_cache=value_cache,
            )
        gene_ids.append(best[0] if best is not None else '')

    gff_feat.loc[:, 'gene_id'] = gene_ids
    out = gff_feat.loc[gff_feat.loc[:, 'gene_id'] != '', :]
    if (out.shape[0]==0):
        print('No match was found.')
    return out

def add_id_column(gff, seq_names, new_col='gene_id'):
    print('Adding gene id column: {}'.format(datetime.datetime.now()), flush=True)
    gff = gff.copy()
    lookup, min_len, max_len = build_search_term_lookup(seq_names)
    value_cache = {}
    gene_ids = []
    for attr_text in gff.loc[:, 'attributes'].fillna('').astype(str).tolist():
        _, _, match = parse_attribute_summary(attr_text, lookup, min_len, max_len, value_cache)
        gene_ids.append(match[0] if match is not None else '')
    gff.loc[:,new_col] = gene_ids
    return gff


def summarize_gene_features(gff, out_cols, id_col='gene_id'):
    gene_ids = gff.loc[:, id_col].fillna('').astype(str).to_numpy()
    if gene_ids.size == 0:
        return pandas.DataFrame(columns=out_cols)
    sequences = gff.loc[:, 'sequence'].to_numpy()
    strands = gff.loc[:, 'strand'].to_numpy()
    starts = gff.loc[:, 'start'].to_numpy()
    ends = gff.loc[:, 'end'].to_numpy()

    order = []
    by_gene = {}
    for gene_id, sequence, strand, start, end in zip(gene_ids, sequences, strands, starts, ends):
        if gene_id == '':
            continue
        start = int(start)
        end = int(end)
        feature_len = end - start + 1
        rec = by_gene.get(gene_id)
        if rec is None:
            by_gene[gene_id] = {
                'gene_id': gene_id,
                'feature_size': feature_len,
                'num_intron': 0,
                'intron_offsets': [],
                'chromosome': sequence,
                'start': start,
                'end': end,
                'strand': strand,
            }
            order.append(gene_id)
            continue
        rec['intron_offsets'].append(rec['feature_size'])
        rec['feature_size'] += feature_len
        rec['num_intron'] += 1
        rec['end'] = end

    if len(order) == 0:
        return pandas.DataFrame(columns=out_cols)
    rows = []
    for gene_id in order:
        rec = by_gene[gene_id]
        intron_offsets = rec['intron_offsets']
        max_intron_pos = int(intron_offsets[-1]) if len(intron_offsets) > 0 else 0
        if (max_intron_pos > rec['feature_size']):
            txt = 'Intron position cannot be greater than feature size: {}, feature_size={}, max intron position = {}'
            raise Exception(txt.format(gene_id, rec['feature_size'], max_intron_pos))
        rows.append({
            'gene_id': gene_id,
            'feature_size': int(rec['feature_size']),
            'num_intron': int(rec['num_intron']),
            'intron_positions': ';'.join(str(int(pos)) for pos in intron_offsets),
            'chromosome': rec['chromosome'],
            'start': int(rec['start']),
            'end': int(rec['end']),
            'strand': rec['strand'],
        })
    return pandas.DataFrame(rows, columns=out_cols)

def add_intron_info(gff, df_all, id_col='gene_id'):
    print('Adding intron information: {}'.format(datetime.datetime.now()), flush=True)
    seq_names = gff.loc[:,id_col].unique()
    df_intron = pandas.DataFrame(None, index=numpy.arange(len(seq_names)), columns=df_all.columns)
    for j,seq_name in enumerate(seq_names):
        gff_gene = gff.loc[(gff.loc[:,id_col]==seq_name),:].reset_index()
        feature_size = ((gff_gene.loc[:,'end']+1).values - gff_gene.loc[:,'start'].values).sum()
        num_intron = gff_gene.shape[0]-1
        intron_positions = list()
        current_offset = 0
        max_intron_pos = 0
        for i in numpy.arange(gff_gene.shape[0]-1):
            feature_block_size = gff_gene.at[i,'end'] - gff_gene.at[i,'start'] + 1
            intron_pos = current_offset + feature_block_size
            max_intron_pos = intron_pos if (intron_pos > max_intron_pos) else max_intron_pos
            intron_positions.append(str(intron_pos))
            current_offset += feature_block_size
        str_intron_positions = ';'.join(intron_positions)
        df_intron.at[j,'gene_id'] = seq_name
        df_intron.at[j,'feature_size'] = feature_size
        df_intron.at[j,'num_intron'] = num_intron
        df_intron.at[j,'intron_positions'] = str_intron_positions
        if (max_intron_pos > feature_size):
            txt = 'Intron position cannot be greater than feature size: {}, feature_size={}, max intron position = {}'
            raise Exception(txt.format(seq_name, feature_size, max_intron_pos))
    df_all = pandas.concat([df_all, df_intron], axis='index', ignore_index=True)
    return df_all

def add_gene_delim(gff, df_all, id_col='gene_id'):
    print('Adding gene location: {}'.format(datetime.datetime.now()), flush=True)
    seq_names = gff.loc[:,id_col].unique()
    for seq_name in seq_names:
        is_seqname = (df_all.loc[:,id_col]==seq_name)
        gff_gene = gff.loc[(gff.loc[:,id_col]==seq_name),:].reset_index(drop=True)
        df_all.loc[is_seqname,'chromosome'] = gff_gene.at[0,'sequence']
        df_all.loc[is_seqname,'strand'] = gff_gene.at[0,'strand']
        df_all.loc[is_seqname,'start'] = gff_gene.at[0,'start']
        df_all.loc[is_seqname,'end'] = gff_gene.loc[:,'end'].tail(1).values[0]
    return df_all

def filename2sciname(file_name):
    tmp = file_name.replace('_','|',1)
    tmp = re.sub(r'[-_.,].*', '', tmp)
    sci_name = tmp.replace('|', '_')
    return sci_name


def process_single_gff(gff_file, dir_gff, seq_sp_values, feature, multiple_hits, gff_cols, out_cols):
    print("{}: Started processing: {}".format(datetime.datetime.now(), gff_file), flush=True)
    gff_path = os.path.join(dir_gff, gff_file)
    if os.stat(gff_path).st_size == 0:
        sys.stderr.write('Empty file: {}\n'.format(gff_path))
        return pandas.DataFrame(columns=out_cols)
    gff = pandas.read_csv(gff_path, sep='\t', header=None, comment='#', low_memory=False, quoting=3)
    gff.columns = gff_cols
    seq_sp = pandas.Series(seq_sp_values)
    gff_id = extract_by_ids(gff=gff, seq_names=seq_sp, feature=feature, multiple_hits=multiple_hits)
    if gff_id.shape[0] == 0:
        return pandas.DataFrame(columns=out_cols)
    print('Summarizing gene features: {}'.format(datetime.datetime.now()), flush=True)
    return summarize_gene_features(gff=gff_id, out_cols=out_cols)

def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    args.ncpu = max(1, int(args.ncpu))
    start_time = time.time()
    print('gff2genestat.py started: {}'.format(datetime.datetime.now()))

    out_cols = ['gene_id', 'feature_size', 'num_intron', 'intron_positions', 'chromosome', 'start', 'end','strand']
    gff_cols = ['sequence','source','feature','start','end','score','strand','phase','attributes']
    seq_names = read_fasta_seqname(file_path=args.seqfile)
    gff_files = [gff for gff in os.listdir(args.dir_gff) if not gff.startswith('.')]
    gff_files = [gff for gff in gff_files if gff.endswith(('.gff', '.gtf', '.gff3', '.gff.gz', '.gtf.gz', '.gff3.gz'))]
    tasks = []
    for gff_file in gff_files:
        sp_ub = filename2sciname(file_name=gff_file)
        seq_sp = seq_names[seq_names.str.startswith(sp_ub)]
        if seq_sp.shape[0] == 0:
            continue
        tasks.append((gff_file, seq_sp.tolist()))

    frames = []
    if args.ncpu == 1 or len(tasks) <= 1:
        for gff_file, seq_sp_values in tasks:
            df_tmp = process_single_gff(
                gff_file=gff_file,
                dir_gff=args.dir_gff,
                seq_sp_values=seq_sp_values,
                feature=args.feature,
                multiple_hits=args.multiple_hits,
                gff_cols=gff_cols,
                out_cols=out_cols,
            )
            if df_tmp.shape[0] > 0:
                frames.append(df_tmp)
    else:
        max_workers = min(args.ncpu, len(tasks))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    process_single_gff,
                    gff_file,
                    args.dir_gff,
                    seq_sp_values,
                    args.feature,
                    args.multiple_hits,
                    gff_cols,
                    out_cols,
                ): (idx, gff_file)
                for idx, (gff_file, seq_sp_values) in enumerate(tasks)
            }
            ordered_results = {}
            for future in as_completed(futures):
                idx, gff_file = futures[future]
                try:
                    df_tmp = future.result()
                    ordered_results[idx] = df_tmp
                except Exception as exc:
                    raise RuntimeError('Failed processing {}: {}'.format(gff_file, exc)) from exc
            for idx in sorted(ordered_results.keys()):
                df_tmp = ordered_results[idx]
                if df_tmp.shape[0] > 0:
                    frames.append(df_tmp)

    if len(frames) > 0:
        df_all = pandas.concat(frames, axis='index', ignore_index=True)
    else:
        df_all = pandas.DataFrame(columns=out_cols)

    df_all = df_all.drop_duplicates()
    num_input = len(seq_names)
    num_output = df_all.shape[0]
    print('Number of input genes: {}'.format(num_input), flush=True)
    print('Number of output entries: {}'.format(num_output), flush=True)
    if (num_input != num_output):
        observed_gene_ids = set(df_all['gene_id'].fillna('').astype(str).tolist())
        for seq_name in sorted(set(seq_names.astype(str).tolist()) - observed_gene_ids):
            print('Missing in the output file: {}'.format(seq_name), flush=True)
    df_all.to_csv(args.outfile, sep='\t', header=True, index=False)
    print('gff2genestat.py completed in {:,} secs: {}'.format(int(time.time() - start_time), datetime.datetime.now()))


if __name__ == '__main__':
    main()
