#!/usr/bin/env python3

import argparse
import gzip
import hashlib
import json
import os
import re
from pathlib import Path

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


def read_uniprot_descriptions(path, accessions=None):
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
                identifier = parts[0]
                acc = uniprot_accession(identifier)
                if accessions is None:
                    keys = (acc,)
                else:
                    keys = {acc, identifier} & accessions
                    if not keys:
                        continue
                text = parts[1] if len(parts) > 1 else ''
                if acc != identifier:
                    # Official UniProt FASTA headers carry accession and entry
                    # name; search tools commonly return only the accession.
                    text = re.split(r'\s(?:OS|OX|GN|PE|SV)=', text, maxsplit=1)[0]
                for key in keys:
                    desc[key] = text
    return desc


def uniprot_accession(identifier):
    """Recognize official UniProt tokens without rewriting other literal IDs."""
    match = re.fullmatch(r'(?:sp|tr)\|([^|\s]+)\|[^|\s]+', identifier)
    return match.group(1) if match else identifier


def read_uniprot_metadata(path, accessions=None):
    if (path == '') or (not os.path.exists(path)):
        return pandas.DataFrame(columns=['accession'] + METADATA_COLUMNS)
    header = pandas.read_csv(path, sep='\t', dtype=str, low_memory=False, keep_default_na=False, nrows=0)
    if 'accession' in header.columns:
        accession_col = 'accession'
    elif 'sprot_best' in header.columns:
        accession_col = 'sprot_best'
    else:
        return pandas.DataFrame(columns=['accession'] + METADATA_COLUMNS)
    usecols = [accession_col] + [col for col in METADATA_COLUMNS if col in header.columns]
    options = dict(sep='\t', dtype=str, low_memory=False, keep_default_na=False, usecols=usecols)
    if accessions is None:
        df = pandas.read_csv(path, **options)
    else:
        # Parse the complete TSV with pandas, retaining only relevant rows.
        # Chunking preserves quoted fields, literal IDs and first-duplicate wins
        # without materializing the entire reference table for each family.
        frames = []
        with pandas.read_csv(path, chunksize=16384, **options) as reader:
            for chunk in reader:
                selected = chunk.loc[chunk[accession_col].isin(accessions)].copy()
                if not selected.empty:
                    frames.append(selected)
        df = pandas.concat(frames, ignore_index=True) if frames else pandas.DataFrame(columns=usecols)
    if accession_col != 'accession':
        df = df.rename(columns={accession_col: 'accession'})
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


def _diamond_identifier(value):
    """Keep identifier text literal while retaining missing empty fields."""
    return value if value else float('nan')


def load_best_diamond_hits(path):
    if (path == '') or (not os.path.exists(path)) or (os.path.getsize(path) == 0):
        return pandas.DataFrame(columns=DIAMOND_COLUMNS + ['coverage'])

    identifiers = ('qseqid', 'sseqid')
    df = pandas.read_csv(
        path, sep='\t', header=None, names=DIAMOND_COLUMNS, low_memory=False,
        dtype={column: str for column in DIAMOND_COLUMNS if column not in identifiers},
        converters=dict.fromkeys(identifiers, _diamond_identifier),
    )
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


def refresh_go_annotations(annotation_paths, metadata_path, output_dir, audit_path=''):
    """Refresh only GO fields from saved best accessions, without a new search."""
    go_columns = ['go_ids', 'go_aspects', 'go_terms', 'go_evidence']
    paths = [Path(path) for path in annotation_paths]
    if len({path.name for path in paths}) != len(paths):
        raise ValueError('Input annotation basenames must be unique')
    frames = [pandas.read_csv(path, sep='\t', dtype=str, keep_default_na=False) for path in paths]
    accessions = set()
    for frame in frames:
        if not {'gene_id', 'sprot_best'}.issubset(frame.columns) or frame['gene_id'].duplicated().any():
            raise ValueError('Existing annotations require unique gene_id and a saved sprot_best column')
        accessions.update(uniprot_accession(value) for value in frame['sprot_best'] if value)
    header = pandas.read_csv(metadata_path, sep='\t', nrows=0)
    if not set(go_columns).issubset(header.columns):
        raise ValueError('GO refresh requires all four aligned GO metadata columns from the source DAT')
    metadata = read_uniprot_metadata(metadata_path, accessions)
    missing = accessions - set(metadata['accession'])
    if missing:
        raise ValueError(f'Metadata does not cover {len(missing)} saved best accessions; use the matching source DAT')
    for row in metadata.itertuples(index=False):
        ids = row.go_ids.split(';') if row.go_ids else []
        aspects = row.go_aspects.split(';') if row.go_aspects else []
        terms = row.go_terms.split('; ') if row.go_terms else []
        evidence = row.go_evidence.split(';') if ids else []
        if len({len(ids), len(aspects), len(terms), len(evidence)}) != 1:
            raise ValueError('Unaligned GO metadata; rebuild it from the source DAT')
        if any(not re.fullmatch(r'GO:\d{7}', go_id) for go_id in ids) or any(not term for term in terms):
            raise ValueError('Invalid GO metadata ID or empty term')
        if any(aspect.strip() not in {'C', 'F', 'P', 'CC', 'MF', 'BP'} for aspect in aspects):
            raise ValueError('Invalid GO metadata aspect')
    lookup = metadata.set_index('accession')
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    destinations = [output / path.name for path in paths]
    if any(path.exists() for path in destinations) or (audit_path and Path(audit_path).exists()):
        raise ValueError('Refreshed annotation outputs already exist')
    records = []
    for path, target, frame in zip(paths, destinations, frames, strict=True):
        result = frame.copy()
        keys = frame['sprot_best'].map(uniprot_accession)
        protected = [column for column in frame.columns if column not in go_columns]
        for column in go_columns:
            result[column] = keys.map(lookup[column]).fillna('')
        if not result[protected].equals(frame[protected]):
            raise ValueError('GO refresh changed accepted hits or another annotation field')
        opener = gzip.open if target.name.endswith('.gz') else open
        with opener(target, 'xt') as handle:
            result.to_csv(handle, sep='\t', index=False)
        records.append({'input': str(path), 'output': str(target), 'genes': len(frame),
                        'GO_annotated_genes': int(result['go_ids'].ne('').sum()),
                        'input_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                        'output_sha256': hashlib.sha256(target.read_bytes()).hexdigest(),
                        'protected_columns_preserved': True})
    audit = {'homology_search_repeated': False, 'saved_best_accessions_preserved': True,
             'metadata_sha256': hashlib.sha256(Path(metadata_path).read_bytes()).hexdigest(), 'records': records}
    if audit_path:
        with Path(audit_path).open('x') as handle:
            json.dump(audit, handle, indent=2)
    return audit


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--diamond_tsv', metavar='PATH', default='', type=str)
    parser.add_argument('--query_fasta', metavar='PATH', default='', type=str)
    parser.add_argument('--uniprot_fasta', metavar='PATH', default='', type=str)
    parser.add_argument('--uniprot_meta_tsv', metavar='PATH', default='', type=str)
    parser.add_argument('--existing_annotation_tsv', metavar='PATH', action='append', default=[])
    parser.add_argument('--outdir', metavar='PATH', default='', type=str)
    parser.add_argument('--audit_json', metavar='PATH', default='', type=str)
    parser.add_argument('--outfile', metavar='PATH', default='', type=str)
    args = parser.parse_args()

    if args.existing_annotation_tsv:
        if not args.uniprot_meta_tsv or not args.outdir or args.outfile or args.diamond_tsv or args.query_fasta or args.uniprot_fasta:
            parser.error('Saved-annotation GO refresh requires --uniprot_meta_tsv and --outdir only')
        refresh_go_annotations(args.existing_annotation_tsv, args.uniprot_meta_tsv, args.outdir, args.audit_json)
        return
    if not args.outfile or args.outdir or args.audit_json:
        parser.error('Search formatting requires --outfile; --outdir/--audit_json are for GO refresh')

    query_ids = read_fasta_ids(args.query_fasta)
    hits = load_best_diamond_hits(args.diamond_tsv)

    if len(query_ids) == 0:
        query_ids = hits['qseqid'].dropna().drop_duplicates().tolist()

    out = pandas.DataFrame(index=pandas.Index(query_ids, name='gene_id'))
    # No-hit searches still join descriptions/metadata through this key.
    out['sprot_best'] = ''

    if hits.shape[0] > 0:
        hit_map = hits.loc[:, ['qseqid', 'sseqid', 'coverage', 'pident', 'evalue']].set_index('qseqid')
        gene_index = out.index.to_series()
        out['sprot_best'] = gene_index.map(hit_map['sseqid'])
        out['sprot_alias'] = out['sprot_best']
        out['sprot_coverage'] = gene_index.map(hit_map['coverage']).map(lambda x: fmt_float(x, precision=2))
        out['sprot_identity'] = gene_index.map(hit_map['pident']).map(lambda x: fmt_float(x, precision=2))
        out['sprot_evalue'] = gene_index.map(hit_map['evalue']).map(fmt_evalue)

    accessions = set(out['sprot_best'].dropna())
    uniprot_desc = read_uniprot_descriptions(args.uniprot_fasta, accessions)
    if len(uniprot_desc) > 0:
        out['sprot_recname'] = out['sprot_best'].map(uniprot_desc).fillna('')

    metadata = read_uniprot_metadata(args.uniprot_meta_tsv, accessions)
    if metadata.shape[0] > 0:
        out = out.join(metadata.set_index('accession'), on='sprot_best')

    out = out.reset_index()
    for col in OUTPUT_COLUMNS:
        if col not in out.columns:
            out[col] = ''
        out[col] = out[col].fillna('')

    out = out.loc[:, OUTPUT_COLUMNS]
    out.to_csv(args.outfile, sep='\t', index=False)


if __name__ == '__main__':
    main()
