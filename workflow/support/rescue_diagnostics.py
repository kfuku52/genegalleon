"""Store rescue summaries and reconstruct detailed diagnostics from model evidence."""
import csv
import gzip
from collections import Counter
from pathlib import Path

AUDIT = ('candidate', 'donor_gene', 'status', 'reasons', 'coverage', 'identity', 'model_id')
QUALITY = ('candidate', 'model_id', 'status', 'start_codon', 'donor_n_terminus_aligned',
           'donor_c_terminus_aligned', 'donor_species', 'flags', 'donor_aligned_query_fraction',
           'donor_internal_unaligned_query_fraction')
SUMMARY = 'diagnostics.json'
UNRESOLVED = 'unresolved_queries.tsv.gz'
COUNTERS = ('status_counts', 'reason_counts', 'quality_flag_counts')


def audit_row(model):
    return (model['query'], model['evidence']['query'], model['status'], ','.join(model['problems']),
            model['coverage'], model['identity'], model.get('model_id', ''))


def quality_row(model):
    evidence = model['quality_evidence']
    terminal, alignment = evidence['terminal_alignment'], evidence['query_alignment']
    return (model['query'], model.get('model_id', ''), model['status'], evidence['start_codon'],
            terminal['n_aligned'], terminal['c_aligned'], ','.join(evidence['donor_species']),
            ','.join(evidence['flags']), alignment['aligned_query_fraction'],
            alignment['internal_unaligned_query_fraction'])


def new_counts():
    return {key: Counter() for key in COUNTERS}


def update(counts, model):
    counts['status_counts'][model['status']] += 1
    counts['reason_counts'].update(model['problems'])
    counts['quality_flag_counts'].update(model['quality_evidence']['flags'])


def summarize(models):
    counts = new_counts()
    for model in models:
        update(counts, model)
    return counts


def replace_counts(summary, before, after):
    updated = dict(summary)
    for key in COUNTERS:
        counts = Counter(summary[key])
        counts.subtract(before[key])
        if any(value < 0 for value in counts.values()):
            raise ValueError('Stored rescue diagnostics disagree with their model records')
        counts.update(after[key])
        updated[key] = +counts
    return updated


def text_open(path, mode):
    path = Path(path)
    if path.suffix == '.gz':
        return gzip.open(path, mode + 't', newline='', encoding='utf-8', compresslevel=1)
    return path.open(mode, newline='', encoding='utf-8')


def table_writer(stream):
    return csv.writer(stream, delimiter='\t', lineterminator='\n')


def write_table(path, header, rows):
    with text_open(path, 'w') as handle:
        writer = table_writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def unresolved_row(candidate, query):
    return candidate, query, 'unresolved', 'no_alignment', '', '', ''


def unresolved_rows(directory):
    path = directory / UNRESOLVED
    if path.is_file():
        with text_open(path, 'r') as handle:
            rows = csv.reader(handle, delimiter='\t')
            next(rows)
            yield from (unresolved_row(candidate, query) for candidate, query in rows)


def detail_name(stem, mode):
    return stem + ('.gz' if mode == 'full' else '')


def legacy_table(directory, stem):
    for name in (stem, stem + '.gz'):
        if (directory / name).is_file():
            return directory / name
    return None


def export_tables(source, destination, models, candidates=()):
    with (destination / 'audit.tsv').open('x', newline='') as audit, (
            destination / 'quality_flags.tsv').open('x', newline='') as quality:
        audit, quality = table_writer(audit), table_writer(quality)
        audit.writerow(AUDIT)
        quality.writerow(QUALITY)
        for model in models:
            if not candidates or model['query'] in candidates:
                audit.writerow(audit_row(model))
                quality.writerow(quality_row(model))
        audit.writerows(row for row in unresolved_rows(source) if not candidates or row[0] in candidates)
