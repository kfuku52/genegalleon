#!/usr/bin/env python3
"""Validate species-tree constraints and select rows without realigning columns."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

from Bio import Phylo, SeqIO
from Bio.Phylo.NewickIO import NewickError


def read_tree(path):
    tree = Phylo.read(path, 'newick')
    labels = [tip.name for tip in tree.get_terminals()]
    if not labels or any(not label for label in labels):
        raise ValueError('Tree tips must have nonempty names')
    repeated = sorted(label for label, n in Counter(labels).items() if n > 1)
    if repeated:
        raise ValueError(f'Duplicate tree tips: {repeated}')
    return tree, set(labels)


def split_set(tree, taxa):
    """Unrooted, nontrivial splits, optionally restricted to a partial backbone."""
    result = set()
    for clade in tree.find_clades():
        left = frozenset(tip.name for tip in clade.get_terminals()) & taxa
        right = frozenset(taxa - left)
        if min(len(left), len(right)) >= 2:
            result.add(min((left, right), key=lambda x: (len(x), sorted(x))))
    return result


def alignment(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        records = list(SeqIO.parse(handle, 'fasta'))
    labels = [record.id for record in records]
    if not records or len(labels) != len(set(labels)):
        raise ValueError('Alignment must have unique, nonempty sequence IDs')
    lengths = {len(record.seq) for record in records}
    if len(lengths) != 1 or not next(iter(lengths)):
        raise ValueError('Alignment rows must have the same positive length')
    return records


def selected_labels(records, taxa_list):
    available = {record.id for record in records}
    if not taxa_list:
        return available
    labels = Path(taxa_list).read_text().splitlines()
    if any(not label or label != label.strip() or any(c.isspace() for c in label) for label in labels):
        raise ValueError('Taxa list must contain one exact, nonempty sequence ID per line')
    if len(labels) != len(set(labels)) or len(labels) < 4:
        raise ValueError('Taxa list must contain at least four unique IDs')
    unknown = set(labels) - available
    if unknown:
        raise ValueError(f'Taxa list IDs absent from alignment: {sorted(unknown)}')
    return set(labels)


def constraint_info(path, taxa, scope):
    tree, labels = read_tree(path)
    if labels - taxa:
        raise ValueError(f'Constraint tips absent from selected alignment: {sorted(labels - taxa)}')
    if scope == 'complete' and labels != taxa:
        raise ValueError(f'Constraint omits alignment tips: {sorted(taxa - labels)}')
    splits = split_set(tree, labels)
    if not splits:
        raise ValueError('Constraint must specify at least one nontrivial unrooted split')
    return labels, splits


def run(args):
    records = alignment(Path(args.alignment))
    taxa = selected_labels(records, args.taxa_list)
    labels, splits = (set(), set())
    if args.constraint:
        labels, splits = constraint_info(args.constraint, taxa, args.scope)
    report = {'schema_version': 1, 'alignment_taxa': len(records), 'selected_taxa': len(taxa),
              'alignment_columns': len(records[0].seq), 'constraint_scope': args.scope,
              'constraint_taxa': sorted(labels),
              'assumed_unrooted_splits': sorted([sorted(s) for s in splits]),
              'root_inferred_from_constraint': False,
              'support_interpretation': 'Constrained splits are prior assumptions, not independent support.'}
    for name, path in [('alignment', args.alignment), ('constraint', args.constraint), ('taxa_list', args.taxa_list)]:
        if path:
            digest = hashlib.sha256()
            with open(path, 'rb') as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b''):
                    digest.update(chunk)
            report[name + '_sha256'] = digest.hexdigest()
    if args.tree:
        tree, observed = read_tree(args.tree)
        if observed != taxa:
            raise ValueError(f'Result tree taxa differ: missing={sorted(taxa-observed)}, extra={sorted(observed-taxa)}')
        missing = splits - split_set(tree, labels)
        if missing:
            raise ValueError(f'Result tree violates {len(missing)} required unrooted constraint splits')
        report['result_constraint_verified'] = True
    if args.output:
        # Validate everything before opening a destination. Never rewrite a source.
        target = Path(args.output)
        for source in (args.alignment, args.constraint, args.taxa_list, args.tree):
            if source and (target.resolve() == Path(source).resolve() or
                           (target.exists() and target.samefile(source))):
                raise ValueError('Output overlaps an input')
        with target.open('w') as handle:
            SeqIO.write((record for record in records if record.id in taxa), handle, 'fasta')
    print(json.dumps(report, sort_keys=True))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--alignment', required=True)
    p.add_argument('--taxa-list', default='')
    p.add_argument('--constraint', default='')
    p.add_argument('--scope', choices=('complete', 'partial'), default='complete')
    p.add_argument('--output')
    p.add_argument('--tree')
    args = p.parse_args()
    try:
        run(args)
    except (ValueError, OSError, NewickError) as error:
        p.exit(2, f'Species-tree constraint error: {error}\n')


if __name__ == '__main__':
    main()
