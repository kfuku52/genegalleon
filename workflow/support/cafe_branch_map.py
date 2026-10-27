#!/usr/bin/env python3
"""Map a CAFE branch between sampled species trees using descendant taxa.

Only a unique non-root clade is eligible. Added species are reported explicitly;
an added sister species can split one old edge into ambiguous candidates.
"""
import argparse
import csv
import hashlib
import io
import json
import re
from pathlib import Path

from Bio import Phylo


def file_sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for data in iter(lambda: handle.read(1024**2), b''):
            value.update(data)
    return value.hexdigest()


def first_newick(path):
    """Read the first CAFE tree without materializing every family tree."""
    parts = []
    nexus = False
    size = 0
    with Path(path).open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.upper().startswith('#NEXUS'):
                nexus = True
                continue
            if nexus and re.match(r'^TRANSLATE\b', stripped, re.I):
                raise ValueError('Translated Nexus labels require an explicit resolved species tree')
            if not parts:
                match = re.match(r'^TREE\s+.+?=\s*(.*)', stripped, re.I)
                if match:
                    line = match[1]
                elif not nexus and stripped.startswith('('):
                    line = stripped
                else:
                    continue
            parts.append(line)
            size += len(line)
            if size > 1024**2:
                raise ValueError('First CAFE tree exceeds the input bound')
            if ';' in line:
                return ''.join(parts).split(';', 1)[0] + ';'
    raise ValueError('No complete CAFE species tree found')


def branch_inventory(path):
    tree = Phylo.read(io.StringIO(first_newick(path)), 'newick')
    descendants = {}
    rows = []
    species = set()
    branch_ids = set()
    for node in tree.find_clades(order='postorder'):
        label = (node.name or '').replace('*', '')
        match = re.fullmatch(r'(.*?)(<\d+>)(?:_\d+)?', label)
        if match is None or match[2] in branch_ids:
            raise ValueError('Missing or duplicate native CAFE branch ID: ' + label)
        taxon, branch_id = match[1], match[2]
        branch_ids.add(branch_id)
        if node.is_terminal():
            if not taxon or taxon in species:
                raise ValueError('Missing or duplicate species identity: ' + taxon)
            species.add(taxon)
            tips = {taxon}
        else:
            if taxon:
                raise ValueError('Unexpected taxon prefix on an internal CAFE branch')
            tips = set().union(*(descendants[child] for child in node.clades))
        descendants[node] = tips
        rows.append({'branch_id': branch_id, 'native_key': taxon + branch_id,
                     'is_root': node is tree.root, 'descendants': sorted(tips)})
    return species, rows


def map_branch(old_path, new_path, target):
    old_species, old_rows = branch_inventory(old_path)
    new_species, new_rows = branch_inventory(new_path)
    selected = [row for row in old_rows if target in (row['branch_id'], row['native_key'])]
    if len(selected) != 1 or selected[0]['is_root']:
        raise ValueError('Old target must identify exactly one non-root CAFE branch')
    old = selected[0]
    wanted = set(old['descendants'])
    common = old_species & new_species
    missing = wanted - new_species
    candidates = []
    if not missing:
        for row in new_rows:
            tips = set(row['descendants'])
            if not row['is_root'] and tips & common == wanted:
                candidates.append({**row, 'added_descendants': sorted(tips - old_species)})
    candidates.sort(key=lambda row: row['native_key'])
    eligible = len(candidates) == 1 and not missing
    state = ('missing_target_species' if missing else 'unmapped_clade' if not candidates
             else 'ambiguous_edge_split' if len(candidates) > 1 else
             'exact_clade' if not candidates[0]['added_descendants'] else 'unique_clade_with_added_species')
    return {'schema_version': 1, 'state': state, 'eligible': eligible,
            'old_tree_sha256': file_sha(old_path), 'new_tree_sha256': file_sha(new_path),
            'old_target': old, 'missing_target_species': sorted(missing),
            'added_species': sorted(new_species - old_species),
            'removed_species': sorted(old_species - new_species),
            'common_species': sorted(common), 'candidates': candidates,
            'new_target': candidates[0]['native_key'] if eligible else None}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--old-tree', type=Path, required=True)
    parser.add_argument('--new-tree', type=Path, required=True)
    parser.add_argument('--old-target', required=True)
    parser.add_argument('--new-branch-table', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    result = map_branch(args.old_tree, args.new_tree, args.old_target)
    if args.new_branch_table is not None:
        with args.new_branch_table.open() as handle:
            columns = next(csv.reader(handle, delimiter='\t'))
        if any(row['native_key'] not in columns for row in result['candidates']):
            raise ValueError('Mapped CAFE keys do not match the new native branch table')
        result['new_branch_table_sha256'] = file_sha(args.new_branch_table)
    text = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output is not None:
        with args.output.open('x') as handle:
            handle.write(text)
    print(text, end='')
    return 0 if result['eligible'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
