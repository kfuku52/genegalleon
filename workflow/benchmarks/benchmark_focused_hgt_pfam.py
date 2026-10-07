#!/usr/bin/env python3
"""Measure focused HGT Pfam filtering, including saved-hit reads and source checks."""

import argparse
import csv
import hashlib
import json
import os
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def fingerprint(records):
    digest = hashlib.sha256()
    for group in records:
        digest.update(b'group\n')
        for row in group:
            digest.update(json.dumps(row, sort_keys=True, separators=(',', ':')).encode())
            digest.update(b'\n')
    return digest.hexdigest()


def fixture(root, families, events_per_family):
    events, links = [], []
    rps = root/'rpsblast'
    rps.mkdir()
    prefix = 'host_scaffold_background_class_'
    for index in range(families):
        family = f'OG{index:07d}'
        # Repeated local gene names deliberately test family-scoped lookups.
        queries = []
        for side in ('donor', 'recipient'):
            for copy in range(3):
                gene = f'{side}_{copy}'
                row = dict(qacc=gene, sacc='', qlen='100', stitle='', qstart='', qend='', evalue='')
                if index % 7 != 0:
                    row.update(sacc='model', stitle='pfam01053, Enzyme, Saved domain',
                               qstart='1', qend='50' if copy < 2 else '49', evalue='1e-8')
                queries.append(row)
        if index % 11 != 0:
            with (rps/f'{family}_rpsblast.tsv').open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(queries[0]), delimiter='\t', lineterminator='\n')
                writer.writeheader()
                writer.writerows(queries)
        for position in range(events_per_family):
            event = dict(event_id=f'{family}:3:{position + 1}', orthogroup=family, branch_id='3',
                         node_name='n3', event_index=str(position + 1), generax_transfer='Y@D@A',
                         generax_donor_node='D', generax_recipient_node='A')
            events.append(event)
            for side in ('donor', 'recipient'):
                for copy in range(3):
                    link = dict(event, side=side, gene_id=f'{side}_{copy}', eligible_for_context='True',
                                lineage_status='retained', host_scaffold_status='measured', host_scaffold_id='s1')
                    link.update({prefix + name: value for name, value in
                                 [('total_count', '20'), ('compatible_count', '9'), ('incompatible_count', '1'),
                                  ('unresolved_count', '10'), ('classified_fraction', '0.5'), ('compatible_fraction', '0.9')]})
                    links.append(link)
                # These descendants must not rescue a rejected pair.
                links.append(dict(link, gene_id=f'{side}_excluded', eligible_for_context='False',
                                  lineage_status='transferred_out'))
    return events, links


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def worker(args):
    sys.path.insert(0, str(args.support_root))
    from focus_hgt_pfam import EVENT_FIELDS, filter_events

    csv.field_size_limit(100_000_000)
    if args.events:
        events = [{k: v for k, v in row.items() if k not in EVENT_FIELDS} for row in read(args.events)]
        links, root = read(args.links), args.family_root
    else:
        events, links = fixture(args.worker, args.families, args.events_per_family)
        root = args.worker
    inputs = fingerprint((events, links))
    started = time.perf_counter()
    selected, audit, pairs, genes, sources = filter_events(events, links, root)
    seconds = time.perf_counter() - started
    assert fingerprint((events, links)) == inputs, 'Benchmark input records were mutated'
    print(json.dumps(dict(seconds=seconds, peak_rss_kib_linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                          input_sha256=inputs, output_sha256=fingerprint((selected, audit, pairs, genes)),
                          source_sha256=hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest(),
                          event_count=len(events), family_count=len({row['orthogroup'] for row in events}),
                          link_count=len(links), passing_event_count=len(selected), pair_count=len(pairs))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1]/'support')
    parser.add_argument('--families', type=int, default=2048)
    parser.add_argument('--events-per-family', type=int, default=4)
    parser.add_argument('--events', type=Path, help='Existing event TSV; saved Pfam decision columns are removed')
    parser.add_argument('--links', type=Path)
    parser.add_argument('--family-root', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if not 1 <= args.families <= 8192 or not 1 <= args.events_per_family <= 16:
        parser.error('Use 1..8192 families and 1..16 events per family')
    if any((args.events, args.links, args.family_root)) and not all((args.events, args.links, args.family_root)):
        parser.error('--events, --links and --family-root must be supplied together')
    if args.worker:
        worker(args)
        return
    if not args.output:
        parser.error('--output is required')
    results = dict(support_root=str(args.support_root), python=sys.version, platform=platform.platform(),
                   command=[sys.executable, *sys.argv],
                   workload='existing_saved_results' if args.events else 'synthetic_multiple_families',
                   configuration=dict(families=args.families, events_per_family=args.events_per_family,
                                      min_shared_pfam_coverage=0.5, allow_both_no_pfam=False),
                   inputs={name: str(path.resolve()) for name, path in
                           [('events', args.events), ('links', args.links), ('family_root', args.family_root)] if path},
                   samples=[])
    environment = dict(os.environ, PYTHONHASHSEED='0')
    with tempfile.TemporaryDirectory(prefix='gg-focused-hgt-pfam-benchmark-') as temporary:
        for trial in range(4):
            root = Path(temporary)/str(trial)
            root.mkdir()
            command = [sys.executable, str(Path(__file__).resolve()), '--worker', str(root),
                       '--support-root', str(args.support_root), '--families', str(args.families),
                       '--events-per-family', str(args.events_per_family)]
            if args.events:
                command += ['--events', str(args.events.resolve()), '--links', str(args.links.resolve()),
                            '--family-root', str(args.family_root.resolve())]
            sample = json.loads(subprocess.check_output(command, text=True, env=environment))
            if trial:
                results['samples'].append(sample)
            print('warmup' if trial == 0 else trial, sample['seconds'], flush=True)
    for field in ('input_sha256', 'output_sha256', 'source_sha256'):
        assert len({sample[field] for sample in results['samples']}) == 1, field
    for metric in ('seconds', 'peak_rss_kib_linux'):
        results['median_' + metric] = statistics.median(sample[metric] for sample in results['samples'])
    args.output.write_text(json.dumps(results, indent=2) + '\n')
    print('medians', results['median_seconds'], results['median_peak_rss_kib_linux'], flush=True)


if __name__ == '__main__':
    main()
