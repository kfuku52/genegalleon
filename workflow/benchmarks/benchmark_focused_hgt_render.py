#!/usr/bin/env python3
"""Compare real focused-HGT tree/context exports using existing saved inputs.

Run inside one GeneGalleon runtime with frozen baseline and candidate support
directories. The exporter uses the shared gg_gene_evolution tree renderer.
No sequence searches, phylogenetic inference, or candidate filtering is run.
"""

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import resource
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

PDF_DATES = re.compile(rb'/(?:CreationDate|ModDate)\s*\(D:[^)]*\)')


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def normalize_settings(settings):
    """Only generated input paths change; retain every scientific argument."""
    if isinstance(settings, dict):
        return {key: normalize_settings(value) for key, value in settings.items()}
    if isinstance(settings, list):
        return [normalize_settings(value) for value in settings]
    if isinstance(settings, str):
        return re.sub(r'/[^,|\s]*/family_inputs(?=/|$)', '{family_inputs}', settings)
    return settings


def pdf_identity(path, pages):
    from pypdf import PdfReader

    raw = path.read_bytes()
    measured = len(PdfReader(path).pages)
    if not raw.startswith(b'%PDF') or measured != pages:
        raise ValueError(f'Expected {pages} PDF pages: {path}')
    # Date removal is intentionally narrow. Images, fonts, content, IDs and
    # structural bytes must remain identical for a successful comparison.
    stable = PDF_DATES.sub(lambda match: b'_' * len(match.group()), raw)
    return dict(pages=measured, size_bytes=len(raw), without_dates_sha256=digest(stable))


def worker(args):
    sys.path.insert(0, str(args.support_root))
    import pypdf.filters
    from focus_hgt_gene_trees import export_gene_trees

    csv.field_size_limit(100_000_000)
    # Keep decompression bounded while accommodating large saved gene trees.
    pypdf.filters.ZLIB_MAX_OUTPUT_LENGTH = 268435456
    events, links, request_sources = [], [], {}
    for family in args.families:
        directory = args.event_plan / family
        for name, destination in [('events.tsv', events), ('links.tsv', links)]:
            path = directory / name
            request_sources[f'{family}/{name}'] = digest(path.read_bytes())
            destination.extend(read_tsv(path))
    args.worker.mkdir(parents=True)
    native = args.worker / 'native'
    native.mkdir()
    run = subprocess.run
    renderer_calls = []

    def capture(command, *positional, **kwargs):
        result = run(command, *positional, **kwargs)
        if command[0] == 'Rscript':
            renderer_calls.append(Path(command[1]).name)
            if result.returncode:
                return result
            if Path(command[1]).name == 'stat_branch2tree_plot.r':
                table = next(value.split('=', 1)[1] for value in command if value.startswith('--stat_branch='))
                family = Path(table).name.split('_focused_stat.branch.tsv')[0]
                shutil.copyfile(Path(kwargs['cwd']) / 'stat_branch2tree_plot.pdf', native / f'{family}.pdf')
            elif Path(command[1]).name == 'tree_plot_batch.r':
                for job in json.loads(Path(command[2]).read_text()):
                    shutil.copyfile(job['output'], native / f"{job['id']}.pdf")
        return result

    subprocess.run = capture
    output = args.worker / 'export'
    os.environ['MPLCONFIGDIR'] = str(args.output / 'matplotlib-cache')
    options = dict(gff_root=args.inputs / 'gff_info',
                   context_annotations=args.inputs / 'context_gene_annotations.tsv',
                   mmseqs2_taxonomy_dir=args.inputs / 'species_cds_mmseqs2taxonomy',
                   scaffold_taxonomy_dir=args.inputs / 'species_scaffold_taxonomy',
                   taxonomy_dbfile=args.inputs / 'context_taxonomy.sqlite', minimum_ufboot=None)
    started = time.perf_counter()
    result = export_gene_trees(output, events, links, args.inputs / 'families', **options)
    seconds = time.perf_counter() - started
    python_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    child_rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    subprocess.run = run
    # Verification is outside the timed interval. Keep all PDFs and audits.
    if result['rendered_family_count'] != len(args.families):
        raise ValueError('Some requested families were not rendered')
    identity = dict(result=result, request_sources=request_sources, artifacts={})
    for path in sorted(output.rglob('*')):
        if not path.is_file():
            continue
        relative = str(path.relative_to(output))
        if path.suffix == '.pdf':
            identity['artifacts'][relative] = pdf_identity(path, 2)
        elif path.name == 'renderer_settings.json':
            normalized = normalize_settings(json.loads(path.read_text()))
            identity['artifacts'][relative] = digest(json.dumps(normalized, sort_keys=True).encode())
        else:
            identity['artifacts'][relative] = digest(path.read_bytes())
    for family in args.families:
        identity['artifacts'][f'native/{family}.pdf'] = pdf_identity(native / f'{family}.pdf', 1)
    identity_sha = digest(json.dumps(identity, sort_keys=True).encode())
    sample = dict(seconds=seconds, python_peak_rss_kib_linux=python_rss,
                  child_peak_rss_kib_linux=child_rss,
                  identity_sha256=identity_sha, renderer_calls=renderer_calls,
                  family_count=len(args.families), event_count=len(events), link_count=len(links))
    (args.worker / 'identity.json').write_text(json.dumps(identity, indent=2) + '\n')
    (args.worker / 'sample.json').write_text(json.dumps(sample, indent=2) + '\n')
    print(json.dumps(sample), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--event-plan', type=Path, required=True, help='Directory containing FAMILY/events.tsv and links.tsv')
    parser.add_argument('--families', nargs='+', required=True, help='Explicit saved family identifiers to compare')
    parser.add_argument('--baseline-root', type=Path)
    parser.add_argument('--candidate-root', type=Path)
    parser.add_argument('--trials', type=int, default=3)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--support-root', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not 1 <= args.trials <= 10:
        parser.error('Use a trial count from 1..10')
    for name in ('inputs', 'event_plan', 'output', 'support_root', 'baseline_root', 'candidate_root', 'worker'):
        value = getattr(args, name)
        if value:
            setattr(args, name, value.resolve())
    if args.worker:
        if not args.support_root:
            parser.error('--worker needs --support-root')
        worker(args)
        return
    if not args.baseline_root or not args.candidate_root:
        parser.error('--baseline-root and --candidate-root are required for a comparison')
    args.output.mkdir(parents=True)
    report = dict(command=[sys.executable, *sys.argv], platform=platform.platform(), python=sys.version,
                  families=args.families, inputs=str(args.inputs),
                  baseline_root=str(args.baseline_root), candidate_root=str(args.candidate_root),
                  rss_definition='Maximum resident set of Python process and largest child, separately; not their sum',
                  samples=[], warmups=[])
    environment = dict(os.environ, PYTHONHASHSEED='0', PYTHONDONTWRITEBYTECODE='1')
    for trial in range(args.trials + 1):
        order = ('baseline', 'candidate') if trial % 2 else ('candidate', 'baseline')
        for variant in order:
            directory = args.output / f'{variant}-{trial}'
            root = getattr(args, variant + '_root')
            command = [sys.executable, str(Path(__file__).resolve()), '--worker', str(directory),
                       '--support-root', str(root), '--inputs', str(args.inputs),
                       '--event-plan', str(args.event_plan), '--families', *args.families,
                       '--output', str(args.output)]
            sample = json.loads(subprocess.check_output(command, text=True, env=environment))
            sample.update(variant=variant, trial=trial)
            report['warmups' if trial == 0 else 'samples'].append(sample)
            print(variant, 'warmup' if trial == 0 else trial, sample['seconds'], flush=True)
            (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    for variant in ('baseline', 'candidate'):
        samples = [sample for sample in report['samples'] if sample['variant'] == variant]
        report[variant] = {f'median_{metric}': statistics.median(sample[metric] for sample in samples)
                           for metric in ('seconds', 'python_peak_rss_kib_linux', 'child_peak_rss_kib_linux')}
    report['equivalent'] = len({sample['identity_sha256'] for sample in report['samples'] + report['warmups']}) == 1
    report['speedup'] = report['baseline']['median_seconds'] / report['candidate']['median_seconds']
    (args.output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: report[key] for key in ('baseline', 'candidate', 'speedup', 'equivalent')}), flush=True)
    if not report['equivalent']:
        raise ValueError('Exported PDFs, audits, settings or consumed-input hashes differ')


if __name__ == '__main__':
    main()
