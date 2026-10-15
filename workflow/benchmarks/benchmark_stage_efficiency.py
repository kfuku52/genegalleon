#!/usr/bin/env python3
"""Compare stage preparation, summary, ordered queues and multi-format plotting."""
import argparse
import contextlib
import hashlib
import io
import json
import os
import platform
import re
import resource
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
CASES = ('preparation', 'orthologs', 'closest', 'ordered-queue', 'plots')


def worker(args):
    sys.path.insert(0, str(ROOT))
    sys.path.insert(0, str(args.support))
    counts = Counter()
    fingerprint = hashlib.sha256()
    with tempfile.TemporaryDirectory(prefix='gg-stage-efficiency-') as temporary:
        root = Path(temporary)
        if args.case == 'preparation':
            import pairwise_synteny as module

            sources = []
            for index in range(4):
                species = f'Species_{index}'
                fasta, gff = root / (species + '.fa'), root / (species + '.gff3')
                with fasta.open('w') as fa, gff.open('w') as annotation:
                    annotation.write('##gff-version 3\n')
                    for number in range(args.genes):
                        name, start = f'g{number:05d}', number * 1000 + 1
                        fa.write(f'>{name}\n' + 'ATG' + 'GCT' * 298 + 'TAA\n')
                        annotation.write(f'chr1\ts\tgene\t{start}\t{start+899}\t.\t+\t.\tID=locus_{name}\n')
                        annotation.write(f'chr1\ts\tmRNA\t{start}\t{start+899}\t.\t+\t.\tID={name};Parent=locus_{name}\n')
                sources.append(dict(species=species, fasta=str(fasta), gff=str(gff), mode='cds',
                                    genetic_code=1, feature='mRNA', attribute='ID'))
            original_records, original_annotation = module.fasta_records, module.annotation_to_genes
            def records(*a, **kw):
                counts['fasta_parses'] += 1
                yield from original_records(*a, **kw)
            def annotation(*a, **kw):
                counts['gff_parses'] += 1
                return original_annotation(*a, **kw)
            module.fasta_records, module.annotation_to_genes = records, annotation
            def native_boundary(command, directory, _commands, label):
                if label == 'mcscan':
                    identifiers = [[line.split('\t')[3] for line in (directory / (side + '.bed')).read_text().splitlines()]
                                   for side in ('target', 'query')]
                    anchors = '# block\n' + ''.join(f'{a}\t{b}\t100\n' for a, b in zip(*identifiers, strict=True))
                    (directory / 'target.query.lifted.anchors').write_text(anchors)
                else:
                    (directory / 'target.query.screened.simple').write_text('# bounded native substitute\n')
            module.run_command = native_boundary
            plan = dict(pairs=[dict(analysis_id=f'pair{i}', target=sources[i % 4], query=sources[(i+1) % 4])
                               for i in range(args.pairs)], tools={}, parameters=dict(minimum_mapping_fraction=1,
                               cscore=.7, min_anchors=4, distance=20))
            output = root / 'analysis'
            started = time.perf_counter()
            module.analyze(plan, output, args.workers)
            seconds = time.perf_counter() - started
            for path in sorted(output.rglob('*')):
                if path.is_file():
                    fingerprint.update(str(path.relative_to(output)).encode())
                    fingerprint.update(path.read_bytes().replace(str(root).encode(), b'<workspace>'))
        elif args.case in ('orthologs', 'closest'):
            import query_gene_orthologs as module
            from Bio import Phylo

            from workflow.tests.test_query_gene_orthologs import weak_duplication_fixture

            seed = root / 'seed'
            seed.mkdir()
            fixture = weak_duplication_fixture.__wrapped__(seed)
            query, families, output = root / 'queries', root / 'families', root / 'output'
            query.mkdir()
            for subdir in ('stat_branch', 'cds_fasta'):
                (families / subdir).mkdir(parents=True)
            for index in range(args.families):
                family = f'FAM{index:05d}'
                shutil.copyfile(fixture['query_dir'] / 'FAM', query / family)
                for subdir, suffix in (('stat_branch', '_stat.branch.tsv'), ('cds_fasta', '_cds.fasta')):
                    shutil.copyfile(fixture['output_root'] / subdir / ('FAM' + suffix), families / subdir / (family + suffix))
            for name in ('read_stat_branch', 'read_family_cds_fasta_ids', 'build_tree_index'):
                # New public wrappers can reuse an input context. Count actual
                # reads/validation, matching the baseline's unwrapped bodies.
                target = '_' + name if hasattr(module, '_' + name) else name
                original = getattr(module, target)
                def observed(*a, _name=name, _original=original, **kw):
                    counts[_name] += 1
                    return _original(*a, **kw)
                setattr(module, target, observed)
            original_phylo_read = Phylo.read
            def tree(*a, **kw):
                counts['species_tree_parses'] += 1
                return original_phylo_read(*a, **kw)
            Phylo.read = tree
            options = dict(basis='query_gene' if args.case == 'closest' else 'reference_species',
                dir_gene_family=str(families), dir_query_gene=str(query), family_file='',
                reference_species='Reference_species', dup_conf_score_threshold=.05,
                **{'out_' + name: str(output / (name + '.tsv')) for name in
                   ('columns', 'glyphs', 'tree', 'synteny', 'ufboot', 'dup_conf', 'query_map', 'selection')})
            if args.case == 'closest':
                options.update(query_selection='closest', species_tree=str(fixture['species_tree']),
                               target_species='Reference_species', selection_species='Reference_species')
            started = time.perf_counter()
            with contextlib.redirect_stdout(io.StringIO()):
                module.run(SimpleNamespace(**options))
            seconds = time.perf_counter() - started
            for path in sorted(output.iterdir()):
                fingerprint.update(path.name.encode() + path.read_bytes())
        elif args.case == 'ordered-queue':
            import parallel_io

            def arguments():
                for index in range(args.tasks):
                    counts['submitted'] += 1
                    yield (index,)
            def task(index):
                if index == 0:
                    time.sleep(.1)
                return index, b'x' * args.payload
            started = time.perf_counter()
            with ThreadPoolExecutor(max_workers=args.workers) as pool:
                if hasattr(parallel_io, 'ordered_bounded_results'):
                    results = parallel_io.ordered_bounded_results(pool, task, arguments(), args.workers * 8)
                else:
                    results = pool.map(task, (item[0] for item in arguments()))
                with contextlib.closing(results):
                    for expected, (index, payload) in enumerate(results):
                        if expected == 0:
                            counts['submitted_before_first_result'] = counts['submitted']
                        assert index == expected
                        fingerprint.update(index.to_bytes(8, 'big') + payload)
                        counts['results'] += 1
            seconds = time.perf_counter() - started
        else:
            directory = root / 'plot'
            directory.mkdir()
            for side, prefix in (('target', 't'), ('query', 'q')):
                (directory / (side + '.bed')).write_text(''.join(f'chr1\t{i*10}\t{i*10+5}\t{prefix}{i}\n' for i in range(6)))
            (directory / 'seqids').write_text('chr1\nchr1\n')
            (directory / 'layout').write_text('0.7,0.12,0.92,0,,Target species,top,target.bed,top\n'
                '0.3,0.12,0.92,0,,Query species,bottom,query.bed,bottom\ne,0,1,colored.simple\n')
            (directory / 'colored.simple').write_text('#88afc4*t0 t2 q1 q4 5 +\n')
            (directory / 'karyotype_colors.json').write_text(json.dumps({'chromosomes': {
                side: {'chr1': '#88afc4'} for side in ('target', 'query')}}))
            config = root / 'matplotlib'
            config.mkdir()
            (config / 'matplotlibrc').write_text('svg.hashsalt: gg-stage-efficiency\n')
            env = dict(os.environ, MPLCONFIGDIR=str(config), SOURCE_DATE_EPOCH='1780000000')
            script = args.support / 'pairwise_synteny_karyotype.py'
            command = [sys.executable, str(script), '--directory', str(directory),
                       '--target-species', 'Target_species', '--query-species', 'Query_species']
            formats = ('pdf', 'png', 'svg')
            batched = '--formats' in script.read_text()
            commands = [command + ['--formats', ','.join(formats)]] if batched else [command + ['--format', fmt] for fmt in formats]
            started = time.perf_counter()
            for argv in commands:
                subprocess.run(argv, check=True, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                counts['plot_processes'] += 1
            seconds = time.perf_counter() - started
            for fmt in formats:
                data = (directory / ('karyotype.' + fmt)).read_bytes()
                # Dates and SVG element IDs are non-scientific metadata. The
                # fixed epoch/hash salt above normally makes this a no-op.
                if fmt == 'pdf':
                    data = re.sub(rb'/CreationDate \([^)]*\)', b'/CreationDate (<fixed>)', data)
                elif fmt == 'svg':
                    data = re.sub(rb'<dc:date>[^<]*</dc:date>', b'<dc:date>fixed</dc:date>', data)
                fingerprint.update(fmt.encode() + data)
            fingerprint.update((directory / 'karyotype_style.json').read_bytes())
        rss = max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
        return dict(seconds=seconds, peak_rss_kib=rss, counts=dict(counts), logical_sha256=fingerprint.hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-support', type=Path)
    parser.add_argument('--baseline-only', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cases', choices=CASES, nargs='+', default=list(CASES))
    for name, default in (('genes', 4000), ('pairs', 14), ('families', 128), ('tasks', 20000),
                          ('payload', 4096), ('workers', 4), ('repeats', 3)):
        parser.add_argument('--' + name, type=int, default=default)
    parser.add_argument('--support', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--case', choices=CASES, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if min(args.genes, args.pairs, args.families, args.tasks, args.payload, args.workers, args.repeats) < 1:
        parser.error('sizes, workers and repeats must be positive')
    if args.support:
        print(json.dumps(worker(args)))
        return
    if args.baseline_only and not args.baseline_support:
        parser.error('--baseline-only requires --baseline-support')
    supports = {} if args.baseline_only else {'current': ROOT / 'workflow/support'}
    if args.baseline_support:
        supports = {'baseline': args.baseline_support, **supports}
    result = dict(environment=dict(platform=platform.platform(), python=platform.python_version()),
        configuration={key: getattr(args, key) for key in ('genes', 'pairs', 'families', 'tasks', 'payload', 'workers', 'repeats')},
        warmups=1, workloads={}, source_sha256={label: {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
            for name in ('pairwise_synteny.py', 'query_gene_orthologs.py', 'parallel_io.py', 'pairwise_synteny_ds.py',
                         'wgd_ssd.py', 'pairwise_synteny_karyotype.py')} for label, root in supports.items()})
    for case in args.cases:
        rows = {label: [] for label in supports}
        for trial in range(args.repeats + 1):
            items = list(supports.items())
            if trial % 2:
                items.reverse()
            for label, root in items:
                command = [sys.executable, str(Path(__file__).resolve()), '--support', str(root), '--case', case,
                           '--output', str(args.output)]
                for name in ('genes', 'pairs', 'families', 'tasks', 'payload', 'workers'):
                    command += ['--' + name, str(getattr(args, name))]
                completed = subprocess.run(command, text=True, capture_output=True, check=True)
                row = json.loads(completed.stdout)
                if trial:
                    rows[label].append(row)
        if len({row['logical_sha256'] for values in rows.values() for row in values}) != 1:
            raise ValueError('Different outputs: ' + case)
        result['workloads'][case] = {label: dict(trials=values, median_seconds=statistics.median(row['seconds'] for row in values),
            median_peak_rss_kib=statistics.median(row['peak_rss_kib'] for row in values)) for label, values in rows.items()}
        print(case, {label: result['workloads'][case][label]['median_seconds'] for label in supports}, flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
