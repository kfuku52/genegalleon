#!/usr/bin/env python3
"""Compare real helpers against a saved support directory, with exact outputs.

Use the GeneGalleon runtime. One warmup and alternating fresh-process trials
exclude fixture preparation/imports; cold table-index construction is included.
"""

import argparse
import contextlib
import hashlib
import importlib.metadata
import io
import json
import platform
import random
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
CASES = ('promoter', 'bootstrap', 'geometry', 'joint-layout', 'traits', 'cds', 'hgt', 'dated', 'kmer', 'rooting', 'r-annotation', 'r-busco')


def frame_payload(frame):
    return json.dumps([(str(c), str(t)) for c, t in frame.dtypes.items()]).encode() + frame.to_csv(sep='\t').encode()


def prepare(fixture):
    import numpy as np
    import pandas as pd

    matrix = np.arange(40000 * 80, dtype=float).reshape(40000, 80) / 7
    matrix[::97, 2] = np.nan
    frame = pd.DataFrame(matrix, columns=[f'tissue{i:02d}' for i in range(80)])
    frame.insert(0, 'gene_id', [f'g{i:05d}' for i in range(len(frame))])
    frame.to_csv(fixture / 'traits.tsv', sep='\t', index=False)
    annotation = fixture / 'annotation'
    annotation.mkdir()
    frame = pd.DataFrame({name: [f'{name}{i}' if name == 'gene_id' or i % 7 else '' for i in range(30000)]
        for name in ('gene_id', 'orthogroup', 'sprot_best', 'cdskit_localize_predicted_class', 'busco_id')})
    for col in range(75):
        frame[f'extra{col}'] = np.arange(len(frame), dtype=float) / 7
    frame.to_csv(annotation / 'Species_a.tsv', sep='\t', index=False)
    for name in ('Species_b', 'Species_c'):
        (annotation / (name + '.tsv')).hardlink_to(annotation / 'Species_a.tsv')


def worker(args):
    sys.path.insert(0, str(args.support.resolve()))
    counts = {}
    with tempfile.TemporaryDirectory(prefix='gg-efficiency-') as temporary:
        scratch = Path(temporary)
        if args.case == 'promoter':
            import get_promoter_fasta as module
            import pandas as pd

            genomes = scratch / 'genomes'
            genomes.mkdir()
            length = 2 * 1024 * 1024
            with (genomes / 'Species_a.fa').open('w') as stream:
                for chrom in range(3):
                    stream.write(f'>chr{chrom} reference chromosome\n')
                    sequence = ('ACGTacgtNN' * (length // 10 + 1))[:length]
                    stream.write('\n'.join(sequence[i:i + 60] for i in range(0, len(sequence), 60)) + '\n')
            rows = [dict(gene_id=f'Species_a_g{i:03d}', chromosome=f'chr{i % 3}',
                         start=1 if i == 0 else (length - 100 if i == 59 else 500 + i * 15000),
                         end=length if i == 59 else 600 + i * 15000,
                         strand='+' if i % 2 == 0 else '-') for i in range(60)]
            table, output = scratch / 'genes.tsv', scratch / 'promoters.fa'
            pd.DataFrame(rows).to_csv(table, sep='\t', index=False)
            sys.argv = ['get_promoter_fasta.py', '--dir_genome', str(genomes), '--geneinfo_tsv', str(table),
                        '--outfile', str(output), '--ncpu', '1', '--promoter_bp', '2000']
            with contextlib.redirect_stdout(io.StringIO()):
                started = time.perf_counter()
                module.main()
                elapsed = time.perf_counter() - started
            payload = output.read_bytes()
            counts.update(genes=60, genome_bytes=(genomes / 'Species_a.fa').stat().st_size)
        elif args.case == 'bootstrap':
            import numpy as np
            import subgenome_dominance as module

            values = np.sin(np.arange(10000) * .125) + .1
            blocks = [f'block{i:05d}' for i in range(len(values))]
            rng, null_rng = np.random.default_rng(123), np.random.default_rng(456)
            started = time.perf_counter()
            result = module.inference(values, blocks, 2000, rng, null_rng=null_rng,
                                      permutation_replicates=100, exact_max_states=1)
            elapsed = time.perf_counter() - started
            payload = json.dumps([result, rng.bit_generator.state, null_rng.bit_generator.state], sort_keys=True).encode()
            counts.update(blocks=10000, bootstrap_replicates=2000, permutation_replicates=100)
        elif args.case in ('geometry', 'joint-layout'):
            import pairwise_synteny_layout as module

            if args.case == 'geometry':
                genes = [SimpleNamespace(gene_id=f'g{i:06d}', seqid=f'chr{i % 512:04d}', start=i // 512) for i in range(40000)]
                seqids = [f'chr{i:04d}' for i in range(512)]
                started = time.perf_counter()
                result = module.track_geometry(seqids, genes)
                elapsed = time.perf_counter() - started
                counts.update(genes=40000, scaffolds=512)
            else:
                genes = [[SimpleNamespace(gene_id=f'{side}_g{i:06d}', seqid=f'chr{i % 5}', start=i // 5)
                          for i in range(20000)] for side in ('a', 'b')]
                selected = [[f'chr{i}' for i in range(5)] for _ in range(2)]
                simple = scratch / 'blocks.simple'
                simple.write_text('\n'.join(f'a_g{i:06d}\ta_g{i+50:06d}\tb_g{i:06d}\tb_g{i+50:06d}\t11\t+'
                                            for i in range(0, 20000 - 50, 250)) + '\n')
                started = time.perf_counter()
                result = module.order_both_by_ribbon_length(selected, genes, simple)
                elapsed = time.perf_counter() - started
                counts.update(genes_per_track=20000, chromosomes=5, track_solves=result[1]['track_solves'])
            payload = json.dumps(result, sort_keys=True).encode()
        elif args.case == 'traits':
            import inspect

            import get_trait_matrix as module

            source = args.fixture / 'traits.tsv'
            queries = [[f'g{(family * 191 + j * 997) % 40000:05d}' for j in range(20)] for family in range(20)]
            kwargs = {'cache_dir': scratch / 'cache'} if 'cache_dir' in inspect.signature(module.process_trait_file).parameters else {}
            started = time.perf_counter()
            results = []
            for i, ids in enumerate(queries):
                results.append(module.process_trait_file(source, set(ids), {g: 'Species_a_' + g for g in ids}, **kwargs))
                if i == 0:
                    counts['first_family_seconds'] = time.perf_counter() - started
            elapsed = time.perf_counter() - started
            payload = b''.join(frame_payload(frame) for frame in results)
            counts.update(genes=40000, traits=80, families=20, selected_per_family=20, source_bytes=source.stat().st_size)
        elif args.case == 'cds':
            import inspect

            import gff2genestat as module
            import pandas as pd

            source = scratch / 'Species_a.traits.tsv'
            sequence = 'ATGAAATAA'
            frame = pd.DataFrame({'gene_id': [f'Species_a_g{i}' for i in range(40000)],
                'cds_sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest(),
                'structure_status': 'sequence_verified', 'cds_first_phase': 0, 'feature_size': 9,
                'phase_status': 'valid', 'feature_blocks': '1-9'})
            frame.to_csv(source, sep='\t', index=False)
            (scratch / 'Species_a.resolution.meta.json').write_text(json.dumps({
                'traits_sha256': hashlib.sha256(source.read_bytes()).hexdigest(), 'contract': {'sources': {}}}))
            kwargs = {'cache_dir': scratch / 'cache'} if 'cache_dir' in inspect.signature(module.apply_cds_resolution).parameters else {}
            results = []
            inputs = [(frame.iloc[i*20:(i+1)*20].copy(), [(g, '', sequence) for g in frame.gene_id[i*20:(i+1)*20]]) for i in range(200)]
            started = time.perf_counter()
            for i, (selected, records) in enumerate(inputs):
                module.apply_cds_resolution(selected, records, scratch, **kwargs)
                results.append(selected)
                if i == 0:
                    counts['first_family_seconds'] = time.perf_counter() - started
            elapsed = time.perf_counter() - started
            payload = b''.join(frame_payload(frame) for frame in results)
            counts.update(genes=40000, families=200, selected_per_family=20, external_source_hashes=0)
        elif args.case == 'hgt':
            import focus_hgt_context as module

            rng = random.Random(31)
            rows = [dict(gene_id=f'g{i}', chromosome='s1', start=str(1+i*100),
                         end=str(100+i*100+rng.randrange(1000)), feature_type='gene') for i in range(12000)]
            for i in range(0, len(rows), 500):
                rows[i]['end'] = '2000000'
            header = list(rows[0])
            (scratch / 'Species_a.gff_info.tsv').write_text('\t'.join(header)+'\n'+''.join('\t'.join(r[k] for k in header)+'\n' for r in rows))
            coordinates = module.GenomeCoordinates(scratch)
            links = [dict(gene_species='Species_a', gene_id=f'g{i}') for i in range(50, 12000, 23)]
            started = time.perf_counter()
            results = []
            for link in links:
                neighborhood = coordinates.neighborhood(link)
                models = coordinates.scaffold_models(link, neighborhood[0])
                results.append((neighborhood, [r['gene_id'] for r in models]))
            coordinates.verify()
            elapsed = time.perf_counter() - started
            payload = json.dumps(results, sort_keys=True).encode()
            counts.update(genes=12000, events=len(links))
        elif args.case == 'dated':
            # Measure only the real tip-positioning block, excluding rendering.
            source = (args.support / 'dated_tree_presentation.py').read_text()
            end = source.index('\n    for node in tree.traverse', source.index('    ys = {leaf:'))
            begin = source.rfind('\n    ranks = ', 0, end)
            if begin < 0:
                begin = source.index('    ys = {leaf:')
            scope = {}
            exec('def position(leaves, order, count):\n' + source[begin:end].lstrip('\n') + '\n    return ys', scope)
            leaves = [type('Leaf', (), {'name': f't{i}'})() for i in range(12000)]
            order = [leaf.name for leaf in leaves][::-1]
            started = time.perf_counter()
            result = scope['position'](leaves, order, len(leaves))
            elapsed = time.perf_counter() - started
            payload = json.dumps(list(result.values())).encode()
            counts.update(tips=len(leaves))
        elif args.case == 'kmer':
            import struct

            executable = scratch / 'kernel'
            subprocess.run(['g++', '-std=c++17', '-O3', '-pthread', str(args.support / 'busco_kmer_distance.cpp'), '-o', str(executable)], check=True)
            sketch = scratch / 'sketch'
            sketch.write_bytes(b'GGKMER01' + struct.pack('<QQQQQ', 5, 16, 2, 0, 0))
            listing, output = scratch / 'list', scratch / 'distances'
            listing.write_text((str(sketch)+'\n')*1600)
            started = time.perf_counter()
            measurement = subprocess.run([sys.executable, '-c',
                'import resource,subprocess,sys;subprocess.run(sys.argv[1:],check=True);print(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)',
                str(executable), 'compare', str(listing), str(output), '2'], check=True, capture_output=True, text=True)
            elapsed = time.perf_counter() - started
            counts['kernel_peak_rss_kib_linux'] = int(measurement.stdout)
            payload = output.read_bytes()
            counts.update(species=1600, markers=2, note='Empty sketches isolate pair storage/output; compilation excluded; kernel RSS measured in separate parent.')
        elif args.case == 'rooting':
            import re

            import species_tree_guided_gene_tree_rooting as module

            source, candidates = scratch / 'in.nwk', scratch / 'candidates'
            source.write_text('((A:7.8,B:0.949):0.244,(C:7.8,D:1.439):2.87);')
            candidates.write_bytes(source.read_bytes())
            options = SimpleNamespace(in_tree=source, candidate_trees=candidates, out_tree=scratch/'selected',
                comparison_table=scratch/'tsv', comparison_plot=scratch/'plot.pdf', species_parser='taxonomic')
            with contextlib.redirect_stdout(io.StringIO()):
                started = time.perf_counter()
                module.select_root(options)
                elapsed = time.perf_counter() - started
            payload = options.out_tree.read_bytes() + options.comparison_table.read_bytes() + re.sub(
                rb'/(CreationDate|ModDate) \([^)]*\)', b'', options.comparison_plot.read_bytes())
            counts.update(tips=4)
        else:
            command = ['Rscript', str(Path(__file__).with_suffix('.R')), str(args.support), args.case, str(args.fixture), str(scratch)]
            completed = subprocess.run(command, capture_output=True, text=True, check=True)
            result = json.loads(completed.stdout.splitlines()[-1])
            return dict(result, peak_rss_kib_linux=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
        return dict(seconds=elapsed, peak_rss_kib_linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                    child_peak_rss_kib_linux=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                    output_sha256=hashlib.sha256(payload).hexdigest(), counts=counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-support', type=Path)
    parser.add_argument('--support', type=Path, default=ROOT / 'workflow/support', help=argparse.SUPPRESS)
    parser.add_argument('--case', choices=CASES, help=argparse.SUPPRESS)
    parser.add_argument('--fixture', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--prepare', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--cases', nargs='+', choices=CASES, default=list(CASES))
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if sys.platform != "linux":
        parser.error("Run in the Linux GeneGalleon runtime; RSS is reported in Linux KiB")
    if args.prepare:
        prepare(args.fixture)
        return
    if args.case:
        # Child tools may log to stdout; reserve the final line for the record.
        print(json.dumps(worker(args)))
        return
    if not args.baseline_support or args.repeats < 1:
        parser.error('--baseline-support and positive --repeats are required')
    sources = ('get_promoter_fasta.py', 'subgenome_dominance.py', 'pairwise_synteny_layout.py',
               'get_trait_matrix.py', 'gff2genestat.py', 'focus_hgt_context.py',
               'dated_tree_presentation.py', 'busco_kmer_distance.cpp', 'species_tree_guided_gene_tree_rooting.py',
               'annotation_summary.r', 'trait_table_index.py')
    fingerprints = {label: {name: hashlib.sha256((support / name).read_bytes()).hexdigest()
                            for name in sources if (support / name).is_file()}
                    for label, support in [('baseline', args.baseline_support), ('current', args.support)]}
    result = dict(environment=dict(platform=platform.platform(), python=sys.version,
        packages={name: importlib.metadata.version(name) for name in ('numpy', 'pandas', 'nwkit', 'biopython')},
        tools={name: subprocess.check_output(command, text=True).splitlines()[0] for name, command in (
            ('seqkit', ['seqkit', 'version']), ('R', ['Rscript', '--version']), ('compiler', ['g++', '--version']))}),
        source_sha256=fingerprints, workloads={})
    with tempfile.TemporaryDirectory(prefix='gg-efficiency-fixtures-') as temporary:
        fixture = Path(temporary)
        subprocess.run([sys.executable, __file__, '--prepare', '--fixture', str(fixture), '--output', str(args.output)], check=True)
        for case in args.cases:
            samples = {label: [] for label in ('baseline', 'current')}
            hashes = set()
            for trial in range(args.repeats + 1):
                for label in (samples if trial % 2 == 0 else reversed(samples)):
                    support = args.baseline_support if label == 'baseline' else args.support
                    completed = subprocess.run([sys.executable, __file__, '--case', case, '--support', str(support),
                        '--fixture', str(fixture), '--output', str(args.output)], capture_output=True, text=True)
                    if completed.returncode:
                        raise RuntimeError(completed.stdout + completed.stderr)
                    row = json.loads(completed.stdout.splitlines()[-1])
                    hashes.add(row['output_sha256'])
                    if trial:
                        samples[label].append(row)
            if len(hashes) != 1:
                raise ValueError('Output mismatch: ' + case)
            result['workloads'][case] = {label: dict(trials=rows,
                median_seconds=statistics.median(row['seconds'] for row in rows),
                median_peak_rss_kib_linux=statistics.median(row['peak_rss_kib_linux'] for row in rows)) for label, rows in samples.items()}
            args.output.write_text(json.dumps(result, indent=2) + '\n')
            print(case, {label: result['workloads'][case][label]['median_seconds'] for label in samples}, flush=True)


if __name__ == '__main__':
    main()
