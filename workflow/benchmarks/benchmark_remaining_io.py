#!/usr/bin/env python3
"""Compare verify, DB, store fingerprint and PDF on deterministic private fixtures.

Use --support-root for a saved baseline. Each sample runs in a fresh process;
one warmup is discarded. Timings include the operation, not fixture creation.
"""
import argparse
import contextlib
import hashlib
import io
import json
import os
import platform
import re
import resource
import sqlite3
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def hash_bytes(data):
    return hashlib.sha256(data).hexdigest()


def db_fingerprint(path):
    digest = hashlib.sha256()
    with sqlite3.connect(path) as conn:
        for name, sql in conn.execute("SELECT name,sql FROM sqlite_master WHERE type='table' ORDER BY name"):
            digest.update(sql.encode())
            cols = [row[1] for row in conn.execute(f'PRAGMA table_info("{name}")')]
            order = ','.join('"' + col + '"' for col in cols)
            for row in conn.execute(f'SELECT * FROM "{name}" ORDER BY {order}'):
                digest.update(json.dumps(row, separators=(',', ':')).encode() + b'\n')
    return digest.hexdigest()


def worker(args):
    sys.path.insert(0, str(args.support_root))
    import artifact_provenance as p
    import workflow_api as api
    from gene_family_output_store import read_only_observation
    from workflow_observation import observe_files
    root = args.worker / 'output/orthogroup'
    root.mkdir(parents=True)
    if args.case in ('verify', 'verify-batch'):
        source = args.worker / 'input/shared.tsv'
        source.parent.mkdir()
        source.write_bytes(b'x' * (64 * 1024 * 1024))
        attempt = args.worker / ('a' * 32)
        attempt.mkdir()
        started = time.time_ns()
        for step in ('iqtree_anc', 'csubst', 'csubst_scan', 'summary_statistics', 'tree_plot'):
            output = root / 'stat_branch/OG0000001_stat.branch.tsv'
            output.parent.mkdir(exist_ok=True)
            output.write_text('branch_id\tnode_name\tnum_sp\tso_event\n0\tg1\t1\tL\n')
            manifest = root / ('artifact_provenance/OG0000001.' + step + '.json')
            argv = ['--workspace-root', str(args.worker), '--logical-root', str(root),
                    '--manifest', str(manifest), '--family-id', 'OG0000001', '--step', step,
                    '--input', 'shared=' + str(source), '--output', 'table=' + str(output)]
            os.environ['GG_OBSERVATION_ATTEMPT_DIR'] = str(attempt)
            with contextlib.redirect_stdout(io.StringIO()):
                assert p.dispatch(['record', *argv]) == 0
        (attempt / 'run.json').write_text(json.dumps({'schema': 'genegalleon-observation-v1',
            'attempt_id': attempt.name, 'workflow': 'gg_gene_evolution', 'execution_state': 'exited',
            'execution_accepted': True, 'exit_code': 0, 'accepted_exit_codes': [0],
            'started_at_ns': started, 'finished_at_ns': time.time_ns()}))
        ns = argparse.Namespace(root=root, workspace_root=args.worker, family_id='OG0000001',
            require_step=['iqtree_anc', 'csubst', 'csubst_scan', 'summary_statistics', 'tree_plot'],
            manifest=[], profile=None, attempt=attempt, recorded_workspace_root=None, include_queue=False)
        requests = [ns]
        if args.case == 'verify-batch':
            for index in range(2, 17):
                family = f'OG{index:07d}'
                other = argparse.Namespace(**{**vars(ns), 'family_id': family, 'attempt': None})
                for step in ns.require_step:
                    payload = json.loads((root / f'artifact_provenance/OG0000001.{step}.json').read_text())
                    payload['family_id'] = family
                    (root / f'artifact_provenance/{family}.{step}.json').write_text(json.dumps(payload))
                requests.append(other)
        begin = time.perf_counter()
        with read_only_observation():
            if args.case == 'verify-batch' and hasattr(api, 'verify_many'):
                results = api.verify_many(requests)
            else:
                results = [api.verify(item) for item in requests]
        result = results[0]
        assert all(item['completion_state'] == 'verified_declared_steps' for item in results)
        seconds = time.perf_counter() - begin
        assert result['completion_state'] == 'verified_declared_steps', result
        output = hash_bytes(json.dumps([(r['step'], r['state'], r['attempt_bound'],
            p.contract_comparison_payload(json.loads((root / r['manifest']).read_text())))
            for r in result['contracts']], sort_keys=True).encode())
    elif args.case in ('database-raw', 'database-zip', 'store-digest'):
        for family in range(256):
            name = f'OG{family:07d}'
            for subdir, suffix, header, rows in (
                ('stat_tree', '_stat.tree.tsv', 'num_branch\tnum_spe\tnum_dup\tnum_sp', ['128\t64\t0\t63']),
                ('stat_branch', '_stat.branch.tsv', 'branch_id\tnode_name\tnum_sp\tso_event\tmetric',
                 [f'{i}\tn{i}\t2\tS\t{family+i}' for i in range(128)])):
                path = root / subdir / (name + suffix)
                path.parent.mkdir(exist_ok=True)
                path.write_text(header + '\n' + '\n'.join(rows) + '\n')
        if args.case == 'database-zip':
            catalog = args.worker / 'families.txt'
            catalog.write_text(''.join(f'OG{i:07d}\n' for i in range(256)))
            with contextlib.redirect_stdout(io.StringIO()):
                subprocess.run([sys.executable, str(args.support_root / 'gene_family_output_store.py'),
                    'convert-storage', '--root', str(root), '--mode', 'orthogroup', '--to', 'zip',
                    '--family-id-file', str(catalog), '--compression', 'store', '--progress-interval', '0'],
                    check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if args.case == 'store-digest':
            for i in range(64):
                path = root / 'alignment' / f'OG{i:07d}.fa'
                path.parent.mkdir(exist_ok=True)
                path.write_bytes(b'>gene\n' + b'ACGT' * (256 * 1024))
            p.configure_digest_cache(None)
            begin = time.perf_counter()
            with read_only_observation(), observe_files():
                result = p.gene_family_store_digest(root)
            seconds = time.perf_counter() - begin
            output = hash_bytes(json.dumps(result).encode())
        else:
            import generate_orthogroup_database as db
            database = args.worker / 'result.db'
            sys.argv = [str(args.support_root / 'generate_orthogroup_database.py'), '--dbpath', str(database),
                '--dir_gene_family', str(root), '--dir_stat_tree', str(root / 'stat_tree'),
                '--dir_stat_branch', str(root / 'stat_branch'), '--ncpu', '2', '--row_threshold', '4096']
            begin = time.perf_counter()
            db.main()
            seconds = time.perf_counter() - begin
            output = db_fingerprint(database)
    else:
        stat = args.worker / 'stat.tsv'
        stat.write_text('branch_id\tparent\tsister\tchild1\tchild2\tnode_name\tbl_rooted\tso_event\tso_event_parent\n'
            '4\t-999\t-999\t2\t3\troot\t0\tS\tS\n2\t4\t3\t0\t1\tn4\t1\tS\tS\n'
            '0\t2\t1\t\t\tg1\t1\tL\tS\n1\t2\t0\t\t\tg2\t1\tL\tS\n3\t4\t2\t\t\tg3\t2\tL\tS\n')
        os.chdir(args.worker)
        panel_args = []
        if args.case == 'pdf-panels':
            domain = args.worker / 'domain.tsv'
            domain.write_text('qacc\tsacc\tstitle\tqlen\tslen\tqstart\tqend\n'
                              'g1\tPF0001\tPF0001,Domain_A\t30\t300\t1\t10\n'
                              'g2\tPF0002\tPF0002,Domain_B\t30\t300\t11\t20\n')
            fasta = args.worker / 'alignment.fa'
            fasta.write_text(''.join(f'>g{i}\n' + 'ACGT' * 21 + 'ACGTAA\n' for i in range(1, 4)))
            panel_args = ['--panel2=domain,' + str(domain),
                          '--panel3=alignment,' + str(fasta) + ',' + str(fasta)]
        begin = time.perf_counter()
        script = args.support_root / 'stat_branch2tree_plot.r'
        # Baseline core uses a separate R dependency probe; new core delegates
        # that same check to the rendering process.
        integrated = 'GG_TREE_PLOT_CHECK_GGIMAGE' in script.read_text()
        if not integrated:
            subprocess.run(['Rscript', '-e', "if (!requireNamespace('ggimage', quietly=TRUE)) quit(status=1)"],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        command = ['Rscript', str(script), '--stat_branch=' + str(stat),
            '--max_delta_intron_present=-0.5', '--panel_widths_mm=tree:60', '--panel1=tree,bl_rooted,no,no,L',
            '--show_branch_id=no', '--event_method=species_overlap', '--species_color_table=PLACEHOLDER',
            '--pie_chart_value_transformation=identity', '--long_branch_display=no', *panel_args]
        outputs = []
        if args.case == 'pdf-batch' and (args.support_root / 'tree_plot_batch.r').exists():
            jobs = [{'id':str(i), 'cwd':str(args.worker), 'args':command[2:],
                     'output':str(args.worker / f'{i}.pdf')} for i in range(8)]
            plan = args.worker / 'plot-plan.json'
            plan.write_text(json.dumps(jobs))
            subprocess.run(['Rscript', str(args.support_root / 'tree_plot_batch.r'), str(plan)], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            outputs = [(args.worker / f'{i}.pdf').read_bytes() for i in range(8)]
        else:
            for _ in range(8 if args.case == 'pdf-batch' else 1):
                subprocess.run(command, check=True,
                    env={**os.environ, 'GG_TREE_PLOT_CHECK_GGIMAGE': '1'}, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                outputs.append(Path('stat_branch2tree_plot.pdf').read_bytes())
        seconds = time.perf_counter() - begin
        # R PDFs differ only in creation/modification dates. Compare all other
        # bytes, including compressed drawing commands, fonts and page sizes.
        normalized = [hash_bytes(re.sub(rb'/(CreationDate|ModDate) \([^)]*\)', b'', raw)) for raw in outputs]
        assert len(set(normalized)) == 1
        output = normalized[0]
    (args.worker / 'measurement.json').write_text(json.dumps({'seconds': seconds, 'output_sha256': output,
        'peak_rss_kib_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'peak_child_rss_kib_linux': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--support-root', type=Path, default=Path(__file__).resolve().parents[1] / 'support')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--case', choices=['verify', 'verify-batch', 'database-raw', 'database-zip', 'store-digest', 'pdf', 'pdf-panels', 'pdf-batch'])
    parser.add_argument('--worker', type=Path)
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    args.support_root = args.support_root.resolve()
    if args.worker:
        return worker(args)
    if not args.output or args.repeats < 2:
        parser.error('--output and at least two measured repeats are required')
    results = {'support_root': str(args.support_root), 'python': sys.version,
               'platform': platform.platform(), 'machine': platform.machine(), 'cases': {}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for case in ([args.case] if args.case else ['verify', 'verify-batch', 'database-raw', 'database-zip', 'store-digest', 'pdf', 'pdf-panels', 'pdf-batch']):
        samples = []
        for i in range(args.repeats + 1):
            with tempfile.TemporaryDirectory(prefix='gg-remaining-benchmark-') as tmp:
                with args.output.with_suffix('.log').open('a') as log:
                    subprocess.run([sys.executable, str(Path(__file__).resolve()), '--worker', tmp,
                        '--case', case, '--support-root', str(args.support_root)], check=True,
                        stdout=log, stderr=subprocess.STDOUT)
                sample = json.loads((Path(tmp) / 'measurement.json').read_text())
            if i:
                samples.append(sample)
        assert len({sample['output_sha256'] for sample in samples}) == 1
        results['cases'][case] = {'samples': samples,
            'median_seconds': statistics.median(s['seconds'] for s in samples)}
        args.output.write_text(json.dumps(results, indent=2) + '\n')
        print(case, results['cases'][case]['median_seconds'], flush=True)


if __name__ == '__main__':
    main()
