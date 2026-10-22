#!/usr/bin/env python3
"""Prepare an isolated upstream OrionGeno installation and run native GPU GFF3 prediction."""
import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import shutil
import subprocess
import tempfile
from pathlib import Path

from annevo_model_selection import protect_outputs, verified_input_sha
from annevo_runtime import gff_counts, publish_pair, selection_record
from annevo_runtime_cache import RuntimeVerification, inspect_tree, private_root
from annevo_runtime_download import atomic_json, canonical_bytes, file_sha
from format_species_writers import input_compression
from genome_sequence_identity import genome_sequence_digest
from oriongeno_model_profiles import PROFILES
from oriongeno_outputs import (
    AUDIT_SUFFIXES,
    UTR_POLICY,
    artifact_record,
    companion_paths,
    compress_audit_jsonl,
    validate_repeat_gff,
    write_explicit_utrs,
)
from shared_namespace_lock import namespace_lock
from transcript_utrs import load_models

REPOSITORY = 'https://github.com/BGIResearch/OrionGeno.git'


def environment(root, threads=1):
    env = dict(os.environ)
    for key in ('PYTHONPATH', 'PYTHONHOME', 'PYTHONPYCACHEPREFIX'):
        env.pop(key, None)
    for key, name in [('HF_HOME', 'huggingface'), ('XDG_CACHE_HOME', 'xdg'),
                      ('TORCH_EXTENSIONS_DIR', 'torch_extensions'), ('TMPDIR', 'tmp'),
                      ('NUMBA_CACHE_DIR', 'numba'), ('TORCH_HOME', 'torch'),
                      ('MAMBA_ROOT_PREFIX', 'mamba'), ('CONDA_PKGS_DIRS', 'conda_pkgs'),
                      ('PIP_CACHE_DIR', 'pip'), ('TRITON_CACHE_DIR', 'triton'),
                      ('HF_HUB_CACHE', 'huggingface/hub'), ('HF_XET_CACHE', 'huggingface/xet')]:
        directory = Path(root) / name
        directory.mkdir(parents=True, exist_ok=True)
        env[key] = str(directory)
    env.update(PYTHONDONTWRITEBYTECODE='1', PYTHONNOUSERSITE='1',
               OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS=str(threads),
               NUMBA_NUM_THREADS=str(threads))
    # GG's production image carries the Conda C toolchain under its target
    # name. Triton otherwise searches only for cc/gcc/clang at first inference.
    if not env.get('CC'):
        compiler = shutil.which('cc') or shutil.which('gcc') or shutil.which('x86_64-conda-linux-gnu-cc')
        if compiler:
            env['CC'] = compiler
    return env


def execute(argv, env, **kwargs):
    return subprocess.run(list(map(str, argv)), env=env, check=True, text=True, **kwargs)


def prepare(args):
    if args.model_profile not in PROFILES:
        raise ValueError('Unsupported OrionGeno model profile')
    root = private_root(Path(args.cache_root))
    receipt = Path(args.receipt)
    with namespace_lock(root / '.prepare.lock', exclusive=True, timeout=86400):
        if receipt.exists():
            value = read_selected(receipt)
            if value['model_profile'] != args.model_profile:
                raise ValueError('OrionGeno model profile changed; use a new runtime receipt')
            if Path(value['cache_root']) != root or (args.source_commit and args.source_commit != value['source_commit']):
                raise ValueError('OrionGeno runtime override differs from its receipt')
            if getattr(args, 'source_ref', '') and args.source_ref != value.get('source_ref'):
                raise ValueError('OrionGeno source branch differs from its receipt')
            for option, key in [('python', 'python'), ('checkpoint', 'model_path')]:
                supplied = getattr(args, option, '')
                if supplied:
                    if option == 'python':
                        current = Path(os.path.abspath(Path(supplied).expanduser()))
                        recorded = Path(value[key])
                        # Keep aliases in one environment (including old
                        # receipts), but distinguish venvs sharing one binary.
                        changed = ((current.parent.resolve(), current.resolve())
                                   != (recorded.parent.resolve(), recorded.resolve()))
                    else:
                        changed = Path(supplied).resolve() != Path(value[key]).resolve()
                    if changed:
                        raise ValueError('OrionGeno runtime override differs from its receipt')
            return value
        env = environment(root)
        # A venv interpreter commonly links to the base Python. Dereferencing
        # that final link discards pyvenv.cfg and imports from the wrong prefix.
        python = Path(os.path.abspath(Path(args.python).expanduser())) if args.python else root / 'environment/bin/python'
        source_commit, source_ref, source = '', getattr(args, 'source_ref', ''), None
        if not args.python:
            if args.no_download:
                raise ValueError('OrionGeno runtime cache is cold and downloads are disabled')
            if platform.system() != 'Linux' or platform.machine() != 'x86_64':
                raise ValueError('The upstream OrionGeno native wheels require Linux x86_64')
            source = root / 'source'
            source_ref = getattr(args, 'source_ref', '')
            if not source_ref:
                refs = Path(__file__).resolve().parents[2] / 'container/source_branches.env'
                if not refs.is_file():
                    raise ValueError('OrionGeno requires --source-ref outside a repository checkout')
                source_ref = dict(line.split('=', 1) for line in refs.read_text().splitlines() if line.startswith('GG_SOURCE_'))['GG_SOURCE_ORIONGENO_REPO_REF']
            if not source.exists():
                execute(['git', 'clone', '--depth', '1', '--branch', source_ref, REPOSITORY, source], env)
            source_commit = execute(['git', '-C', source, 'rev-parse', 'HEAD'], env, capture_output=True).stdout.strip()
            if (root/'installation.ready').exists() and (root/'installation.ready').read_text().strip() != source_commit:
                raise ValueError('OrionGeno source changed after environment installation')
            manager = shutil.which('conda') or shutil.which('micromamba')
            if not manager:
                raise ValueError('OrionGeno bootstrap requires conda, or an explicit --python installation')
            if not (root/'installation.ready').exists():
                if not python.exists():
                    execute([manager, 'env', 'create', '-y', '-p', python.parent.parent, '-f', source / 'environment.min.yml'], env)
                execute([python, '-B', '-m', 'pip', 'install', '--no-deps', '-r', source / 'requirements.native-cu126.txt'], env)
                execute([python, '-B', '-m', 'pip', 'install', '--no-deps', source], env)
                execute([python, '-B', '-s', '-c', 'import torch,oriongeno,mamba_ssm,causal_conv1d,numba'], env)
                (root/'installation.ready').write_text(source_commit+'\n')
        if not python.is_file():
            raise ValueError('OrionGeno Python is missing')
        info = json.loads(execute([python, '-B', '-s', '-c',
            "import json,pathlib,oriongeno;print(json.dumps({'package':str(pathlib.Path(oriongeno.__file__).parent)}))"],
            env, capture_output=True).stdout.strip())
        package = Path(info['package'])
        if not source_commit:
            # An external installation must identify its owning source snapshot.
            source_commit = args.source_commit
            if not source_commit or len(source_commit) != 40 or any(c not in '0123456789abcdef' for c in source_commit):
                raise ValueError('An external OrionGeno installation requires --source-commit provenance')
        checkpoint = Path(args.checkpoint).resolve() if args.checkpoint else root / 'models' / args.model_profile
        if not args.checkpoint:
            script = """from huggingface_hub import snapshot_download
import sys
snapshot_download('BGI-Research/OrionGeno', allow_patterns=['checkpoints/oriongeno_'+sys.argv[1]+'/*'], local_dir=sys.argv[2], local_files_only=sys.argv[3]=='1')
"""
            download_root = root / 'weights'
            execute([python, '-B', '-s', '-c', script, args.model_profile, download_root, '1' if args.no_download else '0'], env)
            checkpoint = download_root / ('checkpoints/oriongeno_'+args.model_profile)
        if not (checkpoint / 'config.json').is_file() or not (
                any(checkpoint.glob('*.safetensors')) or any(checkpoint.glob('pytorch_model*.bin'))):
            raise ValueError('OrionGeno checkpoint lacks config.json or model weights')
        trees = {str(p): inspect_tree(p) for p in (package, checkpoint)}
        model_digest = hashlib.sha256(canonical_bytes({k: v['sha256'] for k, v in trees[str(checkpoint)].items() if v['type'] == 'file'})).hexdigest()
        packages = execute([python, '-B', '-m', 'pip', 'freeze', '--all'], env, capture_output=True).stdout
        value = dict(schema_version=1, kind='genegalleon-oriongeno-installed-runtime',
                     python=str(python), python_sha256=file_sha(python), model_path=str(checkpoint),
                     model_profile=args.model_profile, source_commit=source_commit, source_ref=source_ref,
                     model_sha256=model_digest, trees=trees, packages=packages, cache_root=str(root))
        value['runtime_digest'] = hashlib.sha256(canonical_bytes(value)).hexdigest()
        value['cache_key'] = value['runtime_digest']
        atomic_json(receipt, value)
        return value


def read_selected(path, *, verification=None):
    verification = verification if verification is not None else RuntimeVerification()
    with verification.guard():
        verification.file_sha(path)
        value = json.loads(Path(path).read_text())
        if value.get('kind') != 'genegalleon-oriongeno-installed-runtime' or value.get('schema_version') != 1:
            raise ValueError('Invalid OrionGeno installed runtime receipt')
        check = {k: v for k, v in value.items() if k not in {'runtime_digest', 'cache_key'}}
        if (hashlib.sha256(canonical_bytes(check)).hexdigest() != value['runtime_digest']
                or value.get('cache_key') != value['runtime_digest']
                or verification.file_sha(value['python']) != value['python_sha256']):
            raise ValueError('OrionGeno runtime receipt or Python changed')
        for root, tree in value['trees'].items():
            verification.verify_tree(Path(root), tree)
        packages = execute([value['python'], '-B', '-m', 'pip', 'freeze', '--all'], environment(value['cache_root']), capture_output=True).stdout
        if packages != value['packages']:
            raise ValueError('OrionGeno installed dependencies changed')
        verification.check()
        return value


def run(args):
    if args.threads < 1 or args.batch_size < 1:
        raise ValueError('OrionGeno threads and batch size must be positive')
    output_repeat = getattr(args, 'output_repeat', 1)
    output_utrs = getattr(args, 'output_utrs', 1)
    if output_repeat not in (0, 1) or output_utrs not in (0, 1):
        raise ValueError('OrionGeno output flags must be 0 or 1')
    verification = RuntimeVerification()
    runtime = read_selected(args.receipt, verification=verification)
    genome, output, receipt = map(Path, (args.genome, args.output, args.output_receipt))
    if not genome.is_file() or not genome.stat().st_size:
        raise ValueError('OrionGeno genome is missing or empty')
    selection_options = dict(engine='oriongeno', model_profile=runtime['model_profile'], verification=verification)
    selection, selection_sha = selection_record(args, genome, **selection_options)
    protected = {genome.resolve(), Path(args.receipt).resolve(), Path(runtime['python']).resolve()}
    if selection is not None:
        protected.add(Path(args.selection_receipt).resolve())
        protected.update(Path(evidence['path']).resolve() for evidence in selection.get('evidence', []))
        if selection.get('taxonomy'):
            protected.add(Path(selection['taxonomy']['path']).resolve())
    if output.resolve() == receipt.resolve() or any(p.resolve() in protected or any(p.resolve().is_relative_to(Path(t).resolve()) for t in runtime['trees']) for p in (output, receipt)):
        raise ValueError('OrionGeno outputs overwrite inputs or runtime')
    env = environment(runtime['cache_root'], args.threads)
    logs = output.with_name(output.name+'.logs')
    companions = companion_paths(output)
    if not output_repeat:
        companions.pop('repeat_gff')
    publications = (output, receipt, *companions.values())
    for index, destination in enumerate(publications):
        protect_outputs([destination], (*protected, *publications[:index]))
    log_destinations = (logs, logs/'native.stdout.txt', logs/'native.stderr.txt')
    protect_outputs(log_destinations, (*protected, *publications))
    destinations = (*publications, *log_destinations)
    if any(p.resolve().is_relative_to(Path(t).resolve()) for p in destinations for t in runtime['trees']):
        raise ValueError('OrionGeno outputs overwrite the installed runtime')
    installed_files = {tuple(member['signature'][:2]) for tree in runtime['trees'].values()
                       for member in tree.values() if member['type'] == 'file'}
    for path in destinations:
        if path.exists():
            info = path.stat()
            if (info.st_dev, info.st_ino) in installed_files:
                raise ValueError('OrionGeno outputs overwrite an installed runtime file through a hardlink')
    gpu_code = """import json,torch
import mamba_ssm,causal_conv1d
if not torch.cuda.is_available(): raise SystemExit('OrionGeno requires an allocated CUDA GPU')
x=torch.ones((32,32),device='cuda'); y=x@x; torch.cuda.synchronize()
assert y.is_cuda and y[0,0].item()==32
print(json.dumps({'gpu_verified':True,'torch_version':torch.__version__,'device_name':torch.cuda.get_device_name()}))
"""
    gpu = json.loads(execute([runtime['python'], '-B', '-s', '-c', gpu_code], env, capture_output=True).stdout.strip().splitlines()[-1])
    before = verified_input_sha(genome, verification)
    contig_lengths = {}
    sequence_identity = genome_sequence_digest(genome, contig_lengths=contig_lengths)
    if verified_input_sha(genome, verification) != before:
        raise ValueError('OrionGeno genome changed during sequence identity reading')
    output.parent.mkdir(parents=True, exist_ok=True)
    receipt.parent.mkdir(parents=True, exist_ok=True)
    with namespace_lock(output.with_suffix(output.suffix+'.lock'), exclusive=True, timeout=86400):
        with tempfile.TemporaryDirectory(prefix='.oriongeno-', dir=output.parent) as temporary:
            predicted = Path(temporary) / 'prediction.gff3'
            staged = Path(temporary) / 'annotated.gff3'
            native_genome = genome.resolve()
            if input_compression(genome) == 'gzip':
                native_genome = Path(temporary)/'genome.fa'
                with gzip.open(genome, 'rb') as source, native_genome.open('wb') as target:
                    shutil.copyfileobj(source, target)
            command = [runtime['python'], '-B', '-s', '-m', 'oriongeno', '--genome', str(native_genome),
                       '--checkpoint', runtime['model_path'], '--output', str(predicted),
                       '--batch-size', str(args.batch_size), '--output-gene', 'True', '--output-repeat', str(bool(output_repeat)),
                       '--species-name', '']
            logs.mkdir(exist_ok=True)
            with (logs/'native.stdout.txt').open('w') as stdout, (logs/'native.stderr.txt').open('w') as stderr:
                execute(command, env, stdout=stdout, stderr=stderr, cwd=temporary)
            gff_counts(predicted)
            utr_export = write_explicit_utrs(predicted, staged, lengths=contig_lengths) if output_utrs else None
            if not output_utrs:
                load_models(predicted)
                shutil.copyfile(predicted, staged)
            counts = gff_counts(staged)
            companion_sources = {'native_gff': predicted}
            if output_repeat:
                repeat = Path(temporary) / 'prediction.repeat.gff3'
                repeat_counts = validate_repeat_gff(repeat, contig_lengths)
                companion_sources['repeat_gff'] = repeat
            audit_metadata = {}
            for name, suffix in AUDIT_SUFFIXES.items():
                compressed = Path(temporary) / companions[name].name
                audit_metadata[name] = compress_audit_jsonl(Path(str(predicted) + suffix), compressed)
                companion_sources[name] = compressed
            if verified_input_sha(genome, verification) != before:
                raise ValueError('OrionGeno genome changed during prediction')
            if selection_sha is not None:
                _, current_selection_sha = selection_record(args, genome, **selection_options)
                if current_selection_sha != selection_sha:
                    raise ValueError('OrionGeno model selection changed during prediction')
            result = {k: runtime[k] for k in ('runtime_digest', 'cache_key', 'source_commit', 'model_profile', 'model_sha256')}
            result.update(schema_version=1, kind='genegalleon-oriongeno-native-prediction', native_returncode=0,
                          genome_sha256=before, genome_sequence_sha256=sequence_identity,
                          GFF_sha256=file_sha(staged), GFF_bytes=staged.stat().st_size, gpu_probe=gpu,
                          inference_settings={'lineage': runtime['model_profile'], 'batch_size': args.batch_size, 'threads': args.threads,
                                              'output_repeat': bool(output_repeat), 'output_utrs': bool(output_utrs),
                                              'audit_log_compression': 'gzip'},
                          argv=command, feature_counts=counts)
            result['annotation_outputs'] = dict(schema_version=2, output_repeat=bool(output_repeat), output_utrs=bool(output_utrs),
                audit_log_compression='gzip',
                utr_policy=UTR_POLICY if output_utrs else 'native', utr_export=utr_export,
                files={key: artifact_record(path, companions[key]) | audit_metadata.get(key, {})
                       for key, path in companion_sources.items()})
            if output_repeat:
                result['annotation_outputs']['repeat_feature_counts'] = repeat_counts
            if selection is not None:
                result.update(model_selection=selection, model_selection_sha256=selection_sha, species=args.species)
            if read_selected(args.receipt, verification=verification) != runtime:
                raise ValueError('OrionGeno runtime changed during prediction')
            # The package-inventory subprocess above can outlive concurrent
            # input writers. Fence all used inputs and SQLite sidecars last.
            if verified_input_sha(genome, verification) != before:
                raise ValueError('OrionGeno genome changed before publication')
            verification.check()
            if selection_sha is not None:
                _, current_selection_sha = selection_record(args, genome, **selection_options)
                if current_selection_sha != selection_sha:
                    raise ValueError('OrionGeno model selection changed before publication')
            publish_pair(staged, output, receipt, result,
                         companions=[(path, companions[key]) for key, path in companion_sources.items()])
    return result


def collect(prediction_root, output):
    """Combine published per-species native tables without changing their sources."""
    rows, seen, hashes = [], set(), {}
    for table in sorted(Path(prediction_root).glob('*/candidate_inputs.tsv')):
        hashes[table] = file_sha(table)
        with table.open(newline='') as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            if reader.fieldnames != ['species', 'gff', 'receipt']:
                raise ValueError('Invalid native candidate table: '+str(table))
            for row in reader:
                if set(row) != set(reader.fieldnames) or any(not v for v in row.values()) or row['species'] in seen:
                    raise ValueError('Duplicate or incomplete native candidate row')
                seen.add(row['species'])
                rows.append([row['species'], *[str((table.parent/row[key]).resolve()) for key in ('gff', 'receipt')]])
    if not rows:
        raise ValueError('No completed OrionGeno candidate tables were found')
    output = Path(output)
    protect_outputs([output], [*hashes, *(p for row in rows for p in row[1:])])
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=output.parent, delete=False, newline='') as stream:
        temp = Path(stream.name)
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['species', 'gff', 'receipt'])
        writer.writerows(rows)
    try:
        if any(file_sha(path) != digest for path, digest in hashes.items()):
            raise ValueError('Native tables changed during collection')
        if not output.exists() or output.read_bytes() != temp.read_bytes():
            temp.replace(output)
    finally:
        temp.unlink(missing_ok=True)
    return dict(species=len(rows), table=str(output))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    collection = sub.add_parser('collect')
    collection.add_argument('--prediction-root', required=True)
    collection.add_argument('--output', required=True)
    prep = sub.add_parser('prepare')
    prep.add_argument('--cache-root', required=True)
    prep.add_argument('--receipt', required=True)
    prep.add_argument('--model-profile', choices=sorted(PROFILES), default='plants')
    prep.add_argument('--python', default='')
    prep.add_argument('--checkpoint', default='')
    prep.add_argument('--source-commit', default='')
    prep.add_argument('--source-ref', default='', help='Moving branch from container/source_branches.env; forwarded by the entrypoint')
    prep.add_argument('--no-download', action='store_true')
    pred = sub.add_parser('run')
    for key in ('receipt', 'genome', 'output', 'output-receipt'):
        pred.add_argument('--'+key, required=True)
    pred.add_argument('--threads', type=int, default=1)
    pred.add_argument('--batch-size', type=int, default=1)
    pred.add_argument('--output-repeat', type=int, choices=(0, 1), default=1)
    pred.add_argument('--output-utrs', type=int, choices=(0, 1), default=1)
    pred.add_argument('--selection-receipt', default='')
    pred.add_argument('--species', default='')
    args = parser.parse_args()
    result = collect(args.prediction_root, args.output) if args.command == 'collect' else prepare(args) if args.command == 'prepare' else run(args)
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
