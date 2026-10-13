"""Bounded workflow tests; fixture native output is not a GPU qualification."""

import gzip
import hashlib
import importlib.util
import json
import os
import re
import shlex
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow/support"
CORE = ROOT / "workflow/core/gg_genome_annotation_core.sh"
ENTRY = ROOT / "workflow/gg_genome_annotation_entrypoint.sh"


def source_function(command, *, env=None):
    environment = {**os.environ, "GG_CORE_SOURCE_ONLY": "1", **(env or {})}
    return subprocess.run(
        ["bash", "-c", f'source {shlex.quote(str(CORE))}; ' + command],
        env=environment, capture_output=True, text=True, timeout=20,
    )


@pytest.mark.parametrize("backend", ["none", "annevo"])
def test_explicit_backend_selection(backend):
    result = source_function("gg_validate_genome_annotation_backend", env={"genome_annotation_backend": backend})
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("backend", ["auto", "braker", "helixer", "ANNEVO", "annevo; echo unsafe"])
def test_unknown_backends_fail_before_runtime_preparation(backend):
    result = source_function("gg_validate_genome_annotation_backend", env={"genome_annotation_backend": backend})
    assert result.returncode != 0
    assert "must be none or annevo" in result.stderr


def test_verified_representative_view_cannot_be_replaced():
    result = source_function("gg_validate_genome_annotation_backend", env={
        "genome_annotation_backend": "annevo", "representative_inputs": "/verified/inputs.tsv",
    })
    assert result.returncode != 0
    assert "representative_inputs" in result.stderr


@pytest.mark.parametrize("name,value", [
    ("annevo_no_download", "2"), ("annevo_batch_size", "0"),
    ("annevo_sequence_threshold", "-1"), ("annevo_batch_size", "32; echo unsafe"),
    ("annevo_model_profile", "../../outside"),
])
def test_selected_backend_options_are_typed(name, value):
    result = source_function("gg_validate_genome_annotation_backend", env={
        "genome_annotation_backend": "annevo", name: value,
    })
    assert result.returncode != 0


def test_disabled_backend_does_not_validate_or_prepare_unused_lazy_resources(tmp_path):
    result = source_function("gg_validate_genome_annotation_backend", env={
        "genome_annotation_backend": "none", "annevo_no_download": "invalid",
        "annevo_runtime_manifest": str(tmp_path / "missing"),
        "annevo_batch_size": "invalid", "representative_inputs": "/verified/inputs.tsv",
    })
    assert result.returncode == 0, result.stderr
    assert list(tmp_path.iterdir()) == []


def test_all_backend_options_are_forwarded():
    result = subprocess.run(
        ["bash", "-c", 'source "$1"; gg_print_entrypoint_config_vars gg_genome_annotation_entrypoint.sh',
         "vars", str(SUPPORT / "gg_entrypoint_config_vars.sh")], capture_output=True, text=True, check=True,
    )
    expected = {"genome_annotation_backend", "annevo_model_profile", "annevo_runtime_manifest",
                "annevo_no_download", "annevo_batch_size", "annevo_sequence_threshold"}
    assert expected <= set(result.stdout.splitlines())
    assert 'genome_annotation_backend="none"' in ENTRY.read_text()


@pytest.mark.parametrize("backend", ["none", "annevo"])
def test_actual_input_discovery_selects_cds_or_genomes_only(backend, tmp_path):
    cds = tmp_path / "species_cds"
    genome = tmp_path / "species_genome"
    cds.mkdir()
    genome.mkdir()
    (cds / "Alpha_species.fa").write_text(">cds\nATG\n")
    (genome / "Beta_species.fa.gz").write_bytes(gzip.compress(b">chr\nATG\n"))
    (genome / ".hidden.fa").write_text("ignored")
    (genome / "notes.txt").write_text("ignored")
    alias = tmp_path / "linked_genomes"
    alias.symlink_to(genome, target_is_directory=True)
    text = CORE.read_text()
    start = text.index('infiles=()\nif [[ "${genome_annotation_backend}"')
    end = text.index('if [[ ! "${GG_ARRAY_TASK_ID}"', start)
    script = text[start:end] + 'printf "%s\\n" "${infiles[@]}"\n'
    result = subprocess.run(["bash", "-c", script], capture_output=True, text=True, env={
        **os.environ, "dir_sp_cds": str(cds), "dir_sp_genome": str(alias), "genome_annotation_backend": backend,
    })
    assert result.returncode == 0, result.stderr
    expected = cds / "Alpha_species.fa" if backend == "none" else alias / "Beta_species.fa.gz"
    assert result.stdout.strip() == str(expected)


def fake_native_runtime(path):
    # Only the runtime/provider boundary is stubbed. The production shell stage,
    # provenance, locks, formatter, CDS extraction and downstream resolution run.
    path.write_text('''import argparse, hashlib, json, os
from pathlib import Path
import sys

p=argparse.ArgumentParser();p.add_argument('command')
for key in ['cache-root','model-profile','runtime-manifest','receipt','genome','output','threads','batch-size','sequence-threshold','output-receipt']:
    p.add_argument('--'+key)
p.add_argument('--no-download',action='store_true');p.add_argument('--require-gpu',action='store_true')
a=p.parse_args()
log=Path(os.environ['ANNEVO_FIXTURE_LOG'])
with log.open('a') as h:h.write(json.dumps(vars(a))+'\\n')
if a.command=='prepare':
    r=Path(a.receipt)
    if a.no_download and not r.exists():raise SystemExit('fixture no-download cache is cold')
    r.parent.mkdir(parents=True,exist_ok=True)
    if not r.exists():r.write_text(json.dumps({'ready':True,'fixture_runtime_digest':'a'*64,'model_profile':a.model_profile},sort_keys=True))
    print(r.read_text())
elif a.command=='run':
    if os.environ.get('ANNEVO_FIXTURE_FAIL')=='1':raise SystemExit('fixture native failure')
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    gff=('##gff-version 3\\n'
         's\\tANNEVO\\tgene\\t1\\t306\\t.\\t+\\t.\\tID=gene1\\n'
         's\\tANNEVO\\tmRNA\\t1\\t306\\t.\\t+\\t.\\tID=gene1.t1;Parent=gene1\\n'
         's\\tANNEVO\\texon\\t1\\t306\\t.\\t+\\t.\\tID=exon1;Parent=gene1.t1\\n'
         's\\tANNEVO\\tCDS\\t1\\t306\\t.\\t+\\t0\\tID=cds1;Parent=gene1.t1\\n')
    out.write_text(gff)
    Path(a.output_receipt).write_text(json.dumps({'fixture_only':True,'genome_sha256':hashlib.sha256(Path(a.genome).read_bytes()).hexdigest(),'gff_sha256':hashlib.sha256(out.read_bytes()).hexdigest()},sort_keys=True))
''')


def run_selected_stage(tmp_path, *, policy="stop", fail=False, no_download=False, manifest=""):
    workspace = tmp_path / "workspace"
    support = tmp_path / "support"
    if not support.exists():
        support.mkdir()
        for path in SUPPORT.iterdir():
            if path.name != "annevo_runtime.py":
                (support / path.name).symlink_to(path, target_is_directory=path.is_dir())
        fake_native_runtime(support / "annevo_runtime.py")
    genome = workspace / "input/species_genome/Test_species.fa"
    genome.parent.mkdir(parents=True, exist_ok=True)
    if not genome.exists():
        genome.write_text(">s\n" + "ATG" + "AAA" * 100 + "TAA" + "C" * 15 + "\n")
    task = workspace / "output/tmp/test"
    task.mkdir(parents=True, exist_ok=True)
    log = tmp_path / "native_calls.jsonl"
    env = {**os.environ, "gg_support_dir": str(support), "gg_workspace_dir": str(workspace),
           "gg_workspace_output_dir": str(workspace / "output"),
           "gg_workspace_downloads_dir": str(workspace / "downloads"),
           "dir_sp_tmp": str(task), "sp_ub": "Test_species", "file_sp_genome": str(genome),
           "genome_annotation_backend": "annevo", "GG_TASK_CPUS": "4", "artifact_stale_policy": policy,
           "ANNEVO_FIXTURE_LOG": str(log), "ANNEVO_FIXTURE_FAIL": "1" if fail else "0",
           "annevo_no_download": "1" if no_download else "0", "annevo_runtime_manifest": manifest}
    command = ('source "$gg_support_dir/gg_util.sh"; gg_validate_genome_annotation_backend; '
               'gg_generate_annevo_annotation; check_species_cds_dir "$dir_sp_cds/."; gg_check_annevo_cds_directory; '
               'python "$gg_support_dir/cds_resolution.py" --cds "$file_sp_cds" --gff "$file_sp_gff" '
               '--genome "$file_sp_genome" --genetic-code "$genetic_code" '
               '--output-dir "$gg_workspace_output_dir/resolved"; '
               'printf "LAYOUT=%s\\t%s\\t%s\\n" "$file_sp_cds" "$file_sp_gff" "$file_sp_genome"')
    result = source_function(command, env=env)
    calls = [json.loads(line) for line in log.read_text().splitlines()] if log.exists() else []
    return result, calls, workspace


def test_real_formatter_and_resume_use_generated_inputs_without_mutating_curated_genome(tmp_path):
    first, calls, workspace = run_selected_stage(tmp_path)
    assert first.returncode == 0, first.stdout + first.stderr
    layout = re.search(r"^LAYOUT=(.+)$", first.stdout, re.M).group(1).split("\t")
    assert len(layout) == 3 and all(Path(path).is_file() for path in layout)
    assert all(str(workspace / "output/genome_prediction/annevo") in path for path in layout)
    assert [c['command'] for c in calls] == ['prepare', 'run']
    assert calls[-1]['threads'] == '4' and calls[-1]['batch_size'] == '32'
    assert calls[-1]['sequence_threshold'] == '25'
    assert gzip.open(layout[0], 'rt').read().startswith('>Test_species_')
    genome = workspace / "input/species_genome/Test_species.fa"
    before = genome.read_bytes()
    receipts = list((workspace / "output/genome_prediction/annevo").rglob('*.provenance.json'))
    assert len(receipts) == 2
    native_gff = workspace / "output/genome_prediction/annevo/Test_species/Test_species.gff3"
    stamps = {path: path.stat().st_mtime_ns for path in [native_gff, *map(Path, layout), *receipts]}
    again, calls, _ = run_selected_stage(tmp_path, no_download=True)
    assert again.returncode == 0, again.stdout + again.stderr
    assert sum(c['command'] == 'run' for c in calls) == 1
    assert 'Skipped: ANNEVO native genome prediction' in again.stdout
    assert genome.read_bytes() == before
    assert all(path.stat().st_mtime_ns == value for path, value in stamps.items())


def test_stale_genome_stops_and_failed_native_rebuild_preserves_previous_bundle(tmp_path):
    first, _, workspace = run_selected_stage(tmp_path)
    assert first.returncode == 0, first.stdout + first.stderr
    root = workspace / "output/genome_prediction/annevo/Test_species"
    preserved = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in root.rglob('*') if p.is_file()}
    genome = workspace / "input/species_genome/Test_species.fa"
    genome.write_text(genome.read_text() + '>other\nACGT\n')
    stale, calls, _ = run_selected_stage(tmp_path)
    assert stale.returncode == 3, stale.stdout + stale.stderr
    assert 'Stale artifact' in stale.stdout + stale.stderr
    assert sum(c['command'] == 'run' for c in calls) == 1
    failed, _, _ = run_selected_stage(tmp_path, policy='rebuild', fail=True)
    assert failed.returncode != 0 and 'fixture native failure' in failed.stderr
    assert all(hashlib.sha256(p.read_bytes()).hexdigest() == digest for p, digest in preserved.items())
    rebuilt, calls, _ = run_selected_stage(tmp_path, policy='rebuild')
    assert rebuilt.returncode == 0, rebuilt.stdout + rebuilt.stderr
    assert len(list((root / 'formatted').iterdir())) == 2
    # The previous complete formatted generation is retained.
    assert all(hashlib.sha256(p.read_bytes()).hexdigest() == digest for p, digest in preserved.items()
               if '/formatted/' in str(p))


def test_no_download_requires_a_warm_verified_cache(tmp_path):
    result, calls, workspace = run_selected_stage(tmp_path, no_download=True)
    assert result.returncode != 0
    assert 'cache is cold' in result.stderr
    assert [c['command'] for c in calls] == ['prepare']
    assert not list((workspace / 'output').rglob('*.gff3'))


def test_core_native_argv_is_accepted_by_actual_runtime_cli(tmp_path):
    result, calls, _ = run_selected_stage(tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    spec = importlib.util.spec_from_file_location('core_api_fixture', ROOT / 'workflow/tests/test_annevo_runtime.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for call in calls:
        argv = [call['command']]
        for name, value in call.items():
            if name != 'command' and value is not None and value is not False:
                argv += ['--' + name.replace('_', '-')]
                if value is not True:
                    argv.append(value)
        parsed = module.runtime.parser().parse_args(argv)
        assert parsed.command == call['command']
        if parsed.command == 'run':
            assert parsed.threads == 4 and parsed.batch_size == 32 and parsed.sequence_threshold == 25


def test_actual_runtime_warm_no_download_core_resume_requires_no_native_execution(tmp_path, monkeypatch):
    first, _, workspace = run_selected_stage(tmp_path)
    assert first.returncode == 0, first.stdout + first.stderr
    # Install a small offline, explicitly synthetic runtime with the existing
    # loader test fixture. Afterwards the actual CLI verifies its warm cache;
    # its synthetic interpreter must never execute on prediction resume.
    spec = importlib.util.spec_from_file_location('warm_core_fixture', ROOT / 'workflow/tests/test_annevo_runtime.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    cache_root = workspace / 'downloads/annotation/annevo'
    selected = cache_root / 'selected/Magnoliopsida.json'
    assets = tmp_path / 'offline_assets'
    assets.mkdir()
    with monkeypatch.context() as setup:
        manifest, _ = module.bundle(assets, setup)
        module.cache.prepare(cache_root, 'Magnoliopsida', str(manifest), receipt_path=str(selected))
    support = tmp_path / 'support'
    (support / 'annevo_runtime.py').unlink()
    (support / 'annevo_runtime.py').write_bytes((SUPPORT / 'annevo_runtime.py').read_bytes())
    prediction = workspace / 'output/genome_prediction/annevo/Test_species'
    # Bind the pre-existing fixture prediction to the actual loader/receipt
    # contract, so the core must reuse it rather than rerun a native predictor.
    arguments = [
        '--input', f'genome={workspace / "input/species_genome/Test_species.fa"}',
        '--input', f'runtime_receipt={selected}', '--input', f'runner={support / "annevo_runtime.py"}',
        '--output', f'gff={prediction / "Test_species.gff3"}',
        '--output', f'native_receipt={prediction / "native_prediction.json"}',
        '--parameter', 'model_profile=Magnoliopsida', '--parameter', 'batch_size=32',
        '--parameter', 'sequence_threshold=25', '--parameter', 'threads=4', '--parameter', 'require_gpu=1',
    ]
    shell = ('source "$gg_support_dir/gg_util.sh"; '
             'gg_artifact_contract_init proof genome_annotation_annevo_prediction Test_species "$manifest"; '
             'gg_artifact_record "${proof[@]}" "$@"')
    bound = subprocess.run(['bash', '-c', shell, 'bind-fixture', *arguments], env={
        **os.environ, 'gg_support_dir': str(support), 'gg_workspace_dir': str(workspace),
        'gg_workspace_output_dir': str(workspace / 'output'),
        'manifest': str(prediction / 'prediction.provenance.json'),
    }, capture_output=True, text=True, timeout=20)
    assert bound.returncode == 0, bound.stdout + bound.stderr
    before = selected.read_bytes(), selected.stat().st_mtime_ns
    resumed, _, _ = run_selected_stage(tmp_path, no_download=True, manifest=str(manifest))
    assert resumed.returncode == 0, resumed.stdout + resumed.stderr
    assert 'Skipped: ANNEVO native genome prediction' in resumed.stdout
    assert (selected.read_bytes(), selected.stat().st_mtime_ns) == before
    assert 'LAYOUT=' in resumed.stdout


@pytest.mark.parametrize('backend,devices,success', [('none', '', True), ('annevo', '', False),
                                                    ('annevo', '0', True), ('annevo', 'all', False)])
@pytest.mark.parametrize('representation', ['array', 'legacy_string'])
def test_entrypoint_requests_gpu_only_for_selected_backend(backend, devices, success, representation):
    text = ENTRY.read_text()
    start = text.index('if [[ "${genome_annotation_backend}" == "annevo" ]]')
    end = text.index('gg_entrypoint_activate_container_runtime', start)
    declaration = ('singularity_command=(fixture exec)' if representation == 'array'
                   else 'singularity_command="fixture exec"')
    command = declaration + ';\n' + text[start:end] + '\nprintf "%s\\n" "${singularity_command[@]}"'
    result = subprocess.run(['bash', '-c', command], capture_output=True, text=True,
                            env={**os.environ, 'genome_annotation_backend': backend,
                                 'CUDA_VISIBLE_DEVICES': devices, 'GG_CONTAINER_DOCKER_GPU_DEVICES': ''})
    assert (result.returncode == 0) == success, result.stderr
    if success:
        expected = ['fixture', 'exec'] if representation == 'array' else ['fixture exec']
        if backend == 'annevo':
            if representation == 'array':
                expected.append('--nv')
            else:
                expected[0] += ' --nv'
        assert result.stdout.splitlines() == expected


def test_completed_cds_validation_stamp_reuses_only_unchanged_input_content(tmp_path):
    cds = tmp_path / 'cds'
    cds.mkdir()
    fasta = cds / 'Test_species.fa'
    fasta.write_text('>Test_species_gene1\nATGAAATAA\n')
    text = CORE.read_text()
    start = text.index('species_cds_validation_signature=$(\n')
    end = text.index('# Preserve previous manifest-bound results', start)
    stage = text[start:end]
    log = tmp_path / 'validation_calls'
    shell = ('set -euo pipefail\nsource "$gg_support_dir/gg_util.sh"\n'
             'check_species_cds_dir() { printf "checked\\n" >> "$validation_log"; }\n'
             'check_if_species_files_unique() { :; }\n'
             'species_cds_validation_inputs=("$validation_input")\n' + stage)
    env = {**os.environ, 'gg_support_dir': str(SUPPORT), 'dir_sp_cds': str(cds),
           'dir_tmp': str(tmp_path), 'validation_input': str(fasta), 'validation_log': str(log),
           'genome_annotation_backend': 'none'}
    for _ in range(2):
        result = subprocess.run(['bash', '-c', shell], env=env, capture_output=True, text=True, timeout=20)
        assert result.returncode == 0, result.stdout + result.stderr
    assert log.read_text().splitlines() == ['checked']
    stamps = list(tmp_path.glob('species_cds_validation.*.ok'))
    assert len(stamps) == 1 and len(stamps[0].read_text().strip()) == 64
    fasta.write_text(fasta.read_text() + '>Test_species_gene2\nATGAAGTAA\n')
    result = subprocess.run(['bash', '-c', shell], env=env, capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stdout + result.stderr
    assert log.read_text().splitlines() == ['checked', 'checked']
    assert len(list(tmp_path.glob('species_cds_validation.*.ok'))) == 2


def test_generated_cds_scope_keeps_audits_but_rejects_additional_fasta(tmp_path):
    cds = tmp_path / 'species_cds'
    cds.mkdir()
    selected = cds / 'Test_species.fa'
    selected.write_text('>Test_species_gene1\nATGAAATAA\n')
    (cds / 'Test_species.fa.cds-normalisation.json').write_text('{}')
    command = 'source "$gg_support_dir/gg_util.sh"; gg_check_annevo_cds_directory'
    environment = {'gg_support_dir': str(SUPPORT), 'dir_sp_cds': str(cds), 'file_sp_cds': str(selected)}
    first = source_function(command, env=environment)
    assert first.returncode == 0, first.stderr
    extra = cds / 'Other_species.fa'
    extra.write_text('>Other_species_gene1\nATGAAATAA\n')
    failed = source_function(command, env=environment)
    assert failed.returncode != 0
    assert 'exactly the selected species FASTA' in failed.stderr
    extra.unlink()
    selected.unlink()
    missing = source_function(command, env=environment)
    assert missing.returncode != 0
