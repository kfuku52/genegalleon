"""Real subprocess scheduling/publication tests; fixtures explicitly substitute native inference."""
import argparse
import csv
import importlib
import json
import os
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / 'support'
sys.path.insert(0, str(SUPPORT))
pipeline = importlib.import_module('annevo_pipeline')


@pytest.fixture
def batch(tmp_path, monkeypatch):
    source = tmp_path / 'native'
    source.mkdir()
    events = tmp_path / 'events.jsonl'
    common = '''import argparse,json,time,os
from pathlib import Path
p=argparse.ArgumentParser()
for n in ['-g','-p','-o','-t','-m','-l','-s','--batch_size','--num_workers','--min_intron_length','--min_prot_length']:p.add_argument(n)
p.add_argument('--show_log',action='store_true');p.add_argument('--overlap_pred',action='store_true');p.add_argument('--comp',action='store_true')
a=p.parse_args();name=Path(a.g).stem
def event(stage):
 with open(EVENTS,'a') as f:f.write(json.dumps({'event':stage,'name':name,'time':time.monotonic(),'cuda':os.environ.get('CUDA_VISIBLE_DEVICES'),'jit':os.environ.get('NUMBA_CACHE_DIR'),'compression':a.comp})+'\\n')
'''.replace('EVENTS', repr(str(events)))
    (source / 'prediction.py').write_text("length_config={'Magnoliopsida':[30720,5120,32]}\n" + common + '''
event('prediction-start');time.sleep(.12)
Path(a.p).write_text('probabilities');event('prediction-end')
''')
    (source / 'decoding.py').write_text(common + '''
event('decode-start');time.sleep(.35)
if name=='fail':raise SystemExit(7)
gff='seq\\tANNEVO\\tgene\\t1\\t9\\t.\\t+\\t.\\tID=g1\\nseq\\tANNEVO\\tmRNA\\t1\\t9\\t.\\t+\\t.\\tID=t1;Parent=g1\\nseq\\tANNEVO\\tCDS\\t1\\t9\\t.\\t+\\t0\\tID=c1;Parent=t1\\n'
Path(a.o).write_text(gff);event('decode-end')
''')
    selected = tmp_path / 'selected.json'
    selected.write_text('{}')
    value = dict(model_profile='Magnoliopsida', python=sys.executable,
                 native_entrypoint=str(source / 'annotation.py'), generation=str(source),
                 model_path=str(source / 'model.pt'), runtime_digest='a'*64, cache_key='b'*64,
                 source_commit='c'*40, model_sha256='d'*64)
    monkeypatch.setattr(pipeline.native, 'read_selected', lambda p: value)
    monkeypatch.setattr(pipeline.native, 'probe_gpu', lambda *a: {'gpu_verified': True})
    monkeypatch.setattr(pipeline, 'runtime_environment', lambda r: dict(__import__('os').environ))
    inputs = tmp_path / 'inputs.tsv'

    def make(names=('one', 'two', 'three')):
        with inputs.open('w', newline='') as stream:
            writer = csv.writer(stream, delimiter='\t')
            writer.writerow(['species', 'genome', 'receipt'])
            for name in names:
                genome = tmp_path / (name + '.fa')
                genome.write_text('>seq\nATGAAATAA\n')
                writer.writerow([name, genome, selected])
        return argparse.Namespace(inputs=inputs, output=tmp_path/'out', decoder_workers=2,
            decode_threads=1, predict_threads=1, batch_size=32, sequence_threshold=25,
            tmp_root=tmp_path, max_tmp_gib=1)
    return make, events, selected


def test_one_gpu_two_cpu_decoders_overlap_and_export_actual_native_receipts(batch):
    make, events, _ = batch
    args = make()
    result = pipeline.run(args)
    rows = [json.loads(line) for line in events.read_text().splitlines()]
    gpu, cpu, max_gpu, max_cpu = 0, 0, 0, 0
    for row in sorted(rows, key=lambda row: row['time']):
        gpu += {'prediction-start': 1, 'prediction-end': -1}.get(row['event'], 0)
        cpu += {'decode-start': 1, 'decode-end': -1}.get(row['event'], 0)
        max_gpu, max_cpu = max(max_gpu, gpu), max(max_cpu, cpu)
        if row['event'] == 'decode-start':
            assert row['cuda'] == ''
    assert max_gpu == 1 and max_cpu == 2 and gpu == cpu == 0
    assert result['max_inflight_species'] == 3
    assert result['peak_reserved_bytes'] <= result['max_tmp_bytes']
    with (args.output/'candidate_inputs.tsv').open() as stream:
        exported = list(csv.DictReader(stream, delimiter='\t'))
    assert [row['species'] for row in exported] == ['one', 'two', 'three']
    import annevo_gene_models
    for row in exported:
        identity = annevo_gene_models.prediction_identity(row['gff'], row['receipt'],
                     pipeline.file_sha(args.tmp_root/(row['species']+'.fa')))
        assert identity['genome_match'] == 'raw_identity'
    assert not list(args.tmp_root.glob('gg-ap-*'))
    before = events.read_bytes()
    assert all(row['reused'] for row in pipeline.run(args)['species'].values())
    assert events.read_bytes() == before


def test_decode_failure_retains_logs_and_does_not_publish_campaign_or_failed_pair(batch):
    make, _, _ = batch
    args = make(('fail', 'two'))
    with pytest.raises(RuntimeError, match='decoding failed'):
        pipeline.run(args)
    assert not (args.output/'candidate_inputs.tsv').exists()
    assert not (args.output/'fail/native_prediction.json').exists()
    assert (args.output/'fail/logs/decoding.stderr.txt').exists()
    assert json.loads((args.output/'pipeline.failed.json').read_text())['status'] == 'failed'
    assert not list(args.tmp_root.glob('gg-ap-*'))


def test_largest_first_keeps_export_order_and_reuses_unchanged_pairs(batch):
    make, events, _ = batch
    args = make()
    for name, copies in [('one', 2), ('two', 4), ('three', 1)]:
        (args.tmp_root/(name+'.fa')).write_text('>seq\n'+'ATGAAATAA'*copies+'\n')
    args.order = 'largest-first'
    result = pipeline.run(args)
    assert result['execution_order']==['two', 'one', 'three']
    trace=[json.loads(line) for line in events.read_text().splitlines()]
    assert [row['name'] for row in trace if row['event']=='prediction-start']==result['execution_order']
    with (args.output/'candidate_inputs.tsv').open() as stream:
        assert [row['species'] for row in csv.DictReader(stream,delimiter='\t')]==['one','two','three']
    before=events.read_bytes()
    args.order='input'
    assert all(row['reused'] for row in pipeline.run(args)['species'].values())
    assert events.read_bytes()==before


def test_compression_is_forwarded_recorded_and_invalidates_an_existing_plain_pair(batch):
    make, events, _ = batch
    args=make(('one',))
    args.compression='lzf'
    pipeline.run(args)
    trace=[json.loads(line) for line in events.read_text().splitlines()]
    assert all(row['compression'] for row in trace if row['event'].startswith('prediction'))
    receipt=json.loads((args.output/'one/native_prediction.json').read_text())
    assert receipt['parameters']['compression']=='lzf' and '--comp' in receipt['argv']
    args.compression='none'
    with pytest.raises(ValueError,match='Stale pipeline'):
        pipeline.run(args)


def test_shared_jit_cache_is_bound_and_persists_without_touching_runtime(batch):
    make, events, _ = batch
    args=make()
    args.jit_cache_root=args.tmp_root/'shared-jit'
    pipeline.run(args)
    trace=[json.loads(line) for line in events.read_text().splitlines()]
    paths={row['jit'] for row in trace if row['event']=='decode-start'}
    assert len(paths)==1
    directory=Path(paths.pop())
    assert directory.is_dir() and (directory/'cache_identity.json').is_file()
    assert not directory.is_relative_to(args.tmp_root/'native')
    for name in ['one','two','three']:
        receipt=json.loads((args.output/name/'native_prediction.json').read_text())
        assert receipt['pipeline']['jit_cache_directory']==str(directory)
        assert receipt['pipeline']['jit_cache_key']==directory.name
    assert not list(args.tmp_root.glob('gg-ap-*'))


def test_schedule_rejects_genome_changes_before_native_prediction(batch,monkeypatch):
    make,events,_=batch
    args=make(('one',))
    args.order='largest-first'
    original=pipeline.genome_bases
    def changing(path):
        value=original(path)
        Path(path).write_text('>seq\nATGCCCTAA\n')
        return value
    monkeypatch.setattr(pipeline,'genome_bases',changing)
    with pytest.raises(ValueError,match='changed while scheduling'):
        pipeline.run(args)
    assert not events.exists()


@pytest.mark.parametrize(('option','value'), [('order','unknown'),('compression','lossy')])
def test_unknown_performance_options_reject_before_prediction(batch,option,value):
    make,events,_=batch
    args=make()
    setattr(args,option,value)
    with pytest.raises(ValueError,match='Unknown ANNEVO'):
        pipeline.run(args)
    assert not events.exists()


def test_storage_budget_rejects_before_predicting(batch):
    make, events, _ = batch
    args = make(('one',))
    args.max_tmp_gib = .001
    with pytest.raises(ValueError, match='storage budget'):
        pipeline.run(args)
    assert not events.exists()


def test_small_storage_budget_backpressures_without_losing_species(batch):
    make, events, _ = batch
    args = make()
    args.max_tmp_gib = .02  # one reserved species fits; two do not
    result = pipeline.run(args)
    assert result['max_inflight_species'] == 1
    assert len(result['species']) == 3
    rows = [json.loads(line) for line in events.read_text().splitlines()]
    active = 0
    for row in sorted(rows, key=lambda r:r['time']):
        if row['event']=='decode-start':
            active += 1
        elif row['event']=='decode-end':
            active -= 1
        elif row['event']=='prediction-start':
            assert active == 0


def test_cpu_budget_rejects_allocation_overcommit(batch, monkeypatch):
    make, events, _ = batch
    monkeypatch.setenv('SLURM_CPUS_PER_TASK', '2')
    with pytest.raises(ValueError, match='CPU budget'):
        pipeline.run(make())
    assert not events.exists()


def test_stale_genome_and_receipt_pairs_are_not_reused(batch):
    make, _, _ = batch
    args = make(('one',))
    pipeline.run(args)
    (args.tmp_root/'one.fa').write_text('>seq\nATGCCCTAA\n')
    with pytest.raises(ValueError, match='Stale pipeline'):
        pipeline.run(args)


def test_table_change_during_execution_is_rejected(batch, monkeypatch):
    make, _, _ = batch
    args = make(('one',))
    original = pipeline.execute
    def change(*a, **kw):
        result = original(*a, **kw)
        args.inputs.write_text(args.inputs.read_text()+'\n')
        return result
    monkeypatch.setattr(pipeline, 'execute', change)
    with pytest.raises(RuntimeError, match='input table changed'):
        pipeline.run(args)
    assert not (args.output/'candidate_inputs.tsv').exists()


def test_redirected_pair_is_rejected_before_native_prediction(batch):
    make, events, _ = batch
    args = make(('one',))
    directory = args.output/'one'
    directory.mkdir(parents=True)
    (directory/'one.gff3').symlink_to(directory/'missing')
    with pytest.raises(ValueError, match='Redirected'):
        pipeline.run(args)
    assert not events.exists()


@pytest.mark.parametrize('input_role', ['genome', 'later_genome', 'receipt', 'table', 'selection_evidence', 'model'])
@pytest.mark.parametrize('log_name', ['prediction.stdout.txt', 'prediction.stderr.txt', 'decoding.stdout.txt', 'decoding.stderr.txt'])
def test_logs_cannot_truncate_batch_inputs_or_selection_evidence(batch, monkeypatch, input_role, log_name):
    make, events, selected = batch
    args = make()
    protected = {'genome': args.tmp_root/'one.fa', 'later_genome': args.tmp_root/'two.fa',
                 'receipt': selected, 'table': args.inputs}.get(input_role)
    if input_role == 'selection_evidence':
        protected = args.tmp_root/'taxids.tsv'
        protected.write_text('species\ttaxid\none\t3398\n')
        choice = args.tmp_root/'selection.json'
        choice.write_text(json.dumps({'evidence': [{'path': str(protected)}]}))
        rows = list(csv.DictReader(args.inputs.open(), delimiter='\t'))
        with args.inputs.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=['species', 'genome', 'receipt', 'selection_receipt'], delimiter='\t')
            writer.writeheader()
            for row in rows:
                writer.writerow({**row, 'selection_receipt': str(choice) if row['species']=='two' else ''})
    elif input_role == 'model':
        protected = args.tmp_root/'native/model.pt'
        protected.write_text('model weights')
    logs = args.output/'one/logs'
    logs.mkdir(parents=True)
    os.link(protected, logs/log_name)
    before = protected.read_bytes()
    with pytest.raises(ValueError, match='outputs must not overwrite inputs'):
        pipeline.run(args)
    assert protected.read_bytes() == before
    assert not events.exists()


def test_log_directory_cannot_redirect_writes_into_runtime(batch):
    make, events, _ = batch
    args = make(('one',))
    directory = args.output/'one'
    directory.mkdir(parents=True)
    (directory/'logs').symlink_to(args.tmp_root/'native', target_is_directory=True)
    with pytest.raises(ValueError, match='installed runtime'):
        pipeline.run(args)
    assert not events.exists()


def test_failure_receipt_cannot_replace_selection_evidence(batch):
    make, events, _ = batch
    args = make(('fail',))
    args.output.mkdir()
    evidence = args.output/'pipeline.failed.json'
    evidence.write_text('selection evidence')
    choice = args.tmp_root/'selection.json'
    choice.write_text(json.dumps({'evidence': [{'path': str(evidence)}]}))
    args.inputs.write_text(args.inputs.read_text().replace('genome\treceipt\n', 'genome\treceipt\tselection_receipt\n').rstrip('\n')+'\t'+str(choice)+'\n')
    before = evidence.read_bytes()
    with pytest.raises(ValueError, match='outputs must not overwrite inputs'):
        pipeline.run(args)
    assert evidence.read_bytes()==before and not events.exists()


@pytest.mark.parametrize('budget', [float('nan'), float('inf')])
def test_nonfinite_storage_budget_rejected(batch, budget):
    make, events, _ = batch
    args = make()
    args.max_tmp_gib = budget
    with pytest.raises(ValueError, match='bounds must be positive'):
        pipeline.run(args)
    assert not events.exists()
