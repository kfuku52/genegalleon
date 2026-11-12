"""Held-out exon/path truth tested with real miniprot on temporary genomes."""

import gzip
import hashlib
import json
import random
import sqlite3
import subprocess
import sys
from importlib import import_module
from pathlib import Path

import pytest
from Bio.Data import CodonTable
from Bio.Seq import Seq

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))

refinement = import_module("gene_model_refinement")
gff2genestat = import_module("gff2genestat")
fasta_records = import_module("fasta_sequence_store").fasta_records
write_tsv = import_module("test_gene_model_refinement").write_tsv

SCRIPT = Path(refinement.__file__)
GFF_COLS = ["sequence", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"]
OUT_COLS = ["gene_id", "feature_size", "num_intron", "gff_transcript_id", "feature_blocks"]


def test_refinement_cli_help_uses_real_runtime_without_outputs(tmp_path):
    result = subprocess.run([sys.executable, str(SCRIPT), "--help"], cwd=tmp_path,
                            text=True, capture_output=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "usage:" in result.stdout.lower() and "verify-inputs" in result.stdout
    assert list(tmp_path.iterdir()) == []


def test_refinement_preserves_missing_models_from_completed_real_rescue(tmp_path, monkeypatch):
    fixture = import_module('test_rescue_gene_models_runtime')
    data = fixture.hidden_models.__wrapped__(tmp_path)
    anchor_root, names, truth = fixture.make_plan(data)
    before = {path: path.read_bytes() for role in ('cds', 'gff', 'genome') for path in (tmp_path / role).iterdir()}
    fixture.cli('run', '--output', anchor_root, '--cpus', 2)
    models = json.loads((anchor_root / 'rescued' / names[0] / 'models.json').read_text())
    accepted = [model for model in models if model['status'] == 'accepted']
    assert len(accepted) == 1 and accepted[0]['sequence'] == truth[8]
    root = tmp_path / 'refinement'
    producer_plan = (anchor_root / 'plan.json').read_bytes()
    identities = refinement.rescue.identities()
    monkeypatch.setattr(refinement.rescue, 'identities', lambda: {
        **identities, 'genome_index_implementation': 'downstream-catalogue-update'})
    with pytest.raises(ValueError, match='Rescue schema/tools changed'):
        refinement.rescue.load(anchor_root)
    # Published producers remain readable after a downstream-only helper update.
    refinement.plan(root, rescue_output=anchor_root, mode='off')
    value = refinement.load(root)
    source = value['request']['sources'][names[0]]
    assert Path(source['fasta']).is_relative_to(anchor_root / 'augmented')
    assert Path(source['gff']).is_relative_to(anchor_root / 'augmented')
    cli('run', '--output', root, '--cpus', 2)
    cli('qc', '--output', root)
    assert all(refinement.inspect_stages(root, refinement.load(root)).values())
    assert (anchor_root / 'plan.json').read_bytes() == producer_plan
    records = {identifier: sequence for identifier, _header, sequence in
               fasta_records(root / 'effective/species_cds' / (names[0] + '.fa'))}
    assert truth[8] in records.values()
    assert len(records) == 16  # Fifteen original annotations plus one missing model.
    assert accepted[0]['model_id'] in (root / 'effective/full_annotation' / (names[0] + '.gff3')).read_text()
    assert all(path.read_bytes() == content for path, content in before.items())
    augmented = anchor_root / 'augmented/species_cds' / (names[0] + '.rescue.cds.fa')
    augmented.write_text(augmented.read_text() + '\n')
    with pytest.raises(ValueError, match='Frozen refinement input changed'):
        refinement.load(root)
    with pytest.raises(ValueError, match='Augmented rescue inputs incomplete or corrupt'):
        refinement.plan(tmp_path / 'new-refinement', rescue_output=anchor_root, mode='off')


@pytest.mark.parametrize(('fragment_count', 'manual'), [(2, False), (3, False), (3, True)])
def test_split_model_merge_from_real_rescue(tmp_path, fragment_count, manual):
    fixture = import_module('test_rescue_gene_models_runtime')
    data = fixture.hidden_models.__wrapped__(tmp_path)
    _, names, truth = data
    target, sequence = names[0], truth[8]
    cds = tmp_path / 'cds' / (target + '.cds.fa')
    gff = tmp_path / 'gff' / (target + '.gff3')
    genome = tmp_path / 'genome' / (target + '.genome.fa')
    start = 8 * (len(sequence) + 60)
    if fragment_count == 2:
        boundaries = [0, 240, len(sequence)]
    else:
        # At codon 106 the last fragment has its own in-frame ATG and stop.
        boundaries = [0, 159, 318 if manual else 321, len(sequence)]
    fragments = [sequence[a:b] for a, b in zip(boundaries, boundaries[1:], strict=False)]
    intron = 'GT' + 'C' * (40 if fragment_count == 2 else 26) + 'AG'
    span = str(Seq(intron.join(fragments)).reverse_complement())
    end = start + len(span)
    dna = genome.read_text().splitlines()[1]
    genome.write_text('>chr1\n' + dna[:start] + span + dna[end:] + '\n')
    split = []
    exon_end = end
    for index, fragment in enumerate(fragments):
        exon_start = exon_end - len(fragment)
        split.append(('g8' + 'abc'[index], [[exon_start, exon_end, 0]], fragment))
        exon_end = exon_start - len(intron)
    source_genes = [gene for gene, _, _ in split]
    with gff.open('a') as annotation_file, cds.open('a') as fasta:
        for gene, blocks, fragment in split:
            annotation_file.write(annotation('chr1', gene, gene + '.t1', blocks, '-'))
            fasta.write(f'>{target}_{gene}\n{fragment}\n')
    before = {path: path.read_bytes() for path in (cds, gff, genome)}

    anchor, _, _ = fixture.make_plan(data, model_storage='compact')
    fixture.cli('run', '--output', anchor, '--cpus', 2)
    worker = anchor / 'rescued' / target
    proposals = json.loads((worker / 'split_model_candidates.json').read_text())
    matches = [proposal for proposal in proposals if proposal['cds_sha256'] == hashlib.sha256(sequence.encode()).hexdigest()]
    assert len(matches) == 1, proposals
    review = matches[0]
    assert review['strand'] == '-'
    assert review['cds'] == [blocks[0] for _, blocks, _ in split]
    assert [row['gene_id'] for row in review['owner_overlap']] == source_genes
    if fragment_count == 3:
        evidence = [row['splice_revalidation'] for row in review['support'] if 'splice_revalidation' in row]
        assert evidence and all(row['original_cds'] == [[start, end, 0]] for row in evidence)
    identifier = review['candidate_id']
    original_ids = {target + '_' + gene for gene in source_genes}
    augmented = anchor / 'augmented' / 'species_cds' / (target + '.rescue.cds.fa')
    assert original_ids <= {gene for gene, _, _ in fasta_records(augmented)}
    assert identifier not in {gene for gene, _, _ in fasta_records(augmented)}
    producer_receipt = (worker / 'receipt.json').read_bytes()
    approvals = tmp_path / 'merge.tsv'
    write_tsv(approvals, ['species', 'candidate_id'], [{'species': target, 'candidate_id': identifier}])
    with pytest.raises(ValueError, match='--mode conservative'):
        refinement.plan(tmp_path / 'audit', rescue_output=anchor, mode='audit', split_model_merges=approvals)

    output = tmp_path / 'refinement'
    approval_args = ('--split-model-merges', approvals) if manual else ()
    cli('plan', '--rescue-output', anchor, '--output', output, *approval_args)
    cli('run', '--output', output, '--cpus', 2)
    cli('qc', '--output', output)
    effective = output / 'effective'
    inputs = next(row for row in refinement.verify_inputs(effective / 'inputs.tsv') if row['species'] == target)
    for role in ('cds', 'analysis_cds', 'protein'):
        records = {gene: value for gene, _, value in fasta_records(Path(inputs[role]))}
        assert records[identifier] == (str(Seq(sequence).translate()).removesuffix('*') if role == 'protein' else sequence)
        assert not original_ids & records.keys()
    for role in ('gff', 'analysis_gff'):
        text = Path(inputs[role]).read_text()
        assert f'ID={identifier}.t1;Parent={identifier}' in text
        assert not any(f'ID={gene};' in text or f'ID={gene}.t1;' in text for gene in source_genes)
    archive = (effective / 'full_annotation' / (target + '.gff3')).read_text()
    assert all(f'ID={gene}.t1;Parent={gene}' in archive for gene in [*source_genes, identifier])
    assert (effective / 'source_annotation' / (target + '.gff3')).read_bytes() == (
        anchor / 'augmented' / 'species_gff' / (target + '.rescue.gff3')).read_bytes()
    mapping = {row['gene_id']: row for row in refinement.read_table(effective / 'representative_map.tsv')
               if row['species'] == target}
    assert mapping[identifier]['status'] == 'merged_split_model'
    assert not original_ids & mapping.keys()
    audit = json.loads((effective / 'split_model_merges.json').read_text())
    assert len(audit) == 1 and set(audit[0]['replaced_gene_ids']) == original_ids
    if not manual:
        assert audit[0]['adoption'] == 'automatic'
    review_data = import_module('plot_gene_model_refinement').collect(output)
    assert review_data['split_model_merges'] == audit
    assert all(path.read_bytes() == content for path, content in before.items())
    assert (worker / 'receipt.json').read_bytes() == producer_receipt
    if manual:
        approvals.write_text(approvals.read_text() + '\n')
        with pytest.raises(ValueError, match='Frozen refinement input changed'):
            refinement.load(output)


@pytest.mark.parametrize('case', ['partial', 'intact', 'reciprocal', 'manual'])
def test_single_locus_repairs_from_real_rescue(tmp_path, case, monkeypatch):
    manual = case == 'manual'
    fixture = import_module('test_rescue_gene_models_runtime')
    data = fixture.hidden_models.__wrapped__(tmp_path)
    _, names, truth = data
    coding_target, phase_target = names[:2]
    coding_cds = tmp_path / 'cds' / (coding_target + '.cds.fa')
    coding_gff = tmp_path / 'gff' / (coding_target + '.gff3')
    phase_gff = tmp_path / 'gff' / (phase_target + '.gff3')
    phase_genome = tmp_path / 'genome' / (phase_target + '.genome.fa')
    period = len(truth[8]) + 60
    start = 8 * period
    contig = 'chr1'
    if case == 'reciprocal':
        coding_genome = tmp_path / 'genome' / (coding_target + '.genome.fa')
        dna = coding_genome.read_text().splitlines()[1]
        dna = dna[:start] + 'N' * len(truth[8]) + dna[start + len(truth[8]):]
        coding_genome.write_text('>chr1\n' + dna + '\n>isolated\n' + truth[8] + '\n')
        contig, start = 'isolated', 0
    if case in {'intact', 'reciprocal'}:
        source_cds = truth[8][60:]
        assert source_cds.startswith('ATG')
        source_blocks = [[start + 60, start + len(truth[8]), 0]]
    else:
        source_cds = truth[8][:-60]
        source_blocks = [[start, start + len(source_cds), 0]]
    with coding_cds.open('a') as handle:
        handle.write(f'>{coding_target}_g8\n{source_cds}\n')
    with coding_gff.open('a') as handle:
        handle.write(annotation(contig, 'g8', 'g8.t1', source_blocks))
    phase_start = 11 * period
    first_exon = 241
    intron = 'GT' + 'C' * 40 + 'AG'
    phase_junction = phase_start + first_exon
    phase_end = phase_start + len(truth[11]) + len(intron)
    dna = phase_genome.read_text().splitlines()[1]
    spliced = truth[11][:first_exon] + intron + truth[11][first_exon:]
    phase_genome.write_text('>chr1\n' + dna[:phase_start] + spliced + dna[phase_end:] + '\n')
    original = [line for line in phase_gff.read_text().splitlines(keepends=True)
                if not any(token in line for token in ('ID=g11;', 'ID=g11\n', 'ID=g11.t1;', 'Parent=g11.t1'))]
    original.append(annotation('chr1', 'g11', 'g11.t1', [[phase_start, phase_junction, 0],
                                                      [phase_junction + len(intron), phase_end, 0]]))
    phase_gff.write_text(''.join(original))
    before = {path: path.read_bytes() for path in (coding_cds, coding_gff, phase_gff, phase_genome)}
    expected = {coding_target: (8, 'coding_path'), phase_target: (11, 'phase_only')}
    rescue_args = ()
    if manual:
        targets = tmp_path / 'targets.tsv'
        write_tsv(targets, ['species', 'gene_id', 'donor', 'query'],
                  [{'species': target, 'gene_id': f'g{i}', 'donor': name, 'query': f'{name}_g{i}'}
                   for target, (i, _) in expected.items() for name in names[2:4]])
        rescue_args = ('--model-repair-targets', targets)
    anchor, _, _ = fixture.make_plan(data, model_storage='compact', extra_args=rescue_args)
    fixture.cli('run', '--output', anchor, '--cpus', 2)
    approved = []
    for target, (i, kind) in expected.items():
        if not manual and kind == 'phase_only':
            continue
        publication = anchor / 'rescued' / target / 'model_repair_candidates.json'
        candidate, = [row for row in json.loads(publication.read_text())
                      if row['owner_overlap'][0]['gene_id'] == f'g{i}']
        assert candidate['kind'] == kind
        if not manual:
            field = 'reciprocal_locus' if case == 'reciprocal' else 'syntenic_locus'
            assert len({hit['donor'] for hit in candidate['support'] if hit.get(field) == 'g8'}) >= 2
        approved.append({'species': target, 'candidate_id': candidate['candidate_id']})
    augmented_cds = anchor / 'augmented' / 'species_cds' / (coding_target + '.rescue.cds.fa')
    assert {gene: seq for gene, _, seq in fasta_records(augmented_cds)}[coding_target + '_g8'] == source_cds

    output = tmp_path / 'refinement'
    refine_args = ()
    if manual:
        approvals = tmp_path / 'repairs.tsv'
        write_tsv(approvals, ['species', 'candidate_id'], approved)
        refine_args = ('--model-repairs', approvals)
    cli('plan', '--rescue-output', anchor, '--output', output, *refine_args)
    cli('run', '--output', output, '--cpus', 2)
    cli('qc', '--output', output)
    effective = output / 'effective'
    inputs = {row['species']: row for row in refinement.verify_inputs(effective / 'inputs.tsv')}
    mapping = {row['gene_id']: row for row in refinement.read_table(effective / 'representative_map.tsv')}
    for target, (i, _) in expected.items():
        for role in ('cds', 'analysis_cds', 'protein'):
            records = {gene: value for gene, _, value in fasta_records(Path(inputs[target][role]))}
            sequence = str(Seq(truth[i]).translate()).removesuffix('*') if role == 'protein' else truth[i]
            assert records[target + f'_g{i}'] == sequence
            assert not any('ggrepair' in gene for gene in records)
        assert mapping[target + f'_g{i}']['status'] == 'repaired_model'
    assert mapping[phase_target + '_g11']['source_transcript_id'] == 'g11.t1'
    def phases(path):
        rows = (line.split('\t') for line in path.read_text().splitlines() if not line.startswith('#'))
        return [row[7] for row in rows if row[2] == 'CDS' and 'Parent=g11.t1' in row[8]]
    assert phases(Path(inputs[phase_target]['gff'])) == ['0', '2']
    assert phases(effective / 'full_annotation' / (phase_target + '.gff3')) == ['0', '2']
    assert phases(effective / 'source_annotation' / (phase_target + '.gff3')) == ['0', '0']
    audit = json.loads((effective / 'model_repairs.json').read_text())
    assert {row['gene_id'] for row in audit} == {coding_target + '_g8', phase_target + '_g11'}
    if case == 'reciprocal':
        rescue = refinement.rescue
        replay = import_module('rescue_admission_replay')
        output_replay = tmp_path / 'replayed'
        replay_plan = replay.build_plan(anchor, SUPPORT.parents[1], output_replay)
        replay.prepare(output_replay, replay_plan, 2)
        with monkeypatch.context() as guard:
            def forbidden(*args, **kwargs):
                raise AssertionError('Admission replay launched a model repair search')
            guard.setattr(rescue, 'search_intervals', forbidden)
            replayed = rescue.rescue(output_replay, replay_plan, coding_target, 2)
        original = anchor / 'rescued' / coding_target
        for member in ('targets.json', rescue.LOCAL_RECORDS, rescue.INTERVAL_ARCHIVE):
            assert (replayed / 'model_repair_search' / member).read_bytes() == (
                original / 'model_repair_search' / member).read_bytes()
        candidates = json.loads((replayed / 'model_repair_candidates.json').read_text())
        match, = (row for row in candidates if row['owner_overlap'][0]['gene_id'] == 'g8')
        assert len({hit['donor'] for hit in match['support'] if hit.get('reciprocal_locus') == 'g8'}) >= 2
    if not manual:
        coding_basis = 'reciprocal_cds_extension' if case == 'reciprocal' else 'syntenic_cds_gain'
        assert {row['basis'] for row in audit} == {coding_basis, 'exact_source_cds'}
        assert {row['adoption'] for row in audit} == {'automatic'}
    review = import_module('plot_gene_model_refinement').collect(output)
    assert review['model_repairs'] == audit
    assert all(path.read_bytes() == content for path, content in before.items())


def native_candidate_table(tmp_path, target, genome, gff_text):
    """Mock only the native predictor boundary, never genomic/donor evidence."""
    gff, receipt, table = (tmp_path / name for name in ('native.gff3', 'native.receipt.json', 'native.tsv'))
    gff.write_text(gff_text)
    value = dict(schema_version=1, kind='genegalleon-annevo-native-prediction', native_returncode=0,
                 gpu_probe={'gpu_verified': True}, model_profile='Magnoliopsida', runtime_digest='a' * 64,
                 model_sha256='b' * 64, source_commit='c' * 40, cache_key='d' * 64,
                 genome_sha256=hashlib.sha256(Path(genome).read_bytes()).hexdigest(),
                 GFF_sha256=hashlib.sha256(gff.read_bytes()).hexdigest(), GFF_bytes=gff.stat().st_size)
    receipt.write_text(json.dumps(value))
    write_tsv(table, ['species', 'gff', 'receipt'], [dict(species=target, gff=str(gff), receipt=str(receipt))])
    return table, gff, receipt


def test_annevo_exon_candidates_pass_real_refinement_finalization_and_source_fences(tmp_path):
    inputs, edges, _rna, sources, truth = truth_fixture(tmp_path)
    target = 'Species_target'
    source = next(s for s in sources if s['species'] == target)
    original = {s[role]: Path(s[role]).read_bytes() for s in sources for role in ('cds', 'gff', 'genome')}
    text = '##gff-version 3\n'
    for kind in ('repair', 'minus', 'intact', 'pseudogene', 'stop', 'gap'):
        path = truth[target, kind]
        text += annotation('chr_' + kind, 'native_' + kind, 'native_' + kind + '.t', path['blocks'], path['strand'])
    table, native_gff, _receipt = native_candidate_table(tmp_path, target, source['genome'], text)
    root = tmp_path / 'native_refinement'
    cli('plan', '--inputs', inputs, '--edges', edges, '--output', root,
        '--annevo-inputs', table, '--padding', 650, '--max-intron', 2000)
    assert refinement.load(root)['request']['annevo']['predictions'][target]['gff'] == str(native_gff)
    cli('run', '--output', root, '--cpus', 2)
    cli('qc', '--output', root)
    predictions = json.loads((root / 'predictions' / target / 'predictions.json').read_text())
    supported = {r['gene_id'].removeprefix(target + '_'): r for r in predictions
                 if r['status'] == 'accepted' and r['candidate'].get('support', {}).get('annevo')}
    assert {'repair', 'minus'} <= supported.keys()
    assert not {'intact', 'pseudogene', 'stop', 'gap'} & supported.keys()
    selected = {i: s for i, _, s in fasta_records(root / 'effective/species_cds' / (target + '.fa'))}
    for kind in ('repair', 'minus', 'intact'):
        assert selected[target + '_' + kind] == truth[target, kind]['cds']
    derived_gff = (root / 'effective/species_gff' / (target + '.gff3')).read_text()
    assert 'prediction_engine=ANNEVO' in derived_gff and 'annevo_native_transcripts=' in derived_gff
    assert all(Path(p).read_bytes() == data for p, data in original.items())
    cli('run', '--output', root, '--cpus', 2)  # Verified prediction/finalization reuse.
    native_gff.write_text(native_gff.read_text() + '\n')
    result = cli('qc', '--output', root, success=False)
    assert 'Frozen refinement input changed' in result.stderr


def test_annevo_missing_gene_candidates_survive_rescue_and_inherited_refinement(tmp_path):
    fixture = import_module('test_rescue_gene_models_runtime')
    data = fixture.hidden_models.__wrapped__(tmp_path)
    directory, names, truth = data
    target = names[0]
    original = {path: path.read_bytes() for role in ('cds', 'gff', 'genome') for path in (directory / role).iterdir()}
    # The complete donor annotation is an explicit mocked native coding path;
    # miniprot/LAST, source DNA and donor-species admission execute normally.
    donor_gff = next((directory / 'gff').glob(names[1] + '*'))
    genome = next((directory / 'genome').glob(target + '*'))
    groups = {}
    for line in donor_gff.read_text().splitlines():
        if line.startswith('#'):
            continue
        f = line.split('\t')
        if len(f) == 9 and f[2] == 'CDS':
            parent = refinement.rescue.attributes(f[8])['Parent']
            groups.setdefault((f[0], f[6], parent), []).append([int(f[3]) - 1, int(f[4]), int(f[7])])
    text = '##gff-version 3\n'
    for i, ((seqid, strand, _parent), blocks) in enumerate(sorted(groups.items())):
        text += annotation(seqid, 'native_' + str(i), 'native_' + str(i) + '.t', blocks, strand)
    table, _native, _receipt = native_candidate_table(tmp_path, target, genome, text)
    root = tmp_path / 'native_rescue'
    fixture.cli('plan', '--cds-dir', directory / 'cds', '--gff-dir', directory / 'gff',
                '--genome-dir', directory / 'genome', '--busco-dir', directory / 'busco',
                '--tree', directory / 'tree.nwk', '--output', root, '--annevo-inputs', table, '--model-storage', 'compact')
    fixture.cli('run', '--output', root, '--cpus', 2)
    models = list(import_module('rescue_model_store').iter_accepted_models(root / 'rescued' / target))
    assert any(m['sequence'] == truth[8] and any(s.get('annevo') for s in m.get('support', [])) for m in models)
    gff = (root / 'augmented/species_gff' / (target + '.rescue.gff3')).read_text()
    assert 'prediction_engine=ANNEVO' in gff
    final = tmp_path / 'native_inherited'
    cli('plan', '--rescue-output', root, '--output', final, '--mode', 'off')
    assert refinement.load(final)['request']['annevo'] == rescue_frozen(root)
    cli('run', '--output', final, '--cpus', 2)
    cli('qc', '--output', final)
    assert truth[8] in {sequence for _, _, sequence in fasta_records(final / 'effective/species_cds' / (target + '.fa'))}
    assert all(path.read_bytes() == raw for path, raw in original.items())
    genome.write_text(genome.read_text() + '>new_contig\nACGT\n')
    with pytest.raises(ValueError, match='ANNEVO candidate genome differs'):
        refinement.plan(tmp_path / 'changed_genome_inherited', rescue_output=root, mode='off')


def rescue_frozen(root):
    return json.loads((root / 'plan.json').read_text())['request']['annevo']


def test_copy_ambiguous_homology_candidate_is_archived_without_active_gff_cds_publication(tmp_path, monkeypatch):
    inputs, edges, _rna, sources, _truth = truth_fixture(tmp_path)
    original = Path(sources[0]['gff']).read_bytes()
    root = tmp_path / 'copy_ambiguity'
    select = refinement.select_from_store
    def unresolved_final_selection(db, links, **kwargs):
        result = select(db, links, **kwargs)
        if Path(db).parent.name == 'catalog_index_final':
            # Mock only the final copy-ambiguity boundary. Genomic prediction,
            # admission, storage and paired GFF/CDS publication execute normally.
            with sqlite3.connect(db) as connection:
                locus = refinement.load_locus(connection, 'Species_target', 'Species_target_repair')
            baseline = refinement._baseline(locus, locus['candidates'])
            for row in result['selections']:
                if row['species'] == 'Species_target' and row['gene_id'] == 'Species_target_repair':
                    if row['candidate_id'] != baseline['candidate_id']:
                        result['metrics']['changed_representatives'] -= 1
                    row.update(candidate_id=baseline['candidate_id'], source_transcript_id=baseline['source_transcript_id'],
                               status='ambiguous_correspondence', reason='unresolved_paralog_copy')
            for row in result['scores']:
                if row['species'] == 'Species_target' and row['gene_id'] == 'Species_target_repair':
                    row['selected'] = row['candidate_id'] == baseline['candidate_id']
        return result
    monkeypatch.setattr(refinement, 'select_from_store', unresolved_final_selection)
    value = refinement.plan(root, inputs=inputs, edges=edges, padding=650, max_intron=2000)
    refinement.finalize(root, value)
    cli('qc', '--output', root)
    predictions = json.loads((root / 'predictions/Species_target/predictions.json').read_text())
    candidates = {r['candidate']['candidate_id'] for r in predictions
                  if r['gene_id'] == 'Species_target_repair' and r['status'] == 'accepted'}
    assert candidates  # Candidate admission and publication remain separate.
    full = (root / 'effective/full_annotation/Species_target.gff3').read_text()
    paired = {name for name, _, _ in fasta_records(root / 'effective/all_candidates/Species_target.fa')}
    audit = json.loads((root / 'effective/unpublished_candidates.json').read_text())
    assert candidates <= {r['candidate_id'] for r in audit if r['reason'] == 'unresolved_paralog_copy'}
    assert not candidates & paired
    assert all('ID=' + name + ';' not in full for name in candidates)
    assert (root / 'effective/source_annotation/Species_target.gff3').read_bytes() == original
    recatalog = refinement.build_catalog('Species_target', root / 'effective/all_candidates/Species_target.fa',
                                         root / 'effective/full_annotation/Species_target.gff3',
                                         root / 'effective/species_genome/Species_target.fa')
    assert not candidates & {c['source_transcript_id'] for g in recatalog['loci'] for c in g['candidates']}


def coding(seed, residues=280):
    rng = random.Random(seed)
    codons = {}
    for codon, amino_acid in sorted(CodonTable.unambiguous_dna_by_id[1].forward_table.items()):
        codons.setdefault(amino_acid, codon)
    amino_acids = sorted(codons)
    return "ATG" + "".join(codons[rng.choice(amino_acids)] for _ in range(residues - 1)) + "TAA"


def spliced_gene(sequence, strand="+", flank=1100):
    pieces = [sequence[:300], sequence[300:540], sequence[540:]]
    intron = "GT" + "AC" * 38 + "AG"
    genomic = "AC" * (flank // 2)
    blocks = []
    for index, piece in enumerate(pieces):
        start = len(genomic)
        genomic += piece
        blocks.append([start, len(genomic), 0])
        if index + 1 < len(pieces):
            genomic += intron
    genomic += "CA" * (flank // 2)
    if strand == "-":
        size = len(genomic)
        blocks = [[size - end, size - start, phase] for start, end, phase in blocks]
        genomic = str(Seq(genomic).reverse_complement())
    return genomic, blocks


def annotation(contig, gene, transcript, blocks, strand="+", attributes=""):
    start, end = min(block[0] for block in blocks), max(block[1] for block in blocks)
    lines = [f"{contig}\tsynthetic\tgene\t{start + 1}\t{end}\t.\t{strand}\t.\tID={gene}{attributes}",
             f"{contig}\tsynthetic\tmRNA\t{start + 1}\t{end}\t.\t{strand}\t.\tID={transcript};Parent={gene}"]
    for number, (a, b, phase) in enumerate(blocks):
        lines += [f"{contig}\tsynthetic\tCDS\t{a + 1}\t{b}\t.\t{strand}\t{phase}\tID={transcript}.c{number};Parent={transcript}",
                  f"{contig}\tsynthetic\texon\t{a + 1}\t{b}\t.\t{strand}\t.\tID={transcript}.e{number};Parent={transcript}"]
    return "\n".join(lines) + "\n"


def truth_fixture(tmp_path, *, rna=False):
    """Target annotations lose exons/paths; DNA truth and donors stay complete."""
    names = ["Species_target", "Species_donor1", "Species_donor2"]
    kinds = ["repair", "minus", "alternate", "intact", "pseudogene", "stop", "gap"]
    truth, sources = {}, []
    for name in names:
        genome, gff, fasta = [], ["##gff-version 3\n"], []
        for number, kind in enumerate(kinds):
            strand = "-" if kind == "minus" else "+"
            sequence = coding(104 + number)
            if name == names[0] and kind == "stop":
                sequence = sequence[:390] + "TGA" + sequence[393:]
            if name == names[0] and kind == "gap":
                sequence = sequence[:390] + "NNN" + sequence[393:]
            dna, full_blocks = spliced_gene(sequence, strand)
            contig = "chr_" + kind
            genome.append(f">{contig}\n{dna}\n")
            target = name == names[0]
            blocks = full_blocks[:2] if target and kind in {"repair", "minus"} else full_blocks
            # Provider parent IDs and formatter locus tokens are deliberately
            # different; original and predicted paths must share both correctly.
            source_gene = "provider_" + kind
            attributes = ";gene_id=" + kind + (";pseudo=true" if target and kind == "pseudogene" else "")
            gff.append(annotation(contig, source_gene, kind + ".long", blocks, strand, attributes))
            supplied = "".join(sequence[a:b] for a, b in ((0, 300), (300, 540))) if target and kind in {"repair", "minus"} else sequence
            fasta.append(f">{kind}.long\n{supplied}\n")
            if kind == "alternate" and not target:
                short = [full_blocks[0], full_blocks[2]]
                # An exon-skipped coding path exists in both donor annotations.
                gff.append(annotation(contig, source_gene, kind + ".short", short, strand).split("\n", 1)[1])
                fasta.append(f">{kind}.short\n{sequence[:300] + sequence[540:]}\n")
            truth[name, kind] = dict(cds=sequence, blocks=full_blocks, strand=strand,
                                    shorter_cds=sequence[:300] + sequence[540:])
        paths = {role: tmp_path / (name + suffix) for role, suffix in
                 (("cds", ".cds.fa"), ("gff", ".gff3"), ("genome", ".genome.fa"))}
        for role, content in (("cds", fasta), ("gff", gff), ("genome", genome)):
            paths[role].write_text("".join(content))
        sources.append(dict(species=name, genetic_code=1, **{role: str(path) for role, path in paths.items()}))
    inputs, edges = tmp_path / "inputs.tsv", tmp_path / "edges.tsv"
    write_tsv(inputs, list(sources[0]), sources)
    links = [dict(species_a=names[i], gene_a=names[i] + "_" + kind,
                  species_b=names[j], gene_b=names[j] + "_" + kind, weight=1)
             for kind in kinds for i, j in ((0, 1), (0, 2), (1, 2))]
    write_tsv(edges, list(links[0]), links)
    rna_path = None
    if rna:
        full = truth[names[0], "alternate"]["blocks"]
        rna_path = tmp_path / "rna.tsv"
        evidence = dict(species=names[0], seqid="chr_alternate", strand="+",
                        cds_blocks=json.dumps([[a, b] for a, b, _ in (full[0], full[2])]),
                        transcript_id="target_short_rna", count=12)
        write_tsv(rna_path, list(evidence), [evidence])
    return inputs, edges, rna_path, sources, truth


def cli(*args, success=True):
    result = subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], text=True, capture_output=True)
    if success:
        assert result.returncode == 0, result.stdout + result.stderr
    else:
        assert result.returncode != 0, result.stdout + result.stderr
    return result


def test_transcriptome_assembly_refinement_uses_protein_correspondence(tmp_path):
    target, donors = 'Species_assembly', ('Species_donor1', 'Species_donor2')
    body = 'ATG' + 'GCTGAAGGTTTCATCCAAGAC' * 11
    matching, longest = body + 'TAA', body + 'TGG' * 20 + 'TAA'
    assembly = tmp_path / 'assembly.fa.gz'
    with gzip.open(assembly, 'wt') as stream:
        for gene in (1, 2):
            stream.write(f'>TRINITY_DN1_c0_g{gene}_i1\n{longest}\n>TRINITY_DN1_c0_g{gene}_i2\n{matching}\n')
    prepared = tmp_path / 'prepared'
    cli('prepare-transcriptome', '--assembly', assembly, '--species', target, '--assembler', 'trinity', '--output', prepared)
    assert (prepared / (target + '_initial_longestCDS.fa')).read_text().count(longest) == 2
    assert (prepared / 'all_isoforms.cds.fa').read_text().count('>') == 4
    genomic_sources = []
    for species in donors:
        cds, gff, genome = (tmp_path / (species + ext) for ext in ('.cds.fa', '.gff3', '.genome.fa'))
        cds.write_text(f'>one\n{matching}\n>two\n{matching}\n')
        genome.write_text(f'>chr1\n{matching}\n>chr2\n{matching}\n')
        gff.write_text('##gff-version 3\n' + annotation('chr1', 'one', 't1', [(0, len(matching), 0)])
                       + annotation('chr2', 'two', 't2', [(0, len(matching), 0)]))
        genomic_sources.append(dict(species=species, cds=cds, gff=gff, genome=genome))
    genomic = tmp_path / 'genomic.tsv'
    write_tsv(genomic, list(genomic_sources[0]), genomic_sources)
    edges = tmp_path / 'edges.tsv'
    gene_one, gene_two = (f'{target}_DN1-c0-g{gene}' for gene in (1, 2))
    edge_rows = [dict(species_a=target, gene_a=gene_one, species_b=species, gene_b=species + '_one', ambiguous=0)
                 for species in donors]
    edge_rows.append(dict(species_a=target, gene_a=gene_two, species_b=donors[0], gene_b=donors[0] + '_two', ambiguous=1))
    write_tsv(edges, list(edge_rows[0]), edge_rows)
    root = tmp_path / 'refinement'
    cli('plan', '--inputs', genomic, '--transcriptome-inputs', prepared / 'inputs.tsv', '--edges', edges, '--output', root)
    cli('run', '--output', root)
    cli('qc', '--output', root)
    effective = root / 'effective'
    assert {row['species'] for row in refinement.verify_inputs(effective / 'inputs.tsv')} == set(donors)
    sequence = effective / 'sequence'
    assert {row['species'] for row in refinement.read_table(sequence / 'inputs.tsv')} == {target, *donors}
    choices = {row['gene_id']: row for row in refinement.read_table(sequence / 'representative_map.tsv') if row['species'] == target}
    assert (choices[gene_one]['source_transcript_id'], choices[gene_one]['baseline_transcript_id'], choices[gene_one]['status']) == (
        'TRINITY_DN1_c0_g1_i2', 'TRINITY_DN1_c0_g1_i1', 'conserved')
    assert (choices[gene_two]['source_transcript_id'], choices[gene_two]['reason']) == (
        'TRINITY_DN1_c0_g2_i1', 'unresolved_paralog_copy')
    assert (sequence / 'species_cds' / (target + '.fa')).read_text() == f'>{gene_one}\n{matching}\n>{gene_two}\n{longest}\n'
    assert not (effective / 'species_gff' / (target + '.gff3')).exists()
    assert json.loads((root / 'predictions' / target / 'predictions.json').read_text()) == []
    assert json.loads((root / 'predictions' / donors[0] / 'predictions.json').read_text()) == []
    with (prepared / 'transcript_gene_map.tsv').open('a') as stream:
        stream.write('\n')
    assert 'Frozen refinement input changed' in cli('qc', '--output', root, success=False).stderr


@pytest.mark.parametrize('mixed', [False, True])
def test_discovered_assemblies_infer_edges_and_abstain_on_duplicate_copies(tmp_path, mixed):
    fixture = import_module('test_rescue_gene_models_runtime')
    data = fixture.hidden_models.__wrapped__(tmp_path)
    _, _, sequences = data
    directory = tmp_path / 'transcriptome_assembly'
    assemblies = directory / 'assembled_transcripts_with_isoforms'
    assemblies.mkdir(parents=True)
    matching = sequences[0]
    baseline = ''.join(sequences[2][i:i + 3] if 60 <= i < 420 and i % 9 == 0 else matching[i:i + 3]
                       for i in range(0, len(matching), 3))
    target = 'Assembly_target'
    donors = ('Assembly_donor1', 'Assembly_donor2')
    for species in (target, *donors):
        with gzip.open(assemblies / (species + '_isoform.fa.gz'), 'wt') as handle:
            if species == target:
                handle.write(f'>TRINITY_DN1_c0_g1_i1\n{baseline}\n>TRINITY_DN1_c0_g1_i2\n{matching}\n')
                handle.write(f'>TRINITY_DN1_c0_g2_i1\n{sequences[1]}\n')
            else:
                handle.write(f'>TRINITY_DN1_c0_g1_i1\n{matching}\n')
                handle.write(f'>TRINITY_DN1_c0_g2_i1\n{sequences[1]}\n>TRINITY_DN1_c0_g3_i1\n{sequences[1]}\n')
    genomic_args = ()
    if mixed:
        species = donors[1]
        cds, gff, genome = (tmp_path / (species + ext) for ext in ('.cds.fa', '.gff3', '.genome.fa'))
        proteins = [matching, sequences[1], sequences[1]]
        cds.write_text(''.join(f'>DN1-c0-g{i}\n{seq}\n' for i, seq in enumerate(proteins, 1)))
        genome.write_text(''.join(f'>chr{i}\n{seq}\n' for i, seq in enumerate(proteins, 1)))
        gff.write_text('##gff-version 3\n' + ''.join(
            annotation(f'chr{i}', f'DN1-c0-g{i}', f't{i}', [[0, len(seq), 0]])
            for i, seq in enumerate(proteins, 1)))
        genomic = tmp_path / 'genomic.tsv'
        write_tsv(genomic, ['species', 'cds', 'gff', 'genome'],
                  [{'species': species, 'cds': cds, 'gff': gff, 'genome': genome}])
        genomic_args = ('--inputs', genomic)
    output = tmp_path / 'automatic'
    cli('plan', '--transcriptome-assembly-dir', directory, '--output', output, *genomic_args)
    cli('run', '--output', output, '--cpus', 2)
    cli('qc', '--output', output)
    edges = json.loads((output / 'correspondence' / 'edges.json').read_text())
    assert sum(target in (row['species_a'], row['species_b']) for row in edges) == 2
    assert all(row['gene_a'].endswith('_DN1-c0-g1') and row['gene_b'].endswith('_DN1-c0-g1') for row in edges)
    choices = {row['gene_id']: row for row in refinement.read_table(output / 'effective/sequence/representative_map.tsv')}
    assert choices[target + '_DN1-c0-g1']['source_transcript_id'] == 'TRINITY_DN1_c0_g1_i2'
    assert choices[target + '_DN1-c0-g2']['source_transcript_id'] == 'TRINITY_DN1_c0_g2_i1'
    if mixed:
        rows = refinement.verify_inputs(output / 'effective/inputs.tsv')
        assert {row['species'] for row in rows} == {donors[1]}


@pytest.mark.parametrize("rna", [False, True, "bam"])
def test_real_miniprot_repairs_heldout_exons_adds_paths_and_preserves_biological_negatives(tmp_path, rna):
    inputs, edges, evidence, sources, truth = truth_fixture(tmp_path, rna=rna is True)
    original = {row[role]: Path(row[role]).read_bytes() for row in sources for role in ("cds", "gff", "genome")}
    root = tmp_path / "run"
    params = ["plan", "--inputs", inputs, "--edges", edges, "--output", root,
              "--padding", 650, "--max-intron", 2000]
    if evidence:
        params += ["--rna", evidence]
    if rna == "bam":
        pysam = import_module('pysam')
        bam_dir = tmp_path / 'rna_bams'
        bam_dir.mkdir()
        genome = {name: seq for name, _, seq in fasta_records(Path(sources[0]['genome']))}
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
                  'SQ': [{'SN': key, 'LN': len(seq)} for key, seq in genome.items()]}
        species = sources[0]['species']
        blocks = truth[species, 'alternate']['blocks']
        a, b, _ = blocks[0]
        c, d, _ = blocks[2]
        with pysam.AlignmentFile(str(bam_dir / (species + '.bam')), 'wb', header=header) as bam:
            read = pysam.AlignedSegment(bam.header)
            read.query_name = 'whole_coding_chain'
            read.query_sequence = genome['chr_alternate'][a:b] + genome['chr_alternate'][c:d]
            read.reference_id = bam.get_tid('chr_alternate')
            read.reference_start = a
            read.mapping_quality = 60
            read.query_name = 'deletion_instead_of_splice'
            read.cigartuples = [(0, b - a), (2, c - b), (0, d - c)]
            bam.write(read)
            read.query_name = 'whole_coding_chain'
            read.cigartuples = [(0, b - a), (3, c - b), (0, d - c)]
            read.set_tag('XS', 50)
            bam.write(read)
        params += ['--rna-bam-dir', bam_dir]
    cli(*params)
    for command in ("catalog", "correspondence", "select", "predict", "finalize", "qc"):
        cli(command, "--output", root, "--cpus", 2)
    target = "Species_target"
    predictions = json.loads((root / "predictions" / target / "predictions.json").read_text())
    accepted = {row["gene_id"].removeprefix(target + "_"): row for row in predictions if row["status"] == "accepted"}
    for kind in ("repair", "minus"):
        assert kind in accepted, [(row["gene_id"], row["status"], row["problems"]) for row in predictions]
        model = accepted[kind]
        assert model["candidate"]["cds"] == truth[target, kind]["cds"]
        assert model["candidate"]["blocks"] == truth[target, kind]["blocks"]
        assert model["candidate"]["strand"] == truth[target, kind]["strand"]
        assert model["change_type"] == "model_revision"
        assert model["donors"] == ["Species_donor1", "Species_donor2"]
    alternative = next(row for row in predictions if row['gene_id'] == target + '_alternate')
    assert ("alternate" in accepted) is bool(rna)
    if not rna:
        assert alternative['status'] == 'proposal'
        assert 'unconfirmed_source_structure_change' in alternative['problems']
    assert alternative["candidate"]["cds"] == truth[target, "alternate"]["shorter_cds"]
    assert alternative["change_type"] == "isoform_addition"
    assert alternative["evidence_class"] == ("rna_path_supported" if rna else "homology_only_predicted")
    assert alternative["candidate"]["quality"]["representative_eligible"] is bool(rna)
    if rna == "bam":
        assert alternative['rna_paths'] == [sources[0]['species'] + '.bam:whole_coding_chain']
    assert not ({"intact", "pseudogene", "stop", "gap"} & set(accepted))
    effective = root / "effective"
    summary = json.loads((root / "catalog_index_final/summary.json").read_text())
    final_db = root / "catalog_index_final/loci.sqlite3"
    candidates = sum(len(locus['candidates']) for locus in refinement.iter_loci(final_db))
    assert Path(summary['database']) == final_db
    assert summary['candidate_count'] == candidates
    assert summary['accepted_prediction_count'] >= len(accepted)
    selected = {identifier: sequence for identifier, _header, sequence in
                fasta_records(effective / "species_cds" / (target + ".fa"))}
    proteins = {identifier: sequence for identifier, _header, sequence in
                fasta_records(effective / "species_protein" / (target + ".fa"))}
    assert not ({target + "_" + kind for kind in ("pseudogene", "stop", "gap")} & set(proteins))
    for kind in ("repair", "minus", "intact", "alternate"):
        assert target + "_" + kind in proteins
        assert proteins[target + "_" + kind] == str(Seq(selected[target + "_" + kind]).translate()).removesuffix("*")
    assert not any("*" in sequence for sequence in proteins.values())
    for kind in ("repair", "minus", "intact", "pseudogene", "stop", "gap"):
        assert selected[target + "_" + kind] == truth[target, kind]["cds"]
    if not rna:
        assert selected[target + "_alternate"] == truth[target, "alternate"]["cds"]
    full = (effective / "full_annotation" / (target + ".gff3")).read_text()
    assert all("ID=" + kind + ".long;" in full for kind in ("repair", "minus", "alternate"))
    assert ("ID=" + alternative["candidate"]["source_transcript_id"] + ";" in full) is bool(rna)
    reimported = refinement.build_catalog(target, effective / "all_candidates" / (target + ".fa"),
                                          effective / "full_annotation" / (target + ".gff3"),
                                          effective / "species_genome" / (target + ".fa"))
    rebound = {candidate["source_transcript_id"]: (locus, candidate) for locus in reimported["loci"]
               for candidate in locus["candidates"]}
    for kind in ("repair", "minus", "alternate"):
        transcripts = [kind + '.long']
        if kind in accepted:
            transcripts.append(accepted[kind]['candidate']['source_transcript_id'])
        for transcript in transcripts:
            locus, candidate = rebound[transcript]
            assert locus["gene_id"] == target + "_" + kind
            assert candidate["gene_token"] == kind
            assert candidate["source_gene_id"] == "provider_" + kind
    assert (effective / "source_annotation" / (target + ".gff3")).read_bytes() == original[sources[0]["gff"]]
    # Downstream reader reports the very transcript/coordinates behind selected CDS.
    choices = {row["gene_id"]: row for row in refinement.read_table(effective / "representative_map.tsv") if row["species"] == target}
    traits = gff2genestat.process_single_gff(target + ".gff3", str(effective / "species_gff"),
                                           list(selected), "CDS", "longest", GFF_COLS, OUT_COLS,
                                           representative_map=effective / "representative_map.tsv")
    for row in traits.itertuples(index=False):
        assert row.gff_transcript_id == choices[row.gene_id]["source_transcript_id"]
        assert row.feature_size == len(selected[row.gene_id])
    assert all(Path(path).read_bytes() == content for path, content in original.items())
    assert not list(tmp_path.glob("*.fai"))
    receipts = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in root.rglob("receipt.json")}
    cli("run", "--output", root, "--cpus", 2)
    assert receipts == {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in root.rglob("receipt.json")}
    cli("qc", "--output", root)


def test_real_miniprot_wrong_neighbor_hits_are_proposals_not_repairs(tmp_path):
    names = ["Species_target", "Species_donor1", "Species_donor2"]
    sequence = coding(800, residues=300)
    donor_dna, donor_blocks = spliced_gene(sequence)
    rows = []
    for name in names:
        cds, gff, genome = [tmp_path / (name + suffix) for suffix in (".cds.fa", ".gff3", ".genome.fa")]
        if name == names[0]:
            # The owned locus is physically truncated. A nearby annotated duplicate
            # carries the matching full protein and lies inside the search window.
            partial = sequence[:540]
            clone_start = 1100
            dna = "AC" * 100 + partial + "AC" * ((clone_start - 200 - len(partial)) // 2) + sequence + "AC" * 100
            genome.write_text(">chr1\n" + dna + "\n")
            cds.write_text(">wrong.t\n" + partial + "\n>neighbor.t\n" + sequence + "\n")
            gff.write_text("##gff-version 3\n" + annotation("chr1", "wrong", "wrong.t", [[200, 740, 0]])
                           + annotation("chr1", "neighbor", "neighbor.t", [[clone_start, clone_start + len(sequence), 0]]))
        else:
            genome.write_text(">chr1\n" + donor_dna + "\n")
            cds.write_text(">wrong.t\n" + sequence + "\n")
            gff.write_text("##gff-version 3\n" + annotation("chr1", "wrong", "wrong.t", donor_blocks))
        rows.append(dict(species=name, cds=str(cds), gff=str(gff), genome=str(genome), genetic_code=1))
    inputs, edges = tmp_path / "inputs.tsv", tmp_path / "edges.tsv"
    write_tsv(inputs, list(rows[0]), rows)
    links = [dict(species_a=names[0], gene_a=names[0] + "_wrong", species_b=donor, gene_b=donor + "_wrong")
             for donor in names[1:]]
    write_tsv(edges, list(links[0]), links)
    root = tmp_path / "run"
    cli("plan", "--inputs", inputs, "--edges", edges, "--output", root, "--padding", 1600)
    cli("predict", "--output", root, "--task-index", 3 if sorted(names).index(names[0]) == 2 else sorted(names).index(names[0]) + 1)
    predictions = json.loads((root / "predictions" / names[0] / "predictions.json").read_text())
    assert predictions
    assert not any(row["status"] == "accepted" for row in predictions)
    assert any("overlap_other_locus" in row["problems"] for row in predictions)
    assert any("no_overlap_owned_locus" in row["problems"] for row in predictions)


@pytest.mark.parametrize('with_assembly', [False, True])
def test_real_synteny_reuse_recovers_correspondence_for_existing_anchor_excluded_model(tmp_path, with_assembly):
    """A masked final CDS base removes an anchor while its two flanks survive."""
    rescue = refinement.rescue
    names = ["Plant_species0", "Plant_species1", "Plant_species2"]
    truth = [coding(1400 + number, residues=240) for number in range(7)]
    for role in ("cds", "gff", "genome", "busco"):
        (tmp_path / role).mkdir()
    for name in names:
        genome, annotations, records = "", ["##gff-version 3\n"], []
        for number, sequence in enumerate(truth):
            start = len(genome)
            genome += sequence + "AC" * 60
            end = start + len(sequence) - (1 if name == names[0] and number == 3 else 0)
            annotations.append(annotation("chr1", "g" + str(number), "g" + str(number) + ".t", [[start, end, 0]]))
            records.append(f">{name}_g{number}\n{sequence[:end - start]}\n")
        (tmp_path / "cds" / (name + ".cds.fa")).write_text("".join(records))
        (tmp_path / "gff" / (name + ".gff3")).write_text("".join(annotations))
        (tmp_path / "genome" / (name + ".genome.fa")).write_text(">chr1\n" + genome + "\n")
        (tmp_path / "busco" / (name + ".busco.short.txt")).write_text(
            "# BUSCO version is: 6.0.0\n# The lineage dataset is: embryophyta_odb12\n"
            "# BUSCO was run in mode: transcriptome\nC:95.0%[S:95.0%,D:0.0%],F:0.0%,M:5.0%,n:100\n")
    tree = tmp_path / "tree.nwk"
    tree.write_text("(Plant_species0:1,(Plant_species1:1,Plant_species2:1):1);\n")
    anchor_root, root = tmp_path / "anchors", tmp_path / "refinement"
    command = [sys.executable, rescue.__file__, "plan", "--cds-dir", str(tmp_path / "cds"),
               "--gff-dir", str(tmp_path / "gff"), "--genome-dir", str(tmp_path / "genome"),
               "--busco-dir", str(tmp_path / "busco"), "--tree", str(tree), "--output", str(anchor_root),
               "--common-references", "2", "--nearest-references", "1", "--min-anchors", "3"]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    anchors = rescue.load(anchor_root)
    jobs = [job for job in anchors["synteny_jobs"] if job["a"] != job["b"]]
    assert len(jobs) <= 3
    for job in jobs:
        rescue.synteny(anchor_root, anchors, job["index"], 2, tmp_path / "comparisons")
    receipts = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in (anchor_root / "synteny").rglob("receipt.json")}
    mapping = rescue.table(anchor_root / "prepared" / names[0] / "genes.id_map.tsv")
    excluded = next(row for row in mapping if row["original_id"] == names[0] + "_g3")
    assert excluded["status"] == "translation_excluded"
    assembly_args = []
    if with_assembly:
        species = 'Plant_transcriptome'
        (tmp_path / 'assembly.fa').write_text(f'>transcript\n{truth[3]}\n')
        (tmp_path / 'map.tsv').write_text(f'transcript_id\tgene_id\ntranscript\t{species}_g3\n')
        (tmp_path / 'assembly.tsv').write_text(f'species\tcds\ttranscript_gene_map\n{species}\tassembly.fa\tmap.tsv\n')
        edges = tmp_path / 'edges.tsv'
        edges.write_text('species_a\tgene_a\tspecies_b\tgene_b\n'
                         + ''.join(f'{species}\t{species}_g3\t{name}\t{name}_g3\n' for name in names))
        assembly_args = ['--transcriptome-inputs', tmp_path / 'assembly.tsv', '--edges', edges]
    cli("plan", "--rescue-output", anchor_root, "--output", root, "--padding", 150, *assembly_args)
    cli("correspondence", "--output", root, "--cpus", 2, "--comparison-cache", tmp_path / "comparisons")
    actual = json.loads((root / "correspondence" / "edges.json").read_text())
    restored = [edge for edge in actual if edge["gene_a"].endswith("_g3") and edge["gene_b"].endswith("_g3")
                and names[0] in {edge["species_a"], edge["species_b"]} and 'comparison_basis' not in edge]
    assert len(restored) == 2, actual
    assert all(edge["evidence"]["kind"] == "two_flanking_anchors" and not edge["ambiguous"] for edge in restored)
    if with_assembly:
        assert sum(edge.get('comparison_basis') == 'protein_alignment' for edge in actual) == len(names)
    assert receipts == {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in (anchor_root / "synteny").rglob("receipt.json")}
    cli("predict", "--output", root, "--task-index", 1, "--cpus", 2)
    predictions = json.loads((root / "predictions" / names[0] / "predictions.json").read_text())
    repaired = [row for row in predictions if row["gene_id"] == names[0] + "_g3" and row["status"] == "accepted"]
    assert len(repaired) == 1, predictions
    assert repaired[0]["candidate"]["cds"] == truth[3]
    assert repaired[0]["donors"] == names[1:]


def test_real_miniprot_repairs_source_cds_mismatch_without_changing_coordinates(tmp_path):
    """Independent donor predictions repair DNA-consistent original coding blocks."""
    inputs, edges, _, sources, truth = truth_fixture(tmp_path)
    target = sources[0]["species"]
    expected = truth[target, "intact"]["cds"]
    supplied = expected[:-3] + "NNN" + "AAA" * 18 + "NNN"
    source_path = Path(sources[0]["cds"])
    records = list(fasta_records(source_path))
    source_path.write_text("".join(
        ">" + identifier + "\n" + (supplied if identifier == "intact.long" else sequence) + "\n"
        for identifier, _header, sequence in records))
    original_source = source_path.read_bytes()
    original_annotation = Path(sources[0]["gff"]).read_bytes()
    root = tmp_path / "refinement"
    cli("plan", "--inputs", inputs, "--edges", edges, "--output", root,
        "--mode", "conservative", "--policy", "conserved",
        "--isoform-adoption", "conservation_supported", "--min-support", 2,
        "--minimum-coverage", 0.99, "--minimum-identity", 0.99)
    cli("run", "--output", root, "--cpus", 1)
    cli("qc", "--output", root)
    cli("verify-inputs", "--inputs", root / "effective" / "inputs.tsv")
    accepted = [row for row in json.loads(
        (root / "predictions" / target / "predictions.json").read_text())
        if row.get("status") == "accepted" and row.get("gene_id") == target + "_intact"]
    assert len(accepted) == 1
    repair = accepted[0]
    assert repair["change_type"] == "model_revision"
    assert repair["candidate"]["blocks"] == truth[target, "intact"]["blocks"]
    assert repair["candidate"]["cds"] == expected
    assert repair["same_coordinate_source_repair"]["independent_external_species"] == [
        source["species"] for source in sources[1:]]
    assert not repair["candidate"]["quality"]["rna_supported"]
    selection = json.loads((root / "selection_final" / "selection.json").read_text())
    chosen = next(row for row in selection["selections"]
                  if row["species"] == target and row["gene_id"] == target + "_intact")
    assert chosen["status"] == "sequence_repaired"
    assert chosen["source_transcript_id"] == "intact.long"
    assert chosen["candidate_id"] == repair["candidate"]["candidate_id"]
    effective = root / "effective"
    assert (effective / "source_cds" / (target + ".fa")).read_bytes() == original_source
    assert (effective / "source_annotation" / (target + ".gff3")).read_bytes() == original_annotation
    proteins = {identifier: sequence for identifier, _header, sequence in
                fasta_records(effective / "species_protein" / (target + ".fa"))}
    assert proteins[target + "_intact"] == str(Seq(expected).translate()).removesuffix("*")
    selected = (effective / "species_gff" / (target + ".gff3")).read_text()
    matching = [line for line in selected.splitlines() if "\tmRNA\t" in line and "ID=intact.long;" in line]
    assert len(matching) == 1
