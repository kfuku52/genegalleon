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

import pysam
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


@pytest.mark.parametrize("model_storage", ["legacy", "compact"])
def test_refinement_preserves_missing_models_from_completed_real_rescue(tmp_path, monkeypatch, model_storage):
    from workflow.support.rescue_model_store import iter_accepted_models

    fixture = import_module('test_rescue_gene_models_runtime')
    data = fixture.hidden_models.__wrapped__(tmp_path)
    anchor_root, names, truth = fixture.make_plan(data, model_storage=model_storage)
    before = {path: path.read_bytes() for role in ('cds', 'gff', 'genome') for path in (tmp_path / role).iterdir()}
    fixture.cli('run', '--output', anchor_root, '--cpus', 2)
    worker = anchor_root / 'rescued' / names[0]
    accepted = list(iter_accepted_models(worker))
    assert len(accepted) == 1 and accepted[0]['sequence'] == truth[8]
    root = tmp_path / 'refinement'
    producer_plan = (anchor_root / 'plan.json').read_bytes()
    identities = refinement.rescue.identities()
    monkeypatch.setattr(refinement.rescue, 'identities', lambda: {
        **identities, 'genome_index_implementation': 'downstream-catalogue-update'})
    with pytest.raises(ValueError, match='Rescue schema/tools changed'):
        refinement.rescue.load(anchor_root)
    # Published producers remain readable after a downstream-only helper update.
    refinement.plan(root, rescue_output=anchor_root, mode='off', cpus=4)
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
    worker_file = worker / ('model_store/accepted.sqlite' if model_storage == 'compact' else 'models.json')
    changed = tmp_path / 'changed-worker'
    original_identities = refinement.dependency_identities

    def replace_verified_worker():
        replacement = worker_file.with_suffix('.replacement')
        replacement.write_bytes(worker_file.read_bytes())
        replacement.replace(worker_file)
        return original_identities()

    with monkeypatch.context() as patch:
        patch.setattr(refinement, 'dependency_identities', replace_verified_worker)
        with pytest.raises(ValueError, match='Stage dependency content changed'):
            refinement.plan(changed, rescue_output=anchor_root, mode='off', cpus=4)
    assert not (changed / 'plan.json').exists()
    augmented = anchor_root / 'augmented/species_cds' / (names[0] + '.rescue.cds.fa')
    augmented.write_text(augmented.read_text() + '\n')
    with pytest.raises(ValueError, match='Frozen refinement input changed'):
        refinement.load(root)
    with pytest.raises(ValueError, match='Augmented rescue inputs incomplete or corrupt'):
        refinement.plan(tmp_path / 'new-refinement', rescue_output=anchor_root, mode='off')


@pytest.mark.parametrize(('fragment_count', 'route'), [
    (2, 'partial'), (3, 'partial'), (3, 'manual'),
    (3, 'intact_without_prediction'), (3, 'intact_with_oriongeno'), (3, 'intact_with_other_path'),
])
def test_split_model_merge_from_real_rescue(tmp_path, fragment_count, route):
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
        boundaries = [0, 159, 321 if route == 'partial' else 318, len(sequence)]
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
    approval_args = ('--split-model-merges', approvals) if route == 'manual' else ()
    native_args = ()
    if route in {'intact_with_oriongeno', 'intact_with_other_path'}:
        blocks = [block[:] for block in review['cds']]
        if route == 'intact_with_other_path':
            blocks[0][1] += 3
        text = '##gff-version 3\n' + annotation('chr1', 'native_merge', 'native_merge.t1', blocks, '-')
        table, _, _ = native_candidate_table(tmp_path, target, genome, text, engine='oriongeno')
        native_args = ('--oriongeno-inputs', table)
    cli('plan', '--rescue-output', anchor, '--output', output, *approval_args, *native_args)
    cli('run', '--output', output, '--cpus', 2)
    cli('qc', '--output', output)
    effective = output / 'effective'
    inputs = next(row for row in refinement.verify_inputs(effective / 'inputs.tsv') if row['species'] == target)
    if route in {'intact_without_prediction', 'intact_with_other_path'}:
        records = {gene for gene, _, _ in fasta_records(Path(inputs['cds']))}
        assert original_ids <= records and identifier not in records
        assert not (effective / 'split_model_merges.json').exists()
        return
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
    if route != 'manual':
        assert audit[0]['adoption'] == 'automatic'
    if route == 'intact_with_oriongeno':
        assert audit[0]['basis'] == 'oriongeno_and_complementary_donors'
        assert audit[0]['oriongeno_transcripts'] == ['native_merge.t1']
        assert len({hit['donor'] for hit in audit[0]['support']}) >= 2
        for hit in audit[0]['support']:
            spans = hit['fragment_donor_intervals']
            assert [span['gene_id'] for span in spans] == source_genes
            assert all(a['end'] <= b['start'] for a, b in zip(spans, spans[1:], strict=False))

    review_data = import_module('plot_gene_model_refinement').collect(output)
    assert review_data['split_model_merges'] == audit
    assert all(path.read_bytes() == content for path, content in before.items())
    assert (worker / 'receipt.json').read_bytes() == producer_receipt
    if route == 'manual':
        approvals.write_text(approvals.read_text() + '\n')
        with pytest.raises(ValueError, match='Frozen refinement input changed'):
            refinement.load(output)


@pytest.mark.parametrize(('strand', 'reference_case'), [
    ('+', 'conserved'), ('-', 'conserved'), ('+', 'divergent_bridge'), ('+', 'unmapped_reference'),
    ('+', 'competing_path'), ('+', 'generated_reference'),
])
def test_single_exon_split_merge_revalidates_same_species_copy(tmp_path, strand, reference_case):
    name = 'Plant_example'
    sequence = coding(107)
    bridge_start, bridge_end = 330, 378
    sequence = sequence[:bridge_end] + 'ATG' + sequence[bridge_end + 3:]
    reference = sequence
    if reference_case == 'divergent_bridge':
        alternative = coding(913)
        reference = sequence[:bridge_start] + alternative[bridge_start:bridge_end] + sequence[bridge_end:]
    # Two intact candidate starts give distinct valid CDSs when both are published.
    sequence = 'ATG' + sequence
    bridge_start, bridge_end = bridge_start + 3, bridge_end + 3

    def coordinates(left, right):
        return [[left, right, 0]] if strand == '+' else [[len(sequence) - right, len(sequence) - left, 0]]

    genome, gff, fasta = (tmp_path / n for n in ('genome.fa', 'source.gff3', 'source.fa'))
    genomic = sequence if strand == '+' else str(Seq(sequence).reverse_complement())
    genome.write_text(f'>split\n{genomic}\n>reference\n{reference}\n')
    copy_gff = annotation('reference', 'copy', 'copy.t1', [[0, len(reference), 0]])
    if reference_case == 'generated_reference':
        copy_gff = copy_gff.replace('\tsynthetic\t', '\tGeneGalleon\t')
    gff.write_text('##gff-version 3\n'
                   + annotation('split', 'upstream', 'upstream.t1', coordinates(12, bridge_start), strand)
                   + annotation('split', 'downstream', 'downstream.t1', coordinates(bridge_end, len(sequence)), strand)
                   + copy_gff)
    fasta.write_text(f'>{name}_upstream\n{sequence[12:bridge_start]}\n'
                     f'>{name}_downstream\n{sequence[bridge_end:]}\n>{name}_copy\n{reference}\n')
    catalog = refinement.build_catalog(name, fasta, gff, genome)
    catalog_dir, db = tmp_path / 'catalog', tmp_path / 'loci.sqlite'
    refinement.write_catalog(catalog, catalog_dir)
    refinement.build_store([catalog_dir], db)
    selected = {locus['gene_id']: {'candidate_id': refinement._baseline(locus, locus['candidates'])['candidate_id']}
                for locus in catalog['loci']}

    prefix = name + '_'
    query = prefix + 'copy'
    with gff.open() as handle:
        owner_index = refinement.rescue.OwnershipIndex(
            refinement.rescue.original_ownership(handle, refinement.rescue.attributes))
    proposals = []
    for left in ([0, 3] if reference_case == 'competing_path' else [0]):
        shape = {'seqid': 'split', 'strand': strand, 'cds': coordinates(left, len(sequence))}
        proposals.append({**shape, 'candidate_id': prefix + f'join{left}',
                          'owner_overlap': refinement.rescue.split_model_geometry(
                              refinement.repair_coordinates(shape), owner_index),
                          'cds_sha256': hashlib.sha256(sequence[left:].encode()).hexdigest(),
                          'support': [{'donor': name, 'query': query, 'coverage': .99, 'identity': .95}]})
        assert len(proposals[-1]['owner_overlap']) == 2

    rescue_root = tmp_path / 'rescue'
    prepared = rescue_root / 'prepared' / name
    prepared.mkdir(parents=True)
    pep, mapping = prepared / 'genes.pep', prepared / 'genes.id_map.tsv'
    pep.write_text(f'>{query}\n{str(Seq(reference).translate()).removesuffix("*")}\n')
    locus = 'downstream' if reference_case == 'unmapped_reference' else 'copy'
    write_tsv(mapping, ['original_id', 'locus_id', 'jcvi_id', 'status'],
              [{'original_id': query, 'locus_id': locus, 'jcvi_id': query, 'status': 'selected'}])
    plan_sha = 'a' * 64
    (prepared / 'receipt.json').write_text(json.dumps({'key': {'plan': plan_sha, 'species': name},
        'files': {path.name: refinement.digest(path) for path in (pep, mapping)}}))
    request = {'parameters': refinement.DEFAULTS.copy(),
               'sources': {name: {'fasta': str(fasta), 'gff': str(gff), 'genome': str(genome), 'genetic_code': 1}},
               'automatic_split_merges': {name: proposals}}
    assert all(refinement.same_species_split_hits(row, name) for row in proposals)
    with refinement.ReceiptSnapshot(refinement._receipt_verifier(1)) as proof:
        proteins, loci = refinement.freeze_split_donors(
            request['automatic_split_merges'], None, rescue_root, plan_sha, request['sources'], proof)
    request.update(split_merge_donor_proteins=proteins, split_merge_same_species_loci=loci)
    new_genes, _, merges, _, _ = refinement.prepare_reviewed_models(request, name, db, selected)
    if reference_case != 'conserved':
        assert new_genes == merges == []
        return
    assert len(new_genes) == len(merges) == 1
    assert new_genes[0]['candidates'][0]['cds'] == sequence
    assert merges[0]['basis'] == 'same_species_copy_and_conserved_coding_bridge'
    assert set(merges[0]['replaced_gene_ids']) == {prefix + 'upstream', prefix + 'downstream'}
    support, = merges[0]['support']
    assert support['reference_gene_id'] == 'copy'
    assert support['coding_bridge']['identical_bridge_residues'] == 16
    assert support['coding_bridge']['candidate_interval'] == [bridge_start // 3, bridge_end // 3]


@pytest.mark.parametrize('case', ['partial', 'intact', 'reciprocal', 'manual', 'transcript', 'existing_transcript'])
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
    rows = annotation(contig, 'g8', 'g8.t1', source_blocks).splitlines()
    if case in {'transcript', 'existing_transcript'}:
        # Two mRNAs encode the same CDS but have different annotated UTRs.
        alt = annotation(contig, 'g8', 'g8.copy', source_blocks).splitlines()[1:]
        extended = []
        for line in [rows[0], *alt]:
            fields = line.split('\t')
            if fields[2] in {'gene', 'mRNA', 'exon'}:
                fields[3] = str(start - 5)
            extended.append('\t'.join(fields))
        rows = [extended[0], *rows[1:], *extended[1:],
                f'{contig}\tsynthetic\tfive_prime_UTR\t{start - 5}\t{start}\t.\t+\t.\tParent=g8.copy']
    with coding_gff.open('a') as handle:
        handle.write('\n'.join(rows) + '\n')
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
    if case == 'transcript':
        # The bounded prediction pass defers this locus. The frozen genome-wide
        # rescue projection still supports the completed CDS of both transcripts.
        refine_args = ('--candidate-limit', '1')
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
        status = 'conserved' if case == 'existing_transcript' and target == coding_target else 'repaired_model'
        assert mapping[target + f'_g{i}']['status'] == status
    assert mapping[phase_target + '_g11']['source_transcript_id'] == 'g11.t1'
    def phases(path):
        rows = (line.split('\t') for line in path.read_text().splitlines() if not line.startswith('#'))
        return [row[7] for row in rows if row[2] == 'CDS' and 'Parent=g11.t1' in row[8]]
    assert phases(Path(inputs[phase_target]['gff'])) == ['0', '2']
    assert phases(effective / 'full_annotation' / (phase_target + '.gff3')) == ['0', '2']
    assert phases(effective / 'source_annotation' / (phase_target + '.gff3')) == ['0', '0']
    audit = json.loads((effective / 'model_repairs.json').read_text())
    expected_audit = {phase_target + '_g11'}
    if case != 'existing_transcript':
        expected_audit.add(coding_target + '_g8')
    assert {row['gene_id'] for row in audit} == expected_audit
    if case == 'existing_transcript':
        tid = mapping[coding_target + '_g8']['source_transcript_id']
        full = (effective/'full_annotation'/(coding_target+'.gff3')).read_text()
        models = import_module('transcript_utrs').read_models(full, strict=False)
        assert {t for t,m in models.items() if m['gene'] == 'g8'} == {'g8.copy', 'g8.t1', tid}
    if case == 'transcript':
        record, = [x for x in audit if x['gene_id'] == coding_target + '_g8']
        assert record['derived_from_transcripts'] == ['g8.copy', 'g8.t1']
        assert record['representative']
        full = (effective/'full_annotation'/(coding_target+'.gff3')).read_text()
        for tid in ('g8.copy', 'g8.t1', record['source_transcript_id']):
            assert 'ID=' + tid + ';Parent=g8;' in full or 'ID=' + tid + ';Parent=g8\n' in full
        assert 'derived_from_transcript=g8.copy,g8.t1' in full
        assert 'five_prime_UTR' in full and 'Parent=g8.copy' in full
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
        if case == 'transcript':
            coding_basis = 'comparative_transcript_repair'
        expected_bases = {'exact_source_cds'} if case == 'existing_transcript' else {coding_basis, 'exact_source_cds'}
        assert {row['basis'] for row in audit} == expected_bases
        assert {row['adoption'] for row in audit} == {'automatic'}
    review = import_module('plot_gene_model_refinement').collect(output)
    assert review['model_repairs'] == audit
    assert all(path.read_bytes() == content for path, content in before.items())


def native_candidate_table(tmp_path, target, genome, gff_text, *, engine='annevo'):
    """Mock only the native predictor boundary, never genomic/donor evidence."""
    gff, receipt, table = (tmp_path / name for name in ('native.gff3', 'native.receipt.json', 'native.tsv'))
    gff.write_text(gff_text)
    value = dict(schema_version=1, kind='genegalleon-annevo-native-prediction', native_returncode=0,
                 gpu_probe={'gpu_verified': True}, model_profile='Magnoliopsida', runtime_digest='a' * 64,
                 model_sha256='b' * 64, source_commit='c' * 40, cache_key='d' * 64,
                 genome_sha256=hashlib.sha256(Path(genome).read_bytes()).hexdigest(),
                 GFF_sha256=hashlib.sha256(gff.read_bytes()).hexdigest(), GFF_bytes=gff.stat().st_size)
    if engine == 'oriongeno':
        value.update(kind='genegalleon-oriongeno-native-prediction', model_profile='plants',
                     inference_settings={'lineage': 'plants'})
    receipt.write_text(json.dumps(value))
    write_tsv(table, ['species', 'gff', 'receipt'], [dict(species=target, gff=str(gff), receipt=str(receipt))])
    return table, gff, receipt


@pytest.mark.parametrize('replay_rescue', [False, True])
def test_native_oriongeno_discovers_unpublished_split_models(tmp_path, replay_rescue):
    """Exercise discovery, paired residue support and publication from the CLI."""
    names = ['Species_target', 'Species_donor1', 'Species_donor2']
    targets = ('recover', 'only_one_fragment')
    sequences = {'recover': coding(313, residues=300),
                 'only_one_fragment': coding(314, residues=300)}
    cutoff = 72
    sequence = sequences['only_one_fragment']
    sequences['only_one_fragment'] = sequence[:cutoff] + 'ATG' + sequence[cutoff + 3:]
    sources, native_gff = [], ['##gff-version 3\n']
    for name in names:
        genome, annotation_rows, fasta = [], ['##gff-version 3\n'], []
        for locus in targets:
            sequence = sequences[locus]
            if locus == 'recover' and name != names[0]:
                # A homologous model can retain a long insertion between the two fragments.
                sequence = sequence[:420] + coding(315, residues=60)[:-3] + sequence[420:]
            dna, blocks = spliced_gene(sequence)
            models = [(locus, blocks, sequence)]
            if locus == 'only_one_fragment':
                right = [[blocks[0][0] + cutoff, blocks[0][1], 0], *blocks[1:]]
                models = [(locus, right, sequence[cutoff:])]
            if name == names[0]:
                if locus == 'recover':
                    models = [(locus + '_a', blocks[:1], sequence[:300]),
                              (locus + '_b', blocks[-1:], sequence[540:])]
                else:
                    first = [[blocks[0][0], blocks[0][0] + cutoff, 0]]
                    models = [(locus + '_a', first, sequence[:cutoff]),
                              (locus + '_b', right, sequence[cutoff:])]
                native_gff.append(annotation(locus, locus + '_native', locus + '_native.t1', blocks))
            for gene, cds, supplied in models:
                annotation_rows.append(annotation(locus, gene, gene + '.t1', cds))
                fasta.append(f'>{gene}.t1\n{supplied}\n')
            genome.append(f'>{locus}\n{dna}\n')
        paths = {key: tmp_path / f'{name}.{key}{suffix}'
                 for key, suffix in (('genome', '.fa'), ('cds', '.fa'), ('gff', '.gff3'))}
        for key, lines in (('genome', genome), ('cds', fasta), ('gff', annotation_rows)):
            paths[key].write_text(''.join(lines))
        sources.append({'species': name, 'genetic_code': 1, **paths})

    inputs, edges = tmp_path / 'inputs.tsv', tmp_path / 'edges.tsv'
    write_tsv(inputs, list(sources[0]), sources)
    write_tsv(edges, ['species_a', 'gene_a', 'species_b', 'gene_b', 'weight'], [])
    native, _, _ = native_candidate_table(tmp_path, names[0], sources[0]['genome'],
                                           ''.join(native_gff), engine='oriongeno')
    output = tmp_path / 'refinement'
    if replay_rescue:
        fixture = import_module('test_rescue_gene_models_runtime')
        for role in ('cds', 'gff', 'genome', 'busco'):
            (tmp_path / role).mkdir()
        for source in sources:
            name = source['species']
            for role, suffix in (('cds', '.cds.fa'), ('gff', '.gff3'), ('genome', '.genome.fa')):
                (tmp_path / role / (name + suffix)).symlink_to(source[role])
            (tmp_path / 'busco' / (name + '.busco.short.txt')).write_text(
                '# BUSCO version is: 6.0.0\n# The lineage dataset is: embryophyta_odb12\n'
                '# BUSCO was run in mode: transcriptome\nC:95.0%[S:95.0%,D:0.0%],F:0.0%,M:5.0%,n:100\n')
        (tmp_path / 'tree.nwk').write_text('((Species_target:1,Species_donor1:1):1,Species_donor2:2);\n')
        rescue_root, _, _ = fixture.make_plan((tmp_path, names, ()), model_storage='compact',
                                              extra_args=('--oriongeno-inputs', native, '--common-references', 3))
        for stage in ('prepare', 'synteny', 'rescue'):
            fixture.cli(stage, '--output', rescue_root, '--cpus', 2)
        name = names[0]
        worker = rescue_root / 'rescued' / name / 'receipt.json'
        before = worker.read_bytes()
        replay = rescue_root / 'native_split_candidates' / name
        fixture.cli('split-candidates', '--output', rescue_root, '--cpus', 2)
        assert worker.read_bytes() == before
        stamp = (replay / 'receipt.json').stat().st_mtime_ns
        assert 'Reused native_split_candidates/' + name in fixture.cli(
            'split-candidates', '--output', rescue_root, '--cpus', 2).stdout
        fixture.cli('finalize', '--output', rescue_root, '--cpus', 2)
        assert (replay / 'receipt.json').stat().st_mtime_ns == stamp
        assert worker.read_bytes() == before
        cli('plan', '--rescue-output', rescue_root, '--edges', edges, '--output', output)
        frozen = refinement.load(output)['request']
        assert frozen['native_split_panel']['donors'] == {name: names[1:]}
        row, = frozen['native_split_publications'][name]
        assert row['oriongeno_transcripts'] == ['recover_native.t1']
    else:
        cli('plan', '--inputs', inputs, '--edges', edges, '--oriongeno-inputs', native, '--output', output)
        frozen = refinement.load(output)['request']
        assert frozen['native_split_panel']['donors'] == {names[0]: names[1:]}
    assert 'automatic_split_merges' not in frozen
    cli('run', '--output', output, '--cpus', 2)
    cli('qc', '--output', output)
    effective = output / 'effective'
    accepted, = json.loads((effective / 'split_model_merges.json').read_text())
    assert accepted['basis'] == 'oriongeno_and_cross_fragment_protein_homology'
    assert accepted['oriongeno_transcripts'] == ['recover_native.t1']
    assert {hit['donor'] for hit in accepted['support']} == set(names[1:])
    assert all(hit['paired_coverage'] < .90 <= hit['coverage'] for hit in accepted['support'])
    assert {x['gene_id'] for x in accepted['owner_overlap']} == {'recover_a', 'recover_b'}
    assert sum(x['overlap_bases'] for x in accepted['owner_overlap']) < .75 * len(sequences['recover'])
    emitted = effective / 'analysis_cds' / (names[0] + '.fa')
    sequences_by_gene = {gene: value for gene, _, value in fasta_records(emitted)}
    assert sequences_by_gene[accepted['gene_id']] == sequences['recover']
    assert {names[0] + '_' + locus for locus in ('only_one_fragment_a', 'only_one_fragment_b')} <= sequences_by_gene.keys()


def test_native_oriongeno_splits_fused_models_only_with_collinear_support(tmp_path):
    target = 'Species_target'
    donors = [f'Species_donor{i:02}' for i in range(1, 13)]
    names = [target, *donors]
    accession_cases = ('accession_supported', 'accession_target_fusion', 'accession_crossing', 'accession_paralog')
    cases = ('recover', 'ten_split_two_fused', 'closest_fusion_only', 'close_donor_rna', 'close_mixed',
             *accession_cases,
             'one_collinear_pair', 'one_sister_branch', 'intact_rna', 'one_rna', 'junction_only', 'weak_junction',
             'duplicate_junction')
    fused_donors = {'ten_split_two_fused': donors[-2:], 'closest_fusion_only': donors[:2],
                    'close_donor_rna': donors[:1], 'close_mixed': donors[:1],
                    **{locus: donors[:1] for locus in accession_cases}}
    linked_donors = {'one_collinear_pair': donors[:1], 'one_sister_branch': donors[:2]}
    desired, source_paths, sources = {}, {}, []
    native_gff = ['##gff-version 3\n']
    for name in names:
        genome, gff, cds = [], ['##gff-version 3\n'], []
        for index, locus in enumerate(cases):
            first = coding(711 + index, residues=223)
            second = coding(815 + index, residues=120)
            intron = 'GT' + 'AC' * 12 + 'AG'
            # The fusion splices past the independent upstream stop exon.
            head = first[:-9]
            tail = first[-9:]
            sequence = 'AC' * 30 + head + intron + tail + intron + second + 'AC' * 30
            a = [60, 60 + len(head), 0]
            b = [a[1] + len(intron), a[1] + len(intron) + len(tail), 0]
            c = [b[1] + len(intron), b[1] + len(intron) + len(second), 0]
            strand = '-' if locus == 'recover' else '+'
            if strand == '-':
                size = len(sequence)
                a, b, c = ([size - end, size - start, 0] for start, end, _ in (a, b, c))
                sequence = str(Seq(sequence).reverse_complement())
            genome.append(f'>{locus}\n{sequence}\n')
            models = [(locus + '_first', [a, b], first), (locus + '_second', [c], second)]
            if name == names[0]:
                for gene, blocks, _ in models:
                    native_gff.append(annotation(locus, 'native_' + gene, 'native_' + gene + '.t', blocks, strand))
                desired[locus] = (first, second)
                source_paths[locus] = [a, c]
            if name == target or name in fused_donors.get(locus, ()):
                models = [(locus, [a, c], head + second)]
            elif locus in linked_donors and name not in linked_donors[locus]:
                # The short homolog has the right protein but is on an unlinked contig.
                scaffold = locus + '_unlinked'
                genome.append(f'>{scaffold}\n{sequence}\n')
                gff.append(annotation(scaffold, locus + '_second', locus + '_second.t', [c], strand))
                cds.append(f'>{locus}_second.t\n{second}\n')
                models = models[:1]
            for gene, blocks, supplied in models:
                gff.append(annotation(locus, gene, gene + '.t', blocks, strand))
                cds.append(f'>{gene}.t\n{supplied}\n')
        paths = {role: tmp_path / f'{name}.{role}' for role in ('genome', 'gff', 'cds')}
        for role, lines in (('genome', genome), ('gff', gff), ('cds', cds)):
            paths[role].write_text(''.join(lines))
        sources.append({'species': name, 'genetic_code': 1, **paths})
    inputs, edges, rna, tree = (tmp_path / label for label in ('inputs.tsv', 'edges.tsv', 'rna.tsv', 'tree.nwk'))
    write_tsv(inputs, list(sources[0]), sources)
    write_tsv(edges, ['species_a', 'gene_a', 'species_b', 'gene_b', 'weight'], [])
    write_tsv(rna, ['species', 'seqid', 'strand', 'cds_blocks', 'count', 'transcript_id'], [{
        'species': donors[0], 'seqid': 'close_donor_rna', 'strand': '+',
        'cds_blocks': json.dumps([block[:2] for block in source_paths['close_donor_rna']]),
        'count': 2, 'transcript_id': 'replicated_donor_fusion',
    }])
    subtrees = ['(' + ','.join(donors[a:b]) + ')' for a, b in ((0, 2), (2, 4), (4, 8), (8, 12))]
    tree.write_text('((((' + target + ',' + '),'.join(subtrees) + ');\n')

    gtf, protein, accessions = (tmp_path / name for name in ('accession.gtf', 'accession.faa', 'accessions.tsv'))
    exons, proteins = [], []
    for locus in accession_cases:
        for side, cds in enumerate(desired[locus]):
            start, end = 100 + side * 1000, 99 + side * 1000 + len(cds)
            gene = f'{locus}_{side}'
            attrs = f'gene_id "{gene}"; transcript_id "{gene}.t1";'
            exons.append(f'{locus}\tIsoSeq\texon\t{start}\t{end}\t.\t+\t.\t{attrs}\n')
            peptide = str(Seq(cds).translate()).removesuffix('*')
            proteins.append(f'>{gene}.t1\n{peptide}\n')
            if locus in ('accession_crossing', 'accession_supported'):
                crossed = f'gene_id "{locus}_both"; transcript_id "{locus}_both.t1";'
                label = 'IsoSeq' if locus == 'accession_crossing' else 'AbInitio'
                exons.append(f'{locus}\t{label}\texon\t{start}\t{end}\t.\t+\t.\t{crossed}\n')
            if locus == 'accession_paralog' and side == 0:
                exons.append(f'other\tIsoSeq\texon\t{start}\t{end}\t.\t+\t.\t'
                             'gene_id "paralog"; transcript_id "paralog.t1";\n')
                proteins.append(f'>paralog.t1\n{peptide}\n')
    gtf.write_text(''.join(exons))
    protein.write_text(''.join(proteins))
    write_tsv(accessions, refinement.refinement_accession_rna.FIELDS, [{
        'species': target, 'target_accession': 'Primary', 'evidence_accession': 'Secondary',
        'gtf': gtf.name, 'protein': protein.name, 'rna_sources': 'IsoSeq'}])
    rna_bams = tmp_path / 'rna_bam'
    rna_bams.mkdir()
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [
        {'SN': contig, 'LN': len(seq)} for contig, _, seq in fasta_records(sources[0]['genome'])]}
    entries = {'accession_target_fusion': ('target_fusion_a', 'target_fusion_b'),
               'intact_rna': ('molecule_a', 'molecule_b'), 'one_rna': ('one_molecule',),
               'junction_only': ('local_read1', 'local_read2', 'pcr_copy'),
               'weak_junction': ('short_anchor1', 'short_anchor2'),
               'duplicate_junction': ('same_extent_a', 'same_extent_b')}
    for library in ('a', 'b'):
        bam = rna_bams / f'{target}.{library}.bam'
        with pysam.AlignmentFile(bam, 'wb', header=header) as stream:
            for locus, ids in entries.items():
                (a, b, _), (c, d, _) = source_paths[locus]
                for i, identifier in enumerate(ids):
                    read = pysam.AlignedSegment(stream.header)
                    read.reference_id = stream.get_tid(locus)
                    left, right = (b - a, d - c)
                    if locus in ('junction_only', 'duplicate_junction'):
                        left, right = 20, 20 + (locus == 'junction_only' and i == 1)
                    elif locus == 'weak_junction':
                        left, right = 12, 8 + i
                    read.reference_start = b - left
                    read.cigartuples = [(0, left), (3, c - b), (0, right)]
                    read.query_name, read.query_sequence = identifier, 'A' * (left + right)
                    read.mapping_quality = 60
                    read.set_tag('XS', '+')
                    stream.write(read)
        pysam.index(str(bam))

    native, _, _ = native_candidate_table(tmp_path, target, sources[0]['genome'],
                                           ''.join(native_gff), engine='oriongeno')
    output = tmp_path / 'refinement'
    cli('plan', '--inputs', inputs, '--edges', edges, '--rna', rna, '--rna-bam-dir', rna_bams,
        '--species-tree', tree, '--oriongeno-inputs', native, '--accession-rna', accessions, '--output', output)
    cli('run', '--output', output, '--cpus', 2)
    cli('qc', '--output', output)
    effective = output / 'effective'
    published_splits = json.loads((effective / 'fused_model_splits.json').read_text())
    by_source = {row['source_gene_id']: row for row in published_splits}
    assert set(by_source) == {'recover', 'ten_split_two_fused', 'weak_junction', 'duplicate_junction',
                              'accession_supported'}
    accepted = by_source['recover']
    assert accepted['replaced_gene_id'] == target + '_recover'
    assert accepted['replaced_splice'] == [source_paths['recover'][1][1], source_paths['recover'][0][0]]
    assert {row['species'] for row in accepted['support']} == set(donors)
    conflict = by_source['ten_split_two_fused']
    assert len(conflict['support']) == 10
    assert {row['donor'] for row in conflict['fusion_annotation_support']} == set(donors[-2:])
    assert conflict['closest_informative_species'] == donors[:2]

    decisions = {row['source_gene_id']: row for row in json.loads(
        (effective / 'fused_model_split_decisions.json').read_text())}
    assert decisions['intact_rna']['reason'] == 'target_fusion_supported_by_rna'
    assert decisions['intact_rna']['rna_witnesses'] == 2
    assert decisions['one_rna']['reason'] == 'possible_target_fusion_one_rna_witness'
    assert decisions['one_rna']['rna_witnesses'] == 1
    assert decisions['junction_only']['reason'] == 'possible_target_fusion_replicated_junction'
    assert decisions['junction_only']['rna_witnesses'] == 0
    assert decisions['junction_only']['fusion_junction_rna_witnesses'] == 2
    assert decisions['duplicate_junction']['fusion_junction_rna_witnesses'] == 1
    assert decisions['close_mixed']['reason'] == 'mixed_closest_donor_architectures'
    assert decisions['accession_target_fusion']['reason'] == 'target_fusion_supported_by_rna'
    for locus in ('accession_crossing', 'accession_paralog'):
        assert decisions[locus]['reason'] == 'mixed_closest_donor_architectures'
    witness, = by_source['accession_supported']['accession_rna_support']
    assert witness['evidence_accession'] == 'Secondary'
    assert witness['transcripts'] == ['accession_supported_0.t1', 'accession_supported_1.t1']
    assert {hit['species'] for hit in by_source['accession_supported']['support']} == set(donors[1:])
    assert decisions['closest_fusion_only']['reason'] == 'closest_informative_donors_support_fusion'
    assert decisions['close_donor_rna']['reason'] == 'closest_donor_fusion_supported_by_rna'
    assert decisions['one_collinear_pair']['reason'] == 'insufficient_collinear_donor_species'
    assert decisions['one_sister_branch']['reason'] == 'split_support_on_one_sister_branch'

    new_ids = {model['gene_id'] for row in published_splits for model in row['models']}
    for view in ('analysis_cds', 'species_cds'):
        published = {gene: value for gene, _, value in fasta_records(effective / view / (target + '.fa'))}
        assert set(published) == new_ids | {target + '_' + locus for locus in set(cases) - by_source.keys()}
        for locus, split in by_source.items():
            assert tuple(published[model['gene_id']] for model in split['models']) == desired[locus]
    gff = (effective / 'species_gff' / (target + '.gff3')).read_text()
    for model, native_label in zip(accepted['models'], ('first', 'second'), strict=True):
        gene = model['gene_id']
        transcript = gene + '.t1'
        assert f'ID={transcript};Parent={gene};' in gff
        assert f'oriongeno_native_transcripts=native_recover_{native_label}.t' in gff
        assert f'Parent={transcript};gene_id={gene}' in gff
    assert 'Parent=recover;' not in gff
    archived = (effective / 'source_annotation' / (target + '.gff3')).read_text()
    assert archived == sources[0]['gff'].read_text()
    review = import_module('plot_gene_model_refinement').collect(output)
    assert review['fused_model_splits'] == published_splits
    assert review['species'][target]['fused_models_split'] == len(by_source)
    with protein.open('a') as handle:
        handle.write('\n')
    with pytest.raises(ValueError, match='Frozen refinement input changed'):
        refinement.load(output)


def test_native_terminal_repairs_and_complete_fragment_merges_use_distinct_evidence(tmp_path):
    target, *donors = ['Species_target', 'Species_donor1', 'Species_donor2']
    extensions = ('partial', 'complete_without_rna', 'complete_with_rna', 'donors_lack_extension')
    merges = ('joined_with_rna', 'joined_one_read', 'joined_competing_rna', 'joined_without_rna')
    sources, native_gff, rna, expected = [], ['##gff-version 3\n'], [], {}
    for name in (target, *donors):
        genome, gff, fasta = [], ['##gff-version 3\n'], []
        for i, locus in enumerate(extensions):
            sequence = coding(510 + i, residues=300)
            cut = 180
            sequence = sequence[:cut] + ('ATG' if locus.startswith('complete') else 'GAT') + sequence[cut + 3:]
            dna, blocks = spliced_gene(sequence)
            short = [[blocks[0][0] + cut, blocks[0][1], 0], *blocks[1:]]
            models = [(locus, short, sequence[cut:])] if name == target or locus == 'donors_lack_extension' else [
                (locus, blocks, sequence)]
            genome.append(f'>{locus}\n{dna}\n')
            for gene, cds, bases in models:
                gff.append(annotation(locus, gene, gene + '.t', cds))
                fasta.append(f'>{gene}.t\n{bases}\n')
            if name == target:
                native_gff.append(annotation(locus, 'native_' + locus, 'native_' + locus + '.t', blocks))
                expected[locus] = sequence
                if locus == 'complete_with_rna':
                    rna.append((locus, blocks, 2, 'extended_rna'))
        for i, locus in enumerate(merges):
            first, second = coding(380 + i, residues=220), coding(480 + i, residues=130)
            head, tail = first[:-9], first[-9:]
            intron = 'GT' + 'AC' * 16 + 'AG'
            dna = 'AC' * 30 + head + intron + tail + intron + second + 'AC' * 30
            a = [60, 60 + len(head), 0]
            b = [a[1] + len(intron), a[1] + len(intron) + len(tail), 0]
            c = [b[1] + len(intron), b[1] + len(intron) + len(second), 0]
            alternative, joined = [a, b], [a, c]
            genome.append(f'>{locus}\n{dna}\n')
            models = ([(locus + '_a', alternative, first), (locus + '_b', [c], second)]
                      if name == target else [(locus, joined, head + second)])
            for gene, cds, bases in models:
                gff.append(annotation(locus, gene, gene + '.t', cds))
                fasta.append(f'>{gene}.t\n{bases}\n')
            if name == target:
                native_gff.append(annotation(locus, 'native_' + locus, 'native_' + locus + '.t', joined))
                expected[locus] = head + second
                if locus != 'joined_without_rna':
                    rna.append((locus, joined, 1 if locus == 'joined_one_read' else 2, 'joined_rna'))
                if locus == 'joined_competing_rna':
                    rna.append((locus, alternative, 2, 'separate_rna'))
        paths = {role: tmp_path / f'{name}.{role}' for role in ('genome', 'gff', 'cds')}
        for role, lines in (('genome', genome), ('gff', gff), ('cds', fasta)):
            paths[role].write_text(''.join(lines))
        sources.append({'species': name, 'genetic_code': 1, **paths})
    inputs, edges, rna_path = (tmp_path / f'{label}.tsv' for label in ('inputs', 'edges', 'rna'))
    write_tsv(inputs, list(sources[0]), sources)
    write_tsv(edges, ['species_a', 'gene_a', 'species_b', 'gene_b', 'weight'], [])
    write_tsv(rna_path, ['species', 'seqid', 'strand', 'cds_blocks', 'count', 'transcript_id'], [
        {'species': target, 'seqid': locus, 'strand': '+', 'cds_blocks': json.dumps([v[:2] for v in blocks]),
         'count': count, 'transcript_id': transcript} for locus, blocks, count, transcript in rna])
    native, _, _ = native_candidate_table(tmp_path, target, sources[0]['genome'],
                                           ''.join(native_gff), engine='oriongeno')
    output = tmp_path / 'refinement'
    cli('plan', '--inputs', inputs, '--edges', edges, '--rna', rna_path, '--oriongeno-inputs', native, '--output', output)
    cli('run', '--output', output, '--cpus', 2)
    cli('qc', '--output', output)
    effective = output / 'effective'
    repairs = {row['owner_overlap'][0]['gene_id']: row for row in json.loads(
        (effective / 'model_repairs.json').read_text())}
    assert set(repairs) == {'partial', 'complete_with_rna'}
    assert repairs['partial']['basis'] == 'oriongeno_partial_cds_extension'
    assert repairs['complete_with_rna']['whole_cds_rna_witnesses'] == 2
    merge, = json.loads((effective / 'split_model_merges.json').read_text())
    assert merge['basis'] == 'oriongeno_cross_fragment_homology_and_target_rna'
    assert merge['whole_cds_rna_witnesses'] == 2
    assert merge['replaced_gene_ids'] == [target + '_joined_with_rna_' + label for label in ('a', 'b')]
    seqs = {name: sequence for name, _, sequence in fasta_records(effective / 'analysis_cds' / (target + '.fa'))}
    for locus in repairs:
        assert seqs[target + '_' + locus] == expected[locus]
    assert seqs[merge['gene_id']] == expected['joined_with_rna']
    assert all(target + '_' + locus + suffix in seqs for locus in merges[1:] for suffix in ('_a', '_b'))
    assert all(seqs[target + '_' + locus] == expected[locus][180:]
               for locus in ('complete_without_rna', 'donors_lack_extension'))


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


def test_fixed_panel_refinement_reuses_donors_and_isolates_genome_and_assembly_targets(tmp_path):
    fixture = import_module('test_rescue_gene_models_runtime')
    data = fixture.hidden_models.__wrapped__(tmp_path)
    _, names, truth = data
    target, *donors = names
    target_files = tmp_path / 'target_files'
    for role in ('cds', 'gff', 'genome', 'busco'):
        dst = target_files / role
        dst.mkdir(parents=True)
        source, = (tmp_path / role).glob(target + '.*')
        source.rename(dst / source.name)
    panel, _, _ = fixture.make_plan(data, directory_name='donors')
    fixture.cli('prepare', '--output', panel)
    rescue_root = tmp_path / 'target_rescue'
    fixture.cli('plan', '--cds-dir', target_files / 'cds', '--gff-dir', target_files / 'gff',
                '--genome-dir', target_files / 'genome', '--busco-dir', target_files / 'busco',
                '--tree', tmp_path / 'tree.nwk', '--output', rescue_root, '--donor-plan', panel)
    prepared = panel / 'prepared' / donors[0]
    receipt_path = prepared / 'receipt.json'
    receipt_before = receipt_path.read_bytes()
    changed = json.loads(receipt_before)
    changed['files']['unexpected.txt'] = hashlib.sha256(b'new').hexdigest()
    (prepared / 'unexpected.txt').write_bytes(b'new')
    refinement.atomic_json(receipt_path, changed)
    result = subprocess.run([sys.executable, str(Path(refinement.rescue.__file__)), 'synteny',
                             '--output', str(rescue_root), '--task-index', '1'],
                            text=True, capture_output=True)
    assert result.returncode and 'receipt changed' in result.stderr, result.stdout + result.stderr
    receipt_path.write_bytes(receipt_before)
    (prepared / 'unexpected.txt').unlink()
    fixture.cli('run', '--output', rescue_root, '--cpus', 2)
    rescue_plan = refinement.rescue.load(rescue_root)
    assert rescue_plan['species'] == [target] and set(rescue_plan['donors'][target]) == set(donors)
    genome_output = tmp_path / 'genomic'
    cli('plan', '--rescue-output', rescue_root, '--output', genome_output)
    cli('run', '--output', genome_output, '--cpus', 2)
    cli('qc', '--output', genome_output)
    genome_inputs = genome_output / 'effective/inputs.tsv'
    assert [r['species'] for r in refinement.verify_inputs(genome_inputs)] == [target]
    busco = import_module('gene_model_refinement_busco')
    pairs = busco.input_pairs(genome_output, tmp_path / 'cds', three_stage=True)
    assert [r['species'] for r in pairs] == [target]
    assert Path(pairs[0]['pre_rescue']).parent == target_files / 'cds'
    utr = tmp_path / 'genomic_utrs'
    cli('reproject-utrs', '--source-output', genome_output,
        '--source-receipt-sha256', refinement.digest(genome_output / 'effective/receipt.json'), '--output', utr)
    assert {p.stem for p in (utr / 'effective/species_gff').glob('*.gff3')} == {target}
    genomic_cds = genome_output / 'effective/species_cds' / (target + '.fa')
    effective_cds = {sequence for _, _, sequence in fasta_records(genomic_cds)}
    assert truth[8] in effective_cds
    assert not effective_cds.intersection({truth[9], truth[10]})
    assert set(p.name for p in (genome_output / 'predictions').iterdir()) == {target}
    receipt = genome_output / 'effective/receipt.json'
    before = (receipt.read_bytes(), receipt.stat().st_mtime_ns)
    cache = tmp_path / 'gene_model_refinement_donor_cache'
    donor_paths = sorted(cache.glob('catalog/*/receipt.json'))
    assert len(donor_paths) == len(donors)
    donor_mtimes = {p: p.stat().st_mtime_ns for p in donor_paths}

    assemblies = tmp_path / 'assembly' / 'assembled_transcripts_with_isoforms'
    assemblies.mkdir(parents=True)
    baseline = ''.join(truth[2][i:i + 3] if 60 <= i < 420 and i % 9 == 0 else truth[0][i:i + 3]
                       for i in range(0, len(truth[0]), 3))
    for name in ('Assembly_first', donors[0]):
        (assemblies / (name + '_isoform.fa')).write_text(
            f'>TRINITY_DN1_c0_g1_i1\n{baseline}\n'
            f'>TRINITY_DN1_c0_g1_i2\n{truth[0]}\n'
        )
    first_databases = None
    for name in ('Assembly_first', donors[0]):
        output = tmp_path / name
        cli('plan', '--donor-plan', panel, '--target', name,
            '--transcriptome-assembly-dir', assemblies.parent, '--output', output)
        recorded = json.loads((output / 'plan.json').read_text())['request']['fixed_donors']
        assert set(recorded) == set(donors) - {name}
        cli('run', '--output', output, '--cpus', 2)
        cli('qc', '--output', output)
        initial = output / 'selection_initial'
        final = output / 'selection_final'
        for filename in ('selection.json', 'representative_map.tsv'):
            assert (final / filename).read_bytes() == (initial / filename).read_bytes()
        assert json.loads((final / 'receipt.json').read_bytes())['key']['dependencies']['initial'] == refinement.digest(initial / 'receipt.json')
        sequence = output / 'effective/sequence'
        assert [r['species'] for r in refinement.read_table(sequence / 'inputs.tsv')] == [name]
        decision, = refinement.read_table(sequence / 'representative_map.tsv')
        assert decision['source_transcript_id'] == 'TRINITY_DN1_c0_g1_i2', decision
        assert not (output / 'effective/inputs.tsv').exists()
        assert set(p.name for p in (output / 'predictions').iterdir()) == {name}
        databases = {p: p.stat().st_mtime_ns for p in cache.glob('diamond/*/proteins.dmnd')}
        assert len(databases) == len(donors)
        if first_databases is not None:
            assert databases == first_databases
        first_databases = databases
    assert len(list(cache.glob('protein_pairs/*/edges.json'))) == len(donors) * 2 - 1
    assert (receipt.read_bytes(), receipt.stat().st_mtime_ns) == before
    assert {p: p.stat().st_mtime_ns for p in donor_paths} == donor_mtimes
    report = cli('run', '--output', genome_output, '--cpus', 2)
    assert 'Reused effective' in report.stdout
    assert (receipt.read_bytes(), receipt.stat().st_mtime_ns) == before


def test_fixed_refinement_reuses_exact_real_miniprot_search(tmp_path, monkeypatch):
    _, _, _, sources, truth = truth_fixture(tmp_path)
    source = sources[0]
    genome_path = Path(source['genome'])
    code = source.get('genetic_code', 1)
    source = {'genome': str(genome_path), 'genetic_code': code}
    query = str(Seq(truth[sources[0]['species'], 'intact']['cds']).translate()).removesuffix('*')
    donor = 'Donor_species'
    proteins = {donor: {'protein_1': query}}
    executable = refinement.shutil.which('miniprot')
    request = {'fixed_donors': [donor], 'donor_plan_path': str(tmp_path / 'panel/plan.json'),
               'files': {str(genome_path): refinement.digest(genome_path)},
               'miniprot': {'path': executable, 'sha256': refinement.digest(executable)},
               'implementation': refinement.implementation()}
    # The enclosing runtime test covers real plan loads and input fences.
    monkeypatch.setattr(refinement, 'load', lambda _root: request)
    real_search = refinement.rescue.search_intervals
    searches = []
    def count(*args, **kwargs):
        searches.append(kwargs['executable'])
        return real_search(*args, **kwargs)
    monkeypatch.setattr(refinement.rescue, 'search_intervals', count)

    with refinement.indexed_genome(genome_path) as genome:
        seqid = 'chr_intact'
        end = genome.get_reference_length(seqid)
        windows = {(seqid, 0, end): [{'id': 'query1', 'query': 'protein_1', 'donor': donor}]}
        params = {'max_intron': 10000}
        cold, warm = (tmp_path / stage for stage in ('cold', 'warm'))
        cold.mkdir()
        warm.mkdir()
        first = refinement.search_refinement_intervals(tmp_path, request, cold, windows, proteins,
                                                       genome, source, params, 1)
        second = refinement.search_refinement_intervals(tmp_path, request, warm, windows, proteins,
                                                        genome, source, params, 2)
    assert first and first == second and searches == [executable]
    for member in (refinement.rescue.INTERVAL_ARCHIVE, refinement.rescue.LOCAL_RECORDS):
        assert (cold / member).read_bytes() == (warm / member).read_bytes()


def test_posthoc_omark_mixed_cohort_reuses_final_proteomes(tmp_path):
    """Real rescue/refinement publications, stubbed external OMArk tools."""
    genome = import_module("test_genome_evolution_protein_mode")
    generation = import_module("test_gg_input_generation_end_to_end")
    fixture = import_module("test_rescue_gene_models_runtime")
    anchor, names, _ = fixture.make_plan(fixture.hidden_models.__wrapped__(tmp_path))
    fixture.cli("run", "--output", anchor, "--cpus", 2)
    assembly = "Plant_assembly"
    mapping = tmp_path / "isoforms.tsv"
    mapping.write_text("transcript_id\tgene_id\tbaseline\nshort\tPlant_assembly_g\t1\nlong\tPlant_assembly_g\t0\n")
    (tmp_path / "assembly.fa").write_text(">short\nATGAAATAA\n>long\nATGAAAAAATAA\n")
    sources = tmp_path / "assembly.tsv"
    sources.write_text(f"species\tcds\ttranscript_gene_map\n{assembly}\tassembly.fa\tisoforms.tsv\n")
    root = tmp_path / "refinement"
    refinement.plan(root, rescue_output=anchor, transcriptome_inputs=sources, mode="off", cpus=2)
    cli("run", "--output", root, "--cpus", 2)

    workspace = tmp_path / "workspace"
    (workspace / "tmp_runtime").mkdir(parents=True)
    tools = genome._prepare_stub_binaries(tmp_path)
    db = tmp_path / "fixture.h5"
    db.write_text("fixture database")
    env = generation._core_env(workspace, None, tools, "refinement_omark")
    env.update(gene_model_refinement_dir=str(root), omark_db_path=str(db))
    captured = tmp_path / "capture/omamer_queries.txt"
    report = root.with_name("refinement.review") / "omark_comparison.tsv"
    # Precompute one genome and the transcriptome even when effective/ is unavailable.
    effective = root / "effective"
    effective.rename(root / "pending")
    plan_species = json.loads((root / "plan.json").read_text())["species"]
    for index, name in enumerate([names[0], assembly], start=1):
        baseline_env = dict(env, input_generation_mode="refinement_omark_baseline",
                            GG_ARRAY_TASK_ID=str(plan_species.index(name) + 1))
        result = subprocess.run(["bash", generation.CORE_PATH], env=baseline_env, cwd=tmp_path, text=True, capture_output=True)
        assert result.returncode == 0, result.stdout + result.stderr
        assert f"{name} input protein_set_sha256=" in result.stdout
        assert len(captured.read_text().splitlines()) == index
    assert not report.exists()
    (root / "pending").rename(effective)
    for _ in range(2):
        result = subprocess.run(["bash", generation.CORE_PATH], env=env, cwd=tmp_path, text=True, capture_output=True)
        assert result.returncode == 0, result.stdout + result.stderr
        rows = {(r["species"], r["stage"]): r for r in refinement.read_table(report)}
        assert len(rows) == 3 * len(names) + 2
        assert rows[names[0], "input"]["sum_file"] != rows[names[0], "rescued"]["sum_file"]
        assert rows[names[0], "rescued"]["sum_file"] == rows[names[0], "refined"]["sum_file"]
        assert rows[names[1], "input"]["sum_file"] == rows[names[1], "rescued"]["sum_file"]
        assert rows[assembly, "input"]["sum_file"] == rows[assembly, "refined"]["sum_file"]
        assert (assembly, "rescued") not in rows
        assert int(rows[assembly, "refined"]["consistent_full"]) == 885
        assert len(captured.read_text().splitlines()) == len(names) + 2

    proteins = workspace / "input/species_protein"
    proteins.mkdir(parents=True)
    for row in refinement.read_table(root / "effective/sequence/inputs.tsv"):
        (proteins / (row["species"] + ".fa")).hardlink_to(root / "effective/sequence" / row["protein"])
    settings = {"run_species_omark": "1", "run_build_species_omark_summary": "1", "run_orthofinder": "0",
                "input_sequence_mode": "protein", "omark_db_path": str(db), "artifact_stale_policy": "rebuild"}
    result = genome._run_core(tmp_path, settings)
    assert result.returncode == 0, result.stdout + result.stderr
    assert len(captured.read_text().splitlines()) == len(names) + 2
    for name in [*names, assembly]:
        downstream = workspace / "output/genome_evolution/omark" / name / (name + ".sum")
        assert downstream.samefile(rows[name, "refined"]["sum_file"])

    changed = proteins / (assembly + ".fa")
    changed.unlink()
    changed.write_text(">Plant_assembly_g\nMTT\n")
    result = genome._run_core(tmp_path, settings)
    assert result.returncode == 0, result.stdout + result.stderr
    assert len(captured.read_text().splitlines()) == len(names) + 3
