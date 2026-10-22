"""Experimental preference must preserve DNA, coordinate proof and one locus."""
import gzip
import hashlib
import json
from pathlib import Path

import pytest
from Bio.Seq import Seq
from test_cloned_cds_models import cloned, fixture, fsi


def projected_fixture(root, strand='+', split=300):
    task, sequence, paths, row = fixture(root, strand=strand, split=split)
    from cloned_cds_models import source_graph
    _, _, coding, _ = source_graph(paths['gff'])
    blocks = sorted([[a, b, 0] for name in ('N_fragment', 'C_fragment') for a, b in coding[name]], reverse=strand == '-')
    cumulative = 0
    for block in blocks:
        block[2] = (3 - cumulative % 3) % 3
        cumulative += block[1] - block[0]
    name = 'experimental/clone:1'
    target = dict(seqid='chr1', strand=strand, cds=blocks, source_gene_ids=['N_fragment', 'C_fragment'])
    row.update(cloned_cds_mode='experimental', cloned_cds_projection_targets=json.dumps({name: target}))
    cloned.attach(task, row, root)
    return task, sequence, paths, row, target


def change_clone(task, paths, row, sequence, *, names=('experimental/clone:1',)):
    paths['clone'].write_text(''.join('>' + name + '\n' + sequence + '\n' for name in names))
    row['cloned_cds_sha256'] = hashlib.sha256(paths['clone'].read_bytes()).hexdigest()
    cloned.attach(task, row, paths['clone'].parent)


def format_pair(task, root):
    root.mkdir()
    cds = fsi.format_cds(task, root, False, False)
    gff = fsi.format_gff(task, root, False, False, formatted_cds_path=cds['output_path'])
    genome = fsi.format_genome(task, root, False, False)
    from experimental_cds_evidence import validate_pair
    assert validate_pair(cds['output_path'], gff['output_path'], genome['output_path'], 'Test_species') == 1
    from validate_cds_gff_mapping import validate_single_species
    result = validate_single_species(dict(index=1, species_prefix='Test_species', cds_file=cds['output_path'],
        gff_file=gff['output_path'], genome_file=genome['output_path'], strict=True), 10)
    assert result['ok'], result
    return cds['output_path'], gff['output_path'], genome['output_path']


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('split', [298, 299, 300])
@pytest.mark.parametrize('difference', ['synonymous', 'nonsynonymous', 'insertion', 'deletion'])
def test_native_format_and_exported_qc_preserve_clone_with_verified_projection(tmp_path, strand, split, difference):
    task, sequence, paths, row, _ = projected_fixture(tmp_path / 'inputs', strand, split)
    if difference == 'synonymous':
        for i in range(3, len(sequence) - 3, 3):
            candidate = sequence[:i+2] + ('C' if sequence[i+2] != 'C' else 'T') + sequence[i+3:]
            if str(Seq(candidate).translate()) == str(Seq(sequence).translate()):
                changed = candidate
                break
    elif difference == 'nonsynonymous':
        changed = sequence[:30] + ('AAA' if sequence[30:33] != 'AAA' else 'GCT') + sequence[33:]
    elif difference == 'insertion':
        changed = sequence[:150] + 'GCT' + sequence[150:]
    else:
        changed = sequence[:150] + sequence[153:]
    change_clone(task, paths, row, changed)
    before = {p: p.read_bytes() for p in paths.values()}
    cds, gff, _ = format_pair(task, tmp_path / 'formatted')
    assert changed in dict(fsi.iter_fasta_records(cds)).values()
    audit = json.loads(Path(task['cloned_model_audit']).read_text())
    assert audit['mode'] == 'experimental' and audit['reference_genome_modified'] is False
    assert audit['replacements'][0]['genome_identical'] is False
    assert audit['replacements'][0]['differences']
    assert audit['removed_CDS_ids'] == ['N_fragment.t1', 'C_fragment.t1']
    with gzip.open(gff, 'rt') as stream:
        text = stream.read()
    assert text.count('\tgene\t') == 2
    assert 'gg_sequence_origin=experimental_clone' in text
    assert all(p.read_bytes() == raw for p, raw in before.items())
    from validate_longest_cds_selection import collect_expected_longest_records
    expected, _ = collect_expected_longest_records(task)
    assert len(expected) == 2


def test_different_clones_for_one_gene_require_priority_and_publish_one_model(tmp_path):
    task, sequence, paths, row, target = projected_fixture(tmp_path / 'inputs')
    other = sequence[:30] + ('AAA' if sequence[30:33] != 'AAA' else 'GCT') + sequence[33:]
    paths['clone'].write_text('>experimental/clone:1\n' + sequence + '\n>preferred_clone\n' + other + '\n')
    row.update(cloned_cds_projection_targets=json.dumps({'experimental/clone:1': target, 'preferred_clone': target}),
               cloned_cds_sha256=hashlib.sha256(paths['clone'].read_bytes()).hexdigest())
    cloned.attach(task, row, tmp_path)
    with pytest.raises(ValueError, match='priority'):
        cloned.apply_to_task(task, tmp_path / 'failed')
    row['cloned_cds_priority_ids'] = '["preferred_clone"]'
    cloned.attach(task, row, tmp_path)
    cds, _, _ = format_pair(task, tmp_path / 'formatted')
    sequences = dict(fsi.iter_fasta_records(cds))
    assert len(sequences) == 2 and other in sequences.values() and sequence not in sequences.values()
    audit = json.loads(Path(task['cloned_model_audit']).read_text())
    assert len(audit['replacements']) == 1
    assert audit['replacements'][0]['aliases'] == ['experimental/clone:1']
    assert audit['skipped'][0]['status'] == 'same_gene_lower_priority'


def test_new_locus_is_explicit_and_does_not_remove_neighbouring_genes(tmp_path):
    task, sequence, paths, row, target = projected_fixture(tmp_path / 'inputs')
    text = paths['gff'].read_text().splitlines()
    paths['gff'].write_text('\n'.join(x for x in text if 'fragment' not in x) + '\n')
    paths['cds'].write_text('>neighbour.t1\nATGAAACCCTAA\n')
    target.update(source_gene_ids=[], new_locus=True)
    row['cloned_cds_projection_targets'] = json.dumps({'experimental/clone:1': target})
    cloned.attach(task, row, tmp_path)
    format_pair(task, tmp_path / 'formatted')
    audit = json.loads(Path(task['cloned_model_audit']).read_text())
    assert audit['replacements'][0]['new_locus'] is True
    assert audit['removed_CDS_ids'] == []


def test_longer_only_still_retains_longer_source_over_short_experimental_clone(tmp_path):
    task, sequence, paths, row, target = projected_fixture(tmp_path / 'inputs')
    paths['clone'].write_text('>experimental/clone:1\n' + sequence[:300] + 'TAA\n')
    row['cloned_cds_sha256'] = hashlib.sha256(paths['clone'].read_bytes()).hexdigest()
    row['cloned_cds_longer_only_ids'] = '["experimental/clone:1"]'
    # The supplied CDS is longer than neither an independently supplied full
    # source nor an assembled N/C locus. A full source length is authoritative.
    paths['cds'].write_text('>N_fragment.t1\n' + sequence + '\n>C_fragment.t1\n' + sequence[300:] + '\n>neighbour.t1\nATGAAACCCTAA\n')
    cloned.attach(task, row, tmp_path)
    cloned.apply_to_task(task, tmp_path / 'curated')
    audit = json.loads(Path(task['cloned_model_audit']).read_text())
    assert not audit['replacements'] and audit['skipped'][0]['status'] == 'not_longer'


@pytest.mark.parametrize('damage', ['coordinate', 'clone', 'reference', 'cigar', 'missing_genome', 'missing_proof', 'contig'])
def test_projection_proof_tampering_is_rejected_after_export(tmp_path, damage):
    task, sequence, paths, row, _ = projected_fixture(tmp_path / 'inputs')
    changed = sequence[:150] + 'GCT' + sequence[150:]
    change_clone(task, paths, row, changed)
    cds, gff, genome = format_pair(task, tmp_path / 'formatted')
    if damage == 'clone':
        records = list(fsi.iter_fasta_records(cds))
        assert changed in dict(records).values()
        with gzip.open(cds, 'wt') as f:
            f.write(''.join('>' + h + '\n' + (sequence if s == changed else s) + '\n' for h, s in records))
    elif damage == 'reference':
        records = list(fsi.iter_fasta_records(genome))
        assert any(sequence[:30] in s for _, s in records)
        with gzip.open(genome, 'wt') as f:
            f.write(''.join('>' + h + '\n' + s.replace(sequence[:30], 'A' * 30) + '\n' for h, s in records))
    elif damage == 'contig':
        for path in (gff, genome):
            with gzip.open(path, 'rt') as stream:
                text = stream.read()
            with gzip.open(path, 'wt') as stream:
                stream.write(text.replace('chr1', 'other_contig'))
    elif damage in {'coordinate', 'cigar', 'missing_proof'}:
        with gzip.open(gff, 'rt') as f:
            text = f.read()
        if damage == 'coordinate':
            text = text.replace('\t101\t400\t', '\t102\t400\t')
        elif damage == 'cigar':
            text = text.replace('gg_projection_cigar=', 'gg_projection_cigar=1D')
        else:
            text = text.replace('gg_projection_clone_sha256=', 'missing_projection_clone_sha256=')
        with gzip.open(gff, 'wt') as f:
            f.write(text)
    from experimental_cds_evidence import validate_pair
    with pytest.raises(ValueError):
        validate_pair(cds, gff, None if damage == 'missing_genome' else genome, 'Test_species')


def test_genomic_catalog_keeps_experimental_differences_separate_from_reference_translation(tmp_path):
    task, sequence, paths, row, _ = projected_fixture(tmp_path / 'inputs')
    changed = sequence[:30] + ('AAA' if sequence[30:33] != 'AAA' else 'GCT') + sequence[33:]
    change_clone(task, paths, row, changed)
    cds, gff, genome = format_pair(task, tmp_path / 'formatted')
    from gene_model_catalog import build_catalog
    result = build_catalog('Test_species', cds, gff, genome)
    candidates = [c for locus in result['loci'] for c in locus['candidates'] if c.get('experimental_cds_projection')]
    assert len(candidates) == 1
    assert candidates[0]['quality']['sequence_mismatch'] is True
    assert candidates[0]['quality']['usable'] is False
    assert candidates[0]['source_cds'][0]['cds'] == changed


def test_shorter_reviewed_partial_clone_preserves_full_reference_evidence(tmp_path):
    task, sequence, paths, row, target = projected_fixture(tmp_path / 'inputs')
    offset = sequence.index('ATG', 3, 290)
    target['cds'][0][0] += offset
    cumulative = 0
    for block in target['cds']:
        block[2] = (3 - cumulative % 3) % 3
        cumulative += block[1] - block[0]
    row['cloned_cds_projection_targets'] = json.dumps({'experimental/clone:1': target})
    change_clone(task, paths, row, sequence[offset:])
    with pytest.raises(ValueError, match='truncates'):
        cloned.apply_to_task(task, tmp_path / 'without_partial')
    target['partial_3prime'] = True
    row['cloned_cds_projection_targets'] = json.dumps({'experimental/clone:1': target})
    cloned.attach(task, row, tmp_path)
    with pytest.raises(ValueError, match='truncates'):
        cloned.apply_to_task(task, tmp_path / 'wrong_partial_end')
    target.pop('partial_3prime')
    target['partial_5prime'] = True
    row['cloned_cds_projection_targets'] = json.dumps({'experimental/clone:1': target})
    cloned.attach(task, row, tmp_path)
    format_pair(task, tmp_path / 'formatted')
    audit_path = Path(task['cloned_model_audit'])
    audit = json.loads(audit_path.read_text())
    assert audit['replacements'][0]['projection']['partial_5prime'] == 'true'
    removed = (audit_path.parent / audit['outputs']['removed_cds']).read_text()
    assert sequence[:300] in removed and sequence[300:] in removed
    assert audit['outputs']['removed_gff']


def test_native_summary_prevents_silent_projection_metadata_removal(tmp_path):
    task, sequence, paths, row, _ = projected_fixture(tmp_path / 'inputs')
    cds, gff, genome = format_pair(task, tmp_path / 'formatted')
    with gzip.open(gff, 'rt') as stream:
        lines = stream.read().splitlines()
    with gzip.open(gff, 'wt') as stream:
        for line in lines:
            if not line.startswith('#'):
                fields = line.split('\t')
                fields[8] = ';'.join(a for a in fields[8].split(';') if not a.startswith(('gg_projection_', 'gg_sequence_origin=')))
                line = '\t'.join(fields)
            stream.write(line + '\n')
    from validate_cds_gff_mapping import validate_single_species
    result = validate_single_species(dict(index=1, species_prefix='Test_species', cds_file=cds, gff_file=gff,
        genome_file=genome, strict=True, summary_row=dict(cloned_cds_mode='experimental',
        cloned_model_audit_path=task['cloned_model_audit'])), 10)
    assert result['ok'] is False


@pytest.mark.parametrize('damage', ['unreviewed', 'bad_phase', 'wrong_gene', 'new_locus_over_gene'])
def test_invalid_reviewed_projection_does_not_publish(tmp_path, damage):
    task, sequence, paths, row, target = projected_fixture(tmp_path / 'inputs')
    if damage == 'unreviewed':
        row['cloned_cds_projection_targets'] = '{}'
    elif damage == 'bad_phase':
        target['cds'][1][2] = 1
    elif damage == 'wrong_gene':
        target['source_gene_ids'] = ['neighbour']
    else:
        target.update(source_gene_ids=[], new_locus=True)
    if damage != 'unreviewed':
        row['cloned_cds_projection_targets'] = json.dumps({'experimental/clone:1': target})
    with pytest.raises(ValueError):
        cloned.attach(task, row, tmp_path)
        cloned.apply_to_task(task, tmp_path / 'curated')


@pytest.mark.parametrize('reference_damage', ['stop', 'extra_base'])
def test_intact_clone_can_replace_disrupted_reference_with_explicit_difference_proof(tmp_path, reference_damage):
    task, sequence, paths, row, target = projected_fixture(tmp_path / 'inputs')
    if reference_damage == 'stop':
        genome = paths['genome'].read_text()
        changed = sequence[:30] + 'TAA' + sequence[33:300]
        paths['genome'].write_text(genome.replace(sequence[:300], changed))
    else:
        # A reference-only base makes the projected CDS non-triplet length.
        target['cds'][-1][1] += 1
        row['cloned_cds_projection_targets'] = json.dumps({'experimental/clone:1': target})
        cloned.attach(task, row, tmp_path)
    format_pair(task, tmp_path / 'formatted')
    audit = json.loads(Path(task['cloned_model_audit']).read_text())
    assert audit['replacements'][0]['differences']
    assert audit['replacements'][0]['genome_identical'] is False
