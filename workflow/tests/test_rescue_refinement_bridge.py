"""Missing-gene nominations retain genomic QC and exact existing-locus trust."""
import copy
import json
import sqlite3
import sys
from importlib import import_module
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / 'support'
sys.path.insert(0, str(SUPPORT))
refinement = import_module('gene_model_refinement')


def fixture(*, intact=False):
    sequence = 'ATG' + 'AAA' * 18 + 'TAA'
    original = {'candidate_id': 'Species_target_t', 'source_transcript_id': 't', 'source_gene_id': 'g',
                'seqid': 'chr1', 'strand': '+', 'blocks': [[100, 130, 0]], 'cds': sequence[:30],
                'origin': 'original', 'quality': {'valid_orf': intact}}
    locus = {'gene_id': 'Species_target_g', 'species': 'Species_target', 'seqid': 'chr1', 'strand': '+',
             'candidates': [original]}
    catalog = {'species': 'Species_target', 'genetic_code': 1, 'loci': [locus]}
    models = [dict(gene_id='Species_target_g', seqid='chr1', strand='+', cds=[[100, 160, 0]],
                   sequence=sequence, problems=[], donor_species=donor, donor_candidate=donor + '_t',
                   identity=1.0, coverage=1.0) for donor in ('Species_donor1', 'Species_donor2')]
    edges = [dict(species_a='Species_target', gene_a='Species_target_g', species_b=donor,
                  gene_b=donor + '_g', weight=1.0, ambiguous=False)
             for donor in ('Species_donor1', 'Species_donor2')]
    return catalog, models, edges


def classify(catalog, models, edges, rna=()):
    return refinement.classify_predictions(models, catalog, edges, refinement.DEFAULTS, list(rna), 'assembly')


class Genome:
    def __init__(self, sequence):
        self.sequence = 'C' * 100 + sequence + 'C' * 100

    def get_reference_length(self, seqid):
        assert seqid == 'chr1'
        return len(self.sequence)

    def fetch(self, seqid, start, end):
        assert seqid == 'chr1'
        return self.sequence[start:end]


def rescue_model(models, *, extra_metrics=True):
    raw = copy.deepcopy(models[0])
    raw.update(frameshift=False, model_id='rescue-shape', revision_owner_ids=['g'],
               problems=['overlap_existing_annotation'], evidence={'donor': 'Species_donor1', 'query': 'q1'})
    extra = {'donor': 'Species_donor2', 'query': 'q2'}
    if extra_metrics:
        extra['alignment'] = {'coverage': 1.0, 'identity': 1.0, 'problems': []}
    raw['support'] = [raw['evidence'], extra]
    return raw


def resolve(donor, query):
    return (donor + '_g', donor + '_t') if query == {'Species_donor1': 'q1', 'Species_donor2': 'q2',
                                                    'Species_donor3': 'q3'}.get(donor) else None


@pytest.mark.parametrize('strand,accepted', [('+', False), ('-', True), ('.', False), ('?', False)])
def test_foreign_coding_overlap_is_strand_aware_and_unknown_is_conservative(strand, accepted):
    catalog, models, edges = fixture()
    neighbor = copy.deepcopy(catalog['loci'][0])
    neighbor.update(gene_id='Species_target_neighbor', strand=strand)
    neighbor['candidates'][0].update(strand=strand, blocks=[[150, 180, 0]], source_gene_id='neighbor')
    catalog['loci'].append(neighbor)
    row = classify(catalog, models, edges)[0]
    assert (row['status'] == 'accepted') is accepted
    assert ('overlap_other_locus' in row['problems']) is not accepted


@pytest.mark.parametrize('strand,accepted', [('+', False), ('-', True), ('.', False), ('?', False)])
def test_noncoding_foreign_span_keeps_strand_and_unknown_protection(tmp_path, strand, accepted):
    catalog, models, edges = fixture()
    gff = tmp_path / 'noncoding.gff3'
    gff.write_text(f'chr1\tp\tgene\t151\t180\t.\t{strand}\t.\tID=nc\n'
                   f'chr1\tp\tncRNA\t151\t180\t.\t{strand}\t.\tID=nct;Parent=nc\n')
    catalog['annotation_spans'] = refinement.annotation_ownership_spans(gff, catalog)
    assert {r['strand'] for r in catalog['annotation_spans']} == {strand if strand in {'+', '-'} else '.'}
    assert all('cds_blocks' not in r for r in catalog['annotation_spans'])
    row = classify(catalog, models, edges)[0]
    assert (row['status'] == 'accepted') is accepted


def test_known_foreign_coding_intron_does_not_claim_an_independent_gene(tmp_path):
    catalog, models, edges = fixture()
    neighbor = copy.deepcopy(catalog['loci'][0])
    neighbor.update(gene_id='Species_target_neighbor')
    neighbor['candidates'][0].update(blocks=[[80, 95, 0], [180, 201, 0]], source_gene_id='neighbor', source_transcript_id='nt')
    catalog['loci'].append(neighbor)
    gff = tmp_path / 'coding.gff3'
    gff.write_text('chr1\tp\tgene\t81\t201\t.\t+\t.\tID=neighbor\n'
                   'chr1\tp\tmRNA\t81\t201\t.\t+\t.\tID=nt;Parent=neighbor\n'
                   'chr1\tp\tCDS\t81\t95\t.\t+\t0\tParent=nt\n'
                   'chr1\tp\tCDS\t181\t201\t.\t+\t0\tParent=nt\n')
    catalog['annotation_spans'] = refinement.annotation_ownership_spans(gff, catalog)
    assert all(r.get('cds_blocks') for r in catalog['annotation_spans'])
    assert classify(catalog, models, edges)[0]['status'] == 'accepted'
    for model in models:
        model['cds'] = [[140, 200, 0]]
    row = classify(catalog, models, edges)[0]
    assert row['status'] == 'proposal' and 'overlap_other_locus' in row['problems']


def test_opposite_strand_predictions_are_not_competing_coding_paths():
    catalog, models, edges = fixture()
    second = copy.deepcopy(catalog['loci'][0])
    second.update(gene_id='Species_target_second', strand='-')
    second['candidates'][0].update(source_gene_id='second', source_transcript_id='st', strand='-', blocks=[[180, 210, 0]])
    catalog['loci'].append(second)
    alternatives = [{**m, 'gene_id': second['gene_id'], 'strand': '-', 'cds': [[150, 210, 0]]} for m in models]
    edges += [{**edge, 'gene_a': second['gene_id']} for edge in edges]
    rows = classify(catalog, models + alternatives, edges)
    assert len(rows) == 2 and all(row['status'] == 'accepted' for row in rows)


def test_imported_two_donor_revision_runs_genomic_validation_and_normal_gates():
    catalog, models, edges = fixture()
    raw = rescue_model(models)
    imported, proposals = refinement.import_rescue_revisions([raw], catalog, edges, refinement.DEFAULTS,
                                                            Genome(raw['sequence']), resolve)
    assert len(imported) == 2 and proposals == []
    assert raw['problems'] == ['overlap_existing_annotation']
    row = classify(catalog, imported, edges)[0]
    assert row['status'] == 'accepted'
    assert row['change_type'] == 'model_revision'
    assert row['donors'] == ['Species_donor1', 'Species_donor2']
    assert row['candidate']['quality']['representative_eligible']
    assert row['gene_id'] == 'Species_target_g'


def test_unmeasured_extra_rescue_support_cannot_inflate_independent_donor_count():
    catalog, models, edges = fixture()
    raw = rescue_model(models, extra_metrics=False)
    imported, proposals = refinement.import_rescue_revisions([raw], catalog, edges, refinement.DEFAULTS,
                                                            Genome(raw['sequence']), resolve)
    assert len(imported) == 1
    assert proposals[0]['reason'] == 'rescue_support_without_alignment_metrics'
    row = classify(catalog, imported, edges)[0]
    assert row['status'] == 'proposal' and row['donors'] == ['Species_donor1']
    assert 'insufficient_independent_support' in row['problems']


def test_weak_primary_cannot_veto_two_independently_measured_rescue_donors():
    catalog, models, edges = fixture()
    raw = rescue_model(models)
    raw.update(coverage=0.5, problems=['overlap_existing_annotation', 'low_coverage'])
    raw['support'].append({'donor': 'Species_donor3', 'query': 'q3',
                           'alignment': {'coverage': 1.0, 'identity': 1.0, 'problems': []}})
    edges.append({**edges[0], 'species_b': 'Species_donor3', 'gene_b': 'Species_donor3_g'})
    imported, proposals = refinement.import_rescue_revisions([raw], catalog, edges, refinement.DEFAULTS,
                                                            Genome(raw['sequence']), resolve)
    row = classify(catalog, imported, edges)[0]
    assert not proposals and row['status'] == 'accepted'
    assert row['donors'] == ['Species_donor2', 'Species_donor3']


@pytest.mark.parametrize('case', ['absent_owner', 'multiple_owners', 'ambiguous_owner'])
def test_import_does_not_choose_or_invent_an_existing_owner(case):
    catalog, models, edges = fixture()
    raw = rescue_model(models)
    if case == 'absent_owner':
        raw['revision_owner_ids'] = ['missing']
    elif case == 'multiple_owners':
        raw['revision_owner_ids'] = ['g', 'other']
    else:
        duplicate = copy.deepcopy(catalog['loci'][0])
        duplicate['gene_id'] = 'Species_target_other'
        catalog['loci'].append(duplicate)
    imported, proposals = refinement.import_rescue_revisions([raw], catalog, edges, refinement.DEFAULTS,
                                                            Genome(raw['sequence']), resolve)
    assert imported == []
    assert proposals[0]['reason'] == 'ambiguous_or_absent_rescue_revision_owner'


def test_donor_species_alone_does_not_establish_exact_query_correspondence():
    catalog, models, edges = fixture()
    raw = rescue_model(models)
    def wrong(donor, query):
        return donor + '_paralog', donor + '_paralog_t'
    imported, proposals = refinement.import_rescue_revisions([raw], catalog, edges, refinement.DEFAULTS,
                                                            Genome(raw['sequence']), wrong)
    assert imported == [] and len(proposals) == 2
    assert all(p['reason'] == 'untrusted_rescue_donor_query_correspondence' for p in proposals)


@pytest.mark.parametrize('problem', ['frameshift', 'outside_search_window', 'noncanonical_splice', 'assembly_gap_or_ambiguity'])
def test_import_preserves_original_hard_failure_even_when_current_dna_is_intact(problem):
    catalog, models, edges = fixture()
    raw = rescue_model(models)
    raw['problems'].append(problem)
    imported, _ = refinement.import_rescue_revisions([raw], catalog, edges, refinement.DEFAULTS,
                                                   Genome(raw['sequence']), resolve)
    row = classify(catalog, imported, edges)[0]
    assert row['status'] == 'proposal' and problem in row['problems']


def test_import_cannot_replace_mismatched_dna_or_protected_annotation():
    catalog, models, edges = fixture()
    raw = rescue_model(models)
    raw['sequence'] = raw['sequence'].replace('AAA', 'CCC')
    imported, _ = refinement.import_rescue_revisions([raw], catalog, edges, refinement.DEFAULTS,
                                                   Genome(models[0]['sequence']), resolve)
    row = classify(catalog, imported, edges)[0]
    assert row['status'] == 'proposal' and 'imported_rescue_sequence_mismatch' in row['problems']
    raw = rescue_model(models)
    catalog['loci'][0]['candidates'][0]['quality']['annotated_exception'] = 'pseudogene'
    imported, _ = refinement.import_rescue_revisions([raw], catalog, edges, refinement.DEFAULTS,
                                                   Genome(raw['sequence']), resolve)
    row = classify(catalog, imported, edges)[0]
    assert row['status'] == 'proposal' and 'protected_annotation_exception' in row['problems']


def test_imported_intact_locus_addition_retains_target_rna_adoption_gate():
    catalog, models, edges = fixture(intact=True)
    raw = rescue_model(models)
    imported, _ = refinement.import_rescue_revisions([raw], catalog, edges, refinement.DEFAULTS,
                                                   Genome(raw['sequence']), resolve)
    row = classify(catalog, imported, edges)[0]
    assert row['status'] == 'accepted' and row['change_type'] == 'isoform_addition'
    assert not row['candidate']['quality']['representative_eligible']
    rna = [{'species': 'Species_target', 'seqid': 'chr1', 'strand': '+', 'cds_blocks': '[[100,160]]',
            'transcript_id': 'whole-path', 'count': 1}]
    assert classify(catalog, imported, edges, rna)[0]['candidate']['quality']['representative_eligible']


def test_gene_level_primary_fasta_does_not_mark_rescue_alternative_as_mismatched(tmp_path):
    genome, gff, cds = [tmp_path / suffix for suffix in ('genome.fa', 'models.gff3', 'cds.fa')]
    genome.write_text('>chr1\nATGAAAGTAAAAAGCCCGGGCCCTAA\n')
    cds.write_text('>rescued\nATGAAACCCGGGCCCTAA\n')
    gff.write_text('##gff-version 3\n'
                   'chr1\tgenegalleon_rescue\tgene\t1\t26\t.\t+\t.\tID=rescued\n'
                   'chr1\tgenegalleon_rescue\tmRNA\t1\t26\t.\t+\t.\tID=primary;Parent=rescued\n'
                   'chr1\tgenegalleon_rescue\tCDS\t1\t6\t.\t+\t0\tParent=primary\n'
                   'chr1\tgenegalleon_rescue\tCDS\t15\t26\t.\t+\t0\tParent=primary\n'
                   'chr1\tgenegalleon_rescue\tmRNA\t1\t26\t.\t+\t.\tID=alternative;Parent=rescued;support=homology_coding_path\n'
                   'chr1\tgenegalleon_rescue\tCDS\t1\t6\t.\t+\t0\tParent=alternative\n'
                   'chr1\tgenegalleon_rescue\tCDS\t21\t26\t.\t+\t0\tParent=alternative\n')
    catalog = refinement.build_catalog('Species_target', cds, gff, genome)
    locus, = catalog['loci']
    assert locus['source_baseline_candidate_id'] == 'Species_target_primary'
    assert len(locus['candidates']) == 2
    assert all(c['quality']['usable'] and not c['quality']['sequence_mismatch'] for c in locus['candidates'])
    assert {c['cds'] for c in locus['candidates']} == {'ATGAAACCCGGGCCCTAA', 'ATGAAACCCTAA'}
    alternative, = [c for c in locus['candidates'] if c.get('rescue_alternative_coding_path')]
    assert alternative['origin'] == 'predicted' and not alternative['quality']['representative_eligible']
    assert catalog['summary']['rescue_alternative_coding_paths'] == 1
    inputs, edges = tmp_path / 'inputs.tsv', tmp_path / 'edges.tsv'
    inputs.write_text(f'species\tcds\tgff\tgenome\tgenetic_code\nSpecies_target\t{cds}\t{gff}\t{genome}\t1\n')
    edges.write_text('species_a\tgene_a\tspecies_b\tgene_b\n')
    rna = tmp_path / 'rna.tsv'
    rna.write_text('species\tseqid\tstrand\tcds_blocks\ttranscript_id\tcount\n'
                   'Species_target\tchr1\t+\t[[0,6],[20,26]]\tfull-rna\t1\n')
    for mode, evidence, admitted in [('rna_required', None, False), ('conservation_supported', None, True),
                                    ('rna_required', rna, True)]:
        root = tmp_path / ('catalog_' + mode + ('_rna' if evidence else ''))
        value = refinement.plan(root, inputs=inputs, edges=edges, rna=evidence, mode='off', isoform_adoption=mode)
        directory = refinement.catalog_species(root, value, 'Species_target')
        row = json.loads((directory / 'loci.jsonl').read_text())
        candidate, = [c for c in row['candidates'] if c.get('rescue_alternative_coding_path')]
        assert candidate['quality']['representative_eligible'] is admitted
        assert (candidate['quality'].get('representative_admission') == 'conservation_supported') is (mode == 'conservation_supported')
        effective = refinement.finalize(root, value)
        full = refinement.build_catalog('Species_target', effective / 'species_cds/Species_target.fa',
                                        effective / 'full_annotation/Species_target.gff3', genome)
        assert len(full['loci'][0]['candidates']) == 2
        assert all(c['quality']['usable'] for c in full['loci'][0]['candidates'])
        assert refinement.verify_inputs(effective / 'inputs.tsv')
    # A marker supplied through explicit inputs cannot waive a partial ORF,
    # even when target RNA follows exactly that incomplete coding path.
    gff.write_text(gff.read_text().replace('CDS\t21\t26', 'CDS\t21\t23'))
    rna.write_text(rna.read_text().replace('[[0,6],[20,26]]', '[[0,6],[20,23]]'))
    for mode, evidence in [('conservation_supported', None), ('rna_required', rna)]:
        root = tmp_path / ('partial_' + mode)
        value = refinement.plan(root, inputs=inputs, edges=edges, rna=evidence, mode='off', isoform_adoption=mode)
        directory = refinement.catalog_species(root, value, 'Species_target')
        locus = json.loads((directory / 'loci.jsonl').read_text())
        candidate, = [c for c in locus['candidates'] if c.get('rescue_alternative_coding_path')]
        assert candidate['quality']['partial'] and not candidate['quality']['valid_orf']
        assert not candidate['quality']['representative_eligible']


def test_refinement_prediction_completes_genomic_terminus_before_normal_acceptance(tmp_path, monkeypatch):
    # Search results are supplied below, so this unit test must not depend on
    # the optional predictor being installed on the fast-test runner.
    predictor = tmp_path / 'miniprot'
    predictor.write_text('#!/bin/sh\nexit 99\n')
    predictor.chmod(0o755)
    original_which = refinement.shutil.which
    monkeypatch.setattr(refinement.shutil, 'which',
                        lambda name: str(predictor) if name == 'miniprot' else original_which(name))
    sequence = 'ATGAAACCCGGGCCCTAA'
    names = ['Species_target', 'Species_donor1', 'Species_donor2']
    inputs, edges = tmp_path / 'inputs.tsv', tmp_path / 'edges.tsv'
    rows = []
    for name in names:
        cds, gff, genome = [tmp_path / (name + suffix) for suffix in ('.cds.fa', '.gff3', '.genome.fa')]
        start = 3 if name == 'Species_target' else 0
        cds.write_text('>g\n' + sequence[start:] + '\n')
        genome.write_text('>chr1\n' + sequence + '\n')
        gff.write_text(f'##gff-version 3\nchr1\ts\tgene\t{start + 1}\t18\t.\t+\t.\tID=g\n'
                       f'chr1\ts\tmRNA\t{start + 1}\t18\t.\t+\t.\tID=t;Parent=g\n'
                       f'chr1\ts\tCDS\t{start + 1}\t18\t.\t+\t0\tParent=t\n')
        rows.append('\t'.join([name, str(cds), str(gff), str(genome), '1']))
    inputs.write_text('species\tcds\tgff\tgenome\tgenetic_code\n' + '\n'.join(rows) + '\n')
    edges.write_text('species_a\tgene_a\tspecies_b\tgene_b\n' + ''.join(
        f'Species_target\tSpecies_target_g\t{donor}\t{donor}_g\n' for donor in names[1:]))
    def local_search(tmp, windows, proteins, genome, code, params, cpus):
        return [dict(query=region['id'], seqid=region['id'], strand='+', cds=[[3, 18, 0]], frameshift=False,
                     coverage=.8, identity=1, query_start=1, query_end=5, query_length=5)
                for regions in windows.values() for region in regions]
    monkeypatch.setattr(refinement.rescue, 'search_intervals', local_search)
    root = tmp_path / 'refinement'
    value = refinement.plan(root, inputs=inputs, edges=edges)
    directory = refinement.predict_species(root, value, 'Species_target')
    predicted, = json.loads((directory / 'predictions.json').read_text())
    assert predicted['status'] == 'accepted' and predicted['change_type'] == 'model_revision'
    assert predicted['candidate']['cds'] == sequence
    assert predicted['candidate']['blocks'] == [[0, 18, 0]]
    assert predicted['candidate']['quality']['representative_eligible']
    assert all(a['terminal_completion']['status'] == 'completed' for a in predicted['alignments'])


def test_rescue_query_resolution_requires_usable_exact_frozen_protein(tmp_path):
    mapping, protein = tmp_path / 'mapping.tsv', tmp_path / 'genes.pep'
    mapping.write_text('jcvi_id\tlocus_id\toriginal_id\tstatus\nq1\tg\tSpecies_donor1_g\tselected\n')
    protein.write_text('>q1\nMK\n')
    request = {'revision_donor_maps': {'Species_donor1': {'mapping': str(mapping), 'protein': str(protein)}}}
    candidate = {'candidate_id': 'Species_donor1_t', 'protein': 'MK', 'quality': {'usable': True},
                 'source_fasta_ids': ['Species_donor1_g'], 'coding_key': 'shape'}
    locus = {'gene_id': 'Species_donor1_g', 'source_baseline_candidate_id': candidate['candidate_id'], 'candidates': [candidate]}
    connection = sqlite3.connect(':memory:')
    connection.execute('PRAGMA user_version=1')
    connection.execute('CREATE TABLE loci(species TEXT,gene_id TEXT,json TEXT)')
    connection.execute('INSERT INTO loci VALUES (?,?,?)', ('Species_donor1', locus['gene_id'], json.dumps(locus)))
    resolver = refinement.rescue_revision_donor_resolver(request, connection)
    assert resolver('Species_donor1', 'q1') == ('Species_donor1_g', 'Species_donor1_t')
    assert resolver('Species_donor1', 'absent') is None
    candidate['quality']['usable'] = False
    connection.execute('UPDATE loci SET json=?', (json.dumps(locus),))
    assert refinement.rescue_revision_donor_resolver(request, connection)('Species_donor1', 'q1') is None
    connection.close()
