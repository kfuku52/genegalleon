import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from workflow.support.gff2genestat import attach_transcript_structure, extract_by_ids, summarize_gene_features
from workflow.support.gff_feature_structure import ordered_annotated_blocks


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('exception,mode', [('exception=ribosomal slippage','ribosomal-slippage'), ('pseudo=true','pseudogene')])
def test_slippage_keeps_overlap_and_scaffold_but_not_global_phase(strand, exception, mode):
    attrs = 'ID=cds1;Parent=g1;'+exception
    cds = pd.DataFrame([
        ['s1', strand, 1, 10, attrs, 0],
        ['s1', strand, 10, 20, attrs, 2],
    ], columns=['sequence', 'strand', 'start', 'end', 'attributes', 'phase'])
    cds = cds.assign(gene_id='Species_a_g1', selected_transcript='cds1', feature='CDS')
    annotated = attach_transcript_structure(cds, cds)
    result = summarize_gene_features(annotated, ['gene_id', 'feature_size', 'chromosome', 'phase_status', 'cds_first_phase', 'num_intron'])
    assert result.iloc[0].feature_size == 21  # not the genomic union (20)
    assert result.iloc[0].chromosome == 's1'
    assert result.iloc[0].phase_status == mode
    assert pd.isna(result.iloc[0].cds_first_phase)
    assert pd.isna(result.iloc[0].num_intron)


@pytest.mark.parametrize('strand', ['+', '-'])
def test_source_audited_overlap_keeps_reused_bases_without_biological_exception(strand):
    attrs = 'ID=cds1;Parent=g1;gg_source_overlap=confirmed'
    cds = pd.DataFrame([
        ['s1', strand, 1, 10, attrs, 0],
        ['s1', strand, 10, 20, attrs, 0],
    ], columns=['sequence', 'strand', 'start', 'end', 'attributes', 'phase'])
    cds = cds.assign(gene_id='Species_a_g1', selected_transcript='cds1', feature='CDS')
    annotated = attach_transcript_structure(cds, cds)
    result = summarize_gene_features(annotated, ['gene_id', 'feature_size', 'chromosome', 'phase_status', 'cds_first_phase', 'num_intron'])
    assert result.iloc[0].feature_size == 21
    assert result.iloc[0].phase_status == 'source-overlap'
    assert pd.isna(result.iloc[0].cds_first_phase)
    assert pd.isna(result.iloc[0].num_intron)


def test_source_audited_overlap_requires_marker_on_every_part():
    with pytest.raises(ValueError, match='source-overlap'):
        ordered_annotated_blocks([
            ('s', '+', 1, 10, 'ID=a;Parent=g;gg_source_overlap=confirmed'),
            ('s', '+', 10, 20, 'ID=a;Parent=g'),
        ], 'g')


def test_source_audited_marker_cannot_bypass_overlap_validation():
    with pytest.raises(ValueError, match='requires overlapping blocks'):
        ordered_annotated_blocks([
            ('s', '+', 1, 10, 'ID=a;Parent=g;gg_source_overlap=confirmed'),
            ('s', '+', 20, 30, 'ID=a;Parent=g;gg_source_overlap=confirmed'),
        ], 'g')


def test_partial_cds_length_validation_uses_fuzzy_termini_and_phase():
    from workflow.support.gff2genestat import validate_cds_lengths
    traits = pd.DataFrame([
        dict(gene_id='start', feature_size=242, cds_first_phase=2,
             splice_mode='cis', cds_partial='5prime'),
        dict(gene_id='end', feature_size=193, cds_first_phase=0,
             splice_mode='cis', cds_partial='3prime'),
        dict(gene_id='pseudo', feature_size=245, cds_first_phase=float('nan'),
             splice_mode='pseudogene', cds_partial='5prime'),
    ])
    records = [
        ('start', 'start', 'A' * 240),
        ('end', 'end', 'A' * 195),
        ('pseudo', 'pseudo', 'A' * 246),
    ]
    validate_cds_lengths(traits, records)


def test_partial_cds_length_validation_does_not_relax_complete_records():
    from workflow.support.gff2genestat import validate_cds_lengths
    traits = pd.DataFrame([dict(gene_id='complete', feature_size=10, cds_partial='none')])
    with pytest.raises(ValueError, match='GFF=10, CDS=11'):
        validate_cds_lengths(traits, [('complete', 'complete', 'A' * 11)])


@pytest.mark.parametrize('attrs', [
    ('ID=a;Parent=g', 'ID=a;Parent=g'),
    ('ID=a;Parent=g;exception=ribosomal slippage', 'ID=b;Parent=g;exception=ribosomal slippage'),
    ('ID=a;Parent=g;exception=ribosomal slippage', 'ID=a;Parent=g'),
])
def test_unexplained_or_mixed_overlap_still_fails(attrs):
    with pytest.raises(ValueError):
        ordered_annotated_blocks([('s', '+', 1, 10, attrs[0]), ('s', '+', 10, 20, attrs[1])], 'g')


def test_direct_gene_parent_does_not_merge_alternative_proteins():
    gff = pd.DataFrame([
        ['s', '+', 1, 90, 'ID=g1', 'gene'],
        ['s', '+', 1, 90, 'ID=long;Parent=g1', 'CDS'],
        ['s', '+', 10, 90, 'ID=short;Parent=g1', 'CDS'],
    ], columns=['sequence', 'strand', 'start', 'end', 'attributes', 'feature'])
    result = extract_by_ids(gff, pd.Series(['Species_a_g1']), 'CDS', 'longest')
    assert result.selected_transcript.tolist() == ['long']
    assert result.start.tolist() == [1]


def test_phase_report_preserves_coordinates_and_strict_still_fails():
    cds = pd.DataFrame([
        ['s', '+', 1, 10, 'ID=c;Parent=t', 0],
        ['s', '+', 20, 30, 'ID=c;Parent=t', 0],
    ], columns=['sequence', 'strand', 'start', 'end', 'attributes', 'phase'])
    cds = cds.assign(gene_id='g', selected_transcript='t', feature='CDS')
    with pytest.raises(ValueError, match='Conflicting CDS phases'):
        attach_transcript_structure(cds, cds)
    result = attach_transcript_structure(cds, cds, phase_policy='report')
    assert result.phase_status.eq('conflicting').all()
    assert result.cds_first_phase.isna().all()
    assert result.start.tolist() == [1, 20]
    cds.loc[1, 'strand'] = '-'
    with pytest.raises(ValueError, match='coordinate systems'):
        attach_transcript_structure(cds, cds, phase_policy='report')


def test_require_matches_does_not_publish_empty_table(tmp_path):
    (tmp_path/'Species_a.gff').write_text('s\tx\tCDS\t1\t9\t.\t+\t0\tID=wrong\n')
    fasta = tmp_path/'cds.fa'
    fasta.write_text('>Species_a_right\nATGAAATAA\n')
    out = tmp_path/'out.tsv'
    script = Path(__file__).resolve().parents[1]/'support/gff2genestat.py'
    result = subprocess.run([sys.executable, str(script), '--dir_gff', str(tmp_path),
                             '--seqfile', str(fasta), '--outfile', str(out), '--require-matches'],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert 'No GFF IDs matched' in result.stderr
    assert not out.exists()


def test_structural_gene_id_precedes_another_genes_display_alias():
    gff = pd.DataFrame([
        ['s', '+', 1, 90, 'ID=gene:z;Name=a;gene_id=z', 'gene'],
        ['s', '+', 1, 90, 'ID=t;Parent=gene:z', 'mRNA'],
        ['s', '+', 1, 90, 'ID=c;Parent=t', 'CDS'],
    ], columns=['sequence', 'strand', 'start', 'end', 'attributes', 'feature'])
    result=extract_by_ids(gff,pd.Series(['Species_a_a','Species_a_z']),'CDS','longest')
    assert result.gene_id.tolist()==['Species_a_z']


def test_numbered_cross_contig_cds_preserves_source_order():
    rows = [('b', '+', 566, 1135, 'ID=c3;Parent=t;number=3'),
            ('a', '+', 248256, 248359, 'ID=c1;Parent=t;number=1'),
            ('a', '+', 248577, 249072, 'ID=c2;Parent=t;number=2')]
    blocks, mode = ordered_annotated_blocks(rows, 'g')
    assert mode == 'ordered-fragments'
    assert blocks == [row[:4] for row in rows[1:] + rows[:1]]
    assert sum(b[3] - b[2] + 1 for b in blocks) == 1170
    cds = pd.DataFrame([(*row, phase) for row, phase in zip(rows, [0, 0, 2], strict=True)],
                       columns=['sequence', 'strand', 'start', 'end', 'attributes', 'phase'])
    cds = cds.assign(gene_id='Species_a_g', selected_transcript='t', feature='CDS')
    with pytest.raises(ValueError, match='Conflicting CDS phases'):
        attach_transcript_structure(cds, cds)
    annotated = attach_transcript_structure(cds, cds, phase_policy='report')
    result = summarize_gene_features(annotated, ['gene_id', 'feature_size', 'chromosome', 'num_intron',
                                                'splice_mode', 'feature_blocks', 'transcript_junction_positions'])
    record = result.iloc[0]
    assert record.feature_size == 1170 and record.chromosome == ''
    assert pd.isna(record.num_intron) and record.splice_mode == 'ordered-fragments'
    assert record.transcript_junction_positions == '104;600'
    for damaged in [
        [(*r[:4], r[4].replace(';number=3', '')) for r in rows],
        [(*r[:4], r[4].replace('number=3', 'number=4')) for r in rows],
        [(*r[:4], r[4].replace('number=3', 'number=2')) for r in rows],
        [(*r[:4], r[4].replace('Parent=t;number=3', 'Parent=other;number=3')) for r in rows],
        [rows[0], (*rows[1][:4], 'ID=c1;Parent=t;number=2'), (*rows[2][:4], 'ID=c2;Parent=t;number=1')],
    ]:
        with pytest.raises(ValueError):
            ordered_annotated_blocks(damaged, 'g')
