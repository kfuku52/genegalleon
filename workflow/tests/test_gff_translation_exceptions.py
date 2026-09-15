import pandas as pd
import pytest
import subprocess
import sys
from pathlib import Path

from workflow.support.gff2genestat import extract_by_ids, attach_transcript_structure, summarize_gene_features
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
