from pathlib import Path
import pytest
from workflow.support.repair_coge_transcripts import repair


def inputs(tmp_path):
    gff=tmp_path/'input.gff';cds=tmp_path/'cds.fa';genome=tmp_path/'genome.fa'
    gff.write_text('chr\tCoGe\tgene\t1\t12\t.\t+\t.\tID=g\n'+''.join(
        f'chr\tCoGe\t{feature}\t{start}\t{end}\t.\t+\t.\tID={identifier};Parent={parent};coge_fid=1\n'
        for feature,start,end,identifier,parent in [
            ('mRNA',1,3,'t1','g'),('exon',1,3,'e1','t1'),('CDS',1,3,'c1','t1'),
            ('mRNA',10,12,'t2','t1'),('exon',10,12,'e2','t2'),('CDS',10,12,'c2','t2')]))
    cds.write_text('>Species_a_g\nATGCCC\n');genome.write_text('>chr\nATGAAAAAACCC\n')
    return dict(gff=gff,cds=cds,genome=genome,species='Species_a',output=tmp_path/'output.gff',audit=tmp_path/'audit.json')


def test_repairs_chain_only_after_sequence_verification(tmp_path):
    args=inputs(tmp_path);before=args['gff'].read_bytes()
    result=repair(**args)
    assert len(result['repairs'])==1
    text=args['output'].read_text()
    assert text.count('\tmRNA\t')==1 and '\tmRNA\t1\t12\t' in text
    assert 'ID=c2;Parent=t1;' in text and 'ID=e2;Parent=t1;' in text
    assert args['gff'].read_bytes()==before


@pytest.mark.parametrize('defect',['sequence','length','fid','branch','all_unknown'])
def test_rejects_unverified_or_ambiguous_repairs(tmp_path,defect):
    args=inputs(tmp_path)
    if defect in {'sequence','length','all_unknown'}:
        args['cds'].write_text('>Species_a_g\n'+{'sequence':'ATGCCA','length':'ATG','all_unknown':'NNNNNN'}[defect]+'\n')
    else:
        text=args['gff'].read_text()
        if defect=='fid':text=text.replace('ID=c2;Parent=t2;coge_fid=1','ID=c2;Parent=t2;coge_fid=2')
        else:text+='chr\tCoGe\tmRNA\t4\t6\t.\t+\t.\tID=t3;Parent=t1;coge_fid=1\n'
        args['gff'].write_text(text)
    with pytest.raises(ValueError):repair(**args)
    assert not args['output'].exists() and not args['audit'].exists()


def test_minus_strand_and_masked_bases_keep_unknown_phase(tmp_path):
    args=inputs(tmp_path)
    args['gff'].write_text(args['gff'].read_text().replace('\t+\t','\t-\t'))
    args['cds'].write_text('>Species_a_g\nGGGCAN\n')
    result=repair(**args)
    assert result['repairs'][0]['masked_bases']==1
    assert all(line.split('\t')[7]=='.' for line in args['output'].read_text().splitlines() if not line.startswith('#'))
