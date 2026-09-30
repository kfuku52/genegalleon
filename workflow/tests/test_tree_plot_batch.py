"""Real R worker isolation and equivalent PDF drawing commands."""
import json
import os
import re
import subprocess
from pathlib import Path

SUPPORT = Path(__file__).resolve().parents[1] / 'support'


def normalized(path):
    return re.sub(rb'/(CreationDate|ModDate) \([^)]*\)', b'', path.read_bytes())


def test_batch_render_matches_cli_and_recovers_from_one_family_error(tmp_path):
    (tmp_path/'identity').write_text('a file named like a scalar option')
    (tmp_path/'tree').mkdir()
    stat = tmp_path / 'stat.tsv'
    stat.write_text('branch_id\tparent\tsister\tchild1\tchild2\tnode_name\tbl_rooted\tso_event\tso_event_parent\n'
        '4\t-999\t-999\t2\t3\troot\t0\tS\tS\n2\t4\t3\t0\t1\tn4\t1\tS\tS\n'
        '0\t2\t1\t\t\tg1\t1\tL\tS\n1\t2\t0\t\t\tg2\t1\tL\tS\n3\t4\t2\t\t\tg3\t2\tL\tS\n')
    args = ['--stat_branch='+str(stat), '--max_delta_intron_present=-0.5', '--panel_widths_mm=tree:60',
            '--panel1=tree,bl_rooted,no,no,L', '--show_branch_id=no', '--event_method=species_overlap',
            '--species_color_table=PLACEHOLDER', '--pie_chart_value_transformation=identity', '--long_branch_display=no']
    baseline = subprocess.run(['Rscript', str(SUPPORT/'stat_branch2tree_plot.r'), *args], cwd=tmp_path,
        env={**os.environ, 'GG_TREE_PLOT_CHECK_GGIMAGE':'1'}, capture_output=True)
    assert baseline.returncode == 0, baseline.stderr.decode()
    expected = normalized(tmp_path/'stat_branch2tree_plot.pdf')
    domain = tmp_path/'domain.tsv'
    domain.write_text('qacc\tsacc\tstitle\tqlen\tslen\tqstart\tqend\n'
                      'g1\tPF0001\tPF0001,Domain_A\t30\t300\t1\t10\n')
    fasta = tmp_path/'alignment.fa'
    fasta.write_text(''.join(f'>g{i}\n' + 'ACGT'*21 + 'ACGTAA\n' for i in range(1,4)))
    rich = [*args,'--panel2=domain,'+str(domain),'--panel3=alignment,'+str(fasta)+','+str(fasta)]
    rendered = subprocess.run(['Rscript',str(SUPPORT/'stat_branch2tree_plot.r'),*rich],cwd=tmp_path,capture_output=True)
    assert rendered.returncode == 0, rendered.stderr.decode()
    expected_rich = normalized(tmp_path/'stat_branch2tree_plot.pdf')
    stale = tmp_path/'bad.pdf'
    stale.write_bytes(b'preserve existing failed output')
    jobs = [{'id':str(i), 'cwd':str(tmp_path), 'args':args, 'output':str(tmp_path/f'{i}.pdf')} for i in (0,2)]
    jobs[0]['args'] = rich
    jobs.insert(1, {'id':'bad', 'cwd':str(tmp_path), 'args':['--stat_branch='+str(tmp_path/'missing.tsv')],
                    'output':str(stale)})
    plan = tmp_path/'plan.json'
    plan.write_text(json.dumps(jobs))
    result = subprocess.run(['Rscript', str(SUPPORT/'tree_plot_batch.r'), str(plan)], capture_output=True)
    assert result.returncode == 1, result.stderr.decode()
    receipts = json.loads(Path(str(plan)+'.results.json').read_text())
    assert receipts['completion_evidence'] is False
    assert [(row['id'],row['exit_code']) for row in receipts['results']] == [('0',0),('bad',1),('2',0)]
    assert stale.read_bytes() == b'preserve existing failed output'
    assert normalized(tmp_path/'0.pdf') == expected_rich
    assert normalized(tmp_path/'2.pdf') == expected
    assert not list(tmp_path.glob('gg-plot-*'))



def test_batch_rejects_renderer_source_change_before_pdf_publication(tmp_path):
    import shutil
    copied = tmp_path/'support'
    copied.mkdir()
    for name in ('stat_branch2tree_plot.r','tree_plot_batch.r'):
        shutil.copy2(SUPPORT/name, copied/name)
    stat = tmp_path/'stat.tsv'
    stat.write_text('branch_id\tparent\tsister\tchild1\tchild2\tnode_name\tbl_rooted\tso_event\tso_event_parent\n'
        '4\t-999\t-999\t2\t3\troot\t0\tS\tS\n2\t4\t3\t0\t1\tn4\t1\tS\tS\n'
        '0\t2\t1\t\t\tg1\t1\tL\tS\n1\t2\t0\t\t\tg2\t1\tL\tS\n3\t4\t2\t\t\tg3\t2\tL\tS\n')
    renderer = copied/'stat_branch2tree_plot.r'
    profile = tmp_path/'profile.R'
    profile.write_text("setHook(packageEvent('cowplot','onLoad'), function(...) {\n"
        " ns=asNamespace('cowplot'); original=get('save_plot',ns); unlockBinding('save_plot',ns)\n"
        " assign('save_plot',function(...) { result=original(...); cat('# changed\\n',file="+
        json.dumps(str(renderer))+",append=TRUE); result },ns); lockBinding('save_plot',ns)\n})\n")
    output = tmp_path/'preserve.pdf'
    output.write_bytes(b'existing PDF')
    plan = tmp_path/'mutation.json'
    plan.write_text(json.dumps([{'id':'family','cwd':str(tmp_path),'output':str(output),
        'args':['--stat_branch='+str(stat),'--max_delta_intron_present=-0.5',
                '--panel_widths_mm=tree:60','--panel1=tree,bl_rooted,no,no,L',
                '--show_branch_id=no','--event_method=species_overlap',
                '--species_color_table=PLACEHOLDER','--pie_chart_value_transformation=identity',
                '--long_branch_display=no']}]))
    result = subprocess.run(['Rscript',str(copied/'tree_plot_batch.r'),str(plan)],capture_output=True,
                            env={**os.environ,'R_PROFILE_USER':str(profile)})
    assert result.returncode == 1, result.stderr.decode()
    receipts = json.loads(Path(str(plan)+'.results.json').read_text())
    assert receipts['results'][0]['exit_code'] == 1
    assert 'Renderer source changed before publication' in receipts['results'][0]['detail']
    assert output.read_bytes() == b'existing PDF'
