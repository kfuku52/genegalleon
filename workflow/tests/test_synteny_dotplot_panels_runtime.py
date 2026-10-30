import copy
import csv
import json
from pathlib import Path

import pytest

from workflow.support import pairwise_synteny_panels as panels
from workflow.tests.test_pairwise_synteny_runtime import run_core, write_genome


@pytest.fixture
def saved_dotplots(tmp_path):
    plots, ds = tmp_path / 'plots', tmp_path / 'ds'
    species = ['Anc_species', 'Dio_species', 'Tri_species']
    cds, beds = {}, {}
    for sp, prefix in zip(species, 'adt', strict=True):
        cds[sp] = tmp_path / (sp + '.fa')
        cds[sp].write_text(''.join(f'>{prefix}{i}\nATGGCTGCT\n' for i in range(8)))
        beds[sp] = ''.join(f'chr{1 + i // 4}\t{i * 10}\t{i * 10 + 5}\t{prefix}{i}\n' for i in range(8))
    pairs = []
    for comparison, a, b, pa, pb in [('ad', 0, 1, 'a', 'd'), ('ta', 2, 0, 't', 'a'), ('td', 2, 1, 't', 'd')]:
        pair = {'analysis_id': comparison, 'target': {'species': species[a]}, 'query': {'species': species[b]},
                'target_species': species[a], 'query_species': species[b],
                'ds': {side: {'fasta': str(cds[species[index]])} for side, index in [('target', a), ('query', b)]}}
        pairs.append(pair)
        directory, distances = plots / comparison, ds / comparison
        directory.mkdir(parents=True)
        distances.mkdir(parents=True)
        for side, index in [('target', a), ('query', b)]:
            (directory / (side + '.bed')).write_text(beds[species[index]])
            (directory / ('dotplot.' + side + '.bed')).write_text(beds[species[index]])
        # Two blocks overlap at the first pair. Its dS is estimated once, and
        # the combined figure must retain one point for every unique pair.
        (directory / 'target.query.lifted.anchors').write_text(
            '###\n' + ''.join(f'{pa}{i}\t{pb}{i}\t100\n' for i in range(8)) + f'###\n{pa}0\t{pb}0\t100\n')
        with (distances / 'ds.tsv').open('w', newline='') as stream:
            writer = csv.writer(stream, delimiter='\t')
            writer.writerow(['pair_id', 'dS', 'status'])
            for i in range(8):
                writer.writerow([f'{pa}{i}|{pb}{i}', '' if i == 0 else (3 if i == 1 else .5),
                                 'saturated' if i == 0 else 'ok'])
    config = {'panel_id': 'three', 'comparisons': ['ad', 'ta', 'td'], 'columns': 2,
              'tracks': [{'species': s, 'seqids': ['chr2', 'chr1'],
                          'labels': {'chr1': '1', 'chr2': '2'}} for s in species]}
    config['tracks'][0]['orientation'] = {'chr2': '-'}
    plan = {'pairs': pairs, 'formats': ['pdf', 'svg', 'png'], 'dotplot_color': 'ds', 'ds_color_max': 2,
            'dotplot_panels': [config]}
    return plan, plots, ds


def test_three_pairwise_panels_keep_missing_high_ds_and_common_layout(saved_dotplots, monkeypatch):
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.axes import Axes
    from PIL import Image
    from pypdf import PdfReader

    plan, plots, ds = saved_dotplots
    before = {p: p.read_bytes() for p in plots.rglob('*') if p.is_file()}
    observed, real = [], Axes.scatter
    def capture(ax, x, y, **kwargs):
        observed.append((list(x), list(y), np.asarray(kwargs['c']), kwargs['norm']))
        return real(ax, x, y, **kwargs)
    monkeypatch.setattr(Axes, 'scatter', capture)
    figures = plt.get_fignums()
    panels.render_workflow_panels(plan, plots, ds)
    output = plots / 'dotplot_panels/three'
    report = json.loads((output / 'dotplot_panels.json').read_text())
    assert report['species_order'] == ['Anc_species', 'Dio_species', 'Tri_species']
    assert report['shared_ds_scale'] == [0, 2]
    assert report['tracks'][0]['orientation'] == {'chr2': '-'}
    assert report['tracks'][1]['seqids'] == ['chr2', 'chr1']
    assert [r['axes_swapped'] for r in report['comparisons']] == [False, True, True]
    for row, (x, y, values, norm) in zip(report['comparisons'], observed, strict=True):
        assert row['displayed_pairs'] == row['unique_anchor_pairs'] == 8
        assert row['missing_ds_pairs'] == row['above_color_limit'] == 1
        assert row['status_counts'] == {'ok': 7, 'saturated': 1}
        assert len(x) == 8 and np.isnan(values).sum() == 1 and max(values[1:]) == 3
        assert (norm.vmin, norm.vmax, norm.clip) == (0, 2, True)
    # chr2 is first and reversed for Anc in both comparisons. The same gene
    # must keep its exact display rank even when Anc is the analysis query.
    assert observed[0][0][4:8] == [3.5, 2.5, 1.5, .5]
    assert observed[1][0] == observed[0][0]
    assert {p: p.read_bytes() for p in before} == before
    assert plt.get_fignums() == figures
    page = PdfReader(output / 'dotplot.pdf').pages[0]
    assert float(page.mediabox.width) == pytest.approx(518.4)
    text = page.extract_text()
    assert text.count('dS') == 1 and text.count('Unestimable or saturated') == 1
    assert text.count('Anc species') == text.count('Dio species') == text.count('Tri species') == 2
    assert '(rev)' in text
    for suffix in ['png', 'svg']:
        assert (output / ('dotplot.' + suffix)).stat().st_size > 1000
    with Image.open(output / 'dotplot.png') as image:
        image.verify()
    assert len(list((output / 'inputs').glob('*/*.bed'))) == 12


@pytest.mark.parametrize('change,message', [
    ('bed', 'different BED'), ('cds', 'different BED'), ('eligible', 'different BED'),
    ('missing-ds', 'exactly every'), ('extra-ds', 'exactly every'),
    ('unknown-seqid', 'physical-length filter'), ('empty-pairs', 'No anchors'),
])
def test_inconsistent_or_incomplete_inputs_do_not_publish(saved_dotplots, change, message):
    plan, plots, ds = saved_dotplots
    if change == 'bed':
        file = plots / 'ta/query.bed'
        file.write_text(file.read_text().replace('\t0\t', '\t1\t', 1))
    elif change == 'cds':
        source = Path(plan['pairs'][1]['ds']['query']['fasta'])
        different = source.with_suffix('.different.fa')
        different.write_text(source.read_text().replace('GCT', 'GCC', 1))
        plan['pairs'][1]['ds']['query']['fasta'] = str(different)
    elif change == 'eligible':
        file = plots / 'ta/dotplot.query.bed'
        file.write_text(''.join(line for line in file.read_text().splitlines(True) if not line.startswith('chr2')))
    elif change == 'missing-ds':
        file = ds / 'ad/ds.tsv'
        file.write_text(''.join(file.read_text().splitlines(True)[:-1]))
    elif change == 'extra-ds':
        with (ds / 'ad/ds.tsv').open('a') as f:
            f.write('other|pair\t1\tok\n')
    elif change == 'unknown-seqid':
        plan['dotplot_panels'][0]['tracks'][0]['seqids'] = ['absent']
    else:
        for track in plan['dotplot_panels'][0]['tracks']:
            track['seqids'] = ['chr1'] if track['species'] == 'Anc_species' else ['chr2']
            track.pop('labels', None)
            track.pop('orientation', None)
    with pytest.raises(ValueError, match=message):
        panels.render_workflow_panels(plan, plots, ds)
    assert not (plots / 'dotplot_panels/three').exists()


def test_source_replacement_during_save_never_publishes(saved_dotplots, monkeypatch):
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure

    plan, plots, ds = saved_dotplots
    real = Figure.savefig
    def replaced(fig, path, **kwargs):
        real(fig, path, **kwargs)
        source = plots / 'ad/target.bed'
        replacement = source.with_suffix('.new')
        replacement.write_bytes(source.read_bytes())
        replacement.replace(source)
    monkeypatch.setattr(Figure, 'savefig', replaced)
    previous = plt.get_fignums()
    with pytest.raises(OSError, match='changed'):
        panels.render_workflow_panels(plan, plots, ds)
    assert plt.get_fignums() == previous
    assert not (plots / 'dotplot_panels/three').exists()


def test_existing_output_is_preserved(saved_dotplots):
    plan, plots, ds = saved_dotplots
    target = plots / 'dotplot_panels/three'
    target.mkdir(parents=True)
    (target / 'keep.txt').write_text('keep')
    with pytest.raises(ValueError, match='must be new'):
        panels.render_workflow_panels(plan, plots, ds)
    assert (target / 'keep.txt').read_text() == 'keep'


@pytest.mark.parametrize('change', ['unknown-comparison', 'duplicate-comparison', 'missing-species',
                                  'extra-track-option', 'duplicate-seqid', 'invalid-orientation',
                                  'duplicate-label', 'boolean-columns', 'duplicate-panel', 'unknown-option'])
def test_panel_configuration_boundary_rejects_invalid_requests(saved_dotplots, tmp_path, change):
    plan, _plots, _ds = saved_dotplots
    config = {'panels': copy.deepcopy(plan['dotplot_panels'])}
    panel = config['panels'][0]
    if change == 'unknown-comparison':
        panel['comparisons'].append('missing')
    elif change == 'duplicate-comparison':
        panel['comparisons'].append('ad')
    elif change == 'missing-species':
        panel['tracks'].pop()
    elif change == 'extra-track-option':
        panel['tracks'][0]['bed'] = '/unexpected/source'
    elif change == 'duplicate-seqid':
        panel['tracks'][0]['seqids'].append('chr1')
    elif change == 'invalid-orientation':
        panel['tracks'][0]['orientation']['chr1'] = 'reverse'
    elif change == 'duplicate-label':
        panel['tracks'][0]['labels']['chr2'] = '1'
    elif change == 'boolean-columns':
        panel['columns'] = True
    elif change == 'duplicate-panel':
        config['panels'].append(panel)
    else:
        panel['unknown'] = 1
    file = tmp_path / 'config.json'
    file.write_text(json.dumps(config))
    with pytest.raises(ValueError):
        panels.load_panels(file, plan['pairs'])


def test_panel_config_changes_only_plot_contract(tmp_path):
    from workflow.support.pairwise_synteny import contract_args

    plan = {'workspace': str(tmp_path), 'pairs': [], 'parameters': {}, 'tools': {'jcvi': 'test'},
            'formats': ['svg'], 'ds_tools': {'source_hashes': {}},
            'dotplot_panels_file': str(tmp_path / 'panels.json'),
            'dotplot_panels': [{'panel_id': 'three', 'columns': 2}]}
    for phase in ['analysis', 'ds', 'plots']:
        before = contract_args(plan, phase)
        plan['dotplot_panels'][0]['columns'] = 1
        after = contract_args(plan, phase)
        assert (before == after) == (phase != 'plots')
        plan['dotplot_panels'][0]['columns'] = 2


def test_native_three_pair_workflow_and_plot_only_reuse(tmp_path):
    workspace = tmp_path / 'workspace'
    for species, reverse in [('Anc_species', False), ('Dio_species', True), ('Tri_species', False)]:
        for mode in ['protein', 'cds']:
            write_genome(workspace, species, mode, [('Chr1', reverse)])
    (workspace / 'input/synteny_pairs.tsv').write_text(
        'analysis_id\ttarget_species\tquery_species\nad\tAnc_species\tDio_species\n'
        'ta\tTri_species\tAnc_species\ntd\tTri_species\tDio_species\n')
    file = workspace / 'input/panels.json'
    config = {'panels': [{'panel_id': 'three', 'comparisons': ['ad', 'ta', 'td'],
                          'tracks': [{'species': s} for s in ['Anc_species', 'Dio_species', 'Tri_species']]}]}
    file.write_text(json.dumps(config))
    result = run_core(workspace, synteny_dotplot_color='ds', synteny_plot_formats='svg',
                      synteny_karyotype_sort='none', synteny_dotplot_panels_file=str(file))
    assert result.returncode == 0, result.stdout + result.stderr
    root = workspace / 'output/genome_evolution/synteny'
    immutable = {p: (p.read_bytes(), p.stat().st_mtime_ns)
                 for phase in ['analysis', 'ds'] for p in (root / phase).rglob('*') if p.is_file()}
    output = root / 'plots/dotplot_panels/three'
    assert len(json.loads((output / 'dotplot_panels.json').read_text())['comparisons']) == 3
    config['panels'][0]['columns'] = 1
    config['panels'][0]['tracks'][1]['orientation'] = {'Chr1': '-'}
    file.write_text(json.dumps(config))
    rerun = run_core(workspace, synteny_dotplot_color='ds', synteny_plot_formats='svg',
                     synteny_karyotype_sort='none', synteny_dotplot_panels_file=str(file),
                     synteny_plot_only='1', artifact_stale_policy='rebuild', synteny_ds_color_max='1.5')
    assert rerun.returncode == 0, rerun.stdout + rerun.stderr
    assert {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in immutable} == immutable
    report = json.loads((output / 'dotplot_panels.json').read_text())
    assert report['shared_ds_scale'] == [0, 1.5]
    assert report['tracks'][1]['orientation'] == {'Chr1': '-'}
    assert json.loads((root / 'plots/ad/dotplot_order.json').read_text())['orientation_changed'] is False
    mtime = (output / 'dotplot.svg').stat().st_mtime_ns
    unchanged = run_core(workspace, synteny_dotplot_color='ds', synteny_plot_formats='svg',
                         synteny_karyotype_sort='none', synteny_dotplot_panels_file=str(file),
                         synteny_plot_only='1', synteny_ds_color_max='1.5')
    assert unchanged.returncode == 0, unchanged.stdout + unchanged.stderr
    assert (output / 'dotplot.svg').stat().st_mtime_ns == mtime


def test_combined_panels_require_explicit_ds_mode_before_analysis(tmp_path):
    workspace = tmp_path / 'workspace'
    for species in ['Anc_species', 'Dio_species']:
        write_genome(workspace, species, 'protein', [('Chr1', False)])
    (workspace / 'input/synteny_pairs.tsv').write_text(
        'analysis_id\ttarget_species\tquery_species\nad\tAnc_species\tDio_species\n')
    file = workspace / 'input/panels.json'
    file.write_text(json.dumps({'panels': [{'panel_id': 'two', 'comparisons': ['ad'],
        'tracks': [{'species': 'Anc_species'}, {'species': 'Dio_species'}]}]}))
    result = run_core(workspace, synteny_dotplot_panels_file=str(file))
    assert result.returncode != 0 and 'require dotplot-color ds' in result.stderr
    assert not (workspace / 'output/genome_evolution/synteny/analysis').exists()


def test_defaults_preserve_the_first_native_filtered_bed_order(saved_dotplots):
    plan, plots, ds = saved_dotplots
    for track in plan['dotplot_panels'][0]['tracks']:
        track.pop('seqids')
        track.pop('orientation', None)
    for file in plots.glob('*/dotplot.*.bed'):
        rows = file.read_text().splitlines(True)
        file.write_text(''.join(rows[4:] + rows[:4]))
    panels.render_workflow_panels(plan, plots, ds)
    report = json.loads((plots / 'dotplot_panels/three/dotplot_panels.json').read_text())
    assert all(t['seqids'] == ['chr2', 'chr1'] for t in report['tracks'])


def test_equal_coordinate_genes_keep_the_native_bed_rank(tmp_path):
    path = tmp_path / 'ties.bed'
    path.write_text('chr2\t0\t5\tg2\nchr1\t10\t15\tg10\nchr1\t10\t15\tg1\n')
    genes = panels.ordered_bed(path)
    assert [g.gene_id for g in genes] == ['g2', 'g10', 'g1']
    ranks, bounds, _labels = panels.gene_ranks(genes, {'seqids': ['chr1', 'chr2']})
    assert ranks == {'g10': .5, 'g1': 1.5, 'g2': 2.5} and bounds == [0, 2, 3]
