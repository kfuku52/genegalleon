import json
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))

from focus_hgt_gene_trees import annotate, background_supported  # noqa: E402
from focus_hgt_traits import generate as generate_filtered  # noqa: E402
from focus_hgt_traits import read_tsv, write_tsv  # noqa: E402
from species_trait_schema import schema_path, schema_payload  # noqa: E402


def generate(*args, **kwargs):
    # Existing trait/context tests isolate that behavior from the additional
    # domain filter. Its enabled/default path has dedicated regressions below.
    kwargs.setdefault('require_shared_pfam', False)
    kwargs.setdefault('require_length_ratio', False)
    return generate_filtered(*args, **kwargs)


def saved_species_taxonomy(path):
    rows = [dict(species=name, resolution_status='resolved', domain=domain, phylum=phylum, **{'class': klass})
            for name, domain, phylum, klass in [('A', 'Eukaryota', 'Arthropoda', 'Insecta'),
                                               ('B', 'Eukaryota', 'Arthropoda', 'Insecta'),
                                               ('C', 'Eukaryota', 'Arthropoda', 'Branchiopoda'),
                                               ('D', 'Eukaryota', 'Streptophyta', 'Bryopsida'),
                                               ('K', 'Bacteria', '', ''),
                                               ('U', 'Eukaryota', '', '')]]
    write_tsv(path, list(rows[0]), rows)
    return path


@pytest.mark.parametrize('donor,recipient,status,reason', [
    ('D', 'A', 'passed', 'all_donor_tips_outside_arthropoda_all_recipient_tips_in_insecta'),
    ('K', 'all_insects', 'passed', 'all_donor_tips_outside_arthropoda_all_recipient_tips_in_insecta'),
    ('C', 'A', 'excluded_direction', 'donor_within_arthropoda'),
    ('mixed_donor', 'A', 'withheld', 'donor_mixed_arthropoda'),
    ('unknown_donor', 'A', 'withheld', 'donor_unknown_arthropoda'),
    ('missing', 'A', 'withheld', 'donor_unmapped_arthropoda'),
    ('D', 'C', 'excluded_direction', 'recipient_outside_insecta'),
    ('D', 'mixed_recipient', 'withheld', 'recipient_mixed_insecta'),
    ('D', 'U', 'withheld', 'recipient_unknown_insecta'),
])
def test_direction_uses_exact_species_branches_and_withholds_unknown_tips(tmp_path, donor, recipient, status, reason):
    from focus_hgt_direction import filter_events
    path = saved_species_taxonomy(tmp_path/'taxonomy.tsv')
    nodes = {name: (name,) for name in ('A', 'B', 'C', 'D', 'K', 'U')}
    nodes.update(all_insects=('A', 'B'), mixed_donor=('C', 'D'), unknown_donor=('D', 'U'), mixed_recipient=('A', 'C'))
    event = dict(event_id='OG1:3:1', orthogroup='OG1', generax_donor_node=donor, generax_recipient_node=recipient,
                 generax_transfer=f'Y@{donor}@{recipient}', donor_phylum='Arthropoda', donor_class='Insecta')
    before = path.read_bytes()
    selected, audit, branches, _ = filter_events([event], nodes, path)
    assert audit[0]['direction_filter_status'] == status and audit[0]['direction_filter_reason'] == reason
    assert bool(selected) == (status == 'passed')
    assert 'direction_filter_status' not in event and path.read_bytes() == before
    mixed = next(r for r in branches if r['species_branch'] == 'mixed_donor')
    assert mixed['clade_tip_labels'] == 'C; D' and mixed['arthropoda_within_tip_count'] == 1


def test_direction_taxonomy_aliases_are_exact_and_malformed_classifications_fail(tmp_path):
    from focus_hgt_direction import filter_events, read_taxonomy
    path = saved_species_taxonomy(tmp_path/'taxonomy.tsv')
    fields, rows = read_tsv(path)
    rows[0].update(species='canonical_A', tree_status='mapped', tree_tip='A')
    fields += ['tree_status', 'tree_tip']
    write_tsv(path, fields, rows)
    aliases, _ = read_taxonomy(path)
    assert aliases['A'] is aliases['canonical_A']
    rows[1].update(tree_status='mapped', tree_tip='A')
    write_tsv(path, fields, rows)
    with pytest.raises(ValueError, match='alias'):
        read_taxonomy(path)
    rows[1].update(tree_tip='B', phylum='Streptophyta')
    write_tsv(path, fields, rows)
    with pytest.raises(ValueError, match='Arthropoda'):
        read_taxonomy(path)
    with pytest.raises(ValueError, match='transfer token'):
        filter_events([dict(generax_donor_node='D', generax_recipient_node='A', generax_transfer='Y@C@A')],
                      {'A': ('A',), 'D': ('D',)}, saved_species_taxonomy(path))


def test_direction_is_evaluated_once_after_pfam_before_multiple_traits(source, monkeypatch):
    import focus_hgt_direction
    _, links = read_tsv(source[1])
    links = [dict(supported_link(r['event_id'], r['side'], r['gene_id']), **r) for r in links]
    write_tsv(source[1], list(links[0]), links)
    root = source[0].parent/'families'
    saved_pfam(root, dict(D_gene=['PF01053'], A_gene=['PF01053'], B_gene=[], C_gene=['PF01053']))
    taxonomy = saved_species_taxonomy(source[0].parent/'taxonomy.tsv')
    calls = []
    original = focus_hgt_direction.filter_events
    def capture(events, *args):
        calls.append({r['event_id'] for r in events})
        assert all(r['pfam_filter_status'] == 'passed' and 'focus_trait' not in r for r in events)
        return original(events, *args)
    monkeypatch.setattr(focus_hgt_direction, 'filter_events', capture)
    before = taxonomy.read_bytes()
    report = generate_filtered(*source, plots=False, gene_family_root=root,
                               direction_filter='non_arthropoda_to_insecta', species_taxonomy=taxonomy)
    assert calls == [{'OG1:3:1', 'OG1:3:2', 'OG1:3:3', 'OG1:3:5'}]
    assert report['shared_pfam_filter']['passed_event_count'] == 4
    assert report['shared_direction_filter']['passed_event_count'] == 2
    assert report['filtering_order'] == ['input_cohort', 'pfam_pair_filter', 'species_branch_direction_filter', 'trait_category1']
    assert len(read_tsv(source[-1]/'direction_event_audit.tsv')[1]) == 4
    assert len(read_tsv(source[-1]/'direction_species_branches.tsv')[1]) == 7
    assert {r['event_id'] for r in read_tsv(source[-1]/'traits/category/all_category1/events.tsv')[1]} == {'OG1:3:1'}
    assert report['pfam_filter']['category']['category1_passed_event_count'] == 2
    assert str(taxonomy.resolve()) in report['inputs_sha256'] and taxonomy.read_bytes() == before
    with pytest.raises(ValueError, match='requires existing species taxonomy'):
        generate_filtered(*source, plots=False, gene_family_root=root, direction_filter='non_arthropoda_to_insecta')


def test_direction_filter_flow_is_post_pfam_and_rejects_cohort_mismatches():
    from focus_hgt_figures import filtering_counts
    events = [dict(event_id=str(i), orthogroup='OG'+str(i)) for i in range(4)]
    counts = filtering_counts(events, events[:1], pfam_selected=events[:3], direction_selected=events[:2])
    assert [r['stage'] for r in counts] == ['Input supported-event cohort', 'Event-gene pair Pfam filter',
                                           'Non-Arthropoda donor & category = 1 recipient']
    assert [r['event_count'] for r in counts] == [4, 3, 1]
    with pytest.raises(ValueError, match='subset'):
        filtering_counts(events, events[:1], pfam_selected=events[:1], direction_selected=events[:2])


def test_combined_flow_keeps_verified_support_grain_and_historical_audit(tmp_path):
    from focus_hgt_figures import export_filtering_flow, filtering_counts
    events = [dict(event_id=str(i), orthogroup='OG'+str(i), donor_classification='outside',
                   recipient_classification='insect', status='accepted', support_used='90',
                   support_source='support_generax_ufboot') for i in range(2)]
    # A numeric support on an excluded event is not independently verified
    # support, so removing its direction row cannot make it a passing event.
    audit = events + [dict(events[0], event_id='other', status='excluded_direction',
                          donor_classification='insect', support_used='99', support_source='NA')]
    path = tmp_path/'audit.tsv'
    write_tsv(path, list(audit[0]), audit)
    original = path.read_bytes()
    counts = export_filtering_flow(tmp_path/'plots', events, events[:1], 'gall', path,
                                   pfam_selected=events, direction_selected=events[:1], support_filter_enabled=True)
    assert [r['event_count'] for r in counts] == [3, 2, 2, 2, 1]
    assert counts[-1]['stage'] == 'Non-Arthropoda donor & gall = 1 recipient'
    assert len(counts) == 5 and counts[2]['stage'] == 'Bilateral scaffold background'
    assert path.read_bytes() == original
    legacy = filtering_counts(events, events[:1], path, pfam_selected=events, support_filter_enabled=True)
    assert legacy[1]['stage'] == 'Non-Insecta to Insecta'
    from pypdf import PdfReader
    text = PdfReader(tmp_path/'plots/filtering_flow.pdf').pages[0].extract_text()
    assert 'Upstream directional cohort' not in text
    assert 'Non-Arthropoda donor & gall = 1 recipient' in text
    assert 'previously verified input cohort' in text


def supported_link(event_id, side, gene):
    return dict(event_id=event_id, orthogroup="OG1", side=side, gene_id=gene,
                eligible_for_context="True", host_scaffold_status="measured", host_scaffold_id="scaffold1",
                host_scaffold_background_class_total_count="20", host_scaffold_background_class_compatible_count="9",
                host_scaffold_background_class_incompatible_count="1", host_scaffold_background_class_unresolved_count="10",
                host_scaffold_background_class_classified_fraction="0.5",
                host_scaffold_background_class_compatible_fraction="0.9")


def saved_pfam(root, gene_domains, family='OG1', legacy=False):
    rows = []
    fields = ['qacc', 'sacc', 'qlen', 'stitle', 'qstart', 'qend', 'evalue']
    for gene, domains in gene_domains.items():
        for domain in domains or ['']:
            rows.append(dict(qacc=gene, sacc='123' if domain else '', qlen='100',
                             stitle='pfam' + domain.removeprefix('PF') + ', Name, Description' if domain else '',
                             qstart='1' if domain else '', qend='90' if domain else '', evalue='1e-20' if domain else ''))
    path = root/'rpsblast'/(family + ('.rpsblast.tsv' if legacy else '_rpsblast.tsv'))
    write_tsv(path, fields, rows)
    return path


@pytest.mark.parametrize('donor,recipient,allow,expected,status', [
    (['PF01053'], ['PF01053'], False, True, 'shared_pfam_detected'),
    (['PF01053'], ['PF00001'], True, False, 'detected_pfam_sets_disjoint'),
    ([], ['PF01053'], True, False, 'one_searched_no_pfam_hit'),
    (['PF01053'], [], True, False, 'one_searched_no_pfam_hit'),
    ([], [], False, False, 'both_searched_no_pfam_hit'),
    ([], [], True, True, 'both_searched_no_pfam_hit'),
    (None, [], True, False, 'annotation_record_unavailable'),
    (None, None, True, False, 'annotation_record_unavailable'),
])
def test_pfam_pair_filter_distinguishes_saved_no_hits_from_missing(tmp_path, donor, recipient, allow, expected, status):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    saved_pfam(tmp_path, {gene: domains for gene, domains in [('D_gene', donor), ('A_gene', recipient)]
                          if domains is not None})
    result, audit, pairs, genes, sources = filter_events(events, links, tmp_path, allow_both_no_pfam=allow)
    assert bool(result) == expected
    assert audit[0]['pfam_filter_status'] == ('passed' if expected else 'withheld')
    assert pairs[0]['pair_status'] == status
    assert pairs[0]['passes_pfam_filter'] == str(expected)
    assert len(genes) == 2 and sources


def test_pfam_filter_is_exact_event_pair_and_supported_gene_specific(tmp_path):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    other = dict(events[0], event_id='OG1:3:2', event_index='2', generax_transfer='Y@D@B')
    links += [supported_link(other['event_id'], 'donor', 'D_other'),
              supported_link(other['event_id'], 'recipient', 'B_gene'),
              dict(supported_link(events[0]['event_id'], 'recipient', 'A_unconfirmed'), host_scaffold_status='unavailable'),
              dict(supported_link(events[0]['event_id'], 'recipient', 'A_transferred_out'), lineage_status='transferred_out')]
    saved_pfam(tmp_path, dict(D_gene=['PF01053'], A_gene=['PF00001'], D_other=['PF00001'],
                              B_gene=['PF00001'], A_unconfirmed=['PF01053'], A_transferred_out=['PF01053'],
                              neighbor=['PF01053']))
    selected, audit, pairs, _, _ = filter_events(events + [other], links, tmp_path)
    assert [r['event_id'] for r in selected] == [other['event_id']]
    assert [r['pfam_passing_pair_count'] for r in audit] == [0, 1]
    assert len(pairs) == 2  # No event-wide/domain-wide unions, unconfirmed genes or neighbor rescue.
    # Additional recipient copy forms its own pair; one passing pair retains the event.
    links.append(supported_link(events[0]['event_id'], 'recipient', 'A_good_copy'))
    saved_pfam(tmp_path, dict(D_gene=['PF01053'], A_gene=['PF00001'], A_good_copy=['PF01053'],
                              D_other=['PF00001'], B_gene=['PF00001']))
    selected, audit, _, _, _ = filter_events(events + [other], links, tmp_path)
    assert len(selected) == 2 and audit[0]['pfam_compared_pair_count'] == 2
    assert audit[0]['pfam_shared_accessions'] == 'PF01053'


def test_pfam_filter_never_borrows_another_family_record_and_supports_legacy_filename(tmp_path):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    saved_pfam(tmp_path, dict(D_gene=['PF01053'], A_gene=['PF01053']), family='OG2')
    assert not filter_events(events, links, tmp_path, allow_both_no_pfam=True)[0]
    saved_pfam(tmp_path, dict(D_gene=['PF01053'], A_gene=['PF01053']), legacy=True)
    assert len(filter_events(events, links, tmp_path)[0]) == 1
    links[0]['orthogroup'] = 'OG2'
    with pytest.raises(ValueError, match='identity mismatch'):
        filter_events(events, links, tmp_path)


def test_pfam_filter_rejects_malformed_no_hit_or_model_records(tmp_path):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    path = saved_pfam(tmp_path, dict(D_gene=[], A_gene=[]))
    fields, rows = read_tsv(path)
    rows[0]['qstart'] = '1'
    write_tsv(path, fields, rows)
    with pytest.raises(ValueError, match='no-hit record'):
        filter_events(events, links, tmp_path, allow_both_no_pfam=True)
    path = saved_pfam(tmp_path, dict(D_gene=['PF01053'], A_gene=['PF01053']))
    fields, rows = read_tsv(path)
    rows[0]['stitle'] = 'Unidentified domain'
    write_tsv(path, fields, rows)
    with pytest.raises(ValueError, match='Unmapped Pfam'):
        filter_events(events, links, tmp_path)


def test_pfam_filter_reads_archived_query_hits_without_materializing_inputs(tmp_path):
    from focus_hgt_pfam import filter_events
    from gene_family_output_store import archive_completed_outputs
    _, events, links = focused_node_source()
    root = tmp_path/'family'
    path = saved_pfam(root, dict(D_gene=['PF01053'], A_gene=['PF01053']))
    # Existing store completion contract includes the rendered family tree.
    for subdir, filename in [('mafft', 'OG1_cds.aln.fa.gz'),
                             ('stat_branch', 'OG1_stat.branch.tsv'), ('stat_tree', 'OG1_stat.tree.tsv'),
                             ('tree_plot', 'OG1_tree_plot.pdf')]:
        target = root/subdir/filename
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text('archive fixture\n')
    original = filter_events(events, links, root)
    archive_completed_outputs(root, 'query2family', ['OG1'], lambda name: 'OG1' if name.startswith('OG1_') else None,
                              min_files=1)
    assert not path.exists()
    before = {p: p.read_bytes() for p in root.rglob('*') if p.is_file()}
    archived = filter_events(events, links, root)
    assert original == archived
    assert before == {p: p.read_bytes() for p in root.rglob('*') if p.is_file()}


def test_pfam_filter_flow_preserves_pre_filter_trait_count_and_validates_subset():
    from focus_hgt_figures import filtering_counts
    events = [dict(event_id=f'e{i}', orthogroup='OG1') for i in range(3)]
    counts = filtering_counts(events, events[:1], prefilter_selected=events[:2])
    assert [r['event_count'] for r in counts] == [3, 2, 1]
    assert counts[-1]['stage'] == 'Event-gene pair Pfam filter'
    with pytest.raises(ValueError, match='not a subset'):
        filtering_counts(events, events[2:], prefilter_selected=events[:2])
    counts = filtering_counts(events, events[:1], pfam_selected=events[:2])
    assert [r['stage'] for r in counts][-2:] == ['Event-gene pair Pfam filter', 'Category = 1 recipients']
    assert [r['event_count'] for r in counts] == [3, 2, 1]
    with pytest.raises(ValueError, match='not a subset'):
        filtering_counts(events, events[2:], pfam_selected=events[:2])


def test_pfam_filter_rejects_reserved_columns_even_without_category1_events(source):
    fields, events = read_tsv(source[0])
    for row in events:
        row['pfam_filter_status'] = 'passed'
    write_tsv(source[0], fields + ['pfam_filter_status'], events)
    source[3].write_text('species\tbinary\nA\t0\nB\t0\nC\t0\nD\t0\n')
    schema_path(source[3]).write_bytes(schema_payload(source[3].read_bytes(), {'binary': 'binary'}))
    with pytest.raises(ValueError, match='Reserved Pfam'):
        generate_filtered(*source, plots=False)


def test_pfam_filter_defaults_apply_to_tables_without_plotting_and_empty_exception_is_opt_in(source):
    events, links, tree, trait, output = source
    _, erows = read_tsv(events)
    _, lrows = read_tsv(links)
    lrows = [dict(supported_link(r['event_id'], r['side'], r['gene_id']), **r) for r in lrows]
    write_tsv(links, list(lrows[0]), lrows)
    root = events.parent/'families'
    saved_pfam(root, dict(D_gene=['PF01053'], A_gene=['PF01053'], B_gene=[], C_gene=[]))
    originals = events.read_bytes(), links.read_bytes()
    report = generate_filtered(*source, plots=False, gene_family_root=root)
    selected = read_tsv(output/'traits/binary/all_category1/events.tsv')[1]
    assert {r['event_id'] for r in selected} == {'OG1:3:1', 'OG1:3:2'}
    assert report['require_shared_pfam'] is True and report['allow_both_no_pfam'] is False
    assert len(read_tsv(output/'traits/binary/pfam_event_audit.tsv')[1]) == 3
    assert read_tsv(output/'traits/binary/tips/B/direct_events.tsv')[1] == []
    assert len(read_tsv(output/'traits/binary/tips/B/ancestral_recipient_events.tsv')[1]) == 1
    assert all(r['pfam_filter_status'] == 'passed' for r in selected)
    saved_pfam(root, dict(D_gene=[], A_gene=[], B_gene=[], C_gene=[]))
    generate_filtered(*source, plots=False, gene_family_root=root)
    assert read_tsv(output/'traits/binary/all_category1/events.tsv')[1] == []
    generate_filtered(*source, plots=False, gene_family_root=root, allow_both_no_pfam=True)
    assert len(read_tsv(output/'traits/binary/all_category1/events.tsv')[1]) == 3
    assert (events.read_bytes(), links.read_bytes()) == originals
    assert len(erows) == 5


@pytest.mark.parametrize('native_status', ['selected', 'withheld'])
def test_filtered_plot_consumers_receive_the_same_cohort_as_filtered_tables(source, monkeypatch, native_status):
    import focus_hgt_figures
    import focus_hgt_gene_trees
    import focus_hgt_traits
    fields, links = read_tsv(source[1])
    links = [dict(supported_link(r['event_id'], r['side'], r['gene_id']), **r) for r in links]
    write_tsv(source[1], list(links[0]), links)
    root = source[0].parent/'families'
    saved_pfam(root, dict(D_gene=['PF01053'], A_gene=['PF01053'], B_gene=[], C_gene=[]))
    original = focus_hgt_traits.export_bundle
    def without_species_pdf(*args, **kwargs):
        args = list(args)
        args[8] = False
        return original(*args, **kwargs)
    trees, figures = [], []
    def capture_tree(directory, events, links, family_root, **kwargs):
        trees.append({r['event_id'] for r in events})
        write_tsv(directory/'event_node_audit.tsv', ['event_id', 'status'],
                  [dict(event_id=r['event_id'], status=native_status) for r in events])
        return dict(rendered_family_count=1)
    def capture_figure(directory, source_events, selected, *args, **kwargs):
        counts = focus_hgt_figures.filtering_counts(source_events, selected, pfam_selected=kwargs['pfam_selected'])
        assert kwargs['analyzed_orthogroups'] is None  # The fixture lacks native branch summaries.
        figures.append([r['event_count'] for r in counts])
        return dict(pdf_count=3)
    monkeypatch.setattr(focus_hgt_traits, 'export_bundle', without_species_pdf)
    monkeypatch.setattr(focus_hgt_gene_trees, 'export_gene_trees', capture_tree)
    monkeypatch.setattr(focus_hgt_figures, 'export_figures', capture_figure)
    generate_filtered(*source, plots=True, gene_family_root=root)
    assert trees == [{'OG1:3:1', 'OG1:3:2'}, {'OG1:3:1'}]
    assert figures == [[5, 2, 2], [5, 2, 1]]
    for trait, expected in zip(('binary', 'category'), trees, strict=True):
        assert {r['event_id'] for r in read_tsv(source[-1]/f'traits/{trait}/all_category1/events.tsv')[1]} == expected


def test_shared_pfam_is_evaluated_once_before_multiple_traits(source, monkeypatch):
    import focus_hgt_pfam
    fields, links = read_tsv(source[1])
    links = [dict(supported_link(r['event_id'], r['side'], r['gene_id']), **r) for r in links]
    write_tsv(source[1], list(links[0]), links)
    root = source[0].parent/'families'
    saved_pfam(root, dict(D_gene=['PF01053'], A_gene=['PF01053'], B_gene=[], C_gene=[]))
    calls=[]
    original=focus_hgt_pfam.filter_events
    def capture(events, *args, **kwargs):
        calls.append({r['event_id'] for r in events})
        assert all('focus_trait' not in r for r in events)
        return original(events, *args, **kwargs)
    monkeypatch.setattr(focus_hgt_pfam, 'filter_events', capture)
    report=generate_filtered(*source, plots=False, gene_family_root=root)
    assert calls == [{f'OG1:3:{i}' for i in range(1,6)}]
    assert report['shared_pfam_filter']['input_event_count']==5
    assert report['shared_pfam_filter']['passed_event_count']==2
    assert report['filtering_order']==['input_cohort','pfam_pair_filter','trait_category1']
    assert len(read_tsv(source[-1]/'pfam_event_audit.tsv')[1])==5
    assert len(read_tsv(source[-1]/'pfam_events.tsv')[1])==2
    assert all(report['pfam_filter'][trait]['passed_event_count']==2 for trait in ('binary','category'))


def test_no_scaffold_supported_genes_do_not_consume_irrelevant_domain_inputs(tmp_path):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    for row in links:
        row['host_scaffold_status']='unavailable'
    path=tmp_path/'rpsblast/OG1_rpsblast.tsv'
    path.parent.mkdir()
    path.write_text('Legacy plotting-only fixture\n')
    selected, audit, pairs, genes, sources=filter_events(events,links,tmp_path)
    assert not selected and not pairs and not genes and not sources
    assert audit[0]['pfam_filter_reason']=='no_bilateral_scaffold_supported_gene_pair'


def focused_node_source():
    event = dict(event_id="OG1:3:1", orthogroup="OG1", gene_tree_branch_id="3", gene_tree_node="n3",
                 event_index="1", generax_transfer="Y@D@A", support_used="90")
    stat = [dict(branch_id="0", node_name="A_gene", child1="-999", child2="-999", generax_transfer="N",
                 support_generax_ufboot=""),
            dict(branch_id="3", node_name="n3", child1="0", child2="1", generax_transfer="Y@D@A;Y@D@B",
                 support_generax_ufboot="90"),
            dict(branch_id="1", node_name="D_gene", child1="-999", child2="-999", generax_transfer="N",
                 support_generax_ufboot="")]
    links = [supported_link(event["event_id"], side, gene) for side, gene in
             (("donor", "D_gene"), ("recipient", "A_gene"))]
    return stat, [event], links


def test_focused_gene_nodes_use_exact_event_branch_and_inclusive_bilateral_profile():
    stat, events, links = focused_node_source()
    output, audit = annotate(stat, events, links)
    assert audit[0]["status"] == "selected" and audit[0]["support_generax_ufboot"] == 90
    assert output[0]["hgtfocus_recipient_flag"] == 1
    assert output[1]["hgtfocus_event_ids"] == "OG1:3:1"
    assert output[1]["hgtfocus_node_label"] == "HGT1 UFB=90"
    assert stat[1].get("hgtfocus_event_ids") is None  # Original table is preserved.
    # Two events sharing a gene-tree node are evaluated separately.
    other = dict(events[0], event_id="OG1:3:2", event_index="2", generax_transfer="Y@D@B")
    extra = supported_link(other["event_id"], "recipient", "B_gene")
    output, audit = annotate(stat, events + [other], links + [extra])
    assert [row["status"] for row in audit] == ["selected", "withheld"]
    assert output[1]["hgtfocus_event_count"] == 1  # Missing donor evidence cannot borrow from another event.


@pytest.mark.parametrize("alteration,reason", [
    ("node", "gene_tree_branch_node_unmapped"),
    ("terminal", "terminal_gene_tree_branch"),
    ("token", "transfer_token_or_event_index_unmapped"),
    ("index", "transfer_token_or_event_index_unmapped"),
    ("missing_support", "generax_ufboot_unavailable"),
    ("low_support", "generax_ufboot_below_threshold"),
    ("recipient_background", "bilateral_class_background_not_supported_or_unavailable"),
    ("context_gene", "supported_context_gene_unmapped_to_family_tip"),
])
def test_focused_gene_nodes_withhold_missing_or_unmapped_evidence(alteration, reason):
    stat, events, links = focused_node_source()
    if alteration == "node":
        events[0]["gene_tree_node"] = "wrong"
    elif alteration == "terminal":
        stat[1]["child1"] = "-999"
    elif alteration == "token":
        events[0]["generax_transfer"] = "Y@X@A"
    elif alteration == "index":
        events[0]["event_index"] = "2"
    elif alteration.endswith("support"):
        stat[1]["support_generax_ufboot"] = "" if alteration == "missing_support" else "89.9"
    elif alteration == "context_gene":
        links[0]['gene_id'] = 'other_family_gene'
    else:
        links[1]["host_scaffold_background_class_total_count"] = ""
    output, audit = annotate(stat, events, links, minimum_ufboot=90)
    assert audit[0]["reason"] == reason
    assert all(row["hgtfocus_event_count"] == 0 for row in output)


def test_focused_gene_nodes_reject_conflicting_support_and_background_fractions():
    stat, events, links = focused_node_source()
    events[0]["support_used"] = "100"
    with pytest.raises(ValueError, match="support disagrees"):
        annotate(stat, events, links)
    links[0]["host_scaffold_background_class_compatible_fraction"] = "1"
    with pytest.raises(ValueError, match="background fraction"):
        background_supported(links[0])


def test_focused_gene_tree_folder_audits_missing_families_and_preserves_inputs(tmp_path):
    from focus_hgt_gene_trees import export_gene_trees

    stat, events, links = focused_node_source()
    root = tmp_path / "families"
    table = root / "stat_branch/OG1_stat.branch.tsv"
    write_tsv(table, list(stat[0]), stat)
    original = table.read_bytes()
    def renderer(annotated, pdf):
        assert read_tsv(annotated)[1][1]["hgtfocus_event_count"] == "1"
        pdf.write_bytes(b"%PDF-native-renderer-fixture")
    output = tmp_path / "tree_plot"
    report = export_gene_trees(output, events, links, root, renderer=renderer)
    assert report["rendered_family_count"] == report["selected_event_count"] == 1
    assert table.read_bytes() == original
    assert (output / "OG1_focused_hgt_tree_plot.pdf").is_file()
    missing = dict(events[0], orthogroup="OG2")
    report = export_gene_trees(tmp_path / "missing", [missing], [], root, renderer=renderer)
    assert report["selected_event_count"] == 0
    assert read_tsv(tmp_path / "missing/event_node_audit.tsv")[1][0]["reason"] == "stat_branch_unavailable"


@pytest.fixture
def source(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("((A:1,B:1)0042:1,(C:1,D:1)mixed:1)root;")
    trait = tmp_path / "traits.tsv"
    trait.write_text("species\tbinary\tcategory\tcontinuous\nA\t1\t1\t1\nB\t1\t2\t0\nC\t0\t1\t0.5\nD\tNA\tNA\t2\n")
    schema_path(trait).write_bytes(schema_payload(trait.read_bytes(),
                                                {"binary": "binary", "category": "categorical", "continuous": "numeric"}))
    events = []
    for number, recipient in enumerate(("A", "0042", "mixed", "B", "C"), 1):
        events.append(dict(event_id=f"OG1:3:{number}", orthogroup="OG1", generax_transfer=f"Y@D@{recipient}",
                           generax_donor_node="D", generax_recipient_node=recipient, mapping_status="matched",
                           support_generax_ufboot="90", existing_annotation="Product name", sequence_quality_status="NA"))
    event = tmp_path / "events.tsv"
    write_tsv(event, list(events[0]), events)
    links = []
    for row in events:
        allowed = {'A': {'A'}, 'B': {'B'}, 'C': {'C'}, '0042': {'A', 'B'}, 'mixed': {'C', 'D'}}
        for species, side in (("D", "donor"), ("A", "recipient"), ("B", "recipient"), ("C", "recipient")):
            links.append(dict(event_id=row["event_id"], orthogroup="OG1", gene_id=species + "_gene", gene_species=species,
                              side=side, eligible_for_context=str(side == 'donor' or species in allowed[row['generax_recipient_node']]),
                              product_name="Protein " + species, synteny_support_score="NA"))
    link = tmp_path / "links.tsv"
    write_tsv(link, list(links[0]), links)
    return event, link, tree, trait, tmp_path / "output"


def test_category_and_binary_focus_preserve_cohort_fields_and_event_identity(source):
    manifest = generate(*source, plots=False)
    root = source[-1]
    _, binary = read_tsv(root / "traits/binary/all_category1/events.tsv")
    assert {row["event_id"] for row in binary} == {"OG1:3:1", "OG1:3:2", "OG1:3:4"}
    _, category = read_tsv(root / "traits/category/all_category1/events.tsv")
    assert {row["event_id"] for row in category} == {"OG1:3:1", "OG1:3:5"}
    assert all(row["support_generax_ufboot"] == "90" and row["existing_annotation"] == "Product name" for row in binary)
    assert next(row for row in binary if row["generax_recipient_node"] == "0042")["recipient_clade_tip_labels"] == "A; B"
    assert all(row["focus_ancestral_state_inferred"] == "0" for row in binary)
    _, direct_a = read_tsv(root / "traits/binary/tips/A/direct_events.tsv")
    _, ancestors_a = read_tsv(root / "traits/binary/tips/A/ancestral_recipient_events.tsv")
    _, ancestors_b = read_tsv(root / "traits/binary/tips/B/ancestral_recipient_events.tsv")
    assert [row["event_id"] for row in direct_a] == ["OG1:3:1"]
    assert [row["event_id"] for row in ancestors_a] == [row["event_id"] for row in ancestors_b] == ["OG1:3:2"]
    assert not (root / "traits/binary/tips/D").exists()
    assert not (root / "traits/continuous").exists()
    assert manifest["source_event_count"] == 5
    assert next(row for row in manifest["trait_selection"] if row["trait"] == "continuous")["status"] == "excluded"
    _, recipient = read_tsv(root / "traits/binary/tips/A/recipient_genes.tsv")
    assert len(recipient) == 1 and recipient[0]["gene_species"] == "A"
    assert recipient[0]["product_name"] == "Protein A" and recipient[0]["synteny_support_score"] == "NA"


@pytest.mark.parametrize('eligible', ['true', 'TRUE', '1'])
def test_focus_gene_tables_use_the_same_eligibility_as_event_validation(source, eligible):
    fields, links = read_tsv(source[1])
    for row in links:
        if row['eligible_for_context'] == 'True':
            row['eligible_for_context'] = eligible
    write_tsv(source[1], fields, links)
    manifest = generate(*source, plots=False)
    target = source[-1] / 'traits/binary/tips/A'
    assert read_tsv(target / 'recipient_genes.tsv')[1][0]['gene_id'] == 'A_gene'
    assert read_tsv(target / 'donor_genes.tsv')[1][0]['gene_id'] == 'D_gene'
    summary = next(row for row in manifest['result_index'] if row['trait'] == 'binary' and row['target'] == 'A')
    assert summary['recipient_gene_count'] == summary['donor_gene_count'] == 1


def test_category1_gene_tree_export_receives_only_aggregate_events(source, monkeypatch):
    import focus_hgt_figures
    import focus_hgt_gene_trees
    import focus_hgt_traits

    calls = []
    real_export = focus_hgt_traits.export_bundle
    def export_without_species_pdf(*args, **kwargs):
        # Keep this wiring test independent of the species-tree renderer.
        positional = list(args)
        positional[8] = False
        return real_export(*positional, **kwargs)
    def capture(directory, events, links, root, gff_root='', context_annotations=''):
        assert context_annotations == str(annotation_path)
        calls.append((directory, {e['event_id'] for e in events}, root))
        write_tsv(directory/'event_node_audit.tsv', ['event_id','status'],
                  [dict(event_id=e['event_id'],status='selected') for e in events])
        return dict(rendered_family_count=1, selected_event_count=len(events))
    monkeypatch.setattr(focus_hgt_traits, 'export_bundle', export_without_species_pdf)
    monkeypatch.setattr(focus_hgt_gene_trees, 'export_gene_trees', capture)
    monkeypatch.setattr(focus_hgt_figures, 'export_figures', lambda *args, **kwargs:dict(pdf_count=3))
    annotation_path = source[0].parent/'existing_annotations.tsv'
    annotation_path.write_text('Existing normalized annotation input\n')
    report = generate(*source, plots=True, gene_family_root='existing-families', context_annotations=str(annotation_path))
    assert len(calls) == 2  # One aggregate per binary/categorical trait, no per-tip rendering.
    assert all(path.name == 'tree_plot' and root == 'existing-families' for path, _, root in calls)
    assert calls[0][1] == {'OG1:3:1', 'OG1:3:2', 'OG1:3:4'}
    assert calls[1][1] == {'OG1:3:1', 'OG1:3:5'}
    assert report['gene_tree_plots']['binary']['selected_event_count'] == 3
    assert report['summary_figures']['binary']['pdf_count'] == 3
    assert str(annotation_path.resolve()) in report['inputs_sha256']
    assert not list(source[-1].rglob('transfer_tree.pdf'))


def test_empty_category1_targets_are_reported(source):
    source[0].write_text(source[0].read_text().replace("Y@D@A", "Y@D@D").replace("\tD\tA\t", "\tD\tD\t"))
    fields, links = read_tsv(source[1])
    for link in links:
        if link['event_id'] == 'OG1:3:1' and link['side'] == 'recipient':
            link['eligible_for_context'] = 'False'
    write_tsv(source[1], fields, links)
    generate(*source, plots=False)
    summary = json.loads((source[-1] / "traits/category/tips/A/summary.json").read_text())
    assert summary["event_count"] == 0
    assert read_tsv(source[-1] / "traits/category/tips/A/events.tsv")[1] == []


@pytest.mark.parametrize("alteration,match", [
    ("duplicate_event", "Duplicate or empty transfer"),
    ("duplicate_link", "Duplicate or invalid"),
    ("family_mismatch", "identity mismatch"),
    ("tree_mismatch", "clade tips disagree"),
    ("tip_count_mismatch", "clade tip count disagrees"),
    ("token_mismatch", "transfer token disagrees"),
    ("schema_mismatch", "schema does not match"),
    ("bad_binary", "Invalid binary"),
    ("alias_duplicate", "unique nonempty species"),
    ("positive_outside_tree", "absent from analysis"),
])
def test_invalid_inputs_do_not_publish_partial_bundle(source, alteration, match):
    event, link, _, trait, out = source
    if alteration == "duplicate_event":
        with event.open("a") as handle:
            handle.write(event.read_text().splitlines()[1] + "\n")
    elif alteration == "duplicate_link":
        with link.open("a") as handle:
            handle.write(link.read_text().splitlines()[1] + "\n")
    elif alteration == "family_mismatch":
        link.write_text(link.read_text().replace("\tOG1\t", "\tOTHER\t", 1))
    elif alteration in {"tree_mismatch", "tip_count_mismatch"}:
        fields, rows = read_tsv(event)
        for row in rows:
            row["recipient_clade_tip_labels"] = "Z" if alteration == "tree_mismatch" else {"A": "A", "B": "B", "C": "C", "0042": "A; B", "mixed": "C; D"}[row["generax_recipient_node"]]
            row["recipient_clade_tip_count"] = "99"
        write_tsv(event, fields + ["recipient_clade_tip_labels", "recipient_clade_tip_count"], rows)
    elif alteration == "token_mismatch":
        event.write_text(event.read_text().replace("Y@D@A", "Y@C@A"))
    else:
        replacement = {"schema_mismatch": ("A\t1", "A\t0"), "bad_binary": ("A\t1", "A\t2"),
                       "alias_duplicate": ("B\t1", "A\t1"), "positive_outside_tree": ("A\t1", "Unknown\t1")}[alteration]
        trait.write_text(trait.read_text().replace(*replacement))
        if alteration != "schema_mismatch":
            schema = json.loads(schema_path(trait).read_text())["traits"]
            schema_path(trait).write_bytes(schema_payload(trait.read_bytes(), schema))
    with pytest.raises(ValueError, match=match):
        generate(*source, plots=False)
    assert not out.exists()
    assert not list(out.parent.glob(".hgt-trait-focus-*"))


def test_schema_free_binary_only_and_missing_retained(source):
    schema_path(source[3]).unlink()
    generate(*source, plots=False)
    root = source[-1]
    assert (root / "traits/binary").exists()
    assert not (root / "traits/category").exists()
    assert not (root / "traits/continuous").exists()
    _, values = read_tsv(root / "traits/binary/species_trait_category1.tsv")
    assert next(row for row in values if row["species"] == "D")["binary"] == ""


def test_focused_export_retains_edge_tables_and_omits_unrequested_species_pdfs(source):
    manifest = generate(*source, plots=True, arrow_alpha=0.4)
    assert manifest["transfer_arrow_alpha"] == 0.4
    assert manifest["plot_scope"] == "trait_aggregate_only"
    root = source[-1] / "traits/binary/all_category1"
    _, edges = read_tsv(root / "transfer_edges.tsv")
    assert sum(int(row["hgt_event_count"]) for row in edges) == 3
    assert not list(source[-1].rglob('transfer_tree.pdf'))
    for target in manifest["result_index"]:
        if target["target_type"] != "aggregate":
            directory = source[-1] / target["relative_path"]
            _, target_edges = read_tsv(directory / "transfer_edges.tsv")
            assert sum(int(row["hgt_event_count"]) for row in target_edges) == target["event_count"]
            assert (directory / "events.tsv").is_file()
    from plot_hgt_summary import read_transfer_traits
    assert list(read_transfer_traits(source[3]).columns) == ["binary", "continuous"]
    # Republishing a managed legacy bundle removes old per-recipient PDFs.
    legacy = source[-1] / "traits/binary/tips/A/transfer_tree.pdf"
    legacy.write_bytes(b"previous recipient figure")
    generate(*source, plots=True, arrow_alpha=0.4)
    assert not legacy.exists()
    assert (legacy.parent / "events.tsv").is_file()
    assert (legacy.parent / "transfer_edges.tsv").is_file()


def test_republication_replaces_managed_bundle_and_failure_preserves_it(source, monkeypatch):
    import focus_hgt_traits
    from focus_hgt_traits import build_focus
    generate(*source, plots=False)
    first = (source[-1] / "manifest.json").read_bytes()
    generate(*source, plots=False)
    assert (source[-1] / "manifest.json").read_bytes() == first
    def fail(*args, **kwargs):
        build_focus(*args, **kwargs)
        raise ValueError("Intentional failure")
    monkeypatch.setattr(focus_hgt_traits, "build_focus", fail)
    with pytest.raises(ValueError, match="Intentional failure"):
        generate(*source, plots=False)
    assert (source[-1] / "manifest.json").read_bytes() == first


def test_unmanaged_output_and_input_containment_are_rejected(source):
    source[-1].mkdir()
    with pytest.raises(ValueError, match="unmanaged"):
        generate(*source, plots=False)
    with pytest.raises(ValueError, match="contain an input"):
        generate(*source[:-1], source[0].parent, plots=False)


@pytest.mark.parametrize('parameter', ['gene_family_root', 'gff_root'])
@pytest.mark.parametrize('aliased', [False, True])
def test_managed_output_cannot_replace_existing_family_or_gff_inputs(source, monkeypatch, parameter, aliased):
    import focus_hgt_traits

    output = source[-1]
    output.mkdir()
    manifest = output / 'manifest.json'
    manifest.write_text(json.dumps({'schema_version': focus_hgt_traits.VERSION}))
    inputs = output / 'existing_inputs'
    inputs.mkdir()
    sentinel = inputs / 'curated.tsv'
    sentinel.write_text('existing input must survive\n')
    root = inputs
    if aliased:
        root = output.parent / 'source_alias'
        root.symlink_to(inputs, target_is_directory=True)
    # Avoid scientific rendering; the publication guard must fire before dispatch.
    def synthetic_build(*args, **kwargs):
        return {'schema_version': focus_hgt_traits.VERSION}
    monkeypatch.setattr(focus_hgt_traits, 'build_focus', synthetic_build)
    before = manifest.read_bytes()
    with pytest.raises(ValueError, match='contain an input'):
        generate(*source, **{parameter: str(root)})
    assert sentinel.read_text() == 'existing input must survive\n'
    assert manifest.read_bytes() == before


def test_observation_columns_are_not_focus_traits(source):
    import hashlib

    from species_trait_contract import CONTRACT_VERSION
    trait = source[3]
    metadata = {"schema_version": CONTRACT_VERSION, "table_sha256": hashlib.sha256(trait.read_bytes()).hexdigest(),
                "traits": {"binary": {"role": "quality", "source": "user"}}}
    Path(str(trait) + ".metadata.json").write_text(json.dumps(metadata))
    generate(*source, plots=False)
    assert not (source[-1] / "traits/binary").exists()


def test_categorical_only_table_exports_focused_tables_without_species_pdf(source):
    from plot_hgt_summary import read_transfer_traits
    trait = source[3]
    trait.write_text("species\tcategory\nA\t1\nB\t2\nC\t1\nD\tNA\n")
    schema_path(trait).write_bytes(schema_payload(trait.read_bytes(), {"category": "categorical"}))
    assert read_transfer_traits(trait).shape == (4, 0)
    generate(*source, plots=True)
    assert (source[-1] / "traits/category/all_category1/events.tsv").is_file()
    assert not (source[-1] / "traits/category/all_category1/transfer_tree.pdf").exists()


def test_changed_input_preserves_previous_bundle(source, monkeypatch):
    import focus_hgt_traits
    generate(*source, plots=False)
    previous = (source[-1] / "manifest.json").read_bytes()
    build = focus_hgt_traits.build_focus
    def change_source(*args, **kwargs):
        result = build(*args, **kwargs)
        with source[0].open("a") as handle:
            handle.write("\n")
        return result
    monkeypatch.setattr(focus_hgt_traits, "build_focus", change_source)
    with pytest.raises(ValueError, match="inputs changed"):
        generate(*source, plots=False)
    assert (source[-1] / "manifest.json").read_bytes() == previous


def test_failed_publication_restores_previous_bundle(source, monkeypatch):
    import focus_hgt_traits
    generate(*source, plots=False)
    previous = (source[-1] / "manifest.json").read_bytes()
    replace = focus_hgt_traits.os.replace
    def fail_stage(path, destination):
        if Path(path).name.startswith(".hgt-trait-focus-") and not Path(path).name.startswith(".hgt-trait-focus-backup-"):
            raise OSError("publication failed")
        return replace(path, destination)
    monkeypatch.setattr(focus_hgt_traits.os, "replace", fail_stage)
    with pytest.raises(OSError, match="publication failed"):
        generate(*source, plots=False)
    assert (source[-1] / "manifest.json").read_bytes() == previous


def test_context_structures_keep_genomic_units_and_utr_boundaries():
    from focus_hgt_context import structure
    row = dict(feature_blocks='300-350;100-150', utr_blocks='50-99;351-400',
               feature_type='CDS', strand='-', start='50', end='400', splice_mode='cis')
    result = structure(row)
    assert result['coding'] == [(100, 150), (300, 350)]
    assert result['utr'] == [(50, 99), (351, 400)]
    assert result['introns'] == [(151, 299)]
    assert structure(dict(row, splice_mode='trans'))['introns'] == []
    assert structure(dict(row, feature_blocks='', utr_blocks=''))['status'] == 'exon_coordinates_unavailable'
    with pytest.raises(ValueError, match='outside'):
        structure(dict(row, feature_blocks='1-150'))


def test_renderer_argument_record_roundtrip_and_exact_alignment_tip_validation(tmp_path):
    from gene_family_output_store import GeneFamilyOutputStore, read_only_observation
    from gene_tree_plot_config import record, replay
    root = tmp_path/'families'
    path = root/'artifact_provenance/OG1.tree_plot.args.json'
    args = [f'--stat_branch={root}/stat_branch/OG1_stat.branch.tsv',
            '--panel1=tree,bl_rooted,support_unrooted,l1ou_regime,L',
            f'--panel2=domain,{root}/rpsblast/OG1_rpsblast.tsv',
            '--panel3=synteny_similarity,missing.tsv,20,15',
            '--panel4=sequence_similarity,missing.fa,cds,15,1',
            '--panel5=synteny,missing.tsv,5',
            '--long_branch_display=auto', '--event_method=auto']
    record(path, root, args, 'taxonomic')
    domain = root/'rpsblast/OG1_rpsblast.tsv'
    domain.parent.mkdir()
    domain.write_text('gene\tdomain\nA_gene\tD1\n')
    rows, _, _ = focused_node_source()
    with read_only_observation():
        store = GeneFamilyOutputStore(root)
        with store.read_snapshot():
            report = replay(store, 'OG1', rows, tmp_path/'materialized', {})
    assert '--panel1=tree,bl_rooted,support_unrooted,l1ou_regime,L' in report['arguments']
    assert '--long_branch_display=auto' in report['arguments']
    assert not any(arg.startswith('--stat_branch=') for arg in report['arguments'])
    assert report['settings_source'] == 'recorded_gg_gene_evolution_arguments'
    assert not any('=sequence_similarity,' in arg or '=synteny_similarity,' in arg for arg in report['arguments'])
    assert '--panel3=synteny,missing.tsv,5' in report['arguments']
    assert report['focused_disabled_panels'] == ['synteny_similarity', 'sequence_similarity']
    assert (tmp_path/'materialized/rpsblast/OG1_rpsblast.tsv').read_bytes() == domain.read_bytes()


def test_context_page_has_shared_scale_and_preserves_missing_coordinates(tmp_path):
    from focus_hgt_context import GenomeCoordinates, render_context
    from pypdf import PdfReader
    stat, events, links = focused_node_source()
    for row in stat:
        row['parent'] = '3' if row['branch_id'] != '3' else '-999'
        row['bl_rooted'] = '.1'
    for link in links:
        link['gene_species'] = link['gene_id'].split('_')[0]
    gff = tmp_path/'gff'
    gff.mkdir()
    row = dict(gene_id='D_gene', feature_size='102', num_intron='1', chromosome='scaffold1',
               start='1', end='30001', strand='+', feature_blocks='1-51;29951-30001',
               utr_blocks='', feature_type='CDS', splice_mode='cis')
    # A nearby gene's very long intron must not shrink the focal display.
    neighbor = dict(row, gene_id='D_nearby', start='30002', end='500000',
                    feature_blocks='30002-30100;499900-500000')
    write_tsv(gff/'D.gff_info.tsv', list(row), [row, neighbor])
    pdf = tmp_path/'context.pdf'
    audit = render_context(pdf, stat, events, links, GenomeCoordinates(gff))
    assert len(PdfReader(pdf).pages) == 1
    assert len({(r['shared_axis_min_kb'],r['shared_axis_max_kb']) for r in audit}) == 1
    assert audit[0]['shared_axis_max_kb'] == 35
    assert 'extends beyond the display window' in PdfReader(pdf).pages[0].extract_text()
    assert audit[0]['structure_status'] == 'coding_exons_utr_unavailable'
    assert audit[1]['structure_status'] == 'gff_gene_unavailable'
    assert audit[1]['intron_count'] == ''  # Unknown is not zero introns.
    context_only = tmp_path/'context_only.pdf'
    simplified = render_context(context_only, stat, events, links, GenomeCoordinates(gff), gene_tree_panel=False)
    text = PdfReader(context_only).pages[0].extract_text()
    assert 'Existing gene-tree branch length' not in text
    assert 'substitution/site' not in text
    assert 'HGT node n3 | UFB 90' in text
    assert 'UFB = Ultrafast bootstrap' in text
    assert 'UF=' not in text
    assert [r['gene_id'] for r in simplified] == [r['gene_id'] for r in audit]
    assert [(r['shared_axis_min_kb'],r['shared_axis_max_kb']) for r in simplified] == [(r['shared_axis_min_kb'],r['shared_axis_max_kb']) for r in audit]


def test_bounded_context_caps_distinct_genes_and_audits_omitted_and_unknown_evidence(tmp_path):
    from focus_hgt_context import GenomeCoordinates, choose_context_genes, render_context
    from pypdf import PdfReader
    stat, events, original = focused_node_source()
    events.append(dict(events[0], event_id='OG1:3:2', event_index='2', generax_transfer='Y@D@B'))
    links = []
    for side, genes in [('donor', ['D_gene', 'D_copy1', 'D_copy2', 'D_copy3']),
                        ('recipient', ['A_gene', 'A_copy1', 'A_copy2', 'A_failed', 'A_unknown'])]:
        for gene in genes:
            for event in events:
                link = dict(supported_link(event['event_id'], side, gene), gene_species=gene[0])
                if gene == 'A_failed':
                    link.update(host_scaffold_background_class_compatible_count='8',
                                host_scaffold_background_class_incompatible_count='2',
                                host_scaffold_background_class_compatible_fraction='0.8')
                elif gene == 'A_unknown':
                    link.update(host_scaffold_status='unavailable', host_scaffold_id='',
                                host_scaffold_background_class_total_count='',
                                host_scaffold_background_class_classified_fraction='',
                                host_scaffold_background_class_compatible_fraction='')
                links.append(link)
            if gene not in {'A_gene', 'D_gene'}:
                stat.append(dict(stat[0], branch_id=str(len(stat)+10), node_name=gene))
    links.append(dict(links[0]))  # Duplicate listings do not count or render twice.
    coordinates = GenomeCoordinates(tmp_path/'gff')
    selected, audit, totals = choose_context_genes(events, links, coordinates)
    assert totals['donor'] == dict(total=4, shown=3, omitted=1, supported=4)
    assert totals['recipient'] == dict(total=5, shown=3, omitted=2, supported=3)
    assert len(selected) == 6 and len(audit) == 18  # Two exact event references per distinct side/gene.
    assert len({(r['side'], r['link']['gene_id']) for r in selected}) == 6
    assert {r['link']['gene_id'] for r in selected if r['side'] == 'recipient'} == {'A_gene', 'A_copy1', 'A_copy2'}
    failed = next(r for r in audit if r['gene_id'] == 'A_failed')
    unknown = next(r for r in audit if r['gene_id'] == 'A_unknown')
    assert failed['scaffold_support_status'] == 'scaffold_thresholds_not_met' and failed['displayed'] == 0
    assert unknown['scaffold_support_status'] == 'scaffold_evidence_unavailable'
    assert unknown['host_scaffold_background_class_classified_fraction'] == ''
    assert unknown['intron_count'] == ''
    pdf = tmp_path/'bounded.pdf'
    audit = render_context(pdf, stat, events, links, coordinates, gene_tree_panel=False, max_genes_per_side=3)
    reader = PdfReader(pdf)
    assert len(reader.pages) == 1
    assert float(reader.pages[0].mediabox.width) == 22 * 72
    text = reader.pages[0].extract_text()
    assert 'DONOR DESCENDANTS' in text and 'RECIPIENT DESCENDANTS' in text
    assert 'Shown 3 of 4 genes | 1 omitted' in text and 'Shown 3 of 5 genes | 2 omitted' in text
    assert 'substitution/site' not in text
    assert {(r['shared_axis_min_kb'], r['shared_axis_max_kb']) for r in audit} == {(-20, 20)}
    # With fewer passing recipients, a measured failure is shown distinctly from unknown evidence.
    smaller = [r for r in links if r['gene_id'] not in {'A_copy1', 'A_copy2'}]
    selected, audit, totals = choose_context_genes(events, smaller, coordinates)
    assert totals['recipient']['shown'] == 3
    assert {r['status'] for r in selected if r['side'] == 'recipient'} == {
        'scaffold_supported', 'scaffold_thresholds_not_met', 'scaffold_evidence_unavailable'}
    gff = tmp_path/'gff'
    gff.mkdir()
    gff_row = dict(gene_id='A_unknown', chromosome='recorded_scaffold', start='100', end='300',
                   strand='+', feature_type='CDS', feature_blocks='100-150;250-300', utr_blocks='', num_intron='1')
    write_tsv(gff/'A.gff_info.tsv', list(gff_row), [gff_row])
    _, audit, _ = choose_context_genes(events, smaller, GenomeCoordinates(gff))
    unknown = next(r for r in audit if r['gene_id'] == 'A_unknown')
    assert unknown['scaffold_support_status'] == 'scaffold_evidence_unavailable'
    assert unknown['scaffold'] == 'recorded_scaffold' and unknown['scaffold_basis'] == 'existing_gff'
    assert unknown['intron_count'] == '1'
    conflicting = dict(original[0], gene_species='D', host_scaffold_id='different_scaffold')
    with pytest.raises(ValueError, match='Conflicting scaffold evidence'):
        choose_context_genes(events, links + [conflicting], coordinates)
    for limit in [0, 4, True, 1.5]:
        with pytest.raises(ValueError, match='display limit'):
            choose_context_genes(events, links, coordinates, limit)
    with pytest.raises(ValueError, match='passing event-linked gene'):
        choose_context_genes(events, [r for r in links if r['side'] == 'donor'], coordinates)


def test_context_annotations_require_exact_gene_family_and_same_best_hit(tmp_path):
    from focus_hgt_context_annotations import RANKS, ContextAnnotations, annotation_cells

    row = dict(gene_id='A_gene', orthogroup='OG1', protein_product_name='Own product',
               protein_product_status='exact_transcript_product', besthit_accession='P12345',
               besthit_organism='Hit species', swissprot_best_hit_protein_name='Different hit product',
               **{f'besthit_{rank}': '' for rank in RANKS})
    path = tmp_path/'annotations.tsv'
    write_tsv(path, list(row), [row])
    annotations = ContextAnnotations(path)
    exact = annotations.get('A_gene', 'OG1')
    cells = annotation_cells(exact, 'A', 'Focal')
    assert 'Different hit product' in cells[1] and 'Own product' not in cells[1]
    assert 'best-hit prediction' in cells[1] and 'GFF' not in cells[1]
    assert annotation_cells(dict(exact, swissprot_best_hit_protein_name=''), 'A', 'Focal')[1] == 'Unavailable'
    assert annotation_cells(dict(exact, besthit_accession=''), 'A', 'Focal')[1] == 'Unavailable'
    assert 'Kingdom: unavailable' in cells[2]  # No name-based taxonomy inference.
    missing = annotations.get('A_neighbor')
    assert missing['besthit_accession'] == '' and missing['protein_product_name'] == ''
    with pytest.raises(ValueError, match='family/gene mapping'):
        annotations.get('A_gene', 'wrong_family')
    with pytest.raises(ValueError, match='best hit disagrees'):
        annotations.get('A_gene', 'OG1', dict(node_name='A_gene', child1='-999', child2='-999', sprot_best='DifferentHit'))
    with pytest.raises(ValueError, match='organism disagrees'):
        annotations.get('A_gene', 'OG1', dict(node_name='A_gene', child1='-999', child2='-999', sprot_best='P12345', organism='Wrong organism'))
    write_tsv(path, list(row), [dict(row, protein_product_name='Changed')])
    with pytest.raises(ValueError, match='changed during rendering'):
        annotations.verify()
    write_tsv(path, list(row), [row, row])
    with pytest.raises(ValueError, match='Duplicate or empty'):
        ContextAnnotations(path)
    write_tsv(path, list(row), [dict(row, besthit_accession='')])
    with pytest.raises(ValueError, match='same hit accession'):
        ContextAnnotations(path)
    slim = {k: v for k, v in row.items() if not k.startswith('protein_product_')}
    write_tsv(path, list(slim), [slim])
    assert ContextAnnotations(path).get('A_gene', 'OG1')['swissprot_best_hit_protein_name'] == 'Different hit product'


def test_context_annotation_page_six_tracks_with_neighbors_and_long_rank_names(tmp_path):
    from focus_hgt_context import GenomeCoordinates, render_context
    from focus_hgt_context_annotations import (
        RANKS,
        TABLE_EDGES,
        TABLE_WIDTH_PT,
        ContextAnnotations,
        annotation_cells,
        text_width,
    )
    from pypdf import PdfReader

    stat, events, _ = focused_node_source()
    gff = tmp_path/'gff'
    gff.mkdir()
    links, records, expected = [], [], set()
    for side, species in [('donor', 'D'), ('recipient', 'A')]:
        gff_rows = []
        for i in range(3):
            gene = species + ('_gene' if i == 0 else f'_copy{i}')
            scaffold = f'scaffold{i}'
            links.append(dict(supported_link(events[0]['event_id'], side, gene), gene_species=species,
                              host_scaffold_id=scaffold))
            if gene not in {r['node_name'] for r in stat}:
                stat.append(dict(stat[0], branch_id=str(len(stat)+10), node_name=gene))
            ids = [f'{gene}_neighbor{j}' for j in range(3)] + [gene] + [f'{gene}_neighbor{j}' for j in range(3, 6)]
            for j, ident in enumerate(ids):
                start = 1000 + j*3000
                gff_rows.append(dict(gene_id=ident, chromosome=scaffold, start=str(start), end=str(start+900),
                                     strand='+', feature_type='CDS', feature_blocks=f'{start}-{start+300};{start+700}-{start+900}',
                                     utr_blocks='', num_intron='1'))
                row = dict(gene_id=ident, orthogroup='OG1' if ident == gene else 'OtherOG',
                           protein_product_name='Transcriptional regulator (NtrC/NifA family)' if species == 'D' else '',
                           protein_product_status='locus_level_product' if species == 'D' else 'unavailable',
                           swissprot_best_hit_protein_name='U3 small nucleolar RNA-associated protein 6',
                           besthit_accession='P12345', besthit_organism='Schizosaccharomyces pombe (strain 972 / ATCC 24843) (Fission yeast)',
                           **{f'besthit_{rank}': value for rank, value in zip(RANKS,
                               ['Fungi', 'Ascomycota', 'Schizosaccharomycetes', 'Schizosaccharomycetales',
                                'Schizosaccharomycetaceae', 'Schizosaccharomyces'], strict=True)})
                records.append(row)
                if j in {1, 2, 3, 4, 5}:
                    expected.add((gene, ident))
        write_tsv(gff/f'{species}.gff_info.tsv', list(gff_rows[0]), gff_rows)
    path = tmp_path/'annotations.tsv'
    write_tsv(path, list(records[0]), records)
    annotations = ContextAnnotations(path)
    pdf = tmp_path/'contexts.pdf'
    render_context(pdf, stat, events, links, GenomeCoordinates(gff), gene_tree_panel=False,
                   max_genes_per_side=3, annotations=annotations)
    reader = PdfReader(pdf)
    assert len(reader.pages) == 1
    text = ' '.join(reader.pages[0].extract_text().split())
    assert text.count('Best-hit taxonomic ranks') == 6
    assert text.count('MMseqs2 classification') == 6
    assert 'best-hit prediction' in text and 'GFF product' not in text
    assert 'Transcriptional regulator (NtrC/NifA family)' not in text
    assert 'Schizosaccharomycetaceae' in text and 'Genus: Schizosaccharomyces' in text
    assert len(annotations.display_audit) == 30
    assert {(r['context_focal_gene_id'], r['gene_id']) for r in annotations.display_audit} == expected
    assert all(r['event_ids'] == events[0]['event_id'] for r in annotations.display_audit)
    assert all(r['orthogroup'] == 'OtherOG' for r in annotations.display_audit if r['context_role'] == 'neighbor')
    for record in records:
        cells = annotation_cells(record, record['gene_id'][0], 'Focal')
        for i, cell in enumerate(cells):
            assert all(text_width(line) <= (TABLE_EDGES[i+1]-TABLE_EDGES[i])*TABLE_WIDTH_PT-10
                       for line in cell.split('\n'))


def test_context_displays_own_mmseqs2_taxonomy_separately_from_swissprot(tmp_path):
    from focus_hgt_context_annotations import RANKS, ContextAnnotations, annotation_cells
    from scaffold_taxonomy import RANKS as HOST_RANKS

    raw, host = tmp_path/'raw', tmp_path/'host'
    raw.mkdir()
    host.mkdir()
    path = raw/'A_mmseqs2taxonomy.tsv'
    path.write_text('A_gene\t4792\tspecies\tPhytophthora nicotianae\n'
                    'A_neighbor\t6656\tphylum\tArthropoda\nA_unknown\t0\tno rank\tunclassified\n')
    rows = [dict(species='A', gene_id=gene, scaffold='scaffold', locus_id=gene, count_unit='gff_locus',
                 rank=rank, host_taxid=str(4792 if rank == 'species' else 1),
                 label='compatible' if gene == 'A_gene' else 'unresolved')
            for gene in ['A_gene', 'A_neighbor'] for rank in HOST_RANKS]
    write_tsv(host/'A_gene_taxonomy.tsv', list(rows[0]), rows)
    annotations = ContextAnnotations(mmseqs2_taxonomy_dir=raw, scaffold_taxonomy_dir=host)
    record = annotations.get('A_gene', species='A')
    record.update(besthit_accession='Q96T49', besthit_organism='Homo sapiens',
                  swissprot_best_hit_protein_name='Human predicted product', **{f'besthit_{rank}': '' for rank in RANKS})
    cells = annotation_cells(record, 'A', 'Focal')
    assert 'Phytophthora nicotianae' in cells[3] and 'Homo sapiens' not in cells[3]
    assert 'Homo sapiens' in cells[2] and 'Phytophthora nicotianae' not in cells[2]
    assert 'Host class: compatible' in cells[3] and 'taxid: 4792' in cells[3]
    neighbor = annotations.get('A_neighbor', species='A')
    assert neighbor['mmseqs2_lca_name'] == 'Arthropoda' and neighbor['mmseqs2_host_class_label'] == 'unresolved'
    assert annotations.get('A_unknown', species='A')['mmseqs2_status'] == 'unclassified'
    missing = annotations.get('A_missing', species='A')
    assert missing['mmseqs2_status'] == 'gene_record_unavailable' and missing['mmseqs2_lca_name'] == ''
    assert annotations.get('B_gene', species='B')['mmseqs2_status'] == 'source_unavailable'
    path.write_text(path.read_text() + 'A_copy\t1\tno rank\troot\n')
    with pytest.raises(ValueError, match='changed during rendering'):
        annotations.verify()
    for text, message in [('A_gene\t9606.1\tspecies\tHuman\n', 'taxid'),
                          ('A_gene\t9606\tspecies\tHuman\nA_gene\t4792\tspecies\tOther\n', 'Duplicate'),
                          ('A_gene\t9606\n', 'Malformed')]:
        path.write_text(text)
        with pytest.raises(ValueError, match=message):
            ContextAnnotations(mmseqs2_taxonomy_dir=raw).get('A_gene', species='A')


def test_context_neighbor_selection_reserves_two_flanks_and_adds_intron_hosting_gene(tmp_path):
    from focus_hgt_context import GenomeCoordinates

    focal = dict(gene_id='A_focal', chromosome='scaffold', start='2244212', end='2244850',
                 feature_type='CDS', feature_blocks='2244212-2244850')
    enclosing = dict(focal, gene_id='A_enclosing', start='2226370', end='2246543',
                     feature_blocks='2226370-2226660;2235466-2235743;2246441-2246543')
    rows = [focal, enclosing]
    rows += [dict(focal, gene_id=f'A_flank{i}', start=str(2240000-i*1000), end=str(2240100-i*1000)) for i in range(6)]
    rows += [dict(focal, gene_id=f'A_right{i}', start=str(2250000+i*1000), end=str(2250100+i*1000)) for i in range(6)]
    write_tsv(tmp_path/'A.gff_info.tsv', list(focal), rows)
    _, neighbors, _ = GenomeCoordinates(tmp_path).neighborhood(dict(gene_species='A', gene_id='A_focal'))
    assert 'A_enclosing' in {r['gene_id'] for r in neighbors}
    assert len(neighbors) == 6  # Focal, two left, two right, one overlapping.
    from focus_hgt_context import neighbor_counts
    counts = neighbor_counts(focal, neighbors)
    assert counts['neighbor_left_gene_count'] == counts['neighbor_right_gene_count'] == 2


def test_all_models_in_view_are_drawn_without_expanding_annotation_table(tmp_path):
    from focus_hgt_context import GenomeCoordinates, model_span, prepare_context_models
    from focus_hgt_context_annotations import ContextAnnotations, context_annotation_rows
    focal = dict(gene_id='A_focal', chromosome='s1', start='10000', end='10100',
                 feature_type='CDS', feature_blocks='10000-10100', utr_blocks='')
    rows = [focal]
    for i, start in enumerate([1000,2000,3000,7000,11000,12000,13000,14000,15000]):
        rows.append(dict(focal, gene_id=f'A_neighbor{i}', start=str(start), end=str(start+100),
                         feature_blocks=f'{start}-{start+100}'))
    rows += [dict(focal, gene_id=f'A_nested{i}', start='10010', end='10050',
                  feature_blocks='10010-10050') for i in range(4)]
    rows += [dict(focal, gene_id='A_far', start='200000', end='200100', feature_blocks='200000-200100'),
             dict(focal, gene_id='A_wrong_scaffold', chromosome='s2')]
    write_tsv(tmp_path/'A.gff_info.tsv', list(focal), rows)
    coords = GenomeCoordinates(tmp_path)
    link = dict(gene_id='A_focal', gene_species='A')
    _, neighbors, _ = coords.neighborhood(link)
    entry = dict(focal=focal, neighbors=neighbors, link=link, side='recipient', event_ids={'e1'})
    extent = prepare_context_models([entry], coords)
    expected = {r['gene_id'] for r in rows if r['chromosome']=='s1'
                and entry['display_coordinates'].point(model_span(r)[0]) <= extent
                and entry['display_coordinates'].point(model_span(r)[1]) >= -extent}
    assert {r['gene_id'] for r in entry['models']} == expected
    assert len(entry['models']) > len(neighbors)
    assert 'A_wrong_scaffold' not in expected
    assert len({entry['model_lanes'][f'A_nested{i}'] for i in range(4)}) == 4
    assert entry['model_height_pt'] > 62
    for row in entry['models']:
        assert entry['display_coordinates'].point(int(row['end'])+1)-entry['display_coordinates'].point(int(row['start'])) == pytest.approx(
            (int(row['end'])+1-int(row['start']))/1000)
    tables = context_annotation_rows(entry, ContextAnnotations(), 'OG1', {})
    assert {r['gene_id'] for r in tables} == {r['gene_id'] for r in neighbors}


def test_mmseqs2_query_rank_names_follow_saved_lineage_not_host_or_best_hit(tmp_path):
    import sqlite3

    from focus_hgt_context_annotations import RANKS, ContextAnnotations, annotation_cells
    db = tmp_path/'taxa.sqlite'
    with sqlite3.connect(db) as conn:
        conn.execute('CREATE TABLE species(taxid INTEGER PRIMARY KEY,spname TEXT,rank TEXT,track TEXT)')
        conn.execute('CREATE TABLE merged(taxid_old INTEGER,taxid_new INTEGER)')
        for i, rank in enumerate(RANKS, 1):
            conn.execute('INSERT INTO species VALUES(?,?,?,?)', (i, 'Query_'+rank, rank, ','.join(map(str, range(i,0,-1)))))
        conn.execute('INSERT INTO species VALUES(7,"Query_species","species","7,6,5,4,3,2,1")')
    (tmp_path/'A_mmseqs2taxonomy.tsv').write_text('A_gene\t7\tspecies\tQuery_species\t1\t1\t1\t1\t1;2;3;4;5;6;7\n'
                                               'A_broad\t2\tphylum\tQuery_phylum\t1\t1\t1\t1\t1;2\n')
    before = db.read_bytes()
    annotations = ContextAnnotations(mmseqs2_taxonomy_dir=tmp_path, taxonomy_dbfile=db)
    row = annotations.get('A_gene', species='A')
    assert all(row['mmseqs2_'+rank] == 'Query_'+rank for rank in RANKS)
    broad = annotations.get('A_broad', species='A')
    assert broad['mmseqs2_kingdom']=='Query_kingdom' and broad['mmseqs2_genus']==''
    assert broad['mmseqs2_lineage_status']=='saved_lineage_resolved'
    cells = annotation_cells(dict(row, besthit_accession='P1', besthit_organism='Other organism',
                                  besthit_genus='Other genus', swissprot_best_hit_protein_name='Product'), 'A', 'Focal')
    assert 'Product' in cells[1] and 'Other organism' in cells[2] and 'Query_genus' in cells[3]
    assert 'Other genus' not in cells[3]
    annotations.verify()
    assert db.read_bytes()==before
    missing = ContextAnnotations(mmseqs2_taxonomy_dir=tmp_path).get('A_gene', species='A')
    assert missing['mmseqs2_lineage_status']=='taxonomy_source_unavailable' and missing['mmseqs2_genus']==''
    with sqlite3.connect(db) as conn:
        conn.execute('UPDATE species SET spname="Changed" WHERE taxid=6')
    with pytest.raises(ValueError,match='changed during rendering'):
        annotations.verify()


def test_context_two_flanks_ignore_distance_cutoff_and_keep_scaffolds_separate(tmp_path):
    from focus_hgt_context import GenomeCoordinates, neighbor_counts
    focal = dict(gene_id='A_focal', chromosome='s1', start='1000000', end='1001000',
                 feature_type='CDS', feature_blocks='1000000-1001000')
    rows = [focal]
    for name, start in [('left_far', 100), ('left_near', 900000), ('right_near', 1100000),
                        ('right_far', 2000000), ('right_extra', 3000000)]:
        rows.append(dict(focal, gene_id='A_'+name, start=str(start), end=str(start+100),
                         feature_blocks=f'{start}-{start+100}'))
    rows.append(dict(rows[1], gene_id='A_other_scaffold', chromosome='s2'))
    write_tsv(tmp_path/'A.gff_info.tsv', list(focal), rows)
    _, neighbors, _ = GenomeCoordinates(tmp_path).neighborhood(dict(gene_species='A', gene_id='A_focal'))
    assert {r['gene_id'] for r in neighbors} == {'A_focal','A_left_far','A_left_near','A_right_near','A_right_far'}
    counts = neighbor_counts(focal, neighbors)
    assert counts['neighbor_left_status'] == counts['neighbor_right_status'] == 'minimum_met'
    sparse = neighbor_counts(focal, [focal, rows[1]])
    assert sparse['neighbor_left_gene_count'] == 1 and sparse['neighbor_right_gene_count'] == 0
    assert sparse['neighbor_right_status'] == 'insufficient_annotated_loci'
    assert neighbor_counts(None, [])['neighbor_left_gene_count'] == ''


def test_gap_compression_preserves_unknown_gene_spans_and_overlapping_loci():
    from focus_hgt_context import GapCompressedCoordinates
    focal = dict(gene_id='A_focal', start='100000', end='130000', utr_blocks='99000-99999')
    enclosing = dict(gene_id='A_enclosing', start='95000', end='135000')
    left = dict(gene_id='A_left', start='1000', end='2000')
    right = dict(gene_id='A_right', start='1000000', end='1001000')
    display = GapCompressedCoordinates(focal, [left, enclosing, focal, right])
    assert display.point(115000) == 0
    for row in [left, enclosing, focal, right]:
        assert display.point(int(row['end'])+1) - display.point(int(row['start'])) == pytest.approx(
            (int(row['end'])+1-int(row['start']))/1000)
    assert display.point(129000)-display.point(101000) == pytest.approx(28)  # Introns remain genomic-length.
    gaps = display.audit()
    assert len(gaps) == 2
    assert all(g['display_end_kb']-g['display_start_kb'] == pytest.approx(2) for g in gaps)
    assert gaps[0]['genomic_start_bp'] == 2001 and gaps[0]['genomic_end_exclusive_bp'] == 95000
    assert gaps[0]['omitted_bp'] == 92999 - 2000
    positions = [1000,2001,50000,95000,115000,135001,500000,1000000,1001001]
    assert [display.point(x) for x in positions] == sorted(display.point(x) for x in positions)


def test_compressed_introns_preserve_nested_focal_exon_and_true_coordinate_audit():
    from focus_hgt_context import GapCompressedCoordinates
    host = dict(gene_id='A_host', start='100', end='100000', feature_type='CDS',
                feature_blocks='100-200;99900-100000')
    focal = dict(gene_id='A_focal', start='50000', end='50638', feature_type='CDS',
                 feature_blocks='50000-50638')
    display = GapCompressedCoordinates(focal, [host, focal])
    assert display.point(50319) == 0
    assert display.point(50639)-display.point(50000) == pytest.approx(.639)
    assert display.point(201)-display.point(100) == pytest.approx(.101)
    assert len(display.audit()) == 2 and all(g['gap_type'] == 'intronic' for g in display.audit())
    assert sum(g['omitted_bp'] for g in display.audit()) == (50000-201)+(99900-50639)-4000
    assert all(g['display_end_kb']-g['display_start_kb'] == pytest.approx(2) for g in display.audit())


def test_filter_flow_validates_event_grain_and_does_not_invent_upstream_counts(tmp_path):
    from focus_hgt_figures import filtering_counts, product_labels
    events = [dict(event_id='e1',orthogroup='OG1'),dict(event_id='e2',orthogroup='OG1')]
    assert filtering_counts(events, events[:1])[0]['orthogroup_count'] == 1
    assert len(filtering_counts(events, events[:1])) == 2
    rows = [dict(event_id='e1',orthogroup='OG1',donor_classification='outside',recipient_classification='insect',status='accepted',support_used='90',support_source='support_generax_ufboot')]
    path = tmp_path/'audit.tsv'
    write_tsv(path,list(rows[0]),rows)
    with pytest.raises(ValueError,match='not a subset'):
        filtering_counts(events,events[:1],path)
    names = product_labels(['OG1'],[dict(orthogroup='OG1',side='recipient',gene_id='g1',protein_product_name='Test protein',protein_product_source='existing GFF')])
    assert names['OG1']['protein_product'] == 'Annotation unavailable'
    assert names['OG1']['annotation_gene_id'] == ''
    fallback = product_labels(['OG1'], [dict(orthogroup='OG1', side='recipient', gene_id='g1',
                              protein_product_name='NA', protein_product_source='NA',
                              best_available_product_label='Existing best-hit name',
                              best_available_product_label_basis='existing best-hit annotation')])['OG1']
    assert fallback['protein_product'] == 'Existing best-hit name'
    assert fallback['all_recipient_product_labels'] == 'Existing best-hit name'
    assert fallback['annotation_basis'] == 'existing best-hit annotation'
    mixed = product_labels(['OG1'], [dict(orthogroup='OG1', side='recipient', gene_id='g1',
                           protein_product_name='Recorded product'),
                           dict(orthogroup='OG1', side='recipient', gene_id='g2', protein_product_name='NA',
                           swissprot_best_hit_protein_name='Existing best-hit name')])['OG1']
    assert mixed['protein_product'] == 'Existing best-hit name'
    assert mixed['annotation_gene_id'] == 'g2'
    assert mixed['annotation_basis'] == 'SwissProt_best_hit_prediction'
    assert mixed['all_recipient_product_labels'] == 'Existing best-hit name'
    unverified = product_labels(['OG1'], [dict(orthogroup='OG1', side='recipient', gene_id='g1',
                                 best_available_product_label='Unspecified annotation')])['OG1']
    assert unverified['protein_product'] == 'Annotation unavailable'
    retained = product_labels(['OG1'], [
        dict(orthogroup='OG1', side='recipient', gene_id='A_transferred_out',
             eligible_for_context='False', swissprot_best_hit_protein_name='Excluded descendant hit'),
        dict(orthogroup='OG1', side='recipient', gene_id='Z_retained',
             eligible_for_context='True', swissprot_best_hit_protein_name='Retained descendant hit')])['OG1']
    assert retained['annotation_gene_id'] == 'Z_retained'
    assert retained['all_recipient_product_labels'] == 'Retained descendant hit'


def test_context_annotations_reject_malformed_rows_and_fractional_taxids(tmp_path):
    from focus_hgt_context_annotations import FIELDS, ContextAnnotations, taxid

    path = tmp_path/'annotations.tsv'
    row = dict.fromkeys(FIELDS, '')
    row.update(gene_id='A_gene', orthogroup='OG1', besthit_accession='P1', besthit_taxid='9606.0')
    write_tsv(path, FIELDS, [row])
    assert ContextAnnotations(path).get('A_gene')['besthit_taxid'] == '9606.0'
    assert taxid('9606.0') == 9606 and taxid('.') is None
    for invalid in ('9606.1', 'Infinity', 'abc', '0', '-1'):
        with pytest.raises(ValueError, match='Invalid best-hit taxid'):
            taxid(invalid)
    for mutation, reason in [(dict(row, orthogroup=''), 'orthogroup'),
                             (dict(row, besthit_taxid='9606.1'), 'taxid')]:
        write_tsv(path, FIELDS, [mutation])
        with pytest.raises(ValueError, match=reason):
            ContextAnnotations(path)
    write_tsv(path, FIELDS, [row])
    valid = path.read_text()
    for malformed in [valid.replace('\n', '\textra\n', 1), valid.rstrip('\n')+'\textra\n',
                      valid.rsplit('\t', 1)[0]+'\n', valid.replace('gene_id\t', 'gene_id\tgene_id\t', 1)]:
        path.write_text(malformed)
        with pytest.raises(ValueError, match='Malformed|Duplicate|columns'):
            ContextAnnotations(path)


def test_context_annotations_verify_neighbor_own_family_and_leaf_identity(tmp_path):
    from focus_hgt_context_annotations import FIELDS, ContextAnnotations, context_annotation_rows
    from gene_family_output_store import GeneFamilyOutputStore, read_only_observation

    root = tmp_path/'families'
    root.joinpath('stat_branch').mkdir(parents=True)
    leaf = dict(node_name='A_neighbor', child1='-999', child2='-999', sprot_best='P1',
                sprot_recname='Correct name', organism='Hit organism', taxid_y='9606.0')
    write_tsv(root/'stat_branch/OG2_stat.branch.tsv', list(leaf), [leaf])
    row = dict.fromkeys(FIELDS, '')
    row.update(gene_id='A_neighbor', orthogroup='OG2', besthit_accession='P1', besthit_taxid='9606',
               swissprot_best_hit_protein_name='Correct name')
    path = tmp_path/'annotations.tsv'
    write_tsv(path, FIELDS, [row])
    with read_only_observation():
        annotations = ContextAnnotations(path, store=GeneFamilyOutputStore(root))
        exact = annotations.get('A_neighbor')
        assert exact['besthit_organism'] == 'Hit organism'
        assert exact['annotation_validation_status'] == 'exact_family_leaf_verified'
        assert 'stat_branch/OG2_stat.branch.tsv' in annotations.family_sources
        annotations.verify()
        with pytest.raises(ValueError, match='exact family leaf'):
            annotations.get('A_neighbor', 'OG2', dict(leaf, node_name='A_different'))
        with pytest.raises(ValueError, match='protein name disagrees'):
            annotations.get('A_neighbor', 'OG2', dict(leaf, sprot_recname='Wrong name'))
        with pytest.raises(ValueError, match='best hit disagrees'):
            annotations.get('A_neighbor', 'OG2', dict(leaf, sprot_best=''))
        # A neighbor in the displayed tree cannot claim a different OG.
        entry = dict(link=dict(gene_id='A_focal', gene_species='A'), side='recipient',
                     neighbors=[dict(gene_id='A_neighbor', start='1')], event_ids={'e1'})
        with pytest.raises(ValueError, match='family/gene mapping'):
            context_annotation_rows(entry, annotations, 'OG1', {'A_neighbor': leaf})
    write_tsv(root/'stat_branch/OG2_stat.branch.tsv', list(leaf), [dict(leaf, node_name='A_other')])
    with read_only_observation():
        with pytest.raises(ValueError, match='absent from its own family'):
            ContextAnnotations(path, store=GeneFamilyOutputStore(root)).get('A_neighbor')


def test_exon_only_structure_is_not_labeled_or_drawn_as_utr():
    from focus_hgt_context import structure

    row = dict(feature_type='exon', feature_blocks='100-150;300-350', utr_blocks='', start='100', end='350')
    result = structure(row)
    assert result['status'] == 'annotated_exons_CDS_UTR_unavailable'
    assert result['coding'] == result['utr'] == []
    assert result['exon'] == [(100, 150), (300, 350)]
    assert result['introns'] == [(151, 299)]
    cds = structure(dict(row, feature_type='CDS', utr_blocks='50-99;351-400'))
    assert cds['utr'] == [(50, 99), (351, 400)]  # Saved CDS spans exclude flanking UTRs.


def test_filtering_cohorts_cannot_borrow_event_identity_or_duplicate_counts():
    from focus_hgt_figures import filtering_counts

    row = dict(event_id='e1', orthogroup='OG1', gene_tree_branch_id='3', generax_transfer='Y@D@A')
    with pytest.raises(ValueError, match='Duplicate'):
        filtering_counts([row, row], [row])
    with pytest.raises(ValueError, match='not a subset'):
        filtering_counts([row], [dict(row, event_id='e2')])
    for field in ('orthogroup', 'gene_tree_branch_id', 'generax_transfer'):
        with pytest.raises(ValueError, match='identity disagrees'):
            filtering_counts([row], [dict(row, **{field: 'wrong'})])
    with pytest.raises(ValueError, match='identity disagrees'):
        filtering_counts([row], [dict(event_id='e1', orthogroup='OG1', branch_id='wrong')])


def test_eligible_transfer_gene_cannot_borrow_another_species_branch(source):
    fields, links = read_tsv(source[1])
    links[2]['eligible_for_context'] = 'True'  # B_gene is outside the first event's A recipient.
    write_tsv(source[1], fields, links)
    with pytest.raises(ValueError, match='outside its species branch'):
        generate(*source, plots=False)


def test_native_gene_tree_column_shows_both_roles_and_unconfirmed_recipient():
    stat, events, links = focused_node_source()
    stat.append(dict(stat[0], node_name='A_copy', branch_id='4'))
    links.append(dict(supported_link(events[0]['event_id'], 'recipient', 'A_copy'), host_scaffold_status='unavailable'))
    annotated, _ = annotate(stat, events, links)
    by_name = {r['node_name']: r for r in annotated}
    assert by_name['D_gene']['hgtfocus_donor_flag'] == 1
    assert by_name['D_gene']['hgtfocus_tip_status'] == 'Scaffold-supported donor descendant'
    assert by_name['A_gene']['hgtfocus_tip_status'] == 'Scaffold-supported recipient descendant'
    assert by_name['A_copy']['hgtfocus_recipient_flag'] == 0
    assert by_name['A_copy']['hgtfocus_tip_status'] == 'Scaffold-unconfirmed recipient descendant'


def test_distribution_figure_uses_the_same_supplemental_hit_prediction(tmp_path):
    from io import StringIO

    from Bio import Phylo
    from focus_hgt_context_annotations import FIELDS
    from focus_hgt_figures import export_figures
    from pypdf import PdfReader

    stat, events, links = focused_node_source()
    for row in stat:
        row.update(sprot_best='P1' if row['node_name'] == 'A_gene' else '', sprot_recname='', organism='')
    root = tmp_path/'families'
    write_tsv(root/'stat_branch/OG1_stat.branch.tsv', list(stat[0]), stat)
    for link, species in zip(links, ['D', 'A'], strict=True):
        link['gene_species'] = species
    row = dict.fromkeys(FIELDS, '')
    row.update(gene_id='A_gene', orthogroup='OG1', besthit_accession='P1', swissprot_best_hit_protein_name='Exact predicted product')
    path = tmp_path/'annotations.tsv'
    write_tsv(path, FIELDS, [row])
    tree = Phylo.read(StringIO('(A:1,D:1)root;'), 'newick')
    for event in events:
        event.update(generax_donor_node='D', generax_recipient_node='A')
    output = tmp_path/'plots'
    report = export_figures(output, events, events, links, tree, {'A': 1, 'D': 0}, root, 'gall', context_annotations=path)
    assert report['pdf_count'] == 3 and str(path.resolve()) in report['context_annotation_source_sha256']
    distribution = read_tsv(output/'orthogroup_species_distribution.tsv')[1]
    assert all(r['protein_product'] == 'Exact predicted product' for r in distribution)
    assert len(list(output.glob('*.pdf'))) == 3
    assert all(len(PdfReader(pdf).pages) == 1 for pdf in output.glob('*.pdf'))


def test_missing_gff_coordinates_remain_unavailable_and_invalid_spans_fail(tmp_path):
    from focus_hgt_context import GenomeCoordinates

    root = tmp_path/'gff'
    row = dict(gene_id='A_gene', chromosome='s1', start='NA', end='100', feature_type='CDS', feature_blocks='')
    write_tsv(root/'A.gff_info.tsv', list(row), [row])
    focal, neighbors, status = GenomeCoordinates(root).neighborhood(dict(gene_id='A_gene', gene_species='A'))
    assert focal is None and neighbors == [] and status == 'gff_coordinates_unavailable'
    for span in [('0', '100'), ('101', '100')]:
        write_tsv(root/'A.gff_info.tsv', list(row), [dict(row, start=span[0], end=span[1])])
        with pytest.raises(ValueError, match='Invalid GFF coordinate span'):
            GenomeCoordinates(root).load('A')


@pytest.mark.parametrize('donor_end,recipient_end,minimum,passed', [
    (50, 50, .5, True), (49, 90, .5, False), (90, 49, .5, False),
    (1, 1, 0, True), (100, 100, 1, True)])
def test_shared_pfam_query_coverage_is_inclusive_on_both_genes(tmp_path, donor_end, recipient_end, minimum, passed):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    path = saved_pfam(tmp_path, dict(D_gene=['PF01053'], A_gene=['PF01053']))
    fields, rows = read_tsv(path)
    for row, end in zip(rows, (donor_end, recipient_end), strict=True):
        row['qend'] = str(end)
    write_tsv(path, fields, rows)
    selected, audit, pairs, _, _ = filter_events(events, links, tmp_path, min_shared_pfam_coverage=minimum)
    assert bool(selected) is passed
    assert pairs[0]['donor_shared_pfam_covered_aa'] == donor_end
    assert pairs[0]['recipient_shared_pfam_query_coverage'] == recipient_end / 100
    assert pairs[0]['coverage_status'] == ('passed' if passed else 'below_minimum')
    if not passed:
        assert audit[0]['pfam_filter_reason'] == 'shared_pfam_below_minimum_query_coverage'


def test_shared_domain_intervals_are_unioned_without_overlap_double_counting(tmp_path):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    path = saved_pfam(tmp_path, dict(D_gene=['PF01053'], A_gene=['PF01053']))
    fields, rows = read_tsv(path)
    hits = []
    for row in rows:
        hits.extend([dict(row, qstart='1', qend='30'), dict(row, qstart='11', qend='40')])
    write_tsv(path, fields, hits)
    selected, _, pairs, _, _ = filter_events(events, links, tmp_path)
    assert not selected and pairs[0]['donor_shared_pfam_covered_aa'] == 40
    hits.append(dict(rows[0], stitle='pfam00001, Other, Description', qstart='41', qend='100'))
    write_tsv(path, fields, hits)
    assert not filter_events(events, links, tmp_path)[0]  # Unshared domain cannot rescue coverage.


def test_coverages_cannot_be_borrowed_from_different_pairs(tmp_path):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    eid = events[0]['event_id']
    links += [supported_link(eid, 'donor', 'D_other'), supported_link(eid, 'recipient', 'A_other')]
    path = saved_pfam(tmp_path, dict(D_gene=['PF01053'], A_gene=['PF01053'],
                                    D_other=['PF00001'], A_other=['PF00001']))
    fields, rows = read_tsv(path)
    for row in rows:
        row['qend'] = '90' if row['qacc'] in {'D_gene', 'A_other'} else '20'
    write_tsv(path, fields, rows)
    selected, audit, pairs, _, _ = filter_events(events, links, tmp_path)
    assert not selected and len(pairs) == 4
    assert min(audit[0]['pfam_best_pair_donor_query_coverage'], audit[0]['pfam_best_pair_recipient_query_coverage']) == .2
    for row in rows:
        if row['qacc'] == 'A_gene':
            row['qend'] = '50'
    write_tsv(path, fields, rows)
    selected, audit, pairs, _, _ = filter_events(events, links, tmp_path)
    assert len(selected) == 1 and audit[0]['pfam_passing_pair_count'] == 1
    assert audit[0]['pfam_best_pair_donor_gene'] == 'D_gene'
    assert audit[0]['pfam_best_pair_recipient_gene'] == 'A_gene'


@pytest.mark.parametrize('minimum', [-.1, 1.1, float('nan'), float('inf'), 'bad', None])
def test_invalid_pfam_coverage_configuration_fails_even_for_empty_cohort(minimum):
    from focus_hgt_pfam import filter_events
    with pytest.raises(ValueError, match='finite fraction'):
        filter_events([], [], '', min_shared_pfam_coverage=minimum)


def test_domain_review_flags_are_not_hard_exclusions(tmp_path):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    path = saved_pfam(tmp_path, dict(D_gene=['PF00023', 'PF00001'], A_gene=['PF00023']))
    fields, rows = read_tsv(path)
    for row in rows:
        row.update(qlen='90', qend='60')
        if row['stitle'].startswith('pfam00023'):
            row['stitle'] = 'pfam00023, Ank, Ankyrin repeat'
    write_tsv(path, fields, rows)
    selected, _, pairs, genes, _ = filter_events(events, links, tmp_path)
    assert selected
    flags = pairs[0]['pair_attention_flags']
    assert 'shared_repeat_or_generic_binding_domain_only' in flags
    assert 'pfam_domain_sets_differ_architecture_review' in flags
    assert 'donor_short_query_protein_lt100aa' in flags and 'recipient_short_query_protein_lt100aa' in flags
    assert all(r['gene_attention_flags'] == 'short_query_protein_lt100aa' for r in genes)
    assert selected[0]['pfam_attention_flags'] == '; '.join(sorted(flags.split('; ')))


def test_explicit_no_domain_exception_has_unmeasured_coverage(tmp_path):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    saved_pfam(tmp_path, dict(D_gene=[], A_gene=[]))
    selected, audit, pairs, _, _ = filter_events(events, links, tmp_path, allow_both_no_pfam=True)
    assert selected and pairs[0]['coverage_status'] == 'explicit_bilateral_no_hit_exception'
    assert pairs[0]['donor_shared_pfam_query_coverage'] == ''
    assert pairs[0]['donor_shared_pfam_covered_aa'] == ''
    assert audit[0]['pfam_coverage_passing_pair_count'] == 0


def test_coverage_parameter_reaches_trait_tables_without_plotting(source):
    _, links = read_tsv(source[1])
    links = [dict(supported_link(r['event_id'], r['side'], r['gene_id']), **r) for r in links]
    write_tsv(source[1], list(links[0]), links)
    root = source[0].parent / 'families'
    saved_pfam(root, dict(D_gene=['PF01053'], A_gene=['PF01053'], B_gene=['PF01053'], C_gene=['PF01053']))
    report = generate_filtered(*source, plots=False, gene_family_root=root, min_shared_pfam_coverage=.95)
    assert report['min_shared_pfam_coverage'] == .95
    assert report['shared_pfam_filter']['passed_event_count'] == 0
    assert report['shared_pfam_filter']['attention_flags_are_exclusion_criteria'] is False
    assert not read_tsv(source[-1]/'traits/category/all_category1/events.tsv')[1]
    assert all(r['coverage_status'] == 'below_minimum' for r in read_tsv(source[-1]/'pfam_pair_audit.tsv')[1])
    report = generate_filtered(*source, plots=False, gene_family_root=root, min_shared_pfam_coverage=.9)
    assert report['shared_pfam_filter']['passed_event_count'] > 0


@pytest.mark.parametrize('support', ['', '0', '89.9', '90'])
def test_native_focus_does_not_reapply_ufb_filter(support):
    stat, events, links = focused_node_source()
    stat[1]['support_generax_ufboot'] = support
    events[0]['support_used'] = support
    output, audit = annotate(stat, events, links)
    assert audit[0]['status'] == 'selected'
    assert output[1]['hgtfocus_node_label'] == 'HGT1 UFB=' + (support or 'NA')
    assert audit[0]['support_generax_ufboot'] == (float(support) if support else '')


@pytest.mark.parametrize('field', ['branch_id', 'gene_tree_branch_id', 'node_name', 'gene_tree_node'])
def test_pfam_links_cannot_borrow_project_or_native_branch_identity(tmp_path, field):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    links[0][field] = 'wrong'
    with pytest.raises(ValueError, match='identity mismatch'):
        filter_events(events, links, tmp_path)


def test_disabled_pfam_still_records_imported_code(source):
    report = generate(*source, plots=False)
    assert {'focus_hgt_pfam.py', 'focus_hgt_gene_trees.py', 'gene_family_output_store.py'} <= set(report['code_sha256'])


def test_no_ufb_flow_has_analyzed_ogs_zero_step_and_preserves_cohort(tmp_path):
    from focus_hgt_figures import export_filtering_flow
    from pypdf import PdfReader
    events = [dict(event_id='e1', orthogroup='OG1', support_generax_ufboot='12',
                   mapping_status='matched', species_tree_mapping_status='matched_external_tree')]
    audit = tmp_path/'audit.tsv'
    write_tsv(audit, list(events[0]), events)
    output = tmp_path/'plots'
    counts = export_filtering_flow(output, events, events, 'gall', audit,
                                  pfam_selected=events, direction_selected=events,
                                  analyzed_orthogroups=['OG1', 'OG_without_transfers'])
    displayed = read_tsv(output/'filtering_flow.tsv')[1]
    assert displayed[0] == dict(step='00', stage='All analyzed orthogroups', event_count='NA', orthogroup_count='2')
    assert displayed[-1]['event_count'] == '1'
    text = PdfReader(output/'filtering_flow.pdf').pages[0].extract_text()
    assert 'Matched event and species branches' not in text and 'UFB >=90' not in text
    assert 'No UFB threshold' in text
    assert len(counts) == 4


def test_query_taxonomic_ranks_cannot_borrow_a_foreign_saved_lineage(tmp_path):
    import sqlite3

    from focus_hgt_context_annotations import RANKS, ContextAnnotations
    db = tmp_path/'taxa.sqlite'
    with sqlite3.connect(db) as conn:
        conn.execute('CREATE TABLE species(taxid INTEGER PRIMARY KEY,spname TEXT,rank TEXT,track TEXT)')
        conn.execute('CREATE TABLE merged(taxid_old INTEGER,taxid_new INTEGER)')
        conn.executemany('INSERT INTO species VALUES(?,?,?,?)', [
            (1, 'root', 'no rank', '1'), (2, 'Own kingdom', 'kingdom', '2,1'),
            (3, 'Own species', 'species', '3,2,1'), (4, 'Foreign kingdom', 'kingdom', '4,1')])
        conn.execute('INSERT INTO merged VALUES(20,2)')
    annotations = ContextAnnotations(taxonomy_dbfile=db)
    invalid = annotations.query_lineage(3, (1,4,3))
    assert invalid['mmseqs2_lineage_status'] == 'saved_lineage_conflicts_existing_database'
    assert all(invalid['mmseqs2_'+rank] == '' for rank in RANKS)
    merged = annotations.query_lineage(3, (1,20,3))
    assert merged['mmseqs2_kingdom'] == 'Own kingdom'
    missing = annotations.query_lineage(999, (1,4,999))
    assert missing['mmseqs2_lineage_status'] == 'lca_taxid_unresolved_in_existing_database'
    assert all(missing['mmseqs2_'+rank] == '' for rank in RANKS)


def test_missing_context_support_label_is_explicit_na():
    from focus_hgt_context import support_label
    assert support_label({'support_generax_ufboot': ''}) == 'NA'
    assert support_label({'support_generax_ufboot': '0'}) == '0'


def test_no_domain_exception_records_the_pair_that_actually_passed(tmp_path):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    eid = events[0]['event_id']
    links += [supported_link(eid, 'donor', 'D_nohit'), supported_link(eid, 'recipient', 'A_nohit')]
    path = saved_pfam(tmp_path, {'D_gene': ['PF01053'], 'A_gene': ['PF01053'], 'D_nohit': [], 'A_nohit': []})
    fields, hits = read_tsv(path)
    for hit in hits:
        if hit['sacc']:
            hit['qend'] = '49'
    write_tsv(path, fields, hits)
    selected, _, pairs, _, _ = filter_events(events, links, tmp_path, allow_both_no_pfam=True)
    assert selected[0]['pfam_passing_pair_count'] == 1
    assert selected[0]['pfam_best_pair_donor_gene'] == 'D_nohit'
    assert selected[0]['pfam_best_pair_recipient_gene'] == 'A_nohit'
    assert selected[0]['pfam_best_pair_donor_query_coverage'] == ''
    assert selected[0]['pfam_best_pair_recipient_query_coverage'] == ''
    assert next(row for row in pairs if row['passes_pfam_filter'] == 'True')['coverage_status'] == 'explicit_bilateral_no_hit_exception'


@pytest.mark.parametrize('enabled,expected', [('0', False), ('false', False), ('1', True), ('true', True)])
def test_no_domain_option_cannot_use_python_string_truthiness(tmp_path, enabled, expected):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    saved_pfam(tmp_path, {'D_gene': [], 'A_gene': []})
    assert bool(filter_events(events, links, tmp_path, allow_both_no_pfam=enabled)[0]) is expected


@pytest.mark.parametrize('case', ['branch_alias', 'duplicate_link', 'transferred_out', 'duplicate_event'])
def test_direct_native_renderer_cannot_borrow_invalid_event_gene_evidence(case):
    stat, events, links = focused_node_source()
    if case == 'branch_alias':
        links[0]['branch_id'] = 'wrong'
    elif case == 'duplicate_link':
        links.append(dict(links[0]))
    elif case == 'transferred_out':
        links[0]['lineage_status'] = 'transferred_out_of_donor_lineage'
    else:
        events.append(dict(events[0]))
    with pytest.raises(ValueError):
        annotate(stat, events, links)


@pytest.mark.parametrize('accession', ['', 'P_stale'])
def test_distribution_link_hit_cannot_override_an_explicit_native_no_hit(tmp_path, accession):
    from io import StringIO

    from Bio import Phylo
    from focus_hgt_figures import export_figures
    stat, events, links = focused_node_source()
    for row in stat:
        row.update(sprot_best='', sprot_recname='', organism='')
    root = tmp_path/'families'
    write_tsv(root/'stat_branch/OG1_stat.branch.tsv', list(stat[0]), stat)
    for link, species in zip(links, ['D', 'A'], strict=True):
        link['gene_species'] = species
    links[1].update(besthit_accession=accession, swissprot_best_hit_protein_name='Stale predicted product')
    events[0].update(generax_donor_node='D', generax_recipient_node='A')
    tree = Phylo.read(StringIO('(A:1,D:1)root;'), 'newick')
    with pytest.raises(ValueError, match='best hit disagrees'):
        export_figures(tmp_path/'plots', events, events, links, tree, {'A': 1, 'D': 0}, root, 'gall')


def test_focus_api_normalizes_disabled_pfam_flags_before_recording_the_manifest(source):
    report = generate_filtered(*source, plots=False, require_shared_pfam='0', allow_both_no_pfam='false',
                               require_length_ratio='false')
    assert report['require_shared_pfam'] is False
    assert report['allow_both_no_pfam'] is False
    assert 'pfam_pair_filter' not in report['filtering_order']
    assert not (source[-1]/'pfam_pair_audit.tsv').exists()
    assert len(read_tsv(source[-1]/'traits/binary/all_category1/events.tsv')[1]) == 3


@pytest.mark.parametrize('name', ['require_shared_pfam', 'allow_both_no_pfam', 'require_length_ratio'])
@pytest.mark.parametrize('value', [None, '', 'yes', 2, float('nan'), []])
def test_focus_api_rejects_ambiguous_pfam_flags_before_writing(source, name, value):
    with pytest.raises(ValueError, match=name):
        generate_filtered(*source, plots=False, **{name: value})
    assert not source[-1].exists()


def test_direct_native_export_validates_links_before_skipping_a_missing_family(tmp_path):
    from focus_hgt_gene_trees import export_gene_trees
    _, events, links = focused_node_source()
    links[0]['orthogroup'] = 'OG_wrong'
    output = tmp_path/'plots'
    with pytest.raises(ValueError, match='identity mismatch'):
        export_gene_trees(output, events, links, tmp_path/'missing')
    assert not output.exists()


def test_native_annotation_never_reuses_one_family_tree_for_another_family():
    stat, events, links = focused_node_source()
    other = dict(events[0], event_id='OG2:3:1', orthogroup='OG2')
    with pytest.raises(ValueError, match='single orthogroup'):
        annotate(stat, events + [other], links)


def test_distribution_only_annotates_the_selected_events(tmp_path):
    from io import StringIO

    from Bio import Phylo
    from focus_hgt_figures import export_figures
    stat, events, links = focused_node_source()
    root = tmp_path/'families'
    write_tsv(root/'stat_branch/OG1_stat.branch.tsv', list(stat[0]), stat)
    for link, species in zip(links, ['D', 'A'], strict=True):
        link['gene_species'] = species
    events[0].update(generax_donor_node='D', generax_recipient_node='A')
    other = dict(events[0], event_id='OG1:3:2', event_index='2', generax_transfer='Y@D@B')
    links.append(dict(links[1], event_id=other['event_id'], gene_id='A_missing',
                      swissprot_best_hit_protein_name='Unselected product'))
    tree = Phylo.read(StringIO('(A:1,D:1)root;'), 'newick')
    output = tmp_path/'plots'
    export_figures(output, events + [other], events, links, tree, {'A': 1, 'D': 0}, root, 'gall')
    assert all(row['protein_product'] == 'Annotation unavailable'
               for row in read_tsv(output/'orthogroup_species_distribution.tsv')[1])


def test_pfam_family_index_keeps_reused_gene_names_and_event_pairs_separate(tmp_path):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    other = dict(events[0], event_id='OG2:3:1', orthogroup='OG2')
    links += [dict(link, event_id=other['event_id'], orthogroup='OG2') for link in links]
    saved_pfam(tmp_path, {'D_gene': ['PF01053'], 'A_gene': ['PF01053'], 'unused': ['PF00001']})
    saved_pfam(tmp_path, {'D_gene': ['PF01053'], 'A_gene': ['PF00001']}, family='OG2')
    selected, audit, pairs, genes, sources = filter_events(events + [other], links, tmp_path)
    assert [row['event_id'] for row in selected] == [events[0]['event_id']]
    assert [row['pfam_filter_status'] for row in audit] == ['passed', 'withheld']
    assert [row['pair_status'] for row in pairs] == ['shared_pfam_detected', 'detected_pfam_sets_disjoint']
    assert [(row['orthogroup'], row['gene_id'], row['pfam_accessions']) for row in genes] == [
        ('OG1', 'D_gene', 'PF01053'), ('OG1', 'A_gene', 'PF01053'),
        ('OG2', 'D_gene', 'PF01053'), ('OG2', 'A_gene', 'PF00001')]
    assert len(sources) == 2


def test_unused_pfam_rows_still_receive_structural_validation(tmp_path):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    path = saved_pfam(tmp_path, {'D_gene': ['PF01053'], 'A_gene': ['PF01053'], 'unused': ['PF00001']})
    with path.open('a') as handle:
        handle.write('unused\tmalformed\n')
    with pytest.raises(ValueError, match='Malformed saved Pfam RPS-BLAST row'):
        filter_events(events, links, tmp_path)


@pytest.mark.parametrize('donor_length,recipient_length,passed', [(200, 100, True), (100, 200, True),
                                                               (201, 100, False), (100, 201, False)])
def test_protein_length_ratio_is_inclusive_and_default_on(tmp_path, donor_length, recipient_length, passed):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    path = saved_pfam(tmp_path, dict(D_gene=['PF01053'], A_gene=['PF01053']))
    fields, rows = read_tsv(path)
    for row, length in zip(rows, (donor_length, recipient_length), strict=True):
        row.update(qlen=str(length), qend=str(length))  # Both shared-domain coverages are 100%.
    write_tsv(path, fields, rows)
    selected, audit, pairs, _, _ = filter_events(events, links, tmp_path)
    assert bool(selected) is passed
    pair = pairs[0]
    assert pair['protein_length_ratio'] == min(donor_length, recipient_length) / max(donor_length, recipient_length)
    assert pair['passes_shared_pfam_filter'] == 'True'
    assert pair['passes_length_ratio_filter'] == str(passed)
    assert pair['passes_pair_filter'] == pair['passes_pfam_filter'] == str(passed)
    assert audit[0]['length_ratio_filter_enabled'] == 'True'
    assert audit[0]['pfam_filter_reason'] == ('' if passed else 'protein_length_ratio_below_minimum_or_unmeasured')
    assert filter_events(events, links, tmp_path, require_length_ratio='0')[0]


def test_length_and_pfam_cannot_be_satisfied_by_different_pairs(tmp_path):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    eid = events[0]['event_id']
    links.append(supported_link(eid, 'recipient', 'A_other'))
    path = saved_pfam(tmp_path, dict(D_gene=['PF01053'], A_gene=['PF01053'], A_other=['PF00001']))
    fields, rows = read_tsv(path)
    for row in rows:
        length = 201 if row['qacc'] == 'D_gene' else 100 if row['qacc'] == 'A_gene' else 201
        row.update(qlen=str(length), qend=str(length))
    write_tsv(path, fields, rows)
    selected, audit, pairs, _, _ = filter_events(events, links, tmp_path)
    assert not selected
    assert audit[0]['pfam_coverage_passing_pair_count'] == audit[0]['length_ratio_passing_pair_count'] == 1
    assert all(row['passes_pair_filter'] == 'False' for row in pairs)
    assert filter_events(events, links, tmp_path, require_length_ratio=False)[0]
    assert filter_events(events, links, tmp_path, require_shared_pfam=False)[0]


def test_no_pfam_requirement_does_not_disable_length_or_allow_missing_lengths(tmp_path):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    # A recorded no-hit protein provides a length; a missing search does not.
    saved_pfam(tmp_path, dict(D_gene=[], A_gene=[]))
    selected, _, pairs, _, _ = filter_events(events, links, tmp_path, require_shared_pfam=False)
    assert selected and pairs[0]['passes_shared_pfam_filter'] == 'False'
    assert pairs[0]['protein_length_ratio'] == 1
    saved_pfam(tmp_path, dict(D_gene=[]))
    selected, _, pairs, _, _ = filter_events(events, links, tmp_path, require_shared_pfam=False, allow_both_no_pfam=True)
    assert not selected and pairs[0]['protein_length_ratio'] == ''
    assert pairs[0]['length_ratio_status'] == 'unmeasured'


def test_no_domain_exception_cannot_bypass_length_threshold(tmp_path):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    path = saved_pfam(tmp_path, dict(D_gene=[], A_gene=[]))
    fields, rows = read_tsv(path)
    rows[0]['qlen'] = '201'
    write_tsv(path, fields, rows)
    selected, _, pairs, _, _ = filter_events(events, links, tmp_path, allow_both_no_pfam=True)
    assert not selected and pairs[0]['passes_shared_pfam_filter'] == 'True'
    assert filter_events(events, links, tmp_path, allow_both_no_pfam=True, require_length_ratio=False)[0]


@pytest.mark.parametrize('minimum', [-.1, 1.1, float('nan'), float('inf'), None, 'bad'])
def test_invalid_length_ratio_configuration_is_rejected_for_empty_cohort(minimum):
    from focus_hgt_pfam import filter_events
    with pytest.raises(ValueError, match='Protein length ratio must be a finite fraction'):
        filter_events([], [], '', min_length_ratio=minimum)


def test_pair_evidence_cache_keeps_each_event_identity_and_rows_independent(tmp_path):
    from focus_hgt_pfam import filter_events
    _, events, links = focused_node_source()
    other = dict(events[0], event_id='OG1:3:2', event_index='2')
    links += [dict(link, event_id=other['event_id'], event_index='2') for link in list(links)]
    saved_pfam(tmp_path, dict(D_gene=['PF01053'], A_gene=['PF01053']))
    _, _, pairs, _, _ = filter_events(events + [other], links, tmp_path)
    assert [row['event_id'] for row in pairs] == [events[0]['event_id'], other['event_id']]
    assert pairs[0] is not pairs[1]
    pairs[0]['pair_attention_flags'] = 'test_mutation'
    assert pairs[1]['pair_attention_flags'] == ''


def test_length_only_mode_reaches_trait_tables_and_manifest(source):
    _, links = read_tsv(source[1])
    links = [dict(supported_link(r['event_id'], r['side'], r['gene_id']), **r) for r in links]
    write_tsv(source[1], list(links[0]), links)
    root = source[0].parent / 'families'
    saved_pfam(root, dict(D_gene=[], A_gene=[], B_gene=[], C_gene=[]))
    report = generate_filtered(*source, plots=False, gene_family_root=root, require_shared_pfam='0')
    assert report['require_length_ratio'] is True
    assert report['pair_filter_criteria'] == ['protein_length_ratio']
    assert report['pair_filter']['passed_event_count'] == 5
    assert 'protein_length_pair_filter' in report['filtering_order']
    assert all(row['passes_pair_filter'] == 'True' and row['passes_shared_pfam_filter'] == 'False'
               for row in read_tsv(source[-1] / 'pfam_pair_audit.tsv')[1])
    focused = read_tsv(source[-1] / 'traits/binary/all_category1/events.tsv')[1]
    assert len(focused) == 3 and all(row['length_ratio_filter_enabled'] == 'True' for row in focused)


def test_pair_flow_label_records_enabled_length_rule(tmp_path):
    from focus_hgt_figures import export_filtering_flow
    events = [dict(event_id='e1', orthogroup='OG1')]
    export_filtering_flow(tmp_path, events, events, 'gall', pfam_selected=events,
                          pair_filter_label='Event-pair length ratio >=0.5')
    stages = read_tsv(tmp_path / 'filtering_flow.tsv')[1]
    assert stages[-2]['stage'] == 'Event-pair length ratio >=0.5'


def test_origin_source_changed_after_assessment_preserves_previous_bundle(source, monkeypatch):
    import focus_hgt_traits
    _, links = read_tsv(source[1])
    links = [dict(supported_link(r['event_id'], r['side'], r['gene_id']), **r) for r in links]
    write_tsv(source[1], list(links[0]), links)
    root = source[0].parent / 'families'
    saved_pfam(root, dict(D_gene=['PF01053'], A_gene=['PF01053'], B_gene=['PF01053'], C_gene=['PF01053']))
    sequence = source[0].parent / 'mmseqs'
    sequence.mkdir()
    raw = sequence / 'D_mmseqs2taxonomy.tsv'
    raw.write_text('D_gene\t11\tclass\tSaved donor class\t1\t1\t1\t1\t1;2;11\n')
    generate_filtered(*source, plots=False, gene_family_root=root, mmseqs2_taxonomy_dir=sequence)
    previous = (source[-1] / 'manifest.json').read_bytes()
    original = focus_hgt_traits.export_bundle
    def change_after_origin(*args, **kwargs):
        result = original(*args, **kwargs)
        with raw.open('a') as handle:
            handle.write('unused\t11\tclass\tSaved donor class\t1\t1\t1\t1\t1;2;11\n')
        return result
    monkeypatch.setattr(focus_hgt_traits, 'export_bundle', change_after_origin)
    with pytest.raises(ValueError, match='Origin taxonomy input changed'):
        generate_filtered(*source, plots=False, gene_family_root=root, mmseqs2_taxonomy_dir=sequence)
    assert (source[-1] / 'manifest.json').read_bytes() == previous
