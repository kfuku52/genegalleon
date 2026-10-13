import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))

from focus_hgt_origin import (  # noqa: E402
    OriginEvidence,
    annotate_origin,
    auxiliary_status,
    endpoint_relation,
    taxid_values,
)
from species_taxonomy import pack_values  # noqa: E402

D = 'Plant_host_many_name_parts'
INSECT = 'Insect_host_many_name_parts'
P = D + '_subsample'
K = 'Bacterium_host'


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


@pytest.fixture
def taxonomy(tmp_path):
    path = tmp_path / 'species_taxonomy.tsv'
    fields = ['species', 'resolution_status', 'domain', 'kingdom', 'phylum', 'class',
              'domain_taxid', 'kingdom_taxid', 'phylum_taxid', 'class_taxid']
    rows = [dict(zip(fields, values, strict=True)) for values in (
        (D, 'resolved', 'Eukaryota', 'Plantae', 'PlantPhylum', 'PlantClass', '2', '10', '110', '11'),
        (INSECT, 'resolved', 'Eukaryota', 'Metazoa', 'Arthropoda', 'Insecta', '2', '20', '210', '21'),
        (P, 'resolved', 'Eukaryota', 'Plantae', 'PlantPhylum', 'OtherPlantClass', '2', '10', '110', '13'),
        (K, 'resolved', 'Bacteria', 'BacterialKingdom', 'BacterialPhylum', 'BacterialClass', '1', '40', '410', '41'),
    )]
    write_tsv(path, fields, rows)
    return path


def mmseq(root, species, rows):
    path = root / (species + '_mmseqs2taxonomy.tsv')
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        for gene, assigned, rank, lineage in rows:
            writer.writerow([gene, assigned, rank, 'Saved LCA organism', '4', '4', '4', '1.0', lineage])
    return path


def scaffold(root, species, gene, labels):
    fields = ['species', 'gene_id', 'scaffold', 'locus_id', 'count_unit', 'rank', 'host_taxid', 'label']
    path = root / (species + '_gene_taxonomy.tsv')
    rows = [dict(zip(fields, (species, gene, 'scaffold1', 'locus1', 'gff_locus', rank, host, label), strict=True))
            for rank, host, label in [
                ('domain', '2', labels.get('domain', 'compatible')), ('phylum', '110', 'unresolved'),
                ('class', '11', labels.get('class', 'compatible')), ('order', '', 'unresolved'),
                ('family', '', 'unresolved'), ('genus', '', 'unresolved'), ('species', '', 'unresolved')]]
    write_tsv(path, fields, rows)
    return path


def event(identity='OG1:4:1', donor=D, recipient=INSECT, **extra):
    return dict(event_id=identity, orthogroup='OG1', event_index='1', generax_donor_node=donor,
                generax_recipient_node=recipient, generax_transfer=f'Y@{donor}@{recipient}', **extra)


def gene_link(e, side, gene, species):
    return dict(event_id=e['event_id'], orthogroup=e['orthogroup'], side=side, gene_id=gene,
                gene_species=species, lineage_status='retained', expression_measured='False',
                intron_supported='', synteny_support_score='')


def gene_rows(links):
    return [{field: row[field] for field in ('event_id', 'orthogroup', 'side', 'gene_id')} for row in links]


def pair(e, donor, recipient, passing=True, **extra):
    return dict(e, donor_gene_id=donor, recipient_gene_id=recipient,
                passes_pfam_filter=str(passing), **extra)


def run(e, links, pairs, taxonomy, tmp_path, **options):
    nodes = {D: (D,), INSECT: (INSECT,), P: (P,), K: (K,), 'n43': (D, INSECT, P), 'n6': (D, P)}
    return annotate_origin([e], pairs, gene_rows(links), nodes, links=links,
                           mmseqs2_taxonomy_dir=tmp_path / 'mmseq', species_taxonomy=taxonomy, **options)


@pytest.mark.parametrize('assigned,rank,lineage,expected,highest', [
    ('112', 'species', '100;2;10;110;11;112', 'compatible', ''),
    ('132', 'species', '100;2;10;110;13;132', 'incompatible', 'class'),
    ('212', 'species', '100;2;20;210;21;212', 'incompatible', 'kingdom'),
    ('412', 'species', '100;1;40;410;41;412', 'incompatible', 'domain'),
    ('110', 'phylum', '100;2;10;110', 'unresolved', ''),
    ('210', 'phylum', '100;2;20;210', 'unresolved', 'kingdom'),
    ('210', 'no rank', '100;2;20;210', 'unresolved', 'kingdom'),
    ('0', 'no rank', '', 'unresolved', ''),
    ('112', 'species', '', 'unresolved', ''),
])
def test_focal_classification_resolved_broad_and_unclassified_are_distinct(taxonomy, tmp_path, assigned, rank, lineage, expected, highest):
    e = event()
    links = [gene_link(e, 'donor', 'opaque_gene_id', D), gene_link(e, 'recipient', 'recipient', INSECT)]
    source = mmseq(tmp_path / 'mmseq', D, [('opaque_gene_id', assigned, rank, lineage)])
    mmseq(tmp_path / 'mmseq', INSECT, [('recipient', '212', 'species', '100;2;20;210;21;212')])
    before = source.read_bytes()
    events, pairs, genes, report = run(e, links, [pair(e, 'opaque_gene_id', 'recipient')], taxonomy, tmp_path)
    assert genes[0]['origin_gene_species'] == D  # The ID itself cannot establish species.
    assert genes[0]['origin_host_class_status'] == expected
    assert genes[0]['origin_highest_mismatch_rank'] == highest
    assert pairs[0]['donor_origin_host_class_status'] == expected
    assert events[0][f'origin_donor_class_{expected}_pair_count'] == 1
    assert report['source_snapshot_verified'] is True
    assert str(source.resolve()) in report['source_sha256'] and source.read_bytes() == before
    assert genes[0]['origin_expression_evidence_status'] == 'not_measured'
    assert genes[0]['origin_intron_evidence_status'] == genes[0]['origin_synteny_evidence_status'] == 'unavailable'


def test_pair_flags_use_same_passing_pair_and_do_not_borrow_a_compatible_copy(taxonomy, tmp_path):
    e = event()
    links = [gene_link(e, 'donor', 'good', D), gene_link(e, 'donor', 'bad', D), gene_link(e, 'recipient', 'r', INSECT)]
    mmseq(tmp_path / 'mmseq', D, [('good', '112', 'species', '100;2;10;110;11;112'),
                                 ('bad', '212', 'species', '100;2;20;210;21;212')])
    mmseq(tmp_path / 'mmseq', INSECT, [('r', '212', 'species', '100;2;20;210;21;212')])
    pairs = [pair(e, 'bad', 'r'), pair(e, 'good', 'r', True, passes_pair_filter='False')]
    events, rows, _, _ = run(e, links, pairs, taxonomy, tmp_path)
    assert events[0]['origin_passing_pair_count'] == 1
    assert events[0]['origin_donor_class_incompatible_pair_count'] == 1
    assert events[0]['origin_donor_class_compatible_pair_count'] == 0
    assert 'all_passing_donor_pairs_class_incompatible' in events[0]['origin_review_flags']
    assert rows[1]['donor_origin_host_class_status'] == 'compatible'
    assert rows[0]['donor_origin_host_class_status'] == 'incompatible'
    e2 = event()
    events, _, _, _ = run(e2, links, [pair(e2, 'bad', 'r'), pair(e2, 'good', 'r')], taxonomy, tmp_path)
    assert 'some_passing_donor_pairs_class_incompatible' in events[0]['origin_review_flags']
    assert 'all_passing_donor_pairs_class_incompatible' not in events[0]['origin_review_flags']
    assert events[0]['origin_passing_pair_count'] == 2


def test_exact_species_identifiers_do_not_alias_shared_underscore_prefixes(taxonomy, tmp_path):
    e = event(donor='n6')
    links = [gene_link(e, 'donor', D + '_subsample_gene', P), gene_link(e, 'recipient', 'r', INSECT)]
    mmseq(tmp_path / 'mmseq', D, [(D + '_subsample_gene', '112', 'species', '100;2;10;110;11;112')])
    mmseq(tmp_path / 'mmseq', P, [(D + '_subsample_gene', '132', 'species', '100;2;10;110;13;132')])
    mmseq(tmp_path / 'mmseq', INSECT, [('r', '212', 'species', '100;2;20;210;21;212')])
    _, _, genes, report = run(e, links, [pair(e, D + '_subsample_gene', 'r')], taxonomy, tmp_path)
    assert genes[0]['origin_gene_species'] == P and genes[0]['origin_host_class_taxid'] == '13'
    assert str((tmp_path / 'mmseq' / (D + '_mmseqs2taxonomy.tsv')).resolve()) not in report['source_sha256']


def test_missing_file_record_species_and_unclassified_remain_explicit(taxonomy, tmp_path):
    e = event()
    links = [gene_link(e, 'donor', 'missing', D), gene_link(e, 'recipient', 'r', INSECT)]
    mmseq(tmp_path / 'mmseq', INSECT, [('other', '212', 'species', '100;2;20;210;21;212')])
    events, _, genes, report = run(e, links, [pair(e, 'missing', 'r')], taxonomy, tmp_path)
    assert genes[0]['origin_mmseqs2_status'] == 'source_file_unavailable'
    assert genes[1]['origin_mmseqs2_status'] == 'gene_record_unavailable'
    assert events[0]['origin_donor_class_missing_pair_count'] == 1
    assert events[0]['origin_recipient_class_missing_pair_count'] == 1
    assert len(report['missing_sources']) == 1
    del links[0]['gene_species']
    events, _, genes, _ = run(event(), links, [pair(e, 'missing', 'r')], taxonomy, tmp_path)
    assert genes[0]['origin_gene_species_status'] == 'gene_species_unavailable'
    assert genes[0]['origin_host_class_status'] == 'missing'


def test_saved_scaffold_label_conflicts_are_flagged_not_promoted_or_excluded(taxonomy, tmp_path):
    e = event()
    links = [gene_link(e, 'donor', 'd', D), gene_link(e, 'recipient', 'r', INSECT)]
    mmseq(tmp_path / 'mmseq', D, [('d', '212', 'species', '100;2;20;210;21;212')])
    mmseq(tmp_path / 'mmseq', INSECT, [('r', '212', 'species', '100;2;20;210;21;212')])
    path = scaffold(tmp_path / 'scaffold', D, 'd', {'class': 'compatible'})
    events, _, genes, report = run(e, links, [pair(e, 'd', 'r')], taxonomy, tmp_path,
                                 scaffold_taxonomy_dir=path.parent)
    assert len(events) == 1
    assert genes[0]['origin_host_class_status'] == 'conflicting'
    assert events[0]['origin_donor_class_conflicting_pair_count'] == 1
    assert 'classification_sources_conflict' in events[0]['origin_review_flags']
    assert report['attention_flags_are_exclusion_criteria'] is False


@pytest.mark.parametrize('field,value,status', [
    ('expression_measured', None, 'unavailable'), ('expression_measured', '', 'unavailable'),
    ('expression_measured', 'False', 'not_measured'), ('expression_measured', False, 'not_measured'),
    ('expression_measured', 'True', 'measured'),
    ('intron_supported', '', 'unavailable'), ('intron_supported', 'False', 'recorded_support_flag'),
    ('intron_supported', True, 'recorded_support_flag'), ('synteny_support_score', '', 'unavailable'),
    ('synteny_support_score', '0', 'recorded_score'), ('synteny_support_score', '0.5', 'recorded_score'),
])
def test_auxiliary_missingness_never_becomes_measured_zero(field, value, status):
    assert auxiliary_status({field: value}, field) == status
    if field == 'intron_supported' and value == 'False':
        assert auxiliary_status({field: value, 'num_intron': '0'}, field) == 'observed_count_available'
        assert auxiliary_status({field: value, 'num_intron': '-999'}, field) == 'recorded_support_flag'


def test_focal_lca_above_class_cannot_inherit_resolved_saved_class_label(taxonomy, tmp_path):
    e = event()
    links = [gene_link(e, 'donor', 'd', D), gene_link(e, 'recipient', 'r', INSECT)]
    mmseq(tmp_path / 'mmseq', D, [('d', '110', 'phylum', '100;2;10;110')])
    path = scaffold(tmp_path / 'scaffold', D, 'd', {'class': 'compatible'})
    _, _, genes, _ = run(e, links, [pair(e, 'd', 'r')], taxonomy, tmp_path, scaffold_taxonomy_dir=path.parent)
    assert genes[0]['origin_host_class_status'] == 'conflicting'


def test_scope_rejects_ambiguous_gene_species_and_a_late_appearing_missing_source(taxonomy, tmp_path):
    e = event()
    links = [gene_link(e, 'donor', 'd', D), gene_link(e, 'recipient', 'r', INSECT)]
    with pytest.raises(ValueError, match='Ambiguous'):
        OriginEvidence(links + [gene_link(e, 'donor', 'd', P)])
    scope = OriginEvidence(links, tmp_path / 'mmseq', species_taxonomy=taxonomy)
    assert scope.gene('d', D)['origin_host_class_status'] == 'missing'
    mmseq(tmp_path / 'mmseq', D, [('d', '112', 'species', '100;2;10;110;11;112')])
    with pytest.raises(ValueError, match='appeared'):
        scope.finalize_report()


@pytest.mark.parametrize('donor,recipient,expected', [
    (('A', 'B', 'C'), ('A', 'B'), 'donor_ancestor_of_recipient'),
    (('A', 'B'), ('A', 'B', 'C'), 'recipient_ancestor_of_donor'),
    (('A',), ('B',), 'incomparable_species_clades'),
    (('A', 'B'), ('A', 'B'), 'same_species_clade'),
    (('A', 'B'), ('B', 'C'), 'overlapping_species_clades'),
    (None, ('B',), 'species_branch_unmapped'),
])
def test_endpoint_relation_tracks_internal_clades(donor, recipient, expected):
    assert endpoint_relation(donor, recipient) == expected


def test_ancestor_and_mixed_unknown_branches_are_review_only_with_no_pair_filter(taxonomy, tmp_path):
    e = event(donor='n43', recipient='n6')
    events, pairs, genes, report = run(e, [], [], taxonomy, tmp_path)
    assert events[0]['origin_endpoint_relation'] == 'donor_ancestor_of_recipient'
    assert events[0]['origin_donor_taxonomy_status'] == 'mixed_kingdom_clade'
    assert events[0]['origin_recipient_taxonomy_status'] == 'mixed_class_clade'
    assert 'donor_ancestor_of_recipient_origin_review' in events[0]['origin_review_flags']
    assert len(events) == 1 and not pairs and not genes and report['passing_pair_count'] == 0
    unknown = event(donor='missing')
    events, _, _, _ = run(unknown, [], [], taxonomy, tmp_path)
    assert events[0]['origin_donor_taxonomy_status'] == 'species_branch_unmapped'


def test_internal_branch_keeps_partial_classification_unknown_without_requiring_absent_ranks(taxonomy, tmp_path):
    with taxonomy.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields, rows = reader.fieldnames, list(reader)
    for row in rows:
        if row['species'] == P:
            row.update({'class': '', 'class_taxid': ''})
        if row['species'] == K:
            row.update({'kingdom': '', 'kingdom_taxid': '', 'class': '', 'class_taxid': ''})
    write_tsv(taxonomy, fields, rows)
    scope = OriginEvidence(species_taxonomy=taxonomy)
    assert scope.branch_taxonomy((D, P)) == 'partially_unresolved_class_clade'
    assert scope.branch_taxonomy((K,)) == 'resolved_uniform_clade'
    e = event(donor='n6')
    links = [gene_link(e, 'donor', 'd', D), gene_link(e, 'recipient', 'r', INSECT)]
    mmseq(tmp_path / 'mmseq', D, [('d', '112', 'species', '100;2;10;110;11;112')])
    mmseq(tmp_path / 'mmseq', INSECT, [('r', '212', 'species', '100;2;20;210;21;212')])
    events, _, _, _ = run(e, links, [pair(e, 'd', 'r')], taxonomy, tmp_path)
    assert events[0]['origin_donor_taxonomy_status'] == 'partially_unresolved_class_clade'
    assert events[0]['origin_review_status'] == 'review_required'


@pytest.mark.parametrize('value,expected', [
    ('2759', (2759,)), ('2759.0', (2759,)), ('', ()), ('[]', ()),
    ('["1", "227290"]', (1, 227290)), ('[1, 227290]', (1, 227290)),
])
def test_taxid_cells_accept_native_singleton_and_repeated_rank_formats(value, expected):
    assert taxid_values(value) == expected


@pytest.mark.parametrize('value', ['["1", "1.5"]', '["1", null]', '[true, 2]', '[0, 1]',
                                  '[[-1], 2]', '["1"', '["1", ""]', '{"taxid": 1}'])
def test_malformed_packed_taxids_are_not_silently_dropped(value):
    with pytest.raises(ValueError):
        taxid_values(value)


def test_native_repeated_no_rank_taxids_preserve_all_host_ancestors(taxonomy, tmp_path):
    with taxonomy.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields, rows = reader.fieldnames, list(reader)
    fields = fields + ['no rank', 'no rank_taxid']
    for row in rows:
        row.update({'no rank': pack_values(('root', 'Agrobacterium/Rhizobium group')),
                    'no rank_taxid': pack_values(('1', '227290'))})
    write_tsv(taxonomy, fields, rows)
    e = event()
    links = [gene_link(e, 'donor', 'd', D), gene_link(e, 'recipient', 'r', INSECT)]
    mmseq(tmp_path / 'mmseq', D, [('d', '112', 'species', '100;2;10;110;11;112')])
    events, _, genes, _ = run(e, links, [pair(e, 'd', 'r')], taxonomy, tmp_path)
    assert genes[0]['origin_host_class_status'] == 'compatible'
    assert events[0]['origin_donor_taxonomy_status'] == 'resolved_uniform_clade'
    scope = OriginEvidence(species_taxonomy=taxonomy)
    assert {1, 227290}.issubset(scope._host(D)[0])


@pytest.mark.parametrize('values,names,expected', [
    (('11', '13'), ('PlantClass', 'OtherPlantClass'), 'unresolved'),
    (('11', '11'), ('PlantClass', 'PlantClass'), 'compatible'),
    (('11',), ('PlantClass', 'OtherPlantClass'), 'unresolved'),
])
def test_packed_target_ranks_never_choose_an_arbitrary_host_rank(taxonomy, tmp_path, values, names, expected):
    with taxonomy.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields, rows = reader.fieldnames, list(reader)
    for row in rows:
        if row['species'] == D:
            row.update({'class': pack_values(names), 'class_taxid': pack_values(values)})
    write_tsv(taxonomy, fields, rows)
    e = event()
    links = [gene_link(e, 'donor', 'd', D), gene_link(e, 'recipient', 'r', INSECT)]
    mmseq(tmp_path / 'mmseq', D, [('d', '112', 'species', '100;2;10;110;11;112')])
    host_source = scaffold(tmp_path / 'scaffold', D, 'd', {'class': 'compatible'})
    events, _, genes, _ = run(e, links, [pair(e, 'd', 'r')], taxonomy, tmp_path, scaffold_taxonomy_dir=host_source.parent)
    assert genes[0]['origin_host_class_status'] == expected
    assert genes[0]['origin_host_class_taxid_candidates'] == ';'.join(dict.fromkeys(values))
    if expected == 'unresolved':
        assert genes[0]['origin_host_class_taxid'] == ''
        assert 'host_class_rank_ambiguous' in genes[0]['origin_gene_attention_flags']
        assert events[0]['origin_donor_taxonomy_status'] == 'ambiguous_class_clade'
        assert events[0]['origin_review_status'] == 'review_required'
    else:
        assert genes[0]['origin_host_class_taxid'] == '11'
        assert events[0]['origin_donor_taxonomy_status'] == 'resolved_uniform_clade'


def test_malformed_non_target_repeated_rank_taxids_fail_when_consumed(taxonomy, tmp_path):
    with taxonomy.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields, rows = reader.fieldnames, list(reader)
    for row in rows:
        row['no rank_taxid'] = '["1", "bad"]'
    write_tsv(taxonomy, fields + ['no rank_taxid'], rows)
    e = event()
    links = [gene_link(e, 'donor', 'd', D), gene_link(e, 'recipient', 'r', INSECT)]
    with pytest.raises(ValueError, match='taxid'):
        run(e, links, [pair(e, 'd', 'r')], taxonomy, tmp_path)


def test_reusable_scope_reads_each_species_once_and_checks_source_changes(taxonomy, tmp_path, monkeypatch):
    first, second = event(), event('OG1:5:1')
    links = [gene_link(e, side, gene, species) for e in (first, second)
             for side, gene, species in [('donor', 'd', D), ('recipient', 'r', INSECT)]]
    path = mmseq(tmp_path / 'mmseq', D, [('d', '112', 'species', '100;2;10;110;11;112')])
    mmseq(tmp_path / 'mmseq', INSECT, [('r', '212', 'species', '100;2;20;210;21;212')])
    scope = OriginEvidence(links, tmp_path / 'mmseq', species_taxonomy=taxonomy)
    loaded, original = [], scope._rows
    def count(source):
        loaded.append(source)
        yield from original(source)
    monkeypatch.setattr(scope, '_rows', count)
    for e in (first, second):
        selected_links = [row for row in links if row['event_id'] == e['event_id']]
        _, _, _, report = run(e, selected_links, [pair(e, 'd', 'r')], taxonomy, tmp_path, evidence=scope)
        assert report['source_snapshot_verified'] is False
    assert len(loaded) == 2 and scope.finalize_report()['loaded_species_count'] == 2
    path.write_text(path.read_text() + '\n')
    with pytest.raises(ValueError, match='changed during generation'):
        scope.finalize_report()


@pytest.mark.parametrize('problem', ['duplicate_link', 'other_event_pair', 'other_gene_pair', 'transferred_out', 'wrong_species', 'invalid_decision'])
def test_broken_exact_event_pair_identities_are_rejected(taxonomy, tmp_path, problem):
    e = event()
    links = [gene_link(e, 'donor', 'd', D), gene_link(e, 'recipient', 'r', INSECT)]
    pairs = [pair(e, 'd', 'r')]
    if problem == 'duplicate_link':
        links.append(dict(links[0]))
    elif problem == 'other_event_pair':
        pairs[0]['event_id'] = 'other'
    elif problem == 'other_gene_pair':
        pairs[0]['donor_gene_id'] = 'other'
    elif problem == 'transferred_out':
        links[0]['lineage_status'] = 'transferred_out'
    elif problem == 'wrong_species':
        links[0]['gene_species'] = P
    else:
        pairs[0]['passes_pfam_filter'] = ''
    with pytest.raises(ValueError):
        run(e, links, pairs, taxonomy, tmp_path)


def test_malformed_requested_saved_lineage_and_scaffold_ranks_fail(taxonomy, tmp_path):
    e = event()
    links = [gene_link(e, 'donor', 'd', D), gene_link(e, 'recipient', 'r', INSECT)]
    mmseq(tmp_path / 'mmseq', D, [('d', '112', 'species', '100;2;11;112;112')])
    with pytest.raises(ValueError, match='lineage'):
        run(e, links, [pair(e, 'd', 'r')], taxonomy, tmp_path)
    mmseq(tmp_path / 'mmseq', D, [('d', '112', 'species', '100;2;10;110;11;112')])
    host_path = scaffold(tmp_path / 'scaffold', D, 'd', {})
    host_path.write_text('\n'.join(host_path.read_text().splitlines()[:-1]) + '\n')
    with pytest.raises(ValueError, match='Incomplete'):
        run(e, links, [pair(e, 'd', 'r')], taxonomy, tmp_path, scaffold_taxonomy_dir=host_path.parent)
