import json
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))

from focus_hgt_gene_trees import annotate, background_supported  # noqa: E402
from focus_hgt_traits import generate, read_tsv, write_tsv  # noqa: E402
from species_trait_schema import schema_path, schema_payload  # noqa: E402


def supported_link(event_id, side, gene):
    return dict(event_id=event_id, orthogroup="OG1", side=side, gene_id=gene,
                eligible_for_context="True", host_scaffold_status="measured", host_scaffold_id="scaffold1",
                host_scaffold_background_class_total_count="20", host_scaffold_background_class_compatible_count="9",
                host_scaffold_background_class_incompatible_count="1", host_scaffold_background_class_unresolved_count="10",
                host_scaffold_background_class_classified_fraction="0.5",
                host_scaffold_background_class_compatible_fraction="0.9")


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
    output, audit = annotate(stat, events, links)
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
    assert 'Kingdom: unavailable' in cells[3]  # No name-based taxonomy inference.
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
    assert 'best-hit prediction' in text and 'GFF product' not in text
    assert 'Transcriptional regulator (NtrC/NifA family)' not in text
    assert 'Schizosaccharomycetaceae' in text and 'Genus: Schizosaccharomyces' in text
    assert len(annotations.display_audit) == 42
    assert {(r['context_focal_gene_id'], r['gene_id']) for r in annotations.display_audit} == expected
    assert all(r['event_ids'] == events[0]['event_id'] for r in annotations.display_audit)
    assert all(r['orthogroup'] == 'OtherOG' for r in annotations.display_audit if r['context_role'] == 'neighbor')
    for record in records:
        cells = annotation_cells(record, record['gene_id'][0], 'Focal')
        for i, cell in enumerate(cells):
            assert all(text_width(line) <= (TABLE_EDGES[i+1]-TABLE_EDGES[i])*TABLE_WIDTH_PT-10
                       for line in cell.split('\n'))


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

    row = dict(event_id='e1', orthogroup='OG1', gene_tree_branch_id='3')
    with pytest.raises(ValueError, match='Duplicate'):
        filtering_counts([row, row], [row])
    with pytest.raises(ValueError, match='not a subset'):
        filtering_counts([row], [dict(row, event_id='e2')])
    for field in ('orthogroup', 'gene_tree_branch_id'):
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
