"""Container-only rendering and gene-summary integration for saved-family plots."""
import os
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas
import pytest

from workflow.support.gene_family_output_store import archive_completed_outputs, family_context
from workflow.tests.test_presence_absence_inputs import saved_manifest

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow/support"


def run_r(args, directory):
    assert shutil.which("Rscript"), "Run this test in a GeneGalleon container"
    return subprocess.run(["Rscript", str(SUPPORT / "plot_query2family_presence_absence.R"), *args],
                          cwd=directory, capture_output=True, text=True, check=True)


def summary_fixture(tmp_path):
    tree = tmp_path / "species.nwk"
    tree.write_text("(Species_one:1,Species_two:1);\n")
    families = [f"family_{i}" for i in range(70)]
    long_table = tmp_path / "long.tsv"
    pandas.DataFrame([
        dict(species=species, species_display=species.replace("_", " "), query=family,
             query_order=i, presence=1, copy_number=1, status="complete")
        for i, family in enumerate(families, 1) for species in ("Species_one", "Species_two")
    ]).to_csv(long_table, sep="\t", index=False)
    return [f"--species_tree={tree}", f"--long_table={long_table}"]


def test_empty_outputs_auto_width_label_map_and_focus_rows(tmp_path):
    args = summary_fixture(tmp_path)
    labels = tmp_path / "labels.tsv"
    display_name = "Readable_family_name_with_wide_WWWWWWWWWWWWWWWWWWWWW"
    labels.write_text(f"kind\tid\tlabel\nfamily\tfamily_0\t{display_name}\n")
    svg = tmp_path / "summary.svg"
    run_r([*args, "--width=auto", f"--label_map={labels}", "--focus_species=Species_one",
           "--out_pdf=", f"--out_svg={svg}"], tmp_path)
    assert not (tmp_path / "1").exists()
    root = ET.parse(svg).getroot()
    width, height = [float(value) for value in root.attrib["viewBox"].split()[2:]]
    assert width == pytest.approx((4.8 + 70 * 0.14) * 72, abs=0.02)
    text = svg.read_text()
    assert display_name in text
    assert ">family_0</text>" not in text
    label_boxes = [(float(y), float(length)) for y, length in re.findall(
        r"translate\([0-9.]+,([0-9.]+)\) rotate\(-90\).*?textLength='([0-9.]+)px'", text)]
    assert len(label_boxes) == 70
    assert max(y + length for y, length in label_boxes) < height - 3
    # One row outline, distinct from the matrix tiles and their white borders.
    assert any(element.tag.endswith("rect") and "stroke-width: 0.64" in element.attrib.get("style", "")
               and "stroke: none" not in element.attrib.get("style", "")
               and float(element.attrib.get("width", "0")) > width * 0.8 for element in root.iter())


def test_empty_svg_and_both_empty_outputs_do_not_write_sentinel_file(tmp_path):
    args = summary_fixture(tmp_path)
    pdf = tmp_path / "summary.pdf"
    run_r([*args, f"--out_pdf={pdf}", "--out_svg="], tmp_path)
    assert pdf.read_bytes().startswith(b"%PDF")
    result = run_r([*args, "--out_pdf=", "--out_svg="], tmp_path)
    assert "Wrote PDF" not in result.stdout and "Wrote SVG" not in result.stdout
    assert not (tmp_path / "1").exists()


@pytest.mark.parametrize("source", ["query2family", "orthogroup"])
@pytest.mark.parametrize("cutoff", [0, 1])
def test_gene_summary_manifest_only_reuses_saved_trees_and_evidence(tmp_path, source, cutoff):
    assert shutil.which("Rscript"), "Run this test in a GeneGalleon container"
    store, manifest = saved_manifest(tmp_path)
    tree = tmp_path / "species.nwk"
    tree.write_text("((Anchor_one:1,Anchor_two:1):1,Missing_species:2);\n")
    workspace = tmp_path / "workspace"
    output = tmp_path / "summary"
    env = os.environ.copy()
    env.update(gg_support_dir=str(SUPPORT), gg_workspace_dir=str(workspace),
               gene_family_source=source, summary_output_dir=str(output),
               run_species_taxonomy="0", run_family_completion_summary="0",
               run_presence_absence_summary="1", presence_absence_species_tree=str(tree),
               presence_absence_species_tree_ci="", presence_absence_species_tree_support="",
               presence_absence_busco_table="", presence_absence_ortholog_basis="query_gene",
               presence_absence_family_manifest=str(manifest), presence_absence_query_label="label",
               presence_absence_dup_conf_score_threshold=str(cutoff),
               presence_absence_focus_species="Anchor_one", presence_absence_plot_width="auto",
               presence_absence_legend_columns="3" if source == "orthogroup" else "auto")
    annotate_loss = source == "query2family" and cutoff == 0
    if annotate_loss:
        (tmp_path / "losses.tsv").write_text(
            "family_id\tspecies\tanchor_cds_fasta_id\nCombined\tMissing_species\tAnchor_one_a\n"
        )
        env["presence_absence_loss_table"] = "losses.tsv"
    if source == "orthogroup":
        def clade(name, event="", species="", children=""):
            events = f'<eventsRec><{event} speciesLocation="{species}"/></eventsRec>' if event else ""
            return f"<clade><name>{name}</name>{events}{children}</clade>"

        species = clade("root", children=clade("Missing_species") + clade(
            "inner", children=clade("Anchor_one") + clade("Anchor_two")))
        gene = clade("dup", "duplication", "root", clade(
            "owner", "speciation", "root", clade("loss", "loss", "Missing_species") + clade(
                "a", "speciation", "inner", clade("Anchor_one_a", "leaf", "Anchor_one")
                + clade("Anchor_two_a", "leaf", "Anchor_two"))) + clade(
            "b", "speciation", "inner", clade("Anchor_one_b", "leaf", "Anchor_one")
            + clade("Anchor_two_b", "leaf", "Anchor_two")))
        (store / "generax_xml").mkdir()
        (store / "generax_xml/OG1_generax.xml").write_text(
            f"<recPhylo><spTree><phylogeny>{species}</phylogeny></spTree>"
            f"<recGeneTree><phylogeny>{gene}</phylogeny></recGeneTree></recPhylo>"
        )
        query_dir = tmp_path / "queries"
        query_dir.mkdir()
        (query_dir / "OG1").write_text(">near\nAA\n")
        for subdir, suffix in (("stat_tree", "_stat.tree.tsv"), ("tree_plot", "_tree_plot.pdf")):
            folder = store / subdir
            folder.mkdir()
            (folder / f"OG1{suffix}").write_text("fixture\n")
        families, from_name = family_context("query2family", query_dir=query_dir)
        archive_completed_outputs(store, "query2family", families, from_name)
        assert not (store / "generax_xml").exists()
        if cutoff:
            env["presence_absence_mark_terminal_losses"] = "0"
    result = subprocess.run(["bash", str(ROOT / "workflow/core/gg_gene_summary_core.sh")],
                            env=env, cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    prefix = output / "query2family_query_gene_orthologs"
    assert Path(str(prefix) + ".pdf").read_bytes().startswith(b"%PDF")
    assert ">HOG1 Anchor_one_a</text>" in Path(str(prefix) + ".svg").read_text()
    columns = pandas.read_csv(str(prefix) + ".columns.tsv", sep="\t")
    assert len(columns) == 2 and set(columns.family_id) == {"Combined"}
    mapping = pandas.read_csv(str(prefix) + ".query_map.tsv", sep="\t")
    assert set(mapping.hog_ids) == {"HOG1", "HOG2"}
    assert Path(str(prefix) + ".selection.tsv").is_file()
    candidates = pandas.read_csv(str(prefix) + ".dup_conf.tsv", sep="\t")
    if cutoff:
        assert len(candidates) == 2 and set(candidates.dup_conf_score) == {1}
        assert set(candidates.mrca_branch_id) == {0}
        assert "Additional ortholog candidate" in Path(str(prefix) + ".svg").read_text()
    else:
        assert candidates.empty
    assert Path(str(prefix) + ".overlap.tsv").is_file()
    assert "Skipping gene-family presence/absence summary because no" not in result.stdout
    svg_path = Path(str(prefix) + ".svg")
    text = svg_path.read_text()
    if annotate_loss or (source == "orthogroup" and cutoff == 0):
        assert text.count(">0†</text>") == 2  # Matrix and legend.
    else:
        assert "†" not in text
    auto_table = Path(str(prefix) + ".terminal_losses.tsv")
    if source == "orthogroup" and cutoff == 0:
        assert auto_table.read_text().splitlines()[1:] == ["Combined\tMissing_species\tAnchor_one_a"]
    elif source == "query2family" and cutoff:
        assert len(auto_table.read_text().splitlines()) == 1  # Missing XML.
    else:
        assert not auto_table.exists()
    label_boxes = [(float(y), float(length)) for y, length in re.findall(
        r"translate\([0-9.]+,([0-9.]+)\) rotate\(-90\).*?textLength='([0-9.]+)px'", text)]
    svg = ET.parse(svg_path).getroot()
    legend = next(element for element in svg.iter()
                  if element.tag.endswith("text") and element.text.startswith("Query-gene orthologs"))
    assert max(y + length for y, length in label_boxes) + 6 < float(legend.attrib["y"])
    if source == "orthogroup":
        legend_text = [element for element in svg.iter()
                       if element.tag.endswith("text") and element.text in (
                           "undetected", "query-anchor-specific", "shared across query anchors", "D#: mapped duplication")]
        assert len({element.attrib["x"] for element in legend_text}) > 1
        for left in legend_text:
            for right in legend_text:
                if left.attrib["y"] == right.attrib["y"] and float(left.attrib["x"]) < float(right.attrib["x"]):
                    assert float(left.attrib["x"]) + float(left.attrib["textLength"].removesuffix("px")) + 6 < float(right.attrib["x"])


def canonical_svg(path):
    text = Path(path).read_text()
    identifiers = {value: f'id{index}' for index, value in enumerate(re.findall(r"\bid=['\"]([^'\"]+)", text))}
    for value, replacement in sorted(identifiers.items(), key=lambda item: -len(item[0])):
        text = text.replace(value, replacement)
    return ET.tostring(ET.fromstring(text))


def test_related_plot_batch_matches_separate_views(tmp_path):
    busco = tmp_path / 'busco'
    busco.mkdir()
    for species in ('Species_one', 'Species_two'):
        (busco / f'{species}.tsv').write_text('# Busco id\tStatus\tSequence\nBUSCO1\tComplete\tgene1\n')
    args = [*summary_fixture(tmp_path), '--width=auto', f'--busco_table={busco}']
    jobs = []
    for index, mode in enumerate(('presence', 'copy_number', 'presence')):
        single = tmp_path / f'single{index}.svg'
        run_r([*args, f'--value={mode}', f'--out_svg={single}'], tmp_path)
        if jobs:
            jobs.append('--next-plot')
        jobs.extend([*args, f'--value={mode}', f'--out_svg={tmp_path / f"batch{index}.svg"}',
                     f'--out_pdf={tmp_path / f"batch{index}.pdf"}'])
    run_r(['--batch', *jobs], tmp_path)
    for index in range(3):
        assert canonical_svg(tmp_path / f'single{index}.svg') == canonical_svg(tmp_path / f'batch{index}.svg')
        assert (tmp_path / f'batch{index}.pdf').read_bytes().startswith(b'%PDF')
    assert not list(tmp_path.glob('.gg-presence-*'))


@pytest.mark.parametrize('failure', ['render', 'duplicate', 'input_alias'])
def test_plot_batch_failure_preserves_all_previous_views(tmp_path, failure):
    args = summary_fixture(tmp_path)
    outputs = [tmp_path / 'first.svg', tmp_path / 'second.svg']
    for index, path in enumerate(outputs):
        path.write_text(f'previous{index}')
    second = [*args, f'--out_svg={outputs[1]}']
    if failure == 'render':
        second.append('--width=invalid')
    elif failure == 'duplicate':
        second.append(f'--out_svg={outputs[0]}')
    else:
        # The first output would replace the second view's input.
        second = [*args, f'--long_table={outputs[0]}', f'--out_svg={outputs[1]}']
    result = subprocess.run(['Rscript', str(SUPPORT / 'plot_query2family_presence_absence.R'), '--batch',
                             *args, f'--out_svg={outputs[0]}', '--next-plot', *second],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode != 0
    assert [path.read_text() for path in outputs] == ['previous0', 'previous1']
    assert not list(tmp_path.glob('.gg-presence-*'))


@pytest.mark.parametrize('restore_failure', [False, True])
def test_plot_batch_publication_failure_rolls_back_or_retains_backup(tmp_path, restore_failure):
    args = summary_fixture(tmp_path)
    outputs = [tmp_path / 'first.svg', tmp_path / 'second.svg']
    for index, path in enumerate(outputs):
        path.write_text(f'previous{index}')
    wrapper = tmp_path / 'publish.R'
    wrapper.write_text(
        f'source("{SUPPORT / "presence_absence_batch.R"}")\n'
        'file.rename <- function(from, to) {\n'
        '  if (basename(to) == "second.svg" && !grepl("backup", from)) return(FALSE)\n'
        + ('  if (grepl("backup", from)) return(FALSE)\n' if restore_failure else '') +
        '  base::file.rename(from, to)\n}\n'
        f'gg_render_presence_batch(commandArgs(TRUE), "{SUPPORT / "plot_query2family_presence_absence.R"}")\n'
    )
    result = subprocess.run(['Rscript', str(wrapper), *args, f'--out_svg={outputs[0]}',
                             '--next-plot', *args, f'--out_svg={outputs[1]}'],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode != 0
    assert outputs[1].read_text() == 'previous1'
    if restore_failure:
        backups = list(tmp_path.glob('.gg-presence-backup-*'))
        assert len(backups) == 1 and backups[0].read_text() == 'previous0'
        assert 'previous bytes retained at:' in result.stderr
    else:
        assert outputs[0].read_text() == 'previous0'
        assert not list(tmp_path.glob('.gg-presence-*'))


@pytest.mark.parametrize('basis', ['reference_species', 'query_gene', 'both'])
def test_gene_summary_batches_live_views_in_the_core(tmp_path, basis):
    from test_query_gene_orthologs import weak_duplication_fixture

    fixture = weak_duplication_fixture.__wrapped__(tmp_path)
    workspace = tmp_path / 'workspace'
    shutil.copytree(fixture['query_dir'], workspace / 'input/query_gene')
    shutil.copytree(fixture['output_root'], workspace / 'output/query2family')
    output = tmp_path / 'summary'
    env = os.environ.copy()
    env.update(gg_support_dir=str(SUPPORT), gg_workspace_dir=str(workspace),
               summary_output_dir=str(output), gene_family_source='query2family',
               run_species_taxonomy='0', run_family_completion_summary='0', run_presence_absence_summary='1',
               presence_absence_species_tree=str(fixture['species_tree']), presence_absence_species_tree_ci='',
               presence_absence_species_tree_support='', presence_absence_busco_table='',
               presence_absence_ortholog_basis=basis, presence_absence_include_incomplete='1',
               GG_COMMON_REFERENCE_SPECIES='Reference_species')
    result = subprocess.run(['bash', str(ROOT / 'workflow/core/gg_gene_summary_core.sh')],
                            env=env, cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (output / 'query2family_presence_absence.svg').is_file()
    for selected in (('reference_gene', 'query_gene') if basis == 'both' else
                     ('reference_gene',) if basis == 'reference_species' else ('query_gene',)):
        prefix = output / f'query2family_{selected}_orthologs'
        assert Path(str(prefix) + '.pdf').read_bytes().startswith(b'%PDF')
        assert len(pandas.read_csv(str(prefix) + '.columns.tsv', sep='\t')) == 2
    assert not list(output.glob('.gg-presence-*'))
