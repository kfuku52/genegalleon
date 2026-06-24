import ast
import gzip
from importlib.util import module_from_spec, spec_from_file_location
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas
import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "csubst_site_wrapper.py"


def load_module():
    spec = spec_from_file_location("csubst_site_wrapper", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_extract_pdb_id():
    source = SCRIPT_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(SCRIPT_PATH))
    extract_node = None
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == "extract_pdb_id":
            extract_node = node
            break
    if extract_node is None:
        raise AssertionError("extract_pdb_id function not found in csubst_site_wrapper.py")
    module_tree = ast.Module(body=[extract_node], type_ignores=[])
    code = compile(module_tree, filename=str(SCRIPT_PATH), mode="exec")
    namespace = {"os": os, "re": re}
    exec(code, namespace)
    return namespace["extract_pdb_id"]


def test_get_matplotlib_imports_plotting_submodules():
    mod = load_module()

    matplotlib = mod._get_matplotlib()

    assert hasattr(matplotlib, "pyplot")
    assert hasattr(matplotlib, "patches")


def test_extract_pdb_id_returns_none_for_missing_directory(tmp_path):
    extract_pdb_id = load_extract_pdb_id()
    missing = tmp_path / "missing"
    assert extract_pdb_id(str(missing)) is None


def test_extract_pdb_id_returns_none_when_no_match(tmp_path):
    extract_pdb_id = load_extract_pdb_id()
    (tmp_path / "README.txt").write_text("x", encoding="utf-8")
    (tmp_path / "csubst_site.tsv").write_text("x", encoding="utf-8")
    assert extract_pdb_id(str(tmp_path)) is None


def test_extract_pdb_id_parses_expected_identifier(tmp_path):
    extract_pdb_id = load_extract_pdb_id()
    (tmp_path / "csubst_sites.2XYZ.fa").write_text(">x\nAA\n", encoding="utf-8")
    assert extract_pdb_id(str(tmp_path)) == "2XYZ"


def test_extract_pdb_id_ignores_hidden_files_and_is_deterministic(tmp_path):
    extract_pdb_id = load_extract_pdb_id()
    (tmp_path / ".csubst_sites.ZZZZ.fa").write_text(">x\nAA\n", encoding="utf-8")
    (tmp_path / "csubst_sites.BBBB.fa").write_text(">x\nAA\n", encoding="utf-8")
    (tmp_path / "csubst_sites.AAAA.fa").write_text(">x\nAA\n", encoding="utf-8")
    assert extract_pdb_id(str(tmp_path)) == "AAAA"


def test_resolve_site_output_dir_prefers_latest_directory_name(tmp_path):
    mod = load_module()
    latest = tmp_path / "OG1_1_2" / "csubst_sites.branch_id1,2"
    legacy = tmp_path / "OG1_1_2" / "csubst_site.branch_id1,2"
    latest.mkdir(parents=True)
    legacy.mkdir(parents=True)

    out = mod.resolve_site_output_dir(str(tmp_path / "OG1_1_2"), "1,2")

    assert out == str(latest)


def test_resolve_site_artifacts_reads_latest_manifest_outputs(tmp_path):
    mod = load_module()
    site_dir = tmp_path / "OG1_1_2" / "csubst_sites.branch_id1,2"
    site_dir.mkdir(parents=True)
    site_tsv = site_dir / "csubst_sites.tsv"
    site_pdf = site_dir / "csubst_sites.pdf"
    pymol_pdf = site_dir / "csubst_sites.2XYZ.pymol.pdf"
    site_tsv.write_text("codon_site_alignment\tOCNany2spe\n1\t1\n", encoding="utf-8")
    site_pdf.write_text("pdf", encoding="utf-8")
    pymol_pdf.write_text("pdf", encoding="utf-8")
    (site_dir / "csubst_sites.2XYZ.fa").write_text(">x\nAA\n", encoding="utf-8")
    pandas.DataFrame(
        [
            {"output_kind": "site_table_tsv", "output_path": str(site_tsv), "file_exists": "Y"},
            {"output_kind": "site_summary_pdf", "output_path": str(site_pdf), "file_exists": "Y"},
            {"output_kind": "pymol_summary_pdf", "output_path": str(pymol_pdf), "file_exists": "Y"},
        ]
    ).to_csv(site_dir / "csubst_sites.outputs.tsv", sep="\t", index=False)

    artifacts = mod.resolve_site_artifacts(str(tmp_path / "OG1_1_2"), "1,2")

    assert artifacts["site_dir"] == str(site_dir)
    assert artifacts["site_table_tsv"] == str(site_tsv)
    assert artifacts["site_summary_pdf"] == str(site_pdf)
    assert artifacts["pymol_summary_pdf"] == str(pymol_pdf)
    assert artifacts["pdb_id"] == "2XYZ"


def test_resolve_site_artifacts_falls_back_to_legacy_names_without_manifest(tmp_path):
    mod = load_module()
    site_dir = tmp_path / "OG1_1_2" / "csubst_site.branch_id1,2"
    site_dir.mkdir(parents=True)
    site_tsv = site_dir / "csubst_site.2XYZ.tsv"
    site_pdf = site_dir / "csubst_site.2XYZ.pdf"
    pymol_pdf = site_dir / "csubst_site.2XYZ.pymol.pdf"
    site_tsv.write_text("codon_site_alignment\tOCNany2spe\n1\t1\n", encoding="utf-8")
    site_pdf.write_text("pdf", encoding="utf-8")
    pymol_pdf.write_text("pdf", encoding="utf-8")
    (site_dir / "csubst_site.2XYZ.fa").write_text(">x\nAA\n", encoding="utf-8")

    artifacts = mod.resolve_site_artifacts(str(tmp_path / "OG1_1_2"), "1,2")

    assert artifacts["site_dir"] == str(site_dir)
    assert artifacts["site_table_tsv"] == str(site_tsv)
    assert artifacts["site_summary_pdf"] == str(site_pdf)
    assert artifacts["pymol_summary_pdf"] == str(pymol_pdf)
    assert artifacts["pdb_id"] == "2XYZ"


def test_get_cb_required_columns_keeps_needed_columns_only():
    mod = load_module()
    cols = [
        "orthogroup",
        "OCNany2spe",
        "ECNany2spe",
        "OCSany2spe",
        "ECSany2spe",
        "omegaCany2spe",
        "OCNCoD",
        "branch_id_1",
        "branch_id_2",
        "is_fg_traitA",
        "branch_num_fg_stem_traitA",
        "unused_col",
    ]

    out = mod.get_cb_required_columns(cols, ["traitA"])

    assert out == [
        "orthogroup",
        "OCNany2spe",
        "ECNany2spe",
        "OCSany2spe",
        "ECSany2spe",
        "omegaCany2spe",
        "OCNCoD",
        "branch_id_1",
        "branch_id_2",
        "is_fg_traitA",
        "branch_num_fg_stem_traitA",
    ]


def test_skip_lower_order_filters_subset_rows_per_orthogroup():
    mod = load_module()
    cb_passed = pandas.DataFrame(
        [
            {"orthogroup": "OG1", "branch_id_1": 1, "branch_id_2": 2},
            {"orthogroup": "OG1", "branch_id_1": 2, "branch_id_2": 3},
            {"orthogroup": "OG2", "branch_id_1": 5, "branch_id_2": 6},
        ]
    )
    already = {
        "traitA": {
            "OG1": [frozenset([1, 2, 9])],
            "OG2": [frozenset([7, 8])],
        }
    }

    out, updated = mod.skip_lower_order(cb_passed, 2, "traitA", already)

    assert out["orthogroup"].tolist() == ["OG1", "OG2"]
    assert out.loc[0, ["branch_id_1", "branch_id_2"]].tolist() == [2, 3]
    assert out.loc[1, ["branch_id_1", "branch_id_2"]].tolist() == [5, 6]
    assert frozenset([2, 3]) in updated["traitA"]["OG1"]
    assert frozenset([5, 6]) in updated["traitA"]["OG2"]


def test_load_annotation_besthits_reads_besthit_columns_only(tmp_path):
    mod = load_module()
    annot_dir = tmp_path / "Orthogroups_filtered"
    annot_dir.mkdir(parents=True)
    infile = annot_dir / "Orthogroups.GeneCount.annotated.tsv"
    infile.write_text(
        "Orthogroup\tTotal\tbesthit_0.5\tbesthit_0.95\n"
        "OG1\t4\thitA\thitB\n",
        encoding="utf-8",
    )

    out = mod.load_annotation_besthits(str(tmp_path))

    assert out.columns.tolist() == ["orthogroup", "besthit_0.5", "besthit_0.95"]
    assert out.loc[0, "orthogroup"] == "OG1"


def test_load_annotation_besthits_requires_filtered_directory(tmp_path):
    mod = load_module()
    legacy_dir = tmp_path / "Orthogroups"
    legacy_dir.mkdir(parents=True)
    (legacy_dir / "Orthogroups.GeneCount.annotated.tsv").write_text(
        "Orthogroup\tTotal\tbesthit_0.5\n"
        "OG1\t4\tlegacy_hit\n",
        encoding="utf-8",
    )

    with pytest.raises(FileNotFoundError) as excinfo:
        mod.load_annotation_besthits(str(tmp_path))

    expected = tmp_path / "Orthogroups_filtered" / "Orthogroups.GeneCount.annotated.tsv"
    assert Path(excinfo.value.filename) == expected


def test_gene_evolution_artifact_paths_use_current_underscore_names(tmp_path):
    mod = load_module()

    dir_og = tmp_path / "orthogroup"
    dir_out_og = tmp_path / "csubst_site" / "OG0001_1_2"

    assert mod.get_iqtree_anc_zip_path(str(dir_og), "OG0001") == str(
        dir_og / "iqtree_anc" / "OG0001_iqtree.anc.zip"
    )
    assert mod.get_iqtree_anc_dir(str(dir_out_og), "OG0001") == str(
        dir_out_og / "OG0001.iqtree.anc"
    )
    assert mod.get_stat_branch_path(str(dir_og), "OG0001") == str(
        dir_og / "stat_branch" / "OG0001_stat.branch.tsv"
    )
    assert mod.get_rpsblast_path(str(dir_og), "OG0001") == str(
        dir_og / "rpsblast" / "OG0001_rpsblast.tsv"
    )


def test_get_alignment_for_tree_plot_reads_current_clipkit_name(tmp_path):
    mod = load_module()
    dir_og = tmp_path / "orthogroup"
    dir_out_og = tmp_path / "csubst_site" / "OG0001_1_2"
    clipkit_dir = dir_og / "clipkit"
    clipkit_dir.mkdir(parents=True)
    dir_out_og.mkdir(parents=True)
    alignment = clipkit_dir / "OG0001_cds.clipkit.fa.gz"
    with gzip.open(alignment, "wt") as handle:
        handle.write(">gene1\nATG\n")

    out = mod.get_alignment_for_tree_plot(str(dir_og), "OG0001", str(dir_out_og))

    assert out == str(dir_out_og / "OG0001_cds.clipkit.plot.fasta")
    assert Path(out).read_text(encoding="utf-8") == ">gene1\nATG\n"


def test_get_untrimmed_alignment_for_tree_plot_reads_current_mafft_name(tmp_path):
    mod = load_module()
    dir_og = tmp_path / "orthogroup"
    dir_out_og = tmp_path / "csubst_site" / "OG0001_1_2"
    mafft_dir = dir_og / "mafft"
    mafft_dir.mkdir(parents=True)
    dir_out_og.mkdir(parents=True)
    alignment = mafft_dir / "OG0001_cds.aln.fa.gz"
    with gzip.open(alignment, "wt") as handle:
        handle.write(">gene1\nATG---AAA\n")

    out = mod.get_untrimmed_alignment_for_tree_plot(str(dir_og), "OG0001", str(dir_out_og))

    assert out == str(dir_out_og / "OG0001_cds.untrimmed.plot.fasta")
    assert Path(out).read_text(encoding="utf-8") == ">gene1\nATG---AAA\n"


def test_get_untrimmed_alignment_for_tree_plot_returns_none_when_missing(tmp_path):
    mod = load_module()
    dir_og = tmp_path / "orthogroup"
    dir_out_og = tmp_path / "csubst_site" / "OG0001_1_2"
    dir_og.mkdir(parents=True)
    dir_out_og.mkdir(parents=True)

    assert mod.get_untrimmed_alignment_for_tree_plot(str(dir_og), "OG0001", str(dir_out_og)) is None


def test_build_alignment_panel_arg_includes_untrimmed_when_available():
    mod = load_module()

    assert mod.build_alignment_panel_arg("trim.fa", "untrim.fa") == "--panel11=alignment,trim.fa,untrim.fa"
    assert mod.build_alignment_panel_arg("trim.fa", None) == "--panel11=alignment,trim.fa"


def test_get_alignment_for_tree_plot_rejects_legacy_dot_clipkit_name(tmp_path):
    mod = load_module()
    dir_og = tmp_path / "orthogroup"
    dir_out_og = tmp_path / "csubst_site" / "OG0001_1_2"
    clipkit_dir = dir_og / "clipkit"
    clipkit_dir.mkdir(parents=True)
    (clipkit_dir / "OG0001.cds.clipkit.fa").write_text(">gene1\nATG\n", encoding="utf-8")

    with pytest.raises(FileNotFoundError):
        mod.get_alignment_for_tree_plot(str(dir_og), "OG0001", str(dir_out_og))


def test_help_has_no_side_effect_files(tmp_path):
    proc = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "--help"],
        cwd=str(tmp_path),
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0
    assert not (tmp_path / "generate_orthogroup_database.log").exists()
    assert not (tmp_path / "mpl").exists()
