import ast
import gzip
import os
import re
import subprocess
import sys
from importlib.util import module_from_spec, spec_from_file_location
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


def test_generate_trait_colors_marks_every_nonbackground_lineage_id(monkeypatch, tmp_path):
    mod = load_module()
    monkeypatch.chdir(tmp_path)
    df_trait = pandas.DataFrame(
        {
            "species": ["sp0", "sp1", "sp2", "sp16", "sp_missing"],
            "trait": [0, 1, 2, 16, None],
        }
    )

    mod.generate_trait_colors(df_trait=df_trait, trait_names=["trait"])

    colors = pandas.read_csv(tmp_path / "trait_trait.color.tsv", sep="\t")
    assert colors["color"].tolist() == ["black", "firebrick", "firebrick", "firebrick", "black"]


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
    (tmp_path / "csubst.2XYZ.fa").write_text(">x\nAA\n", encoding="utf-8")
    assert extract_pdb_id(str(tmp_path)) == "2XYZ"


def test_extract_pdb_id_ignores_hidden_files_and_is_deterministic(tmp_path):
    extract_pdb_id = load_extract_pdb_id()
    (tmp_path / ".csubst.ZZZZ.fa").write_text(">x\nAA\n", encoding="utf-8")
    (tmp_path / "csubst.BBBB.fa").write_text(">x\nAA\n", encoding="utf-8")
    (tmp_path / "csubst.AAAA.fa").write_text(">x\nAA\n", encoding="utf-8")
    assert extract_pdb_id(str(tmp_path)) == "AAAA"


def test_resolve_site_output_dir_uses_current_csubst_namespace(tmp_path):
    mod = load_module()
    expected = (
        tmp_path
        / "OG0009394_2_66_70_79_99_116"
        / "csubst_sites"
        / "csubst.branch_id2,66,70,79,99,116"
    )

    out = mod.resolve_site_output_dir(
        str(tmp_path / "OG0009394_2_66_70_79_99_116"),
        "2,66,70,79,99,116",
    )

    assert out == str(expected)


def test_resolve_site_artifacts_reads_current_manifest_outputs(tmp_path):
    mod = load_module()
    site_dir = tmp_path / "OG1_1_2" / "csubst_sites" / "csubst.branch_id1,2"
    site_dir.mkdir(parents=True)
    site_tsv = site_dir / "csubst.2XYZ.tsv"
    site_pdf = site_dir / "csubst.2XYZ.pdf"
    pymol_pdf = site_dir / "csubst.2XYZ.pymol.pdf"
    site_tsv.write_text("codon_site_alignment\tOCNany2spe\n1\t1\n", encoding="utf-8")
    site_pdf.write_text("pdf", encoding="utf-8")
    pymol_pdf.write_text("pdf", encoding="utf-8")
    (site_dir / "csubst.2XYZ.fa").write_text(">x\nAA\n", encoding="utf-8")
    pandas.DataFrame(
        [
            {"output_kind": "site_table_tsv", "output_path": str(site_tsv), "file_exists": "Y"},
            {"output_kind": "site_summary_pdf", "output_path": str(site_pdf), "file_exists": "Y"},
            {"output_kind": "pymol_summary_pdf", "output_path": str(pymol_pdf), "file_exists": "Y"},
        ]
    ).to_csv(site_dir / "csubst.outputs.tsv", sep="\t", index=False)

    artifacts = mod.resolve_site_artifacts(str(tmp_path / "OG1_1_2"), "1,2")

    assert artifacts["site_dir"] == str(site_dir)
    assert artifacts["site_table_tsv"] == str(site_tsv)
    assert artifacts["site_summary_pdf"] == str(site_pdf)
    assert artifacts["pymol_summary_pdf"] == str(pymol_pdf)
    assert artifacts["pdb_id"] == "2XYZ"


def test_resolve_site_artifacts_falls_back_to_current_names_without_manifest(tmp_path):
    mod = load_module()
    site_dir = tmp_path / "OG1_1_2" / "csubst_sites" / "csubst.branch_id1,2"
    site_dir.mkdir(parents=True)
    site_tsv = site_dir / "csubst.2XYZ.tsv"
    site_pdf = site_dir / "csubst.2XYZ.pdf"
    pymol_pdf = site_dir / "csubst.2XYZ.pymol.pdf"
    site_tsv.write_text("codon_site_alignment\tOCNany2spe\n1\t1\n", encoding="utf-8")
    site_pdf.write_text("pdf", encoding="utf-8")
    pymol_pdf.write_text("pdf", encoding="utf-8")
    (site_dir / "csubst.2XYZ.fa").write_text(">x\nAA\n", encoding="utf-8")

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
        "omegaCany2any",
        "omegaCdif2spe",
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
        "omegaCany2any",
        "omegaCdif2spe",
        "OCNCoD",
        "branch_id_1",
        "branch_id_2",
        "is_fg_traitA",
        "branch_num_fg_stem_traitA",
    ]


def test_get_annotation_text_includes_observed_branch_statistics():
    mod = load_module()

    text = mod.get_annotation_text(
        og="OG0001",
        arity=2,
        branch_id_str="1,2",
        trait="traitA",
        min_OCNany2spe=1.8,
        min_omegaCany2spe=3.0,
        min_OCNCoD=0.0,
        besthit_values=["a", "b", "c", "d", "e"],
        observed_values={
            "OCNany2spe": 2.4,
            "ECNany2spe": 0.8,
            "OCSany2spe": 1.0,
            "ECSany2spe": 0.5,
            "omegaCany2spe": 5.6,
            "omegaCany2any": 4.7,
            "omegaCdif2spe": 2.3,
            "OCNCoD": 1.2,
        },
    )

    assert "Threshold\n\nOCNany2spe: 1.8" in text
    assert "omegaCany2spe: 3.0" in text
    assert "OCNCoD: 0.0" in text
    assert "Observed\n\nOCNany2spe: 2.4" in text
    assert "ECNany2spe: 0.8" in text
    assert "OCSany2spe: 1.0" in text
    assert "ECSany2spe: 0.5" in text
    assert "omegaCany2spe: 5.6" in text
    assert "omegaCany2any: 4.7" in text
    assert "omegaCdif2spe: 2.3" in text
    assert "OCNCoD: 1.2" in text


def test_append_csubst_sites_command_quotes_command_in_annotation_text():
    mod = load_module()

    text = mod.append_csubst_sites_command(
        "Annotation body",
        ["csubst", "sites", "--branch_id", "1,2", "--alignment_file", "dir with spaces/csubst.fasta"],
    )

    assert "Annotation body" in text
    assert "CSUBST sites command" in text
    assert "csubst sites --branch_id 1,2 --alignment_file 'dir with spaces/csubst.fasta'" in text


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


def test_write_recoded_site_alignment_uses_csubst_state_symbols(tmp_path):
    mod = load_module()
    recoding_table = tmp_path / "csubst_nonsyn_recoding.tsv"
    rows = [
        ("1", "AGPST", "A"),
        ("1", "AGPST", "G"),
        ("1", "AGPST", "P"),
        ("1", "AGPST", "S"),
        ("1", "AGPST", "T"),
        ("2", "DENQ", "D"),
        ("2", "DENQ", "E"),
        ("2", "DENQ", "N"),
        ("2", "DENQ", "Q"),
        ("3", "HKR", "H"),
        ("3", "HKR", "K"),
        ("3", "HKR", "R"),
        ("4", "ILMV", "I"),
        ("4", "ILMV", "L"),
        ("4", "ILMV", "M"),
        ("4", "ILMV", "V"),
        ("5", "FWY", "F"),
        ("5", "FWY", "W"),
        ("5", "FWY", "Y"),
        ("6", "C", "C"),
    ]
    recoding_table.write_text(
        "state_id\tstate_label\tamino_acid\n"
        + "\n".join("\t".join(row) for row in rows)
        + "\n",
        encoding="utf-8",
    )
    codon_alignment = tmp_path / "csubst.fasta"
    codon_alignment.write_text(
        ">gene1\nGCTGATAAAATTTTTTGT\n"
        ">gene2\nTCTGAACATGTGTATTGC\n",
        encoding="utf-8",
    )
    out = tmp_path / "recoded.fasta"

    result = mod.write_recoded_site_alignment(str(codon_alignment), str(recoding_table), str(out))

    assert result == str(out)
    assert out.read_text(encoding="utf-8") == ">gene1\nABCDEF\n>gene2\nABCDEF\n"


def test_build_tree_plot_panel_args_adds_recoded_state_panel_when_available():
    mod = load_module()

    out = mod.build_tree_plot_panel_args(
        file_og_rpsblast="rps.tsv",
        file_csubst_input_fasta="csubst.fasta",
        convergent_site_str="2:5",
        file_og_alignment="trim.fa",
        file_og_untrimmed_alignment="untrim.fa",
        recoded_site_alignment="dayhoff6.fa",
        csubst_nonsyn_recode="dayhoff6",
    )

    assert "--panel10=amino_acid_site,1,2:5,csubst.fasta" in out
    assert "--panel11=site_state,site_state_dayhoff6,2:5,dayhoff6.fa,Recoded state (dayhoff6)" in out
    assert "--panel12=alignment,trim.fa,untrim.fa" in out
    assert "--panel13=fimo,2000,0.05" in out


def test_build_tree_plot_panel_args_keeps_existing_numbering_without_recoded_panel():
    mod = load_module()

    out = mod.build_tree_plot_panel_args(
        file_og_rpsblast="rps.tsv",
        file_csubst_input_fasta="csubst.fasta",
        convergent_site_str="2:5",
        file_og_alignment="trim.fa",
        file_og_untrimmed_alignment=None,
        recoded_site_alignment=None,
        csubst_nonsyn_recode="no",
    )

    assert "--panel10=amino_acid_site,1,2:5,csubst.fasta" in out
    assert "--panel11=alignment,trim.fa" in out
    assert "--panel12=fimo,2000,0.05" in out
    assert not any(arg.startswith("--panel11=site_state") for arg in out)


def test_csubst_nonsyn_recode_output_suffix_preserves_default_name():
    mod = load_module()

    assert mod.csubst_nonsyn_recode_output_suffix("no") == ""
    assert mod.csubst_nonsyn_recode_output_suffix("Dayhoff6") == "_nonsynRecode-dayhoff6"


def test_build_csubst_sites_command_passes_nondefault_nonsyn_recode_only():
    mod = load_module()

    default_cmd = mod.build_csubst_sites_command("OG0001.iqtree.anc", "/tmp/OG0001.iqtree.anc", "1,2", 2, "no")
    recoded_cmd = mod.build_csubst_sites_command("OG0001.iqtree.anc", "/tmp/OG0001.iqtree.anc", "1,2", 2, "dayhoff6")

    assert "--nonsyn_recode" not in default_cmd
    assert recoded_cmd[recoded_cmd.index("--nonsyn_recode") + 1] == "dayhoff6"
    assert default_cmd[default_cmd.index("--outdir") + 1] == "csubst_sites"
    assert default_cmd[default_cmd.index("--output_prefix") + 1] == "csubst"
    assert default_cmd[default_cmd.index("--pdb") + 1] == "besthit"


def test_build_csubst_sites_command_can_disable_pdb_search():
    mod = load_module()

    command = mod.build_csubst_sites_command(
        "OG0001.iqtree.anc",
        "/tmp/OG0001.iqtree.anc",
        "1,2",
        2,
        "no",
        pdb="none",
    )

    assert "--pdb" not in command


def test_run_stat_branch2tree_plot_accepts_one_explicit_site(monkeypatch, tmp_path):
    mod = load_module()
    monkeypatch.chdir(tmp_path)
    site_table = tmp_path / "csubst.tsv"
    site_table.write_text(
        "codon_site_alignment\tOCNany2spe\n2\t2.0\n5\t2.0\n",
        encoding="utf-8",
    )
    iqtree_dir = tmp_path / "OG0001.iqtree.anc"
    iqtree_dir.mkdir()
    (iqtree_dir / "csubst.fasta").write_text(">a\nATG\n", encoding="utf-8")
    monkeypatch.setattr(mod, "get_stat_branch_path", lambda **kwargs: "stat.tsv")
    monkeypatch.setattr(mod, "get_rpsblast_path", lambda **kwargs: "rps.tsv")
    monkeypatch.setattr(mod, "get_alignment_for_tree_plot", lambda **kwargs: "trim.fa")
    monkeypatch.setattr(mod, "get_untrimmed_alignment_for_tree_plot", lambda **kwargs: None)
    monkeypatch.setattr(mod, "get_iqtree_anc_dir", lambda **kwargs: str(iqtree_dir))
    monkeypatch.setattr(
        mod,
        "resolve_site_artifacts",
        lambda **kwargs: {"site_table_tsv": str(site_table), "site_dir": str(tmp_path)},
    )
    monkeypatch.setattr(mod, "prepare_recoded_site_alignment", lambda **kwargs: None)
    commands = []

    def fake_run(command, check):
        commands.append(command)
        (tmp_path / "stat_branch2tree_plot.pdf").write_bytes(b"focused")

    monkeypatch.setattr(mod.subprocess, "run", fake_run)
    output = tmp_path / "focused.pdf"

    mod.run_stat_branch2tree_plot(
        og="OG0001",
        branch_id_str="1,2",
        file_trait_color="trait.tsv",
        dir_out_og=str(tmp_path),
        dir_og=str(tmp_path),
        convergent_sites=[5],
        file_tree_plot_out=str(output),
    )

    assert output.read_bytes() == b"focused"
    assert "--panel10=amino_acid_site,1,5," + str(iqtree_dir / "csubst.fasta") in commands[0]
    assert not any("amino_acid_site,1,2:5," in token for token in commands[0])

    site_table.write_text(
        "codon_site_alignment\tOCNany2spe\n2\t0.0\n5\t0.0\n",
        encoding="utf-8",
    )
    legacy_output = tmp_path / "legacy-empty-sites.pdf"
    mod.run_stat_branch2tree_plot(
        og="OG0001",
        branch_id_str="1,2",
        file_trait_color="trait.tsv",
        dir_out_og=str(tmp_path),
        dir_og=str(tmp_path),
        file_tree_plot_out=str(legacy_output),
    )
    assert legacy_output.read_bytes() == b"focused"

    with pytest.raises(ValueError, match="No convergent sites"):
        mod.run_stat_branch2tree_plot(
            og="OG0001",
            branch_id_str="1,2",
            file_trait_color="trait.tsv",
            dir_out_og=str(tmp_path),
            dir_og=str(tmp_path),
            convergent_sites=[],
            file_tree_plot_out=str(tmp_path / "explicit-empty-sites.pdf"),
        )


def test_get_alignment_for_tree_plot_rejects_legacy_dot_clipkit_name(tmp_path):
    mod = load_module()
    dir_og = tmp_path / "orthogroup"
    dir_out_og = tmp_path / "csubst_site" / "OG0001_1_2"
    clipkit_dir = dir_og / "clipkit"
    clipkit_dir.mkdir(parents=True)
    (clipkit_dir / "OG0001.cds.clipkit.fa").write_text(">gene1\nATG\n", encoding="utf-8")

    with pytest.raises(FileNotFoundError):
        mod.get_alignment_for_tree_plot(str(dir_og), "OG0001", str(dir_out_og))


def test_materialize_csubst_site_inputs_reads_only_requested_family_from_zip(tmp_path):
    from workflow.support.gene_family_output_store import archive_completed_outputs, family_context

    mod = load_module()
    dir_og = tmp_path / "orthogroup"
    genecount = tmp_path / "Orthogroups.GeneCount.selected.tsv"
    genecount.write_text("Orthogroup\tTotal\nOG0001\t1\nOG0002\t1\n", encoding="utf-8")
    for family_id in ("OG0001", "OG0002"):
        for subdir, suffix in (
            ("stat_branch", "_stat.branch.tsv"),
            ("stat_tree", "_stat.tree.tsv"),
            ("tree_plot", "_tree_plot.pdf"),
            ("iqtree_anc", "_iqtree.anc.zip"),
            ("rpsblast", "_rpsblast.tsv"),
            ("clipkit", "_cds.clipkit.fa.gz"),
        ):
            path = dir_og / subdir / f"{family_id}{suffix}"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(f"{family_id}:{subdir}\n".encode())
    family_ids, family_from_name = family_context("orthogroup", genecount=genecount)
    archive_completed_outputs(dir_og, "orthogroup", family_ids, family_from_name)

    destination = tmp_path / "materialized"
    materialized = mod.materialize_csubst_site_inputs(
        str(dir_og),
        "OG0001",
        str(destination),
    )

    assert materialized
    assert (destination / "iqtree_anc" / "OG0001_iqtree.anc.zip").is_file()
    assert (destination / "clipkit" / "OG0001_cds.clipkit.fa.gz").is_file()
    assert not (destination / "stat_branch" / "OG0002_stat.branch.tsv").exists()


def test_process_og_batch_materializes_each_family_only_once(tmp_path, monkeypatch):
    mod = load_module()
    dir_og = tmp_path / "orthogroup"
    (dir_og / ".gg_store").mkdir(parents=True)
    materialize_calls = []
    process_calls = []

    def fake_materialize(dir_og, og, destination_root):
        materialize_calls.append((dir_og, og, destination_root))
        return []

    def fake_process(
        og,
        branch_id_str,
        dir_out,
        effective_dir_og,
        file_trait_color,
        ncpu,
        csubst_nonsyn_recode,
        annotation_text,
    ):
        assert Path(effective_dir_og).is_dir()
        process_calls.append((og, branch_id_str, effective_dir_og, annotation_text))
        return og, None

    monkeypatch.setattr(mod, "materialize_csubst_site_inputs", fake_materialize)
    monkeypatch.setattr(mod, "process_index", fake_process)

    results = mod.process_og_batch(
        "OG0001",
        ["1,2", "3,4"],
        str(tmp_path / "out"),
        str(dir_og),
        "trait.tsv",
        2,
        "no",
        ["first", "second"],
    )

    assert results == [("OG0001", None), ("OG0001", None)]
    assert len(materialize_calls) == 1
    assert [call[1] for call in process_calls] == ["1,2", "3,4"]
    assert process_calls[0][2] == process_calls[1][2]
    assert not Path(process_calls[0][2]).exists()
    mod.cleanup_materialization_metadata(tmp_path / "out")
    assert not (tmp_path / "out" / ".gg_materialized.lock").exists()
    assert (tmp_path / ".gg_materialized.lock").is_file()


def test_locked_materialization_directory_reclaims_stale_runs(tmp_path):
    mod = load_module()
    output = tmp_path / "out"
    stale = output / ".gg_materialized" / "stale"
    stale.mkdir(parents=True)
    (stale / ".run.lock").write_bytes(b"")
    (stale / "large-input.bin").write_bytes(b"x" * 1024)

    managed = mod.LockedMaterializationDirectory(output, "OG0001")
    try:
        assert not stale.exists()
        assert Path(managed.name).is_dir()
    finally:
        managed.cleanup()
        mod.cleanup_materialization_metadata(output)

    assert not (output / ".gg_materialized").exists()
    assert not (output / ".gg_materialized.lock").exists()
    assert (tmp_path / ".gg_materialized.lock").is_file()


def test_locked_materialization_directory_preserves_concurrent_runs(tmp_path):
    mod = load_module()
    output = tmp_path / "out"

    first = mod.LockedMaterializationDirectory(output, "OG0001")
    second = mod.LockedMaterializationDirectory(output, "OG0002")
    try:
        assert Path(first.name).is_dir()
        assert Path(second.name).is_dir()
        first.cleanup()
        assert not Path(first.name).exists()
        assert Path(second.name).is_dir()
    finally:
        first.cleanup()
        second.cleanup()

    assert not (output / ".gg_materialized").exists()
    assert (tmp_path / ".gg_materialized.lock").is_file()


def test_process_index_restores_cwd_when_summary_already_exists(tmp_path):
    mod = load_module()
    output_root = tmp_path / "out"
    output = output_root / "OG0001_1_2"
    output.mkdir(parents=True)
    (output / "summary.OG0001_branch_id1,2.pdf").write_bytes(b"pdf")
    original_cwd = os.getcwd()

    assert mod.process_index(
        "OG0001",
        "1,2",
        str(output_root),
        str(tmp_path / "orthogroup"),
        "trait.tsv",
        1,
        "no",
        "annotation",
    ) == ("OG0001", None)
    assert os.getcwd() == original_cwd


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


def test_raise_on_processing_failures_propagates_worker_errors():
    mod = load_module()

    with pytest.raises(RuntimeError, match="OG0002: failed site analysis"):
        mod.raise_on_processing_failures(
            [("OG0002", RuntimeError("failed site analysis"))]
        )


def test_raise_on_processing_failures_accepts_success():
    mod = load_module()

    assert mod.raise_on_processing_failures([]) is None
