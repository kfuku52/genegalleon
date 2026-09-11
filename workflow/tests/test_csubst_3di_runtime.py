"""Real predictor test; opt in with GG_TEST_CSUBST_3DI=1 in a GeneGalleon runtime."""
import copy
import json
import os
import shutil
import subprocess
import zipfile
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from workflow.support import csubst_input_bundle as bundle
from workflow.support import csubst_scan_candidate_sites as candidates
from workflow.support import csubst_site_wrapper as sites

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow/support"

pytestmark = pytest.mark.skipif(os.environ.get("GG_TEST_CSUBST_3DI") != "1", reason="Set GG_TEST_CSUBST_3DI=1 for real 3Di predictor validation")


def run(command, cwd, env, name):
    result = subprocess.run(command, cwd=cwd, env=env, capture_output=True, text=True, timeout=300)
    (cwd / (name + ".log")).write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    return result


@pytest.mark.parametrize("code", [1, 2])
def test_full_bundle_search_scan_and_sites_with_real_predictor(tmp_path, code, monkeypatch):
    core = (ROOT / "workflow/core/gg_gene_evolution_core.sh").read_text()
    fixture = ROOT / "workflow/tests/data/csubst_scan_inference"
    shutil.copyfile(fixture / "input.fa", tmp_path / "full.fa")
    records = list(SeqIO.parse(tmp_path / "full.fa", "fasta"))
    if code == 2:
        # Use sense codons under the mitochondrial code, including its TGA Trp.
        # Reinterpreting the standard-code fixture leaves AGA/AGG stops and no
        # TGA observations, which makes its empirical-frequency fit degenerate.
        for index, record in enumerate(records):
            codons = [str(record.seq)[i:i + 3] for i in range(0, len(record.seq), 3)]
            replacements = {"AGA": "CGA", "AGG": "CGG"}
            if index % 2 == 0:
                replacements["TGG"] = "TGA"
            record.seq = Seq("".join(replacements.get(codon, codon) for codon in codons))
        records[1].seq = Seq(str(records[1].seq[:267]) + "---" + str(records[1].seq[270:]))
        records[2].seq = Seq(str(records[2].seq[:282]) + "NNN" + str(records[2].seq[285:]))
        extra = copy.deepcopy(records[0])
        extra.id, extra.description = "excluded_tip", ""
        records.append(extra)
        SeqIO.write(records, tmp_path / "full.fa", "fasta")
    # Removing ten codons makes trimmed site coordinates demonstrably different.
    for record in records:
        record.seq = record.seq[30:]
    SeqIO.write(records, tmp_path / "trimmed.fa", "fasta")
    shutil.copyfile(fixture / "tree.nwk", tmp_path / "rooted.nwk")
    codon_model = "ECMK07+F+R4" if code == 1 else "GY+F+R4"
    env = dict(os.environ, csubst_nonsyn_recode="3di20", genetic_code=str(code),
               GG_TASK_CPUS="2", codon_model=codon_model, gg_support_dir=str(SUPPORT), og_id="OG3DI",
               file_og_trimmed_aln_analysis=str(tmp_path / "trimmed.fa"),
               file_og_untrimmed_aln_analysis=str(tmp_path / "full.fa"),
               file_og_rooted_tree_analysis=str(tmp_path / "rooted.nwk"),
               csubst_max_arity="2", csubst_exhaustive_until="2", csubst_cutoff_stat="OCNany2spe,0",
               csubst_max_combination="10000", csubst_fg_exclude_wg="no", csubst_fg_stem_only="no",
               csubst_scan_unit_mode="clade", csubst_scan_match="any2spe", csubst_scan_min_event_pp="0.5",
               csubst_scan_min_support="2", csubst_scan_other_scope="all", csubst_scan_site_plot="yes",
               csubst_scan_tree_site_plot_format="pdf", csubst_scan_tree_site_plot_max_sites="30")
    start = core.index('  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "tmp.csubst.fasta"')
    end = core.index('    if [[ -e "${og_id}.iqtree.anc"', start)
    script = 'build_iqtree_mem_args() { IQTREE_MEM_ARGS=(); }\n' + core[start:end] + '\nfi\n'
    run(["bash", "-euo", "pipefail", "-c", script], tmp_path, env, "bundle")
    origin = tmp_path / "original"
    origin.mkdir()
    for path in tmp_path.glob("csubst.*"):
        shutil.move(path, origin / path.name)
    # Archive and extract to a new location before every analysis, like real sites.
    with zipfile.ZipFile(tmp_path / "bundle.zip", "w") as archive:
        for path in origin.rglob("*"):
            if path.is_file():
                archive.write(path, Path("OG3DI.iqtree.anc") / path.relative_to(origin))
    shutil.rmtree(origin)
    with zipfile.ZipFile(tmp_path / "bundle.zip") as archive:
        archive.extractall(tmp_path)
    base = tmp_path / "OG3DI.iqtree.anc"
    full_dir = bundle.structural_directory(base)
    assert len(next(SeqIO.parse(full_dir / "csubst.fasta", "fasta")).seq) == 300
    full_records = {record.id: str(record.seq) for record in SeqIO.parse(full_dir / "csubst.fasta", "fasta")}
    assert len(full_records) == 8 and "excluded_tip" not in full_records
    if code == 2:
        assert full_records["b"][267:270] == "---" and full_records["c"][282:285] == "NNN"
    assert len(next(SeqIO.parse(base / "csubst.fasta", "fasta")).seq) == 270
    assert sites.get_csubst_branch_clade_signatures(base) == sites.get_csubst_branch_clade_signatures(full_dir)
    (tmp_path / "foreground.tsv").write_text("name\ttrait\na\t1\nb\t1\nc\t0\nd\t0\ne\t1\nf\t1\ng\t0\nh\t0\n")
    # Cached model resources must suffice after bundle generation.
    env.update(HF_HUB_OFFLINE="1", TRANSFORMERS_OFFLINE="1")
    start = core.index('  csubst_input_base="./${og_id}.iqtree.anc/csubst"')
    end = core.index('\n  csubst_b_src=', start)
    result = run(["bash", "-euo", "pipefail", "-c", 'foreground_params=()\n' + core[start:end]], tmp_path, env, "search")
    assert "Loaded 3Di state cache" in result.stdout
    scan_start = core.index('  csubst_input_base="./${og_id}.iqtree.anc/csubst"', end)
    scan_end = core.index('\n  if [[ -s "${csubst_scan_dir}/csubst_scan.tsv"', scan_start)
    result = run(["bash", "-euo", "pipefail", "-c", core[scan_start:scan_end]], tmp_path, env, "scan")
    assert "Loaded 3Di state cache" in result.stdout
    scan = pd.read_csv(tmp_path / "csubst_scan/csubst_scan.tsv", sep="\t")
    if code == 1:
        assert not scan.empty
    if not scan.empty:
        assert set(scan.nonsyn_recode) == {"3di20"}
        assert (tmp_path / "csubst_scan/csubst_scan.tree_site.pdf").stat().st_size > 1000
    # Independent direct call uses the same public fit inputs, without the GG shell.
    common = ["--nonsyn_recode", "3di20", "--full_cds_alignment_file", str(full_dir / "csubst.fasta"),
              "--rooted_tree_file", str(full_dir / "csubst.nwk"), "--genetic_code", str(code),
              "--iqtree_model", codon_model, "--sa_state_cache", "yes", "--sa_state_cache_file", str(full_dir / "csubst_3di_state_cache.npz")]
    for suffix in ("treefile", "state", "rate", "iqtree", "log"):
        common += ["--iqtree_" + suffix, str(full_dir / ("csubst." + suffix))]
    run(["csubst", "search", *common, "--max_arity", "2", "--exhaustive_until", "2",
         "--cutoff_stat", "OCNany2spe,0", "--max_combination", "10000", "--fg_exclude_wg", "no",
         "--fg_stem_only", "no", "--outdir", "direct_search"], tmp_path, env, "direct_search")
    pd.testing.assert_frame_equal(pd.read_csv(tmp_path / "csubst_search/csubst_cb_2.tsv", sep="\t"),
                                  pd.read_csv(tmp_path / "direct_search/csubst_cb_2.tsv", sep="\t"))
    run(["csubst", "scan", *common, "--foreground", "foreground.tsv", "--fg_format", "2",
         "--scan_pvalue_calibration", "none", "--scan_n_permutations", "0", "--outdir", "direct_scan"], tmp_path, env, "direct_scan")
    pd.testing.assert_frame_equal(scan, pd.read_csv(tmp_path / "direct_scan/csubst_scan.tsv", sep="\t"))
    branches = str(scan.iloc[0].support_branch_ids) if not scan.empty else "1,2"
    for label in ("normal_sites", "candidate_sites"):
        directory = tmp_path / label
        directory.mkdir()
        command = sites.build_csubst_sites_command(str(base), str(base), branches, 2, "3di20", code, pdb="none")
        result = run(command, directory, env, label)
        assert "Loaded 3Di state cache" in result.stdout
        artifacts = sites.resolve_site_artifacts(str(directory), branches)
        assert artifacts["site_summary_pdf"] and Path(artifacts["site_summary_pdf"]).stat().st_size > 1000
        frame = pd.read_csv(artifacts["site_table_tsv"], sep="\t")
        assert set(frame.codon_site_alignment) == set(range(1, 101))
        if label == "normal_sites":
            reference = frame
        else:
            pd.testing.assert_frame_equal(reference, frame)
        states = directory / "states.fa"
        bundle.write_structural_tip_alignment(full_dir / "csubst.fasta", states)
        assert all(len(record.seq) == 100 for record in SeqIO.parse(states, "fasta"))
    assert (scan.codon_site_alignment == scan.site + 1).all()
    if code == 1:
        family, colors = prepare_report_family(tmp_path, base)
        # Synthetic tip IDs have no public structure accession. Keep real sites inference.
        real_builder = sites.build_csubst_sites_command
        def without_structure_lookup(*args, **kwargs):
            kwargs["pdb"] = "none"
            return real_builder(*args, **kwargs)
        monkeypatch.setattr(sites, "build_csubst_sites_command", without_structure_lookup)
        for key, value in env.items():
            monkeypatch.setenv(key, value)
        output = tmp_path / "normal_report"
        output.mkdir()
        result = sites.process_index("OG3DI", branches, str(output), str(family), str(colors), 2, "3di20", "Real 3Di validation")
        assert result[1] is None, result[1]
        assert list(output.rglob("summary.*.pdf"))
        selected = scan.iloc[:1].copy()
        selected["orthogroup"] = "OG3DI"
        selected["_canonical_support_branch_ids"] = selected.support_branch_ids.astype(str)
        record = candidates.assign_candidate_ids(selected, "3di20", "none").iloc[0].to_dict()
        record["_required_input_signature"] = "runtime-fixture"
        result = candidates.analyze_candidate(record, tmp_path / "candidate_report", str(family),
                                             {"trait": str(colors)}, "3di20", "none")
        assert result["status"] == "completed", result
        assert list((tmp_path / "candidate_report").rglob("*.focused_tree_site.pdf"))
    metadata = json.loads((base / "csubst.input.json").read_text())
    assert metadata["genetic_code"] == code


def prepare_report_family(root, base):
    from csubst import ete, parser_misc, tree
    family = root / "report_family"
    for name in ("iqtree_anc", "stat_branch", "clipkit", "mafft", "rpsblast"):
        (family / name).mkdir(parents=True)
    shutil.copyfile(root / "bundle.zip", family / "iqtree_anc/OG3DI_iqtree.anc.zip")
    shutil.copyfile(root / "trimmed.fa", family / "clipkit/OG3DI_cds.clipkit.fa")
    shutil.copyfile(root / "full.fa", family / "mafft/OG3DI_cds.aln.fa")
    context = parser_misc.annotate_tree(tree.read_treefile({"rooted_tree_file": str(base / "csubst.nwk"),
                                                          "iqtree_treefile": str(base / "csubst.treefile")}))
    def index(node):
        return int(ete.get_prop(node, "numerical_label")) if node is not None else -999
    rows = []
    for node in context["tree"].traverse():
        children = list(node.children)
        rows.append(dict(branch_id=index(node), numerical_label=index(node), node_name=node.name,
                         parent=index(node.up), child1=index(children[0]) if children else -999,
                         child2=index(children[1]) if children else -999, bl_rooted=node.dist,
                         support_unrooted=float("nan") if ete.is_root(node) else 100,
                         gene_labels="; ".join(ete.get_leaf_names(node)),
                         so_event="L" if ete.is_leaf(node) else "S", so_event_parent="S",
                         taxon=node.name if ete.is_leaf(node) else "", start=1, end=300,
                         chromosome="chr1", expression_trait=index(node)+1, num_intron=0))
    pd.DataFrame(rows).to_csv(family / "stat_branch/OG3DI_stat.branch.tsv", sep="\t", index=False)
    colors = family / "color.tsv"
    colors.write_text("species\tcolor\n" + "".join(f"{name}\t#336699\n" for name in "abcdefgh"))
    return family, colors
