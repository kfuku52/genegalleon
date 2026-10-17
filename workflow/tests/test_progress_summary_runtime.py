"""Execute the real progress core in the GeneGalleon runtime, independent of cwd."""
import os
import subprocess
from pathlib import Path

import pandas


def test_progress_core_keeps_all_summary_tables_in_relocated_workspace(tmp_path):
    workspace = tmp_path / "relocated workspace"
    query = workspace / "input/query_gene"
    query.mkdir(parents=True)
    (query / "YABBY").write_text("gene1\n")
    family = workspace / "output/query2family/tree_plot"
    family.mkdir(parents=True)
    (family / "YABBY_tree_plot.pdf").write_text("plot\n")
    orthogroup = workspace / "output/orthogroup/cds_fasta"
    orthogroup.mkdir(parents=True)
    (orthogroup / "HOG0000010_cds.fa").write_text(">gene1\nATG\n")
    counts = workspace / "output/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv"
    counts.parent.mkdir(parents=True)
    counts.write_text("Orthogroup\tTotal\nHOG0000010\t1\n")
    (workspace / "input/species_rnaseq/Arabidopsis_thaliana").mkdir(parents=True)
    transcripts = workspace / "output/transcriptome_assembly/longestcds"
    transcripts.mkdir(parents=True)
    (transcripts / "Arabidopsis_thaliana.fa").write_text(">gene1\nATG\n")
    caller = tmp_path / "unrelated working directory"
    caller.mkdir()
    core = Path(__file__).resolve().parents[1] / "core/gg_progress_summary_core.sh"
    result = subprocess.run(["bash", str(core)], cwd=caller, text=True,
        capture_output=True, timeout=120, env={**os.environ,
            "gg_workspace_dir": str(workspace), "GG_COMMON_TMP_ROOT": "workspace",
            "gene_family_output_storage": "files", "ncpu_progress_summary": "1",
            "mode_transcriptome_assembly": "fastq"})
    assert result.returncode == 0, result.stdout + result.stderr
    expected = {
        "query2family_summary.tsv": "YABBY",
        "orthogroup_summary.tsv": "HOG0000010",
        "orthogroup_genecount.alignment_stats.tsv": "HOG0000010",
        "transcriptome_assembly_summary.tsv": "Arabidopsis_thaliana",
    }
    for filename, row_id in expected.items():
        assert row_id in pandas.read_csv(workspace / filename, sep="\t", index_col=0).index
        assert not (caller / filename).exists()
    query_summary = pandas.read_csv(workspace / "query2family_summary.tsv", sep="\t", index_col=0)
    assert query_summary.loc["YABBY", "tree_plot"] == 1
    assert counts.read_text() == "Orthogroup\tTotal\nHOG0000010\t1\n"

