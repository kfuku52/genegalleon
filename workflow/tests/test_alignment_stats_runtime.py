"""Run the actual gene-evolution statistics blocks and their summary consumers."""
import gzip
import importlib.util
import shlex
import subprocess
from pathlib import Path

import pandas
import pytest

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize("stage", ["original", "cleaned"])
@pytest.mark.parametrize("seq_type", ["dna", "aa"])
def test_alignment_stats_stage_and_summary_readers(tmp_path, stage, seq_type):
    family = "HOG0000010"
    source = tmp_path / "alignment.fa.gz"
    sequence = "ACGTNN--" if seq_type == "dna" else "MPEPX.--"
    with gzip.open(source, "wt") as handle:
        handle.write(f">a\n{sequence}\n>b\n{sequence}\n")
    output_dir = tmp_path / f"alignment_stats_{stage}"
    output_dir.mkdir()
    output = output_dir / f"{family}_alignment_stats.{stage}.tsv"
    core = (ROOT / "workflow/core/gg_gene_evolution_core.sh").read_text()
    start = core.index(f'task="cdskit stats for {stage} alignment"')
    end = core.index('\nfi', start) + len('\nfi')
    block = core[start:end]
    harness = f"""set -euo pipefail
cd {shlex.quote(str(tmp_path))}
og_id={family}
GG_TASK_CPUS=1
alignment_stats_data_type={seq_type}
dir_output_active={shlex.quote(str(tmp_path))}
gg_workspace_dir=$dir_output_active
file_og_untrimmed_aln_analysis={shlex.quote(str(source))}
file_og_trimmed_aln_analysis=$file_og_untrimmed_aln_analysis
file_og_alignment_stats_{stage}={shlex.quote(str(output))}
run_alignment_stats_{stage}=1
disable_if_no_input_file() {{ :; }}
gg_artifact_prepare_stage() {{ printf -v "$1" '%s' 1; }}
gg_step_start() {{ :; }}
gg_step_skip() {{ :; }}
mv_out() {{ mv -- "$1" "$2"; }}
gg_artifact_record() {{ printf '%s\\n' "$@" > recorded.txt; }}
{block}
"""
    subprocess.run(["bash", "-c", harness], check=True, capture_output=True, text=True)
    table = pandas.read_csv(output, sep="\t")
    assert table.loc[0, "No_of_taxa"] == 2
    assert table.loc[0, "Alignment_length"] == 8
    assert table.loc[0, "Missing_percent"] == 50
    assert table.loc[0, "No_variable_sites"] == 0
    assert table.loc[0, "Alignment_name"] == f"{family}.alignment_stats.{stage}.input.fasta"
    assert not (tmp_path / f"{family}.alignment_stats.{stage}.input.fasta").exists()
    assert "statistics_engine=cdskit-stats-alignment" in (tmp_path / "recorded.txt").read_text()
    for reader in ["orthogroup_output_summary", "query2family_output_summary"]:
        spec = importlib.util.spec_from_file_location(reader, ROOT / f"workflow/support/{reader}.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        frame = pandas.DataFrame({"Total": [2]}, index=[family])
        kwargs = {"query_id_matchers": [family]} if reader == "query2family_output_summary" else {}
        result = module.get_alignment_stats(frame, str(output_dir), stage, ncpu=1, **kwargs)
        assert result.loc[family, f"No_of_taxa_{stage}"] == 2
        gc = result.loc[family, f"GC_content_{stage}"]
        assert pandas.isna(gc) if seq_type == "aa" else gc == 0.5


@pytest.mark.parametrize("mode", ["orthogroup", "query2family"])
@pytest.mark.parametrize("archived", [False, True])
@pytest.mark.parametrize("layout", ["underscore", "dot"])
@pytest.mark.parametrize("seq_type", ["dna", "aa"])
def test_one_time_statistics_migration(tmp_path, archived, layout, seq_type, mode):
    import json
    import sys
    import zipfile

    from workflow.support.gene_family_output_store import (
        GeneFamilyOutputStore,
        convert_storage_to_zip,
        family_context,
    )

    root = tmp_path / mode
    family = "HOG0000010" if mode == "orthogroup" else "query_alpha"
    genecount = tmp_path / "counts.tsv"
    genecount.write_text(f"Orthogroup\tTotal\n{family}\t2\n")
    sequence = "ACGTNN--" if seq_type == "dna" else "MPEPX.--"
    source = root / "mafft" / f"{family}_cds.aln.fa.gz"
    source.parent.mkdir(parents=True)
    with gzip.open(source, "wt") as handle:
        handle.write(f">a\n{sequence}\n>b\n{sequence}\n")
    query_dir = tmp_path / "query_gene"
    query_dir.mkdir()
    (query_dir / family).write_text(">query\nACGT\n")
    original_bytes = source.read_bytes()
    for stage in ["original", "cleaned"]:
        subdir = f"amas_{stage}" if layout == "underscore" else f"amas.{stage}"
        prefix = f"{family}_" if layout == "underscore" else f"{family}."
        old = root / subdir / f"{prefix}amas.{stage}.tsv"
        old.parent.mkdir()
        old.write_text("No_of_taxa\n99\n")
        manifest = root / "artifact_provenance" / f"{family}.amas_{stage}.json"
        manifest.parent.mkdir(exist_ok=True)
        manifest.write_text(json.dumps({
            "family_id": family, "step": f"amas_{stage}",
            "inputs": [{"label": "alignment", "scope": "logical", "path": str(source.relative_to(root))}],
            "parameters": {"data_type": seq_type},
        }))
    if archived:
        ids, identify = family_context(mode, genecount=genecount, query_dir=query_dir)
        convert_storage_to_zip(root, mode, ids, identify)
        assert not source.exists()
    summary = tmp_path / "summary.tsv"
    old_aggregate = tmp_path / "orthogroup_genecount.amas.tsv"
    if mode == "orthogroup":
        old_aggregate.write_text("obsolete\n")
    command = [sys.executable, str(ROOT / "workflow/migrations/migrate_alignment_statistics.py"),
               "--root", str(root), "--workspace-root", str(tmp_path), "--mode", mode,
               "--genecount", str(genecount), "--dir-query-gene", str(query_dir), "--summary-out", str(summary)]
    if layout == "dot":
        mapping = tmp_path / "mapping.tsv"
        mapping.write_text("family_id\tstage\talignment\tseq_type\n"
                           f"{family}\toriginal\t{source.relative_to(root)}\t{seq_type}\n")
        command += ["--alignments", str(mapping)]
    subprocess.run(command, check=True, text=True, timeout=40)
    store = GeneFamilyOutputStore(root)
    if mode == "orthogroup":
        assert not old_aggregate.exists()
    for stage in ["original", "cleaned"]:
        assert store.logical_exists(f"alignment_stats_{stage}/{family}_alignment_stats.{stage}.tsv")
        assert not store.logical_exists(f"artifact_provenance/{family}.amas_{stage}.json")
    out = pandas.read_csv(summary, sep="\t")
    assert out.loc[0, "No_of_taxa_original"] == 2
    assert out.loc[0, "No_of_taxa_clean"] == 2
    assert not any("amas" in col for col in out.columns)
    with store.open_binary("mafft", source.name) as handle:
        assert handle.read() == original_bytes
    for archive in root.rglob("*.zip"):
        with zipfile.ZipFile(archive) as handle:
            assert not any("amas" in name for name in handle.namelist())
    subprocess.run(command, check=True, text=True, timeout=40)


def test_migration_preserves_legacy_outputs_when_alignment_is_invalid(tmp_path):
    import sys

    root = tmp_path / "orthogroup"
    old = root / "amas_original" / "HOG0000010_amas.original.tsv"
    old.parent.mkdir(parents=True)
    old.write_text("No_of_taxa\n99\n")
    source = root / "invalid.fa"
    source.write_text(">a\nACGT\n>b\nA\n")
    mapping = tmp_path / "mapping.tsv"
    mapping.write_text("family_id\tstage\talignment\tseq_type\nHOG0000010\toriginal\tinvalid.fa\tdna\n")
    genecount = tmp_path / "counts.tsv"
    genecount.write_text("Orthogroup\tTotal\nHOG0000010\t2\n")
    result = subprocess.run([
        sys.executable, str(ROOT / "workflow/migrations/migrate_alignment_statistics.py"),
        "--root", str(root), "--workspace-root", str(tmp_path), "--mode", "orthogroup",
        "--genecount", str(genecount), "--summary-out", str(tmp_path / "summary.tsv"),
        "--alignments", str(mapping),
    ], capture_output=True, text=True, timeout=40)
    assert result.returncode != 0
    assert old.read_text() == "No_of_taxa\n99\n"
    assert not (root / "alignment_stats_original").exists()
