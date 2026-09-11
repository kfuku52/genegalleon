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
    output_dir = tmp_path / f"amas_{stage}"
    output_dir.mkdir()
    output = output_dir / f"{family}_amas.{stage}.tsv"
    core = (ROOT / "workflow/core/gg_gene_evolution_core.sh").read_text()
    start = core.index(f'task="cdskit stats for {stage} alignment"')
    end = core.index('\nfi', start) + len('\nfi')
    block = core[start:end]
    harness = f"""set -euo pipefail
cd {shlex.quote(str(tmp_path))}
og_id={family}
GG_TASK_CPUS=1
amas_data_type={seq_type}
dir_output_active={shlex.quote(str(tmp_path))}
gg_workspace_dir=$dir_output_active
file_og_untrimmed_aln_analysis={shlex.quote(str(source))}
file_og_trimmed_aln_analysis=$file_og_untrimmed_aln_analysis
file_og_amas_{stage}={shlex.quote(str(output))}
run_amas_{stage}=1
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
    assert table.loc[0, "Alignment_name"] == f"{family}.amas.{stage}.input.fasta"
    assert not (tmp_path / f"{family}.amas.{stage}.input.fasta").exists()
    assert "statistics_engine=cdskit-stats-alignment" in (tmp_path / "recorded.txt").read_text()
    for reader in ["orthogroup_output_summary", "query2family_output_summary"]:
        spec = importlib.util.spec_from_file_location(reader, ROOT / f"workflow/support/{reader}.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        frame = pandas.DataFrame({"Total": [2]}, index=[family])
        kwargs = {"query_id_matchers": [family]} if reader == "query2family_output_summary" else {}
        result = module.get_amas_stats(frame, str(output_dir), stage, ncpu=1, **kwargs)
        assert result.loc[family, f"No_of_taxa_{stage}"] == 2
        gc = result.loc[family, f"GC_content_{stage}"]
        assert pandas.isna(gc) if seq_type == "aa" else gc == 0.5
