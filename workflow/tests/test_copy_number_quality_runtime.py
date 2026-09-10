"""Execute the actual core quality stage and standalone plotting CLI in GG runtime."""

import csv
import hashlib
import json
import os
import shlex
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow/support"


def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fixture_inputs(tmp_path):
    tree = tmp_path / "dated tree.nwk"
    tree.write_text("(((sp1:1,sp2:1):1,(sp3:1,sp4:1):1):1,((sp5:1,sp6:1):1,(sp7:1,sp8:1):1):1);")
    traits = tmp_path / "traits.tsv"
    traits.write_text("species\theight\thabitat\n" + "".join(f"sp{i}\t{i}\tforest\n" for i in range(1, 9)))
    Path(str(traits) + ".schema.json").write_text(json.dumps({
        "schema_version": 1, "table_sha256": hashlib.sha256(traits.read_bytes()).hexdigest(),
        "traits": {"height": "numeric", "habitat": "categorical"},
    }))
    copies = tmp_path / "copies.tsv"
    copies.write_text("Description\tOrthogroup\t" + "\t".join(f"sp{i}" for i in range(1, 9)) + "\n"
                      + "test\tOG1\t1\t2\t3\t4\t5\t6\t7\t8\n")
    short = tmp_path / "short summaries"
    short.mkdir()
    for i in range(1, 9):
        (short / f"sp{i}.busco.short.txt").write_text(
            "# BUSCO version is: 6.0.0\n# The lineage dataset is: eukaryota_odb12\n"
            f"# BUSCO was run in mode: proteins\nC:{80+i}%[S:{79+i}%,D:1%],F:1%,M:{19-i}%,n:100\n")
    return tree, traits, copies, short


def test_matrix_only_cli_includes_busco_and_pairwise_n(tmp_path):
    tree, traits, _, short = fixture_inputs(tmp_path)
    output = tmp_path / "matrix only"
    result = subprocess.run(["Rscript", str(SUPPORT / "copy_number_quality_diagnostics.r"),
                             f"--file_sptree={tree}", f"--file_trait={traits}",
                             f"--busco_short_dir={short}", f"--outdir={output}"], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    rows = read_table(output / "trait_correlations.tsv")
    pair = next(r for r in rows if r["trait_x"] == "height" and r["trait_y"] == "busco_complete_pct")
    assert float(pair["correlation"]) == 1
    assert pair["n_species"] == "8"
    assert (output / "trait_correlations.pdf").read_bytes().startswith(b"%PDF")
    assert not (output / "family_busco_associations.tsv").exists()
    assert read_table(output / "excluded_traits.tsv")[0]["trait"] == "habitat"


def test_core_stage_calls_real_diagnostics_and_records_all_inputs(tmp_path):
    tree, traits, copies, short = fixture_inputs(tmp_path)
    cafe = tmp_path / "cafe"
    cafe.mkdir()
    (cafe / "Gamma_family_results.txt").write_text("#FamilyID\tpvalue\tSignificant at 0.05\nOG1\t0.001\ty\n")
    core = (ROOT / "workflow/core/gg_genome_evolution_core.sh").read_text()
    stage = core.split('task="Copy-number BUSCO quality diagnostics and trait correlations"', 1)[1]
    stage = stage.split('task="GO enrichment analysis"', 1)[0]
    env = dict(os.environ)
    env.update({
        "dir_orthogroup_copy_number": str(tmp_path / "copy output"),
        "genome_evolution_provenance_dir": str(tmp_path / "provenance"),
        "file_orthogroup_copy_number": str(copies), "file_dated_species_tree": str(tree),
        "gg_support_dir": str(SUPPORT), "genome_nwkit_identity": "runtime-test",
        "copy_number_quality_busco_table": "auto", "orthogroup_copy_number_trait": "all",
        "orthogroup_copy_number_trait_min_species": "4", "orthogroup_copy_number_trait_response_families": "", "orthogroup_copy_number_trait_family_ids": "",
        "orthogroup_copy_number_trait_family_file": "", "orthogroup_copy_number_trait_max_families": "all",
        "orthogroup_copy_number_trait_alpha": "0.05", "copy_number_quality_high_completeness": "95",
        "copy_number_quality_correlation_method": "spearman", "copy_number_quality_sensitivity": "0",
        "dir_species_busco_short": str(short), "gg_workspace_dir": str(tmp_path), "file_trait": str(traits),
        "dir_cafe_output": str(cafe), "run_cafe": "1", "run_copy_number_quality_diagnostics": "1", "task": "quality",
        "run_orthogroup_copy_number_trait_pgls": "0",
    })
    audit = tmp_path / "artifact-args.txt"
    # Scheduling is isolated; the helper, parser, NWKIT fits and publication are real.
    prelude = """set -euo pipefail
gg_artifact_contract_init() { :; }
gg_artifact_prepare_stage() { printf -v "$1" '%s' 1; }
gg_step_start() { :; }
gg_step_skip() { :; }
"""
    # Use the actual optional-input function to test directory provenance wiring.
    validation = ROOT / "workflow/support/gg_util/06_workspace_validation.sh"
    prelude += f"source {shlex.quote(str(validation))}\n"
    # Sourcing the library defines unrelated functions; scheduler stubs remain authoritative.
    prelude += 'gg_artifact_contract_init() { :; }\ngg_artifact_prepare_stage() { printf -v "$1" \'%s\' 1; }\n'
    prelude += f"gg_artifact_record() {{ printf '%s\\n' \"$@\" > {shlex.quote(str(audit))}; }}\n"
    result = subprocess.run(["bash", "-c", prelude + stage], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    output = tmp_path / "copy output/quality_diagnostics"
    row = read_table(output / "cafe_family_quality.tsv")[0]
    assert row["#FamilyID"] == "OG1" and row["Significant at 0.05"] == "y"
    assert row["busco_quality_flag"] in {"quality_associated", "no_quality_association_detected"}
    contract = audit.read_text()
    assert f"busco_short_directory={short}" in contract
    assert f"cafe_family_results={cafe}/Gamma_family_results.txt" in contract
    assert "--output-logical-directory" in contract
    assert "nwkit_identity=runtime-test" in contract
    assert f"trait_schema={traits}.schema.json" in contract
    assert "trait_contract=" in contract


def test_declared_numeric_codes_and_observation_metadata_are_respected(tmp_path):
    tree, _, _, short = fixture_inputs(tmp_path)
    traits = tmp_path / "typed.tsv"
    traits.write_text("species\theight\thabitat\tgbif_observed_latitude_mean\n" +
                      "".join(f"sp{i}\t{i}\t{i % 3}\t{i + 20}\n" for i in range(1, 9)))
    schema = Path(str(traits) + ".schema.json")
    schema.write_text(json.dumps({"schema_version": 1,
        "table_sha256": hashlib.sha256(traits.read_bytes()).hexdigest(),
        "traits": {"height": "numeric", "habitat": "categorical", "gbif_observed_latitude_mean": "numeric"}}))
    output = tmp_path / "out"
    command = ["Rscript", str(SUPPORT / "copy_number_quality_diagnostics.r"),
               f"--file_sptree={tree}", f"--file_trait={traits}",
               f"--busco_short_dir={short}", f"--outdir={output}"]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert {r["trait_x"] for r in read_table(output / "trait_correlations.tsv")} == {"height", "busco_complete_pct"}
    excluded = {r["trait"]: r["reason"] for r in read_table(output / "excluded_traits.tsv")}
    assert excluded == {"habitat": "declared_categorical", "gbif_observed_latitude_mean": "observation_contract_excluded"}
    before = (output / "manifest.tsv").read_bytes()
    for trait, expected in [("habitat", "declared categorical"),
                            ("gbif_observed_latitude_mean", "hash-bound")]:
        failed = subprocess.run(command + [f"--trait={trait}"], capture_output=True, text=True)
        assert failed.returncode != 0 and expected in failed.stderr
        assert (output / "manifest.tsv").read_bytes() == before
    schema.write_text(schema.read_text().replace(hashlib.sha256(traits.read_bytes()).hexdigest(), "0" * 64))
    failed = subprocess.run(command, capture_output=True, text=True)
    assert failed.returncode != 0 and "does not match" in failed.stderr


def test_explicit_gbif_selection_survives_and_masks_incomplete_observations(tmp_path):
    tree, _, _, short = fixture_inputs(tmp_path)
    metric = "gbif_observed_northern_limit_lat"
    traits = tmp_path / "observed.tsv"
    traits.write_text(f"species\t{metric}\n" + "".join(f"sp{i}\t{i}\n" for i in range(1, 9)))
    metadata = {"schema_version": 1, "table_sha256": hashlib.sha256(traits.read_bytes()).hexdigest(),
        "traits": {metric: {"role": "observation", "source": "gbif", "source_column": metric}},
        "gbif": {"quality": [{"species": f"sp{i}", "status": "complete_download", "analysis_eligible": i != 8}
                              for i in range(1, 9)]}}
    # Use the canonical completion value from the owning acquisition contract.
    import sys
    sys.path.insert(0, str(SUPPORT))
    from gbif_observations import COMPLETE_STATES
    for row in metadata["gbif"]["quality"]:
        row["status"] = sorted(COMPLETE_STATES)[0]
    Path(str(traits) + ".metadata.json").write_text(json.dumps(metadata))
    output = tmp_path / "out"
    result = subprocess.run(["Rscript", str(SUPPORT / "copy_number_quality_diagnostics.r"),
        f"--file_sptree={tree}", f"--file_trait={traits}", f"--trait={metric}",
        f"--busco_short_dir={short}", f"--outdir={output}"], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    pair = next(r for r in read_table(output / "trait_correlations.tsv") if r["trait_x"] == metric and r["trait_y"] == "busco_complete_pct")
    assert pair["n_species"] == "7" and float(pair["correlation"]) == 1
    assert "GBIF values describe retained occurrence records" in (output / "trait_correlations.svg").read_text()
    assert json.loads((output / "species_trait_input.json").read_text())["masked_species"][metric][0]["species"] == "sp8"
