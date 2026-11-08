"""Regressions at the shared helpers' real caller boundaries."""

import importlib.util
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

import pandas
import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"


def load(name):
    spec = importlib.util.spec_from_file_location(name, SUPPORT / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_numbered_tree_json_bridge_protects_its_request(tmp_path):
    request, response = tmp_path / 'request.json', tmp_path / 'response.json'
    payload = dict(operation='clade_keys', tree={
        'tip.label': ['A', 'B'], 'Nnode': 1, 'edge': [[3, 1], [3, 2]]})
    request.write_text(json.dumps(payload))
    command = [sys.executable, str(SUPPORT / 'transfer_tree_labels.py'), str(request)]
    failed = subprocess.run(command + [str(request)], capture_output=True, text=True)
    assert failed.returncode != 0 and 'overwrite' in failed.stderr
    assert json.loads(request.read_text()) == payload
    subprocess.run(command + [str(response)], check=True)
    values = json.loads(response.read_text())
    assert values[:2] == [None, None]
    assert values[2].startswith('clade-bits-v1:') and values[2].endswith(':3')


def test_catalog_snapshot_roundtrip_preserves_species_keys_and_optional_fields(tmp_path):
    snapshot = load("build_id_options_snapshot")
    manifest = load("build_download_manifest")
    entries = [
        {"id": "shared-source", "species": "Target species", "species_key": key,
         "cds_url": f"https://example.org/{key}.fa", "fernbase_confidence_mode": "high-low combined"}
        for key in ("Target_species_a", "Target_species_b")
    ]
    source = tmp_path / "snapshot.json"
    source.write_text(json.dumps({"providers": {"direct": entries}}))
    first = snapshot.load_snapshot_direct_catalog(source)
    source.write_text(json.dumps({"providers": {"direct": first}}))
    second = manifest.load_direct_catalog_snapshot(source)
    assert first == second
    assert [row["species_key"] for row in second] == ["Target_species_a", "Target_species_b"]
    assert all(row["fernbase_confidence_mode"] == "high-low combined" for row in second)
    assert [row["cds_url"] for row in second] == [entry["cds_url"] for entry in entries]
    with pytest.raises(ValueError, match="catalog choices must be unique"):
        manifest.write_manifest_xlsx([], tmp_path / "ambiguous.xlsx", direct_catalog_entries=second)


def test_catalog_reader_rejects_ambiguous_ids_but_accepts_full_choices():
    manifest = load("format_species_manifest")
    records = [{"id": "shared", "choice": f"shared ({name})", "species_key": name,
                "cds_url": f"https://example.org/{name}.fa"} for name in ("Species_a", "Species_b")]
    lookup = manifest.build_direct_catalog_lookup_from_xlsx(records)
    with pytest.raises(ValueError, match="Ambiguous direct catalog selection"):
        manifest.apply_direct_catalog_xlsx_defaults([{"provider": "direct", "id": "shared"}], lookup)
    rows = manifest.apply_direct_catalog_xlsx_defaults(
        [{"provider": "direct", "id": "shared (Species_b)"}], lookup)
    assert rows[0]["species_key"] == "Species_b"
    assert rows[0]["cds_url"] == "https://example.org/Species_b.fa"
    explicit = {"provider": "direct", "id": "shared", "cds_url": "https://example.org/explicit.fa"}
    assert manifest.apply_direct_catalog_xlsx_defaults([explicit], lookup) == [explicit]


@pytest.mark.parametrize("use_shared", [False, True])
def test_genetic_code_helpers_keep_group_subgroup_and_hybrid_overrides(tmp_path, use_shared):
    labels = ["Oryza_sativa_Indica_Group", "Oryza_sativa_tropical_japonica_subgroup",
              "Cenchrus_americanus_x_Cenchrus_purpureus"]
    cds = tmp_path / "cds"
    cds.mkdir()
    for label in labels:
        (cds / f"{label}.cds.fa.gz").write_bytes(b"filename-only fixture")
    table = tmp_path / "codes.tsv"
    table.write_text("species\tgenetic_code\n" + "".join(f"{label}\t4\n" for label in labels))
    output = tmp_path / "resolved.tsv"
    program = f"set -euo pipefail\nsource {shlex.quote(str(SUPPORT / 'gg_util/06_workspace_validation.sh'))}\n"
    function = "gg_prepare_species_genetic_code_table"
    if not use_shared:
        core = (SUPPORT.parent / "core/gg_genome_evolution_core.sh").read_text()
        program += core[core.index("species_genetic_code_table_path() {"):core.index("lookup_species_genetic_code() {")]
        function = "prepare_species_genetic_code_table"
    program += f'{function} "$1" 1 "$2" "$3"\n'
    result = subprocess.run(["bash", "-c", program, "bash", str(cds), str(output), str(table)],
                            env={**os.environ, "gg_support_dir": str(SUPPORT)}, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    actual = pandas.read_csv(output, sep="\t", dtype=str)
    assert actual["species"].tolist() == sorted(labels)
    assert actual["genetic_code"].tolist() == ["4"] * 3
    assert actual["source"].tolist() == ["species_genetic_code.tsv"] * 3
    assert not result.stderr


def test_summary_missing_directory_policies_are_preserved(tmp_path):
    og = load("orthogroup_output_summary")
    query = load("query2family_output_summary")
    source = pandas.DataFrame({"Total": [1]}, index=["OG000001"])
    missing = str(tmp_path / "missing")
    group = og.get_alignment_stats(source.copy(), missing, "original", 2)
    family = query.get_alignment_stats(source.copy(), missing, "original", ["OG000001"], 2)
    assert group.shape == (1, 9)
    assert group.iloc[:, 1:].isna().all().all()
    pandas.testing.assert_frame_equal(family, source)


@pytest.mark.parametrize("name", ["validate_cds_gff_mapping", "validate_longest_cds_selection"])
def test_validator_parallel_errors_preserve_species_and_result_order(monkeypatch, name):
    module = load(name)
    tasks = [{"index": i, "species_prefix": f"Species_{i}", "summary_row": {"species_prefix": f"Species_{i}"},
              "cds_file": Path(f"Species_{i}.fa"), "gff_file": Path(f"Species_{i}.gff")} for i in (2, 0, 1)]

    def validate(task, missing_limit):
        if task["index"] == 1:
            raise RuntimeError("executor failure")
        return {"index": task["index"], "ok": True}

    monkeypatch.setattr(module, "validate_single_species", validate)
    results = module.run_validation_tasks(tasks, 5, 2)
    assert [row["index"] for row in results] == [0, 1, 2]
    assert results[1]["ok"] is False
    assert "Species_1" in results[1]["error"] and "executor failure" in results[1]["error"]
    with pytest.raises(RuntimeError, match="executor failure"):
        module.run_validation_tasks(tasks, 5, 1)


def test_boolean_contracts_keep_unknown_table_values_and_reject_unknown_cli_values():
    values = load("boolean_values")
    assert values.parse_boolish("yes") is True
    assert values.parse_boolish("0") is False
    assert values.parse_boolish("unknown") is pandas.NA
    assert values.parse_boolish("on") is pandas.NA
    assert values.strtobool("on") == 1
    with pytest.raises(ValueError, match="invalid truth value"):
        values.strtobool("unknown")


def test_package_identity_is_fresh_and_includes_the_recorded_revision(tmp_path):
    package = tmp_path / "identity_fixture"
    package.mkdir()
    source = package / "__init__.py"
    source.write_text("value = 1\n")
    revisions = tmp_path / "revisions.tsv"
    revisions.write_text("source\trevision\nidentity_fixture\t" + "a" * 40 + "\n")
    command = [sys.executable, str(SUPPORT / "python_source_identity.py"), "identity_fixture", str(revisions)]
    environment = {**os.environ, "PYTHONPATH": str(tmp_path)}
    first = subprocess.check_output(command, env=environment, text=True).strip()
    source.write_text("value = 2\n")
    second = subprocess.check_output(command, env=environment, text=True).strip()
    assert first != second
    assert first.startswith("identity_fixture unpackaged source_sha256=")
    assert first.endswith(" " + "a" * 40) and second.endswith(" " + "a" * 40)
    resource = package / "coordinates.R"
    resource.write_text("value <- 1\n")
    third = subprocess.check_output(command, env=environment, text=True).strip()
    resource.write_text("value <- 2\n")
    fourth = subprocess.check_output(command, env=environment, text=True).strip()
    assert len({second, third, fourth}) == 3


@pytest.mark.parametrize("method", ["blastp", "mmseqs2"])
@pytest.mark.parametrize("fail", [False, True])
def test_shared_uniprot_search_preserves_flags_and_propagates_failures(tmp_path, method, fail):
    script = r'''
set -euo pipefail
source "$1/gg_util/10_reference_databases.sh"
gg_support_dir=$1
gg_memory_fraction_gb() { printf 3; }
ensure_uniprot_sprot_blast_db() { printf database; }
ensure_uniprot_sprot_mmseqs_db() { printf database; }
validate_uniprot_sprot_db_prefix() { [[ "$1" == database ]]; }
ensure_uniprot_sprot_metadata_tsv() { printf metadata.tsv; }
blastp() {
  printf '%s\n' "$*" >> search-args
  [[ "$fail" == 0 ]] || return 23
}
mmseqs() {
  printf '%s\n' "$*" >> search-args
  [[ "$fail" == 0 ]] || return 23
}
python() { printf '%s\n' "$*" > formatter-args; }
if gg_annotate_uniprot_proteins 'query with spaces.fa' "$PWD" "$method" 2 4 annotation.tsv; then
  [[ "$fail" == 0 ]]
else
  status=$?
  [[ "$fail" == 1 && "$status" == 23 ]]
fi
'''
    result = subprocess.run(["bash", "-c", script, "bash", str(SUPPORT)], cwd=tmp_path,
                            env={**os.environ, "method": method, "fail": str(int(fail))},
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    args = (tmp_path / "search-args").read_text()
    if fail:
        assert not (tmp_path / "formatter-args").exists()
        assert len(args.splitlines()) == 1
    else:
        formatter = (tmp_path / "formatter-args").read_text()
        assert "--query_fasta query with spaces.fa" in formatter
        assert "--outfile annotation.tsv" in formatter
        assert "--uniprot_meta_tsv metadata.tsv" in formatter
        if method == "blastp":
            assert "-num_threads 2" in args and "-max_target_seqs 1 -evalue 1e-2" in args
        else:
            assert "--split-memory-limit 3G --max-seqs 1 -e 1e-2 -s 7.5" in args
            assert "--format-output query,target,pident,alnlen,evalue,bits,qlen" in args


@pytest.mark.parametrize("storage", ["zip", "raw", "files"])
def test_shared_storage_configuration_preserves_normalization_and_arguments(tmp_path, storage):
    script = r'''
set -euo pipefail
source "$1/gg_util/01_runtime_config.sh"
gg_normalize_gene_family_storage_config 100
printf '%s\n' "$gene_family_output_storage" "${gene_family_archive_write_args[@]}"
'''
    result = subprocess.run(["bash", "-c", script, "bash", str(SUPPORT)],
                            env={**os.environ, "gene_family_output_storage": storage,
                                 "gene_family_zip_compression": "STORE", "gene_family_zip_workers": "2"},
                            text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.splitlines() == ["zip" if storage == "zip" else "files", "--compression", "store",
                                          "--compression-level", "6", "--workers", "2", "--max-final-zip-bytes", "0"]
