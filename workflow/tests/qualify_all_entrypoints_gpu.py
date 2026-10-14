#!/usr/bin/env python3
"""Bounded qualification of all eight real entrypoints in one GPU allocation.

Use a registered executor with a sealed source/input inventory and a fresh root.
The checked-in seven-species inputs are copied/consumed read-only; all generated
data and logs belong to the qualification. Optional database-heavy analyses are
disabled explicitly. Main enabled producers use their ordinary scientific
parameters, and success requires their outputs as well as wrapper exit zero.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
from pathlib import Path

from qualify_annevo_entrypoint_gpu import invoke, qualify, require, sha, write_json

PREFIXES = {
    "input_generation": "GG_INPUT_", "transcriptome_generation": "GG_TRANSCRIPTOME_",
    "fractionation_bias": "GG_FRACTIONATION_BIAS_", "genome_evolution": "GG_GENOME_EVOLUTION_",
    "gene_evolution": "GG_GENE_EVOLUTION_", "gene_summary": "GG_GENE_SUMMARY_",
    "progress_summary": "GG_PROGRESS_SUMMARY_",
}


def registered(source, stage):
    return subprocess.check_output(["bash", "-c", 'source "$1"; gg_print_entrypoint_config_vars "$2"',
        "config", str(source / "workflow/support/gg_entrypoint_config_vars.sh"),
        "gg_" + stage + "_entrypoint.sh"], text=True).splitlines()


def run_case(manifest, stage, workspace, settings, outputs):
    root, source = Path(manifest["root"]), Path(manifest["source"]["root"])
    names = registered(source, stage)
    options = {name: "0" for name in names if name.startswith("run_")}
    options.update(settings)
    require(set(options) <= set(names) | {"delete_tmp_dir"}, "unregistered qualification option: " + stage)
    env = {k: v for k, v in os.environ.items() if k.startswith("SLURM_") or k in
        {"HOME", "USER", "LOGNAME", "LANG", "LC_ALL", "CUDA_VISIBLE_DEVICES"}}
    env.update(PATH="/usr/local/bin:/usr/bin:/bin", PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1",
        gg_workspace_dir=str(workspace), gg_container_image_path=manifest["image"],
        GG_CONTAINER_PROJECT_ROOT_BIND=str(root) + ":" + str(root), GG_COMMON_TMP_ROOT="workspace",
        GG_COMMON_REFERENCE_SPECIES="Arabidopsis_thaliana", GG_COMMON_GENETIC_CODE="1",
        GG_COMMON_BUSCO_LINEAGE="embryophyta_odb12", SLURM_ARRAY_TASK_ID="1",
        APPTAINER_BINDPATH=str(source) + ":" + str(source) + ":ro",
        SINGULARITY_BINDPATH=str(source) + ":" + str(source) + ":ro")
    env.update({PREFIXES[stage] + name.upper(): str(value) for name, value in options.items()})
    entrypoint = source / "workflow" / ("gg_" + stage + "_entrypoint.sh")
    invoke(root, "all_" + stage, ["bash", str(entrypoint)], env, timeout=3600)
    require(": Ending" in (root / "logs" / ("all_" + stage + ".out")).read_text(), "wrapper did not end: " + stage)
    paths = [workspace / p for p in outputs]
    require(all(p.is_file() and p.stat().st_size > 0 for p in paths), "required main-stage output missing: " + stage)
    require(any(workspace.glob("output/versions/container.*.versions.log")), "version report missing: " + stage)
    return {"entrypoint": entrypoint.name, "entrypoint_sha256": sha(entrypoint), "workspace": str(workspace),
        "settings": options, "outputs": {str(p): sha(p) for p in paths}, "verified": True}


def all_entrypoints(manifest):
    root, source = Path(manifest["root"]), Path(manifest["source"]["root"])
    require(os.environ.get("SLURM_CPUS_PER_TASK") == "4" and os.environ.get("SLURM_MEM_PER_NODE") == "32768"
        and os.environ.get("CUDA_VISIBLE_DEVICES"), "allocated 4 CPU / 32 GiB / CUDA required")
    seed = source / "workspace/input"
    require(seed.is_dir(), "sealed checked-in qualification inputs missing")
    inventory = json.loads(Path(manifest["source"]["inventory"]).read_text())["files"]
    seed_files = {name: digest for name, digest in inventory.items() if name.startswith("workspace/input/")}
    require(seed_files and all(sha(source / name) == digest for name, digest in seed_files.items()), "input seed seal changed")
    workspace = root / "comparative_workspace"
    workspace.mkdir(mode=0o700)
    inputs = workspace / "input"
    inputs.mkdir()
    # Entrypoint snapshots have their own project mount. Keep local download
    # sources inside that reviewed mount rather than relying on inherited
    # Apptainer bind environment variables outside the project.
    assets = root / "input_assets"
    assets.mkdir(mode=0o700)
    for role in ("species_cds", "species_gff", "species_genome"):
        shutil.copytree(seed / role, assets / role)
        for path in (assets / role).iterdir():
            require(sha(path) == seed_files["workspace/input/" + role + "/" + path.name], "copied local input changed")
            path.chmod(0o400)
    query = inputs / "query_gene"
    query.mkdir()
    shutil.copyfile(seed / "query_gene/YABBY", query / "YABBY")
    shutil.copytree(seed / "species_rnaseq/Amborella_trichopoda", inputs / "species_rnaseq/Amborella_trichopoda")
    fixtures = json.loads(Path(manifest["model_qualification"]).read_text())
    taxonomy = workspace / "downloads/ete_taxonomy/taxa.sqlite"
    taxonomy.parent.mkdir(parents=True)
    shutil.copyfile(fixtures["taxonomy_db"], taxonomy)
    taxonomy.chmod(0o400)
    require(sha(taxonomy) == fixtures["taxonomy_sha256"], "copied taxonomy changed")
    manifest_path = inputs / "input_generation/local.tsv"
    manifest_path.parent.mkdir()
    cds_files = sorted((assets / "species_cds").glob("*.fa"))
    require(len(cds_files) == 7, "exact seven-species biological fixture required")
    species = []
    with manifest_path.open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["provider", "id", "species_key", "local_cds_path", "local_gff_path", "local_genome_path",
                         "cds_filename", "gff_filename", "genome_filename"])
        for cds in cds_files:
            name = "_".join(cds.name.split("_")[:2])
            gff = list((assets / "species_gff").glob(name + "*"))
            genome = list((assets / "species_genome").glob(name + "*"))
            require(len(gff) == len(genome) == 1, "paired input missing: " + name)
            species.append(name)
            # Local discovery uses role names; a source genome named simply
            # <species>.fa must not be classified as provided CDS.
            writer.writerow(["local", name, name, str(cds), str(gff[0]), str(genome[0]),
                             name + ".cds.fa", name + ".gff3", name + ".genome.fa"])
    results = [run_case(manifest, "input_generation", workspace, {
        "provider": "local", "download_manifest": str(manifest_path), "strict": "1",
        "require_cds": "1", "require_gff": "1", "require_genome": "1",
        "species_cds_dir": str(inputs / "species_cds"), "species_gff_dir": str(inputs / "species_gff"),
        "species_genome_dir": str(inputs / "species_genome"), "run_format_inputs": "1",
        "run_validate_inputs": "1", "run_cds_fx2tab": "1", "delete_tmp_dir": "0"},
        ["output/input_generation/gg_input_generation_species.tsv"])]
    require(len(list((inputs / "species_cds").glob("*.fa*"))) >= 7, "formatted species CDS missing")
    results.append(run_case(manifest, "transcriptome_generation", workspace, {
        "mode_transcriptome_assembly": "fastq", "assembly_method": "rnaSPAdes", "amalgkit_rrna_filter": "no",
        "run_amalgkit_metadata_or_integrate": "1", "run_amalgkit_getfastq": "1", "run_assembly": "1",
        "run_longestcds": "1", "run_longestcds_fx2tab": "1", "run_assembly_stat": "1",
        "run_amalgkit_quant": "1", "run_amalgkit_merge": "1", "delete_tmp_dir": "0"},
        ["output/transcriptome_assembly/assembled_transcripts_with_isoforms/Amborella_trichopoda_isoform.fa.gz",
         "output/transcriptome_assembly/longest_cds/Amborella_trichopoda_longestCDS.fa.gz",
         "output/transcriptome_assembly/assembly_stat/Amborella_trichopoda_assembly_stat.tsv"]))
    fractionation = root / "fractionation_workspace"
    fractionation.mkdir(mode=0o700)
    fixture_code = (
        "import importlib.util;from pathlib import Path;"
        "s=importlib.util.spec_from_file_location('fixture'," + repr(str(source / "workflow/tests/test_fractionation_bias_integration.py")) + ");"
        "m=importlib.util.module_from_spec(s);s.loader.exec_module(m);p=Path(" + repr(str(fractionation / "input")) + ");"
        "m._write_genome(p,'Target_species','t',('target_chr',));"
        "m._write_genome(p,'Query_species','q',('query_a','query_b'));"
        "(p/'fractionation_bias_pairs.tsv').write_text(" +
        repr("analysis_id\ttarget_species\tquery_species\tquota\twindow_size\nsynthetic\tTarget_species\tQuery_species\t1:2\t4\n") + ")")
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    invoke(root, "all_fractionation_fixture", [manifest["container"], "exec", "--containall", "--no-home", "--bind",
        str(root) + ":" + str(root) + "," + str(source) + ":" + str(source) + ":ro", manifest["image"],
        "python", "-B", "-c", fixture_code], env, timeout=60)
    results.append(run_case(manifest, "fractionation_bias", fractionation, {"run_kffractbias": "1", "delete_tmp_dir": "0"},
        ["output/kffractbias/synthetic/synthetic.genes.tsv", "output/kffractbias/synthetic/synthetic.summary.json",
         "output/kffractbias/synthetic/synthetic.synteny.zip"]))
    results.append(run_case(manifest, "genome_evolution", workspace, {"genome_evolution_mode": "orthogroups",
        "run_cds_translation": "1", "run_orthofinder": "1", "run_og_selection": "1", "delete_tmp_dir": "0"},
        ["output/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv"]))
    results.append(run_case(manifest, "gene_evolution", workspace, {"mode_gene_evolution": "query2family",
        "gene_family_output_storage": "files",
        "run_extract_query_fasta": "1", "run_query_blast": "1", "run_extract_primary_fasta": "1", "run_mafft": "1",
        "run_alignment_stats_original": "1", "run_clipkit": "1", "run_alignment_stats_cleaned": "1",
        "run_iqtree": "1", "run_tree_root": "1", "run_summary": "1", "run_tree_plot": "1", "delete_tmp_dir": "0"},
        ["output/query2family/iqtree_tree/YABBY_iqtree.nwk"]))
    results.append(run_case(manifest, "gene_summary", workspace, {"gene_family_source": "query2family",
        "run_family_completion_summary": "1", "run_presence_absence_summary": "1"},
        ["output/gene_summary/query2family/query2family_summary.tsv"]))
    results.append(run_case(manifest, "progress_summary", workspace, {"mode_transcriptome_assembly": "fastq",
        "gene_family_output_storage": "files"},
        ["output/query2family/query2family_summary.tsv"]))
    require(all(sha(source / name) == digest for name, digest in seed_files.items()), "sealed biological inputs changed")
    native = qualify(manifest)
    results.append({"entrypoint": "gg_genome_annotation_entrypoint.sh", "verified": True, "profiles_passed": 6,
        "entrypoint_runs_passed": 12, "evidence": "receipts/entrypoint_qualification.json"})
    expected = {"gg_" + stage + "_entrypoint.sh" for stage in PREFIXES} | {"gg_genome_annotation_entrypoint.sh"}
    require({row["entrypoint"] for row in results} == expected and all(row["verified"] for row in results), "all eight entrypoints required")
    proof = {"all_entrypoints_verified": True, "entrypoints_passed": 8, "entrypoint_runs_passed": 19,
        "profiles_passed": native["profiles_passed"], "cases": results, "fixture_species": species,
        "scientific_completion_verified": False, "scope": "eight-real-entrypoints-with-enabled-main-producers-and-six-GPU-profiles"}
    write_json(root / "receipts/all_entrypoints_qualification.json", proof)
    return proof


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    all_entrypoints(json.loads(args.manifest.read_text()))
