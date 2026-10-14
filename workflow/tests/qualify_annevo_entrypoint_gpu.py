#!/usr/bin/env python3
"""Bounded, allocated six-profile qualification of the real annotation entrypoint.

The site's registered executor verifies the sealed package and allocation before
calling this harness. Run on the allocated host, so the entrypoint itself starts
the SIF. Outputs are created only in a fresh package child; existing runs stay put.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

PROFILES = {"Mammalia", "Aves", "Actinopteri", "Insecta", "Fungi", "Magnoliopsida"}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, sort_keys=True, indent=2)


def invoke(root, label, argv, env, *, timeout=5400):
    started = time.monotonic()
    print("Starting " + label, flush=True)
    with (root / "logs" / (label + ".out")).open("x") as out, (root / "logs" / (label + ".err")).open("x") as err:
        result = subprocess.run(argv, cwd=root, env=env, stdout=out, stderr=err, timeout=timeout)
    write_json(root / "receipts" / (label + ".exit.json"), {
        "argv": argv, "returncode": result.returncode, "wall_seconds": time.monotonic() - started})
    require(result.returncode == 0, label + " failed with exit " + str(result.returncode))
    print("Passed " + label, flush=True)


def qualify(manifest):
    root, source = Path(manifest["root"]), Path(manifest["source"]["root"])
    fixtures = json.loads(Path(manifest["model_qualification"]).read_text())
    require(set(fixtures["fixtures"]) == PROFILES, "exact six profiles required")
    require(os.environ.get("SLURM_CPUS_PER_TASK") == "4" and os.environ.get("SLURM_MEM_PER_NODE") == "32768"
            and os.environ.get("CUDA_VISIBLE_DEVICES"), "allocated 4 CPU / 32 GiB / CUDA required")
    workspace, assets = root / "entrypoint_workspace", root / "entrypoint_assets"
    cache = workspace / "downloads/annotation/annevo"
    # The qualified conda-pack payload has 255-byte binary prefix fields. Its
    # final prefix is fixed before unpacking; installed generations cannot move.
    prefix = cache / "runtimes" / ("a" * 32) / ("g-" + "a" * 16) / "environment"
    require(len(os.fsencode(prefix)) <= 255,
            "qualification output path exceeds the qualified conda-pack 255-byte prefix limit; choose a shorter fresh package root")
    workspace.mkdir(mode=0o700)
    assets.mkdir(mode=0o700)
    genomes = workspace / "input/species_genome"
    genomes.mkdir(parents=True, mode=0o700)
    database = assets / "taxa.sqlite"
    shutil.copyfile(fixtures["taxonomy_db"], database)
    database.chmod(0o400)
    require(sha(database) == fixtures["taxonomy_sha256"], "copied taxonomy changed")
    ordered = sorted(fixtures["fixtures"].values(), key=lambda row: row["species"])
    metadata = assets / "reviewed_genomes.tsv"
    with metadata.open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["species", "taxid", "genome_sha256"])
        for row in ordered:
            require(0 < Path(row["genome"]).stat().st_size <= 1100000 and sha(row["genome"]) == row["sha256"], "small fixture changed")
            target = genomes / (row["species"] + ".fa")
            shutil.copyfile(row["genome"], target)
            target.chmod(0o400)
            writer.writerow([row["species"], row["taxid"], row["sha256"]])
    metadata.chmod(0o400)
    env = {key: value for key, value in os.environ.items() if key.startswith("SLURM_") or key in
           {"HOME", "USER", "LOGNAME", "LANG", "LC_ALL", "CUDA_VISIBLE_DEVICES"}}
    env.update(PATH="/usr/local/bin:/usr/bin:/bin", PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1",
               gg_workspace_dir=str(workspace), gg_container_image_path=manifest["image"],
               GG_CONTAINER_PROJECT_ROOT_BIND=str(root) + ":" + str(root), GG_COMMON_TMP_ROOT="workspace")
    overrides = {"GENOME_ANNOTATION_BACKEND": "annevo", "ANNEVO_MODEL_PROFILE": "auto",
                 "ANNEVO_GENOME_METADATA": str(metadata), "ANNEVO_TAXONOMY_DB": str(database),
                 "ANNEVO_RUNTIME_MANIFEST": manifest["runtime_locator"], "ANNEVO_BATCH_SIZE": "32",
                 "ANNEVO_SEQUENCE_THRESHOLD": "25", "RUN_COLLECT_GFF_INFO": "1", "RUN_CDS_FX2TAB": "1",
                 "RUN_GENOME_FX2TAB": "1", "RUN_CDSKIT_LOCALIZE": "0", "RUN_SCAFFOLD_TAXONOMY": "0",
                 "RUN_MULTISPECIES_SUMMARY": "0", "DELETE_TMP_DIR": "0"}
    env.update({"GG_GENOME_ANNOTATION_" + key: value for key, value in overrides.items()})
    entrypoint = source / "workflow/gg_genome_annotation_entrypoint.sh"
    results, generations, acquired = [], set(), set()
    for index, row in enumerate(ordered, 1):
        species, profile = row["species"], row["profile"]
        env.update(SLURM_ARRAY_TASK_ID=str(index), GG_GENOME_ANNOTATION_ANNEVO_NO_DOWNLOAD="0")
        label = "entrypoint_cold_" + profile
        invoke(root, label, ["bash", str(entrypoint)], env)
        require(": Ending" in (root / "logs" / (label + ".out")).read_text(), "entrypoint did not finish")
        prediction = workspace / "output/genome_prediction/annevo" / species
        receipt = json.loads((prediction / "native_prediction.json").read_text())
        runtime = json.loads((prediction / "runtime.json").read_text())
        selection = json.loads((prediction / "model_selection.json").read_text())
        require(receipt["model_profile"] == selection["model_profile"] == runtime["model_profile"] == profile
                and receipt["gpu_probe"]["gpu_verified"] is True and receipt["genome_sha256"] == row["sha256"], "selected CUDA prediction identity changed")
        generations.add(runtime["generation"])
        acquired.add(profile)
        require({p.name for p in cache.rglob("ANNEVO_*.pt")} == {"ANNEVO_" + p + ".pt" for p in acquired}, "unselected model downloaded")
        candidates = prediction / "candidate_inputs.tsv"
        with candidates.open() as stream:
            require(len(list(csv.DictReader(stream, delimiter="\t"))) == 1, "candidate export missing")
        required = [prediction / (species + ".gff3"), prediction / "native_prediction.json", candidates,
                    workspace / ("output/species_gff_info/" + species + "_gff_info.tsv"),
                    workspace / ("output/species_cds_fx2tab/" + species + "_fx2tab_cds.tsv"),
                    workspace / ("output/species_genome_fx2tab/" + species + "_fx2tab_genome.tsv")]
        require(all(p.is_file() and p.stat().st_size > 0 for p in required), "enabled downstream output missing")
        before = {str(p): sha(p) for p in required}
        env["GG_GENOME_ANNOTATION_ANNEVO_NO_DOWNLOAD"] = "1"
        warm_label = "entrypoint_warm_" + profile
        invoke(root, warm_label, ["bash", str(entrypoint)], env)
        warm_log = (root / "logs" / (warm_label + ".out")).read_text()
        require(": Ending" in warm_log and "Skipped: ANNEVO native genome prediction" in warm_log, "unchanged native prediction was not reused")
        require(before == {str(p): sha(p) for p in required}, "warm entrypoint changed completed artifacts")
        formats = list((prediction / "formatted").iterdir())
        require(len(formats) == 1, "unexpected formatter generation")
        formatted = formats[0]
        prefix = [manifest["container"], "exec", "--containall", "--no-home", "--bind",
                  str(root) + ":" + str(root) + "," + str(source) + ":" + str(source) + ":ro", manifest["image"], "python", "-B"]
        common = ["--species-cds-dir", str(formatted / "species_cds"), "--species-summary", str(formatted / "species.tsv"), "--nthreads", "1"]
        mapping, longest = root / "receipts" / (profile + ".mapping.json"), root / "receipts" / (profile + ".longest.json")
        invoke(root, "entrypoint_mapping_" + profile, prefix + [str(source / "workflow/support/validate_cds_gff_mapping.py"), *common,
               "--species-gff-dir", str(formatted / "species_gff"), "--species-genome-dir", str(formatted / "species_genome"),
               "--stats-output", str(mapping)], env, timeout=180)
        invoke(root, "entrypoint_longest_" + profile, prefix + [str(source / "workflow/support/validate_longest_cds_selection.py"), *common,
               "--summary-outputs-only", "--stats-output", str(longest)], env, timeout=180)
        mapped, chosen = json.loads(mapping.read_text()), json.loads(longest.read_text())
        require(mapped["species_checked"] == mapped["species_passed"] == chosen["species_checked"] == chosen["species_passed"] == 1
                and mapped["phase_conflicts_total"] == mapped["fallback_ids_total"] == 0
                and mapped["cds_ids_total"] == mapped["mapped_ids_total"] > 0, "CDS/GFF/longest QC failed")
        results.append({"profile": profile, "species": species, "taxid": row["taxid"], "mapped_cds": mapped["mapped_ids_total"],
                        "cold_entrypoint_verified": True, "warm_no_download_verified": True, "artifact_sha256": before,
                        "model_sha256": receipt["model_sha256"], "gpu_probe": receipt["gpu_probe"]})
    require(len(generations) == 1, "six profiles did not share one runtime")
    require(any(workspace.glob("output/versions/container.*.versions.log")), "successful version dump missing")
    require(all(sha(genomes / (row["species"] + ".fa")) == row["sha256"] for row in ordered), "curated qualification inputs changed")
    result = {"core_entrypoint_verified": True, "profiles_passed": 6, "entrypoint_runs_passed": 12, "cases": results,
              "runtime_generation": generations.pop(), "overrides": overrides, "scientific_completion_verified": False}
    write_json(root / "receipts/entrypoint_qualification.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    os.umask(0o077)
    qualify(json.loads(args.manifest.read_text()))
