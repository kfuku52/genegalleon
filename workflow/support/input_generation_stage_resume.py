"""Content-verified format/validation checkpoints for native species workers."""

import argparse
import contextlib
import csv
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from input_generation_array_state import (
    FreshDigestBatch,
    atomic_json,
    claim_workspace,
    digest,
    digest_paths,
    load_plan,
    namespace_path,
    prepared,
    verify_receipt,
)
from performance_metrics import count, measure
from shared_namespace_lock import inspect_lock, namespace_lock

FORMAT_PARAMETERS = ("provider", "gene_grouping_mode", "gff_repair_mode", "strict", "genetic_code",
                     "require_cds", "require_gff", "require_genome")
PHASE_READER_TIMEOUT = 30


def parameters(settings, stage, format_contract_version):
    result = {key: str(settings.get(key, "0" if key.startswith("require_") else "1" if key == "genetic_code" else ""))
              for key in FORMAT_PARAMETERS}
    result["format_contract_version"] = str(format_contract_version)
    if stage == "validate":
        result.update(validation_contract_version="4", run_validate_inputs=str(settings["run_validate_inputs"]))
    return result


def checkpoint_path(root, species, stage):
    return root / "tmp/stage_checkpoints" / f"{stage}.{species}.json"


def context(plan_path, index, root, stage, *, namespace_root=None):
    plan = load_plan(plan_path)
    task = plan["tasks"][index - 1]
    settings = json.loads(Path(str(plan_path) + ".settings.json").read_text())
    meta = json.loads((root / "tmp/task_meta_shards" / f"{index}.json").read_text())
    if meta["species_prefix"] != task["species_prefix"] or meta["task_index"] != index:
        raise ValueError("Resume metadata belongs to another species/task")
    paths = {"cds": meta["cds_output_path"], "stats": str(root / "tmp/task_stats_shards" / f"{index}.json"),
             "summary": str(root / "tmp/species_summary_shards" / f"{index}.tsv")}
    if meta.get("gff_output_path"):
        paths["gff"] = meta["gff_output_path"]
    if stage == "format":
        if meta.get("genome_output_path"):
            paths["genome"] = meta["genome_output_path"]
        paths.update({key: meta[key] for key in ("cds_path", "gff_path", "genome_path", "gbff_path") if meta.get(key)})
    else:
        paths["ownership_qc"] = str(root / "tmp/task_stats_shards" / f"{index}.longest.json")
        if "gff" in paths:
            paths["mapping_qc"] = str(root / "tmp/task_stats_shards" / f"{index}.mapping.json")
    paths = {label: str(namespace_path(path, namespace_root)) for label, path in paths.items()}
    return task, settings, meta, paths


def snapshot(paths, *, batch=None):
    if any(not Path(path).is_file() or not Path(path).stat().st_size for path in paths.values()):
        raise ValueError("Stage resume requires all declared nonempty files")
    hashes = (batch.read if batch is not None else digest_paths)(paths.values())
    return {label: {"path": path, "sha256": hashes[path]} for label, path in paths.items()}


def record(plan_path, index, root, stage, format_contract_version, *, expected_hashes=None, batch=None):
    task, settings, _, paths = context(plan_path, index, root, stage)
    if stage == "validate" and settings["run_validate_inputs"] != "1":
        raise ValueError("Disabled validation cannot produce a successful validation checkpoint")
    files = snapshot(paths, batch=batch)
    if expected_hashes is not None and any(files[label]["sha256"] != value for label, value in expected_hashes.items()):
        raise ValueError("Copied stage output differs from its verified source: " + task["species_prefix"])
    if batch is not None:
        batch.check()
    atomic_json(checkpoint_path(root, task["species_prefix"], stage), {
        "schema_version": 1, "stage": stage, "species": task["species_prefix"],
        "parameters": parameters(settings, stage, format_contract_version), "files": files,
    })


def verified_snapshot(plan_path, index, root, stage, format_contract_version, *, namespace_root=None, batch=None):
    try:
        task, settings, _, paths = context(plan_path, index, root, stage, namespace_root=namespace_root)
        saved = json.loads(checkpoint_path(root, task["species_prefix"], stage).read_text())
        compatible = (saved.get("schema_version") == 1 and saved.get("stage") == stage
                and saved.get("species") == task["species_prefix"]
                and saved.get("parameters") == parameters(settings, stage, format_contract_version))
        if not compatible:
            return None
        observed = snapshot(paths, batch=batch)
        expected = {label: {**item, "path": str(namespace_path(item["path"], namespace_root))}
                    for label, item in saved["files"].items()}
        return observed if expected == observed else None
    except (OSError, ValueError, KeyError, TypeError, IndexError):
        return None


def valid(plan_path, index, root, stage, format_contract_version, *, namespace_root=None):
    return verified_snapshot(plan_path, index, root, stage, format_contract_version,
                             namespace_root=namespace_root) is not None


def copy_atomic(source, destination, *, expected_sha256=None):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        if expected_sha256 is None:
            shutil.copyfile(source, temporary)
        else:
            hashed = hashlib.sha256()
            with Path(source).open("rb") as src, temporary.open("wb") as dst:
                for chunk in iter(lambda: src.read(1024 * 1024), b""):
                    dst.write(chunk)
                    hashed.update(chunk)
                    count("copy_bytes", len(chunk))
                    count("sha256_bytes", len(chunk))
            count("sha256_reads")
            if hashed.hexdigest() != expected_sha256:
                raise ValueError("Copied format output differs from its verified source: " + str(source))
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def reject_output_overlap(destinations, protected):
    protected = [Path(path) for path in protected]
    locations = {path.resolve() for path in protected}
    identities = set()
    for path in protected:
        try:
            info = path.stat()
        except FileNotFoundError:
            continue
        identities.add((info.st_dev, info.st_ino))
    for path in destinations:
        destination = Path(path)
        try:
            info = destination.stat()
        except FileNotFoundError:
            identity = None
        else:
            identity = (info.st_dev, info.st_ino)
        if destination.resolve() in locations or identity in identities:
            raise ValueError("Imported output paths overlap donor or raw inputs: " + str(destination))


def read_verified_summary(path, *, expected_sha256):
    data = Path(path).read_bytes()
    count("sha256_bytes", len(data))
    count("sha256_reads")
    if hashlib.sha256(data).hexdigest() != expected_sha256:
        raise ValueError("Summary differs from its verified source: " + str(path))
    return data


def import_fx2tab(source_root, target_root, species, source_cds, target_cds, source_settings, target_settings):
    if target_settings.get("run_cds_fx2tab") != "1":
        return
    manifest = source_root / "artifact_provenance" / f"fx2tab.{species}.json"
    if not manifest.is_file():
        return
    payload = json.loads(manifest.read_text())
    expected_parameters = dict(length="yes", name="yes", gc="yes", gc_skew="yes", only_id="yes")
    if (payload.get("schema_version") != 1 or payload.get("step") != "input_generation_cds_fx2tab"
            or str(payload.get("family_id")) != species or payload.get("parameters") != expected_parameters):
        raise ValueError("Unsupported fx2tab provenance: " + species)
    source = Path(source_settings["species_cds_fx2tab_dir"]) / f"{species}_fx2tab_cds.tsv"
    destination = Path(target_settings["species_cds_fx2tab_dir"]) / source.name
    reject_output_overlap([destination], [source, source_cds, target_cds])
    entries = {item["label"]: item for group in ("inputs", "outputs") for item in payload[group]}
    hashes = digest_paths([str(source_cds), str(source)])
    if (hashes[str(source_cds)] != entries["species_cds"]["sha256"]
            or hashes[str(source)] != entries["fx2tab"]["sha256"]):
        raise ValueError("fx2tab output/input differs from its successful provenance: " + species)
    copy_atomic(source, destination)
    if digest(destination) != hashes[str(source)]:
        raise ValueError("fx2tab output changed during import: " + species)
    command = [sys.executable, str(Path(__file__).with_name("artifact_provenance.py")), "record",
               "--manifest", str(target_root / "artifact_provenance" / manifest.name),
               "--step", "input_generation_cds_fx2tab", "--family-id", species,
               "--logical-root", str(target_root.parent / ".gg_global_artifacts"),
               "--workspace-root", str(target_root.parent.parent),
               "--input", "species_cds=" + str(target_cds), "--output", "fx2tab=" + str(destination)]
    for key, value in expected_parameters.items():
        command.extend(["--parameter", key + "=" + value])
    subprocess.run(command, check=True)


def import_busco(source_root, target_root, species, source_cds, target_cds, source_settings, target_settings):
    """Reuse only the identical lineage/mode contract and freshly verified CDS."""
    if (source_settings.get("run_species_busco") != "1" or target_settings.get("run_species_busco") != "1"
            or source_settings.get("busco_lineage") != target_settings.get("busco_lineage")):
        return False
    manifest = source_root / "artifact_provenance" / f"busco.{species}.json"
    if not manifest.is_file():
        return False
    batch = FreshDigestBatch()
    batch.read([manifest, source_root / "tmp/busco_lineage.resolved.txt", target_root / "tmp/busco_lineage.resolved.txt"])
    source_lineage = (source_root / "tmp/busco_lineage.resolved.txt").read_text().strip()
    target_lineage = (target_root / "tmp/busco_lineage.resolved.txt").read_text().strip()
    if not source_lineage or source_lineage != target_lineage:
        return False
    expected = dict(busco_lineage_request=target_settings["busco_lineage"], busco_lineage_resolved=target_lineage,
                    busco_mode="transcriptome", evalue="1e-03", limit="20")
    payload = json.loads(manifest.read_text())
    if (payload.get("schema_version") != 1 or payload.get("step") != "input_generation_species_busco"
            or payload.get("family_id") != species or payload.get("parameters") != expected):
        raise ValueError("Unsupported BUSCO provenance: " + species)
    inputs, outputs = payload["inputs"], payload["outputs"]
    output_roles = {item["label"] for item in outputs}
    if ([item["label"] for item in inputs] != ["species_cds"]
            or output_roles not in ({"busco_full", "busco_short"}, {"busco_full", "busco_short", "busco_single_copy"})
            or len(outputs) != len(output_roles)):
        raise ValueError("Unsupported BUSCO artifact roles: " + species)
    sources = {label: Path(source_settings["species_busco_" + label + "_dir"]) /
               (species + ".busco." + ("full.tsv" if label == "full" else "short.txt")) for label in ("full", "short")}
    destinations = {label: Path(target_settings["species_busco_" + label + "_dir"]) / path.name
                    for label, path in sources.items()}
    if "busco_single_copy" in output_roles:
        sources["single_copy"] = Path(source_settings["species_busco_full_dir"]) / "single_copy" / (species + ".json.gz")
        destinations["single_copy"] = Path(target_settings["species_busco_full_dir"]) / "single_copy" / (species + ".json.gz")
    reject_output_overlap(destinations.values(), [*sources.values(), source_cds, target_cds])
    hashes = batch.read([source_cds, target_cds, *sources.values()])
    output_hashes = {item["label"]: item["sha256"] for item in outputs}
    if (hashes[str(source_cds)] != inputs[0]["sha256"] or hashes[str(target_cds)] != inputs[0]["sha256"]
            or any(hashes[str(path)] != output_hashes["busco_" + label] for label, path in sources.items())):
        raise ValueError("BUSCO output/input differs from its successful provenance: " + species)
    batch.check()
    for label, path in sources.items():
        copy_atomic(path, destinations[label], expected_sha256=hashes[str(path)])
    batch.check()
    if any(digest(destinations[label]) != hashes[str(path)] for label, path in sources.items()):
        raise ValueError("BUSCO output changed during import: " + species)
    command = [sys.executable, str(Path(__file__).with_name("artifact_provenance.py")), "record",
               "--manifest", str(target_root / "artifact_provenance" / manifest.name),
               "--step", "input_generation_species_busco", "--family-id", species,
               "--logical-root", str(target_root.parent / ".gg_global_artifacts"),
               "--workspace-root", str(target_root.parent.parent), "--input", "species_cds=" + str(target_cds)]
    for label, path in destinations.items():
        command.extend(["--output", "busco_" + label + "=" + str(path)])
    for key, value in expected.items():
        command.extend(["--parameter", key + "=" + value])
    subprocess.run(command, check=True)
    batch.check()
    return True


def import_stages(args):
    with measure("stage_import"):
        if getattr(args, "source_only", False):
            return _import_stages(args)
        # Shared owners briefly hold the registration gate too. Retry that
        # contention while continuing to exclude prepare/finalize writers.
        with namespace_lock(args.root / ".array-phase.lock", exclusive=False,
                            timeout=PHASE_READER_TIMEOUT) as acquired:
            if not acquired:
                raise ValueError("Target input-generation workspace has active prepare/finalize")
            return _import_stages(args)


def _import_stages(args):
    """Import native proofs, never infer successful validation from output existence."""
    source_plan = args.source_plan.resolve(strict=True)
    if source_plan == args.task_plan.resolve():
        raise ValueError("Stage import requires a separate frozen plan/output workspace")
    source_root = args.source_root.resolve(strict=True)
    source_namespace = source_root.parent.parent
    if source_root == args.root.resolve():
        raise ValueError("Stage import cannot share the donor output workspace")
    source_only = getattr(args, "source_only", False)
    with namespace_lock(source_root / ".array-phase.lock", exclusive=source_only,
                        nonblocking=source_only, timeout=PHASE_READER_TIMEOUT) as acquired:
        if not acquired:
            raise ValueError("Source input-generation workspace has active workers/shared stages")
        if digest(source_plan) != args.source_plan_sha256:
            raise ValueError("Source plan differs from the sealed resume SHA-256")
        source = load_plan(source_plan)
        target = load_plan(args.task_plan)
        donor_indices = {task["species_prefix"]: i for i, task in enumerate(source["tasks"], 1)}
        selected = getattr(args, "task_index", None)
        donor_index = (donor_indices.get(target["tasks"][selected - 1]["species_prefix"])
                       if selected is not None and not source_only else None)
        if not prepared(source_plan, namespace_root=source_namespace, task_index=donor_index):
            raise ValueError("Source prepare/settings evidence is missing or stale")
        claim_workspace(source_plan, source_root, namespace_root=source_namespace)
        claim_workspace(args.task_plan, args.root)
        source_settings = json.loads(Path(str(source_plan) + ".settings.json").read_text())
        source_settings = {key: str(namespace_path(value, source_namespace))
                           if isinstance(value, str) and (value == "/workspace" or value.startswith("/workspace/")) else value
                           for key, value in source_settings.items()}
        target_settings = json.loads(Path(str(args.task_plan) + ".settings.json").read_text())
        if parameters(source_settings, "format", args.format_contract_version) != parameters(target_settings, "format", args.format_contract_version):
            raise ValueError("Source/target formatting or required-output parameters differ")
        if source_only:
            print(json.dumps({"source_verified": True, "source_plan_sha256": args.source_plan_sha256}))
            return
        donor_indices = {task["species_prefix"]: i for i, task in enumerate(source["tasks"], 1)}
        imported = []
        skipped = []
        for index, task in enumerate(target["tasks"], 1):
            if getattr(args, "task_index", None) is not None and index != args.task_index:
                continue
            species = task["species_prefix"]
            old_index = donor_indices.get(species)
            if old_index is None:
                skipped.append(species)
                continue
            donor_task = source["tasks"][old_index - 1]
            if any(task.get(key) != donor_task.get(key) for key in ("provider", "species_key")):
                raise ValueError("Source/target species identity differs: " + species)
            # A donor worker holds this same task lock while writing. Root
            # shared ownership excludes prepare/finalize; unrelated species
            # and other import readers may progress concurrently.
            with namespace_lock(Path(str(source_plan) + ".locks") / f"{old_index}.lock",
                                exclusive=False, nonblocking=True) as task_acquired:
                if not task_acquired:
                    raise ValueError("Source input-generation workspace has active workers for " + species)
                target_lock = Path(str(args.task_plan) + ".locks") / f"{index}.lock"
                token = getattr(args, "target_lock_token", None)
                if token:
                    owner = inspect_lock(target_lock)["exclusive"]
                    if owner is None or owner["token"] != token:
                        raise ValueError("Target task lock ownership changed")
                with (contextlib.nullcontext(True) if token else
                      namespace_lock(target_lock, exclusive=True, nonblocking=True)) as target_acquired:
                    if not target_acquired:
                        raise ValueError("Target species task already has an active owner")
                    _import_one(args, source_plan, source_root, source_namespace, source_settings,
                                target_settings, index, task, donor_task, old_index, imported, skipped)
        if digest(source_plan) != args.source_plan_sha256:
            raise ValueError("Source plan changed during stage import")
        report = {"imported": imported, "without_verified_format": skipped}
        print(json.dumps(report))
        return report


def start_worker(args):
    """Resolve the worker once, sharing only fresh checks before any writes."""
    target_lock = Path(str(args.task_plan) + ".locks") / f"{args.task_index}.lock"
    owner = inspect_lock(target_lock)["exclusive"]
    if owner is None or owner["token"] != args.target_lock_token:
        raise ValueError("Target task lock ownership changed")
    batch = FreshDigestBatch()
    current = (not args.overwrite and verified_snapshot(args.task_plan, args.task_index, args.root,
               "format", args.format_contract_version, batch=batch) is not None)
    if not current and not args.overwrite and args.source_plan is not None:
        sources = [(args.source_plan, args.source_root, args.source_plan_sha256)]
        if getattr(args, "fallback_source_plan", None) is not None:
            sources.append((args.fallback_source_plan, args.fallback_source_root, args.fallback_source_plan_sha256))
        for plan, root, sha in sources:
            donor_args = argparse.Namespace(**vars(args))
            donor_args.source_plan, donor_args.source_root, donor_args.source_plan_sha256 = plan, root, sha
            donor_args.source_only = False
            result = import_stages(donor_args)
            if result["imported"]:
                return
    from run_input_generation_task import build_arg_parser, describe_task
    settings = json.loads(Path(str(args.task_plan) + ".settings.json").read_text())
    command = ["--task-plan", str(args.task_plan), "--task-index", str(args.task_index),
               "--task-meta-output", str(args.root / "tmp/task_meta_shards" / f"{args.task_index}.json"),
               "--download-timeout", str(args.download_timeout)]
    for key in ("species_cds_dir", "species_gff_dir", "species_genome_dir"):
        command.extend(["--" + key.replace("_", "-"), settings[key]])
    for header in args.http_header:
        command.extend(["--http-header", header])
    command.extend(["--auth-bearer-token-env", args.auth_bearer_token_env])
    # An incompatible checkpoint may have performed partial checks. Discard
    # them rather than carrying their fences across an attempted import.
    describe_task(build_arg_parser().parse_args(command), batch=batch if current else FreshDigestBatch())


def _import_one(args, source_plan, source_root, source_namespace, source_settings,
                target_settings, index, task, donor_task, old_index, imported, skipped):
    species = task["species_prefix"]
    before_copy = FreshDigestBatch()
    original_files = verified_snapshot(source_plan, old_index, source_root, "format", args.format_contract_version,
                                       namespace_root=source_namespace, batch=before_copy)
    format_valid = original_files is not None
    # Completion is a legacy alternative, not a prerequisite for a
    # current format checkpoint (and may contain large BUSCO outputs).
    complete = (not format_valid and verify_receipt(source_plan, old_index, args.source_plan_sha256,
                                                   namespace_root=source_namespace))
    if not complete and not format_valid:
        skipped.append(species)
        return
    _, _, _, old_paths = context(source_plan, old_index, source_root, "format", namespace_root=source_namespace)
    receipt = json.loads(Path(str(source_plan) + f".completed/{old_index}.json").read_text()) if complete else {}
    if receipt:
        receipt["files"] = {str(namespace_path(path, source_namespace)): value for path, value in receipt["files"].items()}
    if complete and not format_valid:
        # Old workers bound their metadata/shards and formatted outputs
        # to completion, and their format manifest declared this contract.
        manifest = json.loads((source_root / "artifact_provenance" / f"format.{species}.json").read_text())
        old_parameters = manifest.get("parameters", {})
        if (manifest.get("step") != "input_generation_format"
                or str(manifest.get("family_id")) != species
                or str(old_parameters.get("format_contract_version")) != args.format_contract_version):
            raise ValueError("Unsupported legacy format provenance: " + species)
        if any(str(old_parameters.get(key, "")) != str(source_settings[key])
               for key in ("provider", "strict", "gene_grouping_mode", "gff_repair_mode")):
            raise ValueError("Legacy format parameters differ from their frozen settings: " + species)
        if any(path not in receipt["files"] for path in old_paths.values()):
            raise ValueError("Legacy completion does not certify every format shard: " + species)
    if original_files is None:
        original_files = snapshot(old_paths, batch=before_copy)
    if format_valid:
        saved_files = json.loads(checkpoint_path(source_root, species, "format").read_text())["files"]
        saved_files = {label: {**item, "path": str(namespace_path(item["path"], source_namespace))}
                       for label, item in saved_files.items()}
        if original_files != saved_files:
            raise ValueError("Source format files changed after verification: " + species)
    elif any(item["sha256"] != receipt["files"][item["path"]] for item in original_files.values()):
        raise ValueError("Source completed files changed after verification: " + species)
    meta_file = args.root / "tmp/task_meta_shards" / f"{index}.json"
    protected_paths = [*old_paths.values(), source_root / "tmp/task_meta_shards" / f"{old_index}.json",
                       source_root / "tmp/task_stats_shards" / f"{old_index}.longest.json",
                       source_root / "tmp/task_stats_shards" / f"{old_index}.mapping.json"]
    reject_output_overlap([meta_file], protected_paths)
    from run_input_generation_task import build_arg_parser, describe_task
    describe_args = build_arg_parser().parse_args([
        "--task-plan", str(args.task_plan), "--task-index", str(index), "--describe-only",
        "--task-meta-output", str(meta_file),
        *[part for key in ("species_cds_dir", "species_gff_dir", "species_genome_dir")
          for part in ("--" + key.replace("_", "-"), target_settings[key])]])
    describe_task(describe_args, batch=before_copy)
    _, _, _, new_paths = context(args.task_plan, index, args.root, "format")
    old_raw = {key: path for key, path in old_paths.items() if key.endswith("_path")}
    new_raw = {key: path for key, path in new_paths.items() if key.endswith("_path")}
    if old_raw.keys() != new_raw.keys():
        raise ValueError("Source/target raw input roles differ: " + species)
    raw_hashes = before_copy.read(new_raw.values())
    if any(original_files[key]["sha256"] != raw_hashes[new_raw[key]] for key in old_raw):
        raise ValueError("Source/target raw input content differs: " + species)
    # Separate roots do not imply separate physical outputs: custom directories,
    # symlinks and hard links must not turn a donor reader into a writer.
    reject_output_overlap([path for label, path in new_paths.items() if not label.endswith("_path")],
                          [*protected_paths, *new_raw.values()])
    validation_files = None
    if source_settings.get("run_validate_inputs") == "1" and target_settings.get("run_validate_inputs") == "1":
        validation_files = verified_snapshot(source_plan, old_index, source_root, "validate",
                                             args.format_contract_version, namespace_root=source_namespace,
                                             batch=before_copy)
    if validation_files is None and "gff" in old_paths and "genome" in old_paths:
        from format_species_annotation.reference import validate_gff_genome_references
        try:
            validate_gff_genome_references(old_paths["gff"], old_paths["genome"])
        except ValueError:
            # A completed formatter is not proof that its reference contract
            # passed. Reformat this failed input under the fixed runtime.
            skipped.append(species)
            return
    if validation_files is not None:
        _, _, _, validation_paths = context(args.task_plan, index, args.root, "validate")
        reject_output_overlap([validation_paths[label] for label in ("ownership_qc", "mapping_qc")
                               if label in validation_files],
                              [*protected_paths, *new_raw.values(),
                               *(item["path"] for item in validation_files.values())])
    before_copy.check()
    # Writes and parsing start a new boundary. No pre-transfer digest is used
    # as authority for the post-transfer content verification below.
    for label in ("cds", "gff", "genome", "stats"):
        if label in old_paths:
            if label not in new_paths:
                raise ValueError("Source/target formatted output roles differ: " + species)
            copy_atomic(old_paths[label], new_paths[label], expected_sha256=original_files[label]["sha256"])
    replacements = {old_paths[label]: new_paths[label] for label in old_paths if label in new_paths}
    summary = read_verified_summary(old_paths["summary"], expected_sha256=original_files["summary"]["sha256"])
    reader = csv.DictReader(io.StringIO(summary.decode("utf-8"), newline=""), delimiter="\t")
    fields = reader.fieldnames
    rows = [{key: replacements.get(str(namespace_path(value, source_namespace)), value)
             for key, value in row.items()} for row in reader]
    if len(rows) != 1 or rows[0].get("species_prefix") != species:
        raise ValueError("Source summary shard does not identify exactly one requested species")
    destination = Path(new_paths["summary"])
    destination.parent.mkdir(parents=True, exist_ok=True)
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    summary_bytes = buffer.getvalue().encode("utf-8")
    destination.write_bytes(summary_bytes)
    if validation_files is not None:
        for label in ("ownership_qc", "mapping_qc"):
            if label in validation_files:
                copy_atomic(validation_files[label]["path"], validation_paths[label],
                            expected_sha256=validation_files[label]["sha256"])
    after_copy = FreshDigestBatch()
    if snapshot(old_paths, batch=after_copy) != original_files:
        raise ValueError("Source format files changed during import: " + species)
    if validation_files is not None:
        _, _, _, old_validation_paths = context(source_plan, old_index, source_root, "validate",
                                               namespace_root=source_namespace)
        if snapshot(old_validation_paths, batch=after_copy) != validation_files:
            raise ValueError("Source validation files changed during import: " + species)
    summary_sha = hashlib.sha256(summary_bytes).hexdigest()
    record(args.task_plan, index, args.root, "format", args.format_contract_version,
           expected_hashes={**{label: item["sha256"] for label, item in original_files.items()},
                            "summary": summary_sha}, batch=after_copy)
    if validation_files is not None:
        record(args.task_plan, index, args.root, "validate", args.format_contract_version,
               expected_hashes={**{label: item["sha256"] for label, item in validation_files.items()},
                                "summary": summary_sha}, batch=after_copy)
    after_copy.check()
    import_fx2tab(source_root, args.root, species, old_paths["cds"], new_paths["cds"],
                  source_settings, target_settings)
    busco = import_busco(source_root, args.root, species, old_paths["cds"], new_paths["cds"],
                         source_settings, target_settings)
    imported.append({"species": species, "validation": validation_files is not None, "busco": busco})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("check", "record", "import", "check-source", "start-worker"))
    parser.add_argument("--task-plan", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--format-contract-version", required=True)
    parser.add_argument("--stage", choices=("format", "validate"))
    parser.add_argument("--source-plan", type=Path)
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--source-plan-sha256")
    parser.add_argument("--fallback-source-plan", type=Path)
    parser.add_argument("--fallback-source-root", type=Path)
    parser.add_argument("--fallback-source-plan-sha256")
    parser.add_argument("--target-lock-token", help="Existing native worker's target task lock ownership token")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--download-timeout", type=float, default=120)
    parser.add_argument("--http-header", action="append", default=[])
    parser.add_argument("--auth-bearer-token-env", default="")
    args = parser.parse_args()
    fallback = (args.fallback_source_plan, args.fallback_source_root, args.fallback_source_plan_sha256)
    if any(fallback) and (args.action != "start-worker" or not all(fallback) or args.source_plan is None
                         or args.fallback_source_plan.resolve() == args.source_plan.resolve()):
        parser.error("Fallback requires start-worker, a distinct source plan, root and sealed SHA-256")
    if args.action == "start-worker":
        if not args.task_index or not 1 <= args.task_index <= load_plan(args.task_plan)["task_count"]:
            parser.error("A valid task index is required")
        if args.source_plan is not None and not all((args.source_root, args.source_plan_sha256)):
            parser.error("Import requires the source root and sealed SHA-256")
        start_worker(args)
    elif args.action in ("import", "check-source"):
        if not all((args.source_plan, args.source_root, args.source_plan_sha256)):
            parser.error("Import requires the source plan, root and sealed SHA-256")
        if args.task_index is not None and not 1 <= args.task_index <= load_plan(args.task_plan)["task_count"]:
            parser.error("A valid task index is required")
        args.source_only = args.action == "check-source"
        import_stages(args)
    else:
        if not args.stage or not args.task_index or not 1 <= args.task_index <= load_plan(args.task_plan)["task_count"]:
            parser.error("A stage and valid task index are required")
        if args.action == "check":
            raise SystemExit(0 if valid(args.task_plan, args.task_index, args.root, args.stage, args.format_contract_version) else 1)
        record(args.task_plan, args.task_index, args.root, args.stage, args.format_contract_version)


if __name__ == "__main__":
    main()
