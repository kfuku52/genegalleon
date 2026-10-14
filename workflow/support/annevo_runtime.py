#!/usr/bin/env python3
"""Prepare optional ANNEVO assets and invoke the unmodified native predictor."""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import uuid
from pathlib import Path

from annevo_runtime_cache import prepare, runtime_environment, verify_ready
from annevo_runtime_download import atomic_json, file_sha
from genome_sequence_identity import genome_sequence_digest
from shared_namespace_lock import namespace_lock


def read_selected(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise RuntimeError("ANNEVO runtime receipt is missing or redirected")
    receipt = json.loads(path.read_text())
    if receipt.get("schema_version") != 1 or receipt.get("kind") != "genegalleon-annevo-installed-runtime":
        raise RuntimeError("Invalid ANNEVO selected runtime receipt")
    verify_ready(receipt)
    model = Path(receipt["model_path"])
    expected_model = receipt["manifest"]["models"][receipt["model_profile"]]
    expected_path = Path(receipt["cache_root"]) / "models" / expected_model["sha256"] / Path(expected_model["path"]).name
    if (model != expected_path or receipt["model_sha256"] != expected_model["sha256"]
            or receipt["model_size"] != expected_model["size"]):
        raise RuntimeError("ANNEVO selected model identity changed")
    if model.is_symlink() or not model.is_file() or model.stat().st_size != receipt["model_size"] or file_sha(model) != receipt["model_sha256"]:
        raise RuntimeError("ANNEVO model checksum changed")
    return receipt


def probe_gpu(receipt, env):
    code = """import json,torch
if not torch.cuda.is_available():
 raise SystemExit('ANNEVO requires the allocated GPU; enable SIF --nv or Docker GPU forwarding')
device=torch.cuda.current_device()
x=torch.ones((64,64),device='cuda');y=x@x;torch.cuda.synchronize()
assert y.is_cuda and y[0,0].item()==64
print(json.dumps({'gpu_verified':True,'torch_version':torch.__version__,
 'compiled_cuda':torch.version.cuda,'device':device,
 'device_name':torch.cuda.get_device_name(device)}))
"""
    result = subprocess.run([receipt["python"], "-E", "-s", "-B", "-c", code],
                            env=env, capture_output=True, text=True, timeout=180)
    if result.returncode:
        raise RuntimeError("ANNEVO allocated GPU probe failed: " + result.stderr[-4000:])
    data = json.loads(result.stdout.strip().splitlines()[-1])
    expected = receipt["manifest"]["cuda_profile"]
    if data.get("gpu_verified") is not True or any(data.get(k) != expected[k] for k in ("torch_version", "compiled_cuda")):
        raise RuntimeError("ANNEVO GPU runtime profile changed")
    return data


def gff_counts(path):
    counts = {}
    identifiers = set()
    parents = set()
    with Path(path).open() as stream:
        for line in stream:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or int(fields[3]) <= 0 or int(fields[4]) < int(fields[3]):
                raise RuntimeError("ANNEVO emitted invalid GFF coordinates")
            if fields[6] not in {"+", "-"} or (fields[2] == "CDS" and fields[7] not in {"0", "1", "2"}):
                raise RuntimeError("ANNEVO emitted invalid GFF strand or CDS phase")
            attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            if attrs.get("ID"):
                if attrs["ID"] in identifiers:
                    raise RuntimeError("ANNEVO emitted duplicate feature identifiers")
                identifiers.add(attrs["ID"])
            parents.update(attrs.get("Parent", "").split(","))
            counts[fields[2]] = counts.get(fields[2], 0) + 1
    if not counts.get("gene") or not counts.get("CDS") or parents - identifiers - {""}:
        raise RuntimeError("ANNEVO GFF lacks genes/CDS or has unresolved parents")
    return counts


def publish_pair(predicted, output, output_receipt, receipt):
    """Stage both files, then restore the previous pair if publication fails."""
    token = uuid.uuid4().hex
    stages = [p.with_name(p.name + ".staging-" + token) for p in (output, output_receipt)]
    backups = [p.with_name(p.name + ".backup-" + token) for p in (output, output_receipt)]
    moved, published = [], []
    try:
        shutil.copyfile(predicted, stages[0])
        atomic_json(stages[1], receipt)
        for destination, backup in zip((output, output_receipt), backups, strict=True):
            if destination.is_symlink() or (destination.exists() and not destination.is_file()):
                raise ValueError("ANNEVO publication target is redirected or nonregular")
            if destination.exists():
                os.replace(destination, backup)
                moved.append((destination, backup))
        for stage, destination in zip(stages, (output, output_receipt), strict=True):
            os.replace(stage, destination)
            published.append(destination)
    except BaseException:
        for destination in published:
            destination.unlink(missing_ok=True)
        for destination, backup in reversed(moved):
            os.replace(backup, destination)
        raise
    finally:
        for temporary in stages:
            temporary.unlink(missing_ok=True)
    for backup in backups:
        backup.unlink(missing_ok=True)


def run(args):
    selected = Path(args.receipt)
    runtime = read_selected(selected)
    if args.model_profile != runtime["model_profile"]:
        raise ValueError("ANNEVO selected model profile changed")
    for value in (args.threads, args.batch_size, args.sequence_threshold):
        if value <= 0:
            raise ValueError("ANNEVO CPU, batch and sequence thresholds must be positive")
    genome, output, output_receipt = map(Path, (args.genome, args.output, args.output_receipt))
    if genome.is_symlink() or not genome.is_file() or genome.stat().st_size == 0:
        raise ValueError("ANNEVO genome input must be a nonempty regular file")
    protected = {genome.resolve(), selected.resolve(), Path(runtime["model_path"]).resolve()}
    generation = Path(runtime["generation"]).resolve()
    if (output.resolve() == output_receipt.resolve()
            or any(p.resolve() in protected or p.resolve().is_relative_to(generation)
                   for p in (output, output_receipt))):
        raise ValueError("ANNEVO outputs must not overwrite inputs or the installed runtime")
    output.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    output_receipt.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    env = runtime_environment(runtime)
    env["OMP_NUM_THREADS"] = str(args.threads)
    env["MKL_NUM_THREADS"] = str(args.threads)
    original_sha, selected_sha = file_sha(genome), file_sha(selected)
    sequence_identity = genome_sequence_digest(genome)
    if file_sha(genome) != original_sha:
        raise RuntimeError("ANNEVO genome changed during sequence identity reading")
    with namespace_lock(output.with_name(output.name + ".publication.lock"), exclusive=True, timeout=86400):
        gpu = probe_gpu(runtime, env)
        # Keep native multiprocessing sockets short, irrespective of workspace path length.
        with tempfile.TemporaryDirectory(prefix="gg-annevo-", dir="/tmp") as temporary:
            work = Path(temporary)
            predicted = work / "prediction.gff3"
            native_genome = genome.absolute()
            compressed = genome.name.endswith(".gz")
            if compressed:
                native_genome = work / "genome.fa"
                with gzip.open(genome, "rb") as source, native_genome.open("xb") as target:
                    shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
            native_genome_sha = file_sha(native_genome)
            argv = [runtime["python"], "-E", "-s", "-B", runtime["native_entrypoint"],
                    "-g", str(native_genome), "-m", runtime["model_path"],
                    "-l", args.model_profile, "-o", str(predicted), "-t", str(args.threads),
                    "-s", str(args.sequence_threshold), "--tmp_path", str(work),
                    "--batch_size", str(args.batch_size), "--num_workers", str(min(2, args.threads)),
                    "--min_intron_length", "20", "--min_prot_length", "100", "--show_log"]
            env["TMPDIR"] = str(work)
            env["NUMBA_CACHE_DIR"] = str(work / "numba-cache")
            env["TORCH_HOME"] = str(work / "torch-cache")
            env["XDG_CACHE_HOME"] = str(work / "xdg-cache")
            log_root = output.with_name(output.name + ".logs")
            log_root.mkdir(mode=0o700, exist_ok=True)
            started = time.monotonic()
            with (log_root / "native.stdout.txt").open("w") as stdout, (log_root / "native.stderr.txt").open("w") as stderr:
                result = subprocess.run(argv, env=env, cwd=work, stdout=stdout, stderr=stderr)
            if result.returncode or not predicted.is_file() or predicted.stat().st_size == 0:
                raise RuntimeError("ANNEVO native prediction failed; retained logs: " + str(log_root))
            counts = gff_counts(predicted)
            read_selected(selected)
            if file_sha(selected) != selected_sha or file_sha(genome) != original_sha:
                raise RuntimeError("ANNEVO inputs changed during prediction")
            receipt = {"schema_version": 1, "kind": "genegalleon-annevo-native-prediction",
                       "native_returncode": result.returncode, "gpu_probe": gpu,
                       "runtime_digest": runtime["runtime_digest"], "cache_key": runtime["cache_key"],
                       "source_commit": runtime["source_commit"], "model_sha256": runtime["model_sha256"],
                       "model_profile": args.model_profile, "genome_sha256": original_sha,
                       "genome_sequence_sha256": sequence_identity,
                       "prepared_receipt_sha256": selected_sha, "argv": argv,
                       "native_genome_sha256": native_genome_sha, "genome_gzip": compressed,
                       "parameters": {"threads": args.threads, "batch_size": args.batch_size,
                                      "sequence_threshold": args.sequence_threshold,
                                      "num_workers": min(2, args.threads), "min_intron_length": 20, "min_prot_length": 100},
                       "GFF_sha256": file_sha(predicted), "GFF_bytes": predicted.stat().st_size,
                       "feature_counts": counts, "native_wall_seconds": time.monotonic() - started,
                       "scientific_completion_verified": False}
            publish_pair(predicted, output, output_receipt, receipt)
    return receipt


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    setup = commands.add_parser("prepare")
    setup.add_argument("--cache-root", required=True)
    setup.add_argument("--runtime-manifest", default="")
    setup.add_argument("--model-profile", default="Magnoliopsida")
    setup.add_argument("--receipt")
    setup.add_argument("--no-download", action="store_true")
    setup.add_argument("--refresh", action="store_true")
    predict = commands.add_parser("run")
    for name in ("receipt", "genome", "output", "output-receipt"):
        predict.add_argument("--" + name, required=True)
    predict.add_argument("--model-profile", default="Magnoliopsida")
    predict.add_argument("--threads", type=int, default=4)
    predict.add_argument("--batch-size", type=int, default=32)
    predict.add_argument("--sequence-threshold", type=int, default=25)
    return result


def main():
    args = parser().parse_args()
    try:
        if args.command == "prepare":
            value = prepare(args.cache_root, args.model_profile, args.runtime_manifest,
                            receipt_path=args.receipt, no_download=args.no_download, refresh=args.refresh)
        else:
            value = run(args)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError, tarfile.TarError) as exc:
        print("ANNEVO: " + str(exc), file=sys.stderr)
        return 1
    print(json.dumps(value, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
